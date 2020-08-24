using DataStructures
using InteractiveUtils
using JuMP
using Logging
# using MathOptInterface
using Printf
using SCIP

include("utils.jl")
include("bounds.jl")

function get_ysol(y, y_indexes, nodes, operA; num_result=1)
    ysol = SortedDict{Any,Any}()

    for n in nodes
        arr_val     = [JuMP.value(y[n,o], result=num_result) for o in operA if (n,o) in y_indexes]
        arr_oper    = [o for o in operA if (n,o) in y_indexes]

        if maximum(arr_val) > 0.1
            ysol[n] = arr_oper[argmax(arr_val)]
        end
    end

    ysol
end

function get_csol(c, nodes; num_result=1)
    SortedDict(n => JuMP.value(c[n]; result=num_result) for n in nodes)
end

function print_y(y, y_indexes, nodes, operA)
    out = "Print y\n"
    out *= "n\\oper\t"
    for o in operA
        out *= "$(o)\t"
    end
    out *= "\n"
    
    for n in nodes
        out *= @sprintf("%6d\t", n)
        for o in operA
            if (n,o) in y_indexes
                if JuMP.value(y[n,o]) < 1e-06
                    out *= ".\t"
                else
                    out *= @sprintf("%.1f\t", JuMP.value(y[n,o]))
                end
            else
                out *= " \t"
            end
        end
        out *= "\n"
    end

    out
end

function print_v(v, nodes, obs_indexes)
    out = "Print v\n"
    out *= "n\\obs\t"
    for i in obs_indexes
        out *= @sprintf("%5d\t", i)
    end
    out *= "\n"
    
    for n in nodes
        out *= @sprintf("%2d\t", n)
        for i in obs_indexes
            out *= @sprintf("%6.3f\t", JuMP.value(v[i,n]))
        end
        out *= "\n"
    end

    out
end

function print_c(c, y, nodes)
    out = "Print c\n"
    for n in nodes
        if JuMP.value(y[n,'C']) > 0.5
            out *= @sprintf("c[%d] = %.6f\n", n, JuMP.value(c[n]))
        end
    end

    out
end

function print_num_oper(num_oper, operA)
    out = "Print num_oper\n"
    for o in operA
        out *= @sprintf("%s: %d\n", o, round(JuMP.value(num_oper[o])))
    end
    out
end

function solve_MINLP(nodes, obs, operators;
                    print_all_solutions=false,
                    formulation="New-NR",
                    integer_constant=false,
                    ysol=nothing, 
                    ysol_dist_max=0, 
                    ysol_dist_min=0, 
                    ysol_fixlevel=-1,
                    csol=nothing,
                    scip_verblevel=3,
                    scip_time=10.0,
                    scip_absgap=0.0,
                    scip_nodes=-1,
                    scip_presolve=-1)
    """
    Solve a MINLP for symbolic regression.
    nodes:      set of nodes that we consider
    obs:        observations in array, assume that the last column is the dependent variable.
    operators:  set of operators that we consider
    
    print_all_solutions: print all solutions found in the log file (both optimal and nonoptimal) 

    formulation:    "Cozad":    Cozad's formulation with no optional cut such as redundancy elimination and no symmetry breacking
                    "New":      New formulation with no optional cut
                    "-CR":      Add Cozad's redundancy elimination cut
                    "-NR":      Add new redundancy elimination cut
        * I have not implemented symmetry breaking and implication constarint yet.

    (matters if ysol != nothing)
    ysol:                           solution of variable y in dict format (ysol[`nodeid`] = `operator`)
    ysol_dist_min, ysol_dist_max:   min/max distance from the current solution y
    ysol_fixlevel:                  if value is 0, fix all y solutions,
                                    if value is 1, fix all nonleaves,
                                    ... as value increases, fix less number of nodes 
    
    (matters if csol != nothing)
    csol:   solution of variable c in dict format (csol[`nodeid`] = `constant value`)
            if it is given, we fix csol

    scip_*: scip parameters 
    """
    

    nodes = OrderedSet(sort(collect(nodes)))
    @info   "Print nodes & operators\n" *
            "nodes:     $(nodes)\n" *
            "operators: $(operators)"

    ## Define useful variables from { nodes, obs, operators }
    leaves		= get_leaves(nodes)
    nleaves		= get_nonleaves(nodes)
    num_obs 	= size(obs, 1)
    num_indvars = size(obs, 2) - 1
    indvars		= OrderedSet(1:num_indvars)
    operA		= union(operators, indvars)
    operB       = intersect(operators, Set("+-*D"))
    operU       = intersect(operators, Set("REL"))
    operBU      = union(operB, operU)
    operL       = union(intersect(operators, Set("CP")), indvars)

    ## Other parameters (bounds, weights, ...)
    v_lb = isempty(intersect(operators, Set("EL"))) ? -1e+03 : -1e+02 
    v_ub = isempty(intersect(operators, Set("EL"))) ? 1e+03 : 1e+02
    v_lb, v_ub = get_vbounds(obs, operators)
    @info "bounds on v" v_lb v_ub

    (c_lb, c_ub, e_lb, e_ub) = (-2, 2, -1e+09, 1e+09)
    e_weight = ones(num_obs)
    lambda = 0
    num_active_ub = 7
    num_cst_ub = 1

    cst_abslb = 0.1
    mult_child_abslb = 0.01
    div_child_abslb = 0.01
    sqrt_rhs_lb = 0.01
    exp_rhs_lb = v_lb > 0 ? max(v_lb, log(v_lb)) : v_lb
    exp_rhs_ub = min(v_ub, log(v_ub))
    log_rhs_lb = max(v_lb, exp(v_lb))
    log_rhs_ub = min(v_ub, exp(v_ub))

    ## Create a model and set the solver's parameters (timelimit, gap, ...)
    @debug "Create a model and set the solver's parameters (timelimit, gap, ...)"
    optimizer = SCIP.Optimizer(
                    display_verblevel=scip_verblevel,           # default = 4, 0:5
                    limits_absgap=scip_absgap,                  # default = 0
                    limits_time=scip_time,  
                    limits_nodes=scip_nodes,   
                    presolving_maxrounds=scip_presolve,         # default = -1                
                    numerics_feastol=1e-07,                     # default = 1e-06
                    numerics_epsilon=1e-09)                     # default = 1e-09
    model = Model(() -> optimizer) 

    ## Decision variables
    @debug "Decision variables"
    # y: binary variable for the assignment of operators and operands
    y_indexes = union(OrderedSet([(n,o) for n in nleaves for o in operA]), 
                    OrderedSet([(n,o) for n in leaves for o in operL]))
    @variable(model, y[n in nodes, o in operA; (n,o) in y_indexes], Bin)

    # v: intermediate value at node for observation
    @variable(model, v_lb <= v[1:num_obs, nodes] <= v_ub)

    # c: constant value at node
    if integer_constant
        @variable(model, c_lb <= c[nodes] <= c_ub, Int)
    else
        @variable(model, c_lb <= c[nodes] <= c_ub)
    end

    # e: error = prediction - true
    @variable(model, e_lb <= e[1:num_obs] <= e_ub)

    ## Expressions
    @debug "Expressions"
    # num_oper: number of operators
    @expression(model, num_oper[o1 in operA], sum(y[n,o] for (n,o) in y_indexes if o==o1))

    ## Objective: Minimize (MSE) + (lambda)*(# of active_nodes)
    @debug "Objectve"
    @objective(model, Min, 1 / num_obs * sum(e.^2) + lambda * sum(y))

    ## Constrs for control
    @debug "Constrs for control"
    ## Constrs for k-neighbors search
    if !isnothing(ysol) 
        @info "Set start value" ysol
        for (n,o) in y_indexes
            if n in keys(ysol) && o == ysol[n]
                JuMP.set_start_value(y[n,o], 1)
            else
                JuMP.set_start_value(y[n,o], 0)
            end
        end

        if ysol_dist_max >= 0
            @info "Constrs for k-neighbors search" ysol ysol_dist_max ysol_dist_min
            @assert ysol_dist_max >= ysol_dist_min
            @assert ysol_dist_min >= 0
            @constraint(model, ysol_neighbors, 
                ysol_dist_min <=
                sum( (1-y[n,ysol[n]]) for n in nodes if n in keys(ysol) ) + 
                sum( sum(y[n,o] for o in operA if (n,o) in y_indexes) for n in nodes if !(n in keys(ysol)) ) <= ysol_dist_max)
        end

        if ysol_fixlevel >= 0
            @info "Fix ysol level = $(ysol_fixlevel)"

            nodes_fixed = keys(ysol)
            for i in 1:ysol_fixlevel
                nodes_fixed = get_nonleaves(nodes_fixed)
            end

            for n in nodes_fixed, o in operA
                if (n,o) in y_indexes
                    JuMP.is_binary(y[n,o]) ? JuMP.unset_binary(y[n,o]) : nothing
                    JuMP.is_fixed(y[n,o]) ? JuMP.unfix(y[n,o]) : nothing
                    JuMP.has_lower_bound(y[n,o]) ? JuMP.delete_lower_bound(y[n,o]) : nothing
                    JuMP.has_upper_bound(y[n,o]) ? JuMP.delete_upper_bound(y[n,o]) : nothing
                    JuMP.fix(y[n,o], o == ysol[n] ? 1 : 0)
                end
            end
        end
    end

    if !isnothing(csol) 
        @info "Fix csol" csol
        for n in nodes
            JuMP.is_integer(c[n]) ? JuMP.unset_integer(c[n]) : nothing
            JuMP.has_lower_bound(c[n]) ? JuMP.delete_lower_bound(c[n]) : nothing
            JuMP.has_upper_bound(c[n]) ? JuMP.delete_upper_bound(c[n]) : nothing
            JuMP.fix(c[n], n in keys(csol) ? csol[n] : 0)
        end
    end

    ## Constrs for defining e
    @debug "Constrs for defining e"
    @constraint(model, edef[i in 1:num_obs], e[i] == v[i,1] - obs[i,end])

    ## Constrs for defining y
    if occursin("Cozad", formulation)
        @debug "Constrs for defining y (Cozad) $(formulation)"
        @constraints(model, begin
            # [Cozad]-(1d) at most one operator is active
            active_at_most[n1 in nodes], sum(y[n,o] for (n,o) in y_indexes if n == n1) <= 1
            # [Cozad]-(1i) at least one independent variable is active
            active_indvar, sum(y[n,o] for (n,o) in y_indexes if o in indvars) >= 1
            # [Cozad]-(1e) if binary/unary is active, then rch is active
            active_rch[n1 in nleaves], 
                sum(y[n1,o] for o in operBU) <=
                sum(y[n,o] for (n,o) in y_indexes if n == 2*n1+1)
            # [Cozad]-(1f) if binary is active, then lch is active
            active_lch[n1 in nleaves], 
                sum(y[n1,o] for o in operB) <=
                sum(y[n,o] for (n,o) in y_indexes if n == 2*n1)
            # [Cozad]-(1h)
            inactive_rch[n1 in nleaves], 
                sum(y[n1,o] for o in operL) <=
                1 - sum(y[n,o] for (n,o) in y_indexes if n == 2*n1+1)
            # [Cozad]-(1g)
            inactive_lch[n1 in nleaves], 
                sum(y[n1,o] for o in union(operU, operL)) <=
                1 - sum(y[n,o] for (n,o) in y_indexes if n == 2*n1)
        end)
    elseif occursin("New", formulation)
        @debug "Constrs for defining y (New) $(formulation)"
        @constraints(model, begin
            # [Cozad]-(1d) at most one operator is active
            active_at_most[n1 in nodes], sum(y[n,o] for (n,o) in y_indexes if n == n1) <= 1
            # [Cozad]-(1i) at least one independent variable is active
            active_indvar, sum(y[n,o] for (n,o) in y_indexes if o in indvars) >= 1
            # improvement of (1e-1h). Also, improvement of [Newman]'s (18-23)
            active_inactive_rch[n1 in nleaves], 
                sum(y[n1,o] for o in operBU) ==
                sum(y[n,o] for (n,o) in y_indexes if n == 2*n1+1)
            active_inactive_lch[n1 in nleaves], 
                sum(y[n1,o] for o in operB) ==
                sum(y[n,o] for (n,o) in y_indexes if n == 2*n1)
        end)
    else
        nothing
    end

    ## Constrs to remove redundancy
    if occursin("-CR", formulation) || occursin("-NR", formulation)
        # no exp(log(x)) or log(exp(x))
        if issubset(Set("EL"), operA)
            @constraints(model, begin
                redun_explog[n in nleaves; 4*n+3 in nodes], y[n,'E'] + y[2*n+1,'L'] <= 1
                redun_logexp[n in nleaves; 4*n+3 in nodes], y[n,'L'] + y[2*n+1,'E'] <= 1
            end)
        end

        if 'C' in operA && !JuMP.is_integer(c[1])
            # no operation of two constants and no 
            @constraint(model, redun_cst_oper[n in nleaves], y[2*n,'C'] + y[2*n+1,'C'] <= 1)

            # no x-cst, x/cst, and no unary operation with a constant
            if occursin("-CR", formulation)
                @debug "Constrs to remove redundancy (Cozad version) $(formulation)"
                @constraints(model, begin
                    redun_cst_unary[n in nleaves], y[2*n+1,'C'] + sum(y[n,o] for (n,o) in y_indexes if o in operU) <= 1
                    redun_cst_minus[n in nleaves], y[2*n+1,'C'] + sum(y[n,o] for (n,o) in y_indexes if o == '-') <= 1
                    redun_cst_div[n in nleaves], y[2*n+1,'C'] + sum(y[n,o] for (n,o) in y_indexes if o == 'D') <= 1
                end)
            elseif occursin("-NR", formulation)
                @debug "Constrs to remove redundancy (New version) $(formulation)"
                @constraint(model, redun_cst_rhs[n in nleaves], y[2*n+1,'C'] <= y[n,'+'] + y[n,'*'])
            end
        end
    end

    ## Constrs for defining domain bounds
    @debug "Constrs for defining domain bounds"
    if 'C' in operA && cst_abslb > 0
        @constraint(model, domain_cst[n in nodes], c[n]^2 >= cst_abslb^2 * y[n,'C'])
    end

    if ('*' in operA && mult_child_abslb > 0) || ('D' in operA && div_child_abslb > 0)
        @constraint(model, domain_lch_abslb[i in 1:num_obs, n in nleaves], 
        v[i,2*n]^2 >= 
            ( '*' in operA ? mult_child_abslb^2 * y[n,'*'] : 0 ) +
            ( 'D' in operA ? div_child_abslb^2 * y[n,'D'] : 0 ) )

        @constraint(model, domain_rch_abslb[i in 1:num_obs, n in nleaves], 
            v[i,2*n+1]^2 >= 
                ( '*' in operA ? mult_child_abslb^2 * y[n,'*'] : 0 ) +
                ( 'D' in operA ? div_child_abslb^2 * y[n,'D'] : 0 ) )
    end

    if ('R' in operA && sqrt_rhs_lb > 0) || ('E' in operA && exp_rhs_lb > v_lb) || 
        ('L' in operA && log_rhs_lb > v_lb)
        @constraint(model, domain_rhs_lb[i in 1:num_obs, n in nleaves], 
            v[i,2*n+1] >= v_lb +
                ( ('R' in operA && sqrt_rhs_lb > 0) ? (sqrt_rhs_lb - v_lb) * y[n,'D'] : 0 ) + 
                ( ('E' in operA && exp_rhs_lb > v_lb) ? (exp_rhs_lb - v_lb) * y[n,'E'] : 0 ) +
                ( ('L' in operA && log_rhs_lb > v_lb) ? (log_rhs_lb - v_lb) * y[n,'L'] : 0 ) )
    end

    if ('E' in operA && exp_rhs_ub < v_ub) || ('L' in operA && log_rhs_ub < v_ub)
        @constraint(model, domain_rhs_ub[i in 1:num_obs, n in nleaves],
            v[i,2*n+1] <= v_ub +
                ( ('E' in operA && exp_rhs_ub < v_ub) ? (exp_rhs_ub - v_ub) * y[n,'E'] : 0 ) +
                ( ('L' in operA && log_rhs_ub < v_ub) ? (log_rhs_ub - v_ub) * y[n,'L'] : 0 ) )
    end


    ## Constrs for defining v
    if occursin("Cozad", formulation) || occursin("New", formulation)
        if occursin("Cozad", formulation)
            @debug "Constrs for defining v (Cozad) $(formulation)"

            @constraints(model, begin
                v_inactive_ub[i in 1:num_obs, n in nodes], v[i,n] <= v_ub * sum(y[n1,o] for (n1,o) in y_indexes if n1 == n)
                v_inactive_lb[i in 1:num_obs, n in nodes], v[i,n] >= v_lb * sum(y[n1,o] for (n1,o) in y_indexes if n1 == n)
            end)

            @constraints(model, begin
                v_indvars_ub[i in 1:num_obs, n in nodes, o in indvars], v[i,n] - obs[i,o] <= (v_ub - obs[i,o]) * (1 - y[n,o])
                v_indvars_lb[i in 1:num_obs, n in nodes, o in indvars], v[i,n] - obs[i,o] >= (v_lb - obs[i,o]) * (1 - y[n,o])
            end)

            if 'P' in operA
                @constraints(model, begin
                    v_pi_ub[i in 1:num_obs, n in nodes], v[i,n] - pi <= (v_ub - v_lb) * (1 - y[n,'P'])
                    v_pi_lb[i in 1:num_obs, n in nodes], v[i,n] - pi >= (v_lb - v_ub) * (1 - y[n,'P'])
                end)
            end

            if 'C' in operA
                @constraints(model, begin
                    v_cst_ub[i in 1:num_obs, n in nodes], v[i,n] - c[n] <= (v_ub - c_lb) * (1 - y[n,'C'])
                    v_cst_lb[i in 1:num_obs, n in nodes], v[i,n] - c[n] >= (v_lb - c_ub) * (1 - y[n,'C'])
                end)
            end

        elseif occursin("New", formulation)
            @debug "Constrs for defining v (New) $(formulation)"

            @constraints(model, begin
                v_indvars_ub[i in 1:num_obs, n1 in nodes],
                    v[i,n1] <= 
                        v_ub * sum(y[n,o] for (n,o) in y_indexes if n==n1 && o in operB) +
                        min(v_ub, sqrt(v_ub)) * ((n1,'R') in y_indexes ? y[n1,'R'] : 0) +
                        min(v_ub, exp(v_ub)) * ((n1,'E') in y_indexes ? y[n1,'E'] : 0) +
                        min(v_ub, log(v_ub)) * ((n1,'L') in y_indexes ? y[n1,'L'] : 0) +
                        c_ub * ((n1,'C') in y_indexes ? y[n1,'C'] : 0) +
                        pi * ((n1,'P') in y_indexes ? y[n1,'P'] : 0) +
                        sum(obs[i,o] * y[n1,o] for o in indvars)
                v_indvars_lb[i in 1:num_obs, n1 in nodes],
                    v[i,n1] >= 
                        v_lb * sum(y[n,o] for (n,o) in y_indexes if n==n1 && o in operB) +
                        max(v_lb, sqrt(sqrt_rhs_lb)) * ((n1,'R') in y_indexes ? y[n1,'R'] : 0) +
                        max(v_lb, exp(exp_rhs_lb)) * ((n1,'E') in y_indexes ? y[n1,'E'] : 0) +
                        max(v_lb, log(log_rhs_lb)) * ((n1,'L') in y_indexes ? y[n1,'L'] : 0) +
                        c_lb * ((n1,'C') in y_indexes ? y[n1,'C'] : 0) +
                        pi * ((n1,'P') in y_indexes ? y[n1,'P'] : 0) +
                        sum(obs[i,o] * y[n1,o] for o in indvars)
            end)

            if 'C' in operA
                @constraints(model, begin
                    v_cst_ub[i in 1:num_obs, n in nodes],
                        v[i,n] - c[n] <= (v_ub - c_lb) * (1 - y[n,'C'])
                    v_cst_lb[i in 1:num_obs, n in nodes],
                        v[i,n] - c[n] >= (v_lb - c_ub) * (1 - y[n,'C'])
                end)
            end
        end


        if '+' in operA
            @constraints(model, begin
                v_plus_ub[i in 1:num_obs, n in nleaves],
                    v[i,n] - (v[i,2*n] + v[i,2*n+1]) <= (v_ub - 2 * v_lb) * (1 - y[n,'+'])
                v_plus_lb[i in 1:num_obs, n in nleaves],
                    v[i,n] - (v[i,2*n] + v[i,2*n+1]) >= (v_lb - 2 * v_ub) * (1 - y[n,'+'])
            end)
        end

        if '-' in operA
            @constraints(model, begin
                v_minus_ub[i in 1:num_obs, n in nleaves],
                    v[i,n] - (v[i,2*n] - v[i,2*n+1]) <= (2 * v_ub - v_lb) * (1 - y[n,'-'])
                v_minus_lb[i in 1:num_obs, n in nleaves],
                    v[i,n] - (v[i,2*n] - v[i,2*n+1]) >= (2 * v_lb - v_ub) * (1 - y[n,'-'])
            end)
        end

        if '*' in operA
            @constraints(model, begin
                v_mult_ub[i in 1:num_obs, n in nleaves],
                    v[i,n] - (v[i,2*n] * v[i,2*n+1]) <= 
                    (v_ub - min(v_lb^2, v_ub^2, v_lb*v_ub)) * (1 - y[n,'*'])
                v_mult_lb[i in 1:num_obs, n in nleaves],
                    v[i,n] - (v[i,2*n] * v[i,2*n+1]) >= 
                    (v_lb - max(v_lb^2, v_ub^2)) * (1 - y[n,'*'])
            end)
        end

        if 'D' in operA
            @constraints(model, begin
                v_div_ub[i in 1:num_obs, n in nleaves],
                    v[i,n] * v[i,2*n+1] - v[i,2*n] <= 
                    (max(v_lb^2, v_ub^2) - v_lb) * (1 - y[n,'D'])
                v_div_lb[i in 1:num_obs, n in nleaves],
                    v[i,n] * v[i,2*n+1] - v[i,2*n] >= 
                    (min(v_lb^2, v_ub^2, v_lb*v_ub) - v_ub) * (1 - y[n,'D'])
            end)
        end

        if 'R' in operA
            @constraints(model, begin
                v_sqrt_ub[i in 1:num_obs, n in nleaves],
                    v[i,n] * v[i,n] - v[i,2*n+1] <= 
                    (max(v_lb^2, v_ub^2) - v_lb) * (1 - y[n,'R'])
                v_sqrt_lb[i in 1:num_obs, n in nleaves],
                    v[i,n] * v[i,n] - v[i,2*n+1] >= 
                    (- v_ub) * (1 - y[n,'R'])
            end)
        end

        if 'E' in operA
            @NLconstraints(model, begin
                v_exp_ub[i in 1:num_obs, n in nleaves],
                    v[i,n] - exp(v[i,2*n+1]) <= v_ub * (1 - y[n,'E'])
                v_exp_lb[i in 1:num_obs, n in nleaves],
                    v[i,n] - exp(v[i,2*n+1]) >= (v_lb - exp(v_ub)) * (1 - y[n,'E'])
            end)
        end

        if 'L' in operA
            @NLconstraints(model, begin
                v_log_ub[i in 1:num_obs, n in nleaves],
                    exp(v[i,n]) - v[i,2*n+1] <= (exp(v_ub) - v_lb) * (1 - y[n,'L'])
                v_log_lb[i in 1:num_obs, n in nleaves],
                    exp(v[i,n]) - v[i,2*n+1] >= (-v_ub) * (1 - y[n,'L'])
            end)
        end
    end   

    @debug "Optimize!"
    optimize!(model)
    time = JuMP.solve_time(model)

    ## Find the index of the best feasible solution
    ## We have this code because the returned solution is sometimes infeasible
    num_results = JuMP.result_count(model)
    if num_results == 0 
        # There is no solution. We stop here.
        @info   "Print result (no solution)\n" *
                @sprintf("time                  = %12.2f\n",    time) *
                @sprintf("termination status    = %s\n",        termination_status(model)) *
                @sprintf("num_results           = %12d\n",      num_results) *
                @sprintf("bound                 = %12.6f\n",    JuMP.objective_bound(model)) *
                @sprintf("relgap                = %12.6f\n",    JuMP.relative_gap(model)) *
                @sprintf("bnbnode               = %12d\n",      JuMP.node_count(model))
        return false, false, time, 1e+09, nothing, nothing, nothing
    end

    arr_obj     = 1e+09 * ones(num_results)
    arr_active  = zeros(num_results)
    for i in num_results:-1:1
        @assert JuMP.has_values(model; result = i)
 
        ysol = get_ysol(y, y_indexes, nodes, operA; num_result=i)
        csol = get_csol(c, nodes; num_result=i)
        vsol, err = get_vsol(ysol, csol, obs)
        arr_obj[i] = err ? 1e+09 : compute_mse(vsol[:,1], obs[:,end])
        arr_active[i] = err ? 0 : length(ysol)
    end

    ## Now we know the best solution. Print and return it.
    # b: the index of the best solution
    b = argmin(arr_obj)

    # Print the best solution
    ysol = get_ysol(y, y_indexes, nodes, operA; num_result=b)
    csol = get_csol(c, nodes; num_result=b)
    vsol, err = get_vsol(ysol, csol, obs)
    @info "Print the best solution\n" * print_tree(get_treesol(ysol, csol))
    @info "Print errors" * print_error(vsol, obs) 
    @debug "ysol" ysol
    @debug "csol" csol
    @debug "treesol" get_treesol(ysol, csol)

    # Print the current solution
    @debug print_y(y, y_indexes, nodes, operA)
    'C' in operA && @debug print_c(c, y, nodes)
    @debug print_v(v, nodes, 1:num_obs)
    @debug print_obs(obs)
    @debug print_num_oper(num_oper, operA)

    # Print summary
    @info   "Print result\n" *
            @sprintf("time                  = %12.2f\n",    time) *
            @sprintf("termination_status    = %s\n",        termination_status(model)) *
            @sprintf("num_results           = %12d\n",      num_results) *
            @sprintf("best_index            = %12d\n",      b) *
            @sprintf("obj (by solver)       = %12.6e\n",    JuMP.objective_value(model)) *
            @sprintf("obj (recomputed)      = %12.6e\n",    arr_obj[b]) *
            @sprintf("objbound              = %12.6e\n",    JuMP.objective_bound(model)) *
            @sprintf("relative_gap          = %12.6f\n",    JuMP.relative_gap(model)) *
            @sprintf("bnbnode_count         = %12d\n",      JuMP.node_count(model)) *
            @sprintf("active_count          = %12d\n",      arr_active[b])
    
    @info "Final result of MINLP()\n" * print_final_table(arr_obj, zeros(num_results), arr_active)

    if print_all_solutions
        for i in 1:num_results
            if arr_active[i] > 0
                @assert JuMP.has_values(model; result = i)
        
                ysol1 = get_ysol(y, y_indexes, nodes, operA; num_result=i)
                csol1 = get_csol(c, nodes; num_result=i)
                @info   "Print solution $i\n" *
                        @sprintf("obj (recomputed)      = %12.6f\n",    arr_obj[i]) *
                        @sprintf("active_count          = %12d\n",      arr_active[i]) *
                        print_tree(get_treesol(ysol1, csol1))   
            end
        end
    end

    relerr_recompute = abs(JuMP.objective_value(model)/arr_obj[b] - 1)

    return true, relerr_recompute < 1e-03, time, arr_obj[b], ysol, csol, vsol
    # return false, nothing, nothing, nothing, nothing, nothing, nothing
end