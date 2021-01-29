function get_vbounds(obs, operators)
    # Initialization
    operators = Set(operators)
    (v_lb, v_ub) = (-1e+09, 1e+09)
    factor1 = 10
    factor2 = 2
    eps = 1e-09

    # if !isempty(intersect(operators, Set("EL")))
    #     factor1 = 5
    #     factor2 = 2
    # end

    ## Compute min/max of indep/dep variables
    # (z_lb, z_ub) = (minimum(obs[:,end]), maximum(obs[:,end]))
    # (x_lb, x_ub) = (minimum(obs[:,1:end-1]), maximum(obs[:,1:end-1]))
    (obs_lb, obs_ub) = (minimum(obs), maximum(obs))

    ## (1) Shrink bounds
    (v_lb, v_ub) = (maximum([v_lb, factor1^(-sign(obs_lb)) * obs_lb]), minimum([v_ub, factor1 * obs_ub]))

    ## (2) Expand bounds if it is too tight
    (v_lb, v_ub) = (minimum([v_lb, factor2^(-sign(obs_lb)) * obs_lb]), maximum([v_ub, factor2 * obs_ub]))

    ## (3) make bounds symmetric and +- epsilon
    magnitude = maximum([abs(v_lb), abs(v_ub)]) + eps
    (v_lb, v_ub) = (-magnitude, magnitude)

    v_lb, v_ub
end

function get_bounds(obs, operators)
    bounds = Dict()
    v_lb, v_ub = get_vbounds(obs, operators)
    merge!(bounds, Dict("v_lb"=>v_lb, "v_ub"=>v_ub))
    merge!(bounds, Dict("c_lb"=>-2.0, "c_ub"=>2.0))

    bounds["cst_abslb"] = 1/(4*pi)
    bounds["mult_child_abslb"] = 0.01
    bounds["div_child_abslb"] = 0.01
    bounds["sqrt_rhs_lb"] = 0.01
    bounds["exp_rhs_lb"] = v_lb > 0 ? max(v_lb, log(v_lb)) : v_lb
    bounds["exp_rhs_ub"] = min(v_ub, log(v_ub))
    bounds["log_rhs_lb"] = max(v_lb, exp(v_lb))
    bounds["log_rhs_ub"] = min(v_ub, exp(v_ub))

    bounds
end

function get_df_vars(nodes, y_indexes, num_obs)
    df_vars = DataFrame(vartype=String[], nodeid=Int[], dataid=Int[], operator=Any[], sol=Float64[])
    for (n,o) in y_indexes
        append!(df_vars, Dict("vartype"=>"y", "nodeid"=>n, "dataid"=>0, "operator"=>o, "sol"=>0.0))
    end

    for n in nodes, i=1:num_obs 
        append!(df_vars, Dict("vartype"=>"v", "nodeid"=>n, "dataid"=>i, "operator"=>0, "sol"=>0.0))
    end

    for n in nodes
        append!(df_vars, Dict("vartype"=>"c", "nodeid"=>n, "dataid"=>0, "operator"=>0, "sol"=>0.0))
    end

    for i=1:num_obs
        append!(df_vars, Dict("vartype"=>"e", "nodeid"=>0, "dataid"=>i, "operator"=>0, "sol"=>0.0))
    end

    df_vars.index = 1:nrow(df_vars)
    select!(df_vars, ["index", "vartype", "nodeid", "dataid", "operator", "sol"])
    df_vars
end

function create_decision_variables(model, nodes, operA, y_indexes, num_obs, bounds)
    @variable(model, y[n in nodes, o in operA; (n,o) in y_indexes], Bin)
    @variable(model, bounds["v_lb"] <= v[1:num_obs, nodes] <= bounds["v_ub"])
    @variable(model, bounds["c_lb"] <= c[nodes] <= bounds["c_ub"])
    @variable(model, -1e+09 <= e[1:num_obs] <= 1e+09)
    df_vars = get_df_vars(nodes, y_indexes, num_obs)
    model, y, c, v, e, df_vars
end

function set_objective_function(model, e, num_obs)
    @objective(model, Min, 1 / num_obs * sum(e.^2))
    model
end

function add_e_deining_constraints(model, v, e, obs, num_obs)
    @constraint(model, edef[i in 1:num_obs], e[i] == v[i,1] - obs[i,end])
    model
end

function add_tree_defining_constraints(model, formulation, y, y_indexes, nodes, nleaves, operBU, operB, operU, operL, indvars)
    if occursin("Coz", formulation)
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
    elseif occursin("Imp", formulation)
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

    model
end

function add_value_defining_constraints(model, formulation, y, c, v, y_indexes, bounds, nodes, nleaves, obs, num_obs, operA, operB, indvars)
    ### Constrs for defining domain bounds
    @debug "Constrs for defining domain bounds"
    
    if 'C' in operA && bounds["cst_abslb"] > 0
        @constraint(model, domain_cst[n in nodes], c[n]^2 >= bounds["cst_abslb"]^2 * y[n,'C'])
    end

    if ('*' in operA && bounds["mult_child_abslb"] > 0) || ('D' in operA && bounds["div_child_abslb"] > 0)
        @constraint(model, domain_lch_abslb[i in 1:num_obs, n in nleaves], 
        v[i,2*n]^2 >= 
            ( '*' in operA ? bounds["mult_child_abslb"]^2 * y[n,'*'] : 0 ) +
            ( 'D' in operA ? bounds["div_child_abslb"]^2 * y[n,'D'] : 0 ) )

        @constraint(model, domain_rch_abslb[i in 1:num_obs, n in nleaves], 
            v[i,2*n+1]^2 >= 
                ( '*' in operA ? bounds["mult_child_abslb"]^2 * y[n,'*'] : 0 ) +
                ( 'D' in operA ? bounds["div_child_abslb"]^2 * y[n,'D'] : 0 ) )
    end

    if ('R' in operA && bounds["sqrt_rhs_lb"] > 0) || ('E' in operA && bounds["exp_rhs_lb"] > v_lb) || 
        ('L' in operA && bounds["log_rhs_lb"] > v_lb)
        @constraint(model, domain_rhs_lb[i in 1:num_obs, n in nleaves], 
            v[i,2*n+1] >= v_lb +
                ( ('R' in operA && bounds["sqrt_rhs_lb"] > 0) ? (bounds["sqrt_rhs_lb"] - v_lb) * y[n,'D'] : 0 ) + 
                ( ('E' in operA && bounds["exp_rhs_lb"] > v_lb) ? (bounds["exp_rhs_lb"] - v_lb) * y[n,'E'] : 0 ) +
                ( ('L' in operA && bounds["log_rhs_lb"] > v_lb) ? (bounds["log_rhs_lb"] - v_lb) * y[n,'L'] : 0 ) )
    end

    if ('E' in operA && bounds["exp_rhs_ub"] < v_ub) || ('L' in operA && bounds["log_rhs_ub"] < v_ub)
        @constraint(model, domain_rhs_ub[i in 1:num_obs, n in nleaves],
            v[i,2*n+1] <= v_ub +
                ( ('E' in operA && bounds["exp_rhs_ub"] < v_ub) ? (bounds["exp_rhs_ub"] - v_ub) * y[n,'E'] : 0 ) +
                ( ('L' in operA && bounds["log_rhs_ub"] < v_ub) ? (bounds["log_rhs_ub"] - v_ub) * y[n,'L'] : 0 ) )
    end

    if 'S' in operA
        @constraints(model, begin
            domain_asincos_node_lb[i in 1:num_obs, n in nleaves], v[i,n] >= v_lb*(1 - y[n,'S'])
            domain_asincos_node_ub[i in 1:num_obs, n in nleaves], v[i,n] <= y[n,'S'] + v_ub*(1 - y[n,'S'])
            domain_asincos_rch_lb[i in 1:num_obs, n in nleaves], v[i,2*n+1] >= -y[n,'S'] + v_lb*(1 - y[n,'S'])
            domain_asincos_rch_ub[i in 1:num_obs, n in nleaves], v[i,2*n+1] <= y[n,'S'] + v_ub*(1 - y[n,'S'])
        end)
    end

    ### Constrs for defining v
    v_lb = bounds["v_lb"]
    v_ub = bounds["v_ub"]
    c_lb = bounds["c_lb"]
    c_ub = bounds["c_ub"]
    if occursin("Coz", formulation)
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
    elseif occursin("Imp", formulation)
        @debug "Constrs for defining v (New) $(formulation)"

        @constraints(model, begin
            v_indvars_ub[i in 1:num_obs, n1 in nodes],
                v[i,n1] <= 
                    v_ub * sum(y[n,o] for (n,o) in y_indexes if n==n1 && o in operB) +
                    min(v_ub, sqrt(v_ub)) * ((n1,'R') in y_indexes ? y[n1,'R'] : 0) +
                    1 * ((n1,'S') in y_indexes ? y[n1,'S'] : 0) +
                    min(v_ub, exp(v_ub)) * ((n1,'E') in y_indexes ? y[n1,'E'] : 0) +
                    min(v_ub, log(v_ub)) * ((n1,'L') in y_indexes ? y[n1,'L'] : 0) +
                    c_ub * ((n1,'C') in y_indexes ? y[n1,'C'] : 0) +
                    pi * ((n1,'P') in y_indexes ? y[n1,'P'] : 0) +
                    sum(obs[i,o] * y[n1,o] for o in indvars)
            v_indvars_lb[i in 1:num_obs, n1 in nodes],
                v[i,n1] >= 
                    v_lb * sum(y[n,o] for (n,o) in y_indexes if n==n1 && o in operB) +
                    max(v_lb, sqrt(bounds["sqrt_rhs_lb"])) * ((n1,'R') in y_indexes ? y[n1,'R'] : 0) +
                    max(v_lb, exp(bounds["exp_rhs_lb"])) * ((n1,'E') in y_indexes ? y[n1,'E'] : 0) +
                    max(v_lb, log(bounds["log_rhs_lb"])) * ((n1,'L') in y_indexes ? y[n1,'L'] : 0) +
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

    if 'S' in operA
        @constraints(model, begin
            v_asincos_ub[i in 1:num_obs, n in nleaves],
                v[i,n] * v[i,n] + v[i,2*n+1] * v[i,2*n+1] <= y[n,'S'] + 2 * max(v_lb^2, v_ub^2) * (1 - y[n,'S'])
            v_asincos_lb[i in 1:num_obs, n in nleaves],
                v[i,n] * v[i,n] + v[i,2*n+1] * v[i,2*n+1] >= y[n,'S']
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

    model 
end

function add_redundancy_eliminating_constraints(model, formulation, y, c, y_indexes, nodes, nleaves, operA, operU)
    ### Nested Operations, e.g., no exp(log(x)) or log(exp(x))
    if issubset(Set("EL"), operA)
        @constraints(model, begin
            redun_explog[n in nleaves; 4*n+3 in nodes], y[n,'E'] + y[2*n+1,'L'] <= 1
            redun_logexp[n in nleaves; 4*n+3 in nodes], y[n,'L'] + y[2*n+1,'E'] <= 1
        end)
    end

    ### Constant Operations
    if 'C' in operA && !JuMP.is_integer(c[1])
        # no operation of two constants and no 
        @constraint(model, redun_cst_oper[n in nleaves], y[2*n,'C'] + y[2*n+1,'C'] <= 1)

        # no x-cst, x/cst, and no unary operation with a constant
        if occursin("Coz", formulation)
            @debug "Constrs to remove redundancy (Cozad version) $(formulation)"
            @constraints(model, begin
                redun_cst_unary[n in nleaves], y[2*n+1,'C'] + sum(y[nn,o] for (nn,o) in y_indexes if nn == n && o in operU) <= 1
                redun_cst_minus[n in nleaves], y[2*n+1,'C'] + sum(y[nn,o] for (nn,o) in y_indexes if nn == n && o == '-') <= 1
                redun_cst_div[n in nleaves], y[2*n+1,'C'] + sum(y[nn,o] for (nn,o) in y_indexes if nn == n && o == 'D') <= 1
            end)
        elseif occursin("Imp", formulation)
            @debug "Constrs to remove redundancy (New version) $(formulation)"
            oper = intersect(operA, "+*")
            @constraint(model, redun_cst_rhs[n in nleaves], y[2*n+1,'C'] <= sum(y[n,o] for o in intersect(operA, "+*")))
        end
    end

    ### Association Properties
    if occursin("Imp", formulation)
        @debug "Constrs to remove redundancy (New version 2) $(formulation)"
        nodes_constr = intersect(get_nodes_complete(nodes), get_nodes_grandparents(nodes))
        if issubset("+-", operA)
            @constraint(model, redun_pm[n in nodes_constr], y[n,'+'] + y[2*n+1,'-'] <= 1)
        end
        if issubset("*D", operA)
            @constraint(model, redun_pd[n in nodes_constr], y[n,'*'] + y[2*n+1,'D'] <= 1)
        end
    end

    model
end

function add_implication_cuts(model)
    model
end

function add_symmetry_breaking_constraints(model, y, v, nodes, nleaves, bounds)
    @debug "Constrs to break symmetry"
    nodes_constr = intersect(get_nodes_complete(nodes), nleaves)
    @constraint(model, sym[n in nodes_constr], v[1,2*n] - v[1,2*n+1] >= (bounds["v_lb"] - bounds["v_ub"]) * (1 - y[n,'+'] - y[n,'*']))
    model
end

function get_MINLP_model(nodes, obs, operators, formulation)
    ### Declare useful parameters 
    leaves		= get_leaves(nodes)
    nleaves		= get_nonleaves(nodes)   
    num_obs 	= size(obs, 1)
    num_indvars = size(obs, 2) - 1
    indvars		= OrderedSet(1:num_indvars)
    operA		= union(operators, indvars)
    operB       = intersect(operators, Set("+-*D"))
    operU       = intersect(operators, Set("RSEL"))
    operBU      = union(operB, operU)
    operL       = union(intersect(operators, Set("CP")), indvars)
    y_indexes   = [(n,o) for n in nodes for o in operA if n in nleaves || o in operL]

    ### Create a model
    model = Model() 
    bounds = get_bounds(obs, operators)
    model, y, c, v, e, df_vars = create_decision_variables(model, nodes, operA, y_indexes, num_obs, bounds)
    model = set_objective_function(model, e, num_obs)
    model = add_e_deining_constraints(model, v, e, obs, num_obs)
    model = add_tree_defining_constraints(model, formulation, y, y_indexes, nodes, nleaves, operBU, operB, operU, operL, indvars)
    model = add_value_defining_constraints(model, formulation, y, c, v, y_indexes, bounds, nodes, nleaves, obs, num_obs, operA, operB, indvars)
    if occursin("-R", formulation) || occursin("-F", formulation) 
        model = add_redundancy_eliminating_constraints(model, formulation, y, c, y_indexes, nodes, nleaves, operA, operU)
    end
    if occursin("-S", formulation) || occursin("-F", formulation) 
        model = add_symmetry_breaking_constraints(model, y, v, nodes, nleaves, bounds)
    end
    # if occursin("Imp", formulation)
    #     model = add_implication_cuts(model)
    # end

    model, df_vars
end

function save_solution(dir, model, df_vars, obs)
    num_results = JuMP.result_count(model)
    df_sol = DataFrame(index=Int[], formula=String[], objval=Float64[], treesize=Int[])
    if num_results > 0
        for i in 1:num_results
            if JuMP.has_values(model; result = i)
                df_vars.sol = JuMP.value.(JuMP.all_variables(model), result=i)
                ysol, csol = get_ysol_csol(df_vars)
                @debug "ysol and csol", ysol, csol
                formula = get_formula(ysol, csol)[1]
                objval = compute_err_formula(obs, formula, "mse"; err2Inf=true)
                append!(df_sol, Dict("index"=>i, "formula"=>formula, "objval"=>objval, "treesize"=>length(ysol)))
            end
        end
        df_sol.time = repeat([JuMP.solve_time(model)], nrow(df_sol))
        df_sol.bnbnodes = repeat([JuMP.node_count(model)], nrow(df_sol))
    end

    # CSV.write(dir * "df_vars.csv", df_vars)
    CSV.write(dir * "df_sol.csv", df_sol)
end

function MINLP(dir, datafile, operators, maxdepth, timelimit, formulation; loglevel="Info")
    io = open(dir * "log.log", "w+")
    if loglevel == "Debug"
        logger = SimpleLogger(io, Logging.Debug)
    else
        logger = SimpleLogger(io, Logging.Info)
    end
    global_logger(logger)

    ### Parameters
    nodes = OrderedSet(1:(2^(maxdepth+1)-1))
    obs = DelimitedFiles.readdlm(dir * datafile)
    
    ### Build a MINLP model
    model, df_vars = get_MINLP_model(nodes, obs, operators, formulation)
    
    ### Set the MINLP solver
    set_optimizer(model, SCIP.Optimizer)
    set_optimizer_attribute(model, "limits/time", timelimit)

    ### Solve
    optimize!(model)
    
    ### Save the result
    save_solution(dir, model, df_vars, obs)

    close(io)
end