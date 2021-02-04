function get_fixed_nodes(ysol, fixlevel)
    if fixlevel < 0
        return nothing
    elseif fixlevel == 0
        return keys(ysol)
    else
        fixed_nodes = keys(ysol)
        for i in 1:fixlevel
            fixed_nodes = get_nonleaves(fixed_nodes)
        end
        return fixed_nodes
    end
end


function reset_variables(model, df_vars, bounds; unused_nodes=nothing, fixed_nodes=nothing)
    """

    """
    variables = JuMP.all_variables(model)

    # ##### (1) reset variables to original setting
    for row in eachrow(df_vars)
        var = variables[row.index]
        if JuMP.is_fixed(var)
            JuMP.unfix(var)
            if row.vartype in ["c"]
                JuMP.set_lower_bound(var, bounds["c_lb"])
                JuMP.set_lower_bound(var, bounds["c_ub"])
            end
        end
    end

    if isnothing(unused_nodes)
        unused_nodes = []
    end

    if isnothing(fixed_nodes)
        fixed_nodes = []
    end

    for row in eachrow(df_vars)
        var = variables[row.index]
        if row.vartype in ["y", "c"]
            if row.nodeid in unused_nodes
                JuMP.fix(var, 0.0, force=true)
            elseif row.nodeid in fixed_nodes
                JuMP.fix(var, row.sol, force=true)
            else
                if row.vartype in ["y"] && row.sol > 0.5
                    JuMP.set_start_value(var, 1.0)
                elseif row.vartype in ["y"] && row.sol <= 0.5
                    JuMP.set_start_value(var, 0.0)
                else
                    JuMP.set_start_value(var, row.sol)
                end
            end
        end
    end

    model
end

# function reset_variables_by_maxdist_fixlevel(model, df_vars, bounds, ysol, maxdist, fixlevel)
#     fixed_nodes = get_fixed_nodes(ysol, fixlevel)
#     nodes_ground = Set(1:maximum(df_vars.nodeid))
#     unused_nodes = get_unused_nodes(ysol, nodes_ground, maxdist)
#     reset_variables(model, df_vars, bounds; unused_nodes=nothing, fixed_nodes=fixed_nodes)
# end

function add_neighbor_constraints(model, df_vars, mindist, maxdist)
    variables = JuMP.all_variables(model)
    df_y_1 = filter(row -> row.vartype == "y" && row.sol > 0.5, df_vars)
    df_y_0 = filter(row -> row.vartype == "y" && row.sol <= 0.5, df_vars)
    active_nodes = unique(df_y_1.nodeid)
    @constraint(model, mindist <= sum(variables[row.index] for row in eachrow(df_y_0) if !(row.nodeid in active_nodes)) 
                + sum(1 - variables[i] for i in df_y_1.index) <= maxdist)

    return model
end

function set_restricted_MINLP(model, df_vars, bounds, ysol, mindist, maxdist, fixlevel)
    fixed_nodes = get_fixed_nodes(ysol, fixlevel)
    nodes_ground = Set(1:maximum(df_vars.nodeid))
    unused_nodes = get_unused_nodes(ysol, nodes_ground, maxdist)
        model = reset_variables(model, df_vars, bounds; unused_nodes=unused_nodes, fixed_nodes=fixed_nodes)
    model = add_neighbor_constraints(model, df_vars, mindist, maxdist)
    model
end

function initialize_STreCH_parameters(params)
    mindist = 0
    maxdist = params["INT_SOLU_DIST_INIT"]
    fixlevel = params["INT_FIXLEVEL_INIT"]
    nodelimit = params["INT_NODELIMIT_INIT"]
    mindist, maxdist, fixlevel, nodelimit
end

function update_STreCH_parameters(params, mindist, maxdist, fixlevel, nodelimit)
    if fixlevel < params["INT_FIXLEVEL_MAX"]
        fixlevel += 1
    else
        if maxdist < params["INT_SOLU_DIST_MAX"]
            fixlevel = params["INT_FIXLEVEL_INIT"]
            mindist = maxdist + 1
            maxdist = maxdist + params["INT_SOLU_DIST_INCR"]
        else
            if nodelimit < params["INT_NODELIMIT_MAX"]
                fixlevel = params["INT_FIXLEVEL_INIT"]
                mindist = 0
                maxdist = params["INT_SOLU_DIST_INIT"]
                nodelimit = 10 * nodelimit
            else
                return nothing, nothing, nothing, nothing
            end
        end
    end

    mindist, maxdist, fixlevel, nodelimit
end

function update_df_sol(df_sol, obs, df_vars, num_iter, time)
    ysol, csol = get_ysol_csol(df_vars)
    formula = get_formula(ysol, csol)[1]
    objval = compute_err_formula(obs, formula, "mse"; err2Inf=true)
    append!(df_sol, Dict("index"=>num_iter, "formula"=>formula, "objval"=>objval, "treesize"=>length(ysol), "time"=>time))
    df_sol
end

function STreCH(dir; datafile="data.txt", paramfile="params.txt")
    EPSILON = 1e-12
    params = read_parameters(dir * paramfile)
    df_sol = DataFrame(index=Int[], formula=String[], objval=Float64[], treesize=Int[])

    io = open_logger(dir * "log.log", params["LOGLEVEL"])

    t_begin = time()
    num_iter = 1

    ### Parameters
    nodes = OrderedSet(1:(2^(params["INT_MAXDEPTH"]+1)-1))
    obs = DelimitedFiles.readdlm(dir * datafile)
    
    ### Solve an initial problem
    model, df_vars, bounds = get_MINLP_model(nodes, obs, params["OPERATORS"], params["FORMULATION"], params)
    model = reset_variables(model, df_vars, bounds, unused_nodes=8:1000)
    set_optimizer(model, SCIP.Optimizer)

    ### Set absgap, timelimit
    remaining_time = params["DBL_TIMELIMIT"] - (time() - t_begin)
    set_optimizer_attribute(model, "limits/time", remaining_time)
    set_optimizer_attribute(model, "display/verblevel", 2)

    ### Solve
    optimize!(model)
    
    ### Save the result
    objval_best, ysol_best, csol_best, df_vars_best = save_solution(dir, model, df_vars, obs, save_df_sol=false)
    df_sol = update_df_sol(df_sol, obs, df_vars_best, num_iter, time() - t_begin)
    CSV.write(dir * "df_sol.csv", df_sol)

    ### Improve
    mindist, maxdist, fixlevel, nodelimit = initialize_STreCH_parameters(params)
    while remaining_time > 0 && objval_best > EPSILON
        num_iter += 1
        @info "num_iter = $num_iter, fixlevel = $fixlevel, mindist=$mindist, maxdist=$maxdist, nodelimit=$nodelimit"
        model, df_vars, bounds = get_MINLP_model(nodes, obs, params["OPERATORS"], params["FORMULATION"], params)
        model = set_restricted_MINLP(model, df_vars_best, bounds, ysol_best, mindist, maxdist, fixlevel)
        
        ### Set the optimizer and its attributes
        set_optimizer(model, SCIP.Optimizer)
        remaining_time = params["DBL_TIMELIMIT"] - (time() - t_begin)
        absgap = objval_best < 1e-06 ? 0.0 : 0.99 * objval_best
        set_optimizer_attribute(model, "limits/time", remaining_time)
        set_optimizer_attribute(model, "limits/nodes", nodelimit)
        set_optimizer_attribute(model, "limits/absgap", absgap)
        set_optimizer_attribute(model, "display/verblevel", 2)
        
        ### Solve
        optimize!(model)
        termination_status = JuMP.termination_status(model)
        @info "num_iter = $num_iter , $termination_status"

        ### Read solution
        objval_new, ysol_new, csol_new, df_vars_new = save_solution(dir, model, df_vars, obs, save_df_sol=false)
        df_sol = update_df_sol(df_sol, obs, df_vars_new, num_iter, time() - t_begin)
        CSV.write(dir * "df_sol.csv", df_sol)

        ### Compare with the best solution and update paramters for the next iteration
        if objval_new < objval_best - EPSILON
            objval_best = objval_new
            ysol_best = ysol_new
            csol_best = csol_new
            df_vars_best = df_vars_new
            mindist, maxdist, fixlevel, nodelimit = initialize_STreCH_parameters(params)
        else
            """
            update parameters and keep the current solution
            """
            mindist, maxdist, fixlevel, nodelimit = update_STreCH_parameters(params, mindist, maxdist, fixlevel, nodelimit)
            if isnothing(mindist)
                break
            end
        end
        
        remaining_time = params["DBL_TIMELIMIT"] - (time() - t_begin)
    end

    # ### Remove files
    # Glob.glob(dir * "*.csv")
    # for file in Glob.glob(dir * "*.csv")
    #     if !occursin("df_sol.csv", file) && !occursin("df_vars.csv", file)
    #         rm(file)
    #     end
    # end 

    close(io)
end