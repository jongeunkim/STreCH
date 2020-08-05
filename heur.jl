include("symbolic.jl")

function update_nodes!(nodes, y, y_indexes, max_depth)
    nodes_current = Set(nodes)
    
    # Expand
    for (n,o) in y_indexes
        if JuMP.value(y[n,o]) > 0.5
            union!(nodes, Set([2*n, 2*n+1]))
        end
    end
    intersect!(nodes, Set(1:(2^(max_depth+1)-1)))

    # Expand if the new set is a subset of the current
    while issubset(nodes, nodes_current)
        for n in nodes
            union!(nodes, Set([2*n, 2*n+1]))
        end
    end
    intersect!(nodes, Set(1:(2^(max_depth+1)-1)))

    OrderedSet(sort(collect(nodes)))
end

function update_num_oper_lb(y, y_indexes)
    num_oper_lb = Dict{Any, Int}()

    for (n,o) in y_indexes
        if JuMP.value(y[n,o]) > 0.5 && o != 'C'
            num_oper_lb[o] = o in keys(num_oper_lb) ? num_oper_lb[o] + 1 : 1
        end
    end

    num_oper_lb
end

function update_num_oper_ub(max_active, num_oper_lb, stepsize)
    num_oper_ub = Dict{Any, Int}()

    for (o, v) in num_oper_lb
        num_oper_ub[o] = num_oper_lb[o] + stepsize
    end

    num_oper_ub['T'] = max_active

    num_oper_ub
end

function print_final_table(arr_obj, arr_time, arr_active)
    out = "    iter\t    time\t      obj\t   active\n"
    
    for i = 1:length(arr_obj)
        out *= @sprintf("%12d\t%12.1f\t%12.6f\t%12d\n", i, arr_time[i], arr_obj[i], arr_active[i])
    end

    out
end

function solve_Heuristic(obs, operators;
                        time_limit=60, max_iter=100, max_depth=4, stepsize=2)

    ## Initialization
    iter            = 1
    nodes           = Set(1:3)
    max_active      = 3
    remaining_time  = time_limit 
    arr_obj         = []
    arr_time        = [] 
    arr_active      = []
    num_oper_lb     = Dict{Any,Int}()
    num_oper_ub     = Dict{Any,Int}('T' => max_active)

    ## Initial solve
    feasible, obj, time, active, y, y_indexes, ysol = 
        solve_MINLP(nodes, obs, operators, TIME_LIMIT=remaining_time)
    time_limit      -= time
    arr_obj         = [arr_obj; obj]
    arr_time        = [arr_time; time]
    arr_active      = [arr_active; active]

    while feasible && iter < max_iter && obj > 1e-04 && remaining_time > 0
        # Update for the next iteration
        iter += 1
        update_nodes!(nodes, y, y_indexes, max_depth)
        num_oper_lb = update_num_oper_lb(y, y_indexes)
        max_active += stepsize
        num_oper_ub = update_num_oper_ub(max_active, num_oper_lb, stepsize)
        
        # @info "Heuristic Iteration $(iter)" nodes max_active num_oper_lb num_oper_ub
        @info "Heuristic Iteration $(iter)" ysol stepsize

        # Solve        
        feasible, obj, time, active, y, y_indexes, ysol = 
            solve_MINLP(nodes, obs, operators, 
                        # num_oper_lb=num_oper_lb, num_oper_ub=num_oper_ub,
                        ysol=ysol, ysol_dist=stepsize,
                        TIME_LIMIT=remaining_time)
        remaining_time -= time       
        arr_obj         = [arr_obj; obj]
        arr_time        = [arr_time; arr_time[end] + time]
        arr_active      = [arr_active; active]
    end

    @info "Final result\n" * print_final_table(arr_obj, arr_time, arr_active)

    arr_obj, arr_time
end

