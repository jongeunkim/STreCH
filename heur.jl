include("symbolic.jl")

function remove_single_child(nodes)
    nodes = Set(nodes)
    updated = true
    
    while updated
        updated = false
        for n in collect(nodes)
            if !isempty(intersect(nodes, Set([2*n, 2*n+1]))) && !issubset(Set([2*n, 2*n+1]), nodes)
                setdiff!(nodes, Set([2*n, 2*n+1]))
                updated = true
                break
            end
        end
    end

    OrderedSet(sort(collect(nodes)))
end

function fill_single_child(nodes)
    nodes = Set(nodes)
    updated = true
    
    while updated
        updated = false
        for n in collect(nodes)
            if !isempty(intersect(nodes, Set([2*n, 2*n+1]))) && !issubset(Set([2*n, 2*n+1]), nodes)
                union!(nodes, Set([2*n, 2*n+1]))
                updated = true
                break
            end
        end
    end

    OrderedSet(sort(collect(nodes)))
end

function expand_nodes(nodes)
    for n in collect(nodes)
        union!(nodes, Set([n, 2*n, 2*n+1]))
    end
    OrderedSet(sort(collect(nodes)))
end

function update_nodes(ysol, max_depth, stepsize)
    nodes = Set(keys(ysol))

    for i in 1:Int(floor((stepsize-1)/2))
        nodes = expand_nodes(nodes)
    end

    nodes = fill_single_child(nodes)

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



function solve_Heuristic(obs, operators; 
                        time_limit=60, stepsize=2, fixlevel=-1, init_solve=1, obj_termination=1e-06,
                        max_iter=100, max_depth=10, min_improvement=0.1)

    ## Initialization
    iter            = 1
    nodes           = remove_single_child(Set(1:max(7, stepsize)))
    # max_active      = 3
    remaining_time  = time_limit 
    arr_obj         = []
    arr_time        = [] 
    arr_active      = []
    arr_ysol        = []
    arr_stepsize    = []
    stepsize_extra  = 0
    # num_oper_lb     = Dict{Any,Int}()
    # num_oper_ub     = Dict{Any,Int}('T' => max_active)

    ## Initial solve
    if init_solve == 1
        # ysol = 0 and ysol_dist=stepsize

        feasible, time, obj, ysol, csol, vsol = 
            solve_MINLP(nodes, obs, operators, ysol=Dict(), ysol_dist=stepsize, TIME_LIMIT=remaining_time)
    elseif init_solve == 2
        # Depth two tree problem
        nodes = Set(1:7)
        feasible, time, obj, ysol, csol, vsol = 
            solve_MINLP(nodes, obs, operators, TIME_LIMIT=remaining_time)
    end
    
    time_limit      -= time
    arr_obj         = [arr_obj; obj]
    arr_time        = [arr_time; time]
    arr_active      = [arr_active; length(ysol)]
    arr_ysol        = [arr_ysol; deepcopy(ysol)]
    arr_stepsize    = [arr_stepsize; stepsize]

    while obj > obj_termination + eps(Float64) && remaining_time > 0 && iter < max_iter  
        # Update for the next iteration
        iter += 1
        nodes       = update_nodes(ysol, max_depth, stepsize+stepsize_extra)
        # num_oper_lb = update_num_oper_lb(y, y_indexes)
        # max_active += stepsize
        # num_oper_ub = update_num_oper_ub(max_active, num_oper_lb, stepsize)

        # @info "Heuristic Iteration $(iter)" nodes max_active num_oper_lb num_oper_ub
        @info "Heuristic Iteration $(iter)" ysol stepsize+stepsize_extra

        # Solve
        ysol_dist       = stepsize + stepsize_extra
        ysol_dist_min   = stepsize_extra > 0 ? stepsize+stepsize_extra : 0
        feasible, time, obj, ysol, csol, vsol  = 
            solve_MINLP(nodes, obs, operators, 
                        # num_oper_lb=num_oper_lb, num_oper_ub=num_oper_ub,
                        ysol=ysol, ysol_dist=ysol_dist, ysol_dist_min=ysol_dist_min, ysol_fix_level=fixlevel,
                        ABS_GAP=(1-min_improvement)*minimum(arr_obj),
                        TIME_LIMIT=remaining_time)
        remaining_time -= time       
        arr_obj         = [arr_obj; obj]
        arr_time        = [arr_time; arr_time[end] + time]
        arr_active      = [arr_active; length(ysol)]
        arr_stepsize    = [arr_stepsize; stepsize + stepsize_extra]
        arr_ysol        = [arr_ysol; deepcopy(ysol)]

        epsilon = 1e-12
        if arr_obj[end] >=  minimum(arr_obj[1:end-1]) + epsilon
            ysol = arr_ysol[argmin(arr_obj[1:end-1])]
            stepsize_extra += 1
        elseif abs(arr_obj[end] - minimum(arr_obj[1:end-1])) < epsilon
            stepsize_extra += 1
        else
            stepsize_extra  = 0
        end
    end

    @info "Final result\n" * print_final_table(arr_obj, arr_time, arr_active)

    @info "Stepsize" arr_stepsize

    arr_obj, arr_time
end
