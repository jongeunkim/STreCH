include("symbolic.jl")
include("utils.jl")

# const PARAMETERSET = [(fixlevel,smin,smax,time) for time in [300, 600, 1200, 2400, 4800, 9600] for (smin, smax) in [(0,3), (4,5), (6,7)] for fixlevel in [1,2]] 
const PARAMETERSET = [(fixlevel,smin,smax,nodelimit) for nodelimit in [1000, 3000, 1e+04, 3e+04, 1e+05, 3e+05, 1e+06, 3e+06] for (smin, smax) in [(0,3), (4,5), (6,7)] for fixlevel in [1,2]] 

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
                        time_limit=60, obj_termination=1e-06,
                        init_solve=1, subsampling=1,
                        max_depth=4, min_improvement=0.01)

    ## Initialization
    iter            = 1
    nodes           = Set([])
    remaining_time  = time_limit
    b               = 1     # the index of the best solution 
    p               = 1     # the index of the current parameter
    arr_obj         = []
    arr_time        = [] 
    arr_active      = []
    arr_ysol        = []
    arr_csol        = []
    arr_p           = []

    ## Initial solve
    nodes = Set(1:(2^(init_solve+1)-1))
    # TIME_LIMIT = remaining_time >= 301 ? 300 : 60
    TIME_LIMIT = remaining_time
    NODE_LIMIT = 1000
    feasible, optfeasible, time, obj, ysol, csol, vsol = solve_MINLP(nodes, obs, operators, TIME_LIMIT=TIME_LIMIT, NODE_LIMIT=NODE_LIMIT)

    
    remaining_time -= time
    arr_obj         = [arr_obj; obj]
    arr_time        = [arr_time; time]
    arr_active      = [arr_active; length(ysol)]
    arr_ysol        = [arr_ysol; deepcopy(ysol)]
    arr_csol        = [arr_csol; deepcopy(csol)]
    arr_p           = [arr_p; p]

    while obj > obj_termination + eps(Float64) && remaining_time > 0 
        # Update for the next iteration
        iter += 1

        # Set parameters
        # (ysol_fix_level, ysol_dist_min, ysol_dist, TIME_LIMIT) = PARAMETERSET[p] 
        # TIME_LIMIT = max(10, min(TIME_LIMIT, remaining_time))
        # NODE_LIMIT = -1
        # @info "Heuristic Iteration $(iter)" ysol ysol_fix_level ysol_dist_min ysol_dist TIME_LIMIT

        (ysol_fix_level, ysol_dist_min, ysol_dist, NODE_LIMIT) = PARAMETERSET[p] 
        TIME_LIMIT = max(10, remaining_time)
        @info "Heuristic Iteration $(iter)" ysol ysol_fix_level ysol_dist_min ysol_dist NODE_LIMIT

        nodes       = update_nodes(ysol, max_depth, ysol_dist)

        # Subsampling
        if length(nodes) > 7
            vsol, err       = get_vsol(arr_ysol[b], arr_csol[b], obs)
            abserr          = abs.(vsol[:,1] - obs[:,end])
            subobs_size     = Int64(round(size(obs,1) * subsampling))
            subobs_indexes  = sortperm(abserr, rev=true)[1:subobs_size]
            subobs          = obs[subobs_indexes,:]
        else
            subobs          = obs
        end


        # Solve
        feasible, optfeasible, time, obj, ysol, csol, vsol  = 
            solve_MINLP(nodes, subobs, operators, 
                        ysol=ysol, ysol_dist=ysol_dist, ysol_dist_min=ysol_dist_min, ysol_fix_level=ysol_fix_level,
                        ABS_GAP=(1-min_improvement)*minimum(arr_obj),
                        TIME_LIMIT=TIME_LIMIT, NODE_LIMIT=NODE_LIMIT,
                        PRESOLVE=(TIME_LIMIT < 300 ? 0 : 1))
        remaining_time -= time

        # Recompute the obj for the whole observations
        if feasible
            vsol, err = get_vsol(ysol, csol, obs)
            obj = compute_mse(vsol[:,1], obs[:,end])
        else
            obj = 1e+09
        end

        arr_obj         = [arr_obj; obj]
        arr_time        = [arr_time; arr_time[end] + time]
        arr_active      = [arr_active; feasible ? length(ysol) : 0]
        arr_p           = [arr_p; p]
        arr_ysol        = [arr_ysol; feasible ? deepcopy(ysol) : nothing]
        arr_csol        = [arr_csol; feasible ? deepcopy(csol) : nothing]

        epsilon = 1e-12
        if !feasible || (arr_obj[end] >=  arr_obj[b] + epsilon)
            ysol            = arr_ysol[b]
            p              += 1
        elseif abs(arr_obj[end] - arr_obj[b]) < epsilon            
            p              += 1
        else
            b               = length(arr_obj)
            p               = 1
        end

        @info "Final result\n" * print_final_table(arr_obj, arr_time, arr_active)
        @info "p" arr_p
    end


    arr_obj, arr_time, arr_active
end
