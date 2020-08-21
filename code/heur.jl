include("symbolic.jl")
include("utils.jl")

function solve_Heuristic(obs, 
                        operators; 
                        time_limit=60, 
                        obj_termination=1e-06,
                        init_solve=1, 
                        subsampling=1,
                        max_depth=4, 
                        min_improvement=0.01)

                        # PARAMETERSET = [(fixlevel,smin,smax,time) for time in [300, 600, 1200, 2400, 4800, 9600] for (smin, smax) in [(0,3), (4,5), (6,7)] for fixlevel in [1,2]]
    PARAMETERSET = [(fixlevel,smin,smax,nodelimit) for nodelimit in [10^i for i in 4:10] for (smin, smax) in [(0,3), (4,5), (6,7)] for fixlevel in [1,2]] 

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
    # (TIME_LIMIT, NODE_LIMIT) = (min(0.1 * remaining_time, 300), -1)
    (TIME_LIMIT, NODE_LIMIT) = (remaining_time, 10 * PARAMETERSET[1][4])
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
