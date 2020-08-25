using Logging

include("obsgen.jl")
include("opttree.jl")
include("minlp.jl")
include("heur.jl")

## Functions

function write_result(FILENAME, obj, active, time, niters, errmsg)
    (DIR, FILE) = splitdir(FILENAME)
    (DIR, INS) = splitdir(DIR)

    ## Write optimal values
    io = open(DIR * "/result.csv", "a+")
    write(io, @sprintf("%40s,\t%20s,\t%20s,\t%12.6e,\t%12d,\t%12.2f,\t%4d,\t%24s\n", DIR, INS, FILE, obj, active, time, niters, errmsg))
    close(io)
end


function main(args)
    i = 0
    FILENAME        = args[i+=1]
    id              = parse(Int, args[i+=1])
    seed            = parse(Int, args[i+=1])
    num_obs         = parse(Int, args[i+=1])
    noise_level     = parse(Float64, args[i+=1])
    time_limit      = parse(Int, args[i+=1])
    model           = args[i+=1]
    model_param1    = length(args) >= i+1 ? parse(Int, args[i+=1]) : 0
    model_param2    = length(args) >= i+1 ? parse(Int, args[i+=1]) : 0
    model_param3    = length(args) >= i+1 ? parse(Int, args[i+=1]) : 0

    ## Dummy solve
    io = open("dummy.log", "w+")
    logger = SimpleLogger(io, Logging.Info)
    global_logger(logger)
    solve_MINLP(OrderedSet(1:3), rand(5,3), "+-*DC"; scip_time=10, scip_verblevel=1)
    close(io)

    ## Parameters
    rel_gap = 0.00

    ## Observations
    obs, obs_info = obs_generator(id, seed, num_obs, noise_level)
    optval, temp = obs_optval(id, seed, num_obs, noise_level)

    ## Operators and operands
    # Binary:   + - * D
    # Unary:    R (sqrt), E (exp), L (log) 
    # Constant: C (integer), P (pi)
    operators = OrderedSet("+-*DRCP")

    ## Solve
    io = open(FILENAME * ".log", "w+")
    logger = SimpleLogger(io, Logging.Info)
    if model == "minlp"
        max_depth = model_param1
        if model_param2 == 1
            formulation = "Cozad"
        elseif model_param2 == 2
            formulation = "Cozad-CR"
        elseif model_param2 == 3
            formulation = "New"
        elseif model_param2 == 4
            formulation = "New-NR"
        end
        nodes = OrderedSet(1:(2^(max_depth+1)-1))

        global_logger(logger)
        feasible, optfeasible, time, obj, ysol, csol, vsol = solve_MINLP(nodes, obs, operators, scip_time=time_limit, formulation=formulation, print_all_solutions=true)
        active = feasible ? length(ysol) : 0
        niters = 0
        errmsg = optfeasible ? "" : "optinfeasible"
        close(io)
    elseif model == "optcheck"
        ysol, csol = opttree(id)
        operators = union(Set("+-*D"), intersect(Set("RELCP"), Set(values(ysol))))

        if !isnothing(ysol)
            nodes = Set(keys(ysol))

            global_logger(logger)
            feasible, optfeasible, time, obj, ysol, csol, vsol = solve_MINLP(nodes, obs, operators, scip_time=time_limit, 
                ysol=ysol, ysol_dist=0, ysol_dist_min=0, ysol_fix_level=0, csol=csol, print_all_solutions=true)
            active = feasible ? length(ysol) : 0
            niters = 0
            errmsg = optfeasible ? "" : "optinfeasible"
            close(io)
        else
            (obj, active, time, niters, errmsg) = (-1, 0, 0, 0, "no optimal solution in opttree.jl")
        end
    elseif model == "heur"
        max_depth = model_param1
        init_solve  = model_param2
        # subsampling = 0.01 * model_param3
        
        global_logger(logger)
        arr_obj, arr_time, arr_active = 
            solve_Heuristic(obs, operators, time_limit=time_limit, max_depth=max_depth,
                            init_solve=init_solve, obj_termination=optval)
        b = argmin(arr_obj)
        obj = arr_obj[b]
        time = arr_time[end]
        active = arr_active[b]
        niters = length(arr_obj)
        errmsg = ""
        close(io)
    end

    write_result(FILENAME, obj, active, time, niters, errmsg)
end

main(ARGS)