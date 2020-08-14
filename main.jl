using Logging

include("obsgen.jl")
include("symbolic.jl")
include("heur.jl")

## Functions
read_obs = (name) -> DelimitedFiles.readdlm(name * ".obs", header=true)

# function solve_MINLP_with_logging(nodes, obs, operators, time_limit,
#                                     logger1, logger2)
#     solve_MINLP(nodes, obs, operators, TIME_LIMIT=time_limit)
# end

function write_result(FILENAME, obj, active, time, niters, errmsg)
    (DIR, FILE) = splitdir(FILENAME)
    (DIR, dummy) = splitdir(DIR)

    ## Write optimal values
    io = open(DIR * "/result.txt", "a+")
    write(io, @sprintf("%80s\t%12.6f\t%12.6e\t%12d\t%12.2f\t%4d\t%24s\n", FILENAME, obj, obj, active, time, niters, errmsg))
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
    solve_MINLP(OrderedSet(1:3), rand(5,3), "+-*DC", TIME_LIMIT=10, DISPLAY_VERBLEVEL=1)
    close(io)

    ## Parameters
    rel_gap = 0.00

    ## Observations
    obs, obs_info = obs_generator(id, seed, num_obs, noise_level)
    optval, temp = obs_optval(id, seed, num_obs, noise_level)

    ## Operators and operands
    # Binary:   + - * D
    # Unary:    R (sqrt), E (exp), L (log) 
    # Constant: C
    operators = OrderedSet("+-*DRC")

    ## Solve
    io = open(FILENAME * ".log", "w+")
    logger = SimpleLogger(io, Logging.Info)
    if model == "minlp"
        max_depth = model_param1
        nodes = OrderedSet(1:(2^(max_depth+1)-1))

        global_logger(logger)
        feasible, optfeasible, time, obj, ysol, csol, vsol = solve_MINLP(nodes, obs, operators, TIME_LIMIT=time_limit, print_all_solutions=true)
        active = feasible ? length(ysol) : 0
        niters = 0
        errmsg = optfeasible ? "" : "optinfeasible"
        close(io)
    elseif model == "heur"
        init_solve  = model_param1
        subsampling = 0.01 * model_param2

        global_logger(logger)
        arr_obj, arr_time, arr_active = 
            solve_Heuristic(obs, operators, time_limit=time_limit, init_solve=init_solve, subsampling=subsampling, obj_termination=optval)
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