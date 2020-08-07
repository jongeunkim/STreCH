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
    solve_MINLP(OrderedSet(1:3), rand(5,3), "+-*DC", TIME_LIMIT=10, DISPLAY_VERBLEVEL=0)
    close(io)

    ## Parameters
    rel_gap = 0.00

    ## Observations
    obs, obs_info = obs_generator(id, seed, num_obs, noise_level)
    optval = obs_optval(id, seed, num_obs, noise_level)

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
        solve_MINLP(nodes, obs, operators, TIME_LIMIT=time_limit, print_all_solutions=true)
        close(io)
    elseif model == "heur"
        init_solve  = model_param1
        stepsize    = model_param2
        fixlevel    = model_param3

        global_logger(logger)
        solve_Heuristic(obs, operators, 
                        time_limit=time_limit, init_solve=init_solve, stepsize=stepsize, obj_termination=optval)
        close(io)
    end
end

main(ARGS)