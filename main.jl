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
    max_depth       = parse(Int, args[i+=1])
    time_limit      = parse(Int, args[i+=1])
    model           = args[i+=1]
    stepsize        = length(args) >= i+1 ? parse(Int, args[i+=1]) : 0

    ## Dummy solve
    io = open("dummy.log", "w+")
    logger = SimpleLogger(io, Logging.Info)
    global_logger(logger)
    solve_MINLP(OrderedSet(1:3), rand(5,3), "+-*DC", TIME_LIMIT=10)
    close(io)

    ## Parameters
    rel_gap = 0.00

    ## Observations
    obs, obs_info = obs_generator(id, seed, num_obs, noise_level) #, FILENAME * ".obs")
    # obs, obs_info = read_obs(FILENAME)

    ## Operators and operands
    # Binary:   + - * D
    # Unary:    R (sqrt), E (exp), L (log) 
    # Constant: C
    operators = OrderedSet("+-*DRC")

    ## Tree
    nodes = OrderedSet(1:(2^(max_depth+1)-1))

    ## Solve
    io = open(FILENAME * ".log", "w+")
    logger = SimpleLogger(io, Logging.Info)
    if model == "minlp"
        global_logger(logger)
        solve_MINLP(nodes, obs, operators, TIME_LIMIT=time_limit)
        close(io)
    elseif model == "heur"
        global_logger(logger)
        solve_Heuristic(obs, operators, time_limit=time_limit, max_depth=max_depth, stepsize=stepsize)
        close(io)
    end
end

main(ARGS)