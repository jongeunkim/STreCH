using Logging

include("obsgen.jl")
include("symbolic.jl")
include("heur.jl")

## Functions
read_obs = (name) -> DelimitedFiles.readdlm(name * ".obs", header=true)

function main(args)
    FILENAME        = args[1]
    id              = parse(Int, args[2])
    seed            = parse(Int, args[3])
    num_obs         = parse(Int, args[4])
    noise_level     = parse(Float64, args[5])
    max_depth       = parse(Int, args[6])
    time_limit      = parse(Int, args[7])
    model           = args[8]
    stepsize        = length(args) >= 9 ? parse(Int, args[9]) : 0

    ## Dummy solve
    io = open("dummy.log", "w+")
    logger = SimpleLogger(io, Logging.Info)
    global_logger(logger)
    solve_MINLP(OrderedSet(1:3), rand(5,3), "+-*DC", TIME_LIMIT=10)
    close(io)

    ## Parameters
    rel_gap = 0.00

    ## Observations
    obs_generator(id, seed, num_obs, noise_level, FILENAME * ".obs")
    obs, obs_info = read_obs(FILENAME)

    ## Operators and operands
    # Binary:   + - * D
    # Unary:    R (sqrt), E (exp), L (log) 
    # Constant: C
    operators = OrderedSet("+-*D")

    ## Tree
    nodes = OrderedSet(1:(2^(max_depth+1)-1))

    ## Solve
    io = open(FILENAME * ".log", "w+")
    logger = SimpleLogger(io, Logging.Info)
    global_logger(logger)
    if model == "minlp"
        solve_MINLP(nodes, obs, operators, TIME_LIMIT=time_limit, REL_GAP=rel_gap)
    elseif model == "heur"
        solve_Heuristic(obs, operators, time_limit=time_limit, max_depth=max_depth, stepsize=stepsize)
    end
    close(io)
end

main(ARGS)