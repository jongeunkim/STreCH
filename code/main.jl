using Logging

include("obsgen.jl")
include("opttree.jl")
include("minlp.jl")
include("heur.jl")

## Write a single line of result to csv file
function write_result(FILENAME, obj, active, time, bnbnodes, niters, errmsg)
    (DIR, FILE) = splitdir(FILENAME)
    (DIR, INS) = splitdir(DIR)

    ## Write optimal values
    io = open(DIR * "/result.csv", "a+")
    write(io, @sprintf("%40s,\t%20s,\t%20s,\t%12.6e,\t%12d,\t%12.2f,\t%8d,\t%4d,\t%24s\n", DIR, INS, FILE, obj, active, time, bnbnodes, niters, errmsg))
    close(io)
end


## Main function
function main(args)
    ## Read argument
    i = 0
    FILENAME        = args[i+=1]
    id              = parse(Int, args[i+=1])
    seed            = parse(Int, args[i+=1])
    num_obs         = parse(Int, args[i+=1])
    noise_level     = parse(Float64, args[i+=1])
    time_limit      = parse(Int, args[i+=1])
    method           = args[i+=1]
    method_params    = [args[i+j] for j in 1:(length(args)-i)]

    ## Run a dummy solve to compile
    io = open("dummy.log", "w+")
    logger = SimpleLogger(io, Logging.Info)
    global_logger(logger)
    solve_MINLP(OrderedSet(1:3), rand(5,3), "+-*DC"; scip_time=10, scip_verblevel=1)
    close(io)

    ## Read observations and optimal value
    obs, obs_info = obs_generator(id, seed, num_obs, noise_level)
    optval, temp = obs_optval(id, seed, num_obs, noise_level)

    ## Assign operators and operands
    # Binary:   + - * D
    # Unary:    R (sqrt), E (exp), L (log) 
    # Constant: C (general), P (pi)
    operators = OrderedSet("+-*DRC")

    ## Solve
    io = open(FILENAME * ".log", "w+")
    logger = SimpleLogger(io, Logging.Debug)
    if method == "minlp"
        ## If method is `minlp` then param1 = max_depth and param2 = formulation
        max_depth = parse(Int, method_params[1])
        formulation = method_params[2]

        ## Create nodes by max_depth
        nodes = OrderedSet(1:(2^(max_depth+1)-1))

        ## Solve a minlp
        global_logger(logger)
        feasible, optfeasible, time, obj, ysol, csol, vsol, bnbnodes = solve_MINLP(
            nodes, obs, operators, optimizer_name="SCIP", scip_time=time_limit, formulation=formulation, print_all_solutions=true)
        active = feasible ? length(ysol) : 0
        niters = 0
        errmsg = optfeasible ? "" : "optinfeasible"
        close(io)
    elseif method == "heur"
        ## If method is `heur` then param1 = max_depth and param2 = init_solve
        max_depth = parse(Int, method_params[1])
        init_solve = parse(Int, method_params[2])

        ## Create nodes by max_depth
        nodes = OrderedSet(1:(2^(max_depth+1)-1))
        
        ## Solve a heuristic
        global_logger(logger)
        arr_obj, arr_time, arr_active = 
            solve_Heuristic(nodes, obs, operators;
                            time_limit=time_limit, init_solve=init_solve, obj_termination=optval)
        b = argmin(arr_obj)
        obj = arr_obj[b]
        time = arr_time[end]
        active = arr_active[b]
        bnbnodes = 0
        niters = length(arr_obj)
        errmsg = ""
        close(io)
    elseif method == "optcheck"
        ## Check the optimal value by reading an optimal tree. Currently, only one optimal tree is available (id=37)
        ysol, csol = opttree(id)
        operators = union(Set("+-*D"), intersect(Set("RELCP"), Set(values(ysol))))

        if !isnothing(ysol)
            nodes = Set(keys(ysol))

            global_logger(logger)
            feasible, optfeasible, time, obj, ysol, csol, vsol, bnbnodes = solve_MINLP(nodes, obs, operators, scip_time=time_limit, formulation="Cozad-CR",
                ysol=ysol, ysol_fixlevel=0, csol=csol, print_all_solutions=true)
            active = feasible ? length(ysol) : 0
            bnbnodes = 0
            niters = 0
            errmsg = optfeasible ? "" : "optinfeasible"
            close(io)
        else
            (obj, active, time, niters, errmsg) = (-1, 0, 0, 0, "no optimal solution in opttree.jl")
        end
    end

    write_result(FILENAME, obj, active, time, bnbnodes, niters, errmsg)
end



main(ARGS)