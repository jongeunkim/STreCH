# Copyright (c) <2021>, <Sven Leyffer & Jongeun Kim>
# All rights reserved.

# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree. 

function get_leaves(nodes)
    OrderedSet([n for n in nodes if isempty(intersect(Set([2*n, 2*n+1]), nodes))])
end

function get_nonleaves(nodes)
    leaves = get_leaves(nodes)
    setdiff(nodes, leaves)
end

function get_nodes_complete(nodes)
    """
    Return a set of nodes that its rooted subtree is a complete binary tree.
    Assume that `nodes' is a proper binary tree.
    """
    # nodes_arr = sort(collect(nodes), rev=true)
    height = Dict()
    nodes_complete = Set()
    for n in sort(collect(nodes), rev=true)
        if 2*n in nodes
            if height[2*n] >= 0 && height[2*n] == height[2*n+1]
                height[n] = height[2*n] + 1  
                push!(nodes_complete, n) 
            else
                height[n] = -1  
            end
        else
            height[n] = 0
            push!(nodes_complete, n) 
        end
    end

    OrderedSet(sort(collect(nodes_complete)))
end

function get_nodes_grandparents(nodes)
    """
    Return a set of nodes that has any of 4n, 4n+1, 4n+2, 4n+3
    """
    # nodes_arr = sort(collect(nodes), rev=true)
    nodes_grandparents = Set()
    for n in nodes
        if !isempty(intersect(Set([4*n 4*n+1 4*n+2 4*n+3]), nodes))
            push!(nodes_grandparents, n) 
        end
    end

    OrderedSet(sort(collect(nodes_grandparents)))
end


function fill_single_child(nodes)
    """
    Add a missing child to make a full (proper) binary tree which every node has 0 or 2 children
    """
    
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
    """
    Attach both children at each leaf
    """

    for n in collect(nodes)
        union!(nodes, Set([n, 2*n, 2*n+1]))
    end
    OrderedSet(sort(collect(nodes)))
end

function get_nodes_by_depth(max_depth)
    OrderedSet(1:(2^(max_depth+1)-1))
end

function get_unused_nodes(ysol, nodes_ground, stepsize)
    """
    Return nodes that is large enough to serach k-neighbor of ysol given max_depth and stepsize (distance)
    """

    nodes = Set(keys(ysol))

    for i in 1:Int(floor((stepsize-1)/2))
        nodes = expand_nodes(nodes)
    end

    nodes = fill_single_child(nodes)

    nodes = setdiff(nodes_ground, nodes)
 
    OrderedSet(sort(collect(nodes)))
end







function get_ysol_csol(df_vars)
    df_vars_y = filter(row -> row.vartype=="y" && row.sol > 0.5, df_vars)
    ysol = Dict(row.nodeid => row.operator for row in eachrow(df_vars_y))
    csol = Dict()
    for (k,v) in ysol
        if v == 'C'
            df = filter(row -> row.vartype=="c" && row.nodeid == k, df_vars)
            csol[k] = df.sol[1]
        end
    end
    ysol, csol
end

function get_formula(ysol, csol)
    """
    ysol: Dict(nodeid -> operator)
    csol: Dict(nodeid -> constant value)
    """

    if !(1 in keys(ysol))
        return SortedDict(1=>"0")
    end

    formula = SortedDict()
    for (n,o) in ysol
        if o isa Integer
            formula[n] = "x$o"
        elseif o == 'C'
            formula[n] = "($(csol[n]))"
        elseif o == 'P'
            formula[n] = "pi"
        else
            formula[n] = nothing
        end 
    end

    nodes = reverse(sort(collect(keys(formula))))
    for n in nodes
        if isnothing(formula[n])
            if ysol[n] in Set("+-*")
                formula[n] = "($(formula[2*n])$(ysol[n])$(formula[2*n+1]))"
            elseif ysol[n] == 'D'
                formula[n] = "($(formula[2*n])/$(formula[2*n+1]))"
            elseif ysol[n] == 'R'
                formula[n] = "sqrt($(formula[2*n+1]))"
            elseif ysol[n] == 'S'
                formula[n] = "sqrt(1-$(formula[2*n+1])^2)"
            elseif ysol[n] == 'E'
                formula[n] = "exp($(formula[2*n+1]))"
            elseif ysol[n] == 'L'
                formula[n] = "log($(formula[2*n+1]))"
            else
                println("tree2formula():: somrthing is wrong")
            end
        end
    end

    return formula
end


function compute_mse(y_pred, y_true; model="mse")
    err = abs.(y_pred - y_true)
    
    if model == "mse"
        ##### sqaured
        return Statistics.mean(err.^2)
    elseif model == "rmse"
        return sqrt(Statistics.mean(err.^2))
    elseif model == "alternative"
        ##### alternative of sqaured
        return Statistics.mean(log2.(1 .+ err.*2^30))
    else
        return 1e+06
    end
end

function compute_err_formula(obs, formula, errmodel; err2Inf=false)
    nrows, ncols = size(obs)
    nvars = ncols - 1

    y_true = Float64[]
    y_pred = Float64[]
    for i = 1:nrows
        expr = join(["x$(j) = $(obs[i,j]); " for j = 1:nvars]) * formula

        try
            push!(y_pred, eval(Meta.parse(expr)))
            push!(y_true, obs[i,end])
        catch DomainError
            nothing
            # return err2Inf ? Inf : "DomainError"
        end
    end

    if length(y_pred) <= 0.5 * nrows
        return err2Inf ? Inf : "DomainError"
    else
        return compute_mse(y_pred, y_true, model=errmodel)
    end
end

function read_parameters(file)
    params = Dict{String,Any}()

    for line in eachline(file)
        if length(line) == 0 || line[1] == '#'
            continue
        end

        chunks = split(line, limit=3, keepempty=false)
        if length(chunks) < 2
            continue
        end

        key = filter(x -> !isspace(x), chunks[1]) 
        value = filter(x -> !isspace(x), chunks[2])

        if length(key) >=4
            if "DBL_" == key[1:4]
                value = parse(Float64, value)
            elseif "INT_" == key[1:4]
                value = parse(Int, value)
            else
                nothing
            end
        end

        params[key] = value
    end

    params
end

function save_parameters(file, params)
    open(file, "w+") do io
        for (k,v) in params
            write(io, "$k\t$v\n")
        end
    end
end



function open_logger(logfile, loglevel)
    io = open(logfile, "a+")
    if loglevel == "Debug"
        logger = SimpleLogger(io, Logging.Debug)
    else
        logger = SimpleLogger(io, Logging.Info)
    end
    global_logger(logger)
    io
end

