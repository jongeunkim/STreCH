"""
utils.jl
Collection of utility functions
"""

using DataStructures
using Printf

function compute_mse(y_pred, y_true)
    err = y_pred - y_true
    mse = sum(err.^2) / length(y_true)
    mse
end

function get_leaves(nodes)
    OrderedSet([n for n in nodes if isempty(intersect(Set([2*n, 2*n+1]), nodes))])
end

function get_nonleaves(nodes)
    leaves = get_leaves(nodes)
    setdiff(nodes, leaves)
end

function get_treesol(ysol, csol)
    treesol = deepcopy(ysol)
    
    for (n,o) in treesol
        if o == 'C'
            treesol[n] = @sprintf("%.0f", csol[n])
        elseif typeof(o) == Int
            treesol[n] = @sprintf("x%d", o)
        end
    end

    treesol
end

function get_vsol(ysol, csol, obs)
    domain_error = false
    nodes_rev = sort(collect(keys(ysol)), rev=true)
    num_obs = size(obs)[1]
    vsol = zeros(num_obs, maximum(nodes_rev))

    for n in nodes_rev
        o = ysol[n]

        if typeof(o) == Int
            vsol[:,n] = obs[:,o]
        elseif o == 'C'
            vsol[:,n] .= csol[n]
        elseif o == 'P'
            vsol[:,n] .= pi
        elseif o in '+'
            for i in 1:num_obs
                vsol[i,n] = vsol[i,2*n] + vsol[i,2*n+1]
            end
        elseif o == '-'
            for i in 1:num_obs
                vsol[i,n] = vsol[i,2*n] - vsol[i,2*n+1]
            end
        elseif o == '*'
            for i in 1:num_obs
                vsol[i,n] = vsol[i,2*n] * vsol[i,2*n+1]
            end
        elseif o == 'D'
            for i in 1:num_obs
                if abs(vsol[i,2*n+1]) < eps(Float64)
                    @info "domain error in compute_vsol!()" n o
                    domain_error = true
                    break
                end
                vsol[i,n] = vsol[i,2*n] / vsol[i,2*n+1]
            end
        elseif o == 'R'
            for i in 1:num_obs
                if vsol[i,2*n+1] < 0
                    @info "domain error in compute_vsol!()" n o
                    domain_error = true
                    break
                end
                vsol[i,n] = sqrt(vsol[i,2*n+1])
            end
        elseif o == 'E'
            for i in 1:num_obs
                vsol[i,n] = exp(vsol[i,2*n+1])
            end
        elseif o == 'L'
            for i in 1:num_obs
                if vsol[i,2*n+1] <= 0
                    @info "domain error in compute_vsol!()" n o
                    domain_error = true
                    break
                end
                vsol[i,n] = log(vsol[i,2*n+1])
            end
        else
            println("Not a valid operator at node $n")
        end
    end

    if domain_error
        vsol .= 0
    end

    vsol, domain_error
end

function print_tree(treesol)
    out = "Print tree\n"

    depth = Int(floor(log2(maximum(keys(treesol)))))
    for d in 0:depth
        out = out * @sprintf("d=%2d\t\t", d)

        for n in 2^d:(2^(d+1)-1)
            if n in keys(treesol)
                out = out * @sprintf("%d:%s  ", n, treesol[n])
            end
        end

        out = out * "\n"
    end

    out
end


function print_final_table(arr_obj, arr_time, arr_active)
    out = "  i/iter\t    time\t      obj\t   obj(e)\t   active\n"
    
    for i = 1:length(arr_obj)
        out *= @sprintf("%12d\t%12.1f\t%12.6f\t%12.6e\t%12d\n", i, arr_time[i], arr_obj[i], arr_obj[i], arr_active[i])
    end

    b = argmin(arr_obj)
    out *= repeat("-", 80) * "\n"
    out *= @sprintf("%12s\t%12.1f\t%12.6f\t%12.6e\t%12d\n", "Best", arr_time[b], arr_obj[b], arr_obj[b], arr_active[b])

    out
end

function print_obs(obs)
    (m,n) = size(obs)

    out = "Print obs\n"
    out *= "i\\var\t"
    for j in 1:n-1
        out *= @sprintf("x%d\t", j)
    end
    out *= @sprintf("z\n")
    
    for i in 1:m
        out *= @sprintf("i=%2d\t", i)
        for j in 1:n
            out *= @sprintf("%.3f\t", obs[i,j])
        end
        out *= @sprintf("\n")
    end   
    
    out
end

function remove_single_child(nodes)
    """
    Remove a sole child to make a full (proper) binary tree which every node has 0 or 2 children
    """

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

    # println("height")
    # for n in sort(collect(nodes))
    #     println("height[$(n)] = $(height[n])")
    # end

    # println("nodes_complete")
    # println(nodes_complete)

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

function update_nodes(ysol, nodes_ground, stepsize)
    """
    Return nodes that is large enough to serach k-neighbor of ysol given max_depth and stepsize (distance)
    """

    nodes = Set(keys(ysol))

    for i in 1:Int(floor((stepsize-1)/2))
        nodes = expand_nodes(nodes)
    end

    nodes = fill_single_child(nodes)

    intersect!(nodes, nodes_ground)
 
    OrderedSet(sort(collect(nodes)))
end

function print_error(vsol, obs)
    """
    Compare predicted values and true values and return string listing all errors
    """

    err     = vsol[:,1] - obs[:,end]
    abserr  = abs.(err)

    out = "i\terr\n"
    for i in sortperm(abserr)
        out *= @sprintf("%2d\t%12.6f\n", i, err[i])
    end

    out
end

