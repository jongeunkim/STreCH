using DataStructures

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
