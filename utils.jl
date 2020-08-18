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

# println(get_nonleaves(1:7))
# println(get_nonleaves(1:5))
# println(get_nonleaves(1:3))
# println(get_nonleaves(1:1))
