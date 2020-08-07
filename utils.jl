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
