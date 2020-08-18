include("obsgen.jl")

function get_vbounds(obs, operators=Set())
    if isnothing(obs)
        return nothing, nothing
    end        

    # Initialization
    operators = Set(operators)
    (v_lb, v_ub) = (-1e+09, 1e+09)
    factor1 = 10
    factor2 = 2
    eps = 1e-09

    # if !isempty(intersect(operators, Set("EL")))
    #     factor1 = 5
    #     factor2 = 2
    # end

    ## Compute min/max of indep/dep variables
    # (z_lb, z_ub) = (minimum(obs[:,end]), maximum(obs[:,end]))
    # (x_lb, x_ub) = (minimum(obs[:,1:end-1]), maximum(obs[:,1:end-1]))
    (obs_lb, obs_ub) = (minimum(obs), maximum(obs))

    ## (1) Shrink bounds
    (v_lb, v_ub) = (maximum([v_lb, factor1^(-sign(obs_lb)) * obs_lb]), minimum([v_ub, factor1 * obs_ub]))

    ## (2) Expand bounds if it is too tight
    (v_lb, v_ub) = (minimum([v_lb, factor2^(-sign(obs_lb)) * obs_lb]), maximum([v_ub, factor2 * obs_ub]))

    ## (3) make bounds symmetric and +- epsilon
    magnitude = maximum([abs(v_lb), abs(v_ub)]) + eps
    (v_lb, v_ub) = (-magnitude, magnitude)

    v_lb, v_ub
end

# seed=1
# num_obs=10
# noise_level=1e-04

# ids = [201 202 203 243 244 248 280 286 287 294]

# for id in ids
#     obs, obs_info = obs_generator(id, seed, num_obs, noise_level)
#     v_lb, v_ub = get_vbounds(obs ,"EL")
#     println("$id\t $v_lb\t $v_ub")
# end