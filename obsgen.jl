using DelimitedFiles
using Printf
using Random


function obs_generator(id, seed, num_obs, noise_level, filename="")
    function initialize_obs(num_var, seed)
        obs = zeros(num_obs, num_var)
        Random.seed!(seed)
        Random.rand!(obs)
        obs
    end
    
    function uniform_lu!(arr, cols, lb, ub)
        # Assume that arr is in Uniform[0,1)
        # Transform each column to Uniform[lb,ub)

        for c in cols
            arr[:,c] .= (arr[:,c] .* (ub - lb)) .+ lb
        end
    end

    # Initialization
    num_var = 0
    obs = Float64[]
    obs_info = ""

    if id == 1
        # depth = 2
        # +-*/C

        obs_info = "x1_[1,2)\tx2_[3,4)\tf=x1*x2+x2+1"
        
        ## Generate independent varaibles in Uniform[0,1)
        num_var = 2
        initialize_obs(num_var, seed)

        ## Affine-transform indep variables
        uniform_lu!(obs, [1], 1, 2)
        uniform_lu!(obs, [2], 3, 4)

        ## Add a dependent varaible
        z = obs[:,1] .* obs[:,2] + obs[:,2] + ones(num_obs)
        obs = [obs reshape(z, num_obs, 1)]
    elseif id == 2
        # depth = 3
        # +-*/

        obs_info = "x1_[1,2)\tx2_[3,4)\tx3_[5,6)\tf=x1*(x2+x3)+2*x2*x3"
        
        ## Generate independent varaibles in Uniform[0,1)
        num_var = 3
        initialize_obs(num_var, seed)

        ## Affine-transform indep variables
        uniform_lu!(obs, [1], 1, 2)
        uniform_lu!(obs, [2], 3, 4)
        uniform_lu!(obs, [3], 5, 6)

        ## Add a dependent varaible
        z = obs[:,1] .* (obs[:,2] + obs[:,3]) +
            2 .* obs[:,2] .* obs[:,3]
        obs = [obs reshape(z, num_obs, 1)]
    elseif id == 3
        # depth = 3
        # +-*/

        obs_info = "x1_[1,3)\tx2_[2,6)\tx3_[6,9)\tx4_[6,8)\tf=(x3+x2-x1)/(x4-x2+x1)"
        
        ## Generate independent varaibles in Uniform[0,1)
        num_var = 4
        initialize_obs(num_var, seed)

        ## Affine-transform indep variables
        uniform_lu!(obs, [1], 1, 3)
        uniform_lu!(obs, [2], 2, 6)
        uniform_lu!(obs, [3], 6, 9)
        uniform_lu!(obs, [4], 6, 8)

        ## Add a dependent varaible
        z = (obs[:,3] .+ obs[:,2] .- obs[:,1]) ./ (obs[:,4] .- obs[:,2] .+ obs[:,1])
        obs = [obs reshape(z, num_obs, 1)]
    elseif id == 4
        # depth = 4
        # +-*/

        obs_info = "x1_[4,6)\tx2_[1,3)\tf=(x1^2-x2^2+x1x2)/(x1x2-x1+x2)"
        
        ## Generate independent varaibles in Uniform[0,1)
        num_var = 2
        initialize_obs(num_var, seed)

        ## Affine-transform indep variables
        uniform_lu!(obs, [1], 4, 6)
        uniform_lu!(obs, [2], 1, 3)

        ## Add a dependent varaible
        z = ( obs[:,1].^2 - obs[:,2].^2 + (obs[:,1] .* obs[:,2]) ) ./ 
            ( (obs[:,1] .* obs[:,2]) - obs[:,1] + obs[:,2])
        obs = [obs reshape(z, num_obs, 1)]
    elseif id == 5
        # depth = 4
        # +-*/

        obs_info = "x1_[1,6)\tx2_[1,3)\tx3_[1,2)\tx4_[3,4)\tx5_[5,6)\tf=(x1^2*x2)/(x3^2+x4^2+x5^2)"
        
        ## Generate independent varaibles in Uniform[0,1)
        num_var = 5
        initialize_obs(num_var, seed)

        ## Affine-transform indep variables
        uniform_lu!(obs, [1], 1, 6)
        uniform_lu!(obs, [2], 1, 3)
        uniform_lu!(obs, [3], 1, 2)
        uniform_lu!(obs, [4], 1, 4)
        uniform_lu!(obs, [5], 1, 6)

        ## Add a dependent varaible
        z = ( obs[:,1].^2 .* obs[:,2] ) ./ 
            ( obs[:,3].^2 + obs[:,4].^2 + obs[:,5].^2 )
        obs = [obs reshape(z, num_obs, 1)]
    elseif id == 6
        # depth = 4
        # +-*/

        obs_info = "x1_[1,6)\tx2_[1,3)\tx3_[1,2)\tx4_[1,4)\tx5_[1,6)\tf=(x3^2+x4^2+x5^2)/(x1*x2)"
        
        ## Generate independent varaibles in Uniform[0,1)
        num_var = 5
        initialize_obs(num_var, seed)

        ## Affine-transform indep variables
        uniform_lu!(obs, [1], 1, 6)
        uniform_lu!(obs, [2], 1, 3)
        uniform_lu!(obs, [3], 1, 2)
        uniform_lu!(obs, [4], 1, 4)
        uniform_lu!(obs, [5], 1, 6)

        ## Add a dependent varaible
        z = ( obs[:,3].^2 + obs[:,4].^2 + obs[:,5].^2 ) ./ 
            ( obs[:,1].^2 .* obs[:,2] )
        obs = [obs reshape(z, num_obs, 1)]
    else
        nothing
    end

    # Add Gaussian noise
    Random.seed!(seed+1000)
    obs[:,end] = obs[:,end] .* (noise_level .* Random.randn((num_obs, 1)) .+ 1)

    if filename != ""
        open(filename, "w") do io
            DelimitedFiles.write(io, obs_info * "\n")
            DelimitedFiles.writedlm(io, obs)
        end
    end

    obs, obs_info
end


## Generate obs
# for id in [1], seed in [1], num_obs in [5], noise_level in [0, 1e-02, 1e-04, 1e-06]
#     obs_generator(id, seed, num_obs, noise_level)
# end


