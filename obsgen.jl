using DelimitedFiles
using Printf
using Random

function initialize_obs(num_obs, num_var, seed)
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

function uniform_lu_array!(arr, lbs, ubs)
    # Assume that arr is in Uniform[0,1)
    # Transform each column to Uniform[lb,ub)

    for c in 1:length(lbs)
        arr[:,c] .= (arr[:,c] .* (ubs[c] - lbs[c])) .+ lbs[c]
    end
end

function simple_functions(id, seed, num_obs)
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
        obs = initialize_obs(num_obs, num_var, seed)

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
        obs = initialize_obs(num_obs, num_var, seed)

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
        obs = initialize_obs(num_obs, num_var, seed)

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
        obs = initialize_obs(num_obs, num_var, seed)

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
        obs = initialize_obs(num_obs, num_var, seed)

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
        obs = initialize_obs(num_obs, num_var, seed)

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

    obs, obs_info
end

function MRA_functions(id, seed, num_obs)
    ## Table A.5 in Multivariate Rational Approximation paper

    if id == 106
        z_info = "MRA-A.5.6"
        num_var = 2
        lbs = [1, 1]
        ubs = [9, 4]

        ## Generate independent varaibles
        obs = initialize_obs(num_obs, num_var, seed)
        uniform_lu_array!(obs, lbs, ubs)

        ## Add a dependent varaible
        z = ( obs[:,1] + obs[:,2].^3 ) ./ ( obs[:,1] .* obs[:,2].^2 .+ 1 )

    elseif id == 107
        z_info = "MRA-A.5.7"
        num_var = 2
        lbs = [2, 2]
        ubs = [6, 6]

        ## Generate independent varaibles
        obs = initialize_obs(num_obs, num_var, seed)
        uniform_lu_array!(obs, lbs, ubs)

        ## Add a dependent varaible
        z = ( obs[:,1].^2 + obs[:,2].^2 + obs[:,1] - obs[:,2] .- 1 ) ./ ( (obs[:,1] .- 1) .* (obs[:,2] .- 1) ) 
    elseif id == 108
        z_info = "MRA-A.5.8"
        num_var = 2
        lbs = [2, 2]
        ubs = [4, 4]

        ## Generate independent varaibles
        obs = initialize_obs(num_obs, num_var, seed)
        uniform_lu_array!(obs, lbs, ubs)

        ## Add a dependent varaible
        z = ( obs[:,1].^4 + obs[:,2].^4 + obs[:,1].^2 .* obs[:,2].^2 + obs[:,1] .* obs[:,2] ) ./ ( (obs[:,1] .- 1) .* (obs[:,2] .- 1) ) 
    elseif id == 109
        z_info = "MRA-A.5.9"
        num_var = 4
        lbs = [2, 2, 3, 3]
        ubs = [6, 6, 6, 6]

        ## Generate independent varaibles
        obs = initialize_obs(num_obs, num_var, seed)
        uniform_lu_array!(obs, lbs, ubs)

        ## Add a dependent varaible
        z = ( obs[:,1].^2 + obs[:,2].^2 + obs[:,1] - obs[:,2] .- 1 ) ./ ( (obs[:,3] .- 2) .* (obs[:,4] .- 2) )
    elseif id == 115
        z_info = "MRA-A.5.15:Breit-Wigner"
        num_var = 3
        lbs = [80, 5, 90]
        ubs = [100, 10, 93]

        ## Generate independent varaibles
        obs = initialize_obs(num_obs, num_var, seed)
        uniform_lu_array!(obs, lbs, ubs)

        ## Add a dependent varaible
        GAMMA   = sqrt.(obs[:,3].^2 .* (obs[:,3].^2 + obs[:,2].^2))
        PI      = 3
        CONST   = 2 * sqrt(2)

        numer   = CONST * obs[:,3] .* obs[:,2] .* GAMMA
        denom1  = PI * sqrt.(obs[:,1].^2 + GAMMA)
        denom2  = (obs[:,1].^2 - obs[:,3].^2).^2 + obs[:,3].^2 .* obs[:,2].^2
        z       = numer ./ ( denom1 .* denom2 )
    else
        return nothing, nothing
    end

    obs         = [obs reshape(z, num_obs, 1)]
    obs_info    = join([@sprintf("x%d[%d,%d)", c, lbs[c], ubs[c]) for c in 1:num_var], "\t") * "\t" * z_info 

    obs, obs_info
end

function AIF_functions(id, seed, num_obs)
    ## From AI-Feynmann papers

    # Initialization
    num_var = 0
    obs = Float64[]
    obs_info = ""


    
    obs, obs_info
end


function obs_generator(id, seed, num_obs, noise_level, filename="")
    # Initialization
    obs = Float64[]
    obs_info = ""

    if id < 100
        obs, obs_info = simple_functions(id, seed, num_obs)
    elseif id < 200
        obs, obs_info = MRA_functions(id, seed, num_obs)
    elseif id < 300
        obs, obs_info = AIF_functions(id, seed, num_obs)
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
for id in 1:6, seed in [1], num_obs in [10], noise_level in [0]
    obs_generator(id, seed, num_obs, noise_level, "obs/i$id.obs")
end

for id in [106, 107, 108, 109, 115], seed in [1], num_obs in [10], noise_level in [0]
    obs_generator(id, seed, num_obs, noise_level, "obs/i$id.obs")
end


