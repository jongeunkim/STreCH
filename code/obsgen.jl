using CSV
using DataFrames
using DelimitedFiles
using Printf
using Random

include("utils.jl")

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
    obs = nothing
    obs_info = nothing
    seed = seed + id

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

    obs, obs_info = nothing, nothing
    seed = seed + id

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

function obsgen_from_csv(filename, number, num_obs=1, seed=1)
    seed = seed + number


    # Read CSV file
    df = DataFrame(CSV.File(filename; silencewarnings=true))
    
    # Drop missing rows
    DataFrames.dropmissing!(df, :Number)
    
    # Find the row index
    j = findfirst(df[!,:Number] .== number)

    # Generate obs matrix of size (num_obs X num_var) with Uniform(lbs, ubs)
    num_var = df[j, "# variables"]
    obs = initialize_obs(num_obs, num_var, seed)
    lbs = [df[j, "v$(v)_low"] for v in 1:num_var]
    ubs = [df[j, "v$(v)_high"] for v in 1:num_var]
    uniform_lu_array!(obs, lbs, ubs)

    # Read formula
    formula = df[j, :Formula]
    formula = replace(formula, "**" => "^")
    formula = replace(formula, "ln" => "log")

    # Check if there is a non-valid expression
    not_valid_expressions = ["sin", "cos", "tan", "arcsin", "arccos", "arctan"]
    if any([!isnothing(findfirst(e, formula)) for e in not_valid_expressions])
        return false, nothing, nothing, formula, nothing, nothing
    end

    z = zeros(num_obs)
    for i in 1:num_obs
        expr = ""
        for v in 1:num_var
            varname     = df[j, "v$(v)_name"]
            value       = obs[i, v]
            expr *= "$(varname) = $(value); "
        end

        expr *= formula
        
        z[i] = eval(Meta.parse(expr))
    end

    true, obs, z, formula, lbs, ubs
end


function AIF_functions(id, seed, num_obs)
    ## From AI-Feynmann papers

    valid, obs, z, z_info, lbs, ubs = obsgen_from_csv("FeynmanEquations.csv", id-200, num_obs, seed)
    
    if !valid
        return nothing, nothing
    end

    num_var = size(obs, 2)

    # if id == 201
    #     z_info = "AIF-I.9.18"
    #     num_var = 9   # [G m1 m2 x1 x2 y1 y2 z1 z2]
    #     lbs =           [1, 1, 1, 1, 5, 1, 5, 1, 5]
    #     ubs =           [2, 4, 8, 4, 8, 4, 8, 4, 8]

    #     ## Generate independent varaibles
    #     obs = initialize_obs(num_obs, num_var, seed)
    #     uniform_lu_array!(obs, lbs, ubs)

    #     ## Add a dependent varaible
    #     z = ( obs[:,1] .* obs[:,2] .* obs[:,3] ) ./ 
    #         ( (obs[:,4] - obs[:,5]).^2 + (obs[:,6] - obs[:,7]).^2 + (obs[:,8] - obs[:,9]).^2 )
    # elseif id == 202
    #     z_info = "AIF-I.15.3t"
    #     num_var = 4   # [x  c  u  t]
    #     lbs =           [1, 3, 1, 1]
    #     ubs =           [5,10, 2, 5]

    #     ## Generate independent varaibles
    #     obs = initialize_obs(num_obs, num_var, seed)
    #     uniform_lu_array!(obs, lbs, ubs)

    #     ## Add a dependent varaible
    #     x = obs[:,1]
    #     c = obs[:,2]
    #     u = obs[:,3]
    #     t = obs[:,4]
    #     z = ( t - u .* x ./ c.^2 ) ./ sqrt.(1 .- u.^2 ./ c.^2)
    # else
    #     return nothing, nothing
    # end

    obs         = [obs reshape(z, num_obs, 1)]
    obs_info    = join([@sprintf("x%d[%d,%d)", c, lbs[c], ubs[c]) for c in 1:num_var], "\t") * "\t" * z_info 

    obs, obs_info
end


function obs_generator(id, seed, num_obs, noise_level, filename="")
    # Initialization
    obs = nothing
    obs_info = nothing

    if id <= 100
        obs, obs_info = simple_functions(id, seed, num_obs)
    elseif id <= 200
        obs, obs_info = MRA_functions(id, seed, num_obs)
    elseif id <= 300
        obs, obs_info = AIF_functions(id, seed, num_obs)
    end

    if isnothing(obs)
        return nothing, nothing, nothing
    end

    # Add Gaussian noise
    if noise_level > 0
        if id in []
            noise_level = max(1e-01, noise_level)
        elseif id in [115]
            noise_level = max(1e-02, noise_level)
        elseif id in [5]
            noise_level = max(1e-03, noise_level)
        end

        Random.seed!(seed+1000)
        z = obs[:,end] .* (noise_level .* Random.randn((num_obs, 1)) .+ 1)
        optval = compute_mse(z, obs[:,end])
        while optval < 1e-08
            noise_level = 10 * noise_level
            Random.seed!(seed+1000)
            z = obs[:,end] .* (noise_level .* Random.randn((num_obs, 1)) .+ 1)
            optval = compute_mse(z, obs[:,end])
        end

        obs[:,end] = z
    end


    if filename != ""
        open(filename, "w") do io
            DelimitedFiles.write(io, obs_info * "\n")
            DelimitedFiles.writedlm(io, obs)
        end
    end

    obs, obs_info, noise_level
end

function obs_optval(id, seed, num_obs, noise_level)
    obs_noiseless, obs_info, temp = obs_generator(id, seed, num_obs, 0)

    if isnothing(obs_noiseless)
        return nothing, nothing
    end

    obs_noise, obs_info, noise_level = obs_generator(id, seed, num_obs, noise_level)

    mse = compute_mse(obs_noise[:,end], obs_noiseless[:,end])

    mse, noise_level
end

function get_ids(ids, condition_number)
    ids_selected = []
    for id in ids
        obs, obs_info = obs_generator(id, 1, 1, 0, "")

        if !isnothing(obs)
            isexplog    = any([!isnothing(findfirst("exp", obs_info)), !isnothing(findfirst("log", obs_info))])
            num_var     = size(obs,2) - 1

            if condition_number == 0
                push!(valid_ids, id)
            elseif condition_number == 1 && !isexplog && num_var <= 3
                push!(valid_ids, id)
            elseif condition_number == 2 && !isexplog && num_var >= 4
                push!(valid_ids, id)
            elseif condition_number == 3 && isexplog
                push!(valid_ids, id)
            else
                nothing
            end
        end
    end
    return valid_ids
end

# ids = []
# # ids = [ids; 1:6]
# # ids = [ids; [106, 107, 108, 109, 115]]
# ids = [ids; 201:300]
# seeds = [1]
# num_obss = [10]

# for i in 91:100
#     success, obs, z, formula = obsgen_from_csv("FeynmanEquations.csv", i, 1, 1)
#     if !success
#         println(i, ":: Failed, ", formula)
#     end
# end

# Generate obs and write optimal solutions
# DIR = "../obs/"
# if !isdir(DIR)
#     mkdir(DIR)
# end
# io = open(DIR * "optval.log", "w+")
# write(io, @sprintf("%8s\t%8s\t%8s\t%16s\t%16s\n", "id", "seed", "num_obs", "noise_level", "optval"))
# for id in 1:300, seed in [1], num_obs in [10]
#     noise_level = 0
#     obs, obs_info = obs_generator(id, seed, num_obs, noise_level, DIR * "i$(id)_n$(num_obs)_z$(noise_level).obs")
#     for noise_level in [1e-04]
#         optval, noise_level = obs_optval(id, seed, num_obs, noise_level)
#         if !isnothing(optval)
#             write(io, @sprintf("%8d\t%8d\t%8d\t%16.1e\t%16.6e\n", id, seed, num_obs, noise_level, optval))
#         end
#     end
# end
# close(io)

# for condition_number in [0,1,2,3]
#     ids = get_ids(201:300, condition_number)
#     cnt = length(ids)
#     println("condition_number=$(condition_number), cnt=$(cnt), ", join(ids, " "))
# end