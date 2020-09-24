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
    # Set the seed
    seed = seed + number

    # Read CSV file
    df = DataFrame(CSV.File(filename; silencewarnings=true))
    
    # Drop missing rows
    DataFrames.dropmissing!(df, :Number)
    
    # Find the row index
    j = findfirst(df[!,:Number] .== number)
    if isnothing(j)
        println("there is no number == $(number)")
        return nothing, nothing, formula
    end

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
        println("there is a invalid operation in the formula, $(formula)")
        return nothing, nothing, formula
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

    obs      = [obs reshape(z, num_obs, 1)]
    obs_info = join([@sprintf("x%d[%d,%d)", c, lbs[c], ubs[c]) for c in 1:num_var], "\t") * "\t" * formula 

    obs, obs_info, formula
end

function obs_generator(id, seed, num_obs, noise_level, filename="", header=true)
    # Initialization
    obs = nothing
    obs_info = nothing

    obs, obs_info, formula = obsgen_from_csv("Equations.csv", id, num_obs, seed)

    if isnothing(obs)
        return nothing, nothing, nothing
    end

    # Add Gaussian noise
    if noise_level > 0
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
            if header
                DelimitedFiles.write(io, obs_info * "\n")
            end
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
                push!(ids_selected, id)
            elseif condition_number == 1 && !isexplog && num_var <= 3
                push!(ids_selected, id)
            elseif condition_number == 2 && !isexplog && num_var >= 4
                push!(ids_selected, id)
            elseif condition_number == 3 && isexplog
                push!(ids_selected, id)
            else
                nothing
            end
        end
    end
    return ids_selected
end


function main(args)
    i = 0
    dir             = args[i+=1]
    id              = parse(Int, args[i+=1])
    seed            = parse(Int, args[i+=1])
    num_obs         = parse(Int, args[i+=1])
    noise_level     = parse(Float64, args[i+=1])
    
    obs_generator(id, seed, num_obs, noise_level,"$(dir)data.obs", false)
end

main(ARGS)
## Generate obs and write optimal solutions
# DIR = "../obs/"ju
# if !isdir(DIR)
#     mkdir(DIR)
# end
# io = open(DIR * "optval.log", "w+")
# write(io, @sprintf("%8s\t%8s\t%8s\t%16s\t%16s\n", "id", "seed", "num_obs", "noise_level", "optval"))
# for id in 1:100, seed in [1], num_obs in [10], noise_level in [0]
#     INS = "i$(id)_n$(num_obs)_z$(noise_level)"

#     obs, obs_info = obs_generator(id, seed, num_obs, noise_level, "$(DIR)$(INS).obs", false)
#     for noise_level in [1e-04]
#         optval, noise_level = obs_optval(id, seed, num_obs, noise_level)
#         if !isnothing(optval)
#             write(io, @sprintf("%8d\t%8d\t%8d\t%16.1e\t%16.6e\n", id, seed, num_obs, noise_level, optval))
#         end
#     end
# end
# close(io)

# for condition_number in [0,1,2,3]
#     ids = get_ids(1:100, condition_number)
#     cnt = length(ids)
#     println("condition_number=$(condition_number), cnt=$(cnt), ", join(ids, " "))
# end