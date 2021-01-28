include("src/header.jl")
include("src/util.jl")
include("src/MINLP.jl")
include("src/STreCH.jl")

# Binary:   + - * D
# Unary:    R (sqrt), E (exp), L (log), S(arccos(sin(x)) = sqrt(1-x^2))
# Constant: C (general), P (pi)

function example()
    dir = "example/"
    datafile = "data.txt"
    operators = "+-*DC"
    maxdepth = 2
    timelimit = 3600
    formulation = "Imp-F"

    MINLP(dir, datafile, operators, maxdepth, timelimit, formulation)
end

example()