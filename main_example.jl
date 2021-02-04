include("src/header.jl")
include("src/util.jl")
include("src/MINLP.jl")
include("src/STreCH.jl")


"""
In `example` directory, there are `data.txt` and `params.txt` files.
The solver will read those two files and solve.
"""

MINLP("example/")
# STreCH("example/")