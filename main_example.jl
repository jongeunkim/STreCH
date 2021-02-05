# Copyright (c) <2021>, <Sven Leyffer & Jongeun Kim>
# All rights reserved.

# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree. 

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