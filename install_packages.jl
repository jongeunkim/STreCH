# Copyright (c) <2021>, <Sven Leyffer & Jongeun Kim>
# All rights reserved.

# This source code is licensed under the BSD-style license found in the
# LICENSE file in the root directory of this source tree. 

using Pkg

function install_packages()
    pkgs = ["CSV", "DataFrames", "DataStructures", "Glob", "JuMP", "MathOptInterface", "SCIP"]
    Pkg.add(pkgs)
    Pkg.status()
end

install_packages()
