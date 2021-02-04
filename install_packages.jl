using Pkg

function install_packages()
    pkgs = ["CSV", "DataFrames", "DataStructures", "Glob", "JuMP", "MathOptInterface", "SCIP"]
    Pkg.add(pkgs)
    Pkg.status()
end

install_packages()
