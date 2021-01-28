using Pkg

function install_packages()
    pkgs = ["CSV", "DataFrames", "DataStructures", "JuMP", "SCIP"]
    Pkg.add(pkgs)
    Pkg.status()
end

install_packages()
