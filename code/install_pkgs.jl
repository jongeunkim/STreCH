using Pkg

pkgs = ["CSV", "DataFrames", "DataStructures", "JuMP", "SCIP"]
Pkg.add(pkgs)

Pkg.status()
