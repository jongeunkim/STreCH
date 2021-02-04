Install Julia. Available at [Julia Download](https://julialang.org/downloads/).

Install SCIPOptSuite. Available at [SCIPOptSuite Download](https://www.scipopt.org/index.php#download). 
Before you download, please check which version is supported by the SCIP package in Julia [SCIP.jl](https://github.com/scipopt/SCIP.jl).
Version 7.0.0 is the most recent supported version (confirmedd on 08/04/2020).

Execute `install_pkgs.jl` to install the required julia packages.

`example/data.txt` is a sample dataset. Each row represents a data point. The code assumes that the last column is the dependent varaible.

`main_example.jl` is an example code to run MINLP and STreCH. 

Once MINLP or STreCH is finished, it will create `df_sol.csv`, which includes the list of formulas and their objective values (MSE).

The paper is working at Overleaf.

