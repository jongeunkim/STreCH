Install Julia. Available at [Julia Download](https://julialang.org/downloads/). This code is compatible with Julia 1.5.3.

Install SCIPOptSuite. Available at [SCIPOptSuite Download](https://www.scipopt.org/index.php#download). 
Before you download, please check which version is supported by the SCIP package in Julia [SCIP.jl News](https://github.com/scipopt/SCIP.jl/blob/master/NEWS.md).
The most recent version SCIP v.7.0.2 is supported (confirmed on 02/04/2021).

Execute `install_pkgs.jl` to install the required julia packages.

`example/data.txt` is a sample dataset. The formula is x1*x2+1. Each row represents a data point. The code assumes that the last column is the dependent varaible.

`main_example.jl` is an example code to run MINLP and STreCH. 

Once MINLP or STreCH is finished, it will create `df_sol.csv`, which includes the list of formulas and their objective values (MSE).

The paper is working at Overleaf.

`data/` folder contains all the datasets used in the paper. Each folder in `data/` includes a info file `info.txt`, two training set with/without noise, one validation set, and one testing set.
