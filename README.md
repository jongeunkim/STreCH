Install SCIPOptSuite and Julia. 
Refer https://docs.google.com/document/d/1DCRQAZugv181wbx-RZYTGezuMRHp14NmeFdsRL8omHY/edit?usp=sharing.

Run install_pkgs.jl to install the required julia packages.

Then, you can execute run.sh to run a test. 

The paper is working at Overleaf.

# Code Scheme

![image of code scheme](https://www.dropbox.com/s/3s2lnkafykl81hs/symbolic_regression_code_scheme.png?dl=1)

# Key functions

## obsgen.jl
- `obs, obs_info, noise_level = obs_generator(id, seed, num_obs, noise_level, filename="")`

    Return a observation. If `filename` != "", save the observation at `filename`.

- `mse, noise_level = obs_optval(id, seed, num_obs, noise_level)`

    Return the optimal value (MSE). If `noise_level` > 0, then `mse` > 0. 

## minlp.jl
- `feasible, optfeasible, time, optval, ysol, csol, vsol = solve_MINLP(nodes, obs, operators;
                                                                        print_all_solutions=false,
                                                                        formulation="New-NR",
                                                                        integer_constant=true,
                                                                        ysol=nothing, 
                                                                        ysol_dist_max=0, 
                                                                        ysol_dist_min=0, 
                                                                        ysol_fixlevel=-1,
                                                                        csol=nothing,
                                                                        scip_verblevel=3,
                                                                        scip_time=10.0,
                                                                        scip_absgap=0.0,
                                                                        scip_nodes=-1,
                                                                        scip_presolve=-1)`

    Solve a MINLP for symbolic regression

    Mandatory arguments are nodes, observations, and operators.
        
        nodes:      set of nodes that we consider        
        obs:        observations in array, assume that the last column is the dependent variable.        
        operators:  set of operators that we consider
    
    Optional arguments include formulation, SCIP parameters, and y-solution for local branching.

        print_all_solutions: print all solutions found in the log file (both optimal and nonoptimal) 

        formulation:    "Cozad":    Cozad's formulation with no optional cut such as redundancy elimination and no symmetry breacking
                        "New":      New formulation with no optional cut
                        "-CR":      Add Cozad's redundancy elimination cut
                        "-NR":      Add new redundancy elimination cut
            * I have not implemented symmetry breaking and implication constarint yet.

        integer_constant: constnat is an integer or not, true/false

        (matters if ysol != nothing)
        ysol:                           solution of variable y in dict format (ysol[`nodeid`] = `operator`)
        ysol_dist_min, ysol_dist_max:   min/max distance from the current solution y
        ysol_fixlevel:                  if value is 0, fix all y solutions,
                                        if value is 1, fix all nonleaves,
                                        ... as value increases, fix less number of nodes 
        
        (matters if csol != nothing)
        csol:   solution of variable c in dict format (csol[`nodeid`] = `constant value`)
                if it is given, we fix csol

        scip_*: scip parameters 

    It returns seven items.
    - `feasible`: true if there is a feasible solution, otherwise false
    - `optfeasible`: true if SCIP's optimal solution is feasible, otherwise false (it may be infeasible because of numerical error)
    - `time`: elapsed time read from the solver (SCIP)
    - `optval`: optimal value, if the solver's optimal solution is not feasible then return the best solution's value
    - `ysol, csol, vsol`: optimal solution of y, c, and v. Keep in mind that this may not be feasible (check `optfeasible`)


## heur.jl
`arr_obj, arr_time, arr_active = solve_Heuristic(obs, operators; 
                                                    time_limit=60, 
                                                    obj_termination=1e-06,
                                                    init_solve=1, 
                                                    subsampling=1,
                                                    max_depth=4, 
                                                    min_improvement=0.01)`
    
Solve a local branching heuristic.

`obs, operators` are mandatory inputs. Optional inputs are the followings.
- `time_limit`: time limit
- `obj_termination`: heuristic ends when it finds a solution whose objective value is lower or equal to `obj_termination`
- `init_solve`: heuristic solves depth-`init_solve` problem as an initial problem 
- `subsampling`: heuristic solves with `100 * subsampling` % of high error observations at every iteration
- `max_depth`: a solution expression tree is limited by depth `max_depth`
- `min_improvement`: at every iteration, the solver stops when it finds a solution whose objective value is `100 * min_improvement` % lower than the current one

It returns three items, arrays of objective values, cumulative times, and # of active nodes.
`arr_obj[i], arr_time[i], arr_active[i]` correspond to the result of the i-th iteration.


 
