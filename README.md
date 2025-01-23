# geothermal-modelling

## Cloning the repository

The repository is structured to contain the PyPSA-Earth model and PyPSA-USA as submodules. 

There are two ways to clone the repository (1) Recursive cloning of the repository including the submodules or (2) Cloning the main and submodules independently

(1) To recursively clone the modules, the flag `--recurse-submodules` is to be used along with the `git clone` command. The repository can be cloned either using the ssh / http clone commands

    git clone --recurse-submodules git@github.com:open-energy-transition/geothermal-modelling.git

OR 

    git clone --recurse-submodules https://github.com/open-energy-transition/geothermal-modelling.git

(2) To independently clone the modules,

Firstly, the repository is cloned using the `git clone` command. The modules can be cloned either using the ssh / http clone commands

    git clone git@github.com:open-energy-transition/geothermal-modelling.git
OR 

    git clone https://github.com/open-energy-transition/geothermal-modelling.git

and then using the `git submodule` command:

    git submodule update --init --recursive

The new commits of the submodule can be fetched and updated using the command:

    git submodule update --remote

## Solver options

The model requires a solver to be installed to perform the optimisations. The following solvers can be used with PyPSA-Eur and can be set in the config files:

- [Gurobi](https://support.gurobi.com/hc/en-us/articles/14799677517585-Getting-Started-with-Gurobi-Optimizer)
- [Cplex](https://www.ibm.com/products/ilog-cplex-optimization-studio)
- [cbc](https://github.com/coin-or/Cbc#DownloadandInstall)
- [HiGHs](https://highs.dev/)
- [GLPK](https://www.gnu.org/software/glpk/)
- [SCIP](https://scipopt.github.io/PySCIPOpt/docs/html/index.html)

Gurobi and Cplex are commercial solvers and will require licenses to run the model.


## Config settings

Config parameter                            | Initial run                            | Subsequent runs           | 
--------------------------------------------|----------------------------------------|---------------|
PyPSA-Earth related parameters                                                                                    |
--------------------------------------------|----------------------------------------|---------------|
enable -> retrieve_databundle               |true                                    |false                       |
enable -> download_osm_data                 |true                                    |false                       |
enable -> build_natura_raster               |false                                   |false                       |
enable -> retrieve_cost_data                |true                                    |true                        |
cluster_options -> alternative_clustering   |true / false                            |true / false                |
electricity -> custom_powerplants           |replace                                 |replace                     |
costs -> version                            |>= 0.10.0 (esp for sector coupled model)|>= 0.10.0                   |
--------------------------------------------|----------------------------------------|---------------|
Geothermal exclusive parameters                                                                 |
--------------------------------------------|----------------------------------------|---------------|
geothermal -> retrieve_geothermal_databundle|true                                    |false                       |
geothermal -> demand_year                   |2021 (baseline run)            |2021 (baseline run)                  |

## Running the workflow

To run the workflow use the following snakemake rule
        snakemake -call summary

This rule triggers all the rules in the geothermal repository as well as PyPSA-Earth imported rules

## Format the code

Note: Use this only if something is being committed to the repository

To format the code run the command:

    pre-commit run --all