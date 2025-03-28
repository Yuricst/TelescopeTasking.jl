# Scripts used to run optimization instances and analyze/plot the results

## Primary directory

- `config` contains problem config JSON files
- `plots` contains plots generated during the analysis
- `solutuion_XYZ` includes solution files; `solution_JASS` are the solutions used in the analysis of the paper, with a max slew rate of 2 deg/s. `solution_JASS_slow` are solutions with a max slew rate of 0.5 deg/s.

## Scripts for plot in JASS manuscript

- Fig 2: `plot_sparsity_STTP.jl`
- Fig 3: `plot_targets.jl`
- Fig 4: `plot_targets_starlink.jl`
- Fig 5/6: `read_gurobi_log.py`
- Fig 7: `plot_STTP_slew.jl`
- Fig 8/9: `plot_MTTP_solution_passes.jl`
- Fig 10: `tabulate_results.jl`