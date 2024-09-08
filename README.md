# TelescopeTasking.jl
Deterministic operation scheduling problem for nightly observation

### Quick start

This repository is developed on Julia v1.10. Download from: [https://julialang.org/downloads/](https://julialang.org/downloads/). 
Once julia is setup, `git clone` and `cd` to the repository, start the julia-repl, then run

```julia
julia> ]
(@v1.10) pkg> activate .
(TelescopeTasking) pkg>
```

Then, run tests to check installation is successful:

```julia
julia> include("tests/runtest.jl")
```

### Dependencies

Dependencies are listed in `Project.toml`. 
Part of the package dependencies are subsets of the [SatelliteToolbox.jl](https://juliahub.com/ui/Packages/General/SatelliteToolbox) developed by the [Brazilian National Institute for Space Research (INPE)](https://www.gov.br/inpe/pt-br), namely:

- [SatelliteToolboxTle.jl](https://github.com/JuliaSpace/SatelliteToolboxTle.jl) for handling TLEs
- [SatelliteToolboxSgp4.jl](https://github.com/JuliaSpace/SatelliteToolboxSgp4.jl) for propagation
- [SatelliteToolboxTransformations.jl](https://github.com/JuliaSpace/SatelliteToolboxTransformations.jl) for frame transformations

The `TLE` object's structure can be found in the [SatelliteToolboxTle.jl docs page](https://juliaspace.github.io/SatelliteToolboxTle.jl/stable/man/tle_structure/). 


### Scripts for AAS 2025

Scripts for AAS 2025 are located in folder `./scripts`;

- `scripts/solve_STTP.jl` : solve an instance of the STTP with $E = 1,2,3$ with chosen target set
- `scripts/solve_MTTP.jl` : solve an instance of the MTTP with $E = 1,2,3$ with chosen target set

