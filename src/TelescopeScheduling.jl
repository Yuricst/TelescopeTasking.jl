module TelescopeScheduling

using GLMakie
using JuMP
using Interpolations
using LinearAlgebra
using Printf: @printf
using SatelliteToolboxSgp4
using SatelliteToolboxTle
using SatelliteToolboxTransformations

include("transformations.jl")
include("pass.jl")
include("process_tles.jl")
include("problem.jl")
include("visualizations.jl")

end # module TelescopeScheduling
