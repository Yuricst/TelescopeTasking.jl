module TelescopeScheduling

using Interpolations
using JuMP
using LinearAlgebra
using Printf: @printf
using SatelliteToolboxTle

include("pass.jl")
include("problem.jl")

end # module TelescopeScheduling
