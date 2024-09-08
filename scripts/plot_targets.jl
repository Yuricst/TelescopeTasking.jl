"""Plot targets' orbits"""

using Colors
using ColorSchemes
using GeometryBasics
using CairoMakie
using GLPK
using Gurobi
using HiGHS
using JSON
using JuMP
using LinearAlgebra
using ProgressMeter: @showprogress
using Printf: @printf
using SatelliteToolboxTle
using SatelliteToolboxSgp4
using SatelliteToolboxTransformations

include(joinpath(@__DIR__, "../src/TelescopeTasking.jl"))

# load config jsons
target_choice = "A"

# load Earth parameters
# eop_iau1980 = fetch_iers_eop()
eop_file = joinpath(@__DIR__, "..", "data", "eop_iau1980", "finals.all.csv")
eop_iau1980 = read_iers_eop(eop_file, Val(:IAU1980))

# load TLE files
tles = read_tles(read(joinpath(@__DIR__, "..", "data", "tles", "AAS25target$(target_choice).txt"), String))
println("There are $(length(tles)) TLEs in the file")

# plot in 3D
fig = Figure(size=(500, 500))
ax3 = Axis3(fig[1,1]; aspect = (1, 1, 1), xlabel="X, km", ylabel="Y, km", zlabel="Z, km")

# propagate TLE
dt_sec = 20
for tle in tles#[1:10]
    period_sec = TelescopeTasking.tle2period(tle)
    _dts_min = (0:dt_sec:period_sec) / 60        # convert seconds to minutes
    rvs = TelescopeTasking.integrate_sgp4(
        tle,
        _dts_min,
        :eci,
        eop_iau1980,
        [0,0,0],
    )
    lines!(ax3, rvs[1,:]./1e3, rvs[2,:]./1e3, rvs[3,:]./1e3, color=:red, linewidth=0.2)
end

# # plot Earth
# mesh!(ax3, Sphere(GeometryBasics.Point3f0(0), 6378.0/1e3), alpha=1)
display(fig)
