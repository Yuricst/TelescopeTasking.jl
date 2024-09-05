"""Developing methods to save instances"""

using Colors
using ColorSchemes
using GeometryBasics
using GLMakie
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

include(joinpath(@__DIR__, "../src/TelescopeScheduling.jl"))

# load Earth parameters
eop_file = joinpath(@__DIR__, "..", "data", "eop_iau1980", "finals.all.csv")
eop_iau1980 = read_iers_eop(eop_file, Val(:IAU1980))

# load tles used
tles = read_tles(read(joinpath(@__DIR__, "tles_used.txt"), String))

# load solutions
solution_dict = JSON.parsefile(joinpath(@__DIR__, "test_solution.json"))
jd0_obs = solution_dict["jd0_obs"]
obs_duration = solution_dict["obs_duration"]
min_elevation = solution_dict["min_elevation"]
min_obs_duration = solution_dict["min_obs_duration"]
exposure_duration = solution_dict["exposure_duration"]
observer_lla = solution_dict["observer_lla"]

passes = [TelescopeScheduling.VisiblePass(pass_dict) for pass_dict in solution_dict["passes_dict"]]
# recreate passes for making sure it is correctly ordered
passes2, _ = TelescopeScheduling.tles_to_passes(
    tles,
    eop_iau1980,
    jd0_obs,
    obs_duration,
    min_elevation,
    min_obs_duration,
    exposure_duration,
    observer_lla,
    dt_sec=10,
)
for (p1,p2) in zip(passes, passes2)
    @assert p1.tle == p2.tle
    @assert p1.times == p2.times
    @assert p1.azimuths == p2.azimuths
    @assert p1.elevations == p2.elevations
end

# get selected passes from solution 
Y = [el > 0.5 for el in solution_dict["Y"]]
selected_passes = [pass for (pass, y) in zip(passes, value.(Y)) if y > 0.5]


# plot of selected passes
fig_sol = Figure(size=(1000,500))
ax_sol = PolarAxis(fig_sol[1:2,1])
TelescopeScheduling.polar_plot_passes!(ax_sol, passes; color=:grey, linewidth=1.0)
TelescopeScheduling.polar_plot_passes!(ax_sol, selected_passes; 
    linewidth=1.5, color_by_target=true, exposure_only=true)

# plot time-history
axes = [Axis(fig_sol[1,2]; xlabel="Time, min", ylabel="Azimuth, deg"),
        Axis(fig_sol[2,2]; xlabel="Time, min", ylabel="Elevation, deg")]
TelescopeScheduling.plot_time_history!(axes, passes; jd_ref=jd0_obs, color=:grey, linewidth=1.0)
TelescopeScheduling.plot_time_history!(axes, selected_passes; jd_ref=jd0_obs, 
    linewidth=1.5, color_by_target=true, exposure_only=true)
#suptitle("E = $num_exposure")
save(joinpath(@__DIR__, "solution_passes.png"), fig_sol)

display(fig_sol)