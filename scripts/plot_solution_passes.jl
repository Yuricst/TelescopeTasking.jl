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

include(joinpath(@__DIR__, "../src/TelescopeTasking.jl"))

# load Earth parameters
eop_file = joinpath(@__DIR__, "..", "data", "eop_iau1980", "finals.all.csv")
eop_iau1980 = read_iers_eop(eop_file, Val(:IAU1980))


# load config jsons
telescope = JSON.parsefile(joinpath(@__DIR__, "configs/config_telescope.json"))
config = JSON.parsefile(joinpath(@__DIR__, "configs/config_STTP.json"))
target_choice = "A"
num_exposure = 1
solver_choice = "Gurobi"
experiment_name = config["name"] * "_target$(target_choice)_E$(num_exposure)"
@printf("Plotting pass from experiment %s\n", experiment_name)

# load tles used
tles = read_tles(read(joinpath(@__DIR__, "../data/tles/AAS25target$(target_choice).txt"), String))

# load solutions
solution_dict = JSON.parsefile(joinpath(@__DIR__, "solutions/solution_$(experiment_name)_$(solver_choice).json"))
jd0_obs = solution_dict["jd0_obs"]
obs_duration = solution_dict["obs_duration"]
min_elevation = solution_dict["min_elevation"]
min_obs_duration = solution_dict["min_obs_duration"]
exposure_duration = solution_dict["exposure_duration"]
observer_lla = solution_dict["observer_lla"]

passes = [TelescopeTasking.VisiblePass(pass_dict) for pass_dict in solution_dict["passes_dict"]]
# recreate passes for making sure it is correctly ordered
passes2, _ = TelescopeTasking.tles_to_passes(
    tles,
    eop_iau1980,
    jd0_obs,
    obs_duration,
    min_elevation,
    min_obs_duration,
    exposure_duration,
    observer_lla,
    dt_sec=10,
    num_exposure = num_exposure,
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
fontsize = 24
fig_sol = Figure(size=(500,500);)
ax_sol = PolarAxis(fig_sol[1:1,1];
    rticks = ([0,20,40,60], ["90","70","50","30"]),
    # rlabelsize=fontsize,
    # thetalabelsize=fontsize,
    rticklabelsize=fontsize-1,
    thetaticklabelsize=fontsize-1,)
TelescopeTasking.polar_plot_passes!(ax_sol, passes; color=:grey50, linewidth=0.5)
TelescopeTasking.polar_plot_passes!(ax_sol, selected_passes; 
    linewidth=1.5, color_by_target=true, exposure_only=true)

# # plot time-history
# axes = [Axis(fig_sol[1,2]; xlabel="Time, min", ylabel="Azimuth, deg"),
#         Axis(fig_sol[2,2]; xlabel="Time, min", ylabel="Elevation, deg")]
# TelescopeTasking.plot_time_history!(axes, passes; jd_ref=jd0_obs, color=:grey, linewidth=1.0)
# TelescopeTasking.plot_time_history!(axes, selected_passes; jd_ref=jd0_obs, 
#     linewidth=1.5, color_by_target=true, exposure_only=true)
save(joinpath(@__DIR__, "plots", "passes_$(experiment_name)_$(solver_choice).png"), fig_sol)


# plot of non-exposure distance between passes
fig_L = Figure(size=(500,500))
ax_L = Axis(fig_L[1,1]; xlabel="Pass index", ylabel="Non-exposure distance, deg")
_, Ls = TelescopeTasking.get_nonexposure_distance(passes, Y)
Ls_nonzero = [L for L in Ls if L > 0]
scatter!(ax_L, 1:length(Ls_nonzero), rad2deg.(Ls_nonzero); color=:grey10, markersize=5)

# display figure
display(fig_L)