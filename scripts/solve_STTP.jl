"""Script to solve Single Telescope Tasking Problem (STTP)"""

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

# load config json
config = JSON.parsefile(joinpath(@__DIR__, "configs/config_STTP.json"))
target_choice = "B"
experiment_name = config["name"] * "_target$(target_choice)"

# initial epoch of local nightfall
jd0_obs = 2.46055755221534e6 + 0.45     # in julian date

# load Earth parameters
# eop_iau1980 = fetch_iers_eop()
eop_file = joinpath(@__DIR__, "..", "data", "eop_iau1980", "finals.all.csv")
eop_iau1980 = read_iers_eop(eop_file, Val(:IAU1980))

# load TLE files
tles = read_tles(read(joinpath(@__DIR__, "..", "data", "tles", "AAS25target$(target_choice).txt"), String))
@show length(tles)

# get passes
obs_duration = 8 * 3600             # in seconds
min_elevation = deg2rad(30)
min_obs_duration = 60               # in seconds
exposure_duration = 45              # in seconds
observer_lat = deg2rad(config["observer"]["latitude"])          # degrees --> radians
observer_lon = deg2rad(config["observer"]["longitude"])         # degrees --> radians
observer_alt = deg2rad(config["observer"]["altitude"])          # meters
observer_lla = [observer_lat, observer_lon, observer_alt]

passes, _ = TelescopeTasking.tles_to_passes(
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
@printf("Detected %d passes\n", length(passes))

# construct problem
num_exposure = 1
slew_rate = deg2rad(2)      # rad/s
buffer_times = [15, 5]      # times in seconds
problem = TelescopeTasking.TelescopeTaskingProblem(
    passes, num_exposure,
    slew_rate;
    buffer_times = buffer_times
)
@show problem;

# solve problem
solver = MOI.OptimizerWithAttributes(Gurobi.Optimizer,
    "TimeLimit" => 1200)
# solver = HiGHS.Optimizer
X, Y, status = TelescopeTasking.solve!(problem, solver)
selected_passes = [pass for (pass, y) in zip(passes, value.(Y)) if y > 0.5]

# save to dictionary
solution_dict = TelescopeTasking.solution_to_dict(
    passes,
    jd0_obs,
    obs_duration,
    min_elevation,
    min_obs_duration,
    exposure_duration,
    observer_lla,
    X,
    Y,
)
open(joinpath(@__DIR__, "solutions/test_solution.json"), "w") do io
    write(io, JSON.json(solution_dict))
end


# plot of selected passes
fig_sol = Figure(size=(1000,500))
ax_sol = PolarAxis(fig_sol[1:2,1])
TelescopeTasking.polar_plot_passes!(ax_sol, passes; color=:grey, linewidth=0.3)
TelescopeTasking.polar_plot_passes!(ax_sol, selected_passes; 
    linewidth=1.5, color_by_target=true, exposure_only=true)

# plot time-history
axes = [Axis(fig_sol[1,2]; xlabel="Time, min", ylabel="Azimuth, deg"),
        Axis(fig_sol[2,2]; xlabel="Time, min", ylabel="Elevation, deg")]
TelescopeTasking.plot_time_history!(axes, passes; jd_ref=jd0_obs, color=:grey, linewidth=0.3)
TelescopeTasking.plot_time_history!(axes, selected_passes; jd_ref=jd0_obs, 
    linewidth=1.5, color_by_target=true, exposure_only=true)
save(joinpath(@__DIR__, "plots/solution_passes.png"), fig_sol)

# plot sparsity of T and A
fig_spy = Figure(size=(800,400))
ax_spy_A = Axis(fig_spy[1,1]; title="A matrix", xlabel="Target", ylabel="Pass", yreversed=true)
spy!(ax_spy_A, transpose(problem.A))
ax_spy_T = Axis(fig_spy[1,2]; title="T matrix", xlabel="Pass j", ylabel="Pass i", yreversed=true)
spy!(ax_spy_T, transpose(problem.T))

display(fig_sol)