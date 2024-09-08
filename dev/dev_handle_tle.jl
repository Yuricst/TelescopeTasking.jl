"""
Handling TLEs and solving optimization problem
"""

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


# initial epoch of local nightfall
jd0_obs = 2.46055755221534e6 + 0.45     # in julian date

# criteria for valid passes
min_elevation = deg2rad(15)
min_duration_day = 60/86400        # minimum duration: 1 minute, in days
exposure_duration_day = 45/86400   # exposure duration: 45 seconds, in days

# load Earth parameters
# eop_iau1980 = fetch_iers_eop()
eop_file = joinpath(@__DIR__, "..", "data", "eop_iau1980", "finals.all.csv")
eop_iau1980 = read_iers_eop(eop_file, Val(:IAU1980))

# load TLE files
# path_to_tles = joinpath(@__DIR__, "..", "data", "tles", "iridium-33-debris.txt")
path_to_tles = joinpath(@__DIR__, "..", "data", "tles", "active.txt")
tles_str = read(path_to_tles, String)

# convert to TLE objects
tles = read_tles(tles_str)
@show length(tles)

# filter them
names_include = ["STARLINK",]#, "GLOBALSTAR", "IRIDIUM"]
tles = TelescopeTasking.filter(tles, names_include = names_include)
@show length(tles)
tles = tles[1:1000]          # only use subset of them

# save to file tles used
TelescopeTasking.tles_to_file(tles, joinpath(@__DIR__, "tles_used.txt"))

# get passes
obs_duration = 8 * 3600             # in seconds
min_elevation = deg2rad(30)
min_obs_duration = 100              # in seconds
exposure_duration = 60              # in seconds
observer_lat = deg2rad(45)
observer_lon = deg2rad(100)
observer_alt = 30.0
observer_lla = [observer_lat, observer_lon, observer_alt]

passes, sph_ENU_list = TelescopeTasking.tles_to_passes(
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
@show length(passes)
smas = [TelescopeTasking.tle2sma(pass.tle) for pass in passes]       # for sanity check

# construct problem
num_exposure = 2
slew_rate = deg2rad(2)      # rad/s
buffer_times = [15, 0]      # times in seconds
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
solution_dict = TelescopeTasking.STTP_solution_to_dict(
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
open(joinpath(@__DIR__, "test_solution.json"), "w") do io
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
save(joinpath(@__DIR__, "solution_passes.png"), fig_sol)

# plot sparsity of T and A
fig_spy = Figure(size=(800,400))
ax_spy_A = Axis(fig_spy[1,1]; title="A matrix", xlabel="Target", ylabel="Pass", yreversed=true)
spy!(ax_spy_A, transpose(problem.A))
ax_spy_T = Axis(fig_spy[1,2]; title="T matrix", xlabel="Pass j", ylabel="Pass i", yreversed=true)
spy!(ax_spy_T, transpose(problem.T))

display(fig_sol)