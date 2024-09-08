"""
Developing multi-telescope tasking problem
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
jd0_ref = 2.46055755221534e6 + 0.45     # in julian date

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
tles = tles[1:200]          # only use subset of them

# get passes
obs_duration_1 = 8 * 3600             # in seconds
min_elevation = deg2rad(30)
min_obs_duration = 100              # in seconds
exposure_duration = 60              # in seconds

observer_lla_1 = [deg2rad(45), deg2rad(13), 30.0]
jds_night_1 = TelescopeTasking.earliest_night(jd0_ref, observer_lla_1, eop_iau1980)
jd0_obs_1 = jds_night_1[1]
obs_duration_1 = 86400 * (jds_night_1[2] - jds_night_1[1])
@assert jd0_ref <= jd0_obs_1

observer_lla_2 = [deg2rad(21), deg2rad(102), 30.0]
jds_night_2 = TelescopeTasking.earliest_night(jd0_ref, observer_lla_2, eop_iau1980)
jd0_obs_2 = jds_night_2[1]
obs_duration_2 = 86400 * (jds_night_2[2] - jds_night_2[1])
@assert jd0_ref <= jd0_obs_2

passes_1, sph_ENU_list = TelescopeTasking.tles_to_passes(
    tles, eop_iau1980, jd0_obs_1, obs_duration_1,
    min_elevation, min_obs_duration, exposure_duration, observer_lla_1;
)
passes_2, sph_ENU_list = TelescopeTasking.tles_to_passes(
    tles, eop_iau1980, jd0_obs_2, obs_duration_2,
    min_elevation, min_obs_duration, exposure_duration, observer_lla_2;
)
passes_per_telescope = [passes_1, passes_2]
@printf("Pass from observer 1: %d\n", length(passes_1))
@printf("Pass from observer 2: %d\n", length(passes_2))
@printf("Total passes: %d\n", length(passes_1) + length(passes_2))

# construct problem
num_exposure = 1
slew_rate = deg2rad(2)      # rad/s
buffer_times = [15, 0]      # times in seconds
problem = TelescopeTasking.MultiTelescopeTaskingProblem(
    passes_per_telescope,
    num_exposure,
    slew_rate;
    buffer_times = buffer_times
)
@show problem;

# solve problem
solver = MOI.OptimizerWithAttributes(Gurobi.Optimizer,
    "TimeLimit" => 1200)
# solver = HiGHS.Optimizer
X, Y, Y_per_telescope, status = TelescopeTasking.solve!(problem, solver)
selected_passes_per_telescope = [
    [pass for (pass, y) in zip(passes, value.(Y)) if y > 0.5]
    for (passes,Y) in zip(passes_per_telescope, Y_per_telescope)
]

# save to dictionary
observer_lla_per_telescope = [observer_lla_1, observer_lla_2]
obs_duration_per_telescope = [obs_duration_1, obs_duration_2]
solution_dict = TelescopeTasking.MTTP_solution_to_dict(
    problem,
    passes_per_telescope,
    jd0_ref,
    obs_duration_per_telescope,
    min_elevation,
    min_obs_duration,
    exposure_duration,
    observer_lla_per_telescope,
    X,
    Y_per_telescope,
)

# plot of selected passes
fig_sol = Figure(size=(1400,800))
ax_polar1 = PolarAxis(fig_sol[1,1])
TelescopeTasking.polar_plot_passes!(ax_polar1, passes_1; color=:grey, linewidth=0.3)
TelescopeTasking.polar_plot_passes!(ax_polar1, selected_passes_per_telescope[1]; 
    linewidth=1.5, color_by_target=true, exposure_only=true)

    ax_polar2 = PolarAxis(fig_sol[2,1])
TelescopeTasking.polar_plot_passes!(ax_polar2, passes_2; color=:grey, linewidth=0.3)
TelescopeTasking.polar_plot_passes!(ax_polar2, selected_passes_per_telescope[2]; 
    linewidth=1.5, color_by_target=true, exposure_only=true)
    
# plot time-history
axes = [Axis(fig_sol[1,2]; xlabel="Time, hour", ylabel="Azimuth, deg"),
        Axis(fig_sol[1,3]; xlabel="Time, hour", ylabel="Elevation, deg")]
TelescopeTasking.plot_time_history!(axes, passes_1; jd_ref=jd0_ref, color=:grey, linewidth=0.3)
TelescopeTasking.plot_time_history!(axes, selected_passes_per_telescope[1]; 
    jd_ref=jd0_ref,  linewidth=1.5, color_by_target=true, exposure_only=true)


axes = [Axis(fig_sol[2,2]; xlabel="Time, hour", ylabel="Azimuth, deg"),
        Axis(fig_sol[2,3]; xlabel="Time, hour", ylabel="Elevation, deg")]
TelescopeTasking.plot_time_history!(axes, passes_2; jd_ref=jd0_ref, color=:grey, linewidth=0.3)
TelescopeTasking.plot_time_history!(axes, selected_passes_per_telescope[2]; 
    jd_ref=jd0_ref,  linewidth=1.5, color_by_target=true, exposure_only=true)

save(joinpath(@__DIR__, "MTTP_solution_passes.png"), fig_sol)
display(fig_sol)
println("Done!")
