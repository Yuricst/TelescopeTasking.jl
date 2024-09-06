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

# load Earth parameters
# eop_iau1980 = fetch_iers_eop()
eop_file = joinpath(@__DIR__, "..", "data", "eop_iau1980", "finals.all.csv")
eop_iau1980 = read_iers_eop(eop_file, Val(:IAU1980))

# load config jsons
telescope = JSON.parsefile(joinpath(@__DIR__, "configs/config_telescope.json"))
config = JSON.parsefile(joinpath(@__DIR__, "configs/config_STTP.json"))
target_choice = "C"
num_exposure = 1
solver_choice = "Gurobi"

# choose solver
if solver_choice == "Gurobi"
    solver = MOI.OptimizerWithAttributes(
        Gurobi.Optimizer,
        "TimeLimit" => 1200,
    )
elseif solver_choice == "GLPK"
    solver = GLPK.Optimizer
elseif solver_choice == "HiGHS"
    solver = HiGHS.Optimizer
else
    error("Solver choice $solver_choice not recognized!")
end

# load TLE files
tles = read_tles(read(joinpath(@__DIR__, "..", "data", "tles", "AAS25target$(target_choice).txt"), String))
println("There are $(length(tles)) TLEs in the file")

# initial epoch of local nightfall
jd0_ref = telescope["jd0_ref"]
@assert maximum([tle_epoch(tle) for tle in tles]) <= jd0_ref "TLEs are later than reference JD!"
jd0_obs = jd0_ref + 0.2

# get passes
slew_rate = deg2rad(telescope["slew_rate"])                         # rad/s
buffer_times = [telescope["buffer_t0"], telescope["buffer_t1"]]     # times in seconds
obs_duration = 8 * 3600             # in seconds
min_elevation = deg2rad(telescope["min_elevation"] )            # in radians
min_obs_duration = telescope["min_obs_duration"]                # in seconds
exposure_duration = telescope["exposure_duration"]              # in seconds
observer_lat = deg2rad(config["observer"]["latitude"])          # degrees --> radians
observer_lon = deg2rad(config["observer"]["longitude"])         # degrees --> radians
observer_alt = deg2rad(config["observer"]["altitude"])          # meters
observer_lla = [observer_lat, observer_lon, observer_alt]


# iterate through num_exposure
num_exposures = [1, 2, 3]
for num_exposure in num_exposures
    experiment_name = config["name"] * "_target$(target_choice)_E$(num_exposure)"
    println(" *************** Experiment name: $experiment_name *************** ")
    passes, _ = TelescopeTasking.tles_to_passes(
        tles,
        eop_iau1980,
        jd0_obs,
        obs_duration,
        min_elevation,
        min_obs_duration,
        exposure_duration,
        observer_lla,
        dt_sec = 10,
        num_exposure = num_exposure,
    )
    @printf("Detected %d passes\n", length(passes))
    if length(passes) == 0
        @printf("No passes detected, skipping experiment with E = %d\n", num_exposure)
        continue
    end

    # construct problem
    problem = TelescopeTasking.TelescopeTaskingProblem(
        passes, num_exposure,
        slew_rate;
        buffer_times = buffer_times
    )
    @printf("Detected %d observable targets with %d exposures\n", problem.m, num_exposure)

    # solve problem
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
    open(joinpath(@__DIR__, "solutions", "solution_$(experiment_name)_$(solver_choice).json"), "w") do io
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
    save(joinpath(@__DIR__, "plots", "passes_$(experiment_name)_$(solver_choice).png"), fig_sol)

    # plot sparsity of T and A
    fig_spy = Figure(size=(800,400))
    ax_spy_A = Axis(fig_spy[1,1]; title="A matrix", xlabel="Target", ylabel="Pass", yreversed=true)
    spy!(ax_spy_A, transpose(problem.A))
    ax_spy_T = Axis(fig_spy[1,2]; title="T matrix", xlabel="Pass j", ylabel="Pass i", yreversed=true)
    spy!(ax_spy_T, transpose(problem.T))
end