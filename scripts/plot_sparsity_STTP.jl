"""Script to solve Single Telescope Tasking Problem (STTP)"""

using CairoMakie
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

CairoMakie.activate!()

include(joinpath(@__DIR__, "../src/TelescopeTasking.jl"))

# load Earth parameters
# eop_iau1980 = fetch_iers_eop()
eop_file = joinpath(@__DIR__, "..", "data", "eop_iau1980", "finals.all.csv")
eop_iau1980 = read_iers_eop(eop_file, Val(:IAU1980))

config_filenames = [
    "config_STTP1.json",
    # "config_STTP2.json",
    # "config_STTP3.json",
]

for config_filename in config_filenames
    # load config jsons
    config_telescope = JSON.parsefile(joinpath(@__DIR__, "configs/config_telescope.json"))
    config = JSON.parsefile(joinpath(@__DIR__, "configs", config_filename))
    target_choice = "A" #"S1"
    solver_choice = "Gurobi"
    num_exposures = [1, 2, 3]

    # choose solver
    if solver_choice == "Gurobi"
        solver = MOI.OptimizerWithAttributes(
            Gurobi.Optimizer,
            "TimeLimit" => 600,
            "Method" => 2,
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

    # get passes
    slew_rate = deg2rad(config_telescope["slew_rate"])                         # rad/s
    buffer_times = [config_telescope["buffer_t0"], config_telescope["buffer_t1"]]     # times in seconds
    min_elevation = deg2rad(config_telescope["min_elevation"] )            # in radians
    min_obs_duration = config_telescope["min_obs_duration"]                # in seconds
    exposure_duration = config_telescope["exposure_duration"]              # in seconds
    observer_lat = deg2rad(config["observer"]["latitude"])          # degrees --> radians
    observer_lon = deg2rad(config["observer"]["longitude"])         # degrees --> radians
    observer_alt = config["observer"]["altitude"]                   # meters
    observer_lla = [observer_lat, observer_lon, observer_alt]

    # initial epoch of local nightfall
    jd0_ref = config_telescope["jd0_ref"]
    @assert maximum([tle_epoch(tle) for tle in tles]) <= jd0_ref "TLEs are later than reference JD!"
    jds_night = TelescopeTasking.earliest_night(jd0_ref, observer_lla, eop_iau1980)
    jd0_obs = jds_night[1]
    obs_duration = 86400 * (jds_night[2] - jds_night[1])

    # iterate through num_exposure
    for num_exposure in num_exposures[1:1]
        _experiment_name = config["name"] * "_target$(target_choice)_E$(num_exposure)"
        println(" *************** Experiment name: $_experiment_name *************** ")
        times_measure = [time(),]
        _passes, _ = TelescopeTasking.tles_to_passes(
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
        push!(times_measure, time())
        @printf("Detected %d passes; interval time: %1.4f sec\n",
            length(_passes), times_measure[end] - times_measure[end-1])
        if length(_passes) == 0
            @printf("No passes detected, skipping experiment with E = %d\n", num_exposure)
            continue
        end

        # construct problem
        # slew_rate = deg2rad(0.05)         # overwrite for manuscript
        problem = TelescopeTasking.TelescopeTaskingProblem(
            _passes, num_exposure,
            slew_rate;
            buffer_times = buffer_times
        )
        
        # plot T
        fontsize = 24
        fig_spy_T = Figure(size=(400,400))
        ax_spy_T = Axis(fig_spy_T[1,1];
            xlabelsize = fontsize,
            ylabelsize = fontsize,
            xticklabelsize = fontsize-2,
            yticklabelsize = fontsize-2,
            xlabel=L"Pass $j$", 
            ylabel=L"Pass $i$", 
            yreversed=true)
        spy!(ax_spy_T, transpose(problem.T))
        save(joinpath(@__DIR__, "plots",
            "spy_T_$_experiment_name.png"), fig_spy_T)

        fig_spy_A = Figure(size=(400,400))
        ax_spy_A = Axis(fig_spy_A[1,1];
            xlabelsize = fontsize,
            ylabelsize = fontsize,
            xticklabelsize = fontsize-2,
            yticklabelsize = fontsize-2,
            # xlabel=L"Pass $i$", 
            # ylabel=L"Target $k$", 
            xlabel=L"Target $k$", 
            ylabel=L"Pass $i$", 
            yreversed=true)
        #spy!(ax_spy_A, problem.A) #transpose(problem.A))
        spy!(ax_spy_A, transpose(problem.A))
        save(joinpath(@__DIR__, "plots",
            "spy_A_$_experiment_name.png"), fig_spy_A)
    end
end