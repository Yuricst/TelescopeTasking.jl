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


"""
`target_choice` must be from "A", "S1", or "S2
"""
function main(;target_choice = "A")
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
        config_telescope = JSON.parsefile(joinpath(@__DIR__, "configs/config_telescope_slow.json"))
        config = JSON.parsefile(joinpath(@__DIR__, "configs", config_filename))
        solver_choice = "Gurobi"    # Gurobi or GPLK or HiGHS
        num_exposures = [1, 2, 3]
        time_limit = 3600           # 600 or 3600, in seconds
        save_dir = "solutions_JASS_slow"

        # iterate through num_exposure
        for num_exposure in num_exposures
            _experiment_name = config["name"] * "_target$(target_choice)_E$(num_exposure)"
            println(" *************** Experiment name: $_experiment_name *************** ")

            filepath_log = "log_$(_experiment_name)_$(solver_choice).log"
            filepath_solution = "solution_$(_experiment_name)_$(solver_choice).json"
            filepath_stats = "solve_stats_$(_experiment_name)_$(solver_choice).json"

            # choose solver
            if solver_choice == "Gurobi"
                solver = MOI.OptimizerWithAttributes(
                    Gurobi.Optimizer,
                    "WorkLimit" => time_limit,
                    "LogFile" => joinpath(@__DIR__, save_dir, "logs", filepath_log),
                    "Method" => 0,
                )
            elseif solver_choice == "GLPK"
                solver = GLPK.Optimizer
            elseif solver_choice == "HiGHS"
                solver = MOI.OptimizerWithAttributes(
                    HiGHS.Optimizer,
                    "TimeLimit" => time_limit,      # FIXME syntax?
                )
            else
                error("Solver choice $solver_choice not recognized!")
            end

            # load TLE files
            tles = read_tles(read(joinpath(@__DIR__, "..", "data", "tles", "AAS25target$(target_choice).txt"), String))
            println("There are $(length(tles)) TLEs in the file")

            # get passes
            @show config_telescope["slew_rate"]
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

        # # iterate through num_exposure
        # for num_exposure in num_exposures
        #     _experiment_name = config["name"] * "_target$(target_choice)_E$(num_exposure)"
        #     println(" *************** Experiment name: $_experiment_name *************** ")
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
            problem = TelescopeTasking.TelescopeTaskingProblem(
                _passes, num_exposure,
                slew_rate;
                buffer_times = buffer_times
            )
            push!(times_measure, time())
            @printf("Constructed problem - Detected %d observable targets with %d exposure; interval time: %1.4f sec\n",
                problem.m, num_exposure, times_measure[end] - times_measure[end-1])

            # solve problem
            _X, _Y, _solve_stats_dict = TelescopeTasking.solve(problem, solver; verbose = true)
            push!(times_measure, time())
            @printf("Finished solve; interval time: %1.4f sec\n", times_measure[end] - times_measure[end-1])

            # save to dictionary
            _solution_dict = TelescopeTasking.STTP_solution_to_dict(
                _passes,
                jd0_obs,
                obs_duration,
                min_elevation,
                min_obs_duration,
                exposure_duration,
                observer_lla,
                _X,
                _Y,
            )
            open(joinpath(@__DIR__, save_dir, filepath_solution), "w") do io
                write(io, JSON.json(_solution_dict))
            end

            open(joinpath(@__DIR__, save_dir, filepath_stats), "w") do io
                write(io, JSON.json(_solve_stats_dict))
            end
        end
    end
end

for _target_choice in ["A", "S1", "S2"]
    main(target_choice = _target_choice)
end