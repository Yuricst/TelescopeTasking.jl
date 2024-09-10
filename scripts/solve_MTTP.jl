"""Script to solve Multi Telescope Tasking Problem (MTTP)"""

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

function main()
    # load Earth parameters
    # eop_iau1980 = fetch_iers_eop()
    eop_file = joinpath(@__DIR__, "..", "data", "eop_iau1980", "finals.all.csv")
    eop_iau1980 = read_iers_eop(eop_file, Val(:IAU1980))

    config_filenames = [
        #"config_MTTP1.json",
        #"config_MTTP2.json",
        #"config_MTTP3.json",
        "config_MTTP4.json",
    ]

    for config_filename in config_filenames

        # load config jsons
        config_telescope = JSON.parsefile(joinpath(@__DIR__, "configs/config_telescope.json"))
        config = JSON.parsefile(joinpath(@__DIR__, "configs", config_filename))
        target_choice = "A"
        solver_choice = "Gurobi"

        # choose solver
        if solver_choice == "Gurobi"
            solver = MOI.OptimizerWithAttributes(
                Gurobi.Optimizer,
                "TimeLimit" => 600,
                "Method" => 0,
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

        observer_lla_per_telescope = Vector[]
        jd0_obs_per_telescope = Real[]
        obs_duration_per_telescope = Real[]
        jd0_ref_per_telescope = Real[]
        for observer in config["observers"]
            # get observer location
            observer_lat = deg2rad(observer["latitude"])          # degrees --> radians
            observer_lon = deg2rad(observer["longitude"])         # degrees --> radians
            observer_alt = observer["altitude"]                   # meters
            observer_lla = [observer_lat, observer_lon, observer_alt]
            push!(observer_lla_per_telescope, observer_lla)

            # initial epoch of local nightfall
            jd0_ref = observer["jd0_ref"]
            @assert maximum([tle_epoch(tle) for tle in tles]) <= jd0_ref "TLEs are later than reference JD!"
            jds_night = TelescopeTasking.earliest_night(jd0_ref, observer_lla, eop_iau1980)
            jd0_obs = jds_night[1]
            obs_duration = 86400 * (jds_night[2] - jds_night[1])
            push!(jd0_obs_per_telescope, jd0_obs)
            push!(obs_duration_per_telescope, obs_duration)
            push!(jd0_ref_per_telescope, jd0_ref)

            @printf("Night for observer in %s starts at MJD %1.3f and lasts %1.2f hours\n", 
                observer["city"], jd0_obs - 2400000.5, obs_duration/3600)
        end

        # create passes (irrespective of number of exposures)
        passes_per_telescope = Vector{Vector{TelescopeTasking.VisiblePass}}()
        for (q, (observer_lla, jd0_obs, obs_duration)) in enumerate(zip(observer_lla_per_telescope,
                                                                    jd0_obs_per_telescope,
                                                                    obs_duration_per_telescope))
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
            )
            push!(passes_per_telescope, _passes)
            @printf("Detected %d passes for telescope %d\n", length(_passes), q)

            if length(_passes) == 0
                @printf("No passes detected, skipping experiment with E = %d\n", num_exposure)
                continue
            end
        end

        # iterate through num_exposure
        num_exposures = [1,]# 2, 3]
        for num_exposure in num_exposures
            _experiment_name = config["name"] * "_target$(target_choice)_E$(num_exposure)"
            println(" *************** Experiment name: $_experiment_name *************** ")
            times_measure = [time(),]

            # construct problem
            _problem = TelescopeTasking.MultiTelescopeTaskingProblem(
                passes_per_telescope, 
                num_exposure,
                slew_rate;
                buffer_times = buffer_times
            )
            push!(times_measure, time())
            @printf("Constructed problem - Detected %d observable targets with %d exposure; interval time: %1.4f sec\n",
                _problem.m, num_exposure, times_measure[end] - times_measure[end-1])
            if _problem.m == 0
                @printf("No passes detected, skipping experiment with E = %d\n", num_exposure)
                continue
            end

            # solve problem
            _X, _Y, _Y_per_telescope, _solve_stats_dict = TelescopeTasking.solve(_problem, solver; verbose = true)
            push!(times_measure, time())
            @printf("Finished solve; interval time: %1.4f sec\n", times_measure[end] - times_measure[end-1])

            # save to dictionary
            _solution_dict = TelescopeTasking.MTTP_solution_to_dict(
                _problem,
                passes_per_telescope,
                jd0_ref_per_telescope,
                obs_duration_per_telescope,
                min_elevation,
                min_obs_duration,
                exposure_duration,
                observer_lla_per_telescope,
                _X,
                _Y_per_telescope,
            )
            open(joinpath(@__DIR__, "solutions", "solution_$(_experiment_name)_$(solver_choice).json"), "w") do io
                write(io, JSON.json(_solution_dict))
            end

            open(joinpath(@__DIR__, "solutions", "solve_stats_$(_experiment_name)_$(solver_choice).json"), "w") do io
                write(io, JSON.json(_solve_stats_dict))
            end
        end
    end
end

main()