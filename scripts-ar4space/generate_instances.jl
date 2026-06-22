"""Systematically generate instances for extensive testing"""

using Colors
using ColorSchemes
using GeometryBasics
using GLMakie
using Gurobi
using HiGHS
using JSON
using JuMP
using LinearAlgebra
using ProgressMeter: @showprogress
using Printf: @printf
using Random
using SatelliteToolboxTle
using SatelliteToolboxSgp4
using SatelliteToolboxTransformations


if !@isdefined(passes_per_telescope__)
    include(joinpath(@__DIR__, "../src/TelescopeTasking.jl"))

    # load Earth parameters
    # eop_iau1980 = fetch_iers_eop()
    eop_file = joinpath(@__DIR__, "..", "data", "eop_iau1980", "finals.all.csv")
    eop_iau1980 = read_iers_eop(eop_file, Val(:IAU1980))

    # choose instance 
    target_choice = "debris"     # A or S1 or S2 or debris
    config_filename = "config_MTTP4.json"
    num_exposure = 1       # 1, 2, or 3

    # load config jsons
    config_telescope = JSON.parsefile(joinpath(@__DIR__, "../scripts/configs/config_telescope.json"))
    config = JSON.parsefile(joinpath(@__DIR__, "../scripts/configs", config_filename))
    solver_choice = "Gurobi"    # Gurobi or GPLK or HiGHS
    time_limit = 3600           # 600 or 3600, in seconds
    save_dir = "solutions_JASS"

    # define names for logging
    _experiment_name = config["name"] * "_target$(target_choice)_E$(num_exposure)"
    filepath_log = "log_$(_experiment_name)_$(solver_choice).log"
    filepath_solution = "solution_$(_experiment_name)_$(solver_choice).json"
    filepath_stats = "solve_stats_$(_experiment_name)_$(solver_choice).json"

    # choose solver
    if solver_choice == "Gurobi"
        solver = MOI.OptimizerWithAttributes(
            Gurobi.Optimizer,
            "TimeLimit" => time_limit,
            # "LogFile" => filepath_log,
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
    if target_choice == "debris"
        tles = read_tles(read(joinpath(@__DIR__, "..", "data", "tles", "$(target_choice).txt"), String))
    else
        tles = read_tles(read(joinpath(@__DIR__, "..", "data", "tles", "AAS25target$(target_choice).txt"), String))
    end
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

    if occursin("STTP", config["name"])
        # get observer location
        observer = config["observer"]
        observer_lat = deg2rad(observer["latitude"])          # degrees --> radians
        observer_lon = deg2rad(observer["longitude"])         # degrees --> radians
        observer_alt = observer["altitude"]                   # meters
        observer_lla = [observer_lat, observer_lon, observer_alt]
        push!(observer_lla_per_telescope, observer_lla)

        # initial epoch of local nightfall
        jd0_ref = config_telescope["jd0_ref"]
        @assert maximum([tle_epoch(tle) for tle in tles]) <= jd0_ref "TLEs are later than reference JD!"
        jds_night = TelescopeTasking.earliest_night(jd0_ref, observer_lla, eop_iau1980)
        jd0_obs = jds_night[1]
        obs_duration = 86400 * (jds_night[2] - jds_night[1])
        push!(jd0_obs_per_telescope, jd0_obs)
        push!(obs_duration_per_telescope, obs_duration)
        push!(jd0_ref_per_telescope, jd0_ref)

        @printf("Night for observer in %s starts at MJD %1.3f and lasts %1.2f hours\n", 
            observer["city"], jd0_obs - 2400000.5, obs_duration/3600)
    else
        for observer in config["observers"]
            # get observer location
            _observer_lat = deg2rad(observer["latitude"])          # degrees --> radians
            _observer_lon = deg2rad(observer["longitude"])         # degrees --> radians
            _observer_alt = observer["altitude"]                   # meters
            _observer_lla = [_observer_lat, _observer_lon, _observer_alt]
            push!(observer_lla_per_telescope, _observer_lla)

            # initial epoch of local nightfall
            _jd0_ref = observer["jd0_ref"]
            @assert maximum([tle_epoch(tle) for tle in tles]) <= _jd0_ref "TLEs are later than reference JD!"
            _jds_night = TelescopeTasking.earliest_night(_jd0_ref, _observer_lla, eop_iau1980)
            _jd0_obs = _jds_night[1]
            _obs_duration = 86400 * (_jds_night[2] - _jds_night[1])
            push!(jd0_obs_per_telescope, _jd0_obs)
            push!(obs_duration_per_telescope, _obs_duration)
            push!(jd0_ref_per_telescope, _jd0_ref)

            @printf("Night for observer in %s starts at MJD %1.3f and lasts %1.2f hours\n", 
                observer["city"], _jd0_obs - 2400000.5, _obs_duration/3600)
        end
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
end

function dropout_passes(passes, ID)
    isempty(passes) && return passes
    Random.seed!(ID)
    N_sample = min(length(passes), Int(ceil(length(passes) * ID / 100)))
    indices = randperm(length(passes))[1:N_sample]
    return TelescopeTasking.sort(passes[indices])
end

save_directory = joinpath(@__DIR__, "problem-data-test", target_choice)
@showprogress for ID in 1:100
    # keep ceil(ID/100) fraction of passes per telescope (unique, time-sorted)
    if ID < 100
        _passes_per_telescope_dropped = [dropout_passes(passes, ID) for passes in passes_per_telescope]
    else
        _passes_per_telescope_dropped = passes_per_telescope
    end

    # construct problem & solve

    # construct problem
    _problem_ID = TelescopeTasking.MultiTelescopeTaskingProblem(
        _passes_per_telescope_dropped, 
        num_exposure,
        slew_rate;
        buffer_times = buffer_times
    )
    if _problem_ID.m == 0
        @printf("No passes detected, skipping instance ID = %d\n", ID)
        continue
    end

    # # get JuMP model & print to file
    # model = TelescopeTasking.solve(_problem, solver; get_model = true)
    # save_instance_file = joinpath(@__DIR__, "model_$(config["name"])" * "_target$(target_choice)_E$(num_exposure).mps")
    # write_to_file(model, save_instance_file)

    # export to JSON file
    _instance_file = Dict(
        "ID" => ID,
        "config_filename" => config_filename,
        "target_choice" => target_choice,
        "m" => _problem_ID.m,
        "n_per_telescope" => _problem_ID.n_per_telescope,
        "n_total" => _problem_ID.n_total,
        "A_index_per_telescope" => [[(i, j) for i in axes(_A, 1), j in axes(_A, 2) if _A[i, j] > 0.5] for _A in _problem_ID.A_per_telescope],
        "T_per_telescope" => _problem_ID.T_per_telescope,
    )

    _save_instance_file = "instance_ID$(ID).json"
    open(joinpath(save_directory, _save_instance_file), "w") do io
        write(io, JSON.json(_instance_file))
    end
    _xz_file = joinpath(save_directory, _save_instance_file * ".xz")
    if isfile(_xz_file)
        rm(_xz_file)
    end
    run(`xz $(joinpath(save_directory, _save_instance_file))`)
    # println("Saved instance file to $(_save_instance_file).xz!")
end
println("Saved instances to $(save_directory)!")