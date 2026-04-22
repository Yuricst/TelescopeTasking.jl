"""
Handling TLEs
"""

using Colors
using ColorSchemes
using GeometryBasics
using GLMakie
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

@testset "Test for solving instance of single telescope scheduling problem" begin
    # initial epoch of local nightfall
    jd0_obs = 2.46055755221534e6 + 0.45     # in julian date

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
    tles = tles[1:100]          # only use subset of them

    # get passes
    obs_duration = 8 * 3600             # in seconds
    min_elevation = deg2rad(15)
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

    # solve problem
    solver = HiGHS.Optimizer
    X, Y, status = TelescopeTasking.solve(problem, solver; verbose=false)
    @test status["relative_gap"] <= 1e-8
end