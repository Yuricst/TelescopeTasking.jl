"""Test with passes IO functions"""

using Test
using SatelliteToolboxTransformations

include(joinpath(@__DIR__, "../src/TelescopeTasking.jl"))

@testset "Test for IO with passes" begin
    # load Earth parameters
    # eop_iau1980 = fetch_iers_eop()
    eop_file = joinpath(@__DIR__, "..", "data", "eop_iau1980", "finals.all.csv")
    eop_iau1980 = read_iers_eop(eop_file, Val(:IAU1980))

    # load tles
    tles = read_tles(read(joinpath(@__DIR__, "../data/tles/iridium-33-debris.txt"), String))
    
    # get passes
    jd0_obs = 2.46055755221534e6 + 0.45     # in julian date
    obs_duration = 8 * 3600             # in seconds
    min_elevation = deg2rad(15)
    min_obs_duration = 100              # in seconds
    exposure_duration = 60              # in seconds
    observer_lat = deg2rad(45)
    observer_lon = deg2rad(100)
    observer_alt = 30.0
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

    # convert them to dictionaries
    passes_dict = [TelescopeTasking.pass_to_dict(pass) for pass in passes]

    # reconstruct passes
    passes_reconstruct = [TelescopeTasking.VisiblePass(pass_dict) for pass_dict in passes_dict]

    # tests
    for (p1,p2) in zip(passes, passes_reconstruct)
        @test p1.tle == p2.tle
        @test p1.times == p2.times
        @test p1.azimuths == p2.azimuths
        @test p1.elevations == p2.elevations
    end
end