"""Tabulate RSO passes for STTP"""

using JSON
using Printf
using SatelliteToolboxTle
using SatelliteToolboxTransformations

include(joinpath(@__DIR__, "../src/TelescopeTasking.jl"))

# load Earth parameters
eop_file = joinpath(@__DIR__, "..", "data", "eop_iau1980", "finals.all.csv")
eop_iau1980 = read_iers_eop(eop_file, Val(:IAU1980))


# load config jsons
config_telescope = JSON.parsefile(joinpath(@__DIR__, "configs/config_telescope.json"))
config = JSON.parsefile(joinpath(@__DIR__, "configs/config_STTP1.json"))
target_choices = ["A", "S1", "S2"]

for target_choice in target_choices
    if target_choice == "A"
        target_set_name = "\$G\$"
    elseif target_choice == "B"
        target_set_name = "Debris"
    elseif target_choice == "C"
        target_set_name = "Starlink"
    elseif target_choice == "S1"
        target_set_name = "\$S_1\$"
    elseif target_choice == "S2"
        target_set_name = "\$S_2\$"
    else
        error("Target choice $target_choice not recognized!")
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
    observer_alt = deg2rad(config["observer"]["altitude"])          # meters
    observer_lla = [observer_lat, observer_lon, observer_alt]

    # initial epoch of local nightfall
    jd0_ref = config_telescope["jd0_ref"]
    @assert maximum([tle_epoch(tle) for tle in tles]) <= jd0_ref "TLEs are later than reference JD!"
    jds_night = TelescopeTasking.earliest_night(jd0_ref, observer_lla, eop_iau1980)
    jd0_obs = jds_night[1]
    obs_duration = 86400 * (jds_night[2] - jds_night[1])

    # get passes
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
        num_exposure = 1,
    )
    _designators = unique([pass.tle.international_designator for pass in _passes])
    m = length(_designators)
    n = length(_passes)
    # @printf("Problem %s, target %s : RSO = %d, passes = %d\n", config["name"], target_set_name, m, n)
    @printf("  %s & \$%d\$ & \$%d\$ & \$%d\$ & \$%d\$ \\\\\n", 
        target_set_name, m, n, m+n, m+n*(n-1)/2)

    # save dictionary of passes
    passes_dict = [TelescopeTasking.pass_to_dict(pass) for pass in _passes]
    open(joinpath(@__DIR__, "passes", "passes_$(config["name"])_target$(target_choice).json"), "w") do io
        write(io, JSON.json(passes_dict))
    end
end

println("Done!")