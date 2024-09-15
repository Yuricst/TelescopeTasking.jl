"""Tabulate RSO passes for MTTP"""

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
config = JSON.parsefile(joinpath(@__DIR__, "configs/config_MTTP4.json"))
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

    observer_lla_per_telescope = Vector[]
    jd0_obs_per_telescope = Real[]
    obs_duration_per_telescope = Real[]
    for observer in config["observers"]
        # get observer location
        observer_lat = deg2rad(observer["latitude"])          # degrees --> radians
        observer_lon = deg2rad(observer["longitude"])         # degrees --> radians
        observer_alt = observer["altitude"]          # meters
        jd0_ref = observer["jd0_ref"]
        observer_lla = [observer_lat, observer_lon, observer_alt]
        push!(observer_lla_per_telescope, observer_lla)

        # initial epoch of local nightfall
        @assert maximum([tle_epoch(tle) for tle in tles]) <= jd0_ref "TLEs are later than reference JD!"
        jds_night = TelescopeTasking.earliest_night(jd0_ref, observer_lla, eop_iau1980)
        jd0_obs = jds_night[1]
        obs_duration = 86400 * (jds_night[2] - jds_night[1])
        push!(jd0_obs_per_telescope, jd0_obs)
        push!(obs_duration_per_telescope, obs_duration)
    end

    # create passes (irrespective of number of exposures)
    designators = String[]
    n_per_telescope = Int[]
    n_constraints = 0
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
        
        _designators = unique([pass.tle.international_designator for pass in _passes])
        #m = length(_designators)
        n = length(_passes)
        push!(n_per_telescope, n)
        append!(designators, _designators)
        n_constraints += n*(n-1)/2
    end
    n_total = sum(n_per_telescope)
    m_total = length(unique(designators))
    n_constraints += m_total
    @printf("  %s & \$%d\$ & \$%d\$ & \$%d\$ & \$%d\$ \\\\\n",
        target_set_name, m_total, n_total, m_total+n_total, n_constraints)
end

println("Done!")