"""
Handling TLEs
"""

using Colors
using ColorSchemes
using GeometryBasics
using GLMakie
using GLPK
using Gurobi
using JuMP
using LinearAlgebra
using ProgressMeter: @showprogress
using Printf: @printf
using SatelliteToolboxTle
using SatelliteToolboxSgp4
using SatelliteToolboxTransformations

include(joinpath(@__DIR__, "../src/TelescopeScheduling.jl"))


function dcm_ned_to_enu()
    T3 = [0 1 0;
          -1 0 0;
          0 0 1];
    T1 = [1 0 0;
          0 -1 0;
          0 0 -1];
    return T1 * T3
end

function cart2sph(rvec)
    az = atan(rvec[2], rvec[1])
    el = atan(rvec[3], sqrt(rvec[1]^2 + rvec[2]^2))
    r = norm(rvec)
    return [az; el; r]
end

function periodic_lines!(ax, xs, ys, diff_max; color=:red, linewidth=0.5)
    i_from = 1
    for i in 1:length(xs)-1
        diff = abs(xs[i] - xs[i+1])
        if diff > diff_max
            lines!(ax, xs[i_from:i], ys[i_from:i], color=color, linewidth=linewidth)
            i_from = i + 1
        end
    end
end

T_NED2ENU = dcm_ned_to_enu()

# observer location in geodetic coordinates
lat = 45.5
lon = 135.3
alt = 100.0         # in meters
r_ECEF_obs = geodetic_to_ecef(lat |> deg2rad, lon |> deg2rad, alt) / 1e3    # in km
r_ENU_obs = T_NED2ENU * ecef_to_ned(r_ECEF_obs, lat |> deg2rad, lon |> deg2rad, alt; translate = false)

# initial epoch of local nightfall
jd0_obs = 2.46055755221534e6 + 0.5      # in julian date
obs_duration_min = 8 * 60               # in minutes

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
names_include = ["GLOBALSTAR", "IRIDIUM"]
TelescopeScheduling.filter!(tles, names_include=names_include)
@show length(tles)

# get passes
min_elevation = deg2rad(15)
min_obs_duration = 10 * 600         # in seconds
exposure_duration = 60              # in seconds
observer_lat = deg2rad(20)
observer_lon = deg2rad(100)
observer_alt = 30.0
observer_lla = [observer_lat, observer_lon, observer_alt]

passes = TelescopeScheduling.tles_to_passes(
    tles[1:100],
    eop_iau1980,
    jd0_obs,
    min_elevation,
    min_obs_duration,
    exposure_duration,
    observer_lla,
    dt_sec=10,
)
println("Done!")