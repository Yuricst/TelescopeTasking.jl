"""
Frame transformation development
"""

using Colors
using ColorSchemes
using GLMakie
using GeometryBasics
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

# criteria for valid pass
min_elevation = deg2rad(15)
min_duration_day = 60/86400        # minimum duration: 1 minute, in days

# load Earth parameters
# eop_iau1980 = fetch_iers_eop()
eop_file = joinpath(@__DIR__, "..", "data", "eop_iau1980", "finals.all.csv")
eop_iau1980 = read_iers_eop(eop_file, Val(:IAU1980))

# load TLE files
path_to_tles = joinpath(@__DIR__, "..", "data", "tles", "iridium-33-debris.txt")
tles_str = read(path_to_tles, String)

# convert to TLE objects
tles = read_tles(tles_str)

# initialize figure
fig = Figure(size=(800, 800))
ax1 = Axis(fig[1,1]; xlabel="Azimuth, deg", ylabel="Elevation, deg")
colors = cgrad(:hawaii, max(2,length(tles)))

visible_arcs = []
passes = TelescopeScheduling.VisiblePass[]

@showprogress for (idx,tle) in enumerate(tles)
    _jd0_tle = tle_epoch(tle)            # initial epoch of TLE in Julian day
    _sgp4d = sgp4_init(tle)
    _dt_min = 0.1
    
    t0_tle = (jd0_obs - _jd0_tle) * 24*60            # time between TLE time and observation start time, in minutes
    dts_min_tle = t0_tle:_dt_min:(t0_tle + 12*60)    # propagate over 6 hours, in minutes

    _rs_TEME, _vs_TEME = zeros(3, length(dts_min_tle)), zeros(3, length(dts_min_tle))
    _rs_ITRF = zeros(3, length(dts_min_tle))
    _rs_ENU = zeros(3, length(dts_min_tle))
    _sph_ENU = zeros(3, length(dts_min_tle))

    @showprogress for (idx,dt_min) in enumerate(dts_min_tle)
        # propagate state in TEME frame
        _rs_TEME[:,idx], _vs_TEME[:,idx] = sgp4!(_sgp4d, dt_min)

        # current julian date
        jd = _jd0_tle + dt_min/60/24

        # get transformation matrix from ECI to ECEF and apply conversion
        T_TEME2ITRF = r_eci_to_ecef(TEME(), ITRF(), jd, eop_iau1980)
        _rs_ITRF[:,idx] = T_TEME2ITRF * _rs_TEME[:,idx]

        # convert to azimuth and elevation angles
        r_NED = ecef_to_ned(_rs_ITRF[:,idx], lat |> deg2rad, lon |> deg2rad, alt; translate = false)
        _rs_ENU[:,idx] = T_NED2ENU * r_NED - r_ENU_obs
        _sph_ENU[:,idx] = cart2sph(_rs_ENU[:,idx])    # [az,el,radius]
    end

    # plot elevation against azimuth
    periodic_lines!(ax1, rad2deg.(_sph_ENU[1,:]), rad2deg.(_sph_ENU[2,:]), 30;
        color=colors[idx], linewidth=0.4)

    # get visible passes
    times_jd = Array(_jd0_tle .+ dts_min_tle/60/24)
    _passes = TelescopeScheduling.azel_history_to_passes(
        tle, times_jd, _sph_ENU[1,:], _sph_ENU[2,:], min_elevation, min_duration_day)
    push!(passes, _passes...)
end

# plot only visible passes
ax2 = Axis(fig[2,1]; xlabel="Azimuth, deg", ylabel="Elevation, deg")
for pass in passes
    periodic_lines!(ax2, rad2deg.(pass.azimuths), rad2deg.(pass.elevations), 30;
        color=:black, linewidth=0.8)
end

# formatting plots
for _ax in [ax1, ax2]
    ylims!(_ax, 0, 90)
    hlines!(_ax, [rad2deg(min_elevation)], color=:black, linestyle=:dash, linewidth=2.5)
end

display(fig)