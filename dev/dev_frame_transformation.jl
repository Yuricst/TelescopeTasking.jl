"""
Frame transformation development
"""

using GLMakie
using GeometryBasics
using LinearAlgebra
using ProgressMeter: @showprogress
using SatelliteToolboxTle
using SatelliteToolboxSgp4
using SatelliteToolboxTransformations


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

# load Earth parameters
# eop_iau1980 = fetch_iers_eop()
eop_file = joinpath(@__DIR__, "..", "data", "eop_iau1980", "finals.all.csv")
eop_iau1980 = read_iers_eop(eop_file, Val(:IAU1980))

# load TLE files
path_to_tles = joinpath(@__DIR__, "..", "data", "tles", "active_GPS.txt")
tles_str = read(path_to_tles, String)

# convert to TLE objects
tles = read_tles(tles_str)

# propagate one of them
jd0 = tle_epoch(tles[1])            # initial epoch of TLE in Julian day

sgp4d = sgp4_init(tles[1])
dt_min = 0.5
dts_min = 0:dt_min:30*24*60        # propagate over 6 hours
rs_TEME, vs_TEME = zeros(3, length(dts_min)), zeros(3, length(dts_min))
rs_ITRF = zeros(3, length(dts_min))
rs_ENU = zeros(3, length(dts_min))
sph_ENU = zeros(3, length(dts_min))

@showprogress for (idx,dt_min) in enumerate(dts_min)
    # propagate state in TEME frame
    rs_TEME[:,idx], vs_TEME[:,idx] = sgp4!(sgp4d, dt_min)

    # current julian date
    jd = jd0 + dt_min/60/24

    # get transformation matrix from ECI to ECEF and apply conversion
    T_TEME2ITRF = r_eci_to_ecef(TEME(), ITRF(), jd, eop_iau1980)
    rs_ITRF[:,idx] = T_TEME2ITRF * rs_TEME[:,idx]

    # convert to azimuth and elevation angles
    r_NED = ecef_to_ned(rs_ITRF[:,idx], lat |> deg2rad, lon |> deg2rad, alt; translate = false)
    rs_ENU[:,idx] = T_NED2ENU * r_NED - r_ENU_obs
    sph_ENU[:,idx] = cart2sph(rs_ENU[:,idx])    # [az,el,radius]
end

# plot in 3D
fig = Figure(size=(1200, 400))
ax1 = Axis3(fig[1,1]; aspect = (1, 1, 1), title="TEME frame")
ax2 = Axis3(fig[1,2]; aspect = (1, 1, 1), title="ITRF frame")
ax3 = Axis3(fig[1,3]; aspect = (1, 1, 1), title="ENU frame")
ax4 = Axis(fig[1,4]; xlabel="Azimuth, deg", ylabel="Elevation, deg")

lines!(ax1, rs_TEME[1,:], rs_TEME[2,:], rs_TEME[3,:], color=:red, linewidth=0.8)
lines!(ax2, rs_ITRF[1,:], rs_ITRF[2,:], rs_ITRF[3,:], color=:red, linewidth=0.8)
scatter!(ax2, [r_ECEF_obs[1]], [r_ECEF_obs[2]], [r_ECEF_obs[3]], color=:blue, marker=:circle)
lines!(ax3, rs_ENU[1,:], rs_ENU[2,:], rs_ENU[3,:], color=:red, linewidth=0.8)
scatter!(ax3, [0], [0], [0], color=:blue, marker=:circle)

mesh!(ax1, Sphere(GeometryBasics.Point3f0(0), 6378.0), alpha=0.5)
mesh!(ax2, Sphere(GeometryBasics.Point3f0(0), 6378.0), alpha=0.5)

periodic_lines!(ax4, rad2deg.(sph_ENU[1,:]), rad2deg.(sph_ENU[2,:]), 30, color=:red, linewidth=0.8)

display(fig)