"""Plot local night time based on Sun location"""


using GLMakie
using GeometryBasics
using LinearAlgebra
using ProgressMeter: @showprogress
using SatelliteToolboxCelestialBodies
using SatelliteToolboxTle
using SatelliteToolboxSgp4
using SatelliteToolboxTransformations

include(joinpath(@__DIR__, "../src/TelescopeTasking.jl"))

T_NED2ENU = TelescopeTasking.dcm_ned_to_enu()

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
path_to_tles = joinpath(@__DIR__, "..", "data", "tles", "planet.txt")
tles_str = read(path_to_tles, String)

# convert to TLE objects
tles = read_tles(tles_str)
jd0 = tle_epoch(tles[1])

jds = LinRange(jd0, jd0+2.0, 5000)

# sun position
observer_lla = [deg2rad(45), 0, 100]
isun = zeros(3, length(jds))
sph_sun = zeros(3, length(jds))
for (i,jd) in enumerate(jds)
    isun[:,i] = TelescopeTasking.sun_direction_enu(jd, observer_lla, eop_iau1980)
    sph_sun[:,i] = TelescopeTasking.cart2sph(isun[:,i])
end

# get night information
jds_night = TelescopeTasking.earliest_night(jd0, observer_lla, eop_iau1980)

# plot sun information
figure = Figure(size=(1000,400))
ax3d = Axis3(figure[1,1]; aspect = (1, 1, 1))
lines!(ax3d, isun[1,:], isun[2,:], isun[3,:], color=:red)

axAzEl = Axis(figure[1,2]; xlabel="Azimuth, deg", ylabel="Elevation, deg")
xlims!(axAzEl, [-180,180])
ylims!(axAzEl, [-90,90])
scatter!(axAzEl, rad2deg.(sph_sun[1,1:1]), rad2deg.(sph_sun[2,1:1]), 
    color = :red,)
# TelescopeTasking.periodic_lines!(axAzEl, rad2deg.(sph_sun[1,:]), rad2deg.(sph_sun[2,:]), 160,
#     color = :red,
#     linewidth = 1.0)
lines!(axAzEl, rad2deg.(sph_sun[1,:]), rad2deg.(sph_sun[2,:]), 
    color = :red,
    linewidth = 1.0)

ax_time = Axis(figure[1,3]; xlabel="Elapsed time, JD", ylabel="Elevation, deg")
lines!(ax_time, jds .- jd0, rad2deg.(sph_sun[2,:]), color=:red)
vlines!(ax_time, jds_night .- jd0, color=:black, linestyle=:dash)
display(figure)