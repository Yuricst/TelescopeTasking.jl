"""Tabulate telescopes"""

using CairoMakie
using JSON
using Printf: @printf
using SatelliteToolboxTle
using SatelliteToolboxTransformations

include(joinpath(@__DIR__, "../src/TelescopeTasking.jl"))

# load Earth parameters
# eop_iau1980 = fetch_iers_eop()
eop_file = joinpath(@__DIR__, "..", "data", "eop_iau1980", "finals.all.csv")
eop_iau1980 = read_iers_eop(eop_file, Val(:IAU1980))

# load general telescope configuration
config_telescope = JSON.parsefile(joinpath(@__DIR__, "configs/config_telescope.json"))

# get telescopes A, C
config = JSON.parsefile(joinpath(@__DIR__, "configs/config_MTTP1.json"))
telescope_A = config["observers"][1]
telescope_C = config["observers"][2]

# get telescope B
config = JSON.parsefile(joinpath(@__DIR__, "configs/config_STTP2.json"))
telescope_B = config["observer"]

# get telescope C
config = JSON.parsefile(joinpath(@__DIR__, "configs/config_MTTP2.json"))
telescope_D = config["observers"][2]

# get telescope D, E
config = JSON.parsefile(joinpath(@__DIR__, "configs/config_MTTP3.json"))
telescope_E = config["observers"][2]
telescope_F = config["observers"][3]

telescopes = [
    telescope_A,
    telescope_B,
    telescope_C,
    telescope_D,
    telescope_E,
    telescope_F,
]
names = ["A", "B", "C", "D", "E", "F"]
for (name,telescope) in zip(names, telescopes)
    # compute night epochs 
    _observer_lla = [
        deg2rad(telescope["latitude"]),
        deg2rad(telescope["longitude"]),
        telescope["altitude"],
    ]
    _jds_night = TelescopeTasking.earliest_night(config_telescope["jd0_ref"], _observer_lla, eop_iau1980)
    _obs_duration = 24 * (_jds_night[2] - _jds_night[1])

    # print for table
    @printf(
        "     %s  & \$%1.3f\$ & \$%1.3f\$ & \$%1.0f\$ & \$%1.3f\$ & \$%1.3f\$ & \$%1.2f\$ & & \\\\\n",
        name,
        telescope["latitude"],
        telescope["longitude"],
        telescope["altitude"],
        _jds_night[1] - 2400000.5, 
        _jds_night[2] - 2400000.5,
        _obs_duration
    )
end

# load TLE files
target_choice = "B"
tles = read_tles(read(joinpath(@__DIR__, "..", "data", "tles", "AAS25target$(target_choice).txt"), String))
println("There are $(length(tles)) TLEs in the file")

# compute passes
config_telescope = JSON.parsefile(joinpath(@__DIR__, "configs/config_telescope.json"))
slew_rate           = deg2rad(config_telescope["slew_rate"])                         # rad/s
buffer_times        = [config_telescope["buffer_t0"], config_telescope["buffer_t1"]]     # times in seconds
min_elevation       = deg2rad(config_telescope["min_elevation"] )            # in radians
min_obs_duration    = config_telescope["min_obs_duration"]                # in seconds
exposure_duration   = config_telescope["exposure_duration"]              # in seconds

for (name,telescope) in zip(names, telescopes)
    # compute night epochs 
    _observer_lla = [
        deg2rad(telescope["latitude"]),
        deg2rad(telescope["longitude"]),
        telescope["altitude"],
    ]
    _jds_night = TelescopeTasking.earliest_night(config_telescope["jd0_ref"], _observer_lla, eop_iau1980)
    _obs_duration = 86400 * (_jds_night[2] - _jds_night[1])
    
    # get passes
    _passes, _ = TelescopeTasking.tles_to_passes(
        tles,
        eop_iau1980,
        _jds_night[1],
        _obs_duration,
        min_elevation,
        min_obs_duration,
        exposure_duration,
        _observer_lla,
        dt_sec = 10,
    )
    _designators = unique([pass.tle.international_designator for pass in _passes])
    m = length(_designators)
    n = length(_passes)
    @printf("telescope %s : RSO = %d, passes = %d\n", name, m, n)
end


# sun position
jd0 = config_telescope["jd0_ref"]
jds = LinRange(jd0, jd0+2.0, 5000)
observer_lla = [
    deg2rad(telescope_E["latitude"]),
    deg2rad(telescope_E["longitude"]),
    telescope_E["altitude"],
]
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