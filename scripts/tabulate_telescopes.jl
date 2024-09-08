"""Tabulate telescopes"""

using JSON

include(joinpath(@__DIR__, "../src/TelescopeTasking.jl"))


config_telescope = JSON.parsefile(joinpath(@__DIR__, "configs/config_telescope.json"))

# get telescopes A, B
config = JSON.parsefile(joinpath(@__DIR__, "configs/config_MTTP1.json"))
telescope_A = config["observers"][1]
telescope_B = config["observers"][2]

# get telescope C
config = JSON.parsefile(joinpath(@__DIR__, "configs/config_MTTP2.json"))
telescope_C = config["observers"][2]

# get telescope D, E
config = JSON.parsefile(joinpath(@__DIR__, "configs/config_MTTP3.json"))
telescope_D = config["observers"][2]
telescope_E = config["observers"][3]

telescopes = [
    telescope_A,
    telescope_B,
    telescope_C,
    telescope_D,
    telescope_E,
]
names = ["A", "B", "C", "D", "E"]
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
        "     %s  & \$%1.3f\$ & \$%1.3f\$ & \$%1.0f\$ & \$%1.3f\$ & \$%1.3f\$ & \$%1.2f\$ \\\\\n",
        name,
        telescope["latitude"],
        telescope["longitude"],
        telescope["altitude"],
        _jds_night[1], 
        _jds_night[2],
        _obs_duration
    )
end

# sun position
jd0 = config_telescope["jd0_ref"]
jds = LinRange(jd0, jd0+2.0, 5000)
observer_lla = [
    deg2rad(telescope_D["latitude"]),
    deg2rad(telescope_D["longitude"]),
    telescope_D["altitude"],
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