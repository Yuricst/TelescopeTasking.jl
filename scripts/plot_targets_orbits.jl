"""Plot targets' orbits"""

using Colors
using ColorSchemes
using GeometryBasics
using CairoMakie
using GLPK
using Gurobi
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

# load TLE files
tles = read_tles(read(joinpath(@__DIR__, "..", "data", "tles", "AAS25targetS0.txt"), String))
println("There are $(length(tles)) TLEs in the file")


# load Earth parameters
# eop_iau1980 = fetch_iers_eop()
eop_file = joinpath(@__DIR__, "..", "data", "eop_iau1980", "finals.all.csv")
eop_iau1980 = read_iers_eop(eop_file, Val(:IAU1980))

# plot in 3D
colors = cgrad(:matter, 5, categorical = true)[2:end]
fontsize = 32
fig_i1 = Figure(size=(500, 500))
ax3_1 = Axis3(fig_i1[1,1]; aspect = :data, xlabel="X, km", ylabel="Y, km", zlabel="Z, km",
    xlabelsize=fontsize,
    ylabelsize=fontsize,
    zlabelsize=fontsize,
    xticklabelsize=fontsize-1,
    yticklabelsize=fontsize-1,
    zticklabelsize=fontsize-1)

fig_i2 = Figure(size=(500, 500))
ax3_2 = Axis3(fig_i2[1,1]; aspect = :data, xlabel="X, km", ylabel="Y, km", zlabel="Z, km",
    xlabelsize=fontsize,
    ylabelsize=fontsize,
    zlabelsize=fontsize,
    xticklabelsize=fontsize-1,
    yticklabelsize=fontsize-1,
    zticklabelsize=fontsize-1)

fig_i3 = Figure(size=(500, 500))
ax3_3 = Axis3(fig_i3[1,1]; aspect = :data, xlabel="X, km", ylabel="Y, km", zlabel="Z, km",
    xlabelsize=fontsize,
    ylabelsize=fontsize,
    zlabelsize=fontsize,
    xticklabelsize=fontsize-1,
    yticklabelsize=fontsize-1,
    zticklabelsize=fontsize-1)

fig_i4 = Figure(size=(500, 500))
ax3_4 = Axis3(fig_i4[1,1]; aspect = :data, xlabel="X, km", ylabel="Y, km", zlabel="Z, km",
    xlabelsize=fontsize,
    ylabelsize=fontsize,
    zlabelsize=fontsize,
    xticklabelsize=fontsize-1,
    yticklabelsize=fontsize-1,
    zticklabelsize=fontsize-1)

# propagate TLE
n_per_inclinations = [0,0,0,0]
dt_sec = 20
@showprogress for tle in tles
    if tle.inclination < 50
        _ax = ax3_1
        _color = colors[1]
        n_per_inclinations[1] += 1
    elseif tle.inclination < 60
        _ax = ax3_2
        _color = colors[2]
        n_per_inclinations[2] += 1
    elseif tle.inclination < 80
        _ax = ax3_3
        _color = colors[3]
        n_per_inclinations[3] += 1
    else
        _ax = ax3_4
        _color = colors[4]
        n_per_inclinations[4] += 1
    end
    period_sec = TelescopeTasking.tle2period(tle)
    _dts_min = (0:dt_sec:period_sec) / 60        # convert seconds to minutes
    rvs = TelescopeTasking.integrate_sgp4(
        tle,
        _dts_min,
        :eci,
        eop_iau1980,
        [0,0,0],
    )
    lines!(_ax, rvs[1,:]./1e3, rvs[2,:]./1e3, rvs[3,:]./1e3, color=_color, linewidth=0.2)
end
@show n_per_inclinations

# # plot Earth
# for ax3 in [ax3_1, ax3_2, ax3_3, ax3_4]
#     mesh!(ax3, Sphere(GeometryBasics.Point3f0(0), 6378.0/1e3), alpha=0.7)
# end

save(joinpath(@__DIR__, "plots", "orbits_Starlink_ECI_i1.png"), fig_i1)
save(joinpath(@__DIR__, "plots", "orbits_Starlink_ECI_i2.png"), fig_i2)
save(joinpath(@__DIR__, "plots", "orbits_Starlink_ECI_i3.png"), fig_i3)
save(joinpath(@__DIR__, "plots", "orbits_Starlink_ECI_i4.png"), fig_i4)
display(fig_i1)