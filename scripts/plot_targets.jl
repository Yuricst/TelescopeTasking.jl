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

# # load config jsons
# target_choice = "A"

# # load Earth parameters
# # eop_iau1980 = fetch_iers_eop()
# eop_file = joinpath(@__DIR__, "..", "data", "eop_iau1980", "finals.all.csv")
# eop_iau1980 = read_iers_eop(eop_file, Val(:IAU1980))

# # load TLE files
# tles = read_tles(read(joinpath(@__DIR__, "..", "data", "tles", "AAS25target$(target_choice).txt"), String))
# println("There are $(length(tles)) TLEs in the file")

# # plot in 3D
# fig = Figure(size=(500, 500))
# ax3 = Axis3(fig[1,1]; aspect = (1, 1, 1), xlabel="X, km", ylabel="Y, km", zlabel="Z, km")

# # propagate TLE
# dt_sec = 20
# for tle in tles#[1:10]
#     period_sec = TelescopeTasking.tle2period(tle)
#     _dts_min = (0:dt_sec:period_sec) / 60        # convert seconds to minutes
#     rvs = TelescopeTasking.integrate_sgp4(
#         tle,
#         _dts_min,
#         :eci,
#         eop_iau1980,
#         [0,0,0],
#     )
#     lines!(ax3, rvs[1,:]./1e3, rvs[2,:]./1e3, rvs[3,:]./1e3, color=:red, linewidth=0.2)
# end

files = ["cosmos-1408-debris.txt", "cosmos-2251-debris.txt", "fengyun-1c-debris.txt", "iridium-33-debris.txt"]
tles_cosmos1408 = read_tles(read(joinpath(@__DIR__, "..", "data", "tles", files[1]), String))
tles_cosmos2251 = read_tles(read(joinpath(@__DIR__, "..", "data", "tles", files[2]), String))
tles_fengyun1c = read_tles(read(joinpath(@__DIR__, "..", "data", "tles", files[3]), String))
tles_iridium33 = read_tles(read(joinpath(@__DIR__, "..", "data", "tles", files[4]), String))

names = ["Kosmos 1408", "Kosmos 2251", "Fengyun 1C", "Iridium 33"]
tles_list = [tles_cosmos1408, tles_cosmos2251, tles_fengyun1c, tles_iridium33]

fontsize = 20

fig_en = Figure(size=(500,400))
ax_en = Axis(fig_en[1,1];
    yscale = log10,
    xlabelsize=fontsize,
    ylabelsize=fontsize,
    xticklabelsize=fontsize-1,
    yticklabelsize=fontsize-1,
    xlabel="Mean motion, rev/day", ylabel="Eccentricity")

fig_iW = Figure(size=(500,400))
ax_iW = Axis(fig_iW[1,1];
    xticks = [0, 90, 180, 270, 360],
    xlabelsize=fontsize,
    ylabelsize=fontsize,
    xticklabelsize=fontsize-1,
    yticklabelsize=fontsize-1,
    xlabel="RAAN, deg", ylabel="Inclination, deg")
colors = cgrad(:matter, 5, categorical = true)[2:end]
markers = [:circle, :diamond, :cross, :rect]
for (idx,(tles, name)) in enumerate(zip(tles_list, names))
    mean_motion = [tle.mean_motion for tle in tles]
    eccentricity = [tle.eccentricity for tle in tles]
    inclination = [tle.inclination for tle in tles]
    raan = [tle.raan for tle in tles]

    name_ = "$name ($(length(tles)))"
    scatter!(ax_en, mean_motion, eccentricity, label=name_, markersize=5,
        color=colors[idx], marker=markers[idx])
    scatter!(ax_iW, raan, inclination, label=name_, markersize=5, 
        color=colors[idx], marker=markers[idx])
end
axislegend(ax_en, position=:lb, labelsize=fontsize-5)

save(joinpath(@__DIR__, "plots", "targets_orbits_en.png"), fig_en)
save(joinpath(@__DIR__, "plots", "targets_orbits_iW.png"), fig_iW)
display(fig_en)
