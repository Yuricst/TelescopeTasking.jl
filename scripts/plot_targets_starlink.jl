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
tles_S1 = read_tles(read(joinpath(@__DIR__, "..", "data", "tles", "AAS25targetS1.txt"), String))
tles_S2 = read_tles(read(joinpath(@__DIR__, "..", "data", "tles", "AAS25targetS2.txt"), String))
tles_S3 = read_tles(read(joinpath(@__DIR__, "..", "data", "tles", "AAS25targetS3.txt"), String))
tles_S4 = read_tles(read(joinpath(@__DIR__, "..", "data", "tles", "AAS25targetS4.txt"), String))
println("There are $(length(tles)) TLEs in the file")

fontsize = 26

fig_en = Figure(size=(500,400))
ax_en = Axis(fig_en[1,1];
    yscale = log10,
    xlabelsize=fontsize,
    ylabelsize=fontsize,
    xticklabelsize=fontsize-1,
    yticklabelsize=fontsize-1,
    xlabel="Semimajor axis, km", ylabel="Eccentricity")

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

TG = 86164.0

colors = cgrad(:matter, 5, categorical = true)[2:end]
names = ["Shell 1", "Shell 2", "Shell 3", "Shell 4"]
for (idx,(tles,set_name)) in enumerate(zip([tles_S1,tles_S2,tles_S3,tles_S4],names))
    smas = [TelescopeTasking.tle2sma(tle) for tle in tles]
    mean_motion = [tle.mean_motion for tle in tles]
    periods = [TelescopeTasking.tle2period(tle) for tle in tles]
    eccentricity = [tle.eccentricity for tle in tles]
    inclination = [tle.inclination for tle in tles]
    raan = [tle.raan for tle in tles]
    DeltaE = rad2deg.(2pi  * periods ./ TG)

    set_name_ = set_name * " ($(length(tles)))"
    scatter!(ax_en, smas, eccentricity, label=set_name_, markersize=5, color=colors[idx])
    scatter!(ax_iW, raan, inclination, label=set_name_, markersize=5, color=colors[idx])
end

axislegend(ax_en, position=:lb, labelsize=fontsize-5)

save(joinpath(@__DIR__, "plots", "targets_orbits_en_STARLINK.png"), fig_en)
save(joinpath(@__DIR__, "plots", "targets_orbits_iW_STARLINK.png"), fig_iW)
display(fig_en)
