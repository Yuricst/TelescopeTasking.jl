"""Plot MTTP solution passes"""

using CairoMakie
# using GLMakie

using Colors
using ColorSchemes
using GeometryBasics
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

# load Earth parameters
eop_file = joinpath(@__DIR__, "..", "data", "eop_iau1980", "finals.all.csv")
eop_iau1980 = read_iers_eop(eop_file, Val(:IAU1980))

# load config jsons
telescope = JSON.parsefile(joinpath(@__DIR__, "configs/config_telescope.json"))
config = JSON.parsefile(joinpath(@__DIR__, "configs/config_MTTP1.json"))
target_choice = "A"
num_exposure = 1
solver_choice = "Gurobi"
experiment_name = config["name"] * "_target$(target_choice)_E$(num_exposure)"
@printf("Plotting pass from experiment %s\n", experiment_name)

# load tles used
tles = read_tles(read(joinpath(@__DIR__, "../data/tles/AAS25target$(target_choice).txt"), String))

# load solutions
solution_dict = JSON.parsefile(joinpath(@__DIR__, "solutions/solution_$(experiment_name)_$(solver_choice).json"))
min_obs_duration = solution_dict["min_obs_duration"]
exposure_duration = solution_dict["exposure_duration"]
obs_duration_per_telescope = solution_dict["obs_duration_per_telescope"]
observer_lla_per_telescope = solution_dict["observer_lla_per_telescope"]
Y_per_telescope = solution_dict["Y_per_telescope"]
passes_per_telescope = solution_dict["passes_per_telescope"]

# passes_per_telescope = [
#     [TelescopeTasking.VisiblePass(pass_dict) for pass_dict in passes_dict_list]
#     for passes_dict_list in passes_per_telescope
# ]
# selected_passes_per_telescope = [
#     [pass for (pass, y) in zip(passes, value.(Y)) if y > 0.5]
#     for (passes,Y) in zip(passes_per_telescope, Y_per_telescope)
# ]
passes_per_telescope = [
    [TelescopeTasking.VisiblePass(pass_dict) for pass_dict in passes_dict_list]
    for passes_dict_list in passes_per_telescope
]

selected_passes_per_telescope = TelescopeTasking.filter([
    [pass for (pass, y) in zip(passes, value.(Y)) if y > 0.5]
    for (passes,Y) in zip(passes_per_telescope, Y_per_telescope)
], num_exposure)


# plot of selected passes
fontsize = 24
for q in 1:length(passes_per_telescope)
    _fig_sol = Figure(size=(500,500))
    _ax_polar = PolarAxis(_fig_sol[1,1];
        rticks = ([0,20,40,60], ["90","70","50","30"]),
        rticklabelsize=fontsize-1,
        thetaticklabelsize=fontsize-1,
    )
    TelescopeTasking.polar_plot_passes!(_ax_polar, passes_per_telescope[q]; color=:grey60, linewidth=0.1)
    TelescopeTasking.polar_plot_passes!(_ax_polar, selected_passes_per_telescope[q]; 
        linewidth=1.5, color_by_target=true, exposure_only=true)

    # saveplot
    filename = "passes_polar_$(experiment_name)_$(solver_choice)_telescope$q.png"
    save(joinpath(@__DIR__, "plots", filename), _fig_sol)
    @printf("Saved %s\n", filename)
end

# plot of selected passes
fontsize = 24
fig_time = Figure(size=(1400,800);)
for (idx, (passes, selected_passes)) in enumerate(zip(passes_per_telescope, selected_passes_per_telescope))        
    # plot time-history
    jd0_obs = solution_dict["jd0_obs"][idx]
    axes = [Axis(fig_time[1,idx]; xlabel="Time, hour", ylabel="Azimuth, deg"),
            Axis(fig_time[2,idx]; xlabel="Time, hour", ylabel="Elevation, deg")]
    TelescopeTasking.plot_time_history!(axes, passes; jd_ref=jd0_obs, color=:grey70, linewidth=0.1)
    TelescopeTasking.plot_time_history!(axes, selected_passes; jd_ref=jd0_obs, 
        linewidth=0.8, color_by_target=true, exposure_only=true)
end


# # plot of non-exposure distance between passes
# fig_L = Figure(size=(500,500))
# ax_L = Axis(fig_L[1,1]; xlabel="Pass index", ylabel="Non-exposure distance, deg")
# _, Ls = TelescopeTasking.get_nonexposure_distance(passes, Y)
# Ls_nonzero = [L for L in Ls if L > 0]
# scatter!(ax_L, 1:length(Ls_nonzero), rad2deg.(Ls_nonzero); color=:grey10, markersize=5)

# # display figure
display(fig_sol)
# display(fig_time)
# # display(fig_L)