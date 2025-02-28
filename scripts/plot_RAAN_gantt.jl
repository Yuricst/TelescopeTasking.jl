"""Plot of Gantt-chart-like plot with vertical axis in terms of RAAN"""

using Colors
using ColorSchemes
using GeometryBasics
using GLMakie
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
target_choice = "S1"
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
fig_gantt = Figure(size=(800,400))
ax = Axis(fig_gantt[1,1];
    xlabelsize = fontsize,
    ylabelsize = fontsize,
    xticklabelsize = fontsize-2,
    yticklabelsize = fontsize-2,
    xlabel="Elapsed time, hours", 
    ylabel="RAAN, deg",
    yreversed=true
)
colors = [:blue, :red]

for q in 1:length(passes_per_telescope)
    selected_passes = selected_passes_per_telescope[q]

    # extract times and raans
    gantt_x = []
    gantt_y = []
    for pass in selected_passes
        push!(gantt_x, [pass.t0_exposure, pass.tf_exposure] .- jd0_obs)
        push!(gantt_y, [pass.tle.raan, pass.tle.raan])
    end

    # create plot
    for (idx,(x,y)) in enumerate(zip(gantt_x, gantt_y))
        if idx == 1
            scatterlines!(ax, x*24, y, color=colors[q], label="Telescope $q")
        else
            scatterlines!(ax, x*24, y, color=colors[q])
        end
    end
end
filename = "gantt_raan_$(experiment_name)_$(solver_choice).png"
# save(joinpath(@__DIR__, "plots", filename), fig_gantt)
display(fig_gantt)