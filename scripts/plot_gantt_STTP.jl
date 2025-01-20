"""Make gantt chart plot for STTP"""


using GLMakie
using JSON
using ProgressMeter: @showprogress
using Printf: @printf, @sprintf
using SatelliteToolboxTle
using SatelliteToolboxSgp4
using SatelliteToolboxTransformations

GLMakie.activate!(;)

include(joinpath(@__DIR__, "../src/TelescopeTasking.jl"))


# load Earth parameters
eop_file = joinpath(@__DIR__, "..", "data", "eop_iau1980", "finals.all.csv")
eop_iau1980 = read_iers_eop(eop_file, Val(:IAU1980))

# load config jsons
telescope = JSON.parsefile(joinpath(@__DIR__, "configs/config_telescope.json"))
config = JSON.parsefile(joinpath(@__DIR__, "configs/config_STTP1.json"))
target_choice = "A"
num_exposure = 1
solver_choice = "Gurobi"
experiment_name = config["name"] * "_target$(target_choice)_E$(num_exposure)"
@printf("Plotting pass from experiment %s\n", experiment_name)

# load tles used
tles = read_tles(read(joinpath(@__DIR__, "../data/tles/AAS25target$(target_choice).txt"), String))

# load solutions
solution_dict = JSON.parsefile(joinpath(@__DIR__, "solutions/solution_$(experiment_name)_$(solver_choice).json"))
jd0_obs = solution_dict["jd0_obs"]
obs_duration = solution_dict["obs_duration"]
min_elevation = solution_dict["min_elevation"]
min_obs_duration = solution_dict["min_obs_duration"]
exposure_duration = solution_dict["exposure_duration"]
observer_lla = solution_dict["observer_lla"]

# get selected passes from solution 
passes = [TelescopeTasking.VisiblePass(pass_dict) for pass_dict in solution_dict["passes_dict"]]
Y = [el > 0.5 for el in solution_dict["Y"]]
selected_passes = [pass for (pass, y) in zip(passes, value.(Y)) if y > 0.5]


function schedule_band(ax, color, alpha=0.3)
    band!(ax, xs, ylow, yhigh, color=color, alpha=alpha)
end