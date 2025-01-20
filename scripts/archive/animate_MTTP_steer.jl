"""Animate MTTP"""

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
using Printf: @printf, @sprintf
using SatelliteToolboxTle
using SatelliteToolboxSgp4
using SatelliteToolboxTransformations

include(joinpath(@__DIR__, "../src/TelescopeTasking.jl"))

# ------------------------------------------------------------------------ #
# load Earth parameters
eop_file = joinpath(@__DIR__, "..", "data", "eop_iau1980", "finals.all.csv")
eop_iau1980 = read_iers_eop(eop_file, Val(:IAU1980))

# ------------------------------------------------------------------------ #
# load config jsons
telescope = JSON.parsefile(joinpath(@__DIR__, "configs/config_telescope.json"))
config = JSON.parsefile(joinpath(@__DIR__, "configs/config_MTTP1.json"))
target_choice = "A"
num_exposure = 1
solver_choice = "Gurobi"
experiment_name = config["name"] * "_target$(target_choice)_E$(num_exposure)"
@printf("Plotting pass from experiment %s\n", experiment_name)

# ------------------------------------------------------------------------ #
# load tles used
tles = read_tles(read(joinpath(@__DIR__, "../data/tles/AAS25target$(target_choice).txt"), String))

# ------------------------------------------------------------------------ #
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



# ------------------------------------------------------------------------ #
# create animator
save_to_file = false
q = 1
timestamps = [0.0, 1.0, 2.0]

# plot of allocated solution
with_theme(theme_light()) do
# with_theme(theme_dark()) do
    time = Observable(0.0)

    # compute time-step idx s.t. time_steps[idx] <= time < time_steps[idx+1]
    function time2timestep(time)
        return max(1, Int(ceil(time/dt)))
    end

    # initialize figure
    fontsize = 24
    fig = Figure(size = (500, 500))
    _ax_polar = PolarAxis(
        fig[1,1];
        title = @lift("Elapsed time = $(@sprintf "%1.2f" $time )"),
        rticks = ([0,20,40,60], ["90","70","50","30"]),
        rticklabelsize=fontsize-1,
        thetaticklabelsize=fontsize-1,
    )
    TelescopeTasking.polar_plot_passes!(_ax_polar, passes_per_telescope[q]; color=:grey60, linewidth=0.1)
    # TelescopeTasking.polar_plot_passes!(_ax_polar, selected_passes_per_telescope[q]; 
    #     linewidth=1.5, color_by_target=true, exposure_only=true)

    # plot passes in order
    exposures_per_pass = []
    for pass in selected_passes_per_telescope[q]
        _, az_exposure, el_exposure = TelescopeTasking.get_exposure_history(pass)
        push!(exposures_per_pass, )
    end

    # # plot static targets
    # static_target_states = hcat([
    #     [el[1][1], el[1][2], el[1][3]] for el in targets_states
    # ]...)
    # for k = 1:nk
    #     color_observable = @lift(
    #         if maximum(M_bin[:,:,time2timestep($time), k] .* x_[:,:,time2timestep($time)]) == 1
    #             RGBf(0.53, 0.81, 0.92)
    #         else
    #             RGBf(1.0, 0.25, 0.25)
    #         end
    #     )
    #     if k == 1
    #         label = "static target"
    #     else
    #         label = nothing
    #     end
    #     scatter!(_ax_polar,
    #             [targets_states[k][1][1]], [targets_states[k][1][2]], [targets_states[k][1][3]],
    #             color=color_observable,
    #             markersize=8,
    #             label = label,
    #             marker=:circle)
    # end

    if save_to_file == true
        # animate to file
        record(
            fig,
            joinpath(
                save_dir,
                "csda_pmedian_$(target_set)_p$(p)M$(Mcut)_$(method).gif"
            ),
            animator.timestamps;
            framerate = animator.framerate) do t
            animator.time[] = t
            sleep(0.05)
        end

    else 
        # live display 
        display(fig)
        for t in timestamps
            time[] = t
            sleep(0.05)
            yield() #-> not required with record
        end
    end
end