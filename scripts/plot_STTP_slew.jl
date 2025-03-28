"""Plot STTP solution slew path"""

using CairoMakie
using Colors
using ColorSchemes
using GeometryBasics
using GLMakie
using JSON
using JuMP
using LinearAlgebra
using ProgressMeter: @showprogress
using Printf: @printf
using SatelliteToolboxTle
using SatelliteToolboxSgp4
using SatelliteToolboxTransformations

CairoMakie.activate!()

include(joinpath(@__DIR__, "../src/TelescopeTasking.jl"))

function main()
    # load Earth parameters
    eop_file = joinpath(@__DIR__, "..", "data", "eop_iau1980", "finals.all.csv")
    eop_iau1980 = read_iers_eop(eop_file, Val(:IAU1980))

    # load config jsons
    telescope = JSON.parsefile(joinpath(@__DIR__, "configs/config_telescope.json"))
    config = JSON.parsefile(joinpath(@__DIR__, "configs/config_STTP1.json"))
    target_choice = "A"    # A or S1
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

    passes = [TelescopeTasking.VisiblePass(pass_dict) for pass_dict in solution_dict["passes_dict"]]
    # # recreate passes for making sure it is correctly ordered
    # passes2, _ = TelescopeTasking.tles_to_passes(
    #     tles,
    #     eop_iau1980,
    #     jd0_obs,
    #     obs_duration,
    #     min_elevation,
    #     min_obs_duration,
    #     exposure_duration,
    #     observer_lla,
    #     dt_sec=10,
    #     num_exposure = num_exposure,
    # )
    # for (p1,p2) in zip(passes, passes2)
    #     @assert p1.tle == p2.tle
    #     @assert p1.times == p2.times
    #     @assert p1.azimuths == p2.azimuths
    #     @assert p1.elevations == p2.elevations
    # end

    # get selected passes from solution 
    Y = [el > 0.5 for el in solution_dict["Y"]]
    selected_passes = [pass for (pass, y) in zip(passes, value.(Y)) if y > 0.5]
    @show length(selected_passes)
    figures = Figure[]

    i_choose_list = [1, 15, 30, 45] #, 60, 75]
    N_passes_to_plot_list = [15,16,16,16]
    for j in 1:length(i_choose_list)
        i_choose = i_choose_list[j]
        N_passes_to_plot = N_passes_to_plot_list[j]
        final_index = min(i_choose + N_passes_to_plot - 1, length(selected_passes))
        selected_passes_plot = selected_passes[i_choose:final_index]
        colors = cgrad(:hawaii, N_passes_to_plot, categorical = true)

        @show i_choose, final_index

        # plot of selected passes
        fontsize = 24
        fig_sol = Figure(size=(600,500);)
        ax_sol = PolarAxis(fig_sol[1,1];
            rticks = ([0,30,60], ["90","60","30"]),
            # rlabelsize=fontsize,
            # thetalabelsize=fontsize,
            rgridwidth = 0.3,
            thetagridwidth = 0.3,
            rticklabelsize=fontsize-1,
            thetaticklabelsize=fontsize-1,)
        TelescopeTasking.polar_plot_passes!(ax_sol, selected_passes_plot; color=:grey50, linewidth=0.25)
        TelescopeTasking.polar_plot_passes!(ax_sol, selected_passes_plot; 
            linewidth=3, color=colors, color_by_target=false, exposure_only=true)
        TelescopeTasking.polar_plot_interpass_slew!(ax_sol, selected_passes_plot, colors; linestyle = :dash)
        Colorbar(fig_sol[1,2], limits = (0.5, N_passes_to_plot+0.5), 
            colormap = cgrad(:hawaii, N_passes_to_plot, categorical = true),
            labelsize = fontsize, ticklabelsize = fontsize,
            ticks = (1:N_passes_to_plot, string.(i_choose:final_index)),
            label = "Pass number")
        filepath = joinpath(@__DIR__, "plots", "slew_STTP_target$(target_choice)_start$j.png")
        CairoMakie.save(filepath, fig_sol)
        @show filepath
        push!(figures, fig_sol)
    end
    return figures
end 

figures = main()
display(figures[1])
println("Done!")