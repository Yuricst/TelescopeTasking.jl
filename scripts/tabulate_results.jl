"""Tabulate STTP results"""

using CairoMakie
using JSON
using Printf: @printf
using SatelliteToolboxTransformations

include(joinpath(@__DIR__, "../src/TelescopeTasking.jl"))

CairoMakie.activate!()

function main()
    # load Earth parameters
    eop_file = joinpath(@__DIR__, "..", "data", "eop_iau1980", "finals.all.csv")
    eop_iau1980 = read_iers_eop(eop_file, Val(:IAU1980))

    # load config jsons
    instance_names = ["STTP1", "MTTP1", "MTTP2", "MTTP3", "MTTP4"]  #"STTP2", "MTTP1", "MTTP2", "MTTP3"]
    @show target_choice = "A"

    experiment_name_dict = Dict(
        "STTP1" => "STTP I",
        "STTP2" => "STTP II",
        "MTTP1" => "MTTP I",
        "MTTP2" => "MTTP II",
        "MTTP3" => "MTTP III",
        "MTTP4" => "MTTP IV",
    )

    for instance_name in instance_names
        config_telescope = JSON.parsefile(joinpath(@__DIR__, "configs/config_telescope.json"))
        config = JSON.parsefile(joinpath(@__DIR__, "configs/config_$(instance_name).json"))
                
        slew_rate = deg2rad(config_telescope["slew_rate"])                         # rad/s
        exposure_duration = config_telescope["exposure_duration"]              # in seconds
        for num_exposure in [1,2,3]
            # if instance_name == "MTTP1" && num_exposure == 4
            #     continue 
            # end
            
            _solver_choice = "Gurobi"
            _experiment_name = config["name"] * "_target$(target_choice)_E$(num_exposure)"

            # load solutions
            solution_dict = JSON.parsefile(joinpath(@__DIR__, "solutions/solution_$(_experiment_name)_$(_solver_choice).json"))
            solve_stats = JSON.parsefile(joinpath(@__DIR__, "solutions/solve_stats_$(_experiment_name)_$(_solver_choice).json"))

            # instance name, with multicolumn
            if num_exposure == 1
                if instance_name == "MTTP4"
                    n_row = 5
                else
                    n_row = 3
                end
                instance_name_disp = "\\multirow{3}{*}{$(experiment_name_dict[instance_name])}"
            else
                instance_name_disp = "                   "
            end

            # print results
            if occursin("STTP", instance_name)
                if (sum(solution_dict["X"]) > 0)
                    passes = [TelescopeTasking.VisiblePass(pass_dict) for pass_dict in solution_dict["passes_dict"]]
                    tau_idle = TelescopeTasking.get_idletimeratio(
                        passes,
                        solution_dict["Y"],
                        solution_dict["obs_duration"],
                        exposure_duration,
                        slew_rate,
                    )
                    @printf("    %s & \$%1.0f\$ & \$%1.3f\$ & \$%1.2f\$ & \$%1.0f\$ & \$%1.2f\$ & \$%1.2f\$ \\\\\n",
                        instance_name_disp,
                        num_exposure,
                        solve_stats["solve_time"],
                        solve_stats["relative_gap"] * 100,
                        sum(solution_dict["X"]),
                        sum(solution_dict["eta"]),
                        # sum(solution_dict["L"]),
                        tau_idle,
                    )
                else
                    @printf("    %s & \$%1.0f\$ & n.a. & n.a. & n.a. & n.a. & n.a. \\\\\n",
                        instance_name_disp,
                        num_exposure,
                    )
                end
            elseif occursin("MTTP", instance_name)
                s = length(solution_dict["eta_per_telescope"])
                eta_per_telescope = solution_dict["eta_per_telescope"]
                L_per_telescope = [sum(Ls) for Ls in solution_dict["L_per_telescope"]]
                passes_per_telescope = [
                    [TelescopeTasking.VisiblePass(pass_dict) for pass_dict in passes_dict_list]
                    for passes_dict_list in solution_dict["passes_per_telescope"]
                ]       
                tau_idle_per_telescope = [
                    tau_idle = TelescopeTasking.get_idletimeratio(
                        passes_per_telescope[q],
                        solution_dict["Y_per_telescope"][q],
                        solution_dict["obs_duration_per_telescope"][q],
                        exposure_duration,
                        slew_rate,
                    )
                    for q = 1:length(eta_per_telescope)
                ]
                if s == 2
                    @printf("    %s & \$%1.0f\$ & \$%1.3f\$ & \$%1.2f\$ & \$%1.0f\$ & \$%1.2f\$, \$%1.2f\$ & \$%1.2f\$, \$%1.2f\$ \\\\\n",
                        instance_name_disp,
                        num_exposure,
                        solve_stats["solve_time"],
                        solve_stats["relative_gap"] * 100,
                        sum(solution_dict["X"]),
                        eta_per_telescope...,
                        # L_per_telescope...,
                        tau_idle_per_telescope...
                    )
                elseif s == 3
                    @printf("    %s & \$%1.0f\$ & \$%1.3f\$ & \$%1.2f\$ & \$%1.0f\$ & \$%1.2f\$, \$%1.2f\$, \$%1.2f\$ & \$%1.2f\$, \$%1.2f\$, \$%1.2f\$ \\\\\n",
                        instance_name_disp,
                        num_exposure,
                        solve_stats["solve_time"],
                        solve_stats["relative_gap"] * 100,
                        sum(solution_dict["X"]),
                        eta_per_telescope...,
                        # L_per_telescope...,
                        tau_idle_per_telescope...
                    )
                elseif s == 4
                    @printf("    %s & \$%1.0f\$ & \$%1.3f\$ & \$%1.2f\$ & \$%1.0f\$ & \\begin{tabular}[c]{@{}l@{}}\$%1.2f\$, \$%1.2f\$,\\\\\$%1.2f\$, \$%1.2f\$\\end{tabular} & \\begin{tabular}[c]{@{}l@{}}\$%1.2f\$, \$%1.2f\$,\\\\\$%1.2f\$, \$%1.2f\$ \\end{tabular} \\\\\n",
                        instance_name_disp,
                        num_exposure,
                        solve_stats["solve_time"],
                        solve_stats["relative_gap"] * 100,
                        sum(solution_dict["X"]),
                        eta_per_telescope...,
                        # L_per_telescope...,
                        tau_idle_per_telescope...
                    )
                end
            else
                @printf("    %s & \$%1.0f\$ & n.a. & n.a. & n.a. & n.a. & n.a. \\\\\n",
                    experiment_name_dict[instance_name],
                    num_exposure,
                )
            end
        end
        if instance_name != instance_names[end]
            @printf("    \\midrule \n")
        else
            @printf("    \\bottomrule \n")
        end
    end

    # create figure
    println("Creating plot...")
    fontsize = 24
    fig_observed = Figure(size=(500,400))
    axis = Axis(fig_observed[1,1];
        xlabel="Number of exposures", ylabel="Number of observed targets",
        xlabelsize=fontsize, ylabelsize=fontsize,
        xticklabelsize=fontsize-2, yticklabelsize=fontsize-2)
    axis.xticks = [1,2,3]
    
    _solver_choice = "Gurobi"

    markers = [:circle, :rect, :diamond, :vline, :hline]
    colors = cgrad(:hawaii, length(instance_names)+1, categorical=true)
    for (idx,instance_name) in enumerate(instance_names)
        config = JSON.parsefile(joinpath(@__DIR__, "configs/config_$(instance_name).json"))
        sum_X_per_E = Real[]
        for num_exposure in [1,2,3]
            _experiment_name = config["name"] * "_target$(target_choice)_E$(num_exposure)"
            solution_dict = JSON.parsefile(joinpath(@__DIR__, "solutions/solution_$(_experiment_name)_$(_solver_choice).json"))
            push!(sum_X_per_E, sum(sum(solution_dict["X"])))
        end
        scatterlines!(axis, [1,2,3], sum_X_per_E; 
            label=experiment_name_dict[instance_name],
            markersize=18,
            marker=markers[idx],
            color=colors[idx],)
    end
    if target_choice == "A"
        axislegend(axis, position=:lb)
    else
        axislegend(axis, position=:rt)
    end
    save(joinpath(@__DIR__, "plots", "NumObserved_vs_E_target$(target_choice).pdf"), fig_observed)
    display(fig_observed)
end

main()