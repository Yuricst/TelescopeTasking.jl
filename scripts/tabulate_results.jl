"""Tabulate STTP results"""

using CairoMakie
using JSON
using Printf: @printf
using SatelliteToolboxTransformations

include(joinpath(@__DIR__, "../src/TelescopeTasking.jl"))

function main()
    # load Earth parameters
    eop_file = joinpath(@__DIR__, "..", "data", "eop_iau1980", "finals.all.csv")
    eop_iau1980 = read_iers_eop(eop_file, Val(:IAU1980))

    # load config jsons
    instance_names = ["STTP1", "STTP2", "MTTP1", "MTTP2", "MTTP3"]
    for instance_name in instance_names
        config_telescope = JSON.parsefile(joinpath(@__DIR__, "configs/config_telescope.json"))
        config = JSON.parsefile(joinpath(@__DIR__, "configs/config_$(instance_name).json"))
        target_choice = "B"
        
        slew_rate = deg2rad(config_telescope["slew_rate"])                         # rad/s
        exposure_duration = config_telescope["exposure_duration"]              # in seconds

        experiment_name_dict = Dict(
            "STTP1" => "STTP I",
            "STTP2" => "STTP II",
            "MTTP1" => "MTTP I",
            "MTTP2" => "MTTP II",
            "MTTP3" => "MTTP III",
        )

        for num_exposure in [1,2,3]
            _solver_choice = "Gurobi"
            _experiment_name = config["name"] * "_target$(target_choice)_E$(num_exposure)"

            # load solutions
            solution_dict = JSON.parsefile(joinpath(@__DIR__, "solutions/solution_$(_experiment_name)_$(_solver_choice).json"))
            solve_stats = JSON.parsefile(joinpath(@__DIR__, "solutions/solve_stats_$(_experiment_name)_$(_solver_choice).json"))

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
                    @printf("    %s & \$%1.0f\$ & \$%1.2f\$ & \$%1.2f\$ & \$%1.0f\$ & \$%1.2f\$ & \$%1.2f\$ \\\\\n",
                        experiment_name_dict[instance_name],
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
                        experiment_name_dict[instance_name],
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
                    @printf("    %s & \$%1.0f\$ & \$%1.2f\$ & \$%1.2f\$ & \$%1.0f\$ & \$%1.2f\$, \$%1.2f\$ & \$%1.2f\$, \$%1.2f\$ \\\\\n",
                        experiment_name_dict[instance_name],
                        num_exposure,
                        solve_stats["solve_time"],
                        solve_stats["relative_gap"] * 100,
                        sum(solution_dict["X"]),
                        eta_per_telescope...,
                        # L_per_telescope...,
                        tau_idle_per_telescope...
                    )
                elseif s == 3
                    @printf("    %s & \$%1.0f\$ & \$%1.2f\$ & \$%1.2f\$ & \$%1.0f\$ & \$%1.2f\$, \$%1.2f\$, \$%1.2f\$ & \$%1.2f\$, \$%1.2f\$, \$%1.2f\$ \\\\\n",
                        experiment_name_dict[instance_name],
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
        @printf("    \\midrule \n")
    end
end

main()