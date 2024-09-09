"""Tabulate STTP results"""

using CairoMakie
using JSON
using Printf: @printf
using SatelliteToolboxTransformations

include(joinpath(@__DIR__, "../src/TelescopeTasking.jl"))

# load Earth parameters
eop_file = joinpath(@__DIR__, "..", "data", "eop_iau1980", "finals.all.csv")
eop_iau1980 = read_iers_eop(eop_file, Val(:IAU1980))

# load config jsons
instance_names = ["STTP1", "STTP2", "MTTP1", "MTTP2", "MTTP3"]


function compute_gap(Z_opt, Z_greedy)
    return (sum(Z_opt) - sum(Z_greedy)) / sum(Z_opt)
end


for instance_name in instance_names
    telescope = JSON.parsefile(joinpath(@__DIR__, "configs/config_telescope.json"))
    config = JSON.parsefile(joinpath(@__DIR__, "configs/config_$(instance_name).json"))
    target_choice = "B"

    experiment_name_dict = Dict(
        "STTP1" => "STTP I",
        "STTP2" => "STTP II",
        "MTTP1" => "MTTP I",
        "MTTP2" => "MTTP II",
        "MTTP3" => "MTTP III",
    )

    for num_exposure in [1]
        _experiment_name = config["name"] * "_target$(target_choice)_E$(num_exposure)"

        # load solutions
        solution_dict = JSON.parsefile(joinpath(@__DIR__, "solutions/solution_$(_experiment_name)_greedy.json"))
        gurobi_solution_dict = JSON.parsefile(joinpath(@__DIR__, "solutions/solution_$(_experiment_name)_gurobi.json"))


        # print results
        if occursin("STTP", instance_name)
            if (sum(solution_dict["X"]) > 0)
                @printf("    %s & \$%1.0f\$ & \$%1.2f\$ & \$%1.2f\$ & \$%1.0f\$ & \$%1.2f\$ & \$%1.1f\$ \\\\\n",
                    experiment_name_dict[instance_name],
                    num_exposure,
                    solution_dict["solve_time"],
                    compute_gap(sum(gurobi_solution_dict["X"]), sum(solution_dict["X"])) * 100,
                    sum(solution_dict["X"]),
                    sum(solution_dict["eta"]),
                    sum(solution_dict["L"]),
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
            if s == 2
                @printf("    %s & \$%1.0f\$ & \$%1.2f\$ & \$%1.2f\$ & \$%1.0f\$ & \$%1.2f\$, \$%1.2f\$ & \$%1.1f\$, \$%1.1f\$ \\\\\n",
                    experiment_name_dict[instance_name],
                    num_exposure,
                    solution_dict["solve_time"],
                    compute_gap(sum(gurobi_solution_dict["X"]), sum(solution_dict["X"])) * 100,
                    sum(solution_dict["X"]),
                    eta_per_telescope...,
                    L_per_telescope...
                )
            elseif s == 3
                @printf("    %s & \$%1.0f\$ & \$%1.2f\$ & \$%1.2f\$ & \$%1.0f\$ & \$%1.2f\$, \$%1.2f\$, \$%1.2f\$ & \$%1.1f\$, \$%1.1f\$, \$%1.1f\$ \\\\\n",
                    experiment_name_dict[instance_name],
                    num_exposure,
                    solution_dict["solve_time"],
                    compute_gap(sum(gurobi_solution_dict["X"]), sum(solution_dict["X"])) * 100,
                    sum(solution_dict["X"]),
                    eta_per_telescope...,
                    L_per_telescope...
                )
            end
        else
            @printf("    %s & \$%1.0f\$ & n.a. & n.a. & n.a. & n.a. & n.a. \\\\\n",
                experiment_name_dict[instance_name],
                num_exposure,
            )
        end
    end
end