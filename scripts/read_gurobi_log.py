"""Ready and plot iteration results from gurobi log file"""

import os
import numpy as np
import matplotlib.pyplot as plt
import re


fontsize = 15
plt.rc('font', size=fontsize)          # controls default text sizes


def process_line(line):
    try:
        nodes_expl = float(line[1:7])
    except:
        return None
    
    # Match rows with numeric values in the Incumbent and BestBd columns
    #match = re.search(r'(\d+\.\d+|-) +(\d+\.\d+)', line)
    match = re.search(r'(\d+\.\d+|-) +(\d+\.\d+) +(\d+\.\d+)%', line)
    if match:
        incumbent = match.group(1)
        best_bd = match.group(2)
        gap = match.group(3)
        # Convert to float if numeric, otherwise None
        incumbent = float(incumbent) if incumbent != '-' else None
        best_bd = float(best_bd)
        gap = float(gap) if gap != '-' else None

    return nodes_expl, incumbent, best_bd, gap


def read_gurobi_log(filename, ylabel=None): #p, N_demand, size_depot, size_roundtrip, case_name=2):
    """Read gurobi log file and return iteration results"""
    with open(filename, "r") as f:
            lines = f.readlines()

    if ylabel is None:
        ylabel = "Objective value"

    # find iteration results
    for i, line in enumerate(lines):
        if " Expl Unexpl |  Obj  Depth IntInf | Incumbent    BestBd   Gap | It/Node Time" in line:
            i_start = i
            break

    print(f"Found iteration results at line {i_start}")
    iter_results = []
    for line in lines[i_start+3:]:
        if "Explored" in line:
            break
        line_res = process_line(line)
        if line_res is None:
            break
        iter_results.append(line_res)
    iter_results = np.array(iter_results)
    _, idx_safe = np.unique(iter_results[:,0], return_index=True)
    iter_results_safe = np.array(iter_results[idx_safe,:])
    
    # fill needs arrays that are not constructed in an 'ugly' way so we create new arrays
    fill_x, fill_y1, fill_y2 = np.zeros(len(iter_results_safe[:,0])), np.zeros(len(iter_results_safe[:,0])), np.zeros(len(iter_results_safe[:,0]))
    for idx in range(len(iter_results_safe[:,0])):
        fill_x[idx] = iter_results_safe[idx,0]
        fill_y1[idx] = iter_results_safe[idx,1]
        fill_y2[idx] = iter_results_safe[idx,2]

    # plot upper & lower bound
    fig_bounds, ax = plt.subplots(figsize=(6,3.2))
    ax.plot(iter_results[:,0]/1e3, iter_results[:,2], marker="+", color="red", label="Best Bound")
    ax.plot(iter_results[:,0]/1e3, iter_results[:,1], marker="x", color="blue", label="Incumbent")
    ax.fill_between(fill_x/1e3, fill_y1, fill_y2,
                    color="grey", alpha=0.25, label="Optimality gap")
    ax.set(xlabel=r"Number of nodes explored, $\times 10^3$", ylabel=ylabel)
    ax.set(yscale="log")
    ax.legend()
    plt.tight_layout()
    # fig.savefig(
    #     os.path.join(os.path.dirname(__file__),
    #                  f"plots/grb_iter_case{case_name}_p{p}_N{N_demand}.png"),
    #     dpi = 300)
    
    # plot gap
    fig_gap, ax = plt.subplots(figsize=(6,3.2))
    ax.grid(True, alpha=0.5)
    ax.plot(iter_results[:,0]/1e3, iter_results[:,3], marker="x", color="black")
    ax.set(xlabel=r"Number of nodes explored, $\times 10^3$", ylabel="Optimality gap, %")
    if max(iter_results[:,3]) >= 80:
        ax.set(yticks=[0,25,50,75,100])
    plt.tight_layout()
    # fig.savefig(
    #     os.path.join(os.path.dirname(__file__),
    #                  f"plots/grb_gap_case{case_name}_p{p}_N{N_demand}.png"),
    #     dpi = 300)
    return fig_bounds, fig_gap

if __name__ == "__main__":
    import argparse
    parser  = argparse.ArgumentParser()
    parser.add_argument('--instance', '-i', type=str, default="MTTP4")
    parser.add_argument('--target', '-t', type=str, default="S2")
    parser.add_argument('--exposure', '-e', type=int, default=2)
    args = parser.parse_args()

    filename = os.path.join(
        os.path.dirname(__file__),
        f"solutions_JASS/logs/",
        f"log_{args.instance}_target{args.target}_E{args.exposure}_Gurobi.log",
    )
    #\scripts\solutions_JASS\logs\log_MTTP4_targetS2_E2_Gurobi.log
    fig_bounds, fig_gap = read_gurobi_log(filename)

    # save to file
    fig_bounds.savefig(
        os.path.join(os.path.dirname(__file__),
                     f"plots/logs/",
                     f"grb_bnds_{args.instance}_target{args.target}_E{args.exposure}.png"),
        dpi = 300)
    fig_gap.savefig(
        os.path.join(os.path.dirname(__file__),
                     f"plots/logs/",
                     f"grb_gap_{args.instance}_target{args.target}_E{args.exposure}.png"),
        dpi = 300)
    plt.show()