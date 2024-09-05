"""
Handling TLEs
"""

using Colors
using ColorSchemes
using GeometryBasics
using GLMakie
using GLPK
using Gurobi
using JuMP
using LinearAlgebra
using ProgressMeter: @showprogress
using Printf: @printf
using SatelliteToolboxTle
using SatelliteToolboxSgp4
using SatelliteToolboxTransformations

include(joinpath(@__DIR__, "../src/TelescopeScheduling.jl"))


function dcm_ned_to_enu()
    T3 = [0 1 0;
          -1 0 0;
          0 0 1];
    T1 = [1 0 0;
          0 -1 0;
          0 0 -1];
    return T1 * T3
end

function cart2sph(rvec)
    az = atan(rvec[2], rvec[1])
    el = atan(rvec[3], sqrt(rvec[1]^2 + rvec[2]^2))
    r = norm(rvec)
    return [az; el; r]
end

function periodic_lines!(ax, xs, ys, diff_max; color=:red, linewidth=0.5)
    i_from = 1
    for i in 1:length(xs)-1
        diff = abs(xs[i] - xs[i+1])
        if diff > diff_max
            lines!(ax, xs[i_from:i], ys[i_from:i], color=color, linewidth=linewidth)
            i_from = i + 1
        end
    end
end

T_NED2ENU = dcm_ned_to_enu()

# observer location in geodetic coordinates
lat = 45.5
lon = 135.3
alt = 100.0         # in meters
r_ECEF_obs = geodetic_to_ecef(lat |> deg2rad, lon |> deg2rad, alt) / 1e3    # in km
r_ENU_obs = T_NED2ENU * ecef_to_ned(r_ECEF_obs, lat |> deg2rad, lon |> deg2rad, alt; translate = false)

# initial epoch of local nightfall
jd0_obs = 2.46055755221534e6 + 0.45     # in julian date

# criteria for valid passes
min_elevation = deg2rad(15)
min_duration_day = 60/86400        # minimum duration: 1 minute, in days
exposure_duration_day = 45/86400   # exposure duration: 45 seconds, in days

# load Earth parameters
# eop_iau1980 = fetch_iers_eop()
eop_file = joinpath(@__DIR__, "..", "data", "eop_iau1980", "finals.all.csv")
eop_iau1980 = read_iers_eop(eop_file, Val(:IAU1980))

# load TLE files
# path_to_tles = joinpath(@__DIR__, "..", "data", "tles", "iridium-33-debris.txt")
path_to_tles = joinpath(@__DIR__, "..", "data", "tles", "active.txt")
tles_str = read(path_to_tles, String)

# convert to TLE objects
tles = read_tles(tles_str)
@show length(tles)

# filter them
names_include = ["GLOBALSTAR", "IRIDIUM"]
TelescopeScheduling.filter!(tles, names_include=names_include)
@show length(tles)

# get passes
obs_duration = 8 * 3600             # in seconds
min_elevation = deg2rad(15)
min_obs_duration = 100              # in seconds
exposure_duration = 60              # in seconds
observer_lat = deg2rad(45)
observer_lon = deg2rad(100)
observer_alt = 30.0
observer_lla = [observer_lat, observer_lon, observer_alt]

passes, sph_ENU_list = TelescopeScheduling.tles_to_passes(
    tles[1:100],
    eop_iau1980,
    jd0_obs,
    obs_duration,
    min_elevation,
    min_obs_duration,
    exposure_duration,
    observer_lla,
    dt_sec=10,
)
@show length(passes)
smas = [TelescopeScheduling.tle2sma(pass.tle) for pass in passes]

# construct problem
num_exposure = 2
slew_rate = deg2rad(2)      # rad/s
buffer_times = [15, 0]      # times in seconds
problem = TelescopeScheduling.TelescopeSchedulingProblem(
    passes, num_exposure,
    slew_rate;
    buffer_times = buffer_times
)
@show problem;

# solve problem
solver = MOI.OptimizerWithAttributes(Gurobi.Optimizer,
    "TimeLimit" => 30)
X, Y = TelescopeScheduling.solve!(problem, solver)
selected_passes = [pass for (pass, y) in zip(passes, value.(Y)) if y > 0.5]
@printf("Number of observed targets: %d\n", sum(value.(X)))

# plot of selected passes
fig_sol = Figure(size=(600,600))
ax_sol = PolarAxis(fig_sol[1,1])
TelescopeScheduling.polar_plot_passes!(ax_sol, passes; color=:grey, linewidth=1.0)
TelescopeScheduling.polar_plot_passes!(ax_sol, selected_passes; linewidth=1.5, color_by_target=true)


# plot time-history
fig_sol2 = Figure(size=(800, 600))
axes = [Axis(fig_sol2[1,1]; xlabel="Time, hour", ylabel="Azimuth, deg"),
        Axis(fig_sol2[2,1]; xlabel="Time, hour", ylabel="Elevation, deg")]
TelescopeScheduling.plot_time_history!(axes, passes; jd_ref=jd0_obs, color=:grey, linewidth=1.0)
TelescopeScheduling.plot_time_history!(axes, selected_passes; jd_ref=jd0_obs, linewidth=1.5, color_by_target=true)

display(fig_sol2)