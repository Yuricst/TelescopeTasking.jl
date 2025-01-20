"""Create animation of steering history of telescopes"""

using Dates
using GLMakie
using JSON
using ProgressMeter: @showprogress
using Printf: @printf, @sprintf
using SatelliteToolboxTle
using SatelliteToolboxSgp4
using SatelliteToolboxTransformations

GLMakie.activate!(;)

include(joinpath(@__DIR__, "../src/TelescopeTasking.jl"))

# ------------------------------------------------------------------------ #
# load Earth parameters
eop_file = joinpath(@__DIR__, "..", "data", "eop_iau1980", "finals.all.csv")
eop_iau1980 = read_iers_eop(eop_file, Val(:IAU1980))

# ------------------------------------------------------------------------ #
# load config jsons
telescope = JSON.parsefile(joinpath(@__DIR__, "configs/config_telescope.json"))
config_name = "config_STTP1"
config = JSON.parsefile(joinpath(@__DIR__, "configs/$(config_name).json"))
target_choice = "A"
num_exposure = 1
solver_choice = "Gurobi"
experiment_name = config["name"] * "_target$(target_choice)_E$(num_exposure)"
@printf("Plotting pass from experiment %s\n", experiment_name)

# load tles used
tles = read_tles(read(joinpath(@__DIR__, "../data/tles/AAS25target$(target_choice).txt"), String))
# tles = tles[1:10]        # FIXME -- reducing for debugging

# load solutions
solution_dict = JSON.parsefile(joinpath(@__DIR__, "solutions/solution_$(experiment_name)_$(solver_choice).json"))
jd0_obs = solution_dict["jd0_obs"]
obs_duration = solution_dict["obs_duration"]
# min_elevation = solution_dict["min_elevation"]
# min_obs_duration = solution_dict["min_obs_duration"]
# exposure_duration = solution_dict["exposure_duration"]
# observer_lla = solution_dict["observer_lla"]

passes = [TelescopeTasking.VisiblePass(pass_dict) for pass_dict in solution_dict["passes_dict"]]
Y = [el > 0.5 for el in solution_dict["Y"]]
selected_passes = [pass for (pass, y) in zip(passes, value.(Y)) if y > 0.5]

# # get exposure azimtuhs & elevations for selected_passes
# selected_exp_azel = []
# for (ipass,pass) in enumerate(selected_passes)
#     _,az_exp,el_exp = TelescopeTasking.get_exposure_history(pass)
#     selected_exp_azel = [az_exp el_exp]
# end

# get list of observed targets and earliest times at which their observation is complete
observed_sat_names = unique([pass.tle.name for pass in selected_passes])
observed_earliest_JD = Dict()
for pass in selected_passes
    observed_earliest_JD[pass.tle.name] = pass.t0
end

# observer(s)
min_elevation = deg2rad(telescope["min_elevation"])
observer_lat = deg2rad(config["observer"]["latitude"])          # degrees --> radians
observer_lon = deg2rad(config["observer"]["longitude"])         # degrees --> radians
observer_alt = config["observer"]["altitude"]                   # meters
observer_lla = [observer_lat, observer_lon, observer_alt]
r_ECEF_obs = geodetic_to_ecef(observer_lla[1], observer_lla[2], observer_lla[3]) / 1e3

# ------------------------------------------------------------------------ #
# duration of visualization
sec2day = function(dt_sec)
    return dt_sec / 60 / 60 / 24
end
# obs_duration = 1 * 3600    # FIXME overwrite duration of animation
dt_sec = 40.0              # step-size for state query, in seconds
jd_eval = LinRange(jd0_obs, jd0_obs+sec2day(obs_duration), Int(ceil(obs_duration/dt_sec)))    # duration in terms of julian date

# get state history
println("Querying TLE state history...")
vis_frame = :ecef
tles_sph = zeros(length(tles),length(jd_eval),3)
elevations = zeros(length(tles),length(jd_eval))
for (i,tle) in enumerate(tles)
    dts_min = [el for el in  (jd_eval .- tle_epoch(tle)) * 24 * 60]  # days --> minutes
    tles_sph[i,:,:] = TelescopeTasking.integrate_sgp4(
        tle, dts_min, :sph, eop_iau1980, observer_lla)'
    elevations[i,:] = tles_sph[i,:,2]
    # tles_sph[i,:,1] *= 180/π                # azimmuths in radians
    tles_sph[i,:,2] = 90 .- rad2deg.(tles_sph[i,:,2])  # evelations as 90 - rad2deg(elv) for plotting
end

# define whether RSO is visible based on line of sight for each target & time-step
tles_visibility = zeros(Int,length(tles),length(jd_eval))
for isat in 1:size(tles_sph,1)
    for itime in 1:size(tles_sph,2)
        if elevations[isat] > rad2deg(min_elevation)
            tles_visibility[isat,itime] = 1
        end
    end
end

# define whether a target has been observed by a given time
tles_observed_by_now = zeros(Int,length(tles),length(jd_eval))
for (isat,tle) in enumerate(tles)
    for itime in 1:size(tles_sph,2)
        if tle.name in observed_sat_names
            if jd_eval[itime] >= observed_earliest_JD[tle.name]
                tles_observed_by_now[isat,itime] = 1
            end
        end
    end
end

# total number of sat observed over time
N_total = length(tles)
total_observed = zeros(Int, length(jd_eval))
for itime in 1:size(tles_sph,2)
    total_observed[itime] = sum(tles_observed_by_now[:,itime])
end

# define whether a target has been observed by a given time
passes_observed_by_now = zeros(Int,length(selected_passes),length(jd_eval))
for (ipass,pass) in enumerate(selected_passes)
    for itime in 1:size(tles_sph,2)
        if pass.tle.name in observed_sat_names
            if jd_eval[itime] >= observed_earliest_JD[pass.tle.name]
                passes_observed_by_now[ipass,itime] = 1
            end
        end
    end
end

# ------------------------------------------------------------------------ #
# create animation
mutable struct PassAnimator
    framerate::Int
    time_step::Observable{Int}
    orbits_states::Array{Float64, 3}   # N_sat x N_t x 3
    lifted_states

    function PassAnimator(
        framerate::Int,
        orbits_states::Array{Float64,3},
    )
        time_step = Observable{Int}(1)  # initial index
        lifted_states = Vector{Observable{Float64}}[]
        for isat in 1:size(orbits_states,1)
            push!(
                lifted_states,
                [@lift(orbits_states[isat,$time_step,1]),   # azimuth
                 @lift(orbits_states[isat,$time_step,2]),]  # elevation
            )
        end
        new(
            framerate,
            time_step,
            orbits_states,
            lifted_states,
        )
    end
end

function Base.show(io::IO, animator::PassAnimator)
    println(io, "Pass animator struct")
    nsat,ntime,_ = size(animator.orbits_states)
    println(io, "    Number of objects    = $(nsat)")
    println(io, "    Number of time-steps = $(ntime)")
end

# create animator object
save_to_file = true
animation_config = :solution     # :problem or :solution
framerate = 30
animator = OrbitAnimator(framerate, tles_sph)
save_dir = joinpath(@__DIR__, "animations")
timestamps = 1:length(jd_eval)
sleep_display = 0.05

function LinRangeObservable(val,obs,steps)
    dval = (val - obs)/steps
    return [val + dval*i for i in 0:steps]
end

# colors with gradual change across time
colors_timeline = cgrad(:winter, length(timestamps), categorical = true)
colors_greys = cgrad(:greys)    # has inidices 1 ~ 256

with_theme(theme_light()) do
    # initilalize figure
    fontsize = 24
    fig = Figure(size = (500,500))
    time_step_local = animator.time_step    # somehow need to redefine in this scope!
    if save_to_file == true
        title = @lift("$(string(julian2datetime(jd_eval[$time_step_local]))[1:19]) (observed $(total_observed[$time_step_local]) / $(N_total))")
    else
        title = @lift("$(string(julian2datetime(jd_eval[$time_step_local]))[1:19]) time-step = $($time_step_local) / $(size(animator.orbits_states,2))")
    end
    ax1 = PolarAxis(fig[1,1];
        rticks = ([0,30,60], ["90","60","30"]),
        title = title,
        # rlabelsize=fontsize,
        # thetalabelsize=fontsize,
        rgridwidth = 0.3,
        thetagridwidth = 0.3,
        rticklabelsize=fontsize-1,
        thetaticklabelsize=fontsize-1,)
    rlims!(ax1, 0, 60)

    # TelescopeTasking.polar_plot_passes!(ax1, selected_passes;
    #     color=:grey50, linewidth=0.25)

    # plot passes
    for (ipass,pass) in enumerate(selected_passes)
        color_observable = @lift(
            if passes_observed_by_now[ipass,$time_step_local] == 1
                # colors_timeline[$time_step_local]
                colors_timeline[findfirst(isequal(1),passes_observed_by_now[ipass,:])]
            else
                colors_greys[200]
            end
        )
        lw_observable = @lift(
            if passes_observed_by_now[ipass,$time_step_local] == 1
                1.0
            else
                0.1
            end
        )

        _, az_exposure, el_exposure = TelescopeTasking.get_exposure_history(pass)

        lines!(ax1, pass.azimuths, 90 .- rad2deg.(pass.elevations), 
            color = colors_greys[200],
            linewidth = lw_observable)

        lines!(ax1, az_exposure, 90 .- rad2deg.(el_exposure), 
                color = color_observable,
                linewidth = lw_observable)
    end

    # # plot objects
    # for (isat,sat_lifted_states) in enumerate(animator.lifted_states)
    #     color_observable = @lift(
    #         if tles_visibility[isat,$time_step_local] == 1
    #             :red
    #         elseif (animation_config == :solution) && (tles_observed_by_now[isat,$time_step_local] == 1)
    #             :lime
    #         else
    #             :grey
    #         end
    #     )
    #     # line of sight between observer & object
    #     alpha_observable = @lift(
    #         if tles_visibility[isat,$time_step_local] == 1
    #             1.0
    #         else
    #             1.0
    #         end
    #     )

    #     # scatter marker of objects in orbit 
    #     scatter!(ax1, sat_lifted_states[1:2]...,
    #              color = color_observable,
    #              alpha = alpha_observable,
    #              label = nothing, markersize = 5, marker = :diamond)
    # end

    # animate to file
    if save_to_file == true
        println("Saving to file...")
        record(
            fig, 
            joinpath(save_dir, "polar_$(string(animation_config))_$(string(vis_frame))_$(experiment_name)_$(solver_choice).gif");
            framerate = animator.framerate,
        ) do io
            @showprogress for t in timestamps
                animator.time_step[] = t
                recordframe!(io)  # record a new frame
            end
        end

    # live display
    else 
        println("Displaying...")
        display(fig)
        @showprogress for t in timestamps
            animator.time_step[] = t
            sleep(sleep_display)
            # ax1.azimuth[] = mod(ax1.azimuth[] - deg2rad(0.5), 2π)
            # yield() #-> not required with record
        end
    end
end
println("Done!")