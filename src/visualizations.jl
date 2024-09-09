"""Visualization of results"""


function periodic_lines!(ax, xs, ys, diff_max; color=:red, linewidth=0.5)
    i_from = 1
    for i in 1:length(xs)-1
        diff = abs(xs[i] - xs[i+1])
        if diff > diff_max          # stop current line here
            lines!(ax, xs[i_from:i], ys[i_from:i], color=color, linewidth=linewidth)
            i_from = i + 1
        end
    end
end



"""
Create polar plot of passes
"""
function polar_plot_passes!(
    ax::Axis,
    passes::Vector{VisiblePass};
    color = :black,
    linewidth = 0.4,
    color_by_target::Bool = false,
    cgrad_designator = :winter,
)
    if color_by_target == true
        designators = unique([pass.tle.international_designator for pass in passes])
        colors = cgrad(cgrad_designator, max(2,length(designators)), categorical = true)
    end

    for pass in passes
        if color_by_target == true
            index = findfirst(x -> x == pass.tle.international_designator, designators)
            _color = colors[index]
        else
            _color = color
        end

        # periodic_lines!(ax, rad2deg.(pass.azimuths), rad2deg.(pass.elevations), 100,
        #     color = _color,
        #     linewidth = linewidth)
        lines!(ax, rad2deg.(pass.azimuths), rad2deg.(pass.elevations), 
            color = _color,
            linewidth = linewidth)
    end
end



"""
Create polar plot of passes
"""
function polar_plot_passes!(
    ax::PolarAxis,
    passes::Vector{VisiblePass};
    color = :black,
    linewidth = 0.4,
    color_by_target::Bool = false,
    cgrad_designator = :winter,
    exposure_only::Bool = false,
)
    if color_by_target == true
        designators = unique([pass.tle.international_designator for pass in passes])
        colors = cgrad(cgrad_designator, max(2,length(designators)), categorical = true)
    end

    for pass in passes
        if color_by_target == true
            index = findfirst(x -> x == pass.tle.international_designator, designators)
            _color = colors[index]
        else
            _color = color
        end
        
        if exposure_only == false
            lines!(ax, pass.azimuths, 90 .- rad2deg.(pass.elevations),
                color = _color,
                linewidth = linewidth)
        else
            _, az_exposure, el_exposure = get_exposure_history(pass)
            lines!(ax, az_exposure, 90 .- rad2deg.(el_exposure),
                color = _color,
                #marker= :circle,
                #markersize = 1,)
                linewidth = linewidth)
        end
    end
end


"""
Plot time-history of azimuth and elevation
"""
function plot_time_history!(
    axes::Vector{Axis},
    passes::Vector{VisiblePass};
    color = :black,
    linewidth = 0.4,
    color_by_target::Bool = false,
    cgrad_designator = :winter,
    exposure_only::Bool = false,
    jd_ref::Union{Float64, Nothing} = nothing,
    time_multiplier::Real = 24,
)
    if color_by_target == true
        designators = unique([pass.tle.international_designator for pass in passes])
        colors = cgrad(cgrad_designator, max(2,length(designators)), categorical = true)
    end

    for pass in passes
        if color_by_target == true
            index = findfirst(x -> x == pass.tle.international_designator, designators)
            _color = colors[index]
        else
            _color = color
        end
        
        if exposure_only == false
            
            if isnothing(jd_ref)
                times = pass.times
            else
                times = (pass.times .- jd_ref) * time_multiplier
            end
            lines!(axes[1], times, rad2deg.(pass.azimuths), color = _color, linewidth = linewidth)
            lines!(axes[2], times, rad2deg.(pass.elevations), color = _color, linewidth = linewidth)
        else
            times_exposure, az_exposure, el_exposure = get_exposure_history(pass)
            if isnothing(jd_ref)
                times_exposure = times_exposure
            else
                times_exposure = (times_exposure .- jd_ref) * time_multiplier
            end
            lines!(axes[1], times_exposure, rad2deg.(az_exposure), color = _color, linewidth = linewidth)
            lines!(axes[2], times_exposure, rad2deg.(el_exposure), color = _color, linewidth = linewidth)
        end
    end
end


"""
Plot slewing maneuver between two passes
"""
function polar_plot_interpass_slew!(
    ax::PolarAxis,
    passes::Vector{VisiblePass};
    color = :black,
    linewidth::Real = 1.0,
    linestyle = :solid,
    steps::Int = 100,
)
    for i in 1:length(passes)-1
        _, az_1, el_1 = get_exposure_history(passes[i])
        _, az_2, el_2 = get_exposure_history(passes[i+1])
        az_slew, el_slew = sphere_shortest_path(az_1[end], el_1[end], az_2[1], el_2[1], steps)
        lines!(ax, az_slew, 90 .- rad2deg.(el_slew),
            color = color,
            linewidth = linewidth,
            linestyle = linestyle)
    end
end