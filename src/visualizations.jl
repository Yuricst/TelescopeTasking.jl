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

        lines!(ax, pass.azimuths, 90 .- rad2deg.(pass.elevations),
            color = _color,
            linewidth = linewidth)
    end
end



function plot_time_history!(
    axes::Vector{Axis},
    passes::Vector{VisiblePass};
    color = :black,
    linewidth = 0.4,
    color_by_target::Bool = false,
    cgrad_designator = :winter,
    jd_ref::Union{Float64, Nothing} = nothing
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
        if isnothing(jd_ref)
            times = pass.times
        else
            times = (pass.times .- jd_ref) * 24
        end
        lines!(axes[1], times, rad2deg.(pass.azimuths), color = _color, linewidth = linewidth)
        lines!(axes[2], times, rad2deg.(pass.elevations), color = _color, linewidth = linewidth)
    end
end