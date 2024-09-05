"""Visualization of results"""


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