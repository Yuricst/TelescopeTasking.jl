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

    for (i,pass) in enumerate(passes)
        if color_by_target == true
            index = findfirst(x -> x == pass.tle.international_designator, designators)
            _color = colors[index]
        else
            if typeof(color) == Symbol
                _color = color
            else
                _color = color[i]
            end
        end
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

    for (i,pass) in enumerate(passes)
        if color_by_target == true
            index = findfirst(x -> x == pass.tle.international_designator, designators)
            _color = colors[index]
        else
            if typeof(color) == Symbol
                _color = color
            else
                _color = color[i]
            end
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
    color::Symbol = :black,
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


"""
Plot slewing maneuver between two passes
"""
function polar_plot_interpass_slew!(
    ax::PolarAxis,
    passes::Vector{VisiblePass},
    colors;
    linewidth::Real = 1.0,
    linestyle = :solid,
    steps::Int = 100,
)
    @assert length(colors) >= length(passes) - 1
    for i in 1:length(passes)-1
        _, az_1, el_1 = get_exposure_history(passes[i])
        _, az_2, el_2 = get_exposure_history(passes[i+1])
        az_slew, el_slew = sphere_shortest_path(az_1[end], el_1[end], az_2[1], el_2[1], steps)
        lines!(ax, az_slew, 90 .- rad2deg.(el_slew),
            color = colors[i],
            linewidth = linewidth,
            linestyle = linestyle)
    end
end



"""
Plot wireframe of a sphere

# Arguments
- `ax::Union{Axis3,LScene}`: Axis3 or LScene object
- `radius::Real`: radius of the sphere
- `center::Vector`: center of the sphere
- `nsph::Int=20`: number of points along latitude and longitude
- `color=:black`: color of the wireframe
- `linewidth=1.0`: linewidth of the wireframe
"""
function plot_sphere_wireframe!(ax::Union{Axis3,LScene}, radius::Real, center::Vector, nsph::Int=20;
    color=:black, linewidth=1.0, label=nothing)
    # Generate spherical coordinates
    θ = range(0, stop=2π, length=nsph)
    ϕ = range(0, stop=π, length=Int(ceil(nsph/2)))

    # Generate points on the sphere
    xsphere = [center[1] + radius * cos(θ[i]) * sin(ϕ[j]) for j in 1:Int(ceil(nsph/2)), i in 1:nsph]
    ysphere = [center[2] + radius * sin(θ[i]) * sin(ϕ[j]) for j in 1:Int(ceil(nsph/2)), i in 1:nsph]
    zsphere = [center[3] + radius * cos(ϕ[j]) for j in 1:Int(ceil(nsph/2)), i in 1:nsph]
    wireframe!(ax, xsphere, ysphere, zsphere, 
               label=label, color=color, linewidth=linewidth)
    return
end