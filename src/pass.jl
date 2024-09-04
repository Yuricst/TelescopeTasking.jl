"""Struct to handle vsible passes"""


"""
Struct holds visible a single visible pass
"""
struct VisiblePass
    tle::SatelliteToolboxTle.TLE
    t0::Number
    tm::Number
    tf::Number
    azel0::Vector{Float64}
    azelm::Vector{Float64}
    azelf::Vector{Float64}

    times::Vector{Float64}
    azimuths::Vector{Float64}
    elevations::Vector{Float64}

    function VisiblePass(tle, pass_times, pass_azimuths, pass_elevations)
        new(
            tle,
            pass_times[1],
            pass_times[div(end,2)],
            pass_times[end],
            [pass_azimuths[1], pass_elevations[1]],
            [pass_azimuths[div(end,2)], pass_elevations[div(end,2)]],
            [pass_azimuths[end], pass_elevations[end]],
            pass_times,
            pass_azimuths,
            pass_elevations
        )
    end
end


"""
Overload method for showing
"""
function Base.show(io::IO, pass::VisiblePass)
    println(io, "Visible Pass of $(pass.tle.name) (number $(pass.tle.satellite_number))")
    println(io, "    Start: $(pass.t0) JD")
    println(io, "    End: $(pass.tf) JD")
    println(io, "    Max elevation: $(rad2deg(pass.azelm[2])) deg")
    println(io, "    Duration (min): $((pass.tf - pass.t0) * 24 * 60) min")
end


"""
Get vector of passes from azimuth and elevation time history
"""
function azel_history_to_passes(
    tle::SatelliteToolboxTle.TLE,
    times::Vector{Float64},
    azimuths::Vector{Float64},
    elevations::Vector{Float64},
    min_elevation::Number,
    min_duration::Number,
)
    passes = VisiblePass[]
    new_pass = true         # initialize boolean for checking if we are in a new pass
    in_pass = false         # initialize boolean for checking if we are in a pass

    pass_times, pass_azimuths, pass_elevations = [], [], []  # initialize pass vectors

    # iterate through time, azimuth, and elevation
    for (t,az,el) in zip(times, azimuths, elevations)
        if (new_pass == true) && (el > min_elevation)       # start a new pass
            pass_times = [t,]               # re-initialize
            pass_azimuths = [az,]           # re-initialize
            pass_elevations = [el,]         # re-initialize
            new_pass = false
            in_pass  = true

        elseif (new_pass == false) && (el > min_elevation)  # continue the pass_azimuths
            push!(pass_times, t)
            push!(pass_azimuths, az)
            push!(pass_elevations, el)
            @assert in_pass == true

        elseif (new_pass == false) && (el <= min_elevation) && (in_pass == true) # end the pass
            if (pass_times[end] - pass_times[1]) > min_duration
                push!(passes, VisiblePass(tle, pass_times, pass_azimuths, pass_elevations))
            end
            new_pass = true             # toggle for next pass detection
            in_pass = false             # toggle for next pass detection
        end
    end
    return passes
end