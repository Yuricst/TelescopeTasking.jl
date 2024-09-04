"""Struct to handle vsible passes"""


"""
Struct holds visible a single visible pass
"""
struct VisiblePass
    tle::SatelliteToolboxTle.TLE
    t0::Number
    tm::Number
    tf::Number
    t0_exposure::Number
    tf_exposure::Number

    azel0::Vector{Float64}
    azelm::Vector{Float64}
    azelf::Vector{Float64}

    azel0_exposure::Vector{Float64}
    azelf_exposure::Vector{Float64}

    times::Vector{Float64}
    azimuths::Vector{Float64}
    elevations::Vector{Float64}

    function VisiblePass(tle, pass_times, pass_azimuths, pass_elevations, exposure_duration::Number)
        # compute exposure times, cented at the middle of the pass
        tm = pass_times[div(end,2)]
        t0_exposure = tm - exposure_duration/2
        tf_exposure = tm + exposure_duration/2

        # spline interpolate azimuth and elevation
        az_itp_cubic = linear_interpolation(pass_times, pass_azimuths)
        el_itp_cubic = linear_interpolation(pass_times, pass_elevations)
        azel0_exposure = [az_itp_cubic(t0_exposure), el_itp_cubic(t0_exposure)]
        azelf_exposure = [az_itp_cubic(tf_exposure), el_itp_cubic(tf_exposure)]

        
        # create instance
        new(
            tle,
            pass_times[1],
            pass_times[div(end,2)],
            pass_times[end],
            t0_exposure,
            tf_exposure,
            [pass_azimuths[1], pass_elevations[1]],
            [pass_azimuths[div(end,2)], pass_elevations[div(end,2)]],
            [pass_azimuths[end], pass_elevations[end]],
            azel0_exposure,
            azelf_exposure,
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
    println(io, "    Exposure start:        $(pass.t0_exposure) JD")
    println(io, "    Exposure end:          $(pass.tf_exposure) JD")
    println(io, "    Max elevation:         $(rad2deg(pass.azelm[2])) deg")
    println(io, "    Pass duration (min):   $((pass.tf - pass.t0) * 24 * 60) min")
    println(io, "    Number of data points: $(length(pass.times))")
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
    exposure_duration::Number,
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
                push!(passes, VisiblePass(tle, pass_times, pass_azimuths, pass_elevations, exposure_duration))
            end
            new_pass = true             # toggle for next pass detection
            in_pass = false             # toggle for next pass detection
        end
    end
    return passes
end


"""
Sort list of visible passes
"""
function sort(passes::Vector{VisiblePass})
    # get all exposure start times
    ts_exposure_start = [pass.t0_exposure for pass in passes]
    # sort passes by exposure start time
    passes_sorted = passes[sortperm(ts_exposure_start)]
    return passes_sorted
end


"""
Filter list of visible passes
"""
function filter(passes::Vector{VisiblePass}, E::Int)
    pass_to_designators = [pass.tle.international_designator for pass in passes]
    designators = unique(pass_to_designators)
end