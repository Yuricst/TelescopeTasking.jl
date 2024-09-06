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
end


"""
Create instance of VisiblePass from TLE and pass data. This constructor is called inside
the function `tles_to_passes()`. 
"""
function VisiblePass(tle, pass_times, pass_azimuths, pass_elevations, exposure_duration::Number)
    # compute exposure times, cented at the middle of the pass
    tm = (pass_times[1] + pass_times[end]) / 2 
    t0_exposure = tm - exposure_duration/86400 / 2
    tf_exposure = tm + exposure_duration/86400 / 2

    # spline interpolate azimuth and elevation
    az_itp_cubic = linear_interpolation(pass_times, pass_azimuths)
    el_itp_cubic = linear_interpolation(pass_times, pass_elevations)
    @assert pass_times[1] <= t0_exposure <= pass_times[end]
    @assert pass_times[1] <= tf_exposure <= pass_times[end]
    azel0_exposure = [az_itp_cubic(t0_exposure), el_itp_cubic(t0_exposure)]
    azelf_exposure = [az_itp_cubic(tf_exposure), el_itp_cubic(tf_exposure)]
    
    # create instance of VisiblePass
    VisiblePass(
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


"""
Construct instance of VisiblePass from dictionary
"""
function VisiblePass(pass::Dict)
    VisiblePass(
        read_tle(pass["tle_string"]),
        pass["t0"],
        pass["tm"],
        pass["tf"],
        pass["t0_exposure"],
        pass["tf_exposure"],
        pass["azel0"],
        pass["azelm"],
        pass["azelf"],
        pass["azel0_exposure"],
        pass["azelf_exposure"],
        pass["times"],
        pass["azimuths"],
        pass["elevations"]
    )
end


"""
Overload method for showing
"""
function Base.show(io::IO, pass::VisiblePass)
    println(io, "Visible Pass of $(pass.tle.name) (designator $(pass.tle.international_designator))")
    println(io, "    Exposure start:        $(pass.t0_exposure) JD")
    println(io, "    Exposure end:          $(pass.tf_exposure) JD")
    println(io, "    Max elevation:         $(rad2deg(pass.azelm[2])) deg")
    println(io, "    Pass duration (min):   $((pass.tf - pass.t0) * 24 * 60) min")
    println(io, "    Number of data points: $(length(pass.times))")
end


"""
    get_exposure_history(pass::VisiblePass)
    
Get time stamps, azimuths, and elevations during exposure
"""
function get_exposure_history(pass::VisiblePass)
    # select times between exposure
    times_after_t0 = pass.times[pass.times .>= pass.t0_exposure]
    times_exposure = times_after_t0[times_after_t0 .<= pass.tf_exposure]

    az_after_t0 = pass.azimuths[pass.times .>= pass.t0_exposure]
    az_exposure = az_after_t0[times_after_t0 .<= pass.tf_exposure]
    el_after_t0 = pass.elevations[pass.times .>= pass.t0_exposure]
    el_exposure = el_after_t0[times_after_t0 .<= pass.tf_exposure]
    return times_exposure, az_exposure, el_exposure
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
    min_obs_duration::Number,
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
            if (pass_times[end] - pass_times[1]) > min_obs_duration / 86400
                push!(passes, VisiblePass(tle, pass_times, pass_azimuths, pass_elevations, exposure_duration))
                #println("Got a pass: $((pass_times[end] - pass_times[1])*86400) sec")
            # else
            #     println("Too short: $((pass_times[end] - pass_times[1])*86400) sec, need $min_obs_duration sec")
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
Remove passes with less than `num_exposure` exposures from list of passes
"""
function filter(passes::Vector{VisiblePass}, num_exposure::Int)
    pass_to_designators = [pass.tle.international_designator for pass in passes]
    designators = unique(pass_to_designators)
    passes_filtered = VisiblePass[]
    for designator in designators
        passes_designator = [pass for pass in passes if pass.tle.international_designator == designator]
        if length(passes_designator) >= num_exposure
            append!(passes_filtered, passes_designator)
        end
    end
    return passes_filtered
end