"""Functions to process input TLEs"""


function remove!(a, item)
    deleteat!(a, findall(x->x==item, a))
end


function tle2sma(tle::SatelliteToolboxTle.TLE)
    mean_motion_rad_per_sec = tle.mean_motion * 2π / 86400
    return (398600.4418 / mean_motion_rad_per_sec^2)^(1/3)
end


function tle2period(tle::SatelliteToolboxTle.TLE)
    return 86400/tle.mean_motion
end



function filter(
    tles::Vector{SatelliteToolboxTle.TLE};
    names_include::Union{Vector{String},Nothing} = nothing,
    sma_max = 1e16,
)
    """Filter TLEs by eccentricity"""
    tles_filtered = SatelliteToolboxTle.TLE[]
    for tle in tles
        # check on SMA
        _sma = tle2sma(tle)

        # check based on name
        if isnothing(names_include)
            name_valid = true
        else
            name_valid = false
            for substring in names_include
                if occursin(substring, tle.name)
                    name_valid = true
                    break
                end
            end
        end
        
        # if not valid, remove from list
        if (_sma > sma_max) || (name_valid == false)
            remove!(tles, tle)
        else
            push!(tles_filtered, tle)
        end
    end
    return tles_filtered
end


"""
    integrate_sgp4(tle, dts_min, return_type, eop_iau1980, observer_lla)

Integrate TLE with SGP4 model.

# Arguments
- `tle::SatelliteToolboxTle.TLE`: TLE object
- `dts_min::Union{Vector,StepRangeLen}`: time grid for evaluating state in minutes
- `return_type::Symbol`: return type, one of `:eci`, `:ecef`, `:enu`, `:sph`
- `eop_iau1980::EopIau1980`: Earth orientation parameters
- `observer_lla::Vector`: observer location in geodetic coordinates
"""
function integrate_sgp4(
    tle::SatelliteToolboxTle.TLE,
    dts_min::Union{Vector,StepRangeLen},
    return_type::Symbol,
    eop_iau1980::EopIau1980,
    observer_lla::Vector,
)
    # initialize integrator and storage
    sgp4d = sgp4_init(tle)
    jd0_tle = tle_epoch(tle)            # initial epoch of TLE in Julian day

    if return_type == :eci
        storage = zeros(6, length(dts_min))
    elseif return_type == :ecef
        storage = zeros(3, length(dts_min))
    elseif return_type == :enu
        storage = zeros(3, length(dts_min))
    elseif return_type == :sph
        storage = zeros(3, length(dts_min))
    else
        error("Invalid return type $return_type")
    end

    # observer location in ENU
    T_NED2ENU = dcm_ned_to_enu()

    # geodetic_to_ecef() expects h in m and returns in m; we convert to in km
    r_ECEF_obs = geodetic_to_ecef(observer_lla[1], observer_lla[2], observer_lla[3]) / 1e3
    r_ENU_obs = T_NED2ENU * ecef_to_ned(r_ECEF_obs, observer_lla[1], observer_lla[2], observer_lla[3]; translate = false)

    for (idx,dt_min) in enumerate(dts_min)
        # propagate state in TEME frame
        r_TEME, v_TEME = sgp4!(sgp4d, dt_min)
        
        # if only ECI coordinates are needed, we don't need the remainder of the for-loop
        if return_type == :eci
            storage[:,idx] = [r_TEME; v_TEME]
            continue
        end

        # current julian date
        jd = jd0_tle + dt_min/60/24

        # get transformation matrix from ECI to ECEF and apply conversion
        T_TEME2ITRF = r_eci_to_ecef(TEME(), ITRF(), jd, eop_iau1980)
        r_ITRF = T_TEME2ITRF * r_TEME

        # convert to azimuth and elevation angles
        r_NED = ecef_to_ned(r_ITRF, observer_lla[1], observer_lla[2], observer_lla[3]; translate = false)
        r_ENU = T_NED2ENU * r_NED - r_ENU_obs
        sph_ENU = cart2sph(r_ENU)    # [az,el,radius]

        # store for non-ECI case
        if return_type == :ecef
            storage[:,idx] = r_ITRF
        elseif return_type == :enu
            storage[:,idx] = r_ENU
        elseif return_type == :sph
            storage[:,idx] = sph_ENU
        end
    end
    return storage
end


"""
    tles_to_passes(
        tles::Vector{SatelliteToolboxTle.TLE},
        eop_iau1980::EopIau1980,
        jd0_obs,
        min_elevation,
        min_obs_duration,
        exposure_duration,
        observer_lla::Vector;
        dt_sec::Number = 10.0,
    )

Construct vector of visible passes from vector of TLEs.
All time arguments except jd0 are in seconds. 
"""
function tles_to_passes(
    tles::Vector{SatelliteToolboxTle.TLE},
    eop_iau1980::EopIau1980,
    jd0_obs,
    obs_duration,
    min_elevation,
    min_obs_duration,
    exposure_duration,
    observer_lla::Vector;
    dt_sec::Number = 10.0,
    num_exposure::Int = 1,
    filter_by_num_exposure::Bool = false,
)
    @assert deg2rad(min_elevation) >= 0 "Minimum elevation must be positive"
    @assert min_obs_duration >= exposure_duration "Minimum observation duration must be greater than exposure duration"

    sph_ENU_list = []

    passes = VisiblePass[]                              # initialize vector of passes   
    dt_min = dt_sec / 60                                # convert to minutes

    for (idx,tle) in enumerate(tles)
        # initialize time grid
        jd0_tle = tle_epoch(tle)            # initial epoch of TLE in Julian day
        
        tprop_extra = 30                                                                # extra time to make sure observation is within time window
        t0_tle = (jd0_obs - jd0_tle) * 24*60                                            # time between TLE time and observation start time, in minutes
        dts_min_tle = t0_tle-tprop_extra:dt_min:(t0_tle + obs_duration/60)+tprop_extra  # propagation time grid in minutes
        times_jd = Array(jd0_tle .+ dts_min_tle/60/24)                                  # time grid in Julian day
        jd_obs_window = [jd0_obs, jd0_obs + obs_duration/86400]                         # observation window in Julian day

        # integrate over time grid and get passes
        sph_ENU = integrate_sgp4(tle, dts_min_tle, :sph, eop_iau1980, observer_lla)
        _passes = azel_history_to_passes(
            tle, times_jd, sph_ENU[1,:], sph_ENU[2,:], 
            min_elevation, min_obs_duration, exposure_duration, jd_obs_window)
        push!(passes, _passes...)
        push!(sph_ENU_list, sph_ENU)
    end

    # filter ones that do not result in enough exopsures
    if (num_exposure > 1) && (filter_by_num_exposure == true)
        passes = filter(passes, num_exposure)
    end

    # sort passes by exposure start time
    passes .= sort(passes)
    return passes, sph_ENU_list
end