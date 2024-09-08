"""Handling sun location"""


"""
Sun direction in ENU frame
"""
function sun_direction_enu(
    jd::Real,
    observer_lla::Vector,
    eop_iau1980::EopIau1980,
)
    # observer location in ENU
    T_NED2ENU = dcm_ned_to_enu()
    
    # query sun position and transform to ENU frame
    i_sun = sun_position_mod(jd)
    T_TEME2ITRF = r_eci_to_ecef(TEME(), ITRF(), jd, eop_iau1980)
    i_NED = ecef_to_ned(T_TEME2ITRF * i_sun/norm(i_sun), 
        observer_lla[1], observer_lla[2], observer_lla[3]; translate = false)
    i_ENU = T_NED2ENU * i_NED
    return i_ENU
end



"""
Compute earliest night time at a given observer location
"""
function earliest_night(
    jd0::Real,
    observer_lla::Vector,
    eop_iau1980::EopIau1980,;
    elevation_threshold::Real = -deg2rad(12),
    steps::Int = 2000,
)
    # epochs at which we query sun position
    jds = LinRange(jd0, jd0 + 2.5, steps)

    # sun position over time
    sph_sun = zeros(3, length(jds))
    for (i,jd) in enumerate(jds)
        isun = sun_direction_enu(jd, observer_lla, eop_iau1980)
        sph_sun[:,i] = TelescopeTasking.cart2sph(isun)
    end

    # detect earliest beginning of night time
    idx_dusk = 1
    for idx in 1:length(sph_sun[2,:])-1
        if (sph_sun[2,idx] < elevation_threshold) && (sph_sun[2,idx] > sph_sun[2,idx+1])
            idx_dusk = idx
            break
        end
    end
    jd_night_start = jds[idx_dusk]

    # remaining time
    jd_remain = jds[idx_dusk:end]
    sph_sun_remain = sph_sun[:,idx_dusk:end]

    # detect earliest end of night time
    idx_dawn = findfirst(sph_sun_remain[2,:] .> elevation_threshold)
    jd_night_end = jd_remain[idx_dawn]
    return [jd_night_start,jd_night_end]
end