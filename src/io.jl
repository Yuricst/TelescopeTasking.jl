"""Input/output functions"""


function tles_to_file(tles::Vector{SatelliteToolboxTle.TLE}, filename)
    open(filename, "w") do file
        write(file, [convert(String, tle) * "\n" for tle in tles]...)
    end
    println("Saved TLEs to $(filename)!")
    return
end


"""Convert struct to dictionary"""
typedict(x::T) where {T} = Dict(fn=>getfield(x, fn) for fn ∈ fieldnames(T))


"""Convert pass to dictionary"""
function pass_to_dict(pass::VisiblePass)
    pass_dict = typedict(pass)
    delete!(pass_dict, :tle)
    pass_dict[:tle_string] = convert(String, pass.tle)
    return pass_dict
end


"""
Save solution to dictionary
"""
function solution_to_dict(
    passes::Vector,
    jd0_obs,
    obs_duration,
    min_elevation,
    min_obs_duration,
    exposure_duration,
    observer_lla::Vector,
    X::Union{BitVector, Vector{Bool}},
    Y::Union{BitVector, Vector{Bool}},
)
    # create list of visible passes dictionaries
    passes_dict = [pass_to_dict(pass) for pass in passes]
    return Dict(
        "passes_dict" => passes_dict,
        "jd0_obs" => jd0_obs,
        "obs_duration" => obs_duration,
        "min_elevation" => min_elevation,
        "min_obs_duration" => min_obs_duration,
        "exposure_duration" => exposure_duration,
        "observer_lla" => observer_lla,
        "X" => X,
        "Y" => Y,
    )
end