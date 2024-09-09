"""Input/output functions"""


function tles_to_file(tles::Vector{SatelliteToolboxTle.TLE}, filename)
    open(filename, "w") do file
        write(file, [convert(String, tle) * "\n" for tle in tles]...)
    end
    println("Saved TLEs to $(filename)!")
    return
end


"""Convert struct to dictionary"""
#typedict(x::T) where {T} = Dict(fn=>getfield(x, fn) for fn ∈ fieldnames(T))
function typedict(x::T) where {T}
    #dict_symbol = Dict(fn=>getfield(x, fn) for fn ∈ fieldnames(T))
    dict_string = Dict(string(fn)=>getfield(x, fn) for fn ∈ fieldnames(T))
    return dict_string
end


"""Convert pass to dictionary"""
function pass_to_dict(pass::VisiblePass)
    pass_dict = typedict(pass)
    delete!(pass_dict, :tle)
    pass_dict["tle_string"] = convert(String, pass.tle)
    return pass_dict
end


"""
Save solution to dictionary
"""
function STTP_solution_to_dict(
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
    # extract problem dimensions
    m = length(X)
    n = length(Y)

    # compute quality quality_metrics
    eta = get_temporal_efficiency(obs_duration, exposure_duration, sum(Y))
    L, Ls = get_nonexposure_distance(passes, Y)

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
        "eta" => eta,
        "Ls" => Ls,
        "L" => L,
        "m" => m,
        "n" => n,
    )
end


"""
Save solution to dictionary
"""
function MTTP_solution_to_dict(
    problem::MultiTelescopeTaskingProblem,
    passes_per_telescope::Vector,
    jd0_obs,
    obs_duration_per_telescope,
    min_elevation,
    min_obs_duration,
    exposure_duration,
    observer_lla_per_telescope::Vector,
    X::Union{BitVector, Vector{Bool}, Vector{Int}},
    Y_per_telescope::Vector{Vector},
)
    # compute quality quality_metrics
    eta_per_telescope = Real[]
    L_per_telescope = Real[]
    Ls_per_telescope = Vector[]
    for (Y,obs_duration,passes) in zip(Y_per_telescope, obs_duration_per_telescope, passes_per_telescope)
        push!(eta_per_telescope, get_temporal_efficiency(obs_duration, exposure_duration, sum(Y)))
        L, Ls = get_nonexposure_distance(passes, Y)
        push!(L_per_telescope, L)
        push!(Ls_per_telescope, Ls)
    end

    # create list of visible passes dictionaries
    passes_dict = [[pass_to_dict(pass) for pass in passes] for passes in passes_per_telescope]
    return Dict(
        "passes_per_telescope" => passes_dict,
        "jd0_obs" => jd0_obs,
        "obs_duration_per_telescope" => obs_duration_per_telescope,
        "min_elevation" => min_elevation,
        "min_obs_duration" => min_obs_duration,
        "exposure_duration" => exposure_duration,
        "observer_lla_per_telescope" => observer_lla_per_telescope,
        "X" => X,
        "Y_per_telescope" => Y_per_telescope,
        "eta_per_telescope" => eta_per_telescope,
        "Ls_per_telescope" => Ls_per_telescope,
        "L_per_telescope" => L_per_telescope,
        "m" => problem.m,
        "n_per_telescope" => problem.n_per_telescope,
        "map_y2Y" => problem.map_y2Y,
    )
end