"""Functions to compute quality metrics"""


"""
Compute temporal efficiency of tasking solution
"""
function get_temporal_efficiency(
    observation_duration::Real,
    exposure_duration::Real,
    num_passes_used::Real,
)
    return exposure_duration / observation_duration * num_passes_used
end


"""
Compute non-exposure distance L of tasking solution
"""
function get_nonexposure_distance(
    passes::Vector,
    Y::Union{BitVector, Vector{Bool}},
)
    return 0
end