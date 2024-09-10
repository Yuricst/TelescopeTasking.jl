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
Compute subsequent non-exposure distances L of tasking solution for a single telescope.

```math
L = sum_{i=1}^{n-1} sum_{j=i}^{n} acos(dot(r_unit_i, r_unit_j)) * Y_i * Y_j
```
"""
function get_nonexposure_distance(
    passes::Vector,
    Y,
)
    # chronological
    used_passes = [pass for (pass, y) in zip(passes, Y) if y > 0.5]
    Ls = Real[]
    for i in 1:length(used_passes)-1
        r_unit_i = sph2cart(vcat(used_passes[i].azelf_exposure, [1]))
        r_unit_j = sph2cart(vcat(used_passes[i+1].azel0_exposure, [1]))
        @assert norm(r_unit_j) ≈ 1.0
        @assert norm(r_unit_i) ≈ 1.0
        L = acos(dot(r_unit_i, r_unit_j))
        push!(Ls, L)
    end
    return sum(Ls), Ls
end


"""
Compute ratio of time where telescope was idle
"""
function get_idletimeratio(
    passes,
    Y,
    obs_duration,
    exposure_duration,
    slew_rate
)
    L, _ = get_nonexposure_distance(passes, Y)
    return 1 - 1/obs_duration*(L/slew_rate + exposure_duration*sum(Y))
end