"""Functions for path on sphere"""


"""
Perform spherical linear interpolation (slerp) between two vectors
"""
function slerp(v0, v1, t)
    omega = acos(dot(v0, v1) / (norm(v0) * norm(v1)))
    so = sin(omega)
    v0 = sin((1 - t) * omega) / so * v0
    v1 = sin(t * omega) / so * v1
    return v0 + v1
end


"""
Compute angle history along shortest path on a sphere
"""
function sphere_shortest_path(az1, el1, az2, el2, num_points)
    # Convert spherical to Cartesian
    p1 = sph2cart([az1, el1, 1])
    p2 = sph2cart([az2, el2, 1])
    
    # Initialize lists to store angles
    azimuths = Float64[]
    elevations = Float64[]
    
    # Interpolate points
    for t in 0:1/(num_points-1):1
        point = slerp(p1, p2, t)
        azelr = cart2sph(point)
        push!(azimuths, azelr[1])
        push!(elevations, azelr[2])
    end
    return azimuths, elevations
end