    """Load TLE files"""

using GLMakie
using GeometryBasics
using LinearAlgebra
using ProgressMeter: @showprogress
using SatelliteToolboxTle
using SatelliteToolboxSgp4

# load TLE files
path_to_tles = joinpath(@__DIR__, "..", "data", "tles", "iridium-33-debris.txt")
tles_str = read(path_to_tles, String)

# convert to TLE objects
tles = read_tles(tles_str)

# propagate one of them
sgp4d = sgp4_init(tles[1])
dt_min = 1.0
r_teme, v_teme = sgp4!(sgp4d, dt_min)

# propagate over 12 hours
dts_min = 0:dt_min:12*60
rs_teme, vs_teme = zeros(3, length(dts_min)), zeros(3, length(dts_min))
for (idx,dt_min) in enumerate(dts_min)
    rs_teme[:,idx], vs_teme[:,idx] = sgp4!(sgp4d, dt_min)
end

# plot in 3D
fig = Figure(size=(1000, 500))
ax3 = Axis3(fig[1,1], aspect = (1, 1, 1))
ax_history = Axis(fig[1,2]; xlabel="Time, min", ylabel="Radius, km")

# plot Earth
# mesh!(ax, Sphere(GeometryBasics.Point3f0(0), 6378.0), alpha=0.3)

# propagate TLEs over 12 hours and plot
println("Propagating TLEs...")
@showprogress for _tle in tles
    _sgp4d = sgp4_init(_tle)
    _rs_teme, _vs_teme = zeros(3, length(dts_min)), zeros(3, length(dts_min))
    for (_idx,_dt_min) in enumerate(dts_min)
        _rs_teme[:,_idx], _vs_teme[:,_idx] = sgp4!(_sgp4d, _dt_min)
    end
    lines!(ax3, _rs_teme[1,:], _vs_teme[2,:], rs_teme[3,:], color=:red, linewidth=0.8)

    rnorms = [norm(_r) for _r in eachcol(_rs_teme)]
    lines!(ax_history, dts_min, rnorms, color=:red, linewidth=0.8)
end

display(fig)