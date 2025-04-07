module TG

using SatelliteToolboxBase
using SatelliteToolboxPropagators
using GLMakie
using Setfield

export plot_orbit

function plot_orbit(orb::KeplerianElements)
    N = 100
    θ = LinRange(0, 2π, N)

    orbit = ((@set orb.f = theta) for theta in θ) .|> kepler_to_rv .|> first |> stack
    current_pos, current_velocity = kepler_to_rv(orb)
    velocity_arrow_base = current_pos
    velocity_arrow_tip = velocity_arrow_base + current_velocity / √sum(current_velocity .^ 2) * 0.4 * √sum(current_pos .^ 2)
    velocity_arrow_data = [velocity_arrow_base velocity_arrow_tip] |> Matrix
    fig = Figure()
    ax3d = Axis3(fig[1, 1])
    lines!(ax3d, orbit[1, :], orbit[2, :], orbit[3, :], color=:red)
    wireframe!(ax3d, Sphere(Point3(0.0), EARTH_EQUATORIAL_RADIUS), color=:cyan, alpha=0.3)
    scatter!(ax3d, current_pos, markersize=20)
    lines!(ax3d, velocity_arrow_data[1, :], velocity_arrow_data[2, :], velocity_arrow_data[3, :])
    fig
end

end # module TG
