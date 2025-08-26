using GLMakie
using SatelliteToolboxBase

function plot_orbit(orbs::KeplerianElements...)
    N = 100
    θ = LinRange(0, 2π, N)

    fig = Figure()
    ax3d = Axis3(fig[1, 1])
    
    for (orb, c) in zip(orbs, Makie.wong_colors())
        orbit = ((@set orb.f = theta) for theta in θ) .|> kepler_to_rv .|> first |> stack
        current_pos, current_velocity = kepler_to_rv(orb)
        velocity_arrow_base = current_pos
        velocity_arrow_tip = velocity_arrow_base + current_velocity / √sum(current_velocity .^ 2) * 0.4 * √sum(current_pos .^ 2)
        velocity_arrow_data = [velocity_arrow_base velocity_arrow_tip] |> Matrix
        lines!(ax3d, orbit[1, :], orbit[2, :], orbit[3, :], color=c)
        scatter!(ax3d, current_pos, markersize=20, color=c)
        lines!(ax3d, velocity_arrow_data[1, :], velocity_arrow_data[2, :], velocity_arrow_data[3, :], color=c)
    end
    
    wireframe!(ax3d, Sphere(Point3(0.0), EARTH_EQUATORIAL_RADIUS), color=:cyan, alpha=0.3)
    fig, ax3d
end

function add_discretized_trajectory!(ax3d, solved_r, color="green")
    scatter!(ax3d, solved_r[1, :], solved_r[2, :], solved_r[3, :], color=color)
end

function save_with_views!(ax3d, f, prefix)
    save(prefix*"_3d.png", f)
    
    az = [pi/2, 0, 0]
    el = [0, 0, pi/2]
    name = ["y+", "x+", "z+"]
    for (a, e, n) in zip(az, el, name)
        ax3d.azimuth = a
        ax3d.elevation = e
        save(prefix*"_"*n*".png", f)
    end
end