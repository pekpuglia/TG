using GLMakie
using SatelliteToolboxBase
using Setfield

function orbit_plot_data(orbs::KeplerianElements...)
    N = 100
    θ = LinRange(0, 2π, N)

    plot_data = Dict[]

    for orb in orbs
        orbit = ((@set orb.f = theta) for theta in θ) .|> kepler_to_rv .|> first |> stack
        current_pos, current_velocity = kepler_to_rv(orb)
        velocity_arrow_base = current_pos
        velocity_arrow_tip = velocity_arrow_base + current_velocity / √sum(current_velocity .^ 2) * 0.4 * √sum(current_pos .^ 2)
        velocity_arrow_data = [velocity_arrow_base velocity_arrow_tip] |> Matrix

        push!(plot_data, Dict(
            :orbit => orbit,
            :curr_pos => current_pos,
            :vel_arrow_data => velocity_arrow_data
        ))
    end
    plot_data
end

function plot_orbit(plot_data::Vector{Dict})
    fig = Figure()
    ax3d = Axis3(fig[1, 1])
    
    orb_lines = []

    fig = Figure()
    ax3d = Axis3(fig[1, 1])

    for (data, c) in zip(plot_data, Makie.wong_colors())
        orbit = data[:orbit]
        current_pos = data[:curr_pos]
        velocity_arrow_data = data[:vel_arrow_data]
        l = lines!(ax3d, orbit[1, :], orbit[2, :], orbit[3, :], color=c)
        push!(orb_lines, l)
        scatter!(ax3d, current_pos, markersize=20, color=c)
        lines!(ax3d, velocity_arrow_data[1, :], velocity_arrow_data[2, :], velocity_arrow_data[3, :], color=c)
    end

    wireframe!(ax3d, Sphere(Point3(0.0), EARTH_EQUATORIAL_RADIUS), color=:cyan, alpha=0.3)
    fig, ax3d, orb_lines
end

function plot_orbit(orbs::KeplerianElements...)
    plot_data = orbit_plot_data(orbs...)

    plot_orbit(plot_data)
end

const VIEW_DICT = Dict(
    :Xp => [2, 3],
    :Ym => [1, 3],
    :Zp => [1, 2]
)

const LABEL_DICT = Dict(
    :Xp => ["y", "z"],
    :Ym => ["x", "z"],
    :Zp => ["x", "y"]
)

function plot_orbit_2d(view, plot_data::Vector{Dict})
    view_inds = VIEW_DICT[view]

    orb_lines = []
    
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = LABEL_DICT[view][1], ylabel=LABEL_DICT[view][2])
    
    for (data, c) in zip(plot_data, Makie.wong_colors())
        orbit = data[:orbit]
        current_pos = data[:curr_pos]
        velocity_arrow_data = data[:vel_arrow_data]
        l = lines!(ax, orbit[view_inds[1], :], orbit[view_inds[2], :], color=c)
        push!(orb_lines, l)
        scatter!(ax, current_pos[view_inds[1]], current_pos[view_inds[2]], markersize=20, color=c)
        lines!(ax, velocity_arrow_data[view_inds[1], :], velocity_arrow_data[view_inds[2], :], color=c)
    end
    
    wireframe!(ax, Circle(Point2(0.0), EARTH_EQUATORIAL_RADIUS), color=:cyan, alpha=0.3)
    
    fig , ax, orb_lines
end

function plot_orbit_2d(view, orbs::KeplerianElements...)
    plot_data = orbit_plot_data(orbs...)

    plot_orbit_2d(view, plot_data)
end

function add_discretized_trajectory!(ax3d, solved_r, color="green")
    scatter!(ax3d, solved_r[1, :], solved_r[2, :], solved_r[3, :], color=color)
end

function add_discretized_trajectory_2d!(ax, view, solved_r, color="green")
    view_inds = VIEW_DICT[view]
    scatter!(ax, solved_r[view_inds[1], :], solved_r[view_inds[2], :], color=color)
end

function transfer_plot_data(solved_transfer::Transfer, scaling=1e3)
    #get first position
    last_r = if solved_transfer.sequence[1] isa Coast
        solved_transfer.sequence[1].rcoast[:, 1]
    else
        solved_transfer.sequence[2].rcoast[:, 1]
    end

    coast_points = []

    imp_data = []

    for el in solved_transfer.sequence
        if el isa Coast
            # c= add_discretized_trajectory!(ax3d, el.rcoast)
            push!(coast_points, el.rcoast)
            last_r = el.rcoast[:, end]
        else
            # scatter!(ax3d, last_r...)
            # #add line on impulse
            dir = el.deltaVdir
            mag = el.deltaVmag
            
            arrow_data = [last_r last_r+mag*dir*scaling]
            # i = lines!(ax3d, arrow_data[1, :], arrow_data[2, :], arrow_data[3, :], color="red")
            # push!(imp_arrows, i)
            push!(imp_data, (last_r, arrow_data))
        end
    end
    (coast_points, imp_data)
end

function add_transfer!(ax3d, tpd::Tuple)
    #get first position
    # last_r = if solved_transfer.sequence[1] isa Coast
    #     solved_transfer.sequence[1].rcoast[:, 1]
    # else
    #     solved_transfer.sequence[2].rcoast[:, 1]
    # end

    coast_points, imp_data = tpd

    coast_plots = []

    for cp = coast_points
        c= add_discretized_trajectory!(ax3d, cp)
        push!(coast_plots, c)
    end

    imp_arrows = []

    for id = imp_data
        r, arrow_data = id
        scatter!(ax3d, r...)
        i = lines!(ax3d, arrow_data[1, :], arrow_data[2, :], arrow_data[3, :], color="red")
        push!(imp_arrows, i)
    end
    coast_plots, imp_arrows
end

function add_transfer!(ax3d, solved_transfer::Transfer, scaling=1e3)
    tpd = transfer_plot_data(solved_transfer, scaling)

    add_transfer!(ax3d, tpd)
end

function add_transfer_2d!(ax, view, tpd::Tuple)
    view_inds = VIEW_DICT[view]

    coast_points, imp_data = tpd

    coast_plots = []

    for cp = coast_points
        c= add_discretized_trajectory_2d!(ax, view, cp)
        push!(coast_plots, c)
    end

    imp_arrows = []

    for id = imp_data
        r, arrow_data = id
        scatter!(ax, r[view_inds]...)
        i = lines!(ax, arrow_data[view_inds[1], :], arrow_data[view_inds[2], :], color="red")
        push!(imp_arrows, i)
    end
    coast_plots, imp_arrows
end

function add_transfer_2d!(ax, view, solved_transfer::Transfer, scaling=1e3)
    tpd = transfer_plot_data(solved_transfer, scaling)

    add_transfer_2d!(ax, view, tpd)
end

function plot_transfer(orb1::KeplerianElements, orb2::KeplerianElements, solved_transfer::Transfer, scaling)
    f, ax3d, orb_lines = plot_orbit(orb1, orb2)
    coast_ps, ia = add_transfer!(ax3d, solved_transfer, scaling)
    Legend(f[1, 2], [orb_lines..., coast_ps[1], ia[1]], ["Initial orbit", "Final orbit", "Coasting arc", "Impulse"], position = (0.8, 0.9))
    f
end

function plot_transfer_2d(view, orb1::KeplerianElements, orb2::KeplerianElements, solved_transfer::Transfer, scaling)
    f, ax, orb_lines = plot_orbit_2d(view, orb1, orb2)
    coast_ps, ia = add_transfer_2d!(ax, view, solved_transfer, scaling)
    # Legend(f[1, 2], [orb_lines..., coast_ps[1], ia[1]], ["Initial orbit", "Final orbit", "Coasting arc", "Impulse"], position = (0.8, 0.9))
    f
end

function save_with_views!(prefix, orb1::KeplerianElements, orb2::KeplerianElements, solved_transfer::Transfer, scaling)
    #3d plot
    f = plot_transfer(orb1, orb2, solved_transfer, scaling)
    save(prefix*"_3d.png", f, px_per_unit = 300/96)
    
    #2d plots
    fxp = plot_transfer_2d(:Xp, HOHMANN_START, HOHMANN_END, solved_transfer, scaling)
    save(prefix*"_x+.png", fxp, px_per_unit = 300/96)
    fym = plot_transfer_2d(:Ym, HOHMANN_START, HOHMANN_END, solved_transfer, scaling)
    save(prefix*"_y-.png", fym, px_per_unit = 300/96)
    fzp = plot_transfer_2d(:Zp, HOHMANN_START, HOHMANN_END, solved_transfer, scaling)
    save(prefix*"_z+.png", fzp, px_per_unit = 300/96)
end

function plot_primer_vector(transfer::Transfer, tspan_ppdot::Tuple; name, label, style)
    
    f = Figure()
    ax1 = Axis(f[1, 1], xlabel = "t (s)", ylabel = "|p|", title="Maneuver: $name, $(transfer_type(transfer))")
    ax2 = Axis(f[2, 1], xlabel = "t (s)", ylabel = L"d |p| / dt")
    
    
    plot_primer_vector!(f, ax1, ax2, transfer, tspan_ppdot; label, style)
end

function plot_primer_vector!(f, ax1, ax2, transfer::Transfer, tspan_ppdot::Tuple; label, style)
    tspan, ppdot = tspan_ppdot
    normp = norm.(eachcol(ppdot[1:3, :]))
    normpdot = [dot(ppdoti[1:3], ppdoti[4:6]) / norm(ppdoti[1:3]) for ppdoti in eachcol(ppdot)]
    
    imp_ts = impulse_times(transfer)
    
    lines!(ax1, tspan, normp, label = label, linestyle=style)
    vlines!(ax1, imp_ts, linestyle=:dash, color=:gray)
    
    lines!(ax2, tspan, normpdot, label = label, linestyle=style)
    vlines!(ax2, imp_ts, linestyle=:dash, color=:gray)

    f, (ax1, ax2)
end