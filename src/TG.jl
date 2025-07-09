module TG

using ForwardDiff
using GLMakie
using Ipopt
using JuMP
using LinearAlgebra
using SatelliteToolboxBase
using SatelliteToolboxPropagators
using Setfield

export plot_orbit, add_discretized_trajectory!, save_with_views!

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

function add_discretized_trajectory!(ax3d, solved_r)
    scatter!(ax3d, solved_r[1, :], solved_r[2, :], solved_r[3, :], color="green")
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

export orbital_period
function orbital_period(orb::KeplerianElements, m0)
    orbital_period(orb.a,m0)
end

orbital_period(a, m0) = 2π√(a^3/m0)

include("legacy.jl")

include("./glandorf.jl")
include("integrators.jl")

export two_body_dyn
function two_body_dyn(X, mu)
    r = X[1:3]
    v = X[4:6]
    [
        v
        -mu/(√(r'*r)^3)*r
    ]
end

export add_coast_segment
function add_coast_segment(model, deltat, N, ind; dyn=(X -> two_body_dyn(X, GM_EARTH)), integrator=RK4)
    #[(x, y, z), time instants]
    r = @variable(model, [1:3, 1:N+1], base_name="r_coast_$ind")
    v = @variable(model, [1:3, 1:N+1], base_name="v_coast_$ind")


    for i = 1:(N+1)
        # @constraint(model, r[:, i]' * r[:, i] >= EARTH_EQUATORIAL_RADIUS^2)
        # @constraint(model, cross(r[:, i], v[:, i])[3] >=0)
    end

    for i = 1:N
        @constraint(model, integrator(dyn, [r[:, i]; v[:, i]], deltat/N) .== [r[:, i+1]; v[:, i+1]])
    end

    r, v
end

export p0dot_tpbvp, ppdot_deltavs, diagnose_ppdot
function p0dot_tpbvp(p0, pf, delta_t, prop)
    Phi = TG.Phi_time(prop, delta_t)

    M = Phi[1:3, 1:3]
    N = Phi[1:3, 4:6]

    if abs(det(N)) <= 1e-10
        @warn "N singular"
        #CHECK THIS IS TRUE
        #only works for coplanar transfers?
        r1, v1 = kepler_to_rv(prop.tbd.orb₀)
        A = N * [r1 v1]
        b = pf - M * p0
        p0dot = [r1 v1] * (A \ b)
    else
        p0dot = N \ (pf - M * p0)
    end

    p0dot
end


function ppdot_deltavs(transfer_propagator, deltav1, deltav2, delta_t, N)
    p0 = deltav1 / norm(deltav1)
    pf = deltav2 / norm(deltav2)
    p0dot = p0dot_tpbvp(p0, pf, delta_t, transfer_propagator)
    tspan = range(0, delta_t, N)
    ppdot = [TG.Phi_time(transfer_propagator, t) * [p0; p0dot] for t in tspan]
    tspan, ppdot
end

#doesn't account for continuity yet
@enum PRIMER_DIAGNOSTIC IC_FC IC_LA ED_FC ED_LA MID OPT
function diagnose_ppdot(normp, normp_dot, rtol = 1e-4)
    normp_dot0 = normp_dot[1]
    normp_dotf = normp_dot[end]
    
    max_normp_dot = maximum(abs.(normp_dot))
    tol = rtol*max_normp_dot

    if normp_dot0 > tol && normp_dotf < -tol
        IC_FC
    elseif normp_dot0 > tol && normp_dotf > tol
        IC_LA
    elseif normp_dot0 < -tol && normp_dotf < -tol
        ED_FC
    elseif normp_dot0 < -tol && normp_dotf > tol
        ED_LA
    elseif any(normp .> 1 + rtol)
        MID
    else
        OPT
    end
end

end # module TG
