module TG

using SatelliteToolboxBase
using SatelliteToolboxPropagators
using GLMakie
using Setfield
using LinearAlgebra

export plot_orbit

function plot_orbit(orbs::KeplerianElements...)
    N = 100
    θ = LinRange(0, 2π, N)

    fig = Figure()
    ax3d = Axis3(fig[1, 1])
    
    for orb in orbs
        orbit = ((@set orb.f = theta) for theta in θ) .|> kepler_to_rv .|> first |> stack
        current_pos, current_velocity = kepler_to_rv(orb)
        velocity_arrow_base = current_pos
        velocity_arrow_tip = velocity_arrow_base + current_velocity / √sum(current_velocity .^ 2) * 0.4 * √sum(current_pos .^ 2)
        velocity_arrow_data = [velocity_arrow_base velocity_arrow_tip] |> Matrix
        lines!(ax3d, orbit[1, :], orbit[2, :], orbit[3, :], color=:red)
        scatter!(ax3d, current_pos, markersize=20)
        lines!(ax3d, velocity_arrow_data[1, :], velocity_arrow_data[2, :], velocity_arrow_data[3, :])
    end
    
    wireframe!(ax3d, Sphere(Point3(0.0), EARTH_EQUATORIAL_RADIUS), color=:cyan, alpha=0.3)
    fig
end

export orbital_period
function orbital_period(orb::KeplerianElements, m0)
    2π√(orb.a^3/m0)
end

export propagate_coast
function propagate_coast(xi, yi, zi, vxi, vyi, vzi, ti, deltat)
    r = [xi; yi; zi]
    v = [vxi; vyi; vzi]
    #readd nan tests? safe, non safe version?
    @assert norm(v) > 1e3 "Velocity too small"
    @assert norm(cross(r, v)) > 1e9 "Angular momentum too small"
    energy = -tbc_m0/norm(r) + (v'*v)/2
    @assert energy < 0 "Energy non negative: $energy, r = $r, v = $v"

    # println(-tbc_m0 / (√(xi^2+yi^2+zi^2)) + (vxi^2+vyi^2+vzi^2)/2)
    orbi = rv_to_kepler([xi, yi, zi], [vxi, vyi, vzi], ti)
    propagator = Propagators.init(Val(:TwoBody), orbi, propagation_type=Real)
    r, v = Propagators.propagate!(propagator, deltat)
    [r; v; propagator.tbd.orbk.t]
end

export memoize
# https://jump.dev/JuMP.jl/stable/tutorials/nonlinear/tips_and_tricks/#User-defined-operators-with-vector-outputs
"""
    memoize(foo::Function, n_outputs::Int)

Take a function `foo` and return a vector of length `n_outputs`, where element
`i` is a function that returns the equivalent of `foo(x...)[i]`.

To avoid duplication of work, cache the most-recent evaluations of `foo`.
Because `foo_i` is auto-differentiated with ForwardDiff, our cache needs to
work when `x` is a `Float64` and a `ForwardDiff.Dual`.
"""
function memoize(foo::Function, n_outputs::Int)
    last_x, last_f = nothing, nothing
    last_dx, last_dfdx = nothing, nothing
    function foo_i(i, x::T...) where {T<:Real}
        if T == Float64
            if x !== last_x
                last_x, last_f = x, foo(x...)
            end
            return last_f[i]::T
        else
            if x !== last_dx
                last_dx, last_dfdx = x, foo(x...)
            end
            return last_dfdx[i]::T
        end
    end
    return [(x...) -> foo_i(i, x...) for i in 1:n_outputs]
end
end # module TG
