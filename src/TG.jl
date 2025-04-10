module TG

using ForwardDiff
using GLMakie
using Ipopt
using JuMP
using LinearAlgebra
using SatelliteToolboxBase
using SatelliteToolboxPropagators
using Setfield

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

export add_coast_operators!
function add_coast_operators!(model)
    memoized_propagate_coast = memoize(propagate_coast, 7)
    
    #for some reason this is required to avoid error no method matching getindex(::Nothing, ::Int64)
    #no f clue why
    #maybe something to do with the memoization initialization?
    Vcirc = [0, √(tbc_m0 / EARTH_EQUATORIAL_RADIUS), 0]
    rcirc = [EARTH_EQUATORIAL_RADIUS, 0, 0]
    Tcirc = orbital_period(rv_to_kepler(rcirc, Vcirc), tbc_m0)
    ForwardDiff.gradient(x -> memoized_propagate_coast[7](x...), [rcirc..., Vcirc..., 0.0, Tcirc])

    coast_position_x = @operator(model, coast_position_x, 8, memoized_propagate_coast[1])
    coast_position_y = @operator(model, coast_position_y, 8, memoized_propagate_coast[2])
    coast_position_z = @operator(model, coast_position_z, 8, memoized_propagate_coast[3])
    coast_velocity_x = @operator(model, coast_velocity_x, 8, memoized_propagate_coast[4])
    coast_velocity_y = @operator(model, coast_velocity_y, 8, memoized_propagate_coast[5])
    coast_velocity_z = @operator(model, coast_velocity_z, 8, memoized_propagate_coast[6])
    coast_time       = @operator(model, coast_time      , 8, memoized_propagate_coast[7])

    return (
        (r, v, t, dt) -> [
            coast_position_x(r..., v..., t, dt)
            coast_position_y(r..., v..., t, dt)
            coast_position_z(r..., v..., t, dt)
        ],
        (r, v, t, dt) -> [
            coast_velocity_x(r..., v..., t, dt)
            coast_velocity_y(r..., v..., t, dt)
            coast_velocity_z(r..., v..., t, dt)
        ],
        (r, v, t, dt) -> coast_time(r..., v..., t, dt)
    )
end

export single_maneuver_model
function single_maneuver_model(orb0, r_final, total_time)
    ## auxiliary parameters
    Vesc = √(2tbc_m0 / EARTH_EQUATORIAL_RADIUS)
    Vmin = 100.0
    model = Model(Ipopt.Optimizer)
    #control variables
    Δt_maneuver = @variable(model, Δt_maneuver, start=0.0)
    @constraint(model, 0 <= Δt_maneuver <= total_time)

    ΔV = @variable(model, -Vesc <= ΔV[i = 1:3] <= Vesc, start=1.0)
    
    coast_r, coast_v, coast_t = add_coast_operators!(model)
    
    r0, v0 = kepler_to_rv(orb0)

    #first coast
    r_maneuver = coast_r(r0, v0, orb0.t, Δt_maneuver)
    v_pre_maneuver = coast_v(r0, v0, orb0.t, Δt_maneuver)
    time_maneuver = coast_t(r0, v0, orb0.t, Δt_maneuver)
    
    v_post_maneuver = v_pre_maneuver + ΔV
    
    @constraint(model, Vesc^2 >= v_post_maneuver' * v_post_maneuver >= Vmin^2)
    @constraint(model, (v_post_maneuver' * v_post_maneuver) * √(r_maneuver' * r_maneuver) / (2tbc_m0) <= 1-1e-6)
    
    #second coast
    rf = coast_r(r_maneuver, v_post_maneuver, time_maneuver, total_time - Δt_maneuver)
    vf = coast_v(r_maneuver, v_post_maneuver, time_maneuver, total_time - Δt_maneuver)
    
    @constraints(model, begin
        rf[1] == r_final[1]
        rf[2] == r_final[2]
        rf[3] == r_final[3]
    end)
    
    @objective(model, MIN_SENSE, √(ΔV' * ΔV))

    model, r_maneuver, v_post_maneuver, rf, vf
end

end # module TG
