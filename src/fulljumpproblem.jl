using SatelliteToolboxBase
using SatelliteToolboxPropagators
using JuMP
using Ipopt
using Setfield
using ForwardDiff
include("TG.jl")
using .TG
using LinearAlgebra
##
#designing single maneuver inversely
orb0 = KeplerianElements(
                  date_to_jd(2023, 1, 1, 0, 0, 0),
                  8000e3,
                  0.101111,
                  20 |> deg2rad,
                  50    |> deg2rad,
                  30     |> deg2rad,
                  70     |> deg2rad
)
T0 = orbital_period(orb0, tbc_m0)
r0, v0 = kepler_to_rv(orb0)
plot_orbit(orb0)
##
r_preman, v_preman, orbp_preman = Propagators.propagate(Val(:TwoBody), T0/4, orb0)
orb_preman = orbp_preman.tbd.orbk
deltaV = [-1000, 0, 0]
v_postman = v_preman + deltaV
orb_postman = rv_to_kepler(r_preman, v_postman, orb_preman.t)
plot_orbit(orb0, orb_postman)
## final position desired
T_postman = orbital_period(orb_postman, tbc_m0)
r_final, v_final, orbp_postman = Propagators.propagate(Val(:TwoBody), 0.3*T_postman, orb_postman)
orb_final = orbp_postman.tbd.orbk
total_time = T0/4+0.3*T_postman
plot_orbit(orb_final)
##
plot_orbit(orb0, orb_postman, orb_final)
## optimize [t; ΔV] for this maneuver
noNaNs(x::Real) = true
noNaNs(x::ForwardDiff.Dual) = !any(isnan, ForwardDiff.partials(x))

function propagate_coast(xi, yi, zi, vxi, vyi, vzi, ti, deltat)
    # println(-tbc_m0 / (√(xi^2+yi^2+zi^2)) + (vxi^2+vyi^2+vzi^2)/2)
    orbi = rv_to_kepler([xi, yi, zi], [vxi, vyi, vzi], ti)
    propagator = Propagators.init(Val(:TwoBody), orbi, propagation_type=Real)
    r, v = Propagators.propagate!(propagator, deltat)
    [r; v; propagator.tbd.orbk.t]
end
##
## https://jump.dev/JuMP.jl/stable/tutorials/nonlinear/tips_and_tricks/#User-defined-operators-with-vector-outputs
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

memoized_propagate_coast = memoize(propagate_coast, 7)
#for some reason this is required 
ForwardDiff.gradient(x -> memoized_propagate_coast[7](x...), [r0..., v0..., orb0.t, T0/4])
##
model = Model(Ipopt.Optimizer)
#control variables
@variable(model, Δt_maneuver, start=0.0)
@variable(model, ΔV[i = 1:3], start=0.0)

@operator(model, coast_position_x, 8, memoized_propagate_coast[1])
@operator(model, coast_position_y, 8, memoized_propagate_coast[2])
@operator(model, coast_position_z, 8, memoized_propagate_coast[3])
@operator(model, coast_velocity_x, 8, memoized_propagate_coast[4])
@operator(model, coast_velocity_y, 8, memoized_propagate_coast[5])
@operator(model, coast_velocity_z, 8, memoized_propagate_coast[6])
@operator(model, coast_time      , 8, memoized_propagate_coast[7])

#first coast
r_maneuver = [
    coast_position_x(r0..., v0..., orb0.t, Δt_maneuver)
    coast_position_y(r0..., v0..., orb0.t, Δt_maneuver)
    coast_position_z(r0..., v0..., orb0.t, Δt_maneuver)
]
v_pre_maneuver = [
    coast_velocity_x(r0..., v0..., orb0.t, Δt_maneuver)
    coast_velocity_y(r0..., v0..., orb0.t, Δt_maneuver)
    coast_velocity_z(r0..., v0..., orb0.t, Δt_maneuver)
]
time_maneuver = coast_time(r0..., v0..., orb0.t, Δt_maneuver)

v_post_maneuver = v_pre_maneuver + ΔV

#second coast
rf = [
    coast_position_x(r_maneuver..., v_post_maneuver..., time_maneuver, total_time - Δt_maneuver)
    coast_position_y(r_maneuver..., v_post_maneuver..., time_maneuver, total_time - Δt_maneuver)
    coast_position_z(r_maneuver..., v_post_maneuver..., time_maneuver, total_time - Δt_maneuver)
]
vf = [
    coast_velocity_x(r_maneuver..., v_post_maneuver..., time_maneuver, total_time - Δt_maneuver)
    coast_velocity_y(r_maneuver..., v_post_maneuver..., time_maneuver, total_time - Δt_maneuver)
    coast_velocity_z(r_maneuver..., v_post_maneuver..., time_maneuver, total_time - Δt_maneuver)
]

@constraints(model, begin
    rf[1] == r_final[1]
    rf[2] == r_final[2]
    rf[3] == r_final[3]
    vf[1] == v_final[1]
    # vf[2] == v_final[2]
    # vf[3] == v_final[3]
end)

model
##
optimize!(model)
##
objective_value(model)
##
value(model[:Δt_maneuver])
##
value.(model[:ΔV])
##
solved_rf = value.(rf)
solved_vf = value.(vf)
##
plot_orbit(orb0, orb_final, rv_to_kepler(solved_rf, solved_vf))