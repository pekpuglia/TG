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

function final_state(maneuver_params, prob_params)
    ri, vi, ti       = prob_params[1]
    total_time = prob_params[2]

    #normalized time
    maneuver_delta_t = maneuver_params[1] * total_time
    maneuver_deltaV = maneuver_params[2:4] * 1000

    first_coast_ret = propagate_coast(ri..., vi..., ti, maneuver_delta_t)
    r_preman, v_preman, tman = first_coast_ret[1:3], first_coast_ret[4:6], first_coast_ret[7]
    v_postman = v_preman + maneuver_deltaV
    
    final_coast_ret = propagate_coast(r_preman..., v_postman..., tman, total_time-maneuver_delta_t)

    r_final, v_final, tfinal = final_coast_ret[1:3], final_coast_ret[4:6], final_coast_ret[7]

    [r_final; v_final]
end
##
prob_params = [
    (r0, v0, orb0.t),
    T0/4+0.3*T_postman,
    r_final, v_final
]
total_time = prob_params[2]
##
prob_answer = [
    T0/4 / (prob_params[2]);
    deltaV / 1000
]
##
final_state(prob_answer, prob_params)
##
man0 = [
    0.0
    0.0
    0.0
    0.0
]
##
final_state(man0, prob_params)
##
ForwardDiff.jacobian(m -> final_state(m, prob_params), man0)
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

memoized_final_state = memoize(
    (t, dVx, dVy, dVz) -> final_state([t;dVx;dVy;dVz], prob_params), 6)
##
#for some reason this is required 
ForwardDiff.gradient(x -> memoized_propagate_coast[7](x...), [r0..., v0..., orb0.t, T0/4])
##
r_norm = √sum(r_final' * r_final)
v_norm = √sum(v_final' * v_final)
# time_scale = total_time
# deltaV_scale = v_norm
# objective_scaling = 10000.0
##
model = Model(Ipopt.Optimizer)

#control variables
@variable(model, Δt_maneuver, start=0.0)
# @constraint(model, 0 <= Δt_maneuver <= 1)

@variable(model, ΔV[i = 1:3], start=0.0)

#relevant state variables
@variable(model, r_maneuver[i=1:3], start=r0[i])
@variable(model, v_pre_maneuver[i=1:3], start=v0[i])
@variable(model, v_post_maneuver[i=1:3], start=v0[i])
@variable(model, time_maneuver, start = orb0.t)

@variable(model, rf[i=1:3], start=r0[i])
@variable(model, vf[i=1:3], start=v0[i])

@operator(model, coast_position_x, 8, memoized_propagate_coast[1])
@operator(model, coast_position_y, 8, memoized_propagate_coast[2])
@operator(model, coast_position_z, 8, memoized_propagate_coast[3])
@operator(model, coast_velocity_x, 8, memoized_propagate_coast[4])
@operator(model, coast_velocity_y, 8, memoized_propagate_coast[5])
@operator(model, coast_velocity_z, 8, memoized_propagate_coast[6])
@operator(model, coast_time      , 8, memoized_propagate_coast[7])

#first coast
@constraints(model, begin
    r_maneuver[1]     == coast_position_x(r0..., v0..., orb0.t, Δt_maneuver)
    r_maneuver[2]     == coast_position_y(r0..., v0..., orb0.t, Δt_maneuver)
    r_maneuver[3]     == coast_position_z(r0..., v0..., orb0.t, Δt_maneuver)
    v_pre_maneuver[1] == coast_velocity_x(r0..., v0..., orb0.t, Δt_maneuver)
    v_pre_maneuver[2] == coast_velocity_y(r0..., v0..., orb0.t, Δt_maneuver)
    v_pre_maneuver[3] == coast_velocity_z(r0..., v0..., orb0.t, Δt_maneuver)
    time_maneuver     == coast_time(      r0..., v0..., orb0.t, Δt_maneuver)
end)

@constraints(model, begin
    v_post_maneuver == v_pre_maneuver + ΔV
end)

# @constraint(model, -tbc_m0 / √(r_maneuver' * r_maneuver) + 1/2 * v_post_maneuver' * v_post_maneuver <= 0)

# h_post_maneuver = cross(r_maneuver/r_norm, v_post_maneuver/v_norm)
# @constraint(model, √(h_post_maneuver' * h_post_maneuver) >= 1e8 / (r_norm*v_norm))

#second coast
@constraints(model, begin
    rf[1] == coast_position_x(r_maneuver..., v_post_maneuver..., time_maneuver, total_time - Δt_maneuver)
    rf[2] == coast_position_y(r_maneuver..., v_post_maneuver..., time_maneuver, total_time - Δt_maneuver)
    rf[3] == coast_position_z(r_maneuver..., v_post_maneuver..., time_maneuver, total_time - Δt_maneuver)
    vf[1] == coast_velocity_x(r_maneuver..., v_post_maneuver..., time_maneuver, total_time - Δt_maneuver)
    vf[2] == coast_velocity_y(r_maneuver..., v_post_maneuver..., time_maneuver, total_time - Δt_maneuver)
    vf[3] == coast_velocity_z(r_maneuver..., v_post_maneuver..., time_maneuver, total_time - Δt_maneuver)
end)

# @constraint(model, 0 <= Δt_maneuver <= 1)
# @constraint(model, 0 <= ΔV' * ΔV <= 8)

@objective(model, MIN_SENSE, 
    sum(((rf- r_final)) .^ 2) + 
    sum(((vf - v_final)) .^ 2))
# @constraint(model, final_state([Δt_maneuver; ΔV], prob_params) == [prob_params[3]; prob_params[4]])
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
solved_rf = value.(model[:rf])
solved_vf = value.(model[:vf])
##
value(model[:time_maneuver])
##
plot_orbit(orb_final, rv_to_kepler(solved_rf, solved_vf))