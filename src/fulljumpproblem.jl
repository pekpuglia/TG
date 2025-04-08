using SatelliteToolboxBase
using SatelliteToolboxPropagators
using JuMP
using Ipopt
using Setfield
using ForwardDiff
include("TG.jl")
using .TG
## ######################################################
#       designing single maneuver inversely
#
#########################################################

orb0 = KeplerianElements(
                  date_to_jd(2023, 1, 1, 0, 0, 0),
                  7190.982e3,
                  0.001111,
                  10 |> deg2rad,
                  0.5    |> deg2rad,
                  0.5     |> deg2rad,
                  0.5     |> deg2rad
)
T0 = orbital_period(orb0, tbc_m0)
r0, v0 = kepler_to_rv(orb0)
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
# SatelliteToolbox stores time as Julian Day
episode_time = (orb_final.t - orb0.t) * 2π/EARTH_ANGULAR_SPEED
plot_orbit(orb_final)
##
prob_answer = [
    T0/4 / (T0/4+0.3*T_postman);
    deltaV
]
#####################################################################################
## optimize [t; ΔV] for this maneuver - put into TG module

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

noNaNs(x::Real) = true
noNaNs(x::ForwardDiff.Dual) = !any(isnan, ForwardDiff.partials(x))

function propagate_state(x0, y0, z0, vx0, vy0, vz0, t0, Δt)
    orb0 = rv_to_kepler([x0;y0;z0], [vx0; vy0; vz0], t0)

    @assert noNaNs(orb0.a) "NaN found"
    @assert noNaNs(orb0.e) "NaN found"
    @assert noNaNs(orb0.f) "NaN found"
    @assert noNaNs(orb0.i) "NaN found"
    @assert noNaNs(orb0.t) "NaN found"
    @assert noNaNs(orb0.Ω) "NaN found"
    @assert noNaNs(orb0.ω) "NaN found"

    propagator0 = Propagators.init(Val(:TwoBody), orb0, propagation_type=Real)

    rf, vf = Propagators.propagate!(propagator0, Δt)
    
    @assert noNaNs(rf[1]) "NaN found"
    @assert noNaNs(rf[2]) "NaN found"
    @assert noNaNs(rf[3]) "NaN found"
    @assert noNaNs(vf[1]) "NaN found"
    @assert noNaNs(vf[2]) "NaN found"
    @assert noNaNs(vf[3]) "NaN found"

    [rf; vf; propagator0.tbd.orbk.t]
end
##
propagate_state(r_preman..., v_postman..., orb_preman.t, 0.3*T_postman)
##
ForwardDiff.jacobian(x -> propagate_state(x..., 0.3*T_postman), [r_preman..., v_postman..., orb_preman.t])
##
memoized_propagate_state = memoize(propagate_state, 7)
##
[f(r_preman..., v_postman..., orb_preman.t, 0.3*T_postman) for f in memoized_propagate_state]
##
function residuals(x, y, z, vx, vy, vz, r_target, v_target)
    r_norm = √sum(r_target' * r_target)
    v_norm = √sum(v_target' * v_target)

    sum(
        ([x;y;z]/r_norm - r_target/r_norm).^2
    )+sum(
        ([vx;vy;vz]/v_norm - v_target/v_norm).^2
    )
end
##
model = Model(Ipopt.Optimizer)

propagated_position_x = @operator(model, propagated_position_x, 8, memoized_propagate_state[1])
propagated_position_y = @operator(model, propagated_position_y, 8, memoized_propagate_state[2])
propagated_position_z = @operator(model, propagated_position_z, 8, memoized_propagate_state[3])
propagated_velocity_x = @operator(model, propagated_velocity_x, 8, memoized_propagate_state[4])
propagated_velocity_y = @operator(model, propagated_velocity_y, 8, memoized_propagate_state[5])
propagated_velocity_z = @operator(model, propagated_velocity_z, 8, memoized_propagate_state[6])
propagated_time       = @operator(model, propagated_time,       8, memoized_propagate_state[7])

Δt_maneuver = @variable(model, Δt_maneuver, start=0.0)
ΔV         = @variable(model, ΔV[i = 1:3], start=0.0)

x_premaneuver  = propagated_position_x(r0..., v0..., orb0.t, Δt_maneuver)
y_premaneuver  = propagated_position_y(r0..., v0..., orb0.t, Δt_maneuver)
z_premaneuver  = propagated_position_z(r0..., v0..., orb0.t, Δt_maneuver)
vx_premaneuver = propagated_velocity_x(r0..., v0..., orb0.t, Δt_maneuver)
vy_premaneuver = propagated_velocity_y(r0..., v0..., orb0.t, Δt_maneuver)
vz_premaneuver = propagated_velocity_z(r0..., v0..., orb0.t, Δt_maneuver)
t_maneuver     = propagated_time(      r0..., v0..., orb0.t, Δt_maneuver)

vx_postmaneuver = vx_premaneuver + ΔV[1]
vy_postmaneuver = vy_premaneuver + ΔV[2]
vz_postmaneuver = vz_premaneuver + ΔV[3]

r_premaneuver  = [x_premaneuver  ; y_premaneuver  ; z_premaneuver]
v_postmaneuver = [vx_postmaneuver; vy_postmaneuver; vz_postmaneuver]

x_final  = propagated_position_x(r_premaneuver..., v_postmaneuver..., t_maneuver, episode_time - Δt_maneuver)
y_final  = propagated_position_y(r_premaneuver..., v_postmaneuver..., t_maneuver, episode_time - Δt_maneuver)
z_final  = propagated_position_z(r_premaneuver..., v_postmaneuver..., t_maneuver, episode_time - Δt_maneuver)
vx_final = propagated_velocity_x(r_premaneuver..., v_postmaneuver..., t_maneuver, episode_time - Δt_maneuver)
vy_final = propagated_velocity_y(r_premaneuver..., v_postmaneuver..., t_maneuver, episode_time - Δt_maneuver)
vz_final = propagated_velocity_z(r_premaneuver..., v_postmaneuver..., t_maneuver, episode_time - Δt_maneuver)

@constraint(model, 0 <= t_maneuver <= 1)
@constraint(model, 0 <= ΔV' * ΔV <= 2e6)

@objective(model, MIN_SENSE, residuals(
    x_final,
    y_final,
    z_final,
    vx_final,
    vy_final,
    vz_final,
    r_final, v_final
))
# @constraint(model, final_state([t_maneuver; ΔV], prob_params) == [prob_params[3]; prob_params[4]])
model
##
optimize!(model)
##
value.(all_variables(model))
##