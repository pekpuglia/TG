using SatelliteToolboxBase
using SatelliteToolboxPropagators
using Optimization
using OptimizationOptimJL
using Optim
using Setfield
using ForwardDiff
include("TG.jl")
using .TG
##
#designing single maneuver inversely
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
## optimize [t; ΔV] for this maneuver
noNaNs(x::Real) = true
noNaNs(x::ForwardDiff.Dual) = !any(isnan, ForwardDiff.partials(x))
function final_state_residuals(maneuver_params, prob_params)
    orb0       = prob_params[1]
    total_time = prob_params[2]
    r_target   = prob_params[3]
    v_target   = prob_params[4]

    r_norm = √sum(r_target' * r_target)
    v_norm = √sum(v_target' * v_target)

    #normalized time
    maneuver_time = maneuver_params[1] * total_time
    maneuver_deltaV = maneuver_params[2:4]


    propagator_preman = Propagators.init(Val(:TwoBody), orb0, propagation_type=Number)

    r_preman, v_preman = Propagators.propagate!(propagator_preman, maneuver_time)
    v_postman = v_preman + maneuver_deltaV

    @assert noNaNs(r_preman[1]) "NaN found"
    @assert noNaNs(r_preman[2]) "NaN found"
    @assert noNaNs(r_preman[3]) "NaN found"
    @assert noNaNs(v_postman[1]) "NaN found"
    @assert noNaNs(v_postman[2]) "NaN found"
    @assert noNaNs(v_postman[3]) "NaN found"

    orb_postman = rv_to_kepler(r_preman, v_postman, propagator_preman.tbd.orbk.t)
    
    @assert noNaNs(orb_postman.a) "NaN found"
    @assert noNaNs(orb_postman.e) "NaN found"
    @assert noNaNs(orb_postman.f) "NaN found"
    @assert noNaNs(orb_postman.i) "NaN found"
    @assert noNaNs(orb_postman.t) "NaN found"
    @assert noNaNs(orb_postman.Ω) "NaN found"
    @assert noNaNs(orb_postman.ω) "NaN found"

    propagator_postman = Propagators.init(Val(:TwoBody), orb_postman, propagation_type=Number)

    r_final, v_final = Propagators.propagate!(propagator_postman, total_time-maneuver_time)
    
    sum((r_final/r_norm - r_target/r_norm).^2) + sum((v_final/v_norm - v_target/v_norm).^2)
end
##
prob_params = [
    orb0,
    T0/4+0.3*T_postman,
    r_final, v_final
]
##
prob_answer = [
    T0/4 / (prob_params[2]);
    deltaV
]
##
final_state_residuals(prob_answer, prob_params)
##
man0 = [
    0.0
    0.0
    0.0
    0.0
]
##
final_state_residuals(man0, prob_params)
##
ForwardDiff.gradient(m -> final_state_residuals(m, prob_params), man0)
##
f = OptimizationFunction(final_state_residuals, Optimization.AutoForwardDiff())
##
prob = OptimizationProblem(f, man0, prob_params)
##
sol = solve(prob, Newton())