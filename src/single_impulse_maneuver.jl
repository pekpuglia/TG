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
## ##########################################
#designing single maneuver inversely
function final_position(orb0, deltaV, period_fraction1, period_fraction2)
    T0 = orbital_period(orb0, tbc_m0)

    r_preman, v_preman, orbp_preman = Propagators.propagate(Val(:TwoBody), period_fraction1*T0, orb0)
    orb_preman = orbp_preman.tbd.orbk

    v_postman = v_preman + deltaV
    orb_postman = rv_to_kepler(r_preman, v_postman, orb_preman.t)
    T_postman = orbital_period(orb_postman, tbc_m0)

    r_final, v_final, orbp_postman = Propagators.propagate(Val(:TwoBody), period_fraction2*T_postman, orb_postman)
    total_time = orbp_preman.tbd.Δt + orbp_postman.tbd.Δt
    r_final, total_time
end
##
orb0 = KeplerianElements(
                  date_to_jd(2023, 1, 1, 0, 0, 0),
                  30000e3,
                  0.101111,
                  20 |> deg2rad,
                  50    |> deg2rad,
                  30     |> deg2rad,
                  70     |> deg2rad
)
deltaV = [1000, 2000, 0]
r_final, total_time = final_position(orb0, deltaV, 0.25, 0.3)
# ##
plot_orbit(orb0, orb_postman, orb_final)
## solving maneuver
model, r_maneuver, v_post_maneuver, rf, vf = single_maneuver_model(orb0, r_final, 1.3total_time);
##
optimize!(model)
##
solved_rf = value.(rf)
solved_vf = value.(vf)

solved_r_maneuver = value.(r_maneuver)
solved_v_post_maneuver = value.(v_post_maneuver)

plot_orbit(
orb0, 
rv_to_kepler(solved_r_maneuver, solved_v_post_maneuver), 
orb_final, 
rv_to_kepler(solved_rf, solved_vf)) |> display

objective_value(model), value(model[:Δt_maneuver])
##
