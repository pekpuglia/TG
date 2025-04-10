using SatelliteToolboxBase
using SatelliteToolboxPropagators
using JuMP
using Ipopt
using Setfield
using ForwardDiff
include("TG.jl")
using .TG
using LinearAlgebra
## ##########################################
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
# plot_orbit(orb0)
##
r_preman, v_preman, orbp_preman = Propagators.propagate(Val(:TwoBody), T0/4, orb0)
orb_preman = orbp_preman.tbd.orbk
deltaV = [1000, 0, 0]
v_postman = v_preman + deltaV
orb_postman = rv_to_kepler(r_preman, v_postman, orb_preman.t)
# plot_orbit(orb0, orb_postman)
## final position desired
T_postman = orbital_period(orb_postman, tbc_m0)
r_final, v_final, orbp_postman = Propagators.propagate(Val(:TwoBody), 0.3*T_postman, orb_postman)
orb_final = orbp_postman.tbd.orbk
total_time = T0/4+0.3*T_postman
# plot_orbit(orb_final)
##
plot_orbit(orb0, orb_postman, orb_final)
#######################################################
## solving maneuver

##
function maneuver_model(orb0, r_final, total_time)
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
##
model, r_maneuver, v_post_maneuver, rf, vf = maneuver_model(orb0, r_final, 2*total_time);
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
