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

@constraint(model, v_post_maneuver' * v_post_maneuver >= 1e4)

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

@objective(model, MIN_SENSE, (ΔV' * ΔV))

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
solved_r_maneuver = value.(r_maneuver)
solved_v_post_maneuver = value.(v_post_maneuver)
##
plot_orbit(
    orb0, 
    rv_to_kepler(solved_r_maneuver, solved_v_post_maneuver), 
    orb_final, 
    rv_to_kepler(solved_rf, solved_vf))