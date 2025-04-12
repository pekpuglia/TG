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
    7000e3,
    0.01,
    1.5 |> deg2rad,
    1.5    |> deg2rad,
    1.5     |> deg2rad,
    1.5     |> deg2rad
)
r0, v0 = kepler_to_rv(orb0)

deltaV = [-1000, 0, 0]

r_final, total_time, orb_postman, orb_final, v_final = final_position(
    orb0, 
    deltaV, 
    0.2, 
    0.0
)

plot_orbit(orb0, orb_postman, orb_final)
## solving maneuver
# model, r_maneuver, v_post_maneuver, rf, vf = single_maneuver_model(
#     orb0, 
#     r_final, 
#     dir_v_final,
#     total_time);
# set_silent(model)
model = Model(
    optimizer_with_attributes(Ipopt.Optimizer,
        "max_iter" => 10000,
        "max_wall_time" => 30.0)
)

Δt_maneuver = @variable(model, Δt_maneuver, start=0.5*total_time)
# @constraint(model, 0 <= Δt_maneuver <= total_time)  #set start value of constraints?

ΔVmag = @variable(model, 0 <= ΔVmag, start=1.0)
    
ΔVdir = @variable(model, -1 <= ΔVdir[1:3] <= 1, start = 0.0)
set_start_value(ΔVdir[1], 1.0)

# @constraint(model, ΔVdir' * ΔVdir == 1)

orbparams_i = add_orbital_elements!(model)

orbparams_pre_maneuver = add_orbital_elements!(model)
# orbparams_post_maneuver = add_orbital_elements!(model)

# orbparams_f = add_orbital_elements!(model)

add_coast_set_boundaries!(model, 
    orbparams_i,
    orbparams_pre_maneuver,
    r0, v0,
    r_final, v_final - ΔVmag * ΔVdir, 
    Δt_maneuver
)

# add_coast_set_boundaries!(model,
#     orbparams_post_maneuver,
#     orbparams_f,
#     orbparams_pre_maneuver.r,
#     orbparams_post_maneuver.v + ΔVmag*ΔVdir,
#     r_final, v_final,
#     total_time - Δt_maneuver
# )

# @objective(model, MIN_SENSE, ΔVmag)

model

##
optimize!(model)
##
solved_rf = value.(orbparams_pre_maneuver.r)
solved_vf = value.(orbparams_pre_maneuver.v)

solved_deltaV = value.(ΔVmag * ΔVdir)

solved_r_maneuver = value.(orbparams_pre_maneuver.r)
solved_v_post_maneuver = value.(orbparams_pre_maneuver.v + solved_deltaV)

plot_orbit(
orb0, 
rv_to_kepler(solved_r_maneuver, solved_v_post_maneuver), 
orb_final, 
rv_to_kepler(solved_rf, solved_vf)) |> display
##
objective_value(model), value(model[:Δt_maneuver])