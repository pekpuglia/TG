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
    0.01,
    1.5 |> deg2rad,
    1.5    |> deg2rad,
    1.5     |> deg2rad,
    1.5     |> deg2rad
)

deltaV = [0, -1000, 0]

r_final, total_time, orb_postman, orb_final = final_position(
    orb0, 
    deltaV, 
    0.5, 
    0.7)

plot_orbit(orb0, orb_postman, orb_final)
##
function single_maneuver_model_fix(orb0, r_final, total_time)
    ## auxiliary parameters
    Vesc = √(2tbc_m0 / EARTH_EQUATORIAL_RADIUS)
    Vmin = 100.0

    moon_distance = 384400.e2 

    model = Model(
        optimizer_with_attributes(Ipopt.Optimizer,
            "max_iter" => 10000)
    )

    #control variables
    Δt_maneuver = @variable(model, Δt_maneuver, start=0.5*total_time)
    @constraint(model, 0 <= Δt_maneuver <= total_time)  #set start value of constraints?

    #colocar no referencial local da v_pre_maneuver
    ΔVmag = @variable(model, 0 <= ΔV <= Vesc, start=1.0)
    
    ΔVdir = @variable(model, -1 <= ΔVdir[1:3] <= 1, start = 0.0)
    set_start_value(ΔVdir[1], 1.0)

    @constraint(model, ΔVdir' * ΔVdir == 1)

    coast_r, coast_v, coast_t = add_coast_operators!(model)
    
    r0, v0 = kepler_to_rv(orb0)

    #first coast
    r_maneuver = coast_r(r0, v0, orb0.t, Δt_maneuver)

    @constraint(model, EARTH_EQUATORIAL_RADIUS^2 <= r_maneuver' * r_maneuver <= moon_distance^2)

    v_pre_maneuver = coast_v(r0, v0, orb0.t, Δt_maneuver)

    @constraints(model, begin
        Vesc >= v_pre_maneuver[1] >= -Vesc
        Vesc >= v_pre_maneuver[2] >= -Vesc 
        Vesc >= v_pre_maneuver[3] >= -Vesc 
    end)
    @constraint(model, Vesc^2 >= v_pre_maneuver' * v_pre_maneuver >= Vmin^2)
    @constraint(model, (v_pre_maneuver' * v_pre_maneuver) * √(r_maneuver' * r_maneuver) / (2tbc_m0) <= 1-1e-6)

    time_maneuver = coast_t(r0, v0, orb0.t, Δt_maneuver)
    
    #delta V is in tangential reference frame
    #x ∥ v
    #z ∥ h = r × v
    #y = z × x

    x_tang = v_pre_maneuver ./ √(v_pre_maneuver' * v_pre_maneuver)

    h = cross(r_maneuver, v_pre_maneuver)

    z_tang = h ./ √(h' * h)

    y_tang = cross(z_tang, x_tang)

    ΔV = ΔVmag * ΔVdir

    v_post_maneuver = v_pre_maneuver + [x_tang y_tang z_tang] * ΔV

    @constraints(model, begin
        Vesc >= v_post_maneuver[1] >= -Vesc
        Vesc >= v_post_maneuver[2] >= -Vesc 
        Vesc >= v_post_maneuver[3] >= -Vesc 
    end)
    @constraint(model, Vesc^2 >= v_post_maneuver' * v_post_maneuver >= Vmin^2)
    @constraint(model, (v_post_maneuver' * v_post_maneuver) * √(r_maneuver' * r_maneuver) / (2tbc_m0) <= 1-1e-6)
    
    #second coast
    rf = coast_r(r_maneuver, v_post_maneuver, time_maneuver, total_time - Δt_maneuver)

    @constraint(model, EARTH_EQUATORIAL_RADIUS^2 <= rf' * rf <= moon_distance^2)

    vf = coast_v(r_maneuver, v_post_maneuver, time_maneuver, total_time - Δt_maneuver)
    
    @constraints(model, begin
        Vesc >= vf[1] >= -Vesc
        Vesc >= vf[2] >= -Vesc 
        Vesc >= vf[3] >= -Vesc 
    end)
    @constraint(model, Vesc^2 >= vf' * vf >= Vmin^2)
    @constraint(model, (vf' * vf) * √(rf' * rf) / (2tbc_m0) <= 1-1e-6)

    
    @constraints(model, begin
        -100.0 <= rf[1] - r_final[1] <= 100.0
        -100.0 <= rf[2] - r_final[2] <= 100.0
        -100.0 <= rf[3] - r_final[3] <= 100.0
    end)
    
    @objective(model, MIN_SENSE, ΔVmag)

    model, r_maneuver, v_post_maneuver, rf, vf
end
## solving maneuver
model, r_maneuver, v_post_maneuver, rf, vf = single_maneuver_model_fix(
    orb0, 
    r_final, 
    total_time);
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
