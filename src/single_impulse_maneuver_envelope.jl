using SatelliteToolboxBase
using SatelliteToolboxPropagators
using JuMP
using Ipopt
using Setfield
using ForwardDiff
include("TG.jl")
using .TG
using LinearAlgebra
using LazyGrids
using Random
## ##########################################
struct InputParams
    a
    e
    f
    dvx
    dvy
    dvz
    fT1
    fT2
end
##
  a_list = 7000e3:1000e3:20000e3
  e_list = [0.001, 0.01, 0.02, 0.1, 0.4]
  f_list = 0:4
dvx_list = [-1200, -1100, -1000, -500, 0, 500, 1000, 1100, 1200]
dvy_list = dvx_list
dvz_list = dvx_list
fT1_list = 0.1:0.2:0.9
fT2_list = 0.1:0.2:0.9
##
input_grid = ndgrid(a_list, e_list, f_list, dvx_list, dvy_list, dvz_list, fT1_list, fT2_list);
input_params_array = InputParams.(input_grid...);
input_params_list = reshape(input_params_array, :);
shuffle!(input_params_list)
##
struct OutputParams
    inp::InputParams
    error
    exec_time
    is_solved_and_feasible
    deltaV
    time_man
end
##
Nsamples = 3
random_inputs = input_params_list[1:Nsamples]

outputs = []

out_file = "./src/single_impulse_out"

open(out_file, "a") do file
    
    for (i, inp) in enumerate(random_inputs)

        if mod(i, round(Nsamples/10)) == 0
            @info "$i of $Nsamples done"
        end

        orb0 = KeplerianElements(
            date_to_jd(2023, 1, 1, 0, 0, 0),
            inp.a,
            inp.e,
            1.5 |> deg2rad,
            1.5    |> deg2rad,
            1.5     |> deg2rad,
            inp.f     |> deg2rad
        )
        if orb0.a * (1 - orb0.e) < EARTH_EQUATORIAL_RADIUS
            continue
        end
        deltaV = [inp.dvx, inp.dvy, inp.dvz]
    
        r_final, total_time, orb_postman, orb_final, dir_v_final = final_position(
            orb0, 
            deltaV, 
            inp.fT1, 
            inp.fT2
        )
    
        if orb_final.a * (1 - orb_final.e) < EARTH_EQUATORIAL_RADIUS
            continue
        end
    
        model, r_maneuver, v_post_maneuver, rf, vf = single_maneuver_model(
            orb0, 
            r_final, 
            dir_v_final,
            total_time);
        set_silent(model)
    
        stats = nothing
        output = nothing
        try
            stats = @timed optimize!(model)    
        catch e
            output = OutputParams(
                inp,
                e,
                nothing,
                nothing,
                nothing,
                nothing
            )
        else
            output = OutputParams(
                inp,
                nothing,
                is_solved_and_feasible(model),
                stats.time,
                value(model[:ΔVmag]) * value.(model[:ΔVdir]),
                value(model[:Δt_maneuver])
            )
        end

        push!(outputs, output)
        println(file, output)
        flush(file)
    end
end

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

deltaV = [1000, 0, 1000]

r_final, total_time, orb_postman, orb_final, dir_v_final = final_position(
    orb0, 
    deltaV, 
    0.5, 
    0.7
)

plot_orbit(orb0, orb_postman, orb_final)
## solving maneuver
model, r_maneuver, v_post_maneuver, rf, vf = single_maneuver_model(
    orb0, 
    r_final, 
    dir_v_final,
    total_time);
# set_silent(model)
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
