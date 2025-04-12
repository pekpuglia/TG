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
  a_list = 7000e3:100e3:8000e3
  e_list = [0.001, 0.01, 0.02, 0.03, 0.05]
  f_list = [0.1, 1, 180, 50]
dvx_list = [-1100, -1000, -500, -200, 0, 200, 500, 1000, 1100]
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
Nsamples = 100
random_inputs = input_params_list[1:Nsamples]

outputs = []

out_file = "./src/LEO_single_impulse_out"

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
