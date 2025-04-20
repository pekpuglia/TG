using SatelliteToolboxBase
using SatelliteToolboxPropagators
using JuMP
using Ipopt
using Setfield
using ForwardDiff
include("TG.jl")
using .TG
using LinearAlgebra
using Random
using Printf
##
a_list = 7000e3:100e3:10000e3
e_list = [0.0, 0.001, 0.01, 0.02, 0.1, 0.2, 0.4, 0.7, 0.9]
f_list = deg2rad.([0, 0.1, 1, (10:10:360)...])
i_list = deg2rad.(0:10:180.0)
Ω_list = deg2rad.(0:10:360.0)
ω_list = deg2rad.(0:10:360.0)
t_list = 0:0.1:3
##
solution_types = ["start constraint", "end constraint"]
##
function output_format(sol_type, exec_time, deltat, givenorb::KeplerianElements, start_orb_params::FullOrbitalParameters, end_orb_params::FullOrbitalParameters; error = nothing, converged)
    given_a = givenorb.a
    given_e = givenorb.e
    given_i = givenorb.i
    given_Ω = givenorb.Ω
    given_ω = givenorb.ω
    given_f = givenorb.f
    given_r, given_v = kepler_to_rv(givenorb)
    
    start_a = value(start_orb_params.a)
    start_e = value(start_orb_params.e)
    start_i = value(start_orb_params.i)
    start_Ω = value(start_orb_params.Ω)
    start_ω = value(start_orb_params.ω)
    start_f = value(start_orb_params.nu)
    start_r = value.(start_orb_params.r)
    start_v = value.(start_orb_params.v)
    
    end_a = value(end_orb_params.a)
    end_e = value(end_orb_params.e)
    end_i = value(end_orb_params.i)
    end_Ω = value(end_orb_params.Ω)
    end_ω = value(end_orb_params.ω)
    end_f = value(end_orb_params.nu)
    end_r = value.(end_orb_params.r)
    end_v = value.(end_orb_params.v)
    
    "Solution type: $sol_type\nSolution time: $exec_time\n"*something(error,"")*"Solved and feasible: $converged\nΔt=$deltat\n"*
    (
        @sprintf "\t       Given orbital elements: a = %5.4e, e = %5.4f, f = %5.4f, i = %5.4f, Ω = %5.4f, ω = %5.4f\n" given_a given_e given_f mod(given_i, 2pi) mod(given_Ω, 2pi) mod(given_ω, 2pi)
    )*(
        @sprintf "\tSolved start orbital elements: a = %5.4e, e = %5.4f, f = %5.4f, i = %5.4f, Ω = %5.4f, ω = %5.4f\n" start_a     start_e     start_f     mod(start_i, 2pi)     mod(start_Ω, 2pi)     mod(start_ω, 2pi)
    )*(
        @sprintf "\t  Solved end orbital elements: a = %5.4e, e = %5.4f, f = %5.4f, i = %5.4f, Ω = %5.4f, ω = %5.4f\n" end_a     end_e     end_f     mod(end_i, 2pi)     mod(end_Ω, 2pi)     mod(end_ω, 2pi)
    )*
    "\t  Given state vector: r = $given_r, v = $given_v\n"*
    "\tInitial state vector: r = $(start_r), v = $(start_v)\n\n"*
    "\t  Final state vector: r = $(end_r), v = $(end_v)\n\n"
end
##
Nsamples = 100

outfile = "./src/LEO_coast_out"
open(outfile, "a") do file
    for i in 1:Nsamples
        orb = KeplerianElements(
            date_to_jd(2023, 1, 1, 0, 0, 0),
            rand(a_list),
            rand(e_list),
            rand(i_list),
            rand(Ω_list),
            rand(ω_list),
            rand(f_list),
        )
        while true
            if orb.a*(1-orb.e) >= EARTH_EQUATORIAL_RADIUS
                break
            end
            orb = @set orb.a = rand(a_list)
            orb = @set orb.e = rand(e_list)
        end

        given_r, given_v = kepler_to_rv(orb)

        model = Model(
            optimizer_with_attributes(Ipopt.Optimizer,
            "max_wall_time" => 30.0)
        )
        set_silent(model)

        orbparams_i = add_orbital_elements!(model)
        r0, v0, a0, e0, i0, Ω0, ω0, nu0, M0, E0 = getfield.(Ref(orbparams_i), fieldnames(FullOrbitalParameters))

        orbparams_f = add_orbital_elements!(model)
        rf, vf, af, ef, i_f, Ωf, ωf, nuf, Mf, Ef = getfield.(Ref(orbparams_f), fieldnames(FullOrbitalParameters))

        @constraint(model, af == a0)
        @constraint(model, ef == e0)
        @constraint(model, i_f == i0)
        @constraint(model, Ωf == Ω0)
        @constraint(model, ωf == ω0)
        
        T = orbital_period(orbparams_i.a, GM_EARTH)
        selected_t = rand(t_list)
        Δt = selected_t * orbital_period(orb, GM_EARTH)

        @constraint(model, Δt == (Mf - M0) / (2π) * T)

        sol_type = rand(solution_types)
        if sol_type == solution_types[1]
            @constraint(model, rf .== given_r)
            @constraint(model, vf .== given_v)
        else
            @constraint(model, r0 .== given_r)
            @constraint(model, v0 .== given_v)
        end

        error = nothing
        stats = nothing
        try
            stats = @timed optimize!(model)
        catch e
            error = e
        end

        output = output_format(sol_type, stats.time, selected_t, orb, orbparams_i, orbparams_f; error=error, converged=is_solved_and_feasible(model))

        println(file, output)
        flush(file)
    end
end