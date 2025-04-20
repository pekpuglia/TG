using SatelliteToolboxBase
using SatelliteToolboxPropagators
using JuMP
using Ipopt
using Setfield
using ForwardDiff
include("../TG.jl")
using .TG
using LinearAlgebra
using Random
using Printf
##
a_list = 7000e3:100e3:40000e3
e_list = [0.0, 0.001, 0.01, 0.02, 0.1, 0.2, 0.4, 0.7, 0.9]
f_list = deg2rad.([0, 0.1, 1, (10:10:360)...])
i_list = deg2rad.(0:10:180.0)
Ω_list = deg2rad.(0:10:360.0)
ω_list = deg2rad.(0:10:360.0)
##
solution_types = ["rv input", "elements input"]
##
function output_format(sol_type, time, orb::KeplerianElements, givenr, givenv, orbparams::FullOrbitalParameters; error = nothing, converged)
    a = value(orbparams.a)
    e = value(orbparams.e)
    i = value(orbparams.i)
    Ω = value(orbparams.Ω)
    ω = value(orbparams.ω)
    f = value(orbparams.nu)
    r = value.(orbparams.r)
    v = value.(orbparams.v)
    "Solution type: $sol_type\nSolution time: $time\n"*something(error,"")*"Solved and feasible: $converged\n"*
    (
        @sprintf "\tGiven orbital elements : a = %5.4e, e = %5.4f, f = %5.4f, i = %5.4f, Ω = %5.4f, ω = %5.4f\n" orb.a orb.e orb.f mod(orb.i, 2pi) mod(orb.Ω, 2pi) mod(orb.ω, 2pi)
    )*(
        @sprintf "\tSolved orbital elements: a = %5.4e, e = %5.4f, f = %5.4f, i = %5.4f, Ω = %5.4f, ω = %5.4f\n" a     e     f     mod(i, 2pi)     mod(Ω, 2pi)     mod(ω, 2pi)
    )*
    "\t Given state vector: r = $givenr, v = $givenv\n"*
    "\tSolved state vector: r = $(r), v = $(v)\n\n"
end
##
Nsamples = 100

outfile = "./src/elements/elements_out"
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

        orbparams = add_orbital_elements!(model)
        r, v, a, e, i, Ω, ω, f, M, E = getfield.(Ref(orbparams), fieldnames(FullOrbitalParameters))
        
        sol_type = rand(solution_types)
        if sol_type == "rv input"
            @constraint(model, r .== given_r)
            @constraint(model, v .== given_v)
        else
            @constraint(model, a == orb.a)
            @constraint(model, e == orb.e)
            @constraint(model, i == orb.i)
            @constraint(model, Ω == orb.Ω)
            @constraint(model, ω == orb.ω)
            @constraint(model, f == orb.f)
        end

        error = nothing
        stats = nothing
        try
            stats = @timed optimize!(model)
        catch e
            error = e
        end

        output = output_format(sol_type, stats.time, orb, given_r, given_v, orbparams; error=error, converged=is_solved_and_feasible(model))

        println(file, output)
        flush(file)
    end
end