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
a_list = 7000e3:1000e3:40000e3
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
        @sprintf "\tGiven orbital elements : a = %5.4e, e = %5.4f, f = %5.4f, i = %5.4f, Ω = %5.4f, ω = %5.4f\n" orb.a orb.e mod(orb.f, 2pi) mod(orb.i, 2pi) mod(orb.Ω, 2pi) mod(orb.ω, 2pi)
    )*(
        @sprintf "\tSolved orbital elements: a = %5.4e, e = %5.4f, f = %5.4f, i = %5.4f, Ω = %5.4f, ω = %5.4f\n" a     e     mod(f, 2pi)     mod(i, 2pi)     mod(Ω, 2pi)     mod(ω, 2pi)
    )*
    "\t Given state vector: r = $givenr, v = $givenv\n"*
    "\tSolved state vector: r = $(r), v = $(v)\n\n"
end
##
function add_orbital_elements_fix!(model)
    Vorb_sup = √(GM_EARTH/EARTH_EQUATORIAL_RADIUS)
    rscaled = @variable(model, [1:3])
    @constraint(model, rscaled' * rscaled >= 1)
    r = EARTH_EQUATORIAL_RADIUS * rscaled
    vscaled = @variable(model, [1:3])
    v = Vorb_sup*vscaled

    a = @variable(model, lower_bound = EARTH_EQUATORIAL_RADIUS, start = 2EARTH_EQUATORIAL_RADIUS)
    e = @variable(model, lower_bound = 0, upper_bound = 1) 
    @constraint(model, a*(1-e) >= EARTH_EQUATORIAL_RADIUS)
    i = @variable(model, lower_bound = 0, upper_bound = π, base_name = "i")
    Ω = @variable(model, base_name = "Ω")
    ω = @variable(model, base_name = "ω")
    nu = @variable(model, lower_bound = -2π, upper_bound = 2π, base_name = "nu")

    #rad!!!
    M = @variable(model, lower_bound = 0.0, base_name = "M")
    E = @variable(model, lower_bound = 0.0, base_name= "E")
    
    R3Omega = [
         cos(Ω) sin(Ω) 0
        -sin(Ω) cos(Ω) 0
        0          0        1
    ]

    R1i = [
        1  0         0
        0  cos(i) sin(i)
        0 -sin(i) cos(i)
    ]

    R3omega = [
        cos(ω)  sin(ω) 0
        -sin(ω) cos(ω) 0
        0          0        1
    ]

    QXxbar = R3omega * R1i * R3Omega

    #curtis chap 4
    #h^2/mu = p = a (1-e^2)
    r_perifocal = a*(1-e^2) * 1/(1+e*cos(nu)) * [cos(nu); sin(nu); 0]
    
    h = √(GM_EARTH*a*(1-e^2))
    v_perifocal = GM_EARTH / h * [-sin(nu); e + cos(nu); 0]


    @constraint(model, r .== QXxbar' * r_perifocal)
    @constraint(model, v .== QXxbar' * v_perifocal)
    
    @constraint(model, E - e*sin(E) == M)

    #curtis page 144 & 145
    @constraint(model, nu == 2 * atan(√(1+e)*sin(E/2), √(1-e)*cos(E/2)))

    FullOrbitalParameters(r, v, a, e, i, Ω, ω, nu, M, E)
end
##
Nsamples = 100

outfile = "./src/elements/test_elements_out"
open(outfile, "w") do file
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

        orbparams = add_orbital_elements_fix!(model)
        r, v, a, e, i, Ω, ω, f, M, E = getfield.(Ref(orbparams), fieldnames(FullOrbitalParameters))
        
        sol_type = "rv input"
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