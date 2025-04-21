using SatelliteToolboxBase
using SatelliteToolboxPropagators
using JuMP
using Ipopt
using Setfield
using ForwardDiff
include("../TG.jl")
using .TG
using LinearAlgebra
using Printf
##
orb = KeplerianElements(
    date_to_jd(2023, 1, 1, 0, 0, 0),
    9.7000e+06,
    0.0010,
    0.1000,
    4.1888,
    1.3963,
    4.7124
)
given_r, given_v = kepler_to_rv(orb)
plot_orbit(orb)
##
# function add_orbital_elements_fix!(model)
model = Model(
    optimizer_with_attributes(Ipopt.Optimizer,
    "max_wall_time" => 30.0)
)

Vorb_sup = √(GM_EARTH/EARTH_EQUATORIAL_RADIUS)
rscaled = @variable(model, [1:3])
@constraint(model, rscaled' * rscaled >= 1)
r = EARTH_EQUATORIAL_RADIUS * rscaled
vscaled = @variable(model, [1:3])
v = Vorb_sup*vscaled

ascaled = @variable(model, lower_bound = 1.0)
a = EARTH_EQUATORIAL_RADIUS * ascaled
e = @variable(model, lower_bound = 0, upper_bound = 1) 
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
r_perifocal_scaled = ascaled*(1-e^2) * 1/(1+e*cos(nu)) * [cos(nu); sin(nu); 0]

h = √(GM_EARTH*a*(1-e^2))
v_perifocal = GM_EARTH / h * [-sin(nu); e + cos(nu); 0]


@constraint(model, rscaled .== QXxbar' * r_perifocal_scaled)
@constraint(model, v .== QXxbar' * v_perifocal)

@constraint(model, E - e*sin(E) == M)

#curtis page 144 & 145
@constraint(model, nu == 2 * atan(√(1+e)*sin(E/2), √(1-e)*cos(E/2)))

@constraint(model, rscaled .== (given_r/EARTH_EQUATORIAL_RADIUS))
Vorb_sup = √(GM_EARTH/EARTH_EQUATORIAL_RADIUS)

@constraint(model, vscaled .== (given_v/Vorb_sup))

model
##
optimize!(model)
##
value(a)
##
value(e)
##
value(i)
##
#oK!
value(Ω)
##
value(ω)
##
value(nu)
##
solved_r = value.(rscaled)*EARTH_EQUATORIAL_RADIUS
solved_v = value.(vscaled)*Vorb_sup

plot_orbit(
    orb,
    KeplerianElements(
        date_to_jd(2023, 1, 1, 0, 0, 0),
        value(a),
        clamp(value(e), 0, 1),
        value(i),
        value(Ω),
        value(ω),
        value(nu)
    ),
    rv_to_kepler(solved_r, solved_v),
)
## ##############################################################################
# monte carlo
a_list = 7000e3:100e3:10000e3
e_list = [0.0, 0.001, 0.01, 0.02, 0.1, 0.2, 0.4, 0.7, 0.9]
f_list = deg2rad.([0, 0.1, 1, (10:10:360)...])
i_list = deg2rad.(0:10:80.0)
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

        rscaled, vscaled, a, e, i, Ω, ω, nu, M, E = add_orbital_elements_fix!(model)                
        
        Vorb_sup = √(GM_EARTH/EARTH_EQUATORIAL_RADIUS)
        orbparams = FullOrbitalParameters(
            rscaled*EARTH_EQUATORIAL_RADIUS,
            vscaled*Vorb_sup,
            a, e, i, Ω, ω, nu, M, E
        )
        
        sol_type = "rv input"
        if sol_type == "rv input"
            @constraint(model, rscaled .== (given_r/EARTH_EQUATORIAL_RADIUS))
            
            @constraint(model, vscaled .== (given_v/Vorb_sup))
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