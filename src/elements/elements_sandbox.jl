using SatelliteToolboxBase
using SatelliteToolboxPropagators
using JuMP
using Ipopt
using Setfield
using ForwardDiff
include("../TG.jl")
using .TG
using LinearAlgebra
##
orb = KeplerianElements(
    date_to_jd(2023, 1, 1, 0, 0, 0),
    7.0000e+06,
    0.000,
    0.6981,
    4.0143,
    0.6981,
    4.1888
)
given_r, given_v = kepler_to_rv(orb)
plot_orbit(orb)
##
function add_orbital_elements_fix!(model)
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
    r_perifocal = a*(1-e^2) * 1/(1+e*cos(nu)) * [cos(nu); sin(nu); 0]
    
    h = √(GM_EARTH*a*(1-e^2))
    v_perifocal = GM_EARTH / h * [-sin(nu); e + cos(nu); 0]


    @constraint(model, rscaled .== QXxbar' * r_perifocal / EARTH_EQUATORIAL_RADIUS)
    @constraint(model, vscaled .== QXxbar' * v_perifocal / Vorb_sup)
    
    @constraint(model, E - e*sin(E) == M)

    #curtis page 144 & 145
    @constraint(model, nu == 2 * atan(√(1+e)*sin(E/2), √(1-e)*cos(E/2)))

    FullOrbitalParameters(r, v, a, e, i, Ω, ω, nu, M, E), rscaled, vscaled
end
##
model = Model(
    optimizer_with_attributes(Ipopt.Optimizer,
    "max_wall_time" => 30.0)
)
orbparams, rsc, vsc = add_orbital_elements_fix!(model)
r, v, a, e, i, Ω, ω, f, M, E = getfield.(Ref(orbparams), fieldnames(FullOrbitalParameters))
        
@constraint(model, rsc .== (given_r/EARTH_EQUATORIAL_RADIUS))
Vorb_sup = √(GM_EARTH/EARTH_EQUATORIAL_RADIUS)

@constraint(model, vsc .== (given_v/Vorb_sup))

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
value(f)
##
solved_r = value.(r)
solved_v = value.(v)

plot_orbit(
    orb,
    KeplerianElements(
        date_to_jd(2023, 1, 1, 0, 0, 0),
        value(a),
        clamp(value(e), 0, 1),
        value(i),
        value(Ω),
        value(ω),
        value(f)
    ),
    rv_to_kepler(solved_r, solved_v),
)