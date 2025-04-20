using SatelliteToolboxBase
using SatelliteToolboxPropagators
using JuMP
using Ipopt
using Setfield
using ForwardDiff
include("../TG.jl")
using .TG
using LinearAlgebra
## example
orb = KeplerianElements(
    date_to_jd(2023, 1, 1, 0, 0, 0),
    9.6000e+06,
    0.0100,
    0.6981,
    0,
    1.2217,
    5.4105
)
given_r, given_v = kepler_to_rv(orb)
plot_orbit(orb)
##
function add_orbital_elements_fix!(model, given_rv = false)
    Vorb_sup = √(GM_EARTH/EARTH_EQUATORIAL_RADIUS)
    rscaled = @variable(model, [1:3], start = 1.0)
    r = EARTH_EQUATORIAL_RADIUS * rscaled
    vscaled = @variable(model, [1:3])
    set_start_value(vscaled[1], 1.0)
    v = Vorb_sup*vscaled

    a = @variable(model, lower_bound = EARTH_EQUATORIAL_RADIUS, start = 2EARTH_EQUATORIAL_RADIUS)
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

    if given_rv
        hvec = cross(r, v)
        h = √(hvec' * hvec)
    else
        h = √(GM_EARTH*a*(1-e^2))
    end

    #curtis chap 4
    #h^2/mu = p = a (1-e^2)
    r_perifocal = a*(1-e^2) * 1/(1+e*cos(nu)) * [cos(nu); sin(nu); 0]

    v_perifocal = GM_EARTH / h * [-sin(nu); e + cos(nu); 0]


    @constraint(model, r .== QXxbar' * r_perifocal)
    @constraint(model, v .== QXxbar' * v_perifocal)
    
    @constraint(model, E - e*sin(E) == M)

    #curtis page 144 & 145
    @constraint(model, nu == 2 * atan(√(1+e)*sin(E/2), √(1-e)*cos(E/2)))

    FullOrbitalParameters(r, v, a, e, i, Ω, ω, nu, M, E)
end
##
model = Model(
    optimizer_with_attributes(Ipopt.Optimizer,
    "max_wall_time" => 30.0)
)
orbparams = add_orbital_elements_fix!(model)
r, v, a, e, i, Ω, ω, f, M, E = getfield.(Ref(orbparams), fieldnames(FullOrbitalParameters))
        
@constraint(model, r .== given_r)
@constraint(model, v .== given_v)
# set_start_value(a, orb.a)
# @constraint(model, a == orb.a)
# @constraint(model, e == orb.e)
# @constraint(model, i == orb.i)
# @constraint(model, Ω == orb.Ω)
# @constraint(model, ω == orb.ω)
# @constraint(model, f == orb.f)
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
    rv_to_kepler(solved_r, solved_v),
)