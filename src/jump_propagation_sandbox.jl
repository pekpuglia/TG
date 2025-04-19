using SatelliteToolboxBase
using SatelliteToolboxPropagators
using JuMP
using Ipopt
using Setfield
using ForwardDiff
include("TG.jl")
using .TG
using LinearAlgebra
## example
rp = (6378+400)*1000.0
ra = (6378+4000)*1000.0
agiven = (rp + ra) / 2
egiven = (ra - rp) / (ra + rp)
orb0 = KeplerianElements(
    date_to_jd(2023, 1, 1, 0, 0, 0),
    agiven,
    egiven,
    30 |> deg2rad,
    15    |> deg2rad,
    60     |> deg2rad,
    180     |> deg2rad
)
r0, v0 = kepler_to_rv(orb0)
##
function add_orbital_elements_fix!(model, given_rv = true)
    Vorb_sup = √(GM_EARTH/EARTH_EQUATORIAL_RADIUS)
    rscaled = @variable(model, [1:3], start = 1.0)
    r = EARTH_EQUATORIAL_RADIUS * rscaled
    vscaled = @variable(model, [1:3])
    set_start_value(vscaled[1], 1.0)
    v = Vorb_sup*vscaled

    a = @variable(model, lower_bound = 0, start = EARTH_EQUATORIAL_RADIUS)
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

orbparams_i = add_orbital_elements_fix!(model)
r, v, a, e, i, Ω, ω, nu, M, E = getfield.(Ref(orbparams_i), fieldnames(FullOrbitalParameters))

@constraint(model, r .== r0)
@constraint(model, v .== v0)


model
##
optimize!(model)
##
#ok!
value(a), agiven
##
#ok!
value(e), egiven
##
#ok!
rad2deg(value(i))
##
#oK!
rad2deg(value(Ω))
##
#ok!
rad2deg(value(ω))
##
rad2deg(value(nu))
##
solved_r = value.(r)
solved_v = value.(v)

plot_orbit(
    orb0,
    rv_to_kepler(solved_r, solved_v),
)