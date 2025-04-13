using SatelliteToolboxBase
using SatelliteToolboxPropagators
using JuMP
using Ipopt
using Setfield
using ForwardDiff
include("TG.jl")
using .TG
using LinearAlgebra
## example 2.7
rp = (6378+400)*1000.0
ra = (6378+4000)*1000.0
agiven = (rp + ra) / 2
egiven = (ra - rp) / (ra + rp)
orb0 = KeplerianElements(
    date_to_jd(2023, 1, 1, 0, 0, 0),
    agiven,
    egiven,
    30 |> deg2rad,
    0    |> deg2rad,
    0     |> deg2rad,
    0     |> deg2rad
)
r0, v0 = kepler_to_rv(orb0)
plot_orbit(orb0)
##
function add_orbital_elements_fix!(model)
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
    nu = @variable(model, base_name = "nu")

    #rad!!!
    M = @variable(model, base_name = "M")
    E = @variable(model, base_name= "E")
    
    hvec = cross(r, v)
    h = √(hvec' * hvec)

    #curtis chap 4
    r_perifocal = h^2/GM_EARTH * 1/(1+e*cos(nu)) * [cos(nu); sin(nu); 0]

    v_perifocal = GM_EARTH / h * [-sin(nu); e + cos(nu); 0]

    R3Omega = [
        cos(Ω) -sin(Ω) 0
        sin(Ω)  cos(Ω) 0
        0           0         1
    ]

    R1i = [
        1 0          0
        0 cos(i) -sin(i)
        0 sin(i)  cos(i)
    ]

    R3omega = [
        cos(Ω) -sin(Ω) 0
        sin(Ω)  cos(Ω) 0
        0           0         1
    ]

    QXxbar = R3omega * R1i * R3Omega

    @constraint(model, r .== QXxbar' * r_perifocal)
    @constraint(model, v .== QXxbar' * v_perifocal)
    
    @constraint(model, E - e*sin(E) == M)

    #curtis page 144 & 145
    @constraint(model, E == 2 * atan(√((1-e)/(1+e))*tan(nu/2)))

    FullOrbitalParameters(r, v, a, e, i, Ω, ω, nu, M, E)
end
##
model = Model(
    optimizer_with_attributes(Ipopt.Optimizer,
    "max_wall_time" => 30.0)
)
orbparams_i = add_orbital_elements_fix!(model)
r, v, a, e, i, Ω, ω, nu, M, E = getfield.(Ref(orbparams_i), fieldnames(FullOrbitalParameters))

# @constraint(model, a*(1-e) == rp)
# @constraint(model, a*(1+e) == ra)

# @constraint(model, i == deg2rad(30))

# @constraint(model, nu == 0)

# @constraint(model, Ω == 0)

# @constraint(model, ω == 0)

@constraint(model, r .== r0)
@constraint(model, v .== v0)

model
##
optimize!(model)
##
value(a), agiven
##
value(e), egiven
##
value(nu)
##
rad2deg(value(i))
##
value(Ω)
##
value(ω)
##
r0, value.(r)
##
v0, value.(v)
##
solved_r = value.(r)
solved_v = value.(v)
plot_orbit(
    orb0,
    rv_to_kepler(solved_r, solved_v),
)