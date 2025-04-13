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
rp = 9600e3
ra = 21000e3
a = (rp + ra) / 2
e = (ra - rp) / (ra + rp)
orbi = KeplerianElements(
    date_to_jd(2023, 1, 1, 0, 0, 0),
    a,
    e,
    30 |> deg2rad,
    0    |> deg2rad,
    0     |> deg2rad,
    0     |> deg2rad
)
ri, vi = kepler_to_rv(orbi)
Δt = 3600.0
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
    nu = @variable(model, base_name = "nu")

    #rad!!!
    M = @variable(model, base_name = "M")
    E = @variable(model, base_name= "E")
    
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
    @constraint(model, nu == 2 * atan(√((1+e)/(1-e))*tan(E/2)))

    FullOrbitalParameters(r, v, a, e, i, Ω, ω, nu, M, E)
end
##
model = Model(
    optimizer_with_attributes(Ipopt.Optimizer,
    "max_wall_time" => 30.0)
)
orbparams_i = add_orbital_elements_fix!(model, true)
r0, v0, a0, e0, i0, Ω0, ω0, nu0, M0, E0 = getfield.(Ref(orbparams_i), fieldnames(FullOrbitalParameters))

orbparams_f = add_orbital_elements_fix!(model, false)
rf, vf, af, ef, i_f, Ωf, ωf, nuf, Mf, Ef = getfield.(Ref(orbparams_f), fieldnames(FullOrbitalParameters))

@constraint(model, r0 .== ri)
@constraint(model, v0 .== vi)

@constraint(model, af == a0)
@constraint(model, ef == e0)
@constraint(model, i_f == i0)
@constraint(model, Ωf == Ω0)
@constraint(model, ωf == ω0)

T = orbital_period(orbi, GM_EARTH)

@constraint(model, Δt == (Mf - M0) / (2π) * T)


model
##
optimize!(model)
##
value(a0), value(af), a
##
value(e0), value(ef), e
##
value(nu0), value(nuf)
##
value(i0), value(i_f)
##
value(Ω0), value(Ωf)
##
value(ω0), value(ωf)
##
value(E0), value(Ef)
##
(value(Mf) - value(M0)) / (2π) * T
#Ef = 3.4794
#nu = 193.2deg
##
solved_r0 = value.(r0)
solved_v0 = value.(v0)

solved_rf = value.(rf)
solved_vf = value.(vf)

plot_orbit(
    orbi,
    rv_to_kepler(solved_r0, solved_v0),
    rv_to_kepler(solved_rf, solved_vf),
)