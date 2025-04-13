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
Δt = 3*3600.0
##
model = Model(
    optimizer_with_attributes(Ipopt.Optimizer,
    "max_wall_time" => 30.0)
)
orbparams_i = add_orbital_elements!(model, true)
r0, v0, a0, e0, i0, Ω0, ω0, nu0, M0, E0 = getfield.(Ref(orbparams_i), fieldnames(FullOrbitalParameters))

orbparams_f = add_orbital_elements!(model, false)
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