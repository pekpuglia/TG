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
orbi = KeplerianElements(
    date_to_jd(2023, 1, 1, 0, 0, 0),
    1.0000e+07,
    0.0100,
    0.1,
    4.5379,
    3.4907,
    2.6180
)
ri, vi = kepler_to_rv(orbi)
T = orbital_period(orbi, GM_EARTH)
Δt = 1.4T
plot_orbit(orbi)
##
model = Model(
    optimizer_with_attributes(Ipopt.Optimizer,
    "max_wall_time" => 30.0)
)
orbparams_i = add_orbital_elements!(model)
r0, v0, a0, e0, i0, Ω0, ω0, nu0, M0, E0 = getfield.(Ref(orbparams_i), fieldnames(FullOrbitalParameters))

orbparams_f = add_orbital_elements!(model)
rf, vf, af, ef, i_f, Ωf, ωf, nuf, Mf, Ef = getfield.(Ref(orbparams_f), fieldnames(FullOrbitalParameters))

@constraint(model, r0 .== ri)
@constraint(model, v0 .== vi)

@constraint(model, af == a0)
@constraint(model, ef == e0)
@constraint(model, i_f == i0)
@constraint(model, Ωf == Ω0)
@constraint(model, ωf == ω0)

@constraint(model, Δt == (Mf - M0) / (2π) * T)
model
##
optimize!(model)
##
value(a0), value(af)
##
value(e0), value(ef)
##
value(i0), value(i_f)
##
#oK!
value(Ω0), value(Ωf)
##
value(ω0), value(ωf)
##
value(nu0), value(nuf)
##
solved_rf = value.(rf)
solved_vf = value.(vf)

solved_r0 = value.(r0)
solved_v0 = value.(v0)


plot_orbit(
    orbi,
    rv_to_kepler(solved_rf, solved_vf),
    rv_to_kepler(solved_r0, solved_v0),
)