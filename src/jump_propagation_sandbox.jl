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
agiven = (rp + ra) / 2
egiven = (ra - rp) / (ra + rp)
orbi = KeplerianElements(
    date_to_jd(2023, 1, 1, 0, 0, 0),
    agiven,
    egiven,
    30 |> deg2rad,
    0    |> deg2rad,
    0     |> deg2rad,
    0     |> deg2rad
)
ri, vi = kepler_to_rv(orbi)
T = orbital_period(orbi, GM_EARTH)
Δt = 1.2T
##
model = Model(
    optimizer_with_attributes(Ipopt.Optimizer,
    "max_wall_time" => 30.0)
)
orbparams_i = add_orbital_elements!(model)
r0, v0, a0, e0, i0, Ω0, ω0, nu0, M0, E0 = getfield.(Ref(orbparams_i), fieldnames(FullOrbitalParameters))

orbparams_f = add_orbital_elements!(model)
rf, vf, af, ef, i_f, Ωf, ωf, nuf, Mf, Ef = getfield.(Ref(orbparams_f), fieldnames(FullOrbitalParameters))

@constraint(model, rf .== ri)
@constraint(model, vf .== vi)

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
#ok!
value(af), agiven
##
#ok!
value(ef), egiven
##
#ok!
rad2deg(value(i_f))
##
#oK!
rad2deg(value(Ωf))
##
#ok!
rad2deg(value(ωf))
##
rad2deg(value(nuf))
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