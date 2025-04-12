using SatelliteToolboxBase
using SatelliteToolboxPropagators
using JuMP
using Ipopt
using Setfield
using ForwardDiff
include("TG.jl")
using .TG
using LinearAlgebra
##
##
## obj: dar r0, v0, deltaT, receber r, v
#example 3.1 curtis
rp = 9600e3
ra = 21000e3
a = (rp + ra) / 2
e = (ra - rp) / (ra + rp)
orb0 = KeplerianElements(
    date_to_jd(2023, 1, 1, 0, 0, 0),
    a,
    e,
    30 |> deg2rad,
    0    |> deg2rad,
    0     |> deg2rad,
    0     |> deg2rad
)
r0, v0 = kepler_to_rv(orb0)
orbf = @set orb0.f = 120 |> deg2rad
r, v = kepler_to_rv(orbf)
T = orbital_period(orb0, GM_EARTH)
plot_orbit(orb0, orbf)
## CURTIS chap 3
model = Model(
    optimizer_with_attributes(Ipopt.Optimizer,
        "max_iter" => 10000,
        "max_wall_time" => 30.0)
)


ri, vi, ai, ei, ii, Ωi, ωi, nui, Mi, Ei = add_orbital_elements!(model)
rf, vf, af, ef, i_f, Ωf, ωf, nuf, Mf, Ef = add_orbital_elements!(model)

@constraint(model, ri .== r0)
@constraint(model, vi .== v0)
@constraint(model, rf .== r)
@constraint(model, vf .== v)

@variable(model, Δt)

@constraint(model, Δt == (Mf - Mi) / (2π) * T)

model
##
optimize!(model)
value.(all_variables(model))
##
value(ei), value(ef)
##
value(ai), value(af)
##
value(Ei), value(Ef)
##
value(Mi), value(Mf)
##
value(model[:Δt])