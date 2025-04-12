using SatelliteToolboxBase
using SatelliteToolboxPropagators
using JuMP
using Ipopt
using Setfield
using ForwardDiff
include("TG.jl")
using .TG
using LinearAlgebra
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


orbparams_i = add_orbital_elements!(model)
orbparams_f = add_orbital_elements!(model)


@constraint(model, orbparams_i.r .== r0)
@constraint(model, orbparams_i.v .== v0)
@constraint(model, orbparams_f.r .== r)
@constraint(model, orbparams_f.v .== v)

@variable(model, Δt)

@constraint(model, Δt == (orbparams_f.M - orbparams_i.M) / (2π) * T)

model
##
optimize!(model)
value.(all_variables(model))
##
value(orbparams_i.e), value(orbparams_f.e)
##
value(orbparams_i.a), value(orbparams_f.a)
##
value(orbparams_i.E), value(orbparams_f.E)
##
value(orbparams_i.M), value(orbparams_f.M)
##
value(model[:Δt])