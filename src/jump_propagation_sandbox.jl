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
#designing single maneuver inversely
orb0 = KeplerianElements(
    date_to_jd(2023, 1, 1, 0, 0, 0),
    7000e3,
    0.01,
    1.5 |> deg2rad,
    1.5    |> deg2rad,
    1.5     |> deg2rad,
    1.5     |> deg2rad
)
r0, v0 = kepler_to_rv(orb0)
T0 = orbital_period(orb0, GM_EARTH)
## auxiliary parameters

model = Model(
    optimizer_with_attributes(Ipopt.Optimizer,
        "max_iter" => 10000,
        "max_wall_time" => 30.0)
)

@variables(model, begin
    a0
    e0
    0 <= i0d <= 180
    Ω0
    ω0
    M0
    M
    E
    nu0
    nu
end)

@variable(model, r[1:3])
@variable(model, v[1:3])

r0norm = norm(r0)
v0norm = norm(v0)
v0r = dot(r0/r0norm, v0)

h0 = cross(r0, v0)

@constraint(model, cosd(i0d) == h0[3]/√(h0' * h0))

model

##
optimize!(model)
##
value.(all_variables(model))
##
