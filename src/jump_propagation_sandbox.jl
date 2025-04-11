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
    91.5    |> deg2rad,
    1.5     |> deg2rad,
    1.5     |> deg2rad
)
r0, v0 = kepler_to_rv(orb0)
T0 = orbital_period(orb0, GM_EARTH)
# plot_orbit(orb0)
## CURTIS chap 3

model = Model(
    optimizer_with_attributes(Ipopt.Optimizer,
        "max_iter" => 10000,
        "max_wall_time" => 30.0)
)

@variables(model, begin
    a
    e
    0 <= id <= 180
    0 <= Ω <= 360
    0 <= ω <= 360
    nu
end)
set_start_value(Ω, 180)
set_start_value(ω, 180)

# @variable(model, r[1:3])
# @variable(model, v[1:3])

rnorm = norm(r0)
vnorm = norm(v0)
vr = dot(r0/rnorm, v0)

h = cross(r0, v0)

@constraint(model, cosd(id) == h[3]/√(h' * h))

N = cross([0;0;1], h)

Nnorm = √(N' * N)

@constraint(model, cosd(Ω) == N[1]/Nnorm)
@constraint(model, sind(Ω) == N[2]/Nnorm)

exc_vec = (vnorm^2 / GM_EARTH - 1 / rnorm) * r0 - rnorm*vr/GM_EARTH .* v0

e = √(exc_vec' * exc_vec)

@constraint(model, cosd(ω) == dot(N, exc_vec) / (Nnorm*e))

N_e_cross = cross(N, exc_vec)
normal_N_e_cross = dot(N_e_cross, h/√(h' * h))

@constraint(model, sind(ω) == normal_N_e_cross / (Nnorm*e))

model

##
optimize!(model)
value(model[:ω])
##
value.(all_variables(model))
##
value(normal_N_e_cross)
