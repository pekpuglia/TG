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
    91.5 |> deg2rad,
    91.5    |> deg2rad,
    271.5     |> deg2rad,
    263     |> deg2rad
)
r0, v0 = kepler_to_rv(orb0)
# v0 = v0 * 10/7
# T0 = orbital_period(orb0, GM_EARTH)
plot_orbit(rv_to_kepler(r0, v0, orb0.t))
## CURTIS chap 3

model = Model(
    optimizer_with_attributes(Ipopt.Optimizer,
        "max_iter" => 10000,
        "max_wall_time" => 30.0)
)

i = @variable(model, lower_bound = 0, upper_bound = 180, base_name = "i")
Ω = @variable(model, base_name = "Ω")
ω = @variable(model, base_name = "ω")
nu = @variable(model, base_name = "nu")
# @variables(model, begin
#     # 0 <= i <= 180
#     # Ω
#     ω
#     nu
# end)
# @variable(model, r[1:3])
# @variable(model, v[1:3])
tol = 1e-3
rnorm = norm(r0)
vnorm = norm(v0)
vr = dot(r0/rnorm, v0)

h = cross(r0, v0)

normal_direction = h / √(h' * h)

@constraint(model, -tol <= cosd(i) - normal_direction[3] <= tol)

N = cross([0;0;1], h)

Nnorm = √(N' * N)

@constraint(model, -tol <= cosd(Ω) - N[1]/Nnorm <= tol)
@constraint(model, -tol <= sind(Ω) - N[2]/Nnorm <= tol)

exc_vec = (vnorm^2 / GM_EARTH - 1 / rnorm) * r0 - rnorm*vr/GM_EARTH .* v0

e = √(exc_vec' * exc_vec)

@constraint(model, -tol <= cosd(ω) - dot(N, exc_vec) / (Nnorm*e) <= tol)

N_e_cross = cross(N, exc_vec)
normal_N_e_cross = dot(N_e_cross, normal_direction)

@constraint(model, -tol <= sind(ω) - normal_N_e_cross / (Nnorm*e) <= tol)

@constraint(model, -tol <= cosd(nu) - dot(exc_vec, r0) / (e*rnorm) <= tol)

exc_r_cross = cross(exc_vec, r0)
normal_exc_r_cross = dot(exc_r_cross, normal_direction)
@constraint(model, -tol <= sind(nu) -  normal_exc_r_cross / (e*rnorm) <= tol)

model

##
optimize!(model)
value.(all_variables(model))
##
value(e)
##
all_variables(model)