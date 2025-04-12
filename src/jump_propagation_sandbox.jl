using SatelliteToolboxBase
using SatelliteToolboxPropagators
using JuMP
using Ipopt
using Setfield
using ForwardDiff
include("TG.jl")
using .TG
using LinearAlgebra
## example 3.2
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

Vorb_sup = √(GM_EARTH/EARTH_EQUATORIAL_RADIUS)

r0 = @variable(model, [1:3], start = EARTH_EQUATORIAL_RADIUS)
v0 = @variable(model, [1:3])
set_start_value(v0[1], Vorb_sup)

r1 = @variable(model, [1:3], start = EARTH_EQUATORIAL_RADIUS)
v1 = @variable(model, [1:3])
set_start_value(v1[1], Vorb_sup)

#adding exc as variable so bounds will always be respected
#then need to put constraint on it and implement E and M
#deg!!!
e = @variable(model, lower_bound = 0, upper_bound = 1) 
i = @variable(model, lower_bound = 0, upper_bound = 180, base_name = "i")
Ω = @variable(model, base_name = "Ω")
ω = @variable(model, base_name = "ω")
nu0 = @variable(model, base_name = "nu0")
nu1 = @variable(model, base_name = "nu1")

#rad!!!
M0 = @variable(model, base_name = "M0")
E0 = @variable(model, base_name= "E0")
M1 = @variable(model, base_name = "M1")
E1 = @variable(model, base_name= "E1")

tol = 1e-9

r0norm = √(r0' * r0)
v0norm = √(v0' * v0)

r1norm = √(r1' * r1)
v1norm = √(v1' * v1)

a = -1 / (- 2 / r0norm + v0norm^2 / GM_EARTH)

vr0 = dot(r0 ./ r0norm, v0)

h0 = cross(r0, v0)

normal_direction = h0 ./ √(h0' * h0)

@constraint(model, cosd(i) == normal_direction[3])

N = cross([0;0;1], h0)

Nnorm = √(N' * N)

@constraint(model, cosd(Ω) == N[1]/Nnorm)
@constraint(model, -tol <= sind(Ω) - N[2]/Nnorm <= tol)

exc_vec = (v0norm^2 / GM_EARTH - 1 / r0norm) * r0 - r0norm*vr0/GM_EARTH .* v0

exc_vec_norm = √(exc_vec' * exc_vec)

@constraint(model, e == exc_vec_norm)

@constraint(model, cosd(ω) == dot(N, exc_vec) / (Nnorm*exc_vec_norm))

N_e_cross = cross(N, exc_vec)
normal_N_e_cross = dot(N_e_cross, normal_direction)

@constraint(model, -tol <= sind(ω) - normal_N_e_cross / (Nnorm*exc_vec_norm) <= tol)

#nu

@constraint(model, cosd(nu0) == dot(exc_vec, r0) / (exc_vec_norm*r0norm))

exc_r_cross0 = cross(exc_vec, r0)
normal_exc_r_cross0 = dot(exc_r_cross0, normal_direction)
@constraint(model, -tol <= sind(nu0) -  normal_exc_r_cross0 / (exc_vec_norm*r0norm) <= tol)


@constraint(model, cosd(nu1) == dot(exc_vec, r1) / (exc_vec_norm*r1norm))

exc_r_cross1 = cross(exc_vec, r1)
normal_exc_r_cross1 = dot(exc_r_cross1, normal_direction)
@constraint(model, -tol <= sind(nu1) -  normal_exc_r_cross1 / (exc_vec_norm*r1norm) <= tol)

# M

@constraint(model, E0 - e*sin(E0) == M0)
@constraint(model, E1 - e*sin(E1) == M1)

# E
#curtis page 144 & 145
@constraint(model, cos(E0) == (e + cosd(nu0)) / (1 + e*cosd(nu0)))
@constraint(model, -tol <= sin(E0) - √(1-e^2)*sind(nu0) / (1 + e*cosd(nu0)) <= tol)

@constraint(model, cos(E1) == (e + cosd(nu1)) / (1 + e*cosd(nu1)))
@constraint(model, -tol <= sin(E1) - √(1-e^2)*sind(nu1) / (1 + e*cosd(nu1)) <= tol)

# @constraint(model, r0 .== ri)
# @constraint(model, v0 .== vi)


##
optimize!(model)
##
value.(r1)
##
value.(v1)
##
plot_orbit(
    orbi,
    rv_to_kepler(value.(r0), value.(v0)),
    rv_to_kepler(value.(r1), value.(v1))
)