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
function add_orbital_elements!(model)
    Vorb_sup = √(GM_EARTH/EARTH_EQUATORIAL_RADIUS)
    r = @variable(model, [1:3], start = EARTH_EQUATORIAL_RADIUS)
    v = @variable(model, [1:3])
    set_start_value(v[1], Vorb_sup)

    #adding exc as variable so bounds will always be respected
    #then need to put constraint on it and implement E and M
    #deg!!!
    e = @variable(model, lower_bound = 0, upper_bound = 1) 
    i = @variable(model, lower_bound = 0, upper_bound = 180, base_name = "i")
    Ω = @variable(model, base_name = "Ω")
    ω = @variable(model, base_name = "ω")
    nu = @variable(model, base_name = "nu")

    #rad!!!
    M = @variable(model, base_name = "M")
    E = @variable(model, base_name= "E")
    
    tol = 1e-9
    rnorm = √(r' * r)
    vnorm = √(v' * v)
    
    a = -1 / (- 2 / rnorm + vnorm^2 / GM_EARTH)
    
    vr = dot(r ./ rnorm, v)
    
    h = cross(r, v)
    
    normal_direction = h ./ √(h' * h)
    
    @constraint(model, cosd(i) == normal_direction[3])
    
    N = cross([0;0;1], h)
    
    Nnorm = √(N' * N)
    
    @constraint(model, cosd(Ω) == N[1]/Nnorm)
    @constraint(model, -tol <= sind(Ω) - N[2]/Nnorm <= tol)
    
    exc_vec = (vnorm^2 / GM_EARTH - 1 / rnorm) * r - rnorm*vr/GM_EARTH .* v
    
    exc_vec_norm = √(exc_vec' * exc_vec)

    @constraint(model, e == exc_vec_norm)
    
    @constraint(model, cosd(ω) == dot(N, exc_vec) / (Nnorm*exc_vec_norm))
    
    N_e_cross = cross(N, exc_vec)
    normal_N_e_cross = dot(N_e_cross, normal_direction)
    
    @constraint(model, -tol <= sind(ω) - normal_N_e_cross / (Nnorm*exc_vec_norm) <= tol)
    
    @constraint(model, cosd(nu) == dot(exc_vec, r) / (exc_vec_norm*rnorm))
    
    exc_r_cross = cross(exc_vec, r)
    normal_exc_r_cross = dot(exc_r_cross, normal_direction)
    @constraint(model, -tol <= sind(nu) -  normal_exc_r_cross / (exc_vec_norm*rnorm) <= tol)
    
    @constraint(model, E - e*sin(E) == M)

    #curtis page 144 & 145
    @constraint(model, cos(E) == (e + cosd(nu)) / (1 + e*cosd(nu)))
    @constraint(model, -tol <= sin(E) - √(1-e^2)*sind(nu) / (1 + e*cosd(nu)) <= tol)

    r, v, a, e, i, Ω, ω, nu, M, E
end
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