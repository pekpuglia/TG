using LinearAlgebra
using JuMP
using Ipopt
using SatelliteToolboxBase
include("TG.jl")
using .TG
## curtis 5.2
r1 = [
    5000
    10000
    2100
]*1e3

r2 = [
    -14600
    2500
    7000
]*1e3

delta_t = 3600
##
function S(z)
    # if z > 0
    (√(z) - sin(√z)) / (z ^ (3/2))
    # elseif z < 0
    #     (sinh(√-z) - √(-z)) / (z ^ (3/2))
    # else
    #     1/6
    # end
end

function C(z)
    # if z > 0
    (1 - cos(√z)) / z
    # elseif z < 0
    #     (cosh(√-z) - 1) / (-z)
    # else
    #     1/2
    # end
end

# function y(z, r1n, r2n, A)
#     r1n + r2n + A*(z*S(z) - 1) / √(C(z))
# end

function solve_for_z(delta_t, A, r1n, r2n)
    model = Model(Ipopt.Optimizer)

    z = @variable(model, z >= 0)
    
    y = r1n + r2n + A*(z*S(z) - 1) / √(C(z))

    @constraint(model, √GM_EARTH*delta_t == (y/C(z))^(3/2)*S(z) + A*sqrt(y))
    
    optimize!(model)

    value(z), value(y)
end

function lambert(r1, r2, deltat)
    r1n = norm(r1)
    r2n = norm(r2)

    c = cross(r1, r2)
    d = dot(r1, r2)

    delta_theta = begin
        dtheta = acos(d / (r1n*r2n))
        
        if c[3] >= 0
            dtheta
        else
            2π - dtheta
        end
    end

    A = sin(delta_theta) * sqrt(r1n*r2n / (1 - cos(delta_theta)))

    z, y = solve_for_z(delta_t, A, r1n, r2n)

    f = 1 - y/r1n
    g = A*√(y/GM_EARTH)
    gdot = 1 - y / r2n

    v1 = 1/g * (r2 - f*r1)
    v2 = 1/g * (gdot * r2 - r1)

    v1, v2
end
## add test for this
v1, v2 = lambert(r1, r2, delta_t)
##
orb1 = rv_to_kepler(r1, v1)
orb2 = rv_to_kepler(r2, v2)
##
plot_orbit(orb1, orb2)
##
#GLANDORF expressions

g_glandorf = ((p/h)*(p+r)*r*sin(θ) - 3*e*p*delta_t)/(1-e^2)
## R, V, r, theta funções do tempo
#computes more than needed for P but ok
function glandorf_precompute(R, V, t)
    r = norm(R)

    H = cross(R, V)
    h = norm(H)
    
    orb = rv_to_kepler(R, V)

    p = orb.a * (1 - orb.e^2)
    θ = orb.f
    e = orb.e

    g_glandorf = ((p/h)*(p+r)*r*sin(θ) - 3*e*p*t)/(1-e^2)
    f_glandorf = - (p/h) * (p+r)*r*cos(θ)

    F3 = - (h/p)*sin(θ)
    F4 = (h/p)*(e + cos(θ))
    F5 = GM_EARTH / r^3

    F = [
        r*cos(θ)
        r*sin(θ)
        F3
        F4
        F5
        3t
        3GM_EARTH*t/r^3
        F3 + g_glandorf*F5
        F4 + f_glandorf*F5
    ]

    σ = cross(R, H)
    w = cross(V, H)
    b1 =    h + F[5]*(F[1]*f_glandorf - F[2]*g_glandorf)
    b2 = F[6]*h + F[3]*f_glandorf - F[4]*g_glandorf
    b3 = F[6]*h + 2*(F[3]*f_glandorf - F[4]*g_glandorf)
    b4 = F[1]*f_glandorf - F[2]*g_glandorf

    b = [b1; b2;b3;b4]

    F, g_glandorf, f_glandorf, σ, w, b
end
##
R = r2
V = v2
H = cross(R, V)
t = delta_t
F, g_glandorf, f_glandorf, σ, w, b = glandorf_precompute(R, V, t)
##
function PRVt(R, V, t)
    F, g_glandorf, f_glandorf, σ, w, b = glandorf_precompute(R, V, t)
    H = cross(R, V)
    
    [
        F[1]*R-g_glandorf*V F[2]*R-f_glandorf*V 2*R-F[6]*V  V    F[1]*H F[2]*H
        F[8]*R-F[1]*V         F[9]*R-F[2]*V         F[7]*R-V   -F[5]*R F[3]*H F[4]*H
    ]
end
##
P = PRVt(R, V, t)
##
function PinvRVt(R, V, t)
    F, g_glandorf, f_glandorf, σ, w, b = glandorf_precompute(R, V, t)
    H = cross(R, V)
    h = norm(H)
    
    1/h^3 * [
        F[2]*F[5]*σ' - F[4]*w'  2*F[4]*σ' - F[2]*w'
        -F[1]*F[5]*σ' + F[3]*w' -2*F[3]*σ' + F[1]*w'
        h*w'              -h*σ'
        -b[1]*σ'+b[2]*w'      -b[3]*σ' + b[4]*w'
        F[4]*H'             -F[2]*H'
        -F[3]*H'             F[1]*H'
    ]
end
##
Pinv = PinvRVt(R, V, t)
##
P*Pinv
##
PRVt(r1, v1, 0)*PinvRVt(r1, v1, 0)