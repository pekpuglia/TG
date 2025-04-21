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
##
v1, v2 = lambert(r1, r2, delta_t)
##
orb1 = rv_to_kepler(r1, v1)
orb2 = rv_to_kepler(r2, v2)
##
plot_orbit(orb1, orb2)
##
orb1 = rv_to_kepler(r1, v1)
orb2 = rv_to_kepler(r2, v2)
TG.plot_orbit(orb1, orb2)
##
#GLANDORF expressions
R0 = r1
V0 = v1

R = r2
V = v2
r = norm(R)

t = delta_t

H = cross(R0, V0)
h = norm(H)

p = orb1.a * (1 - orb1.e^2)
θ = orb2.f
e = orb1.e
##
g_glandorf = ((p/h)*(p+r)*r*sin(θ) - 3*e*p*delta_t)/(1-e^2)
f_glandorf = - (p/h) * (p+r)*r*cos(θ)
## R, V, r, theta funções do tempo
F1 = r*cos(θ)
F2 = r*sin(θ)
F3 = - (h/p)*sin(θ)
F4 = (h/p)*(e + cos(θ))
##
F5 = GM_EARTH / r^3
F6 = 3t
F7 = 3GM_EARTH*t/r^3
F8 = F3 + g_glandorf*F5
F9 = F4 + f_glandorf*F5
##
P = [
    F1*R-g_glandorf*V F2*R-f_glandorf*V 2*R-F6*V  V    F1*H F2*H
    F8*R-F1*V         F9*R-F2*V         F7*R-V   -F5*R F3*H F4*H
]
##
σ = cross(R, H)
w = cross(V, H)
b1 =    h + F5*(F1*f_glandorf - F2*g_glandorf)
b2 = F6*h + F3*f_glandorf - F4*g_glandorf
b3 = F6*h + 2*(F3*f_glandorf - F4*g_glandorf)
b4 = F1*f_glandorf - F2*g_glandorf
Pinv = 1/h^3 * [
     F2*F5*σ' - F4*w'  2*F4*σ' - F2*w'
    -F1*F5*σ' + F3*w' -2*F3*σ' + F1*w'
    h*w'              -h*σ'
    -b1*σ'+b2*w'      -b3*σ' + b4*w'
    F4*H'             -F2*H'
    -F3*H'             F1*H'
]
##
P*Pinv