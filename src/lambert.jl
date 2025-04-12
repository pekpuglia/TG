using LinearAlgebra
using JuMP
using Ipopt
using SatelliteToolboxBase
##
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
#prograde trajectories

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
##
A = sin(delta_theta) * sqrt(r1n*r2n / (1 - cos(delta_theta)))
##
#z = alpha xi^2
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

y(z) = r1n + r2n + A*(z*S(z) - 1) / √(C(z))
##
function solve_for_z(delta_t, A)
    model = Model(Ipopt.Optimizer)

    @variable(model, z >= 0)
    
    @constraint(model, √GM_EARTH*delta_t == (y(z)/C(z))^(3/2)*S(z) + A*sqrt(y(z)))
    
    optimize!(model)

    value(z)
end
##
z = solve_for_z(delta_t, A)
yz = y(z)
##
f = 1 - yz/r1n
g = A*√(yz/GM_EARTH)
gdot = 1 - yz / r2n
##
v1 = 1/g * (r2 - f*r1)
v2 = 1/g * (gdot * r2 - r1)
