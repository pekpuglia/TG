
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

function solve_for_z(delta_t, A, r1n, r2n)
    model = Model(Ipopt.Optimizer)

    z = @variable(model, z >= 0)
    
    y = r1n + r2n + A*(z*S(z) - 1) / √(C(z))

    @constraint(model, √GM_EARTH*delta_t == (y/C(z))^(3/2)*S(z) + A*sqrt(y))
    
    optimize!(model)

    value(z), value(y)
end
export lambert
#curtis chap 5
function lambert(r1, r2, deltat)
    r1n = norm(r1)
    r2n = norm(r2)

    c = cross(r1, r2)
    d = dot(r1, r2)
    #prograde trajectories
    delta_theta = begin
        dtheta = acos(clamp(d / (r1n*r2n), -1, 1))
        
        if c[3] >= 0
            dtheta
        else
            2π - dtheta
        end
    end

    A = sin(delta_theta) * sqrt(r1n*r2n / (1 - cos(delta_theta)))

    z, y = solve_for_z(deltat, A, r1n, r2n)

    #curtis 5.46
    f = 1 - y/r1n
    g = A*√(y/GM_EARTH)
    gdot = 1 - y / r2n

    v1 = 1/g * (r2 - f*r1)
    v2 = 1/g * (gdot * r2 - r1)

    v1, v2
end