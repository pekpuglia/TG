c0(x) = cos(√x)
c1(x) = sin(√x)/√x
c2(x) = (1-cos(√x))/x
c3(x) = (√x - sin(√x))/(x*√x)

u(x, rho) = √(1 - rho*c1(x)/√c2(x))

export LambertResult, lambert
#sukhanov 7
struct LambertResult
    is_elliptic
    sigma
    sigma_par
    r1
    v1
    r2
    v2
    t
    propagator
end

#sukhanov 7
function lambert(r1, r2, t; prograde=true, RAAN = nothing, i = nothing, setsilent=true)
    r1n = norm(r1)
    r2n = norm(r2)


    c = cross(r1, r2)
    d = dot(r1, r2)
    #prograde trajectories
    phi = begin
        dphi = acos(clamp(d / (r1n*r2n), -1, 1))
        
        if prograde
            if c[3] >= 0
                dphi
            else
                2π - dphi
            end
        else
            if c[3] < 0
                dphi
            else
                2π - dphi
            end
        end
    end

    sigma = √GM_EARTH / (r1n + r2n)^(3/2) * t

    rho = √(2*r1n*r2n) / (r1n + r2n) * cos(phi/2)

    sigma_par = 1/3*(√2 + rho)*√(1-√2*rho)

    if sigma <= sigma_par
        return LambertResult(false, sigma, sigma_par, r1, nothing, r2, nothing, t, nothing)
    end

    model = Model(Ipopt.Optimizer)
    if setsilent
        set_silent(model)
    end

    #7.52 - x initial condition
    um = √(1-√2*abs(rho))

    eps0 = (π/(2/3*um^3+sigma-rho*um))^(1//3)*um

    z0 = π-eps0

    x0 = 2*z0^2

    #sukhanov 7.32 says x>0 elliptic orbit
    x = @variable(model, lower_bound = 0, start = x0)

    @constraint(model, c3(x)/c2(x)^(3/2)*u(x, rho)^3 + rho*u(x, rho) == sigma)

    optimize!(model)

    xsol = value(x)
    ssol = value(√((r1n+r2n)/(GM_EARTH*c2(x)))*u(x, rho))
    
    if norm(c) / (r1n*r2n) > 1e-6
        f = 1 - GM_EARTH*ssol^2*c2(xsol)/r1n
        g = t - GM_EARTH*ssol^3*c3(xsol)
        gdot = 1 - GM_EARTH*ssol^2*c2(xsol)/r2n

        v1 = 1/g * (r2 - f*r1)
        v2 = 1/g * (gdot * r2 - r1)
    else
        #colinear case
        h = - GM_EARTH*c2(xsol)*xsol / ((r1n+r2n)*u(xsol, value(rho))^2)

        #c1 should be very close to 0 around x = pi^2
        #to ensure correctness vr is set to 0 when x-pi^2<1e-6
        if abs(xsol-pi^2) < 1e-6
            vr1 = 0
        else #doesn't work
            vr1 = -r1n*c1(xsol) / (ssol*c2(xsol))
        end

        try
            vn1 = √(h - vr1^2 + 2GM_EARTH/r1n) #should change prograde/retrograde
        catch e
            if e isa DomainError
                return LambertResult(false, sigma, sigma_par, r1, nothing, r2, nothing, t, nothing)
            end
        end

        orbit_normal = [
            sin(RAAN)*sin(i)
            -cos(RAAN)*sin(i)
            cos(i)
        ]

        r1dir = r1 / r1n
        ndir1 = cross(orbit_normal, r1dir)

        v1 = r1dir*vr1 + vn1*ndir1
        
        vn2 = r1n*vn1/r2n

        vr2_squared = (h - vn2^2 + 2GM_EARTH/r2n)

        if abs(vr2_squared) < 1e-6
            vr2_squared = 0.0
        end

        vr2 = √vr2_squared

        r2dir = r2/r2n
        ndirf = cross(orbit_normal, r2dir)

        v2 = r2dir*vr2 + vn2*ndirf
    end

    try
        propagator = Propagators.init(Val(:TwoBody), rv_to_kepler(r1, v1))
        rv_to_kepler(r2, v2) # sometimes fails - don't know why
        LambertResult(true, sigma, sigma_par, r1, v1, r2, v2, t, propagator)
    catch ArgumentError #no idea where this comes from
        LambertResult(false, sigma, sigma_par, r1, v1, r2, v2, t, nothing)
    end
end