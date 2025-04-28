using SatelliteToolboxBase
using SatelliteToolboxPropagators
include("TG.jl")
using .TG
using Setfield
using LinearAlgebra
using GLMakie
using JuMP
using Ipopt
## estimate of transfer_time
extra_phase = 0
hohmann_start_phase_frac = 0
a1 = 7000e3
a2 = 8000e3

hohmann_time = orbital_period((a1+a2)/2, GM_EARTH) / 2
initial_coast_time = max(hohmann_start_phase_frac * extra_phase / 360 * orbital_period(a1, GM_EARTH), 0)
terminal_coast_time = max((extra_phase*(1 - hohmann_start_phase_frac)) / 360 * orbital_period(a2, GM_EARTH), 0)

transfer_time = initial_coast_time + hohmann_time + terminal_coast_time

##
orb1 = KeplerianElements(
    date_to_jd(2023, 1, 1, 0, 0, 0),
    a1,
    0.0,
    51 |> deg2rad,
    0    |> deg2rad,
    0     |> deg2rad,
    0     |> deg2rad
)

orb2 = KeplerianElements(
    orb1.t + transfer_time / 86400,
    a2,
    0.0,
    orb1.i,
    orb1.Ω,
    orb1.ω,
    180+extra_phase     |> deg2rad
)
##
using TG: c1, c2, c3, u
function sukhanov_lambert(r1, r2, t; RAAN = nothing, i = nothing)
    r1n = norm(r1)
    r2n = norm(r2)


    c = cross(r1, r2)
    d = dot(r1, r2)
    #prograde trajectories
    phi = begin
        dphi = acos(clamp(d / (r1n*r2n), -1, 1))
        
        if c[3] >= 0
            dphi
        else
            2π - dphi
        end
    end

    sigma = √GM_EARTH / (r1n + r2n)^(3/2) * t

    rho = √(2*r1n*r2n) / (r1n + r2n) * cos(phi/2)

    model = Model(Ipopt.Optimizer)
    
    #sukhanov 7.32 says x>0 elliptic orbit
    x = @variable(model, lower_bound = 0)

    @constraint(model, c3(x)/c2(x)^(3/2)*u(x, rho)^3 + rho*u(x, rho) == sigma)

    optimize!(model)

    xsol = value(x)
    println(xsol)
    ssol = value(√((r1n+r2n)/(GM_EARTH*c2(x)))*u(x, rho))
    println(ssol)
    
    if norm(c) / (r1n*r2n) > 1e-6
        println("non colinear")

        f = 1 - GM_EARTH*ssol^2*c2(xsol)/r1n
        g = t - GM_EARTH*ssol^3*c3(xsol)
        gdot = 1 - GM_EARTH*ssol^2*c2(xsol)/r2n

        v1 = 1/g * (r2 - f*r1)
        v2 = 1/g * (gdot * r2 - r1)
    else
        println("colinear")
        #colinear case
        h = - GM_EARTH*c2(xsol)*xsol / ((r1n+r2n)*u(xsol, value(rho))^2)

        #c1 should be very close to 0 around x = pi^2
        #to ensure correctness vr is set to 0 when x-pi^2<1e-6
        if abs(xsol-pi^2) < 1e-6
            vr1 = 0
        else
            vr1 = -r1n*c1(xsol) / (ssol*c2(xsol))
        end

        vn1 = √(h - vr1^2 + 2GM_EARTH/r1n)

        orbit_normal = [
            sin(RAAN)*sin(i)
            -cos(RAAN)*sin(i)
            cos(i)
        ]

        r1dir = r1 / r1n
        ndir1 = cross(orbit_normal, r1dir)

        v1 = r1dir*vr1 + vn1*ndir1
        
        vn2 = r1n*vn1/r2n

        vr2 = √(h - vn2^2 + 2GM_EARTH/r2n)

        r2dir = r2/r2n
        ndirf = cross(orbit_normal, r2dir)

        v2 = r2dir*vr2 + vn2*ndirf
    end
    v1, v2
end
##
r1, v1             = kepler_to_rv(orb1)
r2, v2             = kepler_to_rv(orb2)
v1sol, v2sol = sukhanov_lambert(r1, r2, (orb2.t - orb1.t)*86400, RAAN = orb1.Ω, i=orb1.i)
dot(v1sol, r1/norm(r1))
##
plot_orbit(
    rv_to_kepler(r1, v1sol),
    rv_to_kepler(r1, v1),
    rv_to_kepler(r2, v2sol),
    rv_to_kepler(r2, v2),
)
##
hohmann_start = KeplerianElements(
    date_to_jd(2023, 1, 1, 0, 0, 0),
    (orb1.a + orb2.a) / 2,
    (orb2.a - orb1.a) / (orb1.a + orb2.a),
    orb1.i,
    orb1.Ω,
    orb1.ω + deg2rad(hohmann_start_phase_frac*extra_phase),
    0 |> deg2rad
)
hohmann_end = @set hohmann_start.f = π
##
function Phi_time(propagator, t)
    r, v = Propagators.propagate!(propagator, t)
    r0, v0 = kepler_to_rv(propagator.tbd.orb₀)
    Phi = P_glandorf(r, v, t) * Pinv_glandorf(r0, v0, 0)
    Phi
end

struct LambertTransfer1
    orb1::KeplerianElements
    orb2::KeplerianElements
    v_transfer_start
    v_transfer_end
    transfer_propagator
    p0_p0dot
end

function LambertTransfer1(orb1::KeplerianElements, orb2::KeplerianElements)
    r1, v1             = kepler_to_rv(orb1)
    r2, v2             = kepler_to_rv(orb2)
    vt_start, vt_end = lambert(r1, r2, (orb2.t - orb1.t)*86400, RAAN=orb1.Ω, i=orb1.i)
    
    orbt1, orbt2 = try
        rv_to_kepler(r1, vt_start), rv_to_kepler(r2, vt_end)
    catch ArgumentError
        return LambertTransfer1(orb1, orb2, nothing, nothing, nothing, nothing)
    end

    propagator = Propagators.init(Val(:TwoBody), orbt1)

    deltaV0 = vt_start - v1
    deltaVf = v2 - vt_end

    p0 = deltaV0 / norm(deltaV0)
    pf = deltaVf / norm(deltaVf)

    Phi = P_glandorf(r2, vt_end, transfer_time) * Pinv_glandorf(r1, vt_start, 0)
    M = Phi[1:3, 1:3]
    N = Phi[1:3, 4:6]

    if abs(det(N)) <= 1e-10
        #CHECK THIS IS TRUE
        #only works for coplanar transfers?
        A = N * [r1 vt_start]
        b = pf - M * p0
        p0dot = [r1 vt_start] * (A \ b)
    else
        p0dot = N \ (pf - M * p0)
    end

    LambertTransfer1(orb1, orb2, vt_start, vt_end, propagator, [p0; p0dot])
end

transfer_start(lt::LambertTransfer1) = rv_to_kepler(kepler_to_rv(orb1)[1], lt.v_transfer_start)
transfer_end(lt::LambertTransfer1) = rv_to_kepler(kepler_to_rv(orb2)[1], lt.v_transfer_end)

function cost(lt::LambertTransfer1)
    if isnothing(lt.v_transfer_start)
        NaN
    else
        norm(lt.v_transfer_start - kepler_to_rv(orb1)[2]) + norm(lt.v_transfer_end - kepler_to_rv(orb2)[2])
    end
end

p(lt::LambertTransfer1, t) = Phi_time(lt.transfer_propagator, t)[1:3, :] * lt.p0_p0dot
pdot(lt::LambertTransfer1, t) = Phi_time(lt.transfer_propagator, t)[4:6, :] * lt.p0_p0dot
##
transfer = LambertTransfer1(orb1, orb2)
##
f = plot_orbit(orb1, orb2,transfer_start(transfer))
##
save("./src/late_hohmann.png", f)
##
t_list = transfer_time * (0.0:0.01:1.0)
p_history = p.(Ref(transfer), t_list)
##
pnorm_history = norm.(p_history)
##
#analisar
#descobrir por que N não é inversível no caso Hohmann e o que fazer
f = lines(t_list, pnorm_history)
##
save("./src/late_hohmann_p_history.png", f)
##
p_x = getindex.(p_history, 1)
p_y = getindex.(p_history, 2)
p_z = getindex.(p_history, 3)
##
lines(p_x, p_y, p_z)
## porkchop plot over 1 synodic period
T1 = orbital_period(a1, GM_EARTH)
T2 = orbital_period(a2, GM_EARTH)
Tsyn = T1*T2 / abs(T2 - T1)
##
N = 100
time_offsets1 =  range(-T1/5, T1/5, N)
time_offsets2 = range(-T2/5, T2/5, N)
orb1_porkchop = [Propagators.propagate(Val(:TwoBody), t, orb1)[3].tbd.orbk for t in time_offsets1]
orb2_porkchop = [Propagators.propagate(Val(:TwoBody), t, orb2)[3].tbd.orbk for t in time_offsets2]
##
porkchop(orb1, orb2) = (orb2.t > orb1.t) ? LambertTransfer1(orb1, orb2) : nothing
transfer_porkchop = porkchop.(orb1_porkchop, permutedims(orb2_porkchop))
##
cost(x::Nothing) = NaN
##
f = Figure()
ax = Axis(f[1, 1], xlabel= "departure", ylabel = "arrival", title="Cost of transfer")
cont = contourf!(ax, time_offsets1, time_offsets2, cost.(transfer_porkchop))
Colorbar(f[1, 2], cont)
f
##
save("./src/porckchop_hohmann2.png", f)