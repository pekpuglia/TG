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
hohmann_start_phase_frac = 1
a1 = 7000e3
a2 = 8000e3

hohmann_time = orbital_period((a1+a2)/2, GM_EARTH) / 2
initial_coast_time = hohmann_start_phase_frac * extra_phase / 360 * orbital_period(a1, GM_EARTH)
terminal_coast_time = (extra_phase*(1 - hohmann_start_phase_frac)) / 360 * orbital_period(a2, GM_EARTH)

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
    deg2rad(180+extra_phase) + orb1.f
)
##
r1, v1             = kepler_to_rv(orb1)
r2, v2             = kepler_to_rv(orb2)
v1sol, v2sol = lambert(r1, r2, (orb2.t - orb1.t)*86400, false, RAAN = orb1.Ω, i=orb1.i)
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

c0(x) = cos(√x)
c1(x) = sin(√x)/√x
c2(x) = (1-cos(√x))/x
c3(x) = (√x - sin(√x))/(x*√x)

u(x, rho) = √(1 - rho*c1(x)/√c2(x))

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
function parab_time_lambert(r1, r2, t; prograde=true, RAAN = nothing, i = nothing, setsilent=true)
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
        LambertResult(true, sigma, sigma_par, r1, v1, r2, v2, t, propagator)
    catch ArgumentError #no idea where this comes from
        LambertResult(false, sigma, sigma_par, r1, v1, r2, v2, t, nothing)
    end
end

struct LambertTransfer2
    orb1::KeplerianElements
    orb2::KeplerianElements
    lambert::LambertResult
    p0_p0dot
end

function LambertTransfer2(orb1::KeplerianElements, orb2::KeplerianElements)
    r1, v1             = kepler_to_rv(orb1)
    r2, v2             = kepler_to_rv(orb2)
    lamb = parab_time_lambert(r1, r2, (orb2.t - orb1.t)*86400, RAAN=orb1.Ω, i=orb1.i)

    if !lamb.is_elliptic
        return LambertTransfer2(orb1, orb2, lamb, nothing)
    end

    propagator = Propagators.init(Val(:TwoBody), lamb.propagator.tbd.orb₀)

    deltaV0 = lamb.v1 - v1
    deltaVf = v2 - lamb.v2

    p0 = deltaV0 / norm(deltaV0)
    pf = deltaVf / norm(deltaVf)

    Phi = P_glandorf(r2, lamb.v2, transfer_time) * Pinv_glandorf(r1, lamb.v1, 0)
    M = Phi[1:3, 1:3]
    N = Phi[1:3, 4:6]

    if abs(det(N)) <= 1e-10
        #CHECK THIS IS TRUE
        #only works for coplanar transfers?
        A = N * [r1 lamb.v1]
        b = pf - M * p0
        p0dot = [r1 lamb.v1] * (A \ b)
    else
        p0dot = N \ (pf - M * p0)
    end

    LambertTransfer2(orb1, orb2, lamb, [p0; p0dot])
end

transfer_start(lt::LambertTransfer2) = lt.lambert.propagator.tbd.orb₀
transfer_end(lt::LambertTransfer2) = lt.lambert.propagator.tbd.orbk

function cost(lt::LambertTransfer2)
    if lt.lambert.is_elliptic
        norm(lt.lambert.v1 - kepler_to_rv(orb1)[2]) + norm(lt.lambert.v2 - kepler_to_rv(orb2)[2])
    else
        NaN
    end
end

function cost(::Nothing)
    NaN
end

p(lt::LambertTransfer2, t) = Phi_time(lt.lambert.propagator, t)[1:3, :] * lt.p0_p0dot
pdot(lt::LambertTransfer2, t) = Phi_time(lt.lambert.propagator, t)[4:6, :] * lt.p0_p0dot
##
transfer = LambertTransfer2(orb1, orb2)
##
f = plot_orbit(orb1, orb2,transfer_start(transfer))
##
# save("./src/late_hohmann.png", f)
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
# save("./src/late_hohmann_p_history.png", f)
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
_, _, prop = Propagators.propagate(Val(:TwoBody), -transfer_time, orb2)
orb2_init = prop.tbd.orbk
##
N = 100
time_offsets1 =  range(0, 2T2, N)
time_offsets2 = range(0, 2T2, N)
orb1_porkchop = [Propagators.propagate(Val(:TwoBody), t, orb1)[3].tbd.orbk for t in time_offsets1]
orb2_porkchop = [Propagators.propagate(Val(:TwoBody), t, orb2_init)[3].tbd.orbk for t in time_offsets2]
##
function porkchop_transfer(orb1, orb2)
    if orb1.t >= orb2.t
        return nothing
    else
        return LambertTransfer2(orb1, orb2)
    end
end
##
transfer_grid = porkchop_transfer.(orb1_porkchop, permutedims(orb2_porkchop))
##
f = Figure()
ax = Axis(f[1, 1], xlabel= "departure", ylabel = "arrival", title="Cost of transfer", aspect=1.0)
cont = contourf!(ax, 1:N, 1:N, cost.(transfer_grid))
hyperbolic = findall(x -> !isnothing(x) && !x.lambert.is_elliptic, transfer_grid)
scatter!(ax, getfield.(hyperbolic, :I) .|> first, getfield.(hyperbolic, :I) .|> last, color=:red)
Colorbar(f[1, 2], cont)
f
## TODO
# add retrograde maneuver?/exclude weird retrograde cases
# optimize initial coast
# better time frames
# check sigma > sigma_par
#adicionar ifs nas funções de stumpff
##
# save("./src/better_porkchop.png", f)