using SatelliteToolboxBase
using SatelliteToolboxPropagators
include("TG.jl")
using .TG
using Setfield
using LinearAlgebra
using GLMakie
## estimate of transfer_time
extra_phase = 0
hohmann_start_phase_frac = 0.5
a1 = 7000e3
a2 = 8000e3

hohmann_time = orbital_period((a1+a2)/2, GM_EARTH) / 2
initial_coast_time = max(hohmann_start_phase_frac * extra_phase / 360 * orbital_period(a1, GM_EARTH), 0)
terminal_coast_time = max((extra_phase*(1 - hohmann_start_phase_frac)) / 360 * orbital_period(a2, GM_EARTH), 0)

transfer_time = initial_coast_time + hohmann_time + terminal_coast_time

##
orb1 = KeplerianElements(
    date_to_jd(2023, 1, 1, 0, 0, 0),
    7000e3,
    0.0,
    51 |> deg2rad,
    0    |> deg2rad,
    0     |> deg2rad,
    0     |> deg2rad
)

orb2 = KeplerianElements(
    orb1.t + transfer_time / 86400,
    8000e3,
    0.0,
    51 |> deg2rad,
    0    |> deg2rad,
    0     |> deg2rad,
    180+extra_phase     |> deg2rad
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
plot_orbit(orb1, orb2, hohmann_start, hohmann_end)
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
    vt_start, vt_end = lambert(r1, r2, (orb2.t - orb1.t)*86400)
    
    propagator = Propagators.init(Val(:TwoBody), rv_to_kepler(r1, vt_start))

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

cost(lt::LambertTransfer1) = norm(lt.v_transfer_start - kepler_to_rv(orb1)[2]) + norm(lt.v_transfer_end - kepler_to_rv(orb2)[2])

p(lt::LambertTransfer1, t) = Phi_time(lt.transfer_propagator, t)[1:3, :] * lt.p0_p0dot
pdot(lt::LambertTransfer1, t) = Phi_time(lt.transfer_propagator, t)[4:6, :] * lt.p0_p0dot
##
transfer = LambertTransfer1(orb1, orb2)
##
plot_orbit(orb1, orb2,transfer_start(transfer))
##
t_list = transfer_time * (0.0:0.01:1.0)
p_history = p.(Ref(transfer), t_list)
##
pnorm_history = norm.(p_history)
##
#analisar
#descobrir por que N não é inversível no caso Hohmann e o que fazer
lines(t_list, pnorm_history)
##
p_x = getindex.(p_history, 1)
p_y = getindex.(p_history, 2)
p_z = getindex.(p_history, 3)
##
lines(p_x, p_y, p_z)