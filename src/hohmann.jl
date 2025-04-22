using SatelliteToolboxBase
using SatelliteToolboxPropagators
include("TG.jl")
using .TG
using Setfield
using LinearAlgebra
##
extra_phase = 0
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
    date_to_jd(2023, 1, 1, 0, 0, 0),
    8000e3,
    0.0,
    51 |> deg2rad,
    0    |> deg2rad,
    0     |> deg2rad,
    180+extra_phase     |> deg2rad
)
##
hohmann_start_phase = extra_phase * 1
hohmann_start = KeplerianElements(
    date_to_jd(2023, 1, 1, 0, 0, 0),
    (orb1.a + orb2.a) / 2,
    (orb2.a - orb1.a) / (orb1.a + orb2.a),
    orb1.i,
    orb1.Ω,
    orb1.ω + deg2rad(hohmann_start_phase),
    0 |> deg2rad
)
hohmann_end = @set hohmann_start.f = π
##
plot_orbit(orb1, orb2, hohmann_start, hohmann_end)
##
hohmann_time = orbital_period(hohmann_start, GM_EARTH) / 2

initial_coast_time = max(hohmann_start_phase / 360 * orbital_period(orb1, GM_EARTH), 0)
terminal_coast_time = max((extra_phase - hohmann_start_phase) / 360 * orbital_period(orb2, GM_EARTH), 0)

transfer_time = initial_coast_time + hohmann_time + terminal_coast_time
##
r1, v1             = kepler_to_rv(orb1)
r2, v2             = kepler_to_rv(orb2)
vt_start, vt_end = lambert(r1, r2, transfer_time)
##
plot_orbit(orb1, orb2, rv_to_kepler(r1, vt_start))
## conway 2.43 → onwards
deltaV0 = vt_start - v1
deltaVf = v2 - vt_end
##
p0 = deltaV0 / norm(deltaV0)
pf = deltaVf / norm(deltaVf)
##
Phi = P_glandorf(r2, vt_end, transfer_time) * Pinv_glandorf(r1, vt_start, 0)
M = Phi[1:3, 1:3]
N = Phi[1:3, 4:6]
S = Phi[4:6, 1:3]
T = Phi[4:6, 4:6]
##
# p0dot = N \ (pf - M * p0)
A = N * [r1 vt_start]
b = pf - M * p0
p0dot = [r1 vt_start] * (A \ b)
## get Phi over time
function Phi_time(propagator, t)
    r, v = Propagators.propagate!(propagator, t)
    r0, v0 = kepler_to_rv(propagator.tbd.orb₀)
    Phi = P_glandorf(r, v, t) * Pinv_glandorf(r0, v0, 0)
    Phi
end
##
propagator = Propagators.init(Val(:TwoBody), rv_to_kepler(r1, vt_start))
t_list = transfer_time * (0.0:0.01:1.0)
p_history = [Phi_time(propagator, t)[1:3, :] * [p0; p0dot] for t in t_list]
##
using GLMakie
##
pnorm_history = norm.(p_history)
##
#analisar
#descobrir por que N não é inversível no caso Hohmann e o que fazer
plot(t_list, pnorm_history)