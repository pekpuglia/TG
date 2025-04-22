using SatelliteToolboxBase
using SatelliteToolboxPropagators
include("TG.jl")
using .TG
using Setfield
using LinearAlgebra
##
extra_phase = 10
orb1 = KeplerianElements(
    date_to_jd(2023, 1, 1, 0, 0, 0),
    7000e3,
    0.0,
    51 |> deg2rad,
    0    |> deg2rad,
    0     |> deg2rad,
    -extra_phase     |> deg2rad
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
## hohmann might cause N to be non invertible
transfer_a = (orb1.a + orb2.a) / 2

transfer_orbit_time = 2π√(transfer_a^3 / GM_EARTH) / 2

min_transfer_time = transfer_orbit_time + extra_phase/360 * orbital_period(orb1, GM_EARTH)
max_transfer_time = transfer_orbit_time + extra_phase/360 * orbital_period(orb2, GM_EARTH)
transfer_time = max_transfer_time
##
r1, v1             = kepler_to_rv(orb1)
r2, v2             = kepler_to_rv(orb2)
vt_start, vt_end = lambert(r1, r2, transfer_time)
##
plot_orbit(orb1, orb2, rv_to_kepler(r1, vt_start))
##
# r1, v1             = kepler_to_rv(orb1)
# rt_start, vt_start = kepler_to_rv(transfer_start)
# rt_end, vt_end     = kepler_to_rv(transfer_end)
# r2, v2             = kepler_to_rv(orb2)
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
p0dot = N \ (pf - M * p0)
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