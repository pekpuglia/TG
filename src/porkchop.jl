using SatelliteToolboxBase
using SatelliteToolboxPropagators
include("TG.jl")
using .TG
using LinearAlgebra
using GLMakie
##
initial_phase_offset = -80
a1 = 7000e3
a2 = 9000e3

hohmann_time = orbital_period((a1+a2)/2, GM_EARTH) / 2
##
orb1 = KeplerianElements(
    date_to_jd(2023, 1, 1, 0, 0, 0),
    a1,
    0.0,
    51 |> deg2rad,
    0    |> deg2rad,
    0     |> deg2rad,
    60     |> deg2rad
)

orb2 = KeplerianElements(
    orb1.t,
    a2,
    orb1.e,
    orb1.i,
    orb1.Ω,
    orb1.ω,
    deg2rad(-80) + orb1.f
)
##
plot_orbit(orb1, orb2)
## fixed departure
N = 100
flight_times = range(0, 2orbital_period(orb2, GM_EARTH), N)
##
i = 50
orb2_arrival = Propagators.propagate(Val(:TwoBody), flight_times[i], orb2)[3].tbd.orbk
##
plot_orbit(orb1, orb2, orb2_arrival)
##
transfer = CoastTransfer(orb1, orb2_arrival)
##
plot_transfer(transfer)
##
orb2_arrival_list = [Propagators.propagate(Val(:TwoBody), time, orb2)[3].tbd.orbk for time in flight_times]
##
transfer_list = CoastTransfer.(Ref(orb1), orb2_arrival_list)
##
plot(flight_times, cost.(transfer_list))