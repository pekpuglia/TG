using SatelliteToolboxBase, SatelliteToolboxPropagators
include("../orb_mech.jl")

NONCOP_START = KeplerianElements(
    date_to_jd(2023, 1, 1, 0, 0, 0),
    6748.1e3,
    0.0,
    42.1 |> deg2rad,
    120.2    |> deg2rad,
    0     |> deg2rad,
    175     |> deg2rad
)

NONCOP_END_0 = KeplerianElements(
    NONCOP_START.t,
    6778.1e3,
    0.0,
    42 |> deg2rad,
    120 |> deg2rad,
    0,
    180 |> deg2rad
)

NONCOP_END = Propagators.propagate(Val(:TwoBody), 2*orbital_period(NONCOP_END_0, GM_EARTH), NONCOP_END_0)[3].tbd.orbk

NAME = "Non-coplanar RV"

TRANSFER_TIME =  2*orbital_period(NONCOP_END, GM_EARTH)
