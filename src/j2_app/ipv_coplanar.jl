using SatelliteToolboxBase
using SatelliteToolboxPropagators
include("../orb_mech.jl")

COPLANAR_START = KeplerianElements(
    date_to_jd(2023, 1, 1, 0, 0, 0),
    6748e3,
    0.0,
    5 |> deg2rad,
    35    |> deg2rad,
    0     |> deg2rad,
    0     |> deg2rad
)

COPLANAR_END_0 = KeplerianElements(
    COPLANAR_START.t,
    6778e3,
    0.0,
    COPLANAR_START.i,
    COPLANAR_START.Ω,
    COPLANAR_START.ω,
    deg2rad(2) + COPLANAR_START.f
)

COPLANAR_END = Propagators.propagate(Val(:TwoBody), 4500, COPLANAR_END_0)[3].tbd.orbk

NAME = "Paper Coplanar Rendez-vous"

TRANSFER_TIME = 4500.0
