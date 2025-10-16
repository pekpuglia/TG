using SatelliteToolboxBase
include("../lib/orb_mech.jl")

PHASING_START = KeplerianElements(
    date_to_jd(2023, 1, 1, 0, 0, 0),
    EARTH_EQUATORIAL_RADIUS + 800e3,
    0.001,
    deg2rad(4),
    deg2rad(0),
    0,
    0
)

PHASING_END = KeplerianElements(
    date_to_jd(2023, 1, 1, 0, 0, 0),
    EARTH_EQUATORIAL_RADIUS + 800e3,
    0.001,
    deg2rad(4),
    deg2rad(0),
    0,
    deg2rad(-30)
)

NAME = "Phasing"

TRANSFER_TIME = 2orbital_period(PHASING_START, GM_EARTH)