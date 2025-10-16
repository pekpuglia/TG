using SatelliteToolboxBase
include("../lib/orb_mech.jl")

LEO_MAINTENANCE_START = KeplerianElements(
    date_to_jd(2023, 1, 1, 0, 0, 0),
    EARTH_EQUATORIAL_RADIUS + 400e3,
    0.003,
    deg2rad(53),
    deg2rad(3),
    0,
    0
)

LEO_MAINTENANCE_END = KeplerianElements(
    date_to_jd(2023, 1, 1, 0, 0, 0),
    EARTH_EQUATORIAL_RADIUS + 600e3,
    0.001,
    deg2rad(51),
    deg2rad(0),
    0,
    deg2rad(180)
)

NAME = "LEO Mainenance"

TRANSFER_TIME = orbital_period(LEO_MAINTENANCE_END, GM_EARTH)