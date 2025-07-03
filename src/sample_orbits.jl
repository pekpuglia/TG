using TG
using SatelliteToolboxBase
##
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
    0
)
##
SSO_MAINTENANCE_START = KeplerianElements(
    date_to_jd(2023, 1, 1, 0, 0, 0),
    EARTH_EQUATORIAL_RADIUS + 700e3,
    0.003,
    deg2rad(95),
    deg2rad(5),
    0,
    0
)

SSO_MAINTENANCE_END = KeplerianElements(
    date_to_jd(2023, 1, 1, 0, 0, 0),
    EARTH_EQUATORIAL_RADIUS + 750e3,
    0.001,
    deg2rad(98),
    deg2rad(0),
    0,
    0
)
##
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
##
GEO_INSERTION_START, GEO_INSERTION_END = let 
    a0 = EARTH_EQUATORIAL_RADIUS + 600e3
    af = 42164e3
    egto = (af - a0) / (af + a0)
    orb1 = KeplerianElements(
        date_to_jd(2023, 1, 1, 0, 0, 0),
        (a0+af)/2,
        egto,
        deg2rad(4),
        deg2rad(0),
        0,
        0
    )

    orb2 = KeplerianElements(
        date_to_jd(2023, 1, 1, 0, 0, 0),
        af,
        0,
        deg2rad(0),
        deg2rad(0),
        0,
        deg2rad(0)
    )
    orb1, orb2
end