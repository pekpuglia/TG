using SatelliteToolboxBase
##
HOHMANN_START = KeplerianElements(
    date_to_jd(2023, 1, 1, 0, 0, 0),
    7000e3,
    0.0,
    51 |> deg2rad,
    0    |> deg2rad,
    0     |> deg2rad,
    0     |> deg2rad
)

HOHMANN_END = KeplerianElements(
    HOHMANN_START.t,
    9000e3,
    0.0,
    HOHMANN_START.i,
    HOHMANN_START.Ω,
    HOHMANN_START.ω,
    deg2rad(180) + HOHMANN_START.f
)
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
    deg2rad(180)
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
    HOHMANN_START = KeplerianElements(
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
    HOHMANN_START, orb2
end

COPLANAR_RDV_START = KeplerianElements(
    date_to_jd(2023, 1, 1, 0, 0, 0),
    6748e3,
    0.0,
    5 |> deg2rad,
    35    |> deg2rad,
    0     |> deg2rad,
    0     |> deg2rad
)

COPLANAR_RDV_END = KeplerianElements(
    HOHMANN_START.t,
    6778e3,
    0.0,
    HOHMANN_START.i,
    HOHMANN_START.Ω,
    HOHMANN_START.ω,
    deg2rad(2) + HOHMANN_START.f
)

ORBIT_STARTS = [
    HOHMANN_START,
    LEO_MAINTENANCE_START,
    SSO_MAINTENANCE_START,
    PHASING_START,
    GEO_INSERTION_START
]

ORBIT_ENDS = [
    HOHMANN_END,
    LEO_MAINTENANCE_END,
    SSO_MAINTENANCE_END,
    PHASING_END,
    GEO_INSERTION_END
]

#add transfer times

PREFIXES = [
    "hohmann"
    "leo_maintenance"
    "sso_maintenance"
    "phasing"
    "geo_insertion"
]