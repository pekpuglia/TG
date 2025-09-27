using SatelliteToolboxBase
##

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