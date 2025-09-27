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

NAME = "Hohmann"

TRANSFER_TIME = orbital_period((HOHMANN_START.a+HOHMANN_END.a)/2, GM_EARTH) / 2