include("../lib.jl")
# include("../orbits/hohmann.jl")
## curtis
#km
h = 600.0
#kg/m^3
rho = 1.137e-13
H = 76.377 #km

DRAG_REF_RADIUS = EARTH_EQUATORIAL_RADIUS + 1000h
DRAG_REF_RHO = 1.137e-13
DRAG_REF_H = 76.377e3

model = J2DragModel(
    GM_EARTH, 
    EGM_2008_J2*GM_EARTH*EARTH_EQUATORIAL_RADIUS^2, 
    EARTH_ANGULAR_SPEED, 
    2.2*pi*0.25/100,
    DRAG_REF_RADIUS,
    DRAG_REF_RHO,
    DRAG_REF_H)
##
orb = KeplerianElements(
    0,
    6955e3,
    0.052049,
    deg2rad(65.1),
    deg2rad(340),
    deg2rad(58),
    deg2rad(332)
)

x0 = vcat(kepler_to_rv(orb)...)

days = 200
xf = final_X(X -> dynamics(X, model), x0, days*24*3600, days*1000, RK8)
norm(xf[1:3]) - EARTH_EQUATORIAL_RADIUS - 100e3