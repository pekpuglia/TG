include("../lib.jl")
include("../orbits/hohmann.jl")
##
model = J2model(GM_EARTH, EGM_2008_J2, EARTH_EQUATORIAL_RADIUS)
# model = TwoBodyModel(GM_EARTH)
t = 2orbital_period(HOHMANN_START, GM_EARTH)
r1, v1 = kepler_to_rv(HOHMANN_START)
##
cart_X = final_X(X -> dynamics(X, model), [r1; v1], t, 100, RK8)
cart_orb = rv_to_kepler(cart_X[1:3], cart_X[4:6])
##
j2osc_orb = Propagators.propagate(sat_toolbox_model(model), t, HOHMANN_START)[3].j2oscd.orbk
j2_X = vcat(kepler_to_rv(j2osc_orb)...)
##
display(norm(cart_X[1:3] - j2_X[1:3]))
plot_orbit(HOHMANN_START, j2osc_orb, cart_orb)[1]
##
t = 200orbital_period(HOHMANN_START, GM_EARTH)

cart_traj = trajectory(X -> dynamics(X, model), [r1; v1], t, 10000, RK8)
f, ax = plot_orbit(HOHMANN_START)
add_discretized_trajectory!(ax, cart_traj)
f
## solve example 10.2
orb = KeplerianElements(
    HOHMANN_START.t,
    8059e3,
    0.17136,
    deg2rad(28),
    deg2rad(45),
    deg2rad(30),
    deg2rad(40)
)
t = 48*3600.0
j2osc_orb = Propagators.propagate(sat_toolbox_model(model), t, orb)[3].j2oscd.orbk
##1 point / 3 min
x0 = vcat(kepler_to_rv(orb)...)
xf = final_X(X -> dynamics(X, model), x0, t, 1000, RK8)
cart_orb = rv_to_kepler(xf[1:3], xf[4:6])