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