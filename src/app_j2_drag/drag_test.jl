include("../lib.jl")
# include("../orbits/hohmann.jl")
#check density vs height
#reproduce 10.1
##
rs = EARTH_EQUATORIAL_RADIUS .+ (1:1:1000)*1e3
rho_if = rho_model_curtis.(rs, [R_TABLE_ATM], [RHO_TABLE_ATM], [H_TABLE_ATM])
rho_smooth = rho_model_smooth.(rs, [R_TABLE_ATM], [RHO_TABLE_ATM], [H_TABLE_ATM], 200)
lines(rs, log.(abs.(rho_if .- rho_smooth) ./ rho_if))
## curtis
model = J2DragModel(
    GM_EARTH, 
    EGM_2008_J2*GM_EARTH*EARTH_EQUATORIAL_RADIUS^2, 
    2.2*pi*0.25/100, 
    EARTH_ANGULAR_SPEED,
    R_TABLE_ATM,
    RHO_TABLE_ATM,
    H_TABLE_ATM)
##
orb = KeplerianElements(
    0,
    6955e3,
    0.052,
    deg2rad(65.1),
    deg2rad(340),
    deg2rad(58),
    deg2rad(332)
)

x0 = vcat(kepler_to_rv(orb)...)

days = 108.9
N_per_day = 2000
xs = trajectory(X -> dynamics(X, model, rho_model_curtis), x0, days*24*3600, round(Int, days*N_per_day), RK8)
orbs = [rv_to_kepler(x[1:3], x[4:6]) for x = eachcol(xs)]
has = [o.a * (1 + o.e) for o = orbs] .- EARTH_EQUATORIAL_RADIUS
hps = [o.a * (1 - o.e) for o = orbs] .- EARTH_EQUATORIAL_RADIUS
hps[end]
##
f = lines(1:(days*N_per_day), has)
lines!(f.axis, 1:(days*N_per_day), hps)
f
##
f, ax3d, _ = plot_orbit(orb)
add_discretized_trajectory!(ax3d, xs)
f
##
lines(norm.(eachcol(xs[1:3, :])) .- EARTH_EQUATORIAL_RADIUS)
## scaling doesn't work!!!!!!!!
L = orb.a
T = orbital_period(orb, GM_EARTH)


rho_model_smooth(norm(x0[1:3]), model.r_table, model.rho_table, model.H_table, 50) * L^3

sc_mod = scale(model, L, T)

rho_model_smooth(norm(x0[1:3]) / L, sc_mod.r_table, sc_mod.rho_table, sc_mod.H_table, 1000)


##
scx0 = [x0[1:3] / L; x0[4:6] * T / L]
t = T/5
xf = final_X(X -> dynamics(X, model), x0, t, 10, RK8)
sc_xf = final_X(X -> dynamics(X, sc_mod), scx0, t / T, 10, RK8)
# unsc_xf = [sc_xf[1:3] * L; sc_xf[4:6] * L / T]