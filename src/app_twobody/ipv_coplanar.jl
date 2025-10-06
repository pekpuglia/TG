using Setfield
include("transfer.jl")
include("./primer_vector.jl")
include("integrators.jl")
include("orb_mech.jl")
include("plotting.jl")
include("casadi_transfer_model.jl")
##
orb1 = KeplerianElements(
    date_to_jd(2023, 1, 1, 0, 0, 0),
    6748e3,
    0.0,
    5 |> deg2rad,
    35    |> deg2rad,
    0     |> deg2rad,
    0     |> deg2rad
)

orb2_0 = KeplerianElements(
    orb1.t,
    6778e3,
    0.0,
    orb1.i,
    orb1.Ω,
    orb1.ω,
    deg2rad(2) + orb1.f
)

orb2 = Propagators.propagate(Val(:TwoBody), 4500, orb2_0)[3].tbd.orbk

r1, v1 = kepler_to_rv(orb1)
r2, v2 = kepler_to_rv(orb2)

#vary tf for each orbit
tf_real = 4500.0

L = (orb1.a+orb2.a)/2
T = tf_real
orb_model = scale(TwoBodyModel(GM_EARTH), L, T)

tfprime = tf_real / T

X1 = [r1 / L; v1 * T / L]
X2 = [r2 / L; v2 * T / L]

##
N = 200
planner, transfer = n_impulse_transfer(orb_model, X1, X2, tfprime, N, 2, false, false);

solver = casadi.nlpsol("S", "ipopt", planner.prob, Dict("ipopt" => Dict(
    "max_iter" => 3000,
    "constr_viol_tol" => 1e-5)));
## initial guess

seq0 = [scale(s, L, T) 
    for s = initial_orb_sequence(orb1, tf_real, N, 2, false, false, [1.0])
    ]

tab0 = vcat(varlist.(seq0)...)

##
sol = solve_planner(solver, planner, tab0);
##
solved_transfer = unscale(sol_to_transfer(sol, transfer), L, T)
total_dV(solved_transfer)
##
f, ax3d = plot_orbit(orb1, orb2)
# add_discretized_trajectory!(ax3d, solved_transfer.sequence[2].rcoast)
add_transfer!(ax3d, solved_transfer)
f
##
save_with_views!(ax3d, f, "results/$(PREFIXES[case_ind])")
##

#automate this - discard early departure/late arrival
tspan_ppdot = primer_vector(solved_transfer, PVTMGlandorf(), 100)
plot_primer_vector(solved_transfer, tspan_ppdot)[1]
##
save("results/"*PREFIXES[case_ind]*"_primer_vector.png", f)
##
# try free impulse time solutions
N = 100
planner_2, transfer_2 = n_impulse_transfer(orb_model, X1, X2, tfprime, N, 2, true, true)

solver_2 = casadi.nlpsol("s", "ipopt", planner.prob, Dict("ipopt" => Dict(
    "max_iter" => 3000,
    "constr_viol_tol" => 1e-5,
    "max_wall_time" => 30)))
##
seq0 = [scale(s, L, T) 
    for s = initial_orb_sequence(orb1, tf_real, N, 2, true, true, [0.3, 0.3, 0.3])
]

tab0 = vcat(varlist.(seq0)...)
##
sol = solve_planner(solver, planner, tab0)
##
solved_transfer_2 = unscale(sol_to_transfer(sol, transfer), L, T)
total_dV(solved_transfer_2)
##
f, ax3d = plot_orbit(orb1, orb2)
add_transfer!(ax3d, solved_transfer_2, 1e3)
f
##
st2_simpl = rediscretize_transfer(solved_transfer_2, 50)
f, ax3d = plot_orbit(orb1, orb2)
add_transfer!(ax3d, st2_simpl, 1e3)
f
## impulse times
tspan_ppdot = primer_vector(solved_transfer, PVTMGlandorf(), 100)
plot_primer_vector(solved_transfer, tspan_ppdot)[1]
## 3 imp - use last solution
#do j2 case!!!
N = 50
planner_3, transfer_3 = n_impulse_transfer(orb_model, X1, X2, tfprime, N, 3, true, true)

solver_3 = casadi.nlpsol("s", "ipopt", planner_3.prob, Dict("ipopt" => Dict(
    "max_iter" => 3000,
    "constr_viol_tol" => 1e-5,
    "max_wall_time" => 60)))
## add null impulse at optimal position
new_transfer = add_null_impulse(st2_simpl, tspan_ppdot)

new_transfer_scaled = scale(new_transfer, L, T)
tab0 = vcat(varlist.(new_transfer_scaled.sequence)...)
##
f, ax3d = plot_orbit(orb1, orb2)
add_transfer!(ax3d, new_transfer, 1e3)
f
##
sol = solve_planner(solver_3, planner_3, tab0)
##
solved_transfer = unscale(sol_to_transfer(sol, transfer_3), L, T)
total_dV(solved_transfer)
##
f, ax3d = plot_orbit(orb1, orb2)
add_transfer!(ax3d, solved_transfer, 1e3)
f
## impulse times
cumtime = cumsum([0; getfield.(filter(x -> x isa Coast, solved_transfer.sequence), :dt)]) #impulse times function
##
tspan_ppdot = primer_vector(solved_transfer, PVTMGlandorf(), 100)
plot_primer_vector(solved_transfer, tspan_ppdot)[1]