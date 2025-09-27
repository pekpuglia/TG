include("../lib.jl")
include("../orbits/hohmann.jl")
##
r1, v1 = kepler_to_rv(HOHMANN_START)
r2, v2 = kepler_to_rv(HOHMANN_END)


#plot scenario

##
L = (HOHMANN_START.a+HOHMANN_END.a)/2
T = 1
tfprime = TRANSFER_TIME / T

orb_model = scale(J2model(GM_EARTH, EGM_2008_J2, EARTH_EQUATORIAL_RADIUS), L, T)

X1 = [r1 / L; v1 * T / L]
X2 = [r2 / L; v2 * T / L]

N = 100
##
planner, transfer = n_impulse_transfer(orb_model, X1, X2, tfprime, N, 2, false, false);

solver = casadi.nlpsol("S", "ipopt", planner.prob, Dict("ipopt" => Dict(
    "max_iter" => 3000,
    "constr_viol_tol" => 1e-5)));
## initial guess

seq0 = [scale(s, L, T) 
    for s = initial_orb_sequence(HOHMANN_START, TRANSFER_TIME, N, 2, false, false, [1])
    ]

tab0 = vcat(varlist.(seq0)...)

##
sol = solve_planner(solver, planner, tab0)
##
solved_transfer = unscale(sol_to_transfer(sol, transfer), L, T)
##
f, ax3d = plot_orbit(HOHMANN_START, HOHMANN_END)
add_transfer!(ax3d, solved_transfer, 1e4)
f
##
# save_with_views!(ax3d, f, "results/$(PREFIXES[case_ind])")
##
tspan_ppdot = primer_vector(solved_transfer, PVTMFromSTM(100, RK8), 100)
#wrooooong!
# tspan_ppdot = primer_vector(solved_transfer, PVTMGlandorf(), 100)
# tspan_ppdot = primer_vector(solved_transfer, PVTMFromODE(100, RK8), 100)
plot_primer_vector(solved_transfer, tspan_ppdot)[1]
## data summary
total_dV(solved_transfer)
##
impulse_times(solved_transfer)
##impulse magnitudes
impulses(solved_transfer) .|> x -> x.deltaVmag
## 2 impulse free solution
planner, transfer = n_impulse_transfer(orb_model, X1, X2, tfprime, N, 2, true, true);

solver = casadi.nlpsol("S", "ipopt", planner.prob, Dict("ipopt" => Dict(
    "max_iter" => 3000,
    "constr_viol_tol" => 1e-5)));
##
seq0 = [scale(s, L, T) 
    for s = initial_orb_sequence(HOHMANN_START, TRANSFER_TIME, N, 2, true, true, [0.3, 0.3, 0.3])
    ]

tab0 = vcat(varlist.(seq0)...)
sol = solve_planner(solver, planner, tab0)
solved_transfer_2 = unscale(sol_to_transfer(sol, transfer), L, T)
##
total_dV(solved_transfer_2)
##
f, ax3d = plot_orbit(HOHMANN_START, HOHMANN_END)
add_transfer!(ax3d, solved_transfer, 1e4)
f
##
tspan_ppdot_2 = primer_vector(solved_transfer, PVTMFromSTM(100, RK8), 100)
plot_primer_vector(solved_transfer, tspan_ppdot_2)[1]
## 3 impulses
planner, transfer = n_impulse_transfer(orb_model, X1, X2, tfprime, N, 3, true, true);

solver = casadi.nlpsol("S", "ipopt", planner.prob, Dict("ipopt" => Dict(
    "max_iter" => 3000,
    "constr_viol_tol" => 1e-5)));
##
null_transfer_2 = add_null_impulse(solved_transfer_2, tspan_ppdot_2)

tab0 = vcat(varlist.(null_transfer_2.sequence)...)
sol = solve_planner(solver, planner, tab0)
solved_transfer_3 = unscale(sol_to_transfer(sol, transfer), L, T)
##
total_dV(solved_transfer_3)
##
f, ax3d = plot_orbit(HOHMANN_START, HOHMANN_END)
add_transfer!(ax3d, solved_transfer_3, 1e4)
f
##
tspan_ppdot_3 = primer_vector(solved_transfer_3, PVTMFromSTM(100, RK8), 100)
plot_primer_vector(solved_transfer_3, tspan_ppdot_3)[1]
##
impulses(solved_transfer_3) .|> x -> x.deltaVmag
##
impulse_times(solved_transfer_3)
