include("../lib.jl")
include("../orbits/hohmann.jl")
##
r1, v1 = kepler_to_rv(HOHMANN_START)
r2, v2 = kepler_to_rv(HOHMANN_END)
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

tolepsilon = 1e-5
solver = casadi.nlpsol("S", "ipopt", planner.prob, Dict("ipopt" => Dict(
    "max_iter" => 3000,
    "constr_viol_tol" => tolepsilon)));
## initial guess

seq0 = [scale(s, L, T) 
    for s = initial_orb_sequence(HOHMANN_START, TRANSFER_TIME, N, 2, false, false, [1])
    ]

tab0 = vcat(varlist.(seq0)...)

##
sol = solve_planner(solver, planner, tab0)
##
solved_transfer2r = unscale(sol_to_transfer(sol, transfer), L, T)
##
f, ax3d, orb_lines = plot_orbit(HOHMANN_START, HOHMANN_END)
imp_scale = 1e3
coast_ps, ia = add_transfer!(ax3d, solved_transfer2r, imp_scale)
Legend(f[1, 2], [orb_lines..., coast_ps[1], ia[1]], ["Initial orbit", "Final orbit", "Coasting arc", "Impulse (×$imp_scale s)"], position = (0.8, 0.9))
f
##
save_with_views("results/j2/hohmann/ICI", HOHMANN_START, HOHMANN_END, solved_transfer2r, imp_scale)
##
tspan_ppdot_stm_2r = primer_vector(solved_transfer2r, PVTMFromSTM(100, RK8), 100)
tspan_ppdot_ode = primer_vector(solved_transfer2r, PVTMFromODE(100, RK8), 100)
##
f, axs = plot_primer_vector(solved_transfer2r, tspan_ppdot_stm_2r, name=NAME, label="STM", style=:dash)
plot_primer_vector!(f, axs[1], axs[2], solved_transfer2r, tspan_ppdot_ode, label="ODE", style=:dashdotdot)
Legend(f[1, 2], axs[1], "Methods")
f
##
save("./results/j2/hohmann/ICI_primer_vector.png", f, px_per_unit = 300/96)
## data summary
hohmann_tab = transfer_table(solved_transfer2r, tspan_ppdot_stm_2r, L, T, tolepsilon, NAME)
dump_table("./report/text/resultados.tex", "% J2 C2C ICI", hohmann_tab)
###########################################################################################################################################
###########################################################################################################################################
###########################################################################################################################################
## 2 impulse free solution
###########################################################################################################################################
###########################################################################################################################################
###########################################################################################################################################

planner, transfer = n_impulse_transfer(orb_model, X1, X2, tfprime, N, 2, true, true);

solver = casadi.nlpsol("S", "ipopt", planner.prob, Dict("ipopt" => Dict(
    "max_iter" => 3000,
    "constr_viol_tol" => tolepsilon)));
##
seq0 = [scale(s, L, T) 
    for s = initial_orb_sequence(HOHMANN_START, TRANSFER_TIME, N, 2, true, true, [0.3, 0.3, 0.3])
    ]

tab0 = vcat(varlist.(seq0)...)
sol = solve_planner(solver, planner, tab0)
solved_transfer2f = unscale(sol_to_transfer(sol, transfer), L, T)
##
f, ax3d, orb_lines = plot_orbit(HOHMANN_START, HOHMANN_END)
imp_scale = 1e4
coast_ps, ia = add_transfer!(ax3d, solved_transfer2f, imp_scale)
Legend(f[1, 2], [orb_lines..., coast_ps[1], ia[1]], ["Initial orbit", "Final orbit", "Coasting arc", "Impulse (×$imp_scale s)"], position = (0.8, 0.9))
f
##
save_with_views("results/j2/hohmann/$(transfer_type(solved_transfer2f))", HOHMANN_START, HOHMANN_END, solved_transfer2f, imp_scale)
##
tspan_ppdot_stm_2f = primer_vector(solved_transfer2f, PVTMFromSTM(100, RK8), 100)
tspan_ppdot_ode = primer_vector(solved_transfer2f, PVTMFromODE(100, RK8), 100)
##
f, axs = plot_primer_vector(solved_transfer2f, tspan_ppdot_stm_2f, name=NAME, label="STM", style=:dash)
plot_primer_vector!(f, axs[1], axs[2], solved_transfer2f, tspan_ppdot_ode, label="ODE", style=:dashdotdot)
Legend(f[1, 2], axs[1], "Methods")
f
##
save("./results/j2/hohmann/$(transfer_type(solved_transfer2f))_primer_vector.png", f, px_per_unit = 300/96)
## data summary
hohmann_tab = transfer_table(solved_transfer2f, tspan_ppdot_stm_2f, L, T, tolepsilon, NAME)
dump_table("./report/text/resultados.tex", "% J2 C2C 2 CICIC", hohmann_tab)
###########################################################################################################################################
###########################################################################################################################################
###########################################################################################################################################
## 3 impulses
###########################################################################################################################################
###########################################################################################################################################
###########################################################################################################################################
planner, transfer = n_impulse_transfer(orb_model, X1, X2, tfprime, N, 3, true, true);

solver = casadi.nlpsol("S", "ipopt", planner.prob, Dict("ipopt" => Dict(
    "max_iter" => 3000,
    "constr_viol_tol" => 1e-5)));
##
null_transfer_2 = add_null_impulse(solved_transfer2f, tspan_ppdot_stm_2f)

tab0 = vcat(varlist.(null_transfer_2.sequence)...)
sol = solve_planner(solver, planner, tab0)
solved_transfer_3 = unscale(sol_to_transfer(sol, transfer), L, T)
##
f, ax3d, orb_lines = plot_orbit(HOHMANN_START, HOHMANN_END)
imp_scale = 1e4
coast_ps, ia = add_transfer!(ax3d, solved_transfer_3, imp_scale)
Legend(f[1, 2], [orb_lines..., coast_ps[1], ia[1]], ["Initial orbit", "Final orbit", "Coasting arc", "Impulse (×$imp_scale s)"], position = (0.8, 0.9))
f
##
save_with_views("results/j2/hohmann/$(transfer_type(solved_transfer_3))", HOHMANN_START, HOHMANN_END, solved_transfer_3, imp_scale)
##
tspan_ppdot_stm_3 = primer_vector(solved_transfer_3, PVTMFromSTM(100, RK8), 100)
tspan_ppdot_ode_3 = primer_vector(solved_transfer_3, PVTMFromODE(100, RK8), 100)
##
f, axs = plot_primer_vector(solved_transfer_3, tspan_ppdot_stm_3, name=NAME, label="STM", style=:dash)
plot_primer_vector!(f, axs[1], axs[2], solved_transfer_3, tspan_ppdot_ode_3, label="ODE", style=:dashdotdot)
Legend(f[1, 2], axs[1], "Methods")
f
##
save("./results/j2/hohmann/$(transfer_type(solved_transfer_3))_primer_vector.png", f, px_per_unit = 300/96)
## data summary
hohmann_tab = transfer_table(solved_transfer_3, tspan_ppdot_stm_3, L, T, tolepsilon, NAME)
dump_table("./report/text/resultados.tex", "% J2 C2C 3 CICICIC", hohmann_tab)
##
sum_tab = summary_table(NAME, solved_transfer2r, solved_transfer2f, solved_transfer_3)
dump_table("./report/text/resultados.tex", "% J2 C2C SUM", sum_tab)
