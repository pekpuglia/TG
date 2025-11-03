include("../lib.jl")
include("../orbits/hohmann.jl")
##
r1, v1 = kepler_to_rv(HOHMANN_START)
r2, v2 = kepler_to_rv(HOHMANN_END)

##plot scenario
fig, ax3d, orb_lines = plot_orbit(HOHMANN_START, HOHMANN_END)
axislegend(ax3d, orb_lines, ["Initial orbit", "Final orbit"], position = (0.8, 0.9))
fig
##
save("./results/two_body/hohmann/scenario.png", fig, px_per_unit = 300/96)
##
L = (HOHMANN_START.a+HOHMANN_END.a)/2
T = 1
tfprime = TRANSFER_TIME / T

orb_model = scale(TwoBodyModel(GM_EARTH), L, T)

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
solved_transfer = unscale(sol_to_transfer(sol, transfer), L, T)
##
plot_transfer(HOHMANN_START, HOHMANN_END, solved_transfer, 1e4)
##
save_with_views!("results/two_body/hohmann/ICI", HOHMANN_START, HOHMANN_END, solved_transfer, 1e4)
##
tspan_ppdot_glandorf = primer_vector(solved_transfer, PVTMGlandorf(), 100)
tspan_ppdot_stm = primer_vector(solved_transfer, PVTMFromSTM(100, RK8), 100)
tspan_ppdot_ode = primer_vector(solved_transfer, PVTMFromODE(100, RK8), 100)
##
f, axs = plot_primer_vector(solved_transfer, tspan_ppdot_glandorf, name=NAME, label="Glandorf", style=:dash)
plot_primer_vector!(f, axs[1], axs[2], solved_transfer, tspan_ppdot_stm, label="STM", style=:dashdot)
plot_primer_vector!(f, axs[1], axs[2], solved_transfer, tspan_ppdot_ode, label="ODE", style=:dashdotdot)
Legend(f[1, 2], axs[1], "Methods")
f
##
save("./results/two_body/hohmann/ICI_primer_vector.png", f, px_per_unit = 300/96)
## data summary
#L T constr_viol_tol deltax_tol  final error
#maxnormp diag
#i = 1 t dV theta_v theta_p, 
hohmann_tab = transfer_table(solved_transfer, tspan_ppdot_glandorf, L, T, tolepsilon, NAME)
dump_table("./report/text/resultados.tex", "% TB C2C ICI", hohmann_tab)
##
hohmann_setup = setup_table(HOHMANN_START, HOHMANN_END, TRANSFER_TIME)
dump_table("./report/text/resultados.tex", "% C2C SETUP", hohmann_setup)