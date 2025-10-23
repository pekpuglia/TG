using Setfield
include("../lib.jl")
##
NAME = "Non-Coplanar rendez-vous"
orb1 = KeplerianElements(
    date_to_jd(2023, 1, 1, 0, 0, 0),
    6748.1e3,
    0.0,
    42.1 |> deg2rad,
    120.2    |> deg2rad,
    0     |> deg2rad,
    175     |> deg2rad
)

orb2_0 = KeplerianElements(
    orb1.t,
    6778.1e3,
    0.0,
    42 |> deg2rad,
    120 |> deg2rad,
    0,
    180 |> deg2rad
)

orb2 = Propagators.propagate(Val(:TwoBody), 2*orbital_period(orb2_0, GM_EARTH), orb2_0)[3].tbd.orbk

r1, v1 = kepler_to_rv(orb1)
x1 = [r1; v1]
r2, v2 = kepler_to_rv(orb2)
x2 = [r2; v2]
##
fig, ax3d, orb_lines = plot_orbit(orb1, orb2)
axislegend(ax3d, orb_lines, ["Initial orbit", "Final orbit"], position = (0.8, 0.9))
fig
##
save_with_views!(ax3d, fig, "./results/two_body/ipv_noncop/scenario")
##
##vary tf for each orbit
tf_real = 2*orbital_period(orb2, GM_EARTH) #slightly different to paper
orb_model = TwoBodyModel(GM_EARTH)
##
ncop_table = setup_table(orb1, orb2, tf_real)
dump_table("./report/text/resultados.tex", "% NCOP SCENARIO", ncop_table)
##
L = (orb1.a+orb2.a)/2
T = 1

tfprime = tf_real / T

X1 = [r1 / L; v1 * T / L]
X2 = [r2 / L; v2 * T / L]
orbm = scale(orb_model, L, T)

N = 200
planner_2r, transfer_2r = n_impulse_transfer(orbm, X1, X2, tfprime, N, 2, false, false);

tolepsilon = 1e-6
solver_2r = casadi.nlpsol("S", "ipopt", planner_2r.prob, Dict("ipopt" => Dict(
    "max_iter" => 3000,
    "constr_viol_tol" => tolepsilon)));
## initial guess
seq0 = [scale(s, L, T) 
    for s = initial_orb_sequence(orb1, tf_real, N, 2, false, false, [1.0])
    ]

tab0 = vcat(varlist.(seq0)...)
##
sol = solve_planner(solver_2r, planner_2r, tab0);
##
solved_transfer_2r = unscale(sol_to_transfer(sol, transfer_2r), L, T)
##
f, ax3d, orb_lines = plot_orbit(orb1, orb2)
imp_sc = 1e3
coast_ps, ia = add_transfer!(ax3d, solved_transfer_2r, imp_sc)
Legend(f[1, 2], [orb_lines..., coast_ps[1], ia[1]], ["Initial orbit", "Final orbit", "Coasting arc", "Impulse x$imp_sc s"], position = (0.8, 0.9))
f
##
save_with_views!(ax3d, f, "results/two_body/ipv_noncop/$(transfer_type(solved_transfer_2r))")
##
#automate this - discard early departure/late arrival
tspan_ppdot_glandorf = primer_vector(solved_transfer_2r, PVTMGlandorf(), 100)
tspan_ppdot_stm      = primer_vector(solved_transfer_2r, PVTMFromSTM(100, RK8), 100)
tspan_ppdot_ode      = primer_vector(solved_transfer_2r, PVTMFromODE(100, RK8), 100)
##
f, axs = plot_primer_vector(solved_transfer_2r, tspan_ppdot_glandorf, name=NAME, label="Glandorf", style=:dash)
plot_primer_vector!(f, axs[1], axs[2], solved_transfer_2r, tspan_ppdot_stm, label="STM", style=:dashdot)
plot_primer_vector!(f, axs[1], axs[2], solved_transfer_2r, tspan_ppdot_ode, label="ODE", style=:dashdotdot)
Legend(f[1, 2], axs[1], "Methods")
f
##
save("./results/two_body/ipv_noncop/$(transfer_type(solved_transfer_2r))_primer_vector.png", f, px_per_unit = 300/96)
##
tableICI = transfer_table(solved_transfer_2r, tspan_ppdot_glandorf, L, T, tolepsilon, NAME)
dump_table("./report/text/resultados.tex", "% TB NCOP ICI", tableICI)
###########################################################################################################################################
###########################################################################################################################################
###########################################################################################################################################
## 2 impulse free solution
###########################################################################################################################################
###########################################################################################################################################
###########################################################################################################################################
# try free impulse time solutions
L = (orb1.a+orb2.a)/2
T = tf_real

tfprime = tf_real / T

X1 = [r1 / L; v1 * T / L]
X2 = [r2 / L; v2 * T / L]
orbm = scale(orb_model, L, T)

N_2 = 100
planner_2f, transfer_2f = n_impulse_transfer(orbm, X1, X2, tfprime, N_2, 2, true, true)

solver_2f = casadi.nlpsol("s", "ipopt", planner_2f.prob, Dict("ipopt" => Dict(
    "max_iter" => 3000,
    "constr_viol_tol" => 1e-5,
    "max_wall_time" => 60)))
## 0.1 0.1 0.8 - 67
# 0.1 0.8 0.1 - 55
seq0 = [scale(s, L, T) 
    for s = initial_orb_sequence(orb1, tf_real, N_2, 2, true, true, [0.7, 0.1, 0.2])
]

tab0 = vcat(varlist.(seq0)...)
##
sol = solve_planner(solver_2f, planner_2f, tab0)

solved_transfer_2f = unscale(sol_to_transfer(sol, transfer_2f), L, T)
##
f, ax3d, orb_lines = plot_orbit(orb1, orb2)
imp_sc = 5e5
coast_ps, ia = add_transfer!(ax3d, solved_transfer_2f, imp_sc)
Legend(f[1, 2], [orb_lines..., coast_ps[1], ia[1]], ["Initial orbit", "Final orbit", "Coasting arc", "Impulse ×" * (@sprintf "%.1e" imp_sc) * " s"], position = (0.8, 0.9))
f
##
save_with_views!(ax3d, f, "results/two_body/ipv_noncop/$(transfer_type(solved_transfer_2f))")
##
tspan_ppdot_glandorf_2f = primer_vector(solved_transfer_2f, PVTMGlandorf(), 100)
tspan_ppdot_stm      = primer_vector(solved_transfer_2f, PVTMFromSTM(100, RK8), 100)
tspan_ppdot_ode      = primer_vector(solved_transfer_2f, PVTMFromODE(100, RK8), 100)
##
f, axs = plot_primer_vector(solved_transfer_2f, tspan_ppdot_glandorf_2f, name=NAME, label="Glandorf", style=:dash)
plot_primer_vector!(f, axs[1], axs[2], solved_transfer_2f, tspan_ppdot_stm, label="STM", style=:dashdot)
plot_primer_vector!(f, axs[1], axs[2], solved_transfer_2f, tspan_ppdot_ode, label="ODE", style=:dashdotdot)
Legend(f[1, 2], axs[1], "Methods")
f
##
save("./results/two_body/ipv_noncop/$(transfer_type(solved_transfer_2f))_primer_vector.png", f, px_per_unit = 300/96)
##
tableCICIC = transfer_table(solved_transfer_2f, tspan_ppdot_glandorf_2f, L, T, tolepsilon, NAME)
dump_table("./report/text/resultados.tex", "% TB NCOP 2 CICIC", tableCICIC)
###########################################################################################################################################
###########################################################################################################################################
###########################################################################################################################################
## 3 imp - use last solution
###########################################################################################################################################
###########################################################################################################################################
###########################################################################################################################################
N_3 = 50
planner_3, transfer_3 = n_impulse_transfer(orbm, X1, X2, tfprime, N_3, 3, true, true)

tolepsilon = 1e-5
solver_3 = casadi.nlpsol("s", "ipopt", planner_3.prob, Dict("ipopt" => Dict(
    "max_iter" => 3000,
    "constr_viol_tol" => tolepsilon,
    "max_wall_time" => 30,
    "max_resto_iter" => 100)))
## 
# best_random_transfer, best_part = random_starts(solver_3, planner_3, transfer_3, orb1, tf_real, L, T, 3)

# total_dV(best_random_transfer)
##
seq0 = [scale(s, L, T) 
    for s = initial_orb_sequence(orb1, tf_real, N_3, 3, true, true, [ 0.22790549731708265
        0.3636515390073254
        0.2337753328534541
        0.17466763082213788])
]

tab0 = vcat(varlist.(seq0)...)

sol = solve_planner(solver_3, planner_3, tab0)
transfer3_manual = unscale(sol_to_transfer(sol, transfer_3), L, T)
total_dV(transfer3_manual)
##
solved_transfer_3 = transfer3_manual
##
f, ax3d, orb_lines = plot_orbit(orb1, orb2)
imp_sc = 5e5
coast_ps, ia = add_transfer!(ax3d, solved_transfer_3, imp_sc)
Legend(f[1, 2], [orb_lines..., coast_ps[1], ia[1]], ["Initial orbit", "Final orbit", "Coasting arc", "Impulse ×" * (@sprintf "%.1e" imp_sc) * " s"], position = (0.8, 0.9))
f
##
save_with_views!(ax3d, f, "results/two_body/ipv_noncop/$(transfer_type(solved_transfer_3))")
##
tspan_ppdot_glandorf_3 = primer_vector(solved_transfer_3, PVTMGlandorf(), 100)
tspan_ppdot_stm      = primer_vector(solved_transfer_3, PVTMFromSTM(100, RK8), 100)
tspan_ppdot_ode      = primer_vector(solved_transfer_3, PVTMFromODE(100, RK8), 100)
##
f, axs = plot_primer_vector(solved_transfer_3, tspan_ppdot_glandorf_3, name=NAME, label="Glandorf", style=:dash)
plot_primer_vector!(f, axs[1], axs[2], solved_transfer_3, tspan_ppdot_stm, label="STM", style=:dashdot)
plot_primer_vector!(f, axs[1], axs[2], solved_transfer_3, tspan_ppdot_ode, label="ODE", style=:dashdotdot)
Legend(f[1, 2], axs[1], "Methods")
f
##
save("./results/two_body/ipv_noncop/$(transfer_type(solved_transfer_3))_primer_vector.png", f, px_per_unit = 300/96)
##
tableCICICIC = transfer_table(solved_transfer_3, tspan_ppdot_glandorf_3, L, T, tolepsilon, NAME)
dump_table("./report/text/resultados.tex", "% TB NCOP 3 CICICIC", tableCICICIC)
##
###########################################################################################################################################
###########################################################################################################################################
###########################################################################################################################################
## 4 impulses
###########################################################################################################################################
###########################################################################################################################################
###########################################################################################################################################
N_4 = 50
planner_4, transfer_4 = n_impulse_transfer(orbm, X1, X2, tfprime, N_4, 4, true, true)

solver_4 = casadi.nlpsol("s", "ipopt", planner_4.prob, Dict("ipopt" => Dict(
    "max_iter" => 3000,
    "constr_viol_tol" => 1e-5,
    "max_wall_time" => 30,
    "max_resto_iter" => 100))
)
##
# best_random_transfer, best_part = random_starts(solver_4, planner_4, transfer_4, orb1, tf_real, L, T, 20)
# total_dV(best_random_transfer)
##
seq0 = [scale(s, L, T) 
    for s = initial_orb_sequence(orb1, tf_real, N_4, 4, true, true, [0.07205383292143142, 0.8437961343872056, 0.03566546650471125, 0.04621096008052483, 0.002273606106126924])
]

tab0 = vcat(varlist.(seq0)...)

sol = solve_planner(solver_4, planner_4, tab0)
transfer4_manual = unscale(sol_to_transfer(sol, transfer_4), L, T)
total_dV(transfer4_manual)
##
solved_transfer_4 = transfer4_manual
##
f, ax3d, orb_lines = plot_orbit(orb1, orb2)
imp_sc = 5e5
coast_ps, ia = add_transfer!(ax3d, solved_transfer_4, imp_sc)
Legend(f[1, 2], [orb_lines..., coast_ps[1], ia[1]], ["Initial orbit", "Final orbit", "Coasting arc", "Impulse ×" * (@sprintf "%.1e" imp_sc) * " s"], position = (0.8, 0.9))
f
##
save_with_views!(ax3d, f, "results/two_body/ipv_noncop/$(transfer_type(solved_transfer_4))")
##
tspan_ppdot_glandorf_4 = primer_vector(solved_transfer_4, PVTMGlandorf(), 100)
tspan_ppdot_stm      = primer_vector(solved_transfer_4, PVTMFromSTM(100, RK8), 100)
tspan_ppdot_ode      = primer_vector(solved_transfer_4, PVTMFromODE(100, RK8), 100)
##
f, axs = plot_primer_vector(solved_transfer_4, tspan_ppdot_glandorf_4, name=NAME, label="Glandorf", style=:dash)
plot_primer_vector!(f, axs[1], axs[2], solved_transfer_4, tspan_ppdot_stm, label="STM", style=:dashdot)
plot_primer_vector!(f, axs[1], axs[2], solved_transfer_4, tspan_ppdot_ode, label="ODE", style=:dashdotdot)
Legend(f[1, 2], axs[1], "Methods")
f
##
save("./results/two_body/ipv_noncop/$(transfer_type(solved_transfer_4))_primer_vector.png", f, px_per_unit = 300/96)
##
tableCICICICIC = transfer_table(solved_transfer_4, tspan_ppdot_glandorf_4, L, T, tolepsilon, NAME)
dump_table("./report/text/resultados.tex", "% TB NCOP 4 CICICICIC", tableCICICICIC)
##
sum_tab = summary_table(NAME, solved_transfer_2r, solved_transfer_2f, solved_transfer_3, solved_transfer_4)
dump_table("./report/text/resultados.tex", "% TB NCOP SUM", sum_tab)
