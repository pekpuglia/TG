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
#vary tf for each orbit
tf_real = 2*orbital_period(orb2, GM_EARTH) #slightly different to paper
orb_model = TwoBodyModel(GM_EARTH)
##
L = (orb1.a+orb2.a)/2
T = 1

tfprime = tf_real / T

X1 = [r1 / L; v1 * T / L]
X2 = [r2 / L; v2 * T / L]
orbm = scale(orb_model, L, T)

N = 200
planner_2r, transfer_2r = n_impulse_transfer(orbm, X1, X2, tfprime, N, 2, false, false);

solver_2r = casadi.nlpsol("S", "ipopt", planner_2r.prob, Dict("ipopt" => Dict(
    "max_iter" => 3000,
    "constr_viol_tol" => 1e-6)));
## initial guess

seq0 = [scale(s, L, T) 
    for s = initial_orb_sequence(orb1, tf_real, N, 2, false, false, [1.0])
    ]

tab0 = vcat(varlist.(seq0)...)

##
sol = solve_planner(solver_2r, planner_2r, tab0);
##
solved_transfer = unscale(sol_to_transfer(sol, transfer_2r), L, T)
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
tspan_ppdot = primer_vector(solved_transfer, PVTMFromSTM(200, RK8), 100)
plot_primer_vector(solved_transfer, tspan_ppdot)[1]
##
save("results/"*PREFIXES[case_ind]*"_primer_vector.png", f)
##
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
tdv2 = total_dV(solved_transfer_2f)
##
f, ax3d = plot_orbit(orb1, orb2)
add_transfer!(ax3d, solved_transfer_2f, 1e3)
f
## impulse times
tspan_ppdot_2f = primer_vector(solved_transfer_2f, PVTMGlandorf(), 100)
plot_primer_vector(solved_transfer_2f, tspan_ppdot_2f)[1]
## 3 imp - use last solution
#do j2 case!!!
N_3 = 50
planner_3, transfer_3 = n_impulse_transfer(orbm, X1, X2, tfprime, N_3, 3, true, true)

# add_inequality!(planner_3, sum(getfield.(impulses(transfer_3), :deltaVmag)), 0, tdv2 * T / L)

solver_3 = casadi.nlpsol("s", "ipopt", planner_3.prob, Dict("ipopt" => Dict(
    "max_iter" => 3000,
    "constr_viol_tol" => 1e-5,
    "max_wall_time" => 180)))
## add null impulse
seq0 = [scale(s, L, T) 
    for s = initial_orb_sequence(orb1, tf_real, N_3, 3, true, true, [0.1, 0.4, 0.3, 0.2])
]
##
function summarize_transf(t::Transfer)
    newseq = []

    for el in t.sequence
        if el isa Coast
            push!(newseq, Coast(
                [el.rcoast[:, 1] el.rcoast[:, end]],
                [el.vcoast[:, 1] el.vcoast[:, end]],
                el.dt
            )) 
        else
            push!(newseq, el)
        end
    end
    Transfer(t.X1, t.X2, t.model, t.transfer_time, t.nimp, t.init_coast, t.final_coast, newseq)
end

function random_starts(solver, planner::CasADiPlanner, sc_t::Transfer, orb1, tf_real, L, T, n_starts)
    _, ncoasts = create_sequence(sc_t.nimp, sc_t.init_coast, sc_t.final_coast)

    N = size(coasts(sc_t)[1].rcoast)[2]

    min_dV = Inf

    best_partition = nothing
    best_transfer = nothing

    for i = 1:n_starts
        partition = rand(ncoasts)
        partition /= sum(partition)
        seq0 = [scale(s, L, T) 
            for s = initial_orb_sequence(orb1, tf_real, N, sc_t.nimp, sc_t.init_coast, sc_t.final_coast, partition)
        ]
        sol = solve_planner(solver, planner, vcat(varlist.(seq0)...))
        solved_transfer = unscale(sol_to_transfer(sol, sc_t), L, T)
        dv = total_dV(solved_transfer)

        if dv <= min_dV
            min_dV = dv
            best_partition = partition
            best_transfer = summarize_transf(solved_transfer)
        end
    end
    best_transfer, best_partition
end
## init cond from sol 2
transfer_2f_redisc = rediscretize_transfer(solved_transfer_2f, N_3)

null_imp_transf_2 = scale(
    add_null_impulse(
        transfer_2f_redisc, 
        tspan_ppdot_2f), L, T)

seq0 = null_imp_transf_2.sequence
##
tab0 = vcat(varlist.(seq0)...)

sol = solve_planner(solver_3, planner_3, tab0)
##
solved_transfer_3 = unscale(sol_to_transfer(sol, transfer_3), L, T)
total_dV(solved_transfer)
## equal to paper with manual init cond
tspan_ppdot_3 = primer_vector(solved_transfer_3, PVTMGlandorf(), 100)
plot_primer_vector(solved_transfer_3, tspan_ppdot_3)[1]
##
f, ax3d = plot_orbit(orb1, orb2)
add_transfer!(ax3d, solved_transfer, 1e3)
f
##
best_transfer, best_part = random_starts(solver_3, planner_3, transfer_3, orb1, tf_real, L, T, 10)
total_dV(best_transfer)
#impulse times: [  1131.0896958664634
#  8773.458597861762
#  11107.157595175846] give 39.8

## impulse times
tspan_ppdot = primer_vector(best_transfer, PVTMGlandorf(), 100)
plot_primer_vector(best_transfer, tspan_ppdot)[1]
## do 4 impulses - still not equal to paper
N_4 = 50
planner_4, transfer_4 = n_impulse_transfer(orbm, X1, X2, tfprime, N_4, 4, true, true)

# add_inequality!(planner_3, sum(getfield.(impulses(transfer_3), :deltaVmag)), 0, tdv2 * T / L)

solver_4 = casadi.nlpsol("s", "ipopt", planner_4.prob, Dict("ipopt" => Dict(
    "max_iter" => 3000,
    "constr_viol_tol" => 1e-5,
    "max_wall_time" => 180))
)

null_imp_transf_3 = scale(
    add_null_impulse(
        solved_transfer_3, 
        tspan_ppdot_3), L, T)

seq0 = null_imp_transf_3.sequence
##
tab0 = vcat(varlist.(seq0)...)

sol = solve_planner(solver_4, planner_4, tab0)
##
solved_transfer_4 = unscale(sol_to_transfer(sol, transfer_4), L, T)
total_dV(solved_transfer)
##
tspan_ppdot = primer_vector(solved_transfer_4, PVTMGlandorf(), 100)
plot_primer_vector(solved_transfer_4, tspan_ppdot)[1]