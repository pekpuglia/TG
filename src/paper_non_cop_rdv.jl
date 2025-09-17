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
    6741.1e3,
    0.0,
    42.1 |> deg2rad,
    120.2    |> deg2rad,
    0     |> deg2rad,
    175     |> deg2rad
)

orb2_0 = KeplerianElements(
    orb1.t,
    6771.1e3,
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

N = 100
planner_2f, transfer_2f = n_impulse_transfer(orbm, X1, X2, tfprime, N, 2, true, true)

solver_2f = casadi.nlpsol("s", "ipopt", planner_2f.prob, Dict("ipopt" => Dict(
    "max_iter" => 3000,
    "constr_viol_tol" => 1e-5,
    "max_wall_time" => 60)))
##
seq0 = [scale(s, L, T) 
    for s = initial_orb_sequence(orb1, tf_real, N, 2, true, true, [0.8, 0.1, 0.1])
]

tab0 = vcat(varlist.(seq0)...)
##
sol = solve_planner(solver_2f, planner_2f, tab0)

solved_transfer_2f = unscale(sol_to_transfer(sol, transfer_2f), L, T)
total_dV(solved_transfer_2f)
##
f, ax3d = plot_orbit(orb1, orb2)
add_transfer!(ax3d, solved_transfer, 1e3)
f
## impulse times
tspan_ppdot = primer_vector(solved_transfer, PVTMGlandorf(), 100)
plot_primer_vector(solved_transfer, tspan_ppdot)[1]
## 3 imp - use last solution
#do j2 case!!!
N = 50
planner_3, transfer_3 = n_impulse_transfer(orbm, X1, X2, tfprime, N, 3, true, true)

solver_3 = casadi.nlpsol("s", "ipopt", planner_3.prob, Dict("ipopt" => Dict(
    "max_iter" => 3000,
    "constr_viol_tol" => 1e-5,
    "max_wall_time" => 60)))
## add null impulse
seq0 = [scale(s, L, T) 
    for s = initial_orb_sequence(orb1, tf_real, N, 3, true, true, [0.5, 0.1, 0.2, 0.1])
]

tab0 = vcat(varlist.(seq0)...)
##
N = 50 #change N later
tspan, ppdot = tspan_ppdot

# normpdot = [dot(ppdoti[1:3], ppdoti[4:6]) / norm(ppdoti[1:3]) for ppdoti in eachcol(ppdot)]

normp = norm.(eachcol(ppdot[1:3, :]))

max_norm_time_ind = findall(
    (prev_el_next) -> prev_el_next[1] <= prev_el_next[2] && prev_el_next[2] >= prev_el_next[3], 
    collect(zip(normp[1:end-2], normp[2:end-1], normp[3:end]))) .+ 1

max_norm_time_ind = max_norm_time_ind[findmax(i -> normp[i], max_norm_time_ind)[2]]

max_norm_time = tspan[max_norm_time_ind]
#find where to insert impulse
imp_ts = impulse_times(solved_transfer)
#last impulse before maxnorm
#if = 0, new impulse comes at the beginning
new_impulse_ind = something(findlast(<=(max_norm_time), imp_ts), 0)

#coast to split index
impulse_indices = findall(el -> el isa Impulse, solved_transfer.sequence)

coast_to_split_sequence_index = (new_impulse_ind == 0) ? 1 : (impulse_indices.+1)[new_impulse_ind]

coast_to_split = solved_transfer.sequence[coast_to_split_sequence_index]

delta_t_split = max_norm_time - [0; imp_ts][new_impulse_ind+1]

xcoast_before = zeros(6, N)

xcoast_before[:, 1] = [coast_to_split.rcoast[:, 1]; coast_to_split.vcoast[:, 1]]

for i = 2:N
    xcoast_before[:, i] =  RK8(X -> dynamics(X, solved_transfer.model), xcoast_before[:, i-1], delta_t_split / (N-1))
end

delta_t_after = coast_to_split.dt - delta_t_split

xcoast_after = zeros(6, N)

xcoast_after[:, 1] = xcoast_before[:, end]

for i = 2:N
    xcoast_after[:, i] =  RK8(X -> dynamics(X, solved_transfer.model), xcoast_after[:, i-1], delta_t_after / (N-1))
end

coast_before = Coast(xcoast_before[1:3, :], xcoast_before[4:6, :], delta_t_split)
coast_after = Coast(xcoast_after[1:3, :], xcoast_after[4:6, :], delta_t_after)


new_impulse = Impulse(0.0, (x -> x / norm(x))(ppdot[1:3, max_norm_time_ind]))

new_seq = [solved_transfer.sequence[1:coast_to_split_sequence_index-1]; coast_before; new_impulse; coast_after; solved_transfer.sequence[coast_to_split_sequence_index+1:end]]
new_transfer = Transfer(solved_transfer.X1, solved_transfer.X2, solved_transfer.model, solved_transfer.transfer_time, new_seq)

new_transfer_scaled = scale(new_transfer, L, T)
tab0 = vcat(varlist.(new_transfer_scaled.sequence)...)
##
f, ax3d = plot_orbit(orb1, orb2)
add_transfer!(ax3d, new_transfer, 1e3)
f
##
sol = solve_planner(solver_3, planner_3, tab0)

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