using Setfield
include("sample_orbits.jl")
include("transfer.jl")
include("./primer_vector.jl")
include("integrators.jl")
include("orb_mech.jl")
include("plotting.jl")
include("casadi_transfer_model.jl")
##
case_ind = 4
orb1, orb2 = ORBIT_STARTS[case_ind], ORBIT_ENDS[case_ind]
r1, v1 = kepler_to_rv(orb1)
r2, v2 = kepler_to_rv(orb2)

#vary tf for each orbit
tf_real = orbital_period((orb1.a+orb2.a)/2, GM_EARTH)

L = (orb1.a+orb2.a)/2
T = 1
tfprime = tf_real / T

MUPRIME = GM_EARTH * T ^ 2 / L ^ 3
f = X -> two_body_dyn(X, MUPRIME)


X1 = [r1 / L; v1 * T / L]
X2 = [r2 / L; v2 * T / L]

N = 200
##
planner, transfer = n_impulse_transfer(X1, X2, tfprime, MUPRIME, N, 2, false, false);
##
solver = casadi.nlpsol("S", "ipopt", planner.prob);
## initial guess

seq0 = [scale(s, L, T) 
    for s = initial_orb_sequence(orb1, tf_real, N, 2, false, false, [1.0])
    ]

tab0 = vcat(varlist.(seq0)...)

##
sol = solve_planner(solver, planner, tab0)
##
solved_model = unscale(sol_to_transfer(sol, transfer), L, T)
solved_orb = rv_to_kepler(solved_model.sequence[2].rcoast[:, 1], solved_model.sequence[2].vcoast[:, 1])
solved_prop = Propagators.init(Val(:TwoBody), solved_orb)
##
f, ax3d = plot_orbit(orb1, orb2)
# add_discretized_trajectory!(ax3d, solved_model.sequence[2].rcoast)
add_transfer!(ax3d, solved_model)
f
##
save_with_views!(ax3d, f, "results/$(PREFIXES[case_ind])")
##

#automate this - discard early departure/late arrival
tspan_ppdot = primer_vector(solved_model, 1000)
tspan, ppdot = tspan_ppdot
normp = norm.(eachcol(ppdot[1:3, :]))
normpdot = [dot(ppdoti[1:3], ppdoti[4:6]) / norm(ppdoti[1:3]) for ppdoti in eachcol(ppdot)]
# ##
pv_diag = diagnose_ppdot(normp, normpdot) #remove?
## automate this
f = Figure()
ax1 = Axis(f[1, 1], xlabel = "t (s)", ylabel = "|p|", title="Diagnostic: "*string(diagnose_ppdot(normp, normpdot)))
lines!(ax1, tspan, normp)
vlines!(ax1, tf_real, linestyle=:dash, color=:gray)
ax2 = Axis(f[2, 1], xlabel = "t (s)", ylabel = L"d |p| / dt")
lines!(ax2, tspan, normpdot)
vlines!(ax2, tf_real, linestyle=:dash, color=:gray)
f
##
save("results/"*PREFIXES[case_ind]*"_primer_vector.png", f)
##
# try free impulse time solutions
N = 50
planner, transfer = n_impulse_transfer(X1, X2, tfprime, MUPRIME, N, 2, true, true)

solver = casadi.nlpsol("s", "ipopt", planner.prob)
##
seq0 = [scale(s, L, T) 
    for s = initial_orb_sequence(orb1, tf_real, N, 2, true, true, [0.5, 0, 0.5])
]

tab0 = vcat(varlist.(seq0)...)
##
sol = solve_planner(solver, planner, tab0)
##
solved_model = unscale(sol_to_transfer(sol, transfer), L, T)
# solved_orbs = [rv_to_kepler(c.rcoast[:, 1], c.vcoast[:, 1]) for c in solved_model.sequence if c isa Coast]
# solved_prop = Propagators.init(Val(:TwoBody), solved_orb)
##
f, ax3d = plot_orbit(orb1, orb2)
add_transfer!(ax3d, solved_model, 1e3)
f
##
tspan_ppdot = primer_vector(solved_model, 100)
tspan, ppdot = tspan_ppdot
normp = norm.(eachcol(ppdot[1:3, :]))
plot(tspan, normp)
##
normpdot = [dot(ppdoti[1:3], ppdoti[4:6]) / norm(ppdoti[1:3]) for ppdoti in eachcol(ppdots)]
plot(tspan, normpdot)