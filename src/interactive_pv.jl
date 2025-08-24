# using SatelliteToolboxBase
# using SatelliteToolboxPropagators
# using JuMP
# using Ipopt
# using GLMakie
# using LinearAlgebra
using Setfield
include("sample_orbits.jl")
include("transfer.jl")
include("./primer_vector.jl")
include("integrators.jl")
include("orb_mech.jl")
include("plotting.jl")
# include("jump_base.jl")
include("casadi_interface.jl")

##
case_ind = 1
orb1, orb2 = ORBIT_STARTS[case_ind], ORBIT_ENDS[case_ind]
r1, v1 = kepler_to_rv(orb1)
r2, v2 = kepler_to_rv(orb2)

#vary tf for each orbit
tf_real = orbital_period((orb1.a+orb2.a)/2, GM_EARTH) / 2

L = (orb1.a+orb2.a)/2
T = tf_real
tfprime = tf_real / T

MUPRIME = GM_EARTH * T ^ 2 / L ^ 3
f = X -> two_body_dyn(X, MUPRIME)


X1 = [r1 / L; v1 * T / L]
X2 = [r2 / L; v2 * T / L]

N = 100
tabX = cvar("X", 6, N)
dVmag = cvar("dVmag", 2)
dVdir = cvar("dVdir", 3, 2)

variables = [tabX..., dVmag..., dVdir...]

planner = CasADiPlanner(variables)

add_bounds!(planner, dVmag[1], 0, Inf)
add_bounds!(planner, dVmag[2], 0, Inf)

add_equality!(planner, dVdir[:, 1]' * dVdir[:, 1], 1)
add_equality!(planner, dVdir[:, 2]' * dVdir[:, 2], 1)

#boundary conditions
add_equality!(planner, tabX[1:3, 1], X1[1:3]);
add_equality!(planner, tabX[4:6, 1] - dVmag[1]*dVdir[:, 1], X1[4:6]);
add_equality!(planner, tabX[1:3, end], X2[1:3]);
add_equality!(planner, tabX[4:6, end] + dVmag[2]*dVdir[:, 2], X2[4:6]);

#integration
for i = 1:(N-1)
    add_equality!(planner, tabX[:, i+1] .- RK8(f, tabX[:, i], tfprime/(N-1)), zeros(6));
end

planner.prob["f"] = sum(sx_iterator(dVmag))
##
solver = casadi.nlpsol("S", "ipopt", planner.prob)
## initial guess
tabx0 = []
dVmag0 = []
dVdir0 = []

dVmags0 = [0; 0];
dVdirs0 = [X1[4:6]/norm(X1[4:6]) X2[4:6]/norm(X2[4:6])];
#integrate motion along initial orbit
tabX0 = repeat(X1, 1, N);

for i = 2:N
    tabX0[:, i] = RK8(f, tabX0[:, i-1], tfprime/(N-1));
end

tab0 = [
    reshape(tabX0, 6*N, 1);
    dVmags0;
    reshape(dVdirs0, 3*2, 1)
];
##
sol = solve_planner(solver, planner, tab0)
##
solved_model = unscale(solved(model_transfer), L, T)
solved_orb = rv_to_kepler(solved_model.sequence[2].rcoast[:, 1], solved_model.sequence[2].vcoast[:, 1])
solved_prop = Propagators.init(Val(:TwoBody), solved_orb)
##
f, ax3d = plot_orbit(orb1, orb2)
add_discretized_trajectory!(ax3d, solved_model.sequence[2].rcoast)
f
##
save_with_views!(ax3d, f, "results/$(PREFIXES[case_ind])")
##
# struct PVTrajectory
#     impulse_times::Vector
#     #list of trajectories for each coasting arc
#     #therefore length(p) = nimp - 1
#     p::Vector{Matrix}
#     diagnostic::PRIMER_DIAGNOSTIC
# end



#automate this - discard early departure/late arrival
tspan_ppdot = primer_vector(solved_model, 1000, planar=true)
tspan, ppdot = tspan_ppdot[1]
normp = norm.(getindex.(ppdot, Ref(1:3)))
normpdot = [dot(ppdoti[1:3], ppdoti[4:6]) / norm(ppdoti[1:3]) for ppdoti in ppdot]
# ##
pv_diag = diagnose_ppdot(normp, normpdot) #remove?
## automate this
f = Figure()
ax1 = Axis(f[1, 1], xlabel = "t (s)", ylabel = "|p|", title="Diagnostic: "*string(diagnose_ppdot(normp, normpdot)))
lines!(ax1, tspan, normp)
vlines!(ax1, tf1, linestyle=:dash, color=:gray)
ax2 = Axis(f[2, 1], xlabel = "t (s)", ylabel = L"d |p| / dt")
lines!(ax2, tspan, normpdot)
vlines!(ax2, tf1, linestyle=:dash, color=:gray)
f
##
save("results/"*PREFIXES[case_ind]*"_primer_vector.png", f)
##
# try free impulse time solutions
model = Model(optimizer_with_attributes(Ipopt.Optimizer,
"max_iter" => 3_000,
# "max_wall_time" => 30.0
))

transfer_time = tf1 #1.5 tf1 gives good orbit for case GEO
model, model_transfer = n_impulse_model(model, X1, X2, transfer_time / T, MUPRIME, 100, 2, true, true)

all_r = cat((model_transfer.sequence[i].rcoast for i = [1, 3, 5])..., dims=2)
all_v = cat((model_transfer.sequence[i].vcoast for i = [1, 3, 5])..., dims=2)

initial_guess_initorb!(all_r, all_v, orb1, transfer_time, L, T)

# #false! need to concatenate initial condition segments
# initial_guess_initorb!(model_transfer.sequence[3].rcoast, model_transfer.sequence[3].vcoast, orb1, tf1, L, T)
# initial_guess_initorb!(model_transfer.sequence[5].rcoast, model_transfer.sequence[5].vcoast, orb1, tf1/2  , L, T)

optimize!(model)
##
solved_model = unscale(solved(model_transfer), L, T)
solved_orbs = [rv_to_kepler(c.rcoast[:, 1], c.vcoast[:, 1]) for c in solved_model.sequence if c isa Coast]
# solved_prop = Propagators.init(Val(:TwoBody), solved_orb)
##
f, ax3d = plot_orbit(orb1, orb2)
add_discretized_trajectory!(ax3d, solved_model.sequence[1].rcoast)
add_discretized_trajectory!(ax3d, solved_model.sequence[3].rcoast)
add_discretized_trajectory!(ax3d, solved_model.sequence[5].rcoast)
f