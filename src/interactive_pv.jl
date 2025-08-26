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
# include("casadi_interface.jl")
include("casadi_transfer_model.jl")

function initial_orb_sequence(orb1, tf, Ndisc, nimp::Int, init_coast::Bool, final_coast::Bool, time_partition)
    type_sequence, ncoasts = create_sequence(nimp, init_coast, final_coast)

    coast_times = tf * time_partition
    # cum_time = cumsum([0, coast_times...])

    # prop0 = Propagators.init(Val(:TwoBody), orb1)

    # props = []

    # for i = 1:ncoasts
    #     Propagators.propagate(prop0, cum_time[i])
    #     orbk = prop0.tbd.orbk
    #     push!(props, Propagators.init(Val(:TwoBody), orbk))
    # end

    sequence = []

    last_r, last_v = kepler_to_rv(orb1)
    coast_ind = 0
    for t in type_sequence
        if t == Impulse #impulse in velocity direction
            push!(sequence, Impulse(0.0, last_v/norm(last_v)))
        elseif t == Coast
            coast_ind += 1
            orb = rv_to_kepler(last_r, last_v)
            prop = Propagators.init(Val(:TwoBody), orb)
            tabr = zeros(3, Ndisc)
            tabv = zeros(3, Ndisc)
            for i = 1:Ndisc
                r, v = Propagators.propagate!(prop, coast_times[coast_ind] * (i-1) / (Ndisc-1))
                tabr[:, i] = r
                tabv[:, i] = v
            end
            last_r = tabr[:, end]
            last_v = tabv[:, end]
            push!(sequence, Coast(tabr, tabv, coast_times[coast_ind]))
        end
    end

    sequence
end
##
case_ind = 3
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

N = 100
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
tspan_ppdot = primer_vector(solved_model, 1000)
tspan, ppdot = tspan_ppdot[1]
normp = norm.(getindex.(ppdot, Ref(1:3)))
normpdot = [dot(ppdoti[1:3], ppdoti[4:6]) / norm(ppdoti[1:3]) for ppdoti in ppdot]
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
add_discretized_trajectory!(ax3d, solved_model.sequence[1].rcoast)
add_discretized_trajectory!(ax3d, solved_model.sequence[3].rcoast)
add_discretized_trajectory!(ax3d, solved_model.sequence[5].rcoast)
f