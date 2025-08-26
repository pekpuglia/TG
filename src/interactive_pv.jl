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

function primer_vector(transfer::Transfer, npoints; tpbvp_kwargs...)
    coast_list = filter(x -> x isa Coast, transfer.sequence)
    impulse_list = filter(x-> x isa Impulse, transfer.sequence)

    dts = getfield.(coast_list, :dt)
    
    impulse_times = cumsum([0; dts])

    #compute primer vector on coasts surrounded by impulses
    two_impulse_coasts = []
    for (e1, e2, e3) in zip(transfer.sequence[1:end-2], transfer.sequence[2:end-1], transfer.sequence[3:end])
        if e1 isa Impulse && e2 isa Coast && e3 isa Impulse
            push!(two_impulse_coasts, (e1, e2, e3))
        end
    end

    # if transfer.sequence[1] isa Coast || transfer.sequence[end] isa Coast
    #     @warn "Unimplemented edge coast case!!!!"
    # end

    tspan_ppdot = []

    for tic in two_impulse_coasts
        i1, c, i2 = tic
        dv1 = i1.deltaVmag * i1.deltaVdir
        dv2 = i2.deltaVmag * i2.deltaVdir

        orbit = rv_to_kepler(c.rcoast[:, 1], c.vcoast[:, 1])

        propagator = Propagators.init(Val(:TwoBody), orbit)

        push!(tspan_ppdot, ppdot_deltavs(propagator, dv1, dv2, c.dt, npoints; tpbvp_kwargs...))
    end

    if transfer.sequence[1] isa Coast
        #ppdot at the first impulse
        #[first coast between 2 impulses][2nd element in (tspan, ppdot)][first ppdot in the coast]
        ppdot_end = tspan_ppdot[1][2][1]

        first_coast_propagator = Propagators.init(Val(:TwoBody), rv_to_kepler(transfer.X1[1:3], transfer.X1[4:6]))
        first_coast_duration = transfer.sequence[1].dt

        #backwards propagation
        tspan = range(0, first_coast_duration, npoints)
        ppdot = []
        for t in tspan
            Phi = TG.Phi_time(first_coast_propagator, t-first_coast_duration)
            push!(ppdot, Phi*ppdot_end)
        end
        tspan_ppdot = [(tspan, ppdot), tspan_ppdot...]
    end

    if transfer.sequence[end] isa Coast
        #ppdot at the last impulse
        #[last coast between 2 impulses][2nd element in (tspan, ppdot)][last ppdot in the coast]
        ppdot_start = tspan_ppdot[end][2][end]

        last_coast_propagator = Propagators.init(Val(:TwoBody), rv_to_kepler(transfer.sequence[end].rcoast[:, 1], transfer.sequence[end].vcoast[:, 1]))
        last_coast_duration = transfer.sequence[end].dt

        #backwards propagation
        tspan = range(0, last_coast_duration, npoints)
        ppdot = []
        for t in tspan
            Phi = TG.Phi_time(last_coast_propagator, t)
            push!(ppdot, Phi*ppdot_start)
        end
        push!(tspan_ppdot, (tspan, ppdot))
    end

    # #merge everything
    # tspans = first.(tspan_ppdot)
    # prev_time = [0, last.(tspans)]
    # tspan = []
    # #accumulate tspan?
    tspan_ppdot
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
function add_transfer!(ax3d, solved_transfer::Transfer, scaling=1e3)
    #get first position
    last_r = if solved_transfer.sequence[1] isa Coast
        solved_transfer.sequence[1].rcoast[:, 1]
    else
        solved_transfer.sequence[2].rcoast[:, 1]
    end

    for el in solved_transfer.sequence
        if el isa Coast
            add_discretized_trajectory!(ax3d, el.rcoast)
            last_r = el.rcoast[:, end]
        else
            scatter!(ax3d, last_r..., color="red")
            #add line on impulse
            dir = el.deltaVdir
            mag = el.deltaVmag

            arrow_data = [last_r last_r+mag*dir*scaling]
            lines!(ax3d, arrow_data[1, :], arrow_data[2, :], arrow_data[3, :], color="red")
        end
    end
end

f, ax3d = plot_orbit(orb1, orb2)
# add_discretized_trajectory!(ax3d, solved_model.sequence[2].rcoast)
add_transfer!(ax3d, solved_model)
f
##
save_with_views!(ax3d, f, "results/$(PREFIXES[case_ind])")
##

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
add_transfer!(ax3d, solved_model, 1e3)
f
##
tspan_ppdot = primer_vector(solved_model, 100)
tspans = first.(tspan_ppdot)
prevtime = cumsum([0; last.(tspans[1:end-1])])
tspan = vcat((ts .+ pv for (ts, pv) in zip(tspans, prevtime))...)
ppdots = hcat(vcat(last.(tspan_ppdot)...)...)
normp = norm.(eachcol(ppdots[1:3, :]))
plot(tspan, normp)
##
normpdot = [dot(ppdoti[1:3], ppdoti[4:6]) / norm(ppdoti[1:3]) for ppdoti in eachcol(ppdots)]
plot(tspan, normpdot)