using SatelliteToolboxBase
using SatelliteToolboxPropagators
using JuMP
using Ipopt
include("TG.jl")
using .TG
using GLMakie
using LinearAlgebra
using Setfield
include("sample_orbits.jl")
##
function initial_guess_initorb!(orb1, tf, r, v)
    prop = Propagators.init(Val(:TwoBody), orb1)

    N = size(r)[2]

    for i = 1:N
        rguess, vguess = Propagators.propagate!(prop, (i-1)*tf/(N-1))
        set_start_value.(r[:, i], rguess)
        set_start_value.(v[:, i], vguess)
    end
end

struct Impulse
    deltaVmag
    deltaVdir
end
solved(i::Impulse) = Impulse(value(i.deltaVmag), value.(i.deltaVdir))
unscale(i::Impulse, L, T) = Impulse(L/T * i.deltaVmag, i.deltaVdir)

struct Coast
    rcoast::Matrix
    vcoast::Matrix
    dt
end
solved(c::Coast) = Coast(value.(c.rcoast), value.(c.vcoast), value(c.dt))
unscale(c::Coast, L, T) = Coast(L * c.rcoast, L/T * c.vcoast, T * c.dt)

struct Transfer
    X1::Vector
    X2::Vector
    mu::Float64
    transfer_time
    # rcoast::Array
    # vcoast::Array
    # deltaVmags::Vector
    # deltaVdirs::Matrix
    # dts::Vector
    sequence::Vector{Union{Impulse, Coast}}
end

#scaling agnostic!
function n_impulse_model(model, X1, X2, tf, mu, Ndisc, nimp::Int, init_coast::Bool, final_coast::Bool)
    ncoasts = nimp - 1 + init_coast + final_coast

    dts = @variable(model, [1:ncoasts], base_name = "dt", lower_bound=0, upper_bound=tf)
    @constraint(model, sum(dts) == tf)

    deltaVmags = @variable(model, [1:nimp], lower_bound = 0, base_name = "dVmag")
    deltaVdirs = @variable(model, [1:3, 1:nimp], base_name = "dVdir")
    @constraint(model, [i=1:nimp], deltaVdirs[:, i]' * deltaVdirs[:, i] == 1)

    impulses = Impulse.(deltaVmags, eachcol(deltaVdirs))

    rvcoasts = [add_coast_segment(
    model, dts[i], Ndisc, i,
    dyn=(X -> two_body_dyn(X, mu))) for i = 1:ncoasts]

    coasts = Coast.(first.(rvcoasts), last.(rvcoasts), dts)
    
    #build sequence - add initial condition to sequence? could simplify logic
    sequence = Vector{Union{Impulse, Coast}}(undef, nimp + ncoasts)


    #add coasts
    for i = 1:ncoasts
        #implicitly takes care of final_coast
        #if initial coast, coasts are on odd indices
        if init_coast
            sequence[2*i-1] = coasts[i]
        else    #coasts on even indices
            sequence[2*i] = coasts[i]
        end
    end

    #add impulses
    for i = 1:nimp
        imp_ind = (init_coast) ? 2*i : 2*i-1
        sequence[imp_ind] = impulses[i]

        #add continuity constraints between coasts
        if 1 < imp_ind < nimp+ncoasts && sequence[imp_ind-1] isa Coast && sequence[imp_ind+1] isa Coast
            last_coast = sequence[imp_ind-1]
            imp = sequence[imp_ind]
            next_coast = sequence[imp_ind+1]
            @constraint(model, last_coast.rcoast[:, end] .== next_coast.rcoast[:, 1])
            @constraint(model, next_coast.vcoast[:, 1] .== last_coast.vcoast[:, 1] + imp.deltaVmag * imp.deltaVdir)
        end
    end

    #initial boundary condition
    if init_coast
        @constraint(model, coasts[1].rcoast[:, 1] .== X1[1:3])
        @constraint(model, coasts[1].vcoast[:, 1] .== X1[4:6])
    else
        #init cond -> imp -> coast
        i = sequence[1]
        c = sequence[2]
        @constraint(model, c.rcoast[:, 1] .== X1[1:3])
        @constraint(model, c.vcoast[:, 1] .== X1[4:6] + i.deltaVmag * i.deltaVdir)
    end

    #final boundary condition
    if final_coast
        @constraint(model, coasts[end].rcoast[:, end] .== X2[1:3])
        @constraint(model, coasts[end].vcoast[:, end] .== X2[4:6])
    else
        #init cond -> imp -> coast
        i = sequence[end]
        c = sequence[end-1]
        @constraint(model, c.rcoast[:, end] .== X2[1:3])
        @constraint(model, X2[4:6] .== c.vcoast[:, end] + i.deltaVmag * i.deltaVdir)
    end

    @objective(model, Min, sum(deltaVmags))

    transfer = Transfer(X1, X2, mu, tf, sequence)

    model, transfer
end

function solved(t::Transfer)
    Transfer(
        t.X1,
        t.X2,
        t.mu,
        t.transfer_time,
        solved.(t.sequence)
    )
end

function unscale(t::Transfer, L, T)
    Transfer(
        diagm([L, L, L, L/T, L/T, L/T]) * t.X1,
        diagm([L, L, L, L/T, L/T, L/T]) * t.X2,
        L^3 / T^2 * t.mu,
        T * t.transfer_time,
        unscale.(t.sequence, L, T)
    )
end

#scale options?
function lambert_transfer_model(X1, X2, deltat, MU, N)
    r1 = X1[1:3]
    v1 = X1[4:6]
    r2 = X2[1:3]
    v2 = X2[4:6]

    model = Model(Ipopt.Optimizer)

    r, v = add_coast_segment(model, deltat, N, "", dyn=(X -> two_body_dyn(X, MU)))

    dV = @variable(model, [1:2], lower_bound=0)
    dVdirs = @variable(model, [1:3, 1:2])
    @constraint(model, [i=1:2], dVdirs[:, i]' * dVdirs[:, i] == 1)

    @constraint(model, r[:, 1] .== r1)
    @constraint(model, r[:, end] .== r2)

    @constraint(model, dV[1] * dVdirs[:, 1] == v[:, 1] - v1)
    @constraint(model, dV[2] * dVdirs[:, 2] == (v[:, end] - v2))

    @objective(model, MIN_SENSE, sum(dV))

    model, Transfer(X1, X2, MU, deltat, [
        Impulse(dV[1], dVdirs[:, 1]),
        Coast(r, v, deltat),
        Impulse(dV[2], dVdirs[:, 2])])
end
##
case_ind = 3
orb1, orb2 = ORBIT_STARTS[case_ind], ORBIT_ENDS[case_ind]
r1, v1 = kepler_to_rv(orb1)
r2, v2 = kepler_to_rv(orb2)

tf1 = orbital_period((orb1.a+orb2.a)/2, GM_EARTH)

L = (orb1.a+orb2.a)/2
T = tf1
MUPRIME = GM_EARTH * T ^ 2 / L ^ 3


X1 = [r1 / L; v1 * T / L]
X2 = [r2 / L; v2 * T / L]


model = Model(Ipopt.Optimizer)
model, model_transfer = n_impulse_model(model, X1, X2, tf1 / T, MUPRIME, 20, 2, false, false)

#this is not scaled!!!
initial_guess_initorb!(orb1, tf1, model_transfer.sequence[2].rcoast, model_transfer.sequence[2].vcoast)
##
optimize!(model)
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

    if transfer.sequence[1] isa Coast || transfer.sequence[end] isa Coast
        @warn "Unimplemented edge coast case!!!!"
    end

    tspan_ppdot = []

    for tic in two_impulse_coasts
        i1, c, i2 = tic
        dv1 = i1.deltaVmag * i1.deltaVdir
        dv2 = i2.deltaVmag * i2.deltaVdir

        orbit = rv_to_kepler(c.rcoast[:, 1], c.vcoast[:, 1])

        propagator = Propagators.init(Val(:TwoBody), orbit)

        push!(tspan_ppdot, ppdot_deltavs(propagator, dv1, dv2, c.dt, npoints; tpbvp_kwargs...))
    end

    # tspan = first.(tspan_ppdot)
    # ppdot = last.(tspan_ppdot)

    #plug diagnostic function, which is incomplete
    tspan_ppdot
end

# deltav1 = solved_model.sequence[2].vcoast[:, 1] - v1
# deltav2 = v2 - solved_model.sequence[2].vcoast[:, end]
# ##
# tspan, ppdot = ppdot_deltavs(solved_prop, deltav1, deltav2, tf1, 100)
tspan_ppdot = primer_vector(solved_model, 100)
tspan, ppdot = tspan_ppdot[1]
normp = norm.(getindex.(ppdot, Ref(1:3)))
normpdot = [dot(ppdoti[1:3], ppdoti[4:6]) / norm(ppdoti[1:3]) for ppdoti in ppdot]
# ##
pv_diag = diagnose_ppdot(normp, normpdot) #remove?
##
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
model = Model(Ipopt.Optimizer)
model, model_transfer = n_impulse_model(model, X1, X2, tf1 / T, MUPRIME, 20, 2, true, true)

initial_guess_initorb!(orb1, 0, model_transfer.sequence[1].rcoast, model_transfer.sequence[1].vcoast)
initial_guess_initorb!(orb1, tf1, model_transfer.sequence[3].rcoast, model_transfer.sequence[3].vcoast)
initial_guess_initorb!(orb1, 0, model_transfer.sequence[5].rcoast, model_transfer.sequence[5].vcoast)

optimize!(model)
##
solved_model = unscale(solved(model_transfer), L, T)
solved_orbs = [rv_to_kepler(c.rcoast[:, 1], c.vcoast[:, 1]) for c in solved_model.sequence if c isa Coast]
# solved_prop = Propagators.init(Val(:TwoBody), solved_orb)
##
f, ax3d = plot_orbit(orb1, orb2)
# add_discretized_trajectory!(ax3d, solved_model.sequence[1].rcoast)
add_discretized_trajectory!(ax3d, solved_model.sequence[3].rcoast)
# add_discretized_trajectory!(ax3d, solved_model.sequence[5].rcoast)
f