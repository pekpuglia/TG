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
## REPRODUCE ORBITS IN INTERACTIVE PRIMER VECTOR
function initial_guess_initorb!(r, v, orb1, tf, L=1, T=1, ret_final_orb=false)
    prop = Propagators.init(Val(:TwoBody), orb1)

    N = size(r)[2]

    for i = 1:N
        rguess, vguess = Propagators.propagate!(prop, (i-1)*tf/(N-1))
        set_start_value.(r[:, i], rguess / L)
        set_start_value.(v[:, i], vguess * T/L)
    end

    if ret_final_orb
        prop.tbd.orbk
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

function create_sequence(nimp::Int, init_coast::Bool, final_coast::Bool)
    ncoasts = nimp - 1 + init_coast + final_coast
    
    sequence = Vector{Type}(undef, nimp+ncoasts)

    if init_coast
        sequence[1:2:end] .= Coast
        sequence[2:2:end] .= Impulse
    else
        sequence[1:2:end] .= Impulse
        sequence[2:2:end] .= Coast
    end

    sequence, ncoasts
end

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
            tabr = zeros(3, Ndisc+1)
            tabv = zeros(3, Ndisc+1)
            for i = 1:Ndisc+1
                r, v = Propagators.propagate!(prop, coast_times[coast_ind] * (i-1) / Ndisc)
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

#scaling agnostic!
function n_impulse_model(model, X1, X2, tf, mu, Ndisc, nimp::Int, init_coast::Bool, final_coast::Bool)
    ncoasts = nimp - 1 + init_coast + final_coast

    dts = @variable(model, [1:ncoasts], base_name = "dt", lower_bound=0)#, upper_bound=tf)
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

##
case_ind = 3
orb1, orb2 = ORBIT_STARTS[case_ind], ORBIT_ENDS[case_ind]
r1, v1 = kepler_to_rv(orb1)
r2, v2 = kepler_to_rv(orb2)

#vary tf for each orbit
tf1 = orbital_period((orb1.a+orb2.a)/2, GM_EARTH)

L = (orb1.a+orb2.a)/2
T = tf1
MUPRIME = GM_EARTH * T ^ 2 / L ^ 3


X1 = [r1 / L; v1 * T / L]
X2 = [r2 / L; v2 * T / L]


model = Model(optimizer_with_attributes(Ipopt.Optimizer,
    "max_iter" => 3_000,
    "max_wall_time" => 30.0
))
model, model_transfer = n_impulse_model(model, X1, X2, tf1 / T, MUPRIME, 100, 2, false, false)

initial_guess_initorb!(model_transfer.sequence[2].rcoast, model_transfer.sequence[2].vcoast, orb1, 0, L, T)
#try different initial conditions - make list and bick best
# initial_guess_initorb!(model_transfer.sequence[2].rcoast, model_transfer.sequence[2].vcoast, orb1, tf1, L, T)
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

#automate this - discard early departure/late arrival
tspan_ppdot = primer_vector(solved_model, 100)
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