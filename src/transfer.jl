using LinearAlgebra
include("orb_mech.jl")
struct Impulse
    deltaVmag
    deltaVdir
end
unscale(i::Impulse, L, T) = Impulse(L/T * i.deltaVmag, i.deltaVdir)
scale(i::Impulse, L, T) = Impulse(T/L * i.deltaVmag, i.deltaVdir)

struct Coast
    rcoast::Matrix
    vcoast::Matrix
    dt
end
unscale(c::Coast, L, T) = Coast(L * c.rcoast, L/T * c.vcoast, T * c.dt)
scale(c::Coast, L, T) = Coast(L \ c.rcoast, T/L * c.vcoast, T \ c.dt)

#add nimp, ncoasts, init_coast, final_coast to struct
struct Transfer
    X1::Vector
    X2::Vector
    model::AbstractOrbitalMechanicsModel
    transfer_time
    # rcoast::Array
    # vcoast::Array
    # deltaVmags::Vector
    # deltaVdirs::Matrix
    # dts::Vector
    sequence::Vector{Union{Impulse, Coast}}
end
total_dV(t::Transfer) = sum(el.deltaVmag for el in t.sequence if el isa Impulse)

impulses(t::Transfer) = [i for i in t.sequence if i isa Impulse]
coasts(t::Transfer) = [c for c in t.sequence if c isa Coast]


function unscale(t::Transfer, L, T)
    Transfer(
        diagm([L, L, L, L/T, L/T, L/T]) * t.X1,
        diagm([L, L, L, L/T, L/T, L/T]) * t.X2,
        unscale(t.model, L, T),
        T * t.transfer_time,
        unscale.(t.sequence, L, T)
    )
end

function scale(t::Transfer, L, T)
    Transfer(
        diagm([L, L, L, L/T, L/T, L/T]) \ t.X1,
        diagm([L, L, L, L/T, L/T, L/T]) \ t.X2,
        scale(t.model, L, T),
        T \ t.transfer_time,
        scale.(t.sequence, L, T)
    )
end

function impulse_times(t::Transfer)
    ts = []

    last_time = 0.0

    for el in t.sequence
        if el isa Coast
            last_time += el.dt
        elseif el isa Impulse
            push!(ts, last_time)
        end
    end
    ts
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

function add_null_impulse(solved_transfer::Transfer, tspan_ppdot)
    tspan, ppdot = tspan_ppdot

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

    N = size(coast_to_split.rcoast)[2]

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

    new_transfer
end

function rediscretize_coast(c::Coast, model::AbstractOrbitalMechanicsModel, N, integrator)
    x = [c.rcoast[:, 1]; c.vcoast[:, 1]]

    rcoast = zeros(3, N)
    rcoast[:, 1] = x[1:3]
    vcoast = zeros(3, N)
    vcoast[:, 1] = x[4:6]


    for i = 2:N
        x = integrator(X -> dynamics(X, model), x, c.dt / (N-1))

        rcoast[:, i] = x[1:3]
        vcoast[:, i] = x[4:6]
    end

    Coast(rcoast, vcoast, c.dt)
end

function rediscretize_transfer(t::Transfer, N, integrator=RK8)
    new_transfer = deepcopy(t)

    coast_indices = findall(el -> el isa Coast, t.sequence)

    for (i, c) = enumerate(coasts(t))
        new_transfer.sequence[coast_indices[i]] = rediscretize_coast(c, t.model, N, integrator)
    end
    new_transfer
end