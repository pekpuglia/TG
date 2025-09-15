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