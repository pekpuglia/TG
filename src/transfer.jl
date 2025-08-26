using LinearAlgebra
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
    mu::Float64
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
        L^3 / T^2 * t.mu,
        T * t.transfer_time,
        unscale.(t.sequence, L, T)
    )
end

function scale(t::Transfer, L, T)
    Transfer(
        diagm([L, L, L, L/T, L/T, L/T]) \ t.X1,
        diagm([L, L, L, L/T, L/T, L/T]) \ t.X2,
        T^2 / L^3 * t.mu,
        T \ t.transfer_time,
        scale.(t.sequence, L, T)
    )
end