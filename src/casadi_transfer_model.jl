using CasADi
include("casadi_interface.jl")
# using SatelliteToolboxBase

varlist(i::Impulse) = [i.deltaVmag, i.deltaVdir...]
varlist(c::Coast) = [c.rcoast..., c.vcoast..., c.dt]
varlist(t::Transfer) = vcat(varlist.(t.sequence)...)

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


#scaling agnostic!
function n_impulse_transfer(X1, X2, tf, mu, Ndisc, nimp::Int, init_coast::Bool, final_coast::Bool)
    type_sequence, ncoasts = create_sequence(nimp, init_coast, final_coast)

    deltaVmags = cvar("dVmag", nimp)
    deltaVdirs = cvar("dVdir", 3, nimp)
    
    dts = cvar("dt", ncoasts)
    rcoasts = [cvar("r$(c)_", 3, Ndisc) for c = 1:ncoasts]
    vcoasts = [cvar("v$(c)_", 3, Ndisc) for c = 1:ncoasts]

    #build Transfer
    imp_ind = 1
    cst_ind = 1
    sequence = []
    coasts = []
    for t = type_sequence
        el = if t == Impulse
            i = Impulse(deltaVmags[imp_ind], deltaVdirs[:, imp_ind])
            imp_ind += 1
            i
        elseif t == Coast
            c = Coast(rcoasts[cst_ind], vcoasts[cst_ind], dts[cst_ind])
            cst_ind += 1
            push!(coasts, c)
            c
        end
        push!(sequence, el)
    end

    transfer = Transfer(X1, X2, mu, tf, sequence)
    vars = varlist(transfer)

    planner = CasADiPlanner(vars)
    
    #bounds
    for i = 1:nimp
        add_bounds!(planner, deltaVmags[i], 0, Inf)
    end

    for i = 1:ncoasts
        add_bounds!(planner, dts[i], 0, Inf)
    end

    #total time
    add_equality!(planner, sum(dts), tf)

    #direction with unit magnitude
    for i = 1:nimp
        add_equality!(planner, deltaVdirs[:, i]' * deltaVdirs[:, i], 1)
    end

    f = X -> two_body_dyn(X, mu)

    #coast integration
    for c = 1:ncoasts
        dtc = dts[c] / (N-1)
        for i = 1:(N-1)
            Xi = [rcoasts[c][:, i]; vcoasts[c][:, i]]
            Xi1 = [rcoasts[c][:, i+1]; vcoasts[c][:, i+1]]
            
            add_equality!(planner, Xi1 - RK8(f, Xi, dtc), zeros(6))
        end
    end

    #initial boundary condition
    if init_coast
        c = sequence[1]
        add_equality!(planner, c.rcoast[:, 1], X1[1:3])
        add_equality!(planner, c.vcoast[:, 1], X1[4:6]) 
    else
        i = sequence[1]
        c = sequence[2]

        add_equality!(planner, c.rcoast[:, 1], X1[1:3])
        add_equality!(planner, c.vcoast[:, 1] - i.deltaVmag * i.deltaVdir, X1[4:6])
    end
    #intermediate boundary conditions (coast, impulse, coast)
    for (c1, i, c2) in zip(sequence[1:end-2], sequence[2:end-1], sequence[3:end])
        if c1 isa Coast && i isa Impulse && c2 isa Coast
            add_equality!(planner, c1.rcoast[:, end] - c2.rcoast[:, 1], zeros(3))
            add_equality!(planner, c1.vcoast[:, end] - (c2.rcoast[:, 1] + i.deltaVmag * i.deltaVdir), zeros(3))
        end
    end

    #final boundary conditions
    if final_coast
        c = sequence[end]
        add_equality!(planner, c.rcoast[:, end], X2[1:3])
        add_equality!(planner, c.vcoast[:, end], X2[4:6]) 
    else
        i = sequence[end]
        c = sequence[end-1]

        add_equality!(planner, c.rcoast[:, end], X2[1:3])
        add_equality!(planner, c.vcoast[:, end] + i.deltaVmag * i.deltaVdir, X2[4:6])
    end

    planner.prob["f"] = sum(deltaVmags)

    planner, transfer
end

#redo with casadi
# solved(i::Impulse) = Impulse(value(i.deltaVmag), value.(i.deltaVdir))
# solved(c::Coast) = Coast(value.(c.rcoast), value.(c.vcoast), value(c.dt))
# function solved(t::Transfer)
#     Transfer(
#         t.X1,
#         t.X2,
#         t.mu,
#         t.transfer_time,
#         solved.(t.sequence)
#     )
# end