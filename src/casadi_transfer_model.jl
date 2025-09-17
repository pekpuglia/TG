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

scale(x::Vector, L, T) = [x[1:3] / L; x[4:6] * T / L]

#scaling agnostic!
function n_impulse_transfer(model::AbstractOrbitalMechanicsModel, X1, X2, tf, Ndisc, nimp::Int, init_coast::Bool, final_coast::Bool)
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

    transfer = Transfer(X1, X2, model, tf, sequence)
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

    f = X -> dynamics(X, model)

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
            add_equality!(planner, c2.vcoast[:, 1] - (c1.vcoast[:, end] + i.deltaVmag * i.deltaVdir), zeros(3))
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

function sol_to_transfer(sol::Dict, casadi_transfer::Transfer)
    x = sol["x"].toarray()
    variables = varlist(casadi_transfer)

    solved_sequence = []

    for el = casadi_transfer.sequence
        if el isa Impulse
            imp_ind = findfirst(v -> casadi.is_equal(v, el.deltaVmag), variables)
            push!(solved_sequence, Impulse(x[imp_ind], x[(imp_ind+1):(imp_ind+3)]))
        elseif el isa Coast
            cst_ind = findfirst(v -> casadi.is_equal(v, el.rcoast[1, 1]), variables)
            cst_array_len = length(el.rcoast)
            cst_array_size = size(el.rcoast)
            push!(solved_sequence, Coast(
                reshape(x[(cst_ind                ):(cst_ind +cst_array_len-1)], cst_array_size),
                reshape(x[(cst_ind+cst_array_len):(cst_ind+2cst_array_len-1)], cst_array_size),
                x[cst_ind+2cst_array_len]
            ))
        end
    end
    Transfer(
        casadi_transfer.X1,
        casadi_transfer.X2,
        casadi_transfer.model,
        casadi_transfer.transfer_time,
        solved_sequence)
end