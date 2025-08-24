using CasADi
include("casadi_interface.jl")
# using SatelliteToolboxBase

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

    # dts = @variable(model, [1:ncoasts], base_name = "dt", lower_bound=0)#, upper_bound=tf)
    dts = cvar("dt", ncoasts)
    # deltaVmags = @variable(model, [1:nimp], lower_bound = 0, base_name = "dVmag")
    # deltaVdirs = @variable(model, [1:3, 1:nimp], base_name = "dVdir")
    deltaVmags = cvar("dVmag", nimp)
    deltaVdirs = cvar("dVdir", 3, nimp)
    
    rcoasts = [cvar("r$(c)_", 3, Ndisc) for c = 1:ncoasts]
    vcoasts = [cvar("v$(c)_", 3, Ndisc) for c = 1:ncoasts]

    #build Transfer
    imp_ind = 1
    cst_ind = 1
    sequence = []
    for t = type_sequence
        el = if t == Impulse
            i = Impulse(deltaVmags[imp_ind], deltaVdirs[:, imp_ind])
            imp_ind += 1
            i
        elseif t == Coast
            c = Coast(rcoasts[cst_ind], vcoasts[cst_ind], dts[cst_ind])
            cst_ind += 1
            c
        end
        push!(sequence, el)
    end
    # rvcoasts = [add_coast_segment(
    # model, dts[i], Ndisc, i,
    # dyn=(X -> two_body_dyn(X, mu))) for i = 1:ncoasts]
    
    # coasts = Coast.(first.(rvcoasts), last.(rvcoasts), dts)
    
    # @constraint(model, sum(dts) >= tf)
    # @constraint(model, [i=1:nimp], deltaVdirs[:, i]' * deltaVdirs[:, i] == 1)

    # impulses = Impulse.(deltaVmags, eachcol(deltaVdirs))
    
    # #build sequence - add initial condition to sequence? could simplify logic
    # sequence = Vector{Union{Impulse, Coast}}(undef, nimp + ncoasts)


    # #add coasts
    # for i = 1:ncoasts
    #     #implicitly takes care of final_coast
    #     #if initial coast, coasts are on odd indices
    #     if init_coast
    #         sequence[2*i-1] = coasts[i]
    #     else    #coasts on even indices
    #         sequence[2*i] = coasts[i]
    #     end
    # end

    # #add impulses
    # for i = 1:nimp
    #     imp_ind = (init_coast) ? 2*i : 2*i-1
    #     sequence[imp_ind] = impulses[i]

    #     #add continuity constraints between coasts
    #     if 1 < imp_ind < nimp+ncoasts && sequence[imp_ind-1] isa Coast && sequence[imp_ind+1] isa Coast
    #         last_coast = sequence[imp_ind-1]
    #         imp = sequence[imp_ind]
    #         next_coast = sequence[imp_ind+1]
    #         @constraint(model, last_coast.rcoast[:, end] .== next_coast.rcoast[:, 1])
    #         @constraint(model, next_coast.vcoast[:, 1] .== last_coast.vcoast[:, 1] + imp.deltaVmag * imp.deltaVdir)
    #     end
    # end

    # #initial boundary condition
    # if init_coast
    #     @constraint(model, coasts[1].rcoast[:, 1] .== X1[1:3])
    #     @constraint(model, coasts[1].vcoast[:, 1] .== X1[4:6])
    # else
    #     #init cond -> imp -> coast
    #     i = sequence[1]
    #     c = sequence[2]
    #     @constraint(model, c.rcoast[:, 1] .== X1[1:3])
    #     @constraint(model, c.vcoast[:, 1] .== X1[4:6] + i.deltaVmag * i.deltaVdir)
    # end

    # #final boundary condition
    # if final_coast
    #     @constraint(model, coasts[end].rcoast[:, end] .== X2[1:3])
    #     @constraint(model, coasts[end].vcoast[:, end] .== X2[4:6])
    # else
    #     #init cond -> imp -> coast
    #     i = sequence[end]
    #     c = sequence[end-1]
    #     @constraint(model, c.rcoast[:, end] .== X2[1:3])
    #     @constraint(model, X2[4:6] .== c.vcoast[:, end] + i.deltaVmag * i.deltaVdir)
    # end

    # @objective(model, Min, sum(deltaVmags))

    # transfer = Transfer(X1, X2, mu, tf, sequence)

    Transfer(X1, X2, mu, tf, sequence)
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