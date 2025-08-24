using JuMP
using SatelliteToolboxBase


function add_coast_segment(model, deltat, N, ind; dyn=(X -> two_body_dyn(X, GM_EARTH)), integrator=RK4)
    #[(x, y, z), time instants]
    r = @variable(model, [1:3, 1:N+1], base_name="r_coast_$ind")
    v = @variable(model, [1:3, 1:N+1], base_name="v_coast_$ind")


    for i = 1:(N+1)
        # @constraint(model, r[:, i]' * r[:, i] >= EARTH_EQUATORIAL_RADIUS^2)
        # @constraint(model, cross(r[:, i], v[:, i])[3] >=0)
    end

    for i = 1:N
        @constraint(model, integrator(dyn, [r[:, i]; v[:, i]], deltat/N) .== [r[:, i+1]; v[:, i+1]])
    end

    r, v
end

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

#scaling agnostic!
function n_impulse_model(model, X1, X2, tf, mu, Ndisc, nimp::Int, init_coast::Bool, final_coast::Bool)
    ncoasts = nimp - 1 + init_coast + final_coast

    dts = @variable(model, [1:ncoasts], base_name = "dt", lower_bound=0)#, upper_bound=tf)
    @constraint(model, sum(dts) >= tf)

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

solved(i::Impulse) = Impulse(value(i.deltaVmag), value.(i.deltaVdir))
solved(c::Coast) = Coast(value.(c.rcoast), value.(c.vcoast), value(c.dt))
function solved(t::Transfer)
    Transfer(
        t.X1,
        t.X2,
        t.mu,
        t.transfer_time,
        solved.(t.sequence)
    )
end