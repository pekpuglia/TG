using LinearAlgebra
using SatelliteToolboxBase
using SatelliteToolboxPropagators
function glandorf_precompute(R, V, t)
    r = norm(R)

    H = cross(R, V)
    h = norm(H)
    
    orb = rv_to_kepler(R, V)

    p = orb.a * (1 - orb.e^2)
    θ = orb.f
    e = orb.e

    g_glandorf = ((p/h)*(p+r)*r*sin(θ) - 3*e*p*t)/(1-e^2)
    f_glandorf = - (p/h) * (p+r)*r*cos(θ)

    F3 = - (h/p)*sin(θ)
    F4 = (h/p)*(e + cos(θ))
    F5 = GM_EARTH / r^3

    F = [
        r*cos(θ)
        r*sin(θ)
        F3
        F4
        F5
        3t
        3GM_EARTH*t/r^3
        F3 + g_glandorf*F5
        F4 + f_glandorf*F5
    ]

    σ = cross(R, H)
    w = cross(V, H)
    b1 =    h + F[5]*(F[1]*f_glandorf - F[2]*g_glandorf)
    b2 = F[6]*h + F[3]*f_glandorf - F[4]*g_glandorf
    b3 = F[6]*h + 2*(F[3]*f_glandorf - F[4]*g_glandorf)
    b4 = F[1]*f_glandorf - F[2]*g_glandorf

    b = [b1; b2;b3;b4]

    F, g_glandorf, f_glandorf, σ, w, b
end

export P_glandorf
function P_glandorf(R, V, t)
    F, g_glandorf, f_glandorf, σ, w, b = glandorf_precompute(R, V, t)
    H = cross(R, V)
    
    [
        F[1]*R-g_glandorf*V F[2]*R-f_glandorf*V 2*R-F[6]*V  V    F[1]*H F[2]*H
        F[8]*R-F[1]*V         F[9]*R-F[2]*V         F[7]*R-V   -F[5]*R F[3]*H F[4]*H
    ]
end

export Pinv_glandorf
function Pinv_glandorf(R, V, t)
    F, g_glandorf, f_glandorf, σ, w, b = glandorf_precompute(R, V, t)
    H = cross(R, V)
    h = norm(H)
    
    1/h^3 * [
        F[2]*F[5]*σ' - F[4]*w'  2*F[4]*σ' - F[2]*w'
        -F[1]*F[5]*σ' + F[3]*w' -2*F[3]*σ' + F[1]*w'
        h*w'              -h*σ'
        -b[1]*σ'+b[2]*w'      -b[3]*σ' + b[4]*w'
        F[4]*H'             -F[2]*H'
        -F[3]*H'             F[1]*H'
    ]
end

function Phi_time(propagator, t)
    r, v = Propagators.propagate!(propagator, t)
    r0, v0 = kepler_to_rv(propagator.tbd.orb₀)
    Phi = P_glandorf(r, v, t) * Pinv_glandorf(r0, v0, 0)
    Phi
end

function p0dot_tpbvp(p0, pf, delta_t, prop; det_tol=1e-6, planar=false)
    Phi = Phi_time(prop, delta_t)

    M = Phi[1:3, 1:3]
    N = Phi[1:3, 4:6]

    if abs(det(N)) <= det_tol || planar
        @warn "N singular"
        #CHECK THIS IS TRUE
        #only works for coplanar transfers?
        r1, v1 = kepler_to_rv(prop.tbd.orb₀)
        A = N * [r1 v1]
        b = pf - M * p0
        p0dot = [r1 v1] * (A \ b)
    else
        p0dot = N \ (pf - M * p0)
    end

    p0dot
end


function ppdot_deltavs(transfer_propagator, deltav1, deltav2, delta_t, N; tpbvp_kwargs...)
    p0 = deltav1 / norm(deltav1)
    pf = deltav2 / norm(deltav2)
    p0dot = p0dot_tpbvp(p0, pf, delta_t, transfer_propagator; tpbvp_kwargs...)
    tspan = range(0, delta_t, N)
    ppdot = [Phi_time(transfer_propagator, t) * [p0; p0dot] for t in tspan]
    tspan, ppdot
end

#doesn't account for continuity yet
export PRIMER_DIAGNOSTIC, IC_FC, IC_LA, ED_FC, ED_LA, MID, OPT
@enum PRIMER_DIAGNOSTIC IC_FC IC_LA ED_FC ED_LA MID OPT
function diagnose_ppdot(normp, normp_dot, rtol = 1e-4)
    normp_dot0 = normp_dot[1]
    normp_dotf = normp_dot[end]
    
    max_normp_dot = maximum(abs.(normp_dot))
    tol = rtol*max_normp_dot

    if normp_dot0 > tol && normp_dotf < -tol
        IC_FC
    elseif normp_dot0 > tol && normp_dotf > tol
        IC_LA
    elseif normp_dot0 < -tol && normp_dotf < -tol
        ED_FC
    elseif normp_dot0 < -tol && normp_dotf > tol
        ED_LA
    elseif any(normp .> 1 + rtol)
        MID
    else
        OPT
    end
end

#change to numerical alg??
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
        final_r = transfer.sequence[1].rcoast[:, end]
        final_v = transfer.sequence[1].vcoast[:, end]
        #wroooooong!
        first_coast_propagator = Propagators.init(Val(:TwoBody), rv_to_kepler(final_r, final_v))
        first_coast_duration = transfer.sequence[1].dt

        #backwards propagation
        tspan = range(0, first_coast_duration, npoints)
        ppdot = []
        for t in tspan
            Phi = Phi_time(first_coast_propagator, t-first_coast_duration)
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

        #forward propagation
        tspan = range(0, last_coast_duration, npoints)
        ppdot = []
        for t in tspan
            Phi = Phi_time(last_coast_propagator, t)
            push!(ppdot, Phi*ppdot_start)
        end
        push!(tspan_ppdot, (tspan, ppdot))
    end

    # #merge everything
    tspans = first.(tspan_ppdot)
    prevtime = cumsum([0; last.(tspans[1:end-1])])
    tspan = vcat((ts .+ pv for (ts, pv) in zip(tspans, prevtime))...)
    ppdots = hcat(vcat(last.(tspan_ppdot)...)...)
    (tspan, ppdots)
end