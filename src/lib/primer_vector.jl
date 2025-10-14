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

abstract type AbstractPVTMAlgo end

struct PVTMGlandorf <: AbstractPVTMAlgo end

function Phi_time(::PVTMGlandorf, ::TwoBodyModel, x0, t)
    r0, v0 = x0[1:3], x0[4:6]
    propagator = Propagators.init(Val(:TwoBody), rv_to_kepler(r0, v0))
    r, v = Propagators.propagate!(propagator, t)
    Phi = P_glandorf(r, v, t) * Pinv_glandorf(r0, v0, 0)
    Phi
end

using ForwardDiff


struct PVTMFromSTM <: AbstractPVTMAlgo
    N
    integrator
end

function Phi_time(pvtm_algo::PVTMFromSTM, orb_model::AbstractOrbitalMechanicsModel, x0, t)
    end_state(initial_state) = final_X(X -> dynamics(X, orb_model), initial_state, t, pvtm_algo.N, pvtm_algo.integrator)
    
    ForwardDiff.jacobian(end_state, x0)
end

# Phi with integration of diff eqs
struct PVTMFromODE <: AbstractPVTMAlgo
    N
    integrator
end

function Phi_time(pvtm_algo::PVTMFromODE, orb_model::AbstractOrbitalMechanicsModel, x0, t)
    f = X -> dynamics(X, orb_model)

    lambdadot_matrix(x) = - ForwardDiff.jacobian(f, x)'
    #MAYBE THIS ONLY WORKS FOR CONSERVATIVE MODEL???
    J = [zeros(3, 3) I(3)
        -I(3) zeros(3,3)]
    pvdot_matrix(x) = J' * lambdadot_matrix(x) * J
    
    pvdot(p, x) = pvdot_matrix(x) * p
    
    full_sys(xp) = [
        f(xp[1:6])
        pvdot(xp[7:12], xp[1:6])
    ]
    
    Phi_tf_t0 = I(6) |> Matrix{Float64}

    #integration per column
    for i = 1:6
        Phi_tf_t0[:, i] = final_X(full_sys, [x0; Phi_tf_t0[:, i]], t, pvtm_algo.N, pvtm_algo.integrator)[7:12]
    end

    Phi_tf_t0
end

function p0dot_tpbvp(p0, pf, delta_t, x0, pvtm_algo::AbstractPVTMAlgo, orb_model::AbstractOrbitalMechanicsModel)
    Phi = Phi_time(pvtm_algo, orb_model, x0, delta_t)

    pdot0 = SX("pdot0", 3)

    prob = Dict(
        "f" => pdot0' * pdot0,
        "x" => pdot0,
        "g" => vcat(pf - Phi[1:3, :] * [p0; sx_iterator(pdot0)...])
    )

    solver = casadi.nlpsol("S", "ipopt", prob, Dict("ipopt" => Dict(
        "max_iter" => 3000,
        "constr_viol_tol" => 1e-5)))

    sol = solver(x0 = [0.0,0,0], lbg=[0.0,0,0], ubg=[0.0,0,0])

    sol["x"].toarray() |> vec
end


function ppdot_deltavs(pvtm::AbstractPVTMAlgo, orb_model::AbstractOrbitalMechanicsModel, x0, deltav1, deltav2, delta_t, N; tpbvp_kwargs...)
    p0 = deltav1 / norm(deltav1)
    pf = deltav2 / norm(deltav2)
    p0dot = p0dot_tpbvp(p0, pf, delta_t, x0, pvtm, orb_model)
    tspan = range(0, delta_t, N)
    ppdot = [Phi_time(pvtm, orb_model, x0, t) * [p0; p0dot] for t in tspan]
    tspan, ppdot
end

#doesn't account for continuity yet
export PRIMER_DIAGNOSTIC, IC_FC, IC_LA, ED_FC, ED_LA, MID, OPT
@enum PRIMER_DIAGNOSTIC IC_FC IC_LA ED_FC ED_LA MID OPT
function diagnose_ppdot(tspan_ppdot, rtol = 1e-4)

    tspan, ppdot = tspan_ppdot
    normp = norm.(eachcol(ppdot[1:3, :]))
    normpdot = [dot(ppdoti[1:3], ppdoti[4:6]) / norm(ppdoti[1:3]) for ppdoti in eachcol(ppdot)]

    normp_dot0 = normpdot[1]
    normp_dotf = normpdot[end]
    
    max_normp_dot = maximum(abs.(normpdot))
    tol = rtol*max_normp_dot

    if normp_dot0 > tol && normp_dotf < -tol
        "Initial + Final coast"
    elseif normp_dot0 > tol && normp_dotf > tol
        "Initial coast"
    elseif normp_dot0 < -tol && normp_dotf < -tol
        "Final coast"
    elseif any(normp .> 1 + rtol)
        "Add impulse"
    # elseif normp_dot0 < -tol && normp_dotf > tol
    #     "More time"
    else
        "Local optimum"
    end
end

#change to numerical alg??
function primer_vector(transfer::Transfer, pvtm::AbstractPVTMAlgo, npoints; tpbvp_kwargs...)
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

        x0 = [c.rcoast[:, 1]; c.vcoast[:, 1]]
        # orbit = rv_to_kepler(c.rcoast[:, 1], c.vcoast[:, 1])

        # propagator = Propagators.init(Val(:TwoBody), orbit)

        push!(tspan_ppdot, ppdot_deltavs(pvtm, transfer.model, x0, dv1, dv2, c.dt, npoints))
    end

    if transfer.sequence[1] isa Coast
        #ppdot at the first impulse
        #[first coast between 2 impulses][2nd element in (tspan, ppdot)][first ppdot in the coast]
        ppdot_end = tspan_ppdot[1][2][1]
        final_r = transfer.sequence[1].rcoast[:, end]
        final_v = transfer.sequence[1].vcoast[:, end]
        final_x = [final_r; final_v]
        first_coast_duration = transfer.sequence[1].dt

        #backwards propagation
        tspan = range(0, first_coast_duration, npoints)
        ppdot = []
        for t in tspan
            Phi = Phi_time(pvtm, transfer.model,  final_x, t-first_coast_duration)
            push!(ppdot, Phi*ppdot_end)
        end
        tspan_ppdot = [(tspan, ppdot), tspan_ppdot...]
    end

    if transfer.sequence[end] isa Coast
        #ppdot at the last impulse
        #[last coast between 2 impulses][2nd element in (tspan, ppdot)][last ppdot in the coast]
        ppdot_start = tspan_ppdot[end][2][end]

        first_x = [transfer.sequence[end].rcoast[:, 1]; transfer.sequence[end].vcoast[:, 1]]
        last_coast_duration = transfer.sequence[end].dt

        #forward propagation
        tspan = range(0, last_coast_duration, npoints)
        ppdot = []
        for t in tspan
            Phi = Phi_time(pvtm, transfer.model, first_x, t)
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