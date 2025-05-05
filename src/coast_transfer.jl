export CoastTransfer
struct CoastTransfer
    orb1::KeplerianElements
    orb2::KeplerianElements
    lambert::LambertResult
    p0_p0dot
end

function CoastTransfer(orb1::KeplerianElements, orb2::KeplerianElements; prograde=true)
    r2, v2             = kepler_to_rv(orb2)
    r1, v1             = kepler_to_rv(orb1)
    lamb = lambert(r1, r2, (orb2.t - orb1.t)*86400, RAAN=orb1.Ω, i=orb1.i, prograde=prograde)

    if !lamb.is_elliptic
        return CoastTransfer(orb1, orb2, lamb, nothing)
    end

    propagator = Propagators.init(Val(:TwoBody), lamb.propagator.tbd.orb₀)

    deltaV0 = lamb.v1 - v1
    deltaVf = v2 - lamb.v2

    p0 = deltaV0 / norm(deltaV0)
    pf = deltaVf / norm(deltaVf)

    Phi = P_glandorf(r2, lamb.v2, transfer_time) * Pinv_glandorf(r1, lamb.v1, 0)
    M = Phi[1:3, 1:3]
    N = Phi[1:3, 4:6]

    if abs(det(N)) <= 1e-10
        #CHECK THIS IS TRUE
        #only works for coplanar transfers?
        A = N * [r1 lamb.v1]
        b = pf - M * p0
        p0dot = [r1 lamb.v1] * (A \ b)
    else
        p0dot = N \ (pf - M * p0)
    end

    CoastTransfer(orb1, orb2, lamb, [p0; p0dot])
end

#min of prograde, retrograde
function CoastTransfer(orb1::KeplerianElements, orb2::KeplerianElements)
    prograde   = CoastTransfer(orb1, orb2, prograde=true)
    retrograde = CoastTransfer(orb1, orb2, prograde=false)

    if prograde.lambert.is_elliptic && retrograde.lambert.is_elliptic
        if cost(prograde) <= cost(retrograde)
            prograde
        else
            retrograde
        end
    elseif prograde.lambert.is_elliptic && !retrograde.lambert.is_elliptic
        prograde
    elseif !prograde.lambert.is_elliptic && retrograde.lambert.is_elliptic
        retrograde
    else
        prograde
    end
end

export transfer_start, transfer_end, cost, p, pdot

transfer_start(lt::CoastTransfer) = lt.lambert.propagator.tbd.orb₀
transfer_end(lt::CoastTransfer) = lt.lambert.propagator.tbd.orbk

function cost(lt::CoastTransfer)
    if lt.lambert.is_elliptic
        norm(lt.lambert.v1 - kepler_to_rv(lt.orb1)[2]) + norm(lt.lambert.v2 - kepler_to_rv(lt.orb2)[2])
    else
        NaN
    end
end

function cost(::Nothing)
    NaN
end


p(lt::CoastTransfer, t) = Phi_time(lt.lambert.propagator, t)[1:3, :] * lt.p0_p0dot
pdot(lt::CoastTransfer, t) = Phi_time(lt.lambert.propagator, t)[4:6, :] * lt.p0_p0dot