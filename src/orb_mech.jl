using SatelliteToolboxBase

function hohmann_cost(a1, a2, mu)
    v1 = sqrt(mu/a1)
    v2 = sqrt(mu/a2)

    a_transf = (a1+a2)/2

    vt_peri = sqrt(2*mu*(1/a1 - 1/(2*a_transf)))
    vt_apo = sqrt(2*mu*(1/a2 - 1/(2*a_transf)))

    v2-vt_apo + vt_peri-v1
end

function orbital_period(orb::KeplerianElements, m0)
    orbital_period(orb.a,m0)
end

orbital_period(a, m0) = 2π√(a^3/m0)

abstract type AbstractOrbitalMechanicsModel end

normalize(m::AbstractOrbitalMechanicsModel, L, T) = error("not implemented")

abstract type AbstractConservativeModel <: AbstractOrbitalMechanicsModel end

struct TwoBodyModel <: AbstractConservativeModel
    mu
end

scale(tbm::TwoBodyModel, L, T) = TwoBodyModel(tbm.mu * T ^ 2 / L ^ 3)
unscale(tbm::TwoBodyModel, L, T) = TwoBodyModel(tbm.mu * L^3 / T ^ 2)

function dynamics(X, model::TwoBodyModel)
    r = X[1:3]
    v = X[4:6]
    [
        v
        -model.mu/(√(r'*r)^3)*r
    ]
end

# struct J2model <: AbstractOrbitalMechanicsModel
#     mu
#     J2
# end