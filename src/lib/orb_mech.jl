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

struct J2model <: AbstractOrbitalMechanicsModel
    mu
    J2_mu_R2
end

J2model(mu, J2, R) = J2model(mu, J2*mu*R^2)

scale(tbm::J2model, L, T) = J2model(tbm.mu * T ^ 2 / L ^ 3, tbm.J2_mu_R2 * T ^ 2 / L ^ 5)
unscale(tbm::J2model, L, T) = J2model(tbm.mu * L^3 / T ^ 2, tbm.J2_mu_R2 * L ^ 5 / T ^ 2)

function dynamics(X, model::J2model)
    r = X[1:3]
    v = X[4:6]

    rnorm = √(r'*r)

    z_r2 = r[3]^2 / (r' * r)

    j2_term = 3*model.J2_mu_R2 / (2*(r'*r)^2) * [
        r[1] / rnorm * (5*z_r2 - 1)
        r[2] / rnorm * (5*z_r2 - 1)
        r[3] / rnorm * (5*z_r2 - 3)
    ]

    [
        v
        -model.mu/(rnorm^3)*r + j2_term
    ]
end