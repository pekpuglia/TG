using SatelliteToolboxBase, CasADi

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

sat_toolbox_model(::TwoBodyModel) = Val(:TwoBody)

function dynamics(X, model::TwoBodyModel)
    r = X[1:3]
    v = X[4:6]
    [
        v
        -model.mu/(√(r'*r)^3)*r
    ]
end

struct J2model <: AbstractConservativeModel
    mu
    J2_mu_R2
end

J2model(mu, J2, R) = J2model(mu, J2*mu*R^2)

scale(tbm::J2model, L, T) = J2model(tbm.mu * T ^ 2 / L ^ 3, tbm.J2_mu_R2 * T ^ 2 / L ^ 5)
unscale(tbm::J2model, L, T) = J2model(tbm.mu * L^3 / T ^ 2, tbm.J2_mu_R2 * L ^ 5 / T ^ 2)

sat_toolbox_model(::J2model) = Val(:J2osc)

function dynamics(X, model::J2model)
    r = X[1:3]
    v = X[4:6]

    rnorm = √(r'*r)

    z_r2 = r[3]^2 / (r' * r)

    j2_term = 3/2 * model.J2_mu_R2 / (r'*r)^2 * [
        r[1] / rnorm * (5*z_r2 - 1)
        r[2] / rnorm * (5*z_r2 - 1)
        r[3] / rnorm * (5*z_r2 - 3)
    ]

    [
        v
        -model.mu/(rnorm^3)*r + j2_term
    ]
end

#m
const R_TABLE_ATM = EARTH_EQUATORIAL_RADIUS .+ 1e3*[0 25 30 40 50 60 70 80 90 100 110 120 130 140 150 180 200 250 300 350 400 450 500 600 700 800 900 1000]
#kg/m^3
const RHO_TABLE_ATM = [1.225 4.008e-2 1.841e-2 3.996e-3 1.027e-3 3.097e-4 8.283e-5 1.846e-5 3.416e-6 5.606e-7 9.708e-8 2.222e-8 8.152e-9 3.831e-9 2.076e-9 5.194e-10 2.541e-10 6.073e-11 1.916e-11 7.014e-12 2.803e-12 1.184e-12 5.215e-13 1.137e-13 3.070e-14 1.136e-14 5.759e-15 3.561e-15]
#m
const H_TABLE_ATM = 1e3*[ 7.310 6.427 6.546 7.360 8.342 7.583 6.661 5.927 5.533 5.703 6.782 9.973 13.243 16.322 21.652 27.974 34.934 43.342 49.755 54.513 58.019 60.980 65.654 76.377 100.587 147.203 208.020]

struct J2DragModel <: AbstractOrbitalMechanicsModel
    mu
    J2_mu_R2
    B
    omegaE
    r_table
    rho_table
    H_table
end

# J2Dragmodel(mu, J2, R) = J2DragModel(mu, J2*mu*R^2)

scale(j2d::J2DragModel, L, T) = J2DragModel(
    j2d.mu * T ^ 2 / L ^ 3, 
    j2d.J2_mu_R2 * T ^ 2 / L ^ 5,
    j2d.B / L^2,
    j2d.omegaE * T,
    j2d.r_table / L,
    j2d.rho_table * L^3,
    j2d.H_table / L
    )

unscale(j2d::J2DragModel, L, T) = J2DragModel(
    j2d.mu * L ^ 3 / T ^ 2, 
    j2d.J2_mu_R2 * L ^ 5 / T ^ 2,
    j2d.B * L^2,
    j2d.omegaE / T,
    j2d.r_table * L,
    j2d.rho_table / L^3,
    j2d.H_table * L
)


# sat_toolbox_model(::J2DragModel) = Val(:J2osc)

function rho_model_curtis(r, r_table, rho_table, H_table)
    
    #km
    height = (r - EARTH_EQUATORIAL_RADIUS)/1e3
    if height > 1000
        height = 1000;
    elseif height < 0
        height = 0;
    end
    #Determine the interpolation interval:
    i = 1;
    for j = 1:27
        if r >= r_table[j] && height < r_table[j+1]
            i = j;
        end
    end
    if height == 1000.0
        i = 27;
    end
    
    #Exponential interpolation:
    rho_table[i]*exp(-(r - r_table[i])/H_table[i]);
end #atmopshere

Base.exp(x::SX) = casadi.exp(x)
Base.tanh(x::SX) = casadi.tanh(x)

heaviside(x, k) = (tanh(x*k) + 1) / 2
rect(x, k) = 1/2 * (tanh(k*x) + tanh(k*(1-x)))


function rho_model_smooth(r, r_table, rho_table, H_table, k=1000)
     #Exponential interpolation:
    sum(
        rho_table[i]*exp(-(r - r_table[i])/H_table[i]) * 
        rect((r - r_table[i]) / (r_table[i+1] - r_table[i]), k) 
        for i = 1:27
    ) + rho_table[1] * heaviside((-r + r_table[1])/(r_table[2] - r_table[1]), k) + rho_table[end] * heaviside((r - r_table[end])/(r_table[end] - r_table[end-1]), k)
end #atmopshere

function drag_acc(X, model::J2DragModel, rho_model = rho_model_smooth)
    r = X[1:3]
    v = X[4:6]

    rnorm = √(r'*r)

    vrel = v - cross([0; 0; model.omegaE], r)
    rho = rho_model(rnorm, model.r_table, model.rho_table, model.H_table)

    vrelnorm = √(vrel' * vrel)

    drag = -(1/2 * rho * model.B) * vrelnorm * vrel

    drag
end

function dynamics(X, model::J2DragModel, rho_model = rho_model_smooth)
    r = X[1:3]
    v = X[4:6]

    rnorm = √(r'*r)

    z_r2 = r[3]^2 / (r' * r)

    j2_term = 3/2 * model.J2_mu_R2 / (r'*r)^2 * [
        r[1] / rnorm * (5*z_r2 - 1)
        r[2] / rnorm * (5*z_r2 - 1)
        r[3] / rnorm * (5*z_r2 - 3)
    ]

    drag = drag_acc(X, model, rho_model)

    [
        v
        -model.mu/(rnorm^3)*r + j2_term + drag
    ]
end