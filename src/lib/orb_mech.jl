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

struct J2DragModel <: AbstractOrbitalMechanicsModel
    mu
    J2_mu_R2
    omegaE
    B
    r0
    rho0
    H0
end

J2Dragmodel(mu, J2, R) = J2DragModel(mu, J2*mu*R^2)

scale(tbm::J2DragModel, L, T) = J2DragModel(tbm.mu * T ^ 2 / L ^ 3, tbm.J2_mu_R2 * T ^ 2 / L ^ 5)
unscale(tbm::J2DragModel, L, T) = J2DragModel(tbm.mu * L^3 / T ^ 2, tbm.J2_mu_R2 * L ^ 5 / T ^ 2)

# sat_toolbox_model(::J2DragModel) = Val(:J2osc)

function rho_model_if(height)
    #
    # ATMOSPHERE calculates density for altitudes from sea level
    # through 1000 km using exponential interpolation.
    # 
    #Geometric altitudes (km):
    h = 
    [0 25 30 40 50 60 70 80 90 100 110 120 130 140 150 180 200 250 300 350 400 450 500 600 700 800 900 1000];
    
    #Corresponding densities (kg/m^3) from USSA76:
    r = 
    [1.225 4.008e-2 1.841e-2 3.996e-3 1.027e-3 3.097e-4 8.283e-5 1.846e-5 3.416e-6 5.606e-7 9.708e-8 2.222e-8 8.152e-9 3.831e-9 2.076e-9 5.194e-10 2.541e-10 6.073e-11 1.916e-11 7.014e-12 2.803e-12 1.184e-12 5.215e-13 1.137e-13 3.070e-14 1.136e-14 5.759e-15 3.561e-15];
    #Scale heights (km):
    H = 
    [ 7.310 6.427 6.546 7.360 8.342 7.583 6.661 5.927 5.533 5.703 6.782 9.973 13.243 16.322 21.652 27.974 34.934 43.342 49.755 54.513 58.019 60.980 65.654 76.377 100.587 147.203 208.020];
    #Handle altitudes outside of the range:
    if height > 1000
        height = 1000;
    elseif height < 0
        height = 0;
    end
    #Determine the interpolation interval:
    i = 1;
    for j = 1:27
        if height >= h[j] && height < h[j+1]
            i = j;
        end
    end
    if height == 1000.0
        i = 27;
    end
    
    #Exponential interpolation:
    r[i]*exp(-(height - h[i])/H[i]);
end #atmopshere


heaviside(x, k) = (tanh(x*k) + 1) / 2
rect(x, k) = 1/2 * (tanh(k*x) + tanh(k*(1-x)))

function rho_model_smooth(height, k=50)
    #
    # ATMOSPHERE calculates density for altitudes from sea level
    # through 1000 km using exponential interpolation.
    # 
    #Geometric altitudes (km):
    h = 
    [0 25 30 40 50 60 70 80 90 100 110 120 130 140 150 180 200 250 300 350 400 450 500 600 700 800 900 1000];
    
    #Corresponding densities (kg/m^3) from USSA76:
    r = 
    [1.225 4.008e-2 1.841e-2 3.996e-3 1.027e-3 3.097e-4 8.283e-5 1.846e-5 3.416e-6 5.606e-7 9.708e-8 2.222e-8 8.152e-9 3.831e-9 2.076e-9 5.194e-10 2.541e-10 6.073e-11 1.916e-11 7.014e-12 2.803e-12 1.184e-12 5.215e-13 1.137e-13 3.070e-14 1.136e-14 5.759e-15 3.561e-15];
    #Scale heights (km):
    H = 
    [ 7.310 6.427 6.546 7.360 8.342 7.583 6.661 5.927 5.533 5.703 6.782 9.973 13.243 16.322 21.652 27.974 34.934 43.342 49.755 54.513 58.019 60.980 65.654 76.377 100.587 147.203 208.020];
    #Handle altitudes outside of the range:
    # if height > 1000
    #     height = 1000;
    # elseif height < 0
    #     height = 0;
    # end
    # #Determine the interpolation interval:
    # i = 1;
    # for j = 1:27
    #     if height >= h[j] && height < h[j+1]
    #         i = j;
    #     end
    # end
    # if height == 1000.0
    #     i = 27;
    # end
    
    #Exponential interpolation:
    sum(
        r[i]*exp(-(height - h[i])/H[i]) * 
        rect((height - h[i]) / (h[i+1] - h[i]), k) 
        for i = 1:27
    ) + r[1] * heaviside(-height, k) + r[end] * heaviside(height - h[end], k)
end #atmopshere


function dynamics(X, model::J2DragModel)
    r = X[1:3]
    v = X[4:6]

    rnorm = √(r'*r)

    z_r2 = r[3]^2 / (r' * r)

    j2_term = 3/2 * model.J2_mu_R2 / (r'*r)^2 * [
        r[1] / rnorm * (5*z_r2 - 1)
        r[2] / rnorm * (5*z_r2 - 1)
        r[3] / rnorm * (5*z_r2 - 3)
    ]

    vrel = v - cross([0; 0; model.omegaE], r)
    rho = rhomodel(model, rnorm)

    vrelnorm = √(vrel' * vrel)

    drag = -1/2 * rho * vrelnorm * model.B * vrel

    [
        v
        -model.mu/(rnorm^3)*r + drag
    ]
end