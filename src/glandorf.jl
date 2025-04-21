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