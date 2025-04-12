using SatelliteToolboxBase
using SatelliteToolboxPropagators
using JuMP
using Ipopt
using Setfield
using ForwardDiff
include("TG.jl")
using .TG
using LinearAlgebra
## example 3.1
rp = 9600e3
ra = 21000e3
a = (rp + ra) / 2
e = (ra - rp) / (ra + rp)
orb0 = KeplerianElements(
    date_to_jd(2023, 1, 1, 0, 0, 0),
    a,
    e,
    30 |> deg2rad,
    0    |> deg2rad,
    0     |> deg2rad,
    0     |> deg2rad
)
r0, v0 = kepler_to_rv(orb0)
orbf = @set orb0.f = 120 |> deg2rad
r, v = kepler_to_rv(orbf)
plot_orbit(orb0, orbf)
##
function add_orbital_elements_fix!(model)
    Vorb_sup = √(GM_EARTH/EARTH_EQUATORIAL_RADIUS)
    r = @variable(model, [1:3], start = EARTH_EQUATORIAL_RADIUS)
    v = @variable(model, [1:3])
    set_start_value(v[1], Vorb_sup)

    #adding exc as variable so bounds will always be respected
    #then need to put constraint on it and implement E and M
    #deg!!!
    e = @variable(model, lower_bound = 0, upper_bound = 1) 
    i = @variable(model, lower_bound = 0, upper_bound = 180, base_name = "i")
    Ω = @variable(model, base_name = "Ω")
    ω = @variable(model, base_name = "ω")
    nu = @variable(model, base_name = "nu")

    #rad!!!
    M = @variable(model, base_name = "M")
    E = @variable(model, base_name= "E")
    
    hvec = cross(r, v)
    h = √(hvec' * hvec)

    #curtis chap 4
    r_perifocal = h^2/GM_EARTH * 1/(1+e*cosd(nu)) * [cosd(nu); sind(nu); 0]

    v_perifocal = GM_EARTH / h * [-sind(nu); e + cosd(nu); 0]

    R3Omega = [
        cosd(Ω) -sind(Ω) 0
        sind(Ω)  cosd(Ω) 0
        0           0         1
    ]

    R1i = [
        1 0          0
        0 cosd(i) -sind(i)
        0 sind(i)  cosd(i)
    ]

    R3omega = [
        cosd(Ω) -sind(Ω) 0
        sind(Ω)  cosd(Ω) 0
        0           0         1    
    ]

    QXxbar = R3omega * R1i * R3Omega

    @constraint(model, r .== QXxbar' * r_perifocal)
    @constraint(model, v .== QXxbar' * v_perifocal)
    
    @constraint(model, E - e*sin(E) == M)

    #curtis page 144 & 145
    @constraint(model, cos(E) == (e + cosd(nu)) / (1 + e*cosd(nu)))
    @constraint(model, -tol <= sin(E) - √(1-e^2)*sind(nu) / (1 + e*cosd(nu)) <= tol)

    FullOrbitalParameters(r, v, a, e, i, Ω, ω, nu, M, E)
end
##
model = Model(
    optimizer_with_attributes(Ipopt.Optimizer,
    "max_wall_time" => 30.0)
)
orbparams_i = add_orbital_elements_fix!(model)
ri, vi, ai, ei, ii, Ωi, ωi, nui, Mi, Ei = getfield.(Ref(orbparams_i), fieldnames(FullOrbitalParameters))
orbparams_f = add_orbital_elements_fix!(model)
rf, vf, af, ef, i_f, Ωf, ωf, nuf, Mf, Ef = getfield.(Ref(orbparams_f), fieldnames(FullOrbitalParameters))

@variable(model, Δt)

add_coast_set_boundaries!(model, orbparams_i, orbparams_f, r0, v0, r, v, Δt)
model
##
optimize!(model)
##
solved_ri = value.(ri)
solved_vi = value.(vi)
solved_rf = value.(rf)
solved_vf = value.(vf)
plot_orbit(
    rv_to_kepler(solved_ri, solved_vi),
    rv_to_kepler(solved_rf, solved_vf)
)