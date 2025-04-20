using SatelliteToolboxBase
using SatelliteToolboxPropagators
using JuMP
using Ipopt
using Setfield
using ForwardDiff
include("../TG.jl")
using .TG
using LinearAlgebra
##
orb = KeplerianElements(
    date_to_jd(2023, 1, 1, 0, 0, 0),
    7.0000e+06,
    0.000,
    0.6981,
    4.0143,
    0.6981,
    4.1888
)
given_r, given_v = kepler_to_rv(orb)
plot_orbit(orb)
##
model = Model(
    optimizer_with_attributes(Ipopt.Optimizer,
    "max_wall_time" => 30.0)
)
orbparams = add_orbital_elements!(model)
r, v, a, e, i, Ω, ω, f, M, E = getfield.(Ref(orbparams), fieldnames(FullOrbitalParameters))
        
@constraint(model, r .== given_r)
@constraint(model, v .== given_v)
model
##
optimize!(model)
##
value(a)
##
value(e)
##
value(i)
##
#oK!
value(Ω)
##
value(ω)
##
value(f)
##
solved_r = value.(r)
solved_v = value.(v)

plot_orbit(
    orb,
    KeplerianElements(
        date_to_jd(2023, 1, 1, 0, 0, 0),
        value(a),
        clamp(value(e), 0, 1),
        value(i),
        value(Ω),
        value(ω),
        value(f)
    ),
    rv_to_kepler(solved_r, solved_v),
)