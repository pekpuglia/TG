using SatelliteToolboxBase
using SatelliteToolboxPropagators
using JuMP
using Ipopt
using Setfield
using ForwardDiff
include("../TG.jl")
using .TG
using LinearAlgebra
## example
orb = KeplerianElements(
    date_to_jd(2023, 1, 1, 0, 0, 0),
    9.6000e+06,
    0.0100,
    0.6981,
    0,
    1.2217,
    5.4105
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

# @constraint(model, a == orb.a)
# @constraint(model, e == orb.e)
# @constraint(model, i == orb.i)
# @constraint(model, Ω == orb.Ω)
# @constraint(model, ω == orb.ω)
# @constraint(model, f == orb.f)
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
    rv_to_kepler(solved_r, solved_v),
)