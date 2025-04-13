using SatelliteToolboxBase
using SatelliteToolboxPropagators
using JuMP
using Ipopt
using Setfield
using ForwardDiff
include("TG.jl")
using .TG
using LinearAlgebra
## example 
##
model = Model(
    optimizer_with_attributes(Ipopt.Optimizer,
    "max_wall_time" => 30.0)
)
orbparams_i = add_orbital_elements!(model, true)
r, v, a, e, i, Ω, ω, nu, M, E = getfield.(Ref(orbparams_i), fieldnames(FullOrbitalParameters))

# @constraint(model, a*(1-e) == rp)
# @constraint(model, a*(1+e) == ra)

# @constraint(model, i == orb0.i)

# @constraint(model, nu == orb0.f)

# @constraint(model, Ω == orb0.Ω)

# @constraint(model, ω == orb0.ω)

@constraint(model, r .== r0)
@constraint(model, v .== v0)

model
##
optimize!(model)
##
value(a), agiven
##
value(e), egiven
##
value(nu)
##
rad2deg(value(i))
##
value(Ω)
##
value(ω)
##
value(E)
##
r0, value.(r)
##
v0, value.(v)
##
solved_r = value.(r)
solved_v = value.(v)
plot_orbit(
    orb0,
    rv_to_kepler(solved_r, solved_v),
    KeplerianElements(
        orb0.t,
        value(a),
        value(e),
        value(i),
        value(Ω),
        value(ω),
        value(nu)
    )
)