using SatelliteToolboxBase
using SatelliteToolboxPropagators
using JuMP
using Ipopt
using TG
using GLMakie
using LinearAlgebra
using Setfield

#TODO
#do entire primer vector trajectory
#implement interactive primer vector
#do multiple impulses
#free initial and final states
##
a1 = 7000e3
a2 = 40000e3

hohmann_time = orbital_period((a1+a2)/2, GM_EARTH) / 2

transfer_time = 2hohmann_time

orb1 = KeplerianElements(
    date_to_jd(2023, 1, 1, 0, 0, 0),
    a1,
    0.0,
    51 |> deg2rad,
    0    |> deg2rad,
    0     |> deg2rad,
    0     |> deg2rad
)
prop1 = Propagators.init(Val(:TwoBody), orb1)

orb2 = KeplerianElements(
    orb1.t + transfer_time / 86400,
    a2,
    0.2,
    orb1.i,
    orb1.Ω + deg2rad(31),
    orb1.ω,
    deg2rad(180) + orb1.f
)

hohmann_orb = KeplerianElements(
    orb1.t,
    (a1+a2)/2,
    (a2-a1) / (a1+a2),
    orb1.i,
    orb1.Ω,
    orb1.ω,

    deg2rad(60)
)

prop2 = Propagators.init(Val(:TwoBody), orb2)
f, ax3d = plot_orbit(orb1, orb2, hohmann_orb)
f
##
save("./report/img/hohmann_condition.png", f)
##
ax3d.azimuth = -π/2
ax3d.elevation = orb1.i
save("./report/img/hohmann_condition_in_plane.png", f)
##
r1, v1 = kepler_to_rv(orb1)
r2, v2 = kepler_to_rv(orb2)

##
#solve lambert problem to use as initial condition
N = 100

model = Model(Ipopt.Optimizer)

r, v = add_coast_segment(model, transfer_time, N, "lamb")

dt10 = 0/3*transfer_time
dt20 = 3/3*transfer_time

r1_lamb, v1_lamb = Propagators.propagate!(prop1, dt10)
r2_lamb = Propagators.propagate!(prop2, dt10+dt20-transfer_time)[1]

@constraint(model, r[:, 1] .== r1_lamb)
@constraint(model, r[:, end] .== r2_lamb)

#follow natural orbit
for i = 1:size(r)[2]
    rlamb_guess, vlamb_guess = Propagators.propagate!(prop1, dt10 + dt20 * (i/(N+1)))
    set_start_value.(r[:, i], rlamb_guess)
    set_start_value.(v[:, i], vlamb_guess)
end
##
optimize!(model)
##
lamb1 = rv_to_kepler(value.(r[:, 1]), value.(v[:, 1]))
lamb_prop = Propagators.init(Val(:TwoBody), lamb1)
f, ax3d = plot_orbit(orb1, orb2, lamb1)
f
# save("./report/img/hohmann_lambert_guess.png", f)
##
ax3d.azimuth = -π/2
ax3d.elevation = orb1.i
save("./report/img/hohmann_lambert_guess_in_plane.png", f)
##
N=50
model = Model(
    optimizer_with_attributes(Ipopt.Optimizer,
    "max_iter" => 10_000
    )
)

dt1 = @variable(model, base_name="dt1", lower_bound=0, upper_bound=transfer_time)
dt2 = @variable(model, base_name="dt2", lower_bound=0, upper_bound=transfer_time)
dt3 = @variable(model, base_name="dt3", lower_bound=0, upper_bound=transfer_time)
@constraint(model, dt1 + dt2 + dt3 == transfer_time)

@constraint(model, dt1 == 0)
@constraint(model, dt3 == 0)

@constraint(model, dt1+dt2 <= transfer_time)

deltaV1mag = @variable(model, lower_bound = 0, base_name="dV1mag")
deltaV2mag = @variable(model, lower_bound = 0, base_name="dV2mag")

deltaV1dir = @variable(model, [1:3], base_name="dV1dir")
@constraint(model, deltaV1dir' * deltaV1dir == 1)
deltaV2dir = @variable(model, [1:3], base_name="dV2dir")
@constraint(model, deltaV2dir' * deltaV2dir == 1)


rcoast1, vcoast1 = add_coast_segment(model, dt1, N, 1)
rcoast2, vcoast2 = add_coast_segment(model, dt2, N, 2)
rcoast3, vcoast3 = add_coast_segment(model, dt3, N, 3)

@constraint(model, rcoast1[:, 1] .== r1)
@constraint(model, vcoast1[:, 1] .== v1)

@constraint(model, rcoast2[:, 1] .== rcoast1[:, end])
deltaV1 = deltaV1mag * deltaV1dir
@constraint(model, vcoast2[:, 1] .== vcoast1[:, end] + deltaV1)


@constraint(model, rcoast3[:, 1] .== rcoast2[:, end])
deltaV2 = deltaV2mag * deltaV2dir
@constraint(model, vcoast3[:, 1] .== vcoast2[:, end] + deltaV2)


@constraint(model, rcoast3[:, end] .== r2)
@constraint(model, vcoast3[:, end] .== v2)

@objective(model, MIN_SENSE, deltaV1mag + deltaV2mag)

#start condition - set to maneuvers along lambert between start and end
#lambert maneuver
# set_start_value(dt1, dt10)
# set_start_value(dt2, dt20)
# set_start_value(dt2, transfer_time - dt10 - dt20)

dv1 = value.(v[:, 1]) - v1
set_start_value(deltaV1mag, norm(dv1))
set_start_value.(deltaV1dir, dv1/norm(dv1))

dv2 = value.(v[:, end]) - v2
set_start_value(deltaV2mag, norm(dv2))
set_start_value.(deltaV2dir, dv2/norm(dv2))

for i = 1:(N+1)
    r1i, v1i = Propagators.propagate!(prop1, dt10*(i-1)/N)
    set_start_value.(rcoast1[:, i], r1i)
    set_start_value.(vcoast1[:, i], v1i)
    
    r2i, v2i = Propagators.propagate!(lamb_prop, dt20*(i-1)/N)
    set_start_value.(rcoast2[:, i], r2i)
    set_start_value.(vcoast2[:, i], v2i)
    
    r3i, v3i = Propagators.propagate!(prop2, (transfer_time-dt10-dt20)*((i-1)/N - 1))
    set_start_value.(rcoast3[:, i], r3i)
    set_start_value.(vcoast3[:, i], v3i)
end
##
optimize!(model)
##
deltaV1_solved = value.(deltaV1)
deltaV2_solved = value.(deltaV2)
##
deltaVmodel = objective_value(model)
##
transf_orb1 = rv_to_kepler(value.(rcoast1[:, 1]), value.(vcoast1[:, 1]))
transf_orb2 = rv_to_kepler(value.(rcoast2[:, 1]), value.(vcoast2[:, 1]))
transf_orb3 = rv_to_kepler(value.(rcoast3[:, 1]), value.(vcoast3[:, 1]))
##
v1n = norm(v1)
v2n = norm(v2)
vtransf1 = √(2*GM_EARTH*(1/a1 - 1/(a1+a2)))
vtransf2 = √(2*GM_EARTH*(1/a2 - 1/(a1+a2)))
deltaV_hohmann = v2n - vtransf2 + vtransf1 - v1n
##
solved_r = [value.(rcoast1) value.(rcoast2) value.(rcoast3)]
##
# f = Figure()
# ax3d = Axis3(f[1, 1])

f, ax3d = plot_orbit(orb1, orb2)

scatter!(ax3d, solved_r[1, :], solved_r[2, :], solved_r[3, :], color="green")
scatter!(ax3d, solved_r[1, (N+1):(N+1):end-1], solved_r[2, (N+1):(N+1):end-1], solved_r[3, (N+1):(N+1):end-1], markersize=20, color=:red)
wireframe!(ax3d, Sphere(Point3(0.0), EARTH_EQUATORIAL_RADIUS), color=:cyan, alpha=0.3)

f
##
save("./report/img/hohmann_solved.png", f)
ax3d.azimuth = -π/2
ax3d.elevation = orb1.i
save("./report/img/hohmann_solved_in_plane.png", f)
##
f, _ = plot_orbit(orb1, orb2, transf_orb1, transf_orb2, transf_orb3)
f
## primer vector 

p0 = deltaV1_solved / norm(deltaV1_solved)
pf = deltaV2_solved / norm(deltaV2_solved)

inter_impulse_prop = Propagators.init(Val(:TwoBody), transf_orb2)

p0dot = p0dot_tpbvp(p0, pf, value(dt2), inter_impulse_prop)
## plot p and pdot
tspan = range(0, value(dt2), 100)
##
ppdot = [TG.Phi_time(inter_impulse_prop, t) * [p0; p0dot] for t in tspan]
##
normpdot = [dot(ppdoti[1:3], ppdoti[4:6]) / norm(ppdoti[1:3]) for ppdoti in ppdot]
##
f = Figure()
ax1 = Axis(f[1, 1], xlabel = "t (s)", ylabel = "|p|")
lines!(ax1, tspan, norm.(getindex.(ppdot, Ref(1:3))))

ax2 = Axis(f[2, 1], xlabel = "t (s)", ylabel = L"d |p| / dt")
lines!(ax2, tspan, normpdot)
f
# save("./report/img/hohmann_primer_vector_history.png", f)