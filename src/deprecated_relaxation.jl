using SatelliteToolboxBase
using SatelliteToolboxPropagators
using JuMP
using Ipopt
include("TG.jl")
using .TG
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

L = (a1+a2)/2
T = transfer_time

MUPRIME = T^2 * GM_EARTH / L ^ 3

##
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

r1 = r1 / L
v1 = v1 * T / L

r2 = r2 / L
v2 = v2 * T / L
##
#solve lambert problem to use as initial condition
N = 100

model = Model(Ipopt.Optimizer)

r, v = add_coast_segment(model, transfer_time / T, N, "lamb",
    dyn=(X -> two_body_dyn(X, MUPRIME)))

dt10 = 0/3*transfer_time
dt20 = 3/3*transfer_time

r1_lamb, v1_lamb = Propagators.propagate!(prop1, dt10)
r2_lamb = Propagators.propagate!(prop2, dt10+dt20-transfer_time)[1]

@constraint(model, r[:, 1] .== r1_lamb / L)
@constraint(model, r[:, end] .== r2_lamb / L)

#follow natural orbit
for i = 1:size(r)[2]
    rlamb_guess, vlamb_guess = Propagators.propagate!(prop1, dt10 + dt20 * (i/(N+1)))
    set_start_value.(r[:, i], rlamb_guess / L)
    set_start_value.(v[:, i], vlamb_guess * T / L)
end
##
optimize!(model)
##
lamb1 = rv_to_kepler(L * value.(r[:, 1]), L/T * value.(v[:, 1]))
lamb_prop = Propagators.init(Val(:TwoBody), lamb1)
f, ax3d = plot_orbit(orb1, orb2, lamb1)
f
# save("./report/img/hohmann_lambert_guess.png", f)
##
ax3d.azimuth = -π/2
ax3d.elevation = orb1.i
save("./report/img/hohmann_lambert_guess_in_plane.png", f)
##

nimp = 2

N=50
model = Model(
    optimizer_with_attributes(Ipopt.Optimizer,
    "max_iter" => 3_000
    )
)

# dt1 = @variable(model, base_name="dt1", lower_bound=0, upper_bound=transfer_time)
# dt2 = @variable(model, base_name="dt2", lower_bound=0, upper_bound=transfer_time)
# dt3 = @variable(model, base_name="dt3", lower_bound=0, upper_bound=transfer_time)
# @constraint(model, dt1 + dt2 + dt3 == transfer_time)

dts = @variable(model, [1:nimp+1], base_name = "dt", lower_bound=0, upper_bound=transfer_time / T)
@constraint(model, sum(dts) == transfer_time / T)

# deltaV1mag = @variable(model, lower_bound = 0, base_name="dV1mag")
# deltaV2mag = @variable(model, lower_bound = 0, base_name="dV2mag")

# deltaV1dir = @variable(model, [1:3], base_name="dV1dir")
# @constraint(model, deltaV1dir' * deltaV1dir == 1)
# deltaV2dir = @variable(model, [1:3], base_name="dV2dir")
# @constraint(model, deltaV2dir' * deltaV2dir == 1)

deltaVmags = @variable(model, [1:nimp], lower_bound = 0, base_name = "dVmag")
deltaVdirs = @variable(model, [1:3, 1:nimp], base_name = "dVdir")
@constraint(model, [i=1:nimp], deltaVdirs[:, i]' * deltaVdirs[:, i] == 1)

# rcoast1, vcoast1 = add_coast_segment(model, dt1, N, 1)
# rcoast2, vcoast2 = add_coast_segment(model, dt2, N, 2)
# rcoast3, vcoast3 = add_coast_segment(model, dt3, N, 3)

rvcoasts = [add_coast_segment(
    model, dts[i], N, i,
    dyn=(X -> two_body_dyn(X, MUPRIME))) for i = 1:nimp+1]

rcoast = cat(first.(rvcoasts)..., dims=3)
vcoast = cat(last.(rvcoasts)..., dims=3)

@constraint(model, rcoast[:, 1, 1] .== r1)
@constraint(model, vcoast[:, 1, 1] .== v1)

# @constraint(model, rcoast[:, 1, 2] .== rcoast[:, end, 1])
# deltaV1 = deltaV1mag * deltaV1dir
# @constraint(model, vcoast2[:, 1] .== vcoast1[:, end] + deltaV1)


# @constraint(model, rcoast3[:, 1] .== rcoast2[:, end])
# deltaV2 = deltaV2mag * deltaV2dir
# @constraint(model, vcoast3[:, 1] .== vcoast2[:, end] + deltaV2)

#continuity constraints
for i = 1:nimp
    @constraint(model, rcoast[:, 1, i+1] .== rcoast[:, end, i])
    deltaV = deltaVmags[i] * deltaVdirs[:, i]
    @constraint(model, vcoast[:, 1, i+1] .== vcoast[:, end, i] + deltaV)
end


@constraint(model, rcoast[:, end, end] .== r2)
@constraint(model, vcoast[:, end, end] .== v2)

@objective(model, MIN_SENSE, sum(deltaVmags))

#start condition - set to maneuvers along lambert between start and end
#lambert maneuver
set_start_value.(dts, [dt10; dt20; transfer_time-dt10-dt20] / T)
# set_start_value(dt2, dt20)
# set_start_value(dt3, transfer_time - dt10 - dt20)

dv1 = value.(v[:, 1]) - v1
set_start_value(deltaVmags[1], norm(dv1))
set_start_value.(deltaVdirs[1], dv1/norm(dv1))

dv2 = value.(v[:, end]) - v2
set_start_value(deltaVmags[2], norm(dv2))
set_start_value.(deltaVdirs[2], dv2/norm(dv2))

for i = 1:(N+1)
    r1i, v1i = Propagators.propagate!(prop1, dt10*(i-1)/N)
    set_start_value.(rcoast[:, i, 1], r1i / L)
    set_start_value.(vcoast[:, i, 1], v1i * T / L)
    
    r2i, v2i = Propagators.propagate!(lamb_prop, dt20*(i-1)/N)
    set_start_value.(rcoast[:, i, 2], r2i / L)
    set_start_value.(vcoast[:, i, 2], v2i * T / L)
    
    r3i, v3i = Propagators.propagate!(prop2, (transfer_time-dt10-dt20)*((i-1)/N - 1))
    set_start_value.(rcoast[:, i, 3], r3i / L)
    set_start_value.(vcoast[:, i, 3], v3i * T / L)
end
##
optimize!(model)
##
deltaVmag_solved = value.(deltaVmags)
deltaVdir_solved = value.(deltaVdirs)

deltaV_solved = deltaVdir_solved .* deltaVmag_solved'
##
deltaVmodel = objective_value(model) * L / T
##
transf_orb1 = rv_to_kepler(L * value.(rcoast[:, 1, 1]), L / T * value.(vcoast[:, 1, 1]))
transf_orb2 = rv_to_kepler(L * value.(rcoast[:, 1, 2]), L / T * value.(vcoast[:, 1, 2]))
transf_orb3 = rv_to_kepler(L * value.(rcoast[:, 1, 3]), L / T * value.(vcoast[:, 1, 3]))
##
v1n = norm(v1)
v2n = norm(v2)
vtransf1 = √(2*GM_EARTH*(1/a1 - 1/(a1+a2)))
vtransf2 = √(2*GM_EARTH*(1/a2 - 1/(a1+a2)))
deltaV_hohmann = v2n - vtransf2 + vtransf1 - v1n
##
solved_r = L * [value.(rcoast[:, :, 1]) value.(rcoast[:, :, 2]) value.(rcoast[:, :, 3])]
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