using SatelliteToolboxBase
using SatelliteToolboxPropagators
using JuMP
using Ipopt
using TG
using GLMakie
using LinearAlgebra
##
a1 = 7000e3
a2 = 9000e3

hohmann_time = orbital_period((a1+a2)/2, GM_EARTH) / 2

transfer_time = hohmann_time

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
    0.0,
    orb1.i,
    orb1.Ω,
    orb1.ω,
    deg2rad(180) + orb1.f
)
prop2 = Propagators.init(Val(:TwoBody), orb2)
f = plot_orbit(orb1, orb2)
save("./report/img/hohmann_condition.png", f)
##
r1, v1 = kepler_to_rv(orb1)
r2, v2 = kepler_to_rv(orb2)
##
function two_body_dyn(X, mu)
    r = X[1:3]
    v = X[4:6]
    [
        v
        -mu/(√(r'*r)^3)*r
    ]
end

function euler(f, X, dt)
    X + f(X)*dt
end

function RK4(f, X, dt)
    k1 = dt*f(X)
    k2 = dt*f(X+k1/2)
    k3 = dt*f(X+k2/2)
    k4 = dt*f(X+k3)
    X + (k1+2k2+2k3+k4)/6
end

function Xs(f, X0, dt, N, integrator=euler)
    sol = [X0]
    for i = 1:N
        push!(sol, integrator(f, sol[i], dt))
    end
    sol
end
##
function add_coast_segment(model, deltat, N, ind)
    #[(x, y, z), time instants]
    r = @variable(model, [1:3, 1:N+1], base_name="r_coast_$ind")
    v = @variable(model, [1:3, 1:N+1], base_name="v_coast_$ind")


    for i = 1:(N+1)
        # @constraint(model, r[:, i]' * r[:, i] >= EARTH_EQUATORIAL_RADIUS^2)
        @constraint(model, cross(r[:, i], v[:, i])[3] >=0)
    end

    for i = 1:N
        @constraint(model, RK4(X -> two_body_dyn(X, GM_EARTH), [r[:, i]; v[:, i]], deltat/N) .== [r[:, i+1]; v[:, i+1]])
    end

    r, v
end
##
#solve lambert problem to use as initial condition
N = 10

model = Model(Ipopt.Optimizer)

r, v = add_coast_segment(model, transfer_time, N, "lamb")

dt10 = 1/3*transfer_time
dt20 = 1/3*transfer_time

r1_lamb, v1_lamb = Propagators.propagate!(prop1, dt10)
r2_lamb = Propagators.propagate!(prop2, dt10+dt20-transfer_time)[1]

@constraint(model, r[:, 1] .== r1_lamb)
@constraint(model, r[:, end] .== r2_lamb)

for i = 1:size(r)[2]
    set_start_value.(r[:, i], r1_lamb)
    set_start_value.(v[:, i], v1_lamb)
end
##
optimize!(model)
##
lamb1 = rv_to_kepler(value.(r[:, 1]), value.(v[:, 1]))
lamb_prop = Propagators.init(Val(:TwoBody), lamb1)
f = plot_orbit(orb1, orb2, lamb1)
save("./report/img/hohmann_lambert_guess.png", f)
##
N=10
model = Model(
    optimizer_with_attributes(Ipopt.Optimizer,
    "max_iter" => 10_000
    )
)

dt1 = @variable(model, base_name="dt1", lower_bound=0, upper_bound=transfer_time)
dt2 = @variable(model, base_name="dt2", lower_bound=0, upper_bound=transfer_time)

@constraint(model, dt1+dt2 <= transfer_time)

deltaV1mag = @variable(model, lower_bound = 0, base_name="dV1mag")
deltaV2mag = @variable(model, lower_bound = 0, base_name="dV2mag")

deltaV1dir = @variable(model, [1:3], base_name="dV1dir")
@constraint(model, deltaV1dir' * deltaV1dir == 1)
deltaV2dir = @variable(model, [1:3], base_name="dV2dir")
@constraint(model, deltaV2dir' * deltaV2dir == 1)


rcoast1, vcoast1 = add_coast_segment(model, dt1, N, 1)
rcoast2, vcoast2 = add_coast_segment(model, dt2, N, 2)
rcoast3, vcoast3 = add_coast_segment(model, transfer_time-dt1-dt2, N, 3)

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

#start condition
#lambert maneuver
set_start_value(dt1, dt10)
set_start_value(dt2, dt20)

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
f = Figure()
ax3d = Axis3(f[1, 1])

scatter!(ax3d, solved_r[1, :], solved_r[2, :], solved_r[3, :])
scatter!(ax3d, solved_r[1, (N+1):(N+1):end], solved_r[2, (N+1):(N+1):end], solved_r[3, (N+1):(N+1):end], markersize=20, color=:red)
wireframe!(ax3d, Sphere(Point3(0.0), EARTH_EQUATORIAL_RADIUS), color=:cyan, alpha=0.3)

f
save("./report/img/hohmann_solved.png", f)
##
plot_orbit(orb1, orb2, transf_orb1, transf_orb2, transf_orb3)
## primer vector 

p0 = deltaV1_solved / norm(deltaV1_solved)
pf = deltaV2_solved / norm(deltaV2_solved)

inter_impulse_prop = Propagators.init(Val(:TwoBody), transf_orb2)
Phi = TG.Phi_time(inter_impulse_prop, value(dt2))

M = Phi[1:3, 1:3]
N = Phi[1:3, 4:6]

if abs(det(N)) <= 1e-10
    #CHECK THIS IS TRUE
    #only works for coplanar transfers?
    A = N * [r1 v1]
    b = pf - M * p0
    p0dot = [r1 v1] * (A \ b)
else
    p0dot = N \ (pf - M * p0)
end
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

save("./report/img/hohmann_primer_vector_history.png", f)