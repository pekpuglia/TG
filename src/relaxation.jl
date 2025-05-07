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
    60     |> deg2rad
)

orb2 = KeplerianElements(
    orb1.t + transfer_time / 86400,
    a2,
    0.0,
    orb1.i,
    orb1.Ω,
    orb1.ω,
    deg2rad(180) + orb1.f
)
plot_orbit(orb1, orb2)
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
        @constraint(model, r[:, i]' * r[:, i] >= EARTH_EQUATORIAL_RADIUS^2)
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

@constraint(model, r[:, 1] .== r1)
@constraint(model, r[:, end] .== r2)

for i = 1:size(r)[2]
    set_start_value.(r[:, i], r1)
end
##
optimize!(model)
##
lamb1 = rv_to_kepler(value.(r[:, 1]), value.(v[:, 1]))
plot_orbit(orb1, orb2, lamb1)
##
N=10
model = Model(
    optimizer_with_attributes(Ipopt.Optimizer,
    "max_iter" => 10_000
    )
)

t1 = @variable(model, base_name="t1", lower_bound=0, upper_bound=transfer_time)
t2 = @variable(model, base_name="t2", lower_bound=0, upper_bound=transfer_time)

@constraint(model, t1+t2 <= transfer_time)

deltaV1mag = @variable(model, lower_bound = 0, base_name="dV1mag")
deltaV2mag = @variable(model, lower_bound = 0, base_name="dV2mag")

deltaV1dir = @variable(model, [1:3], base_name="dV1dir")
@constraint(model, deltaV1dir' * deltaV1dir == 1)
deltaV2dir = @variable(model, [1:3], base_name="dV2dir")
@constraint(model, deltaV2dir' * deltaV2dir == 1)


rcoast1, vcoast1 = add_coast_segment(model, t1, N, 1)
rcoast2, vcoast2 = add_coast_segment(model, t2, N, 2)
rcoast3, vcoast3 = add_coast_segment(model, transfer_time-t1-t2, N, 3)

for i = 1:(N+1)
    #generate initial condition w lambert?
    set_start_value.(rcoast1[:, i], r1)
    set_start_value.(rcoast2[:, i], r1)
    set_start_value.(rcoast3[:, i], r2)
end


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
##
optimize!(model)
##
solved_r = [value.(rcoast1) value.(rcoast2) value.(rcoast3)]
##
f = Figure()
ax3d = Axis3(f[1, 1])

scatter!(ax3d, solved_r[1, :], solved_r[2, :], solved_r[3, :])
scatter!(ax3d, solved_r[1, (N+1):(N+1):end], solved_r[2, (N+1):(N+1):end], solved_r[3, (N+1):(N+1):end], markersize=20, color=:red)

f