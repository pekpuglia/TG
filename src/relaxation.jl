using SatelliteToolboxBase
using SatelliteToolboxPropagators
using JuMP
using Ipopt
using TG
using GLMakie
using LinearAlgebra
##
a1 = 7000e3

orb1 = KeplerianElements(
    date_to_jd(2023, 1, 1, 0, 0, 0),
    a1,
    0.0,
    51 |> deg2rad,
    0    |> deg2rad,
    0     |> deg2rad,
    60     |> deg2rad
)
##
r1, v1 = kepler_to_rv(orb1)
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
tmax = orbital_period(orb1, GM_EARTH)/2
N = 10
dt = tmax/N
sol = Xs(X -> two_body_dyn(X, GM_EARTH), [r1;v1], dt, N, RK4)
##
r = getindex.(sol, Ref(1:3))
r_history = hcat(r...)
##
scatter(r_history[1, :], r_history[2, :], r_history[3, :])
##
rf = r[end]
##
model = Model(Ipopt.Optimizer)

r = @variable(model, [1:3, 1:N+1], base_name="r")

for i = 1:(N+1)
    set_start_value.(r[:, i], r1)
    # @constraint(model, r[:, i]' * r[:, i] >= EARTH_EQUATORIAL_RADIUS^2)
end

v = @variable(model, [1:3, 1:N+1], base_name="v")

@constraint(model, r[:, 1] .== r1)
@constraint(model, r[:, end] .== rf)

for i = 1:N
    # @constraint(model, cross(r[:, i], v[:, i])[3] >= )
    @constraint(model, RK4(X -> two_body_dyn(X, GM_EARTH), [r[:, i]; v[:, i]], dt) .== [r[:, i+1]; v[:, i+1]])
end
##
optimize!(model)
##
solved_r = value.(r)
##
scatter(solved_r[1, :], solved_r[2, :], solved_r[3, :])
