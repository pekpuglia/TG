using SatelliteToolboxBase
using SatelliteToolboxPropagators
using ForwardDiff
using GLMakie
using Setfield
##
orb = KeplerianElements(
                  date_to_jd(2023, 1, 1, 0, 0, 0),
                  7190.982e3,
                  0.001111,
                  98.405 |> deg2rad,
                  0    |> deg2rad,
                  90     |> deg2rad,
                  19     |> deg2rad
)
##
N = 100
θ = LinRange(0, 2π, N)

orbit = ((@set orb.f = theta) for theta in θ) .|> kepler_to_rv .|> first |> stack
current_pos, current_velocity = kepler_to_rv(orb)
velocity_arrow_base = current_pos
velocity_arrow_tip = velocity_arrow_base + current_velocity / √sum(current_velocity .^ 2) * 0.4 * √sum(current_pos .^ 2)
velocity_arrow_data = [velocity_arrow_base velocity_arrow_tip] |> Matrix
##
fig = Figure()
ax3d = Axis3(fig[1, 1])
lines!(ax3d, orbit[1, :], orbit[2, :], orbit[3, :], color=:red)
wireframe!(ax3d, Sphere(Point3(0.0), EARTH_EQUATORIAL_RADIUS), color=:cyan, alpha=0.3)
scatter!(ax3d, current_pos, markersize=20)
lines!(ax3d, velocity_arrow_data[1, :], velocity_arrow_data[2, :], velocity_arrow_data[3, :])
fig
##
orbp = Propagators.init(Val(:TwoBody), orb)
##
Propagators.propagate!(orbp, 50.0, OrbitStateVector)
orbp.tbd.orbk
##
# v_t(t) = Propagators.propagate!(orbp, t, OrbitStateVector).v
# ##
# ForwardDiff.derivative(v_t, 5.0)
##
function f_els(el_vector)
    r, e, i, RAAN, omega, f, time = el_vector
    # time = el_vector[1]

    orb = KeplerianElements(
        date_to_jd(2023, 1, 1, 0, 0, 0),
        r,
        e,
        i,
        RAAN,
        omega,
        f
    )
    orbp = Propagators.init(Val(:TwoBody), orb; propagation_type=Number)
    Propagators.propagate!(orbp, time)
    [
        orbp.tbd.orbk.a
        orbp.tbd.orbk.e
        orbp.tbd.orbk.i
        orbp.tbd.orbk.Ω
        orbp.tbd.orbk.ω
        orbp.tbd.orbk.f
    ]
end
##
ForwardDiff.jacobian(f_els, [
    7190.982e3,
    0.001111,
    98.405 |> deg2rad,
    100    |> deg2rad,
    90     |> deg2rad,
    19     |> deg2rad,
    5.0])
