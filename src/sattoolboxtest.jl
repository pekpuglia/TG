using SatelliteToolboxBase
using SatelliteToolboxPropagators
using ForwardDiff
##
orb = KeplerianElements(
                  date_to_jd(2023, 1, 1, 0, 0, 0),
                  7190.982e3,
                  0.001111,
                  98.405 |> deg2rad,
                  100    |> deg2rad,
                  90     |> deg2rad,
                  19     |> deg2rad
)
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
