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
Propagators.propagate!(orbp, 5.0, OrbitStateVector)
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
    orbp = Propagators.init(Val(:TwoBody), orb)
    Propagators.propagate!(orbp, time)
    orbp.tbd.orbk.f
end
##
ForwardDiff.gradient(f_els, [
    7190.982e3,
    0.001111,
    98.405 |> deg2rad,
    100    |> deg2rad,
    90     |> deg2rad,
    19     |> deg2rad,
    5.0])
