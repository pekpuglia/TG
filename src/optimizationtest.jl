using SatelliteToolboxBase
using SatelliteToolboxPropagators
using Optimization
using ForwardDiff
##
#example: find time so that 100s after f = 60deg
orb = KeplerianElements(
                  date_to_jd(2023, 1, 1, 0, 0, 0),
                  7190.982e3,
                  0.001111,
                  98.405 |> deg2rad,
                  100    |> deg2rad,
                  90     |> deg2rad,
                  19     |> deg2rad
)

function f_after_time(t, p)
    orb = p[1]
    target_f = p[2]
    orbp = Propagators.init(Val(:TwoBody), orb, propagation_type=Number)
    Propagators.propagate!(orbp, t[1])
    (orbp.tbd.orbk.f - target_f)^2
end
##
f = OptimizationFunction(f_after_time, Optimization.AutoForwardDiff())
##
t0 = [1.0]
p = [orb, 60.0 |> deg2rad]
##
prob = OptimizationProblem(f, t0, p)
##
using OptimizationOptimJL
##
sol = solve(prob, GradientDescent())