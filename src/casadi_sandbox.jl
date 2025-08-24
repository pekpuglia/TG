include("casadi_interface.jl")
##


x = SX("x")
y = SX("y")

planner = CasADiPlanner([x; y])

add_equality!(planner, x - y, 0.1)
add_inequality!(planner, x, 0, Inf)
add_bounds!(planner, y, 0, 1)

α = 1
b = 100
f = (α - x)^2 + b*(y - x^2)^2

planner.prob["f"] = f

S = casadi.nlpsol("S", "ipopt", planner.prob)

sol, solvetime = solve_planner(S, planner, [0.0,0])

println("Optimal solution: x = ", sol["x"].toarray()[1], ", y = ", sol["x"].toarray()[2])