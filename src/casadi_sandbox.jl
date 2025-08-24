using CasADi
##
struct CasADiPlanner
    prob::Dict
    tab_lbg
    tab_ubg
    tab_lbx
    tab_ubx
    function CasADiPlanner(variables)
        new(Dict(
            "x" => vcat(variables),
            "f" => 0,
            "g" => []
        ),
        [],
        [],
        -Inf*ones(size(variables)),
        Inf*ones(size(variables))
        )
    end
end

function add_equality!(planner::CasADiPlanner, expr, val)
    planner.prob["g"] = vcat([planner.prob["g"]..., expr])
    push!(planner.tab_lbg, val)
    push!(planner.tab_ubg, val)
end

function add_inequality!(planner::CasADiPlanner, expr, lb, ub)
    planner.prob["g"] = vcat([planner.prob["g"]..., expr])
    push!(planner.tab_lbg, lb)
    push!(planner.tab_ubg, ub)
end

sx_iterator(sx) = (sx[i] for i = 1:length(sx))

function add_bounds!(planner::CasADiPlanner, var, lb, ub)
    var_index = findfirst(v -> casadi.is_equal(var, v), sx_iterator(planner.prob["x"]))
    planner.tab_lbx[var_index] = lb
    planner.tab_ubx[var_index] = ub
end

function solve_planner(nlpsol, planner::CasADiPlanner, guess)
    stats = @timed nlpsol(x0 = guess, lbg=planner.tab_lbg, ubg=planner.tab_ubg, lbx=planner.tab_lbx, ubx=planner.tab_ubx)
    stats.value, stats.time
end

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