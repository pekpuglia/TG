function cvar(name, sizes...)
    sxvar = SX(name, sizes...)
    tabvar = reshape([sx_iterator(sxvar)...], sizes...)
end

using CasADi
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
        Float64[],
        Float64[],
        -Inf*ones(size(variables)),
        Inf*ones(size(variables))
        )
    end
end

function add_equality!(planner::CasADiPlanner, expr::SX, val::Number)
    planner.prob["g"] = vcat([planner.prob["g"]..., expr])
    push!(planner.tab_lbg, val)
    push!(planner.tab_ubg, val)
end

function add_equality!(planner::CasADiPlanner, expr, val::AbstractVector)
    for i = 1:length(val)
        add_equality!(planner, expr[i], val[i])
    end
end

# function add_equality!(planner::CasADiPlanner, expr::Array, val::AbstractVector)
#     #implement
#     planner.prob["g"] = vcat([planner.prob["g"]..., expr...])
#     append!(planner.tab_lbg, val)
#     append!(planner.tab_ubg, val)
# end

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
    # stats = @timed nlpsol(x0 = guess, lbg=planner.tab_lbg, ubg=planner.tab_ubg, lbx=planner.tab_lbx, ubx=planner.tab_ubx)
    # stats.value, stats.time
    nlpsol(x0 = guess, lbg=planner.tab_lbg, ubg=planner.tab_ubg, lbx=planner.tab_lbx, ubx=planner.tab_ubx)
end