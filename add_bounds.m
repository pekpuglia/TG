function planner = add_bounds(planner, var,lb, ub)
%ADD_BOUNDS Add variable bounds to specific variable
% Can only do one variable at a time
%   var: CasADi variable to bound (should be contained in planner.prob.x)

    var_index = find(arrayfun(@(x) is_equal(x, var), planner.prob.x));

    planner.tab_lbx(var_index) = lb;
    planner.tab_ubx(var_index) = ub;
end

