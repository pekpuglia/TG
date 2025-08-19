function planner = add_equality(planner, newg, newval)
%ADD_INEQUALITY Summary of this function goes here
%   Detailed explanation goes here

    assert(all(size(newg) == size(newval)))

    planner.prob.g = [planner.prob.g; newg];
    planner.tab_lbg = [planner.tab_lbg;newval];
    planner.tab_ubg = [planner.tab_ubg; newval];
end

