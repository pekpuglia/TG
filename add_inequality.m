function planner = add_inequality(planner, newg, newlb, newub)
%ADD_INEQUALITY Summary of this function goes here
%   Detailed explanation goes here

    assert(all(size(newg) == size(newlb)))
    assert(all(size(newg) == size(newub)))

    planner.prob.g = [planner.prob.g; newg];
    planner.tab_lbg = [planner.tab_lbg;newlb];
    planner.tab_ubg = [planner.tab_ubg; newub];
end

