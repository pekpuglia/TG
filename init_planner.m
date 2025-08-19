function planner = init_planner(tab_variables, n, N)
%INIT_PLANNER Summary of this function goes here
%   Detailed explanation goes here
    planner.prob.x = tab_variables;
    planner.prob.g = [];
    planner.prob.p = [];

    planner.tab_lbg = [];
    planner.tab_ubg = [];

    planner.tab_lbx = -Inf*ones(size(tab_variables));
    planner.tab_ubx = Inf*ones(size(tab_variables));

    planner.n = n;
    planner.N = N;
end