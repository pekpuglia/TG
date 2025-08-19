function [sol, solvetime] = solve_planner(planner,x0, options)
%solve_planner Optimize problem described by planner
%   planner: struct created by multiple shooting functions
%   x0: solver initial guess for all variables in vector form
%   options: option struct
    import casadi.*

    solver = nlpsol('solver', 'ipopt', planner.prob,options);

    tic
    sol = solver('x0', x0, 'lbg', planner.tab_lbg, 'ubg', planner.tab_ubg, 'lbx', planner.tab_lbx, 'ubx', planner.tab_ubx);
    solvetime = toc;
end

