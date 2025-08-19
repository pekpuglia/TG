function [planner] = fixed_time_planning(n,m,N,T,f_dyn)
%n: dimension of the state
%m: dimension of the control
%N: number of time samples
%T: total time
%R_input: running weight on the inputs
%Q_state: running weight on the state
%Q_f_state weights on the final state
%f_dyn: continuous time dynamics as a function on the state and controls

% This function creates the optimization problem and all of the internal
% constraints, i.e. integration constraints. Boundary conditions are to be
% specified separately with the add_equality or add_inequality functions.

import casadi.*

%control trajectory
tab_u=SX.sym('tab_u',m,N-1);
%state  trajectory
tab_X= SX.sym('tab_X',n,N);

tab_X_U=[
    reshape(tab_X,n*N,1);
    reshape(tab_u,m*(N-1),1)
];%vectorized version of the optimisation variable

planner = init_planner(tab_X_U, n, m, N);

cost_dyn = @(xaug, u) [ ...
    f_dyn(xaug(1:n), u); ...
    (1/2)*(u'*R_input*u)+ (1/2)*(xaug(1:n)'*Q_state*xaug(1:n)) ...
];

% state and cost integrator
F = RK4_multistep(cost_dyn, n+1, m, 4);

g=[];% general inequality and equality constraints
J = 0;%total cost
% Formulate the NLP
for k=1:N-1
    % Integrate till the end of the interval
    Xk=tab_X(:,k);
    Fk = F('x0',[Xk; J], 'dt', T/(N-1),'u0',tab_u(:,k));
    
    Xk_end = Fk.xf(1:n);%prediction of the state
    J = Fk.xf(end);%new cumulated cost
    
    
    % Add equality constraint
    g = [
        g; 
        Xk_end-tab_X(:,k+1)
    ];
end
planner = add_equality(planner, g, zeros(n*(N-1), 1));

J = J + tab_X(:,N)' * Q_f_state * tab_X(:,N);

%Casadi problem has fields f, x, g, p
planner.prob.f = J;
planner.T = T;
planner.fdyn = f_dyn;
planner.R = R_input;
planner.Q = Q_state;
planner.Qf = Q_f_state;
planner.tabX = tab_X;
planner.tabu = tab_u;
end

