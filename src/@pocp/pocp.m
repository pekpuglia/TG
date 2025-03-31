function ocp = pocp(varargin)
%@POCP/POCP - Create a polynomial optimal control problem
%
%   Consider the cost function
%   
%   J(0, T, x(0), u(t)) = Integral 0->T {h(t, x(t), u(t))} dt + H(x(T))
%   
%   and the continuous-time system
%
%   dx(t)/dt = f(t, x, u)
%
%   Two different kinds of optimal control problems can be modeled by
%   defining a POCP object. The first one is
%
%      min       J(0, T, x(0), u(t))
%    x(0),u(t)
%                dx(t)/dt = f(t, x, u)
%                x(0) in Ci
%                (t, x(t), u(t)) in Ct
%                x(T) in Cf
%
%   To model this problem with POCP, the functions f, h, and H must be
%   polynomial functions and the sets Ci, Ct, and Cf must be semialgebraic.
%   Another class of problems which can be considered are of the form
%
%      min       Integral {J(0, T, x0, u(t, x0))} dmu(x0)
%    u(t,x0)
%                dx(t)/dt = f(t, x, u)
%                (t, x(t), u(t)) in Ct
%                x(T) in Cf
%   
%   where mu(x0) is an assigned probability measures.
%   Although in the problems above the horizon T is given, also problems
%   with free horizon can be considered. For more details on what problems
%   can be modeled with POCP check the PDF documentation which is available
%   with this software.
%
%   The syntax of the pocp command is the following:
%
%   [OCP]=POCP() creates an empty POCP object.
%
%   [OCP]=POCP('ParameterName1', 'Argument1', 'Argument2', ...,
%      'ParameterName2', ...) creates a POCP object with the specified
%      data. The syntax is the same as the function pocp/set. For more
%      information on how the specify the parameters type 'help pocp/set'
%      and press Enter.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2008 by Carlo Savorgnan 
% 
% This program is free software; you can redistribute it and/or modify 
% it under the terms of the GNU General Public License version 2 as 
% published by the Free Software Foundation. 
% 
% This program is distributed in the hope that it will be useful, 
% but WITHOUT ANY WARRANTY; without even the implied warranty of 
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
% GNU General Public License for more details. 
% 
% You should have received a copy of the GNU General Public License 
% along with this program; if not, write to the 
% Free Software Foundation, Inc., 
% 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% data entered by the user
ocp.state = [];    % state variables
ocp.input = [];    % input variables
ocp.time = [];     % time variable
ocp.timeset = false;  % indicates if the user set the time variable or not
ocp.dynamics = []; % system dynamics
ocp.horizon = 0;  % ocp horizon (zero when the horizon is not fixed)
ocp.discrete = false;   % false if the dynamics is continuous, true if it is discrete
ocp.testtime = [];  % [] default behaviour, false do not test time, true test time

ocp.icon = [];  % opc constraints on the initial condition
ocp.icon.ineq = [];
ocp.icon.dirac.var = [];
ocp.icon.dirac.value = [];
ocp.icon.dirac.weight = [];
ocp.icon.unif.var = [];
ocp.icon.unif.interval = [];

ocp.fcon = [];  % opc constraints on the final condition
ocp.fcon.ineq = [];
ocp.fcon.dirac.var = [];
ocp.fcon.dirac.value = [];
ocp.fcon.dirac.weight = [];

ocp.tcon = [];  % opc constraints on the trajectories
ocp.tcon.ineq = [];    % constraints on the istantaneous value
ocp.tcon.int = [];    % integral constraint (of type momcon)

ocp.intcost = [];  % opc integral cost
ocp.fcost = [];  % opc final cost

% variables introduced to solve the ocp (not entered by the user)
ocp.istate = [];   % state at initial time (remains empty if the statistics of the initial state is completely assigned)
ocp.fstate = [];   % state at final time (remains empty if the statistics of the final state is completely assigned)
ocp.ftime = [];    % time at final time (equal to ocp.horizon when the horizon is fixed) a mpol variable otherwise

% variables which contain the measures
ocp.imeas = []; % occupation measure at initial time (remains empty if the statistics of the initial state is completely assigned)
ocp.fmeas = []; % occupation measure at final time (remains empty if the statistics of the initial state is completely assigned)
ocp.tmeas = []; % occupation measure for the trajectory

% variables used to define the gloptipoly problem
ocp.supcon = [];
ocp.momcon = [];
ocp.gloptipoly = [];

% other
ocp.ntestfun = [];  %number of test function used
ocp.iassigned = []; %true is the statistics of the variables at time 0 is completely assigned
ocp.fassigned = []; %true is the statistics of the variables at final time is completely assigned
ocp.cost = [];
ocp.testfun = []; %vector containing the test functions

ocp.sol.cost = [];

ocp = class(ocp, 'pocp');

if nargin~=0
    cmd = 'ocp = set(ocp';
    for index = 1:nargin
        cmd = sprintf('%s, varargin{%d}', cmd, index);
    end
    cmd = sprintf('%s);', cmd);
    eval(cmd);
end

