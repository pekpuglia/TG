function [status, cost, measures, varargout] = solvepocp(ocp, varargin)
%@POCP/SOLVEPOCP - Solve polynomial optimal control problem
%
%   [STATUS, COST, MEASURES] = SOLVEPOCP(PROB, 'mo', DEG) solves the polynomial
%      optimal control problem PROB. Moments up to degree DEG (if DEG is even, or
%      the next even number DEG it is odd) will be used to represent the trajectory
%      occupation measure.
%   [STATUS, COST, MEASURES] = SOLVEPOCP(PROB, 'tf', DEG) solves the polynomial
%      optimal control problem PROB. Test functions up to degree DEG will
%      be used.
%   [STATUS, COST, MEASURES] = SOLVEPOCP(PROB, DEG) when no 'mo' or 'tf' is
%      specified the default behaviour is used ('om').
%
%   Output parameters:
%      -STATUS can take the following values: -1, 0, +1. In the former case
%       the problem is unfeasible or could not be solved. When STATUS is 0
%       the problem could be solved. When STATUS is +1, the problem could be solved
%       and GloptiPoly found an optimal solution to the moment problem
%       formulated (which doesn't imply that the optimal control problem
%       has been solved exactly).
%      -COST suboptimal value of the cost of the optimal control problem.
%      -MEASURES is a structure containing the measures used by Gloptipoly
%       to solve the moment problem greated. If the initial or final
%       conditions have been assigned completely, the corresponding measure
%       is empty.
%
%   [STATUS, COST, MEASURES, VF] = SOLVEPOCP(PROB, ...) when four output
%      arguments are specified a polynomial subsolution of the Hamilton
%      Jacobi Bellman equation is also calculated. The polynomial will be
%      stored in VF (value function). If the option 'tf' is used, DEG
%      coresponds to the maximal degree of VF.
%
%   NOTE: when four output arguments are specified Gloptipoly is used in
%   combination with YALMIP (which should be of course installed in this case).

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

global MMM;

% checking if state and dynamics have been set
if isempty(ocp.state) || isempty(ocp.dynamics)
    error('One or more mandatory fields of the POCP object are empty.');
end

% determining the degree of the test function to be used and parsing the
% input
if nargin==2
    degree = varargin{1};
    if ~isnumeric(degree) || (degree~=round(degree)) || degree<1
        error('When two input arguments are specified, the second should be a positive integer.');
    end
    degree = 2*ceil(degree/2);
    degree = degree-deg(ocp.dynamics)+1;
    if degree<=0
        error('The moment order specified is too low to construct the moment problem.');
    end
elseif nargin==3
    mode = varargin{1};
    if ~ischar(mode) || (~strcmp(mode, 'tf') && ~strcmp(mode, 'mo'))
        error('When three input arguments are specified the second one must be either ''tf'' or ''mo''.');
    end
    degree = varargin{2};
    if ~isnumeric(degree) || (degree~=round(degree)) || degree<1
        error('When three input arguments are specified, the third one must be a positive integer.');
    end
    if strcmp(mode, 'mo')
        degree = 2*ceil(degree/2);
        degree = degree-deg(ocp.dynamics)+1;
        if degree<=0
            error('The moment order specified is too low to construct the moment problem.');
        end
    end 
else
    error('Two or three input arguments must be specified.');
end

% checking for some common mistakes and setting some dafaults variables
if isempty(ocp.icon.ineq) && isempty(ocp.icon.dirac.var) && isempty(ocp.icon.unif.var) && ocp.horizon==0
    disp('Warning: No constraints on the initial condition have been specified. Since the');
    disp('         horizon is free, this would cause the solver to find a solution corresponding');
    disp('         to a zero trajectory.');
    disp('         To avoid this situation the initial condition is assigned automatically');
    disp('         (a uniform distribution in the interval [-1, 1] is set for every state');
    disp('         variable).');
    disp('         THE ASSIGNMENT IS ONLY TEMPORARILY.');
    pause(1);
    nstate=length(ocp.state);
    ocp=setiuniform(ocp, ocp.state, [-ones(nstate, 1), ones(nstate, 1)]);
end

if isempty(ocp.fcon.ineq) && isempty(ocp.fcon.dirac.var) && ocp.horizon==0 && isempty(ocp.fcost)
    disp('Warning: No constraints on the final condition have been specified. Since the');
    disp('         horizon is free and there is no final cost, this would cause the solver');
    disp('         to find a solution corresponding to a zero trajectory.');
    disp('         To avoid this situation the final condition is assigned automatically');
    disp('         (it is set to zero).');
    disp('         THE ASSIGNMENT IS ONLY TEMPORARILY.');
    pause(1);
    nstate=length(ocp.state);
    ocp=setfdirac(ocp, ocp.state, zeros(nstate,1));
end

if isempty(ocp.intcost) && isempty(ocp.fcost)
    disp('Warning: No cost has been specified. The cost is assigned automatically in order to');
    disp('         minimize the trace of the trajectory moment matrix.');
    disp('         THE ASSIGNMENT IS ONLY TEMPORARILY.');
    pause(1);
    %the assignment takes place when the function build is called.
end

%building problem
ocp = build(ocp, degree);

%solving problem
if nargout < 4
    [status,cost] = msol(ocp.gloptipoly);
    measures.initial = ocp.imeas;
    measures.final = ocp.fmeas;
    measures.trajectory = ocp.tmeas;
else
    mset('yalmip',true);
    [status,cost] = msol(ocp.gloptipoly);
    measures.initial = ocp.imeas;
    measures.final = ocp.fmeas;
    measures.trajectory = ocp.tmeas;
    coefficients=double(dual(MMM.yalmipdata.F(1)));
    varargout={ocp.testfun'*coefficients(1:length(ocp.testfun))};
end

if status==1
    disp(' ');
    disp('ATTENTION: global optimality of the solution to the moment problem does not');
    disp('     imply that the solution of the optimal control problem has been found.');
end