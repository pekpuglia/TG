function [ocp] = set(ocp, varargin)
%@POCP/SET - Set polynomial optimal control problem data
%
%   [PROB]=SET(PROB, 'ParameterName1', 'Argument1', 'Argument2', ...,
%      'ParamenterName2', ...) sets the data of the polynomial optimal
%      control problem PROB.
%
%   The value of 'ParameterName' can be chosen between the following
%   values (partial matching is performed):
%   GENERAL DATA
%   -'state'. Mandatory field. One argument specifying the vector of state
%    variables must be entered. The state variables must be of type mpol.
%      Example
%       mpol x 2;
%       prob = pocp();
%       prob = set(prob, 'state', x);
%   -'input'. One argument specifying the vector of input variables must be
%    entered. The input variables must be of type mpol.
%      Example (continued)
%       mpol u;
%       prob = set(prob, 'input', u);
%   -'time'. One argument specifying the time variable must be
%    entered.
%      Example (continued)
%       mpol t;
%       prob = set(prob, 'time', t);
%   -'dynamics'. Mandatory field. One argument specifying the system dynamics
%    must be entered. The dynamics must be a vector of type mpol of the same
%    size of the specified state vector.
%      Example (continued)
%       A = [0 1; 1 1];
%       B = [0; 1];
%       prob = set(prob, 'dynamics', A*x+B*u);
%   -'horizon'. One argument specifying the horizon of the optimal control
%    problem must be entered. If set to 0 the horizon is free. The default
%    value of this field is 0.
%      Example (continued)
%       prob = set(prob, 'horizon', 1);
%   -'fcost'. One argument specifying the final
%    cost (or terminal cost) must be entered.
%      Example (continued)
%       prob = set(prob, 'fcost', x(1)^2);
%   -'scost'. One argument specifying the
%    integral cost (or running cost) must be specified.
%      Example (continued)
%       prob = set(prob, 'scost', x'*x);
%   -'testtime'. One argument must be specified. If the argument true (false),
%    test functions will (not) depend on time. If [] is specified the
%    default behaviour is set (the test functions don't depend on time if
%    the dynamics doesn't depend on time and the horizon is free).
%    When the default behaviour is such that the test functions don't
%    depend on time, setting 'testtime' to true increases CPU time.
%      Example (continued)
%       prob = set(prob, 'testtime', false);
%   CONSTRAINTS ON THE INITIAL STATE
%   -'iconstraint'. Used to set polynomial constraints on the value of the state
%    variables at time 0. One argument specifying the constraints must be
%    entered. The argument must be a scalar (or a vector) of type supcon.
%      Example (continued)
%       prob = set(prob, iconstraint, x(1)<=1);
%   -'idirac'. Used when one or more state variables at time 0 take their value
%    inside a discrete set of points with a specified probability. Two or three
%    arguments must be specified. The first argument is a vector containing the
%    state variables involved. The second argument is a matrix whose columns
%    contain the value of the variables. The third argument (optional when the
%    set contains only one point) is a row vector specifying the probability of
%    each point).
%      Example (continued)
%       % the value of the state x at time 0 is x=(0, 1) with probability 0.8
%       % and x=(1, 1) with probability 0.2.
%       prob = set(prob, 'idirac', x, [0 1; 1 1]; [0.8 0.2]);
%   -'iuniform'. Used when one or more state variables at time 0 take their
%    value inside an interval with uniform probability. Two argument must be
%    specified. The first argument is a vector containing the state variables
%    involved. The second argument is a matrix with two columns and whose rows
%    contain the two values which specify the beginning and end of the interval.
%      Example (continued)
%       % the variable x(1) is uniformly distriburted in [-1, 1]
%       prob = set(prob, 'iuniform', x(1), [-1, 1]);
%   CONSTRAINTS ON THE FINAL STATE
%   -'fconstraint'. Used to set polynomial constraints on the value of the state
%    variables at the end of the horizon. Same usage as 'iconstraint'.
%   -'fdirac' (long name 'finaldirac'): used when one or more state
%    variables at the end of the horizon take their value inside a discrete
%    set of points with a specified probability. Same usage as 'idirac'.
%   CONSTRAINTS ON THE VARIABLES ALONG THE TRAJECTORY
%   -'tconstraint'. One argument specifying the polynomial constraints on the
%    variables must be entered.
%      Example (continued)
%       prob = set(prob, 'tconstraint', [u<=1; x(1)+t+u<=1]);
%   -'sconstraint'. Used when the integral of one or more polynomial function
%    of the system variables is constrained. One argument containing the
%    constraints of type momcon must be specified.
%      Example (continued)
%       % suppose the integral for the time t going from 0 to the end of the 
%       % horizon of u(t)^2 must be smaller than 1
%       prob = set(prob, mom(u^2)<=1);
% 
%   NOTE 1: the current version of POCP doesn't check for conflicts between
%    different kinds of constraints.
%
%   NOTE 2: the system variables (state vector, input vector, time) must be
%    set before they are used. E.g., constraints on the initial state cannot
%    be enetered if the state vector has not been set.
%
%   Full example
%    mpol x;
%    mpol u;
%    prob = pocp();
%    prob = set(prob, ...
%           'state', x, ...
%           'input', u, ...
%           'dynamics', x+u, ...
%           'horizon', 0, ... % free horizon
%           'idirac', x, 1, ... % the initial condition is x(0)=1
%           'fdirac', x, 0, ... % the state should reach the origin
%           'tcon', u^2<=1, ... % constraint on the input varaible
%           'scost', 1); % the running cost is set to 1 to obtain a
%                          % minimum time optimal control problem

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

% The ParameterName 'discrete' is not documented at the moment

% list containing:
% function name which can be found in private, parameter name, minimum number
% of characters for partial matching, minimum number of arguments, maximal number of arguments
list = { ...
    'setstate',     'state',            1, 1, 1;
    'setinput',     'input',            2, 1, 1;
    'settime',      'time',             2, 1, 1;
    'setdynamics',  'dynamics',         2, 1, 1;
    'sethorizon',   'horizon',          1, 1, 1;
    'setfcon',      'fconstraint',      4, 1, 1;
    'setfcost',     'fcost',            4, 1, 1;
    'setfdirac',    'fdirac',           2, 2, 3;
    'seticon',      'iconstraint',      2, 1, 1;
    'setidirac',    'idirac',           2, 2, 3;
    'setintcon',    'sconstraint',      4, 1, 1;
    'setintcost',   'scost',            4, 1, 1;
    'setiuniform',  'iuniform',         2, 2, 2;
    'settcon',      'tconstraint',      1, 1, 1;
    'settesttime',  'testtime',         2, 1, 1;
    'setdiscrete',  'discrete',         2, 0, 1; ...
    };
    
narg = nargin-1;

index = 1;

while index<=narg
    argument = varargin{index};
    index = index+1;
    %looking for argument in list
    idx = match(list, argument);
    if idx==0
        error('invalid option "%s"', argument);
    end
    %building the string containing the command to be executed
    cmd = sprintf('ocp = %s(ocp', list{idx, 1});
    maxarg = list{idx, 5};
    minarg = list{idx, 4};
    finished = false;
    for aind = 1:min(maxarg, narg-index+1)
        if aind<=minarg
            cmd = sprintf('%s, varargin{%d}', cmd, index);
            index = index+1;
        else
            %if next varargin is a valid name of property to be set pwd then the parsing
            %of the current argument is terminated
            if ischar(varargin{index})
                if ~match(list, varargin{index})
                    cmd = sprintf('%s, varargin{%d}', cmd, index);
                    index = index+1;
                end
            else
                cmd = sprintf('%s, varargin{%d}', cmd, index);
                index = index+1;
            end
        end
        if finished
            break
        end
    end
    cmd = sprintf('%s);', cmd);
    eval(cmd);
end

% function which checks a possible matching of a parameter name inside the
% talbe list
function [idx] = match(list, string)
idx = 0;
if ischar(string)
    for sind = 1:size(list, 1)
        if length(string)>=list{sind, 3}
            if strncmp(string, list{sind, 2}, length(string))
                idx = sind;
            end
        end
    end
end
