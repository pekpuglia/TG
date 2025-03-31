function [ocp] = setintcost(ocp, cost)
%@POCP/PRIVATE/SETINTCOST - Internal use only. 

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

% Sets the integral cost.

%checking input
if ~((isvector(cost) && isa(cost, 'mpol')) || isempty(cost) || isnumeric(cost))
    error('the first argument must be a scalar of type mpol, [] or a scalar numeric');
end

[nr, nc] = size(cost);

if nr>1 || nc>1
    error('the first argument must be a scalar of type mpol, [] or a scalar numeric');
end

cost = mpol(cost);

if ~containsonly(ocp, cost, ocp.state, ocp.input, ocp.time)
    error('one or more variables utilized is not valid');
end

%setting integral cost
ocp.intcost = cost;
