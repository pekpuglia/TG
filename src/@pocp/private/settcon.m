function [ocp] = settcon(ocp, ineq)
%@POCP/PRIVATE/SETTCON - Internal use only.

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

% Sets the inequality constraints for the trajectory.

%checking input
if ~((isvector(ineq) && isa(ineq, 'supcon')) || isempty(ineq))
    error('the first argument must be a vector of type supcon');
end

if ~containsonly(ocp, [left(ineq) right(ineq)], ocp.state, ocp.input, ocp.time)
    error('one or more variables used are not valid');
end

%setting constraint
ocp.tcon.ineq = ineq;