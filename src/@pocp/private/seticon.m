function [ocp] = seticon(ocp, ineq)
%POCP/SETICON set the inequality constraints for the state variable at time 0
%   [ocp] = seticon(ocp, ineq) sets the inequality constraints for the
%      state at time 0 to ineq. ineq must be a vector of type supcon

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

%checking input
msg1 = 'the second argument must be a vector of type supcon';
msg2 = 'one or more variables used is not valid';

if ~((isvector(ineq) && isa(ineq, 'supcon')) || isempty(ineq))
    error(msg1);
end

if ~containsonly(ocp, [left(ineq) right(ineq)], ocp.state)
    error(msg2);
end

%setting constraint
ocp.icon.ineq = ineq;
