function [ocp] = sethorizon(ocp, horizon)
%POCP/SETHORIZON sets the OCP horizon
%   [ocp] = sethorizon(ocp, horizon) sets the horizon of the optimal
%   control problem to horizon. horizon should be a scalar grater than or equal
%   to 0 (0 means that the horizon is free)

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
msg='the second argument must scalar greater than or equal to zero';

if ~isnumeric(horizon) || ~(horizon>=0)
    error(msg);
end

%setting horizon
ocp.horizon=horizon;