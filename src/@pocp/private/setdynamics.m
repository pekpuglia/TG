function [ocp] = setdynamics(ocp, dynamics)
%@POCP/PRIVATE/SETDYNAMICS - Internal use only.

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

% Sets the system dynamics for the OCP.

%TODO checking if the variables is already used

if ~isvector(dynamics) || ~isa(dynamics, 'mpol')
    error('the first argument must be a column vector of type mpol');
end

[nrows, ncols]=size(dynamics);

if ncols~=1 || nrows==0
    error('the first argument must be a column vector of type mpol');
end

if ~valdynamics(ocp, dynamics)
    error('one or more variables used are not declared variables');
end

if length(dynamics)~=length(ocp.state)
    error('the dynamics specified has wrong dimension');
end


%setting dynamics
ocp.dynamics=dynamics;