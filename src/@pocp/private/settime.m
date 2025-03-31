function [ocp] = settime(ocp, time)
%@POCP/PRIVATE/SETTIME - Internal use only.

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

% Sets the time variable.

%TODO checking if the variables is already used

%checking input
if ~((isa(time, 'mpol') && isvar(time)) || isempty(time))
    error('the first argument must be a variable of type mpol');
end

if ocp.timeset
    error('the time variable cannot be changed');
end

%setting time
ocp.time = time;
ocp.timeset = true;