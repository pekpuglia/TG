function [ocp] = setinput(ocp, input)
%@POCP/PRIVATE/SETINPUT - Internal use only. 

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

% Sets the input variables.

%checking input
if ~((isvector(input) && isa(input, 'mpol')) || isempty(input))
    error('the first argument must be a vector of type mpol');
end

for index = 1:length(input)
    if ~isvar(input(index))
        error('the first argument must be a vector of type mpol');
    end
end

if inuse(ocp, input, ocp.state, ocp.time)
    error('one or more of the variables entered are already in use in the current problem');
end

if ~isempty(ocp.input)
    error('the input vector cannot be changed');
end
    
%setting input (always saved as a column vector)
if size(input, 2)>1
    ocp.input = input';
else
    ocp.input = input;
end
