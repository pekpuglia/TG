function [ocp] = setstate(ocp, state)
%@POCP/PRIVATE/SETSTATE - Internal use only.

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

% Sets the state vector.

%checking input 
if ~isvector(state) || ~isa(state, 'mpol')
    error('the first argument must be a vector of type mpol');
end

[nrows, ncols]=size(state);

if ncols==0 || nrows==0
    error('the first argument must be a vector of type mpol');
end

for index = 1:length(state)
    if ~isvar(state(index))
        error('the first argument must be a vector of type mpol');
    end
end

if inuse(ocp, state, ocp.input, ocp.time)
    error('one or more of the variables entered are already in use in the current problem');
end

if ~isempty(ocp.state)
    error('the state vector cannot be changed');
end

%setting state (always saved as a column vector)
if size(state, 2)>1
    ocp.state=state';
else
    ocp.state=state;
end