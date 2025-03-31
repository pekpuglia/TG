function [ocp] = setiuniform(ocp, var, inter)
%@POCP/PRIVATE/SETIUNIFORM - Internal use only.

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

% Sets the statistics of the specified variables to be uniform at time 0.
%   [ocp] = setiunif(ocp, var, inter) sets the statistics of the variables
%      var to be uniform in the interval inter at time 0. var must be a
%      vector of variables of type mpol and inter must be a matrix with two
%      columns (coresponding to the beginning end the end of the interval)
%      and a number of rows equal to the length of var

%checking input
if ~((isvector(var) && isa(var, 'mpol')) || isempty(ineq))
    error('the first argument must be a vector of variables of type mpol');
end

nvar = length(var);
for index = 1:nvar
    if ~isvar(var(index))
        error('the first argument must be a vector of variables of type mpol');
    end
end

[nrows, ncols] = size(inter);

if ~(isnumeric(inter) && nrows==nvar && ((ncols==0 && nvar==0) || ncols==2)) 
    error('the second argument must be a numeric matrix with two columns and a number of rows equal to the length of the first argument');
end

if ~containsonly(ocp, var, ocp.state)
    error('one or more variables used is not valid');
end

if nvar>0 && min(inter(:,2)-inter(:,1))<0
    error('the second argument is not a valid interval');
end
    
%setting constraint
ocp.icon.unif.var = var;
ocp.icon.unif.interval = inter;
    