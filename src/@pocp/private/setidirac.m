function [ocp] = setidirac(ocp, var, value, varargin)
%POCP/SETIDIRAC set the statistics of the specified variables to be a
%   combination of Dirac's deltas at time 0
%   [ocp] = setiunif(ocp, var, value, weight) sets the statistics of the
%       variables var to be a combination of diracs. Every column of value
%       must specify the support of a dirac's delta. the Dirac's deltas are
%       weighted according to weight. If weight is not specified, equal
%       weights are used

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
msg1 = 'the second argument must be a vector of variables of type mpol';
msg2 = 'the third argument must be a numeric matrix with an arbitrary number of columns and a number of rows equal to the length of the second argument';
msg3 = 'one or more variables used is not valid';
msg4 = 'the forth argument must be a vector of positive numbers whose length is equal to the number of columns of the third argument';

if ~((isvector(var) && isa(var, 'mpol')) || isempty(ineq))
    error(msg1);
end

nvar = length(var);
for index = 1:nvar
    if ~isvar(var(index))
        error(msg1);
    end
end

[nrows, ncols] = size(value);

if ~(isnumeric(value) && nrows==nvar && (ncols>0 || (ncols==0 && nvar==0))) 
    error(msg2);
end

if ~containsonly(ocp, var, ocp.state)
    error(msg3);
end

if nargin==4
    weight = varargin{1};
    if ~(isnumeric(weight) && isvector(weight) && min(weight)>=0 && sum(weight)~=0)
        error(msg4);
    end
    if length(weight)~=ncols
        error(msg4);
    end
    weight = weight/sum(weight);
else
    weight = ones(ncols,1)/ncols;
end

%setting constraint
ocp.icon.dirac.var = var;
ocp.icon.dirac.value = value;
ocp.icon.dirac.weight = weight;
    