function [out] = containsonly(ocp, tocheck, varargin)
%@POCP/PRIVATE/CONTAINSONLY - Internal use only.

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

% Detects if the polynomial vector tocheck contains only variables
% specified in varargin


if nargin<3
    error('not enough input arguments');
end

tocheck=listvar(tocheck);
nexisting=nargin-2; %number of vectors of existing variables
out = true;

for nind = 1:length(tocheck)
    found = false;
    for eind = 1:nexisting
        existing=varargin{eind};
        for evind = 1:length(existing)
            if isequal(tocheck(nind), existing(evind))
                found = true;
                break
            end
        end
    end
    if found==false
        out = false;
        return
    end
end