function [out] = inuse(ocp, tocheck, varargin)
%@POCP/PRIVATE/INUSE - Internal use only.

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

% Detects if variables in tocheck are already used in varargin.

if nargin<3
    error('not enough input arguments');
end

nexisting=nargin-2; %number of vectors of existing variables

for eind = 1:nexisting
    existing=varargin{eind};
    for evind = 1:length(existing)
        for nind = 1:length(tocheck)
            if isequal(tocheck(nind), existing(evind))
                out = true;
                return
            end
        end
    end
end

out = false;