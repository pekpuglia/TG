function [out] = valdynamics(ocp, dynamics)
%@POCP/PRIVATE/VALDYNAMICS - Internal use only.

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

% Detects if the entered dynamics is valid

list=listvar(dynamics);

out = true;

for vind = 1:length(list)
    var = list(vind);
    evar = false;   %existing variable
    for eind = 1:length(ocp.state)
        if isequal(var, ocp.state(eind))
            evar = true;
        end
    end
    for eind = 1:length(ocp.input)
        if isequal(var, ocp.input(eind))
            evar = true;
        end
    end
    if isequal(var, ocp.time)
        evar = true;
    end
    
    if ~evar
        out = false;
        break
    end
end