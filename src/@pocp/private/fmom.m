function [fmon] = fmom(ocp, mon)
%@POCP/PRIVATE/FMON - Internal use only.

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

% Substitutes the variables in the vector of monomials mon.
% E.g. suppose mom is the monomial x(1)^2*x(2) where x is the state
% vector. Furthermore, assume that the value of x(1) at final time has
% been assigned, while x(2) has to be determined. In the returned fmon
% x(1)^2 will be substituted with the numeric value of the second order
% moment of the variable. x(2) will be substituted with the corresponding
% variable in the occupation measure at the end of the horizon.

fmon = mon;

for mind = 1:length(fmon)
    list = listvar(fmon(mind));
    for vind = 1:length(list)
        ex = deg(fmon(mind), list(vind));
        if isequal(list(vind), ocp.time)
            %time variable
            fmon(mind) = subs(fmon(mind), list(vind), ocp.ftime);
        elseif ~isempty(locate(ocp.fcon.dirac.var, list(vind)))
            %Dirac's delta
            cind = locate(ocp.fcon.dirac.var, list(vind));
            tmp = 0;
            for pind = 1:length(ocp.fcon.dirac.weight)
                tmp = tmp + ocp.fcon.dirac.weight(pind)*ocp.fcon.dirac.value(cind, pind)^ex;
            end
            fmon(mind) = subs(fmon(mind), list(vind), 1)*tmp;
        else
            %the statistics of the variable hasn't been assigned
            fmon(mind) = subs(fmon(mind), list(vind), ocp.fstate(locate(ocp.state, list(vind))));
        end
    end
end