function [ocp] = build(ocp, tfdeg)
%@POCP/PRIVATE/BUILD - Internal use only

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

% Builds the problem to be solved by gloptipoly.
% tfdeg is the maximal degree of the test function used.

nvar = length(ocp.state); %number of state variables

%determining if the test functions should depend on time or not
testtime = false;

if ocp.testtime==true
    testtime = true;
end

if ocp.horizon~=0
    testtime = true;
    ocp.ftime = ocp.horizon;
end

if ocp.testtime==false
    testtime = false;
end

%determining if the statistics of the state vector has been completely
%assigned at time zero and at the end of the horizon
if (length(ocp.icon.dirac.var)+length(ocp.icon.unif.var)) < nvar
    ocp.iassigned = false;
else
    ocp.iassigned = true;
end

if (length(ocp.fcon.dirac.var)) < nvar
    ocp.fassigned = false;
else
    ocp.fassigned = true;
end

%creating needed variables to define the occupation measures
if ~ocp.iassigned
    eval(sprintf('mpol pocpinitialstate %d', nvar));
    ocp.istate = pocpinitialstate;
end

if ~ocp.fassigned
    eval(sprintf('mpol pocpfinalstate %d', nvar));
    ocp.fstate = pocpfinalstate;
end

if testtime && isempty(ocp.time)
    mpol pocptime;
    ocp.time = pocptime;
end

if testtime && ocp.horizon==0
    ocp.ftime = ocp.time;
end

%defining measures
if ~ocp.iassigned
    ocp.imeas = meas(ocp.istate);
end

if ~ocp.fassigned
    ocp.fmeas = meas(ocp.fstate);
end

if ~testtime
    ocp.tmeas = meas([ocp.state; ocp.input]);
else
    ocp.tmeas = meas([ocp.state; ocp.input; ocp.time]);
end

%setting the cost
cost = 0;
if ~isempty(ocp.intcost)
    if pow(ocp.intcost)==0
        %the integral cost is a scalar
        cost = cost + mass(ocp.tmeas)*double(ocp.intcost);
    else
        cost = cost + mom(ocp.intcost);
    end
end    
if ~isempty(ocp.fcost)
    cost = cost + mom(subs(ocp.fcost, ocp.state, ocp.fstate));
end
if isempty(ocp.intcost) && isempty(ocp.fcost)
    %cost is assigned automatically
    tmp = ocp.dynamics;
    maxdeg = ceil((tfdeg + deg(tmp) - 1)/2);
    if testtime
        vars = [ocp.state; ocp.input; ocp.time];
    else
        vars = [ocp.state; ocp.input];
    end
    cost = mom(sum(mmon(vars, 0, maxdeg).^2));
end

%generating test functions
if ~testtime
    if isempty(ocp.imeas) && isempty(ocp.fmeas)
        tfun = mmon(ocp.state, 1, tfdeg);
    else
        tfun = mmon(ocp.state, 0, tfdeg);
    end
else
    if isempty(ocp.imeas) && isempty(ocp.fmeas)
        tfun = mmon([ocp.state; ocp.time], 1, tfdeg);
    else
        tfun = mmon([ocp.state; ocp.time], 0, tfdeg);
    end
end
ocp.testfun = tfun;
ocp.ntestfun = length(tfun);

%moment constraints coming from the trajectory
if ~ocp.discrete 
    if ~testtime
        ocp.momcon = (0 == mom(fmom(ocp, tfun)) - mom(imom(ocp, tfun)) - mom(diff(tfun, ocp.state)*ocp.dynamics));
    else
        ocp.momcon = (0 == mom(fmom(ocp, tfun)) - mom(imom(ocp, tfun)) - mom(diff(tfun, ocp.time) + diff(tfun, ocp.state)*ocp.dynamics));
    end
else
    if ~testtime
        ocp.momcon = (0 == mom(fmom(ocp, tfun)) - mom(imom(ocp, tfun)) - mom(subs(tfun, ocp.state, ocp.dynamics)-tfun));
    else
        ocp.momcon = (0 == mom(fmom(ocp, tfun)) - mom(imom(ocp, tfun)) - mom(subs(tfun, [ocp.state; ocp.time], [ocp.dynamics; ocp.time+1])-tfun));
    end
end

%generating the moment constraints on the initial and final moments deriving from
%the fact that the statistics of (part) of the variables has been assigned
mons = mmon(ocp.state, 1, tfdeg);
%initial state
if ~isempty(ocp.imeas)
    left = subs(mons, ocp.state, ocp.istate);
    right = imom(ocp, mons);
    %eliminating useless elements
    for index = length(left):-1:1
        if isequal(left(index), right(index))
            left = left([1:index-1 index+1:end]);
            right = right([1:index-1 index+1:end]);
        end
    end
    if ~isempty(left)
        ocp.momcon = [ocp.momcon; 0 == mom(left) - mom(right)];
    end
    ocp.momcon = [ocp.momcon; mass(ocp.imeas)==1];
end

%final state
if ~isempty(ocp.fmeas)
    left = subs(mons, ocp.state, ocp.fstate);
    right = fmom(ocp, mons);
    %eliminating useless elements
    for index = length(left):-1:1
        if isequal(left(index), right(index))
            left = left([1:index-1 index+1:end]);
            right = right([1:index-1 index+1:end]);
        end
    end
    if ~isempty(left)
        ocp.momcon = [ocp.momcon; 0 == mom(left) - mom(right)];
    end
    ocp.momcon = [ocp.momcon; mass(ocp.fmeas)==1];
end

%adding moment constraints coming from the integral constraints on the
%trajectory
ocp.momcon = [ocp.momcon; ocp.tcon.int];

%adding support constraints
ocp.supcon = ocp.tcon.ineq;

if ~isempty(ocp.imeas)
    ocp.supcon = [ocp.supcon; subs(ocp.icon.ineq, ocp.state, ocp.istate)];
end

if ~isempty(ocp.fmeas)
    ocp.supcon = [ocp.supcon; subs(ocp.fcon.ineq, ocp.state, ocp.fstate)];
end

%buiding gloptipoly problem
if isempty(ocp.supcon)
    ocp.gloptipoly = msdp(min(cost), ocp.momcon);
else
    ocp.gloptipoly = msdp(min(cost), ocp.momcon, ocp.supcon);
end
