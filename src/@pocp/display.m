function [] = display(ocp)
%POCP/DISPLAY - Display the polynomial optimal control problem data
%
%   DISPLAY(OCP) display the data of OCP

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

nstate = length(ocp.state);
ninput = length(ocp.input);

disp(' ');
disp('POCP - Polynomial Optimal Control Problem');

disp(' ');
disp('state variables');
if isempty(ocp.state)
    disp('  not defined');
else
    for index=1:nstate
        [tmp, state] = display(ocp.state(index));
        disp(['  ' state{1,1}])
    end
end

if ~isempty(ocp.input)
    disp(' ');
    disp('input variables');
    for index=1:ninput
        [tmp, input] = display(ocp.input(index));
        disp(['  ' input{1,1}])
    end
end

if ~isempty(ocp.time)
    disp(' ');
    disp('time variable');
    [tmp, time] = display(ocp.time);
    disp(['  ' time{1,1}])
end

disp(' ');
disp('dynamics');
if isempty(ocp.dynamics)
    disp('  not defined');
else
    if ocp.discrete
        disp('  discrete-time system');
    else
        disp('  continuous-time system');
    end
    if isempty(ocp.time)
        time = '(time)';
    else
        [tmp, time] = display(ocp.time);
        time = time{1,1};
    end
    for index=1:nstate
        [tmp, state] = display(ocp.state(index));
        [tmp, dynamics] = display(ocp.dynamics(index));
        disp(['  D(' state{1,1} ') = ' dynamics{1,1}])
    end
end

disp(' ');
disp('horizon');
if isempty(ocp.horizon)
    disp('  not defined')
elseif ocp.horizon==0
    disp('  free');
else
    disp(['  ' num2str(ocp.horizon)]);
end

if ~isempty(ocp.intcost)
    disp(' ');
    disp('integral cost');
    [tmp, intcost] = display(ocp.intcost);
    disp(['  ' intcost{1,1}])
end

if ~isempty(ocp.fcost)
    disp(' ');
    disp('final cost');
    [tmp, fcost] = display(ocp.fcost);
    disp(['  ' fcost{1,1}])
end

if ~isempty(ocp.tcon.int) || ~isempty(ocp.tcon.ineq)
    disp(' ');
    disp('-- constraints on the trajectory --');
end
    
if ~isempty(ocp.tcon.int)
    disp(' ');
    disp('integral constraints');
    for index = 1:length(ocp.tcon.int)
        [tmp, pol] = display(ocp.tcon.int);
        disp(['  ' pol])
    end
end

if ~isempty(ocp.tcon.ineq)
    disp(' ');
    disp('inequality constraints');
    for index = 1:length(ocp.tcon.ineq)
        [tmp, ineq] = display(ocp.tcon.ineq(index));
        disp(['  ' ineq{1,1}]);
    end
end

if ~isempty(ocp.icon.ineq) || ~isempty(ocp.icon.dirac.var) || ~isempty(ocp.icon.unif.var)
    disp(' ');
    disp('-- constraints at initial time --');
end

if ~isempty(ocp.icon.ineq)
    disp(' ');
    disp('inequality constraints');
    for index = 1:length(ocp.icon.ineq)
        [tmp, ineq] = display(ocp.icon.ineq(index));
        disp(['  ' ineq{1,1}]);
    end
end

if ~isempty(ocp.icon.dirac.var)
    disp(' ');
    disp('Dirac''s delta distributed variables');
    disp(' ');
    st = '          ';
    for vind = 1:length(ocp.icon.dirac.weight)
        st = [st sprintf('     point %3i', vind)];
    end
    disp(st);
    disp(' ');
    for vind = 1:length(ocp.icon.dirac.var)
        [tmp, variable] = display(ocp.icon.dirac.var(vind));
        st = sprintf('  %8s', variable{1,1});
        for pind = 1:size(ocp.icon.dirac.value, 2)
            st = [st sprintf('%14.4g', ocp.icon.dirac.value(vind, pind))];
        end
        disp(st);
    end
    disp(' ');
    st = '    weigth';
    for wind = 1:length(ocp.icon.dirac.weight)
        st = [st sprintf('%14.4g', ocp.icon.dirac.weight(wind))];
    end
    disp(st);
end

if ~isempty(ocp.icon.unif.var)
    disp(' ');
    disp('Uniformly distributed variables');
    disp(' ');
    for vind = 1:length(ocp.icon.unif.var)
        [tmp, variable] = display(ocp.icon.unif.var(vind));
        disp([sprintf('  %8s', variable{1,1}) '     [ ' num2str(ocp.icon.unif.interval(vind,1)) ' , ' num2str(ocp.icon.unif.interval(vind,2)) ' ]']);
    end
end

if ~isempty(ocp.fcon.ineq) || ~isempty(ocp.fcon.dirac.var)
    disp(' ');
    disp('-- constraints at final time --');
end

if ~isempty(ocp.fcon.ineq)
    disp(' ');
    disp('inequality constraints');
    for index = 1:length(ocp.fcon.ineq)
        [tmp, ineq] = display(ocp.fcon.ineq(index));
        disp(['  ' ineq{1,1}]);
    end
end

if ~isempty(ocp.fcon.dirac.var)
    disp(' ');
    disp('Dirac''s delta distributed variables');
    disp(' ');
    st = '          ';
    for vind = 1:length(ocp.fcon.dirac.weight)
        st = [st sprintf('     point %3i', vind)];
    end
    disp(st);
    disp(' ');
    for vind = 1:length(ocp.fcon.dirac.var)
        [tmp, variable] = display(ocp.fcon.dirac.var(vind));
        st = sprintf('  %8s', variable{1,1});
        for pind = 1:size(ocp.fcon.dirac.value, 2)
            st = [st sprintf('%14.4g', ocp.fcon.dirac.value(vind, pind))];
        end
        disp(st);
    end
    disp(' ');
    st = '    weigth';
    for wind = 1:length(ocp.fcon.dirac.weight)
        st = [st sprintf('%14.4g', ocp.fcon.dirac.weight(wind))];
    end
    disp(st);
end

disp(' ');