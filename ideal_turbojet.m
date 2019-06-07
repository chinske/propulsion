function [f_mdot,s,inputs] = ideal_turbojet(varargin)
% IDEAL_TURBOJET Ideal Cycle Analysis - Ideal Turbojet
% Copyright 2014 Christopher Chinske
% This file is part of GTPT.
% 
% GTPT is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% GTPT is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with GTPT.  If not, see <https://www.gnu.org/licenses/>.
% 
% [F_MDOT, S, INPUTS] = IDEAL_TURBOJET
% [F_MDOT, S, INPUTS] = IDEAL_TURBOJET(INPUTS)

switch nargin
  case 0
    disp(' ')
    inputs.t0 = input('T0 (K): ');
    inputs.gam = input('Ratio of specific heats: ');
    inputs.h = input('Fuel heating value (J/kg): ');
    inputs.cp = input('Specific heat of air Cp (J/(kg*K)): ');
    inputs.tau_lam = input('tau_lam: ');
    inputs.pi_c = input('Compressor pressure ratio: ');
    inputs.m0 = input('Flight Mach number: ');
    disp(' ')
  case 1
    inputs = varargin{1};
  otherwise
    error('Too many input arguments')
end

% reassign input variables
t0 = inputs.t0;
gam = inputs.gam;
h = inputs.h;
cp = inputs.cp;
tau_lam = inputs.tau_lam;
pi_c = inputs.pi_c;
m0 = inputs.m0;

% equations
r = ( (gam-1)./gam ).*cp;
a0 = sqrt(gam.*r.*t0);
tau_r = 1 + ((gam-1)./2).*m0.^2;
tau_c = pi_c.^((gam-1)./gam);
f_mdot = a0.*(( (2.*tau_r./(gam-1)) .* ...
		((tau_lam./(tau_r.*tau_c)) - 1) .* ...
		(tau_c - 1) + ...
		(tau_lam./(tau_r.*tau_c)).*m0.^2).^(1./2) - m0);
f = (cp.*t0./h).*(tau_lam - tau_r.*tau_c);
s = (f./f_mdot).*1E6;

end
