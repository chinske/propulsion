function [f_mdot,s,inputs] = ideal_turbojet_ab(varargin)
% IDEAL_TURBOJET_AB Ideal Turbojet with Afterburning
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
% [F_MDOT, S, INPUTS] = IDEAL_TURBOJET_AB
% [F_MDOT, S, INPUTS] = IDEAL_TURBOJET_AB(INPUTS)

switch nargin
  case 0
    disp(' ')
    inputs.t0 = input('T0 (K): ');
    inputs.gam = input('Ratio of specific heats: ');
    inputs.h = input('Fuel heating value (J/kg): ');
    inputs.cp = input('Specific heat of air Cp (J/(kg*K)): ');
    inputs.tau_lam = input('tau_lam: ');
    inputs.tau_lam_ab = input('tau_lam_ab: ');
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
tau_lam_ab = inputs.tau_lam_ab;
pi_c = inputs.pi_c;
m0 = inputs.m0;

% equations
r = ( (gam-1)./gam ).*cp;
a0 = sqrt(gam.*r.*t0);
tau_r = 1 + ((gam-1)./2).*m0.^2;
tau_c = pi_c.^((gam-1)./gam);
tau_t = 1 - (tau_r./tau_lam).*(tau_c - 1);

if sum(tau_lam_ab < (tau_lam.*tau_t)) ~= 0
  disp('tau_lam_ab must be >= tau_lam.*tau_t')
  disp('Returning f_mdot = 0')
  disp('Returning s = 0')
  f_mdot = 0;
  s = 0;
  return
end

f_mdot = a0.*(( (2./(gam-1)) .* ...
		(tau_lam_ab./(tau_r.*tau_c.*tau_t)) .* ...
		(tau_r.*tau_c.*tau_t - 1)).^(1./2)  - m0);

f = (cp.*t0./h).*(tau_lam_ab - tau_r);
s = (f./f_mdot).*1E6;

end
