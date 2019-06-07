function [f_mdot,s,inputs] = ideal_turbofan_ab(varargin)
% IDEAL_TURBOFAN_AB Ideal Turbofan with Afterburning
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
% [F_MDOT, S, INPUTS] = IDEAL_TURBOFAN_AB
% [F_MDOT, S, INPUTS] = IDEAL_TURBOFAN_AB(INPUTS)

disp(' ')
switch nargin
  case 0
    inputs.t0 = input('T0 (K): ');
    inputs.gam = input('Ratio of specific heats: ');
    inputs.h = input('Fuel heating value (J/kg): ');
    inputs.cp = input('Specific heat of air Cp (J/(kg*K)): ');
    inputs.tau_lam = input('tau_lam: ');

    disp(' ')
    disp('Afterburner Performance')
    disp('tau_lam_ab = 0  --> No primary stream burning')
    disp('tau_lam_ab2 = 0 --> No secondary stream burning')
    disp(' ')
    inputs.tau_lam_ab = input('tau_lam_ab: ');
    inputs.tau_lam_ab2 = input('tau_lam_ab2: ');
    disp(' ')

    inputs.pi_c = input('Compressor pressure ratio: ');
    inputs.pi_c2 = input('Fan pressure ratio: ');
    inputs.m0 = input('Flight Mach number: ');
    inputs.alpha = input('Bypass ratio: ');
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
tau_lam_ab2 = inputs.tau_lam_ab2;
pi_c = inputs.pi_c;
pi_c2 = inputs.pi_c2;
m0 = inputs.m0;
alpha = inputs.alpha;

% equations
r = ( (gam-1)./gam ).*cp;
a0 = sqrt(gam.*r.*t0);
tau_r = 1 + ((gam-1)./2).*m0.^2;
tau_c = pi_c.^((gam-1)./gam);
tau_c2 = pi_c2.^((gam-1)./gam);
tau_t = 1 - (tau_r./tau_lam).*( (tau_c - 1) + alpha.*(tau_c2 - 1) );

% check tau_lam_ab and tau_lam_ab2 inputs, or handle if zero
c1 = (tau_lam_ab > 0) & (tau_lam_ab < tau_lam.*tau_t);
c2 = (tau_lam_ab2 > 0) & (tau_lam_ab2 < tau_r.*tau_c2);
if any(c1) || any(c2)
  disp('tau_lam_ab must be >= tau_lam.*tau_t')
  disp('tau_lam_ab2 must be >= tau_r.*tau_c2')
  disp('Returning f_mdot = 0')
  disp('Returning s = 0')
  f_mdot = 0;
  s = 0;
  return
end

if any(tau_lam_ab == 0) && ~isscalar(tau_lam_ab)
  disp('tau_lam_ab = 0 must be a scalar')
  disp('tau_lam_ab2 = 0 must be a scalar')
  disp('Returning f_mdot = 0')
  disp('Returning s = 0')
  f_mdot = 0;
  s = 0;
  return
end

if any(tau_lam_ab2 == 0) && ~isscalar(tau_lam_ab2)
  disp('tau_lam_ab = 0 must be a scalar')
  disp('tau_lam_ab2 = 0 must be a scalar')
  disp('Returning f_mdot = 0')
  disp('Returning s = 0')
  f_mdot = 0;
  s = 0;
  return
end

if tau_lam_ab == 0
  tau_lam_ab = tau_lam.*tau_t;
end

if tau_lam_ab2 == 0
  tau_lam_ab2 = tau_r.*tau_c2;
end

% specific thrust and specific fuel consumption

ip1 = (2./(gam-1)) .* (tau_lam_ab./(tau_r.*tau_c.*tau_t)) .* ...
      (tau_r.*tau_c.*tau_t - 1);

ip2 = (2./(gam-1)) .* (tau_lam_ab2./(tau_r.*tau_c2)) .* ...
      (tau_r.*tau_c2 - 1);

f_mdot = (a0./(1+alpha)).*(ip1.^(1./2) - ...
			   m0 + alpha.*(ip2.^(1./2) - m0));

disp('f_mdot = F / (mdot_c + mdot_F)')
disp(' ')

f = (cp.*t0./h).*(tau_lam_ab + alpha.*tau_lam_ab2 - (1+alpha).*tau_r);
s = f.*1E6./( (1+alpha).*f_mdot );

end
