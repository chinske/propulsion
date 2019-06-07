function [f_mdot,s,inputs] = ideal_turbofan_mixed(varargin)
% IDEAL_TURBOFAN_MIXED Ideal Turbofan with Mixed Exhaust Streams
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
% [F_MDOT, S, INPUTS] = IDEAL_TURBOFAN_MIXED
% [F_MDOT, S, INPUTS] = IDEAL_TURBOFAN_MIXED(INPUTS)

disp(' ')
switch nargin
  case 0
    inputs.t0 = input('T0 (K): ');
    inputs.gam = input('Ratio of specific heats: ');
    inputs.h = input('Fuel heating value (J/kg): ');
    inputs.cp = input('Specific heat of air Cp (J/(kg*K)): ');
    inputs.tau_lam = input('tau_lam: ');
    inputs.pi_c = input('Compressor pressure ratio: ');
    inputs.pi_c2 = input('Fan pressure ratio: ');
    inputs.m0 = input('Flight Mach number: ');
    inputs.alpha = input('Bypass ratio: ');
    inputs.m5 = input('Mach number at turbine exit: ');
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
pi_c2 = inputs.pi_c2;
m0 = inputs.m0;
alpha = inputs.alpha;
m5 = inputs.m5;

% equations
r = ( (gam-1)./gam ).*cp;
a0 = sqrt(gam.*r.*t0);
tau_r = 1 + ((gam-1)./2).*m0.^2;
tau_c = pi_c.^((gam-1)./gam);
tau_c2 = pi_c2.^((gam-1)./gam);
tau_t = 1 - (tau_r./tau_lam).*( (tau_c - 1) + alpha.*(tau_c2 - 1) );
pi_t = tau_t.^(gam./(gam-1));

% equations for the ideal constant-area mixer
pt32pt5 = pi_c2./(pi_c.*pi_t);
tt32tt5 = tau_r.*tau_c2./(tau_lam.*tau_t);

[pt7pt5, m7, m32] = ideal_ca_mixer(pt32pt5,tt32tt5,m5,alpha,gam);

pt9pt5 = pt7pt5;
pt9pt32 = pt9pt5.*pi_c.*pi_t./pi_c2;

% specific thrust and specific fuel consumption

term1 = (2./(gam-1))./(1+alpha);
term2 = 1 - 1./(tau_r.*tau_c2.*pt9pt32.^((gam-1)./gam));
term3 = tau_lam.*tau_t + alpha.*tau_r.*tau_c2;

f_mdot = a0.*( (term1.*term2.*term3).^(1./2) - m0 );

disp('f_mdot = F / (mdot_c + mdot_F)')
disp(' ')

f = (cp.*t0./h).*(tau_lam - tau_r.*tau_c);
s = f.*1E6./( (1+alpha).*f_mdot );

end
