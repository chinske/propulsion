function [tau_t_opt,alpha_opt,f_mdot_opt,s_opt,inputs] = ...
	 ideal_turbofan_opt(varargin)
% IDEAL_TURBOFAN Ideal Cycle Analysis - Ideal Turbofan Optimal
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
% [TAU_T_OPT, ALPHA_OPT, F_MDOT_OPT, S_OPT, INPUTS] = ...
% IDEAL_TURBOFAN_OPT
% 
% [TAU_T_OPT, ALPHA_OPT, F_MDOT_OPT, S_OPT, INPUTS] = ...
% IDEAL_TURBOFAN_OPT(INPUTS)

disp(' ')
switch nargin
  case 0
    disp('Constants')
    disp('---------')
    inputs.gam = input('Ratio of specific heats: ');
    inputs.h = input('Fuel heating value (J/kg): ');
    inputs.cp = input('Specific heat of air Cp (J/(kg*K)): ');
    disp(' ')
    disp('Operating Conditions')
    disp('--------------------')
    inputs.t0 = input('T0 (K): ');
    inputs.m0 = input('M0: ');
    disp(' ')
    disp('Design Limits')
    disp('-------------')
    inputs.tau_lam = input('tau_lam: ');
    disp(' ')
    disp('Design Choices')
    disp('--------------')
    inputs.pi_c = input('Compressor pressure ratio: ');
    inputs.pi_c2 = input('Fan pressure ratio: ');
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

% --------------------------------------------------
% equations
r = ( (gam-1)./gam ).*cp;
a0 = sqrt(gam.*r.*t0);
tau_r = 1 + ((gam-1)./2).*m0.^2;
tau_c = pi_c.^((gam-1)./gam);
tau_c2 = pi_c2.^((gam-1)./gam);

% --------------------------------------------------
% turbine temperature ratio and bypass ratio
% optimal case

tau_t = (1./(tau_r.*tau_c)) + ...
	(1./(4.*tau_lam)).*((tau_r.*tau_c2 - 1).^(1./2) + ...
			    (tau_r - 1).^(1./2)).^2;

alpha = (1./(tau_c2 - 1)).*((tau_lam./tau_r).*(1 - tau_t) - ...
			    (tau_c - 1));

% --------------------------------------------------
% equations cont'd with optimal values

ip1 = (2./(gam-1)) .* (tau_lam./(tau_r.*tau_c)) .* ...
      (tau_r.*tau_c.*tau_t - 1);

ip2 = (2./(gam-1)) .* (tau_r.*tau_c2 - 1);

f_mdot = (a0./(1+alpha)).*(ip1.^(1./2) - ...
			   m0 + alpha.*(ip2.^(1./2) - m0));

disp('f_mdot = F / (mdot_c + mdot_F)')
disp(' ')

f = (cp.*t0./h).*(tau_lam - tau_r.*tau_c);
s = f.*1E6./( (1+alpha).*f_mdot );

% --------------------------------------------------
% reassign variables to output
tau_t_opt = tau_t;
alpha_opt = alpha;
f_mdot_opt = f_mdot;
s_opt = s;

end
