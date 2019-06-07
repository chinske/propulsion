function [f_mdot,s,inputs] = nonideal_turbojet(varargin)
% NONIDEAL_TURBOJET Turbojet with Losses
% Copyright 2015 Christopher Chinske
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
% [F_MDOT, S, INPUTS] = NONIDEAL_TURBOJET
% [F_MDOT, S, INPUTS] = NONIDEAL_TURBOJET(INPUTS)

switch nargin
  case 0
    disp(' ')

    inputs.ab = input('Afterburning flag: ');
    if inputs.ab == 0
      disp('Afterburning is not present')
    elseif inputs.ab == 1
      disp('Afterburning is present')
    else
      disp('Afterburning flag must be 0 or 1')
      disp('Returning f_mdot = 0')
      disp('Returning s = 0')
      f_mdot = 0;
      s = 0;
      return
    end

    inputs.t0 = input('T0 (K): ');

    inputs.gam_c = input('Ratio of specific heats, compressor: ');
    inputs.gam_t = input('Ratio of specific heats, turbine: ');
    if inputs.ab == 1
      inputs.gam_ab = input('Ratio of specific heats, afterburner: ');
    end

    inputs.cp_c = ...
    input('Specific heat of air Cp (J/(kg*K)), compressor: ');
    inputs.cp_t = ...
    input('Specific heat of air Cp (J/(kg*K)), turbine: ');
    if inputs.ab == 1
      inputs.cp_ab = ...
      input('Specific heat of air Cp (J/(kg*K)), afterburner: ');
    end

    inputs.h = input('Fuel heating value (J/kg): ');

    inputs.pi_d = input('Pressure ratio, inlet: ');
    inputs.pi_b = input('Pressure ratio, burner: ');
    inputs.pi_n = input('Pressure ratio, nozzle: ');

    inputs.eta_b = input('Efficiency, burner: ');
    if inputs.ab == 1
      inputs.eta_ab = input('Efficiency, afterburner: ');
    end
    inputs.eta_m = input('Efficiency, mechanical: ');

    inputs.e_c = input('Polytropic efficiency, compressor: ');
    inputs.e_t = input('Polytropic efficiency, turbine: ');

    inputs.p9p0 = input('p9/p0: ');

    inputs.tau_lam = input('tau_lam: ');
    if inputs.ab == 1
      inputs.tau_lam_ab = input('tau_lam_ab: ');
    end

    inputs.pi_c = input('Compressor pressure ratio: ');
    inputs.m0 = input('Flight Mach number: ');

    disp(' ')
  case 1
    inputs = varargin{1};
  otherwise
    error('Too many input arguments')
end

% reassign input variables
if inputs.ab == 0
  disp('Afterburning is not present')
  ab = inputs.ab;
elseif inputs.ab == 1
  disp('Afterburning is present')
  ab = inputs.ab;
else
  disp('Afterburning flag must be 0 or 1')
  disp('Returning f_mdot = 0')
  disp('Returning s = 0')
  f_mdot = 0;
  s = 0;
  return
end

t0 = inputs.t0;
gam_c = inputs.gam_c;
gam_t = inputs.gam_t;
cp_c = inputs.cp_c;
cp_t = inputs.cp_t;
h = inputs.h;
pi_d = inputs.pi_d;
pi_b = inputs.pi_b;
pi_n = inputs.pi_n;
eta_b = inputs.eta_b;
eta_m = inputs.eta_m;
e_c = inputs.e_c;
e_t = inputs.e_t;
p9p0 = inputs.p9p0;
tau_lam = inputs.tau_lam;
pi_c = inputs.pi_c;
m0 = inputs.m0;

if ab == 1
  gam_ab = inputs.gam_ab;
  cp_ab = inputs.cp_ab;
  eta_ab = inputs.eta_ab;
  tau_lam_ab = inputs.tau_lam_ab;
else
  cp_ab = cp_t;
  gam_ab = gam_t;
  eta_ab = 0;
end

% equations
r_c = ((gam_c-1)./gam_c).*cp_c;
a0 = sqrt(gam_c.*r_c.*t0);
tau_r = 1 + ((gam_c-1)./2).*m0.^2;
pi_r = tau_r.^(gam_c./(gam_c-1));
tau_c = pi_c.^((gam_c-1)./(gam_c.*e_c));
f = (tau_lam - tau_r.*tau_c)./((h.*eta_b./(cp_c.*t0)) - tau_lam);
tau_t = 1 - (1./(eta_m.*(1+f))).*(tau_r./tau_lam).*(tau_c-1);
pi_t = tau_t.^(gam_t./((gam_t-1).*e_t));
pt9p9 = (1./p9p0).*pi_r.*pi_d.*pi_c.*pi_b.*pi_t.*pi_n;

if ab == 0
  tau_lam_ab = tau_lam.*tau_t;
end

t9t0 = ((cp_c./cp_ab).*tau_lam_ab)./(pt9p9.^((gam_ab-1)./gam_ab));

f_ab = (1+f).*(tau_lam_ab - tau_lam.*tau_t)./ ...
       ((h.*eta_ab./(cp_c.*t0)) - tau_lam_ab);

term1 = 2./(gam_c-1);
term2 = tau_lam_ab;
term3 = 1 - pt9p9.^(-((gam_ab-1)./gam_ab));

m0u9u0 = (term1.*term2.*term3).^(1./2);

term4 = 1 + f + f_ab;

f_mdot = a0.*(term4.*m0u9u0 - m0 + ...
	      (term4.*t9t0./(gam_c.*m0u9u0)).*(1-(1./p9p0)));

s = ((f+f_ab)./(f_mdot)).*10.^6;

end
