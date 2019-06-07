function [f_mdot,s,inputs] = nonideal_turbofan2(inputs,verbose)
% NONIDEAL_TURBOFAN2 Turbofan with Losses, Convergent Exit Nozzles
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
% [F_MDOT, S, INPUTS] = NONIDEAL_TURBOFAN2(INPUTS,VERBOSE)
% 
% VERBOSE == 1 : Displays additional information that is required for
% off-design performance computations.
% 
% Application
% -----------
% Subsonic transport
% Separate stream turbofan engines
% No afterburning
% Pressure ratio across nozzles is small
% 
% Only convergent nozzles are employed
% 

% set afterburning flags to 0
disp('Setting afterburning flags to 0')
inputs.ab = 0;
inputs.ab2 = 0;

% reassign input variables
if inputs.ab == 0
  disp('Primary stream afterburning is not present')
  ab = inputs.ab;
elseif inputs.ab == 1
  disp('Primary stream afterburning is present')
  ab = inputs.ab;
else
  disp('Afterburning flag must be 0 or 1')
  disp('Returning f_mdot = 0')
  disp('Returning s = 0')
  f_mdot = 0;
  s = 0;
  return
end

if inputs.ab2 == 0
  disp('Secondary stream afterburning is not present')
  ab2 = inputs.ab2;
elseif inputs.ab2 == 1
  disp('Secondary stream afterburning is present')
  ab2 = inputs.ab2;
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
pi_n2 = inputs.pi_n2;
eta_b = inputs.eta_b;
eta_m = inputs.eta_m;
e_c = inputs.e_c;
e_c2 = inputs.e_c2;
e_t = inputs.e_t;
% p9p0 = inputs.p9p0;
% p92p0 = inputs.p92p0;
tau_lam = inputs.tau_lam;
pi_c = inputs.pi_c;
pi_c2 = inputs.pi_c2;
m0 = inputs.m0;
alpha = inputs.alpha;

% check Mach number
if m0 >= 1
  disp('Mach number must be less than 1')
  disp('Returning f_mdot = 0')
  disp('Returning s = 0')
  f_mdot = 0;
  s = 0;
  return
end

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

if ab2 == 1
  gam_ab2 = inputs.gam_ab2;
  cp_ab2 = inputs.cp_ab2;
  eta_ab2 = inputs.eta_ab2;
  tau_lam_ab2 = inputs.tau_lam_ab2;
else
  cp_ab2 = cp_c;
  gam_ab2 = gam_c;
  eta_ab2 = 0;
end

% equations
r_c = ((gam_c-1)./gam_c).*cp_c;
a0 = sqrt(gam_c.*r_c.*t0);
tau_r = 1 + ((gam_c-1)./2).*m0.^2;
pi_r = tau_r.^(gam_c./(gam_c-1));
tau_c = pi_c.^((gam_c-1)./(gam_c.*e_c));
tau_c2 = (pi_c2).^((gam_c-1)./(gam_c.*e_c2));
f = (tau_lam - tau_r.*tau_c)./((h.*eta_b./(cp_c.*t0)) - tau_lam);
tau_t = 1 - (1./(eta_m.*(1+f))).*(tau_r./tau_lam).* ...
	    ((tau_c-1) + alpha.*(tau_c2-1));
pi_t = tau_t.^(gam_t./((gam_t-1).*e_t));

% operation with convergent exit nozzles
p0p92 = (((gam_c+1)./2).^(gam_c./(gam_c-1)))./ ...
	(pi_r.*pi_d.*pi_c2.*pi_n2);

if p0p92 > 1
  disp('Secondary nozzle not choked')
  p92p0 = 1;
  p0p92 = 1;
else
  p92p0 = 1./p0p92;
end

p0p9 = (((gam_t+1)./2).^(gam_t./(gam_t-1)))./ ...
       (pi_r.*pi_d.*pi_c.*pi_b.*pi_t.*pi_n);

if p0p9 > 1
  disp('Primary nozzle not choked')
  p9p0 = 1;
  p0p9 = 1;
else
  p9p0 = 1./p0p9;
end

% equations continued
pt9p9 = (1./p9p0).*pi_r.*pi_d.*pi_c.*pi_b.*pi_t.*pi_n;

if ab == 0
  tau_lam_ab = tau_lam.*tau_t;
end

t9t0 = ((cp_c./cp_ab).*tau_lam_ab)./ ...
       ((pt9p9).^((gam_ab-1)./gam_ab));

term1 = 2./(gam_c-1);
term2 = tau_lam_ab;
term3 = 1 - (pt9p9).^(-((gam_ab-1)./gam_ab));

m0u9u0 = (term1.*term2.*term3).^(1./2);

pt92p92 = (1./p92p0).*pi_r.*pi_d.*pi_c2.*pi_n2;

if ab2 == 0
  tau_lam_ab2 = tau_r.*tau_c2;
end

t92t0 = ((cp_c./cp_ab2).*tau_lam_ab2)./ ...
	((pt92p92).^((gam_ab2-1)./gam_ab2));

term4 = 2./(gam_c-1);
term5 = tau_lam_ab2;
term6 = 1 - (pt92p92).^(-((gam_ab2-1)./gam_ab2));

m0u92u0 = (term4.*term5.*term6).^(1./2);

f_ab = (1+f).* ...
       ((tau_lam_ab - tau_lam.*tau_t)./ ...
	((h.*eta_ab./(cp_c.*t0)) - tau_lam_ab));

f_ab2 = (tau_lam_ab2 - tau_r.*tau_c2)./ ...
	((h.*eta_ab2./(cp_c.*t0)) - tau_lam_ab2);

f_mdot = (a0./(1+alpha)).* ...
	 ((1+f+f_ab).*(m0u9u0) - ...
	  m0 + ...
	  (1+f+f_ab).*(1./(gam_c.*(m0u9u0))).* ...
	  (t9t0).*(1-(1./p9p0)) + ...
	  alpha.* ...
	  ((1+f_ab2).*(m0u92u0) - ...
	   m0 + ...
	   (1+f_ab2).*(1./(gam_c.*(m0u92u0))).* ...
	   (t92t0).*(1-(1./p92p0))));

s = ((f + f_ab + alpha.*f_ab2)./...
     ((1+alpha).*f_mdot)).*(10.^6);

if verbose == 1
  pi_ch = pi_c./pi_c2;
  tau_ch = tau_c./tau_c2;

  eta_c2 = (pi_c2.^((gam_c-1)./gam_c) - 1)./(tau_c2 - 1);
  eta_ch = (pi_ch.^((gam_c-1)./gam_c) - 1)./(tau_ch - 1);
  eta_t = (1 - tau_t)./(1 - pi_t.^((gam_t-1)./gam_t));

  disp(' ')
  disp(['Etta ch: ',num2str(eta_ch)])
  disp(['Etta c'': ',num2str(eta_c2)])

  disp(' ')
  disp(['M0*u9/u0: ',num2str(m0u9u0)])
  disp(['M0*u9''/u0: ',num2str(m0u92u0)])

  disp(' ')
  disp(['Turbine temperature ratio: ',num2str(tau_t)])
  disp(['Turbine efficiency: ',num2str(eta_t)])
  disp(['Fan temperature ratio: ',num2str(tau_c2)])
  disp(['Fuel-to-air ratio: ',num2str(f)])

  disp(' ')
end

end
