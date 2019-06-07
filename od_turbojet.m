function [f_rat,s_rat,a9_rat,mdot_c_rat,pi_c_rat] = ...
	 od_turbojet(inputs)
% OD_TURBOJET Off-Design Fixed-Area Turbojet
% Copyright 2016 Christopher Chinske
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
% [F_RAT,S_RAT,A9_RAT,MDOT_C_RAT,PI_C_RAT] = ...
% OD_TURBOJET(INPUTS)

% reassign input variables
gam_c = inputs.gam_c;
gam_t = inputs.gam_t;
p0p0r = inputs.p0p0r;
t0t0r = inputs.t0t0r;
eta_m = inputs.eta_m;
eta_c = inputs.eta_c;
eta_t = inputs.eta_t;
pi_cr = inputs.pi_cr;
pi_dr = inputs.pi_dr;
pi_br = inputs.pi_br;
pi_nr = inputs.pi_nr;
tau_lamr = inputs.tau_lamr;
m0r = inputs.m0r;
pi_d = inputs.pi_d;
tau_lam = inputs.tau_lam;
m0 = inputs.m0;

% equations
tau_rr = 1 + ((gam_c-1)./2).*m0r.^2;
pi_rr = tau_rr.^(gam_c./(gam_c-1));

tau_r = 1 + ((gam_c-1)./2).*m0.^2;
pi_r = tau_r.^(gam_c./(gam_c-1));

tau_cr = 1 + (1./eta_c).*(pi_cr.^((gam_c-1)./gam_c) - 1);
tau_c = 1 + (tau_cr - 1).*(tau_lam./tau_lamr).*(tau_rr./tau_r);
pi_c = (1 + eta_c.*(tau_c - 1)).^(gam_c./(gam_c-1));

tau_tr = 1 - (1./eta_m).*(tau_rr./tau_lamr).*(tau_cr - 1);
tau_t = tau_tr;
pi_tr = (1 - (1./eta_t).*(1-tau_tr)).^(gam_t./(gam_t-1));
pi_t = pi_tr;

pi_b = pi_br;
pi_n = pi_nr;

pt9p0r = pi_b.*pi_t.*pi_n.*pi_rr.*pi_dr.*pi_cr;

pt9p9 = pt9p0r.*(pi_r.*pi_d.*pi_c)./(pi_rr.*pi_dr.*pi_cr);

m0u9u0 = ((2./(gam_c-1)).*tau_lam.*tau_t.* ...
	  (1 - (pt9p9).^(-(gam_t-1)./gam_t))).^(1./2);

% assume p9 = p0 =>
% (pt9/p9)_R = (pt9/p0)_R
% pt9/p9 = pt9/p0
pt9p9r = pt9p0r;
pt9p0 = pt9p9;

m0u9u0r = ((2./(gam_c-1)).*tau_lamr.*tau_tr.* ...
	  (1 - (pt9p9r).^(-(gam_t-1)./gam_t))).^(1./2);

% outputs
f_rat = (pi_r.*pi_d.*pi_c./(pi_rr.*pi_dr.*pi_cr)) .* ...
	p0p0r .* ...
	((tau_lamr./tau_lam).^(1./2)) .* ...
	((m0u9u0 - m0)./(m0u9u0r - m0r));

s_rat = (t0t0r.^(1./2)) .* ...
	((m0u9u0r - m0r)./(m0u9u0 - m0)) .* ...
	((tau_lam - tau_r.*tau_c)./(tau_lamr - tau_rr.*tau_cr));

a9_rat = ((pt9p0./pt9p0r).^((gam_t+1)./(2.*gam_t))) .* ...
	 ((pt9p0r.^((gam_t-1)./gam_t) - 1) ./ ...
	  (pt9p0.^((gam_t-1)./gam_t) - 1)).^(1./2);

mdot_c_rat = (pi_c./pi_cr).*((tau_cr-1)./(tau_c-1)).^(1./2);

pi_c_rat = pi_c./pi_cr;

end
