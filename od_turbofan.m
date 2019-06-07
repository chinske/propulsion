function [s,f_rat] = od_turbofan(inputs)
% OD_TURBOFAN Off-Design Turbofan
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
% [S,F_RAT] = OD_TURBOFAN(INPUTS)

% reassign input variables
gam_t = inputs.gam_t;
cp_c = inputs.cp_c;
cp_t = inputs.cp_t;
eta_ch = inputs.eta_ch;
eta_c2 = inputs.eta_c2;
eta_b = inputs.eta_b;
pi_d = inputs.pi_d;
pi_b = inputs.pi_b;
pi_n = inputs.pi_n;
pi_n2 = inputs.pi_n2;
t0 = inputs.t0;
p0p0r = inputs.p0p0r;
h = inputs.h;
tau_lam = inputs.tau_lam;
m0 = inputs.m0;
m0r = inputs.m0r;
tau_lamr = inputs.tau_lamr;
pi_c2r = inputs.pi_c2r;
pi_cr = inputs.pi_cr;
alphar = inputs.alphar;
eta_chr = inputs.eta_chr;
eta_c2r = inputs.eta_c2r;
eta_br = inputs.eta_br;
pi_dr = inputs.pi_dr;
pi_br = inputs.pi_br;
pi_nr = inputs.pi_nr;
pi_n2r = inputs.pi_n2r;
t0r = inputs.t0r;
m0u9u0r = inputs.m0u9u0r;
m0u92u0r = inputs.m0u92u0r;
tau_tr = inputs.tau_tr;
eta_tr = inputs.eta_tr;
tau_c2r = inputs.tau_c2r;
fr = inputs.fr;

% equations
tau_r = 1 + (m0.^2)./5;
pi_r = tau_r.^3.5;
tau_rr = 1 + (m0r.^2)./5;
pi_rr = tau_rr.^3.5;
a0 = 20.04.*sqrt(t0);
a0r = 20.04.*sqrt(t0r);

% power balance
iter = 0;
es = 0.00005;
ea = 100;
pi_chr = pi_cr./pi_c2r;
tau_c2 = tau_c2r;
while ea > es
  solold = tau_c2;

  pi_ch = (1 + (eta_ch./eta_chr) .* ...
	       ((pi_cr./pi_c2r).^(1./3.5) - 1) .* ...
	       ((tau_lam./(tau_r.*tau_c2)) ./ ...
		(tau_lamr./(tau_rr.*tau_c2r)))).^3.5;

  alpha = alphar.*(pi_chr./pi_ch).* ...
	  ((tau_lam./(tau_r.*tau_c2)) ./ ...
	   (tau_lamr./(tau_rr.*tau_c2r))).^(1./2);

  tau_c2 = 1 + (tau_c2r - 1).*((alphar+1)./(alpha+1)).* ...
	       ((tau_lam./tau_r)./(tau_lamr./tau_rr));

  iter = iter + 1;
  ea = abs((tau_c2 - solold)./tau_c2).*100;
end

disp(['Bypass ratio (off-design): ',num2str(alpha)])
disp(['Fan temperature ratio (off-design): ',num2str(tau_c2)])
disp(['Number of iterations: ',num2str(iter)])

% equations continued
tau_t = tau_tr;

pi_tr = (1 + (1./eta_tr).*(tau_tr - 1)).^(gam_t./(gam_t-1));
pi_t = pi_tr;

tau_c = tau_c2.*(1 + (1./eta_ch).*(pi_ch.^(1./3.5) - 1));

pi_c2 = (1 + eta_c2.*(tau_c2 - 1)).^3.5;

f = (tau_lam - tau_r.*tau_c)./((h.*eta_b./(cp_c.*t0)) - tau_lam);

pt9p9 = pi_r.*pi_d.*pi_c2.*pi_ch.*pi_b.*pi_t.*pi_n;

m0u9u0 = (5 .* tau_lam .* tau_t .* ...
	  (1 - pt9p9.^(-((gam_t-1)./gam_t)) )).^(1./2);

pt92p9 = pi_r.*pi_d.*pi_c2.*pi_n2;

m0u92u0 = (5 .* tau_r .* tau_c2 .* ...
	   (1 - pt92p9.^(-1./3.5))).^(1./2);

f_mdot = (a0./(1+alpha)) .* ...
	 ((1+f).*m0u9u0 - m0 + alpha.*(m0u92u0 - m0));

s = (f./((1+alpha).*f_mdot)).*10.^6;

term1 = alphar./alpha;
term2 = (pi_r.*pi_d.*pi_c2.*pi_n2)./ ...
	(pi_rr.*pi_dr.*pi_c2r.*pi_n2r);
term3 = ((tau_rr.*tau_c2r)./(tau_r.*tau_c2)).^(1./2);
term4 = p0p0r;
term5 = ((1+f).*m0u9u0 - m0 + alpha.*(m0u92u0 - m0))./ ...
	((1+fr).*m0u9u0r - m0r + alphar.*(m0u92u0r - m0r));

f_rat = term1.*term2.*term3.*term4.*term5;

end
