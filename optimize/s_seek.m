% S_SEEK Script
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
% This script finds an engine (set of parameters) to match a given
% specific fuel consumption.  Free parameters are defined within this
% script.  Fixed parameters and the given specific fuel consumption
% are defined in the objective function PHI.
% 
% This script can be used to find a feasible set of parameters that
% match the known, on-design performance of a particular engine.
%

% ------------------------------------------------------------
% Parameter Bounds
% Define lower bound and upper bound for each free parameter.
% See Table 6.3 in Oates for typical ranges of parameters.
% ------------------------------------------------------------

lb_pi_d = 0.98;
ub_pi_d = 0.998;

lb_pi_b = 0.93;
ub_pi_b = 0.98;

lb_pi_n = 0.99;
ub_pi_n = 0.998;

lb_pi_n2 = 0.99;
ub_pi_n2 = 0.998;

lb_eta_b = 0.96;
ub_eta_b = 0.998;

lb_e_c = 0.86;
ub_e_c = 0.94;

lb_e_c2 = 0.85;
ub_e_c2 = 0.92;

lb_e_t = 0.85;
ub_e_t = 0.92;

lb_tau_lam = 4.4;
ub_tau_lam = 5.4;

lb_pi_c2 = 1.1;
ub_pi_c2 = 2.0;

% build lower bounds (LB) and upper bounds (UB) vectors
lb = [lb_pi_d; ...
      lb_pi_b; ...
      lb_pi_n; ...
      lb_pi_n2; ...
      lb_eta_b; ...
      lb_e_c; ...
      lb_e_c2; ...
      lb_e_t; ...
      lb_tau_lam; ...
      lb_pi_c2];

ub = [ub_pi_d; ...
      ub_pi_b; ...
      ub_pi_n; ...
      ub_pi_n2; ...
      ub_eta_b; ...
      ub_e_c; ...
      ub_e_c2; ...
      ub_e_t; ...
      ub_tau_lam; ...
      ub_pi_c2];

% define initial guess
x0 = (lb+ub)./2;

% ------------------------------------------------------------
% Call SQP
% ------------------------------------------------------------

x = sqp(x0,@phi,[],[],lb,ub);
