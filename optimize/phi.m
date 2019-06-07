function s_err = phi(x)
% PHI Objective Function for S_SEEK
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
% S_ERR = PHI(X) returns S_ERR, the difference between the specific
% fuel consumption computed with the parameters X and the desired
% specific fuel consumption.

% desired specific fuel consumption
s_goal = 9.3473;

% fixed parameters
inputs.t0 = 303.15;
inputs.gam_c = 1.4;
inputs.gam_t = 1.35;
inputs.cp_c = 1004.9;
inputs.cp_t = 1004.9;
inputs.h = 4.4194E7;
inputs.m0 = 0;

inputs.eta_m = 1;
inputs.pi_c = 25.70;
inputs.alpha = 4.80;

% build inputs structure array
inputs.pi_d = x(1);
inputs.pi_b = x(2);
inputs.pi_n = x(3);
inputs.pi_n2 = x(4);
inputs.eta_b = x(5);
inputs.e_c = x(6);
inputs.e_c2 = x(7);
inputs.e_t = x(8);
inputs.tau_lam = x(9);
inputs.pi_c2 = x(10);

% call cycle analysis function
[f_mdot,s,inputs] = nonideal_turbofan2(inputs,0);

% compute error
s_err = abs(s - s_goal);

end
