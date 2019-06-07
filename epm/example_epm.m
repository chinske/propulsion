function [s,f] = example_epm(atmos,altitude,mach,throttle_fraction)
% EXAMPLE_EPM Example Engine Performance Modeling
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
% [S, F] = EXAMPLE_EPM(ATMOS, ALTITUDE, MACH, THROTTLE_FRACTION)
% 
% Engine Performance Modeling, CFM 56-3C1
% 
% Inputs
% ------
% ATMOS is a matrix (lookup-table) with a reference atmosphere.
% Column 1: Altitude (km)
% Column 3: Pressure (Pa)
% Column 4: Temperature (K)
% 
% ALTITUDE (ft)
% MACH
% 
% THROTTLE_FRACTION is a value between 0 and 1, where 1 implies
% maximum TT4.
% 
% Outputs
% -------
% S: Specific Fuel Consumption [ (lbm fuel/h)/(lbf thrust) ]
% F: Thrust [lbf]
% 

% assign atmos data
atmos_h = atmos(:,1);
atmos_p = atmos(:,3);
atmos_t = atmos(:,4);

% convert altitude to km
altitude = altitude.*0.3048E-3;

% build inputs
inputs.gam_t = 1.35;
inputs.cp_c = 1004.9;
inputs.cp_t = 1004.9;
inputs.eta_ch = 0.87953;
inputs.eta_c2 = 0.88515;
inputs.eta_b = 0.98608;
inputs.pi_d = 0.99206;
inputs.pi_b = 0.96201;
inputs.pi_n = 0.99535;
inputs.pi_n2 = 0.99520;
inputs.h = 4.4194E7;
inputs.m0r = 0;
inputs.tau_lamr = 4.8931;
inputs.pi_c2r = 1.5383;
inputs.pi_cr = 25.700;
inputs.alphar = 4.80;
inputs.eta_chr = 0.87953;
inputs.eta_c2r = 0.88515;
inputs.eta_br = 0.98608;
inputs.pi_dr = 0.99206;
inputs.pi_br = 0.96201;
inputs.pi_nr = 0.99535;
inputs.pi_n2r = 0.99520;
inputs.t0r = 303.15;
inputs.m0u9u0r = 0.85925;
inputs.m0u92u0r = 0.80373;
inputs.tau_tr = 0.50274;
inputs.eta_tr = 0.92886;
inputs.tau_c2r = 1.1479;
inputs.fr = 0.015435;

% performance conditions
f_ref = 23500;
tt4_max = 1483.3;
tt4 = tt4_max.*throttle_fraction;
inputs.t0 = interp1(atmos_h,atmos_t,altitude);
inputs.p0p0r = interp1(atmos_h,atmos_p,altitude)./101325;
inputs.tau_lam = tt4./inputs.t0;
inputs.m0 = mach;

% call od_turbofan
[s,f_rat] = od_turbofan(inputs);

% outputs
s = s./28.325;
f = f_rat.*f_ref;

end
