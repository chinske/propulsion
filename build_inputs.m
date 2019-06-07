function inputs = build_inputs()
% BUILD_INPUTS Build Inputs Structure
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
% INPUTS = BUILD_INPUTS

disp(' ')

inputs.ab = input('Primary stream afterburning flag: ');
inputs.ab2 = input('Secondary stream afterburning flag: ');

inputs.t0 = input('T0 (K): ');

inputs.gam_c = input('Ratio of specific heats, compressor: ');
inputs.gam_t = input('Ratio of specific heats, turbine: ');
if inputs.ab == 1
  inputs.gam_ab = ...
  input('Ratio of specific heats, PRI afterburner: ');
end
if inputs.ab2 == 1
  inputs.gam_ab2 = ...
  input('Ratio of specific heats, SEC afterburner: ');
end

inputs.cp_c = ...
input('Specific heat of air Cp (J/(kg*K)), compressor: ');
inputs.cp_t = ...
input('Specific heat of air Cp (J/(kg*K)), turbine: ');
if inputs.ab == 1
  inputs.cp_ab = ...
  input('Specific heat of air Cp (J/(kg*K)), PRI afterburner: ');
end
if inputs.ab2 == 1
  inputs.cp_ab2 = ...
  input('Specific heat of air Cp (J/(kg*K)), SEC afterburner: ');
end

inputs.h = input('Fuel heating value (J/kg): ');

inputs.pi_d = input('Pressure ratio, inlet: ');
inputs.pi_b = input('Pressure ratio, burner: ');
inputs.pi_n = input('Pressure ratio, PRI nozzle: ');
inputs.pi_n2 = input('Pressure ratio, SEC nozzle: ');

inputs.eta_b = input('Efficiency, burner: ');
if inputs.ab == 1
  inputs.eta_ab = input('Efficiency, PRI afterburner: ');
end
if inputs.ab2 == 1
  inputs.eta_ab2 = input('Efficiency, SEC afterburner: ');
end
inputs.eta_m = input('Efficiency, mechanical: ');

inputs.e_c = input('Polytropic efficiency, PRI compressor: ');
inputs.e_c2 = input('Polytropic efficiency, SEC compressor: ');
inputs.e_t = input('Polytropic efficiency, turbine: ');

inputs.p9p0 = input('p9/p0, PRI: ');
inputs.p92p0 = input('p9/p0, SEC: ');

inputs.tau_lam = input('tau_lam: ');
if inputs.ab == 1
  inputs.tau_lam_ab = input('PRI tau_lam_ab: ');
end
if inputs.ab2 == 1
  inputs.tau_lam_ab2 = input('SEC tau_lam_ab: ');
end

inputs.pi_c = input('Compressor pressure ratio: ');
inputs.pi_c2 = input('Fan pressure ratio: ');

inputs.m0 = input('Flight Mach number: ');
inputs.alpha = input('Bypass ratio: ');

end
