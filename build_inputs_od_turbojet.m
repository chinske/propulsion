function inputs = build_inputs_od_turbojet()
% BUILD_INPUTS_OD_TURBOJET Build Inputs Off-Design Turbojet
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
% INPUTS = BUILD_INPUTS_OD_TURBOJET

disp(' ')

inputs.gam_c = input('Ratio of specific heats, compressor: ');
inputs.gam_t = input('Ratio of specific heats, turbine: ');

inputs.p0p0r = input('p0/p0_R: ');
inputs.t0t0r = input('t0/t0_R: ');

inputs.eta_m = input('Efficiency, mechanical: ');
inputs.eta_c = input('Efficiency, compressor: ');
inputs.eta_t = input('Efficiency, turbine: ');

inputs.pi_cr = input('Pressure ratio, compressor (Ref.): ');
inputs.pi_dr = input('Pressure ratio, inlet (Ref.): ');
inputs.pi_br = input('Pressure ratio, burner (Ref.): ');
inputs.pi_nr = input('Pressure ratio, nozzle (Ref.): ');
inputs.tau_lamr = input('tau_lam (Ref.): ');
inputs.m0r = input('M0 (Ref.): ');

inputs.pi_d = input('Pressure ratio, inlet: ');
inputs.tau_lam = input('tau_lam: ');

inputs.m0 = input('Flight Mach number: ');

end
