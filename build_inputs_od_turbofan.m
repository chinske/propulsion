function inputs = build_inputs_od_turbofan()
% BUILD_INPUTS_OD_TURBOFAN Build Inputs Off-Design Turbofan
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
% INPUTS = BUILD_INPUTS_OD_TURBOFAN

disp(' ')

inputs.gam_t = input('Ratio of specific heats, turbine: ');

inputs.cp_c = ...
input('Specific heat of air Cp (J/(kg*K)), compressor: ');
inputs.cp_t = ...
input('Specific heat of air Cp (J/(kg*K)), turbine: ');

disp(' ')
disp('Off-Design Values')
disp('-----------------')

inputs.eta_ch = input('Etta ch (high-pressure compressor): ');
inputs.eta_c2 = input('Etta c'' (fan): ');
inputs.eta_b = input('Etta b (burner): ');

inputs.pi_d = input('Pressure ratio, inlet: ');
inputs.pi_b = input('Pressure ratio, burner: ');
inputs.pi_n = input('Pressure ratio, PRI nozzle: ');
inputs.pi_n2 = input('Pressure ratio, SEC nozzle: ');

inputs.t0 = input('T0 (K): ');
inputs.p0p0r = input('P0/P0_R: ');

inputs.h = input('Fuel heating value (J/kg): ');
inputs.tau_lam = input('tau_lam: ');

inputs.m0 = input('Flight Mach number: ');

disp(' ')
disp('Reference Values')
disp('----------------')

inputs.m0r = input('M0 (Ref.): ');

disp(' ')
inputs.tau_lamr = input('tau_lam (Ref.): ');
inputs.pi_c2r = input('Pressure ratio, fan (Ref.): ');
inputs.pi_cr = input('Pressure ratio, compressor (Ref.): ');

inputs.alphar = input('Bypass ratio (Ref.): ');

disp(' ')
inputs.eta_chr = input('Etta ch (high-pressure compressor) (Ref.): ');
inputs.eta_c2r = input('Etta c'' (fan) (Ref.): ');
inputs.eta_br = input('Etta b (burner) (Ref.): ');

disp(' ')
inputs.pi_dr = input('Pressure ratio, inlet (Ref.): ');
inputs.pi_br = input('Pressure ratio, burner (Ref.): ');
inputs.pi_nr = input('Pressure ratio, PRI nozzle (Ref.): ');
inputs.pi_n2r = input('Pressure ratio, SEC nozzle (Ref.): ');

disp(' ')
inputs.t0r = input('T0 (Ref.) (K): ');

inputs.m0u9u0r = input('M0*u9/u0 (Ref.): ');
inputs.m0u92u0r = input('MO*u9''/u0 (Ref.): ');

disp(' ')
inputs.tau_tr = input('Turbine temperature ratio (Ref.): ');
inputs.eta_tr = input('Turbine efficiency (Ref.): ');
inputs.tau_c2r = input('Fan temperature ratio (Ref.): ');
inputs.fr = input('Fuel-to-air ratio (Ref.): ');

disp(' ')

end
