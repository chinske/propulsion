function [pt7pt5,m7,m32] = ...
	 ideal_ca_mixer(pt32pt5,tt32tt5,m5,alpha,gam)
% IDEAL_CA_MIXER Ideal Constant-Area Mixer
% Copyright 2014 Christopher Chinske
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
% [PT7PT5, M7, M32] = ...
% IDEAL_CA_MIXER(PT32PT5, TT32TT5, M5, ALPHA, GAM)
% 
% PT32PT5 = p_t_3' / p_t_5
% TT32TT5 = T_t_3' / T_t_5
% 

term1 = 2./(gam-1);
term2 = pt32pt5.^((gam-1)./gam);
term3 = (1 + ((gam-1)./2).*m5.^2);

m32 = sqrt( term1 .* (term2.*term3 - 1) );

a32a5 = alpha.*sqrt(tt32tt5).*(m5./m32).* ...
	pt32pt5.^-((gam-1)./(2.*gam));

f = (1 + alpha) .* ...
    (1 + alpha.*tt32tt5) .* ...
    (1./sqrt(fn(m5,gam)) + ...
     alpha.*sqrt(tt32tt5)./sqrt(fn(m32,gam))).^(-2);

m7 = (2.*f./( 1 - 2.*gam.*f + sqrt(1 - 2.*(gam+1).*f) )).^(1./2);

term1 = ((1+alpha).*(1 + alpha.*tt32tt5)).^(1./2);
term2 = 1 + a32a5;
term3 = m5./m7;
term4 = 1 + ((gam-1)./2).*m7.^2;
term5 = 1 + ((gam-1)./2).*m5.^2;

pt7pt5 = (term1./term2).*term3.*(term4./term5).^ ...
				((gam+1)./(2.*(gam-1)));

end

function f = fn(M,gam)
f = (M.^2) .* (1 + ((gam-1)./2).*M.^2) .* (1 + gam.*M.^2).^(-2);
end
