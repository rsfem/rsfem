% Copyright (C) 2022 M. A. Yalçın, İ. Temizer
% This file is distributed under the terms of the GNU General Public License.
% See the file COPYING for license details.

function [rho_c d_rho_c] = coreDensity(Z,Zion,r)
[c_core r_core] = nlccParameter(Z,Zion);
rho_c = c_core*(Z-Zion)/(sqrt(2*pi)*r_core)^3*exp((-r^2/2/r_core^2));
d_rho_c = -(2^(1/2)*c_core*r*exp(-r^2/(2*r_core^2))*(Z - Zion))/(4*pi^(3/2)*r_core^5);

