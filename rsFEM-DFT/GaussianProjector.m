% Copyright (C) 2022 M. A. Yalçın, İ. Temizer
% This file is distributed under the terms of the GNU General Public License.
% See the file COPYING for license details.

function [p] = GaussianProjector(r,l,i,ppP)
rl = ppP((l+1)*2,1);
g = gamma(l + (4*i-1)/2);
p = sqrt(2)*r^(l+2*(i-1))*exp(-r^2/2/rl^2)/(rl^(l+(4*i-1)/2)*sqrt(g));
end
