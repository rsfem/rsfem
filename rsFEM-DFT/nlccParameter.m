% Copyright (C) 2022 M. A. Yalçın, İ. Temizer
% This file is distributed under the terms of the GNU General Public License.
% See the file COPYING for license details.

function [c_core r_core] = nlccParameter(Z,Zion)
if Z == 3 && Zion == 1 % Li
    r_core = 0.7680116744624;
    c_core =  0.28730453041;
elseif Z == 13 && Zion == 3 % Al
    r_core =   0.4302905414;
    c_core =   0.3583342894;
elseif Z ==1 && Zion == 1 % H  
    r_core =   0.0;
    c_core =   0.0;
end
