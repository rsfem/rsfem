% Copyright (C) 2022 M. A. Yalçın, İ. Temizer
% This file is distributed under the terms of the GNU General Public License.
% See the file COPYING for license details.

function [num]= integrator_two_matrices(K1,K2,E,GQ,  GQ_weights,value_PHI, value_diff_PHI,connectivity,coord,NumElmNods,LLQ,LLQ_weights,LLQ_points)
% spherical integral calculator using Lebedev Laikov grid
num = 0;
for element = 1:E  
    for k = 1:GQ 
        K1_sphere = K1(GQ*(element-1)+k,:);
        K2_sphere = K2(GQ*(element-1)+k,:);
        r = 0 ;
        F_s = 0;
        for i = 1:NumElmNods
          r = r + value_PHI(i,k) * coord(connectivity(element,i));
          F_s = F_s + value_diff_PHI(i,k) * coord(connectivity(element,i));
        end
        for l = 1:LLQ
            value_k = K1_sphere(1,l)*K2_sphere(1,l);
            num = num + value_k *  F_s * GQ_weights(k) * LLQ_weights(l)*4*pi*r^2;
        end
    
    end
end

