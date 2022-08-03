% Copyright (C) 2022 M. A. Yalçın, İ. Temizer
% This file is distributed under the terms of the GNU General Public License.
% See the file COPYING for license details.

function [rho_on_sphere rho_c d_rho_c] = InitialRhoGuess(E,GQ,connectivity,value_PHI,coord,NumElmNods,  GQ_weights,value_diff_PHI,LLQ,LLQ_points,LLQ_weights,Z,Zion,NLCC_flag,exc_flag,pp_flag)
%spherical integral rho*dr
num_ip = 0;
rho_on_sphere = [];
if NLCC_flag == 1 
    num_desired = Zion;
    rho_c = zeros(GQ*E,LLQ);
    d_rho_c = zeros(GQ*E,LLQ,3);
    d_rho_c_u = zeros(GQ*E,LLQ);
else
    if pp_flag == 1 
        num_desired = Zion;
        rho_c = [];
        d_rho_c = [];
    elseif pp_flag == 0 
        num_desired = Z;
        rho_c = [];
        d_rho_c = [];
    end
end
for element = 1:E
    for k = 1:GQ            
        r = 0 ;
        F_s_r = 0;
        for i = 1:NumElmNods
          r = r + value_PHI(i,k) * coord(connectivity(element,i));
          F_s_r = F_s_r + value_diff_PHI(i,k) * coord(connectivity(element,i));
        end
        for l =1:LLQ
            if NLCC_flag == 1 && pp_flag == 1  
                [rho_c(GQ*(element-1)+k,l) d_rho_c_u(GQ*(element-1)+k,l)]= coreDensity(Z,Zion,r);
                d_rho_c(GQ*(element-1)+k,l,1) = d_rho_c(GQ*(element-1)+k,l,1) + d_rho_c_u(GQ*(element-1)+k,l)*LLQ_points(l,1); 
                d_rho_c(GQ*(element-1)+k,l,2) = d_rho_c(GQ*(element-1)+k,l,2) + d_rho_c_u(GQ*(element-1)+k,l)*LLQ_points(l,2);
                d_rho_c(GQ*(element-1)+k,l,3) = d_rho_c(GQ*(element-1)+k,l,3) + d_rho_c_u(GQ*(element-1)+k,l)*LLQ_points(l,3);
            end
            value_k = exp(-2*r);
            rho_on_ip(l) = value_k;
            num_ip = num_ip + value_k * F_s_r *LLQ_weights(l)*4*pi * GQ_weights(k)*r^2; 
        end
        rho_on_sphere=[rho_on_sphere ;rho_on_ip];
    end
end
rho_on_sphere = rho_on_sphere ./num_ip * num_desired; 
