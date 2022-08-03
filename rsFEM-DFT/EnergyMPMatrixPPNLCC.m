% Copyright (C) 2022 M. A. Yalçın, İ. Temizer
% This file is distributed under the terms of the GNU General Public License.
% See the file COPYING for license details.

function [E num_exc num_vh  num_T num_f2w eps_ent]= EnergyMPMatrixPPNLCC(eps,v_h_full, v_h_in_full,rho, rho_c, E, GQ,rho_in_old, rho_c_in_old, GQ_weights,value_PHI, value_diff_PHI,connectivity,coord,NumElmNods,Z,Zion,fractional_occupancy,LLQ,LLQ_points,LLQ_weights,L_cutoff,N,grad_rho_v,grad_rho_old,grad_rho_c,NLCC_flag,exc_flag,pp_flag,Ss,Clm,eps_f,L_poisson) %% all inputs should be vectors; fo as input
% E(psi) = (2*SUM( fo *eps_i)-<rho*(v_h_old + v_xc_old)>) + 1/2*<rho*v_h>+<rho*eps_xc> 
% v_eff = v_h + v_xc ;
% E_self == E_ext if E_nn == 0; otherwise E_ext = E_self + E_nn

num_exc = 0;
num_Vext = 0;
num_f2xc = 0;
num_f2w = 0;
num_f2vh = 0;
num_f3 = 0;
num_f4 = 0;
f1 = 0;

for element = 1:E
	for k = 1:GQ            
        r = 0 ;
        F_s = 0;
        for i = 1:NumElmNods
            r = r + value_PHI(i,k) * coord(connectivity(element,i));
            F_s = F_s + value_diff_PHI(i,k) * coord(connectivity(element,i));
        end
        orb_v_h_in = zeros(L_poisson^2);
        orb_v_h = zeros(L_poisson^2);
        for orb = 1:L_poisson^2
            v_h_in = v_h_in_full((orb-1)*N+1:orb*N);
            v_h =    v_h_full((orb-1)*N+1:orb*N);  
            for i = 1:NumElmNods
                orb_v_h_in(orb) = orb_v_h_in(orb) + value_PHI(i,k) * v_h_in(connectivity(element,i));
                orb_v_h(orb)    = orb_v_h(orb)    + value_PHI(i,k) * v_h(connectivity(element,i));
            end
        end
        for l = 1:LLQ
            rho_one_sphere = rho(:,l);
            rho_on_point = rho_one_sphere(GQ*(element-1)+k);
            new_v_h_in = 0;
            new_v_h = 0;
            for orb = 1:L_poisson^2
                sh = evalSphericalHarmonic_hb(orb,0,r*LLQ_points(l,1), r*LLQ_points(l,2), r*LLQ_points(l,3),"counts","wo_diff","cartesian");
                new_v_h_in = new_v_h_in + orb_v_h_in(orb) * sh; 
                new_v_h    = new_v_h + orb_v_h(orb)       * sh;
            end
            value_f2_vh = rho_on_point * new_v_h_in ;
            value_f3 = 1/2 * rho_on_point * new_v_h;
            num_f3 = num_f3 + value_f3 * F_s * GQ_weights(k)*LLQ_weights(l)*4*pi*r^2;
            num_f2vh = num_f2vh + value_f2_vh * F_s * GQ_weights(k)*LLQ_weights(l)*4*pi*r^2;
        end
    end
end
    
if  exc_flag == 1
    [~,exc,~] = exchangeCorrelation(rho,rho_c,grad_rho_v,grad_rho_c,exc_flag,pp_flag,NLCC_flag);
    [vcx_old,~,~] = exchangeCorrelation(rho_in_old,rho_c,grad_rho_v,grad_rho_c,exc_flag,pp_flag,NLCC_flag);
elseif exc_flag == 2 
    [~,exc,~] = exchangeCorrelation(rho,rho_c,grad_rho_v,grad_rho_c,exc_flag,pp_flag,NLCC_flag);
    [vcx_old,~,wxc_old] = exchangeCorrelation(rho_in_old,rho_c,grad_rho_old,grad_rho_c,exc_flag,pp_flag,NLCC_flag);
    if pp_flag == 0
        grad_rho_c = zeros(E*GQ,LLQ);
        grad_rho = grad_rho_v + grad_rho_c;
    else
        grad_rho = grad_rho_v + grad_rho_c;
    end
end

for element = 1:E
    for k = 1:GQ            
        r = 0 ;
        F_s = 0;
        for i = 1:NumElmNods
          r = r + value_PHI(i,k) * coord(connectivity(element,i));
          F_s = F_s + value_diff_PHI(i,k) * coord(connectivity(element,i));
        end
        for l = 1:LLQ
            rho_point = rho(GQ*(element-1)+k,l);
            if pp_flag == 1
                rho_c_point = rho_c(GQ*(element-1)+k,l);
                if exc_flag == 2
                    grad_rho_point = grad_rho_v(GQ*(element-1)+k,l,:);
                    value_gga_wxc = dot(grad_rho_point,wxc_old(GQ*(element-1)+k,l,:));
                    num_f2w = num_f2w + value_gga_wxc * F_s * GQ_weights(k)*LLQ_weights(l)*4*pi*r^2;
                end
                value_f4 = (rho_point+rho_c_point)*(exc(GQ*(element-1)+k,l)) ; 
                Z_in = Zion;
            else
                if exc_flag == 2
                    grad_rho_point = grad_rho_v(GQ*(element-1)+k,l,:);
                    value_gga_wxc = dot(grad_rho_point,wxc_old(GQ*(element-1)+k,l,:));
                    num_f2w = num_f2w + value_gga_wxc * F_s * GQ_weights(k)*LLQ_weights(l)*4*pi*r^2;
                end
                value_f4 = rho_point*(exc(GQ*(element-1)+k,l)) ;    
                Z_in = Z;
            end
            value_f2_xc = rho_point*(vcx_old(GQ*(element-1)+k,l)) ;
            num_f4 = num_f4 + value_f4 * F_s * GQ_weights(k)*LLQ_weights(l)*4*pi*r^2;         
            num_f2xc = num_f2xc + value_f2_xc * F_s * GQ_weights(k)*LLQ_weights(l)*4*pi*r^2;  
            num_exc = num_exc + value_f4 * F_s * GQ_weights(k)*LLQ_weights(l)*4*pi*(r^2);
            num_Vext   = num_Vext   + rho_point * -Z_in/r * F_s * GQ_weights(k)*LLQ_weights(l)*4*pi*r^2;
        end
    end
end

% num_f2vh = rho *vh
% num_f3 = 1/2 rho*vh
% num_f4 = rho * exc
num = -(num_f2xc + num_f2vh + num_f2w) + num_f3 + num_f4;
for i = 1:length(fractional_occupancy)
    f1 = f1 + fractional_occupancy(i)*eps(i,i);
end
num_f2 = num_f2xc + num_f2vh + num_f2w;
num_T = 2*f1-num_f2 -num_Vext;
num_vh = num_f3;
E = 2*f1 + num;
eps_ent = E -Ss;
end
