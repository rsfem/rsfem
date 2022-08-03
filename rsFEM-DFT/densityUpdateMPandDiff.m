% Copyright (C) 2022 M. A. Yalçın, İ. Temizer
% This file is distributed under the terms of the GNU General Public License.
% See the file COPYING for license details.

function [density grad_density r_on_ip] = densityUpdateMPandDiff(E,N,GQ,value_PHI,value_diff_PHI,coord,Clm,connectivity,fractional_occupancy,NumElmNods,LLQ,LLQ_points,L_cutoff,exc_flag)
%% (ELEMENT LOOP, GAUSS LOOP, PSI_İ * N_İ , PSI_İ * dN_İ, 2*SUM(F_İ *PSİ_İ^2))
density = zeros(E*GQ,LLQ);
if exc_flag == 2
    grad_density = zeros(E*GQ,LLQ,3);
else
    grad_density = [];
end
for e = 1:E
    for j=1:GQ
        for k = 1:length(fractional_occupancy)
            if exc_flag == 2
                new_d_psi = zeros(1,L_cutoff^2);
            end
            F_s = 0;
            r = 0;
            for i=1:NumElmNods
                F_s = F_s + value_diff_PHI(i,j) * coord(connectivity(e,i));
                r = r + value_PHI(i,j) * coord(connectivity(e,i));
            end
            r_on_ip((e-1)*GQ+j) = r;
            new_psi = zeros(1,L_cutoff^2);
            for orb = 1:L_cutoff^2
                Clm_so = Clm((orb-1)*N+1:orb*N,k);
                for i=1:NumElmNods
                    new_psi(orb)= new_psi(orb) + Clm_so(connectivity(e,i)) * value_PHI(i,j);
                    if exc_flag == 2
                        derPhi = value_diff_PHI(i,j) / F_s ;
                        new_d_psi(orb) = new_d_psi(orb) + Clm_so(connectivity(e,i)) * derPhi; 
                    end
                end
            end
            for l = 1:LLQ
                %% rho = 2*fi*(psi1*Yl1m1*psi1*Yl1m1)
                %% d_rho = 2*fi*2*(d_psi*Ylm+psi*d_Ylm)*psi*Ylm
                new_psi_llq = 0;
                new_grad_psi_llq= zeros(1,3);
                for orb = 1:L_cutoff^2
                    ylm = evalSphericalHarmonic_hb(orb,0,r*LLQ_points(l,1),r*LLQ_points(l,2),r*LLQ_points(l,3),"counts","wo_diff","cartesian");
                    new_psi_llq = new_psi_llq + ylm * new_psi(orb);
                    if exc_flag == 2
                        [~,dYlm_dx, dYlm_dy, dYlm_dz] = evalSphericalHarmonic_hb(orb,0,r*LLQ_points(l,1),r*LLQ_points(l,2),r*LLQ_points(l,3),"counts","w_diff","cartesian");
                        new_grad_psi_llq(1) =  new_grad_psi_llq(1) + (ylm * LLQ_points(l,1)*new_d_psi(orb) + dYlm_dx * new_psi(orb));
                        new_grad_psi_llq(2) =  new_grad_psi_llq(2) + (ylm * LLQ_points(l,2)*new_d_psi(orb) + dYlm_dy * new_psi(orb));
                        new_grad_psi_llq(3) =  new_grad_psi_llq(3) + (ylm * LLQ_points(l,3)*new_d_psi(orb) + dYlm_dz * new_psi(orb));                                     
                    end
                end
                density(GQ*(e-1)+j,l) = density(GQ*(e-1)+j,l) + 2 * fractional_occupancy(k) * new_psi_llq^2; 
                if exc_flag == 2
                    grad_density(GQ*(e-1)+j,l,1) = grad_density(GQ*(e-1)+j,l,1) + 2*fractional_occupancy(k)*2*new_psi_llq * new_grad_psi_llq(1);
                    grad_density(GQ*(e-1)+j,l,2) = grad_density(GQ*(e-1)+j,l,2) + 2*fractional_occupancy(k)*2*new_psi_llq * new_grad_psi_llq(2);
                    grad_density(GQ*(e-1)+j,l,3) = grad_density(GQ*(e-1)+j,l,3) + 2*fractional_occupancy(k)*2*new_psi_llq * new_grad_psi_llq(3);
                end                
            end % end of LLQ loop
        end % end of fractional occupancy loop
    end % end of GQ loop
end % end of element loop
