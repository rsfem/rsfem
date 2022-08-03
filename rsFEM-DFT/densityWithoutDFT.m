% Copyright (C) 2022 M. A. Yalçın, İ. Temizer
% This file is distributed under the terms of the GNU General Public License.
% See the file COPYING for license details.

function [density r_on_ip thetas phis] = densityWithoutDFT(E,N,GQ,value_PHI,value_diff_PHI,coord,Clm,connectivity,NumElmNods,LLQ,L_cutoff,fractional_occupancy)
%% (ELEMENT LOOP, GAUSS LOOP, PSI_İ * N_İ , PSI_İ * dN_İ, 2*SUM(F_İ *PSİ_İ^2))
phis = [pi/2*ones(1,(N+1)/2) 3*pi/2*ones(1,(N-1)/2)];
theta_first_half = [linspace(0,pi,(N+1)/2)];
theta_Second_half = [linspace(-pi,0,(N+1)/2)];
thetas = [theta_first_half theta_Second_half(2:end)];
density = zeros(E*GQ,length(thetas));
for e = 1:E
    for j=1:GQ
        for k = 1:length(fractional_occupancy)
            r = 0;
            for i=1:NumElmNods
                r = r + value_PHI(i,j) * coord(connectivity(e,i));
            end
            r_on_ip((e-1)*GQ+j) = r;
            new_psi = zeros(1,L_cutoff^2);
            for orb = 1:L_cutoff^2
                Clm_so = Clm((orb-1)*N+1:orb*N,k);
                for i=1:NumElmNods
                    new_psi(orb)= new_psi(orb) + Clm_so(connectivity(e,i)) * value_PHI(i,j);
                end
            end
            for tet_i = 1:length(thetas)
                new_psi_llq = 0;
                for orb = 1:L_cutoff^2
                    ylm = evalSphericalHarmonic_hb(orb,0,r,thetas(tet_i),phis(tet_i),"counts","wo_diff","spherical");
                    new_psi_llq = new_psi_llq + ylm * new_psi(orb);
                end
                density(GQ*(e-1)+j,tet_i) = density(GQ*(e-1)+j,tet_i) + 2 * fractional_occupancy(k) * new_psi_llq^2;    
            end
        end % end of fractional occupancy loop
    end % end of GQ loop
end % end of element loop
