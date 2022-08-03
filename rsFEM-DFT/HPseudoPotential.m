% Copyright (C) 2022 M. A. Yalçın, İ. Temizer
% This file is distributed under the terms of the GNU General Public License.
% See the file COPYING for license details.

function [Hnl_blk] = HPseudoPotential(E,GQ,connectivity,value_PHI,value_diff_PHI,coord,GQ_weights,Z,L_cutoff,L_potential,Vh_in,rho_in_n,N,NumElmNods,ppP,NLCC_flag,exc_flag,pp_flag)
Lq = [0 1 2 3 4 5 6 7 8 9];
i_cutoff = 3; % investigated orbitals cutoff
sz1 = size(ppP);
for i = 2:sz1-1
    if isempty(nonzeros(ppP(i:i+1,:))) == 1
        Lc_limit = (i-1)/2;
        break
    end
end
if L_cutoff<= Lc_limit
    L_cutoff_nl = L_cutoff; 
else
    L_cutoff_nl = Lc_limit; 
end
if Z == 1
    L_cutoff_nl = 0;
end
%r_local    C1      C2      C3      C4
%r_0        h^0_11  h^0_22  h^0_33
%           0       0       0            
%r_1        h^1_11  h^1_22  h^1_33
%	        k^1_11  k^1_22  k^1_33
%r_2        h^2_11  h^2_22  h^2_33
%	        k^2_11  k^2_22  k^2_33
%L_cutoff_nl = 2;
nlP = nlParameter(ppP,L_cutoff_nl,i_cutoff,Z,NLCC_flag,exc_flag,pp_flag);
Hnl_cell_matrix1  = cell(L_cutoff_nl^2,1); 
for il = 1:L_cutoff_nl
    l = Lq(il);
    im_vector = [-l:l];                       
    for im = 1:length(im_vector)
        m = im_vector(im);
        Hnl_cell = cell(i_cutoff,1); 
        for p_i = 1:i_cutoff
                Hnl_vector = zeros(N,1); 
                for element = 1:E
                    Hnl_e = zeros(NumElmNods,1);
                    for k = 1:GQ 
                        r = 0;
                        F_s = 0;
                        for i = 1:NumElmNods
                            r = r + value_PHI(i,k) * coord(connectivity(element,i));
                            F_s = F_s + value_diff_PHI(i,k) * coord(connectivity(element,i));
                        end
                        [value_p1]= GaussianProjector(r,l,p_i,ppP);
                        for i = 1:NumElmNods 
                        	Hnl_e(i) = Hnl_e(i) + ((value_p1)* F_s   * value_PHI(i,k)  ) * GQ_weights(k)*r^2;%*4*pi;?
                        end        
                    end
                    for j1 = 1:NumElmNods
                        gn1 =	connectivity(element,j1);      
                        Hnl_vector(gn1,1)= Hnl_vector(gn1, 1) + Hnl_e(j1,1);   
                    end
                end % end of element loop
            Hnl_cell{p_i} = Hnl_vector;
        end % end of investigated orbitals loop          
        Hnl_matrix = zeros(N,N); 
        for p_i = 1:i_cutoff
            for p_j = 1:i_cutoff
                [l,p_i,p_j,nlP{l+1}(p_i,p_j)];
               Hnl_matrix = Hnl_matrix + nlP{l+1}(p_i,p_j) * Hnl_cell{p_i} * transpose(Hnl_cell{p_j});
            end
        end
        Hnl_matrix = (Hnl_matrix + transpose(Hnl_matrix))/2;
        Hnl_cell_matrix1{l^2+l+1+m} = Hnl_matrix;  
    end % end of m loop
end % end of l loop
Hnl_blk  = [];
for i = 1:L_cutoff_nl^2
    Hnl_blk = blkdiag(Hnl_blk,Hnl_cell_matrix1{i});
end
if L_cutoff>Lc_limit
    f= L_cutoff^2-L_cutoff_nl^2;
    for i=1:f
        Hnl_blk = blkdiag(Hnl_blk,zeros(N,N));
    end
end     
end

