% Copyright (C) 2022 M. A. Yalçın, İ. Temizer
% This file is distributed under the terms of the GNU General Public License.
% See the file COPYING for license details.

function [h_cell] = nlParameter(ppP,L_cutoff_nl,i_cutoff,Z,NLCC_flag,exc_flag,pp_flag)
%% ppP
%r_local    C1      C2      C3      C4
%r_0        h^0_11  h^0_22  h^0_33  0
%           0       0       0       0            
%r_1        h^1_11  h^1_22  h^1_33  0
%	        k^1_11  k^1_22  k^1_33  0
%r_2        h^2_11  h^2_22  h^2_33  0
%	        k^2_11  k^2_22  k^2_33  0

%% h^l_i,j
h_cell = cell(L_cutoff_nl+1,1);
h = zeros(i_cutoff);

%% h_cell
%h^0_11 h^0_12 h^0_13 
%h^0_21 h^0_22 h^0_23  
%h^0_31 h^0_32 h^0_33  

%%
for ic = 1:(L_cutoff_nl+1)
    l = ic - 1;
    h(1,1) = ppP((l+1)*2,2);
    h(2,2) = ppP((l+1)*2,3);
    h(3,3) = ppP((l+1)*2,4);
    if exc_flag == 2 || NLCC_flag == 1 %This condition is for GGA strictly       
        h(1,2) = ppP((l+1)*2,5);
        h(2,1) = h(1,2); 
        h(1,3) = ppP((l+1)*2,6);
        h(3,1) = h(1,3); 
        h(2,3) = ppP((l+1)*2,7);
        h(3,2) = h(2,3);               
    else % LDA without NLCC
        h(1,2) = -1/2*sqrt(((l+1)*2+1)/((l+1)*2+3))*h(2,2);
        h(2,1) = h(1,2);
        h(1,3) = 1/2*sqrt(((l+1)*2+3)*((l+1)*2+5/((l+1)*2+7)/((l+1)*2+9)))*h(3,3);
        h(3,1) = h(1,3);
        h(2,3) = -1/2*2*((l+1)*2+5)/sqrt(((l+1)*2+7)*((l+1)*2+9))*h(3,3);
        h(3,2) = h(2,3);
    end
    h_cell{ic} = h; 
end
end
