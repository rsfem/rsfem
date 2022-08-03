% Copyright (C) 2022 M. A. Yalçın, İ. Temizer
% This file is distributed under the terms of the GNU General Public License.
% See the file COPYING for license details.

function [c]= MixingConstantCalculator(n,rho_out,rho_in,E,GQ,  GQ_weights,value_PHI, value_diff_PHI,connectivity,coord,NumElmNods,LLQ,LLQ_points,LLQ_weights)
% warning('off','all')
n_minus_1 = n - 1;
F_n= (rho_out{n}-rho_in{n}); 
F = [];
K = [];
for m = 1:(n_minus_1)
    F_n_minus_m = (rho_out{n-m}-rho_in{n-m}); 
    for k = 1:(n_minus_1)
        F_n_minus_k = (rho_out{n-k}-rho_in{n-k}); 
        Fnnm = (F_n- F_n_minus_m);
        Fnnk = (F_n - F_n_minus_k);    
        K(m,k) = integrator_two_matrices(Fnnm,Fnnk,E,GQ,GQ_weights,value_PHI, value_diff_PHI,connectivity,coord,NumElmNods,LLQ,LLQ_weights,LLQ_points);
    end
    F(m,1) = integrator_two_matrices(Fnnm,F_n,E,GQ,GQ_weights,value_PHI, value_diff_PHI,connectivity,coord,NumElmNods,LLQ,LLQ_weights,LLQ_points);
end

c = K\F;

