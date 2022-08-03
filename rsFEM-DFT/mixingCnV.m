% Copyright (C) 2022 M. A. Yalçın, İ. Temizer
% This file is distributed under the terms of the GNU General Public License.
% See the file COPYING for license details.

function [rho_in_mixed grad_rho_mixed]= mixingCnV(n_plus_1,rho_in,rho_out,grad_rho_in,grad_rho_out,E,GQ,  GQ_weights,value_PHI, value_diff_PHI,connectivity,coord,NumElmNods,LLQ,LLQ_points,LLQ_weights,NLCC_flag,exc_flag,pp_flag)
gamma = 0.5;
if n_plus_1 == 2 % sc_count = 1
    rho_in_mixed = gamma*rho_in{1} + gamma*rho_out{1}; 
    if exc_flag == 2 
        grad_rho_mixed = gamma*grad_rho_in{1} + gamma*grad_rho_out{1}; 
    else
        grad_rho_mixed = [];
    end
    cnn = 1;
else
n = n_plus_1 - 1;
n_minus_1 = n_plus_1 - 2;
rho_out_n = rho_out{n};
rho_out_n_1 = rho_out{n_minus_1};
c_n = MixingConstantCalculator(n,rho_out(1,1:n),rho_in(1,1:n),E,GQ,  GQ_weights,value_PHI, value_diff_PHI,connectivity,coord,NumElmNods,LLQ,LLQ_points,LLQ_weights);
c_n = [c_n; 1-sum(c_n)];
rho_hat_in_rem = 0;
rho_hat_out_rem = 0;
if exc_flag == 2 
    grad_rho_hat_in_rem = 0;
    grad_rho_hat_out_rem = 0;
end
for k = 1:(n-1)
    rho_hat_in_rem = rho_hat_in_rem + c_n(k)*(rho_in{n-k}); % ck can take any index ?
    rho_hat_out_rem = rho_hat_out_rem + c_n(k)*(rho_out{n-k});
    if exc_flag == 2 
        grad_rho_hat_in_rem = grad_rho_hat_in_rem + c_n(k)*(grad_rho_in{n-k});
        grad_rho_hat_out_rem = grad_rho_hat_out_rem + c_n(k)*(grad_rho_out{n-k});
    end
end
rho_hat_in = c_n(n) * rho_in{n} + rho_hat_in_rem;
rho_hat_out = c_n(n) * rho_out{n} + rho_hat_out_rem;
rho_in_mixed = gamma*rho_hat_out + (1 - gamma)*rho_hat_in;
rho_in_mixed(rho_in_mixed <= 0) = 1e-24;
if exc_flag == 2 
    grad_rho_hat_in = c_n(n) .* grad_rho_in{n} + grad_rho_hat_in_rem;
    grad_rho_hat_out = c_n(n) .* grad_rho_out{n} + grad_rho_hat_out_rem;
    grad_rho_mixed = gamma*grad_rho_hat_out + (1 - gamma)*grad_rho_hat_in;
else
    grad_rho_mixed = [];
end
cnn = c_n;
end
