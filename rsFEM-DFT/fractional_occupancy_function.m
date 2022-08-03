% Copyright (C) 2022 M. A. Yalçın, İ. Temizer
% This file is distributed under the terms of the GNU General Public License.
% See the file COPYING for license details.

function [occ eps_f ent] = fractional_occupancy_function(Z,Zion,eps,pp_flag)
eps = nonzeros(eps);
max_iter = 1000;
k_B = 3.1668114e-6;
%temp = 100;
temp = 1e-3/k_B;
sigma = k_B*temp;
numeig0 = 20; 
tol = 1e-9;

%% min
k0 = 0;

%% max
if numeig0 > length(eps)
    numeig0=length(eps);
end
if pp_flag == 1 
   numsap0 = Zion/2;
else
    numsap0 = Z/2;
end
eps_min = 0;
eps_max = 0;
eps_min = min(eps((k0+1):(k0+numeig0)))-1; % res <0
eps_max = max(eps((k0+1):(k0+numeig0)))+1; % res >0

%% minimum residual
res = 0;
for k = (k0+1):(k0+numeig0)
val = (eps(k)-eps_min) / sigma;
val = exp(val);
occ(k) = 1/(1+val);
res = res + occ(k);
end
res = res -numsap0;
if res > 0
    error("res_min>0")    
end

%% maximum residual
res = 0;
for k = (k0+1):(k0+numeig0)
    aa = eps(k);
val = (eps(k)-eps_max) / sigma;
val = exp(val);
occ(k) = 1/(1+val);
res = res + occ(k);
end
res = res -numsap0;
if res < 0
    error("res_max<0")    
end

%% bisection
iter= 1;
res = 1;
while iter<=max_iter 
    eps_f = (eps_min+eps_max)/2;
    res = 0;
    for k = (k0+1):(k0+numeig0)
        val = (eps(k)-eps_f) / sigma;
        val = exp(val);
        occ(k) = 1/(1+val);
        res = res + occ(k);
    end
    res = res -numsap0;
    if abs(res) < tol
        break
    end
    if res < 0 
        eps_min = eps_f;
    else
        eps_max = eps_f;
    end
    iter = iter + 1;
end

%% enthropy
occ(occ < 1e-30)=0;
occ = nonzeros(occ);
numeig0 = length(occ);
ent = 0;
for k =1:numeig0
    occ0=occ(k);
    if occ(k)<tol
        occ0 = tol;
    end
    if occ(k)> (1-tol)
        occ0 = 1-tol;
    end
    ent = ent - sigma * (occ0*log(occ0) + (1-occ0)*log(1-occ0));
end
ent = 2*ent;

