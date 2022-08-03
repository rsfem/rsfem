% Copyright (C) 2022 M. A. Yalçın, İ. Temizer
% This file is distributed under the terms of the GNU General Public License.
% See the file COPYING for license details.

function [vxc,exc,wxc_vector] = exchangeCorrelation(rho_v,rho_c,grad_rho_v,grad_rho_c,exc_flag,pp_flag,NLCC_flag)
sz = size(rho_v);
if NLCC_flag == 0
    rho_c = zeros(sz(1),sz(2));
    if exc_flag == 2
        grad_rho_c = zeros(sz(1),sz(2),3);
    end
end
rho= rho_v + rho_c;
rho = reshape(rho,sz(1)*sz(2),1);
if exc_flag == 2
    grad_rho = grad_rho_v + grad_rho_c;
    sigma = (grad_rho(:,:,1).*grad_rho(:,:,1)+grad_rho(:,:,2).*grad_rho(:,:,2)+grad_rho(:,:,3).*grad_rho(:,:,3)); 
    wxc_vector = zeros(sz(1),sz(2),3);
    sigma_store=sigma;
    sigma = reshape(sigma,sz(1)*sz(2),1);
    [exc,vxc,wxc]=libxcBridgeGeneral(rho,sigma,exc_flag,pp_flag);
    wxc = reshape(wxc,sz(1),sz(2));
    wxc_vector(:,:,1) =  2*wxc.*grad_rho(:,:,1);
    wxc_vector(:,:,2) =  2*wxc.*grad_rho(:,:,2);
    wxc_vector(:,:,3) =  2*wxc.*grad_rho(:,:,3);
else
    sigma = [];
    [exc,vxc,~]=libxcBridgeGeneral(rho,sigma,exc_flag,pp_flag);
    wxc_vector = [];
end
exc = reshape(exc,sz(1),sz(2));
vxc = reshape(vxc,sz(1),sz(2));
end

