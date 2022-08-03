% Copyright (C) 2022 M. A. Yalçın, İ. Temizer
% This file is distributed under the terms of the GNU General Public License.
% See the file COPYING for license details.

function [exc,vxc,wxc]=libxcBridgeGeneral(rhos,sigmas,exc_flag,pp_flag)
if exc_flag == 1 % LDA
    if pp_flag == 0
        mex('-R2018a','-lxc','lda_slater_x.c','-silent'); 
        mex('-R2018a','-lxc','lda_vwn_c.c','-silent'); 
        [ec,vc]=lda_slater_x(rhos);
        [ex,vx]=lda_vwn_c(rhos);
        exc = ec + ex;
        vxc = vc + vx;
        wxc = [];
    else
        mex('-R2018a','-lxc','lda_teter93_xc.c','-silent'); 
        [exc,vxc]=lda_teter93_xc(rhos);
        wxc = [];
    end       
elseif exc_flag == 2 % GGA
    mex('-R2018a','-lxc','gga_pbe_c.c','-silent'); 
    mex('-R2018a','-lxc','gga_pbe_x.c','-silent'); 
    [ec,vc,wc]=gga_pbe_c(rhos,sigmas);
    [ex,vx,wx]=gga_pbe_x(rhos,sigmas);
    exc = ec + ex;
    vxc = vc + vx;
    wxc = wc + wx;
end


