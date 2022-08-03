% Copyright (C) 2022 M. A. Yalçın, İ. Temizer
% This file is distributed under the terms of the GNU General Public License.
% See the file COPYING for license details.

function [M_blk H_blk variables1] = schrodingerEqnWithVhandVxcCombined(E,GQ,connectivity,value_PHI,value_diff_PHI,coord,GQ_weights,Z,Zion,L_cutoff,L_potential,Vh,rho,rho_c,N,NumElmNods,LLQ,LLQ_weights,LLQ_points,ppP,grad_rho,grad_rho_c,NLCC_flag,exc_flag,pp_flag,external_field_flag,L_poisson,DCoefficients,potential_multiplier)
Lq = [0 1 2 3 4 5 6 7 8 9]; % quantum number
M = zeros(N,N); 
A = zeros(N,N);
B = zeros(N,N);
B1_cell = cell(1,L_poisson^2);
C = zeros(N,N);
D = zeros(N,N);

A_blk  = []; 
B1_blk = [];
B2_blk = [];
D_blk  = [];
M_blk  = [];
H_blk  = [];

for element = 1:E
    A_e = zeros(NumElmNods,NumElmNods);
    B_e = zeros(NumElmNods,NumElmNods);
    C_e = zeros(NumElmNods,NumElmNods);
    M_e = zeros(NumElmNods,NumElmNods);  
    D_e = zeros(NumElmNods,NumElmNods);  
    for k = 1:GQ            
        r = 0 ;
        F_s = 0;
        for i = 1:NumElmNods
          r = r + value_PHI(i,k) * coord(connectivity(element,i));
          F_s = F_s + value_diff_PHI(i,k) * coord(connectivity(element,i));
        end
        value_a = 1/2; % kinetic energy coefficient
        value_c = 1/(2*r^2); %l dependent term   / centrifugal term
        if pp_flag == 1 
            value_d_int = r/ppP(1,1);
            value_d = -potential_multiplier/r*erf(value_d_int/sqrt(2))+exp(-1/2*(value_d_int)^2)*(ppP(1,2)+ppP(1,3)*(value_d_int)^2+ppP(1,4)*(value_d_int)^4+ppP(1,5)*(value_d_int)^6); % external potential            
        
            value_b3_int = r/ppP(1,1);
            value_b3 = -Zion/r*erf(value_b3_int/sqrt(2))+exp(-1/2*(value_b3_int)^2)*(ppP(1,2)+ppP(1,3)*(value_b3_int)^2+ppP(1,4)*(value_b3_int)^4+ppP(1,5)*(value_b3_int)^6); % nucleus-electron potential
        else 
            value_d = -potential_multiplier/r ; % external potential  
            value_b3 = -Z/r; % nucleus-electron potential
        end        
        for i = 1:NumElmNods  
            for j = 1:NumElmNods
                A_e(i,j) = A_e(i,j) + (value_a * 1/F_s * (value_diff_PHI(i,k) * value_diff_PHI(j,k))) * GQ_weights(k)*r^2;
                B_e(i,j) = B_e(i,j) + ((value_b3 )* F_s   * (value_PHI(i,k)      * value_PHI(j,k))     ) * GQ_weights(k)*r^2;                                    
                C_e(i,j) = C_e(i,j) + (value_c * F_s   * (value_PHI(i,k)      * value_PHI(j,k))     ) * GQ_weights(k)*r^2;
                D_e(i,j) = D_e(i,j) + (value_d * F_s   * (value_PHI(i,k)      * value_PHI(j,k))     ) * GQ_weights(k)*r^2;
                                     
                M_e(i,j) = M_e(i,j) + (F_s * (value_PHI(i,k) * value_PHI(j,k))) * GQ_weights(k)*r^2;
             end
        end
    end % end of GQ loop
    for j1 = 1:NumElmNods
        gn1 =        connectivity(element,j1);   
        for j2 = 1:NumElmNods
            gn2  =        connectivity(element,j2);     
            A(gn1, gn2) = A(gn1, gn2) + A_e(j1, j2);                
            B(gn1, gn2) = B(gn1, gn2) + B_e(j1, j2);
            C(gn1, gn2) = C(gn1, gn2) + C_e(j1, j2);
            D(gn1, gn2) = D(gn1, gn2) + D_e(j1, j2);
            M(gn1, gn2) = M(gn1, gn2) + M_e(j1, j2);                                 
        end
    end               
end % end of element loop 

C_dep = C;
for il1 = 1:L_cutoff  
    C = C_dep*Lq(il1)*(Lq(il1)+1);
    mq=(-Lq(il1):Lq(il1));
    for im1=1:length(mq) 
        A_blk=blkdiag(A_blk , A + C + B);
        M_blk=blkdiag(M_blk,M);      
    end % end of m loop
     H_blk = A_blk;
end % end of l loop

%% Vxc
if  exc_flag == 1         
    [vxcMat] = exchangeCorrelation(rho,rho_c,grad_rho,grad_rho_c,exc_flag,pp_flag,NLCC_flag);
elseif exc_flag == 2 
    [vxcMat,~,wxcMat] = exchangeCorrelation(rho,rho_c,grad_rho,grad_rho_c,exc_flag,pp_flag,NLCC_flag);
end
% addressing order  [1 2 3]
%                   [4 5 6]
%                   [7 8 9]   
for il1 = 1:L_cutoff
    l1 = Lq(il1);
    im1_vector = [-l1:l1];
    for im1 = 1:length(im1_vector)
        m1 = im1_vector(im1);
        for il2 = 1:L_cutoff
            l2 = Lq(il2);
            im2_vector = [-l2:l2];                       
            for im2 = 1:length(im2_vector)
                m2 = im2_vector(im2);
                B2 = zeros(N,N);
                if exc_flag == 2 
                    B_gga1 = zeros(N,N);
                    B_gga2 = zeros(N,N);
                end
                for element = 1:E
                    B_e_2 = zeros(NumElmNods,NumElmNods);
                    if exc_flag == 2 
                        B_e_gga1 = zeros(NumElmNods,NumElmNods);
                        B_e_gga2 = zeros(NumElmNods,NumElmNods);
                    end   
                    for k = 1:GQ 
                        if exc_flag == 2 
                            derPhi = zeros(NumElmNods,GQ);
                        end 
                        r = 0;
                        F_s = 0;
                        for i = 1:NumElmNods
                            F_s = F_s + value_diff_PHI(i,k) * coord(connectivity(element,i));
                       end
                       for i = 1:NumElmNods
                            r = r + value_PHI(i,k) * coord(connectivity(element,i));
                            if exc_flag == 2 
                                derPhi(i,j) =  derPhi(i,j) + value_diff_PHI(i,j) / F_s ;
                            end
                       end
                        for l = 1:LLQ
                            sh1 =  evalSphericalHarmonic_hb(l1,m1,r*LLQ_points(l,1),r*LLQ_points(l,2),r*LLQ_points(l,3),"orbs","wo_diff","cartesian");
                            sh2 =  evalSphericalHarmonic_hb(l2,m2,r*LLQ_points(l,1),r*LLQ_points(l,2),r*LLQ_points(l,3),"orbs","wo_diff","cartesian");
                            value_vxc = vxcMat(GQ*(element-1)+k,l)*sh1*sh2;%
                            if exc_flag == 2                    
                                [~, d_Ylm1_d_x, d_Ylm1_d_y, d_Ylm1_d_z] =  evalSphericalHarmonic_hb(l1,m1,r*LLQ_points(l,1),r*LLQ_points(l,2),r*LLQ_points(l,3),"orbs","w_diff","cartesian");
                                [~, d_Ylm2_d_x, d_Ylm2_d_y, d_Ylm2_d_z] =  evalSphericalHarmonic_hb(l2,m2,r*LLQ_points(l,1),r*LLQ_points(l,2),r*LLQ_points(l,3),"orbs","w_diff","cartesian");
                            end
                            for i = 1:NumElmNods %% element calc for B2 only (rho dependent)
                                for j = 1:NumElmNods
                                    B_e_2(i,j) = B_e_2(i,j) + value_vxc * F_s   * (value_PHI(i,k)      * value_PHI(j,k)) * GQ_weights(k)*LLQ_weights(l)*r^2*4*pi;            
                                     if exc_flag == 2    
                                        value_b1 =  wxcMat(GQ*(element-1)+k,l,1) * value_PHI(j,k) * sh1 * (value_diff_PHI(i,k)*LLQ_points(l,1) * sh2 + value_PHI(i,k)*d_Ylm2_d_x*F_s) + ...
                                                    wxcMat(GQ*(element-1)+k,l,2) * value_PHI(j,k) * sh1 * (value_diff_PHI(i,k)*LLQ_points(l,2) * sh2 + value_PHI(i,k)*d_Ylm2_d_y*F_s) + ...
                                                    wxcMat(GQ*(element-1)+k,l,3) * value_PHI(j,k) * sh1 * (value_diff_PHI(i,k)*LLQ_points(l,3) * sh2 + value_PHI(i,k)*d_Ylm2_d_z*F_s);                                                    
                                        value_b2 =  wxcMat(GQ*(element-1)+k,l,1) * value_PHI(i,k) * sh2 * (value_diff_PHI(j,k)*LLQ_points(l,1) * sh1 + value_PHI(j,k)*d_Ylm1_d_x*F_s) + ...
                                                    wxcMat(GQ*(element-1)+k,l,2) * value_PHI(i,k) * sh2 * (value_diff_PHI(j,k)*LLQ_points(l,2) * sh1 + value_PHI(j,k)*d_Ylm1_d_y*F_s) + ...
                                                    wxcMat(GQ*(element-1)+k,l,3) * value_PHI(i,k) * sh2 * (value_diff_PHI(j,k)*LLQ_points(l,3) * sh1 + value_PHI(j,k)*d_Ylm1_d_z*F_s);                                                   
                                        B_e_gga1(i,j) = B_e_gga1(i,j) + value_b1 * GQ_weights(k)*LLQ_weights(l)*r^2*4*pi; 
                                        B_e_gga2(i,j) = B_e_gga2(i,j) + value_b2 * GQ_weights(k)*LLQ_weights(l)*r^2*4*pi;                           
                                    end
                                end
                            end             
                        end % end of LLQ loop
                    end % end of GQ loop
                     for j1 = 1:NumElmNods
                        gn1 =        connectivity(element,j1);    
                        for j2 = 1:NumElmNods
                            gn2  =        connectivity(element,j2);   
                            B2(gn1, gn2)= B2(gn1, gn2) + B_e_2(j1, j2); 
                            if exc_flag == 2 
                                B_gga1(gn1, gn2)= B_gga1(gn1, gn2) + B_e_gga1(j1, j2);
                                B_gga2(gn1, gn2)= B_gga2(gn1, gn2) + B_e_gga2(j1, j2);
                            else
                            B_gga1 = 0;
                            B_gga2 = 0;
                            end
                        end
                    end
                end % end of element loop
                B2_blk((l1^2+(m1+l1+1)-1)*N+1:(l1^2+(m1+l1+1))*N,(l2^2+(m2+l2+1)-1)*N+1:(l2^2+(m2+l2+1))*N) = B2+B_gga1+B_gga2;
            end % end of m2 loop
        end % end of l2 loop
    end % end of m1 loop
end % end of l1 loop

%% Vh 
% % addressing order  [1 2 3]
% %                   [4 6 7]
% %                   [5 8 9]
for orb = 1:L_poisson^2
    B1 = zeros(N,N);
    Vh_orb = Vh((orb-1)*N+1:orb*N);
    for element = 1:E
        B_e_1 = zeros(NumElmNods,NumElmNods);
        for k = 1:GQ 
            r = 0 ;
            F_s = 0;
            new_Vh = 0;
            for i = 1:NumElmNods
                r = r + value_PHI(i,k) * coord(connectivity(element,i));
                F_s = F_s + value_diff_PHI(i,k) * coord(connectivity(element,i));
                new_Vh = new_Vh + value_PHI(i,k) * Vh_orb(connectivity(element,i)); 
            end
                value_b1    = new_Vh;
                for i = 1:NumElmNods
                    for j = 1:NumElmNods
                        B_e_1(i,j) = B_e_1(i,j) + ((value_b1 )* F_s   * (value_PHI(i,k)      * value_PHI(j,k))     )  *r^2* GQ_weights(k);                            
                    end
                end
        end % end of GQ loop
        for j1 = 1:NumElmNods
            gn1 =        connectivity(element,j1);    
            for j2 = 1:NumElmNods
                gn2  =        connectivity(element,j2);   
                B1(gn1, gn2)= B1(gn1, gn2) + B_e_1(j1, j2);         
            end
        end       
    end % end of element loop 
    B1_cell{orb} = B1;     
end % end of orbital loop

%% B1_blk matrix setup 
for il1 = 1:L_cutoff  
    mq1=(-Lq(il1):Lq(il1));
    mql1=length(mq1);
    l1 = Lq(il1);
    for il2=1:L_cutoff 
        l2 = Lq(il2);
        mq2=(-Lq(il2):Lq(il2));
        mql2=length(mq2);
        B1_macro  = zeros(length(mq1)*N,length(mq2)*N);   
        for im1=1:length(mq1) 
            m1 = mq1(im1);
            for im2=1:length(mq2) 
                m2 = mq2(im2);
                CoefficientIndex = 0;
                B1_micro = zeros(N,N); 
                for il3 = 1:L_poisson 
                    l3 = Lq(il3);
                    mq3=(-Lq(il3):Lq(il3));
                    for im3 = 1:length(mq3) 
                        m3 = mq3(im3);
                        CoefficientIndex = CoefficientIndex + 1;
                        q = 0;
                        for l=1:LLQ                  
                            sh1 =  evalSphericalHarmonic_hb(l1,m1,r*LLQ_points(l,1),r*LLQ_points(l,2),r*LLQ_points(l,3),"orbs","wo_diff","cartesian");
                            sh2 =  evalSphericalHarmonic_hb(l2,m2,r*LLQ_points(l,1),r*LLQ_points(l,2),r*LLQ_points(l,3),"orbs","wo_diff","cartesian");
                            sh3 =  evalSphericalHarmonic_hb(l3,m3,r*LLQ_points(l,1),r*LLQ_points(l,2),r*LLQ_points(l,3),"orbs","wo_diff","cartesian");
                            q = q + sh1*sh2*sh3*LLQ_weights(l)*4*pi;
                        end                    
                        B1_micro = B1_micro + B1_cell{CoefficientIndex}*q ;
                    end  % end of m3 loop
                end % end of l3 loop
                B1_macro((im1-1)*N+1:(im1)*N,(im2-1)*N+1:(im2)*N)  =  B1_micro;
            end % end of m2 loop 
        end % end of m1 loop
            B1_blk(Lq(il1)^2*N+1:(Lq(il1)^2+mql1)*N  , Lq(il2)^2*N+1:(Lq(il2)^2+mql2)*N) = B1_macro; 
    end % end of l2 loop
end % end of l1 loop

%% D_blk matrix formation
if external_field_flag == 1
    for il1 = 1:L_cutoff  
        mq1=(-Lq(il1):Lq(il1));
        mql1=length(mq1);
        l1 = Lq(il1);
        for il2=1:L_cutoff 
        l2 = Lq(il2);
            mq2=(-Lq(il2):Lq(il2));
            mql2=length(mq2);
            D_macro  = zeros(length(mq1)*N,length(mq2)*N);
            for im1=1:length(mq1)
                m1 = mq1(im1);
                for im2=1:length(mq2) 
                    m2 = mq2(im2);
                    CoefficientIndexD = 0;
                    D_micro = zeros(N,N); 
                    if external_field_flag == 1
                        for il3 = 1:L_potential
                            l3 = Lq(il3);
                            mq3=(-Lq(il3):Lq(il3));
                            for im3 = 1:length(mq3)  
                                m3 = mq3(im3);
                                CoefficientIndexD = CoefficientIndexD + 1;
                                q = 0;
                                for l=1:LLQ                            
                                    sh1 =  evalSphericalHarmonic_hb(l1,m1,r*LLQ_points(l,1),r*LLQ_points(l,2),r*LLQ_points(l,3),"orbs","wo_diff","cartesian");
                                    sh2 =  evalSphericalHarmonic_hb(l2,m2,r*LLQ_points(l,1),r*LLQ_points(l,2),r*LLQ_points(l,3),"orbs","wo_diff","cartesian");
                                    sh3 =  evalSphericalHarmonic_hb(l3,m3,r*LLQ_points(l,1),r*LLQ_points(l,2),r*LLQ_points(l,3),"orbs","wo_diff","cartesian");
                                    q = q + (sh1*sh2*sh3)*LLQ_weights(l)*4*pi;
                                end
                                D_micro = D_micro + D.*q*DCoefficients(CoefficientIndexD); 
                            end % end of m3 loop   
                        end % end of l3 loop
                        D_macro((im1-1)*N+1:(im1)*N,(im2-1)*N+1:(im2)*N)  =  D_micro;
                    end
                end % end of m2 loop
            end % end of m1 loop
            D_blk(Lq(il1)^2*N+1:(Lq(il1)^2+mql1)*N  , Lq(il2)^2*N+1:(Lq(il2)^2+mql2)*N) = D_macro;
        end % end of l2 loop 
    end % end of l1 loop
end

H_blk = H_blk + (B1_blk + transpose(B1_blk))/2;
H_blk = H_blk + (B2_blk + transpose(B2_blk))/2;
if external_field_flag == 1
    H_blk = H_blk + (D_blk + transpose(D_blk))/2;
end

variables1 = {E,GQ,LLQ,NumElmNods,coord,connectivity,LLQ_points,value_PHI};
