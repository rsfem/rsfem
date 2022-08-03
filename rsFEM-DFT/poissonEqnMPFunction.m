% Copyright (C) 2022 M. A. Yalçın, İ. Temizer
% This file is distributed under the terms of the GNU General Public License.
% See the file COPYING for license details.

function [Vh Vh00] = poissonEqnMPFunction(connectivity,coord,E,GQ,GQ_weights,L,L_poisson,N,NumElmNods,p,value_diff_PHI,value_PHI,Z,Zion,LLQ,LLQ_points,LLQ_weights,flag,density,density00,NLCC_flag,exc_flag,pp_flag)
%flag == 0 righthand side is rho. solution is vh.
%flag == 3 vh is given twice for two different rho matrices.
if pp_flag == 0
    Z_in = Z;
else
    Z_in = Zion;
end
Lq = [0 1 2 3 4 5 6 7 8 9];
Dlm2 = [];
Dlm1 = [];
Dlm = [];

%% element loop for left-hand side of poisson problem
Dlm_micro1 = zeros(N,N);
Dlm_micro2 = zeros(N,N);
for element = 1:E 
    Dlm_micro_e1 = zeros(NumElmNods,NumElmNods);
    Dlm_micro_e2 = zeros(NumElmNods,NumElmNods);
    for k = 1:GQ           
        r = 0 ;
        F_s_r = 0;
        for i = 1:NumElmNods
            r = r + value_PHI(i,k) * coord(connectivity(element,i));
            F_s_r = F_s_r + value_diff_PHI(i,k) * coord(connectivity(element,i));
        end
        value_lhs1 = 1/(4*pi);
        value_lhs2 = 1/r^2/(4*pi); 
        for j = 1:NumElmNods
            for i = 1:NumElmNods
                Dlm_micro_e1(i,j)= Dlm_micro_e1(i,j) + value_lhs1 * 1/F_s_r * value_diff_PHI(j,k) * value_diff_PHI(i,k) * GQ_weights(k) * r^2 ; 
                Dlm_micro_e2(i,j)= Dlm_micro_e2(i,j) + value_lhs2 * F_s_r * value_PHI(j,k) * value_PHI(i,k) * GQ_weights(k) * r^2 ;
            end
        end
    end        
% Assembly Process
    for j1 = 1:NumElmNods
        gn1 = connectivity(element,j1);                               
        for j2 = 1:NumElmNods
            gn2 = connectivity(element, j2);
            Dlm_micro1(gn1, gn2) = Dlm_micro1(gn1, gn2) + Dlm_micro_e1(j1, j2);                       
            Dlm_micro2(gn1, gn2) = Dlm_micro2(gn1, gn2) + Dlm_micro_e2(j1, j2);  
        end
    end
end % element loop ends
Dlm_store = Dlm_micro2;
for in = 1:L_poisson      
    il = Lq(in); 
    Dlm_micro2 = Dlm_store*(il)*(il+1);
    for bl=1:(2*Lq(in)+1)
        Dlm1 = blkdiag(Dlm1 ,Dlm_micro1);
        Dlm2 = blkdiag(Dlm2 ,Dlm_micro2);
    end
end % Lq loop ends
Dlm = Dlm1 + Dlm2;

if flag == 0
%% element loop for right-hand side of poisson problem
Clm_micro = zeros(N,L_poisson^2);
for element = 1:E
    Clm_micro_e = zeros(NumElmNods,L_poisson^2);
    for k = 1:GQ           
        r = 0 ;
        F_s_r = 0;
        for i = 1:NumElmNods
            r = r + value_PHI(i,k) * coord(connectivity(element,i));
            F_s_r = F_s_r + value_diff_PHI(i,k) * coord(connectivity(element,i));
        end
        for orb= 1:L_poisson^2
            Clm_micro_e_orb = zeros(NumElmNods,1);            
            for l = 1:LLQ
                density_one_sphere = density(:,l);
                sh = evalSphericalHarmonic_hb(orb,0,r*LLQ_points(l,1),r*LLQ_points(l,2),r*LLQ_points(l,3),"counts","wo_diff","cartesian");
                value_rhs = density_one_sphere(GQ*(element-1)+k)*sh;                
                for j = 1:NumElmNods
                    Clm_micro_e_orb(j)= Clm_micro_e_orb(j) + value_rhs * F_s_r *value_PHI(j,k) * GQ_weights(k)*LLQ_weights(l)*4*pi*r^2; % 
                end  
            end % LLQ loop ends 
            Clm_micro_e(1:j,orb) = Clm_micro_e(1:j,orb) + Clm_micro_e_orb;
        end % L_cutoff loop ends
    end % GQ loop ends  
    for j1 = 1:NumElmNods
        gn1 = connectivity(element,j1);                               
        for orb = 1:L_poisson^2
            Clm_micro(gn1, orb) = Clm_micro(gn1, orb) + Clm_micro_e(j1, orb);   
        end
    end
end % element loop ends
Clm = Clm_micro;
Clm_bc =Clm(:,1);
Dlm_bc = Dlm(1:N,1:N);

%% boundary conditions 
Clm_bc(:) = Clm_bc(:) - Dlm_bc(:,N) * (+Z_in/L*(sqrt(4*pi)));   

Dlm_bc(N,:) = [];        Clm_bc(N) = [];
Dlm_bc(:,N) = [];

Vh_1s = Dlm_bc\Clm_bc;
Vh_1s(N) = +Z_in/L*(sqrt(4*pi));

Vh = [];
if L_poisson >= 2
    for i = 2:L_poisson^2
        Clm_no_bc =Clm(:,i);
        Dlm_no_bc = Dlm((i-1)*N+1:i*N,(i-1)*N+1:i*N);
        Dlm_no_bc(N,:) = [];        Clm_no_bc(N) = [];
        Dlm_no_bc(:,N) = [];
        Vh_after_1s = Dlm_no_bc\Clm_no_bc;
        Vh_after_1s(N) = 0;
        Vh = [Vh;Vh_after_1s];
    end
    Vh = [Vh_1s ; Vh];
else
    Vh = Vh_1s;
end
Vh00 = [];

elseif flag == 3
%% element loop for right-hand side of poisson problem
Clm_micro = zeros(N,L_poisson^2);
Clm_micro00 = zeros(N,L_poisson^2);
for element = 1:E % 
    Clm_micro_e = zeros(NumElmNods,L_poisson^2);
    Clm_micro_e00 = zeros(NumElmNods,L_poisson^2);
    for k = 1:GQ           
        r = 0 ;
        F_s_r = 0;
        for i = 1:NumElmNods
            r = r + value_PHI(i,k) * coord(connectivity(element,i));
            F_s_r = F_s_r + value_diff_PHI(i,k) * coord(connectivity(element,i));
        end
        for orb= 1:L_poisson^2
            Clm_micro_e_orb = zeros(NumElmNods,1);
            Clm_micro_e_orb00 = zeros(NumElmNods,1);     
            for l = 1:LLQ
                density_one_sphere = density(:,l);
                density_one_sphere00 = density00(:,l);
                sh = evalSphericalHarmonic_hb(orb,0,r*LLQ_points(l,1),r*LLQ_points(l,2),r*LLQ_points(l,3),"counts","wo_diff","cartesian");
                value_rhs = density_one_sphere(GQ*(element-1)+k)*sh;
                value_rhs00 = density_one_sphere00(GQ*(element-1)+k)*sh;                  
                for j = 1:NumElmNods
                    Clm_micro_e_orb(j)= Clm_micro_e_orb(j) + value_rhs * F_s_r *value_PHI(j,k) * GQ_weights(k)*LLQ_weights(l)*4*pi*r^2; 
                    Clm_micro_e_orb00(j)= Clm_micro_e_orb00(j) + value_rhs00 * F_s_r *value_PHI(j,k) * GQ_weights(k)*LLQ_weights(l)*4*pi*r^2;
                end  
            end % LLQ loop ends 
            Clm_micro_e(1:j,orb) = Clm_micro_e(1:j,orb) + Clm_micro_e_orb;
            Clm_micro_e00(1:j,orb) = Clm_micro_e00(1:j,orb) + Clm_micro_e_orb00;
        end % L_cutoff loop ends 
    end % GQ loop ends     
    for j1 = 1:NumElmNods
        gn1 = connectivity(element,j1);                               
        for orb = 1:L_poisson^2
            Clm_micro(gn1, orb) = Clm_micro(gn1, orb) + Clm_micro_e(j1, orb);   
            Clm_micro00(gn1, orb) = Clm_micro00(gn1, orb) + Clm_micro_e00(j1, orb); 
        end
    end
end % element loop ends
Clm = Clm_micro;
Clm00 = Clm_micro00;
Clm_bc =Clm(:,1);
Clm_bc00 =Clm00(:,1);
Dlm_bc = Dlm(1:N,1:N);

%% boundary conditions
Clm_bc(:) = Clm_bc(:) - Dlm_bc(:,N) * (+Z_in/L*(sqrt(4*pi)));   
Clm_bc00(:) = Clm_bc00(:) - Dlm_bc(:,N) * (+Z_in/L*(sqrt(4*pi)));  

Dlm_bc(N,:) = [];        Clm_bc(N) = []; Clm_bc00(N) = [];
Dlm_bc(:,N) = [];
Vh_1s = Dlm_bc\Clm_bc;
Vh_1s(N) = +Z_in/L*(sqrt(4*pi));

Vh_1s00 = Dlm_bc\Clm_bc00;
Vh_1s00(N) = +Z_in/L*(sqrt(4*pi));

Vh = [];
Vh00 = [];
if L_poisson >= 2
    for i = 2:L_poisson^2
        Clm_no_bc =Clm(:,i);
        Clm_no_bc00 =Clm00(:,i);
        Dlm_no_bc = Dlm((i-1)*N+1:i*N,(i-1)*N+1:i*N);
        Dlm_no_bc(N,:) = [];        
                            Clm_no_bc(N) = [];
        Dlm_no_bc(:,N) = [];        
                            Clm_no_bc00(N) = [];
        Vh_after_1s = Dlm_no_bc\Clm_no_bc;
        Vh_after_1s00 = Dlm_no_bc\Clm_no_bc00;
        Vh_after_1s(N) = 0;
        Vh_after_1s00(N) = 0;
        Vh = [Vh;Vh_after_1s];
        Vh00 = [Vh00;Vh_after_1s00];
    end
    Vh = [Vh_1s ; Vh];
    Vh00 = [Vh_1s00 ; Vh00];
else 
    Vh = Vh_1s;
    Vh00 = Vh_1s00;
end
  
end

