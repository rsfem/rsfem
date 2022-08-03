% Copyright (C) 2022 M. A. Yalçın, İ. Temizer
% This file is distributed under the terms of the GNU General Public License.
% See the file COPYING for license details.

function [Res Eps density r_on_ip thetas phis] = HYDROGEN_STARK_FUNCTION(p_l, E_l,L_cutoff_l,potential_multiplier_l,density_plot_flag,alpha)
tic
density = {};
r_on_ip = {};
thetas = {};
phis = {};
for m_ind = 1:length(p_l)
fprintf(" Test =\n\n       %d\n\n ", m_ind)
%% Modifiable Variables
E           = E_l(m_ind)    ;      
p           = p_l(m_ind)    ;      
GQ          = p+3           ;      
n           = 2             ;      
L           = 30            ;      
LLQ         = 110           ;      

L_cutoff    = L_cutoff_l(m_ind);
L_potential = 2;

r1          = 10    ;   
r2          = r1+9  ;   

potential_multiplier = potential_multiplier_l(m_ind);
potential_orbital = 3; 
                  
%% Other Stored Variables
DCoefficients =  zeros(1,L_potential^2);
DCoefficients(1,potential_orbital) = potential_multiplier ;
N   = E * p + 1;        % Number of nodes
NumElmNods = p + 1;     % Number of element nodes
x_min = 0;  x_max = L;  % Boundaries


%% element size calculation
y_coord=zeros(1,E+1);
x_coord=zeros(1,E+1);
y_coord(1,1)=0;
for i=1:E
    y_coord(1,i+1)=y_coord(1,i)+1/E;  
end
x_coord(1,1)=0;
for i=1:E
    x_coord(1,i+1)=L*(y_coord(1,i+1))^n; 
end
for i=1:E
    he_lg(i)=x_coord(1,i+1)-x_coord(1,i);
end

%% Coordinate System Configuration
coord=[];
coord(1)=x_min;
no = 1;
for j_c = 1:E
    for k_c = 1:p
      no = no + 1;
      coord(no) = coord(no-1) + he_lg(j_c)/p;
    end
end

%% Connectivity Array
connectivity = zeros(E,NumElmNods);
counter = 1;
for i = 1:E
    for j = 1:NumElmNods
        connectivity(i,j) = counter;
        counter = counter + 1;
    end
    counter = counter - 1;
end

%% Lebedev Laikov Quadrature
[LLQ_points LLQ_weights] = LebedevLaikovQuadrature(LLQ);
%% Lagrange Element Polynomials
[GQ_points,  GQ_weights] = GaussianQuadrature(GQ); % g, w
%% Shape Function Calculations
[value_PHI value_diff_PHI] = ShapeFunctionandDiff(p,GQ_points);  


%% Generalized Eigenvalue Problem Matrices
M = zeros(N,N); 
A = zeros(N,N);
B = zeros(N,N);
C = zeros(N,N);
D = zeros(N,N);

A_blk=[]; 
B_blk=[];
C_blk=[];
M_blk=[];
H_blk=[];
D_blk=[];


for element = 1:E
    % Initializing the element level matrices for each element
    H_e = zeros(NumElmNods,NumElmNods);
    A_e = zeros(NumElmNods,NumElmNods);
    B_e = zeros(NumElmNods,NumElmNods);
    C_e = zeros(NumElmNods,NumElmNods);
    M_e = zeros(NumElmNods,NumElmNods);  
    D_e = zeros(NumElmNods,NumElmNods);
   
   
    for k = 1:GQ            
        % Calculate x and Jacobian
        x = 0 ;
        F_s = 0;
        for i = 1:NumElmNods
          x = x + value_PHI(i,k) * coord(connectivity(element,i));
          F_s = F_s + value_diff_PHI(i,k) * coord(connectivity(element,i));
        end
        
        % Coefficients
        r = x;
        value_a = 1/2;
        value_b = -1/r;
        value_c = 1/(2*r^2);     
        
        % mollifier
        y = (r^2-r1^2)/(r2^2-r1^2);
        m = 1;
        if alpha == 1
            if r<r1
                m = 1;
            elseif r>r2
                m=0;
            else
                m = 1 + y^2 * (2*y-3);
            end
        end
        value_d = -(sqrt(4*pi/3)*r^alpha)*m;    
        
        
        for i = 1:NumElmNods
            for j = 1:NumElmNods
                A_e(i,j) = A_e(i,j) + (value_a * 1/F_s * (value_diff_PHI(i,k) * value_diff_PHI(j,k))) * GQ_weights(k)*r^2;
                B_e(i,j) = B_e(i,j) + (value_b * F_s   * (value_PHI(i,k)      * value_PHI(j,k))     ) * GQ_weights(k)*r^2;
                C_e(i,j) = C_e(i,j) + (value_c * F_s   * (value_PHI(i,k)      * value_PHI(j,k))     ) * GQ_weights(k)*r^2;
                D_e(i,j) = D_e(i,j) + (value_d * F_s   * (value_PHI(i,k)      * value_PHI(j,k))     ) * GQ_weights(k)*r^2;
                
                M_e(i,j) = M_e(i,j) + (F_s * (value_PHI(i,k) * value_PHI(j,k))) * GQ_weights(k)*r^2;
            end
        end
    end
%% Assembly Process
        for j1 = 1:NumElmNods
            gn1 =        connectivity(element,j1);    
            for j2 = 1:NumElmNods
                gn2 =        connectivity(element,j2);   
                A(gn1, gn2) = A(gn1, gn2) + A_e(j1, j2);
                B(gn1, gn2) = B(gn1, gn2) + B_e(j1, j2);
                C(gn1, gn2) = C(gn1, gn2) + C_e(j1, j2);
                D(gn1, gn2) = D(gn1, gn2) + D_e(j1, j2);
                
                M(gn1, gn2) = M(gn1, gn2) + M_e(j1, j2);
            end
        end    
end % element loop 

%%
Lq = [0 1 2 3 4 5 6 7 8 9 10]; % l orbital quantum number
cnt_s = 0; % counter for sparse matrix
C_dep = C;

for il = 1:L_cutoff % loop  for orbital quantum number           
    C = C_dep*Lq(il)*(Lq(il)+1);
    mq=(-Lq(il):Lq(il));
    for iL=1:length(mq) % loop for magnetic quantum number
        cnt_s = cnt_s + 1;
        A_blk=blkdiag(A_blk,A);
        B_blk=blkdiag(B_blk,B);
        C_blk=blkdiag(C_blk,C);
        
        M_blk=blkdiag(M_blk,M);
    end
    H_blk = A_blk + B_blk + C_blk;
end 

%% D_blk matrix formation
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
            end % end of m2 loop
        end % end of m1 loop
        D_blk(Lq(il1)^2*N+1:(Lq(il1)^2+mql1)*N  , Lq(il2)^2*N+1:(Lq(il2)^2+mql2)*N) = D_macro;
    end % end of l2 loop 
end % end of l1 loop

%% Eigenvalue Problem Solution
H_blk = H_blk + 1/2*(D_blk+D_blk');
[Clm,Eps]=eigs(H_blk,M_blk,N*L_cutoff^2,-0.503771);      

%% Post proccessing
fprintf("Eigenvalue =\n\n%21.12f\n\n ", Eps(1,1))
Eps_stc=nonzeros(Eps);
Eps_sort = sort(Eps_stc);

%% Density calculation for plotting purposes
if density_plot_flag == 1
    Clm_appended = [Eps_stc transpose(Clm)];
    Clm_appended_sorted = transpose(sortrows(Clm_appended,1));
    Clm_sorted = Clm_appended_sorted(2:end,:);
    [density{m_ind} r_on_ip{m_ind} thetas{m_ind} phis{m_ind}] = densityWithoutDFT(E,N,GQ,value_PHI,value_diff_PHI,coord,Clm_sorted,connectivity,NumElmNods,LLQ,L_cutoff);
    save('data_for_plot.mat');
end

%% DElivery of outputs
Res(m_ind,1) = Eps_sort(1,1);
Res(m_ind,2) = N;
Res(m_ind,3) = p;
Res(m_ind,4) = E;
Res(m_ind,5) = L_cutoff;
Res(m_ind,6) = r1;
Res(m_ind,7) = r2;
Res(m_ind,9) = potential_multiplier;
T1 = toc;
fprintf("Time =\n\n%21.12f\n\n ", T1)
Res(m_ind,8) = T1;
end
