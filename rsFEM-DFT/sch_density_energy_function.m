% Copyright (C) 2022 M. A. Yalçın, İ. Temizer
% This file is distributed under the terms of the GNU General Public License.
% See the file COPYING for license details.

function [Res2 Clm_m density thetas phis]= sch_density_energy_function(E_l,L_cutoff_l,F_l,p_l,mp_l,LLQ_l,conf,density_plot_flag)
density = cell(length(E_l));
r_on_ip = cell(length(E_l));
thetas = cell(length(E_l));
phis = cell(length(E_l));
for m= 1:length(E_l)
    fprintf(" Test =\n\n       %d\n\n ", m)
        tic
    E           = E_l(m);           
    p           = p_l(m);              
    GQ          = p+3;           
    LLQ         = LLQ_l(m);            
    n           = 2 ;            
    L = 30;     % 

    %% Atom properties
    Z =13;      % atomic number
    Zion = 3;   % valence electrons number

    %% orbital parameters
    L_cutoff = L_cutoff_l(m);
    L_poisson = L_cutoff+1;
    L_potential = 4;

    %% configuration
    % exchange-correlation functional: LDA (1) or GGA (2);
    exc_flag =conf(1);
    % pseudopotential flag: off (0) of on (1). 0 means all-electron
    pp_flag = conf(2);
    % nlcc_flag: off (0) of on (1), requires pp_flag = 1;
    NLCC_flag = conf(3);

    flagCheck(exc_flag,pp_flag,NLCC_flag);

    mixing_flag = true;

    if pp_flag == 0
        Zion = Z;
    end

    %% external potential initiation
    % external_field_flag = 0 -> no external potential  
    % external_field_flag = 1 -> external potential in multipole form
    % external_field_flag = 2 -> external potential in analytical form
    F = F_l(m);
    potential_multiplier =F;
    external_field_flag = 1; % to shut off external field effect completely, select 0. 
    DCoefficients = zeros(1,L_potential^2);
    if external_field_flag == 1      
        DCoefficients(1,3) = mp_l(1,m);
        DCoefficients(1,6) = mp_l(2,m);
        DCoefficients(1,16) = mp_l(3,m);        
    end

    %% Stored Variables
    N   = E * p + 1;        % Number of nodes
    NumElmNods = p + 1;     % Number of element nodes

    %% FEM building up
    [coord connectivity]= BuildUp(E,n,p,NumElmNods,L);

    %% Lagrange Element Polynomials
    [GQ_points  GQ_weights] = GaussianQuadrature(GQ); 
    [LLQ_points LLQ_weights] = LebedevLaikovQuadrature(LLQ);

    %% Shape Function Calculations
    [value_PHI value_diff_PHI] = ShapeFunctionandDiff(p,GQ_points); 

    %% rho guess, initialization of rho vector from an exponentially decaying function
    [rho rho_c grad_rho_c] = InitialRhoGuess(E,GQ,connectivity,value_PHI,coord,NumElmNods,  GQ_weights,value_diff_PHI,LLQ,LLQ_points,LLQ_weights,Z,Zion,NLCC_flag,exc_flag,pp_flag);
    rho_in_n = rho;
    rho_in{1} = rho_in_n;
    if NLCC_flag == 0 || (Z == 1)
        rho_c = zeros(E*GQ,LLQ);
        grad_rho_c = zeros(E*GQ,LLQ,3);
    end

    %% poisson solution
    [Vh] = poissonEqnMPFunction(connectivity,coord,E,GQ,GQ_weights,L,L_poisson,N,NumElmNods,p,value_diff_PHI,value_PHI,Z,Zion,LLQ,LLQ_points,LLQ_weights,0,rho_in_n,rho_in_n,NLCC_flag,exc_flag,pp_flag);
    Vh_in = Vh;

    %% convergence tolerances
    acceptable_error_density = 1e-9;
    acceptable_error_energy = 1e-9;
    error_density = 1;
    error_energy = 1;
    sc_count = 1;
    Emer_previous = 0;

    %% pseudo-potential contribution
    if pp_flag == 1
        if exc_flag == 2
            grad_rho_in_n = zeros(E*GQ,LLQ,3);
            grad_rho_in{1} = grad_rho_in_n;
        else
            grad_rho_in = [];
            grad_rho_in_n = [];
        end
        ppP = ppParameters(Z,NLCC_flag,exc_flag,pp_flag);
        [Hnl_blk]= HPseudoPotential(E,GQ,connectivity,value_PHI,value_diff_PHI,coord,GQ_weights,Z,L_cutoff,L_potential,Vh_in,rho_in_n,N,NumElmNods,ppP,NLCC_flag,exc_flag,pp_flag);
    else
        if exc_flag == 2
            grad_rho_in_n = zeros(E*GQ,LLQ,3);
            grad_rho_in{1} = grad_rho_in_n;
        else
            grad_rho_in = [];
            grad_rho_in_n = [];
        end
        Hnl_blk = zeros(N*L_cutoff^2,N*L_cutoff^2);
        ppP = [];
    end

    %% self consistent field loop start
    while  sc_count < 40 && (error_energy > acceptable_error_energy || error_density > acceptable_error_density)
        error_density_mat = [];   

        %%   schrodinger equation solver  
         [M_blk,H_blk]= schrodingerEqnWithVhandVxcCombined(E,GQ,connectivity,value_PHI,value_diff_PHI,coord,GQ_weights,Z,Zion,L_cutoff,L_potential,Vh_in,rho_in_n,rho_c,N,NumElmNods,LLQ,LLQ_weights,LLQ_points,ppP,grad_rho_in_n,grad_rho_c,NLCC_flag,exc_flag,pp_flag,external_field_flag,L_poisson,DCoefficients,potential_multiplier);    
         [Clm,Eps]=eigs(H_blk+Hnl_blk,M_blk,N*L_cutoff^2,-999);      
         [fractional_occupancy, eps_f, ent] = fractional_occupancy_function(Z,Zion,Eps,pp_flag);

        %% density update
        [rho_out_n,grad_rho_out_n] = densityUpdateMPandDiff(E,N,GQ,value_PHI,value_diff_PHI,coord,Clm,connectivity,fractional_occupancy,NumElmNods,LLQ,LLQ_points,L_cutoff,exc_flag);
        rho_out{sc_count} = rho_out_n;
        if exc_flag == 2 
            grad_rho_out{sc_count} = grad_rho_out_n;
            grad_rho_in_old = grad_rho_in_n;
        else
            grad_rho_in_old = [];
            grad_rho_out = [];
        end

       %% mixing   
        rho_in_old = rho_in_n;
        sc_count = sc_count + 1;
        [rho_in_n,grad_rho_in_n]= mixingCnV(sc_count,rho_in,rho_out,grad_rho_in,grad_rho_out,E,GQ,  GQ_weights,value_PHI, value_diff_PHI,connectivity,coord,NumElmNods,LLQ,LLQ_points,LLQ_weights,NLCC_flag,exc_flag,pp_flag);
        rho_in{sc_count} = rho_in_n; 
        if exc_flag == 2 
            grad_rho_in{sc_count} = grad_rho_in_n;
        end

        %% poisson equation solver with updated density    
        [Vh,Vh_out] = poissonEqnMPFunction(connectivity,coord,E,GQ,GQ_weights,L,L_poisson,N,NumElmNods,p,value_diff_PHI,value_PHI,Z,Zion,LLQ,LLQ_points,LLQ_weights,3,rho_in_n,rho_out_n,NLCC_flag,exc_flag,pp_flag);

        %% Energy   
        [Etotal,Exc,Ecoul,Ekin,~,Emer]= EnergyMPMatrixPPNLCC(Eps,Vh_out, Vh_in,rho_out_n,rho_c, E, GQ,rho_in_old, rho_c, GQ_weights,value_PHI, value_diff_PHI,connectivity,coord,NumElmNods,Z,Zion,fractional_occupancy,LLQ,LLQ_points,LLQ_weights,L_cutoff,N,grad_rho_out_n,grad_rho_in_old,grad_rho_c,NLCC_flag,exc_flag,pp_flag,ent,Clm,eps_f,L_poisson);
        Vh_in = Vh;

        %% Error calculation
        if sc_count < 3
            error_density_mat = ones(size(rho_out_n));
        else
            for i = 1:size(rho_out_n,1)
                for j = 1:size(rho_out_n,2)
                    if mixing_flag == true
                        error_density_mat(i,j) = rho_in_old(i,j) - rho_in_n(i,j);  %
                    else
                        error_density_mat(i,j) = rho_in_old(i,j) - rho_out_n(i,j);  %
                    end
                end
            end 
        end   
        
        fprintf("\n%-10s \n\n      %d\n\n\n ", "Iteration = ",sc_count-1)
        fprintf("%-10s \n\n%21.12f\n\n\n ","Total Energy = ",Etotal)
        error_energy = abs(Emer - Emer_previous) ;
        error_density = (integrator_two_matrices(error_density_mat,error_density_mat,E,GQ,  GQ_weights,value_PHI, value_diff_PHI,connectivity,coord,NumElmNods,LLQ,LLQ_weights,LLQ_points))^(1/2);
        fprintf("%-10s \n\n%21.12f\n\n\n ", "Error Energy = ",error_energy)
        fprintf("%-10s \n\n%21.12f\n\n\n ","Error Density = ",error_density)
        
        
        Emer_previous = Emer;
        if sc_count == 40
            display("Problem unconverged.")
        end
    end
    if density_plot_flag == 1
        [density{m},r_on_ip{m},thetas{m},phis{m}] = densityWithoutDFT(E,N,GQ,value_PHI,value_diff_PHI,coord,Clm,connectivity,NumElmNods,LLQ,L_cutoff,fractional_occupancy);
        save('data_for_plot.mat');
    end
    display("Problem converged.")
    fprintf("%-10s \n\n%21.12f\n\n\n ","Total Energy = ",Etotal)
    fprintf("%-10s \n\n%21.12f\n\n\n ","Kinetic Energy = ",Ekin)
    fprintf("%-10s \n\n%21.12f\n\n\n ","Coluomb Energy = ",Ecoul)
    fprintf("%-10s \n\n%21.12f\n\n\n ","Exchange Correlation Energy = ",Exc)
    Clm_m{m} = Clm;
    Res2(m,1)=Z;
    Res2(m,2)= exc_flag*100 + pp_flag*10 + NLCC_flag;
    Res2(m,3)=L;
    Res2(m,4)=L_cutoff;
    Res2(m,5)=L_poisson;
    Res2(m,6)=L_potential;
    Res2(m,7)=E;
    Res2(m,8)=p;
    Res2(m,9)=LLQ;
    Res2(m,10)=acceptable_error_density;
    Res2(m,11)=F_l(m);
    Res2(m,12)=Etotal;
    Res2(m,13)=Eps(1,1);
    Res2(m,14) = sc_count;
    One_s = Eps(1,1);
    fprintf("%-10s \n\n%21.12f\n\n\n ","1s = ",One_s)
    if Z > 2 
    Two_s = Eps(2,2);
    fprintf("%-10s \n\n%21.12f\n\n\n ","2s = ",Two_s)
    end
    soln = [Etotal Ekin Ecoul Exc One_s];
    T1 = toc;
    fprintf("Time =\n\n%21.12f\n\n\n ", T1)
    Res2(m,15)=T1;
end
end
