% Copyright (C) 2022 M. A. Yalçın, İ. Temizer
% This file is distributed under the terms of the GNU General Public License.
% See the file COPYING for license details.

function ppP = ppParameters(Z,NLCC_flag,exc_flag,pp_flag)
%% 
% r0 C1   C2   C3   C4   0    0
% rs Hs11 Hs22 Hs33 Hs12 Hs13 Hs23
% rs Ks11 Ks22 Ks33 Ks12 Ks13 Ks23 no spin dependence, probably will always be null
% rp Hp11 Hp22 Hp33 Hp12 Hp13 Hp23
% rp Kp11 Kp22 Kp33 Kp12 Kp13 Kp23
% rd Hd11 Hd22 Hd33 Hd12 Hd13 Hd23
% rd Kd11 Kd22 Kd33 Kd12 Kd13 Kd23

ppP = zeros(7,7);
if Z == 1  % hydrogen
    if NLCC_flag == 0
        if exc_flag == 1
            ppP(1,1) = 0.2;
            ppP(1,2) = -4.180237    ;
            ppP(1,3) = 0.725075     ;
        elseif exc_flag == 2
            ppP(1,1) = 0.2;
            ppP(1,2) = -4.17890044    ;
            ppP(1,3) = 0.72446331     ;
        end   
    else %NLCC = 1 strictly for GGA
        ppP(1,1) = 0.228727444559;
        ppP(1,2) = -3.722900963668    ;
        ppP(1,3) = 0.657596429179    ;
    end       
elseif Z == 3 % lithium
    %% Local part
    if NLCC_flag == 0
        if exc_flag == 1 % LDA
            ppP(1,1) = 0.787553     ; 
            ppP(1,2) = -1.892612    ;
            ppP(1,3) = 0.286060     ;
        else            % GGA cant work without NLCC
            error("GGA cant work without NLCC for Li")
        end
    else % NLCC on, strictly for GGA        
		  ppP(1,1) =   0.811627593845;
	      ppP(1,2) = - 1.07248735835 ;
	      ppP(1,3) =   0.201649152445;
    end
    
    %% Nonlocal part
    if NLCC_flag == 0
        if exc_flag == 1 % LDA
            ppP(2,1) = 0.666375     ; 
            ppP(2,2) = 1.858811     ;
		    ppP(4,1) = 1.079306     ; 
		    ppP(4,2) = -0.005895    ; 
        else % GGA
            % only local 
        end
    else % NLCC on, strictly for GGA
        ppP(2,1) = 0.6975989912;
        ppP(2,2) = 1.55506211655;  
    end
    
elseif Z == 13 % aluminum
    %% Local Part
    if NLCC_flag == 0 
        if exc_flag == 1 % LDA
    	ppP(1,1) =  0.45   ; 
        ppP(1,2) =  -8.491351   ;
        ppP(1,3) =  0    ;
        else % GGA
        ppP(1,1) =   0.45000000;
		ppP(1,2) = - 7.55476126;  
        ppP(1,3) =  0    ;
        end
    else % NLCC on, strictly for GGA
		ppP(1,1) =   0.60224003597;
		ppP(1,2) = - 5.18802089974;
		ppP(1,3) =   0.6961574711;
    
    end
    
    %% Nonlocal Part
    if NLCC_flag == 0
        if exc_flag == 1 % LDA
        ppP(2,1) =  0.460104    ;  
        ppP(2,2) =  5.088340    ; 
        ppP(2,3) =  2.6797    ;
     	ppP(3,1) = 0            ; 
        ppP(3,2) = 0            ;
        ppP(3,3) = 0            ;
     	ppP(4,1) = 0.536744            ; 
        ppP(4,2) = 2.193438            ;
        ppP(4,3) = 0            ;
        else % GGA
          ppP(2,1) =   0.48743529;
		  ppP(2,2) =   6.95993832;
		  ppP(2,3) =   2.43847659;
		  ppP(2,5) = - 1.88883584;
	  
		  ppP(4,1) =   0.56218949;
		  ppP(4,2) =   1.86529857;
		  ppP(4,3) = 0;
        end
    else % NLCC on, strictly for GGA
          ppP(2,1) =   0.451219678;
		  ppP(2,2) =   6.175029169180;
		  ppP(2,3) =   3.6822889039;
		  ppP(2,5) = - 1.91493592966;		
	  
		  ppP(4,1) =   0.53746997212;
		  ppP(4,2) =   2.272958057;
		  ppP(4,3) =   0;
    end
end
