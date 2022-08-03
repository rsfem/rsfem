% Copyright (C) 2022 M. A. Yalçın, İ. Temizer
% This file is distributed under the terms of the GNU General Public License.
% See the file COPYING for license details.

function flagCheck(exc_flag,pp_flag,NLCC_flag)
if exc_flag == 1
    if pp_flag == 0
        string1 = "All-electron & LDA selected. Slater exchange and Vosko, Wilk & Nusair correlation functionals will be used";
    	if NLCC_flag == 0
        	string2 = " \n";   
        elseif NLCC_flag == 1
        	error("Nonlinear core correction cannot be used with All-electron setting.\n");   
        end
    elseif pp_flag == 1
        string1 = "Pseudopotential & LDA selected. Teter 93 exchange and correlation functional will be used. Make sure to enter all PP parameters to following functions:";
        if NLCC_flag == 0
        	string2 = "Nonlinear core correction will not be considered.\n";   
        elseif NLCC_flag == 1
        	string2 = "Nonlinear core correction will be considered. Make sure to enter all NLCC parameters to following functions:\nnlParameter.m\nnlccParameter.m";   
        end
    end
elseif exc_flag == 2
    if pp_flag == 0
        string1 = "All-electron & GGA selected. Perdew, Burke & Ernzerhof exchange and correlation functional will be used";
        if NLCC_flag == 0
            string2 = " \n";   
        elseif NLCC_flag == 1
            error("Nonlinear core correction cannot be used with All-electron setting.");   
        end
    elseif pp_flag == 1
        string1 = "Pseudopotential & GGA selected. Perdew, Burke & Ernzerhof exchange and correlation functional will be used. Make sure to enter all PP parameters to following functions:\nppParameters.m";
        if NLCC_flag == 0
            string2 = "Nonlinear core correction will not be considered.\n";  
        elseif NLCC_flag == 1
        	string2 = "Nonlinear core correction will be considered. Make sure to enter all NLCC parameters to following functions:\nnlParameter.m\nnlccParameter.m";   
        end
    end
end
strB = string1 + "\n" + string2 + "\n";
fprintf(strB)

    
