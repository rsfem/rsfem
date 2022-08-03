% Copyright (C) 2022 M. A. Yalçın, İ. Temizer
% This file is distributed under the terms of the GNU General Public License.
% See the file COPYING for license details.

clc
clear
close all

%% inputs
% for help: call inputs_help from command window 
E_l = [8*ones(1,1)];
p_l = [6*ones(1,1)];
F_l = [3 ];
L_cutoff_l = [3*ones(1,1)];
LLQ_l = [110*ones(1,1) ];
mp_l = [1 ;0 ;0 ];
conf = [2 1 1];
density_plot_flag=1;

%% output
[Res Clm density thetas phis] = sch_density_energy_function(E_l,L_cutoff_l,F_l,p_l,mp_l,LLQ_l,conf,density_plot_flag);
