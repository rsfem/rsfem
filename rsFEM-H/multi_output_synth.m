% Copyright (C) 2022 M. A. Yalçın, İ. Temizer
% This file is distributed under the terms of the GNU General Public License.
% See the file COPYING for license details.

clc
clear
close all

%% inputs
% for help: call inputs_help from command window 
% specify test parameters in a string
potential_multiplier_l = [0.3];
E_l = [8];
p_l = [6];
L_cutoff_l = [4];
density_plot_flag = 1;
alpha = -1;

%% outputs
[Res Eps density r_on_ip thetas phis] = HYDROGEN_STARK_FUNCTION(p_l, E_l,L_cutoff_l,potential_multiplier_l,density_plot_flag,alpha);
