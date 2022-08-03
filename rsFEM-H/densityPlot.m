% Copyright (C) 2022 M. A. Yalçın, İ. Temizer
% This file is distributed under the terms of the GNU General Public License.
% See the file COPYING for license details.

clear
clc
close all

load('data_for_plot.mat')

yy = transpose(r_on_ip{1}) * (sin(thetas{1}));
zz = transpose(r_on_ip{1}) * cos(thetas{1});

yl = [-4 4];
zl = [-3.5 4.5];
lims = [0 0.2];
colormap(jet)
tiledlayout(1,4)


contour(yy,zz,density{1,1},[0.002:0.04:0.2 ])
xlabel('y')
ylabel('z','rotation',0,'VerticalAlignment','middle')
xlim(yl)
ylim(zl)

daspect([	1,1,1])
yline(0,'--')

