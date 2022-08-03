% Copyright (C) 2022 M. A. Yalçın, İ. Temizer
% This file is distributed under the terms of the GNU General Public License.
% See the file COPYING for license details.

clear
clc
close all

load('data_for_plot.mat')

yy = transpose(r_on_ip{1}) * (sin(thetas{1}));
zz = transpose(r_on_ip{1}) * cos(thetas{1});
yl = [-5.5 5.5];
yl = [-5 6];
lims = [0 0.2];
colormap(jet);

contour(yy,zz,density{1,1},[0.002:0.02:0.2 ]);

xlabel('y')
ylabel('z','rotation',0,'VerticalAlignment','middle')
xlim(yl)
ylim(yl)
cc1= colorbar('location','northoutside','Ticks',[0 0.1  0.2],'AxisLocation','out');
c1.CLim = lims;
daspect([1,1,1]);
yline(0,'--');
set(cc1,'FontSize',7);
Y=get(cc1,'Position');
Y(4)=Y(4)*0.6; 
set(cc1,'Position',Y)