% Copyright (C) 2022 M. A. Yalçın, İ. Temizer
% This file is distributed under the terms of the GNU General Public License.
% See the file COPYING for license details.

function [PHI diffPHI] = ShapeFunctionandDiff(p,GQ_points)
x_l = linspace(-1,1,p+1);
PHI=ones(p+1,length(GQ_points));
PHI2=ones(p+1,length(GQ_points));
diffPHI=zeros(p+1,length(GQ_points));
for k=1:length(GQ_points)
    z = GQ_points(k);
    for i=1:p+1
        for j=1:p+1
            if i ~= j               
                PHI(i,k) = PHI(i,k) * (z - x_l(j))/(x_l(i) - x_l(j));
                if x_l(j) == 0
                    PHI2(i,k) = PHI2(i,k) * (1)/(x_l(i) - x_l(j));
                else
                    PHI2(i,k) = PHI2(i,k) * (- x_l(j))/(x_l(i) - x_l(j));
                end
            end
        end   
    end
end

for k=1:length(GQ_points)
    z = GQ_points(k);
    for i=1:p+1
        for j=1:p+1
            if i ~= j   
                if isnan(PHI(i,k) / (z - x_l(j))) == 1
                    diffPHI(i,k) = diffPHI(i,k) + PHI2(i,k);
                else
                    diffPHI(i,k) = diffPHI(i,k) + PHI(i,k) / (z - x_l(j)); 
                end
            end
        end   
    end
end
end
