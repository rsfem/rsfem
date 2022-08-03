% Copyright (C) 2022 M. A. Yalçın, İ. Temizer
% This file is distributed under the terms of the GNU General Public License.
% See the file COPYING for license details.

function  [coord connectivity]= BuildUp(E,n,p,NumElmNods,L)
% Domain Boundaries
x_min = 0;  x_max = L;  

% element size calculation
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

% coordinates
coord=[];
coord(1)=x_min;
no = 1;
for j_c = 1:E
    for k_c = 1:p
      no = no + 1;
      coord(no) = coord(no-1) + he_lg(j_c)/p;
    end
end

% Connectivity Array
connectivity = zeros(E,NumElmNods);
counter = 1;
for i = 1:E
    for j = 1:NumElmNods
        connectivity(i,j) = counter;
        counter = counter + 1;
    end
    counter = counter - 1;
end
