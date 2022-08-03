% Copyright (C) 2022 M. A. Yalçın, İ. Temizer
% This file is distributed under the terms of the GNU General Public License.
% See the file COPYING for license details.

function [GQ_points GQ_weights] = GaussianQuadrature(n)

%% Gaussian Quadrature Points and Weights for n Guassian Points
if n == 2
    GQ_points  = zeros(1,2);      GQ_weights = zeros(1,2);
    GQ_points(1) = -0.5773502691896257;     GQ_weights(1) = 1;
    GQ_points(2) = 0.5773502691896257;      GQ_weights(2) = 1;
end

if n == 3
    GQ_points  = zeros(1,3);      GQ_weights = zeros(1,3);
    GQ_points(1) = -sqrt(3/5);      GQ_points(2) = 0;       GQ_points(3) = sqrt(3/5);
    GQ_weights(1) = 5/9;            GQ_weights(2) = 8/9;    GQ_weights(3) = 5/9;
end
if n == 4
    GQ_points  = zeros(1,4);      GQ_weights = zeros(1,4);
    GQ_points(2) = -0.3399810435848563;                       GQ_weights(2) = 0.6521451548625461;
    GQ_points(3) = 0.3399810435848563;     GQ_weights(3) = 0.6521451548625461;
    GQ_points(1) = -0.8611363115940526;      GQ_weights(1) = 0.3478548451374538;
    GQ_points(4) = 0.8611363115940526;     GQ_weights(4) = 0.3478548451374538;
    
end
if n == 5
    GQ_points  = zeros(1,5);      GQ_weights = zeros(1,5);
    GQ_points(3) = 0;                       GQ_weights(3) = 0.5688888888888889;
    GQ_points(2) = -0.5384693101056831;     GQ_weights(2) = 0.4786286704993665;
    GQ_points(4) = 0.5384693101056831;      GQ_weights(4) = 0.4786286704993665;
    GQ_points(1) = -0.9061798459386640;     GQ_weights(1) = 0.2369268850561891;
    GQ_points(5) = 0.9061798459386640;      GQ_weights(5) = 0.2369268850561891;
end
if n == 6
    GQ_points  = zeros(1,6);      GQ_weights = zeros(1,6);
    GQ_weights(5) = 0.3607615730481386;     GQ_points(5) = 0.6612093864662645;
    GQ_weights(2) = 0.3607615730481386;     GQ_points(2) = -0.6612093864662645;
    GQ_weights(3) = 0.4679139345726910;     GQ_points(3) = -0.2386191860831969;
    GQ_weights(4) = 0.4679139345726910;     GQ_points(4) = 0.2386191860831969;
    GQ_weights(1) = 0.1713244923791704;     GQ_points(1) = -0.9324695142031521;  
    GQ_weights(6) = 0.1713244923791704;     GQ_points(6) = 0.9324695142031521;
end
if n == 7
    GQ_points  = zeros(1,7);      GQ_weights = zeros(1,7);
    GQ_weights(4) = 0.4179591836734694;     GQ_points(4) = 0;
    GQ_weights(5) = 0.3818300505051189;     GQ_points(5) = 0.4058451513773972;
    GQ_weights(3) = 0.3818300505051189;     GQ_points(3) = -0.4058451513773972;
    GQ_weights(2) = 0.2797053914892766;     GQ_points(2) = -0.7415311855993945;
    GQ_weights(6) = 0.2797053914892766;     GQ_points(6) = 0.7415311855993945;  
    GQ_weights(1) = 0.1294849661688697;     GQ_points(1) = -0.9491079123427585;
    GQ_weights(7) = 0.1294849661688697;     GQ_points(7) = 0.9491079123427585;
end
if n == 8
    GQ_points  = zeros(1,8);      GQ_weights = zeros(1,8);
    GQ_weights(4) = 0.3626837833783620;     GQ_points(4) = -0.1834346424956498;
    GQ_weights(5) = 0.3626837833783620;     GQ_points(5) = 0.1834346424956498;
    GQ_weights(3) = 0.3137066458778873;     GQ_points(3) = -0.5255324099163290;
    GQ_weights(6) = 0.3137066458778873;     GQ_points(6) = 0.5255324099163290;
    GQ_weights(2) = 0.2223810344533745;     GQ_points(2) = -0.7966664774136267;  
    GQ_weights(7) = 0.2223810344533745;     GQ_points(7) = 0.7966664774136267;
    GQ_weights(1) = 0.1012285362903763;     GQ_points(1) = -0.9602898564975363;
    GQ_weights(8) = 0.1012285362903763;     GQ_points(8) = 0.9602898564975363;
end
if n == 9
    GQ_points  = zeros(1,9);      GQ_weights = zeros(1,9);
    GQ_weights(5) = 0.3302393550012598;     GQ_points(5) = 0;
    GQ_weights(2) = 0.1806481606948574;     GQ_points(2) = -0.8360311073266358;
    GQ_weights(8) = 0.1806481606948574;     GQ_points(8) = 0.8360311073266358;
    GQ_weights(1) = 0.0812743883615744;     GQ_points(1) = -0.9681602395076261;
    GQ_weights(9) = 0.0812743883615744;     GQ_points(9) = 0.9681602395076261;  
    GQ_weights(4) = 0.3123470770400029;     GQ_points(4) = -0.3242534234038089;
    GQ_weights(6) = 0.3123470770400029;     GQ_points(6) = 0.3242534234038089;
    GQ_weights(3) = 0.2606106964029354;     GQ_points(3) = -0.6133714327005904;
    GQ_weights(7) = 0.2606106964029354;     GQ_points(7) = 0.6133714327005904;     
end

