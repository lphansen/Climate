function [x,w] = quad_points_legendre(n)
%
% Generate nodes and weights for
% Legendre-Gauss quadrature on [-1,1].
% Note that x is a column vector
% and w is a row vector.
% Taken from John Gubner (2014)
%
u = sqrt(1./(4-1./[1:n-1].^2)); % upper diag.
[V,Lambda] = eig(diag(u,1)+diag(u,-1));
[x,i] = sort(diag(Lambda));
Vtop = V(1,:);
Vtop = Vtop(i);
w = 2*Vtop.^2;
