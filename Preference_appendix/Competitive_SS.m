%%%%% This file generates simulation results for competitive model for 
%%%%% preference setting, as well as competitive emissions 
%%%%% data for Figure G.2 in the appendix.

% Authors: Mike Barnett, Jieyao Wang
% Last update: Nov 25, 2019

close all
clear all
clc

%% competitive simulation result

% steady states
kappa = 0.032;
i = 0.09;
q = 2.5;
e = 0.015;
rho = 0.02;
y = log(0.98);
delta = 0.01;

nu = delta.*kappa./e;
c = delta.*(1-kappa).*q./(1-nu);
phi_prod = 1;
phi_1 = (phi_prod.*q-1)./i;
phi_0 = phi_prod./phi_1;
mu_k = rho-phi_0.*log(1+phi_1.*i);
log_j = (log(nu)-log(delta.*(1-kappa)./c)+log(e-delta));
j = exp(log_j);
psi_1 = (e-delta)./(rho+e);
psi_0 = exp(log(rho+e)+psi_1.*y-psi_1.*log_j);
alpha = i+j+c;

% simulation
filename2 = [pwd,'/SS_competitive_Sims'];

T = 100; % 100 years
pers = 4*T; % quarterly
dt = T/pers;
nDims = 3;
its = 1;

R_max = exp(9);
K_max = exp(18);
R_min = exp(0);
K_min = exp(0);
F_max = 4000;
F_min = 0;

R_0 = 650;
K_0 = 80/alpha;
F_0 = (870-580);

muR = @(x) -e+psi_0.*(j.*x(:,2)./x(:,1)).^psi_1;
muK = @(x) (mu_k + phi_0.*log(1+i.*phi_1));
muF = @(x) e.*x(:,1);

sigmaR = @(x) [zeros(size(x(:,1:3)))];
sigmaK = @(x) [zeros(size(x(:,1:3)))];
sigmaF = @(x) [zeros(size(x(:,1:3)))];


upperBounds = [R_max,K_max,F_max];
lowerBounds = [R_min,K_min,F_min];

hists = zeros(pers,nDims,its);
hists2 = hists;
e_hists = zeros(pers,its);
e_hists2 = e_hists;
j_hists = zeros(pers,its);
j_hists2 = j_hists;
i_k_hists = zeros(pers,its);
i_k_hists2 = i_k_hists;

for iters = 1:its
    
hist2 = zeros(pers,nDims);
e_hist2 = zeros(pers,1);
i_k_hist2 = zeros(pers,1);
j_hist2 = zeros(pers,1);


hist2(1,:) = [R_0,K_0,F_0];
e_hist2(1) =  e.*hist2(1,1);
i_k_hist2(1) =  i.*hist2(1,2);
j_hist2(1) =  j.*hist2(1,2);

for jj = 2:pers
shock = normrnd(0,sqrt(dt), 1, nDims);
hist2(jj,1) = max(min(hist2(jj-1,1).*exp((muR(hist2(jj-1,:))-0.5.*sum((sigmaR(hist2(jj-1,:))).^2) )* dt ...
                                  +sigmaR(hist2(jj-1,:))* shock'), upperBounds(:,1)), lowerBounds(:,1));
hist2(jj,2) = max(min(hist2(jj-1,2).*exp((muK(hist2(jj-1,:))-0.5.*sum((sigmaK(hist2(jj-1,:))).^2) )* dt ...
                                  +sigmaK(hist2(jj-1,:))* shock'), upperBounds(:,2)), lowerBounds(:,2));                              
hist2(jj,3) = max(min(hist2(jj-1,3) + muF(hist2(jj-1,:)) * dt + sigmaF(hist2(jj-1,:))* shock', upperBounds(:,3)), lowerBounds(:,3)); 

e_hist2(jj) = e.*hist2(jj-1,1);
i_k_hist2(jj) = i.*hist2(jj-1,2);
j_hist2(jj) =  j.*hist2(jj-1,2);

end

hists2(:,:,iters) = hist2;
e_hists2(:,iters) = e_hist2;
i_k_hists2(:,iters) = i_k_hist2;
j_hists2(:,iters) = j_hist2; 

end

% save results
save(filename2);

%% Generate Competitive emissions data for Figure G.2 in online appendix
e_values = mean(e_hists2,2);
s_e = strcat('e_comp.mat');
save(s_e,'e_values');
