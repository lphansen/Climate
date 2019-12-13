%%%%% This file generates simulation result for the Ambiguity Neutral model.
% Authors: Mike Barnett, Jieyao Wang
% Last update: Nov 20,2019

close all
clear all
clc

load([pwd, '/HJB_Growth_Neutral_1e14']);
filename = 'HJB_Growth_Neutral_1e14';

filename2 = [filename,'_Sims_JY'];

T = 100; 
pers = 4*T; 
dt = T/pers;
nDims = 4;
its = 1;

R_0 = 650;
K_0 = 80./A_O;
T_0 = (870-580);

efunc = griddedInterpolant(r_mat,t_mat,k_mat,e,'linear');
ffunc = griddedInterpolant(r_mat,t_mat,k_mat,f,'linear');
i_kfunc = griddedInterpolant(r_mat,t_mat,k_mat,i_k,'linear');

v_drfunc = griddedInterpolant(r_mat,t_mat,k_mat,v0_dr,'linear');
v_dtfunc = griddedInterpolant(r_mat,t_mat,k_mat,v0_dt,'linear');
v_dkfunc = griddedInterpolant(r_mat,t_mat,k_mat,v0_dk,'linear');
v_func = griddedInterpolant(r_mat,t_mat,k_mat,out_comp,'linear');

pi_tilde_1func = griddedInterpolant(r_mat,t_mat,k_mat,pi_tilde_1_norm,'linear');
pi_tilde_2func = griddedInterpolant(r_mat,t_mat,k_mat,pi_tilde_2_norm,'linear');
pi_tilde_3func = griddedInterpolant(r_mat,t_mat,k_mat,pi_tilde_3_norm,'linear');
pi_tilde_4func = griddedInterpolant(r_mat,t_mat,k_mat,pi_tilde_4_norm,'linear');
pi_tilde_5func = griddedInterpolant(r_mat,t_mat,k_mat,pi_tilde_5_norm,'linear');
pi_tilde_6func = griddedInterpolant(r_mat,t_mat,k_mat,pi_tilde_6_norm,'linear');
pi_tilde_7func = griddedInterpolant(r_mat,t_mat,k_mat,pi_tilde_7_norm,'linear');
pi_tilde_8func = griddedInterpolant(r_mat,t_mat,k_mat,pi_tilde_8_norm,'linear');
pi_tilde_9func = griddedInterpolant(r_mat,t_mat,k_mat,pi_tilde_9_norm,'linear');

e_func = @(x) efunc(log(x(:,1)),x(:,3),log(x(:,2)));
f_func = @(x) max(ffunc(log(x(:,1)),x(:,3),log(x(:,2))),0);
i_k_func = @(x) i_kfunc(log(x(:,1)),x(:,3),log(x(:,2)));

v_dr_func = @(x) v_drfunc(log(x(:,1)),x(:,3),log(x(:,2)));
v_dt_func = @(x) v_dtfunc(log(x(:,1)),x(:,3),log(x(:,2)));
v_dk_func = @(x) v_dkfunc(log(x(:,1)),x(:,3),log(x(:,2)));
v_func = @(x) v_func(log(x(:,1)),x(:,3),log(x(:,2)));

pi_tilde_1_func = @(x) pi_tilde_1func(log(x(:,1)),x(:,3),log(x(:,2)));
pi_tilde_2_func = @(x) pi_tilde_2func(log(x(:,1)),x(:,3),log(x(:,2)));
pi_tilde_3_func = @(x) pi_tilde_3func(log(x(:,1)),x(:,3),log(x(:,2)));
pi_tilde_4_func = @(x) pi_tilde_4func(log(x(:,1)),x(:,3),log(x(:,2)));
pi_tilde_5_func = @(x) pi_tilde_5func(log(x(:,1)),x(:,3),log(x(:,2)));
pi_tilde_6_func = @(x) pi_tilde_6func(log(x(:,1)),x(:,3),log(x(:,2)));
pi_tilde_7_func = @(x) pi_tilde_7func(log(x(:,1)),x(:,3),log(x(:,2)));
pi_tilde_8_func = @(x) pi_tilde_8func(log(x(:,1)),x(:,3),log(x(:,2)));
pi_tilde_9_func = @(x) pi_tilde_9func(log(x(:,1)),x(:,3),log(x(:,2)));

REfunc = griddedInterpolant(r_mat,t_mat,k_mat,RE,'linear');
RE_func = @(x) REfunc(log(x(:,1)),x(:,3),log(x(:,2)));

a1 = (gamma0(1)+gamma1(1).*t_bar+0.5.*gamma2(1).*(t_bar.^2));
b1 = t_mat.*(gamma1(1)+gamma2(1).*t_bar);
c1 = gamma2(1).*(t_mat.^2);

a2 = (gamma0(2)+gamma1(2).*t_bar+0.5.*gamma2(2).*(t_bar.^2));
b2 = t_mat.*(gamma1(2)+gamma2(2).*t_bar);
c2 = gamma2(2).*(t_mat.^2);

a3 = (gamma0(3)+gamma1(3).*t_bar+0.5.*gamma2(3).*(t_bar.^2));
b3 = t_mat.*(gamma1(3)+gamma2(3).*t_bar);
c3 = gamma2(3).*(t_mat.^2);

a4 = (gamma0(4)+gamma1(4).*t_bar+0.5.*gamma2(4).*(t_bar.^2));
b4 = t_mat.*(gamma1(4)+gamma2(4).*t_bar);
c4 = gamma2(4).*(t_mat.^2);

a5 = (gamma0(5)+gamma1(5).*t_bar+0.5.*gamma2(5).*(t_bar.^2));
b5 = t_mat.*(gamma1(5)+gamma2(5).*t_bar);
c5 = gamma2(5).*(t_mat.^2);

a6 = (gamma0(6)+gamma1(6).*t_bar+0.5.*gamma2(6).*(t_bar.^2));
b6 = t_mat.*(gamma1(6)+gamma2(6).*t_bar);
c6 = gamma2(6).*(t_mat.^2);

a7 = (gamma0(7)+gamma1(7).*t_bar+0.5.*gamma2(7).*(t_bar.^2));
b7 = t_mat.*(gamma1(7)+gamma2(7).*t_bar);
c7 = gamma2(7).*(t_mat.^2);

a8 = (gamma0(8)+gamma1(8).*t_bar+0.5.*gamma2(8).*(t_bar.^2));
b8 = t_mat.*(gamma1(8)+gamma2(8).*t_bar);
c8 = gamma2(8).*(t_mat.^2);

a9 = (gamma0(9)+gamma1(9).*t_bar+0.5.*gamma2(9).*(t_bar.^2));
b9 = t_mat.*(gamma1(9)+gamma2(9).*t_bar);
c9 = gamma2(9).*(t_mat.^2);

dmg_tilt_1 = a1+b1.*beta_tilde_1+0.5.*c1.*(beta_tilde_1.^2)+0.5.*c1./lambda_tilde_1;    
dmg_tilt_2 = a2+b2.*beta_tilde_2+0.5.*c2.*(beta_tilde_2.^2)+0.5.*c2./lambda_tilde_2;    
dmg_tilt_3 = a3+b3.*beta_tilde_3+0.5.*c3.*(beta_tilde_3.^2)+0.5.*c3./lambda_tilde_3;    
dmg_tilt_4 = a4+b4.*beta_tilde_4+0.5.*c4.*(beta_tilde_4.^2)+0.5.*c4./lambda_tilde_4;    
dmg_tilt_5 = a5+b5.*beta_tilde_5+0.5.*c5.*(beta_tilde_5.^2)+0.5.*c5./lambda_tilde_5;    
dmg_tilt_6 = a6+b6.*beta_tilde_6+0.5.*c6.*(beta_tilde_6.^2)+0.5.*c6./lambda_tilde_6;    
dmg_tilt_7 = a7+b7.*beta_tilde_7+0.5.*c7.*(beta_tilde_7.^2)+0.5.*c7./lambda_tilde_7;    
dmg_tilt_8 = a8+b8.*beta_tilde_8+0.5.*c8.*(beta_tilde_8.^2)+0.5.*c8./lambda_tilde_8;    
dmg_tilt_9 = a9+b9.*beta_tilde_9+0.5.*c9.*(beta_tilde_9.^2)+0.5.*c9./lambda_tilde_9;    

dmg_1 = a1+b1.*beta_f+0.5.*c1.*(beta_f.^2)+0.5.*c1./lambda;    
dmg_2 = a2+b2.*beta_f+0.5.*c2.*(beta_f.^2)+0.5.*c2./lambda;    
dmg_3 = a3+b3.*beta_f+0.5.*c3.*(beta_f.^2)+0.5.*c3./lambda;    
dmg_4 = a4+b4.*beta_f+0.5.*c4.*(beta_f.^2)+0.5.*c4./lambda;    
dmg_5 = a5+b5.*beta_f+0.5.*c5.*(beta_f.^2)+0.5.*c5./lambda;    
dmg_6 = a6+b6.*beta_f+0.5.*c6.*(beta_f.^2)+0.5.*c6./lambda;    
dmg_7 = a7+b7.*beta_f+0.5.*c7.*(beta_f.^2)+0.5.*c7./lambda;    
dmg_8 = a8+b8.*beta_f+0.5.*c8.*(beta_f.^2)+0.5.*c8./lambda;    
dmg_9 = a9+b9.*beta_f+0.5.*c9.*(beta_f.^2)+0.5.*c9./lambda;       


base_driftKfunc1 = griddedInterpolant(r_mat,t_mat,k_mat,dmg_1,'linear');
base_driftKfunc2 = griddedInterpolant(r_mat,t_mat,k_mat,dmg_2,'linear');
base_driftKfunc3 = griddedInterpolant(r_mat,t_mat,k_mat,dmg_3,'linear');
base_driftKfunc4 = griddedInterpolant(r_mat,t_mat,k_mat,dmg_4,'linear');
base_driftKfunc5 = griddedInterpolant(r_mat,t_mat,k_mat,dmg_5,'linear');
base_driftKfunc6 = griddedInterpolant(r_mat,t_mat,k_mat,dmg_6,'linear');
base_driftKfunc7 = griddedInterpolant(r_mat,t_mat,k_mat,dmg_7,'linear');
base_driftKfunc8 = griddedInterpolant(r_mat,t_mat,k_mat,dmg_8,'linear');
base_driftKfunc9 = griddedInterpolant(r_mat,t_mat,k_mat,dmg_9,'linear');
base_driftK_func1 = @(x) base_driftKfunc1(log(x(:,1)),x(:,3),log(x(:,2)));
base_driftK_func2 = @(x) base_driftKfunc2(log(x(:,1)),x(:,3),log(x(:,2)));
base_driftK_func3 = @(x) base_driftKfunc3(log(x(:,1)),x(:,3),log(x(:,2)));
base_driftK_func4 = @(x) base_driftKfunc4(log(x(:,1)),x(:,3),log(x(:,2)));
base_driftK_func5 = @(x) base_driftKfunc5(log(x(:,1)),x(:,3),log(x(:,2)));
base_driftK_func6 = @(x) base_driftKfunc6(log(x(:,1)),x(:,3),log(x(:,2)));
base_driftK_func7 = @(x) base_driftKfunc7(log(x(:,1)),x(:,3),log(x(:,2)));
base_driftK_func8 = @(x) base_driftKfunc8(log(x(:,1)),x(:,3),log(x(:,2)));
base_driftK_func9 = @(x) base_driftKfunc9(log(x(:,1)),x(:,3),log(x(:,2)));

tilt_driftKfunc1 = griddedInterpolant(r_mat,t_mat,k_mat,dmg_tilt_1,'linear');
tilt_driftKfunc2 = griddedInterpolant(r_mat,t_mat,k_mat,dmg_tilt_2,'linear');
tilt_driftKfunc3 = griddedInterpolant(r_mat,t_mat,k_mat,dmg_tilt_3,'linear');
tilt_driftKfunc4 = griddedInterpolant(r_mat,t_mat,k_mat,dmg_tilt_4,'linear');
tilt_driftKfunc5 = griddedInterpolant(r_mat,t_mat,k_mat,dmg_tilt_5,'linear');
tilt_driftKfunc6 = griddedInterpolant(r_mat,t_mat,k_mat,dmg_tilt_6,'linear');
tilt_driftKfunc7 = griddedInterpolant(r_mat,t_mat,k_mat,dmg_tilt_7,'linear');
tilt_driftKfunc8 = griddedInterpolant(r_mat,t_mat,k_mat,dmg_tilt_8,'linear');
tilt_driftKfunc9 = griddedInterpolant(r_mat,t_mat,k_mat,dmg_tilt_9,'linear');
tilt_driftK_func1 = @(x) tilt_driftKfunc1(log(x(:,1)),x(:,3),log(x(:,2)));
tilt_driftK_func2 = @(x) tilt_driftKfunc2(log(x(:,1)),x(:,3),log(x(:,2)));
tilt_driftK_func3 = @(x) tilt_driftKfunc3(log(x(:,1)),x(:,3),log(x(:,2)));
tilt_driftK_func4 = @(x) tilt_driftKfunc4(log(x(:,1)),x(:,3),log(x(:,2)));
tilt_driftK_func5 = @(x) tilt_driftKfunc5(log(x(:,1)),x(:,3),log(x(:,2)));
tilt_driftK_func6 = @(x) tilt_driftKfunc6(log(x(:,1)),x(:,3),log(x(:,2)));
tilt_driftK_func7 = @(x) tilt_driftKfunc7(log(x(:,1)),x(:,3),log(x(:,2)));
tilt_driftK_func8 = @(x) tilt_driftKfunc8(log(x(:,1)),x(:,3),log(x(:,2)));
tilt_driftK_func9 = @(x) tilt_driftKfunc9(log(x(:,1)),x(:,3),log(x(:,2)));


Gamma_base = @(x) weight1.*base_driftK_func1(x)+ weight2.*base_driftK_func2(x)...
              + weight3.*base_driftK_func3(x)+ weight4.*base_driftK_func4(x)...
              + weight5.*base_driftK_func5(x)+ weight6.*base_driftK_func6(x)...
              + weight7.*base_driftK_func7(x)+ weight8.*base_driftK_func8(x)...
              + weight9.*base_driftK_func9(x);

Gamma_tilted = @(x) pi_tilde_1_func(x).*tilt_driftK_func1(x)+pi_tilde_2_func(x).*tilt_driftK_func2(x)...
                 +pi_tilde_3_func(x).*tilt_driftK_func3(x)+pi_tilde_4_func(x).*tilt_driftK_func4(x)...
                 +pi_tilde_5_func(x).*tilt_driftK_func5(x)+pi_tilde_6_func(x).*tilt_driftK_func6(x)...
                 +pi_tilde_7_func(x).*tilt_driftK_func7(x)+pi_tilde_8_func(x).*tilt_driftK_func8(x)...
                 +pi_tilde_9_func(x).*tilt_driftK_func9(x);

muR = @(x) -e_func(x)+Gamma_r.*(f_func(x).*x(:,2)./x(:,1)).^Theta_r;
muK_base = @(x) (Alpha + Gamma.*log(1+i_k_func(x)./Theta)-Gamma_base(x));
muT = @(x) e_func(x).*x(:,1);
muK_tilted = @(x) (Alpha + Gamma.*log(1+i_k_func(x)./Theta)-Gamma_tilted(x));

sigmaR = @(x) [zeros(size(x(:,1:4)))];
sigmaK = @(x) [zeros(size(x(:,1:4)))];
sigmaT = @(x) [zeros(size(x(:,1:4)))];

R_max = exp(r_max);
K_max = exp(k_max);
T_max = t_max;

R_min = exp(r_min);
K_min = exp(k_min);
T_min = t_min;

upperBounds = [R_max,K_max,T_max,K_max];
lowerBounds = [R_min,K_min,T_min,K_min];

hists = zeros(pers,nDims,its);
hists2 = hists;
e_hists = zeros(pers,its);
e_hists2 = e_hists;
f_hists = zeros(pers,its);
f_hists2 = f_hists;
i_k_hists = zeros(pers,its);
i_k_hists2 = i_k_hists;

for iters = 1:its

hist2 = zeros(pers,nDims);
e_hist2 = zeros(pers,1);
i_k_hist2 = zeros(pers,1);
f_hist2 = zeros(pers,1);

RE_hist2 = zeros(pers,1);
pi_tilde_1_hist2 = zeros(pers,1);
pi_tilde_2_hist2 = zeros(pers,1);
pi_tilde_3_hist2 = zeros(pers,1);
pi_tilde_4_hist2 = zeros(pers,1);
pi_tilde_5_hist2 = zeros(pers,1);
pi_tilde_6_hist2 = zeros(pers,1);
pi_tilde_7_hist2 = zeros(pers,1);
pi_tilde_8_hist2 = zeros(pers,1);
pi_tilde_9_hist2 = zeros(pers,1);

hist2(1,:) = [R_0,K_0,T_0,K_0];
e_hist2(1) =  e_func(hist2(1,:)).*hist2(1,1);
i_k_hist2(1) =  i_k_func(hist2(1,:)).*hist2(1,2);
f_hist2(1) =  f_func(hist2(1,:)).*hist2(1,2);
RE_hist2 = RE_func(hist2(1,:));
pi_tilde_1_hist2 = pi_tilde_1_func(hist2(1,:));
pi_tilde_2_hist2 = pi_tilde_2_func(hist2(1,:));
pi_tilde_3_hist2 = pi_tilde_3_func(hist2(1,:));
pi_tilde_4_hist2 = pi_tilde_4_func(hist2(1,:));
pi_tilde_5_hist2 = pi_tilde_5_func(hist2(1,:));
pi_tilde_6_hist2 = pi_tilde_6_func(hist2(1,:));
pi_tilde_7_hist2 = pi_tilde_7_func(hist2(1,:));
pi_tilde_8_hist2 = pi_tilde_8_func(hist2(1,:));
pi_tilde_9_hist2 = pi_tilde_9_func(hist2(1,:));

for j = 2:pers

shock = normrnd(0,sqrt(dt), 1, nDims);
hist2(j,1) = max(min(hist2(j-1,1).*exp((muR(hist2(j-1,:))-0.5.*sum((sigmaR(hist2(j-1,:))).^2) )* dt ...
                                  +sigmaR(hist2(j-1,:))* shock'), upperBounds(:,1)), lowerBounds(:,1));
hist2(j,2) = max(min(hist2(j-1,2).*exp((muK_base(hist2(j-1,:))-0.5.*sum((sigmaK(hist2(j-1,:))).^2) )* dt ...
                                  +sigmaK(hist2(j-1,:))* shock'), upperBounds(:,2)), lowerBounds(:,2));                              
hist2(j,3) = max(min(hist2(j-1,3) + muT(hist2(j-1,:)) * dt + sigmaT(hist2(j-1,:))* shock', upperBounds(:,3)), lowerBounds(:,3)); 
hist2(j,4) = max(min(hist2(j-1,4).*exp((muK_tilted(hist2(j-1,:))-0.5.*sum((sigmaK(hist2(j-1,:))).^2) )* dt ...
                                  +sigmaK(hist2(j-1,:))* shock'), upperBounds(:,4)), lowerBounds(:,4));  
                              
e_hist2(j) = e_func(hist2(j-1,:)).*hist2(j-1,1);
i_k_hist2(j) = i_k_func(hist2(j-1,:)).*hist2(j-1,2);
f_hist2(j) =  f_func(hist2(j-1,:)).*hist2(j-1,2);

RE_hist2(j) = RE_func(hist2(j-1,:));
pi_tilde_1_hist2(j) = pi_tilde_1_func(hist2(j-1,:));
pi_tilde_2_hist2(j) = pi_tilde_2_func(hist2(j-1,:));
pi_tilde_3_hist2(j) = pi_tilde_3_func(hist2(j-1,:));
pi_tilde_4_hist2(j) = pi_tilde_4_func(hist2(j-1,:));
pi_tilde_5_hist2(j) = pi_tilde_5_func(hist2(j-1,:));
pi_tilde_6_hist2(j) = pi_tilde_6_func(hist2(j-1,:));
pi_tilde_7_hist2(j) = pi_tilde_7_func(hist2(j-1,:));
pi_tilde_8_hist2(j) = pi_tilde_8_func(hist2(j-1,:));
pi_tilde_9_hist2(j) = pi_tilde_9_func(hist2(j-1,:));
end

hists2(:,:,iters) = hist2;
e_hists2(:,iters) = e_hist2;
i_k_hists2(:,iters) = i_k_hist2;
f_hists2(:,iters) = f_hist2; 

RE_hists2(:,iters) =  RE_hist2;
pi_tilde_1_hists2(:,iters) =  pi_tilde_1_hist2;
pi_tilde_2_hists2(:,iters) =  pi_tilde_2_hist2;
pi_tilde_3_hists2(:,iters) =  pi_tilde_3_hist2;
pi_tilde_4_hists2(:,iters) =  pi_tilde_4_hist2;
pi_tilde_5_hists2(:,iters) =  pi_tilde_5_hist2;
pi_tilde_6_hists2(:,iters) =  pi_tilde_6_hist2;
pi_tilde_7_hists2(:,iters) =  pi_tilde_7_hist2;
pi_tilde_8_hists2(:,iters) =  pi_tilde_8_hist2;
pi_tilde_9_hists2(:,iters) =  pi_tilde_9_hist2;
end

save(filename2);
