%%%%% This file generates HJB results for the Ambiguity Neutral model.
% Authors: Mike Barnett, Jieyao Wang
% Last update: Dec 3,2019

close all
clear all
clc

%% Step 0: set up solver
addpath('/mnt/ide0/home/wangjieyao/Climate/FT/')
addpath('/home/wangjieyao/FT/')
addpath('/Volumes/homes/FT/')

%% Step 1: Set up parameters
McD = csvread('TCRE_MacDougallEtAl2017_update.csv');
McD = McD./1000;
par_lambda_McD = McD(:,1);
beta_f = mean(par_lambda_McD); 
var_beta_f = var(par_lambda_McD);
lambda = 1./var_beta_f;
delta = 0.01;
kappa= 0.032;
sigma_g = 0.02;
sigma_k = 0.0161;
sigma_r = 0.0339;
alpha = 0.115000000000000;
phi_0 = 0.0600;
phi_1 = 16.666666666666668;
mu_k = -0.034977443912449;
psi_0 = 0.112733407891680;
psi_1 = 0.142857142857143;
t_bar = 13;

% ambiguity parameter
xi_p = 1000;

%% Step 2: solve HJB
r_min = 0;
r_max = 9;
t_min = 0;
t_max = 750;
k_min = 0;
k_max = 18;

hr = 0.05;
ht = 25;
hk = 0.15;

r = r_min:hr:r_max;
t = t_min:ht:t_max;
k = k_min:hk:k_max;

[r_mat,F_mat,k_mat] = ndgrid(r,t,k); 

dt = 0.1;
tol = 1e-12;
v0 = (kappa).*r_mat+(1-kappa).*k_mat;
v1_initial = v0.*ones(size(r_mat));
out = v0;
vold = v0 .* ones(size(v0));
iter = 1;

mu_1 = 1.272e-02;
mu_2 = -4.871e-04;
sigma_1 = 3.248e-03;
sigma_2 = 1.029e-04;
rho_12 = -2.859133e-07;

gamma_1 = -mu_1;
gamma_2 = -mu_2*2;
sigma_2 = sigma_2*2;
rho_12 = rho_12*2;
mu = [gamma_1,gamma_2];

sigma = [sigma_1.^2, rho_12; rho_12, sigma_2.^2];
Sigma = [var_beta_f,0,0;0,sigma_1.^2, rho_12; 0, rho_12, sigma_2.^2];

% dmg
[gam1,w1] = quad_points_hermite(3) ;
gamm1 = sqrt(2) * 1 * gam1 + 0 ;
[gam2,w2]= quad_points_hermite(3);
gamm2 = sqrt(2) * 1 * gam2 + 0 ;
At=(chol(sigma))';
x(:,1) = [gamma_1;gamma_2] + At*[gamm1(1);gamm2(1)];
x(:,2) = [gamma_1;gamma_2] + At*[gamm1(1);gamm2(2)];
x(:,3) = [gamma_1;gamma_2] + At*[gamm1(1);gamm2(3)];
x(:,4) = [gamma_1;gamma_2] + At*[gamm1(2);gamm2(1)];
x(:,5) = [gamma_1;gamma_2] + At*[gamm1(2);gamm2(2)];
x(:,6) = [gamma_1;gamma_2] + At*[gamm1(2);gamm2(3)];
x(:,7) = [gamma_1;gamma_2] + At*[gamm1(3);gamm2(1)];
x(:,8) = [gamma_1;gamma_2] + At*[gamm1(3);gamm2(2)];
x(:,9) = [gamma_1;gamma_2] + At*[gamm1(3);gamm2(3)];
w(:,1) = [w1(1);w2(1)];
w(:,2) = [w1(1);w2(2)];
w(:,3) = [w1(1);w2(3)];
w(:,4) = [w1(2);w2(1)];
w(:,5) = [w1(2);w2(2)];
w(:,6) = [w1(2);w2(3)];
w(:,7) = [w1(3);w2(1)];
w(:,8) = [w1(3);w2(2)];
w(:,9) = [w1(3);w2(3)];
gamma1 = x(1,:);
gamma2 = x(2,:);
wgt1 = w(1,:);
wgt2 = w(2,:);
vals = linspace(0,30,100);
weight1 = wgt1(1).*wgt2(1);
weight2 = wgt1(2).*wgt2(2);
weight3 = wgt1(3).*wgt2(3);
weight4 = wgt1(4).*wgt2(4);
weight5 = wgt1(5).*wgt2(5);
weight6 = wgt1(6).*wgt2(6);
weight7 = wgt1(7).*wgt2(7);
weight8 = wgt1(8).*wgt2(8);
weight9 = wgt1(9).*wgt2(9);
weight_total = weight1+weight2+weight3+weight4+weight5+weight6+weight7+weight8+weight9;
weight1 = weight1./weight_total;
weight2 = weight2./weight_total;
weight3 = weight3./weight_total;
weight4 = weight4./weight_total;
weight5 = weight5./weight_total;
weight6 = weight6./weight_total;
weight7 = weight7./weight_total;
weight8 = weight8./weight_total;
weight9 = weight9./weight_total;
gamma0(1) = max(-(gamma1(1).*vals+0.5.*gamma2(1).*vals.^2));
gamma0(2) = max(-(gamma1(2).*vals+0.5.*gamma2(2).*vals.^2));
gamma0(3) = max(-(gamma1(3).*vals+0.5.*gamma2(3).*vals.^2));
gamma0(4) = max(-(gamma1(4).*vals+0.5.*gamma2(4).*vals.^2));
gamma0(5) = max(-(gamma1(5).*vals+0.5.*gamma2(5).*vals.^2));
gamma0(6) = max(-(gamma1(6).*vals+0.5.*gamma2(6).*vals.^2));
gamma0(7) = max(-(gamma1(7).*vals+0.5.*gamma2(7).*vals.^2));
gamma0(8) = max(-(gamma1(8).*vals+0.5.*gamma2(8).*vals.^2));
gamma0(9) = max(-(gamma1(9).*vals+0.5.*gamma2(9).*vals.^2));

v0 = v0.*ones(size(v0));
v0 = reshape(out,size(v0));

v0_dt = zeros(size(v0));
v0_dt(:,2:end-1,:) = (1./(2.*ht)).*(v0(:,3:end,:)-v0(:,1:end-2,:));
v0_dt(:,end,:) = (1./ht).*(v0(:,end,:)-v0(:,end-1,:));
v0_dt(:,1,:) = (1./ht).*(v0(:,2,:)-v0(:,1,:));

v0_dr = zeros(size(v0));
v0_dr(2:end-1,:,:) = (1./(2.*hr)).*(v0(3:end,:,:)-v0(1:end-2,:,:));
v0_dr(end,:,:) = (1./hr).*(v0(end,:,:)-v0(end-1,:,:));
v0_dr(1,:,:) = (1./hr).*(v0(2,:,:)-v0(1,:,:));

v0_dk = zeros(size(v0));
v0_dk(:,:,2:end-1) = (1./(2.*hk)).*(v0(:,:,3:end)-v0(:,:,1:end-2));
v0_dk(:,:,end) = (1./hk).*(v0(:,:,end)-v0(:,:,end-1));
v0_dk(:,:,1) = (1./hk).*(v0(:,:,2)-v0(:,:,1));

v0_drr = zeros(size(v0));
v0_drr(2:end-1,:,:) = (1./(hr.^2)).*(v0(3:end,:,:)+v0(1:end-2,:,:)-2.*v0(2:end-1,:,:));
v0_drr(end,:,:) = (1./(hr.^2)).*(v0(end,:,:)+v0(end-2,:,:)-2.*v0(end-1,:,:));
v0_drr(1,:,:) = (1./(hr.^2)).*(v0(3,:,:)+v0(1,:,:)-2.*v0(2,:,:));

v0_drr(v0_dr<1e-16) = 0;
v0_dr(v0_dr<1e-16) = 1e-16;

v0_dtt = zeros(size(v0));
v0_dtt(:,2:end-1) = (1./(ht.^2)).*(v0(:,3:end)+v0(:,1:end-2)-2.*v0(:,2:end-1));
v0_dtt(:,end) = (1./(ht.^2)).*(v0(:,end)+v0(:,end-2)-2.*v0(:,end-1));
v0_dtt(:,1) = (1./(ht.^2)).*(v0(:,3)+v0(:,1)-2.*v0(:,2));

v0_dkk = zeros(size(v0));
v0_dkk(:,:,2:end-1) = (1./(hk.^2)).*(v0(:,:,3:end)+v0(:,:,1:end-2)-2.*v0(:,:,2:end-1));
v0_dkk(:,:,end) = (1./(hk.^2)).*(v0(:,:,end)+v0(:,:,end-2)-2.*v0(:,:,end-1));
v0_dkk(:,:,1) = (1./(hk.^2)).*(v0(:,:,3)+v0(:,:,1)-2.*v0(:,:,2));
   
% FOC  
B1 = v0_dr-v0_dt.*exp(r_mat);  
C1 = delta.*kappa;
e = C1./B1;
e_hat = e;
e_star = e_hat;

Acoeff = ones(size(r_mat));
Bcoeff = ((delta.*(1-kappa).*phi_1+phi_0.*phi_1.*v0_dk).*delta.*(1-kappa)./(v0_dr.*psi_0.*0.5)...
    .*exp(0.5.*(r_mat-k_mat)))./(delta.*(1-kappa).*phi_1);
Ccoeff = -alpha - 1./phi_1;
j = ((-Bcoeff+sqrt(Bcoeff.^2 - 4.*Acoeff.*Ccoeff))./(2.*Acoeff)).^(2);
i_k = alpha-j-(delta.*(1-kappa))./(v0_dr.*psi_0.*0.5).*j.^0.5.*exp(0.5.*(r_mat-k_mat));
     
% % % % % %
a_1 = -v0_dk.*(gamma0(1)+gamma1(1).*t_bar+0.5.*gamma2(1).*(t_bar.^2));
b_1 = -v0_dk.*F_mat.*(gamma1(1)+gamma2(1).*t_bar);
c_1 = -v0_dk.*gamma2(1).*(F_mat.^2);
lambda_tilde_1 = lambda+c_1./xi_p;
beta_tilde_1 = beta_f-c_1./xi_p./lambda_tilde_1.*beta_f-1./xi_p./lambda_tilde_1.*b_1;
I_1 = a_1-0.5.*log(lambda).*xi_p + 0.5.*log(lambda_tilde_1).*xi_p ...
    +0.5.*lambda.*beta_f.^2.*xi_p -0.5.*lambda_tilde_1.*(beta_tilde_1).^2.*xi_p;
R_1 = 1./xi_p.*(I_1-(a_1+b_1.*beta_tilde_1+0.5.*c_1.*(beta_tilde_1).^2+0.5.*c_1./lambda_tilde_1));
pi_tilde_1 = weight1.*exp(-1./xi_p.*I_1);    
    
% % % % % %     
a_2 = -v0_dk.*(gamma0(2)+gamma1(2).*t_bar+0.5.*gamma2(2).*(t_bar.^2));
b_2 = -v0_dk.*F_mat.*(gamma1(2)+gamma2(2).*t_bar);
c_2 = -v0_dk.*gamma2(2).*(F_mat.^2);
lambda_tilde_2 = lambda+c_2./xi_p;
beta_tilde_2 = beta_f-c_2./xi_p./lambda_tilde_2.*beta_f-1./xi_p./lambda_tilde_2.*b_2;
I_2 = a_2-0.5.*log(lambda).*xi_p + 0.5.*log(lambda_tilde_2).*xi_p ...
    +0.5.*lambda.*beta_f.^2.*xi_p -0.5.*lambda_tilde_2.*(beta_tilde_2).^2.*xi_p;
R_2 = 1./xi_p.*(I_2-(a_2+b_2.*beta_tilde_2+0.5.*c_2.*(beta_tilde_2).^2+0.5.*c_2./lambda_tilde_2));
pi_tilde_2 = weight2.*exp(-1./xi_p.*I_2);   
    
% % % % % % 
a_3 = -v0_dk.*(gamma0(3)+gamma1(3).*t_bar+0.5.*gamma2(3).*(t_bar.^2));
b_3 = -v0_dk.*F_mat.*(gamma1(3)+gamma2(3).*t_bar);
c_3 = -v0_dk.*gamma2(3).*(F_mat.^2);
lambda_tilde_3 = lambda+c_3./xi_p;
beta_tilde_3 = beta_f-c_3./xi_p./lambda_tilde_3.*beta_f-1./xi_p./lambda_tilde_3.*b_3;
I_3 = a_3-0.5.*log(lambda).*xi_p + 0.5.*log(lambda_tilde_3).*xi_p ...
    +0.5.*lambda.*beta_f.^2.*xi_p -0.5.*lambda_tilde_3.*(beta_tilde_3).^2.*xi_p;
R_3 = 1./xi_p.*(I_3-(a_3+b_3.*beta_tilde_3+0.5.*c_3.*(beta_tilde_3).^2+0.5.*c_3./lambda_tilde_3));
pi_tilde_3 = weight3.*exp(-1./xi_p.*I_3);   
    
% % % % % % 
a_4 = -v0_dk.*(gamma0(4)+gamma1(4).*t_bar+0.5.*gamma2(4).*(t_bar.^2));
b_4 = -v0_dk.*F_mat.*(gamma1(4)+gamma2(4).*t_bar);
c_4 = -v0_dk.*gamma2(4).*(F_mat.^2);
lambda_tilde_4 = lambda+c_4./xi_p;
beta_tilde_4 = beta_f-c_4./xi_p./lambda_tilde_4.*beta_f-1./xi_p./lambda_tilde_4.*b_4;
I_4 = a_4-0.5.*log(lambda).*xi_p + 0.5.*log(lambda_tilde_4).*xi_p ...
    +0.5.*lambda.*beta_f.^2.*xi_p -0.5.*lambda_tilde_4.*(beta_tilde_4).^2.*xi_p;
R_4 = 1./xi_p.*(I_4-(a_4+b_4.*beta_tilde_4+0.5.*c_4.*(beta_tilde_4).^2+0.5.*c_4./lambda_tilde_4));
pi_tilde_4 = weight4.*exp(-1./xi_p.*I_4); 
    
% % % % % % 
a_5 = -v0_dk.*(gamma0(5)+gamma1(5).*t_bar+0.5.*gamma2(5).*(t_bar.^2));
b_5 = -v0_dk.*F_mat.*(gamma1(5)+gamma2(5).*t_bar);
c_5 = -v0_dk.*gamma2(5).*(F_mat.^2);
lambda_tilde_5 = lambda+c_5./xi_p;
beta_tilde_5 = beta_f-c_5./xi_p./lambda_tilde_5.*beta_f-1./xi_p./lambda_tilde_5.*b_5;
I_5 = a_5-0.5.*log(lambda).*xi_p + 0.5.*log(lambda_tilde_5).*xi_p ...
    +0.5.*lambda.*beta_f.^2.*xi_p -0.5.*lambda_tilde_5.*(beta_tilde_5).^2.*xi_p;
R_5 = 1./xi_p.*(I_5-(a_5+b_5.*beta_tilde_5+0.5.*c_5.*(beta_tilde_5).^2+0.5.*c_5./lambda_tilde_5));
pi_tilde_5 = weight5.*exp(-1./xi_p.*I_5); 
    
% % % % % % 
a_6 = -v0_dk.*(gamma0(6)+gamma1(6).*t_bar+0.5.*gamma2(6).*(t_bar.^2));
b_6 = -v0_dk.*F_mat.*(gamma1(6)+gamma2(6).*t_bar);
c_6 = -v0_dk.*gamma2(6).*(F_mat.^2);
lambda_tilde_6 = lambda+c_6./xi_p;
beta_tilde_6 = beta_f-c_6./xi_p./lambda_tilde_6.*beta_f-1./xi_p./lambda_tilde_6.*b_6;
I_6 = a_6-0.5.*log(lambda).*xi_p + 0.5.*log(lambda_tilde_6).*xi_p ...
    +0.5.*lambda.*beta_f.^2.*xi_p -0.5.*lambda_tilde_6.*(beta_tilde_6).^2.*xi_p;
R_6 = 1./xi_p.*(I_6-(a_6+b_6.*beta_tilde_6+0.5.*c_6.*(beta_tilde_6).^2+0.5.*c_6./lambda_tilde_6));
pi_tilde_6 = weight6.*exp(-1./xi_p.*I_6); 
    
% % % % % % 
a_7 = -v0_dk.*(gamma0(7)+gamma1(7).*t_bar+0.5.*gamma2(7).*(t_bar.^2));
b_7 = -v0_dk.*F_mat.*(gamma1(7)+gamma2(7).*t_bar);
c_7 = -v0_dk.*gamma2(7).*(F_mat.^2);
lambda_tilde_7 = lambda+c_7./xi_p;
beta_tilde_7 = beta_f-c_7./xi_p./lambda_tilde_7.*beta_f-1./xi_p./lambda_tilde_7.*b_7;
I_7 = a_7-0.5.*log(lambda).*xi_p + 0.5.*log(lambda_tilde_7).*xi_p ...
    +0.5.*lambda.*beta_f.^2.*xi_p -0.5.*lambda_tilde_7.*(beta_tilde_7).^2.*xi_p;
R_7 = 1./xi_p.*(I_7-(a_7+b_7.*beta_tilde_7+0.5.*c_7.*(beta_tilde_7).^2+0.5.*c_7./lambda_tilde_7));
pi_tilde_7 = weight7.*exp(-1./xi_p.*I_7); 
    
% % % % % % 
a_8 = -v0_dk.*(gamma0(8)+gamma1(8).*t_bar+0.5.*gamma2(8).*(t_bar.^2));
b_8 = -v0_dk.*F_mat.*(gamma1(8)+gamma2(8).*t_bar);
c_8 = -v0_dk.*gamma2(8).*(F_mat.^2);
lambda_tilde_8 = lambda+c_8./xi_p;
beta_tilde_8 = beta_f-c_8./xi_p./lambda_tilde_8.*beta_f-1./xi_p./lambda_tilde_8.*b_8;
I_8 = a_8-0.5.*log(lambda).*xi_p + 0.5.*log(lambda_tilde_8).*xi_p ...
    +0.5.*lambda.*beta_f.^2.*xi_p -0.5.*lambda_tilde_8.*(beta_tilde_8).^2.*xi_p;
R_8 = 1./xi_p.*(I_8-(a_8+b_8.*beta_tilde_8+0.5.*c_8.*(beta_tilde_8).^2+0.5.*c_8./lambda_tilde_8));
pi_tilde_8 = weight8.*exp(-1./xi_p.*I_8); 
    
% % % % % % 
a_9 = -v0_dk.*(gamma0(9)+gamma1(9).*t_bar+0.5.*gamma2(9).*(t_bar.^2));
b_9 = -v0_dk.*F_mat.*(gamma1(9)+gamma2(9).*t_bar);
c_9 = -v0_dk.*gamma2(9).*(F_mat.^2);
lambda_tilde_9 = lambda+c_9./xi_p;
beta_tilde_9 = beta_f-c_9./xi_p./lambda_tilde_9.*beta_f-1./xi_p./lambda_tilde_9.*b_9;
I_9 = a_9-0.5.*log(lambda).*xi_p + 0.5.*log(lambda_tilde_9).*xi_p ...
    +0.5.*lambda.*beta_f.^2.*xi_p -0.5.*lambda_tilde_9.*(beta_tilde_9).^2.*xi_p;
R_9 = 1./xi_p.*(I_9-(a_9+b_9.*beta_tilde_9+0.5.*c_9.*(beta_tilde_9).^2+0.5.*c_9./lambda_tilde_9));
pi_tilde_9 = weight9.*exp(-1./xi_p.*I_9); 
    
% % % % % %
pi_tilde_total = pi_tilde_1+pi_tilde_2+pi_tilde_3+pi_tilde_4...
    +pi_tilde_5+pi_tilde_6+pi_tilde_7+pi_tilde_8+pi_tilde_9;
    
pi_tilde_1_norm = pi_tilde_1./pi_tilde_total;
pi_tilde_2_norm = pi_tilde_2./pi_tilde_total;
pi_tilde_3_norm = pi_tilde_3./pi_tilde_total;
pi_tilde_4_norm = pi_tilde_4./pi_tilde_total;
pi_tilde_5_norm = pi_tilde_5./pi_tilde_total;
pi_tilde_6_norm = pi_tilde_6./pi_tilde_total;
pi_tilde_7_norm = pi_tilde_7./pi_tilde_total;
pi_tilde_8_norm = pi_tilde_8./pi_tilde_total;
pi_tilde_9_norm = pi_tilde_9./pi_tilde_total;

J_1 = a_1+b_1.*beta_tilde_1+0.5.*c_1.*(beta_tilde_1.^2)+0.5.*c_1./lambda_tilde_1;    
J_2 = a_2+b_2.*beta_tilde_2+0.5.*c_2.*(beta_tilde_2.^2)+0.5.*c_2./lambda_tilde_2;    
J_3 = a_3+b_3.*beta_tilde_3+0.5.*c_3.*(beta_tilde_3.^2)+0.5.*c_3./lambda_tilde_3;    
J_4 = a_4+b_4.*beta_tilde_4+0.5.*c_4.*(beta_tilde_4.^2)+0.5.*c_4./lambda_tilde_4;    
J_5 = a_5+b_5.*beta_tilde_5+0.5.*c_5.*(beta_tilde_5.^2)+0.5.*c_5./lambda_tilde_5;    
J_6 = a_6+b_6.*beta_tilde_6+0.5.*c_6.*(beta_tilde_6.^2)+0.5.*c_6./lambda_tilde_6;    
J_7 = a_7+b_7.*beta_tilde_7+0.5.*c_7.*(beta_tilde_7.^2)+0.5.*c_7./lambda_tilde_7;    
J_8 = a_8+b_8.*beta_tilde_8+0.5.*c_8.*(beta_tilde_8.^2)+0.5.*c_8./lambda_tilde_8;    
J_9 = a_9+b_9.*beta_tilde_9+0.5.*c_9.*(beta_tilde_9.^2)+0.5.*c_9./lambda_tilde_9;    

I_term = -1.*xi_p.*log(pi_tilde_1+pi_tilde_2+pi_tilde_3+pi_tilde_4...
    +pi_tilde_5+pi_tilde_6+pi_tilde_7+pi_tilde_8+pi_tilde_9);

R_1 = 1./xi_p.*(I_1 - J_1);
R_2 = 1./xi_p.*(I_2 - J_2);
R_3 = 1./xi_p.*(I_3 - J_3);
R_4 = 1./xi_p.*(I_4 - J_4);
R_5 = 1./xi_p.*(I_5 - J_5);
R_6 = 1./xi_p.*(I_6 - J_6);
R_7 = 1./xi_p.*(I_7 - J_7);
R_8 = 1./xi_p.*(I_8 - J_8);
R_9 = 1./xi_p.*(I_9 - J_9);

drift_distort = (pi_tilde_1_norm.*J_1 ...
    +pi_tilde_2_norm.*J_2...
    +pi_tilde_3_norm.*J_3...
    +pi_tilde_4_norm.*J_4...
    +pi_tilde_5_norm.*J_5...
    +pi_tilde_6_norm.*J_6...
    +pi_tilde_7_norm.*J_7...
    +pi_tilde_8_norm.*J_8...
    +pi_tilde_9_norm.*J_9);
    
RE = pi_tilde_1_norm.*R_1 + pi_tilde_2_norm.*R_2 + pi_tilde_3_norm.*R_3+ pi_tilde_4_norm.*R_4 ...
    +pi_tilde_5_norm.*R_5 + pi_tilde_6_norm.*R_6 + pi_tilde_7_norm.*R_7+ pi_tilde_8_norm.*R_8+ pi_tilde_9_norm.*R_9 ... ...
    +pi_tilde_1_norm.*(log(pi_tilde_1_norm./weight1)) ...
    +pi_tilde_2_norm.*(log(pi_tilde_2_norm./weight2))...
    +pi_tilde_3_norm.*(log(pi_tilde_3_norm./weight3))...
    +pi_tilde_4_norm.*(log(pi_tilde_4_norm./weight4))...
    +pi_tilde_5_norm.*(log(pi_tilde_5_norm./weight5)) ...
    +pi_tilde_6_norm.*(log(pi_tilde_6_norm./weight6))...
    +pi_tilde_7_norm.*(log(pi_tilde_7_norm./weight7))...
    +pi_tilde_8_norm.*(log(pi_tilde_8_norm./weight8))...
    +pi_tilde_9_norm.*(log(pi_tilde_9_norm./weight9));
RE_total = 1.*xi_p.*RE;

% inputs for solver
A = -delta.*ones(size(r_mat));
B_r = -e_star+psi_0.*(j.^psi_1).*exp(psi_1.*(k_mat-r_mat))-0.5.*(sigma_r.^2);
B_k = mu_k+phi_0.*log(1+i_k.*phi_1)-0.5.*(sigma_k.^2);
B_t = e_star.*exp(r_mat);
C_rr = 0.5.*sigma_r.^2.*ones(size(r_mat));
C_kk = 0.5.*sigma_k.^2.*ones(size(r_mat));
C_tt = zeros(size(r_mat));
D = delta.*kappa.*log(e_star)+delta.*kappa.*r_mat ...
    + delta.*(1-kappa).*(log(alpha-i_k-j)+(k_mat)) ...
    +I_term;

stateSpace = [r_mat(:), F_mat(:), k_mat(:)]; 
model      = {};
model.A    = A(:); 
model.B    = [B_r(:), B_t(:), B_k(:)];
model.C    = [C_rr(:), C_tt(:), C_kk(:)];
model.D    = D(:);
model.v0   = v0(:);
model.dt   = dt;

out = solveCGNatural_1e16(stateSpace, model);
out_comp = reshape(out,size(v0)).*ones(size(r_mat));
    
disp(['Error: ', num2str(max(max(max(abs(out_comp - v1_initial)))))])
v0 = v0.*ones(size(v0));
v0 = reshape(out,size(v0));
q = delta.*(1-kappa)./(alpha -i_k-j);
eta = 0.05;

lhs_error = max(max(max(abs(out_comp - v1_initial))));
lhs_err(iter) = lhs_error;

while (max(max(max(max(abs(out_comp - vold)))))) > tol % check for convergence
    tic
    vold = v0 .* ones(size(v0));
    
    stateSpace = [r_mat(:), F_mat(:), k_mat(:)]; 
    model      = {};
    model.A    = A(:); 
    model.B    = [B_r(:), B_t(:), B_k(:)];
    model.C    = [C_rr(:), C_tt(:), C_kk(:)];
    model.D    = D(:);
    model.v0   = v0(:);
    model.dt   = dt;

    out = solveCGNatural_1e16(stateSpace, model);
    out_comp = reshape(out,size(v0)).*ones(size(r_mat));
    
    iter = iter+1;
    disp(['Diff: ', num2str(max(max(max(abs(out_comp - vold))))),' Number of Iters: ',num2str(iter)])

    v0 = v0.*ones(size(v0));
    v0 = reshape(out,size(v0));

    v0_dt = zeros(size(v0));
    v0_dt(:,2:end-1,:) = (1./(2.*ht)).*(v0(:,3:end,:)-v0(:,1:end-2,:));
    v0_dt(:,end,:) = (1./ht).*(v0(:,end,:)-v0(:,end-1,:));
    v0_dt(:,1,:) = (1./ht).*(v0(:,2,:)-v0(:,1,:));
    
    v0_dr = zeros(size(v0));
    v0_dr(2:end-1,:,:) = (1./(2.*hr)).*(v0(3:end,:,:)-v0(1:end-2,:,:));
    v0_dr(end,:,:) = (1./hr).*(v0(end,:,:)-v0(end-1,:,:));
    v0_dr(1,:,:) = (1./hr).*(v0(2,:,:)-v0(1,:,:));

    v0_dk = zeros(size(v0));
    v0_dk(:,:,2:end-1) = (1./(2.*hk)).*(v0(:,:,3:end)-v0(:,:,1:end-2));
    v0_dk(:,:,end) = (1./hk).*(v0(:,:,end)-v0(:,:,end-1));
    v0_dk(:,:,1) = (1./hk).*(v0(:,:,2)-v0(:,:,1));
       
    v0_drr = zeros(size(v0));
    v0_drr(2:end-1,:,:) = (1./(hr.^2)).*(v0(3:end,:,:)+v0(1:end-2,:,:)-2.*v0(2:end-1,:,:));
    v0_drr(end,:,:) = (1./(hr.^2)).*(v0(end,:,:)+v0(end-2,:,:)-2.*v0(end-1,:,:));
    v0_drr(1,:,:) = (1./(hr.^2)).*(v0(3,:,:)+v0(1,:,:)-2.*v0(2,:,:));
    
    v0_drr(v0_dr<1e-16) = 0;
    v0_dr(v0_dr<1e-16) = 1e-16;
    
    v0_dtt = zeros(size(v0));
    v0_dtt(:,2:end-1) = (1./(ht.^2)).*(v0(:,3:end)+v0(:,1:end-2)-2.*v0(:,2:end-1));
    v0_dtt(:,end) = (1./(ht.^2)).*(v0(:,end)+v0(:,end-2)-2.*v0(:,end-1));
    v0_dtt(:,1) = (1./(ht.^2)).*(v0(:,3)+v0(:,1)-2.*v0(:,2));

    v0_dkk = zeros(size(v0));
    v0_dkk(:,:,2:end-1) = (1./(hk.^2)).*(v0(:,:,3:end)+v0(:,:,1:end-2)-2.*v0(:,:,2:end-1));
    v0_dkk(:,:,end) = (1./(hk.^2)).*(v0(:,:,end)+v0(:,:,end-2)-2.*v0(:,:,end-1));
    v0_dkk(:,:,1) = (1./(hk.^2)).*(v0(:,:,3)+v0(:,:,1)-2.*v0(:,:,2));  
    
    e_hat = e_star;
    B1 = v0_dr-v0_dt.*exp(r_mat);  
    C1 = delta.*kappa;
    e = C1./B1;
    e_star = e;

    Converged = 0;
    nums = 0;
    while Converged == 0
        istar = (phi_0.*phi_1.*v0_dk./q - 1)./phi_1;
        jstar = (q.*exp(psi_1.*(r_mat-k_mat))./((v0_dr).*psi_0.*psi_1)).^(1./(psi_1 - 1));
        if alpha > (istar+jstar)
            qstar = eta.*delta.*(1-kappa)./(alpha-istar-jstar)+(1-eta).*q;
        else
            qstar = 2.*q;
        end

        if (max(max(max(max(abs((qstar-q)./eta)))))<=1e-10)
            Converged = 1; 
        end
        q = qstar;
        i_k = istar;
        j = jstar;
        nums = nums+1;
    end

    % % % % % %
    a_1 = -v0_dk.*(gamma0(1)+gamma1(1).*t_bar+0.5.*gamma2(1).*(t_bar.^2));
    b_1 = -v0_dk.*F_mat.*(gamma1(1)+gamma2(1).*t_bar);
    c_1 = -v0_dk.*gamma2(1).*(F_mat.^2);
    lambda_tilde_1 = lambda+c_1./xi_p;
    beta_tilde_1 = beta_f-c_1./xi_p./lambda_tilde_1.*beta_f-1./xi_p./lambda_tilde_1.*b_1;
    I_1 = a_1-0.5.*log(lambda).*xi_p + 0.5.*log(lambda_tilde_1).*xi_p ...
        +0.5.*lambda.*beta_f.^2.*xi_p -0.5.*lambda_tilde_1.*(beta_tilde_1).^2.*xi_p;
    R_1 = 1./xi_p.*(I_1-(a_1+b_1.*beta_tilde_1+0.5.*c_1.*(beta_tilde_1).^2+0.5.*c_1./lambda_tilde_1));
    pi_tilde_1 = weight1.*exp(-1./xi_p.*I_1);    

    % % % % % %     
    a_2 = -v0_dk.*(gamma0(2)+gamma1(2).*t_bar+0.5.*gamma2(2).*(t_bar.^2));
    b_2 = -v0_dk.*F_mat.*(gamma1(2)+gamma2(2).*t_bar);
    c_2 = -v0_dk.*gamma2(2).*(F_mat.^2);
    lambda_tilde_2 = lambda+c_2./xi_p;
    beta_tilde_2 = beta_f-c_2./xi_p./lambda_tilde_2.*beta_f-1./xi_p./lambda_tilde_2.*b_2;
    I_2 = a_2-0.5.*log(lambda).*xi_p + 0.5.*log(lambda_tilde_2).*xi_p ...
        +0.5.*lambda.*beta_f.^2.*xi_p -0.5.*lambda_tilde_2.*(beta_tilde_2).^2.*xi_p;
    R_2 = 1./xi_p.*(I_2-(a_2+b_2.*beta_tilde_2+0.5.*c_2.*(beta_tilde_2).^2+0.5.*c_2./lambda_tilde_2));
    pi_tilde_2 = weight2.*exp(-1./xi_p.*I_2);   

    % % % % % % 
    a_3 = -v0_dk.*(gamma0(3)+gamma1(3).*t_bar+0.5.*gamma2(3).*(t_bar.^2));
    b_3 = -v0_dk.*F_mat.*(gamma1(3)+gamma2(3).*t_bar);
    c_3 = -v0_dk.*gamma2(3).*(F_mat.^2);
    lambda_tilde_3 = lambda+c_3./xi_p;
    beta_tilde_3 = beta_f-c_3./xi_p./lambda_tilde_3.*beta_f-1./xi_p./lambda_tilde_3.*b_3;
    I_3 = a_3-0.5.*log(lambda).*xi_p + 0.5.*log(lambda_tilde_3).*xi_p ...
        +0.5.*lambda.*beta_f.^2.*xi_p -0.5.*lambda_tilde_3.*(beta_tilde_3).^2.*xi_p;
    R_3 = 1./xi_p.*(I_3-(a_3+b_3.*beta_tilde_3+0.5.*c_3.*(beta_tilde_3).^2+0.5.*c_3./lambda_tilde_3));
    pi_tilde_3 = weight3.*exp(-1./xi_p.*I_3);   

    % % % % % % 
    a_4 = -v0_dk.*(gamma0(4)+gamma1(4).*t_bar+0.5.*gamma2(4).*(t_bar.^2));
    b_4 = -v0_dk.*F_mat.*(gamma1(4)+gamma2(4).*t_bar);
    c_4 = -v0_dk.*gamma2(4).*(F_mat.^2);
    lambda_tilde_4 = lambda+c_4./xi_p;
    beta_tilde_4 = beta_f-c_4./xi_p./lambda_tilde_4.*beta_f-1./xi_p./lambda_tilde_4.*b_4;
    I_4 = a_4-0.5.*log(lambda).*xi_p + 0.5.*log(lambda_tilde_4).*xi_p ...
        +0.5.*lambda.*beta_f.^2.*xi_p -0.5.*lambda_tilde_4.*(beta_tilde_4).^2.*xi_p;
    R_4 = 1./xi_p.*(I_4-(a_4+b_4.*beta_tilde_4+0.5.*c_4.*(beta_tilde_4).^2+0.5.*c_4./lambda_tilde_4));
    pi_tilde_4 = weight4.*exp(-1./xi_p.*I_4); 

    % % % % % % 
    a_5 = -v0_dk.*(gamma0(5)+gamma1(5).*t_bar+0.5.*gamma2(5).*(t_bar.^2));
    b_5 = -v0_dk.*F_mat.*(gamma1(5)+gamma2(5).*t_bar);
    c_5 = -v0_dk.*gamma2(5).*(F_mat.^2);
    lambda_tilde_5 = lambda+c_5./xi_p;
    beta_tilde_5 = beta_f-c_5./xi_p./lambda_tilde_5.*beta_f-1./xi_p./lambda_tilde_5.*b_5;
    I_5 = a_5-0.5.*log(lambda).*xi_p + 0.5.*log(lambda_tilde_5).*xi_p ...
        +0.5.*lambda.*beta_f.^2.*xi_p -0.5.*lambda_tilde_5.*(beta_tilde_5).^2.*xi_p;
    R_5 = 1./xi_p.*(I_5-(a_5+b_5.*beta_tilde_5+0.5.*c_5.*(beta_tilde_5).^2+0.5.*c_5./lambda_tilde_5));
    pi_tilde_5 = weight5.*exp(-1./xi_p.*I_5); 

    % % % % % % 
    a_6 = -v0_dk.*(gamma0(6)+gamma1(6).*t_bar+0.5.*gamma2(6).*(t_bar.^2));
    b_6 = -v0_dk.*F_mat.*(gamma1(6)+gamma2(6).*t_bar);
    c_6 = -v0_dk.*gamma2(6).*(F_mat.^2);
    lambda_tilde_6 = lambda+c_6./xi_p;
    beta_tilde_6 = beta_f-c_6./xi_p./lambda_tilde_6.*beta_f-1./xi_p./lambda_tilde_6.*b_6;
    I_6 = a_6-0.5.*log(lambda).*xi_p + 0.5.*log(lambda_tilde_6).*xi_p ...
        +0.5.*lambda.*beta_f.^2.*xi_p -0.5.*lambda_tilde_6.*(beta_tilde_6).^2.*xi_p;
    R_6 = 1./xi_p.*(I_6-(a_6+b_6.*beta_tilde_6+0.5.*c_6.*(beta_tilde_6).^2+0.5.*c_6./lambda_tilde_6));
    pi_tilde_6 = weight6.*exp(-1./xi_p.*I_6); 

    % % % % % % 
    a_7 = -v0_dk.*(gamma0(7)+gamma1(7).*t_bar+0.5.*gamma2(7).*(t_bar.^2));
    b_7 = -v0_dk.*F_mat.*(gamma1(7)+gamma2(7).*t_bar);
    c_7 = -v0_dk.*gamma2(7).*(F_mat.^2);
    lambda_tilde_7 = lambda+c_7./xi_p;
    beta_tilde_7 = beta_f-c_7./xi_p./lambda_tilde_7.*beta_f-1./xi_p./lambda_tilde_7.*b_7;
    I_7 = a_7-0.5.*log(lambda).*xi_p + 0.5.*log(lambda_tilde_7).*xi_p ...
        +0.5.*lambda.*beta_f.^2.*xi_p -0.5.*lambda_tilde_7.*(beta_tilde_7).^2.*xi_p;
    R_7 = 1./xi_p.*(I_7-(a_7+b_7.*beta_tilde_7+0.5.*c_7.*(beta_tilde_7).^2+0.5.*c_7./lambda_tilde_7));
    pi_tilde_7 = weight7.*exp(-1./xi_p.*I_7); 

    % % % % % % 
    a_8 = -v0_dk.*(gamma0(8)+gamma1(8).*t_bar+0.5.*gamma2(8).*(t_bar.^2));
    b_8 = -v0_dk.*F_mat.*(gamma1(8)+gamma2(8).*t_bar);
    c_8 = -v0_dk.*gamma2(8).*(F_mat.^2);
    lambda_tilde_8 = lambda+c_8./xi_p;
    beta_tilde_8 = beta_f-c_8./xi_p./lambda_tilde_8.*beta_f-1./xi_p./lambda_tilde_8.*b_8;
    I_8 = a_8-0.5.*log(lambda).*xi_p + 0.5.*log(lambda_tilde_8).*xi_p ...
        +0.5.*lambda.*beta_f.^2.*xi_p -0.5.*lambda_tilde_8.*(beta_tilde_8).^2.*xi_p;
    R_8 = 1./xi_p.*(I_8-(a_8+b_8.*beta_tilde_8+0.5.*c_8.*(beta_tilde_8).^2+0.5.*c_8./lambda_tilde_8));
    pi_tilde_8 = weight8.*exp(-1./xi_p.*I_8); 

    % % % % % % 
    a_9 = -v0_dk.*(gamma0(9)+gamma1(9).*t_bar+0.5.*gamma2(9).*(t_bar.^2));
    b_9 = -v0_dk.*F_mat.*(gamma1(9)+gamma2(9).*t_bar);
    c_9 = -v0_dk.*gamma2(9).*(F_mat.^2);
    lambda_tilde_9 = lambda+c_9./xi_p;
    beta_tilde_9 = beta_f-c_9./xi_p./lambda_tilde_9.*beta_f-1./xi_p./lambda_tilde_9.*b_9;
    I_9 = a_9-0.5.*log(lambda).*xi_p + 0.5.*log(lambda_tilde_9).*xi_p ...
        +0.5.*lambda.*beta_f.^2.*xi_p -0.5.*lambda_tilde_9.*(beta_tilde_9).^2.*xi_p;
    R_9 = 1./xi_p.*(I_9-(a_9+b_9.*beta_tilde_9+0.5.*c_9.*(beta_tilde_9).^2+0.5.*c_9./lambda_tilde_9));
    pi_tilde_9 = weight9.*exp(-1./xi_p.*I_9); 

    % % % % % %
    pi_tilde_total = pi_tilde_1+pi_tilde_2+pi_tilde_3+pi_tilde_4...
        +pi_tilde_5+pi_tilde_6+pi_tilde_7+pi_tilde_8+pi_tilde_9;

    pi_tilde_1_norm = pi_tilde_1./pi_tilde_total;
    pi_tilde_2_norm = pi_tilde_2./pi_tilde_total;
    pi_tilde_3_norm = pi_tilde_3./pi_tilde_total;
    pi_tilde_4_norm = pi_tilde_4./pi_tilde_total;
    pi_tilde_5_norm = pi_tilde_5./pi_tilde_total;
    pi_tilde_6_norm = pi_tilde_6./pi_tilde_total;
    pi_tilde_7_norm = pi_tilde_7./pi_tilde_total;
    pi_tilde_8_norm = pi_tilde_8./pi_tilde_total;
    pi_tilde_9_norm = pi_tilde_9./pi_tilde_total;

    J_1 = a_1+b_1.*beta_tilde_1+0.5.*c_1.*(beta_tilde_1.^2)+0.5.*c_1./lambda_tilde_1;    
    J_2 = a_2+b_2.*beta_tilde_2+0.5.*c_2.*(beta_tilde_2.^2)+0.5.*c_2./lambda_tilde_2;    
    J_3 = a_3+b_3.*beta_tilde_3+0.5.*c_3.*(beta_tilde_3.^2)+0.5.*c_3./lambda_tilde_3;    
    J_4 = a_4+b_4.*beta_tilde_4+0.5.*c_4.*(beta_tilde_4.^2)+0.5.*c_4./lambda_tilde_4;    
    J_5 = a_5+b_5.*beta_tilde_5+0.5.*c_5.*(beta_tilde_5.^2)+0.5.*c_5./lambda_tilde_5;    
    J_6 = a_6+b_6.*beta_tilde_6+0.5.*c_6.*(beta_tilde_6.^2)+0.5.*c_6./lambda_tilde_6;    
    J_7 = a_7+b_7.*beta_tilde_7+0.5.*c_7.*(beta_tilde_7.^2)+0.5.*c_7./lambda_tilde_7;    
    J_8 = a_8+b_8.*beta_tilde_8+0.5.*c_8.*(beta_tilde_8.^2)+0.5.*c_8./lambda_tilde_8;    
    J_9 = a_9+b_9.*beta_tilde_9+0.5.*c_9.*(beta_tilde_9.^2)+0.5.*c_9./lambda_tilde_9;    

    I_term = -1.*xi_p.*log(pi_tilde_1+pi_tilde_2+pi_tilde_3+pi_tilde_4...
        +pi_tilde_5+pi_tilde_6+pi_tilde_7+pi_tilde_8+pi_tilde_9);

    R_1 = 1./xi_p.*(I_1 - J_1);
    R_2 = 1./xi_p.*(I_2 - J_2);
    R_3 = 1./xi_p.*(I_3 - J_3);
    R_4 = 1./xi_p.*(I_4 - J_4);
    R_5 = 1./xi_p.*(I_5 - J_5);
    R_6 = 1./xi_p.*(I_6 - J_6);
    R_7 = 1./xi_p.*(I_7 - J_7);
    R_8 = 1./xi_p.*(I_8 - J_8);
    R_9 = 1./xi_p.*(I_9 - J_9);

    drift_distort = (pi_tilde_1_norm.*J_1 ...
        +pi_tilde_2_norm.*J_2...
        +pi_tilde_3_norm.*J_3...
        +pi_tilde_4_norm.*J_4...
        +pi_tilde_5_norm.*J_5...
        +pi_tilde_6_norm.*J_6...
        +pi_tilde_7_norm.*J_7...
        +pi_tilde_8_norm.*J_8...
        +pi_tilde_9_norm.*J_9);

    RE = pi_tilde_1_norm.*R_1 + pi_tilde_2_norm.*R_2 + pi_tilde_3_norm.*R_3+ pi_tilde_4_norm.*R_4 ...
        +pi_tilde_5_norm.*R_5 + pi_tilde_6_norm.*R_6 + pi_tilde_7_norm.*R_7+ pi_tilde_8_norm.*R_8+ pi_tilde_9_norm.*R_9 ... ...
        +pi_tilde_1_norm.*(log(pi_tilde_1_norm./weight1)) ...
        +pi_tilde_2_norm.*(log(pi_tilde_2_norm./weight2))...
        +pi_tilde_3_norm.*(log(pi_tilde_3_norm./weight3))...
        +pi_tilde_4_norm.*(log(pi_tilde_4_norm./weight4))...
        +pi_tilde_5_norm.*(log(pi_tilde_5_norm./weight5)) ...
        +pi_tilde_6_norm.*(log(pi_tilde_6_norm./weight6))...
        +pi_tilde_7_norm.*(log(pi_tilde_7_norm./weight7))...
        +pi_tilde_8_norm.*(log(pi_tilde_8_norm./weight8))...
        +pi_tilde_9_norm.*(log(pi_tilde_9_norm./weight9));
    RE_total = 1.*xi_p.*RE;

    % inputs for solver
    A = -delta.*ones(size(r_mat));
    B_r = -e_star+psi_0.*(j.^psi_1).*exp(psi_1.*(k_mat-r_mat))-0.5.*(sigma_r.^2);
    B_k = mu_k+phi_0.*log(1+i_k.*phi_1)-0.5.*(sigma_k.^2);
    B_t = e_star.*exp(r_mat);
    C_rr = 0.5.*sigma_r.^2.*ones(size(r_mat));
    C_kk = 0.5.*sigma_k.^2.*ones(size(r_mat));
    C_tt = zeros(size(r_mat));
    D = delta.*kappa.*log(e_star)+delta.*kappa.*r_mat ...
        + delta.*(1-kappa).*(log(alpha-i_k-j)+(k_mat)) ...
        +I_term;
    
    % check equation
    pde_error_new = A.*v0+B_r.*v0_dr+B_t.*v0_dt+B_k.*v0_dk+C_rr.*v0_drr+C_kk.*v0_dkk+C_tt.*v0_dtt+D;
    pde_error = pde_error_new;
    disp(['PDE Error: ', num2str(max(max(max(abs(pde_error)))))])
    
    % check convergence
    lhs_error = max(max(max(abs(out_comp - vold))));
    lhs_err(iter) = lhs_error;
    
    if lhs_error < tol
        fprintf('PDE Converges');
    end

    toc
end

s0 = '/HJB_Growth_Neutral';
filename = [pwd, s0];
save(filename);
