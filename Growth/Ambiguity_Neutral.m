%%%%% This file generates results for the Ambiguity Neutral model.
% Authors: Mike Barnett, Jieyao Wang
% Last update: Oct 07,2019

close all
clear all
clc

%% Step 0: set up solver
%%%%% change to location of the solver
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
alpha= 0.032;
sigma_g = 0.02;
sigma_k = 0.0161;
sigma_r = 0.0339;
A_O = 0.115000000000000;
Gamma = 0.0600;
Theta = 1./16.666666666666668;
Alpha = -0.034977443912449;
Gamma_r = 0.112733407891680;
Theta_r = 0.142857142857143;
t_bar = 13;

% ambiguity parameter
theta = 0.001; 

%% Step 2: solve HJB
r_min = 0;
r_max = 9;
t_min = 0;
t_max = 750;
k_min = 0;
k_max = 9;
nr = 30;
nt = 30;
nk = 25;
r = linspace(r_min,r_max,nr)';
t = linspace(t_min,t_max,nt)';
k = linspace(k_min,k_max,nk)';
hr = r(2) - r(1);
hk = k(2) - k(1);
ht = t(2) - t(1);
[r_mat,t_mat,k_mat] = ndgrid(r,t,k); 

dt = 0.5;
tol = 1e-16;
v0 = (alpha).*r_mat+(1-alpha).*k_mat;
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
v0_dr(v0_dr<1e-8) = 1e-8;

v0_dk = zeros(size(v0));
v0_dk(:,:,2:end-1) = (1./(2.*hk)).*(v0(:,:,3:end)-v0(:,:,1:end-2));
v0_dk(:,:,end) = (1./hk).*(v0(:,:,end)-v0(:,:,end-1));
v0_dk(:,:,1) = (1./hk).*(v0(:,:,2)-v0(:,:,1));

v0_drr = zeros(size(v0));
v0_drr(2:end-1,:,:) = (1./(hr.^2)).*(v0(3:end,:,:)+v0(1:end-2,:,:)-2.*v0(2:end-1,:,:));
v0_drr(end,:,:) = (1./(hr.^2)).*(v0(end,:,:)+v0(end-2,:,:)-2.*v0(end-1,:,:));
v0_drr(1,:,:) = (1./(hr.^2)).*(v0(3,:,:)+v0(1,:,:)-2.*v0(2,:,:));

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
C1 = delta.*alpha;
e = C1./B1;
e_hat = e;
e_star = e_hat;

Acoeff = ones(size(r_mat));
Bcoeff = ((delta.*(1-alpha)./Theta+Gamma./Theta.*v0_dk).*delta.*(1-alpha)./(v0_dr.*Gamma_r.*0.5)...
    .*exp(0.5.*(r_mat-k_mat)))./(delta.*(1-alpha)./Theta);
Ccoeff = -A_O - Theta;
f = ((-Bcoeff+sqrt(Bcoeff.^2 - 4.*Acoeff.*Ccoeff))./(2.*Acoeff)).^(2);
i_k = A_O-f-(delta.*(1-alpha))./(v0_dr.*Gamma_r.*0.5).*f.^0.5.*exp(0.5.*(r_mat-k_mat));
     
% % % % % %
a_1 = -v0_dk.*(gamma0(1)+gamma1(1).*t_bar+0.5.*gamma2(1).*(t_bar.^2));
b_1 = -v0_dk.*t_mat.*(gamma1(1)+gamma2(1).*t_bar);
c_1 = -v0_dk.*gamma2(1).*(t_mat.^2);
lambda_tilde_1 = lambda+c_1.*theta;
beta_tilde_1 = beta_f-c_1.*theta./lambda_tilde_1.*beta_f-theta./lambda_tilde_1.*b_1;
I_1 = a_1-0.5.*log(lambda)./theta + 0.5.*log(lambda_tilde_1)./theta ...
    +0.5.*lambda.*beta_f.^2./theta -0.5.*lambda_tilde_1.*(beta_tilde_1).^2./theta;
R_1 = theta.*(I_1-(a_1+b_1.*beta_tilde_1+0.5.*c_1.*(beta_tilde_1).^2+0.5.*c_1./lambda_tilde_1));
pi_tilde_1 = weight1.*exp(-theta.*I_1);    
    
% % % % % %     
a_2 = -v0_dk.*(gamma0(2)+gamma1(2).*t_bar+0.5.*gamma2(2).*(t_bar.^2));
b_2 = -v0_dk.*t_mat.*(gamma1(2)+gamma2(2).*t_bar);
c_2 = -v0_dk.*gamma2(2).*(t_mat.^2);
lambda_tilde_2 = lambda+c_2.*theta;
beta_tilde_2 = beta_f-c_2.*theta./lambda_tilde_2.*beta_f-theta./lambda_tilde_2.*b_2;
I_2 = a_2-0.5.*log(lambda)./theta + 0.5.*log(lambda_tilde_2)./theta ...
    +0.5.*lambda.*beta_f.^2./theta -0.5.*lambda_tilde_2.*(beta_tilde_2).^2./theta;
R_2 = theta.*(I_2-(a_2+b_2.*beta_tilde_2+0.5.*c_2.*(beta_tilde_2).^2+0.5.*c_2./lambda_tilde_2));
pi_tilde_2 = weight2.*exp(-theta.*I_2);   
    
% % % % % % 
a_3 = -v0_dk.*(gamma0(3)+gamma1(3).*t_bar+0.5.*gamma2(3).*(t_bar.^2));
b_3 = -v0_dk.*t_mat.*(gamma1(3)+gamma2(3).*t_bar);
c_3 = -v0_dk.*gamma2(3).*(t_mat.^2);
lambda_tilde_3 = lambda+c_3.*theta;
beta_tilde_3 = beta_f-c_3.*theta./lambda_tilde_3.*beta_f-theta./lambda_tilde_3.*b_3;
I_3 = a_3-0.5.*log(lambda)./theta + 0.5.*log(lambda_tilde_3)./theta ...
    +0.5.*lambda.*beta_f.^2./theta -0.5.*lambda_tilde_3.*(beta_tilde_3).^2./theta;
R_3 = theta.*(I_3-(a_3+b_3.*beta_tilde_3+0.5.*c_3.*(beta_tilde_3).^2+0.5.*c_3./lambda_tilde_3));
pi_tilde_3 = weight3.*exp(-theta.*I_3);   
    
% % % % % % 
a_4 = -v0_dk.*(gamma0(4)+gamma1(4).*t_bar+0.5.*gamma2(4).*(t_bar.^2));
b_4 = -v0_dk.*t_mat.*(gamma1(4)+gamma2(4).*t_bar);
c_4 = -v0_dk.*gamma2(4).*(t_mat.^2);
lambda_tilde_4 = lambda+c_4.*theta;
beta_tilde_4 = beta_f-c_4.*theta./lambda_tilde_4.*beta_f-theta./lambda_tilde_4.*b_4;
I_4 = a_4-0.5.*log(lambda)./theta + 0.5.*log(lambda_tilde_4)./theta ...
    +0.5.*lambda.*beta_f.^2./theta -0.5.*lambda_tilde_4.*(beta_tilde_4).^2./theta;
R_4 = theta.*(I_4-(a_4+b_4.*beta_tilde_4+0.5.*c_4.*(beta_tilde_4).^2+0.5.*c_4./lambda_tilde_4));
pi_tilde_4 = weight4.*exp(-theta.*I_4); 
    
% % % % % % 
a_5 = -v0_dk.*(gamma0(5)+gamma1(5).*t_bar+0.5.*gamma2(5).*(t_bar.^2));
b_5 = -v0_dk.*t_mat.*(gamma1(5)+gamma2(5).*t_bar);
c_5 = -v0_dk.*gamma2(5).*(t_mat.^2);
lambda_tilde_5 = lambda+c_5.*theta;
beta_tilde_5 = beta_f-c_5.*theta./lambda_tilde_5.*beta_f-theta./lambda_tilde_5.*b_5;
I_5 = a_5-0.5.*log(lambda)./theta + 0.5.*log(lambda_tilde_5)./theta ...
    +0.5.*lambda.*beta_f.^2./theta -0.5.*lambda_tilde_5.*(beta_tilde_5).^2./theta;
R_5 = theta.*(I_5-(a_5+b_5.*beta_tilde_5+0.5.*c_5.*(beta_tilde_5).^2+0.5.*c_5./lambda_tilde_5));
pi_tilde_5 = weight5.*exp(-theta.*I_5); 
    
% % % % % % 
a_6 = -v0_dk.*(gamma0(6)+gamma1(6).*t_bar+0.5.*gamma2(6).*(t_bar.^2));
b_6 = -v0_dk.*t_mat.*(gamma1(6)+gamma2(6).*t_bar);
c_6 = -v0_dk.*gamma2(6).*(t_mat.^2);
lambda_tilde_6 = lambda+c_6.*theta;
beta_tilde_6 = beta_f-c_6.*theta./lambda_tilde_6.*beta_f-theta./lambda_tilde_6.*b_6;
I_6 = a_6-0.5.*log(lambda)./theta + 0.5.*log(lambda_tilde_6)./theta ...
    +0.5.*lambda.*beta_f.^2./theta -0.5.*lambda_tilde_6.*(beta_tilde_6).^2./theta;
R_6 = theta.*(I_6-(a_6+b_6.*beta_tilde_6+0.5.*c_6.*(beta_tilde_6).^2+0.5.*c_6./lambda_tilde_6));
pi_tilde_6 = weight6.*exp(-theta.*I_6); 
    
% % % % % % 
a_7 = -v0_dk.*(gamma0(7)+gamma1(7).*t_bar+0.5.*gamma2(7).*(t_bar.^2));
b_7 = -v0_dk.*t_mat.*(gamma1(7)+gamma2(7).*t_bar);
c_7 = -v0_dk.*gamma2(7).*(t_mat.^2);
lambda_tilde_7 = lambda+c_7.*theta;
beta_tilde_7 = beta_f-c_7.*theta./lambda_tilde_7.*beta_f-theta./lambda_tilde_7.*b_7;
I_7 = a_7-0.5.*log(lambda)./theta + 0.5.*log(lambda_tilde_7)./theta ...
    +0.5.*lambda.*beta_f.^2./theta -0.5.*lambda_tilde_7.*(beta_tilde_7).^2./theta;
R_7 = theta.*(I_7-(a_7+b_7.*beta_tilde_7+0.5.*c_7.*(beta_tilde_7).^2+0.5.*c_7./lambda_tilde_7));
pi_tilde_7 = weight7.*exp(-theta.*I_7); 
    
% % % % % % 
a_8 = -v0_dk.*(gamma0(8)+gamma1(8).*t_bar+0.5.*gamma2(8).*(t_bar.^2));
b_8 = -v0_dk.*t_mat.*(gamma1(8)+gamma2(8).*t_bar);
c_8 = -v0_dk.*gamma2(8).*(t_mat.^2);
lambda_tilde_8 = lambda+c_8.*theta;
beta_tilde_8 = beta_f-c_8.*theta./lambda_tilde_8.*beta_f-theta./lambda_tilde_8.*b_8;
I_8 = a_8-0.5.*log(lambda)./theta + 0.5.*log(lambda_tilde_8)./theta ...
    +0.5.*lambda.*beta_f.^2./theta -0.5.*lambda_tilde_8.*(beta_tilde_8).^2./theta;
R_8 = theta.*(I_8-(a_8+b_8.*beta_tilde_8+0.5.*c_8.*(beta_tilde_8).^2+0.5.*c_8./lambda_tilde_8));
pi_tilde_8 = weight8.*exp(-theta.*I_8); 
    
% % % % % % 
a_9 = -v0_dk.*(gamma0(9)+gamma1(9).*t_bar+0.5.*gamma2(9).*(t_bar.^2));
b_9 = -v0_dk.*t_mat.*(gamma1(9)+gamma2(9).*t_bar);
c_9 = -v0_dk.*gamma2(9).*(t_mat.^2);
lambda_tilde_9 = lambda+c_9.*theta;
beta_tilde_9 = beta_f-c_9.*theta./lambda_tilde_9.*beta_f-theta./lambda_tilde_9.*b_9;
I_9 = a_9-0.5.*log(lambda)./theta + 0.5.*log(lambda_tilde_9)./theta ...
    +0.5.*lambda.*beta_f.^2./theta -0.5.*lambda_tilde_9.*(beta_tilde_9).^2./theta;
R_9 = theta.*(I_9-(a_9+b_9.*beta_tilde_9+0.5.*c_9.*(beta_tilde_9).^2+0.5.*c_9./lambda_tilde_9));
pi_tilde_9 = weight9.*exp(-theta.*I_9); 
    
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

I_term = -1./theta.*log(pi_tilde_1+pi_tilde_2+pi_tilde_3+pi_tilde_4...
    +pi_tilde_5+pi_tilde_6+pi_tilde_7+pi_tilde_8+pi_tilde_9);

R_1 = theta.*(I_1 - J_1);
R_2 = theta.*(I_2 - J_2);
R_3 = theta.*(I_3 - J_3);
R_4 = theta.*(I_4 - J_4);
R_5 = theta.*(I_5 - J_5);
R_6 = theta.*(I_6 - J_6);
R_7 = theta.*(I_7 - J_7);
R_8 = theta.*(I_8 - J_8);
R_9 = theta.*(I_9 - J_9);

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
RE_total = 1./theta.*RE;

% inputs for solver
A = -delta.*ones(size(r_mat));
B_r = -e_star+Gamma_r.*(f.^Theta_r).*exp(Theta_r.*(k_mat-r_mat))-0.5.*(sigma_r.^2);
B_k = Alpha+Gamma.*log(1+i_k./Theta)-0.5.*(sigma_k.^2);
B_t = e_star.*exp(r_mat);
C_rr = 0.5.*sigma_r.^2.*ones(size(r_mat));
C_kk = 0.5.*sigma_k.^2.*ones(size(r_mat));
C_tt = zeros(size(r_mat));
D = delta.*alpha.*log(e_star)+delta.*alpha.*r_mat ...
    + delta.*(1-alpha).*(log(A_O-i_k-f)+(k_mat)) ...
    +I_term;

stateSpace = [r_mat(:), t_mat(:), k_mat(:)]; 
model      = {};
model.A    = A(:); 
model.B    = [B_r(:), B_t(:), B_k(:)];
model.C    = [C_rr(:), C_tt(:), C_kk(:)];
model.D    = D(:);
model.v0   = v0(:);
model.dt   = dt;

out = solveCGNatural(stateSpace, model);
out_comp = reshape(out,size(v0)).*ones(size(r_mat));
    
disp(['Error: ', num2str(max(max(max(abs(out_comp - v1_initial)))))])
v0 = v0.*ones(size(v0));
v0 = reshape(out,size(v0));

eta = 0.01;

while (max(max(max(max(abs(out_comp - vold)))))) > tol % check for convergence
    tic
    vold = v0 .* ones(size(v0));
    
    v0_dt = zeros(size(v0));
    v0_dt(:,2:end-1,:) = (1./(2.*ht)).*(v0(:,3:end,:)-v0(:,1:end-2,:));
    v0_dt(:,end,:) = (1./ht).*(v0(:,end,:)-v0(:,end-1,:));
    v0_dt(:,1,:) = (1./ht).*(v0(:,2,:)-v0(:,1,:));
    
    v0_dr = zeros(size(v0));
    v0_dr(2:end-1,:,:) = (1./(2.*hr)).*(v0(3:end,:,:)-v0(1:end-2,:,:));
    v0_dr(end,:,:) = (1./hr).*(v0(end,:,:)-v0(end-1,:,:));
    v0_dr(1,:,:) = (1./hr).*(v0(2,:,:)-v0(1,:,:));
    v0_dr(v0_dr<1e-8) = 1e-8;

    v0_dk = zeros(size(v0));
    v0_dk(:,:,2:end-1) = (1./(2.*hk)).*(v0(:,:,3:end)-v0(:,:,1:end-2));
    v0_dk(:,:,end) = (1./hk).*(v0(:,:,end)-v0(:,:,end-1));
    v0_dk(:,:,1) = (1./hk).*(v0(:,:,2)-v0(:,:,1));
       
    v0_drr = zeros(size(v0));
    v0_drr(2:end-1,:,:) = (1./(hr.^2)).*(v0(3:end,:,:)+v0(1:end-2,:,:)-2.*v0(2:end-1,:,:));
    v0_drr(end,:,:) = (1./(hr.^2)).*(v0(end,:,:)+v0(end-2,:,:)-2.*v0(end-1,:,:));
    v0_drr(1,:,:) = (1./(hr.^2)).*(v0(3,:,:)+v0(1,:,:)-2.*v0(2,:,:));
    
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
    C1 = delta.*alpha;
    e = C1./B1;
    e_star = e;
    
    q0 = delta.*(1-alpha)./(A_O-i_k-f);
    q = q0;

    Converged = 0;
    
    while Converged == 0
        istar = (Gamma./Theta.*v0_dk./q - 1).*Theta;
        jstar = (q.*exp(Theta_r.*(r_mat-k_mat))./((v0_dr).*Gamma_r.*Theta_r)).^(1./(Theta_r - 1));
        if A_O > (istar+jstar)
            qstar = eta.*delta.*(1-alpha)./(A_O-istar-jstar)+(1-eta).*q;
        else
            qstar = 2.*q;
        end

        if (max(max(max(max(abs(istar-i_k)))))<=1e-8) && (max(max(max(max(abs(jstar-f)))))<=1e-8)
            Converged = 1; 
        end
        q = qstar;
        i_k = istar.*(v0_dr>1e-8)+(v0_dr<=1e-8).*(v0_dk.*Gamma.*A_O - delta.*(1-alpha).*Theta)./(delta.*(1-alpha)+v0_dk.*Gamma);
        f = jstar.*(v0_dr>1e-8);
    end
  
    % % % % % %
    a_1 = -v0_dk.*(gamma0(1)+gamma1(1).*t_bar+0.5.*gamma2(1).*(t_bar.^2));
    b_1 = -v0_dk.*t_mat.*(gamma1(1)+gamma2(1).*t_bar);
    c_1 = -v0_dk.*gamma2(1).*(t_mat.^2);
    lambda_tilde_1 = lambda+c_1.*theta;
    beta_tilde_1 = beta_f-c_1.*theta./lambda_tilde_1.*beta_f-theta./lambda_tilde_1.*b_1;
    I_1 = a_1-0.5.*log(lambda)./theta + 0.5.*log(lambda_tilde_1)./theta ...
        +0.5.*lambda.*beta_f.^2./theta -0.5.*lambda_tilde_1.*(beta_tilde_1).^2./theta;
    R_1 = theta.*(I_1-(a_1+b_1.*beta_tilde_1+0.5.*c_1.*(beta_tilde_1).^2+0.5.*c_1./lambda_tilde_1));
    pi_tilde_1 = weight1.*exp(-theta.*I_1);    

    % % % % % %     
    a_2 = -v0_dk.*(gamma0(2)+gamma1(2).*t_bar+0.5.*gamma2(2).*(t_bar.^2));
    b_2 = -v0_dk.*t_mat.*(gamma1(2)+gamma2(2).*t_bar);
    c_2 = -v0_dk.*gamma2(2).*(t_mat.^2);
    lambda_tilde_2 = lambda+c_2.*theta;
    beta_tilde_2 = beta_f-c_2.*theta./lambda_tilde_2.*beta_f-theta./lambda_tilde_2.*b_2;
    I_2 = a_2-0.5.*log(lambda)./theta + 0.5.*log(lambda_tilde_2)./theta ...
        +0.5.*lambda.*beta_f.^2./theta -0.5.*lambda_tilde_2.*(beta_tilde_2).^2./theta;
    R_2 = theta.*(I_2-(a_2+b_2.*beta_tilde_2+0.5.*c_2.*(beta_tilde_2).^2+0.5.*c_2./lambda_tilde_2));
    pi_tilde_2 = weight2.*exp(-theta.*I_2);   

    % % % % % % 
    a_3 = -v0_dk.*(gamma0(3)+gamma1(3).*t_bar+0.5.*gamma2(3).*(t_bar.^2));
    b_3 = -v0_dk.*t_mat.*(gamma1(3)+gamma2(3).*t_bar);
    c_3 = -v0_dk.*gamma2(3).*(t_mat.^2);
    lambda_tilde_3 = lambda+c_3.*theta;
    beta_tilde_3 = beta_f-c_3.*theta./lambda_tilde_3.*beta_f-theta./lambda_tilde_3.*b_3;
    I_3 = a_3-0.5.*log(lambda)./theta + 0.5.*log(lambda_tilde_3)./theta ...
        +0.5.*lambda.*beta_f.^2./theta -0.5.*lambda_tilde_3.*(beta_tilde_3).^2./theta;
    R_3 = theta.*(I_3-(a_3+b_3.*beta_tilde_3+0.5.*c_3.*(beta_tilde_3).^2+0.5.*c_3./lambda_tilde_3));
    pi_tilde_3 = weight3.*exp(-theta.*I_3);   

    % % % % % % 
    a_4 = -v0_dk.*(gamma0(4)+gamma1(4).*t_bar+0.5.*gamma2(4).*(t_bar.^2));
    b_4 = -v0_dk.*t_mat.*(gamma1(4)+gamma2(4).*t_bar);
    c_4 = -v0_dk.*gamma2(4).*(t_mat.^2);
    lambda_tilde_4 = lambda+c_4.*theta;
    beta_tilde_4 = beta_f-c_4.*theta./lambda_tilde_4.*beta_f-theta./lambda_tilde_4.*b_4;
    I_4 = a_4-0.5.*log(lambda)./theta + 0.5.*log(lambda_tilde_4)./theta ...
        +0.5.*lambda.*beta_f.^2./theta -0.5.*lambda_tilde_4.*(beta_tilde_4).^2./theta;
    R_4 = theta.*(I_4-(a_4+b_4.*beta_tilde_4+0.5.*c_4.*(beta_tilde_4).^2+0.5.*c_4./lambda_tilde_4));
    pi_tilde_4 = weight4.*exp(-theta.*I_4); 

    % % % % % % 
    a_5 = -v0_dk.*(gamma0(5)+gamma1(5).*t_bar+0.5.*gamma2(5).*(t_bar.^2));
    b_5 = -v0_dk.*t_mat.*(gamma1(5)+gamma2(5).*t_bar);
    c_5 = -v0_dk.*gamma2(5).*(t_mat.^2);
    lambda_tilde_5 = lambda+c_5.*theta;
    beta_tilde_5 = beta_f-c_5.*theta./lambda_tilde_5.*beta_f-theta./lambda_tilde_5.*b_5;
    I_5 = a_5-0.5.*log(lambda)./theta + 0.5.*log(lambda_tilde_5)./theta ...
        +0.5.*lambda.*beta_f.^2./theta -0.5.*lambda_tilde_5.*(beta_tilde_5).^2./theta;
    R_5 = theta.*(I_5-(a_5+b_5.*beta_tilde_5+0.5.*c_5.*(beta_tilde_5).^2+0.5.*c_5./lambda_tilde_5));
    pi_tilde_5 = weight5.*exp(-theta.*I_5); 

    % % % % % % 
    a_6 = -v0_dk.*(gamma0(6)+gamma1(6).*t_bar+0.5.*gamma2(6).*(t_bar.^2));
    b_6 = -v0_dk.*t_mat.*(gamma1(6)+gamma2(6).*t_bar);
    c_6 = -v0_dk.*gamma2(6).*(t_mat.^2);
    lambda_tilde_6 = lambda+c_6.*theta;
    beta_tilde_6 = beta_f-c_6.*theta./lambda_tilde_6.*beta_f-theta./lambda_tilde_6.*b_6;
    I_6 = a_6-0.5.*log(lambda)./theta + 0.5.*log(lambda_tilde_6)./theta ...
        +0.5.*lambda.*beta_f.^2./theta -0.5.*lambda_tilde_6.*(beta_tilde_6).^2./theta;
    R_6 = theta.*(I_6-(a_6+b_6.*beta_tilde_6+0.5.*c_6.*(beta_tilde_6).^2+0.5.*c_6./lambda_tilde_6));
    pi_tilde_6 = weight6.*exp(-theta.*I_6); 

    % % % % % % 
    a_7 = -v0_dk.*(gamma0(7)+gamma1(7).*t_bar+0.5.*gamma2(7).*(t_bar.^2));
    b_7 = -v0_dk.*t_mat.*(gamma1(7)+gamma2(7).*t_bar);
    c_7 = -v0_dk.*gamma2(7).*(t_mat.^2);
    lambda_tilde_7 = lambda+c_7.*theta;
    beta_tilde_7 = beta_f-c_7.*theta./lambda_tilde_7.*beta_f-theta./lambda_tilde_7.*b_7;
    I_7 = a_7-0.5.*log(lambda)./theta + 0.5.*log(lambda_tilde_7)./theta ...
        +0.5.*lambda.*beta_f.^2./theta -0.5.*lambda_tilde_7.*(beta_tilde_7).^2./theta;
    R_7 = theta.*(I_7-(a_7+b_7.*beta_tilde_7+0.5.*c_7.*(beta_tilde_7).^2+0.5.*c_7./lambda_tilde_7));
    pi_tilde_7 = weight7.*exp(-theta.*I_7); 

    % % % % % % 
    a_8 = -v0_dk.*(gamma0(8)+gamma1(8).*t_bar+0.5.*gamma2(8).*(t_bar.^2));
    b_8 = -v0_dk.*t_mat.*(gamma1(8)+gamma2(8).*t_bar);
    c_8 = -v0_dk.*gamma2(8).*(t_mat.^2);
    lambda_tilde_8 = lambda+c_8.*theta;
    beta_tilde_8 = beta_f-c_8.*theta./lambda_tilde_8.*beta_f-theta./lambda_tilde_8.*b_8;
    I_8 = a_8-0.5.*log(lambda)./theta + 0.5.*log(lambda_tilde_8)./theta ...
        +0.5.*lambda.*beta_f.^2./theta -0.5.*lambda_tilde_8.*(beta_tilde_8).^2./theta;
    R_8 = theta.*(I_8-(a_8+b_8.*beta_tilde_8+0.5.*c_8.*(beta_tilde_8).^2+0.5.*c_8./lambda_tilde_8));
    pi_tilde_8 = weight8.*exp(-theta.*I_8); 

    % % % % % % 
    a_9 = -v0_dk.*(gamma0(9)+gamma1(9).*t_bar+0.5.*gamma2(9).*(t_bar.^2));
    b_9 = -v0_dk.*t_mat.*(gamma1(9)+gamma2(9).*t_bar);
    c_9 = -v0_dk.*gamma2(9).*(t_mat.^2);
    lambda_tilde_9 = lambda+c_9.*theta;
    beta_tilde_9 = beta_f-c_9.*theta./lambda_tilde_9.*beta_f-theta./lambda_tilde_9.*b_9;
    I_9 = a_9-0.5.*log(lambda)./theta + 0.5.*log(lambda_tilde_9)./theta ...
        +0.5.*lambda.*beta_f.^2./theta -0.5.*lambda_tilde_9.*(beta_tilde_9).^2./theta;
    R_9 = theta.*(I_9-(a_9+b_9.*beta_tilde_9+0.5.*c_9.*(beta_tilde_9).^2+0.5.*c_9./lambda_tilde_9));
    pi_tilde_9 = weight9.*exp(-theta.*I_9); 

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

    I_term = -1./theta.*log(pi_tilde_1+pi_tilde_2+pi_tilde_3+pi_tilde_4...
        +pi_tilde_5+pi_tilde_6+pi_tilde_7+pi_tilde_8+pi_tilde_9);

    R_1 = theta.*(I_1 - J_1);
    R_2 = theta.*(I_2 - J_2);
    R_3 = theta.*(I_3 - J_3);
    R_4 = theta.*(I_4 - J_4);
    R_5 = theta.*(I_5 - J_5);
    R_6 = theta.*(I_6 - J_6);
    R_7 = theta.*(I_7 - J_7);
    R_8 = theta.*(I_8 - J_8);
    R_9 = theta.*(I_9 - J_9);

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
    RE_total = 1./theta.*RE;

    % inputs for solver
    A = -delta.*ones(size(r_mat));
    B_r = -e_star+Gamma_r.*(f.^Theta_r).*exp(Theta_r.*(k_mat-r_mat))-0.5.*(sigma_r.^2);
    B_k = Alpha+Gamma.*log(1+i_k./Theta)-0.5.*(sigma_k.^2);
    B_t = e_star.*exp(r_mat);
    C_rr = 0.5.*sigma_r.^2.*ones(size(r_mat));
    C_kk = 0.5.*sigma_k.^2.*ones(size(r_mat));
    C_tt = zeros(size(r_mat));
    D = delta.*alpha.*log(e_star)+delta.*alpha.*r_mat ...
        + delta.*(1-alpha).*(log(A_O-i_k-f)+(k_mat)) ...
        +I_term;

    stateSpace = [r_mat(:), t_mat(:), k_mat(:)]; 
    model      = {};
    model.A    = A(:); 
    model.B    = [B_r(:), B_t(:), B_k(:)];
    model.C    = [C_rr(:), C_tt(:), C_kk(:)];
    model.D    = D(:);
    model.v0   = v0(:);
    model.dt   = dt;

    out = solveCGNatural(stateSpace, model);
    out_comp = reshape(out,size(v0)).*ones(size(r_mat));

    iter = iter+1;
    disp(['Diff: ', num2str(max(max(max(abs(out_comp - vold))))),' Number of Iters: ',num2str(iter)])

    v0 = v0.*ones(size(v0));
    v0 = reshape(out,size(v0));  

    pde_error = A.*v0+B_r.*v0_dr+B_t.*v0_dt+B_k.*v0_dk+C_rr.*v0_drr+C_kk.*v0_dkk+C_tt.*v0_dtt+D;
    disp(['PDE Error: ', num2str(max(max(max(abs(pde_error)))))])

    if (max(max(abs(out_comp - vold)))) < tol
        fprintf('PDE Converges');
    end

    toc
end

s0 = '/HJB_NonLinGrowth_NoAmb';
filename = [pwd, s0];
save(filename); 

%% Step 3: simulation

filename2 = [filename,'_Sims'];

T = 100; 
pers = 4*T; 
dt = T/pers;
nDims = 4;
its = 1;

R_0 = 650;
K_0 = 80./A_O;
T_0 = (870-580);

efunc = griddedInterpolant(r_mat,t_mat,k_mat,e,'spline');
ffunc = griddedInterpolant(r_mat,t_mat,k_mat,f,'linear');
i_kfunc = griddedInterpolant(r_mat,t_mat,k_mat,i_k,'spline');

v_drfunc = griddedInterpolant(r_mat,t_mat,k_mat,v0_dr,'linear');
v_dtfunc = griddedInterpolant(r_mat,t_mat,k_mat,v0_dt,'spline');
v_dkfunc = griddedInterpolant(r_mat,t_mat,k_mat,v0_dk,'spline');
v_func = griddedInterpolant(r_mat,t_mat,k_mat,out_comp,'spline');

pi_tilde_1func = griddedInterpolant(r_mat,t_mat,k_mat,pi_tilde_1_norm,'spline');
pi_tilde_2func = griddedInterpolant(r_mat,t_mat,k_mat,pi_tilde_2_norm,'spline');
pi_tilde_3func = griddedInterpolant(r_mat,t_mat,k_mat,pi_tilde_3_norm,'spline');
pi_tilde_4func = griddedInterpolant(r_mat,t_mat,k_mat,pi_tilde_4_norm,'spline');
pi_tilde_5func = griddedInterpolant(r_mat,t_mat,k_mat,pi_tilde_5_norm,'spline');
pi_tilde_6func = griddedInterpolant(r_mat,t_mat,k_mat,pi_tilde_6_norm,'spline');
pi_tilde_7func = griddedInterpolant(r_mat,t_mat,k_mat,pi_tilde_7_norm,'spline');
pi_tilde_8func = griddedInterpolant(r_mat,t_mat,k_mat,pi_tilde_8_norm,'spline');
pi_tilde_9func = griddedInterpolant(r_mat,t_mat,k_mat,pi_tilde_9_norm,'spline');

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

REfunc = griddedInterpolant(r_mat,t_mat,k_mat,RE,'spline');
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

base_driftKfunc1 = griddedInterpolant(r_mat,t_mat,k_mat,dmg_1,'spline');
base_driftKfunc2 = griddedInterpolant(r_mat,t_mat,k_mat,dmg_2,'spline');
base_driftKfunc3 = griddedInterpolant(r_mat,t_mat,k_mat,dmg_3,'spline');
base_driftKfunc4 = griddedInterpolant(r_mat,t_mat,k_mat,dmg_4,'spline');
base_driftKfunc5 = griddedInterpolant(r_mat,t_mat,k_mat,dmg_5,'spline');
base_driftKfunc6 = griddedInterpolant(r_mat,t_mat,k_mat,dmg_6,'spline');
base_driftKfunc7 = griddedInterpolant(r_mat,t_mat,k_mat,dmg_7,'spline');
base_driftKfunc8 = griddedInterpolant(r_mat,t_mat,k_mat,dmg_8,'spline');
base_driftKfunc9 = griddedInterpolant(r_mat,t_mat,k_mat,dmg_9,'spline');
base_driftK_func1 = @(x) base_driftKfunc1(log(x(:,1)),x(:,3),log(x(:,2)));
base_driftK_func2 = @(x) base_driftKfunc2(log(x(:,1)),x(:,3),log(x(:,2)));
base_driftK_func3 = @(x) base_driftKfunc3(log(x(:,1)),x(:,3),log(x(:,2)));
base_driftK_func4 = @(x) base_driftKfunc4(log(x(:,1)),x(:,3),log(x(:,2)));
base_driftK_func5 = @(x) base_driftKfunc5(log(x(:,1)),x(:,3),log(x(:,2)));
base_driftK_func6 = @(x) base_driftKfunc6(log(x(:,1)),x(:,3),log(x(:,2)));
base_driftK_func7 = @(x) base_driftKfunc7(log(x(:,1)),x(:,3),log(x(:,2)));
base_driftK_func8 = @(x) base_driftKfunc8(log(x(:,1)),x(:,3),log(x(:,2)));
base_driftK_func9 = @(x) base_driftKfunc9(log(x(:,1)),x(:,3),log(x(:,2)));

tilt_driftKfunc1 = griddedInterpolant(r_mat,t_mat,k_mat,dmg_tilt_1,'spline');
tilt_driftKfunc2 = griddedInterpolant(r_mat,t_mat,k_mat,dmg_tilt_2,'spline');
tilt_driftKfunc3 = griddedInterpolant(r_mat,t_mat,k_mat,dmg_tilt_3,'spline');
tilt_driftKfunc4 = griddedInterpolant(r_mat,t_mat,k_mat,dmg_tilt_4,'spline');
tilt_driftKfunc5 = griddedInterpolant(r_mat,t_mat,k_mat,dmg_tilt_5,'spline');
tilt_driftKfunc6 = griddedInterpolant(r_mat,t_mat,k_mat,dmg_tilt_6,'spline');
tilt_driftKfunc7 = griddedInterpolant(r_mat,t_mat,k_mat,dmg_tilt_7,'spline');
tilt_driftKfunc8 = griddedInterpolant(r_mat,t_mat,k_mat,dmg_tilt_8,'spline');
tilt_driftKfunc9 = griddedInterpolant(r_mat,t_mat,k_mat,dmg_tilt_9,'spline');
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
muK_tilted = @(x) (Alpha + Gamma.*log(1+i_k_func(x)./Theta)-Gamma_tilted(x));
muT = @(x) e_func(x).*x(:,1);
muK_base = @(x) (Alpha + Gamma.*log(1+i_k_func(x)./Theta)-Gamma_base(x));

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
hist2(j,2) = max(min(hist2(j-1,2).*exp((muK_tilted(hist2(j-1,:))-0.5.*sum((sigmaK(hist2(j-1,:))).^2) )* dt ...
                                  +sigmaK(hist2(j-1,:))* shock'), upperBounds(:,2)), lowerBounds(:,2));                              
hist2(j,3) = max(min(hist2(j-1,3) + muT(hist2(j-1,:)) * dt + sigmaT(hist2(j-1,:))* shock', upperBounds(:,3)), lowerBounds(:,3)); 
hist2(j,4) = max(min(hist2(j-1,4).*exp((muK_base(hist2(j-1,:))-0.5.*sum((sigmaK(hist2(j-1,:))).^2) )* dt ...
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


% save results 
save(filename2);
e_values = mean(e_hists2,2);
s1 = num2str(1./theta,4);
s1 = strrep(s1,'.','');
s_e = strcat('e_A=',s1,'.mat');
save(s_e,'e_values');

%% Step 4: SCC calculation
close all
clear all
clc

load([pwd, '/HJB_NonLinGrowth_NoAmb']);

a1 = v0_dk.*zeros(size(t_mat));
b1 = v0_dk.*(gamma1(1)+gamma2(1).*t_bar);
c1 = 2.*v0_dk.* gamma2(1).*(t_mat);

a2 = v0_dk.*zeros(size(t_mat));
b2 = v0_dk.*(gamma1(2)+gamma2(2).*t_bar);
c2 = 2.*v0_dk.*gamma2(2).*(t_mat);

a3 = v0_dk.*zeros(size(t_mat));
b3 = v0_dk.*(gamma1(3)+gamma2(3).*t_bar);
c3 = 2.*v0_dk.*gamma2(3).*(t_mat);

a4 = v0_dk.*zeros(size(t_mat));
b4 = v0_dk.*(gamma1(4)+gamma2(4).*t_bar);
c4 = 2.*v0_dk.*gamma2(4).*(t_mat);

a5 = v0_dk.*zeros(size(t_mat));
b5 = v0_dk.*(gamma1(5)+gamma2(5).*t_bar);
c5 = 2.*v0_dk.*gamma2(5).*(t_mat);

a6 = v0_dk.*zeros(size(t_mat));
b6 = v0_dk.*(gamma1(6)+gamma2(6).*t_bar);
c6 = 2.*v0_dk.*gamma2(6).*(t_mat);

a7 = v0_dk.*zeros(size(t_mat));
b7 = v0_dk.*(gamma1(7)+gamma2(7).*t_bar);
c7 = 2.*v0_dk.*gamma2(7).*(t_mat);

a8 = v0_dk.*zeros(size(t_mat));
b8 = v0_dk.*(gamma1(8)+gamma2(8).*t_bar);
c8 = 2.*v0_dk.*gamma2(8).*(t_mat);

a9 = v0_dk.*zeros(size(t_mat));
b9 = v0_dk.*(gamma1(9)+gamma2(9).*t_bar);
c9 = 2.*v0_dk.*gamma2(9).*(t_mat); 

dmg_1 = a1+b1.*beta_f+0.5.*c1.*(beta_f.^2)+0.5.*c1./lambda;    
dmg_2 = a2+b2.*beta_f+0.5.*c2.*(beta_f.^2)+0.5.*c2./lambda;    
dmg_3 = a3+b3.*beta_f+0.5.*c3.*(beta_f.^2)+0.5.*c3./lambda;    
dmg_4 = a4+b4.*beta_f+0.5.*c4.*(beta_f.^2)+0.5.*c4./lambda;    
dmg_5 = a5+b5.*beta_f+0.5.*c5.*(beta_f.^2)+0.5.*c5./lambda;    
dmg_6 = a6+b6.*beta_f+0.5.*c6.*(beta_f.^2)+0.5.*c6./lambda;    
dmg_7 = a7+b7.*beta_f+0.5.*c7.*(beta_f.^2)+0.5.*c7./lambda;    
dmg_8 = a8+b8.*beta_f+0.5.*c8.*(beta_f.^2)+0.5.*c8./lambda;    
dmg_9 = a9+b9.*beta_f+0.5.*c9.*(beta_f.^2)+0.5.*c9./lambda;    

flow_base =  weight1.*dmg_1+ weight2.*dmg_2...
      + weight3.*dmg_3+ weight4.*dmg_4...
      + weight5.*dmg_5+ weight6.*dmg_6...
      + weight7.*dmg_7+ weight8.*dmg_8...
      + weight9.*dmg_9;

a1 = v0_dk.*(gamma0(1)+gamma1(1).*t_bar+0.5.*gamma2(1).*(t_bar.^2));
b1 = v0_dk.*t_mat.*(gamma1(1)+gamma2(1).*t_bar);
c1 = v0_dk.*gamma2(1).*(t_mat.^2);

a2 = v0_dk.*(gamma0(2)+gamma1(2).*t_bar+0.5.*gamma2(2).*(t_bar.^2));
b2 = v0_dk.*t_mat.*(gamma1(2)+gamma2(2).*t_bar);
c2 = v0_dk.*gamma2(2).*(t_mat.^2);

a3 = v0_dk.*(gamma0(3)+gamma1(3).*t_bar+0.5.*gamma2(3).*(t_bar.^2));
b3 = v0_dk.*t_mat.*(gamma1(3)+gamma2(3).*t_bar);
c3 = v0_dk.*gamma2(3).*(t_mat.^2);

a4 = v0_dk.*(gamma0(4)+gamma1(4).*t_bar+0.5.*gamma2(4).*(t_bar.^2));
b4 = v0_dk.*t_mat.*(gamma1(4)+gamma2(4).*t_bar);
c4 = v0_dk.*gamma2(4).*(t_mat.^2);

a5 = v0_dk.*(gamma0(5)+gamma1(5).*t_bar+0.5.*gamma2(5).*(t_bar.^2));
b5 = v0_dk.*t_mat.*(gamma1(5)+gamma2(5).*t_bar);
c5 = v0_dk.*gamma2(5).*(t_mat.^2);

a6 = v0_dk.*(gamma0(6)+gamma1(6).*t_bar+0.5.*gamma2(6).*(t_bar.^2));
b6 = v0_dk.*t_mat.*(gamma1(6)+gamma2(6).*t_bar);
c6 = v0_dk.*gamma2(6).*(t_mat.^2);

a7 = v0_dk.*(gamma0(7)+gamma1(7).*t_bar+0.5.*gamma2(7).*(t_bar.^2));
b7 = v0_dk.*t_mat.*(gamma1(7)+gamma2(7).*t_bar);
c7 = v0_dk.*gamma2(7).*(t_mat.^2);

a8 = v0_dk.*(gamma0(8)+gamma1(8).*t_bar+0.5.*gamma2(8).*(t_bar.^2));
b8 = v0_dk.*t_mat.*(gamma1(8)+gamma2(8).*t_bar);
c8 = v0_dk.*gamma2(8).*(t_mat.^2);

a9 = v0_dk.*(gamma0(9)+gamma1(9).*t_bar+0.5.*gamma2(9).*(t_bar.^2));
b9 = v0_dk.*t_mat.*(gamma1(9)+gamma2(9).*t_bar);
c9 = v0_dk.*gamma2(9).*(t_mat.^2);  

dmg_1 = a1+b1.*beta_f+0.5.*c1.*(beta_f.^2)+0.5.*c1./lambda;    
dmg_2 = a2+b2.*beta_f+0.5.*c2.*(beta_f.^2)+0.5.*c2./lambda;    
dmg_3 = a3+b3.*beta_f+0.5.*c3.*(beta_f.^2)+0.5.*c3./lambda;    
dmg_4 = a4+b4.*beta_f+0.5.*c4.*(beta_f.^2)+0.5.*c4./lambda;    
dmg_5 = a5+b5.*beta_f+0.5.*c5.*(beta_f.^2)+0.5.*c5./lambda;    
dmg_6 = a6+b6.*beta_f+0.5.*c6.*(beta_f.^2)+0.5.*c6./lambda;    
dmg_7 = a7+b7.*beta_f+0.5.*c7.*(beta_f.^2)+0.5.*c7./lambda;    
dmg_8 = a8+b8.*beta_f+0.5.*c8.*(beta_f.^2)+0.5.*c8./lambda;    
dmg_9 = a9+b9.*beta_f+0.5.*c9.*(beta_f.^2)+0.5.*c9./lambda;    

Gamma_base =  weight1.*dmg_1+ weight2.*dmg_2...
      + weight3.*dmg_3+ weight4.*dmg_4...
      + weight5.*dmg_5+ weight6.*dmg_6...
      + weight7.*dmg_7+ weight8.*dmg_8...
      + weight9.*dmg_9;

A = -delta.*ones(size(r_mat));
B_r = -e+Gamma_r.*(f.^Theta_r).*exp(Theta_r.*(k_mat-r_mat))-0.5.*(sigma_r.^2);
B_k = Alpha+Gamma.*log(1+i_k./Theta)-0.5.*(sigma_k.^2)-Gamma_base;
B_t = e_star.*exp(r_mat);
C_rr = 0.5.*sigma_r.^2.*ones(size(r_mat));
C_kk = 0.5.*sigma_k.^2.*ones(size(r_mat));
C_tt = zeros(size(r_mat));

D = flow_base;

stateSpace = [r_mat(:), t_mat(:), k_mat(:)]; 
model      = {};
model.A    = A(:); 
model.B    = [B_r(:), B_t(:), B_k(:)];
model.C    = [C_rr(:), C_tt(:), C_kk(:)];
model.D    = D(:);
model.v0   = v0(:).*0;
model.dt   = 1.0;

out = solveCGNatural_1(stateSpace, model);
v0 = reshape(out,size(r_mat));
filename2 = [pwd, '/SCC_mat_Cumu_base_GrowthNoAmb'];

save(filename2)
%%%%%%%%%
close all
clear all
clc

load([pwd, '/HJB_NonLinGrowth_NoAmb']);

a1 = v0_dk.*zeros(size(t_mat));
b1 = v0_dk.*(gamma1(1)+gamma2(1).*t_bar);
c1 = 2.*v0_dk.* gamma2(1).*(t_mat);

a2 = v0_dk.*zeros(size(t_mat));
b2 = v0_dk.*(gamma1(2)+gamma2(2).*t_bar);
c2 = 2.*v0_dk.*gamma2(2).*(t_mat);

a3 = v0_dk.*zeros(size(t_mat));
b3 = v0_dk.*(gamma1(3)+gamma2(3).*t_bar);
c3 = 2.*v0_dk.*gamma2(3).*(t_mat);

a4 = v0_dk.*zeros(size(t_mat));
b4 = v0_dk.*(gamma1(4)+gamma2(4).*t_bar);
c4 = 2.*v0_dk.*gamma2(4).*(t_mat);

a5 = v0_dk.*zeros(size(t_mat));
b5 = v0_dk.*(gamma1(5)+gamma2(5).*t_bar);
c5 = 2.*v0_dk.*gamma2(5).*(t_mat);

a6 = v0_dk.*zeros(size(t_mat));
b6 = v0_dk.*(gamma1(6)+gamma2(6).*t_bar);
c6 = 2.*v0_dk.*gamma2(6).*(t_mat);

a7 = v0_dk.*zeros(size(t_mat));
b7 = v0_dk.*(gamma1(7)+gamma2(7).*t_bar);
c7 = 2.*v0_dk.*gamma2(7).*(t_mat);

a8 = v0_dk.*zeros(size(t_mat));
b8 = v0_dk.*(gamma1(8)+gamma2(8).*t_bar);
c8 = 2.*v0_dk.*gamma2(8).*(t_mat);

a9 = v0_dk.*zeros(size(t_mat));
b9 = v0_dk.*(gamma1(9)+gamma2(9).*t_bar);
c9 = 2.*v0_dk.*gamma2(9).*(t_mat);

dmg_tilt_1 = a1+b1.*beta_tilde_1+0.5.*c1.*(beta_tilde_1.^2)+0.5.*c1./lambda_tilde_1;    
dmg_tilt_2 = a2+b2.*beta_tilde_2+0.5.*c2.*(beta_tilde_2.^2)+0.5.*c2./lambda_tilde_2;    
dmg_tilt_3 = a3+b3.*beta_tilde_3+0.5.*c3.*(beta_tilde_3.^2)+0.5.*c3./lambda_tilde_3;    
dmg_tilt_4 = a4+b4.*beta_tilde_4+0.5.*c4.*(beta_tilde_4.^2)+0.5.*c4./lambda_tilde_4;    
dmg_tilt_5 = a5+b5.*beta_tilde_5+0.5.*c5.*(beta_tilde_5.^2)+0.5.*c5./lambda_tilde_5;    
dmg_tilt_6 = a6+b6.*beta_tilde_6+0.5.*c6.*(beta_tilde_6.^2)+0.5.*c6./lambda_tilde_6;    
dmg_tilt_7 = a7+b7.*beta_tilde_7+0.5.*c7.*(beta_tilde_7.^2)+0.5.*c7./lambda_tilde_7;    
dmg_tilt_8 = a8+b8.*beta_tilde_8+0.5.*c8.*(beta_tilde_8.^2)+0.5.*c8./lambda_tilde_8;    
dmg_tilt_9 = a9+b9.*beta_tilde_9+0.5.*c9.*(beta_tilde_9.^2)+0.5.*c9./lambda_tilde_9;    

flow_tilted =    pi_tilde_1_norm.*dmg_tilt_1+pi_tilde_2_norm.*dmg_tilt_2...
     +pi_tilde_3_norm.*dmg_tilt_3+pi_tilde_4_norm.*dmg_tilt_4...
     +pi_tilde_5_norm.*dmg_tilt_5+pi_tilde_6_norm.*dmg_tilt_6...
     +pi_tilde_7_norm.*dmg_tilt_7+pi_tilde_8_norm.*dmg_tilt_8...
     +pi_tilde_9_norm.*dmg_tilt_9;
 
% % % % % % % % % % % % % % % % % % % % % % % 
a1 = v0_dk.*(gamma0(1)+gamma1(1).*t_bar+0.5.*gamma2(1).*(t_bar.^2));
b1 = v0_dk.*t_mat.*(gamma1(1)+gamma2(1).*t_bar);
c1 = v0_dk.*gamma2(1).*(t_mat.^2);

a2 = v0_dk.*(gamma0(2)+gamma1(2).*t_bar+0.5.*gamma2(2).*(t_bar.^2));
b2 = v0_dk.*t_mat.*(gamma1(2)+gamma2(2).*t_bar);
c2 = v0_dk.*gamma2(2).*(t_mat.^2);

a3 = v0_dk.*(gamma0(3)+gamma1(3).*t_bar+0.5.*gamma2(3).*(t_bar.^2));
b3 = v0_dk.*t_mat.*(gamma1(3)+gamma2(3).*t_bar);
c3 = v0_dk.*gamma2(3).*(t_mat.^2);

a4 = v0_dk.*(gamma0(4)+gamma1(4).*t_bar+0.5.*gamma2(4).*(t_bar.^2));
b4 = v0_dk.*t_mat.*(gamma1(4)+gamma2(4).*t_bar);
c4 = v0_dk.*gamma2(4).*(t_mat.^2);

a5 = v0_dk.*(gamma0(5)+gamma1(5).*t_bar+0.5.*gamma2(5).*(t_bar.^2));
b5 = v0_dk.*t_mat.*(gamma1(5)+gamma2(5).*t_bar);
c5 = v0_dk.*gamma2(5).*(t_mat.^2);

a6 = v0_dk.*(gamma0(6)+gamma1(6).*t_bar+0.5.*gamma2(6).*(t_bar.^2));
b6 = v0_dk.*t_mat.*(gamma1(6)+gamma2(6).*t_bar);
c6 = v0_dk.*gamma2(6).*(t_mat.^2);

a7 = v0_dk.*(gamma0(7)+gamma1(7).*t_bar+0.5.*gamma2(7).*(t_bar.^2));
b7 = v0_dk.*t_mat.*(gamma1(7)+gamma2(7).*t_bar);
c7 = v0_dk.*gamma2(7).*(t_mat.^2);

a8 = v0_dk.*(gamma0(8)+gamma1(8).*t_bar+0.5.*gamma2(8).*(t_bar.^2));
b8 = v0_dk.*t_mat.*(gamma1(8)+gamma2(8).*t_bar);
c8 = v0_dk.*gamma2(8).*(t_mat.^2);

a9 = v0_dk.*(gamma0(9)+gamma1(9).*t_bar+0.5.*gamma2(9).*(t_bar.^2));
b9 = v0_dk.*t_mat.*(gamma1(9)+gamma2(9).*t_bar);
c9 = v0_dk.*gamma2(9).*(t_mat.^2);

dmg_tilt_1 = a1+b1.*beta_tilde_1+0.5.*c1.*(beta_tilde_1.^2)+0.5.*c1./lambda_tilde_1;    
dmg_tilt_2 = a2+b2.*beta_tilde_2+0.5.*c2.*(beta_tilde_2.^2)+0.5.*c2./lambda_tilde_2;    
dmg_tilt_3 = a3+b3.*beta_tilde_3+0.5.*c3.*(beta_tilde_3.^2)+0.5.*c3./lambda_tilde_3;    
dmg_tilt_4 = a4+b4.*beta_tilde_4+0.5.*c4.*(beta_tilde_4.^2)+0.5.*c4./lambda_tilde_4;    
dmg_tilt_5 = a5+b5.*beta_tilde_5+0.5.*c5.*(beta_tilde_5.^2)+0.5.*c5./lambda_tilde_5;    
dmg_tilt_6 = a6+b6.*beta_tilde_6+0.5.*c6.*(beta_tilde_6.^2)+0.5.*c6./lambda_tilde_6;    
dmg_tilt_7 = a7+b7.*beta_tilde_7+0.5.*c7.*(beta_tilde_7.^2)+0.5.*c7./lambda_tilde_7;    
dmg_tilt_8 = a8+b8.*beta_tilde_8+0.5.*c8.*(beta_tilde_8.^2)+0.5.*c8./lambda_tilde_8;    
dmg_tilt_9 = a9+b9.*beta_tilde_9+0.5.*c9.*(beta_tilde_9.^2)+0.5.*c9./lambda_tilde_9;    


Gamma_tilted =    pi_tilde_1_norm.*dmg_tilt_1+pi_tilde_2_norm.*dmg_tilt_2...
     +pi_tilde_3_norm.*dmg_tilt_3+pi_tilde_4_norm.*dmg_tilt_4...
     +pi_tilde_5_norm.*dmg_tilt_5+pi_tilde_6_norm.*dmg_tilt_6...
     +pi_tilde_7_norm.*dmg_tilt_7+pi_tilde_8_norm.*dmg_tilt_8...
     +pi_tilde_9_norm.*dmg_tilt_9;

A = -delta.*ones(size(r_mat));
B_r = -e+Gamma_r.*(f.^Theta_r).*exp(Theta_r.*(k_mat-r_mat))-0.5.*(sigma_r.^2);
B_k = Alpha+Gamma.*log(1+i_k./Theta)-0.5.*(sigma_k.^2)-Gamma_tilted;
B_t = e_star.*exp(r_mat);
C_rr = 0.5.*sigma_r.^2.*ones(size(r_mat));
C_kk = 0.5.*sigma_k.^2.*ones(size(r_mat));
C_tt = zeros(size(r_mat));

D = flow_tilted;

stateSpace = [r_mat(:), t_mat(:), k_mat(:)]; 
model      = {};
model.A    = A(:); 
model.B    = [B_r(:), B_t(:), B_k(:)];
model.C    = [C_rr(:), C_tt(:), C_kk(:)];
model.D    = D(:);
model.v0   = v0(:).*0;
model.dt   = 1.0;

out = solveCGNatural_1(stateSpace, model);
v0 = reshape(out,size(r_mat));
filename2 = [pwd, '/SCC_mat_Cumu_worst_GrowthNoAmb'];
save(filename2)

% SCC plots
close all;
clc;
clear all;

file2 = [pwd, '/SCC_mat_Cumu_base_GrowthNoAmb'];

Model2 = load(file2,'v0','r_mat','k_mat','t_mat','i_k','f','e','v0_dk','v0_dr','theta',...
    'A_O','alpha','delta','Theta','Gamma','gamma_1','gamma_2','t_bar',...
    'beta_f','var_beta_f');

external_v0 = Model2.v0;
r_mat_1 = Model2.r_mat;
k_mat_1 = Model2.k_mat;
t_mat_1 = Model2.t_mat;
gamma_1 = Model2.gamma_1;
gamma_2 = Model2.gamma_2;
t_bar = Model2.t_bar;
theta = Model2.theta;
A_O = Model2.A_O;
alpha = Model2.alpha;
delta = Model2.delta;
Gamma = Model2.Gamma;
Theta = Model2.Theta;
beta_f = Model2.beta_f;
var_beta_f = Model2.var_beta_f;
e_1 = Model2.e;
f_1 = Model2.f;
i_k_1 = Model2.i_k;

time_vec = linspace(0,100,400);

file2 = [pwd, '/HJB_NonLinGrowth_NoAmb'];

Model2 = load(file2,'v0_dr','v0_dt');
v0_dr_1 = Model2.v0_dr;
v0_dt_1 = Model2.v0_dt;

MC = delta.*(1-alpha)./(A_O.*exp(k_mat_1)-i_k_1.*exp(k_mat_1)-f_1.*exp(k_mat_1));
ME = delta.*alpha./(e_1.*exp(r_mat_1));
SCC = 1000*ME./MC;
SCC_func = griddedInterpolant(r_mat_1,t_mat_1,k_mat_1,SCC,'spline');

ME1 = (v0_dr_1.*exp(-r_mat_1));
SCC1 = 1000*ME1./MC; 
SCC1_func = griddedInterpolant(r_mat_1,t_mat_1,k_mat_1,SCC1,'linear');

ME2_base =  external_v0;
SCC2_base = 1000*ME2_base./MC;
SCC2_base_func = griddedInterpolant(r_mat_1,t_mat_1,k_mat_1,SCC2_base,'spline');

ME2_base_a = -v0_dt_1;
SCC2_base_a = 1000*ME2_base_a./MC;
SCC2_base_a_func = griddedInterpolant(r_mat_1,t_mat_1,k_mat_1,SCC2_base_a,'spline');

file2 = [pwd, '/SCC_mat_Cumu_worst_GrowthNoAmb'];
Model2 = load(file2,'v0');
external_v0_worst = Model2.v0;

ME2_tilt =  external_v0_worst;
SCC2_tilt = 1000*ME2_tilt./MC;
SCC2_tilt_func = griddedInterpolant(r_mat_1,t_mat_1,k_mat_1,SCC2_tilt,'spline');

file1 = [pwd, '/HJB_NonLinGrowth_NoAmb_Sims'];
Model1 = load(file1,'hists2','e_hists2','RE_hists2');

hists2_A = Model1.hists2;
RE_hist = Model1.RE_hists2;
e_hist = Model1.e_hists2;

for time=1:400
    for path=1:1
    SCC_values(time,path) = SCC_func(log((hists2_A(time,1,path))),(hists2_A(time,3,path)),log((hists2_A(time,4,path))));
    SCC1_values(time,path) = SCC1_func(log((hists2_A(time,1,path))),(hists2_A(time,3,path)),log((hists2_A(time,4,path))));
    SCC2_base_values(time,path) = SCC2_base_func(log((hists2_A(time,1,path))),(hists2_A(time,3,path)),log((hists2_A(time,4,path))));
    SCC2_base_a_values(time,path) = SCC2_base_a_func(log((hists2_A(time,1,path))),(hists2_A(time,3,path)),log((hists2_A(time,4,path))));
    SCC2_tilt_values(time,path) = SCC2_tilt_func(log((hists2_A(time,1,path))),(hists2_A(time,3,path)),log((hists2_A(time,4,path))));
    end
end

SCC_total = mean(SCC_values,2);
SCC_private = mean(SCC1_values,2);
SCC2_FK_base = mean(SCC2_base_values,2);
SCC2_V_f = mean(SCC2_base_a_values,2);
SCC2_FK_tilt = mean(SCC2_tilt_values,2);
external = mean(SCC2_V_f,2); 
uncert = mean(SCC2_FK_tilt-SCC2_FK_base,2);


fileID = fopen('Growth_numbers_NoAmb.txt','w');
fprintf(fileID,'external_0yr: %.6f \n',external(1));
fprintf(fileID,'external_50yr: %.6f \n',external(end./2));
fprintf(fileID,'external_100yr: %.6f \n',external(end));
fprintf(fileID,'uncertainty_0yr: %.6f \n',uncert(1));
fprintf(fileID,'uncertainty_50yr: %.6f \n',uncert(end./2));
fprintf(fileID,'uncertainty_100yr: %.6f \n',uncert(end));
fprintf(fileID,'Emissions_0yr: %.6f \n',e_hist(1));
fprintf(fileID,'Emissions_50yr: %.6f \n',e_hist(end./2));
fprintf(fileID,'Emissions_100yr: %.6f \n',e_hist(end));
fprintf(fileID,'Relative Entropy_0yr: %.6f \n',RE_hist(1));
fprintf(fileID,'Relative Entropy_50yr: %.6f \n',RE_hist(end./2));
fprintf(fileID,'Relative Entropy_100yr: %.6f \n',RE_hist(end));
fclose(fileID);


