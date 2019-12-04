%%%%% This file generates results for the Ambiguity Averse model.
% Authors: Mike Barnett, Jieyao Wang
% Last update: Oct 09, 2019
close all
clear all
clc

%% Step 0: set up solver
%%%%% change to location of the solver
% addpath('/mnt/ide0/home/wangjieyao/Climate/FT/')
% addpath('/home/wangjieyao/FT/')
% addpath('/Volumes/homes/FT/')
% addpath('/Volumes/wangjieyao/Climate/FT')

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
Gamma_r = 0.112733407891680;
psi_1 = 0.142857142857143;
power = 2;
gamma_1 =  0.00017675;
gamma_2 =  2.*0.0022;
gamma_2_plus = 2.*0.0197;
sigma_1 = 0;
sigma_2 = 0;
rho_12 = 0;
gamma_bar = 2;
crit = 2;
F0 = 1;
sigma = [sigma_1.^2, rho_12; rho_12, sigma_2.^2];
Sigma = [var_beta_f,0,0;0,sigma_1.^2, rho_12; 0, rho_12, sigma_2.^2];
dee = [gamma_1 + gamma_2.*F0+gamma_2_plus.*(F0-gamma_bar).^2.*(F0>=2),beta_f,beta_f.*F0];
sigma_d = sqrt(dee*Sigma*dee');
xi_d = -1.*(1-kappa);

weight = 1; % weight on nordhaus
bar_gamma_2_plus = weight.*0+(1-weight).*gamma_2_plus;

xi_p = 1.; % ambiguity parameter
nums_count = [];
%% Step 2: solve HJB
r_min = 0;
% r_1 = 5;
% r_2 = 9;
r_max = 9; %r = logR
F_min = 0;
% F_1 = 200;
% F_2 = 1200;
F_max = 4000;
% % nr = 500;
% nr = 300;
% nF = 40;
% r_step_1 = 0.1;
% r_step_2 = 0.025;
% r_step_3 = 0.1;
% F_step_1 = 100;
% F_step_2 = 25;
% F_step_3 = 100;
% 
% r_range_1 = r_min:r_step_1:r_1;
% r_range_2 = r_1:r_step_2:r_2;
% r_range_3 = r_2:r_step_3:r_max;
% 
% F_range_1 = F_min:F_step_1:F_1;
% F_range_2 = F_1:F_step_2:F_2;
% F_range_3 = F_2:F_step_3:F_max;
% 
% r = [r_range_1(1:end-1)'; r_range_2(1:end-1)';r_range_3'];
% t = [F_range_1(1:end-1)'; F_range_2(1:end-1)';F_range_3'];

hr = 0.05;
ht = 25;

r = r_min:hr:r_max;
t = F_min:ht:F_max;
[r_mat,F_mat] = ndgrid(r,t); 

[nr x] = size(r);
[nt x] = size(t);

hr = round(r(2:end)-r(1:end-1), 4);
ht = round(t(2:end)-t(1:end-1), 4);

quadrature = 'legendre';
n = 30;
a = beta_f-5.*sqrt(var_beta_f);
b = beta_f+5.*sqrt(var_beta_f);

tol = 1e-8; % tol
dt = 0.3; % epsilon
v0 = (kappa).*r_mat-beta_f.*F_mat; % initial guess
v1_initial = v0.*ones(size(r_mat));
out = v0;
vold = v0 .* ones(size(v0));
iter = 1;

v0 = v0.*ones(size(v0));
v0 = reshape(out,size(v0));

v0_dt = finiteDifferenceinterp(v0,2,1,ht,'central');
v0_dr = finiteDifferenceinterp(v0,1,1,hr,'central',1e-16);
v0_dtt = finiteDifferenceinterp(v0_dt,2,1,ht,'central');
v0_drr = finiteDifferenceinterp(v0_dr,1,1,hr,'central');
v0_drr(v0_dr<=1e-16)=0;

v0_dk = 1-kappa;


% FOC
B1 = v0_dr-xi_d.*(gamma_1(1)+gamma_2(1)*(F_mat).*beta_f+gamma_2_plus(1).*(F_mat.*beta_f-gamma_bar).^(power-1).*(F_mat>=(crit./beta_f))) ...
    .*beta_f.*exp(r_mat)-v0_dt.*exp(r_mat);  
C1 = -delta.*kappa;
e = -C1./B1;
e_hat = e;

i_k = (alpha.*phi_0.*phi_1.*v0_dk-delta.*(1-kappa))./(delta.*(1-kappa).*phi_1+phi_0.*phi_1.*v0_dk);

% update prob.
a_1 = zeros(size(r_mat));
b_1 = xi_d.*e_hat.*exp(r_mat).*gamma_1;
c_1 = 2.*xi_d.*e_hat.*exp(r_mat).*F_mat.*gamma_2;
lambda_tilde_1 = lambda+c_1./xi_p;
beta_tilde_1 = beta_f-c_1./xi_p./lambda_tilde_1.*beta_f-1./xi_p./lambda_tilde_1.*b_1;
I_1 = a_1-0.5.*log(lambda).*xi_p + 0.5.*log(lambda_tilde_1).*xi_p ...
    +0.5.*lambda.*beta_f.^2.*xi_p -0.5.*lambda_tilde_1.*(beta_tilde_1).^2.*xi_p;
J_1_without_e = xi_d.*(gamma_1.*beta_tilde_1 ...
    +gamma_2.*F_mat.*(beta_tilde_1.^2+1./lambda_tilde_1)).*exp(r_mat) ;
pi_tilde_1 = weight.*exp(-1./xi_p.*I_1);

scale_2_fnc = @(x) exp(-1./xi_p.*xi_d.*(gamma_1.*x ...
    +gamma_2.*x.^2.*F_mat ...
    +gamma_2_plus.*x.*(x.*F_mat-gamma_bar).^(power-1).*((x.*F_mat-gamma_bar)>=0)).*exp(r_mat).*e_hat) ...
    .*normpdf(x,beta_f,sqrt(var_beta_f));
scale_2 = quad_int(scale_2_fnc, [a], [b], n, 'legendre'); % quadrature
q2_tilde_fnc = @(x) exp(-1./xi_p.*xi_d.*(gamma_1.*x ...
    +gamma_2.*x.^2.*F_mat ...
    +gamma_2_plus.*x.*(x.*F_mat-gamma_bar).^(power-1).*((x.*F_mat-gamma_bar)>=0)).*exp(r_mat).*e_hat)./scale_2;
I_2 = -1.*xi_p.*log(scale_2);
J_2_without_e_fnc = @(x) xi_d.*exp(r_mat).*...
    q2_tilde_fnc(x) ...
    .*(gamma_1.*x +gamma_2.*F_mat.*x.^2 ...
    +gamma_2_plus.*x.*(x.*F_mat-gamma_bar).^(power-1).*((x.*F_mat-gamma_bar)>=0)) ...
    .*normpdf(x,beta_f,sqrt(var_beta_f));
J_2_without_e = quad_int(J_2_without_e_fnc, [a], [b], n,'legendre'); % quadrature
J_2_with_e = J_2_without_e.*e_hat;
pi_tilde_2 = (1-weight).*exp(-1./xi_p.*I_2);
pi_tilde_1_norm = pi_tilde_1./(pi_tilde_1+pi_tilde_2);
pi_tilde_2_norm = 1-pi_tilde_1_norm;  
    
expec_e_sum = (pi_tilde_1_norm.*(J_1_without_e) ...
    +pi_tilde_2_norm.*(J_2_without_e));
B1 = v0_dr-v0_dt.*exp(r_mat)-expec_e_sum;  
C1 = -delta.*kappa;
e = -C1./B1; % FOC 
e_star = e; 

J_1 = J_1_without_e.*e_star;
J_2 = J_2_without_e.*e_star;
I_term = -1.*xi_p.*log(pi_tilde_1+pi_tilde_2);
R_1 = 1./xi_p.*(I_1 - J_1);     
R_2 = 1./xi_p.*(I_2 - J_2);
drift_distort = (pi_tilde_1_norm.*J_1 ...
    +pi_tilde_2_norm.*J_2);   
if (weight == 0 || weight == 1)
    RE = pi_tilde_1_norm.*R_1 + pi_tilde_2_norm.*R_2;
else
    RE = pi_tilde_1_norm.*R_1 + pi_tilde_2_norm.*R_2 ...
        +pi_tilde_1_norm.*(log(pi_tilde_1_norm./weight)) ...
        +pi_tilde_2_norm.*(log(pi_tilde_2_norm./(1-weight)));        
end
RE_total = 1.*xi_p.*RE;

% inputs for solver
A = -delta.*ones(size(r_mat));
B_r = -e_star-0.5.*(sigma_r.^2);
B_t = e_star.*exp(r_mat);
C_rr = 0.5.*sigma_r.^2.*ones(size(r_mat));
C_tt = zeros(size(r_mat));

D = delta.*kappa.*log(e_star)+delta.*kappa.*r_mat ...
    + delta.*(1-kappa).*(log(alpha-i_k)) ...
    + drift_distort + RE_total;    

stateSpace = [r_mat(:), F_mat(:)]; 
model      = {};
model.A    = A(:); 
model.B    = [B_r(:), B_t(:)];
model.C    = [C_rr(:), C_tt(:)];
model.D    = D(:);
model.v0   = v0(:);
model.dt   = dt;


out = solveCGNatural(stateSpace, model);
out_comp = reshape(out,size(v0)).*ones(size(r_mat));
    
disp(['Error: ', num2str(max(max(max(abs(out_comp - v1_initial)))))])
v0 = v0.*ones(size(v0));
v0 = reshape(out,size(v0));


pde_error_old = (out_comp - v1_initial);
pde_error = ones(size(r_mat));

s0 = '/HJB_2D_VG_low_150';
filename = [pwd, s0];
load(filename, 'v0'); % save HJB solution
dt = 0.2
% iter = 0
xi_p = 1./175
% vold = out_comp + 100;

while (max(max(max(abs(out_comp - vold))))) > tol % check for convergence
   tic
   vold = v0 .* ones(size(v0));
   pde_error_old = pde_error;
   
    stateSpace = [r_mat(:), F_mat(:)]; 
    model      = {};
    model.A    = A(:); 
    model.B    = [B_r(:), B_t(:)];
    model.C    = [C_rr(:), C_tt(:)];
    model.D    = D(:);
    model.v0   = v0(:);
    model.dt   = dt;
    
    out = solveCGNatural(stateSpace, model);
    out_comp = reshape(out,size(v0)).*ones(size(r_mat));

    iter = iter+1;
    if mod(iter, 2000) == 0
        s0 = '/HJB_2D_VG_low_';
        filename = [pwd, s0, num2str(1/xi_p), '_temp'];
        save(filename); % save HJB solution
    end
    Diff = max(max(max(abs(out_comp - vold))));
    disp(['Diff: ', num2str(max(max(max(abs(out_comp - vold))))),' Number of Iters: ',num2str(iter)])

    v0 = v0.*ones(size(v0));
    v0 = reshape(out,size(v0));
    

% v0_dt = finiteDifference(v0,2,1,ht,'central');
% v0_dr = finiteDifference(v0,1,1,hr,'central',1e-16);
% v0_dtt = finiteDifference(v0,2,2,ht,'central');
% v0_drr = finiteDifference(v0,1,2,hr,'central');
% v0_drr(v0_dr<=1e-16)=0;

v0_dt = finiteDifferenceinterp(v0,2,1,ht,'central');
v0_dr = finiteDifferenceinterp(v0,1,1,hr,'central',1e-16);
v0_dtt = finiteDifferenceinterp(v0_dt,2,1,ht,'central');
v0_drr = finiteDifferenceinterp(v0_dr,1,1,hr,'central');
v0_drr(v0_dr<=1e-16)=0;


v0_dk = 1-kappa;

    e_hat = e_star;

i_k = (alpha.*phi_0.*phi_1.*v0_dk-delta.*(1-kappa))./(delta.*(1-kappa).*phi_1+phi_0.*phi_1.*v0_dk);

    a_1 = zeros(size(r_mat));
    b_1 = xi_d.*e_hat.*exp(r_mat).*gamma_1;
    c_1 = 2.*xi_d.*e_hat.*exp(r_mat).*F_mat.*gamma_2;
    lambda_tilde_1 = lambda+c_1./xi_p;
    beta_tilde_1 = beta_f-c_1./xi_p./lambda_tilde_1.*beta_f-1./xi_p./lambda_tilde_1.*b_1;
    I_1 = a_1-0.5.*log(lambda).*xi_p + 0.5.*log(lambda_tilde_1).*xi_p ...
        +0.5.*lambda.*beta_f.^2.*xi_p -0.5.*lambda_tilde_1.*(beta_tilde_1).^2.*xi_p;
    R_1 = 1./xi_p.*(I_1-(a_1+b_1.*beta_tilde_1+0.5.*c_1.*(beta_tilde_1).^2+0.5.*c_1./lambda_tilde_1));
    J_1_without_e = xi_d.*(gamma_1.*beta_tilde_1 ...
        +gamma_2.*F_mat.*(beta_tilde_1.^2+1./lambda_tilde_1)).*exp(r_mat) ;

    pi_tilde_1 = weight.*exp(-1./xi_p.*I_1);


    scale_2_fnc = @(x) exp(-1./xi_p.*xi_d.*(gamma_1.*x ...
        +gamma_2.*x.^2.*F_mat ...
        +gamma_2_plus.*x.*(x.*F_mat-gamma_bar).^(power-1).*((x.*F_mat-gamma_bar)>=0)).*exp(r_mat).*e_hat) ...
        .*normpdf(x,beta_f,sqrt(var_beta_f));
    scale_2 = quad_int(scale_2_fnc, [a], [b], n, 'legendre');
    q2_tilde_fnc = @(x) exp(-1./xi_p.*xi_d.*(gamma_1.*x ...
        +gamma_2.*x.^2.*F_mat ...
        +gamma_2_plus.*x.*(x.*F_mat-gamma_bar).^(power-1).*((x.*F_mat-gamma_bar)>=0)).*exp(r_mat).*e_hat)./scale_2;
    I_2 = -1.*xi_p.*log(scale_2);
    J_2_without_e_fnc = @(x) xi_d.*exp(r_mat).*...
        q2_tilde_fnc(x) ...
        .*(gamma_1.*x +gamma_2.*F_mat.*x.^2 ...
        +gamma_2_plus.*x.*(x.*F_mat-gamma_bar).^(power-1).*((x.*F_mat-gamma_bar)>=0)) ...
        .*normpdf(x,beta_f,sqrt(var_beta_f));
    J_2_without_e = quad_int(J_2_without_e_fnc, [a], [b], n,'legendre');
    J_2_with_e = J_2_without_e.*e_hat;   
    R_2 = 1./xi_p.*(I_2 - J_2_with_e);
    pi_tilde_2 = (1-weight).*exp(-1./xi_p.*I_2);
    pi_tilde_1_norm = pi_tilde_1./(pi_tilde_1+pi_tilde_2);
    pi_tilde_2_norm = 1-pi_tilde_1_norm;  

    expec_e_sum = (pi_tilde_1_norm.*(J_1_without_e) ...
        +pi_tilde_2_norm.*(J_2_without_e));
    B1 = v0_dr-v0_dt.*exp(r_mat)-expec_e_sum;  
    C1 = -delta.*kappa;
    e = -C1./B1;
    e_star = e; 
    
    J_1 = J_1_without_e.*e_star;
    J_2 = J_2_without_e.*e_star;
    I_term = -1.*xi_p.*log(pi_tilde_1+pi_tilde_2);
    drift_distort = (pi_tilde_1_norm.*J_1 ...
        +pi_tilde_2_norm.*J_2);
    R_1 = 1./xi_p.*(I_1 - J_1);
    R_2 = 1./xi_p.*(I_2 - J_2);
    if (weight == 0 || weight == 1)
        RE = pi_tilde_1_norm.*R_1 + pi_tilde_2_norm.*R_2;
    else
        RE = pi_tilde_1_norm.*R_1 + pi_tilde_2_norm.*R_2 ...
            +pi_tilde_1_norm.*(log(pi_tilde_1_norm./weight)) ...
            +pi_tilde_2_norm.*(log(pi_tilde_2_norm./(1-weight)));        
    end
    RE_total = 1.*xi_p.*RE;

    A = -delta.*ones(size(r_mat));
    B_r = -e_star-0.5.*(sigma_r.^2);
    B_t = e_star.*exp(r_mat);
    C_rr = 0.5.*sigma_r.^2.*ones(size(r_mat));
    C_tt = zeros(size(r_mat));

    D = delta.*kappa.*log(e_star)+delta.*kappa.*r_mat ...
        + delta.*(1-kappa).*(log(alpha-i_k)) ...
        + drift_distort + RE_total;   



    pde_error = A.*v0+B_r.*v0_dr+B_t.*v0_dt+C_rr.*v0_drr+C_tt.*v0_dtt+D;
    disp(['PDE Error: ', num2str(max(max(max(abs(pde_error)))))]) % check equation
    
    disp(['PDE Error Diff: ', num2str(max(max(max(abs(pde_error-pde_error_old)))))]) % check equation


    if (max(max(abs(out_comp - vold)))) < tol
        fprintf('PDE Converges');
    end
    
    toc
end

s0 = '/HJB_2D_VG_low_';
filename = [pwd, s0, num2str(1/xi_p)];
save(filename); % save HJB solution

% stop
% load(filename)
% 
% figure()
% plot(nums_count,'LineWidth',2)
% ylabel('number of inner loops')
% xlabel('Nth iteration of outer loop')
% saveas(gcf,'Iterations.png')
% 
% Every_ten_th = [];
% for i = 1:length(nums_count)
%     if mod(i,10)==0
%         Every_ten_th = [Every_ten_th;nums_count(i)];
%     end
% end

% figure()
% plot(Every_ten_th,'LineWidth',2)
% ylabel('number of inner loops')
% xlabel('Every 10th iteration of outer loop')
% saveas(gcf,'Iterations_10.png')



% figure()
% plot(r,e(:,end./2).*(exp(r_mat(:,end./2))),'LineWidth',2)
% hold on
% plot(r,e(:,1).*(exp(r_mat(:,1))),'LineWidth',2)
% % ylabel('number of inner loops')
% % xlabel('Every 10th iteration of outer loop')
% saveas(gcf,'Iterations_10.png')
% 
% 
% figure()
% plot(r,v0_dr(:,end./2),'LineWidth',2)
% % ylabel('number of inner loops')
% % xlabel('Every 10th iteration of outer loop')
% saveas(gcf,'Iterations_10.png')
% 
% 
% figure()
% plot(t,v0_dr(end./2,:),'LineWidth',2)
% % ylabel('number of inner loops')
% % xlabel('Every 10th iteration of outer loop')
% saveas(gcf,'Iterations_10.png')
% 
% 
