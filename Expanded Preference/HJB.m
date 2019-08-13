%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Main program %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all
clc
% stop
directory = pwd;
load([directory, '/climate_pars']);
theta = 4500;

weight = 0.25;% on nordhaus

s1 = num2str(1./theta,4);
s1 = strrep(s1,'.','');
s2 = num2str(kappa,4);
s2 = strrep(s2,'.','');
s = strcat('/HJB_NonLinPref_A=',s1,'_M=',s2);
s0 = '/HJB_NonLinPref_Cumu_check_exp';
filename2 = [pwd, s0];

addpath('/Users/johnwilson/Google Drive/BBH Climate Work/Draft_ResultsAndCode_04062019/Code/06052019/20190724')

sigma_o = sigma_k;
indic = 0;

%%% beta_f
beta_f = mean(par_lambda_McD);
var_beta_f = var(par_lambda_McD);
par_beta_f = par_lambda_McD;

lambda = 1./var_beta_f;

%%%%%%%%%%%%%%%%

% stop
xi_d = -1.*(1-alpha);


%%%%%%%%%%%%%%%%
r_min = 0;
r_max = 9;
t_min = 0;
t_max = 4000;
k_min = 0;
k_max = 9;

% eps = 1;
tol = 1e-10;

% nr = 30;
% nt = 100;
% nk = 15;
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

%%% time step
dt = 0.5;

n = 30;

gamma1 = -[-0.00017674741764803637, -0.00017674741764803637, ...
    -0.00017674741764803637, -0.00017674741764803637];
gamma2 = -[-0.0021810708664979964, -0.00277323, ...
    -0.00411616, -0.0071769];

a = beta_f-5.*sqrt(var_beta_f);
b = beta_f+5.*sqrt(var_beta_f);

step_val1 = round(length(r)./5);
step_val2 = round(length(t)./5);
step_val3 = round(length(k)./5);


sigma_1 = 0;
sigma_2 = 0;
rho_12 = 0;
F0 = 1;
sigma = [sigma_1.^2, rho_12; rho_12, sigma_2.^2];
Sigma = [var_beta_f,0,0;0,sigma_1.^2, rho_12; 0, rho_12, sigma_2.^2];
dee = [mean(gamma1) + mean(gamma2).*F0,beta_f,beta_f.*F0];
sigma_d = sqrt(dee*Sigma*dee');


v0 = (alpha).*r_mat+(1-alpha).*k_mat-beta_f.*t_mat; %% initial guess
v1_initial = v0.*ones(size(r_mat)); %% initial guess
out = v0;

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
    v0_dtt(:,2:end-1,:) = (1./(ht.^2)).*(v0(:,3:end,:)+v0(:,1:end-2,:)-2.*v0(:,2:end-1,:));
    v0_dtt(:,end,:) = (1./(ht.^2)).*(v0(:,end,:)+v0(:,end-2,:)-2.*v0(:,end-1,:));
    v0_dtt(:,1,:) = (1./(ht.^2)).*(v0(:,3,:)+v0(:,1,:)-2.*v0(:,2,:));

    v0_dkk = zeros(size(v0));
    v0_dkk(:,:,2:end-1) = (1./(hk.^2)).*(v0(:,:,3:end)+v0(:,:,1:end-2)-2.*v0(:,:,2:end-1));
    v0_dkk(:,:,end) = (1./(hk.^2)).*(v0(:,:,end)+v0(:,:,end-2)-2.*v0(:,:,end-1));
    v0_dkk(:,:,1) = (1./(hk.^2)).*(v0(:,:,3)+v0(:,:,1)-2.*v0(:,:,2));

    %%%%%%% FOC_all combinations    
  
    B1 = v0_dr-xi_d.*(gamma1(4)+gamma2(4)*(t_mat).*beta_f) ...
        .*beta_f.*exp(r_mat)-v0_dt.*exp(r_mat);  
    C1 = -delta.*alpha;
          e = -C1./B1;
          e_hat = e;

    Acoeff = exp(r_mat-k_mat);
    Bcoeff = delta.*(1-alpha)./(exp(-r_mat+k_mat).*v0_dr.*Gamma_r.*0.5)...
                        +v0_dk.*Gamma./(exp(-r_mat+k_mat).*v0_dr.*Gamma_r.*0.5);
    Ccoeff = -A_O - Theta;
    f = ((-Bcoeff+sqrt(Bcoeff.^2 - 4.*Acoeff.*Ccoeff))./(2.*Acoeff)).^(2);
    
    i_k = (v0_dk.*Gamma./(exp(-r_mat+k_mat).*v0_dr.*Gamma_r.*0.5)).*(f.^(0.5))-Theta;

    
    %%% lower damage model 1:
    a_1 = zeros(size(r_mat));
    b_1 = xi_d.*e_hat.*exp(r_mat).*gamma1(1);
    c_1 = 2.*xi_d.*e_hat.*exp(r_mat).*t_mat.*gamma2(1);
    lambda_tilde_1 = lambda+c_1.*theta;
    beta_tilde_1 = beta_f-c_1.*theta./lambda_tilde_1.*beta_f-theta./lambda_tilde_1.*b_1;
    I_1 = a_1-0.5.*log(lambda)./theta + 0.5.*log(lambda_tilde_1)./theta ...
        +0.5.*lambda.*beta_f.^2./theta -0.5.*lambda_tilde_1.*(beta_tilde_1).^2./theta;
    R_1 = theta.*(I_1-(a_1+b_1.*beta_tilde_1+0.5.*c_1.*(beta_tilde_1).^2+0.5.*c_1./lambda_tilde_1));
    J_1_without_e = xi_d.*(gamma1(1).*beta_tilde_1 ...
        +gamma2(1).*t_mat.*(beta_tilde_1.^2+1./lambda_tilde_1)).*exp(r_mat) ;
    
    pi_tilde_1 = weight.*exp(-theta.*I_1);
    
    %%% lower damage model 2:
    a_2 = zeros(size(r_mat));
    b_2 = xi_d.*e_hat.*exp(r_mat).*gamma1(2);
    c_2 = 2.*xi_d.*e_hat.*exp(r_mat).*t_mat.*gamma2(2);
    lambda_tilde_2 = lambda+c_2.*theta;
    beta_tilde_2 = beta_f-c_2.*theta./lambda_tilde_2.*beta_f-theta./lambda_tilde_2.*b_2;
    I_2 = a_2-0.5.*log(lambda)./theta + 0.5.*log(lambda_tilde_2)./theta ...
        +0.5.*lambda.*beta_f.^2./theta -0.5.*lambda_tilde_2.*(beta_tilde_2).^2./theta;
    R_2 = theta.*(I_2-(a_2+b_2.*beta_tilde_2+0.5.*c_2.*(beta_tilde_2).^2+0.5.*c_2./lambda_tilde_2));
    J_2_without_e = xi_d.*(gamma1(2).*beta_tilde_2 ...
        +gamma2(2).*t_mat.*(beta_tilde_2.^2+1./lambda_tilde_2)).*exp(r_mat) ;
    
    pi_tilde_2 = weight.*exp(-theta.*I_2);
    
    
    %%% lower damage model 3:
    a_3 = zeros(size(r_mat));
    b_3 = xi_d.*e_hat.*exp(r_mat).*gamma1(3);
    c_3 = 2.*xi_d.*e_hat.*exp(r_mat).*t_mat.*gamma2(3);
    lambda_tilde_3 = lambda+c_3.*theta;
    beta_tilde_3 = beta_f-c_3.*theta./lambda_tilde_3.*beta_f-theta./lambda_tilde_3.*b_3;
    I_3 = a_3-0.5.*log(lambda)./theta + 0.5.*log(lambda_tilde_3)./theta ...
        +0.5.*lambda.*beta_f.^2./theta -0.5.*lambda_tilde_3.*(beta_tilde_3).^2./theta;
    R_3 = theta.*(I_3-(a_3+b_3.*beta_tilde_3+0.5.*c_3.*(beta_tilde_3).^2+0.5.*c_3./lambda_tilde_3));
    J_3_without_e = xi_d.*(gamma1(3).*beta_tilde_3 ...
        +gamma2(3).*t_mat.*(beta_tilde_3.^2+1./lambda_tilde_3)).*exp(r_mat) ;
    
    pi_tilde_3 = weight.*exp(-theta.*I_3);
    
    
    %%% lower damage model 4:
    a_4 = zeros(size(r_mat));
    b_4 = xi_d.*e_hat.*exp(r_mat).*gamma1(4);
    c_4 = 2.*xi_d.*e_hat.*exp(r_mat).*t_mat.*gamma2(4);
    lambda_tilde_4 = lambda+c_4.*theta;
    beta_tilde_4 = beta_f-c_4.*theta./lambda_tilde_4.*beta_f-theta./lambda_tilde_4.*b_4;
    I_4 = a_4-0.5.*log(lambda)./theta + 0.5.*log(lambda_tilde_4)./theta ...
        +0.5.*lambda.*beta_f.^2./theta -0.5.*lambda_tilde_4.*(beta_tilde_4).^2./theta;
    R_4 = theta.*(I_4-(a_4+b_4.*beta_tilde_4+0.5.*c_4.*(beta_tilde_4).^2+0.5.*c_4./lambda_tilde_4));
    J_4_without_e = xi_d.*(gamma1(4).*beta_tilde_4 ...
        +gamma2(4).*t_mat.*(beta_tilde_4.^2+1./lambda_tilde_4)).*exp(r_mat) ;
    
    pi_tilde_4 = weight.*exp(-theta.*I_4);

    pi_tilde_sum = pi_tilde_1 + pi_tilde_2 + pi_tilde_3 + pi_tilde_4;
    
    pi_tilde_1_norm = pi_tilde_1./pi_tilde_sum;
    pi_tilde_2_norm = pi_tilde_2./pi_tilde_sum;
    pi_tilde_3_norm = pi_tilde_3./pi_tilde_sum;
    pi_tilde_4_norm = pi_tilde_4./pi_tilde_sum;
    
    expec_e_sum = pi_tilde_1_norm.*(J_1_without_e) ...
        +pi_tilde_2_norm.*(J_2_without_e) ...
        +pi_tilde_3_norm.*(J_3_without_e) ...
        +pi_tilde_4_norm.*(J_4_without_e);
%     
    B1 = v0_dr-v0_dt.*exp(r_mat)-expec_e_sum;  
    C1 = -delta.*alpha;
    e = -C1./B1;
    e_star = e; 

    J_1 = J_1_without_e.*e_star;
    J_2 = J_2_without_e.*e_star;
    J_3 = J_3_without_e.*e_star;
    J_4 = J_4_without_e.*e_star;
    

    I_term = -1./theta.*log(pi_tilde_sum);
    

    R_1 = theta.*(I_1 - J_1);     
    R_2 = theta.*(I_2 - J_2);
    R_3 = theta.*(I_3 - J_3);     
    R_4 = theta.*(I_4 - J_4);
    
    drift_distort = (pi_tilde_1_norm.*J_1 ...
        +pi_tilde_2_norm.*J_2 ...
        +pi_tilde_3_norm.*J_3 ...
        +pi_tilde_4_norm.*J_4);
    
    RE = pi_tilde_1_norm.*R_1 + pi_tilde_2_norm.*R_2 + pi_tilde_3_norm.*R_3 + pi_tilde_4_norm.*R_4 ...
        +pi_tilde_1_norm.*(log(pi_tilde_1_norm./weight)) ...
        +pi_tilde_2_norm.*(log(pi_tilde_2_norm./weight)) ...
        +pi_tilde_3_norm.*(log(pi_tilde_3_norm./weight)) ...
        +pi_tilde_4_norm.*(log(pi_tilde_4_norm./weight));
    RE_total = 1./theta.*RE;

    %%%%%%%% PDE inputs %%%%%%%%%%%

    A = -delta.*ones(size(r_mat));
    B_r = -e_star+Gamma_r.*(f.^Theta_r)-0.5.*(sigma_r.^2);
    B_k = Alpha+Gamma.*log(1+i_k./Theta)-0.5.*(sigma_k.^2);
    B_t = e_star.*exp(r_mat);
    C_rr = 0.5.*sigma_r.^2.*ones(size(r_mat));
    C_kk = 0.5.*sigma_k.^2.*ones(size(r_mat));
    C_tt = zeros(size(r_mat));
        
    D = delta.*alpha.*log(e_star)+delta.*alpha.*r_mat ...
        + delta.*(1-alpha).*(log(A_O-i_k-f.*exp(r_mat-k_mat))+(k_mat)) ...
       + I_term;
   
     %%% create stateSpace and model
    stateSpace = [r_mat(:), t_mat(:), k_mat(:)]; 
    model      = {};
    model.A    = A(:); 
    model.B    = [B_r(:), B_t(:), B_k(:)];
    model.C    = [C_rr(:), C_tt(:), C_kk(:)];
    model.D    = D(:);
    model.v0   = v0(:);
    model.dt   = dt;

    out = solveLinearSys2nd0_CG3(stateSpace, model);
    out_comp = reshape(out,size(v0)).*ones(size(r_mat));
    
    disp(['Error: ', num2str(max(max(max(abs(out_comp - v1_initial) / dt))))])
    iter = 1;
    
    vold = v0 .* ones(size(v0));
    
timeDerivative = [];


while (max(max(max(abs(out_comp - vold) / dt)))) > tol
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
    v0_dtt(:,2:end-1,:) = (1./(ht.^2)).*(v0(:,3:end,:)+v0(:,1:end-2,:)-2.*v0(:,2:end-1,:));
    v0_dtt(:,end,:) = (1./(ht.^2)).*(v0(:,end,:)+v0(:,end-2,:)-2.*v0(:,end-1,:));
    v0_dtt(:,1,:) = (1./(ht.^2)).*(v0(:,3,:)+v0(:,1,:)-2.*v0(:,2,:));

    v0_dkk = zeros(size(v0));
    v0_dkk(:,:,2:end-1) = (1./(hk.^2)).*(v0(:,:,3:end)+v0(:,:,1:end-2)-2.*v0(:,:,2:end-1));
    v0_dkk(:,:,end) = (1./(hk.^2)).*(v0(:,:,end)+v0(:,:,end-2)-2.*v0(:,:,end-1));
    v0_dkk(:,:,1) = (1./(hk.^2)).*(v0(:,:,3)+v0(:,:,1)-2.*v0(:,:,2));
    
    e_hat = e_star;
    
    
    f = ((A_O + Theta).*exp(-r_mat+k_mat).*(v0_dr.*Gamma_r.*Theta_r)...
       ./((v0_dr.*Gamma_r.*Theta_r).*f.^(Theta_r)+...
       (delta.*(1-alpha)+v0_dk.*Gamma))).^(1./(1-Theta_r));

    f = f.*(v0_dr>1e-8);

    i_k = ((v0_dk.*Gamma./(exp(-r_mat+k_mat).*v0_dr.*Gamma_r.*Theta_r)).*(f.^(1-Theta_r))-Theta).*(v0_dr>1e-8)...
        + (v0_dr<=1e-8).*(v0_dk.*Gamma.*A_O - delta.*(1-alpha).*Theta)./(delta.*(1-alpha)+v0_dk.*Gamma);
   
    %%% lower damage model 1:
    a_1 = zeros(size(r_mat));
    b_1 = xi_d.*e_hat.*exp(r_mat).*gamma1(1);
    c_1 = 2.*xi_d.*e_hat.*exp(r_mat).*t_mat.*gamma2(1);
    lambda_tilde_1 = lambda+c_1.*theta;
    beta_tilde_1 = beta_f-c_1.*theta./lambda_tilde_1.*beta_f-theta./lambda_tilde_1.*b_1;
    I_1 = a_1-0.5.*log(lambda)./theta + 0.5.*log(lambda_tilde_1)./theta ...
        +0.5.*lambda.*beta_f.^2./theta -0.5.*lambda_tilde_1.*(beta_tilde_1).^2./theta;
    R_1 = theta.*(I_1-(a_1+b_1.*beta_tilde_1+0.5.*c_1.*(beta_tilde_1).^2+0.5.*c_1./lambda_tilde_1));
    J_1_without_e = xi_d.*(gamma1(1).*beta_tilde_1 ...
        +gamma2(1).*t_mat.*(beta_tilde_1.^2+1./lambda_tilde_1)).*exp(r_mat) ;
    
    pi_tilde_1 = weight.*exp(-theta.*I_1);
    
    %%% lower damage model 2:
    a_2 = zeros(size(r_mat));
    b_2 = xi_d.*e_hat.*exp(r_mat).*gamma1(2);
    c_2 = 2.*xi_d.*e_hat.*exp(r_mat).*t_mat.*gamma2(2);
    lambda_tilde_2 = lambda+c_2.*theta;
    beta_tilde_2 = beta_f-c_2.*theta./lambda_tilde_2.*beta_f-theta./lambda_tilde_2.*b_2;
    I_2 = a_2-0.5.*log(lambda)./theta + 0.5.*log(lambda_tilde_2)./theta ...
        +0.5.*lambda.*beta_f.^2./theta -0.5.*lambda_tilde_2.*(beta_tilde_2).^2./theta;
    R_2 = theta.*(I_2-(a_2+b_2.*beta_tilde_2+0.5.*c_2.*(beta_tilde_2).^2+0.5.*c_2./lambda_tilde_2));
    J_2_without_e = xi_d.*(gamma1(2).*beta_tilde_2 ...
        +gamma2(2).*t_mat.*(beta_tilde_2.^2+1./lambda_tilde_2)).*exp(r_mat) ;
    
    pi_tilde_2 = weight.*exp(-theta.*I_2);
    
    
    %%% lower damage model 3:
    a_3 = zeros(size(r_mat));
    b_3 = xi_d.*e_hat.*exp(r_mat).*gamma1(3);
    c_3 = 2.*xi_d.*e_hat.*exp(r_mat).*t_mat.*gamma2(3);
    lambda_tilde_3 = lambda+c_3.*theta;
    beta_tilde_3 = beta_f-c_3.*theta./lambda_tilde_3.*beta_f-theta./lambda_tilde_3.*b_3;
    I_3 = a_3-0.5.*log(lambda)./theta + 0.5.*log(lambda_tilde_3)./theta ...
        +0.5.*lambda.*beta_f.^2./theta -0.5.*lambda_tilde_3.*(beta_tilde_3).^2./theta;
    R_3 = theta.*(I_3-(a_3+b_3.*beta_tilde_3+0.5.*c_3.*(beta_tilde_3).^2+0.5.*c_3./lambda_tilde_3));
    J_3_without_e = xi_d.*(gamma1(3).*beta_tilde_3 ...
        +gamma2(3).*t_mat.*(beta_tilde_3.^2+1./lambda_tilde_3)).*exp(r_mat) ;
    
    pi_tilde_3 = weight.*exp(-theta.*I_3);
    
    
    %%% lower damage model 4:
    a_4 = zeros(size(r_mat));
    b_4 = xi_d.*e_hat.*exp(r_mat).*gamma1(4);
    c_4 = 2.*xi_d.*e_hat.*exp(r_mat).*t_mat.*gamma2(4);
    lambda_tilde_4 = lambda+c_4.*theta;
    beta_tilde_4 = beta_f-c_4.*theta./lambda_tilde_4.*beta_f-theta./lambda_tilde_4.*b_4;
    I_4 = a_4-0.5.*log(lambda)./theta + 0.5.*log(lambda_tilde_4)./theta ...
        +0.5.*lambda.*beta_f.^2./theta -0.5.*lambda_tilde_4.*(beta_tilde_4).^2./theta;
    R_4 = theta.*(I_4-(a_4+b_4.*beta_tilde_4+0.5.*c_4.*(beta_tilde_4).^2+0.5.*c_4./lambda_tilde_4));
    J_4_without_e = xi_d.*(gamma1(4).*beta_tilde_4 ...
        +gamma2(4).*t_mat.*(beta_tilde_4.^2+1./lambda_tilde_4)).*exp(r_mat) ;
    
    pi_tilde_4 = weight.*exp(-theta.*I_4);

    
    
    
    
    pi_tilde_sum = pi_tilde_1 + pi_tilde_2 + pi_tilde_3 + pi_tilde_4;
    
    pi_tilde_1_norm = pi_tilde_1./pi_tilde_sum;
    pi_tilde_2_norm = pi_tilde_2./pi_tilde_sum;
    pi_tilde_3_norm = pi_tilde_3./pi_tilde_sum;
    pi_tilde_4_norm = pi_tilde_4./pi_tilde_sum;
    
    expec_e_sum = pi_tilde_1_norm.*(J_1_without_e) ...
        +pi_tilde_2_norm.*(J_2_without_e) ...
        +pi_tilde_3_norm.*(J_3_without_e) ...
        +pi_tilde_4_norm.*(J_4_without_e);
  
    B1 = v0_dr-v0_dt.*exp(r_mat)-expec_e_sum;  
    C1 = -delta.*alpha;
    e = -C1./B1;
    e_star = e; 

    J_1 = J_1_without_e.*e_star;
    J_2 = J_2_without_e.*e_star;
    J_3 = J_3_without_e.*e_star;
    J_4 = J_4_without_e.*e_star;
    
    I_term = -1./theta.*log(pi_tilde_sum);
    
    drift_distort = (pi_tilde_1_norm.*J_1 ...
        +pi_tilde_2_norm.*J_2);

    drift_distort = (pi_tilde_1_norm.*J_1 ...
        +pi_tilde_2_norm.*J_2 ...
        +pi_tilde_3_norm.*J_3 ...
        +pi_tilde_4_norm.*J_4);
    
    RE = pi_tilde_1_norm.*R_1 + pi_tilde_2_norm.*R_2 + ...
         pi_tilde_3_norm.*R_3 + pi_tilde_4_norm.*R_4 ...
        +pi_tilde_1_norm.*(log(pi_tilde_1_norm./weight)) ...
        +pi_tilde_2_norm.*(log(pi_tilde_2_norm./weight)) ...
        +pi_tilde_3_norm.*(log(pi_tilde_3_norm./weight)) ...
        +pi_tilde_4_norm.*(log(pi_tilde_4_norm./weight));
    RE_total = 1./theta.*RE;

    %%%%%%%% PDE inputs %%%%%%%%%%%

    A = -delta.*ones(size(r_mat));
    B_r = -e_star+Gamma_r.*(f.^Theta_r)-0.5.*(sigma_r.^2);
    B_k = Alpha+Gamma.*log(1+i_k./Theta)-0.5.*(sigma_k.^2);
    B_t = e_star.*exp(r_mat);
    C_rr = 0.5.*sigma_r.^2.*ones(size(r_mat));
    C_kk = 0.5.*sigma_k.^2.*ones(size(r_mat));
    C_tt = zeros(size(r_mat));
        
        D = delta.*alpha.*log(e_star)+delta.*alpha.*r_mat ...
        + delta.*(1-alpha).*(log(A_O-i_k-f.*exp(r_mat-k_mat))+(k_mat)) ...
       + I_term;
    

     %%% create stateSpace and model
    stateSpace = [r_mat(:), t_mat(:), k_mat(:)]; 
    model      = {};
    model.A    = A(:); 
    model.B    = [B_r(:), B_t(:), B_k(:)];
    model.C    = [C_rr(:), C_tt(:), C_kk(:)];
    model.D    = D(:);
    model.v0   = v0(:);
    model.dt   = dt;

    out = solveLinearSys2nd0_CG3(stateSpace, model);
    out_comp = reshape(out,size(v0)).*ones(size(r_mat));
    
    iter = iter+1;
    disp(['Diff: ', num2str(max(max(max(abs(out_comp - vold) / dt)))),' Number of Iters: ',num2str(iter)])
    timeDerivative = [timeDerivative max(max(max(abs(out_comp - vold) / dt)))];
         
    v0 = v0.*ones(size(v0));
    v0 = reshape(out,size(v0));
    
    pde_error = A.*v0+B_r.*v0_dr+B_t.*v0_dt+B_k.*v0_dk+C_rr.*v0_drr+C_kk.*v0_dkk+C_tt.*v0_dtt+D;
    disp(['PDE Error: ', num2str(max(max(max(abs(pde_error)))))])
    
    if (max(max(abs(out_comp - vold) / dt))) < tol
        fprintf('PDE Converges');
    end
    
    toc
end


save(filename2);


