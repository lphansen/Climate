%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Main program %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all
clc
% stop
directory = pwd;
addpath('/Users/johnwilson/Google Drive/BBH Climate Work/Draft_ResultsAndCode_04062019/Code/06052019/20190724')
load('climate_pars');

s1 = num2str(1./theta,4);
s1 = strrep(s1,'.','');
s2 = num2str(kappa,4);
s2 = strrep(s2,'.','');
s = strcat('/HJB_NonLinGrowth_A=',s1,'_M=',s2);
s0 = '/HJB_NonLinGrowth_NoAmb';
filename2 = [pwd, s0];

% theta = 0.001;
theta = 200;
RobDrift = 0;

sigma_o = sigma_k;
indic = 0;

beta_f = mean(par_lambda_McD);
var_beta_f = var(par_lambda_McD);
par_beta_f = par_lambda_McD;

t_bar = 13;

lambda = 1./var_beta_f;

%%%%%%%%%%%%%%%%
%%% 2 state variables: r,t
r_min = 0;
r_max = 9;
t_min = 0;
t_max = 750;
k_min = 0;
k_max = 9;

tol = 1e-8;


%--------------------------Kushner-Dupuis Set-Up---------------------------
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
sigs = 5;

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

% stop




% Defined domain of integration in z
a = beta_f-sigs.*sqrt(var_beta_f);
b = beta_f+sigs.*sqrt(var_beta_f);
[gam1,w1] = quad_points_hermite(3) ;
gamm1 = sqrt(2) * 1 * gam1 + 0 ;
[gam2,w2]=quad_points_hermite(3);
gamm2 = sqrt(2) * 1 * gam2 + 0 ;
w1 = w1 ./ sum(w1);
w2 = w2 ./ sum(w2);

[A,B] = meshgrid(gamm1,gamm2);
c=cat(2,A',B');
gs=reshape(c,[],2);

[A,B] = meshgrid(w1,w2);
c=cat(2,A',B');
ws=reshape(c,[],2);
ws = prod(ws, 2);

weight1 = ws(1);
weight2 = ws(2);
weight3 = ws(3);
weight4 = ws(4);
weight5 = ws(5);
weight6 = ws(6);
weight7 = ws(7);
weight8 = ws(8);
weight9 = ws(9);
% sum(ws)

A=(chol(sigma))';

gs = [gamma_1;gamma_2] + A*gs';
gs = gs';

arg_maxs = -gs(:, 1) ./ gs(:, 2);
maxs = -gs(:, 1) .* arg_maxs - .5*gs(:, 2) .* arg_maxs.^2;

gammas = zeros(9, 3);
gammas(:,1) = -maxs;
gammas(:,2:end) = -gs;

% stop

dee = [gamma_1 + gamma_2,beta_f,beta_f];
sigma_d = sqrt(dee*Sigma*dee');

step_val1 = round(length(r)./5);
step_val2 = round(length(t)./5);
step_val3 = round(length(k)./5);


v0 = (alpha).*r_mat+(1-alpha).*k_mat;

% v0 = (alpha).*r_mat+(1-alpha).*k_mat;
% -t_mat.^2;
v1_initial = v0.*ones(size(r_mat));
out = v0;
vold = v0 .* ones(size(v0));
%%% for first iter
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
    
    v0_dtt = zeros(size(v0));
    v0_dtt(:,2:end-1) = (1./(ht.^2)).*(v0(:,3:end)+v0(:,1:end-2)-2.*v0(:,2:end-1));
    v0_dtt(:,end) = (1./(ht.^2)).*(v0(:,end)+v0(:,end-2)-2.*v0(:,end-1));
    v0_dtt(:,1) = (1./(ht.^2)).*(v0(:,3)+v0(:,1)-2.*v0(:,2));

    v0_dkk = zeros(size(v0));
    v0_dkk(:,:,2:end-1) = (1./(hk.^2)).*(v0(:,:,3:end)+v0(:,:,1:end-2)-2.*v0(:,:,2:end-1));
    v0_dkk(:,:,end) = (1./(hk.^2)).*(v0(:,:,end)+v0(:,:,end-2)-2.*v0(:,:,end-1));
    v0_dkk(:,:,1) = (1./(hk.^2)).*(v0(:,:,3)+v0(:,:,1)-2.*v0(:,:,2));
   
   
    B1 = v0_dr-v0_dt.*exp(r_mat);  
    C1 = delta.*alpha;
    e = C1./B1;
    e_hat = e;

    Acoeff = exp(r_mat-k_mat);
    Bcoeff = delta.*(1-alpha)./(exp(-r_mat+k_mat).*v0_dr.*Gamma_r.*0.5)...
                        +v0_dk.*Gamma./(exp(-r_mat+k_mat).*v0_dr.*Gamma_r.*0.5);
    Ccoeff = -A_O - Theta;
    f = ((-Bcoeff+sqrt(Bcoeff.^2 - 4.*Acoeff.*Ccoeff))./(2.*Acoeff)).^(2);
    
    i_k = (v0_dk.*Gamma./(exp(-r_mat+k_mat).*v0_dr.*Gamma_r.*0.5)).*(f.^(0.5))-Theta;



gamma0 = -gammas(:,1);
gamma1 = -gammas(:,2);
gamma2 = -gammas(:,3);



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
    
    
    
% % % % % % % % % % % % % % % % % % % % % % % % % %     
% % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % 

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
     
    
    B1 = v0_dr-v0_dt.*exp(r_mat);  
    C1 = delta.*alpha;
    e = C1./B1;
    e_star = e;

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
    
    %%%
    
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


    %%%%%%%% PDE inputs %%%%%%%%%%%
    
    A = -delta.*ones(size(r_mat));
    B_r = -e_star+Gamma_r.*(f.^Theta_r)-0.5.*(sigma_r.^2);
    B_k = Alpha+Gamma.*log(1+i_k./Theta)-0.5.*(sigma_o.^2);
    B_t = e_star.*exp(r_mat);
    C_rr = 0.5.*sigma_r.^2.*ones(size(r_mat));
    C_kk = 0.5.*sigma_k.^2.*ones(size(r_mat));
    C_tt = zeros(size(r_mat));

    D = delta.*alpha.*log(e_star)+delta.*alpha.*r_mat ...
        + delta.*(1-alpha).*(log(A_O-i_k-f.*exp(r_mat-k_mat))+(k_mat)) ...
        +I_term;
            
     %%% create stateSpace and model
    stateSpace = [r_mat(:), t_mat(:), k_mat(:)]; %% vectorizing state space; this would correspond as the domain argument in the mex file
    model      = {};
    model.A    = A(:); 
    model.B    = [B_r(:), B_t(:), B_k(:)];
    model.C    = [C_rr(:), C_tt(:), C_kk(:)];
    model.D    = D(:);
    model.v0   = v0(:);
    model.dt   = dt;
    %%% solve linear through mex fnction
    out = solveLinearSys2nd0_CG3(stateSpace, model);
    out_comp = reshape(out,size(v0)).*ones(size(r_mat));
    
    %%% check convergence
    
    
    disp(['Error: ', num2str(max(max(max(max(abs(out_comp - v1_initial) / dt)))))])
%     error_plot = zeros(size(r),20);
    iter = 1;
    
        v0 = v0.*ones(size(v0));
    v0 = reshape(out,size(v0));
    
   
    
timeDerivative = [];
while (max(max(max(max(abs(out_comp - vold) / dt))))) > tol
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
    
         f = ((A_O + Theta).*exp(-r_mat+k_mat).*(v0_dr.*Gamma_r.*Theta_r)...
       ./((v0_dr.*Gamma_r.*Theta_r).*f.^(Theta_r)+...
       (delta.*(1-alpha)+v0_dk.*Gamma))).^(1./(1-Theta_r));
    f = f.*(v0_dr>1e-8);
    
    i_k = ((v0_dk.*Gamma./(exp(-r_mat+k_mat).*v0_dr.*Gamma_r.*Theta_r)).*(f.^(1-Theta_r))-Theta).*(v0_dr>1e-8)...
        + (v0_dr<=1e-8).*(v0_dk.*Gamma.*A_O - delta.*(1-alpha).*Theta)./(delta.*(1-alpha)+v0_dk.*Gamma);
  
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
    
    
    
% % % % % % % % % % % % % % % % % % % % % % % % % %     
% % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % 

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
     
    
    B1 = v0_dr-v0_dt.*exp(r_mat);  
    C1 = delta.*alpha;
    e = C1./B1;
    e_star = e;

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
    
    %%%
    
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


    %%%%%%%% PDE inputs %%%%%%%%%%%
    
    A = -delta.*ones(size(r_mat));
    B_r = -e_star+Gamma_r.*(f.^Theta_r)-0.5.*(sigma_r.^2);
    B_k = Alpha+Gamma.*log(1+i_k./Theta)-0.5.*(sigma_o.^2);
    B_t = e_star.*exp(r_mat);
    C_rr = 0.5.*sigma_r.^2.*ones(size(r_mat));
    C_kk = 0.5.*sigma_k.^2.*ones(size(r_mat));
    C_tt = zeros(size(r_mat));
  
    D = delta.*alpha.*log(e_star)+delta.*alpha.*r_mat ...
        + delta.*(1-alpha).*(log(A_O-i_k-f.*exp(r_mat-k_mat))+(k_mat)) ...
        +I_term;
    
    
      %%% create stateSpace and model
    stateSpace = [r_mat(:), t_mat(:), k_mat(:)]; %% vectorizing state space; this would correspond as the domain argument in the mex file
    model      = {};
    model.A    = A(:); 
    model.B    = [B_r(:), B_t(:), B_k(:)];
    model.C    = [C_rr(:), C_tt(:), C_kk(:)];
    model.D    = D(:);
    model.v0   = v0(:);
    model.dt   = dt;
    %%% solve linear through mex fnction
     out = solveLinearSys2nd0_CG3(stateSpace, model);
    out_comp = reshape(out,size(v0)).*ones(size(r_mat));
    
    
    %%% check convergence
    iter = iter+1;
    disp(['Diff: ', num2str(max(max(max(max(abs(out_comp - vold) / dt))))),' Number of Iters: ',num2str(iter)])
    timeDerivative = [timeDerivative max(max(max(max(abs(out_comp - vold) / dt))))];
    
    v0 = v0.*ones(size(v0));
    v0 = reshape(out,size(v0));
    

    pde_error = A.*v0+B_r.*v0_dr+B_t.*v0_dt+B_k.*v0_dk ...
                     +C_rr.*v0_drr+C_kk.*v0_dkk+C_tt.*v0_dtt+D;
    disp(['PDE Error: ', num2str(max(max(max(max(abs(pde_error))))))])
    if (max(max(max(max(abs(out_comp - vold) / dt))))) < tol
        fprintf('PDE Converges');
    end
    
    toc
end

save(filename2);

figure
plot(exp(k),exp(k).*squeeze(i_k(1.*step_val1,3*step_val2,:)),'-b','LineWidth',1)
hold on
plot(exp(k),exp(k).*squeeze(i_k(3.*step_val1,3*step_val2,:)),'-r','LineWidth',1)
plot(exp(k),exp(k).*squeeze(i_k(5.*step_val1-1,3*step_val2,:)),'-g','LineWidth',1)
title(sprintf('Invest'))
xlabel('K')
ylabel('I_K')
% ylim([0.0 2.0])
% xlim([0.0 5.0])
legend(sprintf('R=%1.1f',exp(r(step_val1))),sprintf('R=%1.1f',exp(r(3*step_val1))),...
    sprintf('R=%1.1f',exp(r(5*step_val1-1))),...
    'location','northwest')
set(findall(gcf,'type','axes'),'fontsize',16,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',16,'fontWeight','bold')
print('I_k_r','-dpng')

figure
plot(exp(k),exp(k).*squeeze(i_k(3.*step_val1,1*step_val2,:)),'-b','LineWidth',1)
hold on
plot(exp(k),exp(k).*squeeze(i_k(3.*step_val1,3*step_val2,:)),'-r','LineWidth',1)
plot(exp(k),exp(k).*squeeze(i_k(3.*step_val1,5*step_val2,:)),'-g','LineWidth',1)
title(sprintf('Invest'))
xlabel('K')
ylabel('I_K')
% ylim([0.0 2.0])
% xlim([0.0 5.0])
legend(sprintf('T=%1.1f',t(step_val2)),sprintf('T=%1.1f',t(3*step_val2)),...
    sprintf('T=%1.1f',t(5*step_val2)),...
    'location','northwest')
set(findall(gcf,'type','axes'),'fontsize',16,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',16,'fontWeight','bold')
print('I_k_t','-dpng')


figure
plot(exp(r),exp(r).*squeeze(f(:,1*step_val2,3.*step_val3)),'-b','LineWidth',1)
hold on
plot(exp(r),exp(r).*squeeze(f(:,3*step_val2,3.*step_val3)),'-r','LineWidth',1)
plot(exp(r),exp(r).*squeeze(f(:,5*step_val2,3.*step_val3)),'-g','LineWidth',1)
hold off
title(sprintf('F'))
xlabel('R')
ylabel('F')
% ylim([0.0 2.0])
% xlim([0.0 5.0])
legend(sprintf('T=%1.1f',t(step_val2)),sprintf('T=%1.1f',t(3*step_val2)),...
    sprintf('T=%1.1f',t(5*step_val2)),...
    'location','northwest')
set(findall(gcf,'type','axes'),'fontsize',16,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',16,'fontWeight','bold')
print('F_r_t','-dpng')

figure
plot(exp(r),exp(r).*squeeze(f(:,3*step_val2,1.*step_val3)),'-b','LineWidth',1)
hold on
plot(exp(r),exp(r).*squeeze(f(:,3*step_val2,3.*step_val3)),'-r','LineWidth',1)
plot(exp(r),exp(r).*squeeze(f(:,3*step_val2,5.*step_val3)),'-g','LineWidth',1)
hold off
title(sprintf('F'))
xlabel('R')
ylabel('F')
% ylim([0.0 2.0])
% xlim([0.0 5.0])
legend(sprintf('K=%1.1f',exp(k(step_val3))),sprintf('K=%1.1f',exp(k(3*step_val3))),...
    sprintf('K=%1.1f',exp(k(5*step_val3))),...
    'location','northeast')
set(findall(gcf,'type','axes'),'fontsize',16,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',16,'fontWeight','bold')
print('F_r_k','-dpng')


figure
plot(exp(r),exp(r).*squeeze(e(:,1*step_val2,3.*step_val3)),'-b','LineWidth',1)
hold on
plot(exp(r),exp(r).*squeeze(e(:,3*step_val2,3.*step_val3)),'-r','LineWidth',1)
plot(exp(r),exp(r).*squeeze(e(:,5*step_val2,3.*step_val3)),'-g','LineWidth',1)
hold off
title(sprintf('E'))
xlabel('R')
ylabel('E')
% ylim([0.0 2.0])
% xlim([0.0 5.0])
legend(sprintf('T=%1.1f',t(step_val2)),sprintf('T=%1.1f',t(3*step_val2)),...
    sprintf('T=%1.1f',t(5*step_val2)),...
    'location','northwest')
set(findall(gcf,'type','axes'),'fontsize',16,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',16,'fontWeight','bold')
print('E_r_t','-dpng')

figure
plot(exp(r),exp(r).*squeeze(e(:,3*step_val2,1.*step_val3)),'-b','LineWidth',1)
hold on
plot(exp(r),exp(r).*squeeze(e(:,3*step_val2,3.*step_val3)),'-r','LineWidth',1)
plot(exp(r),exp(r).*squeeze(e(:,3*step_val2,5.*step_val3)),'-g','LineWidth',1)
hold off
title(sprintf('E'))
xlabel('R')
ylabel('E')
% ylim([0.0 2.0])
% xlim([0.0 5.0])
legend(sprintf('K=%1.1f',exp(k(step_val3))),sprintf('K=%1.1f',exp(k(3*step_val3))),...
    sprintf('K=%1.1f',exp(k(5*step_val3))),...
    'location','northwest')
set(findall(gcf,'type','axes'),'fontsize',16,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',16,'fontWeight','bold')
print('E_r_k','-dpng')








