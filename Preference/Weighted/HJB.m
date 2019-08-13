%%%%%%%%%%%%%%%%%%%%%%%%%%%% Value Function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all
clc

%% Model set-up
%%%%% Step 0: add path to location of the solver
addpath('/home/wangjieyao/FT/')
%%%%% Step 1: Set up model-specific parameters
load([pwd, '/climate_pars']); % common parameters
theta = 4500; % ambiguity parameter
power = 2;
gamma_1 =  0.00017675;
gamma_2 =  2.*0.0022;
gamma_2_plus = 2.*0.0197;
sigma_1 = 0;
sigma_2 = 0;
rho_12 = 0;
weight = 0.5; % weight on nordhaus
bar_gamma_2_plus = weight.*0+(1-weight).*gamma_2_plus;
beta_f = mean(par_lambda_McD); % from TCRE_MacDougall Et Al (2017)
var_beta_f = var(par_lambda_McD);
par_beta_f = par_lambda_McD;
lambda = 1./var_beta_f;
f_bar = 2;
crit = 2;
F0 = 1;
sigma = [sigma_1.^2, rho_12; rho_12, sigma_2.^2];
Sigma = [var_beta_f,0,0;0,sigma_1.^2, rho_12; 0, rho_12, sigma_2.^2];
dee = [gamma_1 + gamma_2.*F0+gamma_2_plus.*(F0-f_bar).^2.*(F0>=2),beta_f,beta_f.*F0];
sigma_d = sqrt(dee*Sigma*dee');
xi_d = -1.*(1-alpha);

%%%%% Step 2: Set up state space
r_min = 0;
r_max = 9; %r = logR
t_min = 0;
t_max = 4000;
k_min = 0;
k_max = 9; %k = logK
nr = 30;
nt = 40;
nk = 25;
r = linspace(r_min,r_max,nr)';
t = linspace(t_min,t_max,nt)';
k = linspace(k_min,k_max,nk)';
hr = r(2) - r(1);
hk = k(2) - k(1);
ht = t(2) - t(1);
[r_mat,t_mat,k_mat] = ndgrid(r,t,k); 

%%%%% Step 3: Set up quadrature
quadrature = 'legendre';
n = 30;
a = beta_f-5.*sqrt(var_beta_f);
b = beta_f+5.*sqrt(var_beta_f);

%%%%% Step 4: Set up PDE structure
tol = 1e-10; % tolerance
dt = 0.5; % epsilon
v0 = (alpha).*r_mat+(1-alpha).*k_mat-beta_f.*t_mat; % initial guess
v1_initial = v0.*ones(size(r_mat));
out = v0;
vold = v0 .* ones(size(v0));
timeDerivative = [];
iter = 1;

%% PDE set-up
%%%%% Step 1: Finite difference scheme
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

%%%%% Step 2: First order conditions   
B1 = v0_dr-xi_d.*(gamma_1(1)+gamma_2(1)*(t_mat).*beta_f+gamma_2_plus(1).*(t_mat.*beta_f-f_bar).^(power-1).*(t_mat>=(crit./beta_f))) ...
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
    
%%%%% Step 3: Ambiguity Adjusted Probabilities
%%% low damage specification
a_1 = zeros(size(r_mat));
b_1 = xi_d.*e_hat.*exp(r_mat).*gamma_1;
c_1 = 2.*xi_d.*e_hat.*exp(r_mat).*t_mat.*gamma_2;
lambda_tilde_1 = lambda+c_1.*theta;
beta_tilde_1 = beta_f-c_1.*theta./lambda_tilde_1.*beta_f-theta./lambda_tilde_1.*b_1;
I_1 = a_1-0.5.*log(lambda)./theta + 0.5.*log(lambda_tilde_1)./theta ...
    +0.5.*lambda.*beta_f.^2./theta -0.5.*lambda_tilde_1.*(beta_tilde_1).^2./theta;
R_1 = theta.*(I_1-(a_1+b_1.*beta_tilde_1+0.5.*c_1.*(beta_tilde_1).^2+0.5.*c_1./lambda_tilde_1));
J_1_without_e = xi_d.*(gamma_1.*beta_tilde_1 ...
    +gamma_2.*t_mat.*(beta_tilde_1.^2+1./lambda_tilde_1)).*exp(r_mat) ;
pi_tilde_1 = weight.*exp(-theta.*I_1);
%%% high damage specification
scale_2_fnc = @(x) exp(-theta.*xi_d.*(gamma_1.*x ...
    +gamma_2.*x.^2.*t_mat ...
    +gamma_2_plus.*x.*(x.*t_mat-f_bar).^(power-1).*((x.*t_mat-f_bar)>=0)).*exp(r).*e_hat) ...
    .*normpdf(x,beta_f,sqrt(var_beta_f));
scale_2 = quad_int(scale_2_fnc, [a], [b], n, 'legendre');
q2_tilde_fnc = @(x) exp(-theta.*xi_d.*(gamma_1.*x ...
    +gamma_2.*x.^2.*t_mat ...
    +gamma_2_plus.*x.*(x.*t_mat-f_bar).^(power-1).*((x.*t_mat-f_bar)>=0)).*exp(r).*e_hat)./scale_2;
I_2 = -1./theta.*log(scale_2);
J_2_without_e_fnc = @(x) xi_d.*exp(r_mat).*...
    q2_tilde_fnc(x) ...
    .*(gamma_1.*x +gamma_2.*t_mat.*x.^2 ...
    +gamma_2_plus.*x.*(x.*t_mat-f_bar).^(power-1).*((x.*t_mat-f_bar)>=0)) ...
    .*normpdf(x,beta_f,sqrt(var_beta_f));
J_2_without_e = quad_int(J_2_without_e_fnc, [a], [b], n,'legendre');
J_2_with_e = J_2_without_e.*e_hat;
R_2 = theta.*(I_2 - J_2_with_e);
pi_tilde_2 = (1-weight).*exp(-theta.*I_2);
pi_tilde_1_norm = pi_tilde_1./(pi_tilde_1+pi_tilde_2);
pi_tilde_2_norm = 1-pi_tilde_1_norm;  
    
%%%%% Step 4: Update first order condition for e
expec_e_sum = (pi_tilde_1_norm.*(J_1_without_e) ...
    +pi_tilde_2_norm.*(J_2_without_e));

B1 = v0_dr-v0_dt.*exp(r_mat)-expec_e_sum;  
C1 = -delta.*alpha;
e = -C1./B1;
e_star = e; 

%%%%% Step 5: Construct drifts
J_1 = J_1_without_e.*e_star;
J_2 = J_2_without_e.*e_star;
I_term = -1./theta.*log(pi_tilde_1+pi_tilde_2);
R_1 = theta.*(I_1 - J_1);     
R_2 = theta.*(I_2 - J_2);
drift_distort = (pi_tilde_1_norm.*J_1 ...
    +pi_tilde_2_norm.*J_2);   
if (weight == 0 || weight == 1)
    RE = pi_tilde_1_norm.*R_1 + pi_tilde_2_norm.*R_2;
else
    RE = pi_tilde_1_norm.*R_1 + pi_tilde_2_norm.*R_2 ...
        +pi_tilde_1_norm.*(log(pi_tilde_1_norm./weight)) ...
        +pi_tilde_2_norm.*(log(pi_tilde_2_norm./(1-weight)));        
end
RE_total = 1./theta.*RE;

%%%%% Step 6: Set up inputs for solver
A = -delta.*ones(size(r_mat));
B_r = -e_star+Gamma_r.*(f.^Theta_r)-0.5.*(sigma_r.^2);
B_k = Alpha+Gamma.*log(1+i_k./Theta)-0.5.*(sigma_k.^2);
B_t = e_star.*exp(r_mat);
C_rr = 0.5.*sigma_r.^2.*ones(size(r_mat));
C_kk = 0.5.*sigma_k.^2.*ones(size(r_mat));
C_tt = zeros(size(r_mat));

D = delta.*alpha.*log(e_star)+delta.*alpha.*r_mat ...
    + delta.*(1-alpha).*(log(A_O-i_k-f.*exp(r_mat-k_mat))+(k_mat)) ...
    + drift_distort + RE_total;    
   
%%%%% Step 7: Create state space and model inputs
stateSpace = [r_mat(:), t_mat(:), k_mat(:)]; 
model      = {};
model.A    = A(:); 
model.B    = [B_r(:), B_t(:), B_k(:)];
model.C    = [C_rr(:), C_tt(:), C_kk(:)];
model.D    = D(:);
model.v0   = v0(:);
model.dt   = dt;

%%%%% Step 8: Solve the model and compare results
out = solveCGNatural(stateSpace, model);
out_comp = reshape(out,size(v0)).*ones(size(r_mat));
    
disp(['Error: ', num2str(max(max(max(abs(out_comp - v1_initial) / dt))))])
v0 = v0.*ones(size(v0));
v0 = reshape(out,size(v0));

%%%%% Step 9: Loop for iterative scheme: repeat the above process
while (max(max(max(abs(out_comp - vold) / dt)))) > tol
   tic
   vold = v0 .* ones(size(v0));
    %%%%% Finite difference scheme
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

    %%%%% First order conditions 
    e_hat = e_star;

    f = ((A_O + Theta).*exp(-r_mat+k_mat).*(v0_dr.*Gamma_r.*Theta_r)...
       ./((v0_dr.*Gamma_r.*Theta_r).*f.^(Theta_r)+...
       (delta.*(1-alpha)+v0_dk.*Gamma))).^(1./(1-Theta_r));
    f = f.*(v0_dr>1e-8);

    i_k = ((v0_dk.*Gamma./(exp(-r_mat+k_mat).*v0_dr.*Gamma_r.*Theta_r)).*(f.^(1-Theta_r))-Theta).*(v0_dr>1e-8)...
        + (v0_dr<=1e-8).*(v0_dk.*Gamma.*A_O - delta.*(1-alpha).*Theta)./(delta.*(1-alpha)+v0_dk.*Gamma);

    %%%%% Ambiguity Adjusted Probabilities
    %%% low damage specification
    a_1 = zeros(size(r_mat));
    b_1 = xi_d.*e_hat.*exp(r_mat).*gamma_1;
    c_1 = 2.*xi_d.*e_hat.*exp(r_mat).*t_mat.*gamma_2;
    lambda_tilde_1 = lambda+c_1.*theta;
    beta_tilde_1 = beta_f-c_1.*theta./lambda_tilde_1.*beta_f-theta./lambda_tilde_1.*b_1;
    I_1 = a_1-0.5.*log(lambda)./theta + 0.5.*log(lambda_tilde_1)./theta ...
        +0.5.*lambda.*beta_f.^2./theta -0.5.*lambda_tilde_1.*(beta_tilde_1).^2./theta;
    R_1 = theta.*(I_1-(a_1+b_1.*beta_tilde_1+0.5.*c_1.*(beta_tilde_1).^2+0.5.*c_1./lambda_tilde_1));
    J_1_without_e = xi_d.*(gamma_1.*beta_tilde_1 ...
        +gamma_2.*t_mat.*(beta_tilde_1.^2+1./lambda_tilde_1)).*exp(r_mat) ;

    pi_tilde_1 = weight.*exp(-theta.*I_1);

    %%% high damage specification
    scale_2_fnc = @(x) exp(-theta.*xi_d.*(gamma_1.*x ...
        +gamma_2.*x.^2.*t_mat ...
        +gamma_2_plus.*x.*(x.*t_mat-f_bar).^(power-1).*((x.*t_mat-f_bar)>=0)).*exp(r).*e_hat) ...
        .*normpdf(x,beta_f,sqrt(var_beta_f));
    scale_2 = quad_int(scale_2_fnc, [a], [b], n, 'legendre');
    q2_tilde_fnc = @(x) exp(-theta.*xi_d.*(gamma_1.*x ...
        +gamma_2.*x.^2.*t_mat ...
        +gamma_2_plus.*x.*(x.*t_mat-f_bar).^(power-1).*((x.*t_mat-f_bar)>=0)).*exp(r).*e_hat)./scale_2;
    I_2 = -1./theta.*log(scale_2);
    J_2_without_e_fnc = @(x) xi_d.*exp(r_mat).*...
        q2_tilde_fnc(x) ...
        .*(gamma_1.*x +gamma_2.*t_mat.*x.^2 ...
        +gamma_2_plus.*x.*(x.*t_mat-f_bar).^(power-1).*((x.*t_mat-f_bar)>=0)) ...
        .*normpdf(x,beta_f,sqrt(var_beta_f));
    J_2_without_e = quad_int(J_2_without_e_fnc, [a], [b], n,'legendre');
    J_2_with_e = J_2_without_e.*e_hat;   
    R_2 = theta.*(I_2 - J_2_with_e);
    pi_tilde_2 = (1-weight).*exp(-theta.*I_2);
    pi_tilde_1_norm = pi_tilde_1./(pi_tilde_1+pi_tilde_2);
    pi_tilde_2_norm = 1-pi_tilde_1_norm;  

    %%%%% Update first order condition for e    
    expec_e_sum = (pi_tilde_1_norm.*(J_1_without_e) ...
        +pi_tilde_2_norm.*(J_2_without_e));
    B1 = v0_dr-v0_dt.*exp(r_mat)-expec_e_sum;  
    C1 = -delta.*alpha;
    e = -C1./B1;
    e_star = e; 

    %%%%% Construct drifts      
    J_1 = J_1_without_e.*e_star;
    J_2 = J_2_without_e.*e_star;
    I_term = -1./theta.*log(pi_tilde_1+pi_tilde_2);
    drift_distort = (pi_tilde_1_norm.*J_1 ...
        +pi_tilde_2_norm.*J_2);
    R_1 = theta.*(I_1 - J_1);
    R_2 = theta.*(I_2 - J_2);
    if (weight == 0 || weight == 1)
        RE = pi_tilde_1_norm.*R_1 + pi_tilde_2_norm.*R_2;
    else
        RE = pi_tilde_1_norm.*R_1 + pi_tilde_2_norm.*R_2 ...
            +pi_tilde_1_norm.*(log(pi_tilde_1_norm./weight)) ...
            +pi_tilde_2_norm.*(log(pi_tilde_2_norm./(1-weight)));        
    end
    RE_total = 1./theta.*RE;

    %%%%% Set up inputs for solver
    A = -delta.*ones(size(r_mat));
    B_r = -e_star+Gamma_r.*(f.^Theta_r)-0.5.*(sigma_r.^2);
    B_k = Alpha+Gamma.*log(1+i_k./Theta)-0.5.*(sigma_k.^2);
    B_t = e_star.*exp(r_mat);
    C_rr = 0.5.*sigma_r.^2.*ones(size(r_mat));
    C_kk = 0.5.*sigma_k.^2.*ones(size(r_mat));
    C_tt = zeros(size(r_mat));

    D = delta.*alpha.*log(e_star)+delta.*alpha.*r_mat ...
        + delta.*(1-alpha).*(log(A_O-i_k-f.*exp(r_mat-k_mat))+(k_mat)) ...
        + drift_distort + RE_total;    

    %%%%% Create state space and model inputs
    stateSpace = [r_mat(:), t_mat(:), k_mat(:)]; 
    model      = {};
    model.A    = A(:); 
    model.B    = [B_r(:), B_t(:), B_k(:)];
    model.C    = [C_rr(:), C_tt(:), C_kk(:)];
    model.D    = D(:);
    model.v0   = v0(:);
    model.dt   = dt;
    
    %%%%% Solve the model and compare results
    out = solveCGNatural(stateSpace, model);
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

%%% save file
s0 = '/HJB_NonLinPref_Cumu';
filename = [pwd, s0];
save(filename);


