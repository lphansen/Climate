%%%%%%%%%%%%%%%%%%%%%% Feyman Kac for SCC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc

%%%%% Step 0: Add path to location of the solver
addpath('/mnt/ide0/home/wangjieyao/Climate/FT/')
addpath('/home/wangjieyao/FT/')
addpath('/Volumes/homes/FT/')

% for baseline model

%%%%% Step 1: Load solution file
load([pwd, '/HJB_NonLinPref_Cumu']);

%%%%% Step 2: Set up state space
d_min = 0;
d_max = 0.5;
nd = 3;
d = linspace(d_min,d_max,nd)';
hd = d(2) - d(1);
[r_mat_1,t_mat_1,k_mat_1,d_mat_1] = ndgrid(r,t,k,d); 

%%%%% Step 3: Construct PDE drifts
e = repmat(e,[1,1,1,nd]);
f = repmat(f,[1,1,1,nd]);
i_k = repmat(i_k,[1,1,1,nd]);
v0_dr = repmat(v0_dr,[1,1,1,nd]);
v0_dt = repmat(v0_dr,[1,1,1,nd]);
v0_dk = repmat(v0_dk,[1,1,1,nd]);

base_model_drift_func = @(x) ...
    (gamma_1.*x +gamma_2.*t_mat_1.*x.^2 ...
    +bar_gamma_2_plus.*x.*(x.*t_mat_1-f_bar).^(power-1).*((x.*t_mat_1-f_bar)>=0))...
    .*normpdf(x,beta_f,sqrt(var_beta_f)); 
base_model_drift = quad_int(base_model_drift_func, a, b, n,'legendre');
muD_base_1 = base_model_drift;

base_model_flow_func = @(x) ...
    (gamma_2.*x.^2 +bar_gamma_2_plus.*x.^2.*((x.*t_mat_1-f_bar)>=0)).*exp(r_mat_1).*e...
    .*normpdf(x,beta_f,sqrt(var_beta_f)); 
base_model_flow = quad_int(base_model_flow_func, a, b, n,'legendre');
flow_base_1 = base_model_flow;

%%%%% Step 4: Set up inputs for solver
A = -delta.*ones(size(r_mat_1));
B_r = -e+Gamma_r.*(f.^Theta_r)-0.5.*(sigma_r.^2);
B_k = Alpha+Gamma.*log(1+i_k./Theta)-0.5.*(sigma_k.^2);
B_t = e.*exp(r_mat_1);
B_d = muD_base_1.*exp(r_mat_1).*e;
C_rr = 0.5.*sigma_r.^2.*ones(size(r_mat_1));
C_kk = 0.5.*sigma_k.^2.*ones(size(r_mat_1));
C_tt = zeros(size(r_mat_1));
C_dd = 0.5.*sigma_d.^2.*e.^2.*exp(2.*r_mat_1);
D = flow_base_1;
            
%%%%% Step 5: Create stateSpace and model
stateSpace = [r_mat_1(:), t_mat_1(:), k_mat_1(:), d_mat_1(:)]; 
model      = {};
model.A    = A(:); 
model.B    = [B_r(:), B_t(:), B_k(:), B_d(:)];
model.C    = [C_rr(:), C_tt(:), C_kk(:), C_dd(:)];
model.D    = D(:);
model.v0   = r_mat_1(:).*0;
model.dt   = 1.0;
    
%%%%% Step 6: Solve model
out = solveCGNatural_1(stateSpace, model);

%%%%% Step 7: Save results
v0 = reshape(out,size(r_mat_1));
filename2 = [pwd, '/SCC_mat_Cumu_base'];
save(filename2)

%%%%% Step 8: Repeat for tilted model:
clear all
close all
clc

%%%%% Load solution file
load([pwd, '/HJB_NonLinPref_Cumu']);

%%%%% Set up state space
d_min = 0;
d_max = 0.5;
nd = 3;
d = linspace(d_min,d_max,nd)';
hd = d(2) - d(1);
[r_mat_1,t_mat_1,k_mat_1,d_mat_1] = ndgrid(r,t,k,d); 

%%%%% Construct PDE drifts
e = repmat(e,[1,1,1,nd]);
f = repmat(f,[1,1,1,nd]);
i_k = repmat(i_k,[1,1,1,nd]);
v0_dr = repmat(v0_dr,[1,1,1,nd]);
v0_dt = repmat(v0_dr,[1,1,1,nd]);
v0_dk = repmat(v0_dk,[1,1,1,nd]);

mean_nordhaus = beta_tilde_1;
lambda_tilde_nordhaus = lambda_tilde_1;
nordhaus_model_drift = (gamma_1.*mean_nordhaus...
    +gamma_2.*t_mat_1.*(1./lambda_tilde_nordhaus+mean_nordhaus.^2));

scale_2_1_fnc = @(x) exp(-theta.*xi_d.*(gamma_1.*x ...
    +gamma_2.*x.^2.*t_mat_1 ...
    +gamma_2_plus.*x.*(x.*t_mat_1-f_bar).^(power-1).*((x.*t_mat_1-f_bar)>=0)).*exp(r_mat_1).*e) ...
    .*normpdf(x,beta_f,sqrt(var_beta_f));
scale_2_1 = quad_int(scale_2_1_fnc, [a], [b], n, 'legendre');
q2_tilde_1_fnc = @(x) exp(-theta.*xi_d.*(gamma_1.*x ...
    +gamma_2.*x.^2.*t_mat_1 ...
    +gamma_2_plus.*x.*(x.*t_mat_1-f_bar).^(power-1).*((x.*t_mat_1-f_bar)>=0)).*exp(r_mat_1).*e)./scale_2_1;
weitzman_model_drift_func = @(x) q2_tilde_1_fnc(x) ...
    .*(gamma_1.*x +gamma_2.*t_mat_1.*x.^2 ...
    +gamma_2_plus.*x.*(x.*t_mat_1 - f_bar).^(power-1).*((x.*t_mat_1-f_bar)>=0)) ...
    .*normpdf(x,beta_f,sqrt(var_beta_f));
weitzman_model_drift = quad_int(weitzman_model_drift_func, [a], [b], n, 'legendre');

nordhaus_model_flow = (gamma_2.*(1./lambda_tilde_nordhaus+mean_nordhaus.^2)).*exp(r_mat_1).*e;

weitzman_model_flow_func = @(x) q2_tilde_1_fnc(x) ...
    .*(gamma_2.*x.^2 +gamma_2_plus.*x.^2.*((x.*t_mat_1-f_bar)>=0)).*exp(r_mat_1).*e ...
    .*normpdf(x,beta_f,sqrt(var_beta_f));
weitzman_model_flow = quad_int(weitzman_model_flow_func, [a], [b], n, 'legendre');
    
I_1 = -0.5.*log(lambda)./theta + 0.5.*log(lambda_tilde_1)./theta ...
    +0.5.*lambda.*beta_f.^2./theta -0.5.*lambda_tilde_1.*(beta_tilde_1).^2./theta;
pi_tilde_1 = weight.*exp(-theta.*I_1);  
I_2 = -1./theta.*log(scale_2_1);
pi_tilde_2 = (1-weight).*exp(-theta.*I_2);
pi_tilde_1_norm = pi_tilde_1./(pi_tilde_1+pi_tilde_2);
pi_tilde_2_norm = 1-pi_tilde_1_norm;

muD_tilted = pi_tilde_1_norm.*nordhaus_model_drift ...
    +(1-pi_tilde_1_norm).*weitzman_model_drift;
muD_tilted_1 = muD_tilted;

flow_tilted = pi_tilde_1_norm.*nordhaus_model_flow ...
    +(1-pi_tilde_1_norm).*weitzman_model_flow;
flow_tilted_1 = flow_tilted;

%%%%% Set up inputs for solver
A = -delta.*ones(size(r_mat_1));
B_r = -e+Gamma_r.*(f.^Theta_r)-0.5.*(sigma_r.^2);
B_k = Alpha+Gamma.*log(1+i_k./Theta)-0.5.*(sigma_k.^2);
B_t = e.*exp(r_mat_1);
B_d = muD_tilted_1.*e.*exp(r_mat_1);
C_rr = 0.5.*sigma_r.^2.*ones(size(r_mat_1));
C_kk = 0.5.*sigma_k.^2.*ones(size(r_mat_1));
C_tt = zeros(size(r_mat_1));
C_dd = 0.5.*sigma_d.^2.*e.^2.*exp(2.*r_mat_1);
D = flow_tilted_1;

%%%%% Create stateSpace and model  
stateSpace = [r_mat_1(:), t_mat_1(:), k_mat_1(:), d_mat_1(:)]; 
model      = {};
model.A    = A(:); 
model.B    = [B_r(:), B_t(:), B_k(:), B_d(:)];
model.C    = [C_rr(:), C_tt(:), C_kk(:), C_dd(:)];
model.D    = D(:);
model.v0   = r_mat_1(:).*0;
model.dt   = 1.0;

%%%%% Solve model
out = solveCGNatural_1(stateSpace, model);

%%%%% Save results
v0 = reshape(out,size(r_mat_1));
filename2 = [pwd, '/SCC_mat_Cumu_worst'];
save(filename2)