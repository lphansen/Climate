%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Main program %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close all
clear all
clc
% stop
directory = pwd;
% filename2 = [directory, '/SCC_mat'];
% load([directory, '/HJB_NonLinPref_Cumu_check']);
load([directory, '/HJB_NonLinPref_Cumu_check_exp']);

addpath('/Users/johnwilson/Google Drive/BBH Climate Work/Draft_ResultsAndCode_04062019/Code/06052019/20190724')

dt = 0.1;
tol = 1e-5;
%%% state space
d_min = 0;
d_max = 0.5;
nd = 3;

d = linspace(d_min,d_max,nd)';

hd = d(2) - d(1);

[r_mat_1,t_mat_1,k_mat_1,d_mat_1] = ndgrid(r,t,k,d); 


e = repmat(e,[1,1,1,nd]);
f = repmat(f,[1,1,1,nd]);
i_k = repmat(i_k,[1,1,1,nd]);
v0_dr = repmat(v0_dr,[1,1,1,nd]);
v0_dt = repmat(v0_dr,[1,1,1,nd]);
v0_dk = repmat(v0_dk,[1,1,1,nd]);

    %%%%%%%%% try different cases %%%%%%%%%%

    mean_base = beta_f;
    lambda_tilde_base = lambda;
    base_model_drift = (mean(gamma1).*mean_base...
        +mean(gamma2).*t_mat_1.*(1./lambda_tilde_base+mean_base.^2));
    muD_base_1 = base_model_drift;
    
    base_model_flow = (mean(gamma2).*(1./lambda_tilde_base+mean_base.^2)).*exp(r_mat_1).*e;
    flow_base_1 = base_model_flow;

v0 = (alpha).*r_mat_1+(1-alpha).*k_mat_1; %% initial guess
v1_initial = v0.*ones(size(r_mat_1)); %% initial guess
out = v0;
    v0 = v0.*ones(size(v0));
    v0 = reshape(out,size(v0));
    %%%%%%%% PDE inputs %%%%%%%%%%%

    %%%%%%%%% baseline model %%%%%%%%%%
    
    A = -delta.*ones(size(r_mat_1));
    B_r = -e+Gamma_r.*(f.^Theta_r)-0.5.*(sigma_r.^2);
    B_k = Alpha+Gamma.*log(1+i_k./Theta)-0.5.*(sigma_o.^2);
    B_t = e.*exp(r_mat_1);
    B_d = muD_base_1.*exp(r_mat_1).*e;
    C_rr = 0.5.*sigma_r.^2.*ones(size(r_mat_1));
    C_kk = 0.5.*sigma_k.^2.*ones(size(r_mat_1));
    C_tt = zeros(size(r_mat_1));
    C_dd = 0.5.*sigma_d.^2.*e.^2.*exp(2.*r_mat_1);
  
    D = flow_base_1;
            
     %%% create stateSpace and model
    stateSpace = [r_mat_1(:), t_mat_1(:), k_mat_1(:), d_mat_1(:)]; 
    model      = {};
    model.A    = A(:); 
    model.B    = [B_r(:), B_t(:), B_k(:), B_d(:)];
    model.C    = [C_rr(:), C_tt(:), C_kk(:), C_dd(:)];
    model.D    = D(:);
    model.v0   = v0(:).*0;
    model.dt   = 1.0;

    out = solveCGNatural_1(stateSpace, model);
    out_comp = reshape(out,size(v0)).*ones(size(r_mat_1));
    

    disp(['Error: ', num2str(max(max(max(max(abs(out_comp - v1_initial) / dt)))))])
%     error_plot = zeros(size(r),20);
    iter = 1;
    
    vold = v0 .* ones(size(v0));
    v0 = v0.*ones(size(v0));
    v0 = reshape(out,size(v0));

directory = pwd;
filename2 = [directory, '/SCC_mat_Cumu_base'];

save(filename2)
% stop
clear all
    
clc
% stop
directory = pwd;

load([directory, '/HJB_NonLinPref_Cumu_check_exp']);

addpath('/Users/johnwilson/Google Drive/BBH Climate Work/Draft_ResultsAndCode_04062019/Code/06052019/20190724')

dt = 0.1;
tol = 1e-5;
%%% state space
d_min = 0;
d_max = 0.5;
nd = 3;

d = linspace(d_min,d_max,nd)';

hd = d(2) - d(1);

[r_mat_1,t_mat_1,k_mat_1,d_mat_1] = ndgrid(r,t,k,d); 

step_val1 = round(length(r)./5);
step_val2 = round(length(t)./5);
step_val3 = round(length(k)./5);
step_val4 = round(length(d)./5);

e = repmat(e,[1,1,1,nd]);
f = repmat(f,[1,1,1,nd]);
i_k = repmat(i_k,[1,1,1,nd]);
beta_tilde_1 = repmat(beta_tilde_1,[1,1,1,nd]);
lambda_tilde_1 = repmat(lambda_tilde_1,[1,1,1,nd]);

beta_tilde_2 = repmat(beta_tilde_2,[1,1,1,nd]);
lambda_tilde_2 = repmat(lambda_tilde_2,[1,1,1,nd]);

beta_tilde_3 = repmat(beta_tilde_3,[1,1,1,nd]);
lambda_tilde_3 = repmat(lambda_tilde_3,[1,1,1,nd]);

beta_tilde_4 = repmat(beta_tilde_4,[1,1,1,nd]);
lambda_tilde_4 = repmat(lambda_tilde_4,[1,1,1,nd]);

    
    %%%%%%%%% try different cases %%%%%%%%%%
    
    mean_nordhaus_1 = beta_tilde_1;
    lambda_tilde_nordhaus_1 = lambda_tilde_1;
    nordhaus_model_drift_1 = (gamma1(1).*mean_nordhaus_1...
        +gamma2(1).*t_mat_1.*(1./lambda_tilde_nordhaus_1+mean_nordhaus_1.^2));
 
    nordhaus_model_flow_1 = (gamma2(1).*(1./lambda_tilde_nordhaus_1+mean_nordhaus_1.^2)).*exp(r_mat_1).*e;
    
    mean_nordhaus_2 = beta_tilde_2;
    lambda_tilde_nordhaus_2 = lambda_tilde_2;
    nordhaus_model_drift_2 = (gamma1(2).*mean_nordhaus_2...
        +gamma2(2).*t_mat_1.*(1./lambda_tilde_nordhaus_2+mean_nordhaus_2.^2));
 
    nordhaus_model_flow_2 = (gamma2(2).*(1./lambda_tilde_nordhaus_2+mean_nordhaus_2.^2)).*exp(r_mat_1).*e;
    
    mean_nordhaus_3 = beta_tilde_3;
    lambda_tilde_nordhaus_3 = lambda_tilde_3;
    nordhaus_model_drift_3 = (gamma1(3).*mean_nordhaus_3...
        +gamma2(3).*t_mat_1.*(1./lambda_tilde_nordhaus_3+mean_nordhaus_3.^2));
 
    nordhaus_model_flow_3 = (gamma2(3).*(1./lambda_tilde_nordhaus_3+mean_nordhaus_3.^2)).*exp(r_mat_1).*e;
    
    mean_nordhaus_4 = beta_tilde_4;
    lambda_tilde_nordhaus_4 = lambda_tilde_4;
    nordhaus_model_drift_4 = (gamma1(4).*mean_nordhaus_4...
        +gamma2(4).*t_mat_1.*(1./lambda_tilde_nordhaus_4+mean_nordhaus_4.^2));
 
    nordhaus_model_flow_4 = (gamma2(4).*(1./lambda_tilde_nordhaus_4+mean_nordhaus_4.^2)).*exp(r_mat_1).*e;

    muD_tilted = pi_tilde_1_norm.*nordhaus_model_drift_1 ...
        +pi_tilde_2_norm.*nordhaus_model_drift_2...
        +pi_tilde_3_norm.*nordhaus_model_drift_3...
        +pi_tilde_4_norm.*nordhaus_model_drift_4;

    muD_tilted_1 = muD_tilted;
    
    flow_tilted = pi_tilde_1_norm.*nordhaus_model_flow_1 ...
        +pi_tilde_2_norm.*nordhaus_model_flow_2 ...
        +pi_tilde_3_norm.*nordhaus_model_flow_3 ...
        +pi_tilde_4_norm.*nordhaus_model_flow_4;
    
    flow_tilted_1 = flow_tilted;


    %%%%%%%%% worst case model %%%%%%%%%%
    v0 = (alpha).*r_mat_1+(1-alpha).*k_mat_1; %% initial guess
v1_initial = v0.*ones(size(r_mat_1)); %% initial guess
out = v0;
    v0 = v0.*ones(size(v0));
    v0 = reshape(out,size(v0));
    
    A = -delta.*ones(size(r_mat_1));
    B_r = -e+Gamma_r.*(f.^Theta_r)-0.5.*(sigma_r.^2);
    B_k = Alpha+Gamma.*log(1+i_k./Theta)-0.5.*(sigma_o.^2);
    B_t = e.*exp(r_mat_1);
    B_d = muD_tilted_1.*e.*exp(r_mat_1);
    C_rr = 0.5.*sigma_r.^2.*ones(size(r_mat_1));
    C_kk = 0.5.*sigma_k.^2.*ones(size(r_mat_1));
    C_tt = zeros(size(r_mat_1));
    C_dd = 0.5.*sigma_d.^2.*e.^2.*exp(2.*r_mat_1);
  
    D = flow_tilted_1;
            
     %%% create stateSpace and model
    stateSpace = [r_mat_1(:), t_mat_1(:), k_mat_1(:), d_mat_1(:)]; 
    model      = {};
    model.A    = A(:); 
    model.B    = [B_r(:), B_t(:), B_k(:), B_d(:)];
    model.C    = [C_rr(:), C_tt(:), C_kk(:), C_dd(:)];
    model.D    = D(:);
    model.v0   = v0(:).*0;
    model.dt   = 1.0;

    out = solveCGNatural_1(stateSpace, model);
    out_comp = reshape(out,size(v0)).*ones(size(r_mat_1));

    disp(['Error: ', num2str(max(max(max(max(abs(out_comp - v1_initial) / dt)))))])

    iter = 1;
    
    vold = v0 .* ones(size(v0));
    v0 = v0.*ones(size(v0));
    v0 = reshape(out,size(v0));

directory = pwd;
filename2 = [directory, '/SCC_mat_Cumu_worst'];

save(filename2)