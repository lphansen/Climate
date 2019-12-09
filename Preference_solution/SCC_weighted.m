load('Averse_weighted.mat');
% clear all;

%% SCC load
% base = load('SCC_base_averse_weighted.mat', 'v0');
% base = load('C:\Users\jiamingwang\Dropbox\share with John\Final Code\Preference\Jieyao_weighted_averse_sol\SCC_mat_Cumu_base_4000.mat', 'v0');
% external_v0 = base.v0;

% worst = load('SCC_worst_averse_weighted.mat', 'v0');
% worst = load('C:\Users\jiamingwang\Dropbox\share with John\Final Code\Preference\Jieyao_weighted_averse_sol\SCC_mat_Cumu_worst_4000.mat', 'v0');
worst_jy = load('SCC_base_averse_weighted.mat', 'v0');
% external_v0_worst_jy = worst_jy.external_v0_worst;
external_v0 = worst_jy.v0;

worst = load('SCC_worst_averse_weighted.mat', 'v0');
external_v0_worst = worst.v0;


% model_SCC = load('SCC_averse_weighted.mat', 'external_v0', 'external_v0_worst');
% external_v0 = model_SCC.external_v0;
% external_v0_worst = model_SCC.external_v0_worst;
% save('SCC_averse_low_guess.mat', 'external_v0_worst', 'external_v0');


%% check SCC base

base_model_flow_func = @(x) ...
    (gamma_2.*x.^2 +bar_gamma_2_plus.*x.^2.*((x.*F_mat-gamma_bar)>=0)).*exp(r_mat).*e...
    .*normpdf(x,beta_f,sqrt(var_beta_f)); 
base_model_flow = quad_int(base_model_flow_func, a, b, n,'legendre');
flow_base = base_model_flow;

% inputs for solver
A = -delta.*ones(size(r_mat));
B_r = -e+psi_0.*(j.^psi_1).*exp(psi_1.*(k_mat-r_mat))-0.5.*(sigma_r.^2);
B_k = mu_k+phi_0.*log(1+i_k.*phi_1)-0.5.*(sigma_k.^2);
B_t = e.*exp(r_mat);
C_rr = 0.5.*sigma_r.^2.*ones(size(r_mat));
C_kk = 0.5.*sigma_k.^2.*ones(size(r_mat));
C_tt = zeros(size(r_mat));

D = flow_base;

%% check SCC worst

mean_nordhaus = beta_tilde_1;
lambda_tilde_nordhaus = lambda_tilde_1;

scale_2_1_fnc = @(x) exp(-1./xi_p.*xi_d.*(gamma_1.*x ...
    +gamma_2.*x.^2.*F_mat ...
    +gamma_2_plus.*x.*(x.*F_mat-gamma_bar).^(power-1).*((x.*F_mat-gamma_bar)>=0)).*exp(r_mat).*e) ...
    .*normpdf(x,beta_f,sqrt(var_beta_f));
scale_2_1 = quad_int(scale_2_1_fnc, [a], [b], n, 'legendre');

q2_tilde_1_fnc = @(x) exp(-1./xi_p.*xi_d.*(gamma_1.*x ...
    +gamma_2.*x.^2.*F_mat ...
    +gamma_2_plus.*x.*(x.*F_mat-gamma_bar).^(power-1).*((x.*F_mat-gamma_bar)>=0)).*exp(r_mat).*e)./scale_2_1;

nordhaus_model_flow = (gamma_2.*(1./lambda_tilde_nordhaus+mean_nordhaus.^2)).*exp(r_mat).*e;
weitzman_model_flow_func = @(x) q2_tilde_1_fnc(x) ...
    .*(gamma_2.*x.^2 +gamma_2_plus.*x.^2.*((x.*F_mat-gamma_bar)>=0)).*exp(r_mat).*e ...
    .*normpdf(x,beta_f,sqrt(var_beta_f));
weitzman_model_flow = quad_int(weitzman_model_flow_func, [a], [b], n, 'legendre');

I_1 = -0.5.*log(lambda).*xi_p + 0.5.*log(lambda_tilde_1).*xi_p ...
    +0.5.*lambda.*beta_f.^2.*xi_p -0.5.*lambda_tilde_1.*(beta_tilde_1).^2.*xi_p;
pi_tilde_1 = weight.*exp(-1./xi_p.*I_1);  
I_2 = -1.*xi_p.*log(scale_2_1);
pi_tilde_2 = (1-weight).*exp(-1./xi_p.*I_2);
pi_tilde_1_norm = pi_tilde_1./(pi_tilde_1+pi_tilde_2);
pi_tilde_2_norm = 1-pi_tilde_1_norm;

flow_tilted = pi_tilde_1_norm.*nordhaus_model_flow ...
    +(1-pi_tilde_1_norm).*weitzman_model_flow;

A = -delta.*ones(size(r_mat));
B_r = -e+psi_0.*(j.^psi_1).*exp(psi_1.*(k_mat-r_mat))-0.5.*(sigma_r.^2);
B_k = mu_k+phi_0.*log(1+i_k.*phi_1)-0.5.*(sigma_k.^2);
B_t = e.*exp(r_mat);
C_rr = 0.5.*sigma_r.^2.*ones(size(r_mat));
C_kk = 0.5.*sigma_k.^2.*ones(size(r_mat));
C_tt = zeros(size(r_mat));

D = flow_tilted;

%% run small test

pass_in_worst = external_v0_worst_jy;

stateSpace = [r_mat(:), F_mat(:), k_mat(:)]; 
model      = {};
model.A    = A(:); 
model.B    = [B_r(:), B_t(:), B_k(:)];
model.C    = [C_rr(:), C_tt(:), C_kk(:)];
model.D    = D(:);
model.v0   = v0(:).*0;
model.v1   = pass_in_worst(:);
model.dt   = 1;

out = solveCGNatural1(stateSpace, model);

external_v0_worst = reshape(out,size(r_mat));

%% check pde error

v0 = external_v0_worst;
v0_dr = finiteDifference(v0,1,1,hr,'central');
v0_dt = finiteDifference(v0,2,1,ht,'central');
v0_dk = finiteDifference(v0,3,1,hk,'central');
v0_drr = finiteDifference(v0,1,2,hr,'central');
v0_dtt = finiteDifference(v0,2,2,ht,'central');
v0_dkk = finiteDifference(v0,3,2,hk,'central');

pde_error = A.*v0+B_r.*v0_dr+B_t.*v0_dt+B_k.*v0_dk+C_rr.*v0_drr+C_kk.*v0_dkk+C_tt.*v0_dtt+D;
disp(max(max(max(abs(pde_error)))));
disp(mean(v0,'all'));

%% 
T = 100; % 100 years
pers = 4*T; % quarterly
dt = T/pers;
nDims = 5;
its = 1;

efunc = griddedInterpolant(r_mat,F_mat,k_mat,e,'linear');
% if weight == 1
jfunc = griddedInterpolant(r_mat,F_mat,k_mat,j,'linear');
% else
% %     jfunc = griddedInterpolant(r_mat,F_mat,k_mat,j,'linear');
%     j_psifunc = griddedInterpolant(r_mat,F_mat,k_mat,j.^psi_1,'nearest');
% end
i_kfunc = griddedInterpolant(r_mat,F_mat,k_mat,i_k,'linear');

v_drfunc = griddedInterpolant(r_mat,F_mat,k_mat,v0_dr,'linear');
v_dtfunc = griddedInterpolant(r_mat,F_mat,k_mat,v0_dt,'linear');
v_dkfunc = griddedInterpolant(r_mat,F_mat,k_mat,v0_dk,'linear');
v_func = griddedInterpolant(r_mat,F_mat,k_mat,out_comp,'linear');

pi_tilde_1func = griddedInterpolant(r_mat,F_mat,k_mat,pi_tilde_1_norm,'linear');
pi_tilde_2func = griddedInterpolant(r_mat,F_mat,k_mat,1-pi_tilde_1_norm,'linear');

e_func = @(x) efunc(log(x(:,1)),x(:,3),log(x(:,2)));
% j_func = @(x) max(jfunc(log(x(:,1)),x(:,3),log(x(:,2))),0);
j_func = @(x) max(jfunc(log(x(:,1)),x(:,3),log(x(:,2))),0);
i_k_func = @(x) i_kfunc(log(x(:,1)),x(:,3),log(x(:,2)));

v_dr_func = @(x) v_drfunc(log(x(:,1)),x(:,3),log(x(:,2)));
v_dt_func = @(x) v_dtfunc(log(x(:,1)),x(:,3),log(x(:,2)));
v_dk_func = @(x) v_dkfunc(log(x(:,1)),x(:,3),log(x(:,2)));
v_func = @(x) v_func(log(x(:,1)),x(:,3),log(x(:,2)));
bar_gamma_2_plus = (1-weight).*gamma_2_plus;
pi_tilde_1_func = @(x) pi_tilde_1func(log(x(:,1)),x(:,3),log(x(:,2)));
pi_tilde_2_func = @(x) pi_tilde_2func(log(x(:,1)),x(:,3),log(x(:,2)));

base_model_drift_func = @(x) exp(r_mat).*e...
        .*(gamma_1.*x +gamma_2.*F_mat.*x.^2 ...
        +bar_gamma_2_plus.*x.*(x.*F_mat-gamma_bar).^(power-1).*((x.*F_mat-gamma_bar)>=0))...
        .*normpdf(x,beta_f,sqrt(var_beta_f)); ...
base_model_drift = quad_int(base_model_drift_func, [a], [b], n,'legendre');

mean_nordhaus = beta_tilde_1;
lambda_tilde_nordhaus = lambda_tilde_1;
nordhaus_model_drift = (gamma_1.*mean_nordhaus...
    +gamma_2.*(1./lambda_tilde_nordhaus+mean_nordhaus.^2).*F_mat).*exp(r_mat).*e;

weitzman_model_drift_func = @(x) exp(r_mat).*e.*...
        q2_tilde_fnc(x) ...
        .*(gamma_1.*x +gamma_2.*F_mat.*x.^2 ...
        +gamma_2_plus.*x.*(x.*F_mat-gamma_bar).^(power-1).*((x.*F_mat-gamma_bar)>=0))...
        .*normpdf(x,beta_f,sqrt(var_beta_f)); ...
weitzman_model_drift = quad_int(weitzman_model_drift_func, [a], [b], n,'legendre');

nordhaus_driftfunc = griddedInterpolant(r_mat,F_mat,k_mat,nordhaus_model_drift,'linear');
weitzman_driftfunc = griddedInterpolant(r_mat,F_mat,k_mat,weitzman_model_drift,'linear');
base_driftfunc = griddedInterpolant(r_mat,F_mat,k_mat,base_model_drift,'linear');
nordhaus_drift_func = @(x) nordhaus_driftfunc(log(x(:,1)),x(:,3),log(x(:,2)));
weitzman_drift_func = @(x) weitzman_driftfunc(log(x(:,1)),x(:,3),log(x(:,2)));
base_drift_func = @(x) base_driftfunc(log(x(:,1)),x(:,3),log(x(:,2)));

% initial points
R_0 = 650;
K_0 = 80/alpha;
F_0 = (870-580);
initial_val = [R_0 K_0 F_0];
D_0_base = base_drift_func(initial_val);
D_0_tilted = pi_tilde_1_func(initial_val).*nordhaus_drift_func(initial_val)...
    +(1-pi_tilde_1_func(initial_val)).*weitzman_drift_func(initial_val);

% function handles
muR = @(x) -e_func(x)+psi_0.*(j_func(x).*x(:,2)./x(:,1)).^psi_1;
muK = @(x) (mu_k + phi_0.*log(1+i_k_func(x).*phi_1));
muF = @(x) e_func(x).*x(:,1);
muD_base = @(x) base_drift_func(x);
muD_tilted = @(x) pi_tilde_1_func(x).*nordhaus_drift_func(x)...
    +(1-pi_tilde_1_func(x)).*weitzman_drift_func(x);

sigmaR = @(x) [zeros(size(x(:,1:5)))];
sigmaK = @(x) [zeros(size(x(:,1:5)))];
sigmaF = @(x) [zeros(size(x(:,1:5)))];
sigmaD = @(x) [zeros(size(x(:,1:5)))];

% set bounds
R_max = exp(r_max);
K_max = exp(k_max);
D_max = 5.0;

R_min = exp(r_min);
K_min = exp(k_min);
D_min = -5; 

upperBounds = [R_max,K_max,F_max,D_max,D_max];
lowerBounds = [R_min,K_min,F_min,D_min,D_min];

hists = zeros(pers,nDims,its);
hists2 = hists;
e_hists = zeros(pers,its);
e_hists2 = e_hists;
j_hists = zeros(pers,its);
j_hists2 = j_hists;
i_k_hists = zeros(pers,its);
i_k_hists2 = i_k_hists;
j_temp = j;

for iters = 1:its
    
hist2 = zeros(pers,nDims);
e_hist2 = zeros(pers,1);
i_k_hist2 = zeros(pers,1);
j_hist2 = zeros(pers,1);

v_dr_hist2 = zeros(pers,1);
v_dt_hist2 = zeros(pers,1);
v_dk_hist2 = zeros(pers,1);
v_hist2 = zeros(pers,1);

hist2(1,:) = [R_0,K_0,F_0,D_0_base,D_0_tilted];
e_hist2(1) =  e_func(hist2(1,:)).*hist2(1,1);
i_k_hist2(1) =  i_k_func(hist2(1,:)).*hist2(1,2);
j_hist2(1) =  j_func(hist2(1,:));
J_hist2(1) =  j_func(hist2(1,:)).*hist2(1,2);
v_dr_hist2(1) =  v_dr_func(hist2(1,:));
v_dt_hist2(1) =  v_dt_func(hist2(1,:));
v_dk_hist2(1) =  v_dk_func(hist2(1,:));
v_hist2(1) =  v_func(hist2(1,:));

for j = 2:pers
shock = normrnd(0,sqrt(dt), 1, nDims);
hist2(j,1) = max(min(hist2(j-1,1).*exp((muR(hist2(j-1,:))-0.5.*sum((sigmaR(hist2(j-1,:))).^2) )* dt ...
                                  +sigmaR(hist2(j-1,:))* shock'), upperBounds(:,1)), lowerBounds(:,1));
hist2(j,2) = max(min(hist2(j-1,2).*exp((muK(hist2(j-1,:))-0.5.*sum((sigmaK(hist2(j-1,:))).^2) )* dt ...
                                  +sigmaK(hist2(j-1,:))* shock'), upperBounds(:,2)), lowerBounds(:,2));                              
hist2(j,3) = max(min(hist2(j-1,3) + muF(hist2(j-1,:)) * dt + sigmaF(hist2(j-1,:))* shock', upperBounds(:,3)), lowerBounds(:,3)); 
hist2(j,4) = max(min(hist2(j-1,4) + muD_base(hist2(j-1,:)) * dt + sigmaD(hist2(j-1,:))* shock', upperBounds(:,4)), lowerBounds(:,4)); 
hist2(j,5) = max(min(hist2(j-1,5) + muD_tilted(hist2(j-1,:)) * dt + sigmaD(hist2(j-1,:))* shock', upperBounds(:,5)), lowerBounds(:,5)); 

e_hist2(j) = e_func(hist2(j-1,:)).*hist2(j-1,1);
i_k_hist2(j) = i_k_func(hist2(j-1,:)).*hist2(j-1,2);
j_hist2(j) =  j_func(hist2(j-1,:));
J_hist2(1) =  j_func(hist2(j-1,:)).*hist2(j-1,2);

v_dr_hist2(j) =  v_dr_func(hist2(j-1,:));
v_dt_hist2(j) =  v_dt_func(hist2(j-1,:));
v_dk_hist2(j) =  v_dk_func(hist2(j-1,:));
v_hist2(j) =  v_func(hist2(j-1,:));
end

hists2(:,:,iters) = hist2;
e_hists2(:,iters) = e_hist2;
i_k_hists2(:,iters) = i_k_hist2;
j_hists2(:,iters) = j_hist2; 

v_dr_hists2(:,iters) =  v_dr_hist2;
v_dt_hists2(:,iters) =  v_dt_hist2;
v_dk_hists2(:,iters) =  v_dk_hist2;
v_hists2(:,iters) =  v_hist2;
end

j = j_temp;

MC = delta.*(1-kappa)./(alpha.*exp(k_mat)-i_k.*exp(k_mat)-j.*exp(k_mat));
ME = delta.*kappa./(e.*exp(r_mat));
SCC = 1000*ME./MC;
SCC_func = griddedInterpolant(r_mat,F_mat,k_mat,SCC,'linear');

ME1 = (v0_dr.*exp(-r_mat));  
SCC1 = 1000*ME1./MC;
SCC1_func = griddedInterpolant(r_mat,F_mat,k_mat,SCC1,'linear');

ME2_base =  (1-kappa).*external_v0;
SCC2_base = 1000*ME2_base./MC;
SCC2_base_func = griddedInterpolant(r_mat,F_mat,k_mat,SCC2_base,'linear');

V_d_baseline_func = @(x) xi_d...
         .*(gamma_1.*x +gamma_2.*F_mat.*x.^2 ...
        +bar_gamma_2_plus.*x.*(x.*F_mat-gamma_bar).^(power-1).*((x.*F_mat-gamma_bar)>=0)) ...
        .*normpdf(x,beta_f,sqrt(var_beta_f));
V_d_baseline = quad_int(V_d_baseline_func, [a], [b], n,'legendre');

ME2b = - V_d_baseline;
SCC2_V_d_baseline = 1000*ME2b./MC;
SCC2_V_d_baseline_func = griddedInterpolant(r_mat,F_mat,k_mat,SCC2_V_d_baseline,'linear');


ME2_tilt =  (1-kappa).*external_v0_worst;
SCC2_tilt = 1000*ME2_tilt./MC;
SCC2_tilt_func = griddedInterpolant(r_mat,F_mat,k_mat,SCC2_tilt,'linear');

ME2b = - expec_e_sum.*exp(-r_mat);
SCC2_V_d_tilt = 1000*ME2b./MC;
SCC2_V_d_tilt_func = griddedInterpolant(r_mat,F_mat,k_mat,SCC2_V_d_tilt,'linear');


hists2_A = hists2;

% calculate SCC
for time=1:400
    for path=1:1
    SCC_values(time,path) = SCC_func(log((hists2_A(time,1,path))),(hists2_A(time,3,path)),log((hists2_A(time,2,path))));
    SCC1_values(time,path) = SCC1_func(log((hists2_A(time,1,path))),(hists2_A(time,3,path)),log((hists2_A(time,2,path))));
    SCC2_base_values(time,path) = SCC2_base_func(log((hists2_A(time,1,path))),(hists2_A(time,3,path)),log((hists2_A(time,2,path))));
    SCC2_tilt_values(time,path) = SCC2_tilt_func(log((hists2_A(time,1,path))),(hists2_A(time,3,path)),log((hists2_A(time,2,path))));
    SCC2_V_d_baseline_values(time,path) = SCC2_V_d_baseline_func(log((hists2_A(time,1,path))),(hists2_A(time,3,path)),log((hists2_A(time,2,path))));
    SCC2_V_d_tilt_values(time,path) = SCC2_V_d_tilt_func(log((hists2_A(time,1,path))),(hists2_A(time,3,path)),log((hists2_A(time,2,path))));
    
    end
end

SCC_total = mean(SCC_values,2);
SCC_private = mean(SCC1_values,2);
SCC2_FK_base = mean(SCC2_base_values,2);
SCC2_FK_tilt = mean(SCC2_tilt_values,2);
SCC2_V_d_baseline = mean(SCC2_V_d_baseline_values,2);
SCC2_V_d_tilt = mean(SCC2_V_d_tilt_values,2);

SCC = SCC_total;
SCC1 = SCC_private;
SCC2 = SCC2_FK_base+SCC2_V_d_baseline;
SCC3 = SCC2_V_d_tilt-SCC2_V_d_baseline...
    +SCC2_FK_tilt-SCC2_FK_base;

% save results
s1 = num2str(1.*xi_p,4);
s1 = strrep(s1,'.','');
s_scc = strcat('SCC_AA_num=',s1,'.mat');
save(s_scc,'SCC','SCC1','SCC2','SCC3');

%%
SCC_base = SCC .* exp(hists2(:,4,:));
SCC1_base = SCC1 .* exp(hists2(:,4,:));
SCC2_base = SCC2 .* exp(hists2(:,4,:));
SCC3_base = SCC3 .* exp(hists2(:,4,:));

SCC_tilt = SCC .* exp(hists2(:,5,:));
SCC1_tilt = SCC1 .* exp(hists2(:,5,:));
SCC2_tilt = SCC2 .* exp(hists2(:,5,:));
SCC3_tilt = SCC3 .* exp(hists2(:,5,:));


%% 

years = 0.25:0.25:100;
plot(years, SCC, 'color', [0 0.4470 0.7410], 'LineWidth', 1)
hold on
plot(years, SCC_base, '--', 'color', [0 0.4470 0.7410], 'LineWidth', 1)
hold on
plot(years, SCC_tilt, '-.', 'color', [0 0.4470 0.7410], 'LineWidth', 1)
hold on

plot(years, SCC1, 'color', [0.8500 0.3250 0.0980], 'LineWidth', 1)
hold on
plot(years, SCC1_base, '--', 'color', [0.8500 0.3250 0.0980], 'LineWidth', 1)
hold on
plot(years, SCC1_tilt, '-.', 'color', [0.8500 0.3250 0.0980], 'LineWidth', 1)
hold on

plot(years, SCC2, 'color', [0.9290 0.6940 0.1250], 'LineWidth', 1)
hold on
plot(years, SCC2_base, '--', 'color', [0.9290 0.6940 0.1250], 'LineWidth', 1)
hold on
plot(years, SCC2_tilt, '-.', 'color', [0.9290 0.6940 0.1250], 'LineWidth', 1)
hold on

plot(years, SCC3, 'color', [0.4940 0.1840 0.5560], 'LineWidth', 1)
hold on
plot(years, SCC3_base, '--', 'color', [0.4940 0.1840 0.5560], 'LineWidth', 1)
hold on
plot(years, SCC3_tilt, '-.', 'color', [0.4940 0.1840 0.5560], 'LineWidth', 1)
hold on