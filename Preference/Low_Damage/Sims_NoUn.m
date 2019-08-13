%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc

%%%%% Step 0: Load solution and save filename
filename = [pwd,'/HJB_NonLinPref_Cumu_NoUn'];
load(filename);
filename2 = [pwd,'/HJB_NonLinPref_Cumu_NoUn'];
filename3 = [filename2,'_Sims'];

%%%%% Step 1: Set up simulation
T = 100; % 100 years
pers = 4*T; % quarterly
dt = T/pers;
nDims = 5;
its = 1;

%%%%% Step 2: Create functions
efunc = griddedInterpolant(r_mat,t_mat,k_mat,e,'spline');
if weight == 1
    ffunc = griddedInterpolant(r_mat,t_mat,k_mat,f,'spline');
else
    ffunc = griddedInterpolant(r_mat,t_mat,k_mat,f,'linear');
end
i_kfunc = griddedInterpolant(r_mat,t_mat,k_mat,i_k,'spline');

v_drfunc = griddedInterpolant(r_mat,t_mat,k_mat,v0_dr,'spline');
v_dtfunc = griddedInterpolant(r_mat,t_mat,k_mat,v0_dt,'spline');
v_dkfunc = griddedInterpolant(r_mat,t_mat,k_mat,v0_dk,'spline');
v_func = griddedInterpolant(r_mat,t_mat,k_mat,out_comp,'spline');

pi_tilde_1func = griddedInterpolant(r_mat,t_mat,k_mat,pi_tilde_1_norm,'spline');
pi_tilde_2func = griddedInterpolant(r_mat,t_mat,k_mat,1-pi_tilde_1_norm,'spline');

e_func = @(x) efunc(log(x(:,1)),x(:,3),log(x(:,2)));
f_func = @(x) max(ffunc(log(x(:,1)),x(:,3),log(x(:,2))),0);
i_k_func = @(x) i_kfunc(log(x(:,1)),x(:,3),log(x(:,2)));

v_dr_func = @(x) v_drfunc(log(x(:,1)),x(:,3),log(x(:,2)));
v_dt_func = @(x) v_dtfunc(log(x(:,1)),x(:,3),log(x(:,2)));
v_dk_func = @(x) v_dkfunc(log(x(:,1)),x(:,3),log(x(:,2)));
v_func = @(x) v_func(log(x(:,1)),x(:,3),log(x(:,2)));
bar_gamma_2_plus = (1-weight).*gamma_2_plus;
pi_tilde_1_func = @(x) pi_tilde_1func(log(x(:,1)),x(:,3),log(x(:,2)));
pi_tilde_2_func = @(x) pi_tilde_2func(log(x(:,1)),x(:,3),log(x(:,2)));

base_model_drift_func = @(x) exp(r_mat).*e...
        .*(gamma_1.*x +gamma_2.*t_mat.*x.^2 ...
        +bar_gamma_2_plus.*x.*(x.*t_mat-f_bar).^(power-1).*((x.*t_mat-f_bar)>=0))...
        .*normpdf(x,beta_f,sqrt(var_beta_f)); ...
base_model_drift = quad_int(base_model_drift_func, [a], [b], n,'legendre');

mean_nordhaus = beta_tilde_1;
lambda_tilde_nordhaus = lambda_tilde_1;
nordhaus_model_drift = (gamma_1.*mean_nordhaus...
    +gamma_2.*(1./lambda_tilde_nordhaus+mean_nordhaus.^2).*t_mat).*exp(r_mat).*e;

weitzman_model_drift_func = @(x) exp(r_mat).*e.*...
        q2_tilde_fnc(x) ...
        .*(gamma_1.*x +gamma_2.*t_mat.*x.^2 ...
        +gamma_2_plus.*x.*(x.*t_mat-f_bar).^(power-1).*((x.*t_mat-f_bar)>=0))...
        .*normpdf(x,beta_f,sqrt(var_beta_f)); ...
weitzman_model_drift = quad_int(weitzman_model_drift_func, [a], [b], n,'legendre');

nordhaus_driftfunc = griddedInterpolant(r_mat,t_mat,k_mat,nordhaus_model_drift,'spline');
weitzman_driftfunc = griddedInterpolant(r_mat,t_mat,k_mat,weitzman_model_drift,'spline');
base_driftfunc = griddedInterpolant(r_mat,t_mat,k_mat,base_model_drift,'spline');
nordhaus_drift_func = @(x) nordhaus_driftfunc(log(x(:,1)),x(:,3),log(x(:,2)));
weitzman_drift_func = @(x) weitzman_driftfunc(log(x(:,1)),x(:,3),log(x(:,2)));
base_drift_func = @(x) base_driftfunc(log(x(:,1)),x(:,3),log(x(:,2)));

%%%%% Step 3: Set up initial values
R_0 = 650;
K_0 = 80/A_O;
T_0 = (870-580);
initial_val = [R_0 K_0 T_0];
D_0_base = base_drift_func(initial_val);
D_0_tilted = pi_tilde_1_func(initial_val).*nordhaus_drift_func(initial_val)...
    +(1-pi_tilde_1_func(initial_val)).*weitzman_drift_func(initial_val);

%%%%% Step 4: Create function handle for drifts and vols

muR = @(x) -e_func(x)+Gamma_r.*f_func(x).^Theta_r;
muK = @(x) (Alpha + Gamma.*log(1+i_k_func(x)./Theta));
muT = @(x) e_func(x).*x(:,1);
muD_base = @(x) base_drift_func(x);
muD_tilted = @(x) pi_tilde_1_func(x).*nordhaus_drift_func(x)...
    +(1-pi_tilde_1_func(x)).*weitzman_drift_func(x);

sigmaR = @(x) [zeros(size(x(:,1:5)))];
sigmaK = @(x) [zeros(size(x(:,1:5)))];
sigmaT = @(x) [zeros(size(x(:,1:5)))];
sigmaD = @(x) [zeros(size(x(:,1:5)))];

%%%%% Step 5: Set bounds
R_max = exp(r_max);
K_max = exp(k_max);
T_max = t_max;
D_max = 5.0;

R_min = exp(r_min);
K_min = exp(k_min);
T_min = t_min;
D_min = -5; 

upperBounds = [R_max,K_max,T_max,D_max,D_max];
lowerBounds = [R_min,K_min,T_min,D_min,D_min];

%%%%% Step 6: Initialize values and simulate trajectories
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

v_dr_hist2 = zeros(pers,1);
v_dt_hist2 = zeros(pers,1);
v_dk_hist2 = zeros(pers,1);
v_hist2 = zeros(pers,1);

hist2(1,:) = [R_0,K_0,T_0,D_0_base,D_0_tilted];
e_hist2(1) =  e_func(hist2(1,:)).*hist2(1,1);
i_k_hist2(1) =  i_k_func(hist2(1,:)).*hist2(1,2);
f_hist2(1) =  f_func(hist2(1,:)).*hist2(1,1);
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
hist2(j,3) = max(min(hist2(j-1,3) + muT(hist2(j-1,:)) * dt + sigmaT(hist2(j-1,:))* shock', upperBounds(:,3)), lowerBounds(:,3)); 
hist2(j,4) = max(min(hist2(j-1,4) + muD_base(hist2(j-1,:)) * dt + sigmaD(hist2(j-1,:))* shock', upperBounds(:,4)), lowerBounds(:,4)); 
hist2(j,5) = max(min(hist2(j-1,5) + muD_tilted(hist2(j-1,:)) * dt + sigmaD(hist2(j-1,:))* shock', upperBounds(:,5)), lowerBounds(:,5)); 

e_hist2(j) = e_func(hist2(j-1,:)).*hist2(j-1,1);
i_k_hist2(j) = i_k_func(hist2(j-1,:)).*hist2(j-1,2);
f_hist2(j) =  f_func(hist2(j-1,:)).*hist2(j-1,1);

v_dr_hist2(j) =  v_dr_func(hist2(j-1,:));
v_dt_hist2(j) =  v_dt_func(hist2(j-1,:));
v_dk_hist2(j) =  v_dk_func(hist2(j-1,:));
v_hist2(j) =  v_func(hist2(j-1,:));
end

hists2(:,:,iters) = hist2;
e_hists2(:,iters) = e_hist2;
i_k_hists2(:,iters) = i_k_hist2;
f_hists2(:,iters) = f_hist2; 

v_dr_hists2(:,iters) =  v_dr_hist2;
v_dt_hists2(:,iters) =  v_dt_hist2;
v_dk_hists2(:,iters) =  v_dk_hist2;
v_hists2(:,iters) =  v_hist2;
end

%%%%% Step 7: save results 
save(filename3);
