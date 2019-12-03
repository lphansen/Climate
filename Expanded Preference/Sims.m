
%%% shut down exogenous forcing

clear all
clc

set(0,'defaulttextinterpreter','tex')
set(0, 'defaultAxesTickLabelInterpreter','tex'); 
set(0, 'defaultLegendInterpreter','tex');

directory = pwd;

filename2 = [directory,'/HJB_NonLinPref_Cumu_check_exp'];

load(filename2);

directory = pwd;

filename3 = [filename2,'_Sims'];

R_0 = 650;
K_0 = 666.67;
T_0 = (870-580);
D_0 = mean(gamma1).*T_0+mean(gamma2).*(T_0).^2;

initial_val = [R_0 K_0 T_0];

%%

T = 100; 
pers = 4*T; 
dt = T/pers;


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
pi_tilde_2func = griddedInterpolant(r_mat,t_mat,k_mat,pi_tilde_2_norm,'spline');
pi_tilde_3func = griddedInterpolant(r_mat,t_mat,k_mat,pi_tilde_3_norm,'spline');
pi_tilde_4func = griddedInterpolant(r_mat,t_mat,k_mat,pi_tilde_4_norm,'spline');


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


mean_base = beta_f;
lambda_tilde_base = lambda;
base_model_drift = (mean(gamma1).*mean_base...
    +mean(gamma2).*t_mat.*(1./lambda_tilde_base+mean_base.^2));
damage_base = (mean(gamma1).*mean_base...
    +mean(gamma2).*t_mat.*(1./lambda_tilde_base+mean_base.^2));

mean_nordhaus_1 = beta_tilde_1;
lambda_tilde_nordhaus_1 = lambda_tilde_1;
nordhaus_model_drift_1 = (gamma1(1).*mean_nordhaus_1...
    +gamma2(1).*(1./lambda_tilde_nordhaus_1+mean_nordhaus_1.^2).*t_mat).*exp(r_mat).*e;
damage_1 = (gamma1(1).*mean_nordhaus_1...
    +gamma2(1).*(1./lambda_tilde_nordhaus_1+mean_nordhaus_1.^2).*t_mat);

mean_nordhaus_2 = beta_tilde_2;
lambda_tilde_nordhaus_2 = lambda_tilde_2;
nordhaus_model_drift_2 = (gamma1(2).*mean_nordhaus_2...
    +gamma2(2).*(1./lambda_tilde_nordhaus_2+mean_nordhaus_2.^2).*t_mat).*exp(r_mat).*e;
damage_2 = (gamma1(2).*mean_nordhaus_2...
    +gamma2(2).*(1./lambda_tilde_nordhaus_2+mean_nordhaus_2.^2).*t_mat);

mean_nordhaus_3 = beta_tilde_3;
lambda_tilde_nordhaus_3 = lambda_tilde_3;
nordhaus_model_drift_3 = (gamma1(3).*mean_nordhaus_3...
    +gamma2(3).*(1./lambda_tilde_nordhaus_3+mean_nordhaus_3.^2).*t_mat).*exp(r_mat).*e;
damage_3 = (gamma1(3).*mean_nordhaus_3...
    +gamma2(3).*(1./lambda_tilde_nordhaus_3+mean_nordhaus_3.^2).*t_mat);

mean_nordhaus_4 = beta_tilde_4;
lambda_tilde_nordhaus_4 = lambda_tilde_4;
nordhaus_model_drift_4 = (gamma1(4).*mean_nordhaus_4...
    +gamma2(4).*(1./lambda_tilde_nordhaus_4+mean_nordhaus_4.^2).*t_mat).*exp(r_mat).*e;
damage_4 = (gamma1(4).*mean_nordhaus_4...
    +gamma2(4).*(1./lambda_tilde_nordhaus_4+mean_nordhaus_4.^2).*t_mat);

figure()
plot(t_mat, damage_4)


%%

base_driftfunc = griddedInterpolant(r_mat,t_mat,k_mat,base_model_drift,'spline');
base_drift_func = @(x) base_driftfunc(log(x(:,1)),x(:,3),log(x(:,2)));

nordhaus_driftfunc_1 = griddedInterpolant(r_mat,t_mat,k_mat,nordhaus_model_drift_1,'spline');
nordhaus_drift_func_1 = @(x) nordhaus_driftfunc_1(log(x(:,1)),x(:,3),log(x(:,2)));

nordhaus_driftfunc_2 = griddedInterpolant(r_mat,t_mat,k_mat,nordhaus_model_drift_2,'spline');
nordhaus_drift_func_2 = @(x) nordhaus_driftfunc_2(log(x(:,1)),x(:,3),log(x(:,2)));

nordhaus_driftfunc_3 = griddedInterpolant(r_mat,t_mat,k_mat,nordhaus_model_drift_3,'spline');
nordhaus_drift_func_3 = @(x) nordhaus_driftfunc_3(log(x(:,1)),x(:,3),log(x(:,2)));

nordhaus_driftfunc_4 = griddedInterpolant(r_mat,t_mat,k_mat,nordhaus_model_drift_4,'spline');
nordhaus_drift_func_4 = @(x) nordhaus_driftfunc_4(log(x(:,1)),x(:,3),log(x(:,2)));

D_0_base = base_drift_func(initial_val);
D_0_tilted = pi_tilde_1_func(initial_val).*nordhaus_drift_func_1(initial_val)...
    +pi_tilde_2_func(initial_val).*nordhaus_drift_func_2(initial_val)...
    +pi_tilde_3_func(initial_val).*nordhaus_drift_func_3(initial_val)...
    +pi_tilde_4_func(initial_val).*nordhaus_drift_func_4(initial_val);


%%%Create function handle for the drifts

muR = @(x) -e_func(x)+Gamma_r.*f_func(x).^Theta_r;
muK = @(x) (Alpha + Gamma.*log(1+i_k_func(x)./Theta));
muT = @(x) e_func(x).*x(:,1);
muD = @(x) (mean(gamma1).*beta_f + mean(gamma2).*(x(:,3)).*beta_f.^2)...
    .*e_func(x).*x(:,1);
muD_base = @(x) base_drift_func(x);


muD_tilted = @(x) pi_tilde_1_func(x).*nordhaus_drift_func_1(x)...
    +pi_tilde_2_func(x).*nordhaus_drift_func_2(x)...
    +pi_tilde_3_func(x).*nordhaus_drift_func_3(x)...
    +pi_tilde_4_func(x).*nordhaus_drift_func_4(x);

%%%Create function handles for the vols
sigmaR = @(x) [zeros(size(x(:,1:6)))];
sigmaK = @(x) [zeros(size(x(:,1:6)))];
sigmaT = @(x) [zeros(size(x(:,1:6)))];
sigmaD = @(x) [zeros(size(x(:,1:6)))];

R_max = exp(r_max);
K_max = exp(k_max);
T_max = t_max;
D_max = 5.0;

R_min = exp(r_min);
K_min = exp(k_min);
T_min = t_min;
D_min = -5; 

upperBounds = [R_max,K_max,T_max,D_max,D_max,D_max];
lowerBounds = [R_min,K_min,T_min,D_min,D_min,D_min];

nDims = 6;
its = 1;
% stop

hists = zeros(pers,nDims,its);
hists2 = hists;
e_hists = zeros(pers,its);
e_hists2 = e_hists;
f_hists = zeros(pers,its);
f_hists2 = f_hists;
i_k_hists = zeros(pers,its);
i_k_hists2 = i_k_hists;

tic

for iters = 1:its
% for iters = 1:its
    
hist2 = zeros(pers,nDims);
e_hist2 = zeros(pers,1);
i_k_hist2 = zeros(pers,1);
f_hist2 = zeros(pers,1);

v_dr_hist2 = zeros(pers,1);
v_dt_hist2 = zeros(pers,1);
v_dk_hist2 = zeros(pers,1);
v_hist2 = zeros(pers,1);

hist2(1,:) = [R_0,K_0,T_0,D_0,D_0_base,D_0_tilted];
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
hist2(j,4) = max(min(hist2(j-1,4) + muD(hist2(j-1,:)) * dt + sigmaD(hist2(j-1,:))* shock', upperBounds(:,4)), lowerBounds(:,4)); 
hist2(j,5) = max(min(hist2(j-1,5) + muD_base(hist2(j-1,:)) * dt + sigmaD(hist2(j-1,:))* shock', upperBounds(:,5)), lowerBounds(:,5)); 
hist2(j,6) = max(min(hist2(j-1,6) + muD_tilted(hist2(j-1,:)) * dt + sigmaD(hist2(j-1,:))* shock', upperBounds(:,6)), lowerBounds(:,6)); 

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


%% save results 

save(filename3);

figure
plot(squeeze(e_hists2))
title('E')

figure
plot(squeeze(hists2(:,1,:)))
title('R')

figure
plot(squeeze(hists2(:,2,:)))
title('K')

figure
plot(squeeze(hists2(:,3,:)))
title('F')

figure
plot(squeeze(hists2(:,5,:)))
hold on
plot(squeeze(hists2(:,6,:)))
title('D')

figure
plot(squeeze(f_hists2))
title('F')


