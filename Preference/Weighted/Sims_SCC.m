% load(XXX);
%%
T = 100; % 100 years
pers = 4*T; % quarterly
dt = T/pers;
nDims = 5;
its = 1;

efunc = griddedInterpolant(r_mat,F_mat,k_mat,e,'linear');
% if weight == 1
jfunc = griddedInterpolant(r_mat,F_mat,k_mat,j.^psi_1,'linear');
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
j_func = @(x) exp(log(max(jfunc(log(x(:,1)),x(:,3),log(x(:,2))),0))/psi_1);
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

% % save results
% save('HJB_NonLinPref_Cumu.mat');
% e_values = mean(e_hists2,2);
% s1 = num2str(1.*xi_p,4);
% s1 = strrep(s1,'.','');
% s_e = strcat('e_A=',s1,'.mat');
% save(s_e,'e_values');

%% Step 4: SCC calculation
% close all
% clear all
% clc

% Feyman Kac
% load([pwd, '/HJB_NonLinPref_Cumu']);
j = j_temp;
base_model_flow_func = @(x) ...
    (gamma_2.*x.^2 +bar_gamma_2_plus.*x.^2.*((x.*F_mat-gamma_bar)>=0)).*exp(r_mat).*e...
    .*normpdf(x,beta_f,sqrt(var_beta_f)); 
base_model_flow = quad_int(base_model_flow_func, a, b, n,'legendre');
flow_base = base_model_flow;

v0SCC = (kappa).*r_mat+(1-kappa).*k_mat-beta_f.*F_mat; % initial guess
vold = ones(size(v0SCC));
converged = 0;
iter = 1;
while converged == 0

v0SCC_dr = finiteDifference(v0SCC,1,1,hr,'central', 1e-8);
v0SCC_dt = finiteDifference(v0SCC,2,1,ht,'central');
v0SCC_dk = finiteDifference(v0SCC,3,1,ht,'central');

v0SCC_drr = finiteDifference(v0SCC,1,2,hr,'central');
v0SCC_drr(v0SCC_dr <= 1e-8) = 0;
v0SCC_dtt = finiteDifference(v0SCC,2,2,ht,'central');
v0SCC_dkk = finiteDifference(v0SCC,3,2,ht,'central');

% inputs for solver
A = -delta.*ones(size(r_mat));
B_r = -e+psi_0.*(j.^psi_1).*exp(psi_1.*(k_mat-r_mat))-0.5.*(sigma_r.^2);
B_k = mu_k+phi_0.*log(1+i_k.*phi_1)-0.5.*(sigma_k.^2);
B_t = e.*exp(r_mat);
C_rr = 0.5.*sigma_r.^2.*ones(size(r_mat));
C_kk = 0.5.*sigma_k.^2.*ones(size(r_mat));
C_tt = zeros(size(r_mat));

D = flow_base;
            
stateSpace = [r_mat(:), F_mat(:), k_mat(:)]; 
model      = {};
model.A    = A(:); 
model.B    = [B_r(:), B_t(:), B_k(:)];
model.C    = [C_rr(:), C_tt(:), C_kk(:)];
model.D    = D(:);
model.v0   = v0SCC(:);
model.dt   = 0.5;
out = solveCGNatural(stateSpace, model);

% save results
external_v0 = reshape(out,size(r_mat));
SCC_err = max(max(max(abs(external_v0 - v0SCC))));
disp(['Diff: ', num2str(SCC_err),' Number of Iters: ',num2str(iter)])

if SCC_err < 1e-8
    converged = 1;
else
    v0SCC = external_v0;
    iter = iter + 1;
end



end
% filename2 = [pwd, '/SCC_mat_Cumu_base'];
% save(filename2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear all
% close all
% clc
%% 
% load([pwd, '/HJB_NonLinPref_Cumu']);

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

v0SCC = (kappa).*r_mat+(1-kappa).*k_mat-beta_f.*F_mat; % initial guess
vold = ones(size(v0SCC));
converged = 0;
iter = 1;

while converged == 0

v0SCC_dr = finiteDifference(v0SCC,1,1,hr,'central', 1e-8);
v0SCC_dt = finiteDifference(v0SCC,2,1,ht,'central');
v0SCC_dk = finiteDifference(v0SCC,3,1,ht,'central');

v0SCC_drr = finiteDifference(v0SCC,1,2,hr,'central');
v0SCC_drr(v0SCC_dr <= 1e-8) = 0;
v0SCC_dtt = finiteDifference(v0SCC,2,2,ht,'central');
v0SCC_dkk = finiteDifference(v0SCC,3,2,ht,'central');
            
stateSpace = [r_mat(:), F_mat(:), k_mat(:)]; 
model      = {};
model.A    = A(:); 
model.B    = [B_r(:), B_t(:), B_k(:)];
model.C    = [C_rr(:), C_tt(:), C_kk(:)];
model.D    = D(:);
model.v0   = v0SCC(:);
model.dt   = 0.5;
out = solveCGNatural(stateSpace, model);

% save results

external_v0_worst = reshape(out,size(r_mat));
SCC_err = max(max(max(abs(external_v0_worst - v0SCC))));
disp(['Diff: ', num2str(SCC_err),' Number of Iters: ',num2str(iter)])

if SCC_err < 1e-8
    converged = 1;
else
    v0SCC = external_v0_worst;
    iter = iter + 1;
end



end

% stateSpace = [r_mat(:), F_mat(:), k_mat(:)]; 
% model      = {};
% model.A    = A(:); 
% model.B    = [B_r(:), B_t(:), B_k(:)];
% model.C    = [C_rr(:), C_tt(:), C_kk(:)];
% model.D    = D(:);
% model.v0   = v0(:).*0;
% model.dt   = 1.0;
% 
% out = solveCGNatural_1(stateSpace, model);
% 
% external_v0_worst = reshape(out,size(r_mat));
% filename2 = [pwd, '/SCC_mat_Cumu_worst'];
% save(filename2)

%%% SCC calculation
% close all;
% clear all;
% clc;

% load results
% file2 = [pwd, '/SCC_mat_Cumu_base'];
% Model2 = load(file2,'v0','r_mat','k_mat','F_mat','i_k','j','e','v0_dk','v0_dr','xi_p',...
%     'alpha','kappa','delta','expec_e_sum','xi_d','a','b','n','gamma_1','gamma_2','gamma_bar',...
%     'bar_gamma_2_plus','beta_f','var_beta_f','power');
% external_v0 = Model2.v0;
% r_mat = Model2.r_mat;
% k_mat = Model2.k_mat;
% F_mat = Model2.F_mat;
% gamma_1 = Model2.gamma_1;
% gamma_2 = Model2.gamma_2;
% gamma_bar = Model2.gamma_bar;
% xi_p = Model2.xi_p;
% alpha = Model2.alpha;
% kappa = Model2.kappa;
% delta = Model2.delta;
% xi_d = Model2.xi_d;
% a = Model2.a;
% b = Model2.b;
% n = Model2.n;
% bar_gamma_2_plus = Model2.bar_gamma_2_plus;
% beta_f = Model2.beta_f;
% var_beta_f = Model2.var_beta_f;
% power = Model2.power;

% file1 = [pwd,'/HJB_NonLinPref_Cumu'];
% Model1 = load(file1,'v0_dr','v0_dk','v0_dt','i_k','j','e','expec_e_sum');
% v0_dk = Model1.v0_dk;
% v0_dr = Model1.v0_dr;
% v0_dt = Model1.v0_dt;
% i_k = Model1.i_k;
% j = Model1.j;
% e = Model1.e;
% expec_e_sum = Model1.expec_e_sum;
%% 
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

% file2 = [pwd, '/SCC_mat_Cumu_worst'];
% Model2 = load(file2,'v0','r_mat','k_mat','F_mat');
% external_v0_worst = Model2.v0;
% r_mat = Model2.r_mat;
% k_mat = Model2.k_mat;
% F_mat = Model2.F_mat;

ME2_tilt =  (1-kappa).*external_v0_worst;
SCC2_tilt = 1000*ME2_tilt./MC;
SCC2_tilt_func = griddedInterpolant(r_mat,F_mat,k_mat,SCC2_tilt,'linear');

ME2b = - expec_e_sum.*exp(-r_mat);
SCC2_V_d_tilt = 1000*ME2b./MC;
SCC2_V_d_tilt_func = griddedInterpolant(r_mat,F_mat,k_mat,SCC2_V_d_tilt,'linear');

% file1 = [pwd, '/HJB_NonLinPref_Cumu_Sims'];
% Model1 = load(file1,'hists2');
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

save(xxx);