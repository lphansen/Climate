%%%%% This file generates Social Cost of Carbon decomposition for the Consumption Damage model.
% Authors: Mike Barnett, Jieyao Wang
% Last update: Dec 9, 2019
close all
clear all
clc

%% Step 0: Set up solver
mex solveCGNatura_l.cpp;

%% Step 1: Specify ambiguity and damage level

% Ambiguity setting: 'averse' or 'neutral'
ambiguity = 'averse';
% Damage setting: 'high', 'low', or 'weighted'
damage_level = 'weighted';

if strcmp(ambiguity,'averse')
    xi_p = 1 ./ 4000; % ambiguity parameter
elseif strcmp(ambiguity,'neutral')
    xi_p = 1000; % ambiguity parameter
else
    disp('Error: please choose between ''averse'' and ''weighted'' for ambiguity level');
end

if strcmp(damage_level,'high')
    weight = 0.0; % weight on nordhaus
elseif strcmp(damage_level,'low')
    weight = 1.0; % weight on nordhaus
elseif strcmp(damage_level,'weighted')
    weight = 0.5; % weight on nordhaus
else
    disp('Error: please choose from ''high'',''weighted'' and ''low'' for damage level')
end

%% Step 2: Solve FK equation
% load HJB solution
load([ambiguity,'_',damage_level]);

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
            
stateSpace = [r_mat(:), F_mat(:), k_mat(:)]; 
model      = {};
model.A    = A(:); 
model.B    = [B_r(:), B_t(:), B_k(:)];
model.C    = [C_rr(:), C_tt(:), C_kk(:)];
model.D    = D(:);
model.v0   = v0(:).*0;
model.dt   = 1.0;
    
out = solveCGNatural_1(stateSpace, model);

% save results
v0 = reshape(out,size(r_mat));
save([ambiguity,'_',damage_level,'_SCC_base']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load HJB solution
load([ambiguity,'_',damage_level]);


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

stateSpace = [r_mat(:), F_mat(:), k_mat(:)]; 
model      = {};
model.A    = A(:); 
model.B    = [B_r(:), B_t(:), B_k(:)];
model.C    = [C_rr(:), C_tt(:), C_kk(:)];
model.D    = D(:);
model.v0   = v0(:).*0;
model.dt   = 1.0;

out = solveCGNatural_1(stateSpace, model);

v0 = reshape(out,size(r_mat));
save([ambiguity,'_',damage_level,'_SCC_worst']);


%% Step 3: Calculate SCC

% load results
file2 = [ambiguity,'_',damage_level,'_SCC_base'];
Model2 = load(file2,'v0','r_mat','k_mat','F_mat','i_k','j','e','v0_dk','v0_dr','xi_p',...
    'alpha','kappa','delta','expec_e_sum','xi_d','a','b','n','gamma_1','gamma_2','gamma_bar',...
    'bar_gamma_2_plus','beta_f','var_beta_f','power');
external_v0 = Model2.v0;
r_mat = Model2.r_mat;
k_mat = Model2.k_mat;
F_mat = Model2.F_mat;
gamma_1 = Model2.gamma_1;
gamma_2 = Model2.gamma_2;
gamma_bar = Model2.gamma_bar;
xi_p = Model2.xi_p;
alpha = Model2.alpha;
kappa = Model2.kappa;
delta = Model2.delta;
xi_d = Model2.xi_d;
a = Model2.a;
b = Model2.b;
n = Model2.n;
bar_gamma_2_plus = Model2.bar_gamma_2_plus;
beta_f = Model2.beta_f;
var_beta_f = Model2.var_beta_f;
power = Model2.power;

file1 = [ambiguity,'_',damage_level];
Model1 = load(file1,'v0_dr','v0_dk','v0_dt','i_k','j','e','expec_e_sum');
v0_dk = Model1.v0_dk;
v0_dr = Model1.v0_dr;
v0_dt = Model1.v0_dt;
i_k = Model1.i_k;
j = Model1.j;
e = Model1.e;
expec_e_sum = Model1.expec_e_sum;

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

file2 = [ambiguity,'_',damage_level,'_SCC_worst'];
Model2 = load(file2,'v0','r_mat','k_mat','F_mat');
external_v0_worst = Model2.v0;
r_mat = Model2.r_mat;
k_mat = Model2.k_mat;
F_mat = Model2.F_mat;

ME2_tilt =  (1-kappa).*external_v0_worst;
SCC2_tilt = 1000*ME2_tilt./MC;
SCC2_tilt_func = griddedInterpolant(r_mat,F_mat,k_mat,SCC2_tilt,'linear');

ME2b = - expec_e_sum.*exp(-r_mat);
SCC2_V_d_tilt = 1000*ME2b./MC;
SCC2_V_d_tilt_func = griddedInterpolant(r_mat,F_mat,k_mat,SCC2_V_d_tilt,'linear');

file1 = [ambiguity,'_',damage_level,'_sim'];
Model1 = load(file1,'hists2');
hists2_A = Model1.hists2;

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
s_scc = [ambiguity,'_',damage_level,'_SCC_A=',s1,'.mat'];
save(s_scc,'SCC','SCC1','SCC2','SCC3');
