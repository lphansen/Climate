%%%%% This file generates data for Social Cost of Carbon competitive 
%%%%% setting for left panel in Figure G.1

% Authors: Mike Barnett, Jieyao Wang
% Last update: Nov 25, 2019

close all;
clear all;
clc;

file2 = [pwd, '/averse_weighted_SCC_base'];
Model2 = load(file2,'v0','r_mat','k_mat','F_mat','i_k','e','v0_dk','v0_dr','xi_p',...
              'alpha','kappa','delta','expec_e_sum','xi_d','a','b','n','gamma_1','gamma_2','gamma_bar',...
              'bar_gamma_2_plus','beta_f','var_beta_f','power');
external_v0 = Model2.v0;
r_mat = Model2.r_mat;
k_mat = Model2.k_mat;
F_mat = Model2.F_mat;
gamma_1 = Model2.gamma_1;
gamma_2 = Model2.gamma_2;
gamma_bar = Model2.gamma_bar;
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

file1 = [pwd,'/averse_weighted'];
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

file2 = [pwd, '/averse_weighted_SCC_worst'];
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

file1 = [pwd, '/SS_competitive_Sims'];
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

% Generate data for left panel in Figure G.1 in online appendix
s_scc = strcat('SCC_comp_SS_averse_weighted.mat');
save(s_scc,'SCC','SCC1','SCC2','SCC3');

