% close all;
clc;
clear all;

warning off all;

directory2 = pwd;


%% SCC for original

file2 = [directory2, '/SCC_mat_Cumu_base'];
Model2 = load(file2,'v0','r_mat_1','k_mat_1','d_mat_1','t_mat_1','i_k','f','e','v0_dk','v0_dr','theta','kappa',...
    'A_O','alpha','delta','Theta','Gamma', 'nd','expec_e_sum','xi_d','a','b','n','gamma1','gamma2','f_bar',...
    'lambda','beta_f','var_beta_f','power');

external_v0 = Model2.v0;
r_mat_1 = Model2.r_mat_1;
k_mat_1 = Model2.k_mat_1;
d_mat_1 = Model2.d_mat_1;
t_mat_1 = Model2.t_mat_1;
% f_1 = Model2.f;
% e_1 = Model2.e;
theta = Model2.theta;
kappa = Model2.kappa;
A_O = Model2.A_O;
alpha = Model2.alpha;
delta = Model2.delta;
lambda = Model2.lambda;
gamma1 = Model2.gamma1;
gamma2 = Model2.gamma2;
Gamma = Model2.Gamma;
Theta = Model2.Theta;
nd = Model2.nd;
% expec_e_sum_1 = Model2.expec_e_sum;
xi_d = Model2.xi_d;
a = Model2.a;
b = Model2.b;
n = Model2.n;
beta_f = Model2.beta_f;
var_beta_f = Model2.var_beta_f;

time_vec = linspace(0,100,400);

%% recompute the plots:

%% 
file1 = [pwd,'/HJB_NonLinPref_Cumu_check_exp'];
Model1 = load(file1,'v0_dr','v0_dk','v0_dt','v0_dtt','i_k','f','e','expec_e_sum');
v0_dk_1 = repmat(Model1.v0_dk,[1,1,1,nd]);
v0_dr_1 = repmat(Model1.v0_dr,[1,1,1,nd]);
v0_dt_1 = repmat(Model1.v0_dt,[1,1,1,nd]);
v0_dtt_1 = repmat(Model1.v0_dtt,[1,1,1,nd]);
i_k_1 = repmat(Model1.i_k,[1,1,1,nd]);
f_1 = repmat(Model1.f,[1,1,1,nd]);
e_1 = repmat(Model1.e,[1,1,1,nd]);
expec_e_sum_1 = repmat(Model1.expec_e_sum,[1,1,1,nd]);


MC = delta.*(1-alpha)./(A_O.*exp(k_mat_1)-i_k_1.*exp(k_mat_1)-f_1.*exp(r_mat_1));

%% 
ME1 = (v0_dr_1.*exp(-r_mat_1));  
%% total
ME = delta.*alpha./(e_1.*exp(r_mat_1));
SCC = 1000*ME./MC;
SCC_func = griddedInterpolant(r_mat_1,t_mat_1,k_mat_1,d_mat_1,SCC,'spline');

%% private
SCC1 = 1000*ME1./MC; %private
SCC1_func = griddedInterpolant(r_mat_1,t_mat_1,k_mat_1,d_mat_1,SCC1,'linear');

%% Feyman Kac under baseline
ME2_base =  (1-alpha).*external_v0;
SCC2_base = 1000*ME2_base./MC;
SCC2_base_func = griddedInterpolant(r_mat_1,t_mat_1,k_mat_1,d_mat_1,SCC2_base,'spline');


%% comparison for -V_f
ME2_base_a = -v0_dt_1;
SCC2_base_a = 1000*ME2_base_a./MC;
SCC2_base_a_func = griddedInterpolant(r_mat_1,t_mat_1,k_mat_1,d_mat_1,SCC2_base_a,'spline');

%% V_d term under baseline

mean_base = beta_f;
lambda_tilde_base = lambda;
base_model_drift = (mean(gamma1).*mean_base...
    +mean(gamma2).*t_mat_1.*(1./lambda_tilde_base+mean_base.^2));
    
    
% V_d_baseline_func = @(x) xi_d...
%          .*(gamma_1.*x +gamma_2.*t_mat_1.*x.^2 ...
%         +bar_gamma_2_plus.*x.*(x.*t_mat_1-f_bar).^(power-1).*((x.*t_mat_1-f_bar)>=0)) ...
%         .*normpdf(x,beta_f,sqrt(var_beta_f));
V_d_baseline = xi_d * base_model_drift;



ME2b = - V_d_baseline;
SCC2_V_d_baseline = 1000*ME2b./MC;
SCC2_V_d_baseline_func = griddedInterpolant(r_mat_1,t_mat_1,k_mat_1,d_mat_1,SCC2_V_d_baseline,'spline');

%%
file2 = [pwd, '/SCC_mat_Cumu_worst'];

Model2 = load(file2,'v0','r_mat_1','k_mat_1','d_mat_1','t_mat_1','i_k','f','e','v0_dk','v0_dr','theta','kappa',...
    'A_O','alpha','delta','v0_dt','expec_e_sum',...
    'xi_d','gamma_1','gamma_2','gamma_2_plus','f_bar','beta_f','var_beta_f',...
    'power','lambda','weight','lambda_tilde_1','beta_tilde_1','a','b','n');

external_v0_worst = Model2.v0;
r_mat_1 = Model2.r_mat_1;
k_mat_1 = Model2.k_mat_1;
d_mat_1 = Model2.d_mat_1;
t_mat_1 = Model2.t_mat_1;


%% Feyman Kac under tilted
ME2_tilt =  (1-alpha).*external_v0_worst;
SCC2_tilt = 1000*ME2_tilt./MC;
SCC2_tilt_func = griddedInterpolant(r_mat_1,t_mat_1,k_mat_1,d_mat_1,SCC2_tilt,'spline');
%% V_d term under tilted

ME2b = - expec_e_sum_1.*exp(-r_mat_1);

SCC2_V_d_tilt = 1000*ME2b./MC;
SCC2_V_d_tilt_func = griddedInterpolant(r_mat_1,t_mat_1,k_mat_1,d_mat_1,SCC2_V_d_tilt,'spline');

%% plug in simulation for original one

file1 = [directory2, '/HJB_NonLinPref_Cumu_check_exp_Sims'];
Model1 = load(file1,'hists2','e_hists2','theta','kappa');

hists2_A = Model1.hists2;
e_hist = Model1.e_hists2;
% parfor time=1:400
for time=1:400
    for path=1:1
    SCC_values(time,path) = ...
        SCC_func(log((hists2_A(time,1,path))),(hists2_A(time,3,path)),log((hists2_A(time,2,path))),...
        (hists2_A(time,6,path)));
    SCC1_values(time,path) = SCC1_func(log((hists2_A(time,1,path))),(hists2_A(time,3,path)),log((hists2_A(time,2,path))),(hists2_A(time,5,path)));
    SCC2_base_values(time,path) = SCC2_base_func(log((hists2_A(time,1,path))),(hists2_A(time,3,path)),log((hists2_A(time,2,path))),(hists2_A(time,5,path)));
    SCC2_base_a_values(time,path) = SCC2_base_a_func(log((hists2_A(time,1,path))),(hists2_A(time,3,path)),log((hists2_A(time,2,path))),(hists2_A(time,5,path)));
    SCC2_tilt_values(time,path) = SCC2_tilt_func(log((hists2_A(time,1,path))),(hists2_A(time,3,path)),log((hists2_A(time,2,path))),(hists2_A(time,6,path)));
    SCC2_V_d_baseline_values(time,path) = SCC2_V_d_baseline_func(log((hists2_A(time,1,path))),(hists2_A(time,3,path)),log((hists2_A(time,2,path))),(hists2_A(time,5,path)));
    SCC2_V_d_tilt_values(time,path) = SCC2_V_d_tilt_func(log((hists2_A(time,1,path))),(hists2_A(time,3,path)),log((hists2_A(time,2,path))),(hists2_A(time,6,path)));
    
    end
end

SCC_total = mean(SCC_values,2);
SCC_private = mean(SCC1_values,2);
SCC2_FK_base = mean(SCC2_base_values,2);
SCC2_V_f = mean(SCC2_base_a_values,2);
SCC2_FK_tilt = mean(SCC2_tilt_values,2);
SCC2_V_d_baseline = mean(SCC2_V_d_baseline_values,2);
SCC2_V_d_tilt = mean(SCC2_V_d_tilt_values,2);

SCC = SCC_total;
SCC1 = SCC_private;
SCC2 = SCC2_FK_base+SCC2_V_d_baseline;
SCC3 = SCC2_V_d_tilt-SCC2_V_d_baseline...
    +SCC2_FK_tilt-SCC2_FK_base;


s1 = num2str(1./theta,4);
s1 = strrep(s1,'.','');
s2 = num2str(kappa,4);
s2 = strrep(s2,'.','');
% s_scc = strcat('SCC_A=',s1,'_M=',s2,'.mat');
s_scc = strcat('SCC_',s1,'.mat');
save(s_scc,'SCC','SCC1','SCC2','SCC3');


figure('pos',[10,10,600,500]);
plot(time_vec,SCC,'LineWidth',2)
hold on
plot(time_vec,SCC1,'LineWidth',2)
plot(time_vec,SCC2,'LineWidth',2)
plot(time_vec,SCC3,'LineWidth',2);

legend('Total SCC','Private','Social','Uncertainty')
hold off 
title('SCC Decomposition, $\xi_a=0.0002$','Interpreter','latex')
xlabel('Years')
ylabel('$/Ton')
set(findall(gcf,'type','axes'),'fontsize',12,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',12,'fontWeight','bold')
print('Simu-Un-compare','-dpng')



% figure('pos',[10,10,600,500]);
% plot(time_vec,SCC2_FK_base,'LineWidth',2)
% hold on
% plot(time_vec,SCC2_FK_tilt,'LineWidth',2)
% legend('Base','Tilt')
% hold off 
% title('SCC FK Pieces','Interpreter','latex')
% xlabel('Years')
% ylabel('$/Ton')
% set(findall(gcf,'type','axes'),'fontsize',12,'fontWeight','bold')
% set(findall(gcf,'type','text'),'FontSize',12,'fontWeight','bold')
% print('Simu-Un-compare','-dpng')

figure('pos',[10,10,600,500]);
plot(time_vec,e_hist(:,1,1),'LineWidth',2)
title(sprintf('Emissions Time Series Path, $\\xi_a =$ %0.5f',1/theta),'Interpreter','latex')
xlabel('Years')
ylabel('GtC')
set(findall(gcf,'type','axes'),'fontsize',12,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',12,'fontWeight','bold')
print('Simu-Un-compare','-dpng')

