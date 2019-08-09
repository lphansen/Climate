close all;
clc;
clear all;

warning off all;

directory2 = pwd;
 

%% SCC for original social planner-deterministic:
% 
file2 = [directory2, '/SCC_mat_Cumu_base_GrowthAmb_v2b'];
% file2 = [directory2, '/SCC_mat_Cumu_base_GrowthNoAmb'];
% file2 = [directory2, '/SCC_mat_Cumu_base2'];
% file2 = [directory2, '/SCC_mat_Cumu_worst'];

Model2 = load(file2,'v0','r_mat_1','k_mat_1','d_mat_1','t_mat_1','i_k','f','e','v0_dk','v0_dr','theta','kappa',...
    'A_O','alpha','delta','Theta','Gamma', 'expec_e_sum','xi_d','a','b','n','gamma_1','gamma_2','t_bar',...
    'beta_f','var_beta_f','power');

external_v0 = Model2.v0;
r_mat_1 = Model2.r_mat_1;
k_mat_1 = Model2.k_mat_1;
t_mat_1 = Model2.t_mat_1;
gamma_1 = Model2.gamma_1;
gamma_2 = Model2.gamma_2;
t_bar = Model2.t_bar;
theta = Model2.theta;
kappa = Model2.kappa;
A_O = Model2.A_O;
alpha = Model2.alpha;
delta = Model2.delta;
Gamma = Model2.Gamma;
Theta = Model2.Theta;
a = Model2.a;
b = Model2.b;
beta_f = Model2.beta_f;
var_beta_f = Model2.var_beta_f;
e_1 = Model2.e;
f_1 = Model2.f;
i_k_1 = Model2.i_k;

time_vec = linspace(0,100,400);


%% 
file2 = [pwd, '/HJB_NonLinGrowth_NoAmb'];
% file2 = [pwd, '/../HJB_NonLinGrowth_sig5_v2b'];

Model2 = load(file2,'v0_dr','v0_dt','quadrature');
v0_dr_1 = Model2.v0_dr;
v0_dt_1 = Model2.v0_dt;


%% total
MC = delta.*(1-alpha)./(A_O.*exp(k_mat_1)-i_k_1.*exp(k_mat_1)-f_1.*exp(r_mat_1));

ME = delta.*alpha./(e_1.*exp(r_mat_1));
SCC = 1000*ME./MC;
SCC_func = griddedInterpolant(r_mat_1,t_mat_1,k_mat_1,SCC,'spline');

%% private
ME1 = (v0_dr_1.*exp(-r_mat_1));
SCC1 = 1000*ME1./MC; %private
SCC1_func = griddedInterpolant(r_mat_1,t_mat_1,k_mat_1,SCC1,'linear');

%% Feyman Kac under baseline
ME2_base =  external_v0;
SCC2_base = 1000*ME2_base./MC;
SCC2_base_func = griddedInterpolant(r_mat_1,t_mat_1,k_mat_1,SCC2_base,'spline');

%% comparison for -V_f
ME2_base_a = -v0_dt_1;
SCC2_base_a = 1000*ME2_base_a./MC;
SCC2_base_a_func = griddedInterpolant(r_mat_1,t_mat_1,k_mat_1,SCC2_base_a,'spline');

%%
file2 = [pwd, '/SCC_mat_Cumu_worst_GrowthAmb_v2b'];
% file2 = [pwd, '/SCC_mat_Cumu_worst_GrowthNoAmb'];


Model2 = load(file2,'v0','r_mat_1','k_mat_1','d_mat_1','t_mat_1','i_k','f','e','v0_dk','v0_dr','theta','kappa',...
    'A_O','alpha','delta','Theta','Gamma', 'nd','expec_e_sum','xi_d','a','b','n','gamma_1','gamma_2','t_bar',...
    'beta_f','var_beta_f','power');

external_v0_worst = Model2.v0;

%% Feyman Kac under tilted
ME2_tilt =  external_v0_worst;
SCC2_tilt = 1000*ME2_tilt./MC;
SCC2_tilt_func = griddedInterpolant(r_mat_1,t_mat_1,k_mat_1,SCC2_tilt,'spline');

%% plug in simulation for original one

file1 = [directory2, '/HJB_NonLinGrowth_NoAmb_Sims_determ'];
% file1 = [directory2, '/../HJB_NonLinGrowth_sig5_v2b_Sims_determ'];
Model1 = load(file1,'hists2','e_hists2','theta','kappa','RE_hists2','pi_tilde_1_hists2'...
    ,'pi_tilde_2_hists2','pi_tilde_3_hists2','pi_tilde_4_hists2','pi_tilde_5_hists2'...
    ,'pi_tilde_6_hists2','pi_tilde_7_hists2','pi_tilde_8_hists2','pi_tilde_9_hists2');

hists2_A = Model1.hists2;
RE_hist = Model1.RE_hists2;
e_hist = Model1.e_hists2;
pi_tilde_1_hist = Model1.pi_tilde_1_hists2;
pi_tilde_2_hist = Model1.pi_tilde_2_hists2;
pi_tilde_3_hist = Model1.pi_tilde_3_hists2;
pi_tilde_4_hist = Model1.pi_tilde_4_hists2;
pi_tilde_5_hist = Model1.pi_tilde_5_hists2;
pi_tilde_6_hist = Model1.pi_tilde_6_hists2;
pi_tilde_7_hist = Model1.pi_tilde_7_hists2;
pi_tilde_8_hist = Model1.pi_tilde_8_hists2;
pi_tilde_9_hist = Model1.pi_tilde_9_hists2;
% parfor time=1:400
for time=1:400
    for path=1:1
    SCC_values(time,path) = SCC_func(log((hists2_A(time,1,path))),(hists2_A(time,3,path)),log((hists2_A(time,2,path))));
    SCC1_values(time,path) = SCC1_func(log((hists2_A(time,1,path))),(hists2_A(time,3,path)),log((hists2_A(time,2,path))));
    SCC2_base_values(time,path) = SCC2_base_func(log((hists2_A(time,1,path))),(hists2_A(time,3,path)),log((hists2_A(time,2,path))));
    SCC2_base_a_values(time,path) = SCC2_base_a_func(log((hists2_A(time,1,path))),(hists2_A(time,3,path)),log((hists2_A(time,2,path))));
    SCC2_tilt_values(time,path) = SCC2_tilt_func(log((hists2_A(time,1,path))),(hists2_A(time,3,path)),log((hists2_A(time,2,path))));
    end
end

SCC_total = mean(SCC_values,2);
SCC_private = mean(SCC1_values,2);
SCC2_FK_base = mean(SCC2_base_values,2);
SCC2_V_f = mean(SCC2_base_a_values,2);
SCC2_FK_tilt = mean(SCC2_tilt_values,2);


figure('pos',[10,10,600,500]);
plot(time_vec,mean(RE_hist,2),'LineWidth',2)
% legend('Relative Entropy','location','northwest')
hold off 
title(sprintf('Relative Entropy, $\\xi_a =$ %0.5f',1/theta),'Interpreter','latex')
xlabel('Years')
% ylabel('$/Ton')
axis([0 100 0 1])
set(findall(gcf,'type','axes'),'fontsize',12,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',12,'fontWeight','bold')
print('Simu-Un-compare','-dpng')


figure('pos',[10,10,600,500]);
hold on
plot(time_vec,mean(pi_tilde_1_hist,2),'LineWidth',2)
plot(time_vec,mean(pi_tilde_2_hist,2),'LineWidth',2)
plot(time_vec,mean(pi_tilde_3_hist,2),'LineWidth',2)
plot(time_vec,mean(pi_tilde_4_hist,2),'LineWidth',2)
plot(time_vec,mean(pi_tilde_5_hist,2),'LineWidth',2)
plot(time_vec,mean(pi_tilde_6_hist,2),'LineWidth',2)
plot(time_vec,mean(pi_tilde_7_hist,2),'LineWidth',2)
plot(time_vec,mean(pi_tilde_8_hist,2),'LineWidth',2)
plot(time_vec,mean(pi_tilde_9_hist,2),'LineWidth',2)
legend('Distorted Weight 1', 'Distorted Weight 2',...
    'Distorted Weight 3', 'Distorted Weight 4','Distorted Weight 5', 'Distorted Weight 6',...
    'Distorted Weight 7', 'Distorted Weight 8', 'Distorted Weight 9','location','northwest')
hold off 
title(sprintf('Distorted Weights, $\\xi_a =$ %0.5f',1/theta),'Interpreter','latex')
xlabel('Years')
% ylabel('$/Ton')
axis([0 100 0 1])
set(findall(gcf,'type','axes'),'fontsize',12,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',12,'fontWeight','bold')
print('Simu-Un-compare','-dpng')


figure('pos',[10,10,600,500]);
hold on
plot(time_vec,mean(SCC2_FK_tilt,2),'LineWidth',2)
plot(time_vec,mean(SCC2_V_f,2),'LineWidth',2)
plot(time_vec,mean(SCC2_FK_base,2),'LineWidth',2)
legend('Tilted F-K','Actual -V_f','Baseline F-K','location','northwest')
hold off 
title(sprintf('SCC -V_F Component, $\\xi_a =$ %0.2f',1/theta),'Interpreter','latex')
xlabel('Years')
ylabel('$/Ton')
set(findall(gcf,'type','axes'),'fontsize',12,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',12,'fontWeight','bold')
print('Simu-Un-compare','-dpng')


figure('pos',[10,10,600,500]);
hold on
plot(time_vec,mean(SCC2_FK_tilt,2),'LineWidth',2)
plot(time_vec,mean(SCC2_FK_base,2),'LineWidth',2)
plot(time_vec,mean(SCC2_FK_tilt-SCC2_FK_base,2),'LineWidth',2)
legend('Tilt','Base','Difference','location','northwest')
hold off 
title(sprintf('SCC Tilted and Baseline Comparison, $\\xi_a =$ %0.2f',1/theta),'Interpreter','latex')
xlabel('Years')
ylabel('$/Ton')
set(findall(gcf,'type','axes'),'fontsize',12,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',12,'fontWeight','bold')
print('Simu-Un-compare','-dpng')

figure('pos',[10,10,600,500]);
plot(time_vec,mean(SCC_total,2),'LineWidth',2)
hold on
% plot(time_vec,mean(SCC_private,2),'LineWidth',2)
% % plot(time_vec,mean(SCC2_V_f,2),'LineWidth',2)
% plot(time_vec,mean(SCC2_FK_base,2),'LineWidth',2)
% plot(time_vec,mean(SCC2_FK_tilt-SCC2_FK_base,2),'LineWidth',2);
% plot(time_vec,mean(SCC2_FK_tilt,2),'LineWidth',2);
plot(time_vec,mean(SCC2_V_f,2),'LineWidth',2);
legend('Total SCC','Private','Social','Uncertainty','location','northwest')
hold off 
title(sprintf('SCC Decomposition, $\\xi_a =$ %0.5f',1/theta),'Interpreter','latex')
xlabel('Years')
ylabel('$/Ton')
set(findall(gcf,'type','axes'),'fontsize',12,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',12,'fontWeight','bold')
print('Simu-Un-compare','-dpng')

% external = mean(SCC2_FK_base,2);
% uncert = mean(SCC2_FK_tilt-SCC2_FK_base,2);
% external = mean(SCC2_FK_base,2);
external = mean(SCC2_V_f,2);
uncert = mean(SCC2_FK_tilt-SCC2_FK_base,2);

external(1)
external(end./2)
external(end)

uncert(1)
uncert(end./2)
uncert(end)

e_hist(1,1,1)
e_hist(end./2,1,1)
e_hist(end,1,1)

RE_hist(1)
RE_hist(end./2)
RE_hist(end)

figure('pos',[10,10,600,500]);
plot(time_vec,hists2_A(:,1,1),'LineWidth',2)
title(sprintf('Reserves Time Series Path, $\\xi_a =$ %0.2f',1/theta),'Interpreter','latex')
xlabel('Years')
ylabel('GtC')
set(findall(gcf,'type','axes'),'fontsize',12,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',12,'fontWeight','bold')
print('Simu-Un-compare','-dpng')

figure('pos',[10,10,600,500]);
plot(time_vec,hists2_A(:,3,1),'LineWidth',2)
title(sprintf('Concentrations Time Series Path, $\\xi_a =$ %0.2f',1/theta),'Interpreter','latex')
xlabel('Years')
ylabel('GtC')
set(findall(gcf,'type','axes'),'fontsize',12,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',12,'fontWeight','bold')
print('Simu-Un-compare','-dpng')

figure('pos',[10,10,600,500]);
plot(time_vec,e_hist(:,1,1),'LineWidth',2)
title(sprintf('Emissions Time Series Path, $\\xi_a =$ %0.5f',1/theta),'Interpreter','latex')
xlabel('Years')
ylabel('GtC')
set(findall(gcf,'type','axes'),'fontsize',12,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',12,'fontWeight','bold')
print('Simu-Un-compare','-dpng')
