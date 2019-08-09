clear all
close all
clc


%% load in data
file1 = [pwd, '/HJB_NonLinPref_Cumu'];
Model1 = load(file1,'R_1','pi_tilde_1_norm','r_mat','k_mat','t_mat',...
    'beta_tilde_1','beta_f','lambda_tilde_1','theta','var_beta_f'); 


RE_1 = Model1.R_1;
pi_tilde_1 = Model1.pi_tilde_1_norm;
beta_tilde_1 = Model1.beta_tilde_1;
beta_f = Model1.beta_f;
lambda_tilde_1 = Model1.lambda_tilde_1;
theta = Model1.theta;
var_beta_f = Model1.var_beta_f;

r_mat = Model1.r_mat;
t_mat = Model1.t_mat;
k_mat = Model1.k_mat;

a = beta_f-10.*sqrt(var_beta_f);
b = beta_f+10.*sqrt(var_beta_f);

RE_1_func = griddedInterpolant(r_mat,t_mat,k_mat,RE_1,'spline');
pi_tilde_1_func = griddedInterpolant(r_mat,t_mat,k_mat,pi_tilde_1,'spline');
beta_tilde_1_func = griddedInterpolant(r_mat,t_mat,k_mat,beta_tilde_1,'spline');
lambda_tilde_1_func = griddedInterpolant(r_mat,t_mat,k_mat,lambda_tilde_1,'spline');

file1 = [pwd, '/HJB_NonLinPref_Cumu_Sims'];
Model1 = load(file1,'hists2'); 
hists2 = Model1.hists2;
T_value = mean(squeeze(hists2(:,3,:)),2);
R_value = mean(squeeze(hists2(:,1,:)),2);
K_value = mean(squeeze(hists2(:,2,:)),2);

time_vec = linspace(0,100,400);

%% Generate Relative Entropy and Weights

for time=1:400
    RE_1_plot(time) = RE_1_func(log(R_value(time,1)),T_value(time,1),log(K_value(time,1)));
    weight_plot(time) = pi_tilde_1_func(log(R_value(time,1)),T_value(time,1),log(K_value(time,1)));
    nordhaus_mean(time) = beta_tilde_1_func(log(R_value(time,1)),T_value(time,1),log(K_value(time,1)));
    sd_nordhaus(time) = 1./sqrt(lambda_tilde_1_func(log(R_value(time,1)),T_value(time,1),log(K_value(time,1))));
end

fileID = fopen('Relative Entropy.txt','w');
fprintf(fileID,'xi_a: %.6f \n',1./theta);
fprintf(fileID,'0 yr: %.6f \n',RE_1_plot(1));
fprintf(fileID,'25 yr: %.6f \n',RE_1_plot(100));
fprintf(fileID,'50 yr: %.6f \n',RE_1_plot(200));
fprintf(fileID,'75 yr: %.6f \n',RE_1_plot(300));
fprintf(fileID,'100 yr: %.6f \n',RE_1_plot(400));
fclose(fileID);

fileID = fopen('Shifted Mean.txt','w');
fprintf(fileID,'xi_a: %.6f \n',1./theta);
fprintf(fileID,'0 yr: %.6f \n',nordhaus_mean(1));
fprintf(fileID,'25 yr: %.6f \n',nordhaus_mean(100));
fprintf(fileID,'50 yr: %.6f \n',nordhaus_mean(200));
fprintf(fileID,'75 yr: %.6f \n',nordhaus_mean(300));
fprintf(fileID,'100 yr: %.6f \n',nordhaus_mean(400));
fclose(fileID);

fileID = fopen('Shifted Standard Deviation.txt','w');
fprintf(fileID,'xi_a: %.6f \n',1./theta);
fprintf(fileID,'0 yr: %.6f \n',sd_nordhaus(1));
fprintf(fileID,'25 yr: %.6f \n',sd_nordhaus(100));
fprintf(fileID,'50 yr: %.6f \n',sd_nordhaus(200));
fprintf(fileID,'75 yr: %.6f \n',sd_nordhaus(300));
fprintf(fileID,'100 yr: %.6f \n',sd_nordhaus(400));
fclose(fileID);

fileID = fopen('Nordhaus Weight.txt','w');
fprintf(fileID,'xi_a: %.6f \n',1./theta);
fprintf(fileID,'0 yr: %.6f \n',weight_plot(1));
fprintf(fileID,'25 yr: %.6f \n',weight_plot(100));
fprintf(fileID,'50 yr: %.6f \n',weight_plot(200));
fprintf(fileID,'75 yr: %.6f \n',weight_plot(300));
fprintf(fileID,'100 yr: %.6f \n',weight_plot(400));
fclose(fileID);

figure('pos',[10,10,800,500]);
yyaxis left
plot(time_vec,RE_1_plot,'-','LineWidth',2.5);
ylim([0 1])
hold on
yyaxis right
ylim([0 1])
plot(time_vec,weight_plot,'-','LineWidth',2.5);
title('Relative Entropy and Weights','Interpreter','latex')
legend('Relative Entropy','normalized weight on Nordhaus','Original','weighted')
xlabel('Year','Interpreter','latex')
set(findall(gcf,'type','axes'),'fontsize',16,'fontWeight','bold')
set(findall(gcf,'type','text'),'FontSize',16,'fontWeight','bold')
print('RE_weights','-dpng')

beta_f_space = linspace(a,b,200);
save('beta_f_space','beta_f_space')

%% T=0
time=1;

R0 = R_value(time,1);
K0 = K_value(time,1);
F0 = T_value(time,1);

mean_distort_nordhaus = beta_tilde_1_func(log(R0),F0,log(K0))-beta_f;
lambda_tilde_nordhaus = lambda_tilde_1_func(log(R0),F0,log(K0));

original_dist = normpdf(beta_f_space,beta_f,sqrt(var_beta_f));
nordhaus_dist_0 = normpdf(beta_f_space,mean_distort_nordhaus+beta_f,1./sqrt(lambda_tilde_nordhaus));

nordhaus = nordhaus_dist_0;
original = original_dist;

save('Dist_0yr','nordhaus','original')


%% T=25
time=100;

R0 = R_value(time,1);
K0 = K_value(time,1);
F0 = T_value(time,1);

mean_distort_nordhaus = beta_tilde_1_func(log(R0),F0,log(K0))-beta_f;
lambda_tilde_nordhaus = lambda_tilde_1_func(log(R0),F0,log(K0));

original_dist = normpdf(beta_f_space,beta_f,sqrt(var_beta_f));
nordhaus_dist_0 = normpdf(beta_f_space,mean_distort_nordhaus+beta_f,1./sqrt(lambda_tilde_nordhaus));

nordhaus = nordhaus_dist_0;
original = original_dist;

save('Dist_25yr','nordhaus','original')


%% T=50
time=200;

R0 = R_value(time,1);
K0 = K_value(time,1);
F0 = T_value(time,1);

mean_distort_nordhaus = beta_tilde_1_func(log(R0),F0,log(K0))-beta_f;
lambda_tilde_nordhaus = lambda_tilde_1_func(log(R0),F0,log(K0));

original_dist = normpdf(beta_f_space,beta_f,sqrt(var_beta_f));
nordhaus_dist_0 = normpdf(beta_f_space,mean_distort_nordhaus+beta_f,1./sqrt(lambda_tilde_nordhaus));

nordhaus = nordhaus_dist_0;
original = original_dist;

save('Dist_50yr','nordhaus','original')


%% T=75
time=300;

R0 = R_value(time,1);
K0 = K_value(time,1);
F0 = T_value(time,1);

mean_distort_nordhaus = beta_tilde_1_func(log(R0),F0,log(K0))-beta_f;
lambda_tilde_nordhaus = lambda_tilde_1_func(log(R0),F0,log(K0));

original_dist = normpdf(beta_f_space,beta_f,sqrt(var_beta_f));
nordhaus_dist_0 = normpdf(beta_f_space,mean_distort_nordhaus+beta_f,1./sqrt(lambda_tilde_nordhaus));

nordhaus = nordhaus_dist_0;
original = original_dist;

save('Dist_75yr','nordhaus','original')

%% T=100
time=400;

R0 = R_value(time,1);
K0 = K_value(time,1);
F0 = T_value(time,1);

mean_distort_nordhaus = beta_tilde_1_func(log(R0),F0,log(K0))-beta_f;
lambda_tilde_nordhaus = lambda_tilde_1_func(log(R0),F0,log(K0));

original_dist = normpdf(beta_f_space,beta_f,sqrt(var_beta_f));
nordhaus_dist_0 = normpdf(beta_f_space,mean_distort_nordhaus+beta_f,1./sqrt(lambda_tilde_nordhaus));

nordhaus = nordhaus_dist_0;
original = original_dist;

save('Dist_100yr','nordhaus','original')



