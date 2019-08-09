clear all
close all
clc
% rng(20190430);
pwd;
directory2 = pwd;
 
file1 = [directory2, '/HJB_NonLinGrowth_NoAmb'];
Model1 = load(file1,'RE','R_1','R_2','R_3','R_4','R_5','R_6','R_7','R_8','R_9',...
    'pi_tilde_1_norm','pi_tilde_2_norm','pi_tilde_3_norm','pi_tilde_4_norm','pi_tilde_5_norm',...
    'pi_tilde_6_norm','pi_tilde_7_norm','pi_tilde_8_norm','pi_tilde_9_norm',...
    'r_mat','k_mat','t_mat','v0_dk',...
    'weight1','weight2','weight3','weight4','weight5','weight6','weight7','weight8','weight9',...
    'beta_f','t_max','t_min','theta','var_beta_f','gamma1','gamma2',...
    'gamma_1','gamma_2','power','e','beta_tilde_1','beta_tilde_2','beta_tilde_3','beta_tilde_4',...
    'beta_tilde_5','beta_tilde_6','beta_tilde_7','beta_tilde_8','beta_tilde_9',...
    'lambda_tilde_1','lambda_tilde_2','lambda_tilde_3','lambda_tilde_4',...
    'lambda_tilde_5','lambda_tilde_6','lambda_tilde_7','lambda_tilde_8','lambda_tilde_9',...
    't_bar','a','b'); 

RE = Model1.RE;
RE_1 = Model1.R_1;
RE_2 = Model1.R_2;
RE_3 = Model1.R_3;
RE_4 = Model1.R_4;
RE_5 = Model1.R_5;
RE_6 = Model1.R_6;
RE_7 = Model1.R_7;
RE_8 = Model1.R_8;
RE_9 = Model1.R_9;
weight_1 = Model1.weight1;
weight_2 = Model1.weight2;
weight_3 = Model1.weight3;
weight_4 = Model1.weight4;
weight_5 = Model1.weight5;
weight_6 = Model1.weight6;
weight_7 = Model1.weight7;
weight_8 = Model1.weight8;
weight_9 = Model1.weight9;
pi_tilde_1 = Model1.pi_tilde_1_norm;
pi_tilde_2 = Model1.pi_tilde_2_norm;
pi_tilde_3 = Model1.pi_tilde_3_norm;
pi_tilde_4 = Model1.pi_tilde_4_norm;
pi_tilde_5 = Model1.pi_tilde_5_norm;
pi_tilde_6 = Model1.pi_tilde_6_norm;
pi_tilde_7 = Model1.pi_tilde_7_norm;
pi_tilde_8 = Model1.pi_tilde_8_norm;
pi_tilde_9 = Model1.pi_tilde_9_norm;
% scale_1 = Model1.scale_1;
% scale_2 = Model1.scale_2;
% scale_3 = Model1.scale_3;
% scale_4 = Model1.scale_4;
% scale_5 = Model1.scale_5;
% scale_6 = Model1.scale_6;
% scale_7 = Model1.scale_7;
% scale_8 = Model1.scale_8;
% scale_9 = Model1.scale_9;
beta_f = Model1.beta_f;
var_beta_f = Model1.var_beta_f;
beta_tilde_1 = Model1.beta_tilde_1;
lambda_tilde_1 = Model1.lambda_tilde_1;
beta_tilde_2 = Model1.beta_tilde_2;
lambda_tilde_2 = Model1.lambda_tilde_2;
beta_tilde_3 = Model1.beta_tilde_3;
lambda_tilde_3 = Model1.lambda_tilde_3;
beta_tilde_4 = Model1.beta_tilde_4;
lambda_tilde_4 = Model1.lambda_tilde_4;
beta_tilde_5 = Model1.beta_tilde_5;
lambda_tilde_5 = Model1.lambda_tilde_5;
beta_tilde_6 = Model1.beta_tilde_6;
lambda_tilde_6 = Model1.lambda_tilde_6;
beta_tilde_7 = Model1.beta_tilde_7;
lambda_tilde_7 = Model1.lambda_tilde_7;
beta_tilde_8 = Model1.beta_tilde_8;
lambda_tilde_8 = Model1.lambda_tilde_8;
beta_tilde_9 = Model1.beta_tilde_9;
lambda_tilde_9 = Model1.lambda_tilde_9;

t_max = Model1.t_max;
t_min = Model1.t_min;
theta = Model1.theta;
gamma_1 = Model1.gamma_1;
gamma_2 = Model1.gamma_2;
gamma1 = Model1.gamma1;
gamma2 = Model1.gamma2;
e = Model1.e;
t_bar = Model1.t_bar;
r_mat = Model1.r_mat;
t_mat = Model1.t_mat;
k_mat = Model1.k_mat;
v0_dk = Model1.v0_dk;

a = beta_f-10.*sqrt(var_beta_f);
b = beta_f+10.*sqrt(var_beta_f);

A = Model1.a;
B = Model1.b;

RE_func = griddedInterpolant(r_mat,t_mat,k_mat,RE,'spline');
RE_1_func = griddedInterpolant(r_mat,t_mat,k_mat,RE_1,'spline');
RE_2_func = griddedInterpolant(r_mat,t_mat,k_mat,RE_2,'spline');
RE_3_func = griddedInterpolant(r_mat,t_mat,k_mat,RE_3,'spline');
RE_4_func = griddedInterpolant(r_mat,t_mat,k_mat,RE_4,'spline');
RE_5_func = griddedInterpolant(r_mat,t_mat,k_mat,RE_5,'spline');
RE_6_func = griddedInterpolant(r_mat,t_mat,k_mat,RE_6,'spline');
RE_7_func = griddedInterpolant(r_mat,t_mat,k_mat,RE_7,'spline');
RE_8_func = griddedInterpolant(r_mat,t_mat,k_mat,RE_8,'spline');
RE_9_func = griddedInterpolant(r_mat,t_mat,k_mat,RE_9,'spline');
pi_tilde_1_func = griddedInterpolant(r_mat,t_mat,k_mat,pi_tilde_1,'spline');
pi_tilde_2_func = griddedInterpolant(r_mat,t_mat,k_mat,pi_tilde_2,'spline');
pi_tilde_3_func = griddedInterpolant(r_mat,t_mat,k_mat,pi_tilde_3,'spline');
pi_tilde_4_func = griddedInterpolant(r_mat,t_mat,k_mat,pi_tilde_4,'spline');
pi_tilde_5_func = griddedInterpolant(r_mat,t_mat,k_mat,pi_tilde_5,'spline');
pi_tilde_6_func = griddedInterpolant(r_mat,t_mat,k_mat,pi_tilde_6,'spline');
pi_tilde_7_func = griddedInterpolant(r_mat,t_mat,k_mat,pi_tilde_7,'spline');
pi_tilde_8_func = griddedInterpolant(r_mat,t_mat,k_mat,pi_tilde_8,'spline');
pi_tilde_9_func = griddedInterpolant(r_mat,t_mat,k_mat,pi_tilde_9,'spline');
e_func = griddedInterpolant(r_mat,t_mat,k_mat,e,'spline');
% scale_1_func = griddedInterpolant(r_mat,t_mat,k_mat,scale_1,'spline');
% scale_2_func = griddedInterpolant(r_mat,t_mat,k_mat,scale_2,'spline');
% scale_3_func = griddedInterpolant(r_mat,t_mat,k_mat,scale_3,'spline');
% scale_4_func = griddedInterpolant(r_mat,t_mat,k_mat,scale_4,'spline');
% scale_5_func = griddedInterpolant(r_mat,t_mat,k_mat,scale_5,'spline');
% scale_6_func = griddedInterpolant(r_mat,t_mat,k_mat,scale_6,'spline');
% scale_7_func = griddedInterpolant(r_mat,t_mat,k_mat,scale_7,'spline');
% scale_8_func = griddedInterpolant(r_mat,t_mat,k_mat,scale_8,'spline');
% scale_9_func = griddedInterpolant(r_mat,t_mat,k_mat,scale_9,'spline');
beta_tilde_1_func = griddedInterpolant(r_mat,t_mat,k_mat,beta_tilde_1,'spline');
lambda_tilde_1_func = griddedInterpolant(r_mat,t_mat,k_mat,lambda_tilde_1,'spline');
beta_tilde_2_func = griddedInterpolant(r_mat,t_mat,k_mat,beta_tilde_2,'spline');
lambda_tilde_2_func = griddedInterpolant(r_mat,t_mat,k_mat,lambda_tilde_2,'spline');
beta_tilde_3_func = griddedInterpolant(r_mat,t_mat,k_mat,beta_tilde_3,'spline');
lambda_tilde_3_func = griddedInterpolant(r_mat,t_mat,k_mat,lambda_tilde_3,'spline');
beta_tilde_4_func = griddedInterpolant(r_mat,t_mat,k_mat,beta_tilde_4,'spline');
lambda_tilde_4_func = griddedInterpolant(r_mat,t_mat,k_mat,lambda_tilde_4,'spline');
beta_tilde_5_func = griddedInterpolant(r_mat,t_mat,k_mat,beta_tilde_5,'spline');
lambda_tilde_5_func = griddedInterpolant(r_mat,t_mat,k_mat,lambda_tilde_5,'spline');
beta_tilde_6_func = griddedInterpolant(r_mat,t_mat,k_mat,beta_tilde_6,'spline');
lambda_tilde_6_func = griddedInterpolant(r_mat,t_mat,k_mat,lambda_tilde_6,'spline');
beta_tilde_7_func = griddedInterpolant(r_mat,t_mat,k_mat,beta_tilde_7,'spline');
lambda_tilde_7_func = griddedInterpolant(r_mat,t_mat,k_mat,lambda_tilde_7,'spline');
beta_tilde_8_func = griddedInterpolant(r_mat,t_mat,k_mat,beta_tilde_8,'spline');
lambda_tilde_8_func = griddedInterpolant(r_mat,t_mat,k_mat,lambda_tilde_8,'spline');
beta_tilde_9_func = griddedInterpolant(r_mat,t_mat,k_mat,beta_tilde_9,'spline');
lambda_tilde_9_func = griddedInterpolant(r_mat,t_mat,k_mat,lambda_tilde_9,'spline');

file1 = [directory2, '/HJB_NonLinGrowth_NoAmb_Sims_determ'];
Model1 = load(file1,'hists2','e_hists2','v_dk_hists2'); 
hists2 = Model1.hists2;
v_dk_hists2 = Model1.v_dk_hists2;
e = Model1.e_hists2;
T_value = mean(squeeze(hists2(:,3,:)),2);
R_value = mean(squeeze(hists2(:,1,:)),2);
K_value = mean(squeeze(hists2(:,2,:)),2);
V_DK_value = mean(squeeze(v_dk_hists2(:,:)),2);
E_value = mean(e,2);

time_vec = linspace(0,100,400);


%% Plot_20190718


for time=1:400
    RE_plot(time) = RE_func(log(R_value(time,1)),T_value(time,1),log(K_value(time,1)));
    weight1_plot(time) = pi_tilde_1_func(log(R_value(time,1)),T_value(time,1),log(K_value(time,1)));
    weight2_plot(time) = pi_tilde_2_func(log(R_value(time,1)),T_value(time,1),log(K_value(time,1)));
    weight3_plot(time) = pi_tilde_3_func(log(R_value(time,1)),T_value(time,1),log(K_value(time,1)));
    weight4_plot(time) = pi_tilde_4_func(log(R_value(time,1)),T_value(time,1),log(K_value(time,1)));
    weight5_plot(time) = pi_tilde_5_func(log(R_value(time,1)),T_value(time,1),log(K_value(time,1)));
    weight6_plot(time) = pi_tilde_6_func(log(R_value(time,1)),T_value(time,1),log(K_value(time,1)));
    weight7_plot(time) = pi_tilde_7_func(log(R_value(time,1)),T_value(time,1),log(K_value(time,1)));
    weight8_plot(time) = pi_tilde_8_func(log(R_value(time,1)),T_value(time,1),log(K_value(time,1)));
    weight9_plot(time) = pi_tilde_9_func(log(R_value(time,1)),T_value(time,1),log(K_value(time,1)));

end


fileID = fopen('Relative Entropy.txt','w');
fprintf(fileID,'xi_a: %.6f \n',1./theta);
fprintf(fileID,'0 yr: %.6f \n',RE_plot(1));
fprintf(fileID,'25 yr: %.6f \n',RE_plot(100));
fprintf(fileID,'50 yr: %.6f \n',RE_plot(200));
fprintf(fileID,'75 yr: %.6f \n',RE_plot(300));
fprintf(fileID,'100 yr: %.6f \n',RE_plot(400));
fclose(fileID);


fileID = fopen('Weights File.txt','w');
fprintf(fileID,'xi_a: %.6f \n',1./theta);
fprintf(fileID,'Tilted Weights: \n',1./theta);
fprintf(fileID,'0 yr weight 1: %.6f \n',weight1_plot(1));
fprintf(fileID,'25 yr weight 1: %.6f \n',weight1_plot(100));
fprintf(fileID,'50 yr weight 1: %.6f \n',weight1_plot(200));
fprintf(fileID,'75 yr weight 1: %.6f \n',weight1_plot(300));
fprintf(fileID,'100 yr weight 1: %.6f \n',weight1_plot(400));
fprintf(fileID,'0 yr weight 2: %.6f \n',weight2_plot(1));
fprintf(fileID,'25 yr weight 2: %.6f \n',weight2_plot(100));
fprintf(fileID,'50 yr weight 2: %.6f \n',weight2_plot(200));
fprintf(fileID,'75 yr weight 2: %.6f \n',weight2_plot(300));
fprintf(fileID,'100 yr weight 2: %.6f \n',weight2_plot(400));
fprintf(fileID,'0 yr weight 3: %.6f \n',weight3_plot(1));
fprintf(fileID,'25 yr weight 3: %.6f \n',weight3_plot(100));
fprintf(fileID,'50 yr weight 3: %.6f \n',weight3_plot(200));
fprintf(fileID,'75 yr weight 3: %.6f \n',weight3_plot(300));
fprintf(fileID,'100 yr weight 3: %.6f \n',weight3_plot(400));
fprintf(fileID,'0 yr weight 4: %.6f \n',weight4_plot(1));
fprintf(fileID,'25 yr weight 4: %.6f \n',weight4_plot(100));
fprintf(fileID,'50 yr weight 4: %.6f \n',weight4_plot(200));
fprintf(fileID,'75 yr weight 4: %.6f \n',weight4_plot(300));
fprintf(fileID,'100 yr weight 4: %.6f \n',weight4_plot(400));
fprintf(fileID,'0 yr weight 5: %.6f \n',weight5_plot(1));
fprintf(fileID,'25 yr weight 5: %.6f \n',weight5_plot(100));
fprintf(fileID,'50 yr weight 5: %.6f \n',weight5_plot(200));
fprintf(fileID,'75 yr weight 5: %.6f \n',weight5_plot(300));
fprintf(fileID,'100 yr weight 5: %.6f \n',weight5_plot(400));
fprintf(fileID,'0 yr weight 6: %.6f \n',weight6_plot(1));
fprintf(fileID,'25 yr weight 6: %.6f \n',weight6_plot(100));
fprintf(fileID,'50 yr weight 6: %.6f \n',weight6_plot(200));
fprintf(fileID,'75 yr weight 6: %.6f \n',weight6_plot(300));
fprintf(fileID,'100 yr weight 6: %.6f \n',weight6_plot(400));
fprintf(fileID,'0 yr weight 7: %.6f \n',weight7_plot(1));
fprintf(fileID,'25 yr weight 7: %.6f \n',weight7_plot(100));
fprintf(fileID,'50 yr weight 7: %.6f \n',weight7_plot(200));
fprintf(fileID,'75 yr weight 7: %.6f \n',weight7_plot(300));
fprintf(fileID,'100 yr weight 7: %.6f \n',weight7_plot(400));
fprintf(fileID,'0 yr weight 8: %.6f \n',weight8_plot(1));
fprintf(fileID,'25 yr weight 8: %.6f \n',weight8_plot(100));
fprintf(fileID,'50 yr weight 8: %.6f \n',weight8_plot(200));
fprintf(fileID,'75 yr weight 8: %.6f \n',weight8_plot(300));
fprintf(fileID,'100 yr weight 8: %.6f \n',weight8_plot(400));
fprintf(fileID,'0 yr weight 9: %.6f \n',weight9_plot(1));
fprintf(fileID,'25 yr weight 9: %.6f \n',weight9_plot(100));
fprintf(fileID,'50 yr weight 9: %.6f \n',weight9_plot(200));
fprintf(fileID,'75 yr weight 9: %.6f \n',weight9_plot(300));
fprintf(fileID,'100 yr weight 9: %.6f \n',weight9_plot(400));
fprintf(fileID,'Baseline Weights: \n');
fprintf(fileID,'baseline weight 1: %.6f \n',weight_1);
fprintf(fileID,'baseline weight 2: %.6f \n',weight_2);
fprintf(fileID,'baseline weight 3: %.6f \n',weight_3);
fprintf(fileID,'baseline weight 4: %.6f \n',weight_4);
fprintf(fileID,'baseline weight 5: %.6f \n',weight_5);
fprintf(fileID,'baseline weight 6: %.6f \n',weight_6);
fprintf(fileID,'baseline weight 7: %.6f \n',weight_7);
fprintf(fileID,'baseline weight 8: %.6f \n',weight_8);
fprintf(fileID,'baseline weight 9: %.6f \n',weight_9);
fclose(fileID);

beta_f_space = linspace(a,b,200);
save('beta_f_space','beta_f_space')

%% T=0
time=1;

R0 = R_value(time,1);
K0 = K_value(time,1);
F0 = T_value(time,1);

weight1 = pi_tilde_1_func(log(R0),F0,log(K0));
weight2 = pi_tilde_2_func(log(R0),F0,log(K0));
weight3 = pi_tilde_3_func(log(R0),F0,log(K0));
weight4 = pi_tilde_4_func(log(R0),F0,log(K0));
weight5 = pi_tilde_5_func(log(R0),F0,log(K0));
weight6 = pi_tilde_6_func(log(R0),F0,log(K0));
weight7 = pi_tilde_7_func(log(R0),F0,log(K0));
weight8 = pi_tilde_8_func(log(R0),F0,log(K0));
weight9 = pi_tilde_9_func(log(R0),F0,log(K0));
    
mean_distort_1 = beta_tilde_1_func(log(R0),F0,log(K0))-beta_f;
lambda_tilde_1 = lambda_tilde_1_func(log(R0),F0,log(K0));
mean_distort_2 = beta_tilde_2_func(log(R0),F0,log(K0))-beta_f;
lambda_tilde_2 = lambda_tilde_2_func(log(R0),F0,log(K0));
mean_distort_3 = beta_tilde_3_func(log(R0),F0,log(K0))-beta_f;
lambda_tilde_3 = lambda_tilde_3_func(log(R0),F0,log(K0));
mean_distort_4 = beta_tilde_4_func(log(R0),F0,log(K0))-beta_f;
lambda_tilde_4 = lambda_tilde_4_func(log(R0),F0,log(K0));
mean_distort_5 = beta_tilde_5_func(log(R0),F0,log(K0))-beta_f;
lambda_tilde_5 = lambda_tilde_5_func(log(R0),F0,log(K0));
mean_distort_6 = beta_tilde_6_func(log(R0),F0,log(K0))-beta_f;
lambda_tilde_6 = lambda_tilde_6_func(log(R0),F0,log(K0));
mean_distort_7 = beta_tilde_7_func(log(R0),F0,log(K0))-beta_f;
lambda_tilde_7 = lambda_tilde_7_func(log(R0),F0,log(K0));
mean_distort_8 = beta_tilde_8_func(log(R0),F0,log(K0))-beta_f;
lambda_tilde_8 = lambda_tilde_8_func(log(R0),F0,log(K0));
mean_distort_9 = beta_tilde_9_func(log(R0),F0,log(K0))-beta_f;
lambda_tilde_9 = lambda_tilde_9_func(log(R0),F0,log(K0));


% weight not right?
original_dist = normpdf(beta_f_space,beta_f,sqrt(var_beta_f));
tilt_dist_1 = normpdf(beta_f_space,mean_distort_1+beta_f,1./sqrt(lambda_tilde_1));
tilt_dist_2 = normpdf(beta_f_space,mean_distort_2+beta_f,1./sqrt(lambda_tilde_2));
tilt_dist_3 = normpdf(beta_f_space,mean_distort_3+beta_f,1./sqrt(lambda_tilde_3));
tilt_dist_4 = normpdf(beta_f_space,mean_distort_4+beta_f,1./sqrt(lambda_tilde_4));
tilt_dist_5 = normpdf(beta_f_space,mean_distort_5+beta_f,1./sqrt(lambda_tilde_5));
tilt_dist_6 = normpdf(beta_f_space,mean_distort_6+beta_f,1./sqrt(lambda_tilde_6));
tilt_dist_7 = normpdf(beta_f_space,mean_distort_7+beta_f,1./sqrt(lambda_tilde_7));
tilt_dist_8 = normpdf(beta_f_space,mean_distort_8+beta_f,1./sqrt(lambda_tilde_8));
tilt_dist_9 = normpdf(beta_f_space,mean_distort_9+beta_f,1./sqrt(lambda_tilde_9));
weight1_0 = weight1;
weight2_0 = weight2;
weight3_0 = weight3;
weight4_0 = weight4;
weight5_0 = weight5;
weight6_0 = weight6;
weight7_0 = weight7;
weight8_0 = weight8;
weight9_0 = weight9;

model1_0 = tilt_dist_1;
model2_0 = tilt_dist_2;
model3_0 = tilt_dist_3;
model4_0 = tilt_dist_4;
model5_0 = tilt_dist_5;
model6_0 = tilt_dist_6;
model7_0 = tilt_dist_7;
model8_0 = tilt_dist_8;
model9_0 = tilt_dist_9;
original = original_dist;
weighted = weight1_0.*model1_0+weight2_0.*model2_0+weight3_0.*model3_0+weight4_0.*model4_0...
    +weight5_0.*model5_0+weight6_0.*model6_0+weight7_0.*model7_0+weight8_0.*model8_0+weight9_0.*model9_0;


save('Dist_0yr','model1_0','model2_0','model3_0','model4_0','model5_0','model6_0','model7_0','model8_0'...
    ,'model9_0','original','weighted')


%% T=25
time=100;


R0 = R_value(time,1);
K0 = K_value(time,1);
F0 = T_value(time,1);

weight1 = pi_tilde_1_func(log(R0),F0,log(K0));
weight2 = pi_tilde_2_func(log(R0),F0,log(K0));
weight3 = pi_tilde_3_func(log(R0),F0,log(K0));
weight4 = pi_tilde_4_func(log(R0),F0,log(K0));
weight5 = pi_tilde_5_func(log(R0),F0,log(K0));
weight6 = pi_tilde_6_func(log(R0),F0,log(K0));
weight7 = pi_tilde_7_func(log(R0),F0,log(K0));
weight8 = pi_tilde_8_func(log(R0),F0,log(K0));
weight9 = pi_tilde_9_func(log(R0),F0,log(K0));
    
mean_distort_1 = beta_tilde_1_func(log(R0),F0,log(K0))-beta_f;
lambda_tilde_1 = lambda_tilde_1_func(log(R0),F0,log(K0));
mean_distort_2 = beta_tilde_2_func(log(R0),F0,log(K0))-beta_f;
lambda_tilde_2 = lambda_tilde_2_func(log(R0),F0,log(K0));
mean_distort_3 = beta_tilde_3_func(log(R0),F0,log(K0))-beta_f;
lambda_tilde_3 = lambda_tilde_3_func(log(R0),F0,log(K0));
mean_distort_4 = beta_tilde_4_func(log(R0),F0,log(K0))-beta_f;
lambda_tilde_4 = lambda_tilde_4_func(log(R0),F0,log(K0));
mean_distort_5 = beta_tilde_5_func(log(R0),F0,log(K0))-beta_f;
lambda_tilde_5 = lambda_tilde_5_func(log(R0),F0,log(K0));
mean_distort_6 = beta_tilde_6_func(log(R0),F0,log(K0))-beta_f;
lambda_tilde_6 = lambda_tilde_6_func(log(R0),F0,log(K0));
mean_distort_7 = beta_tilde_7_func(log(R0),F0,log(K0))-beta_f;
lambda_tilde_7 = lambda_tilde_7_func(log(R0),F0,log(K0));
mean_distort_8 = beta_tilde_8_func(log(R0),F0,log(K0))-beta_f;
lambda_tilde_8 = lambda_tilde_8_func(log(R0),F0,log(K0));
mean_distort_9 = beta_tilde_9_func(log(R0),F0,log(K0))-beta_f;
lambda_tilde_9 = lambda_tilde_9_func(log(R0),F0,log(K0));


% weight not right?
original_dist = normpdf(beta_f_space,beta_f,sqrt(var_beta_f));
tilt_dist_1 = normpdf(beta_f_space,mean_distort_1+beta_f,1./sqrt(lambda_tilde_1));
tilt_dist_2 = normpdf(beta_f_space,mean_distort_2+beta_f,1./sqrt(lambda_tilde_2));
tilt_dist_3 = normpdf(beta_f_space,mean_distort_3+beta_f,1./sqrt(lambda_tilde_3));
tilt_dist_4 = normpdf(beta_f_space,mean_distort_4+beta_f,1./sqrt(lambda_tilde_4));
tilt_dist_5 = normpdf(beta_f_space,mean_distort_5+beta_f,1./sqrt(lambda_tilde_5));
tilt_dist_6 = normpdf(beta_f_space,mean_distort_6+beta_f,1./sqrt(lambda_tilde_6));
tilt_dist_7 = normpdf(beta_f_space,mean_distort_7+beta_f,1./sqrt(lambda_tilde_7));
tilt_dist_8 = normpdf(beta_f_space,mean_distort_8+beta_f,1./sqrt(lambda_tilde_8));
tilt_dist_9 = normpdf(beta_f_space,mean_distort_9+beta_f,1./sqrt(lambda_tilde_9));
weight1_0 = weight1;
weight2_0 = weight2;
weight3_0 = weight3;
weight4_0 = weight4;
weight5_0 = weight5;
weight6_0 = weight6;
weight7_0 = weight7;
weight8_0 = weight8;
weight9_0 = weight9;

model1_0 = tilt_dist_1;
model2_0 = tilt_dist_2;
model3_0 = tilt_dist_3;
model4_0 = tilt_dist_4;
model5_0 = tilt_dist_5;
model6_0 = tilt_dist_6;
model7_0 = tilt_dist_7;
model8_0 = tilt_dist_8;
model9_0 = tilt_dist_9;
original = original_dist;
weighted = weight1_0.*model1_0+weight2_0.*model2_0+weight3_0.*model3_0+weight4_0.*model4_0...
    +weight5_0.*model5_0+weight6_0.*model6_0+weight7_0.*model7_0+weight8_0.*model8_0+weight9_0.*model9_0;


save('Dist_25yr','model1_0','model2_0','model3_0','model4_0','model5_0','model6_0','model7_0','model8_0'...
    ,'model9_0','original','weighted')

%% T=50
time=200;

R0 = R_value(time,1);
K0 = K_value(time,1);
F0 = T_value(time,1);

weight1 = pi_tilde_1_func(log(R0),F0,log(K0));
weight2 = pi_tilde_2_func(log(R0),F0,log(K0));
weight3 = pi_tilde_3_func(log(R0),F0,log(K0));
weight4 = pi_tilde_4_func(log(R0),F0,log(K0));
weight5 = pi_tilde_5_func(log(R0),F0,log(K0));
weight6 = pi_tilde_6_func(log(R0),F0,log(K0));
weight7 = pi_tilde_7_func(log(R0),F0,log(K0));
weight8 = pi_tilde_8_func(log(R0),F0,log(K0));
weight9 = pi_tilde_9_func(log(R0),F0,log(K0));
    
mean_distort_1 = beta_tilde_1_func(log(R0),F0,log(K0))-beta_f;
lambda_tilde_1 = lambda_tilde_1_func(log(R0),F0,log(K0));
mean_distort_2 = beta_tilde_2_func(log(R0),F0,log(K0))-beta_f;
lambda_tilde_2 = lambda_tilde_2_func(log(R0),F0,log(K0));
mean_distort_3 = beta_tilde_3_func(log(R0),F0,log(K0))-beta_f;
lambda_tilde_3 = lambda_tilde_3_func(log(R0),F0,log(K0));
mean_distort_4 = beta_tilde_4_func(log(R0),F0,log(K0))-beta_f;
lambda_tilde_4 = lambda_tilde_4_func(log(R0),F0,log(K0));
mean_distort_5 = beta_tilde_5_func(log(R0),F0,log(K0))-beta_f;
lambda_tilde_5 = lambda_tilde_5_func(log(R0),F0,log(K0));
mean_distort_6 = beta_tilde_6_func(log(R0),F0,log(K0))-beta_f;
lambda_tilde_6 = lambda_tilde_6_func(log(R0),F0,log(K0));
mean_distort_7 = beta_tilde_7_func(log(R0),F0,log(K0))-beta_f;
lambda_tilde_7 = lambda_tilde_7_func(log(R0),F0,log(K0));
mean_distort_8 = beta_tilde_8_func(log(R0),F0,log(K0))-beta_f;
lambda_tilde_8 = lambda_tilde_8_func(log(R0),F0,log(K0));
mean_distort_9 = beta_tilde_9_func(log(R0),F0,log(K0))-beta_f;
lambda_tilde_9 = lambda_tilde_9_func(log(R0),F0,log(K0));


% weight not right?
original_dist = normpdf(beta_f_space,beta_f,sqrt(var_beta_f));
tilt_dist_1 = normpdf(beta_f_space,mean_distort_1+beta_f,1./sqrt(lambda_tilde_1));
tilt_dist_2 = normpdf(beta_f_space,mean_distort_2+beta_f,1./sqrt(lambda_tilde_2));
tilt_dist_3 = normpdf(beta_f_space,mean_distort_3+beta_f,1./sqrt(lambda_tilde_3));
tilt_dist_4 = normpdf(beta_f_space,mean_distort_4+beta_f,1./sqrt(lambda_tilde_4));
tilt_dist_5 = normpdf(beta_f_space,mean_distort_5+beta_f,1./sqrt(lambda_tilde_5));
tilt_dist_6 = normpdf(beta_f_space,mean_distort_6+beta_f,1./sqrt(lambda_tilde_6));
tilt_dist_7 = normpdf(beta_f_space,mean_distort_7+beta_f,1./sqrt(lambda_tilde_7));
tilt_dist_8 = normpdf(beta_f_space,mean_distort_8+beta_f,1./sqrt(lambda_tilde_8));
tilt_dist_9 = normpdf(beta_f_space,mean_distort_9+beta_f,1./sqrt(lambda_tilde_9));
weight1_0 = weight1;
weight2_0 = weight2;
weight3_0 = weight3;
weight4_0 = weight4;
weight5_0 = weight5;
weight6_0 = weight6;
weight7_0 = weight7;
weight8_0 = weight8;
weight9_0 = weight9;

model1_0 = tilt_dist_1;
model2_0 = tilt_dist_2;
model3_0 = tilt_dist_3;
model4_0 = tilt_dist_4;
model5_0 = tilt_dist_5;
model6_0 = tilt_dist_6;
model7_0 = tilt_dist_7;
model8_0 = tilt_dist_8;
model9_0 = tilt_dist_9;
original = original_dist;
weighted = weight1_0.*model1_0+weight2_0.*model2_0+weight3_0.*model3_0+weight4_0.*model4_0...
    +weight5_0.*model5_0+weight6_0.*model6_0+weight7_0.*model7_0+weight8_0.*model8_0+weight9_0.*model9_0;


save('Dist_50yr','model1_0','model2_0','model3_0','model4_0','model5_0','model6_0','model7_0','model8_0'...
    ,'model9_0','original','weighted')

%% T=75
time=300;

R0 = R_value(time,1);
K0 = K_value(time,1);
F0 = T_value(time,1);

weight1 = pi_tilde_1_func(log(R0),F0,log(K0));
weight2 = pi_tilde_2_func(log(R0),F0,log(K0));
weight3 = pi_tilde_3_func(log(R0),F0,log(K0));
weight4 = pi_tilde_4_func(log(R0),F0,log(K0));
weight5 = pi_tilde_5_func(log(R0),F0,log(K0));
weight6 = pi_tilde_6_func(log(R0),F0,log(K0));
weight7 = pi_tilde_7_func(log(R0),F0,log(K0));
weight8 = pi_tilde_8_func(log(R0),F0,log(K0));
weight9 = pi_tilde_9_func(log(R0),F0,log(K0));
    
mean_distort_1 = beta_tilde_1_func(log(R0),F0,log(K0))-beta_f;
lambda_tilde_1 = lambda_tilde_1_func(log(R0),F0,log(K0));
mean_distort_2 = beta_tilde_2_func(log(R0),F0,log(K0))-beta_f;
lambda_tilde_2 = lambda_tilde_2_func(log(R0),F0,log(K0));
mean_distort_3 = beta_tilde_3_func(log(R0),F0,log(K0))-beta_f;
lambda_tilde_3 = lambda_tilde_3_func(log(R0),F0,log(K0));
mean_distort_4 = beta_tilde_4_func(log(R0),F0,log(K0))-beta_f;
lambda_tilde_4 = lambda_tilde_4_func(log(R0),F0,log(K0));
mean_distort_5 = beta_tilde_5_func(log(R0),F0,log(K0))-beta_f;
lambda_tilde_5 = lambda_tilde_5_func(log(R0),F0,log(K0));
mean_distort_6 = beta_tilde_6_func(log(R0),F0,log(K0))-beta_f;
lambda_tilde_6 = lambda_tilde_6_func(log(R0),F0,log(K0));
mean_distort_7 = beta_tilde_7_func(log(R0),F0,log(K0))-beta_f;
lambda_tilde_7 = lambda_tilde_7_func(log(R0),F0,log(K0));
mean_distort_8 = beta_tilde_8_func(log(R0),F0,log(K0))-beta_f;
lambda_tilde_8 = lambda_tilde_8_func(log(R0),F0,log(K0));
mean_distort_9 = beta_tilde_9_func(log(R0),F0,log(K0))-beta_f;
lambda_tilde_9 = lambda_tilde_9_func(log(R0),F0,log(K0));


% weight not right?
original_dist = normpdf(beta_f_space,beta_f,sqrt(var_beta_f));
tilt_dist_1 = normpdf(beta_f_space,mean_distort_1+beta_f,1./sqrt(lambda_tilde_1));
tilt_dist_2 = normpdf(beta_f_space,mean_distort_2+beta_f,1./sqrt(lambda_tilde_2));
tilt_dist_3 = normpdf(beta_f_space,mean_distort_3+beta_f,1./sqrt(lambda_tilde_3));
tilt_dist_4 = normpdf(beta_f_space,mean_distort_4+beta_f,1./sqrt(lambda_tilde_4));
tilt_dist_5 = normpdf(beta_f_space,mean_distort_5+beta_f,1./sqrt(lambda_tilde_5));
tilt_dist_6 = normpdf(beta_f_space,mean_distort_6+beta_f,1./sqrt(lambda_tilde_6));
tilt_dist_7 = normpdf(beta_f_space,mean_distort_7+beta_f,1./sqrt(lambda_tilde_7));
tilt_dist_8 = normpdf(beta_f_space,mean_distort_8+beta_f,1./sqrt(lambda_tilde_8));
tilt_dist_9 = normpdf(beta_f_space,mean_distort_9+beta_f,1./sqrt(lambda_tilde_9));
weight1_0 = weight1;
weight2_0 = weight2;
weight3_0 = weight3;
weight4_0 = weight4;
weight5_0 = weight5;
weight6_0 = weight6;
weight7_0 = weight7;
weight8_0 = weight8;
weight9_0 = weight9;

model1_0 = tilt_dist_1;
model2_0 = tilt_dist_2;
model3_0 = tilt_dist_3;
model4_0 = tilt_dist_4;
model5_0 = tilt_dist_5;
model6_0 = tilt_dist_6;
model7_0 = tilt_dist_7;
model8_0 = tilt_dist_8;
model9_0 = tilt_dist_9;
original = original_dist;
weighted = weight1_0.*model1_0+weight2_0.*model2_0+weight3_0.*model3_0+weight4_0.*model4_0...
    +weight5_0.*model5_0+weight6_0.*model6_0+weight7_0.*model7_0+weight8_0.*model8_0+weight9_0.*model9_0;


save('Dist_75yr','model1_0','model2_0','model3_0','model4_0','model5_0','model6_0','model7_0','model8_0'...
    ,'model9_0','original','weighted')

%% T=100
time=400;

R0 = R_value(time,1);
K0 = K_value(time,1);
F0 = T_value(time,1);

weight1 = pi_tilde_1_func(log(R0),F0,log(K0));
weight2 = pi_tilde_2_func(log(R0),F0,log(K0));
weight3 = pi_tilde_3_func(log(R0),F0,log(K0));
weight4 = pi_tilde_4_func(log(R0),F0,log(K0));
weight5 = pi_tilde_5_func(log(R0),F0,log(K0));
weight6 = pi_tilde_6_func(log(R0),F0,log(K0));
weight7 = pi_tilde_7_func(log(R0),F0,log(K0));
weight8 = pi_tilde_8_func(log(R0),F0,log(K0));
weight9 = pi_tilde_9_func(log(R0),F0,log(K0));
    
mean_distort_1 = beta_tilde_1_func(log(R0),F0,log(K0))-beta_f;
lambda_tilde_1 = lambda_tilde_1_func(log(R0),F0,log(K0));
mean_distort_2 = beta_tilde_2_func(log(R0),F0,log(K0))-beta_f;
lambda_tilde_2 = lambda_tilde_2_func(log(R0),F0,log(K0));
mean_distort_3 = beta_tilde_3_func(log(R0),F0,log(K0))-beta_f;
lambda_tilde_3 = lambda_tilde_3_func(log(R0),F0,log(K0));
mean_distort_4 = beta_tilde_4_func(log(R0),F0,log(K0))-beta_f;
lambda_tilde_4 = lambda_tilde_4_func(log(R0),F0,log(K0));
mean_distort_5 = beta_tilde_5_func(log(R0),F0,log(K0))-beta_f;
lambda_tilde_5 = lambda_tilde_5_func(log(R0),F0,log(K0));
mean_distort_6 = beta_tilde_6_func(log(R0),F0,log(K0))-beta_f;
lambda_tilde_6 = lambda_tilde_6_func(log(R0),F0,log(K0));
mean_distort_7 = beta_tilde_7_func(log(R0),F0,log(K0))-beta_f;
lambda_tilde_7 = lambda_tilde_7_func(log(R0),F0,log(K0));
mean_distort_8 = beta_tilde_8_func(log(R0),F0,log(K0))-beta_f;
lambda_tilde_8 = lambda_tilde_8_func(log(R0),F0,log(K0));
mean_distort_9 = beta_tilde_9_func(log(R0),F0,log(K0))-beta_f;
lambda_tilde_9 = lambda_tilde_9_func(log(R0),F0,log(K0));


% weight not right?
original_dist = normpdf(beta_f_space,beta_f,sqrt(var_beta_f));
tilt_dist_1 = normpdf(beta_f_space,mean_distort_1+beta_f,1./sqrt(lambda_tilde_1));
tilt_dist_2 = normpdf(beta_f_space,mean_distort_2+beta_f,1./sqrt(lambda_tilde_2));
tilt_dist_3 = normpdf(beta_f_space,mean_distort_3+beta_f,1./sqrt(lambda_tilde_3));
tilt_dist_4 = normpdf(beta_f_space,mean_distort_4+beta_f,1./sqrt(lambda_tilde_4));
tilt_dist_5 = normpdf(beta_f_space,mean_distort_5+beta_f,1./sqrt(lambda_tilde_5));
tilt_dist_6 = normpdf(beta_f_space,mean_distort_6+beta_f,1./sqrt(lambda_tilde_6));
tilt_dist_7 = normpdf(beta_f_space,mean_distort_7+beta_f,1./sqrt(lambda_tilde_7));
tilt_dist_8 = normpdf(beta_f_space,mean_distort_8+beta_f,1./sqrt(lambda_tilde_8));
tilt_dist_9 = normpdf(beta_f_space,mean_distort_9+beta_f,1./sqrt(lambda_tilde_9));
weight1_0 = weight1;
weight2_0 = weight2;
weight3_0 = weight3;
weight4_0 = weight4;
weight5_0 = weight5;
weight6_0 = weight6;
weight7_0 = weight7;
weight8_0 = weight8;
weight9_0 = weight9;

model1_0 = tilt_dist_1;
model2_0 = tilt_dist_2;
model3_0 = tilt_dist_3;
model4_0 = tilt_dist_4;
model5_0 = tilt_dist_5;
model6_0 = tilt_dist_6;
model7_0 = tilt_dist_7;
model8_0 = tilt_dist_8;
model9_0 = tilt_dist_9;
original = original_dist;
weighted = weight1_0.*model1_0+weight2_0.*model2_0+weight3_0.*model3_0+weight4_0.*model4_0...
    +weight5_0.*model5_0+weight6_0.*model6_0+weight7_0.*model7_0+weight8_0.*model8_0+weight9_0.*model9_0;


save('Dist_100yr','model1_0','model2_0','model3_0','model4_0','model5_0','model6_0','model7_0','model8_0'...
    ,'model9_0','original','weighted')
