%%%%% This file generates numbers in left column from Table G.1 
%%%%% for preference setting in the online appendix.

% Authors: Mike Barnett, Jieyao Wang
% Last update: Nov 25, 2019

close all
clear all
clc

% load results
file1 = [pwd, '/averse_weighted'];
Model1 = load(file1,'RE','pi_tilde_1_norm','pi_tilde_2_norm','r_mat','k_mat','F_mat',...
    'beta_tilde_1','beta_f','lambda_tilde_1','xi_p','var_beta_f','xi_d',...
    'gamma_1','gamma_2','gamma_2_plus','power','gamma_bar','e','n','a','b'); 

RE = Model1.RE;
pi_tilde_1 = Model1.pi_tilde_1_norm;
pi_tilde_2 = Model1.pi_tilde_2_norm;
beta_tilde_1 = Model1.beta_tilde_1;
beta_f = Model1.beta_f;
lambda_tilde_1 = Model1.lambda_tilde_1;
xi_p = Model1.xi_p;
var_beta_f = Model1.var_beta_f;
xi_d = Model1.xi_d;
gamma_1 = Model1.gamma_1;
gamma_2 = Model1.gamma_2;
gamma_2_plus = Model1.gamma_2_plus;
power = Model1.power;
gamma_bar = Model1.gamma_bar;
e = Model1.e;
n = Model1.n;
r_mat = Model1.r_mat;
F_mat = Model1.F_mat;
k_mat = Model1.k_mat;
a = beta_f-10.*sqrt(var_beta_f);
b = beta_f+10.*sqrt(var_beta_f);
A = Model1.a;
B = Model1.b;


RE_func = griddedInterpolant(r_mat,F_mat,k_mat,RE,'linear');
pi_tilde_1_func = griddedInterpolant(r_mat,F_mat,k_mat,pi_tilde_1,'linear');
pi_tilde_2_func = griddedInterpolant(r_mat,F_mat,k_mat,pi_tilde_2,'linear');
beta_tilde_1_func = griddedInterpolant(r_mat,F_mat,k_mat,beta_tilde_1,'linear');
lambda_tilde_1_func = griddedInterpolant(r_mat,F_mat,k_mat,lambda_tilde_1,'linear');
e_func = griddedInterpolant(r_mat,F_mat,k_mat,e,'linear');


file1 = [pwd, '/SS_competitive_Sims'];
Model1 = load(file1,'hists2'); 
hists2 = Model1.hists2;
T_value = mean(squeeze(hists2(:,3,:)),2);
R_value = mean(squeeze(hists2(:,1,:)),2);
K_value = mean(squeeze(hists2(:,2,:)),2);

for time=1:400
    RE_plot(time) = RE_func(log(R_value(time,1)),T_value(time,1),log(K_value(time,1)));
    weight_plot(time) = pi_tilde_1_func(log(R_value(time,1)),T_value(time,1),log(K_value(time,1)));

end

% Generate results of left column of Table G.1

fileID = fopen('Relative Entropy_we.txt','w');
fprintf(fileID,'xi_a: %.6f \n',xi_p);
fprintf(fileID,'0 yr: %.6f \n',RE_plot(1));
fprintf(fileID,'30 yr: %.6f \n',RE_plot(120));
fprintf(fileID,'60 yr: %.6f \n',RE_plot(240));
fclose(fileID);

fileID = fopen('Nordhaus Weight_we.txt','w');
fprintf(fileID,'xi_a: %.6f \n',xi_p);
fprintf(fileID,'0 yr: %.6f \n',weight_plot(1));
fprintf(fileID,'30 yr: %.6f \n',weight_plot(120));
fprintf(fileID,'60 yr: %.6f \n',weight_plot(240));
fclose(fileID);

