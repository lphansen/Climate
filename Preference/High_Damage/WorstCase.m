%%%%%%%%%%%%%%%%%%%%%%%%% Worst case probabilities %%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc

%%%%% Step 0: Load file
file1 = [pwd, '/HJB_NonLinPref_Cumu'];
Model1 = load(file1,'R_2','pi_tilde_1_norm','r_mat','k_mat','t_mat',...
    'beta_tilde_1','beta_f','lambda_tilde_1','theta','var_beta_f','xi_d',...
    'gamma_1','gamma_2','gamma_2_plus','power','e','f_bar','n','a','b'); 

RE_2 = Model1.R_2;
pi_tilde_1 = Model1.pi_tilde_1_norm;
beta_tilde_1 = Model1.beta_tilde_1;
beta_f = Model1.beta_f;
lambda_tilde_1 = Model1.lambda_tilde_1;
theta = Model1.theta;
var_beta_f = Model1.var_beta_f;
xi_d = Model1.xi_d;
gamma_1 = Model1.gamma_1;
gamma_2 = Model1.gamma_2;
gamma_2_plus = Model1.gamma_2_plus;
power = Model1.power;
e = Model1.e;
f_bar = Model1.f_bar;
n = Model1.n;
r_mat = Model1.r_mat;
t_mat = Model1.t_mat;
k_mat = Model1.k_mat;
a = beta_f-10.*sqrt(var_beta_f);
b = beta_f+10.*sqrt(var_beta_f);
A = Model1.a;
B = Model1.b;

%%%%% Step 1: Create function handles
RE_2_func = griddedInterpolant(r_mat,t_mat,k_mat,RE_2,'spline');
pi_tilde_1_func = griddedInterpolant(r_mat,t_mat,k_mat,pi_tilde_1,'spline');
beta_tilde_1_func = griddedInterpolant(r_mat,t_mat,k_mat,beta_tilde_1,'spline');
lambda_tilde_1_func = griddedInterpolant(r_mat,t_mat,k_mat,lambda_tilde_1,'spline');
e_func = griddedInterpolant(r_mat,t_mat,k_mat,e,'spline');

%%%%% Step 2: Load in simulated trajectory
file1 = [pwd, '/HJB_NonLinPref_Cumu_Sims'];
Model1 = load(file1,'hists2'); 
hists2 = Model1.hists2;
T_value = mean(squeeze(hists2(:,3,:)),2);
R_value = mean(squeeze(hists2(:,1,:)),2);
K_value = mean(squeeze(hists2(:,2,:)),2);

%%%%% Step 3: Generate relative entropy and weights
for time=1:400
    RE_2_plot(time) = RE_2_func(log(R_value(time,1)),T_value(time,1),log(K_value(time,1)));
    weight_plot(time) = pi_tilde_1_func(log(R_value(time,1)),T_value(time,1),log(K_value(time,1)));

end

%%%%% Step 4: Save values for 0, 25, 50, 75, 100 years
fileID = fopen('Relative Entropy.txt','w');
fprintf(fileID,'xi_a: %.6f \n',1./theta);
fprintf(fileID,'0 yr: %.6f \n',RE_2_plot(1));
fprintf(fileID,'25 yr: %.6f \n',RE_2_plot(100));
fprintf(fileID,'50 yr: %.6f \n',RE_2_plot(200));
fprintf(fileID,'75 yr: %.6f \n',RE_2_plot(300));
fprintf(fileID,'100 yr: %.6f \n',RE_2_plot(400));
fclose(fileID);

fileID = fopen('Nordhaus Weight.txt','w');
fprintf(fileID,'xi_a: %.6f \n',1./theta);
fprintf(fileID,'0 yr: %.6f \n',weight_plot(1));
fprintf(fileID,'25 yr: %.6f \n',weight_plot(100));
fprintf(fileID,'50 yr: %.6f \n',weight_plot(200));
fprintf(fileID,'75 yr: %.6f \n',weight_plot(300));
fprintf(fileID,'100 yr: %.6f \n',weight_plot(400));
fclose(fileID);

%%%%% Step 5: Save plotting space
beta_f_space = linspace(a,b,200);
save('beta_f_space','beta_f_space')

%%%%% Step 6: Generate tilted probabilities for 0, 25, 50, 75, 100 years

% year 0
time=1;
R0 = R_value(time,1);
K0 = K_value(time,1);
F0 = T_value(time,1);

scale_2_fnc = @(x) exp(-theta.*xi_d.*(gamma_1.*x ...
    +gamma_2.*x.^2.*F0 ...
    +gamma_2_plus.*x.*(x.*F0-f_bar).^(power-1).*((x.*F0-f_bar)>=0)).*R0.*e_func(log(R0),F0,log(K0)))...
    .*normpdf(x, beta_f, sqrt(var_beta_f)); ...
scale_2 = quad_int(scale_2_fnc, A, B, n,'legendre');
q2_tilde_fnc_0 =  exp(-theta.*xi_d.*(gamma_1.*beta_f_space ...
        +gamma_2.*beta_f_space.^2.*F0 ...
        +gamma_2_plus.*beta_f_space.*(beta_f_space.*F0-f_bar).^(power-1).*((beta_f_space.*F0-f_bar)>=0)).*R0.*e_func(log(R0),F0,log(K0)))...
        ./scale_2.*normpdf(beta_f_space,beta_f,sqrt(var_beta_f));

original_dist = normpdf(beta_f_space,beta_f,sqrt(var_beta_f));
weitzman = q2_tilde_fnc_0;
original = original_dist;
% save
save('Dist_0yr','weitzman','original')

% year 25
time=100;
R0 = R_value(time,1);
K0 = K_value(time,1);
F0 = T_value(time,1);

scale_2_fnc = @(x) exp(-theta.*xi_d.*(gamma_1.*x ...
    +gamma_2.*x.^2.*F0 ...
    +gamma_2_plus.*x.*(x.*F0-f_bar).^(power-1).*((x.*F0-f_bar)>=0)).*R0.*e_func(log(R0),F0,log(K0)))...
    .*normpdf(x, beta_f, sqrt(var_beta_f)); ...
scale_2 = quad_int(scale_2_fnc, A, B, n,'legendre');
q2_tilde_fnc_0 =  exp(-theta.*xi_d.*(gamma_1.*beta_f_space ...
        +gamma_2.*beta_f_space.^2.*F0 ...
        +gamma_2_plus.*beta_f_space.*(beta_f_space.*F0-f_bar).^(power-1).*((beta_f_space.*F0-f_bar)>=0)).*R0.*e_func(log(R0),F0,log(K0)))...
        ./scale_2.*normpdf(beta_f_space,beta_f,sqrt(var_beta_f));

original_dist = normpdf(beta_f_space,beta_f,sqrt(var_beta_f));
weitzman = q2_tilde_fnc_0;
original = original_dist;
% save
save('Dist_25yr','weitzman','original')

% year 50
time=200;
R0 = R_value(time,1);
K0 = K_value(time,1);
F0 = T_value(time,1);

scale_2_fnc = @(x) exp(-theta.*xi_d.*(gamma_1.*x ...
    +gamma_2.*x.^2.*F0 ...
    +gamma_2_plus.*x.*(x.*F0-f_bar).^(power-1).*((x.*F0-f_bar)>=0)).*R0.*e_func(log(R0),F0,log(K0)))...
    .*normpdf(x, beta_f, sqrt(var_beta_f)); ...
scale_2 = quad_int(scale_2_fnc, A, B, n,'legendre');
q2_tilde_fnc_0 =  exp(-theta.*xi_d.*(gamma_1.*beta_f_space ...
        +gamma_2.*beta_f_space.^2.*F0 ...
        +gamma_2_plus.*beta_f_space.*(beta_f_space.*F0-f_bar).^(power-1).*((beta_f_space.*F0-f_bar)>=0)).*R0.*e_func(log(R0),F0,log(K0)))...
        ./scale_2.*normpdf(beta_f_space,beta_f,sqrt(var_beta_f));

original_dist = normpdf(beta_f_space,beta_f,sqrt(var_beta_f));
weitzman = q2_tilde_fnc_0;
original = original_dist;
% save
save('Dist_50yr','weitzman','original')


% year 75
time=300;

R0 = R_value(time,1);
K0 = K_value(time,1);
F0 = T_value(time,1);

scale_2_fnc = @(x) exp(-theta.*xi_d.*(gamma_1.*x ...
    +gamma_2.*x.^2.*F0 ...
    +gamma_2_plus.*x.*(x.*F0-f_bar).^(power-1).*((x.*F0-f_bar)>=0)).*R0.*e_func(log(R0),F0,log(K0)))...
    .*normpdf(x, beta_f, sqrt(var_beta_f)); ...
scale_2 = quad_int(scale_2_fnc, A, B, n,'legendre');
q2_tilde_fnc_0 =  exp(-theta.*xi_d.*(gamma_1.*beta_f_space ...
        +gamma_2.*beta_f_space.^2.*F0 ...
        +gamma_2_plus.*beta_f_space.*(beta_f_space.*F0-f_bar).^(power-1).*((beta_f_space.*F0-f_bar)>=0)).*R0.*e_func(log(R0),F0,log(K0)))...
        ./scale_2.*normpdf(beta_f_space,beta_f,sqrt(var_beta_f));

original_dist = normpdf(beta_f_space,beta_f,sqrt(var_beta_f));
weitzman = q2_tilde_fnc_0;
original = original_dist;
% save
save('Dist_75yr','weitzman','original')

% year 100
time=400;
R0 = R_value(time,1);
K0 = K_value(time,1);
F0 = T_value(time,1);

scale_2_fnc = @(x) exp(-theta.*xi_d.*(gamma_1.*x ...
    +gamma_2.*x.^2.*F0 ...
    +gamma_2_plus.*x.*(x.*F0-f_bar).^(power-1).*((x.*F0-f_bar)>=0)).*R0.*e_func(log(R0),F0,log(K0)))...
    .*normpdf(x, beta_f, sqrt(var_beta_f)); ...
scale_2 = quad_int(scale_2_fnc, A, B, n,'legendre');
q2_tilde_fnc_0 =  exp(-theta.*xi_d.*(gamma_1.*beta_f_space ...
        +gamma_2.*beta_f_space.^2.*F0 ...
        +gamma_2_plus.*beta_f_space.*(beta_f_space.*F0-f_bar).^(power-1).*((beta_f_space.*F0-f_bar)>=0)).*R0.*e_func(log(R0),F0,log(K0)))...
        ./scale_2.*normpdf(beta_f_space,beta_f,sqrt(var_beta_f));

original_dist = normpdf(beta_f_space,beta_f,sqrt(var_beta_f));
weitzman = q2_tilde_fnc_0;
original = original_dist;
% save
save('Dist_100yr','weitzman','original')



