%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parameter file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc

%% Uncertainty Parameters
theta = 0.001;
kappa = 1000;
if kappa>10
    RobDrift = 0;
else
    RobDrift = 1./(2.*kappa);
end

%% Emission Parameters

McD = csvread('TCRE_MacDougallEtAl2017_update.csv');
McD = McD./1000;
lambda_McD = mean(McD(:,1));
var_lambda_McD = var(McD(:,1));
lambda_max_McD = max(McD(:,1));
skew_lambda_McD = skewness(McD(:,1)).*std(McD(:,1)).^3;
kur_lambda_McD = kurtosis(McD(:,1)).*std(McD(:,1)).^4;
par_lambda_McD = McD(:,1);

%% Economic/capital
delta = 0.01;
alphaO= 0.032;
xi_m = alphaO;
xi_m_3state = 1;
xi_o = alphaO;
xi_g = 1-alphaO;
sigma_g = 0.02;
sigma_k = 0.0161;
sigma_r = 0.0339;
sigma_m = 0.01;
sigma_n = 0.01;
A_O = 0.115000000000000;
A_G = 0.1;
alpha = alphaO;
Gamma = 0.0600;
Theta = 1./16.666666666666668;
Alpha = -0.034977443912449;
Gamma_r = 0.112733407891680;
Theta_r = 0.142857142857143;
xi_k = 1-alpha;

%% save filename
s0 = '/climate_pars';
filename = [pwd, s0];
save(filename);

