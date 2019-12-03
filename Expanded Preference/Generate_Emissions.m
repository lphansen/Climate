clear all
% close all
clc

file1 = [pwd, '/HJB_NonLinPref_Cumu_check_exp_Sims'];
Model1 = load(file1,'e_hists2','theta','kappa','hists2');
hists2 = Model1.hists2;
e_hists2 = Model1.e_hists2;
e_values_NoUn = mean(e_hists2,2);
theta = Model1.theta;
kappa = Model1.kappa;
s1 = num2str(1./theta,4);
s1 = strrep(s1,'.','');
s2 = num2str(kappa,4);
s2 = strrep(s2,'.','');

% s_e = strcat('e_A=',s1,'_M=',s2,'.mat');
s_e = strcat('emissions_',s1,'.mat');

save(s_e,'e_values_NoUn');


%%

file1 = [pwd, '/HJB_NonLinPref_Cumu_check_exp_Sims'];
Model1 = load(file1,'hists2','e_hists2','theta','kappa');

hists2 = Model1.hists2;
e_hists2 = Model1.e_hists2;
e_values = mean(e_hists2,2);
theta = Model1.theta;
kappa = Model1.kappa;
s1 = num2str(1./theta,4);
s1 = strrep(s1,'.','');
s2 = num2str(kappa,4);
s2 = strrep(s2,'.','');

% s_e = strcat('e_A=',s1,'_M=',s2,'.mat');
s_e = strcat('emissions_',s1,'.mat');

save(s_e,'e_values');

