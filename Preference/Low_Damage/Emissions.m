%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Emissions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc

%%%%% Step 0: Load simulated trajectory for no uncertainty
file1 = [pwd, '/HJB_NonLinPref_Cumu_NoUn_Sims'];
Model1 = load(file1,'e_hists2','theta','kappa','hists2');
e_hists2 = Model1.e_hists2;
e_values_NoUn = mean(e_hists2,2);
theta = Model1.theta;
kappa = Model1.kappa;

%%%%% Step 1: Save file
s1 = num2str(1./theta,4);
s1 = strrep(s1,'.','');
s2 = num2str(kappa,4);
s2 = strrep(s2,'.','');
s_e = strcat('e_A=',s1,'_M=',s2,'.mat');
save(s_e,'e_values_NoUn');


%%%%% Step 2: Repeat for with uncertainty simulated trajectory

%%%%% Load simulated trajectory for no uncertainty
file1 = [pwd, '/HJB_NonLinPref_Cumu_Sims'];
Model1 = load(file1,'hists2','e_hists2','theta','kappa');
e_hists2 = Model1.e_hists2;
e_values = mean(e_hists2,2);
theta = Model1.theta;
kappa = Model1.kappa;

%%%%% Save file
s1 = num2str(1./theta,4);
s1 = strrep(s1,'.','');
s2 = num2str(kappa,4);
s2 = strrep(s2,'.','');
s_e = strcat('e_A=',s1,'_M=',s2,'.mat');
save(s_e,'e_values');

