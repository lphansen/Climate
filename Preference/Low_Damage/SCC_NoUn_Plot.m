%%%%%%%%%%%%%%%%%%%%%%%%% Social Cost of Carbon %%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear all;
close all;

%%%%% Step 0: Load simulated trajectory
file1 = [pwd, '/HJB_NonLinPref_Cumu_NoUn_Sims'];
Model1 = load(file1,'hists2','i_k_hists2','e_hists2','f_hists2','A_O','alpha','theta','kappa');
data1 = Model1.hists2;
i_k = Model1.i_k_hists2;
e = Model1.e_hists2;
f = Model1.f_hists2;
A_O = Model1.A_O;
alpha = Model1.alpha;
theta = Model1.theta;
kappa = Model1.kappa;

%%%%% Step 1: Compute SCC
SCC = 1000*mean(((alpha./(1-alpha)).*(A_O.*squeeze(data1(:,2,:))-i_k-f)./e),2);

%%%%% Step 2: Save file
s1 = num2str(1./theta,4);
s1 = strrep(s1,'.','');
s2 = num2str(kappa,4);
s2 = strrep(s2,'.','');
s = strcat('/HJB_NonLinPref_A=',s1,'_M=',s2);
s0 = '/HJB_NonLinPref';
filename2 = [pwd, s0];
s_scc = strcat('SCC_A=',s1,'_M=',s2,'.mat');
% s_scc = strcat('SCC_',s1,'.mat');
save(s_scc,'SCC');

