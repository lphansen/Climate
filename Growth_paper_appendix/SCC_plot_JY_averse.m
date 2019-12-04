%%%%% This file generates figures for Table 4 and Table G.3 for the Ambiguity Averse model.
% Authors: Mike Barnett, Jieyao Wang
% Last update: Nov 20,2019

%% result for Table 4
close all;
clc;
clear all;

file2 = [pwd, '/SCC_mat_Cumu_base_GrowthAmb_1e9'];

Model2 = load(file2,'v0','r_mat','k_mat','t_mat','i_k','f','e','v0_dk','v0_dr',...
             'A_O','alpha','delta','Theta','Gamma','gamma_1','gamma_2','t_bar');

external_v0 = Model2.v0;
r_mat_1 = Model2.r_mat;
k_mat_1 = Model2.k_mat;
t_mat_1 = Model2.t_mat;
gamma_1 = Model2.gamma_1;
gamma_2 = Model2.gamma_2;
t_bar = Model2.t_bar;
A_O = Model2.A_O;
alpha = Model2.alpha;
delta = Model2.delta;
Gamma = Model2.Gamma;
Theta = Model2.Theta;
e_1 = Model2.e;
f_1 = Model2.f;
i_k_1 = Model2.i_k;

file2 = [pwd, '/HJB_Growth_Averse'];
Model2 = load(file2,'v0_dr','v0_dt');
v0_dr_1 = Model2.v0_dr;
v0_dt_1 = Model2.v0_dt;

MC = delta.*(1-alpha)./(A_O.*exp(k_mat_1)-i_k_1.*exp(k_mat_1)-f_1.*exp(k_mat_1));
ME = delta.*alpha./(e_1.*exp(r_mat_1));
SCC = 1000*ME./MC;
SCC_func = griddedInterpolant(r_mat_1,t_mat_1,k_mat_1,SCC,'linear');

ME1 = (v0_dr_1.*exp(-r_mat_1));
SCC1 = 1000*ME1./MC;
SCC1_func = griddedInterpolant(r_mat_1,t_mat_1,k_mat_1,SCC1,'linear');

ME2_base =  external_v0;
SCC2_base = 1000*ME2_base./MC;
SCC2_base_func = griddedInterpolant(r_mat_1,t_mat_1,k_mat_1,SCC2_base,'linear');

ME2_base_a = -v0_dt_1;
SCC2_base_a = 1000*ME2_base_a./MC;
SCC2_base_a_func = griddedInterpolant(r_mat_1,t_mat_1,k_mat_1,SCC2_base_a,'linear');


file2 = [pwd, '/SCC_mat_Cumu_worst_GrowthAmb_1e9'];
Model2 = load(file2,'v0');
external_v0_worst = Model2.v0;

ME2_tilt =  external_v0_worst;
SCC2_tilt = 1000*ME2_tilt./MC;
SCC2_tilt_func = griddedInterpolant(r_mat_1,t_mat_1,k_mat_1,SCC2_tilt,'linear');

file1 = [pwd, '/HJB_Growth_Averse_Sims_JY'];
Model1 = load(file1,'hists2','e_hists2','theta','RE_hists2');

hists2_A = Model1.hists2;
RE_hist = Model1.RE_hists2;
e_hist = Model1.e_hists2;

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
SCC2_FK_tilt = mean(SCC2_tilt_values,2);
external = mean(SCC2_FK_base,2);
uncert_diff_in_FK = mean(SCC2_FK_tilt-SCC2_FK_base,2);

% generate results of Table 4 in the paper
fileID = fopen('Growth_numbers_Amb.txt','w');
fprintf(fileID,'total_0yr: %.6f \n',SCC_total(1));
fprintf(fileID,'total_50yr: %.6f \n',SCC_total(end./2));
fprintf(fileID,'total_100yr: %.6f \n',SCC_total(end));
fprintf(fileID,'uncertainty_0yr_diff_in_FK: %.6f \n',uncert_diff_in_FK(1));
fprintf(fileID,'uncertainty_50yr_diff_in_FK: %.6f \n',uncert_diff_in_FK(end./2));
fprintf(fileID,'uncertainty_100yr_diff_in_FK: %.6f \n',uncert_diff_in_FK(end));
fprintf(fileID,'Emissions_0yr: %.6f \n',e_hist(1));
fprintf(fileID,'Emissions_50yr: %.6f \n',e_hist(end./2));
fprintf(fileID,'Emissions_100yr: %.6f \n',e_hist(end));
fprintf(fileID,'Relative Entropy_0yr: %.6f \n',RE_hist(1));
fprintf(fileID,'Relative Entropy_50yr: %.6f \n',RE_hist(end./2));
fprintf(fileID,'Relative Entropy_100yr: %.6f \n',RE_hist(end));
fclose(fileID);


%% Result for Table G.3 from the Online Appendix
for time=1:400
    for path=1:1
    SCC_values(time,path) = SCC_func(log((hists2_A(time,1,path))),(hists2_A(time,3,path)),log((hists2_A(time,4,path))));
    SCC1_values(time,path) = SCC1_func(log((hists2_A(time,1,path))),(hists2_A(time,3,path)),log((hists2_A(time,4,path))));
    SCC2_base_values(time,path) = SCC2_base_func(log((hists2_A(time,1,path))),(hists2_A(time,3,path)),log((hists2_A(time,4,path))));
    SCC2_base_a_values(time,path) = SCC2_base_a_func(log((hists2_A(time,1,path))),(hists2_A(time,3,path)),log((hists2_A(time,4,path))));
    SCC2_tilt_values(time,path) = SCC2_tilt_func(log((hists2_A(time,1,path))),(hists2_A(time,3,path)),log((hists2_A(time,4,path))));

    end
end

SCC_total = mean(SCC_values,2);
SCC_private = mean(SCC1_values,2);
SCC2_FK_base = mean(SCC2_base_values,2);
SCC2_FK_tilt = mean(SCC2_tilt_values,2);
external = mean(SCC2_FK_base,2);
uncert_diff_in_FK = mean(SCC2_FK_tilt-SCC2_FK_base,2);

% generate results of Table G.3 in the paper
fileID = fopen('Growth_numbers_Amb_Numeraire.txt','w');
fprintf(fileID,'total_0yr: %.6f \n',SCC_total(1));
fprintf(fileID,'total_50yr: %.6f \n',SCC_total(end./2));
fprintf(fileID,'total_100yr: %.6f \n',SCC_total(end));
fprintf(fileID,'uncertainty_0yr_diff_in_FK: %.6f \n',uncert_diff_in_FK(1));
fprintf(fileID,'uncertainty_50yr_diff_in_FK: %.6f \n',uncert_diff_in_FK(end./2));
fprintf(fileID,'uncertainty_100yr_diff_in_FK: %.6f \n',uncert_diff_in_FK(end));
fprintf(fileID,'Emissions_0yr: %.6f \n',e_hist(1));
fprintf(fileID,'Emissions_50yr: %.6f \n',e_hist(end./2));
fprintf(fileID,'Emissions_100yr: %.6f \n',e_hist(end));
fprintf(fileID,'Relative Entropy_0yr: %.6f \n',RE_hist(1));
fprintf(fileID,'Relative Entropy_50yr: %.6f \n',RE_hist(end./2));
fprintf(fileID,'Relative Entropy_100yr: %.6f \n',RE_hist(end));
fclose(fileID);



