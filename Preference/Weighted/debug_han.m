% a_dic = load("HJB_NonLinPref_Cumu_Sims.mat")
% a = a_dic.j_hists2
% plot(a)
% hold on
% b_dic = load("HJB_NonLinPref_Cumu_old_Sims.mat")
% b = b_dic.j_hists2
% plot(b)
% hold off
% hold on

% figure()
% c_dic = load("HJB_NonLinPref_Cumu_Sims.mat")
% c = c_dic.v_dr_hists2
% plot(c)
% hold on
% d_dic = load("HJB_NonLinPref_Cumu_old_Sims.mat")
% d = d_dic.v_dr_hists2
% plot(d)
% hold off
% a_dic = load("HJB_NonLinPref_Cumu_Sims.mat")
% hists2 = a_dic.hists2
% plot(hists2(:,1))
% hold on
% b_dic = load("HJB_NonLinPref_Cumu_old_Sims.mat")
% hists2 = b_dic.hists2
% plot(hists2(:,1))
% hold off
% hold on

% figure()
% c_dic = load("HJB_NonLinPref_Cumu.mat")
% plot(c_dic.v0_dr(20,:,18))
% 
% hold on
% d_dic = load("HJB_NonLinPref_Cumu_old.mat")
% plot(d_dic.v0_dr(20,:,18))
% 
% legend('new','old')
% 
% hold off
% 
% a_result = load("HJB_NonLinPref_Cumu.mat")
% b_result = load("HJB_NonLinPref_Cumu_old.mat")
% 
% max(abs(a_result.j - b_result.j),[],"all")
% min(abs(a_result.j),[],"all")
% 
% max(abs(a_result.v0_dk-b_result.v0_dk),[],"all")
% min(a_result.v0_dk,[],"all")
% 
% max(abs(a_result.v0_dr-b_result.v0_dr),[],"all")
% min(a_result.v0_dr,[],"all")


% max(abs(a_result.v0 - b_result.v0),[],'all')
% sum(abs(a_result.v0 - b_result.v0),'all')


% j_result = load("HJB_NonLinPref_Cumu_jieyao.mat")
% m_result = load("HJB_NonLinPref_Cumu_mike.mat")

% baseline
base_coarse = load("C:\\Users\\hanxuh\\Google Drive\\MFR\\github\\Climate\\Preference\\Weighted\\HJB_NonLinPref_Cumu_Sims_F4000_R90_K90.mat")
base_fine = load("C:\\Users\\hanxuh\\Google Drive\\MFR\\github\\Climate\\Preference\\Weighted\\HJB_NonLinPref_Cumu_Sims_F4000_R90_K90_D150.mat")

% 50% more length
long_coarse = load("C:\\Users\\hanxuh\\Google Drive\\MFR\\github\\Climate\\Preference\\Weighted\\HJB_NonLinPref_Cumu_Sims_F6000_R135_K135.mat")
long_fine = load("C:\\Users\\hanxuh\\Google Drive\\MFR\\github\\Climate\\Preference\\Weighted\\HJB_NonLinPref_Cumu_Sims_F6000_R135_K135_D150.mat")

% figure()
% plot(base_coarse.j_hists2)
% hold on
% plot(base_fine.j_hists2)
% hold off


figure()
plot(long_coarse.j_hists2)
hold on
plot(long_fine.j_hists2)
hold off

