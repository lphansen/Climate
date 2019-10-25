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


% % baseline
% base_coarse = load("C:\\Users\\hanxuh\\Google Drive\\MFR\\github\\Climate\\Preference\\Weighted\\HJB_NonLinPref_Cumu_Sims_F4000_R90_K90.mat")
% base_fine = load("C:\\Users\\hanxuh\\Google Drive\\MFR\\github\\Climate\\Preference\\Weighted\\HJB_NonLinPref_Cumu_Sims_F4000_R90_K90_D150.mat")
% 
% % 50% more length
% long_coarse = load("C:\\Users\\hanxuh\\Google Drive\\MFR\\github\\Climate\\Preference\\Weighted\\HJB_NonLinPref_Cumu_Sims_F6000_R135_K135.mat")
% long_fine = load("C:\\Users\\hanxuh\\Google Drive\\MFR\\github\\Climate\\Preference\\Weighted\\HJB_NonLinPref_Cumu_Sims_F6000_R135_K135_D150.mat")

% figure()
% plot(base_coarse.j_hists2)
% hold on
% plot(base_fine.j_hists2)
% hold off

% base_coarse = load("C:\\Users\\hanxuh\\Google Drive\\MFR\\github\\Climate\\Preference\\Weighted\\HJB_NonLinPref_Cumu_Sims_F6000_R135_K135_D150.mat")
% base_fine = load("C:\\Users\\hanxuh\\Google Drive\\MFR\\github\\Climate\\Preference\\Weighted\\HJB_NonLinPref_Cumu_Sims_F6000_R135_K135_D200.mat")
% 
% figure()
% plot(base_coarse.v_dr_hists2)
% hold on
% plot(base_fine.v_dr_hists2)
% hold off
% 
% legend('v_r (150% grid density)','v_r (200% grid density)')
% 
% figure()
% plot(base_coarse.j_hists2)
% hold on
% plot(base_fine.j_hists2)
% hold off
% 
% legend('j (150% grid density)','j (200% grid density)')
% 
% figure()
% plot(base_coarse.j_hists2.^psi_1)
% hold on
% plot(base_fine.j_hists2.^psi_1)
% hold off
% 
% legend('j''contribution to drift (150% grid density)','j''contribution to drift (200% grid density)')
% 
% figure()
% plot(log(base_coarse.hists2(:,1)))
% hold on
% plot(log(base_fine.hists2(:,1)))
% hold off
% 
% legend('logR (150% grid density)','logR (200% grid density)')
% 
% figure()
% plot(base_coarse.e_hists2)
% hold on
% plot(base_fine.e_hists2)
% hold off
% 
% legend('emission (150% grid density)','emission (200% grid density)')

% andrew = load("HJB_NonLinPref_Cumu_Sims.mat");
% base = load("HJB_NonLinPref_Cumu_Sims_F4000_R90_K90");
% 
% figure()
% plot(andrew.j_hists2)
% hold on
% plot(base.j_hists2)
% hold off


% baseline = load("C:\\Users\\hanxuh\\Google Drive\\MFR\\github\\Climate\\Preference\\Weighted\\HJB_NonLinPref_Cumu_Sims_F6000_R135_K135_D200.mat")
% fine_F = load("C:\\Users\\hanxuh\\Google Drive\\MFR\\github\\Climate\\Preference\\Weighted\\HJB_NonLinPref_Cumu_Sims_F6000_R135_K135_4xF_2xRK.mat")
% 
% 
% figure()
% plot(baseline.v_dr_hists2)
% hold on
% plot(fine_F.v_dr_hists2)
% hold off
% 
% legend('v_r (2x F,R,K)','v_r (4x F, 2x R,K)')
% 
% figure()
% plot(baseline.j_hists2)
% hold on
% plot(fine_F.j_hists2)
% hold off
% 
% legend('j (2x F,R,K)','j (4x F, 2x R,K)')
% 
% figure()
% plot(baseline.j_hists2.^baseline.psi_1)
% hold on
% plot(fine_F.j_hists2.^baseline.psi_1)
% hold off
% 
% legend('j''contribution to drift (2x F,R,K)','j''contribution to drift (4x F, 2x R,K)')
% 
% figure()
% plot(log(baseline.hists2(:,1)))
% hold on
% plot(log(fine_F.hists2(:,1)))
% hold off
% 
% legend('logR (2x F,R,K)','logR (4x F, 2x R,K)')
% 
% figure()
% plot(baseline.e_hists2)
% hold on
% plot(fine_F.e_hists2)
% hold off

% legend('emission (2x F,R,K)','emission (4x F, 2x R,K)')

baseline = load("C:\\Users\\hanxuh\\Google Drive\\MFR\\github\\Climate\\Preference\\Weighted\\HJB_NonLinPref_Cumu_Sims_F6000_R135_K135.mat")
fine_F = load("C:\\Users\\hanxuh\\Google Drive\\MFR\\github\\Climate\\Preference\\Weighted\\HJB_NonLinPref_Cumu_Sims_F6000_R135_K135_4xRF_1xK_VG.mat")


figure()
plot(baseline.v_dr_hists2)
hold on
plot(fine_F.v_dr_hists2)
hold off

legend('v_r (baseline)','v_r (4x FR, 1x K)')

figure()
plot(baseline.j_hists2)
hold on
plot(fine_F.j_hists2)
hold off

legend('j (baseline)','j (4x FR, 1x K)')

figure()
plot(baseline.j_hists2.^baseline.psi_1)
hold on
plot(fine_F.j_hists2.^baseline.psi_1)
hold off

legend('j''contribution to drift (baseline)','j''contribution to drift (4x FR, 1x K)')

figure()
plot(log(baseline.hists2(:,1)))
hold on
plot(log(fine_F.hists2(:,1)))
hold off

legend('logR (baseline)','logR (4x FR, 1x K)')

figure()
plot(baseline.e_hists2)
hold on
plot(fine_F.e_hists2)
hold off

legend('emission (baseline)','emission (4x FR, 1x K)')