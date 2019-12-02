%%%%% This file solves Feyman Kac for the Ambiguity Averse model.
% Authors: Mike Barnett, Jieyao Wang
% Last update: Nov 20,2019

close all
clear all
clc

%% Step 0: Load HJB results
load([pwd,'/HJB_Growth_Averse']);

%% Step 1: set up solver
addpath('/mnt/ide0/home/wangjieyao/Climate/FT/')
addpath('/home/wangjieyao/FT/')
addpath('/Volumes/homes/FT/')
%% Step 2: Solve Feyman-Kac
[r_mat_1,t_mat_1,k_mat_1] = ndgrid(r,t,k); 

    a1 = v0_dk.*zeros(size(t_mat));
    b1 = v0_dk.*(gamma1(1)+gamma2(1).*t_bar);
    c1 = 2.*v0_dk.* gamma2(1).*(t_mat);
    
    a2 = v0_dk.*zeros(size(t_mat));
    b2 = v0_dk.*(gamma1(2)+gamma2(2).*t_bar);
    c2 = 2.*v0_dk.*gamma2(2).*(t_mat);
    
    a3 = v0_dk.*zeros(size(t_mat));
    b3 = v0_dk.*(gamma1(3)+gamma2(3).*t_bar);
    c3 = 2.*v0_dk.*gamma2(3).*(t_mat);
    
    a4 = v0_dk.*zeros(size(t_mat));
    b4 = v0_dk.*(gamma1(4)+gamma2(4).*t_bar);
    c4 = 2.*v0_dk.*gamma2(4).*(t_mat);
    
    a5 = v0_dk.*zeros(size(t_mat));
    b5 = v0_dk.*(gamma1(5)+gamma2(5).*t_bar);
    c5 = 2.*v0_dk.*gamma2(5).*(t_mat);
    
    a6 = v0_dk.*zeros(size(t_mat));
    b6 = v0_dk.*(gamma1(6)+gamma2(6).*t_bar);
    c6 = 2.*v0_dk.*gamma2(6).*(t_mat);
    
    a7 = v0_dk.*zeros(size(t_mat));
    b7 = v0_dk.*(gamma1(7)+gamma2(7).*t_bar);
    c7 = 2.*v0_dk.*gamma2(7).*(t_mat);
    
    a8 = v0_dk.*zeros(size(t_mat));
    b8 = v0_dk.*(gamma1(8)+gamma2(8).*t_bar);
    c8 = 2.*v0_dk.*gamma2(8).*(t_mat);
    
    a9 = v0_dk.*zeros(size(t_mat));
    b9 = v0_dk.*(gamma1(9)+gamma2(9).*t_bar);
    c9 = 2.*v0_dk.*gamma2(9).*(t_mat);
    
    dmg_tilt_1 = a1+b1.*beta_tilde_1+0.5.*c1.*(beta_tilde_1.^2)+0.5.*c1./lambda_tilde_1;    
    dmg_tilt_2 = a2+b2.*beta_tilde_2+0.5.*c2.*(beta_tilde_2.^2)+0.5.*c2./lambda_tilde_2;    
    dmg_tilt_3 = a3+b3.*beta_tilde_3+0.5.*c3.*(beta_tilde_3.^2)+0.5.*c3./lambda_tilde_3;    
    dmg_tilt_4 = a4+b4.*beta_tilde_4+0.5.*c4.*(beta_tilde_4.^2)+0.5.*c4./lambda_tilde_4;    
    dmg_tilt_5 = a5+b5.*beta_tilde_5+0.5.*c5.*(beta_tilde_5.^2)+0.5.*c5./lambda_tilde_5;    
    dmg_tilt_6 = a6+b6.*beta_tilde_6+0.5.*c6.*(beta_tilde_6.^2)+0.5.*c6./lambda_tilde_6;    
    dmg_tilt_7 = a7+b7.*beta_tilde_7+0.5.*c7.*(beta_tilde_7.^2)+0.5.*c7./lambda_tilde_7;    
    dmg_tilt_8 = a8+b8.*beta_tilde_8+0.5.*c8.*(beta_tilde_8.^2)+0.5.*c8./lambda_tilde_8;    
    dmg_tilt_9 = a9+b9.*beta_tilde_9+0.5.*c9.*(beta_tilde_9.^2)+0.5.*c9./lambda_tilde_9;    
        
    dmg_1 = a1+b1.*beta_f+0.5.*c1.*(beta_f.^2)+0.5.*c1./lambda;    
    dmg_2 = a2+b2.*beta_f+0.5.*c2.*(beta_f.^2)+0.5.*c2./lambda;    
    dmg_3 = a3+b3.*beta_f+0.5.*c3.*(beta_f.^2)+0.5.*c3./lambda;    
    dmg_4 = a4+b4.*beta_f+0.5.*c4.*(beta_f.^2)+0.5.*c4./lambda;    
    dmg_5 = a5+b5.*beta_f+0.5.*c5.*(beta_f.^2)+0.5.*c5./lambda;    
    dmg_6 = a6+b6.*beta_f+0.5.*c6.*(beta_f.^2)+0.5.*c6./lambda;    
    dmg_7 = a7+b7.*beta_f+0.5.*c7.*(beta_f.^2)+0.5.*c7./lambda;    
    dmg_8 = a8+b8.*beta_f+0.5.*c8.*(beta_f.^2)+0.5.*c8./lambda;    
    dmg_9 = a9+b9.*beta_f+0.5.*c9.*(beta_f.^2)+0.5.*c9./lambda;    
    

    flow_base =  weight1.*dmg_1+ weight2.*dmg_2...
          + weight3.*dmg_3+ weight4.*dmg_4...
          + weight5.*dmg_5+ weight6.*dmg_6...
          + weight7.*dmg_7+ weight8.*dmg_8...
          + weight9.*dmg_9;

    flow_tilted =    pi_tilde_1_norm.*dmg_tilt_1+pi_tilde_2_norm.*dmg_tilt_2...
         +pi_tilde_3_norm.*dmg_tilt_3+pi_tilde_4_norm.*dmg_tilt_4...
         +pi_tilde_5_norm.*dmg_tilt_5+pi_tilde_6_norm.*dmg_tilt_6...
         +pi_tilde_7_norm.*dmg_tilt_7+pi_tilde_8_norm.*dmg_tilt_8...
         +pi_tilde_9_norm.*dmg_tilt_9;
% % % % % % % % % % % % % % % % % % % % % % % 
    a1 = v0_dk.*(gamma0(1)+gamma1(1).*t_bar+0.5.*gamma2(1).*(t_bar.^2));
    b1 = v0_dk.*t_mat.*(gamma1(1)+gamma2(1).*t_bar);
    c1 = v0_dk.*gamma2(1).*(t_mat.^2);
    
    a2 = v0_dk.*(gamma0(2)+gamma1(2).*t_bar+0.5.*gamma2(2).*(t_bar.^2));
    b2 = v0_dk.*t_mat.*(gamma1(2)+gamma2(2).*t_bar);
    c2 = v0_dk.*gamma2(2).*(t_mat.^2);
    
    a3 = v0_dk.*(gamma0(3)+gamma1(3).*t_bar+0.5.*gamma2(3).*(t_bar.^2));
    b3 = v0_dk.*t_mat.*(gamma1(3)+gamma2(3).*t_bar);
    c3 = v0_dk.*gamma2(3).*(t_mat.^2);
    
    a4 = v0_dk.*(gamma0(4)+gamma1(4).*t_bar+0.5.*gamma2(4).*(t_bar.^2));
    b4 = v0_dk.*t_mat.*(gamma1(4)+gamma2(4).*t_bar);
    c4 = v0_dk.*gamma2(4).*(t_mat.^2);
    
    a5 = v0_dk.*(gamma0(5)+gamma1(5).*t_bar+0.5.*gamma2(5).*(t_bar.^2));
    b5 = v0_dk.*t_mat.*(gamma1(5)+gamma2(5).*t_bar);
    c5 = v0_dk.*gamma2(5).*(t_mat.^2);
    
    a6 = v0_dk.*(gamma0(6)+gamma1(6).*t_bar+0.5.*gamma2(6).*(t_bar.^2));
    b6 = v0_dk.*t_mat.*(gamma1(6)+gamma2(6).*t_bar);
    c6 = v0_dk.*gamma2(6).*(t_mat.^2);
    
    a7 = v0_dk.*(gamma0(7)+gamma1(7).*t_bar+0.5.*gamma2(7).*(t_bar.^2));
    b7 = v0_dk.*t_mat.*(gamma1(7)+gamma2(7).*t_bar);
    c7 = v0_dk.*gamma2(7).*(t_mat.^2);
    
    a8 = v0_dk.*(gamma0(8)+gamma1(8).*t_bar+0.5.*gamma2(8).*(t_bar.^2));
    b8 = v0_dk.*t_mat.*(gamma1(8)+gamma2(8).*t_bar);
    c8 = v0_dk.*gamma2(8).*(t_mat.^2);
    
    a9 = v0_dk.*(gamma0(9)+gamma1(9).*t_bar+0.5.*gamma2(9).*(t_bar.^2));
    b9 = v0_dk.*t_mat.*(gamma1(9)+gamma2(9).*t_bar);
    c9 = v0_dk.*gamma2(9).*(t_mat.^2);
    
    dmg_tilt_1 = a1+b1.*beta_tilde_1+0.5.*c1.*(beta_tilde_1.^2)+0.5.*c1./lambda_tilde_1;    
    dmg_tilt_2 = a2+b2.*beta_tilde_2+0.5.*c2.*(beta_tilde_2.^2)+0.5.*c2./lambda_tilde_2;    
    dmg_tilt_3 = a3+b3.*beta_tilde_3+0.5.*c3.*(beta_tilde_3.^2)+0.5.*c3./lambda_tilde_3;    
    dmg_tilt_4 = a4+b4.*beta_tilde_4+0.5.*c4.*(beta_tilde_4.^2)+0.5.*c4./lambda_tilde_4;    
    dmg_tilt_5 = a5+b5.*beta_tilde_5+0.5.*c5.*(beta_tilde_5.^2)+0.5.*c5./lambda_tilde_5;    
    dmg_tilt_6 = a6+b6.*beta_tilde_6+0.5.*c6.*(beta_tilde_6.^2)+0.5.*c6./lambda_tilde_6;    
    dmg_tilt_7 = a7+b7.*beta_tilde_7+0.5.*c7.*(beta_tilde_7.^2)+0.5.*c7./lambda_tilde_7;    
    dmg_tilt_8 = a8+b8.*beta_tilde_8+0.5.*c8.*(beta_tilde_8.^2)+0.5.*c8./lambda_tilde_8;    
    dmg_tilt_9 = a9+b9.*beta_tilde_9+0.5.*c9.*(beta_tilde_9.^2)+0.5.*c9./lambda_tilde_9;    
        
    dmg_1 = a1+b1.*beta_f+0.5.*c1.*(beta_f.^2)+0.5.*c1./lambda;    
    dmg_2 = a2+b2.*beta_f+0.5.*c2.*(beta_f.^2)+0.5.*c2./lambda;    
    dmg_3 = a3+b3.*beta_f+0.5.*c3.*(beta_f.^2)+0.5.*c3./lambda;    
    dmg_4 = a4+b4.*beta_f+0.5.*c4.*(beta_f.^2)+0.5.*c4./lambda;    
    dmg_5 = a5+b5.*beta_f+0.5.*c5.*(beta_f.^2)+0.5.*c5./lambda;    
    dmg_6 = a6+b6.*beta_f+0.5.*c6.*(beta_f.^2)+0.5.*c6./lambda;    
    dmg_7 = a7+b7.*beta_f+0.5.*c7.*(beta_f.^2)+0.5.*c7./lambda;    
    dmg_8 = a8+b8.*beta_f+0.5.*c8.*(beta_f.^2)+0.5.*c8./lambda;    
    dmg_9 = a9+b9.*beta_f+0.5.*c9.*(beta_f.^2)+0.5.*c9./lambda;    
    

    Gamma_base =  weight1.*dmg_1+ weight2.*dmg_2...
          + weight3.*dmg_3+ weight4.*dmg_4...
          + weight5.*dmg_5+ weight6.*dmg_6...
          + weight7.*dmg_7+ weight8.*dmg_8...
          + weight9.*dmg_9;

    Gamma_tilted =    pi_tilde_1_norm.*dmg_tilt_1+pi_tilde_2_norm.*dmg_tilt_2...
             +pi_tilde_3_norm.*dmg_tilt_3+pi_tilde_4_norm.*dmg_tilt_4...
             +pi_tilde_5_norm.*dmg_tilt_5+pi_tilde_6_norm.*dmg_tilt_6...
             +pi_tilde_7_norm.*dmg_tilt_7+pi_tilde_8_norm.*dmg_tilt_8...
             +pi_tilde_9_norm.*dmg_tilt_9;


    v0 = (alpha).*r_mat_1+(1-alpha).*k_mat_1;
    v1_initial = v0.*ones(size(r_mat_1)); 

    %%%%%%%%% baseline model %%%%%%%%%%
    A = -delta.*ones(size(r_mat_1));
    B_r = -e_star+Gamma_r.*(f.^Theta_r).*exp(Theta_r.*(k_mat - r_mat))-0.5.*(sigma_r.^2);
    B_k = Alpha+Gamma.*log(1+i_k./Theta)-0.5.*(sigma_k.^2)-Gamma_base;
    B_t = e_star.*exp(r_mat_1);
    C_rr = 0.5.*sigma_r.^2.*ones(size(r_mat_1));
    C_kk = 0.5.*sigma_k.^2.*ones(size(r_mat_1));
    C_tt = zeros(size(r_mat));
    D = flow_base;
    %%%%%%%% PDE inputs %%%%%%%%%%%
    stateSpace = [r_mat_1(:), t_mat_1(:), k_mat_1(:)]; 
    model      = {};
    model.A    = A(:); 
    model.B    = [B_r(:), B_t(:), B_k(:)];
    model.C    = [C_rr(:), C_tt(:), C_kk(:)];
    model.D    = D(:);
    model.v0   = v0(:).*0;
    model.dt   = 1.0;
    %%% solve system of equations
    out = solveCGNatural_1_1e9(stateSpace, model);
    out_comp = reshape(out,size(v0)).*ones(size(r_mat_1));

    disp(['Error: ', num2str(max(max(max(max(abs(out_comp - v1_initial))))))])
    iter = 1;
    
    v0 = v0.*ones(size(v0));
    v0 = reshape(out,size(v0));
 
filename2 = [pwd, '/SCC_mat_Cumu_base_GrowthAmb_1e9'];

save(filename2)

