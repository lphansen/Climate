%%%%% This file generates results for the ambiguity adjusted probabilities of Consumption Damage model.
% Authors: Mike Barnett, Jieyao Wang
% Last update: Dec 9, 2019
close all
clear all
clc

%% Step 1: Specify ambiguity and damage level

% Ambiguity setting: 'averse' or 'neutral'
ambiguity = 'averse';
% Damage setting: 'high', 'low', or 'weighted'
damage_level = 'weighted';

if strcmp(ambiguity,'averse')
    xi_p = 1 ./ 4000; % ambiguity parameter
elseif strcmp(ambiguity,'neutral')
    xi_p = 1000; % ambiguity parameter
else
    disp('Error: please choose between ''averse'' and ''weighted'' for ambiguity level');
end

if strcmp(damage_level,'high')
    weight = 0.0; % weight on nordhaus
elseif strcmp(damage_level,'low')
    weight = 1.0; % weight on nordhaus
elseif strcmp(damage_level,'weighted')
    weight = 0.5; % weight on nordhaus
else
    disp('Error: please choose from ''high'',''weighted'' and ''low'' for damage level')
end


%% Step 2: Calculate probabilities
% Case 1: High damage
if strcmp(damage_level,'high')
    % load results
    file1 = [ambiguity,'_',damage_level];
    Model1 = load(file1,'R_1','pi_tilde_1_norm','r_mat','k_mat','t_mat',...
        'beta_tilde_1','beta_f','lambda_tilde_1','theta','var_beta_f'); 

    RE_1 = Model1.R_1;
    pi_tilde_1 = Model1.pi_tilde_1_norm;
    beta_tilde_1 = Model1.beta_tilde_1;
    beta_f = Model1.beta_f;
    lambda_tilde_1 = Model1.lambda_tilde_1;
    theta = Model1.theta;
    var_beta_f = Model1.var_beta_f;
    r_mat = Model1.r_mat;
    t_mat = Model1.t_mat;
    k_mat = Model1.k_mat;
    a = beta_f-10.*sqrt(var_beta_f);
    b = beta_f+10.*sqrt(var_beta_f);

    RE_1_func = griddedInterpolant(r_mat,t_mat,k_mat,RE_1,'linear');
    pi_tilde_1_func = griddedInterpolant(r_mat,t_mat,k_mat,pi_tilde_1,'linear');
    beta_tilde_1_func = griddedInterpolant(r_mat,t_mat,k_mat,beta_tilde_1,'linear');
    lambda_tilde_1_func = griddedInterpolant(r_mat,t_mat,k_mat,lambda_tilde_1,'linear');

    file1 = [ambiguity,'_',damage_level,'_sim'];
    Model1 = load(file1,'hists2'); 
    hists2 = Model1.hists2;
    T_value = mean(squeeze(hists2(:,3,:)),2);
    R_value = mean(squeeze(hists2(:,1,:)),2);
    K_value = mean(squeeze(hists2(:,2,:)),2);

    % relative entropy and weights
    for time=1:400
        RE_1_plot(time) = RE_1_func(log(R_value(time,1)),T_value(time,1),log(K_value(time,1)));
        weight_plot(time) = pi_tilde_1_func(log(R_value(time,1)),T_value(time,1),log(K_value(time,1)));
        nordhaus_mean(time) = beta_tilde_1_func(log(R_value(time,1)),T_value(time,1),log(K_value(time,1)));
        sd_nordhaus(time) = 1./sqrt(lambda_tilde_1_func(log(R_value(time,1)),T_value(time,1),log(K_value(time,1))));
    end

    fileID = fopen([ambiguity,'_',damage_level,'_Relative Entropy.txt'],'w');
    fprintf(fileID,'xi_a: %.6f \n',1./theta);
    fprintf(fileID,'0 yr: %.6f \n',RE_1_plot(1));
    fprintf(fileID,'25 yr: %.6f \n',RE_1_plot(100));
    fprintf(fileID,'50 yr: %.6f \n',RE_1_plot(200));
    fprintf(fileID,'75 yr: %.6f \n',RE_1_plot(300));
    fprintf(fileID,'100 yr: %.6f \n',RE_1_plot(400));
    fclose(fileID);

    fileID = fopen([ambiguity,'_',damage_level,'_Shifted Mean.txt'],'w');
    fprintf(fileID,'xi_a: %.6f \n',1./theta);
    fprintf(fileID,'0 yr: %.6f \n',nordhaus_mean(1));
    fprintf(fileID,'25 yr: %.6f \n',nordhaus_mean(100));
    fprintf(fileID,'50 yr: %.6f \n',nordhaus_mean(200));
    fprintf(fileID,'75 yr: %.6f \n',nordhaus_mean(300));
    fprintf(fileID,'100 yr: %.6f \n',nordhaus_mean(400));
    fclose(fileID);

    fileID = fopen([ambiguity,'_',damage_level,'_Shifted Standard Deviation.txt'],'w');
    fprintf(fileID,'xi_a: %.6f \n',1./theta);
    fprintf(fileID,'0 yr: %.6f \n',sd_nordhaus(1));
    fprintf(fileID,'25 yr: %.6f \n',sd_nordhaus(100));
    fprintf(fileID,'50 yr: %.6f \n',sd_nordhaus(200));
    fprintf(fileID,'75 yr: %.6f \n',sd_nordhaus(300));
    fprintf(fileID,'100 yr: %.6f \n',sd_nordhaus(400));
    fclose(fileID);

    fileID = fopen([ambiguity,'_',damage_level,'_Nordhaus Weight.txt'],'w');
    fprintf(fileID,'xi_a: %.6f \n',1./theta);
    fprintf(fileID,'0 yr: %.6f \n',weight_plot(1));
    fprintf(fileID,'25 yr: %.6f \n',weight_plot(100));
    fprintf(fileID,'50 yr: %.6f \n',weight_plot(200));
    fprintf(fileID,'75 yr: %.6f \n',weight_plot(300));
    fprintf(fileID,'100 yr: %.6f \n',weight_plot(400));
    fclose(fileID);

    beta_f_space = linspace(a,b,200); % space for beta_f
    save([ambiguity,'_',damage_level,'_beta_f_space'],'beta_f_space')

    %%% probabilitiess

    % year 0
    time=1;
    R0 = R_value(time,1);
    K0 = K_value(time,1);
    F0 = T_value(time,1);

    mean_distort_nordhaus = beta_tilde_1_func(log(R0),F0,log(K0))-beta_f;
    lambda_tilde_nordhaus = lambda_tilde_1_func(log(R0),F0,log(K0));
    original_dist = normpdf(beta_f_space,beta_f,sqrt(var_beta_f));
    nordhaus_dist_0 = normpdf(beta_f_space,mean_distort_nordhaus+beta_f,1./sqrt(lambda_tilde_nordhaus));
    nordhaus = nordhaus_dist_0;
    original = original_dist;
    save([ambiguity,'_',damage_level,'_Dist_0yr'],'nordhaus','original')

    % year 25
    time=100;
    R0 = R_value(time,1);
    K0 = K_value(time,1);
    F0 = T_value(time,1);

    mean_distort_nordhaus = beta_tilde_1_func(log(R0),F0,log(K0))-beta_f;
    lambda_tilde_nordhaus = lambda_tilde_1_func(log(R0),F0,log(K0));
    original_dist = normpdf(beta_f_space,beta_f,sqrt(var_beta_f));
    nordhaus_dist_0 = normpdf(beta_f_space,mean_distort_nordhaus+beta_f,1./sqrt(lambda_tilde_nordhaus));
    nordhaus = nordhaus_dist_0;
    original = original_dist;
    save([ambiguity,'_',damage_level,'_Dist_25yr'],'nordhaus','original')

    % year 50
    time=200;
    R0 = R_value(time,1);
    K0 = K_value(time,1);
    F0 = T_value(time,1);

    mean_distort_nordhaus = beta_tilde_1_func(log(R0),F0,log(K0))-beta_f;
    lambda_tilde_nordhaus = lambda_tilde_1_func(log(R0),F0,log(K0));
    original_dist = normpdf(beta_f_space,beta_f,sqrt(var_beta_f));
    nordhaus_dist_0 = normpdf(beta_f_space,mean_distort_nordhaus+beta_f,1./sqrt(lambda_tilde_nordhaus));
    nordhaus = nordhaus_dist_0;
    original = original_dist;
    save([ambiguity,'_',damage_level,'_Dist_50yr'],'nordhaus','original')

    % year 75
    time=300;
    R0 = R_value(time,1);
    K0 = K_value(time,1);
    F0 = T_value(time,1);

    mean_distort_nordhaus = beta_tilde_1_func(log(R0),F0,log(K0))-beta_f;
    lambda_tilde_nordhaus = lambda_tilde_1_func(log(R0),F0,log(K0));
    original_dist = normpdf(beta_f_space,beta_f,sqrt(var_beta_f));
    nordhaus_dist_0 = normpdf(beta_f_space,mean_distort_nordhaus+beta_f,1./sqrt(lambda_tilde_nordhaus));
    nordhaus = nordhaus_dist_0;
    original = original_dist;
    save([ambiguity,'_',damage_level,'_Dist_75yr'],'nordhaus','original')

    % year 100
    time=400;
    R0 = R_value(time,1);
    K0 = K_value(time,1);
    F0 = T_value(time,1);

    mean_distort_nordhaus = beta_tilde_1_func(log(R0),F0,log(K0))-beta_f;
    lambda_tilde_nordhaus = lambda_tilde_1_func(log(R0),F0,log(K0));
    original_dist = normpdf(beta_f_space,beta_f,sqrt(var_beta_f));
    nordhaus_dist_0 = normpdf(beta_f_space,mean_distort_nordhaus+beta_f,1./sqrt(lambda_tilde_nordhaus));
    nordhaus = nordhaus_dist_0;
    original = original_dist;
    save([ambiguity,'_',damage_level,'_Dist_100yr'],'nordhaus','original');

% Case 2: Low damage
elseif strcmp(damage_level,'low')
    % load results
    file1 = [ambiguity,'_',damage_level];
    Model1 = load(file1,'R_1','pi_tilde_1_norm','r_mat','k_mat','t_mat',...
        'beta_tilde_1','beta_f','lambda_tilde_1','theta','var_beta_f'); 

    RE_1 = Model1.R_1;
    pi_tilde_1 = Model1.pi_tilde_1_norm;
    beta_tilde_1 = Model1.beta_tilde_1;
    beta_f = Model1.beta_f;
    lambda_tilde_1 = Model1.lambda_tilde_1;
    theta = Model1.theta;
    var_beta_f = Model1.var_beta_f;
    r_mat = Model1.r_mat;
    t_mat = Model1.t_mat;
    k_mat = Model1.k_mat;
    a = beta_f-10.*sqrt(var_beta_f);
    b = beta_f+10.*sqrt(var_beta_f);

    RE_1_func = griddedInterpolant(r_mat,t_mat,k_mat,RE_1,'linear');
    pi_tilde_1_func = griddedInterpolant(r_mat,t_mat,k_mat,pi_tilde_1,'linear');
    beta_tilde_1_func = griddedInterpolant(r_mat,t_mat,k_mat,beta_tilde_1,'linear');
    lambda_tilde_1_func = griddedInterpolant(r_mat,t_mat,k_mat,lambda_tilde_1,'linear');

    file1 = [ambiguity,'_',damage_level,'_sim'];
    Model1 = load(file1,'hists2'); 
    hists2 = Model1.hists2;
    T_value = mean(squeeze(hists2(:,3,:)),2);
    R_value = mean(squeeze(hists2(:,1,:)),2);
    K_value = mean(squeeze(hists2(:,2,:)),2);

    % relative entropy and weights
    for time=1:400
        RE_1_plot(time) = RE_1_func(log(R_value(time,1)),T_value(time,1),log(K_value(time,1)));
        weight_plot(time) = pi_tilde_1_func(log(R_value(time,1)),T_value(time,1),log(K_value(time,1)));
        nordhaus_mean(time) = beta_tilde_1_func(log(R_value(time,1)),T_value(time,1),log(K_value(time,1)));
        sd_nordhaus(time) = 1./sqrt(lambda_tilde_1_func(log(R_value(time,1)),T_value(time,1),log(K_value(time,1))));
    end

    fileID = fopen([ambiguity,'_',damage_level,'_Relative Entropy.txt'],'w');
    fprintf(fileID,'xi_a: %.6f \n',1./theta);
    fprintf(fileID,'0 yr: %.6f \n',RE_1_plot(1));
    fprintf(fileID,'25 yr: %.6f \n',RE_1_plot(100));
    fprintf(fileID,'50 yr: %.6f \n',RE_1_plot(200));
    fprintf(fileID,'75 yr: %.6f \n',RE_1_plot(300));
    fprintf(fileID,'100 yr: %.6f \n',RE_1_plot(400));
    fclose(fileID);

    fileID = fopen([ambiguity,'_',damage_level,'_Shifted Mean.txt'],'w');
    fprintf(fileID,'xi_a: %.6f \n',1./theta);
    fprintf(fileID,'0 yr: %.6f \n',nordhaus_mean(1));
    fprintf(fileID,'25 yr: %.6f \n',nordhaus_mean(100));
    fprintf(fileID,'50 yr: %.6f \n',nordhaus_mean(200));
    fprintf(fileID,'75 yr: %.6f \n',nordhaus_mean(300));
    fprintf(fileID,'100 yr: %.6f \n',nordhaus_mean(400));
    fclose(fileID);

    fileID = fopen([ambiguity,'_',damage_level,'Shifted Standard Deviation.txt'],'w');
    fprintf(fileID,'xi_a: %.6f \n',1./theta);
    fprintf(fileID,'0 yr: %.6f \n',sd_nordhaus(1));
    fprintf(fileID,'25 yr: %.6f \n',sd_nordhaus(100));
    fprintf(fileID,'50 yr: %.6f \n',sd_nordhaus(200));
    fprintf(fileID,'75 yr: %.6f \n',sd_nordhaus(300));
    fprintf(fileID,'100 yr: %.6f \n',sd_nordhaus(400));
    fclose(fileID);

    fileID = fopen([ambiguity,'_',damage_level,'_Nordhaus Weight.txt'],'w');
    fprintf(fileID,'xi_a: %.6f \n',1./theta);
    fprintf(fileID,'0 yr: %.6f \n',weight_plot(1));
    fprintf(fileID,'25 yr: %.6f \n',weight_plot(100));
    fprintf(fileID,'50 yr: %.6f \n',weight_plot(200));
    fprintf(fileID,'75 yr: %.6f \n',weight_plot(300));
    fprintf(fileID,'100 yr: %.6f \n',weight_plot(400));
    fclose(fileID);

    beta_f_space = linspace(a,b,200); % space for beta_f
    save([ambiguity,'_',damage_level,'_beta_f_space'],'beta_f_space')

    %%% probabilitiess

    % year 0
    time=1;
    R0 = R_value(time,1);
    K0 = K_value(time,1);
    F0 = T_value(time,1);

    mean_distort_nordhaus = beta_tilde_1_func(log(R0),F0,log(K0))-beta_f;
    lambda_tilde_nordhaus = lambda_tilde_1_func(log(R0),F0,log(K0));
    original_dist = normpdf(beta_f_space,beta_f,sqrt(var_beta_f));
    nordhaus_dist_0 = normpdf(beta_f_space,mean_distort_nordhaus+beta_f,1./sqrt(lambda_tilde_nordhaus));
    nordhaus = nordhaus_dist_0;
    original = original_dist;
    save([ambiguity,'_',damage_level,'_Dist_0yr'],'nordhaus','original')

    % year 25
    time=100;
    R0 = R_value(time,1);
    K0 = K_value(time,1);
    F0 = T_value(time,1);

    mean_distort_nordhaus = beta_tilde_1_func(log(R0),F0,log(K0))-beta_f;
    lambda_tilde_nordhaus = lambda_tilde_1_func(log(R0),F0,log(K0));
    original_dist = normpdf(beta_f_space,beta_f,sqrt(var_beta_f));
    nordhaus_dist_0 = normpdf(beta_f_space,mean_distort_nordhaus+beta_f,1./sqrt(lambda_tilde_nordhaus));
    nordhaus = nordhaus_dist_0;
    original = original_dist;
    save([ambiguity,'_',damage_level,'_Dist_25yr'],'nordhaus','original')

    % year 50
    time=200;
    R0 = R_value(time,1);
    K0 = K_value(time,1);
    F0 = T_value(time,1);

    mean_distort_nordhaus = beta_tilde_1_func(log(R0),F0,log(K0))-beta_f;
    lambda_tilde_nordhaus = lambda_tilde_1_func(log(R0),F0,log(K0));
    original_dist = normpdf(beta_f_space,beta_f,sqrt(var_beta_f));
    nordhaus_dist_0 = normpdf(beta_f_space,mean_distort_nordhaus+beta_f,1./sqrt(lambda_tilde_nordhaus));
    nordhaus = nordhaus_dist_0;
    original = original_dist;
    save([ambiguity,'_',damage_level,'_Dist_50yr'],'nordhaus','original')

    % year 75
    time=300;
    R0 = R_value(time,1);
    K0 = K_value(time,1);
    F0 = T_value(time,1);

    mean_distort_nordhaus = beta_tilde_1_func(log(R0),F0,log(K0))-beta_f;
    lambda_tilde_nordhaus = lambda_tilde_1_func(log(R0),F0,log(K0));
    original_dist = normpdf(beta_f_space,beta_f,sqrt(var_beta_f));
    nordhaus_dist_0 = normpdf(beta_f_space,mean_distort_nordhaus+beta_f,1./sqrt(lambda_tilde_nordhaus));
    nordhaus = nordhaus_dist_0;
    original = original_dist;
    save([ambiguity,'_',damage_level,'_Dist_75yr'],'nordhaus','original')

    % year 100
    time=400;
    R0 = R_value(time,1);
    K0 = K_value(time,1);
    F0 = T_value(time,1);

    mean_distort_nordhaus = beta_tilde_1_func(log(R0),F0,log(K0))-beta_f;
    lambda_tilde_nordhaus = lambda_tilde_1_func(log(R0),F0,log(K0));
    original_dist = normpdf(beta_f_space,beta_f,sqrt(var_beta_f));
    nordhaus_dist_0 = normpdf(beta_f_space,mean_distort_nordhaus+beta_f,1./sqrt(lambda_tilde_nordhaus));
    nordhaus = nordhaus_dist_0;
    original = original_dist;
    save([ambiguity,'_',damage_level,'_Dist_100yr'],'nordhaus','original');

% Case 3: Weighted damage    
else
    % load results
    file1 = [ambiguity,'_',damage_level];
    Model1 = load(file1,'R_1','pi_tilde_1_norm','r_mat','k_mat','t_mat',...
        'beta_tilde_1','beta_f','lambda_tilde_1','theta','var_beta_f'); 

    RE_1 = Model1.R_1;
    pi_tilde_1 = Model1.pi_tilde_1_norm;
    beta_tilde_1 = Model1.beta_tilde_1;
    beta_f = Model1.beta_f;
    lambda_tilde_1 = Model1.lambda_tilde_1;
    theta = Model1.theta;
    var_beta_f = Model1.var_beta_f;
    r_mat = Model1.r_mat;
    t_mat = Model1.t_mat;
    k_mat = Model1.k_mat;
    a = beta_f-10.*sqrt(var_beta_f);
    b = beta_f+10.*sqrt(var_beta_f);

    RE_1_func = griddedInterpolant(r_mat,t_mat,k_mat,RE_1,'linear');
    pi_tilde_1_func = griddedInterpolant(r_mat,t_mat,k_mat,pi_tilde_1,'linear');
    beta_tilde_1_func = griddedInterpolant(r_mat,t_mat,k_mat,beta_tilde_1,'linear');
    lambda_tilde_1_func = griddedInterpolant(r_mat,t_mat,k_mat,lambda_tilde_1,'linear');

    file1 = [ambiguity,'_',damage_level,'_sim'];
    Model1 = load(file1,'hists2'); 
    hists2 = Model1.hists2;
    T_value = mean(squeeze(hists2(:,3,:)),2);
    R_value = mean(squeeze(hists2(:,1,:)),2);
    K_value = mean(squeeze(hists2(:,2,:)),2);

    % relative entropy and weights
    for time=1:400
        RE_1_plot(time) = RE_1_func(log(R_value(time,1)),T_value(time,1),log(K_value(time,1)));
        weight_plot(time) = pi_tilde_1_func(log(R_value(time,1)),T_value(time,1),log(K_value(time,1)));
        nordhaus_mean(time) = beta_tilde_1_func(log(R_value(time,1)),T_value(time,1),log(K_value(time,1)));
        sd_nordhaus(time) = 1./sqrt(lambda_tilde_1_func(log(R_value(time,1)),T_value(time,1),log(K_value(time,1))));
    end

    fileID = fopen([ambiguity,'_',damage_level,'_Relative Entropy.txt'],'w');
    fprintf(fileID,'xi_a: %.6f \n',1./theta);
    fprintf(fileID,'0 yr: %.6f \n',RE_1_plot(1));
    fprintf(fileID,'25 yr: %.6f \n',RE_1_plot(100));
    fprintf(fileID,'50 yr: %.6f \n',RE_1_plot(200));
    fprintf(fileID,'75 yr: %.6f \n',RE_1_plot(300));
    fprintf(fileID,'100 yr: %.6f \n',RE_1_plot(400));
    fclose(fileID);

    fileID = fopen([ambiguity,'_',damage_level,'_Shifted Mean.txt'],'w');
    fprintf(fileID,'xi_a: %.6f \n',1./theta);
    fprintf(fileID,'0 yr: %.6f \n',nordhaus_mean(1));
    fprintf(fileID,'25 yr: %.6f \n',nordhaus_mean(100));
    fprintf(fileID,'50 yr: %.6f \n',nordhaus_mean(200));
    fprintf(fileID,'75 yr: %.6f \n',nordhaus_mean(300));
    fprintf(fileID,'100 yr: %.6f \n',nordhaus_mean(400));
    fclose(fileID);

    fileID = fopen([ambiguity,'_',damage_level,'_Shifted Standard Deviation.txt'],'w');
    fprintf(fileID,'xi_a: %.6f \n',1./theta);
    fprintf(fileID,'0 yr: %.6f \n',sd_nordhaus(1));
    fprintf(fileID,'25 yr: %.6f \n',sd_nordhaus(100));
    fprintf(fileID,'50 yr: %.6f \n',sd_nordhaus(200));
    fprintf(fileID,'75 yr: %.6f \n',sd_nordhaus(300));
    fprintf(fileID,'100 yr: %.6f \n',sd_nordhaus(400));
    fclose(fileID);

    fileID = fopen([ambiguity,'_',damage_level,'_Nordhaus Weight.txt'],'w');
    fprintf(fileID,'xi_a: %.6f \n',1./theta);
    fprintf(fileID,'0 yr: %.6f \n',weight_plot(1));
    fprintf(fileID,'25 yr: %.6f \n',weight_plot(100));
    fprintf(fileID,'50 yr: %.6f \n',weight_plot(200));
    fprintf(fileID,'75 yr: %.6f \n',weight_plot(300));
    fprintf(fileID,'100 yr: %.6f \n',weight_plot(400));
    fclose(fileID);

    beta_f_space = linspace(a,b,200); % space for beta_f
    save([ambiguity,'_',damage_level,'_beta_f_space'],'beta_f_space')

    %%% probabilitiess

    % year 0
    time=1;
    R0 = R_value(time,1);
    K0 = K_value(time,1);
    F0 = T_value(time,1);

    mean_distort_nordhaus = beta_tilde_1_func(log(R0),F0,log(K0))-beta_f;
    lambda_tilde_nordhaus = lambda_tilde_1_func(log(R0),F0,log(K0));
    original_dist = normpdf(beta_f_space,beta_f,sqrt(var_beta_f));
    nordhaus_dist_0 = normpdf(beta_f_space,mean_distort_nordhaus+beta_f,1./sqrt(lambda_tilde_nordhaus));
    nordhaus = nordhaus_dist_0;
    original = original_dist;
    save([ambiguity,'_',damage_level,'_Dist_0yr'],'nordhaus','original')

    % year 25
    time=100;
    R0 = R_value(time,1);
    K0 = K_value(time,1);
    F0 = T_value(time,1);

    mean_distort_nordhaus = beta_tilde_1_func(log(R0),F0,log(K0))-beta_f;
    lambda_tilde_nordhaus = lambda_tilde_1_func(log(R0),F0,log(K0));
    original_dist = normpdf(beta_f_space,beta_f,sqrt(var_beta_f));
    nordhaus_dist_0 = normpdf(beta_f_space,mean_distort_nordhaus+beta_f,1./sqrt(lambda_tilde_nordhaus));
    nordhaus = nordhaus_dist_0;
    original = original_dist;
    save([ambiguity,'_',damage_level,'_Dist_25yr'],'nordhaus','original')

    % year 50
    time=200;
    R0 = R_value(time,1);
    K0 = K_value(time,1);
    F0 = T_value(time,1);

    mean_distort_nordhaus = beta_tilde_1_func(log(R0),F0,log(K0))-beta_f;
    lambda_tilde_nordhaus = lambda_tilde_1_func(log(R0),F0,log(K0));
    original_dist = normpdf(beta_f_space,beta_f,sqrt(var_beta_f));
    nordhaus_dist_0 = normpdf(beta_f_space,mean_distort_nordhaus+beta_f,1./sqrt(lambda_tilde_nordhaus));
    nordhaus = nordhaus_dist_0;
    original = original_dist;
    save([ambiguity,'_',damage_level,'_Dist_50yr'],'nordhaus','original')

    % year 75
    time=300;
    R0 = R_value(time,1);
    K0 = K_value(time,1);
    F0 = T_value(time,1);

    mean_distort_nordhaus = beta_tilde_1_func(log(R0),F0,log(K0))-beta_f;
    lambda_tilde_nordhaus = lambda_tilde_1_func(log(R0),F0,log(K0));
    original_dist = normpdf(beta_f_space,beta_f,sqrt(var_beta_f));
    nordhaus_dist_0 = normpdf(beta_f_space,mean_distort_nordhaus+beta_f,1./sqrt(lambda_tilde_nordhaus));
    nordhaus = nordhaus_dist_0;
    original = original_dist;
    save([ambiguity,'_',damage_level,'_Dist_75yr'],'nordhaus','original')

    % year 100
    time=400;
    R0 = R_value(time,1);
    K0 = K_value(time,1);
    F0 = T_value(time,1);

    mean_distort_nordhaus = beta_tilde_1_func(log(R0),F0,log(K0))-beta_f;
    lambda_tilde_nordhaus = lambda_tilde_1_func(log(R0),F0,log(K0));
    original_dist = normpdf(beta_f_space,beta_f,sqrt(var_beta_f));
    nordhaus_dist_0 = normpdf(beta_f_space,mean_distort_nordhaus+beta_f,1./sqrt(lambda_tilde_nordhaus));
    nordhaus = nordhaus_dist_0;
    original = original_dist;
    save([ambiguity,'_',damage_level,'_Dist_100yr'],'nordhaus','original')
end
