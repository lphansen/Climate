
   %%%%% This file calls all other scripts in order to generate the result
   %%%%% in Table 4 in the paper and online appendix for the ambiguity
   %%%%% neutral model
   
   close all;
   clc;
   clear all;

   run('Growth_Neutral_approach_1')
   % run('Growth_Neutral_approach_2')
   run('SCC_neutral_base_1e14')
   run('SCC_neutral_worst_1e14')
   run('Sims_Jieyao_Neutral')
   run('SCC_plot_JY_neutral')
   
   
