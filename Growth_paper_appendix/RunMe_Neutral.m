
   %%%%% This file calls all other scripts in order to generate the result
   %%%%% in Table 4 in the paper and online appendix for the ambiguity
   %%%%% neutral model
   
   close all;
   clc;
   clear all;

   run('Growth_Neutral')
   run('SCC_neutral_base_1e9')
   run('SCC_neutral_worst_1e9')
   run('Sims_Jieyao_Neutral')
   run('SCC_plot_JY_neutral')
   
   