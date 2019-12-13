# Preference Damage Model

## Prerequisites

Please follow instructions under the parent folder, and copy "eigen3" folder and MarketIO.h file from the "Solver" folder into this folder.

## Code Explanation:
This directory contains code used to solve the model in the case of consumption damages. A few simple files generate the results reported in the paper:

1. HJB_Solution.m: Solves the model HJB and simulates a time-path trajectory.
2. SCC_Solution.m: Solves the Feynman-Kac equations for the Social Cost of Carbon calculation and simulates a time-path trajectory. 
3. Ambiguity_adjusted_probability.m: Computes the ambiguity-adjusted probability distributions at different time horizons based on the simulated trajectory.

To generate paper results, first set up the solvers (see instructions for solvers). Next run items 1-3 above in sequence. Users can specify alternative ambiguity and damage levels in item 1. Refer back to the notebook for further details. 
