# Expanded Preference Code

This code explores using 4 damage functions with quasi-analytical solutions, as opposed to the regular preference code which uses one model with quasi-analytical solutions and one model which requires numerical integration to evaluate.

Run the code in the following order:
1) HJB.m: Code solving the HJB equation
2) Sims.m: Simulate an economy based on the HJB solution
3) WorstCase.m: Generate plots for the worst case distributions
4) SCC_EX.m: Calculate the external component of SCC using a Feynman-Kac equation
5) Generate_Emissions.m: Plot the emissions
6) Generate_SCC.m: Plot the various components of SCC
