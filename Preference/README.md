# Preference Damage Model

## Code Explanation:
1. RunMe.m: run this file to generate results.
2. Parameters_config1.m: This is the file which contains parameters documented in the model
3. HJB.m: This is the main file for solving the model's HJB equation. 
4. HJB_NoUn.m: This is the main file for solving the model's HJB equation with no ambiguity. 
5. Sims.m: This code simulates the model based upon the solution to the HJB equation.
6. Sims_NoUn.m: This code simulates the model based upon the solution to the HJB equation with no ambiguity.
7. SCC_EX.m: Solves the Feynman-Kac equation associated with the model, yielding the external component of the social cost of carbon.
8. SCC_Plot.m: This pieces together relevant components from the HJB and F-K equations to construct and plot different pieces of the social cost of carbon.
9. SCC_NoUn_Plot.m: This pieces together relevant components from the HJB and F-K equations to construct and plot different pieces of the social cost of carbon for no ambiguity setting.
10. WorstCase.m: This calculates the worst case distributions under each damage equation.
11. quad_int.m: This calculates the intergral for a given function using quadrature rules.
12. quad_points_legendre: This generates weights and points for legendre quadrature
13. quad_points_hermite: This generates weights and points for hermite quadrature
