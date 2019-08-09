# Growth Damage Model

## Sequence for running code:
 1. HJB.m: This is the main file for solving the model's HJB equation. 
 2. Sims.m: This code simulates the model based upon the solution to the HJB equation.
 3. WorstCase.m: This calculates the worst case distributions under each damage equation.
 4. SCC_EX.m: Solves the Feynman-Kac equation associated with the model, yielding the external component of the social cost of carbon.
 5. Graph_SCC.m: This pieces together relevant components from the HJB and F-K equations to construct and plot different pieces of the social cost of carbon.
