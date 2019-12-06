# Preference Damage Model

## Code Explanation:
The directory contains three folders corresponding to three cases of the model from the paper. 
Under each folder there are 5 files:
1. AmbiguityAverse.m : This generates results for the ambiguity averse model.
2. AmbiguityNeutral.m : This generates results for ambiguity neutral model
3. quad_int.m : This calculates the intergral for a given function using quadrature rules.
4. quad_points_legendre : This generates weights and points for legendre quadrature
5. TCRE_MacDougallEtAl2017_update.csv : Parameter estimates from MacDougal et. al (2017)

To generate paper results, first set up the solvers (see instructions for solvers). Next run AmbiguityAverse.m and AmbiguityNeutral.m. Refer back to the notebook for further details.
