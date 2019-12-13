# Pricing Uncertainty Induced by Climate Change

This folder contains the Matlab codes to replicate the numerical results for the growth damage model for the paper

## Prerequisites

Please follow instructions under the parent folder, and copy "eigen3" folder and MarketIO.h file from the "Solver" folder into this folder.

### Installing and activating the environment

For both Mac and Windows Users, please navigate to the folder, and open any .m script file in Matlab.
Then, in the command prompt, type in the following command:

```
mex solveCGNatural_1_1e14.cpp
mex solveCGNatural_1e10.cpp
mex solveCGNatural_1e14.cpp
mex solveCGNatural_1e16.cpp
```

To ensure success of the mex process for the solver, make sure that the eigen3 folder, MarketIO.h file and cpp files above are under the same folder.

## Generating results

### Running code to get data

To run the code, simply run the following two .m scripts in Matlab

```
RunMe_Averse.m
RunMe_Neutral.m
```
These 2 scripts will call other scripts in the folder, and will generate numerical results for the growth damage model. The numbers go into the Table 4 in the paper, and Table G.3 in the online appendix.
Users may also modify these two files and choose which numerical approach they wish to take in order to generate numerical results. The default is approach one; users may comment out approach one and use approach two instead by putting a percentage mark in front of line 10 and deleting the percentage mark in front of line 11 to achieve this.
For detailed description of two approaches users may refer to online Jupyter notebook.

## Authors

* **Michael Barnett** - *Michael.D.Barnett@asu.edu*
* **Jieyao Wang** - *wangjieyao@uchicago.edu*

## Acknowledgments

* We thank Joseph Huang for providing assistance on software implementation, and John Wilson for providing assistance on quadrature methods.