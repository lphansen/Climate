# Pricing Uncertainty Induced by Climate Change

This folder contains the Matlab codes to replicate the numerical results for the growth damage model for the paper

## Getting Started

To copy the code to your machine, you may either download it from the Github website directly or you may clone the repository in read-only mode.

### Prerequisites

This project requires users have access to Matlab software. For license, please visit https://www.mathworks.com/

### Installing and activating the environment

For both Mac and Windows Users, please navigate to the folder, and open any .m script file in Matlab.
Then, in the command prompt, type in the following command:

```
mex solveCGNatural_1_1e9.cpp
mex solveCGNatural_1e10.cpp
```

## Generating results

### Running code to get data

To run the code, simply run the following two .m scripts in Matlab

```
RunMe_Averse.m
RunMe_Neutral.m
```
These 2 scripts will call other scripts in the folder, and will generate numerical results for the growth damage model. The numbers go into the Table 4 in the paper, and Table G.3 in the online appendix.

## Authors

* **Michael Barnett** - *Michael.D.Barnett@asu.edu*
* **Jieyao Wang** - *wangjieyao@uchicago.edu*

## Acknowledgments

* We thank Joseph Huang for providing assistance on software implementation, and John Wilson for providing assistance on quadrature methods.
