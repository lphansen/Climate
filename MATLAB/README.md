## Acessing our MATLAB codes
To copy the code to your machine, you may either download it from the Github website directly or you may clone the repository in read-only mode.
Jupyter notebooks located in the seperate subfolder provided a clean version of the code we used while we still upload our MATLAB source codes are provided here for user to play with. They're structured as below:

### Growth_paper_appendix
This folder contains Matlab scripts that generate the results related to the growth damage model.

### Preference_solution
This folder contains Matlab scripts that generate the paper results related to the consumption damage model.

### Preference_appendix
This folder contains Matlab scripts that generate the appendix results related to the consumption damage model.

### Solver
This folder contains eigen3 folder and MarketIO.h file.
eigen3 is a C++ template library for linear algebra: matrices, vectors, numerical solvers, and related algorithms. (Source:http://eigen.tuxfamily.org/)

### Prerequisites
This project requires users to have access to Matlab software. For license, please visit https://www.mathworks.com/

### Preliminary steps
As an example, if users wish to generate results for the growth model, first navigate to "Solver" folder and unzip the "eigen3" folder. Then, copy the "eigen3" folder and MarketIO.h file into "Growth_paper_appendix" folder. Follow the instructions under "Growth_paper_appendix" folder for more details.
