
# Pricing Uncertainty Induced by Climate Change
This repository contains codes and jupyter notebooks which estimates and demonstrates results of the empirical model in "Pricing Uncertainty Induced by Climate Change" by [Michael Barnett][id3], [William Brock][id2] and [Lars Peter Hansen][id1]. Latest version could be found [here][id4].

[id1]: https://larspeterhansen.org/
[id2]: https://www.ssc.wisc.edu/~wbrock/
[id3]: https://sites.google.com/site/michaelduglasbarnett/home
[id4]: https://larspeterhansen.org/research/papers/

## Acessing our jupyter notebook
To access our notebook, please follow steps below:
1.	Open a Windows command prompt or Mac command terminal and change into the folder you would like to store the files. 
    - You can do this using the command __cd__ in the command prompt.    
    - For example, running "cd “C:\Users\jiamingwang\Documents\Python Projects” " (don’t forget “” around the path name to use an absolute path) would lead me to my designated folder.
```
cd [folder path name]” 
```
2.	Clone the github repository for the paper 
    - If you don’t have github installed, try installing it from this page: https://git-scm.com/download/win.
    - You can do this by running in the command prompt. 
```
git clone https://github.com/lphansen/Climate
```
3.	To download the data you'll need, go to https://drive.google.com/drive/folders/1TVVqditjA489NNSD8H32u0e4eIJkoj0r?usp=sharing and copy the .pickle files into the subfolder __data__ under the git repo __Climate__ you just cloned.
4.	Run __install.bat__ for those who used windows or run __install.sh__ for those who used Mac OS by double clicking the file.
5.	Go back to command line prompt, change directories into the ‘Climate’ folder and open jupyter notebook by running below in command prompt
```
cd Climate
jupyter notebook
```
6.	"PaperResultsIllustration.ipynb" demonstrates interactive figures shown in the paper. “PreferenceModel.ipynb” demonstrates the computational details and codes on how we solve the model under consumption damage settings while “GrowthModel.ipynb” demonstrates growth damage settings. Click on the files to access them.

## Acessing our MATLAB codes
To copy the code to your machine, you may either download it from the Github website directly or you may clone the repository in read-only mode.
Above jupyter notebooks provided a clean version of the code we used while we still upload our MATLAB source codes are provided here for user to play with. They're structured as below:

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

