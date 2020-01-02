/* ****************** */
/* Include packages   */
/* ****************** */
#include <math.h>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <algorithm>
#include <random>
#include <iostream>
#include <time.h>
#include <stack>
#include <assert.h>
//include Eigen
// #include "eigen3/Eigen/Dense"
// #include "eigen3/Eigen/Sparse"
// #include "eigen3/Eigen/SparseLU"
// #include "eigen3/Eigen/SparseQR"
// #include "eigen3/Eigen/SparseCholesky"
// #include "eigen3/Eigen/IterativeLinearSolvers"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/SparseQR>
#include <Eigen/SparseCholesky>
#include <Eigen/IterativeLinearSolvers>
typedef Eigen::SparseMatrix<double > SpMat;
typedef Eigen::Triplet<double> T;
#include <typeinfo>
#include "MarketIO.h"
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>

//include MEX related files
// #include "mex.h"
// #include "matrix.h"

/* ******************** */
/* State Variable Class */
/* ******************** */

Eigen::MatrixXd empty;
Eigen::ArrayXd emptyAry;



namespace py = pybind11;
using namespace std;
using MatrixXdR = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;


/* Timer functions                                    */
/******************************************************/

std::stack<clock_t> tictoc_stack;

void tic() {
    tictoc_stack.push(clock());
}

void toc() {
    std::cout << "Time elapsed: "
    << ((double)(clock() - tictoc_stack.top())) / CLOCKS_PER_SEC
    << std::endl;
    tictoc_stack.pop();
}

class stateVars {
    
public:
    Eigen::MatrixXd stateMat; //matrix to store state variables
    Eigen::MatrixXd stateMatNorm; //matrix to store normalized state variables [-1,1]
    Eigen::ArrayXd increVec; //vector to record steps
    Eigen::ArrayXd dVec; //vector to record steps
    int N; // num of dimensions
    int S; // number of rows for the grid
    Eigen::ArrayXd upperLims;
    Eigen::ArrayXd lowerLims;
    Eigen::ArrayXd gridSizes;

    stateVars (Eigen::ArrayXd, Eigen::ArrayXd, Eigen::ArrayXd); //constructors with arrays of upper/lower bounds and gridsizes 
    stateVars (Eigen::MatrixXd); //constructors by loading in data

};


stateVars::stateVars (Eigen::ArrayXd upper, Eigen::ArrayXd lower, Eigen::ArrayXd gridSizes) {
    
    upperLims = upper;
    lowerLims = lower;
    N = upperLims.size();
    S = gridSizes.prod();
    stateMat.resize(S,N);
    dVec.resize(N);
    increVec.resize(N);
    increVec(0) = 1;
        
    //fill in the state object; similar to the ndgrid function in MATLAB
    
    for (int n = 0; n < N; ++n) {
            
        if (n != 0) {
            increVec(n) = gridSizes(n - 1) * increVec(n - 1);
        }
        dVec(n) = (upper(n) - lower(n)) / (gridSizes(n) - 1);
        
        for (int i = 0; i < S; ++i) {
            stateMat(i,n) = lower(n) + dVec(n) * ( int(i /  increVec(n) ) % int( gridSizes(n) ) );
        }
            
    }
    
}

stateVars::stateVars (Eigen::MatrixXd preLoad) {

    //fill in stateMat based on data loaded
    N = preLoad.cols();
    S = preLoad.rows();
    stateMat.resize(S,N);
    dVec.resize(N); dVec.setZero();
    increVec.resize(N); increVec.setZero();
    upperLims.resize(N); lowerLims.resize(N);
    for (int j = 0; j < preLoad.cols(); ++j) {
        upperLims(j) = preLoad.col(j).maxCoeff();
        lowerLims(j) = preLoad.col(j).minCoeff();
    }
    
    stateMat = preLoad;
    
    //figure out dVec and increVec
    for (int i = 1; i < S; ++i) {
        for (int n = 0; n < N; ++n ) {
            double diff = stateMat(i,n) - stateMat(i-1,n);
            if (diff > 0 && dVec(n) == 0 && increVec(n) == 0) {
                dVec(n) = diff;
                increVec(n) = i;
            }
        }
        
    }
    


}

struct bc {
    double a0;
    double a0S;
    bool natural;
    Eigen::ArrayXd level;
    Eigen::ArrayXd first;
    Eigen::ArrayXd second;
    
    bc(int d) {
        level.resize(d); first.resize(d); second.resize(d);
    }
};

struct elas {
    Eigen::MatrixXd elas1sc;
    Eigen::MatrixXd elas1c; //exposure elas
    Eigen::MatrixXd elas2sc;
    Eigen::MatrixXd elas2c; //exposure elas
    Eigen::MatrixXd elas1p; //price elas
    Eigen::MatrixXd elas2p;  //price elas
    
    elas(int T, int S) {
        elas1sc.resize(S,T); elas1c.resize(S,T); elas1p.resize(S,T);
        elas2sc.resize(S,T); elas2c.resize(S,T); elas2p.resize(S,T);
    }
};




class linearSysVars {
    
public:
    double dt;
    int k;
    Eigen::MatrixXd A; 
    Eigen::MatrixXd B;
    Eigen::MatrixXd C;
    Eigen::MatrixXd D;

    Eigen::ArrayXd atBoundIndicators;

    std::vector<T> matList; 
    std::string solverType;
    SpMat Le;

    //member functions
    
    //constructor
    linearSysVars(stateVars & state_vars, Eigen::MatrixXd A, Eigen::MatrixXd B, Eigen::MatrixXd C, Eigen::MatrixXd D, double dt);
    
    //function to construt matrix
    
    void constructMatFT(stateVars & state_vars);
    void constructMatFK(stateVars & state_vars);

};

linearSysVars::linearSysVars(stateVars & state_vars, Eigen::MatrixXd AInput, Eigen::MatrixXd BInput, Eigen::MatrixXd CInput, Eigen::MatrixXd DInput, double dtInput) {
        
    Le.resize(state_vars.S,state_vars.S);
    A.resize(state_vars.S,1); B.resize(state_vars.S,state_vars.N);
    C.resize(state_vars.S,state_vars.N); D.resize(state_vars.S,1);
    A = AInput; B = BInput; C = CInput; D = DInput;
    dt = dtInput;

}



void linearSysVars::constructMatFT(stateVars & state_vars) {
    matList.clear();
    matList.reserve(10 * state_vars.S);
    atBoundIndicators.resize(state_vars.N);
    double atBound = -1;
    double upperBound = -1;
    //construct matrix

    for (int i = 0; i < state_vars.S; ++i) {
        //level and time deriv
        
        atBound = -1;
        //check boundaries
        
        matList.push_back(T(i,i, (1.0 - dt * A(i,0))  ));
        
        for (int n  = (state_vars.N - 1); n >=0; --n ) {
            
            atBoundIndicators(n) = -1.0;
            
            double firstCoefE = B(i,n);
            
            double secondCoefE = C(i,n);
            
            //check whether it's at upper or lower boundary
            if ( std::abs(state_vars.stateMat(i,n) - state_vars.upperLims(n)) < state_vars.dVec(n)/2.0 ) {  //upper boundary
                atBoundIndicators(n) = 1.0;
        
                atBound = 1.0;
                upperBound = 1.0;
                /* Uncomment this section if you want natural boundaries */
                
                 matList.push_back(T(i, i, - dt * ( firstCoefE/state_vars.dVec(n) + secondCoefE / pow(state_vars.dVec(n), 2) ) ) );
                 matList.push_back(T(i, i - state_vars.increVec(n), - dt * ( - firstCoefE/state_vars.dVec(n) - 2 * secondCoefE / pow(state_vars.dVec(n), 2) ) ));
                 matList.push_back(T(i, i - 2*state_vars.increVec(n), - dt * ( secondCoefE / pow(state_vars.dVec(n), 2) ) ) );
                
                /* Uncomment this section if you want first derivatives = constant  
                 matList.push_back(T(i,i, - (1.0 - dt * A(i,0) ) ));
                /*
                 matList.push_back(T(i, i, - dt * ( 1.0/state_vars.dVec(n)  ) ) );
                 matList.push_back(T(i, i - state_vars.increVec(n), - dt * ( - 1.0/state_vars.dVec(n)  ) ));*/
                 /*
                if ((n == 0) && atBoundIndicators(1) > 0 ) {
                 matList.push_back(T(i, i,  dt * ( firstCoefE/state_vars.dVec(n)  ) ) );
                 matList.push_back(T(i, i - state_vars.increVec(n),  dt * ( - firstCoefE/state_vars.dVec(n)  ) ));
                }*/
                /* Uncomment this section if you want second derivatives = constant  */
                //matList.push_back(T(i,i, - (1.0 - dt * A(i,0) )  ));   
                /*
                matList.push_back(T(i, i, - dt * (  secondCoefE / pow(state_vars.dVec(n), 2) ) ) );
                matList.push_back(T(i, i - state_vars.increVec(n), - dt * ( - 2 * secondCoefE / pow(state_vars.dVec(n), 2) ) ));
                matList.push_back(T(i, i - 2*state_vars.increVec(n), - dt * ( secondCoefE / pow(state_vars.dVec(n), 2) ) ) );
                */
                /*
                matList.push_back(T(i, i, - dt * (  1.0 / pow(state_vars.dVec(n), 2) ) ) );
                matList.push_back(T(i, i - state_vars.increVec(n), - dt * ( - 2 *  1.0 / pow(state_vars.dVec(n), 2) ) ));
                matList.push_back(T(i, i - 2*state_vars.increVec(n), - dt * ( 1.0 / pow(state_vars.dVec(n), 2) ) ) );
                */
            } else if ( std::abs(state_vars.stateMat(i,n) - state_vars.lowerLims(n)) < state_vars.dVec(n)/2.0 ) { //lower boundary
                
                atBoundIndicators(n) = 1.0;
                atBound = 1.0;

                ///* Uncomment this section if you want natural boundaries
                
                 matList.push_back(T(i, i, - dt * ( - firstCoefE/state_vars.dVec(n) + secondCoefE / pow(state_vars.dVec(n), 2) ) ) );
                 matList.push_back(T(i, i + state_vars.increVec(n), - dt * ( firstCoefE/state_vars.dVec(n) - 2 * secondCoefE / pow(state_vars.dVec(n), 2) ) ));
                 matList.push_back(T(i, i + 2*state_vars.increVec(n), - dt * ( secondCoefE / pow(state_vars.dVec(n), 2) ) ) );
                
                //*/
                /* Uncomment this section if you want first derivatives = constant
                 */
                // matList.push_back(T(i,i, - (1.0 - dt * A(i,0) ) ));
                // matList.push_back(T(i, i, - dt * ( - 1/state_vars.dVec(n)  ) ) );
                // matList.push_back(T(i, i + state_vars.increVec(n), - dt * ( 1/state_vars.dVec(n) ) ));
                 /*
                if ((n == 0) && atBoundIndicators(1) > 0 ) {
                 matList.push_back(T(i, i,  dt * ( - firstCoefE/state_vars.dVec(n)  ) ) );
                 matList.push_back(T(i, i + state_vars.increVec(n),  dt * ( firstCoefE/state_vars.dVec(n) ) ));
                }*/
                /* Uncomment this section if you want second derivatives = constant  */
                //matList.push_back(T(i,i, - (1.0 - dt * A(i,0) )/ state_vars.N ));
                /*
                matList.push_back(T(i, i, - dt * ( secondCoefE / pow(state_vars.dVec(n), 2) ) ) );
                matList.push_back(T(i, i + state_vars.increVec(n), - dt * (  - 2 * secondCoefE / pow(state_vars.dVec(n), 2) ) ));
                matList.push_back(T(i, i + 2*state_vars.increVec(n), - dt * ( secondCoefE / pow(state_vars.dVec(n), 2) ) ) );
                */
                /*
                    matList.push_back(T(i, i, - dt * ( 1.0 / pow(state_vars.dVec(n), 2) ) ) );
                    matList.push_back(T(i, i + state_vars.increVec(n), - dt * (  - 2 * 1.0 / pow(state_vars.dVec(n), 2) ) ));
                    matList.push_back(T(i, i + 2*state_vars.increVec(n), - dt * ( 1.0 / pow(state_vars.dVec(n), 2) ) ) );
                */
                
            }



        }
        
             
        if (atBound < 0 ) {
          // matList.push_back(T(i,i, (1.0 - dt * A(i,0))  ));
        }
        for (int n = (state_vars.N - 1); n >= 0; --n) {
            
            //add elements to the vector of triplets for matrix construction
            if ( atBoundIndicators(n) < 0) {
                double firstCoefE = B(i,n);
                double secondCoefE = C(i,n);
                
                //first derivative
                 matList.push_back(T(i,i, - dt * ( -firstCoefE * ( firstCoefE > 0) + firstCoefE * ( firstCoefE < 0) ) / state_vars.dVec(n)  ) );
                 matList.push_back(T(i,i + state_vars.increVec(n), - dt * firstCoefE * ( firstCoefE > 0) / state_vars.dVec(n) ));
                 matList.push_back(T(i,i - state_vars.increVec(n), - dt *  - firstCoefE * ( firstCoefE < 0) / state_vars.dVec(n) ));
                
                    matList.push_back(T(i, i, - dt * -2 * secondCoefE / ( pow(state_vars.dVec(n), 2) ) ));
                    matList.push_back(T(i, i + state_vars.increVec(n), - dt * secondCoefE / ( pow(state_vars.dVec(n), 2) ) ));
                    matList.push_back(T(i, i - state_vars.increVec(n), - dt * secondCoefE / ( pow(state_vars.dVec(n), 2) ) ));
                
            }

        }


    }
    //form matrices

    Le.setFromTriplets(matList.begin(), matList.end());

    //compress
    Le.makeCompressed(); 
}

void linearSysVars::constructMatFK(stateVars & state_vars) {
    matList.clear();
    matList.reserve(10 * state_vars.S);
    atBoundIndicators.resize(state_vars.N);
    double atBound = -1;
    double upperBound = -1;
    //construct matrix

    for (int i = 0; i < state_vars.S; ++i) {
        //level and time deriv
        
        atBound = -1;
        //check boundaries
        
        matList.push_back(T(i,i, (0.0 - dt * A(i,0))  ));
        
        for (int n  = (state_vars.N - 1); n >=0; --n ) {
            
            atBoundIndicators(n) = -1.0;
            
            double firstCoefE = B(i,n);
            
            double secondCoefE = C(i,n);
            
            //check whether it's at upper or lower boundary
            if ( std::abs(state_vars.stateMat(i,n) - state_vars.upperLims(n)) < state_vars.dVec(n)/2.0 ) {  //upper boundary
                atBoundIndicators(n) = 1.0;
        
                atBound = 1.0;
                upperBound = 1.0;
                /* Uncomment this section if you want natural boundaries */
                
                 matList.push_back(T(i, i, - dt * ( firstCoefE/state_vars.dVec(n) + secondCoefE / pow(state_vars.dVec(n), 2) ) ) );
                 matList.push_back(T(i, i - state_vars.increVec(n), - dt * ( - firstCoefE/state_vars.dVec(n) - 2 * secondCoefE / pow(state_vars.dVec(n), 2) ) ));
                 matList.push_back(T(i, i - 2*state_vars.increVec(n), - dt * ( secondCoefE / pow(state_vars.dVec(n), 2) ) ) );
                
                /* Uncomment this section if you want first derivatives = constant  
                 matList.push_back(T(i,i, - (1.0 - dt * A(i,0) ) ));
                /*
                 matList.push_back(T(i, i, - dt * ( 1.0/state_vars.dVec(n)  ) ) );
                 matList.push_back(T(i, i - state_vars.increVec(n), - dt * ( - 1.0/state_vars.dVec(n)  ) ));*/
                 /*
                if ((n == 0) && atBoundIndicators(1) > 0 ) {
                 matList.push_back(T(i, i,  dt * ( firstCoefE/state_vars.dVec(n)  ) ) );
                 matList.push_back(T(i, i - state_vars.increVec(n),  dt * ( - firstCoefE/state_vars.dVec(n)  ) ));
                }*/
                /* Uncomment this section if you want second derivatives = constant  */
                //matList.push_back(T(i,i, - (1.0 - dt * A(i,0) )  ));   
                /*
                matList.push_back(T(i, i, - dt * (  secondCoefE / pow(state_vars.dVec(n), 2) ) ) );
                matList.push_back(T(i, i - state_vars.increVec(n), - dt * ( - 2 * secondCoefE / pow(state_vars.dVec(n), 2) ) ));
                matList.push_back(T(i, i - 2*state_vars.increVec(n), - dt * ( secondCoefE / pow(state_vars.dVec(n), 2) ) ) );
                */
                /*
                matList.push_back(T(i, i, - dt * (  1.0 / pow(state_vars.dVec(n), 2) ) ) );
                matList.push_back(T(i, i - state_vars.increVec(n), - dt * ( - 2 *  1.0 / pow(state_vars.dVec(n), 2) ) ));
                matList.push_back(T(i, i - 2*state_vars.increVec(n), - dt * ( 1.0 / pow(state_vars.dVec(n), 2) ) ) );
                */
            } else if ( std::abs(state_vars.stateMat(i,n) - state_vars.lowerLims(n)) < state_vars.dVec(n)/2.0 ) { //lower boundary
                
                atBoundIndicators(n) = 1.0;
                atBound = 1.0;

                ///* Uncomment this section if you want natural boundaries
                
                 matList.push_back(T(i, i, - dt * ( - firstCoefE/state_vars.dVec(n) + secondCoefE / pow(state_vars.dVec(n), 2) ) ) );
                 matList.push_back(T(i, i + state_vars.increVec(n), - dt * ( firstCoefE/state_vars.dVec(n) - 2 * secondCoefE / pow(state_vars.dVec(n), 2) ) ));
                 matList.push_back(T(i, i + 2*state_vars.increVec(n), - dt * ( secondCoefE / pow(state_vars.dVec(n), 2) ) ) );
                
                //*/
                /* Uncomment this section if you want first derivatives = constant
                 */
                // matList.push_back(T(i,i, - (1.0 - dt * A(i,0) ) ));
                // matList.push_back(T(i, i, - dt * ( - 1/state_vars.dVec(n)  ) ) );
                // matList.push_back(T(i, i + state_vars.increVec(n), - dt * ( 1/state_vars.dVec(n) ) ));
                 /*
                if ((n == 0) && atBoundIndicators(1) > 0 ) {
                 matList.push_back(T(i, i,  dt * ( - firstCoefE/state_vars.dVec(n)  ) ) );
                 matList.push_back(T(i, i + state_vars.increVec(n),  dt * ( firstCoefE/state_vars.dVec(n) ) ));
                }*/
                /* Uncomment this section if you want second derivatives = constant  */
                //matList.push_back(T(i,i, - (1.0 - dt * A(i,0) )/ state_vars.N ));
                /*
                matList.push_back(T(i, i, - dt * ( secondCoefE / pow(state_vars.dVec(n), 2) ) ) );
                matList.push_back(T(i, i + state_vars.increVec(n), - dt * (  - 2 * secondCoefE / pow(state_vars.dVec(n), 2) ) ));
                matList.push_back(T(i, i + 2*state_vars.increVec(n), - dt * ( secondCoefE / pow(state_vars.dVec(n), 2) ) ) );
                */
                /*
                    matList.push_back(T(i, i, - dt * ( 1.0 / pow(state_vars.dVec(n), 2) ) ) );
                    matList.push_back(T(i, i + state_vars.increVec(n), - dt * (  - 2 * 1.0 / pow(state_vars.dVec(n), 2) ) ));
                    matList.push_back(T(i, i + 2*state_vars.increVec(n), - dt * ( 1.0 / pow(state_vars.dVec(n), 2) ) ) );
                */
                
            }



        }
        
             
        if (atBound < 0 ) {
          // matList.push_back(T(i,i, (1.0 - dt * A(i,0))  ));
        }
        for (int n = (state_vars.N - 1); n >= 0; --n) {
            
            //add elements to the vector of triplets for matrix construction
            if ( atBoundIndicators(n) < 0) {
                double firstCoefE = B(i,n);
                double secondCoefE = C(i,n);
                
                //first derivative
                 matList.push_back(T(i,i, - dt * ( -firstCoefE * ( firstCoefE > 0) + firstCoefE * ( firstCoefE < 0) ) / state_vars.dVec(n)  ) );
                 matList.push_back(T(i,i + state_vars.increVec(n), - dt * firstCoefE * ( firstCoefE > 0) / state_vars.dVec(n) ));
                 matList.push_back(T(i,i - state_vars.increVec(n), - dt *  - firstCoefE * ( firstCoefE < 0) / state_vars.dVec(n) ));
                
                    matList.push_back(T(i, i, - dt * -2 * secondCoefE / ( pow(state_vars.dVec(n), 2) ) ));
                    matList.push_back(T(i, i + state_vars.increVec(n), - dt * secondCoefE / ( pow(state_vars.dVec(n), 2) ) ));
                    matList.push_back(T(i, i - state_vars.increVec(n), - dt * secondCoefE / ( pow(state_vars.dVec(n), 2) ) ));
                
            }

        }


    }
    //form matrices

    Le.setFromTriplets(matList.begin(), matList.end());

    //compress
    Le.makeCompressed(); 
}

py::tuple solveFT(Eigen::Ref<MatrixXdR> preLoadMat, Eigen::Ref<MatrixXdR> A, Eigen::Ref<MatrixXdR> B, Eigen::Ref<MatrixXdR> C,  Eigen::Ref<MatrixXdR> D, Eigen::Ref<MatrixXdR> v0, double dt, int tol)
{
    py::tuple data(3);
    stateVars stateSpace(preLoadMat);

    linearSysVars linearSys_vars(stateSpace, A,B,C,D,dt);
    linearSys_vars.constructMatFT(stateSpace);

    Eigen::VectorXd rhs; 

    rhs = v0.array() + dt * D.array(); // transform v0 into rhs
    /*********************************************/
    /* Change RHS to reflect boundary conditions */
    /*********************************************/

    //construct matrix
    /* uncomment this section if you want to set the boundary conditions to a constant */
    // for (int i = 0; i < stateSpace.S; ++i) {

    //     for (int n = (stateSpace.N - 1); n >=0; --n ) {
            
    //         //check whether it's at upper or lower boundary
    //         if ( std::abs(stateSpace.stateMat(i,n) - stateSpace.upperLims(n)) < stateSpace.dVec(n)/2 ) {  //upper boundary
    //          //   v0(i) = 0.0001;
    //         } else if ( std::abs( stateSpace.stateMat(i,n) - stateSpace.lowerLims(n)) < stateSpace.dVec(n)/2 ) { //lower boundary
    //          //v0(i) = 0.0001;            
    //         }
    //     }
    // }
     
    /* Initialize Eigen's cg solver */
 
    Eigen::VectorXd XiEVector;
    Eigen::LeastSquaresConjugateGradient<SpMat > cgE;
    // cgE.setMaxIterations(10000);
    cgE.setTolerance( pow(10,tol) );
    cgE.compute(linearSys_vars.Le);

    XiEVector = cgE.solveWithGuess(rhs, v0);
    data[0] = int(cgE.iterations());
    data[1] = cgE.error();
    data[2] = XiEVector;
    return data;    

}

py::tuple solveFK(Eigen::Ref<MatrixXdR> preLoadMat, Eigen::Ref<MatrixXdR> A, Eigen::Ref<MatrixXdR> B, Eigen::Ref<MatrixXdR> C,  Eigen::Ref<MatrixXdR> D, Eigen::Ref<MatrixXdR> v0, int iters)
{
    py::tuple data(3);
    stateVars stateSpace(preLoadMat);
    double dt(1.0);
    linearSysVars linearSys_vars(stateSpace, A,B,C,D,dt);
    linearSys_vars.constructMatFK(stateSpace);

    Eigen::VectorXd rhs;
    rhs =  dt * D.array(); // transform v0 into rhs
    /*********************************************/
    /* Change RHS to reflect boundary conditions */
    /*********************************************/

    //construct matrix
    /* uncomment this section if you want to set the boundary conditions to a constant */
    for (int i = 0; i < stateSpace.S; ++i) {

        for (int n = (stateSpace.N - 1); n >=0; --n ) {
            
            //check whether it's at upper or lower boundary
            if ( std::abs(stateSpace.stateMat(i,n) - stateSpace.upperLims(n)) < stateSpace.dVec(n)/2 ) {  //upper boundary
             //   v0(i) = 0.0001;
            } else if ( std::abs( stateSpace.stateMat(i,n) - stateSpace.lowerLims(n)) < stateSpace.dVec(n)/2 ) { //lower boundary
             //v0(i) = 0.0001;            
            }
        }
    }
     
    /* Initialize Eigen's cg solver */
    Eigen::VectorXd XiEVector;
    Eigen::LeastSquaresConjugateGradient<SpMat > cgE;
    cgE.setMaxIterations(iters);
    cgE.setTolerance( 0.000001 );
    cgE.compute(linearSys_vars.Le);  // update with Sparse matrix A
    XiEVector = cgE.solveWithGuess(rhs,v0);  // (rhs, guess)
    data[0] = int(cgE.iterations());
    data[1] = cgE.error();
    data[2] = XiEVector;

    return data;    

}
/*************************************/
/* Using pybind11 to interface       */
/* with python                       */
/*************************************/

PYBIND11_MODULE(SolveLinSys,m){
    m.doc() = "PDE Solver in cpp";

    m.def("solveFT", &solveFT, py::arg("stateSpace"),
        py::arg("A"), py::arg("B"), py::arg("C"), py::arg("D"),
        py::arg("v0"), py::arg("dt"), py::arg("tol"));

    m.def("solveFK", &solveFK, py::arg("stateSpace"),
        py::arg("A"), py::arg("B"), py::arg("C"), py::arg("D"),
        py::arg("v0"), py::arg("iters"));

}