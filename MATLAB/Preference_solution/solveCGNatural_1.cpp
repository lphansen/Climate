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

//include Eigen
#include "eigen3/Eigen/Dense"
#include "eigen3/Eigen/Sparse"
#include "eigen3/Eigen/SparseLU"
#include "eigen3/Eigen/SparseQR"
#include "eigen3/Eigen/SparseCholesky"
#include "eigen3/Eigen/IterativeLinearSolvers"
typedef Eigen::SparseMatrix<double > SpMat;
typedef Eigen::Triplet<double> T;
#include <typeinfo>
#include "MarketIO.h"
//include MEX related files
#include "mex.h"
#include "matrix.h"

/* ******************** */
/* State Variable Class */
/* ******************** */

Eigen::MatrixXd empty;
Eigen::ArrayXd emptyAry;
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
    SpMat Le;

    //member functions
    
    //constructor
    linearSysVars(stateVars & state_vars, Eigen::MatrixXd A, Eigen::MatrixXd B, Eigen::MatrixXd C, Eigen::MatrixXd D, double dt);
    
    //function to construt matrix
    
    void constructMat(stateVars & state_vars);
    
};

linearSysVars::linearSysVars(stateVars & state_vars, Eigen::MatrixXd AInput, Eigen::MatrixXd BInput, Eigen::MatrixXd CInput, Eigen::MatrixXd DInput, double dtInput) {
        
    Le.resize(state_vars.S,state_vars.S);
    A.resize(state_vars.S,1); B.resize(state_vars.S,state_vars.N);
    C.resize(state_vars.S,state_vars.N); D.resize(state_vars.S,1);
    A = AInput; B = BInput; C = CInput; D = DInput;
    dt = dtInput;

}



void linearSysVars::constructMat(stateVars & state_vars) {
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
        
        //matList.push_back(T(i,i, (1.0 - dt * A(i,0))  ));
        matList.push_back(T(i,i, (0.0 - dt * A(i,0))  ));
        
        for (int n = (state_vars.N - 1); n >=0; --n ) {
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
                
                /* Uncomment this section if you want first derivatives = constant  */
//                 matList.push_back(T(i,i, - (1.0 - dt * A(i,0) ) ));
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
//                 matList.push_back(T(i,i, - (1.0 - dt * A(i,0) ) ));
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
//     saveMarket(Le,"Le_local_dt.dat");

}


void
mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
    int        i;

    
    /* Translating state space inputs into C++ EIGEN */
    
    /* Create state space */

    int nRows = mxGetM(prhs[0]); int nCols = mxGetN(prhs[0]);
    Eigen::Map<Eigen::MatrixXd> preLoadMat((double *)mxGetPr(prhs[0]),nRows,nCols);
    
    stateVars stateSpace(preLoadMat);    


    /* Load matrix coefficients */
    mxArray *mxValue; 
    mxValue = mxGetField(prhs[1], 0, "A");
    nRows = mxGetM(mxValue); nCols = mxGetN(mxValue);
    Eigen::Map<Eigen::MatrixXd> A((double *)mxGetPr(mxValue),nRows,nCols);

    
    mxValue = mxGetField(prhs[1], 0, "B");      
    nRows = mxGetM(mxValue); nCols = mxGetN(mxValue);
    Eigen::Map<Eigen::MatrixXd> B((double *)mxGetPr(mxValue),nRows,nCols);
    
    mxValue = mxGetField(prhs[1], 0, "C");      
    nRows = mxGetM(mxValue); nCols = mxGetN(mxValue);
    Eigen::Map<Eigen::MatrixXd> C((double *)mxGetPr(mxValue),nRows,nCols);
    
    mxValue = mxGetField(prhs[1], 0, "D");      
    nRows = mxGetM(mxValue); nCols = mxGetN(mxValue);
    Eigen::Map<Eigen::MatrixXd> D((double *)mxGetPr(mxValue),nRows,nCols);

    mxValue = mxGetField(prhs[1], 0, "v0");      
    nRows = mxGetM(mxValue); nCols = mxGetN(mxValue);
    Eigen::Map<Eigen::MatrixXd> v0((double *)mxGetPr(mxValue),nRows,nCols);
    
    mxValue = mxGetField(prhs[1], 0, "v1");      
    nRows = mxGetM(mxValue); nCols = mxGetN(mxValue);
    Eigen::Map<Eigen::MatrixXd> v1((double *)mxGetPr(mxValue),nRows,nCols);

    /* Construct linear system */
    mxValue = mxGetField(prhs[1], 0, "dt"); 
    double dt = mxGetPr(mxValue)[0];  //Load in dt;
    //saveMarket(dt,"Le_local_dt.dat");
    linearSysVars linearSys_vars(stateSpace, A,B,C,D,dt);
    
    linearSys_vars.constructMat(stateSpace);
    
    // Eigen::VectorXd v1; v1.resize(stateSpace.S, stateSpace.N); v1 = v0; // smart guess
    v0 = v0.array() + dt * D.array(); // transform v0 into rhs
    
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
     
    /***************************************/
    /* Solve the system and send to MATLAB */  
    /***************************************/
    
    /* Initialize Eigen's cg solver */
    
    Eigen::VectorXd XiEVector;
    Eigen::LeastSquaresConjugateGradient<SpMat > cgE;
    cgE.setMaxIterations(400000);
    cgE.setTolerance( 0.000001 );
    cgE.compute(linearSys_vars.Le);
    XiEVector = cgE.solveWithGuess(v0,v1);
    mexPrintf("CONJUGATE GRADIENT TOOK (number of iterations): %3i% \n",  cgE.iterations() );
    mexPrintf("CONJUGATE GRADIENT error: %3f% \n",  cgE.error());
    mexEvalString("drawnow;");
    v1 = XiEVector;
    
    /* Try to use LU decomp*/
    /*
    Eigen::SparseLU<SpMat > solver; 
    solver.analyzePattern(linearSys_vars.Le);
    solver.factorize(linearSys_vars.Le);
    v1 = solver.solve(v0);
    */
    /* Send solution to MATLAB */
    
    mwSize rows = stateSpace.S;
    mwSize cols = 1;
    
    //mexPrintf("Conjugate gradient took iterations:%3i%\n", cgE.iterations());

    plhs[0] = mxCreateDoubleMatrix(rows, cols, mxREAL); // Create MATLAB array of same size
    Eigen::Map<Eigen::MatrixXd> map(mxGetPr(plhs[0]), rows, cols); // Map the array
    map = v1;
    
//     plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL); // Create MATLAB array of same size
//     Eigen::Map<Eigen::MatrixXd> map1(mxGetPr(plhs[1]), 1, 1); // Map the array
//     map1 = cgE.iterations();
//     
//     plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL); // Create MATLAB array of same size
//     Eigen::Map<Eigen::MatrixXd> map2(mxGetPr(plhs[2]), 1, 1); // Map the array
//     map2 = cgE.error();
//     
}

