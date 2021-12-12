from petsc4py.PETSc cimport Mat, PetscMat

from petsc4py.PETSc import Error
cimport numpy as np # this allows access to the data member

cdef extern from "petsclinearsystem.h":
    ctypedef struct Params:
        double lambda_
    int FormLinearSystem_C(double *R, double *F, double *K, double *A, double *B_r, double *B_f_k, double *B_k, double *C_rr, double *C_ff, double *C_kk, double dt, double *lowerLims, double *upperLims, double *dVec, int *incVec, int n, PetscMat petsc_mat)

def formLinearSystem(np.ndarray R, np.ndarray F, np.ndarray K, np.ndarray A, np.ndarray B_r, np.ndarray B_f, np.ndarray B_k, np.ndarray C_rr, np.ndarray C_ff, np.ndarray C_kk, double dt, np.ndarray lowerLims, np.ndarray upperLims, np.ndarray dVec, np.ndarray incVec, Mat pymat):
    cdef int ierr = 0
    cdef int n = len(A)

    ierr = FormLinearSystem_C(<double*>R.data,<double*>F.data, <double*>K.data, <double*>A.data, <double*>B_r.data, <double*>B_f.data, <double*>B_k.data, <double*>C_rr.data, <double*>C_ff.data, <double*>C_kk.data, dt, <double*>lowerLims.data, <double*>upperLims.data, <double*>dVec.data, <int*>incVec.data, n, pymat.mat)
    if ierr != 0: raise Error(ierr)
    return pymat