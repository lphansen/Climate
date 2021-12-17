#include "petsclinearsystem.h"

static inline void fill_mat_values(PetscScalar *State, PetscInt i, PetscInt center, PetscInt j, PetscScalar *lowerLims, PetscScalar *upperLims, PetscScalar *dVec, PetscInt *incVec, PetscInt maxcols,PetscScalar *B, PetscScalar *C, PetscScalar dt, PetscInt *cols, PetscScalar *vals)
{ 
  PetscScalar firstCoefE = B[i],secondCoefE = C[i];

  //check whether it's at upper or lower boundary
  if (PetscAbs(State[i]-upperLims[j]) < dVec[j]/2.0) {//upper boundary
    vals[center] += -dt*(firstCoefE/dVec[j] + secondCoefE/PetscPowReal(dVec[j],2));
    vals[center-1-2*j] = -dt*(-firstCoefE/dVec[j] - 2.0*secondCoefE/PetscPowReal(dVec[j],2));
    vals[center-2-2*j] = -dt*(secondCoefE/PetscPowReal(dVec[j],2));
    cols[center-1-2*j] = i - incVec[j];
    cols[center-2-2*j] = i - 2*incVec[j];
  } else if (PetscAbs(State[i]-lowerLims[j]) < dVec[j]/2.0) {//lower boundary
    vals[center] += -dt*(-firstCoefE/dVec[j] + secondCoefE/PetscPowReal(dVec[j],2));
    vals[center+1+2*j] = -dt*(firstCoefE/dVec[j] - 2.0*secondCoefE/PetscPowReal(dVec[j],2));
    vals[center+2+2*j] = -dt*(secondCoefE/PetscPowReal(dVec[j],2));
    if (i + incVec[j] < maxcols) cols[center+1+2*j] = i + incVec[j]; // ignore out of bound entries
    if (i + 2*incVec[j] < maxcols) cols[center+2+2*j] = i + 2*incVec[j];
  } else {
    //first derivative
    vals[center] += -dt*(-firstCoefE*(firstCoefE > 0) + firstCoefE*( firstCoefE<0))/dVec[j] - dt*(-2)*secondCoefE/(PetscPowReal(dVec[j],2));
    vals[center+1+2*j] = -dt*firstCoefE*(firstCoefE > 0)/dVec[j] - dt* secondCoefE/(PetscPowReal(dVec[j],2));
    vals[center-1-2*j] = -dt*-firstCoefE*(firstCoefE < 0)/dVec[j] - dt*secondCoefE/(PetscPowReal(dVec[j],2));
    cols[center-1-2*j] = i - incVec[j];
    if (i + incVec[j] < maxcols) cols[center+1+2*j] = i + incVec[j]; // ignore out of bound entries
  }
}

PetscErrorCode FormLinearSystem_C(PetscScalar *R, PetscScalar *F, PetscScalar *K, PetscScalar *A, PetscScalar *B_r, PetscScalar *B_f, PetscScalar *B_k, PetscScalar *C_rr, PetscScalar *C_ff, PetscScalar *C_kk, PetscScalar dt, PetscScalar *lowerLims, PetscScalar *upperLims, PetscScalar *dVec, PetscInt *incVec, PetscInt n, Mat petsc_mat)
{
  PetscErrorCode ierr;
  PetscInt       i, center;
  PetscInt       cols[13];
  PetscScalar    vals[13];

  PetscFunctionBegin;
  for (i = 0; i < n; ++i) {
    center = 3*4/2;
    memset(vals,0,13*sizeof(PetscScalar));
    memset(cols,-1,13*sizeof(PetscInt));
    cols[center] = i;
    vals[center] = 1.0 - dt * A[i];
    fill_mat_values(R,i,center,0,lowerLims,upperLims,dVec,incVec,n,B_r,C_rr,dt,cols,vals);
    fill_mat_values(F,i,center,1,lowerLims,upperLims,dVec,incVec,n,B_f,C_ff,dt,cols,vals);
    fill_mat_values(K,i,center,2,lowerLims,upperLims,dVec,incVec,n,B_k,C_kk,dt,cols,vals);
    ierr = MatSetValues(petsc_mat,1,&i,3*4+1,cols,vals,INSERT_VALUES);CHKERRQ(ierr);
  }
  ierr = MatAssemblyBegin(petsc_mat,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(petsc_mat,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
