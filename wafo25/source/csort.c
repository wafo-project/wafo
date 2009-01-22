#include "mex.h"
/*
 * csort - Counting sort
 *
 *
 * This is a MEX-file for MATLAB.
 * by Per Andreas Brodtkorb 15.08.2001
 */ 


int csort(double *a, double *len,double *b, double *ind, int n,int k,double amin)
{ int i,iy,ix;
/*CSORT Counting sort 
% b   = csort(a) sorted vector 
% ind = index vector defined by b = a(I);
% len = temporary work array
% a   = vector of integers to sort (length N).
% 
% CSORT assumes that each of the N input elements is an integer. Let K
% denote the range of the integers. When K = O(N) the sorting
% runs in O(N) time.*/
 for (i=0; i<n; i++) {   
   iy = ((int) (*(a + i)-amin));
   *(len+iy) =  *(len+iy) +1.0;
 }
 /* len(i) now contains the number of elements equal to i */
 for (i=1; i<k; i++) {
   *(len+i) =  *(len+i) + *(len+i-1);
 }
 /* len(i) now contains the number of elements less than or equal to i */
 /* Cumulative run lengths to indices */
 for (i=n; i>0; i--) {   
   iy = ((int) (*(a + i-1)-amin));
   ix = ((int) *(len + iy))-1;
   *(ind + ix) = i;
   *(b+ix) =  *(a + i-1);
   *(len+iy)=ix; 
 }
 return i;
}

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  mxArray *aptr, *aptr1;
  double *a,*len, *ind, *b, amin,amax;
  int length, mrows,ncols, ix,n,k;
  
  /*  check for proper number of arguments */
  if(nrhs<1) {
     mexErrMsgTxt("Not enough input arguments.");
  }
  if(nrhs>4) {
    mexErrMsgTxt("Too many input arguments.");
  }
  if(nlhs<1 || 2<nlhs) 
    mexErrMsgTxt("One or two outputs required.");
  
  /*  create a pointer to the input matrix a */
  a = mxGetPr(prhs[0]);

/*  get the dimensions of the matrix input a */
  mrows = mxGetM(prhs[0]);
  ncols = mxGetN(prhs[0]);
  if (ncols>mrows){
	  n = ncols;
  }
  else{
     n = mrows;
  }
  /*n  = max(ncols,mrows);*/
 
  if((ncols*mrows != n ) || (ncols*mrows == 0)){
    mexErrMsgTxt("Inputs must be row vectors or column vectors");
  }

  if (nrhs<2 || mxIsEmpty(prhs[1])){ 
    /* find minimum value*/
    amin = *(a+0);
    for (ix=1;ix<n;ix++){
      if (*(a+ix)<amin ) {
	amin = *(a+ix);
      }
    }
  }
  else{
    amin = mxGetScalar(prhs[1]);   
  }
  if (nrhs<3 || mxIsEmpty(prhs[2])){
    /* find maximum value*/
    amax = *(a+0);
    for (ix=1;ix<n;ix++) {
	if (*(a+ix)>amax ){
	  amax = *(a+ix);
	}
    }
  }
  else{
    amax = mxGetScalar(prhs[2]) ;
  }

  k = ((int) (amax-amin+1.0));
  /*  allocate enough memory to hold the resultant matrix */
  plhs[0]=mxCreateDoubleMatrix(n,1, mxREAL);
 
  aptr1 = mxCreateDoubleMatrix(k,1, mxREAL);
  /* and Create a pointers to the output matrix */
  b   =  mxGetPr(plhs[0]);  
  
  len =  mxGetPr(aptr1); 
  if (nlhs>1){
    plhs[1] = mxCreateDoubleMatrix(n,1, mxREAL);
    ind =  mxGetPr(plhs[1]); 
  }
  else{
    aptr = mxCreateDoubleMatrix(n,1, mxREAL);
    ind =  mxGetPr(aptr); 
  }

  /*  call the C subroutine */
  length = csort(a, len,b, ind, n,k,amin);

  /* Deallocate the old data */
  if (nlhs<2){
    mxDestroyArray(aptr);
  }
  mxDestroyArray(aptr1);

}
