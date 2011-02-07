#include "mex.h"

/*
 * findcross.c - 
 *
 *  Returns indices to level l crossings of argument vector
 *
 * This is a MEX-file for MATLAB.
 * 1998 by Per Andreas Brodtkorb. last modified 23.06-98
 */ 




int findcross(double *y, double h, double *ind, int n)
{ int i,start, ix=0,dcross=0;

 if  ( *(y +0)< h){
    dcross=-1; /* first is a up-crossing*/ 
 }
 if  ( *(y +0)> h){
    dcross=1;  /* first is a down-crossing*/ 
 }
 start=0;
 if  ( *(y +0)== h){
    /* Find out what type of crossing we have next time.. */
    for (i=1; i<n; i++) {
       start=i;
       if  ( *(y +i)< h){
	  *(ind + ix) = i; /* first crossing is a down crossing*/ 
	  ix++; 
	  dcross=-1; /* The next crossing is a up-crossing*/ 
	  break;
       }
       if  ( *(y +i)> h){
	  *(ind + ix) = i; /* first crossing is a up-crossing*/ 
	  ix++; 
	  dcross=1;  /*The next crossing is a down-crossing*/ 
	  break;
       }
    }
 }
 
 for (i=start; i<n-1; i++) {
    if (( (dcross==-1) && (*(y +i)<=h) && (*(y+i+1) > h)  )  || ((dcross==1 ) && (*(y +i)>=h) && (*(y+i+1) < h) ) )  { 
      
       *(ind + ix) = i+1 ;
       ix++;
       dcross=-dcross;
    }
    
    
    
 }
 return ix;
}

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  mxArray *array_ptr;
  double *y, *ind;
  double  h;
  int length, status,mrows,ncols, ix;
  
  /*  check for proper number of arguments */
  if(nrhs>2) 
     mexErrMsgTxt("One or two inputs required.");
  if(nlhs!=1) 
    mexErrMsgTxt("One output required.");

  /*  create a pointer to the input matrix y */
  y = mxGetPr(prhs[0]);
  if(nrhs==2){ 
     /* check to make sure the second input argument is a scalar */
     if( !mxIsNumeric(prhs[1]) || !mxIsDouble(prhs[1]) ||
	 mxIsEmpty(prhs[1])    || mxIsComplex(prhs[1]) ||
	 mxGetN(prhs[1])*mxGetM(prhs[1])!=1 ) {
	mexErrMsgTxt("Input x must be a scalar.");
     } 
  
     /*  get the scalar input x */
     h = mxGetScalar(prhs[1]); }
  else
     h=0.0;
  
  
  /*  get the dimensions of the matrix input y */
  mrows = mxGetM(prhs[0]);
  ncols = mxGetN(prhs[0]);
  
  if(ncols<mrows){
     if( ncols!=1  )
	mexErrMsgTxt("Inputs must be row vectors or column vectors");
  }
  else {
  if( mrows!=1  )
    mexErrMsgTxt("Inputs must be row vectors or column vectors");
  mrows=ncols;
  }
    
  /*  set pointer to the resultant matrix */
  array_ptr=mxCreateDoubleMatrix((mrows-1),1, mxREAL);
  /* allocate enough memory to hold the indices
    ie. make a copy of the output matrix */
  ind=  mxGetPr(array_ptr);   
  
  /*  call the C subroutine */
  length = findcross(y,h,ind,mrows);

  if (length==0){
    /* returning the empty matrix, no TP found  */
     plhs[0]=mxCreateDoubleMatrix(0,0, mxREAL);
     mxSetN(plhs[0],0); 
     mxSetM(plhs[0],0); 
     /* Deallocate the old data */
     mxDestroyArray(array_ptr);
  }
  else{ 
     mxSetM(array_ptr,length);
      /* Redo ind to respond to the new size */ 
     plhs[0]=mxCreateDoubleMatrix(length,1, mxREAL);
     /* copy the elements */
     memcpy(mxGetPr(plhs[0]),ind, length*(sizeof(double)));
     /* Deallocate the old data */
     mxDestroyArray(array_ptr);
  }


}
