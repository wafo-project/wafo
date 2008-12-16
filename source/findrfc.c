#include "mex.h"
#include "math.h"
/*
 * findrfc.c - 
 *
 *  Returns indices to RFC turningpoints of a vector
 *  of turningpoints
 *
 * This is a MEX-file for MATLAB.
 * 1998 by Per Andreas Brodtkorb.
 */

int rfc_filt(double *y1,double hmin, double *ind, int n) {
   double xminus,xplus,Tpl,Tmi,*y,Tstart;
   int i,j,ix=0,NC,iy;
   
   

   if (*(y1+0)> *(y1+1)){ /* if first is a max*/
      y=&(*(y1+1));  /* ignore the first max*/
      NC=floor((n-1)/2);
      Tstart=2;
   }
   else {
      y=y1;
      NC=floor(n/2);
      Tstart=1;
   }
    
   if (NC<1){
      return ix; /* No RFC cycles*/
   }
   

   if (( *(y+0) > *(y+1)) && ( *(y+1) > *(y+2)) ){
      ix=-1;
      return ix; /*This is not a sequence of turningpoints, exit */
   }
   if ((*(y+0) < *(y+1)) && (*(y+1)< *(y+2))){
      ix=-1;
      return ix; /*This is not a sequence of turningpoints, exit */
   }
   
   
   for (i=0; i<NC; i++) {
      
      Tmi=Tstart+2*i;
      Tpl=Tstart+2*i+2;
      xminus=*(y+2*i);
      xplus=*(y+2*i+2);

      if(i!=0){
	 j=i-1;
	 while((j>=0) && (*(y+2*j+1)<=*(y+2*i+1))){
	    if( (*(y+2*j)<xminus) ){
	       xminus=*(y+2*j);
	       Tmi=Tstart+2*j;
	    } /*if */
	    j--;
	 } /*while j*/
      } /*if i */
      if ( xminus >= xplus){
	 if ( (*(y+2*i+1)-xminus) >= hmin){
	    *(ind+ix)=Tmi;
	    ix++;
	    *(ind+ix)=(Tstart+2*i+1);
	    ix++;
	 } /*if*/
	 goto L180;
      } 
      
      j=i+1;
      while((j<NC) ) {
	 if (*(y+2*j+1) >= *(y+2*i+1)) goto L170;
	 if( (*(y+2*j+2) <= xplus) ){
	    xplus=*(y+2*j+2);
	    Tpl=(Tstart+2*j+2);
	 }/*if*/
	    j++;
      } /*while*/

      
      if ( (*(y+2*i+1)-xminus) >= hmin) {
	 *(ind+ix)=Tmi;
	 ix++;	
	 *(ind+ix)=(Tstart+2*i+1);
	 ix++;
	 
      } /*if*/
      goto L180;
   L170: 
      if (xplus <= xminus ) {
	 if ( (*(y+2*i+1)-xminus) >= hmin){
	    *(ind+ix)=Tmi;
	    ix++;
	    *(ind+ix)=(Tstart+2*i+1);
	    ix++;
	 } /*if*/
	 /*goto L180;*/
      }
      else{	    
	 if ( (*(y+2*i+1)-xplus) >= hmin) {
	    *(ind+ix)=(Tstart+2*i+1);
	    ix++;	
	    *(ind+ix)=Tpl;
	    ix++;	
	 } /*if*/
      } /*elseif*/
   L180:
     iy=i;
   }  /* for i */
  return ix;
}

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  mxArray *output_array[1],*array_ptr;
  double *y, *ind,*tmp;
  double  h;
  int length, mrows,ncols;
  
  /*  check for proper number of arguments */
  if(nrhs!=2) 
     mexErrMsgTxt("Two inputs required.");
  if(nlhs!=1) 
    mexErrMsgTxt("One output required.");
  
  
  /* check to make sure the second input argument is a scalar */
  if( !mxIsNumeric(prhs[1]) || !mxIsDouble(prhs[1]) ||
      mxIsEmpty(prhs[1])    || mxIsComplex(prhs[1]) ||
      mxGetN(prhs[1])*mxGetM(prhs[1])!=1 ) {
     mexErrMsgTxt("Input x must be a scalar.");
  } 
  
  
  /*  get the scalar input x */
  h = mxGetScalar(prhs[1]);
  if( h<=0) {
     mexErrMsgTxt("Input h must be larger than zero.");
  }
  
  /*  get the dimensions of the matrix input  */
  mrows = mxGetM(prhs[0]);
  ncols = mxGetN(prhs[0]);

  if(ncols<mrows){
     if( ncols!=1  )
	mexErrMsgTxt("Input must be row vectors or column vectors");
  }
  else {
     if( mrows!=1  )
	mexErrMsgTxt("Input must be row vectors or column vectors");
     mrows=ncols;
  }

  /*  create a pointer to the input matrix y*//* temporary matrix*/ 
  tmp = mxGetPr(prhs[0]);
  
 
  y= (double *)mxCalloc(mrows+1, sizeof(double));
  /* copy and reshape the y matrix so that the last element is zero*/
  memcpy(y,(const void *)tmp, mrows*(sizeof(double)));
  
  array_ptr=mxCreateDoubleMatrix((mrows+1),1, mxREAL);
  ind=  mxGetPr(array_ptr); 
  length = rfc_filt(y,h,ind,mrows);
  mxFree(y);


  if (length==-1){
     mexErrMsgTxt("The input is not a sequence of turningpoints");
  }else{
     if (length<0){
	mexErrMsgTxt("Something is wrong! I don't know what.");
     }
  }

  
  if (length==0){
    /* returning the empty matrix, no rfc cycles found  */
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
     /* Sorting the indices*/
     mexCallMATLAB(1,output_array, 1,plhs,"sort");
      /*  load the new matrix data into plhs[0] */
    memcpy(mxGetPr(plhs[0]),mxGetPr(output_array[0]), length*(sizeof(double)));
    mxDestroyArray(output_array[0]);
  }
   
}



