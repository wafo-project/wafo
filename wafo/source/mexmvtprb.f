!* This is a MEX-file for MATLAB.
* and contains a mex-interface to Alan Genz's program,MVNDST a subroutine
* for computing  non-central multivariate t probabilities.
!     This subroutine uses an algorithm (QRSVN) described in the paper
!     "Methods for the Computation of Multivariate t-Probabilities",
!        by Alan Genz and Frank Bretz
!
!          Alan Genz 
!          Department of Mathematics
!          Washington State University 
!          Pullman, WA 99164-3113
!          Email : AlanGenz@wsu.edu
!
* The mex-interface was written by
*			Per Andreas Brodtkorb
*			Norwegian Defence Research Establishment
*			P.O. Box 115
*			N-3191 Horten
*			Norway
*			Email: Per.Brodtkorb@ffi.no
!
! CALL [value,error,inform] =
!               mexmvtprb(cov,A,constr,B,delta,nu,maxpts,abseps,releps)
!
!    COV   REAL NxN covariance matrix; only the diagonal and lower left 
!             triangle are used. COVRNC matrix must be positive semidefinite.
!     NU     INTEGER, number of degrees of freedom; if NU < 1, then
!             a multivariate normal computation is completed. 
!     M      INTEGER, number of linear constraints for integration region.
!     A       REAL size M array of lower integration limits.
!     CONSTR REAL MxN constraint matrix; integration region is all X with
!                     LOWER < DELTA + CONSTR*X < UPPER.
!     B      REAL, array of upper integration limits.
!     DELTA  REAL array of non-centrality parameters.
!     MAXPTS INTEGER, maximum number of function values allowed. This 
!            parameter can be used to limit the time. A sensible 
!            strategy is to start with MAXPTS = 1000*N, and then
!            increase MAXPTS if ERROR is too large.
!     ABSEPS REAL absolute error tolerance.
!     RELEPS REAL relative error tolerance.
!     ERROR  REAL estimated absolute error, with 99% confidence level.
!     VALUE  REAL estimated value for the integral
!     INFORM INTEGER, termination status parameter:
!            if INFORM = 0, normal completion with ERROR < EPS;
!            if INFORM = 1, completion with ERROR > EPS and MAXPTS 
!                           function vaules used; increase MAXPTS to 
!                           decrease ERROR;
!            if INFORM = 2, N > NL or N < 1.
!            if INFORM = 3, covariance matrix not positive semidefinite.
!


* This file was successfully compiled for matlab 5.3, 6.1
* using Compaq Visual Fortran 6.1, and Windows 2000, XP.
* The example here uses Fortran90 source.
* First, you will need to modify your mexopts.bat file.
* To find it, issue the command prefdir(1) from the Matlab command line,
* the directory it answers with will contain your mexopts.bat file.
* Open it for editing. The first section will look like:
*
*rem ********************************************************************
*rem General parameters
*rem ********************************************************************
*set MATLAB=%MATLAB%
*set DF_ROOT=C:\Program Files\Microsoft Visual Studio
*set VCDir=%DF_ROOT%\VC98
*set MSDevDir=%DF_ROOT%\Common\msdev98
*set DFDir=%DF_ROOT%\DF98
*set PATH=%MSDevDir%\bin;%DFDir%\BIN;%VCDir%\BIN;%PATH%
*set INCLUDE=%DFDir%\INCLUDE;%DFDir%\IMSL\INCLUDE;%INCLUDE%
*set LIB=%DFDir%\LIB;%VCDir%\LIB
*
* Now, basically the problem is that mex.bat doesn't recognize files
* with extension .f90 as being source files.
! To fix it do:
! 
!  for matlab R11: 
!     (1) open c:\matlab\bin\mex.bat using an editor 
!     (2) line 2338 is:
!           if ($EXTENSION =~ /(c|f|cc|cxx|cpp|for)$/i ) {
!         change it to
!           if ($EXTENSION =~ /(c|f|cc|cxx|cpp|for|f90)$/i ) {
!         (notice the "|f90"?)
!
! for matlab R12>:
!    make the changes as indicated above to either your 
!      c:\matlab\bin\mex.bat
!    or
!      c:\matlab\bin\mex.pl
* 
* then you are ready to compile this file at the matlab prompt using the
* following command:
!
*  mex -O -output mexmvtprb  mvdist.f90 mexmvtprb.f
      subroutine mexFunction(nlhs, plhs, nrhs, prhs)
      USE PRECISION_MODEL
      USE MVSTAT
      IMPLICIT NONE
C-----------------------------------------------------------------------
C     (integer) Replace integer by integer*8 on the DEC Alpha and the
C     SGI 64-bit platforms
C
      integer plhs(*), prhs(*)
      integer mxCreateFull, mxGetPr
      integer V_pr, E_pr, IF_pr, C_pr, A_pr, CT_pr,B_pr,D_pr
C-----------------------------------------------------------------------
C
      integer nlhs, nrhs
      integer mxGetM, mxGetN, mxIsNumeric
      real*8  mxGetScalar
      real*8 :: rIFT
      INTEGER :: K,N1,N3,M1,M2,M3,M4,M5
      INTEGER :: NU, MAXPTS, IFT,IVLS
      REAL(KIND=STND) :: ABSEPS, RELEPS, VAL, ERR, INFINITY     
      PARAMETER (INFINITY = 37.0_STND)
      REAL(KIND=STND) ,ALLOCATABLE,  DIMENSION(:,:)  :: COV,CSTR 
      REAL(KIND=STND) ,ALLOCATABLE,  DIMENSION(:)  :: LOW, UP, D 
      INTEGER, ALLOCATABLE, DIMENSION(:) :: INFIN
C     Check for proper number of arguments. 
      if(nrhs .ne. 9) then
         call mexErrMsgTxt('9 inputs required.')
      elseif(nlhs .ne. 3) then
         call mexErrMsgTxt('3 outputs required.')
      endif
      DO K=1,9
         if(.not.mxIsNumeric(prhs(K))) THEN
            call mexErrMsgTxt('Inputs must be numeric!')
         endif
      ENDDO

C     Get the size of the input array.
      M1 = mxGetM(prhs(1))
      N1 = mxGetN(prhs(1))
C     Column * row should be smaller than 5*990
      if(M1.NE.N1) then
         call mexErrMsgTxt('#Rows and #columns must be equal.')
      endif
      M2 = mxGetN(prhs(2))*mxGetM(prhs(2))
      N3 = mxGetN(prhs(3))
      M3 = mxGetM(prhs(3))
      M4 = mxGetN(prhs(4))*mxGetM(prhs(4))
      M5 = mxGetN(prhs(5))*mxGetM(prhs(5))
      if(M2.ne.M3 .or. M3.ne.M4 .or.M4.ne.M5.OR.N3.NE.N1) then
      call mexErrMsgTxt('Inconsistent size of Cov, A, CSTR, B and d.')
      endif      
      
      NU      = mxGetScalar(prhs(6))
      ABSEPS  = mxGetScalar(prhs(7))
      RELEPS  = mxGetScalar(prhs(8))
      MAXPTS  = mxGetScalar(prhs(9))

C     Create matrix for the return argument.
      plhs(1) = mxCreateFull(1,1,0)
      plhs(2) = mxCreateFull(1,1,0)
      plhs(3) = mxCreateFull(1,1,0)

      V_pr  = mxGetPr(plhs(1))
      E_pr  = mxGetPr(plhs(2))
      IF_pr = mxGetPr(plhs(3))

      ALLOCATE(COV(N1,N1))
      ALLOCATE(LOW(M2))
      ALLOCATE(CSTR(M2,N1))
      ALLOCATE(UP(M2))
      ALLOCATE(D(M2))
      ALLOCATE(INFIN(M2))

      C_pr = mxGetPr(prhs(1))
      A_pr = mxGetPr(prhs(2))
      CT_pr = mxGetPr(prhs(3))
      B_pr = mxGetPr(prhs(4))
      D_pr = mxGetPr(prhs(5))
      call mxCopyPtrToReal8(C_pr,COV,N1*N1)
      call mxCopyPtrToReal8(A_pr,LOW,M2)
      call mxCopyPtrToReal8(CT_pr,CSTR,M2*N1)
      call mxCopyPtrToReal8(B_pr,UP,M2)
      call mxCopyPtrToReal8(D_pr,D,M2)
*     
* Set INFIN  INTEGER, array of integration limits flags:
*            if INFIN(I) < 0, Ith limits are (-infinity, infinity);
*            if INFIN(I) = 0, Ith limits are (-infinity, UPPER(I)];
*            if INFIN(I) = 1, Ith limits are [LOWER(I), infinity);
*            if INFIN(I) = 2, Ith limits are [LOWER(I), UPPER(I)].		

	DO k = 1,M2
           INFIN(K) = 2
           if (LOW(K).LE.-INFINITY) THEN
              if (UP(K) .GE. INFINITY) THEN
                 INFIN(K) = -1
              else
                 INFIN(K) = 0
              endif
           else if (UP(K).GE.INFINITY) THEN 
              INFIN(K) = 1
           endif	
	ENDDO

C	call mxCopyPtrToInteger4(M_pr,MAXPTS,1)
      CALL MVDIST( N1, COV, NU, M2, LOW, CSTR, UP, INFIN, D,
     &     MAXPTS, ABSEPS, RELEPS, ERR, VAL, IVLS, IFT )
      
      rIFT = IFT
C     Load the data into V_pr, which is the output to MATLAB
      call mxCopyReal8ToPtr(VAL,V_pr,1)     
      call mxCopyReal8ToPtr(ERR,E_pr,1)
      call mxCopyReal8ToPtr(rIFT,IF_pr,1)
      DEALLOCATE(COV)
      DEALLOCATE(LOW)
      DEALLOCATE(CSTR)
      DEALLOCATE(UP)
      DEALLOCATE(D)
      DEALLOCATE(INFIN)
      RETURN
      END SUBROUTINE MEXFUNCTION
