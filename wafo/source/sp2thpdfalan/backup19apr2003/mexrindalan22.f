! This is a MEX-file for MATLAB.
! This file contains a mex-interface to RIND a subroutine
! for computing multivariate normal expectations.
! The file is self contained and should compile without errors on (Fortran90) 
! standard Fortran compilers.
!
! The mex-interface was written by
!			Per Andreas Brodtkorb
!			Norwegian Defence Research Establishment
!			P.O. Box 115
!			N-3191 Horten
!			Norway
!			Email: Per.Brodtkorb@ffi.no
!
!          Alan Genz
!          Department of Mathematics
!          Washington State University
!          Pullman, WA 99164-3113
!          Email : alangenz@wsu.edu
!
! 
! MEXRIND Computes multivariate normal expectations
!
!  E[Jacobian*Indicator|Condition ]*f_{Xc}(xc(:,ix)) 
!  where
!      "Indicator" = I{ H_lo(i) < X(i) < H_up(i), i=1:N_t+N_d }
!      "Jacobian"  = J(X(Nt+1),...,X(Nt+Nd+Nc)), special case is 
!      "Jacobian"  = |X(Nt+1)*...*X(Nt+Nd)|=|Xd(1)*Xd(2)..Xd(Nd)|
!      "condition" = Xc=xc(:,ix),  ix=1,...,Nx.
!      X = [Xt; Xd; Xc], a stochastic vector of Multivariate Gaussian 
!          variables where Xt,Xd and Xc have the length Nt, Nd and Nc,
!          respectively. (Recommended limitations Nx,Nt<=100, Nd<=6 and Nc<=10) 
! 
!  CALL: [value,error,inform]=mexrindalan22(S,m,indI,Blo,Bup,INFIN,xc,
!         Nt,SCIS,XcScale,ABSEPS,RELEPS,COVEPS,MAXPTS,MINPTS,seed);
!
!
!    VALUE  = estimated value for the expectation as explained above size 1 x Nx
!    ERROR  = estimated absolute error, with 99% confidence level.   size 1 x Nx
!   INFORM  = INTEGER, termination status parameter: (not implemented yet)
!            if INFORM = 0, normal completion with ERROR < EPS;
!            if INFORM = 1, completion with ERROR > EPS and MAXPTS 
!                           function vaules used; increase MAXPTS to 
!                           decrease ERROR;
!            if INFORM = 2, N > 100 or N < 1.
! 
!         S = Covariance matrix of X=[Xt;Xd;Xc] size Ntdc x Ntdc (Ntdc=Nt+Nd+Nc)
!         m = the expectation of X=[Xt;Xd;Xc]   size N x 1
!      indI = vector of indices to the different barriers in the  
!            indicator function,  length NI, where   NI = Nb+1 
!             (NB! restriction  indI(1)=0, indI(NI)=Nt+Nd )
! B_lo,B_up = Lower and upper barriers used to compute the integration 
!             limits, Hlo and Hup, respectively. size  Mb x Nb 
!    INFIN  = INTEGER, array of integration limits flags:  size 1 x Nb   (in)
!             if INFIN(I) < 0, Ith limits are (-infinity, infinity);
!             if INFIN(I) = 0, Ith limits are (-infinity, Hup(I)];
!             if INFIN(I) = 1, Ith limits are [Hlo(I), infinity);
!             if INFIN(I) = 2, Ith limits are [Hlo(I), Hup(I)].
!        xc = values to condition on            size Nc x Nx
!        Nt = size of Xt
!      SCIS = Integer defining integration method
!             1 Integrate all by SADAPT for Ndim<9 and by KRBVRC otherwise 
!             2 Integrate all by SADAPT by Genz (1992) (Fast)
!             3 Integrate all by KRBVRC by Genz (1993) (Fast)
!             4 Integrate all by KROBOV by Genz (1992) (Fast)
!             5 Integrate all by RCRUDE by Genz (1992)
!   XcScale = REAL to scale the conditinal probability density, i.e.,
!              f_{Xc} = exp(-0.5*Xc*inv(Sxc)*Xc + XcScale)
!    ABSEPS = REAL absolute error tolerance.
!    RELEPS = REAL relative error tolerance.
!    COVEPS = REAL error in cholesky factorization
!    MAXPTS = INTEGER, maximum number of function values allowed. This 
!             parameter can be used to limit the time. A sensible 
!             strategy is to start with MAXPTS = 1000*N, and then
!             increase MAXPTS if ERROR is too large.
!    MINPTS = INTEGER, minimum number of function values allowed
!    SEED   = INTEGER, seed to the random generator used in the integrations
!    NIT    = INTEGER, maximum number of Xt variables to integrate
!   xCutOff = REAL upper/lower truncation limit of the marginal normal CDF 
!
! 
!   If  Mb<Nc+1 then B_lo(Mb+1:Nc+1,:) is assumed to be zero.
!   The relation to the integration limits Hlo and Hup are as follows
!    IF INFIN(j)>=0,
!      IF INFIN(j)~=0,  Hlo(i)=Blo(1,j)+Blo(2:Mb,j).'*xc(1:Mb-1,ix), 
!      IF INFIN(j)~=1,  Hup(i)=Bup(1,j)+Bup(2:Mb,j).'*xc(1:Mb-1,ix), 
!
!   where i=indI(j-1)+1:indI(j), j=2:NI, ix=1:Nx
!
! This file was successfully compiled for matlab 5.3
! using Compaq Visual Fortran 6.1, and Windows 2000 and windows XP.
! The example here uses Fortran90 source.
! First, you will need to modify your mexopts.bat file.
! To find it, issue the command prefdir(1) from the Matlab command line,
! the directory it answers with will contain your mexopts.bat file.
! Open it for editing. The first section will look like:
!
!rem ********************************************************************
!rem General parameters
!rem ********************************************************************
!set MATLAB=%MATLAB%
!set DF_ROOT=C:\Program Files\Microsoft Visual Studio
!set VCDir=%DF_ROOT%\VC98
!set MSDevDir=%DF_ROOT%\Common\msdev98
!set DFDir=%DF_ROOT%\DF98
!set PATH=%MSDevDir%\bin;%DFDir%\BIN;%VCDir%\BIN;%PATH%
!set INCLUDE=%DFDir%\INCLUDE;%DFDir%\IMSL\INCLUDE;%INCLUDE%
!set LIB=%DFDir%\LIB;%VCDir%\LIB
!
! then you are ready to compile this file at the matlab prompt using the following command:
!     
!   mex -O mexrindalan22.f intmodule.f  jacobmod.f alanpab22.f 
!
!           

      SUBROUTINE mexFunction(nlhs, plhs, nrhs, prhs)
      USE GLOBALDATA
      USE RIND
      IMPLICIT NONE
C-----------------------------------------------------------------------
C     (integer) Replace integer by integer*8 on the DEC Alpha and the
C     SGI 64-bit platforms
C
      INTEGER :: plhs(*), prhs(*)
      INTEGER :: mxCreateFull, mxGetPr
      INTEGER :: S_pr, Ex_pr, IN_pr, BL_pr, BU_pr, INF_pr, Xc_pr
      INTEGER :: V_pr, E_pr     ! output pointers
C-----------------------------------------------------------------------
C
      INTEGER :: nlhs, nrhs
      INTEGER :: mxGetM, mxGetN, mxIsNumeric
      REAL*8  :: mxGetScalar
      INTEGER :: Nx,Nt,Nj,Nc,Ntd,Ntdc,Ni,Nb,Mb,K,I
C	DOUBLE PRECISION :: ABSEPS, RELEPS, VAL, ERR
      DOUBLE PRECISION :: XcScale
      INTEGER :: speed1, seed1
      DOUBLE PRECISION, ALLOCATABLE :: BIG(:,:),Xc(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: Blo(:,:),Bup(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: Ex(:),VALS(:),rINDI(:)
      DOUBLE PRECISION, ALLOCATABLE :: ERR(:), rINFI(:)
      INTEGER, ALLOCATABLE  :: IndI(:),INFIN(:)
      INTEGER, ALLOCATABLE  :: seed(:)
      INTEGER               :: seed_size
C     Check for proper number of arguments. 
      IF(nrhs .ne. 18) THEN
         CALL mexErrMsgTxt('18 inputs required.')
      ELSEIF(nlhs .lt.1 .or. nlhs .gt. 2 ) then
         CALL mexErrMsgTxt('1 or 2 outputs required.')
      END IF
C     Check to ensure the array is numeric (not strings).
      DO K=1,18
         IF(mxIsNumeric(prhs(K)) .EQ. 0 ) THEN
            CALL mexErrMsgTxt('Inputs must be numeric !')
         ENDIF
      END DO
C     Get the size of the input array.
      Ntdc = mxGetM(prhs(1))
      IF (Ntdc .NE. mxGetN(prhs(1)) .OR.
     &     Ntdc .NE. mxGetN(prhs(2))*mxGetM(prhs(2)) ) THEN
         CALL mexErrMsgTxt('Inconsistent size of Big and Ex !')
      END IF

      Nc = mxGetM(prhs(7))
      Nx = mxGetN(prhs(7))
      IF (Nx.LT.1) THEN
         Nx = 1;	
      END IF
      Ni = mxGetM(prhs(3))*mxGetN(prhs(3))
      Mb = mxGetM(prhs(4))
      Nb = mxGetN(prhs(4))
      IF (Nb+1.NE. Ni) THEN
         call mexErrMsgTxt('Nb+1 must equal Ni!')
      END IF
      IF (Mb .NE. mxGetM(prhs(5)) .or.
     &     Nb .NE. mxGetN(prhs(5))) THEN
         call mexErrMsgTxt('Size of Blo and Bup must be the same!')
      END IF
      IF (Nb .NE. mxGetM(prhs(6))*mxGetN(prhs(6))) THEN
         call mexErrMsgTxt('Size of INFI must equal Nb of Blo and Bup!')
      END IF
      Nt      = mxGetScalar(prhs(8))
! Set variables in the globaldata module
      SCIS    = mxGetScalar(prhs(9))
      XcScale = mxGetScalar(prhs(10))
      ABSEPS  = mxGetScalar(prhs(11))
      RELEPS  = mxGetScalar(prhs(12))
      EPS2    = mxGetScalar(prhs(13))
      EPS     = SQRT(EPS2)
      NSIMmax = mxGetScalar(prhs(14))
      NSIMmin = mxGetScalar(prhs(15))
      seed1   = mxGetScalar(prhs(16))
      NIT     = mxGetScalar(prhs(17))
      xCutOff = mxGetScalar(prhs(18))

      ALLOCATE(VALS(Nx))
      ALLOCATE(ERR(Nx))
      ALLOCATE(BIG(Ntdc,Ntdc))
      ALLOCATE(Ex(Ntdc))
      ALLOCATE(INDI(Ni))
      ALLOCATE(rINDI(Ni))
      ALLOCATE(Blo(Mb,Nb))
      ALLOCATE(Bup(Mb,Nb))
      ALLOCATE(rINFI(Nb))
      ALLOCATE(INFIN(Nb))
      ALLOCATE(Xc(Nc,Nx))
      
      IF (SCIS.LT.0) THEN
         SCIS = MAX(ABS(SCIS),1)
      ENDIF
      IF (SCIS.gt.0) THEN
         CALL random_seed(SIZE=seed_size) 
         ALLOCATE(seed(seed_size))
                                !print *,'rindinterface seed', seed1
         CALL random_seed(GET=seed(1:seed_size)) ! get current state
         seed(1:seed_size)=seed1 ! change seed
         CALL random_seed(PUT=seed(1:seed_size)) 
         CALL random_seed(GET=seed(1:seed_size)) ! get current state
                                !print *,'rindinterface seed', seed
         DEALLOCATE(seed)
      ENDIF

      Ntd = Ntdc - Nc;
!	Nd  = Ntd - Nt
	

C     Create matrix for the return argument.
      plhs(1) = mxCreateFull(1,Nx,0)
      V_pr    = mxGetPr(plhs(1))
      IF (nlhs.gt. 1) THEN
         plhs(2) = mxCreateFull(1,Nx,0)
         E_pr    = mxGetPr(plhs(2))
      ENDIF  


      S_pr   = mxGetPr(prhs(1))
      Ex_pr  = mxGetPr(prhs(2))
      IN_pr  = mxGetPr(prhs(3))
      Bl_pr  = mxGetPr(prhs(4))
      Bu_pr  = mxGetPr(prhs(5))
      INF_pr = mxGetPr(prhs(6))
      Xc_pr  = mxGetPr(prhs(7))

      call mxCopyPtrToReal8(S_pr,BIG,Ntdc*Ntdc)
      call mxCopyPtrToReal8(Ex_pr,Ex,Ntdc)
      call mxCopyPtrToReal8(IN_pr,rINDI,Ni)	
      call mxCopyPtrToReal8(Bl_pr,Blo,Mb*Nb)
      call mxCopyPtrToReal8(Bu_pr,Bup,Mb*Nb)
      call mxCopyPtrToReal8(INF_pr,rINFI,Nb)
      call mxCopyPtrToReal8(Xc_pr,Xc,Nc*Nx)
      DO K=1,Nb
         INDI(K)  = NINT(rINDI(K))
         INFIN(K) = NINT(rINFI(K))
      END DO
      INDI(Ni) = NINT(rINDI(Ni))
      
      DEALLOCATE(rINFI)
      DEALLOCATE(rINDI)
      
      IF (Ntd.EQ.INDI(NI)) THEN	
C     Call the computational subroutine.
         CALL RINDD(VALS,ERR,Big,Ex,Xc,Nt,INDI,Blo,Bup,INFIN,XcScale)
      ELSE
         DEALLOCATE(VALS)
         DEALLOCATE(ERR)
         DEALLOCATE(BIG)
         DEALLOCATE(Ex)
         DEALLOCATE(Xc)
         DEALLOCATE(INDI)
         DEALLOCATE(Blo)
         DEALLOCATE(Bup)
         DEALLOCATE(INFIN)
         CALL mexErrMsgTxt('INDI(Ni) must equal Nt+Nd !')
      ENDIF
      
C     Load the data into V_pr, which is the output to MATLAB
      CALL mxCopyReal8ToPtr(VALS,V_pr,Nx)   
      IF (nlhs.gt. 1) THEN
         CALL mxCopyReal8ToPtr(ERR,E_pr,Nx)
      ENDIF  
      
      DEALLOCATE(VALS)
      DEALLOCATE(ERR)
      DEALLOCATE(BIG)
      DEALLOCATE(Ex)
      DEALLOCATE(Xc)
      DEALLOCATE(INDI)
      DEALLOCATE(Blo)
      DEALLOCATE(Bup)
      DEALLOCATE(INFIN)
      RETURN
      END SUBROUTINE MEXFUNCTION

