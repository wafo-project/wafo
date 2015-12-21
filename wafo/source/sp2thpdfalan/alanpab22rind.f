!****************************************************************************
!NB! if compilation complains about too many continuation lines extend it.
! 
!
!  modules:   GLOBALDATA, GLOBCONST, FUNCMOD , RIND    Version 1.0   
!
! Programs available in module RIND : 
! (NB! the GLOBALDATA  module is also used to transport the inputs)  
!
!
! SETDATA initializes global constants explicitly:
!   
! CALL SETDATA(RELEPS,ABSEPS,EPS2,xCutOff) 
!
!                   GLOBALDATA module :
!   RELEPS,ABSEPS = relative and absolute requested error 
!        EPS2  = if conditional variance is less it is considered as zero
!                i.e., the variable is considered deterministic 
!      xCutOff = 5 (standard deviations by default)
!     FxCutOff = FI(xCutOff) - FI( -xCutOff), truncation of the normal CDF.  
!
! INITDATA initializes global constants implicitly:
!
! CALL INITDATA (speed)    
!
!        speed = 1,2,...,9 (1=slowest and most accurate,9=fastest, 
!                           but less accurate)
!
! see the GLOBALDATA, FUNCMOD, GLOBALCONST 
!
!
!RINDD computes  E[Jacobian*Indicator|Condition]*f_{Xc}(xc(:,ix)) 
!
! where
!     "Indicator" = I{ H_lo(i) < X(i) < H_up(i), i=1:Nt+Nd }
!     "Jacobian"  = J(X(Nt+1),...,X(Nt+Nd+Nc)), special case is 
!     "Jacobian"  = |X(Nt+1)*...*X(Nt+Nd)|=|Xd(1)*Xd(2)..Xd(Nd)|
!     "condition" = Xc=xc(:,ix),  ix=1,...,Nx.
!     X = [Xt; Xd ;Xc], a stochastic vector of Multivariate Gaussian 
!         variables where Xt,Xd and Xc have the length Nt, Nd and Nc,
!         respectively. 
!         (Recommended limitations Nx, Nt<101, Nd<7 and NIT,Nc<11) 
! (RIND = Random Integration N Dimensions) 
!
!CALL RINDD(E,S,m,xc,Nt,indI,Blo,Bup,INFIN);
!
!        E = expectation/density as explained above size 1 x Nx        (out)
!        S = Covariance matrix of X=[Xt;Xd;Xc] size N x N (N=Nt+Nd+Nc) (in)
!        m = the expectation of X=[Xt;Xd;Xc]   size N x 1              (in)
!       xc = values to condition on            size Nc x Nx            (in)
!     indI = vector of indices to the different barriers in the        (in) 
!            indicator function,  length NI, where   NI = Nb+1 
!            (NB! restriction  indI(1)=0, indI(NI)=Nt+Nd )
!  Blo,Bup = Lower and upper barrier coefficients used to compute the   (in)
!            integration limits A and B, respectively. 
!            size  Mb x Nb. If  Mb<Nc+1 then
!            Blo(Mb+1:Nc+1,:) is assumed to be zero. 
!    INFIN = INTEGER, array of integration limits flags:  size 1 x Nb   (in)
!            if INFIN(I) < 0, Ith limits are (-infinity, infinity);
!            if INFIN(I) = 0, Ith limits are (-infinity, Hup(I)];
!            if INFIN(I) = 1, Ith limits are [Hlo(I), infinity);
!            if INFIN(I) = 2, Ith limits are [Hlo(I), Hup(I)].
! 
! The relation to the integration limits Hlo and Hup are as follows
!    IF INFIN(j)>=0,
!      IF INFIN(j)~=0,  A(i)=Blo(1,j)+Blo(2:Mb,j).'*xc(1:Mb-1,ix), 
!      IF INFIN(j)~=1,  B(i)=Bup(1,j)+Bup(2:Mb,j).'*xc(1:Mb-1,ix), 
!
!            where i=indI(j-1)+1:indI(j), j=1:NI-1, ix=1:Nx
!            Thus the integration limits may change with the conditional
!            variables.
!Example: 
! The indices, indI=[0 3 5 6], and coefficients Blo=[0 0 -1], Bup=[0 0 5],INFIN=[0 1 2] 
! means that   A = [-inf -inf -inf 0 0 -1]  B = [0 0 0 inf inf 5] 
!
! The GLOBALDATA and FUNCMOD modules are used to transport the inputs: 
!     SCIS = 1 Integrate all by SADAPT if NDIM<9 otherwise by KRBVRC (default)
!            2 Integrate all by SADAPT by Genz (1992) (Fast and reliable)
!            3 Integrate all by KRBVRC by Genz (1998) (Fast and reliable)
!            4 Integrate all by KROBOV by Genz (1992) (Fast and reliable)
!            5 Integrate all by RCRUDE by Genz (1992) (Reliable)
!
! (Recommended limitations Nx,Nt<101, Nd<7 and Nc<11)
! Also note that the size information have to be transferred to RINDD through the
! input arguments E,S,m,Nt,IndI,Blo,Bup and INFIN 
!
! if SCIS > 0 then you must initialize the random generator before you 
!  call rindd by the following lines:
!
!      call random_seed(SIZE=seed_size) 
!      allocate(seed(seed_size)) 
!      call random_seed(GET=seed(1:seed_size))  ! get current seed
!      seed(1)=seed1                            ! change seed
!      call random_seed(PUT=seed(1:seed_size)) 
!      deallocate(seed)
!
! For further description see the modules 
!
! References
! Podgorski et. al. (1999)
! "Exact distributions for apparent waves in irregular seas"
! Ocean Engineering                                                    (RINDD)
!
! R. Ambartzumian, A. Der Kiureghian, V. Ohanian and H.
! Sukiasian (1998)
! "Multinormal probabilities by sequential conditioned 
!  importance sampling: theory and application"                        (MVNFUN)
! Probabilistic Engineering Mechanics, Vol. 13, No 4. pp 299-308  
!
! Alan Genz (1992)
! 'Numerical Computation of Multivariate Normal Probabilites'          (MVNFUN)
! J. computational Graphical Statistics, Vol.1, pp 141--149
!
! Alan Genz and Koon-Shing Kwong (1999/2000?)
! 'Numerical Evaluation of Singular Multivariate Normal Distributions' (MVNFUN,COVSRT)
! Submitted to Computational Statistics and Data analysis
!
!

! Tested on:  DIGITAL UNIX Fortran90 compiler
!             PC pentium II with Lahey Fortran90 compiler
!             Solaris with SunSoft F90 compiler Version 1.0.1.0  (21229283) 
! History:
! revised pab 20.05.2003
!  - fixed a bug in covsrt
!  -replaced swapre and swapint with SWAPMOD
! revised pab 10.04.2003
! -fixed bugs in covsrt1 and covsrt
! revised pab 23.03.2003
! - fixed bug in rcswap
! - added erfcoremod
! revised pab 22.02.2003
! -removed some code and cleaned up some code
! TODO: there is a bug in rcswap
! revised pab 18.02.2003
! -commented out all print statements.
! -added funcmod1
! -added mvnlimits to handle differences in probabilities
!              in upper tails much more accurately, 
! revised pab 05.05.2000
!  - found a bug in funcmod2 and funcmod when Nd>0
! revised pab 11.04.2000
!  - found a bug in mvnfun: The values of xd(I) was not computed correctly 
! revised pab 04.04.2000
!  This is based on Alan Genz (1992) and 
!  Alan Genz and Koon-Shing Kwong (1999/2000?)
!  Program for calculating Singular Multivariate probabilites.
!  Currently extended/updated for calculating multivariate Singular
!  Expectations
!  Also updated from FORTRAN77 to FORTRAN90
!  - ADDED THL, EXLMS, EXINV
!*********************************************************************
      MODULE GLOBALCONST        ! global constants	
      IMPLICIT NONE
      DOUBLE PRECISION, PARAMETER :: gSQTWPI1= 0.39894228040143D0  !=1/sqrt(2*pi)
      DOUBLE PRECISION, PARAMETER :: gSQPI1  = 0.56418958354776D0  !=1/sqrt(pi)
      DOUBLE PRECISION, PARAMETER :: gSQPI   = 1.77245385090552D0  !=sqrt(pi)
      DOUBLE PRECISION, PARAMETER :: gSQTW   = 1.41421356237310D0  !=sqrt(2)
      DOUBLE PRECISION, PARAMETER :: gSQTW1  = 0.70710678118655D0  !=1/sqrt(2)
      DOUBLE PRECISION, PARAMETER :: gPI1    = 0.31830988618379D0  !=1/pi
      DOUBLE PRECISION, PARAMETER :: gPI     = 3.14159265358979D0  !=pi
      DOUBLE PRECISION, PARAMETER :: gTWPI   = 6.28318530717958D0  !=2*pi
      DOUBLE PRECISION, PARAMETER :: gSQTWPI = 2.50662827463100D0  !=sqrt(2*pi)
      DOUBLE PRECISION, PARAMETER :: gONE    = 1.D0
      DOUBLE PRECISION, PARAMETER :: gTWO    = 2.D0
      DOUBLE PRECISION, PARAMETER :: gHALF   = 0.5D0
      DOUBLE PRECISION, PARAMETER :: gZERO   = 0.D0
      DOUBLE PRECISION, PARAMETER :: gINFINITY = 37.D0 ! SQRT(-gTWO*LOG(1.D+12*TINY(gONE)))
!     Set gINFINITY (infinity).
!     Such that EXP(-2.x^2) > 10^(12) times TINY
!     SAVE gINFINITY
      END MODULE GLOBALCONST
      MODULE ERFCOREMOD
      IMPLICIT NONE

      INTERFACE CALERF
      MODULE PROCEDURE CALERF
      END INTERFACE 

      INTERFACE DERF
      MODULE PROCEDURE DERF
      END INTERFACE 

      INTERFACE DERFC
      MODULE PROCEDURE DERFC
      END INTERFACE 

      INTERFACE DERFCX
      MODULE PROCEDURE DERFCX
      END INTERFACE 
      CONTAINS
C--------------------------------------------------------------------
C
C DERF subprogram computes approximate values for erf(x).
C   (see comments heading CALERF).
C
C   Author/date: W. J. Cody, January 8, 1985
C
C--------------------------------------------------------------------
      FUNCTION DERF( X ) RESULT (VALUE)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN)  :: X
      DOUBLE PRECISION   :: VALUE
      INTEGER, PARAMETER :: JINT = 0
      CALL CALERF(X,VALUE,JINT)
      RETURN
      END FUNCTION DERF
C--------------------------------------------------------------------
C
C DERFC subprogram computes approximate values for erfc(x).
C   (see comments heading CALERF).
C
C   Author/date: W. J. Cody, January 8, 1985
C
C--------------------------------------------------------------------
      FUNCTION DERFC( X ) RESULT (VALUE)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN)  :: X
      DOUBLE PRECISION :: VALUE
      INTEGER, PARAMETER :: JINT = 1
      CALL CALERF(X,VALUE,JINT)
      RETURN
      END FUNCTION DERFC
C------------------------------------------------------------------
C
C DERFCX subprogram computes approximate values for exp(x*x) * erfc(x).
C   (see comments heading CALERF).
C
C   Author/date: W. J. Cody, March 30, 1987
C
C------------------------------------------------------------------
      FUNCTION DERFCX( X ) RESULT (VALUE)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN)  :: X
      DOUBLE PRECISION   :: VALUE
      INTEGER, PARAMETER :: JINT = 2
      CALL CALERF(X,VALUE,JINT)
      RETURN
      END FUNCTION DERFCX

      SUBROUTINE CALERF(ARG,RESULT,JINT)
      IMPLICIT NONE
C------------------------------------------------------------------
C
C CALERF packet evaluates  erf(x),  erfc(x),  and  exp(x*x)*erfc(x)
C   for a real argument  x.  It contains three FUNCTION type
C   subprograms: ERF, ERFC, and ERFCX (or DERF, DERFC, and DERFCX),
C   and one SUBROUTINE type subprogram, CALERF.  The calling
C   statements for the primary entries are:
C
C                   Y=ERF(X)     (or   Y=DERF(X)),
C
C                   Y=ERFC(X)    (or   Y=DERFC(X)),
C   and
C                   Y=ERFCX(X)   (or   Y=DERFCX(X)).
C
C   The routine  CALERF  is intended for internal packet use only,
C   all computations within the packet being concentrated in this
C   routine.  The function subprograms invoke  CALERF  with the
C   statement
C
C          CALL CALERF(ARG,RESULT,JINT)
C
C   where the parameter usage is as follows
C
C      Function                     Parameters for CALERF
C       call              ARG                  Result          JINT
C
C     ERF(ARG)      ANY REAL ARGUMENT         ERF(ARG)          0
C     ERFC(ARG)     ABS(ARG) .LT. XBIG        ERFC(ARG)         1
C     ERFCX(ARG)    XNEG .LT. ARG .LT. XMAX   ERFCX(ARG)        2
C
C   The main computation evaluates near-minimax approximations
C   from "Rational Chebyshev approximations for the error function"
C   by W. J. Cody, Math. Comp., 1969, PP. 631-638.  This
C   transportable program uses rational functions that theoretically
C   approximate  erf(x)  and  erfc(x)  to at least 18 significant
C   decimal digits.  The accuracy achieved depends on the arithmetic
C   system, the compiler, the intrinsic functions, and proper
C   selection of the machine-dependent constants.
C
C*******************************************************************
C*******************************************************************
C
C Explanation of machine-dependent constants
C
C   XMIN   = the smallest positive floating-point number.
C   XINF   = the largest positive finite floating-point number.
C   XNEG   = the largest negative argument acceptable to ERFCX;
C            the negative of the solution to the equation
C            2*exp(x*x) = XINF.
C   XSMALL = argument below which erf(x) may be represented by
C            2*x/sqrt(pi)  and above which  x*x  will not underflow.
C            A conservative value is the largest machine number X
C            such that   1.0 + X = 1.0   to machine precision.
C   XBIG   = largest argument acceptable to ERFC;  solution to
C            the equation:  W(x) * (1-0.5/x**2) = XMIN,  where
C            W(x) = exp(-x*x)/[x*sqrt(pi)].
C   XHUGE  = argument above which  1.0 - 1/(2*x*x) = 1.0  to
C            machine precision.  A conservative value is
C            1/[2*sqrt(XSMALL)]
C   XMAX   = largest acceptable argument to ERFCX; the minimum
C            of XINF and 1/[sqrt(pi)*XMIN].
C
C   Approximate values for some important machines are:
C
C                          XMIN       XINF        XNEG     XSMALL
C
C    C 7600      (S.P.)  3.13E-294   1.26E+322   -27.220  7.11E-15
C  CRAY-1        (S.P.)  4.58E-2467  5.45E+2465  -75.345  7.11E-15
C  IEEE (IBM/XT,
C    SUN, etc.)  (S.P.)  1.18E-38    3.40E+38     -9.382  5.96E-8
C  IEEE (IBM/XT,
C    SUN, etc.)  (D.P.)  2.23D-308   1.79D+308   -26.628  1.11D-16
C  IBM 195       (D.P.)  5.40D-79    7.23E+75    -13.190  1.39D-17
C  UNIVAC 1108   (D.P.)  2.78D-309   8.98D+307   -26.615  1.73D-18
C  VAX D-Format  (D.P.)  2.94D-39    1.70D+38     -9.345  1.39D-17
C  VAX G-Format  (D.P.)  5.56D-309   8.98D+307   -26.615  1.11D-16
C
C
C                          XBIG       XHUGE       XMAX
C
C    C 7600      (S.P.)  25.922      8.39E+6     1.80X+293
C  CRAY-1        (S.P.)  75.326      8.39E+6     5.45E+2465
C  IEEE (IBM/XT,
C    SUN, etc.)  (S.P.)   9.194      2.90E+3     4.79E+37
C  IEEE (IBM/XT,
C    SUN, etc.)  (D.P.)  26.543      6.71D+7     2.53D+307
C  IBM 195       (D.P.)  13.306      1.90D+8     7.23E+75
C  UNIVAC 1108   (D.P.)  26.582      5.37D+8     8.98D+307
C  VAX D-Format  (D.P.)   9.269      1.90D+8     1.70D+38
C  VAX G-Format  (D.P.)  26.569      6.71D+7     8.98D+307
C
C*******************************************************************
C*******************************************************************
C
C Error returns
C
C  The program returns  ERFC = 0      for  ARG .GE. XBIG;
C
C                       ERFCX = XINF  for  ARG .LT. XNEG;
C      and
C                       ERFCX = 0     for  ARG .GE. XMAX.
C
C
C Intrinsic functions required are:
C
C     ABS, AINT, EXP
C
C
C  Author: W. J. Cody
C          Mathematics and Computer Science Division
C          Argonne National Laboratory
C          Argonne, IL 60439
C
C  Latest modification: March 19, 1990
C  Updated to F90 by pab 23.03.2003
C
C------------------------------------------------------------------
      DOUBLE PRECISION, INTENT(IN) :: ARG
      INTEGER, INTENT(IN)          :: JINT
      DOUBLE PRECISION, INTENT(INOUT):: RESULT
! Local variables
      INTEGER :: I
      DOUBLE PRECISION :: DEL,X,XDEN,XNUM,Y,YSQ
C------------------------------------------------------------------
C  Mathematical constants
C------------------------------------------------------------------
      DOUBLE PRECISION, PARAMETER :: ZERO   = 0.0D0
      DOUBLE PRECISION, PARAMETER :: HALF   = 0.05D0
      DOUBLE PRECISION, PARAMETER :: ONE    = 1.0D0
      DOUBLE PRECISION, PARAMETER :: TWO    = 2.0D0
      DOUBLE PRECISION, PARAMETER :: FOUR   = 4.0D0
      DOUBLE PRECISION, PARAMETER :: SIXTEN = 16.0D0
      DOUBLE PRECISION, PARAMETER :: SQRPI  = 5.6418958354775628695D-1
      DOUBLE PRECISION, PARAMETER :: THRESH = 0.46875D0
C------------------------------------------------------------------
C  Machine-dependent constants
C------------------------------------------------------------------
      DOUBLE PRECISION, PARAMETER :: XNEG   = -26.628D0
      DOUBLE PRECISION, PARAMETER :: XSMALL = 1.11D-16
      DOUBLE PRECISION, PARAMETER :: XBIG   = 26.543D0
      DOUBLE PRECISION, PARAMETER :: XHUGE  = 6.71D7
      DOUBLE PRECISION, PARAMETER :: XMAX   = 2.53D307
      DOUBLE PRECISION, PARAMETER :: XINF   = 1.79D308  
!---------------------------------------------------------------
!     Coefficents to the rational polynomials
!--------------------------------------------------------------
      DOUBLE PRECISION, DIMENSION(5) :: A, Q
      DOUBLE PRECISION, DIMENSION(4) :: B
      DOUBLE PRECISION, DIMENSION(9) :: C
      DOUBLE PRECISION, DIMENSION(8) :: D
      DOUBLE PRECISION, DIMENSION(6) :: P
C------------------------------------------------------------------
C  Coefficients for approximation to  erf  in first interval
C------------------------------------------------------------------
      PARAMETER (A = (/ 3.16112374387056560D00,
     &     1.13864154151050156D02,3.77485237685302021D02,
     &     3.20937758913846947D03, 1.85777706184603153D-1/))
      PARAMETER ( B = (/2.36012909523441209D01,2.44024637934444173D02,
     &       1.28261652607737228D03,2.84423683343917062D03/))
C------------------------------------------------------------------
C  Coefficients for approximation to  erfc  in second interval
C------------------------------------------------------------------
      PARAMETER ( C=(/5.64188496988670089D-1,8.88314979438837594D0,
     1       6.61191906371416295D01,2.98635138197400131D02,
     2       8.81952221241769090D02,1.71204761263407058D03,
     3       2.05107837782607147D03,1.23033935479799725D03,
     4       2.15311535474403846D-8/))
      PARAMETER ( D =(/1.57449261107098347D01,1.17693950891312499D02, 
     1       5.37181101862009858D02,1.62138957456669019D03,
     2       3.29079923573345963D03,4.36261909014324716D03,
     3       3.43936767414372164D03,1.23033935480374942D03/))
C------------------------------------------------------------------
C  Coefficients for approximation to  erfc  in third interval
C------------------------------------------------------------------
      PARAMETER ( P =(/3.05326634961232344D-1,3.60344899949804439D-1,
     1       1.25781726111229246D-1,1.60837851487422766D-2,
     2       6.58749161529837803D-4,1.63153871373020978D-2/))
      PARAMETER (Q =(/2.56852019228982242D00,1.87295284992346047D00,
     1       5.27905102951428412D-1,6.05183413124413191D-2,
     2       2.33520497626869185D-3/))
C------------------------------------------------------------------
      X = ARG
      Y = ABS(X)
      IF (Y .LE. THRESH) THEN
C------------------------------------------------------------------
C  Evaluate  erf  for  |X| <= 0.46875
C------------------------------------------------------------------
         !YSQ = ZERO
         IF (Y .GT. XSMALL) THEN
            YSQ = Y * Y
            XNUM = A(5)*YSQ
            XDEN = YSQ
            DO  I = 1, 3
               XNUM = (XNUM + A(I)) * YSQ
               XDEN = (XDEN + B(I)) * YSQ
            END DO
            RESULT = X * (XNUM + A(4)) / (XDEN + B(4))
         ELSE
            RESULT = X *  A(4) / B(4)
         ENDIF
         IF (JINT .NE. 0) RESULT = ONE - RESULT
         IF (JINT .EQ. 2) RESULT = EXP(YSQ) * RESULT
         GO TO 800
C------------------------------------------------------------------
C     Evaluate  erfc  for 0.46875 <= |X| <= 4.0
C------------------------------------------------------------------
      ELSE IF (Y .LE. FOUR) THEN
         XNUM = C(9)*Y
         XDEN = Y
         DO I = 1, 7
            XNUM = (XNUM + C(I)) * Y
            XDEN = (XDEN + D(I)) * Y
         END DO
         RESULT = (XNUM + C(8)) / (XDEN + D(8))
         IF (JINT .NE. 2) THEN
            YSQ = AINT(Y*SIXTEN)/SIXTEN
            DEL = (Y-YSQ)*(Y+YSQ)
            RESULT = EXP(-YSQ*YSQ) * EXP(-DEL) * RESULT
         END IF
C------------------------------------------------------------------
C     Evaluate  erfc  for |X| > 4.0
C------------------------------------------------------------------
      ELSE
         RESULT = ZERO
         IF (Y .GE. XBIG) THEN
            IF ((JINT .NE. 2) .OR. (Y .GE. XMAX)) GO TO 300
            IF (Y .GE. XHUGE) THEN
               RESULT = SQRPI / Y
               GO TO 300
            END IF
         END IF
         YSQ = ONE / (Y * Y)
         XNUM = P(6)*YSQ
         XDEN = YSQ
         DO I = 1, 4
            XNUM = (XNUM + P(I)) * YSQ
            XDEN = (XDEN + Q(I)) * YSQ
         ENDDO
         RESULT = YSQ *(XNUM + P(5)) / (XDEN + Q(5))
         RESULT = (SQRPI -  RESULT) / Y
         IF (JINT .NE. 2) THEN
            YSQ = AINT(Y*SIXTEN)/SIXTEN
            DEL = (Y-YSQ)*(Y+YSQ)
            RESULT = EXP(-YSQ*YSQ) * EXP(-DEL) * RESULT
         END IF
      END IF
C------------------------------------------------------------------
C     Fix up for negative argument, erf, etc.
C------------------------------------------------------------------
 300  IF (JINT .EQ. 0) THEN
         RESULT = (HALF - RESULT) + HALF
         IF (X .LT. ZERO) RESULT = -RESULT
      ELSE IF (JINT .EQ. 1) THEN
         IF (X .LT. ZERO) RESULT = TWO - RESULT
      ELSE
         IF (X .LT. ZERO) THEN
            IF (X .LT. XNEG) THEN
               RESULT = XINF
            ELSE
               YSQ = AINT(X*SIXTEN)/SIXTEN
               DEL = (X-YSQ)*(X+YSQ)
               Y = EXP(YSQ*YSQ) * EXP(DEL)
               RESULT = (Y+Y) - RESULT
            END IF
         END IF
      END IF
 800  RETURN
      END SUBROUTINE CALERF
      END MODULE ERFCOREMOD
      
      MODULE GLOBALDATA
      IMPLICIT NONE    
      PUBLIC :: INITDATA, SETDATA
      PRIVATE :: FIINV, FI

      INTERFACE FIINV
      MODULE PROCEDURE FIINV
      END INTERFACE
      
      INTERFACE FI
      MODULE PROCEDURE FI
      END INTERFACE 

      INTERFACE INITDATA
      MODULE PROCEDURE INITDATA
      END INTERFACE
      
      INTERFACE SETDATA
      MODULE PROCEDURE SETDATA
      END INTERFACE 
!     Constants determining accuracy of integration
!-----------------------------------------------
!     if the conditional variance are less than EPS2, the variable is 
!     considered deterministic : 
      DOUBLE PRECISION,SAVE :: EPS2   = 1.D-10 
      DOUBLE PRECISION,SAVE :: EPS    = 1.D-5  ! = SQRT(EPS2)
      DOUBLE PRECISION,SAVE :: XCEPS2 = 1.D-16 ! if Var(Xc) is less return NaN
      DOUBLE PRECISION,SAVE :: ABSEPS = 1.D-3  ! requested absolute error 
      DOUBLE PRECISION,SAVE :: RELEPS = 1.D-3  ! requested Relative error, i.e. if 
                                ! 3.0*STD(XIND)/XIND is less we accept the estimate
      DOUBLE PRECISION,SAVE :: fxcEpss = 1.D-20 ! if less do not compute E(...|Xc)
      DOUBLE PRECISION,SAVE :: xCutOff = 5.D0  ! upper/lower truncation limit of the 
                                       ! normal CDF 
      DOUBLE PRECISION,SAVE :: FxCutOff  = 0.99999942669686D0 
      DOUBLE PRECISION,SAVE :: CFxCutOff = 5.733031438470704D-7  ! 1-FxCutOff, 
!      DOUBLE PRECISION,SAVE :: XSMALL    = 4.2D-16 !cut off parameter to FI
      DOUBLE PRECISION,SAVE :: XMAX    = 8.D0 !5.d0 !8.29287554336168D0 !cut off parameter to FI
    
!     Parameters controlling the performance and integration method
      INTEGER, SAVE :: SCIS=1   !=1 Integrate by SADAPT if Ndim<9 KRBVRC otherwise 
                                !=2 Integrate all by SADAPT  (Fast and reliable)
                                !=3 Integrate all by KRBVRC  (Fastest and reliable)
                                !=4 Integrate all by KROBOV  (Fast and reliable)
                                !=5 Integrate all by RCRUDE  (Reliable)
      INTEGER, SAVE :: NIT     = 1000  ! maximum number Xt variables to integrate
      INTEGER, SAVE :: NSIMmax = 10000 ! maximum number of simulations
      INTEGER, SAVE :: NSIMmin = 0		! minimum number of
                                ! simulations
      CONTAINS
      FUNCTION FIINV(P) RESULT (VAL)
!      use GLOBALDATA, only: XMAX,CFxCutOff
      IMPLICIT NONE
*
*	ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3
*
*	Produces the normal deviate Z corresponding to a given lower
*	tail area of P.
*       Absolute error less than 1e-13
*       Relative error less than 1e-15 for abs(VAL)>0.1
*
*	The hash sums below are the sums of the mantissas of the
*	coefficients.   They are included for use in checking
*	transcription.
*
      DOUBLE PRECISION, INTENT(in) :: P
      DOUBLE PRECISION :: VAL
!local variables
      DOUBLE PRECISION SPLIT1, SPLIT2, CONST1, CONST2, ONE, ZERO, HALF, 
     &     A0, A1, A2, A3, A4, A5, A6, A7, B1, B2, B3, B4, B5, B6, B7, 
     &     C0, C1, C2, C3, C4, C5, C6, C7, D1, D2, D3, D4, D5, D6, D7, 
     &     E0, E1, E2, E3, E4, E5, E6, E7, F1, F2, F3, F4, F5, F6, F7, 
     &     Q, R
      PARAMETER ( SPLIT1 = 0.425D0, SPLIT2 = 5.D0,
     &            CONST1 = 0.180625D0, CONST2 = 1.6D0,
     &				ONE = 1.D0, ZERO = 0.D0, HALF = 0.5D0 )
*     
*     Coefficients for P close to 0.5
*     
      PARAMETER (
     *     A0 = 3.38713 28727 96366 6080D0,
     *     A1 = 1.33141 66789 17843 7745D+2,
     *     A2 = 1.97159 09503 06551 4427D+3,
     *     A3 = 1.37316 93765 50946 1125D+4,
     *     A4 = 4.59219 53931 54987 1457D+4,
     *     A5 = 6.72657 70927 00870 0853D+4,
     *     A6 = 3.34305 75583 58812 8105D+4,
     *     A7 = 2.50908 09287 30122 6727D+3,
     *     B1 = 4.23133 30701 60091 1252D+1,
     *     B2 = 6.87187 00749 20579 0830D+2,
     *     B3 = 5.39419 60214 24751 1077D+3,
     *     B4 = 2.12137 94301 58659 5867D+4,
     *     B5 = 3.93078 95800 09271 0610D+4,
     *     B6 = 2.87290 85735 72194 2674D+4,
     *     B7 = 5.22649 52788 52854 5610D+3 )
*     HASH SUM AB    55.88319 28806 14901 4439
*     
*     Coefficients for P not close to 0, 0.5 or 1.
*     
      PARAMETER (
     *     C0 = 1.42343 71107 49683 57734D0,
     *     C1 = 4.63033 78461 56545 29590D0,
     *     C2 = 5.76949 72214 60691 40550D0,
     *     C3 = 3.64784 83247 63204 60504D0,
     *     C4 = 1.27045 82524 52368 38258D0,
     *     C5 = 2.41780 72517 74506 11770D-1,
     *     C6 = 2.27238 44989 26918 45833D-2,
     *     C7 = 7.74545 01427 83414 07640D-4,
     *     D1 = 2.05319 16266 37758 82187D0,
     *     D2 = 1.67638 48301 83803 84940D0,
     *     D3 = 6.89767 33498 51000 04550D-1,
     *     D4 = 1.48103 97642 74800 74590D-1,
     *     D5 = 1.51986 66563 61645 71966D-2,
     *     D6 = 5.47593 80849 95344 94600D-4,
     *     D7 = 1.05075 00716 44416 84324D-9 )
*     HASH SUM CD    49.33206 50330 16102 89036
*
*	Coefficients for P near 0 or 1.
*
      PARAMETER (
     *     E0 = 6.65790 46435 01103 77720D0,
     *     E1 = 5.46378 49111 64114 36990D0,
     *     E2 = 1.78482 65399 17291 33580D0,
     *     E3 = 2.96560 57182 85048 91230D-1,
     *     E4 = 2.65321 89526 57612 30930D-2,
     *     E5 = 1.24266 09473 88078 43860D-3,
     *     E6 = 2.71155 55687 43487 57815D-5,
     *     E7 = 2.01033 43992 92288 13265D-7,
     *     F1 = 5.99832 20655 58879 37690D-1,
     *     F2 = 1.36929 88092 27358 05310D-1,
     *     F3 = 1.48753 61290 85061 48525D-2,
     *     F4 = 7.86869 13114 56132 59100D-4,
     *     F5 = 1.84631 83175 10054 68180D-5,
     *     F6 = 1.42151 17583 16445 88870D-7,
     *     F7 = 2.04426 31033 89939 78564D-15 )
*     HASH SUM EF    47.52583 31754 92896 71629
*     
      Q = ( P - HALF)
      IF ( ABS(Q) .LE. SPLIT1 ) THEN ! Central range.
         R = CONST1 - Q*Q
         VAL = Q*( ( ( ((((A7*R + A6)*R + A5)*R + A4)*R + A3)
     *                  *R + A2 )*R + A1 )*R + A0 )
     *            /( ( ( ((((B7*R + B6)*R + B5)*R + B4)*R + B3)
     *                  *R + B2 )*R + B1 )*R + ONE)
      ELSE ! near the endpoints
         R = MIN( P, ONE - P )
         IF  (R .GT.ZERO) THEN ! ( 2.d0*R .GT. CFxCutOff) THEN ! R .GT.0.d0
            R = SQRT( -LOG(R) )
            IF ( R .LE. SPLIT2 ) THEN
               R = R - CONST2
               VAL = ( ( ( ((((C7*R + C6)*R + C5)*R + C4)*R + C3)
     *                      *R + C2 )*R + C1 )*R + C0 ) 
     *                /( ( ( ((((D7*R + D6)*R + D5)*R + D4)*R + D3)
     *                      *R + D2 )*R + D1 )*R + ONE )
            ELSE
               R = R - SPLIT2
               VAL = ( ( ( ((((E7*R + E6)*R + E5)*R + E4)*R + E3)
     *                      *R + E2 )*R + E1 )*R + E0 )
     *                /( ( ( ((((F7*R + F6)*R + F5)*R + F4)*R + F3)
     *                      *R + F2 )*R + F1 )*R + ONE )
            END IF
         ELSE
            VAL = 37.D0 !XMAX 9.d0
         END IF
         IF ( Q .LT. ZERO ) VAL = - VAL
      END IF
      RETURN
      END FUNCTION FIINV     
      FUNCTION FI( Z ) RESULT (VALUE)
      USE ERFCOREMOD
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in) :: Z
      DOUBLE PRECISION :: VALUE
! Local variables
      DOUBLE PRECISION, PARAMETER:: SQ2M1 = 0.70710678118655D0 !     1/SQRT(2)
      DOUBLE PRECISION, PARAMETER:: HALF = 0.5D0
      VALUE = DERFC(-Z*SQ2M1)*HALF
      RETURN
      END FUNCTION FI
      

      SUBROUTINE INITDATA (speed)                                   
!      USE GLOBALDATA
!      USE FIMOD
      IMPLICIT NONE
      INTEGER , INTENT(in) :: speed
      INTEGER :: TMP
      TMP= max(11-speed,1)
      IF (.TRUE.) THEN
         NsimMax=10000
         SELECT case (speed)
         CASE (11:)
            ABSEPS = 1d-1 
         CASE (10)
            ABSEPS = 1d-2 
         CASE (7:9) 
            ABSEPS = 1d-2 
         CASE (4:6)
            NsimMax=20000
            ABSEPS = 1d-3
         CASE (:3)
            NsimMax=30000
            ABSEPS = 1d-4 
         END SELECT 
       
         RELEPS = MIN(ABSEPS ,1.d-2)
         TMP    = MOD(TMP+1,3)+1
         EPS2   = ABSEPS*(1D-1**TMP)
         !EPS2=MIN(EPS2,1D-6)
         !EPS2=1D-10
         xCutOff = ABS(FIINV(ABSEPS*1D-1*(1D-1**TMP)))
         xCutOff = MAX(xCutOff+0.5d0,4.d0)
      ELSE         
         NsimMax = 10000*TMP
         NsimMin = 0
         ABSEPS  = (1.0D-1)**TMP
         
         
         IF (.FALSE.) THEN
            EPS2    = 1.0D-10
            xCutOff = XMAX 
         ELSE
            
            xCutOff = ABS(FIINV(ABSEPS))
            xCutOff = MAX(xCutOff+0.5d0,4.d0)
                                !xCutOff= MIN(xCutOff,5.d0)
            EPS2   = ABSEPS*1.d-2
            ABSEPS = MIN(ABSEPS*(10**3),0.1D0)
         ENDIF
      ENDIF

      RELEPS = MIN(ABSEPS*1.d-1 ,1.d-2)
      EPS    = SQRT(EPS2)

      CFxCutOff = FI(-xCutOff)*2.0D0 !1.d0-FxCutOff  
      FxCutOff  = 1.0D0 - CFxCutOff  !FI(xCutOff) - FI( -xCutOff) 
     
      XMAX      = xCutOff
      RETURN
      PRINT *, 'SP2THPDFALAN.exe CALLED'
      print *,'INITDATA: Requested parameters :'
      SELECT CASE (SCIS)
      CASE (:1)
         PRINT *,'SCIS = 1 SADAPT if NDIM<9 otherwise by KRBVRC'
      CASE (2) 
         PRINT *,'SCIS = 2 SADAPT'
      CASE (3)
         PRINT *,'SCIS = 3 KRBVRC'
      CASE (4)
         PRINT *,'SCIS = 4 KROBOV'
      CASE (5) 
         PRINT *,'SCIS = 5 RCRUDE'
      CASE (6) 
         PRINT *,'SCIS = 6 RLHCRUDE using MLHD center of cell' 
      CASE (7)
         PRINT *,'SCIS = 7 RLHCRUDE using LHD and center of cell'
      CASE (8)
         PRINT *,'SCIS = 8 RLHCRUDE using MLHD and random point'
      CASE (9:)
         PRINT *,'SCIS = 9 RLHCRUDE using LHD and random point'
      END SELECT
	PRINT *, 'ABSEPS = ', ABSEPS, ' RELEPS = ', RELEPS
	PRINT *, 'EPS2 = ', EPS2, ' xCutOff = ', xCutOff
	PRINT *, 'NsimMax = ', NsimMax
      RETURN                                                             
      END SUBROUTINE INITDATA
      SUBROUTINE SETDATA (dREPS,dAEPS,dEPS2,dXc) 
!      USE GLOBALDATA
!      USE FIMOD, ONLY : FI
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in) :: dREPS,dAEPS,dEPS2,dXc 
      RELEPS   = dREPS          ! Constants controlling 
      ABSEPS   = dAEPS          ! accuracy of integration 
      EPS2     = dEPS2          
      EPS      = SQRT(EPS2)

      xCutOff  = dXc
      CFxCutOff= FI(-xCutOff)*2.0D0 !1.d0-FxCutOff  
      FxCutOff = 1.0D0 - CFxCutOff  !FI(xCutOff) - FI( -xCutOff) 
                                  
      XMAX  = xCutOff
      IF (ABSEPS.LE.1e-4)      NsimMax=20000
      IF (ABSEPS.LE.1e-5)      NsimMax=40000
      IF (ABSEPS.LE.1e-6)      NsimMax=80000
      RETURN
      print *,'SETDATA: Requested parameters :'
      SELECT CASE (SCIS)
      CASE (:1)
         PRINT *,'SCIS = 1 SADAPT if NDIM<9 otherwise by KRBVRC'
      CASE (2) 
         PRINT *,'SCIS = 2 SADAPT'
      CASE (3)  
         PRINT *,'SCIS = 3 KRBVRC'
      CASE (4) 
         PRINT *,'SCIS = 4 KROBOV'
      CASE (5)  
         PRINT *,'SCIS = 5 RCRUDE'
      CASE (6) 
         PRINT *,'SCIS = 6 RLHCRUDE using MLHD center of cell' 
      CASE (7) 
         PRINT *,'SCIS = 7 RLHCRUDE using LHD and center of cell'
      CASE (8) 
         PRINT *,'SCIS = 8 RLHCRUDE using MLHD and random point'
      CASE (9:) 
         PRINT *,'SCIS = 9 RLHCRUDE using LHD and random point'
      END SELECT
      print *,'ABSEPS,RELEPS,EPS2,xCutOff,NsimMax',
     &     ABSEPS,RELEPS,EPS2,xCutOff,  NsimMax
      RETURN            
      END SUBROUTINE SETDATA
      END MODULE GLOBALDATA
      
!
! FIMOD contains functions for calculating 1D and 2D Normal probabilites
!       and  1D expectations
      MODULE FIMOD
      IMPLICIT NONE
      PRIVATE
      PUBLIC :: NORMPRB, FI, FIINV, MVNLIMITS, MVNLMS, BVU,BVNMVN
      PUBLIC :: NORM2DPRB, THL, GAUSINT, GAUSINT2, EXLMS, EXINV

      INTERFACE NORMPRB
      MODULE PROCEDURE NORMPRB
      END INTERFACE 

      INTERFACE FI
      MODULE PROCEDURE FI
      END INTERFACE 
      
      INTERFACE FI2
      MODULE PROCEDURE FI2
      END INTERFACE 

      INTERFACE FIINV
      MODULE PROCEDURE FIINV
      END INTERFACE
      
      INTERFACE  MVNLIMITS
      MODULE PROCEDURE  MVNLIMITS
      END INTERFACE

      INTERFACE  MVNLMS
      MODULE PROCEDURE  MVNLMS
      END INTERFACE

      INTERFACE BVU
      MODULE PROCEDURE BVU
      END INTERFACE

      INTERFACE BVNMVN
      MODULE PROCEDURE BVNMVN
      END INTERFACE
      
      INTERFACE NORM2DPRB
      MODULE PROCEDURE NORM2DPRB
      END INTERFACE
      
      INTERFACE THL
      MODULE PROCEDURE THL
      END INTERFACE

      INTERFACE GAUSINT
      MODULE PROCEDURE GAUSINT
      END INTERFACE

      INTERFACE GAUSINT2
      MODULE PROCEDURE GAUSINT2
      END INTERFACE
      
      INTERFACE EXLMS
      MODULE PROCEDURE EXLMS
      END INTERFACE

      INTERFACE EXINV
      MODULE PROCEDURE EXINV
      END INTERFACE

      CONTAINS
      FUNCTION FIINV(P) RESULT (VAL)
      use GLOBALDATA, only: XMAX,CFxCutOff
      IMPLICIT NONE
*
*	ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3
*
*	Produces the normal deviate Z corresponding to a given lower
*	tail area of P.
*       Absolute error less than 1e-13
*       Relative error less than 1e-15 for abs(VAL)>0.1
*
*	The hash sums below are the sums of the mantissas of the
*	coefficients.   They are included for use in checking
*	transcription.
*
      DOUBLE PRECISION, INTENT(in) :: P
      DOUBLE PRECISION :: VAL
!local variables
      DOUBLE PRECISION SPLIT1, SPLIT2, CONST1, CONST2, ONE, ZERO, HALF, 
     &     A0, A1, A2, A3, A4, A5, A6, A7, B1, B2, B3, B4, B5, B6, B7, 
     &     C0, C1, C2, C3, C4, C5, C6, C7, D1, D2, D3, D4, D5, D6, D7, 
     &     E0, E1, E2, E3, E4, E5, E6, E7, F1, F2, F3, F4, F5, F6, F7, 
     &     Q, R
      PARAMETER ( SPLIT1 = 0.425D0, SPLIT2 = 5.D0,
     &     CONST1 = 0.180625D0, CONST2 = 1.6D0,
     &     ONE = 1.D0, ZERO = 0.D0, HALF = 0.5D0 )
*     
*     Coefficients for P close to 0.5
*     
      PARAMETER (
     *     A0 = 3.38713 28727 96366 6080D0,
     *     A1 = 1.33141 66789 17843 7745D+2,
     *     A2 = 1.97159 09503 06551 4427D+3,
     *     A3 = 1.37316 93765 50946 1125D+4,
     *     A4 = 4.59219 53931 54987 1457D+4,
     *     A5 = 6.72657 70927 00870 0853D+4,
     *     A6 = 3.34305 75583 58812 8105D+4,
     *     A7 = 2.50908 09287 30122 6727D+3,
     *     B1 = 4.23133 30701 60091 1252D+1,
     *     B2 = 6.87187 00749 20579 0830D+2,
     *     B3 = 5.39419 60214 24751 1077D+3,
     *     B4 = 2.12137 94301 58659 5867D+4,
     *     B5 = 3.93078 95800 09271 0610D+4,
     *     B6 = 2.87290 85735 72194 2674D+4,
     *     B7 = 5.22649 52788 52854 5610D+3 )
*     HASH SUM AB    55.88319 28806 14901 4439
*     
*     Coefficients for P not close to 0, 0.5 or 1.
*     
      PARAMETER (
     *     C0 = 1.42343 71107 49683 57734D0,
     *     C1 = 4.63033 78461 56545 29590D0,
     *     C2 = 5.76949 72214 60691 40550D0,
     *     C3 = 3.64784 83247 63204 60504D0,
     *     C4 = 1.27045 82524 52368 38258D0,
     *     C5 = 2.41780 72517 74506 11770D-1,
     *     C6 = 2.27238 44989 26918 45833D-2,
     *     C7 = 7.74545 01427 83414 07640D-4,
     *     D1 = 2.05319 16266 37758 82187D0,
     *     D2 = 1.67638 48301 83803 84940D0,
     *     D3 = 6.89767 33498 51000 04550D-1,
     *     D4 = 1.48103 97642 74800 74590D-1,
     *     D5 = 1.51986 66563 61645 71966D-2,
     *     D6 = 5.47593 80849 95344 94600D-4,
     *     D7 = 1.05075 00716 44416 84324D-9 )
*     HASH SUM CD    49.33206 50330 16102 89036
*
*	Coefficients for P near 0 or 1.
*
      PARAMETER (
     *     E0 = 6.65790 46435 01103 77720D0,
     *     E1 = 5.46378 49111 64114 36990D0,
     *     E2 = 1.78482 65399 17291 33580D0,
     *     E3 = 2.96560 57182 85048 91230D-1,
     *     E4 = 2.65321 89526 57612 30930D-2,
     *     E5 = 1.24266 09473 88078 43860D-3,
     *     E6 = 2.71155 55687 43487 57815D-5,
     *     E7 = 2.01033 43992 92288 13265D-7,
     *     F1 = 5.99832 20655 58879 37690D-1,
     *     F2 = 1.36929 88092 27358 05310D-1,
     *     F3 = 1.48753 61290 85061 48525D-2,
     *     F4 = 7.86869 13114 56132 59100D-4,
     *     F5 = 1.84631 83175 10054 68180D-5,
     *     F6 = 1.42151 17583 16445 88870D-7,
     *     F7 = 2.04426 31033 89939 78564D-15 )
*     HASH SUM EF    47.52583 31754 92896 71629
*     
      Q = ( P - HALF)
      IF ( ABS(Q) .LE. SPLIT1 ) THEN ! Central range.
         R = CONST1 - Q*Q
         VAL = Q*( ( ( ((((A7*R + A6)*R + A5)*R + A4)*R + A3)
     *                  *R + A2 )*R + A1 )*R + A0 )
     *            /( ( ( ((((B7*R + B6)*R + B5)*R + B4)*R + B3)
     *                  *R + B2 )*R + B1 )*R + ONE)
      ELSE ! near the endpoints
         R = MIN( P, ONE - P )
         IF  (R .GT.ZERO) THEN ! ( 2.d0*R .GT. CFxCutOff) THEN ! R .GT.0.d0
            R = SQRT( -LOG(R) )
            IF ( R .LE. SPLIT2 ) THEN
               R = R - CONST2
               VAL = ( ( ( ((((C7*R + C6)*R + C5)*R + C4)*R + C3)
     *                      *R + C2 )*R + C1 )*R + C0 ) 
     *                /( ( ( ((((D7*R + D6)*R + D5)*R + D4)*R + D3)
     *                      *R + D2 )*R + D1 )*R + ONE )
            ELSE
               R = R - SPLIT2
               VAL = ( ( ( ((((E7*R + E6)*R + E5)*R + E4)*R + E3)
     *                      *R + E2 )*R + E1 )*R + E0 )
     *                /( ( ( ((((F7*R + F6)*R + F5)*R + F4)*R + F3)
     *                      *R + F2 )*R + F1 )*R + ONE )
            END IF
         ELSE
            VAL = 37.0d0 !9.D0 !XMAX 9.d0
         END IF
         IF ( Q .LT. ZERO ) VAL = - VAL
      END IF
      RETURN
      END FUNCTION FIINV     
                                ! *********************************     
      SUBROUTINE NORMPRB(Z, P, Q)
      USE ERFCOREMOD
      USE GLOBALDATA, ONLY : XMAX      
! Normal distribution probabilities accurate to 18 digits between 
! -XMAX and XMAX
!
! Z    = no. of standard deviations from the mean.
! P, Q = probabilities to the left & right of Z.   P + Q = 1.
!
! by pab 23.03.2003
!
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN)            :: Z
      DOUBLE PRECISION, INTENT(OUT)           :: P
      DOUBLE PRECISION, INTENT(OUT), OPTIONAL ::  Q
!Local variables
      DOUBLE PRECISION            :: PP, QQ, ZABS
      DOUBLE PRECISION, PARAMETER :: ZERO  = 0.0D0
      DOUBLE PRECISION, PARAMETER :: ONE   = 1.0D0
      DOUBLE PRECISION, PARAMETER :: SQ2M1 = 0.70710678118655D0 !     1/SQRT(2)
      DOUBLE PRECISION, PARAMETER :: HALF  = 0.5D0
      ZABS = ABS(Z)     
!
!     |Z| > 37  (or XMAX)
!
      IF ( ZABS .GT. XMAX ) THEN
!         IF (PRESENT(PDF)) PDF = ZERO
         IF (Z > ZERO) THEN
            !IF (PRESENT(P)) 
            P = ONE
            IF (PRESENT(Q)) Q = ZERO
         ELSE
            !IF (PRESENT(P)) 
            P = ZERO
            IF (PRESENT(Q)) Q = ONE
         END IF
      ELSE
!     
!     |Z| <= 37
!     
         PP = DERFC(ZABS*SQ2M1)*HALF

         IF (Z .LT. ZERO) THEN
!     IF (PRESENT(P)) 
            P = PP
            IF (PRESENT(Q)) Q = ONE - PP
         ELSE
                                !IF (PRESENT(P)) 
            P = ONE - PP
            IF (PRESENT(Q)) Q = PP
         END IF
      END IF

      RETURN
      END SUBROUTINE NORMPRB
      FUNCTION FI( Z ) RESULT (VALUE)
      USE ERFCOREMOD
      USE GLOBALDATA, ONLY : XMAX
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in) :: Z
      DOUBLE PRECISION :: VALUE
! Local variables
      DOUBLE PRECISION :: ZABS
      DOUBLE PRECISION, PARAMETER:: SQ2M1 = 0.70710678118655D0 !     1/SQRT(2)
      DOUBLE PRECISION, PARAMETER:: HALF = 0.5D0
      ZABS = ABS(Z)
*     
*     |Z| > 37  (or XMAX)
*     
      IF ( ZABS .GT. XMAX ) THEN
         IF (Z .LT. 0.0D0) THEN
            VALUE = 0.0D0
         ELSE
            VALUE = 1.0D0
         ENDIF
      ELSE
         VALUE = DERFC(-Z*SQ2M1)*HALF
      ENDIF
      RETURN
      END FUNCTION FI
      
      FUNCTION FI2( Z ) RESULT (VALUE)
      USE GLOBALDATA, ONLY : XMAX
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in) :: Z
      DOUBLE PRECISION :: VALUE
*     
*     Normal distribution probabilities accurate to 1.e-15.
*     relative error less than 1e-8;
*     Z = no. of standard deviations from the mean.
*     
*     Based upon algorithm 5666 for the error function, from:
*     Hart, J.F. et al, 'Computer Approximations', Wiley 1968
*     
*     Programmer: Alan Miller
*     
*     Latest revision - 30 March 1986
*     
      DOUBLE PRECISION :: P0, P1, P2, P3, P4, P5, P6, 
     *     Q0, Q1, Q2, Q3, Q4, Q5, Q6, Q7,
     *     P, EXPNTL, CUTOFF, ROOTPI, ZABS
      PARAMETER(
     *     P0 = 220.20 68679 12376 1D0,
     *     P1 = 221.21 35961 69931 1D0, 
     *     P2 = 112.07 92914 97870 9D0,
     *     P3 = 33.912 86607 83830 0D0,
     *     P4 = 6.3739 62203 53165 0D0,
     *     P5 = 0.70038 30644 43688 1D0, 
     *     P6 = 0.035262 49659 98910 9D0 )
      PARAMETER(
     *     Q0 = 440.41 37358 24752 2D0,
     *     Q1 = 793.82 65125 19948 4D0, 
     *     Q2 = 637.33 36333 78831 1D0,
     *     Q3 = 296.56 42487 79673 7D0, 
     *     Q4 = 86.780 73220 29460 8D0,
     *     Q5 = 16.064 17757 92069 5D0, 
     *     Q6 = 1.7556 67163 18264 2D0,
     *     Q7 = 0.088388 34764 83184 4D0 )
      PARAMETER( ROOTPI = 2.5066 28274 63100 1D0 )
      PARAMETER( CUTOFF = 7.0710 67811 86547 5D0 )
*     
      ZABS = ABS(Z)
*     
*     |Z| > 37  (or XMAX)
*     
      IF ( ZABS .GT. XMAX ) THEN
         P = 0.d0
      ELSE
*     
*     |Z| <= 37
*     
         EXPNTL = EXP( -ZABS**2/2D0 )
*     
*     |Z| < CUTOFF = 10/SQRT(2)
*     
         IF ( ZABS .LT. CUTOFF ) THEN
            P = EXPNTL*( (((((P6*ZABS + P5)*ZABS + P4)*ZABS + P3)*ZABS
     *           + P2)*ZABS + P1)*ZABS + P0)/(((((((Q7*ZABS + Q6)*ZABS
     *           + Q5)*ZABS + Q4)*ZABS + Q3)*ZABS + Q2)*ZABS + Q1)*ZABS
     *           + Q0 )
*     
*     |Z| >= CUTOFF.
*     
         ELSE
            P = EXPNTL/( ZABS + 1.d0/( ZABS + 2.d0/( ZABS + 3.d0/( ZABS 
     *                        + 4.d0/( ZABS + 0.65D0 ) ) ) ) )/ROOTPI
         END IF
      END IF
      IF ( Z .GT. 0.d0 ) P = 1.d0 - P
      VALUE = P
      RETURN
      END FUNCTION FI2

      SUBROUTINE MVNLIMITS( A, B, INFIN, AP, PRB,AQ)
! RETURN probabilities for being between A and B
!  WHERE 
!  AP = FI(A), AQ = 1 - FI(A)
!  BP = FI(B), BQ = 1 - FI(B)
!  PRB = BP-AP IF BP+AP<1
!      = AQ-BQ OTHERWISE
!
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in) :: A, B
      DOUBLE PRECISION, INTENT(out) :: AP
      DOUBLE PRECISION, INTENT(out),OPTIONAL :: PRB,AQ
      INTEGER,INTENT(in) :: INFIN
!     LOCAL VARIABLES
      DOUBLE PRECISION :: BP,AQQ, BQQ
      DOUBLE PRECISION, PARAMETER :: ONE=1.D0, ZERO = 0.D0
      
      SELECT CASE (infin)
      CASE (:-1)
         AP  = ZERO
!     BP  = ONE
         IF (PRESENT(PRB)) PRB = ONE
         IF (PRESENT(AQ)) AQ = ONE
!     IF (PRESENT(BQ)) BQ = ZERO
      CASE (0)
         AP  = ZERO
         CALL NORMPRB(B,BP)     !,BQQ)
         IF (PRESENT(PRB)) PRB = BP
         IF (PRESENT(AQ)) AQ = ONE
!     IF (PRESENT(BQ)) BQ = BQQ
      CASE (1)
!     BP = ONE
         CALL NORMPRB(A,AP,AQQ)
         IF (PRESENT(PRB)) PRB = AQQ
         IF (PRESENT(AQ)) AQ = AQQ
!     IF (PRESENT(BQ)) BQ = ZERO
      CASE (2:)
         CALL NORMPRB(A,AP,AQQ)
         CALL NORMPRB(B,BP,BQQ)
         IF (PRESENT(PRB)) THEN
            IF (AP+BP .LT. ONE) THEN
               PRB = BP - AP
            ELSE
               PRB = AQQ - BQQ
            END IF
         ENDIF
         IF (PRESENT(AQ)) AQ = AQQ
!     IF (PRESENT(BQ)) BQ = BQQ
      END SELECT
      RETURN
      END SUBROUTINE MVNLIMITS 


      SUBROUTINE MVNLMS( A, B, INFIN, LOWER, UPPER )
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in) :: A, B
      DOUBLE PRECISION, INTENT(out) :: LOWER, UPPER
      INTEGER,INTENT(in) :: INFIN

      LOWER = 0.0D0
      UPPER = 1.0D0
      IF ( INFIN .LT. 0 ) RETURN
      IF ( INFIN .NE. 0 ) LOWER = FI(A)
      IF ( INFIN .NE. 1 ) UPPER = FI(B)
      RETURN
      END SUBROUTINE MVNLMS 

      
      FUNCTION BVNMVN( LOWER, UPPER, INFIN, CORREL ) RESULT (VAL)
      IMPLICIT NONE
*
*     A function for computing bivariate normal probabilities.
*
*  Parameters
*
*     LOWER  REAL, array of lower integration limits.
*     UPPER  REAL, array of upper integration limits.
*     INFIN  INTEGER, array of integration limits flags:
*            if INFIN(I) = 0, Ith limits are (-infinity, UPPER(I)];
*            if INFIN(I) = 1, Ith limits are [LOWER(I), infinity);
*            if INFIN(I) = 2, Ith limits are [LOWER(I), UPPER(I)].
*     CORREL REAL, correlation coefficient.
*
      
      DOUBLE PRECISION, DIMENSION(:), INTENT (IN) :: LOWER, UPPER
      DOUBLE PRECISION,               INTENT (IN) :: CORREL 
      INTEGER,          DIMENSION(:), INTENT (IN) :: INFIN
      DOUBLE PRECISION :: VAL
      IF ( INFIN(1) .EQ. 2  .AND. INFIN(2) .EQ. 2 ) THEN
         VAL =  BVU ( LOWER(1), LOWER(2), CORREL )
     +           - BVU ( UPPER(1), LOWER(2), CORREL )
     +           - BVU ( LOWER(1), UPPER(2), CORREL )
     +           + BVU ( UPPER(1), UPPER(2), CORREL )
      ELSE IF ( INFIN(1) .EQ. 2  .AND. INFIN(2) .EQ. 1 ) THEN
         VAL =  BVU ( LOWER(1), LOWER(2), CORREL )
     +           - BVU ( UPPER(1), LOWER(2), CORREL )
      ELSE IF ( INFIN(1) .EQ. 1  .AND. INFIN(2) .EQ. 2 ) THEN
         VAL =  BVU ( LOWER(1), LOWER(2), CORREL )
     +           - BVU ( LOWER(1), UPPER(2), CORREL )
      ELSE IF ( INFIN(1) .EQ. 2  .AND. INFIN(2) .EQ. 0 ) THEN
         VAL =  BVU ( -UPPER(1), -UPPER(2), CORREL )
     +           - BVU ( -LOWER(1), -UPPER(2), CORREL )
      ELSE IF ( INFIN(1) .EQ. 0  .AND. INFIN(2) .EQ. 2 ) THEN
         VAL =  BVU ( -UPPER(1), -UPPER(2), CORREL )
     +           - BVU ( -UPPER(1), -LOWER(2), CORREL )
      ELSE IF ( INFIN(1) .EQ. 1  .AND. INFIN(2) .EQ. 0 ) THEN
         VAL =  BVU ( LOWER(1), -UPPER(2), -CORREL )
      ELSE IF ( INFIN(1) .EQ. 0  .AND. INFIN(2) .EQ. 1 ) THEN
         VAL =  BVU ( -UPPER(1), LOWER(2), -CORREL )
      ELSE IF ( INFIN(1) .EQ. 1  .AND. INFIN(2) .EQ. 1 ) THEN
         VAL =  BVU ( LOWER(1), LOWER(2), CORREL )
      ELSE IF ( INFIN(1) .EQ. 0  .AND. INFIN(2) .EQ. 0 ) THEN
         VAL =  BVU ( -UPPER(1), -UPPER(2), CORREL )
      END IF
      END  FUNCTION BVNMVN
      FUNCTION BVU( SH, SK, R ) RESULT (VAL)
      IMPLICIT NONE
*
*     A function for computing bivariate normal probabilities.
*
*       Yihong Ge
*       Department of Computer Science and Electrical Engineering
*       Washington State University
*       Pullman, WA 99164-2752
*     and
*       Alan Genz
*       Department of Mathematics
*       Washington State University
*       Pullman, WA 99164-3113
*       Email : alangenz@wsu.edu
*
* BVN - calculate the probability that X is larger than SH and Y is
*       larger than SK. (to accuracy of 1e-16?)
*
* Parameters
*
*   SH  REAL, integration limit
*   SK  REAL, integration limit
*   R   REAL, correlation coefficient
*   LG  INTEGER, number of Gauss Rule Points and Weights
*
      DOUBLE PRECISION, INTENT(IN) :: SH, SK, R
      DOUBLE PRECISION  :: VAL
! Local variables
      DOUBLE PRECISION :: ZERO,HALF,ONE,TWO,TWOPI 
      INTEGER :: I, LG, NG
      PARAMETER ( ZERO = 0.D0,ONE=1D0,TWO=2D0,HALF=0.5D0)
      PARAMETER (TWOPI = 6.2831 85307 179586D0 ) 
      DOUBLE PRECISION, DIMENSION(10,3) :: X, W
      DOUBLE PRECISION :: AS, A, B, C, D, RS, XS
      DOUBLE PRECISION :: SN, ASR, H, K, BS, HS, HK
*     Gauss Legendre Points and Weights, N =  6
      DATA  ( W(I,1), X(I,1), I = 1,3) /
     *  0.1713244923791705D+00,-0.9324695142031522D+00,
     *  0.3607615730481384D+00,-0.6612093864662647D+00,
     *  0.4679139345726904D+00,-0.2386191860831970D+00/
*     Gauss Legendre Points and Weights, N = 12
      DATA ( W(I,2), X(I,2), I = 1,6) /
     *  0.4717533638651177D-01,-0.9815606342467191D+00,
     *  0.1069393259953183D+00,-0.9041172563704750D+00,
     *  0.1600783285433464D+00,-0.7699026741943050D+00,
     *  0.2031674267230659D+00,-0.5873179542866171D+00,
     *  0.2334925365383547D+00,-0.3678314989981802D+00,
     *  0.2491470458134029D+00,-0.1252334085114692D+00/
*     Gauss Legendre Points and Weights, N = 20
      DATA ( W(I,3), X(I,3), I = 1,10) /
     *  0.1761400713915212D-01,-0.9931285991850949D+00,
     *  0.4060142980038694D-01,-0.9639719272779138D+00,
     *  0.6267204833410906D-01,-0.9122344282513259D+00,
     *  0.8327674157670475D-01,-0.8391169718222188D+00,
     *  0.1019301198172404D+00,-0.7463319064601508D+00,
     *  0.1181945319615184D+00,-0.6360536807265150D+00,
     *  0.1316886384491766D+00,-0.5108670019508271D+00,
     *  0.1420961093183821D+00,-0.3737060887154196D+00,
     *  0.1491729864726037D+00,-0.2277858511416451D+00,
     *  0.1527533871307259D+00,-0.7652652113349733D-01/
      SAVE W, X
      IF ( ABS(R) .LT. 0.3D0 ) THEN
         NG = 1
         LG = 3
      ELSE IF ( ABS(R) .LT. 0.75D0 ) THEN
         NG = 2
         LG = 6
      ELSE 
         NG = 3
         LG = 10
      ENDIF
      H = SH
      K = SK 
      HK = H*K
      VAL = ZERO
      IF ( ABS(R) .LT. 0.925D0 ) THEN
         HS = ( H*H + K*K )*HALF
         ASR = ASIN(R)
         DO I = 1, LG
            SN  = SIN(ASR*( X(I,NG)+ONE )*HALF)
            VAL = VAL + W(I,NG)*EXP( ( SN*HK - HS )/( ONE - SN*SN ) )
            SN  = SIN(ASR*(-X(I,NG)+ONE )*HALF)
            VAL = VAL + W(I,NG)*EXP( ( SN*HK - HS )/( ONE - SN*SN ) )
         END DO
         VAL = VAL*ASR*HALF/TWOPI + FI(-H)*FI(-K) 
      ELSE
         IF ( R .LT. ZERO ) THEN
            K  = -K
            HK = -HK
         ENDIF
         IF ( ABS(R) .LT. ONE ) THEN
            AS  = ( ONE - R )*( ONE + R )
            A   = SQRT(AS)
            BS  = ( H - K )**2
            C   = ( 4D0 - HK )/8D0 
            D   = ( 12D0 - HK )/16D0
            VAL = A*EXP( -(BS/AS + HK)*HALF )
     +       *( ONE - C*(BS - AS)*(ONE - D*BS/5D0)/3D0 + C*D*AS*AS/5D0 )
            IF ( HK .GT. -160D0 ) THEN
               B   = SQRT(BS)
               VAL = VAL - EXP(-HK/2D0)*SQRT(TWOPI)*FI(-B/A)*B
     +                    *( ONE - C*BS*( ONE - D*BS/5D0 )/3D0 ) 
            ENDIF
            A = A/2D0
            DO I = 1, LG
               XS  = ( A*(X(I,NG)+ONE) )**2
               RS  = SQRT( ONE - XS )
               VAL = VAL + A*W(I,NG)*
     +             ( EXP( -BS/(2D0*XS) - HK/(ONE+RS) )/RS 
     +            - EXP( -(BS/XS+HK)/2D0 )*( ONE + C*XS*( ONE + D*XS )))
               XS  = AS*(ONE-X(I,NG))**2/4D0
               RS  = SQRT( ONE - XS )
               VAL = VAL + A*W(I,NG)*EXP( -(BS/XS + HK)/2D0 )
     +                    *( EXP( -HK*(ONE-RS)/(2D0*(ONE+RS)) )/RS 
     +                       - ( ONE + C*XS*( ONE + D*XS ) ) )
            END DO
            VAL = -VAL/TWOPI
         ENDIF
         IF ( R .GT. ZERO ) VAL =  VAL + FI( -MAX( H, K ) )
         IF ( R .LT. ZERO ) VAL = -VAL + MAX( ZERO, FI(-H)-FI(-K) )     
      ENDIF
      RETURN
      END FUNCTION BVU

      
      
      FUNCTION THL(H1, L1, GH) RESULT (T)                          
      USE GLOBALDATA, ONLY : ABSEPS
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in) :: H1, L1, GH
      DOUBLE PRECISION             :: T                                     
!                                                                       
!     Local constants                                                   
!     
      DOUBLE PRECISION, PARAMETER :: SQPI = 1.77245385090552d0   !=sqrt(pi)
      DOUBLE PRECISION, PARAMETER :: SQTWO= 1.41421356237310d0    !=sqrt(2)
      DOUBLE PRECISION, PARAMETER :: PI1  = 0.31830988618379d0     !=1/pi
      DOUBLE PRECISION, PARAMETER :: PI   = 3.14159265358979D0     !=pi
      DOUBLE PRECISION, PARAMETER :: ZERO = 0.D0
      DOUBLE PRECISION, PARAMETER :: QUART= 0.25D0
      DOUBLE PRECISION, PARAMETER :: HALF = 0.5D0
      DOUBLE PRECISION, PARAMETER :: ONE  = 1.D0
      DOUBLE PRECISION, PARAMETER :: TWO  = 2.D0
      DOUBLE PRECISION, PARAMETER :: EXPLIM=80.D0
      INTEGER,          PARAMETER :: NCk=7
      DOUBLE PRECISION, DIMENSION(NCk) :: Ck
! The following parameters gives absolute error <= 1e-7
      PARAMETER ( Ck =(/ 0.99999124672847d0, -0.99913342611804d0,
     &     0.98540829635096d0, -0.90305751774767d0, 0.66776957888038d0,
     &     -0.32268060187487d0, 0.07170680041655d0/))
! The following parameters gives absolute error <= 1e-11.
!      PARAMETER ( Ck =(/ 0.99999995579043d0, -0.99999111387630d0, 
!     &     0.99969933071259d0,-0.99596450623251d0, 0.97168814716964d0,  
!     &    -0.88105640681133d0, 0.67507517897079d0,-0.38534334229836d0,
!     &     0.13907128134639d0,-0.02317854687612d0/))
    
! Local variables
      INTEGER :: k
      DOUBLE PRECISION :: EPSL,H,L, H2, CON, CN, 
     &   L2, EX, W2, AP, S2, SP, S1, SN, SGN, H2L, Ik

! THL computes the integral
!             L
!   T(h,L) = int exp(h^2/2*(1-x^2))/(1+x^2) dx /(2*pi)
!             0
! 
!   GH = FI(h)
!
!     TLH is a controlled precision Fortran function 
!     which may be used to calculate    
!     bivariate normal probabilities, e.g. the CDF  for   
!     two normal variates X and Y whose correlation is R:
!
!     Prob(X<h,Y<k;R)=(FI(h)+FI(k))/2-T(h,g(h,k,R))-T(k,g(k,h,R))-J(h,k)/2
!  where 
!       J(h,k)   = 0     if h*k>0 or hk=0 and h+k>=0
!                  1     otherwise
!       g(h,k,R) = (h/k-R)/sqrt(1-R^2)
!
!     The accuracy is specified with ABSEPS.  
!
! Some properties of T(H,L):
! T(h,L) = 1/4-(FI(h)-1/2)*(FI(h*L)-1/2)-T(h*L,1/L)
! T(h,L) = T(-h,L)=-T(h,-L)
! T(h,0) = 0
! T(0,L) = atan(L)/(2*pi)
! T(h,1) = FI(h)*(1-FI(h))/2

! References
!  Jagdish K. Patel and Campbell B. Read (1982)
! "Handbook of the normal distribution",
! marcel dekker inc,New York - Basel, Vol. 40, pp 293--300       
!  
! Tested on: Matlab 5.1                                              
! History:
! revised pab 19.04.2000
!   - found a bug when L<-1, now fixed
! pab 20.03.2000
! - fixed some errors, optimized the splitting between Owen and Borth algorithms.
! - added Daley (1974) approximation
! pab 03.03.2000
!   - forgot to pass the calculated value for T on to output 
!     => this is probably why  THL gave floating invalid earlier
! by pab 19.08.199

!      print *, 'THL enter,L,H',L1,H1
      T    = ZERO
      EPSL = ABSEPS*1.d-1
      L    = L1
      H    = H1
      if (abs(L).LE.EPSL) RETURN ! T(H,0)=0
      
      if (abs(H).LE.EPSL) THEN   !T(0,L)
         T=ATAN(L)*HALF*PI1
         return
      endif
      if (abs(abs(L)-ONE).LE.EPSL) THEN ! T(H,1)
         T=HALF*GH*(ONE-GH)*SIGN(ONE,L)
         return
      end if   
      !print *, 'THL still here'
      SGN = ONE
      if (ABS(L).GT. ONE) THEN  ! make sure abs(L)<=1 to avoid numerical problems
                                ! using the identities:
                                ! T(h,L)=-T(h,-L)
                                ! T(h,L)=1/4+(1/2-FI(h))*(FI(h*L)-1/2)-T(h*L,1/L)
         IF (L.GE.0) SGN =-ONE
         H   = H * L                                                       
         L   = ONE /ABS( L)
         T   = (HALF-GH)*(FI(H)-HALF)-SGN*QUART
      endif

      H2 = H * H * HALF                                                                                                                                                                    
      IF (H2 .GE. EXPLIM) RETURN  ! THL is found to the desired accuracy
   
      IF (EPSL.GE.5.5d-5) THEN !Daley (1974) approximation
          SN=ATAN(L)
          H2L=L*L*H2*TWO
          H2L=H2L*H2L
          T = T + SGN*HALF*PI1*SN*EXP(-L*H2/SN)*(ONE+0.00868d0*H2L)
          RETURN
      ENDIF
      EX = EXP(-H2) 
      CON = abs(TWO*PI*EPSL ) 

! The following division between Owens and Borths algorithm is optimal for nearly all cases of EPSL>=1e-8
! This division is done in order to avoid slow convergence when L approx +/-1
!
!           Owens (1956) Algorithm:
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IF (abs(H).LE.TWO.AND. 
     &     abs(L).LE.MAX(-0.72d0*abs(H)+1.936d0,HALF)) THEN

         L2 = L * L
         W2 = H2 * EX                                                       
         AP = ONE; SP=ONE
         S1 = ONE - EX;  S2=S1;    CN=S1
                                                      
         ! Alternating series: check if we can stop before SP=30 (usually SP is less than 8)         
         do while ((abs(CN*L).GT.CON).AND.( SP.LT.30.d0))
            SN =  SP                                                            
            SP =  SP + ONE                                                      
            S2 =  S2 - W2                                                       
            W2 =  W2 * H2 / SP                                                  
            AP = -AP * L2                                                        
            CN =  AP * S2 / (SN + SP)                                           
            S1 =  S1 + CN
         end do
!         print *,'THL, SP=',SP
                                !*, 'THL leaving 1.6' 
         T = T + SGN*(atan(L) - L*S1) *HALF*PI1                                                   
         return
      endif
! Otherwise
!
!           Borth (1974) Algorithm:
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     
! approximate 1/(1+x^2) with a Tchebyshev polynomial of degree 2K
! and then integrate. Ck contains the coefficients of the polynomial.
! The absolute error is less than 1e-7.
        
      H2L=H2*L*L
      SN=SQTWO/H
      W2=exp(-H2L)*H*L/SQTWO ! safer when L approx 0
      Ik=SQPI*(FI(H*L)-HALF)
      
   
      CN=Ck(1)*Ik*SN
      S1=CN
      k=1
      !Alternating series: Check if we can stop before k=Nck
      do while ((abs(CN*EX).GT.CON) .AND. (k.LT.NCk)) 
         Ik=HALF*(DBLE(2*k-1)*Ik-W2)
         W2=W2*H2L
         SN=SN/H2
         k=k+1
         CN=Ck(k)*Ik*SN
         S1=S1+CN   
      end do	
!      print *,'THL K=',K                      
      T = T + SGN*EX*S1*HALF*PI1
      !print *, 'THL leaving last'
      RETURN                                                             
      END FUNCTION THL                                      

      FUNCTION NORM2DPRB(a1,b1,a2,b2,R) RESULT (prb)
      USE GLOBALDATA, ONLY : ABSEPS,CFxCutoff,xCutoff
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in) :: a1,b1,a2,b2,R
      DOUBLE PRECISION             :: prb
! Local variables
      DOUBLE PRECISION, PARAMETER :: PI1  = 0.31830988618379d0 !=1/pi
      DOUBLE PRECISION :: L1,L2,L3,L4,Ga1,Ga2,Gb1,Gb2,SQR,
     &     T11,T12,T21,T22,T31,T32,T41,T42,EPSL

!NORM2DPRB calculates the probability Prob(a1<X1<b1,a2<x2<b2)
!          for a bivariate Gaussian distribution
!          with correlation R
!
! References
!  Jagdish K. Patel and Campbell B. Read (1982)
! "Handbook of the normal distribution",
! marcel dekker inc,New York - Basel, Vol. 40, pp 293--300
! Tested on: Matlab 5.1


! History:
! pab 03.03.2000
!   - found  the error in THL which gave floating invalid
! by pab 19.08.1999
! NOTE!!! this does not work correctly yet. For some mysterious reason
! it crashes from time to time giving the message floating invalid
      EPSL=ABSEPS*1.d-1 
      if ((a1+EPSL.GE.b1).OR.(a2+EPSL.GE.b2)) THEN
         prb=0.d0
         return
      endif	
      Ga1=FI(a1)
      Gb1=FI(b1)
                                ! See if area is symmetric
                                ! The following may be optimized further
                                ! by checking against xCutOff
      if (abs(a1-a2).LE.EPSL) THEN
         Ga2=Ga1
      elseif (abs(a2-b1).LE.EPSL) THEN
         Ga2=Gb1
      else   
         Ga2=FI(a2)
      endif

      if (abs(b1-b2).LE.EPSL) THEN
         Gb2=Gb1
      elseif (abs(b2-a1).LE.EPSL) THEN
         Gb2=Ga1
      else   
         Gb2=FI(b2)
      endif	
      

      if (abs(R).LE.EPSL) then ! R=0
         prb=(Gb1-Ga1)*(Gb2-Ga2)
         GO TO 200
      endif

      if (abs(R)+EPSL.GE.1.d0) THEN
!         IF (abs(R).GT.1.d0+EPSL) PRINT *,'Warning: correlation, R=',R
         if (R.GT.0.d0) then    ! R = 1
            prb=min(Gb1,Gb2)+min(Ga1,Ga2)-min(Ga1,Gb2)-min(Gb1,Ga2)
         else   ! R =-1
            prb=max(Gb1+Gb2-1.d0,0.d0)+max(Ga1+Ga2-1.d0,0.d0)-
     &           max(Ga1+Gb2-1.d0,0.d0)-max(Gb1+Ga2-1.d0,0.d0)
         endif
         GO TO 200
      endif

      SQR = sqrt(1.d0-R*R)
        !print *, 'norm2d enter',ga1,gb1,ga2,gb2
      !print *, 'norm2d enter',a1,b1,a2,b2,R
      ! L1=PSI(b1,b2,R),! L3=PSI(b1,a2,R)
      if (abs(b1).LE.EPSL) then
         if (abs(b2).LE.EPSL) then
            L1=0.25d0+asin(R)*0.5d0*PI1
         else
            T11=SIGN(0.25d0,b2)
            T12=THL(b2,-R/SQR,Gb2)
            L1=0.5d0*Gb2+0.25d0-T11-T12+min(0.d0,SIGN(0.5d0,b2))
         endif
         if (abs(a2).LE.EPSL) then
            L3=0.25d0+asin(R)*0.5d0*PI1
         else
            T31=SIGN(0.25d0,a2)
            T32=THL(a2,-R/SQR,Ga2)
            L3=0.5d0*Ga2+0.25d0-T31-T32+min(0.d0,SIGN(0.5d0,a2))
         endif
      else
         !print *, 'norm2d'
         if (abs(b2).LE.EPSL) then
            T11=THL(b1,-R/SQR,Gb1)
            T12=SIGN(0.25d0,b1)
            L1=0.5d0*Gb1+0.25d0-T11-T12+0.5d0*min(0.d0,SIGN(1.d0,b1))
         else
            !print *, 'norm2d,b1,r,gb1',b1,(b2/b1-R)/SQR,gb1
            T11=THL(b1,(b2/b1-R)/SQR,Gb1)
            !print *, 'norm2d,T11',T11
            if (abs(b2-b1).LE.EPSL) then
               L1=Gb1-2.d0*T11
               T12=T11
            else  
               !print *,'norm2d T12'
               if (abs(b2+b1).LE.EPSL) then
                  T12=T11
               else
                  T12=THL(b2,(b1/b2-R)/SQR,Gb2)
               endif
               !print *,'norm2d T12',T12
               L1=0.5d0*(Gb1+Gb2)-T11-T12+min(0.d0,SIGN(0.5d0,b1*b2))
            endif
            !print *, 'norm2d,L1',L1
         endif
         !print *, 'norm2d'
         if (abs(a2).LE.EPSL ) then
            T32=SIGN(0.25d0,b1)
            T31=THL(b1,-R/SQR,Gb1)
            L3=0.5d0*Gb1+0.25d0-T31-T32+min(0.d0,SIGN(0.5d0,b1))
         else
           ! print *, 'norm2d,b1,r,gb1',b1,(a2/b1-R)/SQR,gb1
            T31=THL(b1,(a2/b1-R)/SQR,Gb1)
            !print *, 'norm2d,T31',T31
            if (abs(a2-b1).LE.EPSL) THEN
               L3=Ga1-2.d0*T31
               T32=T31
            else  
               IF (abs(a2+b1).LE.EPSL) THEN
                  T32=T31
               ELSE
                  T32=THL(a2,(b1/a2-R)/SQR,Ga2)
               ENDIF
               if ((abs(a2+b2).LE.EPSL).AND.(abs(a1+b1).LE.EPSL)) then
                  prb=2.d0*(T31+T32-T11-T12)+1.d0 ! OK
                                !prb=2.d0*(L1-L3)-1.d0
                  goto 200
               endif
               L3=0.5d0*(Gb1+Ga2)-T31-T32+min(0.d0,SIGN(0.5d0,b1*a2))        
            endif	
         endif   
      endif
      !print *, 'norm2d L1,L3',L1,L3
      !L2=PSI(a1,a2,R) L4=PSI(a1,b2,R)
      if (abs(a1).LE.EPSL) then
         if (abs(b2).LE.EPSL) then
            L4=0.25d0+asin(R)*0.5d0*PI1
         else
            T41=SIGN(0.25d0,b2)
            T42=THL(b2,-R/SQR,Gb2)
            L4=0.5d0*Gb2+0.25d0-T41-T42+min(0.d0,SIGN(0.5d0,b2))
         endif
         if (abs(a2).LE.EPSL) then
            L2=0.25d0+asin(R)*0.5d0*PI1
         else
            T21=SIGN(0.25d0,a2)
            T22=THL(a2,-R/SQR,Ga2)
            L2=0.5d0*Ga2+0.25d0-T21-T22+min(0.d0,SIGN(0.5d0,a2))
         endif
      else
         if (abs(b2).LE.EPSL ) then
            T41=THL(a1,-R/SQR,Ga1)
            T42=SIGN(0.25d0,a1)
            L4=0.5d0*Ga1+0.25d0-T41-T42+min(0.d0,SIGN(0.5d0,a1))
         else
            !print *, 'norm2d T41, L',(b2/a1-R)/SQR
            T41=THL(a1,(b2/a1-R)/SQR,Ga1)
            if (abs(b2-a1).LE.EPSL) then
               L4=Gb1-2.d0*T41
               T42=T41
            else    
               if (abs(b2+a1).LE.EPSL) then
                  T42=T41
               else
                  T42=THL(b2,(a1/b2-R)/SQR,Gb2)
               endif
               L4=0.5d0*(Ga1+Gb2)-T41-T42+min(0.d0,SIGN(0.5d0,a1*b2))
            endif	
         endif
         !print *, 'norm2d L4',L4
         if (abs(a2).LE.EPSL) then
            T22=SIGN(0.25d0,a1)
            T21=THL(a2,-R/SQR,Ga1)
            L2=0.5d0*Ga1+0.25d0-T21-T22+min(0.d0,SIGN(0.5d0,a1))
         else
            T21=THL(a1,(a2/a1-R)/SQR,Ga1)
            if (abs(a2-a1).LE.EPSL) then
               L2=Ga1-2.d0*T21
               T22=T21
            else   
               if (abs(a2+a1).LE.EPSL) then
                  T22=T21
               else
                  T22=THL(a2,(a1/a2-R)/SQR,Ga2)
               endif
               L2=0.5d0*(Ga1+Ga2)-T21-T22+min(0.d0,SIGN(0.5d0,a1*a2))
            endif	
         endif
         !print *, 'norm2d L2',L2
      endif

                                !L1,L2,L3,L4
      prb=L1+L2-L3-L4      
 200  if (prb > 1.d0) then ! fix up round off error
         prb=1.d0
      elseif (prb<0.d0) then 
            prb=0.d0
      end if
      !print *, 'norm2d leaving'
      return
      END FUNCTION NORM2DPRB
      
      FUNCTION GAUSINT (X1, X2, A, B, C, D) RESULT (value)               
      USE GLOBALDATA,ONLY:  xCutOff
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in) :: X1,X2,A,B,C,D
      DOUBLE PRECISION             :: value
!Local variables
      DOUBLE PRECISION             :: Y1,Y2,Y3
       DOUBLE PRECISION, PARAMETER :: SQTWOPI1=3.9894228040143d-1 !=1/sqrt(2*pi)
      ! Let  X  be standardized Gaussian variable, 
      ! i.e., X=N(0,1). The function calculate the
      !  following integral E[I(X1<X<X2)(A+BX)(C+DX)
      ! where I(X1<X<X2) is an indicator function of
      ! the set {X1<X<X2}. 
      IF (X1.GE.X2) THEN                                                 
         value = 0.d0                                                    
         RETURN                                                          
      ENDIF
      IF (ABS (X1) .GT.xCutOff) THEN                                           
         Y1 = 0.d0                                                         
      ELSE                                                               
         Y1 = (A * D+B * C + X1 * B * D) * EXP ( - 0.5d0 * X1 * X1)        
      ENDIF
      IF (ABS (X2) .GT.xCutOff) THEN                                           
         Y2 = 0.d0                                                         
      ELSE                                                               
         Y2 = (A * D+B * C + X2 * B * D) * EXP ( - 0.5d0 * X2 * X2)        
      ENDIF
      Y3 = (A * C + B * D) * (FI (X2) - FI (X1) )                        
      value = Y3 + SQTWOPI1 * (Y1 - Y2)                                      
      RETURN                                                             
      END FUNCTION GAUSINT
      
      
      FUNCTION GAUSINT2 (X1, X2, A, B) RESULT (value)               
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in) :: X1,X2,A,B
      DOUBLE PRECISION             :: value
! Local variables
      DOUBLE PRECISION             :: X0,Y0,Y1,Y2
      DOUBLE PRECISION, PARAMETER :: SQTWOPI1=3.9894228040143d-1 !=1/sqrt(2*pi)
!     Let  X  be standardized Gaussian variable, 
!     i.e., X=N(0,1). The function calculate the
!     following integral E[I(X1<X<X2)ABS(A+BX)]
!     where I(X1<X<X2) is an indicator function of
!     the set {X1<X<X2}. 
      IF (X1.GE.X2) THEN                                                 
         value = 0.d0                                                    
         RETURN                                                          
      ENDIF
      IF (ABS(B).EQ.0.d0) THEN
         value = ABS(A)*(FI(X2)-FI(X1))
         RETURN
      ENDIF
     
      Y1 = -A*FI(X1)+SQTWOPI1*B*EXP(-0.5d0*X1*X1)
      Y2 = A*FI(X2)-SQTWOPI1*B*EXP(-0.5d0*X2*X2)
      IF ((B*X1.LT.-A).AND.(-A.LT.B*X2))THEN
         X0 = -A/B  
         Y0 = 2.d0*(A*FI(X0)-SQTWOPI1*B*EXP(-0.5d0*X0*X0))
         value=ABS(Y2-Y1-Y0)
      ELSE
         value=ABS(Y1+Y2) 
      ENDIF                                       
      RETURN                                                             
      END FUNCTION GAUSINT2

      SUBROUTINE EXLMS(A, X1, X2, INFIN, LOWER, UPPER, Ca,Pa)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in) :: A, X1, X2
      DOUBLE PRECISION, INTENT(out) :: LOWER, UPPER,Ca,Pa
      INTEGER,INTENT(in) :: INFIN
      ! Local variables
      DOUBLE PRECISION :: P1
      DOUBLE PRECISION, PARAMETER :: AMAX = 5D15
      
      ! Let  X  be standardized Gaussian variable, 
      ! i.e., X=N(0,1). The function calculate the
      !  following integral E[I(X1<X<X2)ABS(A+X)]/C(A) 
      ! where I(X1<X<X2) is an indicator function of
      ! the set {X1<X<X2} and C(A) is a normalization factor
      ! i.e. E(I(-inf<X<inf)*ABS(A+X)) = C(A) 
      ! 
      ! Pa = probability at the inflection point of the CDF
      ! Ca = C(A) normalization parameter
      ! A  = location parameter of the inflection point
      ! X1,X2 = integration limits < infinity
      LOWER = 0.d0
      UPPER = 1.d0
      
      P1 = EXFUN(-A,A)
      Ca = A+2D0*P1
      Pa = P1/Ca 
      IF ( INFIN .LT. 0 ) RETURN
      IF ( INFIN .NE. 0 ) THEN
         IF (X1.LE.-A) THEN
            LOWER = EXFUN(X1,A)/Ca
         ELSE
            LOWER = 1D0-EXFUN(-X1,-A)/Ca
         ENDIF
      ENDIF
      IF ( INFIN .NE. 1 ) THEN
         IF (X2.LE.-A) THEN
            UPPER = EXFUN(X2,A)/Ca
         ELSE
            UPPER = 1D0-EXFUN(-X2,-A)/Ca
         ENDIF
      ENDIF
      
      RETURN
      CONTAINS
      FUNCTION EXFUN(X,A) RESULT (P)
      DOUBLE PRECISION, INTENT(in) :: X,A
      DOUBLE PRECISION :: P
      DOUBLE PRECISION, PARAMETER :: SQTWOPI1=3.9894228040143d-1 !=1/sqrt(2*pi)
      P = EXP(-X*X*0.5D0)*SQTWOPI1-A*FI(X)
      END FUNCTION EXFUN
      END SUBROUTINE EXLMS 


      FUNCTION EXINV(P,A,Ca,Pa) RESULT (VAL)
!EXINV calculates the inverse of the CDF of abs(x+A)*exp(-x^2/2)/SQRT(2*pi)/C(A)
!      where C(A) is a normalization parameter depending on A
!
!   CALL:  val   = exinv(P,A,Ca,Pa)
! 
!     val  = quantiles
!     p    = probabilites
!     A    = location parameter
!     Ca   = normalization parameter
!     Pa   = excdf(-A,A) probability at the inflection point of the CDF
      double precision, intent(in) :: P,A,Ca,Pa
      double precision :: val
! local variables
      double precision, parameter :: amax = 5.D15
      double precision, parameter :: epsl = 5.D-15
      double precision, parameter :: xmax = 8
      double precision :: P1,Xk,Ak,Zk,SGN

!      if (P<0.D0.OR.P.GT.1.D0) PRINT *,'warning P<0 or P>1'

! The inverse cdf of 0 is -inf, and the inverse cdf of 1 is inf.  
      if (P.LE.EPSL.OR.P+EPSL.GE.1.D0) THEN
         VAL = SIGN(xmax,P-0.5D0)
         return
      endif
      Ak = ABS(A)
      if (EPSL.LT.Ak .AND. Ak.LT.amax) THEN
         IF (ABS(p-Pa).LE.EPSL) THEN
            VAL = SIGN(MIN(Ak,xmax),-A)
            RETURN
         ENDIF
         IF (Ak.LT.1D-2) THEN   ! starting guess always less than 0.2 from the true value
            IF (P.GE.0.5D0) THEN 
               xk = SQRT(-2D0*log(2D0*(1D0-P)))
            ELSE
               xk = -SQRT(-2D0*log(2D0*P))
            ENDIF
         ELSE
            xk = FIINV(P)       ! starting guess always less than 0.8 from the true value
            ! Modify starting guess if possible in order to speed up Newtons method
            IF (1D-3.LE.P.AND. P.LE.0.99D0.AND. 
     &           3.5.LE.Ak.AND.Ak.LE.1D3 ) THEN
               SGN = SIGN(1.d0,-A)
               Zk = xk*SGN
               xk = SGN*(Zk+((1D0/(64.9495D0*Ak-178.3191D0)-0.02D0/Ak)*
     &              Zk+1D0/(-0.99679234298211D0*Ak-0.07195350071872D0))/
     &            (Zk/(-1.48430620263825D0*Ak-0.33340759016175D0)+1D0))
            ELSEIF ((P.LT.1D-3.AND.A.LE.-3.5D0).OR. 
     &              (3.5D0.LE.A.AND.P.GT.0.99D0)) THEN
               SGN = SIGN(1.d0,-A)
               Zk = xk*SGN
               P1 = -2.00126182192701D0*Ak-2.57306603933111D0
               xk = SGN*Zk*(1D0+
     &              P1/((-0.99179258785909D0*Ak-0.21359746002397D0)*
     &              (Zk+P1)))
            ENDIF
         ENDIF
         ! Check if the starting guess is on the correct side of the inflection point
         IF (xk.LE.-A .AND. P.GT.Pa) xk = 1.D-2-A 
         IF (xk.GE.-A .AND. P.LT.Pa) xk = -1.D-2-A

      
         IF (P.LT.Pa) THEN
            VAL = funca(xk,A,P*Ca);
         ELSE  ! exploit the symmetry of the CDF
            VAL = -funca(-xk,-A,(1.D0-P)*Ca)
         ENDIF
         
      ELSEIF (ABS(A).LE.EPSL) THEN
         IF (P>=0.5D0) THEN
            VAL = SQRT(-2D0*log(2D0*(1.D0-P)))
         ELSE
            VAL = -SQRT(-2D0*log(2D0*P));
         ENDIF
      ELSE  ! ABS(A) > AMAX
         VAL = FIINV(P)
      ENDIF
      !CALL EXLMS(A,0.d0,VAL,0,ak,P1,zk,sgn) 
      !If (ABS(p-P1).GT.0.0001) PRINT *,'excdf(x,a)-p',p-P1 
      RETURN
      
      CONTAINS

      function funca(xk0,ak,CaP) RESULT (xk)
      double precision, intent(in) :: xk0,ak,CaP ! =Ca*P
      DOUBLE PRECISION :: xk
!Local variables                                
      INTEGER,          PARAMETER :: ixmax = 25
      double precision, parameter :: crit = 7.1D-08 ! = sqrt(1e-15)
      double precision, parameter :: SQTWOPI1 = 0.39894228040143D0 !=1/SQRT(2*pi)
      double precision, parameter :: SQTWOPI = 2.50662827463100D0 !=SQRT(2*pi)
      INTEGER :: IX
      DOUBLE PRECISION :: H,H1,tmp0,tmp1,XNEW
      ! Newton's Method or Fixed point iteration to find the inverse of the EXCDF.
      ! Assumption: xk0 < -ak and xk < -ak      
      ! Permit no more than IXMAX iterations.
      IX = 0
      H  = 1.D0
      xk = xk0    ! starting guess for the iteration
     
     
                                ! Break out of the iteration loop for the following:
                                !  1) The last update is very small (compared to x).
                                !  2) The last update is very small (compared to sqrt(eps)=crit).
                                !  3) There are more than 15 iterations. This should NEVER happen. 
      IF (.TRUE..OR.ABS(ak).LT.1.D-2) THEN 
      ! Newton's method
      !~~~~~~~~~~~~~~~~~
      DO WHILE( ABS(H).GT.MIN(crit*ABS(xk),crit).AND.IX.LT.IXMAX) 
         
         IX = IX+1   
                                !print *,'Iteration ',IX
         
         tmp0  = FI(xk)
         tmp1  = EXP(-xk*xk*0.5D0)*SQTWOPI1 ! =normpdf(x)
         H1 = (tmp1-ak*tmp0-CaP)/(ABS(xk+ak)*tmp1)
         H  = DSIGN(MIN(ABS(H1),0.7D0/DBLE(IX)),H1) ! Only allow smaller and smaller steps
         
         xnew = xk - H
                                ! Make sure that the current guess is less than -a.
                                ! When Newton's Method suggests steps that lead to -a guesses
                                ! take a step 9/10ths of the way to -a:
         IF (xnew.GT.-ak-crit) THEN
            xnew = (xk - 9.D0*ak)*1D-1;
            H    = xnew - xk
         ENDIF
         xk = xnew
      END DO
      ELSE                      ! FIXED POINT iteration
                                !~~~~~~~~~~~~~~~~~~~~~~~
         DO WHILE (ABS(H).GT.MIN(crit*ABS(xk),crit).AND.IX.LT.IXMAX)                                
            IX   = IX+1   
            tmp0 = SQTWOPI1*EXP(-xk*xk*0.5D0)/FI(xk)
            tmp1 = -2.D0*LOG(SQTWOPI*CaP*tmp0/(tmp0-ak))
            SGN  = sign(1.D0,tmp1)
            xnew = -SQRT(SGN*tmp1)*SGN
                                ! Make sure that the current guess is less than -a.
                                ! When this method suggests steps that lead to -a guesses
                                ! take a step 9/10ths of the way to -a:
            IF (xnew.GT.-ak-crit) xnew = (xk - 9.D0*ak)*1.D-1

            H  = xnew - xk
            xk = xnew
         END DO
      ENDIF

      !print *,'EXINV total number of iterations ',IX
      if (IX.GE.IXMAX) THEN 
!         print *, 'Warning: EXINV did not converge. Cap=',Cap
!         print *, 'The last step was:  ', h, ' value=,',xk,' ak=',ak
      endif
      return
      END FUNCTION FUNCA
      END FUNCTION EXINV
      END MODULE FIMOD

      MODULE SWAPMOD
      INTERFACE SWAP
      MODULE PROCEDURE SWAP_R, SWAP_I, SWAP_C
      END INTERFACE
      CONTAINS

      SUBROUTINE SWAP_R(A,B)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT (INOUT) :: A, B
      DOUBLE PRECISION                 :: TEMP
      TEMP = A; A = B; B = TEMP
      END SUBROUTINE SWAP_R

      SUBROUTINE SWAP_I(A,B)
      IMPLICIT NONE
      INTEGER, INTENT (INOUT) :: A, B
      INTEGER                 :: TEMP
      TEMP = A ; A = B ; B = TEMP
      END SUBROUTINE SWAP_I

      SUBROUTINE SWAP_C(A,B)
      IMPLICIT NONE
      CHARACTER, INTENT (INOUT) :: A, B
      CHARACTER                 :: TEMP
      TEMP = A ; A = B ; B = TEMP
      END SUBROUTINE SWAP_C
      END MODULE SWAPMOD

      MODULE COVSORTMOD
      IMPLICIT NONE

      INTERFACE RCSCALE
      MODULE PROCEDURE  RCSCALE
      END INTERFACE

      INTERFACE RCSWAP
      MODULE PROCEDURE RCSWAP
      END INTERFACE

      INTERFACE CVSRTXC
      MODULE PROCEDURE CVSRTXC
      END INTERFACE

      INTERFACE COVSRT
      MODULE PROCEDURE COVSRT
      END INTERFACE

      INTERFACE COVSRT1
      MODULE PROCEDURE COVSRT1
      END INTERFACE

      INTERFACE C1C2
      MODULE PROCEDURE C1C2
      END INTERFACE
      
      INTERFACE PRINTCOF
      MODULE PROCEDURE PRINTCOF
      END INTERFACE 

	INTERFACE mexPrintf
	MODULE PROCEDURE mexPrintf
      END INTERFACE 

      INTERFACE PRINTVAR
      MODULE PROCEDURE PRINTVARD, PRINTVARI
      END INTERFACE 

      INTERFACE PRINTVEC
      MODULE PROCEDURE PRINTVECD, PRINTVECI
      END INTERFACE 

      INTERFACE ADJLIMITS
      MODULE PROCEDURE ADJLIMITS
      END INTERFACE 
      CONTAINS
      SUBROUTINE ADJLIMITS(A,B, infi)
      USE GLOBALDATA, ONLY : xCutOff
      IMPLICIT NONE
!     Adjust INFI when integration limits A and/or B is too far out in the tail
      DOUBLE PRECISION, INTENT(IN)     :: A,B
      INTEGER,          INTENT(IN OUT) :: infi
!      DOUBLE PRECISION, PARAMETER :: xCutOff = 37.D0
      IF (infi.GE.0) THEN
         IF (infi.NE.0)THEN
            IF (A .LT. -xCutOff) THEN
               infi = infi-2
               !CALL mexprintf('ADJ A')
            ENDIF
         ENDIF
         IF (ABS(infi).NE.1) THEN
            IF (xCutOff .LT. B) THEN
               infi = infi-1
               !CALL mexprintf('ADJ B')
            ENDIF
         END IF
      END IF
      RETURN
      END SUBROUTINE ADJLIMITS

	SUBROUTINE mexPrintf(string1)
      CHARACTER*80,INTENT(IN) :: string1
      PRINT  string1 
      RETURN
      END SUBROUTINE mexPrintf


      SUBROUTINE PRINTVARI(I,TXT)
      INTEGER,  INTENT(IN) :: I
      CHARACTER*80, OPTIONAL, INTENT(IN) :: TXT
      CHARACTER*80 :: string1
      IF (PRESENT(TXT)) THEN
       !  CALL mexprintf(TXT)
       !  CALL mexprintf(':')
      ENDIF
      WRITE(string1,125) I
 125  FORMAT (' ',i5) 
      !CALL mexprintf(string1//CHAR(10)) 
      RETURN
      END SUBROUTINE PRINTVARI
      SUBROUTINE PRINTVARD(D,TXT)
      DOUBLE PRECISION,  INTENT(IN) :: D
      CHARACTER*80, OPTIONAL, INTENT(IN) :: TXT
      CHARACTER*80 :: string1
      IF (PRESENT(TXT)) THEN
       !  CALL mexprintf(TXT)
        ! CALL mexprintf(':')
      ENDIF
      WRITE(string1,115) D
 115  FORMAT (' ',10F10.5) 
      !CALL mexprintf(string1) 
      !CALL mexprintf(CHAR(10))
      RETURN
      END SUBROUTINE PRINTVARD

      SUBROUTINE PRINTVECD(CDI,TXT)
      DOUBLE PRECISION, DIMENSION(:),INTENT(in) :: CDI
      CHARACTER*80, OPTIONAL, INTENT(IN) :: TXT
      INTEGER :: I
      CHARACTER*80 :: string
      IF (PRESENT(TXT)) THEN
       !  CALL mexprintf(TXT) 
        ! CALL mexprintf(':')
      ENDIF
      DO I = 1, SIZE(CDI,1)
         WRITE(string,115) CDI(I)
 115     FORMAT (' ',10F10.5) 
         !CALL mexprintf(string) 
      ENDDO
      !CALL mexprintf(CHAR(10)) 
      RETURN
      END SUBROUTINE PRINTVECD
      SUBROUTINE PRINTVECI(CDI,TXT)
      INTEGER, DIMENSION(:),INTENT(in) :: CDI
      CHARACTER*80, OPTIONAL, INTENT(IN) :: TXT
      INTEGER :: I
      CHARACTER*80 :: string
      IF (PRESENT(TXT)) THEN
       !  CALL mexprintf(TXT) 
        ! CALL mexprintf(':')
      ENDIF
      DO I = 1, SIZE(CDI,1)
         WRITE(string,115) CDI(I)
 115     FORMAT (' ',i5) 
         !CALL mexprintf(string) 
      ENDDO
      !CALL mexprintf(CHAR(10)) 
      RETURN
      END SUBROUTINE PRINTVECI

      SUBROUTINE PRINTCOF(Ntd,A,B,INFI,COF,INDEX1)
      IMPLICIT NONE
      INTEGER ,INTENT(in) :: Ntd
      DOUBLE PRECISION, DIMENSION(:),INTENT(in) :: A,B
      DOUBLE PRECISION, DIMENSION(:,:),INTENT(in) :: COF
      INTEGER , DIMENSION(:), INTENT(in) :: INFI,INDEX1
! Local variables
      INTEGER :: I, J, K
      CHARACTER*80 :: string
      CALL mexprintf('  Lower  Upper  Cholesky Matrix '//CHAR(10))
      
      J = Ntd !MIN(Ntd,5)
      DO I = 1, Ntd   
         IF ( INFI(I) .LT. 0 ) THEN 
            WRITE(string,111) INDEX1(I)
 111        FORMAT (i2,'  -inf    inf ') 
       !     CALL mexprintf(string) 
         ELSE IF ( INFI(I) .EQ. 0 ) THEN 
            WRITE(string,112) INDEX1(I),B(I)
 112        FORMAT (i2,'  -inf ', 5F5.2)
        !    CALL mexprintf(string) 
         ELSE IF ( INFI(I) .EQ. 1 ) THEN 
            WRITE(string,113) INDEX1(I),A(I)
 113        FORMAT (i2,' ', 5F5.2 ,'   inf ')
         !   CALL mexprintf(string) 
          !  CALL mexprintf('   inf ')
         ELSE 
            WRITE(string,114) INDEX1(I),A(I),B(I)
 114        FORMAT (i2,' ', 5F5.2 ,' ', 5F5.2)
           ! CALL mexprintf(string) 
         END IF
         
         DO K = 1,J
            WRITE(string,115) COF(I,K)
 115        FORMAT (' ',10F10.5) 
            !CALL mexprintf(string) 
         ENDDO
         !CALL mexprintf(CHAR(10)) 
      END DO
      RETURN
      END SUBROUTINE PRINTCOF
      SUBROUTINE C1C2(I0,I1,A,B,INFIN, Cm, B1, SQ, AJ, BJ, FINA,FINB)  
! The regression equation for the conditional distr. of Y given X=x
! is equal  to the conditional expectation of Y given X=x, i.e.,
! 
!       E(Y|X=x) = E(Y) + Cov(Y,X)/Var(X)[x-E(X)]
!
!  Let x1=(x-E(X))/SQRT(Var(X)) be zero mean, C1< x1 <C2, B1(I)=COV(Y(I),X)/SQRT(Var(X)). 
!  Then the process  Y(I) with mean Cm(I) can be written as 
!
!       y(I) = Cm(I) + x1*B1(I) + Delta(I) for  I=1,...,N.
!
!  where SQ(I) = sqrt(Var(Y|X)) is the standard deviation of Delta(I). 
!
!  Since we are truncating all Gaussian  variables to                   
!  the interval [-C,C], then if for any I                               
!                                                                       
!  a) Cm(I)+x1*B1(I)-C*SQ(I)>B(I)  or                             
!                                                                       
!  b) Cm(I)+x1*B1(I)+C*SQ(I)<A(I)  then                           
!                                                                       
!  the (XIND|Xn=xn) = 0 !!!!!!!!!                                               
!                                                                       
!  Consequently, for increasing the accuracy (by excluding possible 
!  discontinuouities) we shall exclude such values for which (XIND|X1=x1) = 0.
!  Hence we assume that if Aj<x<Bj any of the previous conditions are 
!  satisfied
!                                                                       
!  OBSERVE!!, Aj, Bj has to be set to (the normalized) lower and upper bounds 
!  of possible values for x1,respectively, i.e.,
!           Aj=max((A-E(X))/SQRT(Var(X)),-C), Bj=min((B-E(X))/SQRT(Var(X)),C) 
!  before calling C1C2 subroutine.                         
!
      USE GLOBALCONST
      USE GLOBALDATA, ONLY : EPS2,EPS ,xCutOff
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: Cm, B1, SQ, A,B
      INTEGER,          DIMENSION(:), INTENT(in) :: INFIN 
      DOUBLE PRECISION,            INTENT(inout) :: AJ,BJ
      INTEGER,                     INTENT(inout) :: FINA, FINB
      INTEGER, INTENT(IN) :: I0, I1
     
!     Local variables 
      DOUBLE PRECISION   :: xCut
      DOUBLE PRECISION, PARAMETER :: TOL = 1.0D-16     
      DOUBLE PRECISION :: AI,BI,CSQ,BdSQ0, LTOL 
      INTEGER :: INFI, I

      xCut = xCutOff
      LTOL = TOL ! EPS2
!      IF (AJ.GE.BJ) GO TO 112
!      CALL PRINTVAR(AJ,TXT='BC1C2: AJ')
!      CALL PRINTVAR(BJ,TXT='BC1C2: BJ')
     
      IF (I1.LT.I0)  RETURN       !Not able to change integration limits
      DO I = I0,I1             
!     C = xCutOff
         INFI = INFIN(I)
         IF (INFI.GE.0) THEN
            BdSQ0 = B1 (I)   
            CSQ   = xCut*SQ(I)
            IF (BdSQ0 .GT. LTOL) THEN
               IF ( INFI .NE. 0 ) THEN
                  IF (FINA.EQ.1) THEN
                     AJ = MAX(AJ,(A(I) - Cm(I) - CSQ)/BdSQ0)
                  ELSE
                     AJ = (A(I) - Cm(I) - CSQ)/BdSQ0
                     FINA = 1
                  ENDIF
                  IF (FINB.GT.0) AJ = MIN(AJ,BJ)
               END IF
               IF ( INFI .NE. 1 ) THEN
                  IF (FINB.EQ.1) THEN
                     BJ = MIN(BJ,(B(I) - Cm(I) + CSQ)/BdSQ0)
                  ELSE
                     BJ = (B(I) - Cm(I) + CSQ)/BdSQ0
                     FINB = 1
                  ENDIF 
                  IF (FINA.GT.0) BJ = MAX(AJ,BJ)             
               END IF
             ELSEIF (BdSQ0 .LT. -LTOL) THEN
                IF ( INFI .NE. 0 ) THEN
                  IF (FINB.EQ.1) THEN
                     BJ = MIN(BJ,(A(I) - Cm(I) - CSQ)/BdSQ0)
                  ELSE
                     BJ = (A(I) - Cm(I) - CSQ)/BdSQ0
                     FINB = 1
                  ENDIF
                  IF (FINA.GT.0) BJ = MAX(AJ,BJ)    
               END IF
               IF ( INFI .NE. 1 ) THEN
                  IF (FINA.EQ.1) THEN
                     AJ = MAX(AJ,(B(I) - Cm(I) + CSQ)/BdSQ0)
                  ELSE
                     AJ = (B(I) - Cm(I) + CSQ)/BdSQ0
                     FINA = 1
                  ENDIF
                  IF (FINB.GT.0) AJ = MIN(AJ,BJ)
               END IF
            END IF               
         ENDIF
      END DO
!      CALL PRINTVAR(AJ,TXT='AC1C2: AJ')
!      CALL PRINTVAR(BJ,TXT='AC1C2: BJ')
      RETURN
      END SUBROUTINE C1C2
      SUBROUTINE CVSRTXC (Nt,Nd,R,index1,INFORM)
      USE GLOBALDATA, ONLY :  XCEPS2
      USE GLOBALCONST
      IMPLICIT NONE
      INTEGER, INTENT(in) :: Nt,Nd 
      DOUBLE PRECISION, DIMENSION(:,:), INTENT(inout) :: R
      INTEGER,          DIMENSION(:  ), INTENT(inout) :: index1 
      INTEGER, INTENT(out) :: INFORM 
! local variables
      DOUBLE PRECISION, DIMENSION(:  ), allocatable   :: SQ
      INTEGER,          DIMENSION(1  )                :: m
      INTEGER :: M1,K,I,J,Ntdc,Ntd,Nc, Nullity,LO
      DOUBLE PRECISION, PARAMETER :: LTOL = 1.0D-16
!CVSRTXC calculate the conditional covariance matrix of Xt and Xd given Xc
! as well as the cholesky factorization for the Xc variable(s)
! The Xc variables are sorted by the largest conditional covariance
! 
! R         = In : Cov(X) where X=[Xt Xd Xc] is stochastic vector
!             Out: sorted Conditional Covar. matrix, i.e.,
!                    [ Cov([Xt,Xd] | Xc)    Shape N X N  (N=Ntdc=Nt+Nd+Nc) 
! index1    = In/Out : permutation vector giving the indices to the variables 
!             original place.   Size  Ntdc
! INFORM    = Out, Returns
!             0 If Normal termination.
!             1 If R is degenerate, i.e., Cov(Xc) is singular.
! 
! R=Cov([Xt,Xd,Xc]) is a covariance matrix of the stochastic 
! vector X=[Xt Xd Xc] where the variables Xt, Xd and Xc have the size
! Nt, Nd and Nc, respectively. 
! Xc are the conditional variables.
! Xd and Xt are the variables to integrate. 
!(Xd,Xt = variables in the jacobian and indicator respectively)
!
! Note: CVSRTXC only works on the upper triangular part of R
      INFORM = 0
      Ntdc = size(R,DIM=1)
      Ntd  = Nt   + Nd
      Nc   = Ntdc - Ntd
      
      IF (Nc.LT.1) RETURN

      ALLOCATE(SQ(1:Ntdc))
      DO I = 1, Ntdc
         SQ(I)  = R(I,I)
      ENDDO
	
	
      LO = 1
      K = Ntdc
      DO I = 1, Nc             ! Condsort Xc
         m  = K+1-MAXLOC(SQ(K:Ntd+1:-1)) 
         M1 = m(1)
         IF (SQ(m1).LE.XCEPS2) THEN
!     PRINT *,'CVSRTXC: Degenerate case of Xc(Nc-J+1) for J=',ix 
            !CALL mexprintf('CVSRTXC: Degenerate case of Xc(Nc-J+1)')
            INFORM = 1
            GOTO 200   ! RETURN    !degenerate case
         ENDIF
         IF (M1.NE.K) THEN
            ! Symmetric row column permuations
            ! Swap row and columns, but only upper triangular part
            CALL RCSWAP( M1, K, Ntdc,Ntd, R,INDEX1,SQ)
         END IF
         R(K,K) = SQRT(SQ(K))
         IF (K .EQ. LO) GOTO 200
         R(LO:K-1,K) = R(LO:K-1,K)/R(K,K)
! Cov(Xi,Xj|Xk,Xk+1,..,Xn)  = ....
!         Cov(Xi,Xj|Xk+1,..,Xn) - Cov(Xi,Xk|Xk+1,..Xn)*Cov(Xj,Xk|Xk+1,..Xn)
         DO J = LO,K-1
                                ! Var(Xj | Xk,Xk+1,...,Xn)
            SQ(J)  =  R(J,J) - R(J,K)*R(J,K)
            IF (SQ(J).LE.LTOL.AND.J.LE.Ntd) THEN
               IF (LO.LT.J) THEN
                  CALL RCSWAP(LO, J, Ntdc,Ntd, R,INDEX1,SQ)
               ENDIF
               R(LO,LO:K-1) = gZERO
               IF (SQ(LO).LT.-SQRT(LTOL)) THEN
                  !R(LO,K) = gZERO
                  !CALL mexprintf('Negative definit BIG!')
               ENDIF
               SQ(LO) = gZERO
               LO = LO + 1
            ELSE
               R(J,J) = SQ(J)
               R(LO:J-1,J) = R(LO:J-1,J) - R(LO:J-1,K)*R(J,K)
            ENDIF
         END DO 
         K = K - 1
      ENDDO
 200  DEALLOCATE(SQ)
      RETURN                                                             
      END SUBROUTINE CVSRTXC

      SUBROUTINE  RCSCALE(K,K0,N1,N,K1,CDI,Cm,R,A,B,INFI,INDEX1,Y)
      USE GLOBALCONST
      USE SWAPMOD
      IMPLICIT NONE
! RCSCALE:  Scale  covariance matrix rows and limits
!  If the conditional covariance matrix diagonal entry is zero, 
!  permute limits and/or rows, if necessary.
      INTEGER, INTENT(IN) :: K, K0
      INTEGER, INTENT(IN) :: N1, N
      INTEGER, INTENT(INOUT) :: K1
      DOUBLE PRECISION, DIMENSION(:),  INTENT(INOUT) :: CDI,A,B,Cm
      DOUBLE PRECISION, DIMENSION(:,:),INTENT(INOUT) :: R
      INTEGER,          DIMENSION(:),  INTENT(INOUT) :: INFI,INDEX1
      DOUBLE PRECISION, DIMENSION(:),OPTIONAL,INTENT(INOUT) :: Y
!Local variables
      DOUBLE PRECISION, PARAMETER :: LTOL = 1.0D-16
      INTEGER :: KK,K00, KKold, I, J, Ntdc
      K00 = K0
      DO WHILE( (0.LT.K00).AND. (ABS(R(K00,K)).LE.LTOL) )
         R(K00,K) = gZERO
         K00      = K00 - 1
      ENDDO
      IF (K00.GT.0) THEN
         CDI(K) = R(K00,K)
         A(K)   = A(K)/CDI(K)
         B(K)   = B(K)/CDI(K)                  
      
         IF ((CDI(K) .LT. gZERO).AND. INFI(K).GE. 0) THEN
            CALL SWAP(A(K),B(K))
            IF (INFI(K).NE. 2) INFI(K) = 1-INFI(K)
         END IF
      
         DO J = 1, K00
            ! conditional covariances
            R(J,K) = R(J,K)/CDI(K) 
            ! conditional standard dev.s used in regression eq.
            R(K,J) = R(K,J)/ABS(CDI(K)) 
         END DO
         DO J = K00 + 1, K
            R(J,K) = gZERO
            R(K,J) = gZERO
         END DO
         
         KKold = K
         KK = K-1
         DO I = K0, K00+1, -1
            DO WHILE ((I.LE.KK) .AND. ABS(R(I,KK)).GT.LTOL)
               DO J = 1,I       !K0 
                  ! SWAP Covariance matrix
                  CALL SWAP(R(J,KK),R(J,KKold))
                  ! SWAP conditional standarddeviations
                  CALL SWAP(R(KK,J),R(KKold,J))
               END DO
               CALL SWAP(CDI(KK),CDI(KKold))
               CALL SWAP(Cm(KK),Cm(KKold))
               CALL SWAP(INDEX1(KK),INDEX1(KKold))
               CALL SWAP(A(KK),A(KKold))
               CALL SWAP(B(KK),B(KKold))
               CALL SWAP(INFI(KK),INFI(KKold))
               IF (PRESENT(Y)) THEN
                   CALL SWAP(Y(KK),Y(KKold))
               ENDIF
               Ntdc = SIZE(R,DIM=1)
               IF (N.LT.Ntdc) THEN
                  ! SWAP Xc entries, i.e, Cov(Xt,Xc) and Cov(Xd,Xc)
                  DO J = N+1, Ntdc
                     CALL SWAP( R(KK,J), R(KKold,J) )
                  END DO
               ENDIF
               KKold = KK
               KK    = KK - 1
            ENDDO
         END DO
         IF (KK.LT.K1) THEN
            K1 = K1 + 1
            !CALL mexPrintf('RCSCALE: updated K1'//CHAR(10))
         END IF
!         CALL PRINTVAR(K,TXT='K')
!         CALL PRINTVAR(KK,TXT='KK')
!         CALL PRINTVAR(K1,TXT='K1')
!         CALL PRINTVAR(K00,TXT='K00')
!         CALL PRINTVAR(K0,TXT='K0')
!         CALL PRINTCOF(N,A,B,INFI,R,INDEX1)
      ELSE 
!     Remove variable it is conditional independent of all other variables
         !CALL mexPrintf('RCSCALE ERROR*********************')
      ENDIF
      END SUBROUTINE RCSCALE
      
      SUBROUTINE COVSRT(BCVSRT, Nt,Nd,R,Cm,A,B,INFI,INDEX1, 
     &     INFIS,INFISD, NDIM, Y, CDI )
      USE FIMOD
      USE SWAPMOD
      USE GLOBALCONST
      USE GLOBALDATA, ONLY : EPS2,NIT,xCutOff
      IMPLICIT NONE
!COVSRT  sort integration limits and determine Cholesky factor.
!
!     Nt, Nd = size info about Xt and Xd variables.
!     R      = Covariance/Cholesky factored matrix for [Xt,Xd,Xc] (in)
!              On input: 
!               1a) the first upper triangular the Nt + Nd times Nt + Nd
!                   block contains COV([Xt,Xd]|Xc)
!                 (conditional covariance matrix for Xt and Xd given Xc)
!               2a) The upper triangular part of the Nt+Nd+Nc times Nc
!                   last block contains the cholesky matrix for Xc,
!                   i.e.,
 
!              On output: 
!               1b) part 2a) mentioned above is unchanged, only necessary
!                 permutations according to INDEX1 is done.
!               2b) part 1a) mentioned above is changed to a special
!                  form of cholesky matrix: (N = Nt+Nd-INFIS-INFISD)
!                  R(1,1) = 1
!                  R(1,2:N) = [COV(X1,X2)/STD(X1),....COV(X1,XN)/STD(X1)]
!                  R(2,2) = 1
!                  R(2,3:N) = [COV(X2,X3)/STD(X2|X1),....COV(X2,XN)/STD(X2|X1)]
!              
!              Note: Only upper triangular part is needed.
!     Cm     = Conditional mean given Xc
!     A,B    = lower and upper integration limits length Nt+Nd
!     INFIN  = INTEGER, array of integration limits flags:  length Nt+Nd   (in)
!             if INFIN(I) < 0, Ith limits are (-infinity, infinity);
!             if INFIN(I) = 0, Ith limits are (-infinity, B(I)];
!             if INFIN(I) = 1, Ith limits are [A(I), infinity);
!             if INFIN(I) = 2, Ith limits are [A(I), B(I)].
!    INDEX1  = permutation index vector, i.e., giving the indices to the
!              variables original place.
!    INFIS   = Number of redundant variables of Xt
!    INFISD  = Number of redundant variables of Xd
!    NDIM    = Number of relevant dimensions to integrate. This is the
!             same as the rank of the submatrix of Cov([Xt,Xd]) minus
!             the INFIS variables of Xt and INFISD variables of Xd.
!    Y       = working array
!    CDI     = Cholesky diagonal elements which contains either
!               CDI(J) = STD(Xj | X1,...,Xj-1,Xc) if Xj is stochastic given
!                        X1,...Xj, Xc
!              or
!               CDI(J) = COV(Xj,Xk | X1,..,Xk-1,Xc  )/STD(Xk | X1,..,Xk-1,Xc) 
!               if Xj is determinstically determined given X1,..,Xk,Xc
!               for some k<j.
!       
!     Subroutine to sort integration limits and determine Cholesky
!     factor.
!     
!     Note: COVSRT works only on the upper triangular part of R
      LOGICAL,                          INTENT(in)    :: BCVSRT
      INTEGER,                          INTENT(in)    :: Nt,Nd
      DOUBLE PRECISION, DIMENSION(:,:), INTENT(inout) :: R
      DOUBLE PRECISION, DIMENSION(:  ), INTENT(inout) :: Cm,A,B
      INTEGER,          DIMENSION(:  ), INTENT(inout) :: INFI
      INTEGER,          DIMENSION(:  ), INTENT(inout) :: INDEX1
      INTEGER,                          INTENT(out) :: INFIS,INFISD,NDIM
      DOUBLE PRECISION, DIMENSION(:  ), INTENT(out)   :: Y, CDI
!     Local variables
      INTEGER :: N,N1,I, J, K, L, M, JMIN,Ndleft
      INTEGER ::  K1, K0, Nullity,INFJ,FINA,FINB
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: AA1,BB1
      DOUBLE PRECISION :: SUMSQ, AJ, BJ, TMP,D, E,EPSL
      DOUBLE PRECISION :: AA, Ca, Pa, APJ, PRBJ
      DOUBLE PRECISION :: CVDIAG, AMIN, BMIN, PRBMIN, YL, YU
      DOUBLE PRECISION :: SQTWPI,LTOL,TOL,ONE,ZERO,HALF 
      PARAMETER ( SQTWPI = 2.506628274631001D0, TOL = 1D-16 )
      PARAMETER (ONE = 1.D0, ZERO = 0.D0, HALF = 0.5D0) 
      
      EPSL   = EPS2
      INFIS  = 0
      INFISD = 0
      Ndim   = 0
      Ndleft = Nd
      N      = Nt+Nd
!      ALLOCATE(AA1(N))
!      ALLOCATE(BB1(N))
!      AA1(:) = ZERO
!      BB1(:) = ZERO


      LTOL = TOL
      TMP  = ZERO
      DO I = 1, N
         IF (R(I,I).GT.TMP) TMP = R(I,I)
      ENDDO
      EPSL = MAX(EPS2,1.0D-15)
      !IF (N.LT.10) EPSL = MIN(1D-10,EPSL)
      !LTOL = EPSL
!      EPSL = MAX(EPS2,LTOL)
      IF (TMP.GT.EPSL) THEN
         DO I = 1, N
            IF ((INFI(I) .LT. 0).OR.(R(I,I).LE.LTOL)) THEN
               IF (INDEX1(I).LE.Nt)  THEN 
                  INFIS = INFIS+1
               ELSEIF (R(I,I).LE.LTOL) THEN
                  INFISD = INFISD+1
               ENDIF
            ENDIF
         END DO
      ELSE
         LTOL       = EPSL
         INFIS      = Nt
         INFISD     = Nd
         R(1:N,1:N) = ZERO
      ENDIF
      
      N1  = N-INFIS-INFISD
      CDI(N1+1:SIZE(CDI)) = gZERO
      !PRINT *,'COVSRT'
      !CALL PRINTCOF(N,A,B,INFI,R,INDEX1)

!     Move any redundant variables of Xd to innermost positions. 
      DO I = N, N-INFISD+1, -1
         IF ( (R(I,I) .GT. LTOL) .OR. (INDEX1(I).LE.Nt)) THEN 
            DO J = 1,I-1
               IF ( (R(J,J) .LE. LTOL) .AND. (INDEX1(J).GT.Nt)) THEN
                  CALL RCSWAP(J, I, N, N, R,INDEX1,Cm, A, B, INFI)
                  GO TO 10
               ENDIF
            END DO
         ENDIF
 10   END DO  
!     
!     Move any doubly infinite limits or any redundant of Xt to the next
!     innermost positions.
!     
      DO I = N-INFISD, N1+1, -1
         IF ( ((INFI(I) .GE. 0).AND.(R(I,I).GT. LTOL))
     &        .OR. (INDEX1(I).GT.Nt)) THEN 
            DO J = 1,I-1
               IF ( (INFI(J) .LT. 0 .OR. (R(J,J).LE. LTOL)) 
     &              .AND. (INDEX1(J).LE.Nt)) THEN
                  CALL RCSWAP( J, I, N,N, R,INDEX1,Cm, A, B, INFI)
                  GO TO 15
               ENDIF
            END DO
         ENDIF
 15   END DO

!      CALL mexprintf('Before sorting')
!      CALL PRINTCOF(N,A,B,INFI,R,INDEX1)
!      CALL PRINTVEC(CDI,'CDI')
!      CALL PRINTVEC(Cm,'Cm')
      
      IF ( N1 .LE. 0 ) GOTO 200
!    
!     Sort remaining limits and determine Cholesky factor.
!     
      Y(:) = gZERO
      K    = 1
      Ndleft  = Nd - INFISD
      Nullity = 0
      DO  WHILE (K .LE. N1) 
     
!     Determine the integration limits for variable with minimum
!     expected probability and interchange that variable with Kth.
     
         K0     = K - Nullity
         PRBMIN = gTWO
         JMIN   = K
         CVDIAG = ZERO
         IF ((Ndleft.GT.0) .OR. (NDIM.LT.Nd+NIT)) THEN
            DO J = K, N1
               IF ( R(J,J) .GT. EPSL) THEN
                  TMP = ZERO    ! =  conditional mean of Y(I) given Y(1:I-1)
                  DO I = 1, K0 - 1
                     TMP = TMP + R(I,J)*Y(I)
                  END DO
                  SUMSQ = SQRT( R(J,J))
                  
                  IF (INFI(J).GE.0) THEN
                                ! May have infinite int. limits if Nd>0
                     IF (INFI(J).NE.0) THEN
                        AJ = ( A(J) - TMP )/SUMSQ
                     ENDIF
                     IF (INFI(J).NE.1) THEN
                        BJ = ( B(J) - TMP )/SUMSQ
                     ENDIF
                  ENDIF
                  IF (INDEX1(J).GT.Nt) THEN
                     AA = (Cm(J)+TMP)/SUMSQ ! inflection point
                     CALL EXLMS(AA,AJ,BJ,INFI(J),D,E,Ca,Pa)
                     PRBJ = E - D
                  ELSE   
                                !CALL MVNLMS( AJ, BJ, INFI(J), D, E )
                     CALL MVNLIMITS(AJ,BJ,INFI(J),APJ,PRBJ)
                    
                  ENDIF
                                !IF ( EMIN + D .GE. E + DMIN ) THEN
                  IF ( PRBJ .LT. PRBMIN ) THEN
                     JMIN = J
                     AMIN = AJ
                     BMIN = BJ
                     PRBMIN = MAX(PRBJ,ZERO)
                     CVDIAG = SUMSQ
                  ENDIF
               ENDIF
            END DO 
         END IF
!     
!     Compute Ith column of Cholesky factor.
!     Compute expected value for Ith integration variable (without
!     considering the jacobian) and
!     scale Ith covariance matrix row and limits.
!     
 40      IF ( CVDIAG.GT.TOL) THEN
            IF (INDEX1(JMIN).GT.Nt) THEN
               Ndleft = Ndleft-1             
            ENDIF
            NDIM = NDIM + 1     !Number of relevant dimensions to integrate

            IF ( K.LT.JMIN ) THEN
               CALL RCSWAP( K, JMIN, N1,N, R,INDEX1,Cm, A, B, INFI)
               ! SWAP conditional standarddeviations
               DO J = 1,K0-1  !MIN(K0, K-1)
                  CALL SWAP(R(K,J),R(JMIN,J))
               END DO
            END IF
                            
            R(K0,K) = CVDIAG 
            CDI(K)  = CVDIAG     ! Store the diagonal element
            DO I = K0+1,K
               R(I,K) = ZERO;
               R(K,I) = ZERO
            END DO
            
            K1  = K	
            I  = K1+1
            DO WHILE (I .LE. N1)              
               TMP = ZERO
               DO J = 1, K0-1
                  !tmp = tmp + L(i,j).*L(k1,j)
                  TMP = TMP + R(J,I)*R(J,K1) 
               END DO
                  ! Cov(Xk,Xi|X1,X2,...Xk-1)/STD(Xk|X1,X2,...Xk-1)
               R(K0,I)  = (R(K1,I) - TMP)/CVDIAG  
                  ! Var(Xi|X1,X2,...Xk)
               R(I,I) = R(I,I) - R(K0,I) * R(K0,I)

               IF (R(I,I).GT.LTOL) THEN
                  R(I,K0) = SQRT(R(I,I)) ! STD(Xi|X1,X2,...Xk)
               ELSE   !!IF (R(I,I) .LE. LTOL) THEN !TOL
                  R(I,K0) = MAX(SQRT(MAX(R(I,I), gZERO)),LTOL)
                  IF (INDEX1(I).GT.Nt) THEN
                     Ndleft = Ndleft - 1
                  ENDIF
                  Nullity = Nullity + 1
                  K  = K + 1
                  IF (K .LT. I) THEN
                     CALL RCSWAP( K, I, N1,N,R,INDEX1,Cm, A, B, INFI)
                     ! SWAP conditional standarddeviations
                     DO J = 1, K0
                        CALL SWAP(R(K,J),R(I,J))
                     END DO
                  ENDIF
                  CALL  RCSCALE(K,K0,N1,N,K1,CDI,Cm,R,A,B,INFI,INDEX1)
               END IF	
               I = I + 1
! 75            CONTINUE
            END DO
            INFJ = INFI(K1)
            IF (K1 .EQ.1) THEN
               FINA = 0
               FINB = 0
               IF (INFJ.GE.0) THEN
                  IF  (INFJ.NE.0) FINA = 1
                  IF  (INFJ.NE.1) FINB = 1
               ENDIF
            
               CALL C1C2(K1+1,N1,A,B,INFI, Y, R(K0,:),R(:,K0),
     &              AMIN, BMIN, FINA,FINB) 
               INFJ = 2*FINA+FINB-1
               CALL MVNLIMITS(AMIN,BMIN,INFJ,APJ,PRBMIN) 
            ENDIF

            IF ( PRBMIN .GT. ZERO) THEN
               YL = ZERO
               YU = ZERO
               IF (INFJ.GE.0) THEN
                  IF (INFJ .NE. 0) YL =-EXP(-HALF*(AMIN**2))/SQTWPI
                  IF (INFJ .NE. 1) YU =-EXP(-HALF*(BMIN**2))/SQTWPI
               ENDIF
               Y(K0) = ( YU - YL )/PRBMIN
            ELSE
               SELECT CASE (INFJ)
               CASE (:-1) 
                  Y(K0) = ZERO
               CASE (0) 
                  Y(K0) = BMIN
               CASE (1)
                  Y(K0) = AMIN
               CASE (2:)
                  Y(K0) = ( AMIN + BMIN )*HALF
               END SELECT 
            END IF
            R(K0,K1) = R(K0,K1)/CVDIAG 
            DO J = 1, K0 - 1
               ! conditional covariances
               R(J,K1) = R(J,K1)/CVDIAG 
               ! conditional standard dev.s used in regression eq.
               R(K1,J) = R(K1,J)/CVDIAG 
            END DO
            
            A(K1) = A(K1)/CVDIAG
            B(K1) = B(K1)/CVDIAG

            K  = K  + 1
! 100        CONTINUE
         ELSE
            R(K:N1,K:N1) = gZERO
!            CALL PRINTCOF(N,A,B,INFI,R,INDEX1)
            DO I = K, N1
!  Scale  covariance matrix rows and limits
!  If the conditional covariance matrix diagonal entry is zero, 
!  permute limits and/or rows, if necessary.
               CALL RCSCALE(I,K0-1,N1,N,K1,CDI,Cm,R,A,B,INFI,INDEX1)
            END DO
            Nullity = N1 - K0 + 1
            GOTO 200  !RETURN	
         END IF
      END DO 
 200  CONTINUE
      IF (Ndim .GT. 0) THEN  ! N1<K
         ! K1 = index to the last stochastic varible to integrate
         IF (INDEX1(K1).LE.Nt) Ndim = Ndim-1
      ENDIF
!      CALL mexprintf('After sorting')
!      CALL PRINTCOF(N,A,B,INFI,R,INDEX1)
!      CALL PRINTVEC(CDI,'CDI')
!      CALL PRINTVEC(Y,'Y')
!      CALL PRINTVEC(AA1,'AA1')
!      CALL PRINTVEC(BB1,'BB1')
      CALL PRINTVAR(NDIM,TXT='NDIM')
!      DEALLOCATE(AA1)
!      DEALLOCATE(BB1)
      RETURN
      END SUBROUTINE COVSRT

      SUBROUTINE COVSRT1(BCVSRT, Nt,Nd,R,Cm,A,B,INFI,INDEX1, 
     &     INFIS,INFISD, NDIM, Y, CDI )
      USE FIMOD
      USE SWAPMOD
      USE GLOBALCONST
      USE GLOBALDATA, ONLY : EPS2,NIT,xCutOff
      IMPLICIT NONE
!COVSRT  sort integration limits and determine Cholesky factor.
!
!     Nt, Nd = size info about Xt and Xd variables.
!     R      = Covariance/Cholesky factored matrix for [Xt,Xd,Xc] (in)
!              On input: 
!               1a) the first upper triangular the Nt + Nd times Nt + Nd
!                   block contains COV([Xt,Xd]|Xc)
!                   (conditional covariance matrix for Xt and Xd given Xc)
!               2a) The upper triangular part of the Nt+Nd+Nc times Nc
!                   last block contains the cholesky matrix for Xc, i.e., 
!
!              On output: 
!               1b) part 2a) mentioned above is unchanged, only necessary
!                   permutations according to INDEX1 is done.
!               2b) part 1a) mentioned above is changed to a special
!                   form of cholesky matrix: (N = Nt+Nd-INFIS-INFISD)
!                  R(1,1) = 1
!                  R(1,2:N) = [COV(X1,X2)/STD(X1),....COV(X1,XN)/STD(X1)]
!                  R(2,2) = 1
!                  R(2,3:N) = [COV(X2,X3)/STD(X2|X1),....COV(X2,XN)/STD(X2|X1)]
!              
!              Note: Only upper triangular part is needed.
!     Cm     = Conditional mean given Xc
!     A,B    = lower and upper integration limits length Nt+Nd
!     INFIN  = INTEGER, array of integration limits flags:  length Nt+Nd   (in)
!             if INFIN(I) < 0, Ith limits are (-infinity, infinity);
!             if INFIN(I) = 0, Ith limits are (-infinity, B(I)];
!             if INFIN(I) = 1, Ith limits are [A(I), infinity);
!             if INFIN(I) = 2, Ith limits are [A(I), B(I)].
!    INDEX1  = permutation index vector
!    INFIS   = Number of redundant variables of Xt
!    INFISD  = Number of redundant variables of Xd
!    NDIM    = Number of relevant dimensions to integrate. This is the
!             same as the rank of the submatrix of Cov([Xt,Xd]) minus
!             the INFIS variables of Xt and INFISD variables of Xd.
!    Y       = working array
!    CDI     = Cholesky diagonal elements which contains either
!               CDI(J) = STD(Xj| X1,,,Xj-1,Xc) if Xj is stochastic given
!                X1,...Xj, Xc
!              or
!               CDI(J) = COV(Xj,Xk|X1,..,Xk-1,Xc  )/STD(Xk| X1,,,Xk-1,Xc) 
!               if Xj is determinstically determined given X1,..,Xk,Xc
!               for some k<j.
!       
!     Subroutine to sort integration limits and determine Cholesky
!     factor.
!     
!     Note: COVSRT1 works only on the upper triangular part of R
      LOGICAL,                          INTENT(in)    :: BCVSRT
      INTEGER,                          INTENT(in)    :: Nt,Nd
      DOUBLE PRECISION, DIMENSION(:,:), INTENT(inout) :: R
      DOUBLE PRECISION, DIMENSION(:  ), INTENT(inout) :: Cm,A,B
      INTEGER,          DIMENSION(:  ), INTENT(inout) :: INFI
      INTEGER,          DIMENSION(:  ), INTENT(inout) :: INDEX1
      INTEGER,                          INTENT(out) :: INFIS,INFISD,NDIM
      DOUBLE PRECISION, DIMENSION(:  ), INTENT(out)   :: Y, CDI
!     Local variables
      INTEGER :: N,N1,I, J, K, L, M, JMIN,Ndleft
      INTEGER ::  K1, K0, Nullity, INFJ, FINA, FINB
      DOUBLE PRECISION :: SUMSQ, AJ, BJ, TMP,D, E,EPSL
      DOUBLE PRECISION :: AA, Ca, Pa, APJ, PRBJ
      DOUBLE PRECISION :: CVDIAG, AMIN, BMIN, PRBMIN, YL, YU
      DOUBLE PRECISION :: SQTWPI,LTOL,TOL,ONE,ZERO,HALF 
      PARAMETER ( SQTWPI = 2.506628274631001D0, TOL = 1D-16 )
      PARAMETER (ONE = 1.D0, ZERO = 0.D0, HALF = 0.5D0) 
      
      EPSL   = EPS2
      INFIS  = 0
      INFISD = 0
      Ndim   = 0
      Ndleft = Nd
      N      = Nt+Nd
      
!     IF (N.LT.10) EPSL = MIN(1D-10,EPS2)
      LTOL = TOL
      TMP  = ZERO
      DO I = 1, N
         IF (R(I,I).GT.TMP) TMP = R(I,I)
      ENDDO
      !      EPSL = MAX(EPS2*TMP,LTOL)
      EPSL = MAX(EPS2,LTOL)
      IF (TMP.GT.EPSL) THEN
         DO I = 1, N
            IF ((INFI(I) .LT. 0).OR.(R(I,I).LE.LTOL)) THEN
               IF (INDEX1(I).LE.Nt)  THEN 
                  INFIS = INFIS+1
               ELSEIF (R(I,I).LE.LTOL) THEN
                  INFISD = INFISD+1
               ENDIF
            ENDIF
         END DO
      ELSE
         !CALL PRINTCOF(N,A,B,INFI,R,INDEX1)
         !CALL PRINTVEC(CDI)
         !CALL PRINTVEC(Cm)
         LTOL   = EPSL
         INFIS  = Nt
         INFISD = Nd
         R(1:N,1:N) = ZERO
      ENDIF

      
      CDI(N-INFIS-INFISD+1:SIZE(CDI)) = ZERO
                                !PRINT *,'COVSRT'
      !CALL PRINTCOF(N,A,B,INFI,R,INDEX1)
         !     Move any redundant variables of Xd to innermost positions. 
      DO I = N, N-INFISD+1, -1
         IF ( R(I,I) .GT. LTOL .OR. INDEX1(I).LE.Nt) THEN 
            DO J = 1,I-1
               IF ( R(J,J) .LE. LTOL .AND. INDEX1(J).GT.Nt) THEN
                  CALL RCSWAP( J, I, N,N, R,INDEX1,Cm, A, B, INFI)
                  GO TO 10
               ENDIF
            END DO
         ENDIF
 10   END DO  
!     
!     Move any doubly infinite limits or any redundant of Xt to the next
!     innermost positions.
!     
      DO I = N-INFISD, N-INFISD-INFIS+1, -1
         IF ( ((INFI(I) .GE. 0).AND. (R(I,I).GT. LTOL) )
     &        .OR. INDEX1(I).GT.Nt) THEN 
            DO J = 1,I-1
               IF ( (INFI(J) .LT. 0 .OR. (R(J,J).LE. LTOL)) 
     &              .AND. INDEX1(J).LE.Nt) THEN
                  CALL RCSWAP( J, I, N,N, R,INDEX1,Cm, A, B, INFI)
                  GO TO 15
               ENDIF
            END DO
         ENDIF
 15   END DO
      
      IF ( INFIS +INFISD .GE. N ) RETURN
!      CALL mexprintf('Before sorting')
!      CALL PRINTCOF(N,A,B,INFI,R,INDEX1)

     
!     Sort remaining limits and determine Cholesky factor.
     
      K  = 1
      N1  = N-INFIS-INFISD
      Ndleft = Nd - INFISD
      Nullity = 0
      Y(1:N1) = ZERO
      DO  WHILE (K .LE. N1) 
     
!     Determine the integration limits for variable with minimum
!     expected probability and interchange that variable with Kth.
     
         K0     = K-Nullity
         PRBMIN = 2.d0
         JMIN   = K
         CVDIAG = ZERO
         IF (Ndleft.GT.0 .OR. NDIM.LT.Nd+NIT) THEN
            DO J = K,N1
               IF ( R(J,J) .GT. EPSL) THEN
                  TMP = Y(J)
                  !TMP = ZERO    ! =  the conditional mean of Y(I) given Y(1:I-1)
                  !DO I = 1, K0-1
                  !   TMP = TMP + R(I,J)*Y(I)
                  !END DO
                  SUMSQ = SQRT( R(J,J))
                  
                  IF (INFI(J).LT.0) GO TO 30 ! May have infinite int. limits if Nd>0
                  IF (INFI(J).NE.0) THEN
                     AJ = ( A(J) - TMP )/SUMSQ
                  ENDIF
                  IF (INFI(J).NE.1) THEN
                     BJ = ( B(J) - TMP )/SUMSQ
                  ENDIF
 30               IF (INDEX1(J).GT.Nt) THEN
                     AA = (Cm(J)+TMP)/SUMSQ ! inflection point
                     CALL EXLMS(AA,AJ,BJ,INFI(J),D,E,Ca,Pa)
                     PRBJ = E-D
                  ELSE
                                !CALL MVNLMS( AJ, BJ, INFI(J), D, E )
                     CALL MVNLIMITS(AJ,BJ,INFI(J),APJ,PRBJ)
                  ENDIF
                                !IF ( EMIN + D .GE. E + DMIN ) THEN
                  IF ( PRBJ .LT. PRBMIN ) THEN
                     JMIN = J
                     AMIN = AJ
                     BMIN = BJ
                     PRBMIN = MAX(PRBJ,ZERO)
                     CVDIAG = SUMSQ
                  ENDIF
               ENDIF
            END DO 
         END IF
!     
!     Compute Ith column of Cholesky factor.
!     Compute expected value for Ith integration variable (without
!     considering the jacobian) and
!     scale Ith covariance matrix row and limits.
!     
 40      IF ( CVDIAG.GT.TOL) THEN
            IF (INDEX1(JMIN).GT.Nt) THEN
               Ndleft = Ndleft-1
            ELSE
               IF (.FALSE..AND.PRBMIN+SQRT(LTOL).GE.gONE) THEN !BCVSRT.EQ.
                  D = gZERO
                  DO I = 1, K0-1
                     D = D + ABS(R(I,JMIN))*xCutOff
                  END DO
                  D = D/CVDIAG
                  M = INFI(JMIN)
                  CALL ADJLIMITS(AMIN+D,BMIN-D,M)
                  IF (M.LT.0) THEN
                     !variable is redundnant
                   IF ( JMIN.LT.N1 ) THEN
                   CALL RCSWAP( JMIN, N1, N1,N, R,INDEX1,Cm, A, B, INFI)
                                ! SWAP conditional standarddeviations
                      DO I = 1,K0-1 
                         CALL SWAP(R(JMIN,I),R(N1,I))
                      END DO
                      CALL SWAP(Y(N1),Y(JMIN))
                   ENDIF
                   INFIS = INFIS+1
                   N1    = N1-1
                   GOTO 100 
                  END IF      
                ENDIF
            ENDIF
            NDIM = NDIM+1     !Number of relevant dimensions to integrate

            IF ( K.LT.JMIN ) THEN
                              
               CALL RCSWAP( K, JMIN, N1,N, R,INDEX1,Cm, A, B, INFI)
               ! SWAP conditional standarddeviations
               DO J=1,K0-1
                  CALL SWAP(R(K,J),R(JMIN,J))
               END DO
               CALL SWAP(Y(K),Y(JMIN))
            END IF
            !IF (INDEX1(K).GT.Nt) Ndleft = Ndleft-1
            
            R(K0,K:N1) = R(K0,K:N1)/CVDIAG
            R(K0,K) = CVDIAG
            CDI(K)  = CVDIAG     ! Store the diagonal element
            DO I = K0+1,K
               R(I,K) = ZERO
               R(K,I) = ZERO
            END DO
            
            K1  = K
            !IF (K .EQ. N1) GOTO 200
            
!  Cov(Xi,Xj|Xk,Xk+1,..,Xn)=Cov(Xi,Xj|Xk+1,..,Xn) - Cov(Xi,Xk|Xk+1,..Xn)*Cov(Xj,Xk|Xk+1,..Xn)
            DO I = K1+1,N1
                                ! Var(Xj | Xk,Xk+1,...,Xn)
               R(I,I)  =  R(I,I) - R(K0,I)*R(K0,I)
               IF (R(I,I).GT.LTOL) THEN
                  R(I,K0)     = SQRT(R(I,I)) ! STD(Xi|X1,X2,...Xk)
                  R(I,I+1:N1) = R(I,I+1:N1) - R(K0,I+1:N1)*R(K0,I)
               ELSE
                  R(I,K0) = MAX(SQRT(MAX(R(I,I), gZERO)),LTOL)
                  !R(I,K0) = SQRT(LTOL) 
                  Nullity = Nullity + 1
                  K  = K + 1
                  IF (K .LT. I) THEN
                     CALL RCSWAP( K, I, N1,N,R,INDEX1,Cm, A, B, INFI)
                     ! SWAP conditional standarddeviations
                     DO J=1,K0
                        CALL SWAP(R(K,J),R(I,J))
                     END DO
                     CALL SWAP(Y(K),Y(I))
                  ENDIF     
                  IF (INDEX1(K).GT.Nt) Ndleft = Ndleft-1
                  CALL  RCSCALE(K,K0,N1,N,K1,CDI,Cm,R,A,B,INFI,INDEX1,Y)
!                  CDI(K) = R(K0,K)
                  !A(K) = (A(K)-xCutOff*R(K,K0))/CDI(K)
                  !B(K) = (B(K)+xCutOff*R(K,K0))/CDI(K)
!                  A(K) = A(K)/CDI(K)
!                  B(K) = B(K)/CDI(K)
!                  IF ((CDI(K) .LT. ZERO).AND. INFI(K).GE. 0) THEN
!                     CALL SWAP(A(K),B(K))
!                     IF (INFI(K).NE. 2) INFI(K) = 1-INFI(K)
!                  END IF
!                  DO J = 1, K0
!                     R(J,K) = R(J,K)/CDI(K)      ! conditional covariances
!                     R(K,J) = R(K,J)/ABS(CDI(K)) ! conditional standard dev.s used in regression eq.
!                  END DO
!                  DO J = K0+1,K
!                     R(J,K) = ZERO
!                     R(K,J) = ZERO
!                  END DO
               END IF	
            END DO 
            INFJ = INFI(K1)
            IF (K1 .EQ.1) THEN
               FINA = 0
               FINB = 0
               IF (INFJ.GE.0) THEN
                  IF  (INFJ.NE.0) FINA = 1
                  IF  (INFJ.NE.1) FINB = 1
               ENDIF
            
               CALL C1C2(K1+1,N1,A,B,INFI, Y, R(K0,:),R(:,K0),
     &              AMIN, BMIN, FINA,FINB) 
               INFJ = 2*FINA+FINB-1
               CALL MVNLIMITS(AMIN,BMIN,INFJ,APJ,PRBMIN) 
            ENDIF
                      
            IF ( PRBMIN .GT. ZERO) THEN
               YL = ZERO
               YU = ZERO
               IF (INFJ.GE.0) THEN
                  IF (INFJ .NE. 0) YL =-EXP(-HALF*(AMIN**2))/SQTWPI
                  IF (INFJ .NE. 1) YU =-EXP(-HALF*(BMIN**2))/SQTWPI
               ENDIF
               Y(K0) = ( YU - YL )/PRBMIN
            ELSE
               SELECT CASE (INFJ)
               CASE (:-1) 
                  Y(K0) = ZERO
               CASE (0) 
                  Y(K0) = BMIN
               CASE (1)
                  Y(K0) = AMIN
               CASE (2:)
                  Y(K0) = ( AMIN + BMIN )*HALF
               END SELECT 
            END IF
                     ! conditional mean (expectation) 
                                ! E(Y(K+1:N)|Y(1),Y(2),...,Y(K)) 
            Y(K+1:N1) = Y(K+1:N1)+Y(K0)*R(K0,K+1:N1)
            R(K0,K1)  = R(K0,K1)/CVDIAG ! conditional covariances
            DO J = 1, K0 - 1
               R(J,K1) = R(J,K1)/CVDIAG ! conditional covariances
               R(K1,J) = R(K1,J)/CVDIAG ! conditional standard dev.s used in regression eq.
            END DO

            A(K1) = A(K1)/CVDIAG
            B(K1) = B(K1)/CVDIAG

            K  = K  + 1
 100        CONTINUE
         ELSE 
            R(K:N1,K:N1) = gZERO
!            CALL PRINTCOF(N,A,B,INFI,R,INDEX1)
            DO I = K, N1
!  Scale  covariance matrix rows and limits
!  If the conditional covariance matrix diagonal entry is zero, 
!  permute limits and/or rows, if necessary.
               CALL RCSCALE(I,K0-1,N1,N,K1,CDI,Cm,R,A,B,INFI,INDEX1)
            END DO
            Nullity = N1 - K0 + 1
            GOTO 200  !RETURN	

            DO I = K,N1
               CDI(I) = R(K0-1,I)
               A(I) = A(I)/CDI(I)
               B(I) = B(I)/CDI(I)
               IF ((CDI(I) .LT. ZERO).AND. INFI(I).GE. 0) THEN
                  CALL SWAP(A(I),B(I))
                  IF (INFI(I).NE. 2) INFI(I) = 1-INFI(I)
               END IF
               DO J = 1,K0-1
                  R(J,I) = R(J,I)/CDI(I)     ! conditional covariances
                  R(I,J) = R(I,J)/ABS(CDI(I))! conditional standard dev.s used in regression eq.
               END DO	
               DO J=K0,I
                  R(J,I) = ZERO ! set the covariance to the rest to
                                ! zero
                  R(I,J) = ZERO
               END DO
            END DO
            Nullity = N1-K0+1
            GOTO 200
                                !RETURN			
         END IF
      END DO 
      
 200  CONTINUE
      IF (Ndim .GT. 0) THEN     ! N1<K
         ! K1 = index to the last stochastic varible to integrate
         IF (INDEX1(K1).LE.Nt) Ndim = Ndim-1
      ENDIF
!      CALL mexprintf('After sorting')
!      CALL PRINTCOF(N,A,B,INFI,R,INDEX1)
!      CALL PRINTVEC(CDI)
!      CALL PRINTVAR(NDIM,TXT='NDIM')
      RETURN
      END SUBROUTINE COVSRT1
	
*
      SUBROUTINE RCSWAP( P, Q, N,Ntd, C,IND,Cm, A, B, INFIN )
      USE SWAPMOD
      IMPLICIT NONE
* RCSWAP  Swaps rows and columns P and Q in situ, with P <= Q.
*
*
*   CALL  RCSWAP( P, Q, N, Ntd, C,IND A, B, INFIN, Cm)
*
*    P, Q  = row/column number to swap P<=Q<=N
*    N     = length of A, B
*    Ntd   = length(Xt)+length(Xd)
*    C     = upper triangular cholesky factor.Cov([Xt,Xd,Xc]) size Ntdc X Ntdc
*    IND   = permutation index vector.
*    Cm    = conditional mean
*    A,B   = lower and upper integration limit, respectively.
*    INFIN = if INFIN(I) < 0, Ith limits are (-infinity, infinity);
*            if INFIN(I) = 0, Ith limits are (-infinity, B(I)];
*            if INFIN(I) = 1, Ith limits are [A(I), infinity);
*            if INFIN(I) = 2, Ith limits are [A(I), B(I)].
*
      DOUBLE PRECISION, DIMENSION(:,:),INTENT(inout) :: C
      INTEGER, DIMENSION(:),INTENT(inout) :: IND
      INTEGER, DIMENSION(:),          OPTIONAL,INTENT(inout) :: INFIN
      DOUBLE PRECISION, DIMENSION(:), OPTIONAL,INTENT(inout) :: A,B,Cm
      INTEGER,INTENT(in) :: P, Q, N, Ntd
! local variable
      INTEGER :: J, Ntdc
      IF (PRESENT(Cm))    CALL SWAP( Cm(P), Cm(Q) )
      IF (PRESENT(A))     CALL SWAP( A(P), A(Q) )
      IF (PRESENT(B))     CALL SWAP( B(P), B(Q) )
      IF (PRESENT(INFIN)) CALL SWAP(INFIN(P),INFIN(Q))
	
      CALL SWAP(IND(P),IND(Q))
      
      CALL SWAP( C(P,P), C(Q,Q) )
      DO J = 1, P-1
         CALL SWAP( C(J,P), C(J,Q) )
      END DO
      DO J = P+1, Q-1
         CALL SWAP( C(P,J), C(J,Q) )
      END DO
      DO J = Q+1, N
         CALL SWAP( C(P,J), C(Q,J) )
      END DO
      Ntdc = SIZE(C,DIM=1)
      IF (N.LT.Ntdc) THEN
         DO J = Ntd+1, Ntdc
            CALL SWAP( C(P,J), C(Q,J) )
         END DO
      ENDIF
      RETURN
      END SUBROUTINE RCSWAP
      END MODULE COVSORTMOD


      MODULE FUNCVARMOD
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: BIG ! ChOlesky Factor/Covariance matrix
      DOUBLE PRECISION, DIMENSION(:  ), ALLOCATABLE :: CDI ! Cholesky DIagonal elements
      DOUBLE PRECISION, DIMENSION(:  ), ALLOCATABLE :: Cm,CDIXd,CmXd
      DOUBLE PRECISION, DIMENSION(:  ), ALLOCATABLE :: Xd,Xc,Y
      DOUBLE PRECISION, DIMENSION(:  ), ALLOCATABLE :: A,B      ! Integration limits
      INTEGER, DIMENSION(:  ), ALLOCATABLE :: INFI,INDEX1      
      INTEGER,SAVE :: Nt,Nd		  ! Size information
      INTEGER,SAVE :: INFIS,INFISD ! Number of redundant variables of Xt and Xd respectively
      END MODULE FUNCVARMOD

      MODULE FUNCMOD0            ! FUNCTION module containing constants transfeered to mvnfun
      USE FUNCVARMOD
      IMPLICIT NONE     
	! This module use FIINV AND regression equation C1C2    
      
      
! variables set in initfun and used in mvnfun:
      INTEGER, PRIVATE :: I0,NdleftN0
      DOUBLE PRECISION, PRIVATE :: E1,D1, VAL0
      
      INTERFACE  INITFUN
      MODULE PROCEDURE INITFUN
      END INTERFACE

      INTERFACE  MVNFUN
      MODULE PROCEDURE MVNFUN
      END INTERFACE

      CONTAINS

      SUBROUTINE INITFUN(VALUE,INFORM)
      USE JACOBMOD
      USE GLOBALCONST
      USE FIMOD
      USE COVSORTMOD, ONLY: C1C2
      USE GLOBALDATA, ONLY: EPS2,EPS,xCutOff
      IMPLICIT NONE
!      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(OUT) :: VALUE
      INTEGER, INTENT(out) :: INFORM
! local variables:
      INTEGER ::  N,NdleftO
      INTEGER :: I, J,  FINA, FINB, INFII
      DOUBLE PRECISION :: AI, BI
      DOUBLE PRECISION :: xCut 
      xCut = 0.0D0 !xCutOff
*     
*     Integrand subroutine
*
!INITFUN initilize the Multivariate Normal integrand function
! COF  - conditional sorted ChOlesky Factor of the covariance matrix (IN)
! CDI  - Cholesky DIagonal elements used to calculate the mean  
! Cm   - conditional mean of Xd and Xt given Xc, E(Xd,Xt|Xc)
! xd   - variables to the jacobian variable, need no initialization size Nd
! xc   - conditional variables (IN)
! INDEX1 - if INDEX1(I)>Nt then variable no. I is one of the Xd
!          variables otherwise it is one of Xt.

      !PRINT *,'Mvnfun,ndim',Ndim
      INFORM = 0
      VALUE  = 0.D0
      VAL0   = 1.D0

      NdleftN0 = Nd               ! Counter for number of Xd variables left
      
      I0   = 0
      FINA = 0
      FINB = 0
      N = Nt+Nd-INFIS-INFISD-1
      IF (INFIS+INFISD .GT. 0) THEN
!     CHCKLIM Check if the conditional mean Cm = E(Xt,Xd|Xc) for the
!     deterministic variables are between the barriers, i.e., A=Hlo-Cm< 0 <B=Hup-Cm
!     INFIN  INTEGER, array of integration limits flags:
!            if INFIN(I) < 0, Ith limits are (-infinity, infinity);
!            if INFIN(I) = 0, Ith limits are (-infinity, B(I)];
!            if INFIN(I) = 1, Ith limits are [A(I), infinity);
!            if INFIN(I) = 2, Ith limits are [A(I), B(I)].
		
         I = N+1
         DO J=1,INFIS+INFISD
            I = I+1
            IF (INFI(I).GE.0) THEN
               IF ((INFI(I).NE.0).AND.(EPS .LT.A(I))) GOTO 200
               IF ((INFI(I).NE.1).AND.(-EPS .GT.B(I))) GOTO 200
            ENDIF
         ENDDO
           
		
         IF (INFISD.GT.0) THEN       
                                ! Redundant variables of Xd: replace Xd with the mean
            I = Nt+Nd !-INFIS
            J = NdleftN0-INFISD
              
            DO WHILE (NdleftN0.GT.J)
               IF (INDEX1(I).GT.Nt) THEN ! isXd
                  xd (NdleftN0) =  Cm (I)
                  NdleftN0=NdleftN0-1                  
               END IF
               I = I-1
            ENDDO
         ENDIF
		
         IF (N+1.LT.1) THEN   
!     Degenerate case, No relevant variables left to integrate  
!     Print *,'rindd ndim1',Ndim1
            IF (Nd.GT.0) THEN 
               VALUE = jacob (xd,xc) ! jacobian of xd,xc
            ELSE
               VALUE = gONE
            END IF
            GOTO 200
         ENDIF
      ENDIF
      
      NdleftO = NdleftN0

      DO I = 1, N+1
         IF (INFI(I).LT.0) GO TO 100            ! May have infinite int. Limits if Nd>0
         IF ( INFI(I) .NE. 0 ) THEN
            IF ( FINA .EQ. 1 ) THEN
               AI = MAX( AI, A(I) - xCut*BIG(I,1))
            ELSE
               AI = A(I)        
               FINA = 1
            END IF
         END IF
         IF ( INFI(I) .NE. 1 ) THEN
            IF ( FINB .EQ. 1 ) THEN
               BI = MIN( BI, B(I) + xCut*BIG(I,1))
            ELSE
               BI = B(I)       
               FINB = 1
            END IF
         END IF
 100     IF (INDEX1(I).GT.Nt) THEN ! Save the mean for Xd
            CmXd(NdleftN0)  = Cm(I) 
            CDIXd(NdleftN0) = CDI(I)  
            NdleftN0 = NdleftN0-1
         END IF
    
         IF (I.EQ.N+1.OR.BIG(2,I+1).GT.gZERO) THEN
            Y = gZERO
            CALL C1C2(I+1,N+1,A,B,INFI, Y, BIG(1,:),BIG(:,1),
     &           AI, BI, FINA,FINB) 
            CALL MVNLMS( AI, BI,2*FINA+FINB-1, D1, E1 )
            IF (D1.GE.E1) GOTO 200
            
            IF (Nd.GT.0.AND. NdleftO.LE.0) VAL0 = JACOB(xd,xc)  
            IF (I.EQ.N+1.AND.NdleftO.LE.0) THEN
               VALUE = (E1-D1)*VAL0
               GO TO 200
            ELSEIF (N.EQ.1.AND.NdleftO.LE.0) THEN
               IF ( ABS( BIG(2,I+1) ) .GT. gZERO ) THEN
                  D1 = SQRT( 1 + BIG(1,I+1)**2 )
                  IF ( INFI(2) .NE. 0 ) A(2) = A(2)/D1
                  IF ( INFI(2) .NE. 1 ) B(2) = B(2)/D1
                  VALUE = BVNMVN( A, B,INFI,BIG(1,I+1)/D1 )*VAL0
               ELSE             ! correlation=+/-1
                  VALUE = (E1-D1)*VAL0
               END IF
               GO TO 200 
            ENDIF
            VAL0 = VAL0*(E1-D1)
            I0   = I
            RETURN
         ENDIF
      ENDDO
      RETURN
 200  INFORM = 1
      
      RETURN
      END SUBROUTINE INITFUN     
*     
*     Integrand subroutine
*
      FUNCTION MVNFUN( Ndim, W ) RESULT (VAL)
      USE JACOBMOD
      USE FIMOD
      USE GLOBALDATA, ONLY: EPS2,EPS,xCutOff
      IMPLICIT NONE
      INTEGER, INTENT (in) :: Ndim
      DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: W
      DOUBLE PRECISION :: VAL
! local variables:
      INTEGER ::  N,I, J, FINA, FINB
      INTEGER ::  NdleftN,NdleftO ,JMX,IK
      DOUBLE PRECISION :: TMP, AI, BI, DI, EI,xCut
!      DOUBLE PRECISION, PARAMETER :: EPSL = 1D-12
      
!MVNFUN Multivariate Normal integrand function
! where the integrand is transformed from an integral 
! having integration limits A and B to an
! integral having  constant integration limits i.e.
!   B                                  1 
! int jacob(xd,xc)*f(xd,xt)dxt dxd = int F2(W) dW
!  A                                  0
!
! W    - new transformed integration variables, valid range 0..1
!        The vector must have the length Ndim returned from Covsrt
! COF  - conditional sorted ChOlesky Factor of the covariance matrix (IN)
! CDI  - Cholesky DIagonal elements used to calculate the mean  
! Cm   - conditional mean of Xd and Xt given Xc, E(Xd,Xt|Xc)
! xd   - variables to the jacobian variable, need no initialization size Nd
! xc   - conditional variables (IN)
! INDEX1 - if INDEX1(I)>Nt then variable No. I is one of the Xd variables otherwise it is one of Xt
      !PRINT *,'Mvnfun,ndim',Ndim
      
      xCut = 0.0D0 ! xCutOff

      N = Nt+Nd-INFIS-INFISD-1
      IK = 1                    ! Counter for Ndim 
      FINA = 0
      FINB = 0

      NdleftN = NdleftN0          ! Counter for number of Xd variables left
      VAL  = VAL0
      NdleftO = Nd-INFISD
      Y(IK) = FIINV( D1 + W(IK)*( E1 - D1 ) )
      IF (NdleftO.GT.NdleftN ) THEN
         xd(NdleftN+1:NdleftO) = CmXd(NdleftN+1:NdleftO)+
     &        Y(IK)*CDIXd(NdleftN+1:NdleftO)
      ENDIF
      NdleftO = NdleftN
      IK = 2                    !=IK+1
            
     
      DO I = I0+1, N+1
         TMP = 0.d0
!         JMX = MAX(I-IK,0)+1
         DO J = 1, IK-1 !I-JMX
            ! E(Y(IK) | Y(1),...Y(IK-1))/STD(Y(IK)|Y(1),,,,Y(IK-1))
            TMP = TMP + BIG(J,I)*Y(J)  
         END DO
         IF (INFI(I).LT.0) GO TO 100            ! May have infinite int. Limits if Nd>0
         IF ( INFI(I) .NE. 0 ) THEN
            IF ( FINA .EQ. 1 ) THEN
               AI = MAX( AI, A(I) - TMP - xCut*BIG(I,IK)) 
            ELSE
               AI = A(I)-TMP 
               FINA = 1
            END IF
            IF (FINB.EQ.1.AND.BI.LE.AI) GOTO 200
         END IF
         IF ( INFI(I) .NE. 1 ) THEN
            IF ( FINB .EQ. 1 ) THEN
               BI = MIN( BI, B(I) - TMP + xCut*BIG(I,IK))	
            ELSE
               BI = B(I)-TMP 
               FINB = 1
            END IF
            IF (FINA.EQ.1.AND.BI.LE.AI) GOTO 200
         END IF
 100     IF (INDEX1(I).GT.Nt) THEN ! Save the mean of xd and Covariance diagonal element 
            CmXd(NdleftN) =  Cm(I)+TMP*CDI(I)    ! Conditional mean E(Xi|X1,..X)
            CDIXd(NdleftN) = CDI(I)              ! Covariance diagonal  
            NdleftN = NdleftN - 1
         END IF
         IF ( I .EQ. N+1 .OR. BIG(IK+1,I+1) .GT.0.D0 ) THEN 
! Must change so that Y = conditional expectation, in order to use C1C2 here:
!             CALL C1C2(I+1,N+1,A,B,INFI, Y, BIG(IK,:),BIG(:,IK),
!     &           AI, BI, FINA,FINB) 
            CALL MVNLMS( AI, BI, 2*FINA+FINB-1, DI, EI )            
            IF ( DI .GE. EI ) GO TO 200
            VAL = VAL*( EI - DI )
            
            IF ( I .LE. N .OR. NdleftN.LT.NdleftO) THEN
               Y(IK) = FIINV( DI + W(IK)*( EI - DI ) )
               IF (NdleftN.LT.NdleftO ) THEN
                  xd(NdleftN+1:NdleftO) = CmXd(NdleftN+1:NdleftO)+
     &                 Y(IK)*CDIXd(NdleftN+1:NdleftO)
                  NdleftO = NdleftN
               ENDIF
            ENDIF            
            IK   = IK + 1
            FINA = 0
            FINB = 0
         END IF
      END DO
      IF (Nd.GT.0) VAL = VAL*jacob(xd,xc)      
      RETURN
 200  VAL = 0.d0
      RETURN
      END FUNCTION MVNFUN
      END MODULE FUNCMOD0 
      MODULE FUNCMOD         ! FUNCTION module containing constants transfeered to mvnfun
      USE FUNCVARMOD
      IMPLICIT NONE     
	! This module use FIINV     
      
      
! variables set in initfun and used in mvnfun:
      INTEGER, PRIVATE :: I0,NdleftN0
      DOUBLE PRECISION, PRIVATE :: E1,D1, VAL0
      
      INTERFACE  INITFUN
      MODULE PROCEDURE INITFUN
      END INTERFACE

      INTERFACE  MVNFUN
      MODULE PROCEDURE MVNFUN
      END INTERFACE

      CONTAINS

      SUBROUTINE INITFUN(VALUE,INFORM)
      USE JACOBMOD
      USE GLOBALCONST
      USE FIMOD
      USE GLOBALDATA, ONLY: EPS2,EPS,xCutOff
      IMPLICIT NONE
!      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(OUT) :: VALUE
      INTEGER, INTENT(out) :: INFORM
! local variables:
      INTEGER ::  N,NdleftO
      INTEGER :: I, J,  FINA, FINB, INFII
      DOUBLE PRECISION :: AI, BI
*     
*     Integrand subroutine
*
!INITFUN initilize the Multivariate Normal integrand function
! COF  - conditional sorted ChOlesky Factor of the covariance matrix (IN)
! CDI  - Cholesky DIagonal elements used to calculate the mean  
! Cm   - conditional mean of Xd and Xt given Xc, E(Xd,Xt|Xc)
! xd   - variables to the jacobian variable, need no initialization size Nd
! xc   - conditional variables (IN)
! INDEX1 - if INDEX1(I)>Nt then variable no. I is one of the Xd variables otherwise it is one of Xt
      !PRINT *,'Mvnfun,ndim',Ndim
      INFORM = 0
      VALUE  = 0.D0
      VAL0   = 1.D0

      NdleftN0 = Nd               ! Counter for number of Xd variables left
      
      I0   = 0
      FINA = 0
      FINB = 0
      N = Nt+Nd-INFIS-INFISD-1
      IF (INFIS+INFISD .GT. 0) THEN
!     CHCKLIM Check if the conditional mean Cm = E(Xt,Xd|Xc) for the
!     deterministic variables are between the barriers, i.e., A=Hlo-Cm< 0 <B=Hup-Cm
!     INFIN  INTEGER, array of integration limits flags:
!            if INFIN(I) < 0, Ith limits are (-infinity, infinity);
!            if INFIN(I) = 0, Ith limits are (-infinity, B(I)];
!            if INFIN(I) = 1, Ith limits are [A(I), infinity);
!            if INFIN(I) = 2, Ith limits are [A(I), B(I)].
		
         I = N+1
         DO J=1,INFIS+INFISD
            I = I+1
            IF (INFI(I).GE.0) THEN
               IF ((INFI(I).NE.0).AND.(EPS .LT.A(I))) GOTO 200
               IF ((INFI(I).NE.1).AND.(-EPS .GT.B(I))) GOTO 200
            ENDIF
         ENDDO
           
		
         IF (INFISD.GT.0) THEN       
                                ! Redundant variables of Xd: replace Xd with the mean
            I = Nt+Nd !-INFIS
            J = NdleftN0-INFISD
              
            DO WHILE (NdleftN0.GT.J)
               IF (INDEX1(I).GT.Nt) THEN ! isXd
                  xd (NdleftN0) =  Cm (I)
                  NdleftN0      = NdleftN0-1                  
               END IF
               I = I-1
            ENDDO
         ENDIF
		
         IF (N+1.LT.1) THEN   
!     Degenerate case, No relevant variables left to integrate  
!     Print *,'rindd ndim1',Ndim1
            IF (Nd.GT.0) THEN 
               VALUE = jacob (xd,xc) ! jacobian of xd,xc
            ELSE
               VALUE = gONE
            END IF
            GOTO 200
         ENDIF
      ENDIF
      
      NdleftO = NdleftN0

      DO I = 1, N+1
         IF (INFI(I).LT.0) GO TO 100            ! May have infinite int. Limits if Nd>0
         IF ( INFI(I) .NE. 0 ) THEN
            IF ( FINA .EQ. 1 ) THEN
               AI = MAX( AI, A(I) - xCutOff*BIG(I,1))
            ELSE
               AI = A(I)        !MAX(A(I),-xCutOff) 
               FINA = 1
            END IF
         END IF
         IF ( INFI(I) .NE. 1 ) THEN
            IF ( FINB .EQ. 1 ) THEN
               BI = MIN( BI, B(I) + xCutOff*BIG(I,1))
            ELSE
               BI = B(I)        !MIN( B(I),xCutOff) 
               FINB = 1
            END IF
         END IF
 100     IF (INDEX1(I).GT.Nt) THEN ! Save the mean for Xd
            CmXd(NdleftN0) =  Cm(I) 
            CDIXd(NdleftN0) = CDI(I)  
            NdleftN0 = NdleftN0-1
         END IF
    
         IF (I.EQ.N+1.OR.BIG(2,I+1).GT.gZERO) THEN
            INFII = 2*FINA+FINB-1
            CALL MVNLMS( AI, BI,INFII , D1, E1 )
            IF (D1.GE.E1) GOTO 200
            
            IF (Nd.GT.0.AND. NdleftO.LE.0) VAL0 = JACOB(xd,xc)  
            IF (I.EQ.N+1.AND.NdleftO.LE.0) THEN
               VALUE = (E1-D1)*VAL0
               GO TO 200
            ELSEIF (N.EQ.1.AND.NdleftO.LE.0) THEN
               IF ( ABS( BIG(2,I+1) ) .GT. 0 ) THEN
                  D1 = SQRT( 1 + BIG(1,I+1)**2 )
                  IF ( INFI(2) .NE. 0 ) A(2) = A(2)/D1
                  IF ( INFI(2) .NE. 1 ) B(2) = B(2)/D1
                  VALUE = BVNMVN( A, B,INFI,BIG(1,I+1)/D1 )*VAL0
               ELSE             ! correlation=+/-1
                  VALUE = (E1-D1)*VAL0
               END IF
               GO TO 200 
            ENDIF
            VAL0 = VAL0*(E1-D1)
            I0   = I
            RETURN
         ENDIF
      ENDDO
      RETURN
 200  INFORM = 1
      
      RETURN
      END SUBROUTINE INITFUN
*     
*     Integrand subroutine
*
      FUNCTION MVNFUN( Ndim, W ) RESULT (VAL)
      USE JACOBMOD
      USE FIMOD
      USE GLOBALDATA, ONLY: EPS2,EPS,xCutOff
      IMPLICIT NONE
      INTEGER, INTENT (in) :: Ndim
      DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: W
      DOUBLE PRECISION :: VAL
! local variables:
      INTEGER ::  N,I, J, FINA, FINB
      INTEGER ::  NdleftN,NdleftO ,JMX,IK
      DOUBLE PRECISION :: TMP, AI, BI, DI, EI
      DOUBLE PRECISION, PARAMETER :: EPSL=1D-12
!MVNFUN Multivariate Normal integrand function
! where the integrand is transformed from an integral 
! having integration limits A and B to an
! integral having  constant integration limits i.e.
!   B                                  1 
! int jacob(xd,xc)*f(xd,xt)dxt dxd = int F2(W) dW
!  A                                  0
!
! W    - new transformed integration variables, valid range 0..1
!        The vector must have the length Ndim returned from Covsrt
! COF  - conditional sorted ChOlesky Factor of the covariance matrix (IN)
! CDI  - Cholesky DIagonal elements used to calculate the mean  
! Cm   - conditional mean of Xd and Xt given Xc, E(Xd,Xt|Xc)
! xd   - variables to the jacobian variable, need no initialization size Nd
! xc   - conditional variables (IN)
! INDEX1 - if INDEX1(I)>Nt then variable No. I is one of the Xd variables otherwise it is one of Xt
      !PRINT *,'Mvnfun,ndim',Ndim
      
      
      N = Nt+Nd-INFIS-INFISD-1
      IK = 1                    ! Counter for Ndim 
      FINA = 0
      FINB = 0

      NdleftN = NdleftN0          ! Counter for number of Xd variables left
      VAL  = VAL0
      NdleftO = Nd-INFISD
      Y(IK) = FIINV( D1 + W(IK)*( E1 - D1 ) )
      IF (NdleftO.GT.NdleftN ) THEN
         xd(NdleftN+1:NdleftO) = CmXd(NdleftN+1:NdleftO)+
     &        Y(IK)*CDIXd(NdleftN+1:NdleftO)
      ENDIF
      NdleftO = NdleftN
      IK = 2                    !=IK+1
            
     
      DO I = I0+1, N+1
         TMP = 0.d0
         JMX = MAX(I-IK,0)+1
         DO J = 1, I-JMX
            TMP = TMP + BIG(J,I)*Y(J)  ! E(Y(IK) | Y(1),...Y(IK-1))/STD(Y(IK)|Y(1),,,,Y(IK-1))
         END DO
         IF (INFI(I).LT.0) GO TO 100            ! May have infinite int. Limits if Nd>0
         IF ( INFI(I) .NE. 0 ) THEN
            IF ( FINA .EQ. 1 ) THEN
               AI = MAX( AI, A(I) - TMP) 
            ELSE
               AI = A(I)-TMP 
				!AI = MAX(A(I) - TMP,-xCutOff) 
               FINA = 1
            END IF
            IF (FINB.EQ.1.AND.BI.LE.AI) GOTO 200
         END IF
         IF ( INFI(I) .NE. 1 ) THEN
            IF ( FINB .EQ. 1 ) THEN
               BI = MIN( BI, B(I) - TMP)	
            ELSE
               BI = B(I)-TMP 
				!BI = MIN(B(I) - TMP,xCutOff)  
               FINB = 1
            END IF
            IF (FINA.EQ.1.AND.BI.LE.AI) GOTO 200
         END IF
 100     IF (INDEX1(I).GT.Nt) THEN ! Save the mean of xd and Covariance diagonal element 
            CmXd(NdleftN) =  Cm(I)+TMP*CDI(I)    ! Conditional mean E(Xi|X1,..X)
            CDIXd(NdleftN) = CDI(I)              ! Covariance diagonal  
            NdleftN = NdleftN-1
         END IF
         IF ( I .EQ. N+1 .OR. BIG(IK+1,I+1) .GT.0.D0 ) THEN   
            CALL MVNLMS( AI, BI, 2*FINA+FINB-1, DI, EI )            
            IF ( DI .GE. EI ) GO TO 200
            VAL = VAL*( EI - DI )
            
            IF ( I .LE. N .OR. NdleftN.LT.NdleftO) THEN
               Y(IK) = FIINV( DI + W(IK)*( EI - DI ) )
               IF (NdleftN.LT.NdleftO ) THEN
                  xd(NdleftN+1:NdleftO) = CmXd(NdleftN+1:NdleftO)+
     &                 Y(IK)*CDIXd(NdleftN+1:NdleftO)
                  NdleftO = NdleftN
               ENDIF
            ENDIF            
            IK   = IK + 1
            FINA = 0
            FINB = 0
         END IF
      END DO
      IF (Nd.GT.0) VAL=VAL*jacob(xd,xc)      
      RETURN
 200  VAL=0.d0
      RETURN
      END FUNCTION MVNFUN
      END MODULE FUNCMOD 

      MODULE FUNCMOD1        ! FUNCTION module containing constants transfeered to mvnfun
      USE FUNCVARMOD
      IMPLICIT NONE     

!     This module use FIINV and MVNLIMITS     
!     MVNLIMITS RETURN probabilities for being between A and B
!     WHERE 
!     AP = FI(A), AQ = 1 - FI(A)
!     BP = FI(B), BQ = 1 - FI(B)
!     PRB = BP-AP IF BP+AP<1
!         = AQ-BQ OTHERWISE
!
!      INTERFACE
!         SUBROUTINE MVNLIMITS( A, B, INFIN, AP, PRB, AQ)
!         DOUBLE PRECISION, INTENT(in) :: A, B
!         DOUBLE PRECISION, INTENT(out) :: AP
!         DOUBLE PRECISION, INTENT(out),OPTIONAL :: PRB,AQ
!         INTEGER,INTENT(in) :: INFIN
!         END SUBROUTINE MVNLIMITS 
!      END INTERFACE
!      INTERFACE
!         FUNCTION FIINV( Z ) RESULT (VALUE)
!         DOUBLE PRECISION, INTENT(in) :: Z
!         DOUBLE PRECISION :: VALUE
!         END FUNCTION FIINV
!      END INTERFACE
      
! Variables set in initfun and used in mvnfun
      INTEGER, PRIVATE,          SAVE 	:: I0,NdleftN0
      DOUBLE PRECISION, PRIVATE, SAVE	:: AP1,AQ1,PRB1,VAL0
      
      INTERFACE  INITFUN
      MODULE PROCEDURE INITFUN
      END INTERFACE

      INTERFACE  MVNFUN
      MODULE PROCEDURE MVNFUN
      END INTERFACE
      

      CONTAINS
     
      SUBROUTINE INITFUN(VALUE,INFORM)
      USE JACOBMOD
      USE FIMOD
      USE GLOBALDATA, ONLY: EPS2,EPS,xCutOff
      IMPLICIT NONE
!      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(OUT) :: VALUE
      INTEGER, INTENT(out) :: INFORM
! local variables:
      INTEGER ::  N,NdleftO
      INTEGER :: I, J,  FINA, FINB
      DOUBLE PRECISION :: AI, BI
	DOUBLE PRECISION,PARAMETER :: ZERO = 0.D0, ONE = 1.D0
*     
*     Integrand subroutine
*
!INITFUN initilize the Multivariate Normal integrand function
! COF  - conditional sorted ChOlesky Factor of the covariance matrix (IN)
! CDI  - Cholesky DIagonal elements used to calculate the mean  
! Cm   - conditional mean of Xd and Xt given Xc, E(Xd,Xt|Xc)
! xd   - variables to the jacobian variable, need no initialization size Nd
! xc   - conditional variables (IN)
! INDEX1 - if INDEX1(I)>Nt then variable no. I is one of the Xd variables otherwise it is one of Xt
      !PRINT *,'Mvnfun,ndim',Ndim
      INFORM = 0
      VALUE  = ZERO
      VAL0   = ONE

      NdleftN0 = Nd             ! Counter for number of Xd variables left
      
      N = Nt+Nd-INFIS-INFISD-1  
      
      IF (INFIS+INFISD .GT. 0) THEN
! CHCKLIM Check if the conditional mean Cm = E(Xt,Xd|Xc) for the deterministic variables
!         is between the barriers, i.e., A=Hlo-Cm< 0 <B=Hup-Cm
!     INFIN  INTEGER, array of integration limits flags:
!            if INFIN(I) < 0, Ith limits are (-infinity, infinity);
!            if INFIN(I) = 0, Ith limits are (-infinity, B(I)];
!            if INFIN(I) = 1, Ith limits are [A(I), infinity);
!            if INFIN(I) = 2, Ith limits are [A(I), B(I)].
		
         I = N+1
         DO J=1,INFIS+INFISD
            I = I+1
            IF (INFI(I).GE.0) THEN
               IF ((INFI(I).NE.0).AND.(EPS .LT.A(I))) GOTO 200
               IF ((INFI(I).NE.1).AND.(-EPS .GT.B(I))) GOTO 200
            ENDIF
         ENDDO
         
		
         IF (INFISD.GT.0) THEN       
                                ! Redundant variables of Xd: replace Xd with the mean
            I = Nt+Nd !-INFIS
            J = NdleftN0-INFISD
            
            DO WHILE (NdleftN0.GT.J)
               IF (INDEX1(I).GT.Nt) THEN ! isXd
                  xd (NdleftN0) =  Cm (I)
                  NdleftN0=NdleftN0-1                  
               END IF
               I = I-1
            ENDDO
         ENDIF
		
         IF (N+1.LT.1) THEN     !degenerate case, No relevant variables left to integrate  
!     Print *,'rindd ndim1',Ndim1
            IF (Nd.GT.0) THEN 
               VALUE = jacob (xd,xc) ! jacobian of xd,xc
            ELSE
               VALUE = ONE
            END IF
            GOTO 200
         ENDIF
      ENDIF

      NdleftO = NdleftN0   
      
      I0   = 0
      FINA = 0
      FINB = 0

      DO I = 1, N+1
         IF (INFI(I).LT.0) GO TO 100            ! May have infinite int. Limits if Nd>0
         IF ( INFI(I) .NE. 0 ) THEN
            IF ( FINA .EQ. 1 ) THEN
               AI = MAX( AI, A(I) )
            ELSE
               AI = A(I) !MAX(A(I),-xCutOff) 
               FINA = 1
            END IF
         END IF
         IF ( INFI(I) .NE. 1 ) THEN
            IF ( FINB .EQ. 1 ) THEN
               BI = MIN( BI, B(I) )
            ELSE
               BI = B(I) !MIN( B(I),xCutOff) 
               FINB = 1
            END IF
         END IF
 100     IF (INDEX1(I).GT.Nt) THEN ! Save the mean for Xd
            CmXd(NdleftN0) =  Cm(I) 
            CDIXd(NdleftN0) = CDI(I)  
            NdleftN0 = NdleftN0-1
         END IF
    
         IF (I.EQ.N+1.OR.BIG(2,I+1).GT.ZERO) THEN
            CALL MVNLIMITS( AI,BI,2*FINA+FINB-1,AP1, PRB1,AQ1)
            IF (PRB1.LE.ZERO) GOTO 200
            
            IF (Nd.GT.0.AND. NdleftO.LE.0) VAL0 = JACOB(xd,xc)  
            IF (I.EQ.N+1.AND.NdleftO.LE.0) THEN
               VALUE = VAL0*PRB1
               GO TO 200
            ELSEIF (N.EQ.1.AND.NdleftO.LE.0) THEN
               IF ( ABS( BIG(2,I+1) ) .GT. ZERO ) THEN
                  AP1 = SQRT( ONE + BIG(1,I+1)**2 )
					INFI(2) =INFI(I+1)
                  IF ( INFI(2) .NE. 0 ) A(2) = A(I+1)/AP1
                  IF ( INFI(2) .NE. 1 ) B(2) = B(I+1)/AP1
                  VALUE = BVNMVN( A, B,INFI,BIG(1,I+1)/AP1 )*VAL0
               ELSE             ! correlation=+/-1
                  VALUE = VAL0*PRB1
               END IF
               GO TO 200 
            ENDIF
            VAL0 = VAL0*PRB1
            I0   = I
            RETURN
         ENDIF
      ENDDO
      RETURN
 200  INFORM = 1
      
      RETURN
      END SUBROUTINE INITFUN

      FUNCTION MVNFUN( Ndim, W ) RESULT (VAL)
      USE JACOBMOD
      USE FIMOD
      USE GLOBALDATA, ONLY: EPS2,EPS,xCutOff
      IMPLICIT NONE
      INTEGER, INTENT (in) :: Ndim
      DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: W
      DOUBLE PRECISION :: VAL
! local variables:
      INTEGER ::  N,I, J, FINA, FINB
      INTEGER ::  NdleftN,NdleftO ,JMX,IK
      DOUBLE PRECISION :: TMP, AI, BI, API, AQI,PRBI
      DOUBLE PRECISION, PARAMETER :: EPSL=1D-12, HALF = 0.5D0
	DOUBLE PRECISION,PARAMETER :: ZERO = 0.D0, ONE = 1.D0
*     
*     Integrand subroutine
*
!MVNFUN Multivariate Normal integrand function
! where the integrand is transformed from an integral 
! having integration limits A and B to an
! integral having  constant integration limits i.e.
!   B                                  1 
! int jacob(xd,xc)*f(xd,xt)dxt dxd = int F2(W) dW
!  A                                  0
!
! W    - new transformed integration variables, valid range 0..1
!        The vector must have the length Ndim returned from Covsrt
! COF  - conditional sorted ChOlesky Factor of the covariance matrix (IN)
! CDI  - Cholesky DIagonal elements used to calculate the mean  
! Cm   - conditional mean of Xd and Xt given Xc, E(Xd,Xt|Xc)
! xd   - variables to the jacobian variable, need no initialization size Nd
! xc   - conditional variables (IN)
! INDEX1 - if INDEX1(I)>Nt then variable No. I is one of the Xd variables otherwise it is one of Xt
      !PRINT *,'Mvnfun,ndim',Ndim
      
      
      N = Nt+Nd-INFIS-INFISD-1
      IK = 1                    ! Counter for Ndim 
      FINA = 0
      FINB = 0

      NdleftN = NdleftN0          ! Counter for number of Xd variables left
      VAL  = VAL0
      NdleftO = Nd-INFISD
      IF (AP1.GT.HALF) THEN
         Y(IK) = -FIINV( AQ1 - W(IK)*PRB1 )
      ELSE
         Y(IK) =  FIINV( AP1 + W(IK)*PRB1 )
      ENDIF
      IF (NdleftO.GT.NdleftN ) THEN
         xd(NdleftN+1:NdleftO) = CmXd(NdleftN+1:NdleftO)+
     &        Y(IK)*CDIXd(NdleftN+1:NdleftO)
      ENDIF
      NdleftO = NdleftN
      IK = 2                    !=IK+1
            
     
      DO I = I0+1, N+1
         TMP = ZERO
         JMX = MAX(I-IK,0)+1
         DO J = 1, I-JMX
            TMP = TMP + BIG(J,I)*Y(J) ! E(Y(IK) | Y(1),...Y(IK-1))/STD(Y(IK)|Y(1),,,,Y(IK-1))
         END DO
         IF (INFI(I).LT.0) GO TO 100            ! May have infinite int. Limits if Nd>0
         IF ( INFI(I) .NE. 0 ) THEN
            IF ( FINA .EQ. 1 ) THEN
!     IF (ABS(CDI(I)).GT.EPS2) THEN
               AI = MAX( AI, A(I) - TMP) 
!     ELSE	
!     PRINT * ,'CDI(',I,') = ' , CDI(I)
!     PRINT * ,'AI = ', AI , ' AI2 =  ', A(I) - TMP
!     ENDIF
            ELSE
               AI = A(I)-TMP 
                                !AI = MAX(A(I) - TMP,-xCutOff) 
               FINA = 1
            END IF
            IF (FINB.EQ.1.AND.BI.LE.AI) GOTO 200
         END IF
         IF ( INFI(I) .NE. 1 ) THEN
            IF ( FINB .EQ. 1 ) THEN
!     IF (ABS(CDI(I)).GT.EPS2) THEN
               BI = MIN( BI, B(I) - TMP)	
!     ELSE
!     PRINT * ,'CDI(',I,') = ' , CDI(I)
!     PRINT * ,'BI = ', BI , ' BI2 =  ', B(I) - TMP
!     ENDIF
            ELSE
               BI = B(I)-TMP 
                                !BI = MIN(B(I) - TMP,xCutOff)  
               FINB = 1
            END IF
            IF (FINA.EQ.1.AND.BI.LE.AI) GOTO 200
         END IF
 100     IF (INDEX1(I).GT.Nt) THEN ! Save the mean of xd and Covariance diagonal element 
            CmXd(NdleftN) =  Cm(I)+TMP*CDI(I)    ! Conditional mean E(Xi|X1,..X)
            CDIXd(NdleftN) = CDI(I)              ! Covariance diagonal  
            NdleftN = NdleftN-1
         END IF
         IF ( I .EQ. N+1 .OR. BIG(IK+1,I+1) .GT.ZERO) THEN   
            CALL MVNLIMITS(AI,BI,2*FINA+FINB-1,API,PRBI,AQI)            
            IF ( PRBI .LE. ZERO ) GO TO 200
            VAL = VAL*PRBI
            
            IF ( I .LE. N .OR. NdleftN.LT.NdleftO) THEN
               IF (API.GT.HALF) THEN
                  Y(IK) = -FIINV( AQI - W(IK)*PRBI )
               ELSE
                  Y(IK) =  FIINV( API + W(IK)*PRBI )
               ENDIF
               IF (NdleftN.LT.NdleftO ) THEN
                  xd(NdleftN+1:NdleftO) = CmXd(NdleftN+1:NdleftO)+
     &                 Y(IK)*CDIXd(NdleftN+1:NdleftO)
                  NdleftO = NdleftN
               ENDIF
            ENDIF            
            IK   = IK + 1
            FINA = 0
            FINB = 0
         END IF
      END DO
      IF (Nd.GT.0) VAL=VAL*jacob(xd,xc)      
      RETURN
 200  VAL=0.d0
      RETURN
      END FUNCTION MVNFUN
      END MODULE FUNCMOD1 


      MODULE FUNCMOD2   ! FUNCTION module containing constants transfeered to mvnfun
	USE FUNCVARMOD
      IMPLICIT NONE     
	! This module use EXINV and FIINV
      
!      INTERFACE
!         SUBROUTINE MVNLMS( A, B, INFIN, LOWER, UPPER )
!         DOUBLE PRECISION, INTENT(in) :: A, B
!         DOUBLE PRECISION, INTENT(out) :: LOWER, UPPER
!         INTEGER,INTENT(in) :: INFIN
!         END SUBROUTINE MVNLMS 
!      END INTERFACE
!
!      INTERFACE
!         FUNCTION FIINV( Z ) RESULT (VALUE)
!         DOUBLE PRECISION, INTENT(in) :: Z
!         DOUBLE PRECISION :: VALUE
!         END FUNCTION FIINV
!      END INTERFACE
      
! Variables set in initfun and used in mvnfun
      INTEGER, PRIVATE, SAVE :: I0,NdleftN0
      DOUBLE PRECISION,PRIVATE, SAVE :: E1,D1
      DOUBLE PRECISION,PRIVATE, SAVE :: VAL0,AA0,Ca0,Pa0
      
      INTERFACE  INITFUN
      MODULE PROCEDURE INITFUN
      END INTERFACE

      INTERFACE  MVNFUN
      MODULE PROCEDURE MVNFUN
      END INTERFACE      

      CONTAINS

      SUBROUTINE INITFUN(VALUE,INFORM)
      USE FIMOD
      USE GLOBALDATA, ONLY: EPS2,EPS,xCutOff
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(out) :: VALUE
      INTEGER, INTENT(out) :: INFORM
! local variables:
      INTEGER ::  N,NdleftO
      INTEGER :: I, J,  FINA, FINB
      DOUBLE PRECISION :: AI, BI
*     
*     Integrand initialiazation subroutine
*
!INITFUN initialize the Multivariate Normal integrand function
!   B                                  1 
! int ABS(xd(1)*xd(2)...xd(Nd))*f(xd,xt)dxt dxd = int F2(W) dW
!  A                                  0
!
! COF  - conditional sorted ChOlesky Factor of the covariance matrix (IN)
! CDI  - Cholesky DIagonal elements used to calculate the mean  
! Cm   - conditional mean of Xd and Xt given Xc, E(Xd,Xt|Xc)
! xd   - variables to the jacobian variable, need no initialization size Nd
! xc   - conditional variables (IN)
! INDEX1 - if INDEX1(I)>Nt then variable no. I is one of the Xd variables otherwise it is one of Xt
      !PRINT *,'Mvnfun,ndim',Ndim
      INFORM = 0
      VALUE  = 0.D0
      VAL0   = 1.d0

      CDIXd  = 0.d0
      CmXd   = 0.d0
      AA0    = 0.d0
      Pa0    = 0.5D0
      Ca0    = 1.D0

      FINA = 0
      FINB = 0
      I0   = 0
      NdleftN0 = Nd               ! Counter for number of Xd variables left
      N    = Nt+Nd-INFIS-INFISD-1
      IF (INFIS+INFISD .GT. 0) THEN
! CHCKLIM Check if the conditional mean Cm = E(Xt,Xd|Xc) for the deterministic variables
!         is between the barriers, i.e., A=Hlo-Cm< 0 <B=Hup-Cm
!     INFIN  INTEGER, array of integration limits flags:
!            if INFIN(I) < 0, Ith limits are (-infinity, infinity);
!            if INFIN(I) = 0, Ith limits are (-infinity, B(I)];
!            if INFIN(I) = 1, Ith limits are [A(I), infinity);
!            if INFIN(I) = 2, Ith limits are [A(I), B(I)].
         I = N+1
         DO J=1,INFIS+INFISD
            I = I+1
            IF (INFI(I).GE.0) THEN
               IF ((INFI(I).NE.0).AND.(EPS .LT.A(I))) GOTO 200
               IF ((INFI(I).NE.1).AND.(-EPS .GT.B(I))) GOTO 200
            ENDIF
         ENDDO
		
         IF (INFISD.GT.0) THEN  ! replace xd with the mean for all deterministic variables
            I = Nt+Nd !-INFIS
            J = Nd-INFISD
                                !PRINT *,'Ndleft0',Ndleft0
            DO WHILE (NdleftN0.GT.J)
               IF (INDEX1(I).GT.Nt) THEN ! isXd
                  VAL0 = VAL0*ABS(Cm (I))
                  NdleftN0=NdleftN0-1                  
               END IF
               I = I-1
            ENDDO
         ENDIF
		
         IF (N+1.LT.1) THEN     !degenerate case, No relevant variables left to integrate  
!     Print *,'rindd ndim1',Ndim1
            VALUE = VAL0
            GOTO 200
         ENDIF
	ENDIF
      

      NdleftO = NdleftN0

      DO I = 1, N+1
         IF (INFI(I).LT.0) GO TO 100            ! May have infinite int. Limits if Nd>0
         IF ( INFI(I) .NE. 0 ) THEN
            IF ( FINA .EQ. 1 ) THEN
               AI = MAX( AI, A(I) ) 
            ELSE
                                !AI = MAX(A(I) ,-xCutOff) 
               AI = A(I)
               FINA = 1
            END IF
         END IF
         IF ( INFI(I) .NE. 1 ) THEN
            IF ( FINB .EQ. 1 ) THEN
               BI = MIN( BI, B(I) ) 
            ELSE
                                !BI = MIN(B(I),xCutOff) 
               BI = B(I)
               FINB = 1
            END IF
         END IF
 100     IF (I.GT.1.AND.INDEX1(I).GT.Nt) THEN ! Save the mean for Xd
            CmXd(NdleftN0)  = Cm(I)
            CDIXd(NdleftN0) = CDI(I)
            NdleftN0 = NdleftN0-1
         END IF
         IF ( I .EQ. N+1 .OR. BIG(2,I+1) .GT. 0.d0 ) THEN 
            
            IF ( INDEX1(1).GT.Nt) THEN
               AA0   = Cm(1)/CDI(1)
               CALL EXLMS(AA0,AI, BI, 2*FINA+FINB-1, D1, E1,Ca0,Pa0)
               VAL0 = VAL0*ABS(CDI(1)*Ca0)
            ELSE
               CALL MVNLMS( AI, BI, 2*FINA+FINB-1, D1, E1 )
            ENDIF
            IF (D1.GE.E1) GOTO 200
            
            IF (I.EQ.N+1.AND.NdleftO.LE.0) THEN
               VALUE = VAL0*(E1-D1)
               GO TO 200
            ELSEIF (N.EQ.1.AND.NdleftO.LE.0) THEN
               IF ( ABS( BIG(2,I+1)) .GT. 0 ) THEN
                  D1 = SQRT( 1 + BIG(I+1,1)**2 )
                  IF ( INFI(I+1) .NE. 0 ) A(2) = A(I+1)/D1
                  IF ( INFI(I+1) .NE. 1 ) B(2) = B(I+1)/D1
                  VALUE = BVNMVN( A, B, INFI,BIG(1,I+1)/D1 )*VAL0
               ELSE             ! correlation=+/-1
                  VALUE = VAL0*(E1-D1)
               END IF
               GO TO 200 
            ENDIF
            VAL0 = VAL0*(E1-D1)
            I0 = I
            !PRINT * ,'INITFUN VAL0,AA0,Ca0,Pa0',VAL0,AA0,Ca0,Pa0
            !PRINT *,'Cmxd',CmXd,'CDIXd',CDIXd
            RETURN
         ENDIF
      END DO
      RETURN
 200  INFORM = 1
      !PRINT * ,'INITFUN VALUE',VALUE
      RETURN
      END SUBROUTINE INITFUN

      FUNCTION MVNFUN( Ndim, W ) RESULT (VAL)
      USE FIMOD
      USE GLOBALDATA, ONLY: EPS2,EPS,xCutOff
      IMPLICIT NONE
      INTEGER, INTENT (in) :: Ndim
      DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: W
      DOUBLE PRECISION :: VAL
! local variables:
      INTEGER ::  N,NdleftN,NdleftO
      INTEGER :: I, J,  FINA, FINB,LK,JMX,IK
      DOUBLE PRECISION :: TMP, AI, BI, DI, EI,TMPOLD
      DOUBLE PRECISION :: AA,Ca ,Pa
      DOUBLE PRECISION, PARAMETER :: EPSL=1D-12
*     
*     Integrand subroutine
*
!MVNFUN Multivariate Normal integrand function
! where the integrand is transformed from an integral 
! having integration limits A and B to an
! integral having  constant integration limits i.e.
!   B                                               1 
! int ABS(xd(1)*xd(2)...xd(Nd))*f(xd,xt)dxt dxd = int F2(W) dW
!  A                                               0
!
! W    - new transformed integration variables, valid range 0..1
!        The vector must have the length Ndim returned from Covsrt
! COF  - conditional sorted ChOlesky Factor of the covariance matrix (IN)
! CDI  - Cholesky DIagonal elements used to calculate the mean  
! Cm   - conditional mean of Xd and Xt given Xc, E(Xd,Xt|Xc)
! xd   - variables to the jacobian variable, need no initialization size Nd
! xc   - conditional variables (IN)
! INDEX1 - if INDEX1(I)>Nt then variable No. I is one of the Xd variables otherwise it is one of Xt
      !PRINT *,'Mvnfun,ndim',Ndim
      VAL = VAL0
      
      N = Nt+Nd-INFIS-INFISD-1
       
      IF (INDEX1(1).GT.Nt ) THEN
         Y(1) = EXINV( D1 + W(1)*( E1 - D1 ),AA0,Ca0,Pa0)
      ELSE
         Y(1) = FIINV( D1 + W(1)*( E1 - D1 ) )   
      ENDIF
      
      NdleftN = NdleftN0
      NdleftO = Nd-INFISD 

      DO I=NdleftN+1,NdleftO
         VAL = VAL*ABS(CmXd(I)+Y(1)*CDIXd(I))
      ENDDO
      !PRINT *,'NdleftO',NdleftO,NdleftN
      NdleftO = NdleftN

      LK   = 0                  ! Counter for LK constraints
      IK   = 2                  ! Counter for Ndim 
      FINA = 0
      FINB = 0
      
     
      DO I = I0+1, N+1
         TMP = 0.d0
         JMX = MAX(I-IK,0)+1
         DO J = 1, I-JMX
            TMP = TMP + BIG(J,I)*Y(J)
         END DO
         IF (INFI(I).LT.0) GO TO 100            ! May have infinite int. Limits if Nd>0
         IF ( INFI(I) .NE. 0 ) THEN
            IF ( FINA .EQ. 1 ) THEN
               AI   = MAX( AI, A(I) - TMP) 
            ELSE
               !AI   = MAX(A(I) - TMP,-xCutOff) 
               AI   = A(I) - TMP 
               FINA = 1
            END IF
            IF (FINB.EQ.1.AND.BI.LE.AI) GOTO 200
         END IF
         IF ( INFI(I) .NE. 1 ) THEN
            IF ( FINB .EQ. 1 ) THEN
               BI   = MIN( BI, B(I) - TMP) 
            ELSE
               !BI   = MIN(B(I) - TMP,xCutOff)  
				BI   = B(I) - TMP  
               FINB = 1
            END IF
            IF (FINA.EQ.1.AND.BI.LE.AI) GOTO 200
         END IF
        
 100     IF (LK.LT.1) THEN
            TMPOLD = TMP
         ELSEIF   (INDEX1(I).GT.Nt) THEN ! Save the mean for Xd
            CmXd(NdleftN)  = Cm(I)+TMP*CDI(I)     ! Mean
            CDIXd(NdleftN) = CDI(I)
            !VAL = VAL*ABS(Cm(I)+TMP*CDI(I)) 
            NdleftN = NdleftN-1
         END IF
         IF ( I .EQ. N+1 .OR. BIG(IK+1,I+1) .GT.0.D0 ) THEN   
            IF (INDEX1(I-LK).GT.Nt) THEN
               AA =  Cm(I-LK)/CDI(I-LK)+TMPOLD  ! location parameter of inflection point
               CALL EXLMS(AA,AI, BI, 2*FINA+FINB-1, DI, EI,Ca,Pa)
               
               IF ( DI .GE. EI ) GO TO 200
               VAL = VAL*ABS(Ca*CDI(I-LK))
               IF ( I .LE. N .OR. NdleftN.LT.NdleftO) THEN
                  Y(IK) = EXINV(DI + W(IK)*( EI - DI ),AA,Ca,Pa)
                  !PRINT * ,'Y',Y(IK)
                  DO J = NdleftN+1,NdleftO
                     VAL = VAL*ABS(CmXd(J)+Y(IK)*CDIXd(J))
                  ENDDO
                  NdleftO = NdleftN
                  !PRINT *,'AA=',AA,' Ca=',Ca,' Pa=',Pa
               ENDIF   
               
            ELSE
               CALL MVNLMS( AI, BI, 2*FINA+FINB-1, DI, EI )            
               IF ( DI .GE. EI ) GO TO 200
               IF ( I .LE. N .OR. NdleftN.LT.NdleftO) THEN
                  Y(IK) = FIINV( DI + W(IK)*( EI - DI ) )
                  DO J = NdleftN+1,NdleftO
                     VAL = VAL*ABS(CmXd(J)+Y(IK)*CDIXd(J))
                  ENDDO
                  NdleftO = NdleftN
               ENDIF
            ENDIF            
            VAL = VAL*( EI - DI )
            !IF (ABS(VAL).GT.1.d0) THEN
            !   PRINT *, 'INDEX1(I-LK)',INDEX1(I-LK),Nt
            !   PRINT *, 'AA Ca,Pa',AA,Ca,Pa,DI,EI
            !   PRINT *, 'CDI,Cm,LK,I', CDI(I-LK),Cm(I-LK),LK,I
            !   PRINT *, 'TMP TMPOLD, Y(IK)',TMP, TMPOLD,Y(IK)
            !ENDIF
            !PRINT *,'I=',I, ' Y=',Y(IK),' val=',val
            LK   = 0
            IK   = IK + 1
            FINA = 0
            FINB = 0
         ELSE
            LK = LK+1
         END IF
      END DO
      RETURN
 200  VAL = 0.d0
      RETURN
      END FUNCTION MVNFUN
      END MODULE FUNCMOD2



!*****************************************************

     
      MODULE RIND
      IMPLICIT NONE
      PRIVATE
      PUBLIC :: RINDD, ECHO
!, INITDATA , SETDATA
 
      INTERFACE
         FUNCTION MVNFUN(N,Z) result (VAL)
         DOUBLE PRECISION,DIMENSION(:), INTENT(IN) :: Z
         INTEGER, INTENT(IN) :: N
         DOUBLE PRECISION :: VAL
         END FUNCTION MVNFUN
      END INTERFACE

      INTERFACE RINDD
      MODULE PROCEDURE RINDD
      END INTERFACE

      INTERFACE  BARRIER  
      MODULE PROCEDURE BARRIER
      END INTERFACE
      
      INTERFACE echo 
      MODULE PROCEDURE echo
      END INTERFACE
     
                       !--------------------------------
      CONTAINS
      

      SUBROUTINE ECHO(array)
      INTEGER ::j,i
      DOUBLE PRECISION,DIMENSION(:,:)::array
      CHARACTER*80 :: string
!      DO j=1,size(array,1)         
!         PRINT 111,j,array(j,:)
!111      FORMAT (i2,':',10F10.5)
!      END DO
      DO j=1,size(array,1)
         WRITE(string,110) j
 110     FORMAT (i2,':')
         !CALL mexprintf(string)
         DO i = 1, size(array,2)
            WRITE(string,111) array(j,i)
 111        FORMAT (' ',10F10.5)
          !  CALL mexprintf(string)
         ENDDO
         !CALL mexprintf(CHAR(10))   !CHAR(10) is a <CR> 
      END DO     
      END SUBROUTINE ECHO

     

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!******************* RINDD - the main program *********************!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE RINDD (VALS,ERR,Big1,Ex,Xc1,Nt1,
     &	              indI,Blo,Bup,INFIN,XcScale)  
      USE FUNCVARMOD
      USE GLOBALDATA 
      USE GLOBALCONST
      USE SWAPMOD
      USE COVSORTMOD
      USE FUNCMOD0              ! use FIINV AND C1C2 (REGRESSION EQ)
!     USE FUNCMOD               ! use FIINV 
!     USE FUNCMOD1              ! use MVNLIMITS and FIINV
!     USE FUNCMOD2              ! use EXINV and FIINV
      USE FIMOD
      USE RCRUDEMOD
      USE KRBVRCMOD
      USE ADAPTMOD
      USE KROBOVMOD
      IMPLICIT NONE  
      DOUBLE PRECISION, DIMENSION(:  ), INTENT(out):: VALS, ERR 
      DOUBLE PRECISION, DIMENSION(:,:), INTENT(in) :: BIG1
      DOUBLE PRECISION, DIMENSION(:,:), INTENT(in) :: xc1 
      DOUBLE PRECISION, DIMENSION(:  ), INTENT(in) :: Ex            
      DOUBLE PRECISION, DIMENSION(:,:), INTENT(in) :: Blo, Bup  
      INTEGER,          DIMENSION(:  ), INTENT(in) :: indI,INFIN
      INTEGER,									 INTENT(in) :: Nt1 
	DOUBLE PRECISION,						 INTENT(in) :: XcScale
! local variables
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: BIG2
      INTEGER, DIMENSION(:),ALLOCATABLE             :: XEDNI,INDEX2
      INTEGER :: K, I, J, ix, INFORM, Nc, Ntdc, Ntd, Nx, NDIM, Ndleft
     &	,MAXPTS,MINPTS
      DOUBLE PRECISION   :: VALUE,SQ0,xx,fxc,quant,ERROR,LABSEPS
	DOUBLE PRECISION   :: XcMin,XcMax,XcPrb, XcMinP
	INTEGER, DIMENSION(1) :: K1
      LOGICAL :: BCVSRT=.TRUE. ! sort covariance before hand if true

! index1, index2    = indices to the variables original place.   Size  Ntdc
! xedni     = indices to the variables new place.        Size  Ntdc
      
! Find the size information
!~~~~~~~~~~~~~~~~~~~~~~~~~~
      Nt   = Nt1
      Nc   = SIZE( Xc1, dim = 1 )
      Nx   = MAX( SIZE( Xc1, dim = 2), 1 )
      Ntdc = SIZE( BIG1, dim = 1 )
      IF (Nt+Nc.GT.Ntdc) Nt = Ntdc - Nc  ! make sure it does not exceed Ntdc-Nc
      Nd = Ntdc-Nt-Nc
      Ntd = Nt+Nd
      IF (Nd.LT.0) THEN
!         PRINT *,'RINDD Nt,Nd,Nc,Ntdc=',Nt,Nd,Nc,Ntdc
         STOP
      ENDIF
      
 !     PRINT *,'Nt Nd Nc Ntd Ntdc,',Nt, Nd, Nc, Ntd, Ntdc
     
! ALLOCATION
!~~~~~~~~~~~~
      IF (Nd.GT.0) THEN
         ALLOCATE(xd(Nd),CmXd(Nd),CDIXd(Nd))
         CmXd(:)  = gZERO
         CDIXd(:) = gZERO
         xd(:)    = gZERO
      END IF
	
      ALLOCATE(BIG(Ntdc,Ntdc),Cm(Ntdc),Y(Ntd))
      ALLOCATE(index1(Ntdc),A(Ntd),B(Ntd),INFI(Ntd),xc(1:Nc)) 
      ALLOCATE(CDI(Ntd),xedni(Ntdc),index2(Ntdc))
 
!     Initialization
!~~~~~~~~~~~~~~~~~~~~~

      BIG    = BIG1           !Copy input matrix
      VALS   = gZERO            ! SET VECTOR TO ZERO
      ERR    = gONE             ! SET Vector to ONE
      index2 = (/(J,J=1,Ntdc)/)
      
!      CALL mexprintf('BIG Before CovsrtXc'//CHAR(10))
!      CALL ECHO(BIG)
                                !sort BIG by decreasing cond. variance for Xc 
      CALL CVSRTXC(Nt,Nd,BIG,index2,INFORM) 
!      CALL mexprintf('BIG after CovsrtXc'//CHAR(10))
!      CALL ECHO(BIG)
      
      IF (INFORM.GT.0) GOTO 110 ! degenerate case exit VALS=0 for all  
                                ! (should perhaps return NaN instead??)


      DO I=Ntdc,1,-1     
         J = index2(I)            ! covariance matrix according to index2
         xedni(J) = I
      END DO

      IF (Nx.EQ.1) THEN
         BCVSRT=.FALSE.
      ELSE
         BCVSRT=.TRUE.
      ENDIF

      IF (BCVSRT) THEN          ! Conditionally sort all before hand       
         Cm = Ex (index2)
         IF (.FALSE.) THEN
            Xc = SUM(Xc1(1:Nc,1:Nx),DIM=2)/DBLE(Nx)
            I = Ntdc 
            DO J = 1, Nc        !Iterative conditioning on the last Nc variables  
               SQ0 = BIG(I,I)   ! SQRT(Var(X(i)|X(i+1),X(i+2),...,X(Ntdc)))
               xx = (Xc(index2(I)-Ntd)-Cm(I))/SQ0
                                ! conditional mean (expectation) 
                                ! E(X(1:i-1)|X(i),X(i+1),...,X(Ntdc)) 
               Cm(1:I-1) = Cm(1:I-1)+xx*BIG (1:I-1,I)
               I = I-1
            ENDDO
         ELSE
            I = Ntdc 
            DO J = 1, Nc        !Iterative conditioning on the last Nc variables  
               SQ0 = BIG(I,I)   ! SQRT(Var(X(i)|X(i+1),X(i+2),...,X(Ntdc)))
               K   = index2(I)-Ntd
               K1  = MINLOC(Xc1(K,:))
               XcMin = (Xc1(K,K1(1))-Cm(I))/SQ0
               K1    = MAXLOC(Xc1(K,:))
               XcMax = (Xc1(K,K1(1))-Cm(I))/SQ0
                                !CALL MVNLIMITS(XcMin,XcMax,2,XcMinP,XcPrb)
               CALL MVNLMS(XcMin,XcMax,2,XcMinP,XcPrb)
               XcPrb = XcPrb-XcMinP
               IF ( XcPrb .GT. gZERO) THEN
                  XcMin = -EXP(-gHALF*(XcMin**2))/gSQTWPI
                  XcMax = -EXP(-gHALF*(XcMax**2))/gSQTWPI
                  xx    = ( XcMax - XcMin )/XcPrb
               ELSE
                  xx = ( XcMin + XcMax )*gHALF
               END IF
               Xc(K) = xx*SQ0 + Cm(I)
                                ! conditional mean (expectation) 
                                ! E(X(1:i-1)|X(i),X(i+1),...,X(Ntdc)) 
               Cm(1:I-1) = Cm(1:I-1)+xx*BIG (1:I-1,I)
               I = I-1
            ENDDO            
         END IF
		
         CALL BARRIER(XEDNI,Cm,Xc,indI,Blo,Bup,INFIN,A,B,INFI) ! compute average integrationlimits   
         CALL COVSRT(BCVSRT,Nt,Nd,BIG,Cm,A,B,INFI,INDEX2, 
     &        INFIS,INFISD,NDIM,Y,CDI)
                                !print *,'xedni',xedni
                                !print *,'index2',index2
         INDEX1 = INDEX2
         DO I = Ntd, 1 , -1     
            J        = index2(I) ! covariance matrix according to index2
            xedni(J) = I
         END DO	
      ELSE
         IF (Nx>1) THEN
            ALLOCATE(BIG2(Ntdc,Ntdc))
            BIG2 = BIG
         ENDIF
      ENDIF
                                ! Now the loop over all different values of
                                ! variables Xc (the one one is conditioning on)
      DO  ix = 1, Nx            ! is started. The density f_{Xc}(xc(:,ix))                            
                                ! will be computed and denoted by  fxc.
         VALUE = gZERO                                                       
         fxc   = gONE
         ERROR = 2.d-16
                                ! Set the original means of the variables
         Cm    = Ex(index2(1:Ntdc)) !   Cm(1:Ntdc)  =Ex (index1(1:Ntdc))
         IF (Nc.GT.0) THEN
            Xc    = xc1(1:Nc,ix)
            QUANT = gZERO
            I = Ntdc 
            DO J = 1, Nc        !Iterative conditioning on the last Nc variables  
               SQ0 = BIG(I,I)   ! SQRT(Var(X(i)|X(i+1),X(i+2),...,X(Ntdc)))
               xx = (Xc(index2(I) - Ntd) - Cm(I))/SQ0
                                !Trick to calculate
                                !fxc = fxc*SQTWPI1*EXP(-0.5*(XX**2))/SQ0   
               QUANT = QUANT - gHALF*xx*xx + LOG(gSQTWPI1) - LOG(SQ0)              
				
                                ! conditional mean (expectation) 
                                ! E(X(1:i-1)|X(i),X(i+1),...,X(Ntdc)) 
               Cm(1:I-1) = Cm(1:I-1) + xx*BIG(1:I-1,I)
               I = I-1
            ENDDO
                                ! Calculating the  
                                ! fxc probability density for i=Ntdc-J+1, 
                                ! fXc=f(X(i)|X(i+1),X(i+2)...X(Ntdc))*
                                !     f(X(i+1)|X(i+2)...X(Ntdc))*..*f(X(Ntdc))
            fxc = EXP(QUANT+XcScale) 
            
            ERROR = ERROR*fxc
            IF (fxc .LT.fxcEpss) GOTO 100 ! Small probability don't bother calculating it
         END IF
                                !Set integration limits A,B and INFI
                                !NOTE: A and B are integration limits with Cm subtracted
         CALL BARRIER(xedni,Cm,xc,indI,Blo,Bup,INFIN,A,B,INFI)
         
         IF (BCVSRT) THEN       ! conditionally sorted before hand just scale integration limits
            WHERE (ABS(CDI).GT.gZERO) ! EPS2)
               A = A/CDI
               B = B/CDI
            END WHERE
                                !The following could be done smarter by changing Blo Bup and INFIN
            DO I = 1,Ntd
               IF ( CDI(I).LT.gZERO.AND.INFI(I).GE.0) THEN
                  CALL SWAP( A(I), B(I) ) 
                  IF ( INFI(I) .NE. 2 ) INFI(I) = 1 - INFI(I)
               END IF
            ENDDO
         ELSE
            INDEX1 = index2   
            IF (ix>1) BIG = BIG2
            CALL COVSRT(BCVSRT,Nt,Nd,BIG,Cm,A,B,INFI, 
     &           INDEX1,INFIS,INFISD,NDIM,Y,CDI)
         ENDIF

         CALL INITFUN(VALUE,INFORM)
!     IF INFORM>0 : degenerate case:
!     Integral can be calculated excactly, ie. 
!     mean of deterministic variables outside the barriers,
!     or NDIM = 1     
         IF (INFORM.GT.0) GO TO 100 
         
         MAXPTS  = NSIMmax
         MINPTS  = NSIMmin
         LABSEPS = ABSEPS       !*fxc
      
         SELECT CASE (SCIS)
         CASE (:1)
            IF (NDIM.LT.9) THEN
               CALL SADAPT(NDIM,MAXPTS,MVNFUN,LABSEPS,
     &              RELEPS,ERROR,VALUE,INFORM)
               VALUE = MAX(VALUE,gZERO)
            ELSE
               CALL KRBVRC(NDIM, MINPTS, MAXPTS, MVNFUN, LABSEPS,RELEPS,
     &              ERROR, VALUE, INFORM )
            ENDIF
         CASE (2)               !        Call the subregion adaptive integration subroutine
            IF ( NDIM .GT. 19.) THEN
!     print *, 'Ndim too large for SADMVN => Calling KRBVRC'
               CALL KRBVRC( NDIM, MINPTS, MAXPTS, MVNFUN, LABSEPS,
     &              RELEPS, ERROR, VALUE, INFORM )
            ELSE
               CALL SADAPT(NDIM,MAXPTS,MVNFUN,LABSEPS,
     &              RELEPS,ERROR,VALUE,INFORM)
               VALUE = MAX(VALUE,gZERO)
            ENDIF
         CASE (3)               !       Call the Lattice rule integration procedure
            CALL KRBVRC( NDIM, MINPTS, MAXPTS, MVNFUN, LABSEPS,
     &           RELEPS, ERROR, VALUE, INFORM )
         CASE (4)               !       Call the Lattice rule integration procedure
            CALL KROBOV( NDIM, MINPTS, MAXPTS, MVNFUN, LABSEPS, RELEPS,
     &           ERROR, VALUE, INFORM )
         CASE (5:)              ! Call Crude Monte Carlo integration procedure
            CALL RANMC( NDIM, MAXPTS, MVNFUN, LABSEPS, 
     &           RELEPS, ERROR, VALUE, INFORM )
         END SELECT   

!     IF (INFORM.gt.0) print *,'RINDD, INFORM,error =',inform,error
 100     VALS(ix) = VALUE*fxc 
         IF (SIZE(ERR, DIM = 1).EQ.Nx) ERR(ix)  = error*fxc       
!     PRINT *,'RINDD MINPTS=',MINPTS
      ENDDO                     !ix

 110  CONTINUE
      IF (ALLOCATED(xc))     DEALLOCATE(xc) 
      IF (ALLOCATED(xd))     DEALLOCATE(xd) 
      IF (ALLOCATED(Cm))     DEALLOCATE(Cm) 
      IF (ALLOCATED(BIG2))   DEALLOCATE(BIG2) 
      IF (ALLOCATED(BIG))    DEALLOCATE(BIG) 
      IF (ALLOCATED(index2)) DEALLOCATE(index2)
      IF (ALLOCATED(index1)) DEALLOCATE(index1)
      IF (ALLOCATED(xedni))  DEALLOCATE(xedni)
      IF (ALLOCATED(A))      DEALLOCATE(A)
      IF (ALLOCATED(B))      DEALLOCATE(B)  
      IF (ALLOCATED(Y))      DEALLOCATE(Y) 
      IF (ALLOCATED(CDI))    DEALLOCATE(CDI)  
      IF (ALLOCATED(CDIXd))  DEALLOCATE(CDIXd)  
      IF (ALLOCATED(CmXd))   DEALLOCATE(CmXd)  
      IF (ALLOCATED(INFI))   DEALLOCATE(INFI)
      RETURN                                                            
      END SUBROUTINE RINDD
      
      SUBROUTINE BARRIER(xedni,Cm,xc,indI,Blo,Bup,INFIN,Hlo,Hup,INFI)     
      USE GLOBALDATA, ONLY : xCutOff
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(:  ), INTENT(in)  :: xc,Cm 
      INTEGER,          DIMENSION(:  ), INTENT(in)  :: xedni,indI,INFIN    
      DOUBLE PRECISION, DIMENSION(:,:), INTENT(in)  :: Blo,Bup   
      DOUBLE PRECISION, DIMENSION(:  ), INTENT(out) :: Hlo,Hup
      INTEGER,          DIMENSION(:  ), INTENT(out) :: INFI  
!Local variables
      INTEGER   :: I, J, K, L,Mb,Nb,NI,Nc           
!this procedure set Hlo,Hup and INFI according to Blo/Bup and INFIN
!
!     INFIN  INTEGER, array of integration limits flags:
!            if INFIN(I) < 0, Ith limits are (-infinity, infinity);
!            if INFIN(I) = 0, Ith limits are (-infinity, Hup(I)];
!            if INFIN(I) = 1, Ith limits are [Hlo(I), infinity);
!            if INFIN(I) = 2, Ith limits are [Hlo(I), Hup(I)].
! Note :
! xedni     = indices to the variables new place.        Size  Ntdc
 
      Mb=size(Blo,DIM=1)
      Nb=size(Blo,DIM=2)
      NI=size(indI,DIM=1)
      Nc=size(xc,DIM=1)
!      IF (Mb.GT.Nc+1) print *,'barrier: Mb,Nc =',Mb,Nc
!      IF (Nb.NE.NI-1) print *,'barrier: Nb,NI =',Nb,NI
      DO J = 2, NI 
         DO I =indI (J - 1) + 1 , indI (J)
            L = xedni(I)  
            INFI(L)=INFIN(J-1)
            Hlo (L) = -xCutOff
            Hup (L) =  xCutOff
            IF (INFI(L).GE.0) THEN
               IF  (INFI(L).NE.0) THEN
                  Hlo (L) = Blo (1, J - 1)-Cm(L)  
                  DO K = 1, Mb-1                                                  
                     Hlo(L) = Hlo(L)+Blo(K+1,J-1)*xc(K)              
                  ENDDO         ! K
               ENDIF
               IF  (INFI(L).NE.1) THEN
                  Hup (L) = Bup (1, J - 1)-Cm(L)  
                  DO K = 1, Mb-1                                                  
                     Hup(L) = Hup(L)+Bup(K+1,J-1)*xc(K)
                  ENDDO  
               ENDIF                 !
            ENDIF            
         ENDDO ! I
      ENDDO ! J
!      print * ,'barrier hup:',size(Hup),Hup(xedni(1:indI(NI)))
!      print * ,'barrier hlo:',size(Hlo),Hlo(xedni(1:indI(NI)))
      RETURN                                                             
      END SUBROUTINE BARRIER                                              
      
      END MODULE RIND         !******************************  
