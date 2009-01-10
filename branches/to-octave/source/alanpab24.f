! Programs available in module RIND : 
!
!   1) setConstants
!   2) RINDD
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
!CALL RINDD(E,err,S,m,xc,Nt,indI,Blo,Bup,INFIN);
!
!        E = expectation/density as explained above size 1 x Nx (out)
!      ERR = estimated absolute error  size 1 x Nx (out)
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
! The indices, indI=[0 3 5 6], and coefficients Blo=[0 0 -1], 
! Bup=[0 0 5], INFIN=[0 1 2] 
! means that   A = [-inf -inf -inf 0 0 -1]  B = [0 0 0 inf inf 5] 
!
!
! (Recommended limitations Nx,Nt<101, Nd<7 and Nc<11)
! Also note that the size information have to be transferred to RINDD
! through the input arguments E,S,m,Nt,IndI,Blo,Bup and INFIN 
!
! For further description see the modules 
!
! References
! Podgorski et al. (2000)
! "Exact distributions for apparent waves in irregular seas"
! Ocean Engineering,  Vol 27, no 1, pp979-1016.                        (RINDD)
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
! Alan Genz and Koon-Shing Kwong (2000?)
! 'Numerical Evaluation of Singular Multivariate Normal Distributions' (MVNFUN,COVSRT)
! Computational Statistics and Data analysis
!
!
! P. A. Brodtkorb (2004),                                 (RINDD, MVNFUN, COVSRT)
! Numerical evaluation of multinormal expectations
! In Lund university report series
! and in the Dr.Ing thesis: 
! The probability of Occurrence of dangerous Wave Situations at Sea.
! Dr.Ing thesis, Norwegian University of Science and Technolgy, NTNU,
! Trondheim, Norway.

! Tested on:  DIGITAL UNIX Fortran90 compiler
!             PC pentium II with Lahey Fortran90 compiler
!             Solaris with SunSoft F90 compiler Version 1.0.1.0  (21229283) 
! History:
! revised pab 23may2004
! RIND module totally rewritten according to the last reference.


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

      MODULE PRINTMOD
      IMPLICIT NONE
!      PUBLIC
! PRINTMOD Module containing print functions to screen in matlab mexfunctions

      INTERFACE echo 
      MODULE PROCEDURE echo
      END INTERFACE
     
      INTERFACE PRINTCOF
      MODULE PROCEDURE PRINTCOF
      END INTERFACE 

      INTERFACE PRINTVAR
      MODULE PROCEDURE PRINTVARD, PRINTVARI
      END INTERFACE 

      INTERFACE PRINTVEC
      MODULE PROCEDURE PRINTVECD, PRINTVECI
      END INTERFACE 


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
         CALL mexprintf(string)
         DO i = 1, size(array,2)
            WRITE(string,111) array(j,i)
 111        FORMAT (' ',10F10.5)
            CALL mexprintf(string)
         ENDDO
         CALL mexprintf(CHAR(10))   !CHAR(10) is a <CR> 
      END DO     
      END SUBROUTINE ECHO

      SUBROUTINE PRINTVARI(I,TXT)
      INTEGER,  INTENT(IN) :: I
      CHARACTER*80, OPTIONAL, INTENT(IN) :: TXT
      CHARACTER*80 :: string1
      IF (PRESENT(TXT)) THEN
         CALL mexprintf(TXT)
         CALL mexprintf(':')
      ENDIF
      WRITE(string1,125) I
 125  FORMAT (' ',i5) 
      CALL mexprintf(string1//CHAR(10)) 
      RETURN
      END SUBROUTINE PRINTVARI
      SUBROUTINE PRINTVARD(D,TXT)
      DOUBLE PRECISION,  INTENT(IN) :: D
      CHARACTER*80, OPTIONAL, INTENT(IN) :: TXT
      CHARACTER*80 :: string1
      IF (PRESENT(TXT)) THEN
         CALL mexprintf(TXT)
         CALL mexprintf(':')
      ENDIF
      WRITE(string1,115) D
 115  FORMAT (' ',10F10.5) 
      CALL mexprintf(string1) 
      CALL mexprintf(CHAR(10))
      RETURN
      END SUBROUTINE PRINTVARD

      SUBROUTINE PRINTVECD(CDI,TXT)
      DOUBLE PRECISION, DIMENSION(:),INTENT(in) :: CDI
      CHARACTER*80, OPTIONAL, INTENT(IN) :: TXT
      INTEGER :: I
      CHARACTER*80 :: string
      IF (PRESENT(TXT)) THEN
         CALL mexprintf(TXT) 
         CALL mexprintf(':')
      ENDIF
      DO I = 1, SIZE(CDI,1)
         WRITE(string,115) CDI(I)
 115     FORMAT (' ',10F10.5) 
         CALL mexprintf(string) 
      ENDDO
      CALL mexprintf(CHAR(10)) 
      RETURN
      END SUBROUTINE PRINTVECD
      SUBROUTINE PRINTVECI(CDI,TXT)
      INTEGER, DIMENSION(:),INTENT(in) :: CDI
      CHARACTER*80, OPTIONAL, INTENT(IN) :: TXT
      INTEGER :: I
      CHARACTER*80 :: string
      IF (PRESENT(TXT)) THEN
         CALL mexprintf(TXT) 
         CALL mexprintf(':')
      ENDIF
      DO I = 1, SIZE(CDI,1)
         WRITE(string,115) CDI(I)
 115     FORMAT (' ',i5) 
         CALL mexprintf(string) 
      ENDDO
      CALL mexprintf(CHAR(10)) 
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
         IF ( INFI(I)  <  0 ) THEN 
            WRITE(string,111) INDEX1(I)
 111        FORMAT (i2,'  -inf    inf '\) 
!            CALL mexprintf(string) 
         ELSE IF ( INFI(I) .EQ. 0 ) THEN 
            WRITE(string,112) INDEX1(I),B(I)
 112        FORMAT (i2,'  -inf ', 5F5.2)
!            CALL mexprintf(string) 
         ELSE IF ( INFI(I) .EQ. 1 ) THEN 
            WRITE(string,113) INDEX1(I),A(I)
 113        FORMAT (i2,' ', 5F5.2 ,'   inf '\)
!            CALL mexprintf(string) 
!            CALL mexprintf('   inf ')
         ELSE 
            WRITE(string,114) INDEX1(I),A(I),B(I)
 114        FORMAT (i2,' ', 5F5.2 ,' ', 5F5.2)
!           CALL mexprintf(string) 
         END IF
         CALL mexprintf(string)
         DO K = 1,J
            WRITE(string,115) COF(I,K)
 115        FORMAT (' ',10F10.5) 
            CALL mexprintf(string) 
         ENDDO
         CALL mexprintf(CHAR(10)) 
      END DO
      RETURN
      END SUBROUTINE PRINTCOF
      end module PRINTMOD
      
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
!C--------------------------------------------------------------------
!C
!C DERF subprogram computes approximate values for erf(x).
!C   (see comments heading CALERF).
!C
!C   Author/date: W. J. Cody, January 8, 1985
!C
!C--------------------------------------------------------------------
      FUNCTION DERF( X ) RESULT (VALUE)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN)  :: X
      DOUBLE PRECISION   :: VALUE
      INTEGER, PARAMETER :: JINT = 0
      CALL CALERF(X,VALUE,JINT)
      RETURN
      END FUNCTION DERF
!C--------------------------------------------------------------------
!C
!C DERFC subprogram computes approximate values for erfc(x).
!C   (see comments heading CALERF).
!C
!C   Author/date: W. J. Cody, January 8, 1985
!C
!C--------------------------------------------------------------------
      FUNCTION DERFC( X ) RESULT (VALUE)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN)  :: X
      DOUBLE PRECISION :: VALUE
      INTEGER, PARAMETER :: JINT = 1
      CALL CALERF(X,VALUE,JINT)
      RETURN
      END FUNCTION DERFC
!C------------------------------------------------------------------
!C
!C DERFCX subprogram computes approximate values for exp(x*x) * erfc(x).
!C   (see comments heading CALERF).
!C
!C   Author/date: W. J. Cody, March 30, 1987
!C
!C------------------------------------------------------------------
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
!C------------------------------------------------------------------
!C
!C CALERF packet evaluates  erf(x),  erfc(x),  and  exp(x*x)*erfc(x)
!C   for a real argument  x.  It contains three FUNCTION type
!C   subprograms: ERF, ERFC, and ERFCX (or DERF, DERFC, and DERFCX),
!C   and one SUBROUTINE type subprogram, CALERF.  The calling
!C   statements for the primary entries are:
!C
!C                   Y=ERF(X)     (or   Y=DERF(X)),
!C
!C                   Y=ERFC(X)    (or   Y=DERFC(X)),
!C   and
!C                   Y=ERFCX(X)   (or   Y=DERFCX(X)).
!C
!C   The routine  CALERF  is intended for internal packet use only,
!C   all computations within the packet being concentrated in this
!C   routine.  The function subprograms invoke  CALERF  with the
!C   statement
!C
!C          CALL CALERF(ARG,RESULT,JINT)
!C
!C   where the parameter usage is as follows
!C
!C      Function                     Parameters for CALERF
!C       call              ARG                  Result          JINT
!C
!C     ERF(ARG)      ANY REAL ARGUMENT         ERF(ARG)          0
!C     ERFC(ARG)     ABS(ARG) .LT. XBIG        ERFC(ARG)         1
!C     ERFCX(ARG)    XNEG .LT. ARG .LT. XMAX   ERFCX(ARG)        2
!C
!C   The main computation evaluates near-minimax approximations
!C   from "Rational Chebyshev approximations for the error function"
!C   by W. J. Cody, Math. Comp., 1969, PP. 631-638.  This
!C   transportable program uses rational functions that theoretically
!C   approximate  erf(x)  and  erfc(x)  to at least 18 significant
!C   decimal digits.  The accuracy achieved depends on the arithmetic
!C   system, the compiler, the intrinsic functions, and proper
!C   selection of the machine-dependent constants.
!C
!C*******************************************************************
!C*******************************************************************
!C
!C Explanation of machine-dependent constants
!C
!C   XMIN   = the smallest positive floating-point number.
!C   XINF   = the largest positive finite floating-point number.
!C   XNEG   = the largest negative argument acceptable to ERFCX;
!C            the negative of the solution to the equation
!C            2*exp(x*x) = XINF.
!C   XSMALL = argument below which erf(x) may be represented by
!C            2*x/sqrt(pi)  and above which  x*x  will not underflow.
!C            A conservative value is the largest machine number X
!C            such that   1.0 + X = 1.0   to machine precision.
!C   XBIG   = largest argument acceptable to ERFC;  solution to
!C            the equation:  W(x) * (1-0.5/x**2) = XMIN,  where
!C            W(x) = exp(-x*x)/[x*sqrt(pi)].
!C   XHUGE  = argument above which  1.0 - 1/(2*x*x) = 1.0  to
!C            machine precision.  A conservative value is
!C            1/[2*sqrt(XSMALL)]
!C   XMAX   = largest acceptable argument to ERFCX; the minimum
!C            of XINF and 1/[sqrt(pi)*XMIN].
!C
!C   Approximate values for some important machines are:
!C
!C                          XMIN       XINF        XNEG     XSMALL
!C
!C    C 7600      (S.P.)  3.13E-294   1.26E+322   -27.220  7.11E-15
!C  CRAY-1        (S.P.)  4.58E-2467  5.45E+2465  -75.345  7.11E-15
!C  IEEE (IBM/XT,
!C    SUN, etc.)  (S.P.)  1.18E-38    3.40E+38     -9.382  5.96E-8
!C  IEEE (IBM/XT,
!C    SUN, etc.)  (D.P.)  2.23D-308   1.79D+308   -26.628  1.11D-16
!C  IBM 195       (D.P.)  5.40D-79    7.23E+75    -13.190  1.39D-17
!C  UNIVAC 1108   (D.P.)  2.78D-309   8.98D+307   -26.615  1.73D-18
!C  VAX D-Format  (D.P.)  2.94D-39    1.70D+38     -9.345  1.39D-17
!C  VAX G-Format  (D.P.)  5.56D-309   8.98D+307   -26.615  1.11D-16
!C
!C
!C                          XBIG       XHUGE       XMAX
!C
!C    C 7600      (S.P.)  25.922      8.39E+6     1.80X+293
!C  CRAY-1        (S.P.)  75.326      8.39E+6     5.45E+2465
!C  IEEE (IBM/XT,
!C    SUN, etc.)  (S.P.)   9.194      2.90E+3     4.79E+37
!C  IEEE (IBM/XT,
!C    SUN, etc.)  (D.P.)  26.543      6.71D+7     2.53D+307
!C  IBM 195       (D.P.)  13.306      1.90D+8     7.23E+75
!C  UNIVAC 1108   (D.P.)  26.582      5.37D+8     8.98D+307
!C  VAX D-Format  (D.P.)   9.269      1.90D+8     1.70D+38
!C  VAX G-Format  (D.P.)  26.569      6.71D+7     8.98D+307
!C
!C*******************************************************************
!C*******************************************************************
!C
!C Error returns
!C
!C  The program returns  ERFC = 0      for  ARG .GE. XBIG;
!C
!C                       ERFCX = XINF  for  ARG .LT. XNEG;
!C      and
!C                       ERFCX = 0     for  ARG .GE. XMAX.
!C
!C
!C Intrinsic functions required are:
!C
!C     ABS, AINT, EXP
!C
!C
!C  Author: W. J. Cody
!C          Mathematics and Computer Science Division
!C          Argonne National Laboratory
!C          Argonne, IL 60439
!C
!C  Latest modification: March 19, 1990
!C  Updated to F90 by pab 23.03.2003
!C
!C------------------------------------------------------------------
      DOUBLE PRECISION, INTENT(IN) :: ARG
      INTEGER, INTENT(IN)          :: JINT
      DOUBLE PRECISION, INTENT(INOUT):: RESULT
! Local variables
      INTEGER :: I
      DOUBLE PRECISION :: DEL,X,XDEN,XNUM,Y,YSQ
!C------------------------------------------------------------------
!C  Mathematical constants
!C------------------------------------------------------------------
      DOUBLE PRECISION, PARAMETER :: ZERO   = 0.0D0
      DOUBLE PRECISION, PARAMETER :: HALF   = 0.05D0
      DOUBLE PRECISION, PARAMETER :: ONE    = 1.0D0
      DOUBLE PRECISION, PARAMETER :: TWO    = 2.0D0
      DOUBLE PRECISION, PARAMETER :: FOUR   = 4.0D0
      DOUBLE PRECISION, PARAMETER :: SIXTEN = 16.0D0
      DOUBLE PRECISION, PARAMETER :: SQRPI  = 5.6418958354775628695D-1
      DOUBLE PRECISION, PARAMETER :: THRESH = 0.46875D0
!C------------------------------------------------------------------
!C  Machine-dependent constants
!C------------------------------------------------------------------
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
!C------------------------------------------------------------------
!C  Coefficients for approximation to  erf  in first interval
!C------------------------------------------------------------------
      PARAMETER (A = (/ 3.16112374387056560D00,
     &     1.13864154151050156D02,3.77485237685302021D02,
     &     3.20937758913846947D03, 1.85777706184603153D-1/))
      PARAMETER ( B = (/2.36012909523441209D01,2.44024637934444173D02,
     &       1.28261652607737228D03,2.84423683343917062D03/))
!C------------------------------------------------------------------
!C  Coefficients for approximation to  erfc  in second interval
!C------------------------------------------------------------------
      PARAMETER ( C=(/5.64188496988670089D-1,8.88314979438837594D0,
     1       6.61191906371416295D01,2.98635138197400131D02,
     2       8.81952221241769090D02,1.71204761263407058D03,
     3       2.05107837782607147D03,1.23033935479799725D03,
     4       2.15311535474403846D-8/))
      PARAMETER ( D =(/1.57449261107098347D01,1.17693950891312499D02, 
     1       5.37181101862009858D02,1.62138957456669019D03,
     2       3.29079923573345963D03,4.36261909014324716D03,
     3       3.43936767414372164D03,1.23033935480374942D03/))
!C------------------------------------------------------------------
!C  Coefficients for approximation to  erfc  in third interval
!C------------------------------------------------------------------
      PARAMETER ( P =(/3.05326634961232344D-1,3.60344899949804439D-1,
     1       1.25781726111229246D-1,1.60837851487422766D-2,
     2       6.58749161529837803D-4,1.63153871373020978D-2/))
      PARAMETER (Q =(/2.56852019228982242D00,1.87295284992346047D00,
     1       5.27905102951428412D-1,6.05183413124413191D-2,
     2       2.33520497626869185D-3/))
!C------------------------------------------------------------------
      X = ARG
      Y = ABS(X)
      IF (Y .LE. THRESH) THEN
!C------------------------------------------------------------------
!C  Evaluate  erf  for  |X| <= 0.46875
!C------------------------------------------------------------------
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
!C------------------------------------------------------------------
!C     Evaluate  erfc  for 0.46875 <= |X| <= 4.0
!C------------------------------------------------------------------
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
!C------------------------------------------------------------------
!C     Evaluate  erfc  for |X| > 4.0
!C------------------------------------------------------------------
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
!C------------------------------------------------------------------
!C     Fix up for negative argument, erf, etc.
!C------------------------------------------------------------------
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
      MODULE TRIVARIATEVAR 
!     Global variables used in calculation of TRIVARIATE
!     normal and TRIVARIATE student T probabilties
      INTEGER :: NU
      DOUBLE PRECISION :: H1, H2, H3, R23, RUA, RUB, AR, RUC
      END MODULE TRIVARIATEVAR 
!
! FIMOD contains functions for calculating 1D and 2D Normal probabilites
!       and  1D expectations
      MODULE FIMOD
      IMPLICIT NONE
      PRIVATE
      PUBLIC :: NORMPRB, FI, FIINV, MVNLIMITS, MVNLMS, BVU,BVNMVN
      PUBLIC :: GAUSINT, GAUSINT2, EXLMS, EXINV
      PUBLIC :: STUDNT, BVTL, TVTL, TVNMVN

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

      INTERFACE STUDNT
      MODULE PROCEDURE STUDNT
      END INTERFACE

      INTERFACE BVTL
      MODULE PROCEDURE BVTL
      END INTERFACE

      INTERFACE TVTL
      MODULE PROCEDURE TVTL
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
         IF ( Q  <  ZERO ) VAL = - VAL
      END IF
      RETURN
      END FUNCTION FIINV     
                                ! *********************************     
      SUBROUTINE NORMPRB(Z, P, Q)
      USE ERFCOREMOD
!      USE GLOBALDATA, ONLY : XMAX      
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
      DOUBLE PRECISION, PARAMETER :: XMAX  = 37D0
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

         IF (Z  <  ZERO) THEN
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
!      USE GLOBALDATA, ONLY : XMAX
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in) :: Z
      DOUBLE PRECISION :: VALUE
! Local variables
      DOUBLE PRECISION :: ZABS
      DOUBLE PRECISION, PARAMETER:: SQ2M1 = 0.70710678118655D0 !     1/SQRT(2)
      DOUBLE PRECISION, PARAMETER:: HALF = 0.5D0
      DOUBLE PRECISION, PARAMETER:: XMAX = 37.D0
      ZABS = ABS(Z)
*     
*     |Z| > 37  (or XMAX)
*     
      IF ( ZABS .GT. XMAX ) THEN
         IF (Z  <  0.0D0) THEN
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
!      USE GLOBALDATA, ONLY : XMAX
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
     *     Q0, Q1, Q2, Q3, Q4, Q5, Q6, Q7,XMAX,
     *     P, EXPNTL, CUTOFF, ROOTPI, ZABS, Z2
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
      PARAMETER( XMAX   = 37.D0 )
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
         Z2 = ZABS * ZABS
         EXPNTL = EXP( -Z2 * 0.5D0 )
*     
*     |Z| < CUTOFF = 10/SQRT(2)
*     
         IF ( ZABS  <  CUTOFF ) THEN
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
            IF (AP+BP  <  ONE) THEN
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
      IF ( INFIN  <  0 ) RETURN
      IF ( INFIN .NE. 0 ) LOWER = FI(A)
      IF ( INFIN .NE. 1 ) UPPER = FI(B)
      RETURN
      END SUBROUTINE MVNLMS 

      FUNCTION TVNMVN(A, B, INFIN, R, EPSI ) RESULT (VAL)
      IMPLICIT NONE
*
*     A function for computing trivariate normal probabilities.
*
*  Parameters
*
*     A  REAL, array of lower integration limits.
*     B  REAL, array of upper integration limits.
!     R  REAL, array of correlation coefficents
!         R = [r12 r13 r23]  
!    EPSI = REAL tolerance
*     INFIN  INTEGER, array of integration limits flags:
*            if INFIN(I) = 0, Ith limits are (-infinity, B(I)];
*            if INFIN(I) = 1, Ith limits are [A(I), infinity);
*            if INFIN(I) = 2, Ith limits are [A(I), B(I)].
      DOUBLE PRECISION, DIMENSION(:), INTENT (IN) :: A, B , R
      DOUBLE PRECISION,               INTENT (IN) :: EPSI
      INTEGER,          DIMENSION(:), INTENT (IN) :: INFIN
      DOUBLE PRECISION :: VAL
      
      IF ( INFIN(1) .EQ. 2 ) THEN
         IF (  INFIN(2) .EQ. 2 ) THEN
            IF (INFIN(3) .EQ. 2 ) THEN !OK
               VAL =  TVNL( B(1), B(2), B(3),R(1),R(2),R(3),EPSI ) 
     &              - TVNL( A(1), B(2), B(3),R(1),R(2),R(3),EPSI )
     &              - TVNL( B(1), A(2), B(3),R(1),R(2),R(3),EPSI )  
     &              - TVNL( B(1), B(2), A(3),R(1),R(2),R(3),EPSI )
     &              + TVNL( A(1), A(2), B(3),R(1),R(2),R(3),EPSI )
     &              + TVNL( A(1), B(2), A(3),R(1),R(2),R(3),EPSI )
     $              + TVNL( B(1), A(2), A(3),R(1),R(2),R(3),EPSI )
     &              - TVNL( A(1), A(2), A(3),R(1),R(2),R(3),EPSI )  
            ELSE IF (INFIN(3) .EQ. 1 ) THEN ! B(3) = inf ok
               VAL = TVNL( B(1), B(2), -A(3),R(1),-R(2),-R(3),EPSI ) 
     &              - TVNL( A(1), B(2), -A(3),R(1),-R(2),-R(3),EPSI )
     &              - TVNL( B(1), A(2), -A(3),R(1),-R(2),-R(3),EPSI )  
     &              + TVNL( A(1), A(2), -A(3),R(1),-R(2),-R(3),EPSI )
            ELSE IF (INFIN(3) .EQ. 0 ) THEN !OK A(3) = -inf
               VAL =  TVNL( B(1), B(2), B(3),R(1),R(2),R(3),EPSI ) 
     &              - TVNL( A(1), B(2), B(3),R(1),R(2),R(3),EPSI )
     &              - TVNL( B(1), A(2), B(3),R(1),R(2),R(3),EPSI )  
     &              + TVNL( A(1), A(2), B(3),R(1),R(2),R(3),EPSI )
            ELSE ! INFIN(1:2)=2
               VAL = BVNMVN( A ,B, INFIN, R(1) )
            ENDIF
         ELSE IF (INFIN(2) .EQ. 1 ) THEN ! B(2) = inf
            IF (INFIN(3) .EQ. 2 ) THEN 
               VAL = TVNL( B(1), -A(2), B(3),-R(1),R(2),-R(3),EPSI ) 
     &              - TVNL( A(1), -A(2), B(3),-R(1),R(2),-R(3),EPSI )
     &              - TVNL( B(1), -A(2), A(3),-R(1),R(2),-R(3),EPSI )  
     &              + TVNL( A(1), -A(2), A(3),-R(1),R(2),-R(3),EPSI )
            ELSE IF (INFIN(3) .EQ. 1 ) THEN
               VAL = TVNL( B(1), -A(2), -A(3),-R(1),-R(2),R(3),EPSI ) 
     $              - TVNL( A(1), -A(2), -A(3),-R(1),-R(2),R(3),EPSI )
            ELSE IF (INFIN(3) .EQ. 0 ) THEN 
               VAL = TVNL( B(1), -A(2), B(3),-R(1),R(2),-R(3),EPSI) 
     $              - TVNL( A(1), -A(2), B(3),-R(1),R(2),-R(3),EPSI)
            ELSE
               VAL = BVNMVN( A ,B, INFIN, R(1) )
            ENDIF
         ELSE IF (INFIN(2) .EQ. 0 ) THEN
            SELECT CASE (INFIN(3))
            CASE (2:)           ! % % A(2)=-INF
               VAL = TVNL( B(1), B(2), B(3),R(1),R(2),R(3),EPSI ) 
     $              - TVNL( A(1), B(2), B(3),R(1),R(2),R(3),EPSI ) 
     $              - TVNL( B(1), B(2), A(3),R(1),R(2),R(3),EPSI ) 
     $              + TVNL( A(1), B(2), A(3),R(1),R(2),R(3),EPSI )
            CASE (1)            !%  % A(2)=-INF B(3) = INF
               VAL = TVNL( B(1), B(2), -A(3),R(1),-R(2),-R(3),EPSI)
     $              - TVNL(A(1), B(2), -A(3),R(1),-R(2),-R(3),EPSI)
            CASE (0)            ! % % A(2)=-INF A(3) = -INF
               VAL = TVNL( B(1), B(2), B(3),R(1),R(2),R(3),EPSI ) 
     $              - TVNL( A(1), B(2), B(3),R(1),R(2),R(3),EPSI)
            CASE DEFAULT
               VAL = BVNMVN(A,B,INFIN,R(1));
            END SELECT
         ELSE
            VAL = BVNMVN(A(1:3:2),B(1:3:2),INFIN(1:3:2),R(2));
         ENDIF
      ELSE IF ( INFIN(1) .EQ. 1 ) THEN
         SELECT CASE (INFIN(2))
         CASE (2) 
            SELECT CASE (INFIN(3))
            CASE (2)            !% B(1) = INF   %OK
               VAL =  TVNL(-A(1), B(2), B(3),-R(1),-R(2),R(3),EPSI )
     $              - TVNL(-A(1), B(2), A(3),-R(1),-R(2),R(3),EPSI )
     $              - TVNL(-A(1), A(2), B(3),-R(1),-R(2),R(3),EPSI )
     $              + TVNL(-A(1), A(2), A(3),-R(1),-R(2),R(3),EPSI )
            CASE (1)            ! % B(1) = INF   B(3) = INF %OK
               VAL = TVNL(-A(1), B(2), -A(3),-R(1),R(2),-R(3),EPSI )
     $              - TVNL(-A(1), A(2), -A(3),-R(1),R(2),-R(3),EPSI)
            CASE (0)            ! % B(1) = INF   A(3) = -INF %OK
               VAL =  TVNL(-A(1), B(2), B(3),-R(1),-R(2),R(3),EPSI )
     $              - TVNL(-A(1), A(2), B(3),-R(1),-R(2),R(3),EPSI)
            CASE (-1)
               VAL = BVNMVN(A,B,INFIN,R(1));
            END SELECT
         CASE (1)               !%B(2) = INF
            SELECT CASE (INFIN(3))
            CASE (2)            ! % B(1) = INF  B(2) = INF % OK
               VAL = TVNL( -A(1), -A(2), B(3),-R(1),R(2),-R(3),EPSI ) 
     &              - TVNL( -A(1), -A(2),A(3),-R(1),R(2),-R(3),EPSI );
            CASE (1)            ! % B(1:3) = INF %OK
               VAL = TVNL( -A(1), -A(2), -A(3),R(1),R(2),R(3),EPSI)
            CASE (0)            !  % B(1:2) = INF A(3) = -INF %OK
               VAL = TVNL( -A(1), -A(2), B(3),R(1),-R(2),-R(3),EPSI )
            CASE (:-1)
               VAL = BVNMVN(A,B,INFIN,R(1))
            END SELECT
         CASE (0) ! A(2) = -INF
            SELECT CASE ( INFIN(3))
            CASE (2)            ! B(1) = INF , A(2) = -INF %OK
               VAL = TVNL( -A(1), B(2), B(3),-R(1),R(2),-R(3),EPSI )
     &              - TVNL( -A(1), B(2),A(3),-R(1),R(2),-R(3),EPSI )
            CASE (1)            ! B(1) = INF , A(2) = -INF  B(3) = INF % OK
               VAL = TVNL( -A(1), B(2), -A(3),-R(1),-R(2),R(3),EPSI)
            CASE (0)            !% B(1) = INF , A(2:3) = -INF 
               VAL = TVNL( -A(1), B(2), B(3),-R(1),-R(2),R(3),EPSI )
            CASE (:-1)
               VAL = BVNMVN(A,B,INFIN,R(1))
            END SELECT
         CASE DEFAULT
            VAL = BVNMVN(A(1:3:2),B(1:3:2),INFIN(1:3:2),R(2))
         END SELECT
      ELSE IF ( INFIN(1) .EQ. 0 ) THEN
         SELECT CASE (INFIN(2))
         CASE (2)
            SELECT CASE (INFIN(3))
            CASE (2:)            ! A(1) = -INF %OK
               VAL =  TVNL( B(1), B(2), B(3),R(1),R(2),R(3),EPSI ) 
     &              - TVNL( B(1), B(2), A(3),R(1),R(2),R(3),EPSI) 
     &              - TVNL( B(1), A(2), B(3),R(1),R(2),R(3),EPSI )
     &              + TVNL( B(1), A(2), A(3),R(1),R(2),R(3),EPSI )
            CASE (1) ! % A(1) = -INF , B(3) = INF %OK
               VAL =  TVNL( B(1), B(2), -A(3),R(1),-R(2),-R(3),EPSI ) 
     $              - TVNL( B(1), A(2), -A(3),R(1),-R(2),-R(3),EPSI )
            CASE (0)            ! A(1) = -INF , A(3) = -INF %OK
               VAL =  TVNL( B(1), B(2), B(3),R(1),R(2),R(3),EPSI ) 
     &              - TVNL( B(1), A(2), B(3),R(1),R(2),R(3),EPSI )
            CASE DEFAULT
               VAL = BVNMVN(A,B,INFIN,R(1));
            END SELECT
         CASE (1)               ! B(2) = INF
            SELECT CASE (INFIN(3))
            CASE (2:)            ! A(1) = -INF B(2) = INF %OK
               VAL =  TVNL( B(1), -A(2), B(3),-R(1),R(2),-R(3),EPSI)
     $              - TVNL( B(1), -A(2), A(3),-R(1),R(2),-R(3),EPSI)
            CASE (1)            ! A(1) = -INF B(2) = INF  B(3) = INF %OK
               VAL =  TVNL( B(1), -A(2), -A(3),-R(1),-R(2),R(3),EPSI)
            CASE (0)            ! % A(1) = -INF B(2) = INF  A(3) = -INF %OK
               VAL = TVNL(B(1), -A(2), B(3),-R(1),R(2),-R(3),EPSI)
            CASE DEFAULT
               VAL = BVNMVN(A,B,INFIN,R(1))
            END SELECT
         CASE (0)               ! A(2) = -INF
            SELECT CASE (INFIN(3))
            CASE (2:)            ! %  A(1:2) = -INF
               VAL =  TVNL( B(1), B(2), B(3),R(1),R(2),R(3),EPSI) 
     $              - TVNL( B(1), B(2), A(3),R(1),R(2),R(3),EPSI)
            CASE (1)            ! A(1:2) = -INF B(3) = INF
               VAL =  TVNL( B(1), B(2), -A(3),R(1),-R(2),-R(3),EPSI)
            CASE (0)            !  % A(1:3) = -INF 
               VAL =  TVNL( B(1), B(2), B(3),R(1),R(2),R(3),EPSI )
            CASE DEFAULT
               VAL = BVNMVN(A,B,INFIN,R(1));
            END  SELECT
         CASE DEFAULT
            VAL = BVNMVN(A(1:3:2),B(1:3:2),INFIN(1:3:2),R(2))
         END SELECT
      ELSE
         VAL = BVNMVN(A(2:3),B(2:3),INFIN(2:3),R(3))
      END IF
      CONTAINS
      DOUBLE PRECISION FUNCTION TVNL(H1,H2,H3, R12,R13,R23, EPSI )
      !Returns Trivariate Normal CDF
      DOUBLE PRECISION, INTENT(IN) :: R12,R13,R23
      DOUBLE PRECISION, INTENT(IN) :: H1,H2,H3, EPSI
!     Locals
      INTEGER, PARAMETER :: NU = 0
      DOUBLE PRECISION,DIMENSION(3) :: H,R
      H(:) = (/ H1, H2, H3 /)
      R(:) = (/ R12, R13, R23 /)
      TVNL = TVTL(NU,H,R,EPSI)
      END FUNCTION TVNL
      END  FUNCTION TVNMVN
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
      DOUBLE PRECISION :: E
      SELECT CASE (INFIN(1))
      CASE (2:)
         SELECT CASE ( INFIN(2) )
         CASE (2:)
            VAL =  BVU ( LOWER(1), LOWER(2), CORREL )
     &           - BVU ( UPPER(1), LOWER(2), CORREL )
     &           - BVU ( LOWER(1), UPPER(2), CORREL )
     &           + BVU ( UPPER(1), UPPER(2), CORREL )

         CASE (1) 
            VAL =  BVU ( LOWER(1), LOWER(2), CORREL )
     &           - BVU ( UPPER(1), LOWER(2), CORREL )
         CASE (0) 
            VAL =  BVU ( -UPPER(1), -UPPER(2), CORREL )
     &           - BVU ( -LOWER(1), -UPPER(2), CORREL )
         CASE DEFAULT
             CALL MVNLIMITS(LOWER(1),UPPER(1),INFIN(1),E,VAL)  
         END SELECT
      CASE (1)
         SELECT CASE ( INFIN(2))
         CASE ( 2: )
            VAL =  BVU ( LOWER(1), LOWER(2), CORREL )
     &           - BVU ( LOWER(1), UPPER(2), CORREL )
         CASE (1) 
            VAL =  BVU ( LOWER(1), LOWER(2), CORREL )
         CASE (0) 
            VAL =  BVU ( LOWER(1), -UPPER(2), -CORREL )
         CASE DEFAULT
            CALL MVNLIMITS(LOWER(2),UPPER(2),INFIN(2),E,VAL)
         END SELECT  
      CASE (0)
         SELECT CASE ( INFIN(2))
         CASE ( 2: ) 
            VAL =  BVU ( -UPPER(1), -UPPER(2), CORREL )
     &           - BVU ( -UPPER(1), -LOWER(2), CORREL )
         CASE ( 1 )
            VAL =  BVU ( -UPPER(1), LOWER(2), -CORREL )
         CASE (0)
            VAL =  BVU ( -UPPER(1), -UPPER(2), CORREL )
         CASE DEFAULT
            CALL MVNLIMITS(LOWER(1),UPPER(1),INFIN(1),E,VAL) 
         END SELECT
      CASE DEFAULT !ELSE  !INFIN(1)<0
         CALL MVNLIMITS(LOWER(2),UPPER(2),INFIN(2),E,VAL)
      END SELECT
      END  FUNCTION BVNMVN
      FUNCTION BVU( SH, SK, R ) RESULT (VAL)
!      USE GLOBALDATA, ONLY: XMAX
      IMPLICIT NONE
*
!     A function for computing bivariate normal probabilities.
!
!       Yihong Ge
!       Department of Computer Science and Electrical Engineering
!       Washington State University
!       Pullman, WA 99164-2752
!     and
!       Alan Genz
!       Department of Mathematics
!       Washington State University
!       Pullman, WA 99164-3113
!       Email : alangenz@wsu.edu
!
!    This function is based on the method described by 
!        Drezner, Z and G.O. Wesolowsky, (1989),
!        On the computation of the bivariate normal integral,
!        Journal of Statist. Comput. Simul. 35, pp. 101-107,
!    with major modifications for double precision, and for |R| close to 1.
!
! BVU - calculate the probability that X > SH and Y > SK. 
!       (to accuracy of 1e-16?)
!
! Parameters
!
!   SH  REAL, lower integration limit
!   SK  REAL, lower integration limit
!   R   REAL, correlation coefficient
!
!   LG  INTEGER, number of Gauss Rule Points and Weights
!
! Revised pab added check on XMAX
      DOUBLE PRECISION, INTENT(IN) :: SH, SK, R
      DOUBLE PRECISION  :: VAL
! Local variables
      DOUBLE PRECISION :: ZERO,ONE,TWO,FOUR
      DOUBLE PRECISION :: TWOPI,SQTWOPI ,TWOPI1,FOURPI1
      DOUBLE PRECISION :: HALF,ONETHIRD,ONEEIGHT,ONESIXTEEN
      DOUBLE PRECISION :: TWELVE, EXPMIN, XMAX
      INTEGER :: I, LG, NG
      PARAMETER ( ZERO = 0.D0,ONE=1.0D0,TWO=2.0D0,HALF=0.5D0) 
      PARAMETER (FOUR = 4.0D0, TWELVE = 12.0D0)
      PARAMETER (EXPMIN =  -100.0D0)
      PARAMETER (ONESIXTEEN = 0.0625D0) !1/16
      PARAMETER (ONEEIGHT = 0.125D0 )  !1/8
      PARAMETER (ONETHIRD = 0.3333333333333333333333D0)
      PARAMETER (TWOPI   = 6.283185307179586D0 ) 
      PARAMETER (TWOPI1  = 0.15915494309190D0 )  !1/(2*pi)
      PARAMETER (FOURPI1 = 0.0795774715459476D0 ) !/1/(4*pi)
      PARAMETER (SQTWOPI = 2.50662827463100D0) ! SQRT(2*pi)
      PARAMETER (XMAX    = 8.3D0)
      DOUBLE PRECISION, DIMENSION(10,3) :: X, W
      DOUBLE PRECISION :: AS, A, B, C, D, RS, XS
      DOUBLE PRECISION :: SN, ASR, H, K, BS, HS, HK
!     Gauss Legendre Points and Weights, N =  6
      DATA  ( W(I,1), X(I,1), I = 1,3) /
     *  0.1713244923791705D+00,-0.9324695142031522D+00,
     *  0.3607615730481384D+00,-0.6612093864662647D+00,
     *  0.4679139345726904D+00,-0.2386191860831970D+00/
!     Gauss Legendre Points and Weights, N = 12
      DATA ( W(I,2), X(I,2), I = 1,6) /
     *  0.4717533638651177D-01,-0.9815606342467191D+00,
     *  0.1069393259953183D+00,-0.9041172563704750D+00,
     *  0.1600783285433464D+00,-0.7699026741943050D+00,
     *  0.2031674267230659D+00,-0.5873179542866171D+00,
     *  0.2334925365383547D+00,-0.3678314989981802D+00,
     *  0.2491470458134029D+00,-0.1252334085114692D+00/
!     Gauss Legendre Points and Weights, N = 20
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
      VAL = ZERO
      HK = MIN(SH,SK)
      IF ( HK  < -XMAX) THEN    ! pab 24.05.2003
         VAL = FI(-MAX(SH,SK))
         RETURN
      ELSE IF ( XMAX  < MAX(SH,SK)) THEN
         RETURN
      ENDIF 
      IF ( ABS(R)  <  0.3D0 ) THEN
         NG = 1
         LG = 3
      ELSE IF ( ABS(R)  <  0.75D0 ) THEN
         NG = 2
         LG = 6
      ELSE 
         NG = 3
         LG = 10
      ENDIF
      H = SH
      K = SK 
      HK = H*K
      
      IF ( ABS(R)  <  0.925D0 ) THEN
         IF (ABS(R) .GT. ZERO ) THEN
            HS = ( H*H + K*K )*HALF
            ASR = ASIN(R)
            DO I = 1, LG
               SN  = SIN(ASR*(ONE + X(I,NG))*HALF)
               VAL = VAL + W(I,NG)*EXP( ( SN*HK - HS )/( ONE - SN*SN ) )
               SN  = SIN(ASR*(ONE - X(I,NG))*HALF)
               VAL = VAL + W(I,NG)*EXP( ( SN*HK - HS )/( ONE - SN*SN ) )
            END DO
            VAL = VAL*ASR*FOURPI1
         ENDIF
         VAL = VAL + FI(-H)*FI(-K) 
      ELSE
         IF ( R  <  ZERO ) THEN
            K  = -K
            HK = -HK
         ENDIF
         IF ( ABS(R)  <  ONE ) THEN
            AS  = ( ONE - R )*( ONE + R )
            A   = SQRT(AS)
            B   = ABS( H - K ) !**2
            BS  = B * B
            C   = ( FOUR - HK ) * ONEEIGHT !/8D0 
            D   = ( TWELVE - HK ) * ONESIXTEEN !/16D0
            ASR =  -(BS/AS + HK)*HALF
            IF (ASR.GT.EXPMIN) THEN
               VAL = A*EXP( ASR ) *
     &              ( ONE - C*(BS - AS)*(ONE - D*BS*0.2D0)*ONETHIRD + 
     &              C*D*AS*AS*0.2D0 )
            ENDIF
            IF ( HK .GT. EXPMIN ) THEN
               VAL = VAL - EXP(-HK*HALF)*SQTWOPI*FI(-B/A)*B
     +              *( ONE - C*BS*( ONE - D*BS*0.2D0 )*ONETHIRD ) 
            ENDIF
            A = A * HALF
            DO I = 1, LG
               XS  = ( A * (ONE + X(I,NG)) ) !**2
               XS  = XS * XS
               RS  = SQRT( ONE - XS )
               ASR = -(BS / XS + HK) * HALF
               IF (ASR.GT.EXPMIN) THEN
                  VAL = VAL + A*W(I,NG)*EXP( ASR )
     &                 * ( EXP( - HALF*HK*( ONE - RS )/( ONE + RS ) )/RS
     $                 -( ONE + C*XS*( ONE + D*XS ) ) )
               ENDIF
               XS  = ( A * (ONE - X(I,NG)) ) !**2
               XS  = XS * XS
               RS  = SQRT( ONE - XS )
               ASR = -(BS / XS + HK) * HALF
               IF (ASR.GT.EXPMIN) THEN
                  VAL = VAL + A*W(I,NG)*EXP( ASR )
     &                 *( EXP( - HALF*HK*( ONE - RS )/( ONE + RS ) )/RS-
     $                 ( ONE + C*XS*( ONE + D*XS ) ) )
               ENDIF
            END DO
            VAL = -VAL*TWOPI1
         ENDIF
         IF ( R .GT. ZERO ) THEN
            VAL =  VAL + FI( -MAX( H, K ) )
         ELSE
            VAL = -VAL 
            IF ( H  <  K )  VAL = VAL +  FI(K)-FI(H) 
         ENDIF
      ENDIF
      RETURN
      END FUNCTION BVU
      DOUBLE PRECISION FUNCTION STUDNT( NU, T )
      IMPLICIT NONE
!
!     Student t Distribution Function
!
!                       T
!         STUDNT = C   I  ( 1 + y*y/NU )**( -(NU+1)/2 ) dy
!                   NU -INF
!
      INTEGER, INTENT(IN) :: NU
      DOUBLE PRECISION, INTENT(IN) :: T
!     Locals
      INTEGER :: J
      DOUBLE PRECISION :: ZRO, ONE
      PARAMETER ( ZRO = 0.0D0, ONE = 1.0D0 )
      DOUBLE PRECISION, PARAMETER :: PI = 3.14159265358979D0
      DOUBLE PRECISION :: CSSTHE, SNTHE, POLYN, TT, TS, RN
    
      IF ( NU  <  1 ) THEN
         STUDNT = FI( T )
      ELSE IF ( NU .EQ. 1 ) THEN
         STUDNT = ( ONE + 2.0D0*ATAN(T)/PI )*0.5D0
      ELSE IF ( NU .EQ. 2 ) THEN
         STUDNT = ( ONE + T/SQRT( 2.0D0 + T*T ))*0.5D0
      ELSE 
         RN = NU ! convert to double
         TT = T * T
         CSSTHE = ONE/( ONE + TT/RN )
         POLYN = 1
         DO J = NU-2, 2, -2
            POLYN = ONE + ( J - 1 )*CSSTHE*POLYN/J
         END DO
         
         IF ( MOD( NU, 2 ) .EQ. 1 ) THEN
            TS = T/SQRT(RN)
            STUDNT = ( ONE + 2.0D0*( ATAN(TS) + 
     &           TS*CSSTHE*POLYN )/PI )*0.5D0
         ELSE
            SNTHE = T/SQRT( RN + TT )
            STUDNT = ( ONE + SNTHE*POLYN )*0.5D0
         END IF
         STUDNT = MAX( ZRO, MIN( STUDNT, ONE ) )
      ENDIF
      END FUNCTION STUDNT
      DOUBLE PRECISION FUNCTION BVTL( NU, DH, DK, R )
      IMPLICIT NONE
!*
!*     A function for computing bivariate t probabilities.
!*
!*       Alan Genz
!*       Department of Mathematics
!*       Washington State University
!*       Pullman, WA 99164-3113
!*       Email : alangenz@wsu.edu
!*
!*    This function is based on the method described by 
!*        Dunnett, C.W. and M. Sobel, (1954),
!*        A bivariate generalization of Student's t-distribution
!*        with tables for certain special cases,
!*        Biometrika 41, pp. 153-169.
!*
!* BVTL - calculate the probability that X < DH and Y < DK. 
!*
!* parameters
!*
!*   NU number of degrees of freedom (NOTE: NU = 0 gives bivariate normal prb)
!*   DH 1st lower integration limit
!*   DK 2nd lower integration limit
!*   R   correlation coefficient
!*
      INTEGER, INTENT(IN) ::NU
      DOUBLE PRECISION, INTENT(IN) :: DH, DK, R
!     Locals
      INTEGER :: J, HS, KS
      DOUBLE PRECISION ::  ORS, HRK, KRH, HRK2, KRH2, BVT
      DOUBLE PRECISION ::  DH2, DK2, SNU ,DNU, DHDK
!, BVND, STUDNT
      DOUBLE PRECISION :: GMPH, GMPK, XNKH, XNHK, QHRK, HKN, HPK, HKRN
      DOUBLE PRECISION :: BTNCKH, BTNCHK, BTPDKH, BTPDHK
      DOUBLE PRECISION :: ZERO, ONE, EPS, PI,TPI
      PARAMETER ( ZERO = 0.0D0, ONE = 1.0D0, EPS = 1.0D-15 )
      PARAMETER (PI =  3.14159265358979D0, TPI = 6.28318530717959D0)
      IF ( NU  <  1 ) THEN
         BVTL = BVU( -DH, -DK, R )
      ELSE IF ( ONE - R .LE. EPS .OR. 1.0D+16<MAX(ABS(DH),ABS(DK))) THEN
            BVTL = STUDNT( NU, MIN( DH, DK ) )
      ELSE IF ( R + ONE  .LE. EPS ) THEN
         IF ( DH .GT. -DK )  THEN
            BVTL = STUDNT( NU, DH ) - STUDNT( NU, -DK )
         ELSE
            BVTL = ZERO
         END IF
      ELSE 
         !PI = ACOS(-ONE)
         !TPI = 2*PI
         DNU = NU           ! convert to double
         SNU = SQRT(DNU)
         ORS = ONE - R * R  
         HRK = DH - R * DK  
         KRH = DK - R * DH  
         DK2 = DK * DK
         DH2 = DH * DH
         IF ( ABS(HRK) + ORS .GT. ZERO ) THEN
            XNHK = HRK**2/( HRK**2 + ORS*( DNU + DK2 ) ) 
            XNKH = KRH**2/( KRH**2 + ORS*( DNU + DH2 ) ) 
         ELSE
            XNHK = ZERO  
            XNKH = ZERO  
         END IF
         HS = SIGN( ONE, HRK) !DH - R*DK )  
         KS = SIGN( ONE, KRH) !DK - R*DH ) 
         IF ( MOD( NU, 2 ) .EQ. 0 ) THEN
            BVT = ATAN2( SQRT(ORS), -R )/TPI 
            GMPH = DH/SQRT( 16.0D0*( DNU + DH2 ) )  
            GMPK = DK/SQRT( 16.0D0*( DNU + DK2 ) )  
            BTNCKH = 2*ATAN2( SQRT( XNKH ), SQRT( ONE - XNKH ) )/PI  
            BTPDKH = 2*SQRT( XNKH*( ONE - XNKH ) )/PI 
            BTNCHK = 2*ATAN2( SQRT( XNHK ), SQRT( ONE - XNHK ) )/PI  
            BTPDHK = 2*SQRT( XNHK*( ONE - XNHK ) )/PI 
            DO J = 1, NU/2
               BVT = BVT + GMPH*( ONE + KS*BTNCKH ) 
               BVT = BVT + GMPK*( ONE + HS*BTNCHK ) 
               BTNCKH = BTNCKH + BTPDKH  
               BTPDKH = 2*J*BTPDKH*( ONE - XNKH )/( 2*J + 1 )  
               BTNCHK = BTNCHK + BTPDHK  
               BTPDHK = 2*J*BTPDHK*( ONE - XNHK )/( 2*J + 1 )  
               GMPH = GMPH*( 2*J - 1 )/( 2*J*( ONE + DH2/DNU ) ) 
               GMPK = GMPK*( 2*J - 1 )/( 2*J*( ONE + DK2/DNU ) ) 
            END DO
         ELSE  ! NU is ODD
            DHDK = DH*DK
            QHRK = SQRT( DH2 + DK2 - 2.0D0*R*DHDK + DNU*ORS )  
            HKRN = DHDK + R*DNU  
            HKN = DHDK - DNU  
            HPK = DH + DK 
            BVT = ATAN2( -SNU*( HKN*QHRK + HPK*HKRN ),
     &                          HKN*HKRN-DNU*HPK*QHRK )/TPI
            IF ( BVT  <  -EPS ) BVT = BVT + ONE
            GMPH = DH/( TPI*SNU*( ONE + DH2/DNU ) )  
            GMPK = DK/( TPI*SNU*( ONE + DK2/DNU ) )  
            BTNCKH = SQRT( XNKH )  
            BTPDKH = BTNCKH 
            BTNCHK = SQRT( XNHK )  
            BTPDHK = BTNCHK  
            DO J = 1, ( NU - 1 )/2
               BVT = BVT + GMPH*( ONE + KS*BTNCKH ) 
               BVT = BVT + GMPK*( ONE + HS*BTNCHK ) 
               BTPDKH = ( 2*J - 1 )*BTPDKH*( ONE - XNKH )/( 2*J )  
               BTNCKH = BTNCKH + BTPDKH  
               BTPDHK = ( 2*J - 1 )*BTPDHK*( ONE - XNHK )/( 2*J )  
               BTNCHK = BTNCHK + BTPDHK  
               GMPH = 2*J*GMPH/( ( 2*J + 1 )*( ONE + DH2/DNU ) ) 
               GMPK = 2*J*GMPK/( ( 2*J + 1 )*( ONE + DK2/DNU ) ) 
            END DO
         END IF
         BVTL = BVT 
      END IF
      END FUNCTION BVTL
     
      DOUBLE PRECISION FUNCTION TVTL( NU1, H, R, EPSI )
      USE TRIVARIATEVAR !Block to transfer variables to TVTMFN
      IMPLICIT NONE
!    
!     A function for computing trivariate normal and t-probabilities.
!     This function uses algorithms developed from the ideas 
!     described in the papers:
!       R.L. Plackett, Biometrika 41(1954), pp. 351-360.
!       Z. Drezner, Math. Comp. 62(1994), pp. 289-294.
!     with adaptive integration from (0,0,1) to (0,0,r23) to R. 
!
!      Calculate the probability that X(I) < H(I), for I = 1,2,3     
!    NU   INTEGER degrees of freedom; use NU = 0 for normal cases.
!    H    REAL array of uppoer limits for probability distribution 
!    R    REAL array of three correlation coefficients, R should 
!         contain the lower left portion of the correlation matrix r. 
!         R should contains the values r21, r31, r23 in that order.
!   EPSI  REAL required absolute accuracy; maximum accuracy for most
!          computations is approximately 1D-14
! 
!    The software is based on work described in the paper
!     "Numerical Computation of Rectangular Bivariate and Trivariate
!      Normal and t Probabilities", by the code author:
!
!       Alan Genz
!       Department of Mathematics
!       Washington State University
!       Pullman, WA 99164-3113
!       Email : alangenz@wsu.edu
!
!      EXTERNAL TVTMFN
      INTEGER,                       INTENT(IN) :: NU1
      DOUBLE PRECISION,DIMENSION(:), INTENT(IN) :: H, R
      DOUBLE PRECISION,              INTENT(IN) ::  EPSI
!Locals
      DOUBLE PRECISION, DIMENSION(3) :: ZROS, HS
      DOUBLE PRECISION, DIMENSION(1) :: HMAX,HMIN
      DOUBLE PRECISION :: R12, R13,  TVT
      DOUBLE PRECISION :: ONE, ZERO, EPS,  PT
!      DOUBLE PRECISION RUA, RUB, AR, RUC,
!        BVTL, PHID, ADONET
      PARAMETER ( ZERO = 0.0D0, ONE = 1.0D0 )
      PARAMETER ( PT = 1.57079632679489661923132169163975D0 ) !pi/2
      
!      COMMON /TVTMBK/ H1, H2, H3, R23, RUA, RUB, AR, RUC, NU
      EPS = MAX( 1.0D-13, EPSI )
!      PT = ASIN(ONE)
      NU = NU1
      H1 = H(1)
      H2 = H(2)
      H3 = H(3)
      R12 = R(1)
      R13 = R(2)
      R23 = R(3)
!
!     Sort R's and check for special cases
!
      IF ( ABS(R12) .GT. ABS(R13) ) THEN
         H2 = H3
         H3 = H(2)
         R12 = R13
         R13 = R(1)
      END IF
      IF ( ABS(R13) .GT. ABS(R23) ) THEN
         H1 = H2
         H2 = H(1)
         R23 = R13
         R13 = R(3)
      END IF
      TVT = 0
      IF ( ABS(H1) + ABS(H2) + ABS(H3)  <  EPS ) THEN 
         TVT = (ONE + (ASIN(R12) + ASIN(R13) + ASIN(R23))/PT )*0.125D0
      ELSE IF ( NU  <  1 .AND. ABS(R12) + ABS(R13)  <  EPS ) THEN 
         TVT = FI(H1)*BVTL( NU, H2, H3, R23 )
      ELSE IF ( NU  <  1 .AND. ABS(R13) + ABS(R23)  <  EPS ) THEN 
         TVT = FI(H3)*BVTL( NU, H1, H2, R12 )
      ELSE IF ( NU  <  1 .AND. ABS(R12) + ABS(R23)  <  EPS ) THEN 
         TVT = FI(H2)*BVTL( NU, H1, H3, R13 )
      ELSE IF ( ONE - R23  <  EPS ) THEN
         TVT = BVTL( NU, H1, MIN( H2, H3 ), R12 )
      ELSE IF ( R23 + ONE  <  EPS ) THEN
         IF  ( H2 .GT. -H3 ) 
     &        TVT = BVTL( NU, H1, H2, R12 ) - BVTL( NU, H1, -H3, R12 )
      ELSE
!
!        Compute singular TVT value
!
         IF ( NU  <  1 ) THEN
            TVT = BVTL( NU, H2, H3, R23 )*FI(H1)
         ELSE IF ( R23 .GE. ZERO ) THEN
            TVT = BVTL( NU, H1, MIN( H2, H3 ), ZERO )
         ELSE IF ( H2 .GT. -H3 ) THEN 
            TVT = BVTL( NU, H1, H2, ZERO ) - BVTL( NU, H1, -H3, ZERO )
         END IF
!
!        Use numerical integration to compute probability
!
!
         RUA = ASIN( R12 )
         RUB = ASIN( R13 )
         AR  = ASIN( R23)
         RUC = SIGN( PT, AR ) - AR
         TVT = TVT + ADONET( TVTMFN, ZERO, ONE, EPS )/( 4.D0*PT )
      END IF
      TVTL = MAX( ZERO, MIN( TVT, ONE ) ) 

      CONTAINS
!
      DOUBLE PRECISION FUNCTION TVTMFN( X )
      USE TRIVARIATEVAR
      IMPLICIT NONE
!
!     Computes Plackett formula integrands
!
      DOUBLE PRECISION, INTENT(IN) :: X
! Locals
      DOUBLE PRECISION R12, RR2, R13, RR3, R, RR, ZRO !, PNTGND
! Parameters transfeered from TRIVARIATEVAR
!      INTEGER :: NU
!      DOUBLE PRECISION :: H1, H2, H3, R23, RUA, RUB, AR, RUC
      PARAMETER ( ZRO = 0.0D0 )    
!      COMMON /TVTMBK/ H1, H2, H3, R23, RUA, RUB, AR, RUC, NU
      TVTMFN = 0.0D0
      CALL SINCS( RUA*X, R12, RR2 )
      CALL SINCS( RUB*X, R13, RR3 )
      IF ( ABS(RUA) .GT. ZRO ) 
     &     TVTMFN = TVTMFN + RUA*PNTGND( NU, H1,H2,H3, R13,R23,R12,RR2 )
      IF ( ABS(RUB) .GT. ZRO )
     &     TVTMFN = TVTMFN + RUB*PNTGND( NU, H1,H3,H2, R12,R23,R13,RR3 )
      IF ( NU .GT. 0 ) THEN
         CALL SINCS( AR + RUC*X, R, RR )
         TVTMFN = TVTMFN - RUC*PNTGND( NU, H2, H3, H1, ZRO, ZRO, R, RR )
      END IF
      END FUNCTION TVTMFN     
!
      SUBROUTINE SINCS( X, SX, CS )
!
!     Computes SIN(X), COS(X)^2, with series approx. for |X| near PI/2
!
      DOUBLE PRECISION, INTENT(IN) :: X
      DOUBLE PRECISION, INTENT(OUT) :: SX, CS
!Locals      
      DOUBLE PRECISION :: EE, PT, KS, KC, ONE, SMALL, HALF, ONETHIRD
      PARAMETER (ONE = 1.0D0, SMALL = 5.0D-5, HALF = 0.5D0 )
      PARAMETER ( PT = 1.57079632679489661923132169163975D0 )
      PARAMETER ( KS = 0.0833333333333333333333333333333D0) !1/12
      PARAMETER ( KC = 0.1333333333333333333333333333333D0) !2/15
      PARAMETER ( ONETHIRD = 0.33333333333333333333333333333333D0) !1/3
      EE = ( PT - ABS(X) )
      EE = EE * EE
      IF ( EE  <  SMALL ) THEN
         SX = SIGN( ONE - EE*( ONE - EE*KS )*HALF, X )
         CS = EE *( ONE - EE*( ONE - EE*KC )*ONETHIRD)
      ELSE
         SX = SIN(X)
         CS = ONE - SX*SX
      END IF
      END SUBROUTINE SINCS
!
      DOUBLE PRECISION FUNCTION PNTGND( NU, BA, BB, BC, RA, RB, R, RR )
      IMPLICIT NONE
!
!     Computes Plackett formula integrand
!
      INTEGER, INTENT(IN) :: NU
      DOUBLE PRECISION, INTENT(IN) :: BA, BB, BC, RA, RB, R, RR
! Locals      
      DOUBLE PRECISION :: DT, FT, BT,RAB2, BARB, rNU!, PHID, STUDNT
      PNTGND = 0.0D0
      FT   = ( RA - RB )
      RAB2 = FT*FT
      DT   = RR*( RR -  RAB2 - 2.0D0*RA*RB*( 1.D0 - R ) )
      IF ( DT .GT. 0.0D0 ) THEN
         BT   = ( BC*RR + BA*( R*RB - RA ) + BB*( R*RA -RB ) )/SQRT(DT) 
         BARB = ( BA - R*BB )
         FT   = ( BARB * BARB ) / RR + BB * BB
         IF ( NU  <  1 ) THEN
            IF ( BT .GT. -10.0D0 .AND. FT  <  100.0D0 ) THEN
               PNTGND = EXP( -FT * 0.5D0 ) * FI( BT )
!               PNTGND = EXP( -FT*0.5D0)
!               IF ( BT  <  10.0D0 ) PNTGND = PNTGND * FI(BT)
            END IF
         ELSE
            rNU = NU
            FT  = SQRT( 1.0D0 + FT/rNU )
            PNTGND = STUDNT( NU, BT/FT )/FT**NU
         END IF
      END IF
      END  FUNCTION PNTGND
!
      DOUBLE PRECISION FUNCTION ADONET( F, A, B, TOL )
      IMPLICIT NONE
!
!     One Dimensional Globally Adaptive Integration Function
!
!      EXTERNAL F
      DOUBLE PRECISION, INTENT(IN) :: A, B, TOL
      INTEGER :: NL, I, IM, IP
      PARAMETER ( NL = 100 )
      DOUBLE PRECISION, DIMENSION(NL) ::  EI, AI, BI, FI
      DOUBLE PRECISION :: FIN, ERR !, KRNRDT
      DOUBLE PRECISION, PARAMETER :: ZERO = 0.0D0
      DOUBLE PRECISION, PARAMETER :: HALF = 0.5D0
      DOUBLE PRECISION, PARAMETER :: ONE  = 1.0D0
      DOUBLE PRECISION, PARAMETER :: FOUR = 4.0D0
      INTERFACE 
         FUNCTION F(Z) RESULT (VAL)
         DOUBLE PRECISION, INTENT(IN) :: Z
         DOUBLE PRECISION :: VAL
         END FUNCTION F
      END INTERFACE
!      COMMON /ABLK/ ERR, IM
      AI(1) = A
      BI(1) = B
      ERR = ONE
      IP = 1
      IM = 1
      DO WHILE ( (ERR .GT. TOL) .AND. (IM  <  NL) ) 
         IM = IM + 1
         BI(IM) = BI(IP)
         AI(IM) = ( AI(IP) + BI(IP) ) * HALF
         BI(IP) = AI(IM)
         FI(IP) = KRNRDT( AI(IP), BI(IP), F, EI(IP) )
         FI(IM) = KRNRDT( AI(IM), BI(IM), F, EI(IM) )
         ERR = ZERO
         FIN = ZERO
         DO I = 1, IM
            IF ( EI(I) .GT. EI(IP) ) IP = I
            FIN = FIN + FI(I)
            ERR = ERR + EI(I)*EI(I)
         END DO
         ERR = FOUR * SQRT( ERR )
      END DO
      ADONET = FIN
      END FUNCTION ADONET
!
      DOUBLE PRECISION FUNCTION KRNRDT( A, B, F, ERR )
!
!     Kronrod Rule
!
      DOUBLE PRECISION A, B, ERR, T, CEN, FC, WID, RESG, RESK
!
!        The abscissae and weights are given for the interval (-1,1);
!        only positive abscissae and corresponding weights are given.
!
!        XGK    - abscissae of the 2N+1-point Kronrod rule: 
!                 XGK(2), XGK(4), ...  N-point Gauss rule abscissae; 
!                 XGK(1), XGK(3), ...  optimally added abscissae.
!        WGK    - weights of the 2N+1-point Kronrod rule.
!        WG     - weights of the N-point Gauss rule.
!
      INTEGER :: J, N
      PARAMETER ( N = 11 )
      DOUBLE PRECISION, PARAMETER :: HALF = 0.5D0
      DOUBLE PRECISION WG(0:(N+1)/2), WGK(0:N), XGK(0:N) 
      SAVE WG, WGK, XGK
      INTERFACE 
          FUNCTION F(Z) RESULT (VAL)
         DOUBLE PRECISION, INTENT(IN) :: Z
         DOUBLE PRECISION :: VAL
         END FUNCTION F
      END INTERFACE
      DATA WG( 0)/ 0.2729250867779007D+00/
      DATA WG( 1)/ 0.5566856711617449D-01/
      DATA WG( 2)/ 0.1255803694649048D+00/
      DATA WG( 3)/ 0.1862902109277352D+00/
      DATA WG( 4)/ 0.2331937645919914D+00/
      DATA WG( 5)/ 0.2628045445102478D+00/
!
      DATA XGK( 0)/ 0.0000000000000000D+00/
      DATA XGK( 1)/ 0.9963696138895427D+00/
      DATA XGK( 2)/ 0.9782286581460570D+00/
      DATA XGK( 3)/ 0.9416771085780681D+00/
      DATA XGK( 4)/ 0.8870625997680953D+00/
      DATA XGK( 5)/ 0.8160574566562211D+00/
      DATA XGK( 6)/ 0.7301520055740492D+00/
      DATA XGK( 7)/ 0.6305995201619651D+00/
      DATA XGK( 8)/ 0.5190961292068118D+00/
      DATA XGK( 9)/ 0.3979441409523776D+00/
      DATA XGK(10)/ 0.2695431559523450D+00/
      DATA XGK(11)/ 0.1361130007993617D+00/
!
      DATA WGK( 0)/ 0.1365777947111183D+00/
      DATA WGK( 1)/ 0.9765441045961290D-02/
      DATA WGK( 2)/ 0.2715655468210443D-01/
      DATA WGK( 3)/ 0.4582937856442671D-01/
      DATA WGK( 4)/ 0.6309742475037484D-01/
      DATA WGK( 5)/ 0.7866457193222764D-01/
      DATA WGK( 6)/ 0.9295309859690074D-01/
      DATA WGK( 7)/ 0.1058720744813894D+00/
      DATA WGK( 8)/ 0.1167395024610472D+00/
      DATA WGK( 9)/ 0.1251587991003195D+00/
      DATA WGK(10)/ 0.1312806842298057D+00/
      DATA WGK(11)/ 0.1351935727998845D+00/
!
!           Major variables
!
!           CEN  - mid point of the interval
!           WID  - half-length of the interval
!           RESG - result of the N-point Gauss formula
!           RESK - result of the 2N+1-point Kronrod formula
!
!           Compute the 2N+1-point Kronrod approximation to
!            the integral, and estimate the absolute error.
!
      WID = ( B - A ) * HALF
      CEN = ( B + A ) * HALF
      FC  = F(CEN)
      RESG = FC * WG(0)
      RESK = FC * WGK(0)
      DO J = 1, N
         T  = WID * XGK(J) 
         FC = F( CEN - T ) + F( CEN + T )
         RESK = RESK + WGK(J) * FC
         IF( MOD( J, 2 ) .EQ. 0 ) RESG = RESG + WG(J/2) * FC
      END DO
      KRNRDT = WID * RESK
      ERR = ABS( WID * ( RESK - RESG ) )
      END FUNCTION KRNRDT
      END FUNCTION TVTL

      FUNCTION GAUSINT (X1, X2, A, B, C, D) RESULT (value)               
!      USE GLOBALDATA,ONLY:  xCutOff
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in) :: X1,X2,A,B,C,D
      DOUBLE PRECISION             :: value
!Local variables
      DOUBLE PRECISION             :: Y1,Y2,Y3
      DOUBLE PRECISION, PARAMETER :: SQTWOPI1=3.9894228040143d-1 !=1/sqrt(2*pi)
      DOUBLE PRECISION, PARAMETER :: XMAX = 37.d0
      ! Let  X  be standardized Gaussian variable, 
      ! i.e., X=N(0,1). The function calculate the
      !  following integral E[I(X1<X<X2)(A+BX)(C+DX)
      ! where I(X1<X<X2) is an indicator function of
      ! the set {X1<X<X2}. 
      IF (X1.GE.X2) THEN                                                 
         value = 0.d0                                                    
         RETURN                                                          
      ENDIF
      IF (ABS (X1) .GT.XMAX) THEN                                           
         Y1 = 0.d0                                                         
      ELSE                                                               
         Y1 = (A * D+B * C + X1 * B * D) * EXP ( - 0.5d0 * X1 * X1)        
      ENDIF
      IF (ABS (X2) .GT.XMAX) THEN                                           
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
      IF ((B*X1 < -A).AND.(-A < B*X2))THEN
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
      IF ( INFIN  <  0 ) RETURN
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
      if (EPSL < Ak .AND. Ak < amax) THEN
         IF (ABS(p-Pa).LE.EPSL) THEN
            VAL = SIGN(MIN(Ak,xmax),-A)
            RETURN
         ENDIF
         IF (Ak < 1D-2) THEN   ! starting guess always less than 0.2 from the true value
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
            ELSEIF ((P < 1D-3.AND.A.LE.-3.5D0).OR. 
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
         IF (xk.GE.-A .AND. P < Pa) xk = -1.D-2-A

      
         IF (P < Pa) THEN
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
     
     
!     Break out of the iteration loop for the following:
!     1) The last update is very small (compared to x).
!     2) The last update is very small (compared to sqrt(eps)=crit).
!     3) There are more than 15 iterations. This should NEVER happen. 
      IF (.TRUE..OR.ABS(ak) < 1.D-2) THEN 
      ! Newton's method
      !~~~~~~~~~~~~~~~~~
      DO WHILE( ABS(H).GT.MIN(crit*ABS(xk),crit).AND.IX < IXMAX) 
         
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
         DO WHILE (ABS(H).GT.MIN(crit*ABS(xk),crit).AND.IX < IXMAX)                                
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

      MODULE RIND
      USE GLOBALCONST
      USE PRINTMOD  ! used for debugging only
      IMPLICIT NONE
      PRIVATE
      PUBLIC :: RINDD, SetConstants
      private :: preInit
      private :: initIntegrand
      private :: initfun,mvnfun,cvsrtxc,covsrt1,covsrt,rcscale,rcswap
      private :: cleanUp
      
! mInfinity = what is considered as infinite value in FI
! mFxcEpss  = if fxc is less, do not compute E(...|Xc)     
! mXcEps2   =  if any Var(Xc(j)|Xc(1),...,Xc(j-1)) <= XCEPS2 then return NAN
      double precision, parameter :: mInfinity = 8.25d0 ! 37.0d0
      double precision, parameter :: mFxcEpss = 1.0D-20 
      double precision, save      :: mXcEps2   = 2.3d-16
!     Constants defining accuracy of integration:
!     mCovEps = termination criteria for Cholesky decomposition
!     mAbsEps = requested absolute tolerance
!     mRelEps = requested relative tolerance
!     mXcutOff = truncation value to c1c2
!     mXcScale = scale factor in the exponential (in order to avoid overflow) 
!     mNc1c2  = number of times to use function c1c2, i.e.,regression
!               equation to restrict integration area.
!     mNIT    = maximum number of Xt variables to integrate
!     mMethod = integration method:
!            1 Integrate all by SADAPT if NDIM<9 otherwise by KRBVRC (default)
!            2 Integrate all by SADAPT by Genz (1992) (Fast and reliable)
!            3 Integrate all by KRBVRC by Genz (1998) (Fast and reliable)
!            4 Integrate all by KROBOV by Genz (1992) (Fast and reliable)
!            5 Integrate all by RCRUDE by Genz (1992) (Reliable)
!            6 integrate all by SOBNIED by Hong and Hickernell
      double precision, save :: mCovEps  = 1.0d-10
      double precision, save :: mAbsEps  = 0.01d0
      double precision, save :: mRelEps  = 0.01d0
      double precision, save :: mXcutOff = 5.d0
      double precision, save :: mXcScale = 0.0d0
      integer, save :: mNc1c2  = 2
      integer, save :: mNIT    = 1000
      integer, save :: mMaxPts = 40000
      integer, save :: mMinPts = 0
      integer, save :: mMethod = 3
      
      
!     Integrand variables:
!     mBIG    = Cholesky Factor/Covariance matrix: 
!               Upper triangular part is the cholesky factor
!               Lower triangular part contains the conditional 
!               standarddeviations 
!               (mBIG2 is only used if mNx>1)
!     mCDI    = Cholesky DIagonal elements
!     mA,mB   = Integration limits
!     mINFI   = integrationi limit flags
!     mCm     = conditional mean
!     mINFIXt, 
!     mINFIXd = # redundant variables of Xt and Xd,
!              respectively      
!     mIndex1, 
!     mIndex2 = indices to the variables original place.   Size  Ntdc
!     xedni   = indices to the variables new place.        Size  Ntdc
!     mNt     = # Xt variables
!     mNd     = # Xd variables
!     mNc     = # Xc variables
!     mNtd    = mNt + mNd
!     mNtdc   = mNt + mNd + mNc
!     mNx     = # different integration limits
      
      double precision,allocatable, dimension(:,:) :: mBIG,mBIG2
      double precision,allocatable, dimension(:)   :: mA,mB,mCDI,mCm 
      INTEGER, DIMENSION(:), ALLOCATABLE :: mInfi,mIndex1,mIndex2,mXedni     
      INTEGER,SAVE :: mNt,mNd,mNc,mNtdc, mNtd, mNx ! Size information
      INTEGER,SAVE :: mInfiXt,mInfiXd 
      logical,save :: mInitIntegrandCalled = .FALSE.

      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: mCDIXd, mCmXd
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: mXd, mXc, mY
      double precision, save :: mSmall = 2.3d-16
      
!     variables set in initfun and used in mvnfun:
      INTEGER, PRIVATE :: mI0,mNdleftN0
      DOUBLE PRECISION, PRIVATE :: mE1,mD1, mVAL0
      
      contains
      subroutine setConstants(method,xcscale,abseps,releps,coveps,
     &     maxpts,minpts,nit,xcutoff,Nc1c2)
      double precision, optional, intent(in) :: xcscale,abseps,releps
     $     ,coveps, xcutoff
      integer, optional,intent(in) :: method,nit,maxpts,minpts,Nc1c2
      double precision, parameter :: one = 1.0d0
      mSmall = spacing(one)
      if (present(method))  mMethod = method
      if (present(xcscale)) mXcScale = xcscale
      if (present(abseps)) mAbsEps = max(abseps,mSmall)
      if (present(releps)) mRelEps = max(releps,0.0d0)
      if (present(coveps)) mCovEps = max(coveps,1d-12)
      if (present(maxpts)) mMaxPts = maxpts
      if (present(minpts)) mMinPts = minpts
      if (present(nit))    mNit    = nit
      if (present(xcutOff)) mXcutOff = xCutOff
      if (present(Nc1c2))  mNc1c2   = max(Nc1c2,1)
      end subroutine setConstants

      subroutine preInit(BIG,Xc,Nt,inform)
      double precision,dimension(:,:), intent(in) :: BIG
      double precision,dimension(:,:), intent(in) :: Xc
      integer, intent(in)  :: Nt
      integer, intent(out) :: inform
!     Local variables
      integer :: I,J
      inform = 0
      mInitIntegrandCalled = .FALSE.
!     Find the size information
!~~~~~~~~~~~~~~~~~~~~~~~~~~
      mNt   = Nt
      mNc   = SIZE( Xc, dim = 1 )
      mNx   = MAX( SIZE( Xc, dim = 2), 1 )
      mNtdc = SIZE( BIG, dim = 1 )
      ! make sure it does not exceed Ntdc-Nc
      IF (mNt+mNc.GT.mNtdc) mNt = mNtdc - mNc  
      mNd  = mNtdc-mNt-mNc
      mNtd = mNt+mNd
      IF (mNd < 0) THEN
!         PRINT *,'RINDD Nt,Nd,Nc,Ntdc=',Nt,Nd,Nc,Ntdc
         ! Size information inconsistent
         inform = 3
         return
      ENDIF
      
 !     PRINT *,'Nt Nd Nc Ntd Ntdc,',Nt, Nd, Nc, Ntd, Ntdc
     
! ALLOCATION
!~~~~~~~~~~~~
      IF (mNd>0) THEN
         ALLOCATE(mXd(mNd),mCmXd(mNd),mCDIXd(mNd))
         mCmXd(:)  = gZERO
         mCDIXd(:) = gZERO
         mxd(:)    = gZERO
      END IF
	
      ALLOCATE(mBIG(mNtdc,mNtdc),mCm(mNtdc),mY(mNtd))
      ALLOCATE(mIndex1(mNtdc),mA(mNtd),mB(mNtd),mINFI(mNtd),mXc(mNc)) 
      ALLOCATE(mCDI(mNtd),mXedni(mNtdc),mIndex2(mNtdc))
 
!     Initialization
!~~~~~~~~~~~~~~~~~~~~~
!     Copy upper triangular of input matrix, only.
      do i = 1,mNtdc
         mBIG(1:i,i)    = BIG(1:i,i) 
      end do
      
      mIndex2 = (/(J,J=1,mNtdc)/)
      
!      CALL mexprintf('BIG Before CovsrtXc'//CHAR(10))
!      CALL ECHO(BIG)
!     sort BIG by decreasing cond. variance for Xc 
      CALL CVSRTXC(mNt,mNd,mBIG,mIndex2,INFORM) 
!      CALL mexprintf('BIG after CovsrtXc'//CHAR(10))
!      CALL ECHO(BIG)
      
      IF (INFORM.GT.0) return ! degenerate case exit VALS=0 for all  
                                ! (should perhaps return NaN instead??)


      DO I=mNtdc,1,-1     
         J = mIndex2(I)            ! covariance matrix according to index2
         mXedni(J) = I
      END DO

      IF (mNx>1) THEN
         ALLOCATE(mBIG2(mNtdc,mNtdc))
         do i = 1,mNtdc
            mBIG2(1:i,i)    = mBIG(1:i,i) !Copy input matrix
         end do
      ENDIF
      return
      end subroutine preInit
      subroutine initIntegrand(ix,Xc,Ex,indI,Blo,Bup,INFIN,
     &     fxc,value,abserr,NDIM,inform)
      integer, intent(in) :: ix ! integrand number
      double precision, dimension(:),intent(in) :: Ex
      double precision, dimension(:,:), intent(in) :: Xc,Blo,Bup
      integer, dimension(:), intent(in) :: indI,INFIN
      double precision, intent(out) :: fxc,value,abserr
      integer, intent(out) :: NDIM, inform
!     Locals
      DOUBLE PRECISION   :: SQ0,xx,quant
      integer :: I,J
      inform = 0
      NDIM   = 0
      VALUE   = gZERO     
      fxc     = gONE
      abserr  = mSmall
         
      IF (mInitIntegrandCalled)  then
         do i = 1,mNtdc
            mBIG(1:i,i)    = mBIG2(1:i,i) !Copy input matrix
         end do
      else
         mInitIntegrandCalled = .TRUE.
      endif
      
                                ! Set the original means of the variables
      mCm(:)  = Ex(mIndex2(1:mNtdc)) !   Cm(1:Ntdc)  =Ex (index1(1:Ntdc))
      IF (mNc>0) THEN
         mXc(:) = Xc(:,ix)
                                !mXc(1:Nc)    = Xc(1:Nc,ix)
         QUANT = DBLE(mNc)*LOG(gSQTWPI1)
         I = mNtdc 
         DO J = 1, mNc        
!     Iterative conditioning on the last Nc variables  
            SQ0 = mBIG(I,I)     ! SQRT(Var(X(i)|X(i+1),X(i+2),...,X(Ntdc)))
            xx = (mXc(mIndex2(I) - mNtd) - mCm(I))/SQ0
                                !Trick to calculate
                                !fxc = fxc*SQTWPI1*EXP(-0.5*(XX**2))/SQ0   
            QUANT = QUANT - gHALF*xx*xx  - LOG(SQ0)            
				
                                ! conditional mean (expectation) 
                                ! E(X(1:i-1)|X(i),X(i+1),...,X(Ntdc)) 
            mCm(1:I-1) = mCm(1:I-1) + xx*mBIG(1:I-1,I)
            I = I-1
         ENDDO
                                ! Calculating the  
                                ! fxc probability density for i=Ntdc-J+1, 
                                ! fXc=f(X(i)|X(i+1),X(i+2)...X(Ntdc))*
                                !     f(X(i+1)|X(i+2)...X(Ntdc))*..*f(X(Ntdc))
         fxc = EXP(QUANT+mXcScale) 
         
                                ! if fxc small:  don't bother 
                                ! calculating it, goto end
         IF (fxc  < mFxcEpss) then
            abserr = gONE
            inform = 1
            return
         endif
      END IF
!     Set integration limits mA,mB and mINFI
!     NOTE: mA and mB are integration limits with mCm subtracted
      CALL setIntLimits(mXc,indI,Blo,Bup,INFIN,inform)
      if (inform>0) return
      mIndex1(:) = mIndex2(:)   
      CALL COVSRT(.FALSE., mNt,mNd,mBIG,mCm,mA,mB,mINFI, 
     &        mINDEX1,mINFIXt,mINFIXd,NDIM,mY,mCDI)

      CALL INITFUN(VALUE,abserr,INFORM)
!     IF INFORM>0 : degenerate case:
!     Integral can be calculated excactly, ie. 
!     mean of deterministic variables outside the barriers,
!     or NDIM = 1 
      return
      end subroutine initIntegrand
      subroutine cleanUp
!     Deallocate all work arrays and vectors      
      IF (ALLOCATED(mXc))     DEALLOCATE(mXc) 
      IF (ALLOCATED(mXd))     DEALLOCATE(mXd) 
      IF (ALLOCATED(mCm))     DEALLOCATE(mCm) 
      IF (ALLOCATED(mBIG2))   DEALLOCATE(mBIG2) 
      IF (ALLOCATED(mBIG))    DEALLOCATE(mBIG) 
      IF (ALLOCATED(mIndex2)) DEALLOCATE(mIndex2)
      IF (ALLOCATED(mIndex1)) DEALLOCATE(mIndex1)
      IF (ALLOCATED(mXedni))  DEALLOCATE(mXedni)
      IF (ALLOCATED(mA))      DEALLOCATE(mA)
      IF (ALLOCATED(mB))      DEALLOCATE(mB)  
      IF (ALLOCATED(mY))      DEALLOCATE(mY) 
      IF (ALLOCATED(mCDI))    DEALLOCATE(mCDI)  
      IF (ALLOCATED(mCDIXd))  DEALLOCATE(mCDIXd)  
      IF (ALLOCATED(mCmXd))   DEALLOCATE(mCmXd)  
      IF (ALLOCATED(mINFI))   DEALLOCATE(mINFI)
      end subroutine cleanUp
      function integrandBound(I0,N,Y,FINY) result (bound1)
      use FIMOD
      integer, intent(in) :: I0,N,FINY
      double precision, intent(in) :: Y
      double precision :: bound1
! locals
      integer :: I,IK,FINA, FINB
      double precision :: AI,BI,D1,E1
      double precision :: upError,loError, TMP
!     Computes the upper bound for the intgrand
      bound1 = gzero
      if (FINY<1) return
      FINA = 0
      FINB = 0
      IK = 2
      DO I = I0, N
             ! E(Y(I) | Y(1))/STD(Y(IK)|Y(1))   
         TMP = mBIG(IK-1,I)*Y 
         IF (mINFI(I) > -1) then
!     May have infinite int. Limits if Nd>0
            IF ( mINFI(I) .NE. 0 ) THEN
               IF ( FINA .EQ. 1 ) THEN
                  AI = MAX( AI, mA(I) - tmp )
               ELSE
                  AI   = mA(I) - tmp        
                  FINA = 1
               END IF
            END IF
            IF ( mINFI(I) .NE. 1 ) THEN
               IF ( FINB .EQ. 1 ) THEN
                  BI = MIN( BI, mB(I) - tmp)
               ELSE
                  BI   = mB(I) - tmp       
                  FINB = 1
               END IF
            END IF
         endif
        
         IF (I.EQ.N.OR.mBIG(IK+1,I+1)>gZERO) THEN
            CALL MVNLMS( AI, BI,2*FINA+FINB-1, D1, E1 )
            IF (D1<E1) bound1 = E1-D1
            return
         ENDIF
      ENDDO
      RETURN
      end function integrandBound
      SUBROUTINE INITFUN(VALUE,abserr,INFORM)
      USE JACOBMOD
      use SWAPMOD
      USE FIMOD
!      USE GLOBALDATA, ONLY: NIT,EPS2,EPS,xCutOff,NC1C2,ABSEPS
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(OUT) :: VALUE,abserr
      INTEGER, INTENT(out) :: INFORM
! local variables:
      INTEGER ::  N,NdleftO
      INTEGER :: I, J,  FINA, FINB, INFII
      DOUBLE PRECISION :: AI, BI, D0,E0,R12, R13, R23
      DOUBLE PRECISION :: xCut , upError,loError,maxTruncError
      LOGICAL :: useC1C2,isXd
     
!     
!     Integrand subroutine
!
!INITFUN initialize the Multivariate Normal integrand function
! COF  - conditional sorted ChOlesky Factor of the covariance matrix (IN)
! CDI  - Cholesky DIagonal elements used to calculate the mean  
! Cm   - conditional mean of Xd and Xt given Xc, E(Xd,Xt|Xc)
! xd   - variables to the jacobian variable, need no initialization size Nd
! xc   - conditional variables (IN)
! INDEX1 - if INDEX1(I)>Nt then variable no. I is one of the Xd
!          variables otherwise it is one of Xt.

      !PRINT *,'Mvnfun,ndim',Ndim
      INFORM = 0
      VALUE  = gZERO
      abserr = max(mCovEps , 6.0d0*mSmall)
      mVAL0   = gONE
     

      mNdleftN0 = mNd               ! Counter for number of Xd variables left
      
      mI0  = 0
      FINA = 0
      FINB = 0
      N = mNt + mNd - mINFIXt - mINFIXd-1
      IF (mINFIXt+mINFIXd > 0) THEN
!     CHCKLIM Check if the conditional mean Cm = E(Xt,Xd|Xc) for the
!     deterministic variables are between the barriers, i.e.,
!     A=Hlo-Cm< 0 <B=Hup-Cm
!     INFIN  INTEGER, array of integration limits flags:
!            if INFIN(I) < 0, Ith limits are (-infinity, infinity);
!            if INFIN(I) = 0, Ith limits are (-infinity, B(I)];
!            if INFIN(I) = 1, Ith limits are [A(I), infinity);
!            if INFIN(I) = 2, Ith limits are [A(I), B(I)].
		
         I = N+1
         DO J=1, mINFIXt + mINFIXd
            I = I + 1
            IF (mINFI(I)>-1) THEN
               IF ((mINFI(I).NE.0).AND.(mAbsEps  < mA(I))) GOTO 200
               IF ((mINFI(I).NE.1).AND.(mB(I) < -mAbsEps )) GOTO 200
            ENDIF
         ENDDO
           
		
         IF (mINFIXd>0) THEN       
            ! Redundant variables of Xd: replace Xd with the mean
            I = mNt + mNd !-INFIS
            J = mNdleftN0-mINFIXd
              
            DO WHILE (mNdleftN0>J)
               isXd = (mNt < mIndex1(I))
               IF (isXd) THEN 
                  mXd (mNdleftN0) =  mCm (I)
                  mNdleftN0 = mNdleftN0-1                  
               END IF
               I = I-1
            ENDDO
         ENDIF
		
         IF (N+1 < 1) THEN   
!     Degenerate case, No relevant variables left to integrate  
!     Print *,'rindd ndim1',Ndim1
            IF (mNd>0) THEN 
               VALUE = jacob (mXd,mXc) ! jacobian of xd,xc
            ELSE
               VALUE = gONE
            END IF
            GOTO 200
         ENDIF
      ENDIF
      IF (mNIT<=100) THEN
         xCut = mXcutOff
      
         J = 1
         DO I = 2, N+1
            IF (mBIG(J+1,I)>gZERO) THEN
               J = J + 1
            ELSE
               ! Add xCut std to deterministic variables to get an upper
               ! bound for integral
              mA(I) =  mA(I) - xCut * mBIG(I,J)
              mB(I) =  mB(I) + xCut * mBIG(I,J)
            ENDIF
         END DO
      ELSE
         xCut = gZERO 
      ENDIF

      NdleftO = mNdleftN0
      useC1C2 = (1<=mNc1c2)
      DO I = 1, N+1
         IF (mINFI(I) > -1) then
!     May have infinite int. Limits if Nd>0
            IF ( mINFI(I) .NE. 0 ) THEN
               IF ( FINA .EQ. 1 ) THEN
                  AI = MAX( AI, mA(I) )
               ELSE
                  AI   = mA(I)        
                  FINA = 1
               END IF
            END IF
            IF ( mINFI(I) .NE. 1 ) THEN
               IF ( FINB .EQ. 1 ) THEN
                  BI = MIN( BI, mB(I) )
               ELSE
                  BI   = mB(I)       
                  FINB = 1
               END IF
            END IF
         endif
         isXd = (mINDEX1(I)>mNt)
         IF (isXd) THEN         ! Save the mean for Xd
            mCmXd(mNdleftN0)  = mCm(I) 
            mCDIXd(mNdleftN0) = mCDI(I)  
            mNdleftN0        = mNdleftN0-1
         END IF
    
         IF (I.EQ.N+1.OR.mBIG(2,I+1)>gZERO) THEN
            IF (useC1C2.AND.I<N) THEN
               mY(:) = gZERO
               
               
               CALL MVNLMS( AI, BI,2*FINA+FINB-1, D0, E0 )
               IF (D0>=E0) GOTO 200

               CALL C1C2(I+1,N+1,1,mA,mB,mINFI,mY,mBIG,AI,BI,FINA,FINB)
               CALL MVNLMS( AI, BI,2*FINA+FINB-1, mD1, mE1 )
               IF (mD1>=mE1) GOTO 200
               maxTruncError = FI(-ABS(mXcutOff))*dble(mNc1c2)
               upError = abs(E0-mE1)
               loError = abs(D0-mD1)
               if (upError>mSmall) then
                  upError = upError*integrandBound(I+1,N+1,BI,FINB) 
               endif             
               if (loError>mSmall) then
                  loError = loError*integrandBound(I+1,N+1,AI,FINA)
               endif
               abserr  = abserr + min(upError + loError,maxTruncError)
               !CALL printvar(log10(loError+upError+msmall),'lo+up-err')
            ELSE
               CALL MVNLMS( AI, BI,2*FINA+FINB-1, mD1, mE1 )
               IF (mD1>=mE1) GOTO 200  
            ENDIF
            !CALL MVNLMS( AI, BI,2*FINA+FINB-1, mD1, mE1 )
            !IF (mD1>=mE1) GOTO 200
            IF ( NdleftO<=0) THEN
               IF (mNd>0) mVAL0 = JACOB(mXd,mXc) 
               SELECT CASE (I-N)
               CASE (1)   !IF (I.EQ.N+1) THEN
                  VALUE  = (mE1-mD1)*mVAL0
                  abserr = abserr*mVAL0
                  GO TO 200
               CASE (0)     !ELSEIF (I.EQ.N) THEN
                                !D1=1/sqrt(1-rho^2)=1/STD(X(I+1)|X(1))
                  mD1 = SQRT( gONE + mBIG(1,I+1)*mBIG(1,I+1) ) 
                  mINFI(2) = mINFI(I+1)
                  mA(1) = AI
                  mB(1) = BI
                  mINFI(1) = 2*FINA+FINB-1
                  IF ( mINFI(2) .NE. 0 ) mA(2) = mA(I+1)/mD1
                  IF ( mINFI(2) .NE. 1 ) mB(2) = mB(I+1)/mD1
                  VALUE = BVNMVN( mA, mB,mINFI,mBIG(1,I+1)/mD1 )*mVAL0
                  abserr = (abserr+1.0d-14)*mVAL0
                  GO TO 200 
               CASE ( -1 )  !ELSEIF (I.EQ.N-1) THEN
                  IF (.FALSE.) THEN
! TODO :this needs further checking! (it should work though) 
                  !1/D1= sqrt(1-r12^2) = STD(X(I+1)|X(1))
                  !1/E1=  STD(X(I+2)|X(1)X(I+1))
                  !D1  = BIG(I+1,1)
                  !E1  = BIG(I+2,2)
                     
                  mD1 = gONE/SQRT( gONE + mBIG(1,I+1)*mBIG(1,I+1) )
                  R12 = mBIG( 1, I+1 ) * mD1
                  if (mBIG(3,I+2)>gZERO) then
                     mE1 = gONE/SQRT( gONE + mBIG(1,I+2)*mBIG(1,I+2) +
     &                    mBIG(2,I+2)*mBIG(2,I+2) )                
                     R13 = mBIG( 1, I+2 ) * mE1
                     R23 = mBIG( 2, I+2 ) * (mE1 * mD1) + R12 * R13
                  else
                     mE1  = mCDI(I+2)
                     R13 = mBIG( 1, I+2 ) * mE1
                     R23 = mE1*mD1 + R12 * R13
                     IF ((mE1  <  gZERO).AND. mINFI(I+2)>-1) THEN
                        CALL SWAP(mA(I+2),mB(I+2))
                        IF (mINFI(I+2).NE. 2) mINFI(I+2) = 1-mINFI(I+2)
                     END IF
                     !R23 = BIG( 2, I+2 ) * (E1 * D1) + R12 * R13
                  endif
                  mINFI(2) = mINFI(I+1)
                  mINFI(3) = mINFI(I+2)
                  mA(1) = AI
                  mB(1) = BI
                  mINFI(1) = 2*FINA+FINB-1
                  IF ( mINFI(2) .NE. 0 ) mA(2) = mA(I+1) * mD1
                  IF ( mINFI(2) .NE. 1 ) mB(2) = mB(I+1) * mD1
                  IF ( mINFI(3) .NE. 0 ) mA(3) = mA(I+2) * mE1
                  IF ( mINFI(3) .NE. 1 ) mB(3) = mB(I+2) * mE1
                  if(.false.) then
                     CALL PRINTVECD((/R12, R13, R23 /),'R12 = ')
                     CALL PRINTVECD((/mD1, mE1 /),'D1 = ')
                     CALL PRINTVECD(mBIG(1,1:3),'BIG(1,1:3) = ')
                     CALL PRINTVECD(mBIG(2,2:3),'BIG(2,2:3) = ')
                     CALL PRINTVECD(mBIG(1:3,1),'BIG(1:3,1) = ')
                     CALL PRINTVECD(mBIG(2:3,2),'BIG(2:3,2) = ')
                     CALL PRINTVECD(mA(1:I+2),'A = ')
                     CALL PRINTVECD(mB(1:I+2),'B = ')
                     CALL PRINTVECI(mINFI(1:I+2),'INFI = ')
                     CALL PRINTVECI(mINDEX1(1:I+2),'index1 = ')
                  endif
                  VALUE = TVNMVN( mA, mB,mINFI,
     &                 (/R12, R13, R23 /),1.0d-13) * mVAL0
                  ABSERR = (ABSERR + 1.0d-13)*mVAL0
                  GOTO 200
                  ENDIF
               END SELECT !ENDIF
            ENDIF
            ABSERR = mVAL0*ABSERR
            mVAL0 = mVAL0 * (mE1-mD1)
            mI0   = I
            RETURN
         ENDIF
      ENDDO
      RETURN
 200  INFORM = 1
      RETURN
      END SUBROUTINE INITFUN     
!     
!     Integrand subroutine
!
      FUNCTION MVNFUN( Ndim, W ) RESULT (VAL)
      USE JACOBMOD
      USE FIMOD
      IMPLICIT NONE
      INTEGER, INTENT (IN) :: Ndim
      DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: W
      DOUBLE PRECISION :: VAL
! local variables:
      INTEGER ::  N,I, J, FINA, FINB
      INTEGER ::  NdleftN, NdleftO ,IK
      DOUBLE PRECISION :: TMP, AI, BI, DI, EI
      LOGICAL :: useC1C2, isXd      
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
! mBIG - conditional sorted ChOlesky Factor of the covariance matrix (IN)
! mCDI - Cholesky DIagonal elements used to calculate the mean  
! mCm  - conditional mean of Xd and Xt given Xc, E(Xd,Xt|Xc)
! mXd  - variables to the jacobian variable, need no initialization size Nd
! mXc  - conditional variables (IN)
! mINDEX1 - if mINDEX1(I)>Nt then variable No. I is one of the Xd
!           variables otherwise it is one of Xt

      !PRINT *,'Mvnfun,ndim',Ndim
      
!      xCut = gZERO ! xCutOff

      N    = mNt+mNd-mINFIXt-mINFIXd-1
      IK   = 1                    ! Counter for Ndim 
      FINA = 0
      FINB = 0

      NdleftN = mNdleftN0          ! Counter for number of Xd variables left
      VAL     = mVAL0
      NdleftO = mNd - mINFIXd
      mY(IK) = FIINV( mD1 + W(IK)*( mE1 - mD1 ) )
      useC1C2 = (IK+1.LE.mNc1c2)
      IF (useC1C2) THEN
         ! Calculate the conditional mean
         ! E(Y(I) | Y(1),...Y(I0))/STD(Y(I)|Y(1),,,,Y(I0))
         mY(mI0+1:N+1) = mBIG(IK, mI0+1:N+1)*mY(IK) 
      ENDIF
      IF (NdleftO.GT.NdleftN ) THEN
         mXd(NdleftN+1:NdleftO) = mCmXd(NdleftN+1:NdleftO)+
     &        mY(IK) * mCDIXd(NdleftN+1:NdleftO)
      ENDIF
      NdleftO = NdleftN
      IK = 2                    !=IK+1
            
      
      DO I = mI0+1, N+1
         IF (useC1C2) THEN
             TMP = mY(I)
          ELSE
            TMP = 0.d0
            DO J = 1, IK-1 
               ! E(Y(I) | Y(1),...Y(IK-1))/STD(Y(IK)|Y(1),,,,Y(IK-1))
               TMP = TMP + mBIG(J,I)*mY(J)  
            END DO
         ENDIF
         IF (mINFI(I) < 0) GO TO 100
            ! May have infinite int. Limits if Nd>0
         IF ( mINFI(I) .NE. 0 ) THEN
            IF ( FINA .EQ. 1 ) THEN
               AI = MAX( AI, mA(I) - TMP) 
            ELSE
               AI = mA(I) - TMP 
               FINA = 1
            END IF
            IF (FINB.EQ.1.AND.BI<=AI) GOTO 200
         END IF
         IF ( mINFI(I) .NE. 1 ) THEN
            IF ( FINB .EQ. 1 ) THEN
               BI = MIN( BI, mB(I) - TMP) 
            ELSE
               BI = mB(I) - TMP
               FINB = 1
            END IF
            IF (FINA.EQ.1.AND.BI<=AI) GOTO 200
         END IF
 100     isXd = (mNt<mINDEX1(I))
         IF (isXd) THEN 
!     Save the mean of xd and Covariance diagonal element 
            ! Conditional mean E(Xi|X1,..X)
            mCmXd(NdleftN)  = mCm(I) + TMP * mCDI(I) 
             ! Covariance diagonal 
            mCDIXd(NdleftN) = mCDI(I)              
            NdleftN        = NdleftN - 1
         END IF 
         IF (I == N+1 .OR. mBIG(IK+1,I+1) > gZERO ) THEN 
            IF (useC1C2) THEN
!     Note: for J =I+1:N+1:  Y(J) = conditional expectation, E(Yj|Y1,...Yk) 
               CALL C1C2(I+1,N+1,IK,mA,mB,mINFI,mY,mBIG,AI,BI,FINA,FINB)
            ENDIF
            CALL MVNLMS( AI, BI, 2*FINA+FINB-1, DI, EI )            
            IF ( DI >= EI ) GO TO 200
            VAL = VAL * ( EI - DI )
            
            IF ( I <= N .OR. (NdleftN < NdleftO)) THEN
               mY(IK) = FIINV( DI + W(IK)*( EI - DI ) )
               IF (NdleftN < NdleftO ) THEN
                  mXd(NdleftN+1:NdleftO) = mCmXd(NdleftN+1:NdleftO)+
     &                 mY(IK) * mCDIXd(NdleftN+1:NdleftO)
                  NdleftO = NdleftN
               ENDIF
               useC1C2 = (IK+1<=mNc1c2)
               IF (useC1C2) THEN
                  
                  ! E(Y(J) | Y(1),...Y(I))/STD(Y(J)|Y(1),,,,Y(I))
                  mY(I+1:N+1) = mY(I+1:N+1) + mBIG(IK, I+1:N+1)*mY(IK)
               ENDIF
            ENDIF            
            IK   = IK + 1
            FINA = 0
            FINB = 0
         END IF
      END DO
      IF (mNd>0) VAL = VAL * jacob(mXd,mXc)      
      RETURN
 200  VAL = gZERO
      RETURN
      END FUNCTION MVNFUN


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!******************* RINDD - the main program *********************!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE RINDD (VALS,ERR,Big,Ex,Xc,Nt,
     &	              indI,Blo,Bup,INFIN)  
      USE RCRUDEMOD
      USE KRBVRCMOD
      USE ADAPTMOD
      USE KROBOVMOD
	USE DKBVRCMOD
      USE SSOBOLMOD
      IMPLICIT NONE  
      DOUBLE PRECISION, DIMENSION(:  ), INTENT(out):: VALS, ERR 
      DOUBLE PRECISION, DIMENSION(:,:), INTENT(in) :: BIG
      DOUBLE PRECISION, DIMENSION(:,:), INTENT(in) :: Xc 
      DOUBLE PRECISION, DIMENSION(:  ), INTENT(in) :: Ex            
      DOUBLE PRECISION, DIMENSION(:,:), INTENT(in) :: Blo, Bup  
      INTEGER,          DIMENSION(:  ), INTENT(in) :: indI,INFIN
      INTEGER,                          INTENT(in) :: Nt 
!      DOUBLE PRECISION,                 INTENT(in) :: XcScale
! local variables
      INTEGER :: ix, INFORM, NDIM, Ndleft ,MAXPTS,MINPTS
      DOUBLE PRECISION :: VALUE,fxc,absERR,absERR2
      double precision :: LABSEPS,LRELEPS
      

      VALS(:) = gZERO    
      ERR(:)  = gONE    
     
      call preInit(BIG,Xc,Nt,inform)      
      IF (INFORM.GT.0) GOTO 110 ! degenerate case exit VALS=0 for all  
                                ! (should perhaps return NaN instead??)

!     Now the loop over all different values of
!     variables Xc (the one one is conditioning on)  
!     is started. The density f_{Xc}(xc(:,ix))
!     will be computed and denoted by  fxc.
      DO  ix = 1, mNx 
         call initIntegrand(ix,Xc,Ex,indI,Blo,Bup,infin,
     &        fxc,value,abserr,NDIM,inform)

        
         IF (INFORM.GT.0) GO TO 100 
         
         MAXPTS  = mMAXPTS
         MINPTS  = mMINPTS
         LABSEPS = max(mABSEPS-abserr,0.2D0*mABSEPS)       !*fxc
         LRELEPS = mRELEPS
         absErr2 = mSmall
         
         SELECT CASE (mMethod)
         CASE (:1)
            IF (NDIM < 9) THEN
               CALL SADAPT(NDIM,MAXPTS,MVNFUN,LABSEPS,
     &              LRELEPS,ABSERR2,VALUE,INFORM)
               VALUE = MAX(VALUE,gZERO)
            ELSE
               CALL KRBVRC(NDIM, MINPTS, MAXPTS, MVNFUN,LABSEPS,LRELEPS,
     &              ABSERR2, VALUE, INFORM )
            ENDIF
         CASE (2)               
!        Call the subregion adaptive integration subroutine
            IF ( NDIM .GT. 19.) THEN
!     print *, 'Ndim too large for SADMVN => Calling KRBVRC'
               CALL KRBVRC( NDIM, MINPTS, MAXPTS, MVNFUN, LABSEPS,
     &              LRELEPS, ABSERR2, VALUE, INFORM )
            ELSE
               CALL SADAPT(NDIM,MAXPTS,MVNFUN,LABSEPS,
     &              LRELEPS,ABSERR2,VALUE,INFORM)
               VALUE = MAX(VALUE,gZERO)
            ENDIF
         CASE (3)               !       Call the Lattice rule integration procedure
            CALL KRBVRC( NDIM, MINPTS, MAXPTS, MVNFUN, LABSEPS,
     &           LRELEPS, ABSERR2, VALUE, INFORM )
         CASE (4)               !       Call the Lattice rule
                                !       integration procedure 
            CALL KROBOV( NDIM, MINPTS, MAXPTS, MVNFUN, LABSEPS,
     &           LRELEPS,ABSERR2, VALUE, INFORM )
         CASE (5)    ! Call Crude Monte Carlo integration procedure
            CALL RANMC( NDIM, MAXPTS, MVNFUN, LABSEPS, 
     &           LRELEPS, ABSERR2, VALUE, INFORM )           
         CASE (6)       !       Call the scrambled Sobol sequence rule integration procedure
           CALL SOBNIED( NDIM, MINPTS, MAXPTS, MVNFUN, LABSEPS, LRELEPS,
     &           ABSERR2, VALUE, INFORM )
	   CASE (7:)
	     CALL DKBVRC( NDIM, MINPTS, MAXPTS, MVNFUN, LABSEPS, LRELEPS,
     &           ABSERR2, VALUE, INFORM )
         END SELECT   

!     IF (INFORM.gt.0) print *,'RINDD, INFORM,error =',inform,error
 100     VALS(ix) = VALUE*fxc 
         IF (SIZE(ERR, DIM = 1).EQ.mNx) ERR(ix)  = (abserr+abserr2)*fxc       
      ENDDO                     !ix

 110  CONTINUE
      call cleanUp
      RETURN                                                            
      END SUBROUTINE RINDD
      
      SUBROUTINE setIntLimits(xc,indI,Blo,Bup,INFIN,inform)     
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(:  ), INTENT(in)  :: xc
      INTEGER,          DIMENSION(:  ), INTENT(in)  :: indI,INFIN    
      DOUBLE PRECISION, DIMENSION(:,:), INTENT(in)  :: Blo,Bup   
      integer, intent(out) :: inform
!Local variables
      INTEGER   :: I, J, K, L,Mb1,Nb,NI,Nc 
      DOUBLE PRECISION :: xCut, SQ0
!this procedure set mA,mB and mInfi according to Blo/Bup and INFIN
!
!     INFIN  INTEGER, array of integration limits flags:
!            if INFIN(I) < 0, Ith limits are (-infinity, infinity);
!            if INFIN(I) = 0, Ith limits are (-infinity, mB(I)];
!            if INFIN(I) = 1, Ith limits are [mA(I), infinity);
!            if INFIN(I) = 2, Ith limits are [mA(I), mB(I)].
! Note on member variables:
! mXedni     = indices to the variables new place after cvsrtXc.  Size Ntdc
! mCm        = E(Xt,Xd|Xc), i.e., conditional mean given Xc
! mBIG(:,1:Ntd) = Cov(Xt,Xd|Xc)

      xCut = ABS(mInfinity)
      Mb1 = size(Blo,DIM=1)-1
      Nb = size(Blo,DIM=2)
      NI = size(indI,DIM=1)
      Nc = size(xc,DIM=1)
      if (Mb1>Nc .or. Nb.NE.NI-1) then
!     size of variables inconsistent         
         inform = 4
         return
      endif
      
!      IF (Mb.GT.Nc+1) print *,'barrier: Mb,Nc =',Mb,Nc
!      IF (Nb.NE.NI-1) print *,'barrier: Nb,NI =',Nb,NI
      DO J = 2, NI 
         DO I = indI (J - 1) + 1 , indI (J)
            L        = mXedni(I)  
            mINFI(L) = INFIN(J-1)
            SQ0      = SQRT(mBIG(L,L))
            mA(L)    = -xCut*SQ0
            mB(L)    =  xCut*SQ0
            IF (mINFI(L).GE.0) THEN
               IF  (mINFI(L).NE.0) THEN
                  mA(L) = Blo (1, J - 1)-mCm(L)  
                  DO K = 1, Mb1
                     mA(L) = mA(L)+Blo(K+1,J-1)*xc(K)              
                  ENDDO         ! K
                  ! This can only be done if 
                  if (mA(L)< -xCut*SQ0) mINFI(L) = mINFI(L)-2
               ENDIF
               IF  (mINFI(L).NE.1) THEN
                  mB(L) = Bup (1, J - 1)-mCm(L)  
                  DO K = 1, Mb1
                     mB(L) = mB(L)+Bup(K+1,J-1)*xc(K)
                  ENDDO  
                  if (xCut*SQ0<mB(L)) mINFI(L) = mINFI(L)-1
               ENDIF            !
            ENDIF            
         ENDDO                  ! I
      ENDDO                     ! J
!     print * ,'barrier hup:',size(Hup),Hup(xedni(1:indI(NI)))
!     print * ,'barrier hlo:',size(Hlo),Hlo(xedni(1:indI(NI)))
      RETURN
      END SUBROUTINE setIntLimits                                              
      

      FUNCTION GETTMEAN(A,B,INFJ,PRB) RESULT (MEAN1)
      USE GLOBALCONST
      IMPLICIT NONE
!     GETTMEAN Returns the expected mean, E(I(a<x<b)*X)
      DOUBLE PRECISION, INTENT(IN) :: A,B,PRB
      INTEGER , INTENT(IN) :: INFJ
      DOUBLE PRECISION :: MEAN1
      DOUBLE PRECISION :: YL,YU
!     DOUBLE PRECISION, PARAMETER:: ZERO = 0.0D0, HALF = 0.5D0
      
      IF ( PRB .GT. gZERO) THEN
         YL = gZERO
         YU = gZERO
         IF (INFJ.GE.0) THEN
            IF (INFJ .NE. 0) YL =-EXP(-gHALF*(A*A))*gSQTWPI1
            IF (INFJ .NE. 1) YU =-EXP(-gHALF*(B*B))*gSQTWPI1
         ENDIF
         MEAN1 = ( YU - YL )/PRB
      ELSE
         SELECT CASE (INFJ)
         CASE (:-1) 
            MEAN1 = gZERO
         CASE (0) 
            MEAN1 = B
         CASE (1)
            MEAN1 = A
         CASE (2:)
            MEAN1 = ( A + B ) * gHALF
         END SELECT
      END IF
      RETURN
      END FUNCTION
      SUBROUTINE ADJLIMITS(A,B, infi)
!      USE GLOBALDATA, ONLY : xCutOff
      IMPLICIT NONE
!     Adjust INFI when integration limits A and/or B is too far out in the tail
      DOUBLE PRECISION, INTENT(IN)     :: A,B
      INTEGER,          INTENT(IN OUT) :: infi
!      DOUBLE PRECISION, PARAMETER :: xCutOff = 8.D0
      IF (infi>-1) THEN
         IF (infi.NE.0)THEN
            IF (A  <  -mXcutOff) THEN
               infi = infi-2
!               CALL mexprintf('ADJ A')
            ENDIF
         ENDIF
         IF (infi.NE.1) THEN
            IF (mXCutOff  <  B) THEN
               infi = infi-1
!               CALL mexprintf('ADJ B')
            ENDIF
         END IF
      END IF
      RETURN
      END SUBROUTINE ADJLIMITS
      SUBROUTINE C1C2(I0,I1,IK,A,B,INFIN, Cm, BIG, AJ, BJ, FINA,FINB)  
! The regression equation for the conditional distr. of Y given X=x
! is equal  to the conditional expectation of Y given X=x, i.e.,
! 
!       E(Y|X=x) = E(Y) + Cov(Y,X)/Var(X)[x-E(X)]
!
!  Let 
!     x1 = (x-E(X))/SQRT(Var(X)) be zero mean, 
!     C1< x1 <C2, 
!     B1(I) = COV(Y(I),X)/SQRT(Var(X)). 
!  Then the process  Y(I) with mean Cm(I) can be written as 
!
!       y(I) = Cm(I) + x1*B1(I) + Delta(I) for  I=I0,...,I1.
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
!      USE GLOBALDATA, ONLY : EPS2,EPS ,xCutOff
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: Cm, A,B !, B1, SQ
      DOUBLE PRECISION, DIMENSION(:,:), INTENT(in) :: BIG 
      INTEGER,          DIMENSION(:), INTENT(in) :: INFIN 
      DOUBLE PRECISION,            INTENT(inout) :: AJ,BJ
      INTEGER,                     INTENT(inout) :: FINA, FINB
      INTEGER, INTENT(IN) :: I0, I1, IK
     
!     Local variables 
      DOUBLE PRECISION   :: xCut
      DOUBLE PRECISION, PARAMETER :: TOL = 1.0D-16     
      DOUBLE PRECISION :: AI,BI,CSQ,BdSQ0, LTOL 
      INTEGER :: INFI, I

      xCut = MIN(ABS(mXcutOff),mInfinity)
      LTOL = mSmall ! EPS2
!      AJ = MAX(AJ,-xCut)
!      BJ = MIN(BJ,xCut)
!      IF (AJ.GE.BJ) GO TO 112
!      CALL PRINTVAR(AJ,TXT='BC1C2: AJ')
!      CALL PRINTVAR(BJ,TXT='BC1C2: BJ')
     
      IF (I1 < I0)  RETURN       !Not able to change integration limits
      DO I = I0,I1             
!     C = xCutOff
         INFI = INFIN(I)
         IF (INFI>-1) THEN
            !BdSQ0 = B1(I)  
            !CSQ   = xCut * SQ(I)
            BdSQ0  = BIG(IK,I)
            CSQ    = xCut * BIG(I,IK)
            IF (BdSQ0 > LTOL) THEN
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
             ELSEIF (BdSQ0  <  -LTOL) THEN
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
!      IF (FINA>0 .AND. FINB>0) THEN
!         IF (AJ<BJ) THEN
!            IF (AJ   <= -xCut) FINA = 0
!            IF (xCut <= BJ   ) FINB = 0
!         ENDIF
!      ENDIF
!      CALL PRINTVAR(AJ,TXT='AC1C2: AJ')
!      CALL PRINTVAR(BJ,TXT='AC1C2: BJ')
      RETURN
      END SUBROUTINE C1C2
      SUBROUTINE CVSRTXC (Nt,Nd,R,index1,INFORM)
!      USE GLOBALDATA, ONLY :  XCEPS2
!      USE GLOBALCONST
      IMPLICIT NONE
      INTEGER, INTENT(in) :: Nt,Nd 
      DOUBLE PRECISION, DIMENSION(:,:), INTENT(inout) :: R
      INTEGER,          DIMENSION(:  ), INTENT(inout) :: index1 
      INTEGER, INTENT(out) :: INFORM 
! local variables
      DOUBLE PRECISION, DIMENSION(:), allocatable   :: SQ
      INTEGER,          DIMENSION(1)                :: m
      INTEGER :: M1,K,I,J,Ntdc,Ntd,Nc, Nullity,LO
      DOUBLE PRECISION :: LTOL, maxSQ
!     if any Var(Xc(j)|Xc(1),...,Xc(j-1)) <= XCEPS2 then return NAN 
      double precision :: XCEPS2 
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
      
      IF (Nc < 1) RETURN
      
      
      
      ALLOCATE(SQ(1:Ntdc))
      maxSQ = gZERO
      DO I = 1, Ntdc
         SQ(I)  = R(I,I)
         if (SQ(I)>maxSQ) maxSQ = SQ(I)
      ENDDO

      XCEPS2 = Ntdc*mSmall*maxSQ
      mXcEps2 = XCEPS2
      LTOL   = mSmall
	
      LO = 1
      K = Ntdc
      DO I = 1, Nc             ! Condsort Xc
         m  = K+1-MAXLOC(SQ(K:Ntd+1:-1)) 
         M1 = m(1)
         IF (SQ(m1)<=XCEPS2) THEN
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
            IF (SQ(J)<=LTOL.AND.J<=Ntd) THEN
               IF (LO < J) THEN
                  CALL RCSWAP(LO, J, Ntdc,Ntd, R,INDEX1,SQ)
               ENDIF
               R(LO,LO:K-1) = gZERO
               IF (SQ(LO) < -10.0D0*SQRT(LTOL)) THEN
                  ! inform = 2
                  !R(LO,K) = gZERO
                 ! CALL mexprintf('Negative definit BIG!'//CHAR(10))
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

      SUBROUTINE  RCSCALE(chkLim,K,K0,N1,N,K1,INFIS,CDI,Cm,
     &     R,A,B,INFI,INDEX1,Y)
      USE GLOBALCONST
      USE SWAPMOD
      IMPLICIT NONE
!RCSCALE:  Scale  covariance matrix and limits
!
!   CALL  RCSCALE( k, k0, N1, N,K1, CDI,Cm,R,A, B, INFIN,index1,Y)
!
!  chkLim  = TRUE  if check if variable K is redundant
!            FALSE 
!    K     = index to variable which is deterministic,i.e.,
!            STD(Xk|X1,...Xr) = 0
!    N1    = Number of significant variables of [Xt,Xd] 
!    N     = length(Xt)+length(Xd)
!    K1    = index to current variable we are conditioning on.
!   CDI    = Cholesky diagonal elements which contains either
!               CDI(J) = STD(Xj | X1,...,Xj-1,Xc) if Xj is stochastic given
!                        X1,...Xj, Xc
!              or
!               CDI(J) = COV(Xj,Xk | X1,..,Xk-1,Xc  )/STD(Xk | X1,..,Xk-1,Xc) 
!               if Xj is determinstically determined given X1,..,Xk,Xc
!               for some k<j.
!   Cm     = conditional mean
!    R     = Matrix containing cholesky factor for
!            X = [Xt,Xd,Xc] in upper triangular part. Lower triangular
!            part contains conditional stdevs. size Ntdc X Ntdc
! INDEX1   = permutation index vector. (index to variables original place).
!    A,B   = lower and upper integration limit, respectively.
!    INFIN = if INFIN(I) < 0, Ith limits are (-infinity, infinity);
!            if INFIN(I) = 0, Ith limits are (-infinity, B(I)];
!            if INFIN(I) = 1, Ith limits are [A(I), infinity);
!            if INFIN(I) = 2, Ith limits are [A(I), B(I)].
!    Y     = work vector
!
! NOTE: RCSWAP works only on the upper triangular part of C
! + check if variable k is redundant
!  If the conditional covariance matrix diagonal entry is zero, 
!  permute limits and/or rows, if necessary.
      LOGICAL, INTENT(IN) :: chkLim
      INTEGER, INTENT(IN) :: K, K0,N
      INTEGER, INTENT(INOUT) :: N1,K1,INFIS
      DOUBLE PRECISION, DIMENSION(:),  INTENT(INOUT) :: CDI,A,B,Cm
      DOUBLE PRECISION, DIMENSION(:,:),INTENT(INOUT) :: R
      INTEGER,          DIMENSION(:),  INTENT(INOUT) :: INFI,INDEX1
      DOUBLE PRECISION, DIMENSION(:),OPTIONAL,INTENT(INOUT) :: Y
!Local variables
      DOUBLE PRECISION, PARAMETER :: LTOL = 1.0D-16
      double precision :: xCut
      DOUBLE PRECISION :: D,AK,BK,CVDIAG
      INTEGER :: KK,K00, KKold, I, J, Ntdc, INFK
      LOGICAL :: isXt
      K00 = K0
      DO WHILE( (0 < K00).AND. (ABS(R(K00,K)).LE.LTOL) )
         R(K00,K) = gZERO
         K00      = K00 - 1
      ENDDO
      IF (K00.GT.0) THEN
      !  CDI(K) = COV(Xk Xj| X1,..,Xj-1,Xc  )/STD(Xj | X1,..,Xj-1,Xc) 
         CDI(K) = R(K00,K)       
         A(K)   = A(K)/CDI(K)
         B(K)   = B(K)/CDI(K)                  
      
         IF ((CDI(K)  <  gZERO).AND. INFI(K).GE. 0) THEN
            CALL SWAP(A(K),B(K))
            IF (INFI(K).NE. 2) INFI(K) = 1-INFI(K)
         END IF
      
         
                                !Scale conditional covariances
         R(1:K00,K) = R(1:K00,K)/CDI(K) 
                                !Scale conditional standard dev.s used in regression eq.
         R(K,1:K00) = R(K,1:K00)/ABS(CDI(K)) 
            
         
         R(K00+1:K,K)   = gZERO
         !R(K,K00+1:K-1) = gZERO ! original
         R(K,K00:K-1) = gZERO    ! check this

         !
         if (chkLim.AND.K00>1) then
            ! Check if variable is redundant
            ! TODO:  this chklim-block does not work correctly yet
            xCut = mInfinity
            I = 1
            Ak = R(I,K)*xCut
            Bk = - (R(I,K))*xCut
            if (INFI(I)>=0) then
               if (INFI(I).ne.0) then
                  Ak = -(R(I,K))*MAX(A(I),-xCut)
               endif
               if (INFI(I).ne.1) then
                  Bk = - (R(I,K))*MIN(B(I),xCut)
               endif
            endif

            if (R(I,K)<gZERO) THEN
               CALL SWAP(Ak,Bk)
            endif
            !call printvar(infi(k),'infi(k)')
            !call printvar(A(k),'AK')
            !call printvar(B(k),'BK')
            !call printvar(Ak,'AK')
            !call printvar(Bk,'BK')
            INFK = INFI(K)
            Ak   = A(K)+Ak
            Bk   = B(K)+Bk
            D = gZERO
            DO I = 2, K00-1
               D = D + ABS(R(I,K))
            END DO
            CVDIAG = abs(R(k,k00))
            !call printvar(cvdiag,'cvdiag')
            Ak = (Ak + (D+cvdiag)*xCut)
            Bk = (Bk - (D+cvdiag)*xCut)
            !call printvar(Ak,'AK')
            !call printvar(Bk,'BK')
            ! If Ak<-xCut and xCut<Bk then variable Xk is redundant 
            CALL ADJLIMITS(Ak,Bk,INFK)
            
! Should change this to check against A(k00) and B(k00)
            IF (INFK < 0) THEN
               !variable is redundnant
               !                     CALL mexPrintf('AdjLim'//CHAR(10))
               IF ( K < N1 ) THEN
                  CALL RCSWAP( K, N1, N1,N, R,INDEX1,Cm, A, B, INFI)
                  
                  ! move conditional standarddeviations
                  R(K,1:K0) = R(N1,1:K0)
                  CDI(K)    = CDI(N1)
                  
                  IF (PRESENT(Y)) THEN
                     Y(K) = Y(N1)
                  ENDIF
               ENDIF
               CDI(N1)    = gZERO
               R(1:N1,N1) = gZERO
               R(N1,1:N1) = gZERO

               INFIS = INFIS+1
               N1    = N1-1
              ! CALL printvar(index1(N1),'index1(n1)')
              ! CALL mexPrintf('RCSCALE: updated N1')
              ! CALL printvar(INFIS,'INFIS ')
               return 
            END IF
         endif
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
               IF (N < Ntdc) THEN
                  ! SWAP Xc entries, i.e, Cov(Xt,Xc) and Cov(Xd,Xc)
                  DO J = N+1, Ntdc
                     CALL SWAP( R(KK,J), R(KKold,J) )
                  END DO
               ENDIF
               KKold = KK
               KK    = KK - 1
            ENDDO
         END DO
         IF (KK < K1) THEN
            K1 = K1 + 1
!            CALL mexPrintf('RCSCALE: updated K1'//CHAR(10))
         END IF
!         CALL PRINTVAR(K,TXT='K')
!         CALL PRINTVAR(KK,TXT='KK')
!         CALL PRINTVAR(K1,TXT='K1')
!         CALL PRINTVAR(K00,TXT='K00')
!         CALL PRINTVAR(K0,TXT='K0')
!         CALL PRINTCOF(N,A,B,INFI,R,INDEX1)
      ELSE 
!     Remove variable if it is conditional independent of all other variables
!         CALL mexPrintf('RCSCALE ERROR*********************'//char(10))
!         call PRINTCOF(N,A,B,INFI,R,INDEX1)
         CALL mexPrintf('RCSCALE ERROR*********************'//char(10))
      ENDIF
!      if (chkLim) then
!         call PRINTCOF(N,A,B,INFI,R,INDEX1)
!      endif
      END SUBROUTINE RCSCALE
      
      SUBROUTINE COVSRT(BCVSRT, Nt,Nd,R,Cm,A,B,INFI,INDEX1, 
     &     INFIS,INFISD, NDIM, Y, CDI )
      USE FIMOD
      USE SWAPMOD
      USE GLOBALCONST
!      USE GLOBALDATA, ONLY : EPS2,NIT,xCutOff
      IMPLICIT NONE
!COVSRT  sort integration limits and determine Cholesky factor.
!
!     Nt, Nd = size info about Xt and Xd variables.
!     R      = Covariance/Cholesky factored matrix for [Xt,Xd,Xc] (in)
!              On input: 
!              Note: Only upper triangular part is needed/used.
!               1a) the first upper triangular the Nt + Nd times Nt + Nd
!                   block contains COV([Xt,Xd]|Xc)
!                 (conditional covariance matrix for Xt and Xd given Xc)
!               2a) The upper triangular part of the Nt+Nd+Nc times Nc
!                   last block contains the cholesky matrix for Xc,
!                   i.e., 
! 
!              On output: 
!               1b) part 2a) mentioned above is unchanged, only necessary
!                 permutations according to INDEX1 is done.
!               2b) part 1a) mentioned above is changed to a special
!                  form of cholesky matrix: (N = Nt+Nd-INFIS-INFISD)
!   C = COVARIANCE
!   R(1,1) = 1
!   R(1,2:N) = [C(X1,X2)/STD(X1)/STD(X2|X1),..,C(X1,XN)/STD(X1)/STD(XN|XN-1,..,X1)]
!   R(2,2) = 1
!   R(2,3:N) =[C(X2,X3|X1)/STD(X2|X1)/STD(X3|X2,X1),..,C(X2,XN|X1)/STD(X2|X1)/STD(XN|XN-1,..,X1)]
!                  ....
!                  etc.
!               3b) The lower triangular part of R contains the
!               normalized conditional standard deviations (which is
!               used in the reqression approximation C1C2), i.e., 
!               R(2:N,1)  = [STD(X2|X1) STD(X3|X1),....,STD(XN|X1) ]/STD(X1) 
!               R(3:N,2)  = [STD(X3|X1,X2),....,STD(XN|X1,X2) ]/STD(X2|X1)  
!               .....
!               etc.
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
      double precision :: covErr
!     Local variables
      INTEGER :: N,N1,I, J, K, L, M, JMIN,Ndleft
      INTEGER ::  K1, K0, Nullity,INFJ,FINA,FINB
      DOUBLE PRECISION :: SUMSQ, AJ, BJ, TMP,D, E,EPSL
      DOUBLE PRECISION :: AA, Ca, Pa, APJ, PRBJ,RMAX
      DOUBLE PRECISION :: CVDIAG, AMIN, BMIN, PRBMIN, YL, YU
      DOUBLE PRECISION :: SQTWPI,LTOL,TOL,ONE,ZERO,HALF,xCut 
      LOGICAL          :: isOK = .TRUE.
      LOGICAL          :: isXd = .FALSE.
      LOGICAL          :: isXt = .FALSE.
      LOGICAL          :: chkLim
      PARAMETER ( SQTWPI = 2.506628274631001D0, TOL = 1D-16 )
      PARAMETER (ONE = 1.D0, ZERO = 0.D0, HALF = 0.5D0) 
      
      xCut = mInfinity 
!     xCut = MIN(ABS(xCutOff),8.0D0)
      INFIS  = 0
      INFISD = 0
      Ndim   = 0
      Ndleft = Nd
      N      = Nt + Nd

      LTOL = mSmall
      TMP  = ZERO
      DO I = 1, N
         IF (R(I,I).GT.TMP) TMP = R(I,I)
      ENDDO
      EPSL = tmp*MAX(mCovEps,N*mSmall) !tmp
      !IF (N < 10) EPSL = MIN(1D-10,EPSL)
      !LTOL = EPSL
!      EPSL = MAX(EPS2,LTOL)
      IF (TMP.GT.EPSL) THEN
         DO I = 1, N
            IF ((INFI(I)  <  0).OR.(R(I,I)<=LTOL)) THEN
               IF (INDEX1(I)<=Nt)  THEN 
                  INFIS = INFIS+1
               ELSEIF (R(I,I)<=LTOL) THEN
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

      covErr = 20.d0*LTOL

      N1  = N-INFIS-INFISD
      CDI(N1+1:N) = gZERO
      !PRINT *,'COVSRT'
      !CALL PRINTCOF(N,A,B,INFI,R,INDEX1)

!     Move any redundant variables of Xd to innermost positions. 
      DO I = N, N-INFISD+1, -1
         isXt = (INDEX1(I)<=Nt)
         IF ( (R(I,I) > LTOL) .OR. (isXt)) THEN 
            DO J = 1,I-1
               isXd = (INDEX1(J)>Nt)
               IF ( (R(J,J) <= LTOL) .AND.isXd) THEN
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
         isXd = (INDEX1(I)>Nt)
         IF ( ((INFI(I) > -1).AND.(R(I,I) > LTOL))
     &        .OR. isXd) THEN 
            DO J = 1,I-1
               isXt = (INDEX1(J)<=Nt)
               IF ( (INFI(J)  <  0 .OR. (R(J,J)<= LTOL)) 
     &              .AND. (isXt)) THEN
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
      
      IF ( N1 <= 0 ) GOTO 200
!    
!     Sort remaining limits and determine Cholesky factor.
!     
      Y(1:N1) = gZERO
      K       = 1
      Ndleft  = Nd - INFISD
      Nullity = 0
      DO  WHILE (K .LE. N1) 
     
!     IF (Ndim.EQ.3) EPSL = MAX(EPS2,1D-10)
!     Determine the integration limits for variable with minimum
!     expected probability and interchange that variable with Kth.
     
         K0     = K - Nullity
         PRBMIN = gTWO
         JMIN   = K
         CVDIAG = ZERO
         RMAX   = ZERO
         IF ((Ndleft>0) .OR. (NDIM < Nd+mNIT)) THEN
            DO J = K, N1
               isXd = (INDEX1(J)>Nt)
               isOK = ((NDIM <= Nd+mNIT).OR.isXd)
               IF ( R(J,J) <= K0*K0*EPSL .OR. (.NOT. isOK)) THEN
                  RMAX = max(RMAX,ABS(R(J,J)))
               ELSE
                  TMP = ZERO    ! =  conditional mean of Y(I) given Y(1:I-1)
                  DO I = 1, K0 - 1
                     TMP = TMP + R(I,J)*Y(I)
                  END DO
                  SUMSQ = SQRT( R(J,J))
                  
                  IF (INFI(J)>-1) THEN
                                ! May have infinite int. limits if Nd>0
                     IF (INFI(J).NE.0) THEN
                        AJ = ( A(J) - TMP )/SUMSQ
                     ENDIF
                     IF (INFI(J).NE.1) THEN
                        BJ = ( B(J) - TMP )/SUMSQ
                     ENDIF
                  ENDIF
                  IF (isXd) THEN
                     AA = (Cm(J)+TMP)/SUMSQ ! inflection point
                     CALL EXLMS(AA,AJ,BJ,INFI(J),D,E,Ca,Pa)
                     PRBJ = E - D
                  ELSE   
                                !CALL MVNLMS( AJ, BJ, INFI(J), D, E )
                     CALL MVNLIMITS(AJ,BJ,INFI(J),APJ,PRBJ)                    
                  ENDIF
                                !IF ( EMIN + D .GE. E + DMIN ) THEN
                  IF ( PRBJ  <  PRBMIN ) THEN
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
            isXd = (INDEX1(JMIN)>Nt)
            IF (isXd) THEN
               Ndleft = Ndleft - 1             
            ELSEIF (BCVSRT.EQ..FALSE..AND.(PRBMIN+LTOL>=gONE)) THEN 
!BCVSRT.EQ.
               J = 1
               AJ = R(J,JMIN)*xCut
               BJ = - (R(J,JMIN))*xCut
               if (INFI(J)>=0) then
                  if (INFI(J).ne.0) then
                     AJ = -(R(J,JMIN))*MAX(A(J),-xCut)
                  endif
                  if (INFI(J).ne.1) then
                     BJ = - (R(J,JMIN))*MIN(B(J),xCut)
                  endif
               endif
               if (R(J,JMIN)<gZERO) THEN
                  CALL SWAP(AJ,BJ)
               endif
               INFJ = INFI(JMIN)
               AJ   = A(JMIN)+AJ
               BJ   = B(JMIN)+BJ
               
               D = gZERO
               DO J = 2, K0-1
                  D = D + ABS(R(J,JMIN))
               END DO
               
               AJ = (AJ + D*xCut)/CVDIAG
               BJ = (BJ - D*xCut)/CVDIAG
               CALL ADJLIMITS(AJ,BJ,INFJ)
               IF (INFJ < 0) THEN
                  !variable is redundnant
                  !                     CALL mexPrintf('AdjLim'//CHAR(10))
                  IF ( JMIN < N1 ) THEN
                   CALL RCSWAP( JMIN,N1,N1,N,R,INDEX1,Cm,A,B,INFI)
                     ! move conditional standarddeviations
                     R(JMIN,1:K0-1) = R(N1,1:K0-1)
  
                     Y(JMIN) = Y(N1)                    
                  ENDIF
                  R(1:N1,N1)     = gZERO
                  R(N1,1:N1)     = gZERO
                  Y(N1)   = gZERO
                  INFIS = INFIS+1
                  N1    = N1-1
                  GOTO 100 
               END IF
            ENDIF
            NDIM = NDIM + 1     !Number of relevant dimensions to integrate

            IF ( K < JMIN ) THEN
               CALL RCSWAP( K, JMIN, N1,N, R,INDEX1,Cm, A, B, INFI)
               ! SWAP conditional standarddeviations
               DO J = 1,K0-1  !MIN(K0, K-1)
                  CALL SWAP(R(K,J),R(JMIN,J))
               END DO
            END IF
                            
            R(K0,K) = CVDIAG 
            CDI(K)  = CVDIAG     ! Store the diagonal element
            DO I = K0+1,K
               R(I,K) = gZERO;
               R(K,I) = gZERO
            END DO

            K1 = K	
            I  = K1 + 1
            DO WHILE (I <= N1)              
               TMP = ZERO
               DO J = 1, K0 - 1
                  !tmp = tmp + L(i,j).*L(k1,j)
                  TMP = TMP + R(J,I)*R(J,K1) 
               END DO
                  ! Cov(Xk,Xi|X1,X2,...Xk-1)/STD(Xk|X1,X2,...Xk-1)
               R(K0,I)  = (R(K1,I) - TMP) /CVDIAG  
                  ! Var(Xi|X1,X2,...Xk)
               R(I,I) = R(I,I) - R(K0,I) * R(K0,I)

               IF (R(I,I).GT.LTOL) THEN
                  R(I,K0) = SQRT(R(I,I)) ! STD(Xi|X1,X2,...Xk)
               ELSE   !!IF (R(I,I) .LE. LTOL) THEN !TOL
                                !CALL mexprintf('Singular')
                  isXd = (index1(I)>Nt)
                  if (isXd) then
                     Ndleft = Ndleft - 1  
                  ELSEIF (BCVSRT.EQ..FALSE.) THEN 
!     BCVSRT.EQ.
                     J = 1
                     AJ = R(J,I)*xCut
                     BJ = - (R(J,I))*xCut
                     if (INFI(J)>=0) then
                        if (INFI(J).ne.0) then
                           AJ = -(R(J,I))*MAX(A(J),-xCut)
                        endif
                        if (INFI(J).ne.1) then
                           BJ = - (R(J,I))*MIN(B(J),xCut)
                        endif
                     endif
                     if (R(J,I)<gZERO) THEN
                        CALL SWAP(AJ,BJ)
                     endif
                     INFJ = INFI(I)
                     AJ   = A(I)+AJ
                     BJ   = B(I)+BJ
               
                     D = gZERO
                     DO J = 2, K0
                        D = D + ABS(R(J,I))
                     END DO
               
                     AJ = (AJ + D*xCut)-mXcutOff
                     BJ = (BJ - D*xCut)+mXcutOff
                     !call printvar(Aj,'Aj')
                     !call printvar(Bj,'Bj')
                     CALL ADJLIMITS(AJ,BJ,INFJ)
                     IF (INFJ < 0) THEN
                                !variable is redundnant
                        !CALL mexPrintf('AdjLim'//CHAR(10))
                        IF ( I < N1 ) THEN
                           CALL RCSWAP( I,N1,N1,N,R,INDEX1,Cm,A,B,INFI)
                                ! move conditional standarddeviations
                           R(I,1:K0-1) = R(N1,1:K0-1)
  
                           Y(I) = Y(N1)                    
                        ENDIF
                        R(1:N1,N1)     = gZERO
                        R(N1,1:N1)     = gZERO
                        Y(N1)   = gZERO
                        INFIS = INFIS+1
                        N1    = N1-1
                        
                        !CALL mexprintf('covsrt updated N1')
                        !call printvar(INFIS,' Infis')
                        GOTO 75 
                     END IF
                  END IF
                  IF (mNIT>100) THEN
                     R(I,K0) = gZERO
                  ELSE
                     R(I,K0) = MAX(SQRT(MAX(R(I,I), gZERO)),LTOL)
                  ENDIF
                  Nullity = Nullity + 1
                  K  = K + 1
                  IF (K  <  I) THEN
                     CALL RCSWAP( K, I, N1,N,R,INDEX1,Cm, A, B, INFI)
                     ! SWAP conditional standarddeviations
                     DO J = 1, K0
                        CALL SWAP(R(K,J),R(I,J))
                     END DO
                  ENDIF
                  chkLim = .FALSE. !((.not.isXd).AND.(BCVSRT.EQ..FALSE.))
                  L = INFIS
                  CALL  RCSCALE(chkLim,K,K0,N1,N,K1,INFIS,CDI,Cm,
     &                 R,A,B,INFI,INDEX1)
                  if (L.ne.INFIS) THEN
                     K = K - 1
                     I = I - 1 
                  ENDIF
               END IF	
               I = I + 1
 75            CONTINUE
            END DO
            INFJ = INFI(K1)

            IF (K1 .EQ.1) THEN
               FINA = 0
               FINB = 0
               IF (INFJ.GE.0) THEN
                  IF  (INFJ.NE.0) FINA = 1
                  IF  (INFJ.NE.1) FINB = 1
               ENDIF
               CALL C1C2(K1+1,N1,K0,A,B,INFI, Y, R, 
     &              AMIN, BMIN, FINA,FINB)
               INFJ = 2*FINA+FINB-1
               CALL MVNLIMITS(AMIN,BMIN,INFJ,APJ,PRBMIN) 
            ENDIF
            
            Y(K0) = gettmean(AMIN,BMIN,INFJ,PRBMIN)

           
            R( K0, K1 ) = R( K0, K1 ) / CVDIAG 
            DO J = 1, K0 - 1
               ! conditional covariances
               R( J, K1 ) = R( J, K1 ) / CVDIAG 
               ! conditional standard dev.s used in regression eq.
               R( K1, J ) = R( K1, J ) / CVDIAG 
            END DO
            
            A( K1 ) = A( K1 )/CVDIAG
            B( K1 ) = B( K1 )/CVDIAG
           
            K  = K  + 1
100         CONTINUE
         ELSE
            covErr = RMAX
            R(K:N1,K:N1) = gZERO
            I = K
            DO WHILE (I <= N1)
!  Scale  covariance matrix rows and limits
!  If the conditional covariance matrix diagonal entry is zero, 
!  permute limits and/or rows, if necessary.
               chkLim = ((index1(I)<=Nt).AND.(BCVSRT.EQ..FALSE.))
               L = INFIS
               CALL RCSCALE(chkLim,I,K0-1,N1,N,K1,INFIS,CDI,Cm,
     &              R,A,B,INFI,INDEX1)
               if (L.EQ.INFIS) I = I + 1
            END DO
            Nullity = N1 - K0 + 1
            GOTO 200  !RETURN	
         END IF
      END DO 
 200  CONTINUE
      IF (Ndim .GT. 0) THEN  ! N1<K
         ! K1 = index to the last stochastic varible to integrate
         ! If last stoch. variable is Xt: reduce dimension of integral by 1
         IF (ALL(INDEX1(K1:N1).LE.Nt)) Ndim = Ndim-1
      ENDIF
!      CALL mexprintf('After sorting')
      
!      CALL PRINTCOF(N,A,B,INFI,R,INDEX1)
!      CALL printvar(A(1),'A1')
!      CALL printvar(B(1),'B1')
!      CALL printvar(INFIS,'INFIS')
!      CALL PRINTVEC(CDI,'CDI')
!      CALL PRINTVEC(Y,'Y')
!      CALL PRINTVEC(AA1,'AA1')
!      CALL PRINTVEC(BB1,'BB1')
!      CALL PRINTVAR(NDIM,TXT='NDIM')
!      CALL PRINTVAR(NIT,TXT='NIT')
!      DEALLOCATE(AA1)
!      DEALLOCATE(BB1)
      RETURN
      END SUBROUTINE COVSRT

      SUBROUTINE COVSRT1(BCVSRT, Nt,Nd,R,Cm,A,B,INFI,INDEX1, 
     &     INFIS,INFISD, NDIM, Y, CDI )
      USE FIMOD
      USE SWAPMOD
!      USE GLOBALCONST
!      USE GLOBALDATA, ONLY : EPS2,NIT,xCutOff,Nc1c2
      IMPLICIT NONE
!COVSRT1  sort integration limits and determine Cholesky factor.
!
!     Nt, Nd = size info about Xt and Xd variables.
!     R      = Covariance/Cholesky factored matrix for [Xt,Xd,Xc] (in)
!              On input: 
!              Note: Only upper triangular part is needed/used.
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
!                  ....
!                  etc.
!               3b) The lower triangular part of R contains the
!               conditional standard deviations, i.e., 
!               R(2:N,1)  = [STD(X2|X1) STD(X3|X1),....,STD(XN|X1) ] 
!               R(3:N,2)  = [STD(X3|X1,X2),....,STD(XN|X1,X2) ]  
!               .....
!               etc.
!              
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
      DOUBLE PRECISION :: AA, Ca, Pa, APJ, PRBJ, RMAX, xCut
      DOUBLE PRECISION :: CVDIAG, AMIN, BMIN, PRBMIN, YL, YU
      DOUBLE PRECISION :: SQTWPI,LTOL,TOL,ONE,ZERO,HALF,EIGHT 
!      INTEGER, PARAMETER :: NMAX = 1500
!      DOUBLE PRECISION, DIMENSION(NMAX) :: AP,BP
!      INTEGER, DIMENSION(NMAX) :: INFP
!      INTEGER :: Nabp
      LOGICAL          :: isOK = .TRUE.
      LOGICAL          :: isXd = .FALSE.
      LOGICAL          :: isXt = .FALSE.
      LOGICAL          :: chkLim
      PARAMETER ( SQTWPI = 2.506628274631001D0, TOL = 1D-16 )
      PARAMETER (ONE = 1.D0, ZERO = 0.D0, HALF = 0.5D0,EIGHT = 8.D0) 
      
      xCut = MIN(ABS(mXcutOff),EIGHT)
      INFIS  = 0
      INFISD = 0
      Ndim   = 0
      Ndleft = Nd
      N      = Nt+Nd
      
!     IF (N < 10) EPSL = MIN(1D-10,EPS2)
      LTOL = TOL
      TMP  = ZERO
      DO I = 1, N
         IF (R(I,I).GT.TMP) TMP = R(I,I)
      ENDDO

      EPSL = MAX(mCovEps,N*TMP*mSmall)      
      IF (TMP.GT.EPSL) THEN
         DO I = 1, N
            IF ((INFI(I)  <  0).OR.(R(I,I).LE.LTOL)) THEN
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
      N1  = N-INFIS-INFISD
      CDI(N1+1:N) = ZERO
!     PRINT *,'COVSRT'
!     CALL PRINTCOF(N,A,B,INFI,R,INDEX1)

!     Move any redundant variables of Xd to innermost positions. 
      DO I = N, N-INFISD+1, -1
         isXt = (INDEX1(I).LE.Nt)
         IF ( R(I,I) .GT. LTOL .OR. isXt) THEN 
            DO J = 1,I-1
               isXd = (INDEX1(J).GT.Nt)
               IF ( R(J,J) .LE. LTOL .AND. isXd) THEN
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
      DO I = N-INFISD, N1+1, -1
         isXd = (INDEX1(I).GT.Nt)
         IF ( ((INFI(I) .GE. 0).AND. (R(I,I).GT. LTOL) )
     &        .OR. isXd) THEN 
            DO J = 1,I-1
               isXt = (INDEX1(J).LE.Nt)
               IF ( (INFI(J)  <  0 .OR. (R(J,J).LE. LTOL)) 
     &              .AND. isXt) THEN
                  CALL RCSWAP( J, I, N,N, R,INDEX1,Cm, A, B, INFI)
                  GO TO 15
               ENDIF
            END DO
         ENDIF
 15   END DO
      
      IF ( N1 .LE. 0 ) RETURN
!      CALL mexprintf('Before sorting')
!      CALL PRINTCOF(N,A,B,INFI,R,INDEX1)

     
!     Sort remaining limits and determine Cholesky factor.
      Y(1:N1) = ZERO
      K  = 1
!      N1  = N-INFIS-INFISD
      Ndleft = Nd - INFISD
      Nullity = 0
      
!      Nabp  = 0
!      AP(1:N1) = ZERO
!      BP(1:N1) = zero
      DO  WHILE (K .LE. N1) 
     
!     Determine the integration limits for variable with minimum
!     expected probability and interchange that variable with Kth.
     
         K0     = K-Nullity
         PRBMIN = 2.d0
         JMIN   = K
         CVDIAG = ZERO
         RMAX   = ZERO
         IF (Ndleft.GT.0 .OR. NDIM < Nd+mNIT) THEN
            DO J = K,N1
               isXd = (INDEX1(J).GT.Nt)
               isOK = ((NDIM <= Nd+mNIT).OR.isXd)
               IF ( R(J,J) .LE. K0*K0*EPSL.OR. (.NOT. isOK)) THEN
                  RMAX = max(RMAX,R(J,J))
               ELSE
                  TMP = Y(J) ! =  the conditional mean of Y(J) given Y(1:J-1)
                  SUMSQ = SQRT( R(J,J))
                  
                  IF (INFI(J) < 0) GO TO 30 ! May have infinite int. limits if Nd>0
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
                  IF ( PRBJ  <  PRBMIN ) THEN
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
               IF (BCVSRT.EQ..FALSE..AND.(PRBMIN+LTOL.GE.gONE)) THEN 
!BCVSRT.EQ.
                  I = 1
                  AJ = R(I,JMIN)*xCut
                  BJ = - (R(I,JMIN))*xCut
                  if (INFI(1)>=0) then
                     if (INFI(1).ne.0) then
                        AJ = -(R(I,JMIN))*MAX(A(I),-xCut)
                     endif
                     if (INFI(1).ne.1) then
                        BJ = - (R(I,JMIN))*MIN(B(I),xCut)
                     endif
                  endif
                  if (R(I,JMIN)<gZERO) THEN
                     CALL SWAP(AJ,BJ)
                  endif
                  INFJ = INFI(JMIN)
                  AJ   = A(JMIN)+AJ
                  BJ   = B(JMIN)+BJ

                  D = gZERO
                  DO I = 2, K0-1
                     D = D + ABS(R(I,JMIN))
                  END DO

                  

                  AJ = (AJ + D*xCut)/CVDIAG
                  BJ = (BJ - D*xCut)/CVDIAG
                  CALL ADJLIMITS(AJ,BJ,INFJ)
                  IF (INFJ < 0) THEN
                     !variable is redundnant
                     CALL mexPrintf('AdjLim'//CHAR(10))
                   IF ( JMIN < N1 ) THEN
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
            NDIM = NDIM + 1     !Number of relevant dimensions to integrate

            IF ( K < JMIN ) THEN
                              
               CALL RCSWAP( K, JMIN, N1,N, R,INDEX1,Cm, A, B, INFI)
               ! SWAP conditional standarddeviations
               DO J=1,K0-1
                  CALL SWAP(R(K,J),R(JMIN,J))
               END DO
               CALL SWAP(Y(K),Y(JMIN))
            END IF
            
            
            R(K0,K:N1) = R(K0,K:N1)/CVDIAG
            R(K0,K) = CVDIAG
            CDI(K)  = CVDIAG     ! Store the diagonal element
            DO I = K0+1,K
               R(I,K) = ZERO
               R(K,I) = ZERO
            END DO
            
            K1  = K
            !IF (K .EQ. N1) GOTO 200
            
!  Cov(Xi,Xj|Xk,Xk+1,..,Xn)=
!             Cov(Xi,Xj|Xk+1,..,Xn) -
 !             Cov(Xi,Xk|Xk+1,..Xn)*Cov(Xj,Xk|Xk+1,..Xn)
            I = K1 +1
            DO WHILE (I <= N1)
                ! Var(Xj | Xk,Xk+1,...,Xn)
               R(I,I)  =  R(I,I) - R(K0,I)*R(K0,I)
               IF (R(I,I).GT.LTOL) THEN
                  R(I,K0)     = SQRT(R(I,I)) ! STD(Xi|X1,X2,...Xk)
                  R(I,I+1:N1) = R(I,I+1:N1) - R(K0,I+1:N1)*R(K0,I)
               ELSE
                  R(I,K0) = MAX(SQRT(MAX(R(I,I), gZERO)),LTOL)
                  Nullity = Nullity + 1
                  K  = K + 1
                  IF (K  <  I) THEN
                     CALL RCSWAP( K, I, N1,N,R,INDEX1,Cm, A, B, INFI)
                     ! SWAP conditional standarddeviations
                     DO J=1,K0
                        CALL SWAP(R(K,J),R(I,J))
                     END DO
                     CALL SWAP(Y(K),Y(I))
                  ENDIF    
                  isXd = (INDEX1(K).GT.Nt)
                  IF (isXd) Ndleft = Ndleft-1
                  chkLim = ((.not.isXd).AND.(BCVSRT.EQ..FALSE.))
                  L = INFIS
                  CALL  RCSCALE(chkLim,K,K0,N1,N,K1,INFIS,CDI,Cm,
     &                 R,A,B,INFI,INDEX1,Y)
                  IF (L.NE.INFIS) I = I - 1
               END IF	
               I = I +1
            END DO 
            INFJ = INFI(K1)
            IF (K0 == 1) THEN
               FINA = 0
               FINB = 0
               IF (INFJ.GE.0) THEN
                  IF  (INFJ.NE.0) FINA = 1
                  IF  (INFJ.NE.1) FINB = 1
               ENDIF
               CALL C1C2(K1+1,N1,K0,A,B,INFI, Y, R, 
     &              AMIN, BMIN, FINA,FINB)
               INFJ = 2*FINA+FINB-1
               CALL MVNLIMITS(AMIN,BMIN,INFJ,APJ,PRBMIN) 
            ENDIF
            Y(K0) = GETTMEAN(AMIN,BMIN,INFJ,PRBMIN)

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
            I = K
            DO WHILE (I <= N1)
!  Scale  covariance matrix rows and limits
!  If the conditional covariance matrix diagonal entry is zero, 
!  permute limits and/or rows, if necessary.
               chkLim = ((index1(I)<=Nt).AND.(BCVSRT.EQ..FALSE.))
               L = INFIS
               CALL RCSCALE(chkLim,I,K0-1,N1,N,K1,INFIS,CDI,Cm,
     &              R,A,B,INFI,INDEX1)
               if (L.EQ.INFIS) I = I + 1
            END DO
            Nullity = N1 - K0 + 1
            GOTO 200  !RETURN	
         END IF
      END DO 
      
 200  CONTINUE
      IF (Ndim .GT. 0) THEN     ! N1<K
         ! K1 = index to the last stochastic varible to integrate
         IF (ALL(INDEX1(K1:N1).LE.Nt)) Ndim = Ndim - 1
      ENDIF
!      CALL mexprintf('After sorting')
!      CALL PRINTCOF(N,A,B,INFI,R,INDEX1)
!      CALL PRINTVEC(CDI)
!      CALL PRINTVAR(NDIM,TXT='NDIM')
      RETURN
      END SUBROUTINE COVSRT1
      
      SUBROUTINE RCSWAP( P, Q, N,Ntd, C,IND,Cm, A, B, INFIN )
      USE SWAPMOD
      IMPLICIT NONE
! RCSWAP  Swaps rows and columns P and Q in situ, with P <= Q.
!
!
!   CALL  RCSWAP( P, Q, N, Ntd, C,IND A, B, INFIN, Cm)
!
!    P, Q  = row/column number to swap P<=Q<=N
!    N     = Number of significant variables of [Xt,Xd] 
!    Ntd   = length(Xt)+length(Xd)
!    C     = upper triangular cholesky factor.Cov([Xt,Xd,Xc]) size Ntdc X Ntdc
!    IND   = permutation index vector. (index to variables original place).
!    Cm    = conditional mean
!    A,B   = lower and upper integration limit, respectively.
!    INFIN = if INFIN(I) < 0, Ith limits are (-infinity, infinity);
!            if INFIN(I) = 0, Ith limits are (-infinity, B(I)];
!            if INFIN(I) = 1, Ith limits are [A(I), infinity);
!            if INFIN(I) = 2, Ith limits are [A(I), B(I)].
!
! NOTE: RCSWAP works only on the upper triangular part of C
      DOUBLE PRECISION, DIMENSION(:,:),INTENT(inout) :: C
      INTEGER, DIMENSION(:),INTENT(inout) :: IND
      INTEGER, DIMENSION(:),          OPTIONAL,INTENT(inout) :: INFIN
      DOUBLE PRECISION, DIMENSION(:), OPTIONAL,INTENT(inout) :: A,B,Cm
      INTEGER,INTENT(in) :: P, Q, N, Ntd
! local variable
      INTEGER :: J, Ntdc
      LOGICAL :: isXc
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
      isXc = (N < Ntdc)
      IF (isXc) THEN
         !Swap row P and Q of Xc variables
         DO J = Ntd+1, Ntdc
            CALL SWAP( C(P,J), C(Q,J) )
         END DO
      ENDIF
      RETURN
      END SUBROUTINE RCSWAP
      end module rind
