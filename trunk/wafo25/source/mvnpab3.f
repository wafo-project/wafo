*
* This file contains a short test program, software for MVN distribution 
* subroutines SADMVN, RANMVN, KROMVN and SPHMVN, plus supporting software.
* The file is self contained and should compile without errors on (77) 
* standard Fortran compilers. The test program demonstrates the use of
* four different methods for computing MVN distribution values for a
* five dimensional example problem, with three different integration limit
* combinations.
*
*          Alan Genz
*          Department of Mathematics
*          Washington State University
*          Pullman, WA 99164-3113
*          Email : alangenz@wsu.edu
*
      PROGRAM TSTNRM
*
*     Test program for SADMVN, RANMVN, KROMVN, SPHMVN
*
      DOUBLE PRECISION :: ABSEPS, RELEPS
      DOUBLE PRECISION :: VALR,ERRR,VALS,ERRS,VALK,ERRK,VALD,ERRD
      INTEGER :: N, NN, I, J, K, IJ, MAXPTS, IFTR, IFTS, IFTK, IFTD
      PARAMETER ( N = 5, NN = ( N - 1 )*N/2, MAXPTS = 5000*N*N*N )
      PARAMETER ( ABSEPS = 0.00005, RELEPS = 0 )
      DOUBLE PRECISION, DIMENSION(NN) :: CORREL
      DOUBLE PRECISION, DIMENSION(N)  ::  LOW, UP
      DOUBLE PRECISION :: CONDIT
      INTEGER, DIMENSION(N) ::  INFIN(N)
*          Chen Problem
      DATA ( UP(I),I=1,N) /0.d0, 1.5198d0, 1.7817d0, 1.4755d0, 1.5949d0/
      DATA (LOW(I), I=1,N) /0.d0,  0.d0 , 1.7817d0, 1.4755d0, 1.5949d0/
      DATA (INFIN(I), I=1,N)/ 1, 2     , 1     , 1     , 0     /
      DATA (CORREL(I),I=1,NN)
     &     /-0.707107d0 ,0.d0, 0.5d0, 0.d0, 2*0.5d0, 0.d0, 3*0.5d0/
      PRINT '(''        Test of RANMVN, SADMVN, KROMVN and SPHMVN'')'
      PRINT '(12X, ''Requested Accuracy '',F8.5)', MAX(ABSEPS,RELEPS)
      PRINT '(''           Number of Dimensions is '',I2)', N
      PRINT '(''     Maximum # of Function Values is '',I7)', MAXPTS
*
      DO K = 1,3
         PRINT '(/'' I     Limits'')'
         PRINT '(4X,''Lower  Upper  Lower Left of Correlation Matrix'')'
         IJ = 0
         DO I = 1,N
            IF ( INFIN(I) .LT. 0 ) THEN 
               PRINT '(I2, '' -infin  infin '', 7F9.5)',
     &              I, ( CORREL(IJ+J), J = 1,I-1 ), 1.0
            ELSE IF ( INFIN(I) .EQ. 0 ) THEN 
               PRINT '(I2, '' -infin'', F7.4, 1X, 7F9.5)',
     &              I, UP(I), ( CORREL(IJ+J), J = 1,I-1 ), 1.0
            ELSE IF ( INFIN(I) .EQ. 1 ) THEN 
               PRINT '(I2, F7.4, ''  infin '', 7F9.5)',
     &              I, LOW(I), ( CORREL(IJ+J), J = 1,I-1 ), 1.0
            ELSE 
               PRINT '(I2, 2F7.4, 1X, 7F9.5)', 
     &              I, LOW(I), UP(I), ( CORREL(IJ+J), J = 1,I-1 ), 1.0
            ENDIF
            IJ = IJ + I-1
         END DO
         PRINT '(20X, ''Condition Number is'', F8.1)', CONDIT(N,CORREL) 
         CALL SADMVN( N, LOW, UP, INFIN, CORREL, 
     &                MAXPTS, ABSEPS, RELEPS, ERRS, VALS, IFTS )
         CALL KROMVN( N, LOW, UP, INFIN, CORREL, 
     &                MAXPTS, ABSEPS, RELEPS, ERRK, VALK, IFTK )
         ERRR=0.d0
         ERRD = 0.d0
         VALR=0.d0
         IFTR=0
         VALD=0.d0
         IFTD=0
!         CALL RANMVN( N, LOW, UP, INFIN, CORREL, 
!     &                MAXPTS, ABSEPS, RELEPS, ERRR, VALR, IFTR )
!         CALL SPHMVN( N, LOW, UP, INFIN, CORREL, 
!     &                MAXPTS, ABSEPS, RELEPS, ERRD, VALD, IFTD )
       PRINT '('' Results for:  RANMVN'',9X,''SADMVN'',9X,''KROMVN'', 
     &                      9X, ''SPHMVN'')'
         PRINT '('' Values:   '',4(F11.6,I4))',
     &        VALR, IFTR, VALS, IFTS, VALK, IFTK, VALD, IFTD
         PRINT '(''Errests:   '',4(2X,''('',F8.6'')'',3X))', 
     &        ERRR, ERRS, ERRK, ERRD
         INFIN(1) = INFIN(1) - 1
      END DO
      END PROGRAM TSTNRM
!      CONTAINS

      DOUBLE PRECISION FUNCTION UNI()
*
*     Uniform (0, 1) random number generator
*
*     Reference:
*     L'Ecuyer, Pierre (1996), 
*     "Combined Multiple Recursive Random Number Generators"
*     Operations Research 44, pp. 816-822.
*
*
      INTEGER :: A12, A13, A21, A23, P12, P13, P21, P23
      INTEGER :: Q12, Q13, Q21, Q23, R12, R13, R21, R23
      INTEGER :: X10, X11, X12, X20, X21, X22, Z, M1, M2, H 
      DOUBLE PRECISION :: INVMP1
      PARAMETER (  M1 = 2147483647,  M2 = 2145483479 )
      PARAMETER ( A12 =   63308,    Q12 = 33921, R12 = 12979 )
      PARAMETER ( A13 = -183326,    Q13 = 11714, R13 =  2883 )
      PARAMETER ( A21 =   86098,    Q21 = 24919, R21 =  7417 )
      PARAMETER ( A23 = -539608,    Q23 =  3976, R23 =  2071 )
      PARAMETER ( INVMP1 = 4.656612873077392578125D-10 ) 
*                 INVMP1 = 1/(M1+1)
      SAVE X10, X11, X12, X20, X21, X22
      DATA       X10,      X11,      X12,      X20,      X21,      X22  
     &    / 11111111, 22222223, 33333335, 44444447, 55555559, 66666661 /
*
*     Component 1
*
      H = X10/Q13
      P13 = -A13*( X10 - H*Q13 ) - H*R13
      H = X11/Q12
      P12 =  A12*( X11 - H*Q12 ) - H*R12
      IF ( P13 .LT. 0 ) P13 = P13 + M1
      IF ( P12 .LT. 0 ) P12 = P12 + M1
      X10 = X11 
      X11 = X12
      X12 = P12 - P13
      IF ( X12 .LT. 0 ) X12 = X12 + M1
*
*     Component 2
*
      H = X20/Q23
      P23 = -A23*( X20 - H*Q23 ) - H*R23
      H = X22/Q21
      P21 =  A21*( X22 - H*Q21 ) - H*R21
      IF ( P23 .LT. 0 ) P23 = P23 + M2
      IF ( P21 .LT. 0 ) P21 = P21 + M2
      X20 = X21 
      X21 = X22
      X22 = P21 - P23
      IF ( X22 .LT. 0 ) X22 = X22 + M2
*
*     Combination
*
      Z = X12 - X22
      IF ( Z .LE. 0 ) Z = Z + M1
      UNI = Z*INVMP1
      RETURN
      END FUNCTION UNI
      SUBROUTINE SADMVN( N, LOWER, UPPER, INFIN, CORREL, MAXPTS,
     &                   ABSEPS, RELEPS, ERROR, VALUE, INFORM )
*
*     A subroutine for computing multivariate normal probabilities.
*     This subroutine uses an algorithm given in the paper
*     "Numerical Computation of Multivariate Normal Probabilities", in
*     J. of Computational and Graphical Stat., 1(1992), pp. 141-149, by
*          Alan Genz 
*          Department of Mathematics
*          Washington State University 
*          Pullman, WA 99164-3113
*          Email : alangenz@wsu.edu
*
*  Parameters
*
*     N      INTEGER, the number of variables.
*     LOWER  REAL, array of lower integration limits.
*     UPPER  REAL, array of upper integration limits.
*     INFIN  INTEGER, array of integration limits flags:
*            if INFIN(I) < 0, Ith limits are (-infinity, infinity);
*            if INFIN(I) = 0, Ith limits are (-infinity, UPPER(I)];
*            if INFIN(I) = 1, Ith limits are [LOWER(I), infinity);
*            if INFIN(I) = 2, Ith limits are [LOWER(I), UPPER(I)].
*     CORREL REAL, array of correlation coefficients; the correlation
*            coefficient in row I column J of the correlation matrix
*            should be stored in CORREL( J + ((I-2)*(I-1))/2 ), for J < I.
*     MAXPTS INTEGER, maximum number of function values allowed. This 
*            parameter can be used to limit the time taken. A 
*            sensible strategy is to start with MAXPTS = 1000*N, and then
*            increase MAXPTS if ERROR is too large.
*     ABSEPS REAL absolute error tolerance.
*     RELEPS REAL relative error tolerance.
*     ERROR  REAL estimated absolute error, with 99% confidence level.
*     VALUE  REAL estimated value for the integral
*     INFORM INTEGER, termination status parameter:
*            if INFORM = 0, normal completion with ERROR < EPS;
*            if INFORM = 1, completion with ERROR > EPS and MAXPTS 
*                           function vaules used; increase MAXPTS to 
*                           decrease ERROR;
*            if INFORM = 2, N > 20 or N < 1.
*
      USE ADAPTMOD
      EXTERNAL MVNFNC
      INTEGER :: N, NL, M, LENWRK, MAXPTS, INFORM, INFIS,
     &     RULCLS, TOTCLS, NEWCLS, MAXCLS
      INTEGER, DIMENSION(:) :: INFIN
      DOUBLE PRECISION, DIMENSION(:) ::  CORREL, LOWER, UPPER
      DOUBLE PRECISION :: ABSEPS, RELEPS, ERROR, VALUE,
     &     OLDVAL, D, E, MVNNIT, MVNFNC
      PARAMETER ( NL = 20 )
      PARAMETER ( LENWRK = 20*NL**2 )
      DOUBLE PRECISION, DIMENSION(LENWRK) :: WORK
      IF ( N .GT. 20 .OR. N .LT. 1 ) THEN
         INFORM = 2
         VALUE = 0.d0
         ERROR = 1.d0
         RETURN
      ENDIF
      INFORM = MVNNIT( N, CORREL, LOWER, UPPER, INFIN, INFIS, D, E )
      M = N - INFIS
      IF ( M .EQ. 0 ) THEN
         VALUE = 1.d0
         ERROR = 0.d0 
      ELSE IF ( M .EQ. 1 ) THEN
         VALUE = E - D
         ERROR = 2E-16
      ELSE
*
*        Call the subregion adaptive integration subroutine
*
         M = M - 1
         RULCLS = 1
         CALL ADAPT( M, RULCLS, 0, MVNFNC, ABSEPS, RELEPS, 
     &               LENWRK, WORK, ERROR, VALUE, INFORM )
         MAXCLS = MIN( 10*RULCLS, MAXPTS )
         TOTCLS = 0
         CALL ADAPT(M, TOTCLS, MAXCLS, MVNFNC, ABSEPS, RELEPS, 
     &        LENWRK, WORK, ERROR, VALUE, INFORM)
         IF ( ERROR .GT. MAX( ABSEPS, RELEPS*ABS(VALUE) ) ) THEN
 10         OLDVAL = VALUE
            MAXCLS = MAX( 2*RULCLS, 
     &           MIN( int(3*MAXCLS/2), MAXPTS - TOTCLS ) )
            NEWCLS = -1
            CALL ADAPT(M, NEWCLS, MAXCLS, MVNFNC, ABSEPS, RELEPS, 
     &           LENWRK, WORK, ERROR, VALUE, INFORM)
            TOTCLS = TOTCLS + NEWCLS
            ERROR = ABS(VALUE-OLDVAL) + SQRT(RULCLS*ERROR**2/TOTCLS)
            IF ( ERROR .GT. MAX( ABSEPS, RELEPS*ABS(VALUE) ) ) THEN
               IF ( MAXPTS - TOTCLS .GT. 2*RULCLS ) GO TO 10
            ELSE 
               INFORM = 0
            END IF
         ENDIF
      ENDIF
      RETURN
      END SUBROUTINE SADMVN
*
      SUBROUTINE KROMVN( N, LOWER, UPPER, INFIN, CORREL, MAXPTS,
     &                   ABSEPS, RELEPS, ERROR, VALUE, INFORM )
*
*     A subroutine for computing multivariate normal probabilities.
*     This subroutine uses an algorithm given in the paper
*     "Numerical Computation of Multivariate Normal Probabilities", in
*     J. of Computational and Graphical Stat., 1(1992), pp. 141-149, by
*          Alan Genz 
*          Department of Mathematics
*          Washington State University 
*          Pullman, WA 99164-3113
*          Email : AlanGenz@wsu.edu
*
*  Parameters
*
*     N      INTEGER, the number of variables.
*     LOWER  REAL, array of lower integration limits.
*     UPPER  REAL, array of upper integration limits.
*     INFIN  INTEGER, array of integration limits flags:
*            if INFIN(I) < 0, Ith limits are (-infinity, infinity);
*            if INFIN(I) = 0, Ith limits are (-infinity, UPPER(I)];
*            if INFIN(I) = 1, Ith limits are [LOWER(I), infinity);
*            if INFIN(I) = 2, Ith limits are [LOWER(I), UPPER(I)].
*     CORREL REAL, array of correlation coefficients; the correlation
*            coefficient in row I column J of the correlation matrix
*            should be stored in CORREL( J + ((I-2)*(I-1))/2 ), for J < I.
*     MAXPTS INTEGER, maximum number of function values allowed. This 
*            parameter can be used to limit the time. A sensible 
*            strategy is to start with MAXPTS = 1000*N, and then
*            increase MAXPTS if ERROR is too large.
*     ABSEPS REAL absolute error tolerance.
*     RELEPS REAL relative error tolerance.
*     ERROR  REAL estimated absolute error, with 99% confidence level.
*     VALUE  REAL estimated value for the integral
*     INFORM INTEGER, termination status parameter:
*            if INFORM = 0, normal completion with ERROR < EPS;
*            if INFORM = 1, completion with ERROR > EPS and MAXPTS 
*                           function vaules used; increase MAXPTS to 
*                           decrease ERROR;
*            if INFORM = 2, N > 100 or N < 1.
*
      USE KROBOVMOD
      EXTERNAL MVNFNC
      INTEGER :: N, MAXPTS, INFORM, INFIS, IVLS
      INTEGER, DIMENSION(:) :: INFIN
      DOUBLE PRECISION, DIMENSION(:) :: CORREL, LOWER, UPPER
      DOUBLE PRECISION :: RELEPS, ABSEPS, ERROR, VALUE, 
     &     E, D, MVNNIT, MVNFNC
      IF ( N .GT. 100 .OR. N .LT. 1 ) THEN
         INFORM = 2
         VALUE = 0.d0
         ERROR = 1.d0
      ELSE
         INFORM = MVNNIT(N, CORREL, LOWER, UPPER, INFIN, INFIS, D, E)
         IF ( N-INFIS .EQ. 0 ) THEN
            VALUE = 1.d0
            ERROR = 0.d0
         ELSE IF ( N-INFIS .EQ. 1 ) THEN
            VALUE = E - D
            ERROR = 2E-16
         ELSE
*
*        Call the lattice rule integration subroutine
*
            IVLS = 0
            CALL KROBOV( N-INFIS-1, IVLS, MAXPTS, MVNFNC, 
     &                   ABSEPS, RELEPS, ERROR, VALUE, INFORM )
         ENDIF
      ENDIF
      RETURN
      END SUBROUTINE KROMVN

*
      DOUBLE PRECISION FUNCTION MVNFNC(N, W)
*     
*     Integrand subroutine
*
      INTEGER :: N, INFIS
      INTEGER, DIMENSION(:) :: INFIN
      DOUBLE PRECISION, DIMENSION(:) :: W, LOWER, UPPER, CORREL
      DOUBLE PRECISION :: ONE
      INTEGER :: NL, IJ, I, J
      PARAMETER ( NL = 100, ONE = 1.d0 )
      DOUBLE PRECISION, DIMENSION((NL*(NL+1))/2) :: COV
      DOUBLE PRECISION, DIMENSION(NL) :: A, B, Y
      INTEGER, DIMENSION(NL) :: INFI
      DOUBLE PRECISION :: PROD, D1, E1, DI, EI, SUM1, PHINV, 
     &     D, E, MVNNIT, BVN
      SAVE D1, E1, A, B, INFI, COV
      DI = D1
      EI = E1
      PROD = E1 - D1
      IJ = 1
      DO I = 1,N
         Y(I) = PHINV( DI + W(I)*(EI-DI) )
         SUM1 = 0.d0
         DO J = 1,I
            IJ = IJ + 1
            SUM1 = SUM1 + COV(IJ)*Y(J)
         END DO
         IJ = IJ + 1
         IF ( COV(IJ) .GT. 0 ) THEN
            CALL LIMITS( A(I+1)-SUM1, B(I+1)-SUM1, INFI(I+1), DI, EI )
         ELSE
            DI = ( 1.d0 + SIGN( ONE, A(I+1)-SUM1 ) )/2.d0
            EI = ( 1.d0 + SIGN( ONE, B(I+1)-SUM1 ) )/2.d0
         ENDIF
         PROD = PROD*(EI-DI)
      END DO
      MVNFNC = PROD
      RETURN
*
*     Entry point for intialization.
*
      ENTRY MVNNIT(N, CORREL, LOWER, UPPER, INFIN, INFIS, D, E)
      MVNNIT = 0
*
*     Initialization and computation of covariance Cholesky factor.
*
      CALL NCVSRT(N, LOWER,UPPER,CORREL,INFIN,Y, INFIS,A,B,INFI,COV,D,E)
      D1 = D
      E1 = E
      IF ( N - INFIS .EQ. 2 ) THEN
         D = SQRT( 1.d0 + COV(2)**2 )
         A(2) = A(2)/D
         B(2) = B(2)/D
         E = BVN( A, B, INFI, COV(2)/D )
         D = 0.d0
         INFIS = INFIS + 1 
      END IF
      END  FUNCTION MVNFNC
      SUBROUTINE LIMITS( A, B, INFIN, LOWER, UPPER )
      DOUBLE PRECISION :: A, B, LOWER, UPPER, PHI
      INTEGER :: INFIN
      LOWER = 0.d0
      UPPER = 1.d0
      IF ( INFIN .GE. 0 ) THEN
         IF ( INFIN .NE. 0 ) LOWER = PHI(A)
         IF ( INFIN .NE. 1 ) UPPER = PHI(B)
      ENDIF
      END  SUBROUTINE LIMITS     
      SUBROUTINE NCVSRT( N, LOWER, UPPER, CORREL, INFIN, Y, INFIS, 
     &                   A, B, INFI, COV, D, E )
*
*     Subroutine to sort integration limits.
*
      INTEGER N, INFI(*), INFIN(*), INFIS
      DOUBLE PRECISION 
     &     A(*), B(*), COV(*), LOWER(*), UPPER(*), CORREL(*), Y(*), D, E
      INTEGER I, J, K, IJ, II, JMIN
      DOUBLE PRECISION SUMSQ, ZERO
      PARAMETER ( ZERO = 0 )
      DOUBLE PRECISION AJ, BJ, SUM1, SQTWPI
      DOUBLE PRECISION CVDIAG, AMIN, BMIN, DMIN, EMIN, YL, YU
      PARAMETER ( SQTWPI = 2.50662 82746 31000 50240 )
      IJ = 0
      II = 0
      INFIS = 0
      DO I = 1,N
         INFI(I) = INFIN(I) 
         IF ( INFI(I) .LT. 0 ) THEN
            INFIS = INFIS + 1
         ELSE 
            A(I) = 0.d0
            B(I) = 0.d0
            IF ( INFI(I) .NE. 0 ) A(I) = LOWER(I)
            IF ( INFI(I) .NE. 1 ) B(I) = UPPER(I)
         ENDIF
         DO J = 1,I-1
            IJ = IJ + 1
            II = II + 1
            COV(IJ) = CORREL(II)
         END DO
         IJ = IJ + 1
         COV(IJ) = 1
      END DO
*
*     First move any doubly infinite limits to innermost positions
*
      IF ( INFIS .LT. N ) THEN
         DO I = N,N-INFIS+1,-1
            IF ( INFI(I) .GE. 0 ) THEN 
               DO J = 1,I-1
                  IF ( INFI(J) .LT. 0 ) THEN
                     CALL RCSWAP(J, I, A, B, INFI, N, COV)
                     GO TO 10
                  ENDIF
               END DO
            ENDIF
 10      END DO
*
*     Sort remaining limits and determine Cholesky decomposition
*
         II = 0
         DO I = 1,N-INFIS
*
*     Determine the integration limits for variable with minimum
*      expected probability and interchange that variable with Ith.
*
            EMIN = 1.d0
            DMIN = 0.d0
            JMIN = I
            CVDIAG = 0.d0
            IJ = II
            DO J = I, N-INFIS
               SUM1 = 0.d0
               SUMSQ = 0.d0
               DO K = 1, I-1
                  SUM1 = SUM1 + COV(IJ+K)*Y(K)
                  SUMSQ = SUMSQ + COV(IJ+K)**2
               END DO
               IJ = IJ + J 
               SUMSQ = SQRT( MAX( COV(IJ)-SUMSQ, ZERO ) )
               IF ( SUMSQ .GT. 0 ) THEN
                  IF ( INFI(J) .NE. 0 ) AJ = ( A(J) - SUM1 )/SUMSQ
                  IF ( INFI(J) .NE. 1 ) BJ = ( B(J) - SUM1 )/SUMSQ
                  CALL LIMITS( AJ, BJ, INFI(J), D, E )
                  IF ( EMIN - DMIN .GE. E - D ) THEN
                     JMIN = J
                     IF ( INFI(J) .NE. 0 ) AMIN = AJ 
                     IF ( INFI(J) .NE. 1 ) BMIN = BJ
                     DMIN = D
                     EMIN = E
                     CVDIAG = SUMSQ
                  ENDIF
               ENDIF
            END DO
            IF ( JMIN .NE. I) CALL RCSWAP(I, JMIN, A, B, INFI, N, COV)
*
*     Compute Ith column of Cholesky factor.
*
            IJ = II + I
            COV(IJ) = CVDIAG
            DO J = I+1, N-INFIS               
               IF ( CVDIAG .GT. 0 ) THEN
                  SUM1 = COV(IJ+I)
                  DO K = 1, I-1
                     SUM1 = SUM1 - COV(II+K)*COV(IJ+K)
                  END DO
                  COV(IJ+I) = SUM1/CVDIAG
               ELSE
                  COV(IJ+I) = 0
               ENDIF
               IJ = IJ + J
            END DO
*
*     Compute expected value for Ith integration variable and
*     scale Ith covariance matrix row and limits.
*
            IF ( CVDIAG .GT. 0 ) THEN
               IF ( EMIN .GT. DMIN + 1D-8 ) THEN
                  YL = 0
                  YU = 0
                  IF ( INFI(I) .NE. 0 ) YL = -EXP( -AMIN**2/2 )/SQTWPI
                  IF ( INFI(I) .NE. 1 ) YU = -EXP( -BMIN**2/2 )/SQTWPI
                  Y(I) = ( YU - YL )/( EMIN - DMIN )
               ELSE
                  IF ( INFI(I) .EQ. 0 ) Y(I) = BMIN
                  IF ( INFI(I) .EQ. 1 ) Y(I) = AMIN
                  IF ( INFI(I) .EQ. 2 ) Y(I) = ( AMIN + BMIN )/2
               END IF
               DO J = 1,I
                  II = II + 1
                  COV(II) = COV(II)/CVDIAG
               END DO
               IF ( INFI(I) .NE. 0 ) A(I) = A(I)/CVDIAG
               IF ( INFI(I) .NE. 1 ) B(I) = B(I)/CVDIAG
            ELSE
               Y(I) = 0.d0
               II = II + I
            ENDIF
         END DO
         CALL LIMITS( A(1), B(1), INFI(1), D, E)
      ENDIF
      RETURN
      END SUBROUTINE NCVSRT
      DOUBLE PRECISION FUNCTION CONDIT( N, SYMIN )
*
*     Computes condition number of symmetric matix in situ
*
      INTEGER NL, N
      PARAMETER ( NL = 100 )
      DOUBLE PRECISION DET, SYMIN(*), SUM1, ROWMX, ROWMXI,
     & SYM(NL*(NL+1)/2)
      INTEGER II, IJ, I, J, IM
      ROWMX = 0
      IJ = 0
      DO I = 1,N
         SUM1 = 0
         IM = (I-2)*(I-1)/2
         DO J = 1,I-1
            IM = IM + 1
            SUM1 = SUM1 + ABS(SYMIN(IM))
            IJ = IJ + 1
            SYM(IJ) = SYMIN(IM)
         END DO
         SUM1 = SUM1 + 1
         IJ = IJ + 1
         SYM(IJ) = 1
         IM = IM + I
         DO J = I,N-1
            SUM1 = SUM1 + ABS(SYMIN(IM))
            IM = IM + J
         END DO
         ROWMX = MAX( SUM1, ROWMX )
      END DO
      CALL SYMINV(N, SYM, DET)
      ROWMXI = 0
      II = 0
      DO I = 1,N
         SUM1 = 0.d0
         IJ = II
         DO J = 1,I
            IJ = IJ + 1
            SUM1 = SUM1 + ABS(SYM(IJ))
         END DO
         DO J = I,N-1
            IJ = IJ + J
            SUM1 = SUM1 + ABS(SYM(IJ))
         END DO
         ROWMXI = MAX( SUM1, ROWMXI )
         II = II + I
      END DO
      CONDIT = ROWMX*ROWMXI
      RETURN
      END FUNCTION CONDIT
      SUBROUTINE SYMINV(N, LOWINV, DET)
*
*     Computes lower symmetric inverse and determinant in situ
*
      INTEGER I, II, N
      DOUBLE PRECISION LOWINV(*), DET
      CALL CHOLSK(N, LOWINV)
      DET = 1
      II = 0
      DO I = 1,N
         II = II + I
         DET = DET*LOWINV(II)
      END DO
      DET = DET*DET
      CALL CHOLNV(N, LOWINV)
      CALL CHOLPI(N, LOWINV)
      END SUBROUTINE SYMINV
      SUBROUTINE CHOLSK(N, CHOFAC)
*
*     Computes Choleski factor in situ
*
      INTEGER I, II, J, JJ, K, N
      DOUBLE PRECISION CHOFAC(*), T
      DOUBLE PRECISION S, ZERO
      PARAMETER ( ZERO = 0.d0 )
      JJ = 0
      DO J = 1,N
         II = JJ
         DO I = J,N
            S = CHOFAC(II+J)
            DO K = 1,J-1
               S = S - CHOFAC(II+K)*CHOFAC(JJ+K)
            END DO
            IF ( I .EQ. J ) THEN
               T = SQRT( MAX( S, ZERO ) )
               CHOFAC(II+J) = T
            ELSE
               CHOFAC(II+J) = S/T
            ENDIF
            II = II + I
         END DO
         JJ = JJ + J
      END DO
      RETURN
      END SUBROUTINE CHOLSK
      SUBROUTINE CHOLNV(N, CHOINV)
*
*     Inverts a lower triangular matrix in situ
*
      INTEGER I, II, J, JJ, K, KK, N
      DOUBLE PRECISION CHOINV(*), T
      DOUBLE PRECISION S
      II = 0
      DO I = 1,N
         T = 1/CHOINV(II+I)
         JJ = 0
         DO J = 1,I-1
            S = 0
            JJ = JJ + J
            KK = JJ
            DO K = J,I-1
               S = S + CHOINV(II+K)*CHOINV(KK)
               KK = KK + K
            END DO
            CHOINV(II+J) = -S*T
         END DO
         II = II + I
         CHOINV(II) = T
      END DO
      END SUBROUTINE CHOLNV
      SUBROUTINE CHOLPI(N, CHOPDI)
*
*     Multiplies Choleski inverse factors in situ
*
      INTEGER I, II, J, JJ, K, KK, N
      DOUBLE PRECISION CHOPDI(*)
      DOUBLE PRECISION S
      II = 0
      DO I = 1,N
         DO J = 1,I
            S = 0
            JJ = II + I
            KK = II + J
            DO K = I,N
               S = S + CHOPDI(KK)*CHOPDI(JJ)
               JJ = JJ + K
               KK = KK + K
            END DO
            CHOPDI(II+J) = S
         END DO
         II = II + I
      END DO
      RETURN
      END SUBROUTINE CHOLPI
      SUBROUTINE RCSWAP(P, Q, A, B, INFIN, N, C)
*
*     Swaps rows and columns P and Q in situ.
*
      DOUBLE PRECISION A(*), B(*), C(*), T
      INTEGER INFIN(*), P, Q, N, I, J, II, JJ
      T = A(P)
      A(P) = A(Q)
      A(Q) = T
      T = B(P)
      B(P) = B(Q)
      B(Q) = T
      J = INFIN(P)
      INFIN(P) = INFIN(Q)
      INFIN(Q) = J
      JJ = (P*(P-1))/2
      II = (Q*(Q-1))/2
      T = C(JJ+P)
      C(JJ+P) = C(II+Q)
      C(II+Q) = T
      DO J = 1, P-1
         T = C(JJ+J)
         C(JJ+J) = C(II+J)
         C(II+J) = T
      END DO
      JJ = JJ + P
      DO I = P+1, Q-1
         T = C(JJ+P)
         C(JJ+P) = C(II+I)
         C(II+I) = T
         JJ = JJ + I
      END DO
      II = II + Q
      DO I = Q+1, N
         T = C(II+P)
         C(II+P) = C(II+Q)
         C(II+Q) = T
         II = II + I
      END DO
      RETURN
      END SUBROUTINE RCSWAP
      DOUBLE PRECISION FUNCTION PHI(Z)
*     
*     Normal distribution probabilities accurate to 1.e-15.
*     Z = no. of standard deviations from the mean.
*     
*     Based upon algorithm 5666 for the error function, from:
*     Hart, J.F. et al, 'Computer Approximations', Wiley 1968
*     
*     Programmer: Alan Miller
*     
*     Latest revision - 30 March 1986
*     
      DOUBLE PRECISION P0, P1, P2, P3, P4, P5, P6, 
     &     Q0, Q1, Q2, Q3, Q4, Q5, Q6, Q7,
     &     Z, P, EXPNTL, CUTOFF, ROOTPI, ZABS
      PARAMETER(
     &     P0 = 220.20 68679 12376 1D0,
     &     P1 = 221.21 35961 69931 1D0, 
     &     P2 = 112.07 92914 97870 9D0,
     &     P3 = 33.912 86607 83830 0D0,
     &     P4 = 6.3739 62203 53165 0D0,
     &     P5 = .70038 30644 43688 1D0, 
     &     P6 = .035262 49659 98910 9D0)
      PARAMETER(
     &     Q0 = 440.41 37358 24752 2D0,
     &     Q1 = 793.82 65125 19948 4D0, 
     &     Q2 = 637.33 36333 78831 1D0,
     &     Q3 = 296.56 42487 79673 7D0, 
     &     Q4 = 86.780 73220 29460 8D0,
     &     Q5 = 16.064 17757 92069 5D0, 
     &     Q6 = 1.7556 67163 18264 2D0,
     &     Q7 = .088388 34764 83184 4D0)
      PARAMETER(ROOTPI = 2.5066 28274 63100 1D0)
      PARAMETER(CUTOFF = 7.0710 67811 86547 5D0)
*     
      ZABS = ABS(Z)
*     
*     |Z| > 37
*     
      IF (ZABS .GT. 37) THEN
         P = 0.d0
      ELSE
*     
*     |Z| <= 37
*     
         EXPNTL = EXP(-ZABS**2/2.d0)
*     
*     |Z| < CUTOFF = 10/SQRT(2)
*     
         IF (ZABS .LT. CUTOFF) THEN
            P = EXPNTL*((((((P6*ZABS + P5)*ZABS + P4)*ZABS + P3)*ZABS
     &           + P2)*ZABS + P1)*ZABS + P0)/(((((((Q7*ZABS + Q6)*ZABS
     &           + Q5)*ZABS + Q4)*ZABS + Q3)*ZABS + Q2)*ZABS + Q1)*ZABS
     &           + Q0)
*     
*     |Z| >= CUTOFF.
*     
         ELSE
            P = EXPNTL/(ZABS + 1.d0/(ZABS + 2.d0/(ZABS + 3.d0/(ZABS + 
     &           4.d0/(ZABS + 0.65D0)))))/ROOTPI
         END IF
      END IF
      IF (Z .GT. 0) P = 1 - P
      PHI = P
      RETURN
      END FUNCTION PHI
      DOUBLE PRECISION FUNCTION PHINV(P)
*
*	ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3
*
*	Produces the normal deviate Z corresponding to a given lower
*	tail area of P.
*
*	The hash sums below are the sums of the mantissas of the
*	coefficients.   They are included for use in checking
*	transcription.
*
      DOUBLE PRECISION SPLIT1, SPLIT2, CONST1, CONST2, 
     &     A0, A1, A2, A3, A4, A5, A6, A7, B1, B2, B3, B4, B5, B6, B7, 
     &     C0, C1, C2, C3, C4, C5, C6, C7, D1, D2, D3, D4, D5, D6, D7, 
     &     E0, E1, E2, E3, E4, E5, E6, E7, F1, F2, F3, F4, F5, F6, F7, 
     &     P, Q, R
      PARAMETER (SPLIT1 = 0.425, SPLIT2 = 5,
     &     CONST1 = 0.180625D0, CONST2 = 1.6D0)
*     
*     Coefficients for P close to 0.5
*     
      PARAMETER (
     &     A0 = 3.38713 28727 96366 6080D0,
     &     A1 = 1.33141 66789 17843 7745D+2,
     &     A2 = 1.97159 09503 06551 4427D+3,
     &     A3 = 1.37316 93765 50946 1125D+4,
     &     A4 = 4.59219 53931 54987 1457D+4,
     &     A5 = 6.72657 70927 00870 0853D+4,
     &     A6 = 3.34305 75583 58812 8105D+4,
     &     A7 = 2.50908 09287 30122 6727D+3,
     &     B1 = 4.23133 30701 60091 1252D+1,
     &     B2 = 6.87187 00749 20579 0830D+2,
     &     B3 = 5.39419 60214 24751 1077D+3,
     &     B4 = 2.12137 94301 58659 5867D+4,
     &     B5 = 3.93078 95800 09271 0610D+4,
     &     B6 = 2.87290 85735 72194 2674D+4,
     &     B7 = 5.22649 52788 52854 5610D+3)
*     HASH SUM AB    55.88319 28806 14901 4439
*     
*     Coefficients for P not close to 0, 0.5 or 1.
*     
      PARAMETER (
     &     C0 = 1.42343 71107 49683 57734D0,
     &     C1 = 4.63033 78461 56545 29590D0,
     &     C2 = 5.76949 72214 60691 40550D0,
     &     C3 = 3.64784 83247 63204 60504D0,
     &     C4 = 1.27045 82524 52368 38258D0,
     &     C5 = 2.41780 72517 74506 11770D-1,
     &     C6 = 2.27238 44989 26918 45833D-2,
     &     C7 = 7.74545 01427 83414 07640D-4,
     &     D1 = 2.05319 16266 37758 82187D0,
     &     D2 = 1.67638 48301 83803 84940D0,
     &     D3 = 6.89767 33498 51000 04550D-1,
     &     D4 = 1.48103 97642 74800 74590D-1,
     &     D5 = 1.51986 66563 61645 71966D-2,
     &     D6 = 5.47593 80849 95344 94600D-4,
     &     D7 = 1.05075 00716 44416 84324D-9)
*     HASH SUM CD    49.33206 50330 16102 89036
*
*	Coefficients for P near 0 or 1.
*
      PARAMETER (
     &     E0 = 6.65790 46435 01103 77720D0,
     &     E1 = 5.46378 49111 64114 36990D0,
     &     E2 = 1.78482 65399 17291 33580D0,
     &     E3 = 2.96560 57182 85048 91230D-1,
     &     E4 = 2.65321 89526 57612 30930D-2,
     &     E5 = 1.24266 09473 88078 43860D-3,
     &     E6 = 2.71155 55687 43487 57815D-5,
     &     E7 = 2.01033 43992 92288 13265D-7,
     &     F1 = 5.99832 20655 58879 37690D-1,
     &     F2 = 1.36929 88092 27358 05310D-1,
     &     F3 = 1.48753 61290 85061 48525D-2,
     &     F4 = 7.86869 13114 56132 59100D-4,
     &     F5 = 1.84631 83175 10054 68180D-5,
     &     F6 = 1.42151 17583 16445 88870D-7,
     &     F7 = 2.04426 31033 89939 78564D-15)
*     HASH SUM EF    47.52583 31754 92896 71629
*     
      Q = ( 2.d0*P - 1.d0 )/2.d0
      IF ( ABS(Q) .LE. SPLIT1 ) THEN
         R = CONST1 - Q*Q
         PHINV = Q*(((((((A7*R + A6)*R + A5)*R + A4)*R + A3)
     &        *R + A2)*R + A1)*R + A0) /
     &        (((((((B7*R + B6)*R + B5)*R + B4)*R + B3)
     &        *R + B2)*R + B1)*R + 1.d0)
      ELSE
         R = MIN( P, 1.d0 - P )
         IF (R .GT. 0) THEN
            R = SQRT( -LOG(R) )
            IF ( R .LE. SPLIT2 ) THEN
               R = R - CONST2
               PHINV = (((((((C7*R + C6)*R + C5)*R + C4)*R + C3)
     &              *R + C2)*R + C1)*R + C0) /
     &              (((((((D7*R + D6)*R + D5)*R + D4)*R + D3)
     &              *R + D2)*R + D1)*R + 1.d0)
            ELSE
               R = R - SPLIT2
               PHINV = (((((((E7*R + E6)*R + E5)*R + E4)*R + E3)
     &              *R + E2)*R + E1)*R + E0) /
     &              (((((((F7*R + F6)*R + F5)*R + F4)*R + F3)
     &              *R + F2)*R + F1)*R + 1.d0)
            END IF
         ELSE
            PHINV = 9
         END IF
         IF ( Q .LT. 0 ) PHINV = - PHINV
      END IF
      END FUNCTION PHINV
      DOUBLE PRECISION FUNCTION BVN ( LOWER, UPPER, INFIN, CORREL )
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
      DOUBLE PRECISION LOWER(*), UPPER(*), CORREL, BVNU
      INTEGER INFIN(*)
      IF ( INFIN(1) .EQ. 2  .AND. INFIN(2) .EQ. 2 ) THEN
         BVN =  BVNU ( LOWER(1), LOWER(2), CORREL )
     +        - BVNU ( UPPER(1), LOWER(2), CORREL )
     +        - BVNU ( LOWER(1), UPPER(2), CORREL )
     +        + BVNU ( UPPER(1), UPPER(2), CORREL )
      ELSE IF ( INFIN(1) .EQ. 2  .AND. INFIN(2) .EQ. 1 ) THEN
         BVN =  BVNU ( LOWER(1), LOWER(2), CORREL )
     +        - BVNU ( UPPER(1), LOWER(2), CORREL )
      ELSE IF ( INFIN(1) .EQ. 1  .AND. INFIN(2) .EQ. 2 ) THEN
         BVN =  BVNU ( LOWER(1), LOWER(2), CORREL )
     +        - BVNU ( LOWER(1), UPPER(2), CORREL )
      ELSE IF ( INFIN(1) .EQ. 2  .AND. INFIN(2) .EQ. 0 ) THEN
         BVN =  BVNU ( -UPPER(1), -UPPER(2), CORREL )
     +        - BVNU ( -LOWER(1), -UPPER(2), CORREL )
      ELSE IF ( INFIN(1) .EQ. 0  .AND. INFIN(2) .EQ. 2 ) THEN
         BVN =  BVNU ( -UPPER(1), -UPPER(2), CORREL )
     +        - BVNU ( -UPPER(1), -LOWER(2), CORREL )
      ELSE IF ( INFIN(1) .EQ. 1  .AND. INFIN(2) .EQ. 0 ) THEN
         BVN =  BVNU ( LOWER(1), -UPPER(2), -CORREL )
      ELSE IF ( INFIN(1) .EQ. 0  .AND. INFIN(2) .EQ. 1 ) THEN
         BVN =  BVNU ( -UPPER(1), LOWER(2), -CORREL )
      ELSE IF ( INFIN(1) .EQ. 1  .AND. INFIN(2) .EQ. 1 ) THEN
         BVN =  BVNU ( LOWER(1), LOWER(2), CORREL )
      ELSE IF ( INFIN(1) .EQ. 0  .AND. INFIN(2) .EQ. 0 ) THEN
         BVN =  BVNU ( -UPPER(1), -UPPER(2), CORREL )
      END IF
      END FUNCTION BVN 
      DOUBLE PRECISION FUNCTION BVNU( SH, SK, R )
*
*     A function for computing bivariate normal probabilities.
*
*       Yihong Ge
*       Department of Computer Science and Electrical Engineering
*       Washington State University
*       Pullman, WA 99164-2752
*       Email : yge@eecs.wsu.edu
*     and
*       Alan Genz
*       Department of Mathematics
*       Washington State University
*       Pullman, WA 99164-3113
*       Email : alangenz@wsu.edu
*
* BVN - calculate the probability that X is larger than SH and Y is
*       larger than SK.
*
* Parameters
*
*   SH  REAL, integration limit
*   SK  REAL, integration limit
*   R   REAL, correlation coefficient
*   LG  INTEGER, number of Gauss Rule Points and Weights
*
      DOUBLE PRECISION BVN, SH, SK, R, ZERO, TWOPI 
      INTEGER I, LG, NG
      PARAMETER ( ZERO = 0, TWOPI = 6.2831 85307 179586 ) 
      DOUBLE PRECISION X(10,3), W(10,3), AS, A, B, C, D, RS, XS
      DOUBLE PRECISION PHI, SN, ASR, H, K, BS, HS, HK
*     Gauss Legendre Points and Weights, N =  6
      DATA ( W(I,1), X(I,1), I = 1,3) /
     &  0.1713244923791705D+00,-0.9324695142031522D+00,
     &  0.3607615730481384D+00,-0.6612093864662647D+00,
     &  0.4679139345726904D+00,-0.2386191860831970D+00/
*     Gauss Legendre Points and Weights, N = 12
      DATA ( W(I,2), X(I,2), I = 1,6) /
     &  0.4717533638651177D-01,-0.9815606342467191D+00,
     &  0.1069393259953183D+00,-0.9041172563704750D+00,
     &  0.1600783285433464D+00,-0.7699026741943050D+00,
     &  0.2031674267230659D+00,-0.5873179542866171D+00,
     &  0.2334925365383547D+00,-0.3678314989981802D+00,
     &  0.2491470458134029D+00,-0.1252334085114692D+00/
*     Gauss Legendre Points and Weights, N = 20
      DATA ( W(I,3), X(I,3), I = 1,10) /
     &  0.1761400713915212D-01,-0.9931285991850949D+00,
     &  0.4060142980038694D-01,-0.9639719272779138D+00,
     &  0.6267204833410906D-01,-0.9122344282513259D+00,
     &  0.8327674157670475D-01,-0.8391169718222188D+00,
     &  0.1019301198172404D+00,-0.7463319064601508D+00,
     &  0.1181945319615184D+00,-0.6360536807265150D+00,
     &  0.1316886384491766D+00,-0.5108670019508271D+00,
     &  0.1420961093183821D+00,-0.3737060887154196D+00,
     &  0.1491729864726037D+00,-0.2277858511416451D+00,
     &  0.1527533871307259D+00,-0.7652652113349733D-01/
      SAVE X, W
      IF ( ABS(R) .LT. 0.3 ) THEN
         NG = 1
         LG = 3
      ELSE IF ( ABS(R) .LT. 0.75 ) THEN
         NG = 2
         LG = 6
      ELSE 
         NG = 3
         LG = 10
      ENDIF
      H = SH
      K = SK 
      HK = H*K
      BVN = 0
      IF ( ABS(R) .LT. 0.925 ) THEN
         HS = ( H*H + K*K )/2
         ASR = ASIN(R)
         DO 10  I = 1, LG
            SN = SIN(ASR*( X(I,NG)+1 )/2)
            BVN = BVN + W(I,NG)*EXP( ( SN*HK - HS )/( 1 - SN*SN ) )
            SN = SIN(ASR*(-X(I,NG)+1 )/2)
            BVN = BVN + W(I,NG)*EXP( ( SN*HK - HS )/( 1 - SN*SN ) )
 10      CONTINUE
         BVN = BVN*ASR/(2*TWOPI) + PHI(-H)*PHI(-K) 
      ELSE
         IF ( R .LT. 0 ) THEN
            K = -K
            HK = -HK
         ENDIF
         IF ( ABS(R) .LT. 1 ) THEN
            AS = ( 1 - R )*( 1 + R )
            A = SQRT(AS)
            BS = ( H - K )**2
            C = ( 4 - HK )/8 
            D = ( 12 - HK )/16
            BVN = A*EXP( -(BS/AS + HK)/2 )
     +             *( 1 - C*(BS - AS)*(1 - D*BS/5)/3 + C*D*AS*AS/5 )
            IF ( HK .GT. -160 ) THEN
               B = SQRT(BS)
               BVN = BVN - EXP(-HK/2)*SQRT(TWOPI)*PHI(-B/A)*B
     +                    *( 1 - C*BS*( 1 - D*BS/5 )/3 ) 
            ENDIF
            A = A/2
            DO 20 I = 1, LG
               XS = ( A*(X(I,NG)+1) )**2
               RS = SQRT( 1 - XS )
               BVN = BVN + A*W(I,NG)*
     +              ( EXP( -BS/(2*XS) - HK/(1+RS) )/RS 
     +              - EXP( -(BS/XS+HK)/2 )*( 1 + C*XS*( 1 + D*XS ) ) )
               XS = AS*(-X(I,NG)+1)**2/4
               RS = SQRT( 1 - XS )
               BVN = BVN + A*W(I,NG)*EXP( -(BS/XS + HK)/2 )
     +                    *( EXP( -HK*(1-RS)/(2*(1+RS)) )/RS 
     +                       - ( 1 + C*XS*( 1 + D*XS ) ) )
 20         CONTINUE 
            BVN = -BVN/TWOPI
         ENDIF
         IF ( R .GT. 0 ) BVN =  BVN + PHI( -MAX( H, K ) )
         IF ( R .LT. 0 ) BVN = -BVN + MAX( ZERO, PHI(-H) - PHI(-K) )
      ENDIF
      BVNU = BVN
      END FUNCTION BVNU
!      END PROGRAM TSTNRM
