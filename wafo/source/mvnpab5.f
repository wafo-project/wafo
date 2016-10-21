      MODULE FUNC
      CONTAINS
*
      DOUBLE PRECISION FUNCTION MVNFNC(N, W)
*     
*     Integrand subroutine
*
      INTEGER :: N
      DOUBLE PRECISION, DIMENSION(:) :: W
      MVNFNC = PRODUCT(W(1:N))
      RETURN
      END  FUNCTION MVNFNC
      END MODULE FUNC

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
      DOUBLE PRECISION :: CONDIT
    
      PRINT '(''        Test of RANMVN, SADMVN, KROMVN and SPHMVN'')'
      PRINT '(12X, ''Requested Accuracy '',F8.5)', MAX(ABSEPS,RELEPS)
      PRINT '(''           Number of Dimensions is '',I2)', N
      PRINT '(''     Maximum # of Function Values is '',I7)', MAXPTS
*
      DO K = 1,3
         CALL KROMVN( N, MAXPTS, ABSEPS, RELEPS, ERRK, VALK, IFTK )
         PRINT *,'VALUE=', VALK
         CALL SADMVN( N, MAXPTS, ABSEPS, RELEPS, ERRS, VALS, IFTS )
        
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
        
      END DO
!      END PROGRAM TSTNRM
      CONTAINS

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
      SUBROUTINE SADMVN( N, MAXPTS,ABSEPS, RELEPS, ERROR, VALUE,INFORM)
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
      USE FUNC
! EXTERNAL MVNFNC
      INTEGER :: N, NL, M, LENWRK, MAXPTS, INFORM, INFIS,
     &     RULCLS, TOTCLS, NEWCLS, MAXCLS
      DOUBLE PRECISION :: ABSEPS, RELEPS, ERROR, VALUE,
     &     OLDVAL, D, E, MVNNIT       !, MVNFNC
      PARAMETER ( NL = 20 )
      PARAMETER ( LENWRK = 20*NL**2 )
      DOUBLE PRECISION, DIMENSION(LENWRK) :: WORK
      IF ( N .GT. 20 .OR. N .LT. 1 ) THEN
         INFORM = 2
         VALUE = 0.d0
         ERROR = 1.d0
         RETURN
      ENDIF
      M = N
      IF ( M .EQ. 0 ) THEN
         VALUE = 1.d0
         ERROR = 0.d0 
      ELSE
*
*        Call the subregion adaptive integration subroutine
*
        
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
      SUBROUTINE KROMVN( N, MAXPTS,ABSEPS,RELEPS,ERROR,VALUE,INFORM)
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
      USE FUNC
!      EXTERNAL MVNFNC
      INTEGER :: N, MAXPTS, INFORM, INFIS, IVLS
      DOUBLE PRECISION :: RELEPS, ABSEPS, ERROR, VALUE, 
     &     E, D, MVNNIT !, MVNFNC
      IF ( N .GT. 100 .OR. N .LT. 1 ) THEN
         INFORM = 2
         VALUE = 0.d0
         ERROR = 1.d0
      ELSE
*
*        Call the lattice rule integration subroutine
*
         IVLS = 0
         CALL KROBOV( N, IVLS, MAXPTS, MVNFNC, 
     &        ABSEPS, RELEPS, ERROR, VALUE, INFORM )
      ENDIF
      RETURN
      END SUBROUTINE KROMVN



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
      END PROGRAM TSTNRM
  
