*
* This file contains a short test program and MVNDST, a subroutine
* for computing multivariate normal distribution function values.
* The file is self contained and should compile without errors on (77) 
* standard Fortran compilers. The test program demonstrates the use of
* MVNDST for computing MVN distribution values for a five dimensional 
* example problem, with three different integration limit combinations.
*
*          Alan Genz
*          Department of Mathematics
*          Washington State University
*          Pullman, WA 99164-3113
*          Email : alangenz@wsu.edu
*
      PROGRAM TSTNRM
*
*     Test program for MVNDST
*
      DOUBLE PRECISION ABSEPS, RELEPS, VAL, ERR
      INTEGER N, NN, I, J, K, IJ, MAXPTS, IFT
      PARAMETER ( N = 5, NN = ( N - 1 )*N/2, MAXPTS = 5000*N*N*N )
      PARAMETER ( ABSEPS = 0.00005, RELEPS = 0 )
      DOUBLE PRECISION CORREL(NN), LOW(N), UP(N)
      INTEGER INFIN(N)
*          Chen Problem
      DATA ( UP(I), I=1,N)  /.0, 1.5198, 1.7817, 1.4755, 1.5949/
      DATA (LOW(I), I=1,N)  /.0,  .0   , 1.7817, 1.4755, 1.5949/
      DATA (INFIN(I), I=1,N)/ 1, 2     , 1     , 1     , 0     /
      DATA (CORREL(I),I=1,NN)/-0.707107,0.0,1.d0,0.0,2*0.5,0.0,3*0.5/
      PRINT '(''               Test of MVNDST'')'
      PRINT '(12X, ''Requested Accuracy '',F8.5)', MAX(ABSEPS,RELEPS)
      PRINT '(''           Number of Dimensions is '',I2)', N
      PRINT '(''     Maximum # of Function Values is '',I7)', MAXPTS
*
      DO K = 1, 3
         PRINT '(/'' I     Limits'')'
         PRINT '(4X,''Lower  Upper  Lower Left of Correlation Matrix'')'
         IJ = 0
         DO I = 1, N
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
         CALL MVNDST( N, LOW, UP, INFIN, CORREL, 
     &                MAXPTS, ABSEPS, RELEPS, ERR, VAL, IFT )
         PRINT '('' Results for:  MVNDST'')'
         PRINT '( ''      Value      :   '', F12.6, I5 )', VAL, IFT
         PRINT '( ''  Error Estimate :   '', ''   ('', F8.6'')'' )', ERR
         INFIN(1) = INFIN(1) - 1
      END DO
      END
*
      SUBROUTINE MVNDST( N, LOWER, UPPER, INFIN, CORREL, MAXPTS,
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
*            THe correlation matrix must be positive semidefinite.
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
      EXTERNAL MVNDFN
      INTEGER N, INFIN(*), MAXPTS, INFORM, INFIS, IVLS
      DOUBLE PRECISION CORREL(*), LOWER(*), UPPER(*), RELEPS, ABSEPS,
     &                 ERROR, VALUE, E, D, MVNDNT, MVNDFN
      COMMON /DKBLCK/IVLS
      IF ( N .GT. 100 .OR. N .LT. 1 ) THEN
         INFORM = 2
         VALUE = 0
         ERROR = 1
      ELSE
         INFORM = MVNDNT(N, CORREL, LOWER, UPPER, INFIN, INFIS, D, E)
         IF ( N-INFIS .EQ. 0 ) THEN
            VALUE = 1
            ERROR = 0
         ELSE IF ( N-INFIS .EQ. 1 ) THEN
            VALUE = E - D
            ERROR = 2E-16
         ELSE
*
*        Call the lattice rule integration subroutine
*
            IVLS = 0
            CALL DKBVRC( N-INFIS-1, IVLS, MAXPTS, MVNDFN, 
     &                   ABSEPS, RELEPS, ERROR, VALUE, INFORM )
         ENDIF
      ENDIF
      END
      DOUBLE PRECISION FUNCTION MVNDFN( N, W )
*     
*     Integrand subroutine
*
      INTEGER N, INFIN(*), INFIS, NL
      DOUBLE PRECISION W(*), LOWER(*), UPPER(*), CORREL(*), D, E
      PARAMETER ( NL = 100 )
      DOUBLE PRECISION COV(NL*(NL+1)/2), A(NL), B(NL), Y(NL)
      INTEGER INFI(NL), I, J, IJ, IK, INFA, INFB,K
      DOUBLE PRECISION SUM, AI, BI, DI, EI, PHINVS, BVNMVN, MVNDNT
      SAVE A, B, INFI, COV
      MVNDFN = 1
      INFA = 0
      INFB = 0
      IK = 1
      IJ = 0
      DO I = 1, N+1
         SUM = 0
         DO J = 1, I-1
            IJ = IJ + 1
            IF ( J .LT. IK ) SUM = SUM + COV(IJ)*Y(J)
         END DO
         IF ( INFI(I) .NE. 0 ) THEN
            IF ( INFA .EQ. 1 ) THEN
               AI = MAX( AI, A(I) - SUM )
            ELSE
               AI = A(I) - SUM 
               INFA = 1
            END IF
         END IF
         IF ( INFI(I) .NE. 1 ) THEN
            IF ( INFB .EQ. 1 ) THEN
               BI = MIN( BI, B(I) - SUM )
            ELSE
               BI = B(I) - SUM 
               INFB = 1
            END IF
         END IF
         IJ = IJ + 1
         IF ( I .EQ. N+1 .OR. COV(IJ+IK+1) .GT. 0 ) THEN 
            CALL MVNLMS( AI, BI, 2*INFA+INFB-1, DI, EI )
            IF ( DI .GE. EI ) THEN
               MVNDFN = 0
               RETURN
            ELSE
               MVNDFN = MVNDFN*( EI - DI )
               IF ( I .LE. N ) Y(IK) = PHINVS( DI + W(IK)*( EI - DI ) )
               IK = IK + 1
               INFA = 0
               INFB = 0
            END IF
         END IF
      END DO
      RETURN
*
*     Entry point for intialization.
*
      ENTRY MVNDNT( N, CORREL, LOWER, UPPER, INFIN, INFIS, D, E )
      MVNDNT = 0
*
*     Initialization and computation of covariance Cholesky factor.
*
      CALL COVSRT( N, LOWER,UPPER,CORREL,INFIN,Y, INFIS,A,B,COV,INFI )
      IF ( N - INFIS .EQ. 1 ) THEN
         CALL MVNLMS( A(1), B(1), INFI(1), D, E ) 
      ELSE IF ( N - INFIS .EQ. 2 ) THEN
         IF ( ABS( COV(3) ) .GT. 0 ) THEN
            D = SQRT( 1 + COV(2)**2 )
            IF ( INFI(2) .NE. 0 ) A(2) = A(2)/D
            IF ( INFI(2) .NE. 1 ) B(2) = B(2)/D
            E = BVNMVN( A, B, INFI, COV(2)/D )
            D = 0
         ELSE
            IF ( INFI(1) .NE. 0 ) THEN
               IF ( INFI(2) .NE. 0 ) A(1) = MAX( A(1), A(2) )
            ELSE
               IF ( INFI(2) .NE. 0 ) A(1) = A(2) 
            END IF
            IF ( INFI(1) .NE. 1 ) THEN
               IF ( INFI(2) .NE. 1 ) B(1) = MIN( B(1), B(2) ) 
            ELSE
               IF ( INFI(2) .NE. 1 ) B(1) = B(2)
            END IF
            IF ( INFI(1) .NE. INFI(2) ) INFI(1) = 2
            CALL MVNLMS( A(1), B(1), INFI(1), D, E ) 
         END IF
         INFIS = INFIS + 1 
      END IF
     
      PRINT '(/'' I     Limits'')'
      PRINT '(4X,''Lower  Upper  Lower Left of Cholesky Matrix'')'
      IJ = 0
      DO I = 1, N
         IF ( INFI(I) .LT. 0 ) THEN 
            PRINT '(I2, '' -infin  infin '', 7F9.5)',
     &           I, ( COV(IJ+J), J = 1,I )
         ELSE IF ( INFI(I) .EQ. 0 ) THEN 
            PRINT '(I2, '' -infin'', F7.4, 1X, 7F9.5)',
     &           I, B(I), ( COV(IJ+J), J = 1,I )
         ELSE IF ( INFI(I) .EQ. 1 ) THEN 
            PRINT '(I2, F7.4, ''  infin '', 7F9.5)',
     &           I, A(I), ( COV(IJ+J), J = 1,I )
         ELSE 
            PRINT '(I2, 2F7.4, 1X, 7F9.5)', 
     &           I, A(I), B(I), ( COV(IJ+J), J = 1,I )
         END IF
         IJ = IJ + I
      END DO
      END
      SUBROUTINE MVNLMS( A, B, INFIN, LOWER, UPPER )
      DOUBLE PRECISION A, B, LOWER, UPPER, MVNPHI
      INTEGER INFIN
      LOWER = 0
      UPPER = 1
      IF ( INFIN .GE. 0 ) THEN
         IF ( INFIN .NE. 0 ) LOWER = MVNPHI(A)
         IF ( INFIN .NE. 1 ) UPPER = MVNPHI(B)
      ENDIF
      END      
      SUBROUTINE COVSRT( N, LOWER, UPPER, CORREL, INFIN, Y, 
     &               INFIS,     A,     B,    COV,  INFI )
*
*     Subroutine to sort integration limits and determine Cholesky factor.
*
      INTEGER N, INFI(*), INFIN(*), INFIS
      DOUBLE PRECISION CDI(5),
     &     A(*), B(*), COV(*), LOWER(*), UPPER(*), CORREL(*), Y(*)
      INTEGER I, J, K, L, M, II, IJ, IL, JMIN,TMP
      integer index1(5)
      DOUBLE PRECISION SUMSQ, AJ, BJ, SUM, SQTWPI, EPS, D, E
      DOUBLE PRECISION CVDIAG, AMIN, BMIN, DMIN, EMIN, YL, YU
      PARAMETER ( SQTWPI = 2.50662 82746 31001D0, EPS = 1D-10 )
      DO I=1,5
       index1(I)=I
      END DO
      IJ = 0
      II = 0
      INFIS = 0
      DO I = 1, N
         A(I) = 0
         B(I) = 0
         INFI(I) = INFIN(I) 
         IF ( INFI(I) .LT. 0 ) THEN
            INFIS = INFIS + 1
         ELSE 
            IF ( INFI(I) .NE. 0 ) A(I) = LOWER(I)
            IF ( INFI(I) .NE. 1 ) B(I) = UPPER(I)
         ENDIF
         DO J = 1, I-1
            IJ = IJ + 1
            II = II + 1
            COV(IJ) = CORREL(II)
         END DO
         IJ = IJ + 1
         COV(IJ) = 1
      END DO
*
*     First move any doubly infinite limits to innermost positions.
*
      IF ( INFIS .LT. N ) THEN
         DO I = N, N-INFIS+1, -1
            IF ( INFI(I) .GE. 0 ) THEN 
               DO J = 1,I-1
                  IF ( INFI(J) .LT. 0 ) THEN
                     TMP=index1(j)
                     index1(j)=index1(i)
                     index1(i)=tmp
                     CALL RCSWP( J, I, A, B, INFI, N, COV )
                     GO TO 10
                  ENDIF
               END DO
            ENDIF
 10      END DO
*
*     Sort remaining limits and determine Cholesky factor.
*
         II = 0
         DO I = 1, N-INFIS
*
*        Determine the integration limits for variable with minimum
*        expected probability and interchange that variable with Ith.
*
            DMIN = 0
            EMIN = 1
            JMIN = I
            CVDIAG = 0
            IJ = II
            DO J = I, N-INFIS
               IF ( COV(IJ+J) .GT. EPS ) THEN
                  SUMSQ = SQRT( COV(IJ+J) )
                  SUM = 0
                  DO K = 1, I-1
                     SUM = SUM + COV(IJ+K)*Y(K)
                  END DO
                  AJ = ( A(J) - SUM )/SUMSQ
                  BJ = ( B(J) - SUM )/SUMSQ
                  CALL MVNLMS( AJ, BJ, INFI(J), D, E )
                  IF ( EMIN + D .GE. E + DMIN ) THEN
                     JMIN = J
                     AMIN = AJ
                     BMIN = BJ
                     DMIN = D
                     EMIN = E
                     CVDIAG = SUMSQ
                  ENDIF
               ENDIF
               IJ = IJ + J 
            END DO
            IF ( JMIN .GT. I ) THEN 
               TMP=index1(jmin)
               index1(jmin)=index1(i)
               index1(i)=tmp
               CALL RCSWP( I, JMIN, A,B, INFI, N, COV )
            END IF
            COV(II+I) = CVDIAG
*
*        Compute Ith column of Cholesky factor.
*        Compute expected value for Ith integration variable and
*         scale Ith covariance matrix row and limits.
*
            IF ( CVDIAG .GT. 0 ) THEN
               CDI(I)=CVDIAG
               IL = II + I
               DO L = I+1, N-INFIS               
                  COV(IL+I) = COV(IL+I)/CVDIAG
                  IJ = II + I
                  DO J = I+1, L
                     COV(IL+J) = COV(IL+J) - COV(IL+I)*COV(IJ+I)
                     IJ = IJ + J
                  END DO
                  IL = IL + L
               END DO
               IF ( EMIN .GT. DMIN + EPS ) THEN
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
               DO J = 1, I
                  II = II + 1
                  COV(II) = COV(II)/CVDIAG
               END DO
               A(I) = A(I)/CVDIAG
               B(I) = B(I)/CVDIAG
            ELSE
               CDI(I)=0.d0
               IL = II + I
               DO L = I+1, N-INFIS                
                  COV(IL+I) = 0
                  IL = IL + L
               END DO
*
*        If the covariance matrix diagonal entry is zero, 
*         permute limits and/or rows, if necessary.
*
*
               DO J = I-1, 1, -1
                  IF ( ABS( COV(II+J) ) .GT. EPS ) THEN
                     A(I) = A(I)/COV(II+J)
                     B(I) = B(I)/COV(II+J)
                     IF ( COV(II+J) .LT. 0 ) THEN
                        CALL DKSWAP( A(I), B(I) ) 
                        IF ( INFI(I) .NE. 2 ) INFI(I) = 1 - INFI(I)
                     END IF
                     DO L = 1, J
                        COV(II+L) = COV(II+L)/COV(II+J)
                     END DO
                     DO L = J+1, I-1 
                        IF( COV((L-1)*L/2+J+1) .GT. 0 ) THEN
                           IJ = II
                           DO K = I-1, L, -1 
                              DO M = 1, K
                                 CALL DKSWAP( COV(IJ-K+M), COV(IJ+M) )
                              END DO
                              CALL DKSWAP( A(K), A(K+1) ) 
                              CALL DKSWAP( B(K), B(K+1) ) 
                              M = INFI(K)
                              INFI(K) = INFI(K+1)
                              INFI(K+1) = M
                              M=INDEX1(K)
                              INDEX1(K)=INDEX1(K+1)
                              INDEX1(K+1)=M
                              IJ = IJ - K 
                           END DO
                           GO TO 20
                        END IF
                     END DO
                     GO TO 20
                  END IF
                  COV(II+J) = 0
               END DO
 20            II = II + I
               Y(I) = 0
            END IF
         END DO
      ENDIF
      PRINT *,index1
      PRINT *,CDI
!      DO I=1,N*(N+1)/2
!         PRINT *,COV(I)
!      END DO
      END
*
      SUBROUTINE DKSWAP( X, Y )
      DOUBLE PRECISION X, Y, T
      T = X
      X = Y
      Y = T
      END
*
      SUBROUTINE RCSWP( P, Q, A, B, INFIN, N, C )
*
*     Swaps rows and columns P and Q in situ, with P <= Q.
*
      DOUBLE PRECISION A(*), B(*), C(*)
      INTEGER INFIN(*), P, Q, N, I, J, II, JJ
      CALL DKSWAP( A(P), A(Q) )
      CALL DKSWAP( B(P), B(Q) )
      J = INFIN(P)
      INFIN(P) = INFIN(Q)
      INFIN(Q) = J
      JJ = ( P*( P - 1 ) )/2
      II = ( Q*( Q - 1 ) )/2
      CALL DKSWAP( C(JJ+P), C(II+Q) )
      DO J = 1, P-1
         CALL DKSWAP( C(JJ+J), C(II+J) )
      END DO
      JJ = JJ + P
      DO I = P+1, Q-1
         CALL DKSWAP( C(JJ+P), C(II+I) )
         JJ = JJ + I
      END DO
      II = II + Q
      DO I = Q+1, N
         CALL DKSWAP( C(II+P), C(II+Q) )
         II = II + I
      END DO
      END
*
      SUBROUTINE DKBVRC( NDIM, MINVLS, MAXVLS, FUNCTN, ABSEPS, RELEPS,
     &                   ABSERR, FINEST, INFORM )
*
*  Automatic Multidimensional Integration Subroutine
*               
*         AUTHOR: Alan Genz
*                 Department of Mathematics
*                 Washington State University
*                 Pulman, WA 99164-3113
*                 Email: AlanGenz@wsu.edu
*
*         Last Change: 5/15/98
*
*  KRBVRC computes an approximation to the integral
*
*      1  1     1
*     I  I ... I       F(X)  dx(NDIM)...dx(2)dx(1)
*      0  0     0
*
*
*  DKBVRC uses randomized Korobov rules for the first 20 variables. 
*  The primary references are
*   "Randomization of Number Theoretic Methods for Multiple Integration"
*    R. Cranley and T.N.L. Patterson, SIAM J Numer Anal, 13, pp. 904-14,
*  and 
*   "Optimal Parameters for Multidimensional Integration", 
*    P. Keast, SIAM J Numer Anal, 10, pp.831-838.
*  If there are more than 20 variables, the remaining variables are
*  integrated using Richtmeyer rules. A reference is
*   "Methods of Numerical Integration", P.J. Davis and P. Rabinowitz, 
*    Academic Press, 1984, pp. 482-483.
*   
***************  Parameters ********************************************
****** Input parameters
*  NDIM    Number of variables, must exceed 1, but not exceed 40
*  MINVLS  Integer minimum number of function evaluations allowed.
*          MINVLS must not exceed MAXVLS.  If MINVLS < 0 then the
*          routine assumes a previous call has been made with 
*          the same integrand and continues that calculation.
*  MAXVLS  Integer maximum number of function evaluations allowed.
*  FUNCTN  EXTERNALly declared user defined function to be integrated.
*          It must have parameters (NDIM,Z), where Z is a real array
*          of dimension NDIM.
*                                     
*  ABSEPS  Required absolute accuracy.
*  RELEPS  Required relative accuracy.
****** Output parameters
*  MINVLS  Actual number of function evaluations used.
*  ABSERR  Estimated absolute accuracy of FINEST.
*  FINEST  Estimated value of integral.
*  INFORM  INFORM = 0 for normal exit, when 
*                     ABSERR <= MAX(ABSEPS, RELEPS*ABS(FINEST))
*                  and 
*                     INTVLS <= MAXCLS.
*          INFORM = 1 If MAXVLS was too small to obtain the required 
*          accuracy. In this case a value FINEST is returned with 
*          estimated absolute accuracy ABSERR.
************************************************************************
      EXTERNAL FUNCTN
      INTEGER NDIM, MINVLS, MAXVLS, INFORM, NP, PLIM, NLIM, KLIM, KLIMI,
     &        SAMPLS, I, INTVLS, MINSMP
      PARAMETER ( PLIM = 25, NLIM = 100, KLIM = 20, MINSMP = 8 )
      INTEGER P(PLIM), C(PLIM,KLIM-1) 
      DOUBLE PRECISION FUNCTN, ABSEPS, RELEPS, FINEST, ABSERR, DIFINT, 
     &                 FINVAL, VARSQR, VAREST, VARPRD, VALUE
      DOUBLE PRECISION X(2*NLIM), VK(KLIM), ONE
      PARAMETER ( ONE = 1 )
      SAVE P, C, SAMPLS, NP, VAREST
      INFORM = 1
      INTVLS = 0
      KLIMI = KLIM
      IF ( MINVLS .GE. 0 ) THEN
         FINEST = 0
         VAREST = 0
         SAMPLS = MINSMP 
         DO I = 1, PLIM
            NP = I
            IF ( MINVLS .LT. 2*SAMPLS*P(I) ) GO TO 10
         END DO
         SAMPLS = MAX( MINSMP, MINVLS/( 2*P(NP) ) )
      ENDIF
 10   VK(1) = ONE/P(NP)
      DO I = 2, MIN( NDIM, KLIM )
         VK(I) = MOD( C(NP, MIN(NDIM-1,KLIM-1))*VK(I-1), ONE )
      END DO
      FINVAL = 0
      VARSQR = 0
      DO I = 1, SAMPLS
         CALL DKSMRC( NDIM, KLIMI, VALUE, P(NP), VK, FUNCTN, X )
         DIFINT = ( VALUE - FINVAL )/I
         FINVAL = FINVAL + DIFINT
         VARSQR = ( I - 2 )*VARSQR/I + DIFINT**2
      END DO
      INTVLS = INTVLS + 2*SAMPLS*P(NP)
      VARPRD = VAREST*VARSQR
      FINEST = FINEST + ( FINVAL - FINEST )/( 1 + VARPRD )
      IF ( VARSQR .GT. 0 ) VAREST = ( 1 + VARPRD )/VARSQR
      ABSERR = 3*SQRT( VARSQR/( 1 + VARPRD ) )
      IF ( ABSERR .GT. MAX( ABSEPS, ABS(FINEST)*RELEPS ) ) THEN
         IF ( NP .LT. PLIM ) THEN
            NP = NP + 1
         ELSE
            SAMPLS = MIN( 3*SAMPLS/2, ( MAXVLS - INTVLS )/( 2*P(NP) ) ) 
            SAMPLS = MAX( MINSMP, SAMPLS )
         ENDIF
         IF ( INTVLS + 2*SAMPLS*P(NP) .LE. MAXVLS ) GO TO 10
      ELSE
         INFORM = 0
      ENDIF
      MINVLS = INTVLS
*
      DATA P( 1), (C( 1,I), I = 1, 19)/     31,     12,      9,      9,
     &      13,     12,     12,     12,     12,     12,     12,     12,
     &      12,      3,      3,      3,     12,      7,      7,     12/
      DATA P( 2), (C( 2,I), I = 1, 19)/     47,     13,     11,     17,
     &      10,     15,     15,     15,     15,     15,     15,     22,
     &      15,     15,      6,      6,      6,     15,     15,      9/
      DATA P( 3), (C( 3,I), I = 1, 19)/     73,     27,     28,     10,
     &      11,     11,     20,     11,     11,     28,     13,     13,
     &      28,     13,     13,     13,     14,     14,     14,     14/
      DATA P( 4), (C( 4,I), I = 1, 19)/    113,     35,     27,     27,
     &      36,     22,     29,     29,     20,     45,      5,      5,
     &       5,     21,     21,     21,     21,     21,     21,     21/
      DATA P( 5), (C( 5,I), I = 1, 19)/    173,     64,     66,     28,
     &      28,     44,     44,     55,     67,     10,     10,     10,
     &      10,     10,     10,     38,     38,     10,     10,     10/
      DATA P( 6), (C( 6,I), I = 1, 19)/    263,    111,     42,     54,
     &     118,     20,     31,     31,     72,     17,     94,     14,
     &      14,     11,     14,     14,     14,     94,     10,     10/
      DATA P( 7), (C( 7,I), I = 1, 19)/    397,    163,    154,     83,
     &      43,     82,     92,    150,     59,     76,     76,     47,
     &      11,     11,    100,    131,    116,    116,    116,    116/
      DATA P( 8), (C( 8,I), I = 1, 19)/    593,    246,    189,    242,
     &     102,    250,    250,    102,    250,    280,    118,    196,
     &     118,    191,    215,    121,    121,     49,     49,     49/
      DATA P( 9), (C( 9,I), I = 1, 19)/    907,    347,    402,    322,
     &     418,    215,    220,    339,    339,    339,    337,    218,
     &     315,    315,    315,    315,    167,    167,    167,    167/
      DATA P(10), (C(10,I), I = 1, 19)/   1361,    505,    220,    601,
     &     644,    612,    160,    206,    206,    206,    422,    134,
     &     518,    134,    134,    518,    652,    382,    206,    158/
      DATA P(11), (C(11,I), I = 1, 19)/   2053,    794,    325,    960,
     &     528,    247,    247,    338,    366,    847,    753,    753,
     &     236,    334,    334,    461,    711,    652,    381,    381/
      DATA P(12), (C(12,I), I = 1, 19)/   3079,   1189,    888,    259,
     &    1082,    725,    811,    636,    965,    497,    497,   1490,
     &    1490,    392,   1291,    508,    508,   1291,   1291,    508/
      DATA P(13), (C(13,I), I = 1, 19)/   4621,   1763,   1018,   1500,
     &     432,   1332,   2203,    126,   2240,   1719,   1284,    878,
     &    1983,    266,    266,    266,    266,    747,    747,    127/
      DATA P(14), (C(14,I), I = 1, 19)/   6947,   2872,   3233,   1534,
     &    2941,   2910,    393,   1796,    919,    446,    919,    919,
     &    1117,    103,    103,    103,    103,    103,    103,    103/
      DATA P(15), (C(15,I), I = 1, 19)/  10427,   4309,   3758,   4034,
     &    1963,    730,    642,   1502,   2246,   3834,   1511,   1102,
     &    1102,   1522,   1522,   3427,   3427,   3928,    915,    915/
      DATA P(16), (C(16,I), I = 1, 19)/  15641,   6610,   6977,   1686,
     &    3819,   2314,   5647,   3953,   3614,   5115,    423,    423,
     &    5408,   7426,    423,    423,    487,   6227,   2660,   6227/
      DATA P(17), (C(17,I), I = 1, 19)/  23473,   9861,   3647,   4073,
     &    2535,   3430,   9865,   2830,   9328,   4320,   5913,  10365,
     &    8272,   3706,   6186,   7806,   7806,   7806,   8610,   2563/
      DATA P(18), (C(18,I), I = 1, 19)/  35221,  10327,   7582,   7124,
     &    8214,   9600,  10271,  10193,  10800,   9086,   2365,   4409,
     &   13812,   5661,   9344,   9344,  10362,   9344,   9344,   8585/
      DATA P(19), (C(19,I), I = 1, 19)/  52837,  19540,  19926,  11582,
     &   11113,  24585,   8726,  17218,    419,   4918,   4918,   4918,
     &   15701,  17710,   4037,   4037,  15808,  11401,  19398,  25950/
      DATA P(20), (C(20,I), I = 1, 19)/  79259,  34566,   9579,  12654,
     &   26856,  37873,  38806,  29501,  17271,   3663,  10763,  18955,
     &    1298,  26560,  17132,  17132,   4753,   4753,   8713,  18624/
      DATA P(21), (C(21,I), I = 1, 19)/ 118891,  31929,  49367,  10982,
     &    3527,  27066,  13226,  56010,  18911,  40574,  20767,  20767,
     &    9686,  47603,  47603,  11736,  11736,  41601,  12888,  32948/
      DATA P(22), (C(22,I), I = 1, 19)/ 178349,  40701,  69087,  77576,
     &   64590,  39397,  33179,  10858,  38935,  43129,  35468,  35468,
     &    2196,  61518,  61518,  27945,  70975,  70975,  86478,  86478/
      DATA P(23), (C(23,I), I = 1, 19)/ 267523, 103650, 125480,  59978,
     &   46875,  77172,  83021, 126904,  14541,  56299,  43636,  11655,
     &   52680,  88549,  29804, 101894, 113675,  48040, 113675,  34987/
      DATA P(24), (C(24,I), I = 1, 19)/ 401287, 165843,  90647,  59925,
     &  189541,  67647,  74795,  68365, 167485, 143918,  74912, 167289,
     &   75517,   8148, 172106, 126159,  35867,  35867,  35867, 121694/
      DATA P(25), (C(25,I), I = 1, 19)/ 601942, 130365, 236711, 110235,
     &  125699,  56483,  93735, 234469,  60549,   1291,  93937, 245291,
     &  196061, 258647, 162489, 176631, 204895,  73353, 172319,  28881/
*
      END
*
      SUBROUTINE DKSMRC( NDIM, KLIM, SUMKRO, PRIME, VK, FUNCTN, X )
      EXTERNAL FUNCTN
      INTEGER NDIM, KLIM, PRIME, K, J, JP, NK
      DOUBLE PRECISION SUMKRO, VK(*), FUNCTN, X(*), ONE, XT, MVNUNI
      PARAMETER ( ONE = 1 )
      SUMKRO = 0
      NK = MIN( NDIM, KLIM )
      DO J = 1, NK-1
         JP = J + MVNUNI()*( NK + 1 - J ) 
         XT = VK(J)
         VK(J) = VK(JP)
         VK(JP) = XT
      END DO
      DO J = 1, NDIM
         X(NDIM+J) = MVNUNI()
      END DO
      DO K = 1, PRIME
         DO J = 1, NK
            X(J) = MOD( K*VK(J), ONE )
         END DO
         IF ( NDIM. GT. KLIM ) CALL DKRCHT( NDIM-KLIM, X(KLIM+1) )
         DO J = 1, NDIM
            XT = X(J) + X(NDIM+J)
            IF ( XT .GT. 1 ) XT = XT - 1
            X(J) = ABS( 2*XT - 1 )
         END DO
         SUMKRO = SUMKRO + ( FUNCTN(NDIM,X) - SUMKRO )/( 2*K - 1 )
         DO J = 1, NDIM
            X(J) = 1 - X(J)
         END DO
         SUMKRO = SUMKRO + ( FUNCTN(NDIM,X) - SUMKRO )/( 2*K )
      END DO
      END
*
      SUBROUTINE DKRCHT( S, QUASI )
*
*     This subroutine generates a new quasi-random Richtmeyer vector. 
*     A reference is
*      "Methods of Numerical Integration", P.J. Davis and P. Rabinowitz, 
*       Academic Press, 1984, pp. 482-483.
*
*
*       INPUTS:
*         S - the number of dimensions; 
*             KRRCHT is initialized for each new S or S < 1.
*
*       OUTPUTS:
*         QUASI - a new quasi-random S-vector
*
      INTEGER MXDIM, MXHSUM, B
      PARAMETER ( MXDIM = 80, MXHSUM = 48, B = 2 )
      INTEGER S, HISUM, I, PRIME(MXDIM), OLDS, N(0:MXHSUM)
      DOUBLE PRECISION QUASI(*), RN, PSQT(MXDIM), ONE
      PARAMETER ( ONE = 1 )
*
      SAVE OLDS, PRIME, PSQT, HISUM, N
      DATA OLDS / 0 /
      IF ( S .NE. OLDS .OR. S .LT. 1 ) THEN
         OLDS = S
         N(0) = 0
         HISUM = 0
         DO I = 1, S
            RN = PRIME(I)
            PSQT(I) = SQRT( RN )
         END DO
      END IF
      DO I = 0, HISUM 
         N(I) = N(I) + 1
         IF ( N(I) .LT. B ) GO TO 10
         N(I) = 0
      END DO
      HISUM = HISUM + 1
      IF ( HISUM .GT. 48 ) HISUM = 0
      N(HISUM) = 1
 10   RN = 0
      DO I = HISUM, 0, -1
         RN = N(I) + B*RN
      END DO
      DO I = 1, S
         QUASI(I) = MOD( RN*PSQT(I), ONE )
      END DO
*
      DATA ( PRIME(I), I =    1,  80 ) / 
     &     2,    3,    5,    7,   11,   13,   17,   19,   23,   29,
     &    31,   37,   41,   43,   47,   53,   59,   61,   67,   71,
     &    73,   79,   83,   89,   97,  101,  103,  107,  109,  113,
     &   127,  131,  137,  139,  149,  151,  157,  163,  167,  173,
     &   179,  181,  191,  193,  197,  199,  211,  223,  227,  229,
     &   233,  239,  241,  251,  257,  263,  269,  271,  277,  281,
     &   283,  293,  307,  311,  313,  317,  331,  337,  347,  349,
     &   353,  359,  367,  373,  379,  383,  389,  397,  401,  409/
      END
*
      DOUBLE PRECISION FUNCTION MVNPHI( Z )
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
     *     Q0, Q1, Q2, Q3, Q4, Q5, Q6, Q7,
     *     Z, P, EXPNTL, CUTOFF, ROOTPI, ZABS
      PARAMETER(
     *     P0 = 220.20 68679 12376 1D0,
     *     P1 = 221.21 35961 69931 1D0, 
     *     P2 = 112.07 92914 97870 9D0,
     *     P3 = 33.912 86607 83830 0D0,
     *     P4 = 6.3739 62203 53165 0D0,
     *     P5 = .70038 30644 43688 1D0, 
     *     P6 = .035262 49659 98910 9D0 )
      PARAMETER(
     *     Q0 = 440.41 37358 24752 2D0,
     *     Q1 = 793.82 65125 19948 4D0, 
     *     Q2 = 637.33 36333 78831 1D0,
     *     Q3 = 296.56 42487 79673 7D0, 
     *     Q4 = 86.780 73220 29460 8D0,
     *     Q5 = 16.064 17757 92069 5D0, 
     *     Q6 = 1.7556 67163 18264 2D0,
     *     Q7 = .088388 34764 83184 4D0 )
      PARAMETER( ROOTPI = 2.5066 28274 63100 1D0 )
      PARAMETER( CUTOFF = 7.0710 67811 86547 5D0 )
*     
      ZABS = ABS(Z)
*     
*     |Z| > 37
*     
      IF ( ZABS .GT. 37 ) THEN
         P = 0
      ELSE
*     
*     |Z| <= 37
*     
         EXPNTL = EXP( -ZABS**2/2 )
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
            P = EXPNTL/( ZABS + 1/( ZABS + 2/( ZABS + 3/( ZABS 
     *                        + 4/( ZABS + 0.65D0 ) ) ) ) )/ROOTPI
         END IF
      END IF
      IF ( Z .GT. 0 ) P = 1 - P
      MVNPHI = P
      END
      DOUBLE PRECISION FUNCTION PHINVS(P)
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
     *     A0, A1, A2, A3, A4, A5, A6, A7, B1, B2, B3, B4, B5, B6, B7, 
     *     C0, C1, C2, C3, C4, C5, C6, C7, D1, D2, D3, D4, D5, D6, D7, 
     *     E0, E1, E2, E3, E4, E5, E6, E7, F1, F2, F3, F4, F5, F6, F7, 
     *     P, Q, R
      PARAMETER ( SPLIT1 = 0.425, SPLIT2 = 5,
     *            CONST1 = 0.180625D0, CONST2 = 1.6D0 )
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
      Q = ( 2*P - 1 )/2
      IF ( ABS(Q) .LE. SPLIT1 ) THEN
         R = CONST1 - Q*Q
         PHINVS = Q*( ( ( ((((A7*R + A6)*R + A5)*R + A4)*R + A3)
     *                  *R + A2 )*R + A1 )*R + A0 )
     *            /( ( ( ((((B7*R + B6)*R + B5)*R + B4)*R + B3)
     *                  *R + B2 )*R + B1 )*R + 1 )
      ELSE
         R = MIN( P, 1 - P )
         IF ( R .GT. 0 ) THEN
            R = SQRT( -LOG(R) )
            IF ( R .LE. SPLIT2 ) THEN
               R = R - CONST2
               PHINVS = ( ( ( ((((C7*R + C6)*R + C5)*R + C4)*R + C3)
     *                      *R + C2 )*R + C1 )*R + C0 ) 
     *                /( ( ( ((((D7*R + D6)*R + D5)*R + D4)*R + D3)
     *                      *R + D2 )*R + D1 )*R + 1 )
            ELSE
               R = R - SPLIT2
               PHINVS = ( ( ( ((((E7*R + E6)*R + E5)*R + E4)*R + E3)
     *                      *R + E2 )*R + E1 )*R + E0 )
     *                /( ( ( ((((F7*R + F6)*R + F5)*R + F4)*R + F3)
     *                      *R + F2 )*R + F1 )*R + 1 )
            END IF
         ELSE
            PHINVS = 9
         END IF
         IF ( Q .LT. 0 ) PHINVS = - PHINVS
      END IF
      END
      DOUBLE PRECISION FUNCTION BVNMVN( LOWER, UPPER, INFIN, CORREL )
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
      DOUBLE PRECISION LOWER(*), UPPER(*), CORREL, BVU
      INTEGER INFIN(*)
      IF ( INFIN(1) .EQ. 2  .AND. INFIN(2) .EQ. 2 ) THEN
         BVNMVN =  BVU ( LOWER(1), LOWER(2), CORREL )
     +           - BVU ( UPPER(1), LOWER(2), CORREL )
     +           - BVU ( LOWER(1), UPPER(2), CORREL )
     +           + BVU ( UPPER(1), UPPER(2), CORREL )
      ELSE IF ( INFIN(1) .EQ. 2  .AND. INFIN(2) .EQ. 1 ) THEN
         BVNMVN =  BVU ( LOWER(1), LOWER(2), CORREL )
     +           - BVU ( UPPER(1), LOWER(2), CORREL )
      ELSE IF ( INFIN(1) .EQ. 1  .AND. INFIN(2) .EQ. 2 ) THEN
         BVNMVN =  BVU ( LOWER(1), LOWER(2), CORREL )
     +           - BVU ( LOWER(1), UPPER(2), CORREL )
      ELSE IF ( INFIN(1) .EQ. 2  .AND. INFIN(2) .EQ. 0 ) THEN
         BVNMVN =  BVU ( -UPPER(1), -UPPER(2), CORREL )
     +           - BVU ( -LOWER(1), -UPPER(2), CORREL )
      ELSE IF ( INFIN(1) .EQ. 0  .AND. INFIN(2) .EQ. 2 ) THEN
         BVNMVN =  BVU ( -UPPER(1), -UPPER(2), CORREL )
     +           - BVU ( -UPPER(1), -LOWER(2), CORREL )
      ELSE IF ( INFIN(1) .EQ. 1  .AND. INFIN(2) .EQ. 0 ) THEN
         BVNMVN =  BVU ( LOWER(1), -UPPER(2), -CORREL )
      ELSE IF ( INFIN(1) .EQ. 0  .AND. INFIN(2) .EQ. 1 ) THEN
         BVNMVN =  BVU ( -UPPER(1), LOWER(2), -CORREL )
      ELSE IF ( INFIN(1) .EQ. 1  .AND. INFIN(2) .EQ. 1 ) THEN
         BVNMVN =  BVU ( LOWER(1), LOWER(2), CORREL )
      ELSE IF ( INFIN(1) .EQ. 0  .AND. INFIN(2) .EQ. 0 ) THEN
         BVNMVN =  BVU ( -UPPER(1), -UPPER(2), CORREL )
      END IF
      END 
      DOUBLE PRECISION FUNCTION BVU( SH, SK, R )
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
      DOUBLE PRECISION MVNPHI, SN, ASR, H, K, BS, HS, HK
*     Gauss Legendre Points and Weights, N =  6
      DATA ( W(I,1), X(I,1), I = 1,3) /
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
         DO I = 1, LG
            SN = SIN(ASR*( X(I,NG)+1 )/2)
            BVN = BVN + W(I,NG)*EXP( ( SN*HK - HS )/( 1 - SN*SN ) )
            SN = SIN(ASR*(-X(I,NG)+1 )/2)
            BVN = BVN + W(I,NG)*EXP( ( SN*HK - HS )/( 1 - SN*SN ) )
         END DO
         BVN = BVN*ASR/(2*TWOPI) + MVNPHI(-H)*MVNPHI(-K) 
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
               BVN = BVN - EXP(-HK/2)*SQRT(TWOPI)*MVNPHI(-B/A)*B
     +                    *( 1 - C*BS*( 1 - D*BS/5 )/3 ) 
            ENDIF
            A = A/2
            DO I = 1, LG
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
            END DO
            BVN = -BVN/TWOPI
         ENDIF
         IF ( R .GT. 0 ) BVN =  BVN + MVNPHI( -MAX( H, K ) )
         IF ( R .LT. 0 ) BVN = -BVN + MAX( ZERO, MVNPHI(-H)-MVNPHI(-K) )     
      ENDIF
      BVU = BVN
      END
      DOUBLE PRECISION FUNCTION MVNUNI()
*
*     Uniform (0,1) random number generator
*
*     Reference:
*     L'Ecuyer, Pierre (1996), 
*     "Combined Multiple Recursive Random Number Generators"
*     Operations Research 44, pp. 816-822.
*
*
      INTEGER A12, A13, A21, A23, P12, P13, P21, P23
      INTEGER Q12, Q13, Q21, Q23, R12, R13, R21, R23
      INTEGER X10, X11, X12, X20, X21, X22, Z, M1, M2, H 
      DOUBLE PRECISION INVMP1
      PARAMETER ( M1 = 2147483647, M2 = 2145483479 )
      PARAMETER ( A12 =   63308, Q12 = 33921, R12 = 12979 )
      PARAMETER ( A13 = -183326, Q13 = 11714, R13 =  2883 )
      PARAMETER ( A21 =   86098, Q21 = 24919, R21 =  7417 )
      PARAMETER ( A23 = -539608, Q23 =  3976, R23 =  2071 )
      PARAMETER ( INVMP1 = 4.656612873077392578125D-10 ) 
*                 INVMP1 = 1/(M1+1)
      SAVE X10, X11, X12, X20, X21, X22
      DATA       X10,      X11,      X12,      X20,      X21,      X22  
     &    / 15485857, 17329489, 36312197, 55911127, 75906931, 96210113 /      
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
      MVNUNI = Z*INVMP1
      END
