*
*   This file contains a short test program, software MVTDST for
*   the MVT distribution, plus supporting software.  The file is
*   self-contained and should compile without errors on standard
*   Fortran(77) compilers. The test program demonstrates the use
*   of MVTDSTfor computing MVT distribution values for a five
*   dimensional example problem, with four different NU values.
*
*          Alan Genz
*          Department of Mathematics
*          Washington State University
*          Pullman, WA 99164-3113
*          Email : AlanGenz@wsu.edu
*
      PROGRAM TSTMVT
*
*     Test program for MVTDST
*
      DOUBLE PRECISION ABSEPS, RELEPS, VALS, ERRS
      INTEGER N, NN, NU, I, J, IJ, MAXPTS, IFTS
      PARAMETER ( N = 5, NN = ((N-1)*N)/2, MAXPTS = 25000 )
      PARAMETER ( ABSEPS = 0, RELEPS = 0.005 )
      DOUBLE PRECISION CORREL(NN), LOW(N), UP(N), DELTA(N)
      INTEGER INFIN(N)
      DATA ( UP(I),     I = 1, N )    /  N*2D0   /
      DATA ( LOW(I),    I = 1, N )    /  N*0D0   /
      DATA ( DELTA(I),  I = 1, N )    /  N*1D0   /
      DATA ( INFIN(I),  I = 1, N )    /  N*0     /
      DATA ( CORREL(I), I = 1, NN )   / NN*0.75D0 /
      PRINT '(''        Test of MVTDST'')'
      PRINT '(5X, ''Requested Accuracy '',F8.5)', MAX(ABSEPS,RELEPS)
      PRINT '(5X,''Number of Dimensions is '',I2)',N
      PRINT '(''     Maximum # of Function Values is '',I7)', MAXPTS
*
      PRINT '(/'' I     Limits''/''    Lower  Upper  Delta'','//    
     &     ' 5X, ''Lower Left of Correlation Matrix'')'
      IJ = 0
      DO I = 1, N
         IF ( INFIN(I) .LT. 0 ) THEN 
            PRINT '(I2, '' -infin  infin '', F7.4, 7F9.4)',
     &           I, DELTA(I), ( CORREL(IJ+J), J = 1,I-1 ), 1.0
         ELSE IF ( INFIN(I) .EQ. 0 ) THEN 
            PRINT '(I2, '' -infin'', 2F7.4, 1X, 6F9.4)',
     &           I, UP(I), DELTA(I), ( CORREL(IJ+J), J = 1,I-1 ), 1.0
         ELSE IF ( INFIN(I) .EQ. 1 ) THEN 
            PRINT '(I2, F7.4, ''  infin '', F7.4, 6F9.4)',
     &           I, LOW(I), DELTA(I), ( CORREL(IJ+J), J = 1,I-1 ), 1.0
         ELSE 
            PRINT '(I2, 3F7.4, 1X, 6F9.4)', I, LOW(I), 
     &           UP(I), DELTA(I), ( CORREL(IJ+J), J = 1,I-1 ), 1.0
         ENDIF
         IJ = IJ + I-1
      END DO
      DO NU = 10, 40, 10
         PRINT '(4X,''Nu is'',I3)', NU
         CALL MVTDST( N, NU, LOW, UP, INFIN, CORREL, DELTA,  
     &        MAXPTS, ABSEPS, RELEPS, ERRS, VALS, IFTS )
         PRINT '('' Results for:  MVTDST'')'
         PRINT '('' Value:    '',2(F11.6,I4))', VALS, IFTS
         PRINT '('' Error:    '',2X,''('',F8.6'')'',3X)', ERRS
      END DO
      END
*
      SUBROUTINE MVTDST( N, NU, LOWER, UPPER, INFIN, CORREL, DELTA, 
     &                   MAXPTS, ABSEPS, RELEPS, ERROR, VALUE, INFORM )       
*
*     A subroutine for computing non-central multivariate t probabilities.
*     This subroutine uses an algorithm (QRSVN) described in the paper
*     "Methods for the  Computation of Multivariate t-Probabilities",
*        by Alan Genz and Frank Bretz
*
*          Alan Genz 
*          Department of Mathematics
*          Washington State University 
*          Pullman, WA 99164-3113
*          Email : AlanGenz@wsu.edu
*
*  Parameters
*
*     N      INTEGER, the number of variables.
*     NU     INTEGER, the number of degrees of freedom.
*            If NU < 1, then an MVN probability is computed.
*     LOWER  DOUBLE PRECISION, array of lower integration limits.
*     UPPER  DOUBLE PRECISION, array of upper integration limits.
*     INFIN  INTEGER, array of integration limits flags:
*             if INFIN(I) < 0, Ith limits are (-infinity, infinity);
*             if INFIN(I) = 0, Ith limits are (-infinity, UPPER(I)];
*             if INFIN(I) = 1, Ith limits are [LOWER(I), infinity);
*             if INFIN(I) = 2, Ith limits are [LOWER(I), UPPER(I)].
*     CORREL DOUBLE PRECISION, array of correlation coefficients; 
*            the correlation coefficient in row I column J of the 
*            correlation matrixshould be stored in 
*               CORREL( J + ((I-2)*(I-1))/2 ), for J < I.
*            The correlation matrix must be positive semi-definite.
*     DELTA  DOUBLE PRECISION, array of non-centrality parameters.
*     MAXPTS INTEGER, maximum number of function values allowed. This 
*            parameter can be used to limit the time. A sensible 
*            strategy is to start with MAXPTS = 1000*N, and then
*            increase MAXPTS if ERROR is too large.
*     ABSEPS DOUBLE PRECISION absolute error tolerance.
*     RELEPS DOUBLE PRECISION relative error tolerance.
*     ERROR  DOUBLE PRECISION estimated absolute error, 
*            with 99% confidence level.
*     VALUE  DOUBLE PRECISION estimated value for the integral
*     INFORM INTEGER, termination status parameter:
*            if INFORM = 0, normal completion with ERROR < EPS;
*            if INFORM = 1, completion with ERROR > EPS and MAXPTS 
*                           function vaules used; increase MAXPTS to 
*                           decrease ERROR;
*            if INFORM = 2, N > 100 or N < 1.
*            if INFORM = 3, correlation matrix not positive semi-definite.
*
      EXTERNAL MVFUNC
      INTEGER N, ND, NU, INFIN(*), MAXPTS, INFORM, IVLS
      DOUBLE PRECISION CORREL(*), LOWER(*), UPPER(*), DELTA(*), RELEPS, 
     &                 ABSEPS, ERROR, VALUE, MVINIT, MVFUNC
      COMMON /PTBLCK/IVLS
      IVLS = 0
      IF ( N .GT. 100 .OR. N .LT. 1 ) THEN
         INFORM = 2
         VALUE = 0
         ERROR = 1
      ELSE
         INFORM = MVINIT( N, NU, CORREL, LOWER, UPPER, DELTA, INFIN,
     &                   ND, VALUE, ERROR )
         IF ( INFORM .EQ. 0 .AND. ND .GT. 1 ) THEN
*
*             Call the lattice rule integration subroutine
*
            IF ( NU .GT. 0 ) THEN
               CALL MVKBRC( ND, IVLS, MAXPTS, MVFUNC, ABSEPS, RELEPS, 
     &                      ERROR, VALUE, INFORM )
            ELSE IF ( ND .GT. 2 ) THEN
               CALL MVKBRC( ND-1, IVLS, MAXPTS, MVFUNC, ABSEPS, RELEPS, 
     &                      ERROR, VALUE, INFORM )
            ENDIF
         ENDIF
      ENDIF
      END
*
      DOUBLE PRECISION FUNCTION MVFUNC( N, W )
*     
*     Integrand subroutine
*
      INTEGER N, NUIN, NU, INFIN(*), ND, NL, INFORM  
      DOUBLE PRECISION W(*), LOWER(*), UPPER(*), CORREL(*), DELTA(*)
      PARAMETER ( NL = 100 )
      DOUBLE PRECISION COV(NL*(NL+1)/2), A(NL),B(NL), DL(NL), Y(NL), SNU
      INTEGER INFI(NL)
      DOUBLE PRECISION MVINIT, MVCHNV, MVFNVL, MVBVN, MVSTDT, R, VL, ER
      SAVE NU, SNU, A, B, DL, INFI, COV
      IF ( NU .LT. 1 ) THEN
         R = 1
         MVFUNC = MVFNVL( N+1, W, R, DL, INFI, A, B, COV, Y )
      ELSE
         R = MVCHNV( NU, W(N) )/SNU
         MVFUNC = MVFNVL( N, W, R, DL, INFI, A, B, COV, Y )
      END IF
      RETURN
*
*     Entry point for intialization.
*
      ENTRY MVINIT( N, NUIN, CORREL, LOWER, UPPER, DELTA, INFIN, 
     &             ND, VL, ER )
*
*     Initialization and computation of covariance Cholesky factor.
*
      CALL MVSORT( N, LOWER, UPPER, DELTA, CORREL, INFIN, Y, 
     &            ND,     A,     B,    DL,    COV,  INFI, INFORM )
      MVINIT = INFORM
      IF ( INFORM .EQ. 0 ) THEN
         NU = NUIN
         VL = 1
         IF ( ND .EQ. 0 ) THEN
            ER = 0
         ELSE IF ( ND .EQ. 2 .AND. NU .LT. 1 ) THEN
            IF ( ABS( COV(3) ) .GT. 0 ) THEN
               R = SQRT( 1 + COV(2)**2 )
               IF ( INFI(2) .NE. 0 ) A(2) = ( A(2) - DL(2) )/R
               IF ( INFI(2) .NE. 1 ) B(2) = ( B(2) - DL(2) )/R
               COV(2) = COV(2)/R
               VL = MVBVN( A, B, INFI, COV(2) )
               ER = 1D-15
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
               ND = 1
            END IF
         END IF
         IF ( ND .EQ. 1 ) THEN
            IF ( INFI(1) .NE. 1 ) VL = MVSTDT( NU, B(1) ) 
            IF ( INFI(1) .NE. 0 ) VL = VL - MVSTDT( NU, A(1) ) 
            ER = 1D-16
         END IF
         IF ( NU .GT. 0 ) SNU = SQRT( DBLE(NU) )
      ELSE
         VL = 0
         ER = 1
      END IF
      END
*
      DOUBLE PRECISION FUNCTION MVFNVL( N, W, R, DL, INFI, A,B, COV, Y )
*     
*     Integrand subroutine
*
      INTEGER N, INFI(*)
      DOUBLE PRECISION W(*), R, DL(*), A(*), B(*), COV(*), Y(*)
      INTEGER I, J, IJ, IK, INFA, INFB
      DOUBLE PRECISION SUM, AI, BI, DI, EI, MVPHNV
      MVFNVL = 1
      INFA = 0
      INFB = 0
      IK = 1
      IJ = 0
      DO I = 1, N
         SUM = DL(I)
         DO J = 1, I-1
            IJ = IJ + 1
            IF ( J .LT. IK ) SUM = SUM + COV(IJ)*Y(J)
         END DO
         IF ( INFI(I) .NE. 0 ) THEN
            IF ( INFA .EQ. 1 ) THEN
               AI = MAX( AI, R*A(I) - SUM )
            ELSE
               AI = R*A(I) - SUM 
               INFA = 1
            END IF
         END IF
         IF ( INFI(I) .NE. 1 ) THEN
            IF ( INFB .EQ. 1 ) THEN
               BI = MIN( BI, R*B(I) - SUM )
            ELSE
               BI = R*B(I) - SUM 
               INFB = 1
            END IF
         END IF
         IJ = IJ + 1
         IF ( I .EQ. N .OR. COV(IJ+IK+1) .GT. 0 ) THEN 
            CALL MVLIMS( AI, BI, INFA + INFA + INFB - 1, DI, EI )
            IF ( DI .GE. EI ) THEN
               MVFNVL = 0
               RETURN
            ELSE
               MVFNVL = MVFNVL*( EI - DI )
               IF ( I .LT. N ) Y(IK) = MVPHNV( DI + W(IK)*( EI - DI ) )
               IK = IK + 1
               INFA = 0
               INFB = 0
            END IF
         END IF
      END DO
      END
*
      SUBROUTINE MVSORT( N, LOWER, UPPER, DELTA, CORREL, INFIN, Y, 
     &                  ND,     A,     B,    DL,    COV,  INFI, INFORM )
*
*     Subroutine to sort integration limits and determine Cholesky factor.
*
      INTEGER N, ND, INFIN(*), INFI(*), INFORM
      DOUBLE PRECISION     A(*),     B(*),    DL(*),    COV(*), 
     &                 LOWER(*), UPPER(*), DELTA(*), CORREL(*), Y(*)
      INTEGER I, J, K, L, M, II, IJ, IL, JMIN
      DOUBLE PRECISION SUMSQ, AJ, BJ, SUM, SQTWPI, EPS, EPSI, D, E
      DOUBLE PRECISION CVDIAG, AMIN, BMIN, DEMIN, YL, YU
      PARAMETER ( SQTWPI = 2.50662 82746 31001D0, EPS = 1D-14 )
      INFORM = 0
      IJ = 0
      II = 0
      ND = N
      DO I = 1, N
         A(I) = 0
         B(I) = 0
         DL(I) = 0
         INFI(I) = INFIN(I) 
         IF ( INFI(I) .LT. 0 ) THEN
            ND = ND - 1
         ELSE 
            IF ( INFI(I) .NE. 0 ) A(I) = LOWER(I)
            IF ( INFI(I) .NE. 1 ) B(I) = UPPER(I)
            DL(I) = DELTA(I)
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
      IF ( ND .GT. 0 ) THEN
         DO I = N, ND + 1, -1
            IF ( INFI(I) .GE. 0 ) THEN 
               DO J = 1, I-1
                  IF ( INFI(J) .LT. 0 ) THEN
                     CALL MVSWAP( J, I, A, B, DL, INFI, N, COV )
                     GO TO 10
                  ENDIF
               END DO
            ENDIF
 10      END DO
*
*     Sort remaining limits and determine Cholesky factor.
*
         II = 0
         DO I = 1, ND
*
*        Determine the integration limits for variable with minimum
*        expected probability and interchange that variable with Ith.
*
            DEMIN = 1
            JMIN = I
            CVDIAG = 0
            IJ = II
            EPSI = EPS*I*I
            DO J = I, ND
               IF ( COV(IJ+J) .GT. EPSI ) THEN
                  SUMSQ = SQRT( COV(IJ+J) )
                  SUM = DL(J) 
                  DO K = 1, I-1
                     SUM = SUM + COV(IJ+K)*Y(K)
                  END DO
                  AJ = ( A(J) - SUM )/SUMSQ
                  BJ = ( B(J) - SUM )/SUMSQ
                  CALL MVLIMS( AJ, BJ, INFI(J), D, E )
                  IF ( DEMIN .GE. E - D ) THEN
                     JMIN = J
                     AMIN = AJ
                     BMIN = BJ
                     DEMIN = E - D
                     CVDIAG = SUMSQ
                  ENDIF
               ENDIF
               IJ = IJ + J 
            END DO
            IF ( JMIN .GT. I ) THEN
               CALL MVSWAP( I, JMIN, A, B, DL, INFI, N, COV )
            END IF
            IF ( COV(II+I) .LT. -EPSI ) THEN
               INFORM = 3
               RETURN
            END IF
            COV(II+I) = CVDIAG
*
*        Compute Ith column of Cholesky factor.
*        Compute expected value for Ith integration variable and
*         scale Ith covariance matrix row and limits.
*
            IF ( CVDIAG .GT. 0 ) THEN
               IL = II + I
               DO L = I+1, ND
                  COV(IL+I) = COV(IL+I)/CVDIAG
                  IJ = II + I
                  DO J = I+1, L
                     COV(IL+J) = COV(IL+J) - COV(IL+I)*COV(IJ+I)
                     IJ = IJ + J
                  END DO
                  IL = IL + L
               END DO
               IF ( DEMIN .GT. EPSI ) THEN
                  YL = 0
                  YU = 0
                  IF ( INFI(I) .NE. 0 ) YL = -EXP( -AMIN**2/2 )/SQTWPI
                  IF ( INFI(I) .NE. 1 ) YU = -EXP( -BMIN**2/2 )/SQTWPI
                  Y(I) = ( YU - YL )/DEMIN
               ELSE
                  IF ( INFI(I) .EQ. 0 ) Y(I) = BMIN
                  IF ( INFI(I) .EQ. 1 ) Y(I) = AMIN
                  IF ( INFI(I) .EQ. 2 ) Y(I) = ( AMIN + BMIN )/2
               END IF
               DO J = 1, I
                  II = II + 1
                  COV(II) = COV(II)/CVDIAG
               END DO
                A(I) =  A(I)/CVDIAG
                B(I) =  B(I)/CVDIAG
               DL(I) = DL(I)/CVDIAG
            ELSE
               IL = II + I
               DO L = I+1, ND
                  COV(IL+I) = 0
                  IL = IL + L
               END DO
*
*        If the covariance matrix diagonal entry is zero, 
*         permute limits and rows, if necessary.
*
*
               DO J = I-1, 1, -1
                  IF ( ABS( COV(II+J) ) .GT. EPSI ) THEN
                      A(I) =  A(I)/COV(II+J)
                      B(I) =  B(I)/COV(II+J)
                     DL(I) = DL(I)/COV(II+J)
                     IF ( COV(II+J) .LT. 0 ) THEN
                        CALL MVSSWP( A(I), B(I) ) 
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
                                 CALL MVSSWP( COV(IJ-K+M), COV(IJ+M) )
                              END DO
                              CALL MVSSWP(  A(K),  A(K+1) ) 
                              CALL MVSSWP(  B(K),  B(K+1) ) 
                              CALL MVSSWP( DL(K), DL(K+1) ) 
                              M = INFI(K)
                              INFI(K) = INFI(K+1)
                              INFI(K+1) = M
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
      END
*
      SUBROUTINE MVLIMS( A, B, INFIN, LOWER, UPPER )
      DOUBLE PRECISION A, B, LOWER, UPPER, MVPHI
      INTEGER INFIN
      LOWER = 0
      UPPER = 1
      IF ( INFIN .GE. 0 ) THEN
         IF ( INFIN .NE. 0 ) LOWER = MVPHI(A)
         IF ( INFIN .NE. 1 ) UPPER = MVPHI(B)
      ENDIF
      UPPER = MAX( UPPER, LOWER )
      END      
*
      SUBROUTINE MVSSWP( X, Y )
      DOUBLE PRECISION X, Y, T
      T = X
      X = Y
      Y = T
      END
*
      SUBROUTINE MVSWAP( P, Q, A, B, D, INFIN, N, C )
*
*     Swaps rows and columns P and Q in situ, with P <= Q.
*
      DOUBLE PRECISION A(*), B(*), C(*), D(*)
      INTEGER INFIN(*), P, Q, N, I, J, II, JJ
      CALL MVSSWP( A(P), A(Q) )
      CALL MVSSWP( B(P), B(Q) )
      CALL MVSSWP( D(P), D(Q) )
      J = INFIN(P)
      INFIN(P) = INFIN(Q)
      INFIN(Q) = J
      JJ = ( P*( P - 1 ) )/2
      II = ( Q*( Q - 1 ) )/2
      CALL MVSSWP( C(JJ+P), C(II+Q) )
      DO J = 1, P-1
         CALL MVSSWP( C(JJ+J), C(II+J) )
      END DO
      JJ = JJ + P
      DO I = P+1, Q-1
         CALL MVSSWP( C(JJ+P), C(II+I) )
         JJ = JJ + I
      END DO
      II = II + Q
      DO I = Q+1, N
         CALL MVSSWP( C(II+P), C(II+Q) )
         II = II + I
      END DO
      END
*
      DOUBLE PRECISION FUNCTION MVPHI( Z )
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
      MVPHI = P
      END
      DOUBLE PRECISION FUNCTION MVPHNV(P)
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
         MVPHNV = Q*( ( ( ((((A7*R + A6)*R + A5)*R + A4)*R + A3)
     *                  *R + A2 )*R + A1 )*R + A0 )
     *            /( ( ( ((((B7*R + B6)*R + B5)*R + B4)*R + B3)
     *                  *R + B2 )*R + B1 )*R + 1 )
      ELSE
         R = MIN( P, 1 - P )
         IF ( R .GT. 0 ) THEN
            R = SQRT( -LOG(R) )
            IF ( R .LE. SPLIT2 ) THEN
               R = R - CONST2
               MVPHNV = ( ( ( ((((C7*R + C6)*R + C5)*R + C4)*R + C3)
     *                      *R + C2 )*R + C1 )*R + C0 ) 
     *                /( ( ( ((((D7*R + D6)*R + D5)*R + D4)*R + D3)
     *                      *R + D2 )*R + D1 )*R + 1 )
            ELSE
               R = R - SPLIT2
               MVPHNV = ( ( ( ((((E7*R + E6)*R + E5)*R + E4)*R + E3)
     *                      *R + E2 )*R + E1 )*R + E0 )
     *                /( ( ( ((((F7*R + F6)*R + F5)*R + F4)*R + F3)
     *                      *R + F2 )*R + F1 )*R + 1 )
            END IF
         ELSE
            MVPHNV = 9
         END IF
         IF ( Q .LT. 0 ) MVPHNV = - MVPHNV
      END IF
      END
      DOUBLE PRECISION FUNCTION MVBVN( LOWER, UPPER, INFIN, CORREL )
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
      DOUBLE PRECISION LOWER(*), UPPER(*), CORREL, MVBVU
      INTEGER INFIN(*)
      IF ( INFIN(1) .EQ. 2  .AND. INFIN(2) .EQ. 2 ) THEN
         MVBVN =  MVBVU ( LOWER(1), LOWER(2), CORREL )
     +           - MVBVU ( UPPER(1), LOWER(2), CORREL )
     +           - MVBVU ( LOWER(1), UPPER(2), CORREL )
     +           + MVBVU ( UPPER(1), UPPER(2), CORREL )
      ELSE IF ( INFIN(1) .EQ. 2  .AND. INFIN(2) .EQ. 1 ) THEN
         MVBVN =  MVBVU ( LOWER(1), LOWER(2), CORREL )
     +           - MVBVU ( UPPER(1), LOWER(2), CORREL )
      ELSE IF ( INFIN(1) .EQ. 1  .AND. INFIN(2) .EQ. 2 ) THEN
         MVBVN =  MVBVU ( LOWER(1), LOWER(2), CORREL )
     +           - MVBVU ( LOWER(1), UPPER(2), CORREL )
      ELSE IF ( INFIN(1) .EQ. 2  .AND. INFIN(2) .EQ. 0 ) THEN
         MVBVN =  MVBVU ( -UPPER(1), -UPPER(2), CORREL )
     +           - MVBVU ( -LOWER(1), -UPPER(2), CORREL )
      ELSE IF ( INFIN(1) .EQ. 0  .AND. INFIN(2) .EQ. 2 ) THEN
         MVBVN =  MVBVU ( -UPPER(1), -UPPER(2), CORREL )
     +           - MVBVU ( -UPPER(1), -LOWER(2), CORREL )
      ELSE IF ( INFIN(1) .EQ. 1  .AND. INFIN(2) .EQ. 0 ) THEN
         MVBVN =  MVBVU ( LOWER(1), -UPPER(2), -CORREL )
      ELSE IF ( INFIN(1) .EQ. 0  .AND. INFIN(2) .EQ. 1 ) THEN
         MVBVN =  MVBVU ( -UPPER(1), LOWER(2), -CORREL )
      ELSE IF ( INFIN(1) .EQ. 1  .AND. INFIN(2) .EQ. 1 ) THEN
         MVBVN =  MVBVU ( LOWER(1), LOWER(2), CORREL )
      ELSE IF ( INFIN(1) .EQ. 0  .AND. INFIN(2) .EQ. 0 ) THEN
         MVBVN =  MVBVU ( -UPPER(1), -UPPER(2), CORREL )
      END IF
      END 
      DOUBLE PRECISION FUNCTION MVBVU( SH, SK, R )
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
      DOUBLE PRECISION MVPHI, SN, ASR, H, K, BS, HS, HK
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
         BVN = BVN*ASR/(2*TWOPI) + MVPHI(-H)*MVPHI(-K) 
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
               BVN = BVN - EXP(-HK/2)*SQRT(TWOPI)*MVPHI(-B/A)*B
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
         IF ( R .GT. 0 ) BVN =  BVN + MVPHI( -MAX( H, K ) )
         IF ( R .LT. 0 ) BVN = -BVN + MAX( ZERO, MVPHI(-H) - MVPHI(-K) )     
      ENDIF
      MVBVU = BVN
      END
*
      DOUBLE PRECISION FUNCTION MVSTDT( NU, T )
*
*     Student t Distribution Function
*
*                       T
*         TSTDNT = C   I  ( 1 + y*y/NU )**( -(NU+1)/2 ) dy
*                   NU -INF
*
      INTEGER NU, J
      DOUBLE PRECISION MVPHI, T, CSTHE, SNTHE, POLYN, TT, TS, RN, PI
      PARAMETER ( PI = 3.141592653589793D0 )
      IF ( NU .LT. 1 ) THEN
         MVSTDT = MVPHI( T )
      ELSE IF ( NU .EQ. 1 ) THEN
         MVSTDT = ( 1 + 2*ATAN( T )/PI )/2
      ELSE IF ( NU .EQ. 2) THEN
         MVSTDT = ( 1 + T/SQRT( 2 + T*T ))/2
      ELSE 
         TT = T*T
         CSTHE = NU/( NU + TT )
         POLYN = 1
         DO J = NU - 2, 2, -2
            POLYN = 1 + ( J - 1 )*CSTHE*POLYN/J
         END DO
         IF ( MOD( NU, 2 ) .EQ. 1 ) THEN
            RN = NU
            TS = T/SQRT(RN)
            MVSTDT = ( 1 + 2*( ATAN( TS ) + TS*CSTHE*POLYN )/PI )/2
         ELSE
            SNTHE = T/SQRT( NU + TT )
            MVSTDT = ( 1 + SNTHE*POLYN )/2
         END IF
         IF ( MVSTDT .LT. 0 ) MVSTDT = 0
      ENDIF
      END
*
      DOUBLE PRECISION FUNCTION MVCHNV( N, P )
*
*                  MVCHNV
*     P =  1 - K  I     exp(-t*t/2) t**(N-1) dt, for N >= 1.
*               N  0
*
      INTEGER I, N, NO
      DOUBLE PRECISION P, TWO, R, RO, LRP, LKN, MVPHNV, MVCHNC, TO
      PARAMETER ( LRP = -.22579135264472743235D0, TWO = 2 )
*                 LRP =   LOG( SQRT( 2/PI ) )
      SAVE NO, LKN
      DATA NO / 0 /
      IF ( N .LE. 1 ) THEN
         R = -MVPHNV( P/2 )
      ELSE IF ( P .LT. 1 ) THEN
         IF ( N .EQ. 2 ) THEN
            R = SQRT( -2*LOG(P) )
         ELSE
            IF ( N .NE. NO ) THEN
               NO = N
               LKN = 0
               DO I = N-2, 2, -2
                  LKN = LKN - LOG( DBLE(I) )
               END DO
               IF ( MOD( N, 2 ) .EQ. 1 ) LKN = LKN + LRP
            END IF
            IF ( N .GE. -5*LOG(1-P)/4 ) THEN
               R = TWO/( 9*N )
               R = N*( -MVPHNV(P)*SQRT(R) + 1 - R )**3
               IF ( R .GT. 2*N+6 ) THEN
                  R = 2*( LKN - LOG(P) ) + ( N - 2 )*LOG(R)
               END IF
            ELSE
               R = EXP( ( LOG( (1-P)*N ) - LKN )*TWO/N )
            END IF
            R = SQRT(R)
            RO = R
            R = MVCHNC( LKN, N, P, R )
            IF ( ABS( R - RO ) .GT. 1D-6 ) THEN
               RO = R
               R = MVCHNC( LKN, N, P, R )
               IF ( ABS( R - RO ) .GT. 1D-6 ) R = MVCHNC( LKN, N, P, R )
            END IF
         END IF
      ELSE
         R = 0
      END IF
      MVCHNV = R
      END
*
      DOUBLE PRECISION FUNCTION MVCHNC( LKN, N, P, R )
*
*     Third order correction to R for MVCHNV
*
      INTEGER N, I
      DOUBLE PRECISION P, R, LKN, DF, RR, RP, PF, MVPHI
      PARAMETER ( RP = 0.79788456080286535588D0 )
      RR = R*R
      PF = 1
      DO I = N - 2, 2, -2
         PF = 1 + RR*PF/I
      END DO
      IF ( MOD( N, 2 ) .EQ. 0 ) THEN
         DF = ( P             - EXP( LOG(      PF ) - RR/2 ) )
      ELSE
         DF = ( P - 2*MVPHI(-R) - EXP( LOG( RP*R*PF ) - RR/2 ) )
      ENDIF
      MVCHNC = R - DF/( EXP(LKN+(N-1)*LOG(R)-RR/2) + DF*(R-(N-1)/R)/2 )   
      END
*
      SUBROUTINE MVKBRC( NDIM, MINVLS, MAXVLS, FUNCTN, ABSEPS, RELEPS,
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
*         Last Change: 10/15/99
*
*  MVKBRC computes an approximation to the integral
*
*      1  1     1
*     I  I ... I       F(X)  dx(NDIM)...dx(2)dx(1)
*      0  0     0
*
*
*  It uses randomized Korobov rules. The primary references are
*   "Randomization of Number Theoretic Methods for Multiple Integration"
*    R. Cranley and T.N.L. Patterson, SIAM J Numer Anal, 13, pp. 904-14,
*  and 
*   "Optimal Parameters for Multidimensional Integration", 
*    P. Keast, SIAM J Numer Anal, 10, pp.831-838.
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
      INTEGER NDIM, MINVLS, MAXVLS, INFORM, NP, PLIM, NLIM,
     &        SAMPLS, I, INTVLS, MINSMP
      PARAMETER ( PLIM = 25, NLIM = 100, MINSMP = 8 )
      INTEGER P(PLIM), C(PLIM,NLIM-1) 
      DOUBLE PRECISION FUNCTN, ABSEPS, RELEPS, FINEST, ABSERR, DIFINT, 
     &                 FINVAL, VARSQR, VAREST, VARPRD, VALUE
      DOUBLE PRECISION X(2*NLIM), VK(NLIM), ONE
      PARAMETER ( ONE = 1 )
      SAVE P, C, SAMPLS, NP, VAREST
      INFORM = 1
      INTVLS = 0
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
      DO I = 2, NDIM
         VK(I) = MOD( C(NP, NDIM-1 )*VK(I-1), ONE )
      END DO
      FINVAL = 0
      VARSQR = 0
      DO I = 1, SAMPLS
         CALL MVKSRC( NDIM, VALUE, P(NP), VK, FUNCTN, X )
         DIFINT = ( VALUE - FINVAL )/I
         FINVAL = FINVAL + DIFINT
         VARSQR = ( I - 2 )*VARSQR/I + DIFINT**2
      END DO
      INTVLS = INTVLS + 2*SAMPLS*P(NP)
      VARPRD = VAREST*VARSQR
      FINEST = FINEST + ( FINVAL - FINEST )/( 1 + VARPRD )
      IF ( VARSQR .GT. 0 ) VAREST = ( 1 + VARPRD )/VARSQR
      ABSERR = 7*SQRT( VARSQR/( 1 + VARPRD ) )/2
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
      DATA P( 1), (C( 1,I), I = 1, 99)/      31,  1*     12,  2*      9,
     &        1*     13,  8*     12,  3*      3,  1*     12,  2*      7,
     &        9*     12,  3*      3,  1*     12,  2*      7,  9*     12,
     &        3*      3,  1*     12,  2*      7,  9*     12,  3*      3,
     &        1*     12,  2*      7,  8*     12,  1*      7,  3*      3,
     &        3*      7, 21*      3/
      DATA P( 2), (C( 2,I), I = 1, 99)/      47,  1*     13,  1*     11,
     &        1*     17,  1*     10,  6*     15,  1*     22,  2*     15,
     &        3*      6,  2*     15,  1*      9,  1*     13,  3*      2,
     &        1*     13,  2*     11,  1*     10,  9*     15,  3*      6,
     &        2*     15,  1*      9,  1*     13,  3*      2,  1*     13,
     &        2*     11,  1*     10,  9*     15,  3*      6,  2*     15,
     &        1*      9,  1*     13,  3*      2,  1*     13,  2*     11,
     &        2*     10,  8*     15,  1*      6,  1*      2,  1*      3,
     &        1*      2,  1*      3, 12*      2/
      DATA P( 3), (C( 3,I), I = 1, 99)/      73,  1*     27,  1*     28,
     &        1*     10,  2*     11,  1*     20,  2*     11,  1*     28,
     &        2*     13,  1*     28,  3*     13, 16*     14,  2*     31,
     &        3*      5,  1*     31,  1*     13,  6*     11,  7*     13,
     &       16*     14,  2*     31,  3*      5,  1*     11,  1*     13,
     &        7*     11,  2*     13,  1*     11,  1*     13,  4*      5,
     &        1*     14,  1*     13,  8*      5/
      DATA P( 4), (C( 4,I), I = 1, 99)/     113,  1*     35,  2*     27,
     &        1*     36,  1*     22,  2*     29,  1*     20,  1*     45,
     &        3*      5, 16*     21,  1*     29, 10*     17, 12*     23,
     &        1*     21,  1*     27,  3*      3,  1*     24,  2*     27,
     &        1*     17,  3*     29,  1*     17,  4*      5, 16*     21,
     &        3*     17,  1*      6,  2*     17,  1*      6,  1*      3,
     &        2*      6,  5*      3/
      DATA P( 5), (C( 5,I), I = 1, 99)/     173,  1*     64,  1*     66,
     &        2*     28,  2*     44,  1*     55,  1*     67,  6*     10,
     &        2*     38,  5*     10, 12*     49,  2*     38,  1*     31,
     &        2*      4,  1*     31,  1*     64,  3*      4,  1*     64,
     &        6*     45, 19*     66,  1*     11,  9*     66,  1*     45,
     &        1*     11,  1*      7,  1*      3,  3*      2,  1*     27,
     &        1*      5,  2*      3,  2*      5,  7*      2/
      DATA P( 6), (C( 6,I), I = 1, 99)/     263,  1*    111,  1*     42,
     &        1*     54,  1*    118,  1*     20,  2*     31,  1*     72,
     &        1*     17,  1*     94,  2*     14,  1*     11,  3*     14,
     &        1*     94,  4*     10,  7*     14,  3*     11,  7*      8,
     &        5*     18,  1*    113,  2*     62,  2*     45, 17*    113,
     &        2*     63,  1*     53,  1*     63, 15*     67,  5*     51,
     &        1*     12,  1*     51,  1*     12,  1*     51,  1*      5,
     &        2*      3,  2*      2,  1*      5/
      DATA P( 7), (C( 7,I), I = 1, 99)/     397,  1*    163,  1*    154,
     &        1*     83,  1*     43,  1*     82,  1*     92,  1*    150,
     &        1*     59,  2*     76,  1*     47,  2*     11,  1*    100,
     &        1*    131,  6*    116,  9*    138, 21*    101,  6*    116,
     &        5*    100,  5*    138, 19*    101,  8*     38,  5*      3/
      DATA P( 8), (C( 8,I), I = 1, 99)/     593,  1*    246,  1*    189,
     &        1*    242,  1*    102,  2*    250,  1*    102,  1*    250,
     &        1*    280,  1*    118,  1*    196,  1*    118,  1*    191,
     &        1*    215,  2*    121, 12*     49, 34*    171,  8*    161,
     &       17*     14,  6*     10,  1*    103,  4*     10,  1*      5/
      DATA P( 9), (C( 9,I), I = 1, 99)/     907,  1*    347,  1*    402,
     &        1*    322,  1*    418,  1*    215,  1*    220,  3*    339,
     &        1*    337,  1*    218,  4*    315,  4*    167,  1*    361,
     &        1*    201, 11*    124,  2*    231, 14*     90,  4*     48,
     &       23*     90, 10*    243,  9*    283,  1*     16,  1*    283,
     &        1*     16,  2*    283/
      DATA P(10), (C(10,I), I = 1, 99)/    1361,  1*    505,  1*    220,
     &        1*    601,  1*    644,  1*    612,  1*    160,  3*    206,
     &        1*    422,  1*    134,  1*    518,  2*    134,  1*    518,
     &        1*    652,  1*    382,  1*    206,  1*    158,  1*    441,
     &        1*    179,  1*    441,  1*     56,  2*    559, 14*     56,
     &        2*    101,  1*     56,  8*    101,  7*    193, 21*    101,
     &       17*    122,  4*    101/
      DATA P(11), (C(11,I), I = 1, 99)/    2053,  1*    794,  1*    325,
     &        1*    960,  1*    528,  2*    247,  1*    338,  1*    366,
     &        1*    847,  2*    753,  1*    236,  2*    334,  1*    461,
     &        1*    711,  1*    652,  3*    381,  1*    652,  7*    381,
     &        1*    226,  7*    326,  1*    126, 10*    326,  2*    195,
     &       19*     55,  7*    195, 11*    132, 13*    387/
      DATA P(12), (C(12,I), I = 1, 99)/    3079,  1*   1189,  1*    888,
     &        1*    259,  1*   1082,  1*    725,  1*    811,  1*    636,
     &        1*    965,  2*    497,  2*   1490,  1*    392,  1*   1291,
     &        2*    508,  2*   1291,  1*    508,  1*   1291,  2*    508,
     &        4*    867,  1*    934,  7*    867,  9*   1284,  4*    563,
     &        3*   1010,  1*    208,  1*    838,  3*    563,  2*    759,
     &        1*    564,  2*    759,  4*    801,  5*    759,  8*    563,
     &       22*    226/
      DATA P(13), (C(13,I), I = 1, 99)/    4621,  1*   1763,  1*   1018,
     &        1*   1500,  1*    432,  1*   1332,  1*   2203,  1*    126,
     &        1*   2240,  1*   1719,  1*   1284,  1*    878,  1*   1983,
     &        4*    266,  2*    747,  2*    127,  1*   2074,  1*    127,
     &        1*   2074,  1*   1400, 10*   1383,  1*   1400,  7*   1383,
     &        1*    507,  4*   1073,  5*   1990,  9*    507, 17*   1073,
     &        6*     22,  1*   1073,  6*    452,  1*    318,  4*    301,
     &        2*     86,  1*     15/
      DATA P(14), (C(14,I), I = 1, 99)/    6947,  1*   2872,  1*   3233,
     &        1*   1534,  1*   2941,  1*   2910,  1*    393,  1*   1796,
     &        1*    919,  1*    446,  2*    919,  1*   1117,  7*    103,
     &        1*   2311,  1*   3117,  1*   1101,  2*   3117,  5*   1101,
     &        8*   2503,  7*    429,  3*   1702,  5*    184, 34*    105,
     &       13*    784/
      DATA P(15), (C(15,I), I = 1, 99)/   10427,  1*   4309,  1*   3758,
     &        1*   4034,  1*   1963,  1*    730,  1*    642,  1*   1502,
     &        1*   2246,  1*   3834,  1*   1511,  2*   1102,  2*   1522,
     &        2*   3427,  1*   3928,  2*    915,  4*   3818,  3*   4782,
     &        1*   3818,  1*   4782,  2*   3818,  7*   1327,  9*   1387,
     &       13*   2339, 18*   3148,  3*   1776,  3*   3354,  1*    925,
     &        2*   3354,  5*    925,  8*   2133/
      DATA P(16), (C(16,I), I = 1, 99)/   15641,  1*   6610,  1*   6977,
     &        1*   1686,  1*   3819,  1*   2314,  1*   5647,  1*   3953,
     &        1*   3614,  1*   5115,  2*    423,  1*   5408,  1*   7426,
     &        2*    423,  1*    487,  1*   6227,  1*   2660,  1*   6227,
     &        1*   1221,  1*   3811,  1*    197,  1*   4367,  1*    351,
     &        1*   1281,  1*   1221,  3*    351,  1*   7245,  1*   1984,
     &        6*   2999,  1*   3995,  4*   2063,  1*   1644,  1*   2063,
     &        1*   2077,  3*   2512,  4*   2077, 19*    754,  2*   1097,
     &        4*    754,  1*    248,  1*    754,  4*   1097,  4*    222,
     &        1*    754, 11*   1982/
      DATA P(17), (C(17,I), I = 1, 99)/   23473,  1*   9861,  1*   3647,
     &        1*   4073,  1*   2535,  1*   3430,  1*   9865,  1*   2830,
     &        1*   9328,  1*   4320,  1*   5913,  1*  10365,  1*   8272,
     &        1*   3706,  1*   6186,  3*   7806,  1*   8610,  1*   2563,
     &        2*  11558,  1*   9421,  1*   1181,  1*   9421,  3*   1181,
     &        1*   9421,  2*   1181,  2*  10574,  5*   3534,  3*   2898,
     &        1*   3450,  7*   2141, 15*   7055,  1*   2831, 24*   8204,
     &        3*   4688,  8*   2831/
      DATA P(18), (C(18,I), I = 1, 99)/   35221,  1*  10327,  1*   7582,
     &        1*   7124,  1*   8214,  1*   9600,  1*  10271,  1*  10193,
     &        1*  10800,  1*   9086,  1*   2365,  1*   4409,  1*  13812,
     &        1*   5661,  2*   9344,  1*  10362,  2*   9344,  1*   8585,
     &        1*  11114,  3*  13080,  1*   6949,  3*   3436,  1*  13213,
     &        2*   6130,  2*   8159,  1*  11595,  1*   8159,  1*   3436,
     &       18*   7096,  1*   4377,  1*   7096,  5*   4377,  2*   5410,
     &       32*   4377,  2*    440,  3*   1199/
      DATA P(19), (C(19,I), I = 1, 99)/   52837,  1*  19540,  1*  19926,
     &        1*  11582,  1*  11113,  1*  24585,  1*   8726,  1*  17218,
     &        1*    419,  3*   4918,  1*  15701,  1*  17710,  2*   4037,
     &        1*  15808,  1*  11401,  1*  19398,  2*  25950,  1*   4454,
     &        1*  24987,  1*  11719,  1*   8697,  5*   1452,  2*   8697,
     &        1*   6436,  1*  21475,  1*   6436,  1*  22913,  1*   6434,
     &        1*  18497,  4*  11089,  2*   3036,  4*  14208,  8*  12906,
     &        4*   7614,  6*   5021, 24*  10145,  6*   4544,  4*   8394/
      DATA P(20), (C(20,I), I = 1, 99)/   79259,  1*  34566,  1*   9579,
     &        1*  12654,  1*  26856,  1*  37873,  1*  38806,  1*  29501,
     &        1*  17271,  1*   3663,  1*  10763,  1*  18955,  1*   1298,
     &        1*  26560,  2*  17132,  2*   4753,  1*   8713,  1*  18624,
     &        1*  13082,  1*   6791,  1*   1122,  1*  19363,  1*  34695,
     &        4*  18770,  1*  15628,  4*  18770,  1*  33766,  6*  20837,
     &        5*   6545, 14*  12138,  5*  30483, 19*  12138,  1*   9305,
     &       13*  11107,  2*   9305/
      DATA P(21), (C(21,I), I = 1, 99)/  118891,  1*  31929,  1*  49367,
     &        1*  10982,  1*   3527,  1*  27066,  1*  13226,  1*  56010,
     &        1*  18911,  1*  40574,  2*  20767,  1*   9686,  2*  47603,
     &        2*  11736,  1*  41601,  1*  12888,  1*  32948,  1*  30801,
     &        1*  44243,  2*  53351,  1*  16016,  2*  35086,  1*  32581,
     &        2*   2464,  1*  49554,  2*   2464,  2*  49554,  1*   2464,
     &        1*     81,  1*  27260,  1*  10681,  7*   2185,  5*  18086,
     &        2*  17631,  3*  18086,  1*  37335,  3*  37774, 13*  26401,
     &        1*  12982,  6*  40398,  3*   3518,  9*  37799,  4*   4721,
     &        4*   7067/
      DATA P(22), (C(22,I), I = 1, 99)/  178349,  1*  40701,  1*  69087,
     &        1*  77576,  1*  64590,  1*  39397,  1*  33179,  1*  10858,
     &        1*  38935,  1*  43129,  2*  35468,  1*   5279,  2*  61518,
     &        1*  27945,  2*  70975,  2*  86478,  2*  20514,  2*  73178,
     &        2*  43098,  1*   4701,  2*  59979,  1*  58556,  1*  69916,
     &        2*  15170,  2*   4832,  1*  43064,  1*  71685,  1*   4832,
     &        3*  15170,  3*  27679,  2*  60826,  2*   6187,  5*   4264,
     &        1*  45567,  4*  32269,  9*  62060, 13*   1803, 12*  51108,
     &        2*  55315,  5*  54140,  1*  13134/
      DATA P(23), (C(23,I), I = 1, 99)/  267523,  1* 103650,  1* 125480,
     &        1*  59978,  1*  46875,  1*  77172,  1*  83021,  1* 126904,
     &        1*  14541,  1*  56299,  1*  43636,  1*  11655,  1*  52680,
     &        1*  88549,  1*  29804,  1* 101894,  1* 113675,  1*  48040,
     &        1* 113675,  1*  34987,  1*  48308,  1*  97926,  1*   5475,
     &        1*  49449,  1*   6850,  2*  62545,  1*   9440,  1*  33242,
     &        1*   9440,  1*  33242,  1*   9440,  1*  33242,  1*   9440,
     &        1*  62850,  3*   9440,  3*  90308,  9*  47904,  7*  41143,
     &        5*  36114,  1*  24997, 14*  65162,  7*  47650,  7*  40586,
     &        4*  38725,  5*  88329/
      DATA P(24), (C(24,I), I = 1, 99)/  401287,  1* 165843,  1*  90647,
     &        1*  59925,  1* 189541,  1*  67647,  1*  74795,  1*  68365,
     &        1* 167485,  1* 143918,  1*  74912,  1* 167289,  1*  75517,
     &        1*   8148,  1* 172106,  1* 126159,  3*  35867,  1* 121694,
     &        1*  52171,  1*  95354,  2* 113969,  1*  76304,  2* 123709,
     &        1* 144615,  1* 123709,  2*  64958,  1*  32377,  2* 193002,
     &        1*  25023,  1*  40017,  1* 141605,  2* 189165,  1* 141605,
     &        2* 189165,  3* 141605,  1* 189165, 20* 127047, 10* 127785,
     &        6*  80822, 16* 131661,  1*   7114,  1* 131661/
      DATA P(25), (C(25,I), I = 1, 99)/  601942,  1* 130365,  1* 236711,
     &        1* 110235,  1* 125699,  1*  56483,  1*  93735,  1* 234469,
     &        1*  60549,  1*   1291,  1*  93937,  1* 245291,  1* 196061,
     &        1* 258647,  1* 162489,  1* 176631,  1* 204895,  1*  73353,
     &        1* 172319,  1*  28881,  1* 136787,  2* 122081,  1* 275993,
     &        1*  64673,  3* 211587,  2* 282859,  1* 211587,  1* 242821,
     &        3* 256865,  1* 122203,  1* 291915,  1* 122203,  2* 291915,
     &        1* 122203,  2*  25639,  1* 291803,  1* 245397,  1* 284047,
     &        7* 245397,  1*  94241,  2*  66575, 19* 217673, 10* 210249,
     &       15*  94453/
*
      END
*
      SUBROUTINE MVKSRC( NDIM, SUMKRO, PRIME, VK, FUNCTN, X )
      EXTERNAL FUNCTN
      INTEGER NDIM, PRIME, K, J, JP
      DOUBLE PRECISION SUMKRO, VK(*), FUNCTN, X(*), ONE, XT, MVUNI, FT
      PARAMETER ( ONE = 1 )
      SUMKRO = 0
      DO J = 1, NDIM-1
         JP = J + MVUNI()*( NDIM + 1 - J ) 
         XT = VK(J)
         VK(J) = VK(JP)
         VK(JP) = XT
      END DO
      DO J = 1, NDIM
         X(NDIM+J) = MVUNI()
      END DO
      DO K = 1, PRIME
         DO J = 1, NDIM
            X(J) = MOD( K*VK(J), ONE )
         END DO
         DO J = 1, NDIM
            XT = X(J) + X(NDIM+J)
            IF ( XT .GT. 1 ) XT = XT - 1
            X(J) = ABS( 2*XT - 1 )
         END DO
         FT =  FUNCTN(NDIM,X)
         DO J = 1, NDIM
            X(J) = 1 - X(J)
         END DO
         FT = ( FT + FUNCTN(NDIM,X) )/2
         SUMKRO = SUMKRO + ( FT - SUMKRO )/K
      END DO
      END
*
      DOUBLE PRECISION FUNCTION MVUNI()
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
*                 INVMP1 = 1/( M1 + 1 )
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
      MVUNI = Z*INVMP1
      END
