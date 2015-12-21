*
* This file contains a short test program, software for MVT distribution 
* subroutines SADMVT, RANMVT, SPHMVT and KROMVT, plus supporting software.
* The file is self contained and should compile without errors on (77) 
* standard Fortran compilers. The test program demonstrates the use of
* four different methods for computing MVT distribution values for a
* five dimensional example problem, with three different NU values.
*
*          Alan Genz
*          Department of Mathematics
*          Washington State University
*          Pullman, WA 99164-3113
*          Email : AlanGenz@wsu.edu
*
      PROGRAM TSTMVT
*
*     Test program for SADMVT, RANMVT and KROMVT
*
      DOUBLE PRECISION ABSEPS,RELEPS, VALR,ERRR, VALS,ERRS, VALK,ERRK
      INTEGER N, NN, NU, I, J, IJ, MAXPTS, IFTR, IFTS, IFTK
      PARAMETER ( N = 5, NN = ((N-1)*N)/2, MAXPTS = 25000 )
      PARAMETER ( ABSEPS = 0, RELEPS = 0.005 )
      DOUBLE PRECISION CORREL(NN), LOW(N), UP(N)
      INTEGER INFIN(N)
      DATA (  UP(I),    I = 1, N )    /  N*2D0   /
      DATA ( LOW(I),    I = 1, N )    /  N*0D0   /
      DATA ( INFIN(I),  I = 1, N )    /  N*0     /
      DATA ( CORREL(I), I = 1, NN )   / NN*0.75D0 /
      PRINT '(''        Test of SADMVT, RANMVT and KROMVT'')'
      PRINT '(12X, ''Requested Accuracy '',F8.5)', MAX(ABSEPS,RELEPS)
      PRINT '(7X,''Number of Dimensions is '',I2)',N
      PRINT '(''     Maximum # of Function Values is '',I7)', MAXPTS
*
      PRINT '(/'' I     Limits''/''    Lower  Upper'','//    
     &     ' 5X, ''Lower Left of Correlation Matrix'')'
      IJ = 0
      DO I = 1, N
         IF ( INFIN(I) .LT. 0 ) THEN 
            PRINT '(I2, '' -infin  infin '', 7F9.4)',
     &           I, ( CORREL(IJ+J), J = 1,I-1 ), 1.0
         ELSE IF ( INFIN(I) .EQ. 0 ) THEN 
            PRINT '(I2, '' -infin'', F7.4, 1X, 6F9.4)',
     &           I, UP(I), ( CORREL(IJ+J), J = 1,I-1 ), 1.0
         ELSE IF ( INFIN(I) .EQ. 1 ) THEN 
            PRINT '(I2, F7.4, ''  infin '', 6F9.4)',
     &           I, LOW(I), ( CORREL(IJ+J), J = 1,I-1 ), 1.0
         ELSE 
            PRINT '(I2, 2F7.4, 1X, 6F9.4)', 
     &           I, LOW(I), UP(I), ( CORREL(IJ+J), J = 1,I-1 ), 1.0
         ENDIF
         IJ = IJ + I-1
      END DO
      DO NU = 10, 40, 10
         PRINT '(4X,''Nu is'',I3)', NU
         CALL SADMVT( N, NU, LOW, UP, INFIN, CORREL, 
     &        MAXPTS, ABSEPS, RELEPS, ERRS, VALS, IFTS )
         CALL KROMVT( N, NU, LOW, UP, INFIN, CORREL, 
     &        MAXPTS, ABSEPS, RELEPS, ERRK, VALK, IFTK )
         CALL RANMVT( N, NU, LOW, UP, INFIN, CORREL, 
     &        MAXPTS, ABSEPS, RELEPS, ERRR, VALR, IFTR )
         PRINT '('' Results for:  RANMVT'',9X,''SADMVT'',9X,''KROMVT'')'
         PRINT '('' Values:   '',3(F11.6,I4))',
     &        VALR, IFTR, VALS, IFTS, VALK, IFTK
         PRINT '('' Errors:   '',3(2X,''('',F8.6'')'',3X))',
     &        ERRR, ERRS, ERRK
      END DO
      END
*
      SUBROUTINE RANMVT(N, NU, LOWER, UPPER, INFIN, CORREL, MAXPTS,
     *      ABSEPS, RELEPS, ERROR, VALUE, INFORM)
*
*     A subroutine for computing multivariate t probabilities.
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
*            parameter can be used to limit the time taken. A sensible 
*            strategy is to start with MAXPTS = 1000*N, and then
*            increase MAXPTS if ERROR is too large.
*     ABSEPS REAL absolute error tolerance.
*     RELEPS REAL relative error tolerance.
*     ERROR  REAL, estimated absolute error, with 99% confidence level.
*     VALUE  REAL, estimated value for the integral
*     INFORM INTEGER, termination status parameter:
*            if INFORM = 0, normal completion with ERROR < EPS;
*            if INFORM = 1, completion with ERROR > EPS and MAXPTS 
*                           function vaules used; increase MAXPTS to 
*                           decrease ERROR;
*            if INFORM = 2, N > 20 or N < 1.
*
      EXTERNAL FNCMVT
      INTEGER N, NU, INFIN(*), MAXPTS, INFORM, INFIS, MPT, IVLS
      DOUBLE PRECISION CORREL(*), LOWER(*), UPPER(*), 
     *     ABSEPS, RELEPS, EPS, ERROR, VALUE, E, D, MVTNIT
      IF ( N .GT. 20 .OR. N .LT. 1 ) THEN
         INFORM = 2
         VALUE = 0
         ERROR = 1
         RETURN
      ENDIF
      INFORM = MVTNIT(N, NU, CORREL, LOWER, UPPER, INFIN, INFIS, D, E)
      IF ( N-INFIS .EQ. 0 ) THEN
         VALUE = 1
         ERROR = 0
      ELSE IF ( N-INFIS .EQ. 1 ) THEN
         VALUE = E - D
         ERROR = 2E-16
      ELSE
*
*        Call then Monte-Carlo integration subroutine
*
         MPT = 25 + 10*N*N
         CALL RCRUDE(N-INFIS-1, MPT, FNCMVT, ERROR, VALUE, 0)
         IVLS = MPT
 10      EPS = MAX( ABSEPS, RELEPS*ABS(VALUE) )
         IF ( ERROR .GT. EPS .AND. IVLS .LT. MAXPTS ) THEN
            MPT = MAX(MIN( INT(MPT*(ERROR/(EPS))**2), MAXPTS-IVLS ), 10)
            CALL RCRUDE(N-INFIS-1, MPT, FNCMVT, ERROR, VALUE, 1)
            IVLS = IVLS + MPT
            GO TO 10
         ENDIF
         IF ( ERROR. GT. EPS .AND. IVLS .GE. MAXPTS ) INFORM = 1
      ENDIF
      END
      SUBROUTINE RCRUDE( NDIM, MAXPTS, FUNCTN, ABSEST, FINEST, IR )
*
*     Crude Monte-Carlo Algorithm with simple antithetic variates
*      and weighted results on restart
*
      EXTERNAL FUNCTN
      INTEGER NDIM, MAXPTS, M, K, IR, NPTS
      DOUBLE PRECISION FINEST, ABSEST, X(100), FUN, FUNCTN, UNI,
     *     VARSQR, VAREST, VARPRD, FINDIF, FINVAL
      SAVE VAREST
      IF ( IR .LE. 0 ) THEN
         VAREST = 0
         FINEST = 0
      ENDIF
      FINVAL = 0
      VARSQR = 0
      NPTS = MAXPTS/2
      DO M = 1,NPTS
         DO K = 1,NDIM
            X(K) = UNI()
         END DO
         FUN = FUNCTN(NDIM, X)
         DO K = 1,NDIM
            X(K) = 1 - X(K)
         END DO
         FUN = ( FUNCTN(NDIM, X) + FUN )/2
         FINDIF = ( FUN - FINVAL )/M
         VARSQR = ( M - 2 )*VARSQR/M + FINDIF**2
         FINVAL = FINVAL + FINDIF
      END DO
      VARPRD = VAREST*VARSQR
      FINEST = FINEST + ( FINVAL - FINEST )/(1 + VARPRD)
      IF ( VARSQR .GT. 0 ) VAREST = (1 + VARPRD)/VARSQR
      ABSEST = 3*SQRT( VARSQR/( 1 + VARPRD ) )
      END
*
      SUBROUTINE KROMVT(N, NU, LOWER, UPPER, INFIN, CORREL, MAXPTS,
     *      ABSEPS, RELEPS, ERROR, VALUE, INFORM)
*
*     A subroutine for computing multivariate t probabilities.
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
*     ABSEPS    REAL absolute error tolerance.
*     RELEPS    REAL relative error tolerance.
*     ERROR  REAL estimated absolute error, with 99% confidence level.
*     VALUE  REAL estimated value for the integral
*     INFORM INTEGER, termination status parameter:
*            if INFORM = 0, normal completion with ERROR < EPS;
*            if INFORM = 1, completion with ERROR > EPS and MAXPTS 
*                           function vaules used; increase MAXPTS to 
*                           decrease ERROR;
*            if INFORM = 2, N > 20 or N < 1.
*
      EXTERNAL FNCMVT
      INTEGER N, NU, INFIN(*), MAXPTS, INFORM, INFIS, IVLS
      DOUBLE PRECISION
     *     CORREL(*), LOWER(*), UPPER(*), RELEPS, ABSEPS,
     *     ERROR, VALUE, E, D, MVTNIT
      IF ( N .GT. 20 .OR. N .LT. 1 ) THEN
         INFORM = 2
         VALUE = 0
         ERROR = 1
         RETURN
      ENDIF
      INFORM = MVTNIT( N, NU, CORREL, LOWER, UPPER, INFIN, INFIS, D, E )
      IF ( N-INFIS .EQ. 0 ) THEN
         VALUE = 1
         ERROR = 0
      ELSE IF ( N-INFIS .EQ. 1 ) THEN
         VALUE = E - D
         ERROR = 2E-16
      ELSE
*
*        Call the lattice rule integration integration subroutine
*
         IVLS = 0
         CALL KROBOV( N-INFIS-1, IVLS, MAXPTS, FNCMVT, ABSEPS, RELEPS, 
     *                ERROR, VALUE, INFORM )
      ENDIF
      END
      SUBROUTINE KROBOV( NDIM, MINVLS, MAXVLS, FUNCTN, ABSEPS, RELEPS,
     *                   ABSERR, FINEST, IFAIL )
*
*  Automatic Multidimensional Integration Subroutine
*               
*         AUTHOR: Alan Genz
*                 Department of Mathematics
*                 Washington State University
*                 Pulman, WA 99164-3113
*
*         Last Change: 5/15/93
*
*  KROBOV computes an approximation to the integral
*
*      1  1     1
*     I  I ... I       F(X)  dx(NDIM)...dx(2)dx(1)
*      0  0     0
*
*
*  KROBOV uses randomized Korobov rules. The primary references are
*  "Randomization of Number Theoretic Methods for Multiple Integration"
*   R. Cranley and T.N.L. Patterson, SIAM J Numer Anal, 13, pp. 904-14,
*  and 
*   "Optimal Parameters for Multidimensional Integration", 
*    P. Keast, SIAM J Numer Anal, 10, pp.831-838.
*   
***************  Parameters for KROBOV  *******************************
****** Input parameters
*  NDIM    Number of variables, must exceed 1, but not exceed 40
*  MINVLS  Integer minimum number of function evaluations allowed.
c          MINVLS must not exceed MAXVLS.  If MINVLS < 0 then the
c          routine assumes a previous call of KROBOV has been made with 
c          the same integrand and continues that calculation.
*  MAXVLS  Integer maximum number of function evaluations allowed.
*  FUNCTN  EXTERNALly declared user defined function to be integrated.
c          It must have parameters (NDIM,Z), where Z is a real array
c          of dimension NDIM.
*  ABSEPS  Required absolute accuracy.
*  RELEPS  Required relative accuracy.
****** Output parameters
*  MINVLS  Actual number of function evaluations used by KROBOV.
*  ABSERR  Estimated absolute accuracy of FINEST.
*  FINEST  Estimated value of integral.
*  IFAIL   IFAIL = 0 for normal exit, when 
*                     ABSERR <= MAX(ABSEPS, RELEPS*ABS(FINEST))
*                  and 
*                     INTVLS <= MAXCLS.
*          IFAIL = 1 If MAXVLS was too small for KROBOV to obtain the
c                  required accuracy. In this case KROBOV returns a
c                  value FINEST with estimated absolute accuracy ABSERR.
************************************************************************
      EXTERNAL FUNCTN
      INTEGER PLIM, NLIM, NDIM, MINVLS, MAXVLS, NP, 
     *     IFAIL, SAMPLS, LPSAMP, I, INTVLS, MINSMP
      PARAMETER ( PLIM = 17, NLIM = 40, MINSMP = 16 )
      DOUBLE PRECISION FUNCTN, ABSEPS, RELEPS, FINEST, ABSERR, DIFINT, 
     *     FINVAL, VARSQR, VAREST, VARPRD, VALUE
      DOUBLE PRECISION ALPHA(NLIM), X(NLIM)
      INTEGER C(PLIM,NLIM), P(PLIM)
      DOUBLE PRECISION VK(NLIM), ONE
      PARAMETER ( ONE = 1 )
      SAVE P, C, SAMPLS, NP
      DATA P( 1),(C( 1,I), I = 1,39) /  173,
     *    73,   34,   57,    9,   12,    2,   16,   30,   30,   42,
     *    70,   86,    2,   53,   53,   30,   30,    5,   42,   42,
     *    70,   42,   53,   42,   42,   53,   42,   53,   53,    2,
     *    86,    2,    2,    2,    2,    2,    2,    2,    2/
      DATA P( 2),(C( 2,I), I = 1,39) /  263,
     *   111,  106,   51,   36,   48,  110,    2,    2,    2,    2,
     *    70,   70,   48,    2,    2,   70,  124,  124,   70,   48,
     *    48,   48,   48,  108,   65,   48,   48,   70,    2,   20,
     *     2,    2,    2,    2,    2,    2,    2,    2,    2/
      DATA P( 3),(C( 3,I), I = 1,39) /  397,
     *   163,  168,  164,  133,   23,   64,    2,    2,  106,   80,
     *    80,  126,   16,   16,   16,   16,   16,   16,  107,   80,
     *     2,    2,    2,   32,   32,   32,   31,   64,   31,   31,
     *     4,    4,    4,  126,   16,   16,   16,   16,   16/
      DATA P( 4),(C( 4,I), I = 1,39) /  593,
     *   229,   40,  268,  240,   31,  119,   71,  296,  130,  199,
     *   149,  149,  149,  149,  149,   31,  130,  149,  149,   79,
     *   119,  119,   31,   82,  130,  122,  122,  122,  122,    2,
     *   130,  130,  130,  130,    2,    2,   82,   82,    2/
      DATA P( 5),(C( 5,I), I = 1,39) /  887,
     *   192,  424,   55,  221,  179,  242,  242,    2,    2,   11,
     *    11,   11,  394,  394,  439,  394,  394,  394,  394,  439,
     *   394,  394,  394,  101,  378,  394,  394,  394,  394,  394,
     *   202,  279,  394,  279,    2,    2,    2,    2,    2/
      DATA P( 6),(C( 6,I), I = 1,39) / 1327,
     *   513,  195,  599,  661,  443,  632,  251,  603,  663,    2,
     *   425,  425,  603,  425,  425,  525,  412,  412,  412,  412,
     *   412,   82,   82,   82,  603,  580,  580,  444,   82,   82,
     *   276,  601,  276,  276,  276,  276,  112,  112,  112/
      DATA P( 7),(C( 7,I), I = 1,39) / 1997,
     *   839,  146,  860,  183,  121,   11,   11,  793,  998,    2,
     *     2,  110,  110,  236,  110,  236,  147,  147,  110,  190,
     *   147,  147,  147,  147,  147,  147,  236,  110,  110,  147,
     *   110,  110,  632,  147,  147,  148,    2,  147,  147/
      DATA P( 8),(C( 8,I), I = 1,39) / 2999,
     *  1148, 1406, 1192, 1094, 1290,  632,  341,  785,  393, 1499,
     *     2,  798,  808,  798,  918,  393,  924,  924,    2,    2,
     *     2, 1499,    2, 1016,  798,  798,  798,  808,  270, 1344,
     *   798,  798,  798,  798,  798,  798,  798,    2,    2/
      DATA P( 9),(C( 9,I), I = 1,39) / 4493,
     *  1360,  383,  842, 2157,   30,  959,    3,  717, 1107,    2,
     *     2,    2,  836,  836, 1134,  836,  836,  426,  898,  898,
     *    65,  836,  836,  836,  836,  216,  104,  300,  836, 1022,
     *  1022, 1022, 1022, 1420, 1478, 1478, 1478,  283, 2246/
      DATA P(10),(C(10,I), I = 1,39) / 6737,
     *  2602, 2818, 3240, 2528, 2260, 3141, 2857, 1484, 2113, 2265,
     *     2,    2, 2207, 2207, 2207,  542,  132,  934,  378,  378,
     *  2099,  934,  225,  225,  225,  169,  378, 2257, 2257, 2257,
     *  2257,  934, 2576,  934,  934,  934,  934,  934, 2257/
      DATA P(11),(C(11,I), I = 1,39) / 10111,
     *  3071, 3114, 1170, 3432, 2726, 1098, 3371,  185,    4, 3143,
     *  5055,    2,    2,    2,    2,    2,    2,  334, 1254, 4146,
     *   617, 1879,    2,    2, 1146,  475, 4725,    2,    2,  475,
     *   475,  475,  475,  475,  638,  638,  638,    2, 3107/
      DATA P(12), (C(12,I), I = 1,39)/  15161,
     *   6280,  2350,  1452,  7009,  4273,  4273,  2538,  3976,  4273,
     *   2716,  2716,  6143,  6143,    91,   107,  2716,     2,     2,
     *      2,     2,     2,     2,     2,     2,     2,     2,     2,
     *      2,     2,     2,     2,     2,     2,     2,     2,     2,
     *      2,     2,     2/
      DATA P(13), (C(13,I), I = 1,39)/ 22751,
     *   9438,  3557,   817,  9167,  7379,   811,   811,  6077,  6205,
     *   2505,  2896,  4055,  4055,  4055,   900,  1333,  1154,     2,
     *      2,     2,     2,     2,     2,     2,     2,     2,     2,
     *      2,     2,     2,     2,     2,     2,     2,     2,     2,
     *      2,     2,     2/
      DATA P(14), (C(14,I), I = 1,39)/ 34127,
     *  12543,  7726,  6142, 11635,  6071,  6071,  6071,  6071,  1886,
     *   1993,   634,   836,  1993,  2095,   142,  6542, 13195,     4,
     *      4,     2,     2,     2,     2,     2,     2,     2,     2,
     *      2,     2,     2,     2,     2,     2,     2,     2,     2,
     *      2,     2,     2/
      DATA P(15), (C(15,I), I = 1,39)/ 51193,
     *  19397, 10691, 20935, 12665,  5132,  5132,  5132,  5132,  5132,
     *   5132,  4347,  1087,  2451,  2134, 10539,  3595,  3892,  6301,
     *   3595,     2,     2,     2,     2,     2,     2,     2,     2,
     *      2,     2,     2,     2,     2,     2,     2,     2,     2,
     *      2,     2,     2/
      DATA P(16), (C(16,I), I = 1,39)/ 76801,
     *  22500, 16187, 33739, 23268, 23268, 23268, 23268, 23268, 23268,
     *  20841, 14733,  5916,  4270, 13443, 13443, 21793, 21793,   871,
     *  18819,     2,     2,     2,     2,     2,     2,     2,     2,
     *      2,     2,     2,     2,     2,     2,     2,     2,     2,
     *      2,     2,     2/
      DATA P(17), (C(17,I), I = 1,39)/115183,
     *  33681, 26530, 42347, 15491, 15491, 15491, 15491, 15491, 24490,
     *   5708, 24490, 34535,   204, 23709,  5474, 14302, 38969, 11298,
     *      2, 22095,  7317,     2,     2,     2,     2,     2,     2,
     *      2,     2,     2,     2,     2,     2,     2,     2,     2,
     *      2,     2,     2/
      IFAIL = 1
      INTVLS = 0
      FINEST = 0
      VAREST = 0
      IF ( MINVLS .GE. 0 ) THEN
         FINEST = 0
         SAMPLS = MINSMP 
         DO I = 1,PLIM
            NP = I
            IF ( MINVLS .LT. SAMPLS*P(I) ) GO TO 10
         END DO
         SAMPLS = MAX( MINSMP, MINVLS/P(NP) )
      ENDIF
 10   VK(1) = ONE/P(NP)
      DO I = 2,NDIM
         VK(I) = MOD( C(NP,NDIM-1)*VK(I-1), ONE )
      END DO
      FINVAL = 0
      VARSQR = 0
      LPSAMP = SAMPLS/2
      DO I = 1, LPSAMP
         CALL KROSUM( NDIM, VALUE, P(NP), VK, FUNCTN, ALPHA, X )
         DIFINT = ( VALUE - FINVAL )/I
         FINVAL = FINVAL + DIFINT
         VARSQR = ( I - 2 )*VARSQR/I + DIFINT**2
      END DO
      INTVLS = INTVLS + 2*LPSAMP*P(NP)
      VARPRD = VAREST*VARSQR
      FINEST = FINEST + ( FINVAL - FINEST )/( 1 + VARPRD )
      IF ( VARSQR .GT. 0 ) VAREST = ( 1 + VARPRD )/VARSQR
      ABSERR = 3*SQRT( VARSQR/( 1 + VARPRD ) )
      IF ( ABSERR .GT. MAX( ABSEPS, ABS(FINEST)*RELEPS ) ) THEN
         IF ( NP .LT. PLIM ) THEN
            NP = NP + 1
         ELSE
            SAMPLS = MAX( MINSMP,MIN(3*SAMPLS/2,(MAXVLS-INTVLS)/P(NP)) ) 
         ENDIF
         IF ( INTVLS+SAMPLS*P(NP) .LE. MAXVLS ) GO TO 10
      ELSE
         IFAIL = 0
      ENDIF
      END
      SUBROUTINE KROSUM( NDIM, SUMKRO, NPTS, VK, FUNCTN, ALPHA, X )
      EXTERNAL FUNCTN
      DOUBLE PRECISION VK(*), ONE
      INTEGER NDIM, NPTS, K, J
      PARAMETER ( ONE = 1 )
      DOUBLE PRECISION ALPHA(*), X(*), SUMKRO
      DOUBLE PRECISION SUMFUN, UNI, FUNCTN
      SUMFUN = 0
      DO J = 1,NDIM
         ALPHA(J) = UNI()
      END DO
      DO K = 1,NPTS
         DO J = 1,NDIM
            X(J) = ABS( 2*MOD( ALPHA(J) + VK(J)*K, ONE ) - 1 )
         END DO
         SUMFUN = SUMFUN + FUNCTN(NDIM,X)/2
         DO J = 1,NDIM
            X(J) = 1 - X(J)
         END DO
         SUMFUN = SUMFUN + FUNCTN(NDIM,X)/2
      END DO
      SUMKRO = SUMFUN/NPTS
      END
*
      SUBROUTINE SADMVT(N, NU, LOWER, UPPER, INFIN, CORREL, MAXPTS,
     *      ABSEPS, RELEPS, ERROR, VALUE, INFORM)
*
*     A subroutine for computing multivariate t probabilities.
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
*            parameter can be used to limit the time taken. A sensible 
*            strategy is to start with MAXPTS = 1000*N, and then
*            increase MAXPTS if ERROR is too large.
*     ABSEPS REAL absolute error tolerance.
*     RELEPS REAL relative error tolerance.
*     ERROR  REAL, estimated absolute error, with 99% confidence level.
*     VALUE  REAL, estimated value for the integral
*     INFORM INTEGER, termination status parameter:
*            if INFORM = 0, normal completion with ERROR < EPS;
*            if INFORM = 1, completion with ERROR > EPS and MAXPTS 
*                           function vaules used; increase MAXPTS to 
*                           decrease ERROR;
*            if INFORM = 2, N > 20 or N < 1.
*
      EXTERNAL FNCMVT
      INTEGER NL, N, NU, M, INFIN(*), LENWRK, MAXPTS, INFORM, INFIS,
     &     RULCLS, TOTCLS, NEWCLS, MAXCLS
      DOUBLE PRECISION CORREL(*), LOWER(*), UPPER(*), ABSEPS, RELEPS, 
     &     ERROR, VALUE, OLDVAL, D, E, MVTNIT
      PARAMETER ( NL = 20 )
      PARAMETER ( LENWRK = 20*NL**2 )
      DOUBLE PRECISION WORK(LENWRK)
      IF ( N .GT. 20 .OR. N .LT. 1 ) THEN
         INFORM = 2
         VALUE = 0
         ERROR = 1
         RETURN
      ENDIF
      INFORM = MVTNIT( N, NU, CORREL, LOWER, UPPER, INFIN, INFIS, D, E )
      M = N - INFIS
      IF ( M .EQ. 0 ) THEN
         VALUE = 1
         ERROR = 0
      ELSE IF ( M .EQ. 1 ) THEN
         VALUE = E - D
         ERROR = 2E-16
      ELSE
*
*        Call the subregion adaptive integration subroutine
*
         M = M - 1
         RULCLS = 1
         CALL ADAPT( M, RULCLS, 0, FNCMVT, ABSEPS, RELEPS,
     *               LENWRK, WORK, ERROR, VALUE, INFORM )
         MAXCLS = MIN( 10*RULCLS, MAXPTS )
         TOTCLS = 0
         CALL ADAPT( M, TOTCLS, MAXCLS, FNCMVT, ABSEPS, RELEPS,
     *               LENWRK, WORK, ERROR, VALUE, INFORM )
         IF ( ERROR .GT. MAX( ABSEPS, RELEPS*ABS(VALUE) ) ) THEN
 10         OLDVAL = VALUE
            MAXCLS = MAX( 2*RULCLS, MIN( 3*MAXCLS/2, MAXPTS - TOTCLS ) )
            NEWCLS = -1
            CALL ADAPT( M, NEWCLS, MAXCLS, FNCMVT, ABSEPS, RELEPS,
     *                  LENWRK, WORK, ERROR, VALUE, INFORM  )
            TOTCLS = TOTCLS + NEWCLS
            ERROR = ABS(VALUE-OLDVAL) + SQRT(RULCLS*ERROR**2/TOTCLS)
            IF ( ERROR .GT. MAX( ABSEPS, RELEPS*ABS(VALUE) ) ) THEN
               IF ( MAXPTS - TOTCLS .GT. 2*RULCLS ) GO TO 10
            ELSE
               INFORM = 0
            END IF
         ENDIF
      ENDIF
      END
*
      SUBROUTINE ADAPT(NDIM, MINCLS, MAXCLS, FUNCTN,
     &     ABSREQ, RELREQ, LENWRK, WORK, ABSEST, FINEST, INFORM)
*
*   Adaptive Multidimensional Integration Subroutine
*
*   Author: Alan Genz
*           Department of Mathematics
*           Washington State University
*           Pullman, WA 99164-3113 USA
*
*  This subroutine computes an approximation to the integral
*
*      1 1     1
*     I I ... I       FUNCTN(NDIM,X)  dx(NDIM)...dx(2)dx(1)
*      0 0     0  
*
***************  Parameters for ADAPT  ********************************
*
****** Input Parameters
*
*  NDIM    Integer number of integration variables.
*  MINCLS  Integer minimum number of FUNCTN calls to be allowed; MINCLS
*          must not exceed MAXCLS. If MINCLS < 0, then ADAPT assumes
*          that a previous call of ADAPT has been made with the same
*          integrand and continues that calculation.
*  MAXCLS  Integer maximum number of FUNCTN calls to be used; MAXCLS
*          must be >= RULCLS, the number of function calls required for
*          one application of the basic integration rule.
*           IF ( NDIM .EQ. 1 ) THEN
*              RULCLS = 11
*           ELSE IF ( NDIM .LT. 15 ) THEN
*              RULCLS = 2**NDIM + 2*NDIM*(NDIM+3) + 1
*           ELSE
*              RULCLS = 1 + NDIM*(24-NDIM*(6-NDIM*4))/3
*           ENDIF
*  FUNCTN  Externally declared real user defined integrand. Its 
*          parameters must be (NDIM, Z), where Z is a real array of
*          length NDIM.
*  ABSREQ  Real required absolute accuracy.
*  RELREQ  Real required relative accuracy.
*  LENWRK  Integer length of real array WORK (working storage); ADAPT
*          needs LENWRK >= 16*NDIM + 27. For maximum efficiency LENWRK
*          should be about 2*NDIM*MAXCLS/RULCLS if MAXCLS FUNCTN
*          calls are needed. If LENWRK is significantly less than this,
*          ADAPT may be less efficient.
*
****** Output Parameters
*
*  MINCLS  Actual number of FUNCTN calls used by ADAPT.
*  WORK    Real array (length LENWRK) of working storage. This contains
*          information that is needed for additional calls of ADAPT
*          using the same integrand (input MINCLS < 0).
*  ABSEST  Real estimated absolute accuracy.
*  FINEST  Real estimated value of integral.
*  INFORM  INFORM = 0 for normal exit, when ABSEST <= ABSREQ or
*                     ABSEST <= |FINEST|*RELREQ with MINCLS <= MAXCLS.
*          INFORM = 1 if MAXCLS was too small for ADAPT to obtain the
*                     result FINEST to within the requested accuracy.
*          INFORM = 2 if MINCLS > MAXCLS, LENWRK < 16*NDIM + 27 or 
*                     RULCLS > MAXCLS.
*
************************************************************************
*
*     Begin driver routine. This routine partitions the working storage 
*      array and then calls the main subroutine ADBASE.
*
      EXTERNAL FUNCTN
      INTEGER NDIM, MINCLS, MAXCLS, LENWRK, INFORM
      DOUBLE PRECISION 
     &     FUNCTN, ABSREQ, RELREQ, WORK(LENWRK), ABSEST, FINEST
      INTEGER SBRGNS, MXRGNS, RULCLS, LENRUL, 
     & INERRS, INVALS, INPTRS, INLWRS, INUPRS, INMSHS, INPNTS, INWGTS, 
     & INLOWR, INUPPR, INWDTH, INMESH, INWORK 
      IF ( NDIM .EQ. 1 ) THEN
         LENRUL = 5
         RULCLS = 9
      ELSE IF ( NDIM .LT. 12 ) THEN
         LENRUL = 6
         RULCLS = 2**NDIM + 2*NDIM*(NDIM+2) + 1
      ELSE
         LENRUL = 6
         RULCLS = 1 + 2*NDIM*(1+2*NDIM)
      ENDIF
      IF ( LENWRK .GE. LENRUL*(NDIM+4) + 10*NDIM + 3 .AND.
     &     RULCLS. LE. MAXCLS .AND. MINCLS .LE. MAXCLS ) THEN
        MXRGNS = ( LENWRK - LENRUL*(NDIM+4) - 7*NDIM )/( 3*NDIM + 3 )
        INERRS = 1
        INVALS = INERRS + MXRGNS
        INPTRS = INVALS + MXRGNS
        INLWRS = INPTRS + MXRGNS
        INUPRS = INLWRS + MXRGNS*NDIM
        INMSHS = INUPRS + MXRGNS*NDIM
        INWGTS = INMSHS + MXRGNS*NDIM
        INPNTS = INWGTS + LENRUL*4
        INLOWR = INPNTS + LENRUL*NDIM
        INUPPR = INLOWR + NDIM
        INWDTH = INUPPR + NDIM
        INMESH = INWDTH + NDIM
        INWORK = INMESH + NDIM
        IF ( MINCLS .LT. 0 ) SBRGNS = WORK(LENWRK)
        CALL ADBASE(NDIM, MINCLS, MAXCLS, FUNCTN, ABSREQ, RELREQ, 
     &       ABSEST, FINEST, SBRGNS, MXRGNS, RULCLS, LENRUL, 
     &       WORK(INERRS), WORK(INVALS), WORK(INPTRS), WORK(INLWRS), 
     &       WORK(INUPRS), WORK(INMSHS), WORK(INWGTS), WORK(INPNTS), 
     &       WORK(INLOWR), WORK(INUPPR), WORK(INWDTH), WORK(INMESH), 
     &       WORK(INWORK), INFORM)
        WORK(LENWRK) = SBRGNS
       ELSE
        INFORM = 2
        MINCLS = RULCLS
      ENDIF
      END
      SUBROUTINE BSINIT(NDIM, W, LENRUL, G)
*
*     For initializing basic rule weights and symmetric sum parameters.
*
      INTEGER NDIM, LENRUL, RULPTS(6), I, J, NUMNUL, SDIM
      PARAMETER ( NUMNUL = 4, SDIM = 12 )
      DOUBLE PRECISION W(LENRUL,4), G(NDIM,LENRUL) 
      DOUBLE PRECISION LAM1, LAM2, LAM3, LAMP, RULCON
*
*     The following code determines rule parameters and weights for a
*      degree 7 rule (W(1,1),...,W(5,1)), two degree 5 comparison rules
*      (W(1,2),...,W(5,2) and W(1,3),...,W(5,3)) and a degree 3 
*      comparison rule (W(1,4),...W(5,4)).
*
*       If NDIM = 1, then LENRUL = 5 and total points = 9.
*       If NDIM < SDIM, then LENRUL = 6 and
*                      total points = 1+2*NDIM*(NDIM+2)+2**NDIM.
*       If NDIM > = SDIM, then LENRUL = 6 and
*                      total points = 1+2*NDIM*(1+2*NDIM).
*
      DO I = 1,LENRUL
         DO J = 1,NDIM
            G(J,I) = 0
         END DO
         DO J = 1,NUMNUL
            W(I,J) = 0
         END DO
      END DO
      RULPTS(5) = 2*NDIM*(NDIM-1)
      RULPTS(4) = 2*NDIM
      RULPTS(3) = 2*NDIM
      RULPTS(2) = 2*NDIM
      RULPTS(1) = 1
      LAMP = 0.85
      LAM3 = 0.4707
      LAM2 = 4/(15 - 5/LAM3)
      W(5,1) = ( 3 - 5*LAM3 )/( 180*(LAM2-LAM3)*LAM2**2 )
      IF ( NDIM .LT. SDIM ) THEN 
         LAM1 = 8*LAM3*(31*LAM3-15)/( (3*LAM3-1)*(5*LAM3-3)*35 )
         W(LENRUL,1) = 1/(3*LAM3)**3/2**NDIM
      ELSE
         LAM1 = ( LAM3*(15 - 21*LAM2) + 35*(NDIM-1)*(LAM2-LAM3)/9 )
     &       /  ( LAM3*(21 - 35*LAM2) + 35*(NDIM-1)*(LAM2/LAM3-1)/9 )
         W(6,1) = 1/(4*(3*LAM3)**3)
      ENDIF
      W(3,1) = ( 15 - 21*(LAM3+LAM1) + 35*LAM3*LAM1 )
     &     /( 210*LAM2*(LAM2-LAM3)*(LAM2-LAM1) ) - 2*(NDIM-1)*W(5,1)
      W(2,1) = ( 15 - 21*(LAM3+LAM2) + 35*LAM3*LAM2 )
     &     /( 210*LAM1*(LAM1-LAM3)*(LAM1-LAM2) )
      IF ( NDIM .LT. SDIM ) THEN
         RULPTS(LENRUL) = 2**NDIM
         LAM3 = SQRT(LAM3)
         DO I = 1,NDIM
            G(I,LENRUL) = LAM3
         END DO
      ELSE
         W(6,1) = 1/(4*(3*LAM3)**3)
         RULPTS(6) = 2*NDIM*(NDIM-1)
         LAM3 = SQRT(LAM3)
         DO I = 1,2
            G(I,6) = LAM3
         END DO
      ENDIF
      IF ( NDIM .GT. 1 ) THEN
         W(5,2) = 1/(6*LAM2)**2 
         W(5,3) = 1/(6*LAM2)**2 
      ENDIF
      W(3,2) = ( 3 - 5*LAM1 )/( 30*LAM2*(LAM2-LAM1) ) 
     &     - 2*(NDIM-1)*W(5,2) 
      W(2,2) = ( 3 - 5*LAM2 )/( 30*LAM1*(LAM1-LAM2) )
      W(4,3) = ( 3 - 5*LAM2 )/( 30*LAMP*(LAMP-LAM2) )
      W(3,3) = ( 3 - 5*LAMP )/( 30*LAM2*(LAM2-LAMP) ) 
     &     - 2*(NDIM-1)*W(5,3)
      W(2,4) = 1/(6*LAM1)
      LAMP = SQRT(LAMP)
      LAM2 = SQRT(LAM2)
      LAM1 = SQRT(LAM1)
      G(1,2) = LAM1
      G(1,3) = LAM2
      G(1,4) = LAMP
      IF ( NDIM .GT. 1 ) THEN
         G(1,5) = LAM2
         G(2,5) = LAM2
      ENDIF
      DO J = 1, NUMNUL
         W(1,J) = 1
         DO I = 2,LENRUL
            W(1,J) = W(1,J) - RULPTS(I)*W(I,J)
         END DO
      END DO
      RULCON = 2
      CALL RULNRM( LENRUL, NUMNUL, RULPTS, W, RULCON )
      END
      SUBROUTINE RULNRM( LENRUL, NUMNUL, RULPTS, W, RULCON )
      INTEGER LENRUL, NUMNUL, I, J, K, RULPTS(*)
      DOUBLE PRECISION ALPHA, NORMCF, NORMNL, W(LENRUL, *), RULCON
*
*     Compute orthonormalized null rules.
*
      NORMCF = 0
      DO I = 1,LENRUL
         NORMCF = NORMCF + RULPTS(I)*W(I,1)*W(I,1)
      END DO
      DO K = 2,NUMNUL
         DO I = 1,LENRUL
            W(I,K) = W(I,K) - W(I,1)
         END DO
         DO J = 2,K-1
            ALPHA = 0
            DO I = 1,LENRUL
               ALPHA = ALPHA + RULPTS(I)*W(I,J)*W(I,K)
            END DO
            ALPHA = -ALPHA/NORMCF
            DO I = 1,LENRUL
               W(I,K) = W(I,K) + ALPHA*W(I,J)
            END DO
         END DO
         NORMNL = 0
         DO I = 1,LENRUL
            NORMNL = NORMNL + RULPTS(I)*W(I,K)*W(I,K)
         END DO
         ALPHA = SQRT(NORMCF/NORMNL)
         DO I = 1,LENRUL
            W(I,K) = ALPHA*W(I,K)
         END DO
      END DO
      DO J = 2, NUMNUL
         DO I = 1,LENRUL
            W(I,J) = W(I,J)/RULCON
         END DO
      END DO
      END
      SUBROUTINE ADBASE(NDIM, MINCLS, MAXCLS, FUNCTN, ABSREQ, RELREQ,
     &     ABSEST, FINEST, SBRGNS, MXRGNS, RULCLS, LENRUL,
     &     ERRORS, VALUES, PONTRS, LOWERS, 
     &     UPPERS, MESHES, WEGHTS, POINTS, 
     &     LOWER, UPPER, WIDTH, MESH, WORK, INFORM)
*
*        Main adaptive integration subroutine
*
      EXTERNAL FUNCTN
      INTEGER I, J, NDIM, MINCLS, MAXCLS, SBRGNS, MXRGNS, 
     &     RULCLS, LENRUL, INFORM, NWRGNS 
      DOUBLE PRECISION FUNCTN, ABSREQ, RELREQ, ABSEST, FINEST,   
     &     ERRORS(*), VALUES(*), PONTRS(*),
     &     LOWERS(NDIM,*), UPPERS(NDIM,*),
     &     MESHES(NDIM,*),WEGHTS(*), POINTS(*),
     &     LOWER(*), UPPER(*), WIDTH(*), MESH(*), WORK(*) 
      INTEGER DIVAXN, TOP, RGNCLS, FUNCLS, DIFCLS
      
*
*     Initialization of subroutine
*
      INFORM = 2
      FUNCLS = 0
      CALL BSINIT(NDIM, WEGHTS, LENRUL, POINTS)
      IF ( MINCLS .GE. 0) THEN
*
*       When MINCLS >= 0 determine initial subdivision of the
*       integration region and apply basic rule to each subregion.
*
         SBRGNS = 0
         DO I = 1,NDIM
            LOWER(I) = 0
            MESH(I) = 1
            WIDTH(I) = 1/(2*MESH(I))
            UPPER(I) = 1
         END DO
         DIVAXN = 0
         RGNCLS = RULCLS
         NWRGNS = 1
 10      CALL DIFFER(NDIM, LOWER, UPPER, WIDTH, WORK, WORK(NDIM+1),  
     &        FUNCTN, DIVAXN, DIFCLS)
         FUNCLS = FUNCLS + DIFCLS
         IF ( FUNCLS + RGNCLS*(MESH(DIVAXN)+1)/MESH(DIVAXN)
     &        .LE. MINCLS ) THEN
            RGNCLS = RGNCLS*(MESH(DIVAXN)+1)/MESH(DIVAXN)
            NWRGNS = NWRGNS*(MESH(DIVAXN)+1)/MESH(DIVAXN)
            MESH(DIVAXN) = MESH(DIVAXN) + 1
            WIDTH(DIVAXN) = 1/( 2*MESH(DIVAXN) )
            GO TO 10
         ENDIF
         IF ( NWRGNS .LE. MXRGNS ) THEN
            DO I = 1,NDIM
               UPPER(I) = LOWER(I) + 2*WIDTH(I)
               MESH(I) = 1
            END DO
         ENDIF
*     
*     Apply basic rule to subregions and store results in heap.
*     
 20      SBRGNS = SBRGNS + 1
         CALL BASRUL(NDIM, LOWER, UPPER, WIDTH, FUNCTN, 
     &        WEGHTS, LENRUL, POINTS, WORK, WORK(NDIM+1), 
     &        ERRORS(SBRGNS),VALUES(SBRGNS))
         CALL TRESTR(SBRGNS, SBRGNS, PONTRS, ERRORS)
         DO I = 1,NDIM
            LOWERS(I,SBRGNS) = LOWER(I)
            UPPERS(I,SBRGNS) = UPPER(I)
            MESHES(I,SBRGNS) = MESH(I)
         END DO
         DO I = 1,NDIM
            LOWER(I) = UPPER(I)
            UPPER(I) = LOWER(I) + 2*WIDTH(I)
            IF ( LOWER(I)+WIDTH(I) .LT. 1 )  GO TO 20
            LOWER(I) = 0
            UPPER(I) = LOWER(I) + 2*WIDTH(I)
         END DO
         FUNCLS = FUNCLS + SBRGNS*RULCLS
      ENDIF
*     
*     Check for termination
*
 30   FINEST = 0
      ABSEST = 0
      DO I = 1, SBRGNS
         FINEST = FINEST + VALUES(I)
         ABSEST = ABSEST + ERRORS(I)
      END DO
      IF ( ABSEST .GT. MAX( ABSREQ, RELREQ*ABS(FINEST) )
     &     .OR. FUNCLS .LT. MINCLS ) THEN  
*     
*     Prepare to apply basic rule in (parts of) subregion with
*     largest error.
*     
         TOP = PONTRS(1)
         RGNCLS = RULCLS
         DO I = 1,NDIM
            LOWER(I) = LOWERS(I,TOP)
            UPPER(I) = UPPERS(I,TOP)
            MESH(I) = MESHES(I,TOP)
            WIDTH(I) = (UPPER(I)-LOWER(I))/(2*MESH(I))
            RGNCLS = RGNCLS*MESH(I)
         END DO
         CALL DIFFER(NDIM, LOWER, UPPER, WIDTH, WORK, WORK(NDIM+1),  
     &        FUNCTN, DIVAXN, DIFCLS)
         FUNCLS = FUNCLS + DIFCLS
         RGNCLS = RGNCLS*(MESH(DIVAXN)+1)/MESH(DIVAXN)
         IF ( FUNCLS + RGNCLS .LE. MAXCLS ) THEN
            IF ( SBRGNS + 1 .LE. MXRGNS ) THEN
*     
*     Prepare to subdivide into two pieces.
*    
               NWRGNS = 1
               WIDTH(DIVAXN) = WIDTH(DIVAXN)/2
            ELSE
               NWRGNS = 0
               WIDTH(DIVAXN) = WIDTH(DIVAXN)
     &                        *MESH(DIVAXN)/( MESH(DIVAXN) + 1 )
               MESHES(DIVAXN,TOP) = MESH(DIVAXN) + 1 
            ENDIF
            IF ( NWRGNS .GT. 0 ) THEN
*     
*     Only allow local subdivision when space is available.
*
               DO J = SBRGNS+1,SBRGNS+NWRGNS
                  DO I = 1,NDIM
                     LOWERS(I,J) = LOWER(I)
                     UPPERS(I,J) = UPPER(I)
                     MESHES(I,J) = MESH(I)
                  END DO
               END DO
               UPPERS(DIVAXN,TOP) = LOWER(DIVAXN) + 2*WIDTH(DIVAXN)
               LOWERS(DIVAXN,SBRGNS+1) = UPPERS(DIVAXN,TOP)
            ENDIF
            FUNCLS = FUNCLS + RGNCLS
            CALL BASRUL(NDIM, LOWERS(1,TOP), UPPERS(1,TOP), WIDTH, 
     &           FUNCTN, WEGHTS, LENRUL, POINTS, WORK, WORK(NDIM+1), 
     &           ERRORS(TOP), VALUES(TOP))
            CALL TRESTR(TOP, SBRGNS, PONTRS, ERRORS)
            DO I = SBRGNS+1, SBRGNS+NWRGNS
*     
*     Apply basic rule and store results in heap.
*     
               CALL BASRUL(NDIM, LOWERS(1,I), UPPERS(1,I), WIDTH,
     &              FUNCTN, WEGHTS, LENRUL, POINTS, WORK, WORK(NDIM+1),  
     &              ERRORS(I), VALUES(I))
               CALL TRESTR(I, I, PONTRS, ERRORS)
            END DO
            SBRGNS = SBRGNS + NWRGNS
            GO TO 30
         ELSE
            INFORM = 1
         ENDIF
      ELSE
         INFORM = 0
      ENDIF
      MINCLS = FUNCLS
      END
      SUBROUTINE BASRUL( NDIM, A, B, WIDTH, FUNCTN, W, LENRUL, G,
     &     CENTER, Z, RGNERT, BASEST )
*
*     For application of basic integration rule
*
      EXTERNAL FUNCTN
      INTEGER I, LENRUL, NDIM
      DOUBLE PRECISION 
     &     A(NDIM), B(NDIM), WIDTH(NDIM), FUNCTN, W(LENRUL,4), 
     &     G(NDIM,LENRUL), CENTER(NDIM), Z(NDIM), RGNERT, BASEST
      DOUBLE PRECISION 
     &     FULSUM, FSYMSM, RGNCMP, RGNVAL, RGNVOL, RGNCPT, RGNERR
*
*     Compute Volume and Center of Subregion
*
      RGNVOL = 1
      DO I = 1,NDIM
         RGNVOL = 2*RGNVOL*WIDTH(I)
         CENTER(I) = A(I) + WIDTH(I)
      END DO
      BASEST = 0
      RGNERT = 0
*
*     Compute basic rule and error
*
 10   RGNVAL = 0
      RGNERR = 0
      RGNCMP = 0
      RGNCPT = 0
      DO I = 1,LENRUL
         FSYMSM = FULSUM(NDIM, CENTER, WIDTH, Z, G(1,I), FUNCTN)
*     Basic Rule
         RGNVAL = RGNVAL + W(I,1)*FSYMSM
*     First comparison rule
         RGNERR = RGNERR + W(I,2)*FSYMSM
*     Second comparison rule
         RGNCMP = RGNCMP + W(I,3)*FSYMSM
*     Third Comparison rule
         RGNCPT = RGNCPT + W(I,4)*FSYMSM
      END DO
*
*     Error estimation
*
      RGNERR = SQRT(RGNCMP**2 + RGNERR**2)
      RGNCMP = SQRT(RGNCPT**2 + RGNCMP**2)
      IF ( 4*RGNERR .LT. RGNCMP ) RGNERR = RGNERR/2
      IF ( 2*RGNERR .GT. RGNCMP ) RGNERR = MAX( RGNERR, RGNCMP )
      RGNERT = RGNERT +  RGNVOL*RGNERR
      BASEST = BASEST +  RGNVOL*RGNVAL
*
*     When subregion has more than one piece, determine next piece and
*      loop back to apply basic rule.
*
      DO I = 1,NDIM
         CENTER(I) = CENTER(I) + 2*WIDTH(I)
         IF ( CENTER(I) .LT. B(I) ) GO TO 10
         CENTER(I) = A(I) + WIDTH(I)
      END DO
      END
      DOUBLE PRECISION FUNCTION FULSUM(S, CENTER, HWIDTH, X, G, F)
*
****  To compute fully symmetric basic rule sum
*
      EXTERNAL F
      INTEGER S, IXCHNG, LXCHNG, I, L
      DOUBLE PRECISION CENTER(S), HWIDTH(S), X(S), G(S), F
      DOUBLE PRECISION INTSUM, GL, GI
      FULSUM = 0
*
*     Compute centrally symmetric sum for permutation of G
*
 10   INTSUM = 0
      DO I = 1,S
         X(I) = CENTER(I) + G(I)*HWIDTH(I)
      END DO
 20   INTSUM = INTSUM + F(S,X)
      DO I = 1,S
         G(I) = -G(I)
         X(I) = CENTER(I) + G(I)*HWIDTH(I)
         IF ( G(I) .LT. 0 ) GO TO 20
      END DO
      FULSUM = FULSUM + INTSUM
*     
*     Find next distinct permuation of G and loop back for next sum
*     
      DO I = 2,S
         IF ( G(I-1) .GT. G(I) ) THEN
            GI = G(I)
            IXCHNG = I - 1
            DO L = 1,(I-1)/2
               GL = G(L)
               G(L) = G(I-L)
               G(I-L) = GL
               IF (  GL  .LE. GI ) IXCHNG = IXCHNG - 1
               IF ( G(L) .GT. GI ) LXCHNG = L
            END DO
            IF ( G(IXCHNG) .LE. GI ) IXCHNG = LXCHNG
            G(I) = G(IXCHNG)
            G(IXCHNG) = GI
            GO TO 10
         ENDIF
      END DO
*     
*     End loop for permutations of G and associated sums
*     
*     Restore original order to G's
*     
      DO I = 1,S/2
         GI = G(I)
         G(I) = G(S+1-I)
         G(S+1-I) = GI 
      END DO
      END
      SUBROUTINE DIFFER(NDIM, A, B, WIDTH, Z, DIF, FUNCTN, 
     &     DIVAXN, DIFCLS)
*
*     Compute fourth differences and subdivision axes
*
      EXTERNAL FUNCTN
      INTEGER I, NDIM, DIVAXN, DIFCLS
      DOUBLE PRECISION 
     &     A(NDIM), B(NDIM), WIDTH(NDIM), Z(NDIM), DIF(NDIM), FUNCTN
      DOUBLE PRECISION FRTHDF, FUNCEN, WIDTHI
      DIFCLS = 0
      DIVAXN = MOD( DIVAXN, NDIM ) + 1
      IF ( NDIM .GT. 1 ) THEN
         DO I = 1,NDIM 
            DIF(I) = 0
            Z(I) = A(I) + WIDTH(I)
         END DO
 10      FUNCEN = FUNCTN(NDIM, Z)
         DO I = 1,NDIM
            WIDTHI = WIDTH(I)/5
            FRTHDF = 6*FUNCEN
            Z(I) = Z(I) - 4*WIDTHI
            FRTHDF = FRTHDF + FUNCTN(NDIM,Z)
            Z(I) = Z(I) + 2*WIDTHI
            FRTHDF = FRTHDF - 4*FUNCTN(NDIM,Z)
            Z(I) = Z(I) + 4*WIDTHI
            FRTHDF = FRTHDF - 4*FUNCTN(NDIM,Z)
            Z(I) = Z(I) + 2*WIDTHI
            FRTHDF = FRTHDF + FUNCTN(NDIM,Z)
*     Do not include differences below roundoff
            IF ( FUNCEN + FRTHDF/8 .NE. FUNCEN ) 
     &           DIF(I) = DIF(I) + ABS(FRTHDF)*WIDTH(I)
            Z(I) = Z(I) - 4*WIDTHI
         END DO
         DIFCLS = DIFCLS + 4*NDIM + 1
         DO I = 1,NDIM
            Z(I) = Z(I) + 2*WIDTH(I)
            IF ( Z(I) .LT. B(I) ) GO TO 10
            Z(I) = A(I) + WIDTH(I)
         END DO
         DO I = 1,NDIM
            IF ( DIF(DIVAXN) .LT. DIF(I) ) DIVAXN = I
         END DO
      ENDIF
      END
      SUBROUTINE TRESTR(POINTR, SBRGNS, PONTRS, RGNERS)
****BEGIN PROLOGUE TRESTR
****PURPOSE TRESTR maintains a heap for subregions.
****DESCRIPTION TRESTR maintains a heap for subregions.
*            The subregions are ordered according to the size of the
*            greatest error estimates of each subregion (RGNERS).
*
*   PARAMETERS
*
*     POINTR Integer.
*            The index for the subregion to be inserted in the heap.
*     SBRGNS Integer.
*            Number of subregions in the heap.
*     PONTRS Real array of dimension SBRGNS.
*            Used to store the indices for the greatest estimated errors
*            for each subregion.
*     RGNERS Real array of dimension SBRGNS.
*            Used to store the greatest estimated errors for each 
*            subregion.
*
****ROUTINES CALLED NONE
****END PROLOGUE TRESTR
*
*   Global variables.
*
      INTEGER POINTR, SBRGNS
      DOUBLE PRECISION PONTRS(*), RGNERS(*)
*
*   Local variables.
*
*   RGNERR Intermediate storage for the greatest error of a subregion.
*   SUBRGN Position of child/parent subregion in the heap.
*   SUBTMP Position of parent/child subregion in the heap.
*
      INTEGER SUBRGN, SUBTMP
      DOUBLE PRECISION RGNERR
*
****FIRST PROCESSING STATEMENT TRESTR
*     
      RGNERR = RGNERS(POINTR)
      IF ( POINTR .EQ. PONTRS(1)) THEN
*
*        Move the new subregion inserted at the top of the heap 
*        to its correct position in the heap.
*
         SUBRGN = 1
 10      SUBTMP = 2*SUBRGN
         IF ( SUBTMP .LE. SBRGNS ) THEN
            IF ( SUBTMP .NE. SBRGNS ) THEN
*     
*              Find maximum of left and right child.
*
               IF ( RGNERS(PONTRS(SUBTMP)) .LT. 
     +              RGNERS(PONTRS(SUBTMP+1)) ) SUBTMP = SUBTMP + 1
            ENDIF
*
*           Compare maximum child with parent.
*           If parent is maximum, then done.
*
            IF ( RGNERR .LT. RGNERS(PONTRS(SUBTMP)) ) THEN
*     
*              Move the pointer at position subtmp up the heap.
*     
               PONTRS(SUBRGN) = PONTRS(SUBTMP)
               SUBRGN = SUBTMP
               GO TO 10
            ENDIF
         ENDIF
      ELSE
*
*        Insert new subregion in the heap.
*
         SUBRGN = SBRGNS
 20      SUBTMP = SUBRGN/2
         IF ( SUBTMP .GE. 1 ) THEN
*
*           Compare child with parent. If parent is maximum, then done.
*     
            IF ( RGNERR .GT. RGNERS(PONTRS(SUBTMP)) ) THEN
*     
*              Move the pointer at position subtmp down the heap.
*
               PONTRS(SUBRGN) = PONTRS(SUBTMP)
               SUBRGN = SUBTMP
               GO TO 20
            ENDIF
         ENDIF
      ENDIF
      PONTRS(SUBRGN) = POINTR
*
****END TRESTR
*
      END
*
      SUBROUTINE SPHMVT( N, NU, LOWER, UPPER, INFIN, CORREL, MAXPTS, 
     *                   ABSEPS, RELEPS, ERROR, VALUE, INFORM )
*
*     A subroutine for computing multivariate t probabilities.
*     This subroutine uses a modified version of the Mont-Carlo 
*     algorithm for multivariatie Normal probabilities in the paper
*       "Three Digit Accurate Multiple Normal Probabilities", 
*          pp. 369-380, Numer. Math. 35(1980), by I. Deak
*
*
*  Parameters
*
*     N      INTEGER, the number of variables.
*     NU     INTEGER, the number of degrees of freedom.
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
*     ERROR  REAL, estimated absolute error, with 99% confidence level.
*     VALUE  REAL, estimated value for the integral
*     INFORM INTEGER, termination status parameter:
*            if INFORM = 0, normal completion with ERROR < EPS;
*            if INFORM = 1, completion with ERROR > EPS and MAXPTS 
*                           function vaules used; increase MAXPTS to 
*                           decrease ERROR;
*            if INFORM = 2, N > 50.
*
      EXTERNAL SPMVTI
      INTEGER N, NU, INFIS, INFIN(*), MAXPTS, MPT, INFORM, NS, IVLS
      DOUBLE PRECISION CORREL(*), LOWER(*), UPPER(*), 
     *     ABSEPS, RELEPS, ERROR, VALUE, D, E, EPS, SPMVTI
      IF ( N .GT. 50 ) THEN
         INFORM = 2
         VALUE = 0
         ERROR = 1
         RETURN
      ENDIF
      INFORM = SPMVTI( N,NU, CORREL, LOWER,UPPER,INFIN, INFIS, D,E, NS )
      IF ( N-INFIS .EQ. 0 ) THEN
         VALUE = 1
         ERROR = 0
      ELSE IF ( N-INFIS .EQ. 1 ) THEN
         VALUE = E - D
         ERROR = 2E-16
      ELSE
*
*        Call the Monte-Carlo integration subroutine
*
         MPT = 25 + NS/N**3
         CALL TCRUDE( N-INFIS, MPT, ERROR, VALUE, 0 )
         IVLS = MPT*NS
 10      EPS = MAX( ABSEPS, RELEPS*ABS(VALUE) )
         IF ( ERROR .GT. EPS .AND. IVLS .LT. MAXPTS ) THEN
            MPT = MAX( MIN( INT( MPT*( ERROR/EPS )**2 ),
     *                      ( MAXPTS - IVLS )/NS ), 10 )
            CALL TCRUDE( N-INFIS, MPT, ERROR, VALUE, 1 )
            IVLS = IVLS + MPT*NS
            GO TO 10
         ENDIF
         IF ( ERROR. GT. EPS .AND. IVLS .GE. MAXPTS ) INFORM = 1
      ENDIF
      END
      DOUBLE PRECISION FUNCTION SPMVT(N)
*     
*     Integrand subroutine
*
      DOUBLE PRECISION LOWER(*), UPPER(*), CORREL(*), D, E, ZERO
      INTEGER N, INFIN(*), INFIS
      INTEGER NL, IJ, I, II, J, K, NS, NSO, ND, NU, NUIN
      PARAMETER ( NL = 50, ND = 2, ZERO = 0 )
      DOUBLE PRECISION A(NL), B(NL), U(NL,NL), Y(NL), COV(NL*(NL+1)/2)
      INTEGER INFI(NL), IS(NL), IC(NL)
      DOUBLE PRECISION RS, TMP, BT, RNOR, SPHLMT, SPMVTI
      SAVE NU, A, B, INFI, U
*
*    First generate U = COV*(random orthogonal matrix)
*
      DO K = N-1, 1, -1
         TMP = 0
         DO J = K, N
            Y(J) = RNOR()
            TMP = TMP + Y(J)**2
         END DO
         TMP = -SQRT(TMP)
         BT = 1/( TMP*( Y(K) + TMP ) )
         Y(K) = Y(K) + TMP
         DO I = 1, N
            TMP = 0
            DO J = K, N
               TMP = TMP + U(I,J)*Y(J)
            END DO
            TMP = BT*TMP
            DO J = K, N
               U(I,J) = U(I,J) - TMP*Y(J)
            END DO
         END DO
      END DO
*
*     Compute integrand average
*
      RS = SQRT( DBLE(ND) )
      DO I = 1,ND
         IC(I) = I
      END DO
      IC(ND+1) = N+1
      SPMVT = 0
      NS = 0
 10   DO I = 1,ND
         IS(I) = -1
      END DO
 20   DO I = 1, N
         TMP = 0
         DO J = 1,ND
            TMP = TMP + IS(J)*U( I, IC(J) )
         END DO
         Y(I) = TMP/RS
      END DO
      NS = NS + 1
      SPMVT = SPMVT + ( SPHLMT( N, NU, A, B, INFI, Y ) - SPMVT )/NS
      DO I = 1, ND
         IS(I) = IS(I) + 2
         IF ( IS(I) .LT. 2 ) GO TO 20
         IS(I) = -1
      END DO
      DO I = 1, ND
         IC(I) = IC(I) + 1
         IF ( IC(I) .LT. IC(I+1)  ) GO TO 10
         IC(I) = I
      END DO
      SPMVT = SPMVT/2
      RETURN
*
      ENTRY SPMVTI( N,NUIN, CORREL, LOWER,UPPER,INFIN, INFIS, D,E, NSO )
      SPMVTI = 0
      NU = NUIN
*
*     Initialisation
*
      II = 0
      IJ = 0
      INFIS = 0
      DO I = 1, N
         INFI(I) = INFIN(I) 
         IF ( INFI(I) .LT. 0 ) THEN
            INFIS = INFIS + 1
         ELSE 
            A(I) = 0
            B(I) = 0
            IF ( INFI(I) .NE. 0 ) A(I) = LOWER(I)
            IF ( INFI(I) .NE. 1 ) B(I) = UPPER(I)
         ENDIF
         DO J = 1, I-1
            II = II + 1
            IJ = IJ + 1
            COV(IJ) = CORREL(II)
         END DO
         IJ = IJ + 1
         COV(IJ) = 1
      END DO
      NSO = 1
      DO I = 1,ND
         NSO = 2*NSO*( N - INFIS - I + 1 )/I
      END DO
*
*     First move any doubly infinite limits to innermost positions
*
      IF ( INFIS .LT. N ) THEN
         DO I = N, N-INFIS+1, -1
            IF ( INFI(I) .GE. 0 ) THEN 
               DO J = 1,I-1
                  IF ( INFI(J) .LT. 0 ) THEN
                     CALL RCSWAP( J, I, A, B, INFI, N, COV )
                     GO TO 30
                  ENDIF
               END DO
            ENDIF
 30      END DO
      ENDIF
      II = 0
      DO I = 1, N-INFIS
         DO J = 1, I
            U(J,I) = 0
            II = II + 1
            U(I,J) = COV(II)
         END DO
      END DO
*
*     Determine Cholesky decomposition
*
      DO J = 1, N-INFIS
         DO I = J, N-INFIS
            TMP = U(I,J)
            DO K = 1, J-1
               TMP = TMP - U(I,K)*U(J,K)
            END DO
            IF ( I .EQ. J ) THEN
               U(J,J) = SQRT( MAX( TMP, ZERO ) )
            ELSE IF ( U(I,I) .GT. 0 ) THEN
               U(I,J) = TMP/U(J,J)
            ELSE
               U(I,J) = 0
            END IF
         END DO
      END DO
      DO I = 1, N-INFIS
         IF ( U(I,I) .GT. 0 ) THEN
            IF ( INFI(I) .NE. 0 ) A(I) = A(I)/U(I,I)
            IF ( INFI(I) .NE. 1 ) B(I) = B(I)/U(I,I)
            DO J = 1,I
               U(I,J) = U(I,J)/U(I,I)
            END DO
         ENDIF
      END DO
      CALL MVTLMS( NU, A(1), B(1), INFI(1), D, E )
      END
      DOUBLE PRECISION FUNCTION SPHLMT( N, NU, A, B, INFI, Y )
      DOUBLE PRECISION A(*), B(*), Y(*), CMN, CMX, SPHNCT
      INTEGER INFI(*), I, N, NU
      CMN = -10*N
      CMX =  10*N
      DO I = 1,N
         IF ( Y(I) .GT. 0 ) THEN
            IF ( INFI(I) .NE. 1 ) CMX = MIN( CMX, B(I)/Y(I) )
            IF ( INFI(I) .NE. 0 ) CMN = MAX( CMN, A(I)/Y(I) )
         ELSE
            IF ( INFI(I) .NE. 1 ) CMN = MAX( CMN, B(I)/Y(I) )
            IF ( INFI(I) .NE. 0 ) CMX = MIN( CMX, A(I)/Y(I) )
         ENDIF
      END DO
      IF ( CMN .LT. CMX ) THEN
         IF ( CMN .GE. 0 .AND. CMX .GE. 0 ) THEN
            SPHLMT = SPHNCT( N, NU,  CMX ) - SPHNCT( N, NU,  CMN )
         ELSEIF ( CMN .LT. 0 .AND. CMX .GE. 0 ) THEN
            SPHLMT = SPHNCT( N, NU, -CMN ) + SPHNCT( N, NU,  CMX )
         ELSE
            SPHLMT = SPHNCT( N, NU, -CMN ) - SPHNCT( N, NU, -CMX )
         ENDIF
      ELSE
         SPHLMT = 0
      ENDIF
      END
      SUBROUTINE TCRUDE( NDIM, MAXPTS, ABSEST, FINEST, IR )
*
*     Crude Monte-Carlo Algorithm for Deak method with
*      weighted results on restart
*
      INTEGER NDIM, MAXPTS, M, K, IR, NPTS
      DOUBLE PRECISION FINEST, ABSEST, X(100), SPMVT, UNI, 
     *     VARSQR, VAREST, VARPRD, FINDIF, FINVAL
      SAVE VAREST
      IF ( IR .LE. 0 ) THEN
         VAREST = 0
         FINEST = 0
      ENDIF
      FINVAL = 0
      VARSQR = 0
      DO M = 1, MAXPTS
         FINDIF = ( SPMVT(NDIM) - FINVAL )/M
         FINVAL = FINVAL + FINDIF
         VARSQR = ( M - 2 )*VARSQR/M + FINDIF**2 
      END DO
      VARPRD = VAREST*VARSQR
      FINEST = FINEST + ( FINVAL - FINEST )/(1 + VARPRD)
      IF ( VARSQR .GT. 0 ) VAREST = (1 + VARPRD)/VARSQR
      ABSEST = 3*SQRT( VARSQR/( 1 + VARPRD ) )
      END
      DOUBLE PRECISION FUNCTION SPHNCT( M, NU, R )
*     
*                   R  
*     SPHNCT =  K  I  ( 1 + t**2/NU )**(-(NU+M)/2 ) t**(M-1) dt, for M > 0.
*                M  0
*     
      INTEGER I, M, NU, NUOLD
      DOUBLE PRECISION R, RR, RT, PI, PF, STUDNT, TCON
      PARAMETER ( PI = 3.14159 26535 89793D0 )
      DATA NUOLD / 0 /
      SAVE NUOLD, TCON
      IF ( R .GT. 0 ) THEN
         IF ( M .LE. 1 ) THEN
            SPHNCT = 2*STUDNT( NU, R ) - 1
         ELSE IF ( M .EQ. 2 ) THEN
            SPHNCT = 1 - 1/SQRT( 1 + R*R/NU )**NU
         ELSE 
            RR = R*R/NU
            RT = RR/( 1 + RR )
            PF = 1
            DO I = M - 2, 2, -2
               PF = 1 + PF*RT*( NU + I - 2 )/I
            END DO
            PF = PF*SQRT( RT/RR )**NU
            IF ( MOD( M, 2 ) .EQ. 0 ) THEN
               SPHNCT = 1 - PF
            ELSE
               IF ( NU .NE. NUOLD ) THEN
                  NUOLD = NU
                  TCON = 1
                  IF ( MOD( NU, 2 ) .EQ. 0 ) THEN 
                     TCON = TCON/2
                  ELSE
                     TCON = TCON/PI
                  END IF
                  DO I = NU-2, 1, -2
                     TCON = ( I + 1 )*TCON/I
                  END DO
               END IF
               SPHNCT = 2*( STUDNT( NU, R ) - TCON*SQRT(RT)*PF ) - 1 
            ENDIF
         ENDIF
      ELSE
         SPHNCT = 0
      ENDIF
      END
      DOUBLE PRECISION FUNCTION RNOR()
*
*       RNOR generates normal random numbers with zero mean and
*       unit standard deviation, often denoted N(0,1).
*       Adapted from RNOR in "Numerical Methods and Software" by
*                D. Kahaner, C. Moler, S. Nash
*                Prentice Hall, 1988
*
      DOUBLE PRECISION AA, B, C, C1, C2, PC, X, Y, XN, V(65), S, VT, UNI
      PARAMETER ( AA = 12.37586, B = 0.4878992, C = 12.67706 ) 
      PARAMETER ( C1 = 0.9689279, C2 = 1.301198, PC = 0.1958303E-1 )  
      PARAMETER ( XN = 2.7769943 )
      INTEGER J
      SAVE V
      DATA V/ .3409450, .4573146, .5397793, .6062427, .6631691,
     & 0.7136975, 0.7596125, 0.8020356, 0.8417227, 0.8792102, 0.9148948,
     & 0.9490791, 0.9820005, 1.0138492, 1.0447810, 1.0749254, 1.1043917,
     & 1.1332738, 1.1616530, 1.1896010, 1.2171815, 1.2444516, 1.2714635,
     & 1.2982650, 1.3249008, 1.3514125, 1.3778399, 1.4042211, 1.4305929,
     & 1.4569915, 1.4834526, 1.5100121, 1.5367061, 1.5635712, 1.5906454,
     & 1.6179680, 1.6455802, 1.6735255, 1.7018503, 1.7306045, 1.7598422,
     & 1.7896223, 1.8200099, 1.8510770, 1.8829044, 1.9155830, 1.9492166,
     & 1.9839239, 2.0198430, 2.0571356, 2.0959930, 2.1366450, 2.1793713,
     & 2.2245175, 2.2725185, 2.3239338, 2.3795007, 2.4402218, 2.5075117,
     & 2.5834658, 2.6713916, 4*XN/
      Y = UNI()
      J = MOD( INT( UNI()*128 ), 64 ) + 1
*
*     Pick sign as Y+Y-1 is positive or negative
*
      VT = V(J+1)
      RNOR = (Y+Y-1)*VT
      IF ( ABS(RNOR) .GT. V(J) ) THEN
         X = ( ABS(RNOR)-V(J) )/( VT-V(J) )
         Y = UNI()
         S = X + Y
         IF ( S .GT. C1 ) THEN
            IF ( S .GT. C2 .OR. Y .GT. C - AA*EXP(-(B-B*X)**2/2) ) THEN
               RNOR = SIGN( B-B*X, RNOR )
            ELSE
               IF( EXP(-VT**2/2) + Y*PC/VT .GT. EXP(-RNOR**2/2) ) THEN
 10               X = 0.3601016*LOG( UNI() )
                  IF ( -2*LOG( UNI() ) .LE. X**2 ) GO TO 10
                  RNOR = SIGN( XN-X, RNOR )
               ENDIF
            ENDIF
         ENDIF
      ENDIF
      END
*
      DOUBLE PRECISION FUNCTION FNCMVT(N, W)
*     
*     Integrand subroutine
*
      INTEGER N, NUIN, INFIN(*), INFIS
      DOUBLE PRECISION W(*), LOWER(*), UPPER(*), CORREL(*), D, E
      INTEGER NL, IJ, I, J, NU
      PARAMETER ( NL = 20 )
      DOUBLE PRECISION COV((NL*(NL+1))/2), A(NL), B(NL), Y(NL)
      INTEGER INFI(NL)
      DOUBLE PRECISION PROD, D1, E1, DI, EI, SUM, STDINV, YD, UI, MVTNIT
      SAVE NU, D1, E1, A, B, INFI, COV
      DI = D1
      EI = E1
      PROD = EI - DI
      IJ = 1
      YD = 1
      DO I = 1, N
         UI = STDINV( NU+I-1, DI + W(I)*( EI - DI ) )
         Y(I) = UI/YD
         YD = YD/SQRT( 1 + ( UI - 1 )*( UI + 1 )/( NU + I ) )
         SUM = 0
         DO J = 1, I
            IJ = IJ + 1
            SUM = SUM + COV(IJ)*Y(J)
         END DO
         IJ = IJ + 1
         CALL MVTLMS( NU+I, ( A(I+1) - SUM )*YD, ( B(I+1) - SUM )*YD, 
     &                INFI(I+1), DI, EI ) 
         PROD = PROD*( EI - DI )
      END DO
      FNCMVT = PROD
      RETURN 
*
*     Entry point for intialization
*
      ENTRY MVTNIT( N, NUIN, CORREL, LOWER, UPPER, INFIN, INFIS, D, E ) 
      MVTNIT = 0
*
*     Initialization and computation of covariance matrix Cholesky factor
*
      CALL MVTSRT( N, NUIN, LOWER, UPPER, CORREL, INFIN, Y, INFIS, 
     &             A, B, INFI, COV, D, E )
      NU = NUIN
      D1 = D
      E1 = E
      END
      SUBROUTINE MVTLMS( NU, A, B, INFIN, LOWER, UPPER )
      DOUBLE PRECISION A, B, LOWER, UPPER, STUDNT
      INTEGER NU, INFIN
      LOWER = 0
      UPPER = 1
      IF ( INFIN .GE. 0 ) THEN
         IF ( INFIN .NE. 0 ) LOWER = STUDNT( NU, A )
         IF ( INFIN .NE. 1 ) UPPER = STUDNT( NU, B )
      ENDIF
      END
      SUBROUTINE MVTSRT( N, NU, LOWER, UPPER, CORREL, INFIN, Y, INFIS, 
     &                   A, B, INFI, COV, D, E )
*
*     Sort limits
*
      INTEGER N, NU, INFI(*), INFIN(*), INFIS
      DOUBLE PRECISION 
     &     A(*), B(*), COV(*), LOWER(*), UPPER(*), CORREL(*), Y(*), D, E
      INTEGER I, J, K, IJ, II, JMIN
      DOUBLE PRECISION SUMSQ, ZERO, TWO, PI, CVDIAG
      DOUBLE PRECISION AI, BI, SUM, YL, YU, YD
      DOUBLE PRECISION AMIN, BMIN, DMIN, EMIN, CON, CONODD, CONEVN
      PARAMETER ( ZERO = 0, TWO = 2, PI = 3.14159 26535 89793 23844 )
      IJ = 0
      II = 0
      INFIS = 0
      DO I = 1, N
         INFI(I) = INFIN(I)
         IF ( INFI(I) .LT. 0 ) THEN
            INFIS = INFIS + 1
         ELSE
            A(I) = 0
            B(I) = 0
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
      CONODD = 1/PI
      CONEVN = 1/TWO
      DO I = 1, NU - 1
         IF ( MOD(I,2) .EQ. 0 ) THEN
            IF ( I .GT. 2 ) CONEVN = CONEVN*(I-1)/(I-2)
         ELSE
            IF ( I .GT. 2 ) CONODD = CONODD*(I-1)/(I-2)
         ENDIF
      END DO
*
*     First move any doubly infinite limits to innermost positions
*
      IF ( INFIS .LT. N ) THEN
         DO I = N, N-INFIS+1, -1
            IF ( INFI(I) .GE. 0 ) THEN
               DO J = 1, I-1
                  IF ( INFI(J) .LT. 0 ) THEN
                     CALL RCSWAP( J, I, A, B, INFI, N, COV )
                     GOTO 10
                  ENDIF
               END DO
            ENDIF
 10      END DO
*
*     Sort remaining limits and determine Cholesky decomposition
*
         II = 0
         YD = 1
         DO I = 1, N-INFIS
*
*     Determine the integration limits for variable with minimum
*      expected probability and interchange that variable with Ith.
*
            EMIN = 1
            DMIN = 0
            JMIN = I
            CVDIAG = 0
            IJ = II
            DO J = I, N-INFIS
               SUM = 0
               SUMSQ = 0
               DO K = 1, I-1
                  SUM = SUM + COV(IJ+K)*Y(K)
                  SUMSQ = SUMSQ + COV(IJ+K)**2
               END DO
               IJ = IJ + J
               SUMSQ = SQRT( MAX( COV(IJ)-SUMSQ, ZERO ) )
               IF ( SUMSQ .GT. 0 ) THEN
                  AI = YD*( A(J) - SUM )/SUMSQ
                  BI = YD*( B(J) - SUM )/SUMSQ
                  CALL MVTLMS( NU+J-1, AI, BI, INFI(J), D, E )
                  IF ( EMIN - DMIN .GE. E - D ) THEN
                     JMIN = J
                     AMIN = AI
                     BMIN = BI
                     DMIN = D
                     EMIN = E
                     CVDIAG = SUMSQ
                  ENDIF
               ENDIF
            END DO
            IF ( JMIN .NE. I ) CALL RCSWAP( I, JMIN, A,B, INFI, N,COV )
*
*     Compute Ith column of Cholesky factor.
*
            IJ = II + I
            COV(IJ) = CVDIAG
            DO J = I+1, N-INFIS
               IF ( CVDIAG .GT. 0 ) THEN
                  SUM = COV(IJ+I)
                  DO K = 1, I-1
                     SUM = SUM - COV(II+K)*COV(IJ+K)
                  END DO
                  COV(IJ+I) = SUM/CVDIAG
               ELSE
                  COV(IJ+I) = 0
               ENDIF
               IJ = IJ + J
            END DO
*
*     Compute expected value for Ith integration variable and
*     scale Ith covariance matrix row and limits.
*
            IF ( MOD(NU+I-1,2) .EQ. 0 ) THEN
               IF ( NU+I-3 .GT. 0 ) CONEVN = CONEVN*(NU+I-2)/(NU+I-3)
               CON = CONEVN
            ELSE
               IF ( NU+I-3 .GT. 0 ) CONODD = CONODD*(NU+I-2)/(NU+I-3)
               CON = CONODD
            ENDIF
            IF ( CVDIAG .GT. 0 ) THEN
               YL = 0
               YU = 0
               IF ( INFI(I) .NE. 0 .AND. NU+I-2 .GT. 0 ) 
     &              YL = -CON*(NU+I-1)/(NU+I-2)
     &              /( 1 + AMIN**2/(NU+I-1) )**( (NU+I-2)/TWO )
               IF ( INFI(I) .NE. 1 .AND. NU+I-2 .GT. 0 ) 
     &              YU = -CON*(NU+I-1)/(NU+I-2)
     &              /( 1 + BMIN**2/(NU+I-1) )**( (NU+I-2)/TWO )
               Y(I) = ( YU - YL )/( EMIN - DMIN )/YD
               DO J = 1,I
                  II = II + 1
                  COV(II) = COV(II)/CVDIAG
               END DO
               IF ( INFI(I) .NE. 0 ) A(I) = A(I)/CVDIAG
               IF ( INFI(I) .NE. 1 ) B(I) = B(I)/CVDIAG
            ELSE
               Y(I) = 0
               II = II + I
            ENDIF
            YD = YD/SQRT( 1 + ( Y(I)*YD + 1 )*( Y(I)*YD - 1 )/(NU+I) )
         END DO
         CALL MVTLMS( NU, A(1), B(1), INFI(1), D, E)
      ENDIF
      END
      SUBROUTINE RCSWAP( P, Q, A, B, INFIN, N, C )
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
      END
      SUBROUTINE RANCOV( N, COV, IK )
*
*     Generates random correlation matrix in situ
*
      INTEGER N, I, II, J, IK
      DOUBLE PRECISION COV(*), SG, T, RHO, UNI, S
      IF ( IK .EQ. 0 ) THEN
         RHO = UNI()
         II = 0
         DO I = 1,N
            DO J = 1,I-1
               COV(II+J) = RHO
            END DO
            II = II + I
            COV(II) = 1
         END DO
      ELSE
         II = 0
         DO I = 1,N
            S = 0
            DO J = 1,I
               T = 2*UNI() - 1
               COV(II+J) = T
               S = S + T*T
            END DO
            SG = 1
            IF ( COV(II + I) .LT. 0 ) SG = -1
            S = SQRT(S)
            DO J = 1,I
               COV(II+J) = SG*COV(II+J)/S
            END DO
            II = II + I
         END DO
         CALL CHOLPD(N, COV)
      ENDIF
      END
      SUBROUTINE CHOLPD(N, CHOPRD)
*
*     Multiplies Choleski factors in situ
*
      INTEGER I, II, J, K, KK, N, NN
      DOUBLE PRECISION CHOPRD(*), S
      NN = (N*(N+1))/2
      KK = NN
      DO K = N,1,-1
         KK = KK - K
         II = NN
         DO I = N,K,-1
            II = II - I
            S = 0
            DO J = 1,K
               S = S + CHOPRD(II+J)*CHOPRD(KK+J)
            END DO
            CHOPRD(II+K) = S
         END DO
      END DO
      END
      SUBROUTINE CHOLSK(N, CHOFAC)
*
*     Computes Choleski factor in situ
*
      INTEGER I, II, J, JJ, K, N
      DOUBLE PRECISION CHOFAC(*), S, T, ZERO
      PARAMETER ( ZERO = 0 )
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
      END
      DOUBLE PRECISION FUNCTION STUDNT( NU, T )
*
*     Student t Distribution Function
*
*                       T
*         STUDNT = C   I  ( 1 + y*y/NU )**( -(NU+1)/2 ) dy
*                   NU -INF
*
      INTEGER NU, J
      DOUBLE PRECISION T, CSSTHE, SNTHE, POLYN, TT, TS, RN, PI, ZERO
      PARAMETER ( PI = 3.14159 26535 89793D0, ZERO = 0 )
      IF ( NU .EQ. 1 ) THEN
         STUDNT = ( 1 + 2*ATAN(T)/PI )/2
      ELSE IF ( NU .EQ. 2) THEN
         STUDNT = ( 1 + T/SQRT( 2 + T*T ))/2
      ELSE 
         TT = T*T
         CSSTHE = 1/( 1 + TT/NU )
         POLYN = 1
         DO J = NU-2, 2, -2
            POLYN = 1 + ( J - 1 )*CSSTHE*POLYN/J
         END DO
         IF ( MOD( NU, 2 ) .EQ. 1 ) THEN
            RN = NU
            TS = T/SQRT(RN)
            STUDNT = ( 1 + 2*( ATAN(TS) + TS*CSSTHE*POLYN )/PI )/2
         ELSE
            SNTHE = T/SQRT( NU + TT )
            STUDNT = ( 1 + SNTHE*POLYN )/2
         END IF
         STUDNT = MAX( ZERO, STUDNT )
      ENDIF
      END
      DOUBLE PRECISION FUNCTION STDINV( N, Z )
*
*     Inverse Student t Distribution Function
*
*                    STDINV
*           Z = C   I      (1 + y*y/N)**(-(N+1)/2) dy
*                N  -INF
*
*      Reference: G.W. Hill, Comm. ACM Algorithm 395
*                 Comm. ACM 13 (1970), pp. 619-620.
*
*      Conversions to double precision and other modifications by
*                 Alan Genz, 1993-4.
*
      INTEGER N
      DOUBLE PRECISION Z, P, PHINV, A, B, C, D, X, Y, PI, TWO
      DOUBLE PRECISION STUDNT, STDJAC
      PARAMETER ( PI = 3.14159 26535 89793D0, TWO = 2  )
      IF ( 0 .LT. Z .AND. Z .LT. 1 ) THEN
         IF ( N .EQ. 1 ) THEN
            STDINV = TAN( PI*( 2*Z - 1 )/2 )
         ELSE IF ( N .EQ. 2) THEN
            STDINV = ( 2*Z - 1 )/SQRT( 2*Z*( 1 - Z ) )
         ELSE 
            IF ( 2*Z .GE. 1 ) THEN 
               P = 2*( 1 - Z )
            ELSE
               P = 2*Z
            END IF
            A = 1/( N - 0.5 )
            B = 48/( A*A )
            C = ( ( 20700*A/B - 98 )*A - 16 )*A + 96.36
            D = ( ( 94.5/( B + C ) - 3 )/B + 1 )*SQRT( A*PI/2 )*N
            X = D*P
            Y = X**( TWO/N )
            IF ( Y .GT. A + 0.05 ) THEN
               X = PHINV( P/2 )
               Y = X*X
               IF ( N .LT. 5 ) C = C + 3*( N - 4.5 )*( 10*X + 6 )/100
               C = ( ( (D*X - 100)*X/20 - 7 )*X - 2 )*X + B + C
               Y = ( ( ( ( (4*Y+63)*Y/10+36 )*Y+94.5 )/C-Y-3 )/B + 1 )*X
               Y = A*Y*Y
               IF ( Y .GT. 0.002 ) THEN
                  Y = EXP(Y) - 1
               ELSE
                  Y = Y*( 1 + Y/2 )
               ENDIF
            ELSE
               Y = ( ( 1/( ( (N+6)/(N*Y) - 0.089*D - 0.822 )*(3*N+6) )
     &              + 0.5/(N+4) )*Y - 1 )*(N+1)/(N+2) + 1/Y
            END IF
            STDINV = SQRT(N*Y)
            IF ( 2*Z .LT. 1 ) STDINV = -STDINV
            IF ( ABS( STDINV ) .GT. 0 ) THEN
*
*     Use one third order correction to the single precision result
*
               X = STDINV
               D = Z - STUDNT(N,X)
               STDINV = X + 2*D/( 2/STDJAC(N,X) - D*(N+1)/(N/X+X) )
            END IF
         END IF
      ELSE
*
*     Use cutoff values for Z near 0 or 1.
*
         STDINV = SQRT( N/( 2D-16*SQRT( 2*PI*N ) )**( TWO/N ) )
         IF ( 2*Z .LT. 1 ) STDINV = -STDINV
      END IF
      END
      DOUBLE PRECISION FUNCTION STDJAC( NU, T )
*
*     Student t Distribution Transformation Jacobean
*
*          T            STDINV(NU,T)
*         I  f(y) dy = I   f(STDINV(NU,Z) STDJAC(NU,STDINV(NU,Z)) dZ
*         -INF          0
*
      INTEGER NU, J
      DOUBLE PRECISION CONST, NUOLD, PI, T, TT
      PARAMETER ( PI = 3.14159 26535 89793D0 )
      SAVE NUOLD, CONST
      DATA NUOLD/ 0D0 /
      IF ( NU .EQ. 1 ) THEN
         STDJAC = PI*( 1 + T*T )
      ELSE IF ( NU .EQ. 2 ) THEN 
         STDJAC = SQRT( 2 + T*T )**3
      ELSE 
         IF ( NU .NE. NUOLD ) THEN
            NUOLD = NU
            IF ( MOD( NU, 2 ) .EQ. 0 ) THEN
               CONST = SQRT(NUOLD)*2
            ELSE
               CONST = SQRT(NUOLD)*PI
            END IF
            DO J = NU-2, 1, -2
               CONST = J*CONST/(J+1)
            END DO
         END IF 
         TT = 1 + T*T/NU
         STDJAC = CONST*TT**( (NU+1)/2 ) 
         IF ( MOD( NU, 2 ) .EQ. 0 ) STDJAC = STDJAC*SQRT( TT )
      END IF
      END
      DOUBLE PRECISION FUNCTION PHI(Z)
*
*	Normal distribution probabilities accurate to 1.e-15.
*	Z = no. of standard deviations from the mean.
*
*       Based upon algorithm 5666 for the error function, from:
*       Hart, J.F. et al, 'Computer Approximations', Wiley 1968
*
*       Programmer: Alan Miller
*
*	Latest revision - 30 March 1986
*
      DOUBLE PRECISION P0, P1, P2, P3, P4, P5, P6, 
     &     Q0, Q1, Q2, Q3, Q4, Q5, Q6, Q7,
     &     Z, P, EXPNTL, CUTOFF, ROOTPI, ZABS
      PARAMETER(  P0 = 220.20 68679 12376 1D0,
     &	          P1 = 221.21 35961 69931 1D0, 
     &            P2 = 112.07 92914 97870 9D0,
     &	          P3 = 33.912 86607 83830 0D0,
     &            P4 = 6.3739 62203 53165 0D0,
     &	          P5 = .70038 30644 43688 1D0, 
     &            P6 = .035262 49659 98910 9D0 )
      PARAMETER(  Q0 = 440.41 37358 24752 2D0,
     &	          Q1 = 793.82 65125 19948 4D0, 
     &            Q2 = 637.33 36333 78831 1D0,
     &	          Q3 = 296.56 42487 79673 7D0, 
     &            Q4 = 86.780 73220 29460 8D0,
     &	          Q5 = 16.064 17757 92069 5D0, 
     &            Q6 = 1.7556 67163 18264 2D0,
     &	          Q7 = .088388 34764 83184 4D0 )
      PARAMETER(  ROOTPI = 2.5066 28274 63100 1D0 )
      PARAMETER(  CUTOFF = 7.0710 67811 86547 5D0 )
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
         EXPNTL = EXP(-ZABS**2/2)
*     
*     |Z| < CUTOFF = 10/SQRT(2)
*     
         IF ( ZABS .LT. CUTOFF ) THEN
            P = EXPNTL*((((((P6*ZABS + P5)*ZABS + P4)*ZABS + P3)*ZABS
     &          + P2)*ZABS + P1)*ZABS + P0)/(((((((Q7*ZABS + Q6)*ZABS
     &          + Q5)*ZABS + Q4)*ZABS + Q3)*ZABS + Q2)*ZABS + Q1)*ZABS
     &          + Q0)
*
*     |Z| >= CUTOFF.
*     
         ELSE
            P = EXPNTL/(ZABS + 1/(ZABS + 2/(ZABS + 3/(ZABS + 4/
     &           (ZABS + 0.65000 00000 00000 0D0)))))/ROOTPI
         END IF
      END IF
      IF (Z .GT. 0) P = 1 - P
      PHI = P
      END
      DOUBLE PRECISION FUNCTION PHINV(P)
*
*     ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3
*     
*     Produces the normal deviate Z corresponding to a given lower
*     tail area of P.
*     
*     The hash sums below are the sums of the mantissas of the
*     coefficients.   They are included for use in checking
*     transcription.
*     
      DOUBLE PRECISION SPLIT1, SPLIT2, CONST1, CONST2,
     &     A0, A1, A2, A3, A4, A5, A6, A7, B1, B2, B3, B4, B5, B6, B7,
     &     C0, C1, C2, C3, C4, C5, C6, C7, D1, D2, D3, D4, D5, D6, D7, 
     &     E0, E1, E2, E3, E4, E5, E6, E7, F1, F2, F3, F4, F5, F6, F7, 
     &     P, Q, R
      PARAMETER (SPLIT1 = 0.425D0, SPLIT2 = 5,
     &     CONST1 = 0.180625D0, CONST2 = 1.6D0)
*     
*     Coefficients for P close to 0.5
*     
      PARAMETER (  A0 = 3.38713 28727 96366 6080D0,
     &		   A1 = 1.33141 66789 17843 7745D+2,
     &		   A2 = 1.97159 09503 06551 4427D+3,
     &		   A3 = 1.37316 93765 50946 1125D+4,
     &		   A4 = 4.59219 53931 54987 1457D+4,
     &		   A5 = 6.72657 70927 00870 0853D+4,
     &		   A6 = 3.34305 75583 58812 8105D+4,
     &		   A7 = 2.50908 09287 30122 6727D+3,
     &		   B1 = 4.23133 30701 60091 1252D+1,
     &		   B2 = 6.87187 00749 20579 0830D+2,
     &		   B3 = 5.39419 60214 24751 1077D+3,
     &		   B4 = 2.12137 94301 58659 5867D+4,
     &		   B5 = 3.93078 95800 09271 0610D+4,
     &		   B6 = 2.87290 85735 72194 2674D+4,
     &		   B7 = 5.22649 52788 52854 5610D+3)
*     HASH SUM AB      55.88319 28806 14901 4439
*     
*     Coefficients for P not close to 0, 0.5 or 1.
*     
      PARAMETER (  C0 = 1.42343 71107 49683 57734D0,
     &             C1 = 4.63033 78461 56545 29590D0,
     &		   C2 = 5.76949 72214 60691 40550D0,
     &		   C3 = 3.64784 83247 63204 60504D0,
     &		   C4 = 1.27045 82524 52368 38258D0,
     &		   C5 = 2.41780 72517 74506 11770D-1,
     &		   C6 = 2.27238 44989 26918 45833D-2,
     &		   C7 = 7.74545 01427 83414 07640D-4,
     &		   D1 = 2.05319 16266 37758 82187D0,
     &		   D2 = 1.67638 48301 83803 84940D0,
     &		   D3 = 6.89767 33498 51000 04550D-1,
     &		   D4 = 1.48103 97642 74800 74590D-1,
     &		   D5 = 1.51986 66563 61645 71966D-2,
     &		   D6 = 5.47593 80849 95344 94600D-4,
     &		   D7 = 1.05075 00716 44416 84324D-9)
*	HASH SUM CD    49.33206 50330 16102 89036
*
*	Coefficients for P near 0 or 1.
*
      PARAMETER (  E0 = 6.65790 46435 01103 77720D0,
     &		   E1 = 5.46378 49111 64114 36990D0,
     &		   E2 = 1.78482 65399 17291 33580D0,
     &		   E3 = 2.96560 57182 85048 91230D-1,
     &		   E4 = 2.65321 89526 57612 30930D-2,
     &		   E5 = 1.24266 09473 88078 43860D-3,
     &		   E6 = 2.71155 55687 43487 57815D-5,
     &		   E7 = 2.01033 43992 92288 13265D-7,
     &		   F1 = 5.99832 20655 58879 37690D-1,
     &		   F2 = 1.36929 88092 27358 05310D-1,
     &		   F3 = 1.48753 61290 85061 48525D-2,
     &		   F4 = 7.86869 13114 56132 59100D-4,
     &		   F5 = 1.84631 83175 10054 68180D-5,
     &		   F6 = 1.42151 17583 16445 88870D-7,
     &		   F7 = 2.04426 31033 89939 78564D-15)
*     HASH SUM EF      47.52583 31754 92896 71629
*     
      Q = ( 2*P - 1 )/2
      IF ( ABS(Q) .LE. SPLIT1 ) THEN
         R = CONST1 - Q*Q
         PHINV = Q*(((((((A7*R + A6)*R + A5)*R + A4)*R + A3)
     &        *R + A2)*R + A1)*R + A0) /
     &        (((((((B7*R + B6)*R + B5)*R + B4)*R + B3)
     &        *R + B2)*R + B1)*R + 1)
      ELSE
         R = MIN( P, 1 - P )
         IF ( R .GT. 0 ) THEN
            R = SQRT(-LOG(R))
            IF ( R .LE. SPLIT2 ) THEN
               R = R - CONST2
               PHINV = (((((((C7*R + C6)*R + C5)*R + C4)*R + C3)
     &              *R + C2)*R + C1)*R + C0) /
     &              (((((((D7*R + D6)*R + D5)*R + D4)*R + D3)
     &              *R + D2)*R + D1)*R + 1)
            ELSE
               R = R - SPLIT2
               PHINV = (((((((E7*R + E6)*R + E5)*R + E4)*R + E3)
     &              *R + E2)*R + E1)*R + E0) /
     &              (((((((F7*R + F6)*R + F5)*R + F4)*R + F3)
     &              *R + F2)*R + F1)*R + 1)
            END IF
         ELSE
            PHINV = 8.5
         END IF
         IF (Q .LT. 0) PHINV = -PHINV
      END IF
      END
*
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
     &    / 11111111, 22222222, 33333333, 44444444, 55555555, 88888888 /
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
      END
