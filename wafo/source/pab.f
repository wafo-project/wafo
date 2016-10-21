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
     &        SAMPLS, I, INTVLS, MINSMP,NK
      PARAMETER ( PLIM = 25, NLIM = 100, KLIM = 20, MINSMP = 8 )
      INTEGER , DIMENSION(PLIM) :: P
      INTEGER , DIMENSION(PLIM,KLIM-1) ::  C 
      DOUBLE PRECISION FUNCTN, ABSEPS, RELEPS, FINEST, ABSERR, DIFINT, 
     &                 FINVAL, VARSQR, VAREST, VARPRD, VALUE, ONE
      PARAMETER ( ONE = 1.D0 )
      DOUBLE PRECISION, DIMENSION(2*NLIM) :: X
      DOUBLE PRECISION, DIMENSION(KLIM  ) :: VK
      
      SAVE P, C, SAMPLS, NP, VAREST
      INFORM = 1
      INTVLS = 0
      KLIMI = KLIM
      IF ( MINVLS .GE. 0 ) THEN
         FINEST = 0.d0
         VAREST = 0.d0
         SAMPLS = MINSMP 
         DO I = 1, PLIM
            NP = I
            IF ( MINVLS .LT. 2*SAMPLS*P(I) ) GO TO 10
         END DO
         SAMPLS = MAX( MINSMP, MINVLS/( 2*P(NP) ) )
      ENDIF
 10   VK(1) = ONE/DBLE(P(NP))
      NK=MIN( NDIM, KLIM )
      DO I = 2, NK
         VK(I) = MOD( DBLE(C(NP,NK-1))*VK(I-1), ONE )
      END DO
      FINVAL = 0.d0
      VARSQR = 0.d0
      DO I = 1, SAMPLS
         CALL DKSMRC( NDIM, KLIMI, VALUE, P(NP), VK, FUNCTN, X )
         DIFINT = ( VALUE - FINVAL )/DBLE(I)
         FINVAL = FINVAL + DIFINT
         VARSQR = DBLE( I - 2 )*VARSQR/DBLE(I) + DIFINT*DIFINT
      END DO
      INTVLS = INTVLS + 2*SAMPLS*P(NP)
      VARPRD = VAREST*VARSQR
      FINEST = FINEST + ( FINVAL - FINEST )/( 1.d0 + VARPRD )
      IF ( VARSQR .GT. 0.d0 ) VAREST = ( 1.d0 + VARPRD )/VARSQR
      ABSERR = 3.d0*SQRT( VARSQR/( 1.d0 + VARPRD ) )
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
      END  SUBROUTINE DKBVRC
*
      SUBROUTINE DKSMRC( NDIM, KLIM, SUMKRO, PRIME, VK, FUNCTN, X )
      EXTERNAL FUNCTN
      INTEGER NDIM, KLIM, PRIME, K, J, JP, NK
      DOUBLE PRECISION SUMKRO, FUNCTN,  ONE, XT, MVNUNI
      DOUBLE PRECISION ,DIMENSION(:) :: VK,X
      PARAMETER ( ONE = 1.d0 )
      SUMKRO = 0.d0
      NK = MIN( NDIM, KLIM )
      DO J = 1, NK-1
         CALL random_number(MVNUNI)
         JP = J + MVNUNI*( NK + 1 - J ) 
         XT = VK(J)
         VK(J) = VK(JP)
         VK(JP) = XT
      END DO
      DO J = 1, NDIM
         CALL random_number(MVNUNI)
         X(NDIM+J) = MVNUNI
      END DO
      DO K = 1, PRIME
         DO J = 1, NK
            X(J) = MOD( DBLE(K)*VK(J), ONE )
         END DO
         IF ( NDIM. GT. KLIM ) CALL DKRCHT( NDIM-KLIM, X(KLIM+1:NDIM) )
         DO J = 1, NDIM
            XT = X(J) + X(NDIM+J)
            IF ( XT .GT. 1 ) XT = XT - 1.d0
            X(J) = ABS( 2.d0*XT - 1.d0 )
         END DO
         SUMKRO = SUMKRO + ( FUNCTN(NDIM,X) - SUMKRO )/DBLE( 2*K - 1 )
         DO J = 1, NDIM
            X(J) = 1.d0 - X(J)
         END DO
         SUMKRO = SUMKRO + ( FUNCTN(NDIM,X) - SUMKRO )/DBLE( 2*K )
      END DO
      END SUBROUTINE DKBVRC
*
      SUBROUTINE DKRCHT( S, QUASI )
*
*     This subroutine generates a new quasi-random Richtmeyer vector. 
*     A reference is
*      "Methods of Numerical Integration", P.J. Davis and P. Rabinowitz, 
*       Academic Press, 1984, pp. 482-483.
*
*       INPUTS:
*         S - the number of dimensions; 
*             KRRCHT is initialized for each new S or S < 1.
*
*       OUTPUTS:
*         QUASI - a new quasi-random S-vector
*
* revised pab 01.11.1999
* updated to fortran 90
      INTEGER MXDIM, MXHSUM, B
      PARAMETER ( MXDIM = 80, MXHSUM = 48, B = 2 )
      INTEGER S, HISUM, I,  OLDS, 
      INTEGER, DIMENSION(MXDIM   )  :: PRIME
      INTEGER, DIMENSION(0:MXHSUM)  :: N
      DOUBLE PRECISION , DIMENSION(:) :: QUASI,  
      DOUBLE PRECISION , DIMENSION(MXDIM) ::PSQT
      DOUBLE PRECISION ::  ONE, RN
      PARAMETER ( ONE = 1.D0 )
      PARAMETER ( PRIME = (/ 
     &     2,    3,    5,    7,   11,   13,   17,   19,   23,   29,
     &    31,   37,   41,   43,   47,   53,   59,   61,   67,   71,
     &    73,   79,   83,   89,   97,  101,  103,  107,  109,  113,
     &   127,  131,  137,  139,  149,  151,  157,  163,  167,  173,
     &   179,  181,  191,  193,  197,  199,  211,  223,  227,  229,
     &   233,  239,  241,  251,  257,  263,  269,  271,  277,  281,
     &   283,  293,  307,  311,  313,  317,  331,  337,  347,  349,
     &   353,  359,  367,  373,  379,  383,  389,  397,  401,  409/))
* Primes to continue
* 419  421   431   433   439   443   449   457   461   463   467   479   487   491   499
* 503   509   521   523   541   547   557   563   569   571   577   587   593   599
      SAVE OLDS, PSQT, HISUM, N
      DATA OLDS / 0 /
      IF ( S .NE. OLDS .OR. S .LT. 1 ) THEN
         OLDS = S
         N(0) = 0
         HISUM = 0
         DO I = 1, S
            RN = DBLE(PRIME(I))
            PSQT(I) = SQRT( RN )
         END DO
      END IF
      DO I = 0, HISUM 
         N(I) = N(I) + 1
         IF ( N(I) .LT. B ) GO TO 10
         N(I) = 0
      END DO
      HISUM = HISUM + 1
      IF ( HISUM .GT. MXHSUM ) HISUM = 0
      N(HISUM) = 1
 10   RN = 0.d0
      DO I = HISUM, 0, -1
         RN = DBLE(N(I)) + DBLE(B)*RN
      END DO
      DO I = 1, S
         QUASI(I) = MOD( RN*PSQT(I), ONE )
      END DO
      END SUBROUTINE DKRCHT
