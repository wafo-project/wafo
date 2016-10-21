MODULE PRECISION_MODEL
   IMPLICIT NONE
   INTEGER, PARAMETER, PUBLIC   :: STND = SELECTED_REAL_KIND(12, 60)
END MODULE PRECISION_MODEL
!
MODULE MVSTAT
!
!
!   This package contains Fortran 90 software for the MVT distribution, 
!     plus supporting software. The main subroutine is MVDIST. 
!   The package is self-contained and should compile without errors on 
!   standard Fortran(90) compilers. 
!
   USE PRECISION_MODEL
   IMPLICIT NONE
   PUBLIC  :: MVDIST,    MVDSTD,    MVCRIT, LWRUPR, LWUPRH, LWUPRT, COVSCL
!                MVT, MVT+Deriv,  Critical, Bounds, Bounds,
   PUBLIC  :: MVKBRV,  UNIFRM
!            KOROBOV,  Uniform
!
   PRIVATE :: MVINIT, MVINTD, MVSORT, MVSWAP, MVFUNC, MVLIMS
   PRIVATE :: MVKRSV, MVSPCL, MVVRNC, MVCHNT
!
   INTEGER,         PRIVATE, PARAMETER        :: ML = 500, NL = 200
   REAL(KIND=STND), PRIVATE                   :: SQTNU
   REAL(KIND=STND), PRIVATE, DIMENSION(ML,NL) :: CNSTRN
   REAL(KIND=STND), PRIVATE, DIMENSION(ML)    :: A, B, D 
   INTEGER,         PRIVATE, DIMENSION(ML)    :: INF
   INTEGER,         PRIVATE, DIMENSION(ML)    :: CONLNG
   INTEGER,         PRIVATE                   :: NP, MP, NUP
!
   PUBLIC  :: MVDNST,    MVUVT, MVSTDT, MVSTNV
!             T-Dens, 1-D Dist, T-Dist,  Inv T
   PUBLIC  :: MVPHI,  MVPHNV,  MVCHNV, MVCHI
!               Phi, Inv Phi, Inv Chi,   Chi
   PUBLIC  :: MVBVT, MVBVTL, MVTVT, MVTVTL, MVBVU
!               BVT,           TVT,          BVNU
   PUBLIC  :: MVSTDC,  MVBVTC, MVTVTC
!        Comp:   1-D,     2-D,    3-D 
   PUBLIC  :: MVADON,         MVKRND
!       1-D Adaptive,   Kronrod Rule
   PRIVATE :: TVTFNC, MVCHNC
   INTEGER,         PRIVATE :: NUF
   REAL(KIND=STND), PRIVATE :: BF1, BF2, BF3, RF31, RF32
CONTAINS
   SUBROUTINE MVDIST( N, COVRNC, NU, M, LOWER, CONSTR, UPPER, INFIN, DELTA,   &
                      MAXPTS, ABSEPS, RELEPS, ERROR, VALUE, NEVALS, INFORM )
!
!     A subroutine for computing non-central multivariate t probabilities.
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
!  Parameters
!
!     N      INTEGER number of variables.
!     COVRNC REAL NxN covariance matrix; only the diagonal and lower left 
!             triangle are used. COVRNC matrix must be positive semidefinite.
!     NU     INTEGER, number of degrees of freedom; if NU < 1, then
!             a multivariate normal computation is completed. 
!     M      INTEGER, number of linear constraints for integration region.
!     LOWER  REAL array of lower integration limits.
!     CONSTR REAL MxN constraint matrix; integration region is all X with
!                     LOWER < DELTA + CONSTR*X < UPPER.
!     UPPER  REAL, array of upper integration limits.
!     INFIN  INTEGER, array of integration limits flags:
!            if INFIN(I) < 0, Ith limits are (-infinity, infinity);
!            if INFIN(I) = 0, Ith limits are (-infinity, UPPER(I)];
!            if INFIN(I) = 1, Ith limits are [LOWER(I), infinity);
!            if INFIN(I) = 2, Ith limits are [LOWER(I), UPPER(I)].
!     DELTA  REAL array of non-centrality parameters.
!     MAXPTS INTEGER, maximum number of function values allowed. This 
!            parameter can be used to limit the time. A sensible 
!            strategy is to start with MAXPTS = 1000*N, and then
!            increase MAXPTS if ERROR is too large.
!     ABSEPS REAL absolute error tolerance.
!     RELEPS REAL relative error tolerance.
!     ERROR  REAL estimated absolute error, with 99% confidence level.
!     VALUE  REAL estimated value for the integral
!     NEVALS Integer number of integrand values used for this call.
!     INFORM INTEGER, termination status parameter:
!            if INFORM = 0, normal completion with ERROR < EPS;
!            if INFORM = 1, completion with ERROR > EPS and MAXPTS 
!                           function vaules used; increase MAXPTS to 
!                           decrease ERROR;
!            if INFORM = 2, N > NL or N < 1.
!            if INFORM = 3, covariance matrix not positive semidefinite.
!
! Parameter Variables
!
      INTEGER,                         INTENT(IN)  :: N, NU, M, MAXPTS
      REAL(KIND=STND), DIMENSION(:),   INTENT(IN)  :: LOWER, UPPER, DELTA
      REAL(KIND=STND), DIMENSION(:,:), INTENT(IN)  :: COVRNC, CONSTR
      INTEGER,    DIMENSION(:),        INTENT(IN)  :: INFIN
      REAL(KIND=STND),                 INTENT(IN)  :: RELEPS, ABSEPS
      INTEGER,                         INTENT(OUT) :: INFORM, NEVALS
      REAL(KIND=STND),                 INTENT(OUT) :: ERROR, VALUE
!
!
!
      REAL(KIND=STND), DIMENSION(1) :: V
     
!
      NEVALS = 0
      CALL MVINIT( N, COVRNC, NU, M, LOWER, CONSTR, UPPER, INFIN, DELTA,      &
                   VALUE, ERROR, INFORM )
      IF ( INFORM == 0 .AND. NP > 0 ) THEN 
         CALL MVKBRV( NP - 1 + MIN(1,NUP), 0, MAXPTS, 1, MVFUNC,              &
                      ABSEPS, RELEPS, ERROR, V, NEVALS, INFORM )
         VALUE = V(1)
      ENDIF
!
   END SUBROUTINE MVDIST
!
   FUNCTION MVFUNC( NF, W ) RESULT(VALUE)
!     
!     Integrand subroutine
!
!
! Global Variables
!
      INTEGER,         INTENT(IN)               :: NF
      REAL(KIND=STND), INTENT(IN), DIMENSION(:) :: W
      REAL(KIND=STND),            DIMENSION(NF) :: VALUE
!
! Local Variables  
!
      REAL(KIND=STND), DIMENSION(NP) :: Y
      INTEGER                        :: I, J, L, LI, FINITA, FINITB
      REAL(KIND=STND)                :: R, S, AI, AT, BI, BT, DI, EI, DE
!
      VALUE(1:NF) = 0
      IF ( NUP > 0 ) THEN
         CALL MVCHNT( NUP, W(NP+NF-1), R, VALUE(1) )
         R = R/SQTNU
      ELSE
         VALUE(1) = 1
         R = 1
      END IF
      LI = 0
      DO I = 1, NP
         FINITA = 0
         FINITB = 0
         DO L = LI + 1, LI + CONLNG(I)         
            S = D(L) + SUM( CNSTRN(L,1:I-1)*Y(1:I-1) ) 
            IF ( INF(L) /= 0 ) THEN
               AI = R*A(L) - S 
               IF ( FINITA > 0 ) AI = MAX( AI, AT )
               FINITA = 1
               AT = AI
            END IF
            IF ( INF(L) /= 1 ) THEN
               BI = R*B(L) - S 
               IF ( FINITB > 0 ) BI = MIN( BI, BT )
               FINITB = 1
               BT = BI
            END IF
         END DO
         CALL MVLIMS( 0, AI, BI, FINITA + FINITA + FINITB - 1, DI, EI, DE )
         VALUE(1) = VALUE(1)*DE
         IF ( DE == 0 ) EXIT
         IF ( I < NP .OR. NF == 2 ) Y(I) = MVPHNV( DI + W(I)*DE )
         LI = LI + CONLNG(I)
      END DO
      IF ( VALUE(1) > 0 .AND. NF == 2 ) VALUE(2) = VALUE(1)*( NP - SUM(Y*Y) )
!
   END FUNCTION MVFUNC
!
!
   SUBROUTINE MVCHNT( N, P, R, S )    
!
!                  R
!     P =  1 - K  I     exp(-t*t/2) t**(N-1) dt, for N >= 1.
!               N  0
!
!   with Normal approximation when N >= NR
!
      INTEGER,         INTENT(IN)  :: N      
      REAL(KIND=STND), INTENT(IN)  :: P
      REAL(KIND=STND), INTENT(OUT) :: R
      REAL(KIND=STND), INTENT(OUT) :: S ! Integrand scale factor for N >= 100
!
      INTEGER,          PARAMETER  :: NR = 1000
      INTEGER                      :: I
!
      REAL(KIND=STND), PARAMETER   :: LRP = -.22579135264472743235E0_STND
!                                     LRP =   LOG( SQRT( 2/PI ) )
      REAL(KIND=STND), PARAMETER   :: LRT =  .91893853320467274177E0_STND
!                                     LRT =   LOG( SQRT( 2*PI ) )

      REAL(KIND=STND)              :: RI
      REAL(KIND=STND),       SAVE  :: LKN
      INTEGER,               SAVE  :: NO = 0
!
      IF ( N < NR ) THEN
         R = MVCHNV( N, P )
         S = 1
      ELSE
!
!     Use Normal approximation for Chi
!
         IF ( N /= NO ) THEN
!
!          First call computes LKN = LOG(K )
!                                         N
           NO = N
            LKN = 0
            DO I = N - 2, 2, -2
               RI = I
               LKN = LKN - LOG( RI )
            END DO
            IF ( MODULO( N, 2 ) == 1 ) LKN = LKN + LRP
         END IF
         RI = N - 1
         RI = SQRT(RI)
         R = RI - MVPHNV(P)
         S = EXP( LRT + LKN + RI*( RI*LOG(R) - R + RI/2 ) )
!
!        S = K  R**(N-1)*EXP( -( R*R - (R-SQRT(N-1))**2 )/2 )*SQRT(2*PI) 
!             N 
      END IF
!
   END SUBROUTINE MVCHNT
!
   SUBROUTINE MVINIT( N, COVRNC, NU, M, LOWER, CONSTR, UPPER, INFIN, DELTA,   &
                      VALUE, ERROR, INFORM )
!
!     Initialization and computation of covariance Cholesky factor.
!
!
! Parameter Variables
!
      INTEGER,                         INTENT(IN)  :: N, NU, M
      REAL(KIND=STND), DIMENSION(N,N), INTENT(IN)  :: COVRNC
      REAL(KIND=STND), DIMENSION(M),   INTENT(IN)  :: LOWER, UPPER, DELTA
      REAL(KIND=STND), DIMENSION(M,N), INTENT(IN)  :: CONSTR
      INTEGER,         DIMENSION(M),   INTENT(IN)  :: INFIN
      REAL(KIND=STND),                 INTENT(OUT) :: VALUE, ERROR
      INTEGER,                         INTENT(OUT) :: INFORM
!
!     Initialization and computation of covariance Cholesky factor.
!
!
      IF ( 0 < N .AND. N <= NL ) THEN
         NUP = MAX( NU, 0 )
         CALL MVSORT( N, COVRNC, M, LOWER, CONSTR, UPPER, INFIN, DELTA, &
                      .TRUE., INFORM )
      ELSE
         INFORM = 2
      END IF
      CALL MVSPCL( NU, VALUE, ERROR, INFORM )
!
   END SUBROUTINE MVINIT
!
   SUBROUTINE MVSPCL( NU, VALUE, ERROR, INFORM )
!
!     Special cases subroutine
!
! Parameter Variables
!
      INTEGER,         INTENT(IN)  :: NU, INFORM
      REAL(KIND=STND), INTENT(OUT) :: VALUE, ERROR
!
! Local Variables  
!
      REAL(KIND=STND) :: DS, R22, R33
!
      IF ( INFORM > 0 ) THEN
         VALUE = 0
         ERROR = 1
      ELSE
         SQTNU = SQRT( REAL( MAX( NUP, 0 ), STND ) ) 
!     
!        Special cases
!
         IF ( NP == 0 ) THEN
            ERROR = 0
            VALUE = 1
         ELSE IF ( NP < 4 ) THEN
            DS = SUM ( ABS( D( 1:CONLNG(1) ) ) )
            IF ( NP == 1 .AND. ( NUP < 1 .OR. DS == 0 ) ) THEN
!     
!              1-d case for normal or central t
!

               VALUE = MVUVT( NUP, A(1)-D(1), B(1)-D(1), INF(1) ) 
               ERROR = 2E-16_STND
               NP = 0
            ELSE IF ( NP == 2 .AND. CONLNG(1) == 1 .AND. CONLNG(2) < 3 ) THEN
               IF ( NUP < 1 .OR. SUM( ABS(D(1:1+CONLNG(2))) ) == 0 ) THEN
                  IF ( INF(1) /= 0 ) A(1) = A(1) - D(1) 
                  IF ( INF(1) /= 1 ) B(1) = B(1) - D(1) 
                  R22 = SQRT( 1 + CNSTRN(2,1)**2 )
                  IF ( INF(2) /= 0 ) A(2) = ( A(2) - D(2) )/R22
                  IF ( INF(2) /= 1 ) B(2) = ( B(2) - D(2) )/R22
                  IF ( CONLNG(2) == 1 ) THEN
!     
!                    2-d case for normal or central t
!
                     CNSTRN(2,1) = CNSTRN(2,1)/R22
                     VALUE = MVBVT( NUP, A, B, INF, CNSTRN(2,1))
                     ERROR = 1E-15_STND
                  ELSE
!     
!                    3-d singular case for normal or central t
!
                     R33 = SQRT( 1 + CNSTRN(3,1)**2 )
                     IF ( INF(3) /= 0 ) A(3) = ( A(3) - D(3) )/R33
                     IF ( INF(3) /= 1 ) B(3) = ( B(3) - D(3) )/R33
                     CNSTRN(3,2) = CNSTRN(3,2) + CNSTRN(2,1)*CNSTRN(3,1)
                     CNSTRN(2,1) = CNSTRN(2,1)/R22
                     CNSTRN(3,1) = CNSTRN(3,1)/R33
                     CNSTRN(3,2) = CNSTRN(3,2)/( R33*R22 )
                     VALUE = MVTVT( NUP, A, B, INF,                          &
                                   (/CNSTRN(2,1),CNSTRN(3,1),CNSTRN(3,2)/) )
                     ERROR = 1E-14_STND
                  END IF
                  NP = 0
               END IF
            ELSE IF ( NP == 3 .AND. CONLNG(1) == 1 .AND. CONLNG(2) == 1      &
                              .AND. CONLNG(3) == 1 ) THEN
               IF ( NUP < 1 .OR. SUM( ABS(D(1:3)) ) == 0 ) THEN
!     
!                 3-d case for normal or central t
!
                  IF ( INF(1) /= 0 ) A(1) = A(1) - D(1) 
                  IF ( INF(1) /= 1 ) B(1) = B(1) - D(1) 
                  R22 = SQRT( 1 + CNSTRN(2,1)**2 )
                  IF ( INF(2) /= 0 ) A(2) = ( A(2) - D(2) )/R22
                  IF ( INF(2) /= 1 ) B(2) = ( B(2) - D(2) )/R22
                  R33 = SQRT( 1 + CNSTRN(3,1)**2 + CNSTRN(3,2)**2 )
                  IF ( INF(3) /= 0 ) A(3) = ( A(3) - D(3) )/R33
                  IF ( INF(3) /= 1 ) B(3) = ( B(3) - D(3) )/R33
                  CNSTRN(3,2) = CNSTRN(3,2) + CNSTRN(2,1)*CNSTRN(3,1)
                  CNSTRN(2,1) = CNSTRN(2,1)/R22
                  CNSTRN(3,1) = CNSTRN(3,1)/R33
                  CNSTRN(3,2) = CNSTRN(3,2)/( R33*R22 )
                  VALUE = MVTVT( NUP, A, B, INF,                            &
                                 (/CNSTRN(2,1),CNSTRN(3,1),CNSTRN(3,2)/) )
                  ERROR = 1E-14_STND
                  NP = 0
               END IF
            END IF
         END IF
      END IF
!
   END SUBROUTINE MVSPCL
!
   SUBROUTINE MVLIMS( NU, A, B, INFIN, LOWER, UPPER, VALUE )
!
! Parameter Variables
!
      INTEGER,         INTENT(IN)  :: NU, INFIN
      REAL(KIND=STND), INTENT(IN)  :: A, B
      REAL(KIND=STND), INTENT(OUT) :: LOWER, UPPER, VALUE
!
! Local Variables  
!
      LOWER = 0
      UPPER = 1
      IF ( INFIN >= 0 ) THEN
         IF ( INFIN /= 0 ) LOWER = MVSTDT( NU, A )
         IF ( INFIN /= 1 ) UPPER = MVSTDT( NU, B )
      END IF
      UPPER = MAX( UPPER, LOWER )
      VALUE = UPPER - LOWER
!
      END SUBROUTINE MVLIMS
!
   SUBROUTINE MVDSTD( N, COVRNC, NU, M, LOWER, CONSTR, UPPER, INFIN, T,       &
                MAXPTS, ABSEPS, RELEPS, ERROR, VALUE, DERIV, NEVALS, INFORM )
!
!     A subroutine for computing multivariate t probability and derivative.
!     This subroutine uses an algorithm (QRSVN) described in the paper
!     "Methods for the  Computation of Multivariate t-Probabilities",
!        by Alan Genz and Frank Bretz
!
!          Alan Genz 
!          Department of Mathematics
!          Washington State University 
!          Pullman, WA 99164-3113
!          Email : AlanGenz@wsu.edu
!
!  Parameters
!
!     N      INTEGER number of variables.
!     COVRNC REAL NxN covariance matrix; only the diagonal and lower left 
!             triangle are used. COVRNC matrix must be positive semidefinite.
!     NU     INTEGER, number of degrees of freedom; if NU < 1, then
!             a multivariate normal computation is completed. 
!     M      INTEGER, number of linear constraints for integration region.
!     LOWER  REAL array of lower integration limits.
!     CONSTR REAL MxN constraint matrix; constraint I is
!                     T*LOWER < DELTA + CONSTR*X < T*UPPER.
!     UPPER  REAL, array of upper integration limits.
!     INFIN  INTEGER, array of integration limits flags:
!            if INFIN(I) < 0, Ith limits are (-infinity, infinity);
!            if INFIN(I) = 0, Ith limits are (-infinity, UPPER(I)];
!            if INFIN(I) = 1, Ith limits are [LOWER(I), infinity);
!            if INFIN(I) = 2, Ith limits are [LOWER(I), UPPER(I)].
!     T      REAL, scale factor for integration limits.
!     MAXPTS INTEGER, maximum number of function values allowed. This 
!            parameter can be used to limit the time. A sensible 
!            strategy is to start with MAXPTS = 1000*N, and then
!            increase MAXPTS if ERROR is too large.
!     ABSEPS REAL absolute error tolerance.
!     RELEPS REAL relative error tolerance.
!     ERROR  REAL estimated absolute error, with 99% confidence level.
!     VALUE  REAL estimated value for the integral
!     DERIV  REAL estimated value for derivative with respect to T
!     NEVALS Integer number of integrand values used for this call.
!     INFORM INTEGER, termination status parameter:
!            if INFORM = 0, normal completion with ERROR < EPS;
!            if INFORM = 1, completion with ERROR > EPS and MAXPTS 
!                           function vaules used; increase MAXPTS to 
!                           decrease ERROR;
!            if INFORM = 2, N > NL or N < 1.
!            if INFORM = 3, covariance matrix not positive semidefinite.
!
!
!
! Parameter Variables
!
      INTEGER,                         INTENT(IN)  :: N, NU, M, MAXPTS
      REAL(KIND=STND), DIMENSION(N,N), INTENT(IN)  :: COVRNC
      REAL(KIND=STND), DIMENSION(M),   INTENT(IN)  :: LOWER, UPPER
      REAL(KIND=STND), DIMENSION(M,N), INTENT(IN)  :: CONSTR
      INTEGER,         DIMENSION(M),   INTENT(IN)  :: INFIN
      REAL(KIND=STND),                 INTENT(IN)  :: T, RELEPS, ABSEPS
      INTEGER,                         INTENT(OUT) :: INFORM, NEVALS
      REAL(KIND=STND),                 INTENT(OUT) :: ERROR, VALUE, DERIV
!
      REAL(KIND=STND), DIMENSION(2) :: V
!
      NEVALS = 0
      CALL MVINTD( N, COVRNC, NU, M, LOWER, CONSTR, UPPER, INFIN, T,          &
                   VALUE, DERIV, ERROR, INFORM )
      IF ( INFORM == 0 .AND. NP > 0 ) THEN 
         CALL MVKBRV( NP + MIN(1,NUP), 0, MAXPTS, 2, MVFUNC,                  &
                      ABSEPS, RELEPS, ERROR, V, NEVALS, INFORM )
         VALUE = V(1)
         DERIV = V(2)/T
      ENDIF
!
      END SUBROUTINE MVDSTD
!
   SUBROUTINE MVINTD( N, COVRNC, NU, M, LOWER, CONSTR, UPPER, INFIN, T,   &
                      VALUE, DERIV, ERROR, INFORM )
!
!    Initialization and computation of covariance Cholesky factor.
!
!
! Parameter Variables
!
      INTEGER,                         INTENT(IN)  :: N, NU, M
      REAL(KIND=STND), DIMENSION(N,N), INTENT(IN)  :: COVRNC
      REAL(KIND=STND), DIMENSION(M),   INTENT(IN)  :: LOWER, UPPER
      REAL(KIND=STND), DIMENSION(M,N), INTENT(IN)  :: CONSTR
      INTEGER,         DIMENSION(M),   INTENT(IN)  :: INFIN
      REAL(KIND=STND),                 INTENT(IN)  :: T
      REAL(KIND=STND),                 INTENT(OUT) :: ERROR, VALUE, DERIV
      INTEGER,                         INTENT(OUT) :: INFORM
!
! Local Variables  
!
      REAL(KIND=STND), DIMENSION(M) :: DELTA, AT, BT
      INTEGER                       :: I
!
!     Initialization and computation of covariance Cholesky factor.
!
!
      IF ( 0 < N .AND. N <= NL ) THEN
         NUP = NU
         DELTA = 0
         DO I = 1, M
            IF ( INFIN(I) /= 0 ) AT(I) = T*LOWER(I)
            IF ( INFIN(I) /= 1 ) BT(I) = T*UPPER(I)
         END DO
         CALL MVSORT( N, COVRNC, M, AT, CONSTR, BT, INFIN, DELTA,           &
                      .FALSE., INFORM )
!
         IF ( INFORM > 0 ) THEN
            VALUE = 0
            DERIV = 0
            ERROR = 1
         ELSE
            SQTNU = SQRT( REAL( MAX( NUP, 0 ), STND ) ) 
            VALUE = 1
            DERIV = 0 
            IF ( NP .EQ. 0 ) THEN
               ERROR = 0
            ELSE IF ( NP .EQ. 1 ) THEN
               IF ( INF(1) /= 1 ) THEN
                  VALUE =         MVSTDT( NU, B(1) ) 
                  DERIV =         B(1)*MVDNST( NU, B(1) )
               END IF
               IF ( INF(1) /= 0 ) THEN
                  VALUE = VALUE - MVSTDT( NU, A(1) ) 
                  DERIV = DERIV - A(1)*MVDNST( NU, A(1) )
               END IF
               ERROR = 1E-16_STND
            END IF
         END IF
      ELSE
         INFORM = 2
      END IF
!
   END SUBROUTINE MVINTD
!
   SUBROUTINE MVSORT( N, COVRNC, M, LOWER, CONSTR, UPPER, INFIN, DELTA,    &
                      PIVOT, INFORM )
!
!     Subroutine to sort integration limits and determine Cholesky factor.
!
!
! Parameter Variables
!
      INTEGER,         INTENT(IN)                    :: N, M
      REAL(KIND=STND), DIMENSION(N,N), INTENT(IN)    :: COVRNC
      REAL(KIND=STND), DIMENSION(M),   INTENT(IN)    :: LOWER, UPPER, DELTA
      REAL(KIND=STND), DIMENSION(M,N), INTENT(IN)    :: CONSTR
      INTEGER,         DIMENSION(M),   INTENT(IN)    :: INFIN
      LOGICAL,                         INTENT(IN)    :: PIVOT
      INTEGER,                         INTENT(INOUT) :: INFORM
!
! Local Variables  
!
      REAL(KIND=STND), DIMENSION(N)   :: V, Y 
      REAL(KIND=STND), DIMENSION(N,N) :: CHLFAC
      INTEGER,         DIMENSION(M)   :: CONLIM
      INTEGER                         :: I, J, K, L, LMN, JMX, INFT, IA, IB
      REAL(KIND=STND)                 :: CVDIAG, S, SUMSQ, AT, BT, DT, EPSI
      REAL(KIND=STND)                 :: AL, BL, VL, VARMN, YL
      REAL(KIND=STND), PARAMETER      :: ONE = 1, EPS = 10000*EPSILON(ONE)
      INFORM = 0
!
!     Copy parameters and drop any doubly infinite constraints.
!
      MP = 0
      DO I = 1, M
         IF ( INFIN(I) >= 0 ) THEN
            MP = MP + 1
            A(MP) = 0
            B(MP) = 0
            INF(MP) = INFIN(I) 
            IF ( INF(MP) /= 0 ) A(MP) = LOWER(I)
            IF ( INF(MP) /= 1 ) B(MP) = UPPER(I)
            CNSTRN(MP,1:N) = CONSTR(I,1:N)
            D(MP) = DELTA(I)
         END IF
      END DO
      CHLFAC = 0
      DO I = 1, N
         CHLFAC(I:N,I) = COVRNC(I:N,I)
      END DO
!
!     Determine Cholesky factor revised contraint matrix.
!
      NP = 0
      DO I = 1, N
         CONLNG(I) = 0
         EPSI = EPS*I*I 
         JMX = I
         DO L = I+1, N
            IF ( CHLFAC(L,L) > CHLFAC(JMX,JMX) ) JMX = L
         END DO
         IF ( JMX > I ) THEN
            CALL MVSWAP( CHLFAC(I,I), CHLFAC(JMX,JMX) )
            DO J = 1, I-1
               CALL MVSWAP( CHLFAC(I,J), CHLFAC(JMX,J) )
            END DO
            DO J = I+1, JMX-1
               CALL MVSWAP( CHLFAC(J,I), CHLFAC(JMX,J) )
            END DO
            DO J = JMX+1, N
               CALL MVSWAP( CHLFAC(J,I), CHLFAC(J,JMX) )
            END DO
            DO L = 1, MP
               CALL MVSWAP( CNSTRN(L,I), CNSTRN(L,JMX) )
            END DO
        END IF
!
!        Compute Ith columns of Cholesky factor and constraint matrix.
!
         IF ( CHLFAC(I,I) < -EPSI ) INFORM = 3
         IF ( CHLFAC(I,I) <  EPSI ) EXIT
         CVDIAG = SQRT( CHLFAC(I,I) )
         CHLFAC(I,I) = CVDIAG
         DO L = I+1, N         
            CHLFAC(L,I) = CHLFAC(L,I)/CVDIAG
            CHLFAC(L,I+1:L) = CHLFAC(L,I+1:L) - CHLFAC(L,I)*CHLFAC(I+1:L,I)
         END DO
         IF( MP > 0 ) THEN
            CNSTRN(1:MP,I) = MATMUL( CNSTRN(1:MP,I:N), CHLFAC(I:N,I) ) 
         END IF
         NP = NP + 1
      END DO
      IF ( INFORM == 0 ) THEN
!
!        Use right reflectors to reduce CNSTRN to lower triangular form
!
         DO I = 1, MIN( NP-1, MP )
            EPSI = EPS*I*I
            IF ( PIVOT ) THEN
!              
!              Permute rows so that smallest variance variables are first
!
               VARMN = 1
               DO L = I, MP
                  V(I:NP) = CNSTRN(L,I:NP)
                  SUMSQ = MAX( SQRT( SUM( V(I:NP)*V(I:NP) ) ), EPSI )
                  S = D(L) + SUM( CNSTRN(L,1:I-1)*Y(1:I-1) )
                  AL = ( A(L) - S )/SUMSQ
                  BL = ( B(L) - S )/SUMSQ
                  CALL MVVRNC( AL, BL, INF(L), EPSI, YL, VL )
                  IF ( VL <= VARMN ) THEN
                     LMN   = L
                     VARMN = VL
                     Y(I) = YL
                  END IF
               END DO
               V(1:NP) = CNSTRN(LMN,1:NP) 
               IF ( LMN > I ) THEN
                  CNSTRN(LMN,1:NP) = CNSTRN(I,1:NP)
                  CNSTRN(  I,1:NP) =        V(1:NP)
                  CALL MVSWAP( A(I), A(LMN) )
                  CALL MVSWAP( B(I), B(LMN) )
                  CALL MVSWAP( D(I), D(LMN) )
                  INFT      = INF(I)
                  INF(I)    = INF(LMN)
                  INF(LMN) = INFT
               END IF
            ELSE
               V(I:NP) = CNSTRN(I,I:NP) 
            END IF
            CNSTRN(I,I+1:NP) = 0
            SUMSQ = SUM( V(I+1:NP)*V(I+1:NP) )
            IF ( SUMSQ > EPSI ) THEN
               SUMSQ = SQRT( SUMSQ + V(I)*V(I) )
               IF ( V(I) < 0 ) SUMSQ = -SUMSQ 
               CNSTRN(I,I) = -SUMSQ
               V(I) = V(I) + SUMSQ
               BT = 1/( SUMSQ*V(I) )
               DO L = I+1, MP
                  AT = BT*SUM( CNSTRN(L,I:NP)*V(I:NP) )
                  CNSTRN(L,I:NP) = CNSTRN(L,I:NP) - AT*V(I:NP)
               END DO
            END IF
         END DO
!
!        Scale and sort constraints
!
         DO I = 1, MP
            Y(1:NP) = CNSTRN(I,1:NP)
            CONLIM(I) = MIN(I,NP)
            JMX = 1
            DO J = 1, CONLIM(I) 
               IF( ABS(Y(J)) > EPS*J*J ) JMX = J
            END DO
            Y(JMX+1:NP) = 0
            CONLNG(JMX) = CONLNG(JMX) + 1 
            AT   = A(I)
            BT   = B(I)
            DT   = D(I)
            INFT = INF(I)
            J = I
            DO L = I-1, 1, -1
               IF( JMX >= CONLIM(L) ) EXIT
               CNSTRN(L+1,1:NP) = CNSTRN(L,1:NP)
               A(L+1)      = A(L)
               B(L+1)      = B(L)
               D(L+1)      = D(L)
               INF(L+1)    = INF(L)
               CONLIM(L+1) = CONLIM(L)
               J = L
            END DO
            CONLIM(J) = JMX
            A(J)   = AT/Y(JMX)
            B(J)   = BT/Y(JMX)
            D(J)   = DT/Y(JMX)
            INF(J) = INFT
            CNSTRN(J,1:NP) = Y(1:NP)/Y(JMX)
            IF ( Y(JMX) < 0 ) THEN
               CALL MVSWAP( A(J), B(J) )
               IF ( INFT /= 2 ) INF(J) = 1 - INFT
            END IF
         END DO
         JMX = 0
         DO I = 1, NP
            IF( CONLNG(I) > 0 ) JMX = I
         END DO
         NP = JMX
         IF ( CONLNG(1) > 1 ) THEN
            IF ( NUP < 1 .OR. SUM( ABS( D(1:CONLNG(1) ) ) ) == 0 ) THEN
!
!              Combine first variable constraints
!
               IA = 0
               IB = 0
               DO L = 1, CONLNG(1)         
                  IF ( INF(L) /= 0 ) THEN
                     AL = A(L) - D(L)
                     IF ( IA > 0 ) AL = MAX( AL, AT )
                     IA = 1
                     AT = AL
                  END IF
                  IF ( INF(L) /= 1 ) THEN
                     BL = B(L) - D(L)
                     IF ( IB > 0 ) BL = MIN( BL, BT )
                     IB = 1
                     BT = BL
                  END IF
               END DO
               A(1)   = AL
               B(1)   = BL
               D(1)   = 0
               INF(1) = IA + IA + IB - 1
               INFT = CONLNG(1) - 1
               DO I = 2, MP - INFT
                  A(I)   = A(I+INFT)
                  B(I)   = B(I+INFT)
                  D(I)   = D(I+INFT)
                  INF(I) = INF(I+INFT)
               END DO
               CONLNG(1) = 1
            END IF
         END IF
      END IF
!
   END SUBROUTINE MVSORT
!
   SUBROUTINE MVVRNC( A, B, INF, EPSI, MEAN, VARANC )
!
!      Compute truncated normal mean and variance
!        mean = -( phi(b) - phi(a) )/( Phi(b) - Phi(a) )
!        variance = 1 -( phi(b)b - phi(a)a )/( Phi(b) - Phi(a) ) - mean^2
!
! Parameter Variables
!
      INTEGER,         INTENT(IN)  :: INF
      REAL(KIND=STND), INTENT(IN)  :: A, B, EPSI
      REAL(KIND=STND), INTENT(OUT) :: MEAN, VARANC
!
! Local Variables  
!
      REAL(KIND=STND) :: DENSA, DENSB, DISTA, DISTB, DISTDF 
!
      CALL MVLIMS( 0, A, B, INF, DISTA, DISTB, DISTDF )
      IF ( INF /= 0 ) DENSA = MVDNST( 0, A )
      IF ( INF /= 1 ) DENSB = MVDNST( 0, B )
      IF ( DISTDF .GT. EPSI ) THEN
         IF ( INF == 0 ) THEN 
            MEAN   =        -   DENSB
            VARANC =        - B*DENSB
         ELSE IF ( INF == 1 ) THEN 
            MEAN   =   DENSA
            VARANC = A*DENSA
         ELSE IF ( INF == 2 ) THEN 
            MEAN   =   DENSA -   DENSB 
            VARANC = A*DENSA - B*DENSB
         END IF
         MEAN   =       MEAN/DISTDF
         VARANC = 1 + VARANC/DISTDF - MEAN**2
      ELSE
         IF ( INF == 0 ) MEAN = B
         IF ( INF == 1 ) MEAN = A
         IF ( INF == 2 ) MEAN = ( A + B )/2
         VARANC = 0
      END IF
!
   END SUBROUTINE MVVRNC
!
   SUBROUTINE MVSWAP( X, Y )
!
! Parameter Variables
!
      REAL(KIND=STND), INTENT(INOUT) :: X, Y
!
! Local Variables  
!
      REAL(KIND=STND) :: T
!
      T = X
      X = Y
      Y = T
!
   END SUBROUTINE MVSWAP
!
   SUBROUTINE MVCRIT( N, COVRNC, NU, M, LOWER, CONSTR, UPPER, INFIN, ALPHA,   &
                       MAXPTS, ABSEPS, ERROR, TALPHA, NVS, INFORM )
!
!     A subroutine for computing critical values. The subroutine tries to
!     compute the critical value TALPHA satisfying P(TALPHA) = 1 - ALPHA.
!     P(T) is a multivariate T probability integral, with covariance 
!     matrix COVRNC, for an integration region determined by 
!             TALPHA*LOWER < DELTA + CONSTR*X < TALPHA*UPPER .
!
!       Reference: Alan Genz and Frank Bretz, "Numerical Computation of 
!         Critical Values for Multiple Comparison Problems", pp. 84-87 
!         in 2000 Proceedings of the Statistical Computing Section, 
!         American Statistical Association, Alexandria, VA. 
!
!          Alan Genz 
!          Department of Mathematics
!          Washington State University 
!          Pullman, WA 99164-3113
!          Email : AlanGenz@wsu.edu
!
!  Parameters
!
!     N      INTEGER number of variables.
!     COVRNC REAL NxN covariance matrix; only the diagonal and lower left 
!             triangle are used. COVRNC matrix must be positive semidefinite.
!     NU     INTEGER, number of degrees of freedom; if NU < 1, then
!             a multivariate normal computation is completed. 
!     M      INTEGER, number of linear constraints for integration region.
!     LOWER  REAL array of lower integration limit scale factors.
!     CONSTR REAL MxN constraint matrix; integration region is all X with
!                     LOWER < CONSTR*X < UPPER.
!     UPPER  REAL, array of upper integration limit scale factors.
!     INFIN  INTEGER, array of integration limits flags:
!            if INFIN(I) = 0, Ith limits are (-infinity, UPPER(I)];
!            if INFIN(I) = 1, Ith limits are [LOWER(I), infinity);
!            if INFIN(I) = 2, Ith limits are [LOWER(I), UPPER(I)].
!     ALPHA  REAL confidence factor.
!     MAXPTS INTEGER, maximum number of function values allowed. This 
!            parameter can be used to limit the time. A sensible 
!            strategy is to start with MAXPTS = 1000*N, and then
!            increase MAXPTS if ERROR is too large.
!     ABSEPS REAL absolute error tolerance for ALPHA.
!     ERROR  REAL estimated absolute error, with 99% confidence level.
!     TALPHA REAL estimated critical value.
!     NVS    Integer number of integrand values used for this call.
!     INFORM INTEGER, termination status parameter:
!            if INFORM = 0, normal completion with ERROR < EPS;
!            if INFORM = 1, completion with ERROR > EPS and MAXPTS 
!                           function values used; increase MAXPTS to 
!                           decrease ERROR;
!            if INFORM = 2  for invalid bounds; MVCRIT requires
!                              LOWER < DELTA < UPPER 
!
!
! Parameter Variables
!
      INTEGER,                         INTENT(IN)  :: N, NU, M, MAXPTS
      REAL(KIND=STND), DIMENSION(N,N), INTENT(IN)  :: COVRNC
      REAL(KIND=STND), DIMENSION(M),   INTENT(IN)  :: LOWER, UPPER
      REAL(KIND=STND), DIMENSION(M,N), INTENT(IN)  :: CONSTR
      INTEGER,         DIMENSION(M),   INTENT(IN)  :: INFIN
      REAL(KIND=STND),                 INTENT(IN)  :: ALPHA, ABSEPS
      INTEGER,                         INTENT(OUT) :: INFORM, NVS
      REAL(KIND=STND),                 INTENT(OUT) :: ERROR, TALPHA
!
! Local Variables
!
      REAL(KIND=STND), DIMENSION(M)   :: DL 
      REAL(KIND=STND), DIMENSION(N)   :: CN 
      REAL(KIND=STND), DIMENSION(2)   :: BD, TD
      REAL(KIND=STND)                 :: AT, BT, ATMN, BTMN, ATMX, BTMX  
      REAL(KIND=STND)                 :: TA, TB, TC, TBF, PA, PB, PC, PCB
      REAL(KIND=STND)                 :: DI, RE, ER, VL, DF
      INTEGER                         :: I, INFI, IMAX, IMIN, NV
      INTEGER                         :: IP = 0 ! > 0 Prints for debugging
!
!     Use 1st order bounds to find interval for TALPHA.
!
      PA = 1
      PB = 0
      INFORM = 0
      NVS = 0
      DO I = 1, M
         CN = CONSTR(I,1:N)
         DI = SQRT( SUM( MATMUL( CN, COVRNC )*CN ) )
         IF ( INFIN(I) /= 0 ) AT = LOWER(I)/DI
         IF ( INFIN(I) /= 1 ) BT = UPPER(I)/DI
         PC = MVSTDC( NU, AT, BT, INFIN(I) )
         IF (  PC <= PA ) THEN
            IMIN = I
            PA = PC
            ATMN = AT
            BTMN = BT
         END IF
         IF (  PC >= PB ) THEN
            IMAX = I
            PB = PC
            ATMX = AT
            BTMX = BT
         END IF
         IF ( INFIN(I) /= 0 .AND. AT >= 0 ) INFORM = 2 
         IF ( INFIN(I) /= 1 .AND. BT <= 0 ) INFORM = 2 
      END DO
      IF ( INFORM == 2 ) THEN
         TALPHA = 0
         ERROR = 1
      ELSE
         IF ( INFIN(IMIN) == 0 ) THEN
            TA =  MVSTNV( NU, 1 - ALPHA )/BTMN
         ELSE IF ( INFIN(IMIN) == 1 ) THEN
            TA = -MVSTNV( NU, 1 - ALPHA )/ATMN
         ELSE
            TA =  MVSTNV( NU, 1 - ALPHA/2 )/MAX( -ATMN, BTMN )
         END IF
         IF ( INFIN(IMAX) == 0 ) THEN
            TB =  MVSTNV( NU, 1 - ALPHA/M )/BTMX
         ELSE IF ( INFIN(IMAX) == 1 ) THEN
            TB = -MVSTNV( NU, 1 - ALPHA/M )/ATMX
         ELSE
            TB =  MVSTNV( NU, 1-ALPHA/(2*M) )/MIN( -ATMX, BTMX )
         END IF
         IF ( IP > 0 ) PRINT '(8X, 2(F8.4,8X,F8.4))', TA, TB
!
!     Use Pegasus method and 2nd order bounds to find better TALPHA interval.
!
!
         DL = 0
         TBF = TB
         DF = 0
         DO I = 1, 2
            CALL LWRUPR( N, COVRNC, NU, M, LOWER, CONSTR, UPPER, INFIN, DL, &
                         TA, BD(2), BD(1) )
            PA = BD(I) - 1 + ALPHA  
            CALL LWRUPR( N, COVRNC, NU, M, LOWER, CONSTR, UPPER, INFIN, DL, &
                         TB, BD(2), BD(1) )
            PB = BD(I) - 1 + ALPHA  
            IF ( IP > 1 ) PRINT '(I8, 2(F8.4,8X,F8.4))', I, TA, TB, PA, PB
!
            DO 
               IF ( ABS( TB - TA ) < ABSEPS/2 ) EXIT
               TC = TB - PB*( TB - TA )/( PB - PA )
               IF ( ( TC - TA )*( TC - TB ) > 0 ) TC = ( TA + TB )/2
               CALL LWRUPR( N, COVRNC, NU, M, LOWER,CONSTR,UPPER,INFIN, DL, &
                            TC, BD(2), BD(1) )
               PC = BD(I) - 1 + ALPHA  
               IF ( IP > 1 ) PRINT '(I8, 6F8.4)', I, TA, TC, TB, PA, PC, PB
               IF ( PC*PB <= 0 ) THEN 
                  TA = TB
                  PA = PB
               ELSE IF ( PA*PC > 0 ) THEN
                  IF ( ABS(PB) < ABS(PA) ) THEN
                     TA = TB
                     PA = PB
                  END IF
               ELSE
                  PA = PA*PB/( PB + PC )
               END IF
               TB = TC
               PB = PC
            END DO
!
            DF = DF + ( ABS( ( PB - PA )/( TB - TA ) ) - DF )/I
            TD(I) = ( TA + TB )/2
            IF ( I == 1 ) THEN
               TA = TD(I)
               TB = TBF
            END IF
         END DO
         TA = TD(1)
         TB = TD(2)
         TC = ( TA + TB )/2 
!
!     Refine TALPHA, if necessary, using Newton iteration.
!
         RE = 0
         PC = 1
         DF = MIN( DF, 2E-1_STND )
         DO 
            CALL LWRUPR( N, COVRNC, NU, M, LOWER, CONSTR, UPPER, INFIN, DL,  &
                      TC, BD(1), BD(2) )
            PCB = MAX( ABS( BD(1) - 1 + ALPHA ), ABS( BD(2) - 1 + ALPHA ) )/2
            IF ( IP > 1 ) PRINT '( I8, 3F8.4, 8X, F10.6, 4X, F6.2, 2F9.6 )',  &
                NVS, TA, TC, TB, PC, DF, BD(1), BD(2)
            ERROR = MIN( ABS( TB - TA )/2, ABS(PC/DF), ABS(PCB/DF) ) 
            IF ( ERROR < ABSEPS ) EXIT
            CALL MVDSTD( N, COVRNC, NU, M, LOWER,CONSTR,UPPER,INFIN, TC,      &
                   MAXPTS-NVS, ABS(DF)*ABSEPS, RE, ER, VL, DF, NV, INFORM )
            NVS = NVS + NV
            IF ( INFORM == 1 ) EXIT
            PC = VL - 1 + ALPHA
            IF ( PC > 0 ) THEN
               TB = ( TC + TB )/2
            ELSE
               TA = ( TA + TC )/2
            END IF
            TC = MIN( MAX( TA, TC - PC/DF ), TB )
         END DO
         TALPHA = TC
      END IF
!
      END SUBROUTINE MVCRIT
!
      SUBROUTINE LWRUPR( N, COVRNC, NU, M, LOWER, CONSTR, UPPER, LF, DELTA,T, &
                         LWRBND, UPRBND )
!
!     Hunter-Worsley and Dawson-Sankoff Bounds
!      References:
!         Hsu, Jason C. (1996), Multiple Comparisons, 
!           Chapman and Hall, London, p. 230;
!         Dawson, D. and Sankoff, A. (1967), An Inequality for Probabilities,
!           Proc. AMS 18, pp. 504-507.
!
!  Parameters
!
!     N      INTEGER number of variables.
!     COVRNC REAL covariance matrix; only the diagonal and lower left 
!             triangle are used. COVRNC matrix must be positive semidefinite.
!     NU     INTEGER, number of degrees of freedom; if NU < 1, then
!             a multivariate normal computation is completed. 
!     M      INTEGER, number of linear constraints for integration region.
!     LOWER  REAL array of lower integration limits.
!     CONSTR REAL MxN constraint matrix; integration region is all X with
!                     T*LOWER < DELTA + CONSTR*X < T*UPPER.
!     UPPER  REAL, array of upper integration limits.
!     LF    INTEGER, array of integration limits flags:
!            if LF(I) < 0, Ith limits are (-infinity, infinity);
!            if LF(I) = 0, Ith limits are (-infinity, UPPER(I)];
!            if LF(I) = 1, Ith limits are [LOWER(I), infinity);
!            if LF(I) = 2, Ith limits are [LOWER(I), UPPER(I)].
!     DELTA  REAL array of non-centrality parameters.
!     T      REAL integration limit scale factor.
!     LWRBND REAL lower bound for probability.
!     UPRBND REAL upper bound for probability.
!
      INTEGER,                         INTENT(IN)  :: N, NU, M
      REAL(KIND=STND), DIMENSION(N,N), INTENT(IN)  :: COVRNC
      REAL(KIND=STND), DIMENSION(M),   INTENT(IN)  :: LOWER, UPPER, DELTA
      REAL(KIND=STND), DIMENSION(M,N), INTENT(IN)  :: CONSTR
      REAL(KIND=STND),                 INTENT(IN)  :: T
      INTEGER,         DIMENSION(M),   INTENT(IN)  :: LF
      REAL(KIND=STND),                 INTENT(OUT) :: LWRBND, UPRBND
!
      INTEGER                         :: I, J, K
      REAL(KIND=STND)                 :: SA, SB, PT
      REAL(KIND=STND), DIMENSION(M)   :: L, U, MXD
      REAL(KIND=STND), DIMENSION(M,M) :: CR
!
!     Initialize maximum spanning tree and compute upper bound
!
      CALL COVSCL( N, COVRNC, M, LOWER, CONSTR, UPPER, LF, DELTA, T, L, U, CR )
      SA = 0
      SB = 0
      DO I = 1, M
         SA = SA + MVSTDC( NU, L(I), U(I), LF(I) ) 
         MXD(I) = 0
      END DO
      LWRBND = 1 - SA
!
!     Compute maximum spanning tree and lower bound
!
      J = 1
      DO K = 1, M
         DO I = 1, M
            IF ( MXD(I) > MXD(J) ) J = I
         END DO
         LWRBND = LWRBND + MXD(J) 
         MXD(J) = -1
         DO I = 1, M
            IF ( MXD(I) >= 0 ) THEN
               PT = MVBVTC( NU, (/ L(I), L(J) /),  (/ U(I), U(J) /),  &
                                (/ LF(I), LF(J) /), CR(I,J) )
               SB = SB + PT
               IF ( PT > MXD(I) ) MXD(I) = PT
            END IF
         END DO
      END DO
      LWRBND = MAX( LWRBND, 0E0_STND )
      K = 1 + 2*SB/SA
      UPRBND = 1 - 2*( SA - SB/K )/( K + 1 )
      UPRBND = MIN( UPRBND, 1E0_STND )
!
      END SUBROUTINE LWRUPR
!
!
      SUBROUTINE LWUPRH( N, COVRNC, NU, M, LOWER, CONSTR, UPPER, LF, DELTA,T, &
                         LWRBND, UPRBND )
!
!     Hybrid Second-Third Order Bounds
!      References:
!         Bukszar, J. and Prekopa, A. (2001), Probability Bounds with Cherry
!           Trees, Math. Oper. Res. 26, pp. 174-192;
!         Tomescu, I. (1986), Hypertrees and Bonferroni Inequalities,
!           Journal of Combinatorial Theory 41, pp. 209-217.
!
!  Parameters
!
!     N      INTEGER number of variables.
!     COVRNC REAL covariance matrix; only the diagonal and lower left 
!             triangle are used. COVRNC matrix must be positive semidefinite.
!     NU     INTEGER, number of degrees of freedom; if NU < 1, then
!             a multivariate normal computation is completed. 
!     M      INTEGER, number of linear constraints for integration region.
!     LOWER  REAL array of lower integration limits.
!     CONSTR REAL MxN constraint matrix; integration region is all X with
!                     T*LOWER < DELTA + CONSTR*X < T*UPPER.
!     UPPER  REAL, array of upper integration limits.
!     LF  INTEGER, array of integration limits flags:
!            if LF(I) < 0, Ith limits are (-infinity, infinity);
!            if LF(I) = 0, Ith limits are (-infinity, UPPER(I)];
!            if LF(I) = 1, Ith limits are [LOWER(I), infinity);
!            if LF(I) = 2, Ith limits are [LOWER(I), UPPER(I)].
!     DELTA  REAL array of non-centrality parameters.
!     T      REAL integration limit scale factor.
!     LWRBND REAL lower bound for probability.
!     UPRBND REAL upper bound for probability.
!
      INTEGER,                         INTENT(IN)  :: N, NU, M
      REAL(KIND=STND), DIMENSION(N,N), INTENT(IN)  :: COVRNC
      REAL(KIND=STND), DIMENSION(M),   INTENT(IN)  :: LOWER, UPPER, DELTA
      REAL(KIND=STND), DIMENSION(M,N), INTENT(IN)  :: CONSTR
      REAL(KIND=STND),                 INTENT(IN)  :: T
      INTEGER,         DIMENSION(M),   INTENT(IN)  :: LF
      REAL(KIND=STND),                 INTENT(OUT) :: LWRBND, UPRBND
!
      INTEGER                         :: I, J, K, L, E, V
      REAL(KIND=STND)                 :: SA, SB, W
      INTEGER,         DIMENSION(M)   :: NE
      REAL(KIND=STND), DIMENSION(M)   :: LT, UT, MXD
      REAL(KIND=STND), DIMENSION(M,M) :: CR, WT
      REAL(KIND=STND), DIMENSION(3)   :: LW, UP, SG
      INTEGER,         DIMENSION(3)   :: INF
      INTEGER                         :: IM, JM, KM
!
!     Compute scaled limits and correlation matrix.
!
      CALL COVSCL( N, COVRNC, M, LOWER, CONSTR, UPPER, LF, DELTA, T, LT,UT,CR )
      SA = 0
      SB = 0
      DO I = 1, M
         SA = SA + MVSTDC( NU, LT(I), UT(I), LF(I) ) 
         DO J = 1, I-1
            LW(1:2) = (/ LT(I), LT(J) /)
            UP(1:2) = (/ UT(I), UT(J) /)
            INF(1:2) = (/ LF(I), LF(J) /)
            WT(I,J) = MVBVTC( NU, LW(1:2), UP(1:2), INF(1:2), CR(I,J) ) 
            WT(J,I) = WT(I,J) 
            SB = SB + WT(I,J)
         END DO
      END DO
!
!     Compute maximum spanning tree lower bound
!
      LWRBND = 1 - SA
      DO I = 1, M
         WT(I,I) = -1
         NE(I) = 1
         MXD(I) = WT(I,1)
      END DO
      J = 2
      L = 2
      DO K = 1, M - 1
         DO I = 2, M
            IF ( MXD(I) > MXD(J) ) J = I
         END DO
         LWRBND = LWRBND + MXD(J) 
         MXD(J) = -1
         IF ( WT( J, NE(J) ) > WT(L, NE(L) ) ) L = J
         DO I = 2, M
            IF ( MXD(I) >= 0 .AND. WT(I,J) > MXD(I) ) THEN
               MXD(I) = WT(I,J)
               NE(I) = J
            END IF
         END DO
      END DO
!
!     Improve lower bound using cherry tree.
!
      NE(1) = 0
      MXD(1) = -1
      MXD(L)     = 1
      MXD(NE(L)) = 1
      DO L = 3, M
         W = -1
         DO K = 1, M
            IF ( MXD(K) < 0 ) THEN
               DO J = 1, M
                  IF ( MXD(J) > 0 .AND. ( J==NE(K) .OR. K==NE(J)  ) ) THEN
                     DO I = 1, M
                        IF ( MXD(I)>0 .AND. ( I==NE(J) .OR. J==NE(I) ) ) THEN 
                           IF( WT(I,K) > W ) THEN
                              W = WT(I,K)
                              IM = I
                              JM = J
                              KM = K
                           END IF
                        END IF
                     END DO
                  END IF
               END DO
            END IF
         END DO
         I = IM
         J = JM
         K = KM
         MXD(K) = 1
         LW = (/ LT(I), LT(J), LT(K) /)
         UP = (/ UT(I), UT(J), UT(K) /)
         INF = (/ LF(I), LF(J), LF(K) /)
         SG = (/ CR(I,J), CR(I,K), CR(J,K) /)
         LWRBND = LWRBND + WT(I,K) - MVTVTC( NU, LW, UP, INF, SG )
      END DO
      LWRBND = MAX( LWRBND, 0E0_STND )
!
!    Compute Tomescu Bound
! 
      UPRBND = 1 - SA + SB
      DO I = 2, M
         DO J = 1, I-1
            IF ( J /= NE(I) .AND. I /= NE(J) ) THEN
               WT(I,J) = -1
               WT(J,I) = -1
            END IF
         END DO
      END DO
      DO L = 3, M
         W = -1
         DO K = 1, M
            E = 0
            DO I = 1, M
               IF ( WT(I,K) >= 0 ) THEN
                  E = E + 1
                  J = I
               END IF
            END DO
            IF ( E == 1 .AND. WT(J,K) > W ) THEN
               W = WT(J,K) 
               V = K
            END IF
         END DO
         K = V
         DO I = 1, M
            IF ( WT(I,K) >= 0 ) J = I
         END DO
         WT(K,J) = -1
         WT(J,K) = -1
         DO I = 2, M
            J = NE(I) 
            IF ( WT(I,J) >= 0 ) THEN 
               LW = (/ LT(I), LT(J), LT(K) /)
               UP = (/ UT(I), UT(J), UT(K) /)
               INF = (/ LF(I), LF(J), LF(K) /)
               SG = (/ CR(I,J), CR(I,K), CR(J,K) /)
               UPRBND = UPRBND - MVTVTC( NU, LW, UP, INF, SG ) 
            END IF
         END DO
      END DO
      UPRBND = MIN( UPRBND, 1E0_STND )
!      K = 1 + 2*SB/SA
!      UPRBND = MIN( 1 - 2*( SA - SB/K )/( K + 1 ), UPRBND )
!
      END SUBROUTINE LWUPRH
!
!
      SUBROUTINE LWUPRT( N, COVRNC, NU, M, LOWER, CONSTR, UPPER, LF, DELTA,T, &
                         LWRBND, UPRBND )
!
!     Third Order Bounds
!      References:
!         Prekopa, A., Stochastic Programming, Kluwer Academic Publishers;  
!
!  Parameters
!
!     N      INTEGER number of variables.
!     COVRNC REAL covariance matrix; only the diagonal and lower left 
!             triangle are used. COVRNC matrix must be positive semidefinite.
!     NU     INTEGER, number of degrees of freedom; if NU < 1, then
!             a multivariate normal computation is completed. 
!     M      INTEGER, number of linear constraints for integration region.
!     LOWER  REAL array of lower integration limits.
!     CONSTR REAL MxN constraint matrix; integration region is all X with
!                     T*LOWER < DELTA + CONSTR*X < T*UPPER.
!     UPPER  REAL, array of upper integration limits.
!     LF  INTEGER, array of integration limits flags:
!            if LF(I) < 0, Ith limits are (-infinity, infinity);
!            if LF(I) = 0, Ith limits are (-infinity, UPPER(I)];
!            if LF(I) = 1, Ith limits are [LOWER(I), infinity);
!            if LF(I) = 2, Ith limits are [LOWER(I), UPPER(I)].
!     DELTA  REAL array of non-centrality parameters.
!     T      REAL integration limit scale factor.
!     LWRBND REAL lower bound for probability.
!     UPRBND REAL upper bound for probability.
!
      INTEGER,                         INTENT(IN)  :: N, NU, M
      REAL(KIND=STND), DIMENSION(N,N), INTENT(IN)  :: COVRNC
      REAL(KIND=STND), DIMENSION(M),   INTENT(IN)  :: LOWER, UPPER, DELTA
      REAL(KIND=STND), DIMENSION(M,N), INTENT(IN)  :: CONSTR
      REAL(KIND=STND),                 INTENT(IN)  :: T
      INTEGER,         DIMENSION(M),   INTENT(IN)  :: LF
      REAL(KIND=STND),                 INTENT(OUT) :: LWRBND, UPRBND
!
      INTEGER                         :: I, J, K
      REAL(KIND=STND)                 :: SA, SB, SC
      REAL(KIND=STND), DIMENSION(M)   :: LT, UT
      REAL(KIND=STND), DIMENSION(M,M) :: CR
      REAL(KIND=STND), DIMENSION(3)   :: LW, UP, SG
      INTEGER,         DIMENSION(3)   :: INF
!
!     Compute scaled limits and correlation matrix.
!
      CALL COVSCL( N, COVRNC, M, LOWER, CONSTR, UPPER, LF, DELTA, T, LT,UT,CR )
      SA = 0
      SB = 0
      SC = 0
      DO I = 1, M
         SA = SA + MVSTDC( NU, LT(I), UT(I), LF(I) ) 
         DO J = 1, I-1
            LW(1:2) = (/ LT(I), LT(J) /)
            UP(1:2) = (/ UT(I), UT(J) /)
            INF(1:2) = (/ LF(I), LF(J) /)
            SB = SB + MVBVTC( NU, LW(1:2), UP(1:2), INF(1:2), CR(I,J) ) 
            DO K = 1, J-1
               LW(3) = LT(K)
               UP(3) = UT(K)
               INF(3) = LF(K)
               SG = (/ CR(I,J), CR(I,K), CR(J,K) /)
               SC = SC + MVTVTC( NU, LW, UP, INF, SG )
            END DO
         END DO
      END DO
      J = 2 + 3*SC/SB
      LWRBND = 1 - ( SA - 2*( (2*J-1)*SB - 3*SC )/( J*(J+1) ) )
      I = 1 + 2*( ( M - 2 )*SB - 3*SC )/( ( M - 1 )*SA - 2*SB )
      UPRBND = 1 - ( (I+2*M-1)*SA - 2*( (2*I+M-2)*SB - 3*SC )/I )/( M*(I+1) ) 
      LWRBND = MAX( LWRBND, 0E0_STND )
      UPRBND = MIN( UPRBND, 1E0_STND )
!
      END SUBROUTINE LWUPRT
!
!
      SUBROUTINE COVSCL( N, COVRNC, M, LOWER, CONSTR, UPPER, LF, DELTA, T, &
                                        AT, BT, CR )
!
!     Compute Correlation Matrix and Scaled Limits
!
!  Parameters
!
!     N      INTEGER number of variables.
!     COVRNC REAL NxN covariance matrix; only the diagonal and lower left 
!             triangle are used. COVRNC must be positive semidefinite.
!     M      INTEGER, number of linear constraints for integration region.
!     LOWER  REAL array of lower integration limits.
!     CONSTR REAL MxN constraint matrix; integration region is all X with
!                     T*LOWER < DELTA + CONSTR*X < T*UPPER.
!     UPPER  REAL, array of upper integration limits.
!     LF  INTEGER, array of integration limits flags:
!            if LF(I) = 0, Ith limits are (-infinity, UPPER(I)];
!            if LF(I) = 1, Ith limits are [LOWER(I), infinity);
!            if LF(I) = 2, Ith limits are [LOWER(I), UPPER(I)].
!     DELTA  REAL array of non-centrality parameters.
!     T      REAL integration limit scale factor.
!     AT      REAL array of scaled lower integration limits.
!     BT      REAL array of scaled lower integration limits.
!     CR REAL MxM correlation matrix.
!
      INTEGER,                         INTENT(IN)  :: N, M
      REAL(KIND=STND), DIMENSION(N,N), INTENT(IN)  :: COVRNC
      REAL(KIND=STND), DIMENSION(M),   INTENT(IN)  :: LOWER, UPPER, DELTA
      REAL(KIND=STND), DIMENSION(M,N), INTENT(IN)  :: CONSTR
      REAL(KIND=STND),                 INTENT(IN)  :: T
      INTEGER,         DIMENSION(M),   INTENT(IN)  :: LF
      REAL(KIND=STND), DIMENSION(M,M), INTENT(OUT) :: CR
      REAL(KIND=STND), DIMENSION(M),   INTENT(OUT) :: AT, BT
!
      INTEGER                         :: I, J
      REAL(KIND=STND)                 :: SIGMAR
      REAL(KIND=STND), DIMENSION(N)   :: COVTMP
!
      DO I = 1, M
         DO J = 1, N
            COVTMP(J) = SUM( COVRNC(J,1:J-1)*CONSTR(I,1:J-1) )                &
                      + SUM( COVRNC(J:N,  J)*CONSTR(I,J:N) )
         END DO
         DO J = I, M
            CR(J,I) = SUM( CONSTR(J,1:N)*COVTMP )
         END DO
         SIGMAR = 1/SQRT( CR(I,I) )
         IF ( LF(I) /= 0 ) AT(I) = SIGMAR*( T*LOWER(I) - DELTA(I) )
         IF ( LF(I) /= 1 ) BT(I) = SIGMAR*( T*UPPER(I) - DELTA(I) )
         CR(I,1:I) = SIGMAR*CR(I,1:I)
         CR(I:M,I) = SIGMAR*CR(I:M,I)
      END DO
      DO I = 1, M
         CR(I,I:M) = CR(I:M,I)
      END DO
!
      END SUBROUTINE COVSCL
!
      SUBROUTINE MVKBRV( N, MINVLS, MAXVLS, NF, FUNCTN, ABSEPS, RELEPS, &
                                        ABSERR, FINEST, INTVLS, INFORM )
!
!  Automatic Multidimensional Integration Subroutine
!               
!       Alan Genz
!       Department of Mathematics
!       Washington State University
!       Pullman, Wa 99164-3113
!       Email : alangenz@wsu.edu
!
!         Last Change: 7/01
!
!  MVKBRV computes an approximation to the integral
!
!      1  1     1
!     I  I ... I   F(X)  dx(N)...dx(2)dx(1)
!      0  0     0
!
!    F(X) is an NF-vector of real integrands.
!
!  MVKBRV uses randomized Korobov rules. 
!  The primary references are
!   "Randomization of Number Theoretic Methods for Multiple Integration"
!    R. Cranley and T.N.L. Patterson, SIAM J Numer Anal, 13, pp. 904-14,
!  and 
!   "Optimal Parameters for Multidimensional Integration", 
!    P. Keast, SIAM J Numer Anal, 10, pp.831-838.
!  If N > 100 the subroutine uses quasi-random Richtmeyer points for 
!  variables with indices > 100. A reference is
!   "Methods of Numerical Integration", P.J. Davis and P. Rabinowitz, 
!    Academic Press, 1984, pp. 482-483.
!
!   
!**************  Parameters ********************************************
!***** Input parameters
!  N       Number of variables, must exceed 1, but not exceed 100
!  MINVLS  Integer minimum number of function evaluations allowed.
!          MINVLS must not exceed MAXVLS.  If MINVLS < 0 then the
!          routine assumes a previous call has been made with 
!          the same integrand and continues that calculation.
!  MAXVLS  Integer maximum number of function evaluations allowed.
!  NF      Integer number of integrands, must exceed 1 and not exceed 5000
!  FUNCTN  user defined function to be integrated.
!          It must have parameters (N,Z), where Z is a real array
!          of dimension N.
!                                     
!  ABSEPS  Required absolute accuracy.
!  RELEPS  Required relative accuracy.
!***** Output parameters
!  MINVLS  Actual number of function evaluations used.
!  ABSERR  Estimated absolute accuracy for FINEST(1).
!  FINEST  Estimated NF-vector of values of integrals.
!  INFORM  INFORM = 0 for normal exit, when 
!                     ABSERR <= MAX(ABSEPS, RELEPS*ABS(FINEST(1)))
!                  and 
!                     INTVLS <= MAXCLS.
!          INFORM = 1 If MAXVLS was too small to obtain the required 
!          accuracy. In this case a value FINEST is returned with 
!          estimated absolute accuracy ABSERR.
!***********************************************************************
!
! Global Variables
!
      INTEGER,                             INTENT(IN) :: N, NF, MINVLS, MAXVLS 
      INTERFACE 
         FUNCTION FUNCTN( NF, X ) RESULT(VALUE)
            USE PRECISION_MODEL
            INTEGER,         INTENT(IN)                :: NF
            REAL(KIND=STND), INTENT(IN), DIMENSION(:)  :: X
            REAL(KIND=STND),             DIMENSION(NF) :: VALUE
         END FUNCTION FUNCTN
      END INTERFACE   
      REAL(KIND=STND),                    INTENT(IN)  :: ABSEPS, RELEPS
      REAL(KIND=STND),     DIMENSION(NF), INTENT(OUT) :: FINEST
      REAL(KIND=STND),                    INTENT(OUT) :: ABSERR
      INTEGER,                            INTENT(OUT) :: INTVLS, INFORM
!
! Local Variables  
!

      INTEGER,                   PARAMETER :: PL = 28, NM = 99, NMP = NM + 1
      INTEGER,                   PARAMETER :: FM = 5000, NMX = 1000
      INTEGER,                   PARAMETER :: MINSMP = 8
      REAL(KIND=STND),           PARAMETER :: ONE = 1
      INTEGER,                        SAVE :: NP, SAMPLS
      REAL(KIND=STND), DIMENSION(FM), SAVE :: VAREST
      INTEGER,      DIMENSION(PL,NM)       :: C
      INTEGER                              :: I, K, KMX
      REAL(KIND=STND)                      :: VARPRD
      REAL(KIND=STND),       DIMENSION(NF) :: DIFINT, FINVAL, VARSQR
      REAL(KIND=STND),        DIMENSION(N) :: VK
      INTEGER, DIMENSION(PL), PARAMETER :: P = (/  31,      47,     73,  113, &
          173,    263,    397,    593,    907,   1361,    2053,   3079, 4621, &
         6947,  10427,  15641,  23473,  35221,  52837,   79259, 118891,       &
       178349, 267523, 401287, 601942, 902933, 1354471, 2031713 /)
      INTEGER, DIMENSION(NM), PARAMETER :: C01 = (/     12,      9,      9, &
                12,     12,     12,     12,     12,     12,     12,     12, &
                12,      3,      3,      3,     12,      7,      7,     12, &
                12,     12,     12,     12,     12,     12,     12,     12, &
                 3,      3,      3,     12,      7,      7,     12,     12, &
                12,     12,     12,     12,     12,     12,     12,      3, &
                 3,      3,     12,      7,      7,     12,     12,     12, &
                12,     12,     12,     12,     12,     12,      3,      3, &
                 3,     12,      7,      7,     12,     12,     12,     12, &
                12,     12,     12,     12,      7,      3,      3,      3, &
                 7,      7,      7,      3,      3,      3,      3,      3, &
                 3,      3,      3,      3,      3,      3,      3,      3, &
                 3,      3,      3,      3,      3,      3,      3,      3 /)
      INTEGER, DIMENSION(NM), PARAMETER :: C02 = (/     13,     11,     11, &
                10,     15,     15,     15,     15,     22,     15,     15, &
                15,     15,      6,      6,      6,     15,     15,      9, &
                13,      2,      2,      2,     13,     11,     11,     10, &
                15,     15,     15,     15,     15,     15,     15,     15, &
                15,      6,      6,      6,     15,     15,      9,     13, &
                 2,      2,      2,     13,     11,     11,     10,     15, &
                15,     15,     15,     15,     15,     15,     15,     15, &
                 6,      6,      6,     15,     15,      9,     13,      2, &
                 2,      2,     13,     11,     11,     10,     10,     15, &
                15,     15,     15,     15,     15,     15,     15,      6, &
                 2,      3,      2,      3,      2,      2,      2,      2, &
                 2,      2,      2,      2,      2,      2,      2,      2 /)
      INTEGER, DIMENSION(NM), PARAMETER :: C03 = (/     27,     28,     10, &
                11,     11,     20,     11,     20,     28,     13,     13, &
                28,     13,     13,     13,     14,     14,     14,     14, &
                14,     14,     14,     14,     14,     14,     14,     14, &
                14,     14,     14,     14,     31,     31,      5,      5, &
                 5,     31,     13,     11,     11,     11,     11,     11, &
                11,     13,     13,     13,     13,     13,     13,     13, &
                14,     14,     14,     14,     14,     14,     14,     14, &
                14,     14,     14,     14,     14,     14,     14,     14, &
                31,     31,      5,      5,      5,     11,     13,     11, &
                11,     11,     11,     11,     11,     11,     13,     13, &
                11,     13,      5,      5,      5,      5,     14,     13, &
                 5,      5,      5,      5,      5,      5,      5,      5 /)
      INTEGER, DIMENSION(NM), PARAMETER :: C04 = (/     35,     27,     27, &
                36,     22,     29,     29,     20,     45,      5,      5, &
                 5,     21,     21,     21,     21,     21,     21,     21, &
                21,     21,     21,     21,     21,     21,     21,     21, &
                21,     29,     17,     17,     17,     17,     17,     17, &
                17,     17,     17,     17,     23,     23,     23,     23, &
                23,     23,     23,     23,     23,     23,     23,     23, &
                21,     27,      3,      3,      3,     24,     27,     27, &
                17,     29,     29,     29,     17,      5,      5,      5, &
                 5,     21,     21,     21,     21,     21,     21,     21, &
                21,     21,     21,     21,     21,     21,     21,     21, &
                21,     17,     17,     17,      6,     17,     17,      6, &
                 3,      6,      6,      3,      3,      3,      3,      3 /)
      INTEGER, DIMENSION(NM), PARAMETER :: C05 = (/     64,     76,     28, &
                28,     44,     44,     55,     31,     10,     52,     10, &
                52,     10,     10,     38,     38,     10,     52,     10, &
                10,     10,     49,     60,     49,     49,     49,     49, &
                49,     49,     49,     49,     49,     49,     38,     38, &
                31,      4,      4,     31,     64,      4,      4,      4, &
                64,     45,     45,     45,     45,     45,     45,     66, &
                66,     66,     66,     66,     66,     66,     66,     66, &
                66,     66,     66,     66,     66,     66,     66,     66, &
                66,     66,     11,     66,     66,     66,     66,     66, &
                66,     66,     66,     66,     45,     11,      7,      3, &
                 2,      2,      2,     27,      5,      3,      3,      5, &
                 5,      2,      2,      2,      2,      2,      2,      2 /)
      INTEGER, DIMENSION(NM), PARAMETER :: C06 = (/    111,     42,     54, &
               118,     20,     31,     31,     72,     17,     94,     14, &
                14,     11,     14,     14,     14,     94,     10,     10, &
                10,     10,     14,     14,     14,     14,     14,     14, &
                14,     11,     11,     11,      8,      8,      8,      8, &
                 8,      8,      8,     18,     18,     18,     18,     18, &
               113,     62,     62,     45,     45,    113,    113,    113, &
               113,    113,    113,    113,    113,    113,    113,    113, &
               113,    113,    113,    113,    113,    113,     63,     63, &
                53,     63,     67,     67,     67,     67,     67,     67, &
                67,     67,     67,     67,     67,     67,     67,     67, &
                67,     51,     51,     51,     51,     51,     12,     51, &
                12,     51,      5,      3,      3,      2,      2,      5 /)
      INTEGER, DIMENSION(NM), PARAMETER :: C07 = (/    163,    154,     83, &
                43,     82,     92,    150,     59,     76,     76,     47, &
                11,     11,    100,    131,    116,    116,    116,    116, &
               116,    116,    138,    138,    138,    138,    138,    138, &
               138,    138,    138,    101,    101,    101,    101,    101, &
               101,    101,    101,    101,    101,    101,    101,    101, &
               101,    101,    101,    101,    101,    101,    101,    101, &
               116,    116,    116,    116,    116,    116,    100,    100, &
               100,    100,    100,    138,    138,    138,    138,    138, &
               101,    101,    101,    101,    101,    101,    101,    101, &
               101,    101,    101,    101,    101,    101,    101,    101, &
               101,    101,    101,     38,     38,     38,     38,     38, &
                38,     38,     38,      3,      3,      3,      3,      3 /)
      INTEGER, DIMENSION(NM), PARAMETER :: C08 = (/    246,    189,    242, &
               250,    250,    250,    102,    250,     36,    118,    196, &
               118,    191,    215,     49,    121,     49,    121,     49, &
                49,     49,     49,    121,     49,     49,     49,     49, &
                49,    171,    171,    171,    171,    171,    171,    171, &
               171,    171,    171,    171,    171,    171,    171,    171, &
               171,    171,    171,    171,    171,    171,    171,    171, &
               171,    171,    171,    171,    171,    171,    171,    171, &
               171,    171,    171,    161,    161,    161,    161,    161, &
               161,    161,    161,     14,     14,     14,     14,     14, &
                14,     14,     14,     14,     14,     14,     14,     14, &
                14,     14,     14,     14,     10,     10,     10,     10, &
                10,     10,    103,     10,     10,     10,     10,      5 /)
      INTEGER, DIMENSION(NM), PARAMETER :: C09 = (/    347,    402,    322, &
               418,    215,    220,    339,    339,    339,    337,    218, &
               167,    315,    315,    315,    167,    167,    167,    167, &
               361,    201,    124,    124,    124,    124,    124,    124, &
               124,    124,    124,    124,    124,    231,    231,     90, &
                90,     90,     90,     90,     90,     90,     90,     90, &
                90,     90,     90,     90,     90,     48,     48,     48, &
                48,     90,     90,     90,     90,     90,     90,     90, &
                90,     90,     90,     90,     90,     90,     90,     90, &
                90,     90,     90,     90,     90,     90,     90,     90, &
               243,    243,    243,    243,    243,    243,    243,    243, &
               243,    243,    283,    283,    283,    283,    283,    283, &
               283,    283,    283,     16,    283,     16,    283,    283 /)
      INTEGER, DIMENSION(NM), PARAMETER :: C10 = (/    505,    220,    601, &
               644,    612,    160,    206,    206,    206,    422,    134, &
               518,    134,    134,    518,    652,    382,    206,    158, &
               441,    179,    441,     56,    559,    559,     56,     56, &
                56,     56,     56,     56,     56,     56,     56,     56, &
                56,     56,     56,     56,    101,    101,     56,    101, &
               101,    101,    101,    101,    101,    101,    101,    193, &
               193,    193,    193,    193,    193,    193,    101,    101, &
               101,    101,    101,    101,    101,    101,    101,    101, &
               101,    101,    101,    101,    101,    101,    101,    101, &
               101,    101,    101,    122,    122,    122,    122,    122, &
               122,    122,    122,    122,    122,    122,    122,    122, &
               122,    122,    122,    122,    101,    101,    101,    101 /)
      INTEGER, DIMENSION(NM), PARAMETER :: C11 = (/    794,    325,    960, &
               528,    247,    247,    338,    366,    847,    753,    753, &
               236,    334,    334,    461,    711,    652,    381,    381, &
               381,    652,    381,    381,    381,    381,    381,    381, &
               381,    226,    326,    326,    326,    326,    326,    326, &
               326,    126,    326,    326,    326,    326,    326,    326, &
               326,    326,    326,    326,    195,    195,     55,     55, &
                55,     55,     55,     55,     55,     55,     55,     55, &
                55,     55,     55,     55,     55,     55,     55,     55, &
                55,    195,    195,    195,    195,    195,    195,    195, &
               132,    132,    132,    132,    132,    132,    132,    132, &
               132,    132,    132,    387,    387,    387,    387,    387, &
               387,    387,    387,    387,    387,    387,    387,    387 /)
      INTEGER, DIMENSION(NM), PARAMETER :: C12 = (/   1189,    888,    259, &
              1082,    725,    811,    636,    965,    497,    497,   1490, &
              1490,    392,   1291,    508,    508,   1291,   1291,    508, &
              1291,    508,    508,    867,    867,    934,    867,    934, &
               867,    867,    867,    867,    867,    867,    867,   1284, &
              1284,   1284,   1284,   1284,   1284,   1284,   1284,   1284, &
               563,    563,    563,    563,   1010,   1010,   1010,    208, &
               838,    563,    563,    563,    759,    759,    564,    759, &
               759,    801,    801,    801,    801,    759,    759,    759, &
               759,    759,    563,    563,    563,    563,    563,    563, &
               563,    563,    226,    226,    226,    226,    226,    226, &
               226,    226,    226,    226,    226,    226,    226,    226, &
               226,    226,    226,    226,    226,    226,    226,    226 /)
      INTEGER, DIMENSION(NM), PARAMETER :: C13 = (/   1763,   1018,   1500, &
               432,   1332,   1159,    126,   2240,   1719,   1284,    878, &
              1983,    266,    266,    266,    266,    747,    747,    127, &
               127,   2074,    127,   2074,   1400,   1383,   1383,   1383, &
              1383,   1383,   1383,   1383,   1383,   1383,   1383,   1400, &
              1383,   1383,   1383,   1383,   1383,   1383,   1383,    507, &
              1073,   1073,   1073,   1073,   1990,   1990,   1990,   1990, &
              1990,    507,    507,    507,    507,    507,    507,    507, &
               507,    507,   1073,   1073,   1073,   1073,   1073,   1073, &
              1073,   1073,   1073,   1073,   1073,   1073,   1073,   1073, &
              1073,   1073,   1073,     22,     22,     22,     22,     22, &
                22,   1073,    452,    452,    452,    452,    452,    452, &
               318,    301,    301,    301,    301,     86,     86,     15 /)
      INTEGER, DIMENSION(NM), PARAMETER :: C14 = (/   2872,   3233,   1534, &
              2941,   2910,    393,   1796,    919,    446,    919,    919, &
              1117,    103,    103,    103,    103,    103,    103,    103, &
              2311,   3117,   1101,   3117,   3117,   1101,   1101,   1101, &
              1101,   1101,   2503,   2503,   2503,   2503,   2503,   2503, &
              2503,   2503,    429,    429,    429,    429,    429,    429, &
               429,   1702,   1702,   1702,    184,    184,    184,    184, &
               184,    105,    105,    105,    105,    105,    105,    105, &
               105,    105,    105,    105,    105,    105,    105,    105, &
               105,    105,    105,    105,    105,    105,    105,    105, &
               105,    105,    105,    105,    105,    105,    105,    105, &
               105,    105,    105,    784,    784,    784,    784,    784, &
               784,    784,    784,    784,    784,    784,    784,    784 /)
      INTEGER, DIMENSION(NM), PARAMETER :: C15 = (/   4309,   3758,   4034, &
              1963,    730,    642,   1502,   2246,   3834,   1511,   1102, &
              1102,   1522,   1522,   3427,   3427,   3928,    915,    915, &
              3818,   3818,   3818,   3818,   4782,   4782,   4782,   3818, &
              4782,   3818,   3818,   1327,   1327,   1327,   1327,   1327, &
              1327,   1327,   1387,   1387,   1387,   1387,   1387,   1387, &
              1387,   1387,   1387,   2339,   2339,   2339,   2339,   2339, &
              2339,   2339,   2339,   2339,   2339,   2339,   2339,   2339, &
              3148,   3148,   3148,   3148,   3148,   3148,   3148,   3148, &
              3148,   3148,   3148,   3148,   3148,   3148,   3148,   3148, &
              3148,   3148,   1776,   1776,   1776,   3354,   3354,   3354, &
               925,   3354,   3354,    925,    925,    925,    925,    925, &
              2133,   2133,   2133,   2133,   2133,   2133,   2133,   2133 /)
      INTEGER, DIMENSION(NM), PARAMETER :: C16 = (/   6610,   6977,   1686, &
              3819,   2314,   5647,   3953,   3614,   5115,    423,    423, &
              5408,   7652,    423,    423,    487,   6227,   2660,   6227, &
              1221,   3811,    197,   4367,   4367,   1281,   1221,    351, &
               351,    351,   1984,   1984,   3995,   2999,   2999,   2999, &
              2999,   2999,   2999,   2063,   2063,   2063,   2063,   1644, &
              2063,   2077,   2512,   2512,   2512,   2077,   2077,   2077, &
              2077,    754,    754,    754,    754,    754,    754,    754, &
               754,    754,    754,    754,    754,    754,    754,    754, &
               754,    754,    754,    754,   1097,   1097,    754,    754, &
               754,    754,    248,    754,   1097,   1097,   1097,   1097, &
               222,    222,    222,    222,    754,   1982,   1982,   1982, &
              1982,   1982,   1982,   1982,   1982,   1982,   1982,   1982 /)  
      INTEGER, DIMENSION(NM), PARAMETER :: C17 = (/   9861,   3647,   4073, &
              2535,   3430,   9865,   2830,   9328,   4320,   5913,  10365, &
              8272,   3706,   6186,   7806,   7806,   7806,   8610,   2563, &
             11558,  11558,   9421,   1181,   9421,   1181,   1181,   1181, &
              9421,   1181,   1181,  10574,  10574,   3534,   3534,   3534, &
              3534,   3534,   2898,   2898,   2898,   3450,   2141,   2141, &
              2141,   2141,   2141,   2141,   2141,   7055,   7055,   7055, &
              7055,   7055,   7055,   7055,   7055,   7055,   7055,   7055, &
              7055,   7055,   7055,   7055,   2831,   8204,   8204,   8204, &
              8204,   8204,   8204,   8204,   8204,   8204,   8204,   8204, &
              8204,   8204,   8204,   8204,   8204,   8204,   8204,   8204, &
              8204,   8204,   8204,   8204,   8204,   4688,   4688,   4688, &
              2831,   2831,   2831,   2831,   2831,   2831,   2831,   2831 /)
      INTEGER, DIMENSION(NM), PARAMETER :: C18 = (/  10327,   7582,   7124, &
              8214,   9600,  10271,  10193,  10800,   9086,   2365,   4409, &
             13812,   5661,   9344,   9344,  10362,   9344,   9344,   8585, &
             11114,  13080,  13080,  13080,   6949,   3436,   3436,   3436, &
             13213,   6130,   6130,   8159,   8159,  11595,   8159,   3436, &
              7096,   7096,   7096,   7096,   7096,   7096,   7096,   7096, &
              7096,   7096,   7096,   7096,   7096,   7096,   7096,   7096, &
              7096,   7096,   4377,   7096,   4377,   4377,   4377,   4377, &
              4377,   5410,   5410,   4377,   4377,   4377,   4377,   4377, &
              4377,   4377,   4377,   4377,   4377,   4377,   4377,   4377, &
              4377,   4377,   4377,   4377,   4377,   4377,   4377,   4377, &
              4377,   4377,   4377,   4377,   4377,   4377,   4377,   4377, &
              4377,   4377,   4377,    440,    440,   1199,   1199,   1199 /)
      INTEGER, DIMENSION(NM), PARAMETER :: C19 = (/  19540,  19926,  11582, &
             11113,  24585,   8726,  17218,    419,   4918,   4918,   4918, &
             15701,  17710,   4037,   4037,  15808,  11401,  19398,  25950, &
             25950,   4454,  24987,  11719,   8697,   1452,   1452,   1452, &
              1452,   1452,   8697,   8697,   6436,  21475,   6436,  22913, &
              6434,  18497,  11089,  11089,  11089,  11089,   3036,   3036, &
             14208,  14208,  14208,  14208,  12906,  12906,  12906,  12906, &
             12906,  12906,  12906,  12906,   7614,   7614,   7614,   7614, &
              5021,   5021,   5021,   5021,   5021,   5021,  10145,  10145, &
             10145,  10145,  10145,  10145,  10145,  10145,  10145,  10145, &
             10145,  10145,  10145,  10145,  10145,  10145,  10145,  10145, &
             10145,  10145,  10145,  10145,  10145,  10145,   4544,   4544, &
              4544,   4544,   4544,   4544,   8394,   8394,   8394,   8394 /)
      INTEGER, DIMENSION(NM), PARAMETER :: C20 = (/  34566,   9579,  12654, &
             26856,  37873,  38806,  29501,  17271,   3663,  10763,  18955, &
              1298,  26560,  17132,  17132,   4753,   4753,   8713,  18624, &
             13082,   6791,   1122,  19363,  34695,  18770,  18770,  18770, &
             18770,  15628,  18770,  18770,  18770,  18770,  33766,  20837, &
             20837,  20837,  20837,  20837,  20837,   6545,   6545,   6545, &
              6545,   6545,  12138,  12138,  12138,  12138,  12138,  12138, &
             12138,  12138,  12138,  12138,  12138,  12138,  12138,  12138, &
             30483,  30483,  30483,  30483,  30483,  12138,  12138,  12138, &
             12138,  12138,  12138,  12138,  12138,  12138,  12138,  12138, &
             12138,  12138,  12138,  12138,  12138,  12138,  12138,  12138, &
              9305,  11107,  11107,  11107,  11107,  11107,  11107,  11107, &
             11107,  11107,  11107,  11107,  11107,  11107,   9305,   9305 /)
      INTEGER, DIMENSION(NM), PARAMETER :: C21 = (/  31929,  49367,  10982, &
              3527,  27066,  13226,  56010,  18911,  40574,  20767,  20767, &
              9686,  47603,  47603,  11736,  11736,  41601,  12888,  32948, &
             30801,  44243,  53351,  53351,  16016,  35086,  35086,  32581, &
              2464,   2464,  49554,   2464,   2464,  49554,  49554,   2464, &
                81,  27260,  10681,   2185,   2185,   2185,   2185,   2185, &
              2185,   2185,  18086,  18086,  18086,  18086,  18086,  17631, &
             17631,  18086,  18086,  18086,  37335,  37774,  37774,  37774, &
             26401,  26401,  26401,  26401,  26401,  26401,  26401,  26401, &
             26401,  26401,  26401,  26401,  26401,  12982,  40398,  40398, &
             40398,  40398,  40398,  40398,   3518,   3518,   3518,  37799, &
             37799,  37799,  37799,  37799,  37799,  37799,  37799,  37799, &
              4721,   4721,   4721,   4721,   7067,   7067,   7067,   7067 /)
      INTEGER, DIMENSION(NM), PARAMETER :: C22 = (/  40701,  69087,  77576, &
             64590,  39397,  33179,  10858,  38935,  43129,  35468,  35468, &
              5279,  61518,  61518,  27945,  70975,  70975,  86478,  86478, &
             20514,  20514,  73178,  73178,  43098,  43098,   4701,  59979, &
             59979,  58556,  69916,  15170,  15170,   4832,   4832,  43064, &
             71685,   4832,  15170,  15170,  15170,  27679,  27679,  27679, &
             60826,  60826,   6187,   6187,   4264,   4264,   4264,   4264, &
              4264,  45567,  32269,  32269,  32269,  32269,  62060,  62060, &
             62060,  62060,  62060,  62060,  62060,  62060,  62060,   1803, &
              1803,   1803,   1803,   1803,   1803,   1803,   1803,   1803, &
              1803,   1803,   1803,   1803,  51108,  51108,  51108,  51108, &
             51108,  51108,  51108,  51108,  51108,  51108,  51108,  51108, &
             55315,  55315,  54140,  54140,  54140,  54140,  54140,  13134 /)
      INTEGER, DIMENSION(NM), PARAMETER :: C23 = (/ 103650, 125480,  59978, &
             46875,  77172,  83021, 126904,  14541,  56299,  43636,  11655, &
             52680,  88549,  29804, 101894, 113675,  48040, 113675,  34987, &
             48308,  97926,   5475,  49449,   6850,  62545,  62545,   9440, &
             33242,   9440,  33242,   9440,  33242,   9440,  62850,   9440, &
              9440,   9440,  90308,  90308,  90308,  47904,  47904,  47904, &
             47904,  47904,  47904,  47904,  47904,  47904,  41143,  41143, &
             41143,  41143,  41143,  41143,  41143,  36114,  36114,  36114, &
             36114,  36114,  24997,  65162,  65162,  65162,  65162,  65162, &
             65162,  65162,  65162,  65162,  65162,  65162,  65162,  65162, &
             65162,  47650,  47650,  47650,  47650,  47650,  47650,  47650, &
             40586,  40586,  40586,  40586,  40586,  40586,  40586,  38725, &
             38725,  38725,  38725,  88329,  88329,  88329,  88329,  88329 /)
      INTEGER, DIMENSION(NM), PARAMETER :: C24 = (/ 165843,  90647,  59925, &
            189541,  67647,  74795,  68365, 167485, 143918,  74912, 167289, &
             75517,   8148, 172106, 126159,  35867,  35867,  35867, 121694, &
             52171,  95354, 113969, 113969,  76304, 123709, 123709, 144615, &
            123709,  64958,  64958,  32377, 193002, 193002,  25023,  40017, &
            141605, 189165, 189165, 141605, 189165, 189165, 141605, 141605, &
            141605, 189165, 127047, 127047, 127047, 127047, 127047, 127047, &
            127047, 127047, 127047, 127047, 127047, 127047, 127047, 127047, &
            127047, 127047, 127047, 127047, 127047, 127047, 127785, 127785, &
            127785, 127785, 127785, 127785, 127785, 127785, 127785, 127785, &
             80822,  80822,  80822,  80822,  80822,  80822, 131661, 131661, &
            131661, 131661, 131661, 131661, 131661, 131661, 131661, 131661, &
            131661, 131661, 131661, 131661, 131661, 131661,   7114, 131661 /)
      INTEGER, DIMENSION(NM), PARAMETER :: C25 = (/ 130365, 236711, 110235, &
            125699,  56483,  93735, 234469,  60549,   1291,  93937, 245291, &
            196061, 258647, 162489, 176631, 204895,  73353, 172319,  28881, &
            136787, 122081, 122081, 275993,  64673, 211587, 211587, 211587, &
            282859, 282859, 211587, 242821, 256865, 256865, 256865, 122203, &
            291915, 122203, 291915, 291915, 122203,  25639,  25639, 291803, &
            245397, 284047, 245397, 245397, 245397, 245397, 245397, 245397, &
            245397,  94241,  66575,  66575, 217673, 217673, 217673, 217673, &
            217673, 217673, 217673, 217673, 217673, 217673, 217673, 217673, &
            217673, 217673, 217673, 217673, 217673, 217673, 217673, 210249, &
            210249, 210249, 210249, 210249, 210249, 210249, 210249, 210249, &
            210249,  94453,  94453,  94453,  94453,  94453,  94453,  94453, &
             94453,  94453,  94453,  94453,  94453,  94453,  94453,  94453 /)
      INTEGER, DIMENSION(NM), PARAMETER :: C26 = (/ 333459, 375354, 102417, &
            383544, 292630,  41147, 374614,  48032, 435453, 281493, 358168, &
            114121, 346892, 238990, 317313, 164158,  35497,  70530,  70530, &
            434839,  24754,  24754,  24754, 393656, 118711, 118711, 148227, &
            271087, 355831,  91034, 417029, 417029,  91034,  91034, 417029, &
             91034, 299843, 299843, 413548, 413548, 308300, 413548, 413548, &
            413548, 308300, 308300, 308300, 413548, 308300, 308300, 308300, &
            308300, 308300,  15311,  15311,  15311,  15311, 176255, 176255, &
             23613,  23613,  23613,  23613,  23613,  23613, 172210, 204328, &
            204328, 204328, 204328, 121626, 121626, 121626, 121626, 121626, &
            200187, 200187, 200187, 200187, 200187, 121551, 121551, 248492, &
            248492, 248492, 248492, 248492, 248492, 248492, 248492, 248492, &
            248492, 248492, 248492,  13942,  13942,  13942,  13942,  13942 /)
      INTEGER, DIMENSION(NM), PARAMETER :: C27 = (/ 500884, 566009, 399251, &
            652979, 355008, 430235, 328722, 670680, 405585, 405585, 424646, &
            670180, 670180, 641587, 215580,  59048, 633320,  81010,  20789, &
            389250, 389250, 638764, 638764, 389250, 389250, 398094,  80846, &
            147776, 147776, 296177, 398094, 398094, 147776, 147776, 396313, &
            578233, 578233, 578233,  19482, 620706, 187095, 620706, 187095, &
            126467, 241663, 241663, 241663, 241663, 241663, 241663, 241663, &
            241663, 241663, 241663, 241663, 241663, 321632,  23210,  23210, &
            394484, 394484, 394484,  78101,  78101,  78101, 542095, 542095, &
            542095, 542095, 542095, 542095, 542095, 542095, 542095, 542095, &
            542095, 542095, 542095, 542095, 542095, 542095, 542095, 542095, &
            542095, 277743, 277743, 277743, 457259, 457259, 457259, 457259, &
            457259, 457259, 457259, 457259, 457259, 457259, 457259, 457259 /)
      INTEGER, DIMENSION(NM), PARAMETER :: C28 = (/ 858339, 918142, 501970, &
            234813, 460565,  31996, 753018, 256150, 199809, 993599, 245149, &
            794183, 121349, 150619, 376952, 809123, 809123, 804319,  67352, &
            969594, 434796, 969594, 804319, 391368, 761041, 754049, 466264, &
            754049, 754049, 466264, 754049, 754049, 282852, 429907, 390017, &
            276645, 994856, 250142, 144595, 907454, 689648, 687580, 687580, &
            687580, 687580, 978368, 687580, 552742, 105195, 942843, 768249, &
            307142, 307142, 307142, 307142, 880619, 880619, 880619, 880619, &
            880619, 880619, 880619, 117185, 117185, 117185, 117185, 117185, &
            117185, 117185, 117185, 117185, 117185, 117185,  60731,  60731, &
             60731,  60731,  60731,  60731,  60731,  60731,  60731,  60731, &
             60731, 178309, 178309, 178309, 178309,  74373,  74373,  74373, &
             74373,  74373,  74373,  74373,  74373, 214965, 214965, 214965 /)
      INTEGER, DIMENSION(NMX), PARAMETER ::   PRIME = (/            &
           2,    3,    5,    7,   11,   13,   17,   19,   23,   29, &
          31,   37,   41,   43,   47,   53,   59,   61,   67,   71, &
          73,   79,   83,   89,   97,  101,  103,  107,  109,  113, &
         127,  131,  137,  139,  149,  151,  157,  163,  167,  173, &
         179,  181,  191,  193,  197,  199,  211,  223,  227,  229, &
         233,  239,  241,  251,  257,  263,  269,  271,  277,  281, &
         283,  293,  307,  311,  313,  317,  331,  337,  347,  349, &
         353,  359,  367,  373,  379,  383,  389,  397,  401,  409, &
         419,  421,  431,  433,  439,  443,  449,  457,  461,  463, &
         467,  479,  487,  491,  499,  503,  509,  521,  523,  541, &
         547,  557,  563,  569,  571,  577,  587,  593,  599,  601, &
         607,  613,  617,  619,  631,  641,  643,  647,  653,  659, &
         661,  673,  677,  683,  691,  701,  709,  719,  727,  733, &
         739,  743,  751,  757,  761,  769,  773,  787,  797,  809, &
         811,  821,  823,  827,  829,  839,  853,  857,  859,  863, &
         877,  881,  883,  887,  907,  911,  919,  929,  937,  941, &
         947,  953,  967,  971,  977,  983,  991,  997, 1009, 1013, &
        1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069, &
        1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151, &
        1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223, &
        1229, 1231, 1237, 1249, 1259, 1277, 1279, 1283, 1289, 1291, &
        1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373, &
        1381, 1399, 1409, 1423, 1427, 1429, 1433, 1439, 1447, 1451, &
        1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511, &
        1523, 1531, 1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583, &
        1597, 1601, 1607, 1609, 1613, 1619, 1621, 1627, 1637, 1657, &
        1663, 1667, 1669, 1693, 1697, 1699, 1709, 1721, 1723, 1733, &
        1741, 1747, 1753, 1759, 1777, 1783, 1787, 1789, 1801, 1811, &
        1823, 1831, 1847, 1861, 1867, 1871, 1873, 1877, 1879, 1889, &
        1901, 1907, 1913, 1931, 1933, 1949, 1951, 1973, 1979, 1987, &
        1993, 1997, 1999, 2003, 2011, 2017, 2027, 2029, 2039, 2053, &
        2063, 2069, 2081, 2083, 2087, 2089, 2099, 2111, 2113, 2129, &
        2131, 2137, 2141, 2143, 2153, 2161, 2179, 2203, 2207, 2213, &
        2221, 2237, 2239, 2243, 2251, 2267, 2269, 2273, 2281, 2287, &
        2293, 2297, 2309, 2311, 2333, 2339, 2341, 2347, 2351, 2357, &
        2371, 2377, 2381, 2383, 2389, 2393, 2399, 2411, 2417, 2423, &
        2437, 2441, 2447, 2459, 2467, 2473, 2477, 2503, 2521, 2531, &
        2539, 2543, 2549, 2551, 2557, 2579, 2591, 2593, 2609, 2617, &
        2621, 2633, 2647, 2657, 2659, 2663, 2671, 2677, 2683, 2687, &
        2689, 2693, 2699, 2707, 2711, 2713, 2719, 2729, 2731, 2741, &
        2749, 2753, 2767, 2777, 2789, 2791, 2797, 2801, 2803, 2819, &
        2833, 2837, 2843, 2851, 2857, 2861, 2879, 2887, 2897, 2903, &
        2909, 2917, 2927, 2939, 2953, 2957, 2963, 2969, 2971, 2999, &
        3001, 3011, 3019, 3023, 3037, 3041, 3049, 3061, 3067, 3079, &
        3083, 3089, 3109, 3119, 3121, 3137, 3163, 3167, 3169, 3181, &
        3187, 3191, 3203, 3209, 3217, 3221, 3229, 3251, 3253, 3257, &
        3259, 3271, 3299, 3301, 3307, 3313, 3319, 3323, 3329, 3331, &
        3343, 3347, 3359, 3361, 3371, 3373, 3389, 3391, 3407, 3413, &
        3433, 3449, 3457, 3461, 3463, 3467, 3469, 3491, 3499, 3511, &
        3517, 3527, 3529, 3533, 3539, 3541, 3547, 3557, 3559, 3571, &
        3581, 3583, 3593, 3607, 3613, 3617, 3623, 3631, 3637, 3643, &
        3659, 3671, 3673, 3677, 3691, 3697, 3701, 3709, 3719, 3727, &
        3733, 3739, 3761, 3767, 3769, 3779, 3793, 3797, 3803, 3821, &
        3823, 3833, 3847, 3851, 3853, 3863, 3877, 3881, 3889, 3907, &
        3911, 3917, 3919, 3923, 3929, 3931, 3943, 3947, 3967, 3989, &
        4001, 4003, 4007, 4013, 4019, 4021, 4027, 4049, 4051, 4057, &
        4073, 4079, 4091, 4093, 4099, 4111, 4127, 4129, 4133, 4139, &
        4153, 4157, 4159, 4177, 4201, 4211, 4217, 4219, 4229, 4231, &
        4241, 4243, 4253, 4259, 4261, 4271, 4273, 4283, 4289, 4297, &
        4327, 4337, 4339, 4349, 4357, 4363, 4373, 4391, 4397, 4409, &
        4421, 4423, 4441, 4447, 4451, 4457, 4463, 4481, 4483, 4493, &
        4507, 4513, 4517, 4519, 4523, 4547, 4549, 4561, 4567, 4583, &
        4591, 4597, 4603, 4621, 4637, 4639, 4643, 4649, 4651, 4657, &
        4663, 4673, 4679, 4691, 4703, 4721, 4723, 4729, 4733, 4751, &
        4759, 4783, 4787, 4789, 4793, 4799, 4801, 4813, 4817, 4831, &
        4861, 4871, 4877, 4889, 4903, 4909, 4919, 4931, 4933, 4937, &
        4943, 4951, 4957, 4967, 4969, 4973, 4987, 4993, 4999, 5003, &
        5009, 5011, 5021, 5023, 5039, 5051, 5059, 5077, 5081, 5087, &
        5099, 5101, 5107, 5113, 5119, 5147, 5153, 5167, 5171, 5179, &
        5189, 5197, 5209, 5227, 5231, 5233, 5237, 5261, 5273, 5279, &
        5281, 5297, 5303, 5309, 5323, 5333, 5347, 5351, 5381, 5387, &
        5393, 5399, 5407, 5413, 5417, 5419, 5431, 5437, 5441, 5443, &
        5449, 5471, 5477, 5479, 5483, 5501, 5503, 5507, 5519, 5521, &
        5527, 5531, 5557, 5563, 5569, 5573, 5581, 5591, 5623, 5639, &
        5641, 5647, 5651, 5653, 5657, 5659, 5669, 5683, 5689, 5693, &
        5701, 5711, 5717, 5737, 5741, 5743, 5749, 5779, 5783, 5791, &
        5801, 5807, 5813, 5821, 5827, 5839, 5843, 5849, 5851, 5857, &
        5861, 5867, 5869, 5879, 5881, 5897, 5903, 5923, 5927, 5939, &
        5953, 5981, 5987, 6007, 6011, 6029, 6037, 6043, 6047, 6053, &
        6067, 6073, 6079, 6089, 6091, 6101, 6113, 6121, 6131, 6133, &
        6143, 6151, 6163, 6173, 6197, 6199, 6203, 6211, 6217, 6221, &
        6229, 6247, 6257, 6263, 6269, 6271, 6277, 6287, 6299, 6301, &
        6311, 6317, 6323, 6329, 6337, 6343, 6353, 6359, 6361, 6367, &
        6373, 6379, 6389, 6397, 6421, 6427, 6449, 6451, 6469, 6473, &
        6481, 6491, 6521, 6529, 6547, 6551, 6553, 6563, 6569, 6571, &
        6577, 6581, 6599, 6607, 6619, 6637, 6653, 6659, 6661, 6673, &
        6679, 6689, 6691, 6701, 6703, 6709, 6719, 6733, 6737, 6761, &
        6763, 6779, 6781, 6791, 6793, 6803, 6823, 6827, 6829, 6833, &
        6841, 6857, 6863, 6869, 6871, 6883, 6899, 6907, 6911, 6917, &
        6947, 6949, 6959, 6961, 6967, 6971, 6977, 6983, 6991, 6997, &
        7001, 7013, 7019, 7027, 7039, 7043, 7057, 7069, 7079, 7103, &
        7109, 7121, 7127, 7129, 7151, 7159, 7177, 7187, 7193, 7207, &
        7211, 7213, 7219, 7229, 7237, 7243, 7247, 7253, 7283, 7297, &
        7307, 7309, 7321, 7331, 7333, 7349, 7351, 7369, 7393, 7411, &
        7417, 7433, 7451, 7457, 7459, 7477, 7481, 7487, 7489, 7499, &
        7507, 7517, 7523, 7529, 7537, 7541, 7547, 7549, 7559, 7561, &
        7573, 7577, 7583, 7589, 7591, 7603, 7607, 7621, 7639, 7643, &
        7649, 7669, 7673, 7681, 7687, 7691, 7699, 7703, 7717, 7723, &
        7727, 7741, 7753, 7757, 7759, 7789, 7793, 7817, 7823, 7829, &
        7841, 7853, 7867, 7873, 7877, 7879, 7883, 7901, 7907, 7919 /)
!
      C = RESHAPE( (/  C01, C02, C03, C04, C05, C06, C07, C08, C09, C10,     &
                       C11, C12, C13, C14, C15, C16, C17, C18, C19, C20,     &
                       C21, C22, C23, C24, C25, C26, C27, C28 /),            & 
                       (/ PL, NM /),  ORDER = (/ 2, 1 /) )
!
      INFORM = 1
      INTVLS = 0
      IF ( MINVLS >= 0 ) THEN
         FINEST = 0
         VAREST = 0
         SAMPLS = MINSMP 
         DO I = MIN( N, 10 ), PL
            NP = I
            IF ( MINVLS < 2*SAMPLS*P(I) ) EXIT
         END DO
         IF ( MINVLS >= 2*SAMPLS*P(PL) ) SAMPLS = MINVLS/( 2*P(PL) ) 
      ENDIF
      DO
         VK(1) = ONE/P(NP)
         DO I = 2, N
            IF ( I <= NMP ) THEN
               VK(I) = MODULO( C(NP,N-1)*VK(I-1), ONE )
            ELSE
               VK(I) = MODULO( SQRT( REAL( PRIME(I-NMP), STND ) ), ONE )
            END IF
         END DO
!
         FINVAL = 0
         VARSQR = 0
         DO I = 1, SAMPLS
            DIFINT = ( MVKRSV( N, P(NP), VK, NF, FUNCTN ) - FINVAL )/I
            FINVAL = FINVAL + DIFINT
            VARSQR = ( I - 2 )*VARSQR/I + DIFINT**2
         END DO
         INTVLS = INTVLS + 2*SAMPLS*P(NP)
         KMX = 1
         DO K = 1, NF
            VARPRD = VAREST(K)*VARSQR(K)
            FINEST(K) = FINEST(K) + ( FINVAL(K) - FINEST(K) )/( 1 + VARPRD )
            IF ( VARSQR(K) > 0 ) VAREST(K) = ( 1 + VARPRD )/VARSQR(K)
         END DO
         ABSERR = 7*SQRT( VARSQR(KMX)/( 1 + VARPRD ) )/2 
         IF ( ABSERR > MAX( ABSEPS, ABS(FINEST(KMX))*RELEPS ) ) THEN
            IF ( NP < PL ) THEN
               NP = NP + 1
            ELSE
               SAMPLS = MIN( 3*SAMPLS/2, ( MAXVLS - INTVLS )/( 2*P(NP) ) ) 
               SAMPLS = MAX( MINSMP, SAMPLS )
            ENDIF
            IF ( INTVLS + 2*SAMPLS*P(NP) > MAXVLS ) EXIT
         ELSE
            INFORM = 0
            EXIT
         ENDIF
      END DO
!
      END SUBROUTINE MVKBRV
!
      FUNCTION MVKRSV( N, PRIME, VK, NF, FUNCTN ) RESULT(VALUES)
!
! Global Variables
!
      INTEGER,                        INTENT(IN)      :: N, NF, PRIME 
      REAL(KIND=STND),  DIMENSION(:), INTENT(INOUT)   :: VK
      INTERFACE 
         FUNCTION FUNCTN( NF, X ) RESULT(VALUE)
            USE PRECISION_MODEL
            INTEGER,                       INTENT(IN) :: NF
            REAL(KIND=STND), DIMENSION(:), INTENT(IN) :: X
            REAL(KIND=STND),            DIMENSION(NF) :: VALUE
         END FUNCTION FUNCTN
      END INTERFACE   
      REAL(KIND=STND),                  DIMENSION(NF) :: VALUES
!
! Local Variables  
!
      INTEGER                       :: K
      REAL(KIND=STND), DIMENSION(N) :: X, R
!
      VALUES = 0
      DO K = 1, N
         R(K) = UNIFRM()
      END DO
      DO K = 1, PRIME
         R = R + VK
         WHERE ( R > 1 ) R = R - 1
         X = ABS ( 2*R - 1 )
         VALUES = VALUES + ( ( FUNCTN(NF,X) + FUNCTN(NF,1-X) )/2 - VALUES )/K
      END DO
!
      END FUNCTION MVKRSV
!
      FUNCTION UNIFRM() RESULT(UNI) 
!
!     Uniform (0,1) random number generator
!
!     Reference:
!     L'Ecuyer, Pierre (1996), 
!     "Combined Multiple Recursive Random Number Generators"
!     Operations Research 44, pp. 816-822.
!
      REAL(KIND=STND) :: UNI
!
      INTEGER         :: Z, H, P12, P13, P21, P23
!
!     Some eight digit primes for seeds
!
      INTEGER,   SAVE :: X10 = 15485857, X11 = 17329489, X12 = 36312197, & 
                         X20 = 55911127, X21 = 75906931, X22 = 96210113       
      REAL(KIND=STND), PARAMETER :: INVMP1 = 4.656612873077392578125E-10_STND 
!                                   INVMP1 = 1/(M1+1)
      INTEGER, PARAMETER :: M1 = 2147483647, M2 = 2145483479 
      INTEGER, PARAMETER :: A12 =   63308, Q12 = 33921, R12 = 12979 
      INTEGER, PARAMETER :: A13 = -183326, Q13 = 11714, R13 =  2883 
      INTEGER, PARAMETER :: A21 =   86098, Q21 = 24919, R21 =  7417 
      INTEGER, PARAMETER :: A23 = -539608, Q23 =  3976, R23 =  2071 
!
!
!     Component 1
!
      H = X10/Q13
      P13 = -A13*( X10 - H*Q13 ) - H*R13
      H = X11/Q12
      P12 =  A12*( X11 - H*Q12 ) - H*R12
      IF ( P13 < 0 ) THEN
         P13 = P13 + M1
      END IF
      IF ( P12 < 0 ) THEN
         P12 = P12 + M1
      END IF
      X10 = X11 
      X11 = X12
      X12 = P12 - P13
      IF ( X12 < 0 ) THEN
         X12 = X12 + M1
      END IF
!
!     Component 2
!
      H = X20/Q23
      P23 = -A23*( X20 - H*Q23 ) - H*R23
      H = X22/Q21
      P21 =  A21*( X22 - H*Q21 ) - H*R21
      IF ( P23 < 0 ) THEN
         P23 = P23 + M2
      END IF
      IF ( P21 < 0 ) THEN
         P21 = P21 + M2
      END IF
      X20 = X21 
      X21 = X22
      X22 = P21 - P23
      IF ( X22 < 0 ) THEN
         X22 = X22 + M2
      END IF
!
!     Combination
!
      Z = X12 - X22
      IF ( Z <= 0 ) THEN
         Z = Z + M1
      END IF
      UNI = Z*INVMP1
!
      END FUNCTION UNIFRM
!
      FUNCTION MVDNST( NU, X ) RESULT(DNSTY)
!
!     Student's T density function
!
      INTEGER,         INTENT(IN) :: NU
      REAL(KIND=STND), INTENT(IN) :: X
      REAL(KIND=STND)             :: DNSTY
!
      REAL(KIND=STND), PARAMETER :: PI = 3.141592653589793E0_STND
      REAL(KIND=STND), PARAMETER :: SQTWPI = 2.506628274631001E0_STND
      REAL(KIND=STND)            :: PROD
      INTEGER :: I
      DNSTY = 0
      IF ( NU < 1 ) THEN
         IF ( ABS(X) < 10 ) DNSTY = EXP( -X*X/2 )/SQTWPI
      ELSE
         PROD = 1/SQRT( REAL( NU, STND ) )
         DO I = NU - 2, 1, -2
            PROD = PROD*( I + 1 )/I
         END DO
         IF ( MODULO( NU, 2 ) == 0 ) THEN
            PROD = PROD/2
         ELSE
            PROD = PROD/PI
         END IF
         DNSTY = PROD/SQRT( 1 + X*X/NU )**( NU + 1 )
      END IF
!
      END FUNCTION MVDNST
!
      FUNCTION MVUVT( NU, A, B, INFIN ) RESULT(VALUE)
!
!     Univariate Distribution Function
!
!     Parameters
!
      INTEGER,         INTENT(IN)  :: NU, INFIN
      REAL(KIND=STND), INTENT(IN)  :: A, B
      REAL(KIND=STND)              :: LOWER, UPPER, VALUE
!
      LOWER = 0
      UPPER = 1
      IF ( INFIN >= 0 ) THEN
         IF ( INFIN /= 0 ) LOWER = MVSTDT( NU, A )
         IF ( INFIN /= 1 ) UPPER = MVSTDT( NU, B )
      END IF
      UPPER = MAX( UPPER, LOWER )
      VALUE = UPPER - LOWER
!
      END FUNCTION MVUVT
!
      FUNCTION MVSTDT( NU, T ) RESULT(DSTRB)
!
!     Student t Distribution Function
!
!                       T
!         MVSTDT = C   I  ( 1 + y*y/NU )**( -(NU+1)/2 ) dy
!                   NU -00
!
      INTEGER,         INTENT(IN) :: NU
      REAL(KIND=STND), INTENT(IN) :: T
      REAL(KIND=STND)             :: DSTRB
!
      REAL(KIND=STND), PARAMETER  :: PI = 3.141592653589793E0_STND 
      REAL(KIND=STND)             :: CSTHE, SNTHE, POLYN, TT, TS, RN
      INTEGER :: J
      IF ( NU < 1 ) THEN
         DSTRB = MVPHI( T )
      ELSE IF ( NU == 1 ) THEN
         DSTRB = ( 1 + 2*ATAN( T )/PI )/2
      ELSE IF ( NU == 2) THEN
         DSTRB = ( 1 + T/SQRT( 2 + T*T ))/2
      ELSE 
         TT = T*T
         CSTHE = NU/( NU + TT )
         POLYN = 1
         DO J = NU - 2, 2, -2
            POLYN = 1 + ( J - 1 )*CSTHE*POLYN/J
         END DO
         IF ( MODULO( NU, 2 ) == 1 ) THEN
            RN = NU
            TS = T/SQRT(RN)
            DSTRB = ( 1 + 2*( ATAN( TS ) + TS*CSTHE*POLYN )/PI )/2
         ELSE
            SNTHE = T/SQRT( NU + TT )
            DSTRB = ( 1 + SNTHE*POLYN )/2
         END IF
         IF ( DSTRB < 0 ) DSTRB = 0
      ENDIF
!
      END FUNCTION MVSTDT
!
      FUNCTION MVSTNV( N, Z ) RESULT(STINV)
!
!     Inverse Student t Distribution Function
!
!                    STINV
!           Z = C   I      (1 + y*y/N)**(-(N+1)/2) dy
!                N  -INF
!
!      Reference: G.W. Hill, Comm. ACM Algorithm 395
!                 Comm. ACM 13 (1970), pp. 619-620.
!
!      Enhancements for double precision and other modifications by
!                 Alan Genz, 1993-2001.
!
      INTEGER,         INTENT(IN) :: N
      REAL(KIND=STND), INTENT(IN) :: Z
      REAL(KIND=STND)             :: STINV
!
      INTEGER                     :: J
      REAL(KIND=STND)             :: P, A, B, C, D, X, XX, Y, CONST, STJAC
      REAL(KIND=STND),       SAVE :: NOLD = 0
      REAL(KIND=STND),  PARAMETER :: PI = 3.141592653589793E0_STND, TWO = 2 
!
      IF ( 0 < Z .AND. Z < 1 ) THEN
         IF ( N < 1 ) THEN
           STINV = MVPHNV( Z )
         ELSE IF ( N == 1 ) THEN
            STINV = TAN( PI*( 2*Z - 1 )/2 )
         ELSE IF ( N == 2) THEN
            STINV = ( 2*Z - 1 )/SQRT( 2*Z*( 1 - Z ) )
         ELSE 
            IF ( 2*Z >= 1 ) THEN 
               P = 2*( 1 - Z )
            ELSE
               P = 2*Z
            END IF
            A = 1/( N - 1/TWO )
            B = 48/( A*A )
            C = ( ( 20700*A/B - 98 )*A - 16 )*A + 96.36E0_STND
            D = ( ( 94.5E0_STND/( B + C ) - 3 )/B + 1 )*SQRT( A*PI/2 )*N
            X = D*P
            Y = X**( TWO/N )
            IF ( Y .GT. A + 0.05E0_STND ) THEN
               X = MVPHNV( P/2 )
               Y = X*X
               IF ( N < 5 ) C = C + 3*( N - 9/TWO )*( 10*X + 6 )/100
               C = ( ( (D*X - 100)*X/20 - 7 )*X - 2 )*X + B + C
               Y = ( ( ( ( (4*Y+63)*Y/10+36 )*Y+94.5E0_STND )/C-Y-3 )/B+1 )*X
               Y = A*Y*Y
               IF ( Y > 0.002E0_STND ) THEN
                  Y = EXP(Y) - 1
               ELSE
                  Y = Y*( 1 + Y/2 )
               ENDIF
            ELSE
               Y = ( ( 1/( ( (N+6)/(N*Y) - 0.089E0_STND*D - 0.822E0_STND )   &
                    *(3*N+6) ) + 0.5E0_STND/(N+4) )*Y - 1 )*(N+1)/(N+2) + 1/Y
            END IF
            STINV = SQRT(N*Y)
            IF ( 2*Z < 1 ) STINV = -STINV
            IF ( ABS( STINV ) > 0 ) THEN
!
!     Use one third order Schroeder correction to single precision result
!
               X = STINV
               D = Z - MVSTDT( N, X )
               IF ( N /= NOLD ) THEN
                  NOLD = N
                  IF ( MODULO( N, 2 ) == 0 ) THEN
                     CONST = SQRT(NOLD)*2
                  ELSE
                     CONST = SQRT(NOLD)*PI
                  END IF
                  DO J = N - 2, 1, -2
                     CONST = J*CONST/( J + 1 )
                  END DO
               END IF
               XX = 1 + X*X/N
               STJAC = CONST*XX**( ( N + 1 )/TWO ) 
               Y = D*STJAC
               STINV = X + Y*( 1 + Y*(N+1)/( 2*( X + N/X ) ) )
            END IF
         END IF
      ELSE
!
!     Use cutoff values for Z near 0 or 1.
!
         STINV = SQRT( N/( 2E-16_STND*SQRT( 2*PI*N ) )**( TWO/N ) )
         IF ( 2*Z < 1 ) STINV = -STINV
      END IF
!
      END FUNCTION MVSTNV
!
      FUNCTION MVSTDC( NU, L, U, INFN ) RESULT(DSTRBC)
!
!     Complementary Student t Distribution Function
!
!     if INFN, then MVSTDC computes
!         < 0        0
!          0     P( X > U )
!          1     P( X < L )
!          2     P( X > U ) + P( X < L )
!
      INTEGER,         INTENT(IN) :: NU, INFN
      REAL(KIND=STND), INTENT(IN) :: L, U
      REAL(KIND=STND)             :: DSTRBC
!
      IF ( INFN == 0 ) THEN
         DSTRBC = MVSTDT( NU, -U )
      ELSE IF ( INFN == 1 ) THEN
         DSTRBC = MVSTDT( NU,  L )
      ELSE IF ( INFN == 2 ) THEN
         DSTRBC = MVSTDT( NU,  L ) + MVSTDT( NU, -U )
      ELSE
         DSTRBC = 0
      END IF
!
      END FUNCTION MVSTDC
!
      FUNCTION MVPHI( Z ) RESULT(P)
!     
!     Normal distribution probabilities accurate to 1.e-15.
!     Z = no. of standard deviations from the mean.
!     
!     Based upon algorithm 5666 for the error function, from:
!     Hart, J.F. et al, 'Computer Approximations', Wiley 1968
!     
!     Programmer: Alan Miller
!     
!     Latest revision - 30 March 1986
!     
      REAL(KIND=STND), INTENT(IN) :: Z
      REAL(KIND=STND)             :: P
!
      REAL(KIND=STND) :: EXPNTL, ZABS
      REAL(KIND=STND), PARAMETER ::         &
           P0 =   220.2068679123761E0_STND, &
           P1 =   221.2135961699311E0_STND, &
           P2 =   112.0792914978709E0_STND, &
           P3 =   33.91286607838300E0_STND, &
           P4 =   6.373962203531650E0_STND, &
           P5 =  0.7003830644436881E0_STND, &
           P6 = 0.03526249659989109E0_STND
      REAL(KIND=STND), PARAMETER ::         &
           Q0 =   440.4137358247522E0_STND, &
           Q1 =   793.8265125199484E0_STND, &
           Q2 =   637.3336333788311E0_STND, &
           Q3 =   296.5642487796737E0_STND, &
           Q4 =   86.78073220294608E0_STND, &
           Q5 =   16.06417757920695E0_STND, &
           Q6 =   1.755667163182642E0_STND, &
           Q7 = 0.08838834764831844E0_STND
!
      REAL(KIND=STND), PARAMETER :: CUTOFF = 7.071067811865475_STND,         &
                                    ROOTPI = 2.506628274631000502415765_STND 
!     
      ZABS = ABS(Z)
!     
!     |Z| > 37
!     
      IF ( ZABS > 37 ) THEN
         P = 0
      ELSE
!     
!     |Z| <= 37
!     
         EXPNTL = EXP( -ZABS**2/2 )
!     
!     |Z| < CUTOFF = 10/SQRT(2)
!     
         IF ( ZABS < CUTOFF ) THEN
            P = EXPNTL*( (((((P6*ZABS + P5)*ZABS + P4)*ZABS + P3)*ZABS      & 
                + P2)*ZABS + P1)*ZABS + P0)/(((((((Q7*ZABS + Q6)*ZABS       & 
                + Q5)*ZABS + Q4)*ZABS + Q3)*ZABS + Q2)*ZABS + Q1)*ZABS + Q0 )
!     
!     |Z| >= CUTOFF.
!     
         ELSE
            P = EXPNTL/( ZABS + 1/( ZABS + 2/( ZABS + 3/( ZABS               &
                              + 4/( ZABS + 0.65E0_STND ) ) ) ) )/ROOTPI
         END IF
      END IF
      IF ( Z > 0 ) THEN
         P = 1 - P
      END IF
!
      END FUNCTION MVPHI
!
      FUNCTION MVPHNV( P ) RESULT(PHINV)
!
!	ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3
!
!	Produces the normal deviate Z corresponding to a given lower
!	tail area of P.
!
      REAL(KIND=STND), INTENT(IN) :: P
      REAL(KIND=STND)             :: PHINV
!
      REAL(KIND=STND) :: Q, R
      REAL(KIND=STND), PARAMETER :: SPLIT1 = 0.425E0_STND, SPLIT2 = 5,      &
                           CONST1 = 0.180625E0_STND, CONST2 = 1.6E0_STND 
!     
!     Coefficients for P close to 0.5
!     
      REAL(KIND=STND), PARAMETER ::            & 
           A0 = 3.3871328727963666080E00_STND, &
           A1 = 1.3314166789178437745E+2_STND, &
           A2 = 1.9715909503065514427E+3_STND, &
           A3 = 1.3731693765509461125E+4_STND, &
           A4 = 4.5921953931549871457E+4_STND, &
           A5 = 6.7265770927008700853E+4_STND, &
           A6 = 3.3430575583588128105E+4_STND, &
           A7 = 2.5090809287301226727E+3_STND, &
           B1 = 4.2313330701600911252E+1_STND, &
           B2 = 6.8718700749205790830E+2_STND, &
           B3 = 5.3941960214247511077E+3_STND, &
           B4 = 2.1213794301586595867E+4_STND, &
           B5 = 3.9307895800092710610E+4_STND, &
           B6 = 2.8729085735721942674E+4_STND, &
           B7 = 5.2264952788528545610E+3_STND  
!     
!     Coefficients for P not close to 0, 0.5 or 1.
!     
      REAL(KIND=STND), PARAMETER ::             & 
           C0 = 1.42343711074968357734E00_STND, &
           C1 = 4.63033784615654529590E00_STND, &
           C2 = 5.76949722146069140550E00_STND, &
           C3 = 3.64784832476320460504E00_STND, & 
           C4 = 1.27045825245236838258E00_STND, &
           C5 = 2.41780725177450611770E-1_STND, &
           C6 = 2.27238449892691845833E-2_STND, &
           C7 = 7.74545014278341407640E-4_STND, &
           D1 = 2.05319162663775882187E00_STND, &
           D2 = 1.67638483018380384940E00_STND, &
           D3 = 6.89767334985100004550E-1_STND, &
           D4 = 1.48103976427480074590E-1_STND, &
           D5 = 1.51986665636164571966E-2_STND, &
           D6 = 5.47593808499534494600E-4_STND, &
           D7 = 1.05075007164441684324E-9_STND  
!
!	Coefficients for P near 0 or 1.
!
      REAL(KIND=STND), PARAMETER ::             &
           E0 = 6.65790464350110377720E00_STND, &
           E1 = 5.46378491116411436990E00_STND, &
           E2 = 1.78482653991729133580E00_STND, &
           E3 = 2.96560571828504891230E-1_STND, &
           E4 = 2.65321895265761230930E-2_STND, &
           E5 = 1.24266094738807843860E-3_STND, &
           E6 = 2.71155556874348757815E-5_STND, &
           E7 = 2.01033439929228813265E-7_STND, &
           F1 = 5.99832206555887937690E-1_STND, &
           F2 = 1.36929880922735805310E-1_STND, &
           F3 = 1.48753612908506148525E-2_STND, &
           F4 = 7.86869131145613259100E-4_STND, &
           F5 = 1.84631831751005468180E-5_STND, &
           F6 = 1.42151175831644588870E-7_STND, &
           F7 = 2.04426310338993978564E-15_STND 
!     
      Q = ( 2*P - 1 )/2
      IF ( ABS(Q) <= SPLIT1 ) THEN
         R = CONST1 - Q*Q
         PHINV = Q*( ( ( ((((A7*R + A6)*R + A5)*R + A4)*R + A3) &
                       *R + A2 )*R + A1 )*R + A0 )              &
                 /( ( ( ((((B7*R + B6)*R + B5)*R + B4)*R + B3)  &
                       *R + B2 )*R + B1 )*R + 1 )
      ELSE
         R = MIN( P, 1 - P )
         IF ( R > 0 ) THEN
            R = SQRT( -LOG(R) )
            IF ( R <= SPLIT2 ) THEN
               R = R - CONST2
               PHINV = ( ( ( ((((C7*R + C6)*R + C5)*R + C4)*R + C3) &
                           *R + C2 )*R + C1 )*R + C0 )              &
                     /( ( ( ((((D7*R + D6)*R + D5)*R + D4)*R + D3)  &
                            *R + D2 )*R + D1 )*R + 1 )
            ELSE
               R = R - SPLIT2
               PHINV = ( ( ( ((((E7*R + E6)*R + E5)*R + E4)*R + E3) &
                           *R + E2 )*R + E1 )*R + E0 )              &
                     /( ( ( ((((F7*R + F6)*R + F5)*R + F4)*R + F3)  &
                           *R + F2 )*R + F1 )*R + 1 ) 
            END IF
         ELSE
            PHINV = 9
         END IF
         IF ( Q < 0 ) THEN
            PHINV = - PHINV
         END IF
      END IF
!
      END FUNCTION MVPHNV
!
      FUNCTION MVCHNV( N, P ) RESULT(R)    
!
!                  MVCHNV
!     P =  1 - K  I     exp(-t*t/2) t**(N-1) dt, for N >= 1.
!               N  0
!
      INTEGER,         INTENT(IN) :: N      
      REAL(KIND=STND), INTENT(IN) :: P
      REAL(KIND=STND)             :: R
!
      INTEGER                     :: I
      REAL(KIND=STND)             :: RI
      REAL(KIND=STND), PARAMETER  :: LRP = -.22579135264472743235E0_STND
!                                    LRP =   LOG( SQRT( 2/PI ) )
      REAL(KIND=STND),       SAVE :: LKN
      INTEGER,               SAVE :: NO = 0
      IF ( N < 1 ) THEN
         R = 1
      ELSE IF ( N == 1 ) THEN
         R = -MVPHNV( P/2 )
      ELSE IF ( P < 1 ) THEN
         IF ( N == 2 ) THEN
            R = SQRT( -2*LOG(P) )
         ELSE
            IF ( N /= NO ) THEN
               NO = N
               LKN = 0
               DO I = N - 2, 2, -2
                  RI = I
                  LKN = LKN - LOG( RI )
               END DO
               IF ( MODULO( N, 2 ) == 1 ) LKN = LKN + LRP
            END IF
            IF ( N >= -5*LOG( 1 - P )/4 ) THEN
               R = 2E0_STND/( 9*N )
               R = N*( -MVPHNV(P)*SQRT(R) + 1 - R )**3
               IF ( R > 2*N+6 ) THEN
                  R = 2*( LKN - LOG(P) ) + ( N - 2 )*LOG(R)
               END IF
            ELSE
               R = EXP( 2*( LOG( ( 1 - P )*N ) - LKN )/N )
            END IF
            R = SQRT(R)
            RI = R
            R = MVCHNC( LKN, N, P, R )
            IF ( ABS( R - RI ) > 2E-1_STND ) THEN
               R = RI
            ELSE 
               DO I = 1, 5
                  IF ( ABS( R - RI ) < 1E-6_STND ) EXIT
                  RI = R
                  R = MVCHNC( LKN, N, P, R )
               END DO
            END IF
         END IF
      ELSE
         R = 0
      END IF
!
      END FUNCTION MVCHNV
!
      FUNCTION MVCHNC( LKN, N, P, R ) RESULT(CHNC)
!
!     Third order Schroeder correction to R for MVCHNV
!
      INTEGER,         INTENT(IN) :: N      
      REAL(KIND=STND), INTENT(IN) :: P, LKN, R
      REAL(KIND=STND)             :: CHNC
!
      REAL(KIND=STND)            :: DF
      DF = P - MVCHI( N, R )
      DF = MAX( P - 1, MIN( DF, P ) )
      DF =  DF/EXP( LKN + ( N - 1 )*LOG(R) - R*R/2 )
      CHNC = R - DF*( 1 - DF*( R - ( N - 1 )/R )/2 )
!
      END FUNCTION MVCHNC
!
      FUNCTION MVCHI( N, R ) RESULT(CHI)
!
!                     R
!      CHI =  1 - K  I     exp(-t*t/2) t**(N-1) dt, for N >= 1.
!                  N  0
!
      INTEGER,         INTENT(IN) :: N      
      REAL(KIND=STND), INTENT(IN) :: R
      REAL(KIND=STND)             :: CHI
!
      INTEGER                     :: I
      REAL(KIND=STND)             :: RN, AL, RI, RR, AI, BI, CI, DI, DL 
      REAL(KIND=STND), PARAMETER  :: RP =  0.79788456080286535588E0_STND
!                                    RP =  SQRT( 2/PI )
      REAL(KIND=STND), PARAMETER  :: LRP = 0.12078223763524522234E0_STND
!                                    LRP =  LOG( 2/SQRT(PI) )
      REAL(KIND=STND), PARAMETER  :: ONE = 1, EPS = EPSILON(ONE)
      REAL(KIND=STND),       SAVE :: LGN
      INTEGER,               SAVE :: NO = 0
      RR = R*R
      IF ( N < 2 ) THEN
         CHI = 2*MVPHI(-R)
      ELSE IF ( N < 100 ) THEN
!
!     Use standard Chi series
!
         RN = 1
         DO I = N - 2, 2, -2
            RN = 1 + RR*RN/I
         END DO
         IF ( MODULO( N, 2 ) == 0 ) THEN
            CHI = EXP( LOG(      RN ) - RR/2 )
         ELSE
            CHI = EXP( LOG( R*RP*RN ) - RR/2 ) + 2*MVPHI(-R)
         ENDIF
      ELSE
         IF ( N /= NO ) THEN
            NO = N
            LGN = 0
            DO I = N - 2, 2, -2
               RI = I
               LGN = LGN + LOG( 2/RI )
            END DO
            IF ( MODULO( N, 2 ) == 1 ) LGN = LGN + LRP
         END IF
         RR = RR/2
         AL = N
         AL = AL/2
         IF ( RR < AL + 1 ) THEN 
!
!           Use Incomplete Gamma series
!
            DL = EXP( -RR + AL*LOG(RR) + LGN )/AL
            CHI = DL
            DO I = 1, 1000
               DL = DL*RR/( AL + I ) 
               CHI = CHI + DL
               IF ( ABS( DL*RR/( AL + I + 1 - RR ) ) < EPS ) EXIT
            END DO
            CHI = 1 - CHI
         ELSE
!
!           Use Incomplete Gamma continued fraction
!
            BI = RR + 1 - AL
            CI = 1/EPS
            DI = BI
            CHI = EXP( -RR + AL*LOG(RR) + LGN )/BI 
            DO I = 1, 200
               AI = I*( AL - I )
               BI = BI + 2
               CI = BI + AI/CI
               IF ( CI == 0 ) CI = EPS 
               DI = BI + AI/DI
               IF ( DI == 0 ) DI = EPS 
               DL = CI/DI
               CHI = CHI*DL
               IF ( ABS( DL - 1 ) < EPS ) EXIT
            END DO
         END IF
      END IF
!
      END FUNCTION MVCHI
!
      FUNCTION MVBVT( NU, LOWER, UPPER, INFIN, CORREL ) RESULT(BVNBVT)
!
!     A function for computing bivariate normal and t probabilities.
!
!  Parameters
!
!     NU     INTEGER degrees of freedom parameter; NU < 1 gives normal case.
!     LOWER  REAL, array of lower integration limits.
!     UPPER  REAL, array of upper integration limits.
!     INFIN  INTEGER, array of integration limits flags:
!            if INFIN(I) = 0, Ith limits are (-infinity, UPPER(I)];
!            if INFIN(I) = 1, Ith limits are [LOWER(I), infinity);
!            if INFIN(I) = 2, Ith limits are [LOWER(I), UPPER(I)].
!     CORREL REAL, correlation coefficient.
!
      INTEGER,                       INTENT(IN) :: NU
      REAL(KIND=STND), DIMENSION(:), INTENT(IN) :: LOWER, UPPER
      INTEGER,         DIMENSION(:), INTENT(IN) :: INFIN
      REAL(KIND=STND),               INTENT(IN) :: CORREL
      REAL(KIND=STND)                           :: BVNBVT
!
      IF ( INFIN(1) < 0 ) THEN
         BVNBVT = MVUVT( NU, LOWER(2), UPPER(2), INFIN(2) )
      ELSE IF ( INFIN(2) < 0 ) THEN
         BVNBVT = MVUVT( NU, LOWER(1), UPPER(1), INFIN(1) )
      ELSE IF ( INFIN(1) == 2  .AND. INFIN(2) == 2 ) THEN
         BVNBVT =  MVBVTL ( NU, UPPER(1), UPPER(2), CORREL )             &
              - MVBVTL ( NU, UPPER(1), LOWER(2), CORREL )                &
              - MVBVTL ( NU, LOWER(1), UPPER(2), CORREL )                &
              + MVBVTL ( NU, LOWER(1), LOWER(2), CORREL )
      ELSE IF ( INFIN(1) == 2  .AND. INFIN(2) == 1 ) THEN
         BVNBVT =  MVBVTL ( NU, -LOWER(1), -LOWER(2), CORREL )           &
              - MVBVTL ( NU, -UPPER(1), -LOWER(2), CORREL )
      ELSE IF ( INFIN(1) == 1  .AND. INFIN(2) == 2 ) THEN
         BVNBVT =  MVBVTL ( NU, -LOWER(1), -LOWER(2), CORREL )           &
              - MVBVTL ( NU, -LOWER(1), -UPPER(2), CORREL )
      ELSE IF ( INFIN(1) == 2  .AND. INFIN(2) == 0 ) THEN
         BVNBVT =  MVBVTL ( NU, UPPER(1), UPPER(2), CORREL )             &
              - MVBVTL ( NU, LOWER(1), UPPER(2), CORREL )
      ELSE IF ( INFIN(1) == 0  .AND. INFIN(2) == 2 ) THEN
         BVNBVT =  MVBVTL ( NU, UPPER(1), UPPER(2), CORREL )             &
              - MVBVTL ( NU, UPPER(1), LOWER(2), CORREL )
      ELSE IF ( INFIN(1) == 1  .AND. INFIN(2) == 0 ) THEN
         BVNBVT =  MVBVTL ( NU, -LOWER(1), UPPER(2), -CORREL )
      ELSE IF ( INFIN(1) == 0  .AND. INFIN(2) == 1 ) THEN
         BVNBVT =  MVBVTL ( NU, UPPER(1), -LOWER(2), -CORREL )
      ELSE IF ( INFIN(1) == 1  .AND. INFIN(2) == 1 ) THEN
         BVNBVT =  MVBVTL ( NU, -LOWER(1), -LOWER(2), CORREL )
      ELSE IF ( INFIN(1) == 0  .AND. INFIN(2) == 0 ) THEN
         BVNBVT =  MVBVTL ( NU, UPPER(1), UPPER(2), CORREL )
      END IF
!
      END FUNCTION MVBVT
!
      FUNCTION MVBVTL( NU, DH, DK, R ) RESULT(BVT)
!
!     a function for computing bivariate t probabilities.
!
!       Alan Genz
!       Department of Mathematics
!       Washington State University
!       Pullman, Wa 99164-3113
!       Email : alangenz@wsu.edu
!
!    this function is based on the method described by 
!        Dunnett, C.W. and M. Sobel, (1954),
!        A bivariate generalization of Student's t-distribution
!        with tables for certain special cases,
!        Biometrika 41, pp. 153-169.
!
! mvbvtl - calculate the probability that x < dh and y < dk. 
!
! parameters
!
!   nu number of degrees of freedom
!   dh 1st lower integration limit
!   dk 2nd lower integration limit
!   r   correlation coefficient
!
      REAL(KIND=STND), INTENT(IN) :: DH, DK, R
      REAL(KIND=STND)             :: BVT
      INTEGER, INTENT(IN)         :: NU      
!
      INTEGER         :: J, HS, KS
      REAL(KIND=STND) :: ORS, HRK, KRH, SQNU, GMPH, GMPK 
      REAL(KIND=STND) :: XNKH, XNHK, QHRK, HKN, HPK, HKRN
      REAL(KIND=STND) :: BTNCKH, BTNCHK, BTPDKH, BTPDHK
      REAL(KIND=STND), PARAMETER :: PI = 3.14159265358979323844E0_STND,      &
                                    TPI = 2*PI, ONE = 1 
      IF ( NU < 1 ) THEN
         BVT = MVBVU( -DH, -DK, R )
      ELSE
         SQNU = NU
         SQNU = SQRT( SQNU )
         ORS = 1 - R*R  
         HRK = DH - R*DK  
         KRH = DK - R*DH  
         IF ( ABS(HRK) + ORS > 0 ) THEN
            XNHK = HRK**2/( HRK**2 + ORS*( NU + DK**2 ) ) 
            XNKH = KRH**2/( KRH**2 + ORS*( NU + DH**2 ) ) 
         ELSE
            XNHK = 0
            XNKH = 0  
         END IF
         HS = SIGN( ONE, DH - R*DK )  
         KS = SIGN( ONE, DK - R*DH ) 
         IF ( MODULO( NU, 2 ) == 0 ) THEN
            BVT = ATAN2( SQRT(ORS), -R )/TPI 
            GMPH = DH/SQRT( 16*( NU + DH**2 ) )  
            GMPK = DK/SQRT( 16*( NU + DK**2 ) )  
            BTNCKH = 2*ATAN2( SQRT( XNKH ), SQRT( 1 - XNKH ) )/PI  
            BTPDKH = 2*SQRT( XNKH*( 1 - XNKH ) )/PI 
            BTNCHK = 2*ATAN2( SQRT( XNHK ), SQRT( 1 - XNHK ) )/PI  
            BTPDHK = 2*SQRT( XNHK*( 1 - XNHK ) )/PI 
            DO J = 1, NU/2
               BVT = BVT + GMPH*( 1 + KS*BTNCKH ) 
               BVT = BVT + GMPK*( 1 + HS*BTNCHK ) 
               BTNCKH = BTNCKH + BTPDKH  
               BTPDKH = 2*J*BTPDKH*( 1 - XNKH )/( 2*J + 1 )  
               BTNCHK = BTNCHK + BTPDHK  
               BTPDHK = 2*J*BTPDHK*( 1 - XNHK )/( 2*J + 1 )  
               GMPH = GMPH*( 2*J - 1 )/( 2*J*( 1 + DH**2/NU ) ) 
               GMPK = GMPK*( 2*J - 1 )/( 2*J*( 1 + DK**2/NU ) ) 
            END DO
         ELSE
            QHRK = SQRT( DH**2 + DK**2 - 2*R*DH*DK + NU*ORS )  
            HKRN = DH*DK + R*NU  
            HKN = DH*DK - NU  
            HPK = DH + DK 
            BVT = ATAN2(-SQNU*(HKN*QHRK+HPK*HKRN),HKN*HKRN-NU*HPK*QHRK)/TPI  
            IF ( BVT < -1E-15_STND ) BVT = BVT + 1
            GMPH = DH/( TPI*SQNU*( 1 + DH**2/NU ) )  
            GMPK = DK/( TPI*SQNU*( 1 + DK**2/NU ) )  
            BTNCKH = SQRT( XNKH )  
            BTPDKH = BTNCKH 
            BTNCHK = SQRT( XNHK )  
            BTPDHK = BTNCHK  
            DO J = 1, ( NU - 1 )/2
               BVT = BVT + GMPH*( 1 + KS*BTNCKH ) 
               BVT = BVT + GMPK*( 1 + HS*BTNCHK ) 
               BTPDKH = ( 2*J - 1 )*BTPDKH*( 1 - XNKH )/( 2*J )  
               BTNCKH = BTNCKH + BTPDKH  
               BTPDHK = ( 2*J - 1 )*BTPDHK*( 1 - XNHK )/( 2*J )  
               BTNCHK = BTNCHK + BTPDHK  
               GMPH = 2*J*GMPH/( ( 2*J + 1 )*( 1 + DH**2/NU ) ) 
               GMPK = 2*J*GMPK/( ( 2*J + 1 )*( 1 + DK**2/NU ) ) 
            END DO
         END IF
      END IF
!
      END FUNCTION MVBVTL
!
      FUNCTION MVBVU( SH, SK, R ) RESULT(BVN)
!
!     A function for computing bivariate normal probabilities.
!
!       Yihong Ge
!     and
!       Alan Genz
!       Department of Mathematics
!       Washington State University
!       Pullman, WA 99164-3113
!       Email : alangenz@wsu.edu
!
!    MVBVU - calculate the probability that X is larger than SH and Y is
!          larger than SK.
!
! Parameters
!
!   SH  REAL, integration limit
!   SK  REAL, integration limit
!   R   REAL, correlation coefficient
!   LG  INTEGER, number of Gauss Rule Points and Weights
!
      REAL(KIND=STND), INTENT(IN) :: SH, SK, R
      REAL(KIND=STND)             :: BVN
!
      INTEGER :: I, LG, NG
      REAL(KIND=STND), PARAMETER :: ZERO = 0,                              & 
                                    TWOPI = 6.2831853071795864769252_STND
      REAL(KIND=STND) :: AS, A, B, C, D, RS, XS, SN, ASR, H, K, BS, HS, HK
      REAL(KIND=STND), DIMENSION(10) :: X, W
!     Gauss Legendre Points and Weights, N =  6
      REAL(KIND=STND), DIMENSION(3), PARAMETER :: W1 = (/                  &
        0.1713244923791705E+00_STND,                                       &
        0.3607615730481384E+00_STND,  0.4679139345726904E+00_STND /)
      REAL(KIND=STND), DIMENSION(3), PARAMETER :: X1 = (/                  &
       -0.9324695142031522E+00_STND,                                       &
       -0.6612093864662647E+00_STND, -0.2386191860831970E+00_STND /)
!     Gauss Legendre Points and Weights, N = 12
      REAL(KIND=STND), DIMENSION(6), PARAMETER :: W2 = (/                  &   
        0.4717533638651177E-01_STND,  0.1069393259953183E+00_STND,         &   
        0.1600783285433464E+00_STND,  0.2031674267230659E+00_STND,         &   
        0.2334925365383547E+00_STND,  0.2491470458134029E+00_STND /)
      REAL(KIND=STND), DIMENSION(6), PARAMETER :: X2 = (/                  &   
       -0.9815606342467191E+00_STND, -0.9041172563704750E+00_STND,         &   
       -0.7699026741943050E+00_STND, -0.5873179542866171E+00_STND,         &   
       -0.3678314989981802E+00_STND, -0.1252334085114692E+00_STND /)
!     Gauss Legendre Points and Weights, N = 20
      REAL(KIND=STND), DIMENSION(10), PARAMETER :: W3 = (/                 &   
        0.1761400713915212E-01_STND,  0.4060142980038694E-01_STND,         &   
        0.6267204833410906E-01_STND,  0.8327674157670475E-01_STND,         &   
        0.1019301198172404E+00_STND,  0.1181945319615184E+00_STND,         &   
        0.1316886384491766E+00_STND,  0.1420961093183821E+00_STND,         &   
        0.1491729864726037E+00_STND,  0.1527533871307259E+00_STND /)
      REAL(KIND=STND), DIMENSION(10), PARAMETER :: X3 = (/                 &   
       -0.9931285991850949E+00_STND, -0.9639719272779138E+00_STND,         &   
       -0.9122344282513259E+00_STND, -0.8391169718222188E+00_STND,         &   
       -0.7463319064601508E+00_STND, -0.6360536807265150E+00_STND,         &   
       -0.5108670019508271E+00_STND, -0.3737060887154196E+00_STND,         &   
       -0.2277858511416451E+00_STND, -0.7652652113349733E-01_STND /)
!
      IF ( ABS(R) < 3E-1_STND ) THEN
         NG = 1
         LG = 3
         X(1:LG) = X1
         W(1:LG) = W1
      ELSE IF ( ABS(R) < 75E-2_STND ) THEN
         NG = 2
         LG = 6
         X(1:LG) = X2
         W(1:LG) = W2
      ELSE 
         NG = 3
         LG = 10
         X(1:LG) = X3
         W(1:LG) = W3
      ENDIF
      H = SH
      K = SK 
      HK = H*K
      BVN = 0
      IF ( ABS(R) < 925E-3_STND ) THEN
         HS = ( H*H + K*K )/2
         ASR = ASIN(R)
         DO I = 1, LG
            SN = SIN( ASR*( 1 + X(I) )/2 )
            BVN = BVN + W(I)*EXP( ( SN*HK - HS )/( 1 - SN*SN ) )
            SN = SIN( ASR*( 1 - X(I) )/2 )
            BVN = BVN + W(I)*EXP( ( SN*HK - HS )/( 1 - SN*SN ) )
         END DO
         BVN = BVN*ASR/(2*TWOPI) + MVPHI(-H)*MVPHI(-K) 
      ELSE
         IF ( R < 0 ) THEN
            K = -K
            HK = -HK
         ENDIF
         IF ( ABS(R) < 1 ) THEN
            AS = ( 1 - R )*( 1 + R )
            A = SQRT(AS)
            BS = ( H - K )**2
            C = ( 4 - HK )/8 
            D = ( 12 - HK )/16
            BVN = A*EXP( -(BS/AS + HK)/2 )                               &
                   *( 1 - C*(BS - AS)*(1 - D*BS/5)/3 + C*D*AS*AS/5 )
            IF ( HK > -160 ) THEN
               B = SQRT(BS)
               BVN = BVN - EXP(-HK/2)*SQRT(TWOPI)*MVPHI(-B/A)*B          &
                                 *( 1 - C*BS*( 1 - D*BS/5 )/3 ) 
            ENDIF
            A = A/2
            DO I = 1, LG
               XS = ( A*( X(I) + 1 ) )**2
               RS = SQRT( 1 - XS )
               BVN = BVN + A*W(I)*                                       &
                    ( EXP( -BS/(2*XS) - HK/(1+RS) )/RS                   &
                    - EXP( -(BS/XS+HK)/2 )*( 1 + C*XS*( 1 + D*XS ) ) )
               XS = AS*( -X(I) + 1 )**2/4
               RS = SQRT( 1 - XS )
               BVN = BVN + A*W(I)*EXP( -(BS/XS + HK)/2 )                 &
                          *( EXP( -HK*(1-RS)/(2*(1+RS)) )/RS             &
                             - ( 1 + C*XS*( 1 + D*XS ) ) )
            END DO
            BVN = -BVN/TWOPI
         ENDIF
         IF ( R > 0 ) THEN
            BVN =  BVN + MVPHI( -MAX( H, K ) )
         ELSE
            BVN = -BVN + MAX( ZERO, MVPHI(-H) - MVPHI(-K) )     
         END IF
      END IF
!
      END FUNCTION MVBVU
!
      FUNCTION MVBVTC( NU, L, U, INFIN, RHO ) RESULT(BVTC)
!
!     A function for computing complementary bivariate normal and t 
!       probabilities.
!
!  Parameters
!
!     NU     INTEGER degrees of freedom parameter.
!     L      REAL, array of lower integration limits.
!     U      REAL, array of upper integration limits.
!     INFIN  INTEGER, array of integration limits flags:
!            if INFIN(1) INFIN(2),        then MVBVTC computes
!                 0         0              P( X>U(1), Y>U(2) ),
!                 1         0              P( X<L(1), Y>U(2) ),
!                 0         1              P( X>U(1), Y<L(2) ),
!                 1         1              P( X<L(1), Y<L(2) ),
!                 2         0      P( X>U(1), Y>U(2) ) + P( X<L(1), Y>U(2) ),
!                 2         1      P( X>U(1), Y<L(2) ) + P( X<L(1), Y<L(2) ),
!                 0         2      P( X>U(1), Y>U(2) ) + P( X>U(1), Y<L(2) ),
!                 1         2      P( X<L(1), Y>U(2) ) + P( X<L(1), Y<L(2) ),
!                 2         2      P( X>U(1), Y<L(2) ) + P( X<L(1), Y<L(2) ),
!                               +  P( X>U(1), Y>U(2) ) + P( X<L(1), Y>U(2) ).
!
!     RHO    REAL, correlation coefficient.
!
      INTEGER,                       INTENT(IN) :: NU
      REAL(KIND=STND), DIMENSION(:), INTENT(IN) :: L, U
      INTEGER,         DIMENSION(:), INTENT(IN) :: INFIN
      REAL(KIND=STND),               INTENT(IN) :: RHO
      REAL(KIND=STND)                           :: BVTC
!
!   Local Variables
!
      REAL(KIND=STND), DIMENSION(2) :: LW, UP
      REAL(KIND=STND)               :: B
      INTEGER,         DIMENSION(2) :: INF
      INTEGER                       :: I
!
      IF ( INFIN(1) < 0 .AND. INFIN(2) < 0 ) THEN
         BVTC = 0
      ELSE IF ( INFIN(1) < 0 ) THEN
         BVTC = MVSTDC( NU, L(2), U(2), INFIN(2) )
      ELSE IF ( INFIN(2) < 0 ) THEN
         BVTC = MVSTDC( NU, L(1), U(1), INFIN(1) )
      ELSE
         DO I = 1, 2
            IF ( MODULO( INFIN(I), 2 ) == 0 ) THEN
               INF(I) = 1
               LW(I) = U(I) 
            ELSE
               INF(I) = 0
               UP(I) = L(I) 
            END IF
         END DO
         B = MVBVT( NU, LW, UP, INF, RHO )
         DO I = 1, 2
            IF ( INFIN(I) == 2 ) THEN
               INF(I) = 0
               UP(I) = L(I) 
               B = B + MVBVT( NU, LW, UP, INF, RHO )
            END IF
         END DO
         IF ( INFIN(1) == 2 .AND. INFIN(2) == 2 ) THEN
            INF(1) = 1
            LW(1) = U(1) 
            B = B + MVBVT( NU, LW, UP, INF, RHO )
         END IF
         BVTC = B
      END IF
!
      END FUNCTION MVBVTC
!
      FUNCTION MVTVT( NU, LOWER, UPPER, INFIN, CR ) RESULT(TVNT)
!
!     A function for computing bivariate normal and t probabilities.
!
!  Parameters
!
!     NU     INTEGER degrees of freedom parameter; NU < 1 gives normal case.
!     LOWER  REAL, array of lower integration limits.
!     UPPER  REAL, array of upper integration limits.
!     INFIN  INTEGER, array of integration limits flags:
!            if INFIN(I) = 0, Ith limits are (-infinity, UPPER(I)];
!            if INFIN(I) = 1, Ith limits are [LOWER(I), infinity);
!            if INFIN(I) = 2, Ith limits are [LOWER(I), UPPER(I)].
!     CR REAL, correlation coefficients.
!
      INTEGER,                       INTENT(IN) :: NU
      REAL(KIND=STND), DIMENSION(:), INTENT(IN) :: LOWER, UPPER
      INTEGER,         DIMENSION(:), INTENT(IN) :: INFIN
      REAL(KIND=STND), DIMENSION(:), INTENT(IN) :: CR
      REAL(KIND=STND)                           :: TVNT
!
      INTEGER,         DIMENSION(2) :: INFT
      REAL(KIND=STND)               :: L1, L2, U1, U2
      REAL(KIND=STND), DIMENSION(2) :: LT, UT
      REAL(KIND=STND), DIMENSION(3) :: LW, UP
!
      IF ( INFIN(1) < 0 ) THEN
! -1 0 0
         LT = LOWER(2:3)
         UT = UPPER(2:3)
         INFT = INFIN(2:3)
         TVNT = MVBVT( NU, LT, UT, INFT, CR(3) ) 
      ELSE IF ( INFIN(2) < 0 ) THEN
!  *-1 *
         LT = (/ LOWER(1), LOWER(3) /) 
         UT = (/ UPPER(1), UPPER(3) /) 
         INFT = (/ INFIN(1), INFIN(3) /)
         TVNT = MVBVT( NU, LT, UT, INFT, CR(2) ) 
      ELSE IF ( INFIN(3) < 0 ) THEN
!  * *-1
         LT = LOWER(1:2)
         UT = UPPER(1:2)
         INFT = INFIN(1:2)
         TVNT = MVBVT( NU, LT, UT, INFT, CR(1) ) 
      ELSE IF ( INFIN(1) == 0 .AND. INFIN(2) == 0 .AND. INFIN(3) == 0 ) THEN
!  0 0 0
         TVNT = MVTVTL( NU, UPPER, CR )  
      ELSE IF ( INFIN(1) == 1 .AND. INFIN(2) == 0 .AND. INFIN(3) == 0 ) THEN
!  1 0 0
         U1 = UPPER(2)
         U2 = UPPER(3)
         UP = UPPER
         UP(1) = LOWER(1)
         TVNT = MVBVTL( NU, U1, U2, CR(3) ) - MVTVTL( NU, UP, CR )  
      ELSE IF ( INFIN(1) == 0 .AND. INFIN(2) == 1 .AND. INFIN(3) == 0 ) THEN
!  0 1 0
         U1 = UPPER(1)
         U2 = UPPER(3)
         UP = UPPER
         UP(2) = LOWER(2)
         TVNT = MVBVTL( NU, U1, U2, CR(2) ) - MVTVTL( NU, UP, CR )  
      ELSE IF ( INFIN(1) == 0 .AND. INFIN(2) == 0 .AND. INFIN(3) == 1 ) THEN
!  0 0 1
         U1 = UPPER(1)
         U2 = UPPER(2)
         UP = UPPER
         UP(3) = LOWER(3)
         TVNT = MVBVTL( NU, U1, U2, CR(1) ) - MVTVTL( NU, UP, CR )  
      ELSE IF ( INFIN(1) == 1 .AND. INFIN(2) == 1 .AND. INFIN(3) == 0 ) THEN
!  1 1 0
         U1 = -LOWER(1)
         U2 = -LOWER(2)
         UP = -LOWER
         UP(3) = -UPPER(3)
         TVNT = MVBVTL( NU, U1, U2, CR(1) ) - MVTVTL( NU, UP, CR )  
      ELSE IF ( INFIN(1) == 1 .AND. INFIN(2) == 0 .AND. INFIN(3) == 1 ) THEN
!  1 0 1
         U1 = -LOWER(1)
         U2 = -LOWER(3)
         UP = -LOWER
         UP(2) = -UPPER(2)
         TVNT = MVBVTL( NU, U1, U2, CR(2) ) - MVTVTL( NU, UP, CR )  
      ELSE IF ( INFIN(1) == 0 .AND. INFIN(2) == 1 .AND. INFIN(3) == 1 ) THEN
!  0 1 1
         U1 = -LOWER(2)
         U2 = -LOWER(3)
         UP = -LOWER
         UP(1) = -UPPER(1)
         TVNT = MVBVTL( NU, U1, U2, CR(3) ) - MVTVTL( NU, UP, CR )  
      ELSE IF ( INFIN(1) == 1 .AND. INFIN(2) == 1 .AND. INFIN(3) == 1 ) THEN
!  1 1 1
         UP = -LOWER
         TVNT = MVTVTL( NU, UP, CR )  
      ELSE IF ( INFIN(1) == 2 .AND. INFIN(2) == 0 .AND. INFIN(3) == 0 ) THEN
!  2 0 0
         UP = UPPER
         UP(1) = LOWER(1)
         TVNT = MVTVTL( NU, UPPER, CR ) - MVTVTL( NU, UP, CR )  
      ELSE IF ( INFIN(1) == 0 .AND. INFIN(2) == 2 .AND. INFIN(3) == 0 ) THEN
!  0 2 0
         UP = UPPER
         UP(2) = LOWER(2)
         TVNT = MVTVTL( NU, UPPER, CR ) - MVTVTL( NU, UP, CR )  
      ELSE IF ( INFIN(1) == 0 .AND. INFIN(2) == 0 .AND. INFIN(3) == 2 ) THEN
!  0 0 2
         UP = UPPER
         UP(3) = LOWER(3)
         TVNT = MVTVTL( NU, UPPER, CR ) - MVTVTL( NU, UP, CR )  
      ELSE IF ( INFIN(1) == 2 .AND. INFIN(2) == 2 .AND. INFIN(3) == 0 ) THEN
!  2 2 0
         UP = UPPER
         UP(1) = LOWER(1)
         TVNT = MVTVTL( NU, UPPER, CR ) - MVTVTL( NU, UP, CR )  
         UP = UPPER
         UP(2) = LOWER(2)
         LW = LOWER
         LW(3) = UPPER(3)
         TVNT = TVNT - MVTVTL( NU, UP, CR ) + MVTVTL( NU, LW, CR )
      ELSE IF ( INFIN(1) == 2 .AND. INFIN(2) == 0 .AND. INFIN(3) == 2 ) THEN
!  2 0 2
         UP = UPPER
         UP(1) = LOWER(1)
         TVNT = MVTVTL( NU, UPPER, CR ) - MVTVTL( NU, UP, CR )  
         UP = UPPER
         UP(3) = LOWER(3)
         LW = LOWER
         LW(2) = UPPER(2)
         TVNT = TVNT - MVTVTL( NU, UP, CR ) + MVTVTL( NU, LW, CR )
      ELSE IF ( INFIN(1) == 0 .AND. INFIN(2) == 2 .AND. INFIN(3) == 2 ) THEN
!  0 2 2
         UP = UPPER
         UP(2) = LOWER(2)
         TVNT = MVTVTL( NU, UPPER, CR ) - MVTVTL( NU, UP, CR )  
         UP = UPPER
         UP(3) = LOWER(3)
         LW = LOWER
         LW(1) = UPPER(1)
         TVNT = TVNT - MVTVTL( NU, UP, CR ) + MVTVTL( NU, LW, CR )
      ELSE IF ( INFIN(1) == 0 .AND. INFIN(2) == 1 .AND. INFIN(3) == 2 ) THEN
!  0 1 2
         U1 = UPPER(1)
         U2 = UPPER(3)
         L1 = UPPER(1)
         L2 = LOWER(3)
         TVNT = MVBVTL( NU, U1, U2, CR(2) ) - MVBVTL( NU, L1, L2, CR(2) ) 
         UP = UPPER
         UP(2) = LOWER(2)
         LW = LOWER
         LW(1) = UPPER(1)
         TVNT = TVNT - MVTVTL( NU, UP, CR ) + MVTVTL( NU, LW, CR ) 
      ELSE IF ( INFIN(1) == 0 .AND. INFIN(2) == 2 .AND. INFIN(3) == 1 ) THEN
!  0 2 1
         U1 = UPPER(1)
         U2 = UPPER(2)
         L1 = UPPER(1)
         L2 = LOWER(2)
         TVNT = MVBVTL( NU, U1, U2, CR(1) ) - MVBVTL( NU, L1, L2, CR(1) ) 
         UP = UPPER
         UP(3) = LOWER(3)
         LW = LOWER
         LW(1) = UPPER(1)
         TVNT = TVNT - MVTVTL( NU, UP, CR ) + MVTVTL( NU, LW, CR ) 
      ELSE IF ( INFIN(1) == 1 .AND. INFIN(2) == 0 .AND. INFIN(3) == 2 ) THEN
!  1 0 2
         U1 = UPPER(2)
         U2 = UPPER(3)
         L1 = UPPER(2)
         L2 = LOWER(3)
         TVNT = MVBVTL( NU, U1, U2, CR(3) ) - MVBVTL( NU, L1, L2, CR(3) ) 
         UP = UPPER
         UP(1) = LOWER(1)
         LW = LOWER
         LW(2) = UPPER(2)
         TVNT = TVNT - MVTVTL( NU, UP, CR ) + MVTVTL( NU, LW, CR ) 
      ELSE IF ( INFIN(1) == 2 .AND. INFIN(2) == 0 .AND. INFIN(3) == 1 ) THEN
!  2 0 1
         U1 = UPPER(1)
         U2 = UPPER(2)
         L1 = LOWER(1)
         L2 = UPPER(2)
         TVNT = MVBVTL( NU, U1, U2, CR(1) ) - MVBVTL( NU, L1, L2, CR(1) ) 
         UP = UPPER
         UP(3) = LOWER(3)
         LW = LOWER
         LW(2) = UPPER(2)
         TVNT = TVNT - MVTVTL( NU, UP, CR ) + MVTVTL( NU, LW, CR ) 
      ELSE IF ( INFIN(1) == 1 .AND. INFIN(2) == 2 .AND. INFIN(3) == 0 ) THEN
!  1 2 0
         U1 = UPPER(2)
         U2 = UPPER(3)
         L1 = LOWER(2)
         L2 = UPPER(3)
         TVNT = MVBVTL( NU, U1, U2, CR(3) ) - MVBVTL( NU, L1, L2, CR(3) ) 
         UP = UPPER
         UP(1) = LOWER(1)
         LW = LOWER
         LW(3) = UPPER(3)
         TVNT = TVNT - MVTVTL( NU, UP, CR ) + MVTVTL( NU, LW, CR ) 
      ELSE IF ( INFIN(1) == 2 .AND. INFIN(2) == 1 .AND. INFIN(3) == 0 ) THEN
!  2 1 0
         U1 = UPPER(1)
         U2 = UPPER(3)
         L1 = LOWER(1)
         L2 = UPPER(3)
         TVNT = MVBVTL( NU, U1, U2, CR(2) ) - MVBVTL( NU, L1, L2, CR(2) ) 
         UP = UPPER
         UP(2) = LOWER(2)
         LW = LOWER
         LW(3) = UPPER(3)
         TVNT = TVNT - MVTVTL( NU, UP, CR ) + MVTVTL( NU, LW, CR ) 
      ELSE IF ( INFIN(1) == 1 .AND. INFIN(2) == 1 .AND. INFIN(3) == 2 ) THEN
!  1 1 2
         UP = -LOWER
         LW = -LOWER
         LW(3) = -UPPER(3)
         TVNT = MVTVTL( NU, UP, CR ) - MVTVTL( NU, LW, CR ) 
      ELSE IF ( INFIN(1) == 1 .AND. INFIN(2) == 2 .AND. INFIN(3) == 1 ) THEN
!  1 2 1
         UP = -LOWER
         LW = -LOWER
         LW(2) = -UPPER(2)
         TVNT = MVTVTL( NU, UP, CR ) - MVTVTL( NU, LW, CR ) 
      ELSE IF ( INFIN(1) == 2 .AND. INFIN(2) == 1 .AND. INFIN(3) == 1 ) THEN
!  2 1 1
         UP = -LOWER
         LW = -LOWER
         LW(1) = -UPPER(1)
         TVNT = MVTVTL( NU, UP, CR ) - MVTVTL( NU, LW, CR ) 
      ELSE IF ( INFIN(1) == 1 .AND. INFIN(2) == 2 .AND. INFIN(3) == 2 ) THEN
!  1 2 2
         LT = LOWER(2:3)
         UT = UPPER(2:3)
         INFT = INFIN(2:3)
         TVNT = MVBVT( NU, LT, UT, INFT, CR(3) ) 
         UP = UPPER
         UP(1) = LOWER(1)
         LW = LOWER
         TVNT = TVNT - MVTVTL( NU, UP, CR ) - MVTVTL( NU, LW, CR ) 
         UP = LOWER
         UP(3) = UPPER(3)
         LW = LOWER
         LW(2) = UPPER(2)
         TVNT = TVNT + MVTVTL( NU, UP, CR ) + MVTVTL( NU, LW, CR ) 
      ELSE IF ( INFIN(1) == 2 .AND. INFIN(2) == 1 .AND. INFIN(3) == 2 ) THEN
!  2 1 2
         LT = (/ LOWER(1), LOWER(3) /) 
         UT = (/ UPPER(1), UPPER(3) /) 
         INFT = (/ INFIN(1), INFIN(3) /)
         TVNT = MVBVT( NU, LT, UT, INFT, CR(2) ) 
         UP = UPPER
         UP(2) = LOWER(2)
         LW = LOWER
         TVNT = TVNT - MVTVTL( NU, UP, CR ) - MVTVTL( NU, LW, CR ) 
         UP = LOWER
         UP(3) = UPPER(3)
         LW = LOWER
         LW(1) = UPPER(1)
         TVNT = TVNT + MVTVTL( NU, UP, CR ) + MVTVTL( NU, LW, CR ) 
      ELSE IF ( INFIN(1) == 2 .AND. INFIN(2) == 2 .AND. INFIN(3) == 1 ) THEN
!  2 2 1
         LT = LOWER(1:2)
         UT = UPPER(1:2)
         INFT = INFIN(1:2)
         TVNT = MVBVT( NU, LT, UT, INFT, CR(1) ) 
         UP = UPPER
         UP(3) = LOWER(3)
         LW = LOWER
         TVNT = TVNT - MVTVTL( NU, UP, CR ) - MVTVTL( NU, LW, CR ) 
         UP = LOWER
         UP(1) = UPPER(1)
         LW = LOWER
         LW(2) = UPPER(2)
         TVNT = TVNT + MVTVTL( NU, UP, CR ) + MVTVTL( NU, LW, CR ) 
      ELSE IF ( INFIN(1) == 2 .AND. INFIN(2) == 2 .AND. INFIN(3) == 2 ) THEN
!  2 2 2
         UP = UPPER
         LW = LOWER
         TVNT = MVTVTL( NU, UP, CR ) -  MVTVTL( NU, LW, CR ) 
         UP = UPPER
         UP(3) = LOWER(3)
         LW = LOWER
         LW(3) = UPPER(3)
         TVNT = TVNT - MVTVTL( NU, UP, CR ) + MVTVTL( NU, LW, CR ) 
         UP = UPPER
         UP(2) = LOWER(2)
         LW = LOWER
         LW(2) = UPPER(2)
         TVNT = TVNT - MVTVTL( NU, UP, CR ) + MVTVTL( NU, LW, CR ) 
         UP = UPPER
         UP(1) = LOWER(1)
         LW = LOWER
         LW(1) = UPPER(1)
         TVNT = TVNT - MVTVTL( NU, UP, CR ) + MVTVTL( NU, LW, CR ) 
      END IF
      TVNT = MIN( MAX( 0E0_STND, TVNT ), 1E0_STND ) 
!
      END FUNCTION MVTVT
!
      FUNCTION MVTVTL( NU, LIMIT, SIGMA ) RESULT(TVT)
!
!     A function for computing trivariate normal and t probabilities.
!     This function uses algorithms developed from the ideas 
!     described in the papers:
!       R.L. Plackett, Biometrika 41(1954), pp. 351-360.
!       Z. Drezner, Math. Comp. 62(1994), pp. 289-294.
!     Adaptive integration is done from SIGMA to singular point.
!
!     Code Author:
!       Alan Genz
!       Department of Mathematics
!       Washington State University
!       Pullman, WA 99164-3113
!       Email : alangenz@wsu.edu
!
!  Calculates the probability that X(I) < LIMIT(I), for I = 1,2,3.
!
! Parameters
!   NU     INTEGER         :: number of degrees of freedom; NU < 1 -> TVN.
!   LIMIT  REAL array of three upper integration limits.
!   SIGMA  REAL array of three correlation coefficients, SIGMA should 
!          contain the lower left portion of the correlation matrix R. 
!          SIGMA(1) = R(2,1), SIGMA(2) = R(3,1),  SIGMA(3) = R(3,2).  
!
      INTEGER,                       INTENT(IN) :: NU
      REAL(KIND=STND), DIMENSION(3), INTENT(IN) :: LIMIT(3), SIGMA(3)
      REAL(KIND=STND)                           :: TVT
!
      REAL(KIND=STND) :: B1, B2, B3, SQ21, SQ31, SQ32
      REAL(KIND=STND) :: R21, R2P, R31, R32, R3P, ASR, XS       
      REAL(KIND=STND), PARAMETER :: ONE = 1, EPS = 1E-13_STND
      REAL(KIND=STND), PARAMETER :: TWOPI = 6.2831853071795864769252_STND
!
      B1 = LIMIT(1)
      B2 = LIMIT(2)
      B3 = LIMIT(3)
      R21 = SIGMA(1)
      R31 = SIGMA(2)
      R32 = SIGMA(3)
      IF ( ABS(B1) + ABS(B2) + ABS(B3) < EPS ) THEN 
         TVT = ( 1 - ( ACOS(R21) + ACOS(R31) + ACOS(R32) )/TWOPI )/2
      ELSE IF ( NU < 1 .AND. ABS(R21) + ABS(R31) < EPS ) THEN 
         TVT = MVPHI(B1)*MVBVU( -B2, -B3, R32 )
      ELSE IF ( NU < 1 .AND. ABS(R31) + ABS(R32) < EPS ) THEN 
         TVT = MVPHI(B3)*MVBVU( -B1, -B2, R21 )
      ELSE IF ( NU < 1 .AND. ABS(R21) + ABS(R32) < EPS ) THEN 
         TVT = MVPHI(B2)*MVBVU( -B1, -B3, R31 )
      ELSE
         SQ21 = SQRT( ( 1 - R21 )*( 1 + R21 ) )
         SQ31 = SQRT( ( 1 - R31 )*( 1 + R31 ) )
         SQ32 = SQRT( ( 1 - R32 )*( 1 + R32 ) )
         R2P = R31*R32 - SIGN( SQ31*SQ32, R31*R32 - R21 )
         R3P = R21*R32 - SIGN( SQ21*SQ32, R21*R32 - R31 )
         IF ( ABS( R31 - R3P ) < ABS( R21 - R2P ) ) THEN
            R31 = R21
            R21 = SIGMA(2)
            R2P = R3P
            B2 = B3
            B3 = LIMIT(2)
         END IF
         R3P = R21*R31 - SIGN( SQ21*SQ31, R21*R31 - R32 )
         IF ( ABS( R32 - R3P ) < ABS( R21 - R2P ) ) THEN
            R32 = R21
            R21 = SIGMA(3)
            R2P = R3P
            B1 = B3
            B3 = LIMIT(1)
         END IF
         R2P = MAX( -ONE, MIN( R2P, ONE ) )
         IF ( ABS( R2P - 1 ) < EPS ) THEN 
            ASR =  TWOPI/4
            TVT = MVBVTL( NU, MIN( B1, B2 ), B3, R31 )
         ELSE IF ( ABS( R2P + 1 ) < EPS ) THEN 
            ASR = -TWOPI/4
            IF ( B1 > -B2 )  THEN
               TVT = MVBVTL( NU, B1, B3, R31 ) - MVBVTL( NU, -B2, B3, R31 )
            ELSE
               TVT = 0
            END IF
         ELSE
            ASR = ASIN(R2P)
            SQ21 = COS(ASR)
            SQ31 = SIGN( SQRT( 1 -  R31*R31  ), R32 - R2P*R31 )
            XS = ( B3*SQ21 - B2*SQ31 )/( R31*SQ21 - R2P*SQ31 ) 
            IF ( SQ31 >= 0 ) THEN
               IF ( XS >= B1 ) THEN
                  IF ( -R2P*SQ31 > -R31*SQ21 ) THEN
                     TVT = MVBVTL( NU, B1, B2, R2P )
                  ELSE
                     TVT = MVBVTL( NU, B1, B3, R31 )
                  END IF
               ELSE
                  TVT = MVBVTL( NU, XS, B2, R2P ) - MVBVTL( NU, XS, B3, R31 )
                  IF ( -R2P*SQ31 > -R31*SQ21 ) THEN
                     TVT = TVT + MVBVTL( NU, B1, B3, R31 )
                  ELSE
                     TVT = MVBVTL( NU, B1, B2, R2P ) - TVT
                  END IF
               END IF
            ELSE
               IF ( XS >= B1 ) THEN
                  TVT = 0
                  IF ( -R2P*SQ31 >= -R31*SQ21 )                               &
                       TVT = MVBVTL( NU,B1,B2,R2P ) - MVBVTL( NU,B1,-B3,-R31 )
               ELSE
                  TVT = MVBVTL( NU, XS, B2, R2P ) - MVBVTL( NU, XS,-B3,-R31 )
                  IF ( -R2P*SQ31 < -R31*SQ21 )                                &
                     TVT = MVBVTL(NU,B1,B2,R2P) - MVBVTL(NU,B1,-B3,-R31) - TVT
               END IF
            END IF
         END IF
         R3P = TVT
         IF ( ABS( R2P - R21 ) > EPS ) THEN
            NUF = NU
            BF1 = B1
            BF2 = B2
            BF3 = B3
            RF31 = R31
            RF32 = R32
            TVT = TVT + MVADON( TVTFNC, ASR, ASIN(R21), EPS )/TWOPI
         END IF
      END IF
!
      TVT = MIN( MAX( 0E0_STND, TVT ), 1E0_STND )
!
      END FUNCTION MVTVTL
!
      FUNCTION TVTFNC( X ) RESULT(TVTSFN)
!
      REAL(KIND=STND), INTENT(IN) :: X
      REAL(KIND=STND)             :: TVTSFN
      REAL(KIND=STND) :: BS, R, RR, EE, FP, HP, CD, RA, RB
      REAL(KIND=STND), PARAMETER :: PT = 1.57079632679489662E0_STND
      TVTSFN = 0 
      EE = ( PT - ABS(X) )**2
      IF ( EE < 5E-5_STND ) THEN
         R = SIGN( 1 - EE*( 1 - EE/12 )/2, X )
         RR = EE*( 1 - EE/3 )
      ELSE
         R = SIN(X)
         RR = 1 - R*R
      END IF
      BS = SIGN( BF2, R*BF2 )
      RA = BF1*RF31 + BF2*RF32
      RB = BF1*RF32 + BF2*RF31
      FP = ( BF1 - BS )**2/RR + 2*BF1*BS/( 1 + ABS(R) )
      HP = BF3*RR + R*RB - RA 
      CD = RR - RF31**2 - RF32**2 + 2*R*RF31*RF32  
      IF ( CD > 0 ) THEN
         HP = HP/SQRT( RR*CD )
         IF ( NUF < 1 ) THEN
            IF ( HP > -10 .AND. FP < 100 ) THEN
               TVTSFN = EXP( -FP/2 )
               IF ( HP < 10 ) TVTSFN = MVPHI(HP)*TVTSFN
            END IF
         ELSE
            FP = SQRT( 1 + FP/NUF )
            TVTSFN = MVSTDT( NUF, HP/FP )/FP**NUF  
         END IF
      ELSE
         IF ( HP > 0 ) THEN
            IF ( NUF < 1 ) THEN
               IF( FP < 100 ) TVTSFN = EXP( -FP/2 )
            ELSE
               TVTSFN = 1/SQRT( 1 + FP/NUF )**NUF
            END IF
         END IF
      END IF
      END FUNCTION TVTFNC
!
      FUNCTION MVADON( FUNCTN, A, B, TOL ) RESULT(FIN)
!
!     One Dimensional Globally Adaptive Integration Function
!
      REAL(KIND=STND),              INTENT(IN)  :: A, B, TOL
      INTERFACE 
         FUNCTION FUNCTN( X ) RESULT(VALUE)
            USE PRECISION_MODEL
            REAL(KIND=STND), INTENT(IN)  :: X
            REAL(KIND=STND)              :: VALUE
         END FUNCTION FUNCTN
      END INTERFACE   
      INTEGER                        :: I, IM, IP
      INTEGER,             PARAMETER :: NL = 25
      REAL(KIND=STND), DIMENSION(NL) :: EI, AI, BI, FI
      REAL(KIND=STND)                :: FIN, ERR
!
      FIN = MVKRND( A, B, FUNCTN, ERR )
      AI(1) = A
      BI(1) = B
      IP = 1
      IM = 1
      DO WHILE ( ERR > TOL .AND. IM < NL ) 
         IM = IM + 1
         BI(IM) = BI(IP)
         AI(IM) = ( AI(IP) + BI(IP) )/2
         BI(IP) = AI(IM)
         FI(IP) = MVKRND( AI(IP), BI(IP), FUNCTN, EI(IP) )
         FI(IM) = MVKRND( AI(IM), BI(IM), FUNCTN, EI(IM) )
         IP = 1
         ERR = 0
         FIN = 0
         DO I = 1, IM
            IF ( EI(I) > EI(IP) ) IP = I
            FIN = FIN + FI(I)
            ERR = ERR + EI(I)**2
         END DO
         ERR = 4*SQRT( ERR )
      END DO
!
      END FUNCTION MVADON
!
      FUNCTION MVKRND( A, B, FUNCTN, ERR ) RESULT(RESK)
!
!     Kronrod Rule
!
      REAL(KIND=STND), INTENT(IN)  :: A, B
      REAL(KIND=STND), INTENT(OUT) :: ERR
      INTERFACE 
         FUNCTION FUNCTN( X ) RESULT(VALUE)
            USE PRECISION_MODEL
            REAL(KIND=STND), INTENT(IN)  :: X
            REAL(KIND=STND)              :: VALUE
         END FUNCTION FUNCTN
      END INTERFACE   
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
      INTEGER,         PARAMETER                  :: N = 11
      REAL(KIND=STND), PARAMETER, DIMENSION(0:5)  :: WG =         &
      (/ 0.2729250867779007E0_STND, 0.5566856711617449E-1_STND,   &
         0.1255803694649048E0_STND, 0.1862902109277352E0_STND,    &
         0.2331937645919914E0_STND, 0.2628045445102478E0_STND /)
!
      REAL(KIND=STND), PARAMETER, DIMENSION(0:11) :: XGK =         &
      (/ 0E0_STND,                  0.9963696138895427E0_STND,     &
         0.9782286581460570E0_STND, 0.9416771085780681E0_STND,     &
         0.8870625997680953E0_STND, 0.8160574566562211E0_STND,     &
         0.7301520055740492E0_STND, 0.6305995201619651E0_STND,     &
         0.5190961292068118E0_STND, 0.3979441409523776E0_STND,     &
         0.2695431559523450E0_STND, 0.1361130007993617E0_STND /)
!
      REAL(KIND=STND), PARAMETER, DIMENSION(0:11) :: WGK =         &
      (/ 0.1365777947111183E0_STND,  0.9765441045961290E-2_STND,   &
         0.2715655468210443E-1_STND, 0.4582937856442671E-1_STND,   &
         0.6309742475037484E-1_STND, 0.7866457193222764E-1_STND,   &
         0.9295309859690074E-1_STND, 0.1058720744813894E0_STND,    &
         0.1167395024610472E0_STND,  0.1251587991003195E0_STND,    &
         0.1312806842298057E0_STND,  0.1351935727998845E0_STND /)
!
      INTEGER         :: J
      REAL(KIND=STND) :: T, FC, CEN, WID, RESG, RESK
!
!           Compute the 2N+1-point Kronrod approximation to
!            the integral, and estimate the absolute error.
!
      WID = ( B - A )/2
      CEN = ( B + A )/2
      FC = FUNCTN(CEN)
      RESG = FC*WG(0)
      RESK = FC*WGK(0)
      DO J = 1, N
         T = WID*XGK(J) 
         FC = FUNCTN( CEN - T ) + FUNCTN( CEN + T )
         RESK = RESK + WGK(J)*FC
         IF( MOD( J, 2 ) == 0 ) RESG = RESG + WG(J/2)*FC
      END DO
      ERR = ABS( WID*( RESK - RESG ) )
      RESK = WID*RESK
!
      END FUNCTION MVKRND
!
      FUNCTION MVTVTC( NU, L, U, INFIN, RHO ) RESULT(TVTC)
!
!     A function for computing complementary trivariate normal and t 
!       probabilities.
!
!  Parameters
!
!     NU     INTEGER degrees of freedom parameter.
!     L      REAL, array of lower integration limits.
!     U      REAL, array of upper integration limits.
!     INFIN  INTEGER, array of integration limits flags:
!            if INFIN, then MVTVTC computes
!             0 0 0   P( X>U(1), Y>U(2), Z>U(3) ),
!             1 0 0   P( X<L(1), Y>U(2), Z>U(3) ), 
!             1 1 0   P( X<L(1), Y<L(2), Z>U(3) ),
!             2 0 0   P( X>U(1), Y>U(2), Z>U(3)) + P( X<L(1), Y>U(2), Z>U(3) ),
!             2 1 0   P( X>U(1), Y<L(2), Z>U(3)) + P( X<L(1), Y<L(2), Z>U(3) ),
!             2 1 1   P( X>U(1), Y<L(2), Z<L(3)) + P( X<L(1), Y<L(2), Z<L(3) ),
!             2 2 0   P( X>U(1), Y>U(2), Z>U(3)) + P( X<L(1), Y>U(2), Z>U(3) ),
!                   + P( X>U(1), Y<L(2), Z>U(3)) + P( X<L(1), Y<L(2), Z>U(3) ),
!             2 2 2   P( X>U(1), Y>U(2), Z<L(3)) + P( X<L(1), Y>U(2), Z<L(3) ),
!                   + P( X>U(1), Y<L(2), Z<L(3)) + P( X<L(1), Y<L(2), Z<L(3) ),
!                     P( X>U(1), Y>U(2), Z>U(3)) + P( X<L(1), Y>U(2), Z>U(3) ),
!                   + P( X>U(1), Y<L(2), Z>U(3)) + P( X<L(1), Y<L(2), Z>U(3) ),
!             with other cases given as permutations of cases given.
!
!     RHO    REAL, correlation coefficient.
!
      INTEGER,                       INTENT(IN) :: NU
      REAL(KIND=STND), DIMENSION(:), INTENT(IN) :: L, U
      INTEGER,         DIMENSION(:), INTENT(IN) :: INFIN
      REAL(KIND=STND), DIMENSION(:), INTENT(IN) :: RHO
      REAL(KIND=STND)                           :: TVTC
!
!   Local Variables
!
      REAL(KIND=STND), DIMENSION(3) :: LW, UP
      REAL(KIND=STND)               :: T
      INTEGER,         DIMENSION(3) :: INF, STRT, CNTR, STEP
      INTEGER                       :: I, J, K
!
      IF ( INFIN(1) < 0 ) THEN
         TVTC = MVBVTC( NU, L(2:3), U(2:3), INFIN(2:3), RHO(3) )
      ELSE IF ( INFIN(2) < 0 ) THEN
         LW(1:2) = (/ L(1), L(3) /) 
         UP(1:2) = (/ U(1), U(3) /) 
         INF(1:2) = (/ INFIN(1), INFIN(3) /) 
         TVTC = MVBVTC( NU, LW(1:2), UP(1:2), INF(1:2), RHO(2) )
      ELSE IF ( INFIN(3) < 0 ) THEN
         TVTC = MVBVTC( NU, L(1:2), U(1:2), INFIN(1:2), RHO(1) )
      ELSE
         T = 0
         STRT = MODULO( INFIN, 2 )
         CNTR = STRT
         STEP = 2 - STRT
         DO
            DO I = 1, 3
               IF ( CNTR(I) == 0 ) THEN
                  INF(I) = 1
                  LW(I) = U(I) 
               ELSE
                  INF(I) = 0
                  UP(I) = L(I) 
               END IF
            END DO
            T = T + MVTVT( NU, LW, UP, INF, RHO )
            DO I = 1, 3
               CNTR(I) = CNTR(I) + STEP(I)
               IF ( CNTR(I) <= INFIN(I) ) EXIT
               CNTR(I) = STRT(I)
            END DO
            IF ( I == 4 ) EXIT
         END DO
         TVTC = MIN( MAX( 0E0_STND, T ), 1E0_STND )
      END IF
!
      END FUNCTION MVTVTC
!
END MODULE MVSTAT

