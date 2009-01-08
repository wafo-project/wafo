   PROGRAM MVTIN
C
C       COMPUTE PROBABILITY INTEGRAL FOR A MULTIVARIATE
C       NORMAL OR A MULTIVARIATE T WITH PRODUCT CORRELATION
C       STRUCTURE USING MVSTUD AND MVNPRD.
C
      DIMENSION A(50),B(50),BPD(50),INF(50),D(50),TMP(2)
      DATA  ZERO, ONE, EPS / 0.0, 1.0, 0.0001 /
      IERC = 0
      HNC = ZERO
C
C       START A NEW PROBLEM:-
C       ENTER N = NO. DIMENSIONS, AND NDF = DEG. OF FREEDOM
C       (FOR INFINITE D.F., ENTER NDF = 0)
C
   10 WRITE (*, 200)
      READ (*, *) N, NDF
      IF (N .LE. 0 .OR. NDF .LT. 0) STOP
      IF (N .GT. 50) THEN
         WRITE (*, 300)
         STOP
      ELSE
         XN = FLOAT(N)
      ENDIF
C
C       SPECIFY THE CORRELATION STRUCTURE:-
C       ENTER IX = 1 IF THERE IS A FIXED RHO
C          OR IX = 2 IF RHO IS DETERMINED FROM SAMPLE SIZES
C          OR IX = 3 IF RHO(I,J) = BPD(I)*BPD(J)
C          OR IX = 4 IF RHO(I,J) = 1 / (1 + SQRT(N))
C
   20 WRITE (*, 210)
      READ (*,*) IX
      IF (IX .LT. 1 .OR. IX .GT. 4) GO TO 10
      IF (IX .EQ. 1 .OR. IX .EQ. 4) THEN
         IF (IX .EQ. 4) THEN
             RHO = ONE / (ONE + SQRT(XN))
         ELSE
   30       WRITE (*, 220)
            READ (*, *) RHO
            IF (RHO .LT. ZERO .OR. RHO .GE. ONE) GO TO 30
         ENDIF
         SQTRHO = ZERO
         IF (RHO .GT. ZERO) SQTRHO = SQRT(RHO)
         DO 40 I = 1, N
            BPD(I) = SQTRHO
   40    CONTINUE
      ELSE
         IF (IX .EQ. 2) THEN
            WRITE (*, 310)
            READ (*, *) X0, (BPD(I), I = 1, N)
            DO 50 I = 1, N
               BPD(I) = ONE / SQRT(ONE + X0 / BPD(I))
   50       CONTINUE
         ELSE
               WRITE (*, 230) N
               READ (*, *) (BPD(I), I = 1, N)
         ENDIF
      ENDIF
C
C       SPECIFY MEANS
C
   80 WRITE (*, 205)
      READ (*, *) IDEL
      IF (IDEL .EQ. 0) THEN
         DO 90 I = 1, N
            D(I) = ZERO
   90    CONTINUE
      ELSE
         NCM = 1
         IF (IDEL .NE. 1) NCM = N
         WRITE (*, 270) NCM
         READ (*, *) (D(I), I = 1, NCM)
         IF (NCM .EQ. N) GO TO 120
         DO 110 I = 1, N
            D(I) = D(1)
  110    CONTINUE
      ENDIF
  120 CONTINUE
      WRITE (*, 240)
      READ (*, *) NCM
         IF (NCM .NE. 1) NCM = N
      NSTRT = 1
  125 CONTINUE
C
C       SPECIFY THE VALUES FOR THE ENDPOINTS
C
            WRITE (*, 245)
            DO 130 I = NSTRT, NCM
  132       IF (NCM .EQ. 1) THEN
                WRITE (*, 250)
            ELSE
                WRITE (*, 260) I
            END IF
            READ (*, *) INF1, TMP(1), (TMP(J), J = 2, INF1)
            IF (INF1 .LT. 0 .OR. INF1 .GT. 2) GO TO 132
            INF(I) = INF1
            IF (INF1 .EQ. 0) THEN
               B(I) = TMP(1)
               A(I) = ZERO
            ELSE
               IF (INF1 .EQ. 1) THEN
                  B(I) = ZERO
                  A(I) = TMP(1)
               ELSE
                  B(I) = TMP(1)
                  A(I) = TMP(2)
               ENDIF
            ENDIF
            IF (NCM .EQ. 1) GO TO 131
  130    CONTINUE
  131    IF (NCM .EQ. 1) THEN
            DO 140 I = 2, N
               A(I) = A(1)
               B(I) = B(1)
               INF(I) = INF(1)
  140       CONTINUE
         ELSE
            CONTINUE
         ENDIF
      NX = N
      IF (N .GT. 4) NX = 4
      WRITE (*, 180) NX
      DO 150 I = 1, NX
         WRITE (*, 190) I, B(I), A(I), INF(I), BPD(I), D(I)
  150 CONTINUE
C
C        COMPUTE PROB
C
      WRITE (*, 290) N, NDF
      CALL MVSTUD(NDF,A,B,BPD,EPS,N,INF,D,IERC,HNC,PROB,BND,IFLT)
      WRITE (*, 280) PROB, BND, IFLT
  160 WRITE (*, 170)
      READ (*, *) NSTRT
      IF (NSTRT .GT. N) GO TO 160
      IF (NSTRT .GE. 1) GO TO 125
      GO TO 10
  170 FORMAT(1X,'Enter 1 to repeat with new values for limits',
     *    /1X,'   or J for new values for J-th to N-th limits, incl.',
     *    /1X '   or 0 for new problem')
  180 FORMAT(1X,'First',I2,' parameter sets:-',
     *       /,1X,' i      L_i      U_i   Type     b_i    \mu_i')
  190 FORMAT(1X, I2, 2F9.3, I6, 2F9.3)
  200 FORMAT(1X, 'Enter No. Dimensions, and D.o.F. (0 for MVN)',
     *  /8X, '(or enter 0, 0 to terminate)')
  205 FORMAT(1X, 'Enter 0 for central case,'/
     *       1X, '   or 1 for vector (D,...,D),'/
     *       1X, '   or 2 for vector (D_1,...,D_N)')
  210 FORMAT(1X, 'Enter 1 for equal-rho case,'/
     *       1X, '   or 2 to compute \rho_ij from sample sizes,'/
     *       1X, '   or 3 if \rho_ij = b_i*b_j,'/
     *       1X, '   or 4 for \rho = 1 / (1+sqrt(N))')
  220 FORMAT(1X, 'Enter common \rho')
  230 FORMAT(1X, 'Enter b_i, i = 1, ',I3)
  240 FORMAT(1X, 'Enter 1 for a common set of limits or 0 if not')
  245 FORMAT(1X,'For each variable enter type of interval and limits,'/
     *       1X,'Type = 0 iff [a,\infty], = 1 iff [-\infty,b]',
     *', = 2 iff [a,b]'/)
  250 FORMAT(1X, 'Enter common values of TYPE and limit(s)')
  260 FORMAT(1X, 'Enter TYPE and limit(s) for I =',I3)
  270 FORMAT(1X, 'Enter',I3,' mean value(s)')
  280 FORMAT(1X, 'PROB  =',F10.6,/,1X,'BOUND =',F10.6,/,1X,
     *  'FAULT =', I10, /)
  290 FORMAT(/,1X,'Computed results for N =',I3,', NDF =',I4,':-')
  300 FORMAT(1X, 'DIMENSION LIMITS EXCEEDED',/)
  310 FORMAT(1X, 'Enter sample sizes (control first)',/)
      END
      SUBROUTINE MVSTUD(NDF,A,B,BPD,ERRB,N,INF,D,IERC,HNC,PROB,
     *    BND,IFLT)
C
C        COMPUTE MULTIVARIATE STUDENT INTEGRAL,
C        USING MVNPRD (DUNNETT, APPL. STAT., 1989)
C        IF RHO(I,J) = BPD(I)*BPD(J).
C
C        IF RHO(I,J) HAS GENERAL STRUCTURE, USE
C        MULNOR (SCHERVISH, APPL. STAT., 1984) AND REPLACE
C        CALL MVNPRD(A,B,BPD,EPS,N,INF,IERC,HNC,PROB,BND,IFLT)
C        BY CALL MULNOR(A,B,SIG,EPS,N,INF,PROB,BND,IFLT).
C
C        AUTHOR: C.W. DUNNETT, MCMASTER UNVERSITY
C
C        BASED ON ADAPTIVE SIMPSON'S RULE ALGORITHM
C        DESCRIBED IN SHAMPINE & ALLEN: "NUMERICAL
C        COMPUTING", (1974), PAGE 240.
C
C        PARAMETERS ARE SAME AS IN ALGORITHM AS 251
C        IN APPL. STAT. (1989), VOL. 38: 564-579
C        WITH THE FOLLOWING ADDITIONS:
C             NDF   INTEGER      INPUT  DEGREES OF FREEDOM
C             D     REAL ARRAY   INPUT  NON-CENTRALITY VECTOR
C        (PUT NDF = 0 FOR INFINITE D.F.)
C
      INTEGER NN
      PARAMETER (NN=50)
      DIMENSION A(*),B(*),BPD(*),INF(*),D(*),F(3),AA(NN),BB(NN)
      DATA ZERO,HALF,TWO,THREE,FOUR / 0.0, 0.5, 2.0, 3.0, 4.0 /
      DO 10 I = 1, N
         AA(I) = A(I) - D(I)
         BB(I) = B(I) - D(I)
   10 CONTINUE
      IF (NDF .LE. 0) THEN
         CALL MVNPRD(AA,BB,BPD,ERRB,N,INF,IERC,HNC,PROB,BND,IFLT)
         RETURN
      ENDIF
      BND   = ZERO
      IFLT  =  0
      MAXDF = 150
      ERB2  = ERRB
C
C        CHECK IF D.F. EXCEED MAXDF; IF YES, THEN PROB
C        IS COMPUTED BY QUADRATIC INTERPOLATION ON 1./D.F.
C
      IF (NDF .LE. MAXDF) GO TO 20
      CALL MVNPRD(AA,BB,BPD,ERB2,N,INF,IERC,HNC,F(1),BND,IFLT)
      NF    =  MAXDF / 2
      CALL SIMPSN(NF,A,B,BPD,ERB2,N,INF,D,IERC,HNC,F(3),BND,IFLT)
      NF    =  NF * 2
      CALL SIMPSN(NF,A,B,BPD,ERB2,N,INF,D,IERC,HNC,F(2),BND,IFLT)
      XX    =  FLOAT(NF) / FLOAT(NDF)
      AX    =  F(3) - F(2)*TWO + F(1)
      BX    =  F(2)*FOUR - F(3) - F(1)*THREE
      PROB  =  F(1) + XX * (AX * XX + BX) * HALF
      RETURN
   20 CALL SIMPSN (NDF,A,B,BPD,ERB2,N,INF,D,IERC,HNC,PROB,BND,IFLT)
      RETURN
      END
