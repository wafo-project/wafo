c                  3-IX-93
C

C
C  The programs RIND computes the following problem:
C
C  Let the process  y(t)=BU(t)+Delta(t), where  Delta(t)  is zero mean
C  Gaussian process and BU(t) be expectation of  y(t). Consider the process  x
C  at the grid T(1),...,T(N), N=0,...,50, (N=0 means the empty grid T).
C  Let  y(I) = BU(T(I)) + Delta(T(I)).  Observe we do not assume that the
C  points T(I) are ordered or from, e.g. T(I) are in R^K.
C
C  The covariance fkn of Delta at the points T(1),...,T(N), are given in
C  the vector R;  Cov( Delta(T(i)), Delta(T(j)) = R(i+(j-1)*N),
C  furter E[y(T(i))] = BU(i). Hence the dimension of R must be N*N=2500.
C  The vector SQ(1), ...., SQ(N) contains standard deviations of the residual
C  Delta(T(I)), e.g. SQ(I) = SQRT (R(I+(I-1)*N)). However the small values
C  of SQ could be corrupted by nummerical errors especially if the covariance
C  structure was computed using the FFT algorithm. IF R(I+(I-1)*N)<EPS
C  SQ(I)=0.    and is used as an indicator that one is not allowed to
C  condition on Delta(T(I)). Further when one have conditioned on the point
C  T(I) the variance is put to zero.
C
C  Consider in addition to y(t) a Gaussian variable  Y, E[Y]=DBUN, Var(Y)=VDER,
C  DB(I)=Cov(Delta(T(I)),Y).
C
C  *** XIND - is the result;  XIND=E[Y^+1{ HH<y(I)<0 for all I, I=1,...,N}] ***
C
C In the speccial case by choosing  DB(I)=0, VDER=0 and DBUN=1, IAC=0,1,
C (if IAC=0 VDER can take any positive value), the output XIND is equal to
C XIND=Prob( HH < y(I) < 0 for all I, I=1,...,N).
C
C
C  Some control variables:
C  INFR=1 means that both R, DB and SQ are the same as in the previous
C  call of RIND subroutine, INFR=0 indicates the new R, DB and SQ. The history
C  of the conditioning is saved in the vectors INF(5), INFO(5): INF(1), ...,
C  INF(5) are the times one has conditioned in the subroutines RINDT1,...,
C  RINDT5, respectively. After conditioning INFO(i)=INF(i). Now if INF=INFO
C  then the conditioning tree has not be changed and one not need to compute
C  the conditonal covariances. This is really time saving trick. We are assume
C  that the program saves all the time the matrices at the same fysical
C  location, e.g. the values of variables are saved during the execution.
C  This has all the time be checked when new compilator will be used.
C
C  The variable ISQ marks which type of conditioning will be used ISQ=0
C  means random time where the probability is minimum, ISQ=1 is the time
C  where the variance of the residual process is minimal(ISQ=1 is faster).
C
C  NIT, IAC are  described in CROSSPACK paper, EPS0 is the accuracy constant
C  used in choosing the number of nodes in numerical integrations
C  (XX1, H1 vectors). The nodes and weights and other parameters are
C  read in the subroutine INITINTEG from files Z.DAT, H.DAT and ACCUR.DAT.
C
C
C    NIT=0, IAC=1 then one uses RIND0 - subroutine, all other cases
C    goes through RIND1, ...,RIND5. NIT=0, here means explicite formula
C    approximation for XIND=E[Y^+1{ HH<BU(I)<0 for all I, I=1,...,N}], where
C    BU(I) is deterministic function.
C
C    NIT=1, leads tp call RIND1, IAC=0 is also explicite form approximation,
C    while IAC=1 leads to maximum one dimensional integral.
C    .......
C    NIT=5, leads tp call RIND5, IAC is maximally 4-dimensional integral,
C    while IAC=1 leads to maximum 5 dimensional integral.
C




      SUBROUTINE RIND(XIND,NIT,R,BU,DBUN,DB,SQ,VDER,EPS0,IAC,N,INFR)
      real*8 SQ
      DIMENSION R(1),BU(1),SQ(1),DB(1),INF(10),INFO(10),HH(101)
      
      COMMON /TBR/ HH
      
      COMMON /INFC/   ISQ,INF,INFO
      COMMON /CHECK1/ III01,III11,III21,III31,III41,III51
     *,III61,III71,III81,III91,III101
      COMMON /EPS/    EPS,EPSS,CEPSS

C
C      III01,III11,... - variables,counts how many times one calls
C      subroutine RIND0,RIND1,..., III*1 are also modified in the
C      subroutines RIND*. This gives us statistics over the complexity of
C      numerical calculations.
C
      XIND=0.
      IF (N.lt.1) go to 99
      
      IF (INFR.EQ.0) THEN
         NNIT=MIN(NIT,N)
         DO 1 I=1,10
            INF(I)=0
            INFO(I)=0
 1       CONTINUE
         III=0
         DO 2 I=1,N
          IF (SQ(I).GT.EPS) then
            III=1
            else
             IF(BU(I).GT.0.) THEN
             RETURN
             END IF
             IF(BU(I).LT.HH(I)) THEN
             RETURN
             END IF
         END IF
 2      CONTINUE
        END IF
        IF (III.eq.0) go to 99

      GO TO (10,20,30,40,50,60,70,80,90,100) NNIT
      if (nnit.gt.10) nnit=10
      
      CALL RIND0(XIND,BU,DBUN,VDER,N)
      iii01=iii01+1
      RETURN
 10   CALL RIND1(XIND,R,BU,DBUN,DB,SQ,VDER,EPS0,IAC,N)
      iii11=iii11+1
      RETURN
 20   CALL RIND2(XIND,R,BU,DBUN,DB,SQ,VDER,EPS0,IAC,N)
      iii21=iii21+1
      RETURN
 30   CALL RIND3(XIND,R,BU,DBUN,DB,SQ,VDER,EPS0,IAC,N)
      iii31=iii31+1
      RETURN
 40   CALL RIND4(XIND,R,BU,DBUN,DB,SQ,VDER,EPS0,IAC,N)
      iii41=iii41+1
      RETURN
 50   CALL RIND5(XIND,R,BU,DBUN,DB,SQ,VDER,EPS0,IAC,N)
      iii51=iii51+1
      RETURN
 60   CALL RIND6(XIND,R,BU,DBUN,DB,SQ,VDER,EPS0,IAC,N)
      iii61=iii61+1
      RETURN
 70   CALL RIND7(XIND,R,BU,DBUN,DB,SQ,VDER,EPS0,IAC,N)
      iii71=iii71+1
      RETURN
 80   CALL RIND8(XIND,R,BU,DBUN,DB,SQ,VDER,EPS0,IAC,N)
      iii81=iii81+1
      RETURN
 90   CALL RIND9(XIND,R,BU,DBUN,DB,SQ,VDER,EPS0,IAC,N)
      iii91=iii91+1
      RETURN
 100  CALL RIND10(XIND,R,BU,DBUN,DB,SQ,VDER,EPS0,IAC,N)
      iii101=iii101+1
      RETURN
 99   continue   
       SDER=0.
       IF(VDER.GT.EPS) SDER=SQRT(VDER)
       XIND=PMEAN(DBUN,SDER)
       return
 
       
       END

      SUBROUTINE RIND0(XIND,BU,DBUN,VDER,N)
      
      DIMENSION BU(1)
      DIMENSION HH(101)

      COMMON /EPS/ EPS,EPSS,CEPSS
      COMMON /TBR/ HH
      
      IF (N.LT.1) GO TO 20
      XIND=0.
      IF(DBUN.LT.0.) THEN
      RETURN
      END IF
      DO 10 I=1,N
      IF(BU(I).GT.0.) THEN
      RETURN
      END IF
      IF(BU(I).LT.HH(I)) THEN
      RETURN
      END IF
10    CONTINUE
20    CONTINUE
      SDER=0.
      IF(VDER.GT.EPS) SDER=SQRT(VDER)
      XIND=PMEAN(DBUN,SDER)
      RETURN
      END

      SUBROUTINE RIND1(XIND,R,BU,DBUN,DB,SQ,VDER,EPS0,IAC,N)
      real*8 XX1,H1,SQ,SQ1
      DIMENSION R(1),BU(1),SQ(1),DB(1),B1(101),SQ1(101)
      DIMENSION XX1(24),H1(24),INF(10),INFO(10)
      DIMENSION HH(101)
      
      COMMON /EPS/    EPS,EPSS,CEPSS
      COMMON /RINT/   C,FC
      COMMON /TBR/    HH
      COMMON /INFC/   ISQ,INF,INFO
      COMMON /CHECK1/ III01,III11,III21,III31,III41,III51
     *,III61,III71,III81,III91,III101
      
c      print *,'Topp of R1:',sq(1),sq(2),sq(3)
      XIND=0.

C Choice of the time for conditioning, two methods
C
C ISQ=1; INF(1)=II0 is the point where the  SQ(I) obtaines its maximum, SQ0
C is the maximal st. deviation of the residual.
C
C ISQ=0; INF(1) is the time point when the probability P(hh<BU(I)+Delta(I)<0)
C obtains its minimum which is denoted by XFF.
C
      XFF=1.
      IF (N.LT.1) GO TO 11
C
C  If N<1 means empty grid we can not condition on any point and hence GO TO 11
C  then XFF=1. and the XIND will be approximated by E[Y^+]=PMEAN(DBUN,SDER).
C  Obs. E[Y]=DBUN and Var(Y)=VDER.
C
      SQ0=0.
      DO 1 I=1,N
      sq1(i)=0.
      IF (SQ(I).LE.eps) GO TO 1
      HHB=HH(I)-BU(I)
C
C Obs. SQ(I)<=EPS idicates that the point is not good for conditioning.
C There can be two reasons for it: 1 Variance of the residual is too small
C or the point was already used before.
C
      if (ISQ.GT.1) then
      SS0=R(I+(I-1)*N)
      DB1N=DB(I)
      VDER1=DB1N*DB1N/SS0
      IF (VDER1.gt.SQ0) Then
      SQ0=VDER1
      II0=I
      end if
      else
      IF (SQ(I).GT.SQ0) THEN
      SQ0=SQ(I)
      II0=I
      END IF
      end if
      X=-BU(I)/SQ(I)
      XH=HHB/SQ(I)
      XF=FI(X)-FI(XH)
      IF(XF.LT.XFF) THEN
      INF(1)=I
      XFF=XF
      END IF
1     CONTINUE
11     CONTINUE
C
C   If the minimum probability XFF is close to 0 !!! then the indicator 1{}
C   can be bounded by EPSS leading to the approximation of XIND=0 and RETURN.
C
      IF(XFF.LT.EPSS) RETURN
C
C  We are stoping because we assume that for all sample pathes I(0,t)=1
C
C   If the minimum probability XFF is close to one then the indicator 1{}
C   can be bounded by 1 leading to the approximation of XIND by E[Y^+]=
C   PMEAN(DBUN,SDER), if IAC=1 or with E[Y]^+=MAX(0,DBUN).
C   This reduces the order of integrations.
C
      IF (XFF.GT.0.999*FC) THEN
      SDER=0.
      IF(VDER.GT.EPS) SDER=SQRT(VDER)
      XIND=MAX(DBUN,0.)
      IF(IAC.LT.1) RETURN
      XIND=PMEAN(DBUN,SDER)
      RETURN
      END IF
C
C  We are conditioning on the point T(INF(1)). If ISQ=1 INF(1)=ii0.
C  Obviously, X(INF(1))=BU(INF(1))+Delta(INF(1)), where SQ0 is a standart
C  deviation of Delta(IN(1)), hence  if X(INF(1))>XMA or X(INF(1))<XM1
C  1{}=0. Hence the values of X(INF(1)) are truncated to the interval [xmi,xma].

      IF(ISQ.EQ.1) INF(1)=II0
      SQ0=SQ(INF(1))
      XMA=-BU(INF(1))/SQ0
      XMI=XMA+HH(INF(1))/SQ0
      XMI=MAX(-C,XMI)
      XMA=MIN(C,XMA)
      IF (XMI.GT.XMA) XMA=-C
C
C   Now we are checking whether INF(I)=INFO(I), I=1,..,5, what indicates
C   that all conditional covariances and variancec are unchanged since the
C   last call of this subroutine (R,SQ,DB) as well as
C
      III=0
      DO 2 I=1,10
      III=III + ABS(INF(I)-INFO(I))
2     CONTINUE
      IF (III.EQ.0) GO TO 99
      DO 3 I=1,N
      B1(I)=R(I+(INF(1)-1)*N)
3     CONTINUE
      SS0=B1(INF(1))
      DB1N=DB(INF(1))
      INFO(1)=INF(1)
      VDER1=VDER-DB1N*DB1N/SS0
      SDER1=0.
      IF(VDER1.GT.EPS) SDER1=SQRT(VDER1)
      DB1N=DB1N/SQ0
      DO 4 I=1,N
4     B1(I)=B1(I)/SQ0
99    CONTINUE
C
C  Here conditioning is done
C
      CALL C1_C2(XMI,XMA,BU,B1,DBUN,DB1N,0.,SQ1,N)
      IF(FI(XMA)-FI(XMI).LT.EPSS) RETURN
C *******************************************************
C
C  In this special case RIND1 if IAC=0 one can explicitly compute XIND
C  and stop.
      IF (IAC.LT.1) THEN
      XIND=GAUSINT(XMI,XMA,DBUN,DB1N,1.,0.)
      RETURN
      END IF
      CALL GAUSS1(N1,H1,XX1,XMI,XMA,EPS0)
c      print *,XMI,XMA,EPS0,N1
c      write(11,*) XMI,XMA,EPS0,N1
      XIND=0.
      DO 10 J=1,N1
      DER=DBUN+XX1(J)*DB1N
      XIND=XIND+PMEAN(DER,SDER1)*H1(J)
      III01=III01+1
c      IF (N.eq.15) then
c      print *,'der,dbun,db1n,sder1',der,dbun,db1n,sder1
c      write(11,*) der,dbun,db1n,sder1
c      end if
10    CONTINUE
c      IF (N.eq.15) then
c      do 999 iii=1,N
c      print *,iii,sq(Iii)
c999   continue
c      write(11,*) XIND,INF(1),INF(2),INF(3),inf(4)
c      pause
c      end if
      RETURN
      END

      SUBROUTINE RIND2(XIND,R,BU,DBUN,DB,SQ,VDER,EPS0,IAC,N)
      real*8 XX1,H1,SQ,SQ1
      DIMENSION R(1),BU(1),SQ(1),DB(1)
      DIMENSION R1(10201),B1(101),DB1(101),BU1(101),SQ1(101)
      DIMENSION XX1(24),H1(24),INF(10),INFO(10)
      DIMENSION HH(101)
      
      COMMON /EPS/    EPS,EPSS,CEPSS
      COMMON /RINT/   C,FC
      COMMON /TBR/    HH
      COMMON /INFC/   ISQ,INF,INFO
c      COMMON/CHECK/III0,III1,III2,III3,III4
      COMMON /CHECK1/ III01,III11,III21,III31,III41,III51
     *,III61,III71,III81,III91,III101
      
c      PRINT *,'Topp of 2:',sq(1),sq(2),sq(3)
      XIND=0.
      XFF=1.
      IF (N.LT.1) GO TO 11
      SQ0=0.
      DO 1 I=1,N
         IF (SQ(I).LE.EPS) THEN
            SQ1(I)=SQ(I)
            GO TO 1
         END IF
c      CSQ=C*SQ(I)
c      IF (BU(I).LT.-CSQ.AND.BU(I).GT.HH(I)+CSQ) GO TO 1
c      IF (BU(I).LT.-CSQ.AND.BU(I).GT.HH(I)+CSQ) SQ1(I)=EPS1
C      IF (BU(I).GT. CSQ.OR.BU(I).LT.HH(I)-CSQ) RETURN
         IF (SQ(I).GT.SQ0) THEN
            SQ0=SQ(I)
            II0=I
         END IF
         X=-BU(I)/SQ(I)
         XH=X+HH(I)/SQ(I)
         XF=FI(X)-FI(XH)
c      IF(XF.GT.CEPSS) SQ1(I)=EPS1
         IF (XF.LT.XFF) THEN
            INF(2)=I
            XFF=XF
         END IF
 1    CONTINUE
11     CONTINUE
      IF(XFF.LT.EPSS) RETURN
      IF (XFF.GT.0.9999*FC) THEN
      SDER=0.
      IF (VDER.GT.EPS) SDER=SQRT(VDER)
      XIND=MAX(DBUN,0.)
      IF(IAC.LT.1) RETURN
      XIND=PMEAN(DBUN,SDER)
      RETURN
      END IF
      IF(ISQ.EQ.1) INF(2)=II0
      SQ0=SQ(INF(2))
      XMA=-BU(INF(2))/SQ0
      XMI=XMA+HH(INF(2))/SQ0
      XMI=MAX(-C,XMI)
      XMA=MIN(C,XMA)
      IF (XMI.GT.XMA) XMA=-C
C *************************************************
C
C  We are conditioning on the point T(INF(2)) and write new model
C  BU(I)+X*B1(I)+Delta1(I), I=1,...,N  (Obs. we do not use I=1,N)
C  SQ1(I) is standard deviation of Delta1 DBUN=BU'(N), DB1N=B1'(N) and X is
C  N(0,1) independent of Delta1, SDER1 is standard deviation of Delta1'(N).
C
      III=0
      DO 2 I=2,10
      III=III + ABS(INF(I)-INFO(I))
2     CONTINUE
      IF (III.EQ.0) GO TO 99
      CALL M_COND(R1,B1,DB1,R,DB,INF(2),N)
C      III1=III1+1
      SS0=B1(INF(2))
      INFO(2)=INF(2)
      DB1N=DB(INF(2))
      VDER1=VDER-DB1N*DB1N/SS0
      SDER1=0.
      IF (VDER1.GT.EPS) SDER1=SQRT(VDER1)
C      SQ0=SQRT(SS0)
      DB1N=DB1N/SQ0
      SQ1(INF(2))=0.
      DO 3 I=1,N
C
C IF (SQ1(I).EQ.EPS1) GO TO 3 - the .EQ. can not be raplaced with .LT. without
C some general changes in the strategy of SQ and SQ1 values. More exactly
C SQ can not be changed in this subroutine when for some  I  we would like to
C put SQ1(I)=EPS1 in the first loop. This SQ1 should not be changed here and
C thus we have GO TO 3 statement. Observe that the other SQ1 values are
C
      IF (SQ(I).LE.EPS) GO TO 3
C      IF (SQ1(I).EQ.EPS1.OR.SQ(I).LE.EPS1) GO TO 3
      XR1=R1(I+(I-1)*N)
c      IF(XR1.LT.0.) CALL ERROR(I,N,-1)
      SQ1(I)=0.
      IF (XR1.GT.EPS) SQ1(I)=SQRT(XR1)
3     CONTINUE
      DO 5 I=1,N
5     B1(I)=B1(I)/SQ0
99    CONTINUE
C
C ***********************************************************
C
C We shall condition on the values of X, XMI<X<XMA, but for some
C X values XIND will be zero leading to reduced accuracy. Hence we try
C to exclude them and narrow the interval [XMI,XMA]
C
c      PRINT *,'2:**  ',XMI,XMA
      CALL C1_C2(XMI,XMA,BU,B1,DBUN,DB1N,SDER1,SQ1,N)
c      PRINT *,'2:****',XMI,XMA
      IF(FI(XMA)-FI(XMI).LT.EPSS)  THEN
c         print *, 'Leaving R2: Exit 4', XIND,XIND1,VDER1
         RETURN
      ENDIF
      CALL GAUSS1(N1,H1,XX1,XMI,XMA,EPS0)
      DO 10 J=1,N1
      DO 20 I=1,N
20    BU1(I)=BU(I)+XX1(J)*B1(I)
      DER=DBUN+XX1(J)*DB1N
c      print *,'R2: before calling R1: (SQ1):',sq1(1),sq1(2),sq1(3)
      CALL RIND1(XIND1,R1,BU1,DER,DB1,SQ1,VDER1,EPS0,IAC,N)
      III11=III11+1
10    XIND=XIND+XIND1*H1(J)
c      print *, 'Leaving R2: Exit 5', XIND,XIND1,VDER1
      RETURN
      END

      SUBROUTINE RIND3(XIND,R,BU,DBUN,DB,SQ,VDER,EPS0,IAC,N)
      real*8 XX1,H1,SQ,SQ1
      DIMENSION R(1),BU(1),SQ(1),DB(1)
      DIMENSION R1(10201),B1(101),DB1(101),BU1(101),SQ1(101)
      DIMENSION XX1(24),H1(24),INF(10),INFO(10)
      DIMENSION HH(101)
      
      COMMON /EPS/    EPS,EPSS,CEPSS
      COMMON /RINT/   C,FC
      COMMON /TBR/    HH
      COMMON /INFC/   ISQ,INF,INFO
C      COMMON/CHECK/III0,III1,III2,III3,III4
      COMMON /CHECK1/ III01,III11,III21,III31,III41,III51
     *,III61,III71,III81,III91,III101
      
c      PRINT *,'Topp of 3:',sq(1),sq(2),sq(3)
      XIND=0.
      XFF=1.
      IF (N.LT.1) GO TO 11
      SQ0=0.
      DO 1 I=1,N
         IF (SQ(I).LE.EPS) THEN
            SQ1(I)=SQ(I)
            GO TO 1
         END IF
         IF (SQ(I).GT.SQ0) THEN
            SQ0=SQ(I)
            II0=I
         END IF
         X=-BU(I)/SQ(I)
         XH=X+HH(I)/SQ(I)
         XF=FI(X)-FI(XH)
         IF(XF.LT.XFF) THEN
            INF(3)=I
            XFF=XF
         END IF
 1    CONTINUE
 11   CONTINUE
      IF (XFF.LT.EPSS) RETURN
      IF (XFF.GT.0.9999*FC) THEN
         SDER=0.
         IF(VDER.GT.EPS) SDER=SQRT(VDER)
         XIND=MAX(DBUN,0.)
         IF (IAC.LT.1) RETURN
         XIND=PMEAN(DBUN,SDER)
         RETURN
      END IF
      IF(ISQ.EQ.1) INF(3)=II0
      SQ0=SQ(INF(3))
      XMA=-BU(INF(3))/SQ0
      XMI=XMA+HH(INF(3))/SQ0
      XMI=MAX(-C,XMI)
      XMA=MIN(C,XMA)
      IF (XMI.GT.XMA) XMA=-C
      III=0
      DO 2 I=3,10
         III=III + ABS(INF(I)-INFO(I))
 2    CONTINUE
      IF (III.EQ.0) GO TO 99
      CALL M_COND(R1,B1,DB1,R,DB,INF(3),N)
      SS0=B1(INF(3))
      DB1N=DB(INF(3))
      VDER1=VDER-DB1N*DB1N/SS0
      SDER1=0.
      IF (VDER1.GT.EPS) SDER1=SQRT(VDER1)
      INFO(3)=INF(3)
      DB1N=DB1N/SQ0
      SQ1(INF(3))=0.
      DO 3 I=1,N
         IF (SQ(I).LE.EPS) GO TO 3
         XR1=R1(I+(I-1)*N)
         SQ1(I)=0.
         IF (XR1.GT.EPS) SQ1(I)=SQRT(XR1)
3     CONTINUE
      DO 5 I=1,N
 5    B1(I)=B1(I)/SQ0
99    CONTINUE
c      PRINT *,'3:**  ',XMI,XMA
      CALL C1_C2(XMI,XMA,BU,B1,DBUN,DB1N,SDER1,SQ1,N)
c      PRINT *,'3:****',XMI,XMA,EPSS
      IF (FI(XMA)-FI(XMI).LT.EPSS) RETURN
      CALL GAUSS1(N1,H1,XX1,XMI,XMA,EPS0)
      DO 10 J=1,N1
         DO 20 I=1,N
20    BU1(I)=BU(I)+XX1(J)*B1(I)
      DER=DBUN+XX1(J)*DB1N
c      print *,'R3: before calling R2: (SQ1):',sq1(1),sq1(2),sq1(3)
      CALL RIND2(XIND1,R1,BU1,DER,DB1,SQ1,VDER1,EPS0,IAC,N)
      III21=III21+1
10    XIND=XIND+XIND1*H1(J)
      RETURN
      END

      SUBROUTINE RIND4(XIND,R,BU,DBUN,DB,SQ,VDER,EPS0,IAC,N)
      real*8 XX1,H1,SQ,SQ1
      DIMENSION R(1),BU(1),SQ(1),DB(1)
      DIMENSION R1(10201),B1(101),DB1(101),BU1(101),SQ1(101)
      DIMENSION XX1(24),H1(24),INF(10),INFO(10)
      DIMENSION HH(101)
      
      COMMON /EPS/    EPS,EPSS,CEPSS
      COMMON /RINT/   C,FC
      COMMON /TBR/    HH
      COMMON /INFC/   ISQ,INF,INFO
C      COMMON/CHECK/III0,III1,III2,III3,III4
      COMMON /CHECK1/ III01,III11,III21,III31,III41,III51
     *,III61,III71,III81,III91,III101
      
c      PRINT *,'Topp of 4:',SQ(1),SQ(2),SQ(3)
      XIND=0.
      XFF=1.
      IF (N.LT.1) GO TO 11
      SQ0=0.
      DO 1 I=1,N
         IF (SQ(I).LE.EPS) THEN
            SQ1(I)=SQ(I)
            GO TO 1
         END IF
         IF (SQ(I).GT.SQ0) THEN
            SQ0=SQ(I)
            II0=I
         END IF
         X=-BU(I)/SQ(I)
         XH=X+HH(I)/SQ(I)
         XF=FI(X)-FI(XH)
         IF (XF.LT.XFF) THEN
            INF(4)=I
            XFF=XF
         END IF
 1    CONTINUE
 11   CONTINUE
      IF (XFF.LT.EPSS) RETURN
      IF (XFF.GT.0.9999*FC) THEN
         SDER=0.
         IF(VDER.GT.EPS) SDER=SQRT(VDER)
         XIND=MAX(DBUN,0.)
         IF (IAC.LT.1) RETURN
         XIND=PMEAN(DBUN,SDER)
         RETURN
      END IF
      IF(ISQ.EQ.1) INF(4)=II0
      SQ0=SQ(INF(4))
      XMA=-BU(INF(4))/SQ0
      XMI=XMA+HH(INF(4))/SQ0
      XMI=MAX(-C,XMI)
      XMA=MIN(C,XMA)
      IF (XMI.GT.XMA) XMA=-C
      III=0
      DO 2 I=4,10
         III=III + ABS(INF(I)-INFO(I))
2     CONTINUE
      IF (III.EQ.0) GO TO 99
      CALL M_COND(R1,B1,DB1,R,DB,INF(4),N)
      SS0=B1(INF(4))
      DB1N=DB(INF(4))
      VDER1=VDER-DB1N*DB1N/SS0
      SDER1=0.
      IF (VDER1.GT.EPS) SDER1=SQRT(VDER1)
      INFO(4)=INF(4)
      DB1N=DB1N/SQ0
      SQ1(INF(4))=0.
      DO 3 I=1,N
         IF (SQ(I).LE.EPS) GO TO 3
         XR1=R1(I+(I-1)*N)
         SQ1(I)=0.
         IF (XR1.GT.EPS) SQ1(I)=SQRT(XR1)
3     CONTINUE
      DO 5 I=1,N
5     B1(I)=B1(I)/SQ0
99    CONTINUE
C      PRINT *,'**',XMI,XMA
      CALL C1_C2(XMI,XMA,BU,B1,DBUN,DB1N,SDER1,SQ1,N)
C      PRINT *,INF(4),XMI,XMA
      IF(FI(XMA)-FI(XMI).LT.EPSS) RETURN
      CALL GAUSS1(N1,H1,XX1,XMI,XMA,EPS0)
      DO 10 J=1,N1
      DO 20 I=1,N
20    BU1(I)=BU(I)+XX1(J)*B1(I)
      DER=DBUN+XX1(J)*DB1N
      CALL RIND3(XIND1,R1,BU1,DER,DB1,SQ1,VDER1,EPS0,IAC,N)
      III31=III31+1
10    XIND=XIND+XIND1*H1(J)
      RETURN
      END

      SUBROUTINE RIND5(XIND,R,BU,DBUN,DB,SQ,VDER,EPS0,IAC,N)
      real*8 XX1,H1,SQ,SQ1

      DIMENSION R(1),BU(1),SQ(1),DB(1)
      DIMENSION R1(10000),B1(100),DB1(100),BU1(100),SQ1(100),HH(101)
      
      DIMENSION XX1(24),H1(24),INF(10),INFO(10)
      
      COMMON /EPS/    EPS,EPSS,CEPSS
      COMMON /RINT/   C,FC
      COMMON /TBR/    HH
      COMMON /INFC/   ISQ,INF,INFO
C      COMMON/CHECK/III0,III1,III2,III3,III4
      COMMON /CHECK1/ III01,III11,III21,III31,III41,III51
     *,III61,III71,III81,III91,III101
      
      XIND=0.
      XFF=1.
      IF (N.LT.1) GO TO 11
      SQ0=0.
      DO 1 I=1,N
         IF (SQ(I).LE.EPS) THEN
            SQ1(I)=SQ(I)
            GO TO 1
         END IF
         IF (SQ(I).GT.SQ0) THEN
            SQ0=SQ(I)
            II0=I
         END IF
         X=-BU(I)/SQ(I)
         XH=X+HH(I)/SQ(I)
         XF=FI(X)-FI(XH)
         IF (XF.LT.XFF) THEN
            INF(5)=I
            XFF=XF
         END IF
 1    CONTINUE
 11   CONTINUE
      IF (XFF.LT.EPSS) RETURN
      IF (XFF.GT.0.9999*FC) THEN
         SDER=0.
         IF(VDER.GT.EPS) SDER=SQRT(VDER)
         XIND=MAX(DBUN,0.)
         IF (IAC.LT.1) RETURN
         XIND=PMEAN(DBUN,SDER)
         RETURN
      END IF
      IF(ISQ.EQ.1) INF(5)=II0
      SQ0=SQ(INF(5))
      XMA=-BU(INF(5))/SQ0
      XMI=XMA+HH(INF(5))/SQ0
      XMI=MAX(-C,XMI)
      XMA=MIN(C,XMA)
      IF (XMI.GT.XMA) XMA=-C
      III=0
      DO 2 I=5,10
         III=III + ABS(INF(I)-INFO(I))
2     CONTINUE
      IF (III.EQ.0) GO TO 99
      CALL M_COND(R1,B1,DB1,R,DB,INF(5),N)
      SS0=B1(INF(5))
      DB1N=DB(INF(5))
      VDER1=VDER-DB1N*DB1N/SS0
      SDER1=0.
      IF (VDER1.GT.EPS) SDER1=SQRT(VDER1)
      INFO(5)=INF(5)
      DB1N=DB1N/SQ0
      SQ1(INF(5))=0.
      DO 3 I=1,N
         IF (SQ(I).LE.EPS) GO TO 3
         XR1=R1(I+(I-1)*N)
         SQ1(I)=0.
         IF (XR1.GT.EPS) SQ1(I)=SQRT(XR1)
3     CONTINUE
      DO 5 I=1,N
5     B1(I)=B1(I)/SQ0
99    CONTINUE
      CALL C1_C2(XMI,XMA,BU,B1,DBUN,DB1N,SDER1,SQ1,N)
      IF(FI(XMA)-FI(XMI).LT.EPSS) RETURN
      CALL GAUSS1(N1,H1,XX1,XMI,XMA,EPS0)
      DO 10 J=1,N1
      DO 20 I=1,N
20    BU1(I)=BU(I)+XX1(J)*B1(I)
      DER=DBUN+XX1(J)*DB1N
      CALL RIND4(XIND1,R1,BU1,DER,DB1,SQ1,VDER1,EPS0,IAC,N)
      III41=III41+1
10    XIND=XIND+XIND1*H1(J)
      RETURN
      END

C      SUBROUTINE RIND5(XIND,R,BU,DBUN,DB,SQ,VDER,EPS0,IAC,N)
C      real*8 XX1,H1,SQ,SQ1
C      DIMENSION R(1),BU(1),SQ(1),DB(1)
C      DIMENSION R1(10201),B1(101),DB1(101),BU1(101),SQ1(101)
C      DIMENSION XX1(24),H1(24),INF(10),INFO(10)
C      DIMENSION HH(101)
C      
C      COMMON /EPS/    EPS,EPSS,CEPSS
C      COMMON /RINT/   C,FC
C      COMMON /TBR/    HH
C      COMMON /INFC/   ISQ,INF,INFO
CC      COMMON/CHECK/III0,III1,III2,III3,III4
C       COMMON /CHECK1/ III01,III11,III21,III31,III41,III51
C     *,III61,III71,III81,III91,III101
     
C      XIND=0.
C      XFF=1.
C      IF (N.LT.1) GO TO 11
C      SQ0=0.
C      DO 1 I=1,N
C         IF (SQ(I).LE.EPS) THEN
C            SQ1(I)=SQ(I)
C            GO TO 1
C         END IF
C         IF (SQ(I).GT.SQ0) THEN
C            SQ0=SQ(I)
C            II0=I
C         END IF
C         X=-BU(I)/SQ(I)
C         XH=X+HH(I)/SQ(I)
C         XF=FI(X)-FI(XH)
C         IF(XF.LT.XFF) THEN
C            INF(5)=I
C            XFF=XF
C         END IF
C 1    CONTINUE
C 11   CONTINUE
C      IF (XFF.LT.EPSS) RETURN
C      IF (XFF.GT.0.9999*FC) THEN
C      SDER=0.
C      IF(VDER.GT.EPS) SDER=SQRT(VDER)
C      XIND=MAX(DBUN,0.)
C      IF(IAC.LT.1) RETURN
C      XIND=PMEAN(DBUN,SDER)
C      RETURN
C      END IF
C      IF (ISQ.EQ.1) INF(5)=II0
C      SQ0=SQ(INF(5))
C      XMA=-BU(INF(5))/SQ0
C      XMI=XMA+HH(INF(5))/SQ0
C      XMI=MAX(-C,XMI)
C      XMA=MIN(C,XMA)
C      IF (XMI.GT.XMA) XMA=-C
C      III=0
C      DO 2 I=5,10
C         III=III + ABS(INF(I)-INFO(I))
C2     CONTINUE
C      IF (III.EQ.0) GO TO 99
C      CALL M_COND(R1,B1,DB1,R,DB,INF(5),N)
C      SS0=B1(INF(5))
C      DB1N=DB(INF(5))
C      VDER1=VDER-DB1N*DB1N/SS0
C      INFO(5)=INF(5)
C      DB1N=DB1N/SQ0
C      SQ1(INF(5))=0.
C      DO 3 I=1,N
C      IF (SQ(I).LE.EPS) GO TO 3
C      XR1=R1(I+(I-1)*N)
C      SQ1(I)=0.
C      IF (XR1.GT.EPS) SQ1(I)=SQRT(XR1)
C3     CONTINUE
C      DO 5 I=1,N
C5     B1(I)=B1(I)/SQ0
C99    CONTINUE
CC      PRINT *,'**',XMI,XMA
C      CALL C1_C2(XMI,XMA,BU,B1,DBUN,DB1N,SDER1,SQ1,N)
CC      PRINT *,INF(5),XMI,XMA
C      IF(FI(XMA)-FI(XMI).LT.EPSS) RETURN
C      CALL GAUSS1(N1,H1,XX1,XMI,XMA,EPS0)
C      DO 10 J=1,N1
C      DO 20 I=1,N
C20    BU1(I)=BU(I)+XX1(J)*B1(I)
C      DER=DBUN+XX1(J)*DB1N
C      CALL RIND4(XIND1,R1,BU1,DER,DB1,SQ1,VDER1,EPS0,IAC,N)
C      III41=III41+1
C10    XIND=XIND+XIND1*H1(J)
C      RETURN
C      END
C
      SUBROUTINE RIND6(XIND,R,BU,DBUN,DB,SQ,VDER,EPS0,IAC,N)
      real*8 XX1,H1,SQ,SQ1
      DIMENSION R(1),BU(1),SQ(1),DB(1)
      DIMENSION R1(10201),B1(101),DB1(101),BU1(101),SQ1(101)
      DIMENSION XX1(24),H1(24),INF(10),INFO(10)
      DIMENSION HH(101)
      
      COMMON /EPS/    EPS,EPSS,CEPSS
      COMMON /RINT/   C,FC
      COMMON /TBR/    HH
      COMMON /INFC/   ISQ,INF,INFO
      COMMON /CHECK1/ III01,III11,III21,III31,III41,III51
     *,III61,III71,III81,III91,III101
      
      XIND=0.
      XFF=1.
      IF (N.LT.1) GO TO 11
      SQ0=0.
      DO 1 I=1,N
         IF (SQ(I).LE.EPS) THEN
            SQ1(I)=SQ(I)
            GO TO 1
         END IF
c      CSQ=C*SQ(I)
c      IF (BU(I).LT.-CSQ.AND.BU(I).GT.HH(I)+CSQ) GO TO 1
c      IF (BU(I).LT.-CSQ.AND.BU(I).GT.HH(I)+CSQ) SQ1(I)=EPS1
C      IF (BU(I).GT. CSQ.OR.BU(I).LT.HH(I)-CSQ) RETURN
         IF (SQ(I).GT.SQ0) THEN
            SQ0=SQ(I)
            II0=I
         END IF
         X=-BU(I)/SQ(I)
         XH=X+HH(I)/SQ(I)
         XF=FI(X)-FI(XH)
c      IF(XF.GT.CEPSS) SQ1(I)=EPS1
         IF (XF.LT.XFF) THEN
            INF(6)=I
            XFF=XF
         END IF
 1    CONTINUE
11     CONTINUE
      IF(XFF.LT.EPSS) RETURN
      IF (XFF.GT.0.9999*FC) THEN
      SDER=0.
      IF (VDER.GT.EPS) SDER=SQRT(VDER)
      XIND=MAX(DBUN,0.)
      IF(IAC.LT.1) RETURN
      XIND=PMEAN(DBUN,SDER)
      RETURN
      END IF
      IF(ISQ.EQ.1) INF(6)=II0
      SQ0=SQ(INF(6))
      XMA=-BU(INF(6))/SQ0
      XMI=XMA+HH(INF(6))/SQ0
      XMI=MAX(-C,XMI)
      XMA=MIN(C,XMA)
      IF (XMI.GT.XMA) XMA=-C
C *************************************************
C
C  We are conditioning on the point T(INF(2)) and write new model
C  BU(I)+X*B1(I)+Delta1(I), I=1,...,N  (Obs. we do not use I=1,N)
C  SQ1(I) is standard deviation of Delta1 DBUN=BU'(N), DB1N=B1'(N) and X is
C  N(0,1) independent of Delta1, SDER1 is standard deviation of Delta1'(N).
C
      III=0
      DO 2 I=6,10
      III=III + ABS(INF(I)-INFO(I))
2     CONTINUE
      IF (III.EQ.0) GO TO 99
      CALL M_COND(R1,B1,DB1,R,DB,INF(6),N)
C      III1=III1+1
      SS0=B1(INF(6))
      INFO(6)=INF(6)
      DB1N=DB(INF(6))
      VDER1=VDER-DB1N*DB1N/SS0
      SDER1=0.
      IF (VDER1.GT.EPS) SDER1=SQRT(VDER1)
C      SQ0=SQRT(SS0)
      DB1N=DB1N/SQ0
      SQ1(INF(6))=0.
      DO 3 I=1,N
C
C IF (SQ1(I).EQ.EPS1) GO TO 3 - the .EQ. can not be raplaced with .LT. without
C some general changes in the strategy of SQ and SQ1 values. More exactly
C SQ can not be changed in this subroutine when for some  I  we would like to
C put SQ1(I)=EPS1 in the first loop. This SQ1 should not be changed here and
C thus we have GO TO 3 statement. Observe that the other SQ1 values are
C
      IF (SQ(I).LE.EPS) GO TO 3
C      IF (SQ1(I).EQ.EPS1.OR.SQ(I).LE.EPS1) GO TO 3
      XR1=R1(I+(I-1)*N)
c      IF(XR1.LT.0.) CALL ERROR(I,N,-1)
      SQ1(I)=0.
      IF (XR1.GT.EPS) SQ1(I)=SQRT(XR1)
3     CONTINUE
      DO 5 I=1,N
5     B1(I)=B1(I)/SQ0
99    CONTINUE
C
C ***********************************************************
C
C We shall condition on the values of X, XMI<X<XMA, but for some
C X values XIND will be zero leading to reduced accuracy. Hence we try
C to exclude them and narrow the interval [XMI,XMA]
C
      CALL C1_C2(XMI,XMA,BU,B1,DBUN,DB1N,SDER1,SQ1,N)
      IF(FI(XMA)-FI(XMI).LT.EPSS)  THEN
         RETURN
      ENDIF
      CALL GAUSS1(N1,H1,XX1,XMI,XMA,EPS0)
      DO 10 J=1,N1
      DO 20 I=1,N
20    BU1(I)=BU(I)+XX1(J)*B1(I)
      DER=DBUN+XX1(J)*DB1N
      CALL RIND5(XIND1,R1,BU1,DER,DB1,SQ1,VDER1,EPS0,IAC,N)
      III51=III51+1
10    XIND=XIND+XIND1*H1(J)
      RETURN
      END

      SUBROUTINE RIND7(XIND,R,BU,DBUN,DB,SQ,VDER,EPS0,IAC,N)
      real*8 XX1,H1,SQ,SQ1
      DIMENSION R(1),BU(1),SQ(1),DB(1)
      DIMENSION R1(10201),B1(101),DB1(101),BU1(101),SQ1(101)
      DIMENSION XX1(24),H1(24),INF(10),INFO(10)
      DIMENSION HH(101)
      
      COMMON /EPS/    EPS,EPSS,CEPSS
      COMMON /RINT/   C,FC
      COMMON /TBR/    HH
      COMMON /INFC/   ISQ,INF,INFO
      COMMON /CHECK1/ III01,III11,III21,III31,III41,III51
     *,III61,III71,III81,III91,III101
      
      XIND=0.
      XFF=1.
      IF (N.LT.1) GO TO 11
      SQ0=0.
      DO 1 I=1,N
         IF (SQ(I).LE.EPS) THEN
            SQ1(I)=SQ(I)
            GO TO 1
         END IF
c      CSQ=C*SQ(I)
c      IF (BU(I).LT.-CSQ.AND.BU(I).GT.HH(I)+CSQ) GO TO 1
c      IF (BU(I).LT.-CSQ.AND.BU(I).GT.HH(I)+CSQ) SQ1(I)=EPS1
C      IF (BU(I).GT. CSQ.OR.BU(I).LT.HH(I)-CSQ) RETURN
         IF (SQ(I).GT.SQ0) THEN
            SQ0=SQ(I)
            II0=I
         END IF
         X=-BU(I)/SQ(I)
         XH=X+HH(I)/SQ(I)
         XF=FI(X)-FI(XH)
c      IF(XF.GT.CEPSS) SQ1(I)=EPS1
         IF (XF.LT.XFF) THEN
            INF(7)=I
            XFF=XF
         END IF
 1    CONTINUE
11     CONTINUE
      IF(XFF.LT.EPSS) RETURN
      IF (XFF.GT.0.9999*FC) THEN
      SDER=0.
      IF (VDER.GT.EPS) SDER=SQRT(VDER)
      XIND=MAX(DBUN,0.)
      IF(IAC.LT.1) RETURN
      XIND=PMEAN(DBUN,SDER)
      RETURN
      END IF
      IF(ISQ.EQ.1) INF(7)=II0
      SQ0=SQ(INF(7))
      XMA=-BU(INF(7))/SQ0
      XMI=XMA+HH(INF(7))/SQ0
      XMI=MAX(-C,XMI)
      XMA=MIN(C,XMA)
      IF (XMI.GT.XMA) XMA=-C
C *************************************************
C
C  We are conditioning on the point T(INF(2)) and write new model
C  BU(I)+X*B1(I)+Delta1(I), I=1,...,N  (Obs. we do not use I=1,N)
C  SQ1(I) is standard deviation of Delta1 DBUN=BU'(N), DB1N=B1'(N) and X is
C  N(0,1) independent of Delta1, SDER1 is standard deviation of Delta1'(N).
C
      III=0
      DO 2 I=7,10
      III=III + ABS(INF(I)-INFO(I))
2     CONTINUE
      IF (III.EQ.0) GO TO 99
      CALL M_COND(R1,B1,DB1,R,DB,INF(7),N)
C      III1=III1+1
      SS0=B1(INF(7))
      INFO(7)=INF(7)
      DB1N=DB(INF(7))
      VDER1=VDER-DB1N*DB1N/SS0
      SDER1=0.
      IF (VDER1.GT.EPS) SDER1=SQRT(VDER1)
C      SQ0=SQRT(SS0)
      DB1N=DB1N/SQ0
      SQ1(INF(7))=0.
      DO 3 I=1,N
C
C IF (SQ1(I).EQ.EPS1) GO TO 3 - the .EQ. can not be raplaced with .LT. without
C some general changes in the strategy of SQ and SQ1 values. More exactly
C SQ can not be changed in this subroutine when for some  I  we would like to
C put SQ1(I)=EPS1 in the first loop. This SQ1 should not be changed here and
C thus we have GO TO 3 statement. Observe that the other SQ1 values are
C
      IF (SQ(I).LE.EPS) GO TO 3
C      IF (SQ1(I).EQ.EPS1.OR.SQ(I).LE.EPS1) GO TO 3
      XR1=R1(I+(I-1)*N)
c      IF(XR1.LT.0.) CALL ERROR(I,N,-1)
      SQ1(I)=0.
      IF (XR1.GT.EPS) SQ1(I)=SQRT(XR1)
3     CONTINUE
      DO 5 I=1,N
5     B1(I)=B1(I)/SQ0
99    CONTINUE
C
C ***********************************************************
C
C We shall condition on the values of X, XMI<X<XMA, but for some
C X values XIND will be zero leading to reduced accuracy. Hence we try
C to exclude them and narrow the interval [XMI,XMA]
C
      CALL C1_C2(XMI,XMA,BU,B1,DBUN,DB1N,SDER1,SQ1,N)
      IF(FI(XMA)-FI(XMI).LT.EPSS)  THEN
         RETURN
      ENDIF
      CALL GAUSS1(N1,H1,XX1,XMI,XMA,EPS0)
      DO 10 J=1,N1
      DO 20 I=1,N
20    BU1(I)=BU(I)+XX1(J)*B1(I)
      DER=DBUN+XX1(J)*DB1N
      CALL RIND6(XIND1,R1,BU1,DER,DB1,SQ1,VDER1,EPS0,IAC,N)
      III61=III61+1
10    XIND=XIND+XIND1*H1(J)
      RETURN
      END

      SUBROUTINE RIND8(XIND,R,BU,DBUN,DB,SQ,VDER,EPS0,IAC,N)
      real*8 XX1,H1,SQ,SQ1
      DIMENSION R(1),BU(1),SQ(1),DB(1)
      DIMENSION R1(10201),B1(101),DB1(101),BU1(101),SQ1(101)
      DIMENSION XX1(24),H1(24),INF(10),INFO(10)
      DIMENSION HH(101)
      
      COMMON /EPS/    EPS,EPSS,CEPSS
      COMMON /RINT/   C,FC
      COMMON /TBR/    HH
      COMMON /INFC/   ISQ,INF,INFO
      COMMON /CHECK1/ III01,III11,III21,III31,III41,III51
     *,III61,III71,III81,III91,III101
      
      XIND=0.
      XFF=1.
      IF (N.LT.1) GO TO 11
      SQ0=0.
      DO 1 I=1,N
         IF (SQ(I).LE.EPS) THEN
            SQ1(I)=SQ(I)
            GO TO 1
         END IF
c      CSQ=C*SQ(I)
c      IF (BU(I).LT.-CSQ.AND.BU(I).GT.HH(I)+CSQ) GO TO 1
c      IF (BU(I).LT.-CSQ.AND.BU(I).GT.HH(I)+CSQ) SQ1(I)=EPS1
C      IF (BU(I).GT. CSQ.OR.BU(I).LT.HH(I)-CSQ) RETURN
         IF (SQ(I).GT.SQ0) THEN
            SQ0=SQ(I)
            II0=I
         END IF
         X=-BU(I)/SQ(I)
         XH=X+HH(I)/SQ(I)
         XF=FI(X)-FI(XH)
c      IF(XF.GT.CEPSS) SQ1(I)=EPS1
         IF (XF.LT.XFF) THEN
            INF(8)=I
            XFF=XF
         END IF
 1    CONTINUE
11     CONTINUE
      IF(XFF.LT.EPSS) RETURN
      IF (XFF.GT.0.9999*FC) THEN
      SDER=0.
      IF (VDER.GT.EPS) SDER=SQRT(VDER)
      XIND=MAX(DBUN,0.)
      IF(IAC.LT.1) RETURN
      XIND=PMEAN(DBUN,SDER)
      RETURN
      END IF
      IF(ISQ.EQ.1) INF(8)=II0
      SQ0=SQ(INF(8))
      XMA=-BU(INF(8))/SQ0
      XMI=XMA+HH(INF(8))/SQ0
      XMI=MAX(-C,XMI)
      XMA=MIN(C,XMA)
      IF (XMI.GT.XMA) XMA=-C
C *************************************************
C
C  We are conditioning on the point T(INF(2)) and write new model
C  BU(I)+X*B1(I)+Delta1(I), I=1,...,N  (Obs. we do not use I=1,N)
C  SQ1(I) is standard deviation of Delta1 DBUN=BU'(N), DB1N=B1'(N) and X is
C  N(0,1) independent of Delta1, SDER1 is standard deviation of Delta1'(N).
C
      III=0
      DO 2 I=8,10
      III=III + ABS(INF(I)-INFO(I))
2     CONTINUE
      IF (III.EQ.0) GO TO 99
      CALL M_COND(R1,B1,DB1,R,DB,INF(8),N)
C      III1=III1+1
      SS0=B1(INF(8))
      INFO(8)=INF(8)
      DB1N=DB(INF(8))
      VDER1=VDER-DB1N*DB1N/SS0
      SDER1=0.
      IF (VDER1.GT.EPS) SDER1=SQRT(VDER1)
C      SQ0=SQRT(SS0)
      DB1N=DB1N/SQ0
      SQ1(INF(8))=0.
      DO 3 I=1,N
C
C IF (SQ1(I).EQ.EPS1) GO TO 3 - the .EQ. can not be raplaced with .LT. without
C some general changes in the strategy of SQ and SQ1 values. More exactly
C SQ can not be changed in this subroutine when for some  I  we would like to
C put SQ1(I)=EPS1 in the first loop. This SQ1 should not be changed here and
C thus we have GO TO 3 statement. Observe that the other SQ1 values are
C
      IF (SQ(I).LE.EPS) GO TO 3
C      IF (SQ1(I).EQ.EPS1.OR.SQ(I).LE.EPS1) GO TO 3
      XR1=R1(I+(I-1)*N)
c      IF(XR1.LT.0.) CALL ERROR(I,N,-1)
      SQ1(I)=0.
      IF (XR1.GT.EPS) SQ1(I)=SQRT(XR1)
3     CONTINUE
      DO 5 I=1,N
5     B1(I)=B1(I)/SQ0
99    CONTINUE
C
C ***********************************************************
C
C We shall condition on the values of X, XMI<X<XMA, but for some
C X values XIND will be zero leading to reduced accuracy. Hence we try
C to exclude them and narrow the interval [XMI,XMA]
C
      CALL C1_C2(XMI,XMA,BU,B1,DBUN,DB1N,SDER1,SQ1,N)
      IF(FI(XMA)-FI(XMI).LT.EPSS)  THEN
         RETURN
      ENDIF
      CALL GAUSS1(N1,H1,XX1,XMI,XMA,EPS0)
      DO 10 J=1,N1
      DO 20 I=1,N
20    BU1(I)=BU(I)+XX1(J)*B1(I)
      DER=DBUN+XX1(J)*DB1N
      CALL RIND7(XIND1,R1,BU1,DER,DB1,SQ1,VDER1,EPS0,IAC,N)
      III71=III71+1
10    XIND=XIND+XIND1*H1(J)
      RETURN
      END

      SUBROUTINE RIND9(XIND,R,BU,DBUN,DB,SQ,VDER,EPS0,IAC,N)
      real*8 XX1,H1,SQ,SQ1
      DIMENSION R(1),BU(1),SQ(1),DB(1)
      DIMENSION R1(10201),B1(101),DB1(101),BU1(101),SQ1(101)
      DIMENSION XX1(24),H1(24),INF(10),INFO(10)
      DIMENSION HH(101)
      
      COMMON /EPS/    EPS,EPSS,CEPSS
      COMMON /RINT/   C,FC
      COMMON /TBR/    HH
      COMMON /INFC/   ISQ,INF,INFO
      COMMON /CHECK1/ III01,III11,III21,III31,III41,III51
     *,III61,III71,III81,III91,III101
      
      XIND=0.
      XFF=1.
      IF (N.LT.1) GO TO 11
      SQ0=0.
      DO 1 I=1,N
         IF (SQ(I).LE.EPS) THEN
            SQ1(I)=SQ(I)
            GO TO 1
         END IF
c      CSQ=C*SQ(I)
c      IF (BU(I).LT.-CSQ.AND.BU(I).GT.HH(I)+CSQ) GO TO 1
c      IF (BU(I).LT.-CSQ.AND.BU(I).GT.HH(I)+CSQ) SQ1(I)=EPS1
C      IF (BU(I).GT. CSQ.OR.BU(I).LT.HH(I)-CSQ) RETURN
         IF (SQ(I).GT.SQ0) THEN
            SQ0=SQ(I)
            II0=I
         END IF
         X=-BU(I)/SQ(I)
         XH=X+HH(I)/SQ(I)
         XF=FI(X)-FI(XH)
c      IF(XF.GT.CEPSS) SQ1(I)=EPS1
         IF (XF.LT.XFF) THEN
            INF(9)=I
            XFF=XF
         END IF
 1    CONTINUE
11     CONTINUE
      IF(XFF.LT.EPSS) RETURN
      IF (XFF.GT.0.9999*FC) THEN
      SDER=0.
      IF (VDER.GT.EPS) SDER=SQRT(VDER)
      XIND=MAX(DBUN,0.)
      IF(IAC.LT.1) RETURN
      XIND=PMEAN(DBUN,SDER)
      RETURN
      END IF
      IF(ISQ.EQ.1) INF(9)=II0
      SQ0=SQ(INF(9))
      XMA=-BU(INF(9))/SQ0
      XMI=XMA+HH(INF(9))/SQ0
      XMI=MAX(-C,XMI)
      XMA=MIN(C,XMA)
      IF (XMI.GT.XMA) XMA=-C
C *************************************************
C
C  We are conditioning on the point T(INF(2)) and write new model
C  BU(I)+X*B1(I)+Delta1(I), I=1,...,N  (Obs. we do not use I=1,N)
C  SQ1(I) is standard deviation of Delta1 DBUN=BU'(N), DB1N=B1'(N) and X is
C  N(0,1) independent of Delta1, SDER1 is standard deviation of Delta1'(N).
C
      III=0
      DO 2 I=9,10
      III=III + ABS(INF(I)-INFO(I))
2     CONTINUE
      IF (III.EQ.0) GO TO 99
      CALL M_COND(R1,B1,DB1,R,DB,INF(9),N)
C      III1=III1+1
      SS0=B1(INF(9))
      INFO(9)=INF(9)
      DB1N=DB(INF(9))
      VDER1=VDER-DB1N*DB1N/SS0
      SDER1=0.
      IF (VDER1.GT.EPS) SDER1=SQRT(VDER1)
C      SQ0=SQRT(SS0)
      DB1N=DB1N/SQ0
      SQ1(INF(9))=0.
      DO 3 I=1,N
C
C IF (SQ1(I).EQ.EPS1) GO TO 3 - the .EQ. can not be raplaced with .LT. without
C some general changes in the strategy of SQ and SQ1 values. More exactly
C SQ can not be changed in this subroutine when for some  I  we would like to
C put SQ1(I)=EPS1 in the first loop. This SQ1 should not be changed here and
C thus we have GO TO 3 statement. Observe that the other SQ1 values are
C
      IF (SQ(I).LE.EPS) GO TO 3
C      IF (SQ1(I).EQ.EPS1.OR.SQ(I).LE.EPS1) GO TO 3
      XR1=R1(I+(I-1)*N)
c      IF(XR1.LT.0.) CALL ERROR(I,N,-1)
      SQ1(I)=0.
      IF (XR1.GT.EPS) SQ1(I)=SQRT(XR1)
3     CONTINUE
      DO 5 I=1,N
5     B1(I)=B1(I)/SQ0
99    CONTINUE
C
C ***********************************************************
C
C We shall condition on the values of X, XMI<X<XMA, but for some
C X values XIND will be zero leading to reduced accuracy. Hence we try
C to exclude them and narrow the interval [XMI,XMA]
C
      CALL C1_C2(XMI,XMA,BU,B1,DBUN,DB1N,SDER1,SQ1,N)
      IF(FI(XMA)-FI(XMI).LT.EPSS)  THEN
         RETURN
      ENDIF
      CALL GAUSS1(N1,H1,XX1,XMI,XMA,EPS0)
      DO 10 J=1,N1
      DO 20 I=1,N
20    BU1(I)=BU(I)+XX1(J)*B1(I)
      DER=DBUN+XX1(J)*DB1N
      CALL RIND8(XIND1,R1,BU1,DER,DB1,SQ1,VDER1,EPS0,IAC,N)
      III81=III81+1
10    XIND=XIND+XIND1*H1(J)
      RETURN
      END

      SUBROUTINE RIND10(XIND,R,BU,DBUN,DB,SQ,VDER,EPS0,IAC,N)
      real*8 XX1,H1,SQ,SQ1
      DIMENSION R(1),BU(1),SQ(1),DB(1)
      DIMENSION R1(10201),B1(101),DB1(101),BU1(101),SQ1(101)
      DIMENSION XX1(24),H1(24),INF(10),INFO(10)
      DIMENSION HH(101)
     
      COMMON /EPS/    EPS,EPSS,CEPSS
      COMMON /RINT/   C,FC
      COMMON /TBR/    HH
      COMMON /INFC/   ISQ,INF,INFO
      COMMON /CHECK1/ III01,III11,III21,III31,III41,III51
     *,III61,III71,III81,III91,III101

      XIND=0.
      XFF=1.
      IF (N.LT.1) GO TO 11
      SQ0=0.
      DO 1 I=1,N
         IF (SQ(I).LE.EPS) THEN
            SQ1(I)=SQ(I)
            GO TO 1
         END IF
c      CSQ=C*SQ(I)
c      IF (BU(I).LT.-CSQ.AND.BU(I).GT.HH(I)+CSQ) GO TO 1
c      IF (BU(I).LT.-CSQ.AND.BU(I).GT.HH(I)+CSQ) SQ1(I)=EPS1
C      IF (BU(I).GT. CSQ.OR.BU(I).LT.HH(I)-CSQ) RETURN
         IF (SQ(I).GT.SQ0) THEN
            SQ0=SQ(I)
            II0=I
         END IF
         X=-BU(I)/SQ(I)
         XH=X+HH(I)/SQ(I)
         XF=FI(X)-FI(XH)
c      IF(XF.GT.CEPSS) SQ1(I)=EPS1
         IF (XF.LT.XFF) THEN
            INF(10)=I
            XFF=XF
         END IF
 1    CONTINUE
11     CONTINUE
      IF(XFF.LT.EPSS) RETURN
      IF (XFF.GT.0.9999*FC) THEN
      SDER=0.
      IF (VDER.GT.EPS) SDER=SQRT(VDER)
      XIND=MAX(DBUN,0.)
      IF(IAC.LT.1) RETURN
      XIND=PMEAN(DBUN,SDER)
      RETURN
      END IF
      IF(ISQ.EQ.1) INF(10)=II0
      SQ0=SQ(INF(10))
      XMA=-BU(INF(10))/SQ0
      XMI=XMA+HH(INF(10))/SQ0
      XMI=MAX(-C,XMI)
      XMA=MIN(C,XMA)
      IF (XMI.GT.XMA) XMA=-C
C *************************************************
C
C  We are conditioning on the point T(INF(2)) and write new model
C  BU(I)+X*B1(I)+Delta1(I), I=1,...,N  (Obs. we do not use I=1,N)
C  SQ1(I) is standard deviation of Delta1 DBUN=BU'(N), DB1N=B1'(N) and X is
C  N(0,1) independent of Delta1, SDER1 is standard deviation of Delta1'(N).
C
      III=0
      DO 2 I=10,10
      III=III + ABS(INF(I)-INFO(I))
2     CONTINUE
      IF (III.EQ.0) GO TO 99
      CALL M_COND(R1,B1,DB1,R,DB,INF(10),N)
C      III1=III1+1
      SS0=B1(INF(10))
      INFO(10)=INF(10)
      DB1N=DB(INF(10))
      VDER1=VDER-DB1N*DB1N/SS0
      SDER1=0.
      IF (VDER1.GT.EPS) SDER1=SQRT(VDER1)
C      SQ0=SQRT(SS0)
      DB1N=DB1N/SQ0
      SQ1(INF(10))=0.
      DO 3 I=1,N
C
C IF (SQ1(I).EQ.EPS1) GO TO 3 - the .EQ. can not be raplaced with .LT. without
C some general changes in the strategy of SQ and SQ1 values. More exactly
C SQ can not be changed in this subroutine when for some  I  we would like to
C put SQ1(I)=EPS1 in the first loop. This SQ1 should not be changed here and
C thus we have GO TO 3 statement. Observe that the other SQ1 values are
cc
      IF (SQ(I).LE.EPS) GO TO 3
C      IF (SQ1(I).EQ.EPS1.OR.SQ(I).LE.EPS1) GO TO 3
      XR1=R1(I+(I-1)*N)
c      IF(XR1.LT.0.) CALL ERROR(I,N,-1)
      SQ1(I)=0.
      IF (XR1.GT.EPS) SQ1(I)=SQRT(XR1)
3     CONTINUE
      DO 5 I=1,N
5     B1(I)=B1(I)/SQ0
99    CONTINUE
C
C ***********************************************************
C
C We shall condition on the values of X, XMI<X<XMA, but for some
C X values XIND will be zero leading to reduced accuracy. Hence we try
C to exclude them and narrow the interval [XMI,XMA]
C
      CALL C1_C2(XMI,XMA,BU,B1,DBUN,DB1N,SDER1,SQ1,N)
      IF(FI(XMA)-FI(XMI).LT.EPSS)  THEN
        RETURN
      ENDIF
      CALL GAUSS1(N1,H1,XX1,XMI,XMA,EPS0)
      DO 10 J=1,N1
      DO 20 I=1,N
20    BU1(I)=BU(I)+XX1(J)*B1(I)
      DER=DBUN+XX1(J)*DB1N
      CALL RIND9(XIND1,R1,BU1,DER,DB1,SQ1,VDER1,EPS0,IAC,N)
      III91=III91+1
10    XIND=XIND+XIND1*H1(J)
      RETURN
      END
      
      SUBROUTINE C1_C2(C1,C2,BU,B1,DBUN,DB1N,SDER,SQ,N)
C
C  We assume that the process  y  is of form y(I)=BU(I)+X*B1(I)+Delta(I),
C  I=1,...,N, SQ(I) is standard deviation of Delta(I), where X is  N(0,1)
C  independent of Delta.  Let Y = DBUN + DB1N*X + Z, where Z is zero-mean
C  Gaussian with standart independent of X (it can depend on Delta(I)) with
C  standart deviation SDER. Since we are truncating all Gaussian  variables to
C  the interval [-C,C], then if for any I
C
C  a) BU(I)+x*B1(I)-C*SQ(I)>0  or
C
C  b) BU(I)+x*B1(I)+C*SQ(I)<HH  then
C
C  XIND|X=x = E[Y^+1{ HH<y(I)<0 for all I, I=1,...,N}|X=x] = 0 !!!!!!!!!
C
C  Further, see discussion in comments to the subroutine PMEAN, by first upper-
C  bounding the indicator 1{} in XIND by 1, XIND|X=x = 0 if
C
C  c) DBUN+x*DB1N+4.5*SDER<0.
C
C  Consequently, for increasing the accuracy (by excluding possible discon-
C  tinuouities) we shall exclude such X=x values for which XIND|X=x = 0.
C  XIND=E([XIND|X]). Hence we assume that if C1<X<C2 any of the previous
C  conditions are satisfied.
C
C  OBSERVE!!, C1, C2 has to be set to upper bounds of possible values, e.g.
C  C1=-C, C2=C before calling C1_C2 subroutine.
C
      real*8 SQ
      DIMENSION BU(1),B1(1),SQ(1)
      DIMENSION HH(101)
      COMMON /EPS/EPS,EPSS,CEPSS
      COMMON/RINT/C,FC
      COMMON/TBR/HH
      DO 10 I=1,N
      CSQ=C*SQ(I)
      HHB=HH(I)-BU(I)
C
C  If ABS(B1(I)) < EPS we can have overflow and hence we consider two cases
C  1) BU(I) is so large or small so we can surely assume that the probability
C     of staying between the barriers is 0, consequently C1=C2=0
C  2) we do not change the original limits.
C
      IF (ABS(B1(I)).LT. EPS) THEN
         IF (BU(I).GT.CSQ.OR.BU(I).LT.HH(I)-CSQ) THEN
            C1=0.
            C2=0.
            RETURN
         END IF
C
C  In other cases this part follows from the description of the problem.
C
       ELSE
         IF (B1(I).GT.EPS) THEN
            CC1=(HHB-CSQ)/B1(I)
            CC2=(-BU(I)+CSQ)/B1(I)
            IF (C1.LT.CC1) C1=CC1
            IF (C2.GT.CC2) C2=CC2
          ELSE
            CC2=(HHB-CSQ)/B1(I)
            CC1=(-BU(I)+CSQ)/B1(I)
            IF (C1.LT.CC1) C1=CC1
            IF (C2.GT.CC2) C2=CC2
          END IF
      END IF
10    CONTINUE
      X=-DBUN-4.5*SDER
      IF(DB1N.GT.EPS.AND.C1.LT.X/DB1N) C1=X/DB1N
      IF(DB1N.LT.-EPS.AND.C2.GT.X/DB1N) C2=X/DB1N
      if(abs(db1n).lt.eps.and.x.gt.0.) then
            C1=0.
            C2=0.
            RETURN
         END IF

c
c  In the following three rows we are cutting C1, C2 to the interval [-C,C].
c  Obs. all tree lines are neccessary.
c
      C1=MAX(-C,C1)
      C2=MIN( C,C2)
      IF (C1.GT.C2) C2=-C
C      PRINT *,2,C1,C2
      RETURN
      END

      FUNCTION GAUSINT(X1,X2,A,B,C,D)
C
C Let  X  be stardized Gaussian variable, i.e. X=N(0,1).
C The function calculate the followin integral E[I(X1<X<X2)(A+BX)(C+DX)],
C where I(X1<X<X2) is an indicator function of the set {X1<X<X2}.
C
      DATA SP/0.398942/
      IF(X1.GE.X2) THEN
      GAUSINT=0.
      RETURN
      END IF
      Y1=(A*D+B*C+X1*B*D)*EXP(-0.5*X1*X1)
      Y2=(A*D+B*C+X2*B*D)*EXP(-0.5*X2*X2)
      Y3=(A*C+B*D)*(FI(X2)-FI(X1))
      GAUSINT=Y3+SP*(Y1-Y2)
      RETURN
      END


      SUBROUTINE GAUSS1(N,H1,XX1,XMI,XMA,EPS0) 
      REAL*8  Z,H,Z1,XX1,H1
      DIMENSION Z(126),H(126),Z1(24),XX1(1),H1(1)
      DIMENSION NN(25)
      COMMON/QUADR/ Z,H,NN,NNW
      COMMON/CHECKQ/ III
      DATA SP/0.398942/
      IF (XMA.LT.XMI) THEN
      PRINT *,'Error XMIN>XMAX in GAUSS1 - stop!'
C     STOP
      END IF
      NNN=0
      DO 10 I=1,NNW
      N=NN(I)
      DO 20 J=1,N
      XX1(J)=0.5*(Z(NNN+J)*(XMA-XMI)+XMA+XMI)
      Z1(J)=XX1(J)*XX1(J)
20    H1(J)=0.5*SP*(XMA-XMI)*H(NNN+J)*EXP(-0.5*Z1(J))
      NNN=NNN+N
      SDOT=GAUSINT(XMI,XMA,0.,1.,0.,1.)
      SDOT1=0.
      DO 30 I1=1,N
30    SDOT1=SDOT1+Z1(I1)*H1(I1)
      DIFF1=ABS(SDOT-SDOT1)
      IF(EPS0.LT.DIFF1) GO TO 10
      III=III+N
C      PRINT *,'N. of nodes',iii
      RETURN
10    CONTINUE
      END

      
      SUBROUTINE M_COND(Syy_cd,Syyii,Syx_cd,Syy,Syx,ii,N)
C
C               INPUT:
C
C   ii     IS THE INDEX OF THE TIME ON WHICH WE ARE CONDITIONING.
C   N      number of variables in covariance matrix  Syy
C
C   Covariance matrix  Syy(I+(J-1)*N)=Cov(Yi,Yj) (is unchanged)
C   Covariance vector  Syx(I)=Cov(Yi,X) (is unchanged)
C
C              OUTPUT:
C
C   Covariance matrix Syy_cd(I+(J-1)*N)=Cov(Xi,Xj|Xii)
C   Covariance vector Syyii(I)=Cov(Xi,Xii)
C   Covariance vector Syx_cd(I)=Cov(Xi,Y|Xii)
C   Variance          Q1=Var(Xii)=Syyii(ii)
c   Obs. If   Q1<EPS there is no conditioning
C
C
      DIMENSION Syy_cd(1),Syyii(1),Syx_cd(1),Syy(1),Syx(1)
      COMMON /EPS/    EPS,EPSS,CEPSS
      IF (II.LE.0.OR.II.GT.N) THEN
       PRINT *,'The conditioning time in M_COND is out of range, stop!'
       STOP
      END IF
C
C   Q1=Var(Xii)=Syyii(ii)
C
      Q1=Syy(II+(II-1)*N)
      IF(Q1.LE.eps) then
      DO 5 I=1,N
      Syyii(I)=0.
 5    CONTINUE
      Q1=1.
      else
      DO 10 I=1,N
      Syyii(I)=Syy(I+(II-1)*N)
10    CONTINUE
      end if
      DO 20 I=1,N
      DO 30 J=1,N
      Syy_cd(I+(J-1)*N)=Syy(I+(J-1)*N)-Syy(II+(J-1)*N)*Syyii(I)/Q1
30    CONTINUE
      Syx_cd(I)=Syx(I)-Syx(II)*Syyii(I)/Q1
20    CONTINUE
      RETURN
      END
      


      FUNCTION PMEAN(XX,SS)
C
C   PMEAN is the positive mean of a Gaussian variable with mean  XX  and
C   standart deviation  SS,  i.e.  PMEAN=SS*FIFUNK(XX/SS), where
C   FIFUNK(x)=f(x)+x*FI(x), f  and  FI  are density and distribution
C   functions of  N(0,1)  variable, respectively. We have modified the
C   algorithm 209 from CACAM for evaluation of  FIFUNK, to avoid the operations
C   of type  SS*XX/SS, which can give numerical errors when  SS  and  XX  are
C   both small.
C
C                     **   NUMERICAL ACCURACY  **
C
C
C  Obs. that, our general assumption is that the process is normalized, i.e.
C  Var (X(t)) = Var (X'(t)) = 1.  Consequently all conditional variances
C  are less than  1, and usualy close to zero, i.e.  SS<1.. Now since  SS<1.
C  and  FIFUNK(x)<0.0000001  for  x<-4.5  we have defined  FIFUNK(x)=0.
C  if  x<-4.5  and  FIFUNK(x)=x  if  x>
C  4.5.  Under we have a table with
C  exact values of  FIFUNK.
C
C    x       FIFUNK(x)
C
C   -5.0    0.00000005233
C   -4.5    0.00000069515
C   -4.0    0.00000711075
C    4.0    4.00000700000
C    4.5    4.50000100000
C
C Obviously the tresholds  -4.5  and  4.5  can be increased.
C
C
      DATA SP/0.398942/
      IF(XX.LT.4.5*SS) GO TO 1
      PMEAN=XX
      RETURN
1     IF(XX.GT.-4.5*SS) GO TO 3
      PMEAN=0.
      RETURN
3     continue
      if (SS .LT. 0.0000001) then
      PMEAN=0.
      RETURN
      end if

      X=XX/SS
      
      IF(X)7,8,7
7     Y=0.5*ABS(X)
      IF(Y-1.) 5,6,6
5     W=Y*Y
      Z=((((((((0.000124818987*W-0.001075204047)*W+0.005198775019)*W
     1         -0.019198292)*W+0.05905403564)*W -0.15196875136)*W
     2         +0.3191529327)*W-0.5319230073)*W+0.7978845606)*Y*2.
      GO TO 100
6     Y=Y-2.
      Z=(((((((((((((-0.000045255659*Y+0.00015252929)*Y-0.000019538132)*
     1Y-0.000676904986)*Y+0.001390604284)*Y-0.000794620820)*Y
     2-0.002034254874)*Y+0.006549791214)*Y-0.010557625006)*Y+0.011630447
     3319)*Y-0.009279453341)*Y+0.005353579108)*Y-0.002141268741)*Y
     4+0.000535310849)*Y+0.9999366575
100   IF(X.GT.0.) PMEAN=SS*SP*EXP(-0.5*X*X)+XX*0.5*(Z+1.)
      IF(X.LT.0.) PMEAN=SS*SP*EXP(-0.5*X*X)+XX*0.5*(1.-Z)
      RETURN
8     PMEAN=SS*SP
      RETURN
      END


      FUNCTION FI(XX)
C
C   Algorithm 209 from CACAM.
C   FI(xx)  is a distribution functions of  N(0,1)  variable.
C
      X=XX
      IF(X) 1,2,1
2     FI=0.5
      RETURN
1     Y=0.5*ABS(X)
      IF(Y-3.) 3,4,4
4     IF(X.GT.0.) FI=1.
      IF(X.LT.0.) FI=0.
      RETURN
3     IF(Y-1.) 5,6,6
5     W=Y*Y
      Z=((((((((0.000124818987*W-0.001075204047)*W+0.005198775019)*W
     1         -0.019198292)*W+0.05905403564)*W -0.15196875136)*W
     2         +0.3191529327)*W-0.5319230073)*W+0.7978845606)*Y*2.
      GO TO 100
6     Y=Y-2.
      Z=(((((((((((((-0.000045255659*Y+0.00015252929)*Y-0.000019538132)*
     1Y-0.000676904986)*Y+0.001390604284)*Y-0.000794620820)*Y
     2-0.002034254874)*Y+0.006549791214)*Y-0.010557625006)*Y+0.011630447
     3319)*Y-0.009279453341)*Y+0.005353579108)*Y-0.002141268741)*Y
     4+0.000535310849)*Y+0.9999366575
100   IF(X.GT.0.) FI=0.5*(Z+1.)
      IF(X.LT.0.) FI=0.5*(1.-Z)
      RETURN
      END



      SUBROUTINE TRANSF(N,T,A,TIME,VALUE,DER) 
C
C N number of data points
C TIME vector of time points
C A a vector of values of a function G(TIME)
C T independent time point
C VALUE is a value of a function at T, i.e. VALUE=G(T).
c DER=G'(t)
C
      DIMENSION A(1),TIME(1)
      
      IF (T.LT.TIME(1))  then
      der=(A(2)-A(1))/(Time(2)-TIME(1))
      T1=T-TIME(1)
      VALUE=A(1)+T1*DER
      return
      end if
      IF (T.GT.TIME(N)) then
      der=(A(N)-A(N-1))/(Time(N)-TIME(N-1))
      T1=T-TIME(N)
      VALUE=A(N)+T1*DER
      return
      end if
      DO 5 I=2,N
      IF (T.LT.TIME(I)) GO TO 10
5     CONTINUE
10    I=I-1
      T1=T-TIME(I)
      DER=(A(I+1)-A(I))/(Time(i+1)-TIME(I))
      VALUE=A(I)+T1*DER
      RETURN
      END
      
      FUNCTION SPLE(N,T,A,TIME)
C
C N number of data points
C TIME vector of time points
C A a vector of values of a function G(TIME)
C T independent time point
C SPLE is a value of a function at T, i.e. SPLE=G(T).
C
      DIMENSION A(1),TIME(1)
      SPLE=-9.9
      IF (T.LT.TIME(1) .OR. T.GT.TIME(N)) RETURN
      DO 5 I=2,N
      IF (T.LT.TIME(I)) GO TO 10
5     CONTINUE
10    I=I-1
      T1=T-TIME(I)
      SPLE=A(I)+T1*(A(I+1)-A(I))/(Time(i+1)-TIME(I))
      RETURN
      END



       SUBROUTINE SVBKSB(U,W,V,M,N,MP,NP,B,X)
C
C   Solves  AX=B  for a vector  X, where  A  is specified by the arrays
C   U, W, V  as returned by SVDCMP.  M  and  N  are the logical
C   dimensions of  A, and will be equal for a square matrices.  MP  and  NP
C   are the phisical dimensions of  A.  B  is the input right-hand side.
C   X  is the output solution vector. No input quantities are destroyed,
C   so the routine may be called sequentialy with different  B's.
C
       PARAMETER (NMAX=100)
C   Maximum anticipated value of N
       DIMENSION U(MP,NP),W(NP),V(NP,NP),B(MP),X(NP),TMP(NMAX)
       DO 12 J=1,N
C   Cumulate U^T*B
         S=0.
         IF (W(J).NE.0.) THEN
C   Nonzero rezult only if  wj  is nonzero
           DO 11 I=1,M
             S=S+U(I,J)*B(I)
11         CONTINUE
           S=S/W(J)
C   This is the divide by  wj
         ENDIF
         TMP(J)=S
12     CONTINUE
       DO 14 J=1,N
         S=0.0
         DO 13 JJ=1,N
           S=S+V(J,JJ)*TMP(JJ)
13       CONTINUE
         X(J)=S
14     CONTINUE
       RETURN
       END



       SUBROUTINE SVDCMP(A,M,N,MP,NP,W,V)
C
C  Given a matrix  A, with logical dimensions  M  by  N  and physical
C  dimensions  MP  by  NP, this routine computes its singular value
C  decomposition,  A=U.W.V^T, see Numerical Recipes, by Press W.,H.
C  Flannery, B. P., Teukolsky S.A. and Vetterling W., T. Cambrige
C  University Press 1986, Chapter 2.9. The matrix  U  replaces A  on
C  output. The diagonal matrix of singular values  W  is ouyput as a vector
C  W. The matrix  V (not the transpose  V^T) is output as  V.  M  must be
C  greater or equal to  N; if it is smaller, then  A  should be filled up
C  to square with zero rows.
C
       PARAMETER (NMAX=100)
C  Maximum anticipated values of  N
       DIMENSION A(MP,NP),W(NP),V(NP,NP),RV1(NMAX)
       IF(M.LT.N) PAUSE 'You must augment  A  with extra zero rows.'
C  Householder reduction to bidiagonal form
       G=0.0
       SCALE=0.0
       ANORM=0.0
       DO 25 I=1,N
          L=I+1
          RV1(I)=SCALE*G
          G=0.0
          S=0.0
          SCALE=0.0
          IF (I.LE.M) THEN
             DO 11 K=I,M
               SCALE=SCALE+ABS(A(K,I))
11           CONTINUE
             IF (SCALE.NE.0.0) THEN
                DO 12 K=I,M
                  A(K,I)=A(K,I)/SCALE
                  S=S+A(K,I)*A(K,I)
12              CONTINUE
                F=A(I,I)
                G=-SIGN(SQRT(S),F)
                H=F*G-S
                A(I,I)=F-G
                IF (I.NE.N) THEN
                   DO 15 J=L,N
                     S=0.0
                     DO 13 K=I,M
                       S=S+A(K,I)*A(K,J)
13                   CONTINUE
                     F=S/H
                     DO 14 K=I,M
                       A(K,J)=A(K,J)+F*A(K,I)
14                   CONTINUE
15              CONTINUE
              ENDIF
              DO 16 K=I,M
                 A(K,I)=SCALE*A(K,I)
16            CONTINUE
           ENDIF
       ENDIF
       W(I)=SCALE*G
       G=0.0
       S=0.0
       SCALE=0.0
       IF ((I.LE.M).AND.(I.NE.N)) THEN
           DO 17 K=L,N
               SCALE=SCALE+ABS(A(I,K))
17           CONTINUE
             IF (SCALE.NE.0.0) THEN
                DO 18 K=L,N
                  A(I,K)=A(I,K)/SCALE
                  S=S+A(I,K)*A(I,K)
18              CONTINUE
                F=A(I,L)
                G=-SIGN(SQRT(S),F)
                H=F*G-S
                A(I,L)=F-G
                DO 19 K=L,N
                  RV1(K)=A(I,K)/H
19              CONTINUE
                IF (I.NE.M) THEN
                   DO 23 J=L,M
                     S=0.0
                     DO 21 K=L,N
                       S=S+A(J,K)*A(I,K)
21                   CONTINUE
                     DO 22 K=L,N
                       A(J,K)=A(J,K)+S*RV1(K)
22                   CONTINUE
23                 CONTINUE
              ENDIF
              DO 24 K=L,N
                 A(I,K)=SCALE*A(I,K)
24            CONTINUE
           ENDIF
       ENDIF
       ANORM=MAX(ANORM,(ABS(W(I))+ABS(RV1(I))))
25     CONTINUE
c        print *,'25'
C   Accumulation of right-hand transformations.
       DO 32 I=N,1,-1
       IF (I.LT.N) THEN
         IF (G.NE.0.0) THEN
           DO 26 J=L,N
             V(J,I)=(A(I,J)/A(I,L))/G
C   Double division to avoid posible underflow.
26          CONTINUE
          DO 29 J=L,N
            S=0.0
            DO 27 K=L,N
              S=S+A(I,K)*V(K,J)
27          CONTINUE
            DO 28 K=L,N
              V(K,J)=V(K,J)+S*V(K,I)
28          CONTINUE
29        CONTINUE
        ENDIF
        DO 31 J=L,N
          V(I,J)=0.0
          V(J,I)=0.0
31      CONTINUE
       ENDIF
       V(I,I)=1.0
       G=RV1(I)
       L=I
32     CONTINUE
c        print *,'32'

C  Accumulation of the left-hang transformation
       DO 39 I=N,1,-1
         L=I+1
         G=W(I)
         IF (I.LT.N) THEN
           DO 33 J=L,N
             A(I,J)=0.0
33         CONTINUE
         ENDIF
         IF (G.NE.0.0) THEN
           G=1.0/G
           IF (I.NE.N) THEN
             DO 36 J=L,N
               S=0.0
               DO 34 K=L,M
                 S=S+A(K,I)*A(K,J)
34             CONTINUE
               F=(S/A(I,I))*G
             DO 35 K=I,M
               A(K,J)=A(K,J)+F*A(K,I)
35           CONTINUE
36         CONTINUE
         ENDIF
        DO 37 J=I,M
          A(J,I)=A(J,I)*G
37      CONTINUE
       ELSE
         DO 38 J=I,M
           A(J,I)=0.0
38       CONTINUE
       ENDIF
       A(I,I)=A(I,I)+1.0
39     CONTINUE
c        print *,'39'

C   Diagonalization of the bidiagonal form
C   Loop over singular values
       DO 49 K=N,1,-1
C   Loop allowed iterations
         DO 48 ITS=1,30
C   Test for spliting
            DO 41 L=K,1,-1
              NM=L-1
C   Note that RV1(1) is always zero
              IF((ABS(RV1(L))+ANORM).EQ.ANORM) GO TO 2
              IF((ABS(W(NM))+ANORM).EQ.ANORM) GO TO 1
41          CONTINUE
c          print *,'41'
1         C=0.0
          S=1.0
          DO 43 I=L,K
            F=S*RV1(I)
            IF ((ABS(F)+ANORM).NE.ANORM) THEN
              G=W(I)
              H=SQRT(F*F+G*G)
              W(I)=H
              H=1.0/H
              C= (G*H)
              S=-(F*H)
              DO 42 J=1,M
                Y=A(J,NM)
                Z=A(J,I)
                A(J,NM)=(Y*C)+(Z*S)
                A(J,I)=-(Y*S)+(Z*C)
42            CONTINUE
            ENDIF
43        CONTINUE
c          print *,'43'
2         Z=W(K)
          IF (L.EQ.K) THEN
C   Convergence
            IF (Z.LT.0.0) THEN
C   Singular values are made nonnegative
              W(K)=-Z
              DO 44 J=1,N
                V(J,K)=-V(J,K)
44            CONTINUE
            ENDIF
            GO TO 3
          ENDIF
          IF (ITS.EQ.30) PAUSE 'No convergence in 30 iterations'
          X=W(L)
          NM=K-1
          Y=W(NM)
          G=RV1(NM)
          H=RV1(K)
          F=((Y-Z)*(Y+Z)+(G-H)*(G+H))/(2.0*H*Y)
          G=SQRT(F*F+1.0)
          F=((X-Z)*(X+Z)+H*((Y/(F+SIGN(G,F)))-H))/X
C   Next  QR  transformation
          C=1.0
          S=1.0
          DO 47 J=L,NM
            I=J+1
            G=RV1(I)
            Y=W(I)
            H=S*G
            G=C*G
            Z=SQRT(F*F+H*H)
            RV1(J)=Z
            C=F/Z
            S=H/Z
            F= (X*C)+(G*S)
            G=-(X*S)+(G*C)
            H=Y*S
            Y=Y*C
            DO 45 NM=1,N
              X=V(NM,J)
              Z=V(NM,I)
              V(NM,J)= (X*C)+(Z*S)
              V(NM,I)=-(X*S)+(Z*C)
45          CONTINUE
c            print *,'45',F,H
            Z=pythag(F,H)
            W(J)=Z
C   Rotation can be arbitrary if  Z=0.
            IF (Z.NE.0.0) THEN
c            print *,1/Z
              Z=1.0/Z
c              print *,'*'
              C=F*Z
              S=H*Z
            ENDIF
            F= (C*G)+(S*Y)
            X=-(S*G)+(C*Y)
            DO 46 NM=1,M
              Y=A(NM,J)
              Z=A(NM,I)
              A(NM,J)= (Y*C)+(Z*S)
              A(NM,I)=-(Y*S)+(Z*C)
46          CONTINUE
c          print *,'46'

47        CONTINUE
c          print *,'47'
          RV1(L)=0.0
          RV1(K)=F
          W(K)=X
48      CONTINUE
3      CONTINUE
49     CONTINUE
c        print *,'49'
       
       RETURN
       END
      
       FUNCTION pythag(a,b)
       REAL a,b,pythag
       REAL absa,absb
       absa=abs(a)
       absb=abs(b)
       IF (absa.GT.absb) THEN
          pythag=absa*SQRT(1.+(absb/absa)**2)
       ELSE
          IF (absb.EQ.0) THEN
             pythag=0.
          ELSE
             pythag=absb*SQRT(1.+(absa/absb)**2)
          ENDIF
       ENDIF
       RETURN
       END
        
      BLOCK DATA inithermite

      REAL*8  Z(126),H(126)
      DIMENSION NN(25)
      
      COMMON /QUADR/  Z,H,NN,NNW
      COMMON /EXPACC/ PMAX
      COMMON /RINT/   C,FC
      
      DATA NNW /13/
      DATA (NN(I),I=1,13)/2,3,4,5,6,7,8,9,10,12,16,20,24/
      DATA PMAX/20./
      DATA C/4.5/
      DATA (H(I),I=1,61)/1.,1.,0.555555555555556,0.888888888888889,
     * 0.555555555555556,0.347854845137454,0.652145154862546,
     * 0.652145154862546,0.347854845137454,0.236926885056189,
     * 0.478628670499366,0.568888888888889,0.478628670499366,
     * 0.236926885056189,0.171324492379170,0.360761573048139,
     * 0.467913934572691,0.467913934572691,0.360761573048139,
     * 0.171324492379170,0.129484966168870,0.279705391489277,
     * 0.381830050505119,0.417959183673469,0.381830050505119,
     * 0.279705391489277,0.129484966168870,0.101228536290376,
     * 0.222381034453374,0.313706645877887,0.362683783378362,
     * 0.362683783378362,0.313706645877887,0.222381034453374,
     * 0.101228536290376,0.081274388361574,0.180648160694857,
     * 0.260610696402935,0.312347077040003,0.330239355001260,
     * 0.312347077040003,0.260610696402935,0.180648160694857,
     * 0.081274388361574,0.066671344308688,0.149451349150581,
     * 0.219086362515982,0.269266719309996,0.295524224714753,
     * 0.295524224714753,0.269266719309996,0.219086362515982,
     * 0.149451349150581,0.066671344308688,0.047175336386512,
     * 0.106939325995318,0.160078328543346,0.203167426723066,
     * 0.233492536538355,0.249147048513403,0.249147048513403/
      DATA (H(I),I=62,101)/0.233492536538355,0.203167426723066,
     * 0.160078328543346,0.106939325995318,
     * 0.047175336386512,0.027152459411754094852,
     * 0.062253523938647892863,0.095158511682492784810,
     * 0.124628971255533872052,0.149595988816576732081,
     * 0.169156519395002538189,0.182603415044923588867, 
     * 0.189450610455068496285,0.189450610455068496285,
     * 0.182603415044923588867,0.169156519395002538189,
     * 0.149595988816576732081,0.124628971255533872052,
     * 0.095158511682492784810,0.062253523938647892863,
     * 0.027152459411754094852,0.017614007139152118312,
     * 0.040601429800386941331,0.062672048334109063570,
     * 0.083276741576704748725,0.101930119817240435037,
     * 0.118194531961518417312,0.131688638449176626898,
     * 0.142096109318382051329,0.149172986472603746788,
     * 0.152753387130725850698,0.152753387130725850698,
     * 0.149172986472603746788,0.142096109318382051329,
     * 0.131688638449176626898,0.118194531961518417312,
     * 0.101930119817240435037,0.083276741576704748725,
     * 0.062672048334109063570,0.040601429800386941331/
      DATA (H(I),I=102,126)/0.017614007139152118312,
     * 0.012341229799987199547, 0.028531388628933663181,
     * 0.044277438817419806169, 0.059298584915436780746,
     * 0.073346481411080305734, 0.086190161531953275917,
     * 0.097618652104113888270, 0.107444270115965634783,
     * 0.115505668053725601353, 0.121670472927803391204,
     * 0.125837456346828296121, 0.127938195346752156974,
     * 0.127938195346752156974, 0.125837456346828296121,
     * 0.121670472927803391204, 0.115505668053725601353,
     * 0.107444270115965634783, 0.097618652104113888270,
     * 0.086190161531953275917, 0.073346481411080305734,
     * 0.059298584915436780746, 0.044277438817419806169,
     * 0.028531388628933663181, 0.012341229799987199547/ 
        
      DATA (Z(I),I=1,58)/-0.577350269189626,0.577350269189626,
     *  -0.774596669241483,
     *  0., 0.774596669241483, -0.861136311594053, -0.339981043584856,
     *      0.339981043584856,  0.861136311594053, -0.906179845938664,
     *     -0.538469310105683,
     *  0., 0.538469310105683,  0.906179845938664, -0.932469514203152,
     *     -0.661209386466265, -0.238619186083197,  0.238619186083197,
     *      0.661209386466265,  0.932469514203152, -0.949107912342759,
     *     -0.741531185599394, -0.405845151377397,
     *  0., 0.405845151377397,  0.741531185599394,  0.949107912342759,
     *     -0.960289856497536, -0.796666477413627, -0.525532409916329,
     *     -0.183434642495650,  0.183434642495650,  0.525532409916329,
     *      0.796666477413627,  0.960289856497536, -0.968160239507626,
     *     -0.836031107326636, -0.613371432700590, -0.324253423403809,
     *  0., 0.324253423403809,  0.613371432700590,  0.836031107326636,
     *      0.968160239507626, -0.973906528517172, -0.865063366688985,
     *     -0.679409568299024, -0.433395394129247, -0.148874338981631,
     *      0.148874338981631,  0.433395394129247,  0.679409568299024,
     *      0.865063366688985,  0.973906528517172, -0.981560634246719,
     *     -0.904117256370475, -0.769902674194305, -0.587317954286617/
      DATA (Z(I),I=59,99)/-0.367831498198180, -0.125233408511469,
     *      0.125233408511469, 0.367831498198180,  0.587317954286617,
     *      0.769902674194305, 0.904117256370475,  0.981560634246719,
     *     -0.989400934991649932596,
     *     -0.944575023073232576078, -0.865631202387831743880,
     *     -0.755404408355003033895, -0.617876244402643748447,
     *     -0.458016777657227386342, -0.281603550779258913230,
     *     -0.095012509837637440185,  0.095012509837637440185,
     *      0.281603550779258913230,  0.458016777657227386342,
     *      0.617876244402643748447,  0.755404408355003033895,
     *      0.865631202387831743880,  0.944575023073232576078,
     *      0.989400934991649932596, -0.993128599185094924786,
     *     -0.963971927277913791268, -0.912234428251325905868,
     *     -0.839116971822218823395, -0.746331906460150792614,
     *     -0.636053680726515025453, -0.510867001950827098004,
     *     -0.373706088715419560673, -0.227785851141645078080,
     *     -0.076526521133497333755,  0.076526521133497333755,
     *      0.227785851141645078080,  0.373706088715419560673,
     *      0.510867001950827098004,  0.636053680726515025453,
     *      0.746331906460150792614,  0.839116971822218823395/
      DATA (Z(I),I=100,126)/0.912234428251325905868,
     *      0.963971927277913791268,  0.993128599185094924786,
     *     -0.995187219997021360180, -0.974728555971309498198,
     *     -0.938274552002732758524, -0.886415527004401034213,
     *     -0.820001985973902921954, -0.740124191578554364244,
     *     -0.648093651936975569252, -0.545421471388839535658,
     *     -0.433793507626045138487, -0.315042679696163374387,
     *     -0.191118867473616309159, -0.064056892862605626085,
     *      0.064056892862605626085,  0.191118867473616309159,
     *      0.315042679696163374387,  0.433793507626045138487,
     *      0.545421471388839535658,  0.648093651936975569252,
     *      0.740124191578554364244,  0.820001985973902921954,
     *      0.886415527004401034213,  0.938274552002732758524,
     *      0.974728555971309498198,  0.995187219997021360180/

      END
