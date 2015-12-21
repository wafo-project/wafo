C     Version  1991-XII-14

C  The MREG programs.
C
C
C  We consider a process X(I)=X(T(I)) at the grid of  N  points  T(1),...,T(N),
C
C          X(I) = -A(I) + Z*A(I+N) + Sum Xj*A(I+(j+1)*N) + Delta(I),
C
C  the sum disepirs if M=1, j=1,...,M-1. We assume that Z,Xj are independend
C  standart Rayleigh, Gaussian distributed rv. and independent of the zero
C  mean Gaussian residual process,  with covariance structure given in  R,
C
C          R(i+(j-1)N) = Cov (Delta(T(i)), Delta(T(j))).
C
C  Aditionally we have a zero mean Gaussian variable XN, 
C  indepensennt of Z,Xj with  covariance structure defined by 
C  B(i)= Cov (Delta(T(i)),XN), i=1,...,N,  B(N+1)=Var(XN). 
C  Furthermore XN and Z,Xj satisfies the following equation system
C
C     (BB + (XN,0,...,0)^T =  AA*(Z,X1,...,Xm-1)^T       (***)
C
C  where AA is (M,M) matrix, BB is M-vector. We rewrite this equation, by
C  introducing a variable X_M=XN/SQRT(XN) and construct   new matrix AA1
c  by adding the collon (SQRT(Var(XN)),0,...,0) and the row with only zeros.
C  The equations (***) writtes
C
C     (BB,0)^T  =  AA1*(Z,X1,...,Xm-1,Xm)^T       (****)
C
C  where AA1 is (M+1,M+1) matrix, We assume that the rank of AA1 is M,
C  otherwise the density is singular and we give a output F=0.CC
C
C  Let Y0 be a zero-mean Gaussian variable indepensennt of Z,Xj 
C  with  covariance structure defined by 
C  DB(i)= Cov (Delta(T(i)),Y0), i=1,...,N,  DB(N+1)=Cov(XN,Y0), Var(Y0)=VDER. 
C  Let Y be defined by
C
C     Y=-DA(1) + Z*DA(2) + Sum Xj*DA(2+j) +Y0,
C
C  j=1,...,M-1. The program computes:
C
C  F = E[ Y^+ *1{ HH<X(I)<0 for all I, I=1,...,N}|Z,X1,...,X_M-1 solves (***)]
C      *f_{Z,X1,....,XM-1}(***).
C
C  In the simplest case NIT=0 we define (Delta(1),...,Delta(N),XN)=0. 
C
C  We renormalize  vectors AN and DA, the covariance fkn R, DB
C  and VDER. Then by renormalization we choose the Gaussian variable X such
C  that F is written in the form
C
C  F = E[(D0(1)+X*D0(2)+Y1)^+*(PC+X*PD)^+*1{HH <A0(I)+X*B0(I)+Delta1(I)<0}]
C
C  Observe, PC+X*PD>0 defines integration region for X.
C  In the simplest case NIT=0 we define (Delta(1),...,Delta(N),Y1)=0. 
C  For NIT=1 only (Delta(1),...,Delta(N))=0, i.e. we have to compute
C  a one dimensional integral. Finally by conditioning on X the problem is
C  put in the format of RIND-problem.
C
C  MN is fysikal dimension for the matrix AA, INF indicates whether one
C  has already called the subroutine before and ONLY! inputs BB, DA or A
C  was changed.
C
C  Observe the limitations are :  N<=100, 0<M<=MN<=5.



      SUBROUTINE TWOG(F,R,B,DB,AA,BB,A,DA,VDER,EPS0,IAC,M,MN,N,NIT,INF)
      real*8 XX1,H1,SQ
      DIMENSION R(1),B(1),DB(1),A(1),DA(1)
      DIMENSION AA(MN,MN),BB(1),AA1(6,6),AO(6)
      DIMENSION W1(6),V1(6,6),DA0(6),A0(202),B0(202)
      DIMENSION R1(10201),SQ(202),B1(202),DB1(202),A1(1414),D0(2)
      DIMENSION XX1(24),H1(24)
      COMMON /EPS/EPS,EPSS,CEPSS
      COMMON/RINT/ EXB,FEXB
C
C  If INF=0 we have to initiate conditioning and renormalization transf.
C
      INF1=0
      F=0.
      IF (N.LT.0) RETURN
      IF (N.EQ.0) GO TO 2
      DO 1 I=1,N
      DO 1 I1=1,M+1
      A1(I+(I1-1)*N)=A(I+(I1-1)*N)
1     CONTINUE
2     CONTINUE
      DO 3 I=1,M+1
      DA0(I)=DA(I)
3     CONTINUE

C
C  Renormalization
C
      IF (INF.EQ.1) GO TO 105
      NNIT=MIN(NIT,N)

      DO 4 i=1,n
             sq(i)=0.
4     continue
      
      DO 5 I=1,M
      DO 5 J=1,M
         AA1(I,J)=AA(I,J)
5     CONTINUE

      NNIT=MIN(NIT,N)
      
      QD=B(N+1)
      DB1N=DB(N+1)
      VDER1=VDER

      IF(QD.le.eps) then

         db1n=0.
         sqd=0.
         nnit=0
         DO 7 I=1,N
            DB1(I)=DB(I)
            A1(I+(M+1)*N)=0.
            SQ(I)=0.
7        CONTINUE

         else

           SQD=SQRT(QD)            
           VDER1=VDER1-DB1N*DB1N/QD
           DB1N=DB1N/SQD
           DO 8 I=1,N
              DB1(I)=DB(I)-DB(N+1)*(B(I)/QD)
              A1(I+(M+1)*N)=B(I)/SQD
8         CONTINUE
      
      end if

      SDER1=0.
      IF (VDER1.GT.EPS) SDER1=SQRT(VDER1)
C      print *,'sqd,SDER1',sqd,SDER1     
      DA0(M+2)=DB1N
      BB(M+1)=0.

      
      AA1(1,M+1)=SQD
      DO 30 I=1,M
       AA1(I+1,M+1)=0.
       AA1(M+1,I)=0.
30    CONTINUE


      CALL SVDCMP(AA1,M+1,M+1,6,6,W1,V1)
      
      DET1=1.
      idet=0
      DO 35 I=1,M+1
       IF ( W1(I).LT.0.00001 ) THEN
          idet=idet+1
          W1(I)=0.
         DO 37 J=1,M+1
            AO(J)=V1(J,I)
37      CONTINUE
        GO TO 35
        END IF
        DET1=DET1*W1(I)
35    CONTINUE
C      print *,'det1',det1

      IF(DET1.LT.0.001) NNIT=0
c
c    Obs. QD can not be small since NNIT>0
c      
      IF (NNIT.GT.1) THEN
      DO 10 I=1,N
         XR1=R(I+(I-1)*N)-B(I)*(B(I)/QD) 
         IF(XR1.GT.EPS) SQ(I)=SQRT(XR1)
10    CONTINUE
           
      DO 15 I=1,N
      DO 15 J=1,N
         R1(J+(I-1)*N)=R(J+(I-1)*N)-B(I)*(B(J)/QD)
15    CONTINUE
      
      
      END IF


105   CONTINUE
      if (idet.gt.1) return
C
C  Renormalization is dome
C
      CALL R_ORT(CC,PC,PD,AA1,V1,W1,AO,BB,A1,A0,B0,DA0,D0,DET1,M+1,6,N)
      IF(CC.LT.0.) RETURN
      XMI=-EXB
      XMA=EXB
      IF(ABS(PD).LE.EPS.AND.PC.LT.0.) RETURN
      IF(ABS(PD).LE.EPS) GO TO 102
      X=-PC/PD
      IF(PD.GT.0..AND.XMI.LT.X) XMI=X
      IF(PD.LT.0..AND.XMA.GT.X) XMA=X
102   CONTINUE
c      PRINT *,'XMI,XMA',XMI,XMA
      IF(NNIT.eq.1.AND.IAC.LT.1.OR.NNIT.eq.0.OR.XMI.GE.XMA) THEN
      CALL C1_C2(XMI,XMA,A0,B0,D0(1),D0(2),0.,SQ,N)
c         PRINT *,'XMI,XMA',XMI,XMA
         F=GAUSINT(XMI,XMA,D0(1),D0(2),PC,PD)*CC
c         print *,'return',f,cc
         RETURN
      END IF
C
C ***********************************************************
C
C We shall condition on the values of X, XMI<X<XMA, but for some
C X values XIND will be zero leading to reduced accuracy. Hence we try
C to exclude them and narrow the interval [XMI,XMA]
c
c      PRINT *,XMI,XMA

      CALL C1_C2(XMI,XMA,A0,B0,D0(1),D0(2),SDER1,SQ,N)
c      PRINT *,XMI,XMA
      IF(FI(XMA)-FI(XMI).LT.EPSS) RETURN
      CALL GAUSS1(N1,H1,XX1,XMI,XMA,EPS0)
      DO 100 I2=1,N1
      FR1=CC*H1(I2)*(PC+XX1(I2)*PD)
      DO 40 I=1,N
      B1(I)=A0(I)+XX1(I2)*B0(I)
40    CONTINUE
      DB0N=D0(1)+XX1(I2)*D0(2)

C
C  INF1=1 means that both R1 and SQ are the same as in the previous
C  call of TWOREG subroutine, INF=0 indicates the new R and SQ.
c      
c      print *,'go in rind'
      CALL RIND(XIND,NNIT-1,R1,B1,DB0N,DB1,SQ,VDER1,EPS0,IAC,N,INF1)
c      if (n.gt.29) print *,XIND,DB0N,VDER1,b1(N),b1(n-1)
      INF1=1
 100  F=F+FR1*XIND
      RETURN
      END


      SUBROUTINE R_ORT(C,PC,PD,AA1,V1,W1,AO,BB,A0,A,B,DA,D0,DET,M,MN,N)
      DIMENSION A(1),B(1),DA(1),AA1(MN,MN),V1(MN,MN),W1(1),BB(1),A0(1)
      DIMENSION XO(6),AO(1),D0(1)
      COMMON/EXPACC/EACC
      DATA SP/0.398942/
      C=-999.
      CALL SVBKSB(AA1,W1,V1,M,M,MN,MN,BB,XO)
      P=0.
      DER0=-DA(1)
      DER1=0.
      DO 10 I=1,M
      P=P+XO(I)*XO(I)
      DER0=DER0+XO(I)*DA(I+1)
      DER1=DER1+AO(I)*DA(I+1)
10    CONTINUE
      IF (P.GT.2.*EACC) RETURN
      C=(SP**(M-2))*EXP(-0.5*P)/ABS(DET)
c      print *,'XO',XO(1),XO(2),XO(3),XO(4)
c      print *,'AO',AO(1),AO(2),AO(3),AO(4)
      if(n.lt.1) go to 100
      DO 20 I=1,N
      A(I)=-A0(I)
      B(I)=0.
      DO 20 J=1,M
      B(I)=B(I)+AO(J)*A0(I+J*N)
      A(I)=A(I)+XO(J)*A0(I+J*N)
20    CONTINUE
100   continue
      D0(1)=DER0
      D0(2)=DER1
      PC=XO(1)
      PD=AO(1)
      RETURN
      END



