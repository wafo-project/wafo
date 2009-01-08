C     Version  1996-I-15 

C  This is a FORTRAN version of mctp2rfc.m and iter.m programs.
C  The maximal dimension of Markov and rainflow matrices is 200.

      PROGRAM MCTP2RFC

      real*8 Frfc(40000),FmM(40000),FMax(200),FMin(200)
      real*8 AA(40000),FMI(200),FMA(200),FNN(200),E(200)
      real*8 BB(40000),FNNX(200),EX(200),PB(40000),PA(40000)
      real*8 XI(40000),FRFC1(40000),FMM1(40000)
      real*8 X,fx,Y,DRFC,FMA1,FMI1,FN1,ENA,epsilon
      real*8 check0,FMaMi,sa,sra,fnorm,check_i,check
      real*8 afx,ax,afnx,aex,aax,aafx,eps0,eps1,check1

      CALL READ_I(FRFC,FMM,epsilon,kkk,N,eps0,eps1)
c      print *,kkk,N
      KITER=0

       check=0.
       do 888 I=1,N*N
       check=check+FRFC(I)
888    continue
       check=check*epsilon

c       DO 3 I=1,N
c	  FMax(I)=0.
c	  FMin(I)=0.
c3      CONTINUE           

c       FMaMI=0.
      
c       DO 4 I=1,N
c	  DO 5 J=1,N
c	     FMax(I)=FMax(I)+FRFC(I+(J-1)*N)
c	     FMin(I)=FMin(I)+FRFC(J+(I-1)*N)
c5         continue
c          IF(FMax(I).GT.FMaMi) FMaMi=FMAX(I)
c          IF(FMin(I).GT.FMaMi) FMaMi=FMIN(I)
c4      CONTINUE            

1      CONTINUE      

       KITER=KITER+1

       DO 2 I=1,N*N
       FMM1(I)=FMM(I)
2      CONTINUE      

      DO 3 II=1,N
	  FMax(II)=0.
	  FMin(II)=0.
3      CONTINUE           

       FMaMI=0.
       DO 4 II=1,N
	  DO 5 JJ=1,N
	     FMax(II)=FMax(II)+FMM1(II+(JJ-1)*N)
	     FMin(II)=FMin(II)+FMM1(JJ+(II-1)*N)
5         continue
          IF(FMax(II).GT.FMaMi) FMaMi=FMAX(II)
          IF(FMin(II).GT.FMaMi) FMaMi=FMIN(II)
4      CONTINUE            


      check0=eps0*FMaMi
      check1=eps1



      DO 10 I=1,N
	 DO 10 J=1,N
	  FRFC1(J+(I-1)*N)=0.
10    continue

      DO 98 k=2,N
	 DO 99 i=2,K
	 NA=I-1
	 SA=0.
	 SRA=0.
	   DO 20 J1=1,NA
	      FMA(J1)=FMAX(N-K+J1)
	      FMI(J1)=FMIN(K-I+J1)
	      DO 20 J2=1,NA
		 AA(J1+(J2-1)*NA)=FMM1(N-K+J1+(K-I+J2-1)*N)
		 SA=SA+FMM1(N-K+J1+(K-I+J2-1)*N)
		 SRA=SRA+FRFC1(N-K+J1+(K-I+J2-1)*N)
20         CONTINUE
        x=aa(1)
        SA=SA-X
C
C Our equation is Y+DRFC=X+F(X), here X is known.
C
        AA(1)=0. 
	DRFC=SRA-SA

        IF(FMA(1).lt.check0.or.FMI(1).lt.check0) go to 99
	
	DO 40 J1=1,NA
	   FNN(J1)=FMA(J1)
	   E(J1)=FMI(NA-J1+1)
40      CONTINUE

	DO 50 J1=1,NA
	     DO 50 J2=1,NA
		BB(J1+(J2-1)*NA)=AA(J2+(NA-J1)*NA)
		FNN(J1)=FNN(J1)-AA(J1+(J2-1)*NA)
		E(J1)=E(J1)-AA(J2+(NA-J1)*NA)
50      CONTINUE
    
       AFX=0.
       AX=0.
       DO 55 J1=2,NA 
	  IF(FNN(J1).GT.check1*FMA(J1)) AFX=1.
	  IF(E(J1-1).GT.check1*FMI(J1-1)) AX=1.
55     continue

       FX=0. 
       AFNX=0.
       AEX=0.
       if (FNN(1)-x.gt.check1*FMA(1)) AFNX=1.
       if (E(NA)-x.gt.check1*FMI(NA)) AEX=1.
c       if (min(AFNX,AEX).lt.0.5) go to 999

c        if (k.eq.N-20) then
c        print *,i,afnx,aex,x,fnn(1)-x,e(na)-x
c        print *,'*',FMI(1),FMA(1)
c        pause
c        end if

       AAX=max(AEX,AX)
       AAFX=MAX(AFNX,AFX)

       if (min(AAX,AAFX).GT.check0) then

	   DO 60 J1=1,NA
	       Fnorm=FMA(j1)
	       if (Fnorm.LT.check0) FNORM=1. 
		  DO 65 J2=1,NA
		     PA(J1+(NA-J2)*NA)=AA(J1+(J2-1)*NA)/FNORM
65                CONTINUE
60      CONTINUE
	   
	   DO 70 J1=1,NA
	       Fnorm=FMI(NA-J1+1)
	       if (Fnorm.lt.CHECK0) FNORM=1.
		  DO 75 J2=1,NA
		    PB(J1+(J2-1)*NA)=BB(J1+(J2-1)*NA)/FNORM
75                 CONTINUE
		    E(j1)=E(j1)/Fnorm
70      CONTINUE
       DO 80 J1=1,NA      
	  FNNx(J1)=FNN(J1)
	  Ex(J1)=E(J1)
80     CONTINUE

       DO 110 I0=1,NA*NA
	  XI(I0)=0.
110     CONTINUE
       IF (NA.GT.1) then
         DO 120 I1=1,NA-1
	  DO 130 J=1,NA-1
	     DO 135 jj=1,NA
      XI(I1+(J-1)*NA)=XI(I1+(J-1)*NA)-PB(I1+(jj-1)*NA)*PA(jj+(J-1)*NA)
135           continue
130        CONTINUE
	  XI(I1+(I1-1)*NA)=1.+XI(I1+(I1-1)*NA)
120     CONTINUE
        end if
        XI(NA*NA)=1.
        FMA1=FMA(1)
        FMI1=FMI(1)
        FN1=FNN(1)
        ENA=E(NA)

      
       CALL M2RFC(FX,x,XI,pA,pB,FNNX,EX,FMA1,FMI1,FN1,eNA,NA)
       end if

999    continue
       Y=max(0.,X+FX-DRFC)
       Frfc1(N-k+1+(k-i)*N)=Y
c        if (k.eq.N-20) then
c        print *,k,i,X
c        pause
c        write(10,600) x,y,fx,drfc
c        end if
c        if (k.gt.N-20) stop
 99   CONTINUE       
      PRINT *,'Number of loops left',max(0,(Kkk-KITER))*N+(N-k+1)
 98   CONTINUE
      
      Check_i=0.

      if (kkk.gt.0) then      

       DO 101 i=1,N*N
         fMM(I)=fMM1(I)+(Frfc(I)-Frfc1(I))
         fMM(I)=max(0.,fMM(I))
         check_i=check_i+abs(fmm(i)-fmm1(i))
101    continue
      end if

      if (check_i.gt.check.and.kiter.lt.kkk) go to 1 
  
       
    
      DO 100 I=1,N
	 DO 100 J=1,N
	   write(3,500) Frfc1(J+(I-1)*N)
	   write(4,500) Fmm(J+(I-1)*N)
100    continue
500    format(2x,e14.8)
600    format(4(2x,e14.8),2x,i4)
      CLOSE(UNIT=3)
      CLOSE(UNIT=4)
      STOP
      END


      SUBROUTINE READ_I(FRFC,FMM,epsilon,k,N,eps0,eps1)
      real*8 FMM(1),FRFC(1)
      real*8 epsilon,eps0,eps1 
      INTEGER N,k

      OPEN(UNIT=1,FILE='rfc_mat.in')
      OPEN(UNIT=2,FILE='mark_mat.in')
      OPEN(UNIT=3,FILE='rfc_mat.out')
      OPEN(UNIT=4,FILE='mark_mat.out')
      OPEN(UNIT=9,FILE='epsilon.in')
      OPEN(UNIT=10,FILE='slask.out')

      READ(9,*) epsilon,eps0,eps1
      read(9,*) k
      CLOSE(UNIT=9)  
      NG=1
 12   READ (1,*,END=11) FRFC(NG)
      NG=NG+1
      GO TO 12
 11   CONTINUE
      NG=NG-1
      IF (NG.GT.40000) THEN
      PRINT *,'rainlow matrix has dimension exciding 200, stop'
      STOP
      END IF
      IF (NG.LT.4) THEN
      PRINT *,'rainlow matrix has dimension smaller than 2, stop'
      STOP
      END IF
      CLOSE(UNIT=1)
      
      NG1=1
 22   READ (2,*,END=21) FMM(NG1)
      NG1=NG1+1
      GO TO 22
 21   CONTINUE
      NG1=NG1-1
      IF (NG1.GT.40000) THEN
      PRINT *,'Markov matrix has dimension exciding 200, stop'
      STOP
      END IF
      IF (NG1.LT.4) THEN
      PRINT *,'Markov matrix has dimension smaller than 2, stop'
      STOP
      END IF
c      print *,epsilon,K    
      CLOSE(UNIT=2)

      IF (NG.NE.NG1) THEN
      PRINT *,'Markov and rainflow matrices have diff. dim.'
      STOP
      END IF

      N=0
20    N=N+1
      IF (N*N.LT.NG) go to 20
      RETURN
      END



       SUBROUTINE M2RFC(FX,X,YI,AA,BB,FNX,EX,FMA1,FMI1,F1,eNA,NA)
       real*8 YI(1),AA(1),BB(1),FNX(1),EX(1)
       real*8 XI(40000),V(40000),W(200),FY(200)
       real*8 FX,X,FMA1,FMI1,F1,ena
       AA(1+(nA-1)*nA)=x/FMA1
       BB(nA)=x/FMI1    
       FNX(1)=F1-x
       ex(nA)=enA-x/FMI1

       FX=0.
       if (nA.EQ.1) THEN
	   fx=FNx(1)*(AA(1)/(YI(1)-BB(1)*AA(1))*ex(1))
	   RETURN
       ENDIF

       DO 10 I=1,NA*NA
	  XI(I)=YI(I)
10     CONTINUE
	  DO 20 J=1,NA
	     DO 25 jj=1,NA
      XI(NA+(J-1)*NA)=XI(NA+(J-1)*NA)-BB(NA+(jj-1)*NA)*AA(jj+(J-1)*NA)
25           continue
20        CONTINUE
	  DO 30 I=1,NA-1
	     DO 35 jj=1,NA
      XI(I+(NA-1)*NA)=XI(I+(NA-1)*NA)-BB(I+(jj-1)*NA)*AA(jj+(NA-1)*NA)
35           continue
30        CONTINUE

       CALL SVDCMP(XI,NA,NA,NA,NA,W,V)
       CALL SVBKSB(XI,W,V,NA,NA,NA,NA,EX,FY)
       
       DO 40 I=1,NA
	  DO 40 J=1,NA
	     FX=FX+FNX(I)*AA(I+(J-1)*NA)*FY(J)
40     CONTINUE
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
       PARAMETER (NMAX=200)
C   Maximum anticipated value of N
       real*8 U(MP,NP),W(NP),V(NP,NP),B(MP),X(NP),TMP(NMAX)
       real*8 s
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
C  output. The diagonal matrix of singular values  V  is output as 
C  a vector W. 
C  The matrix  V (not the transpose  V^T) is output as  V.  M  must be
C  greater or equal to  N; if it is smaller, then  A  should be filled up
C  to square with zero rows.
C
       PARAMETER (NMAX=200)
C  Maximum anticipated values of  N
       real*8 A(MP,NP),W(NP),V(NP,NP),RV1(NMAX)
       real*8 g,scale,anorm,s,f,h,c,y,z,x
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
	  IF (ITS.EQ.300) PAUSE 'No convergence in 300 iterations'
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
       REAL*8 a,b,pythag
       REAL*8 absa,absb
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
