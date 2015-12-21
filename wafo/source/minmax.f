C     Version  1994-X-18 

C     This is a new version of WAMP program computing crest-trough wavelength 
C     and amplitude  density.

      PROGRAM Maxmin
      real*8 SQ0,SQ1,SQ
      DIMENSION UVdens(101,101),HHT(101)
      DIMENSION T(101),Ulev(101),Vlev(101)
      DIMENSION VT(101),UT(101),Vdd(101),Udd(101)
      DIMENSION COV(505),R(10201),R1(10201),R2(10201),R3(10201)
      DIMENSION AA(3,3),BB(4),DAI(5),AI(707)

C
C   The program computes the joint density of maximum the following minimum
C   and the distance between Max and min for a zero-mean stationary 
C   Gaussian process with covariance function defined explicitely with 4 
C   derivatives. The process should be normalized so that the first and 
C   the second spectral moments are equal to 1. The values of Max are taken
C   as the nodes at Hermite-Quadrature and then integrated out so that
C   the output is a joint density of wavelength T and amplitude H=Max-min.
C   The Max values are defined by subroutine Gauss_M with the accuracy
C   input  epsu. The principle is that the integral of the marginal density 
C   of f_Max is computed with sufficient accuracy.
C
      DIMENSION B0(101),DB0(101),DDB0(101)
      DIMENSION B1(101),DB1(101),DDB1(101)
      DIMENSION DB2(101),DDB2(101)
      DIMENSION Q(101),SQ(101),VDER(101),DBI(101),BI(101)
      COMMON/CHECK1/III01,III11,III21,III31,III41,III51
     *,III61,III71,III81,III91,III101 
      COMMON/CHECKQ/III
      COMMON /EPS/  EPS,EPSS,CEPSS

C
C  Initiation of all constants and integration nodes 'INITINTEG'
C
      CALL INITINTEG(EPS0,IAC,NIT)
c
c   OBS. we are using the variables AI,R,R1,R2 as a temporary storage 
C   for transformation  g  of the process.



c
      CALL INITLEVELS(Ulev,NU,Vlev,NV,T,HHT,N,R1,R2,NG)
      IF( R1(1) .gt. R1(ng))  then
      do 13 i=1,ng
      AI(i)=R1(i)
      R(I)=R2(i)
13    continue
      do 17 i=1,ng
      R1(i)=AI(ng-i+1)
      R2(i)=R(ng-i+1)
17    continue
      end if
      if(abs(R1(ng)-R1(1))*abs(R2(ng)-R2(1)).lt.0.01) then
      print *,'The transformation  g is singular, stop'
      stop
      end if
         DO 14 IV=1,Nv
                V=Vlev(IV)
                CALL TRANSF(NG,V,R2,R1,VALUE,DER) 
                VT(IV)=VALUE
                Vdd(Iv)=DER
14       continue
          DO 16 IU=1,Nu
                U=Ulev(IU)
                CALL TRANSF(NG,U,R2,R1,VALUE,DER) 
                UT(IU)=VALUE
                Udd(IU)=DER
                do 16 iv=1,Nv
                UVdens(IU,IV)=0.
16        CONTINUE
      
      
      CALL COVG(XL0,XL2,XL4,COV,R1,R2,R3,T,N)
      
      
      Q0=XL4
      IF (Q0.le.1.+eps) then
      Print *,'Covariance structure is singular, stop.'
      stop
      end if
      SQ0=SQRT(Q0)
      Q1=XL0-XL2*XL2/XL4
      IF (Q1.le.eps) then
      Print *,'Covariance structure is singular, stop.'
      stop
      end if
      SQ1=SQRT(Q1)
      DO 10 I=1,N
      B0(I)=-COV(I+2*N)
      DB0(I)=-COV(I+3*N)
      DDB0(I)=-COV(I+4*N)
      
      B1(I)=COV(I)+COV(I+2*N)*(XL2/XL4)
      DB1(I)=COV(I+N)+COV(I+3*N)*(XL2/XL4)
      DDB1(I)=COV(I+2*N)+XL2*(COV(I+4*N)/XL4)
C
C       Q(I) contains Var(X(T(i))|X'(0),X''(0),X(0))
C    VDER(I) contains Var(X''(T(i))|X'(0),X''(0),X(0))
C
      Q(I)=XL0 - COV(I+N)*(COV(I+N)/XL2) - B0(I)*(B0(I)/Q0)
     1     -B1(I)*(B1(I)/Q1)
      VDER(I)=XL4 - (COV(I+3*N)*COV(I+3*N))/XL2 - (DDB0(I)*DDB0(I))/Q0
     1     - (DDB1(I)*DDB1(I))/Q1 
      

C
C       DDB2(I) contains Cov(X''(T(i)),X(T(i))|X'(0),X''(0),X(0))
C
      DDB2(I)=-XL2 - (COV(I+N)*COV(I+3*N))/XL2 - DDB0(I)*(B0(I)/Q0)
     1     -DDB1(I)*(B1(I)/Q1)
      IF(Q(I).LE.eps) then
      sq(i)=0.
      ddb2(i)=0.
      else
      SQ(I)=SQRT(Q(I))
C
C       VDER(I) contains Var(X''(T(i))|X'(0),X''(0),X(0),X(T(i))
C
      
      VDER(I)=VDER(I) - (DDB2(I)*DDB2(I))/Q(I)
      end if
      
10    CONTINUE
	
      DO 15 I=1,N
      DO 15 J=1,N
C      
C   R1 contains Cov(X(T(I)),X'(T(J))|X'(0),X''(0),X(0)) 
C
      R1(J+(I-1)*N)=R1(J+(I-1)*N) -  COV(I+N)*(COV(J+2*N)/XL2) 
     1 -  (B0(I)*DB0(J)/Q0) -  (B1(I)*DB1(J)/Q1) 
      
C      
C   R2 contains Cov(X'(T(I)),X'(T(J))|X'(0),X''(0),X(0))
C
      R2(J+(I-1)*N) = -R2(J+(I-1)*N) - COV(I+2*N)*(COV(J+2*N)/XL2) 
     1   - DB0(I)*DB0(J)/Q0  - DB1(I)*(DB1(J)/Q1) 
C      
C   R3 contains Cov(X''(T(I)),X'(T(J))|X'(0),X''(0),X(0))
C
      R3(J+(I-1)*N) = R3(J+(I-1)*N) - COV(I+3*N)*(COV(J+2*N)/XL2) 
     1   - DB0(J)*(DDB0(I)/Q0)  - DDB1(I)*(DB1(J)/Q1) 
15    CONTINUE

C  The initiations are finished and we are beginnings with 3 loops
C  on T=T(I), U=Ulevels(IU), V=Ulevels(IV), U>V.
      

      DO 20 I=1,N

           NNIT=NIT
           IF (Q(I).LE.EPS) GO TO 20

         DO 30 I1=1,I
            DB2(I1)=R1(I1+(I-1)*N) 

C     Cov(X'(T(I1)),X(T(i))|X'(0),X''(0),X(0)) 
C     DDB2(I) contains Cov(X''(T(i)),X(T(i))|X'(0),X''(0),X(0))

 30      CONTINUE
         DO 50 I3=1,I
            DBI(I3) = R3(I3+(I-1)*N) - (DDB2(I)*DB2(I3)/Q(I)) 
            BI(I3)  = R2(I3+(I-1)*N) - (DB2(I)*DB2(I3)/Q(I))
 50      CONTINUE
         DO 51 I3=1,I-1
            AI(I3)=0.
            AI(I3+I-1)=DB0(I3)/SQ0 
            AI(I3+2*(I-1))=DB1(I3)/SQ1
            AI(I3+3*(I-1))=DB2(I3)/SQ(I)
 51      CONTINUE
         VDERI=VDER(I)
         DAI(1)=0.
         DAI(2)=DDB0(I)/SQ0
         DAI(3)=DDB1(I)/SQ1 
         DAI(4)=DDB2(I)/SQ(I) 
         AA(1,1)=DB0(I)/SQ0 
         AA(1,2)=DB1(I)/SQ1 
         AA(1,3)=DB2(I)/SQ(I) 
         AA(2,1)=XL2/SQ0
         AA(2,2)=SQ1
         AA(2,3)=0.
         AA(3,1)=B0(I)/SQ0 
         AA(3,2)=B1(I)/SQ1 
         AA(3,3)=SQ(I)
         IF (BI(I).LE.EPS) NNIT=0
         IF (NNIT.GT.1) THEN
            IF(I.LT.1) GO TO 41
            DO 40 I1=1,I-1
               DO 40 I2=1,I-1

C   R contains Cov(X'(T(I1)),X'(T(I2))|X'(0),X''(0),X(0),X(I)) 

      R(I2+(I1-1)*(I-1))=R2(I2+(I1-1)*N)-(DB2(I1)*DB2(I2)/Q(I)) 

 40         CONTINUE
 41         CONTINUE
         END IF

C  Here the covariance of the problem would be innitiated

            INF=0
            Print *,'   Laps to go:',N-I+1
         DO 80 IV=1,Nv
            V=VT(IV)
            IF (ABS(V).GT.5.) GO TO 80                  
            IF (Vdd(IV).LT.EPS0) GO TO 80
            DO 60 IU=1,Nu
                  U=UT(IU)
                  IF (ABS(U).GT.5.) GO TO 60                  
                  IF (Udd(IU).LT.EPS0) GO TO 60
                  BB(1)=0.
                  BB(2)=U
                  BB(3)=V
                  if(u.le.v) go to 60
		
! 		if (IV.EQ.33.AND.IU.EQ.6.AND.I.EQ.11) THEN
!			fffff = 10
!		endif
      CALL TWOG(F,R,BI,DBI,AA,BB,AI,DAI,VDERI,EPS0,IAC,3,3,I-1,NNIT,INF)
                  INF=1
            UVdens(IU,IV) = UVdens(IU,IV) + Udd(IU)*Vdd(IV)*HHT(I)*F
	if (F>0) THEN
!	if (N-I+1 .eq. 38.and.IV.EQ.26.AND.IU.EQ.16) THEN	
!	if (IV.EQ.26.AND.IU.EQ.16) THEN		
          PRINT * ,' R:', R(1:I)
	    PRINT * ,' BI:', BI(1:I)
		PRINT * ,' DBI:', DBI(1:I)
		PRINT * ,' DB2:', DB2(1:I)
		PRINT * ,' DB0(1):', DB0(1)
		PRINT * ,' DB1(1):', DB1(1)
          PRINT * ,' DAI:', DAI
		PRINT * ,' BB:', BB
	    PRINT * ,' VDERI:', VDERI
          PRINT * ,' F    :', F
		PRINT * ,' UVDENS :',  UVdens(IU,IV)
	    fffff = 10
	endif
 60         CONTINUE
 80       continue
 20   CONTINUE
       hhhh=0.
       do 90 Iu=1,Nu
       do 90 Iv=1,Nv
       WRITE(10,300) Ulev(iu),Vlev(iv),UVdens(iu,iv)
       hhhh=hhhh+UVdens(iu,iv)
 90    continue
      if (nu.gt.1.and.nv.gt.1) then
      write(11,*) 'SumSum f_uv *du*dv='
     1,(Ulev(2)-Ulev(1))*(Vlev(2)-Vlev(1))*hhhh
      end if
      
      sder=sqrt(XL4-XL2*XL2/XL0)
      cder=-XL2/sqrt(XL0)
      const=1/sqrt(XL0*XL4)
      DO 95 IU=1,NU     
        U=UT(IU)
        FM=Udd(IU)*const*exp(-0.5*U*U/XL0)*PMEAN(-cder*U,sder)
        WRITE(9,300) Ulev(IU),FM
 95   continue      
      DO 105 IV=1,NV     
        V=VT(IV)
        VV=cder*V
        Fm=Vdd(IV)*const*exp(-0.5*V*V/XL0)*PMEAN(VV,sder)
        WRITE(8,300) Vlev(IV),Fm
 105   continue 
      if (iii.eq.0) iii=1

      write(11,*) 'Rate of calls RINDT0:',float(iii01)/float(iii)
      write(11,*) 'Rate of calls RINDT1:',float(iii11)/float(iii)
      write(11,*) 'Rate of calls RINDT2:',float(iii21)/float(iii)
      write(11,*) 'Rate of calls RINDT3:',float(iii31)/float(iii)
      write(11,*) 'Rate of calls RINDT4:',float(iii41)/float(iii)
      write(11,*) 'Rate of calls RINDT5:',float(iii51)/float(iii)
      write(11,*) 'Rate of calls RINDT6:',float(iii61)/float(iii)
      write(11,*) 'Rate of calls RINDT7:',float(iii71)/float(iii)
      write(11,*) 'Rate of calls RINDT8:',float(iii81)/float(iii)
      write(11,*) 'Rate of calls RINDT9:',float(iii91)/float(iii)
      write(11,*) 'Rate of calls RINDT10:',float(iii101)/float(iii)
      write(11,*) 'Number of calls of RINDT*',iii

      CLOSE(UNIT=8)
      CLOSE(UNIT=9)
      CLOSE(UNIT=10)
      CLOSE(UNIT=11)

 300  FORMAT(4(3X,F10.6))
      STOP
      END

      SUBROUTINE INITLEVELS(ULEVELS,NU,Vlevels,Nv,T,HT,N,TG,XG,NG)
      DIMENSION ULEVELS(1),Vlevels(1),T(1),HT(1),TG(1),XG(1),HH(101)
      COMMON/TBR/HH
      OPEN(UNIT=2,FILE='transf.in')
      OPEN(UNIT=4,FILE='Mm.in')
      OPEN(UNIT=3,FILE='t.in')
      

      NG=1
 12   READ (2,*,END=11) TG(NG),XG(NG)
      NG=NG+1
      GO TO 12
 11   CONTINUE
      NG=NG-1
      IF (NG.GT.501) THEN
      PRINT *,'Vector defining transformation of data > 501, stop'
      STOP
      END IF


      N=1
 32   READ (3,*,END=31) T(N)
      N=N+1
      GO TO 32
 31   CONTINUE
      N=N-1      
      
      CLOSE(UNIT=3)
      
      IF(N.gt.101) then
      print *,'The number of wavelength points >101, stop'
      stop
      end if
      IF(N.lt.2) then
      print *,'The number of wavelength points < 2, stop'
      stop
      end if

      HT(1)=0.5*(T(2)-T(1))
      HT(N)=0.5*(T(N)-T(N-1))
      HH(1)=-100.
      HH(N)=-100.
      DO 10 I=2,N-1
         HT(I)=0.5*(T(I+1)-T(I-1))
         HH(I)=-100.
10    CONTINUE
      

      
      READ(4,*) Umin,Umax,NU
      READ(4,*) Vmin,Vmax,NV

      IF(NU.gt.101) then
      print *,'The number of maxima >101, stop'
      stop
      end if
      IF(NV.gt.101) then
      print *,'The number of minima >101, stop'
      stop
      end if

      IF(NU.LT.1) Then
      print *,'The number of maxima < 1, stop'
      stop
      end if
      IF(NV.LT.1) Then
      print *,'The number of minima < 1, stop'
      stop
      end if

      Ulevels(1)=Umax
      IF (NU.lt.2) go to 25
      HU=(Umax-Umin)/FLOAT(NU-1)
      DO 20 I=1,NU-1
         ULEVELS(I+1)=Umax-FLOAT(I)*HU
20    CONTINUE

 25   continue
      Vlevels(1)=Vmax
      IF (NV.lt.2) go to 35
      HV=(Vmax-Vmin)/FLOAT(NV-1)
      DO 30 I=1,Nv-1
         VLEVELS(I+1)=Vmax-FLOAT(I)*HV
30    CONTINUE
35    continue
      CLOSE(UNIT=4)
      RETURN
      END




      SUBROUTINE COVG(XL0,XL2,XL4,COV,COV1,COV2,COV3,T,N)
C
C  COVG  evaluates: 
C
C  XL0,XL2,XL4 - spectral moments.
C
C  Covariance function and its four derivatives for a vector T of length N. 
C  It is saved in a vector COV; COV(1,...,N)=r(T), COV(N+1,...,2N)=r'(T), etc.
C  The vector COV should be of the length 5*N.
C
C  Covariance matrices COV1=r'(T-T), COV2=r''(T-T) and COV3=r'''(T-T) 
C  Dimension of  COV1, COV2  should be   N*N.
C
      DIMENSION T(1),A(501),TIME(501)
      DIMENSION COV(1),COV1(1),COV2(1),COV3(1)
      OPEN(UNIT=32,FILE='Cd0.in')
      OPEN(UNIT=33,FILE='Cd1.in')
      OPEN(UNIT=34,FILE='Cd2.in')
      OPEN(UNIT=35,FILE='Cd3.in')
      OPEN(UNIT=36,FILE='Cd4.in')
C
C   COV(Y(T),Y(0))
C

      NT=1
 12   READ (32,*,END=11) TIME(NT),A(NT)
      NT=NT+1
      GO TO 12
 11   CONTINUE
      NT=NT-1
      
      
      XL0=SPLE(NT,0,A,TIME)

      DO 10 I=1,N
      COV(I)=SPLE(NT,T(I),A,TIME)
10    CONTINUE

C     
C    DERIVATIVE  COV(Y(T),Y(0))
C

      NT=1
 22   READ (33,*,END=21) TIME(NT),A(NT)
      NT=NT+1
      GO TO 22
 21   CONTINUE
      NT=NT-1      

      II=0
      DO 20 I=1,N
      COV(I+N)=SPLE(NT,T(I),A,TIME)
      DO 20 J=1,N
      II=II+1
      T0=T(J)-T(I)
      TT=ABS(T0)
      COV1(II)=SPLE(NT,TT,A,TIME)
      IF (T0.LT.0.) COV1(II)=-COV1(II)
20    CONTINUE

C    2-DERIVATIVE  COV(Y(T),Y(0))

      NT=1
 32   READ (34,*,END=31) TIME(NT),A(NT)
      NT=NT+1
      GO TO 32
 31   CONTINUE
      NT=NT-1      

      II=0
      XL2=-SPLE(NT,0,A,TIME)

      DO 30 I=1,N
      COV(I+2*N)=SPLE(NT,T(I),A,TIME)
      DO 30 J=1,N
      II=II+1
      TT=ABS(T(J)-T(I))
      COV2(II)=SPLE(NT,TT,A,TIME)
30    CONTINUE

C    3-DERIVATIVE  COV(Y(T),Y(0))

            NT=1
 42   READ (35,*,END=41) TIME(NT),A(NT)
      NT=NT+1
      GO TO 42
 41   CONTINUE
      NT=NT-1      

      
      II=0
      DO 40 I=1,N
      COV(I+3*N)=SPLE(NT,T(I),A,TIME)
      DO 40 J=1,N
      II=II+1
      T0=T(J)-T(I)
      TT=ABS(T0)
      COV3(II)=SPLE(NT,TT,A,TIME)
      IF (T0.LT.0.) COV3(II)=-COV3(II)
40    CONTINUE



C    4-DERIVATIVE  COV(Y(T),Y(0))

      NT=1
 52   READ (36,*,END=51) TIME(NT),A(NT)
      NT=NT+1
      GO TO 52
 51   CONTINUE
      NT=NT-1
      
      XL4=SPLE(NT,0.,A,TIME)

      DO 50 I=1,N
      COV(I+4*N)=SPLE(NT,T(I),A,TIME)
50    CONTINUE
      CLOSE(UNIT=32)
      CLOSE(UNIT=33)
      CLOSE(UNIT=34)
      CLOSE(UNIT=35)
      CLOSE(UNIT=36)
      RETURN
      END

      SUBROUTINE INITINTEG(EPS0,IAC,NIT)
      dimension  INF(10),INFO(10)
      
      COMMON /RINT/   C,FC
      COMMON /EPS/    EPS,EPSS,CEPSS
      COMMON /INFC/   ISQ,INF,INFO
      OPEN(UNIT=1,FILE='accur.in')
      OPEN(UNIT=8,FILE='min.out')
      OPEN(UNIT=9,FILE='Max.out')
      OPEN(UNIT=10,FILE='Maxmin.out')
      OPEN(UNIT=11,FILE='Maxmin.log')

      READ(1,*) NIT,IAC,ISQ
      READ(1,*) EPS,EPSS,EPS0
      
      CLOSE (UNIT=1)
      
      FC=FI(C)-FI(-C)
      CEPSS=1.-EPSS

      RETURN
      END

