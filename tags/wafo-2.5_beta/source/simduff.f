      PROGRAM oscillatorD
      REAL    dt,z,a,b,alf
      INTEGER number_of_points, seed
C
C     Open files.
C      
      OPEN(UNIT=1,FILE='simduff.in')
      OPEN(UNIT=4,FILE='out.dat')
C
      READ(1,*) number_of_points,dt,z,a,b,alf
      READ(1,*) seed
C

      cst=sqrt(2.)
      s  = sqrt(4*z)*exp(1./alf*alog(dt))

      if (alf.gt.1.99999) then
      X1    = gauss(seed)
      if (a.gt.0.01) X1=X1/sqrt(a)
      Xdot1 = gauss(seed)
       kp=1
       DO 10 i=2,number_of_points
         random=gauss(seed)
         X2    = X1+Xdot1*dt
         Xdot2 = Xdot1+(-b*X1-a*X1**3-2*z*Xdot1)
     1   *dt + s*random
         X1=X2
         Xdot1=Xdot2
         if (i/100.eq.kp) then
         write(4,*) i*dt,x2,xdot2
         kp=kp+1
         end if
 10   CONTINUE
         else

      X1    = stable(seed,alf)/cst
      if (a.gt.0.01) X1=X1/sqrt(a)
      Xdot1 = stable(seed,alf)/cst
      kp=1
      DO 50 i=2,number_of_points
         random=stable(seed,alf)/cst
         X2    = X1+Xdot1*dt
         Xdot2 = Xdot1+(-b*X1-a*X1**3-2*z*Xdot1)
     1   *dt + s*random
         X1=X2
         Xdot1=Xdot2
         if (i/100.eq.kp) then
         write(4,*) i*dt,x2,xdot2
         kp=kp+1
         end if

 50   CONTINUE
         endif
C
C     Close files.
C      
      CLOSE(UNIT=1)
      CLOSE(UNIT=4)
C
      STOP
      END
C
C      
      FUNCTION gauss(x)
      INTEGER x
      DATA pi/3.141592654/
      u1=rangen(x)
      u2=rangen(x)
      gauss=sqrt(-2*alog(u1))*cos(2*pi*u2)
      RETURN
      END
C     
C     Fuction to generate uniform random numbers. Note that the
C     seed  X  should be a 9 digit integer.
C     
      FUNCTION rangen(X)
      INTEGER A, P, X, B15, B16, XHI, XALO, LEFTLO, FHI, K
      DATA A/16807/,B15/32768/,B16/65536/,P/2147483647/
      XHI=X/B16
      XALO=(X-XHI*B16)*A
      LEFTLO=XALO/B16
      FHI=XHI*A+LEFTLO
      K=FHI/B15
      X=(((XALO-LEFTLO*B16)-P)+(FHI-K*B15)*B16)+K
      IF (X.LT.0) X=X+P
      rangen=FLOAT(X)*4.656612875E-10
      RETURN
      END


      function stable(ix,alf)
	REAL alf, xk, hh,w,v
	INTEGER ix
      data pi/3.141592654/
      v=pi*(rangen(ix)-0.5)
      xk=alf*v
      w=-alog(rangen(ix))/cos(v-xk)
      hh=w*cos(v)
      stable=sin(xk)*exp(-alog(hh)/alf)*w
      return
      end
