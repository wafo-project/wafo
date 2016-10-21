      PROGRAM down2cc

c     Author: Mats Frendahl
c     Date: Tue Jun 8 -93

c     Copyright, 1993, Mats Frendahl, Dept. of Math. Stat.
c     University of Lund, SWEDEN.
     
c     =======================================================
      
c     Calculates the maximum cycle count given the number of
c     downcrossings/downcrossing intensity.

c     The program reads crossint.dat which must contain the
c     levels  u  and the corresponding crossing intensities
c     in two column format.  u  must be increasing. The number
c     of slice levels are given in the file slice.dat.

c     Output is a file, out.dat, containing the slice levels,
c     the cycle maxima and cycle minima.

c     =======================================================
      
      DIMENSION u(50000),mu(50000),index(50000),out(500000,3)
      
      REAL u,mu,max_mu,out,delta
      
      OPEN(UNIT=1,FILE='crossint.dat')
      OPEN(UNIT=2,FILE='slice.dat')
      OPEN(UNIT=3,FILE='out.dat')

c     Read the number of slice levels.
      
      READ(2,*) number_of_slices

c     Close the file, will not be needed any more.
      
      CLOSE(UNIT=2)

c     Initiate, k, a counter for the number of points defining
c     the crossing intensities and max_mu, the maximum value
c     of the crossing intensity.
      
      k=1
      max_mu=0

c     Read the crossing intensity data, the levels  u  and the
c     crossing intensities  mu.
      
 1    READ(1,*,END=9) u(k),mu(k)
      IF (mu(k).GE.max_mu) max_mu=mu(k)
      k=k+1
      GO TO 1
 9    CONTINUE

c     Close the file, will not be needed any more.
      
      CLOSE(UNIT=1)

c     The number of crossing intensity points are one less.
      
      number_of_u=k-1

c     The slice distance is calculated.
      
      delta=max_mu/number_of_slices

c     Initiate, nc, the number of cycles.
      
      nc=0

c     Fix a slice level.
      
      DO 10 i=1,number_of_slices

c     Calculate the height of the slice level. OBS: We lower the
c     slice levels by delta/2 in order to a) really get a meaning-
c     full cycle if the number of slice levels = 1 and b) to get
c     closer to the bottom where the "tail-cycles" can be hard to get.
         
         h=i*delta-0.5*delta

c     Initiate, n, the number of crossings of the slice level.

         n=0

c     Run trough all  u  and find where the crossing intensity
c     crosses the slice level.
         
         DO 20 j=1,number_of_u-1

c     If the crossing intensity crosses the slice level (up or
c     down) add one to the number of crossings and mark the
c     index for the crossing.
            
            IF ((mu(j).LE.h.AND.h.LE.mu(j+1))
     *           .OR.(mu(j).GE.h.AND.h.GE.mu(j+1))) THEN
               n=n+1
               index(n)=j
            ENDIF
 20      CONTINUE

c     If the number of crossings are greater than 0, extract
c     the cycles by combining a cycle using the two closest
c     u's that cross the slice level. We take the mean of the
c     two u's near the crossing. We also save the slice level.
         
         IF (n.GT.0) THEN
            DO 30 k=1,n,2
               nc=nc+1
               out(nc,1)=h
               out(nc,2)=(u(index(k+1))+u(index(k+1)+1))/2
               out(nc,3)=(u(index(k))+u(index(k)+1))/2
 30         CONTINUE
         ENDIF

c     Pick a new slice level.
         
 10   CONTINUE

c     Write the slice levels and the cycle count to file.
      
      DO 40 l=1,nc
         WRITE(3,*) out(l,1),out(l,2),out(l,3)
 40   CONTINUE

c     Close remaining files
      
      CLOSE(UNIT=3)
      
      END
      
