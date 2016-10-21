!****************************************************************************
!NB! if compilation complains about too many continuation lines extend it.
! 
!
!  modules:   GLOBALDATA, GLOBCONST, FUNCMOD , RIND    Version 1.0   
!
! Programs available in module RIND : 
! (NB! the GLOBALDATA  module is also used to transport the inputs)  
!
!
! SETDATA initializes global constants explicitly:
!   
! CALL SETDATA(RELEPS,ABSEPS,EPS2,xCutOff) 
!
!                   GLOBALDATA module :
!   RELEPS,ABSEPS = relative and absolute requested error 
!        EPS2  = if conditional variance is less it is considered as zero
!                i.e., the variable is considered deterministic 
!      xCutOff = 5 (standard deviations by default)
!     FxCutOff = FI(xCutOff) - FI( -xCutOff), truncation of the normal CDF.  
!
! INITDATA initializes global constants implicitly:
!
! CALL INITDATA (speed)    
!
!        speed = 1,2,...,9 (1=slowest and most accurate,9=fastest, 
!                           but less accurate)
!
! see the GLOBALDATA, FUNCMOD, GLOBALCONST 
!
!
!RINDD computes  E[Jacobian*Indicator|Condition]*f_{Xc}(xc(:,ix)) 
!
! where
!     "Indicator" = I{ H_lo(i) < X(i) < H_up(i), i=1:Nt+Nd }
!     "Jacobian"  = J(X(Nt+1),...,X(Nt+Nd+Nc)), special case is 
!     "Jacobian"  = |X(Nt+1)*...*X(Nt+Nd)|=|Xd(1)*Xd(2)..Xd(Nd)|
!     "condition" = Xc=xc(:,ix),  ix=1,...,Nx.
!     X = [Xt; Xd ;Xc], a stochastic vector of Multivariate Gaussian 
!         variables where Xt,Xd and Xc have the length Nt, Nd and Nc,
!         respectively. 
!         (Recommended limitations Nx, Nt<101, Nd<7 and NIT,Nc<11) 
! (RIND = Random Integration N Dimensions) 
!
!CALL RINDD(E,S,m,xc,Nt,indI,Blo,Bup,INFIN);
!
!        E = expectation/density as explained above size 1 x Nx        (out)
!        S = Covariance matrix of X=[Xt;Xd;Xc] size N x N (N=Nt+Nd+Nc) (in)
!        m = the expectation of X=[Xt;Xd;Xc]   size N x 1              (in)
!       xc = values to condition on            size Nc x Nx            (in)
!     indI = vector of indices to the different barriers in the        (in) 
!            indicator function,  length NI, where   NI = Nb+1 
!            (NB! restriction  indI(1)=0, indI(NI)=Nt+Nd )
!  Blo,Bup = Lower and upper barrier coefficients used to compute the   (in)
!            integration limits A and B, respectively. 
!            size  Mb x Nb. If  Mb<Nc+1 then
!            Blo(Mb+1:Nc+1,:) is assumed to be zero. 
!    INFIN = INTEGER, array of integration limits flags:  size 1 x Nb   (in)
!            if INFIN(I) < 0, Ith limits are (-infinity, infinity);
!            if INFIN(I) = 0, Ith limits are (-infinity, Hup(I)];
!            if INFIN(I) = 1, Ith limits are [Hlo(I), infinity);
!            if INFIN(I) = 2, Ith limits are [Hlo(I), Hup(I)].
! 
! The relation to the integration limits Hlo and Hup are as follows
!    IF INFIN(j)>=0,
!      IF INFIN(j)~=0,  A(i)=Blo(1,j)+Blo(2:Mb,j).'*xc(1:Mb-1,ix), 
!      IF INFIN(j)~=1,  B(i)=Bup(1,j)+Bup(2:Mb,j).'*xc(1:Mb-1,ix), 
!
!            where i=indI(j-1)+1:indI(j), j=1:NI-1, ix=1:Nx
!            Thus the integration limits may change with the conditional
!            variables.
!Example: 
! The indices, indI=[0 3 5 6], and coefficients Blo=[0 0 -1], Bup=[0 0 5],INFIN=[0 1 2] 
! means that   A = [-inf -inf -inf 0 0 -1]  B = [0 0 0 inf inf 5] 
!
! The GLOBALDATA and FUNCMOD modules are used to transport the inputs: 
!     SCIS = 1 Integrate all by SADAPT if NDIM<9 otherwise by KRBVRC (default)
!            2 Integrate all by SADAPT by Genz (1992) (Fast and reliable)
!            3 Integrate all by KRBVRC by Genz (1998) (Fast and reliable)
!            4 Integrate all by KROBOV by Genz (1992) (Fast and reliable)
!            5 Integrate all by RCRUDE by Genz (1992) (Reliable)
!            6 Integrate all by RLHCRUDE using MLHD and center of cell 
!            7 Integrate all by RLHCRUDE using LHD and center of cell
!            8 Integrate all by RLHCRUDE using MLHD and random point within the cell
!            9 Integrate all by RLHCRUDE using LHD and random point within the cell 
!
!jacobdef = defines the jacobian used (default jacobdef=0)
!
! (Recommended limitations Nx,Nt<101, Nd<7 and Nc<11)
! Also note that the size information have to be transferred to RINDD through the
! input arguments E,S,m,Nt,IndI,Blo,Bup and INFIN 
!
! if SCIS > 0 then you must initialize the random generator before you 
!  call rindd by the following lines:
!
!      call random_seed(SIZE=seed_size) 
!      allocate(seed(seed_size)) 
!      call random_seed(GET=seed(1:seed_size))  ! get current seed
!      seed(1)=seed1                            ! change seed
!      call random_seed(PUT=seed(1:seed_size)) 
!      deallocate(seed)
!
! For further description see the modules 
!
!
! References
! Podgorski et. al. (1999)
! "Exact distributions for apparent waves in irregular seas"
! Ocean Engineering                                                    (RINDD)
!
! R. Ambartzumian, A. Der Kiureghian, V. Ohanian and H.
! Sukiasian (1998)
! "Multinormal probabilities by sequential conditioned 
!  importance sampling: theory and application"                         (MVNFUN)
! Probabilistic Engineering Mechanics, Vol. 13, No 4. pp 299-308  
!
! Alan Genz (1992)
! 'Numerical Computation of Multivariate Normal Probabilites'           (MVNFUN)
! J. computational Graphical Statistics, Vol.1, pp 141--149
!
! Alan Genz and Koon-Shing Kwong (1999/2000?)
! 'Numerical Evaluation of Singular Multivariate Normal Distributions'  (MVNFUN,COVSRT)
! Submitted to Computational Statistics and Data analysis
!
!

! Tested on:  DIGITAL UNIX Fortran90 compiler
!             PC pentium II with Lahey Fortran90 compiler
!             Solaris with SunSoft F90 compiler Version 1.0.1.0  (21229283) 
! History:
! revised pab 18.02.2003
! -commented out all print statements.
! There is a bug in rcswap!!! TODO fix this.
! revised pab 05.05.2000
!  - found a bug in funcmod1 and funcmod when Nd>0
! revised pab 11.04.2000
!  - found a bug in mvnfun: The values of xd(I) was not computed correctly 
! revised pab 04.04.2000
!  This is based on Alan Genz (1992) and Alan Genz and Koon-Shing Kwong (1999/2000?)
!  Program for calculating Singular Multivariate probabilites.
!  Currently extended/updated for calculating multivariate Singular Expectations
!  Also updated from FORTRAN77 to FORTRAN90
!  - ADDED THL, EXLMS, EXINV
!*********************************************************************


      MODULE GLOBALDATA
      IMPLICIT NONE             
                      ! Constants determining accuracy of integration
                      !-----------------------------------------------
                      !if the conditional variance are less than: 
      DOUBLE PRECISION :: EPS2=1.d-10   !- EPS2, the variable is 
                                        !  considered deterministic 
      DOUBLE PRECISION :: EPS  = 1.d-5    ! = SQRT(EPS2)
      DOUBLE PRECISION :: EPS0 = 1.d-10  
      DOUBLE PRECISION :: XCEPS2=1.d-16 ! if Var(Xc) is less return NaN
      DOUBLE PRECISION :: AbsEps = 1.d-3 ! requested absolute error 
      DOUBLE PRECISION :: RelEps = 1.d-10 ! requested Relative error, i.e. if 
                                ! 3.0*STD(XIND)/XIND is less we accept the estimate
                                ! The following may be allocated outside RINDD
                                ! if one wants the coefficient of variation, i.e.
                                ! STDEV(XIND)/XIND when SCIS=2.  (NB: size Nx)  
      DOUBLE PRECISION :: fxcEpss=1.d-20 ! if less do not compute E(...|Xc)
      DOUBLE PRECISION :: xCutOff=5.d0  ! upper/lower truncation limit of the 
                                       ! normal CDF 
      DOUBLE PRECISION :: FxCutOff  = 0.99999942669686d0 
      DOUBLE PRECISION :: CFxCutOff = 5.733031438470704d-7       ! 1-FxCutOff, 
      DOUBLE PRECISION :: XSMALL=4.2D-16, XMAX=8.29287554336168D0 !cut off parameters to FI
    
!parameters controlling the performance and integration method
      INTEGER :: SCIS=1         !=1 Integrate by SADAPT if Ndim<9 KRBVRC otherwise 
                                !=2 Integrate all by SADAPT  (Fast and reliable)
                                !=3 Integrate all by KRBVRC  (Fastest and reliable)
                                !=4 Integrate all by KROBOV  (Fast and reliable)
                                !=5 Integrate all by RCRUDE  (Reliable)
                                !=6 Integrate all by RLHCRUDE using MLHD and center of cell 
                                !=7 Integrate all by RLHCRUDE using  LHD and center of cell
                                !=8 Integrate all by RLHCRUDE using MLHD and random point within the cell
                                !=9 Integrate all by RLHCRUDE using LHD and random point within the cell 
      INTEGER :: NSIMmax = 10000 ! maximum number of simulations
      INTEGER :: NSIMmin = 0    ! minimum number of simulations
      INTEGER :: rateLHD = 10   ! rateLhd*Nstoc = size of LHD matrix
      INTEGER :: jacobdef = 0   ! jacobdef determines the jacobian used
                                  
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: COV
      LOGICAL :: MLHD=.true.    ! Modified Latin Hypercube Design when sampling
                                ! WITH SCIS used in rindscis and mnormprb
      LOGICAL :: useMIDP=.true. ! use midpoint of cell instead randomly within cell
                                ! used on rindscis and mnormprb

      END MODULE GLOBALDATA

      MODULE GLOBALCONST        ! global constants
      DOUBLE PRECISION, PARAMETER :: SQTWPI1=3.9894228040143d-1 !=1/sqrt(2*pi)
      DOUBLE PRECISION, PARAMETER :: SQPI1=5.6418958354776d-1   !=1/sqrt(pi)
      DOUBLE PRECISION, PARAMETER :: SQPI= 1.77245385090552d0   !=sqrt(pi)
      DOUBLE PRECISION, PARAMETER :: SQTW=1.41421356237310d0    !=sqrt(2)
      DOUBLE PRECISION, PARAMETER :: SQTW1=0.70710678118655d0   !=1/sqrt(2)
      DOUBLE PRECISION, PARAMETER :: PI1=0.31830988618379d0     !=1/pi
      DOUBLE PRECISION, PARAMETER :: PI= 3.14159265358979D0     !=pi
      DOUBLE PRECISION, PARAMETER :: TWPI=6.28318530717958D0    !=2*pi
      DOUBLE PRECISION, PARAMETER :: SQTWPI=2.50662827463100D0  !=sqrt(2*pi)
      DOUBLE PRECISION, PARAMETER :: ONE=1.d0
      END MODULE GLOBALCONST

!
! FIMOD contains functions for calculating 1D and 2D Normal probabilites
!       and  1D expectations
      MODULE FIMOD
      IMPLICIT NONE

      INTERFACE FI
      MODULE PROCEDURE FI
      END INTERFACE 

      INTERFACE FIINV
      MODULE PROCEDURE FIINV
      END INTERFACE

      INTERFACE FI2
      MODULE PROCEDURE FI2
      END INTERFACE 

      INTERFACE FIINV2
      MODULE PROCEDURE FIINV2
      END INTERFACE
      
      INTERFACE  MVNLMS
      MODULE PROCEDURE  MVNLMS
      END INTERFACE

      INTERFACE BVU
      MODULE PROCEDURE BVU
      END INTERFACE

      INTERFACE BVNMVN
      MODULE PROCEDURE BVNMVN
      END INTERFACE
      
      INTERFACE NORM2DPRB
      MODULE PROCEDURE NORM2DPRB
      END INTERFACE
      
      INTERFACE THL
      MODULE PROCEDURE THL
      END INTERFACE

      INTERFACE GAUSINT
      MODULE PROCEDURE GAUSINT
      END INTERFACE

      INTERFACE GAUSINT2
      MODULE PROCEDURE GAUSINT2
      END INTERFACE
      
      INTERFACE EXLMS
      MODULE PROCEDURE EXLMS
      END INTERFACE

      INTERFACE EXINV
      MODULE PROCEDURE EXINV
      END INTERFACE


      CONTAINS

      FUNCTION FIINV(P) RESULT (VAL)
      use GLOBALDATA, only: XMAX,CFxCutOff
      IMPLICIT NONE
*
*	ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3
*
*	Produces the normal deviate Z corresponding to a given lower
*	tail area of P.
*
*	The hash sums below are the sums of the mantissas of the
*	coefficients.   They are included for use in checking
*	transcription.
*
      DOUBLE PRECISION, INTENT(in) :: P
      DOUBLE PRECISION :: VAL
!local variables
      DOUBLE PRECISION SPLIT1, SPLIT2, CONST1, CONST2, 
     *     A0, A1, A2, A3, A4, A5, A6, A7, B1, B2, B3, B4, B5, B6, B7, 
     *     C0, C1, C2, C3, C4, C5, C6, C7, D1, D2, D3, D4, D5, D6, D7, 
     *     E0, E1, E2, E3, E4, E5, E6, E7, F1, F2, F3, F4, F5, F6, F7, 
     *     Q, R
      PARAMETER ( SPLIT1 = 0.425d0, SPLIT2 = 5.d0,
     *            CONST1 = 0.180625D0, CONST2 = 1.6D0 )
*     
*     Coefficients for P close to 0.5
*     
      PARAMETER (
     *     A0 = 3.38713 28727 96366 6080D0,
     *     A1 = 1.33141 66789 17843 7745D+2,
     *     A2 = 1.97159 09503 06551 4427D+3,
     *     A3 = 1.37316 93765 50946 1125D+4,
     *     A4 = 4.59219 53931 54987 1457D+4,
     *     A5 = 6.72657 70927 00870 0853D+4,
     *     A6 = 3.34305 75583 58812 8105D+4,
     *     A7 = 2.50908 09287 30122 6727D+3,
     *     B1 = 4.23133 30701 60091 1252D+1,
     *     B2 = 6.87187 00749 20579 0830D+2,
     *     B3 = 5.39419 60214 24751 1077D+3,
     *     B4 = 2.12137 94301 58659 5867D+4,
     *     B5 = 3.93078 95800 09271 0610D+4,
     *     B6 = 2.87290 85735 72194 2674D+4,
     *     B7 = 5.22649 52788 52854 5610D+3 )
*     HASH SUM AB    55.88319 28806 14901 4439
*     
*     Coefficients for P not close to 0, 0.5 or 1.
*     
      PARAMETER (
     *     C0 = 1.42343 71107 49683 57734D0,
     *     C1 = 4.63033 78461 56545 29590D0,
     *     C2 = 5.76949 72214 60691 40550D0,
     *     C3 = 3.64784 83247 63204 60504D0,
     *     C4 = 1.27045 82524 52368 38258D0,
     *     C5 = 2.41780 72517 74506 11770D-1,
     *     C6 = 2.27238 44989 26918 45833D-2,
     *     C7 = 7.74545 01427 83414 07640D-4,
     *     D1 = 2.05319 16266 37758 82187D0,
     *     D2 = 1.67638 48301 83803 84940D0,
     *     D3 = 6.89767 33498 51000 04550D-1,
     *     D4 = 1.48103 97642 74800 74590D-1,
     *     D5 = 1.51986 66563 61645 71966D-2,
     *     D6 = 5.47593 80849 95344 94600D-4,
     *     D7 = 1.05075 00716 44416 84324D-9 )
*     HASH SUM CD    49.33206 50330 16102 89036
*
*	Coefficients for P near 0 or 1.
*
      PARAMETER (
     *     E0 = 6.65790 46435 01103 77720D0,
     *     E1 = 5.46378 49111 64114 36990D0,
     *     E2 = 1.78482 65399 17291 33580D0,
     *     E3 = 2.96560 57182 85048 91230D-1,
     *     E4 = 2.65321 89526 57612 30930D-2,
     *     E5 = 1.24266 09473 88078 43860D-3,
     *     E6 = 2.71155 55687 43487 57815D-5,
     *     E7 = 2.01033 43992 92288 13265D-7,
     *     F1 = 5.99832 20655 58879 37690D-1,
     *     F2 = 1.36929 88092 27358 05310D-1,
     *     F3 = 1.48753 61290 85061 48525D-2,
     *     F4 = 7.86869 13114 56132 59100D-4,
     *     F5 = 1.84631 83175 10054 68180D-5,
     *     F6 = 1.42151 17583 16445 88870D-7,
     *     F7 = 2.04426 31033 89939 78564D-15 )
*     HASH SUM EF    47.52583 31754 92896 71629
*     
      Q = ( P - 0.5d0)
      IF ( ABS(Q) .LE. SPLIT1 ) THEN ! Central range.
         R = CONST1 - Q*Q
         VAL = Q*( ( ( ((((A7*R + A6)*R + A5)*R + A4)*R + A3)
     *                  *R + A2 )*R + A1 )*R + A0 )
     *            /( ( ( ((((B7*R + B6)*R + B5)*R + B4)*R + B3)
     *                  *R + B2 )*R + B1 )*R + 1.d0 )
      ELSE ! near the endpoints
         R = MIN( P, 1.d0 - P )
         IF  (R .GT.0.d0) THEN ! ( 2.d0*R .GT. CFxCutOff) THEN ! R .GT.0.d0
            R = SQRT( -LOG(R) )
            IF ( R .LE. SPLIT2 ) THEN
               R = R - CONST2
               VAL = ( ( ( ((((C7*R + C6)*R + C5)*R + C4)*R + C3)
     *                      *R + C2 )*R + C1 )*R + C0 ) 
     *                /( ( ( ((((D7*R + D6)*R + D5)*R + D4)*R + D3)
     *                      *R + D2 )*R + D1 )*R + 1.d0 )
            ELSE
               R = R - SPLIT2
               VAL = ( ( ( ((((E7*R + E6)*R + E5)*R + E4)*R + E3)
     *                      *R + E2 )*R + E1 )*R + E0 )
     *                /( ( ( ((((F7*R + F6)*R + F5)*R + F4)*R + F3)
     *                      *R + F2 )*R + F1 )*R + 1.d0 )
            END IF
         ELSE
            VAL = 9.D0 !XMAX9.d0
         END IF
         IF ( Q .LT. 0.d0 ) VAL = - VAL
      END IF
      RETURN
      END FUNCTION FIINV     

      function FIINV2(p1) RESULT(value)
      use GLOBALDATA, only: XMAX,CFxCutOff
      DOUBLE PRECISION, INTENT(IN) :: P1
      DOUBLE PRECISION :: VALUE
! Local variables
      DOUBLE PRECISION, DIMENSION(4) :: A,B,C
      DOUBLE PRECISION, DIMENSION(2) :: D
      DOUBLE PRECISION :: P,Z,X=0.d0
      DOUBLE PRECISION, PARAMETER :: P0=0.7D0
       DOUBLE PRECISION, PARAMETER :: SQTWO=1.41421356237310d0    !=sqrt(2)
                                ! Coefficients in rational approximations.
      PARAMETER (A=(/ 0.886226899D0, -1.645349621D0, 0.914624893D0,
     &     -0.140543331D0 /))
      PARAMETER (B=(/  -2.118377725D0, 1.442710462D0, -0.329097515D0,  
     &     0.012229801D0/))
      PARAMETER (C=(/ -1.970840454D0, -1.624906493D0,  3.429567803D0, 
     &     1.641345311D0/))
      PARAMETER (D=(/  3.543889200D0,  1.637067800D0 /))
      !This function returns the inverse of the FI function
      ! It is translated from MATLAB where it
      ! originally evaluated the inverse of the error 
      ! function

      P = 2.d0*p1-1.d0
      if (abs(p) .LE. p0) THEN  ! Central range.
         z = p*p;
         x = p * (((a(4)*z+a(3))*z+a(2))*z+a(1)) / 
     &        ((((b(4)*z+b(3))*z+b(2))*z+b(1))*z+1.d0)
                                
      elseif (( p0 .LT. abs(p) ) .AND. 
     &        (abs(p) .LT. 1.d0-0.5d0*CFxCutOff)) THEN ! Near end points of range.
         z = sqrt(-log((1.d0-abs(p))/2.d0))
         x = sign((((c(4)*z+c(3))*z+c(2))*z+c(1))
     &        /((d(2)*z+d(1))*z+1.d0),p)
      
      else  !if (abs(p).GE.1.d0-.5d0*CFxCutOff) then  !    Exceptional case.
         value = sign(XMAX,p)
         return 
                                !elseif (abs(p).GT.1.d0) then
                                !x = NaN
                                !print *,'Error, FIINV p>1, ',p
                                !stop
      endif

! Two steps of Newton-Raphson correction to full accuracy.
! Without these steps, erfinv(y) would be about 3 times
! faster to compute, but accurate to only about 6 digits.

!x = x - (2.d0*FI(x*sqrt(2.d0))-1.d0 - p) / (2/sqrt(pi) * exp(-x^2));
!x = x - (2.d0*FI(x*sqrt(2.d0))-1.d0 - p) / (2/sqrt(pi) * exp(-x^2));

      value = SQTWO*x 
      return
      END FUNCTION FIINV2


                                ! *********************************     
      FUNCTION FI2(X) RESULT(value)
      USE GLOBALDATA, ONLY: XSMALL,XMAX
      IMPLICIT NONE
      DOUBLE PRECISION :: X,XX,value,z,y,del
      DOUBLE PRECISION :: xbreak=0.46875d0
!  MACHINE-DEPENDENT PARAMETERS  put into the GLOBALDATA module                                        
!------------------------------------------------------------------     
!     DOUBLE PRECISION, PARAMETER :: XSMALL=4.2D-16, XMAX=9.269D0        
!   This is a modification of a FORTRAN program by W. J. Cody,
!   Argonne National Laboratory, NETLIB/SPECFUN, March 19, 1990.
!   The original computation evaluated near-minimax approximations
!   from "Rational Chebyshev approximations for the error function"
!   by W. J. Cody, Math. Comp., 1969, PP. 631-638.
!   optimized for speed
!   For XX between -8 and 8
!   the maximum relative error is 2.3e-13 
!   and is obtained at XX=-8 
!   The absolute error is less than 1e-15
   
        
      XX=X*0.70710678118655d0
      y = ABS(XX)
      IF (y.LE.xbreak) THEN    
                                !  FI  for  |x| <= 0.46875*sqrt(2)

         IF (y .GT. XSMALL) THEN
            z = y * y
            value = 0.5d0+0.5d0*XX * ((((0.185777706184603153d0*z+
     &           3.16112374387056560d0)*z+1.13864154151050156d+2)*z+
     &           3.77485237685302021d+2)*z + 3.20937758913846947d+3)/ 
     &           ((((z+2.36012909523441209d+1)*z+2.44024637934444173d+2)
     &           *z+1.28261652607737228d+3)*z + 2.84423683343917062d+3)
         ELSE
            !z = 0.D0
            value=.5d0 !+0.56418958354776d0*XX
         ENDIF
       
         RETURN   
      END IF
   
      IF (y .LE. 4.d0) THEN
                                ! FI  for 0.46875*sqrt(2) <= |x| <= 4.0*sqrt(2)
         value = 0.5d0*((((((((2.15311535474403846d-8*y + 
     &        0.564188496988670089d0)*y+8.88314979438837594d0)*y+
     &        6.61191906371416295d+1)*y+2.98635138197400131d+2)*y+
     &        8.81952221241769090d+2)*y+1.71204761263407058d+3)*y+ 
     &        2.05107837782607147d+3)*y + 1.23033935479799725d+3)/ 
     &        ((((((((y+1.57449261107098347d+1)*y+
     &        1.17693950891312499d+2)*y+5.37181101862009858d+2)*y+
     &        1.62138957456669019d+3)*y+3.29079923573345963d+3)*y+
     &        4.36261909014324716d+3)*y+3.43936767414372164d+3)*y+ 
     &        1.23033935480374942d+3)
            
      ELSE ! FI  for |x| > 4.0*sqrt(2)
         IF (Y .GE. XMAX) THEN
            value=0.d0
            GOTO 300
         ELSE    
           z = 1.d0 / (y * y)
           value=z*(((((1.63153871373020978d-2*z+3.05326634961232344d-1)
     &           *z+3.60344899949804439d-1)*z+1.25781726111229246d-1)*z+
     &           1.60837851487422766d-2)*z + 6.58749161529837803d-4) / 
     &           (((((z+ 2.56852019228982242d0)*z+1.87295284992346047d0)
     &           *z+5.27905102951428412d-1)*z+6.05183413124413191d-2)*z+ 
     &           2.33520497626869185d-3)
           value = .5d0*(0.56418958354776d0 -  value) / y

         END IF
      END IF
      z = DBLE(FLOOR(y*16d0))/16d0
      del = (y-z)*(y+z)
      value = EXP(-z*z-del) * value

                                !   fix up for positive argument, FI
300   IF (XX .GT. xbreak) value=1.d0-value
      RETURN
                                !if (value > 1d0) then
                                !   value=1.d0
                                !else 
                                ! if (value<0.d0) then 
                                !    value=0.d0
                                ! end if
                                !end if
           
      END FUNCTION FI2
      FUNCTION FI( Z ) RESULT (VALUE)
      USE GLOBALDATA, ONLY : XMAX
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in) :: Z
      DOUBLE PRECISION :: VALUE
*     
*     Normal distribution probabilities accurate to 1.e-15.
*     Z = no. of standard deviations from the mean.
*     
*     Based upon algorithm 5666 for the error function, from:
*     Hart, J.F. et al, 'Computer Approximations', Wiley 1968
*     
*     Programmer: Alan Miller
*     
*     Latest revision - 30 March 1986
*     
      DOUBLE PRECISION :: P0, P1, P2, P3, P4, P5, P6, 
     *     Q0, Q1, Q2, Q3, Q4, Q5, Q6, Q7,
     *     P, EXPNTL, CUTOFF, ROOTPI, ZABS
      PARAMETER(
     *     P0 = 220.20 68679 12376 1D0,
     *     P1 = 221.21 35961 69931 1D0, 
     *     P2 = 112.07 92914 97870 9D0,
     *     P3 = 33.912 86607 83830 0D0,
     *     P4 = 6.3739 62203 53165 0D0,
     *     P5 = .70038 30644 43688 1D0, 
     *     P6 = .035262 49659 98910 9D0 )
      PARAMETER(
     *     Q0 = 440.41 37358 24752 2D0,
     *     Q1 = 793.82 65125 19948 4D0, 
     *     Q2 = 637.33 36333 78831 1D0,
     *     Q3 = 296.56 42487 79673 7D0, 
     *     Q4 = 86.780 73220 29460 8D0,
     *     Q5 = 16.064 17757 92069 5D0, 
     *     Q6 = 1.7556 67163 18264 2D0,
     *     Q7 = .088388 34764 83184 4D0 )
      PARAMETER( ROOTPI = 2.5066 28274 63100 1D0 )
      PARAMETER( CUTOFF = 7.0710 67811 86547 5D0 )
*     
      ZABS = ABS(Z)
*     
*     |Z| > 37  (or XMAX)
*     
      IF ( ZABS .GT. 37.D0 ) THEN
         P = 0.d0
      ELSE
*     
*     |Z| <= 37
*     
         EXPNTL = EXP( -ZABS**2/2D0 )
*     
*     |Z| < CUTOFF = 10/SQRT(2)
*     
         IF ( ZABS .LT. CUTOFF ) THEN
            P = EXPNTL*( (((((P6*ZABS + P5)*ZABS + P4)*ZABS + P3)*ZABS
     *           + P2)*ZABS + P1)*ZABS + P0)/(((((((Q7*ZABS + Q6)*ZABS
     *           + Q5)*ZABS + Q4)*ZABS + Q3)*ZABS + Q2)*ZABS + Q1)*ZABS
     *           + Q0 )
*     
*     |Z| >= CUTOFF.
*     
         ELSE
            P = EXPNTL/( ZABS + 1.d0/( ZABS + 2.d0/( ZABS + 3.d0/( ZABS 
     *                        + 4.d0/( ZABS + 0.65D0 ) ) ) ) )/ROOTPI
         END IF
      END IF
      IF ( Z .GT. 0.d0 ) P = 1.d0 - P
      VALUE = P
      RETURN
      END FUNCTION FI

      SUBROUTINE MVNLMS( A, B, INFIN, LOWER, UPPER )
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in) :: A, B
      DOUBLE PRECISION, INTENT(out) :: LOWER, UPPER
      INTEGER,INTENT(in) :: INFIN

      LOWER = 0.d0
      UPPER = 1.d0
      IF ( INFIN .LT. 0 ) RETURN
      IF ( INFIN .NE. 0 ) LOWER = FI(A)
      IF ( INFIN .NE. 1 ) UPPER = FI(B)
      RETURN
      END SUBROUTINE MVNLMS 
      
      FUNCTION BVNMVN( LOWER, UPPER, INFIN, CORREL ) RESULT (VAL)
      IMPLICIT NONE
*
*     A function for computing bivariate normal probabilities.
*
*  Parameters
*
*     LOWER  REAL, array of lower integration limits.
*     UPPER  REAL, array of upper integration limits.
*     INFIN  INTEGER, array of integration limits flags:
*            if INFIN(I) = 0, Ith limits are (-infinity, UPPER(I)];
*            if INFIN(I) = 1, Ith limits are [LOWER(I), infinity);
*            if INFIN(I) = 2, Ith limits are [LOWER(I), UPPER(I)].
*     CORREL REAL, correlation coefficient.
*
      
      DOUBLE PRECISION, DIMENSION(:), INTENT (IN) :: LOWER, UPPER
      DOUBLE PRECISION,               INTENT (IN) :: CORREL 
      INTEGER,          DIMENSION(:), INTENT (IN) :: INFIN
      DOUBLE PRECISION :: VAL
      IF ( INFIN(1) .EQ. 2  .AND. INFIN(2) .EQ. 2 ) THEN
         VAL =  BVU ( LOWER(1), LOWER(2), CORREL )
     +           - BVU ( UPPER(1), LOWER(2), CORREL )
     +           - BVU ( LOWER(1), UPPER(2), CORREL )
     +           + BVU ( UPPER(1), UPPER(2), CORREL )
      ELSE IF ( INFIN(1) .EQ. 2  .AND. INFIN(2) .EQ. 1 ) THEN
         VAL =  BVU ( LOWER(1), LOWER(2), CORREL )
     +           - BVU ( UPPER(1), LOWER(2), CORREL )
      ELSE IF ( INFIN(1) .EQ. 1  .AND. INFIN(2) .EQ. 2 ) THEN
         VAL =  BVU ( LOWER(1), LOWER(2), CORREL )
     +           - BVU ( LOWER(1), UPPER(2), CORREL )
      ELSE IF ( INFIN(1) .EQ. 2  .AND. INFIN(2) .EQ. 0 ) THEN
         VAL =  BVU ( -UPPER(1), -UPPER(2), CORREL )
     +           - BVU ( -LOWER(1), -UPPER(2), CORREL )
      ELSE IF ( INFIN(1) .EQ. 0  .AND. INFIN(2) .EQ. 2 ) THEN
         VAL =  BVU ( -UPPER(1), -UPPER(2), CORREL )
     +           - BVU ( -UPPER(1), -LOWER(2), CORREL )
      ELSE IF ( INFIN(1) .EQ. 1  .AND. INFIN(2) .EQ. 0 ) THEN
         VAL =  BVU ( LOWER(1), -UPPER(2), -CORREL )
      ELSE IF ( INFIN(1) .EQ. 0  .AND. INFIN(2) .EQ. 1 ) THEN
         VAL =  BVU ( -UPPER(1), LOWER(2), -CORREL )
      ELSE IF ( INFIN(1) .EQ. 1  .AND. INFIN(2) .EQ. 1 ) THEN
         VAL =  BVU ( LOWER(1), LOWER(2), CORREL )
      ELSE IF ( INFIN(1) .EQ. 0  .AND. INFIN(2) .EQ. 0 ) THEN
         VAL =  BVU ( -UPPER(1), -UPPER(2), CORREL )
      END IF
      END  FUNCTION BVNMVN
      FUNCTION BVU( SH, SK, R ) RESULT (VAL)
      IMPLICIT NONE
*
*     A function for computing bivariate normal probabilities.
*
*       Yihong Ge
*       Department of Computer Science and Electrical Engineering
*       Washington State University
*       Pullman, WA 99164-2752
*     and
*       Alan Genz
*       Department of Mathematics
*       Washington State University
*       Pullman, WA 99164-3113
*       Email : alangenz@wsu.edu
*
* BVN - calculate the probability that X is larger than SH and Y is
*       larger than SK. (to accuracy of 1e-16?)
*
* Parameters
*
*   SH  REAL, integration limit
*   SK  REAL, integration limit
*   R   REAL, correlation coefficient
*   LG  INTEGER, number of Gauss Rule Points and Weights
*
      DOUBLE PRECISION, INTENT(IN) :: SH, SK, R
      DOUBLE PRECISION  :: VAL
! Local variables
      DOUBLE PRECISION :: ZERO,HALF,ONE,TWO,TWOPI 
      INTEGER :: I, LG, NG
      PARAMETER ( ZERO = 0.D0,ONE=1D0,TWO=2D0,HALF=0.5D0)
      PARAMETER (TWOPI = 6.2831 85307 179586D0 ) 
      DOUBLE PRECISION, DIMENSION(10,3) :: X, W
      DOUBLE PRECISION :: AS, A, B, C, D, RS, XS
      DOUBLE PRECISION :: SN, ASR, H, K, BS, HS, HK
*     Gauss Legendre Points and Weights, N =  6
      DATA  ( W(I,1), X(I,1), I = 1,3) /
     *  0.1713244923791705D+00,-0.9324695142031522D+00,
     *  0.3607615730481384D+00,-0.6612093864662647D+00,
     *  0.4679139345726904D+00,-0.2386191860831970D+00/
*     Gauss Legendre Points and Weights, N = 12
      DATA ( W(I,2), X(I,2), I = 1,6) /
     *  0.4717533638651177D-01,-0.9815606342467191D+00,
     *  0.1069393259953183D+00,-0.9041172563704750D+00,
     *  0.1600783285433464D+00,-0.7699026741943050D+00,
     *  0.2031674267230659D+00,-0.5873179542866171D+00,
     *  0.2334925365383547D+00,-0.3678314989981802D+00,
     *  0.2491470458134029D+00,-0.1252334085114692D+00/
*     Gauss Legendre Points and Weights, N = 20
      DATA ( W(I,3), X(I,3), I = 1,10) /
     *  0.1761400713915212D-01,-0.9931285991850949D+00,
     *  0.4060142980038694D-01,-0.9639719272779138D+00,
     *  0.6267204833410906D-01,-0.9122344282513259D+00,
     *  0.8327674157670475D-01,-0.8391169718222188D+00,
     *  0.1019301198172404D+00,-0.7463319064601508D+00,
     *  0.1181945319615184D+00,-0.6360536807265150D+00,
     *  0.1316886384491766D+00,-0.5108670019508271D+00,
     *  0.1420961093183821D+00,-0.3737060887154196D+00,
     *  0.1491729864726037D+00,-0.2277858511416451D+00,
     *  0.1527533871307259D+00,-0.7652652113349733D-01/
      SAVE W, X
      IF ( ABS(R) .LT. 0.3D0 ) THEN
         NG = 1
         LG = 3
      ELSE IF ( ABS(R) .LT. 0.75D0 ) THEN
         NG = 2
         LG = 6
      ELSE 
         NG = 3
         LG = 10
      ENDIF
      H = SH
      K = SK 
      HK = H*K
      VAL = ZERO
      IF ( ABS(R) .LT. 0.925D0 ) THEN
         HS = ( H*H + K*K )*HALF
         ASR = ASIN(R)
         DO I = 1, LG
            SN  = SIN(ASR*( X(I,NG)+ONE )*HALF)
            VAL = VAL + W(I,NG)*EXP( ( SN*HK - HS )/( ONE - SN*SN ) )
            SN  = SIN(ASR*(-X(I,NG)+ONE )*HALF)
            VAL = VAL + W(I,NG)*EXP( ( SN*HK - HS )/( ONE - SN*SN ) )
         END DO
         VAL = VAL*ASR*HALF/TWOPI + FI(-H)*FI(-K) 
      ELSE
         IF ( R .LT. ZERO ) THEN
            K  = -K
            HK = -HK
         ENDIF
         IF ( ABS(R) .LT. ONE ) THEN
            AS  = ( ONE - R )*( ONE + R )
            A   = SQRT(AS)
            BS  = ( H - K )**2
            C   = ( 4D0 - HK )/8D0 
            D   = ( 12D0 - HK )/16D0
            VAL = A*EXP( -(BS/AS + HK)*HALF )
     +       *( ONE - C*(BS - AS)*(ONE - D*BS/5D0)/3D0 + C*D*AS*AS/5D0 )
            IF ( HK .GT. -160D0 ) THEN
               B   = SQRT(BS)
               VAL = VAL - EXP(-HK/2D0)*SQRT(TWOPI)*FI(-B/A)*B
     +                    *( ONE - C*BS*( ONE - D*BS/5D0 )/3D0 ) 
            ENDIF
            A = A/2D0
            DO I = 1, LG
               XS  = ( A*(X(I,NG)+ONE) )**2
               RS  = SQRT( ONE - XS )
               VAL = VAL + A*W(I,NG)*
     +             ( EXP( -BS/(2D0*XS) - HK/(ONE+RS) )/RS 
     +            - EXP( -(BS/XS+HK)/2D0 )*( ONE + C*XS*( ONE + D*XS )))
               XS  = AS*(ONE-X(I,NG))**2/4D0
               RS  = SQRT( ONE - XS )
               VAL = VAL + A*W(I,NG)*EXP( -(BS/XS + HK)/2D0 )
     +                    *( EXP( -HK*(ONE-RS)/(2D0*(ONE+RS)) )/RS 
     +                       - ( ONE + C*XS*( ONE + D*XS ) ) )
            END DO
            VAL = -VAL/TWOPI
         ENDIF
         IF ( R .GT. ZERO ) VAL =  VAL + FI( -MAX( H, K ) )
         IF ( R .LT. ZERO ) VAL = -VAL + MAX( ZERO, FI(-H)-FI(-K) )     
      ENDIF
      RETURN
      END FUNCTION BVU

      
      
      FUNCTION THL(H1, L1, GH) RESULT (T)                          
      USE GLOBALDATA, ONLY : ABSEPS
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in) :: H1, L1, GH
      DOUBLE PRECISION             :: T                                     
!                                                                       
!     Local constants                                                   
!     
      DOUBLE PRECISION, PARAMETER :: SQPI = 1.77245385090552d0   !=sqrt(pi)
      DOUBLE PRECISION, PARAMETER :: SQTWO= 1.41421356237310d0    !=sqrt(2)
      DOUBLE PRECISION, PARAMETER :: PI1  = 0.31830988618379d0     !=1/pi
      DOUBLE PRECISION, PARAMETER :: PI   = 3.14159265358979D0     !=pi
      DOUBLE PRECISION, PARAMETER :: ZERO = 0.D0
      DOUBLE PRECISION, PARAMETER :: QUART= 0.25D0
      DOUBLE PRECISION, PARAMETER :: HALF = 0.5D0
      DOUBLE PRECISION, PARAMETER :: ONE  = 1.D0
      DOUBLE PRECISION, PARAMETER :: TWO  = 2.D0
      DOUBLE PRECISION, PARAMETER :: EXPLIM=80.D0
      INTEGER,          PARAMETER :: NCk=7
      DOUBLE PRECISION, DIMENSION(NCk) :: Ck
! The following parameters gives absolute error <= 1e-7
      PARAMETER ( Ck =(/ 0.99999124672847d0, -0.99913342611804d0,
     &     0.98540829635096d0, -0.90305751774767d0, 0.66776957888038d0,
     &     -0.32268060187487d0, 0.07170680041655d0/))
! The following parameters gives absolute error <= 1e-11.
!      PARAMETER ( Ck =(/ 0.99999995579043d0, -0.99999111387630d0, 
!     &     0.99969933071259d0,-0.99596450623251d0, 0.97168814716964d0,  
!     &    -0.88105640681133d0, 0.67507517897079d0,-0.38534334229836d0,
!     &     0.13907128134639d0,-0.02317854687612d0/))
    
! Local variables
      INTEGER :: k
      DOUBLE PRECISION :: EPSL,H,L, H2, CON, CN, 
     &   L2, EX, W2, AP, S2, SP, S1, SN, SGN, H2L, Ik

! THL computes the integral
!             L
!   T(h,L) = int exp(h^2/2*(1-x^2))/(1+x^2) dx /(2*pi)
!             0
! 
!   GH = FI(h)
!
!     TLH is a controlled precision Fortran function 
!     which may be used to calculate    
!     bivariate normal probabilities, e.g. the CDF  for   
!     two normal variates X and Y whose correlation is R:
!
!     Prob(X<h,Y<k;R)=(FI(h)+FI(k))/2-T(h,g(h,k,R))-T(k,g(k,h,R))-J(h,k)/2
!  where 
!       J(h,k)   = 0     if h*k>0 or hk=0 and h+k>=0
!                  1     otherwise
!       g(h,k,R) = (h/k-R)/sqrt(1-R^2)
!
!     The accuracy is specified with ABSEPS.  
!
! Some properties of T(H,L):
! T(h,L) = 1/4-(FI(h)-1/2)*(FI(h*L)-1/2)-T(h*L,1/L)
! T(h,L) = T(-h,L)=-T(h,-L)
! T(h,0) = 0
! T(0,L) = atan(L)/(2*pi)
! T(h,1) = FI(h)*(1-FI(h))/2

! References
!  Jagdish K. Patel and Campbell B. Read (1982)
! "Handbook of the normal distribution",
! marcel dekker inc,New York - Basel, Vol. 40, pp 293--300       
!  
! Tested on: Matlab 5.1                                              
! History:
! revised pab 19.04.2000
!   - found a bug when L<-1, now fixed
! pab 20.03.2000
! - fixed some errors, optimized the splitting between Owen and Borth algorithms.
! - added Daley (1974) approximation
! pab 03.03.2000
!   - forgot to pass the calculated value for T on to output 
!     => this is probably why  THL gave floating invalid earlier
! by pab 19.08.199

!      print *, 'THL enter,L,H',L1,H1
      T    = ZERO
      EPSL = ABSEPS*1.d-1
      L    = L1
      H    = H1
      if (abs(L).LE.EPSL) RETURN ! T(H,0)=0
      
      if (abs(H).LE.EPSL) THEN   !T(0,L)
         T=ATAN(L)*HALF*PI1
         return
      endif
      if (abs(abs(L)-ONE).LE.EPSL) THEN ! T(H,1)
         T=HALF*GH*(ONE-GH)*SIGN(ONE,L)
         return
      end if   
      !print *, 'THL still here'
      SGN = ONE
      if (ABS(L).GT. ONE) THEN  ! make sure abs(L)<=1 to avoid numerical problems
                                ! using the identities:
                                ! T(h,L)=-T(h,-L)
                                ! T(h,L)=1/4+(1/2-FI(h))*(FI(h*L)-1/2)-T(h*L,1/L)
         IF (L.GE.0) SGN =-ONE
         H   = H * L                                                       
         L   = ONE /ABS( L)
         T   = (HALF-GH)*(FI(H)-HALF)-SGN*QUART
      endif

      H2 = H * H * HALF                                                                                                                                                                    
      IF (H2 .GE. EXPLIM) RETURN  ! THL is found to the desired accuracy
   
      IF (EPSL.GE.5.5d-5) THEN !Daley (1974) approximation
          SN=ATAN(L)
          H2L=L*L*H2*TWO
          H2L=H2L*H2L
          T = T + SGN*HALF*PI1*SN*EXP(-L*H2/SN)*(ONE+0.00868d0*H2L)
          RETURN
      ENDIF
      EX = EXP(-H2) 
      CON = abs(TWO*PI*EPSL ) 

! The following division between Owens and Borths algorithm is optimal for nearly all cases of EPSL>=1e-8
! This division is done in order to avoid slow convergence when L approx +/-1
!
!           Owens (1956) Algorithm:
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      IF (abs(H).LE.TWO.AND. 
     &     abs(L).LE.MAX(-0.72d0*abs(H)+1.936d0,HALF)) THEN

         L2 = L * L
         W2 = H2 * EX                                                       
         AP = ONE; SP=ONE
         S1 = ONE - EX;  S2=S1;    CN=S1
                                                      
         ! Alternating series: check if we can stop before SP=30 (usually SP is less than 8)         
         do while ((abs(CN*L).GT.CON).AND.( SP.LT.30.d0))
            SN =  SP                                                            
            SP =  SP + ONE                                                      
            S2 =  S2 - W2                                                       
            W2 =  W2 * H2 / SP                                                  
            AP = -AP * L2                                                        
            CN =  AP * S2 / (SN + SP)                                           
            S1 =  S1 + CN
         end do
!         print *,'THL, SP=',SP
                                !*, 'THL leaving 1.6' 
         T = T + SGN*(atan(L) - L*S1) *HALF*PI1                                                   
         return
      endif
! Otherwise
!
!           Borth (1974) Algorithm:
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~     
! approximate 1/(1+x^2) with a Tchebyshev polynomial of degree 2K
! and then integrate. Ck contains the coefficients of the polynomial.
! The absolute error is less than 1e-7.
        
      H2L=H2*L*L
      SN=SQTWO/H
      W2=exp(-H2L)*H*L/SQTWO ! safer when L approx 0
      Ik=SQPI*(FI(H*L)-HALF)
      
   
      CN=Ck(1)*Ik*SN
      S1=CN
      k=1
      !Alternating series: Check if we can stop before k=Nck
      do while ((abs(CN*EX).GT.CON) .AND. (k.LT.NCk)) 
         Ik=HALF*(DBLE(2*k-1)*Ik-W2)
         W2=W2*H2L
         SN=SN/H2
         k=k+1
         CN=Ck(k)*Ik*SN
         S1=S1+CN   
      end do	
!      print *,'THL K=',K                      
      T = T + SGN*EX*S1*HALF*PI1
      !print *, 'THL leaving last'
      RETURN                                                             
      END FUNCTION THL                                      

      FUNCTION NORM2DPRB(a1,b1,a2,b2,R) RESULT (prb)
      USE GLOBALDATA, ONLY : ABSEPS,CFxCutoff,xCutoff
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in) :: a1,b1,a2,b2,R
      DOUBLE PRECISION             :: prb
! Local variables
      DOUBLE PRECISION, PARAMETER :: PI1  = 0.31830988618379d0 !=1/pi
      DOUBLE PRECISION :: L1,L2,L3,L4,Ga1,Ga2,Gb1,Gb2,SQR,
     &     T11,T12,T21,T22,T31,T32,T41,T42,EPSL

!NORM2DPRB calculates the probability Prob(a1<X1<b1,a2<x2<b2)
!          for a bivariate Gaussian distribution
!          with correlation R
!
! References
!  Jagdish K. Patel and Campbell B. Read (1982)
! "Handbook of the normal distribution",
! marcel dekker inc,New York - Basel, Vol. 40, pp 293--300
! Tested on: Matlab 5.1


! History:
! pab 03.03.2000
!   - found  the error in THL which gave floating invalid
! by pab 19.08.1999
! NOTE!!! this does not work correctly yet. For some mysterious reason
! it crashes from time to time giving the message floating invalid
      EPSL=ABSEPS*1.d-1 
      if ((a1+EPSL.GE.b1).OR.(a2+EPSL.GE.b2)) THEN
         prb=0.d0
         return
      endif	
      Ga1=FI(a1)
      Gb1=FI(b1)
                                ! See if area is symmetric
                                ! The following may be optimized further
                                ! by checking against xCutOff
      if (abs(a1-a2).LE.EPSL) THEN
         Ga2=Ga1
      elseif (abs(a2-b1).LE.EPSL) THEN
         Ga2=Gb1
      else   
         Ga2=FI(a2)
      endif

      if (abs(b1-b2).LE.EPSL) THEN
         Gb2=Gb1
      elseif (abs(b2-a1).LE.EPSL) THEN
         Gb2=Ga1
      else   
         Gb2=FI(b2)
      endif	
      

      if (abs(R).LE.EPSL) then ! R=0
         prb=(Gb1-Ga1)*(Gb2-Ga2)
         GO TO 200
      endif

      if (abs(R)+EPSL.GE.1.d0) THEN
!         IF (abs(R).GT.1.d0+EPSL) PRINT *,'Warning: correlation, R=',R
         if (R.GT.0.d0) then    ! R = 1
            prb=min(Gb1,Gb2)+min(Ga1,Ga2)-min(Ga1,Gb2)-min(Gb1,Ga2)
         else   ! R =-1
            prb=max(Gb1+Gb2-1.d0,0.d0)+max(Ga1+Ga2-1.d0,0.d0)-
     &           max(Ga1+Gb2-1.d0,0.d0)-max(Gb1+Ga2-1.d0,0.d0)
         endif
         GO TO 200
      endif

      SQR = sqrt(1.d0-R*R)
        !print *, 'norm2d enter',ga1,gb1,ga2,gb2
      !print *, 'norm2d enter',a1,b1,a2,b2,R
      ! L1=PSI(b1,b2,R),! L3=PSI(b1,a2,R)
      if (abs(b1).LE.EPSL) then
         if (abs(b2).LE.EPSL) then
            L1=0.25d0+asin(R)*0.5d0*PI1
         else
            T11=SIGN(0.25d0,b2)
            T12=THL(b2,-R/SQR,Gb2)
            L1=0.5d0*Gb2+0.25d0-T11-T12+min(0.d0,SIGN(0.5d0,b2))
         endif
         if (abs(a2).LE.EPSL) then
            L3=0.25d0+asin(R)*0.5d0*PI1
         else
            T31=SIGN(0.25d0,a2)
            T32=THL(a2,-R/SQR,Ga2)
            L3=0.5d0*Ga2+0.25d0-T31-T32+min(0.d0,SIGN(0.5d0,a2))
         endif
      else
         !print *, 'norm2d'
         if (abs(b2).LE.EPSL) then
            T11=THL(b1,-R/SQR,Gb1)
            T12=SIGN(0.25d0,b1)
            L1=0.5d0*Gb1+0.25d0-T11-T12+0.5d0*min(0.d0,SIGN(1.d0,b1))
         else
            !print *, 'norm2d,b1,r,gb1',b1,(b2/b1-R)/SQR,gb1
            T11=THL(b1,(b2/b1-R)/SQR,Gb1)
            !print *, 'norm2d,T11',T11
            if (abs(b2-b1).LE.EPSL) then
               L1=Gb1-2.d0*T11
               T12=T11
            else  
               !print *,'norm2d T12'
               if (abs(b2+b1).LE.EPSL) then
                  T12=T11
               else
                  T12=THL(b2,(b1/b2-R)/SQR,Gb2)
               endif
               !print *,'norm2d T12',T12
               L1=0.5d0*(Gb1+Gb2)-T11-T12+min(0.d0,SIGN(0.5d0,b1*b2))
            endif
            !print *, 'norm2d,L1',L1
         endif
         !print *, 'norm2d'
         if (abs(a2).LE.EPSL ) then
            T32=SIGN(0.25d0,b1)
            T31=THL(b1,-R/SQR,Gb1)
            L3=0.5d0*Gb1+0.25d0-T31-T32+min(0.d0,SIGN(0.5d0,b1))
         else
           ! print *, 'norm2d,b1,r,gb1',b1,(a2/b1-R)/SQR,gb1
            T31=THL(b1,(a2/b1-R)/SQR,Gb1)
            !print *, 'norm2d,T31',T31
            if (abs(a2-b1).LE.EPSL) THEN
               L3=Ga1-2.d0*T31
               T32=T31
            else  
               IF (abs(a2+b1).LE.EPSL) THEN
                  T32=T31
               ELSE
                  T32=THL(a2,(b1/a2-R)/SQR,Ga2)
               ENDIF
               if ((abs(a2+b2).LE.EPSL).AND.(abs(a1+b1).LE.EPSL)) then
                  prb=2.d0*(T31+T32-T11-T12)+1.d0 ! OK
                                !prb=2.d0*(L1-L3)-1.d0
                  goto 200
               endif
               L3=0.5d0*(Gb1+Ga2)-T31-T32+min(0.d0,SIGN(0.5d0,b1*a2))        
            endif	
         endif   
      endif
      !print *, 'norm2d L1,L3',L1,L3
      !L2=PSI(a1,a2,R) L4=PSI(a1,b2,R)
      if (abs(a1).LE.EPSL) then
         if (abs(b2).LE.EPSL) then
            L4=0.25d0+asin(R)*0.5d0*PI1
         else
            T41=SIGN(0.25d0,b2)
            T42=THL(b2,-R/SQR,Gb2)
            L4=0.5d0*Gb2+0.25d0-T41-T42+min(0.d0,SIGN(0.5d0,b2))
         endif
         if (abs(a2).LE.EPSL) then
            L2=0.25d0+asin(R)*0.5d0*PI1
         else
            T21=SIGN(0.25d0,a2)
            T22=THL(a2,-R/SQR,Ga2)
            L2=0.5d0*Ga2+0.25d0-T21-T22+min(0.d0,SIGN(0.5d0,a2))
         endif
      else
         if (abs(b2).LE.EPSL ) then
            T41=THL(a1,-R/SQR,Ga1)
            T42=SIGN(0.25d0,a1)
            L4=0.5d0*Ga1+0.25d0-T41-T42+min(0.d0,SIGN(0.5d0,a1))
         else
            !print *, 'norm2d T41, L',(b2/a1-R)/SQR
            T41=THL(a1,(b2/a1-R)/SQR,Ga1)
            if (abs(b2-a1).LE.EPSL) then
               L4=Gb1-2.d0*T41
               T42=T41
            else    
               if (abs(b2+a1).LE.EPSL) then
                  T42=T41
               else
                  T42=THL(b2,(a1/b2-R)/SQR,Gb2)
               endif
               L4=0.5d0*(Ga1+Gb2)-T41-T42+min(0.d0,SIGN(0.5d0,a1*b2))
            endif	
         endif
         !print *, 'norm2d L4',L4
         if (abs(a2).LE.EPSL) then
            T22=SIGN(0.25d0,a1)
            T21=THL(a2,-R/SQR,Ga1)
            L2=0.5d0*Ga1+0.25d0-T21-T22+min(0.d0,SIGN(0.5d0,a1))
         else
            T21=THL(a1,(a2/a1-R)/SQR,Ga1)
            if (abs(a2-a1).LE.EPSL) then
               L2=Ga1-2.d0*T21
               T22=T21
            else   
               if (abs(a2+a1).LE.EPSL) then
                  T22=T21
               else
                  T22=THL(a2,(a1/a2-R)/SQR,Ga2)
               endif
               L2=0.5d0*(Ga1+Ga2)-T21-T22+min(0.d0,SIGN(0.5d0,a1*a2))
            endif	
         endif
         !print *, 'norm2d L2',L2
      endif

                                !L1,L2,L3,L4
      prb=L1+L2-L3-L4      
 200  if (prb > 1.d0) then ! fix up round off error
         prb=1.d0
      elseif (prb<0.d0) then 
            prb=0.d0
      end if
      !print *, 'norm2d leaving'
      return
      END FUNCTION NORM2DPRB
      
      FUNCTION GAUSINT (X1, X2, A, B, C, D) RESULT (value)               
      USE GLOBALDATA,ONLY:  xCutOff
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in) :: X1,X2,A,B,C,D
      DOUBLE PRECISION             :: value
!Local variables
      DOUBLE PRECISION             :: Y1,Y2,Y3
       DOUBLE PRECISION, PARAMETER :: SQTWOPI1=3.9894228040143d-1 !=1/sqrt(2*pi)
      ! Let  X  be standardized Gaussian variable, 
      ! i.e., X=N(0,1). The function calculate the
      !  following integral E[I(X1<X<X2)(A+BX)(C+DX)
      ! where I(X1<X<X2) is an indicator function of
      ! the set {X1<X<X2}. 
      IF (X1.GE.X2) THEN                                                 
         value = 0.d0                                                    
         RETURN                                                          
      ENDIF
      IF (ABS (X1) .GT.xCutOff) THEN                                           
         Y1 = 0.d0                                                         
      ELSE                                                               
         Y1 = (A * D+B * C + X1 * B * D) * EXP ( - 0.5d0 * X1 * X1)        
      ENDIF
      IF (ABS (X2) .GT.xCutOff) THEN                                           
         Y2 = 0.d0                                                         
      ELSE                                                               
         Y2 = (A * D+B * C + X2 * B * D) * EXP ( - 0.5d0 * X2 * X2)        
      ENDIF
      Y3 = (A * C + B * D) * (FI (X2) - FI (X1) )                        
      value = Y3 + SQTWOPI1 * (Y1 - Y2)                                      
      RETURN                                                             
      END FUNCTION GAUSINT
      
      
      FUNCTION GAUSINT2 (X1, X2, A, B) RESULT (value)               
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in) :: X1,X2,A,B
      DOUBLE PRECISION             :: value
! Local variables
      DOUBLE PRECISION             :: X0,Y0,Y1,Y2
      DOUBLE PRECISION, PARAMETER :: SQTWOPI1=3.9894228040143d-1 !=1/sqrt(2*pi)
      ! Let  X  be standardized Gaussian variable, 
      ! i.e., X=N(0,1). The function calculate the
      !  following integral E[I(X1<X<X2)ABS(A+BX)]
      ! where I(X1<X<X2) is an indicator function of
      ! the set {X1<X<X2}. 
      IF (X1.GE.X2) THEN                                                 
         value = 0.d0                                                    
         RETURN                                                          
      ENDIF
      IF (ABS(B).EQ.0.d0) THEN
         value = ABS(A)*(FI(X2)-FI(X1))
         RETURN
      ENDIF
     
      Y1 = -A*FI(X1)+SQTWOPI1*B*EXP(-0.5d0*X1*X1)
      Y2 = A*FI(X2)-SQTWOPI1*B*EXP(-0.5d0*X2*X2)
      IF ((B*X1.LT.-A).AND.(-A.LT.B*X2))THEN
         X0 = -A/B  
         Y0 = 2.d0*(A*FI(X0)-SQTWOPI1*B*EXP(-0.5d0*X0*X0))
         value=ABS(Y2-Y1-Y0)
      ELSE
         value=ABS(Y1+Y2) 
      ENDIF                                       
      RETURN                                                             
      END FUNCTION GAUSINT2

      SUBROUTINE EXLMS(A, X1, X2, INFIN, LOWER, UPPER, Ca,Pa)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in) :: A, X1, X2
      DOUBLE PRECISION, INTENT(out) :: LOWER, UPPER,Ca,Pa
      INTEGER,INTENT(in) :: INFIN
      ! Local variables
      DOUBLE PRECISION :: P1
      DOUBLE PRECISION, PARAMETER :: AMAX = 5D15
      
      ! Let  X  be standardized Gaussian variable, 
      ! i.e., X=N(0,1). The function calculate the
      !  following integral E[I(X1<X<X2)ABS(A+X)]/C(A) 
      ! where I(X1<X<X2) is an indicator function of
      ! the set {X1<X<X2} and C(A) is a normalization factor
      ! i.e. E(I(-inf<X<inf)*ABS(A+X)) = C(A) 
      ! 
      ! Pa = probability at the inflection point of the CDF
      ! Ca = C(A) normalization parameter
      ! A  = location parameter of the inflection point
      ! X1,X2 = integration limits < infinity
      LOWER = 0.d0
      UPPER = 1.d0
      
      P1 = EXFUN(-A,A)
      Ca = A+2D0*P1
      Pa = P1/Ca 
      IF ( INFIN .LT. 0 ) RETURN
      IF ( INFIN .NE. 0 ) THEN
         IF (X1.LE.-A) THEN
            LOWER = EXFUN(X1,A)/Ca
         ELSE
            LOWER = 1D0-EXFUN(-X1,-A)/Ca
         ENDIF
      ENDIF
      IF ( INFIN .NE. 1 ) THEN
         IF (X2.LE.-A) THEN
            UPPER = EXFUN(X2,A)/Ca
         ELSE
            UPPER = 1D0-EXFUN(-X2,-A)/Ca
         ENDIF
      ENDIF
      
      RETURN
      CONTAINS
      FUNCTION EXFUN(X,A) RESULT (P)
      DOUBLE PRECISION, INTENT(in) :: X,A
      DOUBLE PRECISION :: P
      DOUBLE PRECISION, PARAMETER :: SQTWOPI1=3.9894228040143d-1 !=1/sqrt(2*pi)
      P = EXP(-X*X*0.5D0)*SQTWOPI1-A*FI(X)
      END FUNCTION EXFUN
      END SUBROUTINE EXLMS 


      FUNCTION EXINV(P,A,Ca,Pa) RESULT (VAL)
!EXINV calculates the inverse of the CDF of abs(x+A)*exp(-x^2/2)/SQRT(2*pi)/C(A)
!      where C(A) is a normalization parameter depending on A
!
!   CALL:  val   = exinv(P,A,Ca,Pa)
! 
!     val  = quantiles
!     p    = probabilites
!     A    = location parameter
!     Ca   = normalization parameter
!     Pa   = excdf(-A,A) probability at the inflection point of the CDF
      double precision, intent(in) :: P,A,Ca,Pa
      double precision :: val
! local variables
      double precision, parameter :: amax = 5.D15
      double precision, parameter :: epsl = 5.D-15
      double precision, parameter :: xmax = 8
      double precision :: P1,Xk,Ak,Zk,SGN

!      if (P<0.D0.OR.P.GT.1.D0) PRINT *,'warning P<0 or P>1'

! The inverse cdf of 0 is -inf, and the inverse cdf of 1 is inf.  
      if (P.LE.EPSL.OR.P+EPSL.GE.1.D0) THEN
         VAL = SIGN(xmax,P-0.5D0)
         return
      endif
      Ak = ABS(A)
      if (EPSL.LT.Ak .AND. Ak.LT.amax) THEN
         IF (ABS(p-Pa).LE.EPSL) THEN
            VAL = SIGN(MIN(Ak,xmax),-A)
            RETURN
         ENDIF
         IF (Ak.LT.1D-2) THEN   ! starting guess always less than 0.2 from the true value
            IF (P.GE.0.5D0) THEN 
               xk = SQRT(-2D0*log(2D0*(1D0-P)))
            ELSE
               xk = -SQRT(-2D0*log(2D0*P))
            ENDIF
         ELSE
            xk = FIINV(P)       ! starting guess always less than 0.8 from the true value
            ! Modify starting guess if possible in order to speed up Newtons method
            IF (1D-3.LE.P.AND. P.LE.0.99D0.AND. 
     &           3.5.LE.Ak.AND.Ak.LE.1D3 ) THEN
               SGN = SIGN(1.d0,-A)
               Zk = xk*SGN
               xk = SGN*(Zk+((1D0/(64.9495D0*Ak-178.3191D0)-0.02D0/Ak)*
     &              Zk+1D0/(-0.99679234298211D0*Ak-0.07195350071872D0))/
     &            (Zk/(-1.48430620263825D0*Ak-0.33340759016175D0)+1D0))
            ELSEIF ((P.LT.1D-3.AND.A.LE.-3.5D0).OR. 
     &              (3.5D0.LE.A.AND.P.GT.0.99D0)) THEN
               SGN = SIGN(1.d0,-A)
               Zk = xk*SGN
               P1 = -2.00126182192701D0*Ak-2.57306603933111D0
               xk = SGN*Zk*(1D0+
     &              P1/((-0.99179258785909D0*Ak-0.21359746002397D0)*
     &              (Zk+P1)))
            ENDIF
         ENDIF
         ! Check if the starting guess is on the correct side of the inflection point
         IF (xk.LE.-A .AND. P.GT.Pa) xk = 1.D-2-A 
         IF (xk.GE.-A .AND. P.LT.Pa) xk = -1.D-2-A

      
         IF (P.LT.Pa) THEN
            VAL = funca(xk,A,P*Ca);
         ELSE  ! exploit the symmetry of the CDF
            VAL = -funca(-xk,-A,(1.D0-P)*Ca)
         ENDIF
         
      ELSEIF (ABS(A).LE.EPSL) THEN
         IF (P>=0.5D0) THEN
            VAL = SQRT(-2D0*log(2D0*(1.D0-P)))
         ELSE
            VAL = -SQRT(-2D0*log(2D0*P));
         ENDIF
      ELSE  ! ABS(A) > AMAX
         VAL = FIINV(P)
      ENDIF
      !CALL EXLMS(A,0.d0,VAL,0,ak,P1,zk,sgn) 
      !If (ABS(p-P1).GT.0.0001) PRINT *,'excdf(x,a)-p',p-P1 
      RETURN
      
      CONTAINS

      function funca(xk0,ak,CaP) RESULT (xk)
      double precision, intent(in) :: xk0,ak,CaP ! =Ca*P
      DOUBLE PRECISION :: xk
!Local variables                                
      INTEGER,          PARAMETER :: ixmax = 25
      double precision, parameter :: crit = 7.1D-08 ! = sqrt(1e-15)
      double precision, parameter :: SQTWOPI1 = 0.39894228040143D0 !=1/SQRT(2*pi)
      double precision, parameter :: SQTWOPI = 2.50662827463100D0 !=SQRT(2*pi)
      INTEGER :: IX
      DOUBLE PRECISION :: H,H1,tmp0,tmp1,XNEW
      ! Newton's Method or Fixed point iteration to find the inverse of the EXCDF.
      ! Assumption: xk0 < -ak and xk < -ak      
      ! Permit no more than IXMAX iterations.
      IX = 0
      H  = 1.D0
      xk = xk0    ! starting guess for the iteration
     
     
                                ! Break out of the iteration loop for the following:
                                !  1) The last update is very small (compared to x).
                                !  2) The last update is very small (compared to sqrt(eps)=crit).
                                !  3) There are more than 15 iterations. This should NEVER happen. 
      IF (.TRUE..OR.ABS(ak).LT.1.D-2) THEN 
      ! Newton's method
      !~~~~~~~~~~~~~~~~~
      DO WHILE( ABS(H).GT.MIN(crit*ABS(xk),crit).AND.IX.LT.IXMAX) 
         
         IX = IX+1   
                                !print *,'Iteration ',IX
         
         tmp0  = FI(xk)
         tmp1  = EXP(-xk*xk*0.5D0)*SQTWOPI1 ! =normpdf(x)
         H1 = (tmp1-ak*tmp0-CaP)/(ABS(xk+ak)*tmp1)
         H  = DSIGN(MIN(ABS(H1),0.7D0/DBLE(IX)),H1) ! Only allow smaller and smaller steps
         
         xnew = xk - H
                                ! Make sure that the current guess is less than -a.
                                ! When Newton's Method suggests steps that lead to -a guesses
                                ! take a step 9/10ths of the way to -a:
         IF (xnew.GT.-ak-crit) THEN
            xnew = (xk - 9.D0*ak)*1D-1;
            H    = xnew - xk
         ENDIF
         xk = xnew
      END DO
      ELSE                      ! FIXED POINT iteration
                                !~~~~~~~~~~~~~~~~~~~~~~~
         DO WHILE (ABS(H).GT.MIN(crit*ABS(xk),crit).AND.IX.LT.IXMAX)                                
            IX   = IX+1   
            tmp0 = SQTWOPI1*EXP(-xk*xk*0.5D0)/FI(xk)
            tmp1 = -2.D0*LOG(SQTWOPI*CaP*tmp0/(tmp0-ak))
            SGN  = sign(1.D0,tmp1)
            xnew = -SQRT(SGN*tmp1)*SGN
                                ! Make sure that the current guess is less than -a.
                                ! When this method suggests steps that lead to -a guesses
                                ! take a step 9/10ths of the way to -a:
            IF (xnew.GT.-ak-crit) xnew = (xk - 9.D0*ak)*1.D-1

            H  = xnew - xk
            xk = xnew
         END DO
      ENDIF

      !print *,'EXINV total number of iterations ',IX
      if (IX.GE.IXMAX) THEN 
!         print *, 'Warning: EXINV did not converge. Cap=',Cap
!         print *, 'The last step was:  ', h, ' value=,',xk,' ak=',ak
      endif
      return
      END FUNCTION FUNCA
      END FUNCTION EXINV

      END MODULE FIMOD


      MODULE FUNCMOD         ! FUNCTION module containing constants transfeered to mvnfun
      IMPLICIT NONE     
      
      INTERFACE
         SUBROUTINE MVNLMS( A, B, INFIN, LOWER, UPPER )
         DOUBLE PRECISION, INTENT(in) :: A, B
         DOUBLE PRECISION, INTENT(out) :: LOWER, UPPER
         INTEGER,INTENT(in) :: INFIN
         END SUBROUTINE MVNLMS 
      END INTERFACE

!      INTERFACE
!         FUNCTION FI( Z ) RESULT (VALUE)
!         DOUBLE PRECISION, INTENT(in) :: Z
!         DOUBLE PRECISION :: VALUE
!         END FUNCTION FI
!      END INTERFACE

      INTERFACE
         FUNCTION FIINV( Z ) RESULT (VALUE)
         DOUBLE PRECISION, INTENT(in) :: Z
         DOUBLE PRECISION :: VALUE
         END FUNCTION FIINV
      END INTERFACE
      
      DOUBLE PRECISION, DIMENSION(:  ), ALLOCATABLE :: COF,CDI  ! ChOlesky Factor and Cholesky DIagonal elements
      DOUBLE PRECISION, DIMENSION(:  ), ALLOCATABLE :: Cm,CDIXd,CmXd
      DOUBLE PRECISION, DIMENSION(:  ), ALLOCATABLE :: xd,xc,Y
      DOUBLE PRECISION, DIMENSION(:  ), ALLOCATABLE :: A,B      ! Integration limits
      INTEGER, DIMENSION(:  ), ALLOCATABLE :: INFI,INDEX1      
      INTEGER :: Nt,Nd,INFIS,INFISD    ! Size information
      INTEGER :: I0,NdleftN0,IJ0
      DOUBLE PRECISION :: E1,D1
! variables transfeered to mvnfun1:
      DOUBLE PRECISION :: VAL0,AA0,Ca0,Pa0
      
      INTERFACE  INITFUN
      MODULE PROCEDURE INITFUN
      END INTERFACE

      INTERFACE  MVNFUN
      MODULE PROCEDURE MVNFUN
      END INTERFACE
      
      INTERFACE  MVNFUN0
      MODULE PROCEDURE MVNFUN0
      END INTERFACE

      INTERFACE JACOB0
      MODULE PROCEDURE JACOB0
      END INTERFACE 
      

      CONTAINS
     

      FUNCTION JACOB0 (xd0,xc0) RESULT (value1) 
      USE GLOBALDATA, ONLY: jacobdef
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(:),INTENT(in) :: xd0 ,xc0
      DOUBLE PRECISION :: value1
      SELECT CASE (jacobdef)
      CASE (0)                  ! default
         value1 = ABS(PRODUCT(xd0))                             
      CASE (1)
         value1 = 1.d0
      CASE (2) 
         value1 = ABS(PRODUCT(xd0)*PRODUCT(xc0))  
      END SELECT                           
      RETURN                                                             
      END FUNCTION JACOB0

      SUBROUTINE INITFUN(VALUE,INFORM)
	USE JACOBMOD
      USE FIMOD
      USE GLOBALDATA, ONLY: EPS2,EPS,xCutOff
      IMPLICIT NONE
!      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(OUT) :: VALUE
      INTEGER, INTENT(out) :: INFORM
! local variables:
      INTEGER ::  N,NdleftO
      INTEGER :: I, J,  INFA, INFB
      DOUBLE PRECISION :: AI, BI
*     
*     Integrand subroutine
*
!INITFUN initilize the Multivariate Normal integrand function
! COF  - conditional sorted ChOlesky Factor of the covariance matrix (IN)
! CDI  - Cholesky DIagonal elements used to calculate the mean  
! Cm   - conditional mean of Xd and Xt given Xc, E(Xd,Xt|Xc)
! xd   - variables to the jacobian variable, need no initialization size Nd
! xc   - conditional variables (IN)
! INDEX1 - if INDEX1(I)>Nt then variable no. I is one of the Xd variables otherwise it is one of Xt
      !PRINT *,'Mvnfun,ndim',Ndim
      INFORM = 0
      VALUE  = 0.D0
      VAL0   = 1.D0

      NdleftN0 = Nd               ! Counter for number of Xd variables left
      
      IJ0  = 0
      I0   = 0
      INFA = 0
      INFB = 0
      N = Nt+Nd-INFIS-INFISD-1
      IF (INFISD.GT.0) THEN         ! replace xd with the mean
         I = Nt+Nd-INFIS
         J = NdleftN0-INFISD
                                !PRINT *,'Ndleft0',Ndleft0
         DO WHILE (NdleftN0.GT.J)
            IF (INDEX1(I).GT.Nt) THEN ! isXd
               xd (NdleftN0) =  Cm (I)
               NdleftN0=NdleftN0-1                  
            END IF
            I = I-1
         ENDDO
      ENDIF
      NdleftO = NdleftN0

      DO I = 1, N+1
         IJ0=IJ0+I
         IF (INFI(I).LT.0) GO TO 100            ! May have infinite int. Limits if Nd>0
         IF ( INFI(I) .NE. 0 ) THEN
            IF ( INFA .EQ. 1 ) THEN
               AI = MAX( AI, A(I) )
            ELSE
               AI = A(I) !MAX(A(I),-xCutOff) 
               INFA = 1
            END IF
         END IF
         IF ( INFI(I) .NE. 1 ) THEN
            IF ( INFB .EQ. 1 ) THEN
               BI = MIN( BI, B(I) )
            ELSE
               BI = B(I) !MIN( B(I),xCutOff) 
               INFB = 1
            END IF
         END IF
 100     IF (INDEX1(I).GT.Nt) THEN ! Save the mean for Xd
            CmXd(NdleftN0) =  Cm(I) 
            CDIXd(NdleftN0) = CDI(I)  
            NdleftN0 = NdleftN0-1
         END IF
    
         IF (I.EQ.N+1.OR.COF(IJ0+2).GT.0) THEN
            CALL MVNLMS( AI, BI, 2*INFA+INFB-1, D1, E1 )
            IF (D1.GE.E1) GOTO 200
            
            IF (Nd.GT.0.AND. NdleftO.LE.0) VAL0 = JACOB(xd,xc)  
            IF (I.EQ.N+1.AND.NdleftO.LE.0) THEN
               VALUE = (E1-D1)*VAL0
               GO TO 200
            ELSEIF (N.EQ.1.AND.NdleftO.LE.0) THEN
               IF ( ABS( COF(3) ) .GT. 0 ) THEN
                  D1 = SQRT( 1 + COF(2)**2 )
                  IF ( INFI(2) .NE. 0 ) A(2) = A(2)/D1
                  IF ( INFI(2) .NE. 1 ) B(2) = B(2)/D1
                  VALUE = BVNMVN( A, B,INFI,COF(2)/D1 )*VAL0
               ELSE             ! correlation=+/-1
                  VALUE = (E1-D1)*VAL0
               END IF
               GO TO 200 
            ENDIF
            VAL0 = VAL0*(E1-D1)
            I0   = I
            RETURN
         ENDIF
      ENDDO
      RETURN
 200  INFORM = 1
      
      RETURN
      END SUBROUTINE INITFUN

      FUNCTION MVNFUN0( Ndim, W ) RESULT (VAL)
	USE JACOBMOD
      USE FIMOD
      USE GLOBALDATA, ONLY: EPS2,EPS,xCutOff
      IMPLICIT NONE
      INTEGER, INTENT (in) :: Ndim
      DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: W
      DOUBLE PRECISION :: VAL
! local variables:
      INTEGER ::  N
      INTEGER :: I, J, IJ,  INFA, INFB,Ndleft,LK,JMX,IK
      DOUBLE PRECISION :: TMP, AI, BI, DI, EI,TMPOLD
      DOUBLE PRECISION, PARAMETER :: EPSL=1D-12
*     
*    The original  Integrand subroutine
*
!MVNFUN Multivariate Normal integrand function
! where the integrand is transformed from an integral 
! having integration limits A and B to an
! integral having  constant integration limits i.e.
!   B                                  1 
! int jacob(xd,xc)*f(xd,xt)dxt dxd = int F2(W) dW
!  A                                  0
!
! W    - new transformed integration variables, valid range 0..1
!        The vector must have the length Ndim returned from Covsrt
! COF  - conditional sorted ChOlesky Factor of the covariance matrix (IN)
! CDI  - Cholesky DIagonal elements used to calculate the mean  
! Cm   - conditional mean of Xd and Xt given Xc, E(Xd,Xt|Xc)
! xd   - variables to the jacobian variable, need no initialization size Nd
! xc   - conditional variables (IN)
! INDEX1 - if INDEX1(I)>Nt then variable No. I is one of the Xd variables otherwise it is one of Xt
      !PRINT *,'Mvnfun,ndim',Ndim
      VAL = 1.d0
      Ndleft = Nd               ! Counter for number of Xd variables left
      N=Nt+Nd-INFIS-1
      INFA = 0
      INFB = 0
      IK = 1                    ! Counter for Ndim 
      IJ = 0
      LK = 0                    ! Counter for LK constraints
      DO I = 1, N+1
         TMP = 0.d0
         JMX = MAX(I-IK,0)+1
         DO J = 1, I-JMX
            IJ = IJ + 1
            TMP = TMP + COF(IJ)*Y(J)
         END DO
         IJ = IJ+JMX
         IF (INFI(I).LT.0) GO TO 100            ! May have infinite int. Limits if Nd>0
         IF ( INFI(I) .NE. 0 ) THEN
            IF ( INFA .EQ. 1 ) THEN
               AI = MAX( AI, A(I) - TMP) ! - xCutOff*EPS)  ! added eps to make the int.lim. fuzzy
            ELSE
               AI = A(I) - TMP 
               INFA = 1
            END IF
         END IF
         IF ( INFI(I) .NE. 1 ) THEN
            IF ( INFB .EQ. 1 ) THEN
               BI = MIN( BI, B(I) - TMP) ! + xCutOff*EPS)  ! added eps to make the int.lim. fuzzy
            ELSE
               BI = B(I) - TMP  
               INFB = 1
            END IF
         END IF
 100     IF (LK.LT.1) THEN
            TMPOLD = TMP
         ELSEIF (INDEX1(I).GT.Nt) THEN ! Deterministic variable => replace xd with the mean
            xd(Ndleft) =  Cm(I)+TMP*CDI(I)  ! this is wrong
            Ndleft = Ndleft-1
         END IF
         IF ( I .EQ. N+1 .OR. COF(IJ+IK+1) .GT. EPSL ) THEN   
            CALL MVNLMS( AI, BI, 2*INFA+INFB-1, DI, EI )            
            IF ( DI .GE. EI ) THEN
               VAL = 0.d0
               RETURN
            ENDIF
            VAL = VAL*( EI - DI )
            IF ( I .LE. N .OR. INDEX1(I-LK).GT.Nt) THEN
               Y(IK) = FIINV( DI + W(IK)*( EI - DI ) )
               IF (INDEX1(I-LK).GT.Nt ) THEN
                  xd(Ndleft) = Cm(I-LK)+(TMPOLD+Y(IK))*CDI(I-LK)
                  Ndleft = Ndleft-1
               ENDIF
            ENDIF            
            LK   = 0
            IK   = IK + 1
            INFA = 0
            INFB = 0
         ELSE
            LK=LK+1
         END IF
      END DO
      
      IF (Nd.GT.0) VAL = VAL*jacob(xd,xc)      
      RETURN
      END FUNCTION MVNFUN0

      FUNCTION MVNFUN( Ndim, W ) RESULT (VAL)
	USE JACOBMOD
      USE FIMOD
      USE GLOBALDATA, ONLY: EPS2,EPS,xCutOff
      IMPLICIT NONE
      INTEGER, INTENT (in) :: Ndim
      DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: W
      DOUBLE PRECISION :: VAL
! local variables:
      INTEGER ::  N,I, J, IJ, INFA, INFB
      INTEGER ::  NdleftN,NdleftO ,JMX,IK
      DOUBLE PRECISION :: TMP, AI, BI, DI, EI
      DOUBLE PRECISION, PARAMETER :: EPSL=1D-12
*     
*     Integrand subroutine
*
!MVNFUN Multivariate Normal integrand function
! where the integrand is transformed from an integral 
! having integration limits A and B to an
! integral having  constant integration limits i.e.
!   B                                  1 
! int jacob(xd,xc)*f(xd,xt)dxt dxd = int F2(W) dW
!  A                                  0
!
! W    - new transformed integration variables, valid range 0..1
!        The vector must have the length Ndim returned from Covsrt
! COF  - conditional sorted ChOlesky Factor of the covariance matrix (IN)
! CDI  - Cholesky DIagonal elements used to calculate the mean  
! Cm   - conditional mean of Xd and Xt given Xc, E(Xd,Xt|Xc)
! xd   - variables to the jacobian variable, need no initialization size Nd
! xc   - conditional variables (IN)
! INDEX1 - if INDEX1(I)>Nt then variable No. I is one of the Xd variables otherwise it is one of Xt
      !PRINT *,'Mvnfun,ndim',Ndim
      
      
      N = Nt+Nd-INFIS-INFISD-1
      IK = 1                    ! Counter for Ndim 
      INFA = 0
      INFB = 0

      NdleftN = NdleftN0          ! Counter for number of Xd variables left
      IJ   = IJ0
      VAL  = VAL0
      NdleftO = Nd-INFISD
      Y(IK) = FIINV( D1 + W(IK)*( E1 - D1 ) )
      IF (NdleftO.GT.NdleftN ) THEN
         xd(NdleftN+1:NdleftO) = CmXd(NdleftN+1:NdleftO)+
     &        Y(IK)*CDIXd(NdleftN+1:NdleftO)
      ENDIF
      NdleftO = NdleftN
      IK = 2                    !=IK+1
            
     
      DO I = I0+1, N+1
         TMP = 0.d0
         JMX = MAX(I-IK,0)+1
         DO J = 1, I-JMX
            IJ = IJ + 1
            TMP = TMP + COF(IJ)*Y(J)  ! E(Y(IK) | Y(1),...Y(IK-1))/STD(Y(IK)|Y(1),,,,Y(IK-1))
         END DO
         IJ = IJ+JMX
         IF (INFI(I).LT.0) GO TO 100            ! May have infinite int. Limits if Nd>0
         IF ( INFI(I) .NE. 0 ) THEN
            IF ( INFA .EQ. 1 ) THEN
				IF (ABS(CDI(I)).GT.EPS2) THEN
					AI = MAX( AI, A(I) - TMP) 
				ELSE	
!					PRINT * ,'CDI(',I,') = ' , CDI(I)
!					PRINT * ,'AI = ', AI , ' AI2 =  ', A(I) - TMP
				ENDIF
            ELSE
               AI = A(I)-TMP 
				!AI = MAX(A(I) - TMP,-xCutOff) 
               INFA = 1
            END IF
            IF (INFB.EQ.1.AND.BI.LE.AI) GOTO 200
         END IF
         IF ( INFI(I) .NE. 1 ) THEN
            IF ( INFB .EQ. 1 ) THEN
				IF (ABS(CDI(I)).GT.EPS2) THEN
					BI = MIN( BI, B(I) - TMP)	
				ELSE
!					PRINT * ,'CDI(',I,') = ' , CDI(I)
!					PRINT * ,'BI = ', BI , ' BI2 =  ', B(I) - TMP
				ENDIF
            ELSE
               BI = B(I)-TMP 
				!BI = MIN(B(I) - TMP,xCutOff)  
               INFB = 1
            END IF
            IF (INFA.EQ.1.AND.BI.LE.AI) GOTO 200
         END IF
 100     IF (INDEX1(I).GT.Nt) THEN ! Save the mean of xd and Covariance diagonal element 
            CmXd(NdleftN) =  Cm(I)+TMP*CDI(I)    ! Conditional mean E(Xi|X1,..X)
            CDIXd(NdleftN) = CDI(I)              ! Covariance diagonal  
            NdleftN = NdleftN-1
         END IF
         IF ( I .EQ. N+1 .OR. COF(IJ+IK+1) .GT.0.D0 ) THEN   
            CALL MVNLMS( AI, BI, 2*INFA+INFB-1, DI, EI )            
            IF ( DI .GE. EI ) GO TO 200
            VAL = VAL*( EI - DI )
            
            IF ( I .LE. N .OR. NdleftN.LT.NdleftO) THEN
               Y(IK) = FIINV( DI + W(IK)*( EI - DI ) )
               IF (NdleftN.LT.NdleftO ) THEN
                  xd(NdleftN+1:NdleftO) = CmXd(NdleftN+1:NdleftO)+
     &                 Y(IK)*CDIXd(NdleftN+1:NdleftO)
                  NdleftO = NdleftN
               ENDIF
            ENDIF            
            IK   = IK + 1
            INFA = 0
            INFB = 0
         END IF
      END DO
      IF (Nd.GT.0) VAL=VAL*jacob(xd,xc)      
      RETURN
 200  VAL=0.d0
      RETURN
      END FUNCTION MVNFUN
      END MODULE FUNCMOD 

      MODULE FUNCMOD1            ! FUNCTION module containing constants transfeered to mvnfun
      IMPLICIT NONE     
      
      INTERFACE
         SUBROUTINE MVNLMS( A, B, INFIN, LOWER, UPPER )
         DOUBLE PRECISION, INTENT(in) :: A, B
         DOUBLE PRECISION, INTENT(out) :: LOWER, UPPER
         INTEGER,INTENT(in) :: INFIN
         END SUBROUTINE MVNLMS 
      END INTERFACE

!      INTERFACE
!         FUNCTION FI( Z ) RESULT (VALUE)
!         DOUBLE PRECISION, INTENT(in) :: Z
!         DOUBLE PRECISION :: VALUE
!         END FUNCTION FI
!      END INTERFACE

      INTERFACE
         FUNCTION FIINV( Z ) RESULT (VALUE)
         DOUBLE PRECISION, INTENT(in) :: Z
         DOUBLE PRECISION :: VALUE
         END FUNCTION FIINV
      END INTERFACE
      
      DOUBLE PRECISION, DIMENSION(:  ), ALLOCATABLE :: COF,CDI  ! ChOlesky Factor and Cholesky DIagonal elements
      DOUBLE PRECISION, DIMENSION(:  ), ALLOCATABLE :: Cm,CDIXd,CmXd
      DOUBLE PRECISION, DIMENSION(:  ), ALLOCATABLE :: xd,xc,Y
      DOUBLE PRECISION, DIMENSION(:  ), ALLOCATABLE :: A,B      ! Integration limits
      INTEGER, DIMENSION(:  ), ALLOCATABLE :: INFI,INDEX1      
      INTEGER :: Nt,Nd,INFIS,INFISD    ! Size information
      INTEGER :: I0,NdleftN0,IJ0
      DOUBLE PRECISION :: E1,D1
! variables transfeered to mvnfun1:
      DOUBLE PRECISION :: VAL0,AA0,Ca0,Pa0
      
      INTERFACE  INITFUN
      MODULE PROCEDURE INITFUN
      END INTERFACE

      INTERFACE  MVNFUN
      MODULE PROCEDURE MVNFUN
      END INTERFACE

      INTERFACE JACOB0
      MODULE PROCEDURE JACOB0
      END INTERFACE 
      

      CONTAINS
     

      FUNCTION JACOB0 (xd0,xc0) RESULT (value1) 
      USE GLOBALDATA, ONLY: jacobdef
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(:),INTENT(in) :: xd0 ,xc0
      DOUBLE PRECISION :: value1
      SELECT CASE (jacobdef)
      CASE (0)                  ! default
         value1 = ABS(PRODUCT(xd0))                             
      CASE (1)
         value1 = 1.d0
      CASE (2) 
         value1 = ABS(PRODUCT(xd0)*PRODUCT(xc0))  
      END SELECT                           
      RETURN                                                             
      END FUNCTION JACOB0

      SUBROUTINE INITFUN(VALUE,INFORM)
      USE FIMOD
	USE JACOBMOD
      USE GLOBALDATA, ONLY: EPS2,EPS,xCutOff
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(out) :: VALUE
      INTEGER, INTENT(out) :: INFORM
! local variables:
      INTEGER ::  N,NdleftO
      INTEGER :: I, J,  INFA, INFB
      DOUBLE PRECISION :: AI, BI
*     
*     Integrand initialiazation subroutine
*
!INITFUN initialize the Multivariate Normal integrand function
!   B                                  1 
! int ABS(xd(1)*xd(2)...xd(Nd))*f(xd,xt)dxt dxd = int F2(W) dW
!  A                                  0
!
! COF  - conditional sorted ChOlesky Factor of the covariance matrix (IN)
! CDI  - Cholesky DIagonal elements used to calculate the mean  
! Cm   - conditional mean of Xd and Xt given Xc, E(Xd,Xt|Xc)
! xd   - variables to the jacobian variable, need no initialization size Nd
! xc   - conditional variables (IN)
! INDEX1 - if INDEX1(I)>Nt then variable no. I is one of the Xd variables otherwise it is one of Xt
      !PRINT *,'Mvnfun,ndim',Ndim
      INFORM = 0
      VALUE  = 0.D0
      VAL0   = 1.d0

      CDIXd  = 0.d0
      CmXd   = 0.d0
      AA0    = 0.d0
      Pa0    = 0.5D0
      Ca0    = 1.D0

      INFA = 0
      INFB = 0
      IJ0  = 0
      I0   = 0
      NdleftN0 = Nd               ! Counter for number of Xd variables left
      N    = Nt+Nd-INFIS-INFISD-1

      IF (INFISD.GT.0) THEN     ! replace xd with the mean for all deterministic variables
         I = Nt+Nd-INFIS
         J = Nd-INFISD
                                !PRINT *,'Ndleft0',Ndleft0
         DO WHILE (NdleftN0.GT.J)
            IF (INDEX1(I).GT.Nt) THEN ! isXd
               VAL0 = VAL0*ABS(Cm (I))
               NdleftN0=NdleftN0-1                  
            END IF
            I=I-1
         ENDDO
      ENDIF

      NdleftO = NdleftN0

      DO I = 1, N+1
         IJ0=IJ0+I
         IF (INFI(I).LT.0) GO TO 100            ! May have infinite int. Limits if Nd>0
         IF ( INFI(I) .NE. 0 ) THEN
            IF ( INFA .EQ. 1 ) THEN
               AI = MAX( AI, A(I) ) 
            ELSE
               AI = MAX(A(I) ,-xCutOff) 
               INFA = 1
            END IF
         END IF
         IF ( INFI(I) .NE. 1 ) THEN
            IF ( INFB .EQ. 1 ) THEN
               BI = MIN( BI, B(I) ) 
            ELSE
               BI = MIN(B(I),xCutOff) 
               INFB = 1
            END IF
         END IF
 100     IF (I.GT.1.AND.INDEX1(I).GT.Nt) THEN ! Save the mean for Xd
            CmXd(NdleftN0)  = Cm(I)
            CDIXd(NdleftN0) = CDI(I)
            NdleftN0 = NdleftN0-1
         END IF
         IF ( I .EQ. N+1 .OR. COF(IJ0+2) .GT. 0.d0 ) THEN 
            IF ( INDEX1(1).GT.Nt) THEN
               AA0   = Cm(1)/CDI(1)
               CALL EXLMS(AA0,AI, BI, 2*INFA+INFB-1, D1, E1,Ca0,Pa0)
               VAL0 = VAL0*ABS(CDI(1)*Ca0)
            ELSE
               CALL MVNLMS( AI, BI, 2*INFA+INFB-1, D1, E1 )
            ENDIF
            IF (D1.GE.E1) GOTO 200
            
            IF (I.EQ.N+1.AND.NdleftO.LE.0) THEN
               VALUE = VAL0*(E1-D1)
               GO TO 200
            ELSEIF (N.EQ.1.AND.NdleftO.LE.0) THEN
               IF ( ABS( COF(3) ) .GT. 0 ) THEN
                  D1 = SQRT( 1 + COF(2)**2 )
                  IF ( INFI(2) .NE. 0 ) A(2) = A(2)/D1
                  IF ( INFI(2) .NE. 1 ) B(2) = B(2)/D1
                  VALUE = BVNMVN( A, B, INFI, COF(2)/D1 )*VAL0
               ELSE             ! correlation=+/-1
                  VALUE = VAL0*(E1-D1)
               END IF
               GO TO 200 
            ENDIF
            VAL0 = VAL0*(E1-D1)
            I0 = I
            !PRINT * ,'INITFUN VAL0,AA0,Ca0,Pa0',VAL0,AA0,Ca0,Pa0
            !PRINT *,'Cmxd',CmXd,'CDIXd',CDIXd
            RETURN
         ENDIF
      END DO
      RETURN
 200  INFORM = 1
      !PRINT * ,'INITFUN VALUE',VALUE
      RETURN
      END SUBROUTINE INITFUN

      FUNCTION MVNFUN( Ndim, W ) RESULT (VAL)
      USE FIMOD
      USE GLOBALDATA, ONLY: EPS2,EPS,xCutOff
      IMPLICIT NONE
      INTEGER, INTENT (in) :: Ndim
      DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: W
      DOUBLE PRECISION :: VAL
! local variables:
      INTEGER ::  N,NdleftN,NdleftO
      INTEGER :: I, J, IJ,  INFA, INFB,LK,JMX,IK
      DOUBLE PRECISION :: TMP, AI, BI, DI, EI,TMPOLD
      DOUBLE PRECISION :: AA,Ca ,Pa
      DOUBLE PRECISION, PARAMETER :: EPSL=1D-12
*     
*     Integrand subroutine
*
!MVNFUN Multivariate Normal integrand function
! where the integrand is transformed from an integral 
! having integration limits A and B to an
! integral having  constant integration limits i.e.
!   B                                               1 
! int ABS(xd(1)*xd(2)...xd(Nd))*f(xd,xt)dxt dxd = int F2(W) dW
!  A                                               0
!
! W    - new transformed integration variables, valid range 0..1
!        The vector must have the length Ndim returned from Covsrt
! COF  - conditional sorted ChOlesky Factor of the covariance matrix (IN)
! CDI  - Cholesky DIagonal elements used to calculate the mean  
! Cm   - conditional mean of Xd and Xt given Xc, E(Xd,Xt|Xc)
! xd   - variables to the jacobian variable, need no initialization size Nd
! xc   - conditional variables (IN)
! INDEX1 - if INDEX1(I)>Nt then variable No. I is one of the Xd variables otherwise it is one of Xt
      !PRINT *,'Mvnfun,ndim',Ndim
      VAL = VAL0
      IJ  = IJ0
      
      N = Nt+Nd-INFIS-INFISD-1
       
      IF (INDEX1(1).GT.Nt ) THEN
         Y(1) = EXINV( D1 + W(1)*( E1 - D1 ),AA0,Ca0,Pa0)
      ELSE
         Y(1) = FIINV( D1 + W(1)*( E1 - D1 ) )   
      ENDIF
      
      NdleftN = NdleftN0
      NdleftO = Nd-INFISD 

      DO I=NdleftN+1,NdleftO
         VAL = VAL*ABS(CmXd(I)+Y(1)*CDIXd(I))
      ENDDO
      !PRINT *,'IJ',IJ,'NdleftO',NdleftO,NdleftN
      NdleftO = NdleftN

      LK   = 0                  ! Counter for LK constraints
      IK   = 2                  ! Counter for Ndim 
      INFA = 0
      INFB = 0
      
     
      DO I = I0+1, N+1
         TMP = 0.d0
         JMX = MAX(I-IK,0)+1
         DO J = 1, I-JMX
            IJ = IJ + 1
            TMP = TMP + COF(IJ)*Y(J)
         END DO
         IJ = IJ+JMX
         IF (INFI(I).LT.0) GO TO 100            ! May have infinite int. Limits if Nd>0
         IF ( INFI(I) .NE. 0 ) THEN
            IF ( INFA .EQ. 1 ) THEN
               AI   = MAX( AI, A(I) - TMP) 
            ELSE
               AI   = MAX(A(I) - TMP,-xCutOff) 
               INFA = 1
            END IF
            IF (INFB.EQ.1.AND.BI.LE.AI) GOTO 200
         END IF
         IF ( INFI(I) .NE. 1 ) THEN
            IF ( INFB .EQ. 1 ) THEN
               BI   = MIN( BI, B(I) - TMP) 
            ELSE
               BI   = MIN(B(I) - TMP,xCutOff)  
               INFB = 1
            END IF
            IF (INFA.EQ.1.AND.BI.LE.AI) GOTO 200
         END IF
        
 100     IF (LK.LT.1) THEN
            TMPOLD = TMP
         ELSEIF   (INDEX1(I).GT.Nt) THEN ! Save the mean for Xd
            CmXd(NdleftN)  = Cm(I)+TMP*CDI(I)     ! Mean
            CDIXd(NdleftN) = CDI(I)
            !VAL = VAL*ABS(Cm(I)+TMP*CDI(I)) 
            NdleftN = NdleftN-1
         END IF
         IF ( I .EQ. N+1 .OR. COF(IJ+IK+1) .GT.0.D0 ) THEN   
            IF (INDEX1(I-LK).GT.Nt) THEN
               AA =  Cm(I-LK)/CDI(I-LK)+TMPOLD  ! location parameter of inflection point
               CALL EXLMS(AA,AI, BI, 2*INFA+INFB-1, DI, EI,Ca,Pa)
               
               IF ( DI .GE. EI ) GO TO 200
               VAL = VAL*ABS(Ca*CDI(I-LK))
               IF ( I .LE. N .OR. NdleftN.LT.NdleftO) THEN
                  Y(IK) = EXINV(DI + W(IK)*( EI - DI ),AA,Ca,Pa)
                  !PRINT * ,'Y',Y(IK)
                  DO J = NdleftN+1,NdleftO
                     VAL = VAL*ABS(CmXd(J)+Y(IK)*CDIXd(J))
                  ENDDO
                  NdleftO = NdleftN
                  !PRINT *,'AA=',AA,' Ca=',Ca,' Pa=',Pa
                  
               ENDIF   
               
            ELSE
               CALL MVNLMS( AI, BI, 2*INFA+INFB-1, DI, EI )            
               IF ( DI .GE. EI ) GO TO 200
               IF ( I .LE. N .OR. NdleftN.LT.NdleftO) THEN
                  Y(IK) = FIINV( DI + W(IK)*( EI - DI ) )
                  DO J = NdleftN+1,NdleftO
                     VAL = VAL*ABS(CmXd(J)+Y(IK)*CDIXd(J))
                  ENDDO
                  NdleftO = NdleftN
               ENDIF
            ENDIF            
            VAL = VAL*( EI - DI )
            !IF (ABS(VAL).GT.1.d0) THEN
            !   PRINT *, 'INDEX1(I-LK)',INDEX1(I-LK),Nt
            !   PRINT *, 'AA Ca,Pa',AA,Ca,Pa,DI,EI
            !   PRINT *, 'CDI,Cm,LK,I', CDI(I-LK),Cm(I-LK),LK,I
            !   PRINT *, 'TMP TMPOLD, Y(IK)',TMP, TMPOLD,Y(IK)
            !ENDIF
            !PRINT *,'I=',I, ' Y=',Y(IK),' val=',val
            LK   = 0
            IK   = IK + 1
            INFA = 0
            INFB = 0
         ELSE
            LK=LK+1
         END IF
      END DO
      RETURN
 200  VAL = 0.d0
      RETURN
      END FUNCTION MVNFUN
      END MODULE FUNCMOD1



!*****************************************************

     
      MODULE RIND
      IMPLICIT NONE
      PRIVATE
      PUBLIC :: RINDD, INITDATA, SETDATA,ECHO
 
    
      INTERFACE
         FUNCTION MVNFUN(N,Z) result (VAL)
         DOUBLE PRECISION,DIMENSION(:), INTENT(IN) :: Z
         INTEGER, INTENT(IN) :: N
         DOUBLE PRECISION :: VAL
         END FUNCTION MVNFUN
      END INTERFACE

    

      INTERFACE RINDD
      MODULE PROCEDURE RINDD
      END INTERFACE
      
      INTERFACE SETDATA
      MODULE PROCEDURE SETDATA
      END INTERFACE

      INTERFACE INITDATA
      MODULE PROCEDURE INITDATA
      END INTERFACE

      INTERFACE  BARRIER  
      MODULE PROCEDURE BARRIER
      END INTERFACE
      

      INTERFACE  CHCKLIM 
      MODULE PROCEDURE CHCKLIM
      END INTERFACE

      INTERFACE echo 
      MODULE PROCEDURE echo
      END INTERFACE
      
      INTERFACE  RCSWP
      MODULE PROCEDURE  RCSWP
      END INTERFACE

	INTERFACE RCSWAP
	MODULE PROCEDURE RCSWAP
	END INTERFACE

      INTERFACE swapRe
      MODULE PROCEDURE swapRe
      END INTERFACE

      INTERFACE swapint
      MODULE PROCEDURE swapint
      END INTERFACE

      INTERFACE CVSRTXC
      MODULE PROCEDURE CVSRTXC
      END INTERFACE

      INTERFACE CONDSORT2
      MODULE PROCEDURE CONDSORT2
      END INTERFACE

      INTERFACE COVSRT
      MODULE PROCEDURE COVSRT
      END INTERFACE

	INTERFACE COVSRT0
      MODULE PROCEDURE COVSRT0
      END INTERFACE
      
      INTERFACE PRINTCOF
      MODULE PROCEDURE PRINTCOF
      END INTERFACE 

     
                                !--------------------------------
      CONTAINS
      SUBROUTINE SETDATA (dREPS,dAEPS,dEPS2,dXc) 
      USE GLOBALDATA
      USE FIMOD
      IMPLICIT NONE
      DOUBLE PRECISION , INTENT(in) :: dREPS,dAEPS,dEPS2,dXc 
      RelEps   = dREPS          ! Constants controlling 
      AbsEps   = dAEPS          ! accuracy of integration 
      EPS2     = dEPS2          
      EPS      = SQRT(EPS2)

      xCutOff  = dXc 
      FxCutOff = FI(xCutOff) - FI( -xCutOff) 
      CFxCutOff= 1.d0-FxCutOff                              
      !XSMALL=ABSEPS
      XMAX=xCutOff
      IF (ABSEPS.LE.1e-4)      NsimMax=20000
      IF (ABSEPS.LE.1e-5)      NsimMax=40000
      IF (ABSEPS.LE.1e-6)      NsimMax=80000
	RETURN
      print *,'SETDATA: Requested parameters :'
      SELECT CASE (SCIS)
      CASE (:1)
         PRINT *,'SCIS = 1 SADAPT if NDIM<9 otherwise by KRBVRC'
      CASE (2) 
         PRINT *,'SCIS = 2 SADAPT'
      CASE (3)  
         PRINT *,'SCIS = 3 KRBVRC'
      CASE (4) 
         PRINT *,'SCIS = 4 KROBOV'
      CASE (5)  
         PRINT *,'SCIS = 5 RCRUDE'
      CASE (6) 
         PRINT *,'SCIS = 6 RLHCRUDE using MLHD center of cell' 
      CASE (7) 
         PRINT *,'SCIS = 7 RLHCRUDE using LHD and center of cell'
      CASE (8) 
         PRINT *,'SCIS = 8 RLHCRUDE using MLHD and random point'
      CASE (9:) 
         PRINT *,'SCIS = 9 RLHCRUDE using LHD and random point'
      END SELECT
      print *,'ABSEPS,RELEPS,EPS2,xCutOff,NsimMax',
     &     ABSEPS,RELEPS,EPS2,xCutOff,  NsimMax
      RETURN            
      END SUBROUTINE SETDATA
      
      SUBROUTINE INITDATA (speed)                                   
      USE GLOBALDATA
      USE FIMOD
      IMPLICIT NONE
      INTEGER , INTENT(in) :: speed
      INTEGER :: TMP
      TMP= max(11-speed,1)
      IF (.TRUE.) THEN
         NsimMax=10000
         SELECT case (speed)
         CASE (11:)
            ABSEPS = 1d-1 
         CASE (10)
            ABSEPS = 1d-2 
         CASE (7:9) 
            ABSEPS = 1d-2 
         CASE (4:6)
            NsimMax=20000
            ABSEPS = 1d-3
         CASE (:3)
            NsimMax=30000
            ABSEPS = 1d-4 
         END SELECT 
       
       
        
         RELEPS = MIN(ABSEPS ,1.d-2)
         TMP=MOD(TMP+1,3)+1
         EPS2=ABSEPS*(1D-1**TMP)
         !EPS2=MIN(EPS2,1D-6)
         !EPS2=1D-10
         xCutOff=ABS(FIINV(ABSEPS*1D-1*(1D-1**TMP)))
         xCutOff= MAX(xCutOff+0.5d0,4.d0)
      ELSE


         
         NsimMax=10000*TMP
         NsimMin=0
         ABSEPS=1.d-1**TMP
         
         
         IF (.FALSE.) THEN
            EPS2=1.d-10
            xCutOff=XMAX 
         ELSE
            
            xCutOff=ABS(FIINV(ABSEPS))
            xCutOff= MAX(xCutOff+0.5d0,4.d0)
                                !xCutOff= MIN(xCutOff,5.d0)
            EPS2 = ABSEPS*1.d-2
            ABSEPS = MIN(ABSEPS*(10**3),0.1D0)
         ENDIF
      ENDIF

      RELEPS = MIN(ABSEPS*1.d-1 ,1.d-2)
      EPS=SQRT(EPS2)
      
      FxCutOff = FI(xCutOff) - FI( -xCutOff)  
    
      CFxCutOff = 1.d0-FxCutOff
      ! Truncation parameters to FI and FIINV
      XSMALL = ABS(FIINV(.5d0+CFxCutOff/2.d0)) 
      XMAX   = xCutOff
	RETURN
	PRINT *, 'SP2THPDFALAN.exe CALLED'
      print *,'INITDATA: Requested parameters :'
      SELECT CASE (SCIS)
      CASE (:1)
         PRINT *,'SCIS = 1 SADAPT if NDIM<9 otherwise by KRBVRC'
      CASE (2) 
         PRINT *,'SCIS = 2 SADAPT'
      CASE (3)
         PRINT *,'SCIS = 3 KRBVRC'
      CASE (4)
         PRINT *,'SCIS = 4 KROBOV'
      CASE (5) 
         PRINT *,'SCIS = 5 RCRUDE'
      CASE (6) 
         PRINT *,'SCIS = 6 RLHCRUDE using MLHD center of cell' 
      CASE (7)
         PRINT *,'SCIS = 7 RLHCRUDE using LHD and center of cell'
      CASE (8)
         PRINT *,'SCIS = 8 RLHCRUDE using MLHD and random point'
      CASE (9:)
         PRINT *,'SCIS = 9 RLHCRUDE using LHD and random point'
      END SELECT
	PRINT *, 'ABSEPS = ', ABSEPS, ' RELEPS = ', RELEPS
	PRINT *, 'EPS2 = ', EPS2, ' xCutOff = ', xCutOff
	PRINT *, 'NsimMax = ', NsimMax
      RETURN                                                             
      END SUBROUTINE INITDATA

      SUBROUTINE ECHO(array)
      INTEGER ::j
      DOUBLE PRECISION,DIMENSION(:,:)::array
      DO j=1,size(array,1)
         PRINT 111,j,array(j,:)
111      FORMAT (i2,':',10F10.5)
      END DO
      END SUBROUTINE ECHO

      SUBROUTINE PRINTCOF(Ntd,A,B,INFI,COF,INDEX1)
      IMPLICIT NONE
      INTEGER ,INTENT(in) :: Ntd
      DOUBLE PRECISION, DIMENSION(:),INTENT(in) :: A,B,COF
      INTEGER , DIMENSION(:), INTENT(in) :: INFI,INDEX1
! Local variables
      INTEGER :: I, J, IJ
      PRINT '(/'' I     Limits'')'
      PRINT '(4X,''Lower  Upper  Lower Left of Cholesky Matrix'')'
      IJ = 0
      DO I = 1, Ntd
         J=MIN(I,10)
         IF ( INFI(I) .LT. 0 ) THEN 
            PRINT '(I2, '' -infin  infin '', 10F10.5)',
     &           INDEX1(I),  COF(IJ+1:IJ+J)
         ELSE IF ( INFI(I) .EQ. 0 ) THEN 
            PRINT '(I2, '' -infin'', F10.2, 1X, 10F10.5)',
     &           INDEX1(I), B(I), COF(IJ+1:IJ+J)
         ELSE IF ( INFI(I) .EQ. 1 ) THEN 
            PRINT '(I2, F10.2, ''  infin '', 10F10.5)',
     &           INDEX1(I), A(I), COF(IJ+1:IJ+J)
         ELSE 
            PRINT '(I2, 2F10.2, 1X, 10F10.5)', 
     &           INDEX1(I), A(I), B(I), COF(IJ+1:IJ+J)
         END IF
         IJ = IJ + I
      END DO
      END SUBROUTINE PRINTCOF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!******************* RINDD - the main program *********************!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE RINDD (fxind,Big1, Ex, xc1,Nt1,indI,Blo,Bup,INFIN)  
      USE GLOBALDATA 
      USE GLOBALCONST
	USE JACOBMOD
      USE FUNCMOD   ! use FIINV 
!	USE FUNCMOD1  ! use EXINV
      USE RCRUDEMOD
      USE KRBVRCMOD
      USE ADAPTMOD
      USE KROBOVMOD
      IMPLICIT NONE  
      DOUBLE PRECISION, DIMENSION(:  ), INTENT(out):: fxind 
      DOUBLE PRECISION, DIMENSION(:,:), INTENT(in) :: BIG1
      DOUBLE PRECISION, DIMENSION(:,:), INTENT(in) :: xc1 
      DOUBLE PRECISION, DIMENSION(:  ), INTENT(in) :: Ex            
      DOUBLE PRECISION, DIMENSION(:,:), INTENT(in) :: Blo, Bup  
      INTEGER,          DIMENSION(:  ), INTENT(in) :: indI,INFIN
      INTEGER, INTENT(in) :: Nt1 
! local variables
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: BIG
      INTEGER, DIMENSION(:),ALLOCATABLE :: XEDNI,INDEX2,IND
      INTEGER                          :: J,ix,Ntdcmj,Nst,Nsd,INFORM,
     &     Nc,Ntdc,Ntd,Nx,NDIM,NDIM1,Ndleft,MAXPTS,MINPTS, NLHD,I
      DOUBLE PRECISION         :: VALUE,SQ0,xx,fxc,quant,ERROR
      LOGICAL :: BCVSRT=.TRUE.    
      
! Find the size information
!~~~~~~~~~~~~~~~~~~~~~~~~~~
      Nt=Nt1
      Nc=size(xc1,dim=1)
      Nx=MAX(size(xc1,dim=2),1)
      Ntdc=size(BIG1,dim=1)
      IF (Nt+Nc.GT.Ntdc) Nt=Ntdc-Nc  ! make sure it does not exceed Ntdc-Nc
      Nd=Ntdc-Nt-Nc
      Ntd=Nt+Nd
      IF (Nd.LT.0) THEN
!         PRINT *,'RINDD Nt,Nd,Nc,Ntdc=',Nt,Nd,Nc,Ntdc
         STOP
      ENDIF
      
 !     PRINT *,'Nt Nd Nc Ntd Ntdc,',Nt, Nd, Nc, Ntd, Ntdc
     
! ALLOCATION
!~~~~~~~~~~~~
      IF (Nd.GT.0) THEN
         ALLOCATE(xd(Nd),CmXd(Nd),CDIXd(Nd))
         CmXd = 0.d0
         CDIXd = 0.D0
         xd = 0.d0 
      END IF
      ALLOCATE(BIG(Ntdc,Ntdc),Cm(Ntdc),Y(Ntd))
      ALLOCATE(index1(Ntdc),A(Ntd),B(Ntd),INFI(Ntd),xc(1:Nc)) 
      ALLOCATE(COF((Ntd+1)*(Ntd+2)/2),CDI(Ntd),xedni(Ntdc),index2(Ntdc))
 
! Initialization
!~~~~~~~~~~~~~~~

      BIG=BIG1
      index2=(/(J,J=1,Ntdc)/)
      xedni=index2
     
      fxind  = 0.d0             ! initialize
      CALL CVSRTXC(Nt,Nd,BIG,index2,Nst,Nsd,INFORM) ! sort index2 by decreasing cond. variance for Xc 
      IF (INFORM.GT.0) GOTO 110 ! degenerate case exit fxind=0 for all  
                                ! (should perhaps return NaN instead??)
      NDIM1 = Ntd-(Nsd-Nst-1)
      ALLOCATE(ind(1:(Nsd-Nst-1)))
      ind = index2(Nst+1:Nsd-1)
      IF (Nx.EQ.1) THEN
         BCVSRT=.FALSE.
      ELSE
         BCVSRT=.TRUE.
      ENDIF
      IF (BCVSRT) THEN          ! Conditionally sort all before hand       
         Cm  = Ex (1:Ntdc)
         xc = SUM(xc1(1:Nc,1:Nx),DIM=2)/DBLE(Nx)
         CALL BARRIER(xedni,Cm,xc,indI,Blo,Bup,INFIN,A,B,INFI) ! compute average integrationlimits   
         !print *,'infin',infin
         !print *,'xedni',xedni
         CALL COVSRT(Nt,Nd,BIG(1:Ntd,1:Ntd),Cm(1:Ntd),A,B,INFI,XEDNI, 
     &        INFIS,INFISD,NDIM,Y,COF,CDI)
                                !print *,'xedni',xedni
                                !print *,'index2',index2
         INDEX2(1:Ntd) = xedni(1:Ntd)
         INDEX1 = INDEX2
         !IF (INDEX1(Ntd-INFIS-INFISD).LE.Nt) NDIM = NDIM-1
         !PRINT *,'RINDD NDIM=', NDIM
      ENDIF
      CALL SORTBIG(Ntdc,BIG,index2,xedni) !sort BIG by index2
      !print *,'BIG'
      !CALL ECHO(BIG(1:Ntdc,1:Ntdc))
      !CALL PRINTCOF(Ntd,A,B,INFI,COF,INDEX1)
      !print *,'ind1',ind
      ind = xedni(ind)
      !print *,'ind2',ind
      
                                ! Now the loop over all different values of
                                ! variables Xc (the one one is conditioning on)
      DO  ix = 1, Nx            ! is started. The density f_{Xc}(xc(:,ix))                            
                                ! will be computed and denoted by  fxc.
         VALUE = 0.d0                                                       
         fxc   = 1.d0
         ERROR = 0.d0
                                ! Set the original means of the variables
         Cm = Ex (index2(1:Ntdc)) !   Cm(1:Ntdc)  =Ex (index1(1:Ntdc))
         
         xc = xc1(1:Nc,ix)
         DO J = 1, Nc           !Iterative conditioning on the last Nc variables  
            Ntdcmj=Ntdc-J
            SQ0=BIG(Ntdcmj+1,Ntdcmj+1) ! SQRT(Var(X(i)|X(i+1),X(i+2),...,X(Ntdc)))
                                ! i=Ntdc-J+1 (J=1 var(X(Ntdc))
            
            xx = (xc(index2(Ntdcmj+1)-Ntd)-Cm(Ntdcmj+1))/SQ0
            quant = 0.5D0 * xx * xx                                        
            IF (quant.GT.90.d0) GOTO 100                                                     
            
                                ! fxc probability density for i=Ntdc-J+1, 
                                ! fXc=f(X(i)|X(i+1),X(i+2)...X(Ntdc))*
                                !     f(X(i+1)|X(i+2)...X(Ntdc))*..*f(X(Ntdc))
            fxc = fxc*SQTWPI1*EXP(-quant)/SQ0   
                                ! conditional mean (expectation) 
                                ! E(X(1:i-1)|X(i),X(i+1),...,X(Ntdc)) 
            Cm(1:Ntdcmj) = Cm(1:Ntdcmj)+xx*BIG (1:Ntdcmj,Ntdcmj+1)
         ENDDO  
                                !print *,'density',fxc                 ! J
                                !PRINT *, 'Rindd, Cm=',Cm(xedni(max(1,Nt-5):Ntdc))
                                !PRINT *, 'Rindd, Cm=',Cm(xedni(1:Ntdc))

         IF (fxc .LT.fxcEpss) GOTO 100 ! Small probability don't bother calculating it

                                !set the integration limits A,B and INFI
                                !NOTE: A and B are integration limits with Cm subtracted
         CALL BARRIER(xedni,Cm,xc,indI,Blo,Bup,INFIN,A,B,INFI)
         CALL CHCKLIM(ind,A,B,INFI,EPS,INFORM)
         IF (INFORM.GT.0) GO TO 100   !degenerate case mean of deterministic variables  
                                      ! outside the barriers     
         
      
         IF (NDIM1.LT.1) THEN !degenerate case, No relevant variables left to integrate  
!            Print *,'rindd ndim1',Ndim1
            IF (Nd.GT.0) THEN ! replace xd with the mean
               Ndleft=Nd
               I=Ntd
               DO WHILE (Ndleft.GT.0)
                  IF (index2(I).GT.Nt) THEN ! isXd
                     xd (Ndleft) =  Cm (I)
                     Ndleft=Ndleft-1                  
                  END IF
                  I=I-1
               ENDDO
               VALUE = jacob (xd,xc) ! jacobian of xd,xc
            ELSE
               VALUE = 1.d0
            END IF
            ERROR = MAX(ABSEPS,RELEPS*VALUE)
            GOTO 100
         ENDIF
        
         IF (BCVSRT) THEN       ! conditionally sorted before hand just scale integration limits
            WHERE (ABS(CDI).GT.0.d0) ! EPS2)
               A = A/CDI
               B = B/CDI
            END WHERE
                                !The following could be done smarter by changing Blo Bup and INFIN
            DO I = 1,Ntd
               IF ( CDI(I) .LT.0.d0.AND.INFI(I).GE.0) THEN
                  CALL SWAPRE( A(I), B(I) ) 
                  IF ( INFI(I) .NE. 2 ) INFI(I) = 1 - INFI(I)
               END IF
            ENDDO
         ELSE
            INDEX1 = index2   
            CALL COVSRT(Nt,Nd,BIG(1:Ntd,1:Ntd),Cm(1:Ntd),A,B,INFI, 
     &           INDEX1,INFIS,INFISD,NDIM,Y,COF,CDI)
            
         ENDIF
         CALL INITFUN(VALUE,INFORM)
         IF (INFORM.GT.0) GO TO 100 !degenerate case: the integral can be calculated excactly, ie. 
                                !mean of deterministic variables outside the barriers, or NDIM = 1     
         
         ! CALL PRINTCOF(Ntd,A,B,INFI,COF,INDEX1)
         MAXPTS=NSIMmax
         MINPTS=NSIMmin
      
         IF ( NDIM .GT. 19. AND. SCIS.EQ.2) THEN
 !           print *, 'Ndim too large for SADMVN => Calling KRBVRC'
            SCIS=3
         ENDIF
         
         !print * ,'RINDD: CDI',MINVAL(ABS(CDI(1:Nt+Nd-INFIS-INFISD)))
         !print * ,'RINDD: Ndim',Ndim
         SELECT CASE (SCIS)
         CASE (:1)
            IF (NDIM.lt.9) THEN
               CALL SADAPT(Ndim,MAXPTS,MVNFUN,ABSEPS,
     &              RELEPS,ERROR,VALUE,INFORM)
            ELSE
               CALL KRBVRC( NDIM, MINPTS, MAXPTS, MVNFUN, ABSEPS,RELEPS,
     &              ERROR, VALUE, INFORM )
            ENDIF
         CASE (2)            !        Call the subregion adaptive integration subroutine
            CALL SADAPT(Ndim,MAXPTS,MVNFUN,ABSEPS,
     &           RELEPS,ERROR,VALUE,INFORM)
            VALUE=MAX(VALUE,0.d0)
         CASE (3)             !       Call the Lattice rule integration procedure
            CALL KRBVRC( NDIM, MINPTS, MAXPTS, MVNFUN, ABSEPS,
     &        RELEPS, ERROR, VALUE, INFORM )
         CASE (4)             !       Call the Lattice rule integration procedure
           CALL KROBOV( NDIM, MINPTS, MAXPTS, MVNFUN, ABSEPS, RELEPS,
     &          ERROR, VALUE, INFORM )
         CASE (5)             ! Call Crude Monte Carlo integration procedure
           CALL RANMC( NDIM, MAXPTS, MVNFUN, ABSEPS, 
     &          RELEPS, ERROR, VALUE, INFORM )
         CASE (6:)          ! Call Crude Latin Hypercube Monte Carlo integration procedure
            Nlhd=max(1,5+rateLHD*Ndim) ! size of LHD
            IF (SCIS.GE.8) THEN
               useMIDP=.FALSE.  ! Do not Use midpoint of Cell, but sample randomly
            ELSE
               useMIDP=.TRUE.   ! Use center of Cell
            END IF
            IF (SCIS.EQ.7.OR.SCIS.EQ.9) THEN
               MLHD=.FALSE.     ! Do not modify LHD
            ELSE
               MLHD=.TRUE.      ! Modify LHD
            ENDIF
            if (rateLHD.lt.1) then ! make sure
               Nlhd=1 
               usemidp=.FALSE.
!               print * ,'Rindscis: only able to useMIDP if rateLHD>0'
            endif 
            CALL  RANLHMC( NDIM, MAXPTS,NLHD,USEMIDP,MLHD,MVNFUN,ABSEPS, 
     &           RELEPS, ERROR, VALUE, INFORM )
         END SELECT   

         if (allocated(COV)) then ! save the coefficient of variation in COV
            if ((VALUE.gt.0.d0))  COV(ix)=ERROR/VALUE/3.0d0
         endif
 !        IF (INFORM.gt.0) print *,'RINDD, INFORM,error =',inform,error
 100     fxind(ix) = VALUE*fxc                                          
         !PRINT *,'RINDD MINPTS=',MINPTS
      ENDDO                     !ix

 110  CONTINUE
      IF (ALLOCATED(xc))     DEALLOCATE(xc) 
      IF (ALLOCATED(xd))     DEALLOCATE(xd) 
      IF (ALLOCATED(Cm))     DEALLOCATE(Cm) 
      IF (ALLOCATED(BIG))    DEALLOCATE(BIG) 
      IF (ALLOCATED(index2)) DEALLOCATE(index2)
      IF (ALLOCATED(index1)) DEALLOCATE(index1)
      IF (ALLOCATED(xedni))  DEALLOCATE(xedni)
      IF (ALLOCATED(A))      DEALLOCATE(A)
      IF (ALLOCATED(B))      DEALLOCATE(B)  
      IF (ALLOCATED(Y))      DEALLOCATE(Y) 
      IF (ALLOCATED(COF))    DEALLOCATE(COF)
      IF (ALLOCATED(CDI))    DEALLOCATE(CDI)  
      IF (ALLOCATED(CDIXd))  DEALLOCATE(CDIXd)  
      IF (ALLOCATED(CmXd))   DEALLOCATE(CmXd)  
      IF (ALLOCATED(INFI))   DEALLOCATE(INFI)
      RETURN                                                            
      END SUBROUTINE RINDD
      
   
      SUBROUTINE CHCKLIM(ind,Hlo,Hup,INFIN,SQEPS,INFORM)
      IMPLICIT NONE
      INTEGER,DIMENSION(:),           INTENT(in)  :: ind 
      DOUBLE PRECISION, DIMENSION(:), INTENT(in)  :: Hlo,Hup
      INTEGER,          DIMENSION(:), INTENT(in)  :: INFIN
      DOUBLE PRECISION,               INTENT(in)  :: SQEPS
      INTEGER,                        INTENT(out) :: INFORM
! LOCAL variable
      INTEGER :: J,K,N
! CHCKLIM Check if the conditional mean Cm = E(Xt,Xd|Xc) for the deterministic variables
!         is between the barriers, i.e., Hlo=A-Cm< 0 <Hup=B-Cm
!     INFIN  INTEGER, array of integration limits flags:
!            if INFIN(I) < 0, Ith limits are (-infinity, infinity);
!            if INFIN(I) = 0, Ith limits are (-infinity, Hup(I)];
!            if INFIN(I) = 1, Ith limits are [Hlo(I), infinity);
!            if INFIN(I) = 2, Ith limits are [Hlo(I), Hup(I)].
      INFORM=1
      N = SIZE(ind)
      DO J=1,N
         K = ind(J)
         IF (INFIN(K).GE.0) THEN
            IF ((INFIN(K).NE.0).AND.(SQEPS .LT.Hlo (K))) RETURN
            IF ((INFIN(K).NE.1).AND.(-SQEPS .GT. Hup(K))) RETURN
         ENDIF
      ENDDO
      INFORM=0
      RETURN
      END SUBROUTINE CHCKLIM
     

      SUBROUTINE BARRIER(xedni,Cm,xc,indI,Blo,Bup,INFIN,Hlo,Hup,INFI)     
      USE GLOBALDATA, ONLY : xCutOff
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(:  ), INTENT(in)  :: xc,Cm 
      INTEGER,          DIMENSION(:  ), INTENT(in)  :: xedni,indI,INFIN    
      DOUBLE PRECISION, DIMENSION(:,:), INTENT(in)  :: Blo,Bup   
      DOUBLE PRECISION, DIMENSION(:  ), INTENT(out) :: Hlo,Hup
      INTEGER,          DIMENSION(:  ), INTENT(out) :: INFI  
!Local variables
      INTEGER   :: I, J, K, L,Mb,Nb,NI,Nc           
!this procedure set Hlo,Hup and INFI according to Blo/Bup and INFIN
!
!     INFIN  INTEGER, array of integration limits flags:
!            if INFIN(I) < 0, Ith limits are (-infinity, infinity);
!            if INFIN(I) = 0, Ith limits are (-infinity, Hup(I)];
!            if INFIN(I) = 1, Ith limits are [Hlo(I), infinity);
!            if INFIN(I) = 2, Ith limits are [Hlo(I), Hup(I)].
 
      Mb=size(Blo,DIM=1)
      Nb=size(Blo,DIM=2)
      NI=size(indI,DIM=1)
      Nc=size(xc,DIM=1)
!      IF (Mb.GT.Nc+1) print *,'barrier: Mb,Nc =',Mb,Nc
!      IF (Nb.NE.NI-1) print *,'barrier: Nb,NI =',Nb,NI
      DO J = 2, NI 
         DO I =indI (J - 1) + 1 , indI (J)
            L=xedni(I)  
            INFI(L)=INFIN(J-1)
            Hlo (L) = -xCutOff
            Hup (L) =  xCutOff
            IF (INFI(L).GE.0) THEN
               IF  (INFI(L).NE.0) THEN
                  Hlo (L) = Blo (1, J - 1)-Cm(L)  
                  DO K = 1, Mb-1                                                  
                     Hlo(L) = Hlo(L)+Blo(K+1,J-1)*xc(K)              
                  ENDDO         ! K
               ENDIF
               IF  (INFI(L).NE.1) THEN
                  Hup (L) = Bup (1, J - 1)-Cm(L)  
                  DO K = 1, Mb-1                                                  
                     Hup(L) = Hup(L)+Bup(K+1,J-1)*xc(K)
                  ENDDO  
               ENDIF                 !
            ENDIF            
         ENDDO ! I
      ENDDO ! J
!      print * ,'barrier hup:',size(Hup),Hup(xedni(1:indI(NI)))
!      print * ,'barrier hlo:',size(Hlo),Hlo(xedni(1:indI(NI)))
      RETURN                                                             
      END SUBROUTINE BARRIER                                              
      

            !********************************************************************

      SUBROUTINE CVSRTXC (Nt,Nd,R,index1,Nst,Nsd,INFORM)
      USE GLOBALDATA, ONLY :  XCEPS2
      IMPLICIT NONE
      INTEGER, INTENT(in) :: Nt,Nd 
      DOUBLE PRECISION, DIMENSION(:,:), INTENT(inout) :: R
      INTEGER,          DIMENSION(:  ), INTENT(inout) :: index1 
      INTEGER, INTENT(out) :: Nst,Nsd,INFORM 
! local variables
      DOUBLE PRECISION, DIMENSION(:  ), allocatable   :: SQ
      INTEGER,          DIMENSION(1  )                :: m
      INTEGER :: m1,r1,ix,Ntdc,Ntd,Nc
! R         = Input: Cov(X) where X=[Xt Xd Xc] is stochastic vector
!            Output: sorted Conditional Covar. matrix   Shape N X N  (N=Ntdc=Nt+Nd+Nc) 
! index1    = indices to the variables original place.   Size  Ntdc
! xedni     = indices to the variables new place.        Size  Ntdc
! Nst       = index to the last stochastic variable
!            among Nt first of Xt after conditioning on Xc.                                 
! Nsd       = index to the first stochastic variable
!            among Xd  after conditioning on Xc.                                 
! 
! R=Cov([Xt,Xd,Xc]) is a covariance matrix of the stochastic vector X=[Xt Xd Xc]
! where the variables Xt, Xd and Xc have the size Nt, Nd and Nc, respectively. 
! Xc is (are) the conditional variable(s).
! Xd and Xt are the variables to integrate. (Xd,Xt = variables in the jacobian and indicator respectively)
      INFORM=0
      Ntdc = size(R,DIM=1)
      Ntd  = Nt+Nd
      Nc   = Ntdc-Ntd
      Nsd  = Nt+1;Nst=Nt
      IF (Nc.LT.1) RETURN

      ALLOCATE(SQ(1:Ntdc))
      DO ix = 1, Ntdc
         R(ix,ix) = R(ix,ix)+1.D-12  ! adding a nugget effect to ensure that CHOLESKY factorization is not corrupted by numerical errors
         SQ(ix)   = R(ix,ix)
      ENDDO
      DO ix = 1, Nc             ! Condsort Xc
         r1=Ntdc-ix
         m=r1+2-MAXLOC(SQ(r1+1:Ntd+1:-1)) 
         IF (SQ(m(1)).LE.XCEPS2) THEN
!            PRINT *,'CVSRTXC: Degenerate case of Xc(Nc-J+1) for J=',ix 
            INFORM=1
            GOTO 200            ! RETURN    !degenerate case
         ENDIF
         m1=index1(m(1))
         CALL swapint(index1(m(1)),index1(r1+1))
         SQ(r1+1) = SQRT(SQ(m(1)))
         R(m1,m1) = SQ(r1+1)
                                ! sort and calculate conditional covariances
         CALL CONDSORT2(Ntd,R,SQ,index1,Nst,Nsd,m1,r1)
      ENDDO
!      r1=0
!      DO ix=Ntd,Nsd,-1
!         r1=r1+1
!         IF (r1.GE.Nsd) GOTO 100 
!         CALL SWAPINT(INDEX1(r1),INDEX1(ix))
!         CALL SWAPRE(SQ(r1),SQ(ix))
!      END DO
                    ! ix 
  
 200  DEALLOCATE(SQ)
      RETURN                                                             
      END SUBROUTINE CVSRTXC
      


      SUBROUTINE CONDSORT2(Ntd,R,SQ,index1,Nstoc,NstoXd,m1,N)
      USE GLOBALDATA, ONLY : EPS2,XCEPS2
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(:,:), INTENT(inout) :: R
      DOUBLE PRECISION, DIMENSION(:),   INTENT(inout) :: SQ 
      INTEGER,          DIMENSION(:  ), INTENT(inout)   :: index1 
      INTEGER, INTENT(inout)   :: Nstoc,NstoXd
      INTEGER, INTENT(in)   :: m1,N,Ntd
! local variables
      DOUBLE PRECISION :: EPSL
      INTEGER       :: Nsold,Ndold
      INTEGER       :: r1,c1,row,col,iy
! save their old values
      Nsold=Nstoc;Ndold=NstoXd
      !EPSL=EPS2
      !IF (Ntd.LT.10) EPSL = 1D-10
	EPSL = XCEPS2

                                ! Calculating conditional variances for the 
                                ! Xc variables. 
      DO row=Ntd+1,N
         r1=index1(row)
         R(r1,m1)=R(r1,m1)/R(m1,m1)
         R(m1,r1)=R(r1,m1)
         SQ(row)=R(r1,r1)-R(r1,m1)*R(m1,r1)      !/R(m1,m1)
         IF (SQ(row).LE.XCEPS2) THEN 
            R(r1,r1)=0.d0
            SQ(row)=0.d0
            RETURN              ! degenerate case XIND should return NaN
         ELSE
            R(r1,r1)=SQ(row)
            DO col=row+1,N
               c1=index1(col)
               R(c1,m1)=R(c1,m1)/R(m1,m1)
               R(m1,c1)=R(c1,m1)
               R(c1,r1)=R(r1,c1)-R(r1,m1)*R(m1,c1)      !/R(m1,m1)
               R(r1,c1)=R(c1,r1)
            ENDDO
         ENDIF
      ENDDO                     ! Calculating conditional variances for the 
                                ! first Nstoc variables.
                                ! variables with variance less than EPS2 
                                ! will be treated as deterministic and not 
                                ! stochastic variables and are therefore moved
                                ! to the end among these Nt-Nj first variables.
                                ! Nstoc is the # of variables we treat 
                                ! stochastically 
      iy=1
      DO WHILE (iy.LE.Nstoc)
         r1=index1(iy)
         R(r1,m1)=R(r1,m1)/R(m1,m1)
         R(m1,r1)=R(r1,m1)
         SQ(iy)=R(r1,r1)-R(r1,m1)*R(m1,r1)              !/R(m1,m1)
         IF (SQ(iy).LE.EPSL) THEN
            IF (iy.LT.Nstoc) THEN               
               r1=index1(Nstoc)
               R(r1,m1)=R(r1,m1)/R(m1,m1)
               R(m1,r1)=R(r1,m1)
               SQ(Nstoc)=R(r1,r1)-R(r1,m1)*R(m1,r1) !/R(m1,m1)
               DO WHILE ((SQ(Nstoc).LE.EPSL).AND.(iy.LT.Nstoc))
                  SQ(Nstoc)=0.d0 !MAX(0.d0,SQ(Nstoc))
                  Nstoc=Nstoc-1
                  r1=index1(Nstoc)
                  SQ(Nstoc)=R(r1,r1)-R(r1,m1)*R(m1,r1) !/R(m1,m1)
               END DO 
               CALL swapint(index1(iy),index1(Nstoc)) ! swap indices
               !CALL swapre(SQ(iy),SQ(Nstoc))          ! swap values
               SQ(iy)=SQ(Nstoc);
            ENDIF
            SQ(Nstoc)=0.d0
            Nstoc=Nstoc-1
         ENDIF
         iy=iy+1
      END DO  
                                ! Calculating conditional variances for the 
                                ! stochastic variables Xd. 
                                ! Variables with conditional variance less than
                                ! EPS2 are moved to the beginning among these.
      DO iy=Ndold,MIN(Ntd,N)
         r1=index1(iy)
         R(r1,m1)=R(r1,m1)/R(m1,m1)
         R(m1,r1)=R(r1,m1)
         SQ(iy)=R(r1,r1)-R(r1,m1)*R(m1,r1)              !/R(m1,m1)
         IF (SQ(iy).LE.EPSL) THEN
            CALL swapint(index1(iy),index1(NstoXd))
            !CALL swapre(SQ(iy),SQ(NstoXd)) !
            SQ(iy)=SQ(NstoXd);SQ(NstoXd)=0.d0
            NstoXd=NstoXd+1
         ENDIF                  ! SQ < EPS2
      ENDDO
      
            
            ! Calculating Covariances for non-deterministic variables
      DO row=1,Nstoc
         r1=index1(row)
         R(r1,r1)=SQ(row)
         DO col=row+1,Nstoc
            c1=index1(col)
            R(c1,r1)=R(r1,c1)-R(r1,m1)*R(m1,c1)         !/R(m1,m1)
            R(r1,c1)=R(c1,r1)
         ENDDO
         DO col=NstoXd,N
            c1=index1(col)
            R(c1,r1)=R(r1,c1)-R(r1,m1)*R(m1,c1)          !/R(m1,m1)
            R(r1,c1)=R(c1,r1)
         ENDDO
      ENDDO            
      DO row=NstoXd,MIN(Ntd,N)
         r1=index1(row)
         R(r1,r1)=SQ(row)         
         DO col=row+1,N
            c1=index1(col)
            R(c1,r1)=R(r1,c1)-R(r1,m1)*R(m1,c1)          !/R(m1,m1)
            R(r1,c1)=R(c1,r1)
         ENDDO
      ENDDO
                                ! Set covariances for Deterministic variables to zero
                                ! in order to avoid numerical problems
      DO row=Ndold,NStoXd-1
         r1=index1(row)
         SQ(row)  = 0.d0 !MAX(SQ(row),0.d0)
         R(r1,r1) = SQ(row)
         DO col=row+1,N
            c1=index1(col)
            R(c1,r1)=0.d0
            R(r1,c1)=0.d0
         ENDDO
         DO col=1,Nsold
            c1=index1(col)
            R(c1,r1)=0.d0
            R(r1,c1)=0.d0
         ENDDO
      ENDDO
      
      DO row=Nstoc+1,Nsold
         r1=index1(row)
         SQ(row)  = 0.d0 !MAX(SQ(row),0.d0)
         R(r1,r1)=SQ(row)
         DO col=1,row-1
            c1=index1(col)
            R(c1,r1)=0.d0
            R(r1,c1)=0.d0
         ENDDO
         DO col=NstoXd,N
            c1=index1(col)
            R(c1,r1)=0.d0
            R(r1,c1)=0.d0
         ENDDO
      ENDDO
      RETURN   
      END SUBROUTINE CONDSORT2

      SUBROUTINE SORTBIG(Ntdc,R,index1,xedni)
      IMPLICIT NONE
      INTEGER, INTENT(in) :: Ntdc 
      DOUBLE PRECISION, DIMENSION(:,:), INTENT(inout) :: R
      INTEGER,          DIMENSION(:  ), INTENT(inout) :: index1 
      INTEGER,          DIMENSION(:  ), INTENT(out  ) :: xedni 
! local variables
      DOUBLE PRECISION, DIMENSION(Ntdc)  :: SQ
      INTEGER :: r1,c1,r2,c2
      LOGICAL :: changed
! R         = Input: Cov(X) where X=[Xt Xd Xc] is stochastic vector
!            Output: sorted Conditional Covar. matrix   Shape N X N  (N=Ntdc=Nt+Nd+Nc) 
! index1    = indices to the variables original place.   Size  Ntdc
! xedni     = indices to the variables new place.        Size  Ntdc
! Nst       = index to the last stochastic variable
!            among Nt first of Xt after conditioning on Xc.                                 
! Nsd       = index to the first stochastic variable
!            among Xd  after conditioning on Xc.                                 
! 
! R=Cov([Xt,Xd,Xc]) is a covariance matrix of the stochastic vector X=[Xt Xd Xc]
! where the variables Xt, Xd and Xc have the size Nt, Nd and Nc, respectively. 
! Xc is (are) the conditional variable(s).
! Xd and Xt are the variables to integrate. (Xd,Xt = variables in the jacobian and indicator respectively)
      changed=.FALSE.
      DO r2=1,Ntdc
         SQ(r2) = R(r2,r2)
      END DO
      DO r2=Ntdc,1,-1       ! sorting the upper triangular of the 
         r1=index1(r2)          ! covariance matrix according to index1
         xedni(r1)=r2
         !PRINT *,'condsort,xedni',xedni
         !PRINT *,'condsort,r1,r2',r1,r2
         IF ((r1.NE.r2).OR.(changed)) THEN
            changed = .TRUE.
            R(r2,r2) = SQ(r1)
            DO c2=r2+1,Ntdc
               c1=index1(c2)
               IF (c1.GT.r1) THEN
                  R(r2,c2)=R(c1,r1)
               ELSE
                  R(r2,c2)=R(r1,c1)
               END IF
            END DO
         END IF              
      END DO
                                ! you may sort the lower triangular according  
                                ! to index1 also, but it is not needed
                                ! since R is symmetric.  Uncomment the
                                ! following if the whole matrix is needed
!      IF (changed.EQ.1) THEN
!         DO c2=1,Ntdc
!            DO r2=c2+1,Ntdc
!               R(r2,c2)=R(c2,r2) ! R symmetric
!            END DO
!         END DO
!      ENDIF
      END SUBROUTINE SORTBIG

      SUBROUTINE COVSRT0( Nt,Nd,R,Cm,A,B,INFI,INDEX1, 
     &              INFIS,INFISD, NDIM, Y, COF,CDI )
      USE FIMOD
      USE GLOBALDATA, ONLY : EPS2
      IMPLICIT NONE
*
*     Subroutine to sort integration limits and determine Cholesky factor.
*
      INTEGER,                          INTENT(in)    :: Nt,Nd
      DOUBLE PRECISION, DIMENSION(:,:), INTENT(in)    :: R
      DOUBLE PRECISION, DIMENSION(:  ), INTENT(inout) :: Cm,A,B
      INTEGER,          DIMENSION(:  ), INTENT(inout) :: INFI
      INTEGER,          DIMENSION(:  ), INTENT(inout) :: INDEX1
      INTEGER,                          INTENT(out) :: INFIS,INFISD,NDIM
      DOUBLE PRECISION, DIMENSION(:  ), INTENT(out)   :: Y, COF, CDI
! Local variables
      INTEGER :: N,I, J, K, L, M, II, IJ, IL, JMIN,Ndleft
      DOUBLE PRECISION :: SUMSQ, AJ, BJ, TMP, SQTWPI, D, E,EPSL 
      DOUBLE PRECISION :: CVDIAG, AMIN, BMIN, DMIN, EMIN, YL, YU
      PARAMETER ( SQTWPI = 2.506628274631001D0) !, EPSL = 1D-10 )
      
      EPSL=EPS2
      
     
      IJ     = 0
      INFIS  = 0
      INFISD = 0
      Ndim   = 0
      Ndleft = Nd
      N      = Nt+Nd
      IF (N.LT.10) EPSL = 1D-10

      DO I = 1, N
         IF (R(I,I).LT.EPSL) INFI(I)=-1
         IF ( INFI(I) .LT. 0.) THEN
            IF (INDEX1(I).LE.Nt)  THEN 
               INFIS = INFIS+1
            ELSEIF (R(I,I).LT.EPSL) THEN
               INFISD = INFISD+1
            ENDIF
         ENDIF
        
         DO J = 1, I
            IJ = IJ + 1
            COF(IJ) = R(J,I)
         END DO
      END DO
      COF(IJ+1:SIZE(COF))=0.d0
      CDI(N-INFIS-INFISD+1:SIZE(CDI))=0.d0
      !PRINT *,'COVSRT'
      !CALL PRINTCOF(N,A,B,INFI,COF,INDEX1)

      IF ( INFIS +INFISD .GE. N ) RETURN
*     
*     First move Xd to the beginning 
!         IJ = MIN(Nd,Nt)
!         DO I = N, N-IJ+1, -1
!            DO J = 1,IJ
!               CALL RCSWP( J, I, A, B, INFI, N, COF )
!               CALL SWAPINT(INDEX1(J),INDEX1(I))
!               CALL SWAPRE(Cm(J),Cm(I))
!            END DO
!         END DO
*
*     Move any doubly infinite limits of Xt to innermost positions.
*
      DO I = N, N-INFIS+1, -1
         IF ( INFI(I) .GE. 0 .OR. INDEX1(I).GT.Nt) THEN 
            DO J = 1,I-1
               IF ( INFI(J) .LT. 0 .AND. INDEX1(J).LE.Nt) THEN
                  CALL RCSWP( J, I, A, B, INFI, N, COF )
                  CALL SWAPINT(INDEX1(J),INDEX1(I))
                  CALL SWAPRE(Cm(J),Cm(I))
                  GO TO 10
               ENDIF
            END DO
         ENDIF
 10   END DO
*     Then any redundant variables of Xd to the next innermost positions. 
         DO I = N-INFIS, N-INFIS-INFISD+1, -1
            K=INDEX1(I)
            IF ( R(K,K) .GE. EPSL) THEN 
               DO J = 1,I-1
                  M=INDEX1(J)
                  IF ( R(M,M) .LT. EPSL) THEN
                     CALL RCSWP( J, I, A, B, INFI, N, COF )
                     CALL SWAPINT(INDEX1(J),INDEX1(I))
                     CALL SWAPRE(Cm(J),Cm(I))
                     GO TO 15
                  ENDIF
               END DO
            ENDIF
 15      END DO

*     Then any doubly infinite limits of Xd to the next innermost positions. 
!         DO I = N-INFIS, N-INFIS-INFISD+1, -1
!            IF ( INFI(I) .GE. 0) THEN 
!               DO J = 1,I-1
!                  IF ( INFI(J) .LT. 0) THEN
!                     CALL RCSWP( J, I, A, B, INFI, N, COF )
!                     CALL SWAPINT(INDEX1(J),INDEX1(I))
!                     CALL SWAPRE(Cm(J),Cm(I))
!                     GO TO 15
!                  ENDIF
!               END DO
!            ENDIF
! 15      END DO
         
*
*     Sort remaining limits and determine Cholesky factor.
*
      II = 0
      DO I = 1, N-INFIS-INFISD
*
*        Determine the integration limits for variable with minimum
*        expected probability and interchange that variable with Ith.
*
         DMIN = 0.d0
         EMIN = 1.d0
         JMIN = I
         CVDIAG = 0.d0
         IJ = II
         DO J = I,N-INFIS-INFISD
            IF ( COF(IJ+J) .GT. EPSL ) THEN
               SUMSQ = SQRT( COF(IJ+J) )
               TMP = 0.d0  ! =  the conditional mean of Y(I) given Y(1:I-1)
               DO K = 1, I-1
                  TMP = TMP + COF(IJ+K)*Y(K)
               END DO
               IF (INFI(J).LT.0) GO TO 30  ! May have infinite int. limits if Nd>0
               IF (INFI(J).NE.0) AJ = ( A(J) - TMP )/SUMSQ
               IF (INFI(J).NE.1) BJ = ( B(J) - TMP )/SUMSQ
 30            CALL MVNLMS( AJ, BJ, INFI(J), D, E )
               IF ( EMIN + D .GE. E + DMIN ) THEN
                  JMIN = J
                  AMIN = AJ
                  BMIN = BJ
                  DMIN = D
                  EMIN = E
                  CVDIAG = SUMSQ
               ENDIF
            ENDIF
            IJ = IJ + J 
                                !               IF (J.GE.Nd.AND.Ndleft.GT.0) THEN
                                !                   Ndleft=Ndleft-1
                                !                   GOTO 40
                                !                ENDIF
         END DO
 40      IF ( JMIN .GT. I ) THEN
            CALL RCSWP( I, JMIN, A,B, INFI, N, COF )
            CALL SWAPINT(INDEX1(I),INDEX1(JMIN))
            CALL SWAPRE(Cm(I),Cm(JMIN))
         END IF
         IL = II + I
         COF(IL) = CVDIAG
*     
*     Compute Ith column of Cholesky factor.
*     Compute expected value for Ith integration variable (without considering the jacobian) and
*     scale Ith covariance matrix row and limits.
*
         IF ( CVDIAG .GT. 1.d-16 ) THEN
            CDI(I) = CVDIAG         ! Store the diagonal element
            Ndim   = Ndim+1         ! counter for Number of relevant dimensions to integrate
            DO L = I+1, N-INFIS-INFISD               
               COF(IL+I) = COF(IL+I)/CVDIAG
               IJ = II + I
               COF(IL+L) = COF(IL+L) - COF(IL+I)*COF(IL+I)
               IF (.TRUE..OR.COF(IL+L).GT.EPSL) THEN  ! Calculate conditional covariances Y(I+1:N) given Y(1:I)
                  DO J = I+1, L-1
                     COF(IL+J) = COF(IL+J) - COF(IL+I)*COF(IJ+I)
                     IJ = IJ + J            ! update column index
                  END DO
               ELSE  ! NEW CALL This is WRONG!!!!!!!!!!!!!!!!!
                  ! Must move the deterministic variables to the end also
                  ! set covariances with the rest to zero to avoid numerical difficulties
                  ! This is important when EPSL is set large to increase speed
                  ! or high accuracy is requested 
                  ! Set row to zero
                  !PRINT *, ' COF(IL+L)=',COF(IL+L)
                  DO J = I+1, L
                     COF(IL+J) = 0.d0
                  END DO
                  ! Set  column to zero
                  IJ = L
                  DO J = L+1, N-INFIS-INFISD
                     IJ = IJ + J -1
                     COF(IL+IJ) = 0.d0
                  END DO
                  !PRINT *,'setting zeros'
                  !CALL PRINTCOF(N,A,B,INFI,COF,INDEX1)
               END IF
               IL = IL + L            ! goto next row i.e. update row index
            END DO
            

            IF ( EMIN .GT. DMIN + EPSL ) THEN
               YL = 0.d0
               YU = 0.d0
               IF (INFI(I).GE.0) THEN
                  IF ( INFI(I) .NE. 0 ) YL =-EXP(-0.5d0*AMIN**2)/SQTWPI
                  IF ( INFI(I) .NE. 1 ) YU =-EXP(-0.5d0*BMIN**2)/SQTWPI
               ENDIF
               Y(I) = ( YU - YL )/( EMIN - DMIN )
            ELSE
               SELECT CASE (INFI(I))
               CASE (:-1) 
                  Y(I) = 0.d0
               CASE (0) 
                  Y(I) = BMIN
               CASE (1)
                  Y(I) = AMIN
               CASE (2:)
                  Y(I) = ( AMIN + BMIN )*0.5d0
               END SELECT 
            END IF
            DO J = 1, I
               II = II + 1
               COF(II) = COF(II)/CVDIAG
            END DO
            A(I) = A(I)/CVDIAG
            B(I) = B(I)/CVDIAG
         ELSE
            DO L = I+1, N-INFIS                
               COF(IL+I) = 0.d0
               IL = IL + L
            END DO
*     
*     If the covariance matrix diagonal entry is zero, 
*     permute limits and/or rows, if necessary.
*
*               
            DO J = I-1, 1, -1
               IF ( ABS( COF(II+J) ) .GT. EPSL ) THEN
                  CDI(I) = COF(II+J)
                  A(I)   = A(I)/COF(II+J)
                  B(I)   = B(I)/COF(II+J)
                  IF ( COF(II+J).LT.0.d0.AND.INFI(I).GE.0) THEN
                     CALL SWAPRE( A(I), B(I) ) 
                     IF ( INFI(I) .NE. 2 ) INFI(I) = 1 - INFI(I)
                  END IF
                  DO L = 1, J
                     COF(II+L) = COF(II+L)/COF(II+J)
                  END DO ! L
                  DO L = J+1, I-1 
                     IF( COF((L-1)*L/2+J+1) .GT. 0.d0 ) THEN
                        IJ = II
                        DO K = I-1, L, -1 
                           DO M = 1, K
                              CALL SWAPRE( COF(IJ-K+M), COF(IJ+M) )
                           END DO
                           CALL SWAPRE( A(K), A(K+1) ) 
                           CALL SWAPRE( B(K), B(K+1) )
                           CALL SWAPRE( Cm(K),Cm(K+1))
                           CALL SWAPRE( CDI(K),CDI(K+1))
                           CALL SWAPINT(INDEX1(K),INDEX1(K+1))
                           CALL SWAPINT(INFI(K),INFI(K+1))
                           IJ = IJ - K 
                        END DO  ! K
                        GO TO 50
                     END IF
                  END DO ! L
                  GO TO 50
               END IF
               COF(II+J) = 0.d0
            END DO ! J
            CDI(I)=0.d0
 50         II = II + I
            Y(I) = 0.d0
         END IF
         !PRINT *,'Iteration I=',I
         !CALL PRINTCOF(N,A,B,INFI,COF,INDEX1)
      END DO ! I
      !PRINT *,'Covsrt NDIM=',NDIM
      !CALL PRINTCOF(N,A,B,INFI,COF,INDEX1)
      RETURN
      END SUBROUTINE COVSRT0
	SUBROUTINE COVSRT( Nt,Nd,R,Cm,A,B,INFI,INDEX1, 
     &              INFIS,INFISD, NDIM, Y, COF,CDI )
      USE FIMOD
      USE GLOBALDATA, ONLY : EPS2
      IMPLICIT NONE
*
*     Subroutine to sort integration limits and determine Cholesky factor.
*
      INTEGER,                          INTENT(in)    :: Nt,Nd
      DOUBLE PRECISION, DIMENSION(:,:), INTENT(in)    :: R
      DOUBLE PRECISION, DIMENSION(:  ), INTENT(inout) :: Cm,A,B
      INTEGER,          DIMENSION(:  ), INTENT(inout) :: INFI
      INTEGER,          DIMENSION(:  ), INTENT(inout) :: INDEX1
      INTEGER,                          INTENT(out) :: INFIS,INFISD,NDIM
      DOUBLE PRECISION, DIMENSION(:  ), INTENT(out)   :: Y, COF, CDI
! Local variables
      INTEGER :: N,N1,I, J, K, L, M, II, IJ, JMIN,Ndleft
	INTEGER :: KK, K1, KK1, K0, Nullity
      DOUBLE PRECISION :: SUMSQ, AJ, BJ, TMP, SQTWPI, D, E,EPSL, TOL
	DOUBLE PRECISION :: AA, Ca, Pa
      DOUBLE PRECISION :: CVDIAG, AMIN, BMIN, DMIN, EMIN, YL, YU
      PARAMETER ( SQTWPI = 2.506628274631001D0, TOL = 1D-16 )
      
      EPSL=EPS2
      II     = 0
      INFIS  = 0
      INFISD = 0
      Ndim   = 0
      Ndleft = Nd
      N      = Nt+Nd
	
      IF (N.LT.10) EPSL = 1D-10
	DO I = 1, N
         IF (R(I,I).LE.TOL) INFI(I)=-1
         IF ( INFI(I) .LT. 0.) THEN
            IF (INDEX1(I).LE.Nt)  THEN 
               INFIS = INFIS+1
            ELSEIF (R(I,I).LE.TOL) THEN
               INFISD = INFISD+1
            ENDIF
         ENDIF
        
         DO J = 1, I-1
            II = II + 1
            COF(II) = R(J,I)
         END DO
			II = II+1
			COF(II) = R(I,I) !+1D-12  ! adding a nugget effect
      END DO
      COF(II+1:SIZE(COF))=0.d0
      CDI(N-INFIS-INFISD+1:SIZE(CDI))=0.d0
      !PRINT *,'COVSRT'
      !CALL PRINTCOF(N,A,B,INFI,COF,INDEX1)

      IF ( INFIS +INFISD .GE. N ) RETURN
*     
*     First move Xd to the beginning 
!         IJ = MIN(Nd,Nt)
!         DO I = N, N-IJ+1, -1
!            DO J = 1,IJ
!               CALL RCSWP( J, I, A, B, INFI, N, COF )
!               CALL SWAPINT(INDEX1(J),INDEX1(I))
!               CALL SWAPRE(Cm(J),Cm(I))
!            END DO
!         END DO
*
*     Move any doubly infinite limits of Xt to innermost positions.
*
	DO I = N, N-INFIS+1, -1
		IF ( INFI(I) .GE. 0 .OR. INDEX1(I).GT.Nt) THEN 
            DO J = 1,I-1
               IF ( INFI(J) .LT. 0 .AND. INDEX1(J).LE.Nt) THEN
                  CALL RCSWAP( J, I,N, A, B, INFI,INDEX1,Cm, COF )
                  GO TO 10
               ENDIF
            END DO
         ENDIF
 10   END DO
*     Then any redundant variables of Xd to the next innermost positions. 
      DO I = N-INFIS, N-INFIS-INFISD+1, -1
	  K=INDEX1(I)
	  IF ( R(K,K) .GT. TOL) THEN 
		 DO J = 1,I-1
			M=INDEX1(J)
			IF ( R(M,M) .LE. TOL) THEN
			   CALL RCSWAP( J, I,N, A, B, INFI,INDEX1,Cm, COF )
			   GO TO 15
			ENDIF
		 END DO
	  ENDIF
15     END DO

*     Then any doubly infinite limits of Xd to the next innermost positions. 
!         DO I = N-INFIS, N-INFIS-INFISD+1, -1
!            IF ( INFI(I) .GE. 0) THEN 
!               DO J = 1,I-1
!                  IF ( INFI(J) .LT. 0) THEN
!                     CALL RCSWP( J, I, A, B, INFI, N, COF )
!                     CALL SWAPINT(INDEX1(J),INDEX1(I))
!                     CALL SWAPRE(Cm(J),Cm(I))
!                     GO TO 15
!                  ENDIF
!               END DO
!            ENDIF
! 15      END DO
         
      

*
*     Sort remaining limits and determine Cholesky factor.
*
      KK = 0
	K  = 1
	N1  = N-INFIS-INFISD
	Nullity = 0
      DO  WHILE (K .LE. N1) 
*
*        Determine the integration limits for variable with minimum
*        expected probability and interchange that variable with Kth.
*
		K0   = K-Nullity
         DMIN = 0.d0
         EMIN = 1.d0
         JMIN = K
         CVDIAG = 0.d0
         II = KK
         DO J = K,N1
            IF ( COF(II+J) .GT. EPSL) THEN

               TMP = 0.d0  ! =  the conditional mean of Y(I) given Y(1:I-1)
               DO I = 1, K0-1
                  TMP = TMP + COF(II+I)*Y(I)
               END DO
				SUMSQ = SQRT( COF(II+J) )

               IF (INFI(J).LT.0) GO TO 30  ! May have infinite int. limits if Nd>0
               IF (INFI(J).NE.0) AJ = ( A(J) - TMP )/SUMSQ
               IF (INFI(J).NE.1) BJ = ( B(J) - TMP )/SUMSQ
 30				IF (INDEX1(J).GT.Nt) THEN
					AA = (Cm(J)+TMP)/SUMSQ  ! inflection point
					CALL EXLMS(AA,AJ,BJ,INFI(J),D,E,Ca,Pa)
				ELSE
	            CALL MVNLMS( AJ, BJ, INFI(J), D, E )
				ENDIF
               IF ( EMIN + D .GE. E + DMIN ) THEN
                  JMIN = J
                  AMIN = AJ
                  BMIN = BJ
                  DMIN = D
                  EMIN = E
                  CVDIAG = SUMSQ
               ENDIF
            ENDIF
            II = II + J
                                !               IF (J.GE.Nd.AND.Ndleft.GT.0) THEN
                                !                   Ndleft=Ndleft-1
                                !                   GOTO 40
                                !                ENDIF
         END DO        
*     
*     Compute Ith column of Cholesky factor.
*     Compute expected value for Ith integration variable (without considering the jacobian) and
*     scale Ith covariance matrix row and limits.
*
  40     IF ( CVDIAG.GT.TOL) THEN
			NDIM   = NDIM+1         ! counter for Number of relevant dimensions to integrate

	      IF ( K.LT.JMIN ) THEN
	         CALL RCSWAP( K,JMIN,N1,A,B,INFI,INDEX1,Cm,COF)
	      END IF
			KK1 = KK
			K1  = K	
			II  = KK + K0;
	      COF(II) = CVDIAG;
            CDI(K)  = CVDIAG         ! Store the diagonal element
            DO I = K0+1,K
				II = II+1;
				COF(II) = 0.d0;
			END DO
			!now II = KK+K


			DO I = K1+1,N1              
				TMP = 0.D0
				DO J = 1, K0-1
					!tmp = tmp+L(i,j).*L(k1,j)
					TMP = TMP + COF(II+J)*COF(KK1+J) 
               END DO
				!Cov(Xk,Xi|X1,X2,...Xk-1)/STD(Xk|X1,X2,...Xk-1)
				COF(II+K0)  = (COF(II+K1) - TMP)/CVDIAG  
				!  Var(Xk|X1,X2,...Xk)
				COF(II+I) = COF(II+I) - COF(II+K0) * COF(II+K0)

               
               IF (COF(II+I) .LE. EPSL) THEN  !TOL
					Nullity  = Nullity + 1
					KK = KK + K
					K  = K + 1
					IF (K .LT. I) THEN
						CALL RCSWAP(K,I,N1,A,B,INFI,INDEX1,Cm,COF)
					ENDIF
					CDI(K) = COF(KK+K0)
					A(K) = A(K)/CDI(K)
					B(K) = B(K)/CDI(K)
					IF ((CDI(K) .LT. 0.d0).AND. INFI(K).GE. 0) THEN
						CALL SWAPRE(A(K),B(K))
						IF (INFI(K).NE. 2) INFI(K) = 1-INFI(K)
					END IF
					DO J = KK+1,KK+K0
						COF(J) = COF(J)/CDI(K)
					END DO
					!CDI(K) = 0.d0
					DO J = KK+K0+1,KK+K
						COF(J) = 0.d0
					END DO
				END IF	
				II = II + I
			END DO
            

			IF ( EMIN .GT. DMIN + EPSL ) THEN
               YL = 0.d0
               YU = 0.d0
               IF (INFI(K1).GE.0) THEN
                  IF ( INFI(K1) .NE. 0 ) YL =-EXP(-0.5d0*AMIN**2)/SQTWPI
                  IF ( INFI(K1) .NE. 1 ) YU =-EXP(-0.5d0*BMIN**2)/SQTWPI
               ENDIF
               Y(K0) = ( YU - YL )/( EMIN - DMIN )
            ELSE
               SELECT CASE (INFI(K1))
               CASE (:-1) 
                  Y(K0) = 0.d0
               CASE (0) 
                  Y(K0) = BMIN
               CASE (1)
                  Y(K0) = AMIN
               CASE (2:)
                  Y(K0) = ( AMIN + BMIN )*0.5d0
               END SELECT 
            END IF

	      DO J = KK1+1, KK1+K0
               COF(J) = COF(J)/CVDIAG
            END DO

			A(K1) = A(K1)/CVDIAG
            B(K1) = B(K1)/CVDIAG

			KK = KK + K
			K  = K  + 1
         ELSE
			II = KK;
			DO I = K,N1                
				CDI(I) = COF(II+K0-1)
				A(I) = A(I)/CDI(I)
				B(I) = B(I)/CDI(I)
				IF ((CDI(I) .LT. 0.d0).AND. INFI(I).GE. 0) THEN
					CALL SWAPRE(A(I),B(I))
					IF (INFI(I).NE. 2) INFI(I) = 1-INFI(I)
				END IF
				DO J = 1,K0-1
					II = II +1
					COF(II) = COF(II)/CDI(I)
				END DO	
				
				DO J=K0,I
					II = II+1;
					COF(II) = 0.d0;  ! set the covariance to the rest to zero
				END DO
			END DO
			Nullity = N1-K0+1
			GOTO 200
			!RETURN			
         END IF
      END DO 

200	CONTINUE
	! I
      !PRINT *,'Covsrt NDIM=',NDIM
      !CALL PRINTCOF(N,A,B,INFI,COF,INDEX1)
      RETURN
      END SUBROUTINE COVSRT
*
      SUBROUTINE RCSWP( P, Q, A, B, INFIN, N, C )
      IMPLICIT NONE
*
*     Swaps rows and columns P and Q in situ, with P <= Q.
*
      DOUBLE PRECISION, DIMENSION(:),INTENT(inout):: A, B, C
      INTEGER, DIMENSION(:),INTENT(inout) :: INFIN
      INTEGER,INTENT(in) :: P, Q, N
! local variables
      INTEGER :: I, J, II, JJ
      CALL SWAPRE( A(P), A(Q) )
      CALL SWAPRE( B(P), B(Q) )
      CALL SWAPINT(INFIN(P),INFIN(Q))
      JJ = ( P*( P - 1 ) )/2
      II = ( Q*( Q - 1 ) )/2
      CALL SWAPRE( C(JJ+P), C(II+Q) )
      DO J = 1, P-1
         CALL SWAPRE( C(JJ+J), C(II+J) )
      END DO
      JJ = JJ + P
      DO I = P+1, Q-1
         CALL SWAPRE( C(JJ+P), C(II+I) )
         JJ = JJ + I
      END DO
      II = II + Q
      DO I = Q+1, N
         CALL SWAPRE( C(II+P), C(II+Q) )
         II = II + I
      END DO
      RETURN
      END SUBROUTINE RCSWP
      SUBROUTINE RCSWAP( P, Q, N, A, B, INFIN, IND, Cm, C )
      IMPLICIT NONE
* RCSWAP  Swaps rows and columns P and Q in situ, with P <= Q.
*
*
*   CALL  RCSWAP( P, Q, N, A, B, INFIN, IND, Cm, C )
*
*    P, Q  = row/column number to swap P<=Q<=N
*    N     = length of A, B
*    A,B   = lower and upper integration limit, respectively.
*    INFIN = if INFIN(I) < 0, Ith limits are (-infinity, infinity);
*            if INFIN(I) = 0, Ith limits are (-infinity, B(I)];
*            if INFIN(I) = 1, Ith limits are [A(I), infinity);
*            if INFIN(I) = 2, Ith limits are [A(I), B(I)].
*    IND   = permutation indiex vector.
*    Cm    = conditional mean
*    C     = lower triangular cholesky factor.
*
      DOUBLE PRECISION, DIMENSION(:),INTENT(inout):: A, B, C, Cm
      INTEGER, DIMENSION(:),INTENT(inout) :: INFIN,IND
      INTEGER,INTENT(in) :: P, Q, N
! local variables
      INTEGER :: I, J, II, JJ
      CALL SWAPRE( A(P), A(Q) )
      CALL SWAPRE( B(P), B(Q) )
	CALL SWAPRE( Cm(P), Cm(Q) )
      CALL SWAPINT(INFIN(P),INFIN(Q))
	CALL SWAPINT(IND(P),IND(Q))
      JJ = ( P*( P - 1 ) )/2  ! index to element A(p-1,p-1)
      II = ( Q*( Q - 1 ) )/2	! index to element A(q-1,q-1)
      CALL SWAPRE( C(JJ+P), C(II+Q) )
      DO J = 1, P-1
         CALL SWAPRE( C(JJ+J), C(II+J) )
      END DO
      JJ = JJ + P
      DO I = P+1, Q-1
         CALL SWAPRE( C(JJ+P), C(II+I) )
         JJ = JJ + I
      END DO
      II = II + Q
      DO I = Q+1, N
         CALL SWAPRE( C(II+P), C(II+Q) )
         II = II + I
      END DO
      RETURN
      END SUBROUTINE RCSWAP
	

      SUBROUTINE swapRe(m,n)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(inout) :: m,n
      DOUBLE PRECISION                :: tmp      
      tmp=m
      m=n
      n=tmp
      RETURN
      END SUBROUTINE swapRe
  
      SUBROUTINE swapint(m,n)
      IMPLICIT NONE
      INTEGER, INTENT(inout) :: m,n
      INTEGER                :: tmp 
      tmp=m
      m=n
      n=tmp
      RETURN
      END SUBROUTINE swapint
      END MODULE RIND         !******************************  
