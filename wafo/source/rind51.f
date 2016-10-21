!****************************************************************************
!NB! if compilation complains about too many continuation lines extend it.
! 
!
!  modules:   GLOBALDATA, QUAD, RIND    Version 1.0   
!
! Programs available in module RIND : 
! (NB! the GLOBALDATA and QUAD  module is also used to transport the inputs)  
!
!
! SETDATA initializes global constants explicitly:
!   
! CALL SETDATA(EPSS,EPS2,xCutOff,NINT1) 
!
!                   GLOBALDATA module :
!   EPSS,CEPSS = 1.d0 - EPSS , controlling the accuracy of indicator function
!        EPS2  = if conditional variance is less it is considered as zero
!                i.e., the variable is considered deterministic 
!      xCutOff = 5 (standard deviations by default)
!     FxCutOff = FI(xCutOff) - FI( -xCutOff), truncation of the normal CDF.  
!
!                   QUAD module:                                             
!     Nint1(i) = quadrature formulae used in integration of Xd(i)
!                implicitly determining # nodes 
!
! INITDATA initializes global constants implicitly:
!
! CALL INITDATA (speed)    
!
!        speed = 1,2,...,9 (1=slowest and most accurate,9=fastest, 
!                           but less accurate)
!
! see the GLOBALDATA and QUAD module for other constants and default values 
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
!CALL RINDD(E,S,m,xc,indI,Blo,Bup);
!
!        E = expectation/density as explained above size 1 x Nx        (out)
!        S = Covariance matrix of X=[Xt;Xd;Xc] size N x N (N=Nt+Nd+Nc) (inout)
!            NB!: out=conditional sorted Covariance matrix
!        m = the expectation of X=[Xt;Xd;Xc]   size N x 1              (in)
!       xc = values to condition on            size Nc x Nx            (in)
!     indI = vector of indices to the different barriers in the        (in) 
!            indicator function,  length NI, where   NI = Nb+1 
!            (NB! restriction  indI(1)=0, indI(NI)=Nt+Nd )
!Blo,Bup = Lower and upper barrier coefficients used to compute the  (in)
!            integration limits Hlo and Hup, respectively. 
!            size  Mb x Nb. If  Mb<Nc+1 then
!            Blo(Mb+1:Nc+1,:) is assumed to be zero. The relation 
!            to the integration limits Hlo and Hup are as follows
!
!              Hlo(i)=Blo(1,j)+Blo(2:Nc+1,j).'*xc(:,ix), 
!              Hup(i)=Bup(1,j)+Bup(2:Nc+1,j).'*xc(:,ix), 
!
!            where i=indI(j-1)+1:indI(j), j=2:NI, ix=1:Nx
!            Thus the integration limits may change with the conditional
!            variables.
!Example: 
! The indices, indI=[0 3 5], and coefficients Blo=[-inf 0], Bup=[0  inf]
!  gives   Hlo = [-inf -inf -inf 0 0]  Hup = [0 0 0 inf inf] 
!
! The GLOBALDATA and QUAD  modules are used to transport the inputs: 
!     SCIS = 0 => integr. all by quadrature (default)
!            1 => integr. Nt-Nj variables in indicator by SCIS
!                 and the rest by quadrature
!            2 => integr.  all by SCIS 
!      NIT = 0,1,2..., maximum # of iterations/integrations done by quadrature 
!            to calculate the indicator function (default NIT=2)  
!jacobdef = defines the jacobian used (default jacobdef=0)
!  NB!  the size information below must be set before calling RINDD 
!       Nx = # different xc
!       Nt = length of Xt
!       Nd = length of Xd 
!       Nc = length of Xc
!      Ntd = Nt+Nd
!     Ntdc = Nt+Nd+Nc
!       Mb
!       NI
!       Nj = # of variables in indicator integrated directly like the
!            variables in the jacobian (default 0)
!            The order of integration between Xd and Nj of  Xt is done in 
!            decreasing order of conditional variance.
!      Njj = # of variables in indicator integrated directly like the
!            variables in the jacobian (default 0)
!            The Njj variables of Xt is integrated after Xd and Nj of Xt  
!            also in decreasing order of conditional variance. (Not implemented yet)
!
! (Recommended limitations Nx,Nt<101, Nd<7 and NIT,Nc<11) 
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
! Ocean Engineering                                                    (RINDXXX)
!
! R. Ambartzumian, A. Der Kiureghian, V. Ohanian and H.
! Sukiasian (1998)
! "Multinormal probabilities by sequential conditioned 
!  importance sampling: theory and application"             (RINDSCIS, MNORMPRB)
! Probabilistic Engineering Mechanics, Vol. 13, No 4. pp 299-308  
!
! Alan Genz (1992)
! 'Numerical Computation of Multivariate Normal Probabilites'
! J. computational Graphical Statistics, Vol.1, pp 141--149
!
! William H. Press, Saul Teukolsky, 
! William T. Wetterling and Brian P. Flannery (1997)
! "Numerical recipes in Fortran 77", Vol. 1, pp 55-63, 299--305  (SVDCMP,SOBSEQ)
!
! Igor Rychlik and Georg Lindgren (1993)
! "Crossreg - A technique for first passage and wave density analysis" (RINDXXX)
! Probability in the Engineering and informational Sciences,
! Vol 7, pp 125--148
!
! Igor Rychlik (1992)
! "Confidence bands for linear regressions"                      (RIND2,RINDNIT)
! Commun. Statist. -simula., Vol 21,No 2, pp 333--352
!
!  Jagdish K. Patel and Campbell B. Read (1982)
! "Handbook of the normal distribution",
! marcel dekker inc,New York - Basel, Vol. 40, pp 293--300      (NORM2DPRB, THL)   
!
! Donald E. Knuth (1973) "The art of computer programming,",
! Vol. 3, pp 84-  (sorting and searching)                               (SORTRE)

! Tested on:  DIGITAL UNIX Fortran90 compiler
!             PC pentium II with Lahey Fortran90 compiler
!             Solaris with SunSoft F90 compiler Version 1.0.1.0  (21229283) 
! History:
!revised by I.R. 27.01.2000, Removed bugs in RINDNIT (There where some returns
!           without deallocating some variables. A misco error in THL, leading
!           to floating invalid on alpha has been repaired by seting value=zero.
!           Probably there is an error somehere making variable "value" to behave badly.
!Revised by IR. 03.01.2000 Bug in C1C2 fixed      
!revised by I.R. 27.12.1999, New name RIND52.f
!          I have changed assumption about deterministic variables. Those have now 
!          variances equal EPS2 not zero and have consequences for C1C2 and on some
!          places in RINDND. The effect is that bariers becomes fuzzy (not sharp)
!          and prevents for discountinuities due to numerical errors of order 1E-16.
!          The program RIND0 is removed making the structure of program simpler. 
!          We have still a problem when variables in indicator become
!          deterministic before conditioning on derivatives in Xd and even Xd and Xc.
!revised by  Igor Rychlik 06.12.1999.    New name RIND50.f 
!       - The main differance to RIND49.f is that the output contains now
!         both the upper- and the lower-bound for the computed intensity.
!         For NIT=0, the program is much slower then using rind49.f. 
!         For higher NIT-values this difference is negligable relative to
!         the increase of the information. If only en estimate is needed use RIND49.f.
!         At present it seems that the NIT=-2 is to recomend for sp2tthpdf program
!         since some irregularity is smoothed out by integration of halfwavelength.
!         The density can be ploted together with the upper and lower bound. 
!revised by  Igor Rychlik 01.12.1999  New name RIND49.f
!       - changed RINDNIT and ARGP0 in order to exclude
!         irrelevant variables (such that probability of beeing 
!         between bariers is 1.) All computations related to NIT
!         are moved to RINDNIT (removing RIND2,RIND3). This caused some changes 
!         in RIND0,RINDDN. Furthermore RINDD1 is removed and moved          
!         some parts of it to RINDDN. This made program few seconds slower. The lower 
!         bound in older ARGP0 programs contained logical error - corrected.
!revised by Per A. Brodtkorb 08.11.1999
!       - fixed a bug in rinddnd 
!          new line: CmNew(Nst+1:Nsd-1)= Cm(Nst+1:Nsd-1)
!revised by Per A. Brodtkorb 28.10.1999
!       - fixed a bug in rinddnd
!       - changed rindscis, mnormprb 
!       - added MVNFUN, MVNFUN2
!       - replaced CVaccept with RelEps 
!revised by Per A. Brodtkorb 27.10.1999
!       - changed NINT to NINT1 due to naming conflict with an intrinsic of the same name
!revised by Per A. Brodtkorb 25.10.1999
!       - added an alternative FIINV for use in rindscis and mnormprb 
!revised by Per A. Brodtkorb 13.10.1999
!       - added useMIDP for use in rindscis and mnormprb 
!
!revised by Per A. Brodtkorb 22.09.1999
!       - removed all underscore letters due to
!         problems with  SunSoft F90 compiler
!         (i.e. changed  GLOBAL_DATA to GLOBALDATA etc.)
!revised by Per A. Brodtkorb 09.09.1999
!       - added sobseq: Sobol sequence (quasi random numbers)
!              an alternative to random_number in RINDSCIS and mnormprb
!revised by Per A. Brodtkorb 07.09.1999
!       - added pythag,svdcmp,sortre
!       - added RINDSCIS: evaluating multinormal integrals by SCIS
!              condsort3: prepares BIG for use with RINDSCIS and mnormprb
!revised by Per A. Brodtkorb 03.09.1999
!       - added mnormprb: evaluating multinormal probabilities by SCIS
!            See globaldata for SCIS
! revised by Per A. Brodtkorb 01.09.1999
!       - increased the default NUGGET from 1.d-12 to 1.d-8
!       - also set NUGGET depending on speed in INITDATA
! revised by Per A. Brodtkorb 27.08.1999
!       - changed rindnit,rind2: 
!         enabled option to do the integration faster/(smarter?).
!         See GLOBALDATA for XSPLT
! revised by Per A. Brodtkorb 17.08.1999
!       - added THL, norm2dprb not taken in to use
!         due to some mysterious floating invalid
!         occuring from time to time in norm2dprb (on DIGITAL unix)
! revised by Per A. Brodtkorb 02.08.1999
!       - updated condsort
!       - enabled the use of C1C2 in rinddnd
! revised by Per A. Brodtkorb 14.05.1999
!       - updated to fortran90
!       - enabled recursive calls
!       - No limitations on size of the inputs
!       - fixed some bugs
!       - added some additonal checks
!       - added Hermite, Laguerre quadratures for alternative integration
!       - rewritten CONDSORT, conditional covariance matrix in upper 
!         triangular. 
!       - RINDXXX routines only work on the upper triangular
!         of the covariance matrix
!       - Added a Nugget effect to the covariance matrix in order  
!         to ensure the conditioning is not corrupted by numerical errors    
!       - added the option to condsort Nj variables of Xt, i.e.,
!         enabling direct integration like the integration of Xd
! by  Igor Rychlik 29.10.1998 (PROGRAM RIND11 --- Version 1.0)
!         which was a revision of program RIND from 3.9.1993 - the program that
!         is used in wave_t and wave_t2 programs.

!*********************************************************************


      MODULE GLOBALDATA
      IMPLICIT NONE             
                      ! Constants determining accuracy of integration
                      !-----------------------------------------------
                      !if the conditional variance are less than: 
      DOUBLE PRECISION :: EPS2=1.d-4    !- EPS2, the variable is 
                                        !  considered deterministic 
      DOUBLE PRECISION :: XCEPS2=1.d-10 ! if Var(Xc) is less return NaN
      DOUBLE PRECISION :: EPSS = 5.d-5  ! accuracy of Indicator 
      DOUBLE PRECISION :: CEPSS=0.99995 ! accuracy of Indicator 
      DOUBLE PRECISION :: EPS0 = 5.d-5 ! used in GAUSSLE1 to implicitly 
                                       ! determ. # nodes  
   
      DOUBLE PRECISION :: fxcEpss=1.d-20 ! if less do not compute E(...|Xc)
      DOUBLE PRECISION :: xCutOff=5.d0  ! upper/lower truncation limit of the 
                                       ! normal CDF 
      DOUBLE PRECISION :: xCutOff2=25.d0 !=xCutOff^2 maximum exponent in 
                                       ! normal PDF, to avoid underflow
                                       ! FxCutOff=FI(xCutOff)-FI(-CxCutOff), 
      DOUBLE PRECISION :: FxCutOff  = 0.99999942669686d0 
      DOUBLE PRECISION :: CFxCutOff = 5.733031438470704d-7       ! 1-FxCutOff, 
      DOUBLE PRECISION :: XSMALL=4.2D-16, XMAX=8.29287554336168D0 !cut off parameters to FI
                                ! Nugget>0: Adds a small value to diagonal 
                                ! elements of the covariance matrix to ensure
                                ! that the inversion is  not corrupted by  
                                ! round off errors. 
                                ! Good choice might be 1e-8 
      DOUBLE PRECISION :: NUGGET=1.d-8 ! Obs NUGGET must be smaller then EPS2
    
!parameters controlling the performance of RINDSCIS and MNORMPRB:
      INTEGER :: SCIS=0         !=0 => integr. all by quadrature 
                                !=1 => integr. Nt-Nj variables in indicator 
                                !	   by SCIS and the rest by quadrature
                                !=2 => integr.  all by SCIS
      INTEGER :: NSIMmax = 10000 ! maximum number of simulations
      INTEGER :: NSIMmin = 10    ! minimum number of simulations
      INTEGER :: rateLHD = 20    ! rateLhd*Nstoc = size of LHD matrix
      INTEGER :: Ntscis  = 0     ! Ntscis=Nt-Nj-Njj when SCIS>0 Ntscis=0 otherwise
      DOUBLE PRECISION :: RelEps = 0.08 ! Relative error, i.e. if 
                                !2.5*STD(XIND)/XIND is less we accept the estimate
                                ! The following may be allocated outside RINDD
                                ! if one wants the coefficient of variation, i.e.
                                ! STDEV(XIND)/XIND when SCIS=2.  (NB: size Nx)                              
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: COV
      integer :: COVix ! counting variable for COV
      LOGICAL :: MLHD=.true. ! Modified Latin Hypercube Design when sampling
                                ! WITH SCIS used in rindscis and mnormprb
      LOGICAL :: useMIDP=.true. ! use midpoint of cell instead randomly within cell
                                 ! used on rindscis and mnormprb
      LOGICAL :: useC1C2=.true. ! use C1C2 in rindscis,mnormprb 
      LOGICAL :: C1C2det=.true. ! use C1C2 only on the variables that becomes 
                                ! deterministic after conditioning on X(N)
                                ! used in rinddnd rindd1 and rindscis mnormprb

!parameters controlling performance of quadrature integration:
                ! if Hup>=xCutOff AND Hlo<-XSPLT OR
                !    Hup>=XSPLT AND Hl0<=-xCutOff then
                !  do a different integration to increase speed
                ! in rind2 and rindnit. This give slightly different 
                ! results
                ! DEFAULT 5 =xCutOff => do the same integration allways 
      DOUBLE PRECISION :: XSPLT = 5.d0 ! DEFAULT XSPLT= 5 =xCutOff 
      INTEGER :: NIT=2       ! NIT=maximum # of iterations/integrations by
                                ! quadrature used to calculate the indicator function 

      INTEGER :: jacobdef=0    ! jacobdef determines the jacobian used

                                ! size information of the covariance matrix BIG
                                ! Nt,Nd,....Ntd,Nx must be set before calling 
                                ! RINDD.  NsXtmj, NsXdj is set in RINDD 
      INTEGER :: Nt,Nd,Nc,Ntdc,Ntd,Nx
                                ! Constants determines how integration is done
      INTEGER :: Nj=0,Njj=0  ! Njj is not implemented yet
                                ! size information of indI, Blo,Bup
                                ! Blo/Bup size Mb x NI-1
                                ! indI vector of length NI
      INTEGER :: NI,Mb       ! must be set before calling RINDD
     
                                ! The following is allocated in RINDD
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: SQ 
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Hlo,Hup 
      INTEGER,          DIMENSION(:), ALLOCATABLE :: index1,xedni,indXtd
      INTEGER,          DIMENSION(:), ALLOCATABLE :: NsXtmj, NsXdj

                                ! global constants
      DOUBLE PRECISION, PARAMETER :: SQTWOPI1=3.9894228040143d-1 !=1/sqrt(2*pi)
      DOUBLE PRECISION, PARAMETER :: SQPI1=5.6418958354776d-1    !=1/sqrt(pi)
      DOUBLE PRECISION, PARAMETER :: SQPI= 1.77245385090552d0    !=sqrt(pi)
      DOUBLE PRECISION, PARAMETER :: SQTWO=1.41421356237310d0    !=sqrt(2)
      DOUBLE PRECISION, PARAMETER :: SQTWO1=0.70710678118655d0   !=1/sqrt(2)
      DOUBLE PRECISION, PARAMETER :: PI1=0.31830988618379d0      !=1/pi
      DOUBLE PRECISION, PARAMETER :: PI= 3.14159265358979D0      !=pi
      DOUBLE PRECISION, PARAMETER :: TWOPI=6.28318530717958D0    !=2*pi
      END MODULE GLOBALDATA


      MODULE QUAD
      IMPLICIT NONE         ! Quadratures available: Legendre,Hermite,Laguerre

      INTEGER, PARAMETER :: PMAX=24       ! maximum # nodes
      INTEGER, PARAMETER :: sizNint=13    ! size of Nint1
      INTEGER                 :: minQNr=1      ! minimum quadrature number
                                                  ! used in GaussLe1, Gaussle2
      INTEGER                 :: Le2QNr=8      ! quadr. number used in  rind2,rindnit
      INTEGER, DIMENSION(sizNint) :: Nint1 ! use quadr. No. Nint1(i) in 
                                                  ! integration of Xd(i)

                                ! # different quadratures stored for :
                                !------------------------------------- 
      INTEGER,PARAMETER  :: NLeW=13 ! Legendre  
      INTEGER,PARAMETER  :: NHeW=13 ! Hermite   
      INTEGER,PARAMETER  :: NLaW=13 ! Laguerre 
                                ! Quadrature Number stored for :
                                !------------------------------------- 
      INTEGER, DIMENSION(NLeW) :: LeQNr    ! Legendre  
      INTEGER, DIMENSION(NHeW) :: HeQNr    ! Hermite  
      INTEGER, DIMENSION(NLaW) :: LaQNr    ! Laguerre 
      PARAMETER (LeQNr=(/ 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 16, 20, 24 /)) 
      PARAMETER (HeQNr=(/ 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 16, 20, 24 /)) 
      PARAMETER (LaQNr=(/ 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 16, 20, 24 /)) 


                            ! The indices to the weights & nodes stored for: 
                            !------------------------------------------------
      INTEGER, DIMENSION(NLeW+1) :: LeIND  !Legendre
      INTEGER, DIMENSION(NHeW+1) :: HeIND  !Hermite
      INTEGER, DIMENSION(NLaW+1) :: LaIND  !Laguerre 
      PARAMETER (LeIND=(/0,2,5,9,14,20,27,35,44,54,66,82,102,126/)) !Legendre
      PARAMETER (HeIND=(/0,2,5,9,14,20,27,35,44,54,66,82,102,126/)) !Hermite
      PARAMETER (LaIND=(/0,2,5,9,14,20,27,35,44,54,66,82,102,126/)) !Laguerre 

                            !------------------------------------------------
      DOUBLE PRECISION,  DIMENSION(126) :: LeBP,LeWF,HeBP,HeWF 
      DOUBLE PRECISION,  DIMENSION(126) :: LaBP0,LaWF0,LaBP5,LaWF5

!The Hermite Quadrature integrates an integral of the form
!        inf                         n
!       Int (exp(-x^2) F(x)) dx  =  Sum  wf(j)*F( bp(j) )
!       -Inf                        j=1   
!The Laguerre Quadrature integrates an integral of the form
!        inf                               n
!       Int (x^alpha exp(-x) F(x)) dx  =  Sum  wf(j)*F( bp(j) )
!         0                               j=1   
! weights stored here are for alpha=0 and alpha=-0.5

                             ! initialize Legendre weights, wf,  and nodes, bp
      PARAMETER ( LeWF = (/ 1.d0, 1.d0, 0.555555555555556d0,             
     *     0.888888888888889d0, 0.555555555555556d0, 
     *     0.347854845137454d0, 0.652145154862546d0,
     *     0.652145154862546d0, 0.347854845137454d0,         
     *     0.236926885056189d0, 0.478628670499366d0, 
     *     0.568888888888889d0, 0.478628670499366d0,
     *     0.236926885056189d0, 0.171324492379170d0,         
     *     0.360761573048139d0, 0.467913934572691d0,
     *     0.467913934572691d0, 0.360761573048139d0,
     *     0.171324492379170d0, 0.129484966168870d0,         
     *     0.279705391489277d0, 0.381830050505119d0,
     *     0.417959183673469d0, 0.381830050505119d0,
     *     0.279705391489277d0, 0.129484966168870d0,         
     *     0.101228536290376d0, 0.222381034453374d0, 
     *     0.313706645877887d0, 0.362683783378362d0,
     *     0.362683783378362d0, 0.313706645877887d0,         
     *     0.222381034453374d0, 0.101228536290376d0,
     *     0.081274388361574d0, 0.180648160694857d0,
     *     0.260610696402935d0, 0.312347077040003d0,         
     *     0.330239355001260d0, 0.312347077040003d0,
     *     0.260610696402935d0, 0.180648160694857d0,
     *     0.081274388361574d0, 0.066671344308688d0,         
     *     0.149451349150581d0, 0.219086362515982d0,
     *     0.269266719309996d0, 0.295524224714753d0,
     *     0.295524224714753d0, 0.269266719309996d0,         
     *     0.219086362515982d0, 0.149451349150581d0,
     *     0.066671344308688d0, 0.047175336386512d0, 
     *     0.106939325995318d0, 0.160078328543346d0,         
     *     0.203167426723066d0, 0.233492536538355d0, 
     *     0.249147048513403d0, 0.249147048513403d0,
     *     0.233492536538355d0, 
     *     0.203167426723066d0,       0.160078328543346d0, 
     *     0.106939325995318d0,       0.047175336386512d0,         
     *     0.027152459411754094852d0, 0.062253523938647892863d0,
     *     0.095158511682492784810d0, 0.124628971255533872052d0,
     *     0.149595988816576732081d0, 0.169156519395002538189d0,
     *     0.182603415044923588867d0, 0.189450610455068496285d0,
     *     0.189450610455068496285d0, 0.182603415044923588867d0,
     *     0.169156519395002538189d0, 0.149595988816576732081d0,
     *     0.124628971255533872052d0, 0.095158511682492784810d0,
     *     0.062253523938647892863d0, 0.027152459411754094852d0,
     *     0.017614007139152118312d0, 0.040601429800386941331d0,
     *     0.062672048334109063570d0, 0.083276741576704748725d0,
     *     0.101930119817240435037d0, 0.118194531961518417312d0,
     *     0.131688638449176626898d0, 0.142096109318382051329d0,
     *     0.149172986472603746788d0, 0.152753387130725850698d0,
     *     0.152753387130725850698d0, 0.149172986472603746788d0,
     *     0.142096109318382051329d0, 0.131688638449176626898d0,
     *     0.118194531961518417312d0, 0.101930119817240435037d0,
     *     0.083276741576704748725d0, 0.062672048334109063570d0,
     *     0.040601429800386941331d0, 0.017614007139152118312d0,            
     *     0.012341229799987199547d0, 0.028531388628933663181d0,
     *     0.044277438817419806169d0, 0.059298584915436780746d0,
     *     0.073346481411080305734d0, 0.086190161531953275917d0,
     *     0.097618652104113888270d0, 0.107444270115965634783d0,
     *     0.115505668053725601353d0, 0.121670472927803391204d0,
     *     0.125837456346828296121d0, 0.127938195346752156974d0,
     *     0.127938195346752156974d0, 0.125837456346828296121d0,
     *     0.121670472927803391204d0, 0.115505668053725601353d0,
     *     0.107444270115965634783d0, 0.097618652104113888270d0,
     *     0.086190161531953275917d0, 0.073346481411080305734d0,
     *     0.059298584915436780746d0, 0.044277438817419806169d0,
     *     0.028531388628933663181d0, 0.012341229799987199547d0 /))

      PARAMETER ( LeBP=(/-0.577350269189626d0,0.577350269189626d0,
     *    -0.774596669241483d0,0.d0,
     *     0.774596669241483d0, -0.861136311594053d0, 
     *    -0.339981043584856d0, 0.339981043584856d0,
     *     0.861136311594053d0, -0.906179845938664d0,
     *    -0.538469310105683d0, 0.d0,
     *     0.538469310105683d0,  0.906179845938664d0, 
     *    -0.932469514203152d0,  -0.661209386466265d0,
     *    -0.238619186083197d0,  0.238619186083197d0,
     *     0.661209386466265d0,  0.932469514203152d0, 
     *    -0.949107912342759d0,-0.741531185599394d0,
     *    -0.405845151377397d0,  0.d0,
     *     0.405845151377397d0,  0.741531185599394d0, 
     *     0.949107912342759d0,   -0.960289856497536d0, 
     *    -0.796666477413627d0, -0.525532409916329d0,
     *    -0.183434642495650d0,  0.183434642495650d0, 
     *     0.525532409916329d0,    0.796666477413627d0,
     *     0.960289856497536d0, -0.968160239507626d0,
     *    -0.836031107326636d0, -0.613371432700590d0,
     *    -0.324253423403809d0,  0.d0,
     *     0.324253423403809d0,  0.613371432700590d0,
     *     0.836031107326636d0,    0.968160239507626d0, 
     *    -0.973906528517172d0, -0.865063366688985d0,
     *    -0.679409568299024d0, -0.433395394129247d0, 
     *    -0.148874338981631d0,    0.148874338981631d0, 
     *     0.433395394129247d0,  0.679409568299024d0,
     *     0.865063366688985d0,  0.973906528517172d0,
     *    -0.981560634246719d0,    -0.904117256370475d0, 
     *    -0.769902674194305d0, -0.587317954286617d0,
     *    -0.367831498198180d0, -0.125233408511469d0,
     *     0.125233408511469d0, 0.367831498198180d0,  
     *     0.587317954286617d0, 0.769902674194305d0,
     *     0.904117256370475d0,  0.981560634246719d0,
     *    -0.989400934991649932596d0,
     *    -0.944575023073232576078d0, -0.865631202387831743880d0,
     *    -0.755404408355003033895d0, -0.617876244402643748447d0,
     *    -0.458016777657227386342d0, -0.281603550779258913230d0,
     *    -0.095012509837637440185d0,  0.095012509837637440185d0,
     *     0.281603550779258913230d0,  0.458016777657227386342d0,
     *     0.617876244402643748447d0,  0.755404408355003033895d0,
     *     0.865631202387831743880d0,  0.944575023073232576078d0,
     *     0.989400934991649932596d0, -0.993128599185094924786d0,
     *    -0.963971927277913791268d0, -0.912234428251325905868d0,
     *    -0.839116971822218823395d0, -0.746331906460150792614d0,
     *    -0.636053680726515025453d0, -0.510867001950827098004d0,
     *    -0.373706088715419560673d0, -0.227785851141645078080d0,
     *     -0.076526521133497333755d0,  0.076526521133497333755d0,
     *      0.227785851141645078080d0,  0.373706088715419560673d0,
     *      0.510867001950827098004d0,  0.636053680726515025453d0,
     *      0.746331906460150792614d0,  0.839116971822218823395d0,
     *      0.912234428251325905868d0,
     *      0.963971927277913791268d0,  0.993128599185094924786d0,
     *     -0.995187219997021360180d0, -0.974728555971309498198d0,
     *     -0.938274552002732758524d0, -0.886415527004401034213d0,
     *     -0.820001985973902921954d0, -0.740124191578554364244d0,
     *     -0.648093651936975569252d0, -0.545421471388839535658d0,
     *     -0.433793507626045138487d0, -0.315042679696163374387d0,
     *     -0.191118867473616309159d0, -0.064056892862605626085d0,
     *      0.064056892862605626085d0,  0.191118867473616309159d0,
     *      0.315042679696163374387d0,  0.433793507626045138487d0,
     *      0.545421471388839535658d0,  0.648093651936975569252d0,
     *      0.740124191578554364244d0,  0.820001985973902921954d0,
     *      0.886415527004401034213d0,  0.938274552002732758524d0,
     *      0.974728555971309498198d0,  0.995187219997021360180d0 /))    

                                ! initialize Hermite weights in HeWF and 
                                ! nodes in HeBP
                                ! NB! the relative error of these numbers  
                                ! are less than 10^-15 
      PARAMETER (HeWF =(/ 8.8622692545275816d-1,
     *     8.8622692545275816d-1,
     *     2.9540897515091930d-1,   1.1816359006036770d0,
     *     2.9540897515091930d-1,   8.1312835447245310d-2,
     *     8.0491409000551251d-1,   8.0491409000551295d-1,
     *     8.1312835447245213d-2,   1.9953242059045910d-2,
     *     3.9361932315224146d-1,   9.4530872048294134d-1,
     *     3.9361932315224102d-1,   1.9953242059045962d-2,
     *     4.5300099055088378d-3,   1.5706732032285636d-1,
     *     7.2462959522439319d-1,   7.2462959522439241d-1,
     *     1.5706732032285681d-1,   4.5300099055088534d-3,
     *     9.7178124509952175d-4,   5.4515582819126975d-2,
     *     4.2560725261012805d-1,   8.1026461755680768d-1,
     *     4.2560725261012783d-1,   5.4515582819126975d-2,
     *     9.7178124509951828d-4,   1.9960407221136729d-4,
     *     1.7077983007413571d-2,   2.0780232581489183d-1,
     *     6.6114701255824082d-1,   6.6114701255824138d-1,
     *     2.0780232581489202d-1,   1.7077983007413498d-2,
     *     1.9960407221136775d-4,   3.9606977263264446d-5,
     *     4.9436242755369411d-3,   8.8474527394376654d-2,
     *     4.3265155900255586d-1,   7.2023521560605108d-1,
     *     4.3265155900255559d-1,   8.8474527394376543d-2,
     *     4.9436242755369350d-3,   3.9606977263264324d-5,
     *     7.6404328552326139d-6,   1.3436457467812229d-3,
     *     3.3874394455481210d-2,   2.4013861108231502d-1,
     *     6.1086263373532623d-1,   6.1086263373532546d-1,
     *     2.4013861108231468d-1,   3.3874394455480884d-2,
     *     1.3436457467812298d-3,   7.6404328552325919d-6,
     *     2.6585516843562997d-7,   8.5736870435879089d-5,
     *     3.9053905846291028d-3,   5.1607985615883860d-2,
     *     2.6049231026416092d-1,   5.7013523626247820d-1,
     *     5.7013523626248030d-1,   2.6049231026416109d-1,
     *     5.1607985615883846d-2,   3.9053905846290530d-3,
     *     8.5736870435878506d-5,   2.6585516843562880d-7,
     *     2.6548074740111735d-10,  2.3209808448651987d-7,
     *     2.7118600925379007d-5,   9.3228400862418819d-4,
     *     1.2880311535509989d-2,   8.3810041398985652d-2,
     *     2.8064745852853318d-1,   5.0792947901661278d-1,
     *     5.0792947901661356d-1,   2.8064745852853334d-1,
     *     8.3810041398985735d-2,   1.2880311535510015d-2,
     *     9.3228400862418407d-4,   2.7118600925378956d-5,
     *     2.3209808448651966d-7,   2.6548074740111787d-10,
     *     2.2293936455342015d-13,  4.3993409922730765d-10,
     *     1.0860693707692910d-7,   7.8025564785320463d-6,
     *     2.2833863601635403d-4,   3.2437733422378719d-3,
     *     2.4810520887463536d-2,   1.0901720602002360d-1,
     *     2.8667550536283382d-1,   4.6224366960061047d-1,
     *     4.6224366960061070d-1,   2.8667550536283398d-1,
     *     1.0901720602002325d-1,   2.4810520887463588d-2,
     *     3.2437733422378649d-3,   2.2833863601635316d-4,
     *     7.8025564785321005d-6,   1.0860693707692749d-7,
     *     4.3993409922731370d-10,  2.2293936455342167d-13,
     *     1.6643684964891124d-16,  6.5846202430781508d-13,
     *     3.0462542699875022d-10,  4.0189711749413878d-8,
     *     2.1582457049023452d-6,   5.6886916364043773d-5,
     *     8.2369248268841073d-4,   7.0483558100726748d-3,
     *     3.7445470503230736d-2,   1.2773962178455966d-1,
     *     2.8617953534644325d-1,   4.2693116386869828d-1,
     *     4.2693116386869912d-1,   2.8617953534644286d-1,
     *     1.2773962178455908d-1,   3.7445470503230875d-2,
     *     7.0483558100726844d-3,   8.2369248268842027d-4,
     *     5.6886916364044037d-5,   2.1582457049023460d-6,
     *     4.0189711749414963d-8,   3.0462542699876118d-10,
     *     6.5846202430782225d-13,  1.6643684964889408d-16 /))
      
                                !hermite nodes
      PARAMETER (HeBP = (/  -7.07106781186547572d-1,
     *     7.0710678118654752d-1,   -1.2247448713915894d0,
     *     0.d0,                     1.2247448713915894d0,
     *     -1.6506801238857845d0,   -5.2464762327529035d-1,
     *     5.2464762327529035d-1,    1.6506801238857845d0,
     *     -2.0201828704560869d0,   -9.5857246461381806d-1,
     *     0.d0,                     9.5857246461381851d-1,
     *     2.0201828704560860d0,    -2.3506049736744918d0,
     *     -1.3358490740136963d0,   -4.3607741192761629d-1,
     *     4.3607741192761657d-1,    1.3358490740136963d0,
     *     2.3506049736744927d0,    -2.6519613568352334d0,
     *     -1.6735516287674728d0,   -8.1628788285896470d-1,
     *     0.d0,                     8.1628788285896470d-1,
     *     1.6735516287674705d0,     2.6519613568352325d0,
     *     -2.9306374202572423d0,   -1.9816567566958434d0,
     *     -1.1571937124467806d0,   -3.8118699020732233d-1,
     *     3.8118699020732211d-1,    1.1571937124467804d0,
     *     1.9816567566958441d0,     2.9306374202572423d0,
     *     -3.1909932017815290d0,   -2.2665805845318436d0,
     *     -1.4685532892166682d0,   -7.2355101875283812d-1,
     *     0.d0,                     7.2355101875283756d-1,
     *     1.4685532892166657d0,     2.2665805845318405d0,
     *     3.1909932017815281d0,    -3.4361591188377387d0,
     *     -2.5327316742327906d0,   -1.7566836492998805d0,
     *     -1.0366108297895140d0,   -3.4290132722370548d-1,
     *     3.4290132722370464d-1,    1.0366108297895136d0,
     *     1.7566836492998834d0,     2.5327316742327857d0,
     *     3.4361591188377396d0,    -3.8897248978697796d0,
     *     -3.0206370251208856d0,   -2.2795070805010567d0,
     *     -1.5976826351526050d0,   -9.4778839124016290d-1,
     *     -3.1424037625435908d-1,   3.1424037625435935d-1,
     *     9.4778839124016356d-1,    1.5976826351526054d0,
     *     2.2795070805010602d0,     3.0206370251208905d0,
     *     3.8897248978697831d0,    -4.6887389393058214d0,
     *     -3.8694479048601251d0,   -3.1769991619799582d0,
     *     -2.5462021578474765d0,   -1.9517879909162541d0,
     *     -1.3802585391988809d0,   -8.2295144914465523d-1,
     *     -2.7348104613815177d-1,   2.7348104613815244d-1,
     *     8.2295144914465579d-1,    1.3802585391988802d0,
     *     1.9517879909162534d0,     2.5462021578474801d0,
     *     3.1769991619799565d0,     3.8694479048601265d0,
     *     4.6887389393058196d0,    -5.3874808900112274d0,   
     *     -4.6036824495507513d0,   -3.9447640401156296d0,   
     *     -3.3478545673832154d0,   -2.7888060584281300d0,
     *     -2.2549740020892721d0,   -1.7385377121165839d0,
     *     -1.2340762153953209d0,   -7.3747372854539361d-1,
     *     -2.4534070830090124d-1,   2.4534070830090149d-1,
     *     7.3747372854539439d-1,    1.2340762153953226d0,
     *     1.7385377121165866d0,     2.2549740020892770d0,
     *     2.7888060584281282d0,     3.3478545673832105d0,
     *     3.9447640401156230d0,     4.6036824495507398d0,
     *     5.3874808900112274d0,    -6.0159255614257390d0,
     *     -5.2593829276680442d0,   -4.6256627564237904d0,
     *     -4.0536644024481472d0,   -3.5200068130345219d0,
     *     -3.0125461375655647d0,   -2.5238810170114276d0,
     *     -2.0490035736616989d0,   -1.5842500109616944d0,
     *     -1.1267608176112460d0,   -6.7417110703721150d-1,
     *     -2.2441454747251538d-1,   2.2441454747251532d-1,
     *     6.7417110703721206d-1,    1.1267608176112454d0,
     *     1.5842500109616939d0,     2.0490035736616958d0,
     *     2.5238810170114281d0,     3.0125461375655687d0,
     *     3.5200068130345232d0,     4.0536644024481499d0,
     *     4.6256627564237816d0,     5.2593829276680353d0,
     *     6.0159255614257550d0 /))
                          !initialize Laguerre weights and nodes (basepoints)
                          ! for alpha=0  
                          ! NB! the relative error of these numbers 
                          ! are less than 10^-15           
      PARAMETER (LaWF0= (/  8.5355339059327351d-1,
     *     1.4644660940672624d-1, 7.1109300992917313d-1, 
     *     2.7851773356924092d-1,  1.0389256501586137d-2, 
     *     6.0315410434163386d-1, 
     *     3.5741869243779956d-1,  3.8887908515005364d-2, 
     *     5.3929470556132730d-4,  5.2175561058280850d-1, 
     *     3.9866681108317570d-1,  7.5942449681707588d-2, 
     *     3.6117586799220489d-3,  2.3369972385776180d-5, 
     *     4.5896467394996360d-1,  4.1700083077212080d-1, 
     *     1.1337338207404497d-1,  1.0399197453149061d-2, 
     *     2.6101720281493249d-4,  8.9854790642961944d-7, 
     *     4.0931895170127397d-1,  4.2183127786171964d-1, 
     *     1.4712634865750537d-1, 
     *     2.0633514468716974d-2,  1.0740101432807480d-3, 
     *     1.5865464348564158d-5,  3.1703154789955724d-8, 
     *     3.6918858934163773d-1,  4.1878678081434328d-1, 
     *     1.7579498663717152d-1,  3.3343492261215649d-2, 
     *     2.7945362352256712d-3,  9.0765087733581999d-5, 
     *     8.4857467162725493d-7,  1.0480011748715038d-9, 
     *     3.3612642179796304d-1,  4.1121398042398466d-1, 
     *     1.9928752537088576d0,   4.7460562765651609d-2, 
     *     5.5996266107945772d-3,  3.0524976709321133d-4, 
     *     6.5921230260753743d-6,  4.1107693303495271d-8, 
     *     3.2908740303506941d-11, 
     *     3.0844111576502009d-1,  4.0111992915527328d-1, 
     *     2.1806828761180935d-1,  6.2087456098677683d-2, 
     *     9.5015169751810902d-3,  7.5300838858753855d-4, 
     *     2.8259233495995652d-5,  4.2493139849626742d-7, 
     *     1.8395648239796174d-9,  9.9118272196090085d-13,
     &     2.6473137105544342d-01,
     &     3.7775927587313773d-01, 2.4408201131987739d-01,
     &     9.0449222211681030d-02, 2.0102381154634138d-02,
     &     2.6639735418653122d-03, 2.0323159266299895d-04,
     &     8.3650558568197802d-06, 1.6684938765409045d-07,
     &     1.3423910305150080d-09, 3.0616016350350437d-12,
     &     8.1480774674261369d-16, 2.0615171495780091d-01,
     &     3.3105785495088480d-01, 2.6579577764421392d-01,
     &     1.3629693429637740d-01, 4.7328928694125222d-02,
     &     1.1299900080339390d-02, 1.8490709435263156d-03,
     &     2.0427191530827761d-04, 1.4844586873981184d-05,
     &     6.8283193308711422d-07, 1.8810248410796518d-08,
     &     2.8623502429738514d-10, 2.1270790332241105d-12,
     &     6.2979670025179594d-15, 5.0504737000353956d-18,
     &     4.1614623703728548d-22, 1.6874680185111446d-01,
     &     2.9125436200606764d-01, 2.6668610286700062d-01,
     &     1.6600245326950708d-01, 7.4826064668792408d-02,
     &     2.4964417309283247d-02, 6.2025508445722223d-03,
     &     1.1449623864769028d-03, 1.5574177302781227d-04,
     &     1.5401440865224898d-05, 1.0864863665179799d-06,
     &     5.3301209095567054d-08, 1.7579811790505857d-09,
     &     3.7255024025122967d-11, 4.7675292515782048d-13,
     &     3.3728442433624315d-15, 1.1550143395004071d-17,
     &     1.5395221405823110d-20, 5.2864427255691140d-24,
     &     1.6564566124989991d-28, 1.4281197333478154d-01,
     &     2.5877410751742391d-01, 2.5880670727286992d-01,
     &     1.8332268897777793d-01, 9.8166272629918963d-02,
     &     4.0732478151408603d-02, 1.3226019405120104d-02,
     &     3.3693490584783083d-03, 6.7216256409355021d-04,
     &     1.0446121465927488d-04, 1.2544721977993268d-05,
     &     1.1513158127372857d-06, 7.9608129591336357d-08,
     &     4.0728589875500037d-09, 1.5070082262925912d-10,
     &     3.9177365150584634d-12, 6.8941810529581520d-14,
     &     7.8198003824593093d-16, 5.3501888130099474d-18,
     &     2.0105174645555229d-20, 3.6057658645531092d-23,
     &     2.4518188458785009d-26, 4.0883015936805334d-30,
     &     5.5753457883284229d-35 /))      
      PARAMETER (LaBP0=(/   5.8578643762690485d-1, 
     *     3.4142135623730949d+00,  4.1577455678347897d-1, 
     *     2.2942803602790409d0,    6.2899450829374803d0, 
     *     3.2254768961939217d-1,  1.7457611011583465d0, 
     *     4.5366202969211287d0,    9.3950709123011364d0, 
     *     2.6356031971814076d-1,  1.4134030591065161d0, 
     *     3.5964257710407206d0,    7.0858100058588356d0, 
     *     1.2640800844275784d+01,  2.2284660417926061d-1, 
     *     1.1889321016726229d0,    2.9927363260593141d+00, 
     *     5.7751435691045128d0,    9.8374674183825839d0, 
     *     1.5982873980601699d+01,  1.9304367656036231d-1, 
     *     1.0266648953391919d0,    2.5678767449507460d0, 
     *     4.9003530845264844d0,    8.1821534445628572d0, 
     *     1.2734180291797809d+01,  1.9395727862262543d+01, 
     *     1.7027963230510107d-1,  9.0370177679938035d-1, 
     *     2.2510866298661316d0,    4.2667001702876597d0, 
     *     7.0459054023934673d0,    1.0758516010180994d+01, 
     *     1.5740678641278004d+01,  2.2863131736889272d+01, 
     *     1.5232222773180798d-1,  8.0722002274225590d-1, 
     *     2.0051351556193473d0,    3.7834739733312328d0, 
     *     6.2049567778766175d0,    9.3729852516875773d0, 
     *     1.3466236911092089d+01,  1.8833597788991703d+01, 
     *     2.6374071890927389d+01,  1.3779347054049221d-1, 
     *     7.2945454950317090d-1,  1.8083429017403163d0, 
     *     3.4014336978548996d0, 
     *     5.5524961400638029d0,    8.3301527467644991d0, 
     *     1.1843785837900066d+01,  1.6279257831378107d+01, 
     *     2.1996585811980765d+01,  2.9920697012273894d+01 ,
     &     1.1572211735802050d-01, 6.1175748451513112d-01,
     &     1.5126102697764183d+00, 2.8337513377435077d+00,
     &     4.5992276394183476d+00, 6.8445254531151809d+00,
     &     9.6213168424568707d+00, 1.3006054993306348d+01,
     &     1.7116855187462260d+01, 2.2151090379397019d+01,
     &     2.8487967250983996d+01, 3.7099121044466933d+01,
     &     8.7649410478926978d-02, 4.6269632891508106d-01,
     &     1.1410577748312269d+00, 2.1292836450983796d+00,
     &     3.4370866338932058d+00, 5.0780186145497677d+00,
     &     7.0703385350482320d+00, 9.4383143363919331d+00,
     &     1.2214223368866158d+01, 1.5441527368781616d+01,
     &     1.9180156856753147d+01, 2.3515905693991915d+01,
     &     2.8578729742882153d+01,
     &     3.4583398702286622d+01, 4.1940452647688396d+01,
     &     5.1701160339543350d+01, 7.0539889691989419d-02,
     &     3.7212681800161185d-01, 9.1658210248327376d-01,
     &     1.7073065310283420d+00, 2.7491992553094309d+00,
     &     4.0489253138508827d+00, 5.6151749708616148d+00,
     &     7.4590174536710663d+00, 9.5943928695810943d+00,
     &     1.2038802546964314d+01, 1.4814293442630738d+01,
     &     1.7948895520519383d+01, 2.1478788240285009d+01,
     &     2.5451702793186907d+01, 2.9932554631700611d+01,
     &     3.5013434240478986d+01, 4.0833057056728535d+01,
     &     4.7619994047346523d+01, 5.5810795750063903d+01,
     &     6.6524416525615763d+01, 5.9019852181507730d-02,
     &     3.1123914619848325d-01, 7.6609690554593646d-01,
     &     1.4255975908036129d+00, 2.2925620586321909d+00,
     &     3.3707742642089964d+00, 4.6650837034671726d+00,
     &     6.1815351187367655d+00, 7.9275392471721489d+00,
     &     9.9120980150777047d+00, 1.2146102711729766d+01,
     &     1.4642732289596671d+01, 1.7417992646508978d+01,
     &     2.0491460082616424d+01, 2.3887329848169724d+01,
     &     2.7635937174332710d+01, 3.1776041352374712d+01,
     &     3.6358405801651635d+01, 4.1451720484870783d+01,
     &     4.7153106445156347d+01, 5.3608574544695017d+01,
     &     6.1058531447218698d+01, 6.9962240035105026d+01,
     &     8.1498279233948850d+01/))

                                !Laguerre nodes for alpha=-0.5
      PARAMETER (LaBP5 = (/ 2.7525512860841095e-01,
     &      2.7247448713915889e+00, 1.9016350919348812e-01,
     &      1.7844927485432514e+00, 5.5253437422632619e+00,
     &      1.4530352150331699e-01, 1.3390972881263605e+00,
     &      3.9269635013582880e+00, 8.5886356890120332e+00,
     &      1.1758132021177792e-01, 1.0745620124369035e+00,
     &      3.0859374437175511e+00, 6.4147297336620337e+00,
     &      1.1807189489971735e+01, 9.8747014068480951e-02,
     &      8.9830283456961701e-01, 2.5525898026681721e+00,
     &      5.1961525300544675e+00, 9.1242480375311814e+00,
     &      1.5129959781108084e+01, 8.5115442997593743e-02,
     &      7.7213792004277715e-01, 2.1805918884504596e+00,
     &      4.3897928867310174e+00, 7.5540913261017897e+00,
     &      1.1989993039823887e+01, 1.8528277495852500e+01,
     &      7.4791882596818141e-02, 6.7724908764928937e-01,
     &      1.9051136350314275e+00, 3.8094763614849056e+00,
     &      6.4831454286271679e+00, 1.0093323675221344e+01,
     &      1.4972627088426393e+01, 2.1984272840962646e+01,
     &      6.6702230958194261e-02, 6.0323635708174905e-01,
     &      1.6923950797931777e+00, 3.3691762702432655e+00,
     &      5.6944233429577471e+00, 8.7697567302685968e+00,
     &      1.2771825354869195e+01, 1.8046505467728977e+01,
     &      2.5485979166099078e+01, 6.0192063149587700e-02,
     &      5.4386750029464592e-01, 1.5229441054044432e+00,
     &      3.0225133764515753e+00, 5.0849077500985240e+00,
     &      7.7774392315254426e+00, 1.1208130204348663e+01,
     &      1.5561163332189356e+01, 2.1193892096301536e+01,
     &      2.9024950340236231e+01, 5.0361889117293709e-02,
     &      4.5450668156378027e-01, 1.2695899401039612e+00,
     &      2.5098480972321284e+00, 4.1984156448784127e+00,
     &      6.3699753880306362e+00, 9.0754342309612088e+00,
     &      1.2390447963809477e+01, 1.6432195087675318e+01,
     &      2.1396755936166095e+01, 2.7661108779846099e+01,
     &      3.6191360360615583e+01, 3.7962914575312985e-02,
     &      3.4220015601094805e-01, 9.5355315539086472e-01,
     &      1.8779315076960728e+00, 3.1246010507021431e+00,
     &      4.7067267076675874e+00, 6.6422151797414388e+00,
     &      8.9550013377233881e+00, 1.1677033673975952e+01,
     &      1.4851431341801243e+01, 1.8537743178606682e+01,
     &      2.2821300693525199e+01, 2.7831438211328681e+01,
     &      3.3781970488226136e+01, 4.1081666525491165e+01,
     &      5.0777223877537075e+01, 3.0463239279482423e-02,
     &      2.7444471579285024e-01, 7.6388755844391365e-01,
     &      1.5018014976681033e+00, 2.4928301451213657e+00,
     &      3.7434180412162927e+00, 5.2620558537883513e+00,
     &      7.0596277357415627e+00, 9.1498983120306470e+00,
     &      1.1550198286442805e+01, 1.4282403685210406e+01,
     &      1.7374366975199074e+01, 2.0862075185437845e+01,
     &      2.4793039892463458e+01, 2.9231910157093431e+01,
     &      3.4270428925039589e+01, 4.0046815790245596e+01,
     &      4.6788846392124952e+01, 5.4931555621020564e+01,
     &      6.5589931990639684e+01, 2.5437996585689085e-02,
     &      2.2910231649262403e-01, 6.3729027873266897e-01,
     &      1.2517406323627462e+00, 2.0751129098523808e+00,
     &      3.1110524551477146e+00, 4.3642830769353065e+00,
     &      5.8407332713236055e+00, 7.5477046800234531e+00,
     &      9.4940953300264859e+00, 1.1690695926056069e+01,
     &      1.4150586187285759e+01, 1.6889671928527100e+01,
     &      1.9927425875242456e+01, 2.3287932824879903e+01,
     &      2.7001406056472355e+01, 3.1106464709046559e+01,
     &      3.5653703516328221e+01, 4.0711598185543110e+01,
     &      4.6376979557540103e+01, 5.2795432527283602e+01,
     &      6.0206666963057259e+01, 6.9068601975304347e+01,
     &      8.0556280819950416e+01/))
      
      PARAMETER (LaWF5 = (/   1.6098281800110255e+00,
     &      1.6262567089449037e-01, 1.4492591904487846e+00,
     &      3.1413464064571323e-01, 9.0600198110176913e-03,
     &      1.3222940251164819e+00, 4.1560465162978422e-01,
     &      3.4155966014826969e-02, 3.9920814442273529e-04,
     &      1.2217252674706509e+00, 4.8027722216462992e-01,
     &      6.7748788910962143e-02, 2.6872914935624635e-03,
     &      1.5280865710465251e-05, 1.1402704725249586e+00,
     &      5.2098462052832328e-01, 1.0321597123176789e-01,
     &      7.8107811692581406e-03, 1.7147374087175731e-04,
     &      5.3171033687126004e-07, 1.0728118194241802e+00,
     &      5.4621121812849427e-01, 1.3701106844693015e-01,
     &      1.5700109452915889e-02, 7.1018522710384658e-04,
     &      9.4329687100378043e-06, 1.7257182336250307e-08,
     &      1.0158589580332265e+00, 5.6129491705706813e-01,
     &      1.6762008279797133e-01, 2.5760623071019968e-02,
     &      1.8645680172483614e-03, 5.4237201850757696e-05,
     &      4.6419616897304271e-07, 5.3096149480223697e-10,
     &      9.6699138945091101e-01, 5.6961457133995952e-01,
     &      1.9460349528263074e-01, 3.7280084775089407e-02,
     &      3.7770452605368474e-03, 1.8362253735858719e-04,
     &      3.6213089621868382e-06, 2.0934411591584102e-08,
     &      1.5656399544231742e-11, 9.2448733920121973e-01,
     &      5.7335101072566907e-01, 2.1803441204004675e-01,
     &      4.9621041774927162e-02, 6.4875466844757246e-03,
     &      4.5667727203270848e-04, 1.5605112957064066e-05,
     &      2.1721387415385585e-07, 8.7986819845463701e-10,
     &      4.4587872910682818e-13, 8.5386232773739834e-01,
     &      5.7235907069288550e-01, 2.5547924356911883e-01,
     &      7.4890941006461639e-02, 1.4096711620145414e-02,
     &      1.6473849653768340e-03, 1.1377383272808749e-04,
     &      4.3164914098046565e-06, 8.0379423498828602e-08,
     &      6.0925085399751771e-10, 1.3169240486156312e-12,
     &      3.3287369929782692e-16, 7.5047670518560539e-01,
     &      5.5491628460505815e-01, 3.0253946815328553e-01,
     &      1.2091626191182542e-01, 3.5106857663146820e-02,
     &      7.3097806533088429e-03, 1.0725367310559510e-03,
     &      1.0833168123639965e-04, 7.3011702591247581e-06,
     &      3.1483355850911864e-07, 8.1976643295418016e-09,
     &      1.1866582926793190e-10, 8.4300204226528705e-13,
     &      2.3946880341857530e-15, 1.8463473073036743e-18,
     &      1.4621352854768128e-22, 6.7728655485117817e-01,
     &      5.3145650375475362e-01, 3.2675746542654360e-01,
     &      1.5694921173080897e-01, 5.8625131072344717e-02,
     &      1.6921776016516312e-02, 3.7429936591959084e-03,
     &      6.2770718908266166e-04, 7.8738679621849850e-05,
     &      7.2631523013860402e-06, 4.8222883273410492e-07,
     &      2.2424721664551585e-08, 7.0512415827308280e-10,
     &      1.4313056105380569e-11, 1.7611415290432366e-13,
     &      1.2016717578981511e-15, 3.9783620242330409e-18,
     &      5.1351867308233644e-21, 1.7088113927550770e-24,
     &      5.1820874276942667e-29, 6.2200206075592535e-01,
     &      5.0792308532951769e-01, 3.3840894389128295e-01,
     &      1.8364459415856996e-01, 8.0959353969207851e-02,
     &      2.8889923149962169e-02, 8.3060098239550965e-03,
     &      1.9127846396388331e-03, 3.5030086360234562e-04,
     &      5.0571980554969836e-05, 5.6945173834697106e-06,
     &      4.9373179873395243e-07, 3.2450282717915824e-08,
     &      1.5860934990330932e-09, 5.6305930756763865e-11,
     &      1.4093865163091798e-12, 2.3951797309583852e-14,
     &      2.6303192453168292e-16, 1.7460319202373756e-18,
     &      6.3767746470103704e-21, 1.1129154937804721e-23,
     &      7.3700721603011131e-27, 1.1969225386627985e-30,
     &      1.5871102921547987e-35 /))
      END MODULE QUAD

!*****************************************************

     
      MODULE RIND

      IMPLICIT NONE
      PRIVATE
      PUBLIC :: RINDD, INITDATA, SETDATA,ECHO
 
      INTERFACE RINDD
      MODULE PROCEDURE RINDD
      END INTERFACE
      
      INTERFACE FI
      MODULE PROCEDURE FI
      END INTERFACE 
     
      INTERFACE FIINV
      MODULE PROCEDURE FIINV
      END INTERFACE
! alternative FIINV slightly faster but less accurate
      INTERFACE FIINVS 
      MODULE PROCEDURE FIINVS
      END INTERFACE

      INTERFACE MNORMPRB
      MODULE PROCEDURE MNORMPRB
      END INTERFACE

      INTERFACE THL
      MODULE PROCEDURE THL
      END INTERFACE 

      INTERFACE NORM2DPRB
      MODULE PROCEDURE NORM2DPRB
      END INTERFACE 

      INTERFACE SETDATA
      MODULE PROCEDURE SETDATA
      END INTERFACE

      INTERFACE INITDATA
      MODULE PROCEDURE INITDATA
      END INTERFACE

      INTERFACE JACOB
      MODULE PROCEDURE JACOB
      END INTERFACE 
      
      INTERFACE ARGP0
      MODULE PROCEDURE ARGP0
      END INTERFACE
            
      INTERFACE  RINDDND
      MODULE PROCEDURE RINDDND
      END INTERFACE

      INTERFACE  RINDSCIS
      MODULE PROCEDURE RINDSCIS
      END INTERFACE

      INTERFACE  MVNDUN
      MODULE PROCEDURE MVNFUN
      END INTERFACE

      INTERFACE  MVNDUN2
      MODULE PROCEDURE MVNFUN2
      END INTERFACE

      INTERFACE RINDNIT
      MODULE  PROCEDURE RINDNIT
      END INTERFACE
      
      INTERFACE  BARRIER  
      MODULE PROCEDURE BARRIER
      END INTERFACE

      INTERFACE GAUSINT
      MODULE PROCEDURE GAUSINT
      END INTERFACE
      
      INTERFACE GAUSINT2
      MODULE PROCEDURE GAUSINT2
      END INTERFACE


      INTERFACE GAUSSLA0
      MODULE PROCEDURE GAUSSLA0
      END INTERFACE

      INTERFACE GAUSSLE0
      MODULE PROCEDURE GAUSSLE0
      END INTERFACE

      INTERFACE GAUSSHE0
      MODULE PROCEDURE GAUSSHE0
      END INTERFACE
      

      INTERFACE GAUSSLE1
      MODULE PROCEDURE GAUSSLE1 
      END INTERFACE

      INTERFACE GAUSSLE2
      MODULE PROCEDURE GAUSSLE2 
      END INTERFACE

      INTERFACE GAUSSQ
      MODULE PROCEDURE GAUSSQ 
      END INTERFACE

      INTERFACE C1C2
      MODULE PROCEDURE C1C2 
      END INTERFACE

      INTERFACE echo 
      MODULE PROCEDURE echo
      END INTERFACE

      INTERFACE swapRe
      MODULE PROCEDURE swapRe
      END INTERFACE

      INTERFACE swapint
      MODULE PROCEDURE swapint
      END INTERFACE

      INTERFACE getdiag
      MODULE PROCEDURE getdiag
      END INTERFACE


      INTERFACE spearcorr
      MODULE PROCEDURE spearcorr
      END INTERFACE
 
      INTERFACE  SVDCMP
      MODULE PROCEDURE SVDCMP
      END INTERFACE

      INTERFACE  pythag
      MODULE PROCEDURE pythag
      END INTERFACE
   
      INTERFACE CONDSORT
      MODULE PROCEDURE CONDSORT
      END INTERFACE

   
      INTERFACE CONDSORT2
      MODULE PROCEDURE CONDSORT2
      END INTERFACE

      INTERFACE CONDSORT3
      MODULE PROCEDURE CONDSORT3
      END INTERFACE

      INTERFACE SORTRE
      MODULE PROCEDURE SORTRE
      END INTERFACE

      INTERFACE SOBSEQ
      MODULE PROCEDURE SOBSEQ
      END INTERFACE
                                !--------------------------------
      CONTAINS
      FUNCTION JACOB ( xd,xc) RESULT (value1)                                
      USE GLOBALDATA , ONLY : jacobdef
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(:),INTENT(in) :: xd ,xc
      DOUBLE PRECISION :: value1
      SELECT CASE (jacobdef)
      ! if you add another jacobian
      ! then rindd1 may be modified for NIT=0
      ! or Nstoc<1 to increase speed
      CASE (0)                  ! default
         value1 =ABS(PRODUCT(xd))                             
      CASE (1)
         value1=1.d0
      case (2) 
         value1 =ABS(PRODUCT(xd)*product(xc))  
      END SELECT
                           
      RETURN                                                             
      END FUNCTION JACOB

      SUBROUTINE SETDATA (dEPSS,dEPS2,dXc,dNINT) 
      USE GLOBALDATA
      USE QUAD, ONLY: sizNint,Nint1
      IMPLICIT NONE
      DOUBLE PRECISION , INTENT(in) :: dEPSs,dEPS2,dXc 
      DOUBLE PRECISION , DIMENSION(:), INTENT(in) :: dNINT
      INTEGER :: N
      
      N=SIZE(dNINT)
      IF (sizNint.LT.N) THEN
         PRINT *,'Error in setdata, Nint too large'
         N=sizNint
      ENDIF
      NINT1(1:N)=dNINT(1:N)  ! quadrature formulae for the Xd variables
      IF (N.LT.sizNint) THEN
         NINT1(N:sizNint)=NINT1(N)  
      END IF                        
  
      EPS2     = dEPS2          ! Constants controlling 
      EPSS     = dEPSS          ! accuracy of integration 
      CEPSS    = 1.d0 - EPSS 
      xCutOff  = dXc
      xCutOff2 = dXc*dXc        !The region we are integrating over.
      FxCutOff = FI(xCutOff) - FI( -xCutOff) 
      CFxCutOff=1.d0-FxCutOff                              
      XSMALL=EPSS
      XMAX=xCutOff
      RETURN                    ! INITINTEGRAL                         
      END SUBROUTINE SETDATA
      
      SUBROUTINE INITDATA (speed)                                   
      USE GLOBALDATA
      USE QUAD, ONLY: sizNint,Nint1,minQnr,Le2Qnr
      IMPLICIT NONE
      INTEGER , INTENT(in) :: speed 
      
           
      SELECT CASE (speed)
      CASE (9:)
         
         Le2Qnr=1
         NINT1 (1) = 1                                                   
         NINT1 (2) = 1                                                   
         NINT1 (3) = 1                                                   
         EPS2 = 1d-3                                                      
         EPSS = 1d-3 
         NUGGET=1d-6                                                  
      CASE (8)
         Le2Qnr=2
         NINT1 (1) = 2                                                   
         NINT1 (2) = 2                                                   
         NINT1 (3) = 2                                                   
         EPS2 = 1d-3                                                      
         EPSS = 1d-4
         NUGGET=1d-6
      CASE (7)
         Le2Qnr=3
         NINT1 (1) = 2                                                   
         NINT1 (2) = 3                                                   
         NINT1 (3) = 4                                                   
         EPS2 = 1d-3                                                      
         EPSS = 1d-4
         NUGGET=1d-6
      CASE (6)
         Le2Qnr=4
         NINT1 (1) = 3                                                   
         NINT1 (2) = 4                                                   
         NINT1 (3) = 5                                                   
         EPS2 = 1d-4                                                      
         EPSS = 1d-4
         NUGGET=1d-6
      CASE (5)
         Le2Qnr=5
         NINT1 (1) = 4                                                   
         NINT1 (2) = 5                                                   
         NINT1 (3) = 6                                                   
         EPS2 = 1d-4                                                      
         EPSS = 1d-4
         NUGGET=1d-6
      CASE (4)     ! quadrature formulae for the Xd variables
         Le2Qnr=6
         NINT1 (1) = 6      ! use quadr. form. No. 6 in integration of Xd(1)
         NINT1 (2) = 7      ! use quadr. form. No. 7 in integration of Xd(2)
         NINT1 (3) = 8      ! use quadr. form. No. 8 in integration of Xd(3)
                                
         EPS2 = 1d-4        ! Constants controlling 
         EPSS = 1d-5        ! accuracy of integration
         NUGGET=1d-6
      CASE (3)
         Le2Qnr=7
         NINT1 (1) = 8                                                   
         NINT1 (2) = 9                                                   
         NINT1 (3) = 10                                                   
         EPS2 = 1d-4                                                      
         EPSS = 1d-5 
         NUGGET=1d-6
      CASE (2)
         Le2Qnr=8
         NINT1 (1) = 9                                                   
         NINT1 (2) = 10                                                  
         NINT1 (3) = 11                                                  
         EPS2 = 1d-5                                                      
         EPSS = 1d-6
         NUGGET=1d-7
      CASE (:1)
         Le2Qnr=9
         NINT1 (1) = 11                                                  
         NINT1 (2) = 12                                                  
         NINT1 (3) = 13                                                  
         EPS2 = 1d-6                                                    
         EPSS = 1d-7
         NUGGET=1d-8
      END SELECT

      NINT1(4:sizNint)=NINT1(3) 
      minQnr=1                  ! minimum qudrature No. used in GaussLe1,Gaussle2
                                !The region we are integrating over.
      FxCutOff = FI(xCutOff) - FI( -xCutOff)  
      xCutOff2=xCutOff*xCutOff
      CFxCutOff=1.d0-FxCutOff 
                                                       
      
      if (SCIS.gt.0) then  
         !NUGGET=1d-8
         C1C2det=.true.
         EPSS=EPSS*1.d2
      else
         C1C2det=.true.
      endif
      CEPSS = 1.d0 - EPSS 
      XSMALL=ABS(FIINV(.5d0+CFxCutOff/2.d0)) 
      XMAX=xCutOff
      RETURN                                                             
                                ! INITINTEGRAL                         
      END SUBROUTINE INITDATA

      SUBROUTINE ECHO(array)
      INTEGER ::j
      DOUBLE PRECISION,DIMENSION(:,:)::array
      DO j=1,size(array,1)
         PRINT 111,j,array(j,:)
111      FORMAT (i2,':',10F10.5)
      END DO
      END SUBROUTINE ECHO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!******************* RINDD - the main program *********************!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE RINDD (fxind,Big, Ex, xc,indI,Blo,Bup)  
      USE GLOBALDATA, ONLY :Nt,Nj,Njj,Nd,Nc,Nx,Ntd,Ntdc,NsXtmj,NsXdj,
     &     indXtd,index1,xedni,SQ,Hlo,Hup,fxcepss,EPS2,XCEPS2,NIT,
     &     SQTWOPI1,xCutOff2,SCIS,Ntscis,COVix
      IMPLICIT NONE     
      DOUBLE PRECISION, DIMENSION(:,:), INTENT(inout) :: BIG
      DOUBLE PRECISION, DIMENSION(:,:), INTENT(in)    ::  xc 
      DOUBLE PRECISION, DIMENSION(:  ), INTENT(in)    ::  Ex            
      DOUBLE PRECISION, DIMENSION(:,:), INTENT(out)   ::  fxind 
      DOUBLE PRECISION, DIMENSION(:,:), INTENT(in)    :: Blo, Bup  
      INTEGER,          DIMENSION(:  ), INTENT(in)    :: indI    
      DOUBLE PRECISION, DIMENSION(:  ), ALLOCATABLE   :: xd, Cm 
      INTEGER                                     :: J, ix,Ntdcmj,Nstoc
      DOUBLE PRECISION, DIMENSION(2)                  :: xind
      DOUBLE PRECISION                                :: SS0,xx
      DOUBLE PRECISION                                :: fxc,quant,Sqeps
 
                                !Initialization
                                !Call Initdata(speed)
      Nj=MIN(Nj,MAX(Nt,0))      ! make sure Nj<=Nt
      Njj=MIN(Njj,MAX(Nt-Nj,0)) ! make sure Njj<=Nt-Nj

      IF (Nd.GT.0) THEN
         ALLOCATE(xd(1:Nd))
         xd=0.d0 
      END IF

      If (SCIS.GT.0) then
         Ntscis=Nt-Nj-Njj
         ALLOCATE(SQ(1:Ntd,1:Ntd)) ! Cond. stdev's  
         ALLOCATE(NsXtmj(1:Ntd+1)) ! indices to stoch. var. See condsort
      else
         Ntscis=0
         ALLOCATE(SQ(1:Ntd,1:max(Njj+Nj+Nd,1)) ) ! Cond. stdev's  
         ALLOCATE(NsXtmj(1:Nd+Nj+Njj+1)) ! indices to stoch. var. See condsort
      endif

      ALLOCATE(Cm(1:Ntdc))      !Cond. mean which has the same order as local
      Cm=0.d0                   !covariance matrices (after sorting) or excluding
                                !irrelevant variables.
      
      ALLOCATE(index1(1:Ntdc))  ! indices to the var. original place in BIG 
      index1=0                  ! (before sort.)
      ALLOCATE(xedni(1:Ntdc))   ! indices to var. new place (after sorting),  
      xedni=0                   ! eg. the point xedni(1) is the original position
                                ! of variable with conditional mean CM(1).
      ALLOCATE(Hlo(1:Ntd))      ! lower and upper integration limits are computed
                                ! in the new order that is the same as CM.
                                ! This convention is expressed in the vector indXTD.
      Hlo=0.d0                  ! However later on some variables will be exluded
                                ! since those are irrelevant and hence CMnew(1)
                                ! does not to be conditional mean of the same variable
                                ! as CM(1) is from the beginning. Consequently  
      ALLOCATE(Hup(1:Ntd))      ! the order of Hup, Hlo will be unchanged. So we need  
      Hup=0.d0                  ! to know where the relevants variables bounds are
                                ! This will be given in the subroutines by a vector indS.

      ALLOCATE(NsXdj(1:Nd+Nj+1)) ! indices to stoch. var. See condsort
      NsXdj=0
      ALLOCATE(indXtd(1:Ntd))    ! indices to Xt and Xd as they are 
      indXtd=(/(J,J=1,Ntd)/)     ! sorted in Hlo and Hup

                                 !conditional covariance matrix BIG            
      CALL CONDSORT (BIG,SQ, index1, xedni, NsXtmj,NsXdj)
      
      !PRINT *, 'index=', index1
      !PRINT *, 'xedni=', xedni
!       PRINT *,BIG(1,1),BIG(2,2),BIG(1,2)

      fxind  = 0.d0             ! initialize
                                ! Now the loop over all different values of
                                ! variables Xc (the one one is conditioning on)
      DO  ix = 1, Nx            ! is started. The density f_{Xc}(xc(:,ix))                            
         COVix = ix             ! will be computed and denoted by  fxc.
         xind = 0.d0                                                       
         fxc = 1.d0                                                         

                                 ! Set the original means of the variables
         Cm  =Ex (index1(1:Ntdc))!   Cm(1:Ntdc)  =Ex (index1(1:Ntdc))
         DO J = 1, Nc            !Recursive conditioning on the last Nc variables  
            Ntdcmj=Ntdc-J
            SS0=BIG(Ntdcmj+1,Ntdcmj+1) ! var(X(i)|X(i+1),X(i+2),...,X(Ntdc))
                                       ! i=Ntdc-J+1 (J=1 var(X(Ntdc))

            
            IF (SS0.LT.XCEPS2) THEN !Degenerated case the density 
                                    !can not be computed
               PRINT *,'Degenerate case of Xc(Nc-J+1) for J=',j 

               GOTO 110         ! degenerate case exit fxind=0 for all  
            ENDIF               ! (should perhaps return NaN instead??)

            xx = xc(index1(Ntdcmj+1)-Ntd,ix)-Cm(Ntdcmj+1)
            quant = 0.5d0 * xx * xx / SS0                                       
            IF (quant.GT.90.d0) GOTO 100                                                     
            
                                ! fxc probability density for i=Ntdc-J+1, 
                                ! fXc=f(X(i)|X(i+1),X(i+2)...X(Ntdc))*
                                !     f(X(i+1)|X(i+2)...X(Ntdc))*..*f(X(Ntdc))
            fxc = fxc*SQTWOPI1*EXP(-quant)/SQRT(SS0)   
                                ! conditional mean (expectation) 
                                ! E(X(1:i-1)|X(i),X(i+1),...,X(Ntdc)) 
            Cm(1:Ntdcmj) = Cm(1:Ntdcmj)+xx*BIG (1:Ntdcmj,Ntdcmj+1)/SS0
           
         ENDDO                  ! J
         !PRINT *, 'Rindd, Cm=',Cm(xedni(max(1,Nt-5):Ntdc))

         IF (fxc .LT.fxcEpss) THEN
            GOTO 100      ! Small probability dont bother calculating it
         ENDIF

                          !set the global integration limits Hlo,Hup 
         CALL BARRIER(xc(:,ix),indI,Blo,Bup)
         xind(1)=1.0d0
         xind(2)=0.0d0


        sqeps=sqrt(eps2)
        Nstoc= NsXtmj(Ntscis+Njj+Nd+Nj+1)                                                 
        IF (any((Cm(Nstoc+1:Nt-Nj) .GT.Hup(Nstoc+1:Nt-Nj)+sqeps ).OR. 
     *     (Cm (Nstoc+1:Nt-Nj) .LT.Hlo (Nstoc+1:Nt-Nj)-sqeps))) THEN  !degenerate case 
           xind(1)=0.d0              !mean of deterministic variable is 
           goto 100                 ! outside the barriers     
        ENDIF
        
      IF (Nstoc.LT.1.and.SCIS.EQ.0) then
        print *,'Warning nr. of var. in idicator = 0 - indicator is 
     &replaced by 1. or 0.'       
      end if

         IF (SCIS.EQ.2) then    ! integrate all by SCIS
            XIND=RINDSCIS(BIG,Cm,xd,xc(:,ix))
            goto 100
         endif

         SELECT CASE (Nd+Nj)
         CASE (:0) 
           IF (SCIS.NE.0) then       ! integrate all by SCIS
             XIND=MNORMPRB(BIG,Cm(1:Nstoc))
           ELSE 
            XIND=RINDNIT(BIG,SQ(1:Nstoc,1),Cm,indXtd(1:Nstoc),NIT)
           END IF          
         CASE (1:)
            xind=RINDDND(BIG,Cm,xd,xc(:,ix),Nd,Nj) 
         END SELECT
 100     continue
         fxind(ix,1)=xind(1)*fxc
         fxind(ix,2)=xind(2)*fxc                                          
      ENDDO                     !ix

 110  IF (ALLOCATED(xd)) THEN
         DEALLOCATE(xd) 
      END IF
     
      DEALLOCATE(SQ)
      DEALLOCATE(Hlo)
      DEALLOCATE(Hup)  
      DEALLOCATE(Cm) 
      DEALLOCATE(index1)
      DEALLOCATE(xedni) 
      DEALLOCATE(NsXtmj)
      DEALLOCATE(NsXdj)
      DEALLOCATE(indXtd) 
      RETURN                                                            
      END SUBROUTINE RINDD

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!*************************** ARGP0 *********************************!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE ARGP0 (I0,I1,P0,Plo,SQ,Cm,indS,ind,Nind)
      USE GLOBALDATA, ONLY : Hlo,Hup,xCutOff,CFxCutOff,EPSS,EPS2
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(:), INTENT(in)  :: SQ , Cm !stdev./mean
      INTEGER,          DIMENSION(:), INTENT(in)  :: indS
      INTEGER,          DIMENSION(:), INTENT(out) :: ind
      DOUBLE PRECISION,               INTENT(out) :: P0,Plo
      INTEGER,                        INTENT(out) :: I0, I1
      INTEGER,                        INTENT(out) :: Nind
      DOUBLE PRECISION                            :: PPROD,P1
      DOUBLE PRECISION                            :: Xup, XLo
      INTEGER                                     :: I, Nstoc
                                ! indS contains the indices to the limits
      Nstoc=SIZE(indS)          ! in Hlo/Hup of variables in the indicator
                                ! ind contains indeces to the relevant
                                ! variables which are Nind<=Nstoc.
                                ! We wish to compute P(Hlo<X<Hup) but
                                ! only have lower and upper bounds Plo,P0, resp.
                                ! I0 is the position of the minimal
                                ! probability in the vector ind, i.e.
                                ! P0=P(Hlo<X(indS(ind(I0)))<Hup)
                                ! I1 is the second minimum.
      P0=2.d0
      I0=1
      I1=1
      Plo=0.d0
      Nind=0
      DO I = 1, Nstoc        
         Xup=xCutOff
         Xlo=-xCutOff
         IF (SQ(I).GE.EPS2) THEN
            Xup = MIN( (Hup (indS(I)) - Cm (I))/ SQ(I),Xup)
            Xlo = MAX( (Hlo (indS(I)) - Cm (I))/ SQ(I),Xlo)
          ELSE 
            IF (Hup(indS(I))-Cm (I).lt.0.d0) Xup=Xlo
            IF (Hlo(indS(I))-Cm (I).gt.0.d0) Xlo=Xup
         END IF 
         IF (Xup.LE.Xlo+EPSS) THEN  ! +EPSS
            P0=0.d0
            Plo=0.d0
            I0=I
            Nind=1
            RETURN                                             
         ENDIF

       IF ((Xup.LT.xCutOff-epss).or.(Xlo.gt.epss-xCutOff)) THEN
         Nind=Nind+1
         ind(Nind)=i
                                ! this procedure calculates 
         P1 =FI(Xup)-FI(Xlo)
         Plo=Plo+P1
         IF (P1.LT.P0) THEN   
            I1 = I0
            I0 = Nind              ! Prob(I0)=Prob(XMA>X(i0)>XMI)= 
            P0 = P1        !     min Prob(Hup(i)> X(i)>Hlo(i))
            IF (P0.LT.EPSS) THEN
               Plo=0.d0;
               RETURN                           
            ENDIF
         ENDIF
       ENDIF 
      ENDDO

      Plo=max(0.d0,1.d0-DBLE(Nind)+Plo)
      P0=min(1.d0,P0)

      RETURN                                                             
      END SUBROUTINE ARGP0
 
       

                                 !Ntmj is the number of elements in indicator
                                 !since Nj points of process valeus (Nt) have
                                 !been moved to the jacobian.
                                 !index1 contains the original 
                                 !positions of variables in the
                                 !covaraince matrix before condsort
                                 !and that why if index(Ntmj+1)>Nt
                                 !it means the variable to conditon on
                                 !is a derivative isXd=1

                                 != # stochastic variables before
                                 !conditioning on X(Ntmj+1). This
                                 !I still not checked why.
                 


! ******************* RINDDND ****************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
      RECURSIVE FUNCTION RINDDND (BIG,Cm,xd,xc,Ndleft,Njleft) 
     &      RESULT (xind) 
      USE GLOBALDATA, ONLY :SQPI1, SQTWOPI1,Hup,Hlo,Nt,Nj,Njj,Nd,
     &     NsXtmj,NsXdj,EPS2,NIT,xCutOff,CFxCutOff,EPSS,CEPSS,index1,
     &     indXtd,SQ,SQTWO,SQTWO1,SCIS,Ntscis,C1C2det
      USE QUAD , ONLY : PMAX,Nint1,sizNint
      IMPLICIT NONE
      INTEGER,INTENT(in) :: Ndleft,Njleft ! # DIMENSIONs to integrate
      DOUBLE PRECISION, DIMENSION(:,:), INTENT(inout) :: BIG
      DOUBLE PRECISION, DIMENSION(:  ), INTENT(in) :: Cm ! conditional mean
      DOUBLE PRECISION, DIMENSION(:  ), INTENT(inout) :: xd ! integr. variables
      DOUBLE PRECISION, DIMENSION(:  ), INTENT(in)    :: xc ! conditional values
!local variables
      DOUBLE PRECISION, DIMENSION(2)                  :: xind
      DOUBLE PRECISION, DIMENSION(2)                  :: xind1
      DOUBLE PRECISION, DIMENSION(PMAX)               :: WXdi, Xdi !weights/nodes 
      DOUBLE PRECISION, DIMENSION(:  ), ALLOCATABLE   :: CmNEW
      INTEGER  :: Nrr, Nr, J, I0, N,Ndleft1,Ndjleft,Ntmj,isXd
      INTEGER :: Nst,Nsd,  NstN,NsdN 
      DOUBLE PRECISION     :: SS0,SQ0,fxd,XMA,XMI,tmp,sqeps

      Ntmj=Nt-Nj
      Ndjleft= Ndleft+Njleft
      N=Ntmj+Ndjleft
      
      IF (index1(N).GT.Nt) THEN
         isXd=1
      ELSE
         isXd=0
      END IF
      XIND=0.d0 
      SS0 = BIG (N, N)              
      Nst     = NsXtmj(Ntscis+Njj+Ndjleft+1)                                

!********************************************************************************
!** Here Starts the degenerated case the remaining variables are deterministic **
!********************************************************************************

      IF (SS0.LT.EPS2) THEN 
            sqeps=sqrt(EPS2)   
                                    !Next is the check for the special situation
                                    !that after conditioning on Xc all derivatives are
                                    !singular and not satisfying the limitations 
                                    !(so something is generally wrong) 
            IF (any((Cm(Nst+1:N).GT.Hup(Nst+1:N) +sqeps).OR.
     &           (Cm(Nst+1:N).LT.Hlo(Nst+1:N)-sqeps ))) THEN
               RETURN               !the mean of Xd or Xt is too extreme          
            ENDIF
                                    !Here we are putting in all conditional expectations
                                    !for the values of the deterministics derivatives.
         IF (Nd.GT.0) THEN
            Ndleft1=Ndleft
            DO WHILE (Ndleft1.GT.0)
               IF (index1(N).GT.Nt) THEN  ! isXd
                  xd (Ndleft1) =  Cm (N)
                  Ndleft1=Ndleft1-1                  
               END IF
               N=N-1
            ENDDO
            fxd = jacob (xd,xc) ! jacobian of xd,xc
         ELSE
            fxd = 1.d0 !     XIND = FxCutOff???
         END IF

         XIND=fxd
         XIND(2)=0.d0 
         IF ((Nst.gt.0).and.(SCIS.ne.0)) then
             XIND=fxd*MNORMPRB(BIG,Cm(1:Nst))
         END IF 
         IF ((Nst.gt.0).and.(SCIS.eq.0)) then
            XIND=fxd*RINDNIT(BIG,SQ(:,Ntscis+Njj+1),
     &         Cm,indXtd(1:Nst),NIT)
         END IF
         RETURN
      ENDIF
      
!***** Here Starts the conditioning on the last variable (nondeterministic) *
!****************************************************************************

      SQ0 = SQ(N,Ntscis+Njj+Ndjleft) !SQRT (SS0)   !        
     
      XMA=MIN((Hup (indXtd(N))-Cm (N))/SQ0, xCutOff)
      XMI=MAX((Hlo (indXtd(N))-Cm (N))/SQ0,-xCutOff)

         ! See if we can narrow down integration range
         
      Nsd  = NsXdj(Ndjleft+1) ! index to first stoch. variable of Xd before conditioning on X(N)
      NstN = NsXtmj(Ntscis+Njj+Ndjleft) ! index to last stoch. variable of Xt after cond. on X(N)
      NsdN = NsXdj(Ndjleft)  ! index to first stoch. variable of Xd after conditioning on X(N)

      CALL C1C2(XMI,XMA,Cm(1:N-1),BIG(1:N-1,N),
     &        SQ(1:N-1,Ntscis+Njj+Ndjleft),SQ0,indXtd(1:N-1))

      IF (XMA.LE.XMI) THEN
         XIND=0.d0
         RETURN
      ENDIF
      Nrr = NINT1 (MIN(Ndjleft,sizNint))     
      Nr=0 ! initialize # of nodes 
      !print *, 'rinddnd Nrr',Nrr
                        !Grid the interval [XMI,XMA] by  GAUSS quadr.
      CALL GAUSSLE2(Nr, WXdi, Xdi,XMI,XMA, Nrr)  
                                !print *, 'Xdi',Xdi
      ALLOCATE(CmNEW(1:N-1)) 
      
      DO J = 1, Nr         
         IF (Wxdi(J).GT.(CFxCutOff)) THEN ! EPSS???
            IF (isXd.EQ.1) THEN
               xd (Ndleft) =  Xdi (J)*SQ0 + Cm (N)
            END IF
            CmNEW (1:(N-1)) = 
     &           Cm(1:(N-1))+Xdi(J)*(BIG(1:(N-1),N)/SQ0)
                                       !  Here we start with the case whent there 
                                       !  some derivatives left to integrate.
            IF  (Ndjleft.GT.1) THEN
              fxd=Wxdi(J)
              XIND1=RINDDND(BIG,CmNEW,xd,xc,
     &              Ndleft-isXd,Njleft-1+isXd)
            ELSE
                                        !  Here all is conditioned on
                                        !  and we wish to compute the
                                        !  conditional probability that 
                                        !  variables in indicator stays between bariers.
              XIND1=1.d0
              IF (Nd.GT.0) THEN         !if there are derivatives we need 
                                        !to compute jacobian jacob(xd,xc) 
                  fxd = WXdi(J)*jacob(xd(1:Nd),xc) 
                ELSE                     !Here there are no derivatives 
                                         !and we assume that jacob(xc)=1
                  fxd = WXdi(J) 
              END IF
            
            IF (NstN.LT.1) GOTO 100       !Here there are no points in indicator
                                          !left to integrate and hence XIND1=1 or 0..
                                           !integrate by Monte Carlo - SCIS           
         IF (SCIS.NE.0) XIND1=MNORMPRB(BIG,CmNEW)
                                         !integrate by quadrature     
         IF (SCIS.EQ.0) XIND1=RINDNIT(BIG,
     &             SQ(:,Ntscis+Njj+1),CmNEW,indXtd(1:NstN),NIT)
         END IF
100     continue                                         
         XIND = XIND+XIND1 * fxd
         END IF
      ENDDO
                           
      DEALLOCATE(CmNEW) 
      RETURN                                                             
      END FUNCTION RINDDND



! ******************* RINDNIT ****************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                                ! old procedure rind2-6  
      RECURSIVE FUNCTION RINDNIT(R,SQ,Cm,indS,NITL) RESULT (xind)
      USE GLOBALDATA, ONLY : Hlo,Hup,EPS2, EPSS,CEPSS,CFxCutOff
     &         ,xCutOff,XSPLT
      USE QUAD, ONLY : PMAX,Le2QNr
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(:,:), INTENT(in)    :: R 
      DOUBLE PRECISION, DIMENSION(:  ), INTENT(in)    :: SQ
      DOUBLE PRECISION, DIMENSION(:  ), INTENT(in)    :: Cm
      DOUBLE PRECISION, DIMENSION(2)                  :: xind  
      INTEGER,          DIMENSION(:  ), INTENT(in)    :: indS
      INTEGER,                          INTENT(in)    :: NITL

      DOUBLE PRECISION, DIMENSION(2)                  :: xind1,xind2  
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE   :: RNEW
      DOUBLE PRECISION, DIMENSION(:  ), ALLOCATABLE   :: B,SQnew
      DOUBLE PRECISION, DIMENSION(:  ), ALLOCATABLE   :: CmNEW
      INTEGER,          DIMENSION(:  ), ALLOCATABLE   :: indSNEW,ind
      INTEGER                                         :: I0,I1

      DOUBLE PRECISION, DIMENSION(PMAX)               :: H1, XX1
      INTEGER, DIMENSION(1)                           :: m1
      DOUBLE PRECISION              :: SQ0,SQ1,SS0,SS1,SS
      DOUBLE PRECISION              :: XMI,XMA,XMI1,XMA1,SGN,P0,Plo,rho
      INTEGER      :: Ns,Nsnew,row,col,r1,r2,J,N1 

! Assumption is that there is at least one variable X in the indicator,
! LNIT nonegative integer.
! If  LNIT=0 or the number of relevant variables is less then 3, the recursion
! stops. It gives exact value if after removing irrelevant variables there
! are maximum 2 variables left in the indicator. The program is not using
! RIND2 function any more. IR. 28 XI 1999 - Indianapolis.
!
! explanation to variables (above):
! R       = cov. matr.
! B       = R(I,I0) I=1:Ns
! SQ      = SQRT(R(I,I)) I=1:Ns
! Cm      = cond. mean
! indS    = indices to the stochastic variables as they are stored in
!           the global variables Hlo and Hup
! Ns      = size of indS =# of variables in indicator before conditioning
! Nsnew   = # of relevant variables in indicator before conditioning
! I0,I1   = indecies to minimum prob. and next minimal, respectively 
! ..NEW   = the var. above after conditioning on X(I0) or used in recursion  
! ind     = temp. variable storing indices

      Ns=SIZE(indS)       !=# stochastic variables before conditioning
!      print *,'NS,NITL=',NS,NITL
      XIND=1.d0
      if (Ns.lt.1) return
      
      ALLOCATE(ind(1:NS))
      CALL ARGP0(I0,I1,P0,Plo,SQ,Cm,indS,ind,NSnew)
                           !The probability of being between bariers is one
!      print *,'NSnew,P0,Plo=',NSnew,P0,Plo
                                !since there are no relevant variables.
      IF (NSnew.lt.1) GOTO 300

      XIND(1)=P0
      XIND(2)=Plo
                         !Lower bound Plo and upper bound P0 are close
                         !or all variables are close to be irrelevant,
                         !e.g. Nsnew=1.
      IF ((P0-Plo.LT.EPSS).OR.(P0.GT.CEPSS)) go to 300

! Now CEPSS>P0>EPSS and there are more then one relevant variables (NSnew>1)
! Those has indecies ind(I0), ind(I1).
! Hence we have nondegenerated case.
      
      SS0 = R (ind(I0) ,ind(I0))
      SQ0 = sqrt(SS0)  
!      print *,'P0-Plo,SS0,Sq0',P0-Plo,SS0,Sq0
      r1=indS(ind(I0))
      XMA = MIN((Hup (r1)-Cm (ind(I0)))/SQ0,xCutOff)
      XMI = MAX((Hlo (r1)-Cm (ind(I0)))/SQ0,-xCutOff) 

!If NSnew=2 then we ccan compute the probability exactly and recurions stops.
      IF (NSnew.eq.2) then
         XIND(2)=XIND(1)
         I1=2
         if (I0.eq.2) I1=1
         SS1 = R (ind(I1) ,ind(I1))
         SQ1=sqrt(SS1)
         rho=(R(ind(1),ind(2))/SS0)*(R(ind(1),ind(2))/SS1)
         IF (ABS(rho).gt.1.d0) go to 300
         rho=sqrt(rho)
         if (R(ind(1),ind(2)).lt.0.d0) rho=-rho
         r2=indS(ind(I1))
         XMA1 = MIN((Hup (r2)-Cm (ind(I1)))/SQ1,xCutOff)
         XMI1 = MAX((Hlo (r2)-Cm (ind(I1)))/SQ1,-xCutOff)
!         print *,XMA1,XMI1,XMA,XMI,rho
         XIND(1)=NORM2DPRB(XMI,XMA,XMI1,XMA1,rho)
         XIND(2)=XIND(1)
!         print *,XIND
         GOTO 300         
       END IF
        !If  NITL=0 which means computations without conditioning on X(ind(I0))
      IF(NITL.lt.1) go to 300

!We have NITL>0 nad at least 3 variables in the indicator, ie.
!we will condition on X(ind(I0)).
!First we check whether one can use XSPLIT variant of integration.

      if ((XMA.GE.xCutOff).AND.(XMI.LT.-XSPLT)) THEN  ! (.FALSE.).AND.
         XMA=XMI
         XMI=-xCutOff
         SGN=-1.d0
      elseif ((XMA.GT.XSPLT).AND.(XMI.LE.-xCutOff)) THEN
         XMI=XMA
         XMA=xCutOff
         SGN=-1.d0
      else
         SGN=1.d0
         XIND2=0.d0
      endif

         ! Must allocate several variables to recursively
         ! transfer them to rindnit: Rnew, SQnew, CMnew, indSnew
         ! The variable B is used in computations of conditional mean and cov.
         ! The size is NSnew-1 (the relevant variables minus X(ind(I0)).

         ALLOCATE(indSNEW(1:NSnew-1))
         ALLOCATE(RNEW(NSnew-1,NSnew-1))
         ALLOCATE(CMnew(1:NSnew-1))
         ALLOCATE(SQnew(1:NSnew-1))
         ALLOCATE(B(1:NSnew-1))
                               !This DO loop is divided in two parts in order 
                               !to only work on the upper triangular of R 
         DO row=1,I0-1
           r1=ind(row)
           Rnew(row,row:I0-1)=R(r1,ind(row:i0-1))
           Rnew(row,I0:NSnew-1)=R(r1,ind(i0+1:NSnew))
           B(row)=R(r1,ind(i0))
         enddo
         DO row=I0+1,NSnew
              r1=ind(row)
              Rnew(row-1,row-1:NSnew-1)=R(r1,ind(row:NSnew))
              B(row-1)=R(ind(i0),r1)
         enddo
         DO row=I0+1,NSnew
              ind(row-1)=ind(row)
         enddo

         CMnew=CM(ind(1:NSnew-1))
         SQnew=SQ(ind(1:NSnew-1))
         indSnew=indS(ind(1:NSnew-1))

                                                     !USE the  XSPLIT variant
         IF (SGN.LT.0.d0) XIND2=RINDNIT(Rnew,SQnew,CMnew,indSnew,NITL-1)

                          ! Perform conditioning on X(I0)
       NSnew=NSnew-1
       DO row = 1, NSnew                           
         Rnew(row,row:NSnew) = Rnew(row,row:NSnew) -
     &         B(row)*(B(row:NSnew)/SS0)
         SS = RNEW(row,row)
         IF (SS.GE.EPS2) then
               SQNEW (row) = SQRT (SS)
            ELSE
               SQNEW(row) =0.d0
         END IF
       ENDDO

                                !See if we can Narrow down the limits
      CALL C1C2(XMI,XMA,CmNew,B,SQNEW,SQ0,indSnew)       
      XIND(1) =(FI (XMA) - FI (XMI))
      XIND(2)=0.d0
      !print *,'rindnit fi',xind
      IF (XIND(1).LT.EPSS) GOTO 200                                               
      
      !print *,'rindnit gaussle2' 	
      N1=0                      !  computing nodes for num. integration. 
      CALL GAUSSLE2 (N1, H1, XX1, XMI, XMA,LE2Qnr) 
                                ! new conditional covariance

      XIND = 0.d0
      !print *,'rindnit for loop'                                                          
      DO   J = 1, N1
         IF (H1(J).GT.CFxCutOff) THEN       
                                !XIND1=0.d0
            CMnew=CM(ind(1:NSnew)) + XX1(J)*(B / SQ0)            
            XIND1=RINDNIT(Rnew,SQnew,CMnew,indSnew,NITL-1)               
            XIND = XIND+XIND1 * H1 (J)
         END IF                                      
      ENDDO
200   CONTINUE
      IF (SGN.lt.0.d0) call swapRe(XIND(1),XIND(2))
      XIND=XIND2+SGN*XIND
                           !fix up round off errors and make sure 0=<xind<=1
      if (XIND(1).GT.1.d0) THEN
         XIND(1)=1.D0
      elseif (XIND(1).LT.0.D0) THEN
         XIND(1)=0.d0
      endif
      if (XIND(2).GT.1.d0) THEN
         XIND=1.D0
      elseif (XIND(2).LT.0.D0) THEN
         XIND(2)=0.d0
      endif
300   CONTINUE
      if (allocated(INDSNEW)) DEALLOCATE(INDSNEW)
      if (allocated(RNEW))  DEALLOCATE(RNEW)
      if (allocated(CmNEW)) DEALLOCATE(CmNEW)
      if (allocated(SQNEW))  DEALLOCATE(SQNEW)
      if (allocated(B))  DEALLOCATE(B)
      if (allocated(ind)) deallocate(ind)
      !print *,'rindnit leaving end' 
      RETURN                                                             
      END FUNCTION  RINDNIT
                                          
                                !*****************************
      SUBROUTINE BARRIER(xc,indI,Blo,Bup)     
      USE GLOBALDATA, ONLY : Hup,Hlo,Mb,NI,xedni,Ntd
      IMPLICIT NONE
      INTEGER,          DIMENSION(:  ), INTENT(in) :: indI               
      DOUBLE PRECISION, DIMENSION(:,:), INTENT(in) :: Blo,Bup   
      DOUBLE PRECISION, DIMENSION(:  ), INTENT(in) :: xc         
      INTEGER                                      :: I, J, K, L           
!this procedure set Hlo,Hup according to Blo/Bup 
    
      DO J = 2, NI 
         DO I =indI (J - 1) + 1 , indI (J)  
            L=xedni(I)
            Hlo (L) = Blo (1, J - 1)  
            Hup (L) = Bup (1, J - 1)  
            DO K = 1, Mb-1                                                  
               Hlo(L) =Hlo(L)+Blo(K+1,J-1)*xc(K)
               Hup(L) =Hup(L)+Bup(K+1,J-1)*xc(K)
            ENDDO ! K
         ENDDO ! I
      ENDDO ! J
      !print * ,'barrier hup:',Hup(xedni(1:Ntd))

      !print * ,'barrier hlo:',Hlo(xedni(1:Ntd))
      
      RETURN                                                             
      END SUBROUTINE BARRIER                                            
  
! ************************************

      SUBROUTINE C1C2(C1, C2, Cm, B1, SQ, SQ0, ind)  
! The regression equation for the conditional distr. of Y given X=x
! is equal  to the conditional expectation of Y given X=x, i.e.,
! 
!       E(Y|X=x)=E(Y)+Cov(Y,X)/Var(X)[x-E(X)]
!
!  Let x1=(x-E(X))/SQ0 be zero mean, C1<x1<C2, B1(I)=COV(Y(I),X) and 
!  SQ0=sqrt(Var(X)). Then the process  Y(I) with mean Cm(I) can be written as 
!
!       y(I)=Cm(I)+x1*B1(I)/SQ0+Delta(I) for  I=1,...,N.
!
!  where SQ(I)=sqrt(Var(Y|X)) is the standard deviation of Delta(I). 
!
!  Since we are truncating all Gaussian  variables to                   
!  the interval [-C,C], then if for any I                               
!                                                                       
!  a) Cm(I)+x1*B1(I)/SQ0-C*SQ(I)>Hup(I)  or                             
!                                                                       
!  b) Cm(I)+x1*B1(I)/SQ0+C*SQ(I)<Hlo(I)  then                           
!                                                                       
!  (XIND|X1=x1) = 0 !!!!!!!!!                                               
!                                                                       
!  Consequently, for increasing the accuracy (by excluding possible 
!  discontinuouities) we shall exclude such values for which (XIND|X1=x1) = 0.
!  Hence we assume that if C1<x<C2 any of the previous conditions are 
!  satisfied
!                                                                       
!  OBSERVE!!, C1, C2 has to be set to (the normalized) lower and upper bounds 
!  of possible values for x1,respectively, i.e.,
!           C1=max((Hlo-E(X))/SQ0,-C), C2=min((Hup-E(X))/SQ0,C) 
!  before calling C1C2 subroutine.                         
!                                                
      USE GLOBALDATA, ONLY : Hup,Hlo,xCutOff,EPS2
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: Cm, B1, SQ 
      INTEGER,          DIMENSION(:), INTENT(in) :: ind 
      DOUBLE PRECISION,            INTENT(inout) :: C1,C2 
      DOUBLE PRECISION,               INTENT(in) :: SQ0   
     
      ! local variables                          
      
      DOUBLE PRECISION :: CC1, CC2,CSQ,HHup,HHlo,BdSQ0 
      INTEGER :: N,I            !,changedLimits=0

                                !ind contains indices to the varibles 
                                !location in  Hlo and Hup 
      IF(C1.GE.C2) THEN 
         C1=0.d0
         C2=0.d0
         RETURN                                                          
      ENDIF

      N=SIZE(ind)
      IF (N.LT.1) THEN                                                   
         RETURN    !Not able to change integration limits
      ENDIF

    
      DO I = N,1,-1            ! C=xCutOff
         CSQ = xCutOff*SQ(I)

         HHup = Hup (ind(I)) - Cm (I)                                       
         HHlo = Hlo (ind(I)) - Cm (I)
                                !  If ABS(B1(I)/SQ0) < EPS2 overflow may occur 
                                !  and hence if
                                !  1) Cm(I) is so large or small so we can 
                                !     surely assume that the probability
                                !     of staying between the barriers is 0, 
                                !     consequently C1=C2=0        
         BdSQ0=B1 (I)/SQ0
!         print *,'I,HHup,HHlo,Bdsq0',I,HHup,HHlo,BdsQ0
         IF (ABS (BdSQ0 ) .LT.EPS2 ) THEN
            if(SQ(I).lt.EPS2) CSQ= xCutOff*sqrt(EPS2)       
               IF (HHlo.GT.CSQ.OR.HHup.LT. - CSQ) THEN                    
               C1 = 0.d0                                                   
               C2 = 0.d0
!               changedLimits=1
               GOTO 112                                                    
            ENDIF 
         ELSE  !  In other cases this part follows 
               !  from the description of the problem.
            IF (BdSQ0.LT.0.d0) THEN 
                CC2 = (HHlo - CSQ) / BdSQ0   
               CC1 = (HHup + CSQ) / BdSQ0 
            ELSE ! BdSQ0>0
               CC1 = (HHlo - CSQ) / BdSQ0                     
               CC2 = (HHup + CSQ) / BdSQ0
            ENDIF
            IF (C1.LT.CC1) THEN
               C1 = CC1         !                  changedLimits=1
               IF (C2.GT.CC2) THEN 
                  C2 = CC2 
               END IF
               IF (C1.GE.C2) THEN                                                 
                  C2 = 0.d0 
                  C1 = 0.d0
                  GOTO 112 
               ENDIF
            ELSEIF (C2.GT.CC2) THEN 
               C2 = CC2      !                  changedLimits=1
               IF (C1.GE.C2) THEN                                                 
                  C2 = 0.d0 
                  C1 = 0.d0
                  GOTO 112 
               ENDIF
            END IF
         ENDIF
      END DO
                                                        
     
 112  continue 
      !IF (changedLimits) THEN
       ! PRINT *,'C1C2=',C1,C2
      !END IF
      RETURN                                                             
      END SUBROUTINE C1C2

! ************************************

      function MVNFUN(W,BIG,CmN,xd,xc,Pl1,Pu1) RESULT (XIND)
      USE GLOBALDATA, ONLY : Hlo,Hup,xCutOff,Nt,Nd,Nj,Ntd,SQ,
     &     NsXtmj, NsXdj,indXtd,index1,useC1C2,C1C2det,EPS2
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(:  ), INTENT(in)    :: W  
      DOUBLE PRECISION, DIMENSION(:  ), INTENT(inout)    :: CmN  ! conditional mean
      DOUBLE PRECISION, DIMENSION(:,:), INTENT(in) :: BIG ! conditional covariance 
      DOUBLE PRECISION, DIMENSION(:  ), INTENT(inout) :: xd  ! integr. variables
      DOUBLE PRECISION, DIMENSION(:  ), INTENT(in)    :: xc  ! conditional values
      DOUBLE PRECISION, INTENT(in) :: Pl1,Pu1 !FI(XMI),FI(XMA) of variable 1
      DOUBLE PRECISION                                :: XIND
!local variables
      DOUBLE PRECISION :: Pl,Pu
      DOUBLE PRECISION :: X,Y,XMI,XMA,SQ0
      INTEGER          :: Nst,NstN,NsdN,Nst0,Nsd,Nsd0,K
      INTEGER          :: Ndleft,Ndjleft,Ntmj

!MVNFUN Multivariate Normal integrand function
! where the integrand is transformed from an integral 
! having integration limits Hl0 and Hup to an
! integral having  constant integration limits i.e.
!   Hup                               1 
! int jacob(xd,xc)*f(xd,xt)dxt dxd = int F2(W) dW
!Hlo                                 0
!
! W   - new transformed integration variables, valid range 0..1
!       The vector must have the length Nst0+Ntd-Nsd0
! BIG - conditional sorted covariance matrix (IN)
! CmN - conditional mean NB! INOUT vector 
!       must be initialized with the conditional 
!       mean E(Xd,Xt|Xc) before calling mvnfun
! xd  - variables to the jacobian INOUT variable 
!        need no initialization
! xc  - conditional variables (IN) 
! Pl1 = FI(XMI) for the first integration variable (IN)
! Pu1 = FI(XMA) ------||-------------------------------
 
      Nst = NsXtmj(Ntd+1) ! index to last stoch variable of Xt before conditioning on X(Ntd)
      Ntmj=Nt-Nj
      Nsd0=NsXdj(1)
      if (Nt.gt.Nj) then
         Nst0=NsXtmj(Ntmj)
      else
         Nst0=0
      endif
      Pl=Pl1
      Pu=Pu1
            
      Y=Pu-Pl
      if (Nd+Nj.EQ.0) then
         SQ0=SQ(1,1)
         goto 200
      endif
      Ndjleft=Nd+Nj
      Nsd = NsXdj(Ndjleft+1)    ! index to last stoch variable of Xd and Nj of Xt  before conditioning on X(Ntd)
      Ndleft=Nd
      SQ0=SQ(Ntd,Ntd)
                                !print *,'rindsci,nst,nsd,nd,nj',nst,nsd,Nd,Nj   
      !print *,'mvn start K loop'
      DO K=Ntd-1,Nsd0,-1
         X=FIINV(Pl+W(Ntd-K)*(Pu-Pl))
         IF (index1(K+1).GT.Nt) THEN ! isXd
            xd (Ndleft) =  CmN(K+1)+X*SQ0
            Ndleft=Ndleft-1                  
         END IF
         Nst    = NsXtmj(K+1)   ! # stoch. var. of Xt before conditioning on X(K)
         CmN(1:Nst)=CmN(1:Nst)+X*BIG(1:Nst,K+1)/SQ0
         CmN(Nsd:K)=CmN(Nsd:K)+X*BIG(Nsd:K,K+1)/SQ0

         Ndjleft = Ndjleft-1
         Nsd      = NsXdj(Ndjleft+1)
         SQ0      = SQ(K,K)
                                 
         XMA = (Hup (K)-CmN(K))/SQ0
         XMI = (Hlo (K)-CmN(K))/SQ0
               
         if (useC1C2) then      ! see if we can narrow down sampling range 
                                !                  XMI=max(XMI,-xCutOff)
                                !                  XMA=min(XMA,xCutOff)
            if (C1C2det) then
               NsdN = NsXdj(Ndjleft) 
               NstN = NsXtmj(K)
               CALL C1C2(XMI,XMA,CmN(Nsd:NsdN-1),
     &              BIG(Nsd:NsdN-1,K),SQ(Nsd:NsdN-1,K),
     &              SQ0,indXtd(Nsd:NsdN-1))
               CALL C1C2(XMI,XMA,CmN(NstN+1:Nst),
     &              BIG(NstN+1:Nst,K),SQ(NstN+1:Nst,K),
     &              SQ0,indXtd(NstN+1:Nst))
            else
               CALL C1C2(XMI,XMA,CmN(Nsd:K-1),BIG(Nsd:K-1,K),
     &              SQ(Nsd:K-1,Ntmj+Ndjleft),SQ0,indXtd(Nsd:K-1))
               CALL C1C2(XMI,XMA,CmN(1:Nst),BIG(1:Nst,K)
     &              ,SQ(1:Nst,Ntmj+Ndjleft),SQ0,indXtd(1:Nst))
            endif     
            IF (XMA.LE.XMI) THEN
               Y=0.d0
               goto 250
            ENDIF
         endif
         Pl=FI(XMI)
         Pu=FI(XMA)
         Y=Y*(Pu-Pl) 
      ENDDO                     ! K LOOP
      X=FIINV(Pl+W(Ntd-Nsd0+1)*(Pu-Pl))
      Nst    = NsXtmj(Nsd0)     ! # stoch. var. of Xt after conditioning on X(Nsd0)
                                ! and before conditioning on X(1)
      CmN(1:Nst)=CmN(1:Nst)+X*(BIG(1:Nst,Nsd0)/SQ0)
      if (Nd.gt.0) then
         CmN(Nsd:Nsd0-1)=CmN(Nsd:Nsd0-1)+X*BIG(Nsd:Nsd0-1,Nsd0)/SQ0
         if (Ndleft.gt.0) then
            if (index1(Nsd0).GT.Nt) then     
               xd (Ndleft) =  CmN(Nsd0)+X*SQ0
               Ndleft=Ndleft-1
            endif
            K=Nsd0-1  
            do while (Ndleft.gt.0)
               if ((index1(K).GT.Nt)) THEN ! isXd
                  xd (Ndleft) =  CmN(K)
                  Ndleft=Ndleft-1                  
               END IF
               K=K-1
            ENDDO
         endif                  ! Ndleft
         Y = Y*jacob ( xd,xc)   ! jacobian of xd,xc
      endif                     ! Nd>0
      if (Nst0.gt.0) then
         SQ0=SQ(1,1)
         XMA = (Hup (1)-CmN(1))/SQ0
         XMI = (Hlo (1)-CmN(1))/SQ0
               
         if (useC1C2) then
                                !                  XMI=max(XMI,-xCutOff)
                                !                  XMA=min(XMA,xCutOff)
            if (C1C2det) then
               NstN = NsXtmj(1) ! # stoch. var. after conditioning 
               CALL C1C2(XMI,XMA,CmN(NstN+1:Nst),
     &              BIG(1,NstN+1:Nst),SQ(NstN+1:Nst,1),
     &              SQ0,indXtd(NstN+1:Nst))
            else
               CALL C1C2(XMI,XMA,CmN(2:Nst),BIG(1,2:Nst),
     &              SQ(2:Nst,1),SQ0,indXtd(2:Nst))
            endif
     
           
         endif
         IF (XMA.LE.XMI) THEN
            Y=0.d0
            goto 250
         ENDIF
         Pl=FI(XMI)
         Pu=FI(XMA)
         Y=Y*(Pu-Pl)
      endif
      !if (COVix.gt.2) then 
      !print *,' mvnfun start K2 loop'
      !endif
 200  do K=2,Nst0
         X=FIINV(Pl+W(Ntd-Nsd0+K)*(Pu-Pl)) 
         Nst = NsXtmj(K-1)      ! index to last stoch. var. before conditioning on X(K)
         CmN(K:Nst)=CmN(K:Nst)+X*BIG(K-1,K:Nst)/SQ0
         SQ0=SQ(K,K)
         if (SQ0.lt.EPS2) then
            print *,'rindscis index1',index1
            print *,'SQ'
            call echo(SQ(1:Ntd,1:Nt))
            print *,'error rindscis SQ0,',SQ0
         endif
         XMA = (Hup (K)-CmN(K))/SQ0
         XMI = (Hlo (K)-CmN(K))/SQ0
     
               
         if (useC1C2) then
                                !              XMI=max(XMI,-xCutOff)
                                !              XMA=min(XMA,xCutOff)
            if (C1C2det) then
               NstN=NsXtmj(K)   ! index to last stoch. var. after conditioning  X(K)
               CALL C1C2(XMI,XMA,CmN(NstN+1:Nst),
     &              BIG(K,NstN+1:Nst),SQ(NstN+1:Nst,K),
     &              SQ0,indXtd(NstN+1:Nst))
            else
               CALL C1C2(XMI,XMA,CmN(K+1:Nst),BIG(K,K+1:Nst),
     &              SQ(K+1:Nst,K),SQ0,indXtd(K+1:Nst))
            endif

         endif
         IF (XMA.LE.XMI) THEN
            Y=0.d0
            goto 250
         ENDIF
         Pl=FI(XMI)
         Pu=FI(XMA)               
         Y=Y*(Pu-Pl)
      enddo                     ! K loop
 250  XIND=Y
      !print *,' mvnfun leaving'
      return
      END FUNCTION MVNFUN

      function MVNFUN2(W,BIG,CmN,Pl1,Pu1) RESULT (XIND)
      USE GLOBALDATA, ONLY : Hlo,Hup,xCutOff,Njj,Nj,Ntscis,Ntd,SQ,
     &     NsXtmj, NsXdj,indXtd,index1,useC1C2,C1C2det,Nt,EPS2
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(:  ), INTENT(in)    :: W  
      DOUBLE PRECISION, DIMENSION(:  ), INTENT(inout) :: CmN  ! conditional mean
      DOUBLE PRECISION, DIMENSION(:,:), INTENT(inout) :: BIG ! conditional covariance 
      DOUBLE PRECISION, INTENT(in) :: Pl1,Pu1 !FI(XMI),FI(XMA) of variable 1
      DOUBLE PRECISION                                :: XIND
!local variables
      DOUBLE PRECISION :: Pl,Pu
      DOUBLE PRECISION :: X,Y,XMI,XMA,SQ0
      INTEGER          :: Nst,NstN,Nst0,K

!MVNFUN2 Multivariate Normal integrand function
! where the integrand is transformed from an integral 
! having integration limits Hl0 and Hup to an
! integral having  constant integration limits i.e.
!   Hup                  1 
! int f(xd,xt)dxt dxd = int F2(W) dW
!Hlo                    0
!
! W   - new transformed integration variables, valid range 0..1
!       The vector must have the size Nst0
! BIG - conditional sorted covariance matrix (IN)
! CmN - conditional mean NB! INOUT vector 
!       must be initialized with the conditional 
!       mean E(Xd,Xt|Xc) before calling mvnfun
! Pl1 = FI(XMI) for the first integration variable
! Pu1 = FI(XMA) ------||-------------------------
 
 
      Nst0 = NsXtmj(Njj+Ntscis)
      if (Njj.GT.0) then
         Nst  = NsXtmj(Njj)
      else
         Nst  = NsXtmj(Ntscis+1)
      endif
      !CmN(1:Nst)=Cm(1:Nst)
            
      Pl=Pl1
      Pu=Pu1
            
      Y=Pu-Pl
      SQ0=SQ(1,1)
      
      do K=2,Nst0
         X=FIINV(Pl+W(K-1)*(Pu-Pl))   
         Nst = NsXtmj(K-1)      ! index to last stoch. var. before conditioning on X(K)
         CmN(K:Nst)=CmN(K:Nst)+X*BIG(K-1,K:Nst)/SQ0
         SQ0=SQ(K,K)
         if (SQ0.lt.EPS2) then
            print *,'rindscis index1',index1
            print *,'SQ'
            call echo(SQ(1:Ntd,1:Nt))
            print *,'error rindscis SQ0,',SQ0
         endif
         XMA = (Hup (K)-CmN(K))/SQ0
         XMI = (Hlo (K)-CmN(K))/SQ0
       
               
         if (useC1C2) then
                                !              XMI=max(XMI,-xCutOff)
                                !              XMA=min(XMA,xCutOff)
            if (C1C2det) then
               NstN=NsXtmj(K)   ! index to last stoch. var. after conditioning on X(K)
               CALL C1C2(XMI,XMA,CmN(NstN+1:Nst),
     &              BIG(K,NstN+1:Nst),SQ(NstN+1:Nst,K),
     &              SQ0,indXtd(NstN+1:Nst))
            else
               CALL C1C2(XMI,XMA,CmN(K+1:Nst),BIG(K,K+1:Nst),
     &              SQ(K+1:Nst,K),SQ0,indXtd(K+1:Nst))
            endif
            
         endif
         IF (XMA.LE.XMI) THEN
            Y=0.d0
            goto 250
         ENDIF
         Pl=FI(XMI)
         Pu=FI(XMA)               
         Y=Y*(Pu-Pl)
      enddo                     ! K loop
 250  XIND=Y
      return
      END FUNCTION MVNFUN2



      function RINDSCIS(BIG,Cm,xd,xc) RESULT (XIND)
      USE GLOBALDATA, ONLY : Hlo,Hup,xCutOff,NUGGET,EPSS,EPS2,
     &     CFxCutOff, RelEps,NSIMmax,NSIMmin,Nt,Nd,Nj,Ntd,SQ,
     &     NsXtmj, NsXdj,indXtd,rateLHD,index1,
     &     useC1C2,useMIDP,C1C2det,MLHD,COV,COVix
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(:  ), INTENT(in)    :: Cm  ! conditional mean
      DOUBLE PRECISION, DIMENSION(:,:), INTENT(inout) :: BIG ! conditional covariance 
      DOUBLE PRECISION, DIMENSION(:  ), INTENT(inout) :: xd  ! integr. variables
      DOUBLE PRECISION, DIMENSION(:  ), INTENT(in)    :: xc  ! conditional values
      DOUBLE PRECISION                                :: XIND
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE  :: PP,V,C
      DOUBLE PRECISION, DIMENSION(:  ), ALLOCATABLE  :: CmN,D
      INTEGER,          DIMENSION(:,:), ALLOCATABLE  :: lhd

      DOUBLE PRECISION :: Pl1,Pu1,intval,intvalold,varsum,varsum2
      DOUBLE PRECISION :: ErrEst,ErrEst2, RelErrEst,RelErrEst2
      DOUBLE PRECISION :: Y,Y2,dL,dN,dNlhd,XMI,XMA,SQ0
      INTEGER          :: Nst,Nst0,Nsd,Nsd0,L,K,M,Nlhd
      INTEGER          :: Nrep,Ndim,Ndleft,Ntmj

   
    

!RINDSCIS Multivariate Normal integrals by SCIS or LHSCIS
!  SCIS   = Sequential conditioned importance sampling
!  LHSCIS = Latin Hypercube Sequential Conditioned Importance Sampling
!
! References 
! R. Ambartzumian, A. Der Kiureghian, V. Ohanian and H.
! Sukiasian (1998)
! Probabilistic Engineering Mechanics, Vol. 13, No 4. pp 299-308
!
! Alan Genz (1992)
! 'Numerical Computation of Multivariate Normal Probabilities'
! J. computational Graphical Statistics, Vol.1, pp 141--149


      !print *,'enter rindscis'
      Nst  = NsXtmj(Ntd+1)
      Ntmj=Nt-Nj
      if (Ntmj.GT.0) then
         Nst0 = NsXtmj(Ntmj)
      else
         Nst0 = 0
      endif
      Nsd  = NsXdj(Nd+Nj+1)
      Nsd0 = NsXdj(1)
      !if (COVix.Ge.33) then
         
      !   print *,'SQ'
      !   call echo(SQ(1:Ntd,1:Nt))
      !   print *,'rindscis index1',index1
      !print *,'Rindscis NsXtmj', NsXtmj
      !print *,'rindscis NsXdj', NsXdj
      !endif

                      !check deterministic means are within barrier
      if (any((Hup(Nst+1:Nsd-1).lt.Cm(Nst+1:Nsd-1)) .OR.
     &     (Hlo(Nst+1:Nsd-1).gt.Cm(Nst+1:Nsd-1)))) then
         XIND=0.d0
         if (allocated(COV)) then ! save the coefficient of variation in COV
            COV(COVix)=0.d0
         endif
         return
      endif
      
      if (Nd+Nj.gt.0) then
         IF ( BIG(Ntd,Ntd).LT.EPS2) THEN  !degenerate case                  
            IF (Nd.GT.0) THEN
               Ndleft=Nd;K=Ntd
               DO WHILE (Ndleft.GT.0)
                  IF (index1(K).GT.Nt) THEN ! isXd
                     xd (Ndleft) =  Cm (K)
                     Ndleft=Ndleft-1                  
                  END IF
                  K=K-1
               ENDDO
               XIND = jacob ( xd,xc) ! jacobian of xd,xc
            ELSE
               XIND = 1.d0      !     XIND = FxCutOff???
            END IF
            !print *,'jacob,xd',xind,xd            
            IF (Nst.LT.1) then
               if (allocated(COV)) then ! save the coefficient of variation in COV
                  COV(COVix)=0.d0
               endif
               RETURN  
            endif
            !print *,'rindsci '
            XIND=XIND*MNORMPRB(BIG(1:Nst,1:Nst),Cm(1:Nst))
            !print *,'leaving rindscis'
            RETURN
         ENDIF
      elseif (Nst.lt.1) then
         if (allocated(COV)) then ! save the coefficient of variation in COV
            COV(COVix)=0.d0
         endif
      
         XIND=1.d0
         return
      endif
      !print *,' rindscis start calculat'
      intval    = 0.d0
      intvalold = 0.d0
      Varsum    = 0.d0
      Varsum2   = 0.d0
      
     


      if (Nd+Nj.gt.0) then
         SQ0=SQ(Ntd,Ntd)
         XMA = (Hup (Ntd)-Cm(Ntd))/SQ0
         XMI = (Hlo (Ntd)-Cm(Ntd))/SQ0
     
         if (useC1C2) then ! see if we can narrow down sampling range 
!            XMI=max(XMI,-xCutOff)
!            XMA=min(XMA,xCutOff)
            CALL C1C2(XMI,XMA,Cm(1:Ntd-1),BIG(1:Ntd-1,Ntd),
     &           SQ(1:Ntd-1,Ntd),SQ0,indXtd(1:Ntd-1))
            
            IF (XMA.LE.XMI) THEN
                                ! XIND=Y=0 for all return 
               goto 300
            ENDIF
         endif
         Pl1=FI(XMI)
         Pu1=FI(XMA)
      else
         SQ0=SQ(1,1)
         XMA = (Hup (1)-Cm(1))/SQ0
         XMI = (Hlo (1)-Cm(1))/SQ0
     
         if (useC1C2) then ! see if we can narrow down sampling range 
!            XMI=max(XMI,-xCutOff)
!            XMA=min(XMA,xCutOff)
            CALL C1C2(XMI,XMA,Cm(2:Nst),BIG(1,2:Nst),
     &           SQ(2:Nst,1),SQ0,indXtd(2:Nst))
            
            IF (XMA.LE.XMI) THEN
                                ! PQ= Y=0 for all return 
               goto 300
            ENDIF
         endif
         Pl1=FI(XMI)
         Pu1=FI(XMA)
      endif
      Ndim=Nst0+Ntd-Nsd0+1      ! # dim. we treat stochastically 
      Nlhd=max(1,rateLHD*Ndim)             ! size of LHD
      if (rateLHD.lt.1) then ! make sure 
         usemidp=.FALSE.
         print * ,'rindscis: only able to use useMIDP if rateLHD>0'
      endif 
      ALLOCATE(CmN(1:Ntd))
      CmN(Nst+1:Nsd)= Cm(Nst+1:Nsd)
      allocate(pp(1:Nlhd,1:Ndim))
      if (nlhd.GT.1) then
         allocate(lhd(1:Nlhd,1:Ndim)) ! allocate LHD  
         allocate(C(1:Ndim,1:Ndim))    
         if (MLHD) then
            allocate(D(1:Ndim))
            allocate(V(1:Ndim,1:Ndim))
         endif
      end if
 
      dNlhd=dble(Nlhd)
      Nrep=ceiling(dble(NSIMmax)/dNlhd) ! # replications
      dL=0.d0
      dN=0.d0     
      !print *,' rindscis start L loop'
      DO L=1,Nrep
         CALL random_number(PP)
         if (Nlhd.gt.1) then
            LHD(:,1)=(/ (K,K=1,Nlhd)/)
            do M=1,3 ! do 3 attemps to construct a LHD with rank=Ndim
               do K=2,Ndim
                  CALL sortre(LHD(:,K),PP(:,K)) ! lhd = latin hypercube design
               enddo
               CALL spearcorr(C,lhd) ! find rankcorrelation between columns
               do K=1,Ndim-1 ! see if rank=Ndim
                  if (any(abs(C(K,K+1:Ndim)).GE.1.d0))  then
                     if (M.EQ.3) goto 30
                     CALL random_number(PP)
                     cycle
                  endif
               enddo
               goto 20
            enddo
 20         if (MLHD) then   !modify lhd by reducing correlation between columns
               DO K=1,Ndim
                  C(K,K)=C(K,K)+NUGGET ! add nugget effect to ensure that 
                                !inversion is not corrupted by round off errors
               enddo
               CALL svdcmp(C,D,V) ! C=U*D*V'=Q*Q'
               do K=1,Ndim
                 V(K,:)=C(:,K)*sqrt(1/D(K)) ! inverting Q=U*sqrt(D)
               enddo
              
               PP=MATMUL(dble(lhd),V)  ! LHD*inv(Q)
               do K=1,Ndim
                  CALL sortre(LHD(:,K),PP(:,K)) ! lhd = latin hypercube design
               enddo
            endif 
 30         if (USEMIDP) then   ! use the center of the cell
               PP=(dble(lhd)-0.5d0)/dNlhd   
            else                ! distribute uniformly within the cell
               CALL random_number(PP)
               PP=(dble(lhd)-PP)/dNlhd 
            endif
         endif

           
         !print *,' rindscis start M loop'
         DO M=1,Nlhd
            Nst = NsXtmj(Ntd+1) 
            CmN(1:Nst)=Cm(1:Nst) ! initialize conditional mean  
            CmN(Nsd:Ntd)=Cm(Nsd:Ntd)                          
            Y=MVNFUN(PP(M,:),BIG,CmN,xd,xc,Pl1,Pu1) ! evaluate the integrand
            CmN(1:Nst)=Cm(1:Nst) ! initialize conditional mean
            CmN(Nsd:Ntd)=Cm(Nsd:Ntd)
            PP(M,:)=1.d0-PP(M,:) ! using antithetic variables to reduce variance
            Y= (Y+MVNFUN(PP(M,:),BIG,CmN,xd,xc,Pl1,Pu1))/2.d0   
      
            dN=dN+1.d0
            varsum= varsum+(dN-1.d0)*(Y-intval)*(Y-intval)/dN 
            intval=intval+(Y-intval)/dN ! integral value
            ErrEst=2.5d0*SQRT(varsum/(dN*dN)) ! error estimate
            if (abs(intval).gt.0.d0) then
               RelErrEst = ErrEst/abs(intval)
            else
               RelErrEst=0.d0
            end if
                       
            IF (((RelErrEst.LT.RelEps).or.(ErrEst.lt.EPSS))
     &              .AND.(dN.GT.NSIMmin))  GOTO 300
         ENDDO                  ! M loop
         if (Nlhd.gt.1) then
            dL=dble(L)
            Y2=dL*(intval-intvalold)+intvalold
            varsum2= varsum2+(dL-1.d0)*(Y2-intval)*(Y2-intval)/dL ! Better estimate of variance
            ErrEst2=2.5d0*SQRT(varsum2/(dL*dL)) ! error estimate between replicates
            if (abs(intval).gt.0.d0) then
               RelErrEst2 = ErrEst/abs(intval)
            else
               RelErrEst2=0.d0
            end if
                       
            IF (((RelErrEst2.LT.RelEps).or.(ErrEst2.lt.EPSS))
     &           .AND.(L.GT.5))  GOTO 300
            intvalold=intval
         end if
       enddo                     ! L loop
 300  XIND=intval
      !print *,'rindscis L,N,PQ,CV,CV2',L,dN,PQ,CV,CV2
      !print *,'rindscis L,N,PQ',L,dN,PQ
      !print *,'leaving rindscis',PQ
      if (allocated(COV)) then ! save the coefficient of variation in COV
         if ((dL.gt.1).and.(Nlhd.gt.1)) then
            COV(COVix)=min(RelErrEst,RelErrEst2)/2.5d0
         else
            COV(COVix)=RelErrEst/2.5d0
         endif
      endif
      if (allocated(PP)) then
         DEALLOCATE(CmN)
         DEallocate(pp)
      endif
      if (allocated(lhd)) then
         DEallocate(lhd)
         DEallocate(C)
         if (allocated(D)) then
            DEallocate(D)
            DEallocate(V)
         endif
      endif
      !print *,'leaving rindscis2'
      return
      END FUNCTION RINDSCIS

!********************************************************************

      SUBROUTINE CONDSORT (R,CSTD,index1,xedni,NsXtmj,NsXdj)
      USE GLOBALDATA, ONLY : Nt,Nj,Njj,Nd,Nc,Ntdc,Ntd,EPS2,Nugget,
     &    XCEPS2,SCIS,Ntscis
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(:,:), INTENT(inout) :: R
      DOUBLE PRECISION, DIMENSION(:,:), INTENT(out)   :: CSTD 
      INTEGER,          DIMENSION(:  ), INTENT(out)   :: index1 
      INTEGER,          DIMENSION(:  ), INTENT(out)   :: xedni 
      INTEGER,          DIMENSION(:  ), INTENT(out)   :: NsXtmj
      INTEGER,          DIMENSION(:  ), INTENT(out)   :: NsXdj
! local variables
      DOUBLE PRECISION, DIMENSION(:  ), allocatable   :: SQ
      DOUBLE PRECISION, DIMENSION(:,:), allocatable   :: CSTD2
      INTEGER,          DIMENSION(1  )                :: m
      INTEGER       :: Nstoc,Ntmp,NstoXd   !,degenerate
      INTEGER       :: changed,m1,r1,c1,row,col,ix,iy,Njleft,Ntmj
     
! R         = Input: Cov(X) where X=[Xt Xd Xc] is stochastic vector
!            Output: sorted Conditional Covar. matrix   Shape N X N  (N=Nt+Nd+Nc) 
! CSTD      = SQRT(Var(X(1:I-1)|X(I:N)))
!            conditional standard deviation.             Shape Ntd X max(Nd+Nj,1)
! index1    = indices to the variables original place.   Size  Ntdc
! xedni     = indices to the variables new place.        Size  Ntdc
! NsXtmj(I) = indices to the last stochastic variable
!            among Nt-Nj first of Xt after conditioning on 
!            X(Nt-Nj+I).                                 Size  Nd+Nj+Njj+Ntscis+1
! NsXdj(I)  = indices to the first stochastic variable
!            among Xd+Nj of Xt after conditioning on 
!            X(Nt-Nj+I).                                 Size  Nd+Nj+1
! 
! R=Cov([Xt,Xd,Xc]) is a covariance matrix of the stochastic vector X=[Xt Xd Xc]
! where the variables Xt, Xd and Xc have the size Nt, Nd and Nc, respectively. 
! Xc is (are) the conditional variable(s).
! Xd and Xt are the variables to integrate.
! Xd + Nj variables of Xt are integrated directly by the RindDXX 
! subroutines in the order of decreasing conditional variance.
! The remaining Nt-Nj variables of Xt are integrated in random
! order by the RindXX subroutines. 
! CONDSORT prepare and rearrange the covariance matrix 
! by decreasing order of conditional variances in a special way
! to accomodate this strategy:
! 
! After conditioning and sorting, the first Nt-Nj x Nt-Nj block of R 
! will contain the conditional covariance matrix
! of Xt(1:Nt-Nj) given Xt(Nt-Nj+1:Nt)  Xd and Xc, i.e., 
! Cov(Xt(1:Nt-Nj),Xt(1:Nt-Nj)|Xt(Nt-Nj+1:Nt), Xd,Xc)     
! NB! for Nj>0 the order of Xd and Xt(Nt-Nj+1:Nt) may be mixed.
! The covariances, Cov(X(1:I-1),X(I)|X(I+1:N)), needed for computation of the 
! conditional expectation, E(X(1:I-1)|X(I:N), are saved in column I of R 
! for I=Nt-Nj+1:Ntdc.
! 
! IF any of the variables have variance less than EPS2. They will be 
! be treated as deterministic and not stochastic variables by the 
! RindXXX subroutines. The deterministic variables are moved to 
! middle in the order they became deterministic in order to 
! keep track of them. Their variance and covariance with
! the remaining stochastic variables are set to zero in 
! order to avoid numerical difficulties.
! 
! NsXtmj(I) is the number of variables  among the Nt-Nj 
! first we treat stochastically after conditioning on X(Nt-Nj+I).  
! The covariance matrix is sorted so that all variables with indices
! from 1 to NsXtmj(I) are stochastic after conditioning
! on X(Nt-Nj+I).  Thus NsXtmj(I) may also be considered 
! as the index to the last stochastic variable after conditioning
! on X(Nt-Nj+I). In other words NsXtmj keeps track of the deterministic
! and stochastic variables among the Nt-Nj first variables in each 
! conditioning step.
!
! Similarly  NsXdj(I)  keeps track of the deterministic and stochastic
! variables among the Nd+Nj following variables in each conditioning step.
! NsXdj(I) is the index to the first stochastic variable
! among the Nd+Nj following variables after conditioning on X(Nt-Nj+I).   
! The covariance matrix is sorted so that all variables with indices
! from NsXdj(I+1) to NsXdj(I)-1 are  deterministic conditioned on
! X(Nt-Nj+I). 
!

! Var(Xc(1))>Var(Xc(2)|Xc(1))>...>Var(Xc(Nc)|Xc(1),Xc(2),...,Xc(Nc)).
! If Nj=0 then
! Var(Xd(1)|Xc)>Var(Xd(2)|Xd(1),Xc)>...>Var(Xd(Nd)|Xd(1),Xd(2),...,Xd(Nd),Xc). 
!        
! NB!! Since R is symmetric, only the upper triangular contains the
! sorted conditional covariance. The whole matrix
! is easily obtained by copying elements of the upper triangle to 
! the lower or by uncommenting some lines in the end of this subroutine


! Using SQ to temporarily store the diagonal of R
! Adding a nugget effect to ensure the the inversion is 
! not corrupted by round off errors 
! good choice for nugget might be 1e-8         
                             !call getdiag(SQ,R)
      ALLOCATE(SQ(1:Ntdc))

      IF (Nd+Nj+Njj+Ntscis.GT.0) THEN
         ALLOCATE(CSTD2(1:Ntd,1:Nd+Nj+Njj+Ntscis))
         CSTD2=0.d0             ! initialize CSTD 
      ENDIF
      !CALL ECHO(R,Ntdc)
      DO ix = 1, Ntdc 
         R(ix,ix)=R(ix,ix)+Nugget   
         SQ(ix)=R(ix,ix)
         index1 (ix) = ix       ! initialize index1 
      ENDDO
     
      Ntmj=Nt-Nj
      !NsXtmj(Njj+Nd+Nj+1)=Ntmj      ! index to last stochastic variable of Nt-Nj of Xt
      !NsXdj(Nd+Nj+1)=Ntmj+1     ! index to first stochastic variable of Xd and Nj of Xt
      !degenerate=0
      Njleft=Nj
      NstoXd=Ntmj+1;Nstoc=Ntmj
      
     
      DO ix = 1, Nc             ! Condsort Xc
         m=Ntdc-ix+2-MAXLOC(SQ(Ntdc-ix+1:Ntd+1:-1)) 
         IF (SQ(m(1)).LT.XCEPS2) THEN
            PRINT *,'Condsort, degenerate Xc'
                                !degenerate=1
            GOTO 200            ! RETURN    !degenerate case
         ENDIF
         m1=index1(m(1));
         CALL swapint(index1(m(1)),index1(Ntdc-ix+1))
         SQ(Ntdc-ix+1)=SQ(m(1))
                                ! sort and calculate conditional covariances
         CALL CONDSORT2(R,SQ,index1,Nstoc,NstoXd,Njleft,m1,Ntdc-ix)
      ENDDO                     ! ix 
        
      NsXdj(Nd+Nj+1)  = NstoXd  ! index to first stochastic variable of Xd and Nj of Xt
      NsXtmj(Nd+Nj+Njj+Ntscis+1) = Nstoc ! index to last stochastic variable of Nt-Nj of Xt
      !print *, 'condsort index1', index1
      !print *, 'condsort Xd' 
      !call echo(R,Ntdc)
     
      DO ix = 1, Nd+Nj          ! Condsort Xd +  Nj of Xt
         IF (Njleft.GT.0) THEN
            m=Ntd-ix+2-MAXLOC(SQ(Ntd-ix+1:1:-1))
            Ntmp=NstoXd+Njleft-1
            IF (((NstoXd.LE.m(1)).AND.(m(1).LE.Ntmp))
     &           .OR.(m(1).LE.Nstoc)) THEN
               CALL swapint(index1(m(1)),index1(Ntmp))
               CALL swapRe(SQ(m(1)),SQ(Ntmp))
               m(1)=Ntmp
               Njleft=Njleft-1
            END IF
         ELSE
            m=Ntd-ix+2-MAXLOC(SQ(Ntd-ix+1:Ntmj+1:-1))
         END IF
         IF (SQ(m(1)).LT.EPS2) THEN          
                                !PRINT *,'Condsort, degenerate Xd'
                                !degenerate=1
            Ntmp=Nd+Nj+1-ix
            NsXtmj(Ntscis+Njj+1:Ntmp+Ntscis+Njj+1)=Nstoc
            NsXdj(1:Ntmp+1)=NstoXd
            IF (ix.EQ.1) THEN
               DO iy=1,Ntd      !sqrt(VAR(X(I)|X(Ntd-ix+1:Ntdc))
                  r1=index1(iy)
                  CSTD2(r1,Ntscis+Njj+1:Ntmp+Ntscis+Njj)=SQRT(SQ(iy)) 
               ENDDO
            ELSE
               DO iy=ix,Nd+Nj
                  CSTD2(:,Nd+Nj+Ntscis+Njj+1-iy)=
     &                 CSTD2(:,Ntmp+Ntscis+Njj+1)
               ENDDO
            ENDIF
            GOTO 200            ! degenerate case
         END IF
         m1=index1(m(1));
         CALL swapint(index1(m(1)),index1(Ntd-ix+1))
         CSTD2(m1,Nd+Nj+Ntscis+Njj+1-ix)=SQRT(SQ(m(1)))
         SQ(Ntd-ix+1)=SQ(m(1))
         
                                ! Calculating conditional variances 
         CALL CONDSORT2(R,SQ,index1,Nstoc,NstoXd,Njleft,m1,Ntd-ix)
                                ! saving indices
         NsXtmj(Nd+Nj+Njj+Ntscis+1-ix)=Nstoc
         NsXdj(Nd+Nj+1-ix)=NstoXd
         
                                ! Calculating standard deviations non-deterministic variables
         DO row=1,Nstoc
            r1=index1(row)
            CSTD2(r1,Nd+Nj+Njj+Ntscis+1-ix)=SQRT(SQ(row)) !sqrt(VAR(X(I)|X(Ntd-ix+1:Ntdc))  
         ENDDO            
         DO row=NstoXd,Ntd-ix
            r1=index1(row)
            CSTD2(r1,Nd+Nj+Ntscis+Njj+1-ix)=SQRT(SQ(row)) !sqrt(VAR(X(I)|X(Ntd-ix+1:Ntdc))  
         ENDDO
      ENDDO                     ! ix 
      
      
 200  IF ((SCIS.GT.0).OR. (Njj.gt.0)) THEN ! check on Njj instead
          ! Calculating conditional variances and sort for Nstoc of Xt
         CALL CONDSORT3(R,CSTD2,SQ,index1,NsXtmj,Nstoc)
         !Nst0=Nstoc
      ENDIF
      IF (Nd+Nj.EQ.0) THEN 
         IF (Nc.EQ.0) THEN 
            ix=1; Nstoc=Ntmj   
            DO WHILE (ix.LE.Nstoc)
               IF (SQ(ix).LT.EPS2) THEN
                  DO WHILE ((SQ(Nstoc).LT.EPS2).AND.(ix.LT.Nstoc))
                     SQ(Nstoc)=0.d0
                     Nstoc=Nstoc-1
                  END DO 
                  CALL swapint(index1(ix),index1(Nstoc)) ! swap indices
                  SQ(ix)=SQ(Nstoc); SQ(Nstoc)=0.d0
                  Nstoc=Nstoc-1
               ENDIF
               ix=ix+1
            END DO         
         ENDIF
         CSTD(1:Nt,1)=SQRT(SQ(1:Nt))
         NsXtmj(1)=Nstoc
      ELSE
         DO row=1,Ntd           ! sorting CSTD according to index1
            r1=index1(row)         
            CSTD(row,:)= CSTD2(r1,:)
         END DO
         DEALLOCATE(CSTD2)
      ENDIF                      

      changed=0                         
      DO row=Ntdc,1,-1          ! sorting the upper triangular of the 
         r1=index1(row)         ! covariance matrix according to index1
         xedni(r1)=row
         !PRINT *,'condsort,xedni',xedni
         !PRINT *,'condsort,r1,row',r1,row
         IF ((r1.NE.row).OR.(changed.EQ.1)) THEN
            changed=1
            R(row,row)=SQ(row)
            DO col=row+1,Ntdc
               c1=index1(col)
               IF (c1.GT.r1) THEN
                  R(row,col)=R(c1,r1)
               ELSE
                  R(row,col)=R(r1,c1)
               END IF
            END DO
         END IF              
      END DO
                                ! you may sort the lower triangular according  
                                ! to index1 also, but it is not needed
                                ! since R is symmetric.  Uncomment the
                                ! following if the whole matrix is needed
!      DO col=1,Ntdc
!         DO row=col+1,Ntdc
!            R(row,col)=R(col,row) ! R symmetric
!         END DO
!      END DO
!      IF (degenerate.EQ.1) THEN
!         PRINT *,'condsort,R='
!         call echo(R,Ntdc)
!         PRINT *,'condsort,SQ='
!         call echo(CSTD,Ntd)
!         PRINT *,'index=',index1
!         PRINT *,'xedni=',xedni
!      ENDIF
!      PRINT * , 'big'
!600   FORMAT(4F8.4)
!      PRINT 600, R
!      PRINT 600, SQ
      DEALLOCATE(SQ)
     
      RETURN                                                             
      END SUBROUTINE CONDSORT


      SUBROUTINE CONDSORT2(R,SQ,index1,Nstoc,NstoXd,Njleft,m1,N)
      USE GLOBALDATA, ONLY : Ntd,EPS2,XCEPS2
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(:,:), INTENT(inout) :: R
      DOUBLE PRECISION, DIMENSION(:),   INTENT(inout) :: SQ 
      INTEGER,          DIMENSION(:  ), INTENT(inout)   :: index1 
      INTEGER, INTENT(inout)   :: Nstoc,NstoXd,Njleft
      INTEGER, INTENT(in)   :: m1,N
! local variables
      INTEGER       :: Nsold,Ndold, Ntmp
      INTEGER       :: r1,c1,row,col,iy

! save their old values
      Nsold=Nstoc;Ndold=NstoXd

                                ! Calculating conditional variances for the 
                                ! Xc variables. 
      DO row=Ntd+1,N
         r1=index1(row)
         SQ(row)=R(r1,r1)-R(r1,m1)*R(m1,r1)/R(m1,m1)
         IF (SQ(row).LT.XCEPS2) THEN 
            R(r1,r1)=0.d0
            SQ(row)=0.d0
            RETURN              ! degenerate case XIND should return NaN
         ELSE
            R(r1,r1)=SQ(row)
            DO col=row+1,N
               c1=index1(col)
               R(c1,r1)=R(r1,c1)-R(r1,m1)*R(m1,c1)/R(m1,m1)
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
         SQ(iy)=R(r1,r1)-R(r1,m1)*R(m1,r1)/R(m1,m1)
         IF (SQ(iy).LT.EPS2) THEN
            r1=index1(Nstoc)
            SQ(Nstoc)=R(r1,r1)-R(r1,m1)*R(m1,r1)/R(m1,m1)
            DO WHILE ((SQ(Nstoc).LT.EPS2).AND.(iy.LT.Nstoc))
               SQ(Nstoc)=0.d0
               Nstoc=Nstoc-1
               r1=index1(Nstoc)
               SQ(Nstoc)=R(r1,r1)-R(r1,m1)*R(m1,r1)/R(m1,m1)
            END DO 
            CALL swapint(index1(iy),index1(Nstoc)) ! swap indices
            SQ(iy)=SQ(Nstoc); SQ(Nstoc)=0.d0 ! swap values
            Nstoc=Nstoc-1
         ENDIF
         iy=iy+1
      END DO  
      
                                ! Calculating conditional variances for the 
                                ! stochastic variables Xd and Njleft of Xt. 
                                ! Variables with conditional variance less than
                                ! EPS2 are moved to the beginning among these
                                ! with only One exception: if it is one of the
                                ! Xt variables and Nstoc>0 then it switch place 
                                ! with Xt(Nstoc)

      DO iy=Ndold,MIN(Ntd,N)
         r1=index1(iy)
         SQ(iy)=R(r1,r1)-R(r1,m1)*R(m1,r1)/R(m1,m1)
         IF (SQ(iy).LT.EPS2) THEN
            IF (Njleft.GT.0) THEN
               Ntmp=NstoXd+Njleft 
               IF (iy.LT.Ntmp) THEN
                  IF (Nstoc.GT.0) THEN !switch place with Xt(Nstoc)
                     CALL swapint(index1(iy),index1(Nstoc))
                     SQ(iy)=SQ(Nstoc);SQ(Nstoc)=0.d0
                     Nstoc=Nstoc-1
                  ELSE
                     CALL swapint(index1(iy),index1(NstoXd))
                     SQ(iy)=SQ(NstoXd);SQ(NstoXd)=0.d0
                     Njleft=Njleft-1
                     NstoXd=NstoXd+1
                  ENDIF
               ELSE
                  CALL swapint(index1(iy),index1(Ntmp))
                  CALL swapint(index1(Ntmp),index1(NstoXd))
                  SQ(iy)=SQ(Ntmp) 
                  SQ(Ntmp)=SQ(NstoXd);SQ(NstoXd)=0.d0
                  NstoXd=NstoXd+1
               ENDIF
            ELSE
               CALL swapint(index1(iy),index1(NstoXd))
               SQ(iy)=SQ(NstoXd);SQ(NstoXd)=0.d0
               NstoXd=NstoXd+1
            ENDIF
         ENDIF                  ! SQ < EPS2
      ENDDO
      
            
            ! Calculating Covariances for non-deterministic variables
      DO row=1,Nstoc
         r1=index1(row)
         R(r1,r1)=SQ(row)
         DO col=row+1,Nstoc
            c1=index1(col)
            R(c1,r1)=R(r1,c1)-R(r1,m1)*R(m1,c1)/R(m1,m1)
            R(r1,c1)=R(c1,r1)
         ENDDO
         DO col=NstoXd,N
            c1=index1(col)
            R(c1,r1)=R(r1,c1)-R(r1,m1)*R(m1,c1)/R(m1,m1)
            R(r1,c1)=R(c1,r1)
         ENDDO
      ENDDO            
      DO row=NstoXd,MIN(Ntd,N)
         r1=index1(row)
         R(r1,r1)=SQ(row)
         
         DO col=row+1,N
            c1=index1(col)
            R(c1,r1)=R(r1,c1)-R(r1,m1)*R(m1,c1)/R(m1,m1)
            R(r1,c1)=R(c1,r1)
         ENDDO
      ENDDO
      
                                ! Set covariances for Deterministic variables to zero
                                ! in order to avoid numerical problems
      
      DO row=Ndold,NStoXd-1
         r1=index1(row)
         R(r1,r1)=0.d0
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
         R(r1,r1)=0.d0
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

      SUBROUTINE CONDSORT3(R,CSTD2,SQ,index1,NsXtmj,Nstoc)
      USE GLOBALDATA, ONLY : EPS2,Njj,Ntscis
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(:,:), INTENT(inout) :: R,CSTD2
      DOUBLE PRECISION, DIMENSION(:),   INTENT(inout) :: SQ ! diag. of R
      INTEGER,          DIMENSION(:  ), INTENT(inout) :: index1,NsXtmj 
      INTEGER,          DIMENSION(1) :: m
      INTEGER, INTENT(inout)   :: Nstoc
! local variables
      INTEGER  :: m1
      INTEGER       :: Nsold
      INTEGER       :: r1,c1,row,col,iy,ix
! This function condsort all the xt variables for use with RINDSCIS and 
! MNORMPRB

      !Nsoold=Nstoc
      ix=1
     
      DO WHILE ((ix.LE.Nstoc).and.(ix.LE.(Ntscis+Njj)))
         m=ix-1+MAXLOC(SQ(ix:Nstoc)) 
         IF (SQ(m(1)).LT.EPS2) THEN
            PRINT *,'Condsort3, error degenerate X'
            NsXtmj(1:Njj+Ntscis)=0
            Nstoc=0       !degenerate=1
            RETURN    !degenerate case
         ENDIF
         m1=index1(m(1));
         CALL swapint(index1(m(1)),index1(ix))
         SQ(ix)=SQ(m(1))
         CSTD2(m1,ix)=SQRT(SQ(ix))
                                ! Calculating conditional variances for the 
                                ! first Nstoc variables.
                                ! variables with variance less than EPS2 
                                ! will be treated as deterministic and not 
                                ! stochastic variables and are therefore moved
                                ! to the end among these variables.
                                ! Nstoc is the # of variables we treat 
                                ! stochastically 
         iy=ix+1;Nsold=Nstoc
         DO WHILE (iy.LE.Nstoc)
            r1=index1(iy)
            SQ(iy)=R(r1,r1)-R(r1,m1)*R(m1,r1)/R(m1,m1)
            IF (SQ(iy).LT.EPS2) THEN
               r1=index1(Nstoc)
               SQ(Nstoc)=R(r1,r1)-R(r1,m1)*R(m1,r1)/R(m1,m1)
               DO WHILE ((SQ(Nstoc).LT.EPS2).AND.(iy.LT.Nstoc))
                  SQ(Nstoc)=0.d0
                  Nstoc=Nstoc-1
                  r1=index1(Nstoc)
                  SQ(Nstoc)=R(r1,r1)-R(r1,m1)*R(m1,r1)/R(m1,m1)
               END DO 
               CALL swapint(index1(iy),index1(Nstoc)) ! swap indices
               SQ(iy)=SQ(Nstoc); SQ(Nstoc)=0.d0 ! swap values
               Nstoc=Nstoc-1
            ENDIF
            iy=iy+1
         END DO 
         NsXtmj(ix)=Nstoc ! saving index to last stoch. var. after conditioning
             ! Calculating Covariances for non-deterministic variables
         DO row=ix+1,Nstoc
            r1=index1(row)
            R(r1,r1)=SQ(row)
            CSTD2(r1,ix)=SQRT(SQ(row)) ! saving stdev after conditioning on ix
            DO col=row+1,Nstoc
               c1=index1(col)
               R(c1,r1)=R(r1,c1)-R(r1,m1)*R(m1,c1)/R(m1,m1)
               R(r1,c1)=R(c1,r1)
            ENDDO
         ENDDO        
            ! similarly for deterministic values  
         DO row=Nstoc+1,Nsold
            r1=index1(row)
            R(r1,r1)=0.d0
            DO col=ix+1,Nsold   !row-1
               c1=index1(col)
               R(c1,r1)=0.d0
               R(r1,c1)=0.d0
            ENDDO
         ENDDO
         ix=ix+1
      ENDDO
      NsXtmj(Nstoc+1:Njj+Ntscis)=Nstoc
      RETURN   
      END SUBROUTINE CONDSORT3

      SUBROUTINE SORTRE(indices,rarray)
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(:), INTENT(inout) :: rarray
      INTEGER,          DIMENSION(:), INTENT(out)   :: indices
! local variables
       INTEGER  :: i,im,j,k,m,n
   
! diminishing increment sort as described by
! Donald E. Knuth (1973) "The art of computer programming,",
! Vol. 3, pp 84-  (sorting and searching)
      n=size(indices)
      indices=(/(i,i=1,n)/)
100   continue
      if (n.le.1) goto 800
      m=1
200   continue
      m=m+m
      if (m.lt.n) goto 200
      m=m-1
300   continue
      m=m/2
      if (m.eq.0) goto 800
      k=n-m
      j=1
400   continue
      i=j
500   continue
      im=i+m
      if (rarray(i).gt.rarray(im)) goto 700          
600   continue
      j=j+1
      if (j.gt.k) goto 300
      goto 400
700   continue
      CALL swapre(rarray(i),rarray(im))
      CALL swapint(indices(i),indices(im))
      i=i-m
      if (i.lt.1) goto 600
      goto 500
800   continue
      RETURN   
      END SUBROUTINE SORTRE

!______________________________________________________

      subroutine sobseq(n,x)
      double precision,dimension(:) :: x
      integer, intent(in) ::n
      integer,parameter ::MAXBIT=30,MAXDIM=6
      integer :: i,im, in,ipp,j,k,l
      integer, dimension(MAXDIM) :: ip,mdeg,ix
      integer, dimension(MAXDIM,MAXBIT) ::iu
      integer, dimension(MAXDIM*MAXBIT) ::iv
      double precision :: fac
      save ip,mdeg,ix,iv,in,fac
      equivalence (iv,iu)       ! to allow both 1D and 2D addressing
! returns sobols sequence of quasi-random numbers between 0 1
! When n is negative, internally initializes a set of MAXBIT
! direction numbers for each of MAXDIM different sobol
! sequences. When n is positive (but < MAXDIM)
! returns as the vector x(1:n) the next values from n of these sequences 
! (n must not be changed between initializations)
!
! This routine is initialised for maximum of n=6 dimensions
! and a word length of 30 bits. These parameter may be increased by 
!changing MAXBIT and MAXDIM and add more initializing data to 
! ip (primitive polynomials), mdeg (their degrees) and iv 
! (the starting value for the recurrence relation)
 
!reference
! William H. Press, Saul Teukolsky, William T. Wetterling and Brian P. Flannery (1997)
! "Numerical recipes in Fortran 77", Vol. 1, pp 299--305
      if (n.lt.0) then          ! initialize, do not return vector

         ix=0
         in=0  ! random startpoint: CALL random_number(P); in=P*2^MAXBIT
               ! and remove warning message below
         if (iv(1).NE.1) return

         ip=(/0,1,1,2,1,4 /)
         mdeg=(/1,2,3,3,4,4 /)
         ix=(/0,0,0,0,0,0 /)
         iv(1:24)=(/1,1,1,1,1,1,3,1,3,3,1,1,5,
     &        7,7,3,3,5,15,11,5,15,13,9/)
         iv(25:MAXDIM*MAXBIT)=0

         
         fac=1.d0/2.d0**MAXBIT
         do k=1,MAXDIM
            do j=1,mdeg(k)      ! stored values need normalization
               iu(k,j)=iu(k,j)*2**(MAXBIT-j)
            enddo
            do j=1,mdeg(k)+1,MAXBIT ! use reccurence to get other values
               ipp=ip(k)
               i=iu(k,j-mdeg(k))
               i=ieor(i,i/2**mdeg(k))
               do l=mdeg(k)-1,1,-1
                  if (iand(ipp,1).ne.0) i=ieor(i,iu(k,j-l))
                  ipp=floor(dble(ipp/2))
               enddo
               iu(k,j)=i
            enddo      
         enddo
      else                      ! calculate the next vector in the sequence
         im=in
         do j=1,MAXBIT          ! find the rightmost zero bit
            if (iand(im,1).eq.0) goto 1
            im=floor(dble(im/2))
         enddo
         print *,'MAXBIT too small in sobseq'
 1       im=(j-1)*MAXDIM
         do k=1,min(n,MAXDIM)   !XOR the 
            ix(k)=ieor(ix(k),iv(im+k))
            x(k)=ix(k)*fac
         enddo
         in=in+1                ! increment counter

      endif
      return
      end subroutine sobseq


      SUBROUTINE spearcorr(C,D)
      IMPLICIT NONE
      DOUBLE PRECISION, dimension(:,:), INTENT(out) :: C
      integer, dimension(:,:),intent(in) :: D ! rank matrix
      double precision, dimension(:,:),allocatable :: DD,DDT
      double precision, dimension(:),allocatable :: tmp
      INTEGER             :: N,M,ix,iy      
      DOUBLE PRECISION    :: dN      
! this procedure calculates spearmans correlation coefficient
! between the columns of D 
     
      N=size(D,dim=1);M=SIZE(D,dim=2)
      dN=dble(N)
      allocate(DD(1:N,1:M))
      DD=dble(D)
      if (.false.) then ! old call
         allocate(DDt(1:M,1:N))
         DDT=transpose(DD)
         C = matmul(DDt,DD)*12.d0/(dn*(dn*dn-1.d0)) 
         C=(C-3.d0*(dn+1.d0)/(dn-1.d0))
         deallocate(DDT)
      else
         allocate(tmp(1:N))
         do  ix=1, m-1
            do iy=ix+1,m
              tmp= DD(1:N,ix)-DD(1:N,iy)         
              C(ix,iy)=1.d0-6.d0*SUM(tmp*tmp)/dn/(dn*dn-1.d0)  
              C(iy,ix)=C(ix,iy)
           enddo
           C(ix,ix) = 1.d0
        enddo
        C(m,m)=1.d0
        deallocate(tmp)
      endif
      deallocate(DD)

      return
      END SUBROUTINE spearcorr
  

      SUBROUTINE swapRe(m,n)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(inout) :: m,n
      DOUBLE PRECISION                :: tmp      
      tmp=m
      m=n
      n=tmp
      END SUBROUTINE swapRe
  
      SUBROUTINE swapint(m,n)
      IMPLICIT NONE
      INTEGER, INTENT(inout) :: m,n
      INTEGER                :: tmp 
      tmp=m
      m=n
      n=tmp
      END SUBROUTINE swapint
      
      SUBROUTINE getdiag(diag,matrix)
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(:  ), INTENT(out) :: diag
      DOUBLE PRECISION, DIMENSION(:,:), INTENT(in)  :: matrix
      DOUBLE PRECISION, DIMENSION(:  ), ALLOCATABLE :: vector   

      ALLOCATE(vector(SIZE(matrix)))                    
      vector=PACK(matrix,.TRUE.)
      diag=vector(1:SIZE(matrix):SIZE(matrix,dim=1)+1)
      DEALLOCATE(vector)
      END SUBROUTINE getdiag

      SUBROUTINE SVDCMP(A,W,V)
      IMPLICIT NONE 
      DOUBLE PRECISION, DIMENSION(:  ), INTENT(out) :: W
      DOUBLE PRECISION, DIMENSION(:,:), INTENT(inout)  :: A
      DOUBLE PRECISION, DIMENSION(:,:), INTENT(OUT) :: V   
!LOCAL VARIABLES
      DOUBLE PRECISION, DIMENSION(:), allocatable :: RV1   
      DOUBLE PRECISION :: G,S,SCALE,ANORM,F,H,C,X,Y,Z   
      INTEGER M,N,NM,I,J,K,L,ITS
      
      !PARAMETER (NMAX=100)
C  Maximum anticipated values of  N

C  DIMENSION A(MP,NP),W(NP),V(NP,NP),RV1(NMAX)
C  Given a matrix  A, with logical dimensions  M  by  N  and physical
C  dimensions  MP  by  NP, this routine computes its singular value
C  decomposition,  A=U.W.V^T, see Numerical Recipes, by Press W.,H.
C  Flannery, B. P., Teukolsky S.A. and Vetterling W., T. Cambrige
C  University Press 1986, Chapter 2.9. The matrix  U  replaces A  on
C  output. The diagonal matrix of singular values  W  is output as a vector
C  W. The matrix  V (not the transpose  V^T) is output as  V.  M  must be
C  greater or equal to  N; if it is smaller, then  A  should be filled up
C  to square with zero rows.
C
      
       M=size(A,dim=1);N=size(A,dim=2)
       !Mp=M;Np=N
       allocate(RV1(1:N))
      IF (M.LT.N) then
         Print *,'SVDCMP: You must augment  A  with extra zero rows.'
      endif
C  Householder reduction to bidiagonal form
       G=0.d0
       SCALE=0.d0
       ANORM=0.d0
       DO 25 I=1,N
          L=I+1
          RV1(I)=SCALE*G
          G=0.D0
          S=0.D0
          SCALE=0.D0
          IF (I.LE.M) THEN
             DO  K=I,M
               SCALE=SCALE+ABS(A(K,I))
             enddo
             IF (SCALE.NE.0.D0) THEN
                DO  K=I,M
                  A(K,I)=A(K,I)/SCALE
                  S=S+A(K,I)*A(K,I)
                enddo
                F=A(I,I)
                G=-SIGN(SQRT(S),F)
                H=F*G-S
                A(I,I)=F-G
                IF (I.NE.N) THEN
                   DO  J=L,N
                     S=0.D0
                     DO  K=I,M
                       S=S+A(K,I)*A(K,J)
                     enddo
                     F=S/H
                     DO  K=I,M
                       A(K,J)=A(K,J)+F*A(K,I)
                    enddo
                enddo
              ENDIF
              DO  K=I,M
                 A(K,I)=SCALE*A(K,I)
              enddo
           ENDIF
       ENDIF
       W(I)=SCALE*G
       G=0.d0
       S=0.d0
       SCALE=0.d0
       IF ((I.LE.M).AND.(I.NE.N)) THEN
           DO  K=L,N
               SCALE=SCALE+ABS(A(I,K))
           enddo
             IF (SCALE.NE.0.0) THEN
                DO  K=L,N
                  A(I,K)=A(I,K)/SCALE
                  S=S+A(I,K)*A(I,K)
                enddo
                F=A(I,L)
                G=-SIGN(SQRT(S),F)
                H=F*G-S
                A(I,L)=F-G
                DO  K=L,N
                  RV1(K)=A(I,K)/H
                enddo
                IF (I.NE.M) THEN
                   DO  J=L,M
                     S=0.D0
                     DO  K=L,N
                       S=S+A(J,K)*A(I,K)
                    enddo
                     DO  K=L,N
                       A(J,K)=A(J,K)+S*RV1(K)
                     enddo
                   enddo
              ENDIF
              DO  K=L,N
                 A(I,K)=SCALE*A(I,K)
              enddo
           ENDIF
       ENDIF
       ANORM=MAX(ANORM,(ABS(W(I))+ABS(RV1(I))))
25     CONTINUE
c        print *,'25'
C   Accumulation of right-hand transformations.
       DO  I=N,1,-1
       IF (I.LT.N) THEN
         IF (G.NE.0.d0) THEN
           DO  J=L,N
             V(J,I)=(A(I,J)/A(I,L))/G
C   Double division to avoid possible underflow.
           enddo
          DO  J=L,N
            S=0.d0
            DO  K=L,N
              S=S+A(I,K)*V(K,J)
            enddo
            DO  K=L,N
              V(K,J)=V(K,J)+S*V(K,I)
            enddo
          enddo
        ENDIF
        DO  J=L,N
          V(I,J)=0.d0
          V(J,I)=0.d0
        enddo
       ENDIF
       V(I,I)=1.d0
       G=RV1(I)
       L=I
       enddo
c        print *,'32'

C  Accumulation of the left-hang transformation
       DO 39 I=N,1,-1
         L=I+1
         G=W(I)
         IF (I.LT.N) THEN
           DO  J=L,N
             A(I,J)=0.d0
           enddo
         ENDIF
         IF (G.NE.0.d0) THEN
           G=1.d0/G
           IF (I.NE.N) THEN
             DO  J=L,N
               S=0.d0
               DO K=L,M
                 S=S+A(K,I)*A(K,J)
               enddo
               F=(S/A(I,I))*G
             DO  K=I,M
               A(K,J)=A(K,J)+F*A(K,I)
             enddo
           enddo
         ENDIF
        DO  J=I,M
          A(J,I)=A(J,I)*G
        enddo
       ELSE
         DO  J=I,M
           A(J,I)=0.d0
         enddo
       ENDIF
       A(I,I)=A(I,I)+1.d0
39     CONTINUE
c        print *,'39'

C   Diagonalization of the bidiagonal form
C   Loop over singular values
       DO 49 K=N,1,-1
C   Loop allowed iterations
         DO 48 ITS=1,30
C   Test for spliting
            DO  L=K,1,-1
              NM=L-1
C   Note that RV1(1) is always zero
! old call which may cause inconsistent results
!              IF((ABS(RV1(L))+ANORM).EQ.ANORM) GO TO 2
!              IF((ABS(W(NM))+ANORM).EQ.ANORM) GO TO 1
! NEW CALL
              IF (((ABS(RV1(L))+ANORM).GE.NEAREST(ANORM,-1.d0)).AND.
     &          ((ABS(RV1(L))+ANORM).LE.NEAREST(ANORM,1.d0)) ) GO TO 2
              IF (((ABS(W(NM))+ANORM).GE.NEAREST(ANORM,-1.d0)).AND.
     &          ((ABS(W(NM))+ANORM).LE.NEAREST(ANORM,1.d0)) ) GO TO 1

            enddo
c          print *,'41'
1         C=0.d0
          S=1.d0
          DO  I=L,K
            F=S*RV1(I)
! old call which may cause inconsistent results

            IF (((ABS(F)+ANORM).LT.ANORM).OR.
     &            ((ABS(F)+ANORM).GT.ANORM)) THEN
              G=W(I)
              H=SQRT(F*F+G*G)
              W(I)=H
              H=1.D0/H
              C= (G*H)
              S=-(F*H)
              DO  J=1,M
                Y=A(J,NM)
                Z=A(J,I)
                A(J,NM)=(Y*C)+(Z*S)
                A(J,I)=-(Y*S)+(Z*C)
              enddo
            ENDIF
          enddo
c          print *,'43'
2         Z=W(K)
          IF (L.EQ.K) THEN
C   Convergence
            IF (Z.LT.0.d0) THEN
C   Singular values are made nonnegative
              W(K)=-Z
              DO  J=1,N
                V(J,K)=-V(J,K)
              enddo
            ENDIF
            GO TO 3
          ENDIF
          IF (ITS.EQ.30) then
             print *,'SVDCMP: No convergence in 30 iterations'
          endif
          X=W(L)
          NM=K-1
          Y=W(NM)
          G=RV1(NM)
          H=RV1(K)
          F=((Y-Z)*(Y+Z)+(G-H)*(G+H))/(2.d0*H*Y)
          G=SQRT(F*F+1.D0)
          F=((X-Z)*(X+Z)+H*((Y/(F+SIGN(G,F)))-H))/X
C   Next  QR  transformation
          C=1.d0
          S=1.d0
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
            DO  NM=1,N
              X=V(NM,J)
              Z=V(NM,I)
              V(NM,J)= (X*C)+(Z*S)
              V(NM,I)=-(X*S)+(Z*C)
            enddo
c            print *,'45',F,H
            Z=pythag(F,H)
            W(J)=Z
C   Rotation can be arbitrary if  Z=0.
            IF (Z.NE.0.d0) THEN
c            print *,1/Z
              Z=1.d0/Z
c              print *,'*'
              C=F*Z
              S=H*Z
            ENDIF
            F= (C*G)+(S*Y)
            X=-(S*G)+(C*Y)
            DO  NM=1,M
              Y=A(NM,J)
              Z=A(NM,I)
              A(NM,J)= (Y*C)+(Z*S)
              A(NM,I)=-(Y*S)+(Z*C)
            enddo
c          print *,'46'

47        CONTINUE
c          print *,'47'
          RV1(L)=0.D0
          RV1(K)=F
          W(K)=X
48      CONTINUE
3      CONTINUE
49     CONTINUE
c        print *,'49'
       deallocate(RV1)
       RETURN
       END SUBROUTINE SVDCMP

       FUNCTION pythag(a,b) RESULT (VALUE)
       DOUBLE PRECISION, INTENT(IN) :: a,b
       DOUBLE PRECISION :: VALUE
       DOUBLE PRECISION :: absa,absb
       absa=abs(a)
       absb=abs(b)
       IF (absa.GT.absb) THEN
          VALUE=absa*SQRT(1.d0+(absb/absa)**2)
       ELSE
          IF (absb.EQ.0) THEN
             VALUE=0.D0
          ELSE
             VALUE=absb*SQRT(1.d0+(absa/absb)**2)
          ENDIF
       ENDIF
       RETURN
       END FUNCTION PYTHAG
       

      FUNCTION THL(H1, L1, GH) RESULT (value)                          
      USE GLOBALDATA, ONLY : EPSS,PI,PI1,SQTWO,SQPI
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in) :: H1, L1, GH
      DOUBLE PRECISION             :: value                                     
!                                                                       
!     Local constants                                                   
!                                                                       
      DOUBLE PRECISION, PARAMETER :: ZERO=0.D0
      DOUBLE PRECISION, PARAMETER :: QUART=0.25D0
      DOUBLE PRECISION, PARAMETER :: HALF=0.5D0
      DOUBLE PRECISION, PARAMETER :: ONE=1.D0
      DOUBLE PRECISION, PARAMETER :: TWO=2.D0
      DOUBLE PRECISION, PARAMETER :: EXPLIM=80.D0
      DOUBLE PRECISION, DIMENSION(10) :: Ck
      PARAMETER ( Ck =(/ 0.99999995579043d0, -0.99999111387630d0, 
     &     0.99969933071259d0,-0.99596450623251d0, 0.97168814716964d0,  
     &    -0.88105640681133d0, 0.67507517897079d0,-0.38534334229836d0,
     &     0.13907128134639d0,-0.02317854687612d0/))
    
! Local variables
      INTEGER :: k
      DOUBLE PRECISION :: H,L, H2, T, CONEX, CN, 
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
!     The accuracy is specified with EPSS.  
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
!      print *, 'THL enter,L,H',L1,H1
      T=ZERO;L=L1;H=H1
      value=zero
      if (abs(L).LT.EPSS) THEN ! T(H,0)=0
         value=T
         return
      endif
      if (abs(H).LT.EPSS) THEN !T(0,L)
         !print *,'THL, atan'
         value=atan(L)*HALF*PI1
         !print *,'THL, atan',value
         return
      endif
      if (abs(abs(L)-ONE).LT.EPSS) THEN ! T(h,1)
         value=HALF*GH*(ONE-GH)*SIGN(ONE,L)
         return
      end if   
!      print *, 'THL still here'
      SGN=ONE
      if (abs(L).GT. ONE) THEN  ! make sure abs(L)<=1 to avoid numerical problems
                                ! using the identities:
                                ! T(h,L)=-T(h,-L)
                                ! T(h,L)=1/4-(FI(h)-1/2)*(FI(h*L)-1/2)-T(h*L,1/L)
         H = H * L                                                       
         L = ONE / L
                                                            
         if (L .LT. ZERO) THEN
            T = (HALF-GH)*(FI(H)-HALF)-QUART
         else  
            SGN=-ONE;
            T = (HALF-GH)*(FI(H)-HALF)+QUART
         endif
      endif

      H2 = H * H * HALF                                                                                                                                                                    
      if (H2 < EXPLIM) THEN 
         EX = exp(-H2)                                                      
      else
         value=T                ! T is found to the desired accuracy
!          print *, 'THL leaving'
         return
      end if 

      if ((H.LT. 1.6d0) .OR. (L.LE. 0.3d0)) then       
         L2 = L * L
         W2 = H2 * EX                                                       
         AP = ONE; SP=ONE
         S1 = ONE - EX;  S2=S1;    CN=S1
         CONEX = abs(PI*EPSS / L)                                              
         
         do while ((abs(CN).GT.CONEX).AND.( SP.LT.30.d0))
            SN =  SP                                                            
            SP =  SP + ONE                                                      
            S2 =  S2 - W2                                                       
            W2 =  W2 * H2 / SP                                                  
            AP = -AP * L2                                                        
            CN =  AP * S2 / (SN + SP)                                           
            S1 =  S1 + CN
         end do

!         print *, 'THL leaving 1.6'
!         print *,T,SGN,atan(L),S1,HALF,PI1
         
         value = T + SGN*(atan(L) - L*S1) *HALF*PI1         
!         print *,'value',value           
         return
      endif
      ! Otherwise
      ! approximate 1/(1+x^2) with a Tchebyshev polynomial of degree 2K
      ! and then integrate. Ck contains the coefficients of the polynomial.
      ! The absolute error is less than 1e-10.
        
      H2L=H2*L*L
      W2=exp(-H2L)*SQTWO/(H*L)
      Ik=SQPI*(FI(H*L)-HALF)
      SN=SQTWO/H
   
      CONEX = abs(PI*EPSS / EX)                                              
      CN=Ck(1)*Ik*SN
      S1=CN
      k=1
      do while ((abs(CN).GT.CONEX) .AND. (k.LT.10))
         W2=W2*H2L
         Ik=HALF*((TWO*DBLE(k)-ONE)*Ik-W2)
         SN=SN/H2
         k=k+1
         CN=Ck(k)*Ik*SN
         S1=S1+CN   
      end do	
                            
      T = T + SGN*EX*S1*HALF*PI1
!      print *, 'THL leaving last'
      RETURN                                                             
      END FUNCTION THL                                      

      function norm2dprb(a1,b1,a2,b2,R) RESULT (prb)
      USE GLOBALDATA, ONLY : EPSS,CFxCutoff,xCutoff,PI1
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in) :: a1,b1,a2,b2,R
      DOUBLE PRECISION             :: prb
      DOUBLE PRECISION :: L1,L2,L3,L4,Ga1,Ga2,Gb1,Gb2,SQR,
     &     T11,T12,T21,T22,T31,T32,T41,T42

!NORM2DPRB calculates the probability Prob(a1<X1<b1,a2<x2<b2)
!          for a bivariate Gaussian distribution
!          with correlation R
!
! References
!  Jagdish K. Patel and Campbell B. Read (1982)
! "Handbook of the normal distribution",
! marcel dekker inc,New York - Basel, Vol. 40, pp 293--300
! Tested on: Matlab 5.1
! NOTE!!! this does not work correctly yet. For some mysterious reason
! it crashes from time to time giving the message floating invalid

! History:
! by pab 19.08.1999

      if ((a1+EPSS.GE.b1).OR.(a2+EPSS.GE.b2)) THEN
         prb=0.d0
         return
      endif	
      

      Ga1=FI(a1)
      Gb1=FI(b1)
                                ! See if area is symmetric
                                ! The following may be optimized further
                                ! by checking against xCutOff
      if (abs(a1-a2).LT.EPSS) THEN
         Ga2=Ga1
      else	
         if (abs(a2-b1).LT.EPSS) THEN
            Ga2=Gb1
         else   
            Ga2=FI(a2)
         endif
      endif

      if (abs(b1-b2).LT.EPSS) THEN
         Gb2=Gb1
      else	
         if (abs(b2-a1).LT.EPSS) THEN
            Gb2=Ga1
         else   
            Gb2=FI(b2)
         endif	
      endif

      if (abs(R).LT.EPSS) then ! R=0
         prb=(Gb1-Ga1)*(Gb2-Ga2)
         return
      endif

      if (abs(R).GT.1.d0) THEN
         PRINT *,'Warning: correlation larger/smaller than +/- 1, R= ',R
         PRINT *,'limits (a1,b1)',a1,b1,'(a2,b2)',a2,b2
         if (R.GT.0.d0) then    ! R = 1
            prb=min(Gb1,Gb2)+min(Ga1,Ga2)-min(Ga1,Gb2)-min(Gb1,Ga2)
         else   ! R =-1
            prb=max(Gb1+Gb2-1.d0,0.d0)+max(Ga1+Ga2-1.d0,0.d0)-
     &           max(Ga1+Gb2-1.d0,0.d0)-max(Gb1+Ga2-1.d0,0.d0)
         endif
         return
      endif

      SQR=sqrt(1.d0-R*R)
      if (SQR.LT.EPSS) then
         if (R.GT.0.d0) then    ! R = 1
            prb=min(Gb1,Gb2)+min(Ga1,Ga2)-min(Ga1,Gb2)-min(Gb1,Ga2)
         else   ! R =-1
            prb=max(Gb1+Gb2-1.d0,0.d0)+max(Ga1+Ga2-1.d0,0.d0)-
     &           max(Ga1+Gb2-1.d0,0.d0)-max(Gb1+Ga2-1.d0,0.d0)
         endif
         return
      endif

      !print *, 'norm2d enter',ga1,gb1,ga2,gb2
      !print *, 'norm2d enter',a1,b1,a2,b2,R
                                ! L1=PSI(b1,b2,R),! L3=PSI(b1,a2,R)
      if (abs(b1).LT.EPSS) then
         if (abs(b2).LT.EPSS) then
            L1=0.25d0+asin(R)*0.5d0*PI1
         else
            T11=SIGN(0.25d0,b2)
            T12=THL(b2,-R/SQR,Gb2)
            L1=0.5d0*Gb2+0.25d0-T11-T12+min(0.d0,SIGN(0.5d0,b2))
         endif
         if (abs(a2).LT.EPSS) then
            L3=0.25d0+asin(R)*0.5d0*PI1
         else
            T31=SIGN(0.25d0,a2)
            T32=THL(a2,-R/SQR,Ga2)
            L3=0.5d0*Ga2+0.25d0-T31-T32+min(0.d0,SIGN(0.5d0,a2))
         endif
      else
         !print *, 'norm2d'
         if (abs(b2).LT.EPSS) then
            T11=THL(b1,-R/SQR,Gb1)
            T12=SIGN(0.25d0,b1)
            L1=0.5d0*Gb1+0.25d0-T11-T12+0.5d0*min(0.d0,SIGN(1.d0,b1))
         else
            !print *, 'norm2d,b1,r,gb1',b1,(b2/b1-R)/SQR,gb1
            T11=THL(b1,(b2/b1-R)/SQR,Gb1)
            !print *, 'norm2d,T11',T11
            if (abs(b2-b1).LT.EPSS) then
               L1=Gb1-2.d0*T11
               T12=T11
            else  
               !print *,'norm2d T12'
               if (abs(b2+b1).LT.EPSS) then
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
         if (abs(a2).LT.EPSS ) then
            T32=SIGN(0.25d0,b1)
            T31=THL(b1,-R/SQR,Gb1)
            L3=0.5d0*Gb1+0.25d0-T31-T32+min(0.d0,SIGN(0.5d0,b1))
         else
           ! print *, 'norm2d,b1,r,gb1',b1,(a2/b1-R)/SQR,gb1
            T31=THL(b1,(a2/b1-R)/SQR,Gb1)
            !print *, 'norm2d,T31',T31
            if (abs(a2-b1).LT.EPSS) THEN
               L3=Ga1-2.d0*T31
               T32=T31
            else  
               IF (abs(a2+b1).LT.EPSS) THEN
                  T32=T31
               ELSE
                  T32=THL(a2,(b1/a2-R)/SQR,Ga2)
               ENDIF
               if ((abs(a2+b2).LT.EPSS).AND.(abs(a1+b1).LT.EPSS)) then
                  prb=2.d0*(T31+T32-T11-T12)+1.d0 ! OK
                                !prb=2.d0*(L1-L3)-1.d0
                  goto 200
               endif
               L3=0.5d0*(Gb1+Ga2)-T31-T32+min(0.d0,SIGN(0.5d0,b1*a2))        
            endif	
         endif   
      endif
      !print *, 'norm2d L1,L3',L1,L3
                                !!L2=PSI(a1,a2,R) L4=PSI(a1,b2,R)
      if (abs(a1).LT.EPSS) then
         if (abs(b2).LT.EPSS) then
            L4=0.25d0+asin(R)*0.5d0*PI1
         else
            T41=SIGN(0.25d0,b2)
            T42=THL(b2,-R/SQR,Gb2)
            L4=0.5d0*Gb2+0.25d0-T41-T42+min(0.d0,SIGN(0.5d0,b2))
         endif
         if (abs(a2).LT.EPSS) then
            L2=0.25d0+asin(R)*0.5d0*PI1
         else
            T21=SIGN(0.25d0,a2)
            T22=THL(a2,-R/SQR,Ga2)
            L2=0.5d0*Ga2+0.25d0-T21-T22+min(0.d0,SIGN(0.5d0,a2))
         endif
      else
         if (abs(b2).LT.EPSS ) then
            T41=THL(a1,-R/SQR,Ga1)
            T42=SIGN(0.25d0,a1)
            L4=0.5d0*Ga1+0.25d0-T41-T42+min(0.d0,SIGN(0.5d0,a1))
         else
            !print *, 'norm2d T41, L',(b2/a1-R)/SQR
            T41=THL(a1,(b2/a1-R)/SQR,Ga1)
            if (abs(b2-a1).LT.EPSS) then
               L4=Gb1-2.d0*T41
               T42=T41
            else    
               if (abs(b2+a1).LT.EPSS) then
                  T42=T41
               else
                  T42=THL(b2,(a1/b2-R)/SQR,Gb2)
               endif
               L4=0.5d0*(Ga1+Gb2)-T41-T42+min(0.d0,SIGN(0.5d0,a1*b2))
            endif	
         endif
         !print *, 'norm2d L4',L4
         if (abs(a2).LT.EPSS) then
            T22=SIGN(0.25d0,a1)
            T21=THL(a2,-R/SQR,Ga1)
            L2=0.5d0*Ga1+0.25d0-T21-T22+min(0.d0,SIGN(0.5d0,a1))
         else
            T21=THL(a1,(a2/a1-R)/SQR,Ga1)
            if (abs(a2-a1).LT.EPSS) then
               L2=Ga1-2.d0*T21
               T22=T21
            else   
               if (abs(a2+a1).LT.EPSS) then
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
      else 
         if (prb<0.d0) then 
            prb=0.d0
         end if
      end if
      !print *, 'norm2d leaving'
      return
      END FUNCTION NORM2DPRB
      


      function MNORMPRB(BIG,Cm) RESULT (XIND)
      USE GLOBALDATA, ONLY : Hlo,Hup,xCutOff,NUGGET,EPSS,EPS2,
     &     CFxCutOff, RelEps,NSIMmax,NSIMmin,Nt,Nd,Nj,Ntd,SQ,
     &     Njj,Ntscis,NsXtmj, indXtd,rateLHD,index1,
     &     useC1C2,useMIDP,C1C2det,MLHD,COV,COVix
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(:  ), INTENT(in)    :: Cm  ! conditional mean
      DOUBLE PRECISION, DIMENSION(:,:), INTENT(inout) :: BIG ! conditional covariance 
      DOUBLE PRECISION                                :: XIND
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE  :: PP,V,C
      DOUBLE PRECISION, DIMENSION(:  ), ALLOCATABLE  :: CmN,D
      INTEGER,          DIMENSION(:,:), ALLOCATABLE  :: lhd

      DOUBLE PRECISION :: Pl1,Pu1,intval,intvalold,varsum,varsum2
      DOUBLE PRECISION :: ErrEst,ErrEst2, RelErrEst,RelErrEst2
      DOUBLE PRECISION :: Y,Y2,dL,dN,dNlhd,XMI,XMA,SQ0
      INTEGER          :: Nst,Nst0,L,K,M,Nlhd
      INTEGER          :: Nrep,Ndim

!MNORMPRB Multivariate Normal integrals by SCIS or LHSCIS
!  SCIS   = Sequential conditioned importance sampling
!  LHSCIS = Latin Hypercube Sequential Conditioned Importance Sampling
!
! !  NB!!: R must be conditional sorted by condsort3   
!        works on the upper triangular part of R
!
! References 
! R. Ambartzumian, A. Der Kiureghian, V. Ohanian and H.
! Sukiasian (1998)
! Probabilistic Engineering Mechanics, Vol. 13, No 4. pp 299-308
!
! Alan Genz (1992)
! 'Numerical Computation of Multivariate Normal Probabilities'
! J. computational Graphical Statistics, Vol.1, pp 141--149


      !print *,'enter mnormprb'
      Nst0 = NsXtmj(Njj+Ntscis)
      if (Njj.GT.0) then
         Nst  = NsXtmj(Njj)
      else
         Nst  = NsXtmj(Ntscis+1)
      endif
      !Nst=size(Cm)
      if (Nst.lt.Njj+1) then
        XIND=1.d0
        if (allocated(COV)) then ! save the coefficient of variation in COV
           COV(COVix)=0.d0
        endif
        return
      endif

      if (Nst.lt.Njj+1) then
         if (allocated(COV)) then ! save the coefficient of variation in COV
            COV(COVix)=0.d0
         endif
      
         XIND=1.d0
         return
      endif
       !print *,' mnormprb start calculat'
      intval=0.d0
      intvalold=0.d0
      Varsum=0.d0
      Varsum2=0.d0
            
      SQ0=SQ(Njj+1,Njj+1)
      XMA = (Hup (Njj+1)-Cm(Njj+1))/SQ0
      XMI = (Hlo (Njj+1)-Cm(Njj+1))/SQ0
      
      if (useC1C2) then         ! see if we can narrow down sampling range 
                                !            XMI=max(XMI,-xCutOff)
                                !            XMA=min(XMA,xCutOff)
         CALL C1C2(XMI,XMA,Cm(Njj+2:Nst),BIG(1,2:Nst),
     &        SQ(2:Nst,1),SQ0,indXtd(2:Nst))
         
         IF (XMA.LE.XMI) THEN
                                ! PQ= Y=0 for all return 
            goto 300
         ENDIF
      endif
      Pl1=FI(XMI)
      Pu1=FI(XMA)
      Ndim=Nst0-Njj
      Nlhd=max(1,rateLHD*Nst0)             ! size of LHD
      if (rateLHD.lt.1) then ! make sure 
         usemidp=.FALSE.
         print * ,'mnormprb: only able to use useMIDP if rateLHD>0'
      endif 
      ALLOCATE(CmN(1:Nst))
!      CmN(Nst:Nsd)=0.d0
      allocate(pp(1:Nlhd,1:Ndim))
      if (nlhd.GT.1) then
         allocate(lhd(1:Nlhd,1:Ndim)) ! allocate LHD  
         allocate(C(1:Ndim,1:Ndim))    
         if (MLHD) then
            allocate(D(1:Ndim))
            allocate(V(1:Ndim,1:Ndim))
         endif
      end if
 
      dNlhd=dble(Nlhd)
      Nrep=ceiling(dble(NSIMmax)/dNlhd) ! # replications
      dL=0.d0
      dN=0.d0     
      !print *,' mnormprb start L loop'
      DO L=1,Nrep
         CALL random_number(PP)
         if (Nlhd.gt.1) then
            LHD(:,1)=(/ (K,K=1,Nlhd)/)
            do M=1,3 ! do 3 attemps to construct a LHD with rank=Ndim
               do K=2,Ndim
                  CALL sortre(LHD(:,K),PP(:,K)) ! lhd = latin hypercube design
               enddo
               CALL spearcorr(C,lhd) ! find rankcorrelation between columns
               do K=1,Ndim-1 ! see if rank=Ndim
                  if (any(abs(C(K,K+1:Ndim)).GE.1.d0))  then
                     if (M.EQ.3) goto 30
                     CALL random_number(PP)
                     cycle
                  endif
               enddo
               goto 20
            enddo
 20         if (MLHD) then   !modify lhd by reducing correlation between columns
               DO K=1,Ndim
                  C(K,K)=C(K,K)+NUGGET ! add nugget effect to ensure that 
                                !inversion is not corrupted by round off errors
               enddo
               CALL svdcmp(C,D,V) ! C=U*D*V'=Q*Q'
               do K=1,Ndim
                 V(K,:)=C(:,K)*sqrt(1/D(K)) ! inverting Q=U*sqrt(D)
               enddo
              
               PP=MATMUL(dble(lhd),V)  ! LHD*inv(Q)
               do K=1,Ndim
                  CALL sortre(LHD(:,K),PP(:,K)) ! lhd = latin hypercube design
               enddo
            endif 
 30         if (USEMIDP) then   ! use the center of the cell
               PP=(dble(lhd)-0.5d0)/dNlhd   
            else                ! distribute uniformly within the cell
               CALL random_number(PP)
               PP=(dble(lhd)-PP)/dNlhd 
            endif
         endif

           
         !print *,' mnormprb start M loop'
         DO M=1,Nlhd
            CmN(1:Nst-Njj)=Cm(Njj+1:Nst) ! initialize conditional mean  
            Y=MVNFUN2(PP(M,:),BIG,CmN,Pl1,Pu1) ! evaluate the integrand
            CmN(1:Nst)=Cm(1:Nst) ! initialize conditional mean
            PP(M,:)=1.d0-PP(M,:) ! using antithetic variables to reduce variance
            Y= (Y+MVNFUN2(PP(M,:),BIG,CmN,Pl1,Pu1))/2.d0   
      
            dN=dN+1.d0
            varsum= varsum+(dN-1.d0)*(Y-intval)*(Y-intval)/dN 
            intval=intval+(Y-intval)/dN ! integral value
            ErrEst=2.5d0*SQRT(varsum/(dN*dN)) ! error estimate
            if (abs(intval).gt.0.d0) then
               RelErrEst = ErrEst/abs(intval)
            else
               RelErrEst=0.d0
            end if
                       
            IF (((RelErrEst.LT.RelEps).or.(ErrEst.lt.EPSS))
     &           .AND.(dN.GT.NSIMmin))  GOTO 300
      ENDDO                     ! M loop
      if (Nlhd.gt.1) then
            dL=dble(L)
            Y2=dL*(intval-intvalold)+intvalold
            varsum2= varsum2+(dL-1.d0)*(Y2-intval)*(Y2-intval)/dL ! Better estimate of variance
            ErrEst2=2.5d0*SQRT(varsum2/(dL*dL)) ! error estimate between replicates
            if (abs(intval).gt.0.d0) then
               RelErrEst2 = ErrEst/abs(intval)
            else
               RelErrEst2=0.d0
            end if
                       
            IF (((RelErrEst2.LT.RelEps).or.(ErrEst2.lt.EPSS))
     &           .AND.(L.GT.5))  GOTO 300
            intvalold=intval
         end if
       enddo                     ! L loop
 300  XIND=intval
      !print *,'mnormprb L,N,PQ,CV,CV2',L,dN,intval,ErrEst,ErrEst2
      !print *,'leaving mnormprb'
      if (allocated(COV)) then ! save the coefficient of variation in COV
         if ((dL.gt.1).and.(Nlhd.gt.1)) then
            COV(COVix)=min(RelErrEst,RelErrEst2)/2.5d0
         else
            COV(COVix)=RelErrEst/2.5d0
         endif
      endif
      if (allocated(PP)) then
         DEALLOCATE(CmN)
         DEallocate(pp)
      endif
      if (allocated(lhd)) then
         DEallocate(lhd)
         DEallocate(C)
         if (allocated(D)) then
            DEallocate(D)
            DEallocate(V)
         endif
      endif
      !print *,'leaving mnormprb'
      return
      END FUNCTION MNORMPRB

      DOUBLE PRECISION FUNCTION FIINV(P)
      use GLOBALDATA, only: XMAX,CFxCutOff
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
      DOUBLE PRECISION SPLIT1, SPLIT2, CONST1, CONST2, 
     *     A0, A1, A2, A3, A4, A5, A6, A7, B1, B2, B3, B4, B5, B6, B7, 
     *     C0, C1, C2, C3, C4, C5, C6, C7, D1, D2, D3, D4, D5, D6, D7, 
     *     E0, E1, E2, E3, E4, E5, E6, E7, F1, F2, F3, F4, F5, F6, F7, 
     *     P, Q, R
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
      Q = ( P - .5d0)
      IF ( ABS(Q) .LE. SPLIT1 ) THEN ! Central range.
         R = CONST1 - Q*Q
         FIINV = Q*( ( ( ((((A7*R + A6)*R + A5)*R + A4)*R + A3)
     *                  *R + A2 )*R + A1 )*R + A0 )
     *            /( ( ( ((((B7*R + B6)*R + B5)*R + B4)*R + B3)
     *                  *R + B2 )*R + B1 )*R + 1.d0 )
      ELSE ! near the endpoints
         R = MIN( P, 1.d0 - P )
         IF ( R .GT. CFxCutOff/2.d0) THEN ! R .GT.0.d0
            R = SQRT( -LOG(R) )
            IF ( R .LE. SPLIT2 ) THEN
               R = R - CONST2
               FIINV = ( ( ( ((((C7*R + C6)*R + C5)*R + C4)*R + C3)
     *                      *R + C2 )*R + C1 )*R + C0 ) 
     *                /( ( ( ((((D7*R + D6)*R + D5)*R + D4)*R + D3)
     *                      *R + D2 )*R + D1 )*R + 1.d0 )
            ELSE
               R = R - SPLIT2
               FIINV = ( ( ( ((((E7*R + E6)*R + E5)*R + E4)*R + E3)
     *                      *R + E2 )*R + E1 )*R + E0 )
     *                /( ( ( ((((F7*R + F6)*R + F5)*R + F4)*R + F3)
     *                      *R + F2 )*R + F1 )*R + 1.d0 )
            END IF
         ELSE
            FIINV = XMAX !9.d0
         END IF
         IF ( Q .LT. 0.d0 ) FIINV = - FIINV
      END IF
      END FUNCTION FIINV     



      function FIINVS(p1) RESULT(value)
      use GLOBALDATA, only: SQTWO,XMAX,CFxCutOff
      DOUBLE PRECISION, INTENT(IN) :: P1
      DOUBLE PRECISION :: VALUE
      DOUBLE PRECISION, DIMENSION(4) :: A,B,C
      DOUBLE PRECISION, DIMENSION(2) :: D
      DOUBLE PRECISION :: P,Z,X=0.d0
      DOUBLE PRECISION, PARAMETER :: P0=0.7D0
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
      END FUNCTION FIINVS


                                ! *********************************     
      FUNCTION FI(X) RESULT(value)
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
300   IF (XX .GT. xbreak) THEN
         value=1.d0-value
      END IF
      RETURN
                                !if (value > 1d0) then
                                !   value=1.d0
                                !else 
                                ! if (value<0.d0) then 
                                !    value=0.d0
                                ! end if
                                !end if
           
      END FUNCTION FI

      FUNCTION GAUSINT (X1, X2, A, B, C, D) RESULT (value)               
      USE GLOBALDATA,ONLY:    SQTWOPI1,xCutOff
      IMPLICIT NONE
      DOUBLE PRECISION             :: value,Y1,Y2,Y3
      DOUBLE PRECISION, INTENT(in) :: X1,X2,A,B,C,D
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
      USE GLOBALDATA,ONLY:    SQTWOPI1,xCutOff
      IMPLICIT NONE
      DOUBLE PRECISION             :: value,X0,Y0,Y1,Y2
      DOUBLE PRECISION, INTENT(in) :: X1,X2,A,B
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
         value=ABS(A)*(FI(X2)-FI(X1))
         RETURN
      ENDIF
     
 
      Y1=-A*FI(X1)+SQTWOPI1*B*EXP(-0.5d0*X1*X1)
      Y2=A*FI(X2)-SQTWOPI1*B*EXP(-0.5d0*X2*X2)
      IF ((B*X1.LT.-A).AND.(-A.LT.B*X2))THEN
          X0=-A/B  
         Y0=2*(A*FI(X0)-SQTWOPI1*B*EXP(-0.5d0*X0*X0))
         value=ABS(Y2-Y1-Y0)
      ELSE
         value=ABS(Y1+Y2) 
      ENDIF                                       
      RETURN                                                             
      END FUNCTION GAUSINT2
      
      
      
                                                                       
      SUBROUTINE GAUSSLE1 (N,WFout,BPOUT,XMI,XMA)
      USE GLOBALDATA,ONLY : SQTWOPI1,EPS0,xCutOff2
      USE QUAD , ONLY: LeBP,LeWF,LeIND,NLeW,minQnr
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(:),INTENT(out) :: BPOUT, WFout
      DOUBLE PRECISION,               INTENT(in) :: XMI,XMA
      INTEGER,                     INTENT(inout) :: N
      ! local variables
      DOUBLE PRECISION                           :: Z1,SDOT, SDOT1,DIFF1 
      INTEGER                                    :: NN,I,J,k    

      ! The subroutine picks the lowest Gauss-Legendre
      ! quadrature needed to integrate the test function
      ! gaussint to the specified accuracy, EPS0.  
      ! The nodes and weights between the  integration
      !  limits XMI and XMA (all normalized) are returned.  
      ! Note that the weights are multiplied with
      ! 1/sqrt(2*pi)*exp(.5*bpout^2) 

      IF (XMA.LE.XMI) THEN                                               
         PRINT * , 'Warning XMIN>=XMAX in GAUSSLE1 !',XMI,XMA                   
         RETURN
      ENDIF
             
      DO  I = minQnr, NLeW 
         NN = N  !initialize
         DO J = LeIND(I)+1, LeIND(I+1)
            BPOUT (NN+1) = 0.5d0*(LeBP(J)*(XMA-XMI)+XMA+XMI)         
            Z1 = BPOUT (NN+1) * BPOUT (NN+1)
                                !IF (Z1.LE.xCutOff2) THEN
            NN=NN+1
            WFout (NN) = 0.5d0 * SQTWOPI1 * (XMA - XMI) *  
     &              LeWF(J) *EXP ( - 0.5d0* Z1  )  
                                !ENDIF
         ENDDO
         
         SDOT = GAUSINT (XMI, XMA, - 2.5d0, 2.d0, 2.5d0, 2.d0)
         SDOT1 = 0.d0                                                      
       
         DO  k = N+1, NN                                                 
            SDOT1 = SDOT1+WFout(k)*(-2.5d0+2.d0*BPOUT(k) )* 
     &           (2.5d0 + 2.d0 * BPOUT (k) ) 
         ENDDO
         DIFF1 = ABS (SDOT - SDOT1)                                      
                            
         IF (EPS0.GT.DIFF1) THEN
            N=NN
            PRINT * ,'gaussle1, XMI,XMA,NN',XMI,XMA,NN
            RETURN                                      
         END IF                                                          
      END DO
      RETURN
      END SUBROUTINE GAUSSLE1
                                                                        
                                                                        
                                !********************************************* 
      SUBROUTINE GAUSSLE0 (N, wfout, bpout, XMI, XMA, N0) 
      USE GLOBALDATA, ONLY : SQTWOPI1,xCutOff2,EPSS
      USE QUAD,        ONLY : LeBP,LeWF,NLeW,LeIND
      IMPLICIT NONE
      INTEGER, INTENT(in)                         :: N0 
      INTEGER, INTENT(inout)                      :: N
      DOUBLE PRECISION, DIMENSION(:), INTENT(out) :: wfout,bpout 
      DOUBLE PRECISION,                INTENT(in) :: XMI,XMA    
! Local variables
      DOUBLE PRECISION                            :: Z1  
      INTEGER                                     :: J
      ! The subroutine computes Gauss-Legendre
      ! nodes and weights between
      ! the (normalized) integration limits XMI and XMA    
      ! Note that the weights are multiplied with
      ! 1/sqrt(2*pi)*exp(.5*bpout^2) so that
      !  b
      ! int f(x)*exp(-x^2/2)/sqrt(2*pi)dx=sum f(bp(j))*wf(j)
      !  a                                 j

      IF (XMA.LE.XMI) THEN                                               
         PRINT * , 'Warning XMIN>=XMAX in GAUSSLE0 !',XMI,XMA                    
         RETURN  ! no more nodes added
      ENDIF
      IF ((XMA-XMI).LT.EPSS) THEN
         N=N+1
         BPout (N) = 0.5d0 * (XMA + XMI)
         Z1 = BPOUT (N) * BPOUT (N)   
         WFout (N) = SQTWOPI1 * (XMA - XMI) *EXP ( - 0.5d0* Z1  )  
         RETURN
      ENDIF
      IF (N0.GT.NLeW) THEN                                                 
         PRINT * , 'error in GAUSSLE0, quadrature not available' 
         STOP                                                            
      ENDIF 
      !print *, 'GAUSSLE0',N0          
                 
      !print *, N                           
      DO  J = LeIND(N0)+1, LeIND(N0+1) 
             
         BPout (N+1) = 0.5d0 * (LeBP(J) * (XMA - XMI) + XMA + XMI)
         Z1 = BPOUT (N+1) * BPOUT (N+1)                                      
                                !         IF (Z1.LE.xCutOff2) THEN
         N=N+1                  ! add a new node and weight
         WFout (N) = 0.5d0 * SQTWOPI1 * (XMA - XMI) *  
     &        LeWF(J) *EXP ( - 0.5d0* Z1  )  
                                !         ENDIF  
      ENDDO
       !print *,BPout
      RETURN                                                             
      END SUBROUTINE GAUSSLE0
                          !********************************************* 
      SUBROUTINE GAUSSLE2 (N, wfout, bpout, XMI, XMA, N0) 
      USE GLOBALDATA, ONLY : SQTWOPI1,xCutOff,xCutOff2,EPSS
      USE QUAD,        ONLY : LeBP,LeWF,NLeW,LeIND,minQNr
      IMPLICIT NONE
      INTEGER, INTENT(in)                         :: N0 
      INTEGER, INTENT(inout)                      :: N
      DOUBLE PRECISION, DIMENSION(:), INTENT(out) :: wfout,bpout 
      DOUBLE PRECISION,                INTENT(in) :: XMI,XMA    
! Local variables
      DOUBLE PRECISION                            :: Z1  
      INTEGER                                     :: J,N1
      ! The subroutine computes Gauss-Legendre
      ! nodes and weights between
      ! the (normalized) integration limits XMI and XMA
      ! This procedure select number of nodes
      ! depending on the length of the integration interval.    
      ! Note that the weights are multiplied with
      ! 1/sqrt(2*pi)*exp(.5*bpout^2) so that
      !  b
      ! int f(x)*exp(-x^2/2)/sqrt(2*pi)dx=sum f(bp(j))*wf(j)
      !  a                                 j

      IF (XMA.LE.XMI) THEN                                               
         PRINT * , 'Warning XMIN>=XMAX in GAUSSLE2 !',XMI,XMA                    
         RETURN  ! no more nodes added
      ENDIF
      IF ((XMA-XMI).LT.EPSS) THEN
         N=N+1
         BPout (N) = 0.5d0 * (XMA + XMI)
         Z1 = BPOUT (N) * BPOUT (N)   
         WFout (N) = SQTWOPI1 * (XMA - XMI) *EXP ( - 0.5d0* Z1  )  
         RETURN
      ENDIF
      IF (N0.GT.NLeW) THEN                                                 
         PRINT * , 'Warning in GAUSSLE2, quadrature not available' 
      ENDIF 
      !print *, 'GAUSSLE0',N0          
                 
      !print *, N                           
      N1=CEILING(0.65d0*(XMA-XMI)*DBLE(N0)/xCutOff)
      N1=MAX(MIN(N1,NLew),minQNr)
      
      DO  J = LeIND(N1)+1, LeIND(N1+1) 
             
         BPout (N+1) = 0.5d0 * (LeBP(J) * (XMA - XMI) + XMA + XMI)
         Z1 = BPOUT (N+1) * BPOUT (N+1)                                      
                                !         IF (Z1.LE.xCutOff2) THEN
         N=N+1                  ! add a new node and weight
         WFout (N) = 0.5d0 * SQTWOPI1 * (XMA - XMI) *  
     &        LeWF(J) *EXP ( - 0.5d0* Z1  )  
                                !         ENDIF  
      ENDDO
      !PRINT * ,'gaussle2, XMI,XMA,N',XMI,XMA,N
       !print *,BPout
      RETURN                                                             
      END SUBROUTINE GAUSSLE2

      SUBROUTINE GAUSSHE0 (N, WFout, BPout, XMI, XMA, N0)
      USE GLOBALDATA, ONLY :  SQPI1,SQTWO               
      USE QUAD, ONLY : HeBP,HeWF,HeIND,NHeW
      IMPLICIT NONE
      INTEGER,                         INTENT(in) :: N0 
      INTEGER,                      INTENT(inout) :: N
      DOUBLE PRECISION, DIMENSION(:), INTENT(out) :: wfout,bpout 
      DOUBLE PRECISION,                INTENT(in) :: XMI,XMA  
      INTEGER                                     :: J
      ! The subroutine returns modified Gauss-Hermite
      ! nodes and weights between
      ! the integration limits XMI and XMA        
      ! for the chosen number of nodes
      ! implicitly assuming that the integrand
      !  goes smoothly towards zero as its approach XMI or XMA
      ! Note that the nodes and weights are modified
      ! according to
      !  Inf
      ! int f(x)*exp(-x^2/2)/sqrt(2*pi)dx=sum f(bp(j))*wf(j)
      ! -Inf                               j
       
      IF (XMA.LE.XMI) THEN                                               
         PRINT * , 'Warning XMIN>=XMAX in GAUSSHE0 !',XMI,XMA     
         RETURN ! no more nodes added
      ENDIF
      IF (N0.GT.NHeW) THEN                                                 
         PRINT * , 'error in GAUSSHE0, quadrature not available'
         STOP                                                            
      ENDIF 
   
      DO  J = HeIND(N0)+1, HeIND(N0+1)
         BPout (N+1) = HeBP (J) * SQTWO  
         IF (BPout (N+1).GT.XMA) THEN
            RETURN
         END IF   
         IF (BPout (N+1).GE.XMI) THEN
            N=N+1  ! add the node
            WFout (N) = HeWF (J) * SQPI1   
         END IF
      ENDDO
      RETURN                                                             
      END SUBROUTINE GAUSSHE0   

      SUBROUTINE GAUSSLA0 (N, WFout, BPout, XMI, XMA, N0)
      USE GLOBALDATA, ONLY :  SQPI1              
      USE QUAD, ONLY : LaBP5,LaWF5,LaIND,NLaW
      IMPLICIT NONE
      INTEGER,                         INTENT(in) :: N0 
      INTEGER,                      INTENT(inout) :: N
      DOUBLE PRECISION, DIMENSION(:), INTENT(out) :: wfout,bpout 
      DOUBLE PRECISION,                INTENT(in) :: XMI,XMA  
      INTEGER                                     :: J
      ! The subroutine returns modified Gauss-Laguerre
      ! nodes and weights for alpha=-0.5 between
      ! the integration limits XMI and XMA        
      ! for the chosen number of nodes
      ! implicitly assuming the integrand
      !  goes smoothly towards zero as its approach XMI or XMA
      ! Note that the nodes and weights are modified
      ! according to
      !  Inf
      ! int f(x)*exp(-x^2/2)/sqrt(2*pi)dx=sum f(bp(j))*wf(j)
      !  0                                 j
       
      IF (XMA.LE.XMI) THEN                                               
         PRINT * , 'Warning XMIN>=XMAX in GAUSSLA0 !',XMI,XMA     
         RETURN !no more nodes added
      ENDIF
      IF (N0.GT.NLaW) THEN                                                 
         PRINT * , 'error in GAUSSLA0, quadrature not available'
         STOP                                                            
      ENDIF 
   
      DO  J = LaIND(N0)+1, LaIND(N0+1)
         IF (XMA.LE.0.d0) THEN
            BPout (N+1) = -SQRT(2.d0*LaBP5(J))
         ELSE
            BPout (N+1) = SQRT(2.d0*LaBP5(J))   
         END IF
         IF (BPout (N+1).GT.XMA) THEN
            RETURN
         END IF   
         IF (BPout (N+1).GE.XMI) THEN
            N=N+1 ! add the node       
            WFout (N) = LaWF5 (J)*0.5d0*SQPI1
         END IF
      ENDDO
      !PRINT *,'gaussla0, bp',LaBP5(LaIND(N0)+1:LaIND(N0+1))
      !PRINT *,'gaussla0, wf',LaWF5(LaIND(N0)+1:LaIND(N0+1))
      RETURN                                                             
      END SUBROUTINE GAUSSLA0 

      SUBROUTINE GAUSSQ(N, WF, BP, XMI, XMA, N0) 
      USE GLOBALDATA, ONLY : xCutOff
      USE QUAD       , ONLY : minQNr
      IMPLICIT NONE
      INTEGER,                         INTENT(in) :: N0 
      INTEGER,                      INTENT(inout) :: N
      DOUBLE PRECISION, DIMENSION(:), INTENT(out) :: wf,bp 
      DOUBLE PRECISION,                INTENT(in) :: XMI,XMA  
      INTEGER                                     :: N1
      ! The subroutine returns 
      ! nodes and weights between
      ! the integration limits XMI and XMA        
      ! for the chosen number of nodes
      ! Note that the nodes and weights are modified
      ! according to
      !  Inf
      ! int f(x)*exp(-x^2/2)/sqrt(2*pi)dx=sum f(bp(j))*wf(j)
      !  0                                 j
       
      !IF (XMA.LE.XMI) THEN                                               
      !   PRINT * , 'Warning XMIN>=XMAX in GAUSSQ !',XMI,XMA     
      !   RETURN !no more nodes added
      !ENDIF 
      CALL GAUSSLE0(N,WF,BP,XMI,XMA,N0) 
      RETURN 
      IF ((XMA.GE.xCutOff).AND.(XMI.LE.-xCutOff)) THEN
         CALL GAUSSHE0(N,WF,BP,XMI,XMA,N0) 
      ELSE
         CALL GAUSSLE2(N,WF,BP,XMI,XMA,N0) 
         RETURN
         IF (((XMA.LT.xCutOff).AND.(XMI.GT.-xCutOff)).OR.(.TRUE.)
     &        .OR.(XMI.GT.0.d0).OR.(XMA.LT.0.d0)) THEN 
                                ! Grid by Gauss-LegENDre quadrature
            CALL GAUSSLE2(N,WF,BP,XMI,XMA,N0) 
         ELSE
            ! this does not work well
            PRINT *,'N0',N0,N
            N1=CEILING(DBLE(N0)/2.d0)
            IF (XMA.GE.xCutOff) THEN
               IF (XMI.LT.0.d0) THEN
                  CALL GAUSSLE2 (N, WF, BP,XMI ,0.d0,N0)
               ENDIF
               CALL GAUSSLA0 (N, WF, BP,0.d0, XMA, N1)
            ELSE
               IF (XMA.GT.0.d0) THEN
                  CALL GAUSSLE2 (N, WF,BP,0.d0,XMA,N0)
               ENDIF
               CALL GAUSSLA0 (N, WF,BP,XMI,0.d0, N1) 
            END IF
         END IF
      ENDIF
      !PRINT *,'gaussq, wf',wf(1:N)
      !PRINT *,'gaussq, bp',bp(1:N)
      RETURN
      END SUBROUTINE GAUSSQ
      END MODULE RIND         !******************************  


