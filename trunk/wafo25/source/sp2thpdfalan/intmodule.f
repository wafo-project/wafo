! INTMODULE contains the modules:
!        - ADAPTMOD
!        - RCRUDEMOD
!        - KROBOVMOD
!        - KRBVRCMOD
! which contains several different Multidimensional Integration Subroutines
!
! See descriptions below
!
*  ADAPTMOD is a module  containing a:
* 
*   Adaptive Multidimensional Integration Subroutine
*
*   Author: Alan Genz
*           Department of Mathematics
*           Washington State University
*           Pullman, WA 99164-3113 USA
*
* Revised pab 21.11.2000
*  A bug found by Igor in dksmrc: VK was not correctly randomized
*  is now fixed
* Revised pab 07.10.2000, 
*    1) Removed LENWRK and WORK from input in ADAPT. 
*    2) Defined LENWRK internally and Put a save statement before WORK instead 
*    3) Bug fix in ADBASE: DIVAXN was undetermined when MINCLS<0. Solution:
*         put a save statement on DIVAXN in order to save/keep its last value.
*    4) MAXDIM is now a global variable defining the maximum number of dimensions
*       it is possible to integrate.
*       
* revised pab 07.09.2000
*  - solaris compiler complained on the DATA statements
*    for the P and C matrices in the krbvrc and krobov routine
*    => made separate DATA statements for P and C and moved them
*      to right after the variable definitions.      
* revised pab 10.03.2000
*   - updated to f90 (i.e. changed to assumed shape arrays + changing integers to DBLE)
*   - put it into a module
*
*  This subroutine computes an approximation to the integral
*
*      1 1     1
*     I I ... I       FUNCTN(NDIM,X)  dx(NDIM)...dx(2)dx(1)
*      0 0     0  
*
***************  Parameters for SADAPT  ********************************
*
********Input  Parameters
*
*     N      INTEGER, the number of variables.
*     MAXPTS INTEGER, maximum number of function values allowed. This 
*            parameter can be used to limit the time taken. A 
*            sensible strategy is to start with MAXPTS = 1000*N, and then
*            increase MAXPTS if ERROR is too large.
*    FUNCTN  Externally declared real user defined integrand. Its 
*            parameters must be (N, Z), where Z is a real array of
*            length N.
*     ABSEPS REAL absolute error tolerance.
*     RELEPS REAL relative error tolerance.
*
*******Output  Parameters
*
*     ERROR  REAL estimated absolute error, with 99% confidence level.
*     VALUE  REAL estimated value for the integral
*     INFORM INTEGER, termination status parameter:
*            if INFORM = 0, normal completion with ERROR < EPS;
*            if INFORM = 1, completion with ERROR > EPS and MAXPTS 
*                           function vaules used; increase MAXPTS to 
*                           decrease ERROR;
*            if INFORM = 2, N > 20 or N < 1.
*
***************  Parameters for ADAPT  ********************************
*
****** Input Parameters
*
*  NDIM    Integer number of integration variables.
*  MINCLS  Integer minimum number of FUNCTN calls to be allowed; MINCLS
*          must not exceed MAXCLS. If MINCLS < 0, then ADAPT assumes
*          that a previous call of ADAPT has been made with the same
*          integrand and continues that calculation.
*  MAXCLS  Integer maximum number of FUNCTN calls to be used; MAXCLS
*          must be >= RULCLS, the number of function calls required for
*          one application of the basic integration rule.
*           IF ( NDIM .EQ. 1 ) THEN
*              RULCLS = 11
*           ELSE IF ( NDIM .LT. 15 ) THEN
*              RULCLS = 2**NDIM + 2*NDIM*(NDIM+3) + 1
*           ELSE
*              RULCLS = 1 + NDIM*(24-NDIM*(6-NDIM*4))/3
*           ENDIF
*  FUNCTN  Externally declared real user defined integrand. Its 
*          parameters must be (NDIM, Z), where Z is a real array of
*          length NDIM.
*  ABSREQ  Real required absolute accuracy.
*  RELREQ  Real required relative accuracy.
*
****** Output Parameters
*
*  MINCLS  Actual number of FUNCTN calls used by ADAPT.
*  ABSEST  Real estimated absolute accuracy.
*  FINEST  Real estimated value of integral.
*  INFORM  INFORM = 0 for normal exit, when ABSEST <= ABSREQ or
*                     ABSEST <= |FINEST|*RELREQ with MINCLS <= MAXCLS.
*          INFORM = 1 if MAXCLS was too small for ADAPT to obtain the
*                     result FINEST to within the requested accuracy.
*          INFORM = 2 if MINCLS > MAXCLS, LENWRK < 16*NDIM + 27 or 
*                     RULCLS > MAXCLS.
*
*
*
* ADAPT revised by pab 07.10.2000, 
*    1) Removed LENWRK and WORK from input. 
*    2) Defined LENWRK internally and Put a save statement before WORK instead 
*
*  WORK    Real array (length LENWRK) of working storage. This contains
*          information that is needed for additional calls of ADAPT
*          using the same integrand (input MINCLS < 0).
*  LENWRK  Integer length of real array WORK (working storage); ADAPT
*          needs LENWRK >= 16*NDIM + 27. For maximum efficiency LENWRK
*          should be about 2*NDIM*MAXCLS/RULCLS if MAXCLS FUNCTN
*          calls are needed. If LENWRK is significantly less than this,
*          ADAPT may be less efficient.
      MODULE ADAPTMOD
      IMPLICIT NONE
	INTEGER,PRIVATE, PARAMETER :: MAXDIM=20
      PRIVATE
      PUBLIC :: ADAPT, SADAPT

      INTERFACE SADAPT
      MODULE PROCEDURE SADAPT
      END INTERFACE
 
      INTERFACE ADAPT
      MODULE PROCEDURE ADAPT
      END INTERFACE
      
      INTERFACE ADBASE
      MODULE PROCEDURE ADBASE
      END INTERFACE
       
      INTERFACE  BSINIT
      MODULE PROCEDURE BSINIT
      END INTERFACE

      INTERFACE  RULNRM
      MODULE PROCEDURE RULNRM
      END INTERFACE

      INTERFACE  DIFFER
      MODULE PROCEDURE DIFFER
      END INTERFACE

      INTERFACE  BASRUL
      MODULE PROCEDURE BASRUL
      END INTERFACE

      INTERFACE FULSUM
      MODULE  PROCEDURE FULSUM
      END INTERFACE
      
      INTERFACE  TRESTR
      MODULE PROCEDURE TRESTR
      END INTERFACE
                                !--------------------------------
      CONTAINS   

!***********************************************************
!    MAIN INTEGRATION ROUTINE SADAPT
!***********************************************************  

      SUBROUTINE SADAPT(N,MAXPTS,FUNCTN,ABSEPS,
     &     RELEPS,ERROR,VALUE,INFORM)
      IMPLICIT NONE
*
*     A subroutine for computing multivariate integrals 
*     This subroutine uses an algorithm given in the paper
*     "Numerical Computation of Multivariate Normal Probabilities", in
*     J. of Computational and Graphical Stat., 1(1992), pp. 141-149, by
*          Alan Genz 
*          Department of Mathematics
*          Washington State University 
*          Pullman, WA 99164-3113
*          Email : alangenz@wsu.edu
*
* revised pab 15.03.2000
*   - changed name from SADMVN to SADAPT
*   - Made it general for any integral not just the multivariate normal integral
*
********Input  Parameters
*
*     N      INTEGER, the number of variables.
*     MAXPTS INTEGER, maximum number of function values allowed. This 
*            parameter can be used to limit the time taken. A 
*            sensible strategy is to start with MAXPTS = 1000*N, and then
*            increase MAXPTS if ERROR is too large.
*    FUNCTN  Externally declared real user defined integrand. Its 
*            parameters must be (N, Z), where Z is a real array of
*            length N.
*     ABSEPS REAL absolute error tolerance.
*     RELEPS REAL relative error tolerance.
*
*******Output  Parameters
*
*     ERROR  REAL estimated absolute error, with 99% confidence level.
*     VALUE  REAL estimated value for the integral
*     INFORM INTEGER, termination status parameter:
*            if INFORM = 0, normal completion with ERROR < EPS;
*            if INFORM = 1, completion with ERROR > EPS and MAXPTS 
*                           function vaules used; increase MAXPTS to 
*                           decrease ERROR;
*            if INFORM = 2, N > 20 or N < 1.
*
      INTEGER, INTENT(IN)  :: N,  MAXPTS
      INTEGER, INTENT(OUT) :: INFORM
      !INTEGER ::  NL, LENWRK, 
	INTEGER :: RULCLS, TOTCLS, NEWCLS, MAXCLS
      DOUBLE PRECISION, INTENT(IN)  :: ABSEPS, RELEPS
      DOUBLE PRECISION, INTENT(OUT) :: ERROR, VALUE
      DOUBLE PRECISION :: OLDVAL
      !PARAMETER ( NL = 20 )
      !PARAMETER ( LENWRK = 20*NL**2 )
      !DOUBLE PRECISION, DIMENSION(LENWRK) :: WORK
      INTERFACE
         DOUBLE PRECISION FUNCTION FUNCTN(N,Z)
         DOUBLE PRECISION,DIMENSION(*), INTENT(IN) :: Z
         INTEGER, INTENT(IN) :: N
         END FUNCTION FUNCTN
      END INTERFACE
      IF ( N .GT. MAXDIM .OR. N .LT. 1 ) THEN
         INFORM = 2
         VALUE = 0.d0
         ERROR = 1.d0
         RETURN
      ENDIF 
      INFORM = 1
*
*     Call the subregion adaptive integration subroutine
*     
      RULCLS = 1
      CALL ADAPT( N, RULCLS, 0, FUNCTN, ABSEPS, RELEPS, 
     &      ERROR, VALUE, INFORM )
      MAXCLS = MIN( 10*RULCLS, MAXPTS )
      TOTCLS = 0
      CALL ADAPT(N, TOTCLS, MAXCLS, FUNCTN, ABSEPS, RELEPS, 
     &     ERROR, VALUE, INFORM)
      IF ( ERROR .GT. MAX( ABSEPS, RELEPS*ABS(VALUE) ) ) THEN
 10      OLDVAL = VALUE
         MAXCLS = MAX( 2*RULCLS,MIN(INT(3*MAXCLS/2),MAXPTS-TOTCLS))
         NEWCLS = -1
         CALL ADAPT(N, NEWCLS, MAXCLS, FUNCTN, ABSEPS, RELEPS, 
     &        ERROR, VALUE, INFORM)
         TOTCLS = TOTCLS + NEWCLS
         ERROR = ABS(VALUE-OLDVAL) + 
     &        SQRT(DBLE(RULCLS)*ERROR**2/DBLE(TOTCLS))
         IF ( ERROR .GT. MAX( ABSEPS, RELEPS*ABS(VALUE) ) ) THEN
            IF ( MAXPTS - TOTCLS .GT. 2*RULCLS ) GO TO 10
         ELSE 
            INFORM = 0
         END IF
      ENDIF
      
      END SUBROUTINE SADAPT



!***********************************************************
!    MAIN INTEGRATION ROUTINE ADAPT
!***********************************************************  


      SUBROUTINE ADAPT(NDIM, MINCLS, MAXCLS, FUNCTN,
     &     ABSREQ, RELREQ, ABSEST, FINEST, INFORM)
      IMPLICIT NONE      
*
*   Adaptive Multidimensional Integration Subroutine
*
*   Author: Alan Genz
*           Department of Mathematics
*           Washington State University
*           Pullman, WA 99164-3113 USA
*
*  This subroutine computes an approximation to the integral
*
*      1 1     1
*     I I ... I       FUNCTN(NDIM,X)  dx(NDIM)...dx(2)dx(1)
*      0 0     0  
*
***************  Parameters for ADAPT  ********************************
*
****** Input Parameters
*
*  NDIM    Integer number of integration variables.
*  MINCLS  Integer minimum number of FUNCTN calls to be allowed; MINCLS
*          must not exceed MAXCLS. If MINCLS < 0, then ADAPT assumes
*          that a previous call of ADAPT has been made with the same
*          integrand and continues that calculation.
*  MAXCLS  Integer maximum number of FUNCTN calls to be used; MAXCLS
*          must be >= RULCLS, the number of function calls required for
*          one application of the basic integration rule.
*           IF ( NDIM .EQ. 1 ) THEN
*              RULCLS = 11
*           ELSE IF ( NDIM .LT. 15 ) THEN
*              RULCLS = 2**NDIM + 2*NDIM*(NDIM+3) + 1
*           ELSE
*              RULCLS = 1 + NDIM*(24-NDIM*(6-NDIM*4))/3
*           ENDIF
*  FUNCTN  Externally declared real user defined integrand. Its 
*          parameters must be (NDIM, Z), where Z is a real array of
*          length NDIM.
*  ABSREQ  Real required absolute accuracy.
*  RELREQ  Real required relative accuracy.
*
****** Output Parameters
*
*  MINCLS  Actual number of FUNCTN calls used by ADAPT.
*  ABSEST  Real estimated absolute accuracy.
*  FINEST  Real estimated value of integral.
*  INFORM  INFORM = 0 for normal exit, when ABSEST <= ABSREQ or
*                     ABSEST <= |FINEST|*RELREQ with MINCLS <= MAXCLS.
*          INFORM = 1 if MAXCLS was too small for ADAPT to obtain the
*                     result FINEST to within the requested accuracy.
*          INFORM = 2 if MINCLS > MAXCLS, LENWRK < 16*NDIM + 27 or 
*                     RULCLS > MAXCLS.
*
************************************************************************
*
*     Begin driver routine. This routine partitions the working storage 
*      array and then calls the main subroutine ADBASE.
*
* Revised pab 07.10.2000, 
*    1) Removed LENWRK and WORK from input. 
*    2) Defined LENWRK internally and Put a save statement before WORK instead 
*                          
*  LENWRK  Integer length of real array WORK (working storage); ADAPT
*          needs LENWRK >= 16*NDIM + 27. For maximum efficiency LENWRK
*          should be about 2*NDIM*MAXCLS/RULCLS if MAXCLS FUNCTN
*          calls are needed. If LENWRK is significantly less than this,
*          ADAPT may be less efficient.
*
*  WORK    Real array (length LENWRK) of working storage. This contains
*          information that is needed for additional calls of ADAPT
*          using the same integrand (input MINCLS < 0).
*
      INTEGER, INTENT(IN)    :: NDIM,  MAXCLS 
      INTEGER, INTENT(INOUT) :: MINCLS
      INTEGER, INTENT(OUT)   :: INFORM
      DOUBLE PRECISION, INTENT(IN)  :: ABSREQ, RELREQ 
      DOUBLE PRECISION, INTENT(OUT) :: ABSEST, FINEST
*     Local variables	 
      INTEGER, PARAMETER :: LENWRK=20*MAXDIM*MAXDIM 
      DOUBLE PRECISION, DIMENSION(LENWRK) :: WORK        ! length lenwrk
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: POINTS,WEGHTS,LUM 
      INTEGER :: SBRGNS, MXRGNS, RULCLS, LENRUL, 
     & INERRS, INVALS, INPTRS, INLWRS, INUPRS, INMSHS, INPNTS, INWGTS, 
     & INLOWR, INUPPR, INWDTH, INMESH, INWORK
      INTERFACE
         DOUBLE PRECISION FUNCTION FUNCTN(N,Z)
         DOUBLE PRECISION,DIMENSION(*), INTENT(IN) :: Z
         INTEGER, INTENT(IN) :: N
         END FUNCTION FUNCTN
      END INTERFACE
	SAVE WORK
!      print *,'adapt, ndim', ndim 
      IF ( NDIM .EQ. 1 ) THEN
         LENRUL = 5
         RULCLS = 9
      ELSE IF ( NDIM .LT. 12 ) THEN
         LENRUL = 6
         RULCLS = 2**NDIM + 2*NDIM*(NDIM+2) + 1
      ELSE
         LENRUL = 6
!         RULCLS = 1 + 2*NDIM*(1+2*NDIM)  ! old call pab 15.03.2003
         RULCLS = 1851 + 2*NDIM*(1+2*NDIM)
      ENDIF
      IF ( LENWRK .GE. LENRUL*(NDIM+4) + 10*NDIM + 3 .AND.
     &     RULCLS. LE. MAXCLS .AND. MINCLS .LE. MAXCLS ) THEN
         MXRGNS = ( LENWRK - LENRUL*(NDIM+4) - 7*NDIM )/( 3*NDIM + 3 )
         INERRS = 1
         INVALS = INERRS + MXRGNS
         INPTRS = INVALS + MXRGNS
         INLWRS = INPTRS + MXRGNS
         INUPRS = INLWRS + MXRGNS*NDIM
         INMSHS = INUPRS + MXRGNS*NDIM
         INWGTS = INMSHS + MXRGNS*NDIM
         INPNTS = INWGTS + LENRUL*4
         INLOWR = INPNTS + LENRUL*NDIM
         INUPPR = INLOWR + NDIM
         INWDTH = INUPPR + NDIM
         INMESH = INWDTH + NDIM
         INWORK = INMESH + NDIM
          
         ALLOCATE(POINTS(NDIM,LENRUL))
         ALLOCATE(WEGHTS(LENRUL,4))
         ALLOCATE(LUM(NDIM,MXRGNS*3))

         IF (	MINCLS .LT. 0 ) THEN
            SBRGNS = WORK(LENWRK)
            LUM    = reshape(WORK(INLWRS:INWGTS-1),(/ NDIM,MXRGNS*3/))
            WEGHTS = reshape(WORK(INWGTS:INPNTS-1),(/ LENRUL , 4 /))
            POINTS = reshape(WORK(INPNTS:INLOWR-1),(/ NDIM, LENRUL/))
	   !ELSE
	   !   WORK=0.D0;LUM=0.D0;WEGHTS=0.D0;POINTS=0.D0;	
	   ENDIF	
         CALL ADBASE(NDIM, MINCLS, MAXCLS, FUNCTN, ABSREQ, RELREQ, 
     &        ABSEST, FINEST, SBRGNS, MXRGNS, RULCLS, LENRUL, 
     &        WORK(INERRS:INVALS-1), WORK(INVALS:INPTRS-1), 
     &        WORK(INPTRS:INLWRS-1), LUM(:,1:MXRGNS), 
     &        LUM(:,MXRGNS+1:2*MXRGNS),LUM(:,2*MXRGNS+1:3*MXRGNS), 
     &        WEGHTS,POINTS,WORK(INLOWR:INUPPR-1),WORK(INUPPR:INWDTH-1), 
     &        WORK(INWDTH:INMESH-1), WORK(INMESH:INWORK-1), 
     &        WORK(INWORK:INWORK+2*NDIM-1), INFORM)
         WORK(LENWRK) = SBRGNS
         WORK(INLWRS:INWGTS-1) = reshape(LUM   ,(/ NDIM*MXRGNS*3/))  ! LOWERS UPPERS MESHES
         WORK(INWGTS:INPNTS-1) = reshape(WEGHTS,(/ LENRUL*4 /))
         WORK(INPNTS:INLOWR-1) = reshape(POINTS,(/ NDIM*LENRUL/))
         DEALLOCATE(POINTS)
         DEALLOCATE(WEGHTS)
         DEALLOCATE(LUM)
      ELSE
         INFORM = 2
         MINCLS = RULCLS
      ENDIF
      RETURN
      END SUBROUTINE ADAPT
      SUBROUTINE BSINIT(NDIM, W, LENRUL, G)
      IMPLICIT NONE
*
*     For initializing basic rule weights and symmetric sum parameters.
*
      INTEGER, INTENT(IN) :: NDIM, LENRUL
      DOUBLE PRECISION , DIMENSION(:,:), INTENT(OUT) :: W, G
*      DOUBLE PRECISION W(LENRUL,4), G(NDIM,LENRUL) 
*    Local variables
	INTEGER :: I, J
	INTEGER, PARAMETER :: NUMNUL=4, SDIM=12
      INTEGER, DIMENSION(6) ::  RULPTS
      DOUBLE PRECISION LAM1, LAM2, LAM3, LAM4, LAMP, RULCON
*
*     The following code determines rule parameters and weights for a
*      degree 7 rule (W(1,1),...,W(5,1)), two degree 5 comparison rules
*      (W(1,2),...,W(5,2) and W(1,3),...,W(5,3)) and a degree 3 
*      comparison rule (W(1,4),...W(5,4)).
*
*       If NDIM = 1, then LENRUL = 5 and total points = 9.
*       If NDIM < SDIM, then LENRUL = 6 and
*                      total points = 1+2*NDIM*(NDIM+2)+2**NDIM.
*       If NDIM > = SDIM, then LENRUL = 6 and
*                      total points = 1+2*NDIM*(1+2*NDIM).
*
!      print *,'BSINIT, ndim', ndim 
!      DO I = 1,LENRUL
!         DO J = 1,NDIM
!            G(J,I) = 0.d0
!         END DO
!         DO J = 1,NUMNUL
!            W(I,J) = 0.d0
!         END DO
!      END DO
      G = 0.D0
      W = 0.D0
      I = 2*NDIM
      RULPTS(5) = I*(NDIM-1)
      RULPTS(4) = I 
      RULPTS(3) = I
      RULPTS(2) = I
      RULPTS(1) = 1
      LAMP = 0.85d0
      LAM3 = 0.4707d0
      LAM2 = 4d0/(15.d0 - 5.d0/LAM3)
	LAM4 = 1.D0/(27.D0*LAM3*LAM3*LAM3)
      W(5,1) = ( 3.d0 - 5.d0*LAM3 )/( 180.d0*(LAM2-LAM3)*LAM2*LAM2)
      IF ( NDIM .LT. SDIM ) THEN 
         RULPTS(LENRUL) = 2**NDIM
         LAM1 = 8.d0*LAM3*(31.d0*LAM3-15.d0)/
     &        ( (3.d0*LAM3-1.d0)*(5.d0*LAM3-3.d0)*35.d0 )
         W(LENRUL,1) = LAM4/DBLE(RULPTS(LENRUL))  
	   
      ELSE
         LAM1 = ( LAM3*(15.d0 - 21.d0*LAM2) + 
     &        35.d0*DBLE(NDIM-1)*(LAM2-LAM3)/9.d0 )
     &       /  ( LAM3*(21.d0 - 35.d0*LAM2) + 
     &        35.d0*DBLE(NDIM-1)*(LAM2/LAM3-1.d0)/9.d0 )
         W(6,1) = LAM4*0.25D0  
   	   RULPTS(6) = 2*NDIM*(NDIM-1)

      ENDIF
      W(3,1) = ( 15.d0 - 21.d0*(LAM3+LAM1) + 35.d0*LAM3*LAM1 )
     &   /(210.d0*LAM2*(LAM2-LAM3)*(LAM2-LAM1))-DBLE(2*(NDIM-1))*W(5,1)
      W(2,1) = ( 15.d0 - 21.d0*(LAM3+LAM2) + 35.d0*LAM3*LAM2 )
     &     /( 210.d0*LAM1*(LAM1-LAM3)*(LAM1-LAM2) )
      LAM3 = SQRT(LAM3)
      IF ( NDIM .LT. SDIM ) THEN             
          G(1:NDIM,LENRUL) = LAM3
      ELSE
          G(1,6) = LAM3
   	    G(2,6) = LAM3
      ENDIF
      IF ( NDIM .GT. 1 ) THEN
         W(5,2) = 1.d0/(6.d0*LAM2)**2 
         W(5,3) = W(5,2) 
      ENDIF
      W(3,2) = ( 3.d0 - 5.d0*LAM1 )/( 30.d0*LAM2*(LAM2-LAM1) ) 
     &     - DBLE(2*(NDIM-1))*W(5,2) 
      W(2,2) = ( 3.d0 - 5.d0*LAM2 )/( 30.d0*LAM1*(LAM1-LAM2) )
      W(4,3) = ( 3.d0 - 5.d0*LAM2 )/( 30.d0*LAMP*(LAMP-LAM2) )
      W(3,3) = ( 3.d0 - 5.d0*LAMP )/( 30.d0*LAM2*(LAM2-LAMP) ) 
     &     - DBLE(2*(NDIM-1))*W(5,3)
      W(2,4) = 1.d0/(6.d0*LAM1)
      LAMP = SQRT(LAMP)
      LAM2 = SQRT(LAM2)
      LAM1 = SQRT(LAM1)
      G(1,2) = LAM1
      G(1,3) = LAM2
      G(1,4) = LAMP
      IF ( NDIM .GT. 1 ) THEN
         G(1,5) = LAM2
         G(2,5) = LAM2
      ENDIF
      DO J = 1, NUMNUL
         W(1,J) = 1.d0
         DO I = 2,LENRUL
            W(1,J) = W(1,J) - DBLE(RULPTS(I))*W(I,J)
         END DO
      END DO
      RULCON = 0.5d0
      CALL RULNRM( LENRUL, NUMNUL, RULPTS, W, RULCON )
      END SUBROUTINE BSINIT
!
!
      SUBROUTINE RULNRM( LENRUL, NUMNUL, RULPTS, W, RULCON )
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: LENRUL, NUMNUL
      INTEGER, DIMENSION(:), INTENT(IN) :: RULPTS
      DOUBLE PRECISION, DIMENSION(:,:), INTENT(INOUT) :: W       !(LENRUL, *),
      DOUBLE PRECISION, INTENT(IN) ::  RULCON
*     Local variables
	INTEGER          :: I, J, K
      DOUBLE PRECISION :: ALPHA, NORMCF, NORMNL
	
      
*
*     Compute orthonormalized null rules.
*
!      print *,'RULNRM, lenrul, numnul', lenrul,NUMNUL 
      NORMCF = 0.d0
      DO I = 1,LENRUL
         NORMCF = NORMCF + DBLE(RULPTS(I))*W(I,1)*W(I,1)
      END DO
      DO K = 2,NUMNUL
         DO I = 1,LENRUL
            W(I,K) = W(I,K) - W(I,1)
         END DO
         DO J = 2,K-1
            ALPHA = 0.d0
            DO I = 1,LENRUL
               ALPHA = ALPHA + DBLE(RULPTS(I))*W(I,J)*W(I,K)
            END DO
            ALPHA = -ALPHA/NORMCF
            DO I = 1,LENRUL
               W(I,K) = W(I,K) + ALPHA*W(I,J)
            END DO
         END DO
         NORMNL = 0.d0
         DO I = 1,LENRUL
            NORMNL = NORMNL + DBLE(RULPTS(I))*W(I,K)*W(I,K)
         END DO
         ALPHA = SQRT(NORMCF/NORMNL)
         DO I = 1,LENRUL
            W(I,K) = ALPHA*W(I,K)
         END DO
      END DO
      DO J = 2, NUMNUL
         DO I = 1,LENRUL
            W(I,J) = W(I,J)*RULCON
         END DO
      END DO
      RETURN
      END SUBROUTINE RULNRM
!
!
      SUBROUTINE ADBASE(NDIM, MINCLS, MAXCLS, FUNCTN, ABSREQ, RELREQ,
     &     ABSEST, FINEST, SBRGNS, MXRGNS, RULCLS, LENRUL,
     &     ERRORS, VALUES, PONTRS, LOWERS, 
     &     UPPERS, MESHES, WEGHTS, POINTS, 
     &     LOWER, UPPER, WIDTH, MESH, WORK, INFORM)
      IMPLICIT NONE
*
*        Main adaptive integration subroutine
*
      INTEGER,INTENT(IN) :: NDIM,  MAXCLS, MXRGNS,LENRUL, RULCLS
	INTEGER, INTENT(INOUT) :: MINCLS, SBRGNS  
      INTEGER, INTENT(OUT) :: INFORM 
      DOUBLE PRECISION, INTENT(IN) :: ABSREQ, RELREQ
	DOUBLE PRECISION, INTENT(OUT) :: ABSEST, FINEST  
      DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) ::  ERRORS, VALUES, 
     &	PONTRS, LOWER, UPPER, WIDTH, MESH, WORK
      DOUBLE PRECISION, DIMENSION(:,:), INTENT(INOUT) :: WEGHTS, POINTS             
	! shape (LENRUL,4) and (NDIM,LENRUL)
      DOUBLE PRECISION, DIMENSION(:,:), INTENT(INOUT) :: LOWERS, UPPERS, 
     &	MESHES     !SHAPE  (NDIM,MXRGNS),   
      INTEGER :: I, J,NWRGNS,  DIVAXN, TOP, RGNCLS, FUNCLS, DIFCLS      
      INTERFACE
         DOUBLE PRECISION FUNCTION FUNCTN(N,Z)
         DOUBLE PRECISION,DIMENSION(*), INTENT(IN) :: Z
         INTEGER, INTENT(IN) :: N
         END FUNCTION FUNCTN
      END INTERFACE
*
*     Initialization of subroutine
*
!      print *,'ADBASE, ndim', ndim, shape(POINTS) 
	SAVE DIVAXN     ! added pab 07.11.2000 (divaxn may have negative values otherwise)
      INFORM = 2
      FUNCLS = 0
      CALL BSINIT(NDIM, WEGHTS, LENRUL, POINTS)
      IF ( MINCLS .GE. 0) THEN
*
*       When MINCLS >= 0 determine initial subdivision of the
*       integration region and apply basic rule to each subregion.
*
         SBRGNS = 0
         DO I = 1,NDIM
            LOWER(I) = 0.d0
            MESH(I) = 1.d0
            WIDTH(I) = 1.d0/(2.d0*MESH(I))
            UPPER(I) = 1.d0
         END DO
         DIVAXN = 0
         RGNCLS = RULCLS
         NWRGNS = 1
 10      CONTINUE
	   !IF (abs(DIVAXN).GT.NDIM)   PRINT *,'adbase DIVAXN1',DIVAXN
         CALL DIFFER(NDIM, LOWER, UPPER, WIDTH, WORK(1:NDIM),  
     &        WORK(NDIM+1:2*NDIM), FUNCTN, DIVAXN, DIFCLS)
         FUNCLS = FUNCLS + DIFCLS
         IF (DBLE(RGNCLS)*(MESH(DIVAXN)+1.d0)/MESH(DIVAXN)
     &        .LE. DBLE(MINCLS-FUNCLS) ) THEN
            RGNCLS = NINT(DBLE(RGNCLS)*(MESH(DIVAXN)+1.d0)/MESH(DIVAXN))
            NWRGNS = NINT(DBLE(NWRGNS)*(MESH(DIVAXN)+1.d0)/MESH(DIVAXN))
            MESH(DIVAXN) = MESH(DIVAXN) + 1.d0
            WIDTH(DIVAXN) = 1.d0/( 2.d0*MESH(DIVAXN) )
            GO TO 10
         ENDIF
         IF ( NWRGNS .LE. MXRGNS ) THEN
            DO I = 1,NDIM
               UPPER(I) = LOWER(I) + 2.d0*WIDTH(I)
               MESH(I) = 1.d0
            END DO
         ENDIF
*     
*     Apply basic rule to subregions and store results in heap.
*     
 20      SBRGNS = SBRGNS + 1
         CALL BASRUL(NDIM, LOWER, UPPER, WIDTH, FUNCTN, 
     &        WEGHTS, LENRUL, POINTS, WORK(1:NDIM), WORK(NDIM+1:2*NDIM), 
     &        ERRORS(SBRGNS),VALUES(SBRGNS))
         CALL TRESTR(SBRGNS, SBRGNS, PONTRS, ERRORS)
         DO I = 1,NDIM
            LOWERS(I,SBRGNS) = LOWER(I)
            UPPERS(I,SBRGNS) = UPPER(I)
            MESHES(I,SBRGNS) = MESH(I)
         END DO
         DO I = 1,NDIM
            LOWER(I) = UPPER(I)
            UPPER(I) = LOWER(I) + 2.d0*WIDTH(I)
            IF (LOWER(I)+WIDTH(I) .LT. 1.D0)  GO TO 20
            LOWER(I) = 0.d0
            UPPER(I) = LOWER(I) + 2.d0*WIDTH(I)
         END DO
         FUNCLS = FUNCLS + SBRGNS*RULCLS
      ENDIF
*     
*     Check for termination
*
 30   FINEST = 0.d0
      ABSEST = 0.d0
      DO I = 1, SBRGNS
         FINEST = FINEST + VALUES(I)
         ABSEST = ABSEST + ERRORS(I)
      END DO
      IF ( ABSEST .GT. MAX( ABSREQ, RELREQ*ABS(FINEST) )
     &     .OR. FUNCLS .LT. MINCLS ) THEN  
*     
*     Prepare to apply basic rule in (parts of) subregion with
*     largest error.
*     
         TOP = PONTRS(1)
         RGNCLS = RULCLS
         DO I = 1,NDIM
            LOWER(I) = LOWERS(I,TOP)
            UPPER(I) = UPPERS(I,TOP)
            MESH(I) = MESHES(I,TOP)
            WIDTH(I) = (UPPER(I)-LOWER(I))/(2.D0*MESH(I))
            RGNCLS = NINT(DBLE(RGNCLS)*MESH(I))
         END DO
	   !IF (abs(DIVAXN).GT.NDIM)   PRINT *,'adbase DIVAXN2',DIVAXN
         CALL DIFFER(NDIM, LOWER, UPPER, WIDTH, WORK(1:NDIM),   
     &       WORK(NDIM+1:2*NDIM), FUNCTN, DIVAXN, DIFCLS)
         FUNCLS = FUNCLS + DIFCLS
         RGNCLS = NINT(DBLE(RGNCLS)*(MESH(DIVAXN)+1.D0))/MESH(DIVAXN)
         IF ( FUNCLS + RGNCLS .LE. MAXCLS ) THEN
            IF ( SBRGNS + 1 .LE. MXRGNS ) THEN
*     
*     Prepare to subdivide into two pieces.
*    
               NWRGNS = 1
               WIDTH(DIVAXN) = 0.5d0*WIDTH(DIVAXN)
            ELSE
               NWRGNS = 0
               WIDTH(DIVAXN) = WIDTH(DIVAXN)
     &                        *MESH(DIVAXN)/( MESH(DIVAXN) + 1.d0 )
               MESHES(DIVAXN,TOP) = MESH(DIVAXN) + 1.d0 
            ENDIF
            IF ( NWRGNS .GT. 0 ) THEN
*     
*     Only allow local subdivision when space is available.
*
               DO J = SBRGNS+1,SBRGNS+NWRGNS
                  DO I = 1,NDIM
                     LOWERS(I,J) = LOWER(I)
                     UPPERS(I,J) = UPPER(I)
                     MESHES(I,J) = MESH(I)
                  END DO
               END DO
               UPPERS(DIVAXN,TOP) = LOWER(DIVAXN) + 2.d0*WIDTH(DIVAXN)
               LOWERS(DIVAXN,SBRGNS+1) = UPPERS(DIVAXN,TOP)
            ENDIF
            FUNCLS = FUNCLS + RGNCLS
            CALL BASRUL(NDIM, LOWERS(:,TOP), UPPERS(:,TOP), WIDTH, 
     &           FUNCTN, WEGHTS, LENRUL, POINTS, WORK(1:NDIM),  
     &           WORK(NDIM+1:2*NDIM),ERRORS(TOP), VALUES(TOP))
            CALL TRESTR(TOP, SBRGNS, PONTRS, ERRORS)
            DO I = SBRGNS+1, SBRGNS+NWRGNS
*     
*     Apply basic rule and store results in heap.
*     
               CALL BASRUL(NDIM, LOWERS(:,I), UPPERS(:,I), WIDTH,
     &              FUNCTN, WEGHTS, LENRUL, POINTS, WORK(1:NDIM),   
     &              WORK(NDIM+1:2*NDIM),ERRORS(I), VALUES(I))
               CALL TRESTR(I, I, PONTRS, ERRORS)
            END DO
            SBRGNS = SBRGNS + NWRGNS
            GO TO 30
         ELSE
            INFORM = 1
         ENDIF
      ELSE
         INFORM = 0
      ENDIF
      MINCLS = FUNCLS
      RETURN 
      END SUBROUTINE ADBASE
      SUBROUTINE BASRUL( NDIM, A, B, WIDTH, FUNCTN, W, LENRUL, G,
     &     CENTER, Z, RGNERT, BASEST )
      IMPLICIT NONE
*
*     For application of basic integration rule
*
      INTEGER, INTENT(IN) :: LENRUL, NDIM
      DOUBLE PRECISION, DIMENSION(:  ), INTENT(IN) :: A, B, WIDTH         !(NDIM)
      DOUBLE PRECISION, DIMENSION(:,:), INTENT(IN) :: W                   !(LENRUL,4), 
      DOUBLE PRECISION, DIMENSION(:,:), INTENT(INOUT) :: G                !(NDIM,LENRUL), 
      DOUBLE PRECISION, DIMENSION(:  ), INTENT(INOUT) :: CENTER, Z        !(NDIM)
      DOUBLE PRECISION, INTENT(OUT) :: RGNERT, BASEST
	INTEGER :: I
      DOUBLE PRECISION :: FSYMSM, RGNCMP, RGNVAL,
     &     RGNVOL, RGNCPT, RGNERR
      INTERFACE
         DOUBLE PRECISION FUNCTION FUNCTN(N,Z)
         DOUBLE PRECISION,DIMENSION(*), INTENT(IN) :: Z
         INTEGER, INTENT(IN) :: N
         END FUNCTION FUNCTN
      END INTERFACE
*
*     Compute Volume and Center of Subregion
*
!      print *,'BASRULE, ndim', ndim 
      RGNVOL = 1.d0
      DO I = 1,NDIM
         RGNVOL = 2.d0*RGNVOL*WIDTH(I)
         CENTER(I) = A(I) + WIDTH(I)
      END DO
      BASEST = 0.d0
      RGNERT = 0.d0
*
*     Compute basic rule and error
*
 10   RGNVAL = 0.d0
      RGNERR = 0.d0
      RGNCMP = 0.d0
      RGNCPT = 0.d0
      DO I = 1,LENRUL
         FSYMSM = FULSUM(NDIM, CENTER, WIDTH, Z, G(:,I), FUNCTN)
*     Basic Rule
         RGNVAL = RGNVAL + W(I,1)*FSYMSM
*     First comparison rule
         RGNERR = RGNERR + W(I,2)*FSYMSM
*     Second comparison rule
         RGNCMP = RGNCMP + W(I,3)*FSYMSM
*     Third Comparison rule
         RGNCPT = RGNCPT + W(I,4)*FSYMSM
      END DO
*
*     Error estimation
*
      RGNERR = SQRT(RGNCMP*RGNCMP + RGNERR*RGNERR)
      RGNCMP = SQRT(RGNCPT*RGNCPT + RGNCMP*RGNCMP)
      IF ( 4.d0*RGNERR .LT. RGNCMP ) RGNERR = 0.5d0*RGNERR
      IF ( 2.d0*RGNERR .GT. RGNCMP ) RGNERR = MAX( RGNERR, RGNCMP )
      RGNERT = RGNERT +  RGNVOL*RGNERR
      BASEST = BASEST +  RGNVOL*RGNVAL
*
*     When subregion has more than one piece, determine next piece and
*      loop back to apply basic rule.
*
      DO I = 1,NDIM
         CENTER(I) = CENTER(I) + 2.d0*WIDTH(I)
         IF ( CENTER(I) .LT. B(I) ) GO TO 10
         CENTER(I) = A(I) + WIDTH(I)
      END DO
      RETURN
      END SUBROUTINE BASRUL
      DOUBLE PRECISION FUNCTION FULSUM(S, CENTER, HWIDTH, X, G, F)
      IMPLICIT NONE
*
****  To compute fully symmetric basic rule sum
* 
      INTEGER, INTENT(IN) :: S
      DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: CENTER, HWIDTH
	DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: X, G                     ! shape S
      INTEGER          :: IXCHNG, LXCHNG, I, L
      DOUBLE PRECISION :: INTSUM, GL, GI
      INTERFACE
         DOUBLE PRECISION FUNCTION F(N,Z)
         DOUBLE PRECISION,DIMENSION(*), INTENT(IN) :: Z
         INTEGER, INTENT(IN) :: N
         END FUNCTION F
      END INTERFACE
!      print *,'FULSUM, S', S, shape(X) 
      FULSUM = 0.d0
*
*     Compute centrally symmetric sum for permutation of G
*
 10   INTSUM = 0.d0
      !DO I = 1,S
      !   X(I) = CENTER(I) + G(I)*HWIDTH(I)
      !END DO
	X = CENTER + G*HWIDTH
 20   INTSUM = INTSUM + F(S,X)
      DO I = 1,S
         G(I) = -G(I)
         X(I) = CENTER(I) + G(I)*HWIDTH(I)
         IF ( G(I) .LT. 0.d0 ) GO TO 20
      END DO
      FULSUM = FULSUM + INTSUM
*     
*     Find next distinct permuation of G and loop back for next sum
*     
      DO I = 2,S
         IF ( G(I-1) .GT. G(I) ) THEN
            GI = G(I)
            IXCHNG = I - 1
            DO L = 1,(I-1)/2
               GL = G(L)
               G(L) = G(I-L)
               G(I-L) = GL
               IF (  GL  .LE. GI ) IXCHNG = IXCHNG - 1
               IF ( G(L) .GT. GI ) LXCHNG = L
            END DO
            IF ( G(IXCHNG) .LE. GI ) IXCHNG = LXCHNG
            G(I) = G(IXCHNG)
            G(IXCHNG) = GI
            GO TO 10
         ENDIF
      END DO
*     
*     End loop for permutations of G and associated sums
*     
*     Restore original order to G's
*     
      DO I = 1,S/2
         GI = G(I)
         G(I) = G(S+1-I)
         G(S+1-I) = GI 
      END DO
      RETURN
      END FUNCTION FULSUM
      SUBROUTINE DIFFER(NDIM, A, B, WIDTH, Z, DIF, FUNCTN, 
     &     DIVAXN, DIFCLS)
      IMPLICIT NONE
*
*     Compute fourth differences and subdivision axes
*
      INTEGER, INTENT(IN)    :: NDIM
	INTEGER, INTENT(INOUT) :: DIVAXN
	INTEGER, INTENT(OUT)   :: DIFCLS
      DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: A, B, WIDTH   ! (NDIM)
	DOUBLE PRECISION, DIMENSION(:),INTENT(OUT) :: Z, DIF        ! (NDIM)
      DOUBLE PRECISION :: FRTHDF, FUNCEN, WIDTHI
	INTEGER :: I
      INTERFACE
         DOUBLE PRECISION FUNCTION FUNCTN(N,Z)
         DOUBLE PRECISION,DIMENSION(*), INTENT(IN) :: Z
         INTEGER, INTENT(IN) :: N
         END FUNCTION FUNCTN
      END INTERFACE
!      print *,'DIFFER, ndim', ndim, shape(Z) 
      DIFCLS = 0
	!IF (abs(DIVAXN).GT.NDIM)   PRINT *,'DIFFER DIVAXN1',DIVAXN
	
      DIVAXN = MOD(DIVAXN, NDIM ) + 1
	!print *,'DIFFER, divaxn2', divaxn
      IF ( NDIM .GT. 1 ) THEN
         !DO I = 1,NDIM 
         !  DIF(I) = 0.d0
         !   Z(I) = A(I) + WIDTH(I)
         !END DO
	   DIF = 0.D0
	   Z(1:NDIM) = A(1:NDIM) + WIDTH(1:NDIM)
!         print *,'Z', Z
 10      FUNCEN = FUNCTN(NDIM, Z)
         DO I = 1,NDIM
            WIDTHI = 0.2d0*WIDTH(I)
            FRTHDF = 6.d0*FUNCEN
            Z(I) = Z(I) - 4.d0*WIDTHI
            FRTHDF = FRTHDF + FUNCTN(NDIM,Z)
            Z(I) = Z(I) + 2.d0*WIDTHI
            FRTHDF = FRTHDF - 4.d0*FUNCTN(NDIM,Z)
            Z(I) = Z(I) + 4.d0*WIDTHI
            FRTHDF = FRTHDF - 4.d0*FUNCTN(NDIM,Z)
            Z(I) = Z(I) + 2.d0*WIDTHI
            FRTHDF = FRTHDF + FUNCTN(NDIM,Z)
*     Do not include differences below roundoff
!            IF ( FUNCEN + FRTHDF/8.d0 .NE. FUNCEN ) 
             IF ( FUNCEN + FRTHDF*0.125D0 .NE. FUNCEN ) 
     &           DIF(I) = DIF(I) + ABS(FRTHDF)*WIDTH(I)
            Z(I) = Z(I) - 4.d0*WIDTHI
         END DO
         DIFCLS = DIFCLS + 4*NDIM + 1
         DO I = 1,NDIM
            Z(I) = Z(I) + 2.D0*WIDTH(I)
            IF ( Z(I) .LT. B(I) ) GO TO 10
            Z(I) = A(I) + WIDTH(I)
         END DO
	   !IF (abs(DIVAXN).GT.NDIM)   PRINT *,'DIFFER DIVAXN',DIVAXN,shape(dif),ndim
         DO I = 1,NDIM	        
            IF ( DIF(DIVAXN) .LT. DIF(I) ) DIVAXN = I
         END DO
      ENDIF
      RETURN
      END SUBROUTINE DIFFER
      SUBROUTINE TRESTR(POINTR, SBRGNS, PONTRS, RGNERS)
      IMPLICIT NONE
****BEGIN PROLOGUE TRESTR
****PURPOSE TRESTR maintains a heap for subregions.
****DESCRIPTION TRESTR maintains a heap for subregions.
*            The subregions are ordered according to the size of the
*            greatest error estimates of each subregion (RGNERS).
*
*   PARAMETERS
*
*     POINTR Integer.
*            The index for the subregion to be inserted in the heap.
*     SBRGNS Integer.
*            Number of subregions in the heap.
*     PONTRS Real array of dimension SBRGNS.
*            Used to store the indices for the greatest estimated errors
*            for each subregion.
*     RGNERS Real array of dimension SBRGNS.
*            Used to store the greatest estimated errors for each 
*            subregion.
*
****ROUTINES CALLED NONE
****END PROLOGUE TRESTR
*
*   Global variables.
*
      INTEGER, INTENT(IN) ::POINTR, SBRGNS
      DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: PONTRS
	DOUBLE PRECISION, DIMENSION(:), INTENT(IN)    :: RGNERS                 !(*)
*
*   Local variables.
*
*   RGNERR Intermediate storage for the greatest error of a subregion.
*   SUBRGN Position of child/parent subregion in the heap.
*   SUBTMP Position of parent/child subregion in the heap.
*
      INTEGER SUBRGN, SUBTMP
      DOUBLE PRECISION RGNERR
*
****FIRST PROCESSING STATEMENT TRESTR
*     
!      print *,'TRESTR' 
      RGNERR = RGNERS(POINTR)
      IF ( POINTR.EQ.NINT(PONTRS(1))) THEN
*
*        Move the new subregion inserted at the top of the heap 
*        to its correct position in the heap.
*
         SUBRGN = 1
 10      SUBTMP = 2*SUBRGN
         IF ( SUBTMP .LE. SBRGNS ) THEN
            IF ( SUBTMP .NE. SBRGNS ) THEN
*     
*              Find maximum of left and right child.
*
               IF ( RGNERS(NINT(PONTRS(SUBTMP))) .LT. 
     &              RGNERS(NINT(PONTRS(SUBTMP+1))) ) SUBTMP = SUBTMP + 1
            ENDIF
*
*           Compare maximum child with parent.
*           If parent is maximum, then done.
*
            IF ( RGNERR .LT. RGNERS(NINT(PONTRS(SUBTMP))) ) THEN
*     
*              Move the pointer at position subtmp up the heap.
*     
               PONTRS(SUBRGN) = PONTRS(SUBTMP)
               SUBRGN = SUBTMP
               GO TO 10
            ENDIF	
         ENDIF
      ELSE
*
*        Insert new subregion in the heap.
*
         SUBRGN = SBRGNS
 20      SUBTMP = SUBRGN/2
         IF ( SUBTMP .GE. 1 ) THEN
*
*           Compare child with parent. If parent is maximum, then done.
*     
            IF ( RGNERR .GT. RGNERS(NINT(PONTRS(SUBTMP))) ) THEN
*     
*              Move the pointer at position subtmp down the heap.
*
               PONTRS(SUBRGN) = PONTRS(SUBTMP)
               SUBRGN = SUBTMP
               GO TO 20
            ENDIF
         ENDIF
      ENDIF
      PONTRS(SUBRGN) = DBLE(POINTR)
*
****END TRESTR
*
      RETURN
      END SUBROUTINE TRESTR
      END MODULE ADAPTMOD       



*  RCRUDEMOD is a module  containing two:
*
*  Automatic Multidimensional Integration Subroutines
*               
*         AUTHOR: Alan Genz
*                 Department of Mathematics
*                 Washington State University
*                 Pulman, WA 99164-3113
*                 Email: AlanGenz@wsu.edu
*
*         Last Change: 5/15/98
* revised pab 10.03.2000
*   - updated to f90 (i.e. changed to assumed shape arrays + changing integers to DBLE)
*   - put it into a module
*   - added ranlhmc
*
*  RCRUDEMOD computes an approximation to the integral
*
*      1  1     1
*     I  I ... I       F(X)  dx(NDIM)...dx(2)dx(1)
*     0  0     0
! References:
! Alan Genz (1992)
! 'Numerical Computation of Multivariate Normal Probabilites'
! J. computational Graphical Statistics, Vol.1, pp 141--149              (RANMC)
!
! William H. Press, Saul Teukolsky, 
! William T. Wetterling and Brian P. Flannery (1997)
! "Numerical recipes in Fortran 77", Vol. 1, pp 55-63            (SVDCMP,PYTHAG)
!
! Donald E. Knuth (1973) "The art of computer programming,",
! Vol. 3, pp 84-  (sorting and searching)                               (SORTRE)


! You may  initialize the random generator before you 
!  call RANLHMC or RANMC by the following lines:
!
!      call random_seed(SIZE=seed_size) 
!      allocate(seed(seed_size)) 
!      call random_seed(GET=seed(1:seed_size))  ! get current seed
!      seed(1)=seed1                            ! change seed
!      call random_seed(PUT=seed(1:seed_size)) 
!      deallocate(seed)



      MODULE RCRUDEMOD
      IMPLICIT NONE
      PRIVATE
      PUBLIC :: RANMC, RANLHMC
      INTEGER ::  NDIMMAX       

      INTERFACE RANMC
      MODULE PROCEDURE RANMC
      END INTERFACE
      
      INTERFACE RCRUDE
      MODULE PROCEDURE RCRUDE
      END INTERFACE

      INTERFACE SVDCMP
      MODULE PROCEDURE SVDCMP
      END INTERFACE

      INTERFACE PYTHAG
      MODULE PROCEDURE PYTHAG
      END INTERFACE

      INTERFACE SPEARCORR
      MODULE PROCEDURE SPEARCORR
      END INTERFACE
       
      INTERFACE SORTRE
      MODULE PROCEDURE SORTRE
      END INTERFACE
      
      INTERFACE BINSORT
      MODULE PROCEDURE BINSORT
      END INTERFACE

      INTERFACE SWAPRE
      MODULE PROCEDURE SWAPRE
      END INTERFACE

      INTERFACE SWAPINT
      MODULE PROCEDURE SWAPINT
      END INTERFACE

      PARAMETER (NDIMMAX=100)
                                    !--------------------------------
      CONTAINS   
      SUBROUTINE RANMC( N, MAXPTS, FUNCTN, ABSEPS, 
     &     RELEPS, ERROR, VALUE, INFORM )
      IMPLICIT NONE
*
*     A subroutine for computing multivariate integrals.
*     This subroutine uses the Monte-Carlo algorithm given in the paper
*     "Numerical Computation of Multivariate Normal Probabilities", in
*     J. of Computational and Graphical Stat., 1(1992), pp. 141-149, by
*          Alan Genz
*          Department of Mathematics
*          Washington State University
*          Pullman, WA 99164-3113
*          Email : alangenz@wsu.edu
*
*  This subroutine computes an approximation to the integral
*
*      1 1     1
*     I I ... I       FUNCTN(NDIM,X)  dx(NDIM)...dx(2)dx(1)
*      0 0     0  
*
***************  Parameters for RANMC  ********************************
*
****** Input Parameters
*
*     N      INTEGER, the number of variables.
*     MAXPTS INTEGER, maximum number of function values allowed. This 
*            parameter can be used to limit the time taken. A 
*            sensible strategy is to start with MAXPTS = 1000*N, and then
*            increase MAXPTS if ERROR is too large.
*     ABSEPS REAL absolute error tolerance.
*     RELEPS REAL relative error tolerance.
*
****** Output Parameters
*
*     ERROR  REAL estimated absolute error, with 99% confidence level.
*     VALUE  REAL estimated value for the integral
*     INFORM INTEGER, termination status parameter:
*            if INFORM = 0, normal completion with ERROR < EPS;
*            if INFORM = 1, completion with ERROR > EPS and MAXPTS 
*                           function vaules used; increase MAXPTS to 
*                           decrease ERROR;
*            if INFORM = 2, N > 100 or N < 1.
*
      INTEGER :: N, MAXPTS, MPT, INFORM, IVLS
      DOUBLE PRECISION :: ABSEPS, RELEPS, ERROR, VALUE, EPS
      INTERFACE
         DOUBLE PRECISION FUNCTION FUNCTN(N,Z)
         DOUBLE PRECISION,DIMENSION(*), INTENT(IN) :: Z
         INTEGER, INTENT(IN) :: N
         END FUNCTION FUNCTN
      END INTERFACE
      INFORM=0
      IF ( N .GT. NDIMMAX .OR. N .LT. 1 ) THEN
         INFORM = 2
         VALUE = 0.d0
         ERROR = 1.d0
         RETURN
      ENDIF
*
*        Call then Monte-Carlo integration subroutine
*
      MPT = 25 + 10*N
      CALL RCRUDE(N, MPT, FUNCTN, ERROR, VALUE, 0)
      IVLS = MPT
 10   EPS = MAX( ABSEPS, RELEPS*ABS(VALUE) )
      IF ( ERROR .GT. EPS .AND. IVLS .LT. MAXPTS ) THEN 
         MPT = MAX( MIN( INT(MPT*(ERROR/(EPS))**2), 
     &        MAXPTS-IVLS ), 10 )
         CALL RCRUDE(N, MPT, FUNCTN, ERROR, VALUE, 1)
         IVLS = IVLS + MPT
         GO TO 10
      ENDIF
      IF ( ERROR. GT. EPS .AND. IVLS .GE. MAXPTS ) INFORM = 1
      !IF (INFORM.EQ.1) print *,'ranmc eps',EPS 
      END SUBROUTINE RANMC
      SUBROUTINE RANLHMC( N, MAXPTS,NLHD,USEMIDP,MLHD,FUNCTN, ABSEPS, 
     &     RELEPS, ERROR, VALUE, INFORM )
      IMPLICIT NONE
*
*     A subroutine for computing multivariate integrals.
*  This subroutine computes an approximation to the integral
*
*      1 1     1
*     I I ... I       FUNCTN(NDIM,X)  dx(NDIM)...dx(2)dx(1)
*      0 0     0  
*
***************  Parameters for RANLHMC  ********************************
*
****** Input Parameters
*
*     N      INTEGER, the number of variables.
*     MAXPTS INTEGER, maximum number of function values allowed. This 
*            parameter can be used to limit the time taken. A 
*            sensible strategy is to start with MAXPTS = 1000*N, and then
*            increase MAXPTS if ERROR is too large.
*     NLHD   INTEGER, size of Latin Hypercube used. A reasonable value is
*            to start with NLHD = 5*N
*     MLHD   LOGICAL, If true : Modify the Latin Hypercube so that the rank correlation
*            between columns are smaller. 
*    useMIDP LOGICAL, If true : use midpoint of cell instead randomly within cell
*     ABSEPS REAL absolute error tolerance.
*     RELEPS REAL relative error tolerance.
*
****** Output Parameters
*
*     ERROR  REAL estimated absolute error, with 99% confidence level.
*     VALUE  REAL estimated value for the integral
*     INFORM INTEGER, termination status parameter:
*            if INFORM = 0, normal completion with ERROR < EPS;
*            if INFORM = 1, completion with ERROR > EPS and MAXPTS 
*                           function vaules used; increase MAXPTS to 
*                           decrease ERROR;
*            if INFORM = 2, N > 100 or N < 1.
*
      INTEGER :: N, MAXPTS, NLHD, MPT, INFORM, IVLS
      DOUBLE PRECISION :: ABSEPS, RELEPS, ERROR, VALUE, EPS
      LOGICAL :: MLHD           ! If true : Modified Latin Hypercube Design 
      LOGICAL :: useMIDP        ! If true : use midpoint of cell instead randomly within cell
      INTERFACE
         DOUBLE PRECISION FUNCTION FUNCTN(N,Z)
         DOUBLE PRECISION,DIMENSION(*), INTENT(IN) :: Z
         INTEGER, INTENT(IN) :: N
         END FUNCTION FUNCTN
      END INTERFACE
      INFORM=0
      IF ( N .GT. NDIMMAX .OR. N .LT. 1 ) THEN
         INFORM = 2
         VALUE = 0.d0
         ERROR = 1.d0
!         PRINT *,'RANLHMC N NDIMMAX', N, NDIMMAX
         RETURN
      ENDIF
*
*        Call then Latin Hypercube Monte-Carlo integration subroutine
*
      MPT = 25+MAX(6*NLHD, 10*N)
!      PRINT *,'RANLHMC MPT,NLHD', MPT,NLHD
      CALL RLHCRUDE(N, MPT,NLHD,FUNCTN, ERROR, VALUE, 0, USEMIDP,MLHD)
!      PRINT *,'RANLHMC VALUE', VALUE
      IVLS = MPT
 10   EPS = MAX( ABSEPS, RELEPS*ABS(VALUE) )
      IF ( ERROR .GT. EPS .AND. IVLS .LT. MAXPTS ) THEN 
         MPT = MAX( MIN( INT(MPT*(ERROR/(EPS))**2), 
     &        MAXPTS-IVLS ), MAX(10,NLHD*2) )
         CALL RLHCRUDE(N, MPT,NLHD,FUNCTN, ERROR, VALUE, 1,
     &        USEMIDP,MLHD)
         IVLS = IVLS + MPT
         GO TO 10
      ENDIF
      IF ( ERROR. GT. EPS .AND. IVLS .GE. MAXPTS ) INFORM = 1
      END SUBROUTINE RANLHMC
      SUBROUTINE RCRUDE(NDIM, MAXPTS, FUNCTN, ABSEST, FINEST, IR)
      IMPLICIT NONE
*
*     Crude Monte-Carlo Algorithm with simple antithetic variates
*      and weighted results on restart
*
      INTEGER :: NDIM, MAXPTS, M,  IR, NPTS
      DOUBLE PRECISION :: FINEST, ABSEST, FUN, 
     &     VARSQR, VAREST, VARPRD, FINDIF, FINVAL
      DOUBLE PRECISION, DIMENSION(NDIMMAX) :: X
      INTERFACE
         DOUBLE PRECISION FUNCTION FUNCTN(N,Z)
         DOUBLE PRECISION,DIMENSION(*), INTENT(IN) :: Z
         INTEGER, INTENT(IN) :: N
         END FUNCTION FUNCTN
      END INTERFACE
      SAVE VAREST
      IF ( IR .LE. 0 ) THEN
         VAREST = 0.d0
         FINEST = 0.d0
      ENDIF
      FINVAL = 0.d0
      VARSQR = 0.d0
      NPTS = INT(MAXPTS/2)
      DO M = 1,NPTS
         CALL random_number(X(1:NDIM))
         FUN = FUNCTN(NDIM, X(1:NDIM))
         X(1:NDIM) = 1.d0 - X(1:NDIM)
         FUN = (FUNCTN(NDIM, X(1:NDIM)) + FUN )*0.5d0
         FINDIF = ( FUN - FINVAL )/DBLE(M)
         VARSQR = DBLE( M - 2 )*VARSQR/DBLE(M) + FINDIF*FINDIF 
         FINVAL = FINVAL + FINDIF
      END DO
      VARPRD = VAREST*VARSQR
      FINEST = FINEST + ( FINVAL - FINEST )/(1.d0 + VARPRD)
      IF ( VARSQR .GT. 0 ) VAREST = (1.d0 + VARPRD)/VARSQR
      ABSEST = 3.d0*SQRT( VARSQR/( 1.d0 + VARPRD ) )
      MAXPTS=2*NPTS
      END SUBROUTINE RCRUDE


      SUBROUTINE RLHCRUDE(NDIM,MAXPTS,NLHD,FUNCTN, ABSEST, FINEST, IR,
     &     USEMIDP,MLHD)
      IMPLICIT NONE
*
*     Crude Latin Hypercube Monte-Carlo Algorithm with simple antithetic variates
*      and weighted results on restart
*
      INTEGER :: NDIM, NLHD, MAXPTS, M, K,IX, IR, NPTS
      DOUBLE PRECISION :: FINEST, ABSEST, FUN, 
     &     VARSQR, VAREST, VARPRD, FINDIF, FINVAL
      LOGICAL :: MLHD           ! Modified Latin Hypercube Design 
      LOGICAL :: useMIDP        ! use midpoint of cell instead randomly within cell
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: X,V,C
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: D
      INTEGER, DIMENSION(:,:) , ALLOCATABLE :: LHD
      INTERFACE
         DOUBLE PRECISION FUNCTION FUNCTN(N,Z)
         DOUBLE PRECISION,DIMENSION(*), INTENT(IN) :: Z
         INTEGER, INTENT(IN) :: N
         END FUNCTION FUNCTN
      END INTERFACE
      SAVE VAREST
      IF ( IR .LE. 0 ) THEN
         VAREST = 0.d0
         FINEST = 0.d0
      ENDIF
      FINVAL = 0.d0
      VARSQR = 0.d0
      IF (2*Nlhd.LE.NDIM*3) MLHD=.FALSE.
      ALLOCATE(X(NLHD,NDIM))
      IF (NLHD.GT.1) THEN
         allocate(LHD(1:Nlhd,1:Ndim)) ! allocate LHD
         LHD(:,1)=(/ (K,K=1,Nlhd)/)
         do K=2,Ndim
            LHD(:,K)=LHD(:,1) 
         enddo
         allocate(C(1:Ndim,1:Ndim))    
         if (MLHD) then
            allocate(D(1:Ndim))
            allocate(V(1:Ndim,1:Ndim))
         endif
      ELSE                      ! make sure 
          MLHD=.FALSE.
                                !print * ,'rindscis: only able to use useMIDP if NLHD>1'
          usemidp=.FALSE.
      end if
   
      
      NPTS = INT(MAXPTS/(2*Nlhd))
!      PRINT *, 'RLHCRUDE NPTS',NPTS
      DO M = 1,NPTS
         CALL random_number(X)
      
         if (Nlhd.gt.1) then
            IF (NDIM.EQ.1) GOTO 30
            do IX=1,3            ! do 3 attemps to construct a LHD with rank=Ndim
               do K=2,Ndim
                  !CALL sortre(LHD(:,K),X(:,K)) ! lhd = latin hypercube design
                  !PRINT *,'RLHCRUDE X',X(:,k)
                  CALL binsort(LHD(:,K),X(:,K)) 
               enddo
               if (IX.EQ.3 .AND..NOT.MLHD) GO TO 30
               CALL spearcorr(C,lhd) ! find rankcorrelation between columns
               if (IX.EQ.3) goto 20
               do K=1,Ndim-1    ! see if rank=Ndim
                  if (any(abs(C(K,K+1:Ndim)).GE.1.d0))  then
                     CALL random_number(X)
                     GO TO 10
                  endif
               enddo
               GO TO 20
 10         enddo
 20         if (MLHD) then      !modify lhd by reducing correlation between columns
               DO K=1,Ndim
                  C(K,K)=C(K,K)+1.D-12 ! add nugget effect to ensure that 
                                !inversion is not corrupted by round off errors
               enddo
               CALL svdcmp(C,D,V) ! C=U*D*V'=Q*Q'
               do K=1,Ndim
                  V(K,:)=C(:,K)*sqrt(1/D(K)) ! inverting Q=U*sqrt(D)
               enddo
               
               X=MATMUL(DBLE(LHD),V) ! LHD*inv(Q)
               do K=1,Ndim
                  LHD(:,K)=(/ (IX,IX=1,Nlhd)/)
                  CALL sortre(LHD(:,K),X(:,K)) ! lhd = latin hypercube design
                  !CALL binsort(LHD(:,K),X(:,K)) 
               enddo
            endif 
 30         IF (USEMIDP) then   ! use the center of the cell
               X=(DBLE(LHD)-0.5d0)/DBLE(Nlhd)   
            else                ! distribute uniformly within the cell
               CALL random_number(X)
               X=(DBLE(LHD)-X)/DBLE(Nlhd) 
            endif
         ENDIF
         FUN=0.D0   
         DO IX = 1,NLHD
            FUN = FUN+FUNCTN(NDIM, X(IX,1:NDIM))
            X(IX,1:NDIM) = 1.d0 - X(IX,1:NDIM)
            FUN = ( FUNCTN(NDIM, X(IX,1:NDIM)) + FUN )
         END DO
         FUN=FUN/DBLE(2*NLHD)
!         PRINT *,'RLHCRUDE X=',X
!         PRINT *,'RLHCRUDE FUN=',FUN
         FINDIF = ( FUN - FINVAL )/DBLE(M)
         VARSQR = DBLE( M - 2 )*VARSQR/DBLE(M) + FINDIF**2.d0 
         FINVAL = FINVAL + FINDIF
      END DO
      
      VARPRD = VAREST*VARSQR
      FINEST = FINEST + ( FINVAL - FINEST )/(1.d0 + VARPRD)
      IF ( VARSQR .GT. 0 ) VAREST = (1.d0 + VARPRD)/VARSQR
      ABSEST = 3.d0*SQRT( VARSQR/( 1.d0 + VARPRD ) )
      
      if (ALLOCATED(X))  DEALLOCATE(X)
      if (allocated(lhd)) then
         DEallocate(lhd)
         DEallocate(C)
         if (allocated(D)) then
            DEallocate(D)
            DEallocate(V)
         endif
      endif
      END SUBROUTINE RLHCRUDE

      SUBROUTINE BINSORT(indices,rarray)
      IMPLICIT NONE
      TYPE ENTRY
         DOUBLE PRECISION, POINTER :: VAL
         INTEGER :: IX
         TYPE( ENTRY), POINTER :: NEXT
      END TYPE ENTRY 
      DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: rarray
      INTEGER,          DIMENSION(:), INTENT(inout)  :: indices
      DOUBLE PRECISION, DIMENSION(SIZE(rarray)),TARGET  :: A
      TYPE(ENTRY), DIMENSION(:), ALLOCATABLE,TARGET  :: B      
      TYPE(ENTRY), POINTER   :: FIRST,CURRENT

! local variables
      INTEGER  :: i,im,n
      DOUBLE PRECISION :: mx, mn
! Bucket sort: 
! This subroutine sorts the indices according to rarray. The Assumption is that rarray consists of 
! uniformly distributed numbers. If the assumption holds it runs in O(n) time
      n=size(indices)
      IF (n.EQ.1) RETURN
      !indices=(/(i,i=1,n)/)
      mx = MAXVAL(rarray)
      mn = MINVAL(rarray)
      A=(rarray-mn)/(mx-mn)  ! make sure the numbers are between 0 and 1
     
      !print *,'binsort ind=',indices
      !print *,'binsort rar=',rarray
      !print *,'binsort rar=',A
      ALLOCATE(B(0:n-1))
      !IF (ASSOCIATED(B(0)%VAL)) print *,'binsort B(0)=',B(0)%VAL
      DO I=0,n-1
         NULLIFY(B(I)%VAL)
         NULLIFY(B(I)%NEXT)
      ENDDO
      
      DO I=1,n
         IM=min(ABS(FLOOR(n*A(I))),N-1)
         IF (ASSOCIATED(B(IM)%VAL)) THEN  ! insert the new item by insertion sorting
            ALLOCATE(CURRENT)
            IF (A(I).LT.B(IM)%VAL) THEN
              CURRENT = B(IM) 
              B(IM)   = ENTRY(A(I),indices(I),CURRENT)
            ELSE
               FIRST => B(IM)
               DO WHILE(ASSOCIATED(FIRST%NEXT).AND.
     &              FIRST%NEXT%VAL.LT.A(I))
                  FIRST=FIRST%NEXT
               END DO
            
               CURRENT = ENTRY(A(I),indices(I),FIRST%NEXT)
               FIRST%NEXT => CURRENT  
            ENDIF   
         ELSE
            B(IM)%VAL => A(I)
            B(IM)%IX  = indices(I)
         ENDIF
      END DO
      IM=0
      I=0 
      DO WHILE (IM.LT.N .AND. I.LT.N)
         IF (ASSOCIATED(B(I)%VAL)) THEN
            IM=IM+1
            indices(IM)=B(I)%IX
            DO WHILE (ASSOCIATED(B(I)%NEXT)) 
               CURRENT => B(I)%NEXT
               B(I)%NEXT => B(I)%NEXT%NEXT
               IM=IM+1
               indices(IM)=CURRENT%IX
               DEALLOCATE(CURRENT)
            END DO
         ENDIF
         I=I+1
      END DO
      DEALLOCATE(B)
      !print *,'binsort ind=',indices
      RETURN
      END SUBROUTINE BINSORT

      SUBROUTINE SORTRE(indices,rarray)
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(:), INTENT(inout) :: rarray
      INTEGER,          DIMENSION(:), INTENT(inout) :: indices
! local variables
       INTEGER  :: i,im,j,k,m,n
   
! diminishing increment sort as described by
! Donald E. Knuth (1973) "The art of computer programming,",
! Vol. 3, pp 84-  (sorting and searching)
      n=size(indices)
      ! if the below is commented out then assume indices are already initialized
      !indices=(/(i,i=1,n)/)
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
!______________________________________________________

      SUBROUTINE spearcorr(C,D)
      IMPLICIT NONE
      DOUBLE PRECISION, dimension(:,:), INTENT(out) :: C
      integer, dimension(:,:),intent(in) :: D ! rank matrix
      double precision, dimension(:,:),allocatable :: DD !,DDT
      double precision, dimension(:),allocatable :: tmp
      INTEGER             :: N,M,ix,iy      
      DOUBLE PRECISION    :: dN      
! this procedure calculates spearmans correlation coefficient
! between the columns of D 
     
      N=size(D,dim=1);M=SIZE(D,dim=2)
      dN=dble(N)
      allocate(DD(1:N,1:M))
      DD=dble(D)
!      if (.false.) then ! old call
!         allocate(DDt(1:M,1:N))
!         DDT=transpose(DD)
!         C = matmul(DDt,DD)*12.d0/(dn*(dn*dn-1.d0)) 
!         C=(C-3.d0*(dn+1.d0)/(dn-1.d0))
!         deallocate(DDT)
!      else
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
!      endif
      deallocate(DD)
      return
      END SUBROUTINE spearcorr
  
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
!         Print *,'SVDCMP: You must augment  A  with extra zero rows.'
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
!             print *,'SVDCMP: No convergence in 30 iterations'
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
       END MODULE RCRUDEMOD








*  KRBVRCMOD is a module  containing a:
*
*  Automatic Multidimensional Integration Subroutine
*               
*         AUTHOR: Alan Genz
*                 Department of Mathematics
*                 Washington State University
*                 Pulman, WA 99164-3113
*                 Email: AlanGenz@wsu.edu
*
*         Last Change: 5/15/98
* revised pab 10.03.2000
*   - updated to f90 (i.e. changed to assumed shape arrays + changing integers to DBLE)
*   - put it into a module
*
*  KRBVRC computes an approximation to the integral
*
*      1  1     1
*     I  I ... I       F(X)  dx(NDIM)...dx(2)dx(1)
*      0  0     0
*
*
*  KRBVRC uses randomized Korobov rules for the first 20 variables. 
*  The primary references are
*   "Randomization of Number Theoretic Methods for Multiple Integration"
*    R. Cranley and T.N.L. Patterson, SIAM J Numer Anal, 13, pp. 904-14,
*  and 
*   "Optimal Parameters for Multidimensional Integration", 
*    P. Keast, SIAM J Numer Anal, 10, pp.831-838.
*  If there are more than 20 variables, the remaining variables are
*  integrated using Richtmeyer rules. A reference is
*   "Methods of Numerical Integration", P.J. Davis and P. Rabinowitz, 
*    Academic Press, 1984, pp. 482-483.
*   
***************  Parameters for KRBVRC ********************************************
****** Input parameters
*  NDIM    Number of variables, must exceed 1, but not exceed 100
*  MINVLS  Integer minimum number of function evaluations allowed.
*          MINVLS must not exceed MAXVLS.  If MINVLS < 0 then the
*          routine assumes a previous call has been made with 
*          the same integrand and continues that calculation.
*  MAXVLS  Integer maximum number of function evaluations allowed.
*  FUNCTN  EXTERNALly declared user defined function to be integrated.
*          It must have parameters (NDIM,Z), where Z is a real array
*          of dimension NDIM.
*                                     
*  ABSEPS  Required absolute accuracy.
*  RELEPS  Required relative accuracy.
*
****** Output parameters
*
*  MINVLS  Actual number of function evaluations used.
*  ABSERR  Estimated absolute accuracy of FINEST.
*  FINEST  Estimated value of integral.
*  INFORM  INFORM = 0 for normal exit, when 
*                     ABSERR <= MAX(ABSEPS, RELEPS*ABS(FINEST))
*                  and 
*                     INTVLS <= MAXCLS.
*          INFORM = 1 If MAXVLS was too small to obtain the required 
*          accuracy. In this case a value FINEST is returned with 
*          estimated absolute accuracy ABSERR.
************************************************************************
! William H. Press, Saul Teukolsky, 
! William T. Wetterling and Brian P. Flannery (1997)
! "Numerical recipes in Fortran 77", Vol. 1, pp 299--305  (SOBSEQ)

! You may  initialize the random generator before you 
!  call KRBVRC by the following lines:
!
!      call random_seed(SIZE=seed_size) 
!      allocate(seed(seed_size)) 
!      call random_seed(GET=seed(1:seed_size))  ! get current seed
!      seed(1)=seed1                            ! change seed
!      call random_seed(PUT=seed(1:seed_size)) 
!      deallocate(seed)
!
      MODULE KRBVRCMOD
      IMPLICIT NONE
      PRIVATE
      PUBLIC :: KRBVRC
! 
      INTERFACE KRBVRC
      MODULE PROCEDURE KRBVRC
      END INTERFACE
!
      INTERFACE DKSMRC
      MODULE PROCEDURE DKSMRC
      END INTERFACE
!      
      INTERFACE  DKRCHT
      MODULE PROCEDURE  DKRCHT
      END INTERFACE
      
      INTERFACE  SOBSEQ
      MODULE PROCEDURE SOBSEQ
      END INTERFACE
!
      CONTAINS   

!***********************************************************
!    MAIN INTEGRATION ROUTINE KRBVRC
!***********************************************************  

      SUBROUTINE KRBVRC( NDIM, MINVLS, MAXVLS, FUNCTN, ABSEPS, RELEPS,
     &                   ABSERR, FINEST, INFORM )
*
*  Automatic Multidimensional Integration Subroutine
*               
*         AUTHOR: Alan Genz
*                 Department of Mathematics
*                 Washington State University
*                 Pulman, WA 99164-3113
*                 Email: AlanGenz@wsu.edu
*
*         Last Change: 5/15/98
*
*  KRBVRC computes an approximation to the integral
*
*      1  1     1
*     I  I ... I       F(X)  dx(NDIM)...dx(2)dx(1)
*      0  0     0
*
*
*  KRBVRC uses randomized Korobov rules for the first 20 variables. 
*  The primary references are
*   "Randomization of Number Theoretic Methods for Multiple Integration"
*    R. Cranley and T.N.L. Patterson, SIAM J Numer Anal, 13, pp. 904-14,
*  and 
*   "Optimal Parameters for Multidimensional Integration", 
*    P. Keast, SIAM J Numer Anal, 10, pp.831-838.
*  If there are more than 20 variables, the remaining variables are
*  integrated using Richtmeyer rules. A reference is
*   "Methods of Numerical Integration", P.J. Davis and P. Rabinowitz, 
*    Academic Press, 1984, pp. 482-483.
*   
***************  Parameters ********************************************
****** Input parameters
*  NDIM    Number of variables, must exceed 1, but not exceed 100
*  MINVLS  Integer minimum number of function evaluations allowed.
*          MINVLS must not exceed MAXVLS.  If MINVLS < 0 then the
*          routine assumes a previous call has been made with 
*          the same integrand and continues that calculation.
*  MAXVLS  Integer maximum number of function evaluations allowed.
*  FUNCTN  EXTERNALly declared user defined function to be integrated.
*          It must have parameters (NDIM,Z), where Z is a real array
*          of dimension NDIM.
*                                     
*  ABSEPS  Required absolute accuracy.
*  RELEPS  Required relative accuracy.
****** Output parameters
*  MINVLS  Actual number of function evaluations used.
*  ABSERR  Estimated absolute accuracy of FINEST.
*  FINEST  Estimated value of integral.
*  INFORM  INFORM = 0 for normal exit, when 
*                     ABSERR <= MAX(ABSEPS, RELEPS*ABS(FINEST))
*                  and 
*                     INTVLS <= MAXCLS.
*          INFORM = 1 If MAXVLS was too small to obtain the required 
*                  accuracy. In this case a value FINEST is returned with 
*                  estimated absolute accuracy ABSERR.
*          INFORM = 2 If NDIM>100 or NDIM<1
************************************************************************
      INTEGER, INTENT(IN)    :: NDIM,  MAXVLS
      INTEGER, INTENT(INOUT) :: MINVLS
      INTEGER, INTENT(OUT)   :: INFORM
      DOUBLE PRECISION, INTENT(IN)  :: ABSEPS, RELEPS
      DOUBLE PRECISION, INTENT(OUT) :: FINEST, ABSERR 
      INTEGER :: NP,PLIM,NLIM,KLIM,KLIMI,SAMPLS,I,INTVLS,MINSMP,NK
      PARAMETER ( PLIM = 25, NLIM = 100, KLIM = 20, MINSMP = 8 )
      INTEGER , DIMENSION(PLIM) :: P
      INTEGER , DIMENSION(PLIM,KLIM-1) ::  C 
      DOUBLE PRECISION :: DIFINT,FINVAL,VARSQR,VAREST,VARPRD,VALUE,ONE
      PARAMETER ( ONE = 1.D0 )
      DOUBLE PRECISION, DIMENSION(2*NLIM) :: X  = 0.d0
      DOUBLE PRECISION, DIMENSION(KLIM  ) :: VK = 0.d0
      INTERFACE
         DOUBLE PRECISION FUNCTION FUNCTN(N,Z)
         DOUBLE PRECISION,DIMENSION(*), INTENT(IN) :: Z
         INTEGER, INTENT(IN) :: N
         END FUNCTION FUNCTN
      END INTERFACE
      DATA P / 31, 47, 73, 113, 173, 263, 397, 593, 907, 1361,
     &     2053, 3079, 4621, 6947, 10427, 15641, 23473, 35221, 
     &     52837, 79259, 118891, 178349, 267523, 401287, 601942/
      DATA (C( 1,I), I = 1, 19)/       12,      9,      9,
     &      13,     12,     12,     12,     12,     12,     12,     12,
     &      12,      3,      3,      3,     12,      7,      7,     12/
      DATA (C( 2,I), I = 1, 19)/        13,     11,     17,
     &      10,     15,     15,     15,     15,     15,     15,     22,
     &      15,     15,      6,      6,      6,     15,     15,      9/
      DATA (C( 3,I), I = 1, 19)/        27,     28,     10,
     &      11,     11,     20,     11,     11,     28,     13,     13,
     &      28,     13,     13,     13,     14,     14,     14,     14/
      DATA (C( 4,I), I = 1, 19)/        35,     27,     27,
     &      36,     22,     29,     29,     20,     45,      5,      5,
     &       5,     21,     21,     21,     21,     21,     21,     21/
      DATA (C( 5,I), I = 1, 19)/        64,     66,     28,
     &      28,     44,     44,     55,     67,     10,     10,     10,
     &      10,     10,     10,     38,     38,     10,     10,     10/
      DATA (C( 6,I), I = 1, 19)/       111,     42,     54,
     &     118,     20,     31,     31,     72,     17,     94,     14,
     &      14,     11,     14,     14,     14,     94,     10,     10/
      DATA (C( 7,I), I = 1, 19)/       163,    154,     83,
     &      43,     82,     92,    150,     59,     76,     76,     47,
     &      11,     11,    100,    131,    116,    116,    116,    116/
      DATA (C( 8,I), I = 1, 19)/      246,    189,    242,
     &     102,    250,    250,    102,    250,    280,    118,    196,
     &     118,    191,    215,    121,    121,     49,     49,     49/
      DATA (C( 9,I), I = 1, 19)/      347,    402,    322,
     &     418,    215,    220,    339,    339,    339,    337,    218,
     &     315,    315,    315,    315,    167,    167,    167,    167/
      DATA (C(10,I), I = 1, 19)/      505,    220,    601,
     &     644,    612,    160,    206,    206,    206,    422,    134,
     &     518,    134,    134,    518,    652,    382,    206,    158/
      DATA (C(11,I), I = 1, 19)/     794,    325,    960,
     &     528,    247,    247,    338,    366,    847,    753,    753,
     &     236,    334,    334,    461,    711,    652,    381,    381/
      DATA (C(12,I), I = 1, 19)/     1189,    888,    259,
     &    1082,    725,    811,    636,    965,    497,    497,   1490,
     &    1490,    392,   1291,    508,    508,   1291,   1291,    508/
      DATA (C(13,I), I = 1, 19)/     1763,   1018,   1500,
     &     432,   1332,   2203,    126,   2240,   1719,   1284,    878,
     &    1983,    266,    266,    266,    266,    747,    747,    127/
      DATA  (C(14,I), I = 1, 19)/     2872,   3233,   1534,
     &    2941,   2910,    393,   1796,    919,    446,    919,    919,
     &    1117,    103,    103,    103,    103,    103,    103,    103/
      DATA  (C(15,I), I = 1, 19)/    4309,   3758,   4034,
     &    1963,    730,    642,   1502,   2246,   3834,   1511,   1102,
     &    1102,   1522,   1522,   3427,   3427,   3928,    915,    915/
      DATA  (C(16,I), I = 1, 19)/     6610,   6977,   1686,
     &    3819,   2314,   5647,   3953,   3614,   5115,    423,    423,
     &    5408,   7426,    423,    423,    487,   6227,   2660,   6227/
      DATA  (C(17,I), I = 1, 19)/     9861,   3647,   4073,
     &    2535,   3430,   9865,   2830,   9328,   4320,   5913,  10365,
     &    8272,   3706,   6186,   7806,   7806,   7806,   8610,   2563/
      DATA  (C(18,I), I = 1, 19)/   10327,   7582,   7124,
     &    8214,   9600,  10271,  10193,  10800,   9086,   2365,   4409,
     &   13812,   5661,   9344,   9344,  10362,   9344,   9344,   8585/
      DATA (C(19,I), I = 1, 19)/   19540,  19926,  11582,
     &   11113,  24585,   8726,  17218,    419,   4918,   4918,   4918,
     &   15701,  17710,   4037,   4037,  15808,  11401,  19398,  25950/
      DATA  (C(20,I), I = 1, 19)/    34566,   9579,  12654,
     &   26856,  37873,  38806,  29501,  17271,   3663,  10763,  18955,
     &    1298,  26560,  17132,  17132,   4753,   4753,   8713,  18624/
      DATA  (C(21,I), I = 1, 19)/   31929,  49367,  10982,
     &    3527,  27066,  13226,  56010,  18911,  40574,  20767,  20767,
     &    9686,  47603,  47603,  11736,  11736,  41601,  12888,  32948/
      DATA (C(22,I), I = 1, 19)/   40701,  69087,  77576,
     &   64590,  39397,  33179,  10858,  38935,  43129,  35468,  35468,
     &    2196,  61518,  61518,  27945,  70975,  70975,  86478,  86478/
      DATA  (C(23,I), I = 1, 19)/  103650, 125480,  59978,
     &   46875,  77172,  83021, 126904,  14541,  56299,  43636,  11655,
     &   52680,  88549,  29804, 101894, 113675,  48040, 113675,  34987/
      DATA (C(24,I), I = 1, 19)/  165843,  90647,  59925,
     &  189541,  67647,  74795,  68365, 167485, 143918,  74912, 167289,
     &   75517,   8148, 172106, 126159,  35867,  35867,  35867, 121694/
      DATA (C(25,I), I = 1, 19)/  130365, 236711, 110235,
     &  125699,  56483,  93735, 234469,  60549,   1291,  93937, 245291,
     &  196061, 258647, 162489, 176631, 204895,  73353, 172319,  28881/
*
      SAVE P, C, SAMPLS, NP, VAREST
      IF ( NDIM .GT. 100 .OR. NDIM .LT. 1 ) THEN
         INFORM = 2
         FINEST = 0.d0
         ABSERR = 1.d0
         RETURN
      ENDIF
      INFORM = 1
      INTVLS = 0
      KLIMI = KLIM
      IF ( MINVLS .GE. 0 ) THEN
         FINEST = 0.d0
         VAREST = 0.d0
         SAMPLS = MINSMP 
         DO I = 1, PLIM
            NP = I
            IF ( MINVLS .LT. 2*SAMPLS*P(I) ) GO TO 10
         END DO
         SAMPLS = MAX( MINSMP, MINVLS/( 2*P(NP) ) )
      ENDIF
 10   VK(1) = ONE/DBLE(P(NP))
      NK = MIN( NDIM, KLIM )
      DO I = 2, NK
         VK(I) = MOD(DBLE(C(NP,NK-1))*VK(I-1), ONE )
      END DO
      FINVAL = 0.d0
      VARSQR = 0.d0
      DO I = 1, SAMPLS
         CALL DKSMRC( NDIM, KLIMI, VALUE, P(NP), VK, FUNCTN, X )
         DIFINT = ( VALUE - FINVAL )/DBLE(I)
         FINVAL = FINVAL + DIFINT
         VARSQR = DBLE( I - 2 )*VARSQR/DBLE(I) + DIFINT*DIFINT
      END DO
      INTVLS = INTVLS + 2*SAMPLS*P(NP)
      VARPRD = VAREST*VARSQR
      FINEST = FINEST + ( FINVAL - FINEST )/( 1.d0 + VARPRD )
      IF ( VARSQR .GT. 0.d0 ) VAREST = ( 1.d0 + VARPRD )/VARSQR
      ABSERR = 3.d0*SQRT( VARSQR/( 1.d0 + VARPRD ) )
      IF ( ABSERR .GT. MAX( ABSEPS, ABS(FINEST)*RELEPS ) ) THEN
         IF ( NP .LT. PLIM ) THEN
            NP = NP + 1
         ELSE
            SAMPLS = MIN( 3*SAMPLS/2, ( MAXVLS - INTVLS )/( 2*P(NP) ) ) 
            SAMPLS = MAX( MINSMP, SAMPLS )
         ENDIF
         IF ( INTVLS + 2*SAMPLS*P(NP) .LE. MAXVLS ) GO TO 10
      ELSE
         INFORM = 0
      ENDIF
      MINVLS = INTVLS
*
      END  SUBROUTINE KRBVRC
*
      SUBROUTINE DKSMRC( NDIM, KLIM, SUMKRO, PRIME, VK, FUNCTN, X )
      INTEGER, INTENT(IN):: NDIM, KLIM, PRIME
	DOUBLE PRECISION, INTENT(OUT) :: SUMKRO
      DOUBLE PRECISION, DIMENSION(*), INTENT(INOUT) :: VK,X
	INTEGER :: K, J, JP, NK
	DOUBLE PRECISION ::  ONE, XT, MVNUNI
      PARAMETER ( ONE = 1.d0 )
      INTERFACE
         DOUBLE PRECISION FUNCTION FUNCTN(N,Z)
         DOUBLE PRECISION,DIMENSION(*), INTENT(IN) :: Z
         INTEGER, INTENT(IN) :: N
         END FUNCTION FUNCTN
      END INTERFACE      
      SUMKRO = 0.d0
*
*     Randomize Variable Order
*
      NK = MIN( NDIM, KLIM )
      DO J = 1, NK-1
         CALL random_number(MVNUNI)
!         JP = J + NINT(MVNUNI*DBLE( NK + 1 - J )) 
         JP = J + NINT(MVNUNI*DBLE( NK - J ))   ! pab 21.11.2000
         XT = VK(J)
         VK(J) = VK(JP)
         VK(JP) = XT
      END DO
*
*     Determine Random Shifts for each Variable
*
      CALL random_number(X(NDIM+1:2*NDIM))
*
*     Compute periodized and symmetrized  lattice rule sum
*
      DO K = 1, PRIME
         X(1:NK) = MOD( DBLE(K)*VK(1:NK), ONE )
         IF ( NDIM. GT. KLIM ) CALL DKRCHT(KLIM, NDIM-KLIM, X) !X(KLIM+1:NDIM) )
         DO J = 1, NDIM
            XT = X(J) + X(NDIM+J)
            IF ( XT .GT. ONE ) XT = XT - 1.d0
            X(J) = ABS( 2.d0*XT - 1.d0 )
         END DO
         SUMKRO = SUMKRO+(FUNCTN(NDIM,X)-SUMKRO)/DBLE(2*K-1)
         X(1:NDIM) = 1.d0 - X(1:NDIM)
         SUMKRO = SUMKRO+(FUNCTN(NDIM,X)-SUMKRO)/DBLE(2*K)
      END DO
      END  SUBROUTINE DKSMRC 
*
      SUBROUTINE DKRCHT(KLIM, S, QUASI )
*
*     This subroutine generates a new quasi-random Richtmeyer vector. 
*     A reference is
*      "Methods of Numerical Integration", P.J. Davis and P. Rabinowitz, 
*       Academic Press, 1984, pp. 482-483.
*
*       INPUTS:
*      KLIM - Lower start value
*         S - the number of dimensions; 
*             DKRCHT is initialized for each new S or S < 1.
*
*       OUTPUTS:
*         QUASI - a new quasi-random S-vector
*
* revised pab 28.05.2003
* - added klim in order to avoid copying of arrays in and out
* revised pab 01.11.1999
* updated to fortran 90
      INTEGER, INTENT(IN) :: S,KLIM
      DOUBLE PRECISION , DIMENSION(*) :: QUASI
      INTEGER :: MXDIM, MXHSUM, B
      PARAMETER ( MXDIM = 80, MXHSUM = 48, B = 2 )
      INTEGER :: HISUM, I,  OLDS 
      DOUBLE PRECISION , DIMENSION(MXDIM) :: PSQT
      INTEGER, DIMENSION(MXDIM   )  :: PRIME
      INTEGER, DIMENSION(0:MXHSUM)  :: N
        
     
      DOUBLE PRECISION ::  ONE, RN
      PARAMETER ( ONE = 1.D0 )
      PARAMETER ( PRIME = (/ 
     &     2,    3,    5,    7,   11,   13,   17,   19,   23,   29,
     &    31,   37,   41,   43,   47,   53,   59,   61,   67,   71,
     &    73,   79,   83,   89,   97,  101,  103,  107,  109,  113,
     &   127,  131,  137,  139,  149,  151,  157,  163,  167,  173,
     &   179,  181,  191,  193,  197,  199,  211,  223,  227,  229,
     &   233,  239,  241,  251,  257,  263,  269,  271,  277,  281,
     &   283,  293,  307,  311,  313,  317,  331,  337,  347,  349,
     &   353,  359,  367,  373,  379,  383,  389,  397,  401,  409/))
* Primes to continue
* 419  421   431   433   439   443   449   457   461   463   467   479   487   491   499
* 503   509   521   523   541   547   557   563   569   571   577   587   593   599
      SAVE OLDS, PSQT, HISUM, N
      DATA OLDS / 0 /
      IF ( S .NE. OLDS .OR. S .LT. 1 ) THEN                          
         OLDS = ABS(S)                             ! pab 14.03.2000
         N(0) = 0
         HISUM = 0
         DO I = 1, OLDS
            RN = DBLE(PRIME(I))
            PSQT(I) = SQRT( RN )
         END DO
      END IF
      DO I = 0, HISUM 
         N(I) = N(I) + 1
         IF ( N(I) .LT. B ) GO TO 10
         N(I) = 0
      END DO
      HISUM = HISUM + 1
      IF ( HISUM .GT. MXHSUM ) HISUM = 0
      N(HISUM) = 1
 10   RN = 0.d0
      DO I = HISUM, 0, -1
         RN = DBLE(N(I)) + DBLE(B)*RN
      END DO
      DO I = 1, OLDS
         QUASI(KLIM+I) = MOD( RN*PSQT(I), ONE )
      END DO
      END SUBROUTINE DKRCHT
!
! SOBSEQ is not taken in to use:
!
      SUBROUTINE SOBSEQ(N,X)
      IMPLICIT NONE
      DOUBLE PRECISION,DIMENSION(*), INTENT(OUT):: X
      INTEGER, INTENT(IN) :: N
      INTEGER,PARAMETER ::MAXBIT=30,MAXDIM=6
      INTEGER :: I,IM, IN,IPP,J,K,L, OLDN
      INTEGER, DIMENSION(MAXDIM) :: IP,MDEG,IX
      INTEGER, DIMENSION(MAXDIM,MAXBIT) ::IU
      INTEGER, DIMENSION(MAXDIM*MAXBIT) ::IV
      DOUBLE PRECISION :: FAC
      SAVE IP,MDEG,IX,IV,IN,FAC, OLDN
      DATA OLDN / 0 /
      DATA IP /0,1,1,2,1,4 /, MDEG /1,2,3,3,4,4 /
      DATA IX /0,0,0,0,0,0 /
      DATA IV /1,1,1,1,1,1,3,1,3,3,1,1,5,
     &     7,7,3,3,5,15,11,5,15,13,9,156*0/
	!(MAXDIM*MAXBIT-24)
      EQUIVALENCE (IV,IU)       ! to allow both 1D and 2D addressing
! returns sobols sequence of quasi-random numbers between 0 1
! When n is new or is negative, internally initializes a set of MAXBIT
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
      
	
      
      IF (N.LT.0 .OR. OLDN.NE.N ) THEN          ! INITIALIZE, DO NOT RETURN VECTOR
         OLDN = ABS(N)
         IX=0
         IN=0  ! RANDOM STARTPOINT: CALL RANDOM_NUMBER(P); IN=P*2^MAXBIT
               ! AND REMOVE WARNING MESSAGE BELOW
         !IF (IV(1).NE.1) RETURN
 
         IF (IV(1).EQ.1) THEN
            FAC=1.D0/2.D0**MAXBIT
            DO K=1,MAXDIM
               DO J=1,MDEG(K)   ! STORED VALUES NEED NORMALIZATION
                  IU(K,J)=IU(K,J)*2**(MAXBIT-J)
               ENDDO
               DO J=1,MDEG(K)+1,MAXBIT ! USE RECCURENCE TO GET OTHER VALUES
                  IPP=IP(K)
                  I=IU(K,J-MDEG(K))
                  I=IEOR(I,I/2**MDEG(K))
                  DO L=MDEG(K)-1,1,-1
                     IF (IAND(IPP,1).NE.0) I=IEOR(I,IU(K,J-L))
                     IPP=IPP/2
                  ENDDO
                  IU(K,J)=I
               ENDDO      
            ENDDO
         ENDIF
      ENDIF                      ! CALCULATE THE NEXT VECTOR IN THE SEQUENCE
         IM=IN
         DO J=1,MAXBIT          ! FIND THE RIGHTMOST ZERO BIT
            IF (IAND(IM,1).EQ.0) GOTO 1
            IM=IM/2
         ENDDO
!         PRINT *,'MAXBIT TOO SMALL IN SOBSEQ'
 1       IM=(J-1)*MAXDIM
         DO K=1,MIN(OLDN,MAXDIM)   !XOR THE 
            IX(K)=IEOR(IX(K),IV(IM+K))
            X(K)=IX(K)*FAC
         ENDDO
         IN=IN+1                ! INCREMENT COUNTER
      
      RETURN
      END SUBROUTINE SOBSEQ

      END MODULE KRBVRCMOD










*  KROBOVMOD is a module  containing a:
*
*  Automatic Multidimensional Integration Subroutine
*               
*         AUTHOR: Alan Genz
*                 Department of Mathematics
*                 Washington State University
*                 Pulman, WA 99164-3113
*                 Email: AlanGenz@wsu.edu
*
*         Last Change: 4/15/98
*
* revised pab 10.03.2000
*   - updated to f90 (i.e. changed to assumed shape arrays + changing integers to DBLE)
*   - put it into a module
*
*  KROBOV computes an approximation to the integral
*
*      1  1     1
*     I  I ... I       F(X)  dx(NDIM)...dx(2)dx(1)
*      0  0     0
*
*
*  KROBOV uses randomized Korobov rules. The primary references are
*  "Randomization of Number Theoretic Methods for Multiple Integration"
*   R. Cranley and T.N.L. Patterson, SIAM J Numer Anal, 13, pp. 904-14,
*  and 
*   "Optimal Parameters for Multidimensional Integration", 
*    P. Keast, SIAM J Numer Anal, 10, pp.831-838.
*   
***************  Parameters ********************************************
****** Input parameters
*  NDIM    Number of variables, must exceed 1, but not exceed 100
*  MINVLS  Integer minimum number of function evaluations allowed.
*          MINVLS must not exceed MAXVLS.  If MINVLS < 0 then the
*          routine assumes a previous call has been made with 
*          the same integrand and continues that calculation.
*  MAXVLS  Integer maximum number of function evaluations allowed.
*  FUNCTN  EXTERNALly declared user defined function to be integrated.
*          It must have parameters (NDIM,Z), where Z is a real array
*          of dimension NDIM.
*  ABSEPS  Required absolute accuracy.
*  RELEPS  Required relative accuracy.
****** Output parameters
*  MINVLS  Actual number of function evaluations used.
*  ABSERR  Estimated absolute accuracy of FINEST.
*  FINEST  Estimated value of integral.
*  INFORM  INFORM = 0 for normal exit, when 
*                     ABSERR <= MAX(ABSEPS, RELEPS*ABS(FINEST))
*                  and 
*                     INTVLS <= MAXCLS.
*          INFORM = 1 If MAXVLS was too small to obtain the required 
*          accuracy. In this case a value FINEST is returned with 
*          estimated absolute accuracy ABSERR.
************************************************************************
! You may  initialize the random generator before you 
!  call KROBOV by the following lines:
!
!      call random_seed(SIZE=seed_size) 
!      allocate(seed(seed_size)) 
!      call random_seed(GET=seed(1:seed_size))  ! get current seed
!      seed(1)=seed1                            ! change seed
!      call random_seed(PUT=seed(1:seed_size)) 
!      deallocate(seed)


      MODULE KROBOVMOD
      IMPLICIT NONE
      PRIVATE
      PUBLIC :: KROBOV
 
      INTERFACE KROBOV
      MODULE PROCEDURE KROBOV
      END INTERFACE

      INTERFACE KROSUM
      MODULE PROCEDURE KROSUM
      END INTERFACE

      CONTAINS   

!***********************************************************
!    MAIN INTEGRATION ROUTINE KROBOV
!***********************************************************  
      SUBROUTINE KROBOV( NDIM, MINVLS, MAXVLS, FUNCTN, ABSEPS, RELEPS,
     &                   ABSERR, FINEST, INFORM )
*
*  Automatic Multidimensional Integration Subroutine
*               
*         AUTHOR: Alan Genz
*                 Department of Mathematics
*                 Washington State University
*                 Pulman, WA 99164-3113
*                 Email: AlanGenz@wsu.edu
*
*         Last Change: 4/15/98
*
*  KROBOV computes an approximation to the integral
*
*      1  1     1
*     I  I ... I       F(X)  dx(NDIM)...dx(2)dx(1)
*      0  0     0
*
*
*  KROBOV uses randomized Korobov rules. The primary references are
*  "Randomization of Number Theoretic Methods for Multiple Integration"
*   R. Cranley and T.N.L. Patterson, SIAM J Numer Anal, 13, pp. 904-14,
*  and 
*   "Optimal Parameters for Multidimensional Integration", 
*    P. Keast, SIAM J Numer Anal, 10, pp.831-838.
*   
***************  Parameters ********************************************
****** Input parameters
*  NDIM    Number of variables, must exceed 1, but not exceed 100
*  MINVLS  Integer minimum number of function evaluations allowed.
*          MINVLS must not exceed MAXVLS.  If MINVLS < 0 then the
*          routine assumes a previous call has been made with 
*          the same integrand and continues that calculation.
*  MAXVLS  Integer maximum number of function evaluations allowed.
*  FUNCTN  EXTERNALly declared user defined function to be integrated.
*          It must have parameters (NDIM,Z), where Z is a real array
*          of dimension NDIM.
*  ABSEPS  Required absolute accuracy.
*  RELEPS  Required relative accuracy.
****** Output parameters
*  MINVLS  Actual number of function evaluations used.
*  ABSERR  Estimated absolute accuracy of FINEST.
*  FINEST  Estimated value of integral.
*  INFORM  INFORM = 0 for normal exit, when 
*                     ABSERR <= MAX(ABSEPS, RELEPS*ABS(FINEST))
*                  and 
*                     INTVLS <= MAXCLS.
*          INFORM = 1 If MAXVLS was too small to obtain the required 
*                   accuracy. In this case a value FINEST is returned with 
*                   estimated absolute accuracy ABSERR.
*          INFORM = 2 If NDIM>100 or NDIM<1
************************************************************************
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, MAXVLS 
	INTEGER, INTENT(INOUT) ::MINVLS
	INTEGER, INTENT(OUT) ::INFORM
      DOUBLE PRECISION, INTENT(IN)  :: ABSEPS, RELEPS
	DOUBLE PRECISION, INTENT(OUT) :: FINEST, ABSERR
! Local variables:	 
	INTEGER :: NP, PLIM, NLIM, SAMPLS, I, INTVLS, MINSMP
      PARAMETER ( PLIM = 20, NLIM = 100, MINSMP = 6 )
      INTEGER, DIMENSION(PLIM,NLIM) :: C
      INTEGER, DIMENSION(PLIM)      :: P
      DOUBLE PRECISION :: DIFINT, FINVAL, VARSQR, VAREST, VARPRD, VALUE  
      DOUBLE PRECISION, DIMENSION(NLIM) :: ALPHA, X, VK
      DOUBLE PRECISION :: ONE
      PARAMETER ( ONE = 1.d0 )
      INTERFACE
         DOUBLE PRECISION FUNCTION FUNCTN(N,Z)
         DOUBLE PRECISION,DIMENSION(*), INTENT(IN) :: Z
         INTEGER, INTENT(IN) :: N
         END FUNCTION FUNCTN
      END INTERFACE
      DATA P /113, 173, 263,397,593,907,1361,2053,3079,4621,6947,
     &     10427, 15641,23473, 35221, 52837, 79259,
     &     118891, 178349, 267523 /
        DATA ( C( 1,I), I = 1, 99 ) /   
     &     42,    54,    55,    32,    13,    26,    26,    13,    26,
     &     14,    13,    26,    35,     2,     2,     2,     2,    56,
     &     28,     7,     7,    28,     4,    49,     4,    40,    48,
     &      5,    35,    27,    16,    16,     2,     2,     7,    28,
     &      4,    49,     4,    56,     8,     2,     2,    56,     7,
     &     16,    28,     7,     7,    28,     4,    49,     4,    37,
     &     55,    21,    33,    40,    16,    16,    28,     7,    16,
     &     28,     4,    49,     4,    56,    35,     2,     2,     2,
     &     16,    16,    28,     4,    16,    28,     4,    49,     4,
     &     40,    40,     5,    42,    27,    16,    16,    28,     4,
     &     16,    28,     4,    49,     4,     8,     8,     2,     2/
      DATA  ( C( 2,I), I = 1, 99 ) /    
     &     64,    34,    57,     9,    72,    86,    16,    75,    75,
     &     70,    42,     2,    86,    62,    62,    30,    30,     5,
     &     42,    70,    70,    70,    53,    70,    70,    53,    42,
     &     62,    53,    53,    53,    69,    75,     5,    53,    86,
     &      2,     5,    30,    75,    59,     2,    69,     5,     5,
     &     63,    62,     5,    69,    30,    44,    30,    86,    86,
     &      2,    69,     5,     5,     2,     2,    61,    69,    17,
     &      2,     2,     2,    53,    69,     2,     2,    86,    69,
     &     13,     2,     2,    37,    43,    65,     2,     2,    30,
     &     86,    45,    16,    32,    18,    86,    86,    86,     9,
     &     63,    63,    11,    76,    76,    76,    63,    60,    70/
      DATA  ( C( 3,I), I = 1, 99 ) /   
     &    111,    67,    98,    36,    48,   110,     2,   131,     2,
     &      2,   124,   124,    48,     2,     2,   124,   124,    70,
     &     70,    48,   126,    48,   126,    56,    65,    48,    48,
     &     70,     2,    92,   124,    92,   126,   131,   124,    70,
     &     70,    70,    20,   105,    70,     2,     2,    27,   108,
     &     27,    39,     2,   131,   131,    92,    92,    48,     2,
     &    126,    20,   126,     2,     2,   131,    38,   117,     2,
     &    131,    68,    58,    38,    90,    38,   108,    38,     2,
     &    131,   131,   131,    68,    14,    94,   131,   131,   131,
     &    108,    18,   131,    56,    85,   117,   117,     9,   131,
     &    131,    55,    92,    92,    92,   131,   131,    48,    48/
      DATA  ( C( 4,I), I = 1, 99 ) /    
     &    151,   168,    46,   197,    69,    64,     2,   198,   191,
     &    134,   134,   167,   124,    16,   124,   124,   124,   124,
     &    141,   134,   128,     2,     2,    32,    32,    32,    31,
     &     31,    64,    64,    99,     4,     4,   167,   124,   124,
     &    124,   124,   124,   124,   107,    85,    79,    85,   111,
     &     85,   128,    31,    31,    31,    31,    64,   167,     4,
     &    107,   167,   124,   124,   124,   124,   124,   124,   107,
     &    183,     2,     2,     2,    62,    32,    31,    31,    31,
     &     31,    31,   167,     4,   107,   167,   124,   124,   124,
     &    124,   124,   124,   107,   142,   184,   184,    65,    65,
     &    183,    31,    31,    31,    31,    31,   167,     4,   107/
      DATA  ( C( 5,I), I = 1, 99 ) /   
     &    229,    40,   268,    42,   153,   294,    71,     2,   130,
     &    199,   199,   199,   149,   199,   149,   153,   130,   149,
     &    149,    15,   119,   294,    31,    82,   260,   122,   209,
     &    209,   122,   296,   130,   130,   260,   260,    30,   206,
     &     94,   209,    94,   122,   209,   209,   122,   122,   209,
     &    130,     2,   130,   130,    38,    38,    79,    82,    94,
     &     82,   122,   122,   209,   209,   122,   122,   168,   220,
     &     62,    60,   168,   282,   282,    82,   209,   122,    94,
     &    209,   122,   122,   122,   122,   258,   148,   286,   256,
     &    256,    62,    62,    82,   122,    82,    82,   122,   122,
     &    122,   209,   122,    15,    79,    79,    79,    79,   168/
      DATA  ( C( 6,I), I = 1, 99 ) /   
     &    264,   402,   406,   147,   452,   153,   224,     2,     2,
     &    224,   224,   449,   101,   182,   449,   101,   451,   181,
     &    181,   101,   101,   377,    85,   453,   453,   453,    85,
     &    197,   451,     2,     2,   101,   449,   449,   449,   173,
     &    173,     2,   453,   453,     2,   426,    66,   367,   426,
     &    101,   453,     2,    32,    32,    32,   101,     2,     2,
     &    453,   223,   147,   449,   290,     2,   453,     2,    83,
     &    223,   101,   453,     2,    83,    83,   147,     2,   453,
     &    147,   147,   147,   147,   147,   147,   147,   453,   153,
     &    153,   147,     2,   224,   290,   320,   453,   147,   431,
     &    383,   290,   290,     2,   162,   162,   147,     2,   162/
      DATA ( C( 7,I), I = 1, 99 ) /   
     &    505,   220,   195,   410,   199,   248,   460,   471,     2,
     &    331,   662,   547,   209,   547,   547,   209,     2,   680,
     &    680,   629,   370,   574,    63,    63,   259,   268,   259,
     &    547,   209,   209,   209,   547,   547,   209,   209,   547,
     &    547,   108,    63,    63,   108,    63,    63,   108,   259,
     &    268,   268,   547,   209,   209,   209,   209,   547,   209,
     &    209,   209,   547,   108,    63,    63,    63,   405,   285,
     &    234,   259,   259,   259,   259,   209,   209,   209,   209,
     &    209,   209,   209,   209,   547,   289,   289,   234,   285,
     &    316,     2,   410,   259,   259,   259,   268,   209,   209,
     &    209,   209,   547,   547,   209,   209,   209,   285,   316/
      DATA ( C( 8,I), I = 1, 99 ) /   
     &    468,   635,   849,   687,   948,    37,  1014,   513,     2,
     &      2,     2,     2,     2,  1026,     2,     2,  1026,   201,
     &    201,     2,  1026,   413,  1026,  1026,     2,     2,   703,
     &    703,     2,     2,   393,   393,   678,   413,  1026,     2,
     &      2,  1026,  1026,     2,   405,   953,     2,  1026,   123,
     &    123,   953,   953,   123,   405,   794,   123,   647,   613,
     &   1026,   647,   768,   953,   405,   953,   405,   918,   918,
     &    123,   953,   953,   918,   953,   536,   405,    70,   124,
     &   1005,   529,   207,   405,   405,   953,   953,   123,   918,
     &    918,   953,   405,   918,   953,   468,   405,   794,   794,
     &    647,   613,   548,   405,   953,   405,   953,   123,   918/
      DATA ( C( 9,I), I = 1, 99 ) /   
     &   1189,  1423,   287,   186,   341,    77,   733,   733,  1116,
     &      2,  1539,     2,     2,     2,     2,     2,  1116,   847,
     &   1174,     2,   827,   713,   910,   944,   139,  1174,  1174,
     &   1539,  1397,  1397,  1174,   370,    33,  1210,     2,   370,
     &   1423,   370,   370,  1423,  1423,  1423,   434,  1423,   901,
     &    139,  1174,   427,   427,   200,  1247,   114,   114,  1441,
     &    139,   728,  1116,  1174,   139,   113,   113,   113,  1406,
     &   1247,   200,   200,   200,   200,  1247,  1247,    27,   427,
     &    427,  1122,  1122,   696,   696,   427,  1539,   435,  1122,
     &    758,  1247,  1247,  1247,   200,   200,   200,  1247,   114,
     &     27,   118,   118,   113,   118,   453,   453,  1084,  1406/
      DATA ( C(10,I), I = 1, 99 ) /   
     &   1764,  1349,  1859,   693,    78,   438,   531,    68,  2234,
     &   2310,  2310,  2310,     2,  2310,  2310,  2102,  2102,   178,
     &    314,   921,  1074,  1074,  1074,  2147,   314,  1869,   178,
     &    178,  1324,  1324,   510,  2309,  1541,  1541,  1541,  1541,
     &    342,  1324,  1324,  1324,  1324,   510,   570,   570,  2197,
     &    173,  1202,   998,  1324,  1324,   178,  1324,  1324,  1541,
     &   1541,  1541,   342,  1541,   886,   178,  1324,  1324,  1324,
     &    510,   784,   784,   501,   652,  1541,  1541,  1324,   178,
     &   1324,   178,  1324,  1541,   342,  1541,  2144,   784,  2132,
     &   1324,  1324,  1324,  1324,   510,   652,  1804,  1541,  1541,
     &   1541,  2132,  1324,  1324,  1324,   178,   510,  1541,   652/
      DATA  ( C(11,I), I = 1, 99 ) /   
     &   2872,  1238,   387,  2135,   235,  1565,   221,  1515,  2950,
     &    486,  3473,     2,  2950,   982,  2950,  3122,  2950,  3172,
     &   2091,  2091,     9,  3449,  3122,  2846,  3122,  3122,  1947,
     &   2846,  3122,   772,  1387,  2895,  1387,     3,     3,     3,
     &   1320,  1320,  2963,  2963,  1320,  1320,  2380,   108,  1284,
     &    702,  1429,   907,  3220,  3125,  1320,  2963,  1320,  1320,
     &   2963,  1320,  1639,  3168,  1660,  2895,  2895,  2895,  2895,
     &   1639,  1297,  1639,   404,  3168,  2963,  2943,  2943,   550,
     &   1387,  1387,  2895,  2895,  2895,  1387,  2895,  1387,  2895,
     &   1320,  1320,  2963,  1320,  1320,  1320,  2963,  1320,     2,
     &   3473,     2,  3473,   772,  2550,     9,  1320,  2963,  1320/
      DATA ( C(12,I), I = 1, 99 ) /  
     &   4309,  2339,  4154,  4480,  4967,   630,  5212,  2592,  4715,
     &   1808,  1808,  5213,     2,   216,  4014,  3499,  3499,  4204,
     &   2701,  2701,  5213,  4157,  1209,  4157,  4460,   335,  4460,
     &   1533,  4575,  4013,  4460,  1881,  2701,  4030,  4030,  1881,
     &   4030,  1738,   249,   335,    57,  2561,  2561,  2561,  1533,
     &   1533,  1533,  4013,  4013,  4013,  4013,  4013,  1533,   856,
     &    856,   468,   468,   468,  2561,   468,  2022,  2022,  2434,
     &    138,  4605,  1100,  2561,  2561,    57,    57,  3249,   468,
     &    468,   468,    57,   468,  1738,   313,   856,     6,  3877,
     &    468,   557,   468,    57,   468,  4605,  2022,     2,  4605,
     &    138,  1100,    57,  2561,    57,    57,  2022,  5213,  3249/
      DATA  ( C(13,I), I = 1, 99 ) /  
     &   6610,  1658,  3022,  2603,  5211,   265,  4985,     3,  4971,
     &   2127,  1877,  1877,     2,  2925,  3175,  3878,  1940,  1940,
     &   1940,  5117,  5117,  5771,  5117,  5117,  5117,  5117,  5117,
     &   5771,  5771,  5117,  3658,  3658,  3658,  3658,  3658,  3658,
     &   5255,  2925,  2619,  1714,  4100,  6718,  6718,  4100,  2322,
     &    842,  4100,  6718,  5119,  4728,  5255,  5771,  5771,  5771,
     &   5117,  5771,  5117,  5117,  5117,  5117,  5117,  5117,  5771,
     &   5771,  1868,  4483,  4728,  3658,  5255,  3658,  5255,  3658,
     &   3658,  5255,  5255,  3658,  6718,  6718,   842,  2322,  6718,
     &   4100,  6718,  4100,  4100,  5117,  5771,  5771,  5117,  5771,
     &   5771,  5771,  5771,  5117,  5117,  5117,  5771,  5771,  1868/
      DATA  ( C(14,I), I = 1, 99 ) /  
     &   9861,  7101,  6257,  7878, 11170, 11638,  7542,  2592,  2591,
     &   6074,  1428,  8925, 11736,  8925,  5623,  5623,  1535,  6759,
     &   9953,  9953, 11459,  9953,  7615,  7615, 11377, 11377,  2762,
     &  11734, 11459,  6892,  1535,  6759,  4695,  1535,  6892,     2,
     &      2,  6892,  6892,  4177,  4177,  6339,  6950,  1226,  1226,
     &   1226,  4177,  6892,  6890,  3640,  3640,  1226, 10590, 10590,
     &   6950,  6950,  6950,  1226,  6950,  6950,  7586,  7586,  7565,
     &   7565,  3640,  3640,  6950,  7565,  6950,  3599,  3599,  3599,
     &   2441,  4885,  4885,  4885,  7565,  7565,  1226,  1226,  1226,
     &   6950,  7586,  1346,  2441,  6339,  3640,  6950, 10590,  6339,
     &   6950,  6950,  6950,  1226,  1226,  6950,   836,  6891,  7565/
      DATA  ( C(15,I), I = 1, 99 ) /  
     &  13482,  5629,  6068, 11974,  4732, 14946, 12097, 17609, 11740,
     &  15170, 10478, 10478, 17610,     2,     2,  7064,  7064,  7064,
     &   5665,  1771,  2947,  4453, 12323, 17610, 14809, 14809,  5665,
     &   5665,  2947,  2947,  2947,  2947, 12323, 12323,  4453,  4453,
     &   2026, 11772,  2026, 11665, 12323, 12323,  3582,  2940,  2940,
     &   6654,  4449,  9254, 11470,   304,   304, 11470,   304, 11470,
     &   6156,  9254, 11772,  6654, 11772,  6156, 11470, 11470, 11772,
     &  11772, 11772, 11470, 11470,   304, 11470, 11470,   304, 11470,
     &    304, 11470,   304,   304,   304,  6654, 11508,   304,   304,
     &   6156,  3582, 11470, 11470, 11470, 17274,  6654,  6654,  6744,
     &   6711,  6654,  6156,  3370,  6654, 12134,  3370,  6654,  3582/
      DATA  ( C(16,I), I = 1, 99 ) /  
     &  13482,  5629,  6068, 11974,  4732, 14946, 12097, 17609, 11740,
     &  15170, 10478, 10478, 17610,     2,     2,  7064,  7064,  7064,
     &   5665,  1771,  2947,  4453, 12323, 17610, 14809, 14809,  5665,
     &   5665,  2947,  2947,  2947,  2947, 12323, 12323,  4453,  4453,
     &   2026, 11772,  2026, 11665, 12323, 12323,  3582,  2940,  2940,
     &   6654,  4449,  9254, 11470,   304,   304, 11470,   304, 11470,
     &   6156,  9254, 11772,  6654, 11772,  6156, 11470, 11470, 11772,
     &  11772, 11772, 11470, 11470,   304, 11470, 11470,   304, 11470,
     &    304, 11470,   304,   304,   304,  6654, 11508,   304,   304,
     &   6156,  3582, 11470, 11470, 11470, 17274,  6654,  6654,  6744,
     &   6711,  6654,  6156,  3370,  6654, 12134,  3370,  6654,  3582/
      DATA  ( C(17,I), I = 1, 99 ) /  
     &  34566, 38838, 23965, 17279, 35325, 33471,   330, 36050, 26419,
     &   3012, 38428, 36430, 36430, 36755, 39629,  5749,  5749, 36755,
     &   5749, 14353, 14353, 14353, 32395, 32395, 32395, 32395, 32396,
     &  32396, 32396, 32396, 27739, 14353, 36430, 36430, 36430, 15727,
     &  38428, 28987, 28987, 27739, 38428, 27739, 18786, 14353, 15727,
     &  28987, 19151, 19757, 19757, 19757, 14353, 22876, 19151, 24737,
     &  24737,  4412, 30567, 30537, 19757, 30537, 19757, 30537, 30537,
     &   4412, 24737, 28987, 19757, 19757, 19757, 30537, 30537, 33186,
     &   4010,  4010,  4010, 17307, 15217, 32789, 37709,  4010,  4010,
     &   4010, 33186, 33186,  4010, 11057, 39388, 33186,  1122, 15089,
     &  39629,     2,     2, 23899, 16466, 16466, 17038,  9477,  9260/
      DATA ( C(18,I), I = 1, 99 ) / 
     &  31929, 40295,  2610,  5177, 17271, 23770,  9140,   952, 39631,
     &      3, 11424, 49719, 38267, 25172,     2,     2, 59445,     2,
     &  59445, 38267, 44358, 14673, 53892, 14674, 14673, 14674, 41368,
     &  17875, 17875, 30190, 20444, 55869, 15644, 25499, 15644, 20983,
     &  44358, 15644, 15644,   485, 41428,   485,   485,   485, 41428,
     &  53798, 50230, 53798, 50253, 50253, 35677, 35677, 17474,  7592,
     &   4098, 17474,   485, 41428,   485, 41428,   485, 41428,   485,
     &  41428, 41428, 41428, 41428, 41428,  9020, 22816,  4098,  4098,
     &   4098,  7592, 42517,   485, 50006, 50006, 22816, 22816,  9020,
     &    485, 41428, 41428, 41428, 41428, 50006,   485, 41428, 41428,
     &  41428, 41428, 22816, 41428, 41428,   485,   485,   485,  9020/
      DATA  ( C(19,I), I = 1, 99 ) / 
     &  73726, 16352, 16297, 74268, 60788,  8555,  1077, 25486, 86595,
     &  59450, 19958, 62205, 62205,  4825,  4825, 89174, 89174, 62205,
     &  19958, 62205, 19958, 27626, 63080, 62205, 62205, 62205, 19958,
     &   8914, 83856, 30760, 47774, 47774, 19958, 62205, 39865, 39865,
     &  74988, 75715, 75715, 74988, 34522, 74988, 74988, 25101, 44621,
     &  44621, 44621, 25101, 25101, 25101, 44621, 47768, 41547, 44621,
     &  10273, 74988, 74988, 74988, 74988, 74988, 74988, 34522, 34522,
     &  67796, 67796, 30208,     2, 67062, 18500, 29251, 29251,     2,
     &  67796, 67062, 38649, 59302,  6225, 67062,  6475,  6225, 46772,
     &  38649, 67062, 46772, 46772, 67062, 46772, 25372, 67062,  6475,
     &  25372, 67062, 67062, 67062,  6225, 67062, 67062, 68247, 80676/
      DATA ( C(20,I), I = 1, 99 )/ 
     & 103650, 50089, 70223, 41805, 74847,112775, 40889, 64866, 44053,
     &   1754,129471, 13630, 53467, 53467, 61378,133761,     2,133761,
     &      2,133761,133761, 65531, 65531, 65531, 38080,133761,133761,
     & 131061,  5431, 65531, 78250, 11397, 38841, 38841,107233,107233,
     & 111286, 19065, 38841, 19065, 19065, 16099,127638, 82411, 96659,
     &  96659, 82411, 96659, 82411, 51986,101677, 39264, 39264,101677,
     &  39264, 39264, 47996, 96659, 82411, 47996, 10971, 10004, 82411,
     &  96659, 82411, 82411, 82411, 96659, 96659, 96659, 82411, 96659,
     &  51986,110913, 51986, 51986,110913, 82411, 54713, 54713, 22360,
     & 117652, 22360, 78250, 78250, 91996, 22360, 91996, 97781, 91996,
     &  97781, 91996, 97781, 97781, 91996, 97781, 97781, 36249, 39779/
      SAVE P, C, SAMPLS, NP, VAREST
      IF ( NDIM .GT. 100 .OR. NDIM .LT. 1 ) THEN
         INFORM = 2
         FINEST = 0.d0
         ABSERR = 1.d0
         RETURN
      ENDIF
      INFORM = 1
      INTVLS = 0
      IF ( MINVLS .GE. 0 ) THEN
         FINEST = 0.d0
         VAREST = 0.d0
         SAMPLS = MINSMP 
         DO I = 1, PLIM
            NP = I
            IF ( MINVLS .LT. 2*SAMPLS*P(I) ) GO TO 10
         END DO
         SAMPLS = MAX( MINSMP, INT(MINVLS/( 2*P(NP)) ) )
      ENDIF
 10   VK(1) = ONE/DBLE(P(NP))
      DO I = 2, NDIM
         VK(I) = MOD( DBLE(C(NP,NDIM-1))*VK(I-1), ONE )
      END DO
      FINVAL = 0.d0
      VARSQR = 0.d0
*
*     Compute mean and standard error for SAMPLS randomized lattice rules
*
      DO I = 1, SAMPLS
         CALL KROSUM( NDIM, VALUE, P(NP), VK, FUNCTN, ALPHA, X )
         DIFINT = ( VALUE - FINVAL )/DBLE(I)
         FINVAL = FINVAL + DIFINT
         VARSQR = DBLE(I - 2)*VARSQR/DBLE(I) + DIFINT*DIFINT
      END DO
      INTVLS = INTVLS + 2*SAMPLS*P(NP)
      VARPRD = VAREST*VARSQR
      FINEST = FINEST + ( FINVAL - FINEST )/( 1.d0 + VARPRD )
      IF ( VARSQR .GT. 0.d0 ) VAREST = ( 1.d0 + VARPRD )/VARSQR
      ABSERR = 3.d0*SQRT( VARSQR/( 1.d0 + VARPRD ) )
      IF ( ABSERR .GT. MAX( ABSEPS, ABS(FINEST)*RELEPS ) ) THEN
         IF ( NP .LT. PLIM ) THEN
            NP = NP + 1
         ELSE
            SAMPLS = MIN( 3*SAMPLS/2, ( MAXVLS - INTVLS )/( 2*P(NP) ) ) 
            SAMPLS = MAX( MINSMP, SAMPLS )
         ENDIF
         IF ( INTVLS + 2*SAMPLS*P(NP) .LE. MAXVLS ) GO TO 10
      ELSE
         INFORM = 0
      ENDIF
      MINVLS = INTVLS
      END SUBROUTINE KROBOV
*
      SUBROUTINE KROSUM( NDIM, SUMKRO, PRIME, VK, FUNCTN, ALPHA, X )
      INTEGER, INTENT(IN):: NDIM, PRIME	
      DOUBLE PRECISION, INTENT(OUT) :: SUMKRO
	DOUBLE PRECISION, DIMENSION(*), INTENT(INOUT) :: ALPHA,X                 ! size NDIM
      INTEGER :: K                 !, J
      DOUBLE PRECISION :: ONE        
      DOUBLE PRECISION, DIMENSION(*), INTENT(IN) :: VK
      INTERFACE
         DOUBLE PRECISION FUNCTION FUNCTN(N,Z)
         DOUBLE PRECISION,DIMENSION(*), INTENT(IN) :: Z
         INTEGER, INTENT(IN) :: N
         END FUNCTION FUNCTN
      END INTERFACE
      PARAMETER ( ONE = 1.d0 )
      SUMKRO = 0.d0
      CALL random_number(ALPHA(1:NDIM))
      DO K = 1, PRIME
         X(1:NDIM) = MOD( DBLE(K)*VK(1:NDIM) + ALPHA(1:NDIM), ONE )
         X(1:NDIM) = ABS( 2.d0*X(1:NDIM) - ONE )
!         PRINT *,'KROSUM W=',X(1:NDIM)
         SUMKRO = SUMKRO+(FUNCTN(NDIM,X)-SUMKRO)/DBLE(2*K-1)
         X(1:NDIM) = ONE - X(1:NDIM)
         SUMKRO = SUMKRO+(FUNCTN(NDIM,X)-SUMKRO)/DBLE(2*K)
      END DO
      END SUBROUTINE KROSUM
      END MODULE KROBOVMOD

