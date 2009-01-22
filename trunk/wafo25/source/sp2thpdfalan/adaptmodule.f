*  ADAPTMOD is a module  containing a:
* 
*   Adaptive Multidimensional Integration Subroutine
*
*   Author: Alan Genz
*           Department of Mathematics
*           Washington State University
*           Pullman, WA 99164-3113 USA
*
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
*  LENWRK  Integer length of real array WORK (working storage); ADAPT
*          needs LENWRK >= 16*NDIM + 27. For maximum efficiency LENWRK
*          should be about 2*NDIM*MAXCLS/RULCLS if MAXCLS FUNCTN
*          calls are needed. If LENWRK is significantly less than this,
*          ADAPT may be less efficient.
*
****** Output Parameters
*
*  MINCLS  Actual number of FUNCTN calls used by ADAPT.
*  WORK    Real array (length LENWRK) of working storage. This contains
*          information that is needed for additional calls of ADAPT
*          using the same integrand (input MINCLS < 0).
*  ABSEST  Real estimated absolute accuracy.
*  FINEST  Real estimated value of integral.
*  INFORM  INFORM = 0 for normal exit, when ABSEST <= ABSREQ or
*                     ABSEST <= |FINEST|*RELREQ with MINCLS <= MAXCLS.
*          INFORM = 1 if MAXCLS was too small for ADAPT to obtain the
*                     result FINEST to within the requested accuracy.
*          INFORM = 2 if MINCLS > MAXCLS, LENWRK < 16*NDIM + 27 or 
*                     RULCLS > MAXCLS.

      MODULE ADAPTMOD
      IMPLICIT NONE
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
      INTEGER ::  NL, LENWRK, RULCLS, TOTCLS, NEWCLS, MAXCLS
      DOUBLE PRECISION, INTENT(IN)  :: ABSEPS, RELEPS
      DOUBLE PRECISION, INTENT(OUT) :: ERROR, VALUE
      DOUBLE PRECISION :: OLDVAL
      PARAMETER ( NL = 20 )
      PARAMETER ( LENWRK = 20*NL**2 )
      DOUBLE PRECISION, DIMENSION(LENWRK) :: WORK
      INTERFACE
         DOUBLE PRECISION FUNCTION FUNCTN(N,Z)
         DOUBLE PRECISION,DIMENSION(:), INTENT(IN) :: Z
         INTEGER, INTENT(IN) :: N
         END FUNCTION FUNCTN
      END INTERFACE
      IF ( N .GT. NL .OR. N .LT. 1 ) THEN
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
     &     LENWRK, WORK, ERROR, VALUE, INFORM )
      MAXCLS = MIN( 10*RULCLS, MAXPTS )
      TOTCLS = 0
      CALL ADAPT(N, TOTCLS, MAXCLS, FUNCTN, ABSEPS, RELEPS, 
     &     LENWRK, WORK, ERROR, VALUE, INFORM)
      IF ( ERROR .GT. MAX( ABSEPS, RELEPS*ABS(VALUE) ) ) THEN
 10      OLDVAL = VALUE
         MAXCLS = MAX( 2*RULCLS,MIN(INT(3*MAXCLS/2),MAXPTS-TOTCLS))
         NEWCLS = -1
         CALL ADAPT(N, NEWCLS, MAXCLS, FUNCTN, ABSEPS, RELEPS, 
     &        LENWRK, WORK, ERROR, VALUE, INFORM)
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
     &     ABSREQ, RELREQ, LENWRK, WORK, ABSEST, FINEST, INFORM)
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
*  LENWRK  Integer length of real array WORK (working storage); ADAPT
*          needs LENWRK >= 16*NDIM + 27. For maximum efficiency LENWRK
*          should be about 2*NDIM*MAXCLS/RULCLS if MAXCLS FUNCTN
*          calls are needed. If LENWRK is significantly less than this,
*          ADAPT may be less efficient.
*
****** Output Parameters
*
*  MINCLS  Actual number of FUNCTN calls used by ADAPT.
*  WORK    Real array (length LENWRK) of working storage. This contains
*          information that is needed for additional calls of ADAPT
*          using the same integrand (input MINCLS < 0).
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
!         EXTERNAL FUNCTN
      INTEGER :: NDIM, MINCLS, MAXCLS, LENWRK, INFORM
      DOUBLE PRECISION  ::  ABSREQ, RELREQ, ABSEST, FINEST      ! ,FUNCTN
      DOUBLE PRECISION, DIMENSION(:) :: WORK                          ! length lenwrk
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: POINTS,WEGHTS,LUM 
      INTEGER :: SBRGNS, MXRGNS, RULCLS, LENRUL, 
     & INERRS, INVALS, INPTRS, INLWRS, INUPRS, INMSHS, INPNTS, INWGTS, 
     & INLOWR, INUPPR, INWDTH, INMESH, INWORK
      INTERFACE
         DOUBLE PRECISION FUNCTION FUNCTN(N,Z)
         DOUBLE PRECISION,DIMENSION(:), INTENT(IN) :: Z
         INTEGER, INTENT(IN) :: N
         END FUNCTION FUNCTN
      END INTERFACE
!      print *,'adapt, ndim', ndim 
      IF ( NDIM .EQ. 1 ) THEN
         LENRUL = 5
         RULCLS = 9
      ELSE IF ( NDIM .LT. 12 ) THEN
         LENRUL = 6
         RULCLS = 2**NDIM + 2*NDIM*(NDIM+2) + 1
      ELSE
         LENRUL = 6
         RULCLS = 1 + 2*NDIM*(1+2*NDIM)
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
          
    
         IF ( MINCLS .LT. 0 ) SBRGNS = WORK(LENWRK)
         ALLOCATE(POINTS(NDIM,LENRUL))
         ALLOCATE(WEGHTS(LENRUL,4))
         ALLOCATE(LUM(NDIM,MXRGNS*3))
         LUM    = reshape(WORK(INLWRS:INWGTS-1),(/ NDIM,MXRGNS*3/))
         WEGHTS = reshape(WORK(INWGTS:INPNTS-1),(/ LENRUL , 4 /))
         POINTS = reshape(WORK(INPNTS:INLOWR-1),(/ NDIM, LENRUL/))
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
      INTEGER NDIM, LENRUL, I, J, NUMNUL, SDIM
      INTEGER, DIMENSION(6) ::  RULPTS
      PARAMETER ( NUMNUL = 4, SDIM = 12 )
      DOUBLE PRECISION , DIMENSION(:,:) :: W
      DOUBLE PRECISION , DIMENSION(:,:) :: G
*      DOUBLE PRECISION W(LENRUL,4), G(NDIM,LENRUL) 
      DOUBLE PRECISION LAM1, LAM2, LAM3, LAMP, RULCON
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
      DO I = 1,LENRUL
         DO J = 1,NDIM
            G(J,I) = 0.d0
         END DO
         DO J = 1,NUMNUL
            W(I,J) = 0.d0
         END DO
      END DO
      RULPTS(5) = 2*NDIM*(NDIM-1)
      RULPTS(4) = 2*NDIM
      RULPTS(3) = 2*NDIM
      RULPTS(2) = 2*NDIM
      RULPTS(1) = 1
      LAMP = 0.85d0
      LAM3 = 0.4707d0
      LAM2 = 4d0/(15.d0 - 5.d0/LAM3)
      W(5,1) = ( 3.d0 - 5.d0*LAM3 )/( 180.d0*(LAM2-LAM3)*LAM2**2 )
      IF ( NDIM .LT. SDIM ) THEN 
         LAM1 = 8.d0*LAM3*(31.d0*LAM3-15.d0)/
     &        ( (3.d0*LAM3-1.d0)*(5.d0*LAM3-3.d0)*35.d0 )
         W(LENRUL,1) = 1.d0/(3.d0*LAM3)**3.d0/2.d0**DBLE(NDIM)
      ELSE
         LAM1 = ( LAM3*(15.d0 - 21.d0*LAM2) + 
     &        35.d0*DBLE(NDIM-1)*(LAM2-LAM3)/9.d0 )
     &       /  ( LAM3*(21.d0 - 35.d0*LAM2) + 
     &        35.d0*DBLE(NDIM-1)*(LAM2/LAM3-1.d0)/9.d0 )
         W(6,1) = 1.d0/(4.d0*(3.d0*LAM3)**3.d0)
      ENDIF
      W(3,1) = ( 15.d0 - 21.d0*(LAM3+LAM1) + 35.d0*LAM3*LAM1 )
     &   /( 210.d0*LAM2*(LAM2-LAM3)*(LAM2-LAM1) ) - 2.d0*(NDIM-1)*W(5,1)
      W(2,1) = ( 15.d0 - 21.d0*(LAM3+LAM2) + 35.d0*LAM3*LAM2 )
     &     /( 210.d0*LAM1*(LAM1-LAM3)*(LAM1-LAM2) )
      IF ( NDIM .LT. SDIM ) THEN
         RULPTS(LENRUL) = 2**NDIM
         LAM3 = SQRT(LAM3)
         DO I = 1,NDIM
            G(I,LENRUL) = LAM3
         END DO
      ELSE
         W(6,1) = 1.d0/(4.d0*(3.d0*LAM3)**3.d0)
         RULPTS(6) = 2*NDIM*(NDIM-1)
         LAM3 = SQRT(LAM3)
         DO I = 1,2
            G(I,6) = LAM3
         END DO
      ENDIF
      IF ( NDIM .GT. 1 ) THEN
         W(5,2) = 1.d0/(6.d0*LAM2)**2 
         W(5,3) = 1.d0/(6.d0*LAM2)**2 
      ENDIF
      W(3,2) = ( 3.d0 - 5.d0*LAM1 )/( 30.d0*LAM2*(LAM2-LAM1) ) 
     &     - 2.d0*(NDIM-1)*W(5,2) 
      W(2,2) = ( 3.d0 - 5.d0*LAM2 )/( 30.d0*LAM1*(LAM1-LAM2) )
      W(4,3) = ( 3.d0 - 5.d0*LAM2 )/( 30.d0*LAMP*(LAMP-LAM2) )
      W(3,3) = ( 3.d0 - 5.d0*LAMP )/( 30.d0*LAM2*(LAM2-LAMP) ) 
     &     - 2.d0*DBLE(NDIM-1)*W(5,3)
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
      RULCON = 2.d0
      CALL RULNRM( LENRUL, NUMNUL, RULPTS, W, RULCON )
      END SUBROUTINE BSINIT
!
!
      SUBROUTINE RULNRM( LENRUL, NUMNUL, RULPTS, W, RULCON )
      IMPLICIT NONE
      INTEGER :: LENRUL, NUMNUL, I, J, K
      INTEGER, DIMENSION(:) :: RULPTS
      DOUBLE PRECISION :: ALPHA, NORMCF, NORMNL,  RULCON
      DOUBLE PRECISION, DIMENSION(:,:) :: W                !(LENRUL, *),
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
            W(I,J) = W(I,J)/RULCON
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
!      EXTERNAL FUNCTN
      INTEGER :: I, J, NDIM, MINCLS, MAXCLS, SBRGNS, MXRGNS, 
     &     RULCLS, LENRUL, INFORM, NWRGNS 
      DOUBLE PRECISION :: ABSREQ, RELREQ, ABSEST, FINEST              !,  FUNCTN,
      DOUBLE PRECISION, DIMENSION(:) ::  ERRORS, VALUES, PONTRS,
     &     LOWER, UPPER, WIDTH, MESH, WORK
      DOUBLE PRECISION, DIMENSION(:,:) ::  WEGHTS, POINTS             ! shape (LENRUL,4) and (Ndim,LENRUL)
      DOUBLE PRECISION, DIMENSION(:,:) ::  LOWERS, UPPERS, MESHES     !SHAPE  (NDIM,MXRGNS), 
      INTEGER :: DIVAXN, TOP, RGNCLS, FUNCLS, DIFCLS      
      INTERFACE
         DOUBLE PRECISION FUNCTION FUNCTN(N,Z)
         DOUBLE PRECISION,DIMENSION(:), INTENT(IN) :: Z
         INTEGER, INTENT(IN) :: N
         END FUNCTION FUNCTN
      END INTERFACE
*
*     Initialization of subroutine
*
!      print *,'ADBASE, ndim', ndim, shape(POINTS) 
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
         CALL DIFFER(NDIM, LOWER, UPPER, WIDTH, WORK(1:NDIM),  
     &        WORK(NDIM+1:2*NDIM), FUNCTN, DIVAXN, DIFCLS)
         FUNCLS = FUNCLS + DIFCLS
         IF ( FUNCLS + RGNCLS*(MESH(DIVAXN)+1.d0)/MESH(DIVAXN)
     &        .LE. MINCLS ) THEN
            RGNCLS = RGNCLS*(MESH(DIVAXN)+1.d0)/MESH(DIVAXN)
            NWRGNS = NWRGNS*(MESH(DIVAXN)+1.d0)/MESH(DIVAXN)
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
            IF ( LOWER(I)+WIDTH(I) .LT. 1 )  GO TO 20
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
            WIDTH(I) = (UPPER(I)-LOWER(I))/(2*MESH(I))
            RGNCLS = RGNCLS*MESH(I)
         END DO
         CALL DIFFER(NDIM, LOWER, UPPER, WIDTH, WORK(1:NDIM),   
     &       WORK(NDIM+1:2*NDIM), FUNCTN, DIVAXN, DIFCLS)
         FUNCLS = FUNCLS + DIFCLS
         RGNCLS = RGNCLS*(MESH(DIVAXN)+1)/MESH(DIVAXN)
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
      INTEGER :: I, LENRUL, NDIM
      DOUBLE PRECISION, DIMENSION(:  ) ::  A, B, WIDTH                       !(NDIM)
      DOUBLE PRECISION, DIMENSION(:,:) :: W                                  !(LENRUL,4), 
      DOUBLE PRECISION, DIMENSION(:,:) :: G                                  !(NDIM,LENRUL), 
      DOUBLE PRECISION, DIMENSION(:  ) :: CENTER, Z                          !(NDIM)
      DOUBLE PRECISION :: RGNERT, BASEST
      DOUBLE PRECISION :: FSYMSM, RGNCMP, RGNVAL,
     &     RGNVOL, RGNCPT, RGNERR
      INTERFACE
         DOUBLE PRECISION FUNCTION FUNCTN(N,Z)
         DOUBLE PRECISION,DIMENSION(:), INTENT(IN) :: Z
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
      RGNERR = SQRT(RGNCMP**2.d0 + RGNERR**2.d0)
      RGNCMP = SQRT(RGNCPT**2.d0 + RGNCMP**2.d0)
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
 
      INTEGER :: S, IXCHNG, LXCHNG, I, L
      DOUBLE PRECISION, DIMENSION(:) :: CENTER, HWIDTH, X, G                     ! shape S
      DOUBLE PRECISION :: INTSUM, GL, GI
      INTERFACE
         DOUBLE PRECISION FUNCTION F(N,Z)
         DOUBLE PRECISION,DIMENSION(:), INTENT(IN) :: Z
         INTEGER, INTENT(IN) :: N
         END FUNCTION F
      END INTERFACE
!      print *,'FULSUM, S', S, shape(X) 
      FULSUM = 0.d0
*
*     Compute centrally symmetric sum for permutation of G
*
 10   INTSUM = 0.d0
      DO I = 1,S
         X(I) = CENTER(I) + G(I)*HWIDTH(I)
      END DO
 20   INTSUM = INTSUM + F(S,X(1:S))
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
      INTEGER :: I, NDIM, DIVAXN, DIFCLS
      DOUBLE PRECISION, DIMENSION(:) :: A, B, WIDTH, Z, DIF        ! (NDIM)
      DOUBLE PRECISION :: FRTHDF, FUNCEN, WIDTHI
      INTERFACE
         DOUBLE PRECISION FUNCTION FUNCTN(N,Z)
         DOUBLE PRECISION,DIMENSION(:), INTENT(IN) :: Z
         INTEGER, INTENT(IN) :: N
         END FUNCTION FUNCTN
      END INTERFACE
!      print *,'DIFFER, ndim', ndim, shape(Z) 
      DIFCLS = 0
      DIVAXN = MOD( DIVAXN, NDIM ) + 1
      IF ( NDIM .GT. 1 ) THEN
         DO I = 1,NDIM 
            DIF(I) = 0.d0
            Z(I) = A(I) + WIDTH(I)
         END DO
!         print *,'Z', Z
 10      FUNCEN = FUNCTN(NDIM, Z(1:NDIM))
         DO I = 1,NDIM
            WIDTHI = 0.2d0*WIDTH(I)
            FRTHDF = 6.d0*FUNCEN
            Z(I) = Z(I) - 4.d0*WIDTHI
            FRTHDF = FRTHDF + FUNCTN(NDIM,Z(1:NDIM))
            Z(I) = Z(I) + 2.d0*WIDTHI
            FRTHDF = FRTHDF - 4.d0*FUNCTN(NDIM,Z(1:NDIM))
            Z(I) = Z(I) + 4*WIDTHI
            FRTHDF = FRTHDF - 4.d0*FUNCTN(NDIM,Z(1:NDIM))
            Z(I) = Z(I) + 2.d0*WIDTHI
            FRTHDF = FRTHDF + FUNCTN(NDIM,Z(1:NDIM))
*     Do not include differences below roundoff
            IF ( FUNCEN + FRTHDF/8.d0 .NE. FUNCEN ) 
     &           DIF(I) = DIF(I) + ABS(FRTHDF)*WIDTH(I)
            Z(I) = Z(I) - 4.d0*WIDTHI
         END DO
         DIFCLS = DIFCLS + 4*NDIM + 1
         DO I = 1,NDIM
            Z(I) = Z(I) + 2*WIDTH(I)
            IF ( Z(I) .LT. B(I) ) GO TO 10
            Z(I) = A(I) + WIDTH(I)
         END DO
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
      INTEGER POINTR, SBRGNS
      DOUBLE PRECISION, DIMENSION(:) :: PONTRS, RGNERS                 !(*)
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
