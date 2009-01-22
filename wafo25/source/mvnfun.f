    SUBROUTINE INITFUN
      USE GLOBALDATA, ONLY: EPS2,EPS,xCutOff
      IMPLICIT NONE
! local variables:
      INTEGER ::  N
      INTEGER :: I, J, IJ,  INFA, INFB,Ndleft,LK,JMX,IK
      DOUBLE PRECISION :: TMP, AI, BI, DI, EI,TMPOLD
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
!        The vector must have the length Ndim returned from condsort
! COF  - conditional sorted ChOlesky Factor of the covariance matrix (IN)
! CDI  - Cholesky DIagonal elements used to calculate the mean  
! Cm   - conditional mean of Xd and Xt given Xc, E(Xd,Xt|Xc)
! xd   - variables to the jacobian variable, need no initialization size Nd
! xc   - conditional variables (IN)
! INDEX1 - if INDEX1(I)>Nt then variable no. I is one of the Xd variables otherwise it is one of Xt
      !PRINT *,'Mvnfun,ndim',Ndim
     
      Ndleft = Nd               ! Counter for number of Xd variables left
      N=Nt+Nd-INFIS-1
      INFA = 0
      INFB = 0
      IK = 1                    ! Counter for Ndim 
      IJ = 0
      LK = 0                    ! Counter for LK constraints
      DO I = 1, N+1
         TMP = 0.d0
         JMX=MAX(I-IK,0)
         DO J = 1, I-1-JMX
            IJ = IJ + 1
            TMP = TMP + COF(IJ)*Y(J)
         END DO
         IJ=IJ+JMX
         IF (INFI(I).LT.0) GO TO 100            ! May have infinite int. Limits if Nd>0
         IF ( INFI(I) .NE. 0 ) THEN
            IF ( INFA .EQ. 1 ) THEN
               AI = MAX( AI, A(I) - TMP) !- xCutOff*EPS)  ! added eps to make the int.lim. fuzzy
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
 100     IJ = IJ + 1
         IF (LK.LT.1) THEN
            TMP0 = TMP
         ELSEIF (INDEX1(I).GT.Nt) THEN ! Deterministic variable => replace xd with the mean
            xd(Ndleft) =  Cm(I)+TMP*CDI(I)
            Ndleft = Ndleft-1
         END IF
         IF ( I .EQ. N+1 .OR. COF(IJ+IK+1) .GT. EPS2 ) THEN 
            CALL MVNLMS( AI, BI, 2*INFA+INFB-1, D1, E1 )
            I0=I
            LK0=Lk
            RETURN
         ELSE
            LK=LK+1
         END IF
      END DO
      END SUBROUTINE INITFUN















      FUNCTION MVNFUN( Ndim, W ) RESULT (VAL)
      USE GLOBALDATA, ONLY: EPS2,EPS,xCutOff
      IMPLICIT NONE
      INTEGER, INTENT (in) :: Ndim
      DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: W
      DOUBLE PRECISION :: VAL
! local variables:
      INTEGER ::  N
      INTEGER :: I, J, IJ,  INFA, INFB,Ndleft,LK,JMX,IK
      DOUBLE PRECISION :: TMP, AI, BI, DI, EI,TMPOLD
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
!        The vector must have the length Ndim returned from condsort
! COF  - conditional sorted ChOlesky Factor of the covariance matrix (IN)
! CDI  - Cholesky DIagonal elements used to calculate the mean  
! Cm   - conditional mean of Xd and Xt given Xc, E(Xd,Xt|Xc)
! xd   - variables to the jacobian variable, need no initialization size Nd
! xc   - conditional variables (IN)
! INDEX1 - if INDEX1(I)>Nt then variable no. I is one of the Xd variables otherwise it is one of Xt
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
         !DO J = 1, I-1
         !   IJ = IJ + 1
         !   IF ( J .LT. IK ) TMP = TMP + COF(IJ)*Y(J)
         !END DO
         JMX=MIN(IK,I)-1 
         DO J = 1, JMX
            IJ = IJ + 1
            TMP = TMP + COF(IJ)*Y(J)
         END DO
         IF (IK.LT.I) IJ=IJ+I-IK 
         IF (INFI(I).LT.0) GO TO 100            ! May have infinite int. Limits if Nd>0
         IF ( INFI(I) .NE. 0 ) THEN
            IF ( INFA .EQ. 1 ) THEN
               AI = MAX( AI, A(I) - TMP) !- xCutOff*EPS)  ! added eps to make the int.lim. fuzzy
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
 100     IJ = IJ + 1
         IF (LK.LT.1) TMPOLD=TMP
         IF ( I .EQ. N+1 .OR. COF(IJ+IK+1) .GT. 0.d0 ) THEN 
            CALL MVNLMS( AI, BI, 2*INFA+INFB-1, DI, EI )
            IF ( DI .GE. EI ) THEN
               VAL = 0.d0
               RETURN
            ELSE
               VAL = VAL*( EI - DI )
               IF ( I .LE. N .OR. INDEX1(I-LK).GT.Nt) THEN
                  Y(IK) = FIINV( DI + W(IK)*( EI - DI ) )
                  IF (INDEX1(I-LK).GT.Nt ) THEN
                     xd(Ndleft) = Cm(I-LK)+(TMPOLD+Y(IK))*CDI(I-LK)
                     Ndleft = Ndleft-1
                  ENDIF
               ENDIF
               IK   = IK + 1
               INFA = 0
               INFB = 0
               LK   = 0
            END IF
         ELSE
            LK=LK+1
            IF (INDEX1(I).GT.Nt.AND.LK.GT.1) THEN ! Deterministic variable => replace xd with the mean
               xd(Ndleft) =  Cm(I)+TMP*CDI(I)
               Ndleft = Ndleft-1
            END IF
         END IF
      END DO
      
      IF (Nd.GT.0) THEN
         IF (Ndleft.GT.0.AND.INDEX1(N+1).GT.Nt) THEN
            xd(Ndleft)=Cm(N+1)+TMP*CDI(N+1)
            !print *,'mvnfun IK,Ndleft',IK,Ndleft-1
         ENDIF
         VAL=VAL*jacob(xd,xc)
      ENDIF
      
      RETURN
      END FUNCTION MVNFUN





     function MVNFUN(Ndim,W) RESULT (XIND)
      USE MVNFUNDATA, ONLY : BIG,Cm,xd,xc,Pl1,Pu1
      USE GLOBALDATA, ONLY : Hlo,Hup,xCutOff,Nt,Nd,Nj,Ntd,SQ,
     &     NsXtmj, NsXdj,indXtd,index1,useC1C2,C1C2det,EPS2
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(:  ), INTENT(in)  :: W  
!      DOUBLE PRECISION, DIMENSION(:  ), INTENT(inout)    :: CmN  ! conditional mean
!      DOUBLE PRECISION, DIMENSION(:,:), INTENT(in) :: BIG ! conditional covariance 
!      DOUBLE PRECISION, DIMENSION(:  ), INTENT(inout) :: xd  ! integr. variables
!      DOUBLE PRECISION, DIMENSION(:  ), INTENT(in)    :: xc  ! conditional values
!      DOUBLE PRECISION, INTENT(in) :: Pl1,Pu1 !FI(XMI),FI(XMA) of variable 1
      DOUBLE PRECISION, INTENT(in) :: Ndim
      DOUBLE PRECISION                              :: XIND
!local variables
      DOUBLE PRECISION, DIMENSION(:  ), ALLOCATABLE :: CmN  ! conditional mean
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
! W    - new transformed integration variables, valid range 0..1
!        The vector must have the length Ndim=Nst0+Ntd-Nsd0
! BIG  - conditional sorted covariance matrix (IN)
! Cm   = conditional mean of Xd and Xt given Xc, E(Xd,Xt|Xc)
! CmN  - local conditional mean 
! xd   - variables to the jacobian variable, need no initialization
! xc   - conditional variables (IN) 
! Pl1  = FI(XMI) for the first integration variable (IN)
! Pu1  = FI(XMA) ------||-------------------------------
      ALLOCATE(CmN(1:Ntd))
      CmN(1:Ntd) = Cm(1:Ntd) ! initialize conditional mean
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
                                !print *,'mvnfun,nst,nsd,nd,nj',nst,nsd,Nd,Nj   
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
         Y = Y*jacob ( xd,xc1)   ! jacobian of xd,xc
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
            print *,'mvnfun index1',index1
            print *,'SQ'
            call echo(SQ(1:Ntd,1:Nt))
            print *,'error mvnfun SQ0,',SQ0
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
      IF (ALLOCATED(CmN)) DEALLOCATE(CmN)
      return
      END FUNCTION MVNFUN
