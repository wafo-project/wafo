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
            print *,'mvnfun2 index1',index1
            print *,'SQ'
            call echo(SQ(1:Ntd,1:Nt))
            print *,'error mvnfun2 SQ0,',SQ0
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
