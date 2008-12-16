      PROGRAM rind_interface
C***********************************************************************
C     This program is a interface between matlab and fortran           *
C     for the RIND_D program  ALAN2.f                                  *
C***********************************************************************
C compilation is done by
C f90 -g2 -C -automatic -o ../wave/alpha/rind.exe rind42.f rind_interface.f
Cf90  -gline -Nl126 -C -o ../exec/lnx86/rind.exe  intmoduole.f  alan2.f rindinterf2.f 

      use globaldata, only: SCIS,rateLHD
      use mvnfundata, only : jacobdef
      use rind
   
      IMPLICIT NONE
      double precision, dimension(:,:), allocatable :: BIG,xc,B_up,B_lo
      double precision, dimension(:  ), allocatable :: fxind,ex
      integer ,dimension(:), allocatable :: indI,INFIN
      integer :: k,speed,seed1,Nt,Nd,Nc,NI,Mb,Nx,Nj,NIT,Ntdc,Ntd 
      INTEGER,          DIMENSION(:  ), ALLOCATABLE  :: seed
      INTEGER                                        :: seed_size
      
      
      CALL INIT_LEVELS(speed,Nt,Nd,Nc,NI,Mb,Nx,NIT,Nj,jacobdef,seed1,
     &     SCIS,rateLHD)
      if (SCIS.gt.0) then
         call random_seed(SIZE=seed_size) 
         allocate(seed(seed_size))
         !print *,'rindinterface seed', seed1
         call random_seed(GET=seed(1:seed_size)) ! get current state
         seed(1:seed_size)=seed1          ! change seed
         call random_seed(PUT=seed(1:seed_size)) 
         call random_seed(GET=seed(1:seed_size)) ! get current state
         !print *,'rindinterface seed', seed
         deallocate(seed)
      endif

      Ntd=Nt+Nd;Ntdc=Ntd+Nc
      !print *,'Nt,Nd,Nc,NI,Mb,Nx',Nt,Nd,Nc,NI,Mb,Nx
      !print *,'Speed,NIT',speed,NIT
      !print *,'rind_intNj',Nj
      if (speed.gt.0) CALL INITDATA(speed)  
      !print * ,'speed,Nt',speed,Nt
      !print *,'big '  
      allocate(BIG(1:Ntdc,1:Ntdc))
      allocate(ex(1:Ntdc))
      allocate(B_lo(1:Mb,1:(NI-1)))
      allocate(B_up(1:Mb,1:(NI-1)))
      allocate(indI(NI),INFIN(NI-1))
      allocate(fxind(1:Nx))
      !if (Nc.EQ.0) then
      !   allocate( xc(1:Nc,1:1))
      !   xc(1,1)=1.d0;
      !else
      allocate( xc(1:Nc,1:Nx))
      !end if
      !print *,'big ', shape(BIG),Ntdc
      !print *,'Ex ', size(ex),Ntdc
      !print *,'B_lo ', size(B_lo),(NI-1)*Mb
      !print *,'B_up ', size(B_up)
      !print *,'indI ', size(indI)
      !print *,'fxind ', size(fxind)
      CALL DATA_INPUT(Ntdc,Nc,Nx,NI,Mb,BIG,ex,xc,indI,B_lo,B_up,INFIN)
      !print *,'big '
      IF (Ntd.EQ.indI(NI)) THEN
      
                                !print *,'Ntdc ', Ntdc       
                                !CALL ECHO(BIG(1:Ntdc,1:Ntdc),Ntdc)
         CALL RINDD(fxind,Big(1:Ntdc,1:Ntdc),ex(1:Ntdc),
     *        xc,Nt,indI,B_lo,B_up,INFIN)
         !call echo(big)  
                                !print *,'Ready: ',pt,' of ',Ntime
                                !print *, 'ansr', fxind
         open (unit=11, file='rind.out',  STATUS='unknown')
         do k=1,Nx
            write(11,*)  fxind(k)
         enddo
         close(11)
      ELSE
         PRINT *,'Ntd must be equal to indI(NI)!!, stop'
      ENDIF   
      deallocate(BIG)
      deallocate(ex)
      deallocate(B_lo)
      deallocate(B_up)
      deallocate(indI)
      deallocate(fxind)
      deallocate( xc)
      


      CONTAINS
           
C***********************************************************************
C***********************************************************************
       
      SUBROUTINE data_input(Ntdc,Nc,Nx,NI,Mb,BIG,Ex,xc,indI,B_lo,B_up
     &     ,INFIN) 
      double precision, dimension(:,:), intent(out) :: BIG,xc,B_lo,B_up
      double precision, dimension(:  ), intent(out) :: Ex
      integer, dimension(:  ), intent(out) :: indI,INFIN
      integer, intent(in) :: Ntdc,Nc,Nx,NI,Mb
      integer :: i,j,INFA,INFB
      !print *,'Load big ', shape(BIG),Ntdc
      open (unit=1, file='BIG.in', STATUS='unknown')
      !print *,'Load Ex',shape(Ex)
      open (unit=2, file='Ex.in', STATUS='unknown')
      !print *,'Load xc',shape(xc),Nc,Nx
      open (unit=3, file='xc.in', STATUS='unknown')
      ! print *,'Load indI',shape(indI),NI
      open (unit=4, file='indI.in', STATUS='unknown')
      !print *,'Load B_lo',shape(B_lo),Mb,NI-1
      open (unit=5, file='B_lo.in', STATUS='unknown')
      !print *,'Load B_up',shape(B_up)
      open (unit=8, file='B_up.in', STATUS='unknown')
      !print *,'big '

      do j=1,Ntdc
         do i=1,Ntdc
            read(1,*) BIG(j,i)
         enddo
         read(2,*) Ex(j)
      enddo
      !print *,'big ',BIG(1,1:Ntdc)
      !print *,'big ',BIG(2,1:Ntdc)
      do j=1,Nc
         do i=1,Nx
            read(3,*) xc(j,i)
         enddo
      enddo

      do i=1,NI
         read(4,*) indI(i)
      enddo

      do j=1,Mb
         do i=1,NI-1
            read(5,*) B_lo(j,i)
            read(8,*) B_up(j,i)
         enddo
      enddo
      do i=1,NI-1
         INFA=0
         INFB=0
         if (B_lo(1,i).GT.-10) infa=1
         if (B_up(1,i).LT.10) infb=1
         INFIN(I)=2*infa+infb-1
      enddo
!      print *,'R ',R(1:Ntime)
      
      close(1)
      close(2)
      close(3)
      close(4)
      close(5)
      close(8)
      return
      END SUBROUTINE DATA_INPUT 
      
C***********************************************************************
C***********************************************************************
      SUBROUTINE INIT_LEVELS (speed,Nt,Nd,Nc,NI,Mb,Nx,NIT,Nj,
     &     jacob_def,seed,SCIS,rateLHD)
      USE rind
      integer, intent(out) :: speed,Nt,Nd,Nc,NI,Mb,Nx,NIT,Nj,
     &     jacob_def,seed,SCIS,rateLHD
      double precision :: XSPLT 
      !local variables
      double precision :: dREPS,dEPSS,dEPS2,dXc
      integer :: dNINT,dminQnr,dLe2Qnr

      OPEN(UNIT=14,FILE='sizeinfo.in',STATUS= 'UNKNOWN')
      READ (14,*) speed 
      READ (14,*) Nt
      READ (14,*) Nd
      READ (14,*) Nc
      READ (14,*) NI
      READ (14,*) Mb
      READ (14,*) Nx
      READ (14,*) NIT
      READ (14,*) Nj
      READ (14,*) jacob_def
      READ (14,*) seed
      READ (14,*) SCIS   
      READ (14,*) rateLHD
      READ (14,*) XSPLT
      IF (speed.LT.1) THEN
         READ (14,*) dEPSS
         READ (14,*) dEPS2
         READ (14,*) dXc
         READ (14,*) dREPS
         READ (14,*) dNINT
         READ (14,*) dminQnr
         READ (14,*) dLe2Qnr
         CALL SETDATA (dREPS,dEPSS,dEPS2,dXc) 
      ENDIF
      CLOSE(UNIT=14)
      
      
      RETURN
      END SUBROUTINE INIT_LEVELS
            
C**********************************************************************

 

      
      END  PROGRAM  rind_interface
