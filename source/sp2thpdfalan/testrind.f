      PROGRAM TESTRIND
      USE GLOBALDATA
      USE RIND
      IMPLICIT NONE
!
!     Test program for RIND
!
      DOUBLE PRECISION :: XcScale 
      INTEGER :: Nx,Nt,Nj,Nc,Ntd,Ntdc,Ni,Nb,Mb,K,I,J
      INTEGER :: speed1, seed1
      DOUBLE PRECISION, ALLOCATABLE :: BIG(:,:),Xc(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: Blo(:,:),Bup(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: Ex(:),VALS(:)
      DOUBLE PRECISION, ALLOCATABLE :: ERR1(:)
      INTEGER, ALLOCATABLE  :: IndI(:),INFIN(:)
      INTEGER, ALLOCATABLE  :: seed(:)
      INTEGER               :: seed_size
      Ntdc = 20
      Nc = 0
      Ntd = Ntdc - Nc;
      Nt   = Ntd
      Nx   = 1
      Ni   = Ntd + 1
      Nb   = Ntd
      Mb =  1

      
! Set variables in the globaldata module
      SCIS    = 3 
      XcScale = 0.0d0
      ABSEPS  = 1.0d-3
      RELEPS  = 0.0d0
      EPS2    = ABSEPS*1.0d-1
      EPS     = SQRT(EPS2)
      NSIMmax = 1000000
      NSIMmin = 0
      seed1   = 1
      NIT     = Ntdc+1
      xCutOff = 4.5d0

      ALLOCATE(VALS(Nx))
      ALLOCATE(ERR1(Nx))
      ALLOCATE(BIG(Ntdc,Ntdc))
      ALLOCATE(Ex(Ntdc))
      ALLOCATE(INDI(Ni))
      ALLOCATE(Blo(Mb,Nb))
      ALLOCATE(Bup(Mb,Nb))
      ALLOCATE(INFIN(Nb))
      ALLOCATE(Xc(Nc,Nx))
      
      IF (SCIS.gt.0) THEN
         CALL random_seed(SIZE=seed_size) 
         ALLOCATE(seed(seed_size))
                                !print *,'rindinterface seed', seed1
         CALL random_seed(GET=seed(1:seed_size)) ! get current state
         seed(1:seed_size)=seed1 ! change seed
         CALL random_seed(PUT=seed(1:seed_size)) 
         CALL random_seed(GET=seed(1:seed_size)) ! get current state
                                !print *,'rindinterface seed', seed
         DEALLOCATE(seed)
      ENDIF
	PRINT '(''        Test of RIND '')'
      PRINT '(12X, ''Requested Accuracy '',F8.5)', MAX(ABSEPS,RELEPS)
      PRINT '(''           Number of Dimensions is '',I2)', Ntdc
      PRINT '(''     Maximum # of Function Values is '',I7)', NSIMmax
	DO K = 1,10
      CALL RANDOM_NUMBER(BIG)
      DO I=1,Ntdc
         BIG(I,I) = 1.0D0
         DO J=I+1,Ntdc
            BIG(I,J) = 2.0D0*BIG(J,I)-1
         ENDDO
      ENDDO
        
      Ex(:) = 0.0D0
      CALL RANDOM_NUMBER(Blo)
      CALL RANDOM_NUMBER(Bup)
      Blo = -Blo*SQRT(DBLE(Ntdc))
      Bup = Bup*SQRT(DBLE(Ntdc))
      INDI(:) = (/0:Ntd/)
      INFIN(:) = 0
      INFIN(1) = 2

     
      DO J = 1,3
         CALL RINDD(VALS,ERR1,Big,Ex,Xc,Nt,INDI,Blo,Bup,INFIN,XcScale)
       
         PRINT '('' Values:   '',(F11.6))', VALS
         PRINT '('' Error :   '',(F11.6))', ERR1
         INFIN(1) = INFIN(1) - 1
      END DO
      ENDDO
      DEALLOCATE(VALS)
      DEALLOCATE(ERR1)
      DEALLOCATE(BIG)
      DEALLOCATE(Ex)
      DEALLOCATE(Xc)
      DEALLOCATE(INDI)
      DEALLOCATE(Blo)
      DEALLOCATE(Bup)
      DEALLOCATE(INFIN)
      END PROGRAM
