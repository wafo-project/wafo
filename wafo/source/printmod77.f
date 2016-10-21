      MODULE PRINTMOD77
      IMPLICIT NONE
!      PUBLIC
! PRINTMOD77 Module containing print functions to screen 

      INTERFACE echo 
      MODULE PROCEDURE echo
      END INTERFACE
     
      INTERFACE PRINTCOF
      MODULE PROCEDURE PRINTCOF
      END INTERFACE 

      INTERFACE PRINTVAR
      MODULE PROCEDURE PRINTVARD, PRINTVARI
      END INTERFACE 

      INTERFACE PRINTVEC
      MODULE PROCEDURE PRINTVECD, PRINTVECI
      END INTERFACE 


      CONTAINS
      SUBROUTINE ECHO(array)
      INTEGER ::j,i
      DOUBLE PRECISION,DIMENSION(:,:)::array
      CHARACTER*80 :: string
      DO j=1,size(array,1)         
         PRINT 111,j,array(j,:)
111      FORMAT (i2,':',10F10.5)
      END DO
C      DO j=1,size(array,1)
C         WRITE(string,110) j
C 110     FORMAT (i2,':')
C         CALL mexprintf(string)
C         DO i = 1, size(array,2)
C            WRITE(string,111) array(j,i)
C 111        FORMAT (' ',10F10.5)
C            CALL mexprintf(string)
C         ENDDO
C         CALL mexprintf(CHAR(10))   !CHAR(10) is a <CR> 
C      END DO     
      END SUBROUTINE ECHO

      SUBROUTINE PRINTVARI(I,TXT)
      INTEGER,  INTENT(IN) :: I
      CHARACTER*80, OPTIONAL, INTENT(IN) :: TXT
      CHARACTER*80 :: string1
      IF (PRESENT(TXT)) THEN
	 WRITE(*,*) TXT
       WRITE(*,*)  ':'
      ENDIF
      WRITE(string1,125) I
 125  FORMAT (' ',i5) 
      WRITE(*,*) string1
	WRITE(*,*)  CHAR(10) 
      RETURN
      END SUBROUTINE PRINTVARI
      SUBROUTINE PRINTVARD(D,TXT)
      DOUBLE PRECISION,  INTENT(IN) :: D
      CHARACTER*80, OPTIONAL, INTENT(IN) :: TXT
      CHARACTER*80 :: string1
      IF (PRESENT(TXT)) THEN
        WRITE(*,*) TXT
        WRITE(*,*)  ':'
      ENDIF
      WRITE(string1,115) D
 115  FORMAT (' ',10F10.5) 
      WRITE(*,*) string1
	WRITE(*,*)  CHAR(10) 
      RETURN
      END SUBROUTINE PRINTVARD

      SUBROUTINE PRINTVECD(CDI,TXT,N)
      DOUBLE PRECISION, DIMENSION(:),INTENT(in) :: CDI
	INTEGER, OPTIONAL, INTENT(IN):: N
      CHARACTER*80, OPTIONAL, INTENT(IN) :: TXT
      INTEGER :: I,N1
      CHARACTER*80 :: string
      IF (PRESENT(TXT)) THEN
        WRITE(*,*) TXT
        WRITE(*,*)  ':'
      ENDIF
	IF (PRESENT(N)) THEN
	 N1 = N
	ELSE
	 N1 = SIZE(CDI,1)
	ENDIF
      DO I = 1, N1
         WRITE(*,115) CDI(I)
 115     FORMAT (' ',10F10.5) 
!         CALL mexprintf(string) 
      ENDDO
      WRITE(*,*) CHAR(10) 
      RETURN
      END SUBROUTINE PRINTVECD
      SUBROUTINE PRINTVECI(CDI,TXT,N)
      INTEGER, DIMENSION(:),INTENT(in) :: CDI
	INTEGER, OPTIONAL, INTENT(IN):: N
      CHARACTER*80, OPTIONAL, INTENT(IN) :: TXT
      INTEGER :: I, N1
      CHARACTER*80 :: string
      IF (PRESENT(TXT)) THEN
        WRITE(*,*) TXT
        WRITE(*,*)  ':'
      ENDIF
	IF (PRESENT(N)) THEN
	 N1 = N
	ELSE
	 N1 = SIZE(CDI,1)
	ENDIF
      DO I = 1, N1
         WRITE(*,115) CDI(I)
 115     FORMAT (' ',i5) 
!         CALL mexprintf(string) 
      ENDDO
      WRITE(*,*) CHAR(10) 
      RETURN
      END SUBROUTINE PRINTVECI

      SUBROUTINE PRINTCOF(Ntd,A,B,INFI,COF,INDEX1)
      IMPLICIT NONE
      INTEGER ,INTENT(in) :: Ntd
      DOUBLE PRECISION, DIMENSION(:),INTENT(in) :: A,B
      DOUBLE PRECISION, DIMENSION(:,:),INTENT(in) :: COF
      INTEGER , DIMENSION(:), INTENT(in) :: INFI,INDEX1
! Local variables
      INTEGER :: I, J, K
      CHARACTER*80 :: string
      CALL mexprintf('  Lower  Upper  Cholesky Matrix '//CHAR(10))
      
      J = Ntd !MIN(Ntd,5)
      DO I = 1, Ntd   
         IF ( INFI(I)  <  0 ) THEN 
            WRITE(string,111) INDEX1(I)
 111        FORMAT (i2,'  -inf    inf '\) 
!            CALL mexprintf(string) 
         ELSE IF ( INFI(I) .EQ. 0 ) THEN 
            WRITE(string,112) INDEX1(I),B(I)
 112        FORMAT (i2,'  -inf ', 5F5.2)
!            CALL mexprintf(string) 
         ELSE IF ( INFI(I) .EQ. 1 ) THEN 
            WRITE(string,113) INDEX1(I),A(I)
 113        FORMAT (i2,' ', 5F5.2 ,'   inf '\)
!            CALL mexprintf(string) 
!            CALL mexprintf('   inf ')
         ELSE 
            WRITE(string,114) INDEX1(I),A(I),B(I)
 114        FORMAT (i2,' ', 5F5.2 ,' ', 5F5.2)
!           CALL mexprintf(string) 
         END IF
         CALL mexprintf(string)
         DO K = 1,J
            WRITE(string,115) COF(I,K)
 115        FORMAT (' ',10F10.5) 
            CALL mexprintf(string) 
         ENDDO
         CALL mexprintf(CHAR(10)) 
      END DO
      RETURN
      END SUBROUTINE PRINTCOF
      end module PRINTMOD77
