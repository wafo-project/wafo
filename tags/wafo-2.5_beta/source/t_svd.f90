PROGRAM Test_SVD

! A simple test of the SVD module.

! The matrix:

!  ( 1  2  3 )
!  ( 4  5  6 )
!  ( 7  8  9 )
!  (10 11 12 )

! has singular values: 25.4624, 1.2907 and zero.

USE SVD
IMPLICIT NONE

REAL (dp)  :: e(3), x(4,3), s(4), u(4,4), v(3,3), value
INTEGER    :: info, row, col

value = 1.0_dp
DO row = 1, 4
  DO col = 1, 3
    x(row,col) = value
    value = value + 1.0_dp
  END DO
END DO

CALL dsvdc(x, 4, 3, s, e, u, v, 11, info)

! Output the singular values

WRITE(*, 900) s(1:4)
900 FORMAT(' The calculated singular values = ', 4f10.4/)

! Output the U matrix

WRITE(*, *) '     The U-matrix'
DO row = 1, 4
  WRITE(*, 910) u(row,1:4)
  910 FORMAT(4f10.4)
END DO

! Now the right-hand or V-matrix

WRITE(*, *)
WRITE(*, *) '     The V-matrix'
DO row = 1, 3
  WRITE(*, 910) v(row,1:3)
END DO

STOP
END PROGRAM Test_SVD
