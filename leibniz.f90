REAL FUNCTION leibniz(a1, a2, a3)

! Function that calculates the determinant of a 3x3 matrix A
! Requires 3 vectors with 3 components (ex: vector a with xyz component)
!
! |  a  b  c  |      |  a11  a12  a13  |
! |  d  e  f  |  ==  |  a21  a22  a23  |
! |  g  h  i  |      |  a31  a32  a33  |
!
! det(A) = aei + bfg + cdh - ceg - bdi - afh
! det(a) = (a11*a22*a33) + (a12*a23*a31) + (a13*a21*a32) - &
! & (a13*a22*a31) - (a12*a21*a33) - (a11*a23*a32)

  IMPLICIT NONE
  REAL*8, DIMENSION(3), INTENT(IN) :: a1, a2, a3

  leibniz = (a1(1)*a2(2)*a3(3)) + (a1(2)*a2(3)*a3(1)) + (a1(3)*a2(1)*a3(2)) - &
   & (a1(3)*a2(2)*a3(1)) - (a1(2)*a2(1)*a3(3)) - (a1(1)*a2(3)*a3(2))

END FUNCTION leibniz
