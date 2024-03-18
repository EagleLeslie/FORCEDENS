SUBROUTINE calc_forces(ntimesteps, nions, time, natom, atype, wmass, acell, xyz, poly, deriv)

    IMPLICIT NONE

    INTEGER :: i,j,k,l, n
    INTEGER, PARAMETER :: mint=3
    INTEGER, INTENT(IN) :: ntimesteps, nions
    INTEGER, DIMENSION(:), INTENT(IN) :: natom(nions)
    REAL*8, DIMENSION(:) :: cofx(3), cofy(3), cofz(3)
    REAL*8, DIMENSION(:), INTENT(IN) :: time(ntimesteps), wmass(nions)
    REAL*8, DIMENSION(:,:,:), INTENT(IN) :: xyz(ntimesteps, SUM(natom),3), acell(ntimesteps,3,3)
    REAL*8, DIMENSION(:,:,:), INTENT(OUT) :: deriv(ntimesteps, SUM(natom),9)
    REAL*8, DIMENSION(:,:,:,:), INTENT(OUT) :: poly(ntimesteps,SUM(natom),3,3)
    CHARACTER(len=3), DIMENSION(:), INTENT(IN) :: atype(nions)
    REAL*8, PARAMETER :: avog = 6.022E23

    DO i = 1,ntimesteps-1
        DO j = 1,SUM(natom(:))
            ! Polynomial interpolation and calculate polynomial interpolation coefficients
            ! polint(xa,ya,n,x,y,dy)
            ! REAL :: dy,x,y,xa(n),ya(n)
            ! INTEGER :: n, NMAX
            ! Given arrays xa and ya, each of length n, and given a value of x, polint
            ! returns a value y, and an error estiamte dy
            !CALL polint(time(i), xyz(i,j,1), mint, x, y, dy)
            CALL polcof(time(i), xyz(i,j,1), mint, cofx)

            !CALL polint(time(i), xyz(i,j,2), mint, x, y, dy)
            CALL polcof(time(i), xyz(i,j,2), mint, cofy)

            !CALL polint(time(i), xyz(i,j,3), mint, x, y, dy)
            CALL polcof(time(i), xyz(i,j,3), mint, cofz)

            poly(i,j,1,:) = cofx
            poly(i,j,2,:) = cofy
            poly(i,j,3,:) = cofz

            ! Calculate velocities; derivative of polynomial interpolation output
            deriv(i,j,1) = poly(i,j,1,2) + (2.*time(i)*poly(i,j,1,3))
            deriv(i,j,2) = poly(i,j,2,2) + (2.*time(i)*poly(i,j,2,3))
            deriv(i,j,3) = poly(i,j,3,2) + (2.*time(i)*poly(i,j,3,3))

            ! Calculate acceleration; second derivative of polynomial interpolation output
            deriv(i,j,4) = (2.*poly(i,j,1,3))
            deriv(i,j,5) = (2.*poly(i,j,2,3))
            deriv(i,j,6) = (2.*poly(i,j,3,3))

        END DO
    END DO

    WRITE(*,*) wmass(:), atype(:)

    l = 0
    n = 1
    DO k = 1,nions
        l = l + natom(k)
        WRITE(*,*) k, n, l, wmass(k)
        DO j = n,l
            DO i = 1,ntimesteps
                deriv(i,j,7) = deriv(i,j,4) * wmass(k)
                deriv(i,j,8) = deriv(i,j,5) * wmass(k)
                deriv(i,j,9) = deriv(i,j,6) * wmass(k)
            END DO
        END DO
        n = n + natom(k)
    END DO
    
    DO i = 1,ntimesteps-1
        DO j = 1,SUM(natom(:))
            WRITE(88,*) i,j, deriv(i,j,7)*1E10*acell(1,1,1)*6.241509E18, deriv(i,j,8)*1E10*acell(1,2,2)*6.241509E18, &
            & deriv(i,j,9)*1E10*acell(1,3,3)*6.241509E18
        END DO
    END DO

END SUBROUTINE calc_forces 