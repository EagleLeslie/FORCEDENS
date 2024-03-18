! subroutine reading in VASP XDATCAR
! Requires linux/bash command line for "grep, awk, and sed" commands
!
! If user does not have linux/bash shell, uncomment commands that prompt
! user for the number of time steps and number of atom types.
! Don't forget to comment out areas that call out "CALL EXECUTE_COMMAND_LINE"

SUBROUTINE read_xdat(ntimesteps,nions,totatoms,acell,volume,natom,xyz,atype,sim_type)

  IMPLICIT NONE

  INTERFACE
    REAL FUNCTION leibniz(a1, a2, a3)
      IMPLICIT NONE
      REAL*8, DIMENSION(3), INTENT(IN) :: a1, a2, a3
    END FUNCTION leibniz
  END INTERFACE

  INTEGER :: i,j,k,l
  INTEGER, INTENT(IN) :: ntimesteps, nions, totatoms
  INTEGER, DIMENSION(:), INTENT(OUT) :: natom(nions)
  REAL*8, DIMENSION(:), INTENT(OUT) :: volume(ntimesteps)
  REAL*8, DIMENSION(:,:,:), INTENT(OUT) :: acell(ntimesteps,3,3), xyz(ntimesteps,totatoms,3)
  REAL*8 :: scale, dx, dy, dz, disp
  CHARACTER(len=50) :: header, isif
  CHARACTER(len=3), DIMENSION(:), INTENT(OUT) :: atype(nions)
  LOGICAL, INTENT(OUT) :: sim_type

  OPEN(UNIT=10,FILE='XDATCAR',STATUS='OLD',ERR=200)

  READ(10,*) header
  PRINT*, 'Header = ', header

  READ(10,*) scale
  READ(10,*) acell(1,1,:)
  READ(10,*) acell(1,2,:)
  READ(10,*) acell(1,3,:)
  READ(10,*) atype(:)
  READ(10,*) natom(:)
  READ(10,*) ! "Direct configuration="

  DO i = 1,3
    acell(1,i,1) = acell(1,i,1)*scale
    acell(1,i,2) = acell(1,i,2)*scale
    acell(1,i,3) = acell(1,i,3)*scale
  END DO

  PRINT*, scale
  PRINT*, acell(1,1,:)
  PRINT*, acell(1,2,:)
  PRINT*, acell(1,3,:)
  PRINT*, atype(:)
  PRINT*, natom(:)

  volume(1) = leibniz(acell(1,1,:), acell(1,2,:), acell(1,3,:))
  PRINT*, "Initial volume: ", volume(1)

  DO j = 1,totatoms
    READ(10,*) xyz(1,j,:)
  END DO

  READ(10,*) isif ! check for NPH or NPT ISIF is toggled in XDATCAR run
  IF (isif == header) THEN
    PRINT*, 'SIMULATION TYPE: NPH or NPT SIMULATION'
    sim_type = .TRUE.
    DO i = 2,ntimesteps
      READ(10,*) ! scale
      READ(10,*) acell(i,1,:)
      READ(10,*) acell(i,2,:)
      READ(10,*) acell(i,3,:)
      READ(10,*) ! atom types
      READ(10,*) ! number of atom types
      READ(10,*) ! Direct configuration=
      DO k = 1,3
          acell(i,k,1) = acell(i,k,1)*scale
          acell(i,k,2) = acell(i,k,2)*scale
          acell(i,k,3) = acell(i,k,3)*scale
      END DO
      volume(i) = leibniz(acell(i,1,:), acell(i,2,:), acell(i,3,:))
      ! PRINT*, volume(i)
      DO j = 1,totatoms
          READ(10,*) xyz(i,j,:)
      END DO
      READ(10,*,END=301)
      END DO
  ELSE
    PRINT*, 'SIMULATION TYPE: NVT SIMULATION'
    sim_type = .FALSE.
    DO i = 2,ntimesteps
      DO j = 1,totatoms
        READ(10,*) xyz(i,j,:)
      END DO
      READ(10,*,END=300)
    END DO
  END IF
  300 CONTINUE
  301 CONTINUE

  200 CONTINUE
  CLOSE(10)

END SUBROUTINE read_xdat
