SUBROUTINE find_forces(ntimesteps, totatoms, forces)

    IMPLICIT NONE

    INTEGER :: i,j,k,l
    INTEGER, INTENT(IN) :: ntimesteps, totatoms
    REAL*8, DIMENSION(:,:,:), INTENT(OUT) :: forces(ntimesteps,totatoms,6)

    PRINT*, "READING IN ATOMIC FORCES"

    ! Get forces from OUTCAR
    ! FORMAT: XX YY ZZ XY YZ ZX
    ! AX where X = total number of atoms + 2
    ! CALL EXECUTE_COMMAND_LINE("grep -A119 'TOTAL-FORCE (eV/Angst)' OUTCAR > forces") ! For FeMgO with 117 atoms
    ! CALL EXECUTE_COMMAND_LINE("grep -A150 'TOTAL-FORCE (eV/Angst)' OUTCAR > forces") ! For FeMgSiO3 with 148 atoms
    ! CALL EXECUTE_COMMAND_LINE("grep -A92 'TOTAL-FORCE (eV/Angst)' OUTCAR > forces")
    CALL EXECUTE_COMMAND_LINE("python grab_forces.py")
    OPEN(11,FILE='forces')

    READ(11,*)
    READ(11,*) 

    DO i = 1,ntimesteps
        DO j = 1,totatoms
            READ(11,*,END=100) forces(i,j,:)
        END DO
        DO k = 1,4
            READ(11,*,END=111)
        END DO
    END DO
    100 CONTINUE
    111 CONTINUE
    CLOSE(11,STATUS='DELETE')

END SUBROUTINE find_forces