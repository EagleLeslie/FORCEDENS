SUBROUTINE find_params(ntimesteps, nions, totatoms)

  IMPLICIT NONE

  INTEGER :: i, j, k 
  INTEGER, ALLOCATABLE, DIMENSION(:) :: atype
  INTEGER, INTENT(OUT) :: ntimesteps, nions, totatoms

  !! Get number of timesteps from XDATCAR
  !! Comment out this block if no bash shell command line and execute:
  ! WRITE(*,*) 'How many total time steps?'
  ! READ(*,*) ntimesteps

  CALL EXECUTE_COMMAND_LINE("grep 'Direct' XDATCAR | tail -1| awk '{print $3}' > ntime")
  OPEN(11,FILE='ntime')
  READ(11,*,END=100) i
  100 CONTINUE
  ntimesteps = i
  CLOSE(11,STATUS='DELETE')

  WRITE(*,*) 'Total timesteps: ', ntimesteps

  !! Get number of ions
  !! Comment out this block if no bash shell command line and execute:
  ! WRITE(*,*) 'How many total ion types?' ! Example: Fe Mg O == 3 ion types
  ! READ(*,*) nions
  CALL EXECUTE_COMMAND_LINE("sed '7q;d' XDATCAR > blah")
  CALL EXECUTE_COMMAND_LINE("awk '{print NF}' blah | sort -nu | tail -n 1 > nions")
  OPEN(12,FILE='nions')
  READ(12,*,END=101) j
  101 CONTINUE
  nions = j
  CLOSE(12,STATUS='DELETE')

  WRITE(*,*) 'Total atom types: ', nions

  ALLOCATE(atype(nions))

  OPEN(13,FILE='blah')
  READ(13,*,END=102) atype(:)
  102 CONTINUE 
  totatoms = SUM(atype)
  CLOSE(13,STATUS='DELETE')

  WRITE(*,*) 'Total atoms: ', totatoms

  DEALLOCATE(atype)

  RETURN
END SUBROUTINE find_params
