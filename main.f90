PROGRAM main

  ! Test out distance calculations from metric tensor and gdot subroutine
  USE element_table_mod, only: get_element, atomic_element

  IMPLICIT NONE

  INTERFACE
    FUNCTION factorial(n)
      INTEGER :: factorial 
      INTEGER, INTENT(IN) :: n
    END FUNCTION factorial

    REAL FUNCTION gdot(a,b,gij)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION a(3), b(3), gij(3,3)
    END FUNCTION gdot
  END INTERFACE

  INTEGER :: i,j,k,l,m,n
  INTEGER :: ntimesteps, nions, totatoms, xmax, ymax, zmax, nCr
  INTEGER, PARAMETER :: r=2
  INTEGER, ALLOCATABLE, DIMENSION(:) :: natom
  REAL*8 :: start, finish
  REAL*8, DIMENSION(:) :: v1(3), v2(3)
  REAL*8, ALLOCATABLE, DIMENSION(:) :: volume, time, wmass
  REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: acell, xyz, forces, deriv, force_density
  REAL*8, ALLOCATABLE, DIMENSION(:,:,:,:) :: poly, gr, force_density_time
  CHARACTER(len=3), ALLOCATABLE, DIMENSION(:) :: atype
  TYPE(atomic_element) :: element
  REAL*8, PARAMETER :: avog = 6.022E23
  LOGICAL :: sim_type, file_exists 

  ! Time how long program runs
  CALL CPU_TIME(start)

  ! Find parameters: number of timesteps, number of ions in system, and total atoms
  CALL find_params(ntimesteps, nions, totatoms)

  ! Allocate arrays according to parameters
  ALLOCATE(natom(nions))
  ALLOCATE(atype(nions))
  ALLOCATE(volume(ntimesteps))
  ALLOCATE(acell(ntimesteps,3,3))
  ALLOCATE(forces(ntimesteps,totatoms,6))
  ALLOCATE(xyz(ntimesteps,totatoms,3))

  ! Read in XDATCAR
  CALL read_xdat(ntimesteps,nions,totatoms,acell,volume,natom,xyz,atype,sim_type)

  OPEN(UNIT=505,FILE='atoms_rho',STATUS='UNKNOWN')

  ! Write out atom names for density figures
  DO n = 1,nions
    WRITE(505,*) atype(n)
  END DO

  CLOSE(505)

  ! True == NPT or NPH simulation
  ! False == NVT simulation
  PRINT*, sim_type

  ! Initialize
  xmax = INT(CEILING(acell(1,1,1)))
  ymax = INT(CEILING(acell(1,2,2)))
  zmax = INT(CEILING(acell(1,3,3)))

  ! DO i = 1, ntimesteps
  !   DO j = 1, SUM(natom)
  !       xyz(i,j,3) = xyz(i,j,3) + 0.25
  !       IF (xyz(i,j,3).gt.1) xyz(i,j,3) = xyz(i,j,3) - 1
  !       IF (xyz(i,j,3).lt.0) xyz(i,j,3) = xyz(i,j,3) + 1
  !   END DO
  ! END DO

  ! Allocate arrays for force calculation
  ALLOCATE(time(ntimesteps))
  ALLOCATE(wmass(nions))
  ALLOCATE(deriv(ntimesteps, SUM(natom),9))
  ALLOCATE(poly(ntimesteps,SUM(natom),3,3))

  !Get atomic masses of the system
  DO i=1,nions
    ! WRITE(*,*) atype(i)
    element = get_element(atype(i))
    wmass(i) = (element%mass)/avog/1000
    ! WRITE(*,*) "Atomic Mass: ", atype(i), " = ", wmass(i), "kg"
    ! WRITE(*,*) "Atomic Mass: ", atype(i), " = ", wmass(i)*avog, "amu"
  END DO 

  ! Make time array into a real number
  DO i = 1,ntimesteps
    time(i) = i
  END DO

  ! Find total number of atom-pair combinations
  m = nions 
  nCr = (factorial(m) / (factorial(r) * factorial(m-r))) + nions 
  WRITE(*,*) 'Total number of atom type combinations: ', nCr

  !CALL calc_forces(ntimesteps, nions, time, natom, atype, wmass, acell, xyz, poly, deriv)

  CALL find_forces(ntimesteps, totatoms, forces)

  ! Allocate Force density arrays
  ALLOCATE(force_density_time(nions,ntimesteps,100,3))
  ALLOCATE(force_density(nions,100,6))
  ALLOCATE(gr(nCr,ntimesteps,250,2))

  CALL calc_1dforcedens(ntimesteps, nions, natom, acell, 100, forces, xyz, force_density_time, force_density)
  ! CALL calc_frdf(sim_type, ntimesteps, nions, natom, atype, acell, nCr, forces, xyz, gr)

  ! Deallocate arrays
  DEALLOCATE(natom, atype, volume, acell, forces, xyz, force_density_time, force_density, gr)
  DEALLOCATE(time, wmass, deriv, poly)

  CALL CPU_TIME(finish)
  PRINT*, "Time = ", finish-start, "seconds ( = ", (finish-start)/60, "minutes.)"

END PROGRAM main