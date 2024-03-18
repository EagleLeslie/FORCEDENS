SUBROUTINE calc_1dforcedens(ntimesteps, nions, natom, acell,  zmax, forces, xyz, force_density_time, force_density)

  IMPLICIT NONE

  INTERFACE
    REAL FUNCTION gdot(a,b,gij)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION a(3), b(3), gij(3,3)
    END FUNCTION gdot

    REAL FUNCTION leibniz(a1, a2, a3)
      IMPLICIT NONE
      REAL*8, DIMENSION(3), INTENT(IN) :: a1, a2, a3
    END FUNCTION leibniz
  END INTERFACE

  INTEGER :: i,j,k,l,n,x,y,z,t0,tf,tdiff,nbins,at1,at2
  INTEGER, INTENT(IN) :: ntimesteps, nions, zmax
  INTEGER, DIMENSION(:), INTENT(IN) :: natom(nions)
  REAL*8, DIMENSION(:) :: v1(3), counts(3), force_sum(3)
  REAL*8, DIMENSION(:,:) :: gij(3,3)
  ! REAL*8, DIMENSION(:,:,:) :: bins(15237,100,6)
  REAL*8, DIMENSION(:,:,:), INTENT(IN) :: acell(ntimesteps,3,3), xyz(ntimesteps,SUM(natom),3), forces(ntimesteps, SUM(natom),6)
  REAL*8, DIMENSION(:,:,:), INTENT(OUT) :: force_density(nions,zmax,6)
  REAL*8, DIMENSION(:,:,:,:), INTENT(OUT) :: force_density_time(nions,ntimesteps,zmax,3)
  REAL*8 :: dummy, poop, temperature, disp, xref, yref, zref, r_box, rmax, dr, dz
  REAL*8 :: rdispx, rdispy, rdispz, vol, rho0, delta_x, var1, cov1, var_avg, cov_avg, lambda_star, est_star
  REAL*8 :: tol ! angstroms
  REAL*8, PARAMETER :: pi = 4.0*ATAN(1.0), ev_kb = 8.617333262E-5, ev_J  = 1.380649E-23

  CALL EXECUTE_COMMAND_LINE( "grep -m 1 'TEBEG  =' OUTCAR | cut -c 1-18 | awk '{print $3}' > temp ") 
  OPEN(30,FILE='temp',STATUS='OLD')
  READ(30,*) dummy
  temperature = dummy
  CLOSE(30,STATUS='DELETE')

  WRITE(*,*) "Temperature = ", temperature
  vol = leibniz(acell(1,1,:), acell(1,2,:), acell(1,3,:))
  rho0 = natom(1)/vol/2
  WRITE(*,*) "rho0 = ", rho0

  PRINT*, "STARTING FORCE DENSITY CALCULATION"

  ! CALL metric(acell(1,:,:),gij(:,:),1.0)

  t0 = ntimesteps - ntimesteps/2
  tf = ntimesteps
  tdiff = tf - t0 + 1

  PRINT*, t0, tf, tdiff

  nbins = 100
  dz = acell(1,3,3)/nbins

  at1 = 0
  at2 = 0

  ! 1D force density calculation with Heaviside function
  DO n = 1,nions
    IF (n.eq.1) THEN
      at1 = 1
      at2 = natom(n)
    ELSE 
      at1 = at2+1
      at2 = at2 + natom(n)
    END IF
    PRINT*, at1, at2
    DO i = t0, tf
        ! counts(:) = 0
        DO z = 1,nbins
          ! zref = REAL(z)/zmax
          zref = ((dz * z) - (dz/2))/acell(1,3,3)
          r_box = SQRT(1.0)
          ! counts(:) = 0
          force_sum(:) = 0

          DO j = at1,at2
            v1(3) = (xyz(i,j,3) - zref) !- ANINT((xyz(i,j,3) - zref))
            v1(2) = (zref - xyz(i,j,3)) !- ANINT((zref - xyz(i,j,3))) 

            ! Heaviside function
            IF (v1(3).gt.0) THEN
              counts(3) = 1
            ELSEIF (v1(3).le.0) THEN
              counts(3) = 0
            END IF

            ! Heaviside function
            IF (v1(2).gt.0) THEN
              counts(2) = 1
            ELSEIF (v1(2).le.0) THEN
              counts(2) = 0
            END IF

            force_sum(2) = force_sum(2) + (counts(2) * forces(i,j,6))
            force_sum(3) = force_sum(3) + (counts(3) * forces(i,j,6))
          END DO

          force_density_time(n,i,z,2) = force_sum(2)
          force_density_time(n,i,z,3) = force_sum(3) 
          ! force_density(1,1,z,3) = SUM(force_density_time(t0:tf,1,1,z,3),DIM=1)/tdiff /(ev_kb * temperature)
        END DO
    END DO
  END DO

  ! OPEN(UNIT=22,FILE='dens1',STATUS='UNKNOWN')

  130  FORMAT(f18.10, 8X, f18.10)

  k = 21
  l = 31
  DO n = 1,nions
    WRITE(k,130)
    WRITE(l,130)
    k = k + 1
    l = l + 1
    DO z = 1,nbins
      zref = (dz * z) - (dz/2)
      DO i = t0,tf
        ! force_density(1,1,z,2) = SUM(force_density_time(t0:tf,1,1,z,2),DIM=1)/tdiff/ &
          ! & (ev_kb * temperature)/acell(1,1,1)/acell(1,2,2) ! Time averaged rho for estimator 1
          
        ! force_density(1,1,z,3) = SUM(force_density_time(t0:tf,1,1,z,3),DIM=1)/tdiff/ &
          ! & (ev_kb * temperature)/acell(1,1,1)/acell(1,2,2)*(-1) ! Time averaged rho for estimator 2

          force_density(n,z,2) = force_density(n,z,2) + force_density_time(n,i,z,2)
          force_density(n,z,3) = force_density(n,z,3) + force_density_time(n,i,z,3)
        ! WRITE(22,*) zref, force_density(1,1,z,2), force_density(1,1,z,3)
      END DO
      force_density(n,z,2) = force_density(n,z,2)/(ev_kb * temperature)/acell(1,1,1)/acell(1,2,2)/(tdiff)
      force_density(n,z,3) = force_density(n,z,3)/(ev_kb * temperature)/acell(1,1,1)/acell(1,2,2)/(tdiff)*(-1)
      IF (n.eq.1) THEN 
        WRITE(k,130) zref, force_density(n,z,2)
        WRITE(l,130) zref, force_density(n,z,3)
      ELSE
        WRITE(k,130) force_density(n,z,2)
        WRITE(l,130) force_density(n,z,3)
      END IF
      ! WRITE(22,*) zref, force_density(n,z,2), force_density(n,z,3)
    END DO
  END DO

  CALL EXECUTE_COMMAND_LINE("rm fort.21")
  CALL EXECUTE_COMMAND_LINE("paste fort.2* > dens1")
  CALL EXECUTE_COMMAND_LINE("rm fort.2*")

  CALL EXECUTE_COMMAND_LINE("rm fort.31")
  CALL EXECUTE_COMMAND_LINE("paste fort.3* > dens2")
  CALL EXECUTE_COMMAND_LINE("rm fort.3*")

  ! CLOSE(22)

  ! OPEN(UNIT=64,FILE='est_star_rho',STATUS='UNKNOWN')
  ! OPEN(UNIT=65,FILE='lambda_star_rho',STATUS='UNKNOWN')

  k = 41
  l = 51
  ! Calculate variances and co-variances 
  DO n = 1,nions
    WRITE(k,130)
    WRITE(l,130)
    k = k + 1
    l = l + 1
    DO z = 1,nbins
      var1 = 0
      cov1 = 0
      zref = (dz * z) - (dz/2)
      DO i = t0,tf
        delta_x = force_density_time(n,i,z,3) - force_density_time(n,i,z,2)
        var1 = var1 + ((force_density_time(n,i,z,2) - force_density(n,z,2))**2.)
        cov1 = cov1 + ((force_density_time(n,i,z,2) - force_density(n,z,2)) * &
          & (delta_x - (force_density(n,z,3) - force_density(n,z,2))))
      END DO
      var_avg = var1/tdiff
      cov_avg = cov1/tdiff

      lambda_star = (-1)*cov_avg/var_avg

      est_star = ((1-lambda_star)*force_density(n,z,2)) + (force_density(n,z,3))

      IF (n.eq.1) THEN 
        WRITE(k,130) zref, est_star
        WRITE(l,130) zref, lambda_star
      ELSE
        WRITE(k,130) est_star
        WRITE(l,130) lambda_star
      END IF
      ! WRITE(64,*) n, zref, est_star
      ! WRITE(65,*) n, zref, lambda_star
    END DO
  END DO

  CALL EXECUTE_COMMAND_LINE("rm fort.41")
  CALL EXECUTE_COMMAND_LINE("paste fort.4* > est_star_rho")
  CALL EXECUTE_COMMAND_LINE("rm fort.4*")

  CALL EXECUTE_COMMAND_LINE("rm fort.51")
  CALL EXECUTE_COMMAND_LINE("paste fort.5* > lambda_star_rho")
  CALL EXECUTE_COMMAND_LINE("rm fort.5*")

  ! CLOSE(64)
  ! CLOSE(65)

END SUBROUTINE calc_1dforcedens