! Calculated RDF using force fields 

SUBROUTINE calc_frdf(sim_type, ntimesteps, nions, natom, atype, acell, nCr, forces, xyz, gr)

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

    INTEGER :: i,j,k,l,m,n, r, nbins, at1, at2
    INTEGER, INTENT(IN) :: ntimesteps, nions, nCr
    INTEGER, DIMENSION(:), INTENT(IN) :: natom(nions)
    INTEGER, DIMENSION(:) :: ntype(nions+1)
    REAL*8, DIMENSION(:) :: counts(2), fv(3)
    REAL*8, DIMENSION(:,:) :: gij(3,3), v(3,3)
    REAL*8, DIMENSION(:,:,:):: gr_avg(nCr,125,2)
    REAL*8, DIMENSION(:,:,:), INTENT(IN) :: acell(ntimesteps,3,3), xyz(ntimesteps,SUM(natom),3), forces(ntimesteps, SUM(natom),6)
    REAL*8, DIMENSION(:,:,:,:), INTENT(OUT) :: gr(nCr,ntimesteps,250,2)
    REAL*8 :: dummy, beta, delta, rmax, vol, rvec, rdisp, fdisp, inner_dot1, atom_sum1, density
    REAL*8 :: temperature, inner_dot2, atom_sum2, delta_x, var1, cov1, var_avg, cov_avg, lambda_star, est_star
    CHARACTER(len=3), DIMENSION(:) :: aname(nions+1)
    CHARACTER(len=3), DIMENSION(:), INTENT(IN) :: atype(nions)
    REAL*8, PARAMETER :: pi = 4.0*ATAN(1.0), ev_kb = 8.617333262E-5, ev_J  = 1.380649E-23
    LOGICAL, INTENT(IN) :: sim_type

    WRITE(*,*) "CALCULATING FORCE RADIAL DISTRIBUTION FUNCTION"
    WRITE(*,*)

    CALL EXECUTE_COMMAND_LINE( "grep -m 1 'TEBEG  =' OUTCAR | cut -c 1-18 | awk '{print $3}' > temp ") 
    OPEN(30,FILE='temp',STATUS='OLD')
    READ(30,*) dummy
    temperature = dummy
    CLOSE(30,STATUS='DELETE')

    ! Initialize
    beta = 1/(ev_kb * temperature)
    vol = leibniz(acell(1,1,:), acell(1,2,:), acell(1,3,:))
    rmax = vol**(1./3.)
    gr(:,:,:,:) = 0

    WRITE(*,*) "Temperature = ", temperature
    WRITE(*,*) "beta = ", beta

    ntype(:) = 0

    ! Assign number of atoms and atom names to arrays: ntype and aname
    DO i = 1,nions
        ntype(i+1) = ntype(i+1) + ntype(i) + natom(i)
        aname(i+1) = atype(i)
        PRINT*, ntype(i+1), aname(i+1)
    END DO
    

    OPEN(UNIT=5,FILE='atoms',STATUS='UNKNOWN')
    ! print out atoms pairs for python plotting legend
    DO m=2,nions+1
        DO l = m,nions+1
            WRITE(5,*) aname(m), aname(l)
        END DO
    END DO
    CLOSE(5)

  ! Determine whether XDATCAR is an NPT/NPH or NVT simulation
  IF (sim_type) THEN 
    PRINT*, sim_type
    PRINT*, '************************************** Calculating NPT or NPH displacements **************************************'
    
    ! Initialize RDF parameters
    delta = 0.05
    WRITE(*,*) 'Delta: ', delta
    nbins = rmax/delta
    
    ! Calculate distance vectors between two points. Then take the dot product of vectors with metric tensor Gij
    n = 0
    DO m = 2,nions+1
      DO l = m,nions+1
        n = n + 1
        WRITE(*,*)
        WRITE(*,*) 'Computing atom ', aname(m), aname(l), ' pair'
        WRITE(*,*)
        at1 = ntype(m-1)+1
        at2 = ntype(l-1)+1

        density = natom(m-1)*natom(l-1)
        ! WRITE(*,*) ab, ntype(m)-1, cd, ntype(l)
        DO i = 1,ntimesteps
          gij(:,:) = 0.0
          CALL metric(acell(i,:,:),gij(:,:),1.0)
          DO j = at1,ntype(m)-1
            DO k = at2,ntype(l)
            
              ! v(1,1) = (xyz(i,j,1) - xyz(i,k,1)) - ANINT((xyz(i,j,1) - xyz(i,k,1)))
              ! v(1,2) = (xyz(i,j,2) - xyz(i,k,2)) - ANINT((xyz(i,j,2) - xyz(i,k,2)))
              ! v(1,3) = (xyz(i,j,3) - xyz(i,k,3)) - ANINT((xyz(i,j,3) - xyz(i,k,3)))
              ! rdisp = SQRT(gdot(v1(:), v1(:), gij(:,:)))

              ! fv(1,1) = (forces(i,k,4) - forces(i,j,4))
              ! fv(1,2) = (forces(i,k,5) - forces(i,j,5))
              ! fv(1,3) = (forces(i,k,6) - forces(i,j,6))
              ! fdisp = SQRT( fv(1)**2. + fv(2)**2. + fv(3)**2. )
            END DO
          END DO
          IF (mod(i,100) .eq. 0) write(6,'(a25,i7,a10)',ADVANCE='NO') 'Computing displacement ',i,'         '//CHAR(13)
        END DO
      END DO
    END DO

  ELSE
    PRINT*, sim_type
    PRINT*, '************************************** Calculating NVT displacements **************************************'

    ! Calculate metric tensor Gij
    CALL metric(acell(1,:,:),gij(:,:),1.0)
    ! print*, (gij(:,:))

    ! Initialize number of bins
    nbins = 250
    delta = rmax/nbins

    ! Calculate distance vectors between two points. Then take the dot product of vectors with metric tensor Gij
    n = 0
    DO m = 2,nions+1
      DO l = m,nions+1
        n = n + 1
        WRITE(*,*)
        WRITE(*,*) 'Computing atom ', aname(m), aname(l), ' pair'
        WRITE(*,*)
        at1 = ntype(m-1)+1
        at2 = ntype(l-1)+1

        density = natom(m-1)*natom(l-1)

        DO i = 1,ntimesteps
          ! atom_sum = 0
          DO r = 1, nbins
            atom_sum1 = 0
            atom_sum2 = 0
            rvec = (r * delta) - (delta/2)
            IF (rvec.gt.rmax/2) EXIT
            DO j = at1,ntype(m)
              ! inner_sum = 0
              DO k = at2,ntype(l)
                IF (j.eq.k) CYCLE
                v(1,1) = (xyz(i,j,1) - xyz(i,k,1)) - ANINT((xyz(i,j,1) - xyz(i,k,1)))
                v(1,2) = (xyz(i,j,2) - xyz(i,k,2)) - ANINT((xyz(i,j,2) - xyz(i,k,2)))
                v(1,3) = (xyz(i,j,3) - xyz(i,k,3)) - ANINT((xyz(i,j,3) - xyz(i,k,3)))
                rdisp = SQRT(gdot(v(1,:), v(1,:), gij(:,:)))
                ! rdisp = SQRT(v1(1)**2. + v1(2)**2. + v1(3)**2.)
                ! PRINT*, rdisp

                ! Heaviside function
                IF ((rdisp - rvec).gt.0) THEN
                  counts(1) = 1
                ELSEIF ((rdisp - rvec).le.0) THEN
                  counts(1) = 0
                END IF

                ! Heaviside function
                IF ((rvec - rdisp).gt.0) THEN
                  counts(2) = 1
                ELSEIF ((rvec - rdisp).le.0) THEN
                  counts(2) = 0
                END IF

                ! lambda = 0 g(r) RHS
                v(2,1) = (v(1,1)*acell(1,1,1)/(rdisp**3.))*counts(1)
                v(2,2) = (v(1,2)*acell(1,2,2)/(rdisp**3.))*counts(1)
                v(2,3) = (v(1,3)*acell(1,3,3)/(rdisp**3.))*counts(1)

                ! lamba = 1 g(r) RHS
                v(3,1) = (v(1,1)*acell(1,1,1)/(rdisp**3.))*counts(2)*(-1)
                v(3,2) = (v(1,2)*acell(1,2,2)/(rdisp**3.))*counts(2)*(-1)
                v(3,3) = (v(1,3)*acell(1,3,3)/(rdisp**3.))*counts(2)*(-1)

                fv(1) = (forces(i,k,4) - forces(i,j,4))
                fv(2) = (forces(i,k,5) - forces(i,j,5))
                fv(3) = (forces(i,k,6) - forces(i,j,6))
                ! fdisp = SQRT( fv(1)**2. + fv(2)**2. + fv(3)**2. )
                fv = fv/2

                inner_dot1 = DOT_PRODUCT(fv(:), v(2,:))
                inner_dot2 = DOT_PRODUCT(fv(:), v(3,:))

                atom_sum1 = atom_sum1 + inner_dot1
                atom_sum2 = atom_sum2 + inner_dot2
              END DO
              ! atom_sum = inner_sum
            END DO
            gr(n,i,r,1) = 1 + ((vol/density) * (beta/(4*pi)) * atom_sum1)
            gr(n,i,r,2) = ((vol/density) * (beta/(4*pi)) * atom_sum2)
          END DO
          IF (mod(i,100) .eq. 0) write(6,'(a25,i7,a10)',ADVANCE='NO') 'Computing displacement ',i,'         '//CHAR(13)
        END DO
      END DO
    END DO

    130  FORMAT(f18.10, 8X, f18.10)

    k = 21
    DO n = 1, nCr
      WRITE(k,130)
      k = k + 1
      DO r = 1,nbins
        rvec = (r * delta) - (delta/2)
        IF (rvec.gt.rmax/2) EXIT
        gr_avg(n,r,1) = SUM(gr(n,:,r,1),DIM=1)/ntimesteps ! Time averaged g(r) for estimator 1
        IF (n.eq.1) THEN 
          WRITE(k,130) rvec, SUM(gr(n,:,r,1),DIM=1)/ntimesteps
        ELSE
          WRITE(k,130) SUM(gr(n,:,r,1),DIM=1)/ntimesteps
        END IF
      END DO
    END DO

    k = 31
    DO n = 1, nCr
      WRITE(k,130)
      k = k + 1
      DO r = 1,nbins
        rvec = (r * delta) - (delta/2)
        IF (rvec.gt.rmax/2) EXIT
        gr_avg(n,r,2) = SUM(gr(n,:,r,2),DIM=1)/ntimesteps ! Time averaged g(r) for estimator 2
        IF (n.eq.1) THEN 
          WRITE(k,130) rvec, SUM(gr(n,:,r,2),DIM=1)/ntimesteps
        ELSE
          WRITE(k,130) SUM(gr(n,:,r,2),DIM=1)/ntimesteps
        END IF
      END DO
    END DO

  CALL EXECUTE_COMMAND_LINE("rm fort.21")
  CALL EXECUTE_COMMAND_LINE("paste fort.2* > rad1")
  CALL EXECUTE_COMMAND_LINE("rm fort.2*")

  CALL EXECUTE_COMMAND_LINE("rm fort.31")
  CALL EXECUTE_COMMAND_LINE("paste fort.3* > rad2")
  CALL EXECUTE_COMMAND_LINE("rm fort.3*")

    ! Calculate variances and co-variances 
    k = 51
    l = 61
    DO n = 1,nCr
      WRITE(k,130)
      k = k + 1
      l = l + 1
      DO r = 1,nbins
        var1 = 0
        cov1 = 0
        rvec = (r * delta) - (delta/2)
        IF (rvec.gt.rmax/2) EXIT
        DO i = 1,ntimesteps
          delta_x = gr(n,i,r,2) - gr(n,i,r,1)
          var1 = var1 + ((gr(n,i,r,1) - gr_avg(n,r,1))**2.)
          cov1 = cov1 + ((gr(n,i,r,1) - gr_avg(n,r,1)) * (delta_x - (gr_avg(n,r,2)-gr_avg(n,r,1))))
        END DO
        var_avg = var1/ntimesteps
        cov_avg = cov1/ntimesteps
        lambda_star = (-1)*cov_avg/var_avg

        est_star = ((1-lambda_star)*gr_avg(n,r,1)) + (lambda_star*gr_avg(n,r,2))
        IF (n.eq.1) THEN 
          WRITE(k,130) rvec, est_star
          WRITE(l,*) rvec, lambda_star
        ELSE
          WRITE(k,130) est_star
          WRITE(l,*) lambda_star
        END IF

      END DO
    END DO

  CALL EXECUTE_COMMAND_LINE("rm fort.51")
  CALL EXECUTE_COMMAND_LINE("paste fort.5* > rad3")
  CALL EXECUTE_COMMAND_LINE("rm fort.5*")

  CALL EXECUTE_COMMAND_LINE("paste fort.6* > lambda_star_gr")
  CALL EXECUTE_COMMAND_LINE("rm fort.6*")

  END IF

END SUBROUTINE calc_frdf