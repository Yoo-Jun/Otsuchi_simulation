#include "macro.h"

SUBROUTINE step
  USE misc
  USE parameters
  USE calendar
  USE geometry
  USE velocity
  USE subgrid
  USE tracers
  USE io
  USE state
  USE seaice
  USE flux
  USE gls
  USE npzd
  USE ecosystem
  USE particles
#ifdef PROFILE
  USE profile
#endif
  IMPLICIT NONE

  REAL(8) :: tmp(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: avol
  INTEGER :: i, j, k
  INTEGER :: kl, km0, km1
  INTEGER :: n

  REAL(8) :: qsum

  LOGICAL :: stat

  n_timestep = n_timestep + 1

  CALL dump

  CALL update_flux
  CALL update_viscosity

  IF (.NOT. rigid_lid) THEN
     ! restore ssh
     CALL checkin('RESTORE_SSH', tmp(:,:,1), stat)

     IF (stat) THEN
        CALL checkin('RESRATE_SSH', tmp(:,:,2), stat)
        CALL assert(stat, "restoring rate for SSH (RESRATE_SSH) is not specified.")

!$OMP PARALLEL DO
        DO j=1-slv, jsize+slv
        DO i=1-slv, isize+slv
           IF (tmp(i,j,1) == UNDEF) CYCLE
           ssh(i,j) = ssh(i,j) + imask2d(i,j)*(tmp(i,j,1)-ssh(i,j))*tmp(i,j,2)*dtime
        END DO
        END DO
     END IF

!conservation of tracers is not satisfied in z-star mode with surface water suplly
     IF (fix_meanssh) THEN
        avol = 0.0
        IF (vrank==0) THEN
           DO i=1, isize
           DO j=1, jsize
              avol = avol + ssh(i,j)*dsz(i,j)*imask2d(i,j)
           END DO
           END DO
        END IF

        CALL gsum(avol, all=.TRUE.)
        avol = avol / total_area

        DO j=1-slv, jsize+slv
        DO i=1-slv, isize+slv
           ssh(i,j) = ssh(i,j) - imask2d(i,j)*avol
        END DO
        END DO
     END IF

     IF (stat_prec) THEN
        DO j=1-slv, jsize+slv
        DO i=1-slv, isize+slv
           ssh(i,j) = ssh(i,j) + imask2d(i,j)*prec(i,j)*dtime
        END DO
        END DO
     END IF

     IF (stat_evap) THEN
        DO j=1-slv, jsize+slv
        DO i=1-slv, isize+slv
           ssh(i,j) = ssh(i,j) - imask2d(i,j)*evap(i,j)*dtime
        END DO
        END DO
     END IF

     CALL update_geometry

     IF (vrank==0) THEN
        DO n=1, n_tracer
           ! T_new = (T_old*DV_old + P*T_prec - E*T_evap)/DV_new, where DV_new = DV_old + (P - E)*dsz*dtime
           !       =  T_old + (T_prec - T_old)*P*dsz*dtime/DV_new -(T_evap - T_old)*E*dsz*dtime/DV_new
           IF (stat_prec_tracer(n)) THEN
              DO j=1-slv, jsize+slv
              DO i=1-slv, isize+slv
                 IF (use_landwater) THEN
                    IF (imask3d(i,j,ksize)==0 .OR. lwdried(i,j)) CYCLE
                 END IF
                 IF (prec_tracer(i,j,n) == UNDEF) CYCLE
                 tracer(i,j,ksize,n) = tracer(i,j,ksize,n) + (prec_tracer(i,j,n)-tracer(i,j,ksize,n))*prec(i,j)*dsz(i,j)*dtime/dvol(i,j,ksize)
              END DO
              END DO
           END IF

           IF (stat_evap_tracer(n)) THEN
              DO j=1-slv, jsize+slv
              DO i=1-slv, isize+slv
                 IF (use_landwater) THEN
                    IF (imask3d(i,j,ksize)==0 .OR. lwdried(i,j)) CYCLE
                 END IF
                 IF (evap_tracer(i,j,n) == UNDEF) CYCLE
                 tracer(i,j,ksize,n) = tracer(i,j,ksize,n) - (evap_tracer(i,j,n)-tracer(i,j,ksize,n))*evap(i,j)*dsz(i,j)*dtime/dvol(i,j,ksize)
              END DO
              END DO
           END IF
        END DO
     END IF
  END IF

!$OMP PARALLEL
!$OMP WORKSHARE
  u_old(:,:,:) = u(:,:,:)
  v_old(:,:,:) = v(:,:,:)
  w_old(:,:,:) = w(:,:,:)
!$OMP END WORKSHARE

!$OMP DO
  DO k=1, ksize
  DO j=1, jsize
  DO i=0, isize
     gx(i,j,k) = -0.5D0*gx_ab(i,j,k)
     gx_ab(i,j,k) = 0.0D0
  END DO
  END DO
  END DO

!$OMP DO
  DO k=1, ksize
  DO j=0, jsize
  DO i=1, isize
     gy(i,j,k) = -0.5D0*gy_ab(i,j,k)
     gy_ab(i,j,k) = 0.0D0
  END DO
  END DO
  END DO

!$OMP DO
  DO k=0, ksize
  DO j=1, jsize
  DO i=1, isize
     gz(i,j,k) = -0.5D0*gz_ab(i,j,k)
     gz_ab(i,j,k) = 0.0D0
  END DO
  END DO
  END DO

  IF (seaice_coupled .AND. vrank==0) THEN
!$OMP DO
     DO j=1, jsize
     DO i=0, isize
        gix(i,j) = -0.5D0*gix_ab(i,j)
        gix_ab(i,j) = 0.0D0
     END DO
     END DO

!$OMP DO
     DO j=0, jsize
     DO i=1, isize
        giy(i,j) = -0.5D0*giy_ab(i,j)
        giy_ab(i,j) = 0.0D0
     END DO
     END DO
  END IF
!$OMP END PARALLEL

PROFILE_BEGIN('vadv')
  CALL gterm('NLIN', ab=.TRUE.)
  CALL nonlinear(gx_ab, gy_ab, gz_ab)
  CALL gterm('NLIN', ab=.TRUE.)
PROFILE_END('vadv')

  CALL gterm('CORI', ab=.TRUE.)
  CALL coriolis(gx_ab, gy_ab, gz_ab)
  CALL gterm('CORI', ab=.TRUE.)

  CALL checkout('GX', gx_ab) ! x-component  orsum of nonlinear and Coriolis term required for perfect restart [m/s^2]
  CALL checkout('GY', gy_ab) ! y-component  orsum of nonlinear and Coriolis term required for perfect restart [m/s^2]
  CALL checkout('GZ', gz_ab) ! z-component  orsum of nonlinear and Coriolis term required for perfect restart [m/s^2]

!$OMP PARALLEL PRIVATE(i, j, k)
!$OMP DO
  DO k=1, ksize
  DO j=1, jsize
  DO i=0, isize
     gx(i,j,k) = gx(i,j,k) + 1.5D0*gx_ab(i,j,k)
  END DO
  END DO
  END DO

!$OMP DO
  DO k=1, ksize
  DO j=0, jsize
  DO i=1, isize
     gy(i,j,k) = gy(i,j,k) + 1.5D0*gy_ab(i,j,k)
  END DO
  END DO
  END DO

!$OMP DO
  DO k=0, ksize
  DO j=1, jsize
  DO i=1, isize
     gz(i,j,k) = gz(i,j,k) + 1.5D0*gz_ab(i,j,k)
  END DO
  END DO
  END DO
!$OMP END PARALLEL

  CALL gterm('BODY')
  CALL body_forcing(gx, gy, gz)
  CALL gterm('BODY')

  CALL gterm('WAVE')
  CALL wave_forcing(gx, gy, gz)
  CALL gterm('WAVE')

  CALL gterm('STRS')
  CALL surface_stress(gx, gy, gix, giy)
  CALL bottom_stress(gx, gy)
  CALL gterm('STRS')

  IF (use_gls) CALL gls_mixing(visc(:,:,:,3), diff(:,:,:,3))

PROFILE_BEGIN('visc')
  CALL gterm('VISC')
  CALL viscosity(gx, gy, gz)
  CALL gterm('VISC')

  CALL biharmonic_filter(gx, gy, gz)
PROFILE_END('visc')

  CALL grad_slp(gx, gy)

  CALL gterm('BUOY')
  CALL buoyancy(gx, gy, gz)
  CALL gterm('BUOY')

  IF (use_landwater .AND. vrank==0) THEN
     DO j=1, jsize
     DO i=0, isize
        IF (lwflag(i,j) .AND. lwflag(i+1,j)) gx(i,j,ksize) = 0.0D0
     END DO
     END DO

     DO j=0, jsize
     DO i=1, isize
        IF (lwflag(i,j) .AND. lwflag(i,j+1)) gy(i,j,ksize) = 0.0D0
     END DO
     END DO
  END IF

  IF (use_particles) CALL particle_driven_forcing(gx, gy, gz)

  CALL gterm('REST')
  CALL restore_velocity(gx, gy, gz)
  CALL gterm('REST')

  CALL gterm('GSSH')
  CALL grad_ssh(gx, gy)
  CALL gterm('GSSH')


!$OMP PARALLEL DO PRIVATE(km0)
  DO k=1, ksize
     km0 = maskindex(k)

     DO j=1, jsize
     DO i=0, isize
        u(i,j,k)  = u(i,j,k) + dtime*gx(i,j,k)
     END DO
     END DO

     DO j=0, jsize
     DO i=1, isize
        v(i,j,k)  = v(i,j,k) + dtime*gy(i,j,k)
     END DO
     END DO
  END DO

  IF (use_landwater .AND. vrank==0) THEN
     DO j=1, jsize
     DO i=0, isize
        u(i,j,ksize) = u(i,j,ksize) - imask3d(i,j,ksize)*imask3d(i+1,j,ksize) * gravity * (landelev(i+1,j)-landelev(i,j))*idx1(i,j) * dtime
        IF (u(i,j,ksize) > 0.0 .AND. lwdried(i,j))   u(i,j,ksize) = 0.0D0
        IF (u(i,j,ksize) < 0.0 .AND. lwdried(i+1,j)) u(i,j,ksize) = 0.0D0
     END DO
     END DO

     DO j=0, jsize
     DO i=1, isize
        v(i,j,ksize) = v(i,j,ksize) - imask3d(i,j,ksize)*imask3d(i,j+1,ksize) * gravity * (landelev(i,j+1)-landelev(i,j))*idy1(i,j) * dtime
        IF (v(i,j,ksize) > 0.0 .AND. lwdried(i,j))   v(i,j,ksize) = 0.0D0
        IF (v(i,j,ksize) < 0.0 .AND. lwdried(i,j+1)) v(i,j,ksize) = 0.0D0
     END DO
     END DO
  END IF

  IF (hydrostatic) THEN
     CALL step_hydrostatic
  ELSE ! Non-hydrostatic
     CALL step_nonhydrostatic
  END IF

  IF (.NOT. rigid_lid) THEN
     ssh_old(:,:) = ssh(:,:)

!     IF (vrank==0) ssh(1:isize,1:jsize) = ssh(1:isize,1:jsize) + 2.0*p2d(1:isize,1:jsize)/(dtime*gravity)
     IF (vrank==0) ssh(1:isize,1:jsize) = ssh(1:isize,1:jsize) + 2.0*p2d(1:isize,1:jsize)/(dtime*gravity)

     CALL vcast(ssh)
!!!!!!YJKIM
    IF (jcoord==jpes-1) THEN
       DO i=1, isize
!          ssh(i,jsize-1)=ssh(i,jsize-2)
          ssh(i,jsize)=ssh(i,jsize-1)
       END DO
    END IF
    IF (jcoord==0) THEN
      DO i=1, isize
         ssh(i,2)=ssh(i,3)
         ssh(i,1)=ssh(i,2)
         ssh(i,0)=ssh(i,1)
      END DO
   END IF       
    IF (icoord==ipes-1) THEN
      DO j=1, jsize
!         ssh(isize-1,j)=ssh(isize-2,j)
         ssh(isize,j)=ssh(isize-1,j)
      END DO   
    END IF


     CALL update_boundary(ssh)

     CALL update_geometry
     CALL checkout_geometry

  END IF

  CALL update_velocity_boundary(u, v, w)

  IF (seaice_coupled) THEN
     IF (vrank==0) CALL coriolis_2d(ui, vi, gix_ab, giy_ab)

     IF (vrank==0) THEN
        CALL seaice_rheology(ui, vi, ai, hi, gix, giy)

        DO j=1, jsize
        DO i=0, isize
           gix(i,j) = gix(i,j) + 1.5D0 * gix_ab(i,j) - gravity * (ssh(i+1,j) - ssh(i,j)) * idx1(i,j)

           ui(i,j) = ui(i,j) + imask2d(i+1,j)*imask2d(i,j) * gix(i,j)*dtime
        END DO
        END DO

        DO j=0, jsize
        DO i=1, isize
           giy(i,j) = giy(i,j) + 1.5D0 * giy_ab(i,j) - gravity * (ssh(i,j+1) - ssh(i,j)) * idy1(i,j)

           vi(i,j) = vi(i,j) + imask2d(i,j)*imask2d(i,j+1) * giy(i,j)*dtime
        END DO
        END DO
     END IF

     CALL update_boundary(ui)
     CALL update_boundary(vi)
  END IF

  IF (iceshelf_coupled) THEN
     CALL iceshelf_thermodynamics
  END IF

  CALL step_tracers

  IF (seaice_coupled) THEN
     CALL seaice_advection(ui, vi, hi, ai)
  END IF

  IF (tracer_index_t /= 0) CALL heatflux(tracer(:,:,:,tracer_index_t))

  IF (seaice_coupled) THEN
     CALL seaice_thermodynamics(tracer(:,:,:,tracer_index_t), tracer(:,:,:,tracer_index_s), hi, ai)
  END IF

  CALL update_density

  IF (use_npzd)      CALL step_npzd
  IF (use_ecosystem) CALL step_ecosystem

PROFILE_BEGIN('particle')
  IF (use_particles)  CALL step_particles
PROFILE_END('particle')

  CALL openboundary_velocity

  IF (assert_nan) CALL assert(check_nan(u) .AND. check_nan(v) .AND. check_nan(w), "detect NaN in prognostic velocity, abort!")

  CALL checkout_velocity

  CALL update_velocity_diagnostic

  CALL checkout_tracers

  IF (seaice_coupled) CALL checkout_seaice

  t_current = t_current + dtime
  current_datetime = format_datetime(t_current)

  CALL flush_io

  IF (use_particles) CALL flush_particles

  CALL report_ssh
  CALL report_tracers
  CALL report_open
  CALL report_momentum
  CALL report_energy
  CALL report_cfl

  IF (REPORT_UNIT > 0 .AND. REPORT_UNIT /= STDOUT_UNIT) FLUSH(REPORT_UNIT)

CONTAINS
  SUBROUTINE step_nonhydrostatic
    USE solver3d

!$OMP PARALLEL DO PRIVATE(km0, km1)
    DO k=0, ksize
       km0=maskindex(k)
       km1=maskindex(k+1)

       DO j=1, jsize
       DO i=1, isize
          w(i,j,k) = w(i,j,k) + dtime * gz(i,j,k)
       END DO
       END DO
    END DO

    IF (cycle_z) CALL update_boundary(w)

!$OMP PARALLEL DO PRIVATE(kl)
    DO k=1, ksize
       kl=dzindex(k)

       DO j=1, jsize
       DO i=1, isize
          q3d(i,j,k) = imask3d(i,j,k)*( u(i,j,k)*dsx(i,j,kl) - u(i-1,j,k)*dsx(i-1,j,kl) &
                                      + v(i,j,k)*dsy(i,j,kl) - v(i,j-1,k)*dsy(i,j-1,kl) &
                                      + w(i,j,k)*dsz(i,j)    - w(i,j,k-1)*dsz(i,j))     &
                     + imask3d(i,j,k)*( u_old(i,j,k)*dsx(i,j,kl) - u_old(i-1,j,k)*dsx(i-1,j,kl) &
                                      + v_old(i,j,k)*dsy(i,j,kl) - v_old(i,j-1,k)*dsy(i,j-1,kl) &
                                      + w_old(i,j,k)*dsz(i,j)    - w_old(i,j,k-1)*dsz(i,j)) * c_ebcn
       END DO
       END DO
    END DO

    IF (report_qsum) THEN
       qsum = sum(q3d(1:isize,1:jsize,1:ksize))
       CALL gsum(qsum)
       IF (rank==0) WRITE(REPORT_UNIT, '(A,ES12.5)') 'qsum = ', qsum
    END IF

    CALL checkout('Q3D', q3d) ! source term for 3D pressure solver [m^3/s]

    IF (bypass_solver) RETURN

    IF (.NOT. rigid_lid) THEN
       IF (require_checkin('SURFACE_ALPHA')) THEN
          tmp(:,:,1) = 1.0D0
          CALL checkin('SURFACE_ALPHA', tmp(:,:,1))
          DO k=0, ksize+1
          DO j=0, jsize+1
          DO i=0, isize+1
                da_solver(i,j,k) = 2.0*dsz(i,j)/(gravity*dtime**2) * tmp(i,j,1) * b_ebcn * imask3d(i,j,k)*cz_star(i,j,k)
          END DO
          END DO
          END DO
       END IF

PROFILE_BEGIN('reset_solver')
       SELECT CASE (pgradient_scheme)
       CASE (1)
          CALL reset_solver3d(dx_h=dx, dy_h=dy, dz_v=dz, dsx=dsx, dsy=dsy, dsz_h=dsz, &
                              dvol=dvol, alpha=da_solver, mask=lmask3d)
       CASE (2,3)
          CALL reset_solver3d(dx_h=dx, dy_h=dy, dz=dz_star, dsx=dsx, dsy=dsy, dsz_h=dsz, &
                              dvol=dvol, alpha=da_solver, mask=lmask3d)
       END SELECT

PROFILE_END('reset_solver')
    END IF

PROFILE_BEGIN('solver3d')
!$OMP PARALLEL WORKSHARE
    p3d(:,:,:) = p3d(:,:,:)*imask3d(0:isize+1,0:jsize+1,0:ksize+1)
!$OMP END PARALLEL WORKSHARE

    CALL solve3d(q3d, p3d)
PROFILE_END('solver3d')
    CALL checkout('P3D', p3d) ! solution for 3D pressure solver [m^2/s]

    p2d(:,:)= 0.0D0
    DO j=0, jsize+1
    DO i=0, isize+1
       IF (surface_flag(i,j)) THEN
          k = surface_k(i,j)
          p2d(i,j) = p3d(i,j,k)
          IF (pgradient_scheme==3) THEN
             p2d(i,j) = p3d(i,j,k) / dvol(i,j,k)
          ELSE
             p2d(i,j) = p3d(i,j,k)
          END IF
       END IF
    END DO
    END DO


    CALL vsum(p2d, all=.TRUE.)
    CALL checkout('P2D', p2d) ! solution for 2D pressure solver [m^2/s]

    SELECT CASE (pgradient_scheme)
    CASE (1)
!$OMP PARALLEL DO PRIVATE(i, j, k, km0, km1)
       DO k=0, ksize
          km0=maskindex(k)
          km1=maskindex(k+1)
          DO j=1, jsize
          DO i=1, isize
             w(i,j,k) = w(i,j,k) - imask3d(i,j,km0)*imask3d(i,j,  km1) * (p3d(i,j,k+1)-p3d(i,j,k))*idz1(k)
          END DO
          END DO

          IF (k==0) CYCLE

          DO j=1, jsize
          DO i=0, isize
             u(i,j,k) = u(i,j,k) - imask3d(i,j,km0)*imask3d(i+1,j,km0) * (p3d(i+1,j,k)-p3d(i,j,k))*idx1(i,j)
          END DO
          END DO

          DO j=0, jsize
          DO i=1, isize
             v(i,j,k) = v(i,j,k) - imask3d(i,j,km0)*imask3d(i,j+1,km0) * (p3d(i,j+1,k)-p3d(i,j,k))*idy1(i,j)
          END DO
          END DO
       END DO
    CASE (2)
!$OMP PARALLEL DO PRIVATE(i, j, k, km0, km1)
       DO k=0, ksize
          km0=maskindex(k)
          km1=maskindex(k+1)

          DO j=1, jsize
          DO i=1, isize
             w(i,j,k) = w(i,j,k) - imask3d(i,j,km0)*imask3d(i,j,km1) &
                  * (p3d(i,j,k+1) - p3d(i,j,k))*dsz(i,j) * 2.0/(dvol(i,j,k)+dvol(i,j,k+1))
          END DO
          END DO

          IF (k==0) CYCLE

          DO j=1, jsize
          DO i=0, isize
             u(i,j,k) = u(i,j,k) - imask3d(i,j,km0)*imask3d(i+1,j,km0) &
                  * (p3d(i+1,j,k)*dy(i+1,j)*dz_star(i+1,j,k) - p3d(i,j,k)*dy(i,j)*dz_star(i,j,k)) &
                  * 2.0/(dvol(i,j,k)+dvol(i+1,j,k))
          END DO
          END DO

          DO j=0, jsize
          DO i=1, isize
             v(i,j,k) = v(i,j,k) - imask3d(i,j,km0)*imask3d(i,j+1,km0) &
                  * (p3d(i,j+1,k)*dx(i,j+1)*dz_star(i,j+1,k) - p3d(i,j,k)*dx(i,j)*dz_star(i,j,k)) &
                  * 2.0/(dvol(i,j,k)+dvol(i,j+1,k))
          END DO
          END DO
       END DO
    CASE (3)
!$OMP PARALLEL DO PRIVATE(i, j, k, km0, km1)
       DO k=0, ksize
          km0=maskindex(k)
          km1=maskindex(k+1)

          DO j=1, jsize
          DO i=1, isize
             w(i,j,k) = w(i,j,k) - imask3d(i,j,km0)*imask3d(i,j,km1) * (p3d(i,j,k+1)/dz_star(i,j,k+1)-p3d(i,j,k)/dz_star(i,j,k)) * 2.0/(dvol(i,j,k)+dvol(i,j,k+1))
          END DO
          END DO

          IF (k==0) CYCLE

          DO j=1, jsize
          DO i=0, isize
             u(i,j,k) = u(i,j,k) - imask3d(i,j,km0)*imask3d(i+1,j,km0) * (p3d(i+1,j,k)*idx0(i+1,j)-p3d(i,j,k)*idx0(i,j)) * 2.0/(dvol(i,j,k)+dvol(i+1,j,k))
          END DO
          END DO

          DO j=0, jsize
          DO i=1, isize
             v(i,j,k) = v(i,j,k) - imask3d(i,j,km0)*imask3d(i,j+1,km0) * (p3d(i,j+1,k)*idy0(i,j+1)-p3d(i,j,k)*idy0(i,j)) * 2.0/(dvol(i,j,k)+dvol(i,j+1,k))
          END DO
          END DO

       END DO
    END SELECT

    IF (require_checkout('GX_GNHP')) THEN
       SELECT CASE (pgradient_scheme)
       CASE (1)
!$OMP PARALLEL DO
          DO k=1, ksize
          DO j=1, jsize
          DO i=0, isize
             tmp(i,j,k) = imask3d(i,j,k)*imask3d(i+1,j,k)*(p3d(i,j,k)-p3d(i+1,j,k))*idx1(i,j)*idtime
          END DO
          END DO
          END DO
       CASE (2)
!$OMP PARALLEL DO
          DO k=1, ksize
          DO j=1, jsize
          DO i=0, isize
             tmp(i,j,k) = imask3d(i,j,k)*imask3d(i+1,j,k)*(p3d(i,j,k)*dy(i,j)*dz_star(i,j,k)-p3d(i+1,j,k)*dy(i+1,j)*dz_star(i+1,j,k)) * 2.0/(dvol(i,j,k)+dvol(i+1,j,k)) * idtime
          END DO
          END DO
          END DO
       CASE (3)
!$OMP PARALLEL DO
          DO k=1, ksize
          DO j=1, jsize
          DO i=0, isize
             tmp(i,j,k) = imask3d(i,j,k)*imask3d(i+1,j,k)*(p3d(i,j,k)*idx0(i,j) - p3d(i+1,j,k)*idx0(i+1,j)) * 2.0/(dvol(i,j,k)+dvol(i+1,j,k)) * idtime
          END DO
          END DO
          END DO
       END SELECT
       CALL checkout('GX_GNHP', tmp(0:isize,1:jsize,1:ksize)) ! x-compoenent for nonhydrostatic pressure gradient [m/s^2]
    END IF
    IF (require_checkout('GY_GNHP')) THEN
       SELECT CASE (pgradient_scheme)
       CASE (1)
!$OMP PARALLEL DO
          DO k=1, ksize
          DO j=0, jsize
          DO i=1, isize
             tmp(i,j,k) = imask3d(i,j,k)*imask3d(i,j+1,k)*(p3d(i,j,k)-p3d(i,j+1,k))*idy1(i,j)*idtime
          END DO
          END DO
          END DO
       CASE (2)
!$OMP PARALLEL DO
          DO k=1, ksize
          DO j=1, jsize
          DO i=0, isize
             tmp(i,j,k) = imask3d(i,j,k)*imask3d(i,j+1,k)*(p3d(i,j,k)*dx(i,j)*dz_star(i,j,k)-p3d(i,j+1,k)*dx(i,j+1)*dz_star(i,j+1,k)) * 2.0/(dvol(i,j,k)+dvol(i,j+1,k)) * idtime
          END DO
          END DO
          END DO
       CASE (3)
!$OMP PARALLEL DO
          DO k=1, ksize
          DO j=1, jsize
          DO i=0, isize
             tmp(i,j,k) = imask3d(i,j,k)*imask3d(i,j+1,k)*(p3d(i,j,k)*idy0(i,j)-p3d(i,j+1,k)*idy0(i,j+1)) * 2.0/(dvol(i,j,k)+dvol(i,j+1,k)) * idtime
          END DO
          END DO
          END DO
       END SELECT
       CALL checkout('GY_GNHP', tmp(1:isize,0:jsize,1:ksize)) ! y-compoenent for nonhydrostatic pressure gradient [m/s^2]
    END IF
    IF (require_checkout('GZ_GNHP')) THEN
       SELECT CASE (pgradient_scheme)
       CASE (1)
!$OMP PARALLEL DO
          DO k=0, ksize
          DO j=1, jsize
          DO i=1, isize
             tmp(i,j,k) = imask3d(i,j,k)*imask3d(i,j,k+1)*(p3d(i,j,k)-p3d(i,j,k+1))*idz1(k)*idtime
          END DO
          END DO
          END DO
       CASE (2)
!$OMP PARALLEL DO
          DO k=1, ksize
          DO j=1, jsize
          DO i=0, isize
             tmp(i,j,k) = imask3d(i,j,k)*imask3d(i,j+1,k)*(p3d(i,j,k)*dx(i,j)*dy(i,j)-p3d(i,j,k+1)*dx(i,j)*dy(i,j)) * 2.0/(dvol(i,j,k)+dvol(i,j,k+1)) * idtime
          END DO
          END DO
          END DO
       CASE (3)
!$OMP PARALLEL DO
          DO k=1, ksize
          DO j=1, jsize
          DO i=0, isize
             tmp(i,j,k) = imask3d(i,j,k)*imask3d(i,j+1,k)*(p3d(i,j,k)/dz_star(i,j,k)-p3d(i,j,k+1)/dz_star(i,j,k+1)) * 2.0/(dvol(i,j,k)+dvol(i,j,k+1)) * idtime
          END DO
          END DO
          END DO
       END SELECT
       CALL checkout('GZ_GNHP', tmp(1:isize,1:jsize,0:ksize)) ! z-compoenent for nonhydrostatic pressure gradient [m/s^2]
    END IF

  END SUBROUTINE step_nonhydrostatic

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE step_hydrostatic
    USE solver2d

    REAL(8) :: u2d(0:isize, 1:jsize)
    REAL(8) :: v2d(1:isize, 0:jsize)
    REAL(8) :: w2d(1:isize, 1:jsize, 0:ksize)

    INTEGER :: i, j, k
    INTEGER :: kk, kl, km

!$OMP PARALLEL PRIVATE(i, j, k, kk, kl, km)
!$OMP WORKSHARE
    u2d(:,:) = 0.0D0
    v2d(:,:) = 0.0D0
!$OMP END WORKSHARE

    DO k=1, ksize
       kk=k
       kl=dzindex(k)

!$OMP DO
       DO j=1, jsize
       DO i=0, isize
          u2d(i,j) = u2d(i,j) + (u(i,j,kk) + u_old(i,j,kk)*c_ebcn) * dsx(i,j,kl)
       END DO
       END DO

!$OMP DO
       DO j=0, jsize
       DO i=1, isize
          v2d(i,j) = v2d(i,j) + (v(i,j,kk) + v_old(i,j,kk)*c_ebcn) * dsy(i,j,kl)
       END DO
       END DO
    END DO
!$OMP END PARALLEL

    CALL vsum(u2d, all=.TRUE.)
    CALL vsum(v2d, all=.TRUE.)

!$OMP PARALLEL DO
    DO j=1, jsize
    DO i=1, isize
       q2d(i,j) = (u2d(i,j) - u2d(i-1,j) + v2d(i,j) - v2d(i,j-1))
    END DO
    END DO

    CALL checkout('Q2D', q2d) ! source term for 2D pressure solver [m^3/s]

    IF (.NOT. bypass_solver) THEN
       IF (vrank==0) THEN
          IF (.NOT. rigid_lid) THEN
             IF (require_checkin('SURFACE_ALPHA')) THEN
                tmp(:,:,1) = 1.0D0
                CALL checkin('SURFACE_ALPHA', tmp(:,:,1))
                DO j=0, jsize+1
                DO i=0, isize+1
                   da_solver(i,j,1) = 2.0*dsz(i,j)/(gravity*dtime**2) * tmp(i,j,1) * b_ebcn
                END DO
                END DO
             END IF

PROFILE_BEGIN('reset_solver')
             CALL reset_solver2d(dx(0:isize+1,0:jsize+1),    &
                                 dy(0:isize+1,0:jsize+1),    &
                                 dsx2d(0:isize+1,0:jsize+1), &
                                 dsy2d(0:isize+1,0:jsize+1), &
                                 da_solver(:,:,1),           &
                                 lmask2d(0:isize+1,0:jsize+1))
PROFILE_END('reset_solver')
          END IF

PROFILE_BEGIN('solver2d')
!$OMP PARALLEL WORKSHARE
          p2d = p2d(:,:)*imask2d(0:isize+1,0:jsize+1)
!$OMP END PARALLEL WORKSHARE
          CALL solve2d(q2d, p2d)
PROFILE_END('solver2d')
       END IF

       CALL vcast(p2d)
       CALL checkout('P2D', p2d) ! solution for 2D puressure solver [m^2/s]
    END IF

!$OMP PARALLEL DO PRIVATE(km0)
    DO k=1, ksize
       km0 = maskindex(k)

       DO j=1, jsize
       DO i=0, isize
          u(i,j,k) = u(i,j,k) - imask3d(i,j,km0)*imask3d(i+1,j,km0)*(p2d(i+1,j)-p2d(i,j))*idx1(i,j)
       END DO
       END DO

       DO j=0, jsize
       DO i=1, isize
          v(i,j,k) = v(i,j,k) - imask3d(i,j,km0)*imask3d(i,j+1,km0)*(p2d(i,j+1)-p2d(i,j))*idy1(i,j)
       END DO
       END DO
    END DO

    CALL update_w(u, v, w)

  END SUBROUTINE step_hydrostatic

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE gterm(name, ab)
    CHARACTER(*), INTENT(IN) :: name
    LOGICAL,      INTENT(IN), OPTIONAL :: ab

    CHARACTER(8), SAVE :: name_  = ''
    LOGICAL,      SAVE :: switch = .FALSE.
    LOGICAL :: ab_

    REAL(8), ALLOCATABLE, SAVE :: gx_tmp(:,:,:)
    REAL(8), ALLOCATABLE, SAVE :: gy_tmp(:,:,:)
    REAL(8), ALLOCATABLE, SAVE :: gz_tmp(:,:,:)

    ab_ = .FALSE.
    IF (present(ab)) ab_ = ab

    IF (.NOT. allocated(gx_tmp)) ALLOCATE(gx_tmp(0:isize,1:jsize,1:ksize))
    IF (.NOT. allocated(gy_tmp)) ALLOCATE(gy_tmp(1:isize,0:jsize,1:ksize))
    IF (.NOT. allocated(gz_tmp)) ALLOCATE(gz_tmp(1:isize,1:jsize,0:ksize))

    IF (switch) THEN
       CALL assert(trim(name) == trim(name_), "illegal call of GTERM")

       IF (require_checkout('GX_' // trim(name_))) THEN
          IF (ab_) THEN
             CALL axpy(1.0D0, gx_ab, gx_tmp)
          ELSE
             CALL axpy(1.0D0, gx,    gx_tmp)
          END IF
          CALL checkout('GX_'//trim(name_), gx_tmp) ! x-component for a certain tendency term [m/s^2]
       END IF

       IF (require_checkout('GY_' // trim(name_))) THEN
          IF (ab_) THEN
             CALL axpy(1.0D0, gy_ab, gy_tmp)
          ELSE
             CALL axpy(1.0D0, gy,    gy_tmp)
          END IF
          CALL checkout('GY_'//trim(name_), gy_tmp) ! y-component for a certain tendency term [m/s^2]
       END IF

       IF (require_checkout('GZ_' // trim(name_))) THEN
          IF (ab_) THEN
             CALL axpy(1.0D0, gz_ab, gz_tmp)
          ELSE
             CALL axpy(1.0D0, gz,    gz_tmp)
          END IF
          CALL checkout('GZ_'//trim(name_), gz_tmp) ! z-component for a certain tendency term [m/s^2]
       END IF

#ifdef DEBUG
       DO k=1, ksize
       DO j=1, jsize
       DO i=0, isize
          IF (imask3d(i,j,k)*imask3d(i+1,j,k)==0.0 .AND. &
               ((ab_ .AND. gx_ab(i,j,k)/=0.0) .OR. (.NOT.ab_ .AND. gx(i,j,k)/=0.0))) THEN
             WRITE(0,*) trim(name_) // ": GX mask inconsistent at ", i, j, k
             STOP
          END IF
          END DO
          END DO
       END DO
       DO k=1, ksize
       DO j=0, jsize
       DO i=1, isize
          IF (imask3d(i,j,k)*imask3d(i,j+1,k)==0.0 .AND. &
               ((ab_ .AND. gy_ab(i,j,k)/=0.0) .OR. (.NOT.ab_ .AND. gy(i,j,k)/=0.0))) THEN
             WRITE(0,*) trim(name_) // "GY mask inconsistent at ", i, j, k
             STOP
          END IF
       END DO
       END DO
       END DO
       DO k=0, ksize
       DO j=1, jsize
       DO i=1, isize
          IF (imask3d(i,j,k)*imask3d(i,j,k+1)==0.0 .AND. &
               ((ab_ .AND. gz_ab(i,j,k)/=0.0) .OR. (.NOT.ab_ .AND. gz(i,j,k)/=0.0))) THEN
             WRITE(0,*) trim(name_) // "GZ mask inconsistent at ", i, j, k
             STOP
          END IF
       END DO
       END DO
       END DO
#endif
    ELSE
       name_ = trim(name)

       IF (require_checkout('GX_' // trim(name_))) THEN
          IF (ab_) THEN
             CALL scal(-1.0D0, gx_ab, gx_tmp)
          ELSE
             CALL scal(-1.0D0, gx,    gx_tmp)
          END IF
       END IF

       IF (require_checkout('GY_' // trim(name_))) THEN
          IF (ab_) THEN
             CALL scal(-1.0D0, gy_ab, gy_tmp)
          ELSE
             CALL scal(-1.0D0, gy,    gy_tmp)
          END IF
       END IF

       IF (require_checkout('GZ_' // trim(name_))) THEN
          IF (ab_) THEN
             CALL scal(-1.0D0, gz_ab, gz_tmp)
          ELSE
             CALL scal(-1.0D0, gz,    gz_tmp)
          END IF
       END IF
    END IF

    switch = .NOT. switch
  END SUBROUTINE gterm

  SUBROUTINE dump
    CALL dump_data(u)
    CALL dump_data(v)
    CALL dump_data(w)
    CALL dump_data(ssh)
    CALL dump_data(tracer)
  END SUBROUTINE dump

END SUBROUTINE step

SUBROUTINE checkout_geometry
  USE misc
  USE geometry
  USE io
  IMPLICIT NONE

  CALL checkout('SSH', ssh) ! sea surface height [m]

  CALL checkout('CEXT', cext) ! phase speed of external gravity wave [m/s]

  IF (require_checkout('DVOL')) CALL checkout_dvol
  IF (require_checkout('DSX'))  CALL checkout_dsx
  IF (require_checkout('DSY'))  CALL checkout_dsy

  CALL checkout('BATHYMETRY', h_bathymetry) ! bathymetry [m]
  CALL checkout('ICESHELF',   h_iceshelf)   ! iceshelf base [m]

  IF (use_landwater .AND. require_checkout('DRYWET')) CALL checkout_drywet
  IF (require_checkout('PMASK'))                      CALL checkout_pmask

CONTAINS
  SUBROUTINE checkout_dvol
    REAL(8) :: tmp(isize,jsize,ksize)
    INTEGER :: i, j, k

!$OMP PARALLEL DO
    DO k=1, ksize
    DO j=1, jsize
    DO i=1, isize
       tmp(i,j,k) = dvol(i,j,k)*imask3d(i,j,k)
    END DO
    END DO
    END DO

    CALL checkout('DVOL', tmp) ! grid volume [m^3]
  END SUBROUTINE checkout_dvol

  SUBROUTINE checkout_dsx
    REAL(8) :: tmp(0:isize,jsize,ksize)
    INTEGER :: i, j, k

!$OMP PARALLEL DO
    DO k=1, ksize
    DO j=1, jsize
    DO i=0, isize
       tmp(i,j,k) = dsx(i,j,k)*imask3d(i,j,k)*imask3d(i+1,j,k)
    END DO
    END DO
    END DO

    IF (open_w .AND. icoord==0) THEN
       DO k=1, ksize
       DO j=1, jsize
          tmp(0,j,k) = dsx(0,j,k)*imask3d(1,j,k)
       END DO
       END DO
    END IF

    IF (open_e .AND. icoord==ipes-1) THEN
       DO k=1, ksize
       DO j=1, jsize
          tmp(isize,j,k) = dsx(isize,j,k)*imask3d(isize,j,k)
       END DO
       END DO
    END IF

    CALL checkout('DSX', tmp) ! grid area of x-section [m^2]
  END SUBROUTINE checkout_dsx

  SUBROUTINE checkout_dsy
    REAL(8) :: tmp(isize,0:jsize,ksize)
    INTEGER :: i, j, k

!$OMP PARALLEL DO
    DO k=1, ksize
    DO j=0, jsize
    DO i=1, isize
       tmp(i,j,k) = dsy(i,j,k)*imask3d(i,j,k)*imask3d(i,j+1,k)
    END DO
    END DO
    END DO

    IF (open_s .AND. jcoord==0) THEN
       DO k=1, ksize
       DO i=1, isize
          tmp(i,0,k) = dsy(i,0,k)*imask3d(i,1,k)
       END DO
       END DO
    END IF

    IF (open_n .AND. jcoord==jpes-1) THEN
       DO k=1, ksize
       DO i=1, isize
          tmp(i,jsize,k) = dsy(i,jsize,k)*imask3d(i,jsize,k)
       END DO
       END DO
    END IF

    CALL checkout('DSY', tmp) ! grid area of y-section [m^2]
  END SUBROUTINE checkout_dsy

  SUBROUTINE checkout_drywet
    REAL(8) :: tmp(isize,jsize)
    INTEGER :: i, j

    IF (vrank==0) THEN
       DO j=1, jsize
       DO i=1, isize
          IF (.NOT. lmask3d(i,j,ksize)) THEN
             tmp(i,j) = -1.0   ! masked cell
          ELSE IF (.NOT. lwflag(i,j)) THEN
             tmp(i,j) =  2.0   ! ocean cell
          ELSE IF (.NOT. lwdried(i,j)) THEN
             tmp(i,j) =  1.0   ! wet land cell
          ELSE
             tmp(i,j) =  0.0   ! dry land cell
          END IF
       END DO
       END DO
    ELSE
       tmp(:,:) = 0.0
    END IF

    CALL checkout('DRYWET', tmp) ! flags for dry/wet grid
  END SUBROUTINE checkout_drywet

  SUBROUTINE checkout_pmask
    REAL(8) :: tmp(isize,jsize,ksize)
    INTEGER :: i, j, k

    DO k=1, ksize
    DO j=1, jsize
    DO i=1, isize
       tmp(i,j,k) = imask3d(i,j,k) * dvol(i,j,k) / (dx(i,j)*dy(i,j)*dz(k))
    END DO
    END DO
    END DO

    CALL checkout('PMASK', tmp) ! mask with partial cell ratio
  END SUBROUTINE checkout_pmask

END SUBROUTINE checkout_geometry

