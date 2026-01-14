#include "macro.h"

SUBROUTINE step_offline
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

  INTEGER :: i, j, k
  INTEGER :: n

  LOGICAL :: stat

  n_timestep = n_timestep + 1

  CALL update_flux
  CALL update_viscosity

!$OMP PARALLEL WORKSHARE
  u_old(:,:,:) = u(:,:,:)
  v_old(:,:,:) = v(:,:,:)
  w_old(:,:,:) = w(:,:,:)
!$OMP END PARALLEL WORKSHARE

  IF (.NOT. rigid_lid) THEN
     ssh_old(:,:) = ssh(:,:)
     CALL checkin('SSH', ssh, stat)
     IF (stat) CALL update_geometry
     CALL checkout_geometry
  END IF

  CALL checkin('U', u, stat)
  CALL assert(stat, "offline mode requires input of 'U'")

  CALL checkin('V', v, stat)
  CALL assert(stat, "offline mode requires input of 'V'")

  CALL checkin('W', w, stat)
  IF (.NOT. stat) THEN
     CALL update_w(u, v, w, topdown=w_topdown)
     CALL update_boundary(w)
  END IF
  CALL update_velocity_external

  CALl step_tracers

  IF (seaice_coupled) THEN
     CALL seaice_advection(ui, vi, hi, ai)
  END IF

  IF (tracer_index_t /= 0) CALL heatflux(tracer(:,:,:,tracer_index_t))

  CALL update_density

  IF (use_npzd)      CALL step_npzd
  IF (use_ecosystem) CALL step_ecosystem

PROFILE_BEGIN('particle')
  IF (use_particles) CALL step_particles
PROFILE_END('particle')

  IF (assert_nan) CALL assert(check_nan(u) .AND. check_nan(v) .AND. check_nan(w), "detect NaN in velocity for offline mode")

  CALL checkout_velocity

  CALL update_velocity_diagnostic

  CALL checkout_tracers

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

!-----------------------------------------------------------------------------------------------------------------------

END SUBROUTINE step_offline
