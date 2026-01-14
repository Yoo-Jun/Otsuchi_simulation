#include "macro.h"

SUBROUTINE iceshelf_thermodynamics
  USE misc
  USE geometry
  USE parameters
  USE state
  USE tracers
  USE io
  IMPLICIT NONE

  REAL(8) :: ismelt(isize,jsize,ksize)

  REAL(8) :: gt, gs
  REAL(8) :: t, s, z, m, ds
  REAL(8) :: tf, tf_s, sb, t_insitu
  REAL(8) :: a, b, c

  INTEGER :: tracer_index_mwtrc

  INTEGER :: i, j, k

  IF (.NOT. iceshelf_melt) RETURN

  IF (iceshelf_gamma_t <= 0.0) RETURN

  CALL assert(tracer_index_t/=0 .AND. tracer_index_s/=0, "tracers 'T' and 'S' are required for iceshelf thermodynamics")

  gt = iceshelf_gamma_t
  gs = iceshelf_gamma_s

!$OMP PARALLEL PRIVATE(i, j, k) &
!$OMP          PRIVATE(t, s, z, m, ds)  &
!$OMP          PRIVATE(tf, tf_s, sb, a, b, c)
!$OMP DO
  DO k=1, ksize
  DO j=1, jsize
  DO i=1, isize
     IF (.NOT. lmask3d(i,j,k)) THEN
        ismelt(i,j,k) = 0.0
        CYCLE
     END IF

     ds = 0.0
     IF (isflag(i,j,k+1)) ds = ds + dsz(i,j)
     IF (isflag(i-1,j,k)) ds = ds + dsx_ref(i-1,j,k)
     IF (isflag(i+1,j,k)) ds = ds + dsx_ref(i,  j,k)
     IF (isflag(i,j-1,k)) ds = ds + dsy_ref(i,j-1,k)
     IF (isflag(i,j+1,k)) ds = ds + dsy_ref(i,j,  k)

     IF (ds == 0.0) THEN
        ismelt(i,j,k) = 0.0
        CYCLE
     END IF

     t = tracer(i,j,k,tracer_index_t)
     s = tracer(i,j,k,tracer_index_s)
     z = depth_offset + depth(k)

     !use linear approximation of freezing_temperature
     tf   = freezing_temperature(s, z)
     tf_s = freezing_temperature(s+1.0, z) - tf
     tf = tf - tf_s*s

     t_insitu = potential_temperature(t, s, 0.0D0, z)

     a = cp * gt * tf_s
     b = cp * gt * (tf - t_insitu) - gs * l_freeze
     c = gs * l_freeze * s

     sb = (- b - sqrt(b**2 - 4*a*c)) / (2*a)
     tf = tf + tf_s * sb
     ismelt(i,j,k) = rho_0*cp*gt*(t_insitu - tf)/l_freeze * ds/dvol(i,j,k)

     m = ismelt(i,j,k)*dtime

     tracer(i,j,k,tracer_index_t) = t - m*l_freeze/(cp*rho_0)
     tracer(i,j,k,tracer_index_s) = s*rho_0/(rho_0+m)
  END DO
  END DO
  END DO
!$OMP END PARALLEL

  CALL checkout('ISMELT', ismelt) ![kg / (s m^3)]

  CALL update_tracer_boundary(tracer_index_t)
  CALL update_tracer_boundary(tracer_index_s)

  tracer_index_mwtrc = tracer_index('MWTRC')
  IF (tracer_index_mwtrc /= 0) THEN
!OMP PARALLEL DO
     DO k=1, ksize
     DO j=1, jsize
     DO i=1, isize
        tracer(i,j,k,tracer_index_mwtrc) = tracer(i,j,k,tracer_index_mwtrc) + ismelt(i,j,k)*dtime
        ![kg /m^3] contents of iceshelf-melt origin freshwater
     END DO
     END DO
     END DO

     CALL update_tracer_boundary(tracer_index_mwtrc)
  END IF


END SUBROUTINE iceshelf_thermodynamics
