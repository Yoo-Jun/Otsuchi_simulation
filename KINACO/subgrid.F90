#include "macro.h"

MODULE subgrid
  USE misc
  USE geometry
  USE io
  USE velocity
  IMPLICIT NONE
  SAVE

  REAL(8), ALLOCATABLE :: visc(:,:,:,:)
!!!!!!!!!YJKIM  
!!!  REAL(8) :: sponge_2d(0:isize+1,0:jsize+1)
!!!!!!!!!!!!!
  INTEGER :: les_report_interval = -1


CONTAINS
  SUBROUTINE init_subgrid
    ALLOCATE(visc(0:isize+1, 0:jsize+1, 0:ksize+1, 3))
    visc(:,:,:,:) = 0.0
!!!!!!!!!YJKIM
!!!    CALL initial_data('SPONGE_2D', sponge_2d)
!!!!!!!!!!
  END SUBROUTINE init_subgrid

  SUBROUTINE finalize_subgrid
    DEALLOCATE(visc)
  END SUBROUTINE finalize_subgrid

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE viscosity(gx, gy, gz)
    USE parameters, ONLY: viscosity_scheme, visch, viscv, c_smagorinsky_visc

    REAL(8), INTENT(INOUT) :: gx(0:isize, 1:jsize, 1:ksize)
    REAL(8), INTENT(INOUT) :: gy(1:isize, 0:jsize, 1:ksize)
    REAL(8), INTENT(INOUT) :: gz(1:isize, 1:jsize, 0:ksize)

    REAL(8) :: eps(1:isize, 1:jsize, 1:ksize)

    IF (require_checkout('DRATE')) THEN
       CALL assert(viscosity_scheme==1, "Sorry, 'DRATE' output supports isotorpic visicosity (VISCOSITY_MODE=1) only ")
!$OMP PARALLEL WORKSHARE
       eps(:,:,:) = 0.0
!$OMP END PARALLEL WORKSHARE
    ELSE
       eps(1,1,1) = UNDEF
    END IF

    SELECT CASE(viscosity_scheme)
    CASE (0)
       ! no viscosity term
    CASE (1)
       CALL viscosity_iso(visc(:,:,:,1), gx, gy, gz, eps)
    CASE (2)
       CALL viscosity_h(visc(:,:,:,1), gx, gy, gz)
       CALL viscosity_v(visc(:,:,:,2), gx, gy, gz)
    CASE (3)
       CALL assert(.FALSE., "sopycnal/diapycnal viscosity is not yet supported.")
    CASE (4)
       CALL viscosity_h(   visc(:,:,:,1), gx, gy, gz)
       CALL viscosity_vimp(visc(:,:,:,2), gx, gy)
    END SELECT

    IF (use_gls) CALL viscosity_vimp(visc(:,:,:,3), gx, gy)

    CALL checkout('DRATE', eps)

  END SUBROUTINE viscosity

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE viscosity_iso(nu, gx, gy, gz, eps)
    USE parameters, ONLY: rho_0
    REAL(8), INTENT(IN) :: nu(0:isize+1,0:jsize+1,0:ksize+1)

    REAL(8), INTENT(INOUT) :: gx(0:isize, 1:jsize, 1:ksize)
    REAL(8), INTENT(INOUT) :: gy(1:isize, 0:jsize, 1:ksize)
    REAL(8), INTENT(INOUT) :: gz(1:isize, 1:jsize, 0:ksize)

    REAL(8), INTENT(INOUT) :: eps(isize, jsize, ksize)

    REAL(8) :: tau_xx(0:isize+1, 1:jsize)
    REAL(8) :: tau_yy(1:isize,   0:jsize+1)
    REAL(8) :: tau_zz(1:isize,   1:jsize, 0:1)

    REAL(8) :: tau_xy(0:isize, 0:jsize)
    REAL(8) :: tau_yz(1:isize, 0:jsize, 0:1)
    REAL(8) :: tau_zx(0:isize, 1:jsize, 0:1)

    REAL(8) :: tmp(0:isize, 0:jsize)
    LOGICAL :: flag_eps

    INTEGER :: i, j, k
    INTEGER :: kstr, kend
    INTEGER :: kl0, kl1
    INTEGER :: km0, km1

    INTEGER :: tid

    tid  = 0
    kstr = 1
    kend = ksize

    flag_eps = (eps(1,1,1) /= UNDEF)

!$OMP PARALLEL PRIVATE(i, j, k, tid)  &
!$OMP          PRIVATE(tau_xx, tau_yy, tau_zz) &
!$OMP          PRIVATE(tau_xy, tau_yz, tau_zx) &
!$OMP          PRIVATE(kstr, kend) &
!$OMP          PRIVATE(kl0, kl1, km0, km1)

!$  tid  = omp_get_thread_num()
!$  kstr = kstr_t(tid)
!$  kend = kend_t(tid)

    tau_zz(:,:,1) = UNDEF
    tau_zx(:,:,1) = UNDEF
    tau_yz(:,:,1) = UNDEF

    DO k=kstr-2, kend
       km0 = maskindex(k)
       km1 = maskindex(k+1)
       kl0 = dzindex(k)
       kl1 = dzindex(k+1)

       DO j=1, jsize
       DO i=1, isize
          tau_zz(i,j,0) = tau_zz(i,j,1)
          tau_zz(i,j,1) = 2.0*nu(i,j,k+1)*dwdz(i,j,k+1)
       END DO
       END DO

       IF (k==kstr-2) CYCLE

       DO j=0, jsize
       DO i=1, isize
          tau_yz(i,j,0) = tau_yz(i,j,1)
          tau_yz(i,j,1) = 0.25*(nu(i,j,k)+nu(i,j+1,k)+nu(i,j,k+1)+nu(i,j+1,k+1)) * (dwdy(i,j,k)+dvdz(i,j,k))
       END DO
       END DO

       DO j=1, jsize
       DO i=0, isize
          tau_zx(i,j,0) = tau_zx(i,j,1)
          tau_zx(i,j,1) = 0.25*(nu(i,j,k)+nu(i+1,j,k)+nu(i,j,k+1)+nu(i+1,j,k+1)) * (dudz(i,j,k)+dwdx(i,j,k))
       END DO
       END DO

!$     IF (tid/=0 .AND. k==kstr-1) CYCLE

       DO j=1, jsize
       DO i=1, isize
          gz(i,j,k) = gz(i,j,k) + imask3d(i,j,km0)*imask3d(i,j,km1) &
               * (  tau_zx(i,j,1)*(dsx(i,j,kl0)+dsx(i,j,kl1)) - tau_zx(i-1,j,  1)*(dsx(i-1,j,  kl0)+dsx(i-1,j,  kl1)) &
                  + tau_yz(i,j,1)*(dsy(i,j,kl0)+dsy(i,j,kl1)) - tau_yz(i,  j-1,1)*(dsy(i,  j-1,kl0)+dsy(i,  j-1,kl1)) &
                  + tau_zz(i,j,1)*2*dsz(i,j)                  - tau_zz(i,  j,  0)*2*dsz(i,j))                         &
               / (dvol(i,j,kl0)+dvol(i,j,kl1))
       END DO
       END DO

       IF (k==kstr-1) CYCLE

       DO j=1, jsize
       DO i=0, isize+1
          tau_xx(i,j) = 2.0 * nu(i,j,k) * dudx(i,j,k)
       END DO
       END DO

       DO j=0, jsize+1
       DO i=1, isize
          tau_yy(i,j) = 2.0 * nu(i,j,k) * dvdy(i,j,k)
       END DO
       END DO

       DO j=0, jsize
       DO i=0, isize
          tau_xy(i,j) = 0.25 * (nu(i,j,k)+nu(i+1,j,k)+nu(i,j+1,k)+nu(i+1,j+1,k)) * (dvdx(i,j,k)+dudy(i,j,k))
       END DO
       END DO

       DO j=1, jsize
       DO i=0, isize
          gx(i,j,k) = gx(i,j,k) + imask3d(i,j,km0)*imask3d(i+1,j,km0) &
               * (  tau_xx(i+1,j)*(dsx(i,j,kl0)+dsx(i+1,j,kl0)) - tau_xx(i,j)  *(dsx(i-1,j,  kl0)+dsx(i,  j,  kl0)) &
                  + tau_xy(i,  j)*(dsy(i,j,kl0)+dsy(i+1,j,kl0)) - tau_xy(i,j-1)*(dsy(i,  j-1,kl0)+dsy(i+1,j-1,kl0)) &
                  + tau_zx(i,j,1)*(dsz(i,j)    +dsz(i+1,j))     - tau_zx(i,j,0)*(dsz(i,j)        +dsz(i+1,j)))      &
               / (dvol(i,j,kl0)+dvol(i+1,j,kl0))

          IF (k == min(surface_k(i,j), surface_k(i+1,j))) taux_sfc(i,j) = taux_sfc(i,j) + rho_0 * tau_zx(i,j,1)
          IF (k == max( bottom_k(i,j),  bottom_k(i+1,j))) taux_btm(i,j) = taux_btm(i,j) - rho_0 * tau_zx(i,j,0)
       END DO
       END DO

       DO j=0, jsize
       DO i=1, isize
          gy(i,j,k) = gy(i,j,k) + imask3d(i,j,km0)*imask3d(i,j+1,km0) &
               * (  tau_xy(i,j)  *(dsx(i,j,kl0)+dsx(i,j+1,kl0)) - tau_xy(i-1,j)*(dsx(i-1,j,  kl0)+dsx(i-1,j+1,kl0)) &
                  + tau_yy(i,j+1)*(dsy(i,j,kl0)+dsy(i,j+1,kl0)) - tau_yy(i,  j)*(dsy(i,  j-1,kl0)+dsy(i,  j,  kl0)) &
                  + tau_yz(i,j,1)*(dsz(i,j)    +dsz(i,j+1))     - tau_yz(i,j,0)*(dsz(i,j)        +dsz(i,j+1)))      &
               / (dvol(i,j,kl0)+dvol(i,j+1,kl0))

          IF (k == min(surface_k(i,j), surface_k(i,j+1))) tauy_sfc(i,j) = tauy_sfc(i,j) + rho_0 * tau_yz(i,j,1)
          IF (k == max( bottom_k(i,j),  bottom_k(i,j+1))) tauy_btm(i,j) = tauy_btm(i,j) - rho_0 * tau_yz(i,j,0)
       END DO
       END DO

       IF (flag_eps) THEN
          DO j=0, jsize
          DO i=0, isize
             tmp(i,j) = 0.25*tau_xy(i,j)*(dudy(i,j,k)+dvdx(i,j,k))
          END DO
          END DO

          DO j=1, jsize
          DO i=1, isize
             eps(i,j,k) = eps(i,j,k) + tau_xx(i,j)*dudx(i,j,k) + tau_yy(i,j)*dvdy(i,j,k) + tau_zz(i,j,0)*dwdz(i,j,k) &
                                     + tmp(i-1,j-1) + tmp(i,j-1) + tmp(i-1,j) + tmp(i,j)
          END DO
          END DO

          DO j=1, jsize
          DO i=0, isize
             tmp(i,j) = 0.25*(tau_zx(i,j,0)*(dudz(i,j,k-1)+dwdx(i,j,k-1)) + tau_zx(i,j,1)*(dudz(i,j,k)+dwdx(i,j,k)))
          END DO
          END DO

          DO j=1, jsize
          DO i=1, isize
             eps(i,j,k) = eps(i,j,k) + tmp(i-1,j) + tmp(i,j)
          END DO
          END DO

          DO j=0, jsize
          DO i=1, isize
             tmp(i,j) = 0.25*(tau_yz(i,j,0)*(dvdz(i,j,k-1)+dwdy(i,j,k-1)) + tau_yz(i,j,1)*(dvdz(i,j,k)+dwdy(i,j,k)))
          END DO
          END DO

          DO j=1, jsize
          DO i=1, isize
             eps(i,j,k) = eps(i,j,k) + tmp(i,j-1) + tmp(i,j)
          END DO
          END DO
       END IF
    END DO
!$OMP END PARALLEL

  END SUBROUTINE viscosity_iso

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE viscosity_h(nu, gx, gy, gz)
    USE parameters, ONLY: rho_0

    REAL(8), INTENT(IN) :: nu(0:isize+1,0:jsize+1,0:ksize+1)

    REAL(8), INTENT(INOUT) :: gx(0:isize, 1:jsize, 1:ksize)
    REAL(8), INTENT(INOUT) :: gy(1:isize, 0:jsize, 1:ksize)
    REAL(8), INTENT(INOUT) :: gz(1:isize, 1:jsize, 0:ksize)

    REAL(8) :: fx(0:isize+1, 0:jsize+1)
    REAL(8) :: fy(0:isize+1, 0:jsize+1)

    INTEGER :: i, j, k
    INTEGER :: kstr, kend
    INTEGER :: km0, km1
    INTEGER :: kl0, kl1

    INTEGER :: tid

    tid  = 0
    kstr = 1
    kend = ksize

!$OMP PARALLEL PRIVATE(i, j, k, tid) &
!$OMP          PRIVATE(fx, fy)       &
!$OMP          PRIVATE(kstr, kend)   &
!$OMP          PRIVATE(km0, km1, kl0, kl1)

!$  tid  = omp_get_thread_num()
!$  kstr = kstr_t(tid)
!$  kend = kend_t(tid)

    DO k=kstr, kend
       kl0 = dzindex(k)
       km0 = maskindex(k)

       DO j=1, jsize
       DO i=0, isize+1
          fx(i,j) = -nu(i,j,k)*dudx(i,j,k)
       END DO
       END DO

       DO j=0, jsize
       DO i=0, isize
          fy(i,j) = -0.25*(nu(i,j,k)+nu(i+1,j,k)+nu(i,j+1,k)+nu(i+1,j+1,k)) * dudy(i,j,k)
       END DO
       END DO

       DO j=1, jsize
       DO i=0, isize
          gx(i,j,k) = gx(i,j,k) - imask3d(i,j,km0)*imask3d(i+1,j,km0) &
               * (  fx(i+1,j)*(dsx(i,j,kl0)+dsx(i+1,j,kl0)) - fx(i,j)  *(dsx(i-1,j,  kl0)+dsx(i,  j,  kl0))  &
                  + fy(i,  j)*(dsy(i,j,kl0)+dsy(i+1,j,kl0)) - fy(i,j-1)*(dsy(i,  j-1,kl0)+dsy(i+1,j-1,kl0))) &
               / (dvol(i,j,kl0)+dvol(i+1,j,kl0))
       END DO
       END DO
    END DO

    DO k=kstr, kend
       kl0 = dzindex(k)
       km0 = maskindex(k)

       DO j=0, jsize
       DO i=0, isize
          fx(i,j) = -0.25*(nu(i,j,k)+nu(i+1,j,k)+nu(i,j+1,k)+nu(i+1,j+1,k)) * dvdx(i,j,k)
       END DO
       END DO

       DO j=0, jsize+1
       DO i=1, isize
          fy(i,j) = -nu(i,j,k)*dvdy(i,j,k)
       END DO
       END DO

       DO j=0, jsize
       DO i=1, isize
          gy(i,j,k) = gy(i,j,k) - imask3d(i,j,km0)*imask3d(i,j+1,km0) &
               * (  fx(i,j)  *(dsx(i,j,kl0)+dsx(i,j+1,kl0)) - fx(i-1,j)*(dsx(i-1,j,  kl0)+dsx(i-1,j+1,kl0))  &
                  + fy(i,j+1)*(dsy(i,j,kl0)+dsy(i,j+1,kl0)) - fy(i,  j)*(dsy(i,  j-1,kl0)+dsy(i,  j,  kl0))) &
               / (dvol(i,j,kl0)+dvol(i,j+1,kl0))
       END DO
       END DO
    END DO

    DO k=kstr-2, kend
       kl0 = dzindex(k)
       kl1 = dzindex(k+1)
       km0 = maskindex(k)
       km1 = maskindex(k+1)

       IF (k==kstr-2) CYCLE
!$     IF (tid/=0 .AND. k==kstr-1) CYCLE

       DO j=1, jsize
       DO i=0, isize
          fx(i,j) = -0.25*(nu(i,j,k)+nu(i+1,j,k)+nu(i,j,k+1)+nu(i+1,j,k+1)) * dwdx(i,j,k)
       END DO
       END DO

       DO j=0, jsize
       DO i=1, isize
          fy(i,j) = -0.25*(nu(i,j,k)+nu(i,j+1,k)+nu(i,j,k+1)+nu(i,j+1,k+1)) * dwdy(i,j,k)
       END DO
       END DO

       DO j=1, jsize
       DO i=1, isize
          gz(i,j,k) = gz(i,j,k) - imask3d(i,j,km0)*imask3d(i,j,km1) &
               * (  fx(i,j)  *(dsx(i,j,kl0)+dsx(i,j,kl1)) - fx(i-1,j)*(dsx(i-1,j,  kl0)+dsx(i-1,j,  kl1))  &
                  + fy(i,j)  *(dsy(i,j,kl0)+dsy(i,j,kl1)) - fy(i,j-1)*(dsy(i,  j-1,kl0)+dsy(i,  j-1,kl1))) &
               / (dvol(i,j,kl0)+dvol(i,j,kl1))
       END DO
       END DO
    END DO
!$OMP END PARALLEL

  END SUBROUTINE viscosity_h

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE viscosity_v(nu, gx, gy, gz)
    USE parameters, ONLY: rho_0

    REAL(8), INTENT(IN) :: nu(0:isize+1,0:jsize+1,0:ksize+1)

    REAL(8), INTENT(INOUT) :: gx(0:isize, 1:jsize, 1:ksize)
    REAL(8), INTENT(INOUT) :: gy(1:isize, 0:jsize, 1:ksize)
    REAL(8), INTENT(INOUT) :: gz(1:isize, 1:jsize, 0:ksize)

    REAL(8) :: fz(0:isize+1, 0:jsize+1, 0:1)

    INTEGER :: i, j, k
    INTEGER :: kstr, kend
    INTEGER :: km0, km1
    INTEGER :: kl0, kl1

    INTEGER :: tid

    tid  = 0
    kstr = 1
    kend = ksize

!$OMP PARALLEL PRIVATE(i, j, k, tid) &
!$OMP          PRIVATE(fz)           &
!$OMP          PRIVATE(kstr, kend)   &
!$OMP          PRIVATE(km0, km1, kl0, kl1)

!$  tid  = omp_get_thread_num()
!$  kstr = kstr_t(tid)
!$  kend = kend_t(tid)

    fz(:,:,1) = UNDEF

    DO k=kstr-1, kend
       kl0 = dzindex(k)
       km0 = maskindex(k)

       DO j=1, jsize
       DO i=0, isize
          fz(i,j,0) = fz(i,j,1)
          fz(i,j,1) = -0.25*(nu(i,j,k)+nu(i+1,j,k)+nu(i,j,k+1)+nu(i+1,j,k+1)) * dudz(i,j,k)
       END DO
       END DO

       IF (k==kstr-1) CYCLE

       DO j=1, jsize
       DO i=0, isize
          gx(i,j,k) = gx(i,j,k) - imask3d(i,j,km0)*imask3d(i+1,j,km0) &
               * (fz(i,j,1)-fz(i,j,0)) * (dsz(i,j)+dsz(i+1,j)) / (dvol(i,j,kl0)+dvol(i+1,j,kl0))

          IF (k == min(surface_k(i,j), surface_k(i+1,j))) taux_sfc(i,j) = taux_sfc(i,j) + rho_0 * fz(i,j,1)
          IF (k == max( bottom_k(i,j),  bottom_k(i+1,j))) taux_btm(i,j) = taux_btm(i,j) - rho_0 * fz(i,j,0)
       END DO
       END DO
    END DO

    DO k=kstr-1, kend
       kl0 = dzindex(k)
       km0 = maskindex(k)

       DO j=0, jsize
       DO i=1, isize
          fz(i,j,0) = fz(i,j,1)
          fz(i,j,1) = -0.25*(nu(i,j,k)+nu(i,j+1,k)+nu(i,j,k+1)+nu(i,j+1,k+1)) * dvdz(i,j,k)
       END DO
       END DO

       IF (k==kstr-1) CYCLE

       DO j=0, jsize
       DO i=1, isize
          gy(i,j,k) = gy(i,j,k) - imask3d(i,j,km0)*imask3d(i,j+1,km0) &
               * (fz(i,j,1)-fz(i,j,0)) * (dsz(i,j)+dsz(i,j+1)) / (dvol(i,j,kl0)+dvol(i,j+1,kl0))

          IF (k == min(surface_k(i,j), surface_k(i,j+1))) tauy_sfc(i,j) = tauy_sfc(i,j) + rho_0 * fz(i,j,1)
          IF (k == max( bottom_k(i,j),  bottom_k(i,j+1))) tauy_btm(i,j) = tauy_btm(i,j) - rho_0 * fz(i,j,0)
       END DO
       END DO
    END DO

    DO k=kstr-2, kend
       kl0 = dzindex(k)
       kl1 = dzindex(k+1)
       km0 = maskindex(k)
       km1 = maskindex(k+1)

       DO j=1, jsize
       DO i=1, isize
          fz(i,j,0) = fz(i,j,1)
          fz(i,j,1) = -nu(i,j,k+1)*dwdz(i,j,k+1)
       END DO
       END DO

       IF (k==kstr-2) CYCLE
!$     IF (tid/=0 .AND. k==kstr-1) CYCLE

       DO j=1, jsize
       DO i=1, isize
          gz(i,j,k) = gz(i,j,k) - imask3d(i,j,km0)*imask3d(i,j,km1) &
               * (fz(i,j,1)-fz(i,j,0)) * dsz(i,j) * 2.0 / (dvol(i,j,kl0)+dvol(i,j,kl1))
       END DO
       END DO
    END DO
!$OMP END PARALLEL

  END SUBROUTINE viscosity_v

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE viscosity_vimp(nu, gx, gy)
    USE parameters, ONLY: rho_0

    REAL(8), INTENT(IN) :: nu(0:isize+1,0:jsize+1,0:ksize+1)

    REAL(8), INTENT(INOUT) :: gx(0:isize, 1:jsize, 1:ksize)
    REAL(8), INTENT(INOUT) :: gy(1:isize, 0:jsize, 1:ksize)

    REAL(8) :: x(0:isize,1:jsize,0:ksize+1)
    REAL(8) :: y(1:isize,0:jsize,0:ksize+1)

    REAL(8) :: a(0:isize,1:jsize,0:ksize)
    REAL(8) :: b(1:isize,0:jsize,0:ksize)

    REAL(8) :: w(0:isize,0:jsize,ksize+1,2)
    REAL(8) :: t(0:isize,0:jsize,2)

    REAL(8) :: dzu(0:isize,1:jsize,1:ksize)
    REAL(8) :: dzv(1:isize,0:jsize,1:ksize)

    INTEGER :: i, j, k
    INTEGER :: km0, km1

    CALL assert(.NOT. cycle_z, "implicit vertical viscosity is not supported for CYCLE_Z")

!$OMP PARALLEL PRIVATE(i,j,k, km0, km1)
!$OMP DO
    DO k=1, ksize
       DO j=1, jsize
       DO i=0, isize
          dzu(i,j,k) = (dvol(i,j,k)+dvol(i+1,j,k)) / (dsz(i,j)+dsz(i+1,j))
       END DO
       END DO

       DO j=0, jsize
       DO i=1, isize
          dzv(i,j,k) = (dvol(i,j,k)+dvol(i,j+1,k)) / (dsz(i,j)+dsz(i,j+1))
       END DO
       END DO
    END DO

!$OMP DO
    DO k=0, ksize
       km0 = maskindex(k)
       km1 = maskindex(k+1)

       DO j=1, jsize
       DO i=0, isize
          a(i,j,k) = 0.25*(nu(i,j,k)+nu(i+1,j,k)+nu(i,j,k+1)+nu(i+1,j,k+1))*idz1(k)&
                         *imask3d(i,j,km0)*imask3d(i+1,j,km0)*imask3d(i,j,km1)*imask3d(i+1,j,km1)*dtime
       END DO
       END DO

       DO j=0, jsize
       DO i=1, isize
          b(i,j,k) = 0.25*(nu(i,j,k)+nu(i,j+1,k)+nu(i,j,k+1)+nu(i,j+1,k+1))*idz1(k)&
                         *imask3d(i,j,km0)*imask3d(i,j+1,km0)*imask3d(i,j,km1)*imask3d(i,j+1,km1)*dtime
       END DO
       END DO
    END DO

!$OMP WORKSHARE
    t(:,:,:) = 1.0

    x(:,:,0) = 0.0
    y(:,:,0) = 0.0

    x(:,:,ksize+1) = 0.0
    y(:,:,ksize+1) = 0.0
!$OMP END WORKSHARE

!$OMP END PARALLEL

    CALL lrecv(t,        tag=1)
    CALL lrecv(x(:,:,0), tag=2)
    CALL lrecv(y(:,:,0), tag=3)

!$OMP PARALLEL PRIVATE(i,j)
    DO k=1, ksize
!$OMP DO COLLAPSE(2)
       DO j=1, jsize
       DO i=0, isize
          w(i,j,k,1) = -a(i,j,k-1) / t(i,j,1)
          t(i,j,1)   =  dzu(i,j,k) + a(i,j,k-1) + a(i,j,k) + a(i,j,k-1)*w(i,j,k,1)
          x(i,j,k)   = (dzu(i,j,k)*u(i,j,k) + a(i,j,k-1)*x(i,j,k-1)) / t(i,j,1)
       END DO
       END DO
!$OMP DO COLLAPSE(2)
       DO j=0, jsize
       DO i=1, isize
          w(i,j,k,2) = -b(i,j,k-1) / t(i,j,2)
          t(i,j,2)   =  dzv(i,j,k) + b(i,j,k-1) + b(i,j,k) + b(i,j,k-1)*w(i,j,k,2)
          y(i,j,k)   = (dzv(i,j,k)*v(i,j,k) + b(i,j,k-1)*y(i,j,k-1)) / t(i,j,2)
       END DO
       END DO
    END DO

!$OMP DO COLLAPSE(2)
    DO j=1, jsize
    DO i=0, isize
       w(i,j,ksize+1,1) = -a(i,j,ksize) / t(i,j,1)
    END DO
    END DO
!$OMP DO COLLAPSE(2)
    DO j=0, jsize
    DO i=1, isize
       w(i,j,ksize+1,2) = -b(i,j,ksize) / t(i,j,2)
    END DO
    END DO
!$OMP END PARALLEL

    CALL usend(t,              tag=1)
    CALL usend(x(:,:,ksize),   tag=2)
    CALL usend(y(:,:,ksize),   tag=3)

    CALL urecv(x(:,:,ksize+1), tag=4)
    CALL urecv(y(:,:,ksize+1), tag=5)

!$OMP PARALLEL PRIVATE(i,j)
    DO k=ksize, 0, -1
       IF (k==ksize .AND. vrank==0) CYCLE

!$OMP DO COLLAPSE(2)
       DO j=1, jsize
       DO i=0, isize
          x(i,j,k)  = x(i,j,k) - w(i,j,k+1,1)*x(i,j,k+1)
       END DO
       END DO
!$OMP DO COLLAPSE(2)
       DO j=0, jsize
       DO i=1, isize
          y(i,j,k)  = y(i,j,k) - w(i,j,k+1,2)*y(i,j,k+1)
       END DO
       END DO
    END DO
!$OMP END PARALLEL

    CALL lsend(x(:,:,1), tag=4)
    CALL lsend(y(:,:,1), tag=5)

!$OMP PARALLEL DO PRIVATE(km0, km1)
    DO k=1, ksize
       km0 = maskindex(k)

       DO j=1, jsize
       DO i=0, isize
          gx(i,j,k) = gx(i,j,k) + (a(i,j,k)*(x(i,j,k+1)-x(i,j,k)) + a(i,j,k-1)*(x(i,j,k-1)-x(i,j,k))) * imask3d(i,j,km0)*imask3d(i+1,j,km0) * idtime/dzu(i,j,k)
       END DO
       END DO

       DO j=0, jsize
       DO i=1, isize
          gy(i,j,k) = gy(i,j,k) + (b(i,j,k)*(y(i,j,k+1)-y(i,j,k)) + b(i,j,k-1)*(y(i,j,k-1)-y(i,j,k))) * imask3d(i,j,km0)*imask3d(i,j+1,km0) * idtime/dzv(i,j,k)
       END DO
       END DO
    END DO

  END SUBROUTINE viscosity_vimp

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE update_viscosity
    INTEGER :: i, j, k

!$OMP PARALLEL WORKSHARE
    visc(:,:,:,1) = visch
    visc(:,:,:,2) = viscv
!$OMP END PARALLEL WORKSHARE

    CALL checkin('VISC_H', visc(:,:,:,1), add=.TRUE.) !for backward_compatibility
    CALL checkin('VISC_V', visc(:,:,:,2), add=.TRUE.) !for backward_compatibility

    CALL checkin('VISCH', visc(:,:,:,1), add=.TRUE.)
    CALL checkin('VISCV', visc(:,:,:,2), add=.TRUE.)

    IF (c_smagorinsky_visc /= 0.0) THEN 
!$OMP PARALLEL DO
       DO k=0, ksize+1
       DO j=0, jsize+1
       DO i=0, isize+1
          visc(i,j,k,1) = visc(i,j,k,1) + c_smagorinsky_visc**2 * delta2(i,j,k) * sqrt2 * srate(i,j,k)
       END DO
       END DO
       END DO
    END IF

    IF (les_scheme /= 0) THEN
       CALL les_viscosity(visc(:,:,:,3))

!$OMP PARALLEL DO
       DO k=0, ksize+1
       DO j=0, jsize+1
       DO i=0, isize+1
          visc(i,j,k,1) = visc(i,j,k,1) + visc(i,j,k,3)
          visc(i,j,k,2) = visc(i,j,k,2) + visc(i,j,k,3)
       END DO
       END DO
       END DO
    END IF

    CALL checkout('VISCH', visc(:,:,:,1)) ! horizontal viscosity [m^2/s]
    CALL checkout('VISCV', visc(:,:,:,2)) ! vertical viscosity   [m^2/s]

  END SUBROUTINE update_viscosity

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE les_viscosity(visc)
    REAL(8), INTENT(OUT) :: visc(0:isize+1, 0:jsize+1, 0:ksize+1)
    REAL(8) :: max_, sum_
    INTEGER :: i, j, k

    SELECT CASE (les_scheme)
    CASE (1) ! simple Smagorinsky model
       CALL les_smagorinsky(visc)

    CASE (2) ! the WALE model (Nicoud and Ducros 1999)
       CALL les_wale(visc)

    CASE (3) ! the FSF model (Ducros et al. 1996)
       CALL les_fsf(visc)
    END SELECT

    IF (limit_les /= UNDEF) THEN
!$OMP PARALLEL DO
       DO k=0, ksize+1
       DO j=0, jsize+1
       DO i=0, isize+1
          visc(i,j,k) = min(visc(i,j,k),limit_les)
       END DO
       END DO
       END DO
    END IF

    IF (les_report_interval > 0 .AND. mod(n_timestep, les_report_interval) == 0) THEN
       max_  = maxval(visc(1:isize,1:jsize,1:ksize), mask=lmask3d(1:isize,1:jsize,1:ksize))
       sum_  = sum(   visc(1:isize,1:jsize,1:ksize), mask=lmask3d(1:isize,1:jsize,1:ksize))

       CALL gmax(max_)
       CALL gsum(sum_)

       IF (rank==0) THEN
          WRITE(REPORT_UNIT, '(A, ": LES viscosity max=", ES10.3, ", mean=", ES10.3)') &
               trim(current_datetime), max_, sum_/n_cells
       END IF
    END IF

  CONTAINS
    SUBROUTINE les_smagorinsky(visc)
      REAL(8), INTENT(OUT) :: visc(0:isize+1, 0:jsize+1, 0:ksize+1)
      INTEGER :: i, j, k

      REAL(4) :: damp(0:isize+1, 0:jsize+1, 0:ksize+1)
      LOGICAL :: stat

      CALL checkin("SMAGDAMP", damp, stat=stat)
      IF (.NOT. stat) THEN
!$OMP PARALLEL WORKSHARE
         damp(:,:,:) = 1.0
!$OMP END PARALLEL WORKSHARE
      END IF

!$OMP PARALLEL DO
       DO k=0, ksize+1
       DO j=0, jsize+1
       DO i=0, isize+1
          visc(i,j,k) = c_smagorinsky_les**2 * delta2(i,j,k)*(damp(i,j,k)**2) * sqrt2 * srate(i,j,k)
       END DO
       END DO
       END DO

     END SUBROUTINE les_smagorinsky

     SUBROUTINE les_wale(visc)
       REAL(8), INTENT(OUT) :: visc(0:isize+1, 0:jsize+1, 0:ksize+1)
       REAL(8) :: s11, s22, s33, s12, s23, s31, r12, r23, r31, ss, rr, dd
       INTEGER :: i, j, k

!$OMP PARALLEL DO PRIVATE(s11, s22, s33, s12, s23, s31, r12, r23, r31, ss, rr, dd)
       DO k=0, ksize+1
       DO j=0, jsize+1
       DO i=0, isize+1
          IF (.NOT. lmask3d(i,j,k)) THEN
             visc(i,j,k) = 0.0
             CYCLE
          END IF

          s11 = dudx(i,j,k)
          s22 = dvdy(i,j,k)
          s33 = dwdz(i,j,k)

          s12 = 0.125*(dudy(i-1,j-1,k)+dudy(i,j-1,k)+dudy(i-1,j,k)+dudy(i,j,k) &
                      +dvdx(i-1,j-1,k)+dvdx(i,j-1,k)+dvdx(i-1,j,k)+dvdx(i,j,k))

          s23 = 0.125*(dvdz(i,j-1,k-1)+dvdz(i,j,k-1)+dvdz(i,j-1,k)+dvdz(i,j,k) &
                      +dwdy(i,j-1,k-1)+dwdy(i,j,k-1)+dwdy(i,j-1,k)+dwdy(i,j,k))

          s31 = 0.125*(dwdx(i-1,j,k-1)+dwdx(i-1,j,k)+dwdx(i,j,k-1)+dwdx(i,j,k) &
                      +dudz(i-1,j,k-1)+dudz(i-1,j,k)+dudz(i,j,k-1)+dudz(i,j,k))

          r12 = 0.125*(dudy(i-1,j-1,k)+dudy(i,j-1,k)+dudy(i-1,j,k)+dudy(i,j,k) &
                      -dvdx(i-1,j-1,k)-dvdx(i,j-1,k)-dvdx(i-1,j,k)-dvdx(i,j,k))

          r23 = 0.125*(dvdz(i,j-1,k-1)+dvdz(i,j,k-1)+dvdz(i,j-1,k)+dvdz(i,j,k) &
                      -dwdy(i,j-1,k-1)-dwdy(i,j,k-1)-dwdy(i,j-1,k)-dwdy(i,j,k))

          r31 = 0.125*(dwdx(i-1,j,k-1)+dwdx(i-1,j,k)+dwdx(i,j,k-1)+dwdx(i,j,k) &
                      -dudz(i-1,j,k-1)-dudz(i-1,j,k)-dudz(i,j,k-1)-dudz(i,j,k))

          ss = s11**2 + s22**2 + s33**2 + 2*(s12**2 + s23**2 + s31**2)
          rr =                            2*(r12**2 + r23**2 + r31**2)

          dd = -(s11**2 + s12**2 + s31**2)*(r12**2+r31**2)  &
               -(s12**2 + s22**2 + s23**2)*(r23**2+r12**2)  &
               -(s31**2 + s23**2 + s33**2)*(r31**2+r23**2)  &
               + 2*(s12*(s11+s22)+s23*s31)*r23*r31 &
               + 2*(s31*(s33+s11)+s12*s23)*r12*r23 &
               + 2*(s23*(s22+s33)+s31*s12)*r31*r12

          dd = 1/6.0*(ss**2 + rr**2) + 2/3.0*ss*rr + 2*dd

          IF (dd < 1.0E-24) THEN
             visc(i,j,k) = 0.0
          ELSE
             visc(i,j,k) = c_wale**2 * delta2(i,j,k) * dd**1.5  / (ss**2.5 + dd**1.25)
          END IF
       END DO
       END DO
       END DO
     END SUBROUTINE les_wale

     SUBROUTINE les_fsf(visc)
       REAL(8), INTENT(OUT) :: visc(0:isize+1, 0:jsize+1, 0:ksize+1)
       INTEGER :: i, j, k

       REAL(8) :: fu(-1:isize+2, -1:jsize+2, -1:ksize+2)
       REAL(8) :: fv(-1:isize+2, -1:jsize+2, -1:ksize+2)
       REAL(8) :: fw(-1:isize+2, -1:jsize+2, -1:ksize+2)
       REAL(8) :: tmp(-1:isize+2, -1:jsize+2, -1:ksize+2)
       REAL(8) :: fsf
       INTEGER :: n
       INTEGER :: km0, km1, km2

!$OMP PARALLEL DO PRIVATE(km1)
       DO k=-1, ksize+2
          km1 = maskindex(k)
       DO j=-1, jsize+2
       DO i=-1, isize+2
          tmp(i,j,k) = imask3d(i,j,km1)*(u(i-1,j,k)+u(i,j,k))*0.5
       END DO
       END DO
       END DO

       CALL lfilter(tmp, fu, 1)
       CALL lfilter(fu, tmp, 0)
       CALL update_boundary(tmp)
       CALL lfilter(tmp, fu, 1)

!$OMP PARALLEL DO PRIVATE(km1)
       DO k=-1, ksize+2
          km1 = maskindex(k)
       DO j=-1, jsize+2
       DO i=-1, isize+2
          tmp(i,j,k) = imask3d(i,j,km1)*(v(i,j-1,k)+v(i,j,k))*0.5
       END DO
       END DO
       END DO
       CALL lfilter(tmp, fv, 1)
       CALL lfilter(fv, tmp, 0)
       CALL update_boundary(tmp)
       CALL lfilter(tmp, fv, 1)

!$OMP PARALLEL DO PRIVATE(km1)
       DO k=-1, ksize+2
          km1 = maskindex(k)
       DO j=-1, jsize+2
       DO i=-1, isize+2
          tmp(i,j,k) = imask3d(i,j,km1)*(w(i,j,k-1)+w(i,j,k))*0.5
       END DO
       END DO
       END DO
       CALL lfilter(tmp, fw, 1)
       CALL lfilter(fw, tmp, 0)
       CALL update_boundary(tmp)
       CALL lfilter(tmp, fw, 1)

!$OMP PARALLEL DO PRIVATE(fsf, n, km0, km1, km2)
       DO k=1, ksize
          km0 = maskindex(k-1)
          km1 = maskindex(k)
          km2 = maskindex(k+1)
       DO j=1, jsize
       DO i=1, isize
          fsf = 0.0
          n   = 0
          IF (.NOT. lmask3d(i,j,km1)) THEN
             visc(i,j,k) = 0.0
             CYCLE
          END IF

          IF (lmask3d(i+1,j,km1)) THEN
             fsf = fsf + ((fu(i+1,j,k)-fu(i,j,k))**2 + (fv(i+1,j,k)-fv(i,j,k))**2 + (fw(i+1,j,k)-fw(i,j,k))**2) * (delta2(i,j,k)*idx0(i,j)**2)**(1.0/3)
             n = n+1
          END IF
          IF (lmask3d(i-1,j,km1)) THEN
             fsf = fsf + ((fu(i-1,j,k)-fu(i,j,k))**2 + (fv(i-1,j,k)-fv(i,j,k))**2 + (fw(i-1,j,k)-fw(i,j,k))**2) * (delta2(i,j,k)*idx0(i,j)**2)**(1.0/3)
             n = n+1
          END IF
          IF (lmask3d(i,j+1,km1)) THEN
             fsf = fsf + ((fu(i,j+1,k)-fu(i,j,k))**2 + (fv(i,j+1,k)-fv(i,j,k))**2 + (fw(i,j+1,k)-fw(i,j,k))**2) * (delta2(i,j,k)*idy0(i,j)**2)**(1.0/3)
             n = n+1
          END IF
          IF (lmask3d(i,j-1,km1)) THEN
             fsf = fsf + ((fu(i,j-1,k)-fu(i,j,k))**2 + (fv(i,j-1,k)-fv(i,j,k))**2 + (fw(i,j-1,k)-fw(i,j,k))**2) * (delta2(i,j,k)*idy0(i,j)**2)**(1.0/3)
             n = n+1
          END IF
          IF (lmask3d(i,j,km2)) THEN
             fsf = fsf + ((fu(i,j,k+1)-fu(i,j,k))**2 + (fv(i,j,k+1)-fv(i,j,k))**2 + (fw(i,j,k+1)-fw(i,j,k))**2) * (delta2(i,j,k)*idz0(k)**2)**(1.0/3)
             n = n+1
          END IF
          IF (lmask3d(i,j,km0)) THEN
             fsf = fsf + ((fu(i,j,k-1)-fu(i,j,k))**2 + (fv(i,j,k-1)-fv(i,j,k))**2 + (fw(i,j,k-1)-fw(i,j,k))**2) * (delta2(i,j,k)*idz0(k)**2)**(1.0/3)
             n = n+1
          END IF

          IF (n>0) THEN
             fsf = fsf/n
             visc(i,j,k) = c_fsf * sqrt(delta2(i,j,k)*fsf)
          ELSE
             visc(i,j,k) = 0.0
          END IF
       END DO
       END DO
       END DO

       CALL update_boundary(visc)

     END SUBROUTINE les_fsf

     SUBROUTINE lfilter(in, out, halo)
       REAL(8), INTENT(IN)  :: in( -1:isize+2, -1:jsize+2, -1:ksize+2)
       REAL(8), INTENT(OUT) :: out(-1:isize+2, -1:jsize+2, -1:ksize+2)
       INTEGER, INTENT(IN)  :: halo

       INTEGER :: i, j, k
       INTEGER :: km0, km1, km2

!$OMP PARALLEL DO PRIVATE(km0, km1, km2)
       DO k=1-halo, ksize+halo
          km0 = maskindex(k-1)
          km1 = maskindex(k)
          km2 = maskindex(k+1)
       DO j=1-halo, jsize+halo
       DO i=1-halo, isize+halo
          out(i,j,k) = imask3d(i,j,km1)*(imask3d(i+1,j,km1)*(in(i+1,j,k)-in(i,j,k)) &
                                        +imask3d(i-1,j,km1)*(in(i-1,j,k)-in(i,j,k)) &
                                        +imask3d(i,j+1,km1)*(in(i,j+1,k)-in(i,j,k)) &
                                        +imask3d(i,j-1,km1)*(in(i,j-1,k)-in(i,j,k)) &
                                        +imask3d(i,j,  km2)*(in(i,j,k+1)-in(i,j,k)) &
                                        +imask3d(i,j,  km0)*(in(i,j,k-1)-in(i,j,k)))
       END DO
       END DO
       END DO

     END SUBROUTINE lfilter

   END SUBROUTINE les_viscosity

!-----------------------------------------------------------------------------------------------------------------------

   SUBROUTINE biharmonic_filter(gx, gy, gz)
    USE parameters, ONLY : eps_bhfilter

    REAL(8), INTENT(INOUT) :: gx(0:isize, 1:jsize, 1:ksize)
    REAL(8), INTENT(INOUT) :: gy(1:isize, 0:jsize, 1:ksize)
    REAL(8), INTENT(INOUT) :: gz(1:isize, 1:jsize, 0:ksize)

    REAL(8) :: l( -1:isize+1, -1:jsize+1)
    REAL(8) :: fx(-1:isize+2, -1:jsize+2)
    REAL(8) :: fy(-1:isize+2, -1:jsize+2)

    REAL(8) :: sqrtnu(-1:isize+2,-1:jsize+2)

    INTEGER :: i, j, k
    INTEGER :: kl, km

    IF (eps_bhfilter <= 0.0) RETURN

!$OMP PARALLEL DO PRIVATE(l, fx, fy, sqrtnu, kl, km)
    DO k= 1, ksize
       kl = dzindex(k)
       km = maskindex(k)

       DO j=-1, jsize+2
       DO i=-1, isize+2
          sqrtnu(i,j) = eps_bhfilter*(dx(i,j)**2 + dy(i,j)**2)/sqrt(dtime)
       END DO
       END DO

       DO j=-1, jsize+1
       DO i=-1, isize+2
          fx(i,j) = imask3d(i,j,km) * sqrtnu(i,j) * (u(i,j,k)-u(i-1,j,k))*idx0(i,j) * dy(i,j)*dz_star(i,j,kl)
       END DO
       END DO

       DO j=-1, jsize+1
       DO i=-1, isize+1
          fy(i,j) = imask3d(i,j,km)*imask3d(i+1,j,km)*imask3d(i,j+1,km)*imask3d(i+1,j+1,km)        &
               * 0.25*(sqrtnu(i,j)+sqrtnu(i+1,j)+sqrtnu(i,j+1)+sqrtnu(i+1,j+1))      &
               * (u(i,j+1,k)-u(i,j,k))*4.0/(dy(i,j)+dy(i+1,j)+dy(i,j+1)+dy(i+1,j+1)) &
               * 0.5*(dsy(i,j,kl)+dsy(i+1,j,kl))
       END DO
       END DO

       DO j= 0, jsize+1
       DO i=-1, isize+1
          l(i,j) = (fx(i+1,j)-fx(i,j)+fy(i,j)-fy(i,j-1)) * 2.0/(dvol(i,j,kl)+dvol(i+1,j,kl))
       END DO
       END DO

       DO j=1, jsize
       DO i=0, isize+1
          fx(i,j) = imask3d(i,j,km) * sqrtnu(i,j) * (l(i,j)-l(i-1,j))*idx0(i,j) * dy(i,j)*dz_star(i,j,kl)
       END DO
       END DO

       DO j=0, jsize
       DO i=0, isize
          fy(i,j) = imask3d(i,j,km)*imask3d(i+1,j,km)*imask3d(i,j+1,km)*imask3d(i+1,j+1,km)    &
               * 0.25*(sqrtnu(i,j)+sqrtnu(i+1,j)+sqrtnu(i,j+1)+sqrtnu(i+1,j+1))  &
               * (l(i,j+1)-l(i,j))*4.0/(dy(i,j)+dy(i+1,j)+dy(i,j+1)+dy(i+1,j+1)) &
               * 0.5*(dsy(i,j,kl)+dsy(i+1,j,kl))
       END DO
       END DO

       DO j=1, jsize
       DO i=0, isize
          gx(i,j,k) =  gx(i,j,k) - imask3d(i,j,km)*imask3d(i+1,j,km) &
               * (fx(i+1,j)-fx(i,j)+fy(i,j)-fy(i,j-1)) * 2.0/(dvol(i,j,kl)+dvol(i+1,j,kl))
       END DO
       END DO

!-----

       DO j=-1, jsize+1
       DO i=-1, isize+1
          fx(i,j) = imask3d(i,j,km)*imask3d(i,j+1,km)*imask3d(i+1,j,km)*imask3d(i+1,j+1,km)        &
               * 0.25*(sqrtnu(i,j)+sqrtnu(i,j+1)+sqrtnu(i+1,j)+sqrtnu(i+1,j+1))      &
               * (v(i+1,j,k)-v(i,j,k))*4.0/(dx(i,j)+dx(i,j+1)+dx(i+1,j)+dx(i+1,j+1)) &
               * 0.5*(dsx(i,j,kl)+dsx(i,j+1,kl))
       END DO
       END DO

       DO j=-1, jsize+2
       DO i=-1, isize+1
          fy(i,j) = imask3d(i,j,km) * sqrtnu(i,j) * (v(i,j,k)-v(i,j-1,k))*idy0(i,j) * dx(i,j)*dz_star(i,j,kl)
       END DO
       END DO

       DO j=-1, jsize+1
       DO i= 0, isize+1
          l(i,j) = (fx(i,j)-fx(i-1,j)+fy(i,j+1)-fy(i,j)) * 2.0/(dvol(i,j,kl)+dvol(i,j+1,kl))
       END DO
       END DO

       DO j=0, jsize
       DO i=0, isize
          fx(i,j) = imask3d(i,j,km)*imask3d(i,j+1,km)*imask3d(i+1,j,km)*imask3d(i+1,j+1,km)    &
               * 0.25*(sqrtnu(i,j)+sqrtnu(i,j+1)+sqrtnu(i+1,j)+sqrtnu(i+1,j+1))  &
               * (l(i+1,j)-l(i,j))*4.0/(dx(i,j)+dx(i,j+1)+dx(i+1,j)+dx(i+1,j+1)) &
               * 0.5*(dsx(i,j,kl)+dsx(i,j+1,kl))
       END DO
       END DO

       DO j=0, jsize+1
       DO i=1, isize
          fy(i,j) = imask3d(i,j,km) * sqrtnu(i,j) * (l(i,j)-l(i,j-1))*idy0(i,j) * dx(i,j)*dz_star(i,j,kl)
       END DO
       END DO

       DO j=0, jsize
       DO i=1, isize
          gy(i,j,k) =  gy(i,j,k) - imask3d(i,j,km)*imask3d(i,j+1,km) &
               * (fx(i,j)-fx(i-1,j)+fy(i,j+1)-fy(i,j)) * 2.0/(dvol(i,j,kl)+dvol(i,j+1,kl))
       END DO
       END DO
    END DO

  END SUBROUTINE biharmonic_filter

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE viscosity_diapycnal(u, v, w, nu1, nu2, n_x, n_y, n_z, gx, gy, gz)
  REAL(8), INTENT(IN) :: u( -slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8), INTENT(IN) :: v(1-slv:isize+slv,  -slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8), INTENT(IN) :: w(1-slv:isize+slv, 1-slv:jsize+slv,  -slv:ksize+slv)

  REAL(8), INTENT(IN) :: nu1(1-slv:isize+slv,1-slv:jsize+slv,1-slv:ksize+slv)
  REAL(8), INTENT(IN) :: nu2(1-slv:isize+slv,1-slv:jsize+slv,1-slv:ksize+slv)

  REAL(8), INTENT(IN) :: n_x( -slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8), INTENT(IN) :: n_y(1-slv:isize+slv,  -slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8), INTENT(IN) :: n_z(1-slv:isize+slv, 1-slv:jsize+slv,  -slv:ksize+slv)

  REAL(8), INTENT(INOUT) :: gx(0:isize, 1:jsize, 1:ksize)
  REAL(8), INTENT(INOUT) :: gy(1:isize, 0:jsize, 1:ksize)
  REAL(8), INTENT(INOUT) :: gz(1:isize, 1:jsize, 0:ksize)

  REAL(8) :: u_y(-1:isize+1, -1:jsize+1, 1:2)
  REAL(8) :: u_z(-1:isize+1,  0:jsize+1, 1:2)

  REAL(8) :: v_x(-1:isize+1, -1:jsize+1, 1:2)
  REAL(8) :: v_z( 0:isize+1, -1:jsize+1, 1:2)

  REAL(8) :: w_x(-1:isize+1,  0:jsize+1, 1:2)
  REAL(8) :: w_y( 0:isize+1, -1:jsize+1, 1:2)

  REAL(8) :: u_n(-1:isize+1,  0:jsize+1, 1:2)
  REAL(8) :: v_n( 0:isize+1, -1:jsize+1, 1:2)
  REAL(8) :: w_n( 0:isize+1,  0:jsize+1, 1:2)

  REAL(8) :: sn(0:isize+1, 0:jsize+1, 1:2)

  REAL(8) :: tau_xx(0:isize+1, 1:jsize)
  REAL(8) :: tau_yy(1:isize,   0:jsize+1)
  REAL(8) :: tau_zz(1:isize,   1:jsize, 1:2)

  REAL(8) :: tau_xy(0:isize, 0:jsize)
  REAL(8) :: tau_yz(1:isize, 0:jsize, 1:2)
  REAL(8) :: tau_zx(0:isize, 1:jsize, 1:2)

  REAL(4) :: sf_u, sf_l, sf_s

  INTEGER :: i, j, k
  INTEGER :: kstr, kend
  INTEGER :: ka, kb, kt
  INTEGER :: kl1, kl2
  INTEGER :: km1, km2, km3

  INTEGER :: tid

  tid  = 0
  kstr = 1
  kend = ksize

  sf_l = slipfactor_bottom
  sf_s = slipfactor_side

!$OMP PARALLEL PRIVATE(i, j, k, tid)  &
!$OMP          PRIVATE(u_y, u_z, u_n) &
!$OMP          PRIVATE(v_z, v_x, v_n) &
!$OMP          PRIVATE(w_x, w_y, w_n) &
!$OMP          PRIVATE(sn)            &
!$OMP          PRIVATE(tau_xx, tau_yy, tau_zz) &
!$OMP          PRIVATE(tau_xy, tau_yz, tau_zx) &
!$OMP          PRIVATE(kstr, kend) &
!$OMP          PRIVATE(ka, kb, kt) &
!$OMP          PRIVATE(kl1, kl2)   &
!$OMP          PRIVATE(sf_u)       &
!$OMP          PRIVATE(km1, km2, km3)

  ka = 1
  kb = 2

!$  tid = omp_get_thread_num()
!$  kstr=kstr_t(tid)
!$  kend=kend_t(tid)

  DO k=kstr-3, kend
     kt=ka
     ka=kb
     kb=kt

     km2 = maskindex(k+1)
     km3 = maskindex(k+2)

     IF (kcoord==kpes-1 .AND. k==ksize) THEN
        sf_u = slipfactor_surface
     ELSE
        sf_u = slipfactor_top
     END IF

     DO j= 0, jsize+1
     DO i=-1, isize+1
        w_x(i,j,kb) = ( (imask3d(i+1,j,km2)*imask3d(i+1,j,km3) - (1-imask3d(i,  j,km2)*imask3d(i,  j,km3))*sf_s) * w(i+1,j,k+1) &
                       -(imask3d(i,  j,km2)*imask3d(i,  j,km3) - (1-imask3d(i+1,j,km2)*imask3d(i+1,j,km3))*sf_s) * w(i,  j,k+1)) * idx1(i,j)
     END DO
     END DO

     DO j=-1, jsize+1
     DO i= 0, isize+1
        w_y(i,j,kb) = ( (imask3d(i,j+1,km2)*imask3d(i,j+1,km3) - (1-imask3d(i,j,  km2)*imask3d(i,j,  km3))*sf_s) * w(i,j+1,k+1) &
                       -(imask3d(i,j,  km2)*imask3d(i,j,  km3) - (1-imask3d(i,j+1,km2)*imask3d(i,j+1,km3))*sf_s) * w(i,j,  k+1)) * idy1(i,j)
     END DO
     END DO

     DO j= 0, jsize+1
     DO i= 0, isize+1
        w_n(i,j,kb) = (0.5D0*(((w(i,j,k+1)-w(i,j,k))*idz0(k+1))*dz(k+2)+((w(i,j,k+2)-w(i,j,k+1))*idz0(k+2))*dz(k+1))*n_z(i,j,k+1) &
             + 0.125D0*((n_x(i-1,j,k+1)+n_x(i,j,k+1))*dz(k+2)+(n_x(i-1,j,k+2)+n_x(i,j,k+2))*dz(k+1))*(w_x(i-1,j,kb)+w_x(i,j,kb)) &
             + 0.125D0*((n_y(i,j-1,k+1)+n_y(i,j,k+1))*dz(k+2)+(n_y(i,j-1,k+2)+n_y(i,j,k+2))*dz(k+1))*(w_y(i,j-1,kb)+w_y(i,j,kb))) * idz1(k+1)
     END DO
     END DO

     DO j= 0, jsize+1
     DO i=-1, isize+1
        u_z(i,j,kb) = ( (imask3d(i,j,km3)*imask3d(i+1,j,km3) - (1-imask3d(i,j,km2)*imask3d(i+1,j,km2))*sf_l) * u(i,j,k+2) &
                       -(imask3d(i,j,km2)*imask3d(i+1,j,km2) - (1-imask3d(i,j,km3)*imask3d(i+1,j,km3))*sf_u) * u(i,j,k+1)) * idz1(k+1)
     END DO
     END DO

     DO j=-1, jsize+1
     DO i= 0, isize+1
        v_z(i,j,kb) = ( (imask3d(i,j,km3)*imask3d(i,j+1,km3) - (1-imask3d(i,j,km2)*imask3d(i,j+1,km2))*sf_l) * v(i,j,k+2) &
                       -(imask3d(i,j,km2)*imask3d(i,j+1,km2) - (1-imask3d(i,j,km3)*imask3d(i,j+1,km3))*sf_u) * v(i,j,k+1)) * idz1(k+1)
     END DO
     END DO

     IF (k==kstr-3) CYCLE

     DO j=-1, jsize+1
     DO i=-1, isize+1
        u_y(i,j,kb) = ( (imask3d(i,j+1,km2)*imask3d(i+1,j+1,km2) - (1-imask3d(i,j,  km2)*imask3d(i+1,j,  km2))*sf_s) * u(i,j+1,k+1)  &
                       -(imask3d(i,j,  km2)*imask3d(i+1,j,  km2) - (1-imask3d(i,j+1,km2)*imask3d(i+1,j+1,km2))*sf_s) * u(i,j,  k+1)) &
                      * 4.0D0 / ((dy(i+1,j+1)*dx(i,j+1)+dy(i,j+1)*dx(i+1,j+1))*idx1(i,j+1) &
                                +(dy(i+1,j)  *dx(i,j)  +dy(i,j)  *dx(i+1,j)  )*idx1(i,j))
     END DO
     END DO

     DO j= 0, jsize+1
     DO i=-1, isize+1
        u_n(i,j,kb) = (0.5D0*(((u(i,j,k+1)-u(i-1,j,k+1))*idx0(i,j))*dx(i+1,j)+((u(i+1,j,k+1)-u(i,j,k+1))*idx0(i+1,j))*dx(i,j))*n_x(i,j,k+1) &
             + 0.125D0*((n_y(i,j-1,k+1)+n_y(i,j,k+1))*dx(i+1,j)+(n_y(i+1,j-1,k+1)+n_y(i+1,j,k+1))*dx(i,j))*(u_y(i,j-1,kb)+u_y(i,j,kb)) &
             + 0.125D0*((n_z(i,j,  k)  +n_z(i,j,k+1))*dx(i+1,j)+(n_z(i+1,j,  k)  +n_z(i+1,j,k+1))*dx(i,j))*(u_z(i,j,  ka)+u_z(i,j,kb))) * idx1(i,j)
     END DO
     END DO

     DO j=-1, jsize+1
     DO i=-1, isize+1
        v_x(i,j,kb) = ( (imask3d(i+1,j,km2)*imask3d(i+1,j+1,km2) - (1-imask3d(i,  j,km2)*imask3d(i,  j+1,km2))*sf_s) * v(i+1,j,k+1)  &
                       -(imask3d(i,  j,km2)*imask3d(i,  j+1,km2) - (1-imask3d(i+1,j,km2)*imask3d(i+1,j+1,km2))*sf_s) * v(i,  j,k+1)) &
                      * 4.0D0 / ((dx(i+1,j+1)*dy(i+1,j)+dx(i+1,j)*dy(i+1,j+1))*idy1(i+1,j) &
                                +(dx(i,  j+1)*dy(i,  j)+dx(i,  j)*dy(i,  j+1))*idy1(i,  j))
     END DO
     END DO

     DO j=-1, jsize+1
     DO i= 0, isize+1
        v_n(i,j,kb) = (0.5D0*(((v(i,j,k+1)-v(i,j-1,k+1))*idy0(i,j))*dy(i,j+1)+((v(i,j+1,k+1)-v(i,j,k+1))*idy0(i,j+1))*dy(i,j))*n_y(i,j,k+1) &
             + 0.125D0*((n_x(i-1,j,k+1)+n_x(i,j,k+1))*dy(i,j+1)+(n_x(i-1,j+1,k+1)+n_y(i,j+1,k+1))*dy(i,j))*(v_x(i-1,j,kb)+v_x(i,j,kb)) &
             + 0.125D0*((n_z(i,  j,k)  +n_z(i,j,k+1))*dy(i,j+1)+(n_z(i,  j+1,k)  +n_z(i,j+1,k+1))*dy(i,j))*(v_z(i,  j,ka)+v_z(i,j,kb))) * idy1(i,j)
     END DO
     END DO

     DO j=0, jsize+1
     DO i=0, isize+1
        sn(i,j,kb) = 0.5D0 * (u_n(i,j,kb)*n_x(i,j,k+1)+u_n(i-1,j,kb)*n_x(i-1,j,k+1) &
                            + v_n(i,j,kb)*n_y(i,j,k+1)+v_n(i,j-1,kb)*n_y(i,j-1,k+1) &
                            + w_n(i,j,kb)*n_z(i,j,k+1)+w_n(i,j,  ka)*n_z(i,j,k))
     END DO
     END DO

     DO j=1, jsize
     DO i=1, isize
        tau_zz(i,j,kb) = nu1(i,j,k+1) * 2.0D0*((w(i,j,k+1)-w(i,j,k))*idz0(k+1)) &
                       + nu2(i,j,k+1) * ((w_n(i,j,kb)*n_z(i,j,k+1)+w_n(i,j,ka)*n_z(i,j,k)) + 0.5D0*(n_z(i,j,k+1)**2+n_z(i,j,k)**2-2.0D0)*sn(i,j,kb))
     END DO
     END DO

     IF (k==kstr-2) CYCLE

     km1 = maskindex(k)
     kl1 = dzindex(k)
     kl2 = dzindex(k+1)

     DO j=0, jsize
     DO i=1, isize
        tau_yz(i,j,ka) = ((nu1(i,j,k)*dy(i,j)*dz(k)+nu1(i,j+1,k)*dy(i,j+1)*dz(k)+nu1(i,j,k+1)*dy(i,j)*dz(k+1)+nu1(i,j+1,k+1)*dy(i,j+1)*dz(k+1)) * (w_y(i,j,ka)+v_z(i,j,ka)) &
                        + (  nu2(i,j,  k)  *dy(i,j)  *dz(k)  *(w_n(i,j,  ka)*n_y(i,j,k)   + v_n(i,j,ka)*n_z(i,j,  k) + n_y(i,j,k)  *n_z(i,j,  k)*sn(i,j,  ka))              &
                           + nu2(i,j+1,k)  *dy(i,j+1)*dz(k)  *(w_n(i,j+1,ka)*n_y(i,j,k)   + v_n(i,j,ka)*n_z(i,j+1,k) + n_y(i,j,k)  *n_z(i,j+1,k)*sn(i,j+1,ka))              &
                           + nu2(i,j,  k+1)*dy(i,j)  *dz(k+1)*(w_n(i,j,  ka)*n_y(i,j,k+1) + v_n(i,j,kb)*n_z(i,j,  k) + n_y(i,j,k+1)*n_z(i,j,  k)*sn(i,j,  kb))              &
                           + nu2(i,j+1,k+1)*dy(i,j+1)*dz(k+1)*(w_n(i,j+1,ka)*n_y(i,j,k+1) + v_n(i,j,kb)*n_z(i,j+1,k) + n_y(i,j,k+1)*n_z(i,j+1,k)*sn(i,j+1,kb))))            &
                        / ((dy(i,j)+dy(i,j+1))*(dz(k)+dz(k+1)))
     END DO
     END DO

     DO j=1, jsize
     DO i=0, isize
        tau_zx(i,j,ka) = ((nu1(i,j,k)*dx(i,j)*dz(k)+nu1(i,j,k+1)*dx(i,j)*dz(k+1)+nu1(i+1,j,k)*dx(i+1,j)*dz(k)+nu1(i+1,j,k+1)*dx(i+1,j)*dz(k+1)) * (u_z(i,j,ka)+w_x(i,j,ka)) &
                        + (  nu2(i,  j,k)  *dx(i,  j)*dz(k)  *(u_n(i,j,ka)*n_z(i,  j,k) + w_n(i,  j,ka)*n_x(i,j,k)   + n_z(i,  j,k)*n_x(i,j,k)  *sn(i,  j,ka))              &
                           + nu2(i,  j,k+1)*dx(i,  j)*dz(k+1)*(u_n(i,j,kb)*n_z(i,  j,k) + w_n(i,  j,ka)*n_x(i,j,k+1) + n_z(i,  j,k)*n_x(i,j,k+1)*sn(i,  j,kb))              &
                           + nu2(i+1,j,k)  *dx(i+1,j)*dz(k)  *(u_n(i,j,ka)*n_z(i+1,j,k) + w_n(i+1,j,ka)*n_x(i,j,k)   + n_z(i+1,j,k)*n_x(i,j,k)  *sn(i+1,j,ka))              &
                           + nu2(i+1,j,k+1)*dx(i+1,j)*dz(k+1)*(u_n(i,j,kb)*n_z(i+1,j,k) + w_n(i+1,j,ka)*n_x(i,j,k+1) + n_z(i+1,j,k)*n_x(i,j,k+1)*sn(i+1,j,kb))))            &
                        / ((dx(i,j)+dx(i+1,j))*(dz(k)+dz(k+1)))
     END DO
     END DO

!$   IF (tid/=0 .AND. k==kstr-1) CYCLE

     DO j=1, jsize
     DO i=1, isize
        gz(i,j,k) = gz(i,j,k) + imask3d(i,j,km1)*imask3d(i,j,km2) &
             * (  tau_zx(i,j,ka)*(dsx(i,j,kl1)+dsx(i,j,kl2)) - tau_zx(i-1,j,  ka)*(dsx(i-1,j,  kl1)+dsx(i-1,j,  kl2))  &
                + tau_yz(i,j,ka)*(dsy(i,j,kl1)+dsy(i,j,kl2)) - tau_yz(i,  j-1,ka)*(dsy(i,  j-1,kl1)+dsy(i,  j-1,kl2))  &
                + tau_zz(i,j,kb)*(dsz(i,j)    +dsz(i,j))     - tau_zz(i,  j,  ka)*(dsz(i,j)        +dsz(i,j)))         &
             / (dvol(i,j,kl1)+dvol(i,j,kl2))
     END DO
     END DO

     IF (k==kstr-1) CYCLE

     DO j=1, jsize
     DO i=0, isize+1
        tau_xx(i,j) = nu1(i,j,k) * 2.0D0*((u(i,j,k)-u(i-1,j,k))*idx0(i,j)) &
                    + nu2(i,j,k) * ((u_n(i,j,ka)*n_x(i,j,k)+u_n(i-1,j,ka)*n_x(i-1,j,k)) + 0.5D0*(n_x(i,j,k)**2+n_x(i-1,j,k)**2-2.0D0)*sn(i,j,ka))
     END DO
     END DO

     DO j=0, jsize+1
     DO i=1, isize
        tau_yy(i,j) = nu1(i,j,k) * 2.0D0*((v(i,j,k)-v(i,j-1,k))*idy0(i,j)) &
                    + nu2(i,j,k) * ((v_n(i,j,ka)*n_y(i,j,k)+v_n(i,j-1,ka)*n_y(i,j-1,k)) + 0.5D0*(n_y(i,j,k)**2+n_y(i,j-1,k)**2-2.0D0)*sn(i,j,ka))
     END DO
     END DO

     DO j=0, jsize
     DO i=0, isize
        tau_xy(i,j) = ((nu1(i,j,k)*dsz(i,j)+nu1(i+1,j,k)*dsz(i+1,j)+nu1(i,j+1,k)*dsz(i,j+1)+nu1(i+1,j+1,k)*dsz(i+1,j+1)) * (v_x(i,j,ka)+u_y(i,j,ka))          &
                      + (  nu2(i,  j,  k)*dsz(i,  j)  *(v_n(i,  j,ka)*n_x(i,j,  k) + u_n(i,j,  ka)*n_y(i,  j,k) + n_x(i,j,  k)*n_y(i,  j,k)*sn(i,  j,  ka))   &
                         + nu2(i+1,j,  k)*dsz(i+1,j)  *(v_n(i+1,j,ka)*n_x(i,j,  k) + u_n(i,j,  ka)*n_y(i+1,j,k) + n_x(i,j,  k)*n_y(i+1,j,k)*sn(i+1,j,  ka))   &
                         + nu2(i,  j+1,k)*dsz(i,  j+1)*(v_n(i,  j,ka)*n_x(i,j+1,k) + u_n(i,j+1,ka)*n_y(i,  j,k) + n_x(i,j+1,k)*n_y(i,  j,k)*sn(i,  j+1,ka))   &
                         + nu2(i+1,j+1,k)*dsz(i+1,j+1)*(v_n(i+1,j,ka)*n_x(i,j+1,k) + u_n(i,j+1,ka)*n_y(i+1,j,k) + n_x(i,j+1,k)*n_y(i+1,j,k)*sn(i+1,j+1,ka)))) &
                      / (dsz(i,j)+dsz(i+1,j)+dsz(i,j+1)+dsz(i+1,j+1))
     END DO
     END DO

     DO j=1, jsize
     DO i=0, isize
        gx(i,j,k) = gx(i,j,k) + imask3d(i,j,km1)*imask3d(i+1,j,km1) &
             * (  tau_xx(i+1,j) *(dsx(i,j,kl1)+dsx(i+1,j,kl1)) - tau_xx(i,j)   *(dsx(i-1,j,  kl1)+dsx(i,  j,  kl1))  &
                + tau_xy(i,  j) *(dsy(i,j,kl1)+dsy(i+1,j,kl1)) - tau_xy(i,j-1) *(dsy(i,  j-1,kl1)+dsy(i+1,j-1,kl1))  &
                + tau_zx(i,j,ka)*(dsz(i,j)    +dsz(i+1,j))     - tau_zx(i,j,kb)*(dsz(i,j)        +dsz(i+1,j)))       &
             / (dvol(i,j,kl1)+dvol(i+1,j,kl1))
     END DO
     END DO

     DO j=0, jsize
     DO i=1, isize
        gy(i,j,k) = gy(i,j,k) + imask3d(i,j,km1)*imask3d(i,j+1,km1) &
             * (  tau_xy(i,j)   *(dsx(i,j,kl1)+dsx(i,j+1,kl1)) - tau_xy(i-1,j) *(dsx(i-1,j,  kl1)+dsx(i-1,j+1,kl1))  &
                + tau_yy(i,j+1) *(dsy(i,j,kl1)+dsy(i,j+1,kl1)) - tau_yy(i,  j) *(dsy(i,  j-1,kl1)+dsy(i,  j,  kl1))  &
                + tau_yz(i,j,ka)*(dsz(i,j)    +dsz(i,j+1))     - tau_yz(i,j,kb)*(dsz(i,j)        +dsz(i,j+1)))       &
             / (dvol(i,j,kl1)+dvol(i,j+1,kl1))
     END DO
     END DO
  END DO
!$OMP END PARALLEL

END SUBROUTINE viscosity_diapycnal

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE viscosity_nottuned(u, v, w, nu1, nu2, n_x, n_y, n_z, gx, gy, gz)
  REAL(8), INTENT(IN) :: u( -slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8), INTENT(IN) :: v(1-slv:isize+slv,  -slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8), INTENT(IN) :: w(1-slv:isize+slv, 1-slv:jsize+slv,  -slv:ksize+slv)

  REAL(8), INTENT(IN) :: nu1(1-slv:isize+slv,1-slv:jsize+slv,1-slv:ksize+slv)
  REAL(8), INTENT(IN) :: nu2(1-slv:isize+slv,1-slv:jsize+slv,1-slv:ksize+slv)

  REAL(8), INTENT(IN) :: n_x( -slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8), INTENT(IN) :: n_y(1-slv:isize+slv,  -slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8), INTENT(IN) :: n_z(1-slv:isize+slv, 1-slv:jsize+slv,  -slv:ksize+slv)

  REAL(8), INTENT(INOUT) :: gx(0:isize, 1:jsize, 1:ksize)
  REAL(8), INTENT(INOUT) :: gy(1:isize, 0:jsize, 1:ksize)
  REAL(8), INTENT(INOUT) :: gz(1:isize, 1:jsize, 0:ksize)

  REAL(8) :: u_x(-1:isize+2, -1:jsize+2, -1:ksize+2)
  REAL(8) :: u_y(-1:isize+1, -1:jsize+1,  0:ksize+1)
  REAL(8) :: u_z(-1:isize+1,  0:jsize+1, -1:ksize+1)

  REAL(8) :: v_x(-1:isize+1, -1:jsize+1,  0:ksize+1)
  REAL(8) :: v_y(-1:isize+2, -1:jsize+2, -1:ksize+2)
  REAL(8) :: v_z( 0:isize+1, -1:jsize+1, -1:ksize+1)

  REAL(8) :: w_x(-1:isize+1,  0:jsize+1, -1:ksize+1)
  REAL(8) :: w_y( 0:isize+1, -1:jsize+1, -1:ksize+1)
  REAL(8) :: w_z(-1:isize+2, -1:jsize+2, -1:ksize+2)

  REAL(8) :: u_n(-1:isize+1,  0:jsize+1,  0:ksize+1)
  REAL(8) :: v_n( 0:isize+1, -1:jsize+1,  0:ksize+1)
  REAL(8) :: w_n( 0:isize+1,  0:jsize+1, -1:ksize+1)

  REAL(8) :: sn(0:isize+1, 0:jsize+1, 0:ksize+1)

  REAL(8) :: tau_xx(0:isize+1, 1:jsize,   1:ksize)
  REAL(8) :: tau_yy(1:isize,   0:jsize+1, 1:ksize)
  REAL(8) :: tau_zz(1:isize,   1:jsize,   0:ksize+1)

  REAL(8) :: tau_xy(0:isize, 0:jsize, 1:ksize)
  REAL(8) :: tau_yz(1:isize, 0:jsize, 0:ksize)
  REAL(8) :: tau_zx(0:isize, 1:jsize, 0:ksize)

  REAL(4) :: sf_u, sf_l, sf_s

  INTEGER :: i, j, k

  sf_l = slipfactor_bottom
  sf_s = slipfactor_side

  DO k=-1, ksize+2
  DO j=-1, jsize+2
  DO i=-1, isize+2
     u_x(i,j,k) = (u(i,j,k)-u(i-1,j,k))*idx0(i,j)
     v_y(i,j,k) = (v(i,j,k)-v(i,j-1,k))*idy0(i,j)
     w_z(i,j,k) = (w(i,j,k)-w(i,j,k-1))*idz0(k)
  END DO
  END DO
  END DO

  DO k= 0, ksize+1
  DO j=-1, jsize+1
  DO i=-1, isize+1
     u_y(i,j,k) = ( (imask3d(i,j+1,k)*imask3d(i+1,j+1,k) - (1-imask3d(i,j,  k)*imask3d(i+1,j,  k))*sf_s) * u(i,j+1,k)  &
                   -(imask3d(i,j,  k)*imask3d(i+1,j,  k) - (1-imask3d(i,j+1,k)*imask3d(i+1,j+1,k))*sf_s) * u(i,j,  k)) &
                 * 4.0D0 / ((dy(i+1,j+1)*dx(i,j+1)+dy(i,j+1)*dx(i+1,j+1))*idx1(i,j+1) &
                           +(dy(i+1,j)  *dx(i,j)  +dy(i,j)  *dx(i+1,j)  )*idx1(i,j))
  END DO
  END DO
  END DO

  DO k= 0, ksize+1
  DO j=-1, jsize+1
  DO i=-1, isize+1
     v_x(i,j,k) = ( (imask3d(i+1,j,k)*imask3d(i+1,j+1,k) - (1-imask3d(i,  j,k)*imask3d(i,  j+1,k))*sf_s) * v(i+1,j,k)  &
                   -(imask3d(i,  j,k)*imask3d(i,  j+1,k) - (1-imask3d(i+1,j,k)*imask3d(i+1,j+1,k))*sf_s) * v(i,  j,k)) &
                 * 4.0D0 / ((dx(i+1,j+1)*dy(i+1,j)+dx(i+1,j)*dy(i+1,j+1))*idy1(i+1,j) &
                           +(dx(i,  j+1)*dy(i,  j)+dx(i,  j)*dy(i,  j+1))*idy1(i,  j))
  END DO
  END DO
  END DO

  DO k=-1, ksize+1

     IF (kcoord==kpes-1 .AND. k==ksize) THEN
        sf_u = slipfactor_surface
     ELSE
        sf_u = slipfactor_top
     END IF

     DO j= 0, jsize+1
     DO i=-1, isize+1

        u_z(i,j,k) = ( (imask3d(i,j,k+1)*imask3d(i+1,j,k+1) - (1-imask3d(i,j,k)  *imask3d(i+1,j,k)  )*sf_l) * u(i,j,k+1) &
                      -(imask3d(i,j,k)  *imask3d(i+1,j,k)   - (1-imask3d(i,j,k+1)*imask3d(i+1,j,k+1))*sf_u) * u(i,j,k)) * idz1(k)
     END DO
     END DO

     DO j=-1, jsize+1
     DO i= 0, isize+1
        v_z(i,j,k) = ( (imask3d(i,j,k+1)*imask3d(i,j+1,k+1) - (1-imask3d(i,j,k)  *imask3d(i,j+1,k)  )*sf_l)  * v(i,j,k+1) &
                      -(imask3d(i,j,k)  *imask3d(i,j+1,k)   - (1-imask3d(i,j,k+1)*imask3d(i,j+1,k+1))*sf_u) * v(i,j,k)) * idz1(k)
     END DO
     END DO
  END DO

  DO k=-1, ksize+1
     DO j= 0, jsize+1
     DO i=-1, isize+1
        w_x(i,j,k) = ( (imask3d(i+1,j,k)*imask3d(i+1,j,k+1) - (1-imask3d(i,  j,k)*imask3d(i,  j,k+1))*sf_s) * w(i+1,j,k) &
                      -(imask3d(i,  j,k)*imask3d(i,  j,k+1) - (1-imask3d(i+1,j,k)*imask3d(i+1,j,k+1))*sf_s) * w(i,  j,k)) * idx1(i,j)
     END DO
     END DO

     DO j=-1, jsize+1
     DO i= 0, isize+1
        w_y(i,j,k) = ( (imask3d(i,j+1,k)*imask3d(i,j+1,k+1) - (1-imask3d(i,j,  k)*imask3d(i,j,  k+1))*sf_s) * w(i,j+1,k) &
                      -(imask3d(i,j,  k)*imask3d(i,j,  k+1) - (1-imask3d(i,j+1,k)*imask3d(i,j+1,k+1))*sf_s) * w(i,j,  k)) * idy1(i,j)
     END DO
     END DO
  END DO

  DO k= 0, ksize+1
     DO j= 0, jsize+1
     DO i=-1, isize+1
        u_n(i,j,k) = (0.5D0*(u_x(i,j,k)*dx(i+1,j)+u_x(i+1,j,k)*dx(i,j))*n_x(i,j,k) &
             + 0.125D0*((n_y(i,j-1,k)+n_y(i,j,k))*dx(i+1,j)+(n_y(i+1,j-1,k)+n_y(i+1,j,k))*dx(i,j))*(u_y(i,j-1,k)+u_y(i,j,k)) &
             + 0.125D0*((n_z(i,j,k-1)+n_z(i,j,k))*dx(i+1,j)+(n_z(i+1,j,k-1)+n_z(i+1,j,k))*dx(i,j))*(u_z(i,j,k-1)+u_z(i,j,k))) * idx1(i,j)
     END DO
     END DO

     DO j=-1, jsize+1
     DO i= 0, isize+1
        v_n(i,j,k) = (0.5D0*(v_y(i,j,k)*dy(i,j+1)+v_y(i,j+1,k)*dy(i,j))*n_y(i,j,k) &
             + 0.125D0*((n_x(i-1,j,k)+n_x(i,j,k))*dy(i,j+1)+(n_x(i-1,j+1,k)+n_y(i,j+1,k))*dy(i,j))*(v_x(i-1,j,k)+v_x(i,j,k)) &
             + 0.125D0*((n_z(i,j,k-1)+n_z(i,j,k))*dy(i,j+1)+(n_z(i,j+1,k-1)+n_z(i,j+1,k))*dy(i,j))*(v_z(i,j,k-1)+v_z(i,j,k))) * idy1(i,j)
     END DO
     END DO
  END DO

  DO k=-1, ksize+1
  DO j= 0, jsize+1
  DO i= 0, isize+1
     w_n(i,j,k) = (0.5D0*(w_z(i,j,k)*dz(k+1)+w_z(i,j,k+1)*dz(k))*n_z(i,j,k) &
          + 0.125D0*((n_x(i-1,j,k)+n_x(i,j,k))*dz(k+1)+(n_x(i-1,j,k+1)+n_x(i,j,k+1))*dz(k))*(w_x(i-1,j,k)+w_x(i,j,k)) &
          + 0.125D0*((n_y(i,j-1,k)+n_y(i,j,k))*dz(k+1)+(n_y(i,j-1,k+1)+n_y(i,j,k+1))*dz(k))*(w_y(i,j-1,k)+w_y(i,j,k))) * idz1(k)
  END DO
  END DO
  END DO

  DO k=0, ksize+1
  DO j=0, jsize+1
  DO i=0, isize+1
     sn(i,j,k) = 0.5D0 * (u_n(i,j,k)*n_x(i,j,k)+u_n(i-1,j,k)*n_x(i-1,j,k) &
                        + v_n(i,j,k)*n_y(i,j,k)+v_n(i,j-1,k)*n_y(i,j-1,k) &
                        + w_n(i,j,k)*n_z(i,j,k)+w_n(i,j,k-1)*n_z(i,j,k-1))
  END DO
  END DO
  END DO

  DO k=1, ksize
  DO j=1, jsize
  DO i=0, isize+1
     tau_xx(i,j,k) = nu1(i,j,k) * 2.0D0*u_x(i,j,k) &
                   + nu2(i,j,k) * ((u_n(i,j,k)*n_x(i,j,k)+u_n(i-1,j,k)*n_x(i-1,j,k)) + 0.5D0*(n_x(i,j,k)**2+n_x(i-1,j,k)**2-2.0D0)*sn(i,j,k))
  END DO
  END DO
  END DO

  DO k=1, ksize
  DO j=0, jsize+1
  DO i=1, isize
     tau_yy(i,j,k) = nu1(i,j,k) * 2.0D0*v_y(i,j,k) &
                   + nu2(i,j,k) * ((v_n(i,j,k)*n_y(i,j,k)+v_n(i,j-1,k)*n_y(i,j-1,k)) + 0.5D0*(n_y(i,j,k)**2+n_y(i,j-1,k)**2-2.0D0)*sn(i,j,k))
  END DO
  END DO
  END DO

  DO k=0, ksize+1
  DO j=1, jsize
  DO i=1, isize
     tau_zz(i,j,k) = nu1(i,j,k) * 2.0D0*w_z(i,j,k) &
                   + nu2(i,j,k) * ((w_n(i,j,k)*n_z(i,j,k)+w_n(i,j,k-1)*n_z(i,j,k-1)) + 0.5D0*(n_z(i,j,k)**2+n_z(i,j,k-1)**2-2.0D0)*sn(i,j,k))
  END DO
  END DO
  END DO

  DO k=1, ksize
  DO j=0, jsize
  DO i=0, isize
     tau_xy(i,j,k) = ((nu1(i,j,k)*dsz(i,j)+nu1(i+1,j,k)*dsz(i+1,j)+nu1(i,j+1,k)*dsz(i,j+1)+nu1(i+1,j+1,k)*dsz(i+1,j+1)) * (v_x(i,j,k)+u_y(i,j,k))          &
                    + (  nu2(i,  j,  k)*dsz(i,  j)  *(v_n(i,  j,k)*n_x(i,j,  k) + u_n(i,j,  k)*n_y(i,  j,k) + n_x(i,j,  k)*n_y(i,  j,k)*sn(i,  j,  k))     &
                       + nu2(i+1,j,  k)*dsz(i+1,j)  *(v_n(i+1,j,k)*n_x(i,j,  k) + u_n(i,j,  k)*n_y(i+1,j,k) + n_x(i,j,  k)*n_y(i+1,j,k)*sn(i+1,j,  k))     &
                       + nu2(i,  j+1,k)*dsz(i,  j+1)*(v_n(i,  j,k)*n_x(i,j+1,k) + u_n(i,j+1,k)*n_y(i,  j,k) + n_x(i,j+1,k)*n_y(i,  j,k)*sn(i,  j+1,k))     &
                       + nu2(i+1,j+1,k)*dsz(i+1,j+1)*(v_n(i+1,j,k)*n_x(i,j+1,k) + u_n(i,j+1,k)*n_y(i+1,j,k) + n_x(i,j+1,k)*n_y(i+1,j,k)*sn(i+1,j+1,k))))   &
                    / (dsz(i,j)+dsz(i+1,j)+dsz(i,j+1)+dsz(i+1,j+1))
  END DO
  END DO
  END DO

  DO k=0, ksize
  DO j=0, jsize
  DO i=1, isize
     tau_yz(i,j,k) = ((nu1(i,j,k)*dy(i,j)*dz(k)+nu1(i,j+1,k)*dy(i,j+1)*dz(k)+nu1(i,j,k+1)*dy(i,j)*dz(k+1)+nu1(i,j+1,k+1)*dy(i,j+1)*dz(k+1)) * (w_y(i,j,k)+v_z(i,j,k)) &
                    + (  nu2(i,j,  k)  *dy(i,j)  *dz(k)  *(w_n(i,j,  k)*n_y(i,j,k)   + v_n(i,j,k)  *n_z(i,j,  k) + n_y(i,j,k)  *n_z(i,j,  k)*sn(i,j,  k))             &
                       + nu2(i,j+1,k)  *dy(i,j+1)*dz(k)  *(w_n(i,j+1,k)*n_y(i,j,k)   + v_n(i,j,k)  *n_z(i,j+1,k) + n_y(i,j,k)  *n_z(i,j+1,k)*sn(i,j+1,k))             &
                       + nu2(i,j,  k+1)*dy(i,j)  *dz(k+1)*(w_n(i,j,  k)*n_y(i,j,k+1) + v_n(i,j,k+1)*n_z(i,j,  k) + n_y(i,j,k+1)*n_z(i,j,  k)*sn(i,j,  k+1))           &
                       + nu2(i,j+1,k+1)*dy(i,j+1)*dz(k+1)*(w_n(i,j+1,k)*n_y(i,j,k+1) + v_n(i,j,k+1)*n_z(i,j+1,k) + n_y(i,j,k+1)*n_z(i,j+1,k)*sn(i,j+1,k+1))))         &
                    / ((dy(i,j)+dy(i,j+1))*(dz(k)+dz(k+1)))
  END DO
  END DO
  END DO

  DO k=0, ksize
  DO j=1, jsize
  DO i=0, isize
     tau_zx(i,j,k) = ((nu1(i,j,k)*dx(i,j)*dz(k)+nu1(i,j,k+1)*dx(i,j)*dz(k+1)+nu1(i+1,j,k)*dx(i+1,j)*dz(k)+nu1(i+1,j,k+1)*dx(i+1,j)*dz(k+1)) * (u_z(i,j,k)+w_x(i,j,k)) &
                    + (  nu2(i,  j,k)  *dx(i,  j)*dz(k)  *(u_n(i,j,k)  *n_z(i,  j,k) + w_n(i,  j,k)*n_x(i,j,k)   + n_z(i,  j,k)*n_x(i,j,k)  *sn(i,  j,k))             &
                       + nu2(i,  j,k+1)*dx(i,  j)*dz(k+1)*(u_n(i,j,k+1)*n_z(i,  j,k) + w_n(i,  j,k)*n_x(i,j,k+1) + n_z(i,  j,k)*n_x(i,j,k+1)*sn(i,  j,k+1))           &
                       + nu2(i+1,j,k)  *dx(i+1,j)*dz(k)  *(u_n(i,j,k)  *n_z(i+1,j,k) + w_n(i+1,j,k)*n_x(i,j,k)   + n_z(i+1,j,k)*n_x(i,j,k)  *sn(i+1,j,k))             &
                       + nu2(i+1,j,k+1)*dx(i+1,j)*dz(k+1)*(u_n(i,j,k+1)*n_z(i+1,j,k) + w_n(i+1,j,k)*n_x(i,j,k+1) + n_z(i+1,j,k)*n_x(i,j,k+1)*sn(i+1,j,k+1))))         &
                    / ((dx(i,j)+dx(i+1,j))*(dz(k)+dz(k+1)))
  END DO
  END DO
  END DO

  DO k=1, ksize
  DO j=1, jsize
  DO i=0, isize
     gx(i,j,k) = gx(i,j,k) + imask3d(i,j,k)*imask3d(i+1,j,k) &
          * (  tau_xx(i+1,j,k)*(dsx(i,j,k)+dsx(i+1,j,k)) - tau_xx(i,j,  k)  *(dsx(i-1,j,  k)+dsx(i,  j,  k))  &
             + tau_xy(i,  j,k)*(dsy(i,j,k)+dsy(i+1,j,k)) - tau_xy(i,j-1,k)  *(dsy(i,  j-1,k)+dsy(i+1,j-1,k))  &
             + tau_zx(i,  j,k)*(dsz(i,j)  +dsz(i+1,j)  ) - tau_zx(i,j,  k-1)*(dsz(i,  j)    +dsz(i+1,j)))     &
          / (dvol(i,j,k)+dvol(i+1,j,k))
  END DO
  END DO
  END DO

  DO k=1, ksize
  DO j=0, jsize
  DO i=1, isize
     gy(i,j,k) = gy(i,j,k) + imask3d(i,j,k)*imask3d(i,j+1,k) &
          * (  tau_xy(i,j,  k)*(dsx(i,j,k)+dsx(i,j+1,k)) - tau_xy(i-1,j,k)  *(dsx(i-1,j,  k)+dsx(i-1,j+1,k))  &
             + tau_yy(i,j+1,k)*(dsy(i,j,k)+dsy(i,j+1,k)) - tau_yy(i,  j,k)  *(dsy(i,  j-1,k)+dsy(i,  j,  k))  &
             + tau_yz(i,j,  k)*(dsz(i,j)  +dsz(i,j+1)  ) - tau_yz(i,  j,k-1)*(dsz(i,  j)    +dsz(i,  j+1)  )) &
          / (dvol(i,j,k)+dvol(i,j+1,k))
  END DO
  END DO
  END DO

  DO k=0, ksize
  DO j=1, jsize
  DO i=1, isize
     gz(i,j,k) = gz(i,j,k) + imask3d(i,j,k)*imask3d(i,j,k+1) &
          * (  tau_zx(i,j,k)  *(dsx(i,j,k)+dsx(i,j,k+1)) - tau_zx(i-1,j,  k)*(dsx(i-1,j,  k)+dsx(i-1,j,  k+1))  &
             + tau_yz(i,j,k)  *(dsy(i,j,k)+dsy(i,j,k+1)) - tau_yz(i,  j-1,k)*(dsy(i,  j-1,k)+dsy(i,  j-1,k+1))  &
             + tau_zz(i,j,k+1)*(dsz(i,j)  +dsz(i,j)  )   - tau_zz(i,  j,  k)*(dsz(i,  j)    +dsz(i,  j)      )) &
          / (dvol(i,j,k)+dvol(i,j,k+1))
  END DO
  END DO
  END DO


END SUBROUTINE viscosity_nottuned

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE diapycnal_vector(t, s, n_x, n_y, n_z)
  REAL(8), INTENT(IN)  :: t(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8), INTENT(IN)  :: s(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)

  REAL(8), INTENT(OUT) :: n_x( -slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8), INTENT(OUT) :: n_y(1-slv:isize+slv,  -slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8), INTENT(OUT) :: n_z(1-slv:isize+slv, 1-slv:jsize+slv,  -slv:ksize+slv)

  REAL(8) :: sigma(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)

  REAL(8) :: tmpx( -slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: tmpy(1-slv:isize+slv,  -slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: tmpz(1-slv:isize+slv, 1-slv:jsize+slv,  -slv:ksize+slv)

  REAL(8) :: absn(-1:isize+1, -1:jsize+1, -1:ksize+1)

  REAL(8), PARAMETER :: eps = 1.0D-6

  INTEGER :: i, j, k

  tmpx = 0.0D0
  tmpy = 0.0D0
  tmpz = 0.0D0

  DO k= 0, ksize+1
  DO j= 0, jsize+1
  DO i=-1, isize+1
     tmpx(i,j,k) = - imask3d(i,j,k)*imask3d(i+1,j,k) * (sigma(i+1,j,k) - sigma(i,j,k))*idx1(i,j)
  END DO
  END DO
  END DO

  DO k= 0, ksize+1
  DO j=-1, jsize+1
  DO i= 0, isize+1
     tmpy(i,j,k) = - imask3d(i,j,k)*imask3d(i,j+1,k) * (sigma(i,j+1,k) - sigma(i,j,k))*idy1(i,j)
  END DO
  END DO
  END DO

  DO k=-1, ksize+1
  DO j= 0, jsize+1
  DO i= 0, isize+1
     tmpz(i,j,k) = - imask3d(i,j,k)*imask3d(i,j,k+1) * (sigma(i,j,k+1) - sigma(i,j,k))*idz1(k)
  END DO
  END DO
  END DO

  n_x(:,:,:) = 0.0D0
  n_y(:,:,:) = 0.0D0
  n_z(:,:,:) = 1.0D0

  DO k= 0, ksize+1
  DO j= 0, jsize+1
  DO i=-1, isize+1
     absn(i,j,k) = sqrt(tmpx(i,j,k)**2 &
          + ((tmpy(i,j-1,k)+tmpy(i+1,j-1,k)+tmpy(i,j,k)+tmpy(i+1,j,k))*0.25D0)**2 &
          + ((tmpz(i,j,k-1)+tmpz(i+1,j,k-1)+tmpz(i,j,k)+tmpz(i+1,j,k))*0.25D0)**2)
  END DO
  END DO
  END DO

  DO k= 0, ksize+1
  DO j= 0, jsize+1
  DO i=-1, isize+1
     IF (absn(i,j,k) > eps)  n_x(i,j,k) = tmpx(i,j,k)/absn(i,j,k)
  END DO
  END DO
  END DO

  DO k= 0, ksize+1
  DO j=-1, jsize+1
  DO i= 0, isize+1
     absn(i,j,k) = sqrt(tmpy(i,j,k)**2 &
          + ((tmpz(i,j,k-1)+tmpz(i,j+1,k-1)+tmpz(i,j,k)+tmpz(i,j+1,k))*0.25D0)**2 &
          + ((tmpx(i-1,j,k)+tmpx(i-1,j+1,k)+tmpx(i,j,k)+tmpx(i,j+1,k))*0.25D0)**2)
  END DO
  END DO
  END DO

  DO k= 0, ksize+1
  DO j=-1, jsize+1
  DO i= 0, isize+1
     IF (absn(i,j,k) > eps) n_y(i,j,k) = tmpy(i,j,k)/absn(i,j,k)
  END DO
  END DO
  END DO

  DO k=-1, ksize+1
  DO j= 0, jsize+1
  DO i= 0, isize+1
     absn(i,j,k) = sqrt(tmpz(i,j,k)**2 &
          + ((tmpx(i-1,j,k)+tmpx(i-1,j,k+1)+tmpx(i,j,k)+tmpx(i,j,k+1))*0.25D0)**2 &
          + ((tmpy(i,j-1,k)+tmpy(i,j-1,k+1)+tmpy(i,j,k)+tmpy(i,j,k+1))*0.25D0)**2)
  END DO
  END DO
  END DO

  DO k=-1, ksize+1
  DO j= 0, jsize+1
  DO i= 0, isize+1
     IF (absn(i,j,k) > eps) n_z(i,j,k) = tmpz(i,j,k)/absn(i,j,k)
  END DO
  END DO
  END DO

  CALL update_boundary(n_x)
  CALL update_boundary(n_y)
  CALL update_boundary(n_z)

END SUBROUTINE diapycnal_vector

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE contraction_vector(ax, ay, az, bx, by, bz, c)
  REAL(8), INTENT(IN) :: ax( -slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8), INTENT(IN) :: ay(1-slv:isize+slv,  -slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8), INTENT(IN) :: az(1-slv:isize+slv, 1-slv:jsize+slv,  -slv:ksize+slv)

  REAL(8), INTENT(IN) :: bx( -slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8), INTENT(IN) :: by(1-slv:isize+slv,  -slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8), INTENT(IN) :: bz(1-slv:isize+slv, 1-slv:jsize+slv,  -slv:ksize+slv)

  REAL(8), INTENT(OUT) :: c(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)

  INTEGER :: i, j, k

  DO k=1, ksize
  DO j=1, jsize
  DO i=1, isize
     c(i,j,k) = ( (ax(i-1,j,k)+ax(i,j,k))*(bx(i-1,j,k)+bx(i,j,k)) + ax(i-1,j,k)*bx(i-1,j,k) + ax(i,j,k)*bx(i,j,k) &
                 +(ay(i,j-1,k)+ay(i,j,k))*(by(i,j-1,k)+by(i,j,k)) + ay(i,j-1,k)*by(i,j-1,k) + ay(i,j,k)*by(i,j,k) &
                 +(az(i,j,k-1)+az(i,j,k))*(bz(i,j,k-1)+bz(i,j,k)) + az(i,j,k-1)*bz(i,j,k-1) + az(i,j,k)*bz(i,j,k)) / 6.0D0
  END DO
  END DO
  END DO

  CALL update_boundary(c)

END SUBROUTINE contraction_vector

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE contraction_tensor(axx, ayy, azz, axy, ayz, azx, bxx, byy, bzz, bxy, byz, bzx, c)
  REAL(8), INTENT(IN) :: axx(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8), INTENT(IN) :: ayy(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8), INTENT(IN) :: azz(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8), INTENT(IN) :: axy( -slv:isize+slv,  -slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8), INTENT(IN) :: ayz(1-slv:isize+slv,  -slv:jsize+slv,  -slv:ksize+slv)
  REAL(8), INTENT(IN) :: azx( -slv:isize+slv, 1-slv:jsize+slv,  -slv:ksize+slv)

  REAL(8), INTENT(IN) :: bxx(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8), INTENT(IN) :: byy(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8), INTENT(IN) :: bzz(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8), INTENT(IN) :: bxy( -slv:isize+slv,  -slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8), INTENT(IN) :: byz(1-slv:isize+slv,  -slv:jsize+slv,  -slv:ksize+slv)
  REAL(8), INTENT(IN) :: bzx( -slv:isize+slv, 1-slv:jsize+slv,  -slv:ksize+slv)

  REAL(8), INTENT(OUT) :: c(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)

  INTEGER :: i, j, k

  DO k=1, ksize
  DO j=1, jsize
  DO i=1, isize
     c(i,j,k) = axx(i,j,k)*bxx(i,j,k) + ayy(i,j,k)*byy(i,j,k) + azz(i,j,k)*bzz(i,j,k) &
          + ( 4.0D0*( axy(i-1,j-1,k)*bxy(i-1,j-1,k)+axy(i-1,j,k)*bxy(i-1,j,k)+axy(i,j-1,k)*bxy(i,j-1,k) + axy(i,j,k)*bxy(i,j,k)  &
                     +ayz(i,j-1,k-1)*byz(i,j-1,k-1)+ayz(i,j-1,k)*byz(i,j-1,k)+ayz(i,j,k-1)*byz(i,j,k-1) + ayz(i,j,k)*byz(i,j,k)  &
                     +azx(i-1,j,k-1)*bzx(i-1,j,k-1)+azx(i,j,k-1)*bzx(i,j,k-1)+azx(i-1,j,k)*bzx(i-1,j,k) + azx(i,j,k)*bzx(i,j,k)) &
             +2.0D0*( (axy(i-1,j-1,k)+axy(i,j,k))*(bxy(i-1,j,k)+bxy(i,j-1,k))+(axy(i-1,j,k)+axy(i,j-1,k))*(bxy(i-1,j-1,k)+bxy(i,j,k))  &
                     +(ayz(i,j-1,k-1)+ayz(i,j,k))*(byz(i,j-1,k)+byz(i,j,k-1))+(ayz(i,j-1,k)+ayz(i,j,k-1))*(byz(i,j-1,k-1)+byz(i,j,k))  &
                     +(azx(i-1,j,k-1)+azx(i,j,k))*(bzx(i,j,k-1)+bzx(i-1,j,k))+(azx(i,j,k-1)+azx(i-1,j,k))*(bzx(i-1,j,k-1)+bzx(i,j,k))) &
             + axy(i-1,j-1,k)*bxy(i,j,k)+axy(i-1,j,k)*bxy(i,j-1,k)+axy(i,j-1,k)*bxy(i-1,j,k)+axy(i,j,k)*bxy(i-1,j-1,k) &
             + ayz(i,j-1,k-1)*byz(i,j,k)+ayz(i,j-1,k)*byz(i,j,k-1)+ayz(i,j,k-1)*byz(i,j-1,k)+ayz(i,j,k)*byz(i,j-1,k-1) &
             + azx(i-1,j,k-1)*bzx(i,j,k)+azx(i,j,k-1)*bzx(i-1,j,k)+azx(i-1,j,k)*bzx(i,j,k-1)+azx(i,j,k)*bzx(i-1,j,k-1) ) / 18.0D0
  END DO
  END DO
  END DO

  CALL update_boundary(c)

END SUBROUTINE contraction_tensor


END MODULE subgrid
