#include "macro.h"

MODULE dles
  USE misc
  USE geometry
  USE subgrid
  IMPLICIT NONE

  REAL(8), ALLOCATABLE :: u_test(:,:,:)
  REAL(8), ALLOCATABLE :: v_test(:,:,:)
  REAL(8), ALLOCATABLE :: w_test(:,:,:)

  REAL(8), ALLOCATABLE :: tvol(:,:,:)
  REAL(8), ALLOCATABLE :: tvol_x(:,:,:)
  REAL(8), ALLOCATABLE :: tvol_y(:,:,:)
  REAL(8), ALLOCATABLE :: tvol_z(:,:,:)
  REAL(8), ALLOCATABLE :: tvol_xy(:,:,:)
  REAL(8), ALLOCATABLE :: tvol_yz(:,:,:)
  REAL(8), ALLOCATABLE :: tvol_zx(:,:,:)

  REAL(8), ALLOCATABLE :: nfreq2_test(:,:,:)
  REAL(8), ALLOCATABLE :: srate_test(:,:,:)

  REAL(8), ALLOCATABLE :: dles_flag(:,:,:)

  REAL(8), PARAMETER :: c_max = 1.0D2

  REAL(8), PARAMETER :: tg_ratio2 = 4.0D0

  REAL(8), ALLOCATABLE :: dles_tmp1(:,:,:)
  REAL(8), ALLOCATABLE :: dles_tmp2(:,:,:)

  REAL(8), PARAMETER :: dles_gamma = 5.0D-2

  REAL(8), PARAMETER :: eps = 1.0D-8

CONTAINS

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE test_filter(in, out)
  REAL(8), INTENT(IN)  :: in( 1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8), INTENT(OUT) :: out(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)

  REAL(8) :: tmp(0:isize+1, 0:jsize+1, 0:ksize+1)

  INTEGER :: i, j, k

  DO k=0, ksize+1
  DO j=0, jsize+1
  DO i=0, isize+1
     tmp(i,j,k) = imask3d(i,j,k) * in(i,j,k) * dvol(i,j,k)
  END DO
  END DO
  END DO

  DO k=1, ksize
  DO j=1, jsize
  DO i=1, isize
     IF (tvol(i,j,k) == 0.0D0) THEN
        out(i,j,k) = 0.0D0
     ELSE
        out(i,j,k) = (tmp(i,j,k) &
             + 0.5   * (tmp(i+1,j,  k)   + tmp(i-1,j,  k)    &
                      + tmp(i,  j+1,k)   + tmp(i,  j-1,k)    &
                      + tmp(i,  j,  k+1) + tmp(i,  j,  k-1)) &
             + 0.25  * (tmp(i+1,j+1,k)   + tmp(i-1,j+1,k)   + tmp(i+1,j-1,k)   + tmp(i-1,j-1,k)     &
                      + tmp(i  ,j+1,k+1) + tmp(i  ,j-1,k+1) + tmp(i  ,j+1,k-1) + tmp(i  ,j-1,k-1)   &
                      + tmp(i+1,j,  k+1) + tmp(i+1,j,  k-1) + tmp(i-1,j,  k+1) + tmp(i-1,j,  k-1))  &
             + 0.125 * (tmp(i+1,j+1,k+1) + tmp(i-1,j+1,k+1) + tmp(i+1,j-1,k+1) + tmp(i-1,j-1,k+1)   &
                      + tmp(i+1,j+1,k-1) + tmp(i-1,j+1,k-1) + tmp(i+1,j-1,k-1) + tmp(i-1,j-1,k-1))) / tvol(i,j,k)
     END IF
  END DO
  END DO
  END DO

END SUBROUTINE test_filter

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE test_filter_x(in, out)
  REAL(8), INTENT(IN)  :: in( -slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8), INTENT(OUT) :: out(-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)

  REAL(8) :: tmp(0:isize+1, 0:jsize+1, 0:ksize+1)

  INTEGER :: i, j, k

  DO k=0, ksize+1
  DO j=0, jsize+1
  DO i=0, isize+1
     tmp(i,j,k) = imask3d(i,j,k) * 0.5*(in(i-1,j,k)+in(i,j,k)) * dvol(i,j,k)
  END DO
  END DO
  END DO

  DO k=1, ksize
  DO j=1, jsize
  DO i=0, isize
     IF (tvol_x(i,j,k) == 0.0) THEN
        out(i,j,k) = 0.0
     ELSE
        out(i,j,k) = (tmp(i,j,k) + tmp(i+1,j,k) &
             + 0.5 * (tmp(i  ,j+1,k)   + tmp(i,  j-1,k)   + tmp(i  ,j,  k+1) + tmp(i,  j,  k-1)  &
                    + tmp(i+1,j+1,k)   + tmp(i+1,j-1,k)   + tmp(i+1,j,  k+1) + tmp(i+1,j,  k-1)) &
             + 0.25 *(tmp(i,  j+1,k+1) + tmp(i,  j-1,k+1) + tmp(i,  j+1,k-1) + tmp(i,  j-1,k-1)  &
                    + tmp(i+1,j+1,k+1) + tmp(i+1,j-1,k+1) + tmp(i+1,j+1,k-1) + tmp(i+1,j-1,k-1))) / tvol_x(i,j,k)
     END IF
  END DO
  END DO
  END DO

END SUBROUTINE test_filter_x

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE test_filter_y(in, out)
  REAL(8), INTENT(IN)  :: in( 1-slv:isize+slv, -slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8), INTENT(OUT) :: out(1-slv:isize+slv, -slv:jsize+slv, 1-slv:ksize+slv)

  REAL(8) :: tmp(0:isize+1, 0:jsize+1, 0:ksize+1)

  INTEGER :: i, j, k

  DO k=0, ksize+1
  DO j=0, jsize+1
  DO i=0, isize+1
     tmp(i,j,k) = imask3d(i,j,k) * 0.5*(in(i,j-1,k)+in(i,j,k)) * dvol(i,j,k)
  END DO
  END DO
  END DO

  DO k=1, ksize
  DO j=0, jsize
  DO i=1, isize
     IF (tvol_y(i,j,k) == 0.0) THEN
        out(i,j,k) = 0.0
     ELSE
        out(i,j,k) = (tmp(i,j,k) + tmp(i,j+1,k) &
             + 0.5 * (tmp(i,  j,  k+1) + tmp(i,  j,  k-1) + tmp(i+1,j,  k)   + tmp(i-1,j,  k)    &
                    + tmp(i,  j+1,k+1) + tmp(i,  j+1,k-1) + tmp(i+1,j+1,k)   + tmp(i-1,j+1,k))   &
             + 0.25 *(tmp(i+1,j,  k+1) + tmp(i+1,j,  k-1) + tmp(i-1,j,  k+1) + tmp(i-1,j,  k-1)  &
                    + tmp(i+1,j+1,k+1) + tmp(i+1,j+1,k-1) + tmp(i-1,j+1,k+1) + tmp(i-1,j+1,k-1))) / tvol_y(i,j,k)
     END IF
  END DO
  END DO
  END DO

END SUBROUTINE test_filter_y

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE test_filter_z(in, out)
  REAL(8), INTENT(IN)  :: in( 1-slv:isize+slv, 1-slv:jsize+slv, -slv:ksize+slv)
  REAL(8), INTENT(OUT) :: out(1-slv:isize+slv, 1-slv:jsize+slv, -slv:ksize+slv)

  REAL(8) :: tmp(0:isize+1, 0:jsize+1, 0:ksize+1)
  REAL(8) :: tvol

  INTEGER :: i, j, k

  DO k=0, ksize+1
  DO j=0, jsize+1
  DO i=0, isize+1
     tmp(i,j,k) = imask3d(i,j,k) * 0.5*(in(i,j,k-1)+in(i,j,k)) * dvol(i,j,k)
  END DO
  END DO
  END DO

  DO k=0, ksize
  DO j=1, jsize
  DO i=1, isize
     IF (tvol_z(i,j,k) == 0.0) THEN
        out(i,j,k) = 0.0
     ELSE
        out(i,j,k) = (tmp(i,j,k) + tmp(i,j,k+1) &
             + 0.5 * (tmp(i+1,j,  k)   + tmp(i-1,j,  k)   + tmp(i,  j+1,k)   + tmp(i,  j-1,k)    &
                    + tmp(i+1,j,  k+1) + tmp(i-1,j,  k+1) + tmp(i,  j+1,k+1) + tmp(i,  j-1,k+1)) &
             + 0.25 *(tmp(i+1,j+1,k)   + tmp(i-1,j+1,k)   + tmp(i+1,j-1,k)   + tmp(i-1,j-1,k)    &
                    + tmp(i+1,j+1,k+1) + tmp(i-1,j+1,k+1) + tmp(i+1,j-1,k+1) + tmp(i-1,j-1,k+1))) / tvol_z(i,j,k)
     END IF

  END DO
  END DO
  END DO

END SUBROUTINE test_filter_z

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE test_filter_xy(in, out)
  REAL(8), INTENT(IN)  :: in( -slv:isize+slv, -slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8), INTENT(OUT) :: out(-slv:isize+slv, -slv:jsize+slv, 1-slv:ksize+slv)

  REAL(8) :: tmp(0:isize+1, 0:jsize+1, 0:ksize+1)

  INTEGER :: i, j, k

  DO k=0, ksize+1
  DO j=0, jsize+1
  DO i=0, isize+1
     tmp(i,j,k) = imask3d(i,j,k) * 0.25*(in(i-1,j-1,k)+in(i,j-1,k)+in(i-1,j,k)+in(i,j,k)) * dvol(i,j,k)
  END DO
  END DO
  END DO

  DO k=1, ksize
  DO j=0, jsize
  DO i=0, isize
     IF (tvol_xy(i,j,k) == 0.0) THEN
        out(i,j,k) = 0.0
     ELSE
        out(i,j,k) = (tmp(i,j,k)   + tmp(i+1,j,k)   + tmp(i,j+1,k)   + tmp(i+1,j+1,k)   &
             + 0.5 * (tmp(i,j,k+1) + tmp(i+1,j,k+1) + tmp(i,j+1,k+1) + tmp(i+1,j+1,k+1) &
                    + tmp(i,j,k-1) + tmp(i+1,j,k-1) + tmp(i,j+1,k-1) + tmp(i+1,j+1,k-1))) / tvol_xy(i,j,k)
     END IF
  END DO
  END DO
  END DO

END SUBROUTINE test_filter_xy

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE test_filter_yz(in, out)
  REAL(8), INTENT(IN)  :: in( 1-slv:isize+slv, -slv:jsize+slv, -slv:ksize+slv)
  REAL(8), INTENT(OUT) :: out(1-slv:isize+slv, -slv:jsize+slv, -slv:ksize+slv)

  REAL(8) :: tmp(0:isize+1, 0:jsize+1, 0:ksize+1)

  INTEGER :: i, j, k

  DO k=0, ksize+1
  DO j=0, jsize+1
  DO i=0, isize+1
     tmp(i,j,k) = imask3d(i,j,k) * 0.25*(in(i,j-1,k-1)+in(i,j,k-1)+in(i,j-1,k)+in(i,j,k)) * dvol(i,j,k)
  END DO
  END DO
  END DO

  DO k=0, ksize
  DO j=0, jsize
  DO i=1, isize
     IF (tvol_yz(i,j,k) == 0.0) THEN
        out(i,j,k) = 0.0
     ELSE
        out(i,j,k) = (tmp(i,  j,k) + tmp(i,  j+1,k) + tmp(i,  j,k+1) + tmp(i,  j+1,k+1) &
             + 0.5 * (tmp(i+1,j,k) + tmp(i+1,j+1,k) + tmp(i+1,j,k+1) + tmp(i+1,j+1,k+1) &
                    + tmp(i-1,j,k) + tmp(i-1,j+1,k) + tmp(i-1,j,k+1) + tmp(i-1,j+1,k+1))) / tvol_yz(i,j,k)
     END IF
  END DO
  END DO
  END DO

END SUBROUTINE test_filter_yz

!-----------------------------------------------------------------------------------------------------------------------


SUBROUTINE test_filter_zx(in, out)
  REAL(8), INTENT(IN)  :: in( -slv:isize+slv, 1-slv:jsize+slv, -slv:ksize+slv)
  REAL(8), INTENT(OUT) :: out(-slv:isize+slv, 1-slv:jsize+slv, -slv:ksize+slv)

  REAL(8) :: tmp(0:isize+1, 0:jsize+1, 0:ksize+1)

  INTEGER :: i, j, k

  DO k=0, ksize+1
  DO j=0, jsize+1
  DO i=0, isize+1
     tmp(i,j,k) = imask3d(i,j,k) * 0.25*(in(i-1,j,k-1)+in(i-1,j,k)+in(i,j,k-1)+in(i,j,k)) * dvol(i,j,k)
  END DO
  END DO
  END DO

  DO k=0, ksize
  DO j=1, jsize
  DO i=0, isize
     IF (tvol_zx(i,j,k) == 0.0) THEN
        out(i,j,k) = 0.0
     ELSE
        out(i,j,k) = (tmp(i,j,  k) + tmp(i,j,  k+1) + tmp(i+1,j,  k) + tmp(i+1,j,  k+1) &
             + 0.5 * (tmp(i,j+1,k) + tmp(i,j+1,k+1) + tmp(i+1,j+1,k) + tmp(i+1,j+1,k+1) &
                    + tmp(i,j-1,k) + tmp(i,j-1,k+1) + tmp(i+1,j-1,k) + tmp(i+1,j-1,k+1))) / tvol_zx(i,j,k)
     END IF
  END DO
  END DO
  END DO

END SUBROUTINE test_filter_zx

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE init_dles
  INTEGER :: i, j, k
  REAL(8) :: dvol_masked(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)

  ALLOCATE(u_test( -slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv))
  ALLOCATE(v_test(1-slv:isize+slv,  -slv:jsize+slv, 1-slv:ksize+slv))
  ALLOCATE(w_test(1-slv:isize+slv, 1-slv:jsize+slv,  -slv:ksize+slv))
  u_test(:,:,:) = 0.0
  v_test(:,:,:) = 0.0
  w_test(:,:,:) = 0.0

  ALLOCATE(srate_test(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv))
  srate_test(:,:,:) = 0.0

  ALLOCATE(dles_flag(1:isize, 1:jsize, 1:ksize))
  dles_flag(:,:,:) = 0.0

  ALLOCATE(nfreq2_test(1-slv:isize+slv, 1-slv:jsize+slv, -slv:ksize+slv))
  nfreq2_test(:,:,:) = 0.0

  ALLOCATE(tvol(   1:isize, 1:jsize, 1:ksize))
  ALLOCATE(tvol_x( 0:isize, 1:jsize, 1:ksize))
  ALLOCATE(tvol_y( 1:isize, 0:jsize, 1:ksize))
  ALLOCATE(tvol_z( 1:isize, 1:jsize, 0:ksize))
  ALLOCATE(tvol_xy(0:isize, 0:jsize, 1:ksize))
  ALLOCATE(tvol_yz(1:isize, 0:jsize, 0:ksize))
  ALLOCATE(tvol_zx(0:isize, 1:jsize, 0:ksize))

  dvol_masked(:,:,:)=dvol(:,:,:)*imask3d(:,:,:)

  DO k=1, ksize
  DO j=1, jsize
  DO i=1, isize
     tvol(i,j,k) = imask3d(i,j,k)*(dvol_masked(i,j,k) &
          + 0.5   * (dvol_masked(i+1,j,  k)   + dvol_masked(i-1,j,  k)    &
                   + dvol_masked(i,  j+1,k)   + dvol_masked(i,  j-1,k)    &
                   + dvol_masked(i,  j,  k+1) + dvol_masked(i,  j,  k-1)) &
          + 0.25  * (dvol_masked(i+1,j+1,k)   + dvol_masked(i-1,j+1,k)   + dvol_masked(i+1,j-1,k)   + dvol_masked(i-1,j-1,k)     &
                   + dvol_masked(i  ,j+1,k+1) + dvol_masked(i  ,j-1,k+1) + dvol_masked(i  ,j+1,k-1) + dvol_masked(i  ,j-1,k-1)   &
                   + dvol_masked(i+1,j,  k+1) + dvol_masked(i+1,j,  k-1) + dvol_masked(i-1,j,  k+1) + dvol_masked(i-1,j,  k-1))  &
          + 0.125 * (dvol_masked(i+1,j+1,k+1) + dvol_masked(i-1,j+1,k+1) + dvol_masked(i+1,j-1,k+1) + dvol_masked(i-1,j-1,k+1)   &
                   + dvol_masked(i+1,j+1,k-1) + dvol_masked(i-1,j+1,k-1) + dvol_masked(i+1,j-1,k-1) + dvol_masked(i-1,j-1,k-1)))
  END DO
  END DO
  END DO

  DO k=1, ksize
  DO j=1, jsize
  DO i=0, isize
     tvol_x(i,j,k) = (dvol_masked(i,j,k) + dvol_masked(i+1,j,k) &
          + 0.5  * (dvol_masked(i  ,j+1,k)   + dvol_masked(i,  j-1,k)   + dvol_masked(i  ,j,  k+1) + dvol_masked(i,  j,  k-1)  &
                  + dvol_masked(i+1,j+1,k)   + dvol_masked(i+1,j-1,k)   + dvol_masked(i+1,j,  k+1) + dvol_masked(i+1,j,  k-1)) &
          + 0.25 * (dvol_masked(i,  j+1,k+1) + dvol_masked(i,  j-1,k+1) + dvol_masked(i,  j+1,k-1) + dvol_masked(i,  j-1,k-1)  &
                  + dvol_masked(i+1,j+1,k+1) + dvol_masked(i+1,j-1,k+1) + dvol_masked(i+1,j+1,k-1) + dvol_masked(i+1,j-1,k-1)))
  END DO
  END DO
  END DO

  DO k=1, ksize
  DO j=1, jsize
  DO i=0, isize
     tvol_x(i,j,k) = (dvol_masked(i,j,k) + dvol_masked(i+1,j,k) &
          + 0.5  * (dvol_masked(i  ,j+1,k)   + dvol_masked(i,  j-1,k)   + dvol_masked(i  ,j,  k+1) + dvol_masked(i,  j,  k-1)  &
                  + dvol_masked(i+1,j+1,k)   + dvol_masked(i+1,j-1,k)   + dvol_masked(i+1,j,  k+1) + dvol_masked(i+1,j,  k-1)) &
          + 0.25 * (dvol_masked(i,  j+1,k+1) + dvol_masked(i,  j-1,k+1) + dvol_masked(i,  j+1,k-1) + dvol_masked(i,  j-1,k-1)  &
                  + dvol_masked(i+1,j+1,k+1) + dvol_masked(i+1,j-1,k+1) + dvol_masked(i+1,j+1,k-1) + dvol_masked(i+1,j-1,k-1)))
  END DO
  END DO
  END DO

  DO k=1, ksize
  DO j=0, jsize
  DO i=1, isize
     tvol_y(i,j,k) = (dvol_masked(i,j,k) + dvol_masked(i,j+1,k) &
          + 0.5  * (dvol_masked(i,  j,  k+1) + dvol_masked(i,  j,  k-1) + dvol_masked(i+1,j,  k)   + dvol_masked(i-1,j,  k)    &
                  + dvol_masked(i,  j+1,k+1) + dvol_masked(i,  j+1,k-1) + dvol_masked(i+1,j+1,k)   + dvol_masked(i-1,j+1,k))   &
          + 0.25 * (dvol_masked(i+1,j,  k+1) + dvol_masked(i+1,j,  k-1) + dvol_masked(i-1,  j,k+1) + dvol_masked(i-1,j,  k-1)  &
                  + dvol_masked(i+1,j+1,k+1) + dvol_masked(i+1,j+1,k-1) + dvol_masked(i-1,j+1,k+1) + dvol_masked(i-1,j+1,k-1)))
  END DO
  END DO
  END DO

  DO k=0, ksize
  DO j=1, jsize
  DO i=1, isize
     tvol_z(i,j,k) = (dvol_masked(i,j,k) + dvol_masked(i,j,k+1) &
          + 0.5  * (dvol_masked(i+1,j,  k)   + dvol_masked(i-1,j,  k)   + dvol_masked(i,  j+1,k)   + dvol_masked(i,  j-1,k)    &
                  + dvol_masked(i+1,j,  k+1) + dvol_masked(i-1,j,  k+1) + dvol_masked(i,  j+1,k+1) + dvol_masked(i,  j-1,k+1)) &
          + 0.25 * (dvol_masked(i+1,j+1,k)   + dvol_masked(i-1,j+1,k)   + dvol_masked(i+1,j-1,k)   + dvol_masked(i-1,j-1,k)    &
                  + dvol_masked(i+1,j+1,k+1) + dvol_masked(i-1,j+1,k+1) + dvol_masked(i+1,j-1,k+1) + dvol_masked(i-1,j-1,k+1)))
  END DO
  END DO
  END DO

  DO k=1, ksize
  DO j=0, jsize
  DO i=0, isize
     tvol_xy(i,j,k) = (dvol_masked(i,j,k)   + dvol_masked(i+1,j,k)   + dvol_masked(i,j+1,k)   + dvol_masked(i+1,j+1,k)    &
          + 0.5 * (    dvol_masked(i,j,k+1) + dvol_masked(i+1,j,k+1) + dvol_masked(i,j+1,k+1) + dvol_masked(i+1,j+1,k+1)  &
                     + dvol_masked(i,j,k-1) + dvol_masked(i+1,j,k-1) + dvol_masked(i,j+1,k-1) + dvol_masked(i+1,j+1,k-1)))
  END DO
  END DO
  END DO

  DO k=0, ksize
  DO j=0, jsize
  DO i=1, isize
     tvol_yz(i,j,k) = (dvol_masked(i,  j,k) + dvol_masked(i,  j+1,k) + dvol_masked(i,  j,k+1) + dvol_masked(i,  j+1,k+1)  &
          + 0.5 * (    dvol_masked(i+1,j,k) + dvol_masked(i+1,j+1,k) + dvol_masked(i+1,j,k+1) + dvol_masked(i+1,j+1,k+1)  &
                     + dvol_masked(i-1,j,k) + dvol_masked(i-1,j+1,k) + dvol_masked(i-1,j,k+1) + dvol_masked(i-1,j+1,k+1)))
  END DO
  END DO
  END DO

  DO k=0, ksize
  DO j=1, jsize
  DO i=0, isize
     tvol_zx(i,j,k) = (dvol_masked(i,j,  k) + dvol_masked(i,j,  k+1) + dvol_masked(i+1,j,  k) + dvol_masked(i+1,j,  k+1)  &
          + 0.5 * (    dvol_masked(i,j+1,k) + dvol_masked(i,j+1,k+1) + dvol_masked(i+1,j+1,k) + dvol_masked(i+1,j+1,k+1)  &
                     + dvol_masked(i,j-1,k) + dvol_masked(i,j-1,k+1) + dvol_masked(i+1,j-1,k) + dvol_masked(i+1,j-1,k+1)))
  END DO
  END DO
  END DO

  ALLOCATE(dles_tmp1(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv))
  ALLOCATE(dles_tmp2(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv))

  dles_tmp1(:,:,:) = 0.0
  dles_tmp2(:,:,:) = 0.0

END SUBROUTINE init_dles

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE reset_test_volume
  REAL(8) :: dvol_masked(1-slv:isize+slv, 1-slv:jsize+slv, ksize-2:ksize+1)

  INTEGER :: i, j, k

  dvol_masked(:,:,ksize-2:ksize+1)=dvol(:,:,ksize-2:ksize+1)*imask3d(:,:,ksize-2:ksize+1)

  DO k=ksize-1, ksize
  DO j=1, jsize
  DO i=1, isize
     tvol(i,j,k) = imask3d(i,j,k)*(dvol_masked(i,j,k) &
          + 0.5   * (dvol_masked(i+1,j,  k)   + dvol_masked(i-1,j,  k)    &
                   + dvol_masked(i,  j+1,k)   + dvol_masked(i,  j-1,k)    &
                   + dvol_masked(i,  j,  k+1) + dvol_masked(i,  j,  k-1)) &
          + 0.25  * (dvol_masked(i+1,j+1,k)   + dvol_masked(i-1,j+1,k)   + dvol_masked(i+1,j-1,k)   + dvol_masked(i-1,j-1,k)     &
                   + dvol_masked(i  ,j+1,k+1) + dvol_masked(i  ,j-1,k+1) + dvol_masked(i  ,j+1,k-1) + dvol_masked(i  ,j-1,k-1)   &
                   + dvol_masked(i+1,j,  k+1) + dvol_masked(i+1,j,  k-1) + dvol_masked(i-1,j,  k+1) + dvol_masked(i-1,j,  k-1))  &
          + 0.125 * (dvol_masked(i+1,j+1,k+1) + dvol_masked(i-1,j+1,k+1) + dvol_masked(i+1,j-1,k+1) + dvol_masked(i-1,j-1,k+1)   &
                   + dvol_masked(i+1,j+1,k-1) + dvol_masked(i-1,j+1,k-1) + dvol_masked(i+1,j-1,k-1) + dvol_masked(i-1,j-1,k-1)))
  END DO
  END DO
  END DO

  DO k=ksize-1, ksize
  DO j=1, jsize
  DO i=0, isize
     tvol_x(i,j,k) = (dvol_masked(i,j,k) + dvol_masked(i+1,j,k) &
          + 0.5  * (dvol_masked(i  ,j+1,k)   + dvol_masked(i,  j-1,k)   + dvol_masked(i  ,j,  k+1) + dvol_masked(i,  j,  k-1)  &
                  + dvol_masked(i+1,j+1,k)   + dvol_masked(i+1,j-1,k)   + dvol_masked(i+1,j,  k+1) + dvol_masked(i+1,j,  k-1)) &
          + 0.25 * (dvol_masked(i,  j+1,k+1) + dvol_masked(i,  j-1,k+1) + dvol_masked(i,  j+1,k-1) + dvol_masked(i,  j-1,k-1)  &
                  + dvol_masked(i+1,j+1,k+1) + dvol_masked(i+1,j-1,k+1) + dvol_masked(i+1,j+1,k-1) + dvol_masked(i+1,j-1,k-1)))
  END DO
  END DO
  END DO

  DO k=ksize-1, ksize
  DO j=1, jsize
  DO i=0, isize
     tvol_x(i,j,k) = (dvol_masked(i,j,k) + dvol_masked(i+1,j,k) &
          + 0.5  * (dvol_masked(i  ,j+1,k)   + dvol_masked(i,  j-1,k)   + dvol_masked(i  ,j,  k+1) + dvol_masked(i,  j,  k-1)  &
                  + dvol_masked(i+1,j+1,k)   + dvol_masked(i+1,j-1,k)   + dvol_masked(i+1,j,  k+1) + dvol_masked(i+1,j,  k-1)) &
          + 0.25 * (dvol_masked(i,  j+1,k+1) + dvol_masked(i,  j-1,k+1) + dvol_masked(i,  j+1,k-1) + dvol_masked(i,  j-1,k-1)  &
                  + dvol_masked(i+1,j+1,k+1) + dvol_masked(i+1,j-1,k+1) + dvol_masked(i+1,j+1,k-1) + dvol_masked(i+1,j-1,k-1)))
  END DO
  END DO
  END DO

  DO k=ksize-1, ksize
  DO j=0, jsize
  DO i=1, isize
     tvol_y(i,j,k) = (dvol_masked(i,j,k) + dvol_masked(i,j+1,k) &
          + 0.5  * (dvol_masked(i,  j,  k+1) + dvol_masked(i,  j,  k-1) + dvol_masked(i+1,j,  k)   + dvol_masked(i-1,j,  k)    &
                  + dvol_masked(i,  j+1,k+1) + dvol_masked(i,  j+1,k-1) + dvol_masked(i+1,j+1,k)   + dvol_masked(i-1,j+1,k))   &
          + 0.25 * (dvol_masked(i+1,j,  k+1) + dvol_masked(i+1,j,  k-1) + dvol_masked(i-1,  j,k+1) + dvol_masked(i-1,j,  k-1)  &
                  + dvol_masked(i+1,j+1,k+1) + dvol_masked(i+1,j+1,k-1) + dvol_masked(i-1,j+1,k+1) + dvol_masked(i-1,j+1,k-1)))
  END DO
  END DO
  END DO

  DO k=ksize-1, ksize
  DO j=1, jsize
  DO i=1, isize
     tvol_z(i,j,k) = (dvol_masked(i,j,k) + dvol_masked(i,j,k+1) &
          + 0.5  * (dvol_masked(i+1,j,  k)   + dvol_masked(i-1,j,  k)   + dvol_masked(i,  j+1,k)   + dvol_masked(i,  j-1,k)    &
                  + dvol_masked(i+1,j,  k+1) + dvol_masked(i-1,j,  k+1) + dvol_masked(i,  j+1,k+1) + dvol_masked(i,  j-1,k+1)) &
          + 0.25 * (dvol_masked(i+1,j+1,k)   + dvol_masked(i-1,j+1,k)   + dvol_masked(i+1,j-1,k)   + dvol_masked(i-1,j-1,k)    &
                  + dvol_masked(i+1,j+1,k+1) + dvol_masked(i-1,j+1,k+1) + dvol_masked(i+1,j-1,k+1) + dvol_masked(i-1,j-1,k+1)))
  END DO
  END DO
  END DO

  DO k=ksize-1, ksize
  DO j=0, jsize
  DO i=0, isize
     tvol_xy(i,j,k) = (dvol_masked(i,j,k)   + dvol_masked(i+1,j,k)   + dvol_masked(i,j+1,k)   + dvol_masked(i+1,j+1,k)    &
          + 0.5 * (    dvol_masked(i,j,k+1) + dvol_masked(i+1,j,k+1) + dvol_masked(i,j+1,k+1) + dvol_masked(i+1,j+1,k+1)  &
                     + dvol_masked(i,j,k-1) + dvol_masked(i+1,j,k-1) + dvol_masked(i,j+1,k-1) + dvol_masked(i+1,j+1,k-1)))
  END DO
  END DO
  END DO

  DO k=ksize-1, ksize
  DO j=0, jsize
  DO i=1, isize
     tvol_yz(i,j,k) = (dvol_masked(i,  j,k) + dvol_masked(i,  j+1,k) + dvol_masked(i,  j,k+1) + dvol_masked(i,  j+1,k+1)  &
          + 0.5 * (    dvol_masked(i+1,j,k) + dvol_masked(i+1,j+1,k) + dvol_masked(i+1,j,k+1) + dvol_masked(i+1,j+1,k+1)  &
                    +  dvol_masked(i-1,j,k) + dvol_masked(i-1,j+1,k) + dvol_masked(i-1,j,k+1) + dvol_masked(i-1,j+1,k+1)))
  END DO
  END DO
  END DO

  DO k=ksize-1, ksize
  DO j=1, jsize
  DO i=0, isize
     tvol_zx(i,j,k) = (dvol_masked(i,j,  k) + dvol_masked(i,j,  k+1) + dvol_masked(i+1,j,  k) + dvol_masked(i+1,j,  k+1)  &
          + 0.5 * (    dvol_masked(i,j+1,k) + dvol_masked(i,j+1,k+1) + dvol_masked(i+1,j+1,k) + dvol_masked(i+1,j+1,k+1)  &
                     + dvol_masked(i,j-1,k) + dvol_masked(i,j-1,k+1) + dvol_masked(i+1,j-1,k) + dvol_masked(i+1,j-1,k+1)))
  END DO
  END DO
  END DO

END SUBROUTINE reset_test_volume

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE finalize_dles
  IF (ALLOCATED(u_test)) DEALLOCATE(u_test)
  IF (ALLOCATED(v_test)) DEALLOCATE(v_test)
  IF (ALLOCATED(w_test)) DEALLOCATE(w_test)

  IF (ALLOCATED(nfreq2_test)) DEALLOCATE(nfreq2_test)
  IF (ALLOCATED(srate_test)) DEALLOCATE(srate_test)

  IF (ALLOCATED(dles_flag)) DEALLOCATE(dles_flag)

  IF (ALLOCATED(tvol))    DEALLOCATE(tvol)
  IF (ALLOCATED(tvol_x))  DEALLOCATE(tvol_x)
  IF (ALLOCATED(tvol_y))  DEALLOCATE(tvol_y)
  IF (ALLOCATED(tvol_z))  DEALLOCATE(tvol_z)
  IF (ALLOCATED(tvol_xy)) DEALLOCATE(tvol_xy)
  IF (ALLOCATED(tvol_yz)) DEALLOCATE(tvol_yz)
  IF (ALLOCATED(tvol_zx)) DEALLOCATE(tvol_zx)

END SUBROUTINE finalize_dles

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE dles_viscosity(u_grid, v_grid, w_grid, u_test, v_test, w_test, nfreq2_grid, nfreq2_test, &
                          n_x, n_y, n_z, nu1, nu2, flag)
  USE parameters, ONLY : subgrid_isotropic

  REAL(8), INTENT(IN) :: u_grid( -slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8), INTENT(IN) :: v_grid(1-slv:isize+slv,  -slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8), INTENT(IN) :: w_grid(1-slv:isize+slv, 1-slv:jsize+slv,  -slv:ksize+slv)

  REAL(8), INTENT(IN) :: u_test( -slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8), INTENT(IN) :: v_test(1-slv:isize+slv,  -slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8), INTENT(IN) :: w_test(1-slv:isize+slv, 1-slv:jsize+slv,  -slv:ksize+slv)

  REAL(8), INTENT(IN) :: nfreq2_grid(1-slv:isize+slv, 1-slv:jsize+slv, -slv:ksize+slv)
  REAL(8), INTENT(IN) :: nfreq2_test(1-slv:isize+slv, 1-slv:jsize+slv, -slv:ksize+slv)

  REAL(8), INTENT(IN) :: n_x( -slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8), INTENT(IN) :: n_y(1-slv:isize+slv,  -slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8), INTENT(IN) :: n_z(1-slv:isize+slv, 1-slv:jsize+slv,  -slv:ksize+slv)

  REAL(8), INTENT(OUT) :: nu1(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8), INTENT(OUT) :: nu2(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)

  REAL(8), OPTIONAL, INTENT(OUT) :: flag(1:isize, 1:jsize, 1:ksize)

  REAL(8) :: srate_grid(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: srate_test(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)

  REAL(8) :: sr_xx(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: sr_yy(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: sr_zz(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: sr_xy( -slv:isize+slv,  -slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: sr_yz(1-slv:isize+slv,  -slv:jsize+slv,  -slv:ksize+slv)
  REAL(8) :: sr_zx( -slv:isize+slv, 1-slv:jsize+slv,  -slv:ksize+slv)

  REAL(8) :: sn_xx(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: sn_yy(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: sn_zz(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: sn_xy( -slv:isize+slv,  -slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: sn_yz(1-slv:isize+slv,  -slv:jsize+slv,  -slv:ksize+slv)
  REAL(8) :: sn_zx( -slv:isize+slv, 1-slv:jsize+slv,  -slv:ksize+slv)

  REAL(8) :: lxx(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: lyy(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: lzz(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: lxy( -slv:isize+slv,  -slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: lyz(1-slv:isize+slv,  -slv:jsize+slv,  -slv:ksize+slv)
  REAL(8) :: lzx( -slv:isize+slv, 1-slv:jsize+slv,  -slv:ksize+slv)

  REAL(8) :: mxx(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: myy(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: mzz(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: mxy( -slv:isize+slv,  -slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: myz(1-slv:isize+slv,  -slv:jsize+slv,  -slv:ksize+slv)
  REAL(8) :: mzx( -slv:isize+slv, 1-slv:jsize+slv,  -slv:ksize+slv)

  REAL(8) :: nxx(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: nyy(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: nzz(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: nxy( -slv:isize+slv,  -slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: nyz(1-slv:isize+slv,  -slv:jsize+slv,  -slv:ksize+slv)
  REAL(8) :: nzx( -slv:isize+slv, 1-slv:jsize+slv,  -slv:ksize+slv)

  REAL(8) :: lm(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: ln(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: ll(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: mm(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: nn(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: mn(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)

  REAL(8) :: det(1:isize, 1:jsize, 1:ksize)
  REAL(8) :: tr(1:isize, 1:jsize, 1:ksize)

  REAL(8) :: tmp_xx(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: tmp_yy(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: tmp_zz(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: tmp_xy( -slv:isize+slv,  -slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: tmp_yz(1-slv:isize+slv,  -slv:jsize+slv,  -slv:ksize+slv)
  REAL(8) :: tmp_zx( -slv:isize+slv, 1-slv:jsize+slv,  -slv:ksize+slv)

  INTEGER :: i, j, k

  tmp_xx(:,:,:) = 0.0
  tmp_yy(:,:,:) = 0.0
  tmp_zz(:,:,:) = 0.0
  tmp_xy(:,:,:) = 0.0
  tmp_yz(:,:,:) = 0.0
  tmp_zx(:,:,:) = 0.0

  DO k=0, ksize+1
  DO j=0, jsize+1
  DO i=0, isize+1
     tmp_xx(i,j,k) = ((u_grid(i-1,j,k)+u_grid(i,j,k))*0.5)**2
     tmp_yy(i,j,k) = ((v_grid(i,j-1,k)+v_grid(i,j,k))*0.5)**2
     tmp_zz(i,j,k) = ((w_grid(i,j,k-1)+w_grid(i,j,k))*0.5)**2
  END DO
  END DO
  END DO

  DO k= 0, ksize+1
  DO j=-1, jsize+1
  DO i=-1, isize+1
     tmp_xy(i,j,k) = (u_grid(i,j,k)+u_grid(i,j+1,k))*0.5 * (v_grid(i,j,k)+v_grid(i+1,j,k))*0.5
  END DO
  END DO
  END DO

  DO k=-1, ksize+1
  DO j=-1, jsize+1
  DO i= 0, isize+1
     tmp_yz(i,j,k) = (v_grid(i,j,k)+v_grid(i,j,k+1))*0.5 * (w_grid(i,j,k)+w_grid(i,j+1,k))*0.5
  END DO
  END DO
  END DO

  DO k=-1, ksize+1
  DO j= 0, jsize+1
  DO i=-1, isize+1
     tmp_zx(i,j,k) = (w_grid(i,j,k)+w_grid(i+1,j,k))*0.5 * (u_grid(i,j,k)+u_grid(i,j,k+1))*0.5
  END DO
  END DO
  END DO

  CALL test_filter(tmp_xx, lxx)
  CALL test_filter(tmp_yy, lyy)
  CALL test_filter(tmp_zz, lzz)
  CALL test_filter_xy(tmp_xy, lxy)
  CALL test_filter_yz(tmp_yz, lyz)
  CALL test_filter_zx(tmp_zx, lzx)

  DO k=1, ksize
  DO j=1, jsize
  DO i=1, isize
     lxx(i,j,k) = lxx(i,j,k) - ((u_test(i-1,j,k)+u_test(i,j,k))*0.5)**2
     lyy(i,j,k) = lyy(i,j,k) - ((v_test(i,j-1,k)+v_test(i,j,k))*0.5)**2
     lzz(i,j,k) = lzz(i,j,k) - ((w_test(i,j,k-1)+w_test(i,j,k))*0.5)**2

     tr(i,j,k) = lxx(i,j,k)+lyy(i,j,k)+lzz(i,j,k)
     lxx(i,j,k) = lxx(i,j,k) - tr(i,j,k)/3.0
     lyy(i,j,k) = lyy(i,j,k) - tr(i,j,k)/3.0
     lzz(i,j,k) = lzz(i,j,k) - tr(i,j,k)/3.0
  END DO
  END DO
  END DO

  DO k=1, ksize
  DO j=0, jsize
  DO i=0, isize
     lxy(i,j,k) = lxy(i,j,k) - (u_test(i,j,k)+u_test(i,j+1,k))*0.5 * (v_test(i,j,k)+v_test(i+1,j,k))*0.5
  END DO
  END DO
  END DO

  DO k=0, ksize
  DO j=0, jsize
  DO i=1, isize
     lyz(i,j,k) = lyz(i,j,k) - (v_test(i,j,k)+v_test(i,j,k+1))*0.5 * (w_test(i,j,k)+w_test(i,j+1,k))*0.5
  END DO
  END DO
  END DO

  DO k=0, ksize
  DO j=1, jsize
  DO i=0, isize
     lzx(i,j,k) = lzx(i,j,k) - (w_test(i,j,k)+w_test(i+1,j,k))*0.5 * (u_test(i,j,k)+u_test(i,j,k+1))*0.5
  END DO
  END DO
  END DO

  CALL sr_tensor(u_grid, v_grid, w_grid, sr_xx, sr_yy, sr_zz, sr_xy, sr_yz, sr_zx)

  CALL sn_tensor(sr_xx, sr_yy, sr_zz, sr_xy, sr_yz, sr_zx, n_x, n_y, n_z, sn_xx, sn_yy, sn_zz, sn_xy, sn_yz, sn_zx)

  srate_grid(:,:,:) = 0.0

  DO k=1, ksize
  DO j=1, jsize
  DO i=1, isize
     srate_grid(i,j,k) = sqrt(2.0*(sr_xx(i,j,k)**2 + sr_yy(i,j,k)**2 + sr_zz(i,j,k)**2) &
          + 4.0*( ((sr_xy(i-1,j-1,k)+sr_xy(i,j-1,k)+sr_xy(i-1,j,k)+sr_xy(i,j,k))*0.25)**2 &
                  + ((sr_yz(i,j-1,k-1)+sr_yz(i,j,k-1)+sr_yz(i,j-1,k)+sr_yz(i,j,k))*0.25)**2 &
                  + ((sr_zx(i-1,j,k-1)+sr_zx(i-1,j,k)+sr_zx(i,j,k-1)+sr_zx(i,j,k))*0.25)**2))
  END DO
  END DO
  END DO

  CALL update_boundary(srate_grid)

  DO k=0, ksize+1
  DO j=0, jsize+1
  DO i=0, isize+1
     tmp_xx(i,j,k) = sr_xx(i,j,k) * srate_grid(i,j,k)
     tmp_yy(i,j,k) = sr_yy(i,j,k) * srate_grid(i,j,k)
     tmp_zz(i,j,k) = sr_zz(i,j,k) * srate_grid(i,j,k)
  END DO
  END DO
  END DO

  DO k= 0, ksize+1
  DO j=-1, jsize+1
  DO i=-1, isize+1
     tmp_xy(i,j,k) = sr_xy(i,j,k) &
          * (srate_grid(i,j,k)+srate_grid(i+1,j,k)+srate_grid(i,j+1,k)+srate_grid(i+1,j+1,k))*0.25
  END DO
  END DO
  END DO

  DO k=-1, ksize+1
  DO j=-1, jsize+1
  DO i= 0, isize+1
     tmp_yz(i,j,k) = sr_yz(i,j,k) &
          * (srate_grid(i,j,k)+srate_grid(i,j+1,k)+srate_grid(i,j,k+1)+srate_grid(i,j+1,k+1))*0.25
  END DO
  END DO
  END DO

  DO k=-1, ksize+1
  DO j= 0, jsize+1
  DO i=-1, isize+1
     tmp_zx(i,j,k) = sr_zx(i,j,k) &
          * (srate_grid(i,j,k)+srate_grid(i,j,k+1)+srate_grid(i+1,j,k)+srate_grid(i+1,j,k+1))*0.25
  END DO
  END DO
  END DO

  CALL test_filter(tmp_xx, mxx)
  CALL test_filter(tmp_yy, myy)
  CALL test_filter(tmp_zz, mzz)
  CALL test_filter_xy(tmp_xy, mxy)
  CALL test_filter_yz(tmp_yz, myz)
  CALL test_filter_zx(tmp_zx, mzx)

  DO k=0, ksize+1
  DO j=0, jsize+1
  DO i=0, isize+1
     tmp_xx(i,j,k) = sn_xx(i,j,k) * nfreq2_grid(i,j,k)
     tmp_yy(i,j,k) = sn_yy(i,j,k) * nfreq2_grid(i,j,k)
     tmp_zz(i,j,k) = sn_zz(i,j,k) * nfreq2_grid(i,j,k)
  END DO
  END DO
  END DO

  DO k= 0, ksize+1
  DO j=-1, jsize+1
  DO i=-1, isize+1
     tmp_xy(i,j,k) = sn_xy(i,j,k) &
          * (nfreq2_grid(i,j,k)+nfreq2_grid(i+1,j,k)+nfreq2_grid(i,j+1,k)+nfreq2_grid(i+1,j+1,k))*0.25
  END DO
  END DO
  END DO

  DO k=-1, ksize+1
  DO j=-1, jsize+1
  DO i= 0, isize+1
     tmp_yz(i,j,k) = sn_yz(i,j,k) &
          * (nfreq2_grid(i,j,k)+nfreq2_grid(i,j+1,k)+nfreq2_grid(i,j,k+1)+nfreq2_grid(i,j+1,k+1))*0.25
  END DO
  END DO
  END DO

  DO k=-1, ksize+1
  DO j= 0, jsize+1
  DO i=-1, isize+1
     tmp_zx(i,j,k) = sn_zx(i,j,k) &
          * (nfreq2_grid(i,j,k)+nfreq2_grid(i,j,k+1)+nfreq2_grid(i+1,j,k)+nfreq2_grid(i+1,j,k+1))*0.25
  END DO
  END DO
  END DO

  CALL test_filter(tmp_xx, nxx)
  CALL test_filter(tmp_yy, nyy)
  CALL test_filter(tmp_zz, nzz)
  CALL test_filter_xy(tmp_xy, nxy)
  CALL test_filter_yz(tmp_yz, nyz)
  CALL test_filter_zx(tmp_zx, nzx)

  CALL sr_tensor(u_test, v_test, w_test, sr_xx, sr_yy, sr_zz, sr_xy, sr_yz, sr_zx)
  CALL sn_tensor(sr_xx, sr_yy, sr_zz, sr_xy, sr_yz, sr_zx, n_x, n_y, n_z, sn_xx, sn_yy, sn_zz, sn_xy, sn_yz, sn_zx)

  DO k=0, ksize+1
  DO j=0, jsize+1
  DO i=0, isize+1
     srate_test(i,j,k) = sqrt(2.0*(sr_xx(i,j,k)**2 + sr_yy(i,j,k)**2 + sr_zz(i,j,k)**2) &
          + 4.0*( ((sr_xy(i-1,j-1,k)+sr_xy(i,j-1,k)+sr_xy(i-1,j,k)+sr_xy(i,j,k))*0.25)**2 &
                  + ((sr_yz(i,j-1,k-1)+sr_yz(i,j,k-1)+sr_yz(i,j-1,k)+sr_yz(i,j,k))*0.25)**2 &
                  + ((sr_zx(i-1,j,k-1)+sr_zx(i-1,j,k)+sr_zx(i,j,k-1)+sr_zx(i,j,k))*0.25)**2))
  END DO
  END DO
  END DO

  DO k=1, ksize
  DO j=1, jsize
  DO i=1, isize
     mxx(i,j,k) = mxx(i,j,k) - tg_ratio2 * sr_xx(i,j,k) * srate_test(i,j,k)
     myy(i,j,k) = myy(i,j,k) - tg_ratio2 * sr_yy(i,j,k) * srate_test(i,j,k)
     mzz(i,j,k) = mzz(i,j,k) - tg_ratio2 * sr_zz(i,j,k) * srate_test(i,j,k)
  END DO
  END DO
  END DO

  DO k=1, ksize
  DO j=0, jsize
  DO i=0, isize
     mxy(i,j,k) = mxy(i,j,k) - tg_ratio2 * sr_xy(i,j,k) &
          * (srate_test(i,j,k)+srate_test(i+1,j,k)+srate_test(i,j+1,k)+srate_test(i+1,j+1,k))*0.25
  END DO
  END DO
  END DO

  DO k=0, ksize
  DO j=0, jsize
  DO i=1, isize
     myz(i,j,k) = myz(i,j,k) - tg_ratio2 * sr_yz(i,j,k) &
          * (srate_test(i,j,k)+srate_test(i,j+1,k)+srate_test(i,j,k+1)+srate_test(i,j+1,k+1))*0.25
  END DO
  END DO
  END DO

  DO k=0, ksize
  DO j=1, jsize
  DO i=0, isize
     mzx(i,j,k) = mzx(i,j,k) - tg_ratio2 * sr_zx(i,j,k) &
          * (srate_test(i,j,k)+srate_test(i,j,k+1)+srate_test(i+1,j,k)+srate_test(i+1,j,k+1))*0.25
  END DO
  END DO
  END DO

  DO k=1, ksize
  DO j=1, jsize
  DO i=1, isize
     nxx(i,j,k) = nxx(i,j,k) - tg_ratio2 * sn_xx(i,j,k) * nfreq2_test(i,j,k)
     nyy(i,j,k) = nyy(i,j,k) - tg_ratio2 * sn_yy(i,j,k) * nfreq2_test(i,j,k)
     nzz(i,j,k) = nzz(i,j,k) - tg_ratio2 * sn_zz(i,j,k) * nfreq2_test(i,j,k)
  END DO
  END DO
  END DO

  DO k=1, ksize
  DO j=0, jsize
  DO i=0, isize
     nxy(i,j,k) = nxy(i,j,k) - tg_ratio2 * sn_xy(i,j,k) &
          * (nfreq2_test(i,j,k)+nfreq2_test(i+1,j,k)+nfreq2_test(i,j+1,k)+nfreq2_test(i+1,j+1,k))*0.25
  END DO
  END DO
  END DO

  DO k=0, ksize
  DO j=0, jsize
  DO i=1, isize
     nyz(i,j,k) = nyz(i,j,k) - tg_ratio2 * sn_yz(i,j,k) &
          * (nfreq2_test(i,j,k)+nfreq2_test(i,j+1,k)+nfreq2_test(i,j,k+1)+nfreq2_test(i,j+1,k+1))*0.25
  END DO
  END DO
  END DO

  DO k=0, ksize
  DO j=1, jsize
  DO i=0, isize
     nzx(i,j,k) = nzx(i,j,k) - tg_ratio2 * sn_zx(i,j,k) &
          * (nfreq2_test(i,j,k)+nfreq2_test(i,j,k+1)+nfreq2_test(i+1,j,k)+nfreq2_test(i+1,j,k+1))*0.25
  END DO
  END DO
  END DO

  lm = 0.0
  ln = 0.0
  ll = 0.0
  mm = 0.0
  nn = 0.0
  mn = 0.0

  CALL contraction_tensor(lxx, lyy, lzz, lxy, lyz, lzx, mxx, myy, mzz, mxy, myz, mzx, lm)
  CALL contraction_tensor(lxx, lyy, lzz, lxy, lyz, lzx, nxx, nyy, nzz, nxy, nyz, nzx, ln)
  CALL contraction_tensor(lxx, lyy, lzz, lxy, lyz, lzx, lxx, lyy, lzz, lxy, lyz, lzx, ll)
  CALL contraction_tensor(mxx, myy, mzz, mxy, myz, mzx, mxx, myy, mzz, mxy, myz, mzx, mm)
  CALL contraction_tensor(nxx, nyy, nzz, nxy, nyz, nzx, nxx, nyy, nzz, nxy, nyz, nzx, nn)
  CALL contraction_tensor(mxx, myy, mzz, mxy, myz, mzx, nxx, nyy, nzz, nxy, nyz, nzx, mn)

  CALL update_boundary(lm)
  CALL update_boundary(ln)
  CALL update_boundary(ll)
  CALL update_boundary(mm)
  CALL update_boundary(nn)
  CALL update_boundary(mn)

  CALL test_filter(lm, lm)
  CALL test_filter(ln, ln)
  CALL test_filter(ll, ll)
  CALL test_filter(mm, mm)
  CALL test_filter(nn, nn)
  CALL test_filter(mn, mn)


  IF (subgrid_isotropic) THEN
     DO k=1, ksize
     DO j=1, jsize
     DO i=1, isize

        IF (mm(i,j,k) * delta2(i,j)**2 <= eps * ll(i,j,k)) THEN
           nu1(i,j,k) = 0.0
           nu2(i,j,k) = 0.0

           flag(i,j,k) = 0.0
        ELSE
           nu1(i,j,k) = lm(i,j,k)/(2.0*mm(i,j,k))
           nu2(i,j,k) = 0.0

           flag(i,j,k) = 1.0
        ENDIF
     END DO
     END DO
     END DO
  ELSE
     DO k=1, ksize
     DO j=1, jsize
     DO i=1, isize

        det(i,j,k) = mm(i,j,k)*nn(i,j,k) - mn(i,j,k)**2

        IF (mm(i,j,k) * delta2(i,j)**2 > eps * ll(i,j,k)) THEN
           IF (nn(i,j,k) * delta2(i,j)**2 > eps * ll(i,j,k)) THEN
              IF (det(i,j,k) > eps * mm(i,j,k)*nn(i,j,k)) THEN
                 nu1(i,j,k) = (lm(i,j,k)*nn(i,j,k) - ln(i,j,k)*mn(i,j,k)) /det(i,j,k)
                 nu2(i,j,k) = (ln(i,j,k)*mm(i,j,k) - lm(i,j,k)*mn(i,j,k)) /det(i,j,k)

                 flag(i,j,k) = 3.0
              ELSE
                 nu1(i,j,k) = lm(i,j,k) /(2.0*mm(i,j,K))
                 nu2(i,j,k) = ln(i,j,k) /(2.0*nn(i,j,K))

                 flag(i,j,k) = 2.0
              END IF
           ELSE
              nu1(i,j,k) = lm(i,j,k) / mm(i,j,K)
              nu2(i,j,k) = 0.0

              flag(i,j,k) = 1.0
           END IF
        ELSE
           nu1(i,j,k) = 0.0
           nu2(i,j,k) = 0.0

           flag(i,j,k) = 0.0
        END IF
     END DO
     END DO
     END DO
  END IF

  DO k=1, ksize
  DO j=1, jsize
  DO i=1, isize
     nu1(i,j,k) = nu1(i,j,k)*srate_grid(i,j,k)
     nu2(i,j,k) = nu2(i,j,k)*nfreq2_grid(i,j,k)
  END DO
  END DO
  END DO

#ifdef DLES_LIMITER
  DO k=1, ksize
  DO j=1, jsize
  DO i=1, isize
     IF (nu1(i,j,k) < 0.0)        nu1(i,j,k) = 0.0
     IF (nu2(i,j,k) < - nu1(i,j,k)) nu2(i,j,k) = -nu1(i,j,k)
  END DO
  END DO
  END DO
#endif

  CALL update_boundary(nu1)
  CALL update_boundary(nu2)

END SUBROUTINE dles_viscosity

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE dles_viscosity_iso(u_grid, v_grid, w_grid, u_test, v_test, w_test, nu, flag)

  REAL(8), INTENT(IN) :: u_grid( -slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8), INTENT(IN) :: v_grid(1-slv:isize+slv,  -slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8), INTENT(IN) :: w_grid(1-slv:isize+slv, 1-slv:jsize+slv,  -slv:ksize+slv)

  REAL(8), INTENT(IN) :: u_test( -slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8), INTENT(IN) :: v_test(1-slv:isize+slv,  -slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8), INTENT(IN) :: w_test(1-slv:isize+slv, 1-slv:jsize+slv,  -slv:ksize+slv)

  REAL(8), INTENT(OUT) :: nu(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)

  REAL(8), OPTIONAL, INTENT(OUT) :: flag(1:isize, 1:jsize, 1:ksize)

  REAL(8) :: srate_grid(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: srate_test(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)

  REAL(8) :: sr_xx(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: sr_yy(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: sr_zz(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: sr_xy( -slv:isize+slv,  -slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: sr_yz(1-slv:isize+slv,  -slv:jsize+slv,  -slv:ksize+slv)
  REAL(8) :: sr_zx( -slv:isize+slv, 1-slv:jsize+slv,  -slv:ksize+slv)

  REAL(8) :: lxx(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: lyy(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: lzz(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: lxy( -slv:isize+slv,  -slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: lyz(1-slv:isize+slv,  -slv:jsize+slv,  -slv:ksize+slv)
  REAL(8) :: lzx( -slv:isize+slv, 1-slv:jsize+slv,  -slv:ksize+slv)

  REAL(8) :: mxx(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: myy(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: mzz(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: mxy( -slv:isize+slv,  -slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: myz(1-slv:isize+slv,  -slv:jsize+slv,  -slv:ksize+slv)
  REAL(8) :: mzx( -slv:isize+slv, 1-slv:jsize+slv,  -slv:ksize+slv)

  REAL(8) :: lm(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: mm(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)

  REAL(8) :: tr(1:isize, 1:jsize, 1:ksize)

  REAL(8) :: tmp_xx(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: tmp_yy(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: tmp_zz(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: tmp_xy( -slv:isize+slv,  -slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: tmp_yz(1-slv:isize+slv,  -slv:jsize+slv,  -slv:ksize+slv)
  REAL(8) :: tmp_zx( -slv:isize+slv, 1-slv:jsize+slv,  -slv:ksize+slv)

  INTEGER :: i, j, k

  tmp_xx(:,:,:) = 0.0
  tmp_yy(:,:,:) = 0.0
  tmp_zz(:,:,:) = 0.0
  tmp_xy(:,:,:) = 0.0
  tmp_yz(:,:,:) = 0.0
  tmp_zx(:,:,:) = 0.0

  DO k=0, ksize+1
  DO j=0, jsize+1
  DO i=0, isize+1
     tmp_xx(i,j,k) = ((u_grid(i-1,j,k)+u_grid(i,j,k))*0.5)**2
     tmp_yy(i,j,k) = ((v_grid(i,j-1,k)+v_grid(i,j,k))*0.5)**2
     tmp_zz(i,j,k) = ((w_grid(i,j,k-1)+w_grid(i,j,k))*0.5)**2
  END DO
  END DO
  END DO

  DO k= 0, ksize+1
  DO j=-1, jsize+1
  DO i=-1, isize+1
     tmp_xy(i,j,k) = (u_grid(i,j,k)+u_grid(i,j+1,k))*0.5 * (v_grid(i,j,k)+v_grid(i+1,j,k))*0.5
  END DO
  END DO
  END DO

  DO k=-1, ksize+1
  DO j=-1, jsize+1
  DO i= 0, isize+1
     tmp_yz(i,j,k) = (v_grid(i,j,k)+v_grid(i,j,k+1))*0.5 * (w_grid(i,j,k)+w_grid(i,j+1,k))*0.5
  END DO
  END DO
  END DO

  DO k=-1, ksize+1
  DO j= 0, jsize+1
  DO i=-1, isize+1
     tmp_zx(i,j,k) = (w_grid(i,j,k)+w_grid(i+1,j,k))*0.5 * (u_grid(i,j,k)+u_grid(i,j,k+1))*0.5
  END DO
  END DO
  END DO

  CALL test_filter(tmp_xx, lxx)
  CALL test_filter(tmp_yy, lyy)
  CALL test_filter(tmp_zz, lzz)
  CALL test_filter_xy(tmp_xy, lxy)
  CALL test_filter_yz(tmp_yz, lyz)
  CALL test_filter_zx(tmp_zx, lzx)

  DO k=1, ksize
  DO j=1, jsize
  DO i=1, isize
     lxx(i,j,k) = lxx(i,j,k) - ((u_test(i-1,j,k)+u_test(i,j,k))*0.5)**2
     lyy(i,j,k) = lyy(i,j,k) - ((v_test(i,j-1,k)+v_test(i,j,k))*0.5)**2
     lzz(i,j,k) = lzz(i,j,k) - ((w_test(i,j,k-1)+w_test(i,j,k))*0.5)**2

     tr(i,j,k) = lxx(i,j,k)+lyy(i,j,k)+lzz(i,j,k)
     lxx(i,j,k) = lxx(i,j,k) - tr(i,j,k)/3.0
     lyy(i,j,k) = lyy(i,j,k) - tr(i,j,k)/3.0
     lzz(i,j,k) = lzz(i,j,k) - tr(i,j,k)/3.0
  END DO
  END DO
  END DO

  DO k=1, ksize
  DO j=0, jsize
  DO i=0, isize
     lxy(i,j,k) = lxy(i,j,k) - (u_test(i,j,k)+u_test(i,j+1,k))*0.5 * (v_test(i,j,k)+v_test(i+1,j,k))*0.5
  END DO
  END DO
  END DO

  DO k=0, ksize
  DO j=0, jsize
  DO i=1, isize
     lyz(i,j,k) = lyz(i,j,k) - (v_test(i,j,k)+v_test(i,j,k+1))*0.5 * (w_test(i,j,k)+w_test(i,j+1,k))*0.5
  END DO
  END DO
  END DO

  DO k=0, ksize
  DO j=1, jsize
  DO i=0, isize
     lzx(i,j,k) = lzx(i,j,k) - (w_test(i,j,k)+w_test(i+1,j,k))*0.5 * (u_test(i,j,k)+u_test(i,j,k+1))*0.5
  END DO
  END DO
  END DO

  CALL sr_tensor(u_grid, v_grid, w_grid, sr_xx, sr_yy, sr_zz, sr_xy, sr_yz, sr_zx)


  srate_grid(:,:,:) = 0.0

  DO k=1, ksize
  DO j=1, jsize
  DO i=1, isize
     srate_grid(i,j,k) = sqrt(2.0*(sr_xx(i,j,k)**2 + sr_yy(i,j,k)**2 + sr_zz(i,j,k)**2) &
          + 4.0*( ((sr_xy(i-1,j-1,k)+sr_xy(i,j-1,k)+sr_xy(i-1,j,k)+sr_xy(i,j,k))*0.25)**2 &
                  + ((sr_yz(i,j-1,k-1)+sr_yz(i,j,k-1)+sr_yz(i,j-1,k)+sr_yz(i,j,k))*0.25)**2 &
                  + ((sr_zx(i-1,j,k-1)+sr_zx(i-1,j,k)+sr_zx(i,j,k-1)+sr_zx(i,j,k))*0.25)**2))
  END DO
  END DO
  END DO

  CALL update_boundary(srate_grid)

  DO k=0, ksize+1
  DO j=0, jsize+1
  DO i=0, isize+1
     tmp_xx(i,j,k) = sr_xx(i,j,k) * srate_grid(i,j,k)
     tmp_yy(i,j,k) = sr_yy(i,j,k) * srate_grid(i,j,k)
     tmp_zz(i,j,k) = sr_zz(i,j,k) * srate_grid(i,j,k)
  END DO
  END DO
  END DO

  DO k= 0, ksize+1
  DO j=-1, jsize+1
  DO i=-1, isize+1
     tmp_xy(i,j,k) = sr_xy(i,j,k) &
          * (srate_grid(i,j,k)+srate_grid(i+1,j,k)+srate_grid(i,j+1,k)+srate_grid(i+1,j+1,k))*0.25
  END DO
  END DO
  END DO

  DO k=-1, ksize+1
  DO j=-1, jsize+1
  DO i= 0, isize+1
     tmp_yz(i,j,k) = sr_yz(i,j,k) &
          * (srate_grid(i,j,k)+srate_grid(i,j+1,k)+srate_grid(i,j,k+1)+srate_grid(i,j+1,k+1))*0.25
  END DO
  END DO
  END DO

  DO k=-1, ksize+1
  DO j= 0, jsize+1
  DO i=-1, isize+1
     tmp_zx(i,j,k) = sr_zx(i,j,k) &
          * (srate_grid(i,j,k)+srate_grid(i,j,k+1)+srate_grid(i+1,j,k)+srate_grid(i+1,j,k+1))*0.25
  END DO
  END DO
  END DO

  CALL test_filter(tmp_xx, mxx)
  CALL test_filter(tmp_yy, myy)
  CALL test_filter(tmp_zz, mzz)
  CALL test_filter_xy(tmp_xy, mxy)
  CALL test_filter_yz(tmp_yz, myz)
  CALL test_filter_zx(tmp_zx, mzx)

  CALL sr_tensor(u_test, v_test, w_test, sr_xx, sr_yy, sr_zz, sr_xy, sr_yz, sr_zx)

  DO k=0, ksize+1
  DO j=0, jsize+1
  DO i=0, isize+1
     srate_test(i,j,k) = sqrt(2.0*(sr_xx(i,j,k)**2 + sr_yy(i,j,k)**2 + sr_zz(i,j,k)**2) &
          + 4.0*( ((sr_xy(i-1,j-1,k)+sr_xy(i,j-1,k)+sr_xy(i-1,j,k)+sr_xy(i,j,k))*0.25)**2 &
                  + ((sr_yz(i,j-1,k-1)+sr_yz(i,j,k-1)+sr_yz(i,j-1,k)+sr_yz(i,j,k))*0.25)**2 &
                  + ((sr_zx(i-1,j,k-1)+sr_zx(i-1,j,k)+sr_zx(i,j,k-1)+sr_zx(i,j,k))*0.25)**2))
  END DO
  END DO
  END DO

  DO k=1, ksize
  DO j=1, jsize
  DO i=1, isize
     mxx(i,j,k) = mxx(i,j,k) - tg_ratio2 * sr_xx(i,j,k) * srate_test(i,j,k)
     myy(i,j,k) = myy(i,j,k) - tg_ratio2 * sr_yy(i,j,k) * srate_test(i,j,k)
     mzz(i,j,k) = mzz(i,j,k) - tg_ratio2 * sr_zz(i,j,k) * srate_test(i,j,k)
  END DO
  END DO
  END DO

  DO k=1, ksize
  DO j=0, jsize
  DO i=0, isize
     mxy(i,j,k) = mxy(i,j,k) - tg_ratio2 * sr_xy(i,j,k) &
          * (srate_test(i,j,k)+srate_test(i+1,j,k)+srate_test(i,j+1,k)+srate_test(i+1,j+1,k))*0.25
  END DO
  END DO
  END DO

  DO k=0, ksize
  DO j=0, jsize
  DO i=1, isize
     myz(i,j,k) = myz(i,j,k) - tg_ratio2 * sr_yz(i,j,k) &
          * (srate_test(i,j,k)+srate_test(i,j+1,k)+srate_test(i,j,k+1)+srate_test(i,j+1,k+1))*0.25
  END DO
  END DO
  END DO

  DO k=0, ksize
  DO j=1, jsize
  DO i=0, isize
     mzx(i,j,k) = mzx(i,j,k) - tg_ratio2 * sr_zx(i,j,k) &
          * (srate_test(i,j,k)+srate_test(i,j,k+1)+srate_test(i+1,j,k)+srate_test(i+1,j,k+1))*0.25
  END DO
  END DO
  END DO

  lm = 0.0
  mm = 0.0

  CALL contraction_tensor(lxx, lyy, lzz, lxy, lyz, lzx, mxx, myy, mzz, mxy, myz, mzx, lm)
  CALL contraction_tensor(mxx, myy, mzz, mxy, myz, mzx, mxx, myy, mzz, mxy, myz, mzx, mm)

!  CALL update_boundary(lm)
!  CALL update_boundary(mm)
!  CALL test_filter(lm, lm)
!  CALL test_filter(mm, mm)

  DO k=1, ksize
  DO j=1, jsize
  DO i=1, isize

     IF ( abs(lm(i,j,k)) >= 2 * c_max * delta2(i,j) * mm(i,j,k)) THEN
        nu(i,j,k) = sign(1.0D0, lm(i,j,k)) * c_max * delta2(i,j) * srate_grid(i,j,k)

!        flag(i,j,k) =  0.0D0
        flag(i,j,k) =  sign(1.0D0, lm(i,j,k)) * c_max
     ELSE
        nu(i,j,k) = lm(i,j,k)/(2.0*mm(i,j,k)) * srate_grid(i,j,k)

!        flag(i,j,k) = 1.0D0
        flag(i,j,k) = lm(i,j,k)/(2.0*mm(i,j,k))/delta2(i,j)
     ENDIF

#ifdef DLES_LIMITER
     IF (nu(i,j,k) < 0.0) nu(i,j,k) = 0.0
#endif
  END DO
  END DO
  END DO

  CALL update_boundary(nu)
  CALL test_filter(nu, nu)
  CALL update_boundary(nu)

END SUBROUTINE dles_viscosity_iso

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE dles_diffusivity_iso(a_grid, u_grid, v_grid, w_grid, u_test, v_test, w_test, &
                            srate_grid, srate_test, kappa, flag)
  REAL(8), INTENT(IN) :: a_grid(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)

  REAL(8), INTENT(IN) :: u_grid( -slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8), INTENT(IN) :: v_grid(1-slv:isize+slv,  -slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8), INTENT(IN) :: w_grid(1-slv:isize+slv, 1-slv:jsize+slv,  -slv:ksize+slv)

  REAL(8), INTENT(IN) :: u_test( -slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8), INTENT(IN) :: v_test(1-slv:isize+slv,  -slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8), INTENT(IN) :: w_test(1-slv:isize+slv, 1-slv:jsize+slv,  -slv:ksize+slv)

  REAL(8), INTENT(IN) :: srate_grid(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8), INTENT(IN) :: srate_test(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)

  REAL(8), INTENT(OUT) :: kappa(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)

  REAL(8), INTENT(OUT) :: flag(1:isize, 1:jsize, 1:ksize)

  REAL(8) :: a_test(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)

  REAL(8) :: lx( -slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: ly(1-slv:isize+slv,  -slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: lz(1-slv:isize+slv, 1-slv:jsize+slv,  -slv:ksize+slv)

  REAL(8) :: mx( -slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: my(1-slv:isize+slv,  -slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: mz(1-slv:isize+slv, 1-slv:jsize+slv,  -slv:ksize+slv)

  REAL(8) :: lm(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: mm(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)

  INTEGER :: i, j, k

  CALL test_filter(a_grid, a_test)
  CALL update_boundary(a_test)

  lx=0.0D0
  ly=0.0D0
  lz=0.0D0

  DO k= 0, ksize+1
  DO j= 0, jsize+1
  DO i=-1, isize+1
     lx(i,j,k) = u_grid(i,j,k) * (a_grid(i,j,k)+a_grid(i+1,j,k))*0.5D0
  END DO
  END DO
  END DO

  DO k= 0, ksize+1
  DO j=-1, jsize+1
  DO i= 0, isize+1
     ly(i,j,k) = v_grid(i,j,k) * (a_grid(i,j,k)+a_grid(i,j+1,k))*0.5D0
  END DO
  END DO
  END DO

  DO k=-1, ksize+1
  DO j= 0, jsize+1
  DO i= 0, isize+1
     lz(i,j,k) = w_grid(i,j,k) * (a_grid(i,j,k)+a_grid(i,j,k+1))*0.5D0
  END DO
  END DO
  END DO

  CALL test_filter_x(lx, lx)
  CALL test_filter_y(ly, ly)
  CALL test_filter_z(lz, lz)

  DO k=1, ksize
  DO j=1, jsize
  DO i=0, isize
     lx(i,j,k) = imask3d(i,j,k)*imask3d(i+1,j,k) * (lx(i,j,k) - u_test(i,j,k) * (a_test(i,j,k)+a_test(i+1,j,k))*0.5D0)
  END DO
  END DO
  END DO

  DO k=1, ksize
  DO j=0, jsize
  DO i=1, isize
     ly(i,j,k) = imask3d(i,j,k)*imask3d(i,j+1,k) * (ly(i,j,k) - v_test(i,j,k) * (a_test(i,j,k)+a_test(i,j+1,k))*0.5D0)
  END DO
  END DO
  END DO

  DO k=0, ksize
  DO j=1, jsize
  DO i=1, isize
     lz(i,j,k) = imask3d(i,j,k)*imask3d(i,j,k+1) * (lz(i,j,k) - w_test(i,j,k) * (a_test(i,j,k)+a_test(i,j,k+1))*0.5D0)
  END DO
  END DO
  END DO

  mx=0.0D0
  my=0.0D0
  mz=0.0D0

  DO k= 0, ksize+1
  DO j= 0, jsize+1
  DO i=-1, isize+1
     mx(i,j,k) = imask3d(i,j,k)*imask3d(i+1,j,k) * (a_grid(i+1,j,k)-a_grid(i,j,k))*idx1(i,j) &
                                                 * (srate_grid(i,j,k)+srate_grid(i+1,j,k))*0.5D0

  END DO
  END DO
  END DO

  DO k= 0, ksize+1
  DO j=-1, jsize+1
  DO i= 0, isize+1
     my(i,j,k) = imask3d(i,j,k)*imask3d(i,j+1,k) * (a_grid(i,j+1,k)-a_grid(i,j,k))*idy1(i,j) &
                                                 * (srate_grid(i,j,k)+srate_grid(i,j+1,k))*0.5D0
  END DO
  END DO
  END DO

  DO k=-1, ksize+1
  DO j= 0, jsize+1
  DO i= 0, isize+1
     mz(i,j,k) = imask3d(i,j,k)*imask3d(i,j,k+1) * (a_grid(i,j,k+1)-a_grid(i,j,k))*idz1(k) &
                                                 * (srate_grid(i,j,k)+srate_grid(i,j,k+1))*0.5D0
  END DO
  END DO
  END DO

  CALL test_filter_x(mx, mx)
  CALL test_filter_y(my, my)
  CALL test_filter_z(mz, mz)

  DO k=1, ksize
  DO j=1, jsize
  DO i=0, isize
     mx(i,j,k) = imask3d(i,j,k)*imask3d(i+1,j,k) &
          * (mx(i,j,k) - tg_ratio2*(a_test(i+1,j,k)-a_test(i,j,k))*idx1(i,j) * (srate_test(i,j,k)+srate_test(i+1,j,k))*0.5D0)
  END DO
  END DO
  END DO


  DO k=1, ksize
  DO j=0, jsize
  DO i=1, isize
     my(i,j,k) = imask3d(i,j,k)*imask3d(i,j+1,k) * (my(i,j,k) - tg_ratio2*(a_test(i,j+1,k)-a_test(i,j,k))*idy1(i,j) &
                                                                   *(srate_test(i,j,k)+srate_test(i,j+1,k))*0.5D0)

  END DO
  END DO
  END DO

  DO k=0, ksize
  DO j=1, jsize
  DO i=1, isize
     mz(i,j,k) = imask3d(i,j,k)*imask3d(i,j,k+1) * (mz(i,j,k) - tg_ratio2*(a_test(i,j,k+1)-a_test(i,j,k))*idz1(k) &
                                                                   *(srate_test(i,j,k)+srate_test(i,j,k+1))*0.5D0)
  END DO
  END DO
  END DO

  lm = 0.0D0
  mm = 0.0D0

  CALL contraction_vector(lx, ly, lz, mx, my, mz, lm)
  CALL contraction_vector(mx, my, mz, mx, my, mz, mm)

  ! CALL update_boundary(lm)
  ! CALL update_boundary(mm)

  ! CALL test_filter(lm, lm)
  ! CALL test_filter(mm, mm)

  DO k=1, ksize
  DO j=1, jsize
  DO i=1, isize

     IF ( abs(lm(i,j,k)) >= c_max * delta2(i,j) * mm(i,j,k)) THEN
        kappa(i,j,k) = sign(1.0D0, lm(i,j,k)) * c_max* delta2(i,j) * srate_grid(i,j,k)

!        flag(i,j,k) = 0.0D0
        flag(i,j,k) = sign(1.0D0, lm(i,j,k)) * c_max
     ELSE
        kappa(i,j,k) = lm(i,j,k)/mm(i,j,k) * srate_grid(i,j,k)

!        flag(i,j,k) = 1.0D0
        flag(i,j,k) = lm(i,j,k)/mm(i,j,k)/delta2(i,j,k)
     ENDIF

#ifdef DLES_LIMITER
     IF (kappa(i,j,k) < 0.0D0) kappa(i,j,k) = 0.0D0
#endif

  END DO
  END DO
  END DO

  CALL update_boundary(kappa)
  CALL test_filter(kappa, kappa)
  CALL update_boundary(kappa)

END SUBROUTINE dles_diffusivity_iso

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE dles_diffusivity_hv(a_grid, u_grid, v_grid, w_grid, u_test, v_test, w_test, &
                               srate_grid, srate_test, nfreq2_grid, nfreq2_test, kappa1, kappa2, flag)
  REAL(8), INTENT(IN) :: a_grid(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)

  REAL(8), INTENT(IN) :: u_grid( -slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8), INTENT(IN) :: v_grid(1-slv:isize+slv,  -slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8), INTENT(IN) :: w_grid(1-slv:isize+slv, 1-slv:jsize+slv,  -slv:ksize+slv)

  REAL(8), INTENT(IN) :: u_test( -slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8), INTENT(IN) :: v_test(1-slv:isize+slv,  -slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8), INTENT(IN) :: w_test(1-slv:isize+slv, 1-slv:jsize+slv,  -slv:ksize+slv)

  REAL(8), INTENT(IN) :: srate_grid(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8), INTENT(IN) :: srate_test(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)

  REAL(8), INTENT(IN) :: nfreq2_grid(1-slv:isize+slv, 1-slv:jsize+slv, -slv:ksize+slv)
  REAL(8), INTENT(IN) :: nfreq2_test(1-slv:isize+slv, 1-slv:jsize+slv, -slv:ksize+slv)

  REAL(8), INTENT(OUT) :: kappa1(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8), INTENT(OUT) :: kappa2(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)

  REAL(8), INTENT(OUT) :: flag(1:isize, 1:jsize, 1:ksize)

  REAL(8) :: a_test(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)

  REAL(8) :: lx( -slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: ly(1-slv:isize+slv,  -slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: lz(1-slv:isize+slv, 1-slv:jsize+slv,  -slv:ksize+slv)

  REAL(8) :: mx( -slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: my(1-slv:isize+slv,  -slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: mz(1-slv:isize+slv, 1-slv:jsize+slv,  -slv:ksize+slv)

  REAL(8) :: nx( -slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: ny(1-slv:isize+slv,  -slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: nz(1-slv:isize+slv, 1-slv:jsize+slv,  -slv:ksize+slv)

  REAL(8) :: dnx( -slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: dny(1-slv:isize+slv,  -slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: dnz(1-slv:isize+slv, 1-slv:jsize+slv,  -slv:ksize+slv)

  REAL(8) :: lm(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: ln(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: ll(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: mn(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: mm(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: nn(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)

  REAL(8) :: det(1:isize, 1:jsize, 1:ksize)

  INTEGER :: i, j, k

  CALL test_filter(a_grid, a_test)
  CALL update_boundary(a_test)

  lx=0.0D0
  ly=0.0D0
  lz=0.0D0

  DO k= 0, ksize+1
  DO j= 0, jsize+1
  DO i=-1, isize+1
     lx(i,j,k) = u_grid(i,j,k) * (a_grid(i,j,k)+a_grid(i+1,j,k))*0.5D0
  END DO
  END DO
  END DO

  DO k= 0, ksize+1
  DO j=-1, jsize+1
  DO i= 0, isize+1
     ly(i,j,k) = v_grid(i,j,k) * (a_grid(i,j,k)+a_grid(i,j+1,k))*0.5D0
  END DO
  END DO
  END DO

  DO k=-1, ksize+1
  DO j= 0, jsize+1
  DO i= 0, isize+1
     lz(i,j,k) = w_grid(i,j,k) * (a_grid(i,j,k)+a_grid(i,j,k+1))*0.5D0
  END DO
  END DO
  END DO

  CALL test_filter_x(lx, lx)
  CALL test_filter_y(ly, ly)
  CALL test_filter_z(lz, lz)

  DO k=1, ksize
  DO j=1, jsize
  DO i=0, isize
     lx(i,j,k) = imask3d(i,j,k)*imask3d(i+1,j,k) * (lx(i,j,k) - u_test(i,j,k) * (a_test(i,j,k)+a_test(i+1,j,k))*0.5D0)
  END DO
  END DO
  END DO

  DO k=1, ksize
  DO j=0, jsize
  DO i=1, isize
     ly(i,j,k) = imask3d(i,j,k)*imask3d(i,j+1,k) * (ly(i,j,k) - v_test(i,j,k) * (a_test(i,j,k)+a_test(i,j+1,k))*0.5D0)
  END DO
  END DO
  END DO

  DO k=0, ksize
  DO j=1, jsize
  DO i=1, isize
     lz(i,j,k) = imask3d(i,j,k)*imask3d(i,j,k+1) * (lz(i,j,k) - w_test(i,j,k) * (a_test(i,j,k)+a_test(i,j,k+1))*0.5D0)
  END DO
  END DO
  END DO

  mx=0.0D0
  my=0.0D0
  mz=0.0D0

  DO k= 0, ksize+1
  DO j= 0, jsize+1
  DO i=-1, isize+1
     mx(i,j,k) = imask3d(i,j,k)*imask3d(i+1,j,k) * (a_grid(i+1,j,k)-a_grid(i,j,k))*idx1(i,j) &
                                           * (srate_grid(i,j,k)+srate_grid(i+1,j,k))*0.5D0

  END DO
  END DO
  END DO

  DO k= 0, ksize+1
  DO j=-1, jsize+1
  DO i= 0, isize+1
     my(i,j,k) = imask3d(i,j,k)*imask3d(i,j+1,k) * (a_grid(i,j+1,k)-a_grid(i,j,k))*idy1(i,j) &
                                           * (srate_grid(i,j,k)+srate_grid(i,j+1,k))*0.5D0
  END DO
  END DO
  END DO

  DO k=-1, ksize+1
  DO j= 0, jsize+1
  DO i= 0, isize+1
     mz(i,j,k) = imask3d(i,j,k)*imask3d(i,j,k+1) * (a_grid(i,j,k+1)-a_grid(i,j,k))*idz1(k) &
                                           * (srate_grid(i,j,k)+srate_grid(i,j,k+1))*0.5D0
  END DO
  END DO
  END DO

  CALL test_filter_x(mx, mx)
  CALL test_filter_y(my, my)
  CALL test_filter_z(mz, mz)

  DO k=1, ksize
  DO j=1, jsize
  DO i=0, isize
     mx(i,j,k) = imask3d(i,j,k)*imask3d(i+1,j,k) * (mx(i,j,k) - tg_ratio2*(a_test(i+1,j,k)-a_test(i,j,k))*idx1(i,j) &
                                                                   *(srate_test(i,j,k)+srate_test(i+1,j,k))*0.5D0)
  END DO
  END DO
  END DO


  DO k=1, ksize
  DO j=0, jsize
  DO i=1, isize
     my(i,j,k) = imask3d(i,j,k)*imask3d(i,j+1,k) * (my(i,j,k) - tg_ratio2*(a_test(i,j+1,k)-a_test(i,j,k))*idy1(i,j) &
                                                                   *(srate_test(i,j,k)+srate_test(i,j+1,k))*0.5D0)

  END DO
  END DO
  END DO

  DO k=0, ksize
  DO j=1, jsize
  DO i=1, isize
     mz(i,j,k) = imask3d(i,j,k)*imask3d(i,j,k+1) * (mz(i,j,k) - tg_ratio2*(a_test(i,j,k+1)-a_test(i,j,k))*idz1(k) &
                                                                   *(srate_test(i,j,k)+srate_test(i,j,k+1))*0.5D0)
  END DO
  END DO
  END DO

  dnx(:,:,:) = 0.0D0
  dny(:,:,:) = 0.0D0
  dnz(:,:,:) = 0.0D0

  DO k=-1, ksize+1
  DO j=-1, jsize+1
  DO i=-1, isize+1
     dnz(i,j,k) = imask3d(i,j,k)*imask3d(i,j,k+1) * (a_grid(i,j,k+1)-a_grid(i,j,k))*idz1(k)
  END DO
  END DO
  END DO

  nx=0.0D0
  ny=0.0D0
  nz=0.0D0

  DO k=-1, ksize+1
  DO j= 0, jsize+1
  DO i= 0, isize+1
     nz(i,j,k) = imask3d(i,j,k)*imask3d(i,j,k+1) * dnz(i,j,k) * sqrt(max(nfreq2_grid(i,j,k), 0.0D0))
  END DO
  END DO
  END DO

  CALL test_filter_z(nz, nz)

  DO k=0, ksize+1
  DO j=0, jsize+1
  DO i=0, isize+1
     dnz(i,j,k) = imask3d(i,j,k)*imask3d(i,j,k+1) * (a_test(i,j,k+1)-a_test(i,j,k))*idz1(k)
  END DO
  END DO
  END DO

  DO k=0, ksize
  DO j=1, jsize
  DO i=1, isize
     nz(i,j,k) = imask3d(i,j,k)*imask3d(i,j,k+1) * (nz(i,j,k) - tg_ratio2 * dnz(i,j,k)*sqrt(max(nfreq2_test(i,j,k),0.0D0)))
  END DO
  END DO
  END DO

  lm = 0.0D0
  ln = 0.0D0
  ll = 0.0D0
  mn = 0.0D0
  mm = 0.0D0
  nn = 0.0D0

  CALL contraction_vector(lx, ly, lz, mx, my, mz, lm)
  CALL contraction_vector(lx, ly, lz, nx, ny, nz, ln)
  CALL contraction_vector(lx, ly, lz, lx, ly, lz, ll)
  CALL contraction_vector(mx, my, mz, nx, ny, nz, mn)
  CALL contraction_vector(mx, my, mz, mx, my, mz, mm)
  CALL contraction_vector(nx, ny, nz, nx, ny, nz, nn)

  DO k=1, ksize
  DO j=1, jsize
  DO i=1, isize

     IF ((nfreq2_grid(i,j,k-1)+nfreq2_grid(i,j,k))*0.5D0 < eps) THEN
        IF (abs(lm(i,j,k)) >= c_max * delta2(i,j) * mm(i,j,k)) THEN
           kappa1(i,j,k) = sign(1.0D0, lm(i,j,k)) * c_max* delta2(i,j) * srate_grid(i,j,k)

           flag(i,j,k) = 1.0
        ELSE
           kappa1(i,j,k) = lm(i,j,k)/mm(i,j,k) * srate_grid(i,j,k)

           flag(i,j,k) = 2.0
        ENDIF

        kappa2(i,j,k) = 0.0D0
     ELSE
        det(i,j,k) = mm(i,j,k)*nn(i,j,k) - mn(i,j,k)**2
        kappa1(i,j,k) = lm(i,j,k)*nn(i,j,k) - ln(i,j,k)*mn(i,j,k)
        kappa2(i,j,k) = ln(i,j,k)*mm(i,j,k) - lm(i,j,k)*mn(i,j,k)

        IF (abs(kappa1(i,j,k)) >= c_max * delta2(i,j) * det(i,j,k)) THEN
           kappa1(i,j,k) = sign(1.0D0, kappa1(i,j,k)) * c_max* delta2(i,j) * srate_grid(i,j,k)
           flag(i,j,k) = 3.0D0
        ELSE
           kappa1(i,j,k) = kappa1(i,j,k)/det(i,j,k) * srate_grid(i,j,k)
           flag(i,j,k) = 4.0D0
        END IF

        IF (abs(kappa2(i,j,k)) >= c_max * delta2(i,j) * det(i,j,k)) THEN
           kappa2(i,j,k) = sign(1.0D0, kappa2(i,j,k)) * c_max* delta2(i,j) * sqrt((nfreq2_grid(i,j,k-1)+nfreq2_grid(i,j,k))*0.5D0)
           flag(i,j,k) = flag(i,j,k) + 2.0
        ELSE
           kappa2(i,j,k) = kappa2(i,j,k)/det(i,j,k) * sqrt((nfreq2_grid(i,j,k-1)+nfreq2_grid(i,j,k))*0.5D0)
        END IF
     END IF

  END DO
  END DO
  END DO

  CALL update_boundary(kappa1)
  CALL update_boundary(kappa2)
  CALL test_filter(kappa1, kappa1)
  CALL test_filter(kappa2, kappa2)
  CALL update_boundary(kappa1)
  CALL update_boundary(kappa2)


END SUBROUTINE dles_diffusivity_hv

!-----------------------------------------------------------------------------------------------------------------------


SUBROUTINE dles_diffusivity(a_grid, u_grid, v_grid, w_grid, u_test, v_test, w_test, &
                            srate_grid, srate_test, nfreq2_grid, nfreq2_test, n_x, n_y, n_z, kappa1, kappa2, flag)
  REAL(8), INTENT(IN) :: a_grid(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)

  REAL(8), INTENT(IN) :: u_grid( -slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8), INTENT(IN) :: v_grid(1-slv:isize+slv,  -slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8), INTENT(IN) :: w_grid(1-slv:isize+slv, 1-slv:jsize+slv,  -slv:ksize+slv)

  REAL(8), INTENT(IN) :: u_test( -slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8), INTENT(IN) :: v_test(1-slv:isize+slv,  -slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8), INTENT(IN) :: w_test(1-slv:isize+slv, 1-slv:jsize+slv,  -slv:ksize+slv)

  REAL(8), INTENT(IN) :: srate_grid(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8), INTENT(IN) :: srate_test(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)

  REAL(8), INTENT(IN) :: nfreq2_grid(1-slv:isize+slv, 1-slv:jsize+slv, -slv:ksize+slv)
  REAL(8), INTENT(IN) :: nfreq2_test(1-slv:isize+slv, 1-slv:jsize+slv, -slv:ksize+slv)

  REAL(8), INTENT(IN) :: n_x( -slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8), INTENT(IN) :: n_y(1-slv:isize+slv,  -slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8), INTENT(IN) :: n_z(1-slv:isize+slv, 1-slv:jsize+slv,  -slv:ksize+slv)

  REAL(8), INTENT(OUT) :: kappa1(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8), INTENT(OUT) :: kappa2(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)

  REAL(8), INTENT(OUT) :: flag(1:isize, 1:jsize, 1:ksize)

  REAL(8) :: a_test(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)

  REAL(8) :: lx( -slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: ly(1-slv:isize+slv,  -slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: lz(1-slv:isize+slv, 1-slv:jsize+slv,  -slv:ksize+slv)

  REAL(8) :: mx( -slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: my(1-slv:isize+slv,  -slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: mz(1-slv:isize+slv, 1-slv:jsize+slv,  -slv:ksize+slv)

  REAL(8) :: nx( -slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: ny(1-slv:isize+slv,  -slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: nz(1-slv:isize+slv, 1-slv:jsize+slv,  -slv:ksize+slv)

  REAL(8) :: dnx( -slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: dny(1-slv:isize+slv,  -slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: dnz(1-slv:isize+slv, 1-slv:jsize+slv,  -slv:ksize+slv)

  REAL(8) :: lm(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: ln(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: ll(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: mn(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: mm(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
  REAL(8) :: nn(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)

  REAL(8) :: det(1:isize, 1:jsize, 1:ksize)

  INTEGER :: i, j, k

  CALL test_filter(a_grid, a_test)
  CALL update_boundary(a_test)

  lx=0.0D0
  ly=0.0D0
  lz=0.0D0

  DO k= 0, ksize+1
  DO j= 0, jsize+1
  DO i=-1, isize+1
     lx(i,j,k) = u_grid(i,j,k) * (a_grid(i,j,k)+a_grid(i+1,j,k))*0.5D0
  END DO
  END DO
  END DO

  DO k= 0, ksize+1
  DO j=-1, jsize+1
  DO i= 0, isize+1
     ly(i,j,k) = v_grid(i,j,k) * (a_grid(i,j,k)+a_grid(i,j+1,k))*0.5D0
  END DO
  END DO
  END DO

  DO k=-1, ksize+1
  DO j= 0, jsize+1
  DO i= 0, isize+1
     lz(i,j,k) = w_grid(i,j,k) * (a_grid(i,j,k)+a_grid(i,j,k+1))*0.5D0
  END DO
  END DO
  END DO

  CALL test_filter_x(lx, lx)
  CALL test_filter_y(ly, ly)
  CALL test_filter_z(lz, lz)

  DO k=1, ksize
  DO j=1, jsize
  DO i=0, isize
     lx(i,j,k) = imask3d(i,j,k)*imask3d(i+1,j,k) * (lx(i,j,k) - u_test(i,j,k) * (a_test(i,j,k)+a_test(i+1,j,k))*0.5D0)
  END DO
  END DO
  END DO

  DO k=1, ksize
  DO j=0, jsize
  DO i=1, isize
     ly(i,j,k) = imask3d(i,j,k)*imask3d(i,j+1,k) * (ly(i,j,k) - v_test(i,j,k) * (a_test(i,j,k)+a_test(i,j+1,k))*0.5D0)
  END DO
  END DO
  END DO

  DO k=0, ksize
  DO j=1, jsize
  DO i=1, isize
     lz(i,j,k) = imask3d(i,j,k)*imask3d(i,j,k+1) * (lz(i,j,k) - w_test(i,j,k) * (a_test(i,j,k)+a_test(i,j,k+1))*0.5D0)
  END DO
  END DO
  END DO

  mx=0.0D0
  my=0.0D0
  mz=0.0D0

  DO k= 0, ksize+1
  DO j= 0, jsize+1
  DO i=-1, isize+1
     mx(i,j,k) = imask3d(i,j,k)*imask3d(i+1,j,k) * (a_grid(i+1,j,k)-a_grid(i,j,k))*idx1(i,j) &
                                           * (srate_grid(i,j,k)+srate_grid(i+1,j,k))*0.5D0

  END DO
  END DO
  END DO

  DO k= 0, ksize+1
  DO j=-1, jsize+1
  DO i= 0, isize+1
     my(i,j,k) = imask3d(i,j,k)*imask3d(i,j+1,k) * (a_grid(i,j+1,k)-a_grid(i,j,k))*idy1(i,j) &
                                           * (srate_grid(i,j,k)+srate_grid(i,j+1,k))*0.5D0
  END DO
  END DO
  END DO

  DO k=-1, ksize+1
  DO j= 0, jsize+1
  DO i= 0, isize+1
     mz(i,j,k) = imask3d(i,j,k)*imask3d(i,j,k+1) * (a_grid(i,j,k+1)-a_grid(i,j,k))*idz1(k) &
                                           * (srate_grid(i,j,k)+srate_grid(i,j,k+1))*0.5D0
  END DO
  END DO
  END DO

  CALL test_filter_x(mx, mx)
  CALL test_filter_y(my, my)
  CALL test_filter_z(mz, mz)

  DO k=1, ksize
  DO j=1, jsize
  DO i=0, isize
     mx(i,j,k) = imask3d(i,j,k)*imask3d(i+1,j,k) * (mx(i,j,k) - tg_ratio2*(a_test(i+1,j,k)-a_test(i,j,k))*idx1(i,j) &
                                                                   *(srate_test(i,j,k)+srate_test(i+1,j,k))*0.5D0)
  END DO
  END DO
  END DO


  DO k=1, ksize
  DO j=0, jsize
  DO i=1, isize
     my(i,j,k) = imask3d(i,j,k)*imask3d(i,j+1,k) * (my(i,j,k) - tg_ratio2*(a_test(i,j+1,k)-a_test(i,j,k))*idy1(i,j) &
                                                                   *(srate_test(i,j,k)+srate_test(i,j+1,k))*0.5D0)

  END DO
  END DO
  END DO

  DO k=0, ksize
  DO j=1, jsize
  DO i=1, isize
     mz(i,j,k) = imask3d(i,j,k)*imask3d(i,j,k+1) * (mz(i,j,k) - tg_ratio2*(a_test(i,j,k+1)-a_test(i,j,k))*idz1(k) &
                                                                   *(srate_test(i,j,k)+srate_test(i,j,k+1))*0.5D0)
  END DO
  END DO
  END DO

  dnx(:,:,:) = 0.0D0
  dny(:,:,:) = 0.0D0
  dnz(:,:,:) = 0.0D0

  DO k=-1, ksize+1
  DO j=-1, jsize+1
  DO i=-1, isize+1
     dnx(i,j,k) = imask3d(i,j,k)*imask3d(i+1,j,k) * n_x(i,j,k) * (a_grid(i+1,j,k)-a_grid(i,j,k))*idx1(i,j)
     dny(i,j,k) = imask3d(i,j,k)*imask3d(i,j+1,k) * n_y(i,j,k) * (a_grid(i,j+1,k)-a_grid(i,j,k))*idy1(i,j)
     dnz(i,j,k) = imask3d(i,j,k)*imask3d(i,j,k+1) * n_z(i,j,k) * (a_grid(i,j,k+1)-a_grid(i,j,k))*idz1(k)
  END DO
  END DO
  END DO

  CALL  update_boundary(dnx)
  CALL  update_boundary(dny)
  CALL  update_boundary(dnz)

  nx=0.0D0
  ny=0.0D0
  nz=0.0D0

  DO k= 0, ksize+1
  DO j= 0, jsize+1
  DO i=-1, isize+1
     nx(i,j,k) = imask3d(i,j,k)*imask3d(i+1,j,k) * n_x(i,j,k) &
          * (dnx(i,j,k) + (dny(i,j-1,k)+dny(i,j,k)+dny(i+1,j-1,k)+dny(i+1,j,k))*0.25D0  &
                        + (dnz(i,j,k-1)+dnz(i,j,k)+dnz(i+1,j,k-1)+dnz(i+1,j,k))*0.25D0) &
          * (nfreq2_grid(i,j,k)+nfreq2_grid(i+1,j,k))*0.5D0
  END DO
  END DO
  END DO

  DO k= 0, ksize+1
  DO j=-1, jsize+1
  DO i= 0, isize+1
     ny(i,j,k) = imask3d(i,j,k)*imask3d(i,j+1,k) * n_y(i,j,k) &
          * (dny(i,j,k) + (dnz(i,j,k-1)+dnz(i,j,k)+dnz(i,j+1,k-1)+dnz(i,j+1,k))*0.25D0  &
                        + (dnx(i-1,j,k)+dnx(i,j,k)+dnx(i-1,j+1,k)+dnx(i,j+1,k))*0.25D0) &
          * (nfreq2_grid(i,j,k)+nfreq2_grid(i,j+1,k))*0.5D0
  END DO
  END DO
  END DO

  DO k=-1, ksize+1
  DO j= 0, jsize+1
  DO i= 0, isize+1
     nz(i,j,k) = imask3d(i,j,k)*imask3d(i,j,k+1) * n_z(i,j,k) &
          * (dnz(i,j,k) + (dnx(i-1,j,k)+dnx(i,j,k)+dnx(i-1,j,k+1)+dnx(i,j,k+1))*0.25D0  &
                        + (dny(i,j-1,k)+dny(i,j,k)+dny(i,j-1,k+1)+dny(i,j,k+1))*0.25D0) &
          * sqrt(max(nfreq2_grid(i,j,k), 0.0D0))
  END DO
  END DO
  END DO

  CALL test_filter_x(nx, nx)
  CALL test_filter_y(ny, ny)
  CALL test_filter_z(nz, nz)

  DO k=0, ksize+1
  DO j=0, jsize+1
  DO i=0, isize+1
     dnx(i,j,k) = imask3d(i,j,k)*imask3d(i+1,j,k) * n_x(i,j,k) * (a_test(i+1,j,k)-a_test(i,j,k))*idx1(i,j)
     dny(i,j,k) = imask3d(i,j,k)*imask3d(i,j+1,k) * n_y(i,j,k) * (a_test(i,j+1,k)-a_test(i,j,k))*idy1(i,j)
     dnz(i,j,k) = imask3d(i,j,k)*imask3d(i,j,k+1) * n_z(i,j,k) * (a_test(i,j,k+1)-a_test(i,j,k))*idz1(k)
  END DO
  END DO
  END DO

  DO k=1, ksize
  DO j=1, jsize
  DO i=0, isize
     nx(i,j,k) = imask3d(i,j,k)*imask3d(i+1,j,k) * (nx(i,j,k) - tg_ratio2 * n_x(i,j,k) &
          * (dnx(i,j,k) + (dny(i,j-1,k)+dny(i,j,k)+dny(i+1,j-1,k)+dny(i+1,j,k))*0.25D0  &
                        + (dnz(i,j,k-1)+dnz(i,j,k)+dnz(i+1,j,k-1)+dnz(i+1,j,k))*0.25D0) &
          * (nfreq2_test(i,j,k)+nfreq2_test(i+1,j,k))*0.5D0)
  END DO
  END DO
  END DO

  DO k=1, ksize
  DO j=0, jsize
  DO i=1, isize
     ny(i,j,k) = imask3d(i,j,k)*imask3d(i,j+1,k) * (ny(i,j,k) - tg_ratio2 * n_y(i,j,k) &
          * (dny(i,j,k) + (dnz(i,j,k-1)+dnz(i,j,k)+dnz(i,j+1,k-1)+dnz(i,j+1,k))*0.25D0  &
                        + (dnx(i-1,j,k)+dnx(i,j,k)+dnx(i-1,j+1,k)+dnx(i,j+1,k))*0.25D0) &
          * (nfreq2_test(i,j,k)+nfreq2_test(i,j+1,k))*0.5D0)
  END DO
  END DO
  END DO

  DO k=0, ksize
  DO j=1, jsize
  DO i=1, isize
     nz(i,j,k) = imask3d(i,j,k)*imask3d(i,j,k+1) * (nz(i,j,k) - tg_ratio2 * n_z(i,j,k) &
          * (dnz(i,j,k) + (dnx(i-1,j,k)+dnx(i,j,k)+dnx(i-1,j,k+1)+dnx(i,j,k+1))*0.25D0  &
                        + (dny(i,j-1,k)+dny(i,j,k)+dny(i,j-1,k+1)+dny(i,j,k+1))*0.25D0) &
          * (sqrt(max(nfreq2_test(i,j,k),0.0D0))))
  END DO
  END DO
  END DO

  lm = 0.0D0
  ln = 0.0D0
  ll = 0.0D0
  mn = 0.0D0
  mm = 0.0D0
  nn = 0.0D0

  CALL contraction_vector(lx, ly, lz, mx, my, mz, lm)
  CALL contraction_vector(lx, ly, lz, nx, ny, nz, ln)
  CALL contraction_vector(lx, ly, lz, lx, ly, lz, ll)
  CALL contraction_vector(mx, my, mz, nx, ny, nz, mn)
  CALL contraction_vector(mx, my, mz, mx, my, mz, mm)
  CALL contraction_vector(nx, ny, nz, nx, ny, nz, nn)

  ! CALL update_boundary(lm)
  ! CALL update_boundary(ln)
  ! CALL update_boundary(ll)
  ! CALL update_boundary(mm)
  ! CALL update_boundary(nn)
  ! CALL update_boundary(mn)

  ! CALL test_filter(lm, lm)
  ! CALL test_filter(ln, ln)
  ! CALL test_filter(ll, ll)
  ! CALL test_filter(mm, mm)
  ! CALL test_filter(nn, nn)
  ! CALL test_filter(mn, mn)

  DO k=1, ksize
  DO j=1, jsize
  DO i=1, isize

     det(i,j,k) = mm(i,j,k)*nn(i,j,k) - mn(i,j,k)**2

     IF (mm(i,j,k) * delta2(i,j)**2 > eps * ll(i,j,k)) THEN
        IF (nn(i,j,k) * delta2(i,j)**2 > eps * ll(i,j,k)) THEN
           IF (det(i,j,k) > eps * mm(i,j,k)*nn(i,j,k)) THEN
              kappa1(i,j,k) = (lm(i,j,k)*nn(i,j,k) - ln(i,j,k)*mn(i,j,k)) /det(i,j,k) * srate_grid(i,j,k)
              kappa2(i,j,k) = (ln(i,j,k)*mm(i,j,k) - lm(i,j,k)*mn(i,j,k)) /det(i,j,k) * sqrt(max((nfreq2_grid(i,j,k-1)+nfreq2_grid(i,j,k))*0.5D0,0.0D0))

              flag(i,j,k) = 3.0D0
           ELSE
              kappa1(i,j,k) = lm(i,j,k) /(2.0*mm(i,j,K)) * srate_grid(i,j,k)
              kappa2(i,j,k) = ln(i,j,k) /(2.0*nn(i,j,K)) * sqrt(max((nfreq2_grid(i,j,k-1)+nfreq2_grid(i,j,k))*0.5D0,0.0D0))

              flag(i,j,k) = 2.0
           END IF
        ELSE
           kappa1(i,j,k) = lm(i,j,k) / mm(i,j,K) * srate_grid(i,j,k)
           kappa2(i,j,k) = 0.0

           flag(i,j,k) = 1.0
        END IF
     ELSE
        kappa1(i,j,k) = 0.0
        kappa2(i,j,k) = 0.0

        flag(i,j,k) = 0.0
     END IF

  END DO
  END DO
  END DO


  CALL update_boundary(kappa1)
  CALL update_boundary(kappa2)

END SUBROUTINE dles_diffusivity

!-----------------------------------------------------------------------------------------------------------------------

END MODULE les
