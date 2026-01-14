#include "macro.h"

MODULE seaice
  USE misc
  USE geometry
  IMPLICIT NONE
  SAVE

  REAL(8) :: ai_max = 0.99D0

  REAL(8), ALLOCATABLE :: ui(:,:)
  REAL(8), ALLOCATABLE :: vi(:,:)
  REAL(8), ALLOCATABLE :: hi(:,:)
  REAL(8), ALLOCATABLE :: ai(:,:)

  REAL(8), ALLOCATABLE :: gix(:,:)
  REAL(8), ALLOCATABLE :: giy(:,:)

  REAL(8), ALLOCATABLE :: gix_ab(:,:)
  REAL(8), ALLOCATABLE :: giy_ab(:,:)

  REAL(8), PARAMETER :: hi_min = 0.01
  REAL(8), PARAMETER :: ai_min = 0.01
CONTAINS

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE init_seaice
    USE io

    ALLOCATE(hi(1-slv:isize+slv,1-slv:jsize+slv))
    ALLOCATE(ai(1-slv:isize+slv,1-slv:jsize+slv))
    ALLOCATE(ui( -slv:isize+slv,1-slv:jsize+slv))
    ALLOCATE(vi(1-slv:isize+slv, -slv:jsize+slv))

    hi(:,:) = 0.0D0
    ai(:,:) = 0.0D0
    ui(:,:) = 0.0D0
    vi(:,:) = 0.0D0

    ALLOCATE(gix(0:isize,1:jsize))
    ALLOCATE(giy(1:isize,0:jsize))
    gix(:,:) = 0.0D0
    giy(:,:) = 0.0D0

    seaice_coupled = .FALSE.

    IF (.NOT. seaice_coupled) RETURN

    CALL assert(kpes==1, "PARALLEL3D (kpes > 2) is not supported for seaice-coupled run")

    CALL initial_data('HI', hi, default=.TRUE.)
    CALL initial_data('AI', ai, default=.TRUE.)
    CALL initial_data('UI', ui, default=.TRUE.)
    CALL initial_data('VI', vi, default=.TRUE.)

    ALLOCATE(gix_ab(0:isize,1:jsize))
    ALLOCATE(giy_ab(1:isize,0:jsize))

    gix_ab(:,:) = 0.0D0
    giy_ab(:,:) = 0.0D0

    CALL coriolis_2d(ui, vi, gix, giy)

    CALL initial_data('GIX', gix_ab, default=perfect_restart)
    CALL initial_data('GIY', giy_ab, default=perfect_restart)

  END SUBROUTINE init_seaice

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE seaice_advection(ui, vi, hi, ai)
    REAL(8), INTENT(IN)    :: ui( -slv:isize+slv, 1-slv:jsize+slv)
    REAL(8), INTENT(IN)    :: vi(1-slv:isize+slv,  -slv:jsize+slv)
    REAL(8), INTENT(INOUT) :: hi(1-slv:isize+slv, 1-slv:jsize+slv)
    REAL(8), INTENT(INOUT) :: ai(1-slv:isize+slv, 1-slv:jsize+slv)

    REAL(8) :: fx(0:isize, 1:jsize)
    REAL(8) :: fy(1:isize, 0:jsize)

    INTEGER :: i, j

    IF (vrank/=0) RETURN

    fx(:,:) = 0.0D0
    fy(:,:) = 0.0D0
    CALL advection_cosmic(hi, ui, vi, dtime, fx, fy)

    DO j=1, jsize
    DO i=1, isize
       hi(i,j) = hi(i,j) + (  fx(i-1,j)*imask3d(i-1,j,ksize)*imask3d(i,j,ksize)*dy1(i-1,j) &
                            - fx(i,  j)*imask3d(i+1,j,ksize)*imask3d(i,j,ksize)*dy1(i,j)   &
                            + fy(i,j-1)*imask3d(i,j-1,ksize)*imask3d(i,j,ksize)*dx1(i,j-1) &
                            - fy(i,j)  *imask3d(i,j+1,ksize)*imask3d(i,j,ksize)*dx1(i,j)) * dtime/dsz(i,j)
    END DO
    END DO


    fx(:,:) = 0.0D0
    fy(:,:) = 0.0D0
    CALL advection_cosmic(ai, ui, vi, dtime, fx, fy)

    DO j=1, jsize
    DO i=1, isize
       ai(i,j) = ai(i,j) + (  fx(i-1,j)*imask3d(i-1,j,ksize)*imask3d(i,j,ksize)*dy1(i-1,j) &
                            - fx(i,  j)*imask3d(i+1,j,ksize)*imask3d(i,j,ksize)*dy1(i,j)   &
                            + fy(i,j-1)*imask3d(i,j-1,ksize)*imask3d(i,j,ksize)*dx1(i,j-1) &
                            - fy(i,j)  *imask3d(i,j+1,ksize)*imask3d(i,j,ksize)*dx1(i,j)) * dtime/dsz(i,j)

       ai(i,j) = min(max(ai(i,j), 0.0D0), ai_max)
    END DO
    END DO

    CALL update_boundary(ai)
    CALL update_boundary(hi)

  END SUBROUTINE seaice_advection

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE advection_cosmic(a, u, v, dtime, fx, fy)
    REAL(8), INTENT(IN) :: a(1-slv:isize+slv, 1-slv:jsize+slv)

    REAL(8), INTENT(IN) :: u( -slv:isize+slv, 1-slv:jsize+slv)
    REAL(8), INTENT(IN) :: v(1-slv:isize+slv,  -slv:jsize+slv)

    REAL(4), INTENT(IN) :: dtime

    REAL(8), INTENT(INOUT) :: fx(0:isize, 1:jsize)
    REAL(8), INTENT(INOUT) :: fy(1:isize, 0:jsize)

    REAL(8) :: a_x( 1:isize,   -1:jsize+2)
    REAL(8) :: a_y(-1:isize+2,  1:jsize)

    REAL(8) :: udt, vdt

    INTEGER :: i, j

    DO j=-1, jsize+2
    DO i= 1, isize
       udt = imask2d(i-1,j)*imask2d(i,j)*imask2d(i+1,j) * (u(i-1,j)+u(i,j)) * dtime / (2.0D0*dx(i,j)+dx(i-1,j)+dx(i+1,j))
       a_x(i,j) = a(i,j) - 0.5D0 * udt * (a(i+1,j)-a(i-1,j))
    END DO
    END DO

    DO j= 1, jsize
    DO i=-1, isize+2
       vdt = imask2d(i,j-1)*imask2d(i,j)*imask2d(i,j+1) * (v(i,j-1)+v(i,j)) * dtime / (2.0D0*dy(i,j)+dy(i,j-1)+dy(i,j+1))
       a_y(i,j) = a(i,j) - 0.5D0 * vdt * (a(i,j+1)-a(i,j-1))
    END DO
    END DO

    DO j=1, jsize
    DO i=0, isize
       fx(i,j) = fx(i,j) + flux_quickest(u(i,j), imask2d(i-1,j)*a_y(i-1,j) + (1-imask2d(i-1,j))*a_y(i,j),   &
                                                 a_y(i,j),   &
                                                 a_y(i+1,j), &
                                                 imask2d(i+2,j)*a_y(i+2,j) + (1-imask2d(i+2,j))*a_y(i+1,j), &
                                         dx(i-1,j), dx(i,j), dx(i+1,j), dx(i+2,j), dtime)
    END DO
    END DO

    DO j=0, jsize
    DO i=1, isize
       fy(i,j) = fy(i,j) + flux_quickest(v(i,j), imask2d(i,j-1)*a_x(i,j-1) + (1-imask2d(i,j-1))*a_x(i,j),   &
                                                 a_x(i,j),   &
                                                 a_x(i,j+1), &
                                                 imask2d(i,j+2)*a_x(i,j+2) + (1-imask2d(i,j+2))*a_x(i,j+1), &
                                         dy(i,j-1), dy(i,j), dy(i,j+1), dy(i,j+2), dtime)
    END DO
    END DO

CONTAINS

!-----------------------------------------------------------------------------------------------------------------------

  REAL(8) PURE FUNCTION flux_quickest(v, a0, a1, a2, a3, d0, d1, d2, d3, dtime)
    REAL(8), INTENT(IN) :: v
    REAL(8), INTENT(IN) :: a0, a1, a2, a3
    REAL(8),        INTENT(IN) :: d0, d1, d2, d3
    REAL(4),        INTENT(IN) :: dtime

    REAL(8) :: au, ac, ad, a_eff
    REAL(8) :: v_abs
    REAL(8)        :: du, dc, dd
    REAL(8) :: c0, c1, c2

    IF (v > 0.0) THEN
       au = a0
       ac = a1
       ad = a2
       du = d0*0.5 + d1
       dc = d1*0.5
       dd = d2*0.5
       v_abs = v
    ELSE
       au = a3
       ac = a2
       ad = a1
       du = d3*0.5 + d2
       dc = d2*0.5
       dd = d1*0.5
       v_abs = -v
    END IF

    c2 = ((au-ac)/(du-dc) + (ad-ac)/(dd+dc))/(du+dd)
    c1 = (ad-ac)/(dc+dd) - c2*(dd-dc)
    c0 = ac - c2*dc**2 + c1*dc

    a_eff = c2 *((v*dtime)**2/3.0 - (dc+dd)**2/12.0) - c1*v_abs*dtime*0.5 + c0

    flux_quickest = v*limiter(a_eff, au, ac, ad, v_abs*dtime/(dc+dd))
  END FUNCTION flux_quickest

  REAL(8) PURE FUNCTION limiter(a, au, ac, ad, alpha)
    REAL(8), INTENT(IN) :: a, au, ac, ad, alpha

    REAL(8) :: ar

    IF (alpha < 1.0E-5) THEN
       limiter = a
    ELSE IF ((ac < au .AND. ac < ad) .OR. (ac > au .AND. ac > ad)) THEN
       limiter = ac
    ELSE
       ar = au + (ac - au)/alpha
       limiter = a
       IF (ad > au) THEN
          IF (a < ac) limiter = ac
          IF (a > ad) limiter = ad
          IF (a > ar) limiter = ar
       ELSE
          IF (a > ac) limiter = ac
          IF (a < ad) limiter = ad
          IF (a < ar) limiter = ar
       END IF
    END IF
  END FUNCTION limiter

  END SUBROUTINE advection_cosmic

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE seaice_rheology(ui, vi, hi, ai, gix, giy)
    USE parameters
    REAL(8), INTENT(IN) :: ui( -slv:isize+slv, 1-slv:jsize+slv)
    REAL(8), INTENT(IN) :: vi(1-slv:isize+slv,  -slv:jsize+slv)
    REAL(8), INTENT(IN) :: hi(1-slv:isize+slv, 1-slv:jsize+slv)
    REAL(8), INTENT(IN) :: ai(1-slv:isize+slv, 1-slv:jsize+slv)
    REAL(8), INTENT(INOUT) :: gix(0:isize, 1:jsize)
    REAL(8), INTENT(INOUT) :: giy(1:isize, 0:jsize)



  END SUBROUTINE seaice_rheology

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE coriolis_2d(u, v, gx, gy)
    USE io

    REAL(8), INTENT(IN) :: u( -slv:isize+slv, 1-slv:jsize+slv)
    REAL(8), INTENT(IN) :: v(1-slv:isize+slv,  -slv:jsize+slv)

    REAL(8), INTENT(INOUT) :: gx(0:isize, 1:jsize)
    REAL(8), INTENT(INOUT) :: gy(1:isize, 0:jsize)

    LOGICAL :: stat
    INTEGER :: i, j

    CALL checkin('CORIOLIS_Z', corz, stat)

!$OMP PARALLEL PRIVATE(i, j)
!$OMP DO
    DO j=1, jsize
    DO i=0, isize
       gx(i,j) = gx(i,j) + imask2d(i,j)*imask2d(i+1,j) &
            * ( corz(i,  j) * (v(i,  j-1)+v(i,  j))*0.5D0 * dsz(i,  j)  &
              + corz(i+1,j) * (v(i+1,j-1)+v(i+1,j))*0.5D0 * dsz(i+1,j)) &
            / (dsz(i,j)+dsz(i+1,j))
    END DO
    END DO

!$OMP DO
    DO j=0, jsize
    DO i=1, isize
       gy(i,j) = gy(i,j) - imask2d(i,j)*imask2d(i,j+1) &
            * ( corz(i,j)   * (u(i-1,j)  +u(i,j))  *0.5D0 * dsz(i,j)    &
              + corz(i,j+1) * (u(i-1,j+1)+u(i,j+1))*0.5D0 * dsz(i,j+1)) &
            / (dsz(i,j)+dsz(i,j+1))
    END DO
    END DO
!$OMP END PARALLEL

  END SUBROUTINE coriolis_2d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE seaice_thermodynamics(t, s, hi, ai)
    USE parameters
    USE tracers, ONLY : TRC_KIND

    REAL(TRC_KIND), INTENT(INOUT) :: t(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
    REAL(TRC_KIND), INTENT(INOUT) :: s(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)

    REAL(8), INTENT(INOUT) :: hi(1-slv:isize+slv, 1-slv:jsize+slv)
    REAL(8), INTENT(INOUT) :: ai(1-slv:isize+slv, 1-slv:jsize+slv)


  END SUBROUTINE seaice_thermodynamics

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE checkout_seaice
    USE io
    CALL checkout('GIX', gix_ab)
    CALL checkout('GIY', gix_ab)
    CALL checkout('UI', ui)
    CALL checkout('VI', vi)
    CALL checkout('HI', hi)
    CALL checkout('AI', ai)
  END SUBROUTINE checkout_seaice

END MODULE
