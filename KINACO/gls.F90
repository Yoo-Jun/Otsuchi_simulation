#include "macro.h"

MODULE gls
  USE misc
  USE geometry
  USE velocity
  USE tracers
  IMPLICIT NONE
  SAVE
  PRIVATE

  PUBLIC init_gls, gls_mixing

  INTEGER :: tracer_index_tke
  INTEGER :: tracer_index_gls

  INTEGER :: gls_scheme = 3       ! 1: k-epsilon model (Rodi  1987)
                                  ! 2: k-omega   model (Wilcox 1988)
                                  ! 3: GEN model (Umlauf and Burchard 2003)
                                  !    parameters following Warner et al. (2005)
                                  !    stability function follows Canuto et al. (2001)/Buarchard and Bolden (2001) CA1

  REAL(4) :: max_visc = 1.0       ! maximum clip viscosity [m^2/s]
  REAL(4) :: max_diff = 1.0       ! maximum clip diffusivity [m^2/s]

  REAL(4) :: min_tke = 1.0E-6     ! minimum tke [m^2/s^2]
  REAL(4) :: min_gls = 1.0E-12    ! minimum gls

  LOGICAL :: report = .FALSE.

  TYPE :: gls_param_struct
     REAL(4) :: p
     REAL(4) :: m
     REAL(4) :: n
     REAL(4) :: c0
     REAL(4) :: c1
     REAL(4) :: c2
     REAL(4) :: c3p
     REAL(4) :: c3m
     REAL(4) :: stke
     REAL(4) :: sgls
     CHARACTER(16) :: name
  END TYPE gls_param_struct

  TYPE(gls_param_struct) :: gls_param

CONTAINS

  SUBROUTINE init_gls
    LOGICAL :: no_advection = .FALSE.

    INTEGER        :: iostat
    CHARACTER(256) :: iomsg

    LOGICAL :: use_turbulence = .FALSE.

    NAMELIST / gls /        &
         use_gls,           &
         gls_scheme,        &
         no_advection,      &
         max_visc,          &
         max_diff,          &
         min_tke,           &
         min_gls,           &
         report

    ! for backward compatibility
    NAMELIST / turbulence / &
         use_turbulence,    &
         gls_scheme,        &
         no_advection,      &
         max_visc,          &
         max_diff,          &
         min_tke,           &
         min_gls,           &
         report

    IF (rank==0) THEN
       REWIND(CONFIG_UNIT)
       READ(CONFIG_UNIT, NML=gls, IOSTAT=iostat, IOMSG=iomsg)
       CALL assert(iostat <= 0, "failed to read GLS namelist", iomsg)

       REWIND(CONFIG_UNIT)
       READ(CONFIG_UNIT, NML=turbulence, IOSTAT=iostat, IOMSG=iomsg)
       CALL assert(iostat <= 0, "failed to read TURBULENCE", iomsg)

       use_gls = use_gls .OR. use_turbulence
    END IF

    CALL bcast(use_gls)
    CALL bcast(gls_scheme)
    CALL bcast(no_advection)
    CALL bcast(max_visc)
    CALL bcast(max_diff)
    CALL bcast(min_tke)
    CALL bcast(min_gls)
    CALL bcast(report)

    IF (.NOT. use_gls) RETURN

    SELECT CASE(gls_scheme)
    CASE (1) !k-epsilon model
       gls_param%name = 'k-epsioon'

       gls_param%p    =  3.0
       gls_param%m    =  3.0/2
       gls_param%n    = -1.0

       gls_param%c0   =  0.5270
       gls_param%c1   =  1.44
       gls_param%c2   =  1.92
       gls_param%c3p  =  1.0
       gls_param%c3m  = -0.63

       gls_param%stke =  1.0
       gls_param%sgls =  1.3

    CASE (2) !k-omega model !MY2.5
       gls_param%name = 'k-omega'

       gls_param%p    = 0.0
       gls_param%m    = 1.0
       gls_param%n    = 1.0

       gls_param%c0   =  0.5544
       gls_param%c1   =  0.9
       gls_param%c2   =  0.52
       gls_param%c3p  =  2.5
       gls_param%c3m  =  1.0

       gls_param%stke =  1.96
       gls_param%sgls =  1.96

    CASE (3) !gen model
       gls_param%name = 'gen'
       gls_param%p    =  2.0
       gls_param%m    =  1.0
       gls_param%n    = -2.0/3

       gls_param%c0   =  0.5270
       gls_param%c1   =  1.0
       gls_param%c2   =  1.22
       gls_param%c3p  =  1.0
       gls_param%c3m  =  0.05

       gls_param%stke =  0.8
       gls_param%sgls =  1.07

      !  gls_param%p    =  2.0
      !  gls_param%m    =  1.0
      !  gls_param%n    = -0.67

      !  gls_param%c0   =  0.5544
      !  gls_param%c1   =  1.0
      !  gls_param%c2   =  1.22
      !  gls_param%c3p  =  1.0
      !  gls_param%c3m  =  0.1

      !  gls_param%stke =  0.8
      !  gls_param%sgls =  1.07

    CASE DEFAULT
       CALL assert(.FALSE., "unsupported GLS_SCHEME")

    END SELECT

    CALL add_tracer('TKE', tracer_index_tke)
    CALL add_tracer('GLS', tracer_index_gls)

    tracer_info(tracer_index_tke)%diffscheme = 0
    tracer_info(tracer_index_gls)%diffscheme = 0

    IF (no_advection) THEN
       tracer_info(tracer_index_tke)%advscheme = 0
       tracer_info(tracer_index_gls)%advscheme = 0
    END IF

    tracer_info(tracer_index_tke)%cgls = -1.0/gls_param%stke
    tracer_info(tracer_index_gls)%cgls = -1.0/gls_param%sgls
    !negative cgls indicates to use KM (visc) instead of KH (diff)

    IF (rank==0) WRITE(REPORT_UNIT, *) "enable GLS turbulence closure model ('"//trim(gls_param%name)//"')"
  END SUBROUTINE init_gls

  SUBROUTINE gls_mixing(visc, diff)
    USE velocity, ONLY: sfreq2
    USE state,    ONLY: bfreq2
    USE io

    REAL(8), INTENT(OUT) :: visc(0:isize+1,0:jsize+1,0:ksize+1)
    REAL(8), INTENT(OUT) :: diff(0:isize+1,0:jsize+1,0:ksize+1)

    REAL(8) :: tlen(0:isize+1,0:jsize+1,0:ksize+1)

    REAL(8) :: vfric_sfc(0:isize+1,0:jsize+1)
    REAL(8) :: vfric_btm(0:isize+1,0:jsize+1)

    REAL(8) :: am, an, sm, sh, km, kh, tmp
    REAL(8) :: gls, tke, sqrttke

    INTEGER :: i, j, k

    REAL(4) :: p, m, n, c0, c1, c2, c3

    p  = gls_param%p
    m  = gls_param%m
    n  = gls_param%n
    c0 = gls_param%c0
    c1 = gls_param%c1
    c2 = gls_param%c2

!$OMP PARALLEL
!$OMP WORKSHARE
    vfric_sfc(:,:) = 0.0
    vfric_btm(:,:) = 0.0
!$OMP END WORKSHARE

!$OMP DO
    DO j=1, jsize
    DO i=1, isize
       vfric_sfc(i,j) = sqrt( sqrt(0.25*(taux_sfc(i-1,j)+taux_sfc(i,j))**2+0.25*(tauy_sfc(i,j-1)+tauy_sfc(i,j))**2)/rho_0)
       vfric_btm(i,j) = sqrt( sqrt(0.25*(taux_btm(i-1,j)+taux_btm(i,j))**2+0.25*(tauy_btm(i,j-1)+tauy_btm(i,j))**2)/rho_0)
    END DO
    END DO
!$OMP END PARALLEL
    CALL vmax(vfric_sfc, all=.TRUE.)
    CALL vmax(vfric_btm, all=.TRUE.)
    CALL update_boundary(vfric_sfc)
    CALL update_boundary(vfric_btm)


!$OMP PARALLEL DO PRIVATE(am, an, sm, sh, km, kh, tke, gls, tmp, c3, sqrttke)
    DO k=0, ksize+1
    DO j=0, jsize+1
    DO i=0, isize+1
       IF (.NOT. lmask3d(i,j,k)) THEN
          visc(i,j,k) = 0.0
          diff(i,j,k) = 0.0
          tlen(i,j,k) = 0.0
          CYCLE
       END IF

       tke = max(tracer(i,j,k,tracer_index_tke), min_tke)
       gls = max(tracer(i,j,k,tracer_index_gls), min_gls)

!!!!!!YJK
!       IF (tracer(i,j,k,tracer_index_tke) < min_tke) THEN
!         tracer(i,j,k,tracer_index_tke) = min_tke
!       END IF
!       IF (tracer(i,j,k,tracer_index_gls) < min_gls) THEN
!         tracer(i,j,k,tracer_index_gls) = min_gls 
!       END IF 
!!!!!!YJK 

       sqrttke = sqrt(tke)

       ! switch c3 parameter by stratification stability
       IF (bfreq2(i,j,k-1)+bfreq2(i,j,k) < 0.0) THEN
          c3 = gls_param%c3p
       ELSE
          c3 = gls_param%c3m
!!!YJK
          IF (n > 0) THEN
             gls = min(gls, (0.56)**(n/2) * c0**p * tke**(m+n/2) * (0.5*(bfreq2(i,j,k-1)+bfreq2(i,j,k)))**(-n/2))
          ELSE
              gls = max(gls, (0.56)**(n/2) * c0**p * tke**(m+n/2) * (0.5*(bfreq2(i,j,k-1)+bfreq2(i,j,k)))**(-n/2))
!             gls = max(gls, ((0.56)*(c0**(p/n)) * (tke**(m/n+1/2)) * (0.5*(bfreq2(i,j,k-1)+bfreq2(i,j,k)))**(-1/2))**(n))
          END IF
       END IF

       tlen(i,j,k) = c0**(-p/n) * tke**(-m/n) * gls**(1.0/n)
!!!!!!!!YJK
!       IF (bfreq2(i,j,k-1)+bfreq2(i,j,k) > 0.0) THEN
!         tlen(i,j,k)=min(c0**(-p/n) * tke**(-m/n) * gls**(1.0/n), sqrt(0.56*tke*0.5*(bfreq2(i,j,k-1)+bfreq2(i,j,k))))
!       ELSE
!         tlen(i,j,k) = c0**(-p/n) * tke**(-m/n) * gls**(1.0/n)
!       End IF
      
       tmp = tlen(i,j,k)**2 / (c0**6*tke)  !(tke**2/epslion**2)
       am =  tmp * 0.5*(sfreq2(i,j,k-1)+sfreq2(i,j,k))
       an =  tmp * 0.5*(bfreq2(i,j,k-1)+bfreq2(i,j,k))

       tmp = 1.0 + 0.2555*an + 0.02872*am + 0.008677*an*an + 0.005222*an*am - 0.0000337*am*am
       sm  = (0.1070 + 0.01741*an - 0.00012*am) / (tmp*c0**3)
       sh  = (0.1120 + 0.00452*an + 0.00088*am) / (tmp*c0**3)

       km = max(min(sm*sqrt2*sqrttke*tlen(i,j,k), max_visc), 0.0)
       kh = max(min(sh*sqrt2*sqrttke*tlen(i,j,k), max_diff), 0.0)

       ! implict scheme for disipation
       tmp = c0**3 * sqrttke * dtime

       tracer(i,j,k,tracer_index_tke) = tke * tlen(i,j,k) / (tlen(i,j,k)+tmp) &
            + dtime * 0.5*( (km+viscv)*(sfreq2(i,j,k-1)+sfreq2(i,j,k)) &
                           -(kh+diffv)*(bfreq2(i,j,k-1)+bfreq2(i,j,k)))

       tracer(i,j,k,tracer_index_gls) = gls * tlen(i,j,k) / (tlen(i,j,k)+tmp*c2) &
            + dtime * 0.5*( (km+viscv)*(sfreq2(i,j,k-1)+sfreq2(i,j,k))*c1 &
                           -(kh+diffv)*(bfreq2(i,j,k-1)+bfreq2(i,j,k))*c3)*gls/tke

       IF (k==surface_k(i,j)) THEN  ! surface fulx
          tracer(i,j,k,tracer_index_tke) = tracer(i,j,k,tracer_index_tke) &
               + dtime * 0.0 * vfric_sfc(i,j)**3*dsz(i,j)/dvol(i,j,k) ! 100 is the empirical param following UB03
          tracer(i,j,k,tracer_index_gls) = tracer(i,j,k,tracer_index_gls) &
               - dtime * n*km/gls_param%sgls* c0**p * tke**m * (karman_const*dz(k)*0.5)**n / (dz(k)*0.5) * dsz(i,j)/dvol(i,j,k)
       END IF

       IF (k==bottom_k(i,j)) THEN  ! bottom flux
          tracer(i,j,k,tracer_index_tke) = tracer(i,j,k,tracer_index_tke) &
               + dtime * 0.0 * vfric_btm(i,j)**3*dsz(i,j)/dvol(i,j,k) ! 100 is the empirical param following UB03
          tracer(i,j,k,tracer_index_gls) = tracer(i,j,k,tracer_index_gls) &
               - dtime * n*km/gls_param%sgls* c0**p * tke**m * (karman_const*dz(k)*0.5)**n / (dz(k)*0.5) * dsz(i,j)/dvol(i,j,k)
       END IF
!!!!!!YJK nothing done to code just to check
!      IF (tracer(i,j,k,tracer_index_tke) < min_tke) THEN
!        tracer(i,j,k,tracer_index_tke) = min_tke
!      END IF
!      IF (tracer(i,j,k,tracer_index_gls) < min_gls) THEN
!        tracer(i,j,k,tracer_index_gls) = min_gls 
!      END IF 
!!!!!!YJK       
       visc(i,j,k) = km
       diff(i,j,k) = kh
    END DO
    END DO
    END DO

    CALL checkout('VISCT', visc) ! vertical viscosity   by GLS turbulent closure [m^2/s]
    CALL checkout('DIFFT', diff) ! vertical diffusivity by GLS turbulent closure [m^2/s]


    CALL update_tracer_boundary(tracer_index_tke)
    CALL update_tracer_boundary(tracer_index_gls)


    CALL checkout('TLEN', tlen)   ! turbulent length scale [m]

    IF (require_checkout('TEPS')) THEN
!$OMP PARALLEL DO
       DO k=1, ksize
       DO j=1, jsize
       DO i=1, isize
          IF (tlen(i,j,k) > 1.0E-30) THEN
             tlen(i,j,k) = c0**3 * tracer(i,j,k,tracer_index_tke)**1.5 / tlen(i,j,k)
          ELSE
             tlen(i,j,k) = 0.0
          END IF
       END DO
       END DO
       END DO

       CALL checkout('TEPS', tlen) ! turbulent dissipasion rate per unit mass [m^2/s^3]
    END IF

    IF (report) CALL report_gls

  END SUBROUTINE gls_mixing


  SUBROUTINE report_gls
    INTEGER :: i, j, k

    INTEGER(8) :: hist(0:6)
    REAL(8)    :: total

    hist(:) = 0

!$OMP PARALLEL DO REDUCTION(+:hist)
    DO k=1, ksize
    DO j=1, jsize
    DO i=1, isize
       IF (.NOT. lmask3d(i,j,k)) CYCLE

       IF (tracer(i,j,k,tracer_index_tke) >= 1.0) THEN
          hist(0) = hist(0) + 1
       ELSE IF (tracer(i,j,k,tracer_index_tke) >= 1.0E-1) THEN
          hist(1) = hist(1) + 1
       ELSE IF (tracer(i,j,k,tracer_index_tke) >= 1.0E-2) THEN
          hist(2) = hist(2) + 1
       ELSE IF (tracer(i,j,k,tracer_index_tke) >= 1.0E-3) THEN
          hist(3) = hist(3) + 1
       ELSE IF (tracer(i,j,k,tracer_index_tke) >= 1.0E-4) THEN
          hist(4) = hist(4) + 1
       ELSE IF (tracer(i,j,k,tracer_index_tke) >= 1.0E-5) THEN
          hist(5) = hist(5) + 1
       ELSE
          hist(6) = hist(6) + 1
       END IF
    END DO
    END DO
    END DO

    CALL gsum(hist)

    IF (rank==0) THEN
       total = REAL(sum(hist),8)

       WRITE(REPORT_UNIT, '(A)', ADVANCE='NO') " TKE histgram:"
       DO i=6, 1, -1
          WRITE(REPORT_UNIT, '("(>1.0E-",I1, ")", X, F6.3, ", ")', ADVANCE='NO') i, hist(i) / total
       END DO
       WRITE(REPORT_UNIT, '("(>1.0E0)",X, F6.3)') hist(0)/total
    END IF

  END SUBROUTINE report_gls

END MODULE gls
