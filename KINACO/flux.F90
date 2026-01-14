#include "macro.h"

MODULE flux
  USE misc
  USE geometry
  USE velocity
  USE tracers
  USE seaice
  USE parameters
  USE io
  IMPLICIT NONE

  REAL(4), ALLOCATABLE :: aieff(:,:)
  REAL(4), ALLOCATABLE :: hieff(:,:)
  REAL(4), ALLOCATABLE :: dtair(:,:)
  REAL(8), ALLOCATABLE :: wind_x(:,:)
  REAL(8), ALLOCATABLE :: wind_y(:,:)
  REAL(4), ALLOCATABLE :: abswind(:,:)
  REAL(8), ALLOCATABLE :: prec(:,:)
  REAL(8), ALLOCATABLE :: prec_tracer(:,:,:)
  REAL(8), ALLOCATABLE :: evap(:,:)
  REAL(8), ALLOCATABLE :: evap_tracer(:,:,:)
  REAL(4), ALLOCATABLE :: c_drag(:,:)
  REAL(4), ALLOCATABLE :: c_late(:,:)
  REAL(4), ALLOCATABLE :: c_sens(:,:)

  LOGICAL, SAVE :: stat_prec
  LOGICAL, SAVE :: stat_evap
  LOGICAL, SAVE :: stat_prec_tracer(max_tracer)
  LOGICAL, SAVE :: stat_evap_tracer(max_tracer)

CONTAINS
  SUBROUTINE update_flux
    REAL(4) :: absv, qsat
    LOGICAL :: stat
    INTEGER :: i, j, k, n

    IF (.NOT. ALLOCATED(wind_x)) ALLOCATE(wind_x(1-slv:isize+slv, 1-slv:jsize+slv))
    IF (.NOT. ALLOCATED(wind_y)) ALLOCATE(wind_y(1-slv:isize+slv, 1-slv:jsize+slv))

    CALL checkin('WIND_X', wind_x, stat)
    IF (.NOT. stat)  wind_x(:,:) = 0.0

    CALL checkin('WIND_Y', wind_y, stat)
    IF (.NOT. stat)  wind_y(:,:) = 0.0

    IF (.NOT. ALLOCATED(abswind)) ALLOCATE(abswind(1-slv:isize+slv, 1-slv:jsize+slv))

    IF (wind_relative) THEN
       IF (vrank==0) THEN
          DO j=1-slv, jsize+slv
          DO i=1-slv, isize+slv
             abswind(i,j) = sqrt( (wind_x(i,j) - 0.5*(u(i-1,j,ksize)+u(i,j,ksize)))**2 &
                                 +(wind_y(i,j) - 0.5*(v(i,j-1,ksize)+v(i,j,ksize)))**2)
          END DO
          END DO
       END IF
!      CALL vcast(abswind)
!      if abswind is used in the lower-grid, uncomment above vcast call
    ELSE
       DO j=1-slv, jsize+slv
       DO i=1-slv, isize+slv
          abswind(i,j) = sqrt(wind_x(i,j)**2 + wind_y(i,j)**2)
       END DO
       END DO
    END IF

    IF (.NOT. ALLOCATED(dtair))  ALLOCATE(dtair(1-slv:isize+slv, 1-slv:jsize+slv))
    dtair(:,:) = 0.0
    CALL checkin('TAIR', dtair, stat)
    IF (stat) THEN
       CALL assert(tracer_index_t/=0, "'TAIR' input requires tracer 'T'")
       IF (vrank==0) THEN
          DO j=1-slv, jsize+slv
          DO i=1-slv, isize+slv
             dtair(i,j) = imask3d(i,j,ksize)*(dtair(i,j) - REAL(tracer(i,j,ksize,tracer_index_t),4))
          END DO
          END DO
       END IF
!      CALL vcast(dtair)
    END IF

    IF (.NOT. ALLOCATED(c_drag)) ALLOCATE(c_drag(1-slv:isize+slv, 1-slv:jsize+slv))
    IF (.NOT. ALLOCATED(c_late)) ALLOCATE(c_late(1-slv:isize+slv, 1-slv:jsize+slv))
    IF (.NOT. ALLOCATED(c_sens)) ALLOCATE(c_sens(1-slv:isize+slv, 1-slv:jsize+slv))

    SELECT CASE(winddrag_scheme)
    CASE (0) !constant C_d
!$OMP PARALLEL WORKSHARE
       c_drag(:,:) = drag_wind
!$OMP END PARALLEL WORKSHARE

    CASE (1) !following Large and Pond (1981), J. Phys. Oceanogr.
!$OMP PARALLEL DO
       DO j=1-slv, jsize+slv
       DO i=1-slv, isize+slv
          c_drag(i,j) = max(drag_wind, (0.49+0.065*abswind(i,j))*1.0E-3)
       END DO
       END DO


    CASE (2) !following Kara et al. (2000), J. Oce. Atm. Tech.
!$OMP PARALLEL DO PRIVATE(absv)
       DO j=1-slv, jsize+slv
       DO i=1-slv, isize+slv
          absv = max(min(REAL(abswind(i,j),4), 32.5), 2.5)
          c_drag(i,j) = (8.620E-4 + (8.80E-5 - 8.900E-7*absv)*absv) &
                       -(1.034E-4 - (6.78E-6 - 1.147E-7*absv)*absv)*dtair(i,j)
          c_drag(i,j) = max(drag_wind, c_drag(i,j))
       END DO
       END DO

    END SELECT

    IF (use_landwater) THEN
!$OMP PARALLEL DO
       DO j=1-slv, jsize+slv
       DO i=1-slv, isize+slv
          IF (lwflag(i,j)) c_drag(i,j) = 0.0
       END DO
       END DO
    END IF

!$OMP PARALLEL DO PRIVATE(absv)
    DO j=1-slv, jsize+slv
    DO i=1-slv, isize+slv
       absv = max(min(REAL(abswind(i,j),4), 27.5), 3.0)
       c_late(i,j) =  ( 9.94E-4 + (6.10E-5 - 1.00E-6*absv)*absv) &
                     -(-2.00E-5 + (6.91E-4 - 8.17E-4/absv)/absv)*dtair(i,j)

       c_sens(i,j) = 0.96*c_late(i,j)
    END DO
    END DO

    IF (.NOT. ALLOCATED(aieff)) THEN
       ALLOCATE(aieff(1-slv:isize+slv, 1-slv:jsize+slv))
       aieff(:,:) = 0.0
    END IF

    IF (.NOT. ALLOCATED(hieff)) THEN
       ALLOCATE(hieff(1-slv:isize+slv, 1-slv:jsize+slv))
       hieff(:,:) = 0.0
    END IF

    IF (seaice_coupled) THEN
       DO j=1-slv, jsize+slv
       DO i=1-slv, isize+slv
          IF (hi(i,j) < hi_min .OR. ai(i,j) < ai_min) THEN
             aieff(i,j) = 0.0
             hieff(i,j) = hi_min
          ELSE
             aieff(i,j) = ai(i,j)
             hieff(i,j) = max(hi(i,j), hi_min)
             !TODO, if multi-category sea ice model is used, hieff hould be the harmonic mean of ice thickness for each category
          END IF
       END DO
       END DO
    END IF

    IF (rigid_lid) RETURN

    IF (.NOT. ALLOCATED(prec)) ALLOCATE(prec(1-slv:isize+slv, 1-slv:jsize+slv))
    IF (.NOT. ALLOCATED(evap)) ALLOCATE(evap(1-slv:isize+slv, 1-slv:jsize+slv))
    IF (n_tracer > 0 .AND. .NOT. ALLOCATED(prec_tracer)) ALLOCATE(prec_tracer(1-slv:isize+slv, 1-slv:jsize+slv, n_tracer))
    IF (n_tracer > 0 .AND. .NOT. ALLOCATED(evap_tracer)) ALLOCATE(evap_tracer(1-slv:isize+slv, 1-slv:jsize+slv, n_tracer))

    CALL checkin('PREC', prec, stat_prec)  ![kg /(m^2 s)] = [mm/s]
    IF (.NOT. stat_prec) CALL checkin('WATER_SUPPLY', prec, stat_prec)

    IF (stat_prec) THEN
       CALL scaleoffset(prec, scale=1.0E-3, offset=0.0) !convert to [m/s]

       DO n=1, n_tracer
          CALL checkin('PREC_' // trim(tracer_name(n)), prec_tracer(:,:,n), stat_prec_tracer(n))
          IF (.NOT. stat_prec_tracer(n)) CALL checkin('WATER_SUPPLY_' // trim(tracer_name(n)), prec_tracer(:,:,n), stat_prec_tracer(n))
          IF (n==tracer_index_s .AND. .NOT. stat_prec_tracer(n)) THEN
             stat_prec_tracer(n) = .TRUE.
             prec_tracer(:,:,n) = 0.0
          END IF
       END DO
    END IF

    CALL checkin('QAIR', evap, stat_evap) ! humidity [kg/kg]
    IF (stat_evap) THEN
       CALL assert(tracer_index_t /= 0, "'QAIR' input requires tracer 'T'")
       IF (vrank==0) THEN
!$OMP PARALLEL DO PRIVATE(qsat)
          DO j=1-slv, jsize+slv
          DO i=1-slv, isize+slv
             IF (imask3d(i,j,ksize)==0) THEN
                evap(i,j) = 0.0
                CYCLE
             END IF
             qsat = vpsat(REAL(tracer(i,j,ksize,tracer_index_t),4)) !saturated vapor pressure at SST in [hPa]
             qsat = 0.622*qsat/(1013.0 - 0.378*qsat)                !saturated humidity at SST [kg/kg] (SLP is fixed to 1013 hPa)

             evap(i,j) = imask3d(i,j,ksize)* rho_air * c_late(i,j)*abswind(i,j)*(0.98*qsat - evap(i,j)) / rho_0 ![m/s]
          END DO
          END DO
       END IF
    END IF

    CALL checkin('EVAP', evap, stat_evap, add=stat_evap.EQV..TRUE.)

    IF (stat_evap) THEN
       DO n=1, n_tracer
          CALL checkin('EVAP_' // trim(tracer_name(n)), evap_tracer(:,:,n), stat_evap_tracer(n))
          IF (n==tracer_index_s .AND. .NOT. stat_evap_tracer(n)) THEN
             stat_evap_tracer(n) = .TRUE.
             evap_tracer(:,:,n) = 0.0
          END IF
       END DO
    END IF

  END SUBROUTINE update_flux

!----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE surface_stress(gx, gy, gix, giy)
    REAL(8), INTENT(INOUT) :: gx(0:isize, 1:jsize, 1:ksize)
    REAL(8), INTENT(INOUT) :: gy(1:isize, 0:jsize, 1:ksize)

    REAL(8), INTENT(INOUT) :: gix(0:isize, 1:jsize)
    REAL(8), INTENT(INOUT) :: giy(1:isize, 0:jsize)

    REAL(8) :: tau_x(1-slv:isize+slv, 1-slv:jsize+slv)
    REAL(8) :: tau_y(1-slv:isize+slv, 1-slv:jsize+slv)

    REAL(8) :: tmpx(0:isize, 1:jsize)
    REAL(8) :: tmpy(1:isize, 0:jsize)

    REAL(8) :: d, wx, wy

    LOGICAL :: stat_x, stat_y, stat

    INTEGER :: i, j, k

!$OMP PARALLEL WORKSHARE
    tmpx(:,:) = 0.0
    tmpy(:,:) = 0.0
!$OMP END PARALLEL WORKSHARE

    CALL checkin('TAU_X', tau_x, stat_x)
    CALL checkin('TAU_Y', tau_y, stat_y)

    CALL checkin('CORIOLIS_Z', corz)

    IF (vrank/=0) RETURN

!$OMP PARALLEL PRIVATE(i, j, d, wx, wy)
    IF (stat_x) THEN
!$OMP DO
       DO j=1, jsize
       DO i=0, isize
          tmpx(i,j) = tmpx(i,j) + imask3d(i,j,ksize)*imask3d(i+1,j,ksize) * (  (1.0-aieff(i,  j)) * tau_x(i,  j) * dsz(i,  j) &
                                                                             + (1.0-aieff(i+1,j)) * tau_x(i+1,j) * dsz(i+1,j)) ! [kg m / s^2] = [N]
       END DO
       END DO

       IF (seaice_coupled) THEN
!$OMP DO
          DO j=1, jsize
          DO i=0, isize
             IF (aieff(i,j)   > 0.0) THEN
                gix(i,j) = gix(i,j) + imask3d(i,j,ksize)*imask3d(i+1,j,ksize) * aieff(i,j)   * tau_x(i,j)   * dsz(i,j)   &
                     / (rho_ice * hieff(i,  j) * (dsz(i,j)+dsz(i+1,j)))     ! [m / s^2]
             END IF

             IF (aieff(i+1,j) > 0.0) THEN
                gix(i,j) = gix(i,j) + imask3d(i,j,ksize)*imask3d(i+1,j,ksize) * aieff(i+1,j) * tau_x(i+1,j) * dsz(i+1,j) &
                     / (rho_ice * hieff(i+1,j) * (dsz(i,j)+dsz(i+1,j)))
             END IF
          END DO
          END DO
       END IF
    END IF

    IF (stat_y) THEN
!$OMP DO
       DO j=0, jsize
       DO i=1, isize
          tmpy(i,j) = tmpy(i,j) + imask3d(i,j,ksize)*imask3d(i,j+1,ksize) * (  (1.0-aieff(i,j))   * tau_y(i,j)   * dsz(i,j) &
                                                                             + (1.0-aieff(i,j+1)) * tau_y(i,j+1) * dsz(i,j+1)) ! [kg m / s^2] = [N]
       END DO
       END DO

       IF (seaice_coupled) THEN
!$OMP DO
          DO j=0, jsize
          DO i=1, isize
             IF (aieff(i,j)   > 0.0) THEN
                giy(i,j) = giy(i,j) + imask3d(i,j,ksize)*imask3d(i,j+1,ksize) * aieff(i,j)   *  tau_y(i,j)  * dsz(i,j)   &
                     / (rho_ice * hieff(i,j)   * (dsz(i,j)+dsz(i,j+1)))     ! [m / s^2]
             END IF

             IF (aieff(i,j+1) > 0.0) THEN
                giy(i,j) = giy(i,j) + imask3d(i,j,ksize)*imask3d(i,j+1,ksize) * aieff(i,j+1) * tau_y(i,j+1) * dsz(i,j+1) &
                     / (rho_ice * hieff(i,j+1) * (dsz(i,j)+dsz(i,j+1)))
             END IF

          END DO
          END DO
       END IF
    END IF

!$OMP DO
    DO j=0, jsize+1
    DO i=0, isize+1
       wx = wind_x(i,j)
       wy = wind_y(i,j)

       IF (wind_relative) THEN
          wx = wx - 0.5*(u(i-1,j,ksize)+u(i,j,ksize)) ! [m/s]
          wy = wy - 0.5*(v(i,j-1,ksize)+v(i,j,ksize))
       END IF

       d = c_drag(i,j) * rho_air * abswind(i,j) ! [kg / (m^2 s)]

       tau_x(i,j) = d*(cos_turn_wind * wx - sign(1.0D0, corz(i,j))*sin_turn_wind * wy) ! [kg / (m s^2)] = [Pa]
       tau_y(i,j) = d*(cos_turn_wind * wy + sign(1.0D0, corz(i,j))*sin_turn_wind * wx)
    END DO
    END DO

!$OMP DO
    DO j=1, jsize
    DO i=0, isize
       tmpx(i,j) = tmpx(i,j) + imask3d(i,j,ksize)*imask3d(i+1,j,ksize) * (  (1.0-aieff(i,  j)) * tau_x(i,  j) * dsz(i,j) &  ! [kg m / s^2] = [N]
                                                                          + (1.0-aieff(i+1,j)) * tau_x(i+1,j) * dsz(i+1,j))
    END DO
    END DO

!$OMP DO
    DO j=0, jsize
    DO i=1, isize
       tmpy(i,j) = tmpy(i,j) + imask3d(i,j,ksize)*imask3d(i,j+1,ksize) * (  (1.0-aieff(i,j))   * tau_y(i,j)   * dsz(i,j) &
                                                                          + (1.0-aieff(i,j+1)) * tau_y(i,j+1) * dsz(i,j+1))
    END DO
    END DO

    IF (seaice_coupled) THEN
!$OMP DO
       DO j=1, jsize
       DO i=0, isize
          IF (aieff(i,j)   > 0.0) THEN
             gix(i,j) = gix(i,j) + imask3d(i,j,ksize)*imask3d(i+1,j,ksize) * aieff(i,j)   * tau_x(i,  j) * dsz(i,j)   &
                  / (rho_ice * hieff(i,  j) * (dsz(i,j)+dsz(i+1,j))) ! [m / s^2]
          END IF

          IF (aieff(i+1,j) > 0.0) THEN
             gix(i,j) = gix(i,j) + imask3d(i,j,ksize)*imask3d(i+1,j,ksize) * aieff(i+1,j) * tau_x(i+1,j) * dsz(i+1,j) &
                  / (rho_ice * hieff(i+1,j) * (dsz(i,j)+dsz(i+1,j)))
          END IF
       END DO
       END DO

!$OMP DO
       DO j=0, jsize
       DO i=1, isize
          IF (aieff(i,j)   > 0.0) THEN
             giy(i,j) = giy(i,j) + imask3d(i,j,ksize)*imask3d(i,j+1,ksize) * aieff(i,j)   * tau_y(i,j)   * dsz(i,j)    &
                  / (rho_ice * hieff(i,j)   * (dsz(i,j)+dsz(i+1,j)))
          END IF

          IF (aieff(i,j+1) > 0.0) THEN
             giy(i,j) = giy(i,j) + imask3d(i,j,ksize)*imask3d(i,j+1,ksize) * aieff(i,j+1) * tau_y(i,j+1) * dsz(i,j+1)  &
                  / (rho_ice * hieff(i,j+1) * (dsz(i,j)+dsz(i,j+1)))
          END IF
       END DO
       END DO

!$OMP DO
       DO j=0, jsize+1
       DO i=0, isize+1
          wx = 0.5*(ui(i-1,j)+ui(i,j)-u(i-1,j,ksize)+u(i,j,ksize)) ![m / s]
          wy = 0.5*(vi(i,j-1)+vi(i,j)-v(i,j-1,ksize)+v(i,j,ksize))

          d = drag_ice * rho_ice * sqrt(wx**2 + wy**2) ! [kg / (m^2 s)]

          tau_x(i,j) = d*(cos_turn_ice * wx - sign(1.0D0, corz(i,j))*sin_turn_ice * wy) ! [kg / (m s^2)] = [Pa]
          tau_y(i,j) = d*(cos_turn_ice * wy + sign(1.0D0, corz(i,j))*sin_turn_ice * wx)
       END DO
       END DO

!$OMP DO
       DO j=1, jsize
       DO i=0, isize
          tmpx(i,j) = tmpx(i,j) + imask3d(i,j,ksize)*imask3d(i+1,j,ksize) * (  aieff(i,  j) * tau_x(i,  j) * dsz(i,  j) &
                                                                             + aieff(i+1,j) * tau_x(i+1,j) * dsz(i+1,j)) ! [kg m / s^2] = [N]
       END DO
       END DO

!$OMP DO
       DO j=0, jsize
       DO i=1, isize
          tmpy(i,j) = tmpy(i,j) + imask3d(i,j,ksize)*imask3d(i,j+1,ksize) * (  aieff(i,j)   * tau_y(i,j)   * dsz(i,j)   &
                                                                             + aieff(i,j+1) * tau_y(i,j+1) * dsz(i,j+1))
       END DO
       END DO

!$OMP DO
       DO j=1, jsize
       DO i=0, isize
          IF (aieff(i,j)   > 0.0) THEN
             gix(i,j)  = gix(i,j)  - imask3d(i,j,ksize)*imask3d(i+1,j,ksize) * aieff(i,j)   * tau_x(i,j)   * dsz(i,j)   &
                  / (rho_ice * hieff(i,j) * (dsz(i,j)+dsz(i+1,j))) ! [m / s^2]
          END IF

          IF (aieff(i+1,j) > 0.0) THEN
             gix(i,j)  = gix(i,j)  - imask3d(i,j,ksize)*imask3d(i+1,j,ksize) * aieff(i+1,j) * tau_x(i+1,j) * dsz(i+1,j) &
                  / (rho_ice * hieff(i+1,j) * (dsz(i,j)+dsz(i+1,j)))
          END IF
       END DO
       END DO

!$OMP DO
       DO j=0, jsize
       DO i=1, isize
          IF (aieff(i,j)   > 0.0) THEN
             giy(i,j)  = giy(i,j)  - imask3d(i,j,ksize)*imask3d(i,j+1,ksize) * aieff(i,j)   * tau_y(i,j)   * dsz(i,j)   &
                  / (rho_ice * hieff(i,j) * (dsz(i,j)+dsz(i,j+1)))
          END IF

          IF (aieff(i,j+1) > 0.0) THEN
             giy(i,j)  = giy(i,j)  - imask3d(i,j,ksize)*imask3d(i,j+1,ksize) * aieff(i,j+1) * tau_y(i,j+1) * dsz(i,j+1) &
                  / (rho_ice * hieff(i,j+1) * (dsz(i,j)+dsz(i,j+1)))
          END IF
       END DO
       END DO
    END IF

!$OMP DO
    DO j=1, jsize
    DO i=0, isize
       gx(i,j,ksize) = gx(i,j,ksize) + imask3d(i-1,j,ksize)*imask3d(i,j,ksize)*tmpx(i,j)/(rho_0*(dvol(i,j,ksize)+dvol(i+1,j,ksize))) ! [m / s^2]
       taux_sfc(i,j) = taux_sfc(i,j) + imask3d(i-1,j,ksize)*imask3d(i,j,ksize)*tmpx(i,j)/(dsz(i,j)+dsz(i+1,j))                       ! [Pa]
    END DO
    END DO

!$OMP DO
    DO j=0, jsize
    DO i=1, isize
       gy(i,j,ksize) = gy(i,j,ksize) + imask3d(i,j-1,ksize)*imask3d(i,j,ksize)*tmpy(i,j)/(rho_0*(dvol(i,j,ksize)+dvol(i,j+1,ksize)))
       tauy_sfc(i,j) = tauy_sfc(i,j) + imask3d(i,j-1,ksize)*imask3d(i,j,ksize)*tmpy(i,j)/(dsz(i,j)+dsz(i,j+1))
    END DO
    END DO
!$OMP END PARALLEL

  END SUBROUTINE surface_stress

!----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE bottom_stress(gx, gy)
    REAL(8), INTENT(INOUT) :: gx(0:isize, 1:jsize, 1:ksize)
    REAL(8), INTENT(INOUT) :: gy(1:isize, 0:jsize, 1:ksize)

    REAL(8) :: tmp

    INTEGER :: i, j, k, kl0, kl1

    IF (.NOT. bulk_bottom_stress) RETURN

!$OMP PARALLEL PRIVATE(tmp, i, j, k)
!$OMP DO
    DO j=1, jsize
    DO i=0, isize
       k = max(bottom_k(i,j), bottom_k(i+1,j))
       IF (k < 1 .OR. k > ksize) CYCLE

       tmp = drag_bottom * u(i,j,k) &
            * ((1-imask3d(i,  j,k-1))*sqrt(u(i,j,k)**2 + 0.25*(v(i,  j-1,k)+v(i,  j,k))**2) * dsz(i,  j) &
             + (1-imask3d(i+1,j,k-1))*sqrt(u(i,j,k)**2 + 0.25*(v(i+1,j-1,k)+v(i+1,j,k))**2) * dsz(i+1,j)) ! [m^4/s^2]

       taux_btm(i,j) = taux_btm(i,j) - rho_0 * tmp / (dsz(i,j)+dsz(i+1,j)) ! [Pa]
       gx(i,j,k) = gx(i,j,k) - tmp / (dvol(i,j,k)+dvol(i+1,j,k))           ! [m/s^2]
    END DO
    END DO

!$OMP DO
    DO j=0, jsize
    DO i=1, isize
       k = max(bottom_k(i,j), bottom_k(i,j+1))
       IF (k < 1 .OR. k > ksize) CYCLE

       tmp = drag_bottom * v(i,j,k) &
            * ((1-imask3d(i,j,  k-1))*sqrt(v(i,j,k)**2 + 0.25*(u(i-1,j,  k)+u(i,j,  k))**2) * dsz(i,j)   &
             + (1-imask3d(i,j+1,k-1))*sqrt(v(i,j,k)**2 + 0.25*(u(i-1,j+1,k)+u(i,j+1,k))**2) * dsz(i,j+1))

       tauy_btm(i,j) = tauy_btm(i,j) - rho_0 * tmp / (dsz(i,j)+dsz(i,j+1))
       gy(i,j,k) = gy(i,j,k) - tmp / (dvol(i,j,k)+dvol(i,j+1,k))
    END DO
    END DO
!$OMP END PARALLEL

  END SUBROUTINE bottom_stress

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE heatflux(t)
    REAL(TRC_KIND), INTENT(INOUT) :: t(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)

    REAL(8) :: tmp(1-slv:isize+slv, 1-slv:jsize+slv)
    REAL(8) :: net(1-slv:isize+slv, 1-slv:jsize+slv)

    net(:,:) = 0.0

    CALL shortwave
    CALL longwave
    CALL sensible
    CALL latent

    CALL checkout('NETHF', net)

  CONTAINS
    SUBROUTINE shortwave
      REAL(8) :: dt, c
      LOGICAL :: stat
      INTEGER :: i, j, k

      CALL checkin('RADSW', tmp, stat) ! downward shortwave radiation [W/m^2]

      IF (.NOT. stat) RETURN

      IF (vrank==0) THEN
         DO j=1-slv, jsize+slv
         DO i=1-slv, isize+slv
            tmp(i,j) = imask3d(i,j,ksize)*(1.0-albedo_ocean)*(1.0-aieff(i,j))*tmp(i,j) ![W/m^2]
            net(i,j) = net(i,j) + tmp(i,j)
         END DO
         END DO
      ELSE
         tmp(:,:) = 0.0
      END IF

      CALL checkout('NETHF_SW', tmp)

      IF (vrank==0) THEN
         DO j=1-slv, jsize+slv
         DO i=1-slv, isize+slv
            tmp(i,j) = tmp(i,j)*dsz(i,j)*dtime/(cp*rho_0)                           ![K m^3]
            t(i,j,ksize) = t(i,j,ksize) + radsw_r*tmp(i,j) / dvol(i,j,ksize)
            tmp(i,j) = (1.0-radsw_r)*tmp(i,j)
         END DO
         END DO
      ELSE
         CALL urecv(tmp)
      END IF

      DO k=ksize, 1, -1
         c = 1.0 - exp(-dz(k)/radsw_z)
!$OMP PARALLEL DO PRIVATE(dt)
         DO j=1-slv, jsize+slv
         DO i=1-slv, isize+slv

            IF (tmp(i,j) <= 0.0) CYCLE

            IF (imask3d(i,j,k)==0) THEN
               dt = 0.0
               tmp(i,j) = 0.0
            ELSE IF (imask3d(i,j,k-1)==0 .OR. depth(k) > 5*radsw_z) THEN
               dt = tmp(i,j)
               tmp(i,j) = 0.0
            ELSE
               dt = tmp(i,j)*c
               tmp(i,j) = tmp(i,j)-dt
            END IF

            t(i,j,k) = t(i,j,k) + imask3d(i,j,k)*dt / dvol(i,j,k)
         END DO
         END DO
      END DO

      CALL usend(tmp)

    END SUBROUTINE shortwave

    SUBROUTINE longwave
      LOGICAL :: stat
      INTEGER :: i, j

      CALL checkin('RADLW', tmp, stat) ! downward longwave radiation [W/m^2]
      IF (.NOT. stat) RETURN

      IF (vrank==0) THEN
         DO j=1-slv, jsize+slv
         DO i=1-slv, isize+slv
            tmp(i,j) = imask3d(i,j,ksize)*(1.0-aieff(i,j))*emissivity_ocean*(tmp(i,j)-sigma_sb*(t(i,j,ksize)+273.15)**4) ![W/m^2]
            net(i,j) = net(i,j) + tmp(i,j)
            t(i,j,ksize) = t(i,j,ksize) + tmp(i,j)*dsz(i,j)*dtime / (cp*rho_0*dvol(i,j,ksize))
         END DO
         END DO
      ELSE
         tmp(:,:) = 0.0
      END IF

      CALL checkout('NETHF_LW', tmp)
    END SUBROUTINE longwave

    SUBROUTINE sensible
      INTEGER :: i, j

      IF (vrank==0) THEN
         DO j=1-slv, jsize+slv
         DO i=1-slv, isize+slv
            tmp(i,j) = imask3d(i,j,ksize)*(1.0-aieff(i,j))*c_sens(i,j)*rho_air*cp_air*abswind(i,j)*dtair(i,j) ![W/m^2]
            net(i,j) = net(i,j) + tmp(i,j)

            t(i,j,ksize) = t(i,j,ksize) + tmp(i,j)*dsz(i,j)*dtime / (cp*rho_0*dvol(i,j,ksize))
         END DO
         END DO
      ELSE
         tmp(:,:) =  0.0
      END IF

      CALL checkout('NETHF_SENS', tmp)

    END SUBROUTINE sensible

    SUBROUTINE latent
      INTEGER :: i, j

      IF (.NOT. stat_evap)    RETURN
      IF (vrank==0) THEN
         DO j=1-slv, jsize+slv
         DO i=1-slv, isize+slv
            tmp(i,j) = -imask3d(i,j,ksize)*evap(i,j) * l_vapor * rho_0     ! [m/s] [J/kg] [kg/m^3] = [W/m^2]
            net(i,j) = net(i,j) + tmp(i,j)
            t(i,j,ksize) = t(i,j,ksize) + tmp(i,j)*dsz(i,j)*dtime / (cp*rho_0*dvol(i,j,ksize))
         END DO
         END DO
      ELSE
         tmp(:,:) = 0.0
      END IF

      CALL checkout('NETHF_LATE', tmp)
    END SUBROUTINE latent


  END SUBROUTINE heatflux

!----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE report_stress
    REAL(8) :: max_sfc, max_btm
    INTEGER :: i, j

    IF (stress_report_interval <= 0) RETURN

    IF (mod(n_timestep, stress_report_interval) == 0) THEN

       max_sfc = 0.0
       max_btm = 0.0

       DO j=1, jsize
       DO i=1, isize
          max_sfc = max(max_sfc, 0.5*sqrt((taux_sfc(i-1,j)+taux_sfc(i,j))**2 + (tauy_sfc(i,j-1)+tauy_sfc(i,j))**2))
          max_btm = max(max_btm, 0.5*sqrt((taux_btm(i-1,j)+taux_btm(i,j))**2 + (tauy_btm(i,j-1)+tauy_btm(i,j))**2))
       END DO
       END DO

       CALL gmax(max_sfc)
       CALL gmax(max_btm)

       IF (rank==0) THEN
          WRITE(REPORT_UNIT, '("maxmum surface and bottom stress = ", ES10.3, ", ", ES10.3, " [Pa]")') max_sfc, max_btm
       END IF
    ENDIF

  END SUBROUTINE report_stress

!----------------------------------------------------------------------------------------------------------------------

  REAL(4) PURE FUNCTION vpsat(t)
    !return saturated vapor pressure in [hPa]
    REAL(4), INTENT(IN) :: t ! temperature in degC

    vpsat = 6.112*exp(17.67*t/(t+243.50)) ! following Bolton 1980 Mon. Wea. Rev.
  END FUNCTION vpsat

END MODULE flux
