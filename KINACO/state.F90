#include "macro.h"

MODULE state
  USE misc
  USE io
  USE tracers
  IMPLICIT NONE

  REAL(8), ALLOCATABLE :: rho(:,:,:)    ! in-situ water density [kg/m^3]
  REAL(8), ALLOCATABLE :: sigma0(:,:,:) ! potental density anomaly (sigma_theta) [kg/m^3]
  REAL(8), ALLOCATABLE :: ri(:,:,:)     ! Richardson Number
  REAL(8), ALLOCATABLE :: bfreq2(:,:,:) ! square of buoyancy frequency [s^-2]
  REAL(8), ALLOCATABLE :: mld(:,:)      ! mixed layer depth [m]

  LOGICAL, SAVE :: calc_cint = .FALSE.
  LOGICAL, SAVE :: calc_mld  = .FALSE.

CONTAINS

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE init_state
    USE geometry
    USE parameters, ONLY: eq_of_state

    ALLOCATE(rho(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv))
    rho(:,:,:) = rho_0

    ALLOCATE(sigma0(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv))
    sigma0(:,:,:) = rho_0 - 1000.0

    ALLOCATE(mld(1-slv:isize+slv, 1-slv:jsize+slv))
    mld(:,:) = 0.0

    IF (eq_of_state == 2) THEN
       IF (tracer_index_t==0) CALL add_tracer('T', tracer_index_t)
       IF (tracer_index_s==0) CALL add_tracer('S', tracer_index_s)
    END IF

    IF (eq_of_state == 1 .AND. tracer_index_t==0 .AND. tracer_index_s==0) eq_of_state = 0

    ALLOCATE(bfreq2(0:isize+1, 0:jsize+1, -1:ksize+1))
    ALLOCATE(ri(    0:isize+1, 0:jsize+1, -1:ksize+1))

    bfreq2(:,:,:) = 0.0
    ri(:,:,:)     = 0.0

    CALL update_density(no_checkout=.TRUE.)

  END SUBROUTINE init_state

!-----------------------------------------------------------------------------------------------------------------------

  REAL(8) PURE FUNCTION insitu_to_potential(t, s, p_db)
    REAL(8), INTENT(IN) :: t ! in-situ temperature [C deg]
    REAL(8), INTENT(IN) :: s ! salinity [psu]
    REAL(8), INTENT(IN) :: p_db ! pressure offset from sea surface [db] (or depth [m])

    REAL(8), PARAMETER :: a100 =  0.36504D-4
    REAL(8), PARAMETER :: a101 =  0.83198D-5
    REAL(8), PARAMETER :: a102 = -0.54065D-7
    REAL(8), PARAMETER :: a103 =  0.40274D-9
    REAL(8), PARAMETER :: a110 =  0.17439D-5
    REAL(8), PARAMETER :: a111 = -0.29778D-7
    REAL(8), PARAMETER :: a210 = -0.41057D-10
    REAL(8), PARAMETER :: a200 =  0.89309D-8
    REAL(8), PARAMETER :: a201 = -0.31628D-9
    REAL(8), PARAMETER :: a202 =  0.21987D-11
    REAL(8), PARAMETER :: a300 = -0.16056D-12
    REAL(8), PARAMETER :: a301 =  0.50484D-14

    insitu_to_potential = t &
         - p_db    * (a100 + t*(a101 + t*(a102  + t*a103)) + (s-35.0)*(a110 + t*a111)) &
         - p_db**2 * (a200 + t*(a201 + t*a202) + (s-35.0D0)*a210) &
         - p_db**3 * (a300 + t*a301)

    ! c.f., Bryden (1973): Deep sea res., 20, 401-408
  END FUNCTION insitu_to_potential

!-----------------------------------------------------------------------------------------------------------------------

  REAL(8) PURE FUNCTION potential_to_insitu(t, s, p_db)
    REAL(8), INTENT(IN) :: t ! ptential temperature (0db).
    REAL(8), INTENT(IN) :: s ! in-situ salinity [psu]
    REAL(8), INTENT(IN) :: p_db ! pressure offset from sea surface [db] (or depth [m])

    potential_to_insitu = potential_temperature(t, s, 0.0D0, p_db)

  END FUNCTION potential_to_insitu

!-----------------------------------------------------------------------------------------------------------------------

  REAL(8) PURE FUNCTION atg(t, s, p_db)
    REAL(8), INTENT(IN) :: t         ! [deg C]
    REAL(8), INTENT(IN) :: s         ! [psu]
    REAL(8), INTENT(IN) :: p_db      ! pressure offset from sea surface [db] (or depth [m])

    !Compute adiabatic temperature gradient [K/dBar], used in ptential_temperature
    !from ptem.pro <http://www.csag.uct.ac.za/~daithi/idl_lib/pro/potem.pro>
    ! based on Bryden (1974)

    atg = (((-2.1687D-16*t + 1.8676D-14)*t - 4.6206D-13)*p_db &
         + ((2.7759D-12*t - 1.1351D-10)*(s-35.0D0) + ((-5.4481D-14*t + 8.733D-12)*t - 6.7795D-10)*t + 1.8741D-8))*p_db &
         + (-4.2393D-8*t + 1.8932D-6)*(s-35.0D0) + ((6.6228D-10*t - 6.836d-8)*t + 8.5258D-6)*t + 3.5803D-5
  END FUNCTION atg

!-----------------------------------------------------------------------------------------------------------------------

  REAL(8) PURE FUNCTION potential_temperature(t, s, p_db, pref_db)
    REAL(8), INTENT(IN) :: t         ! [deg C]
    REAL(8), INTENT(IN) :: s         ! [psu]
    REAL(8), INTENT(IN) :: p_db      ! pressure offset from sea surface [db] (or depth [m])
    REAL(8), INTENT(IN) :: pref_db   ! refference pressure offset [db] (or depth [m])

    REAL(8) :: tt
    REAL(8) :: ss
    REAL(8) :: pp
    REAL(8) :: ppr
    REAL(8) :: h, xk, q

    tt = t
    ss = s
    pp = p_db
    ppr = pref_db

    h = ppr-pp
    xk = h*ATG(tt,ss,pp)
    tt = tt + 0.5D0*xk
    q = xk
    pp = pp + 0.5D0*h
    xk = h*atg(tt,ss,pp)
    tt = tt + 0.29289322D0*(xk-q)
    q = 0.58578644D0*xk + 0.121320344D0*q
    xk = h*ATG(tt,ss,pp)
    tt = tt + 1.707106781*(xk-q)
    q = 3.414213562D0*xk - 4.121320344D0*q
    pp = pp + 0.5D0*h
    xk = h*ATG(tt,ss,pp)
    potential_temperature = tt + (xk-2.0D0*q)/6.0D0

    !from ptem.pro <http://www.csag.uct.ac.za/~daithi/idl_lib/pro/potem.pro>
    ! checkvalues: s (psu)  t (deg C)  p (dB)  pr (dB)  theta (deg C)
    !               34.75      1.0      4500      0         0.640
    !               34.75      1.0      4500     4000       0.944
    !               34.95      2.5      3500      0         2.207
    !               34.95      2.5      3500     4000       2.558
  END FUNCTION potential_temperature

!-----------------------------------------------------------------------------------------------------------------------

  REAL(8) PURE FUNCTION freezing_temperature(s, p_db)
    REAL(8), INTENT(IN) :: s ! salinity [psu]
    REAL(8), INTENT(IN) :: p_db ! pressure offset from sea surface [db] (or depth [m])

    REAL(8), PARAMETER :: a =  1.710523D-3
    REAL(8), PARAMETER :: b = -2.154996D-4
    REAL(8), PARAMETER :: c = -0.0575
    REAL(8), PARAMETER :: d = -7.53D-4

    freezing_temperature = a*(max(s,0.0)**1.5) + b*s*s + c*s + d*p_db

    ! c.f., Gill (1982)
  END FUNCTION freezing_temperature

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE update_density(no_checkout)
    USE geometry
    USE parameters, ONLY: gravity, rho_0, eq_of_state, eos_t_ref, eos_s_ref, depth_offset
    USE velocity,   ONLY: dudz, dvdz, sfreq2

    LOGICAL, INTENT(IN), OPTIONAL :: no_checkout

    REAL(8) :: tmp(isize,jsize,ksize)

    REAL(8) :: gr
    INTEGER :: i, j, k, n, km0, km1

    gr = gravity / rho_0

    SELECT CASE(eq_of_state)
    CASE (1) ! Linear EOS
!$OMP PARALLEL
       IF (tracer_index_t/=0 .AND. tracer_index_s/=0) THEN
!$OMP DO
          DO k=1-slv, ksize+slv
          DO j=1-slv, jsize+slv
          DO i=1-slv, isize+slv
             rho(i,j,k) = eos_linear(real(tracer(i,j,k,tracer_index_t),8), real(tracer(i,j,k,tracer_index_s),8))
          END DO
          END DO
          END DO
       ELSE IF (tracer_index_t/=0) THEN
!$OMP DO
          DO k=1-slv, ksize+slv
          DO j=1-slv, jsize+slv
          DO i=1-slv, isize+slv
             rho(i,j,k) = eos_linear(real(tracer(i,j,k,tracer_index_t),8), real(eos_s_ref,8))
          END DO
          END DO
          END DO
       ELSE IF (tracer_index_s/=0) THEN
          DO k=1-slv, ksize+slv
          DO j=1-slv, jsize+slv
          DO i=1-slv, isize+slv
             rho(i,j,k) = eos_linear(real(eos_t_ref,8), real(tracer(i,j,k,tracer_index_s),8))
          END DO
          END DO
          END DO
       END IF
!$OMP DO
       DO k=1-slv, ksize+slv
       DO j=1-slv, jsize+slv
       DO i=1-slv, isize+slv
          sigma0(i,j,k) = rho(i,j,k) - 1000.0
       END DO
       END DO
       END DO
!$OMP DO
       DO k=-1, ksize+1
          km0 = maskindex(k)
          km1 = maskindex(k+1)

          DO j=0, jsize+1
          DO i=0, isize+1
             bfreq2(i,j,k) = -imask3d(i,j,km0)*imask3d(i,j,km1)*(rho(i,j,k+1)-rho(i,j,k))*idz1(k) * gr
          END DO
          END DO
       END DO
!$OMP END PARALLEL

    CASE (2) ! UNESCO EOS-80 following Mellor (1991).
!$OMP PARALLEL
!$OMP DO
       DO k=1-slv, ksize+slv
       DO j=1-slv, jsize+slv
       DO i=1-slv, isize+slv
          sigma0(i,j,k) = sigma_theta(tracer(i,j,k,tracer_index_t), tracer(i,j,k,tracer_index_s))
          rho(i,j,k) = eos_mellor(sigma0(i,j,k), tracer(i,j,k,tracer_index_t), tracer(i,j,k,tracer_index_s), depth_offset+depth(k))
       END DO
       END DO
       END DO

!$OMP DO
       DO k=-1, ksize+1
          km0 = maskindex(k)
          km1 = maskindex(k+1)

          DO j=0, jsize+1
          DO i=0, isize+1
             bfreq2(i,j,k) = eos_mellor(sigma0(i,j,k+1), tracer(i,j,k+1,tracer_index_t), tracer(i,j,k+1,tracer_index_s), depth_offset+depth(k))
          END DO
          END DO

          DO j=0, jsize+1
          DO i=0, isize+1
             bfreq2(i,j,k) = -imask3d(i,j,km0)*imask3d(i,j,km1)*(bfreq2(i,j,k)-rho(i,j,k))*idz1(k) * gr

             ri(i,j,k) = bfreq2(i,j,k) / max(1.0D-12, sfreq2(i,j,k))
          END DO
          END DO
       END DO
!$OMP END PARALLEL
    END SELECT

    IF (require_checkout('CINT') .OR. calc_cint) THEN
       cint(:,:) = 0.0
       DO k=1, ksize
          DO j=1, jsize
          DO i=1, isize
             cint(i,j) = cint(i,j) + bfreq2(i,j,k)*0.5*(dz_ref(i,j,k)+dz_ref(i,j,k+1))
          END DO
          END DO
       END DO
       CALL vsum(cint, ALL=.TRUE.)
       DO j=1, jsize
       DO i=1, isize
          cint(i,j) = sqrt(max(cint(i,j) * wct_ref(i,j), 0.0))
       END DO
       END DO
       CALL update_boundary(cint)
    END IF

    IF (require_checkout('MLD') .OR. calc_mld) CALL mixed_layer_depth

    IF (present(no_checkout)) THEN
       IF (no_checkout) RETURN
    END IF

    CALL checkout('RHO',    rho)      ! in-situ density [kg/m^3]
    CALL checkout('SIGMA0', sigma0)   ! potential density anomaly (sigma_theta) [kg/m^3]
    CALL checkout('BFREQ2', bfreq2)   ! square of buoyancy frequency [1/s^2]
    CALL checkout('RICHARDSON', ri)   ! Richardson number
    CALL checkout('CINT',   cint)     ! phase speed of the 1st mode internal gravity wave [m/s]
    CALL checkout('MLD',    mld)      ! mixed layer depth [m]

    IF (require_checkout('SIGMA1')) THEN
!$OMP PARALLEL DO
       DO k=1, ksize
       DO j=1, jsize
       DO i=1, isize
          tmp(i,j,k) = eos_mellor(sigma0(i,j,k), tracer(i,j,k,tracer_index_t), tracer(i,j,k,tracer_index_s), 1000.0D0) - 1000.0
       END DO
       END DO
       END DO
       CALL checkout('SIGMA1', tmp) ! potential density with reference pressure 1000 m [kg/m^3]
    END IF

    IF (require_checkout('SIGMA2')) THEN
!$OMP PARALLEL DO
       DO k=1, ksize
       DO j=1, jsize
       DO i=1, isize
          tmp(i,j,k) = eos_mellor(sigma0(i,j,k), tracer(i,j,k,tracer_index_t), tracer(i,j,k,tracer_index_s), 2000.0D0) - 1000.0
       END DO
       END DO
       END DO
       CALL checkout('SIGMA2', tmp) ! potential density with reference pressure 2000 m [kg/m^3]
    END IF

    IF (require_checkout('SIGMA3')) THEN
!$OMP PARALLEL DO
       DO k=1, ksize
       DO j=1, jsize
       DO i=1, isize
          tmp(i,j,k) = eos_mellor(sigma0(i,j,k), tracer(i,j,k,tracer_index_t), tracer(i,j,k,tracer_index_s), 3000.0D0) - 1000.0
       END DO
       END DO
       END DO
       CALL checkout('SIGMA3', tmp) ! potential density with reference pressure 3000 m [kg/m^3]
    END IF

    IF (require_checkout('SIGMA4')) THEN
!$OMP PARALLEL DO
       DO k=1, ksize
       DO j=1, jsize
       DO i=1, isize
          tmp(i,j,k) = eos_mellor(sigma0(i,j,k), tracer(i,j,k,tracer_index_t), tracer(i,j,k,tracer_index_s), 4000.0D0) - 1000.0
       END DO
       END DO
       END DO
       CALL checkout('SIGMA4', tmp) ! potential density with reference pressure 4000 m [kg/m^3]
    END IF

    IF (require_checkout('SPICINESS')) THEN
       CALL assert(tracer_index_t/=0 .AND. tracer_index_s/=0, "Both of tracers 'T' and 'S' should be explicitly solved for 'SPICINESS' output")
       CALL assert(eq_of_state /=0,                           "'SPICINESS' is not defined for the constant density (EQ_OF_STATE=0)")

       IF (eq_of_state==1) THEN
!$OMP PARALLEL DO
          DO k=1, ksize
          DO j=1, jsize
          DO i=1, isize
             tmp(i,j,k) = spiciness_linear(tracer(i,j,k,tracer_index_t), tracer(i,j,k,tracer_index_s))
          END DO
          END DO
          END DO
       ELSE
!$OMP PARALLEL DO
          DO k=1, ksize
          DO j=1, jsize
          DO i=1, isize
             tmp(i,j,k) = spiciness_flament(tracer(i,j,k,tracer_index_t), tracer(i,j,k,tracer_index_s))
          END DO
          END DO
          END DO
       END IF

       CALL checkout('SPICINESS', tmp) ! seawater spiciness, see Flament 2002 Prog. Oceanogr.
    END IF

  END SUBROUTINE update_density

!-----------------------------------------------------------------------------------------------------------------------

!should be inline-expanded!
  REAL(8) PURE FUNCTION eos(t, s, p_db)
    USE parameters, ONLY: rho_0, eq_of_state

    REAL(8), INTENT(IN) :: t    ! potential temperature [ C degree]
    REAL(8), INTENT(IN) :: s    ! salinity [ psu ]
    REAL(8), INTENT(IN) :: p_db ! pressure offset from sea surface [db] (or depth [m])

    SELECT CASE (eq_of_state)
    CASE (0) ! Constant
       eos = rho_0
    CASE (1) ! Linear EOS
       eos = eos_linear(real(t,8), real(s,8))
    CASE (2) ! UNESCO EOS-80 following Mellor (1991).
       eos = eos_mellor(sigma_theta(t,s), t, s, p_db)
    END SELECT

  END FUNCTION eos

!-----------------------------------------------------------------------------------------------------------------------

  REAL(8) PURE FUNCTION eos_linear(t, s)
    USE parameters, ONLY: rho_0, alpha=>eos_lin_alpha, beta=>eos_lin_beta, t0=>eos_t_ref, s0=>eos_s_ref

    REAL(8), INTENT(IN) :: t
    REAL(8), INTENT(IN) :: s

    eos_linear = rho_0 * (1.0D0 - alpha * (t - t0) + beta * (s - s0))

  END FUNCTION eos_linear

!-----------------------------------------------------------------------------------------------------------------------

  REAL(8) PURE FUNCTION sigma_theta(t, s)

    REAL(8), INTENT(IN) :: t
    REAL(8), INTENT(IN) :: s

    REAL(8), PARAMETER :: a0 = -1.57406D-1
    REAL(8), PARAMETER :: a1 = +6.793952D-2
    REAL(8), PARAMETER :: a2 = -9.095290D-3
    REAL(8), PARAMETER :: a3 = +1.001685D-4
    REAL(8), PARAMETER :: a4 = -1.120083D-6
    REAL(8), PARAMETER :: a5 = +6.5336332D-9

    REAL(8), PARAMETER :: b0 = +8.24493D-1
    REAL(8), PARAMETER :: b1 = -4.0899D-3
    REAL(8), PARAMETER :: b2 = +7.6438D-5
    REAL(8), PARAMETER :: b3 = -8.2467D-7
    REAL(8), PARAMETER :: b4 = +5.3875D-9

    REAL(8), PARAMETER :: c0 = -5.72466D-3
    REAL(8), PARAMETER :: c1 = +1.0227D-4
    REAL(8), PARAMETER :: c2 = -1.6546D-6

    REAL(8), PARAMETER :: d0 = +4.8314D-4

    sigma_theta = (a0 + (a1 + (a2 + (a3 + (a4 + a5*t)*t)*t)*t)*t )  &
                + (b0 + (b1 + (b2 + (b3 + b4*t)*t)*t)*t) * s        &
                + (c0 + (c1 + c2*t)*t) * (max(s,0.0)**1.5)                 &
                +  d0 * s*s

  END FUNCTION sigma_theta

!-----------------------------------------------------------------------------------------------------------------------

  REAL(8) PURE FUNCTION eos_mellor(sigma0, t, s, p_db)
    REAL(8), INTENT(IN) :: sigma0 ! potential density anomaly  [kg / m^3]
    REAL(8), INTENT(IN) :: t ! potential temperature [ C degree]
    REAL(8), INTENT(IN) :: s ! salinity [ psu ]
    REAL(8), INTENT(IN) :: p_db ! pressure offset from sea surface [db] (or depth [m])

    REAL(8) :: c    ! empirical quantity related to the speed of sound
    REAL(8) :: alpha

    c = 1449.2D0 + 0.00821D0*p_db + (4.55D0 - 0.045D0*t)*t + 1.34D0*(s - 35.0D0) + 15.0D-9 * p_db**2
    alpha = p_db/(c*c)

    eos_mellor = 1000.0D0 + sigma0 + 1.0D4 * ( 1.0D0 - 0.20D0*alpha) * alpha

  END FUNCTION eos_mellor

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE buoyancy(gx, gy, gz)
    USE geometry
    USE parameters, ONLY: buoyancy_scheme, rho_0, gravity, cos_gravity, sin_gravity

    REAL(8), INTENT(INOUT) :: gx(0:isize, 1:jsize, 1:ksize)
    REAL(8), INTENT(INOUT) :: gy(1:isize, 0:jsize, 1:ksize)
    REAL(8), INTENT(INOUT) :: gz(1:isize, 1:jsize, 0:ksize)

    REAL(8) :: hsp(0:isize,   0:jsize, 1:2)

    REAL(8) :: rho_bar(0:ksize+1)
    REAL(8) :: gz_bar(0:ksize)
    REAL(8) :: sz_bar(0:ksize)

    REAL(8) :: tmp
    REAL(8) :: gr

    INTEGER :: i, j, k
    INTEGER :: km, kn

    LOGICAL :: stat

    gr = gravity / rho_0

    CALL assert(buoyancy_scheme <= 1 .OR. .NOT. hydrostatic, "BUOYANCY_SCHEME>=2 is not supported in HYDROSTATIC mode")

    SELECT CASE(buoyancy_scheme)
    CASE (1) !gradient of hydrostatic-pressure
       hsp(:,:,:)=0.0

       CALL urecv(hsp)
!$OMP PARALLEL PRIVATE(km, tmp)
       DO k=ksize, 1, -1
          km = maskindex(k)
!$OMP DO
          DO j=1, jsize
          DO i=0, isize
             tmp = imask3d(i,j,km)*imask3d(i+1,j,km) * (rho(i+1,j,k)-rho(i,j,k))*idx1(i,j) * dz(k) * gr
             gx(i,j,k) = gx(i,j,k) - imask3d(i,j,km)*imask3d(i+1,j,km) * (hsp(i,j,1) + tmp*0.5D0)
             hsp(i,j,1) = imask3d(i,j,km)*imask3d(i+1,j,km) * (hsp(i,j,1) + tmp)
          END DO
          END DO

!$OMP DO
          DO j=0, jsize
          DO i=1, isize
             tmp = imask3d(i,j,km)*imask3d(i,j+1,km) * (rho(i,j+1,k)-rho(i,j,k))*idy1(i,j) * dz(k) * gr
             gy(i,j,k) = gy(i,j,k) - imask3d(i,j,km)*imask3d(i,j+1,km) * (hsp(i,j,2) + tmp*0.5D0)
             hsp(i,j,2) = imask3d(i,j,km)*imask3d(i,j+1,km) * (hsp(i,j,2) + tmp)
          END DO
          END DO
       END DO
!$OMP END PARALLEL
       CALL lsend(hsp)

    CASE (2) !explicit gravity forcing with reference density rho_0
!$OMP PARALLEL DO PRIVATE(km, kn)
       DO k=0, ksize
          km = maskindex(k)
          kn = maskindex(k+1)
          DO j=1, jsize
          DO i=1, isize
             gz(i,j,k) = gz(i,j,k) - imask3d(i,j,km)*imask3d(i,j,kn)   * (0.5*(rho(i,j,k) + rho(i,j,k+1))-rho_0) * gr * cos_gravity
          END DO
          END DO

          IF (k==0) CYCLE

          DO j=1, jsize
          DO i=0, isize
             gx(i,j,k) = gx(i,j,k) - imask3d(i,j,km)*imask3d(i+1,j,km) * (0.5*(rho(i,j,k) + rho(i+1,j,k))-rho_0) * gr * sin_gravity
          END DO
          END DO
       END DO

    CASE (3) !explicit gravity forcing with reference density rho_bar(k)
!$OMP PARALLEL DO PRIVATE(km)
       DO k=0, ksize+1
          km = maskindex(k)
          rho_bar(k) = sum(imask3d(1:isize,1:jsize,km)*dsz(1:isize,1:jsize)*rho(1:isize,1:jsize,k))
       END DO
       CALL hsum(rho_bar, all=.TRUE.)
       DO k=0, ksize+1
          IF (layer_area(k) == 0.0) THEN
             rho_bar(k) = rho_0
          ELSE
             rho_bar(k) = rho_bar(k) / layer_area(k)
          END IF
       END DO

!$OMP PARALLEL DO PRIVATE(km, kn)
       DO k=0, ksize
          km = maskindex(k)
          kn = maskindex(k+1)
          DO j=1, jsize
          DO i=1, isize
             gz(i,j,k) = gz(i,j,k) - imask3d(i,j,km)*imask3d(i,j,kn) * 0.5*(rho(i,j,k)-rho_bar(k) + rho(i,j,k+1)-rho_bar(k+1)) * gr
          END DO
          END DO
       END DO

    CASE (4) !explicit gravity forcing with reference density rho_ref(i,j,k)
       CALL buoyancy_explicit(gx, gy, gz)

    END SELECT

    IF (cycle_z) THEN
!$OMP PARALLEL DO
       DO k=0, ksize
          gz_bar(k) = sum(gz(1:isize,1:jsize,k)*dsz(1:isize,1:jsize)*imask3d(1:isize,1:jsize,k)*imask3d(1:isize,1:jsize,k+1))
          sz_bar(k) = sum(                      dsz(1:isize,1:jsize)*imask3d(1:isize,1:jsize,k)*imask3d(1:isize,1:jsize,k+1))
       END DO
       CALL hsum(gz_bar, all=.TRUE.)
       CALL hsum(sz_bar, all=.TRUE.)
!$OMP PARALLEL DO
       DO k=0, ksize
          gz_bar(k) =  gz_bar(k)/sz_bar(k)
          DO j=1, jsize
          DO i=1, isize
             gz(i,j,k) = (gz(i,j,k) - gz_bar(k))*imask3d(i,j,k)*imask3d(i,j,k+1)
          END DO
          END DO
       END DO
    END IF

  END SUBROUTINE buoyancy

  SUBROUTINE buoyancy_explicit(gx, gy, gz)
    USE geometry
    USE parameters, ONLY: buoyancy_scheme, rho_0, gravity, cos_gravity, sin_gravity

    REAL(8), INTENT(INOUT) :: gx(0:isize, 1:jsize, 1:ksize)
    REAL(8), INTENT(INOUT) :: gy(1:isize, 0:jsize, 1:ksize)
    REAL(8), INTENT(INOUT) :: gz(1:isize, 1:jsize, 0:ksize)

    REAL(8) :: drho(0:isize+1,0:jsize+1,0:ksize+1)
    REAL(8) :: tref(0:isize+1,0:jsize+1,0:ksize+1)
    REAL(8) :: sref(0:isize+1,0:jsize+1,0:ksize+1)

    REAL(8) :: gr

    INTEGER :: i, j, k, km, kn
    LOGICAL :: stat

    gr = gravity / rho_0

    CALL checkin('RHO_REF', drho, stat)
    IF (.NOT. stat) THEN
       CALL assert(eq_of_state==1, "BUOYANCY_SCHEME=4 requires 'RHO_REF' input")

       IF (tracer_index_t /= 0) THEN
          CALL checkin('T_REF', tref, stat)
          CALL assert(stat, "BUOYANCY_SCHEME=4 requires 'RHO_REF' (or 'T_REF' and 'S_REF') input")

          IF (cycle_x .AND. tracer_info(tracer_index_t)%offsetx/=0.0) THEN
             IF (icoord == 0)      tref(0,      :,:) = tref(0,      :,:) - tracer_info(tracer_index_t)%offsetx
             IF (icoord == ipes-1) tref(isize+1,:,:) = tref(isize+1,:,:) + tracer_info(tracer_index_t)%offsetx
          END IF
          IF (cycle_y .AND. tracer_info(tracer_index_t)%offsety/=0.0) THEN
             IF (jcoord == 0)      tref(:,0,      :) = tref(:,0,      :) - tracer_info(tracer_index_t)%offsety
             IF (jcoord == jpes-1) tref(:,jsize+1,:) = tref(:,jsize+1,:) + tracer_info(tracer_index_t)%offsety
          END IF
          IF (cycle_z .AND. tracer_info(tracer_index_t)%offsetz/=0.0) THEN
             IF (kcoord == 0)      tref(:,:,0)       = tref(:,:,0)       - tracer_info(tracer_index_t)%offsetz
             IF (kcoord == kpes-1) tref(:,:,ksize+1) = tref(:,:,ksize+1) + tracer_info(tracer_index_t)%offsetz
          END IF
       ELSE
          tref(:,:,:) = eos_t_ref
       END IF

       IF (tracer_index_s /= 0) THEN
          CALL checkin('S_REF', sref, stat)
          CALL assert(stat, "BUOYANCY_SCHEME=4 requires 'RHO_REF' (or 'T_REF' and 'S_REF') input")

          IF (tracer_info(tracer_index_s)%offsetx/=0.0) THEN
             IF (icoord == 0)      sref(0,      :,:) = sref(0,      :,:) - tracer_info(tracer_index_s)%offsetx
             IF (icoord == ipes-1) sref(isize+1,:,:) = sref(isize+1,:,:) + tracer_info(tracer_index_s)%offsetx
          END IF
          IF (tracer_info(tracer_index_s)%offsety/=0.0) THEN
             IF (jcoord == 0)      sref(:,0,      :) = sref(:,0,      :) - tracer_info(tracer_index_s)%offsety
             IF (jcoord == jpes-1) sref(:,jsize+1,:) = sref(:,jsize+1,:) + tracer_info(tracer_index_s)%offsety
          END IF
          IF (tracer_info(tracer_index_s)%offsetz/=0.0) THEN
             IF (kcoord == 0)      sref(:,:,0)       = sref(:,:,0)       - tracer_info(tracer_index_s)%offsetz
             IF (kcoord == kpes-1) sref(:,:,ksize+1) = sref(:,:,ksize+1) + tracer_info(tracer_index_s)%offsetz
          END IF
       ELSE
          sref(:,:,:) = eos_s_ref
       END IF

       DO k=0, ksize+1
       DO j=0, jsize+1
       DO i=0, isize+1
          drho(i,j,k) = eos_linear(tref(i,j,k), sref(i,j,k))
       END DO
       END DO
       END DO
    END IF
    CALL assert(stat, "BUOYANCY_SCHEME=4 requires 'RHO_REF' (or 'T_REF' and 'S_REF') input")

!$OMP PARALLEL WORKSHARE
    drho(:,:,:) = rho(0:isize+1,0:jsize+1,0:ksize+1) - drho(:,:,:)
!$OMP END PARALLEL WORKSHARE
    CALL checkout("DRHO", drho)

!$OMP PARALLEL DO PRIVATE(km, kn)
    DO k=0, ksize
       km = maskindex(k)
       kn = maskindex(k+1)
       DO j=1, jsize
       DO i=1, isize
          gz(i,j,k) = gz(i,j,k) - imask3d(i,j,km)*imask3d(i,j,kn)   * 0.5*(drho(i,j,k) + drho(i,j,k+1)) * gr * cos_gravity
       END DO
       END DO

       IF (k==0) CYCLE

       DO j=1, jsize
       DO i=0, isize
          gx(i,j,k) = gx(i,j,k) - imask3d(i,j,km)*imask3d(i+1,j,km) * 0.5*(drho(i,j,k) + drho(i+1,j,k)) * gr * sin_gravity
       END DO
       END DO
    END DO

  END SUBROUTINE buoyancy_explicit

!-----------------------------------------------------------------------------------------------------------------------

  REAL(8) PURE FUNCTION spiciness_flament(t, s)
    ! sea water spiciness by Flament, P., 2002, Prog. Oceanogr.

    REAL(8), INTENT(IN) :: t ! potential temperature
    REAL(8), INTENT(IN) :: s ! salinity

    REAL(8) :: t_(1:5), s_(1:4)

    REAL(8), PARAMETER :: b(0:5,0:4) = RESHAPE( &
         (/  0.0000D-0,  5.1655D-2,  6.64783D-3, -5.40230D-5,  3.9490D-7, -6.3600D-10, &
             7.7442D-1,  2.0340D-3, -2.46810D-4,  7.32600D-6, -3.0290D-8, -1.3090D-9 , &
            -5.8500D-3, -2.7420D-4, -1.42800D-5,  7.00360D-6, -3.8209D-7,  6.0480D-9 , &
            -9.8400D-4, -8.5000D-6,  3.33700D-5, -3.04120D-6,  1.0012D-7, -1.1409D-9 , &
            -2.0600D-4,  1.3600D-5,  7.89400D-6, -1.08530D-6,  4.7133D-8, -6.6760D-10 /), (/6,5/))

    t_(1) = t
    t_(2) = t_(1)*t_(1)
    t_(3) = t_(2)*t_(1)
    t_(4) = t_(3)*t_(1)
    t_(5) = t_(4)*t_(1)

    s_(1) = s-35.0
    s_(2) = s_(1)*s_(1)
    s_(3) = s_(2)*s_(1)
    s_(4) = s_(3)*s_(1)

    spiciness_flament = (b(0,0) + b(1,0)*t_(1) + b(2,0)*t_(2) + b(3,0)*t_(3) + b(4,0)*t_(4) + b(5,0)*t_(5))       &
                      + (b(0,1) + b(1,1)*t_(1) + b(2,1)*t_(2) + b(3,1)*t_(3) + b(4,1)*t_(4) + b(5,1)*t_(5))*s_(1) &
                      + (b(0,2) + b(1,2)*t_(1) + b(2,2)*t_(2) + b(3,2)*t_(3) + b(4,2)*t_(4) + b(5,2)*t_(5))*s_(2) &
                      + (b(0,3) + b(1,3)*t_(1) + b(2,3)*t_(2) + b(3,3)*t_(3) + b(4,3)*t_(4) + b(5,3)*t_(5))*s_(3) &
                      + (b(0,4) + b(1,4)*t_(1) + b(2,4)*t_(2) + b(3,4)*t_(3) + b(4,4)*t_(4) + b(5,4)*t_(5))*s_(4)
  END FUNCTION spiciness_flament

  REAL(8) PURE FUNCTION spiciness_linear(t, s)
    USE parameters, ONLY: rho_0, alpha=>eos_lin_alpha, beta=>eos_lin_beta, t0=>eos_t_ref, s0=>eos_s_ref

    REAL(8), INTENT(IN) :: t
    REAL(8), INTENT(IN) :: s

    spiciness_linear = rho_0 * (alpha * (t-t0) + beta * (s - s0))

  END FUNCTION spiciness_linear

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE mixed_layer_depth
    USE geometry
    USE parameters, ONLY: mld_dsigma

    REAL(4) :: csigma(1-slv:isize+slv,1-slv:jsize+slv)
    REAL(4) :: r

    INTEGER :: i, j, k

    IF (mld_dsigma == UNDEF) RETURN

    IF (kcoord==kpes-1) THEN
!$OMP PARALLEL DO
       DO j=1-slv, jsize+slv
       DO i=1-slv, isize+slv
          csigma(i,j) = imask3d(i,j,ksize)*(sigma0(i,j,ksize)+mld_dsigma)
       END DO
       END DO
    END IF

    CALL vcast(csigma)

    mld(:,:) = 0.0

!$OMP PARALLEL
    DO k=1, ksize
!$OMP DO
       DO j=1-slv, jsize+slv
       DO i=1-slv, isize+slv
          IF (mld(i,j) > depth(k-1) + dz(k-1)*0.5) CYCLE
          IF (.NOT. lmask3d(i,j,k)) CYCLE
          IF (sigma0(i,j,k) > csigma(i,j)) CYCLE

          IF (lmask3d(i,j,k-1) .AND. sigma0(i,j,k-1) > csigma(i,j)) THEN
             r = (csigma(i,j) - sigma0(i,j,k)) / (sigma0(i,j,k-1) - sigma0(i,j,k))
             r = min(max(0.0, r), 1.0)

             IF (r < 0.5) THEN
                mld(i,j) = max(mld(i,j), depth(k-1) - dz(k)*(0.5-r))
             ELSE
                mld(i,j) = max(mld(i,j), depth(k-1) + dz(k-1)*(r-0.5))
             END IF
          ELSE
             mld(i,j) = max(mld(i,j), depth(k-1))
          END IF
       END DO
       END DO
    END DO
!$OMP END PARALLEL

    CALL vmax(mld, all=.TRUE.)

!$OMP PARALLEL DO
    DO j=1-slv, jsize+slv
    DO i=1-slv, isize+slv
       mld(i,j) = min(mld(i,j), -h_bathymetry(i,j))
    END DO
    END DO

  END SUBROUTINE mixed_layer_depth

END MODULE state
