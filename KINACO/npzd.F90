#include "macro.h"

MODULE npzd
  USE misc
  USE parameters, ONLY: use_npzd, day_sec
  IMPLICIT NONE
  PRIVATE
  PUBLIC init_npzd, step_npzd
  PUBLIC use_convrate, convrate_n2p, convrate_d2n, &
                       convrate_p2n, convrate_p2z, convrate_p2d, &
                       convrate_z2n, convrate_z2d
  PUBLIC PON_SVn
  PUBLIC tracer_index_nut, tracer_index_phy, tracer_index_zoo, tracer_index_det

  REAL(8), SAVE :: PHY_Vmax  = UNDEF
  REAL(8), SAVE :: PHY_KNO3  = UNDEF
  REAL(8), SAVE :: PHY_KGpp  = UNDEF
  REAL(8), SAVE :: PHY_Mor0  = UNDEF
  REAL(8), SAVE :: PHY_KMor  = UNDEF
  REAL(8), SAVE :: PHY_Res0  = UNDEF
  REAL(8), SAVE :: PHY_KRes  = UNDEF
  REAL(8), SAVE :: PHY_Iopt  = UNDEF
  REAL(8), SAVE :: ZOO_GRmax1= UNDEF
  REAL(8), SAVE :: ZOO_KGraz = UNDEF
  REAL(8), SAVE :: ZOO_Lmd   = UNDEF
  REAL(8), SAVE :: ZOO_Star1 = UNDEF
  REAL(8), SAVE :: ZOO_Alpha = UNDEF
  REAL(8), SAVE :: ZOO_Beta  = UNDEF
  REAL(8), SAVE :: ZOO_Mor0  = UNDEF
  REAL(8), SAVE :: ZOO_KMor  = UNDEF
  REAL(8), SAVE :: PON_VP2N0 = UNDEF
  REAL(8), SAVE :: PON_KP2N  = UNDEF
  REAL(8), SAVE :: PON_SVn   = UNDEF
  REAL(8), SAVE :: Lc_alpha1 = UNDEF
  REAL(8), SAVE :: Lc_alpha2 = UNDEF

  REAL(4), SAVE :: PARfrac   = 0.45 !scale factor to convert SW to PAR

  INTEGER :: tracer_index_nut = 0
  INTEGER :: tracer_index_phy = 0
  INTEGER :: tracer_index_zoo = 0
  INTEGER :: tracer_index_det = 0

  REAL(8), ALLOCATABLE :: convrate_n2p(:,:,:)
  REAL(8), ALLOCATABLE :: convrate_d2n(:,:,:)
  REAL(8), ALLOCATABLE :: convrate_p2n(:,:,:)
  REAL(8), ALLOCATABLE :: convrate_p2z(:,:,:)
  REAL(8), ALLOCATABLE :: convrate_p2d(:,:,:)
  REAL(8), ALLOCATABLE :: convrate_z2n(:,:,:)
  REAL(8), ALLOCATABLE :: convrate_z2d(:,:,:)

  LOGICAL, SAVE :: use_convrate = .FALSE.

CONTAINS
  REAL(8) PURE FUNCTION Td(a, b, temp)
    REAL(8), INTENT(IN) :: a, b, temp
    Td     = a * exp(b*Temp)
  END FUNCTION Td

  REAL(8) PURE FUNCTION GraF(a, b, c)
    REAL(8), INTENT(IN) :: a, b, c
    GraF   = max(0.0D0, 1.0D0 - exp(a * (b - c)))
  END FUNCTION GraF

  REAL(8) PURE FUNCTION Mich(a,b)
    REAL(8), INTENT(IN) :: a, b
    Mich   = b / ( a + b )
  END FUNCTION Mich

!-------------------------------------------------------------------

  SUBROUTINE init_npzd
    USE geometry
    USE tracers

    INTEGER        :: iostat
    CHARACTER(256) :: iomsg

! paramset=1 from Onitsuka and Yanagi (2005) for the Sea of Japan (DEFAULT)
    TYPE npzd_params_struct
       CHARACTER(128) :: reference = "Onitsuka and Yanagi, 2005"
   !--- phytoplankkton parameters -----------------------------------------------------------------
       REAL(8) :: PHY_Vmax   = 0.6D0/day_sec   ! Maximum Photosynthetic rate @0degC   [/s]
       REAL(8) :: PHY_KNO3   = 1.5D-6          ! Half satuation constant for Nitrate  [molN/l]
       REAL(8) :: PHY_KGpp   = 6.93D-2         ! Temp. Coeff. for Photosynthetic Rate [/degC]
       REAL(8) :: PHY_Mor0   = 7.0D4/day_sec   ! Mortality Rate @0degC                [l/molN/s]
       REAL(8) :: PHY_KMor   = 6.93D-2         ! Temp. Coeff. for Mortality           [/degC]
       REAL(8) :: PHY_Res0   = 0.03D0/day_sec  ! Respiration Rate at @0degC           [/s]
       REAL(8) :: PHY_KRes   = 0.0519D0        ! Temp. Coeff. for Respiration         [/degC]
       REAL(8) :: PHY_Iopt   = 70.D0/697.675   ! Optimum Light Intensity              [ly/min]

    !--- zooplankton parameters -----------------------------------------------------------------
       REAL(8) :: ZOO_GRmax1 = 0.30D0/day_sec  ! Maximum Rate of Grazing PS @0degC    [/s]
       REAL(8) :: ZOO_KGraz  = 6.93D-2         ! Temp. Coeff. for Grazing             [/degC]
       REAL(8) :: ZOO_Lmd    = 1.40D6          ! Ivlev constant                       [l/molN]
       REAL(8) :: ZOO_Star1  = 4.300D-8        ! Threshold Value for Grazing PS       [molN/l]
       REAL(8) :: ZOO_Alpha  = 0.70D0          ! Assimilation Efficiency              [(nodim)]
       REAL(8) :: ZOO_Beta   = 0.30D0          ! Growth Efficiency                    [(nodim)]
       REAL(8) :: ZOO_Mor0   = 7.00D4/day_sec  ! Mortality Rate @0degC                [l/molN/s]
       REAL(8) :: ZOO_KMor   = 0.0693D0        ! Temp. Coeff. for Mortality           [/degC]

    !--- sinking-particle parameters ------------------------------------------------------------
       REAL(8) :: PON_VP2N0  = 0.05D0/day_sec ! PON Decomp. Rate @0degC               [/s]
       REAL(8) :: PON_KP2N   = 6.93D-2        ! PON Temp. Coeff. for Decomp.          [/degC]
       REAL(8) :: PON_SVn    = 10.0D0/day_sec ! Settling velocity of PON              [m/s]

       REAL(8) :: Lc_alpha1 = 5.0D-2          ! Light Dissipation coefficient         [/m]
       REAL(8) :: Lc_alpha2 = 6.0D4           ! PS Selfshading coefficient            [l/molN/m]
    END TYPE npzd_params_struct

    INTEGER :: paramset = 1
    TYPE(npzd_params_struct) :: params(3)

    NAMELIST / npzd / &
         use_npzd,  &
         paramset,  &
         PHY_Vmax,  &
         PHY_KNO3,  &
         PHY_KGpp,  &
         PHY_Mor0,  &
         PHY_KMor,  &
         PHY_Res0,  &
         PHY_KRes,  &
         PHY_Iopt,  &
         ZOO_GRmax1,&
         ZOO_KGraz, &
         ZOO_Lmd,   &
         ZOO_Star1, &
         ZOO_Alpha, &
         ZOO_Beta,  &
         ZOO_Mor0,  &
         ZOO_KMor,  &
         PON_VP2N0, &
         PON_KP2N,  &
         PON_SVn,   &
         Lc_alpha1, &
         Lc_alpha2, &
         parfrac

! paramset=2 from Yoshikawa et al. (2005) for the subarctic region.
    params(2)%reference = "Yoshikawa et al., 2005"
    params(2)%PHY_Vmax  = 0.5D0/day_sec
    params(2)%PHY_KNO3  = 2.0D-6
    params(2)%PHY_Mor0  = 4.3785D4/day_sec
    params(2)%PHY_Iopt  = 104.7D0/697.675
    params(2)%ZOO_Star1 = 4.000D-8
    params(2)%ZOO_Mor0  = 5.85D4/day_sec
    params(2)%PON_VP2N0 = 0.10D0/day_sec
    params(2)%PON_SVn   = 20.0D0/day_sec
    params(2)%Lc_alpha1 = 4.0D-2
    params(2)%Lc_alpha2 = 4.0D4

! paramset=3 from Kawamiya et al. (2000) for global run.
    params(3)%reference = "Kawamiya et al., 2000"
    params(3)%PHY_Vmax  = 1.2D0/day_sec
    params(3)%PHY_KNO3  = 0.03D-6
    params(3)%PHY_KGpp  = 6.3D-2
    params(3)%PHY_Mor0  = 2.81D4/day_sec
    params(3)%PHY_KMor  = 6.9D-2
    params(3)%PHY_Iopt  = 0.07D0
    params(3)%ZOO_Mor0  = 5.85D4/day_sec
    params(3)%PON_SVn   = 20.0D0/day_sec
    params(3)%Lc_alpha1 = 3.5D-2
    params(3)%Lc_alpha2 = 2.81D4


    IF (rank==0) THEN
       REWIND(CONFIG_UNIT)
       READ(CONFIG_UNIT, NML=npzd, IOSTAT=iostat, IOMSG=iomsg)
       CALL assert(iostat <= 0, "failed to read NPZD namelist", iomsg)
    END IF

    CALL bcast(use_npzd)
    IF (.NOT. use_npzd) RETURN

    IF (rank==0) THEN
       CALL assert(paramset >= 1 .AND. paramset <= size(params), "Invalid NPZD PARAMSET")

       IF (PHY_Vmax   == UNDEF) PHY_Vmax   = params(paramset)%PHY_Vmax
       IF (PHY_KNO3   == UNDEF) PHY_KNO3   = params(paramset)%PHY_KNO3
       IF (PHY_KGpp   == UNDEF) PHY_KGpp   = params(paramset)%PHY_KGpp
       IF (PHY_Mor0   == UNDEF) PHY_Mor0   = params(paramset)%PHY_Mor0
       IF (PHY_KMor   == UNDEF) PHY_KMor   = params(paramset)%PHY_KMor
       IF (PHY_Res0   == UNDEF) PHY_Res0   = params(paramset)%PHY_Res0
       IF (PHY_KRes   == UNDEF) PHY_KRes   = params(paramset)%PHY_KRes
       IF (PHY_Iopt   == UNDEF) PHY_Iopt   = params(paramset)%PHY_Iopt
       IF (ZOO_GRmax1 == UNDEF) ZOO_GRmax1 = params(paramset)%ZOO_GRmax1
       IF (ZOO_KGraz  == UNDEF) ZOO_KGraz  = params(paramset)%ZOO_KGraz
       IF (ZOO_Lmd    == UNDEF) ZOO_Lmd    = params(paramset)%ZOO_Lmd
       IF (ZOO_Star1  == UNDEF) ZOO_Star1  = params(paramset)%ZOO_Star1
       IF (ZOO_Alpha  == UNDEF) ZOO_Alpha  = params(paramset)%ZOO_Alpha
       IF (ZOO_Beta   == UNDEF) ZOO_Beta   = params(paramset)%ZOO_Beta
       IF (ZOO_Mor0   == UNDEF) ZOO_Mor0   = params(paramset)%ZOO_Mor0
       IF (ZOO_KMor   == UNDEF) ZOO_KMor   = params(paramset)%ZOO_KMor
       IF (PON_VP2N0  == UNDEF) PON_VP2N0  = params(paramset)%PON_VP2N0
       IF (PON_KP2N   == UNDEF) PON_KP2N   = params(paramset)%PON_KP2N
       IF (PON_SVn    == UNDEF) PON_SVn    = params(paramset)%PON_SVn
       IF (Lc_alpha1  == UNDEF) Lc_alpha1  = params(paramset)%Lc_alpha1
       IF (Lc_alpha2  == UNDEF) Lc_alpha2  = params(paramset)%Lc_alpha2

       WRITE(REPORT_UNIT, '(A,I0,A)') "NPZD module is enabled: paramset=",paramset," ("//trim(params(paramset)%reference)//")"
#ifdef DEBUG
       WRITE(REPORT_UNIT, '(A,ES10.3)') "   PHY_Vmax   = ", PHY_Vmax
       WRITE(REPORT_UNIT, '(A,ES10.3)') "   PHY_KNO3   = ", PHY_KNO3
       WRITE(REPORT_UNIT, '(A,ES10.3)') "   PHY_KGpp   = ", PHY_KGpp
       WRITE(REPORT_UNIT, '(A,ES10.3)') "   PHY_Mor0   = ", PHY_Mor0
       WRITE(REPORT_UNIT, '(A,ES10.3)') "   PHY_KMor   = ", PHY_KMor
       WRITE(REPORT_UNIT, '(A,ES10.3)') "   PHY_Res0   = ", PHY_Res0
       WRITE(REPORT_UNIT, '(A,ES10.3)') "   PHY_KRes   = ", PHY_KRes
       WRITE(REPORT_UNIT, '(A,ES10.3)') "   PHY_Iopt   = ", PHY_Iopt
       WRITE(REPORT_UNIT, '(A,ES10.3)') "   ZOO_GRmax1 = ", ZOO_GRmax1
       WRITE(REPORT_UNIT, '(A,ES10.3)') "   ZOO_KGraz  = ", ZOO_KGraz
       WRITE(REPORT_UNIT, '(A,ES10.3)') "   ZOO_Lmd    = ", ZOO_Lmd
       WRITE(REPORT_UNIT, '(A,ES10.3)') "   ZOO_Star1  = ", ZOO_Star1
       WRITE(REPORT_UNIT, '(A,ES10.3)') "   ZOO_Alpha  = ", ZOO_Alpha
       WRITE(REPORT_UNIT, '(A,ES10.3)') "   ZOO_Beta   = ", ZOO_Beta
       WRITE(REPORT_UNIT, '(A,ES10.3)') "   ZOO_Mor0   = ", ZOO_Mor0
       WRITE(REPORT_UNIT, '(A,ES10.3)') "   ZOO_KMor   = ", ZOO_KMor
       WRITE(REPORT_UNIT, '(A,ES10.3)') "   PON_VP2N0  = ", PON_VP2N0
       WRITE(REPORT_UNIT, '(A,ES10.3)') "   PON_KP2N   = ", PON_KP2N
       WRITE(REPORT_UNIT, '(A,ES10.3)') "   PON_SVn    = ", PON_SVn
       WRITE(REPORT_UNIT, '(A,ES10.3)') "   Lc_alpha1  = ", Lc_alpha1
       WRITE(REPORT_UNIT, '(A,ES10.3)') "   Lc_alpha2  = ", Lc_alpha2
       WRITE(REPORT_UNIT, '(A,ES10.3)') "   PARfrac    = ", PARfrac
#endif
    END IF

    CALL bcast(PHY_Vmax)
    CALL bcast(PHY_KNO3)
    CALL bcast(PHY_KGpp)
    CALL bcast(PHY_Mor0)
    CALL bcast(PHY_KMor)
    CALL bcast(PHY_Res0)
    CALL bcast(PHY_KRes)
    CALL bcast(PHY_Iopt)
    CALL bcast(ZOO_GRmax1)
    CALL bcast(ZOO_KGraz)
    CALL bcast(ZOO_Lmd)
    CALL bcast(ZOO_Star1)
    CALL bcast(ZOO_Alpha)
    CALL bcast(ZOO_Beta)
    CALL bcast(ZOO_Mor0)
    CALL bcast(ZOO_KMor)
    CALL bcast(PON_VP2N0)
    CALL bcast(PON_KP2N)
    CALL bcast(PON_SVn)
    CALL bcast(Lc_alpha1)
    CALL bcast(Lc_alpha2)
    CALL bcast(PARfrac)

    CALL add_tracer('NUT')
    CALL add_tracer('PHY')
    CALL add_tracer('ZOO')
    CALL add_tracer('DET')

    tracer_index_nut = tracer_index('NUT')
    tracer_index_phy = tracer_index('PHY')
    tracer_index_zoo = tracer_index('ZOO')
    tracer_index_det = tracer_index('DET')

    CALL assert(tracer_index_nut /= 0, "failed to introduce tracer 'NUT' for the NPZD module")
    CALL assert(tracer_index_phy /= 0, "failed to introduce tracer 'PHY' for the NPZD module")
    CALL assert(tracer_index_zoo /= 0, "failed to introduce tracer 'ZOO' for the NPZD module")
    CALL assert(tracer_index_det /= 0, "failed to introduce tracer 'DET' for the NPZD module")

  END SUBROUTINE init_npzd

!-------------------------------------------------------------------

  SUBROUTINE step_npzd
    USE geometry
    USE tracers
    USE io

    REAL(8) :: tmp(1:isize, 1:jsize)

    REAL(8) :: dnut(1:isize, 1:jsize, 1:ksize)
    REAL(8) :: dphy(1:isize, 1:jsize, 1:ksize)
    REAL(8) :: dzoo(1:isize, 1:jsize, 1:ksize)
    REAL(8) :: ddet(1:isize, 1:jsize, 1:ksize)

!     ..... OUTPUT parameters .....
    REAL(8) :: NtoPZD(1:isize, 1:jsize, 1:ksize)
    REAL(8) :: Grazing(1:isize, 1:jsize, 1:ksize)

!     ...... Light condition parameters .......
    REAL(8) :: Lint(1:isize, 1:jsize, 0:ksize)    ! Light Intencity [ly/min]
    REAL(8) :: LfcS(1:isize, 1:jsize, 1:ksize)    ! Light factor for PHY

!     ...... MASK for NPZD cycle .......
    REAL(4) :: npzdmask(1:isize, 1:jsize, 1:ksize)

!     ...... Biological parameters ......
    REAL(8) :: GppPSn,  GppNPS
    REAL(8) :: ResPSn
    REAL(8) :: MorPSn
    REAL(8) :: GraPS2ZSn
    REAL(8) :: EgeZSn, MorZSn, ExcZSn
    REAL(8) :: DecP2N

!     ...... Internal flux condition parameters .......
    REAL(8) :: tnut, tphy, tzoo, tpon
    REAL(8) :: temp, kappa

    REAL(8) :: l
    INTEGER :: i, j, k
    LOGICAL :: stat

    IF (use_convrate) THEN
       IF (.NOT. ALLOCATED(convrate_n2p)) THEN
          ALLOCATE(convrate_n2p(isize,jsize,ksize))
          convrate_n2p(:,:,:) = 0.0
       END IF

       IF (.NOT. ALLOCATED(convrate_p2n)) THEN
          ALLOCATE(convrate_p2n(isize,jsize,ksize))
          convrate_p2n(:,:,:) = 0.0
       END IF
       IF (.NOT. ALLOCATED(convrate_p2z)) THEN
          ALLOCATE(convrate_p2z(isize,jsize,ksize))
          convrate_p2z(:,:,:) = 0.0
       END IF
       IF (.NOT. ALLOCATED(convrate_p2d)) THEN
          ALLOCATE(convrate_p2d(isize,jsize,ksize))
          convrate_p2d(:,:,:) = 0.0
       END IF
       IF (.NOT. ALLOCATED(convrate_z2n)) THEN
          ALLOCATE(convrate_z2n(isize,jsize,ksize))
          convrate_z2n(:,:,:) = 0.0
       END IF
       IF (.NOT. ALLOCATED(convrate_z2d)) THEN
          ALLOCATE(convrate_z2d(isize,jsize,ksize))
          convrate_z2d(:,:,:) = 0.0
       END IF
       IF (.NOT. ALLOCATED(convrate_d2n)) THEN
          ALLOCATE(convrate_d2n(isize,jsize,ksize))
          convrate_d2n(:,:,:) = 0.0
       END IF
    END IF

!******************************************************************************
!     Tendency of conpartments of ecosystem
!******************************************************************************

! Solar Radiation @ Sea Surface
! Conversion from Shortwave Rad.to PAR (from Okunishi"s code)
! ALO(L)   = AL(L,ntime) * 0.45
!
! ATTENTION: swab data unit [W/m^2] : regional model
!                    [erg/cm^2/sec] : Global model
!
! The parameters are tuned so as to use [ly/min] unit.
! [ly/min] = [cal/cm^2/min]
!          = 4.1865*10^4/60 [J/s/m^2]
!          = 697.675 [W/m^2]
!
    tmp(:,:) = 0.0D0
    CALL checkin("RADPA", tmp, stat) !Photosynthetically avairable radiation (PAR) [W/m^2]

    IF (.NOT. stat) THEN
       CALL checkin('SFC_LINT', tmp, stat) !for backward compatibility
       CALL scaleoffset(tmp, scale=parfrac, offset=0.0)
    END IF

    IF (.NOT. stat) THEN
       CALL checkin("RADSW", tmp, stat)
       CALL scaleoffset(tmp, scale=(1.0-albedo_ocean)*parfrac, offset=0.0) ! converted from RADSW if RADPA is not given
    END IF

    CALL assert(stat, "NPZD is enabled but Photosynthetically avairable radiation (PAR) is not given")

    CALL checkin('NPZD_MASK', npzdmask, stat)
    IF (.NOT. stat) npzdmask(:,:,:) = 1.0


    IF (vrank==0) THEN
       DO j = 1, jsize
       DO i = 1, isize
          ! Unit conversion
          Lint(i,j,ksize) = max(tmp(i,j), 0.0D0) / 697.675d0       ! [W/m^2] ==> [ly/min], clip negative to zero
       END DO
       END DO
    ELSE
       CALL urecv(Lint(:,:,ksize))
    END IF

    DO k = ksize, 1, -1
       DO j = 1, jsize
       DO i = 1, isize
          Kappa = Lc_alpha1 + Lc_alpha2 * tracer(i,j,k, tracer_index_phy)*1.0D-3 *imask3d(I, J, K) !(mol/m3)=>(mol/l)
          Lint(i,j,k-1) = Lint(i,j,k)*exp(-Kappa*dz(k))

          l = 0.5*(Lint(i,j,k)+Lint(i,j,k-1))
          LfcS(i,j,k) = l/PHY_Iopt * exp(1.0D0 - l/PHY_Iopt )
       END DO
       END DO
    END DO

    CALL lsend(Lint(:,:,0))

!$OMP DO PRIVATE(tNUT, tPHY, tZOO, tPON, temp, GppNPS, GppPSn, ResPSn, MorPSn, GraPS2ZSn, ExcZsn, EgeZSn, MorZSn, DeCP2N)
    DO k = 1, ksize
    DO j = 1, jsize
    DO i = 1, isize
       IF (.NOT. lmask3d(i,j,k) .OR. npzdmask(i,j,k)==0.0) THEN
          !DO NOTHING
          dnut(i,j,k)    = 0.0
          dphy(i,j,k)    = 0.0
          dzoo(i,j,k)    = 0.0
          ddet(i,j,k)    = 0.0
          NtoPZD(i,j,k)  = 0.0
          Grazing(i,j,k) = 0.0
          CYCLE
       END IF

       ! unit conversion (mol/m^3)=>(mol/l)
       tNUT  =  max(tracer(i,j,k,tracer_index_nut)*1.0D-3, 0.0D0)
       tPHY  =  max(tracer(i,j,k,tracer_index_phy)*1.0D-3, 0.0D0)
       tZOO  =  max(tracer(i,j,k,tracer_index_zoo)*1.0D-3, 0.0D0)
       tPON  =  max(tracer(i,j,k,tracer_index_det)*1.0D-3, 0.0D0)

       temp = tracer(i,j,k,tracer_index_t)

       !   ..... Photosynthesis of PHY .....
       GppNPS    = Mich( PHY_KNO3, tNUT )
       GppPSn    = Td( PHY_Vmax, PHY_KGpp, temp) * LfcS(I, J, K) * tPHY * GppNPS
       ResPSn    = Td( PHY_Res0, PHY_KRes, temp) * tPHY
       MorPSn    = Td( PHY_Mor0, PHY_KMor, temp) * tPHY * tPHY

       !        ..... Grazing PHY --> ZOO .....
       GraPS2ZSn = Td( ZOO_GRmax1,ZOO_KGraz, temp ) * GraF(ZOO_Lmd, ZOO_Star1, tPHY) * tZOO

       !        ..... Mortality, Excration, Egestion for ZOO
       ExcZSn = (ZOO_Alpha- ZOO_Beta )  *  GraPS2ZSn
       EgeZSn = (1.0      - ZOO_Alpha)  *  GraPS2ZSn
       MorZSn = Td( ZOO_Mor0, ZOO_KMor, Temp ) * tZOO * tZOO

       !        ..... Decomposition DET ---> NUT .....
       DecP2N    = Td(PON_VP2N0, PON_KP2N, Temp ) * tPON

       IF (use_convrate) THEN
          IF (tNUT > 0.0) THEN
             convrate_n2p(i,j,k) = min(GppPSn / tNUT, 1.0D0) * dtime
          ELSE
             convrate_n2p(i,j,k) = 0.0D0
          END IF

          IF (tPON > 0.0) THEN
             convrate_d2n(i,j,k) = min(DecP2N / tPON, 1.0D0) * dtime
          ELSE
             convrate_d2n(i,j,k) = 0.0D0
          END IF

          IF (tPHY > 0.0) THEN
             convrate_p2n(i,j,k) = min(ResPSn            / tPHY, 1.0D0) * dtime
             convrate_p2z(i,j,k) = min(GraPS2ZSn         / tPHY, 1.0D0) * dtime
             convrate_p2d(i,j,k) = min(MorPSn            / tPHY, 1.0D0) * dtime
          ELSE
             convrate_p2n(i,j,k) = 0.0D0
             convrate_p2z(i,j,k) = 0.0D0
             convrate_p2d(i,j,k) = 0.0D0
          END IF

          IF (tZOO > 0.0) THEN
             convrate_z2n(i,j,k) = min(ExcZSn            / tZOO, 1.0D0) * dtime
             convrate_z2d(i,j,k) = min((MorZSn + EgeZsn) / tZOO, 1.0D0) * dtime
          ELSE
             convrate_z2n(i,j,k) = 0.0D0
             convrate_z2d(i,j,k) = 0.0D0
          END IF
       END IF

       !unit conversion (mol/l) => (mol/m^3)
       GppPSn = GppPSn * 1.0D3
       ResPSn = ResPSn * 1.0D3
       DecP2N = DecP2N * 1.0D3
       ExcZSn = ExcZSn * 1.0D3
       MorPSn = MorPSn * 1.0D3
       MorZSn = MorZSn * 1.0D3
       EgeZSn = EgeZSn * 1.0D3
       GraPS2ZSn = GraPS2ZSn * 1.0D3

       !.... Tendency terms for biological processes .....
       dnut(i,j,k) = -(GppPSn - ResPSn) + DecP2N + ExcZSn
       dphy(i,j,k) = GppPSn - ResPSn - MorPSn - GraPS2ZSn
       dzoo(i,j,k) = GraPS2ZSn - MorZSn - ExcZSn - EgeZSn
       ddet(i,j,k) = MorPSn + MorZSn + EgeZSn - DecP2N

       !-----------Net primary production from NUT---------------------------
       NtoPZD(i,j,k) = GppPSn - ResPSn

       !-----------Grazing---------------------------------------------------
       Grazing(i,j,k) = GraPS2ZSn


    END DO
    END DO
    END DO

!******************************************************************************
!     Sinking of PON
!******************************************************************************

    !        ..... settling fluxes to under grid .....

    IF (tracer_info(tracer_index_det)%w==0.0 .AND. .NOT. require_checkin('W_DET')) THEN

        tmp(:,:) = 0.0D0

        CALL urecv(tmp)

        DO k = ksize, 1, -1
!$OMP PARALLEL DO
           DO j = 1, jsize
           DO i = 1, isize
              ddet(i,j,k) = ddet(i,j,k) + tmp(i,j) / dz_star(i,j,k)

              tmp(i,j) = PON_SVn * max(tracer(i,j,k,tracer_index_det), 0.0D0) * imask3d(i,j,k)*imask3d(i,j,k-1)

              ddet(i,j,k) = ddet(i,j,k) - tmp(i,j) / dz_star(i,j,k)
           END DO
           END DO
        END DO

        CALL lsend(tmp)
    END IF

!$OMP PARALLEL DO
    DO k=1, ksize
    DO j=1, jsize
    DO i=1, isize
       tracer(i,j,k,tracer_index_nut) = tracer(i,j,k,tracer_index_nut) + dnut(i,j,k) * dtime
       tracer(i,j,k,tracer_index_phy) = tracer(i,j,k,tracer_index_phy) + dphy(i,j,k) * dtime
       tracer(i,j,k,tracer_index_zoo) = tracer(i,j,k,tracer_index_zoo) + dzoo(i,j,k) * dtime
       tracer(i,j,k,tracer_index_det) = tracer(i,j,k,tracer_index_det) + ddet(i,j,k) * dtime
    END DO
    END DO
    END DO

    CALL update_tracer_boundary(tracer_index_nut)
    CALL update_tracer_boundary(tracer_index_phy)
    CALL update_tracer_boundary(tracer_index_zoo)
    CALL update_tracer_boundary(tracer_index_det)

    CALL checkout("DNUT", dnut)
    CALL checkout("DPHY", dphy)
    CALL checkout("DZOO", dzoo)
    CALL checkout("DDET", ddet)

    CALL checkout("Grazing", Grazing)
    CALL checkout("NPP",     NtoPZD)
    CALL checkout("Lintensity", Lint)

    IF (use_convrate) THEN
       CALL checkout("CONVRATE_N2P", convrate_n2p)
       CALL checkout("CONVRATE_P2N", convrate_p2n)
       CALL checkout("CONVRATE_P2Z", convrate_p2z)
       CALL checkout("CONVRATE_P2D", convrate_p2d)
       CALL checkout("CONVRATE_Z2N", convrate_z2n)
       CALL checkout("CONVRATE_Z2D", convrate_z2d)
       CALL checkout("CONVRATE_D2N", convrate_d2n)
    END IF

  END SUBROUTINE step_npzd

END MODULE npzd
