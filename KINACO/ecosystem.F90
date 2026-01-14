#include "macro.h"

! marine ecosystem model based on NEMURO/MEM (Yamanaka et al. 2004 and Shigemitsu et al. 2012)

MODULE ecosystem
  USE misc
  USE parameters, ONLY: use_ecosystem, day_sec
  IMPLICIT NONE
  PRIVATE
  PUBLIC init_ecosystem, step_ecosystem

  INTEGER, PARAMETER :: n_phy = 2
  INTEGER, PARAMETER :: n_zoo = 3

  INTEGER, SAVE :: tracer_index_phy(n_phy)    ! phytoplankton (2 categories, 1: small, 2: diatom)
  INTEGER, SAVE :: tracer_index_zoo(n_zoo)    ! zooplankton   (3 categories, 1: micro, 2: meso, 3: predatory)
  INTEGER, SAVE :: tracer_index_no3           ! nitrate, NO3
  INTEGER, SAVE :: tracer_index_nh4           ! ammonium, NH3
  INTEGER, SAVE :: tracer_index_don           ! disolved organic nitrogen
  INTEGER, SAVE :: tracer_index_pon           ! particulate organic nitrogen
  INTEGER, SAVE :: tracer_index_dsi           ! disolved silicate, Si(OH)4
  INTEGER, SAVE :: tracer_index_psi           ! biogenic silica, opal (not include living diatoms)
  INTEGER, SAVE :: tracer_index_dfe           ! disolved iron
  INTEGER, SAVE :: tracer_index_pfe           ! particulate iron

  TYPE PHYTOPLANKTON_PARAMETER_STRUCT
     REAL(4) :: alpha       = 0.013    / day_sec  ! initial slope of phtosynthesis curve [m^2 / (W s)]
     REAL(4) :: beta        = 1.4E-15  / day_sec  ! photoinhibitation index              [m^2 / (W s)]
     REAL(4) :: gamma       = 0.135    / day_sec  ! ratio of extracellular excertion
     REAL(4) :: pmax        = 0.4      / day_sec  ! maximum light-saturated p-s rate     [1 / s]
     REAL(4) :: vmax        = 0.6      / day_sec  ! maximum growth rate  at 0degC        [1 / s]
     REAL(4) :: mort        = 0.0585E3 / day_sec  ! mortality rate at 0degC              [m^3 / (molN s)]
     REAL(4) :: resp        = 0.03     / day_sec  ! respiration rate at 0degC
     REAL(4) :: k_phot      = 0.0693              ! temp. coefficient for photosysthesis [1 / K]
     REAL(4) :: k_mort      = 0.0693              ! temp. coefficient for mortality      [1 / K]
     REAL(4) :: k_resp      = 0.0519              ! temp. coefficient for respilation    [1 / K]
     REAL(4) :: k_no3       = 1.0E-3              ! half-satulation const. for NO3       [molN  / m^3]
     REAL(4) :: k_nh4       = 0.1E-3              ! half-satulation const. for NH4       [molN  / m^3]
     REAL(4) :: k_dsi       = -1.0                ! half-satulation const. for Si(OH)4   [molSi / m^3]
     REAL(4) :: k_dfe       = -1.0                ! half-satulation const. for disol. Fe [molFe / m^3]
     REAL(4) :: a_no3       = 280                 ! maximum affinity for NO3             [m^3 / (molN s)]
     REAL(4) :: a_nh4                             ! maximum affinity for NH4             [m^3 / (molN s)]
     REAL(4) :: a_dsi                             ! maximum affinity for Si(OH)4         [m^3 / (molSi s)]
     REAL(4) :: a_dfe                             ! maximum affinity for disolved Fe     [m^3 / (molFe s)]
     REAL(4) :: lf0                               ! coefficient for light-limitation
  END type PHYTOPLANKTON_PARAMETER_STRUCT

  TYPE ZOOPLANKTON_PARAMETER_STRUCT
     REAL(4) :: alpha       = 0.7                 ! assimilation efficiency
     REAL(4) :: beta        = 0.3                 ! growth efficiency
     REAL(4) :: lambda      = 1.4E3               ! Ivlev constant                  [1 / molN]
     REAL(4) :: graz(n_phy) = 0.0                 ! maximum grazing rate at 0decC   [1 / s]
     REAL(4) :: pred(n_zoo) = 0.0                 ! maximum pradation rate at 0decC [1 / s]
     REAL(4) :: mort        = 0.0535E3 / day_sec  ! mortality rate at 0degC         [m^3 / (molN s)]
     REAL(4) :: k_gp        = 0.0693              ! temp. coefficient for grazing/pradation [1 / K]
     REAL(4) :: k_mort      = 0.0693              ! temp. coefficient for mortality [1 / K]
  END TYPE ZOOPLANKTON_PARAMETER_STRUCT

  TYPE(PHYTOPLANKTON_PARAMETER_STRUCT), SAVE :: phyparam(n_phy)
  TYPE(  ZOOPLANKTON_PARAMETER_STRUCT), SAVE :: zooparam(n_zoo)

  REAL(4) :: lint_a1 = 0.04E0            ! light dissipation coefficient   [1 / m]
  REAL(4) :: lint_a2 = 0.04E6            ! light shading by phytoplankton  [1 / m / molN]

  REAL(4) :: v_nitrif  = 0.03 / day_sec  ! nitrification rate at 0degC         [1 / s]
  REAL(4) :: k_nitrif  = 0.0693          ! temp. coefficient for nitrification [1 / k]

  REAL(4) :: v_donrem  = 0.15 / day_sec  ! DON remineralization rate at 0degC         [1 / s]
  REAL(4) :: k_donrem  = 0.0693          ! temp. coefficient for DON reminiralization [1 / k]

  REAL(4) :: v_ponrem  = 0.08 / day_sec  ! PON remineralization rate at 0degC         [1 / s]
  REAL(4) :: k_ponrem  = 0.0693          ! temp. coefficient for PON reminiralization [1 / k]

  REAL(4) :: v_pondec  = 0.08 / day_sec  ! PON decomposition rate at 0degC            [1 / s]
  REAL(4) :: k_pondec  = 0.0693          ! temp. coefficient for PON decomposition    [1 / k]

  REAL(4) :: v_psidis = 0.16 / day_sec   ! OPAL dissolution rate at 0degC             [1 / s]
  REAL(4) :: k_psidis = 0.0693           ! temp. coefficient for OPAL dissolution     [1 / k]

  REAL(4) :: gp_threshold = 0.043E-3     ! threshold for grazing/predation [molN / m^3]

  REAL(4) :: w_pon = -40.0  / day_sec    ! sinking velocity of PON  [m /s]
  REAL(4) :: w_psi = -40.0  / day_sec    ! sinking velocity of PSI  [m /s]
  REAL(4) :: w_pfe = -0.001 / day_sec    ! sinking velocity of PFE  [m /s]

  REAL(4) :: l_fedes = 0.003 / day_sec   ! Fe desorption rate at 30 degC
  REAL(4) :: a_fedes = 4000.0            ! slope of Arrhenius relation of Fe desorption [K]

  REAL(4) :: r_cn   = 6.625              ! stoichiometry of C  to N [molC  / molN]
  REAL(4) :: r_fen  = 1.7E-5             ! stoichiometry of Fe to N [molFe / molN]
  REAL(4) :: r_sinh = 1.0                ! stoichiometry of Si to N of diatoms in Fe-replete   condition [molSi / molN]
  REAL(4) :: r_sinl = 3.6                ! stoichiometry of Si to N of diatoms in Fe-deficient condition [molSi / molN]
  REAL(4) :: dfe_threshold = 0.03E-6     ! threshold of disolved Fe for the shift of r_shih/r_sinl

  REAL(4) :: c_dustfe = 0.035            ! ratio of iron content in dust
  REAL(4) :: a_dustfe = 0.04             ! solubility of iron in dust at the sea surface
  REAL(4) :: delta_harddust = 4000       ! dissolution length scale of hard dust
  REAL(4) :: delta_softdust =  600       ! dissolution length scale of soft dust
  REAL(4) :: f_harddust     = 0.97       ! fraction of hard dust
  REAL(4) :: lambda_fescav  = 18.5E6   ! base coefficient scavenging of Fe [m^2 / kg]
  REAL(4) :: gamma_fescav   = 4.40E3 / day_sec  ! propertionality constant for scavenging of Fe [m^3 / (molFe s)]

  REAL(4) :: c_ligand = 0.6E-6           ! total ligand concentration [mol / m^3]

CONTAINS

!-------------------------------------------------------------------

  SUBROUTINE init_ecosystem
    USE geometry
    USE parameters, ONLY: day_sec
    USE tracers

    INTEGER :: n

    INTEGER :: iostat
    CHARACTER(256) :: iomsg

    NAMELIST / ecosystem / &
         use_ecosystem

    IF (rank==0) THEN
       REWIND(CONFIG_UNIT)
       READ(CONFIG_UNIT, NML=ecosystem, IOSTAT=iostat, IOMSG=iomsg)
       CALL assert(iostat <= 0, "failed to read ECOSYSTEM namelist", iomsg)
    END IF

    CALL bcast(use_ecosystem)

    IF (.NOT. use_ecosystem) RETURN

    DO n=1, n_phy
       CALL add_tracer('PHY' // trim(format(n, 'I1')), tracer_index_phy(n))
    END DO
    DO n=1, n_zoo
       CALL add_tracer('ZOO' // trim(format(n, 'I1')), tracer_index_zoo(n))
    END DO

    CALL add_tracer('NO3',  tracer_index_no3)
    CALL add_tracer('NH4',  tracer_index_nh4)
    CALL add_tracer('DON',  tracer_index_don)
    CALL add_tracer('PON',  tracer_index_pon)
    CALL add_tracer('DSi',  tracer_index_dsi)
    CALL add_tracer('PSi',  tracer_index_psi)
    CALL add_tracer('DFe',  tracer_index_dfe)
    CALL add_tracer('PFe',  tracer_index_pfe)

    tracer_info(tracer_index_pon)%w = w_pon
    tracer_info(tracer_index_psi)%w = w_psi
    tracer_info(tracer_index_pfe)%w = w_pfe

! set parameters

! phy(1): small non-diatom phytoplankton
    phyparam(1)%alpha   = 0.013    / day_sec
    phyparam(1)%pmax    = 0.4      / day_sec
    phyparam(1)%vmax    = 0.6      / day_sec
    phyparam(1)%mort    = 0.0585E3 / day_sec
    phyparam(1)%resp    = 0.03     / day_sec
    phyparam(1)%k_no3   = 1.0E-3
    phyparam(1)%k_nh4   = 0.1E-3
    phyparam(1)%k_dfe   = 0.05E-6
    phyparam(1)%a_no3   = 282E-3

! phy(2): diatom
    phyparam(2)%alpha   = 0.045    / day_sec
    phyparam(2)%pmax    = 1.4      / day_sec
    phyparam(2)%vmax    = 0.8      / day_sec
    phyparam(2)%mort    = 0.029E3  / day_sec
    phyparam(1)%resp    = 0.03     / day_sec
    phyparam(2)%k_no3   = 3.0E-3
    phyparam(2)%k_nh4   = 0.3E-3
    phyparam(2)%k_dsi   = 6.0E-3
    phyparam(2)%k_dfe   = 0.1E-6
    phyparam(2)%a_no3   = 252E-3

    DO n=1, n_phy
       phyparam(n)%a_nh4 = phyparam(n)%a_no3 * phyparam(n)%k_no3/phyparam(n)%k_nh4
       phyparam(n)%a_dsi = phyparam(n)%a_no3 * phyparam(n)%k_no3/phyparam(n)%k_dsi
       phyparam(n)%a_dfe = phyparam(n)%a_no3 * phyparam(n)%k_no3/phyparam(n)%k_dfe

       phyparam(n)%lf0 = 1.0 / ( (phyparam(n)%alpha / (phyparam(n)%alpha + phyparam(n)%beta)) &
                                *(phyparam(n)%beta  / (phyparam(n)%alpha + phyparam(n)%beta))**(phyparam(n)%beta/phyparam(n)%alpha))
    END DO

! zoo(1) : micro-zooplankton
    zooparam(1)%graz(1) = 0.4 / day_sec

! zoo(2) : meso-zooplankto
    zooparam(2)%graz(1) = 0.1 / day_sec
    zooparam(2)%graz(2) = 0.4 / day_sec
    zooparam(2)%pred(1) = 0.4 / day_sec

! zoo(3) : predatory-zooplankton
    zooparam(3)%graz(2) = 0.2 / day_sec
    zooparam(3)%pred(1) = 0.2 / day_sec
    zooparam(3)%pred(2) = 0.4 / day_sec

  END SUBROUTINE init_ecosystem

!-------------------------------------------------------------------

  SUBROUTINE step_ecosystem
    USE geometry
    USE tracers
    USE io

    REAL(8) :: input(1:isize, 1:jsize, 2)
    REAL(8) :: lint(1:isize, 1:jsize, 0:ksize)
    REAL(8) :: fz_dustfe(1:isize, 1:jsize, 0:ksize)

    REAL(8) :: netpp(1:isize, 1:jsize, 1:ksize, n_phy)
    REAL(8) :: rnew( 1:isize, 1:jsize, 1:ksize, n_phy)

    REAL(8) :: phy(n_phy), zoo(n_zoo), no3, nh4, dsi, psi, dfe, pfe, don, pon, t
    REAL(8) :: phot(n_phy), fa, mu1, mu2, pmort(n_phy), resp(n_phy)
    REAL(8) :: graz(n_phy, n_zoo), pred(n_zoo, n_zoo), zmort(n_zoo), sumgp(n_zoo), expkt
    REAL(8) :: dno3, dnh4, dpon, ddon, dpfe, ddfe
    REAL(8) :: tmp
    REAL(4) :: r_sin

    INTEGER :: i, j, k, n, m

    !preference factor for ZOO3 [ m^3 / (molN) ]
    REAL(4) :: pref1 = 3.01E3
    REAL(4) :: pref2 = 4.605E3

    lint(:,:,ksize) = 0.0D0
    CALL urecv(lint(:,:,ksize))

!$OMP PARALLEL PRIVATE(i,j,k,n)
    DO k=ksize-1, 0, -1
!$OMP WORKSHARE
       lint(:,:,k) = lint(:,:,k+1) + lint_a1 * dz(k+1)
!$OMP END WORKSHARE
       DO n=1, n_phy
!$OMP DO
          DO j=1, jsize
          DO i=1, isize
             lint(i,j,k) = lint(i,j,k) + lint_a2 * tracer(i,j,k+1,tracer_index_phy(n)) * dz(k+1)
          END DO
          END DO
       END DO
    END DO
!$OMP END PARALLEL

    CALL lsend(lint(:,:,0))

    input(:,:,:) = 0.0D0
    CALL checkin('SFC_LINT', input(:,:,1)) ! surface light intensity [W / m^2]
    CALL checkin('SFC_DUST', input(:,:,2)) ! surface dust fall       [kg / (m^2 s)]

    input(:,:,2) = input(:,:,2) * c_dustfe / (molmass_fe)  ! iron in surface dust flux [molFe / (m^2 s)]

!$OMP PARALLEL DO
    DO k = 0, ksize
       DO j = 1, jsize
       DO i = 1, isize
          lint(i,j,k) = input(i,j,1)*exp(-lint(i,j,k))

          fz_dustfe(i,j,k) = input(i,j,2)*(1.0-a_dustfe) &
               * (f_harddust*exp(-depth(k)/delta_harddust) + (1.0-f_harddust)*exp(-depth(k)/delta_softdust))
       END DO
       END DO
    END DO

!$OMP PARALLEL DO &
!$OMP PRIVATE(phy, zoo, no3, nh4)               &
!&OMP PRIVATE(dsi, psi, dfe, pfe, don, pon, t)  &
!$OMP PRIVATE(phot, resp, pmort)                &
!$OMP PRIVATE(graz, pred, zmort)                &
!$OMP PRIVATE(fa, mu1, mu2, lint, expkt, sumgp) &
!$OMP PRIVATE(dno3, dnh4, dpon, ddon)           &
!$OMP PRIVATE(r_sin)
    DO k = 1, ksize
    DO j = 1, jsize
    DO i = 1, isize
       IF (.NOT. lmask3d(i,j,k)) CYCLE

       IF (use_landwater) THEN
          IF (lwflag(i,j)) CYCLE
       END IF

       no3 = tracer(i,j,k,tracer_index_no3)
       nh4 = tracer(i,j,k,tracer_index_nh4)
       dsi = tracer(i,j,k,tracer_index_dsi)
       dfe = tracer(i,j,k,tracer_index_dfe)
       don = tracer(i,j,k,tracer_index_pon)
       pon = tracer(i,j,k,tracer_index_don)
       pfe = tracer(i,j,k,tracer_index_pfe)
       psi = tracer(i,j,k,tracer_index_psi)
       t   = tracer(i,j,k,tracer_index_t)

       DO n=1, n_phy
          phy(n) = tracer(i,j,k,tracer_index_phy(n))

          fa = 1.0 / (1.0 + sqrt(max(phyparam(n)%a_no3*no3, phyparam(n)%a_nh4*nh4)/phyparam(n)%vmax))
          IF (phyparam(n)%a_dfe > 0) fa = max(fa, 1.0 / (1.0 + sqrt(phyparam(n)%a_dfe*dfe/phyparam(n)%vmax)))

          mu1 = phyparam(n)%vmax * no3 / (no3 / (1.0-fa) + phyparam(n)%vmax/(fa*phyparam(n)%a_no3)) * (1.0 - nh4/(nh4+phyparam(n)%k_nh4))
          mu2 = phyparam(n)%vmax * nh4 / (nh4 / (1.0-fa) + phyparam(n)%vmax/(fa*phyparam(n)%a_nh4))

          phot(n) = mu1 + mu2

          rnew(i,j,k,n) = mu1 / (mu1 + mu2)

          IF (phyparam(n)%k_dsi > 0.0) phot(n) = min(phot(n), phyparam(n)%vmax * dsi / (dsi / (1.0-fa) + phyparam(n)%vmax/(fa*phyparam(n)%a_dsi)))
          IF (phyparam(n)%k_dfe > 0.0) phot(n) = min(phot(n), phyparam(n)%vmax * dfe / (dfe / (1.0-fa) + phyparam(n)%vmax/(fa*phyparam(n)%a_dfe)))

          phot(n) = phot(n) * phyparam(n)%lf0 * (1.0 - exp(-phyparam(n)%alpha*lint(i,j,k)/phyparam(n)%pmax))&
                                              *        exp(-phyparam(n)%beta *lint(i,j,k)/phyparam(n)%pmax) &
                            * exp(phyparam(n)%k_phot * t) * phy(n)

          resp(n)  = phyparam(n)%resp * exp(phyparam(n)%k_resp * t) * phy(n)
          pmort(n) = phyparam(n)%mort * exp(phyparam(n)%k_mort * t) * (phy(n)**2)

          netpp(i,j,k,n) = (phot(n) - resp(n)) * r_cn
       END DO


       DO n=1, n_zoo
          zoo(n) = tracer(i,j,k,tracer_index_zoo(n))

          expkt = exp(zooparam(n)%k_gp * t)

          graz(:,n) = 0.0
          DO m=1, n_phy
             IF (zooparam(n)%graz(m) > 0.0 .AND. phy(m) > gp_threshold) &
                  graz(m, n) = zooparam(n)%graz(m)*(1.0 - exp( -zooparam(n)%lambda*(phy(m)-gp_threshold))) * expkt * zoo(n)
          END DO

          pred(:,n) = 0.0
          DO m=1, n_zoo
             IF (zooparam(n)%pred(m) > 0.0 .AND. zoo(m) > gp_threshold) &
                  pred(m, n) = zooparam(n)%pred(m)*(1.0 - exp( -zooparam(n)%lambda*(zoo(m)-gp_threshold))) * expkt * zoo(n)
          END DO

          zmort(n) = zooparam(n)%mort * exp(zooparam(n)%k_mort * t) * (zoo(n)**2)

          IF (n==3) THEN
          !preference factor for zoo3
             pred(1,n) = pred(1,n) * exp(-pref1 * zoo(2))
             graz(2,n) = graz(2,n) * exp(-pref2 * (zoo(1) + zoo(2)))
          END IF

          sumgp(n) = sum(graz(:,n)) + sum(pred(:,n))
       END DO

       DO n=1, n_phy
          tracer(i,j,k,tracer_index_phy(n)) = tracer(i,j,k,tracer_index_phy(n)) &
               + ((1.0 - phyparam(n)%gamma) * phot(n) - resp(n) - pmort(n) - sum(graz(n,:))) * dtime
       END DO

       DO n=1, n_zoo
          tracer(i,j,k,tracer_index_zoo(n)) = tracer(i,j,k,tracer_index_zoo(n)) &
               + (zooparam(n)%beta*sumgp(n) - zmort(n) - sum(pred(n,:))) * dtime
       END DO

       ! nitrification
       dno3 = v_nitrif * exp(k_nitrif * t) * nh4 * dtime
       dnh4 = - dno3

       ddon = 0.0
       dpon = 0.0

       DO n=1, n_phy
          ! photosynthesis - respiration
          dno3 = dno3 - (phot(n) - resp(n)) * rnew(i,j,k,n)       * dtime
          dnh4 = dnh4 - (phot(n) - resp(n)) * (1.0-rnew(i,j,k,n)) * dtime

          ! extracellular excertion
          ddon = ddon + phyparam(n)%gamma * phot(n) * dtime

          ! mortality of phytoplankton
          dpon = dpon + pmort(n) * dtime
       END DO

       DO n=1, n_zoo
          ! excretion
          dnh4 = dnh4 + (zooparam(n)%alpha - zooparam(n)%beta) * sumgp(n) * dtime

          ! egestion
          dpon = dpon + (1.0 - zooparam(n)%alpha) * sumgp(n) * dtime

          ! mortality of zooplankton
          dpon = dpon + zmort(n) * dtime
       END DO

       ! re-mineraization of DON
       tmp = v_donrem * exp(k_donrem*t) * don * dtime
       dnh4 = dnh4 + tmp
       ddon = ddon - tmp

       ! re-mineraization of PON
       tmp = v_ponrem * exp(k_ponrem*t) * pon * dtime
       dnh4 = dnh4 + tmp
       dpon = dpon - tmp

       ! PON decomposition to DON
       tmp = v_pondec * exp(k_pondec*t) * pon * dtime
       ddon = don + tmp
       dpon = pon - tmp

       tracer(i,j,k,tracer_index_no3) = tracer(i,j,k,tracer_index_no3) + dno3
       tracer(i,j,k,tracer_index_nh4) = tracer(i,j,k,tracer_index_nh4) + dnh4
       tracer(i,j,k,tracer_index_don) = tracer(i,j,k,tracer_index_don) + ddon
       tracer(i,j,k,tracer_index_pon) = tracer(i,j,k,tracer_index_pon) + dpon

       IF (dfe > dfe_threshold) THEN
          r_sin = r_sinh
       ELSE
          r_sin = r_sinl
       END IF

       ! OPAL disolution
       tmp = v_psidis * exp(k_psidis*t) * psi * dtime

       tracer(i,j,k,tracer_index_dsi) = tracer(i,j,k,tracer_index_dsi) + tmp &
            - ((1.0-phyparam(2)%gamma) * phot(2) - resp(2)) * r_sin * dtime

       tracer(i,j,k,tracer_index_psi) = tracer(i,j,k,tracer_index_psi) - tmp &
            + (pmort(2) + sum(graz(2,:))) * r_sin * dtime

       ddfe = (dno3 + dnh4)*r_fen + (fz_dustfe(i,j,k) - fz_dustfe(i,j,k-1))/dz_star(i,j,k) * dtime
       IF (k==surface_k(i,j)) ddfe = ddfe + a_dustfe*input(i,j,2) / dz_star(i,j,k) * dtime

       dpfe = 0.0

       ! Particulate Fe desorption
       tmp = l_fedes * exp(-a_fedes * (1.0D0/(273.15+t)) - 1.0D0/(303.15)) * pfe * dtime

       ! disolved Fe scavanging
       tmp = tmp - lambda_fescav * ((pon*abs(w_pon)) * r_cn )   * dfe * dtime
       IF (dfe > c_ligand) tmp = tmp - gamma_fescav * (dfe - c_ligand) * dfe * dtime

       tracer(i,j,k,tracer_index_dfe) = tracer(i,j,k,tracer_index_dfe) + tmp
       tracer(i,j,k,tracer_index_pfe) = tracer(i,j,k,tracer_index_pfe) - tmp
    END DO
    END DO
    END DO


    DO n=1, n_phy
       CALL update_tracer_boundary(tracer_index_phy(n))

       CALL checkout('NETPP' // trim(format(n, 'I0')), netpp(:,:,:,n))  ! net primary production for each phytoplankton category [molC / ( m^3 s)]
       CALL checkout('RNEW'  // trim(format(n, 'I0')), rnew( :,:,:,n))  ! f-ratio for each phytoplankton category
    END DO
    IF (require_checkout('NETPP')) THEN
       CALL checkout('NETPP', sum(netpp,4))! net primary production for total phytoplanktons [molC / ( m^3 s)]
    END IF

    DO n=1, n_zoo
       CALL update_tracer_boundary(tracer_index_zoo(n))
    END DO

    CALL update_tracer_boundary(tracer_index_no3)
    CALL update_tracer_boundary(tracer_index_nh4)
    CALL update_tracer_boundary(tracer_index_don)
    CALL update_tracer_boundary(tracer_index_pon)
    CALL update_tracer_boundary(tracer_index_dsi)
    CALL update_tracer_boundary(tracer_index_psi)
    CALL update_tracer_boundary(tracer_index_dfe)
    CALL update_tracer_boundary(tracer_index_pfe)

  END SUBROUTINE step_ecosystem

END MODULE ecosystem
