#include "macro.h"

MODULE parameters
  USE misc
  USE parallel
  IMPLICIT NONE
  SAVE

  LOGICAL :: hydrostatic  = .FALSE.
  LOGICAL :: rigid_lid    = .FALSE.
  LOGICAL :: fix_meanssh  = .FALSE.
  LOGICAL :: vcoord_zstar = .FALSE.
  LOGICAL :: iceshelf_melt= .FALSE.
  LOGICAL :: offline = .FALSE.

  LOGICAL :: use_landwater  = .FALSE.
  LOGICAL :: use_gls        = .FALSE.
  LOGICAL :: use_npzd       = .FALSE.
  LOGICAL :: use_ecosystem  = .FALSE.
  LOGICAL :: use_particles  = .FALSE.

  REAL(4) :: gravity      = 9.80665    ![m/s2]
  REAL(4) :: gravity_angle_deg = 0.0   !angle of the gravity vector in xz-plane, degree
  REAL(4) :: gravity_angle_tan = 0.0   !angle of the gravity vector in xz-plane, tangent

  REAL(4) :: earth_radius = 6378E3     ![m]

  REAL(4) :: rho_0        = 1025    ![kg/m3]
  REAL(4) :: rho_ice      =  916    ![kg/m3]
  REAL(4) :: rho_air      = 1.22    ![kg/m3]

  REAL(4) :: rho_sediment = 2650    ![kg/m3]
  REAL(4) :: rho_frazil   =  916    ![kg/m3]

  REAL(4) :: rho_bubble   =  1.3    ![kg/m3]

  !for microplastics
  REAL(4) :: rho_pla_pe   =  965    ![kg/m3] polyethylene
  REAL(4) :: rho_pla_pp   =  910    ![kg/m3] polypropylene
  REAL(4) :: rho_pla_pet  = 1390    ![kg/m3] polyethlen-terephthalate
  REAL(4) :: rho_pla_ps   = 1050    ![kg/m3] polystyrene
  REAL(4) :: rho_pla_pa   = 1140    ![kg/m3] polyamide
  REAL(4) :: rho_pla_pmma = 1190    ![kg/m3] polymethyl-methacrylate


  REAL(4) :: cp       = 3.974E3   ![J / (kg K)]
  REAL(4) :: l_freeze = 3.34E5    ![J / kg]
  REAL(4) :: l_vapor  = 2.26E6    ![J / kg]

  REAL(4) :: cp_air   = 1.005E3   ![J / (kg K)]

  REAL(4) :: s_ref    = 35.0      ![psu]
  REAL(4) :: s_ice    = 10.0      ![psu]
  REAL(4) :: s_frazil =  0.0      ![psu]

  REAL(4) :: visch  = 0.0         ![m2/s]
  REAL(4) :: viscv  = 0.0         ![m2/s]

  REAL(4) :: diffh  = 0.0         ![m2/s]
  REAL(4) :: diffv  = 0.0         ![m2/s]

  REAL(4) :: c_smagorinsky_visc = 0.0
  REAL(4) :: c_smagorinsky_diff = 0.0

  LOGICAL :: delta2_sml93 = .FALSE.
  LOGICAL :: delta2_cube  = .FALSE.

  REAL(4) :: eps_bhfilter = 0.03125

  REAL(4) :: userparam(8) = UNDEF

  INTEGER :: tracer_advection_scheme = 3    ! 0: no advection
                                            ! 1: O(1) upwind
                                            ! 2: O(3) upwind (QUICKESS) with flux limiter
                                            ! 3: COSMIC/QUICKEST

  INTEGER :: tracer_diffusion_scheme = 2    ! 0: no diffusion
                                            ! 1: isotropic
                                            ! 2: horizontal/vertical diffusion
                                            ! 3: isopycnal/diapycnal diffusion

  LOGICAL :: tracer_correct_div  = .TRUE.
  LOGICAL :: tracer_mask_sfcflux = .FALSE.

  REAL(4) :: nu_mol = 1.0E-6                ![m2/s]
  REAL(4) :: prandtl_mol = 7.0E0
  REAL(4) :: schmidt_mol = 7.0E2

  REAL(4) :: albedo_ocean = 0.066 ! ocean surface albedo
  REAL(4) :: albedo_ice   = 0.80  ! sea ice albedo

  REAL(4) :: emissivity_ocean = 0.95 ! ocean surface emissivity
  REAL(4) :: emissivity_ice   = 0.95 ! sea ice emissivity

  REAL(4) :: radsw_r = 0.58 ! ratio of shortwave absorbed in the top layer
  REAL(4) :: radsw_z = 0.35 ! scale hight of the shortwave absorption [m]

  REAL(4), PARAMETER :: gas_const   = 8.314     ![J / (mol K)] ! universal gas constant
  REAL(4), PARAMETER :: molmass_co2 = 44.01E-3  ![kg / mol]
  REAL(4), PARAMETER :: molvol_co2  = 33.00E-6  ![m3 / mol]
  REAL(4), PARAMETER :: atm_pa      = 101325    ![Pa]
  REAL(4), PARAMETER :: db_pa       = 1.0E4     ![Pa]

  REAL(4), PARAMETER :: molmass_fe  = 55.85E-3  ![kg / mol]

  REAL(4), PARAMETER :: day_sec     = 86400.0E0 ![s]

  REAL(4), PARAMETER :: sigma_sb = 5.67E-8      ![W K^-4  m^-2]  Stefen-Boltzmann constant

  REAL(4) :: a_ebcn = 0.5625 ! non-dimensional control parameter for the Euler-backward/Crank-Nicolson scheme
  REAL(8) :: b_ebcn, c_ebcn

  LOGICAL :: slip_surface =  .TRUE.
  LOGICAL :: slip_top     =  .FALSE.
  LOGICAL :: slip_bottom  =  .FALSE.
  LOGICAL :: slip_side    =  .TRUE.

  REAL(4) :: slipfactor_surface  ! 1.0 for slip, -1.0 for non-slip
  REAL(4) :: slipfactor_top
  REAL(4) :: slipfactor_bottom
  REAL(4) :: slipfactor_side

  REAL(4) :: landwater_hdry  =  0.01 ![m]
  REAL(4) :: landwater_limit = 100.0 ![m]

  LOGICAL :: wind_relative = .FALSE.

  INTEGER :: winddrag_scheme = 1    ! 0: const C_d, default 1.2E-3
                                    ! 1: Large and Pond (1981), JPO
                                    ! 2: Kara et al. (2000), J. Oce. Atm. Tech.
  REAL(4) :: drag_wind    = 1.2E-3
  REAL(4) :: drag_ice     = 5.5E-3
  REAL(4) :: drag_bottom  = 5.5E-3

  REAL(4) :: turn_wind    =  0.0    ! [deg] turning angle for wind stress
  REAL(4) :: turn_ice     = 25.0    ! [deg] turning angle for ocean-ice stress
  REAL(4) :: turn_bottom  =  0.0    ! [deg] turning angle for bottom stress

  REAL(8) :: cos_turn_wind
  REAL(8) :: sin_turn_wind

  REAL(8) :: cos_turn_ice
  REAL(8) :: sin_turn_ice

  REAL(8) :: cos_turn_bottom
  REAL(8) :: sin_turn_bottom

  REAL(8) :: cos_gravity
  REAL(8) :: sin_gravity

  REAL(4) :: k_ice     = 2.04 ![W/m/K]

  REAL(4) :: hice_min  = 0.05 ![m]

  REAL(4) :: karman_const = 0.41

  INTEGER :: eq_of_state = 2  ! 0 : constant rho(T,S,p)=rho_0
                              ! 1 : linear eq. of state (no dependency on depth)
                              ! 2 : EOS-80 (UNESCO Tech. report, 1980)
                              ! TEOS-10 is not yet supported.

  REAL(4) :: eos_t_ref      = 10.0
  REAL(4) :: eos_s_ref      = 35.0
  REAL(4) :: eos_lin_alpha  =  2.0E-4
  REAL(4) :: eos_lin_beta   =  8.0E-4

  REAL(4) :: mld_dsigma     = UNDEF ! the dnesity anomaly criteria for the mixd layer depth definition [kg/m^3]

  INTEGER :: viscosity_scheme = 2  ! 0 : skip viscosity term
                                   ! 1 : isotropic visocity
                                   ! 2 : horizontal/vertical viscosity
                                   ! 3 : isopycnal/diapycnal viscosity (currently not working)

  INTEGER :: nonlinear_scheme = 3  ! 0 : skip nonlinear term
                                   ! 1 : O(2)-central
                                   ! 2 : O(4)-central
                                   ! 3 : energy-conservative vector-invariant form

  INTEGER :: pgradient_scheme = 1  ! 1 : default
                                   ! 2 : momentum-conservative ***EXPERIMENTAL***

  INTEGER :: buoyancy_scheme  = 1  ! 0 : no buoyancy force
                                   ! 1 : default, calculate hydrstatic pressure gradient
                                   ! 2 : explicit vertical acceralation for anomaly from the constant reference density rho_0
                                   ! 3 : explicit vertical acceralation for anomaly from the horizontally-averaged density rho_bar(k)

  INTEGER :: les_scheme  = 0       ! 0 : no LES subgrid viscosity/diffusivity model
                                   ! 1 : simple Smagorinsky model (Smagorinsky 1963)
                                   ! 2 : wall-adoptive local eddy-viscosity (WALE) model (Nicoud and Ducros 1999)
                                   ! 3 : filetered structure function (FSF) model (Ducros, Comte and Lesieur 1996)

  REAL(4) :: prandtl_les      = 1.0
  REAL(4) :: c_smagorinsky_les= 0.15 ! 1/pi (3/2 C_K)^{-3/4} where the Kolmogorf constant C_K=1.4
  REAL(4) :: c_fsf  = 8.45E-4        ! 0.0014 C_K^{-3/2}
  REAL(4) :: c_wale = 0.5
  REAL(4) :: limit_les        = UNDEF

  REAL(4) :: p_hibler  =  5.0E3
  REAL(4) :: c_hibler  = 20.0
  REAL(4) :: e_hibler  =  2.0

  REAL(4) :: e_hunke   = 0.25


  INTEGER :: ice_split      = 30

  REAL(4) :: depth_offset = 0.0

  REAL(4) :: iceshelf_gamma_t = 1.0E-4    ! turbulent thermal exchenge coefficient for iceshelf melt [m/s]
  REAL(4) :: iceshelf_gamma_s = UNDEF     ! turbulent salt    exchenge coefficient for iceshelf melt [m/s]

  LOGICAL :: bulk_bottom_stress = .FALSE.

  LOGICAL :: particle_advection_rk4     = .TRUE. ! 4th-order Runge-Kutta scheme for particle advection, defalt is TRUE
  INTEGER :: particle_advection_interp  = 1      ! interpolation of velocity field perpedicular to each component
                                                 ! 2: Quadratic, 1: Linear, 0: W/O interporation, default quadratic
  REAL(4) :: particle_repulsion_surface = 1.0E-2 ! default repulsion coefficients for particles
  REAL(4) :: particle_repulsion_bottom  = 1.0E-2 !
  REAL(4) :: particle_repulsion_side    = 1.0E-2 !

  REAL(4) :: radi_nudge(2)  = (/1.0E-4, 1.0E-2/)  ! nudging timesteps on the radiative boundaries for outward and inward propagation, respectively)

  REAL(4) :: radi_nudge_e(2) = UNDEF
  REAL(4) :: radi_nudge_w(2) = UNDEF
  REAL(4) :: radi_nudge_n(2) = UNDEF
  REAL(4) :: radi_nudge_s(2) = UNDEF

  LOGICAL :: radi_tracers   = .TRUE.
  LOGICAL :: radi_tangential= .TRUE.
  LOGICAL :: radi_oblique   = .TRUE.

  LOGICAL :: w_topdown     = .FALSE.

  REAL(8) :: u_cycleoffset_x = 0.0
  REAL(8) :: u_cycleoffset_y = 0.0
  REAL(8) :: u_cycleoffset_z = 0.0
  REAL(8) :: v_cycleoffset_x = 0.0
  REAL(8) :: v_cycleoffset_y = 0.0
  REAL(8) :: v_cycleoffset_z = 0.0

  CONTAINS
    SUBROUTINE init_parameters
      INTEGER :: i, j, k

      INTEGER   :: iostat
      CHARACTER :: iomsg

      REAL(4) :: visc_h, visc_v, diff_h, diff_v       ! aliases for backward compatibility
      INTEGER :: advection_scheme, diffusion_scheme   ! 

      NAMELIST / parameters / &
           gravity,      &
           gravity_angle_deg, &
           gravity_angle_tan, &
           earth_radius, &
           rho_0,   &
           rho_ice, &
           rho_air, &
           rho_sediment, &
           rho_frazil,   &
           rho_bubble,   &
           visch, viscv, &
           diffh, diffv, &
           visc_h, visc_v, diff_h, diff_v, & !for backward-compatibility
           c_smagorinsky_visc, &
           c_smagorinsky_diff, &
           delta2_sml93,       &
           delta2_cube,        &
           eps_bhfilter, &
           tracer_advection_scheme, advection_scheme, &
           tracer_diffusion_scheme, diffusion_scheme, &
           tracer_correct_div,      &
           tracer_mask_sfcflux,     &
           nu_mol, &
           prandtl_mol, &
           schmidt_mol, &
           viscosity_scheme, &
           nonlinear_scheme, &
           pgradient_scheme, &
           buoyancy_scheme,  &
           les_scheme,       &
           prandtl_les,      &
           limit_les,        &
           c_smagorinsky_les,&
           wind_relative, &
           winddrag_scheme, &
           drag_wind,   &
           drag_ice,    &
           drag_bottom, &
           turn_wind,   &
           turn_ice,    &
           turn_bottom, &
           landwater_hdry,  &
           landwater_limit, &
           k_ice, &
           eq_of_state, &
           eos_t_ref, eos_s_ref, eos_lin_alpha, eos_lin_beta,  &
           mld_dsigma, &
           cp, &
           l_freeze, l_vapor, &
           s_ref, s_ice, s_frazil, &
           e_hibler, p_hibler, c_hibler, &
           e_hunke, &
           albedo_ocean, &
           albedo_ice,   &
           emissivity_ocean, &
           emissivity_ice,   &
           radsw_r, radsw_z, &
           ice_split, &
           depth_offset, &
           slip_surface, slip_top, slip_bottom, slip_side, &
           bulk_bottom_stress, &
           iceshelf_gamma_t,   &
           iceshelf_gamma_s,   &
           a_ebcn,    &
           radi_nudge,       &
           radi_nudge_e,     &
           radi_nudge_w,     &
           radi_nudge_n,     &
           radi_nudge_s,     &
           radi_tracers,     &
           radi_tangential,  &
           radi_oblique,     &
           w_topdown, &
           u_cycleoffset_x, &
           u_cycleoffset_y, &
           u_cycleoffset_z, &
           v_cycleoffset_x, &
           v_cycleoffset_y, &
           v_cycleoffset_z, &
           particle_advection_rk4,     &
           particle_advection_interp,  &
           particle_repulsion_surface, &
           particle_repulsion_bottom,  &
           particle_repulsion_side,    &
           userparam

      !aliases for backward compatibility
      visc_h = UNDEF
      visc_v = UNDEF
      diff_h = UNDEF
      diff_v = UNDEF
      advection_scheme = -1
      diffusion_scheme = -1

      IF (rank==0) THEN
         REWIND(CONFIG_UNIT)
         READ(CONFIG_UNIT, NML=parameters, IOSTAT=iostat, IOMSG=iomsg)
         CALL assert(iostat <= 0, "failed to read PARAMETERS namelist", iomsg)
         IF (advection_scheme /= -1) tracer_advection_scheme = advection_scheme
         IF (diffusion_scheme /= -1) tracer_diffusion_scheme = diffusion_scheme
         IF (visc_h /= UNDEF) visch = visc_h
         IF (visc_v /= UNDEF) viscv = visc_v
         IF (diff_h /= UNDEF) diffh = diff_h
         IF (diff_v /= UNDEF) diffv = diff_v
      END IF

      CALL bcast(gravity)
      CALL bcast(gravity_angle_deg)
      CALL bcast(gravity_angle_tan)

      CALL bcast(earth_radius)

      CALL bcast(rho_0)
      CALL bcast(rho_ice)
      CALL bcast(rho_air)
      CALL bcast(rho_sediment)
      CALL bcast(rho_frazil)
      CALL bcast(rho_bubble)

      CALL bcast(visch)
      CALL bcast(viscv)

      CALL bcast(diffh)
      CALL bcast(diffv)

      CALL bcast(c_smagorinsky_visc)
      CALL bcast(c_smagorinsky_diff)
      CALL bcast(delta2_sml93)
      CALL bcast(delta2_cube)

      CALL bcast(eps_bhfilter)

      CALL bcast(nu_mol)
      CALL bcast(prandtl_mol)
      CALL bcast(schmidt_mol)

      CALL bcast(viscosity_scheme)
      CALL bcast(nonlinear_scheme)
      CALL bcast(pgradient_scheme)
      CALL bcast(buoyancy_scheme)
      CALL bcast(les_scheme)

      CALL bcast(prandtl_les)
      CALL bcast(limit_les)
      CALL bcast(c_smagorinsky_les)
      CALL bcast(c_wale)
      CALL bcast(c_fsf)

      CALL bcast(landwater_hdry)
      CALL bcast(landwater_limit)

      CALL bcast(wind_relative)

      CALL bcast(winddrag_scheme)

      CALL bcast(drag_wind)
      CALL bcast(drag_ice)
      CALL bcast(drag_bottom)

      CALL bcast(turn_wind)
      CALL bcast(turn_ice)
      CALL bcast(turn_bottom)

      CALL bcast(k_ice)

      CALL bcast(eq_of_state)
      CALL bcast(eos_t_ref)
      CALL bcast(eos_s_ref)
      CALL bcast(eos_lin_alpha)
      CALL bcast(eos_lin_beta)

      CALL bcast(mld_dsigma)

      CALL bcast(cp)
      CALL bcast(l_freeze)
      CALL bcast(l_vapor)

      CALL bcast(s_ref)
      CALL bcast(s_ice)
      CALL bcast(s_frazil)

      CALL bcast(hice_min)

      CALL bcast(p_hibler)
      CALL bcast(c_hibler)
      CALL bcast(e_hibler)

      CALL bcast(e_hunke)

      CALL bcast(albedo_ocean)
      CALL bcast(albedo_ice)

      CALL bcast(emissivity_ocean)
      CALL bcast(emissivity_ice)

      CALL bcast(radsw_r)
      CALL bcast(radsw_z)

      CALL bcast(ice_split)

      CALL bcast(depth_offset)

      CALL bcast(slip_surface)
      CALL bcast(slip_top)
      CALL bcast(slip_bottom)
      CALL bcast(slip_side)
      CALL bcast(bulk_bottom_stress)

      CALL bcast(tracer_advection_scheme)
      CALL bcast(tracer_diffusion_scheme)
      CALL bcast(tracer_correct_div)
      CALL bcast(tracer_mask_sfcflux)

      CALL bcast(iceshelf_gamma_t)
      CALL bcast(iceshelf_gamma_s)

      IF (iceshelf_gamma_s == UNDEF) iceshelf_gamma_s = 5.05D-3 * iceshelf_gamma_t !following Holland and Jenkins (1999)

      CALL bcast(a_ebcn)
      b_ebcn = 1.0D0/a_ebcn
      c_ebcn = (1.0D0-a_ebcn)/a_ebcn

      CALL bcast(radi_nudge)
      CALL bcast(radi_nudge_e)
      CALL bcast(radi_nudge_w)
      CALL bcast(radi_nudge_n)
      CALL bcast(radi_nudge_s)

      DO i=1, size(radi_nudge)
         IF (radi_nudge_e(i) == UNDEF) radi_nudge_e(i) = radi_nudge(i)
         IF (radi_nudge_w(i) == UNDEF) radi_nudge_w(i) = radi_nudge(i)
         IF (radi_nudge_n(i) == UNDEF) radi_nudge_n(i) = radi_nudge(i)
         IF (radi_nudge_s(i) == UNDEF) radi_nudge_s(i) = radi_nudge(i)

         radi_nudge(i)   = min(max(radi_nudge(i),   0.0), 1.0)
         radi_nudge_e(i) = min(max(radi_nudge_e(i), 0.0), 1.0)
         radi_nudge_w(i) = min(max(radi_nudge_w(i), 0.0), 1.0)
         radi_nudge_n(i) = min(max(radi_nudge_n(i), 0.0), 1.0)
         radi_nudge_s(i) = min(max(radi_nudge_s(i), 0.0), 1.0)
      END DO

      CALL bcast(radi_tracers)
      CALL bcast(radi_tangential)
      CALL bcast(radi_oblique)

      CALL bcast(w_topdown)

      CALL bcast(u_cycleoffset_x)
      CALL bcast(u_cycleoffset_y)
      CALL bcast(u_cycleoffset_z)
      CALL bcast(v_cycleoffset_x)
      CALL bcast(v_cycleoffset_y)
      CALL bcast(v_cycleoffset_z)

      CALL bcast(particle_advection_rk4)
      CALL bcast(particle_advection_interp)
      CALL bcast(particle_repulsion_surface)
      CALL bcast(particle_repulsion_bottom)
      CALL bcast(particle_repulsion_side)

      CALL bcast(userparam)

      IF (slip_surface) THEN
         slipfactor_surface =  1.0
      ELSE
         slipfactor_surface = -1.0
      END IF

      IF (slip_top) THEN
         slipfactor_top =  1.0
      ELSE
         slipfactor_top = -1.0
      END IF

      IF (slip_bottom) THEN
         slipfactor_bottom =  1.0
      ELSE
         slipfactor_bottom = -1.0
      END IF

      IF (slip_side) THEN
         slipfactor_side =  1.0
      ELSE
         slipfactor_side = -1.0
      END IF

      cos_turn_wind   = cos(pi*turn_wind   / 180.0)
      sin_turn_wind   = sin(pi*turn_wind   / 180.0)

      cos_turn_ice    = cos(pi*turn_ice    / 180.0)
      sin_turn_ice    = sin(pi*turn_ice    / 180.0)

      cos_turn_bottom = cos(pi*turn_bottom / 180.0)
      sin_turn_bottom = sin(pi*turn_bottom / 180.0)

      IF (gravity_angle_tan /= 0.0) THEN
         cos_gravity = 1.0 / sqrt(1+gravity_angle_tan**2)
         sin_gravity = gravity_angle_tan * cos_gravity
      ELSE
         cos_gravity = cos(pi*gravity_angle_deg / 180.0)
         sin_gravity = sin(pi*gravity_angle_deg / 180.0)
      END IF

      IF (les_scheme /= 0) THEN
         delta2_cube = .TRUE.
      END IF


      IF (rank==0) THEN
         WRITE(REPORT_UNIT, '(A)', ADVANCE='NO') 'Equation of state: '
         SELECT CASE (eq_of_state)
         CASE (0)
            WRITE(REPORT_UNIT, '(A)') 'Constant'
         CASE (1)
            WRITE(REPORT_UNIT, '(A)') 'Linear'
         CASE (2)
            WRITE(REPORT_UNIT, '(A)') 'EOS80'
         CASE DEFAULT
            CALL assert(.FALSE., "unsupported EQ_OF_STATE")
         END SELECT

         WRITE(REPORT_UNIT, '(A)', ADVANCE='NO') 'Non-Linear term: '
         SELECT CASE (nonlinear_scheme)
         CASE (0)
            WRITE(REPORT_UNIT, '(A)') 'skiped'
         CASE (1)
            WRITE(REPORT_UNIT, '(A)') 'O(2)-central'
         CASE (2)
            WRITE(REPORT_UNIT, '(A)') 'O(4)-central'
         CASE (3)
            WRITE(REPORT_UNIT, '(A)') 'vector-invariant form'
         CASE DEFAULT
            CALL assert(.FALSE., "unsupported NONLINEAR_SCHEME")
         END SELECT

         WRITE(REPORT_UNIT, '(A)', ADVANCE='NO') 'Viscosity term: '
         SELECT CASE (viscosity_scheme)
         CASE (0)
            WRITE(REPORT_UNIT, '(A)') 'skiped'
         CASE (1)
            WRITE(REPORT_UNIT, '(A)') 'isotropic'
         CASE (2)
            WRITE(REPORT_UNIT, '(A)') 'horizontal/vertical'
         CASE (3)
            WRITE(REPORT_UNIT, '(A)') 'isopycnal/diapycnal'
         CASE (4)
            WRITE(REPORT_UNIT, '(A)') 'horizontal/vertical (implicit)'
         CASE DEFAULT
            CALL assert(.FALSE., "unsupported VISCOSITY_SCHEME")
         END SELECT

         WRITE(REPORT_UNIT, '(A)', ADVANCE='NO') 'Pressure-Gradient term: '
         SELECT CASE (pgradient_scheme)
         CASE (1)
            WRITE(REPORT_UNIT, '(A)') 'default'
         CASE (2)
            WRITE(REPORT_UNIT, '(A)') 'momentum-conservative ***EXPERIMENTAL***'
         CASE (3)
            WRITE(REPORT_UNIT, '(A)') 'solve the pressure as volume-integral  ***EXPERIMENTAL***'
         CASE DEFAULT
            CALL assert(.FALSE., "unsupported PGRADIENT_SCHEME")
         END SELECT

         WRITE(REPORT_UNIT, '(A)', ADVANCE='NO') 'Buoyancy term: '
         SELECT CASE (buoyancy_scheme)
         CASE (0)
            WRITE(REPORT_UNIT, '(A)') 'ignored'
         CASE (1)
            WRITE(REPORT_UNIT, '(A)') 'hydrostatic pressure gradient'
         CASE (2)
            WRITE(REPORT_UNIT, '(A)') 'explicit vertical acceralation for anomaly from the constant reference density rho_0'
         CASE (3)
            WRITE(REPORT_UNIT, '(A)') 'explicit vertical acceralation for anomaly from the horizontally-averaged density rho_bar(k)'
         CASE (4)
            WRITE(REPORT_UNIT, '(A)') 'explicit vertical acceralation for anomaly from the reference density field rho_ref(i,j,k)'
         CASE DEFAULT
            CALL assert(.FALSE., "unsupported BUOYANCY_SCHEME")
         END SELECT

         CALL assert(sin_gravity==0.0 .OR. (buoyancy_scheme==2 .OR. buoyancy_scheme==4), "setting GRAVITY_ANGLE requires BUOYANCY_SCHEME=2 or 4")

         IF (les_scheme /= 0) THEN
            CALL assert(c_smagorinsky_visc == 0.0 .AND. c_smagorinsky_diff == 0.0, &
                 "C_SMAGORINSKY_VISC and C_SMAGORINSKY_DIFF are not compatible with LES_SCHEME, use C_SMAGORINSKY_LES and PRANDTL_LES instead.")

            WRITE(REPORT_UNIT, '(A)', ADVANCE='NO') "LES sugrid-model is enabled (force isotropic viscosity/diffusivity): "

            SELECT CASE (les_scheme)
            CASE (1)
               WRITE(REPORT_UNIT, '(A)') "Smagorinsky model"
            CASE (2)
               WRITE(REPORT_UNIT, '(A)') "WALE (wall-adopted local eddy-viscosity) model"
            CASE (3)
               WRITE(REPORT_UNIT, '(A)') "FSF (filterd structured-function) model"
            CASE DEFAULT
               CALL assert(.FALSE., "unsupported LES_SCHEME")
            END SELECT
         END IF

      END IF

   END SUBROUTINE init_parameters

END MODULE parameters
