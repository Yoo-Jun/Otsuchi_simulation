#include "macro.h"

MODULE geometry
  USE misc
  USE parallel
  USE calendar
  USE parameters
#ifdef PROFILE
  USE profile
#endif
  IMPLICIT NONE

  INTEGER, SAVE :: n_timestep = 0

  INTEGER, SAVE :: dimx = 1
  INTEGER, SAVE :: dimy = 1
  INTEGER, SAVE :: dimz = 1

  INTEGER, SAVE :: isize = 1
  INTEGER, SAVE :: jsize = 1
  INTEGER, SAVE :: ksize = 1

  INTEGER, PARAMETER :: slv=2
  INTEGER, PARAMETER :: maxdim = 32768

  REAL(8), ALLOCATABLE :: ssh(:,:)
  REAL(8), ALLOCATABLE :: ssh_old(:,:)

  LOGICAL, ALLOCATABLE :: lwflag(:,:)
  LOGICAL, ALLOCATABLE :: lwdried(:,:)
  REAL(8), ALLOCATABLE :: landelev(:,:)

  INTEGER, ALLOCATABLE :: surface_k(:,:)
  LOGICAL, ALLOCATABLE :: surface_flag(:,:)

  INTEGER, ALLOCATABLE :: bottom_k(:,:)
  LOGICAL, ALLOCATABLE :: bottom_flag(:,:)

  INTEGER(1), ALLOCATABLE :: imask2d(:,:)
  INTEGER(1), ALLOCATABLE :: imask3d(:,:,:)
  LOGICAL(1), ALLOCATABLE :: lmask2d(:,:)
  LOGICAL(1), ALLOCATABLE :: lmask3d(:,:,:)

  INTEGER, ALLOCATABLE :: maskindex(:)

  REAL(8), ALLOCATABLE :: cext(:,:)
  REAL(8), ALLOCATABLE :: cint(:,:)

  REAL(8), ALLOCATABLE :: dx(:,:)
  REAL(8), ALLOCATABLE :: dy(:,:)
  REAL(8), ALLOCATABLE :: dz(:)
  REAL(8), ALLOCATABLE :: dz_ref(:,:,:)
  REAL(8), ALLOCATABLE :: dz_star(:,:,:)
  INTEGER, ALLOCATABLE :: dzindex(:)

  REAL(8) :: dz_max, dz_min

  REAL(8), ALLOCATABLE :: cz_star(:,:,:)

  REAL(8), ALLOCATABLE :: dx1(:,:)
  REAL(8), ALLOCATABLE :: dy1(:,:)

  REAL(8), ALLOCATABLE :: idx0(:,:)
  REAL(8), ALLOCATABLE :: idy0(:,:)
  REAL(8), ALLOCATABLE :: idz0(:)

  REAL(8), ALLOCATABLE :: idx1(:,:)
  REAL(8), ALLOCATABLE :: idy1(:,:)
  REAL(8), ALLOCATABLE :: idz1(:)

  REAL(8), ALLOCATABLE :: idx2(:,:)
  REAL(8), ALLOCATABLE :: idy2(:,:)
  REAL(8), ALLOCATABLE :: idz2(:)

  REAL(8), ALLOCATABLE :: dvol(:,:,:)
  REAL(8), ALLOCATABLE :: dvol_ref(:,:,:)
  REAL(8), ALLOCATABLE :: dvol_old(:,:,:)
  REAL(8), SAVE :: total_vol
  REAL(8), SAVE :: total_area
  REAL(8), ALLOCATABLE :: layer_area(:)

  INTEGER(8), SAVE :: n_cells

  REAL(8), ALLOCATABLE :: dsx(:,:,:)
  REAL(8), ALLOCATABLE :: dsy(:,:,:)
  REAL(8), ALLOCATABLE :: dsx_ref(:,:,:)
  REAL(8), ALLOCATABLE :: dsy_ref(:,:,:)
  REAL(8), ALLOCATABLE :: dsx_old(:,:,:)
  REAL(8), ALLOCATABLE :: dsy_old(:,:,:)

  REAL(8), ALLOCATABLE :: dsz(:,:)

  REAL(8), ALLOCATABLE :: dsx2d(:,:)
  REAL(8), ALLOCATABLE :: dsy2d(:,:)

  REAL(8), ALLOCATABLE :: metxy(:,:)
  REAL(8), ALLOCATABLE :: metyx(:,:)

  REAL(8), ALLOCATABLE :: corx(:,:)
  REAL(8), ALLOCATABLE :: cory(:,:)
  REAL(8), ALLOCATABLE :: corz(:,:)

  REAL(8), ALLOCATABLE :: depth(:)

  REAL(8), ALLOCATABLE :: h_bathymetry(:,:)
  REAL(8), ALLOCATABLE :: h_iceshelf(:,:)

  REAL(8), ALLOCATABLE :: wct_ref(:,:)

  LOGICAL, ALLOCATABLE :: isflag(:,:,:)

  REAL(8), ALLOCATABLE :: delta2(:,:,:)

  REAL(8), ALLOCATABLE :: da_solver(:,:,:)

  LOGICAL, SAVE :: cycle_x = .FALSE.
  LOGICAL, SAVE :: cycle_y = .FALSE.
  LOGICAL, SAVE :: cycle_z = .FALSE.

  LOGICAL, SAVE :: open_e = .FALSE.
  LOGICAL, SAVE :: open_w = .FALSE.
  LOGICAL, SAVE :: open_n = .FALSE.
  LOGICAL, SAVE :: open_s = .FALSE.

  LOGICAL, SAVE :: tripolar = .FALSE.

  REAL(8), SAVE :: open_alpha_e = UNDEF
  REAL(8), SAVE :: open_alpha_w = UNDEF
  REAL(8), SAVE :: open_alpha_n = UNDEF
  REAL(8), SAVE :: open_alpha_s = UNDEF

  LOGICAL, SAVE :: radi_e = .FALSE.
  LOGICAL, SAVE :: radi_w = .FALSE.
  LOGICAL, SAVE :: radi_n = .FALSE.
  LOGICAL, SAVE :: radi_s = .FALSE.

  INTEGER, SAVE :: cfl_report_interval      = 30
  INTEGER, SAVE :: ssh_report_interval      =  0
  INTEGER, SAVE :: open_report_interval     =  0
  INTEGER, SAVE :: momentum_report_interval =  0
  INTEGER, SAVE :: energy_report_interval   =  0
  INTEGER, SAVE :: tracer_report_interval   = 30
  INTEGER, SAVE :: stress_report_interval   =  0

  LOGICAL, SAVE :: perfect_restart = .FALSE.

  LOGICAL, SAVE :: partial_step    = .TRUE.
  REAL(4), SAVE :: partial_step_threshold = 0.125
  LOGICAL, SAVE :: fill_chokepoint = .FALSE.

  LOGICAL, SAVE :: iceshelf_coupled  = .FALSE.
  LOGICAL, SAVE :: seaice_coupled    = .FALSE.

#ifdef PARALLEL_MPI
  REAL(4), ALLOCATABLE, PRIVATE :: sendbuf_r4_e(:)
  REAL(4), ALLOCATABLE, PRIVATE :: sendbuf_r4_w(:)
  REAL(4), ALLOCATABLE, PRIVATE :: sendbuf_r4_n(:)
  REAL(4), ALLOCATABLE, PRIVATE :: sendbuf_r4_s(:)
  REAL(4), ALLOCATABLE, PRIVATE :: sendbuf_r4_ne(:)
  REAL(4), ALLOCATABLE, PRIVATE :: sendbuf_r4_nw(:)
  REAL(4), ALLOCATABLE, PRIVATE :: sendbuf_r4_se(:)
  REAL(4), ALLOCATABLE, PRIVATE :: sendbuf_r4_sw(:)

  REAL(8), ALLOCATABLE, PRIVATE :: sendbuf_r8_e(:)
  REAL(8), ALLOCATABLE, PRIVATE :: sendbuf_r8_w(:)
  REAL(8), ALLOCATABLE, PRIVATE :: sendbuf_r8_n(:)
  REAL(8), ALLOCATABLE, PRIVATE :: sendbuf_r8_s(:)
  REAL(8), ALLOCATABLE, PRIVATE :: sendbuf_r8_ne(:)
  REAL(8), ALLOCATABLE, PRIVATE :: sendbuf_r8_nw(:)
  REAL(8), ALLOCATABLE, PRIVATE :: sendbuf_r8_se(:)
  REAL(8), ALLOCATABLE, PRIVATE :: sendbuf_r8_sw(:)

  REAL(4), ALLOCATABLE, PRIVATE :: recvbuf_r4_e(:)
  REAL(4), ALLOCATABLE, PRIVATE :: recvbuf_r4_w(:)
  REAL(4), ALLOCATABLE, PRIVATE :: recvbuf_r4_n(:)
  REAL(4), ALLOCATABLE, PRIVATE :: recvbuf_r4_s(:)
  REAL(4), ALLOCATABLE, PRIVATE :: recvbuf_r4_ne(:)
  REAL(4), ALLOCATABLE, PRIVATE :: recvbuf_r4_nw(:)
  REAL(4), ALLOCATABLE, PRIVATE :: recvbuf_r4_se(:)
  REAL(4), ALLOCATABLE, PRIVATE :: recvbuf_r4_sw(:)

  REAL(8), ALLOCATABLE, PRIVATE :: recvbuf_r8_e(:)
  REAL(8), ALLOCATABLE, PRIVATE :: recvbuf_r8_w(:)
  REAL(8), ALLOCATABLE, PRIVATE :: recvbuf_r8_n(:)
  REAL(8), ALLOCATABLE, PRIVATE :: recvbuf_r8_s(:)
  REAL(8), ALLOCATABLE, PRIVATE :: recvbuf_r8_ne(:)
  REAL(8), ALLOCATABLE, PRIVATE :: recvbuf_r8_nw(:)
  REAL(8), ALLOCATABLE, PRIVATE :: recvbuf_r8_se(:)
  REAL(8), ALLOCATABLE, PRIVATE :: recvbuf_r8_sw(:)

#ifdef PARALLEL3D
  REAL(4), ALLOCATABLE, PRIVATE :: sendbuf_r4_u(:)
  REAL(4), ALLOCATABLE, PRIVATE :: sendbuf_r4_l(:)
  REAL(4), ALLOCATABLE, PRIVATE :: sendbuf_r4_ue(:)
  REAL(4), ALLOCATABLE, PRIVATE :: sendbuf_r4_uw(:)
  REAL(4), ALLOCATABLE, PRIVATE :: sendbuf_r4_un(:)
  REAL(4), ALLOCATABLE, PRIVATE :: sendbuf_r4_us(:)
  REAL(4), ALLOCATABLE, PRIVATE :: sendbuf_r4_le(:)
  REAL(4), ALLOCATABLE, PRIVATE :: sendbuf_r4_lw(:)
  REAL(4), ALLOCATABLE, PRIVATE :: sendbuf_r4_ln(:)
  REAL(4), ALLOCATABLE, PRIVATE :: sendbuf_r4_ls(:)
  REAL(4), ALLOCATABLE, PRIVATE :: sendbuf_r4_une(:)
  REAL(4), ALLOCATABLE, PRIVATE :: sendbuf_r4_unw(:)
  REAL(4), ALLOCATABLE, PRIVATE :: sendbuf_r4_use(:)
  REAL(4), ALLOCATABLE, PRIVATE :: sendbuf_r4_usw(:)
  REAL(4), ALLOCATABLE, PRIVATE :: sendbuf_r4_lne(:)
  REAL(4), ALLOCATABLE, PRIVATE :: sendbuf_r4_lnw(:)
  REAL(4), ALLOCATABLE, PRIVATE :: sendbuf_r4_lse(:)
  REAL(4), ALLOCATABLE, PRIVATE :: sendbuf_r4_lsw(:)

  REAL(8), ALLOCATABLE, PRIVATE :: sendbuf_r8_u(:)
  REAL(8), ALLOCATABLE, PRIVATE :: sendbuf_r8_l(:)
  REAL(8), ALLOCATABLE, PRIVATE :: sendbuf_r8_ue(:)
  REAL(8), ALLOCATABLE, PRIVATE :: sendbuf_r8_uw(:)
  REAL(8), ALLOCATABLE, PRIVATE :: sendbuf_r8_un(:)
  REAL(8), ALLOCATABLE, PRIVATE :: sendbuf_r8_us(:)
  REAL(8), ALLOCATABLE, PRIVATE :: sendbuf_r8_le(:)
  REAL(8), ALLOCATABLE, PRIVATE :: sendbuf_r8_lw(:)
  REAL(8), ALLOCATABLE, PRIVATE :: sendbuf_r8_ln(:)
  REAL(8), ALLOCATABLE, PRIVATE :: sendbuf_r8_ls(:)
  REAL(8), ALLOCATABLE, PRIVATE :: sendbuf_r8_une(:)
  REAL(8), ALLOCATABLE, PRIVATE :: sendbuf_r8_unw(:)
  REAL(8), ALLOCATABLE, PRIVATE :: sendbuf_r8_use(:)
  REAL(8), ALLOCATABLE, PRIVATE :: sendbuf_r8_usw(:)
  REAL(8), ALLOCATABLE, PRIVATE :: sendbuf_r8_lne(:)
  REAL(8), ALLOCATABLE, PRIVATE :: sendbuf_r8_lnw(:)
  REAL(8), ALLOCATABLE, PRIVATE :: sendbuf_r8_lse(:)
  REAL(8), ALLOCATABLE, PRIVATE :: sendbuf_r8_lsw(:)

  REAL(4), ALLOCATABLE, PRIVATE :: recvbuf_r4_u(:)
  REAL(4), ALLOCATABLE, PRIVATE :: recvbuf_r4_l(:)
  REAL(4), ALLOCATABLE, PRIVATE :: recvbuf_r4_ue(:)
  REAL(4), ALLOCATABLE, PRIVATE :: recvbuf_r4_uw(:)
  REAL(4), ALLOCATABLE, PRIVATE :: recvbuf_r4_un(:)
  REAL(4), ALLOCATABLE, PRIVATE :: recvbuf_r4_us(:)
  REAL(4), ALLOCATABLE, PRIVATE :: recvbuf_r4_le(:)
  REAL(4), ALLOCATABLE, PRIVATE :: recvbuf_r4_lw(:)
  REAL(4), ALLOCATABLE, PRIVATE :: recvbuf_r4_ln(:)
  REAL(4), ALLOCATABLE, PRIVATE :: recvbuf_r4_ls(:)
  REAL(4), ALLOCATABLE, PRIVATE :: recvbuf_r4_une(:)
  REAL(4), ALLOCATABLE, PRIVATE :: recvbuf_r4_unw(:)
  REAL(4), ALLOCATABLE, PRIVATE :: recvbuf_r4_use(:)
  REAL(4), ALLOCATABLE, PRIVATE :: recvbuf_r4_usw(:)
  REAL(4), ALLOCATABLE, PRIVATE :: recvbuf_r4_lne(:)
  REAL(4), ALLOCATABLE, PRIVATE :: recvbuf_r4_lnw(:)
  REAL(4), ALLOCATABLE, PRIVATE :: recvbuf_r4_lse(:)
  REAL(4), ALLOCATABLE, PRIVATE :: recvbuf_r4_lsw(:)

  REAL(8), ALLOCATABLE, PRIVATE :: recvbuf_r8_u(:)
  REAL(8), ALLOCATABLE, PRIVATE :: recvbuf_r8_l(:)
  REAL(8), ALLOCATABLE, PRIVATE :: recvbuf_r8_ue(:)
  REAL(8), ALLOCATABLE, PRIVATE :: recvbuf_r8_uw(:)
  REAL(8), ALLOCATABLE, PRIVATE :: recvbuf_r8_un(:)
  REAL(8), ALLOCATABLE, PRIVATE :: recvbuf_r8_us(:)
  REAL(8), ALLOCATABLE, PRIVATE :: recvbuf_r8_le(:)
  REAL(8), ALLOCATABLE, PRIVATE :: recvbuf_r8_lw(:)
  REAL(8), ALLOCATABLE, PRIVATE :: recvbuf_r8_ln(:)
  REAL(8), ALLOCATABLE, PRIVATE :: recvbuf_r8_ls(:)
  REAL(8), ALLOCATABLE, PRIVATE :: recvbuf_r8_une(:)
  REAL(8), ALLOCATABLE, PRIVATE :: recvbuf_r8_unw(:)
  REAL(8), ALLOCATABLE, PRIVATE :: recvbuf_r8_use(:)
  REAL(8), ALLOCATABLE, PRIVATE :: recvbuf_r8_usw(:)
  REAL(8), ALLOCATABLE, PRIVATE :: recvbuf_r8_lne(:)
  REAL(8), ALLOCATABLE, PRIVATE :: recvbuf_r8_lnw(:)
  REAL(8), ALLOCATABLE, PRIVATE :: recvbuf_r8_lse(:)
  REAL(8), ALLOCATABLE, PRIVATE :: recvbuf_r8_lsw(:)
#endif

  INTEGER, PRIVATE :: mpireq_2dx_r4(4,0:5)
  INTEGER, PRIVATE :: mpireq_2dx_r8(4,0:5)
  INTEGER, PRIVATE :: mpireq_2dy_r4(4,0:5)
  INTEGER, PRIVATE :: mpireq_2dy_r8(4,0:5)
  INTEGER, PRIVATE :: mpireq_3dx_r4(4,0:5)
  INTEGER, PRIVATE :: mpireq_3dx_r8(4,0:5)
  INTEGER, PRIVATE :: mpireq_3dy_r4(4,0:5)
  INTEGER, PRIVATE :: mpireq_3dy_r8(4,0:5)
  INTEGER, PRIVATE :: mpireq_3dz_r4(4,0:5)
  INTEGER, PRIVATE :: mpireq_3dz_r8(4,0:5)

  CHARACTER(16) :: input_method
  CHARACTER(16) :: output_method
  CHARACTER(16) :: sendrecv_method

#ifdef MPIIO
  TYPE fileview_struct
     LOGICAL :: defined = .FALSE.
     INTEGER :: subarray_2d_i1
     INTEGER :: subarray_2d_r4
     INTEGER :: subarray_2d_r8
     INTEGER :: subarray_3d_i1
     INTEGER :: subarray_3d_r4
     INTEGER :: subarray_3d_r8
     INTEGER :: subarray_3d_i1_desc
     INTEGER :: subarray_3d_r4_desc
     INTEGER :: subarray_3d_r8_desc
     INTEGER(MPI_OFFSET_KIND) :: offset
  END TYPE fileview_struct

  INTEGER, PARAMETER :: MAX_FILEVIEW = 16
  TYPE(fileview_struct), SAVE, PRIVATE :: views(0:MAX_FILEVIEW)

  TYPE file_struct
     INTEGER(1), POINTER :: buffer(:)
     INTEGER :: handle
  END TYPE file_struct

  TYPE(file_struct), SAVE, PRIVATE :: file
#ifdef MPIIO_ASYNCHRONOUS
  INTEGER, PRIVATE, PARAMETER :: n_async = 16
  TYPE(file_struct), SAVE, PRIVATE :: file_async(1:n_async)
#endif

  LOGICAL, PRIVATE :: mpiio_zerofill = .TRUE.
#endif
#endif

  INTEGER, SAVE :: kstr_t(0:max_threads-1)
  INTEGER, SAVE :: kend_t(0:max_threads-1)
  INTEGER, SAVE :: ksize_t(0:max_threads-1)

  LOGICAL, SAVE :: bypass_solver = .FALSE.  ! for debug use only.
  LOGICAL, SAVE :: report_qsum   = .FALSE.  ! for debug use only.

  TYPE subregion_struct
     LOGICAL :: defined
     INTEGER :: x_start
     INTEGER :: y_start
     INTEGER :: z_start
     INTEGER :: x_end
     INTEGER :: y_end
     INTEGER :: z_end
     INTEGER :: x_size
     INTEGER :: y_size
     INTEGER :: z_size
  END type subregion_struct

  INTEGER, PARAMETER :: MAX_SUBREGION = 64
  TYPE(subregion_struct), SAVE :: regions(0:MAX_SUBREGION)

  INTERFACE update_boundary
     MODULE PROCEDURE update_boundary_2d_r4
     MODULE PROCEDURE update_boundary_2d_r8
     MODULE PROCEDURE update_boundary_2d_logical
     MODULE PROCEDURE update_boundary_3d_r4
     MODULE PROCEDURE update_boundary_3d_r8
     MODULE PROCEDURE update_boundary_3d_logical
  END INTERFACE

  INTERFACE update_boundary_2d
     MODULE PROCEDURE update_boundary_2d_r4
     MODULE PROCEDURE update_boundary_2d_r8
     MODULE PROCEDURE update_boundary_2d_logical
  END INTERFACE

  INTERFACE update_boundary_3d
     MODULE PROCEDURE update_boundary_3d_r4
     MODULE PROCEDURE update_boundary_3d_r8
     MODULE PROCEDURE update_boundary_3d_logical
  END INTERFACE

  INTERFACE update_boundary_xz
     MODULE PROCEDURE update_boundary_xz_r4
     MODULE PROCEDURE update_boundary_xz_r8
     MODULE PROCEDURE update_boundary_xz_logical
  END INTERFACE

  INTERFACE update_boundary_yz
     MODULE PROCEDURE update_boundary_yz_r4
     MODULE PROCEDURE update_boundary_yz_r8
     MODULE PROCEDURE update_boundary_yz_logical
  END INTERFACE

  INTERFACE update_boundary_x
     MODULE PROCEDURE update_boundary_x_r4
     MODULE PROCEDURE update_boundary_x_r8
  END INTERFACE

  INTERFACE update_boundary_y
     MODULE PROCEDURE update_boundary_y_r4
     MODULE PROCEDURE update_boundary_y_r8
  END INTERFACE

  INTERFACE update_boundary_z
     MODULE PROCEDURE update_boundary_z_r4
     MODULE PROCEDURE update_boundary_z_r8
  END INTERFACE

  INTERFACE read_data
     MODULE PROCEDURE read_data_2d_r4
     MODULE PROCEDURE read_data_2d_r8
     MODULE PROCEDURE read_data_3d_r4
     MODULE PROCEDURE read_data_3d_r8
  END INTERFACE

  INTERFACE read_data_2d
     MODULE PROCEDURE read_data_2d_r4
     MODULE PROCEDURE read_data_2d_r8
  END INTERFACE

  INTERFACE read_data_3d
     MODULE PROCEDURE read_data_3d_r4
     MODULE PROCEDURE read_data_3d_r8
  END INTERFACE

  INTERFACE read_data_xz
     MODULE PROCEDURE read_data_xz_r4
     MODULE PROCEDURE read_data_xz_r8
  END INTERFACE

  INTERFACE read_data_yz
     MODULE PROCEDURE read_data_yz_r4
     MODULE PROCEDURE read_data_yz_r8
  END INTERFACE

  INTERFACE read_data_x
     MODULE PROCEDURE read_data_x_r4
     MODULE PROCEDURE read_data_x_r8
  END INTERFACE

  INTERFACE read_data_y
     MODULE PROCEDURE read_data_y_r4
     MODULE PROCEDURE read_data_y_r8
  END INTERFACE

  INTERFACE read_data_z
     MODULE PROCEDURE read_data_z_r4
     MODULE PROCEDURE read_data_z_r8
  END INTERFACE

  INTERFACE write_data
     MODULE PROCEDURE write_data_2d_r4
     MODULE PROCEDURE write_data_2d_r8
     MODULE PROCEDURE write_data_3d_r4
     MODULE PROCEDURE write_data_3d_r8
  END INTERFACE

  INTERFACE write_data_2d
     MODULE PROCEDURE write_data_2d_r4
     MODULE PROCEDURE write_data_2d_r8
  END INTERFACE

  INTERFACE write_data_3d
     MODULE PROCEDURE write_data_3d_r4
     MODULE PROCEDURE write_data_3d_r8
  END INTERFACE

  INTERFACE write_data_xz
     MODULE PROCEDURE write_data_xz_r4
     MODULE PROCEDURE write_data_xz_r8
  END INTERFACE

  INTERFACE write_data_yz
     MODULE PROCEDURE write_data_yz_r4
     MODULE PROCEDURE write_data_yz_r8
  END INTERFACE

  INTERFACE write_data_x
     MODULE PROCEDURE write_data_x_r4
     MODULE PROCEDURE write_data_x_r8
  END INTERFACE write_data_x

  INTERFACE write_data_y
     MODULE PROCEDURE write_data_y_r4
     MODULE PROCEDURE write_data_y_r8
  END INTERFACE write_data_y

  INTERFACE write_data_z
     MODULE PROCEDURE write_data_z_r4
     MODULE PROCEDURE write_data_z_r8
  END INTERFACE write_data_z

  INTERFACE in_region_3d
     MODULE PROCEDURE in_region_3d_bycode
     MODULE PROCEDURE in_region_3d_bytype
  END INTERFACE in_region_3d

  INTERFACE in_region_2d
     MODULE PROCEDURE in_region_2d_bycode
     MODULE PROCEDURE in_region_2d_bytype
  END INTERFACE in_region_2d

  INTERFACE in_region_xz
     MODULE PROCEDURE in_region_xz_bycode
     MODULE PROCEDURE in_region_xz_bytype
  END INTERFACE in_region_xz

  INTERFACE in_region_yz
     MODULE PROCEDURE in_region_yz_bycode
     MODULE PROCEDURE in_region_yz_bytype
  END INTERFACE in_region_yz

  INTERFACE in_region_x
     MODULE PROCEDURE in_region_x_bycode
     MODULE PROCEDURE in_region_x_bytype
  END INTERFACE in_region_x

  INTERFACE in_region_y
     MODULE PROCEDURE in_region_y_bycode
     MODULE PROCEDURE in_region_y_bytype
  END INTERFACE in_region_y

  INTERFACE in_region_z
     MODULE PROCEDURE in_region_z_bycode
     MODULE PROCEDURE in_region_z_bytype
  END INTERFACE in_region_z

  LOGICAL, SAVE :: assert_nan = .FALSE.

CONTAINS
  SUBROUTINE init_geometry(default_inputdir, default_outputdir)
    CHARACTER(*), INTENT(IN), OPTIONAL :: default_inputdir
    CHARACTER(*), INTENT(IN), OPTIONAL :: default_outputdir

    CHARACTER(512) :: inputdir
    CHARACTER(512) :: outputdir

    REAL(8) :: delta_x(maxdim)
    REAL(8) :: delta_y(maxdim)
    REAL(8) :: delta_z(maxdim)
    REAL(8) :: tmp(maxdim)

    CHARACTER(128) :: metric_x
    CHARACTER(128) :: metric_y
    CHARACTER(128) :: metric_z
    CHARACTER(128) :: metric_xy
    CHARACTER(128) :: metric_yx

    CHARACTER(128) :: coriolis_x
    CHARACTER(128) :: coriolis_y
    CHARACTER(128) :: coriolis_z

    CHARACTER(128) :: bathymetry
    CHARACTER(128) :: iceshelf
    CHARACTER(128) :: topomask
    CHARACTER(128) :: topography ! for backward compatibility

    LOGICAL :: gridwise_bathymetry
    LOGICAL :: gridwise_iceshelf

    LOGICAL :: output_dx
    LOGICAL :: output_dy
    LOGICAL :: output_dsx
    LOGICAL :: output_dsy
    LOGICAL :: output_dsz
    LOGICAL :: output_dvol
    LOGICAL :: output_mask
    LOGICAL :: output_mask_r8

    CHARACTER(12) :: precision
    CHARACTER(20) :: coordinate

    REAL(8) :: topomin

    INTEGER :: i, j, k
    INTEGER :: n

    REAL(8) :: d1, d2, d3
    REAL(8) :: lat

    INTEGER :: io_kind
    INTEGER :: ierr

    INTEGER        :: iostat
    CHARACTER(256) :: iomsg

    REAL(4) :: origin_x, origin_y, origin_z !defined for Python module, igonred in the model run
    INTEGER :: origin_i, origin_j, origin_k

    NAMELIST / geometry / &
         dimx, dimy, dimz, &
         cycle_x, cycle_y, cycle_z, &
         open_e, open_w, open_n, open_s, &
         open_alpha_e, &
         open_alpha_w, &
         open_alpha_n, &
         open_alpha_s, &
         tripolar,     &
         delta_x, delta_y, delta_z, &
         metric_x,   &
         metric_y,   &
         metric_z,   &
         metric_xy,  &
         metric_yx,  &
         coriolis_x, &
         coriolis_y, &
         coriolis_z, &
         bathymetry, &
         iceshelf,   &
         topomask,   &
         topography, & ! for backward compatibility
         gridwise_bathymetry, &
         gridwise_iceshelf,   &
         partial_step,    &
         partial_step_threshold, &
         fill_chokepoint, &
         output_dx,   &
         output_dy,   &
         output_dsx,  &
         output_dsy,  &
         output_dsz,  &
         output_dvol, &
         output_mask, &
         output_mask_r8,&
         inputdir,    &
         outputdir,   &
         precision,   &
         origin_x, origin_y, origin_z, &
         origin_i, origin_j, origin_k, &
         coordinate

    bathymetry  = ''
    iceshelf    = ''
    topomask    = ''
    topography  = ''

    gridwise_bathymetry = .FALSE.
    gridwise_iceshelf   = .FALSE.

    metric_x    = ''
    metric_y    = ''
    metric_z    = ''
    metric_xy   = ''
    metric_yx   = ''

    coriolis_x  = ''
    coriolis_y  = ''
    coriolis_z  = ''

    inputdir    = ''
    outputdir   = ''

    precision  = 'REAL8'

    coordinate = ''

    output_dx   = .FALSE.
    output_dy   = .FALSE.
    output_dsx  = .FALSE.
    output_dsy  = .FALSE.
    output_dsz  = .FALSE.
    output_dvol = .FALSE.
    output_mask = .TRUE.
    output_mask_r8 = .FALSE.

    delta_x(:) = 0.0
    delta_y(:) = 0.0
    delta_z(:) = 0.0

    IF (rank==0) THEN

       REWIND(CONFIG_UNIT)
       READ(CONFIG_UNIT, NML=geometry, IOSTAT=iostat, IOMSG=iomsg)

       CALL assert(iostat >= 0, "GEOMETRY namelist is mandatory")
       CALL assert(iostat == 0, "failed to read GEOMETRY namelist", iomsg)

       CALL assert(dimx <= maxdim, "DIMX should be less than MAXDIM")
       CALL assert(dimy <= maxdim, "DIMY should be less than MAXDIM")
       CALL assert(dimz <= maxdim, "DIMZ should be less than MAXDIM")

       CALL assert(dimz > 1 .OR. hydrostatic .OR. offline, "non-hydrostatic mode cannot run with DIMZ=1, use HYDROSTATIC mode")

       IF (delta_x(1)==0.0) delta_x(1) = 1.0
       DO i=2, dimx
          IF (delta_x(i) == 0.0) delta_x(i) = delta_x(i-1)
       END DO

       IF (delta_y(1)==0.0) delta_y(1) = 1.0
       DO j=2, dimy
          IF (delta_y(j) == 0.0) delta_y(j) = delta_y(j-1)
       END DO

       IF (delta_z(1)==0.0) delta_z(1) = 1.0
       DO k=2, dimz
          IF (delta_z(k) == 0.0) delta_z(k) = delta_z(k-1)
       END DO

       IF (cycle_x) CALL assert(.NOT. (open_e .OR. open_w), "CYCLE_X and OPEN_E/W cannot be used simultaneously")
       IF (cycle_y) CALL assert(.NOT. (open_n .OR. open_s), "CYCLE_Y and OPEN_N/S cannot be used simultaneously")

       IF (tripolar) CALL assert(.NOT. open_n .AND. .NOT. cycle_y, "TRIPOLAR and CYCLE_Y or OPEN_N cannot be used simultaneously")

       IF (inputdir  == '' .AND. present(default_inputdir))  inputdir  = default_inputdir
       IF (outputdir == '' .AND. present(default_outputdir)) outputdir = default_outputdir

       SELECT CASE (trim(precision))
       CASE ('SINGLE', 'single', 'REAL4', 'real4', 'r4', 'R4', '4')
          io_kind = 4
       CASE ('DOUBLE', 'double', 'REAL8', 'real8', 'r8', 'R8', '8')
          io_kind = 8
       CASE DEFAULT
          CALL assert(.FALSE., "unsupported PRECISION '"//trim(precision)//"' GEOMETRY namelist")
       END SELECT

       IF (metric_z /= '') THEN
          CALL read_file(tmp(1:dimz), path(inputdir, metric_z), io_kind)
          DO k=1, dimz
             delta_z(k) = delta_z(k)*tmp(k)
          END DO
       END IF

       IF (bathymetry == '' .AND. topography /= '') bathymetry = topography ! for backward compatibility
    END IF


    CALL bcast(io_kind)

#ifdef PARALLEL_MPI
    CALL init_parallel_geometry(path(inputdir, bathymetry), io_kind, gridwise_bathymetry, delta_z)
#else
    isize = dimx
    jsize = dimy
    ksize = dimz
#endif

    CALL bcast(delta_x(1:dimx))
    CALL bcast(delta_y(1:dimy))
    CALL bcast(delta_z(1:dimz))

    CALL replace('$RUNNAME', trim(runname), inputdir)
    CALL replace('$RUNNAME', trim(runname), outputdir)

    CALL bcast(inputdir)
    CALL bcast(outputdir)

    CALL bcast(metric_x)
    CALL bcast(metric_y)
    CALL bcast(metric_z)
    CALL bcast(metric_xy)
    CALL bcast(metric_yx)

    CALL bcast(coriolis_x)
    CALL bcast(coriolis_y)
    CALL bcast(coriolis_z)

    CALL bcast(bathymetry)
    CALL bcast(iceshelf)
    CALL bcast(topomask)

    CALL bcast(gridwise_bathymetry)
    CALL bcast(gridwise_iceshelf)

    CALL bcast(output_dx)
    CALL bcast(output_dy)
    CALL bcast(output_dsx)
    CALL bcast(output_dsy)
    CALL bcast(output_dsz)
    CALL bcast(output_dvol)
    CALL bcast(output_mask)
    CALL bcast(output_mask_r8)

    CALL bcast(coordinate)

    IF (cycle_z) THEN
       IF (rank==0) WRITE(REPORT_UNIT, *) " vertically-cyclic condition is applied"
       CALL assert(rigid_lid, "CYCLE_Z can be applied only for the RIGID_LID mode")
       CALL assert(buoyancy_scheme/=1, "CYCLE_Z requires BUOYANCY_SCHEME = 0 or 2")
    END IF

    regions(0)%defined = .TRUE.
    regions(0)%x_start = 1
    regions(0)%y_start = 1
    regions(0)%z_start = 1
    regions(0)%x_end   = dimx
    regions(0)%y_end   = dimy
    regions(0)%z_end   = dimz
    regions(0)%x_size  = dimx
    regions(0)%y_size  = dimy
    regions(0)%z_size  = dimz

    CALL read_subregion_namelist

#ifdef MPIIO
    CALL read_fileview_namelist
#endif

!$ IF (rank==0) WRITE(REPORT_UNIT, *) "OMP_NUM_THREADS=", nthreads

    ksize_t(0:nthreads-1) = ksize/nthreads

    DO n=1, ksize - ksize_t(0)*nthreads
       ksize_t(nthreads-n) = ksize_t(nthreads-n) + 1
    END DO
    CALL assert(sum(ksize_t(0:nthreads-1))==ksize, "invalid KSIZE_T")

    DO n=0, nthreads-1
       kstr_t(n) = ksize - sum(ksize_t(n:nthreads-1)) + 1
       kend_t(n) = sum(ksize_t(0:n))
    END DO

#ifdef DEBUG
!$  IF (rank==0) THEN
!$OMP PARALLEL PRIVATE(n)
!$    n = omp_get_thread_num()
!$    WRITE(REPORT_UNIT, *) "thred num=",n,":  k=", kstr_t(n),":", kend_t(n)
!$OMP END PARALLEL
!$  END IF
#endif

    ALLOCATE(dx(1-slv:isize+slv, 1-slv:jsize+slv))
    ALLOCATE(dy(1-slv:isize+slv, 1-slv:jsize+slv))
    ALLOCATE(dz(1-slv:ksize+slv))

    dx(:,:)   = 1.0
    dy(:,:)   = 1.0
    dz(:)     = 1.0

    ALLOCATE(corx(1-slv:isize+slv, 1-slv:jsize+slv))
    ALLOCATE(cory(1-slv:isize+slv, 1-slv:jsize+slv))
    ALLOCATE(corz(1-slv:isize+slv, 1-slv:jsize+slv))
    corx(:,:) = 0.0
    cory(:,:) = 0.0
    corz(:,:) = 0.0

    IF (coordinate(:7)=='LONLATC' .OR. coordinate(:7)=='lonlatc') THEN
       lat = read_literal(coordinate(8:))
       DO j=1, dimy
          delta_y(j) = delta_y(j)*cos(pi*(lat+0.5*delta_y(j))/180.0)
          lat = lat + delta_y(j)
       END DO
       lat = read_literal(coordinate(8:))
       IF (jcoord > 0) lat = lat + sum(delta_y(1:jsize*jcoord))
       DO j=1, jsize
          lat = lat + 0.5*delta_y(jsize*jcoord+j)
          dx(:,j) = earth_radius * cos(pi*lat/180.0) * pi/180.0

          corz(:,j) = 4*pi * (366/365) / day_sec * sin(pi*lat/180.0)
          lat = lat + 0.5*delta_y(jsize*jcoord+j)
       END DO
       dy(:,:) = earth_radius * pi/180.0
    ELSE IF (coordinate(:6)=='LONLAT' .OR. coordinate(:6)=='lonlat') THEN
       lat = read_literal(coordinate(7:))
       IF (jcoord > 0) lat = lat + sum(delta_y(1:jsize*jcoord))
       DO j=1, jsize
          lat = lat + 0.5*delta_y(jsize*jcoord+j)
          dx(:,j) = earth_radius * cos(pi*lat/180.0) * pi/180.0

          corz(:,j) = 4*pi * (366/365) / day_sec * sin(pi*lat/180.0)
          lat = lat + 0.5*delta_y(jsize*jcoord+j)
       END DO
       dy(:,:) = earth_radius * pi/180.0
    ELSE IF (coordinate /= '') THEN
       CALL assert(.FALSE., "unsupported COORDINATE '"//trim(coordinate)//"'")
    END IF

    IF (trim(metric_x) /= '') THEN
       CALL read_data_2d(dx(1:isize, 1:jsize), path(inputdir, metric_x), KIND=io_kind)
    END IF
    DO i=1, isize
       dx(i,:) = dx(i,:) * delta_x(icoord*isize+i)
    END DO

    IF (trim(metric_y) /= '') THEN
       CALL read_data_2d(dy(1:isize, 1:jsize), path(inputdir, metric_y), KIND=io_kind)
    END IF
    DO j=1, jsize
       dy(:,j) = dy(:,j) * delta_y(jcoord*jsize+j)
    END DO

    DO k=1-slv, ksize+slv
       dz(k) = delta_z(max(1,min(dimz,kcoord*ksize+k)))
    END DO

    dz_min = minval(dz(1:ksize))
    dz_max = maxval(dz(1:ksize))
    CALL vmin(dz_min, all=.TRUE.)
    CALL vmax(dz_max, all=.TRUE.)

    DO i=1, slv
       dx(1-i,:) = dx(1,:)
       dy(1-i,:) = dy(1,:)
       dx(isize+i,:) = dx(isize,:)
       dy(isize+i,:) = dy(isize,:)
    END DO

    DO j=1, slv
       dx(:,1-j) = dx(:,1)
       dy(:,1-j) = dy(:,1)
       dx(:,jsize+j) = dx(:,jsize)
       dy(:,jsize+j) = dy(:,jsize)

       corz(:,1-j)     = corz(:,1)
       corz(:,jsize+j) = corz(:,jsize)
    END DO

    CALL update_boundary(dx, method='STD')
    CALL update_boundary(dy, method='STD')

    IF (trim(coriolis_x) /= '') CALL read_data(corx(1:isize, 1:jsize), path(inputdir, coriolis_x), KIND=io_kind)
    IF (trim(coriolis_y) /= '') CALL read_data(cory(1:isize, 1:jsize), path(inputdir, coriolis_y), KIND=io_kind)
    IF (trim(coriolis_z) /= '') CALL read_data(corz(1:isize, 1:jsize), path(inputdir, coriolis_z), KIND=io_kind)

    CALL update_boundary(corx, method='STD')
    CALL update_boundary(cory, method='STD')
    CALL update_boundary(corz, method='STD')

    ALLOCATE(depth(-slv:ksize+slv))
    DO k=-slv, ksize+slv
       IF (kcoord*ksize+k >= dimz) THEN
          depth(k) = 0.0
       ELSE IF (kcoord*ksize+k <= 0) THEN
          depth(k) = sum(delta_z(1:dimz))
       ELSE
          depth(k) = sum(delta_z(kcoord*ksize+k+1:dimz))
       END IF
    END DO

    CALL assert(.NOT. (use_landwater .AND. trim(bathymetry) == ''), "BATHYMETRY is mandatory when LANDWATER is enabled")
    CALL assert(.NOT. (use_landwater .AND. trim(iceshelf)   /= ''), "ICESHELF is not allowed when LANDWATER is enabled")

    ALLOCATE(h_bathymetry(1-slv:isize+slv, 1-slv:jsize+slv))
    h_bathymetry(:,:) = 0.0D0
    h_bathymetry(1:isize, 1:jsize) = -sum(delta_z(1:dimz))

    IF (trim(bathymetry) /= '') THEN
       CALL read_data_2d(h_bathymetry(1:isize, 1:jsize), path(inputdir, bathymetry), KIND=io_kind)

       topomin = minval(h_bathymetry(1:isize,1:jsize))
       CALL gmin(topomin, all=.TRUE.)

       IF (gridwise_bathymetry .OR. topomin >= 0.0) THEN
          DO j=1, jsize
          DO i=1, isize
             h_bathymetry(i,j) = h_gridwise(h_bathymetry(i,j))
          END DO
          END DO
       END IF
    END IF

    CALL update_boundary(h_bathymetry)

    ALLOCATE(h_iceshelf(1-slv:isize+slv, 1-slv:jsize+slv))
    h_iceshelf(:,:) = UNDEF

    IF (trim(iceshelf) /= '') THEN
       CALL read_data_2d(h_iceshelf(1:isize, 1:jsize), path(inputdir, iceshelf), KIND=io_kind)

       topomin = minval(h_iceshelf(1:isize,1:jsize))
       CALL gmin(topomin, all=.TRUE.)

       IF (gridwise_iceshelf .OR. topomin >= 0.0) THEN
          DO j=1, jsize
          DO i=1, isize
             h_iceshelf(i,j) = h_gridwise(h_iceshelf(i,j))
          END DO
          END DO
       END IF
       iceshelf_coupled = .TRUE.
    END IF

    CALL update_boundary(h_iceshelf)

    ALLOCATE(imask2d(-1:isize+2,-1:jsize+2))
    ALLOCATE(imask3d(-1:isize+2,-1:jsize+2,-1:ksize+2))

    ALLOCATE(lmask2d(-1:isize+2,-1:jsize+2))
    ALLOCATE(lmask3d(-1:isize+2,-1:jsize+2,-1:ksize+2))

    ALLOCATE(maskindex(-1:ksize+2))

    ALLOCATE(surface_k(1-slv:isize+slv, 1-slv:jsize+slv))
    ALLOCATE(surface_flag(1-slv:isize+slv, 1-slv:jsize+slv))

    ALLOCATE(bottom_k(1-slv:isize+slv, 1-slv:jsize+slv))
    ALLOCATE(bottom_flag(1-slv:isize+slv, 1-slv:jsize+slv))

    ALLOCATE(dz_ref(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv))

    ALLOCATE(dz_star(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv))
    ALLOCATE(cz_star(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv))

    IF (use_landwater) THEN
       ALLOCATE(lwflag(1-slv:isize+slv, 1-slv:jsize+slv))
       ALLOCATE(lwdried(1-slv:isize+slv, 1-slv:jsize+slv))
       ALLOCATE(landelev(1-slv:isize+slv, 1-slv:jsize+slv))
    END IF

    ALLOCATE(delta2(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv))

    CALL make_mask

    IF (output_mask .OR. output_mask_r8) THEN
       IF (output_mask_r8) THEN
#ifdef PARALLEL_MPI
          IF (remove_masked_pe .AND. rank==0) CALL dump_file(0.0, path(outputdir, 'MASK'), dimx*dimy*dimz, stride=dimx*dimy, kind=8)
#endif
          CALL write_data_3d(REAL(imask3d(1:isize,1:jsize,1:ksize)), path(outputdir, 'MASK'), kind=8)
       ELSE
#ifdef PARALLEL_MPI
          IF (remove_masked_pe .AND. rank==0)  CALL dump_file(0.0, path(outputdir, 'MASK'), dimx*dimy*dimz, stride=dimx*dimy, kind=1)
#endif
          CALL write_data_3d(REAL(imask3d(1:isize,1:jsize,1:ksize)), path(outputdir, 'MASK'), kind=1)
       END IF
    END IF

    ALLOCATE(dx1(1-slv:isize+slv, 1-slv:jsize+slv-1))
    DO j=1-slv, jsize+slv-1
    DO i=1-slv, isize+slv
       dx1(i,j) = (dx(i,j)*dy(i,j+1) + dx(i,j+1)*dy(i,j))/(dy(i,j)+dy(i,j+1))
    END DO
    END DO

    ALLOCATE(dy1(1-slv:isize+slv-1, 1-slv:jsize+slv))
    DO j=1-slv, jsize+slv
    DO i=1-slv, isize+slv-1
       dy1(i,j) = (dy(i,j)*dx(i+1,j) + dy(i+1,j)*dx(i,j))/(dx(i,j)+dx(i+1,j))
    END DO
    END DO

    ALLOCATE(idx0(1-slv:isize+slv, 1-slv:jsize+slv))
    DO j=1-slv, jsize+slv
    DO i=1-slv, isize+slv
       idx0(i,j) = 1.0/dx(i,j)
    END DO
    END DO

    ALLOCATE(idy0(1-slv:isize+slv,  1-slv:jsize+slv))
    DO j=1-slv, jsize+slv
    DO i=1-slv, isize+slv
       idy0(i,j) = 1.0/dy(i,j)
    END DO
    END DO

    ALLOCATE(idz0(1-slv:ksize+slv))
    DO k=1-slv, ksize+slv
       idz0(k) = 1.0/dz(k)
    END DO

    ALLOCATE(idx1(1-slv:isize+slv-1, 1-slv:jsize+slv))
    DO j=1-slv, jsize+slv
    DO i=1-slv, isize+slv-1
       idx1(i,j) = 2.0/(dx(i,j)+dx(i+1,j))
    END DO
    END DO

    ALLOCATE(idy1(1-slv:isize+slv,  1-slv:jsize+slv-1))
    DO j=1-slv, jsize+slv-1
    DO i=1-slv, isize+slv
       idy1(i,j) = 2.0/(dy(i,j)+dy(i,j+1))
    END DO
    END DO

    ALLOCATE(idz1(1-slv:ksize+slv-1))
    DO k=1-slv, ksize+slv-1
       idz1(k) = 2.0/(dz(k)+dz(k+1))
    END DO

    ALLOCATE(idx2(0:isize+1, -1:jsize+2))
    DO j=-1, jsize+2
    DO i= 0, isize+1
       idx2(i,j) = 1.0/(dx(i,j)+0.5*dx(i+1,j)+0.5*dx(i-1,j))
    END DO
    END DO

    ALLOCATE(idy2(-1:isize+2, 0:jsize+1))
    DO j= 0, jsize+1
    DO i=-1, isize+2
       idy2(i,j) = 1.0/(dy(i,j)+0.5*dy(i,j+1)+0.5*dy(i,j-1))
    END DO
    END DO

    ALLOCATE(idz2(0:ksize+1))
    DO k=0, ksize+1
       idz2(k) = 1.0/(dz(k)+0.5*dz(k+1)+0.5*dz(k-1))
    END DO

    ALLOCATE(dsz(1-slv:isize+slv,1-slv:jsize+slv))
    DO j=1-slv, jsize+slv
    DO i=1-slv, isize+slv
       dsz(i,j) = dx(i,j)*dy(i,j)
    END DO
    END DO

    ALLOCATE(wct_ref(1-slv:isize+slv, 1-slv:jsize+slv))
    wct_ref(:,:) = 0.0D0

    DO k=1, ksize
    DO j=1, jsize
    DO i=1, isize
       wct_ref(i,j) = wct_ref(i,j) + dz_ref(i,j,k)*imask3d(i,j,k)
    END DO
    END DO
    END DO

    CALL vsum(wct_ref, ALL=.TRUE.)
    CALL update_boundary(wct_ref)

    ALLOCATE(metxy(1-slv:isize+slv, 1-slv:jsize+slv))
    ALLOCATE(metyx(1-slv:isize+slv, 1-slv:jsize+slv))
    metxy(:,:) = 0.0
    metyx(:,:) = 0.0

    IF (trim(metric_xy) /= '') THEN
       CALL read_data(metxy(1:isize, 1:jsize), path(inputdir, metric_xy), KIND=io_kind)
    ELSE
       DO j=0, jsize+1
       DO i=0, isize+1
          metxy(i,j) = 0.5 * (log(dx(i,j+1)/dx(i,j))*idy1(i,j) + log(dx(i+1,j+1)/dx(i+1,j))*idy1(i+1,j))
       END DO
       END DO
    END IF

    IF (trim(metric_yx) /= '') THEN
       CALL read_data(metyx(1:isize, 1:jsize), path(inputdir, metric_yx), KIND=io_kind)
    ELSE
       DO j=0, jsize+1
       DO i=0, isize+1
          metyx(i,j) = 0.5 * (log(dy(i+1,j)/dy(i,j))*idx1(i,j) + log(dy(i+1,j+1)/dy(i,j+1))*idx1(i,j+1))
       END DO
       END DO
    END IF

    CALL update_boundary(metxy, method='STD')
    CALL update_boundary(metyx, method='STD')

    ALLOCATE(dzindex(-1:ksize+2))
    dzindex(-1) = -1
    dzindex(0)  = 0
    dzindex(1)  = 1
    DO k=2, ksize
       dzindex(k) = dzindex(k-1)

       DO j=1-slv, jsize+slv
       DO i=1-slv, isize+slv
          IF (surface_k(i,j)==k .OR. surface_k(i,j)==k-1) dzindex(k) = k
          IF (dz_ref(i,j,k) /= dz_ref(i,j,k-1))     dzindex(k) = k
       END DO
       END DO
    END DO
    dzindex(ksize+1) = ksize+1
    dzindex(ksize+2) = ksize+2

    ALLOCATE(ssh(1-slv:isize+slv, 1-slv:jsize+slv))
    ssh(:,:) = 0.0

    ALLOCATE(ssh_old(1-slv:isize+slv, 1-slv:jsize+slv))
    ssh_old(:,:) = 0.0

    ALLOCATE(dvol( 1-slv:isize+slv,   1-slv:jsize+slv,   1-slv:ksize+slv))
    ALLOCATE(dsx(  1-slv:isize+slv-1, 1-slv:jsize+slv,   1-slv:ksize+slv))
    ALLOCATE(dsy(  1-slv:isize+slv,   1-slv:jsize+slv-1, 1-slv:ksize+slv))

    ALLOCATE(dvol_ref( 1-slv:isize+slv,   1-slv:jsize+slv,   1-slv:ksize+slv))
    ALLOCATE(dvol_old( 1-slv:isize+slv,   1-slv:jsize+slv,   1-slv:ksize+slv))

    ALLOCATE(dsx_ref(  1-slv:isize+slv-1, 1-slv:jsize+slv,   1-slv:ksize+slv))
    ALLOCATE(dsy_ref(  1-slv:isize+slv,   1-slv:jsize+slv-1, 1-slv:ksize+slv))

    ALLOCATE(dsx_old(  1-slv:isize+slv-1, 1-slv:jsize+slv,   1-slv:ksize+slv))
    ALLOCATE(dsy_old(  1-slv:isize+slv,   1-slv:jsize+slv-1, 1-slv:ksize+slv))

    ALLOCATE(dsx2d(1-slv:isize+slv-1, 1-slv:jsize+slv))
    ALLOCATE(dsy2d(1-slv:isize+slv,   1-slv:jsize+slv-1))

    dvol(:,:,:) = 0.0
    dsx(:,:,:)  = 0.0
    dsy(:,:,:)  = 0.0

    ALLOCATE(cext(1-slv:isize+slv, 1-slv:jsize+slv))
    cext(:,:) = 0.0

    IF (delta2_sml93) THEN
! square of the length scale used in the Smagorinsky model considering grid-anisotoropy,
! c.f. Scotti-Meneveau-Lilly (1993), Phys. Fluids.
       DO k=1-slv, ksize+slv
       DO j=1-slv, jsize+slv
       DO i=1-slv, isize+slv
          d1 = max(dx(i,j), dy(i,j), dz(k))
          d3 = min(dx(i,j), dy(i,j), dz(k))
          d2 = dx(i,j) + dy(i,j) + dz(k) - d1 - d3

          delta2(i,j,k) = ((d1*d2*d3)**(1.0D0/3) &
               * cosh(sqrt(4.0D0/27 * (log(d2/d1)**2 + log(d3/d1)**2 - log(d2/d1)*log(d3/d1)))))**2
       END DO
       END DO
       END DO
    ELSE IF (delta2_cube) THEN
       DO k=1-slv, ksize+slv
       DO j=1-slv, jsize+slv
       DO i=1-slv, isize+slv
          delta2(i,j,k) = (dx(i,j)*dy(i,j)*dz(k))**(2.0/3)
       END DO
       END DO
       END DO
    ELSE
       DO k=1-slv, ksize+slv
       DO j=1-slv, jsize+slv
       DO i=1-slv, isize+slv
          delta2(i,j,k) = (dx(i,j)**2 + dy(i,j)**2) / 2.0
       END DO
       END DO
       END DO
    END IF

    CALL update_geometry

    dvol_ref(:,:,:) = dvol(:,:,:)
    dvol_old(:,:,:) = dvol(:,:,:)

    dsx_ref(:,:,:)  = dsx(:,:,:)
    dsx_old(:,:,:)  = dsx(:,:,:)

    dsy_ref(:,:,:)  = dsy(:,:,:)
    dsy_old(:,:,:)  = dsy(:,:,:)

    total_area = 0.0
    IF (vrank==0) THEN
       DO j=1, jsize
       DO i=1, isize
          total_area = total_area + dsz(i,j)*imask2d(i,j)
       END DO
       END DO
    END IF
    CALL gsum(total_area, all=.TRUE.)

    ALLOCATE(layer_area(1-slv:ksize+slv))
    DO k=1-slv, ksize+slv
       layer_area(k) = sum(imask3d(1:isize,1:jsize,k)*dsz(1:isize,1:jsize))
    END DO
    CALL hsum(layer_area, all=.TRUE.)

    ALLOCATE(cint(1-slv:isize+slv, 1-slv:jsize+slv))
    cint(:,:) = 0.0

    IF (rank==0) WRITE(REPORT_UNIT, '(A, ES10.3, A)') "total sea surface area = ", total_area, " [m^2]"

    CALL report_openarea

    IF (output_dx)   CALL write_data_2d(  dx(1:isize,1:jsize), path(outputdir, 'DX'))
    IF (output_dy)   CALL write_data_2d(  dy(1:isize,1:jsize), path(outputdir, 'DY'))
    IF (output_dsx)  CALL write_data_3d( dsx(1:isize,1:jsize,1:ksize)*imask3d(1:isize,1:jsize,1:ksize)*imask3d(2:isize+1,1:jsize,1:ksize), path(outputdir, 'DSX'))
    IF (output_dsy)  CALL write_data_3d( dsy(1:isize,1:jsize,1:ksize)*imask3d(1:isize,1:jsize,1:ksize)*imask3d(1:isize,2:jsize+1,1:ksize), path(outputdir, 'DSY'))
    IF (output_dsz)  CALL write_data_2d( dsz(1:isize,1:jsize)        *imask2d(1:isize,1:jsize),         path(outputdir, 'DSZ'))
    IF (output_dvol) CALL write_data_3d(dvol(1:isize,1:jsize,1:ksize)*imask3d(1:isize,1:jsize,1:ksize), path(outputdir, 'DVOL'))

  CONTAINS
    REAL(8) PURE FUNCTION h_gridwise(x)
      REAL(8), INTENT(IN) :: x

      INTEGER :: itmp
      REAL(8) :: rtmp
      INTEGER :: k

      h_gridwise = min(max(x, 0.0D0), 1.0D0*dimz)

      itmp = INT(h_gridwise)
      rtmp = h_gridwise - itmp

      h_gridwise = 0.0
      DO k=dimz, itmp+1, -1
         h_gridwise = h_gridwise - delta_z(k)
      END DO
      IF (partial_step) h_gridwise = h_gridwise + rtmp * delta_z(itmp+1)
    END FUNCTION h_gridwise


    SUBROUTINE report_openarea
      REAL(8) :: a(4)

      a(:) = 0.0

      IF (icoord==0) THEN
         DO k=1, ksize
         DO j=1, jsize
            a(1) = a(1) + dsx_ref(0,j,k)*imask3d(1,j,k)
         END DO
         END DO
      END IF

      IF (icoord==ipes-1) THEN
         DO k=1, ksize
         DO j=1, jsize
            a(2) = a(2) + dsx_ref(isize,j,k)*imask3d(isize,j,k)
         END DO
         END DO
      END IF

      IF (jcoord==0) THEN
         DO k=1, ksize
         DO i=1, isize
            a(3) = a(3) + dsy_ref(i,0,k)*imask3d(i,1,k)
         END DO
         END DO
      END IF

      IF (jcoord==jpes-1) THEN
         DO k=1, ksize
         DO i=1, isize
            a(4) = a(4) + dsy_ref(i,jsize,k)*imask3d(i,jsize,k)
         END DO
         END DO
      END IF

      CALL gsum(a)

      IF (rank==0 .AND. open_w) WRITE(REPORT_UNIT, '(A,ES10.3,A)') " West open boundary cross-section area = ", a(1), " [m^2]"
      IF (rank==0 .AND. open_e) WRITE(REPORT_UNIT, '(A,ES10.3,A)') " East open boundary cross-section area = ", a(2), " [m^2]"
      IF (rank==0 .AND. open_s) WRITE(REPORT_UNIT, '(A,ES10.3,A)') "South open boundary cross-section area = ", a(3), " [m^2]"
      IF (rank==0 .AND. open_n) WRITE(REPORT_UNIT, '(A,ES10.3,A)') "North open boundary cross-section area = ", a(4), " [m^2]"

    END SUBROUTINE report_openarea

    SUBROUTINE read_subregion_namelist
      INTEGER :: x_start, x_end
      INTEGER :: y_start, y_end
      INTEGER :: z_start, z_end
      INTEGER :: id

      INTEGER :: n, iostat

      NAMELIST / subregion / &
           id,             &
           x_start, x_end, &
           y_start, y_end, &
           z_start, z_end

      DO n=1, MAX_SUBREGION
         regions(n)%defined = .FALSE.
      END DO

      IF (rank==0) REWIND(CONFIG_UNIT)
      DO
         IF (rank==0) THEN
            id = 0
            x_start = 1
            x_end   = dimx
            y_start = 1
            y_end   = dimy
            z_start = 1
            z_end   = dimz

            READ(CONFIG_UNIT, NML=subregion, IOSTAT=iostat, IOMSG=iomsg)
         END IF

         CALL bcast(iostat)
         IF (iostat < 0) EXIT

         CALL assert(iostat == 0, "failed to read SUBREGION namelist", iomsg)

         CALL bcast(id)
         CALL bcast(x_start)
         CALL bcast(x_end)
         CALL bcast(y_start)
         CALL bcast(y_end)
         CALL bcast(z_start)
         CALL bcast(z_end)

         CALL assert(id > 0 .AND. id <= MAX_SUBREGION,  "invalid SUGREGION ID")
         CALL assert(x_start > 0 .AND. x_end <= dimx .AND. x_start <= x_end, "invalid SUBREGION X_START, X_END")
         CALL assert(y_start > 0 .AND. y_end <= dimy .AND. y_start <= y_end, "invalid SUBREGION Y_START, Y_END")
         CALL assert(z_start > 0 .AND. z_end <= dimz .AND. z_start <= z_end, "invalid SUBREGION Z_START, Z_END")
         CALL assert(.NOT. regions(id)%defined, "SUBREGION ID="//trim(format(id))//" is already defined")

         regions(id)%defined = .TRUE.
         regions(id)%x_start = x_start
         regions(id)%y_start = y_start
         regions(id)%z_start = z_start
         regions(id)%x_end   = x_end
         regions(id)%y_end   = y_end
         regions(id)%z_end   = z_end
         regions(id)%x_size  = x_end - x_start + 1
         regions(id)%y_size  = y_end - y_start + 1
         regions(id)%z_size  = z_end - z_start + 1

         IF (rank==0) THEN
            WRITE(REPORT_UNIT, '(A,I0, A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A)') &
                 "SUBREGION ", id, " is defined as [x=", x_start, ":", x_end, &
                                                 ", y=", y_start, ":", y_end, &
                                                 ", z=", z_start, ":", z_end, "]"
         END IF
      END DO
    END SUBROUTINE read_subregion_namelist

!-----------------------------------------------------------------------------------------------------------------------

#ifdef MPIIO
    SUBROUTINE read_fileview_namelist
      INTEGER :: id
      INTEGER :: x_size
      INTEGER :: y_size
      INTEGER :: x_start
      INTEGER :: y_start
      INTEGER :: z_start
      INTEGER(MPI_OFFSET_KIND) :: offset

      INTEGER :: n, iostat

      NAMELIST / fileview / &
           id,              &
           offset,          &
           x_size,  y_size, &
           x_start, y_start

      DO n=1, MAX_FILEVIEW
         regions(n)%defined = .FALSE.
      END DO

      IF (rank==0) REWIND(CONFIG_UNIT)
      DO
         IF (rank==0) THEN
            id = 0
            offset  = 0_MPI_OFFSET_KIND
            x_start = 1
            y_start = 1
            x_size  = dimx
            y_size  = dimy

            READ(CONFIG_UNIT, NML=fileview, IOSTAT=iostat, IOMSG=iomsg)
            CALL assert(iostat<=0, "failed to read FILEVIEW namelist", iomsg)
         END IF

         CALL bcast(iostat)
         IF (iostat < 0) EXIT

         CALL bcast(id)
         CALL bcast(x_size)
         CALL bcast(y_size)
         CALL bcast(x_start)
         CALL bcast(y_start)
         CALL bcast(offset)

         CALL assert(id > 0 .AND. id <= MAX_FILEVIEW,             "invalid FILEVIEW ID")
         CALL assert(offset  >= 0,                                "invalid OFFSET for FILEVIEW ID="//format(id))
         CALL assert(x_start >= 0 .AND. x_start+dimx-1 <= x_size, "invalid X_START and/or X_SIZE for FILEVIEW ID="//format(id))
         CALL assert(y_start >= 0 .AND. y_start+dimy-1 <= y_size, "invalid Y_START and/or Y_SIZE for FILEVIEW ID="//format(id))

         CALL assert(.NOT. views(id)%defined, "FILEVIEW ID="//trim(format(id))//" is already defined")

         CALL init_fileview(id, offset, x_size, y_size, x_start, y_start)

      END DO
    END SUBROUTINE read_fileview_namelist
#endif

!-----------------------------------------------------------------------------------------------------------------------

    SUBROUTINE make_mask
      INTEGER :: itmp
      REAL(8) :: rtmp
      LOGICAL :: flag
      INTEGER :: i, j, k
      INTEGER :: ierr

      REAL(4) :: mask(-1:isize+2, -1:jsize+2, -1:ksize+2)

      REAL(4) :: th

!$OMP PARALLEL WORKSHARE
      mask(:,:,:) = 0.0
      mask(1:isize,1:jsize,1:ksize) = 1.0
!$OMP END PARALLEL WORKSHARE

      IF (topomask /= '') THEN
         CALL read_data_3d(mask(1:isize,1:jsize,1:ksize), path(inputdir, topomask), kind=1)
!$OMP PARALLEL DO
         DO k=1, ksize
         DO j=1, jsize
         DO i=1, isize
            IF (mask(i,j,k)/=0.0) mask(i,j,k) = 1.0
         END DO
         END DO
         END DO
      END IF

      DO k=1-slv, ksize+slv
         dz_ref(:,:,k) = dz(k)
      END DO

      DO j=1-slv, jsize+slv
      DO i=1-slv, isize+slv
         DO k=-1, ksize+2
            itmp = k
            IF (-depth(k) > h_bathymetry(i,j)) EXIT
         END DO

         IF (vrank==0 .AND. itmp>=ksize) THEN
            rtmp = 1.0D0
         ELSE
            rtmp = (- min(h_bathymetry(i,j),0.0D0) - depth(itmp))/dz(itmp)
         END IF

         DO k=1, ksize
            IF (k < itmp) mask(i,j,k) = 0.0
         END DO

         IF (itmp >= 0 .AND. itmp <= ksize+1) THEN
            IF (partial_step) THEN
               th = max(partial_step_threshold, 0.125*dz_min/dz(itmp))

               IF (rtmp < th) THEN
                  mask(i,j,itmp) = 0.0
                  dz_ref(i,j,itmp+1) = dz_ref(i,j,itmp+1) + dz(itmp)*rtmp
                  dz_ref(i,j,itmp)   = dz_ref(i,j,itmp)   - dz(itmp)*rtmp
               ELSE
                  dz_ref(i,j,itmp)   = dz_ref(i,j,itmp)   - dz(itmp)*(1.0-rtmp)
                  dz_ref(i,j,itmp-1) = dz_ref(i,j,itmp-1) + dz(itmp)*(1.0-rtmp)
               END IF
            ELSE
               IF (rtmp < 1.0/4.0) THEN
                  mask(i,j,itmp) = 0.0
               END IF
            END IF
         END IF

         IF (h_iceshelf(i,j) == UNDEF) CYCLE

         h_iceshelf(i,j) =  min(h_iceshelf(i,j), 0.0)

         DO k=ksize+2, -1, -1
            itmp = k
            IF (-depth(k-1) < h_iceshelf(i,j)) EXIT
         END DO

         rtmp = (depth(itmp-1) + h_iceshelf(i,j))/dz(itmp)

         DO k=1, ksize
            IF (k > itmp) mask(i,j,k) = 0.0
         END DO

         IF (itmp >= 0 .AND. itmp <= ksize+1) THEN
            IF (partial_step) THEN
               th = max(partial_step_threshold, 0.125*dz_min/dz(itmp))

               IF (rtmp < th) THEN
                  mask(i,j,itmp) = 0.0
                  dz_ref(i,j,itmp-1) = dz_ref(i,j,itmp-1) + dz(itmp)*rtmp
                  dz_ref(i,j,itmp)   = dz_ref(i,j,itmp)   - dz(itmp)*rtmp
               ELSE
                  dz_ref(i,j,itmp)   = dz_ref(i,j,itmp)   - dz(itmp)*(1.0-rtmp)
                  dz_ref(i,j,itmp+1) = dz_ref(i,j,itmp+1) + dz(itmp)*(1.0-rtmp)
                  IF (dz_ref(i,j,itmp) < dz(itmp)*0.25) mask(i,j,itmp) = 0.0
               END IF
            ELSE
               IF (rtmp < 1.0/4.0) THEN
                  mask(i,j,itmp) = 0.0
               END IF
            END IF
         END IF

      END DO
      END DO

      IF (iceshelf_coupled) THEN
         ALLOCATE(isflag(0:isize+1,0:jsize+1,0:ksize+1))
         isflag(:,:,:) = .FALSE.
         DO k=1, ksize
         DO j=1, jsize
         DO i=1, isize
            isflag(i,j,k) = (mask(i,j,k)==0.0 .AND. -depth(k) > h_iceshelf(i,j))
         END DO
         END DO
         END DO
         CALL update_boundary(isflag)
      END IF

      IF (use_landwater) THEN
         IF (vrank==0) THEN
            DO j=1, jsize
            DO i=1, isize
               IF (h_bathymetry(i,j) < landwater_limit) THEN
                  mask(i,j,ksize) = 1.0
               ELSE
                  mask(i,j,ksize) = 0.0
               END IF
            END DO
            END DO
         END IF

         lwflag(:,:)  = .FALSE.
         lwdried(:,:) = .FALSE.

         landelev(:,:)  = 0.0D0

         IF (vrank==0) THEN
            DO j=1-slv, jsize+slv
            DO i=1-slv, isize+slv
               landelev(i,j) = max(h_bathymetry(i,j), 0.0D0)

               IF (mask(i,j,ksize)==0.0) CYCLE

               IF (h_bathymetry(i,j) > -dz(ksize)) THEN
                  lwflag(i,j) = .TRUE.
                  dz_ref(i,j,ksize) = max(-h_bathymetry(i,j), 0.0)
               END IF

               lwdried(i,j) = (dz_ref(i,j,ksize) < landwater_hdry)
            END DO
            END DO
         END IF
      END IF

      CALL update_boundary(mask)

      IF (fill_chokepoint) THEN
         DO k=1, ksize
            flag = .TRUE.
            DO WHILE (flag)
               flag = .FALSE.
               DO j=1, jsize
               DO i=1, isize
                  IF (mask(i,j,k)==1.0 .AND.  mask(i+1,j,k)+mask(i-1,j,k)+mask(i,j+1,k)+mask(i,j-1,k) <= 1.0) THEN
                     mask(i,j,k) = 0.0
                     flag = .TRUE.
                  END IF
               END DO
               END DO
#ifdef PARALLEL_MPI
               CALL mpi_allreduce(MPI_IN_PLACE, flag, 1, MPI_LOGICAL, MPI_LOR, comm, ierr)
#endif
               CALL update_boundary(mask(:,:,k))
            END DO
         END DO
         CALL update_boundary(mask)
      END IF

      IF (open_w .AND. icoord==0) THEN
         DO i=1, slv
            dz_ref(1-i,:,:) = dz_ref(1,:,:)
         END DO
      END IF

      IF (open_e .AND. icoord==ipes-1) THEN
         DO i=1, slv
            dz_ref(isize+i,:,:) = dz_ref(isize,:,:)
         END DO
      END IF

      IF (open_s .AND. jcoord==0) THEN
         DO j=1, slv
            dz_ref(:,1-j,:) = dz_ref(:,1,:)
         END DO
      END IF

      IF (open_n .AND. jcoord==jpes-1) THEN
         DO j=1, slv
            dz_ref(:,jsize+j,:) = dz_ref(:,jsize,:)
         END DO
      END IF

      CALL update_boundary(dz_ref)


      maskindex(-1) = -1
      maskindex( 0) = 0
      maskindex( 1) = 1
      DO k=2, ksize+2
         IF (all(mask(:,:,k) == mask(:,:,k-1))) THEN
            maskindex(k) = maskindex(k-1)
         ELSE
            maskindex(k) = k
         END IF
      END DO

      surface_k(:,:) = 0
      DO k=1, ksize
      DO j=1-slv, jsize+slv
      DO i=1-slv, isize+slv
         IF (mask(i,j,k+1)==0.0 .AND. mask(i,j,k)==1.0) THEN
            surface_k(i,j) = ksize*kcoord + k
         END IF
      END DO
      END DO
      END DO
      CALL vmax(surface_k, all=.TRUE.)

      IF (cycle_z) THEN
         DO j=1-slv, jsize+slv
         DO i=1-slv, isize+slv
            IF (surface_k(i,j)==0) surface_k(i,j) = dimz
         END DO
         END DO
      END IF

      surface_flag(:,:) = .FALSE.
      DO j=1-slv, jsize+slv
      DO i=1-slv, isize+slv
         surface_k(i,j) = surface_k(i,j) - ksize*kcoord
         surface_flag(i,j) = (surface_k(i,j) >= 1 .AND. surface_k(i,j) <= ksize)
      END DO
      END DO

      bottom_k(:,:) = dimz+1
      DO k=ksize, 1, -1
      DO j=1-slv, jsize+slv
      DO i=1-slv, isize+slv
         IF (mask(i,j,k-1)==0.0 .AND. mask(i,j,k)==1.0) THEN
            bottom_k(i,j) = ksize*kcoord + k
         END IF
      END DO
      END DO
      END DO
      CALL vmin(bottom_k, all=.TRUE.)

      IF (cycle_z) THEN
         DO j=1-slv, jsize+slv
         DO i=1-slv, isize+slv
            IF (bottom_k(i,j)==dimz+1) bottom_k(i,j) = 1
         END DO
         END DO
      END IF

      bottom_flag(:,:) = .FALSE.
      DO j=1-slv, jsize+slv
      DO i=1-slv, isize+slv
         bottom_k(i,j) = bottom_k(i,j) - ksize*kcoord
         bottom_flag(i,j) = (bottom_k(i,j) >= 1 .AND. bottom_k(i,j) <= ksize)
      END DO
      END DO


      imask3d(:,:,:) = INT(mask(-1:isize+2, -1:jsize+2, -1:ksize+2))

      DO k=-1, ksize+2
      DO j=-1, jsize+2
      DO i=-1, isize+2
         CALL assert(imask3d(i,j,k)==0 .OR. imask3d(i,j,k)==1, "IMASK3D is wrong")
         lmask3d(i,j,k) = (imask3d(i,j,k)==1)
      END DO
      END DO
      END DO

      imask2d(:,:)  = maxval(imask3d, dim=3)

#ifdef PARALLEL3D
      CALL mpi_allreduce(MPI_IN_PLACE, imask2d, size(imask2d), MPI_BYTE, MPI_BOR, vcomm, ierr)
#endif

      DO j=-1, jsize+2
      DO i=-1, isize+2
         CALL assert(imask2d(i,j)==0 .OR. imask2d(i,j)==1, "IMASK2D is wrong")
         lmask2d(i,j) = (imask2d(i,j)==1)
      END DO
      END DO


      n_cells = sum(INT(imask3d(1:isize,1:jsize,1:ksize)))
#ifdef PARALLEL_MPI
      CALL mpi_allreduce(MPI_IN_PLACE, n_cells, 1, MPI_INTEGER8, MPI_SUM, comm, ierr)
#endif

#ifdef DEBUG
      DO j=1, jsize
      DO i=1, isize
         itmp = 0
         IF (surface_flag(i,j)) itmp = 1
         CALL mpi_allreduce(MPI_IN_PLACE, itmp, 1, MPI_INTEGER, MPI_SUM, vcomm, ierr)
         CALL assert(itmp <= 1, "SURFACE_FLAG was set in two different layers")
         CALL assert(logical(lmask2d(i,j)) .EQV. itmp==1, "SURFACE_FLAG does not match LMASK2D")

         itmp = 0
         IF (bottom_flag(i,j)) itmp = 1
         CALL mpi_allreduce(MPI_IN_PLACE, itmp, 1, MPI_INTEGER, MPI_SUM, vcomm, ierr)
         CALL assert(itmp <= 1, "BOTTOM_FLAG was set in two different layers")
         CALL assert(logical(lmask2d(i,j)) .EQV. itmp==1, "BOTTOM_FLAG does not match LMASK2D")
      END DO
      END DO
#endif

      dz_star(:,:,:) = dz_ref(:,:,:)
      cz_star(:,:,:) = 0.0

      IF (vcoord_zstar) THEN
         DO j=1-slv, jsize+slv
         DO i=1-slv, isize+slv
            rtmp = sum(dz(1:ksize)*mask(i,j,1:ksize))
            CALL vsum(rtmp, all=.TRUE.)
            IF (rtmp == 0.0) CYCLE

            DO k=1-slv, ksize+slv
               cz_star(i,j,k) = dz(k) / rtmp
            END DO
         END DO
         END DO

      ELSE
         DO j=1-slv, jsize+slv
         DO i=1-slv, isize+slv
            IF (surface_flag(i,j)) THEN
               cz_star(i,j,surface_k(i,j))   =  1.0
               cz_star(i,j,surface_k(i,j)+1) = -1.0
            END IF
         END DO
         END DO
      END IF

    END SUBROUTINE make_mask
  END SUBROUTINE init_geometry

!-----------------------------------------------------------------------------------------------------------------------

#ifdef PARALLEL_MPI
  SUBROUTINE init_parallel_geometry(bathymetry, kind, gridwise, delta_z)
    CHARACTER(*), INTENT(IN) :: bathymetry
    INTEGER,      INTENT(IN) :: kind
    LOGICAL,      INTENT(IN) :: gridwise
    REAL(8),      INTENT(IN) :: delta_z(maxdim)

    INTEGER :: coords(3), ierr, itmp

    REAL(8), ALLOCATABLE :: ratio(:)
    REAL(8), ALLOCATABLE :: topo(:,:)
    REAL(8), ALLOCATABLE :: depth(:)

    REAL(8) :: topomin

    INTEGER :: slv0, slv1
    INTEGER :: i, j, k
    INTEGER :: m, n
    LOGICAL :: flag

    INTEGER        :: iostat
    CHARACTER(256) :: iomsg


    NAMELIST / parallel_info / &
         ipes, jpes, kpes,     &
         input_method,         &
         output_method,        &
         sendrecv_method,      &
#ifdef MPIIO
         mpiio_zerofill,       &
#endif
         node_dim,             &
         remove_masked_pe

    CALL mpi_comm_size(MPI_COMM_WORLD, npes, ierr)

    CALL bcast(dimx)
    CALL bcast(dimy)
    CALL bcast(dimz)

    CALL bcast(open_e)
    CALL bcast(open_w)
    CALL bcast(open_n)
    CALL bcast(open_s)

    CALL bcast(open_alpha_e)
    CALL bcast(open_alpha_w)
    CALL bcast(open_alpha_n)
    CALL bcast(open_alpha_s)

    radi_e = open_e .AND. (open_alpha_e == UNDEF)
    radi_w = open_w .AND. (open_alpha_w == UNDEF)
    radi_n = open_n .AND. (open_alpha_n == UNDEF)
    radi_s = open_s .AND. (open_alpha_s == UNDEF)

    CALL bcast(cycle_x)
    CALL bcast(cycle_y)
    CALL bcast(cycle_z)

    CALL bcast(tripolar)

    CALL bcast(partial_step)
    CALL bcast(partial_step_threshold)
    CALL bcast(fill_chokepoint)

    IF (rank==0) THEN
       ipes = -1
       jpes = -1
       kpes =  1

#ifdef MPIIO
       input_method  = 'MPI-IO'
       output_method = 'MPI-IO'
#else
       input_method  = 'LAYERED'
       output_method = 'LAYERED'
#endif
       sendrecv_method = 'STD'

       REWIND(CONFIG_UNIT)
       READ(CONFIG_UNIT, NML=parallel_info, IOSTAT=iostat, IOMSG=iomsg)
       CALL assert(iostat <= 0, "failed to read PARALLEL_INFO namelist", iomsg)
    END IF

    CALL bcast(ipes)
    CALL bcast(jpes)
    CALL bcast(kpes)

    CALL bcast(input_method)
    CALL bcast(output_method)
    CALL bcast(sendrecv_method)

    CALL bcast(node_dim)

    CALL bcast(remove_masked_pe)
    CALL bcast(node_dim)

#ifdef MPIIO
    CALL bcast(mpiio_zerofill)
#endif


    IF (remove_masked_pe) THEN
       CALL assert((ipes >= 1 .AND. jpes >= 1 .AND. kpes >= 1), "IPES and JPES should be explicitly specified when REMOVE_MASKED_PE is enabled")
    ELSE
       CALL assert(mod(npes, kpes)==0, "KEPS should be a divisor of NPES")

       ALLOCATE(ratio(npes))

       IF (ipes < 1 .AND. jpes < 1) THEN
          ratio(:) = 0.0
          ipes = 0

          DO WHILE (ipes < npes/kpes)
             ipes = ipes +1

             IF (mod(npes/kpes, ipes) /= 0) CYCLE
             jpes = npes/kpes/ipes
             IF (mod(dimx, ipes) /= 0) CYCLE
             IF (mod(dimy, jpes) /= 0) CYCLE

             isize = dimx/ipes
             jsize = dimy/jpes

             ratio(ipes) = min(isize*1.0/jsize, jsize*1.0/isize)
          END DO
          ipes = maxloc(ratio, dim=1)
          jpes = npes/kpes/ipes
       ELSE IF (ipes < 1) THEN
          IF (rank==0) CALL assert(jpes <= npes, "JPES should not be greater than NPES")
          ipes = npes/kpes/jpes
       ELSE IF (jpes < 1) THEN
          IF (rank==0) CALL assert(ipes <= npes, "IPES should not be greater than NPES")
          jpes = npes/kpes/ipes
       END IF

       DEALLOCATE(ratio)
    END IF

    isize = dimx/ipes
    jsize = dimy/jpes
    ksize = dimz/kpes

    CALL assert(isize*ipes == dimx, "IPES should be a divisor of DIMX")
    CALL assert(jsize*jpes == dimy, "JPES should be a divisor of DIMY")
    CALL assert(ksize*kpes == dimz, "KPES should be a divisor of DIMZ")

    IF (maxval(node_dim) > 1) THEN
       CALL assert(mod(ipes, node_dim(1))==0, "NODE_DIM(1) should be a diviser of IPES")
       CALL assert(mod(jpes, node_dim(2))==0, "NODE_DIM(2) should be a diviser of JPES")
       CALL assert(mod(kpes, node_dim(3))==0, "NODE_DIM(3) should be a diviser of KPES")

       spes  = product(node_dim(:))
       srank = mod(rank, spes)
    END IF

    ALLOCATE(icoords(0:npes-1))
    ALLOCATE(jcoords(0:npes-1))
    ALLOCATE(kcoords(0:npes-1))

    IF (rank==0) THEN
       WRITE(REPORT_UNIT, *) "MPI PRALLEL JOB"
       WRITE(REPORT_UNIT, '(A,I6,A,I4,A,I4,A,I4,A,I4,A,I4,A,I4,A)') &
            "Total ", npes, " PEs, dimension (", ipes, ",", jpes,",", kpes, "),  ", isize, " x",jsize, " x",  ksize, " grids for each PE."
       WRITE(REPORT_UNIT, *) "     input method: ", trim(input_method)
       WRITE(REPORT_UNIT, *) "    output method: ", trim(output_method)
       WRITE(REPORT_UNIT, *) "  sendrecv method: ", trim(sendrecv_method)

    END IF

    IF (.NOT. remove_masked_pe) THEN
       CALL assert(kpes <= npes,              "KPES should be less than or equal to NPES")
       CALL assert(mod(npes,      kpes) == 0, "KPES should be a divisor of NPES")
       CALL assert(mod(npes/kpes, ipes) == 0, "IPES should be a divisor of NPES/KPES")
       CALL assert(mod(npes/kpes, jpes) == 0, "JPES should be a divisor of NPES/KPES")
       CALL assert(ipes*jpes*kpes == npes,    "IPES*JPES*KPES is not equal to NPES")

       IF (spes > 1) THEN
          CALL cart_coords((/ ipes/node_dim(1), jpes/node_dim(2), kpes/node_dim(3) /), rank/spes, coords)
          icoord = coords(1)*node_dim(1)
          jcoord = coords(2)*node_dim(2)
          kcoord = coords(3)*node_dim(3)

          CALL cart_coords(node_dim, srank, coords)
          icoord = icoord + coords(1)
          jcoord = jcoord + coords(2)
          kcoord = kcoord + coords(3)
       ELSE
          CALL cart_coords((/ipes, jpes, kpes/), rank, coords)
          icoord = coords(1)
          jcoord = coords(2)
          kcoord = coords(3)
       END IF

       CALL mpi_allgather(icoord, 1, MPI_INTEGER, icoords, 1, MPI_INTEGER, comm, ierr)
       CALL mpi_allgather(jcoord, 1, MPI_INTEGER, jcoords, 1, MPI_INTEGER, comm, ierr)
       CALL mpi_allgather(kcoord, 1, MPI_INTEGER, kcoords, 1, MPI_INTEGER, comm, ierr)

    ELSE
       IF (rank==0) THEN
          ALLOCATE(topo(1:dimx, 1:dimy))
          ALLOCATE(depth(0:dimz))

          topo(:,:) = 1.0*dimz
          topo(1:dimx,1:dimy) = 0.0
          IF (bathymetry /= '') THEN
             CALL read_file(topo(1:dimx,1:dimy), trim(bathymetry), KIND=kind)
          END IF

          topomin = minval(topo(1:dimx,1:dimy))

          DO k=0, dimz-1
             depth(k) = sum(delta_z(k+1:dimz))
          END DO
          depth(dimz) = 0.0

          flag = gridwise .OR. (topomin >= 0.0)

          n = 0
          DO i = ipes/node_dim(1)-1, 0, -1
          DO j = jpes/node_dim(2)-1, 0, -1
          DO k = kpes/node_dim(3)-1, 0, -1
             topomin = minval(topo(i*isize*node_dim(1)+1:(i+1)*isize*node_dim(1), j*jsize*node_dim(1)+1:(j+1)*jsize*node_dim(2)))
             IF (flag) THEN
                IF (topomin >= (k+1)*ksize*node_dim(3))         CYCLE
             ELSE
                IF (topomin >= -depth((k+1)*ksize*node_dim(3))) CYCLE
             END IF

             CALL assert(n*spes < npes, "more PEs are required!")

             DO m=0, spes-1
                CALL cart_coords(node_dim, m, coords)

                icoords(n+m) = i*node_dim(1) + coords(1)
                jcoords(n+m) = j*node_dim(2) + coords(2)
                kcoords(n+m) = k*node_dim(3) + coords(3)
             END DO
             n = n+spes
          END DO
          END DO
          END DO

          IF (n < npes)  WRITE(REPORT_UNIT, '(A,I0,A,I0,A)') "Only ", n, " PEs are required, ", npes-n, " unused process(es) will be killed."
          DEALLOCATE(topo)
          DEALLOCATE(depth)
       END IF

       CALL bcast(n)

       IF (rank < n) THEN
          CALL mpi_comm_split(MPI_COMM_WORLD, 0,             rank, comm, ierr)
       ELSE
          CALL mpi_comm_split(MPI_COMM_WORLD, MPI_UNDEFINED, rank, comm, ierr)
          CALL assert(comm==MPI_COMM_NULL, "something goes wrong in MPI_COMM_SPLIT???")
          CALL mpi_finalize(ierr)
          STOP
       END IF
       CALL mpi_barrier(comm, ierr)

       CALL mpi_comm_size(comm, npes, ierr)

       CALL bcast(icoords(0:npes-1))
       CALL bcast(jcoords(0:npes-1))
       CALL bcast(kcoords(0:npes-1))

       icoord = icoords(rank)
       jcoord = jcoords(rank)
       kcoord = kcoords(rank)

#ifndef PARALLEL3D
       CALL assert(kcoord==0, "KEPS should be 1")
#endif
    END IF

    IF (spes > 1) THEN
       CALL mpi_comm_split(MPI_COMM_WORLD, rank/spes,  rank, scomm, ierr)
       CALL mpi_comm_rank(scomm, itmp, ierr)
       CALL assert(srank==itmp, "invalid SRANK")
       CALL mpi_comm_size(scomm, itmp, ierr)
       CALL assert(spes==itmp, "invalid SPES")

       IF (vrank==0 .AND. srank==0) THEN
          CALL mpi_comm_split(comm, 0,             rank, gcomm, ierr)
          CALL mpi_comm_rank(gcomm, grank, ierr)
          CALL mpi_comm_size(gcomm, gpes,  ierr)
       ELSE
          CALL mpi_comm_split(comm, MPI_UNDEFINED, rank, gcomm, ierr)
          gpes  = 0
          grank = MPI_PROC_NULL
       END IF

    ELSE
       IF (vrank==0) THEN
          gcomm = hcomm
          gpes  = hpes
          grank = hrank
       ELSE
          gcomm = MPI_COMM_NULL
          gpes  = 0
          grank = MPI_PROC_NULL
       END IF
    END IF


    ALLOCATE(ranks(-1:ipes, -1:jpes, -1:kpes))
    ranks(:,:,:) = MPI_PROC_NULL

    DO n=0, npes-1
       ranks(icoords(n),jcoords(n),kcoords(n)) = n
    END DO

    IF (cycle_x) THEN
       ranks(-1,:,:)   = ranks(ipes-1,:,:)
       ranks(ipes,:,:) = ranks(0,:,:)
    END IF

    IF (cycle_y) THEN
       ranks(:,-1,:)   = ranks(:,jpes-1,:)
       ranks(:,jpes,:) = ranks(:,0,:)
    END IF

    IF (cycle_z) THEN
       ranks(:,:,-1)   = ranks(:,:,kpes-1)
       ranks(:,:,kpes) = ranks(:,:,0)
    END IF

    rank_e   = ranks(icoord+1, jcoord,   kcoord  )
    rank_w   = ranks(icoord-1, jcoord,   kcoord  )
    rank_n   = ranks(icoord,   jcoord+1, kcoord  )
    rank_s   = ranks(icoord,   jcoord-1, kcoord  )
    rank_u   = ranks(icoord,   jcoord,   kcoord+1)
    rank_l   = ranks(icoord,   jcoord,   kcoord-1)

    rank_ne  = ranks(icoord+1, jcoord+1, kcoord  )
    rank_nw  = ranks(icoord-1, jcoord+1, kcoord  )
    rank_se  = ranks(icoord+1, jcoord-1, kcoord  )
    rank_sw  = ranks(icoord-1, jcoord-1, kcoord  )
    rank_ue  = ranks(icoord+1, jcoord,   kcoord+1)
    rank_uw  = ranks(icoord-1, jcoord,   kcoord+1)
    rank_un  = ranks(icoord  , jcoord+1, kcoord+1)
    rank_us  = ranks(icoord  , jcoord-1, kcoord+1)
    rank_le  = ranks(icoord+1, jcoord,   kcoord-1)
    rank_lw  = ranks(icoord-1, jcoord,   kcoord-1)
    rank_ln  = ranks(icoord  , jcoord+1, kcoord-1)
    rank_ls  = ranks(icoord  , jcoord-1, kcoord-1)
    rank_une = ranks(icoord+1, jcoord+1, kcoord+1)
    rank_unw = ranks(icoord-1, jcoord+1, kcoord+1)
    rank_use = ranks(icoord+1, jcoord-1, kcoord+1)
    rank_usw = ranks(icoord-1, jcoord-1, kcoord+1)
    rank_lne = ranks(icoord+1, jcoord+1, kcoord-1)
    rank_lnw = ranks(icoord-1, jcoord+1, kcoord-1)
    rank_lse = ranks(icoord+1, jcoord-1, kcoord-1)
    rank_lsw = ranks(icoord-1, jcoord-1, kcoord-1)

    sendrank(tag_e)   = rank_e
    sendrank(tag_w)   = rank_w
    sendrank(tag_n)   = rank_n
    sendrank(tag_s)   = rank_s
    sendrank(tag_ne)  = rank_ne
    sendrank(tag_nw)  = rank_nw
    sendrank(tag_se)  = rank_se
    sendrank(tag_sw)  = rank_sw
#ifdef PARALLEL3D
    sendrank(tag_u)   = rank_u
    sendrank(tag_l)   = rank_l
    sendrank(tag_ue)  = rank_ue
    sendrank(tag_uw)  = rank_uw
    sendrank(tag_un)  = rank_un
    sendrank(tag_us)  = rank_us
    sendrank(tag_le)  = rank_le
    sendrank(tag_lw)  = rank_lw
    sendrank(tag_ln)  = rank_ln
    sendrank(tag_ls)  = rank_ls
    sendrank(tag_une) = rank_une
    sendrank(tag_unw) = rank_unw
    sendrank(tag_use) = rank_use
    sendrank(tag_usw) = rank_usw
    sendrank(tag_lne) = rank_lne
    sendrank(tag_lnw) = rank_lnw
    sendrank(tag_lse) = rank_lse
    sendrank(tag_lsw) = rank_lsw
#endif

    recvrank(tag_e)   = rank_w
    recvrank(tag_w)   = rank_e
    recvrank(tag_n)   = rank_s
    recvrank(tag_s)   = rank_n
    recvrank(tag_ne)  = rank_sw
    recvrank(tag_nw)  = rank_se
    recvrank(tag_se)  = rank_nw
    recvrank(tag_sw)  = rank_ne
#ifdef PARALLEL3D
    recvrank(tag_u)   = rank_l
    recvrank(tag_l)   = rank_u
    recvrank(tag_ue)  = rank_lw
    recvrank(tag_uw)  = rank_le
    recvrank(tag_un)  = rank_ls
    recvrank(tag_us)  = rank_ln
    recvrank(tag_le)  = rank_uw
    recvrank(tag_lw)  = rank_ue
    recvrank(tag_ln)  = rank_us
    recvrank(tag_ls)  = rank_un
    recvrank(tag_une) = rank_lsw
    recvrank(tag_unw) = rank_lse
    recvrank(tag_use) = rank_lnw
    recvrank(tag_usw) = rank_lne
    recvrank(tag_lne) = rank_usw
    recvrank(tag_lnw) = rank_use
    recvrank(tag_lse) = rank_unw
    recvrank(tag_lsw) = rank_une
#endif

    rank_tp(:,:) = MPI_PROC_NULL
    IF (tripolar .AND. jcoord==jpes-1) THEN
       CALL assert(mod(ipes, 2)==0, "TRIPOLAR coordinate requires IPES is an even number")

       rank_tp( 0, -1:1) = ranks(ipes-icoord-1, jcoord, kcoord-1:kcoord+1)
       rank_tp( 1, -1:1) = ranks(ipes-icoord-2, jcoord, kcoord-1:kcoord+1)
       rank_tp(-1, -1:1) = ranks(ipes-icoord,   jcoord, kcoord-1:kcoord+1)

       IF (icoord==ipes/2-1) rank_tp( 1,:) = MPI_PROC_NULL
       IF (icoord==ipes/2)   rank_tp(-1,:) = MPI_PROC_NULL
    END IF

#ifdef PARALLEL3D
    CALL mpi_comm_split(comm, jcoord*ipes + icoord, kpes-kcoord, vcomm, ierr)
    CALL mpi_comm_size(vcomm, vpes,  ierr)
    CALL mpi_comm_rank(vcomm, vrank, ierr)

    CALL mpi_comm_split(comm, kcoord, rank, hcomm, ierr)
    CALL mpi_comm_size(hcomm, hpes,  ierr)
    CALL mpi_comm_rank(hcomm, hrank, ierr)

    IF (vrank==0 .AND. hrank==0) CALL assert(rank==0, "unexpected root-rank assignment")
#else
    hcomm = comm
    hrank = rank
    hpes  = npes
    vcomm = MPI_COMM_NULL
    vpes  = 1
    vrank = 0
#endif

    ALLOCATE(icoords_h(0:hpes-1))
    ALLOCATE(jcoords_h(0:hpes-1))
    ALLOCATE(kcoords_v(0:vpes-1))

    CALL mpi_allgather(icoord, 1, MPI_INTEGER, icoords_h, 1, MPI_INTEGER, hcomm, ierr)
    CALL mpi_allgather(jcoord, 1, MPI_INTEGER, jcoords_h, 1, MPI_INTEGER, hcomm, ierr)

    ALLOCATE(hranks(-1:ipes, -1:jpes))

    hranks(:,:) = MPI_PROC_NULL
    DO n=0, hpes-1
       hranks(icoords_h(n), jcoords_h(n)) = n
    END DO

    IF (cycle_x) THEN
       hranks(-1,:)   = hranks(ipes-1,:)
       hranks(ipes,:) = hranks(0,:)
    END IF

    IF (cycle_y) THEN
       hranks(:,-1)   = hranks(:,jpes-1)
       hranks(:,jpes) = hranks(:,0)
    END IF

    CALL mpi_allgather(kcoord, 1, MPI_INTEGER, kcoords_v, 1, MPI_INTEGER, vcomm, ierr)

    ALLOCATE(vranks(-1:kpes))
    vranks(:) = MPI_PROC_NULL
    DO n=0, vpes-1
       vranks(kcoords_v(n)) = n
    END DO
    IF (cycle_z) THEN
       vranks(-1)   = vranks(kpes-1)
       vranks(kpes) = vranks(0)
    END IF

#ifdef DEBUG
    IF (rank==0) WRITE(STDERR_UNIT, *) "rank  (icoord,jcoord,kcoord) : hrank vrank"
    DO n=0, npes-1
       CALL mpi_barrier(comm, ierr)
       IF (n==rank) WRITE(STDERR_UNIT, '(I5,X,"(",I4,",",I4,",",I4,")",X,":",I4,X,I4)') rank, icoord, jcoord, kcoord, hrank, vrank
    END DO
#endif

    ALLOCATE(sendbuf_r4_e(jsize*ksize*(slv+1)))
    ALLOCATE(sendbuf_r8_e(jsize*ksize*(slv+1)))
    ALLOCATE(recvbuf_r4_e(jsize*ksize*(slv+1)))
    ALLOCATE(recvbuf_r8_e(jsize*ksize*(slv+1)))

    ALLOCATE(sendbuf_r4_w(jsize*ksize*(slv+1)))
    ALLOCATE(sendbuf_r8_w(jsize*ksize*(slv+1)))
    ALLOCATE(recvbuf_r4_w(jsize*ksize*(slv+1)))
    ALLOCATE(recvbuf_r8_w(jsize*ksize*(slv+1)))

    ALLOCATE(sendbuf_r4_n((isize+2*slv+1)*ksize*(slv+1)))
    ALLOCATE(sendbuf_r8_n((isize+2*slv+1)*ksize*(slv+1)))
    ALLOCATE(recvbuf_r4_n((isize+2*slv+1)*ksize*(slv+1)))
    ALLOCATE(recvbuf_r8_n((isize+2*slv+1)*ksize*(slv+1)))

    ALLOCATE(sendbuf_r4_s((isize+2*slv+1)*ksize*(slv+1)))
    ALLOCATE(sendbuf_r8_s((isize+2*slv+1)*ksize*(slv+1)))
    ALLOCATE(recvbuf_r4_s((isize+2*slv+1)*ksize*(slv+1)))
    ALLOCATE(recvbuf_r8_s((isize+2*slv+1)*ksize*(slv+1)))

    ALLOCATE(sendbuf_r4_ne(ksize*(slv+1)*(slv+1)))
    ALLOCATE(sendbuf_r8_ne(ksize*(slv+1)*(slv+1)))
    ALLOCATE(recvbuf_r4_ne(ksize*(slv+1)*(slv+1)))
    ALLOCATE(recvbuf_r8_ne(ksize*(slv+1)*(slv+1)))

    ALLOCATE(sendbuf_r4_nw(ksize*(slv+1)*(slv+1)))
    ALLOCATE(sendbuf_r8_nw(ksize*(slv+1)*(slv+1)))
    ALLOCATE(recvbuf_r4_nw(ksize*(slv+1)*(slv+1)))
    ALLOCATE(recvbuf_r8_nw(ksize*(slv+1)*(slv+1)))

    ALLOCATE(sendbuf_r4_se(ksize*(slv+1)*(slv+1)))
    ALLOCATE(sendbuf_r8_se(ksize*(slv+1)*(slv+1)))
    ALLOCATE(recvbuf_r4_se(ksize*(slv+1)*(slv+1)))
    ALLOCATE(recvbuf_r8_se(ksize*(slv+1)*(slv+1)))

    ALLOCATE(sendbuf_r4_sw(ksize*(slv+1)*(slv+1)))
    ALLOCATE(sendbuf_r8_sw(ksize*(slv+1)*(slv+1)))
    ALLOCATE(recvbuf_r4_sw(ksize*(slv+1)*(slv+1)))
    ALLOCATE(recvbuf_r8_sw(ksize*(slv+1)*(slv+1)))

    sendbuf_r4_e(:)  = UNDEF
    sendbuf_r4_w(:)  = UNDEF
    sendbuf_r4_n(:)  = UNDEF
    sendbuf_r4_s(:)  = UNDEF
    sendbuf_r4_ne(:) = UNDEF
    sendbuf_r4_nw(:) = UNDEF
    sendbuf_r4_se(:) = UNDEF
    sendbuf_r4_sw(:) = UNDEF

    sendbuf_r8_e(:)  = UNDEF
    sendbuf_r8_w(:)  = UNDEF
    sendbuf_r8_n(:)  = UNDEF
    sendbuf_r8_s(:)  = UNDEF
    sendbuf_r8_ne(:) = UNDEF
    sendbuf_r8_nw(:) = UNDEF
    sendbuf_r8_se(:) = UNDEF
    sendbuf_r8_sw(:) = UNDEF

    recvbuf_r4_e(:)  = UNDEF
    recvbuf_r4_w(:)  = UNDEF
    recvbuf_r4_n(:)  = UNDEF
    recvbuf_r4_s(:)  = UNDEF
    recvbuf_r4_ne(:) = UNDEF
    recvbuf_r4_nw(:) = UNDEF
    recvbuf_r4_se(:) = UNDEF
    recvbuf_r4_sw(:) = UNDEF

    recvbuf_r8_e(:)  = UNDEF
    recvbuf_r8_w(:)  = UNDEF
    recvbuf_r8_n(:)  = UNDEF
    recvbuf_r8_s(:)  = UNDEF
    recvbuf_r8_ne(:) = UNDEF
    recvbuf_r8_nw(:) = UNDEF
    recvbuf_r8_se(:) = UNDEF
    recvbuf_r8_sw(:) = UNDEF


#ifdef PARALLEL3D
    ALLOCATE(sendbuf_r4_u((isize+2*slv+1)*(jsize+2*slv+1)*(slv+1)))
    ALLOCATE(sendbuf_r8_u((isize+2*slv+1)*(jsize+2*slv+1)*(slv+1)))
    ALLOCATE(recvbuf_r4_u((isize+2*slv+1)*(jsize+2*slv+1)*(slv+1)))
    ALLOCATE(recvbuf_r8_u((isize+2*slv+1)*(jsize+2*slv+1)*(slv+1)))

    ALLOCATE(sendbuf_r4_l((isize+2*slv+1)*(jsize+2*slv+1)*(slv+1)))
    ALLOCATE(sendbuf_r8_l((isize+2*slv+1)*(jsize+2*slv+1)*(slv+1)))
    ALLOCATE(recvbuf_r4_l((isize+2*slv+1)*(jsize+2*slv+1)*(slv+1)))
    ALLOCATE(recvbuf_r8_l((isize+2*slv+1)*(jsize+2*slv+1)*(slv+1)))

    ALLOCATE(sendbuf_r4_ue(jsize*(slv+1)*(slv+1)))
    ALLOCATE(sendbuf_r8_ue(jsize*(slv+1)*(slv+1)))
    ALLOCATE(recvbuf_r4_ue(jsize*(slv+1)*(slv+1)))
    ALLOCATE(recvbuf_r8_ue(jsize*(slv+1)*(slv+1)))

    ALLOCATE(sendbuf_r4_uw(jsize*(slv+1)*(slv+1)))
    ALLOCATE(sendbuf_r8_uw(jsize*(slv+1)*(slv+1)))
    ALLOCATE(recvbuf_r4_uw(jsize*(slv+1)*(slv+1)))
    ALLOCATE(recvbuf_r8_uw(jsize*(slv+1)*(slv+1)))

    ALLOCATE(sendbuf_r4_un(isize*(slv+1)*(slv+1)))
    ALLOCATE(sendbuf_r8_un(isize*(slv+1)*(slv+1)))
    ALLOCATE(recvbuf_r4_un(isize*(slv+1)*(slv+1)))
    ALLOCATE(recvbuf_r8_un(isize*(slv+1)*(slv+1)))

    ALLOCATE(sendbuf_r4_us(isize*(slv+1)*(slv+1)))
    ALLOCATE(sendbuf_r8_us(isize*(slv+1)*(slv+1)))
    ALLOCATE(recvbuf_r4_us(isize*(slv+1)*(slv+1)))
    ALLOCATE(recvbuf_r8_us(isize*(slv+1)*(slv+1)))

    ALLOCATE(sendbuf_r4_le(jsize*(slv+1)*(slv+1)))
    ALLOCATE(sendbuf_r8_le(jsize*(slv+1)*(slv+1)))
    ALLOCATE(recvbuf_r4_le(jsize*(slv+1)*(slv+1)))
    ALLOCATE(recvbuf_r8_le(jsize*(slv+1)*(slv+1)))

    ALLOCATE(sendbuf_r4_lw(jsize*(slv+1)*(slv+1)))
    ALLOCATE(sendbuf_r8_lw(jsize*(slv+1)*(slv+1)))
    ALLOCATE(recvbuf_r4_lw(jsize*(slv+1)*(slv+1)))
    ALLOCATE(recvbuf_r8_lw(jsize*(slv+1)*(slv+1)))

    ALLOCATE(sendbuf_r4_ln(isize*(slv+1)*(slv+1)))
    ALLOCATE(sendbuf_r8_ln(isize*(slv+1)*(slv+1)))
    ALLOCATE(recvbuf_r4_ln(isize*(slv+1)*(slv+1)))
    ALLOCATE(recvbuf_r8_ln(isize*(slv+1)*(slv+1)))

    ALLOCATE(sendbuf_r4_ls(isize*(slv+1)*(slv+1)))
    ALLOCATE(sendbuf_r8_ls(isize*(slv+1)*(slv+1)))
    ALLOCATE(recvbuf_r4_ls(isize*(slv+1)*(slv+1)))
    ALLOCATE(recvbuf_r8_ls(isize*(slv+1)*(slv+1)))

    ALLOCATE(sendbuf_r4_une((slv+1)*(slv+1)*(slv+1)))
    ALLOCATE(sendbuf_r8_une((slv+1)*(slv+1)*(slv+1)))
    ALLOCATE(recvbuf_r4_une((slv+1)*(slv+1)*(slv+1)))
    ALLOCATE(recvbuf_r8_une((slv+1)*(slv+1)*(slv+1)))

    ALLOCATE(sendbuf_r4_unw((slv+1)*(slv+1)*(slv+1)))
    ALLOCATE(sendbuf_r8_unw((slv+1)*(slv+1)*(slv+1)))
    ALLOCATE(recvbuf_r4_unw((slv+1)*(slv+1)*(slv+1)))
    ALLOCATE(recvbuf_r8_unw((slv+1)*(slv+1)*(slv+1)))

    ALLOCATE(sendbuf_r4_use((slv+1)*(slv+1)*(slv+1)))
    ALLOCATE(sendbuf_r8_use((slv+1)*(slv+1)*(slv+1)))
    ALLOCATE(recvbuf_r4_use((slv+1)*(slv+1)*(slv+1)))
    ALLOCATE(recvbuf_r8_use((slv+1)*(slv+1)*(slv+1)))

    ALLOCATE(sendbuf_r4_usw((slv+1)*(slv+1)*(slv+1)))
    ALLOCATE(sendbuf_r8_usw((slv+1)*(slv+1)*(slv+1)))
    ALLOCATE(recvbuf_r4_usw((slv+1)*(slv+1)*(slv+1)))
    ALLOCATE(recvbuf_r8_usw((slv+1)*(slv+1)*(slv+1)))

    ALLOCATE(sendbuf_r4_lne((slv+1)*(slv+1)*(slv+1)))
    ALLOCATE(sendbuf_r8_lne((slv+1)*(slv+1)*(slv+1)))
    ALLOCATE(recvbuf_r4_lne((slv+1)*(slv+1)*(slv+1)))
    ALLOCATE(recvbuf_r8_lne((slv+1)*(slv+1)*(slv+1)))

    ALLOCATE(sendbuf_r4_lnw((slv+1)*(slv+1)*(slv+1)))
    ALLOCATE(sendbuf_r8_lnw((slv+1)*(slv+1)*(slv+1)))
    ALLOCATE(recvbuf_r4_lnw((slv+1)*(slv+1)*(slv+1)))
    ALLOCATE(recvbuf_r8_lnw((slv+1)*(slv+1)*(slv+1)))

    ALLOCATE(sendbuf_r4_lse((slv+1)*(slv+1)*(slv+1)))
    ALLOCATE(sendbuf_r8_lse((slv+1)*(slv+1)*(slv+1)))
    ALLOCATE(recvbuf_r4_lse((slv+1)*(slv+1)*(slv+1)))
    ALLOCATE(recvbuf_r8_lse((slv+1)*(slv+1)*(slv+1)))

    ALLOCATE(sendbuf_r4_lsw((slv+1)*(slv+1)*(slv+1)))
    ALLOCATE(sendbuf_r8_lsw((slv+1)*(slv+1)*(slv+1)))
    ALLOCATE(recvbuf_r4_lsw((slv+1)*(slv+1)*(slv+1)))
    ALLOCATE(recvbuf_r8_lsw((slv+1)*(slv+1)*(slv+1)))

    sendbuf_r8_u(:)   = UNDEF
    sendbuf_r8_l(:)   = UNDEF
    sendbuf_r8_ue(:)  = UNDEF
    sendbuf_r8_uw(:)  = UNDEF
    sendbuf_r8_le(:)  = UNDEF
    sendbuf_r8_lw(:)  = UNDEF
    sendbuf_r8_une(:) = UNDEF
    sendbuf_r8_unw(:) = UNDEF
    sendbuf_r8_use(:) = UNDEF
    sendbuf_r8_usw(:) = UNDEF
    sendbuf_r8_lne(:) = UNDEF
    sendbuf_r8_lnw(:) = UNDEF
    sendbuf_r8_lse(:) = UNDEF
    sendbuf_r8_lsw(:) = UNDEF

    recvbuf_r8_u(:)   = UNDEF
    recvbuf_r8_l(:)   = UNDEF
    recvbuf_r8_ue(:)  = UNDEF
    recvbuf_r8_uw(:)  = UNDEF
    recvbuf_r8_le(:)  = UNDEF
    recvbuf_r8_lw(:)  = UNDEF
    recvbuf_r8_une(:) = UNDEF
    recvbuf_r8_unw(:) = UNDEF
    recvbuf_r8_use(:) = UNDEF
    recvbuf_r8_usw(:) = UNDEF
    recvbuf_r8_lne(:) = UNDEF
    recvbuf_r8_lnw(:) = UNDEF
    recvbuf_r8_lse(:) = UNDEF
    recvbuf_r8_lsw(:) = UNDEF
#endif

    mpireq_2dx_r4(:,:) = MPI_REQUEST_NULL
    mpireq_2dx_r8(:,:) = MPI_REQUEST_NULL
    mpireq_2dy_r4(:,:) = MPI_REQUEST_NULL
    mpireq_2dy_r8(:,:) = MPI_REQUEST_NULL
    mpireq_3dx_r4(:,:) = MPI_REQUEST_NULL
    mpireq_3dx_r8(:,:) = MPI_REQUEST_NULL
    mpireq_3dy_r4(:,:) = MPI_REQUEST_NULL
    mpireq_3dy_r8(:,:) = MPI_REQUEST_NULL
    mpireq_3dz_r4(:,:) = MPI_REQUEST_NULL
    mpireq_3dz_r8(:,:) = MPI_REQUEST_NULL

    DO n=1, 5
       slv0 = (n+1)/2
       slv1 = n/2

       IF (slv0 > 0) THEN
          CALL mpi_send_init(sendbuf_r4_e, slv0*jsize,               MPI_REAL4, rank_e, 1, comm, mpireq_2dx_r4(1,n), ierr)
          CALL mpi_send_init(sendbuf_r8_e, slv0*jsize,               MPI_REAL8, rank_e, 1, comm, mpireq_2dx_r8(1,n), ierr)

          CALL mpi_recv_init(recvbuf_r4_w, slv0*jsize,               MPI_REAL4, rank_w, 1, comm, mpireq_2dx_r4(2,n), ierr)
          CALL mpi_recv_init(recvbuf_r8_w, slv0*jsize,               MPI_REAL8, rank_w, 1, comm, mpireq_2dx_r8(2,n), ierr)

          CALL mpi_send_init(sendbuf_r4_n, slv0*(isize+5),           MPI_REAL4, rank_n, 3, comm, mpireq_2dy_r4(1,n), ierr)
          CALL mpi_send_init(sendbuf_r8_n, slv0*(isize+5),           MPI_REAL8, rank_n, 3, comm, mpireq_2dy_r8(1,n), ierr)

          CALL mpi_recv_init(recvbuf_r4_s, slv0*(isize+5),           MPI_REAL4, rank_s, 3, comm, mpireq_2dy_r4(2,n), ierr)
          CALL mpi_recv_init(recvbuf_r8_s, slv0*(isize+5),           MPI_REAL8, rank_s, 3, comm, mpireq_2dy_r8(2,n), ierr)

          CALL mpi_send_init(sendbuf_r4_e, slv0*jsize*ksize,         MPI_REAL4, rank_e, 1, comm, mpireq_3dx_r4(1,n), ierr)
          CALL mpi_send_init(sendbuf_r8_e, slv0*jsize*ksize,         MPI_REAL8, rank_e, 1, comm, mpireq_3dx_r8(1,n), ierr)

          CALL mpi_recv_init(recvbuf_r4_w, slv0*jsize*ksize,         MPI_REAL4, rank_w, 1, comm, mpireq_3dx_r4(2,n), ierr)
          CALL mpi_recv_init(recvbuf_r8_w, slv0*jsize*ksize,         MPI_REAL8, rank_w, 1, comm, mpireq_3dx_r8(2,n), ierr)

          CALL mpi_send_init(sendbuf_r4_n, slv0*(isize+5)*ksize,     MPI_REAL4, rank_n, 3, comm, mpireq_3dy_r4(1,n), ierr)
          CALL mpi_send_init(sendbuf_r8_n, slv0*(isize+5)*ksize,     MPI_REAL8, rank_n, 3, comm, mpireq_3dy_r8(1,n), ierr)

          CALL mpi_recv_init(recvbuf_r4_s, slv0*(isize+5)*ksize,     MPI_REAL4, rank_s, 3, comm, mpireq_3dy_r4(2,n), ierr)
          CALL mpi_recv_init(recvbuf_r8_s, slv0*(isize+5)*ksize,     MPI_REAL8, rank_s, 3, comm, mpireq_3dy_r8(2,n), ierr)

#ifdef PARALLEL3D
          CALL mpi_send_init(sendbuf_r4_u, slv0*(isize+5)*(jsize+5), MPI_REAL4, rank_u, 5, comm, mpireq_3dz_r4(1,n), ierr)
          CALL mpi_send_init(sendbuf_r8_u, slv0*(isize+5)*(jsize+5), MPI_REAL8, rank_u, 5, comm, mpireq_3dz_r8(1,n), ierr)

          CALL mpi_recv_init(recvbuf_r4_l, slv0*(isize+5)*(jsize+5), MPI_REAL4, rank_l, 5, comm, mpireq_3dz_r4(2,n), ierr)
          CALL mpi_recv_init(recvbuf_r8_l, slv0*(isize+5)*(jsize+5), MPI_REAL8, rank_l, 5, comm, mpireq_3dz_r8(2,n), ierr)
#endif
       END IF

       IF (slv1 > 0) THEN
          CALL mpi_send_init(sendbuf_r4_w, slv1*jsize,               MPI_REAL4, rank_w, 0, comm, mpireq_2dx_r4(3,n), ierr)
          CALL mpi_send_init(sendbuf_r8_w, slv1*jsize,               MPI_REAL8, rank_w, 0, comm, mpireq_2dx_r8(3,n), ierr)

          CALL mpi_recv_init(recvbuf_r4_e, slv1*jsize,               MPI_REAL4, rank_e, 0, comm, mpireq_2dx_r4(4,n), ierr)
          CALL mpi_recv_init(recvbuf_r8_e, slv1*jsize,               MPI_REAL8, rank_e, 0, comm, mpireq_2dx_r8(4,n), ierr)

          CALL mpi_send_init(sendbuf_r4_s, slv1*(isize+5),           MPI_REAL4, rank_s, 2, comm, mpireq_2dy_r4(3,n), ierr)
          CALL mpi_send_init(sendbuf_r8_s, slv1*(isize+5),           MPI_REAL8, rank_s, 2, comm, mpireq_2dy_r8(3,n), ierr)

          CALL mpi_recv_init(recvbuf_r4_n, slv1*(isize+5),           MPI_REAL4, rank_n, 2, comm, mpireq_2dy_r4(4,n), ierr)
          CALL mpi_recv_init(recvbuf_r8_n, slv1*(isize+5),           MPI_REAL8, rank_n, 2, comm, mpireq_2dy_r8(4,n), ierr)

          CALL mpi_send_init(sendbuf_r4_w, slv1*jsize*ksize,         MPI_REAL4, rank_w, 0, comm, mpireq_3dx_r4(3,n), ierr)
          CALL mpi_send_init(sendbuf_r8_w, slv1*jsize*ksize,         MPI_REAL8, rank_w, 0, comm, mpireq_3dx_r8(3,n), ierr)

          CALL mpi_recv_init(recvbuf_r4_e, slv1*jsize*ksize,         MPI_REAL4, rank_e, 0, comm, mpireq_3dx_r4(4,n), ierr)
          CALL mpi_recv_init(recvbuf_r8_e, slv1*jsize*ksize,         MPI_REAL8, rank_e, 0, comm, mpireq_3dx_r8(4,n), ierr)

          CALL mpi_send_init(sendbuf_r4_s, slv1*(isize+5)*ksize,     MPI_REAL4, rank_s, 2, comm, mpireq_3dy_r4(3,n), ierr)
          CALL mpi_send_init(sendbuf_r8_s, slv1*(isize+5)*ksize,     MPI_REAL8, rank_s, 2, comm, mpireq_3dy_r8(3,n), ierr)

          CALL mpi_recv_init(recvbuf_r4_n, slv1*(isize+5)*ksize,     MPI_REAL4, rank_n, 2, comm, mpireq_3dy_r4(4,n), ierr)
          CALL mpi_recv_init(recvbuf_r8_n, slv1*(isize+5)*ksize,     MPI_REAL8, rank_n, 2, comm, mpireq_3dy_r8(4,n), ierr)

#ifdef PARALLEL3D
          CALL mpi_send_init(sendbuf_r4_l, slv1*(isize+5)*(jsize+5), MPI_REAL4, rank_l, 4, comm, mpireq_3dz_r4(3,n), ierr)
          CALL mpi_send_init(sendbuf_r8_l, slv1*(isize+5)*(jsize+5), MPI_REAL8, rank_l, 4, comm, mpireq_3dz_r8(3,n), ierr)

          CALL mpi_recv_init(recvbuf_r4_u, slv1*(isize+5)*(jsize+5), MPI_REAL4, rank_u, 4, comm, mpireq_3dz_r4(4,n), ierr)
          CALL mpi_recv_init(recvbuf_r8_u, slv1*(isize+5)*(jsize+5), MPI_REAL8, rank_u, 4, comm, mpireq_3dz_r8(4,n), ierr)
#endif
       END IF

    END DO

#ifdef MPIIO
    file%handle = MPI_FILE_NULL
    ALLOCATE(file%buffer(isize*jsize*ksize*8))
    file%buffer(:) = 0_1
    
#ifdef MPIIO_ASYNCHRONOUS
    DO n=1, n_async
       file_async(n)%handle = MPI_FILE_NULL
       ALLOCATE(file_async(n)%buffer(isize*jsize*ksize*8))
       file_async(n)%buffer(:) = 0_1
    END DO
#endif

    CALL init_fileview(0, 0_MPI_OFFSET_KIND, dimx, dimy, 1, 1)

    mpiio_zerofill = mpiio_zerofill .AND. remove_masked_pe
#endif
  CONTAINS
    SUBROUTINE cart_coords(sizes, rank, coords)
    !row-major, reverse-order
      INTEGER, INTENT(IN)  :: sizes(3), rank
      INTEGER, INTENT(OUT) :: coords(3)

      coords(1) = rank / (sizes(2)*sizes(3))
      coords(2) = (rank - coords(1)*sizes(2)*sizes(3)) / sizes(3)
      coords(3) = rank - coords(1)*sizes(2)*sizes(3) - coords(2)*sizes(3)
      coords(1) = sizes(1) - coords(1) - 1
      coords(2) = sizes(2) - coords(2) - 1
      coords(3) = sizes(3) - coords(3) - 1

    END SUBROUTINE cart_coords

  END SUBROUTINE init_parallel_geometry


#ifdef MPIIO
  SUBROUTINE init_fileview(id, offset, x_size, y_size, x_start, y_start)
    INTEGER, INTENT(IN) :: id
    INTEGER, INTENT(IN) :: x_size
    INTEGER, INTENT(IN) :: y_size
    INTEGER, INTENT(IN) :: x_start
    INTEGER, INTENT(IN) :: y_start
    INTEGER(MPI_OFFSET_KIND), INTENT(IN) :: offset

    INTEGER :: ierr

    CALL assert(.NOT. views(id)%defined, "FIELVIEW ID="//format(id)//" is already defined")

    CALL mpi_type_create_subarray(2, (/x_size*1, y_size/), (/isize*1, jsize/),               &
                                     (/(x_start+icoord*isize-1)*1, y_start+jcoord*jsize-1/), &
                                     MPI_ORDER_FORTRAN, MPI_BYTE, views(id)%subarray_2d_i1, ierr)

    CALL mpi_type_create_subarray(2, (/x_size*4, y_size/), (/isize*4, jsize/),               &
                                     (/(x_start+icoord*isize-1)*4, y_start+jcoord*jsize-1/), &
                                     MPI_ORDER_FORTRAN, MPI_BYTE, views(id)%subarray_2d_r4, ierr)

    CALL mpi_type_create_subarray(2, (/x_size*8, y_size/), (/isize*8, jsize/),               &
                                     (/(x_start+icoord*isize-1)*8, y_start+jcoord*jsize-1/), &
                                     MPI_ORDER_FORTRAN, MPI_BYTE, views(id)%subarray_2d_r8, ierr)

    CALL mpi_type_create_subarray(3, (/x_size*1, y_size, dimz/), (/isize*1, jsize, ksize/),                &
                                     (/(x_start+icoord*isize-1)*1, y_start+jcoord*jsize-1, kcoord*ksize/), &
                                     MPI_ORDER_FORTRAN, MPI_BYTE, views(id)%subarray_3d_i1, ierr)

    CALL mpi_type_create_subarray(3, (/x_size*4, y_size, dimz/), (/isize*4, jsize, ksize/),                &
                                     (/(x_start+icoord*isize-1)*4, y_start+jcoord*jsize-1, kcoord*ksize/), &
                                     MPI_ORDER_FORTRAN, MPI_BYTE, views(id)%subarray_3d_r4, ierr)

    CALL mpi_type_create_subarray(3, (/x_size*8, y_size, dimz/), (/isize*8, jsize, ksize/),                &
                                     (/(x_start+icoord*isize-1)*8, y_start+jcoord*jsize-1, kcoord*ksize/), &
                                     MPI_ORDER_FORTRAN, MPI_BYTE, views(id)%subarray_3d_r8, ierr)

    CALL mpi_type_create_subarray(3, (/x_size*1, y_size, dimz/), (/isize*1, jsize, ksize/),                         &
                                     (/(x_start+icoord*isize-1)*1, y_start+jcoord*jsize-1, (kpes-kcoord-1)*ksize/), &
                                     MPI_ORDER_FORTRAN, MPI_BYTE, views(id)%subarray_3d_i1_desc, ierr)

    CALL mpi_type_create_subarray(3, (/x_size*4, y_size, dimz/), (/isize*4, jsize, ksize/),                         &
                                     (/(x_start+icoord*isize-1)*4, y_start+jcoord*jsize-1, (kpes-kcoord-1)*ksize/), &
                                     MPI_ORDER_FORTRAN, MPI_BYTE, views(id)%subarray_3d_r4_desc, ierr)

    CALL mpi_type_create_subarray(3, (/x_size*8, y_size, dimz/), (/isize*8, jsize, ksize/),                         &
                                     (/(x_start+icoord*isize-1)*8, y_start+jcoord*jsize-1, (kpes-kcoord-1)*ksize/), &
                                     MPI_ORDER_FORTRAN, MPI_BYTE, views(id)%subarray_3d_r8_desc, ierr)

    CALL mpi_type_commit(views(id)%subarray_3d_i1, ierr)
    CALL mpi_type_commit(views(id)%subarray_3d_r4, ierr)
    CALL mpi_type_commit(views(id)%subarray_3d_r8, ierr)
    CALL mpi_type_commit(views(id)%subarray_2d_i1, ierr)
    CALL mpi_type_commit(views(id)%subarray_2d_r4, ierr)
    CALL mpi_type_commit(views(id)%subarray_2d_r8, ierr)
    CALL mpi_type_commit(views(id)%subarray_3d_i1_desc, ierr)
    CALL mpi_type_commit(views(id)%subarray_3d_r4_desc, ierr)
    CALL mpi_type_commit(views(id)%subarray_3d_r8_desc, ierr)

    views(id)%offset = offset

    views(id)%defined = .TRUE.

  END SUBROUTINE init_fileview
#endif
#endif

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE finalize_geometry
    INTEGER :: n
    INTEGER :: ierr

#ifdef MPIIO_ASYNCHRONOUS
    DO n=1, n_async
       IF (file_async(n)%handle /= MPI_FILE_NULL) THEN
          CALL mpi_file_write_all_end(file_async(n)%handle, file_async(n)%buffer, MPI_STATUS_IGNORE, ierr)
          CALL mpi_file_close(file_async(n)%handle, ierr)
          DEALLOCATE(file_async(n)%buffer)
       END IF
    END DO
#endif
    DEALLOCATE(imask2d)
    DEALLOCATE(imask3d)
    DEALLOCATE(lmask2d)
    DEALLOCATE(lmask3d)

    DEALLOCATE(dx)
    DEALLOCATE(dy)
    DEALLOCATE(dz)
    DEALLOCATE(dz_ref)
    DEALLOCATE(dz_star)

    DEALLOCATE(dzindex)

    DEALLOCATE(idx0)
    DEALLOCATE(idy0)
    DEALLOCATE(idz0)

    DEALLOCATE(idx1)
    DEALLOCATE(idy1)
    DEALLOCATE(idz1)

    DEALLOCATE(idx2)
    DEALLOCATE(idy2)
    DEALLOCATE(idz2)

    DEALLOCATE(dvol)
    DEALLOCATE(dvol_old)

    DEALLOCATE(dsx)
    DEALLOCATE(dsy)
    DEALLOCATE(dsz)

    DEALLOCATE(dsx_old)
    DEALLOCATE(dsy_old)
    DEALLOCATE(dsx_ref)
    DEALLOCATE(dsy_ref)

    DEALLOCATE(dsx2d)
    DEALLOCATE(dsy2d)

    DEALLOCATE(metxy)
    DEALLOCATE(metyx)

    DEALLOCATE(corx)
    DEALLOCATE(cory)
    DEALLOCATE(corz)

    DEALLOCATE(depth)

    DEALLOCATE(delta2)

    DEALLOCATE(ssh)

#ifdef PARALLEL_MPI
    DEALLOCATE(sendbuf_r8_e)
    DEALLOCATE(sendbuf_r8_w)
    DEALLOCATE(sendbuf_r8_n)
    DEALLOCATE(sendbuf_r8_s)
    DEALLOCATE(sendbuf_r8_ne)
    DEALLOCATE(sendbuf_r8_nw)
    DEALLOCATE(sendbuf_r8_se)
    DEALLOCATE(sendbuf_r8_sw)

    DEALLOCATE(recvbuf_r8_e)
    DEALLOCATE(recvbuf_r8_w)
    DEALLOCATE(recvbuf_r8_n)
    DEALLOCATE(recvbuf_r8_s)
    DEALLOCATE(recvbuf_r8_ne)
    DEALLOCATE(recvbuf_r8_nw)
    DEALLOCATE(recvbuf_r8_se)
    DEALLOCATE(recvbuf_r8_sw)

#ifdef PARALLEL3D
    DEALLOCATE(sendbuf_r8_u)
    DEALLOCATE(sendbuf_r8_l)
    DEALLOCATE(sendbuf_r8_ue)
    DEALLOCATE(sendbuf_r8_uw)
    DEALLOCATE(sendbuf_r8_un)
    DEALLOCATE(sendbuf_r8_us)
    DEALLOCATE(sendbuf_r8_le)
    DEALLOCATE(sendbuf_r8_lw)
    DEALLOCATE(sendbuf_r8_ln)
    DEALLOCATE(sendbuf_r8_ls)
    DEALLOCATE(sendbuf_r8_une)
    DEALLOCATE(sendbuf_r8_unw)
    DEALLOCATE(sendbuf_r8_use)
    DEALLOCATE(sendbuf_r8_usw)
    DEALLOCATE(sendbuf_r8_lne)
    DEALLOCATE(sendbuf_r8_lnw)
    DEALLOCATE(sendbuf_r8_lse)
    DEALLOCATE(sendbuf_r8_lsw)

    DEALLOCATE(recvbuf_r8_u)
    DEALLOCATE(recvbuf_r8_l)
    DEALLOCATE(recvbuf_r8_ue)
    DEALLOCATE(recvbuf_r8_uw)
    DEALLOCATE(recvbuf_r8_un)
    DEALLOCATE(recvbuf_r8_us)
    DEALLOCATE(recvbuf_r8_le)
    DEALLOCATE(recvbuf_r8_lw)
    DEALLOCATE(recvbuf_r8_ln)
    DEALLOCATE(recvbuf_r8_ls)
    DEALLOCATE(recvbuf_r8_une)
    DEALLOCATE(recvbuf_r8_unw)
    DEALLOCATE(recvbuf_r8_use)
    DEALLOCATE(recvbuf_r8_usw)
    DEALLOCATE(recvbuf_r8_lne)
    DEALLOCATE(recvbuf_r8_lnw)
    DEALLOCATE(recvbuf_r8_lse)
    DEALLOCATE(recvbuf_r8_lsw)
#endif
#endif

  END SUBROUTINE finalize_geometry

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE update_geometry
#ifdef PARALLEL_MPI
    REAL(8) :: sendbuf((isize+2)*(jsize+2)*2)
    INTEGER :: ierr
#endif

    INTEGER :: i, j, k

!$OMP PARALLEL
!$OMP DO
    DO j=1-slv, jsize+slv
    DO i=1-slv, isize+slv
       IF (.NOT. lmask2d(i,j)) ssh(i,j) = 0.0D0
    END DO
    END DO

!$OMP DO
    DO k=1-slv, ksize+slv
    DO j=1-slv, jsize+slv
    DO i=1-slv, isize+slv
       dvol_old(i,j,k) = dvol(i,j,k)
       dz_star(i,j,k) = dz_ref(i,j,k) + ssh(i,j) * cz_star(i,j,k)
       dvol(i,j,k) = dx(i,j)*dy(i,j)*dz_star(i,j,k)
    END DO
    END DO
    END DO
!$OMP END PARALLEL

    total_vol = sum(dvol(1:isize, 1:jsize, 1:ksize)*imask3d(1:isize,1:jsize,1:ksize))
    CALL gsum(total_vol, all=.TRUE.)

    IF (use_landwater .AND. vrank==0) THEN
       DO j=1-slv, jsize+slv
       DO i=1-slv, isize+slv
          lwdried(i,j) = lwflag(i,j) .AND. (dz_star(i,j,ksize) < landwater_hdry)
          dvol(i,j,ksize)  = max(dvol(i,j,ksize), dx(i,j)*dy(i,j) * landwater_hdry)
       END DO
       END DO
    END IF

!$OMP PARALLEL
!$OMP DO
    DO k=1-slv, ksize+slv
       DO j=1-slv, jsize+slv
       DO i=1-slv, isize+slv-1
          dsx_old(i,j,k) = dsx(i,j,k)
!         dsx(i,j,k) = (dy(i,  j)*dz_star(i,  j,k)*dx(i+1,j) &
!                     + dy(i+1,j)*dz_star(i+1,j,k)*dx(i,  j)) / (dx(i,j)+dx(i+1,j))
          dsx(i,j,k) = 0.5*(dvol(i,j,k)+dvol(i+1,j,k))*idx1(i,j)
       END DO
       END DO

       DO j=1-slv, jsize+slv-1
       DO i=1-slv, isize+slv
          dsy_old(i,j,k) = dsy(i,j,k)
!         dsy(i,j,k) = (dx(i,j  )*dz_star(i,j,  k)*dy(i,j+1) &
!                     + dx(i,j+1)*dz_star(i,j+1,k)*dy(i,j  )) / (dy(i,j)+dy(i,j+1))
          dsy(i,j,k) = 0.5*(dvol(i,j,k)+dvol(i,j+1,k))*idy1(i,j)
       END DO
       END DO
    END DO

!$OMP DO
    DO j=1-slv, jsize+slv
    DO i=1-slv, isize+slv-1
       dsx2d(i,j) = sum(dsx(i,j,1:ksize)*imask3d(i,j,1:ksize)*imask3d(i+1,j,1:ksize))
    END DO
    END DO

!$OMP DO
    DO j=1-slv, jsize+slv-1
    DO i=1-slv, isize+slv
       dsy2d(i,j) = sum(dsy(i,j,1:ksize)*imask3d(i,j,1:ksize)*imask3d(i,j+1,1:ksize))
    END DO
    END DO

!$OMP END PARALLEL

    IF (open_w .AND. icoord==0) THEN
       DO j=1-slv, jsize+slv
          dsx2d(0,j) = sum(dsx(0,j,1:ksize)*imask3d(1,j,1:ksize))
       END DO
    END IF

    IF (open_e .AND. icoord==ipes-1) THEN
       DO j=1-slv, jsize+slv
          dsx2d(isize,j) = sum(dsx(isize,j,1:ksize)*imask3d(isize,j,1:ksize))
       END DO
    END IF

    IF (open_s .AND. jcoord==0) THEN
       DO i=1-slv, isize+slv
          dsy2d(i,0) = sum(dsy(i,0,1:ksize)*imask3d(i,1,1:ksize))
       END DO
    END IF

    IF (open_n .AND. jcoord==jpes-1) THEN
       DO i=1-slv, isize+slv
          dsy2d(i,jsize) = sum(dsy(i,jsize,1:ksize)*imask3d(i,jsize,1:ksize))
       END DO
    END IF

    CALL vsum(dsx2d, all=.TRUE.)
    CALL vsum(dsy2d, all=.TRUE.)

    DO j=1-slv, jsize+slv
    DO i=1-slv, isize+slv
       cext(i,j) = imask2d(i,j)*sqrt((wct_ref(i,j)+ssh(i,j))*gravity)
    END DO
    END DO

  END SUBROUTINE update_geometry

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE update_boundary_2d_r4(a, fill, method)
    REAL(4),      INTENT(INOUT)        :: a(:,:)
    REAL(4),      INTENT(IN), OPTIONAL :: fill
    CHARACTER(*), INTENT(IN), OPTIONAL :: method
#ifdef F2008
    CONTIGUOUS a
#endif

    INTEGER :: nx, ny
    INTEGER :: slv_e, slv_w, slv_n, slv_s

    CHARACTER(16) :: method_

    nx = size(a,1) - isize
    ny = size(a,2) - jsize

    CALL assert(nx >= 0 .AND. nx <= 5, "unsported dimension in UPDATE_BOUNDARY_2D")
    CALL assert(ny >= 0 .AND. ny <= 5, "unsported dimension in UPDATE_BOUNDARY_2D")

    slv_e = nx/2
    slv_n = ny/2
    slv_w = (nx+1)/2
    slv_s = (ny+1)/2

#ifdef PARALLEL_MPI
    method_ = sendrecv_method
    IF (present(method)) method_ = trim(method)

    SELECT CASE (trim(method_))
    CASE ('STD')
       CALL parallel_std
    CASE ('NB')
       CALL parallel_nb
    CASE DEFAULT
       CALL assert(.FALSE., "unsupported SENDRECV_METHOD")
    END SELECT
#else
    CALL single
#endif

  CONTAINS

  SUBROUTINE single
    INTEGER :: i, j

    IF (cycle_x) THEN
       DO i=1, slv_w
          a(i,:) = a(dimx+i,:)
       END DO
       DO i=1, slv_e
          a(slv_w+dimx+i,:) = a(slv_w+i,:)
       END DO
    ELSE IF (present(fill)) THEN
       DO i=1, slv_w
          a(i,:) = fill
       END DO
       DO i=1, slv_e
          a(slv_w+dimx+i,:) = fill
       END DO
    END IF

    IF (cycle_y) THEN
       DO j=1, slv_s
          a(:,j) = a(:,dimy+j)
       END DO
       DO j=1, slv_n
          a(:,slv_s+dimy+j) = a(:,slv_s+j)
       END DO
    ELSE IF (present(fill)) THEN
       DO j=1, slv_s
          a(:,j) = fill
       END DO
       DO j=1, slv_n
          a(:,slv_s+dimy+j) = fill
       END DO
    END IF

    IF (tripolar) THEN
       DO j=1, slv_n
       DO i=1, dimx+slv_w+slv_e
          a(i,dimy+slv_s+j) = a(dimx+slv_w+slv_e+1-i, dimy+slv_n+1-j)
          ! using slv_n instead of slev_s in the RHS is correct since it takes care of cases where ny is odd number (slv_s = slv_n+1)
       END DO
       END DO
    END IF

  END SUBROUTINE single

#ifdef PARALLEL_MPI
  SUBROUTINE parallel_std
    INTEGER :: ierr
    INTEGER :: i, j

    DO j=1, jsize
       DO i=1, slv_e
          sendbuf_r4_w(slv_e*(j-1)+i) = a(slv_w+i,slv_s+j)
       END DO

       DO i=1, slv_w
          sendbuf_r4_e(slv_w*(j-1)+i) = a(slv_w+isize-slv_w+i,slv_s+j)
       END DO
    END DO

    IF (nx == 1) THEN
       CALL mpi_startall(2, mpireq_2dx_r4(:,nx), ierr)
    ELSE IF (nx > 1) THEN
       CALL mpi_startall(4, mpireq_2dx_r4(:,nx), ierr)
    END IF
    CALL mpi_waitall(4, mpireq_2dx_r4(:,nx), MPI_STATUSES_IGNORE, ierr)

    IF (rank_w /= MPI_PROC_NULL) THEN
       DO j=1, jsize
       DO i=1, slv_w
          a(i,slv_s+j) = recvbuf_r4_w(slv_w*(j-1)+i)
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO j=1, jsize
       DO i=1, slv_w
          a(i,slv_s+j) = fill
       END DO
       END DO
    END IF

    IF (rank_e /= MPI_PROC_NULL) THEN
       DO j=1, jsize
       DO i=1, slv_e
          a(slv_w+isize+i,slv_s+j) = recvbuf_r4_e(slv_e*(j-1)+i)
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO j=1, jsize
       DO i=1, slv_e
          a(slv_w+isize+i,slv_s+j) = fill
       END DO
       END DO
    END IF

    DO j=1, slv_n
    DO i=1, isize+slv_w+slv_e
       sendbuf_r4_s((isize+slv_w+slv_e)*(j-1)+i) = a(i,slv_s+j)
    END DO
    END DO

    DO j=1, slv_s
    DO i=1, isize+slv_w+slv_e
       sendbuf_r4_n((isize+slv_w+slv_e)*(j-1)+i) = a(i,jsize+j)
    END DO
    END DO

    IF (ny == 1) THEN
       CALL mpi_startall(2, mpireq_2dy_r4(:,ny), ierr)
    ELSE IF (ny > 1) THEN
       CALL mpi_startall(4, mpireq_2dy_r4(:,ny), ierr)
    END IF
    CALL mpi_waitall(4, mpireq_2dy_r4(:,ny), MPI_STATUSES_IGNORE, ierr)

    IF (rank_s /= MPI_PROC_NULL) THEN
       DO j=1, slv_s
       DO i=1, isize+slv_w+slv_e
          a(i,j) = recvbuf_r4_s((isize+slv_w+slv_e)*(j-1)+i)
      END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO j=1, slv_s
          a(:,j) = fill
       END DO
    END IF

    IF (rank_n /= MPI_PROC_NULL) THEN
       DO j=1, slv_n
       DO i=1, isize+slv_w+slv_e
          a(i,slv_s+jsize+j) = recvbuf_r4_n((isize+slv_w+slv_e)*(j-1)+i)
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO j=1, slv_n
          a(:,slv_s+jsize+j) = fill
       END DO
    END IF

    IF (rank_tp(0,0)/=MPI_PROC_NULL) THEN
       DO j=1, slv_n
       DO i=1, isize+slv_w+slv_e
          sendbuf_r4_n((isize+slv_w+slv_e)*(j-1)+i) = a(isize+slv_w+slv_e+1-i,jsize+slv_n+1-j)
       END DO
       END DO

       CALL mpi_sendrecv(sendbuf_r4_n, (isize+slv_w+slv_e)*slv_n, MPI_REAL4, rank_tp(0,0), 30, &
                         recvbuf_r4_n, (isize+slv_w+slv_e)*slv_n, MPI_REAL4, rank_tp(0,0), 30, comm, MPI_STATUS_IGNORE, ierr)

       DO j=1, slv_n
       DO i=1, isize+slv_w+slv_e
          a(i,jsize+slv_s+j) = recvbuf_r4_n((isize+slv_w+slv_e)*(j-1)+i)
       END DO
       END DO
    END IF

  END SUBROUTINE parallel_std

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE parallel_nb
    INTEGER :: sendreq(8), recvreq(8), ierr
    INTEGER :: i, j

    sendreq(:) = MPI_REQUEST_NULL
    recvreq(:) = MPI_REQUEST_NULL

!$OMP PARALLEL PRIVATE(i, j, ierr)
!$OMP SECTIONS
!$OMP SECTION
    DO j=1, slv_n
    DO i=1, slv_e
       sendbuf_r4_sw(slv_e*(j-1)+i) = a(slv_w+i,slv_s+j)
    END DO
    END DO
!$OMP CRITICAL
    CALL mpi_isend(sendbuf_r4_sw, slv_e*slv_n, MPI_REAL4, rank_sw, tag_sw, comm, sendreq(tag_sw), ierr)
    CALL mpi_irecv(recvbuf_r4_ne, slv_e*slv_n, MPI_REAL4, rank_ne, tag_sw, comm, recvreq(tag_sw), ierr)
!$OMP END CRITICAL

!$OMP SECTION
    DO j=1, slv_n
    DO i=1, slv_w
       sendbuf_r4_se(slv_w*(j-1)+i) = a(isize+i,slv_s+j)
    END DO
    END DO
!$OMP CRITICAL
    CALL mpi_isend(sendbuf_r4_se, slv_w*slv_n, MPI_REAL4, rank_se, tag_se, comm, sendreq(tag_se), ierr)
    CALL mpi_irecv(recvbuf_r4_nw, slv_w*slv_n, MPI_REAL4, rank_nw, tag_se, comm, recvreq(tag_se), ierr)
!$OMP END CRITICAL

!$OMP SECTION
    DO j=1, slv_s
    DO i=1, slv_e
       sendbuf_r4_nw(slv_e*(j-1)+i) = a(slv_w+i,jsize+j)
    END DO
    END DO
!$OMP CRITICAL
    CALL mpi_isend(sendbuf_r4_nw, slv_e*slv_s, MPI_REAL4, rank_nw, tag_nw, comm, sendreq(tag_nw), ierr)
    CALL mpi_irecv(recvbuf_r4_se, slv_e*slv_s, MPI_REAL4, rank_se, tag_nw, comm, recvreq(tag_nw), ierr)
!$OMP END CRITICAL

!$OMP SECTION
    DO j=1, slv_s
    DO i=1, slv_w
       sendbuf_r4_ne(slv_w*(j-1)+i) = a(isize+i,jsize+j)
    END DO
    END DO
!$OMP CRITICAL
    CALL mpi_isend(sendbuf_r4_ne, slv_w*slv_s, MPI_REAL4, rank_ne, tag_ne, comm, sendreq(tag_ne), ierr)
    CALL mpi_irecv(recvbuf_r4_sw, slv_w*slv_s, MPI_REAL4, rank_sw, tag_ne, comm, recvreq(tag_ne), ierr)
!$OMP END CRITICAL

!$OMP SECTION
    DO j=1, jsize
    DO i=1, slv_e
       sendbuf_r4_w(slv_e*(j-1)+i) = a(slv_w+i,slv_s+j)
    END DO
    END DO
!$OMP CRITICAL
    CALL mpi_isend(sendbuf_r4_w, jsize*slv_e, MPI_REAL4, rank_w, tag_w, comm, sendreq(tag_w), ierr)
    CALL mpi_irecv(recvbuf_r4_e, jsize*slv_e, MPI_REAL4, rank_e, tag_w, comm, recvreq(tag_w), ierr)
!$OMP END CRITICAL

!$OMP SECTION
    DO j=1, jsize
    DO i=1, slv_w
       sendbuf_r4_e(slv_w*(j-1)+i) = a(isize+i,slv_s+j)
    END DO
    END DO
!$OMP CRITICAL
    CALL mpi_isend(sendbuf_r4_e, jsize*slv_w, MPI_REAL4, rank_e, tag_e, comm, sendreq(tag_e), ierr)
    CALL mpi_irecv(recvbuf_r4_w, jsize*slv_w, MPI_REAL4, rank_w, tag_e, comm, recvreq(tag_e), ierr)
!$OMP END CRITICAL

!$OMP SECTION
    DO j=1, slv_n
    DO i=1, isize
       sendbuf_r4_s(isize*(j-1)+i) = a(slv_w+i,slv_s+j)
    END DO
    END DO
!$OMP CRITICAL
    CALL mpi_isend(sendbuf_r4_s, isize*slv_n, MPI_REAL4, rank_s, tag_s, comm, sendreq(tag_s), ierr)
    CALL mpi_irecv(recvbuf_r4_n, isize*slv_n, MPI_REAL4, rank_n, tag_s, comm, recvreq(tag_s), ierr)
!$OMP END CRITICAL

!$OMP SECTION
    DO j=1, slv_s
    DO i=1, isize
       sendbuf_r4_n(isize*(j-1)+i) = a(slv_w+i,jsize+j)
    END DO
    END DO
!$OMP CRITICAL
    CALL mpi_isend(sendbuf_r4_n, isize*slv_s, MPI_REAL4, rank_n, tag_n, comm, sendreq(tag_n), ierr)
    CALL mpi_irecv(recvbuf_r4_s, isize*slv_s, MPI_REAL4, rank_s, tag_n, comm, recvreq(tag_n), ierr)
!$OMP END CRITICAL
!$OMP END SECTIONS

!$OMP SECTIONS
!$OMP SECTION
!$OMP CRITICAL
    CALL mpi_wait(recvreq(tag_sw), MPI_STATUS_IGNORE, ierr)
!$OMP END CRITICAL
    IF (rank_ne /= MPI_PROC_NULL) THEN
       DO j=1, slv_n
       DO i=1, slv_e
          a(slv_w+isize+i, slv_s+jsize+j) = recvbuf_r4_ne(slv_e*(j-1)+i)
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO j=1, slv_n
       DO i=1, slv_e
          a(slv_w+isize+i, slv_s+jsize+j) = fill
       END DO
       END DO
    END IF

!$OMP SECTION
!$OMP CRITICAL
    CALL mpi_wait(recvreq(tag_se), MPI_STATUS_IGNORE, ierr)
!$OMP END CRITICAL
    IF (rank_nw /= MPI_PROC_NULL) THEN
       DO j=1, slv_n
       DO i=1, slv_w
          a(i, slv_s+jsize+j) = recvbuf_r4_nw(slv_w*(j-1)+i)
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO j=1, slv_n
       DO i=1, slv_w
          a(i, slv_s+jsize+j) = fill
       END DO
       END DO
    END IF

!$OMP SECTION
!$OMP CRITICAL
    CALL mpi_wait(recvreq(tag_nw), MPI_STATUS_IGNORE, ierr)
!$OMP END CRITICAL
    IF (rank_se /= MPI_PROC_NULL) THEN
       DO j=1, slv_s
       DO i=1, slv_e
          a(slv_w+isize+i,j) = recvbuf_r4_se(slv_e*(j-1)+i)
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO j=1, slv_s
       DO i=1, slv_e
          a(slv_w+isize+i,j) = fill
       END DO
       END DO
    END IF

!$OMP SECTION
!$OMP CRITICAL
    CALL mpi_wait(recvreq(tag_ne), MPI_STATUS_IGNORE, ierr)
!$OMP END CRITICAL
    IF (rank_sw /= MPI_PROC_NULL) THEN
       DO j=1, slv_s
       DO i=1, slv_w
          a(i, j) = recvbuf_r4_sw(slv_w*(j-1)+i)
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO j=1, slv_s
       DO i=1, slv_w
          a(i, j) = fill
       END DO
       END DO
    END IF

!$OMP SECTION
!$OMP CRITICAL
    CALL mpi_wait(recvreq(tag_w), MPI_STATUS_IGNORE, ierr)
!$OMP END CRITICAL
    IF (rank_e /= MPI_PROC_NULL) THEN
       DO j=1, jsize
       DO i=1, slv_e
          a(slv_w+isize+i,slv_s+j) = recvbuf_r4_e(slv_e*(j-1)+i)
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO j=1, jsize
       DO i=1, slv_e
          a(slv_w+isize+i,slv_s+j) = fill
       END DO
       END DO
    END IF

!$OMP SECTION
!$OMP CRITICAL
    CALL mpi_wait(recvreq(tag_e), MPI_STATUS_IGNORE, ierr)
!$OMP END CRITICAL
    IF (rank_w /= MPI_PROC_NULL) THEN
       DO j=1, jsize
       DO i=1, slv_w
          a(i,slv_s+j) = recvbuf_r4_w(slv_w*(j-1)+i)
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO j=1, jsize
       DO i=1, slv_w
          a(i,slv_s+j) = fill
       END DO
       END DO
    END IF

!$OMP SECTION
!$OMP CRITICAL
    CALL mpi_wait(recvreq(tag_s), MPI_STATUS_IGNORE, ierr)
!$OMP END CRITICAL
    IF (rank_n /= MPI_PROC_NULL) THEN
       DO j=1, slv_n
       DO i=1, isize
          a(slv_w+i,slv_s+jsize+j) = recvbuf_r4_n(isize*(j-1)+i)
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO j=1, slv_n
       DO i=1, isize
          a(slv_w+i,slv_s+jsize+j) = fill
       END DO
       END DO
    END IF

!$OMP SECTION
!$OMP CRITICAL
    CALL mpi_wait(recvreq(tag_n), MPI_STATUS_IGNORE, ierr)
!$OMP END CRITICAL
    IF (rank_s /= MPI_PROC_NULL) THEN
       DO j=1, slv_s
       DO i=1, isize
          a(slv_w+i,j) = recvbuf_r4_s(isize*(j-1)+i)
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO j=1, slv_s
       DO i=1, isize
          a(slv_w+i,j) = fill
       END DO
       END DO
    END IF
!$OMP END SECTIONS
!$OMP END PARALLEL

    CALL mpi_waitall(8, sendreq, MPI_STATUSES_IGNORE, ierr)

    IF (rank_tp(0,0)/=MPI_PROC_NULL) THEN
       DO j=1, slv_n
       DO i=1, isize+slv_w+slv_e
          sendbuf_r4_n((isize+slv_w+slv_e)*(j-1)+i) = a(isize+slv_w+slv_e+1-i,jsize+slv_n+1-j)
       END DO
       END DO

       CALL mpi_sendrecv(sendbuf_r4_n, (isize+slv_w+slv_e)*slv_n, MPI_REAL4, rank_tp(0,0), 30, &
                         recvbuf_r4_n, (isize+slv_w+slv_e)*slv_n, MPI_REAL4, rank_tp(0,0), 30, comm, MPI_STATUS_IGNORE, ierr)

       DO j=1, slv_n
       DO i=1, isize+slv_w+slv_e
          a(i,jsize+slv_s+j) = recvbuf_r4_n((isize+slv_w+slv_e)*(j-1)+i)
       END DO
       END DO
    END IF

  END SUBROUTINE parallel_nb
#endif

  END SUBROUTINE update_boundary_2d_r4

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE update_boundary_2d_r8(a, fill, method)
    REAL(8),      INTENT(INOUT)        :: a(:,:)
    REAL(8),      INTENT(IN), OPTIONAL :: fill
    CHARACTER(*), INTENT(IN), OPTIONAL :: method
#ifdef F2008
    CONTIGUOUS a
#endif

    INTEGER :: nx, ny
    INTEGER :: slv_e, slv_w, slv_n, slv_s

    CHARACTER(16) :: method_

    nx = size(a,1) - isize
    ny = size(a,2) - jsize

    CALL assert(nx >= 0 .AND. nx <= 5, "unsported dimension in UPDATE_BOUNDARY_2D")
    CALL assert(ny >= 0 .AND. ny <= 5, "unsported dimension in UPDATE_BOUNDARY_2D")

    slv_e = nx/2
    slv_n = ny/2
    slv_w = (nx+1)/2
    slv_s = (ny+1)/2

#ifdef PARALLEL_MPI
    method_ = sendrecv_method
    IF (present(method)) method_ = trim(method)

    SELECT CASE (trim(method_))
    CASE ('STD')
       CALL parallel_std
    CASE ('NB')
       CALL parallel_nb
    CASE DEFAULT
       CALL assert(.FALSE., "unsupported SENDRECV_METHOD")
    END SELECT
#else
    CALL single
#endif

  CONTAINS

  SUBROUTINE single
    INTEGER :: i, j

    IF (cycle_x) THEN
       DO i=1, slv_w
          a(i,:) = a(dimx+i,:)
       END DO
       DO i=1, slv_e
          a(slv_w+dimx+i,:) = a(slv_w+i,:)
       END DO
    ELSE IF (present(fill)) THEN
       DO i=1, slv_w
          a(i,:) = fill
       END DO
       DO i=1, slv_e
          a(slv_w+dimx+i,:) = fill
       END DO
    END IF

    IF (cycle_y) THEN
       DO j=1, slv_s
          a(:,j) = a(:,dimy+j)
       END DO
       DO j=1, slv_n
          a(:,slv_s+dimy+j) = a(:,slv_s+j)
       END DO
    ELSE IF (present(fill)) THEN
       DO j=1, slv_s
          a(:,j) = fill
       END DO
       DO j=1, slv_n
          a(:,slv_s+dimy+j) = fill
       END DO
    END IF

    IF (tripolar) THEN
       DO j=1, slv_n
       DO i=1, dimx+slv_w+slv_e
          a(i,dimy+slv_s+j) = a(dimx+slv_w+slv_e+1-i, dimy+slv_n+1-j)
       END DO
       END DO
    END IF

  END SUBROUTINE single

#ifdef PARALLEL_MPI
  SUBROUTINE parallel_std
    INTEGER :: ierr
    INTEGER :: i, j

    DO j=1, jsize
       DO i=1, slv_e
          sendbuf_r8_w(slv_e*(j-1)+i) = a(slv_w+i,slv_s+j)
       END DO

       DO i=1, slv_w
          sendbuf_r8_e(slv_w*(j-1)+i) = a(slv_w+isize-slv_w+i,slv_s+j)
       END DO
    END DO

    IF (nx == 1) THEN
       CALL mpi_startall(2, mpireq_2dx_r8(:,nx), ierr)
    ELSE IF (nx > 1) THEN
       CALL mpi_startall(4, mpireq_2dx_r8(:,nx), ierr)
    END IF
    CALL mpi_waitall(4, mpireq_2dx_r8(:,nx), MPI_STATUSES_IGNORE, ierr)

    IF (rank_w /= MPI_PROC_NULL) THEN
       DO j=1, jsize
       DO i=1, slv_w
          a(i,slv_s+j) = recvbuf_r8_w(slv_w*(j-1)+i)
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO j=1, jsize
       DO i=1, slv_w
          a(i,slv_s+j) = fill
       END DO
       END DO
    END IF

    IF (rank_e /= MPI_PROC_NULL) THEN
       DO j=1, jsize
       DO i=1, slv_e
          a(slv_w+isize+i,slv_s+j) = recvbuf_r8_e(slv_e*(j-1)+i)
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO j=1, jsize
       DO i=1, slv_e
          a(slv_w+isize+i,slv_s+j) = fill
       END DO
       END DO
    END IF

    DO j=1, slv_n
    DO i=1, isize+slv_w+slv_e
       sendbuf_r8_s((isize+slv_w+slv_e)*(j-1)+i) = a(i,slv_s+j)
    END DO
    END DO

    DO j=1, slv_s
    DO i=1, isize+slv_w+slv_e
       sendbuf_r8_n((isize+slv_w+slv_e)*(j-1)+i) = a(i,jsize+j)
    END DO
    END DO

    IF (ny == 1) THEN
       CALL mpi_startall(2, mpireq_2dy_r8(:,ny), ierr)
    ELSE IF (ny > 1) THEN
       CALL mpi_startall(4, mpireq_2dy_r8(:,ny), ierr)
    END IF
    CALL mpi_waitall(4, mpireq_2dy_r8(:,ny), MPI_STATUSES_IGNORE, ierr)

    IF (rank_s /= MPI_PROC_NULL) THEN
       DO j=1, slv_s
       DO i=1, isize+slv_w+slv_e
          a(i,j) = recvbuf_r8_s((isize+slv_w+slv_e)*(j-1)+i)
      END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO j=1, slv_s
          a(:,j) = fill
       END DO
    END IF

    IF (rank_n /= MPI_PROC_NULL) THEN
       DO j=1, slv_n
       DO i=1, isize+slv_w+slv_e
          a(i,slv_s+jsize+j) = recvbuf_r8_n((isize+slv_w+slv_e)*(j-1)+i)
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO j=1, slv_n
          a(:,slv_s+jsize+j) = fill
       END DO
    END IF

    IF (rank_tp(0,0)/=MPI_PROC_NULL) THEN
       DO j=1, slv_n
       DO i=1, isize+slv_w+slv_e
          sendbuf_r8_n((isize+slv_w+slv_e)*(j-1)+i) = a(isize+slv_w+slv_e+1-i,jsize+slv_n+1-j)
       END DO
       END DO

       CALL mpi_sendrecv(sendbuf_r8_n, (isize+slv_w+slv_e)*slv_n, MPI_REAL8, rank_tp(0,0), 30, &
                         recvbuf_r8_n, (isize+slv_w+slv_e)*slv_n, MPI_REAL8, rank_tp(0,0), 30, comm, MPI_STATUS_IGNORE, ierr)

       DO j=1, slv_n
       DO i=1, isize+slv_w+slv_e
          a(i,jsize+slv_s+j) = recvbuf_r8_n((isize+slv_w+slv_e)*(j-1)+i)
       END DO
       END DO
    END IF

  END SUBROUTINE parallel_std

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE parallel_nb
    INTEGER :: sendreq(8), recvreq(8), ierr
    INTEGER :: i, j

    sendreq(:) = MPI_REQUEST_NULL
    recvreq(:) = MPI_REQUEST_NULL

!$OMP PARALLEL PRIVATE(i, j, ierr)
!$OMP SECTIONS
!$OMP SECTION
    DO j=1, slv_n
    DO i=1, slv_e
       sendbuf_r8_sw(slv_e*(j-1)+i) = a(slv_w+i,slv_s+j)
    END DO
    END DO
!$OMP CRITICAL
    CALL mpi_isend(sendbuf_r8_sw, slv_e*slv_n, MPI_REAL8, rank_sw, tag_sw, comm, sendreq(tag_sw), ierr)
    CALL mpi_irecv(recvbuf_r8_ne, slv_e*slv_n, MPI_REAL8, rank_ne, tag_sw, comm, recvreq(tag_sw), ierr)
!$OMP END CRITICAL

!$OMP SECTION
    DO j=1, slv_n
    DO i=1, slv_w
       sendbuf_r8_se(slv_w*(j-1)+i) = a(isize+i,slv_s+j)
    END DO
    END DO
!$OMP CRITICAL
    CALL mpi_isend(sendbuf_r8_se, slv_w*slv_n, MPI_REAL8, rank_se, tag_se, comm, sendreq(tag_se), ierr)
    CALL mpi_irecv(recvbuf_r8_nw, slv_w*slv_n, MPI_REAL8, rank_nw, tag_se, comm, recvreq(tag_se), ierr)
!$OMP END CRITICAL

!$OMP SECTION
    DO j=1, slv_s
    DO i=1, slv_e
       sendbuf_r8_nw(slv_e*(j-1)+i) = a(slv_w+i,jsize+j)
    END DO
    END DO
!$OMP CRITICAL
    CALL mpi_isend(sendbuf_r8_nw, slv_e*slv_s, MPI_REAL8, rank_nw, tag_nw, comm, sendreq(tag_nw), ierr)
    CALL mpi_irecv(recvbuf_r8_se, slv_e*slv_s, MPI_REAL8, rank_se, tag_nw, comm, recvreq(tag_nw), ierr)
!$OMP END CRITICAL

!$OMP SECTION
    DO j=1, slv_s
    DO i=1, slv_w
       sendbuf_r8_ne(slv_w*(j-1)+i) = a(isize+i,jsize+j)
    END DO
    END DO
!$OMP CRITICAL
    CALL mpi_isend(sendbuf_r8_ne, slv_w*slv_s, MPI_REAL8, rank_ne, tag_ne, comm, sendreq(tag_ne), ierr)
    CALL mpi_irecv(recvbuf_r8_sw, slv_w*slv_s, MPI_REAL8, rank_sw, tag_ne, comm, recvreq(tag_ne), ierr)
!$OMP END CRITICAL

!$OMP SECTION
    DO j=1, jsize
    DO i=1, slv_e
       sendbuf_r8_w(slv_e*(j-1)+i) = a(slv_w+i,slv_s+j)
    END DO
    END DO
!$OMP CRITICAL
    CALL mpi_isend(sendbuf_r8_w, jsize*slv_e, MPI_REAL8, rank_w, tag_w, comm, sendreq(tag_w), ierr)
    CALL mpi_irecv(recvbuf_r8_e, jsize*slv_e, MPI_REAL8, rank_e, tag_w, comm, recvreq(tag_w), ierr)
!$OMP END CRITICAL

!$OMP SECTION
    DO j=1, jsize
    DO i=1, slv_w
       sendbuf_r8_e(slv_w*(j-1)+i) = a(isize+i,slv_s+j)
    END DO
    END DO
!$OMP CRITICAL
    CALL mpi_isend(sendbuf_r8_e, jsize*slv_w, MPI_REAL8, rank_e, tag_e, comm, sendreq(tag_e), ierr)
    CALL mpi_irecv(recvbuf_r8_w, jsize*slv_w, MPI_REAL8, rank_w, tag_e, comm, recvreq(tag_e), ierr)
!$OMP END CRITICAL

!$OMP SECTION
    DO j=1, slv_n
    DO i=1, isize
       sendbuf_r8_s(isize*(j-1)+i) = a(slv_w+i,slv_s+j)
    END DO
    END DO
!$OMP CRITICAL
    CALL mpi_isend(sendbuf_r8_s, isize*slv_n, MPI_REAL8, rank_s, tag_s, comm, sendreq(tag_s), ierr)
    CALL mpi_irecv(recvbuf_r8_n, isize*slv_n, MPI_REAL8, rank_n, tag_s, comm, recvreq(tag_s), ierr)
!$OMP END CRITICAL

!$OMP SECTION
    DO j=1, slv_s
    DO i=1, isize
       sendbuf_r8_n(isize*(j-1)+i) = a(slv_w+i,jsize+j)
    END DO
    END DO
!$OMP CRITICAL
    CALL mpi_isend(sendbuf_r8_n, isize*slv_s, MPI_REAL8, rank_n, tag_n, comm, sendreq(tag_n), ierr)
    CALL mpi_irecv(recvbuf_r8_s, isize*slv_s, MPI_REAL8, rank_s, tag_n, comm, recvreq(tag_n), ierr)
!$OMP END CRITICAL
!$OMP END SECTIONS

!$OMP SECTIONS
!$OMP SECTION
!$OMP CRITICAL
    CALL mpi_wait(recvreq(tag_sw), MPI_STATUS_IGNORE, ierr)
!$OMP END CRITICAL
    IF (rank_ne /= MPI_PROC_NULL) THEN
       DO j=1, slv_n
       DO i=1, slv_e
          a(slv_w+isize+i, slv_s+jsize+j) = recvbuf_r8_ne(slv_e*(j-1)+i)
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO j=1, slv_n
       DO i=1, slv_e
          a(slv_w+isize+i, slv_s+jsize+j) = fill
       END DO
       END DO
    END IF

!$OMP SECTION
!$OMP CRITICAL
    CALL mpi_wait(recvreq(tag_se), MPI_STATUS_IGNORE, ierr)
!$OMP END CRITICAL
    IF (rank_nw /= MPI_PROC_NULL) THEN
       DO j=1, slv_n
       DO i=1, slv_w
          a(i, slv_s+jsize+j) = recvbuf_r8_nw(slv_w*(j-1)+i)
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO j=1, slv_n
       DO i=1, slv_w
          a(i, slv_s+jsize+j) = fill
       END DO
       END DO
    END IF

!$OMP SECTION
!$OMP CRITICAL
    CALL mpi_wait(recvreq(tag_nw), MPI_STATUS_IGNORE, ierr)
!$OMP END CRITICAL
    IF (rank_se /= MPI_PROC_NULL) THEN
       DO j=1, slv_s
       DO i=1, slv_e
          a(slv_w+isize+i,j) = recvbuf_r8_se(slv_e*(j-1)+i)
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO j=1, slv_s
       DO i=1, slv_e
          a(slv_w+isize+i,j) = fill
       END DO
       END DO
    END IF

!$OMP SECTION
!$OMP CRITICAL
    CALL mpi_wait(recvreq(tag_ne), MPI_STATUS_IGNORE, ierr)
!$OMP END CRITICAL
    IF (rank_sw /= MPI_PROC_NULL) THEN
       DO j=1, slv_s
       DO i=1, slv_w
          a(i, j) = recvbuf_r8_sw(slv_w*(j-1)+i)
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO j=1, slv_s
       DO i=1, slv_w
          a(i, j) = fill
       END DO
       END DO
    END IF

!$OMP SECTION
!$OMP CRITICAL
    CALL mpi_wait(recvreq(tag_w), MPI_STATUS_IGNORE, ierr)
!$OMP END CRITICAL
    IF (rank_e /= MPI_PROC_NULL) THEN
       DO j=1, jsize
       DO i=1, slv_e
          a(slv_w+isize+i,slv_s+j) = recvbuf_r8_e(slv_e*(j-1)+i)
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO j=1, jsize
       DO i=1, slv_e
          a(slv_w+isize+i,slv_s+j) = fill
       END DO
       END DO
    END IF

!$OMP SECTION
!$OMP CRITICAL
    CALL mpi_wait(recvreq(tag_e), MPI_STATUS_IGNORE, ierr)
!$OMP END CRITICAL
    IF (rank_w /= MPI_PROC_NULL) THEN
       DO j=1, jsize
       DO i=1, slv_w
          a(i,slv_s+j) = recvbuf_r8_w(slv_w*(j-1)+i)
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO j=1, jsize
       DO i=1, slv_w
          a(i,slv_s+j) = fill
       END DO
       END DO
    END IF

!$OMP SECTION
!$OMP CRITICAL
    CALL mpi_wait(recvreq(tag_s), MPI_STATUS_IGNORE, ierr)
!$OMP END CRITICAL
    IF (rank_n /= MPI_PROC_NULL) THEN
       DO j=1, slv_n
       DO i=1, isize
          a(slv_w+i,slv_s+jsize+j) = recvbuf_r8_n(isize*(j-1)+i)
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO j=1, slv_n
       DO i=1, isize
          a(slv_w+i,slv_s+jsize+j) = fill
       END DO
       END DO
    END IF

!$OMP SECTION
!$OMP CRITICAL
    CALL mpi_wait(recvreq(tag_n), MPI_STATUS_IGNORE, ierr)
!$OMP END CRITICAL
    IF (rank_s /= MPI_PROC_NULL) THEN
       DO j=1, slv_s
       DO i=1, isize
          a(slv_w+i,j) = recvbuf_r8_s(isize*(j-1)+i)
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO j=1, slv_s
       DO i=1, isize
          a(slv_w+i,j) = fill
       END DO
       END DO
    END IF
!$OMP END SECTIONS
!$OMP END PARALLEL

    CALL mpi_waitall(8, sendreq, MPI_STATUSES_IGNORE, ierr)

    IF (rank_tp(0,0)/=MPI_PROC_NULL) THEN
       DO j=1, slv_n
       DO i=1, isize+slv_w+slv_e
          sendbuf_r8_n((isize+slv_w+slv_e)*(j-1)+i) = a(isize+slv_w+slv_e+1-i,jsize+slv_n+1-j)
       END DO
       END DO

       CALL mpi_sendrecv(sendbuf_r8_n, (isize+slv_w+slv_e)*slv_n, MPI_REAL8, rank_tp(0,0), 30, &
                         recvbuf_r8_n, (isize+slv_w+slv_e)*slv_n, MPI_REAL8, rank_tp(0,0), 30, comm, MPI_STATUS_IGNORE, ierr)

       DO j=1, slv_n
       DO i=1, isize+slv_w+slv_e
          a(i,jsize+slv_s+j) = recvbuf_r8_n((isize+slv_w+slv_e)*(j-1)+i)
       END DO
       END DO
    END IF

  END SUBROUTINE parallel_nb
#endif

  END SUBROUTINE update_boundary_2d_r8

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE update_boundary_3d_r4(a, fill, method)
    REAL(4), INTENT(INOUT) :: a(:,:,:)
    REAL(4), INTENT(IN), OPTIONAL :: fill
    CHARACTER(*), INTENT(IN), OPTIONAL :: method

    INTEGER :: nx, ny, nz
    INTEGER :: slv_e, slv_w, slv_n, slv_s, slv_u, slv_l

    CHARACTER(16) :: method_

    nx = size(a,1) - isize
    ny = size(a,2) - jsize
    nz = size(a,3) - ksize

    CALL assert(nx >= 0 .AND. nx <= 5, "unsported dimension in UPDATE_BOUNDARY_3D")
    CALL assert(ny >= 0 .AND. ny <= 5, "unsported dimension in UPDATE_BOUNDARY_3D")
    CALL assert(nz >= 0 .AND. nz <= 5, "unsported dimension in UPDATE_BOUNDARY_3D")

    slv_e = nx/2
    slv_n = ny/2
    slv_u = nz/2
    slv_w = (nx+1)/2
    slv_s = (ny+1)/2
    slv_l = (nz+1)/2

#ifdef PARALLEL_MPI
    method_ = sendrecv_method
    IF (present(method)) method_ = trim(method)

    SELECT CASE (trim(method_))
    CASE ('STD')
       CALL parallel_std
    CASE ('NB')
       CALL parallel_nb
    CASE DEFAULT
       CALL assert(.FALSE., "unsupported SENDRECV_METHOD")
    END SELECT
#else
    CALL single
#endif

  CONTAINS

  SUBROUTINE single
    INTEGER :: i, j, k

    IF (cycle_x) THEN
       DO i=1, slv_w
          a(i,:,:) = a(dimx+i,:,:)
       END DO
       DO i=1, slv_e
          a(slv_w+dimx+i,:,:) = a(slv_w+i,:,:)
       END DO
    ELSE IF (present(fill)) THEN
       DO i=1, slv_w
          a(i,:,:) = fill
       END DO
       DO i=1, slv_e
          a(slv_w+dimx+i,:,:) = fill
       END DO
    END IF

    IF (cycle_y) THEN
       DO j=1, slv_s
          a(:,j,:) = a(:,dimy+j,:)
       END DO
       DO j=1, slv_n
          a(:,slv_s+dimy+j,:) = a(:,slv_s+j,:)
       END DO
    ELSE IF (present(fill)) THEN
       DO j=1, slv_s
          a(:,j,:) = fill
       END DO
       DO j=1, slv_n
          a(:,slv_s+dimy+j,:) = fill
       END DO
    END IF

    IF (cycle_z) THEN
       DO k=1, slv_l
          a(:,:,k) = a(:,:,dimz+k)
       END DO
       DO k=1, slv_u
          a(:,:,slv_l+dimz+k) = a(:,:,slv_l+k)
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, slv_l
          a(:,:,k) = fill
       END DO
       DO k=1, slv_u
          a(:,:,slv_l+dimz+k) = fill
       END DO
    END IF

    IF (tripolar) THEN
       DO j=1, slv_n
       DO i=1, dimx+slv_w+slv_e
          a(i,dimy+slv_s+j,:) = a(dimx+slv_w+slv_e+1-i, dimy+slv_n+1-j,:)
       END DO
       END DO
    END IF

  END SUBROUTINE single

#ifdef PARALLEL_MPI
  SUBROUTINE parallel_std
    INTEGER :: ierr
    INTEGER :: i, j, k

    DO k=1, ksize
    DO j=1, jsize
       DO i=1, slv_e
          sendbuf_r4_w(jsize*slv_e*(k-1)+slv_e*(j-1)+i) = a(slv_w+i,slv_s+j,slv_l+k)
       END DO

       DO i=1, slv_w
          sendbuf_r4_e(jsize*slv_w*(k-1)+slv_w*(j-1)+i) = a(isize+i,slv_s+j,slv_l+k)
       END DO
    END DO
    END DO

    IF (nx == 1) THEN
       CALL mpi_startall(2, mpireq_3dx_r4(:,nx), ierr)
    ELSE IF (nx > 1) THEN
       CALL mpi_startall(4, mpireq_3dx_r4(:,nx), ierr)
    END IF
    CALL mpi_waitall(4, mpireq_3dx_r4(:,nx), MPI_STATUSES_IGNORE, ierr)

    IF (rank_w /= MPI_PROC_NULL) THEN
       DO k=1, ksize
       DO j=1, jsize
       DO i=1, slv_w
          a(i,slv_s+j,slv_l+k) = recvbuf_r4_w(jsize*slv_w*(k-1)+slv_w*(j-1)+i)
       END DO
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, ksize
       DO j=1, jsize
       DO i=1, slv_w
          a(i,slv_s+j,slv_l+k) = fill
       END DO
       END DO
       END DO
    END IF

    IF (rank_e /= MPI_PROC_NULL) THEN
       DO k=1, ksize
       DO j=1, jsize
       DO i=1, slv_e
          a(slv_w+isize+i,slv_s+j,slv_l+k) = recvbuf_r4_e(jsize*slv_e*(k-1)+slv_e*(j-1)+i)
       END DO
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, ksize
       DO j=1, jsize
       DO i=1, slv_e
          a(slv_w+isize+i,slv_s+j,slv_l+k) = fill
       END DO
       END DO
       END DO
    END IF

    DO k=1, ksize
       DO j=1, slv_n
       DO i=1, isize+slv_w+slv_e
          sendbuf_r4_s((isize+slv_w+slv_e)*slv_n*(k-1)+(isize+slv_w+slv_e)*(j-1)+i) = a(i,slv_s+j,slv_l+k)
       END DO
       END DO

       DO j=1, slv_s
       DO i=1, isize+slv_w+slv_e
          sendbuf_r4_n((isize+slv_w+slv_e)*slv_s*(k-1)+(isize+slv_w+slv_e)*(j-1)+i) = a(i,jsize+j,slv_l+k)
       END DO
       END DO
    END DO

    IF (ny == 1) THEN
       CALL mpi_startall(2, mpireq_3dy_r4(:,ny), ierr)
    ELSE IF (ny > 1) THEN
       CALL mpi_startall(4, mpireq_3dy_r4(:,ny), ierr)
    END IF
    CALL mpi_waitall(4, mpireq_3dy_r4(:,ny), MPI_STATUSES_IGNORE, ierr)

    IF (rank_s /= MPI_PROC_NULL) THEN
       DO k=1, ksize
       DO j=1, slv_s
       DO i=1, isize+slv_w+slv_e
          a(i,j,slv_l+k) = recvbuf_r4_s((isize+slv_w+slv_e)*slv_s*(k-1)+(isize+slv_w+slv_e)*(j-1)+i)
       END DO
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, ksize
       DO j=1, slv_s
          a(:,j,slv_l+k) = fill
       END DO
       END DO
    END IF

    IF (rank_n /= MPI_PROC_NULL) THEN
       DO k=1, ksize
       DO j=1, slv_n
       DO i=1, isize+slv_w+slv_e
          a(i,slv_s+jsize+j,slv_l+k) = recvbuf_r4_n((isize+slv_w+slv_e)*slv_n*(k-1)+(isize+slv_w+slv_e)*(j-1)+i)
       END DO
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, ksize
       DO j=1, slv_n
          a(:,slv_s+jsize+j,slv_l+k) = fill
       END DO
       END DO
    END IF

#ifdef PARALLEL3D
    DO k=1, slv_u
       DO j=1, jsize+slv_s+slv_n
       DO i=1, isize+slv_w+slv_e
          sendbuf_r4_l((isize+slv_w+slv_e)*(jsize+slv_s+slv_n)*(k-1)+(isize+slv_w+slv_e)*(j-1)+i) = a(i,j,slv_l+k)
       END DO
       END DO
    END DO

    DO k=1, slv_l
       DO j=1, jsize+slv_s+slv_n
       DO i=1, isize+slv_w+slv_e
          sendbuf_r4_u((isize+slv_w+slv_e)*(jsize+slv_s+slv_n)*(k-1)+(isize+slv_w+slv_e)*(j-1)+i) = a(i,j,ksize+k)
       END DO
       END DO
    END DO

    IF (nz == 1) THEN
       CALL mpi_startall(2, mpireq_3dz_r4(:,nz), ierr)
    ELSE IF (nz > 1) THEN
       CALL mpi_startall(4, mpireq_3dz_r4(:,nz), ierr)
    END IF
    CALL mpi_waitall(4, mpireq_3dz_r4(:,nz), MPI_STATUSES_IGNORE, ierr)

    IF (rank_l /= MPI_PROC_NULL) THEN
       DO k=1, slv_l
       DO j=1, jsize+slv_s+slv_n
       DO i=1, isize+slv_w+slv_e
          a(i,j,k) = recvbuf_r4_l((isize+slv_w+slv_e)*(jsize+slv_s+slv_n)*(k-1)+(isize+slv_w+slv_e)*(j-1)+i)
       END DO
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, slv_l
          a(:,:,k) = fill
       END DO
    END IF

    IF (rank_u /= MPI_PROC_NULL) THEN
       DO k=1, slv_u
       DO j=1, jsize+slv_s+slv_n
       DO i=1, isize+slv_w+slv_e
          a(i,j,slv_l+ksize+k) = recvbuf_r4_u((isize+slv_w+slv_e)*(jsize+slv_s+slv_n)*(k-1)+(isize+slv_e+slv_w)*(j-1)+i)
       END DO
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, slv_u
          a(:,:,slv_l+ksize+k) = fill
       END DO
    END IF
#else
    IF (present(fill)) THEN
       DO k=1, slv_l
          a(:,:,k) = fill
       END DO
       DO k=1, slv_u
          a(:,:,slv_l+ksize+k) = fill
       END DO
    END IF
#endif

    IF (rank_tp(0,0)/=MPI_PROC_NULL) THEN
       DO k=1, ksize+slv_l+slv_u
       DO j=1, slv_n
       DO i=1, isize+slv_w+slv_e
          sendbuf_r4_n((isize+slv_w+slv_e)*slv_n*(k-1)+(isize+slv_w+slv_e)*(j-1)+i) = a(isize+slv_w+slv_e+1-i,jsize+slv_n+1-j,k)
       END DO
       END DO
       END DO

       CALL mpi_sendrecv(sendbuf_r4_n, (isize+slv_w+slv_e)*slv_n*(ksize+slv_l+slv_u), MPI_REAL4, rank_tp(0,0), 30, &
                         recvbuf_r4_n, (isize+slv_w+slv_e)*slv_n*(ksize+slv_l+slv_u), MPI_REAL4, rank_tp(0,0), 30, comm, MPI_STATUS_IGNORE, ierr)

       DO k=1, ksize+slv_l+slv_u
       DO j=1, slv_n
       DO i=1, isize+slv_w+slv_e
          a(i,jsize+slv_s+j,k) = recvbuf_r4_n((isize+slv_w+slv_e)*slv_n*(k-1)+(isize+slv_w+slv_e)*(j-1)+i)
       END DO
       END DO
       END DO
    END IF

  END SUBROUTINE parallel_std

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE parallel_nb
    INTEGER :: sendreq(n_sendrecv), recvreq(n_sendrecv), ierr
    INTEGER :: i, j, k

    sendreq(:) = MPI_REQUEST_NULL
    recvreq(:) = MPI_REQUEST_NULL

!$OMP PARALLEL PRIVATE(i, j, k, ierr)
!$OMP SECTIONS
!$OMP SECTION
    DO k=1, ksize
    DO j=1, slv_n
    DO i=1, slv_e
       sendbuf_r4_sw(slv_e*slv_n*(k-1)+slv_e*(j-1)+i) = a(slv_w+i,slv_s+j,slv_l+k)
    END DO
    END DO
    END DO
!$OMP CRITICAL
    CALL mpi_isend(sendbuf_r4_sw, slv_e*slv_n*ksize, MPI_REAL4, rank_sw, tag_sw, comm, sendreq(tag_sw), ierr)
    CALL mpi_irecv(recvbuf_r4_ne, slv_e*slv_n*ksize, MPI_REAL4, rank_ne, tag_sw, comm, recvreq(tag_sw), ierr)
!$OMP END CRITICAL

!$OMP SECTION
    DO k=1, ksize
    DO j=1, slv_n
    DO i=1, slv_w
       sendbuf_r4_se(slv_w*slv_n*(k-1)+slv_w*(j-1)+i) = a(isize+i,slv_s+j,slv_l+k)
    END DO
    END DO
    END DO
!$OMP CRITICAL
    CALL mpi_isend(sendbuf_r4_se, slv_w*slv_n*ksize, MPI_REAL4, rank_se, tag_se, comm, sendreq(tag_se), ierr)
    CALL mpi_irecv(recvbuf_r4_nw, slv_w*slv_n*ksize, MPI_REAL4, rank_nw, tag_se, comm, recvreq(tag_se), ierr)
!$OMP END CRITICAL

!$OMP SECTION
    DO k=1, ksize
    DO j=1, slv_s
    DO i=1, slv_e
       sendbuf_r4_nw(slv_e*slv_s*(k-1)+slv_e*(j-1)+i) = a(slv_w+i,jsize+j,slv_l+k)
    END DO
    END DO
    END DO
!$OMP CRITICAL
    CALL mpi_isend(sendbuf_r4_nw, slv_e*slv_s*ksize, MPI_REAL4, rank_nw, tag_nw, comm, sendreq(tag_nw), ierr)
    CALL mpi_irecv(recvbuf_r4_se, slv_e*slv_s*ksize, MPI_REAL4, rank_se, tag_nw, comm, recvreq(tag_nw), ierr)
!$OMP END CRITICAL

!$OMP SECTION
    DO k=1, ksize
    DO j=1, slv_s
    DO i=1, slv_w
       sendbuf_r4_ne(slv_w*slv_s*(k-1)+slv_w*(j-1)+i) = a(isize+i,jsize+j,slv_l+k)
    END DO
    END DO
    END DO
!$OMP CRITICAL
    CALL mpi_isend(sendbuf_r4_ne, slv_w*slv_s*ksize, MPI_REAL4, rank_ne, tag_ne, comm, sendreq(tag_ne), ierr)
    CALL mpi_irecv(recvbuf_r4_sw, slv_w*slv_s*ksize, MPI_REAL4, rank_sw, tag_ne, comm, recvreq(tag_ne), ierr)
!$OMP END CRITICAL

!$OMP SECTION
    DO k=1, ksize
    DO j=1, jsize
    DO i=1, slv_e
       sendbuf_r4_w(slv_e*jsize*(k-1)+slv_e*(j-1)+i) = a(slv_w+i,slv_s+j,slv_l+k)
    END DO
    END DO
    END DO
!$OMP CRITICAL
    CALL mpi_isend(sendbuf_r4_w, slv_e*jsize*ksize, MPI_REAL4, rank_w, tag_w, comm, sendreq(tag_w), ierr)
    CALL mpi_irecv(recvbuf_r4_e, slv_e*jsize*ksize, MPI_REAL4, rank_e, tag_w, comm, recvreq(tag_w), ierr)
!$OMP END CRITICAL

!$OMP SECTION
    DO k=1, ksize
    DO j=1, jsize
    DO i=1, slv_w
       sendbuf_r4_e(slv_w*jsize*(k-1)+slv_w*(j-1)+i) = a(isize+i,slv_s+j,slv_l+k)
    END DO
    END DO
    END DO
!$OMP CRITICAL
    CALL mpi_isend(sendbuf_r4_e, slv_w*jsize*ksize, MPI_REAL4, rank_e, tag_e, comm, sendreq(tag_e), ierr)
    CALL mpi_irecv(recvbuf_r4_w, slv_w*jsize*ksize, MPI_REAL4, rank_w, tag_e, comm, recvreq(tag_e), ierr)
!$OMP END CRITICAL

!$OMP SECTION
    DO k=1, ksize
    DO j=1, slv_n
    DO i=1, isize
       sendbuf_r4_s(isize*slv_n*(k-1)+isize*(j-1)+i) = a(slv_w+i,slv_s+j,slv_l+k)
    END DO
    END DO
    END DO
!$OMP CRITICAL
    CALL mpi_isend(sendbuf_r4_s, isize*slv_n*ksize, MPI_REAL4, rank_s, tag_s, comm, sendreq(tag_s), ierr)
    CALL mpi_irecv(recvbuf_r4_n, isize*slv_n*ksize, MPI_REAL4, rank_n, tag_s, comm, recvreq(tag_s), ierr)
!$OMP END CRITICAL

!$OMP SECTION
    DO k=1, ksize
    DO j=1, slv_s
    DO i=1, isize
       sendbuf_r4_n(isize*slv_s*(k-1)+isize*(j-1)+i) = a(slv_w+i,jsize+j,slv_l+k)
    END DO
    END DO
    END DO
!$OMP CRITICAL
    CALL mpi_isend(sendbuf_r4_n, isize*slv_s*ksize, MPI_REAL4, rank_n, tag_n, comm, sendreq(tag_n), ierr)
    CALL mpi_irecv(recvbuf_r4_s, isize*slv_s*ksize, MPI_REAL4, rank_s, tag_n, comm, recvreq(tag_n), ierr)
!$OMP END CRITICAL

#ifdef PARALLEL3D
!$OMP SECTION
    DO k=1, slv_u
    DO j=1, slv_n
    DO i=1, slv_e
       sendbuf_r4_lsw(slv_e*slv_n*(k-1)+slv_e*(j-1)+i) = a(slv_w+i,slv_s+j,slv_l+k)
    END DO
    END DO
    END DO
!$OMP CRITICAL
    CALL mpi_isend(sendbuf_r4_lsw, slv_e*slv_n*slv_u, MPI_REAL4, rank_lsw, tag_lsw, comm, sendreq(tag_lsw), ierr)
    CALL mpi_irecv(recvbuf_r4_une, slv_e*slv_n*slv_u, MPI_REAL4, rank_une, tag_lsw, comm, recvreq(tag_lsw), ierr)
!$OMP END CRITICAL

!$OMP SECTION
    DO k=1, slv_u
    DO j=1, slv_n
    DO i=1, slv_w
       sendbuf_r4_lse(slv_w*slv_n*(k-1)+slv_w*(j-1)+i) = a(isize+i,slv_s+j,slv_l+k)
    END DO
    END DO
    END DO
!$OMP CRITICAL
    CALL mpi_isend(sendbuf_r4_lse, slv_w*slv_n*slv_u, MPI_REAL4, rank_lse, tag_lse, comm, sendreq(tag_lse), ierr)
    CALL mpi_irecv(recvbuf_r4_unw, slv_w*slv_n*slv_u, MPI_REAL4, rank_unw, tag_lse, comm, recvreq(tag_lse), ierr)
!$OMP END CRITICAL

!$OMP SECTION
    DO k=1, slv_u
    DO j=1, slv_s
    DO i=1, slv_e
       sendbuf_r4_lnw(slv_e*slv_s*(k-1)+slv_e*(j-1)+i) = a(slv_w+i,jsize+j,slv_l+k)
    END DO
    END DO
    END DO
!$OMP CRITICAL
    CALL mpi_isend(sendbuf_r4_lnw, slv_e*slv_s*slv_u, MPI_REAL4, rank_lnw, tag_lnw, comm, sendreq(tag_lnw), ierr)
    CALL mpi_irecv(recvbuf_r4_use, slv_e*slv_s*slv_u, MPI_REAL4, rank_use, tag_lnw, comm, recvreq(tag_lnw), ierr)
!$OMP END CRITICAL

!$OMP SECTION
    DO k=1, slv_u
    DO j=1, slv_s
    DO i=1, slv_w
       sendbuf_r4_lne(slv_w*slv_s*(k-1)+slv_w*(j-1)+i) = a(isize+i,jsize+j,slv_l+k)
    END DO
    END DO
    END DO
!$OMP CRITICAL
    CALL mpi_isend(sendbuf_r4_lne, slv_w*slv_s*slv_u, MPI_REAL4, rank_lne, tag_lne, comm, sendreq(tag_lne), ierr)
    CALL mpi_irecv(recvbuf_r4_usw, slv_w*slv_s*slv_u, MPI_REAL4, rank_usw, tag_lne, comm, recvreq(tag_lne), ierr)
!$OMP END CRITICAL

!$OMP SECTION
    DO k=1, slv_l
    DO j=1, slv_n
    DO i=1, slv_e
       sendbuf_r4_usw(slv_e*slv_n*(k-1)+slv_e*(j-1)+i) = a(slv_w+i,slv_s+j,ksize+k)
    END DO
    END DO
    END DO
!$OMP CRITICAL
    CALL mpi_isend(sendbuf_r4_usw, slv_e*slv_n*slv_l, MPI_REAL4, rank_usw, tag_usw, comm, sendreq(tag_usw), ierr)
    CALL mpi_irecv(recvbuf_r4_lne, slv_e*slv_n*slv_l, MPI_REAL4, rank_lne, tag_usw, comm, recvreq(tag_usw), ierr)
!$OMP END CRITICAL

!$OMP SECTION
    DO k=1, slv_l
    DO j=1, slv_n
    DO i=1, slv_w
       sendbuf_r4_use(slv_w*slv_n*(k-1)+slv_w*(j-1)+i) = a(isize+i,slv_s+j,ksize+k)
    END DO
    END DO
    END DO
!$OMP CRITICAL
    CALL mpi_isend(sendbuf_r4_use, slv_w*slv_n*slv_l, MPI_REAL4, rank_use, tag_use, comm, sendreq(tag_use), ierr)
    CALL mpi_irecv(recvbuf_r4_lnw, slv_w*slv_n*slv_l, MPI_REAL4, rank_lnw, tag_use, comm, recvreq(tag_use), ierr)
!$OMP END CRITICAL

!$OMP SECTION
    DO k=1, slv_l
    DO j=1, slv_s
    DO i=1, slv_e
       sendbuf_r4_unw(slv_e*slv_s*(k-1)+slv_e*(j-1)+i) = a(slv_w+i,jsize+j,ksize+k)
    END DO
    END DO
    END DO
!$OMP CRITICAL
    CALL mpi_isend(sendbuf_r4_unw, slv_e*slv_s*slv_l, MPI_REAL4, rank_unw, tag_unw, comm, sendreq(tag_unw), ierr)
    CALL mpi_irecv(recvbuf_r4_lse, slv_e*slv_s*slv_l, MPI_REAL4, rank_lse, tag_unw, comm, recvreq(tag_unw), ierr)
!$OMP END CRITICAL

!$OMP SECTION
    DO k=1, slv_l
    DO j=1, slv_s
    DO i=1, slv_w
       sendbuf_r4_une(slv_w*slv_s*(k-1)+slv_w*(j-1)+i) = a(isize+i,jsize+j,ksize+k)
    END DO
    END DO
    END DO
!$OMP CRITICAL
    CALL mpi_isend(sendbuf_r4_une, slv_w*slv_s*slv_l, MPI_REAL4, rank_une, tag_une, comm, sendreq(tag_une), ierr)
    CALL mpi_irecv(recvbuf_r4_lsw, slv_w*slv_s*slv_l, MPI_REAL4, rank_lsw, tag_une, comm, recvreq(tag_une), ierr)
!$OMP END CRITICAL

!$OMP SECTION
    DO k=1, slv_u
    DO j=1, jsize
    DO i=1, slv_e
       sendbuf_r4_lw(slv_e*jsize*(k-1)+slv_e*(j-1)+i) = a(slv_w+i,slv_s+j,slv_l+k)
    END DO
    END DO
    END DO
!$OMP CRITICAL
    CALL mpi_isend(sendbuf_r4_lw, slv_e*jsize*slv_u, MPI_REAL4, rank_lw, tag_lw, comm, sendreq(tag_lw), ierr)
    CALL mpi_irecv(recvbuf_r4_ue, slv_e*jsize*slv_u, MPI_REAL4, rank_ue, tag_lw, comm, recvreq(tag_lw), ierr)
!$OMP END CRITICAL

!$OMP SECTION
    DO k=1, slv_u
    DO j=1, jsize
    DO i=1, slv_w
       sendbuf_r4_le(slv_w*jsize*(k-1)+slv_w*(j-1)+i) = a(isize+i,slv_s+j,slv_l+k)
    END DO
    END DO
    END DO
!$OMP CRITICAL
    CALL mpi_isend(sendbuf_r4_le, slv_w*jsize*slv_u, MPI_REAL4, rank_le, tag_le, comm, sendreq(tag_le), ierr)
    CALL mpi_irecv(recvbuf_r4_uw, slv_w*jsize*slv_u, MPI_REAL4, rank_uw, tag_le, comm, recvreq(tag_le), ierr)
!$OMP END CRITICAL

!$OMP SECTION
    DO k=1, slv_u
    DO j=1, slv_n
    DO i=1, isize
       sendbuf_r4_ls(isize*slv_n*(k-1)+isize*(j-1)+i) = a(slv_w+i,slv_s+j,slv_l+k)
    END DO
    END DO
    END DO
!$OMP CRITICAL
    CALL mpi_isend(sendbuf_r4_ls, isize*slv_n*slv_u, MPI_REAL4, rank_ls, tag_ls, comm, sendreq(tag_ls), ierr)
    CALL mpi_irecv(recvbuf_r4_un, isize*slv_n*slv_u, MPI_REAL4, rank_un, tag_ls, comm, recvreq(tag_ls), ierr)
!$OMP END CRITICAL

!$OMP SECTION
    DO k=1, slv_u
    DO j=1, slv_s
    DO i=1, isize
       sendbuf_r4_ln(isize*slv_s*(k-1)+isize*(j-1)+i) = a(slv_w+i,jsize+j,slv_l+k)
    END DO
    END DO
    END DO
!$OMP CRITICAL
    CALL mpi_isend(sendbuf_r4_ln, isize*slv_s*slv_u, MPI_REAL4, rank_ln, tag_ln, comm, sendreq(tag_ln), ierr)
    CALL mpi_irecv(recvbuf_r4_us, isize*slv_s*slv_u, MPI_REAL4, rank_us, tag_ln, comm, recvreq(tag_ln), ierr)
!$OMP END CRITICAL

!$OMP SECTION
    DO k=1, slv_l
    DO j=1, jsize
    DO i=1, slv_e
       sendbuf_r4_uw(slv_e*jsize*(k-1)+slv_e*(j-1)+i) = a(slv_w+i,slv_s+j,ksize+k)
    END DO
    END DO
    END DO
!$OMP CRITICAL
    CALL mpi_isend(sendbuf_r4_uw, slv_e*jsize*slv_l, MPI_REAL4, rank_uw, tag_uw, comm, sendreq(tag_uw), ierr)
    CALL mpi_irecv(recvbuf_r4_le, slv_e*jsize*slv_l, MPI_REAL4, rank_le, tag_uw, comm, recvreq(tag_uw), ierr)
!$OMP END CRITICAL

!$OMP SECTION
    DO k=1, slv_l
    DO j=1, jsize
    DO i=1, slv_w
       sendbuf_r4_ue(slv_w*jsize*(k-1)+slv_w*(j-1)+i) = a(isize+i,slv_s+j,ksize+k)
    END DO
    END DO
    END DO
!$OMP CRITICAL
    CALL mpi_isend(sendbuf_r4_ue, slv_w*jsize*slv_l, MPI_REAL4, rank_ue, tag_ue, comm, sendreq(tag_ue), ierr)
    CALL mpi_irecv(recvbuf_r4_lw, slv_w*jsize*slv_l, MPI_REAL4, rank_lw, tag_ue, comm, recvreq(tag_ue), ierr)
!$OMP END CRITICAL

!$OMP SECTION
    DO k=1, slv_l
    DO j=1, slv_n
    DO i=1, isize
       sendbuf_r4_us(isize*slv_n*(k-1)+isize*(j-1)+i) = a(slv_w+i,slv_s+j,ksize+k)
    END DO
    END DO
    END DO
!$OMP CRITICAL
    CALL mpi_isend(sendbuf_r4_us, isize*slv_n*slv_l, MPI_REAL4, rank_us, tag_us, comm, sendreq(tag_us), ierr)
    CALL mpi_irecv(recvbuf_r4_ln, isize*slv_n*slv_l, MPI_REAL4, rank_ln, tag_us, comm, recvreq(tag_us), ierr)
!$OMP END CRITICAL

!$OMP SECTION
    DO k=1, slv_l
    DO j=1, slv_s
    DO i=1, isize
       sendbuf_r4_un(isize*slv_s*(k-1)+isize*(j-1)+i) = a(slv_w+i,jsize+j,ksize+k)
    END DO
    END DO
    END DO
!$OMP CRITICAL
    CALL mpi_isend(sendbuf_r4_un, isize*slv_s*slv_l, MPI_REAL4, rank_un, tag_un, comm, sendreq(tag_un), ierr)
    CALL mpi_irecv(recvbuf_r4_ls, isize*slv_s*slv_l, MPI_REAL4, rank_ls, tag_un, comm, recvreq(tag_un), ierr)
!$OMP END CRITICAL

!$OMP SECTION
    DO k=1, slv_u
    DO j=1, jsize
    DO i=1, isize
       sendbuf_r4_l(isize*jsize*(k-1)+isize*(j-1)+i) = a(slv_w+i,slv_s+j,slv_l+k)
    END DO
    END DO
    END DO
!$OMP CRITICAL
    CALL mpi_isend(sendbuf_r4_l, isize*jsize*slv_u, MPI_REAL4, rank_l, tag_l, comm, sendreq(tag_l), ierr)
    CALL mpi_irecv(recvbuf_r4_u, isize*jsize*slv_u, MPI_REAL4, rank_u, tag_l, comm, recvreq(tag_l), ierr)
!$OMP END CRITICAL

!$OMP SECTION
    DO k=1, slv_l
    DO j=1, jsize
    DO i=1, isize
       sendbuf_r4_u(isize*jsize*(k-1)+isize*(j-1)+i) = a(slv_w+i,slv_s+j,ksize+k)
    END DO
    END DO
    END DO
!$OMP CRITICAL
    CALL mpi_isend(sendbuf_r4_u, isize*jsize*slv_l, MPI_REAL4, rank_u, tag_u, comm, sendreq(tag_u), ierr)
    CALL mpi_irecv(recvbuf_r4_l, isize*jsize*slv_l, MPI_REAL4, rank_l, tag_u, comm, recvreq(tag_u), ierr)
!$OMP END CRITICAL
#endif
!$OMP END SECTIONS


!$OMP SECTIONS
!$OMP SECTION
!$OMP CRITICAL
    CALL mpi_wait(recvreq(tag_sw), MPI_STATUS_IGNORE, ierr)
!$OMP END CRITICAL
    IF (rank_ne /= MPI_PROC_NULL) THEN
       DO k=1, ksize
       DO j=1, slv_n
       DO i=1, slv_e
          a(slv_w+isize+i,slv_s+jsize+j,slv_l+k) = recvbuf_r4_ne(slv_e*slv_n*(k-1)+slv_e*(j-1)+i)
       END DO
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, ksize
       DO j=1, slv_n
       DO i=1, slv_e
          a(slv_w+isize+i,slv_s+jsize+j,slv_l+k) = fill
       END DO
       END DO
       END DO
    END IF

!$OMP SECTION
!$OMP CRITICAL
    CALL mpi_wait(recvreq(tag_se), MPI_STATUS_IGNORE, ierr)
!$OMP END CRITICAL
    IF (rank_nw /= MPI_PROC_NULL) THEN
       DO k=1, ksize
       DO j=1, slv_n
       DO i=1, slv_w
          a(i,slv_s+jsize+j,slv_l+k) = recvbuf_r4_nw(slv_w*slv_n*(k-1)+slv_w*(j-1)+i)
       END DO
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, ksize
       DO j=1, slv_n
       DO i=1, slv_w
          a(i,slv_s+jsize+j,slv_l+k) = fill
       END DO
       END DO
       END DO
    END IF

!$OMP SECTION
!$OMP CRITICAL
    CALL mpi_wait(recvreq(tag_nw), MPI_STATUS_IGNORE, ierr)
!$OMP END CRITICAL
    IF (rank_se /= MPI_PROC_NULL) THEN
       DO k=1, ksize
       DO j=1, slv_s
       DO i=1, slv_e
          a(slv_w+isize+i,j,slv_l+k) = recvbuf_r4_se(slv_e*slv_s*(k-1)+slv_e*(j-1)+i)
       END DO
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, ksize
       DO j=1, slv_s
       DO i=1, slv_e
          a(slv_w+isize+i,j,slv_l+k) = fill
       END DO
       END DO
       END DO
    END IF

!$OMP SECTION
!$OMP CRITICAL
    CALL mpi_wait(recvreq(tag_ne), MPI_STATUS_IGNORE, ierr)
!$OMP END CRITICAL
    IF (rank_sw /= MPI_PROC_NULL) THEN
       DO k=1, ksize
       DO j=1, slv_s
       DO i=1, slv_w
          a(i, j,slv_l+k) = recvbuf_r4_sw(slv_w*slv_s*(k-1)+slv_w*(j-1)+i)
       END DO
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, ksize
       DO j=1, slv_s
       DO i=1, slv_w
          a(i, j,slv_l+k) = fill
       END DO
       END DO
       END DO
    END IF

!$OMP SECTION
!$OMP CRITICAL
    CALL mpi_wait(recvreq(tag_w), MPI_STATUS_IGNORE, ierr)
!$OMP END CRITICAL
    IF (rank_e /= MPI_PROC_NULL) THEN
       DO k=1, ksize
       DO j=1, jsize
       DO i=1, slv_e
          a(slv_w+isize+i,slv_s+j,slv_l+k) = recvbuf_r4_e(slv_e*jsize*(k-1)+slv_e*(j-1)+i)
       END DO
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, ksize
       DO j=1, jsize
       DO i=1, slv_e
          a(slv_w+isize+i,slv_s+j,slv_l+k) = fill
       END DO
       END DO
       END DO
    END IF

!$OMP SECTION
!$OMP CRITICAL
    CALL mpi_wait(recvreq(tag_e), MPI_STATUS_IGNORE, ierr)
!$OMP END CRITICAL
    IF (rank_w /= MPI_PROC_NULL) THEN
       DO k=1, ksize
       DO j=1, jsize
       DO i=1, slv_w
          a(i,slv_s+j,slv_l+k) = recvbuf_r4_w(slv_w*jsize*(k-1)+slv_w*(j-1)+i)
       END DO
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, ksize
       DO j=1, jsize
       DO i=1, slv_w
          a(i,slv_s+j,slv_l+k) = fill
       END DO
       END DO
       END DO
    END IF

!$OMP SECTION
!$OMP CRITICAL
    CALL mpi_wait(recvreq(tag_s), MPI_STATUS_IGNORE, ierr)
!$OMP END CRITICAL
    IF (rank_n /= MPI_PROC_NULL) THEN
       DO k=1, ksize
       DO j=1, slv_n
       DO i=1, isize
          a(slv_w+i,slv_s+jsize+j,slv_l+k) = recvbuf_r4_n(isize*slv_n*(k-1)+isize*(j-1)+i)
       END DO
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, ksize
       DO j=1, slv_n
       DO i=1, isize
          a(slv_w+i,slv_s+jsize+j,slv_l+k) = fill
       END DO
       END DO
       END DO
    END IF

!$OMP SECTION
!$OMP CRITICAL
    CALL mpi_wait(recvreq(tag_n), MPI_STATUS_IGNORE, ierr)
!$OMP END CRITICAL
    IF (rank_s /= MPI_PROC_NULL) THEN
       DO k=1, ksize
       DO j=1, slv_s
       DO i=1, isize
          a(slv_w+i,j,slv_l+k) = recvbuf_r4_s(isize*slv_s*(k-1)+isize*(j-1)+i)
       END DO
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, ksize
       DO j=1, slv_s
       DO i=1, isize
          a(slv_w+i,j,slv_l+k) = fill
       END DO
       END DO
       END DO
    END IF

#ifdef PARALLEL3D
!$OMP SECTION
!$OMP CRITICAL
    CALL mpi_wait(recvreq(tag_lsw), MPI_STATUS_IGNORE, ierr)
!$OMP END CRITICAL
    IF (rank_une /= MPI_PROC_NULL) THEN
       DO k=1, slv_u
       DO j=1, slv_n
       DO i=1, slv_e
          a(slv_w+isize+i,slv_s+jsize+j,slv_l+ksize+k) = recvbuf_r4_une(slv_e*slv_n*(k-1)+slv_e*(j-1)+i)
       END DO
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, slv_u
       DO j=1, slv_n
       DO i=1, slv_e
          a(slv_w+isize+i,slv_s+jsize+j,slv_l+ksize+k) = fill
       END DO
       END DO
       END DO
    END IF

!$OMP SECTION
!$OMP CRITICAL
    CALL mpi_wait(recvreq(tag_lse), MPI_STATUS_IGNORE, ierr)
!$OMP END CRITICAL
    IF (rank_unw /= MPI_PROC_NULL) THEN
       DO k=1, slv_u
       DO j=1, slv_n
       DO i=1, slv_w
          a(i,slv_s+jsize+j,slv_l+ksize+k) = recvbuf_r4_unw(slv_w*slv_n*(k-1)+slv_w*(j-1)+i)
       END DO
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, slv_u
       DO j=1, slv_n
       DO i=1, slv_w
          a(i,slv_s+jsize+j,slv_l+ksize+k) = fill
       END DO
       END DO
       END DO
    END IF

!$OMP SECTION
!$OMP CRITICAL
    CALL mpi_wait(recvreq(tag_lnw), MPI_STATUS_IGNORE, ierr)
!$OMP END CRITICAL
    IF (rank_use /= MPI_PROC_NULL) THEN
       DO k=1, slv_u
       DO j=1, slv_s
       DO i=1, slv_e
          a(slv_w+isize+i,j,slv_l+ksize+k) = recvbuf_r4_use(slv_e*slv_s*(k-1)+slv_e*(j-1)+i)
       END DO
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, slv_u
       DO j=1, slv_s
       DO i=1, slv_e
          a(slv_w+isize+i,j,slv_l+ksize+k) = fill
       END DO
       END DO
       END DO
    END IF

!$OMP SECTION
!$OMP CRITICAL
    CALL mpi_wait(recvreq(tag_lne), MPI_STATUS_IGNORE, ierr)
!$OMP END CRITICAL
    IF (rank_usw /= MPI_PROC_NULL) THEN
       DO k=1, slv_u
       DO j=1, slv_s
       DO i=1, slv_w
          a(i, j,slv_l+ksize+k) = recvbuf_r4_usw(slv_w*slv_s*(k-1)+slv_w*(j-1)+i)
       END DO
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, slv_u
       DO j=1, slv_s
       DO i=1, slv_w
          a(i, j,slv_l+ksize+k) = fill
       END DO
       END DO
       END DO
    END IF

!$OMP SECTION
!$OMP CRITICAL
    CALL mpi_wait(recvreq(tag_usw), MPI_STATUS_IGNORE, ierr)
!$OMP END CRITICAL
    IF (rank_lne /= MPI_PROC_NULL) THEN
       DO k=1, slv_l
       DO j=1, slv_n
       DO i=1, slv_e
          a(slv_w+isize+i,slv_s+jsize+j,k) = recvbuf_r4_lne(slv_e*slv_n*(k-1)+slv_e*(j-1)+i)
       END DO
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, slv_l
       DO j=1, slv_n
       DO i=1, slv_e
          a(slv_w+isize+i,slv_s+jsize+j,k) = fill
       END DO
       END DO
       END DO
    END IF

!$OMP SECTION
!$OMP CRITICAL
    CALL mpi_wait(recvreq(tag_use), MPI_STATUS_IGNORE, ierr)
!$OMP END CRITICAL
    IF (rank_lnw /= MPI_PROC_NULL) THEN
       DO k=1, slv_l
       DO j=1, slv_n
       DO i=1, slv_w
          a(i,slv_s+jsize+j,k) = recvbuf_r4_lnw(slv_w*slv_n*(k-1)+slv_w*(j-1)+i)
       END DO
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, slv_l
       DO j=1, slv_n
       DO i=1, slv_w
          a(i,slv_s+jsize+j,k) = fill
       END DO
       END DO
       END DO
    END IF

!$OMP SECTION
!$OMP CRITICAL
    CALL mpi_wait(recvreq(tag_unw), MPI_STATUS_IGNORE, ierr)
!$OMP END CRITICAL
    IF (rank_lse /= MPI_PROC_NULL) THEN
       DO k=1, slv_l
       DO j=1, slv_s
       DO i=1, slv_e
          a(slv_w+isize+i,j,k) = recvbuf_r4_lse(slv_e*slv_s*(k-1)+slv_e*(j-1)+i)
       END DO
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, slv_l
       DO j=1, slv_s
       DO i=1, slv_e
          a(slv_w+isize+i,j,k) = fill
       END DO
       END DO
       END DO
    END IF

!$OMP SECTION
!$OMP CRITICAL
    CALL mpi_wait(recvreq(tag_une), MPI_STATUS_IGNORE, ierr)
!$OMP END CRITICAL
    IF (rank_lsw /= MPI_PROC_NULL) THEN
       DO k=1, slv_l
       DO j=1, slv_s
       DO i=1, slv_w
          a(i,j,k) = recvbuf_r4_lsw(slv_w*slv_s*(k-1)+slv_w*(j-1)+i)
       END DO
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, slv_l
       DO j=1, slv_s
       DO i=1, slv_w
          a(i,j,k) = fill
       END DO
       END DO
       END DO
    END IF

!$OMP SECTION
!$OMP CRITICAL
    CALL mpi_wait(recvreq(tag_lw), MPI_STATUS_IGNORE, ierr)
!$OMP END CRITICAL
    IF (rank_ue /= MPI_PROC_NULL) THEN
       DO k=1, slv_u
       DO j=1, jsize
       DO i=1, slv_e
          a(slv_w+isize+i,slv_s+j,slv_l+ksize+k) = recvbuf_r4_ue(slv_e*jsize*(k-1)+slv_e*(j-1)+i)
       END DO
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, slv_u
       DO j=1, jsize
       DO i=1, slv_e
          a(slv_w+isize+i,slv_s+j,slv_l+ksize+k) = fill
       END DO
       END DO
       END DO
    END IF

!$OMP SECTION
!$OMP CRITICAL
    CALL mpi_wait(recvreq(tag_le), MPI_STATUS_IGNORE, ierr)
!$OMP END CRITICAL
    IF (rank_uw /= MPI_PROC_NULL) THEN
       DO k=1, slv_u
       DO j=1, jsize
       DO i=1, slv_w
          a(i,slv_s+j,slv_l+ksize+k) = recvbuf_r4_uw(slv_w*jsize*(k-1)+slv_w*(j-1)+i)
       END DO
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, slv_u
       DO j=1, jsize
       DO i=1, slv_w
          a(i,slv_s+j,slv_l+ksize+k) = fill
       END DO
       END DO
       END DO
    END IF

!$OMP SECTION
!$OMP CRITICAL
    CALL mpi_wait(recvreq(tag_ls), MPI_STATUS_IGNORE, ierr)
!$OMP END CRITICAL
    IF (rank_un /= MPI_PROC_NULL) THEN
       DO k=1, slv_u
       DO j=1, slv_n
       DO i=1, isize
          a(slv_w+i,slv_s+jsize+j,slv_l+ksize+k) = recvbuf_r4_un(isize*slv_n*(k-1)+isize*(j-1)+i)
       END DO
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, slv_u
       DO j=1, slv_n
       DO i=1, isize
          a(slv_w+i,slv_s+jsize+j,slv_l+ksize+k) = fill
       END DO
       END DO
       END DO
    END IF

!$OMP SECTION
!$OMP CRITICAL
    CALL mpi_wait(recvreq(tag_ln), MPI_STATUS_IGNORE, ierr)
!$OMP END CRITICAL
    IF (rank_us /= MPI_PROC_NULL) THEN
       DO k=1, slv_u
       DO j=1, slv_s
       DO i=1, isize
          a(slv_w+i,j,slv_l+ksize+k) = recvbuf_r4_us(isize*slv_s*(k-1)+isize*(j-1)+i)
       END DO
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, slv_u
       DO j=1, slv_s
       DO i=1, isize
          a(slv_w+i,j,slv_l+ksize+k) = fill
       END DO
       END DO
       END DO
    END IF

!$OMP SECTION
!$OMP CRITICAL
    CALL mpi_wait(recvreq(tag_uw), MPI_STATUS_IGNORE, ierr)
!$OMP END CRITICAL
    IF (rank_le /= MPI_PROC_NULL) THEN
       DO k=1, slv_l
       DO j=1, jsize
       DO i=1, slv_e
          a(slv_w+isize+i,slv_s+j,k) = recvbuf_r4_le(slv_e*jsize*(k-1)+slv_e*(j-1)+i)
       END DO
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, slv_l
       DO j=1, jsize
       DO i=1, slv_e
          a(slv_w+isize+i,slv_s+j,k) = fill
       END DO
       END DO
       END DO
    END IF

!$OMP SECTION
!$OMP CRITICAL
    CALL mpi_wait(recvreq(tag_ue), MPI_STATUS_IGNORE, ierr)
!$OMP END CRITICAL
    IF (rank_lw /= MPI_PROC_NULL) THEN
       DO k=1, slv_l
       DO j=1, jsize
       DO i=1, slv_w
          a(i,slv_s+j,k) = recvbuf_r4_lw(slv_w*jsize*(k-1)+slv_w*(j-1)+i)
       END DO
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, slv_l
       DO j=1, jsize
       DO i=1, slv_w
          a(i,slv_s+j,k) = fill
       END DO
       END DO
       END DO
    END IF

!$OMP SECTION
!$OMP CRITICAL
    CALL mpi_wait(recvreq(tag_us), MPI_STATUS_IGNORE, ierr)
!$OMP END CRITICAL
    IF (rank_ln /= MPI_PROC_NULL) THEN
       DO k=1, slv_l
       DO j=1, slv_n
       DO i=1, isize
          a(slv_w+i,slv_s+jsize+j,k) = recvbuf_r4_ln(isize*slv_n*(k-1)+isize*(j-1)+i)
       END DO
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, slv_l
       DO j=1, slv_n
       DO i=1, isize
          a(slv_w+i,slv_s+jsize+j,k) = fill
       END DO
       END DO
       END DO
    END IF

!$OMP SECTION
!$OMP CRITICAL
    CALL mpi_wait(recvreq(tag_un), MPI_STATUS_IGNORE, ierr)
!$OMP END CRITICAL
    IF (rank_ls /= MPI_PROC_NULL) THEN
       DO k=1, slv_l
       DO j=1, slv_s
       DO i=1, isize
          a(slv_w+i,j,k) = recvbuf_r4_ls(isize*slv_s*(k-1)+isize*(j-1)+i)
       END DO
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, slv_l
       DO j=1, slv_s
       DO i=1, isize
          a(slv_w+i,j,k) = fill
       END DO
       END DO
       END DO
    END IF

!$OMP SECTION
!$OMP CRITICAL
    CALL mpi_wait(recvreq(tag_l), MPI_STATUS_IGNORE, ierr)
!$OMP END CRITICAL
    IF (rank_u /= MPI_PROC_NULL) THEN
       DO k=1, slv_u
       DO j=1, jsize
       DO i=1, isize
          a(slv_w+i,slv_s+j,slv_l+ksize+k) = recvbuf_r4_u(isize*jsize*(k-1)+isize*(j-1)+i)
       END DO
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, slv_u
       DO j=1, jsize
       DO i=1, isize
          a(slv_w+i,slv_s+j,slv_l+ksize+k) = fill
       END DO
       END DO
       END DO
    END IF

!$OMP SECTION
!$OMP CRITICAL
    CALL mpi_wait(recvreq(tag_u), MPI_STATUS_IGNORE, ierr)
!$OMP END CRITICAL
    IF (rank_l /= MPI_PROC_NULL) THEN
       DO k=1, slv_l
       DO j=1, jsize
       DO i=1, isize
          a(slv_w+i,slv_s+j,k) = recvbuf_r4_l(isize*jsize*(k-1)+isize*(j-1)+i)
       END DO
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, slv_l
       DO j=1, jsize
       DO i=1, isize
          a(slv_w+i,slv_s+j,k) = fill
       END DO
       END DO
       END DO
    END IF
!$OMP END SECTIONS
!$OMP END PARALLEL

#else
!$OMP END SECTIONS
!$OMP END PARALLEL

    IF (present(fill)) THEN
       DO k=1, slv_l
          a(:,:,k) = fill
       END DO
       DO k=1, slv_u
          a(:,:,slv_l+ksize+k) = fill
       END DO
    END IF
#endif

    CALL mpi_waitall(n_sendrecv, sendreq, MPI_STATUSES_IGNORE, ierr)

    IF (rank_tp(0,0)/=MPI_PROC_NULL) THEN
       DO k=1, ksize+slv_l+slv_u
       DO j=1, slv_n
       DO i=1, isize+slv_w+slv_e
          sendbuf_r4_n((isize+slv_w+slv_e)*slv_n*(k-1)+(isize+slv_w+slv_e)*(j-1)+i) = a(isize+slv_w+slv_e+1-i,jsize+slv_n+1-j,k)
       END DO
       END DO
       END DO

       CALL mpi_sendrecv(sendbuf_r4_n, (isize+slv_w+slv_e)*slv_n*(ksize+slv_l+slv_u), MPI_REAL4, rank_tp(0,0), 30, &
                         recvbuf_r4_n, (isize+slv_w+slv_e)*slv_n*(ksize+slv_l+slv_u), MPI_REAL4, rank_tp(0,0), 30, comm, MPI_STATUS_IGNORE, ierr)

       DO k=1, ksize+slv_l+slv_u
       DO j=1, slv_n
       DO i=1, isize+slv_w+slv_e
          a(i,jsize+slv_s+j,k) = recvbuf_r4_n((isize+slv_w+slv_e)*slv_n*(k-1)+(isize+slv_w+slv_e)*(j-1)+i)
       END DO
       END DO
       END DO
    END IF

  END SUBROUTINE parallel_nb
#endif

  END SUBROUTINE update_boundary_3d_r4

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE update_boundary_3d_r8(a, fill, method)
    REAL(8), INTENT(INOUT) :: a(:,:,:)
    REAL(8), INTENT(IN), OPTIONAL :: fill
    CHARACTER(*), INTENT(IN), OPTIONAL :: method

    INTEGER :: nx, ny, nz
    INTEGER :: slv_e, slv_w, slv_n, slv_s, slv_u, slv_l

    CHARACTER(16) :: method_

    nx = size(a,1) - isize
    ny = size(a,2) - jsize
    nz = size(a,3) - ksize

    CALL assert(nx >= 0 .AND. nx <= 5, "unsported dimension in UPDATE_BOUNDARY_3D")
    CALL assert(ny >= 0 .AND. ny <= 5, "unsported dimension in UPDATE_BOUNDARY_3D")
    CALL assert(nz >= 0 .AND. nz <= 5, "unsported dimension in UPDATE_BOUNDARY_3D")

    slv_e = nx/2
    slv_n = ny/2
    slv_u = nz/2
    slv_w = (nx+1)/2
    slv_s = (ny+1)/2
    slv_l = (nz+1)/2

#ifdef PARALLEL_MPI
    method_ = sendrecv_method
    IF (present(method)) method_ = trim(method)

    SELECT CASE (trim(method_))
    CASE ('STD')
       CALL parallel_std
    CASE ('NB')
       CALL parallel_nb
    CASE DEFAULT
       CALL assert(.FALSE., "unsupported SENDRECV_METHOD")
    END SELECT
#else
    CALL single
#endif

  CONTAINS

  SUBROUTINE single
    INTEGER :: i, j, k

    IF (cycle_x) THEN
       DO i=1, slv_w
          a(i,:,:) = a(dimx+i,:,:)
       END DO
       DO i=1, slv_e
          a(slv_w+dimx+i,:,:) = a(slv_w+i,:,:)
       END DO
    ELSE IF (present(fill)) THEN
       DO i=1, slv_w
          a(i,:,:) = fill
       END DO
       DO i=1, slv_e
          a(slv_w+dimx+i,:,:) = fill
       END DO
    END IF

    IF (cycle_y) THEN
       DO j=1, slv_s
          a(:,j,:) = a(:,dimy+j,:)
       END DO
       DO j=1, slv_n
          a(:,slv_s+dimy+j,:) = a(:,slv_s+j,:)
       END DO
    ELSE IF (present(fill)) THEN
       DO j=1, slv_s
          a(:,j,:) = fill
       END DO
       DO j=1, slv_n
          a(:,slv_s+dimy+j,:) = fill
       END DO
    END IF

    IF (cycle_z) THEN
       DO k=1, slv_l
          a(:,:,k) = a(:,:,dimz+k)
       END DO
       DO k=1, slv_u
          a(:,:,slv_l+dimz+k) = a(:,:,slv_l+k)
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, slv_l
          a(:,:,k) = fill
       END DO
       DO k=1, slv_u
          a(:,:,slv_l+dimz+k) = fill
       END DO
    END IF

    IF (tripolar) THEN
       DO j=1, slv_n
       DO i=1, dimx+slv_w+slv_e
          a(i,dimy+slv_s+j,:) = a(dimx+slv_w+slv_e+1-i, dimy+slv_n+1-j,:)
       END DO
       END DO
    END IF

  END SUBROUTINE single

#ifdef PARALLEL_MPI
  SUBROUTINE parallel_std
    INTEGER :: ierr
    INTEGER :: i, j, k

    DO k=1, ksize
    DO j=1, jsize
       DO i=1, slv_e
          sendbuf_r8_w(jsize*slv_e*(k-1)+slv_e*(j-1)+i) = a(slv_w+i,slv_s+j,slv_l+k)
       END DO

       DO i=1, slv_w
          sendbuf_r8_e(jsize*slv_w*(k-1)+slv_w*(j-1)+i) = a(isize+i,slv_s+j,slv_l+k)
       END DO
    END DO
    END DO

    IF (nx == 1) THEN
       CALL mpi_startall(2, mpireq_3dx_r8(:,nx), ierr)
    ELSE IF (nx > 1) THEN
       CALL mpi_startall(4, mpireq_3dx_r8(:,nx), ierr)
    END IF
    CALL mpi_waitall(4, mpireq_3dx_r8(:,nx), MPI_STATUSES_IGNORE, ierr)

    IF (rank_w /= MPI_PROC_NULL) THEN
       DO k=1, ksize
       DO j=1, jsize
       DO i=1, slv_w
          a(i,slv_s+j,slv_l+k) = recvbuf_r8_w(jsize*slv_w*(k-1)+slv_w*(j-1)+i)
       END DO
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, ksize
       DO j=1, jsize
       DO i=1, slv_w
          a(i,slv_s+j,slv_l+k) = fill
       END DO
       END DO
       END DO
    END IF

    IF (rank_e /= MPI_PROC_NULL) THEN
       DO k=1, ksize
       DO j=1, jsize
       DO i=1, slv_e
          a(slv_w+isize+i,slv_s+j,slv_l+k) = recvbuf_r8_e(jsize*slv_e*(k-1)+slv_e*(j-1)+i)
       END DO
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, ksize
       DO j=1, jsize
       DO i=1, slv_e
          a(slv_w+isize+i,slv_s+j,slv_l+k) = fill
       END DO
       END DO
       END DO
    END IF

    DO k=1, ksize
       DO j=1, slv_n
       DO i=1, isize+slv_w+slv_e
          sendbuf_r8_s((isize+slv_w+slv_e)*slv_n*(k-1)+(isize+slv_w+slv_e)*(j-1)+i) = a(i,slv_s+j,slv_l+k)
       END DO
       END DO

       DO j=1, slv_s
       DO i=1, isize+slv_w+slv_e
          sendbuf_r8_n((isize+slv_w+slv_e)*slv_s*(k-1)+(isize+slv_w+slv_e)*(j-1)+i) = a(i,jsize+j,slv_l+k)
       END DO
       END DO
    END DO

    IF (ny == 1) THEN
       CALL mpi_startall(2, mpireq_3dy_r8(:,ny), ierr)
    ELSE IF (ny > 1) THEN
       CALL mpi_startall(4, mpireq_3dy_r8(:,ny), ierr)
    END IF
    CALL mpi_waitall(4, mpireq_3dy_r8(:,ny), MPI_STATUSES_IGNORE, ierr)

    IF (rank_s /= MPI_PROC_NULL) THEN
       DO k=1, ksize
       DO j=1, slv_s
       DO i=1, isize+slv_w+slv_e
          a(i,j,slv_l+k) = recvbuf_r8_s((isize+slv_w+slv_e)*slv_s*(k-1)+(isize+slv_w+slv_e)*(j-1)+i)
       END DO
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, ksize
       DO j=1, slv_s
          a(:,j,slv_l+k) = fill
       END DO
       END DO
    END IF

    IF (rank_n /= MPI_PROC_NULL) THEN
       DO k=1, ksize
       DO j=1, slv_n
       DO i=1, isize+slv_w+slv_e
          a(i,slv_s+jsize+j,slv_l+k) = recvbuf_r8_n((isize+slv_w+slv_e)*slv_n*(k-1)+(isize+slv_w+slv_e)*(j-1)+i)
       END DO
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, ksize
       DO j=1, slv_n
          a(:,slv_s+jsize+j,slv_l+k) = fill
       END DO
       END DO
    END IF

#ifdef PARALLEL3D
    DO k=1, slv_u
       DO j=1, jsize+slv_s+slv_n
       DO i=1, isize+slv_w+slv_e
          sendbuf_r8_l((isize+slv_w+slv_e)*(jsize+slv_s+slv_n)*(k-1)+(isize+slv_w+slv_e)*(j-1)+i) = a(i,j,slv_l+k)
       END DO
       END DO
    END DO

    DO k=1, slv_l
       DO j=1, jsize+slv_s+slv_n
       DO i=1, isize+slv_w+slv_e
          sendbuf_r8_u((isize+slv_w+slv_e)*(jsize+slv_s+slv_n)*(k-1)+(isize+slv_w+slv_e)*(j-1)+i) = a(i,j,ksize+k)
       END DO
       END DO
    END DO

    IF (nz == 1) THEN
       CALL mpi_startall(2, mpireq_3dz_r8(:,nz), ierr)
    ELSE IF (nz > 1) THEN
       CALL mpi_startall(4, mpireq_3dz_r8(:,nz), ierr)
    END IF
    CALL mpi_waitall(4, mpireq_3dz_r8(:,nz), MPI_STATUSES_IGNORE, ierr)

    IF (rank_l /= MPI_PROC_NULL) THEN
       DO k=1, slv_l
       DO j=1, jsize+slv_s+slv_n
       DO i=1, isize+slv_w+slv_e
          a(i,j,k) = recvbuf_r8_l((isize+slv_w+slv_e)*(jsize+slv_s+slv_n)*(k-1)+(isize+slv_w+slv_e)*(j-1)+i)
       END DO
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, slv_l
          a(:,:,k) = fill
       END DO
    END IF

    IF (rank_u /= MPI_PROC_NULL) THEN
       DO k=1, slv_u
       DO j=1, jsize+slv_s+slv_n
       DO i=1, isize+slv_w+slv_e
          a(i,j,slv_l+ksize+k) = recvbuf_r8_u((isize+slv_w+slv_e)*(jsize+slv_s+slv_n)*(k-1)+(isize+slv_e+slv_w)*(j-1)+i)
       END DO
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, slv_u
          a(:,:,slv_l+ksize+k) = fill
       END DO
    END IF
#else
    IF (present(fill)) THEN
       DO k=1, slv_l
          a(:,:,k) = fill
       END DO
       DO k=1, slv_u
          a(:,:,slv_l+ksize+k) = fill
       END DO
    END IF
#endif

    IF (rank_tp(0,0)/=MPI_PROC_NULL) THEN
       DO k=1, ksize+slv_l+slv_u
       DO j=1, slv_n
       DO i=1, isize+slv_w+slv_e
          sendbuf_r8_n((isize+slv_w+slv_e)*slv_n*(k-1)+(isize+slv_w+slv_e)*(j-1)+i) = a(isize+slv_w+slv_e+1-i,jsize+slv_n+1-j,k)
       END DO
       END DO
       END DO

       CALL mpi_sendrecv(sendbuf_r8_n, (isize+slv_w+slv_e)*slv_n*(ksize+slv_l+slv_u), MPI_REAL8, rank_tp(0,0), 30, &
                         recvbuf_r8_n, (isize+slv_w+slv_e)*slv_n*(ksize+slv_l+slv_u), MPI_REAL8, rank_tp(0,0), 30, comm, MPI_STATUS_IGNORE, ierr)

       DO k=1, ksize+slv_l+slv_u
       DO j=1, slv_n
       DO i=1, isize+slv_w+slv_e
          a(i,jsize+slv_s+j,k) = recvbuf_r8_n((isize+slv_w+slv_e)*slv_n*(k-1)+(isize+slv_w+slv_e)*(j-1)+i)
       END DO
       END DO
       END DO
    END IF

  END SUBROUTINE parallel_std

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE parallel_nb
    INTEGER :: sendreq(n_sendrecv), recvreq(n_sendrecv), ierr
    INTEGER :: i, j, k

    sendreq(:) = MPI_REQUEST_NULL
    recvreq(:) = MPI_REQUEST_NULL

!$OMP PARALLEL PRIVATE(i, j, k, ierr)
!$OMP SECTIONS
!$OMP SECTION
    DO k=1, ksize
    DO j=1, slv_n
    DO i=1, slv_e
       sendbuf_r8_sw(slv_e*slv_n*(k-1)+slv_e*(j-1)+i) = a(slv_w+i,slv_s+j,slv_l+k)
    END DO
    END DO
    END DO
!$OMP CRITICAL
    CALL mpi_isend(sendbuf_r8_sw, slv_e*slv_n*ksize, MPI_REAL8, rank_sw, tag_sw, comm, sendreq(tag_sw), ierr)
    CALL mpi_irecv(recvbuf_r8_ne, slv_e*slv_n*ksize, MPI_REAL8, rank_ne, tag_sw, comm, recvreq(tag_sw), ierr)
!$OMP END CRITICAL

!$OMP SECTION
    DO k=1, ksize
    DO j=1, slv_n
    DO i=1, slv_w
       sendbuf_r8_se(slv_w*slv_n*(k-1)+slv_w*(j-1)+i) = a(isize+i,slv_s+j,slv_l+k)
    END DO
    END DO
    END DO
!$OMP CRITICAL
    CALL mpi_isend(sendbuf_r8_se, slv_w*slv_n*ksize, MPI_REAL8, rank_se, tag_se, comm, sendreq(tag_se), ierr)
    CALL mpi_irecv(recvbuf_r8_nw, slv_w*slv_n*ksize, MPI_REAL8, rank_nw, tag_se, comm, recvreq(tag_se), ierr)
!$OMP END CRITICAL

!$OMP SECTION
    DO k=1, ksize
    DO j=1, slv_s
    DO i=1, slv_e
       sendbuf_r8_nw(slv_e*slv_s*(k-1)+slv_e*(j-1)+i) = a(slv_w+i,jsize+j,slv_l+k)
    END DO
    END DO
    END DO
!$OMP CRITICAL
    CALL mpi_isend(sendbuf_r8_nw, slv_e*slv_s*ksize, MPI_REAL8, rank_nw, tag_nw, comm, sendreq(tag_nw), ierr)
    CALL mpi_irecv(recvbuf_r8_se, slv_e*slv_s*ksize, MPI_REAL8, rank_se, tag_nw, comm, recvreq(tag_nw), ierr)
!$OMP END CRITICAL

!$OMP SECTION
    DO k=1, ksize
    DO j=1, slv_s
    DO i=1, slv_w
       sendbuf_r8_ne(slv_w*slv_s*(k-1)+slv_w*(j-1)+i) = a(isize+i,jsize+j,slv_l+k)
    END DO
    END DO
    END DO
!$OMP CRITICAL
    CALL mpi_isend(sendbuf_r8_ne, slv_w*slv_s*ksize, MPI_REAL8, rank_ne, tag_ne, comm, sendreq(tag_ne), ierr)
    CALL mpi_irecv(recvbuf_r8_sw, slv_w*slv_s*ksize, MPI_REAL8, rank_sw, tag_ne, comm, recvreq(tag_ne), ierr)
!$OMP END CRITICAL

!$OMP SECTION
    DO k=1, ksize
    DO j=1, jsize
    DO i=1, slv_e
       sendbuf_r8_w(slv_e*jsize*(k-1)+slv_e*(j-1)+i) = a(slv_w+i,slv_s+j,slv_l+k)
    END DO
    END DO
    END DO
!$OMP CRITICAL
    CALL mpi_isend(sendbuf_r8_w, slv_e*jsize*ksize, MPI_REAL8, rank_w, tag_w, comm, sendreq(tag_w), ierr)
    CALL mpi_irecv(recvbuf_r8_e, slv_e*jsize*ksize, MPI_REAL8, rank_e, tag_w, comm, recvreq(tag_w), ierr)
!$OMP END CRITICAL

!$OMP SECTION
    DO k=1, ksize
    DO j=1, jsize
    DO i=1, slv_w
       sendbuf_r8_e(slv_w*jsize*(k-1)+slv_w*(j-1)+i) = a(isize+i,slv_s+j,slv_l+k)
    END DO
    END DO
    END DO
!$OMP CRITICAL
    CALL mpi_isend(sendbuf_r8_e, slv_w*jsize*ksize, MPI_REAL8, rank_e, tag_e, comm, sendreq(tag_e), ierr)
    CALL mpi_irecv(recvbuf_r8_w, slv_w*jsize*ksize, MPI_REAL8, rank_w, tag_e, comm, recvreq(tag_e), ierr)
!$OMP END CRITICAL

!$OMP SECTION
    DO k=1, ksize
    DO j=1, slv_n
    DO i=1, isize
       sendbuf_r8_s(isize*slv_n*(k-1)+isize*(j-1)+i) = a(slv_w+i,slv_s+j,slv_l+k)
    END DO
    END DO
    END DO
!$OMP CRITICAL
    CALL mpi_isend(sendbuf_r8_s, isize*slv_n*ksize, MPI_REAL8, rank_s, tag_s, comm, sendreq(tag_s), ierr)
    CALL mpi_irecv(recvbuf_r8_n, isize*slv_n*ksize, MPI_REAL8, rank_n, tag_s, comm, recvreq(tag_s), ierr)
!$OMP END CRITICAL

!$OMP SECTION
    DO k=1, ksize
    DO j=1, slv_s
    DO i=1, isize
       sendbuf_r8_n(isize*slv_s*(k-1)+isize*(j-1)+i) = a(slv_w+i,jsize+j,slv_l+k)
    END DO
    END DO
    END DO
!$OMP CRITICAL
    CALL mpi_isend(sendbuf_r8_n, isize*slv_s*ksize, MPI_REAL8, rank_n, tag_n, comm, sendreq(tag_n), ierr)
    CALL mpi_irecv(recvbuf_r8_s, isize*slv_s*ksize, MPI_REAL8, rank_s, tag_n, comm, recvreq(tag_n), ierr)
!$OMP END CRITICAL

#ifdef PARALLEL3D
!$OMP SECTION
    DO k=1, slv_u
    DO j=1, slv_n
    DO i=1, slv_e
       sendbuf_r8_lsw(slv_e*slv_n*(k-1)+slv_e*(j-1)+i) = a(slv_w+i,slv_s+j,slv_l+k)
    END DO
    END DO
    END DO
!$OMP CRITICAL
    CALL mpi_isend(sendbuf_r8_lsw, slv_e*slv_n*slv_u, MPI_REAL8, rank_lsw, tag_lsw, comm, sendreq(tag_lsw), ierr)
    CALL mpi_irecv(recvbuf_r8_une, slv_e*slv_n*slv_u, MPI_REAL8, rank_une, tag_lsw, comm, recvreq(tag_lsw), ierr)
!$OMP END CRITICAL

!$OMP SECTION
    DO k=1, slv_u
    DO j=1, slv_n
    DO i=1, slv_w
       sendbuf_r8_lse(slv_w*slv_n*(k-1)+slv_w*(j-1)+i) = a(isize+i,slv_s+j,slv_l+k)
    END DO
    END DO
    END DO
!$OMP CRITICAL
    CALL mpi_isend(sendbuf_r8_lse, slv_w*slv_n*slv_u, MPI_REAL8, rank_lse, tag_lse, comm, sendreq(tag_lse), ierr)
    CALL mpi_irecv(recvbuf_r8_unw, slv_w*slv_n*slv_u, MPI_REAL8, rank_unw, tag_lse, comm, recvreq(tag_lse), ierr)
!$OMP END CRITICAL

!$OMP SECTION
    DO k=1, slv_u
    DO j=1, slv_s
    DO i=1, slv_e
       sendbuf_r8_lnw(slv_e*slv_s*(k-1)+slv_e*(j-1)+i) = a(slv_w+i,jsize+j,slv_l+k)
    END DO
    END DO
    END DO
!$OMP CRITICAL
    CALL mpi_isend(sendbuf_r8_lnw, slv_e*slv_s*slv_u, MPI_REAL8, rank_lnw, tag_lnw, comm, sendreq(tag_lnw), ierr)
    CALL mpi_irecv(recvbuf_r8_use, slv_e*slv_s*slv_u, MPI_REAL8, rank_use, tag_lnw, comm, recvreq(tag_lnw), ierr)
!$OMP END CRITICAL

!$OMP SECTION
    DO k=1, slv_u
    DO j=1, slv_s
    DO i=1, slv_w
       sendbuf_r8_lne(slv_w*slv_s*(k-1)+slv_w*(j-1)+i) = a(isize+i,jsize+j,slv_l+k)
    END DO
    END DO
    END DO
!$OMP CRITICAL
    CALL mpi_isend(sendbuf_r8_lne, slv_w*slv_s*slv_u, MPI_REAL8, rank_lne, tag_lne, comm, sendreq(tag_lne), ierr)
    CALL mpi_irecv(recvbuf_r8_usw, slv_w*slv_s*slv_u, MPI_REAL8, rank_usw, tag_lne, comm, recvreq(tag_lne), ierr)
!$OMP END CRITICAL

!$OMP SECTION
    DO k=1, slv_l
    DO j=1, slv_n
    DO i=1, slv_e
       sendbuf_r8_usw(slv_e*slv_n*(k-1)+slv_e*(j-1)+i) = a(slv_w+i,slv_s+j,ksize+k)
    END DO
    END DO
    END DO
!$OMP CRITICAL
    CALL mpi_isend(sendbuf_r8_usw, slv_e*slv_n*slv_l, MPI_REAL8, rank_usw, tag_usw, comm, sendreq(tag_usw), ierr)
    CALL mpi_irecv(recvbuf_r8_lne, slv_e*slv_n*slv_l, MPI_REAL8, rank_lne, tag_usw, comm, recvreq(tag_usw), ierr)
!$OMP END CRITICAL

!$OMP SECTION
    DO k=1, slv_l
    DO j=1, slv_n
    DO i=1, slv_w
       sendbuf_r8_use(slv_w*slv_n*(k-1)+slv_w*(j-1)+i) = a(isize+i,slv_s+j,ksize+k)
    END DO
    END DO
    END DO
!$OMP CRITICAL
    CALL mpi_isend(sendbuf_r8_use, slv_w*slv_n*slv_l, MPI_REAL8, rank_use, tag_use, comm, sendreq(tag_use), ierr)
    CALL mpi_irecv(recvbuf_r8_lnw, slv_w*slv_n*slv_l, MPI_REAL8, rank_lnw, tag_use, comm, recvreq(tag_use), ierr)
!$OMP END CRITICAL

!$OMP SECTION
    DO k=1, slv_l
    DO j=1, slv_s
    DO i=1, slv_e
       sendbuf_r8_unw(slv_e*slv_s*(k-1)+slv_e*(j-1)+i) = a(slv_w+i,jsize+j,ksize+k)
    END DO
    END DO
    END DO
!$OMP CRITICAL
    CALL mpi_isend(sendbuf_r8_unw, slv_e*slv_s*slv_l, MPI_REAL8, rank_unw, tag_unw, comm, sendreq(tag_unw), ierr)
    CALL mpi_irecv(recvbuf_r8_lse, slv_e*slv_s*slv_l, MPI_REAL8, rank_lse, tag_unw, comm, recvreq(tag_unw), ierr)
!$OMP END CRITICAL

!$OMP SECTION
    DO k=1, slv_l
    DO j=1, slv_s
    DO i=1, slv_w
       sendbuf_r8_une(slv_w*slv_s*(k-1)+slv_w*(j-1)+i) = a(isize+i,jsize+j,ksize+k)
    END DO
    END DO
    END DO
!$OMP CRITICAL
    CALL mpi_isend(sendbuf_r8_une, slv_w*slv_s*slv_l, MPI_REAL8, rank_une, tag_une, comm, sendreq(tag_une), ierr)
    CALL mpi_irecv(recvbuf_r8_lsw, slv_w*slv_s*slv_l, MPI_REAL8, rank_lsw, tag_une, comm, recvreq(tag_une), ierr)
!$OMP END CRITICAL

!$OMP SECTION
    DO k=1, slv_u
    DO j=1, jsize
    DO i=1, slv_e
       sendbuf_r8_lw(slv_e*jsize*(k-1)+slv_e*(j-1)+i) = a(slv_w+i,slv_s+j,slv_l+k)
    END DO
    END DO
    END DO
!$OMP CRITICAL
    CALL mpi_isend(sendbuf_r8_lw, slv_e*jsize*slv_u, MPI_REAL8, rank_lw, tag_lw, comm, sendreq(tag_lw), ierr)
    CALL mpi_irecv(recvbuf_r8_ue, slv_e*jsize*slv_u, MPI_REAL8, rank_ue, tag_lw, comm, recvreq(tag_lw), ierr)
!$OMP END CRITICAL

!$OMP SECTION
    DO k=1, slv_u
    DO j=1, jsize
    DO i=1, slv_w
       sendbuf_r8_le(slv_w*jsize*(k-1)+slv_w*(j-1)+i) = a(isize+i,slv_s+j,slv_l+k)
    END DO
    END DO
    END DO
!$OMP CRITICAL
    CALL mpi_isend(sendbuf_r8_le, slv_w*jsize*slv_u, MPI_REAL8, rank_le, tag_le, comm, sendreq(tag_le), ierr)
    CALL mpi_irecv(recvbuf_r8_uw, slv_w*jsize*slv_u, MPI_REAL8, rank_uw, tag_le, comm, recvreq(tag_le), ierr)
!$OMP END CRITICAL

!$OMP SECTION
    DO k=1, slv_u
    DO j=1, slv_n
    DO i=1, isize
       sendbuf_r8_ls(isize*slv_n*(k-1)+isize*(j-1)+i) = a(slv_w+i,slv_s+j,slv_l+k)
    END DO
    END DO
    END DO
!$OMP CRITICAL
    CALL mpi_isend(sendbuf_r8_ls, isize*slv_n*slv_u, MPI_REAL8, rank_ls, tag_ls, comm, sendreq(tag_ls), ierr)
    CALL mpi_irecv(recvbuf_r8_un, isize*slv_n*slv_u, MPI_REAL8, rank_un, tag_ls, comm, recvreq(tag_ls), ierr)
!$OMP END CRITICAL

!$OMP SECTION
    DO k=1, slv_u
    DO j=1, slv_s
    DO i=1, isize
       sendbuf_r8_ln(isize*slv_s*(k-1)+isize*(j-1)+i) = a(slv_w+i,jsize+j,slv_l+k)
    END DO
    END DO
    END DO
!$OMP CRITICAL
    CALL mpi_isend(sendbuf_r8_ln, isize*slv_s*slv_u, MPI_REAL8, rank_ln, tag_ln, comm, sendreq(tag_ln), ierr)
    CALL mpi_irecv(recvbuf_r8_us, isize*slv_s*slv_u, MPI_REAL8, rank_us, tag_ln, comm, recvreq(tag_ln), ierr)
!$OMP END CRITICAL

!$OMP SECTION
    DO k=1, slv_l
    DO j=1, jsize
    DO i=1, slv_e
       sendbuf_r8_uw(slv_e*jsize*(k-1)+slv_e*(j-1)+i) = a(slv_w+i,slv_s+j,ksize+k)
    END DO
    END DO
    END DO
!$OMP CRITICAL
    CALL mpi_isend(sendbuf_r8_uw, slv_e*jsize*slv_l, MPI_REAL8, rank_uw, tag_uw, comm, sendreq(tag_uw), ierr)
    CALL mpi_irecv(recvbuf_r8_le, slv_e*jsize*slv_l, MPI_REAL8, rank_le, tag_uw, comm, recvreq(tag_uw), ierr)
!$OMP END CRITICAL

!$OMP SECTION
    DO k=1, slv_l
    DO j=1, jsize
    DO i=1, slv_w
       sendbuf_r8_ue(slv_w*jsize*(k-1)+slv_w*(j-1)+i) = a(isize+i,slv_s+j,ksize+k)
    END DO
    END DO
    END DO
!$OMP CRITICAL
    CALL mpi_isend(sendbuf_r8_ue, slv_w*jsize*slv_l, MPI_REAL8, rank_ue, tag_ue, comm, sendreq(tag_ue), ierr)
    CALL mpi_irecv(recvbuf_r8_lw, slv_w*jsize*slv_l, MPI_REAL8, rank_lw, tag_ue, comm, recvreq(tag_ue), ierr)
!$OMP END CRITICAL

!$OMP SECTION
    DO k=1, slv_l
    DO j=1, slv_n
    DO i=1, isize
       sendbuf_r8_us(isize*slv_n*(k-1)+isize*(j-1)+i) = a(slv_w+i,slv_s+j,ksize+k)
    END DO
    END DO
    END DO
!$OMP CRITICAL
    CALL mpi_isend(sendbuf_r8_us, isize*slv_n*slv_l, MPI_REAL8, rank_us, tag_us, comm, sendreq(tag_us), ierr)
    CALL mpi_irecv(recvbuf_r8_ln, isize*slv_n*slv_l, MPI_REAL8, rank_ln, tag_us, comm, recvreq(tag_us), ierr)
!$OMP END CRITICAL

!$OMP SECTION
    DO k=1, slv_l
    DO j=1, slv_s
    DO i=1, isize
       sendbuf_r8_un(isize*slv_s*(k-1)+isize*(j-1)+i) = a(slv_w+i,jsize+j,ksize+k)
    END DO
    END DO
    END DO
!$OMP CRITICAL
    CALL mpi_isend(sendbuf_r8_un, isize*slv_s*slv_l, MPI_REAL8, rank_un, tag_un, comm, sendreq(tag_un), ierr)
    CALL mpi_irecv(recvbuf_r8_ls, isize*slv_s*slv_l, MPI_REAL8, rank_ls, tag_un, comm, recvreq(tag_un), ierr)
!$OMP END CRITICAL

!$OMP SECTION
    DO k=1, slv_u
    DO j=1, jsize
    DO i=1, isize
       sendbuf_r8_l(isize*jsize*(k-1)+isize*(j-1)+i) = a(slv_w+i,slv_s+j,slv_l+k)
    END DO
    END DO
    END DO
!$OMP CRITICAL
    CALL mpi_isend(sendbuf_r8_l, isize*jsize*slv_u, MPI_REAL8, rank_l, tag_l, comm, sendreq(tag_l), ierr)
    CALL mpi_irecv(recvbuf_r8_u, isize*jsize*slv_u, MPI_REAL8, rank_u, tag_l, comm, recvreq(tag_l), ierr)
!$OMP END CRITICAL

!$OMP SECTION
    DO k=1, slv_l
    DO j=1, jsize
    DO i=1, isize
       sendbuf_r8_u(isize*jsize*(k-1)+isize*(j-1)+i) = a(slv_w+i,slv_s+j,ksize+k)
    END DO
    END DO
    END DO
!$OMP CRITICAL
    CALL mpi_isend(sendbuf_r8_u, isize*jsize*slv_l, MPI_REAL8, rank_u, tag_u, comm, sendreq(tag_u), ierr)
    CALL mpi_irecv(recvbuf_r8_l, isize*jsize*slv_l, MPI_REAL8, rank_l, tag_u, comm, recvreq(tag_u), ierr)
!$OMP END CRITICAL
#endif
!$OMP END SECTIONS


!$OMP SECTIONS
!$OMP SECTION
!$OMP CRITICAL
    CALL mpi_wait(recvreq(tag_sw), MPI_STATUS_IGNORE, ierr)
!$OMP END CRITICAL
    IF (rank_ne /= MPI_PROC_NULL) THEN
       DO k=1, ksize
       DO j=1, slv_n
       DO i=1, slv_e
          a(slv_w+isize+i,slv_s+jsize+j,slv_l+k) = recvbuf_r8_ne(slv_e*slv_n*(k-1)+slv_e*(j-1)+i)
       END DO
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, ksize
       DO j=1, slv_n
       DO i=1, slv_e
          a(slv_w+isize+i,slv_s+jsize+j,slv_l+k) = fill
       END DO
       END DO
       END DO
    END IF

!$OMP SECTION
!$OMP CRITICAL
    CALL mpi_wait(recvreq(tag_se), MPI_STATUS_IGNORE, ierr)
!$OMP END CRITICAL
    IF (rank_nw /= MPI_PROC_NULL) THEN
       DO k=1, ksize
       DO j=1, slv_n
       DO i=1, slv_w
          a(i,slv_s+jsize+j,slv_l+k) = recvbuf_r8_nw(slv_w*slv_n*(k-1)+slv_w*(j-1)+i)
       END DO
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, ksize
       DO j=1, slv_n
       DO i=1, slv_w
          a(i,slv_s+jsize+j,slv_l+k) = fill
       END DO
       END DO
       END DO
    END IF

!$OMP SECTION
!$OMP CRITICAL
    CALL mpi_wait(recvreq(tag_nw), MPI_STATUS_IGNORE, ierr)
!$OMP END CRITICAL
    IF (rank_se /= MPI_PROC_NULL) THEN
       DO k=1, ksize
       DO j=1, slv_s
       DO i=1, slv_e
          a(slv_w+isize+i,j,slv_l+k) = recvbuf_r8_se(slv_e*slv_s*(k-1)+slv_e*(j-1)+i)
       END DO
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, ksize
       DO j=1, slv_s
       DO i=1, slv_e
          a(slv_w+isize+i,j,slv_l+k) = fill
       END DO
       END DO
       END DO
    END IF

!$OMP SECTION
!$OMP CRITICAL
    CALL mpi_wait(recvreq(tag_ne), MPI_STATUS_IGNORE, ierr)
!$OMP END CRITICAL
    IF (rank_sw /= MPI_PROC_NULL) THEN
       DO k=1, ksize
       DO j=1, slv_s
       DO i=1, slv_w
          a(i, j,slv_l+k) = recvbuf_r8_sw(slv_w*slv_s*(k-1)+slv_w*(j-1)+i)
       END DO
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, ksize
       DO j=1, slv_s
       DO i=1, slv_w
          a(i, j,slv_l+k) = fill
       END DO
       END DO
       END DO
    END IF

!$OMP SECTION
!$OMP CRITICAL
    CALL mpi_wait(recvreq(tag_w), MPI_STATUS_IGNORE, ierr)
!$OMP END CRITICAL
    IF (rank_e /= MPI_PROC_NULL) THEN
       DO k=1, ksize
       DO j=1, jsize
       DO i=1, slv_e
          a(slv_w+isize+i,slv_s+j,slv_l+k) = recvbuf_r8_e(slv_e*jsize*(k-1)+slv_e*(j-1)+i)
       END DO
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, ksize
       DO j=1, jsize
       DO i=1, slv_e
          a(slv_w+isize+i,slv_s+j,slv_l+k) = fill
       END DO
       END DO
       END DO
    END IF

!$OMP SECTION
!$OMP CRITICAL
    CALL mpi_wait(recvreq(tag_e), MPI_STATUS_IGNORE, ierr)
!$OMP END CRITICAL
    IF (rank_w /= MPI_PROC_NULL) THEN
       DO k=1, ksize
       DO j=1, jsize
       DO i=1, slv_w
          a(i,slv_s+j,slv_l+k) = recvbuf_r8_w(slv_w*jsize*(k-1)+slv_w*(j-1)+i)
       END DO
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, ksize
       DO j=1, jsize
       DO i=1, slv_w
          a(i,slv_s+j,slv_l+k) = fill
       END DO
       END DO
       END DO
    END IF

!$OMP SECTION
!$OMP CRITICAL
    CALL mpi_wait(recvreq(tag_s), MPI_STATUS_IGNORE, ierr)
!$OMP END CRITICAL
    IF (rank_n /= MPI_PROC_NULL) THEN
       DO k=1, ksize
       DO j=1, slv_n
       DO i=1, isize
          a(slv_w+i,slv_s+jsize+j,slv_l+k) = recvbuf_r8_n(isize*slv_n*(k-1)+isize*(j-1)+i)
       END DO
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, ksize
       DO j=1, slv_n
       DO i=1, isize
          a(slv_w+i,slv_s+jsize+j,slv_l+k) = fill
       END DO
       END DO
       END DO
    END IF

!$OMP SECTION
!$OMP CRITICAL
    CALL mpi_wait(recvreq(tag_n), MPI_STATUS_IGNORE, ierr)
!$OMP END CRITICAL
    IF (rank_s /= MPI_PROC_NULL) THEN
       DO k=1, ksize
       DO j=1, slv_s
       DO i=1, isize
          a(slv_w+i,j,slv_l+k) = recvbuf_r8_s(isize*slv_s*(k-1)+isize*(j-1)+i)
       END DO
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, ksize
       DO j=1, slv_s
       DO i=1, isize
          a(slv_w+i,j,slv_l+k) = fill
       END DO
       END DO
       END DO
    END IF

#ifdef PARALLEL3D
!$OMP SECTION
!$OMP CRITICAL
    CALL mpi_wait(recvreq(tag_lsw), MPI_STATUS_IGNORE, ierr)
!$OMP END CRITICAL
    IF (rank_une /= MPI_PROC_NULL) THEN
       DO k=1, slv_u
       DO j=1, slv_n
       DO i=1, slv_e
          a(slv_w+isize+i,slv_s+jsize+j,slv_l+ksize+k) = recvbuf_r8_une(slv_e*slv_n*(k-1)+slv_e*(j-1)+i)
       END DO
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, slv_u
       DO j=1, slv_n
       DO i=1, slv_e
          a(slv_w+isize+i,slv_s+jsize+j,slv_l+ksize+k) = fill
       END DO
       END DO
       END DO
    END IF

!$OMP SECTION
!$OMP CRITICAL
    CALL mpi_wait(recvreq(tag_lse), MPI_STATUS_IGNORE, ierr)
!$OMP END CRITICAL
    IF (rank_unw /= MPI_PROC_NULL) THEN
       DO k=1, slv_u
       DO j=1, slv_n
       DO i=1, slv_w
          a(i,slv_s+jsize+j,slv_l+ksize+k) = recvbuf_r8_unw(slv_w*slv_n*(k-1)+slv_w*(j-1)+i)
       END DO
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, slv_u
       DO j=1, slv_n
       DO i=1, slv_w
          a(i,slv_s+jsize+j,slv_l+ksize+k) = fill
       END DO
       END DO
       END DO
    END IF

!$OMP SECTION
!$OMP CRITICAL
    CALL mpi_wait(recvreq(tag_lnw), MPI_STATUS_IGNORE, ierr)
!$OMP END CRITICAL
    IF (rank_use /= MPI_PROC_NULL) THEN
       DO k=1, slv_u
       DO j=1, slv_s
       DO i=1, slv_e
          a(slv_w+isize+i,j,slv_l+ksize+k) = recvbuf_r8_use(slv_e*slv_s*(k-1)+slv_e*(j-1)+i)
       END DO
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, slv_u
       DO j=1, slv_s
       DO i=1, slv_e
          a(slv_w+isize+i,j,slv_l+ksize+k) = fill
       END DO
       END DO
       END DO
    END IF

!$OMP SECTION
!$OMP CRITICAL
    CALL mpi_wait(recvreq(tag_lne), MPI_STATUS_IGNORE, ierr)
!$OMP END CRITICAL
    IF (rank_usw /= MPI_PROC_NULL) THEN
       DO k=1, slv_u
       DO j=1, slv_s
       DO i=1, slv_w
          a(i, j,slv_l+ksize+k) = recvbuf_r8_usw(slv_w*slv_s*(k-1)+slv_w*(j-1)+i)
       END DO
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, slv_u
       DO j=1, slv_s
       DO i=1, slv_w
          a(i, j,slv_l+ksize+k) = fill
       END DO
       END DO
       END DO
    END IF

!$OMP SECTION
!$OMP CRITICAL
    CALL mpi_wait(recvreq(tag_usw), MPI_STATUS_IGNORE, ierr)
!$OMP END CRITICAL
    IF (rank_lne /= MPI_PROC_NULL) THEN
       DO k=1, slv_l
       DO j=1, slv_n
       DO i=1, slv_e
          a(slv_w+isize+i,slv_s+jsize+j,k) = recvbuf_r8_lne(slv_e*slv_n*(k-1)+slv_e*(j-1)+i)
       END DO
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, slv_l
       DO j=1, slv_n
       DO i=1, slv_e
          a(slv_w+isize+i,slv_s+jsize+j,k) = fill
       END DO
       END DO
       END DO
    END IF

!$OMP SECTION
!$OMP CRITICAL
    CALL mpi_wait(recvreq(tag_use), MPI_STATUS_IGNORE, ierr)
!$OMP END CRITICAL
    IF (rank_lnw /= MPI_PROC_NULL) THEN
       DO k=1, slv_l
       DO j=1, slv_n
       DO i=1, slv_w
          a(i,slv_s+jsize+j,k) = recvbuf_r8_lnw(slv_w*slv_n*(k-1)+slv_w*(j-1)+i)
       END DO
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, slv_l
       DO j=1, slv_n
       DO i=1, slv_w
          a(i,slv_s+jsize+j,k) = fill
       END DO
       END DO
       END DO
    END IF

!$OMP SECTION
!$OMP CRITICAL
    CALL mpi_wait(recvreq(tag_unw), MPI_STATUS_IGNORE, ierr)
!$OMP END CRITICAL
    IF (rank_lse /= MPI_PROC_NULL) THEN
       DO k=1, slv_l
       DO j=1, slv_s
       DO i=1, slv_e
          a(slv_w+isize+i,j,k) = recvbuf_r8_lse(slv_e*slv_s*(k-1)+slv_e*(j-1)+i)
       END DO
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, slv_l
       DO j=1, slv_s
       DO i=1, slv_e
          a(slv_w+isize+i,j,k) = fill
       END DO
       END DO
       END DO
    END IF

!$OMP SECTION
!$OMP CRITICAL
    CALL mpi_wait(recvreq(tag_une), MPI_STATUS_IGNORE, ierr)
!$OMP END CRITICAL
    IF (rank_lsw /= MPI_PROC_NULL) THEN
       DO k=1, slv_l
       DO j=1, slv_s
       DO i=1, slv_w
          a(i,j,k) = recvbuf_r8_lsw(slv_w*slv_s*(k-1)+slv_w*(j-1)+i)
       END DO
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, slv_l
       DO j=1, slv_s
       DO i=1, slv_w
          a(i,j,k) = fill
       END DO
       END DO
       END DO
    END IF

!$OMP SECTION
!$OMP CRITICAL
    CALL mpi_wait(recvreq(tag_lw), MPI_STATUS_IGNORE, ierr)
!$OMP END CRITICAL
    IF (rank_ue /= MPI_PROC_NULL) THEN
       DO k=1, slv_u
       DO j=1, jsize
       DO i=1, slv_e
          a(slv_w+isize+i,slv_s+j,slv_l+ksize+k) = recvbuf_r8_ue(slv_e*jsize*(k-1)+slv_e*(j-1)+i)
       END DO
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, slv_u
       DO j=1, jsize
       DO i=1, slv_e
          a(slv_w+isize+i,slv_s+j,slv_l+ksize+k) = fill
       END DO
       END DO
       END DO
    END IF

!$OMP SECTION
!$OMP CRITICAL
    CALL mpi_wait(recvreq(tag_le), MPI_STATUS_IGNORE, ierr)
!$OMP END CRITICAL
    IF (rank_uw /= MPI_PROC_NULL) THEN
       DO k=1, slv_u
       DO j=1, jsize
       DO i=1, slv_w
          a(i,slv_s+j,slv_l+ksize+k) = recvbuf_r8_uw(slv_w*jsize*(k-1)+slv_w*(j-1)+i)
       END DO
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, slv_u
       DO j=1, jsize
       DO i=1, slv_w
          a(i,slv_s+j,slv_l+ksize+k) = fill
       END DO
       END DO
       END DO
    END IF

!$OMP SECTION
!$OMP CRITICAL
    CALL mpi_wait(recvreq(tag_ls), MPI_STATUS_IGNORE, ierr)
!$OMP END CRITICAL
    IF (rank_un /= MPI_PROC_NULL) THEN
       DO k=1, slv_u
       DO j=1, slv_n
       DO i=1, isize
          a(slv_w+i,slv_s+jsize+j,slv_l+ksize+k) = recvbuf_r8_un(isize*slv_n*(k-1)+isize*(j-1)+i)
       END DO
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, slv_u
       DO j=1, slv_n
       DO i=1, isize
          a(slv_w+i,slv_s+jsize+j,slv_l+ksize+k) = fill
       END DO
       END DO
       END DO
    END IF

!$OMP SECTION
!$OMP CRITICAL
    CALL mpi_wait(recvreq(tag_ln), MPI_STATUS_IGNORE, ierr)
!$OMP END CRITICAL
    IF (rank_us /= MPI_PROC_NULL) THEN
       DO k=1, slv_u
       DO j=1, slv_s
       DO i=1, isize
          a(slv_w+i,j,slv_l+ksize+k) = recvbuf_r8_us(isize*slv_s*(k-1)+isize*(j-1)+i)
       END DO
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, slv_u
       DO j=1, slv_s
       DO i=1, isize
          a(slv_w+i,j,slv_l+ksize+k) = fill
       END DO
       END DO
       END DO
    END IF

!$OMP SECTION
!$OMP CRITICAL
    CALL mpi_wait(recvreq(tag_uw), MPI_STATUS_IGNORE, ierr)
!$OMP END CRITICAL
    IF (rank_le /= MPI_PROC_NULL) THEN
       DO k=1, slv_l
       DO j=1, jsize
       DO i=1, slv_e
          a(slv_w+isize+i,slv_s+j,k) = recvbuf_r8_le(slv_e*jsize*(k-1)+slv_e*(j-1)+i)
       END DO
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, slv_l
       DO j=1, jsize
       DO i=1, slv_e
          a(slv_w+isize+i,slv_s+j,k) = fill
       END DO
       END DO
       END DO
    END IF

!$OMP SECTION
!$OMP CRITICAL
    CALL mpi_wait(recvreq(tag_ue), MPI_STATUS_IGNORE, ierr)
!$OMP END CRITICAL
    IF (rank_lw /= MPI_PROC_NULL) THEN
       DO k=1, slv_l
       DO j=1, jsize
       DO i=1, slv_w
          a(i,slv_s+j,k) = recvbuf_r8_lw(slv_w*jsize*(k-1)+slv_w*(j-1)+i)
       END DO
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, slv_l
       DO j=1, jsize
       DO i=1, slv_w
          a(i,slv_s+j,k) = fill
       END DO
       END DO
       END DO
    END IF

!$OMP SECTION
!$OMP CRITICAL
    CALL mpi_wait(recvreq(tag_us), MPI_STATUS_IGNORE, ierr)
!$OMP END CRITICAL
    IF (rank_ln /= MPI_PROC_NULL) THEN
       DO k=1, slv_l
       DO j=1, slv_n
       DO i=1, isize
          a(slv_w+i,slv_s+jsize+j,k) = recvbuf_r8_ln(isize*slv_n*(k-1)+isize*(j-1)+i)
       END DO
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, slv_l
       DO j=1, slv_n
       DO i=1, isize
          a(slv_w+i,slv_s+jsize+j,k) = fill
       END DO
       END DO
       END DO
    END IF

!$OMP SECTION
!$OMP CRITICAL
    CALL mpi_wait(recvreq(tag_un), MPI_STATUS_IGNORE, ierr)
!$OMP END CRITICAL
    IF (rank_ls /= MPI_PROC_NULL) THEN
       DO k=1, slv_l
       DO j=1, slv_s
       DO i=1, isize
          a(slv_w+i,j,k) = recvbuf_r8_ls(isize*slv_s*(k-1)+isize*(j-1)+i)
       END DO
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, slv_l
       DO j=1, slv_s
       DO i=1, isize
          a(slv_w+i,j,k) = fill
       END DO
       END DO
       END DO
    END IF

!$OMP SECTION
!$OMP CRITICAL
    CALL mpi_wait(recvreq(tag_l), MPI_STATUS_IGNORE, ierr)
!$OMP END CRITICAL
    IF (rank_u /= MPI_PROC_NULL) THEN
       DO k=1, slv_u
       DO j=1, jsize
       DO i=1, isize
          a(slv_w+i,slv_s+j,slv_l+ksize+k) = recvbuf_r8_u(isize*jsize*(k-1)+isize*(j-1)+i)
       END DO
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, slv_u
       DO j=1, jsize
       DO i=1, isize
          a(slv_w+i,slv_s+j,slv_l+ksize+k) = fill
       END DO
       END DO
       END DO
    END IF

!$OMP SECTION
!$OMP CRITICAL
    CALL mpi_wait(recvreq(tag_u), MPI_STATUS_IGNORE, ierr)
!$OMP END CRITICAL
    IF (rank_l /= MPI_PROC_NULL) THEN
       DO k=1, slv_l
       DO j=1, jsize
       DO i=1, isize
          a(slv_w+i,slv_s+j,k) = recvbuf_r8_l(isize*jsize*(k-1)+isize*(j-1)+i)
       END DO
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, slv_l
       DO j=1, jsize
       DO i=1, isize
          a(slv_w+i,slv_s+j,k) = fill
       END DO
       END DO
       END DO
    END IF
!$OMP END SECTIONS
!$OMP END PARALLEL

#else
!$OMP END SECTIONS
!$OMP END PARALLEL

    IF (present(fill)) THEN
       DO k=1, slv_l
          a(:,:,k) = fill
       END DO
       DO k=1, slv_u
          a(:,:,slv_l+ksize+k) = fill
       END DO
    END IF

#endif
    CALL mpi_waitall(n_sendrecv, sendreq, MPI_STATUSES_IGNORE, ierr)

    IF (rank_tp(0,0)/=MPI_PROC_NULL) THEN
       DO k=1, ksize+slv_l+slv_u
       DO j=1, slv_n
       DO i=1, isize+slv_w+slv_e
          sendbuf_r8_n((isize+slv_w+slv_e)*slv_n*(k-1)+(isize+slv_w+slv_e)*(j-1)+i) = a(isize+slv_w+slv_e+1-i,jsize+slv_n+1-j,k)
       END DO
       END DO
       END DO

       CALL mpi_sendrecv(sendbuf_r8_n, (isize+slv_w+slv_e)*slv_n*(ksize+slv_l+slv_u), MPI_REAL8, rank_tp(0,0), 30, &
                         recvbuf_r8_n, (isize+slv_w+slv_e)*slv_n*(ksize+slv_l+slv_u), MPI_REAL8, rank_tp(0,0), 30, comm, MPI_STATUS_IGNORE, ierr)

       DO k=1, ksize+slv_l+slv_u
       DO j=1, slv_n
       DO i=1, isize+slv_w+slv_e
          a(i,jsize+slv_s+j,k) = recvbuf_r8_n((isize+slv_w+slv_e)*slv_n*(k-1)+(isize+slv_w+slv_e)*(j-1)+i)
       END DO
       END DO
       END DO
    END IF

  END SUBROUTINE parallel_nb
#endif

  END SUBROUTINE update_boundary_3d_r8

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE update_boundary_xz_r4(a, fill)
    REAL(4), INTENT(INOUT) :: a(:,:)
    REAL(4), INTENT(IN), OPTIONAL :: fill

    INTEGER :: nx, nz
    INTEGER :: slv_e, slv_w, slv_u, slv_l

    nx = size(a,1) - isize
    nz = size(a,2) - ksize

    CALL assert(nx >= 0 .AND. nx <= 5, "unsported dimension in UPDATE_BOUNDARY_XZ")
    CALL assert(nz >= 0 .AND. nz <= 5, "unsported dimension in UPDATE_BOUNDARY_XZ")

    slv_e = nx/2
    slv_u = nz/2
    slv_w = (nx+1)/2
    slv_l = (nz+1)/2

#ifdef PARALLEL_MPI
    CALL parallel
#else
    CALL single
#endif

  CONTAINS

  SUBROUTINE single
    INTEGER :: i, k

    IF (cycle_x) THEN
       DO i=1, slv_w
          a(i,:) = a(dimx+i,:)
       END DO
       DO i=1, slv_e
          a(slv_w+dimx+i,:) = a(slv_w+i,:)
       END DO
    ELSE IF (present(fill)) THEN
       DO i=1, slv_w
          a(i,:) = fill
       END DO
       DO i=1, slv_e
          a(slv_w+dimx+i,:) = fill
       END DO
    END IF

    IF (cycle_z) THEN
       DO k=1, slv_l
          a(:,k) = a(:,dimz+k)
       END DO
       DO k=1, slv_u
          a(:,slv_l+dimz+k) = a(:,slv_l+k)
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, slv_l
          a(:,k) = fill
       END DO
       DO k=1, slv_u
          a(:,slv_l+dimz+k) = fill
       END DO
    END IF

  END SUBROUTINE single

#ifdef PARALLEL_MPI
  SUBROUTINE parallel
    INTEGER :: req(4), ierr
    INTEGER :: i, k

    CALL mpi_irecv(recvbuf_r4_e, slv_e*ksize, MPI_REAL4, rank_e, 1, comm, req(1), ierr)
    CALL mpi_irecv(recvbuf_r4_w, slv_w*ksize, MPI_REAL4, rank_w, 2, comm, req(2), ierr)

    DO k=1, ksize
       DO i=1, slv_e
          sendbuf_r4_w(slv_e*(k-1)+i) = a(slv_w+i,slv_l+k)
       END DO

       DO i=1, slv_w
          sendbuf_r4_e(slv_w*(k-1)+i) = a(slv_w+isize-slv_w+i,slv_l+k)
       END DO
    END DO

    CALL mpi_isend(sendbuf_r4_w, slv_e*ksize, MPI_REAL4, rank_w, 1, comm, req(3), ierr)
    CALL mpi_isend(sendbuf_r4_e, slv_w*ksize, MPI_REAL4, rank_e, 2, comm, req(4), ierr)

    CALL mpi_waitall(4, req, MPI_STATUSES_IGNORE, ierr)

    IF (rank_w /= MPI_PROC_NULL) THEN
       DO k=1, ksize
       DO i=1, slv_w
          a(i,slv_l+k) = recvbuf_r4_w(slv_w*(k-1)+i)
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, ksize
       DO i=1, slv_w
          a(i,slv_l+k) = fill
       END DO
       END DO
    END IF

    IF (rank_e /= MPI_PROC_NULL) THEN
       DO k=1, ksize
       DO i=1, slv_e
          a(slv_w+isize+i,slv_l+k) = recvbuf_r4_e(slv_e*(k-1)+i)
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, ksize
       DO i=1, slv_e
          a(slv_w+isize+i,slv_l+k) = fill
       END DO
       END DO
    END IF

    CALL mpi_irecv(recvbuf_r4_u, (isize+slv_w+slv_e)*slv_u, MPI_REAL4, rank_u, 3, comm, req(1), ierr)
    CALL mpi_irecv(recvbuf_r4_l, (isize+slv_w+slv_e)*slv_l, MPI_REAL4, rank_l, 4, comm, req(2), ierr)

    DO k=1, slv_u
    DO i=1, isize+slv_w+slv_e
       sendbuf_r4_l((isize+slv_w+slv_e)*(k-1)+i) = a(i,slv_l+k)
    END DO
    END DO

    DO k=1, slv_l
    DO i=1, isize+slv_w+slv_e
       sendbuf_r4_u((isize+slv_w+slv_e)*(k-1)+i) = a(i,ksize+k)
    END DO
    END DO

    CALL mpi_isend(sendbuf_r4_l, (isize+slv_w+slv_e)*slv_u, MPI_REAL4, rank_l, 3, comm, req(3), ierr)
    CALL mpi_isend(sendbuf_r4_u, (isize+slv_w+slv_e)*slv_l, MPI_REAL4, rank_u, 4, comm, req(4), ierr)

    CALL mpi_waitall( 4, req, MPI_STATUSES_IGNORE, ierr)

    IF (rank_l /= MPI_PROC_NULL) THEN
       DO k=1, slv_l
       DO i=1, isize+slv_w+slv_e
          a(i,k) = recvbuf_r4_l((isize+slv_w+slv_e)*(k-1)+i)
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, slv_l
          a(:,k) = fill
       END DO
    END IF

    IF (rank_u /= MPI_PROC_NULL) THEN
       DO k=1, slv_u
       DO i=1, isize+slv_w+slv_e
          a(i,slv_l+ksize+k) = recvbuf_r4_u((isize+slv_w+slv_e)*(k-1)+i)
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, slv_u
          a(:,slv_l+ksize+k) = fill
       END DO
    END IF
  END SUBROUTINE parallel
#endif
  END SUBROUTINE update_boundary_xz_r4

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE update_boundary_xz_r8(a, fill)
    REAL(8), INTENT(INOUT) :: a(:,:)
    REAL(8), INTENT(IN), OPTIONAL :: fill

    INTEGER :: nx, nz
    INTEGER :: slv_e, slv_w, slv_u, slv_l

    nx = size(a,1) - isize
    nz = size(a,2) - ksize

    CALL assert(nx >= 0 .AND. nx <= 5, "unsported dimension in UPDATE_BOUNDARY_XZ")
    CALL assert(nz >= 0 .AND. nz <= 5, "unsported dimension in UPDATE_BOUNDARY_XZ")

    slv_e = nx/2
    slv_u = nz/2
    slv_w = (nx+1)/2
    slv_l = (nz+1)/2

#ifdef PARALLEL_MPI
    CALL parallel
#else
    CALL single
#endif

  CONTAINS

  SUBROUTINE single
    INTEGER :: i, k

    IF (cycle_x) THEN
       DO i=1, slv_w
          a(i,:) = a(dimx+i,:)
       END DO
       DO i=1, slv_e
          a(slv_w+dimx+i,:) = a(slv_w+i,:)
       END DO
    ELSE IF (present(fill)) THEN
       DO i=1, slv_w
          a(i,:) = fill
       END DO
       DO i=1, slv_e
          a(slv_w+dimx+i,:) = fill
       END DO
    END IF

    IF (cycle_z) THEN
       DO k=1, slv_l
          a(:,k) = a(:,dimz+k)
       END DO
       DO k=1, slv_u
          a(:,slv_l+dimz+k) = a(:,slv_l+k)
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, slv_l
          a(:,k) = fill
       END DO
       DO k=1, slv_u
          a(:,slv_l+dimz+k) = fill
       END DO
    END IF
  END SUBROUTINE single

#ifdef PARALLEL_MPI
  SUBROUTINE parallel
    INTEGER :: req(4), ierr
    INTEGER :: i, k

    CALL mpi_irecv(recvbuf_r8_e, slv_e*ksize, MPI_REAL8, rank_e, 1, comm, req(1), ierr)
    CALL mpi_irecv(recvbuf_r8_w, slv_w*ksize, MPI_REAL8, rank_w, 2, comm, req(2), ierr)

    DO k=1, ksize
       DO i=1, slv_e
          sendbuf_r8_w(slv_e*(k-1)+i) = a(slv_w+i,slv_l+k)
       END DO

       DO i=1, slv_w
          sendbuf_r8_e(slv_w*(k-1)+i) = a(slv_w+isize-slv_w+i,slv_l+k)
       END DO
    END DO

    CALL mpi_isend(sendbuf_r8_w, slv_e*ksize, MPI_REAL8, rank_w, 1, comm, req(3), ierr)
    CALL mpi_isend(sendbuf_r8_e, slv_w*ksize, MPI_REAL8, rank_e, 2, comm, req(4), ierr)

    CALL mpi_waitall(4, req, MPI_STATUSES_IGNORE, ierr)

    IF (rank_w /= MPI_PROC_NULL) THEN
       DO k=1, ksize
       DO i=1, slv_w
          a(i,slv_l+k) = recvbuf_r8_w(slv_w*(k-1)+i)
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, ksize
       DO i=1, slv_w
          a(i,slv_l+k) = fill
       END DO
       END DO
    END IF

    IF (rank_e /= MPI_PROC_NULL) THEN
       DO k=1, ksize
       DO i=1, slv_e
          a(slv_w+isize+i,slv_l+k) = recvbuf_r8_e(slv_e*(k-1)+i)
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, ksize
       DO i=1, slv_e
          a(slv_w+isize+i,slv_l+k) = fill
       END DO
       END DO
    END IF

    CALL mpi_irecv(recvbuf_r8_u, (isize+slv_w+slv_e)*slv_u, MPI_REAL8, rank_u, 3, comm, req(1), ierr)
    CALL mpi_irecv(recvbuf_r8_l, (isize+slv_w+slv_e)*slv_l, MPI_REAL8, rank_l, 4, comm, req(2), ierr)

    DO k=1, slv_u
    DO i=1, isize+slv_w+slv_e
       sendbuf_r8_l((isize+slv_w+slv_e)*(k-1)+i) = a(i,slv_l+k)
    END DO
    END DO

    DO k=1, slv_l
    DO i=1, isize+slv_w+slv_e
       sendbuf_r8_u((isize+slv_w+slv_e)*(k-1)+i) = a(i,ksize+k)
    END DO
    END DO

    CALL mpi_isend(sendbuf_r8_l, (isize+slv_w+slv_e)*slv_u, MPI_REAL8, rank_l, 3, comm, req(3), ierr)
    CALL mpi_isend(sendbuf_r8_u, (isize+slv_w+slv_e)*slv_l, MPI_REAL8, rank_u, 4, comm, req(4), ierr)

    CALL mpi_waitall( 4, req, MPI_STATUSES_IGNORE, ierr)

    IF (rank_l /= MPI_PROC_NULL) THEN
       DO k=1, slv_l
       DO i=1, isize+slv_w+slv_e
          a(i,k) = recvbuf_r8_l((isize+slv_w+slv_e)*(k-1)+i)
      END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, slv_l
          a(:,k) = fill
       END DO
    END IF

    IF (rank_u /= MPI_PROC_NULL) THEN
       DO k=1, slv_u
       DO i=1, isize+slv_w+slv_e
          a(i,slv_l+ksize+k) = recvbuf_r8_u((isize+slv_w+slv_e)*(k-1)+i)
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, slv_u
          a(:,slv_l+ksize+k) = fill
       END DO
    END IF
  END SUBROUTINE parallel
#endif
  END SUBROUTINE update_boundary_xz_r8

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE update_boundary_yz_r4(a, fill)
    REAL(4), INTENT(INOUT) :: a(:,:)
    REAL(4), INTENT(IN), OPTIONAL :: fill

    INTEGER :: ny, nz
    INTEGER :: slv_n, slv_s, slv_u, slv_l

    ny = size(a,1) - jsize
    nz = size(a,2) - ksize

    CALL assert(ny >= 0 .AND. ny <= 5, "unsported dimension in UPDATE_BOUNDARY_YZ")
    CALL assert(nz >= 0 .AND. nz <= 5, "unsported dimension in UPDATE_BOUNDARY_YZ")

    slv_n = ny/2
    slv_u = nz/2
    slv_s = (ny+1)/2
    slv_l = (nz+1)/2

#ifdef PARALLEL_MPI
    CALL parallel
#else
    CALL single
#endif

  CONTAINS

  SUBROUTINE single
    INTEGER :: j, k

    IF (cycle_y) THEN
       DO j=1, slv_s
          a(j,:) = a(dimy+j,:)
       END DO
       DO j=1, slv_n
          a(slv_s+dimy+j,:) = a(slv_s+j,:)
       END DO
    ELSE IF (present(fill)) THEN
       DO j=1, slv_s
          a(j,:) = fill
       END DO
       DO j=1, slv_n
          a(slv_s+dimy+j,:) = fill
       END DO
    END IF

    IF (cycle_z) THEN
       DO k=1, slv_l
          a(:,k) = a(:,dimz+k)
       END DO
       DO k=1, slv_u
          a(:,slv_l+dimz+k) = a(:,slv_l+k)
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, slv_l
          a(:,k) = fill
       END DO
       DO k=1, slv_u
          a(:,slv_l+dimz+k) = fill
       END DO
    END IF
  END SUBROUTINE single

#ifdef PARALLEL_MPI
  SUBROUTINE parallel
    INTEGER :: req(4), ierr
    INTEGER :: j, k

    CALL mpi_irecv(recvbuf_r4_n, slv_n*ksize, MPI_REAL4, rank_n, 1, comm, req(1), ierr)
    CALL mpi_irecv(recvbuf_r4_s, slv_s*ksize, MPI_REAL4, rank_s, 2, comm, req(2), ierr)

    DO k=1, ksize
       DO j=1, slv_n
          sendbuf_r4_s(slv_n*(k-1)+j) = a(slv_s+j,slv_l+k)
       END DO

       DO j=1, slv_s
          sendbuf_r4_n(slv_s*(k-1)+j) = a(slv_s+jsize-slv_s+j,slv_l+k)
       END DO
    END DO

    CALL mpi_isend(sendbuf_r4_s, slv_n*ksize, MPI_REAL4, rank_s, 1, comm, req(3), ierr)
    CALL mpi_isend(sendbuf_r4_n, slv_s*ksize, MPI_REAL4, rank_n, 2, comm, req(4), ierr)

    CALL mpi_waitall(4, req, MPI_STATUSES_IGNORE, ierr)

    IF (rank_s /= MPI_PROC_NULL) THEN
       DO k=1, ksize
       DO j=1, slv_s
          a(j,slv_l+k) = recvbuf_r4_s(slv_s*(k-1)+j)
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, ksize
       DO j=1, slv_s
          a(j,slv_l+k) = fill
       END DO
       END DO
    END IF

    IF (rank_n /= MPI_PROC_NULL) THEN
       DO k=1, ksize
       DO j=1, slv_n
          a(slv_s+jsize+j,slv_l+k) = recvbuf_r4_n(slv_n*(k-1)+j)
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, ksize
       DO j=1, slv_n
          a(slv_s+jsize+j,slv_l+k) = fill
       END DO
       END DO
    END IF

    CALL mpi_irecv(recvbuf_r4_u, (jsize+slv_s+slv_n)*slv_u, MPI_REAL4, rank_u, 3, comm, req(1), ierr)
    CALL mpi_irecv(recvbuf_r4_l, (jsize+slv_s+slv_n)*slv_l, MPI_REAL4, rank_l, 4, comm, req(2), ierr)

    DO k=1, slv_u
    DO j=1, jsize+slv_s+slv_n
       sendbuf_r4_l((jsize+slv_s+slv_n)*(k-1)+j) = a(j,slv_l+k)
    END DO
    END DO

    DO k=1, slv_l
    DO j=1, jsize+slv_s+slv_n
       sendbuf_r4_u((jsize+slv_s+slv_n)*(k-1)+j) = a(j,ksize+k)
    END DO
    END DO

    CALL mpi_isend(sendbuf_r4_l, (jsize+slv_s+slv_n)*slv_u, MPI_REAL4, rank_l, 3, comm, req(3), ierr)
    CALL mpi_isend(sendbuf_r4_u, (jsize+slv_s+slv_n)*slv_l, MPI_REAL4, rank_u, 4, comm, req(4), ierr)

    CALL mpi_waitall( 4, req, MPI_STATUSES_IGNORE, ierr)

    IF (rank_l /= MPI_PROC_NULL) THEN
       DO k=1, slv_l
       DO j=1, jsize+slv_s+slv_n
          a(j,k) = recvbuf_r4_l((jsize+slv_s+slv_n)*(k-1)+j)
      END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, slv_l
          a(:,k) = fill
       END DO
    END IF

    IF (rank_u /= MPI_PROC_NULL) THEN
       DO k=1, slv_u
       DO j=1, jsize+slv_s+slv_n
          a(j,slv_l+ksize+k) = recvbuf_r4_u((jsize+slv_s+slv_n)*(k-1)+j)
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, slv_u
          a(:,slv_l+ksize+k) = fill
       END DO
    END IF
  END SUBROUTINE parallel
#endif
  END SUBROUTINE update_boundary_yz_r4

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE update_boundary_yz_r8(a, fill)
    REAL(8), INTENT(INOUT) :: a(:,:)
    REAL(8), INTENT(IN), OPTIONAL :: fill

    INTEGER :: ny, nz
    INTEGER :: slv_n, slv_s, slv_u, slv_l

    ny = size(a,1) - jsize
    nz = size(a,2) - ksize

    CALL assert(ny >= 0 .AND. ny <= 5, "unsported dimension in UPDATE_BOUNDARY_YZ")
    CALL assert(nz >= 0 .AND. nz <= 5, "unsported dimension in UPDATE_BOUNDARY_YZ")

    slv_n = ny/2
    slv_u = nz/2
    slv_s = (ny+1)/2
    slv_l = (nz+1)/2

#ifdef PARALLEL_MPI
    CALL parallel
#else
    CALL single
#endif

  CONTAINS

  SUBROUTINE single
    INTEGER :: j, k

    IF (cycle_y) THEN
       DO j=1, slv_s
          a(j,:) = a(dimy+j,:)
       END DO
       DO j=1, slv_n
          a(slv_s+dimy+j,:) = a(slv_s+j,:)
       END DO
    ELSE IF (present(fill)) THEN
       DO j=1, slv_s
          a(j,:) = fill
       END DO
       DO j=1, slv_n
          a(slv_s+dimy+j,:) = fill
       END DO
    END IF

    IF (cycle_z) THEN
       DO k=1, slv_l
          a(:,k) = a(:,dimz+k)
       END DO
       DO k=1, slv_u
          a(:,slv_l+dimz+k) = a(:,slv_l+k)
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, slv_l
          a(:,k) = fill
       END DO
       DO k=1, slv_u
          a(:,slv_l+dimz+k) = fill
       END DO
    END IF
  END SUBROUTINE single

#ifdef PARALLEL_MPI
  SUBROUTINE parallel
    INTEGER :: req(4), ierr
    INTEGER :: j, k

    CALL mpi_irecv(recvbuf_r8_n, slv_n*ksize, MPI_REAL8, rank_n, 1, comm, req(1), ierr)
    CALL mpi_irecv(recvbuf_r8_s, slv_s*ksize, MPI_REAL8, rank_s, 2, comm, req(2), ierr)

    DO k=1, ksize
       DO j=1, slv_n
          sendbuf_r8_s(slv_n*(k-1)+j) = a(slv_s+j,slv_l+k)
       END DO

       DO j=1, slv_s
          sendbuf_r8_n(slv_s*(k-1)+j) = a(slv_s+jsize-slv_s+j,slv_l+k)
       END DO
    END DO

    CALL mpi_isend(sendbuf_r8_s, slv_n*ksize, MPI_REAL8, rank_s, 1, comm, req(3), ierr)
    CALL mpi_isend(sendbuf_r8_n, slv_s*ksize, MPI_REAL8, rank_n, 2, comm, req(4), ierr)

    CALL mpi_waitall(4, req, MPI_STATUSES_IGNORE, ierr)

    IF (rank_s /= MPI_PROC_NULL) THEN
       DO k=1, ksize
       DO j=1, slv_s
          a(j,slv_l+k) = recvbuf_r8_s(slv_s*(k-1)+j)
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, ksize
       DO j=1, slv_s
          a(j,slv_l+k) = fill
       END DO
       END DO
    END IF

    IF (rank_n /= MPI_PROC_NULL) THEN
       DO k=1, ksize
       DO j=1, slv_n
          a(slv_s+jsize+j,slv_l+k) = recvbuf_r8_n(slv_n*(k-1)+j)
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, ksize
       DO j=1, slv_n
          a(slv_s+jsize+j,slv_l+k) = fill
       END DO
       END DO
    END IF

    CALL mpi_irecv(recvbuf_r8_u, (jsize+slv_s+slv_n)*slv_u, MPI_REAL8, rank_u, 3, comm, req(1), ierr)
    CALL mpi_irecv(recvbuf_r8_l, (jsize+slv_s+slv_n)*slv_l, MPI_REAL8, rank_l, 4, comm, req(2), ierr)

    DO k=1, slv_u
    DO j=1, jsize+slv_s+slv_n
       sendbuf_r8_l((jsize+slv_s+slv_n)*(k-1)+j) = a(j,slv_l+k)
    END DO
    END DO

    DO k=1, slv_l
    DO j=1, jsize+slv_s+slv_n
       sendbuf_r8_u((jsize+slv_s+slv_n)*(k-1)+j) = a(j,ksize+k)
    END DO
    END DO

    CALL mpi_isend(sendbuf_r8_l, (jsize+slv_s+slv_n)*slv_u, MPI_REAL8, rank_l, 3, comm, req(3), ierr)
    CALL mpi_isend(sendbuf_r8_u, (jsize+slv_s+slv_n)*slv_l, MPI_REAL8, rank_u, 4, comm, req(4), ierr)

    CALL mpi_waitall( 4, req, MPI_STATUSES_IGNORE, ierr)

    IF (rank_l /= MPI_PROC_NULL) THEN
       DO k=1, slv_l
       DO j=1, jsize+slv_s+slv_n
          a(j,k) = recvbuf_r8_l((jsize+slv_s+slv_n)*(k-1)+j)
      END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, slv_l
          a(:,k) = fill
       END DO
    END IF

    IF (rank_u /= MPI_PROC_NULL) THEN
       DO k=1, slv_u
       DO j=1, jsize+slv_s+slv_n
          a(j,slv_l+ksize+k) = recvbuf_r8_u((jsize+slv_s+slv_n)*(k-1)+j)
       END DO
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, slv_u
          a(:,slv_l+ksize+k) = fill
       END DO
    END IF
  END SUBROUTINE parallel
#endif
  END SUBROUTINE update_boundary_yz_r8

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE update_boundary_x_r4(a, fill)
    REAL(4), INTENT(INOUT) :: a(:)
    REAL(4), INTENT(IN), OPTIONAL :: fill

    INTEGER :: nx
    INTEGER :: slv_e, slv_w

    nx = size(a) - isize

    CALL assert(nx >= 0 .AND. nx <= 5, "unsported dimension in UPDATE_BOUNDARY_X")

    slv_e = nx/2
    slv_w = (nx+1)/2

#ifdef PARALLEL_MPI
    CALL parallel
#else
    CALL single
#endif

  CONTAINS

  SUBROUTINE single
    INTEGER :: i

    IF (cycle_x) THEN
       DO i=1, slv_w
          a(i) = a(dimx+i)
       END DO
       DO i=1, slv_e
          a(slv_w+dimx+i) = a(slv_w+i)
       END DO
    ELSE IF (present(fill)) THEN
       DO i=1, slv_w
          a(i) = fill
       END DO
       DO i=1, slv_e
          a(slv_w+dimx+i) = fill
       END DO
    END IF
  END SUBROUTINE single

#ifdef PARALLEL_MPI
  SUBROUTINE parallel
    INTEGER :: req(4), ierr
    INTEGER :: i

    CALL mpi_irecv(recvbuf_r4_e, slv_e, MPI_REAL4, rank_e, 1, comm, req(1), ierr)
    CALL mpi_irecv(recvbuf_r4_w, slv_w, MPI_REAL4, rank_w, 2, comm, req(2), ierr)

    DO i=1, slv_e
       sendbuf_r4_w(i) = a(slv_w+i)
    END DO

    DO i=1, slv_w
       sendbuf_r4_e(i) = a(isize+i)
    END DO

    CALL mpi_isend(sendbuf_r4_w, slv_e, MPI_REAL4, rank_w, 1, comm, req(3), ierr)
    CALL mpi_isend(sendbuf_r4_e, slv_w, MPI_REAL4, rank_e, 2, comm, req(4), ierr)

    CALL mpi_waitall(4, req, MPI_STATUSES_IGNORE, ierr)

    IF (rank_w /= MPI_PROC_NULL) THEN
       DO i=1, slv_w
          a(i) = recvbuf_r4_w(i)
       END DO
    ELSE IF (present(fill)) THEN
       DO i=1, slv_w
          a(i) = fill
       END DO
    END IF

    IF (rank_e /= MPI_PROC_NULL) THEN
       DO i=1, slv_e
          a(slv_w+isize+i) = recvbuf_r4_e(i)
       END DO
    ELSE IF (present(fill)) THEN
       DO i=1, slv_e
          a(slv_w+isize+i) = fill
       END DO
    END IF

  END SUBROUTINE parallel
#endif
  END SUBROUTINE update_boundary_x_r4

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE update_boundary_x_r8(a, fill)
    REAL(8), INTENT(INOUT) :: a(:)
    REAL(8), INTENT(IN), OPTIONAL :: fill

    INTEGER :: nx
    INTEGER :: slv_e, slv_w

    nx = size(a) - isize

    CALL assert(nx >= 0 .AND. nx <= 5, "unsported dimension in UPDATE_BOUNDARY_X")

    slv_e = nx/2
    slv_w = (nx+1)/2

#ifdef PARALLEL_MPI
    CALL parallel
#else
    CALL single
#endif

  CONTAINS

  SUBROUTINE single
    INTEGER :: i

    IF (cycle_x) THEN
       DO i=1, slv_w
          a(i) = a(dimx+i)
       END DO
       DO i=1, slv_e
          a(slv_w+dimx+i) = a(slv_w+i)
       END DO
    ELSE IF (present(fill)) THEN
       DO i=1, slv_w
          a(i) = fill
       END DO
       DO i=1, slv_e
          a(slv_w+dimx+i) = fill
       END DO
    END IF
  END SUBROUTINE single

#ifdef PARALLEL_MPI
  SUBROUTINE parallel
    INTEGER :: req(4), ierr
    INTEGER :: i

    CALL mpi_irecv(recvbuf_r8_e, slv_e, MPI_REAL8, rank_e, 1, comm, req(1), ierr)
    CALL mpi_irecv(recvbuf_r8_w, slv_w, MPI_REAL8, rank_w, 2, comm, req(2), ierr)

    DO i=1, slv_e
       sendbuf_r8_w(i) = a(slv_w+i)
    END DO

    DO i=1, slv_w
       sendbuf_r8_e(i) = a(isize+i)
    END DO

    CALL mpi_isend(sendbuf_r8_w, slv_e, MPI_REAL8, rank_w, 1, comm, req(3), ierr)
    CALL mpi_isend(sendbuf_r8_e, slv_w, MPI_REAL8, rank_e, 2, comm, req(4), ierr)

    CALL mpi_waitall(4, req, MPI_STATUSES_IGNORE, ierr)

    IF (rank_w /= MPI_PROC_NULL) THEN
       DO i=1, slv_w
          a(i) = recvbuf_r8_w(i)
       END DO
    ELSE IF (present(fill)) THEN
       DO i=1, slv_w
          a(i) = fill
       END DO
    END IF

    IF (rank_e /= MPI_PROC_NULL) THEN
       DO i=1, slv_e
          a(slv_w+isize+i) = recvbuf_r8_e(i)
       END DO
    ELSE IF (present(fill)) THEN
       DO i=1, slv_e
          a(slv_w+isize+i) = fill
       END DO
    END IF

  END SUBROUTINE parallel
#endif
  END SUBROUTINE update_boundary_x_r8

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE update_boundary_y_r4(a, fill)
    REAL(4), INTENT(INOUT) :: a(:)
    REAL(4), INTENT(IN), OPTIONAL :: fill

    INTEGER :: ny
    INTEGER :: slv_n, slv_s

    ny = size(a) - jsize

    CALL assert(ny >= 0 .AND. ny <= 5, "unsported dimension in UPDATE_BOUNDARY_Y")

    slv_n = ny/2
    slv_s = (ny+1)/2

#ifdef PARALLEL_MPI
    CALL parallel
#else
    CALL single
#endif

  CONTAINS

  SUBROUTINE single
    INTEGER :: j

    IF (cycle_y) THEN
       DO j=1, slv_s
          a(j) = a(dimy+j)
       END DO
       DO j=1, slv_n
          a(slv_s+dimy+j) = a(slv_s+j)
       END DO
    ELSE IF (present(fill)) THEN
       DO j=1, slv_s
          a(j) = fill
       END DO
       DO j=1, slv_n
          a(slv_s+dimy+j) = fill
       END DO
    END IF
  END SUBROUTINE single

#ifdef PARALLEL_MPI
  SUBROUTINE parallel
    INTEGER :: req(4), ierr
    INTEGER :: j

    CALL mpi_irecv(recvbuf_r4_n, slv_n, MPI_REAL4, rank_n, 1, comm, req(1), ierr)
    CALL mpi_irecv(recvbuf_r4_s, slv_s, MPI_REAL4, rank_s, 2, comm, req(2), ierr)

    DO j=1, slv_n
       sendbuf_r4_s(j) = a(slv_s+j)
    END DO

    DO j=1, slv_s
       sendbuf_r4_n(j) = a(jsize+j)
    END DO

    CALL mpi_isend(sendbuf_r4_s, slv_n, MPI_REAL4, rank_s, 1, comm, req(3), ierr)
    CALL mpi_isend(sendbuf_r4_n, slv_s, MPI_REAL4, rank_n, 2, comm, req(4), ierr)

    CALL mpi_waitall(4, req, MPI_STATUSES_IGNORE, ierr)

    IF (rank_s /= MPI_PROC_NULL) THEN
       DO j=1, slv_s
          a(j) = recvbuf_r4_s(j)
       END DO
    ELSE IF (present(fill)) THEN
       DO j=1, slv_s
          a(j) = fill
       END DO
    END IF

    IF (rank_n /= MPI_PROC_NULL) THEN
       DO j=1, slv_n
          a(slv_s+jsize+j) = recvbuf_r4_n(j)
       END DO
    ELSE IF (present(fill)) THEN
       DO j=1, slv_n
          a(slv_s+jsize+j) = fill
       END DO
    END IF

  END SUBROUTINE parallel
#endif
  END SUBROUTINE update_boundary_y_r4

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE update_boundary_y_r8(a, fill)
    REAL(8), INTENT(INOUT) :: a(:)
    REAL(8), INTENT(IN), OPTIONAL :: fill

    INTEGER :: ny
    INTEGER :: slv_n, slv_s

    ny = size(a) - jsize

    CALL assert(ny >= 0 .AND. ny <= 5, "unsported dimension in UPDATE_BOUNDARY_Y")

    slv_n = ny/2
    slv_s = (ny+1)/2

#ifdef PARALLEL_MPI
    CALL parallel
#else
    CALL single
#endif

  CONTAINS

  SUBROUTINE single
    INTEGER :: j

    IF (cycle_y) THEN
       DO j=1, slv_s
          a(j) = a(dimy+j)
       END DO
       DO j=1, slv_n
          a(slv_s+dimy+j) = a(slv_s+j)
       END DO
    ELSE IF (present(fill)) THEN
       DO j=1, slv_s
          a(j) = fill
       END DO
       DO j=1, slv_n
          a(slv_s+dimy+j) = fill
       END DO
    END IF
  END SUBROUTINE single

#ifdef PARALLEL_MPI
  SUBROUTINE parallel
    INTEGER :: req(4), ierr
    INTEGER :: j

    CALL mpi_irecv(recvbuf_r8_n, slv_n, MPI_REAL8, rank_n, 1, comm, req(1), ierr)
    CALL mpi_irecv(recvbuf_r8_s, slv_s, MPI_REAL8, rank_s, 2, comm, req(2), ierr)

    DO j=1, slv_n
       sendbuf_r8_s(j) = a(slv_s+j)
    END DO

    DO j=1, slv_s
       sendbuf_r8_n(j) = a(jsize+j)
    END DO

    CALL mpi_isend(sendbuf_r8_s, slv_n, MPI_REAL8, rank_s, 1, comm, req(3), ierr)
    CALL mpi_isend(sendbuf_r8_n, slv_s, MPI_REAL8, rank_n, 2, comm, req(4), ierr)

    CALL mpi_waitall(4, req, MPI_STATUSES_IGNORE, ierr)

    IF (rank_s /= MPI_PROC_NULL) THEN
       DO j=1, slv_s
          a(j) = recvbuf_r8_s(j)
       END DO
    ELSE IF (present(fill)) THEN
       DO j=1, slv_s
          a(j) = fill
       END DO
    END IF

    IF (rank_n /= MPI_PROC_NULL) THEN
       DO j=1, slv_n
          a(slv_s+jsize+j) = recvbuf_r8_n(j)
       END DO
    ELSE IF (present(fill)) THEN
       DO j=1, slv_n
          a(slv_s+jsize+j) = fill
       END DO
    END IF

  END SUBROUTINE parallel
#endif
  END SUBROUTINE update_boundary_y_r8

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE update_boundary_z_r4(a, fill)
    REAL(4), INTENT(INOUT) :: a(:)
    REAL(4), INTENT(IN), OPTIONAL :: fill

    INTEGER :: nz
    INTEGER :: slv_u, slv_l

    nz = size(a) - ksize

    CALL assert(nz >= 0 .AND. nz <= 5, "unsported dimension in UPDATE_BOUNDARY_Z")

    slv_u = nz/2
    slv_l = (nz+1)/2

#ifdef PARALLEL_MPI
    CALL parallel
#else
    CALL single
#endif

  CONTAINS

  SUBROUTINE single
    INTEGER :: k

    IF (cycle_z) THEN
       DO k=1, slv_l
          a(k) = a(dimz+k)
       END DO
       DO k=1, slv_u
          a(slv_l+dimz+k) = a(slv_l+k)
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, slv_l
          a(k) = fill
       END DO
       DO k=1, slv_u
          a(slv_l+dimz+k) = fill
       END DO
    END IF
  END SUBROUTINE single

#ifdef PARALLEL_MPI
  SUBROUTINE parallel
    INTEGER :: req(4), ierr
    INTEGER :: k

    CALL mpi_irecv(recvbuf_r4_u, slv_u, MPI_REAL4, rank_u, 1, comm, req(1), ierr)
    CALL mpi_irecv(recvbuf_r4_l, slv_l, MPI_REAL4, rank_l, 2, comm, req(2), ierr)

    DO k=1, slv_u
       sendbuf_r4_l(k) = a(slv_l+k)
    END DO

    DO k=1, slv_l
       sendbuf_r4_u(k) = a(ksize-slv_l+k)
    END DO

    CALL mpi_isend(sendbuf_r4_l, slv_u, MPI_REAL4, rank_l, 1, comm, req(3), ierr)
    CALL mpi_isend(sendbuf_r4_u, slv_l, MPI_REAL4, rank_u, 2, comm, req(4), ierr)

    CALL mpi_waitall(4, req, MPI_STATUSES_IGNORE, ierr)

    IF (rank_l /= MPI_PROC_NULL) THEN
       DO k=1, slv_l
          a(k) = recvbuf_r4_l(k)
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, slv_l
          a(k) = fill
       END DO
    END IF

    IF (rank_u /= MPI_PROC_NULL) THEN
       DO k=1, slv_u
          a(slv_l+ksize+k) = recvbuf_r4_u(k)
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, slv_u
          a(slv_l+ksize+k) = fill
       END DO
    END IF

  END SUBROUTINE parallel
#endif
  END SUBROUTINE update_boundary_z_r4

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE update_boundary_z_r8(a, fill)
    REAL(8), INTENT(INOUT) :: a(:)
    REAL(8), INTENT(IN), OPTIONAL :: fill

    INTEGER :: nz
    INTEGER :: slv_u, slv_l

    nz = size(a) - ksize

    CALL assert(nz >= 0 .AND. nz <= 5,  "unsported dimension in UPDATE_BOUNDARY_Z")

    slv_u = nz/2
    slv_l = (nz+1)/2

#ifdef PARALLEL_MPI
    CALL parallel
#else
    CALL single
#endif

  CONTAINS

  SUBROUTINE single
    INTEGER :: k

    IF (cycle_z) THEN
       DO k=1, slv_l
          a(k) = a(dimz+k)
       END DO
       DO k=1, slv_u
          a(slv_l+dimz+k) = a(slv_l+k)
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, slv_l
          a(k) = fill
       END DO
       DO k=1, slv_u
          a(slv_l+dimz+k) = fill
       END DO
    END IF
  END SUBROUTINE single

#ifdef PARALLEL_MPI
  SUBROUTINE parallel
    INTEGER :: req(4), ierr
    INTEGER :: k

    CALL mpi_irecv(recvbuf_r8_u, slv_u, MPI_REAL8, rank_u, 1, comm, req(1), ierr)
    CALL mpi_irecv(recvbuf_r8_l, slv_l, MPI_REAL8, rank_l, 2, comm, req(2), ierr)

    DO k=1, slv_u
       sendbuf_r8_l(k) = a(slv_l+k)
    END DO

    DO k=1, slv_l
       sendbuf_r8_u(k) = a(ksize-slv_l+k)
    END DO

    CALL mpi_isend(sendbuf_r8_l, slv_u, MPI_REAL8, rank_l, 1, comm, req(3), ierr)
    CALL mpi_isend(sendbuf_r8_u, slv_l, MPI_REAL8, rank_u, 2, comm, req(4), ierr)

    CALL mpi_waitall(4, req, MPI_STATUSES_IGNORE, ierr)

    IF (rank_l /= MPI_PROC_NULL) THEN
       DO k=1, slv_l
          a(k) = recvbuf_r8_l(k)
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, slv_l
          a(k) = fill
       END DO
    END IF

    IF (rank_u /= MPI_PROC_NULL) THEN
       DO k=1, slv_u
          a(slv_l+ksize+k) = recvbuf_r8_u(k)
       END DO
    ELSE IF (present(fill)) THEN
       DO k=1, slv_u
          a(slv_l+ksize+k) = fill
       END DO
    END IF

  END SUBROUTINE parallel
#endif
  END SUBROUTINE update_boundary_z_r8

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE update_boundary_3d_logical(a)
    LOGICAL, INTENT(INOUT) :: a(:,:,:)

    REAL(4) :: x(size(a,1),size(a,2),size(a,3))
    INTEGER :: i, j, k

!$OMP PARALLEL DO
    DO k=1, size(a,3)
    DO j=1, size(a,2)
    DO i=1, size(a,1)
       IF (a(i,j,k)) THEN
          x(i,j,k) = 1.0
       ELSE
          x(i,j,k) = 0.0
       END IF
    END DO
    END DO
    END DO

    CALL update_boundary_3d(x)

!$OMP PARALLEL DO
    DO k=1, size(a,3)
    DO j=1, size(a,2)
    DO i=1, size(a,1)
       a(i,j,k) = (x(i,j,k) == 1.0)
    END DO
    END DO
    END DO

  END SUBROUTINE update_boundary_3d_logical

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE update_boundary_2d_logical(a)
    LOGICAL, INTENT(INOUT) :: a(:,:)

    REAL(4) :: x(size(a,1),size(a,2))
    INTEGER :: i, j

!$OMP PARALLEL DO
    DO j=1, size(a,2)
    DO i=1, size(a,1)
       IF (a(i,j)) THEN
          x(i,j) = 1.0
       ELSE
          x(i,j) = 0.0
       END IF
    END DO
    END DO

    CALL update_boundary_2d(x)

!$OMP PARALLEL DO
    DO j=1, size(a,2)
    DO i=1, size(a,1)
       a(i,j) = (x(i,j) == 1.0)
    END DO
    END DO

  END SUBROUTINE update_boundary_2d_logical

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE update_boundary_xz_logical(a)
    LOGICAL, INTENT(INOUT) :: a(:,:)

    REAL(4) :: x(size(a,1),size(a,2))
    INTEGER :: i, k

!$OMP PARALLEL DO
    DO k=1, size(a,2)
    DO i=1, size(a,1)
       IF (a(i,k)) THEN
          x(i,k) = 1.0
       ELSE
          x(i,k) = 0.0
       END IF
    END DO
    END DO

    CALL update_boundary_xz(x)

!$OMP PARALLEL DO
    DO k=1, size(a,2)
    DO i=1, size(a,1)
       a(i,k) = (x(i,k) == 1.0)
    END DO
    END DO

  END SUBROUTINE update_boundary_xz_logical

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE update_boundary_yz_logical(a)
    LOGICAL, INTENT(INOUT) :: a(:,:)

    REAL(4) :: x(size(a,1),size(a,2))
    INTEGER :: j, k

!$OMP PARALLEL DO
    DO k=1, size(a,2)
    DO j=1, size(a,1)
       IF (a(j,k)) THEN
          x(j,k) = 1.0
       ELSE
          x(j,k) = 0.0
       END IF
    END DO
    END DO

    CALL update_boundary_yz(x)

!$OMP PARALLEL DO
    DO k=1, size(a,2)
    DO j=1, size(a,1)
       a(j,k) = (x(j,k) == 1.0)
    END DO
    END DO

  END SUBROUTINE update_boundary_yz_logical

!-----------------------------------------------------------------------------------------------------------------------

  LOGICAL PURE FUNCTION check_dollar(filename)
    CHARACTER(*), INTENT(IN)  :: filename
    INTEGER :: i

    DO i=len_trim(filename), 1, -1
       IF (filename(i:i) == '/') EXIT
    END DO

    check_dollar = filename(i+1:i+1) == '$'
  END FUNCTION check_dollar

!-----------------------------------------------------------------------------------------------------------------------

  REAL(8) FUNCTION read_literal(str)
    CHARACTER(*), INTENT(IN)  :: str
    INTEGER :: i, iostat

    IF (trim(str) == 'UNDEF') THEN
       read_literal = UNDEF
       RETURN
    END IF

    READ(str, *, IOSTAT=iostat) read_literal
    CALL assert(iostat == 0, "invalid numeric lieteral '"//trim(str)//"'")

    IF(index(str, "E")/=0 .OR. index(str, "e")/=0) read_literal = REAL(read_literal, 4) !precision conversion

  END FUNCTION read_literal

!-----------------------------------------------------------------------------------------------------------------------

  CHARACTER(1024) PURE FUNCTION truncate_atmark(filename)
    CHARACTER(*), INTENT(IN)  :: filename
    INTEGER :: i, j

    i = index(filename, '@', back=.TRUE.)
    IF (i > 1) THEN
       truncate_atmark = trim(filename(1:i-1))
       j = index(filename(i:), '.')
       IF (j /= 0) truncate_atmark = trim(truncate_atmark) // trim(filename(j:))
    ELSE
       truncate_atmark = trim(filename)
    END IF
  END FUNCTION truncate_atmark

  CHARACTER(4) PURE FUNCTION atmark_code(filename)
    CHARACTER(*), INTENT(IN)  :: filename
    INTEGER :: i, j

    i = index(filename, '@', back=.TRUE.)
    IF (i == 0) THEN
       atmark_code = ''
       RETURN
    END IF

    j = index(filename(i+1:), '.')
    IF (j==0) THEN
         atmark_code = filename(i:)
      ELSE
         atmark_code = filename(i:i+j-1)
      END IF
    END FUNCTION atmark_code

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE read_data_2d_r4(data, filename, kind, region, view)
    REAL(4),      INTENT(OUT) :: data(isize, jsize)
    CHARACTER(*), INTENT(IN)  :: filename
    INTEGER,      INTENT(IN), OPTIONAL :: kind
    INTEGER,      INTENT(IN), OPTIONAL :: region
    INTEGER,      INTENT(IN), OPTIONAL :: view

    REAL(8) :: tmp(isize, jsize)

    CALL read_data_2d_r8(tmp, filename, kind, region, view)

    data(:,:) = REAL(tmp, KIND=4)

  END SUBROUTINE read_data_2d_r4

  SUBROUTINE read_data_3d_r4(data, filename, kind, region, view, descend)
    REAL(4),      INTENT(OUT) :: data(isize, jsize, ksize)
    CHARACTER(*), INTENT(IN)  :: filename
    INTEGER,      INTENT(IN), OPTIONAL :: kind
    INTEGER,      INTENT(IN), OPTIONAL :: region
    INTEGER,      INTENT(IN), OPTIONAL :: view
    LOGICAL,      INTENT(IN), OPTIONAL :: descend

    REAL(8) :: tmp(isize, jsize, ksize)

    CALL read_data_3d_r8(tmp, filename, kind, region, view, descend)

    data(:,:,:) = REAL(tmp, KIND=4)

  END SUBROUTINE read_data_3d_r4

  RECURSIVE SUBROUTINE read_data_xz_r4(data, filename, kind, region, descend)
    REAL(4),      INTENT(OUT) :: data(isize, ksize)
    CHARACTER(*), INTENT(IN)  :: filename
    INTEGER,      INTENT(IN), OPTIONAL :: kind
    INTEGER,      INTENT(IN), OPTIONAL :: region
    LOGICAL,      INTENT(IN), OPTIONAL :: descend

    REAL(8) :: tmp(isize, ksize)

    CALL read_data_xz_r8(tmp, filename, kind, region, descend)

    data(:,:) = REAL(tmp, KIND=4)

  END SUBROUTINE read_data_xz_r4

  RECURSIVE SUBROUTINE read_data_yz_r4(data, filename, kind, region, descend)
    REAL(4),      INTENT(OUT) :: data(jsize, ksize)
    CHARACTER(*), INTENT(IN)  :: filename
    INTEGER,      INTENT(IN), OPTIONAL :: kind
    INTEGER,      INTENT(IN), OPTIONAL :: region
    LOGICAL,      INTENT(IN), OPTIONAL :: descend

    REAL(8) :: tmp(jsize, ksize)

    CALL read_data_yz_r8(tmp, filename, kind, region, descend)

    data(:,:) = REAL(tmp, KIND=4)

  END SUBROUTINE read_data_yz_r4

  RECURSIVE SUBROUTINE read_data_x_r4(data, filename, kind, region)
    REAL(4),      INTENT(OUT) :: data(isize)
    CHARACTER(*), INTENT(IN)  :: filename
    INTEGER,      INTENT(IN), OPTIONAL :: kind
    INTEGER,      INTENT(IN), OPTIONAL :: region

    REAL(8) :: tmp(isize)

    CALL read_data_x_r8(tmp, filename, kind, region)

    data(:) = REAL(tmp, KIND=4)

  END SUBROUTINE read_data_x_r4

  RECURSIVE SUBROUTINE read_data_y_r4(data, filename, kind, region)
    REAL(4),      INTENT(OUT) :: data(jsize)
    CHARACTER(*), INTENT(IN)  :: filename
    INTEGER,      INTENT(IN), OPTIONAL :: kind
    INTEGER,      INTENT(IN), OPTIONAL :: region

    REAL(8) :: tmp(jsize)

    CALL read_data_y_r8(tmp, filename, kind, region)

    data(:) = REAL(tmp, KIND=4)

  END SUBROUTINE read_data_y_r4

  RECURSIVE SUBROUTINE read_data_z_r4(data, filename, kind, region, descend)
    REAL(4),      INTENT(OUT) :: data(ksize)
    CHARACTER(*), INTENT(IN)  :: filename
    INTEGER,      INTENT(IN), OPTIONAL :: kind
    INTEGER,      INTENT(IN), OPTIONAL :: region
    LOGICAL,      INTENT(IN), OPTIONAL :: descend

    REAL(8) :: tmp(ksize)

    CALL read_data_z_r8(tmp, filename, kind, region, descend)

    data(:) = REAL(tmp, KIND=4)

  END SUBROUTINE read_data_z_r4

!-----------------------------------------------------------------------------------------------------------------------

  RECURSIVE SUBROUTINE read_data_2d_r8(data, filename, kind, region, view)
    REAL(8),      INTENT(OUT) :: data(isize, jsize)
    CHARACTER(*), INTENT(IN)  :: filename
    INTEGER,      INTENT(IN), OPTIONAL :: kind
    INTEGER,      INTENT(IN), OPTIONAL :: region
    INTEGER,      INTENT(IN), OPTIONAL :: view

    INTEGER :: istr, iend
    INTEGER :: jstr, jend

    INTEGER :: kind_, region_, view_

    REAL(8) :: d

    kind_ = 8
    IF (present(kind)) kind_ = kind
    CALL assert(kind_==1 .OR. kind_==4 .OR. kind_==8, "unsupported KIND in READ_DATA_2D")

    region_ = 0
    IF (present(region)) region_ = region
    CALL assert(regions(region_)%defined, "undefined SUBREGION ID in READ_DATA_2D")

    view_ = 0
    IF (present(view)) view_ = view

#ifdef MPIIO
    CALL assert(views(view_)%defined,     "undefined FILEVIEW ID in READ_DATA_2D")
    CALL assert(view_==0 .OR. region_==0, "FILEVIEW is not supported for SUBREGION input")
#else
    CALL assert(view_==0,                 "FILEVIEW is supported only on MPIIO environment")
#endif

    istr = max(1, regions(region_)%x_start-icoord*isize)
    jstr = max(1, regions(region_)%y_start-jcoord*jsize)
    iend = min(isize, regions(region_)%x_end-icoord*isize)
    jend = min(jsize, regions(region_)%y_end-jcoord*jsize)

    IF (region_ /= 0) THEN
!$OMP PARALLEL WORKSHARE
       data(:,:) = UNDEF
!$OMP END PARALLEL WORKSHARE
    END IF

    IF (check_dollar(filename)) THEN
       IF (.NOT. in_region_2d(regions(region_))) RETURN

       d = read_literal(trim(filename(index(filename,'$',back=.TRUE.)+1:)))
!$OMP PARALLEL WORKSHARE
       data(istr:iend,jstr:jend) = d
!$OMP END PARALLEL WORKSHARE
       RETURN
    END IF

    IF (index(basename(filename), '@') /= 0) THEN
       CALL read_atmark_data
       RETURN
    END IF

PROFILE_BEGIN('ioread')
#ifdef PARALLEL_MPI
    IF (region_ /= 0) THEN
       CALL read_data_2d_subregion(data, filename, kind_, region_)
    ELSE
       SELECT CASE(trim(input_method))
#ifdef MPIIO
       CASE ('MPI-IO', 'MPIIO', 'mpi-io', 'mpiio')
          CALL read_data_2d_mpiio(data, filename, kind_, view_)
#endif
       CASE ('LAYERED', 'layered')
          CALL read_data_2d_layered(data, filename, kind_)
       CASE ('SERIAL', 'serial')
          CALL read_data_2d_serial(data, filename, kind_)
       CASE DEFAULT
          CALL assert(.FALSE., "parallel input method '"//trim(input_method)//"' is not supported")
       END SELECT

       IF (assert_nan) CALL assert(check_nan(data), "NaN detected in file '"//trim(filename)//"'")
    END IF
#else
    CALL read_file(data(regions(region_)%x_start:regions(region_)%x_end, &
                        regions(region_)%y_start:regions(region_)%y_end), filename, kind_)
#endif
PROFILE_END('ioread')

    IF (rank==0) WRITE(REPORT_UNIT, *) "read '"//trim(filename)//"'"

  CONTAINS
    SUBROUTINE read_atmark_data
      REAL(8), ALLOCATABLE :: tmp(:)
      CHARACTER(4) :: a
      INTEGER :: i, j

      a = atmark_code(filename)
      CALL assert(trim(a) /= '', "illegal call of READ_ATMARK_DATA")

      SELECT CASE (trim(a))
      CASE ('@x', '@X')
         ALLOCATE(tmp(regions(region_)%x_start:regions(region_)%x_end))
         IF (rank==0) CALL read_file(tmp, truncate_atmark(filename), kind_)
         CALL bcast(tmp)

         IF (in_region_2d(regions(region_))) THEN
!$OMP PARALLEL DO
            DO j=jstr, jend
            DO i=istr, iend
               data(i,j) = tmp(icoord*isize + i)
            END DO
            END DO
         END IF
         DEALLOCATE(tmp)
         IF (rank==0) WRITE(REPORT_UNIT, *) "read '"//trim(truncate_atmark(filename))//"'"

      CASE ('@y', '@Y')
         ALLOCATE(tmp(regions(region_)%y_start:regions(region_)%y_end))
         IF (rank==0) CALL read_file(tmp, truncate_atmark(filename), kind_)
         CALL bcast(tmp)

         IF (in_region_2d(regions(region_))) THEN
!$OMP PARALLEL DO
            DO j=jstr, jend
            DO i=istr, iend
               data(i,j) = tmp(jcoord*jsize + j)
            END DO
            END DO
         END IF
         DEALLOCATE(tmp)
         IF (rank==0) WRITE(REPORT_UNIT, *) "read '"//trim(truncate_atmark(filename))//"'"

      CASE ('@xy', '@XY')
         CALL read_data_2d(data, truncate_atmark(filename), kind_, region_, view_)

      CASE DEFAULT
         CALL assert(.FALSE., "invalid at-mark usage for '"//trim(filename)//"' in READ_DATA_2D")

      END SELECT

    END SUBROUTINE read_atmark_data

  END SUBROUTINE read_data_2d_r8

!-----------------------------------------------------------------------------------------------------------------------

  RECURSIVE SUBROUTINE read_data_3d_r8(data, filename, kind, region, view, descend)
    REAL(8),      INTENT(OUT) :: data(isize, jsize, ksize)
    CHARACTER(*), INTENT(IN)  :: filename
    INTEGER,      INTENT(IN), OPTIONAL :: kind
    INTEGER,      INTENT(IN), OPTIONAL :: region
    INTEGER,      INTENT(IN), OPTIONAL :: view
    LOGICAL,      INTENT(IN), OPTIONAL :: descend

    INTEGER :: kind_, region_, view_
    LOGICAL :: descend_

    INTEGER :: istr, iend
    INTEGER :: jstr, jend
    INTEGER :: kstr, kend

    REAL(8) :: d

    kind_ = 8
    IF (present(kind)) kind_ = kind
    CALL assert(kind_==1 .OR. kind_==4 .OR. kind_==8, "unsupported KIND in READ_DATA_3D")

    region_ = 0
    IF (present(region)) region_ = region
    CALL assert(regions(region_)%defined, "undefined subregion-code in READ_DATA_3D")

    view_ = 0
    IF (present(view)) view_ = view

    descend_ = .FALSE.
    IF (present(descend)) descend_ = descend


#ifdef MPIIO
    CALL assert(views(view_)%defined,     "undefined FILEVIEW ID in READ_DATA_2D")
    CALL assert(view_==0 .OR. region_==0, "FILEVIEW is not supported for SUBREGION input")
#else
    CALL assert(view_==0,                 "FILEVIEW is supported only on MPIIO environment")
#endif

    istr = max(1, regions(region_)%x_start-icoord*isize)
    jstr = max(1, regions(region_)%y_start-jcoord*jsize)
    kstr = max(1, regions(region_)%z_start-kcoord*ksize)
    iend = min(isize, regions(region_)%x_end-icoord*isize)
    jend = min(jsize, regions(region_)%y_end-jcoord*jsize)
    kend = min(ksize, regions(region_)%z_end-kcoord*ksize)

    IF (region_ /= 0) THEN
!$OMP PARALLEL WORKSHARE
       data(:,:,:) = UNDEF
!$OMP END PARALLEL WORKSHARE
    END IF

    IF (check_dollar(filename)) THEN
       IF (.NOT. in_region_3d(region_)) RETURN

       d = read_literal(trim(filename(index(filename,'$',back=.TRUE.)+1:)))
!$OMP PARALLEL WORKSHARE
       data(istr:iend,jstr:jend,kstr:kend) = d
!$OMP END PARALLEL WORKSHARE
       RETURN
    END IF

    IF (index(basename(filename), '@') /= 0) THEN
       CALL read_atmark_data
       RETURN
    END IF

PROFILE_BEGIN('ioread')

#ifdef PARALLEL_MPI
    IF (region_ /= 0) THEN
       CALL read_data_3d_subregion(data, filename, kind_, region_, descend_)
    ELSE
       SELECT CASE(trim(input_method))
#ifdef MPIIO
       CASE ('MPI-IO', 'MPIIO', 'mpi-io', 'mpiio')
          CALL read_data_3d_mpiio(data, filename, kind_, view_, descend_)
#endif
       CASE ('LAYERED', 'layered')
          CALL read_data_3d_layered(data, filename, kind_, descend_)
       CASE ('SERIAL', 'serial')
          CALL read_data_3d_serial(data, filename, kind_, descend_)
       CASE DEFAULT
          CALL assert(.FALSE., "parallel input method '"//trim(input_method)//"' is not supported")
       END SELECT

       IF (assert_nan) CALL assert(check_nan(data), "NaN detected in file '"//trim(filename)//"'")
    END IF
#else
    CALL read_file(data, filename, kind_)
    IF (descend_) CALL kreverse(data)
#endif
PROFILE_END('ioread')

    IF (rank==0) WRITE(REPORT_UNIT, *) "read '"//trim(filename)//"'"

  CONTAINS
    SUBROUTINE read_atmark_data
      REAL(8), ALLOCATABLE :: tmp(:,:)
      CHARACTER(4) :: a
      INTEGER :: i, j, k

      a = atmark_code(filename)
      CALL assert(trim(a) /= '', "illegal call of READ_ATMARK_DATA")

      SELECT CASE (trim(a))
      CASE ('@x', '@X')
         ALLOCATE(tmp(regions(region_)%x_start:regions(region_)%x_end,1))

         IF (rank==0) THEN
            CALL read_file(tmp(:,1), truncate_atmark(filename), kind_)
            IF (descend_) CALL kreverse(tmp(:,1))
         END IF
         CALL bcast(tmp(:,1))

         IF (in_region_3d(region_)) THEN
!$OMP PARALLEL DO
            DO k=kstr, kend
            DO j=jstr, jend
            DO i=istr, iend
               data(i,j,k) = tmp(icoord*isize + i,1)
            END DO
            END DO
            END DO
         END IF
         DEALLOCATE(tmp)
         IF (rank==0) WRITE(REPORT_UNIT, *) "read '"//trim(truncate_atmark(filename))//"'"

      CASE ('@y', '@Y')
         ALLOCATE(tmp(regions(region_)%y_start:regions(region_)%y_end,1))
         IF (rank==0) CALL read_file(tmp(:,1), truncate_atmark(filename), kind_)
         CALL bcast(tmp(:,1))

         IF (in_region_3d(region_)) THEN
!$OMP PARALLEL DO
            DO k=kstr, kend
            DO j=jstr, jend
            DO i=istr, iend
               data(i,j,k) = tmp(jcoord*jsize + j,1)
            END DO
            END DO
            END DO
         END IF
         DEALLOCATE(tmp)
         IF (rank==0) WRITE(REPORT_UNIT, *) "read '"//trim(truncate_atmark(filename))//"'"

      CASE ('@z', '@Z')
         ALLOCATE(tmp(regions(region_)%z_start:regions(region_)%z_end,1))

         IF (rank==0) THEN
            CALL read_file(tmp, truncate_atmark(filename), kind_)
            IF (descend_) CALL kreverse(tmp(:,1))
         END IF
         CALL bcast(tmp(:,1))

         IF (in_region_3d(region_)) THEN
!$OMP PARALLEL DO
            DO k=kstr, kend
            DO j=jstr, jend
            DO i=istr, iend
               data(i,j,k) = tmp(kcoord*ksize + k,1)
            END DO
            END DO
            END DO
         END IF
         DEALLOCATE(tmp)
         IF (rank==0) WRITE(REPORT_UNIT, *) "read '"//trim(truncate_atmark(filename))//"'"

      CASE ('@xy', '@XY')
         ALLOCATE(tmp(isize,jsize))
         CALL read_data_2d(tmp, truncate_atmark(filename), kind_, region_, view_)

         IF (in_region_3d(region_)) THEN
!$OMP PARALLEL DO
            DO k=kstr, kend
            DO j=jstr, jend
            DO i=istr, iend
               data(i,j,k) = tmp(i,j)
            END DO
            END DO
            END DO
         END IF
         DEALLOCATE(tmp)

      CASE ('@xz', '@XZ')
         ALLOCATE(tmp(isize,ksize))
         CALL read_data_xz(tmp, truncate_atmark(filename), kind_, region_, descend_)

         IF (in_region_3d(region_)) THEN
!$OMP PARALLEL DO
            DO k=kstr, kend
            DO j=jstr, jend
            DO i=istr, iend
               data(i,j,k) = tmp(i,k)
            END DO
            END DO
            END DO
         END IF
         DEALLOCATE(tmp)

      CASE ('@yz', '@YZ')
         ALLOCATE(tmp(jsize,ksize))
         CALL read_data_yz(tmp, truncate_atmark(filename), kind_, region_, descend_)

         IF (in_region_3d(region_)) THEN
!$OMP PARALLEL DO
            DO k=kstr, kend
            DO j=jstr, jend
            DO i=istr, iend
               data(i,j,k) = tmp(j,k)
            END DO
            END DO
            END DO
         END IF
         DEALLOCATE(tmp)

      CASE ('@xyz', '@XYZ')
         CALL read_data_3d(data, truncate_atmark(filename), kind_, region_, view_, descend_)

      CASE DEFAULT
         CALL assert(.FALSE., "invalid at-mark usage for '"//trim(filename)//"' in READ_DATA_3D")
      END SELECT

      IF (ALLOCATED(tmp)) DEALLOCATE(tmp)

    END SUBROUTINE read_atmark_data

  END SUBROUTINE read_data_3d_r8

!-----------------------------------------------------------------------------------------------------------------------

#ifdef PARALLEL_MPI
#ifdef MPIIO

  RECURSIVE SUBROUTINE read_data_2d_mpiio(data, filename, kind, view)
    REAL(8),      INTENT(OUT) :: data(isize, jsize)
    CHARACTER(*), INTENT(IN)  :: filename
    INTEGER,      INTENT(IN)  :: kind
    INTEGER,      INTENT(IN)  :: view

    LOGICAL :: exist
    INTEGER :: subarray, ierr

    CALL assert(kind==1 .OR. kind==4 .OR. kind==8, "unsupported KIND in READ_DATA_2D_MPIIO")

    IF (rank==0) THEN
       INQUIRE(FILE=trim(filename), EXIST=exist)
       CALL assert(exist, "file '"//trim(filename)//"' does not exist")
    END IF


    IF (vrank==0) THEN
       CALL mpi_file_open(hcomm, trim(filename), MPI_MODE_RDONLY, MPI_INFO_NULL, file%handle, ierr)
       CALL assert(ierr==MPI_SUCCESS, "failed to open '"//trim(filename)//"'")

       SELECT CASE (kind)
       CASE (1)
          subarray = views(view)%subarray_2d_i1
       CASE (4)
          subarray = views(view)%subarray_2d_r4
       CASE (8)
          subarray = views(view)%subarray_2d_r8
       END SELECT

       CALL mpi_file_set_view(file%handle, views(view)%offset, MPI_BYTE, subarray, 'native', MPI_INFO_NULL, ierr)
       CALL assert(ierr==MPI_SUCCESS, "MPI_FILE_SET_VIES failed")

       ! SX-Aurora@RIAM stalls at the following MPI_FILE_READ_ALL  (2022/3/17, Y. Matsumura and H. Tsuji)
       CALL mpi_file_read_all(file%handle, file%buffer, isize*jsize*kind, MPI_BYTE, MPI_STATUS_IGNORE, ierr)
       CALL assert(ierr==MPI_SUCCESS, "MPI_FILE_READ_ALL failed")

       CALL mpi_file_close(file%handle, ierr)

       IF (.NOT. check_endian()) CALL convert_endian(file%buffer, isize*jsize, kind)

       SELECT CASE (kind)
       CASE (1)
          data(1:isize,1:jsize) = REAL(reshape(transfer(file%buffer,   0_1, SIZE=isize*jsize), (/isize, jsize/)), KIND=8)
       CASE (4)
          data(1:isize,1:jsize) = REAL(reshape(transfer(file%buffer, 0.0_4, SIZE=isize*jsize), (/isize, jsize/)), KIND=8)
       CASE (8)
          data(1:isize,1:jsize) = REAL(reshape(transfer(file%buffer, 0.0_8, SIZE=isize*jsize), (/isize, jsize/)), KIND=8)
       END SELECT
    END IF

    CALL vcast(data)

  END SUBROUTINE read_data_2d_mpiio

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE read_data_3d_mpiio(data, filename, kind, view, descend)
    REAL(8),      INTENT(OUT) :: data(isize, jsize, ksize)
    CHARACTER(*), INTENT(IN)  :: filename
    INTEGER,      INTENT(IN)  :: kind
    INTEGER,      INTENT(IN)  :: view
    LOGICAL,      INTENT(IN)  :: descend

    LOGICAL :: exist
    INTEGER :: subarray, ierr

    CALL assert(kind==1 .OR. kind==4 .OR. kind==8, "unsupported KIND in READ_DATA_3D_MPIIO")

    IF (rank==0) THEN
       INQUIRE(FILE=trim(filename), EXIST=exist)
       CALL assert(exist, "file '"//trim(filename)//"' does not exist")
    END IF

    CALL mpi_file_open(comm, trim(filename), MPI_MODE_RDONLY, MPI_INFO_NULL, file%handle, ierr)
    CALL assert(ierr==MPI_SUCCESS, "failed to open '"//trim(filename)//"'")

    SELECT CASE (kind)
    CASE (1)
       IF (descend) THEN
          subarray = views(view)%subarray_3d_i1_desc
       ELSE
          subarray = views(view)%subarray_3d_i1
       END IF
    CASE (4)
       IF (descend) THEN
          subarray = views(view)%subarray_3d_r4_desc
       ELSE
          subarray = views(view)%subarray_3d_r4
       END IF
    CASE (8)
       IF (descend) THEN
          subarray = views(view)%subarray_3d_r8_desc
       ELSE
          subarray = views(view)%subarray_3d_r8
       END IF
    END SELECT

    CALL mpi_file_set_view(file%handle, views(view)%offset, MPI_BYTE, subarray, 'native', MPI_INFO_NULL, ierr)
    CALL assert(ierr==MPI_SUCCESS, "MPI_FILE_SET_VIEW failed")

    CALL mpi_file_read_all(file%handle, file%buffer, isize*jsize*ksize*kind, MPI_BYTE, MPI_STATUS_IGNORE, ierr)
    CALL assert(ierr==MPI_SUCCESS, "MPI_FILE_READ_ALL failed")

    CALL mpi_file_close(file%handle, ierr)

    IF (.NOT. check_endian()) CALL convert_endian(file%buffer, isize*jsize*ksize, kind)

    SELECT CASE (kind)
    CASE (1)
       data(:,:,:) = REAL(reshape(transfer(file%buffer,   0_1, SIZE=isize*jsize*ksize), (/isize,jsize,ksize/)), KIND=8)
    CASE (4)
       data(:,:,:) = REAL(reshape(transfer(file%buffer, 0.0_4, SIZE=isize*jsize*ksize), (/isize,jsize,ksize/)), KIND=8)
    CASE (8)
       data(:,:,:) = REAL(reshape(transfer(file%buffer, 0.0_8, SIZE=isize*jsize*ksize), (/isize,jsize,ksize/)), KIND=8)
    END SELECT

    IF (descend) CALL kreverse(data)

  END SUBROUTINE read_data_3d_mpiio
#endif

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE read_data_2d_serial(data, filename, kind)
    REAL(8),      INTENT(OUT) :: data(isize, jsize)
    CHARACTER(*), INTENT(IN)  :: filename
    INTEGER,      INTENT(IN)  :: kind

    REAL(8) :: buf(isize, jsize)

    INTEGER(1) :: tmp1(isize)
    REAL(4)    :: tmp4(isize)
    REAL(8)    :: tmp8(isize)

    INTEGER :: iostat, ierr, rec

    INTEGER :: n, i, j

    CALL assert(kind==1 .OR. kind==4 .OR. kind==8, "unsupported KIND in READ_DATA_2D_SERIAL")

    IF (vrank==0) THEN
       IF (hrank==0) THEN
          OPEN(UNIT   = TMP_UNIT,      &
               FILE   = trim(filename),&
               FORM   = 'UNFORMATTED', &
               ACCESS = 'DIRECT',      &
               STATUS = 'OLD',         &
               ACTION = 'READ',        &
               RECL   = isize*kind,    &
               IOSTAT = iostat)
          CALL assert(iostat==0, "failed to open '"//trim(filename)//"'")

          DO n=0, hpes-1
             DO j=1, jsize
                rec = (jcoords_h(n)*jsize+j-1)*ipes + icoords_h(n) + 1

                SELECT CASE (kind)
                CASE (1)
                   READ(TMP_UNIT, REC=rec) tmp1(:)
                   buf(:,j) = REAL(tmp1(:), KIND=8)
                CASE (4)
                   READ(TMP_UNIT, REC=rec) tmp4(:)
                   buf(:,j) = REAL(tmp4(:), KIND=8)
                CASE (8)
                   READ(TMP_UNIT, REC=rec) tmp8(:)
                   buf(:,j) = REAL(tmp8(:), KIND=8)
                END SELECT
             END DO

             IF (n > 0)  THEN
                CALL mpi_send(buf, isize*jsize, MPI_REAL8, n, 0, hcomm, ierr)
             ELSE
                data(:,:) = buf(:,:)
             END IF
          END DO
          CLOSE(TMP_UNIT)
       ELSE
          CALL mpi_recv(data, isize*jsize, MPI_REAL8, 0, 0, hcomm, MPI_STATUS_IGNORE, ierr)
       END IF
    END IF

    CALL vcast(data)

  END SUBROUTINE read_data_2d_serial

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE read_data_3d_serial(data, filename, kind, descend)
    REAL(8),      INTENT(OUT) :: data(isize, jsize, ksize)
    CHARACTER(*), INTENT(IN)  :: filename
    INTEGER,      INTENT(IN)  :: kind
    LOGICAL,      INTENT(IN)  :: descend

    REAL(8) :: buf(isize, jsize, ksize)

    INTEGER(1) :: tmp1(isize)
    REAL(4)    :: tmp4(isize)
    REAL(8)    :: tmp8(isize)

    INTEGER :: iostat, ierr, rec

    INTEGER :: n, i, j, k

    CALL assert(kind==1 .OR. kind==4 .OR. kind==8, "unsupported KIND in READ_DATA_3D_SERIAL")

    IF (rank==0) THEN
       OPEN(UNIT   = TMP_UNIT,      &
            FILE   = trim(filename),&
            FORM   = 'UNFORMATTED', &
            ACCESS = 'DIRECT',      &
            STATUS = 'OLD',         &
            ACTION = 'READ',        &
            RECL   = isize*kind,    &
            IOSTAT = iostat)

       CALL assert(iostat==0, "failed to open '"//trim(filename)//"'")

       DO n=0, npes-1
          DO k=1, ksize
          DO j=1, jsize
             IF (descend) THEN
                rec = dimz - (kcoords(n)*ksize+k)
             ELSE
                rec = kcoords(n)*ksize+k - 1
             END IF
             rec = rec*dimy*ipes + (jcoords(n)*jsize+j-1)*ipes + icoords(n) + 1

             SELECT CASE (kind)
             CASE (1)
                READ(TMP_UNIT, REC=rec) tmp1(:)
                buf(:,j,k) = REAL(tmp1(:), KIND=8)
             CASE (4)
                READ(TMP_UNIT, REC=rec) tmp4(:)
                buf(:,j,k) = REAL(tmp4(:), KIND=8)
             CASE (8)
                READ(TMP_UNIT, REC=rec) tmp8(:)
                buf(:,j,k) = REAL(tmp8(:), KIND=8)
             END SELECT
          END DO
          END DO

          IF (n > 0) THEN
             CALL mpi_send(buf, isize*jsize*ksize, MPI_REAL8, n, 0, comm, ierr)
          ELSE
             data(:,:,:) = buf(:,:,:)
          END IF
       END DO
       CLOSE(TMP_UNIT)
    ELSE
       CALL mpi_recv(data, isize*jsize*ksize, MPI_REAL8, 0, 0, comm, MPI_STATUS_IGNORE, ierr)
    END IF

    CALL mpi_barrier(comm, ierr)

  END SUBROUTINE read_data_3d_serial

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE read_data_2d_layered(data, filename, kind)
    REAL(8),      INTENT(OUT) :: data(isize, jsize)
    CHARACTER(*), INTENT(IN)  :: filename
    INTEGER,      INTENT(IN)  :: kind

    REAL(8) :: buf(isize, jsize, 0:hpes-1)

    INTEGER(1) :: tmp1(isize*ipes, jsize*jpes)
    REAL(4)    :: tmp4(isize*ipes, jsize*jpes)
    REAL(8)    :: tmp8(isize*ipes, jsize*jpes)

    INTEGER :: iostat, ierr, rec

    INTEGER :: n, i, j

    CALL assert(kind==1 .OR. kind==4 .OR. kind==8, "unsupported KIND in READ_DATA_2D_LAYERED")

    IF (vrank==0) THEN
       IF (hrank==0) THEN
          OPEN(UNIT   = TMP_UNIT,      &
               FILE   = trim(filename),&
               FORM   = 'UNFORMATTED', &
               ACCESS = 'DIRECT',      &
               STATUS = 'OLD',         &
               ACTION = 'READ',        &
               RECL   = isize*ipes*jsize*jpes*kind, &
               IOSTAT = iostat)

          CALL assert(iostat==0, "failed to open '"//trim(filename)//"'")

          SELECT CASE (kind)
          CASE (1)
             READ(TMP_UNIT, REC=1) tmp1
             DO n=0, hpes-1
                buf(:,:,n) = REAL(tmp1(icoords_h(n)*isize+1:(icoords_h(n)+1)*isize, &
                                       jcoords_h(n)*jsize+1:(jcoords_h(n)+1)*jsize), KIND=8)
             END DO
          CASE (4)
             READ(TMP_UNIT, REC=1) tmp4
             DO n=0, hpes-1
                buf(:,:,n) = REAL(tmp4(icoords_h(n)*isize+1:(icoords_h(n)+1)*isize, &
                                       jcoords_h(n)*jsize+1:(jcoords_h(n)+1)*jsize), KIND=8)
             END DO
          CASE (8)
             READ(TMP_UNIT, REC=1) tmp8
             DO n=0, hpes-1
                buf(:,:,n) = REAL(tmp8(icoords_h(n)*isize+1:(icoords_h(n)+1)*isize, &
                                       jcoords_h(n)*jsize+1:(jcoords_h(n)+1)*jsize), KIND=8)
             END DO
          END SELECT
          CLOSE(TMP_UNIT)
       END IF

       CALL mpi_scatter(buf, isize*jsize, MPI_REAL8, data, isize*jsize, MPI_REAL8, 0, hcomm, ierr)
    END IF

    CALL vcast(data)

  END SUBROUTINE read_data_2d_layered

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE read_data_2d_subregion(data, filename, kind, region)
    REAL(8),      INTENT(OUT) :: data(isize, jsize)
    CHARACTER(*), INTENT(IN)  :: filename
    INTEGER,      INTENT(IN)  :: kind
    INTEGER,      INTENT(IN)  :: region

    INTEGER(1), ALLOCATABLE :: tmp1(:,:)
    REAL(4),    ALLOCATABLE :: tmp4(:,:)
    REAL(8),    ALLOCATABLE :: tmp8(:,:)

    REAL(8) :: buf(isize, jsize, 0:hpes-1)

    TYPE(subregion_struct) :: rgn

    INTEGER :: iostat, ierr, rec

    INTEGER :: istr, iend
    INTEGER :: jstr, jend

    INTEGER :: n, i, j

#ifdef PARALLEL_MPI
    INTEGER :: req(0:hpes-1)
#endif

    rgn = regions(region)
    CALL assert(rgn%defined, "undefined subregion-code in READ_DATA_2D_SUBREGION")

    CALL assert(kind==1 .OR. kind==4 .OR. kind==8, "unsupported KIND in READ_DATA_2D_LAYERED")

    IF (vrank==0) THEN
       IF (hrank==0) THEN
          ALLOCATE(tmp8(rgn%x_start:rgn%x_end, rgn%y_start:rgn%y_end))

          OPEN(UNIT   = TMP_UNIT,      &
               FILE   = trim(filename),&
               FORM   = 'UNFORMATTED', &
               ACCESS = 'DIRECT',      &
               STATUS = 'OLD',         &
               ACTION = 'READ',        &
               RECL   = rgn%x_size*rgn%y_size*kind, &
               IOSTAT = iostat)

          CALL assert(iostat==0, "failed to open '"//trim(filename)//"'")

          SELECT CASE (kind)
          CASE (1)
             ALLOCATE(tmp1(rgn%x_start:rgn%x_end, rgn%y_start:rgn%y_end))
             READ(TMP_UNIT, REC=1) tmp1
             tmp8(:,:) = REAL(tmp1(:,:), KIND=8)
             DEALLOCATE(tmp1)
          CASE (4)
             ALLOCATE(tmp4(rgn%x_start:rgn%x_end, rgn%y_start:rgn%y_end))
             READ(TMP_UNIT, REC=1) tmp4
             tmp8(:,:) = REAL(tmp4(:,:), KIND=8)
             DEALLOCATE(tmp4)
          CASE (8)
             READ(TMP_UNIT, REC=1) tmp8
          END SELECT

          CLOSE(TMP_UNIT)

          req(:) = MPI_REQUEST_NULL

          DO n=0, hpes-1
             IF (.NOT. in_region_h(rgn, n)) CYCLE

             istr = max(1,     rgn%x_start-icoords_h(n)*isize)
             jstr = max(1,     rgn%y_start-jcoords_h(n)*jsize)
             iend = min(isize, rgn%x_end  -icoords_h(n)*isize)
             jend = min(jsize, rgn%y_end  -jcoords_h(n)*jsize)

!$OMP PARALLEL PRIVATE(i,j)
!$OMP WORKSHARE
             buf(:,:,n) = UNDEF
!$OMP END WORKSHARE

!$OMP DO
             DO j=jstr, jend
             DO i=istr, iend
                buf(i,j,n) = tmp8(isize*icoords_h(n)+i,jsize*jcoords_h(n)+j)
             END DO
             END DO
!$OMP END PARALLEL
             CALL mpi_isend(buf(:,:,n), isize*jsize, MPI_REAL8, n, 0, comm, req(n), ierr)
          END DO
       END IF

       IF (in_region_h(rgn)) CALL mpi_recv(data, isize*jsize, MPI_REAL8, 0, 0, comm, MPI_STATUS_IGNORE, ierr)
       IF (hrank==0) CALL mpi_waitall(hpes, req, MPI_STATUSES_IGNORE, ierr)
    END IF

    CALL vcast(data)

  END SUBROUTINE read_data_2d_subregion

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE read_data_3d_layered(data, filename, kind, descend)
    REAL(8),      INTENT(OUT) :: data(isize, jsize, ksize)
    CHARACTER(*), INTENT(IN)  :: filename
    INTEGER,      INTENT(IN)  :: kind
    LOGICAL,      INTENT(IN)  :: descend

    REAL(8) :: buf1(isize, jsize, 0:hpes-1)
    REAL(8) :: buf2(isize, jsize, ksize*kpes)
    REAL(8) :: buf3(isize, jsize, ksize, 0:vpes-1)

    INTEGER(1) :: tmp1(isize*ipes, jsize*jpes)
    REAL(4)    :: tmp4(isize*ipes, jsize*jpes)
    REAL(8)    :: tmp8(isize*ipes, jsize*jpes)

    INTEGER :: iostat, ierr, rec
    INTEGER :: n, i, j, k

    CALL assert(kind==1 .OR. kind==4 .OR. kind==8, "unsupported KIND in READ_DATA_3D_LAYERED")

    IF (vrank==0) THEN
       IF (hrank==0) THEN
          OPEN(UNIT   = TMP_UNIT,      &
               FILE   = trim(filename),&
               FORM   = 'UNFORMATTED', &
               ACCESS = 'DIRECT',      &
               STATUS = 'OLD',         &
               ACTION = 'READ',        &
               RECL   = isize*ipes*jsize*jpes*kind, &
               IOSTAT = iostat)
          CALL assert(iostat==0, "failed to open '"//trim(filename)//"'")
       END IF

       DO k=1, dimz
          IF (descend) THEN
             rec = dimz-k+1
          ELSE
             rec = k
          END IF

          IF (hrank==0) THEN

             SELECT CASE (kind)
             CASE (1)
                READ(TMP_UNIT, REC=rec) tmp1
                DO n=0, hpes-1
                   buf1(:,:,n) = REAL(tmp1(icoords_h(n)*isize+1:(icoords_h(n)+1)*isize, &
                                           jcoords_h(n)*jsize+1:(jcoords_h(n)+1)*jsize), KIND=8)
                END DO
             CASE (4)
                READ(TMP_UNIT, REC=rec) tmp4
                DO n=0, hpes-1
                   buf1(:,:,n) = REAL(tmp4(icoords_h(n)*isize+1:(icoords_h(n)+1)*isize, &
                                           jcoords_h(n)*jsize+1:(jcoords_h(n)+1)*jsize), KIND=8)
                END DO
             CASE (8)
                READ(TMP_UNIT, REC=rec) tmp8
                DO n=0, hpes-1
                   buf1(:,:,n) = tmp8(icoords_h(n)*isize+1:(icoords_h(n)+1)*isize, &
                                      jcoords_h(n)*jsize+1:(jcoords_h(n)+1)*jsize)
                END DO
             END SELECT
          END IF
          CALL mpi_scatter(buf1, isize*jsize, MPI_REAL8, buf2(:,:,k), isize*jsize, MPI_REAL8, 0, hcomm, ierr)
       END DO

       IF (hrank==0) CLOSE(TMP_UNIT)
    END IF

    IF (vrank==0) THEN
       DO n=0, vpes-1
          buf3(:,:,:,n) = buf2(:,:,kcoords_v(n)*ksize+1:(kcoords_v(n)+1)*ksize)
       END DO
    END IF
    CALL mpi_scatter(buf3, isize*jsize*ksize, MPI_REAL8, data(:,:,:), isize*jsize*ksize, MPI_REAL8, 0, vcomm, ierr)

  END SUBROUTINE read_data_3d_layered

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE read_data_3d_subregion(data, filename, kind, region, descend)
    REAL(8),      INTENT(OUT) :: data(isize, jsize, ksize)
    CHARACTER(*), INTENT(IN)  :: filename
    INTEGER,      INTENT(IN)  :: kind
    INTEGER,      INTENT(IN)  :: region
    LOGICAL,      INTENT(IN)  :: descend

    INTEGER(1), ALLOCATABLE :: tmp1(:,:)
    REAL(4),    ALLOCATABLE :: tmp4(:,:)
    REAL(8),    ALLOCATABLE :: tmp8(:,:)

    REAL(8) :: buf(isize, jsize, 0:hpes-1)

    TYPE(subregion_struct) :: rgn

    INTEGER :: iostat, ierr, rec

    INTEGER :: istr, iend
    INTEGER :: jstr, jend
    INTEGER :: kstr, kend

    INTEGER :: n, i, j, k

#ifdef PARALLEL_MPI
    INTEGER :: req(0:npes-1)
#endif
    rgn = regions(region)
    CALL assert(rgn%defined, "undefined subregion-code in READ_DATA_3D_SUBREGION")

    IF (rank==0) THEN
       ALLOCATE(tmp8(rgn%x_start:rgn%x_end, rgn%y_start:rgn%y_end))

       IF (kind==1) ALLOCATE(tmp1(rgn%x_size, rgn%y_size))
       IF (kind==4) ALLOCATE(tmp4(rgn%x_size, rgn%y_size))

       OPEN(UNIT   = TMP_UNIT,      &
            FILE   = trim(filename),&
            FORM   = 'UNFORMATTED', &
            ACCESS = 'DIRECT',      &
            STATUS = 'OLD',         &
            ACTION = 'READ',        &
            RECL   = rgn%x_size*rgn%y_size*kind, &
            IOSTAT = iostat)

       CALL assert(iostat==0, "failed to open '"//trim(filename)//"'")
    END IF

    DO k=rgn%z_start, rgn%z_end
       IF (rank==0) THEN
          IF (descend) THEN
             rec = rgn%z_end - k + 1
          ELSE
             rec = k - rgn%z_start + 1
          END IF

          SELECT CASE (kind)
          CASE (1)
             READ(TMP_UNIT, REC=rec) tmp1
             tmp8(:,:) = REAL(tmp1(:,:), KIND=8)
          CASE (4)
             READ(TMP_UNIT, REC=rec) tmp4
             tmp8(:,:) = REAL(tmp4(:,:), KIND=8)
          CASE (8)
             READ(TMP_UNIT, REC=rec) tmp8
          END SELECT

          req(:) = MPI_REQUEST_NULL

          DO n=0, npes-1
             IF (.NOT. in_region_3d(region, n)) CYCLE

             istr = max(1,     rgn%x_start-icoords(n)*isize)
             jstr = max(1,     rgn%y_start-jcoords(n)*jsize)
             iend = min(isize, rgn%x_end  -icoords(n)*isize)
             jend = min(jsize, rgn%y_end  -jcoords(n)*jsize)

!$OMP PARALLEL PRIVATE(i,j)
!$OMP WORKSHARE
             buf(:,:,n) = UNDEF
!$OMP END WORKSHARE

!$OMP DO
             DO j=jstr, jend
             DO i=istr, iend
                buf(i,j,n) = tmp8(isize*icoords(n)+i,jsize*jcoords(n)+j)
             END DO
             END DO
!$OMP END PARALLEL

             CALL mpi_isend(buf(:,:,n), isize*jsize, MPI_REAL8, n, k, comm, req(n), ierr)
          END DO
       END IF

       IF (in_region_3d(region)) CALL mpi_recv(data(:,:,k - kcoord*ksize), isize*jsize, MPI_REAL8, 0, k, comm, MPI_STATUS_IGNORE, ierr)
       IF (rank==0) CALL mpi_waitall(npes, req, MPI_STATUSES_IGNORE, ierr)
    END DO

    IF (rank==0) THEN
       CLOSE(TMP_UNIT)
       DEALLOCATE(tmp8)
       IF (kind==1) DEALLOCATE(tmp1)
       IF (kind==4) DEALLOCATE(tmp4)
    END IF

  END SUBROUTINE read_data_3d_subregion
#endif
!-----------------------------------------------------------------------------------------------------------------------

  RECURSIVE SUBROUTINE read_data_xz_r8(data, filename, kind, region, descend)
    REAL(8),      INTENT(OUT) :: data(isize, ksize)
    CHARACTER(*), INTENT(IN)  :: filename
    INTEGER,      INTENT(IN), OPTIONAL :: kind
    INTEGER,      INTENT(IN), OPTIONAL :: region
    LOGICAL,      INTENT(IN), OPTIONAL :: descend

    TYPE(subregion_struct) :: rgn
    INTEGER :: istr, iend
    INTEGER :: kstr, kend

    REAL(8), ALLOCATABLE :: tmp(:,:)

#ifdef PARALLEL_MPI
    REAL(8) :: buf(isize, ksize, 0:npes-1)
    INTEGER :: req(0:npes-1)
    INTEGER :: ierr
#endif

    LOGICAL :: descend_

    REAL(8) :: d

    INTEGER :: i, k, n

    IF (present(region)) THEN
       rgn = regions(region)
       CALL assert(rgn%defined, "undefined subregion-code in READ_DATA_XZ")

       IF (region /= 0) THEN
!$OMP PARALLEL WORKSHARE
          data(:,:) = UNDEF
!$OMP END PARALLEL WORKSHARE
       END IF
    ELSE
       rgn = regions(0)
    END IF

    descend_ = .FALSE.
    IF (present(descend)) descend_ = descend

    istr = max(1, rgn%x_start-icoord*isize)
    kstr = max(1, rgn%z_start-kcoord*ksize)
    iend = min(isize, rgn%x_end-icoord*isize)
    kend = min(ksize, rgn%z_end-kcoord*ksize)

    IF (check_dollar(filename)) THEN
       IF (.NOT. in_region_xz(rgn)) RETURN

       d = read_literal(trim(filename(index(filename,'$',back=.TRUE.)+1:)))
!$OMP PARALLEL WORKSHARE
       data(istr:iend,kstr:kend) = d
!$OMP END PARALLEL WORKSHARE
       RETURN
    END IF

    IF (index(basename(filename), '@') /= 0) THEN
       CALL read_atmark_data
       RETURN
    END IF

    ALLOCATE(tmp(rgn%x_start:rgn%x_end, rgn%z_start:rgn%z_size))

    IF (rank==0) THEN
       CALL read_file(tmp, filename, kind)
       IF (descend_) CALL kreverse(tmp)
    END IF

#ifdef PARALLEL_MPI
    IF (rank==0) THEN
       req(:) = MPI_REQUEST_NULL

       DO n=0, npes-1
          IF (.NOT. in_region_xz(rgn, n)) CYCLE

          istr = max(1, rgn%x_start-icoords(n)*isize)
          kstr = max(1, rgn%z_start-kcoords(n)*ksize)
          iend = min(isize, rgn%x_end-icoords(n)*isize)
          kend = min(ksize, rgn%z_end-kcoords(n)*ksize)

!$OMP PARALLEL PRIVATE(i,k)
!$OMP WORKSHARE
          buf(:,:,n) = UNDEF
!$OMP END WORKSHARE

!$OMP DO
          DO k=kstr, kend
          DO i=istr, iend
             buf(i,k,n) = tmp(isize*icoords(n)+i,ksize*kcoords(n)+k)
          END DO
          END DO
!$OMP END PARALLEL
          CALL mpi_isend(buf(:,:,n), isize*ksize, MPI_REAL8, n, 0, comm, req(n), ierr)
       END DO
    END IF

    IF (in_region_xz(rgn)) CALL mpi_recv(data, isize*ksize, MPI_REAL8, 0, 0, comm, MPI_STATUS_IGNORE, ierr)
    IF (rank==0) CALL mpi_waitall(npes, req, MPI_STATUSES_IGNORE, ierr)
#else
!$OMP PARALLEL DO
    DO k=rgn%z_start, rgn%z_end
    DO i=rgn%x_start, rgn%x_end
       data(i,k) = tmp(i,k)
    END DO
    END DO
#endif
    DEALLOCATE(tmp)

    IF (rank==0) WRITE(REPORT_UNIT, *) "read '"//trim(filename)//"'"

  CONTAINS
    SUBROUTINE read_atmark_data
      REAL(8), ALLOCATABLE :: tmp(:)
      CHARACTER(4) :: a
      INTEGER :: i, j

      a = atmark_code(filename)
      CALL assert(trim(a) /= '', "illegal call of READ_ATMARK_DATA")

      SELECT CASE (trim(a))
      CASE ('@x', '@X')
         ALLOCATE(tmp(rgn%x_start:rgn%x_end))
         IF (rank==0) CALL read_file(tmp, truncate_atmark(filename), kind)
         CALL bcast(tmp)

         IF (in_region_xz(rgn)) THEN
!$OMP PARALLEL DO
            DO k=kstr, kend
            DO i=istr, iend
               data(i,k) = tmp(icoord*isize + i)
            END DO
            END DO
         END IF
         DEALLOCATE(tmp)
         IF (rank==0) WRITE(REPORT_UNIT, *) "read '"//trim(truncate_atmark(filename))//"'"

      CASE ('@z', '@Z')
         ALLOCATE(tmp(rgn%z_start:rgn%z_end))
         IF (rank==0) THEN
            CALL read_file(tmp, truncate_atmark(filename), kind)
            IF (descend_) CALL kreverse(tmp)
         END IF
         CALL bcast(tmp)

         IF (in_region_xz(rgn)) THEN
!$OMP PARALLEL DO
            DO k=kstr, kend
            DO i=istr, iend
               data(i,k) = tmp(kcoord*ksize + k)
            END DO
            END DO
         END IF
         DEALLOCATE(tmp)
         IF (rank==0) WRITE(REPORT_UNIT, *) "read '"//trim(truncate_atmark(filename))//"'"

      CASE ('@xz', '@XZ')
         CALL read_data_xz(data, truncate_atmark(filename), kind, region, descend_)

      CASE DEFAULT
         CALL assert(.FALSE., "invalid at-mark usage for '"//trim(filename)//"' in READ_DATA_XZ")

      END SELECT

    END SUBROUTINE read_atmark_data

  END SUBROUTINE read_data_xz_r8

!-----------------------------------------------------------------------------------------------------------------------

  RECURSIVE SUBROUTINE read_data_yz_r8(data, filename, kind, region, descend)
    REAL(8),      INTENT(OUT) :: data(jsize, ksize)
    CHARACTER(*), INTENT(IN)  :: filename
    INTEGER,      INTENT(IN), OPTIONAL :: kind
    INTEGER,      INTENT(IN), OPTIONAL :: region
    LOGICAL,      INTENT(IN), OPTIONAL :: descend

    TYPE(subregion_struct) :: rgn
    INTEGER :: jstr, jend
    INTEGER :: kstr, kend

    REAL(8), ALLOCATABLE :: tmp(:,:)

#ifdef PARALLEL_MPI
    REAL(8) :: buf(jsize, ksize, 0:npes-1)
    INTEGER :: req(0:npes-1)
    INTEGER :: ierr
#endif

    LOGICAL :: descend_

    REAL(8) :: d
    INTEGER :: j, k, n

    IF (present(region)) THEN
       rgn = regions(region)
       CALL assert(rgn%defined, "undefined subregion-code in READ_DATA_YZ")

       IF (region /= 0) THEN
!$OMP PARALLEL WORKSHARE
          data(:,:) = UNDEF
!$OMP END PARALLEL WORKSHARE
       END IF
    ELSE
       rgn = regions(0)
    END IF

    descend_ = .FALSE.
    IF (present(descend)) descend_ = descend

    jstr = max(1, rgn%y_start-jcoord*jsize)
    kstr = max(1, rgn%z_start-kcoord*ksize)
    jend = min(jsize, rgn%y_end-jcoord*jsize)
    kend = min(ksize, rgn%z_end-kcoord*ksize)

    IF (check_dollar(filename)) THEN
       IF (.NOT. in_region_yz(rgn)) RETURN

       d = read_literal(trim(filename(index(filename,'$',back=.TRUE.)+1:)))
!$OMP PARALLEL WORKSHARE
       data(jstr:jend,kstr:kend) = d
!$OMP END PARALLEL WORKSHARE
       RETURN
    END IF

    IF (index(basename(filename), '@') /= 0) THEN
       CALL read_atmark_data
       RETURN
    END IF

    ALLOCATE(tmp(rgn%y_start:rgn%y_end, rgn%z_start:rgn%z_size))

    IF (rank==0) THEN
       CALL read_file(tmp, filename, kind)
       IF (descend_) CALL kreverse(tmp)
    END IF

#ifdef PARALLEL_MPI
    IF (rank==0) THEN
       req(:) = MPI_REQUEST_NULL

       DO n=0, npes-1
          IF (.NOT. in_region_yz(rgn, n)) CYCLE

          jstr = max(1, rgn%y_start-jcoords(n)*jsize)
          kstr = max(1, rgn%z_start-kcoords(n)*ksize)
          jend = min(jsize, rgn%y_end-jcoords(n)*jsize)
          kend = min(ksize, rgn%z_end-kcoords(n)*ksize)

!$OMP PARALLEL PRIVATE(j,k)
!$OMP WORKSHARE
          buf(:,:,n) = UNDEF
!$OMP END WORKSHARE

!$OMP DO
          DO k=kstr, kend
          DO j=jstr, jend
             buf(j,k,n) = tmp(jsize*jcoords(n)+j,ksize*kcoords(n)+k)
          END DO
          END DO
!$OMP END PARALLEL
          CALL mpi_isend(buf(:,:,n), jsize*ksize, MPI_REAL8, n, 0, comm, req(n), ierr)
       END DO
    END IF

    IF (in_region_yz(rgn)) CALL mpi_recv(data, jsize*ksize, MPI_REAL8, 0, 0, comm, MPI_STATUS_IGNORE, ierr)
    IF (rank==0) CALL mpi_waitall(npes, req, MPI_STATUSES_IGNORE, ierr)
#else

!$OMP PARALLEL DO
    DO k=rgn%z_start, rgn%z_end
    DO j=rgn%y_start, rgn%y_end
       data(j,k) = tmp(j,k)
    END DO
    END DO
#endif
    DEALLOCATE(tmp)

    IF (rank==0) WRITE(REPORT_UNIT, *) "read '"//trim(filename)//"'"

  CONTAINS
    SUBROUTINE read_atmark_data
      REAL(8), ALLOCATABLE :: tmp(:)
      CHARACTER(4) :: a
      INTEGER :: j, k

      a = atmark_code(filename)
      CALL assert(trim(a) /= '', "illegal call of READ_ATMARK_DATA")

      SELECT CASE (trim(a))
      CASE ('@y', '@Y')
         ALLOCATE(tmp(rgn%y_start:rgn%y_end))
         IF (rank==0) THEN
            CALL read_file(tmp, truncate_atmark(filename), kind)
         END IF
         CALL bcast(tmp)

         IF (in_region_yz(rgn)) THEN
!$OMP PARALLEL DO
            DO k=kstr, kend
            DO j=jstr, jend
               data(j,k) = tmp(jcoord*jsize + j)
            END DO
            END DO
         END IF
         DEALLOCATE(tmp)
         IF (rank==0) WRITE(REPORT_UNIT, *) "read '"//trim(truncate_atmark(filename))//"'"

      CASE ('@z', '@Z')
         ALLOCATE(tmp(rgn%z_start:rgn%z_end))
         IF (rank==0) THEN
            CALL read_file(tmp, truncate_atmark(filename), kind)
            IF (descend_) CALL kreverse(tmp)
         END IF
         CALL bcast(tmp)

         IF (in_region_yz(rgn)) THEN
!$OMP PARALLEL DO
            DO k=kstr, kend
            DO j=jstr, jend
               data(j,k) = tmp(kcoord*ksize + k)
            END DO
            END DO
         END IF
         DEALLOCATE(tmp)
         IF (rank==0) WRITE(REPORT_UNIT, *) "read '"//trim(truncate_atmark(filename))//"'"

      CASE ('@yz', '@YZ')
         CALL read_data_yz(data, truncate_atmark(filename), kind, region, descend_)

      CASE DEFAULT
         CALL assert(.FALSE., "invalid at-mark usage for '"//trim(filename)//"' in READ_DATA_YZ")

      END SELECT

    END SUBROUTINE read_atmark_data

  END SUBROUTINE read_data_yz_r8

!-----------------------------------------------------------------------------------------------------------------------

  RECURSIVE SUBROUTINE read_data_x_r8(data, filename, kind, region)
    REAL(8),      INTENT(OUT) :: data(isize)
    CHARACTER(*), INTENT(IN)  :: filename
    INTEGER,      INTENT(IN), OPTIONAL :: kind
    INTEGER,      INTENT(IN), OPTIONAL :: region

    TYPE(subregion_struct) :: rgn
    INTEGER :: istr, iend

    REAL(8), ALLOCATABLE :: tmp(:)

#ifdef PARALLEL_MPI
    REAL(8) :: buf(isize, 0:npes-1)
    INTEGER :: req(0:npes-1)
    INTEGER :: ierr
#endif

    REAL(8) :: d
    INTEGER :: i, n

    IF (present(region)) THEN
       rgn = regions(region)

       IF (region /= 0) THEN
          data(:) = UNDEF
       END IF
    ELSE
       rgn = regions(0)
    END IF

    istr = max(1, rgn%x_start-icoord*isize)
    iend = min(isize, rgn%x_end-icoord*isize)

    IF (check_dollar(filename)) THEN
       IF (.NOT. in_region_x(rgn)) RETURN

       d = read_literal(trim(filename(index(filename,'$',back=.TRUE.)+1:)))
       data(istr:iend) = d
       RETURN
    END IF

    ALLOCATE(tmp(rgn%x_start:rgn%x_end))

#ifdef PARALLEL_MPI
    IF (rank==0) THEN
       CALL read_file(tmp, filename, kind)

       req(:) = MPI_REQUEST_NULL

       DO n=0, npes-1
          IF (.NOT. in_region_x(rgn, n)) CYCLE

          istr = max(1, rgn%x_start-icoords(n)*isize)
          iend = min(isize, rgn%x_end-icoords(n)*isize)

          buf(:,n) = UNDEF

          DO i=istr, iend
             buf(i,n) = tmp(isize*icoords(n)+i)
          END DO
          CALL mpi_isend(buf(:,n), isize, MPI_REAL8, n, 0, comm, req(n), ierr)
       END DO
    END IF

    IF (in_region_x(rgn)) CALL mpi_recv(data, isize, MPI_REAL8, 0, 0, comm, MPI_STATUS_IGNORE, ierr)
    IF (rank==0) CALL mpi_waitall(npes, req, MPI_STATUSES_IGNORE, ierr)
#else
    CALL read_file(tmp, filename, kind)

    DO i=rgn%x_start, rgn%x_end
       data(i) = tmp(i)
    END DO
#endif
    DEALLOCATE(tmp)

    IF (rank==0) WRITE(REPORT_UNIT, *) "read '"//trim(filename)//"'"

  END SUBROUTINE read_data_x_r8

!-----------------------------------------------------------------------------------------------------------------------

  RECURSIVE SUBROUTINE read_data_y_r8(data, filename, kind, region)
    REAL(8),      INTENT(OUT) :: data(jsize)
    CHARACTER(*), INTENT(IN)  :: filename
    INTEGER,      INTENT(IN), OPTIONAL :: kind
    INTEGER,      INTENT(IN), OPTIONAL :: region

    TYPE(subregion_struct) :: rgn
    INTEGER :: jstr, jend

    REAL(8), ALLOCATABLE :: tmp(:)

#ifdef PARALLEL_MPI
    REAL(8) :: buf(jsize, 0:npes-1)
    INTEGER :: req(0:npes-1)
    INTEGER :: ierr
#endif

    REAL(8) :: d
    INTEGER :: j, n

    IF (present(region)) THEN
       rgn = regions(region)

       IF (region /= 0) THEN
          data(:) = UNDEF
       END IF
    ELSE
       rgn = regions(0)
    END IF

    jstr = max(1, rgn%y_start-jcoord*jsize)
    jend = min(jsize, rgn%y_end-jcoord*jsize)

    IF (check_dollar(filename)) THEN
       IF (.NOT. in_region_y(rgn)) RETURN

       d = read_literal(trim(filename(index(filename,'$',back=.TRUE.)+1:)))
       data(jstr:jend) = d
       RETURN
    END IF

    ALLOCATE(tmp(rgn%y_start:rgn%y_end))

#ifdef PARALLEL_MPI
    IF (rank==0) THEN
       CALL read_file(tmp, filename, kind)

       req(:) = MPI_REQUEST_NULL

       DO n=0, npes-1
          IF (.NOT. in_region_y(rgn, n)) CYCLE

          jstr = max(1, rgn%y_start-jcoords(n)*jsize)
          jend = min(jsize, rgn%y_end-jcoords(n)*jsize)

          buf(:,n) = UNDEF

          DO j=jstr, jend
             buf(j,n) = tmp(jsize*jcoords(n)+j)
          END DO
          CALL mpi_isend(buf(:,n), jsize, MPI_REAL8, n, 0, comm, req(n), ierr)
       END DO
    END IF

    IF (in_region_y(rgn)) CALL mpi_recv(data, jsize, MPI_REAL8, 0, 0, comm, MPI_STATUS_IGNORE, ierr)
    IF (rank==0) CALL mpi_waitall(npes, req, MPI_STATUSES_IGNORE, ierr)
#else
    CALL read_file(tmp, filename, kind)

    DO j=rgn%y_start, rgn%y_end
       data(j) = tmp(j)
    END DO
#endif
    DEALLOCATE(tmp)

    IF (rank==0) WRITE(REPORT_UNIT, *) "read '"//trim(filename)//"'"

  END SUBROUTINE read_data_y_r8

!-----------------------------------------------------------------------------------------------------------------------

  RECURSIVE SUBROUTINE read_data_z_r8(data, filename, kind, region, descend)
    REAL(8),      INTENT(OUT) :: data(ksize)
    CHARACTER(*), INTENT(IN)  :: filename
    INTEGER,      INTENT(IN), OPTIONAL :: kind
    INTEGER,      INTENT(IN), OPTIONAL :: region
    LOGICAL,      INTENT(IN), OPTIONAL :: descend

    TYPE(subregion_struct) :: rgn
    INTEGER :: kstr, kend

    REAL(8), ALLOCATABLE :: tmp(:)

#ifdef PARALLEL_MPI
    REAL(8) :: buf(ksize, 0:npes-1)
    INTEGER :: req(0:npes-1)
    INTEGER :: ierr
#endif

    LOGICAL :: descend_

    REAL(8) :: d
    INTEGER :: k, n

    IF (present(region)) THEN
       rgn = regions(region)

       IF (region /= 0) THEN
          data(:) = UNDEF
       END IF
    ELSE
       rgn = regions(0)
    END IF

    descend_ = .FALSE.
    IF (present(descend)) descend_ = descend

    kstr = max(1, rgn%z_start-kcoord*ksize)
    kend = min(ksize, rgn%z_end-kcoord*ksize)

    IF (check_dollar(filename)) THEN
       IF (.NOT. in_region_z(rgn)) RETURN

       d = read_literal(trim(filename(index(filename,'$',back=.TRUE.)+1:)))
       data(kstr:kend) = d
       RETURN
    END IF

    ALLOCATE(tmp(rgn%z_start:rgn%z_end))
    IF (rank==0) THEN
       CALL read_file(tmp, filename, kind)
       IF (descend_) CALL kreverse(tmp)
    END IF

#ifdef PARALLEL_MPI
    IF (rank==0) THEN
       req(:) = MPI_REQUEST_NULL

       DO n=0, npes-1
          IF (.NOT. in_region_z(rgn, n)) CYCLE

          kstr = max(1, rgn%z_start-kcoords(n)*ksize)
          kend = min(ksize, rgn%z_end-kcoords(n)*ksize)

          buf(:,n) = UNDEF

          DO k=kstr, kend
             buf(k,n) = tmp(ksize*kcoords(n)+k)
          END DO
          CALL mpi_isend(buf(:,n), ksize, MPI_REAL8, n, 0, comm, req(n), ierr)
       END DO
    END IF

    IF (in_region_z(rgn)) CALL mpi_recv(data, ksize, MPI_REAL8, 0, 0, comm, MPI_STATUS_IGNORE, ierr)
    IF (rank==0) CALL mpi_waitall(npes, req, MPI_STATUSES_IGNORE, ierr)
#else
    DO k=rgn%z_start, rgn%z_end
       data(k) = tmp(k)
    END DO
#endif
    DEALLOCATE(tmp)

    IF (rank==0) WRITE(REPORT_UNIT, *) "read '"//trim(filename)//"'"

  END SUBROUTINE read_data_z_r8

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE write_data_2d_r4(data, filename, kind, region)
    REAL(4),      INTENT(IN) :: data(isize, jsize)
    CHARACTER(*), INTENT(IN) :: filename
    INTEGER,      INTENT(IN), OPTIONAL :: kind
    INTEGER,      INTENT(IN), OPTIONAL :: region

    REAL(8) :: tmp(isize, jsize)

    tmp(:,:) = REAL(data, KIND=8)

    CALL write_data_2d_r8(tmp, filename, kind, region)

  END SUBROUTINE write_data_2d_r4

  SUBROUTINE write_data_3d_r4(data, filename, kind, region)
    REAL(4),      INTENT(IN) :: data(isize, jsize, ksize)
    CHARACTER(*), INTENT(IN):: filename
    INTEGER,      INTENT(IN), OPTIONAL :: kind
    INTEGER,      INTENT(IN), OPTIONAL :: region

    REAL(8) :: tmp(isize, jsize, ksize)

    tmp(:,:,:) = REAL(data, KIND=8)

    CALL write_data_3d_r8(tmp, filename, kind, region)

  END SUBROUTINE write_data_3d_r4

  SUBROUTINE write_data_xz_r4(data, filename, kind, region)
    REAL(4),      INTENT(IN) :: data(isize, ksize)
    CHARACTER(*), INTENT(IN) :: filename
    INTEGER,      INTENT(IN), OPTIONAL :: kind
    INTEGER,      INTENT(IN), OPTIONAL :: region

    REAL(8) :: tmp(isize, ksize)

    tmp(:,:) = REAL(data, KIND=8)

    CALL write_data_xz_r8(tmp, filename, kind, region)

  END SUBROUTINE write_data_xz_r4

  SUBROUTINE write_data_yz_r4(data, filename, kind, region)
    REAL(4),      INTENT(IN) :: data(jsize, ksize)
    CHARACTER(*), INTENT(IN) :: filename
    INTEGER,      INTENT(IN), OPTIONAL :: kind
    INTEGER,      INTENT(IN), OPTIONAL :: region

    REAL(8) :: tmp(jsize, ksize)

    tmp(:,:) = REAL(data, KIND=8)

    CALL write_data_yz_r8(tmp, filename, kind, region)

  END SUBROUTINE write_data_yz_r4

  SUBROUTINE write_data_x_r4(data, filename, kind, region)
    REAL(4),      INTENT(IN) :: data(isize)
    CHARACTER(*), INTENT(IN) :: filename
    INTEGER,      INTENT(IN), OPTIONAL :: kind
    INTEGER,      INTENT(IN), OPTIONAL :: region

    REAL(8) :: tmp(isize)

    tmp(:) = REAL(data, KIND=8)

    CALL write_data_x_r8(tmp, filename, kind, region)

  END SUBROUTINE write_data_x_r4

  SUBROUTINE write_data_y_r4(data, filename, kind, region)
    REAL(4),      INTENT(IN) :: data(jsize)
    CHARACTER(*), INTENT(IN) :: filename
    INTEGER,      INTENT(IN), OPTIONAL :: kind
    INTEGER,      INTENT(IN), OPTIONAL :: region

    REAL(8) :: tmp(jsize)

    tmp(:) = REAL(data, KIND=8)

    CALL write_data_x_r8(tmp, filename, kind, region)

  END SUBROUTINE write_data_y_r4

  SUBROUTINE write_data_z_r4(data, filename, kind, region)
    REAL(4),      INTENT(IN) :: data(ksize)
    CHARACTER(*), INTENT(IN) :: filename
    INTEGER,      INTENT(IN), OPTIONAL :: kind
    INTEGER,      INTENT(IN), OPTIONAL :: region

    REAL(8) :: tmp(ksize)

    tmp(:) = REAL(data, KIND=8)

    CALL write_data_z_r8(tmp, filename, kind, region)

  END SUBROUTINE write_data_z_r4

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE write_data_2d_r8(data, filename, kind, region)
    REAL(8),      INTENT(IN) :: data(isize, jsize)
    CHARACTER(*), INTENT(IN) :: filename
    INTEGER,      INTENT(IN), OPTIONAL :: kind
    INTEGER,      INTENT(IN), OPTIONAL :: region

    INTEGER :: kind_, region_

    INTEGER :: i, j

PROFILE_BEGIN('iowrite')

    kind_ = 8
    IF (present(kind)) kind_ = kind
    CALL assert(kind_==1 .OR. kind_==4 .OR. kind_==8, "unsupported KIND in WRITE_DATA_2D")

    region_ = 0
    IF (present(region)) region_ = region

#ifdef PARALLEL_MPI
    IF (region_ /= 0) THEN
       CALL assert(regions(region)%defined, "undefined subregion-code in WRITE_DATA_2D")

       ! sub-region output is supported only at the layered routine
       CALL write_data_2d_layered(data, filename, kind_, region_)
       IF (rank==0) WRITE(REPORT_UNIT, *) "write '"//trim(filename)//"'"
       RETURN
    END IF

    SELECT CASE(trim(output_method))
#ifdef MPIIO
    CASE ('MPI-IO', 'MPIIO', 'mpi-io', 'mpiio')
       CALL write_data_2d_mpiio(data, filename, kind_)
#endif
    CASE ('LAYERED', 'layered')
       CALL write_data_2d_layered(data, filename, kind_)
    CASE ('SERIAL', 'serial')
       CALL write_data_2d_serial(data, filename, kind_)
    CASE DEFAULT
       CALL assert(.FALSE., "parallel output method '"//trim(output_method)//"' is not supported")
    END SELECT
#else
    CALL write_file(data(regions(region_)%x_start:regions(region_)%x_end,  &
                         regions(region_)%y_start:regions(region_)%y_end), filename, kind_)
#endif

PROFILE_END('iowrite')

    IF (rank==0) WRITE(REPORT_UNIT, *) "write '"//trim(filename)//"'"
  END SUBROUTINE write_data_2d_r8

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE write_data_3d_r8(data, filename, kind, region)
    REAL(8),      INTENT(IN) :: data(isize, jsize, ksize)
    CHARACTER(*), INTENT(IN) :: filename
    INTEGER,      INTENT(IN), OPTIONAL :: kind
    INTEGER,      INTENT(IN), OPTIONAL :: region

    INTEGER :: kind_, region_

    INTEGER :: i, j, k

PROFILE_BEGIN('iowrite')

    kind_ = 8
    IF (present(kind)) kind_ = kind
    CALL assert(kind_==1 .OR. kind_==4 .OR. kind_==8, "unsupported KIND in WRITE_DATA_3D")

    region_ = 0
    IF (present(region)) region_ = region

#ifdef PARALLEL_MPI
    IF (region_ /= 0) THEN
       CALL assert(regions(region)%defined, "undefined subregion-code in WRITE_DATA_3D")
       ! sub-region output is supported only at the layered routine
       CALL write_data_3d_layered(data, filename, kind_, region_)
       IF (rank==0) WRITE(REPORT_UNIT, *) "write '"//trim(filename)//"'"
       RETURN
    END IF

    SELECT CASE(trim(output_method))
#ifdef MPIIO
    CASE ('MPI-IO', 'MPIIO', 'mpi-io', 'mpiio')
       CALL write_data_3d_mpiio(data, filename, kind_)
#endif
    CASE ('LAYERED', 'layered')
       CALL write_data_3d_layered(data, filename, kind_)
    CASE ('SERIAL', 'serial')
       CALL write_data_3d_serial(data, filename, kind_)
    CASE DEFAULT
       CALL assert(.FALSE., "parallel output method '"//trim(output_method)//"' is not supported")
    END SELECT
#else
    CALL write_file(data(regions(region_)%x_start:regions(region_)%x_end,  &
                         regions(region_)%y_start:regions(region_)%y_end,  &
                         regions(region_)%z_start:regions(region_)%z_end), filename, kind_)
#endif

PROFILE_END('iowrite')

    IF (rank==0) WRITE(REPORT_UNIT, *) "write '"//trim(filename)//"'"
  END SUBROUTINE write_data_3d_r8

!-----------------------------------------------------------------------------------------------------------------------

#ifdef PARALLEL_MPI
#ifdef MPIIO

  SUBROUTINE write_data_2d_mpiio(data, filename, kind)
    REAL(8),      INTENT(IN) :: data(isize, jsize)
    CHARACTER(*), INTENT(IN) :: filename
    INTEGER,      INTENT(IN) :: kind

    INTEGER :: subarray, ierr

    IF (vrank /= 0) RETURN

    IF (hrank==0 .AND. mpiio_zerofill) CALL dump_file(0.0, filename, dimx*dimy, kind=kind)

    CALL assert(kind==1 .OR. kind==4 .OR. kind==8, "unsupported KIND in WRITE_DATA_2D_MPIIO")

    CALL mpi_file_open(hcomm, trim(filename), MPI_MODE_CREATE + MPI_MODE_WRONLY, MPI_INFO_NULL, file%handle, ierr)
    CALL assert(ierr==MPI_SUCCESS, "failed to open '"//trim(filename)//"'")

    SELECT CASE (kind)
    CASE (1)
       file%buffer(1:isize*jsize*kind) = transfer(INT( data, KIND=1), file%buffer)
       subarray = views(0)%subarray_2d_i1
    CASE (4)
       file%buffer(1:isize*jsize*kind) = transfer(REAL(data, KIND=4), file%buffer)
       subarray = views(0)%subarray_2d_r4
    CASE (8)
       file%buffer(1:isize*jsize*kind) = transfer(REAL(data, KIND=8), file%buffer)
       subarray = views(0)%subarray_2d_r8
    END SELECT

    IF (.NOT. check_endian()) CALL convert_endian(file%buffer, isize*jsize, kind)

    CALL mpi_file_set_view(file%handle, views(0)%offset, MPI_BYTE, subarray, 'native', MPI_INFO_NULL, ierr)

    CALL mpi_file_write_all(file%handle, file%buffer, isize*jsize*kind, MPI_BYTE, MPI_STATUS_IGNORE, ierr)

    CALL mpi_file_close(file%handle, ierr)

  END SUBROUTINE write_data_2d_mpiio

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE write_data_3d_mpiio(data, filename, kind)
    REAL(8),      INTENT(IN) :: data(isize, jsize, ksize)
    CHARACTER(*), INTENT(IN) :: filename
    INTEGER,      INTENT(IN) :: kind

    INTEGER :: n
    INTEGER :: ierr

    CALL assert(kind==1 .OR. kind==4 .OR. kind==8, "unsupported KIND in WRITE_DATA_3D_MPIIO")
#ifdef MPIIO_ASYNCHRONOUS
    CALL assert(.NOT. mpiio_zerofill, "MPIIO_ZEROFILL is not supported for ASYNCHRONOUS OUTPUT")

    n = query_async_write()

    CALL c_(file_async(n))
    CALL mpi_file_write_all_begin(file_async(n)%handle, file_async(n)%buffer, isize*jsize*ksize*kind, MPI_BYTE, ierr)
#else
    IF (rank==0 .AND. mpiio_zerofill) CALL dump_file(0.0, filename, dimx*dimy*dimz, stride=dimx*dimy, kind=kind)

    CALL c_(file)

    CALL mpi_file_write_all(file%handle, file%buffer, isize*jsize*ksize*kind, MPI_BYTE, MPI_STATUS_IGNORE, ierr)
    CALL mpi_file_close(file%handle, ierr)
#endif

  CONTAINS
    SUBROUTINE c_(f)
      TYPE(file_struct), INTENT(INOUT) :: f
      INTEGER :: subarray, ierr

      CALL mpi_file_open(comm, trim(filename), MPI_MODE_CREATE + MPI_MODE_WRONLY, MPI_INFO_NULL, f%handle, ierr)
      CALL assert(ierr==MPI_SUCCESS, "failed to open '"//trim(filename)//"'")

      SELECT CASE (kind)
      CASE (1)
         f%buffer(1:isize*jsize*ksize*kind) = transfer(INT( data, KIND=1), f%buffer)
         subarray = views(0)%subarray_3d_i1
      CASE (4)
         f%buffer(1:isize*jsize*ksize*kind) = transfer(REAL(data, KIND=4), f%buffer)
         subarray = views(0)%subarray_3d_r4
      CASE (8)
         f%buffer(1:isize*jsize*ksize*kind) = transfer(REAL(data, KIND=8), f%buffer)
         subarray = views(0)%subarray_3d_r8
      END SELECT
      
      IF (.NOT. check_endian()) CALL convert_endian(f%buffer, isize*jsize*ksize, kind)
      
      CALL mpi_file_set_view(f%handle, views(0)%offset, MPI_BYTE, subarray, 'native', MPI_INFO_NULL, ierr)

    END SUBROUTINE c_

  END SUBROUTINE write_data_3d_mpiio

#endif
!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE write_data_2d_serial(data, filename, kind)
    REAL(8),      INTENT(IN) :: data(isize, jsize)
    CHARACTER(*), INTENT(IN) :: filename
    INTEGER,      INTENT(IN) :: kind

    REAL(8) :: buf(isize, jsize)

    INTEGER :: n, i, j
    INTEGER :: iostat, ierr, rec

    IF (vrank /= 0) RETURN

    CALL assert(kind==1 .OR. kind==4 .OR. kind==8, "unsupported KIND in WRITE_DATA_2D_SERIAL")

    IF (hrank==0) THEN
       OPEN(UNIT   = TMP_UNIT,      &
            FILE   = trim(filename),&
            FORM   = 'UNFORMATTED', &
            ACCESS = 'DIRECT',      &
            STATUS = 'REPLACE',     &
            ACTION = 'WRITE',       &
            RECL   = isize*kind,    &
            IOSTAT = iostat)

       CALL assert(iostat==0, "failed to open '"//trim(filename)//"'")

       DO n=0, hpes-1
          IF (n > 0) THEN
             CALL mpi_recv(buf, isize*jsize, MPI_REAL8, n, 0, hcomm, MPI_STATUS_IGNORE, ierr)
          ELSE
             buf(:,:) = data(:,:)
          END IF

          DO j=1, jsize
             rec = (jcoords_h(n)*jsize+j-1)*ipes + icoords_h(n) + 1
             SELECT CASE (kind)
             CASE (1)
                WRITE(TMP_UNIT, REC=rec) INT( buf(:,j), KIND=1)
             CASE (4)
                WRITE(TMP_UNIT, REC=rec) REAL(buf(:,j), KIND=4)
             CASE (8)
                WRITE(TMP_UNIT, REC=rec) REAL(buf(:,j), KIND=8)
             END SELECT
          END DO
       END DO
       CLOSE(TMP_UNIT)
    ELSE
       CALL mpi_send(data, isize*jsize, MPI_REAL8, 0, 0, hcomm, ierr)
    END IF

    CALL mpi_barrier(hcomm, ierr)

  END SUBROUTINE write_data_2d_serial

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE write_data_3d_serial(data, filename, kind)
    REAL(8),      INTENT(IN) :: data(isize, jsize, ksize)
    CHARACTER(*), INTENT(IN) :: filename
    INTEGER,      INTENT(IN) :: kind

    REAL(8) :: buf(isize, jsize, ksize)

    INTEGER :: n, i, j, k
    INTEGER :: iostat, ierr, rec

    CALL assert(kind==1 .OR. kind==4 .OR. kind==8, "unsupported KIND in WRITE_DATA_3D_SERIAL")

    IF (rank==0) THEN
       OPEN(UNIT   = TMP_UNIT,      &
            FILE   = trim(filename),&
            FORM   = 'UNFORMATTED', &
            ACCESS = 'DIRECT',      &
            STATUS = 'REPLACE',     &
            ACTION = 'WRITE',       &
            RECL   = isize*kind,    &
            IOSTAT = iostat)

       CALL assert(iostat==0, "failed to open '"//trim(filename)//"'")

       DO n=0, npes-1
          IF (n > 0) THEN
             CALL mpi_recv(buf, isize*jsize*ksize, MPI_REAL8, n, 0, comm, MPI_STATUS_IGNORE, ierr)
          ELSE
             buf = data
          END IF

          DO k=1, ksize
          DO j=1, jsize
             rec = (kcoords(n)*ksize+k-1)*dimy*ipes + (jcoords(n)*jsize+j-1)*ipes + icoords(n) + 1

             SELECT CASE (kind)
             CASE (1)
                WRITE(TMP_UNIT, REC=rec) INT( buf(:,j,k), KIND=1)
             CASE (4)
                WRITE(TMP_UNIT, REC=rec) REAL(buf(:,j,k), KIND=4)
             CASE (8)
                WRITE(TMP_UNIT, REC=rec) REAL(buf(:,j,k), KIND=8)
             END SELECT
          END DO
          END DO
       END DO

       CLOSE(TMP_UNIT)
    ELSE
       CALL mpi_send(data, isize*jsize*ksize, MPI_REAL8, 0, 0, comm, ierr)
    END IF
    CALL mpi_barrier(comm, ierr)

  END SUBROUTINE write_data_3d_serial

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE write_data_2d_layered(data, filename, kind, region)
    REAL(8),      INTENT(IN) :: data(isize, jsize)
    CHARACTER(*), INTENT(IN) :: filename
    INTEGER,      INTENT(IN) :: kind
    INTEGER,      INTENT(IN), OPTIONAL :: region

    REAL(8) :: buf(isize, jsize, 0:hpes-1)
    REAL(8) :: tmp(isize*ipes, jsize*jpes)

    INTEGER :: n, i, j
    INTEGER :: iostat, ierr, rec

    INTEGER :: req(1:hpes-1)

    TYPE(subregion_struct) :: rgn

    IF (present(region)) THEN
       CALL assert(regions(region)%defined, "undefined subregion-code in WRITE_DATA_2D_LAYERED")
       rgn = regions(region)
    ELSE
       rgn = regions(0)
    END IF

    IF (vrank /= 0) RETURN

    CALL assert(kind==1 .OR. kind==4 .OR. kind==8, "unsupported KIND in WRITE_DATA_2D_LAYERED")

    IF (hrank==0) THEN
       OPEN(UNIT   = TMP_UNIT,                   &
            FILE   = trim(filename),             &
            FORM   = 'UNFORMATTED',              &
            ACCESS = 'DIRECT',                   &
            STATUS = 'REPLACE',                  &
            ACTION = 'WRITE',                    &
            RECL   = rgn%x_size*rgn%y_size*kind, &
            IOSTAT = iostat)

       CALL assert(iostat==0, "failed to open '"//trim(filename)//"'")

       DO n=1, hpes-1
          IF (rgn%x_start <= (icoords_h(n)+1)*isize .AND. rgn%x_end >= icoords_h(n)*isize+1 .AND. &
              rgn%y_start <= (jcoords_h(n)+1)*jsize .AND. rgn%y_end >= jcoords_h(n)*jsize+1) THEN
             CALL mpi_irecv(buf(:,:,n), isize*jsize, MPI_REAL8, n, 0, hcomm, req(n), ierr)
          ELSE
             req(n) = MPI_REQUEST_NULL
             buf(:,:,n) = 0.0D0
          END IF
       END DO

       tmp(:,:) = 0.0

       tmp(icoords_h(0)*isize+1:(icoords_h(0)+1)*isize, &
           jcoords_h(0)*jsize+1:(jcoords_h(0)+1)*jsize) = data(:,:)

       DO n=1, hpes-1
          CALL mpi_wait(req(n), MPI_STATUS_IGNORE, ierr)
          tmp(icoords_h(n)*isize+1:(icoords_h(n)+1)*isize, &
              jcoords_h(n)*jsize+1:(jcoords_h(n)+1)*jsize) = buf(:,:,n)
       END DO

       SELECT CASE (kind)
       CASE (1)
          WRITE(TMP_UNIT, REC=1) INT( tmp(rgn%x_start:rgn%x_end, rgn%y_start:rgn%y_end), KIND=1)
       CASE (4)
          WRITE(TMP_UNIT, REC=1) REAL(tmp(rgn%x_start:rgn%x_end, rgn%y_start:rgn%y_end), KIND=4)
       CASE (8)
          WRITE(TMP_UNIT, REC=1) REAL(tmp(rgn%x_start:rgn%x_end, rgn%y_start:rgn%y_end), KIND=8)
       END SELECT

       CLOSE(TMP_UNIT)
    ELSE
       IF (rgn%x_start <= (icoord+1)*isize .AND. rgn%x_end >= icoord*isize+1 .AND. &
           rgn%y_start <= (jcoord+1)*jsize .AND. rgn%y_end >= jcoord*jsize+1) THEN
          CALL mpi_send(data, isize*jsize, MPI_REAL8, 0, 0, hcomm, ierr)
       END IF
    END IF

  END SUBROUTINE write_data_2d_layered

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE write_data_3d_layered(data, filename, kind, region)
    REAL(8),      INTENT(IN) :: data(isize, jsize, ksize)
    CHARACTER(*), INTENT(IN) :: filename
    INTEGER,      INTENT(IN) :: kind
    INTEGER,      INTENT(IN), OPTIONAL :: region

    REAL(8) :: buf(isize, jsize, 0:npes-1)
    REAL(8) :: tmp(isize*ipes, jsize*jpes)

    INTEGER :: n, i, j, k
    INTEGER :: iostat, ierr, rec

    INTEGER :: req(1:npes-1)

    TYPE(subregion_struct) :: rgn

    IF (present(region)) THEN
       CALL assert(regions(region)%defined, "undefined subregion-code in WRITE_DATA_3D_LAYERED")
       rgn = regions(region)
    ELSE
       rgn = regions(0)
    END IF

    CALL assert(kind==1 .OR. kind==4 .OR. kind==8, "unsupported KIND in WRITE_DATA_3D_LAYERED")

    IF (rank==0) THEN
       OPEN(UNIT   = TMP_UNIT,                   &
            FILE   = trim(filename),             &
            FORM   = 'UNFORMATTED',              &
            ACCESS = 'DIRECT',                   &
            STATUS = 'REPLACE',                  &
            ACTION = 'WRITE',                    &
            RECL   = rgn%x_size*rgn%y_size*kind, &
            IOSTAT = iostat)

       CALL assert(iostat==0, "failed to open '"//trim(filename)//"'")

       DO k=rgn%z_start, rgn%z_end
          DO n=1, npes-1
             IF (rgn%x_start <= (icoords(n)+1)*isize .AND. rgn%x_end >= icoords(n)*isize+1 .AND. &
                 rgn%y_start <= (jcoords(n)+1)*jsize .AND. rgn%y_end >= jcoords(n)*jsize+1 .AND. &
                 k <= (kcoords(n)+1)*ksize .AND. k >= kcoords(n)*ksize+1) THEN
                CALL mpi_irecv(buf(:,:,n), isize*jsize, MPI_REAL8, n, k, comm, req(n), ierr)
             ELSE
                req(n) = MPI_REQUEST_NULL
             END IF
          END DO

          tmp(:,:) = 0.0

          IF (rgn%x_start <= (icoord+1)*isize .AND. rgn%x_end >= icoord*isize+1 .AND. &
              rgn%y_start <= (jcoord+1)*jsize .AND. rgn%y_end >= jcoord*jsize+1 .AND. &
              k <= (kcoord+1)*ksize .AND. k >= kcoord*ksize+1) THEN
             tmp(icoord*isize+1:(icoord+1)*isize, jcoord*jsize+1:(jcoord+1)*jsize) = data(:,:,k-kcoord*ksize)
          END IF

          DO n=1, npes-1
             IF (req(n) == MPI_REQUEST_NULL) CYCLE
             CALL mpi_wait(req(n), MPI_STATUS_IGNORE, ierr)
             tmp(icoords(n)*isize+1:(icoords(n)+1)*isize, &
                 jcoords(n)*jsize+1:(jcoords(n)+1)*jsize) = buf(:,:,n)
          END DO

          SELECT CASE (kind)
          CASE (1)
             WRITE(TMP_UNIT, REC=k-rgn%z_start+1) INT( tmp(rgn%x_start:rgn%x_end, rgn%y_start:rgn%y_end), KIND=1)
          CASE (4)
             WRITE(TMP_UNIT, REC=k-rgn%z_start+1) REAL(tmp(rgn%x_start:rgn%x_end, rgn%y_start:rgn%y_end), KIND=4)
          CASE (8)
             WRITE(TMP_UNIT, REC=k-rgn%z_start+1) REAL(tmp(rgn%x_start:rgn%x_end, rgn%y_start:rgn%y_end), KIND=8)
          END SELECT
       END DO
       CLOSE(TMP_UNIT)
    ELSE
       IF (rgn%x_start <= (icoord+1)*isize .AND. rgn%x_end >= icoord*isize+1 .AND. &
           rgn%y_start <= (jcoord+1)*jsize .AND. rgn%y_end >= jcoord*jsize+1) THEN
          DO k=rgn%z_start, rgn%z_end
             IF (k <= (kcoord+1)*ksize .AND. k >= kcoord*ksize+1) THEN
                CALL mpi_send(data(:,:,k-kcoord*ksize), isize*jsize, MPI_REAL8, 0, k, comm, ierr)
             END IF
          END DO
       END IF
    END IF

  END SUBROUTINE write_data_3d_layered

!-----------------------------------------------------------------------------------------------------------------------

#ifdef MPIIO_ASYNCHRONOUS
  INTEGER FUNCTION query_async_write()
    INTEGER, SAVE :: n = 0
    INTEGER :: ierr

    n = mod(n,n_async)+1

    IF (file_async(n)%handle /= MPI_FILE_NULL) THEN
       CALL mpi_file_write_all_end(file_async(n)%handle, file_async(n)%buffer, MPI_STATUS_IGNORE, ierr)
       CALL mpi_file_close(file_async(n)%handle, ierr)
       CALL assert(file_async(n)%handle == MPI_FILE_NULL, "error in closing file")
    END IF

    query_async_write = n

  END FUNCTION query_async_write
#endif
#endif

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE write_data_xz_r8(data, filename, kind, region)
    REAL(8),      INTENT(IN) :: data(isize, ksize)
    CHARACTER(*), INTENT(IN) :: filename
    INTEGER,      INTENT(IN), OPTIONAL :: kind
    INTEGER,      INTENT(IN), OPTIONAL :: region

#ifdef PARALLEL_MPI
    REAL(8) :: buf1(dimx, dimz)
    REAL(8) :: buf2(isize, ksize, 0:npes-1)
    LOGICAL :: flag(0:ipes-1, 0:kpes-1)

    INTEGER :: n, i, k
    INTEGER :: ierr
#endif

    INTEGER :: x_start, x_end
    INTEGER :: z_start, z_end

    IF (PRESENT(region)) THEN
       CALL assert(regions(region)%defined, "undefined subregion-code in WRITE_DATA_XZ")
       x_start = regions(region)%x_start
       z_start = regions(region)%z_start
       x_end   = regions(region)%x_end
       z_end   = regions(region)%z_end
    ELSE
       x_start = 1
       z_start = 1
       x_end   = dimx
       z_end   = dimz
    END IF

#ifdef PARALLEL_MPI
    CALL mpi_gather(data, isize*ksize, MPI_REAL8, &
                    buf2, isize*ksize, MPI_REAL8, 0, comm, ierr)

    IF (rank /= 0) RETURN

    flag(:,:) = .FALSE.
    buf1(:,:) = 0.0

    DO n=0, npes-1
       DO k=1, ksize
       DO i=1, isize
          buf1(isize*icoords(n)+i,ksize*kcoords(n)+k) = buf1(isize*icoords(n)+i,ksize*kcoords(n)+k) + buf2(i,k,n)
       END DO
       END DO
       flag(icoords(n),kcoords(n)) = .TRUE.
    END DO

    DO k=0, kpes-1
    DO i=0, ipes-1
       IF (.NOT. flag(i,k)) buf1(isize*i+1:isize*(i+1),ksize*k+1:ksize*(k+1)) = 0.0
    END DO
    END DO

    CALL write_file(buf1(x_start:x_end, z_start:z_end), filename, kind)
#else
    CALL write_file(data(x_start:x_end, z_start:z_end), filename, kind)
#endif

    IF (rank==0) WRITE(REPORT_UNIT, *) "write '"//trim(filename)//"'"
  END SUBROUTINE write_data_xz_r8

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE write_data_yz_r8(data, filename, kind, region)
    REAL(8),      INTENT(IN) :: data(jsize, ksize)
    CHARACTER(*), INTENT(IN) :: filename
    INTEGER,      INTENT(IN), OPTIONAL :: kind
    INTEGER,      INTENT(IN), OPTIONAL :: region

#ifdef PARALLEL_MPI
    REAL(8) :: buf1(dimy, dimz)
    REAL(8) :: buf2(jsize, ksize, 0:npes-1)
    LOGICAL :: flag(0:jpes-1, 0:kpes-1)

    INTEGER :: n, j, k
    INTEGER :: ierr
#endif
    INTEGER :: y_start, y_end
    INTEGER :: z_start, z_end

    IF (PRESENT(region)) THEN
       CALL assert(regions(region)%defined, "undefined subregion-code in WRITE_DATA_YZ")
       y_start = regions(region)%y_start
       z_start = regions(region)%z_start
       y_end   = regions(region)%y_end
       z_end   = regions(region)%z_end
    ELSE
       y_start = 1
       z_start = 1
       y_end   = dimy
       z_end   = dimz
    END IF

#ifdef PARALLEL_MPI
    CALL mpi_gather(data, jsize*ksize, MPI_REAL8, &
                    buf2, jsize*ksize, MPI_REAL8, 0, comm, ierr)

    IF (rank /= 0) RETURN

    flag(:,:) = .FALSE.
    buf1(:,:) = 0.0

    DO n=0, npes-1
       DO k=1, ksize
       DO j=1, jsize
          buf1(jsize*jcoords(n)+j,ksize*kcoords(n)+k) = buf1(jsize*jcoords(n)+j,ksize*kcoords(n)+k) + buf2(j,k,n)
       END DO
       END DO
       flag(jcoords(n),kcoords(n)) = .TRUE.
    END DO

    DO k=0, kpes-1
    DO j=0, jpes-1
       IF (flag(j,k)) CYCLE
       buf1(jsize*j+1:jsize*(j+1),ksize*k+1:ksize*(k+1)) = 0.0
    END DO
    END DO

    CALL write_file(buf1(y_start:y_end, z_start:z_end), filename, kind)
#else
    CALL write_file(data(y_start:y_end, z_start:z_end), filename, kind)
#endif

    IF (rank==0) WRITE(REPORT_UNIT, *) "write '"//trim(filename)//"'"
  END SUBROUTINE write_data_yz_r8

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE write_data_x_r8(data, filename, kind, region)
    REAL(8),      INTENT(IN) :: data(isize)
    CHARACTER(*), INTENT(IN) :: filename
    INTEGER,      INTENT(IN), OPTIONAL :: kind
    INTEGER,      INTENT(IN), OPTIONAL :: region

#ifdef PARALLEL_MPI
    REAL(8) :: buf1(dimx)
    REAL(8) :: buf2(isize, 0:npes-1)
    LOGICAL :: flag(0:ipes-1)

    INTEGER :: n, i
    INTEGER :: ierr
#endif

    INTEGER :: x_start, x_end

    IF (PRESENT(region)) THEN
       CALL assert(regions(region)%defined, "undefined subregion-code in WRITE_DATA_X")
       x_start = regions(region)%x_start
       x_end   = regions(region)%x_end
    ELSE
       x_start = 1
       x_end   = dimx
    END IF

#ifdef PARALLEL_MPI
    CALL mpi_gather(data, isize, MPI_REAL8, &
                    buf2, isize, MPI_REAL8, 0, comm, ierr)

    IF (rank /= 0) RETURN

    flag(:) = .FALSE.
    buf1(:) = 0.0

    DO n=0, npes-1
       DO i=1, isize
          buf1(isize*icoords(n)+i) = buf1(isize*icoords(n)+i) + buf2(i,n)
       END DO
       flag(icoords(n)) = .TRUE.
    END DO

    DO i=0, ipes-1
       IF (.NOT. flag(i)) buf1(isize*i+1:isize*(i+1)) = 0.0
    END DO

    CALL write_file(buf1(x_start:x_end), filename, kind)
#else
    CALL write_file(data(x_start:x_end), filename, kind)
#endif

    IF (rank==0) WRITE(REPORT_UNIT, *) "write '"//trim(filename)//"'"
  END SUBROUTINE write_data_x_r8

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE write_data_y_r8(data, filename, kind, region)
    REAL(8),      INTENT(IN) :: data(jsize)
    CHARACTER(*), INTENT(IN) :: filename
    INTEGER,      INTENT(IN), OPTIONAL :: kind
    INTEGER,      INTENT(IN), OPTIONAL :: region

#ifdef PARALLEL_MPI
    REAL(8) :: buf1(dimy)
    REAL(8) :: buf2(jsize, 0:npes-1)
    LOGICAL :: flag(0:jpes-1)

    INTEGER :: n, j
    INTEGER :: ierr
#endif

    INTEGER :: y_start, y_end

    IF (PRESENT(region)) THEN
       CALL assert(regions(region)%defined, "undefined subregion-code in WRITE_DATA_Y")
       y_start = regions(region)%y_start
       y_end   = regions(region)%y_end
    ELSE
       y_start = 1
       y_end   = dimy
    END IF

#ifdef PARALLEL_MPI
    CALL mpi_gather(data, jsize, MPI_REAL8, &
                    buf2, jsize, MPI_REAL8, 0, comm, ierr)

    IF (rank /= 0) RETURN

    flag(:) = .FALSE.
    buf1(:) = 0.0

    DO n=0, npes-1
       DO j=1, jsize
          buf1(jsize*jcoords(n)+j) = buf1(jsize*jcoords(n)+j) + buf2(j,n)
       END DO
       flag(jcoords(n)) = .TRUE.
    END DO

    DO j=0, jpes-1
       IF (.NOT. flag(j)) buf1(jsize*j+1:jsize*(j+1)) = 0.0
    END DO

    CALL write_file(buf1(y_start:y_end), filename, kind)
#else
    CALL write_file(data(y_start:y_end), filename, kind)
#endif

    IF (rank==0) WRITE(REPORT_UNIT, *) "write '"//trim(filename)//"'"
  END SUBROUTINE write_data_y_r8

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE write_data_z_r8(data, filename, kind, region)
    REAL(8),      INTENT(IN) :: data(ksize)
    CHARACTER(*), INTENT(IN) :: filename
    INTEGER,      INTENT(IN), OPTIONAL :: kind
    INTEGER,      INTENT(IN), OPTIONAL :: region

#ifdef PARALLEL_MPI
    REAL(8) :: buf1(dimz)
    REAL(8) :: buf2(ksize, 0:npes-1)
    LOGICAL :: flag(0:kpes-1)

    INTEGER :: n, k
    INTEGER :: ierr
#endif

    INTEGER :: z_start, z_end

    IF (PRESENT(region)) THEN
       CALL assert(regions(region)%defined, "undefined subregion-code in WRITE_DATA_Z")
       z_start = regions(region)%z_start
       z_end   = regions(region)%z_end
    ELSE
       z_start = 1
       z_end   = dimz
    END IF

#ifdef PARALLEL_MPI
    CALL mpi_gather(data, ksize, MPI_REAL8, &
                    buf2, ksize, MPI_REAL8, 0, comm, ierr)

    IF (rank /= 0) RETURN

    flag(:) = .FALSE.
    buf1(:) = 0.0

    DO n=0, npes-1
       DO k=1, ksize
          buf1(ksize*kcoords(n)+k) = buf1(ksize*kcoords(n)+k) + buf2(k,n)
       END DO
       flag(kcoords(n)) = .TRUE.
    END DO

    DO k=0, kpes-1
       IF (.NOT. flag(k)) buf1(ksize*k+1:ksize*(k+1)) = 0.0
    END DO

    CALL write_file(buf1(z_start:z_end), filename, kind)
#else
    CALL write_file(data(z_start:z_end), filename, kind)
#endif

    IF (rank==0) WRITE(REPORT_UNIT, *) "write '"//trim(filename)//"'"
  END SUBROUTINE write_data_z_r8

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE read_geometry(unit, metric, topo)
    INTEGER, OPTIONAL :: unit
    LOGICAL, OPTIONAL :: metric
    LOGICAL, OPTIONAL :: topo

    REAL(8) :: delta_x(maxdim)
    REAL(8) :: delta_y(maxdim)
    REAL(8) :: delta_z(maxdim)

    CHARACTER(512) :: inputdir
    CHARACTER(128) :: metric_x
    CHARACTER(128) :: metric_y
    CHARACTER(128) :: bathymetry
    CHARACTER(128) :: iceshelf

    LOGICAL :: read_metric
    LOGICAL :: read_topo

    INTEGER :: i, j, k

    INTEGER :: iostat

    NAMELIST / geometry /           &
         dimx, dimy, dimz,          &
         cycle_x, cycle_y, cycle_z, &
         delta_x, delta_y, delta_z, &
         metric_x,  metric_y,       &
         inputdir,                  &
         bathymetry, iceshelf

    dimx = 1
    dimy = 1
    dimz = 1

    cycle_x = .FALSE.
    cycle_y = .FALSE.
    cycle_z = .FALSE.

    inputdir   = ''
    metric_x   = ''
    metric_y   = ''
    bathymetry = ''
    iceshelf   = ''

    delta_x(:) = 0.0
    delta_y(:) = 0.0
    delta_z(:) = 0.0

    IF (present(unit)) THEN
       REWIND(unit)
       READ(unit, NML=geometry, IOSTAT=iostat)
    ELSE
       READ(*,    NML=geometry, IOSTAT=iostat)
    END IF

    CALL assert(iostat>=0, "failed to read GEOMETRY namelist")

    DO i=2, dimx
       IF (delta_x(i) == 0.0) delta_x(i) = delta_x(i-1)
    END DO

    DO j=2, dimy
       IF (delta_y(j) == 0.0) delta_y(j) = delta_y(j-1)
    END DO

    DO k=2, dimz
       IF (delta_z(k) == 0.0) delta_z(k) = delta_z(k-1)
    END DO

    ALLOCATE(dx(1-slv:dimx+slv, 1-slv:dimy+slv))
    ALLOCATE(dy(1-slv:dimx+slv, 1-slv:dimy+slv))
    ALLOCATE(dz(1-slv:dimz+slv))

    dx(:,:)   = 1.0
    dy(:,:)   = 1.0
    dz(:)     = 1.0

    read_metric = .TRUE.
    IF (present(metric)) read_metric = metric

    IF (read_metric .AND. trim(metric_x) /= '') THEN
       CALL read_file(dx(1:dimx, 1:dimy), path(inputdir, metric_x))
    END IF
    DO i=1, dimx
       dx(i,:) = dx(i,:) * delta_x(i)
    END DO

    IF (read_metric .AND. trim(metric_y) /= '') THEN
       CALL read_file(dy(1:dimx, 1:dimy), path(inputdir, metric_y))
    END IF

    DO j=1, dimy
       dy(:,j) = dy(:,j) * delta_y(j)
    END DO

    IF (cycle_x) THEN
       DO i=1, slv
          dx(dimx+i,:) = dx(i,:)
          dy(dimx+i,:) = dy(i,:)

          dx(i-slv,:)  = dx(dimx+i-slv,:)
          dy(i-slv,:)  = dy(dimx+i-slv,:)
       END DO
    ELSE
       DO i=1, slv
          dx(dimx+i,:) = dx(dimx,:)
          dy(dimx+i,:) = dy(dimx,:)

          dx(i-slv,:)  = dx(1,:)
          dy(i-slv,:)  = dy(1,:)
       END DO
    END IF

    IF (cycle_y) THEN
       DO j=1, slv
          dx(:,dimy+j) = dx(:,j)
          dy(:,dimy+j) = dy(:,j)

          dx(:,j-slv)  = dx(:,dimy+j-slv)
          dy(:,j-slv)  = dy(:,dimy+j-slv)
       END DO
    ELSE
       DO j=1, slv
          dx(:,dimy+j) = dx(:,dimy)
          dy(:,dimy+j) = dy(:,dimy)

          dx(:,j-slv)  = dx(:,1)
          dy(:,j-slv)  = dy(:,1)
       END DO
    END IF

    DO k=1, dimz
       dz(k) = delta_z(k)
    END DO

    IF (cycle_z) THEN
       DO k=1, slv
          dz(dimz+k) = dz(k)
          dz(k-slv)  = dz(dimz+k-slv)
       END DO
    ELSE
       dz(1-slv:0)         = dz(1)
       dz(dimz+1:dimz+slv) = dz(dimz)
    END IF

    ALLOCATE(depth(0:dimz))
    DO k=0, dimz-1
       depth(k) = sum(delta_z(k+1:dimz))
    END DO
    depth(dimz) = 0.0

  END SUBROUTINE read_geometry

!-----------------------------------------------------------------------------------------------------------------------

  LOGICAL PURE FUNCTION in_region_3d_bytype(rgn, r)
    TYPE(subregion_struct), INTENT(IN) :: rgn
    INTEGER, OPTIONAL,   INTENT(IN) :: r

#ifdef PARALLEL_MPI
    INTEGER :: rank_

    rank_ = rank
    IF (present(r)) rank_ = r

    in_region_3d_bytype = rgn%x_start <= (icoords(rank_)+1)*isize .AND. rgn%x_end >= icoords(rank_)*isize+1 .AND. &
                          rgn%y_start <= (jcoords(rank_)+1)*jsize .AND. rgn%y_end >= jcoords(rank_)*jsize+1 .AND. &
                          rgn%z_start <= (kcoords(rank_)+1)*ksize .AND. rgn%z_end >= kcoords(rank_)*ksize+1
#else
    in_region_3d_bytype = rgn%x_start <= dimx .AND. rgn%x_end >= 1 .AND. &
                          rgn%y_start <= dimy .AND. rgn%y_end >= 1 .AND. &
                          rgn%z_start <= dimz .AND. rgn%z_end >= 1
#endif
  END FUNCTION in_region_3d_bytype

  LOGICAL PURE FUNCTION in_region_3d_bycode(code, r)
    INTEGER, INTENT(IN) :: code
    INTEGER, OPTIONAL, INTENT(IN) :: r

    in_region_3d_bycode = in_region_3d_bytype(regions(code), r)
  END FUNCTION in_region_3d_bycode

  LOGICAL PURE FUNCTION in_region_2d_bytype(rgn, r)
    TYPE(subregion_struct), INTENT(IN) :: rgn
    INTEGER,             INTENT(IN), OPTIONAL :: r

#ifdef PARALLEL_MPI
    INTEGER :: rank_

    rank_ = rank
    IF (present(r)) rank_ = r

    in_region_2d_bytype = rgn%x_start <= (icoords(rank_)+1)*isize .AND. rgn%x_end >= icoords(rank_)*isize+1 .AND. &
                          rgn%y_start <= (jcoords(rank_)+1)*jsize .AND. rgn%y_end >= jcoords(rank_)*jsize+1
#else
    in_region_2d_bytype = rgn%x_start <= dimx .AND. rgn%x_end >= 1 .AND. &
                          rgn%y_start <= dimy .AND. rgn%y_end >= 1
#endif
  END FUNCTION in_region_2d_bytype

  LOGICAL PURE FUNCTION in_region_2d_bycode(code, r)
    INTEGER, INTENT(IN) :: code
    INTEGER, OPTIONAL, INTENT(IN) :: r

    in_region_2d_bycode = in_region_2d_bytype(regions(code), r)
  END FUNCTION in_region_2d_bycode

  LOGICAL PURE FUNCTION in_region_h(rgn, r)
    TYPE(subregion_struct), INTENT(IN) :: rgn
    INTEGER,             INTENT(IN), OPTIONAL :: r

#ifdef PARALLEL_MPI
    INTEGER :: rank_

    rank_ = hrank
    IF (present(r)) rank_ = r

    in_region_h = rgn%x_start <= (icoords_h(rank_)+1)*isize .AND. rgn%x_end >= icoords_h(rank_)*isize+1 .AND. &
                  rgn%y_start <= (jcoords_h(rank_)+1)*jsize .AND. rgn%y_end >= jcoords_h(rank_)*jsize+1
#else
    in_region_h = rgn%x_start <= dimx .AND. rgn%x_end >= 1 .AND. &
                  rgn%y_start <= dimy .AND. rgn%y_end >= 1
#endif
  END FUNCTION in_region_h

  LOGICAL PURE FUNCTION in_region_xz_bytype(rgn, r)
    TYPE(subregion_struct), INTENT(IN) :: rgn
    INTEGER,             INTENT(IN), OPTIONAL :: r

#ifdef PARALLEL_MPI
    INTEGER :: rank_

    rank_ = rank
    IF (present(r)) rank_ = r

    in_region_xz_bytype = rgn%x_start <= (icoords(rank_)+1)*isize .AND. rgn%x_end >= icoords(rank_)*isize+1 .AND. &
                          rgn%z_start <= (kcoords(rank_)+1)*ksize .AND. rgn%z_end >= kcoords(rank_)*ksize+1
#else
    in_region_xz_bytype = rgn%x_start <= dimx .AND. rgn%x_end >= 1 .AND. &
                          rgn%z_start <= dimz .AND. rgn%z_end >= 1
#endif
  END FUNCTION in_region_xz_bytype

  LOGICAL PURE FUNCTION in_region_xz_bycode(code, r)
    INTEGER, INTENT(IN) :: code
    INTEGER, OPTIONAL, INTENT(IN) :: r

    in_region_xz_bycode = in_region_xz_bytype(regions(code), r)
  END FUNCTION in_region_xz_bycode

  LOGICAL PURE FUNCTION in_region_yz_bytype(rgn, r)
    TYPE(subregion_struct), INTENT(IN) :: rgn
    INTEGER,             INTENT(IN), OPTIONAL :: r

#ifdef PARALLEL_MPI
    INTEGER :: rank_

    rank_ = rank
    IF (present(r)) rank_ = r

    in_region_yz_bytype = rgn%y_start <= (jcoords(rank_)+1)*jsize .AND. rgn%y_end >= jcoords(rank_)*jsize+1 .AND. &
                          rgn%z_start <= (kcoords(rank_)+1)*ksize .AND. rgn%z_end >= kcoords(rank_)*ksize+1
#else
    in_region_yz_bytype = rgn%y_start <= dimy .AND. rgn%y_end >= 1 .AND. &
                          rgn%z_start <= dimz .AND. rgn%z_end >= 1
#endif
  END FUNCTION in_region_yz_bytype

  LOGICAL PURE FUNCTION in_region_yz_bycode(code, r)
    INTEGER, INTENT(IN) :: code
    INTEGER, OPTIONAL, INTENT(IN) :: r

    in_region_yz_bycode = in_region_yz_bytype(regions(code), r)
  END FUNCTION in_region_yz_bycode

  LOGICAL PURE FUNCTION in_region_x_bytype(rgn, r)
    TYPE(subregion_struct), INTENT(IN) :: rgn
    INTEGER,             INTENT(IN), OPTIONAL :: r

#ifdef PARALLEL_MPI
    INTEGER :: rank_

    rank_ = rank
    IF (present(r)) rank_ = r

    in_region_x_bytype = rgn%x_start <= (icoords(rank_)+1)*isize .AND. rgn%x_end >= icoords(rank_)*isize+1
#else
    in_region_x_bytype = rgn%x_start <= dimx .AND. rgn%x_end >= 1
#endif
  END FUNCTION in_region_x_bytype

  LOGICAL PURE FUNCTION in_region_x_bycode(code, r)
    INTEGER, INTENT(IN) :: code
    INTEGER, OPTIONAL, INTENT(IN) :: r

    in_region_x_bycode = in_region_x_bytype(regions(code), r)
  END FUNCTION in_region_x_bycode

  LOGICAL PURE FUNCTION in_region_y_bytype(rgn, r)
    TYPE(subregion_struct), INTENT(IN) :: rgn
    INTEGER,             INTENT(IN), OPTIONAL :: r

#ifdef PARALLEL_MPI
    INTEGER :: rank_

    rank_ = rank
    IF (present(r)) rank_ = r

    in_region_y_bytype = rgn%y_start <= (jcoords(rank_)+1)*jsize .AND. rgn%y_end >= jcoords(rank_)*jsize+1
#else
    in_region_y_bytype = rgn%y_start <= dimy .AND. rgn%y_end >= 1
#endif
  END FUNCTION in_region_y_bytype

  LOGICAL PURE FUNCTION in_region_y_bycode(code, r)
    INTEGER, INTENT(IN) :: code
    INTEGER, OPTIONAL, INTENT(IN) :: r

    in_region_y_bycode = in_region_y_bytype(regions(code), r)
  END FUNCTION in_region_y_bycode

  LOGICAL PURE FUNCTION in_region_z_bytype(rgn, r)
    TYPE(subregion_struct), INTENT(IN) :: rgn
    INTEGER,             INTENT(IN), OPTIONAL :: r

#ifdef PARALLEL_MPI
    INTEGER :: rank_

    rank_ = rank
    IF (present(r)) rank_ = r

    in_region_z_bytype = rgn%z_start <= (kcoords(rank_)+1)*ksize .AND. rgn%z_end >= kcoords(rank_)*ksize+1
#else
    in_region_z_bytype = rgn%z_start <= dimz .AND. rgn%z_end >= 1
#endif
  END FUNCTION in_region_z_bytype

  LOGICAL PURE FUNCTION in_region_z_bycode(code, r)
    INTEGER, INTENT(IN) :: code
    INTEGER, OPTIONAL, INTENT(IN) :: r

    in_region_z_bycode = in_region_z_bytype(regions(code), r)
  END FUNCTION in_region_z_bycode

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE checkpoint(n)
    INTEGER, INTENT(IN) :: n

    CALL barrier

    IF (rank==0) WRITE(REPORT_UNIT, *) "----- Check point No.", n, " -------------"

    CALL barrier

  END SUBROUTINE checkpoint

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE report_ssh(force)
    LOGICAL, INTENT(IN), OPTIONAL :: force

    REAL(8) :: avol, tmp, max, min

    INTEGER :: i, j
    INTEGER :: ierr

    IF (rigid_lid) RETURN

    IF (.NOT. present(force)) THEN
       IF (ssh_report_interval <= 0)              RETURN
       IF (mod(n_timestep, ssh_report_interval)/=0) RETURN
    ELSE
       IF (.NOT. force) RETURN
    END IF

    avol = 0.0
    IF (vrank==0) THEN
       DO i=1, isize
       DO j=1, jsize
          avol = avol + ssh(i,j)*dsz(i,j)*imask2d(i,j)
       END DO
       END DO
    END IF

    CALL gsum(avol)

    max = maxval(ssh(1:isize, 1:jsize), MASK=lmask2d(1:isize,1:jsize))
    min = minval(ssh(1:isize, 1:jsize), MASK=lmask2d(1:isize,1:jsize))

    CALL gmax(max)
    CALL gmin(min)

    IF (rank==0) THEN
       WRITE(REPORT_UNIT, '(A, ": SSH max=", ES10.3, ", min=", ES10.3, ", ave=", ES10.3, " [m], volume anomaly=", ES13.6, "[m^3]")') &
            trim(current_datetime), max, min, avol/total_area, avol
    ENDIF

  END SUBROUTINE report_ssh

END MODULE geometry
