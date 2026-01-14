#include "macro.h"
#ifndef PARALLEL_MPI
#undef MPIIO
#endif

MODULE particles
  USE misc
  USE geometry
  USE io
  USE parameters
  USE calendar
  IMPLICIT NONE
  SAVE
  PRIVATE
  PUBLIC particle_property_index, particle_category_index
  PUBLIC init_particles, finalize_particles, step_particles, flush_particles
  PUBLIC create_particle, delete_particle, delete_particles, search_particle
  PUBLIC particle_struct, n_prop
  PUBLIC write_particles, read_particles, decode_particles, encode_particles
  PUBLIC get_particle, particle_handle, particle_handle_next, set_particle_property, get_particle_property, get_particle_category
  PUBLIC get_particle_id, get_particle_xpos, get_particle_ypos, get_particle_zpos
  PUBLIC get_default_property
  PUBLIC set_particle_pos_in_grid, set_particle_xpos, set_particle_ypos, set_particle_zpos
  PUBLIC check_particles
  PUBLIC particle_driven_forcing
  PUBLIC terminal_velocity
  PUBLIC size_of_particle
  PUBLIC particle_density
  PUBLIC particle_histogram
  PUBLIC PROP_KIND

  INTEGER, PARAMETER :: max_categories = 256
  INTEGER, PARAMETER :: max_input      = 256
  INTEGER, PARAMETER :: max_tracks     = 1024
  INTEGER, PARAMETER :: max_initialfile= 8

  INTEGER, PARAMETER :: PROP_KIND = 4
  INTEGER, PARAMETER :: n_prop    = 6
  INTEGER, PARAMETER :: size_of_particle = 8 + 4*3 + 4*3 + 4 + PROP_KIND*n_prop + 4

  INTEGER :: n_categories = 0
  INTEGER :: n_tracks     = 0
  INTEGER :: n_input      = 0

  INTEGER, ALLOCATABLE :: ncreate(:,:,:)

  LOGICAL :: mpiio = .TRUE.

  TYPE particle_struct
     INTEGER(8) :: id = 0

     INTEGER(4) :: ipos = 0
     INTEGER(4) :: jpos = 0
     INTEGER(4) :: kpos = 0

     REAL(4)    :: xpos = 0.0
     REAL(4)    :: ypos = 0.0
     REAL(4)    :: zpos = 0.0

     INTEGER(4) :: category = 0

     REAL(PROP_KIND) :: property(n_prop) = UNDEF

     INTEGER(4) :: seed = 0

     INTEGER :: track_index = -1
  END TYPE particle_struct

  INTEGER(8), ALLOCATABLE :: p_id(:)
  INTEGER, ALLOCATABLE    :: p_ijk(:,:)
  REAL(4), ALLOCATABLE    :: p_xyz(:,:)
  INTEGER, ALLOCATABLE    :: p_category(:)
  REAL(PROP_KIND), ALLOCATABLE :: p_property(:,:)
  INTEGER, ALLOCATABLE    :: p_seed(:)
  INTEGER, ALLOCATABLE    :: p_track_index(:)
  INTEGER, ALLOCATABLE    :: p_next(:)
  INTEGER, ALLOCATABLE    :: p_prev(:)

  INTEGER, PARAMETER :: n_probe = 8
  CHARACTER(6), PARAMETER :: probe_vars(n_probe) = (/"U     ", "V     ", "W     ", "KE    ", "SIGMA0", "RHO   ", "DIVH  ", "ZETA  " /)

  TYPE category_info_struct
     CHARACTER(32)   :: name
     INTEGER         :: code
     CHARACTER(32)   :: property_name(n_prop)
     REAL(PROP_KIND) :: property_default(n_prop)
     LOGICAL         :: fix_x
     LOGICAL         :: fix_y
     LOGICAL         :: fix_z
     REAL(8)         :: lifetime
     REAL(4)         :: rho
     REAL(4)         :: repul
     REAL(4)         :: repul_sfc
     REAL(4)         :: repul_btm
     LOGICAL         :: rmv_sfc
     LOGICAL         :: rmv_btm
     LOGICAL         :: rmv_opn
     REAL(4)         :: disph
     REAL(4)         :: dispv
     REAL(4)         :: dispr(3)
     LOGICAL         :: dispgrad
     REAL(4)         :: fish_kinesis(4)
     REAL(4)         :: tkesource
     LOGICAL         :: rk4
     INTEGER         :: interp
     LOGICAL         :: record_grave
     LOGICAL         :: record_birth
     LOGICAL         :: record_watch
     LOGICAL         :: buoyancy
     INTEGER         :: settling
     REAL(4)         :: hist_bins(32)
     REAL(8)         :: hist_interval
     INTEGER         :: hist_prop
     LOGICAL         :: stdrift
     INTEGER         :: propindex_u
     INTEGER         :: propindex_v
     INTEGER         :: propindex_w
     INTEGER         :: propindex_diam
     INTEGER         :: propindex_mass
     INTEGER         :: propindex_probe(n_probe)
     INTEGER         :: propindex_probe_tracer(2,n_prop)
     INTEGER         :: propindex_cumul_tracer(2,n_prop)
     REAL(8), POINTER:: aspawn(:,:,:)
  END TYPE category_info_struct

  TYPE(category_info_struct) :: category_info(0:max_categories-1)

  INTEGER, ALLOCATABLE :: particle_ptr(:,:,:,:)
  INTEGER, ALLOCATABLE :: particle_cnt(:,:,:,:)

  INTEGER :: heap_ptr
  INTEGER :: heap_cnt

  INTEGER :: max_particles = 0

  INTEGER(8) :: max_id = 0
  INTEGER(8) :: cur_id = 0
  INTEGER(8) :: lim_id = 0


#ifdef PARALLEL_MPI
  INTEGER :: buffersize = 16384

  INTEGER(1), ALLOCATABLE :: sendbuf(:,:)
  INTEGER(1), ALLOCATABLE :: recvbuf(:,:)

  INTEGER :: win_max_id
#endif

  INTEGER, PARAMETER :: chunksize  = 65536
  INTEGER :: ompchunk = -1

  TYPE record_struct
     LOGICAL :: initialized = .FALSE.
     INTEGER :: n_write     = 0
     INTEGER :: count       = 0
#ifdef MPIIO
     INTEGER :: filehandle
#endif
     CHARACTER(1024) :: filepath
#ifdef NO_ALLOCATABLE_IN_TYPE
     INTEGER(1), POINTER     :: buffer(:,:)
#else
     INTEGER(1), ALLOCATABLE :: buffer(:,:)
#endif
     INTEGER :: bufsize
     LOGICAL :: csv      = .FALSE.
     REAL(8) :: start    = UNDEF
     REAL(8) :: end      = UNDEF
     REAL(8) :: interval = UNDEF
     REAL(8) :: lastsync = UNDEF
  END TYPE record_struct

  INTEGER(8) :: track_id( max_tracks) = 0

  TYPE(record_struct) :: track_record(max_tracks)

  TYPE(record_struct) :: grave_record
  TYPE(record_struct) :: birth_record
  TYPE(record_struct) :: watch_record

  INTEGER(8) :: min_track_id = 0
  INTEGER(8) :: max_track_id = 0

  ! control parameters for avection and repulsion, default values are defined in the parameters.F90
  LOGICAL :: default_advection_rk4
  REAL(4) :: default_repulsion_surface
  REAL(4) :: default_repulsion_bottom
  REAL(4) :: default_repulsion_side

  REAL(4) :: npzd_permoln ! number of NPZD particle (molN^-1) used for initial and boundary conditions
  REAL(4) :: npzd_svsigma ! normalized standard deviation for sinking velocity of DET
  REAL(4) :: npzd_dfactor ! factor for deployment of NPZD particle (inverse of the restore-timescale by timesteps)

  INTERFACE particle_category_index
     MODULE PROCEDURE category_index_byname
     MODULE PROCEDURE category_index_bycode
  END INTERFACE particle_category_index

  INTERFACE category_index
     MODULE PROCEDURE category_index_byname
     MODULE PROCEDURE category_index_bycode
  END INTERFACE category_index

  INTERFACE particle_property_index
     MODULE PROCEDURE property_index
  END INTERFACE particle_property_index

  TYPE output_info_struct
     REAL(8) :: start
     REAL(8) :: end
     REAL(8) :: interval
     REAL(8) :: lastwrite
     CHARACTER(1024) :: filepath
     LOGICAL :: omit_time
  END TYPE output_info_struct

  TYPE(output_info_struct) :: output_info

  TYPE input_info_struct
     CHARACTER(1024):: filepath
     LOGICAL        :: csv
     REAL(8)        :: start
     REAL(8)        :: end
     REAL(8)        :: interval
     REAL(8)        :: lastread
     REAL(4)        :: rate(0:1023)
     REAL(8)        :: residual
     INTEGER        :: mode
     INTEGER        :: periods
     INTEGER(8)     :: count
     INTEGER(4)     :: length
     LOGICAL        :: cyclic
     INTEGER(8)     :: cyclic_id
     LOGICAL        :: report
     LOGICAL        :: initialized
     LOGICAL        :: rate_rel
     INTEGER(8)     :: offset_id
     LOGICAL        :: ignore_id
     INTEGER        :: def_cat
#ifdef NO_ALLOCATABLE_IN_TYPE
     TYPE(particle_struct), POINTER     :: particles(:)
#else
     TYPE(particle_struct), ALLOCATABLE :: particles(:)
#endif
  END TYPE input_info_struct

  TYPE(input_info_struct) :: input_info(max_input)

  LOGICAL :: native_bigendian

  INTEGER :: nlimit(3) = 0.0

  INTEGER, PARAMETER :: cat_timestamp = -1

  INTEGER, PARAMETER :: id_nshift = 28

CONTAINS

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE init_particles
    INTEGER :: ptr

    INTEGER(1) :: buf(size_of_particle)

    CHARACTER(512) :: defaultdir
    CHARACTER(512) :: outputdir

    CHARACTER(512) :: output_name
    CHARACTER(16)  :: output_start
    CHARACTER(16)  :: output_end
    CHARACTER(16)  :: output_interval

    CHARACTER(512) :: initialdir
    CHARACTER(512) :: initialfile(max_initialfile)
    LOGICAL        :: initialfile_csv(max_initialfile)
    INTEGER        :: initialfile_default_category(max_initialfile)
    INTEGER(8)     :: initialfile_offset_id(max_initialfile)

    LOGICAL        :: gravefile
    LOGICAL        :: gravefile_csv
    LOGICAL        :: gravefile_append
    CHARACTER(512) :: gravefile_name
    CHARACTER(16)  :: gravefile_start
    CHARACTER(16)  :: gravefile_end
    CHARACTER(16)  :: gravefile_interval

    LOGICAL        :: birthfile
    LOGICAL        :: birthfile_csv
    LOGICAL        :: birthfile_append
    CHARACTER(512) :: birthfile_name
    CHARACTER(16)  :: birthfile_start
    CHARACTER(16)  :: birthfile_end
    CHARACTER(16)  :: birthfile_interval

    LOGICAL        :: watchfile
    LOGICAL        :: watchfile_csv
    LOGICAL        :: watchfile_append
    CHARACTER(512) :: watchfile_name
    CHARACTER(16)  :: watchfile_start
    CHARACTER(16)  :: watchfile_end
    CHARACTER(16)  :: watchfile_interval

    LOGICAL        :: npzd_init    ! switch for initialzing NPZD partieles

    LOGICAL        :: omit_time

    INTEGER :: i, j, k
    INTEGER :: ierr
    INTEGER :: iostat
    CHARACTER(256) :: iomsg

    NAMELIST / particles / &
         use_particles,    &
         outputdir,        &
         output_name,      &
         output_start,     &
         output_end,       &
         output_interval,  &
         initialdir,       &
         initialfile,      &
         initialfile_csv,  &
         initialfile_default_category, &
         initialfile_offset_id, &
         gravefile,        &
         gravefile_csv,    &
         gravefile_append, &
         gravefile_name,   &
         gravefile_start,  &
         gravefile_end,    &
         gravefile_interval,&
         birthfile,        &
         birthfile_csv,    &
         birthfile_append, &
         birthfile_name,   &
         birthfile_interval,&
         watchfile,        &
         watchfile_csv,    &
         watchfile_append, &
         watchfile_name,   &
         watchfile_start,  &
         watchfile_end,    &
         watchfile_interval,&
         mpiio,            &
         npzd_permoln,     &
         npzd_svsigma,     &
         npzd_dfactor,     &
         npzd_init,        &
         omit_time

    native_bigendian = check_endian()

!$  ompchunk = INT(isize*jsize / (nthreads*16))

    CALL default_output( &
         default_outputdir = defaultdir,      &
         default_start     = output_start,    &
         default_end       = output_end,      &
         default_interval  = output_interval)

    outputdir        = defaultdir
    output_name      = 'PARTICLES'

    gravefile        = .FALSE.
    gravefile_name   = 'PARTICLE_GRAVE'
    gravefile_csv    = .FALSE.
    gravefile_append = .TRUE.
    gravefile_start    = output_start
    gravefile_end      = output_end
    gravefile_interval = output_interval

    birthfile        = .FALSE.
    birthfile_name   = 'PARTICLE_BIRTH'
    birthfile_csv    = .FALSE.
    birthfile_append = .TRUE.
    birthfile_start    = output_start
    birthfile_end      = output_end
    birthfile_interval = output_interval


    watchfile        = .FALSE.
    watchfile_name   = 'PARTICLE_WATCH'
    watchfile_csv    = .FALSE.
    watchfile_append = .TRUE.
    watchfile_start    = output_start
    watchfile_end      = output_end
    watchfile_interval = output_interval

    initialdir       = global_initialdir
    initialfile(:)     = ''
    initialfile_csv(:) = .FALSE.
    initialfile_default_category(:) = 0
    initialfile_offset_id(:)        = 0_8

    npzd_init        = .FALSE.
    npzd_permoln     = 100
    npzd_svsigma     = 0.2
    npzd_dfactor     = 0.1

    omit_time        = .FALSE.

#ifdef MPIIO
    mpiio            = .TRUE.
#else
    mpiio            = .FALSE.
#endif

    IF (rank==0) THEN
       REWIND(CONFIG_UNIT)
       READ(CONFIG_UNIT, NML=particles, IOSTAT=iostat, IOMSG=iomsg)
       CALL assert(iostat <= 0, "failed to read PARTICLE namelist", iomsg)
    END IF

    CALL bcast(use_particles)
    IF (.NOT. use_particles) RETURN

    CALL bcast(outputdir)

    CALL replace_vars(outputdir, default=defaultdir)

    CALL bcast(output_name)

    output_info%filepath = path(outputdir, trim(output_name))

    CALL bcast(output_start)
    CALL bcast(output_end)
    CALL bcast(output_interval)

    output_info%start = datetime_seconds(output_start)
    IF (output_end(1:1)=='+') THEN
       output_info%end = output_info%start + interval_seconds(output_end(2:)//' ')
    ELSE
       output_info%end = datetime_seconds(output_end)
    END IF
    output_info%interval  = interval_seconds(output_interval)
    output_info%lastwrite = output_info%start
    IF (omit_time) CALL assert(interval_seconds(output_interval) >= 86400.0, "Cannot specify True for OMIT_TIME (particles) when output interval is less than one day")
    output_info%omit_time = omit_time

    CALL bcast(gravefile)
    CALL bcast(gravefile_csv)
    CALL bcast(gravefile_append)
    CALL bcast(gravefile_name)
    CALL bcast(gravefile_start)
    CALL bcast(gravefile_end)
    CALL bcast(gravefile_interval)
    CALL replace_vars(gravefile_name, default=defaultdir)
    IF (gravefile_csv .AND. .NOT. is_csvfile(gravefile_name)) gravefile_name = trim(gravefile_name) // ".csv"

    CALL bcast(birthfile)
    CALL bcast(birthfile_csv)
    CALL bcast(birthfile_append)
    CALL bcast(birthfile_name)
    CALL bcast(birthfile_start)
    CALL bcast(birthfile_end)
    CALL bcast(birthfile_interval)
    CALL replace_vars(birthfile_name, default=defaultdir)
    IF (birthfile_csv .AND. .NOT. is_csvfile(birthfile_name)) birthfile_name = trim(birthfile_name) // ".csv"

    CALL bcast(watchfile)
    CALL bcast(watchfile_csv)
    CALL bcast(watchfile_append)
    CALL bcast(watchfile_name)
    CALL bcast(watchfile_start)
    CALL bcast(watchfile_end)
    CALL bcast(watchfile_interval)
    CALL replace_vars(watchfile_name, default=defaultdir)
    IF (watchfile_csv .AND. .NOT. is_csvfile(watchfile_name)) watchfile_name = trim(watchfile_name) // ".csv"


    IF (gravefile) CALL init_record(path(outputdir, gravefile_name), gravefile_csv, gravefile_append, gravefile_start, gravefile_end, gravefile_interval, grave_record)
    IF (birthfile) CALL init_record(path(outputdir, birthfile_name), birthfile_csv, birthfile_append, birthfile_start, birthfile_end, birthfile_interval, birth_record)
    IF (watchfile) CALL init_record(path(outputdir, watchfile_name), watchfile_csv, watchfile_append, watchfile_start, watchfile_end, watchfile_interval, watch_record)

    CALL bcast(npzd_permoln)
    CALL bcast(npzd_svsigma)
    CALL bcast(npzd_dfactor)
    CALL bcast(npzd_init)

    CALL bcast(initialdir)
    CALL replace_vars(initialdir,  default=defaultdir)

    DO i=1, max_initialfile
       CALL bcast(initialfile(i))
       CALL bcast(initialfile_csv(i))
       CALL bcast(initialfile_offset_id(i))
       CALL bcast(initialfile_default_category(i))
       CALL replace_vars(initialfile(i), default=defaultdir)
    END DO

    CALL bcast(omit_time)

    CALL bcast(mpiio)
#ifndef MPIIO
    CALL assert(.NOT. mpiio, "MPI-IO is not supported!")
#endif

    CALL read_category_namelist

    CALL read_track_namelist

    ALLOCATE(particle_ptr(0:isize+1, 0:jsize+1, 0:ksize+1, n_categories))
    particle_ptr(:,:,:,:) = NULL

    ALLOCATE(particle_cnt(0:isize+1, 0:jsize+1, 0:ksize+1, n_categories))
    particle_cnt(:,:,:,:) = 0

    heap_ptr = NULL
    heap_cnt = 0

    CALL increase(isize*jsize*ksize)

    CALL init_ncreate

    IF (trim(initialfile(1)) == '') THEN
       IF (initialdir=='') initialdir=defaultdir
       initialfile(1) = trim(output_name) //'.'// format_datetime(t_current, omit_time=omit_time)
    END IF

    DO i=1, max_initialfile
       IF (initialfile(i) /= "") THEN
          CALL read_particles(path(initialdir, initialfile(i)), csv=initialfile_csv(i), &
                              default_category=initialfile_default_category(i), offset_id=initialfile_offset_id(i))
       END IF
    END DO

    CALL check_particles

#ifdef PARALLEL_MPI
    ALLOCATE(sendbuf(size_of_particle,buffersize))
    ALLOCATE(recvbuf(size_of_particle,buffersize))

    CALL mpi_barrier(comm, ierr)
#endif

    CALL read_input_namelist

    IF (npzd_init) CALL init_npzd_particles

    IF (gravefile) CALL flush_record(grave_record, all=.TRUE., timestamp=.TRUE.)
    IF (birthfile) CALL flush_record(birth_record, all=.TRUE., timestamp=.TRUE.)
    IF (watchfile) CALL flush_record(watch_record, all=.TRUE., timestamp=.TRUE.)

  CONTAINS

    SUBROUTINE read_track_namelist
      INTEGER(8)      :: id
      CHARACTER(512)  :: outputdir
      CHARACTER(128)  :: filename
      LOGICAL         :: csv
      LOGICAL         :: append
      LOGICAL         :: initialize
      INTEGER         :: num
      INTEGER         :: stride
      CHARACTER(16)   :: start
      CHARACTER(16)   :: end

      CHARACTER(512)  :: default_outputdir
      CHARACTER(16)   :: default_start
      CHARACTER(16)   :: default_end

      INTEGER         :: n
      CHARACTER(20)   :: str_id

      NAMELIST / particle_track / &
           id, num, stride, outputdir, filename, csv, append, start, end,initialize

      CALL default_output(default_outputdir=default_outputdir, default_start=default_start, default_end=default_end)

      IF (rank==0) REWIND(CONFIG_UNIT)
      DO
         IF (rank==0) THEN
            id = 0
            num    = 1
            stride = 1
            outputdir  = default_outputdir
            filename   = ''
            initialize = .FALSE.
            csv        = .FALSE.
            append     = .TRUE.
            start      = default_start
            end        = default_end

            READ(CONFIG_UNIT, NML=particle_track, IOSTAT=iostat, IOMSG=iomsg)
            CALL assert(iostat <= 0, "failed to read PARTICLE_TRACK namelist", iomsg)

            csv    = csv .OR. is_csvfile(filename)
            append = append .AND. .NOT. initialize
         END IF

         CALL bcast(iostat)
         IF (iostat < 0) EXIT

         CALL bcast(id)
         CALL bcast(num)
         CALL bcast(stride)
         CALL bcast(outputdir)
         CALL bcast(filename)
         CALL bcast(csv)
         CALL bcast(append)
         CALL bcast(start)
         CALL bcast(end)

         CALL assert(id>0, "ID in PARTICLE_TRACK should be a positive integer")

         IF (trim(filename) == '') THEN
            DO n=0, num-1
               WRITE(str_id, '(I0)') id + n*stride
               filename = 'PARTICLE_TRACK.' // trim(str_id)

               CALL register_track(id + n*stride, trim(path(outputdir, filename)), csv, append, start, end)
            END DO
         ELSE
            CALL assert(num==1, "assinging FILENAME for PARTICLE_TRACK with NUM > 1 is not allowed")
            CALL register_track(id, trim(path(outputdir, filename)), csv, append, start, end)
         END IF

      END DO

    END SUBROUTINE read_track_namelist

    SUBROUTINE read_category_namelist
      USE tracers,  ONLY: tracer, tracer_name, n_tracer, tracer_index
      USE npzd,     ONLY: use_convrate
      CHARACTER(32)   :: category_name
      INTEGER         :: category_code
      CHARACTER(32)   :: property_name(n_prop)
      REAL(PROP_KIND) :: property_default(n_prop)
      LOGICAL         :: fix_x
      LOGICAL         :: fix_y
      LOGICAL         :: fix_z
      REAL(4)         :: rho
      CHARACTER(16)   :: lifetime
      LOGICAL         :: advection_rk4
      INTEGER         :: advection_interp
      REAL(4)         :: repulsion_surface
      REAL(4)         :: repulsion_bottom
      REAL(4)         :: repulsion_side
      LOGICAL         :: remove_surface
      LOGICAL         :: remove_bottom
      LOGICAL         :: remove_open
      REAL(4)         :: disph, dispersion_kh
      REAL(4)         :: dispv, dispersion_kv
      LOGICAL         :: dispgrad
      REAL(4)         :: fish_kinesis(4) ! kinesis migration scheme (Humston 2000), update interval[s], prefered_temp_lwoer/upper [degC], swim_speed [m/s]
      REAL(4)         :: dispersion_sc(3), dispersion_pr(3)
      REAL(4)         :: tkesource
      LOGICAL         :: record_gravefile
      LOGICAL         :: record_birthfile
      LOGICAL         :: record_watchfile
      LOGICAL         :: buoyancy_coupling
      INTEGER         :: settling_scheme  !0 : no settling velocity
                                          !1 : combination of Newton's, Allen's, and Stokes's formulas (selected by Raynolds number range)
                                          !2 : following Rubey's formula for suspended sediment matter (c.f. Rubey et al. 1933)
                                          !3 : following Stokes's law (proportional to d^2 for low Raynolds number)
                                          !4 : following Allen's law (proportional to d for mid Raynolds number)
                                          !5 : following Newton's law (proportional to d^{1/2} for high Raynolds number)
                                          !6 : following Kaiser's empilical formula for micro-plastic particles (c.f. Kaiser et al. 2019)
      LOGICAL         :: settling_rubey1933 ! for backward compatibility

      INTEGER         :: histogram_property  ! for which property histogram is calculated
      CHARACTER(16)   :: histogram_interval  ! interval time for histogram calculation
      REAL(4)         :: histogram_bins(32)  ! bin thresholds for histogram (max 32bins)
      LOGICAL         :: stokes_drift

      INTEGER :: cat
      INTEGER :: n, m, i
      INTEGER :: iostat

      NAMELIST / particle_category / &
           category_name,     &
           category_code,     &
           rho,               &
           fix_x,             &
           fix_y,             &
           fix_z,             &
           lifetime,          &
           advection_rk4,     &
           advection_interp,  &
           repulsion_surface, &
           repulsion_bottom,  &
           repulsion_side,    &
           remove_surface,    &
           remove_bottom,     &
           remove_open,       &
           dispersion_kh,     &
           dispersion_kv,     &
           disph,             &
           dispv,             &
           dispgrad,          &
           fish_kinesis,      &
           dispersion_pr,     &
           dispersion_sc,     &
           tkesource,         &
           property_name,     &
           property_default,  &
           record_gravefile,  &
           record_birthfile,  &
           record_watchfile,  &
           buoyancy_coupling, &
           settling_scheme,   &
           settling_rubey1933,&
           histogram_property,&
           histogram_interval,&
           histogram_bins,    &
           stokes_drift

      NAMELIST / particle_category_default / &
           category_name,     &
           rho,               &
           fix_x,             &
           fix_y,             &
           fix_z,             &
           lifetime,          &
           advection_rk4,     &
           advection_interp,  &
           repulsion_surface, &
           repulsion_bottom,  &
           repulsion_side,    &
           remove_surface,    &
           remove_bottom,     &
           remove_open,       &
           dispersion_kh,     &
           dispersion_kv,     &
           disph,             &
           dispv,             &
           dispgrad,          &
           fish_kinesis,      &
           dispersion_pr,     &
           dispersion_sc,     &
           tkesource,         &
           property_name,     &
           property_default,  &
           record_gravefile,  &
           record_birthfile,  &
           record_watchfile,  &
           buoyancy_coupling, &
           settling_scheme,   &
           stokes_drift

      IF (rank==0) THEN
         REWIND(CONFIG_UNIT)
         property_name(:)    = ''
         property_default(:) = 0.0_PROP_KIND
         rho                 = UNDEF
         lifetime            = ''
         fix_x               = .FALSE.
         fix_y               = .FALSE.
         fix_z               = .FALSE.
         advection_rk4       = particle_advection_rk4
         advection_interp    = particle_advection_interp
         repulsion_surface   = particle_repulsion_surface
         repulsion_bottom    = particle_repulsion_bottom
         repulsion_side      = particle_repulsion_side
         remove_surface      = .FALSE.
         remove_bottom       = .FALSE.
         remove_open         = .TRUE.
         dispersion_kh       = UNDEF ! for backward compatibility
         dispersion_kv       = UNDEF ! for backward compatibility
         disph               = 0.0
         dispv               = 0.0
         dispgrad            = .TRUE.
         fish_kinesis        = (/UNDEF, UNDEF, UNDEF, UNDEF/)
         dispersion_pr       = (/UNDEF, UNDEF, UNDEF/)
         dispersion_sc       = (/0.0, UNDEF, 0.0/)
         tkesource           = 0.0
         record_gravefile    = .TRUE.
         record_birthfile    = .TRUE.
         record_watchfile    = .TRUE.
         buoyancy_coupling   = .TRUE.
         settling_scheme     = 0
         settling_rubey1933  = .FALSE.
         stokes_drift        = .TRUE.

         READ(CONFIG_UNIT, NML=particle_category_default, IOSTAT=iostat, IOMSG=iomsg)
         CALL assert(iostat <= 0, "failed to read PARTICLE_CATEGORY_DEFAULT namelist", iomsg)

         category_info(0)%code                = 0
         category_info(0)%name                = "PTCL"
         category_info(0)%property_name(:)    = property_name(:)
         category_info(0)%property_default(:) = property_default(:)
         category_info(0)%fix_x               = fix_x
         category_info(0)%fix_y               = fix_y
         category_info(0)%fix_z               = fix_z
         category_info(0)%rho                 = rho
         category_info(0)%rk4                 = advection_rk4
         category_info(0)%interp              = advection_interp
         category_info(0)%repul               = min(max(repulsion_side,    0.0), 0.5)
         category_info(0)%repul_sfc           = min(max(repulsion_surface, 0.0), 0.5)
         category_info(0)%repul_btm           = min(max(repulsion_bottom,  0.0), 0.5)
         category_info(0)%rmv_sfc             = remove_surface
         category_info(0)%rmv_btm             = remove_bottom
         category_info(0)%rmv_opn             = remove_open
         IF (dispersion_kh /= UNDEF) disph    = dispersion_kh
         IF (dispersion_kv /= UNDEF) dispv    = dispersion_kv
         category_info(0)%disph               = disph
         category_info(0)%dispv               = dispv
         category_info(0)%dispgrad            = dispgrad
         category_info(0)%fish_kinesis(:)     = fish_kinesis(:)

         IF (dispersion_pr(1) /= UNDEF) dispersion_sc(1) = dispersion_pr(1)
         IF (dispersion_pr(2) /= UNDEF) dispersion_sc(2) = dispersion_pr(2)
         IF (dispersion_pr(3) /= UNDEF) dispersion_sc(3) = dispersion_pr(3)

         !category_info%dispr represents foctor of disspersion coefficient with respect to eddy viscosity (i.e. inverse of the Schmidt number)
         IF (dispersion_sc(1) == 0.0) THEN
            category_info(0)%dispr(1) = 0.0
         ELSE
            category_info(0)%dispr(1) = 1.0 / dispersion_sc(1)
         END IF
         IF (dispersion_sc(2) == UNDEF) THEN
            category_info(0)%dispr(2) = category_info(0)%dispr(1)
         ELSE IF (dispersion_sc(2) == 0.0) THEN
            category_info(0)%dispr(2) = 0.0
         ELSE
            category_info(0)%dispr(2) = 1.0/dispersion_sc(2)
         END IF
         IF (dispersion_sc(3) == 0.0) THEN
            category_info(0)%dispr(3) = 0.0
         ELSE
            category_info(0)%dispr(3) = 1.0/dispersion_sc(3)
         END IF

         category_info(0)%tkesource           = tkesource
         category_info(0)%record_grave        = record_gravefile
         category_info(0)%record_birth        = record_birthfile
         category_info(0)%record_watch        = record_watchfile
         category_info(0)%buoyancy            = buoyancy_coupling
         category_info(0)%settling            = settling_scheme
         category_info(0)%stdrift             = stokes_drift

         IF (lifetime=='') THEN
            category_info(0)%lifetime = UNDEF
         ELSE
            category_info(0)%lifetime = interval_seconds(lifetime)
         END IF

         DO m=1, n_prop
            IF (trim(property_name(m)) /= '') CALL assert(valid_name(property_name(m)), "PROPERTY_NAME '"//trim(property_name(m))//"' is invalid")
         END DO

         REWIND(CONFIG_UNIT)
         DO
            category_name       = ''
            category_code       = 0
            property_name(:)    = category_info(0)%property_name(:)
            property_default(:) = category_info(0)%property_default(:)
            rho                 = category_info(0)%rho
            lifetime            = ''
            fix_x               = category_info(0)%fix_x
            fix_y               = category_info(0)%fix_y
            fix_z               = category_info(0)%fix_z
            advection_rk4       = category_info(0)%rk4
            advection_interp    = category_info(0)%interp
            repulsion_side      = category_info(0)%repul
            repulsion_surface   = category_info(0)%repul_sfc
            repulsion_bottom    = category_info(0)%repul_btm
            remove_surface      = category_info(0)%rmv_sfc
            remove_bottom       = category_info(0)%rmv_btm
            remove_open         = category_info(0)%rmv_opn
            dispersion_kh       = UNDEF ! for backward compatibility
            dispersion_kv       = UNDEF ! for backward compatibility
            disph               = category_info(0)%disph
            dispv               = category_info(0)%dispv
            dispersion_pr(:)    = UNDEF
            dispersion_sc(:)    = UNDEF
            dispgrad            = category_info(0)%dispgrad
            fish_kinesis(:)     = category_info(0)%fish_kinesis(:)
            tkesource           = category_info(0)%tkesource
            record_gravefile    = category_info(0)%record_grave
            record_birthfile    = category_info(0)%record_birth
            record_watchfile    = category_info(0)%record_watch
            buoyancy_coupling   = category_info(0)%buoyancy
            settling_scheme     = category_info(0)%settling
            settling_rubey1933  = .FALSE.
            stokes_drift        = category_info(0)%stdrift
            histogram_property  = 0
            histogram_interval  = "1D"
            histogram_bins(:)   = UNDEF

            READ(CONFIG_UNIT, NML=particle_category, IOSTAT=iostat, IOMSG=iomsg)
            CALL assert(iostat <= 0, "failed to read PARTICLE_CATEGORY namelist", iomsg)

            IF (iostat < 0) EXIT

            CALL assert(category_code >  0,               "CATEGORY_CODE should be a positive integer")
            CALL assert(category_index(category_code)==0, "CATEGORY_CODE for '"//trim(category_name)//"' is already used")

            DO m=1, n_prop
               IF (trim(property_name(m)) /= '') CALL assert(valid_name(property_name(m)), "PROPERTY_NAME '"//trim(property_name(m))//"' is invalid")
            END DO

            n_categories = n_categories + 1
            CALL assert(n_categories < max_categories, "number of particle categories should be less than MAX_CATEGORIES")

            cat = n_categories



            IF (category_name=="")  category_name = trim(category_info(0)%name) // trim(format(category_code, 'I0'))

            CALL assert(valid_name(    category_name),    "CATEGORY_NAME '"//trim(category_name)//"' is invalid")
            CALL assert(category_index(category_name)==0, "CATEGORY_NAME '"//trim(category_name)//"' is defined multiple times")

            category_info(cat)%name = category_name
            category_info(cat)%code = category_code


            IF (rho == UNDEF) THEN
               IF (category_name(1:len('SEDIMENT'))=='SEDIMENT') THEN
                  rho = rho_sediment
               ELSE IF (category_name(1:len('BUBBLE'))=='BUBBLE') THEN
                  rho = rho_bubble
               ELSE IF (category_name(1:len('FRAZIL'))=='FRAZIL') THEN
                  rho = rho_frazil
               ELSE
                  rho = rho_0
               END IF
            END IF

            IF (trim(category_name) == 'NPZD') THEN
               CALL assert(use_npzd, "NPZD particle requires USE_NPZD = .TRUE.")
               use_convrate = .TRUE.

               CALL assert(property_name(1) == "" .OR. property_name(1) == 'NPZD_STATUS', "property(1) for NPZD particle should be 'NPZD_STATUS'")
               CALL assert(property_name(2) == "" .OR. property_name(2) == 'COUNT_N2PZD', "property(2) for NPZD particle should be 'COUNT_N2PZD'")
               CALL assert(property_name(3) == "" .OR. property_name(3) == 'W',           "property(3) for NPZD particle should be 'W'")

               property_name(1) = 'NPZD_STATUS'
               property_name(2) = 'COUNT_N2PZD'
               property_name(3) = 'W'
            END IF

            IF (lifetime=='') THEN
               category_info(cat)%lifetime = category_info(0)%lifetime
            ELSE
               category_info(cat)%lifetime = interval_seconds(lifetime)
            END IF

            category_info(cat)%property_name(:)    = property_name(:)
            category_info(cat)%property_default(:) = property_default(:)
            category_info(cat)%fix_x               = fix_x
            category_info(cat)%fix_y               = fix_y
            category_info(cat)%fix_z               = fix_z
            category_info(cat)%rho                 = max(rho, 0.0)
            category_info(cat)%rk4                 = advection_rk4
            category_info(cat)%interp              = advection_interp
            category_info(cat)%repul               = min(max(repulsion_side,    0.0), 0.5)
            category_info(cat)%repul_sfc           = min(max(repulsion_surface, 0.0), 0.5)
            category_info(cat)%repul_btm           = min(max(repulsion_bottom,  0.0), 0.5)
            category_info(cat)%rmv_sfc             = remove_surface
            category_info(cat)%rmv_btm             = remove_bottom
            category_info(cat)%rmv_opn             = remove_open

            IF (dispersion_kh /= UNDEF) disph      = dispersion_kh
            IF (dispersion_kv /= UNDEF) dispv      = dispersion_kv
            category_info(cat)%disph               = disph
            category_info(cat)%dispv               = dispv
            category_info(cat)%dispgrad            = dispgrad
            category_info(cat)%fish_kinesis(:)     = fish_kinesis(:)
            IF (dispersion_pr(1) /= UNDEF) dispersion_sc(1) = dispersion_pr(1)
            IF (dispersion_pr(2) /= UNDEF) dispersion_sc(2) = dispersion_pr(2)
            IF (dispersion_pr(3) /= UNDEF) dispersion_sc(3) = dispersion_pr(3)

            !category_info%dispr represents foctor of disspersion coefficient with respect to eddy viscosity (i.e. inverse of the Schmidt number)
            IF (dispersion_sc(1) == UNDEF) THEN
               category_info(cat)%dispr(1) = category_info(0)%dispr(1)
            ELSE IF (dispersion_sc(1) == 0.0) THEN
               category_info(cat)%dispr(1) = 0.0
            ELSE
               category_info(cat)%dispr(1) = 1.0 / dispersion_sc(1)
            END IF
            IF (dispersion_sc(2) == UNDEF) THEN
               category_info(cat)%dispr(2) = category_info(cat)%dispr(1)
            ELSE IF (dispersion_sc(2) == 0.0) THEN
               category_info(cat)%dispr(2) = 0.0
            ELSE
               category_info(cat)%dispr(2) = 1.0/dispersion_sc(2)
            END IF
            IF (dispersion_sc(3) == UNDEF) THEN
               category_info(cat)%dispr(3) = category_info(0)%dispr(3)
            ELSE IF (dispersion_sc(3) == 0.0) THEN
               category_info(cat)%dispr(3) = 0.0
            ELSE
               category_info(cat)%dispr(3) = 1.0/dispersion_sc(3)
            END IF

            category_info(cat)%tkesource           = tkesource
            category_info(cat)%record_grave        = record_gravefile
            category_info(cat)%record_birth        = record_birthfile
            category_info(cat)%record_watch        = record_watchfile
            category_info(cat)%buoyancy            = buoyancy_coupling
            category_info(cat)%settling            = settling_scheme
            IF (settling_rubey1933) category_info(cat)%settling = 2
            category_info(cat)%hist_prop           = histogram_property
            category_info(cat)%hist_interval       = interval_seconds(histogram_interval)
            category_info(cat)%hist_bins           = histogram_bins(:)

            category_info(cat)%stdrift             = stokes_drift
         END DO
      END IF

      CALL bcast(n_categories)

      CALL assert(n_categories > 0, "no particle category is specified")

      DO cat=1, n_categories
         CALL bcast(category_info(cat)%name)
         CALL bcast(category_info(cat)%code)

         DO n=1, n_prop
            CALL bcast(category_info(cat)%property_name(n))
            CALL bcast(category_info(cat)%property_default(n))
         END DO

         CALL bcast(category_info(cat)%fix_x)
         CALL bcast(category_info(cat)%fix_y)
         CALL bcast(category_info(cat)%fix_z)

         CALL bcast(category_info(cat)%rk4)
         CALL bcast(category_info(cat)%interp)

         CALL bcast(category_info(cat)%lifetime)

         CALL bcast(category_info(cat)%rho)

         CALL bcast(category_info(cat)%repul)
         CALL bcast(category_info(cat)%repul_sfc)
         CALL bcast(category_info(cat)%repul_btm)

         CALL bcast(category_info(cat)%rmv_sfc)
         CALL bcast(category_info(cat)%rmv_btm)
         CALL bcast(category_info(cat)%rmv_opn)

         CALL bcast(category_info(cat)%disph)
         CALL bcast(category_info(cat)%dispv)
         CALL bcast(category_info(cat)%dispr)
         CALL bcast(category_info(cat)%dispgrad)

         CALL bcast(category_info(cat)%fish_kinesis)

         CALL bcast(category_info(cat)%tkesource)

         CALL bcast(category_info(cat)%record_grave)
         CALL bcast(category_info(cat)%record_birth)
         CALL bcast(category_info(cat)%record_watch)

         CALL bcast(category_info(cat)%buoyancy)
         CALL bcast(category_info(cat)%settling)

         CALL bcast(category_info(cat)%hist_prop)
         CALL bcast(category_info(cat)%hist_interval)
         CALL bcast(category_info(cat)%hist_bins)

         CALL bcast(category_info(cat)%stdrift)

         category_info(cat)%propindex_u = property_index(cat, 'U')
         category_info(cat)%propindex_v = property_index(cat, 'V')
         category_info(cat)%propindex_w = property_index(cat, 'W')

         category_info(cat)%propindex_diam = property_index(cat, 'DIAMETER')
         category_info(cat)%propindex_mass = property_index(cat, 'MASS')

         DO n=1, n_probe
            category_info(cat)%propindex_probe(n) = property_index(cat, 'PROBE_'//probe_vars(n))
         END DO

         m = 1
         DO n=1, n_tracer
            i = property_index(cat, 'PROBE_' // trim(tracer_name(n)))
            IF (i /= 0) THEN
               category_info(cat)%propindex_probe_tracer(1,m) = i
               category_info(cat)%propindex_probe_tracer(2,m) = n
               m = m+1
            END IF
         END DO
         category_info(cat)%propindex_probe_tracer(:,m:) = 0

         m = 1
         DO n=1, n_tracer
            i = property_index(cat, 'CUMUL_' // trim(tracer_name(n)))
            IF (i /= 0) THEN
               CALL assert(category_info(cat)%property_name(i+1) == '', "property CUMUL_* requires 2 slots")
               category_info(cat)%propindex_cumul_tracer(1,m) = i
               category_info(cat)%propindex_cumul_tracer(2,m) = n
               m = m+1
            END IF
         END DO
         category_info(cat)%propindex_cumul_tracer(:,m:) = 0

         IF (category_info(cat)%settling /= 0) THEN
            CALL assert (category_info(cat)%propindex_w/=0 .AND. category_info(cat)%propindex_diam/=0, "setting SETTLING_SCHEME requires 'W' and 'DIAMETER' properties")
            IF (sin_gravity /= 0.0) THEN
               CALL assert (category_info(cat)%propindex_u/=0, "setting SETTLING_SCHEME requires 'U' propertiy when non-zero GRAVITY_ANGLE is assigned")
            END IF
         END IF
      END DO

      IF (rank==0) THEN
         WRITE(REPORT_UNIT, *)
         WRITE(REPORT_UNIT,'(A,I3,A)') 'using ', n_categories, ' particle-categories.'
         WRITE(REPORT_UNIT, *)  ' category_code category_name property_name...'
         DO cat=1, n_categories
            WRITE(REPORT_UNIT, '(I4)',    ADVANCE='NO') category_info(cat)%code
            WRITE(REPORT_UNIT, '(X,A16)', ADVANCE='NO') category_info(cat)%name
            DO n=1, n_prop
               IF (category_info(cat)%property_name(n) == '') EXIT
               WRITE(REPORT_UNIT, '(X,A16)', ADVANCE='NO') category_info(cat)%property_name(n)
            END DO
            WRITE(REPORT_UNIT, *)
         END DO
         WRITE(REPORT_UNIT, *)
      END IF

      CALL bcast(use_convrate)

    END SUBROUTINE read_category_namelist

    SUBROUTINE read_input_namelist

      CHARACTER(128) :: filename
      CHARACTER(512) :: inputdir
      LOGICAL        :: csv
      CHARACTER(16)  :: start
      CHARACTER(16)  :: end
      CHARACTER(16)  :: interval
      INTEGER        :: calformat
      REAL(4)        :: rate(0:1023)
      CHARACTER(12)  :: mode
      INTEGER        :: periods
      LOGICAL        :: cyclic
      INTEGER(8)     :: cyclic_id
      LOGICAL        :: report
      INTEGER(8)     :: count
      INTEGER(8)     :: offset_id
      LOGICAL        :: ignore_id
      INTEGER        :: default_category
      LOGICAL        :: rate_relative, rate_propto_filesize

      CHARACTER(512) :: default_inputdir
      CHARACTER(16)  :: default_start
      CHARACTER(16)  :: default_end
      CHARACTER(16)  :: default_interval
      INTEGER        :: default_calformat

      INTEGER :: imode

      INTEGER :: n, m
      INTEGER :: iostat
      CHARACTER(256) :: iomsg

      NAMELIST / particle_input / &
           filename,       &
           inputdir,       &
           csv,            &
           rate,           &
           rate_relative,  &
           mode,           &
           periods,        &
           start,          &
           end,            &
           interval,       &
           calformat,      &
           cyclic,         &
           cyclic_id,      &
           report,         &
           offset_id,      &
           ignore_id,      &
           default_category, &
           rate_propto_filesize

      CALL default_input(default_inputdir = default_inputdir, &
                         default_start    = default_start,    &
                         default_end      = default_end,      &
                         default_interval = default_interval)

      IF (rank==0) REWIND(CONFIG_UNIT)
      DO
         filename      = ''
         csv           = .FALSE.
         inputdir      = default_inputdir
         start         = default_start
         end           = default_end
         interval      = default_interval
         mode          = 'CONST'
         periods       = 1
         rate(:)       = 0.0
         rate_relative = .FALSE.
         cyclic        = .FALSE.
         cyclic_id     = 0_8
         report        = .FALSE.
         offset_id     = 0_8
         ignore_id     = .FALSE.
         default_category = 0
         rate_propto_filesize = .FALSE.

         IF (rank==0) READ(CONFIG_UNIT, NML=particle_input, IOSTAT=iostat, IOMSG=iomsg)

         CALL bcast(iostat)
         IF (iostat < 0) EXIT

         rate_relative = (rate_relative .OR. rate_propto_filesize)

         CALL bcast(filename)
         CALL bcast(csv)
         CALL bcast(inputdir)
         CALL bcast(start)
         CALL bcast(end)
         CALL bcast(interval)
         CALL bcast(mode)
         CALL bcast(periods)
         CALL bcast(rate)
         CALL bcast(rate_relative)
         CALL bcast(cyclic)
         CALL bcast(cyclic_id)
         CALL bcast(report)
         CALL bcast(offset_id)
         CALL bcast(ignore_id)
         CALL bcast(default_category)


         CALL assert(iostat == 0,    "failed to read PARTICLE_INPUT namelist", iomsg)
         CALL assert(filename /= '', "FILENAME is not specified in PARTICLE_INPUT namelist")

         CALL assert(periods >  0 .AND. periods  <= 1024, "PERIODS in PARTICLE_INPUT is in invalid range")

         SELECT CASE (trim(mode))
         CASE ('CONST', 'CONSTANT', 'const', 'constant')
            imode = 2
         CASE ('LINEAR', 'LIN', 'linear', 'lin')
            imode = 3
         CASE ('BURST',  'burst')
            imode = 4
         CASE ('SIN', 'SINE', 'SINE-CURVE', 'sin', 'sine', 'sine-curve')
            imode = 5
         CASE ('HISTORICAL', 'HIST', 'historical', 'hist')
            imode = 10
         CASE ('TIMESTAMP', 'TS', 'timestamp', 'ts')
            imode = 20
         CASE DEFAULT
            CALL assert(.FALSE.,  "unsupported MODE in PARTICLE_INPUT")
         END SELECT

         n_input = n_input + 1
         CALL assert(n_input <= max_input, "number of PARTICLE_INPUT should be less than MAX_INPUT")

         n = n_input

         input_info(n)%filepath = path(inputdir, filename)
         input_info(n)%csv      = csv
         input_info(n)%mode     = imode
         input_info(n)%periods  = periods
         input_info(n)%start    = datetime_seconds(start)
         IF (end(1:1)=='+') THEN
            input_info(n)%end   = input_info(n)%start + interval_seconds(end(2:))
         ELSE
            input_info(n)%end   = datetime_seconds(end)
         END IF
         input_info(n)%interval = interval_seconds(interval)

         input_info(n)%rate     = max(0.0, rate)
         input_info(n)%rate_rel = rate_relative
         input_info(n)%cyclic   = cyclic
         input_info(n)%cyclic_id= cyclic_id
         input_info(n)%report   = report

         input_info(n)%offset_id     = offset_id
         input_info(n)%ignore_id     = ignore_id

         input_info(n)%def_cat       = default_category

         input_info(n)%length   = 0
         input_info(n)%residual = 0.0

         input_info(n)%initialized = .FALSE.
      END DO

    END SUBROUTINE read_input_namelist

  END SUBROUTINE init_particles

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE init_particle_input(info)
    TYPE(input_info_struct), INTENT(INOUT) :: info

    INTEGER(1), ALLOCATABLE :: readbuf(:,:)
    INTEGER(1), ALLOCATABLE :: tmpbuf(:,:)

    INTEGER :: bufsize

    INTEGER :: m
    INTEGER :: iostat

    REAL(8) :: time

    IF (info%initialized) THEN
#ifdef DEBUG
       WRITE(REPORT_UNIT,*) "particle_input is allready initialized"
#endif
       RETURN
    END IF

    bufsize = 1024
    ALLOCATE(readbuf(size_of_particle, bufsize))
    readbuf(:,:) = 0_1

    IF (rank==0) THEN
       IF (is_csvfile(info%filepath) .OR. info%csv) THEN
          OPEN(UNIT   = TMP_UNIT,         &
               FORM   = 'UNFORMATTED',    &
               ACCESS = 'DIRECT',         &
               STATUS = 'SCRATCH',        &
               ACTION = 'READWRITE',      &
               RECL   = size_of_particle)

          OPEN(UNIT   = TMP_UNIT+1,       &
               FILE   = info%filepath,    &
               FORM   = 'FORMATTED',      &
               ACCESS = 'SEQUENTIAL',     &
               STATUS = 'OLD',            &
               ACTION = 'READ',           &
               IOSTAT = iostat)

          CALL encode_particles(TMP_UNIT+1, TMP_UNIT)
          CLOSE(TMP_UNIT+1)
       ELSE
          OPEN(UNIT   = TMP_UNIT,         &
               FILE   = info%filepath,    &
               FORM   = 'UNFORMATTED',    &
               ACCESS = 'DIRECT',         &
               STATUS = 'OLD',            &
               ACTION = 'READ',           &
               RECL   = size_of_particle, &
               IOSTAT = iostat)
       END IF

       CALL assert(iostat==0, "failed to open file '"//trim(info%filepath)//"'")

       m = 0
       DO
          IF (m+1 > bufsize) THEN
             ALLOCATE(tmpbuf(size_of_particle, bufsize))
             tmpbuf(:,1:bufsize) = readbuf(:,1:bufsize)
             DEALLOCATE(readbuf)
             ALLOCATE(readbuf(size_of_particle,2*bufsize))
             readbuf(:,1:bufsize)  = tmpbuf(:,1:bufsize)
             readbuf(:,bufsize+1:) = 0_1
             DEALLOCATE(tmpbuf)
             bufsize = size(readbuf,2)
          END IF
          READ(TMP_UNIT, REC=m+1, IOSTAT=iostat) readbuf(:,m+1)
          IF (iostat /= 0) EXIT
          m = m+1
       END DO

       CLOSE(TMP_UNIT)

       CALL assert(m > 0, "particle_input file '"//trim(info%filepath)//"' is empty")
       str_format = format(m)
       WRITE(REPORT_UNIT, *) "read "//trim(str_format)//" particles from '"//trim(info%filepath)//"'"

       info%length = m
    END IF


#ifdef PARALLEL_MPI
    CALL bcast(info%length)

    IF (rank /= 0 .AND. size(readbuf,2) < info%length) THEN
       DEALLOCATE(readbuf)
       ALLOCATE(readbuf(size_of_particle, info%length))
    END IF

    CALL bcast(readbuf(:,1:info%length))
#endif

    ALLOCATE(info%particles(info%length))

    DO m=1, info%length
       CALL restore_particle(readbuf(:,m), info%particles(m), default_category=info%def_cat)
       IF (info%particles(m)%category < 0) CYCLE

       IF (info%ignore_id)            info%particles(m)%id = 0
       IF (info%particles(m)%id /= 0) info%particles(m)%id = info%particles(m)%id + info%offset_id
       IF (info%particles(m)%id > 0) max_id = max(max_id, info%particles(m)%id)

    END DO

    CALL gmax(max_id, all=.TRUE.)

    IF (allocated(readbuf)) DEALLOCATE(readbuf)

    IF (info%mode==20) THEN
       info%cyclic = .FALSE.
       info%residual = 0.0

       CALL check_timestamp(info%particles(:))

       DO m=1, info%length
          IF (info%particles(m)%category /= cat_timestamp) CYCLE

          IF (info%particles(m)%id >= max(t_current, info%start)) EXIT
       END DO
       info%count = m-1
    ELSE
       IF (info%mode==10) THEN
          time = info%lastread
       ELSE
          time = info%start
       END IF

       info%count = 0
       info%residual = 0.0
       IF (t_current > time) THEN
          DO WHILE (time < t_current)

             IF (.NOT. info%cyclic .AND. info%count >= info%length) EXIT

             IF (mod(t_current - info%start,info%interval) < dtime) info%residual = 0.0

             info%residual = info%residual + get_input_count(info, time)

             info%count    = info%count    + int(floor(info%residual))
             info%residual = info%residual - floor(info%residual)

             time  = time  + dtime
          END DO
       END IF
    END IF

#ifdef DEBUG
    IF (rank==0)  WRITE(REPORT_UNIT,*) "initial particle_input count and residual=", info%count, info%residual
#endif

    info%initialized = .TRUE.

  CONTAINS
    SUBROUTINE check_timestamp(p)
      TYPE(particle_struct), INTENT(IN) :: p(:)
      INTEGER    :: i
      INTEGER(8) :: t

      CALL assert(p(1)%category==cat_timestamp, "the 1st record for 'TIMESTAMP'-mode particle input sould be a timestamp")
      t = p(1)%id

      DO i=2, size(p)
         IF (p(i)%category /= cat_timestamp) CYCLE
         CALL assert(p(i)%id >= t, "TIMESTAMPs for particle imput is not monotonically increasing")
      END DO

    END SUBROUTINE check_timestamp

  END SUBROUTINE init_particle_input

!-----------------------------------------------------------------------------------------------------------------------


  SUBROUTINE finalize_particle_input(info)
    TYPE(input_info_struct), INTENT(INOUT) :: info

    IF (.NOT. info%initialized) RETURN

#ifdef DEBUG
    IF (info%report .AND. rank==0) WRITE(REPORT_UNIT,*) "free the particle input of '"//trim(info%filepath)//"'"
#endif

    DEALLOCATE(info%particles)

    IF (rank==0) THEN
       str_format = format(info%count)
       WRITE(REPORT_UNIT, *) "total "//trim(str_format)//" particles are inputted from '"//trim(info%filepath)//"'"
    END IF

    info%residual = 0.0
    info%count    = 0_8
    info%length   = 0
    info%initialized = .FALSE.

  END SUBROUTINE finalize_particle_input

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE increase(n)
    INTEGER, OPTIONAL, INTENT(IN) :: n

    INTEGER :: i
    INTEGER :: increase_size

    IF (PRESENT(n)) THEN
       increase_size = n
    ELSE
       increase_size = max_particles
    END IF

    CALL increase_buffer(p_id,   increase_size)

    CALL increase_buffer(p_ijk, increase_size, l=3)
    CALL increase_buffer(p_xyz, increase_size, l=3)

    CALL increase_buffer(p_category, increase_size)
    CALL increase_buffer(p_property, increase_size, l=n_prop)

    CALL increase_buffer(p_seed, increase_size)

    CALL increase_buffer(p_track_index, increase_size)

    CALL increase_buffer(p_next, increase_size)
    CALL increase_buffer(p_prev, increase_size)



!$OMP PARALLEL DO
    DO i=max_particles+1, max_particles+increase_size
       p_id(i) = 0

       p_ijk(1,i) = 0
       p_ijk(2,i) = 0
       p_ijk(3,i) = 0

       p_xyz(1,i) = 0.0
       p_xyz(2,i) = 0.0
       p_xyz(3,i) = 0.0

       p_category(i) = 0

       p_property(:,i) = UNDEF

       p_seed(i) = 0

       p_next(i) = i+1
       p_prev(i) = i-1

       p_track_index(i) = -1
    END DO

    p_next(max_particles+increase_size) = heap_ptr

    IF (heap_ptr /= NULL) p_prev(heap_ptr) = max_particles+increase_size

    heap_ptr = max_particles+1

    p_prev(heap_ptr) = NULL

    heap_cnt = heap_cnt + increase_size

    max_particles = max_particles + increase_size

  END SUBROUTINE increase

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE finalize_particles
    INTEGER :: i

#ifdef PARALLEL_MPI
    INTEGER :: ierr
#endif

    DO i=1, n_tracks
       IF (track_id(i) > 0) CALL close_record(track_record(i))
    END DO

    IF (grave_record%initialized)  CALL close_record(grave_record)
    IF (birth_record%initialized)  CALL close_record(birth_record)
    IF (watch_record%initialized)  CALL close_record(watch_record)

    DO i=1, n_input
       IF (input_info(i)%initialized) CALL finalize_particle_input(input_info(i))
    END DO

    IF (allocated(p_ijk))        DEALLOCATE(p_ijk)
    IF (allocated(p_xyz))        DEALLOCATE(p_xyz)
    IF (allocated(p_property))   DEALLOCATE(p_property)
    IF (allocated(p_category))   DEALLOCATE(p_category)
    IF (allocated(p_seed))       DEALLOCATE(p_seed)

    IF (allocated(particle_ptr)) DEALLOCATE(particle_ptr)
    IF (allocated(particle_cnt)) DEALLOCATE(particle_cnt)

#ifdef PARALLEL_MPI
    IF (allocated(sendbuf))      DEALLOCATE(recvbuf)
    IF (allocated(recvbuf))      DEALLOCATE(sendbuf)

    DEALLOCATE(ncreate)

#ifdef PROFILE
    IF (rank==0) THEN
       WRITE(REPORT_UNIT, *) "--profile of particle advection ---"
       WRITE(REPORT_UNIT, '("advection region 1 :",EN12.3,"[sec]")') profile_time('padv1')
       WRITE(REPORT_UNIT, '("advection region 2 :",EN12.3,"[sec]")') profile_time('padv2')
       WRITE(REPORT_UNIT, '("advection region 3 :",EN12.3,"[sec]")') profile_time('padv3')
!      WRITE(REPORT_UNIT, '("checkout :",EN12.3,"[sec]")') profile_time('pout')
!      WRITE(REPORT_UNIT, '("input    :",EN12.3,"[sec]")') profile_time('pinput')
!      WRITE(REPORT_UNIT,'("track    :",EN12.3,"[sec]")') profile_time('ptrack')
       WRITE(REPORT_UNIT, *) "-----------------------------------"
    END IF
#endif
#endif

  END SUBROUTINE finalize_particles

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE create_particle(src, record, ptr)
    TYPE(particle_struct), INTENT(IN) :: src
    LOGICAL, INTENT(IN),  OPTIONAL :: record
    INTEGER, INTENT(OUT), OPTIONAL :: ptr

    TYPE(particle_struct) :: tmp_particle
    LOGICAL :: record_

    INTEGER :: p
    INTEGER :: i, j, k, cat
    INTEGER :: n

    i   = src%ipos
    j   = src%jpos
    k   = src%kpos
    cat = src%category

    IF (cat < 1 .OR. cat > n_categories) THEN
       IF (present(ptr)) ptr = NULL
       RETURN
    END IF

    IF (i < 1 .OR. i > isize .OR. &
        j < 1 .OR. j > jsize .OR. &
        k < 1 .OR. k > ksize) THEN
       !WRITE(STDERR_UNIT, *) "try to create particle (ID=",src%id,") outside the domain, ignored"
       IF (present(ptr)) ptr = NULL
       RETURN
    END IF

    IF (.NOT. lmask3d(i,j,k)) THEN
       !WRITE(STDERR_UNIT, *) "try to create particle (ID=",src%id,") in topography mask, ignored"
       IF (present(ptr)) ptr = NULL
       RETURN
    END IF

    IF (heap_cnt == 0) CALL increase

    CALL assert (heap_cnt > 0, "no more particles in heap!")

    p = heap_ptr
    heap_ptr = p_next(p)
    IF (heap_ptr /= NULL) p_prev(heap_ptr) = NULL
    heap_cnt = heap_cnt - 1

    IF (particle_ptr(i,j,k,cat) /= NULL) p_prev(particle_ptr(i,j,k,cat)) = p

    p_next(p) = particle_ptr(i,j,k,cat)
    p_prev(p) = NULL

    particle_ptr(i,j,k,cat) = p
    particle_cnt(i,j,k,cat) = particle_cnt(i,j,k,cat) + 1

    IF (src%id == 0) THEN
       p_id(p) = gen_id(src%ipos, src%jpos, src%kpos)
    ELSE
       p_id(p) = src%id
    END IF

    p_ijk(1,p) = src%ipos
    p_ijk(2,p) = src%jpos
    p_ijk(3,p) = src%kpos

    p_seed(p) = src%seed
    IF (p_seed(p)==0) p_seed(p) = mkseed(r8=t_current, i8=p_id(p))

    p_xyz(1,p) = src%xpos
    p_xyz(2,p) = src%ypos
    p_xyz(3,p) = src%zpos
    IF (p_xyz(1,p)==UNDEF) p_xyz(1,p) = abs(urand(p_seed(p)))
    IF (p_xyz(2,p)==UNDEF) p_xyz(2,p) = abs(urand(p_seed(p)))
    IF (p_xyz(3,p)==UNDEF) p_xyz(3,p) = abs(urand(p_seed(p)))

    p_category(p) = src%category

    DO n=1, n_prop
       IF (src%property(n) == REAL(UNDEF, PROP_KIND)) THEN
          p_property(n,p) = category_info(cat)%property_default(n)
       ELSE
          p_property(n,p) = src%property(n)
       END IF
    END DO

    IF (category_info(cat)%settling/=0) CALL set_settling_w(p)


    IF (src%track_index < 0) THEN
       p_track_index(p) = get_track_index(p_id(p))
    ELSE
       p_track_index(p) = src%track_index
    END IF

    record_ = birth_record%initialized .AND. t_current >= birth_record%start .AND. t_current < birth_record%end
    record_ = record_ .AND. category_info(p_category(p))%record_birth
    IF (present(record)) record_ = record_ .AND. record

    IF (record_) THEN
       tmp_particle = get_particle(p, seed=n_timestep)

       CALL probe(category_info(cat), tmp_particle%ipos, tmp_particle%jpos, tmp_particle%kpos, &
                                      tmp_particle%xpos, tmp_particle%ypos, tmp_particle%zpos, tmp_particle%property)
       CALL append_record(tmp_particle, birth_record)
    END IF

    IF (present(ptr)) ptr = p

  END SUBROUTINE create_particle

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE set_settling_w(ptr)
    INTEGER, INTENT(IN) :: ptr

    INTEGER :: cat
    REAL(4) :: rho, d
    REAL(4) :: default_w

    cat = p_category(ptr)

    rho = category_info(cat)%rho
    d   = p_property(category_info(cat)%propindex_diam, ptr)

    p_property(category_info(cat)%propindex_w, ptr) = terminal_velocity(d, rho, scheme=category_info(cat)%settling)

    IF (sin_gravity /= 0.0) THEN
       p_property(category_info(cat)%propindex_u, ptr) = p_property(category_info(cat)%propindex_w, ptr)*sin_gravity
       p_property(category_info(cat)%propindex_w, ptr) = p_property(category_info(cat)%propindex_w, ptr)*cos_gravity
    END IF

    default_w = category_info(cat)%property_default(category_info(cat)%propindex_w)
    IF (default_w /= UNDEF) p_property(category_info(cat)%propindex_w, ptr) = p_property(category_info(cat)%propindex_w, ptr) + default_w

  END SUBROUTINE set_settling_w

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE delete_particle(ptr, record)
    INTEGER, INTENT(INOUT) :: ptr
    LOGICAL, INTENT(IN), OPTIONAL :: record

    LOGICAL :: record_
    INTEGER :: next, prev
    INTEGER :: i, j, k, cat

    IF (ptr==NULL) RETURN

    i   = p_ijk(1,ptr)
    j   = p_ijk(2,ptr)
    k   = p_ijk(3,ptr)
    cat = p_category(ptr)

    record_ = grave_record%initialized .AND. t_current >= grave_record%start .AND. t_current < grave_record%end
    record_ = record_ .AND. category_info(cat)%record_grave
    IF (present(record)) record_ = record_ .AND. record

#ifdef DEBUG
    CALL assert(p_id(ptr) > 0, "try to delete a particle in heap")
    CALL assert(p_ijk(1,ptr)==i .AND. &
                p_ijk(2,ptr)==j .AND. &
                p_ijk(3,ptr)==k .AND. &
                p_category(ptr)==cat, "failed to delete multiple particles")
#endif

    IF (p_track_index(ptr) > 0) CALL flush_record(track_record(p_track_index(ptr)))

    next = p_next(ptr)
    prev = p_prev(ptr)

    IF (prev == NULL) THEN
       particle_ptr(i,j,k,cat) = next
    ELSE
       p_next(prev) = next
    END IF

    IF (next /= NULL) p_prev(next) = p_prev(ptr)

    IF (heap_ptr /= NULL) p_prev(heap_ptr) = ptr

    p_next(ptr) = heap_ptr
    p_prev(ptr) = NULL

    heap_ptr = ptr
    heap_cnt = heap_cnt + 1

    IF (record_) CALL append_record(get_particle(ptr, seed=n_timestep), grave_record)

    p_id(ptr)   = 0
    p_category(ptr) = 0
    p_track_index(ptr)=-1

    ptr = next

    particle_cnt(i,j,k,cat) = particle_cnt(i,j,k,cat) - 1

#ifdef DEBUG
    IF (particle_cnt(i,j,k,cat) == 0) CALL assert(particle_ptr(i,j,k,cat) == NULL, "particle_cnt==0 but particle_ptr/=NULL")
#endif


  END SUBROUTINE delete_particle

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE delete_particles(i, j, k, cat, record)
    INTEGER, INTENT(IN) :: i, j, k, cat
    LOGICAL, INTENT(IN), OPTIONAL :: record

    INTEGER :: ptr
    INTEGER :: n
    LOGICAL :: record_

    n = particle_cnt(i, j, k, cat)
    IF (n==0) RETURN

    record_ = grave_record%initialized .AND. category_info(cat)%record_grave
    IF (present(record)) record_ = record_ .AND. record

    ptr = particle_ptr(i,j,k,cat)
    DO
       IF (record_) CALL append_record(get_particle(ptr, seed=n_timestep), grave_record)

       IF (p_track_index(ptr) > 0) CALL flush_record(track_record(p_track_index(ptr)))

       p_id(ptr) = 0
       p_category(ptr) = 0
       p_track_index(ptr)=-1

       IF (p_next(ptr) == NULL) THEN
          p_next(ptr) = heap_ptr
          IF (heap_ptr/=NULL) p_prev(heap_ptr) = ptr
          heap_ptr = particle_ptr(i,j,k,cat)
          heap_cnt = heap_cnt + n
          particle_ptr(i,j,k,cat) = NULL
          particle_cnt(i,j,k,cat) = 0
          EXIT
       END IF

       ptr = p_next(ptr)
    END DO

!    ptr = particle_ptr(i, j, k, cat)
!    DO WHILE (ptr /= NULL)
!       CALL delete_particle(ptr, record)
!    END DO

  END SUBROUTINE delete_particles

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE append_record(particle, record, flush)
    TYPE(particle_struct), INTENT(IN)    :: particle
    TYPE(record_struct),   INTENT(INOUT) :: record
    LOGICAL,     OPTIONAL, INTENT(IN)    :: flush
    LOGICAL :: flush_

    flush_ = .FALSE.
    IF (present(flush)) flush_ = flush

    IF (record%n_write == size(record%buffer, 2)) THEN
       IF (flush_) THEN
          CALL flush_record(record)
       ELSE
          CALL increase_recordbuffer(record)
       END IF
    END IF

    record%n_write = record%n_write + 1
    CALL serialize(particle, record%buffer(:,record%n_write), convert_gridpos=.TRUE., convert_categorycode=.TRUE., convert_bigendian=.TRUE.)

  END SUBROUTINE append_record

!-----------------------------------------------------------------------------------------------------------------------

  PURE FUNCTION get_timestamp(time) RESULT(ts)
    REAL(8), INTENT(IN) :: time
    TYPE(particle_struct) :: ts

    INTEGER :: year, mon, doy, dom, hour, min, sec
    REAL(8) :: rem

    CALL seconds2date(time, year, doy)
    CALL seconds2time(time, hour, min, sec, rem)
    CALL doy2dom(year, doy, mon, dom)

    ts%id = INT(time, 8) ! first record is time in seconds by INTEGER(8)
    ts%ipos = -999
    ts%jpos = -999
    ts%kpos = -999
    ts%xpos = 0.0
    ts%ypos = 0.0
    ts%zpos = 0.0

    ts%category    = cat_timestamp

    ts%property(1) = REAL(year)
    ts%property(2) = REAL(mon)
    ts%property(3) = REAL(dom)
    ts%property(4) = REAL(hour)
    ts%property(5) = REAL(min)
    ts%property(6) = REAL(sec+rem)

    ts%seed        = 0
    ts%track_index = 0

  END FUNCTION get_timestamp

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE increase_recordbuffer(record)
    TYPE(record_struct), INTENT(INOUT) :: record

    INTEGER(1), ALLOCATABLE :: tmp(:,:)

    ALLOCATE(tmp(size_of_particle, record%bufsize))
    tmp(:,:) = record%buffer(:,:)
    DEALLOCATE(record%buffer)
    ALLOCATE(record%buffer(size_of_particle, record%bufsize+chunksize))
    record%buffer(:,:record%bufsize)    = tmp(:,:)
    record%buffer(:,record%bufsize+1:)  = 0_1
    record%bufsize = record%bufsize + chunksize
    DEALLOCATE(tmp)

  END SUBROUTINE increase_recordbuffer

!-----------------------------------------------------------------------------------------------------------------------

  FUNCTION search_particle(id) RESULT(ptr)
    INTEGER(8), INTENT(IN) :: id
    INTEGER :: ptr

    INTEGER :: i

    DO ptr=1, max_particles
       IF (p_id(ptr) == id) RETURN
    END DO
    ptr = NULL

  END FUNCTION search_particle

!-----------------------------------------------------------------------------------------------------------------------

  INTEGER(8) FUNCTION gen_id(i, j, k)
    INTEGER, INTENT(IN) :: i, j, k

    ncreate(i,j,k) = ncreate(i,j,k)+1
    gen_id = 2_8**62 + (((kcoord*ksize+k-1_8)*dimy + jcoord*jsize+j-1_8)*dimx + icoord*isize+i-1_8)*(2_8**id_nshift) + ncreate(i,j,k)

  END FUNCTION gen_id

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE serialize(src, buffer, convert_gridpos, convert_categorycode, convert_bigendian)
    TYPE(particle_struct), INTENT(IN)  :: src
    INTEGER(1),            INTENT(OUT) :: buffer(size_of_particle)
    LOGICAL, OPTIONAL,     INTENT(IN)  :: convert_gridpos       ! whether or not convert local grid indices to global, default .FALSE.
    LOGICAL, OPTIONAL,     INTENT(IN)  :: convert_categorycode  ! whether or not convert the category code for output, default .FALSE.
    LOGICAL, OPTIONAL,     INTENT(IN)  :: convert_bigendian     ! whether or not convert data as bigendian for output, default .FALSE.

    INTEGER(1) :: tmp(1:8)
    INTEGER(4) :: itmp, jtmp, ktmp, cat, seed
    INTEGER :: i, j, n, offset
    LOGICAL :: flag

    buffer( 1: 8) = transfer(src%id,   0_1, 8)

    cat = src%category
    flag = .FALSE.
    IF (present(convert_categorycode)) flag = convert_categorycode
    IF (cat > 0) THEN
       IF (flag) THEN
          cat = category_info(cat)%code
       ELSE
          IF (src%track_index > 0) cat = cat + src%track_index * max_categories
       END IF
    END IF

    flag = .FALSE.
    IF (present(convert_gridpos)) flag = convert_gridpos

    IF (flag .AND. cat>0) THEN
       itmp = src%ipos+isize*icoord
       jtmp = src%jpos+jsize*jcoord
       ktmp = src%kpos+ksize*kcoord

       IF (cycle_x) itmp = mod(dimx + itmp - 1, dimx) + 1
       IF (cycle_y) jtmp = mod(dimy + jtmp - 1, dimy) + 1
       IF (cycle_z) ktmp = mod(dimz + ktmp - 1, dimz) + 1
    ELSE
       itmp = src%ipos
       jtmp = src%jpos
       ktmp = src%kpos
    END IF

    buffer( 9:12) = transfer(itmp, 0_1, 4)
    buffer(13:16) = transfer(jtmp, 0_1, 4)
    buffer(17:20) = transfer(ktmp, 0_1, 4)

    buffer(21:24) = transfer(src%xpos, 0_1, 4)
    buffer(25:28) = transfer(src%ypos, 0_1, 4)
    buffer(29:32) = transfer(src%zpos, 0_1, 4)

    buffer(33:36) = transfer(cat, 0_1, 4)

    buffer(37:36+PROP_KIND*n_prop) = transfer(src%property(1:n_prop), 0_1, PROP_KIND*n_prop)

    buffer(size_of_particle-3:size_of_particle) = transfer(src%seed, 0_1, 4)

    flag = .FALSE.
    IF (present(convert_bigendian))  flag = convert_bigendian .AND. .NOT. native_bigendian
    IF (.NOT. flag) RETURN

    tmp(1:8) = buffer(1:8)
    DO i=1, 8
       buffer(i) = tmp(9-i)
    END DO
    offset = 8

    DO j=1, 7
       tmp(1:4) = buffer(offset+1:offset+4)
       DO i=1, 4
          buffer(offset+i) = tmp(5-i)
       END DO
       offset = offset + 4
    END DO

    DO j=1, n_prop
       tmp(1:PROP_KIND) = buffer(offset+1:offset+PROP_KIND)
       DO i=1, PROP_KIND
          buffer(offset+i) = tmp(PROP_KIND+1-i)
       END DO
       offset = offset + PROP_KIND
    END DO

    tmp(1:4) = buffer(offset+1:offset+4)
    DO i=1, 4
       buffer(offset+i) = tmp(5-i)
    END DO

  END SUBROUTINE serialize

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE deserialize(buffer, tgt, convert_gridpos, convert_categorycode, convert_bigendian)
    INTEGER(1),            INTENT(IN)  :: buffer(size_of_particle)
    TYPE(particle_struct), INTENT(OUT) :: tgt
    LOGICAL, OPTIONAL,     INTENT(IN)  :: convert_gridpos      !whether or not convert global grid indices to local ones,  default .FALSE.
    LOGICAL, OPTIONAL,     INTENT(IN)  :: convert_categorycode !whether or not convert the category code for internal use, default .FALSE.
    LOGICAL, OPTIONAL,     INTENT(IN)  :: convert_bigendian    !whether or not convert big endian buffer to native endian, default .FALSE.

    INTEGER(1) :: tmp(size_of_particle)
    LOGICAL :: flag

    INTEGER :: i, j, offset

    flag = .FALSE.
    IF (present(convert_bigendian)) flag = convert_bigendian .AND. .NOT. native_bigendian

    IF (flag) THEN
       DO i=1, 8
          tmp(i)   = buffer(9-i)
       END DO
       offset = 8

       DO j=1, 7
          DO i=1, 4
             tmp(offset+i) = buffer(offset+5-i)
          END DO
          offset = offset + 4
       END DO

       DO j=1, n_prop
          DO i=1, PROP_KIND
             tmp(offset+i) = buffer(offset+PROP_KIND+1-i)
          END DO
          offset = offset + PROP_KIND
       END DO

       DO i=1, 4
          tmp(offset+i) = buffer(offset+5-i)
       END DO
    ELSE
       tmp(:) = buffer(:)
    END IF

    tgt%id = transfer(tmp(1:8), 0_8)

    tgt%ipos = transfer(tmp( 9:12), 0_4)
    tgt%jpos = transfer(tmp(13:16), 0_4)
    tgt%kpos = transfer(tmp(17:20), 0_4)

    tgt%xpos = transfer(tmp(21:24), 0.0_4)
    tgt%ypos = transfer(tmp(25:28), 0.0_4)
    tgt%zpos = transfer(tmp(29:32), 0.0_4)

    flag = .FALSE.
    IF (present(convert_gridpos)) flag = convert_gridpos

    IF (flag) THEN
       IF (cycle_x) tgt%ipos = modulo(tgt%ipos-1, dimx)+1
       IF (cycle_y) tgt%jpos = modulo(tgt%jpos-1, dimy)+1
       IF (cycle_z .OR. tgt%kpos <= 0) tgt%kpos = modulo(tgt%kpos-1, dimz)+1

       tgt%ipos = tgt%ipos - isize*icoord
       tgt%jpos = tgt%jpos - jsize*jcoord
       tgt%kpos = tgt%kpos - ksize*kcoord
    END IF

    tgt%category = transfer(tmp(33:36), 0_4)

    flag = .FALSE.
    IF (present(convert_categorycode)) flag = convert_categorycode

    IF (flag) THEN
       tgt%track_index = -1
       tgt%category = category_index_bycode(tgt%category)
    ELSE
       tgt%track_index = tgt%category / max_categories
       tgt%category = mod(tgt%category, max_categories)
    END IF

    tgt%property(1:n_prop) = transfer(tmp(37:36+PROP_KIND*n_prop), 0.0_PROP_KIND, n_prop)

    tgt%seed = transfer(tmp(size_of_particle-3:size_of_particle), 0_4)

  END SUBROUTINE deserialize

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE main
    USE velocity,   ONLY: u_new=>u, v_new=>v, w_new=>w, u_old, v_old, w_old, taux_btm, tauy_btm
    USE subgrid,    ONLY: visc
    USE tracers,    ONLY: tracer, tracer_name, n_tracer, tracer_index, tracer_index_t
    USE state,      ONLY: mld

    INTEGER :: ptr, next, prev

    INTEGER :: i, j, k, ijk, m, n, t, cat
    INTEGER :: ipos, jpos, kpos
    REAL(4) :: xpos, ypos, zpos
    REAL(4) :: lx, ly, lz, ll

    INTEGER :: nprop, nt
    INTEGER :: propindex_u, propindex_v, propindex_w
    INTEGER :: propindex_stay
    INTEGER :: propindex_mass
    INTEGER :: propindex_ages, propindex_aged
    INTEGER :: propindex_traj, propindex_traj_100km
    INTEGER :: propindex_uvexp
    INTEGER :: propindex_cx, propindex_cy, propindex_cz
    INTEGER :: propindex_zmin, propindex_zmax
    INTEGER :: propindex_passflag

    REAL(4) :: utmp, vtmp, wtmp

    INTEGER :: propindex_ml_ages, propindex_ml_aged
    INTEGER :: propindex_sfc_ages, propindex_sfc_aged

    REAL(4) :: u_(0:1,0:4,0:2) ! dim1: (W,E), dim2: (O,S,N,L,O), dim3: (OLD,NOW,NEW)
    REAL(4) :: v_(0:1,0:4,0:2) ! dim1: (S,N), dim2: (O,L,O,W,E), dim3: (OLD,NOW,NEW)
    REAL(4) :: w_(0:1,0:4,0:2) ! dim1: (L,U), dim2: (O,W,E,S,N), dim3: (OLD,NOW,NEW)

    REAL(4) :: rkdx(4), rkdy(4), rkdz(4)
    REAL(4) :: repul, repul_sfc, repul_btm
    REAL(4), PARAMETER :: eps = 1.0E-6

    REAL(4) :: disp_h(0:isize+1,0:jsize+1,0:ksize+1)
    REAL(4) :: disp_v(0:isize+1,0:jsize+1,0:ksize+1)

    REAL(4) :: disp_(0:1,1:3)
    LOGICAL :: stat_disp

    REAL(4) :: passlabel(1:isize,1:jsize+1,1:ksize)
    LOGICAL :: stat_passlabel


    REAL(4) :: u_st(0:ksize), v_st(0:ksize)
    LOGICAL :: stat_st

    REAL(8) :: age
    INTEGER :: uswitch, vswitch, wswitch


#ifdef PARALLEL_MPI
    INTEGER :: sendreq(0:5)
    INTEGER :: recvreq(0:5)

    INTEGER :: ierr
#endif
    INTEGER :: tid

    INTEGER :: moved_ptr(0:nthreads-1, 0:isize+1, 0:jsize+1, 0:ksize+1)
    INTEGER :: moved_cnt(0:nthreads-1, 0:isize+1, 0:jsize+1, 0:ksize+1)

    LOGICAL :: stat_wfactor

    REAL(4) :: wfactor(isize,jsize)

    REAL(4) :: dispfactor(0:isize+1,0:jsize+1,0:ksize+1)

    LOGICAL :: rmvmask(isize,jsize)

    LOGICAL :: stat

    INTEGER :: tracer_index_tke
    LOGICAL :: flag_tke
    REAL(4) :: gp

    REAL(4) :: conv_n2p(isize,jsize,ksize)
    REAL(4) :: conv_p2n(isize,jsize,ksize)
    REAL(4) :: conv_p2z(isize,jsize,ksize)
    REAL(4) :: conv_p2d(isize,jsize,ksize)
    REAL(4) :: conv_z2n(isize,jsize,ksize)
    REAL(4) :: conv_z2d(isize,jsize,ksize)
    REAL(4) :: conv_d2n(isize,jsize,ksize)

    tracer_index_tke = tracer_index('TKE')

    u_st(:) = 0.0
    v_st(:) = 0.0
    CALL checkin('U_STOKES', u_st, axis='Z', stat=stat_st)
    CALL checkin('V_STOKES', v_st, axis='Z', stat=stat_st, add=.TRUE.)


    cat = category_index('NPZD')
    IF (cat /= 0) THEN
!$OMP PARALLEL WORKSHARE
       conv_n2p(:,:,:) = 0.0
       conv_p2n(:,:,:) = 0.0
       conv_p2z(:,:,:) = 0.0
       conv_p2d(:,:,:) = 0.0
       conv_z2n(:,:,:) = 0.0
       conv_z2d(:,:,:) = 0.0
       conv_d2n(:,:,:) = 0.0
!$OMP END PARALLEL WORKSHARE

!$OMP PARALLEL DO SCHEDULE(dynamic, ompchunk) PRIVATE(i, j, k, ptr)
       DO ijk=0, isize*jsize*ksize-1
          k = ijk/(isize*jsize) + 1
          j = (ijk - (k-1)*(isize*jsize))/isize + 1
          i = ijk - (k-1)*(isize*jsize) - (j-1)*isize + 1

          ptr = particle_ptr(i,j,k,cat)
          DO WHILE (ptr/=NULL)
             CALL update_npzd(ptr)
             ptr = p_next(ptr)
          END DO
       END DO

       CALL checkout("CONVCOUNT_N2P", conv_n2p)
       CALL checkout("CONVCOUNT_P2N", conv_p2n)
       CALL checkout("CONVCOUNT_P2Z", conv_p2z)
       CALL checkout("CONVCOUNT_P2D", conv_p2d)
       CALL checkout("CONVCOUNT_Z2N", conv_z2n)
       CALL checkout("CONVCOUNT_Z2D", conv_z2d)
       CALL checkout("CONVCOUNT_D2N", conv_d2n)
    END IF

    DO cat=1, n_categories
       CALL checkin('WFACTOR_'//trim(category_info(cat)%name),   wfactor, stat_wfactor)
       IF (.NOT. stat_wfactor) CALL checkin('WFACTOR_PARTICLES', wfactor, stat_wfactor)
       IF (.NOT. stat_wfactor) wfactor(:,:) = 1.0

       DO nprop = n_prop, 1, -1
          IF (trim(category_info(cat)%property_name(nprop)) /= '') EXIT
       END DO

       propindex_u = category_info(cat)%propindex_u
       propindex_v = category_info(cat)%propindex_v
       propindex_w = category_info(cat)%propindex_w

       propindex_mass = category_info(cat)%propindex_mass

       propindex_stay = property_index(cat, 'STAY_SIGMA0')

       propindex_ages = property_index(cat, 'AGES')
       propindex_aged = property_index(cat, 'AGED')

       propindex_traj       = property_index(cat, 'TRAJ')
       propindex_traj_100km = property_index(cat, 'TRAJ_100KM')

       propindex_uvexp = property_index(cat, 'UVEXP')

       propindex_ml_ages = property_index(cat, 'ML_AGES')
       propindex_ml_aged = property_index(cat, 'ML_AGED')
       propindex_sfc_ages = property_index(cat, 'SFC_AGES')
       propindex_sfc_aged = property_index(cat, 'SFC_AGED')

       propindex_zmin = property_index(cat, 'ZMIN')
       propindex_zmax = property_index(cat, 'ZMAX')

       repul      = category_info(cat)%repul
       repul_sfc  = category_info(cat)%repul_sfc
       repul_btm  = category_info(cat)%repul_btm

       uswitch = 1
       vswitch = 1
       wswitch = 1

       IF (category_info(cat)%fix_x) uswitch = 0
       IF (category_info(cat)%fix_y) vswitch = 0
       IF (category_info(cat)%fix_z) wswitch = 0

       propindex_passflag = property_index(cat, 'PASSFLAG')
       IF (propindex_passflag /= 0) THEN
          CALL checkin('PASSLABEL_'//trim(category_info(cat)%name), passlabel, stat=stat_passlabel)
          IF (.NOT. stat_passlabel) CALL checkin('PASSLABEL', passlabel, stat=stat_passlabel)
       END IF

!$OMP PARALLEL WORKSHARE
       disp_h(:,:,:) = category_info(cat)%disph
       disp_v(:,:,:) = category_info(cat)%dispv
!$OMP END PARALLEL WORKSHARE

       stat_disp = category_info(cat)%disph /= 0.0 .OR. category_info(cat)%dispv /= 0.0 .OR. maxval(category_info(cat)%dispr) /= 0.0
       IF (stat_disp .AND. maxval(category_info(cat)%dispr) > 0.0) THEN
!$OMP PARALLEL DO
          DO k=0, ksize+1
          DO j=0, jsize+1
          DO i=0, isize+1
             disp_h(i,j,k) = disp_h(i,j,k) + visc(i,j,k,1)*category_info(cat)%dispr(1)
             IF (viscosity_scheme == 1) THEN
                disp_v(i,j,K) = disp_v(i,j,k) + visc(i,j,k,1)*category_info(cat)%dispr(2)
             ELSE
                disp_v(i,j,K) = disp_v(i,j,k) + visc(i,j,k,2)*category_info(cat)%dispr(2)
             END IF
             IF (category_info(cat)%dispr(3)/=0.0) disp_v(i,j,K) = disp_v(i,j,k) + visc(i,j,k,3)*category_info(cat)%dispr(3)
          END DO
          END DO
          END DO
       END IF


       CALL checkin(trim(category_info(cat)%name) // '_KH', disp_h, add=.TRUE., stat=stat_disp) !for backward compatibility
       CALL checkin(trim(category_info(cat)%name) // '_KV', disp_v, add=.TRUE., stat=stat_disp) !
       CALL checkin('PARTICLES_KH', disp_h, add=.TRUE., stat=stat_disp)                         !
       CALL checkin('PARTICLES_KV', disp_v, add=.TRUE., stat=stat_disp)                         !

       CALL checkin('DISPH_'//trim(category_info(cat)%name), disp_h, add=.TRUE., stat=stat_disp)
       CALL checkin('DISPV_'//trim(category_info(cat)%name), disp_v, add=.TRUE., stat=stat_disp)
       CALL checkin('DISPH', disp_h, add=.TRUE., stat=stat_disp)
       CALL checkin('DISPV', disp_v, add=.TRUE., stat=stat_disp)

       CALL checkin('DISPFACTOR_'//trim(category_info(cat)%name), dispfactor, stat)
       IF (.NOT. stat) CALL checkin('DISPFACTOR', dispfactor, stat)

       IF (stat) THEN
!$OMP PARALLEL DO
          DO k=0, ksize+1
          DO j=0, jsize+1
          DO i=0, isize+1
             disp_h(i,j,k) = disp_h(i,j,k) * dispfactor(i,j,k)
             disp_v(i,j,k) = disp_v(i,j,k) * dispfactor(i,j,k)
          END DO
          END DO
          END DO
       END IF

       CALL checkout('DISPH_'//trim(category_info(cat)%name), disp_h)
       CALL checkout('DISPV_'//trim(category_info(cat)%name), disp_v)

       flag_tke = (propindex_mass/=0) .AND. (propindex_w/=0) .AND. (category_info(cat)%tkesource > 0.0) .AND. (tracer_index_tke/=0)
       gp = gravity * (rho_0 - category_info(cat)%rho) / category_info(cat)%rho ! buoayncy force per unit mass, positive if particle is less-dense than water

       CALL checkin('RMVMASK_'//trim(category_info(cat)%name), rmvmask, stat=stat)
       IF (.NOT. stat) CALL checkin('RMVMASK', rmvmask, stat=stat)
       IF (stat) THEN
!$OMP PARALLEL DO SCHEDULE(dynamic, ompchunk) PRIVATE(i, j, k)
          DO ijk=0, isize*jsize*ksize-1
             k = ijk/(isize*jsize) + 1
             j = (ijk - (k-1)*(isize*jsize))/isize + 1
             i = ijk - (k-1)*(isize*jsize) - (j-1)*isize + 1

             IF (.NOT. rmvmask(i,j)) CYCLE

             IF (particle_cnt(i,j,k,cat) > 0) THEN
!$OMP CRITICAL
                CALL delete_particles(i, j, k, cat)
!$OMP END CRITICAL
             END IF
          END DO
       END IF

       IF (category_info(cat)%rmv_opn) THEN
          IF (open_w .AND. icoord==0) THEN
!$OMP PARALLEL DO SCHEDULE(dynamic, jsize) COLLAPSE(2)
             DO k=1, ksize
             DO j=1, jsize
                IF (particle_cnt(1,j,k,cat) > 0) THEN
!$OMP CRITICAL
                   CALL delete_particles(1, j, k, cat)
!$OMP END CRITICAL
                END IF
             END DO
             END DO
          END IF

          IF (open_e .AND. icoord==ipes-1) THEN
!$OMP PARALLEL DO SCHEDULE(dynamic, jsize) COLLAPSE(2)
             DO k=1, ksize
             DO j=1, jsize
                IF (particle_cnt(isize,j,k,cat) > 0) THEN
!$OMP CRITICAL
                   CALL delete_particles(isize, j, k, cat)
!$OMP END CRITICAL
                END IF
             END DO
             END DO
          END IF

          IF (open_s .AND. jcoord==0) THEN
!$OMP PARALLEL DO SCHEDULE(dynamic, isize) COLLAPSE(2)
             DO k=1, ksize
             DO i=1, isize
                IF (particle_cnt(i,1,k,cat) > 0) THEN
!$OMP CRITICAL
                   CALL delete_particles(i, 1, k, cat)
!$OMP END CRITICAL
                END IF
             END DO
             END DO
          END IF

          IF (open_n .AND. jcoord==jpes-1) THEN
!$OMP PARALLEL DO SCHEDULE(dynamic, isize) COLLAPSE(2)
             DO k=1, ksize
             DO i=1, isize
                IF (particle_cnt(i,jsize,k,cat) > 0) THEN
!$OMP CRITICAL
                   CALL delete_particles(i, jsize, k, cat)
!$OMP END CRITICAL
                END IF
             END DO
             END DO
          END IF
       END IF

       IF (category_info(cat)%rmv_sfc) THEN
          CALL checkin('SFCRMVMASK_'//trim(category_info(cat)%name), rmvmask, stat=stat)
          IF (.NOT. stat) CALL checkin('SFCRMVMASK', rmvmask, stat=stat)
          IF (.NOT. stat) rmvmask(:,:) = .TRUE.

!$OMP PARALLEL DO SCHEDULE(dynamic, ompchunk) PRIVATE(i, j, k, zpos, wtmp, ptr, next)
          DO ijk=0, isize*jsize-1
             j = ijk / isize + 1
             i = ijk - (j-1)*isize + 1

             IF (.NOT. rmvmask(i,j))      CYCLE
             IF (.NOT. surface_flag(i,j)) CYCLE
             k=surface_k(i,j)

             ptr = particle_ptr(i,j,k,cat)
             DO WHILE (ptr/=NULL)
                next = p_next(ptr)

                zpos = p_xyz(3,ptr)
                wtmp = wswitch*(0.5*(w_new(i,j,k-1)+w_old(i,j,k-1))*(1.0-zpos) + 0.5*(w_new(i,j,k)+w_old(i,j,k))*zpos)
                IF (propindex_w/=0)    wtmp = wtmp + p_property(propindex_w, ptr)*wfactor(i,j)
                IF (propindex_stay/=0) wtmp = wtmp + w_sigma(i, j, k, zpos, p_property(propindex_stay, ptr))

                IF (zpos*dz(k) + wtmp*dtime > (1.0-repul_sfc)*dz(k) - max(dz(k)-dz_ref(i,j,k),0.0)) THEN
                   p_xyz(3,ptr) = 1.0 - repul_sfc
!$OMP CRITICAL
                   CALL delete_particle(ptr)
!$OMP END CRITICAL
                END IF

                ptr = next
             END DO
          END DO
       END IF

       IF (category_info(cat)%rmv_btm) THEN
          CALL checkin('BTMRMVMASK_'//trim(category_info(cat)%name), rmvmask, stat=stat)
          IF (.NOT. stat) CALL checkin('BTMRMVMASK', rmvmask, stat=stat)
          IF (.NOT. stat) rmvmask(:,:) = .TRUE.
!$OMP PARALLEL DO SCHEDULE(dynamic, ompchunk) PRIVATE(i, j, k, zpos, wtmp, ptr, next)
          DO ijk=0, isize*jsize-1
             j = ijk / isize + 1
             i = ijk - (j-1)*isize + 1

             IF (.NOT. rmvmask(i,j))     CYCLE
             IF (.NOT. bottom_flag(i,j)) CYCLE
             k=bottom_k(i,j)

             ptr = particle_ptr(i,j,k,cat)
             DO WHILE (ptr/=NULL)
                next = p_next(ptr)

                zpos = p_xyz(3,ptr)
                wtmp = wswitch*(0.5*(w_new(i,j,k-1)+w_old(i,j,k-1))*(1.0-zpos) + 0.5*(w_new(i,j,k)+w_old(i,j,k))*zpos)
                IF (propindex_w/=0)    wtmp = wtmp + p_property(propindex_w,ptr)*wfactor(i,j)
                IF (propindex_stay/=0) wtmp = wtmp + w_sigma(i, j, k, zpos, p_property(propindex_stay,ptr))

                IF (zpos*dz(k) + wtmp*dtime < repul_btm*dz(k) + max(dz(k)-dz_ref(i,j,k),0.0)) THEN
                   p_xyz(3,ptr) = repul_btm
!$OMP CRITICAL
                   CALL delete_particle(ptr)
!$OMP END CRITICAL
                END IF

                ptr = next
             END DO
          END DO
       END IF


PROFILE_BEGIN('padv1')

       tid = 0

!$OMP PARALLEL PRIVATE(i, j, k, m, n, ijk, tid) &
!$OMP          PRIVATE(ipos, jpos, kpos)        &
!$OMP          PRIVATE(xpos, ypos, zpos)        &
!$OMP          PRIVATE(lx, ly, lz, ll)          &
!$OMP          PRIVATE(utmp, vtmp, wtmp)        &
!$OMP          PRIVATE(u_, v_, w_)              &
!$OMP          PRIVATE(disp_)                   &
!$OMP          PRIVATE(rkdx, rkdy, rkdz)        &
!$OMP          PRIVATE(ptr,  prev, next)        &
!$OMP          PRIVATE(age)

!$     tid = omp_get_thread_num()

!$OMP WORKSHARE
       moved_ptr(:,:,:,:) = NULL
       moved_cnt(:,:,:,:) = 0
!$OMP END WORKSHARE

!$OMP DO SCHEDULE(dynamic, ompchunk) REDUCTION(+:nlimit)
       DO ijk=0, isize*jsize*ksize-1
          k = ijk/(isize*jsize) + 1
          j = (ijk - (k-1)*(isize*jsize))/isize + 1
          i = ijk - (k-1)*(isize*jsize) - (j-1)*isize + 1

          IF (particle_cnt(i,j,k,cat) == 0) CYCLE

          u_(0:1,0,0) = u_old(i-1:i,j,k)
          u_(0:1,0,2) = u_new(i-1:i,j,k)

          v_(0:1,0,0) = v_old(i,j-1:j,k)
          v_(0:1,0,2) = v_new(i,j-1:j,k)

          w_(0:1,0,0) = w_old(i,j,k-1:k)
          w_(0:1,0,2) = w_new(i,j,k-1:k)

          IF (category_info(cat)%interp > 0) THEN
             u_(0:1,1,0) = u_(0:1,0,0) + (u_old(i-1:i,j-1,k)-u_(0:1,0,0))*imask3d(i,j-1,k)
             u_(0:1,2,0) = u_(0:1,0,0) + (u_old(i-1:i,j+1,k)-u_(0:1,0,0))*imask3d(i,j+1,k)
             u_(0:1,3,0) = u_(0:1,0,0) + (u_old(i-1:i,j,k-1)-u_(0:1,0,0))*imask3d(i,j,k-1)
             u_(0:1,4,0) = u_(0:1,0,0) + (u_old(i-1:i,j,k+1)-u_(0:1,0,0))*imask3d(i,j,k+1)

             u_(0:1,1,2) = u_(0:1,0,2) + (u_new(i-1:i,j-1,k)-u_(0:1,0,2))*imask3d(i,j-1,k)
             u_(0:1,2,2) = u_(0:1,0,2) + (u_new(i-1:i,j+1,k)-u_(0:1,0,2))*imask3d(i,j+1,k)
             u_(0:1,3,2) = u_(0:1,0,2) + (u_new(i-1:i,j,k-1)-u_(0:1,0,2))*imask3d(i,j,k-1)
             u_(0:1,4,2) = u_(0:1,0,2) + (u_new(i-1:i,j,k+1)-u_(0:1,0,2))*imask3d(i,j,k+1)

             v_(0:1,1,0) = v_(0:1,0,0) + (v_old(i,j-1:j,k-1)-v_(0:1,0,0))*imask3d(i,j,k-1)
             v_(0:1,2,0) = v_(0:1,0,0) + (v_old(i,j-1:j,k+1)-v_(0:1,0,0))*imask3d(i,j,k+1)
             v_(0:1,3,0) = v_(0:1,0,0) + (v_old(i-1,j-1:j,k)-v_(0:1,0,0))*imask3d(i-1,j,k)
             v_(0:1,4,0) = v_(0:1,0,0) + (v_old(i+1,j-1:j,k)-v_(0:1,0,0))*imask3d(i+1,j,k)

             v_(0:1,1,2) = v_(0:1,0,2) + (v_new(i,j-1:j,k-1)-v_(0:1,0,2))*imask3d(i,j,k-1)
             v_(0:1,2,2) = v_(0:1,0,2) + (v_new(i,j-1:j,k+1)-v_(0:1,0,2))*imask3d(i,j,k+1)
             v_(0:1,3,2) = v_(0:1,0,2) + (v_new(i-1,j-1:j,k)-v_(0:1,0,2))*imask3d(i-1,j,k)
             v_(0:1,4,2) = v_(0:1,0,2) + (v_new(i+1,j-1:j,k)-v_(0:1,0,2))*imask3d(i+1,j,k)

             w_(0:1,1,0) = w_(0:1,0,0) + (w_old(i-1,j,k-1:k)-w_(0:1,0,0))*imask3d(i-1,j,k)
             w_(0:1,2,0) = w_(0:1,0,0) + (w_old(i+1,j,k-1:k)-w_(0:1,0,0))*imask3d(i+1,j,k)
             w_(0:1,3,0) = w_(0:1,0,0) + (w_old(i,j-1,k-1:k)-w_(0:1,0,0))*imask3d(i,j-1,k)
             w_(0:1,4,0) = w_(0:1,0,0) + (w_old(i,j+1,k-1:k)-w_(0:1,0,0))*imask3d(i,j+1,k)

             w_(0:1,1,2) = w_(0:1,0,2) + (w_new(i-1,j,k-1:k)-w_(0:1,0,2))*imask3d(i-1,j,k)
             w_(0:1,2,2) = w_(0:1,0,2) + (w_new(i+1,j,k-1:k)-w_(0:1,0,2))*imask3d(i+1,j,k)
             w_(0:1,3,2) = w_(0:1,0,2) + (w_new(i,j-1,k-1:k)-w_(0:1,0,2))*imask3d(i,j-1,k)
             w_(0:1,4,2) = w_(0:1,0,2) + (w_new(i,j+1,k-1:k)-w_(0:1,0,2))*imask3d(i,j+1,k)
          END IF

          u_(:,:,1) = 0.5*(u_(:,:,0)+u_(:,:,2))
          v_(:,:,1) = 0.5*(v_(:,:,0)+v_(:,:,2))
          w_(:,:,1) = 0.5*(w_(:,:,0)+w_(:,:,2))

          IF (stat_disp) THEN
             disp_(0,1) = 0.5*(disp_h(i,j,k)+disp_h(i-1,j,k))*imask3d(i-1,j,k)
             disp_(1,1) = 0.5*(disp_h(i,j,k)+disp_h(i+1,j,k))*imask3d(i+1,j,k)
             disp_(0,2) = 0.5*(disp_h(i,j,k)+disp_h(i,j-1,k))*imask3d(i,j-1,k)
             disp_(1,2) = 0.5*(disp_h(i,j,k)+disp_h(i,j+1,k))*imask3d(i,j+1,k)
             disp_(0,3) = 0.5*(disp_v(i,j,k)+disp_v(i,j,k-1))*imask3d(i,j,k-1)
             disp_(1,3) = 0.5*(disp_v(i,j,k)+disp_v(i,j,k+1))*imask3d(i,j,k+1)
          END IF

          ptr = particle_ptr(i,j,k,cat)

          DO WHILE (ptr /= NULL)
             next =  p_next(ptr)
             prev =  p_prev(ptr)

             xpos = p_xyz(1,ptr)
             ypos = p_xyz(2,ptr)
             zpos = p_xyz(3,ptr)

             lx = xpos
             ly = ypos
             lz = zpos

             IF (propindex_traj /= 0) THEN
                ll = sqrt(lx*lx*dx(i,j)*dx(i,j) + ly*ly*dy(i,j)*dy(i,j) + lz*lz*dz(k)*dz(k))
                p_property(propindex_traj,ptr) = p_property(propindex_traj,ptr) + REAL(ll, PROP_KIND)
                IF (propindex_traj_100km /= 0) THEN
                   p_property(propindex_traj_100km,ptr) = p_property(propindex_traj_100km,ptr) &
                        + INT(p_property(propindex_traj,ptr)/1.0E5)
                   p_property(propindex_traj,ptr) = p_property(propindex_traj,ptr) &
                        - INT(p_property(propindex_traj,ptr)/1.0E5)*1.0E5
                END IF
             END IF

             IF (propindex_ages /= 0) THEN
                p_property(propindex_ages,ptr) = p_property(propindex_ages,ptr) + REAL(dtime, PROP_KIND)
                IF (propindex_aged /= 0) THEN
                   p_property(propindex_aged,ptr) = p_property(propindex_aged,ptr) &
                        + INT(p_property(propindex_ages,ptr)/86400.0)
                   p_property(propindex_ages,ptr) = p_property(propindex_ages,ptr) &
                        - INT(p_property(propindex_ages,ptr)/86400.0)*86400.0
                END IF
             ELSE IF (propindex_aged /= 0) THEN
                p_property(propindex_aged,ptr) = p_property(propindex_aged,ptr) + REAL(dtime/86400.0, PROP_KIND)
             END IF

             IF (propindex_uvexp /= 0) THEN
                p_property(propindex_uvexp,ptr) = p_property(propindex_uvexp,ptr) &
                     + REAL(dtime, PROP_KIND)*exp(-0.04*(depth(k)+(1.0-zpos)*dz(k)))
             END IF

             IF (propindex_zmin /= 0) THEN
                IF (p_property(propindex_zmin,ptr) == UNDEF) THEN
                   p_property(propindex_zmin,ptr) = ksize*kcoord + p_ijk(3,ptr)-1+p_xyz(3,ptr)
                ELSE
                   p_property(propindex_zmin,ptr) = min(p_property(propindex_zmin,ptr), ksize*kcoord + p_ijk(3,ptr)-1+p_xyz(3,ptr))
                END IF
             END IF

             IF (propindex_zmax /= 0) THEN
                IF (p_property(propindex_zmax,ptr) == UNDEF) THEN
                   p_property(propindex_zmax,ptr) = ksize*kcoord + p_ijk(3,ptr)-1+p_xyz(3,ptr)
                ELSE
                   p_property(propindex_zmax,ptr) = max(p_property(propindex_zmax,ptr), ksize*kcoord + p_ijk(3,ptr)-1+p_xyz(3,ptr))
                END IF
             END IF


             IF (propindex_ml_ages /= 0) THEN
                IF (mld(i,j) > depth(k)+(1.0-zpos)*dz(k)) THEN
                   p_property(propindex_ml_ages,ptr) = p_property(propindex_ml_ages,ptr) + REAL(dtime, PROP_KIND)
                   IF (propindex_ml_aged /= 0) THEN
                      p_property(propindex_ml_aged,ptr) = p_property(propindex_ml_aged,ptr) &
                           + INT(p_property(propindex_ml_ages,ptr)/86400.0)
                      p_property(propindex_ml_ages,ptr) = p_property(propindex_ml_ages,ptr) &
                           - INT(p_property(propindex_ml_ages,ptr)/86400.0)*86400.0
                   END IF
                END IF
             ELSE IF (propindex_ml_aged /= 0) THEN
                IF (mld(i,j) > depth(k)+(1.0-zpos)*dz(k)) THEN
                   p_property(propindex_ml_aged,ptr) = p_property(propindex_ml_aged,ptr) + REAL(dtime/86400.0, PROP_KIND)
                END IF
             END IF

             IF (propindex_sfc_ages /= 0) THEN
                IF (kcoord==kpes-1 .AND. k == ksize) THEN
                   p_property(propindex_sfc_ages,ptr) = p_property(propindex_sfc_ages,ptr) + REAL(dtime, PROP_KIND)
                   IF (propindex_sfc_aged /= 0) THEN
                      p_property(propindex_sfc_aged,ptr) = p_property(propindex_sfc_aged,ptr) &
                           + INT(p_property(propindex_sfc_ages,ptr)/86400.0)
                      p_property(propindex_sfc_ages,ptr) = p_property(propindex_sfc_ages,ptr) &
                           - INT(p_property(propindex_sfc_ages,ptr)/86400.0)*86400.0
                   END IF
                END IF
             ELSE IF (propindex_sfc_aged /= 0) THEN
                IF (kcoord==kpes-1 .AND. k == ksize) THEN
                   p_property(propindex_sfc_aged,ptr) = p_property(propindex_sfc_aged,ptr) + REAL(dtime/86400.0, PROP_KIND)
                END IF
             END IF

             IF (category_info(cat)%lifetime /= UNDEF) THEN
                age = 0.0
                IF (propindex_aged /= 0) age = age + p_property(propindex_aged,ptr)*86400.0
                IF (propindex_ages /= 0) age = age + p_property(propindex_ages,ptr)

                IF (category_info(cat)%lifetime < age) THEN
                   CALL delete_particle(ptr)
                   CYCLE
                END IF
             END IF

             CALL probe(category_info(cat), i, j, k, xpos, ypos, zpos, p_property(:,ptr))
             CALL cumul(category_info(cat), i, j, k, xpos, ypos, zpos, p_property(:,ptr))

             IF (propindex_passflag /= 0 .AND. stat_passlabel) CALL set_passflag(passlabel(i,j,k), p_property(propindex_passflag,ptr))

             utmp = 0.0
             vtmp = 0.0
             wtmp = 0.0

             IF (stat_disp) THEN
                ! effective diffusivity coefficient K = sigma^2(r) /(2*dtime) = sigma^2(v)*dtime/2 => sigma(v) = sqrt(2K/dtime)
                ! using uniform distribution random numbers of [-sqrt(3),sqrt(3)] range
                utmp = utmp + urand(p_seed(ptr)) * sqrt(6*interp1d(disp_(:,1), xpos)*idtime_r4)
                vtmp = vtmp + urand(p_seed(ptr)) * sqrt(6*interp1d(disp_(:,2), ypos)*idtime_r4)
                wtmp = wtmp + urand(p_seed(ptr)) * sqrt(6*interp1d(disp_(:,3), zpos)*idtime_r4)
                !utmp = utmp + nrand(p_seed(ptr)) * sqrt(2*interp1d(disp_(:,1), xpos)*idtime)
                !vtmp = vtmp + nrand(p_seed(ptr)) * sqrt(2*interp1d(disp_(:,2), ypos)*idtime)
                !wtmp = wtmp + nrand(p_seed(ptr)) * sqrt(2*interp1d(disp_(:,3), zpos)*idtime)
                IF (category_info(cat)%dispgrad) THEN
                   utmp = utmp + (disp_(1,1)-disp_(0,1)) * idx0(i,j)
                   vtmp = vtmp + (disp_(1,2)-disp_(0,2)) * idy0(i,j)
                   wtmp = wtmp + (disp_(1,3)-disp_(0,3)) * idz0(k)
                END IF
             END IF

             IF ((category_info(cat)%fish_kinesis(1) /= UNDEF) .AND. (mod(t_current, category_info(cat)%fish_kinesis(1)) < dtime)) THEN
                CALL update_fish_kinesis(ptr, category_info(cat)%fish_kinesis(2), category_info(cat)%fish_kinesis(3), category_info(cat)%fish_kinesis(4))
             END IF

             IF (propindex_u /= 0) utmp = utmp + p_property(propindex_u,ptr)
             IF (propindex_v /= 0) vtmp = vtmp + p_property(propindex_v,ptr)
             IF (propindex_w /= 0) wtmp = wtmp + p_property(propindex_w,ptr)*wfactor(i,j)

             IF (category_info(cat)%stdrift .AND. stat_st) THEN
                utmp = utmp +  u_st(k-1)*(1.0-zpos) + u_st(k)*zpos
                vtmp = vtmp +  v_st(k-1)*(1.0-zpos) + v_st(k)*zpos
             END IF

             IF (propindex_stay/=0) wtmp = wtmp + w_sigma(i, j, k, p_xyz(3,ptr), p_property(propindex_stay,ptr))

             IF (flag_tke) tracer(i,j,k,tracer_index_tke) = tracer(i,j,k,tracer_index_tke) &
                  + max(0.0, category_info(cat)%tkesource * wtmp*p_property(propindex_mass,ptr)*gp*dtime) / (rho_0 * dvol(i,j,k))

             IF (category_info(cat)%rk4) THEN
                SELECT CASE (category_info(cat)%interp)
                CASE (2)
                   rkdx(1) = (utmp + uswitch*interp_q(u_(:,:,0), xpos, ypos, zpos))*dtime*idx0(i,j)
                   rkdy(1) = (vtmp + vswitch*interp_q(v_(:,:,0), ypos, zpos, xpos))*dtime*idy0(i,j)
                   rkdz(1) = (wtmp + wswitch*interp_q(w_(:,:,0), zpos, xpos, ypos))*dtime*idz0(k)

                   rkdx(2) = (utmp + uswitch*interp_q(u_(:,:,1), xpos+0.5*rkdx(1), ypos+0.5*rkdy(1), zpos+0.5*rkdz(1)))*dtime*idx0(i,j)
                   rkdy(2) = (vtmp + vswitch*interp_q(v_(:,:,1), ypos+0.5*rkdy(1), zpos+0.5*rkdz(1), xpos+0.5*rkdx(1)))*dtime*idy0(i,j)
                   rkdz(2) = (wtmp + wswitch*interp_q(w_(:,:,1), zpos+0.5*rkdz(1), xpos+0.5*rkdx(1), ypos+0.5*rkdy(1)))*dtime*idz0(k)

                   rkdx(3) = (utmp + uswitch*interp_q(u_(:,:,1), xpos+0.5*rkdx(2), ypos+0.5*rkdy(2), zpos+0.5*rkdz(2)))*dtime*idx0(i,j)
                   rkdy(3) = (vtmp + vswitch*interp_q(v_(:,:,1), ypos+0.5*rkdy(2), zpos+0.5*rkdz(2), xpos+0.5*rkdx(2)))*dtime*idy0(i,j)
                   rkdz(3) = (wtmp + wswitch*interp_q(w_(:,:,1), zpos+0.5*rkdz(2), xpos+0.5*rkdx(2), ypos+0.5*rkdy(2)))*dtime*idz0(k)

                   rkdx(4) = (utmp + uswitch*interp_q(u_(:,:,2), xpos+rkdx(3), ypos+rkdy(3), zpos+rkdz(3)))*dtime*idx0(i,j)
                   rkdy(4) = (vtmp + vswitch*interp_q(v_(:,:,2), ypos+rkdy(3), zpos+rkdz(3), xpos+rkdx(3)))*dtime*idy0(i,j)
                   rkdz(4) = (wtmp + wswitch*interp_q(w_(:,:,2), zpos+rkdz(3), xpos+rkdx(3), ypos+rkdy(3)))*dtime*idz0(k)
                CASE (1)
                   rkdx(1) = (utmp + uswitch*interp_l(u_(:,:,0), xpos, ypos, zpos))*dtime*idx0(i,j)
                   rkdy(1) = (vtmp + vswitch*interp_l(v_(:,:,0), ypos, zpos, xpos))*dtime*idy0(i,j)
                   rkdz(1) = (wtmp + wswitch*interp_l(w_(:,:,0), zpos, xpos, ypos))*dtime*idz0(k)

                   rkdx(2) = (utmp + uswitch*interp_l(u_(:,:,1), xpos+0.5*rkdx(1), ypos+0.5*rkdy(1), zpos+0.5*rkdz(1)))*dtime*idx0(i,j)
                   rkdy(2) = (vtmp + vswitch*interp_l(v_(:,:,1), ypos+0.5*rkdy(1), zpos+0.5*rkdz(1), xpos+0.5*rkdx(1)))*dtime*idy0(i,j)
                   rkdz(2) = (wtmp + wswitch*interp_l(w_(:,:,1), zpos+0.5*rkdz(1), xpos+0.5*rkdx(1), ypos+0.5*rkdy(1)))*dtime*idz0(k)

                   rkdx(3) = (utmp + uswitch*interp_l(u_(:,:,1), xpos+0.5*rkdx(2), ypos+0.5*rkdy(2), zpos+0.5*rkdz(2)))*dtime*idx0(i,j)
                   rkdy(3) = (vtmp + vswitch*interp_l(v_(:,:,1), ypos+0.5*rkdy(2), zpos+0.5*rkdz(2), xpos+0.5*rkdx(2)))*dtime*idy0(i,j)
                   rkdz(3) = (wtmp + wswitch*interp_l(w_(:,:,1), zpos+0.5*rkdz(2), xpos+0.5*rkdx(2), ypos+0.5*rkdy(2)))*dtime*idz0(k)

                   rkdx(4) = (utmp + uswitch*interp_l(u_(:,:,2), xpos+rkdx(3), ypos+rkdy(3), zpos+rkdz(3)))*dtime*idx0(i,j)
                   rkdy(4) = (vtmp + vswitch*interp_l(v_(:,:,2), ypos+rkdy(3), zpos+rkdz(3), xpos+rkdx(3)))*dtime*idy0(i,j)
                   rkdz(4) = (wtmp + wswitch*interp_l(w_(:,:,2), zpos+rkdz(3), xpos+rkdx(3), ypos+rkdy(3)))*dtime*idz0(k)
                CASE (0)
                   rkdx(1) = (utmp + uswitch*interp1d(u_(:,0,0), xpos))*dtime*idx0(i,j)
                   rkdy(1) = (vtmp + vswitch*interp1d(v_(:,0,0), ypos))*dtime*idy0(i,j)
                   rkdz(1) = (wtmp + wswitch*interp1d(w_(:,0,0), zpos))*dtime*idz0(k)

                   rkdx(2) = (utmp + uswitch*interp1d(u_(:,0,1), xpos+0.5*rkdx(1)))*dtime*idx0(i,j)
                   rkdy(2) = (vtmp + vswitch*interp1d(v_(:,0,1), ypos+0.5*rkdy(1)))*dtime*idy0(i,j)
                   rkdz(2) = (wtmp + wswitch*interp1d(w_(:,0,1), zpos+0.5*rkdz(1)))*dtime*idz0(k)

                   rkdx(3) = (utmp + uswitch*interp1d(u_(:,0,1), xpos+0.5*rkdx(2)))*dtime*idx0(i,j)
                   rkdy(3) = (vtmp + vswitch*interp1d(v_(:,0,1), ypos+0.5*rkdy(2)))*dtime*idy0(i,j)
                   rkdz(3) = (wtmp + wswitch*interp1d(w_(:,0,1), zpos+0.5*rkdz(2)))*dtime*idz0(k)

                   rkdx(4) = (utmp + uswitch*interp1d(u_(:,0,2), xpos+rkdx(3)))*dtime*idx0(i,j)
                   rkdy(4) = (vtmp + vswitch*interp1d(v_(:,0,2), ypos+rkdy(3)))*dtime*idy0(i,j)
                   rkdz(4) = (wtmp + wswitch*interp1d(w_(:,0,2), zpos+rkdz(3)))*dtime*idz0(k)
                END SELECT

                xpos = xpos + (rkdx(1)+2*rkdx(2)+2*rkdx(3)+rkdx(4))/6.0
                ypos = ypos + (rkdy(1)+2*rkdy(2)+2*rkdy(3)+rkdy(4))/6.0
                zpos = zpos + (rkdz(1)+2*rkdz(2)+2*rkdz(3)+rkdz(4))/6.0
             ELSE
                SELECT CASE (category_info(cat)%interp)
                CASE (2)
                   xpos = xpos + (utmp + uswitch*interp_q(u_(:,:,1), xpos, ypos, zpos))*dtime*idx0(i,j)
                   ypos = ypos + (vtmp + vswitch*interp_q(v_(:,:,1), ypos, zpos, xpos))*dtime*idy0(i,j)
                   zpos = zpos + (wtmp + wswitch*interp_q(w_(:,:,1), zpos, xpos, ypos))*dtime*idz0(k)
                CASE (1)
                   xpos = xpos + (utmp + uswitch*interp_l(u_(:,:,1), xpos, ypos, zpos))*dtime*idx0(i,j)
                   ypos = ypos + (vtmp + vswitch*interp_l(v_(:,:,1), ypos, zpos, xpos))*dtime*idy0(i,j)
                   zpos = zpos + (wtmp + wswitch*interp_l(w_(:,:,1), zpos, xpos, ypos))*dtime*idz0(k)
                CASE (0)
                   xpos = xpos + (utmp + uswitch*interp1d(u_(:,0,1), xpos))*dtime*idx0(i,j)
                   ypos = ypos + (vtmp + vswitch*interp1d(v_(:,0,1), ypos))*dtime*idy0(i,j)
                   zpos = zpos + (wtmp + wswitch*interp1d(w_(:,0,1), zpos))*dtime*idz0(k)
                END SELECT
             END IF

             lx = abs(lx - xpos)
             ly = abs(ly - ypos)
             lz = abs(lz - zpos)

             ipos = i
             jpos = j
             kpos = k

             IF (lx >= ly .AND. lx >= lz) THEN
                CALL cross_x(ipos, jpos, kpos, xpos)

                IF (ly >= lz) THEN
                   CALL cross_y(ipos, jpos, kpos, ypos)
                   CALL cross_z(ipos, jpos, kpos, zpos)
                ELSE
                   CALL cross_z(ipos, jpos, kpos, zpos)
                   CALL cross_y(ipos, jpos, kpos, ypos)
                END IF
             ELSE IF (ly >= lz .AND. ly >= lx) THEN
                CALL cross_y(ipos, jpos, kpos, ypos)
                IF (lz >= lx) THEN
                   CALL cross_z(ipos, jpos, kpos, zpos)
                   CALL cross_x(ipos, jpos, kpos, xpos)
                ELSE
                   CALL cross_x(ipos, jpos, kpos, xpos)
                   CALL cross_z(ipos, jpos, kpos, zpos)
                END IF
             ELSE
                CALL cross_z(ipos, jpos, kpos, zpos)
                IF (lx >= ly) THEN
                   CALL cross_x(ipos, jpos, kpos, xpos)
                   CALL cross_y(ipos, jpos, kpos, ypos)
                ELSE
                   CALL cross_y(ipos, jpos, kpos, ypos)
                   CALL cross_x(ipos, jpos, kpos, xpos)
                END IF
             END IF

             p_xyz(1,ptr) = xpos
             p_xyz(2,ptr) = ypos
             p_xyz(3,ptr) = zpos

             IF (i/= ipos .OR. j/=jpos .OR. k/=kpos) THEN
                p_ijk(1,ptr) = ipos
                p_ijk(2,ptr) = jpos
                p_ijk(3,ptr) = kpos

                IF ((ipos==i+1 .AND. xpos>=1.0-eps) .OR. (ipos==i-1 .AND. xpos<=eps)) nlimit(1) = nlimit(1) + 1
                IF ((jpos==j+1 .AND. ypos>=1.0-eps) .OR. (jpos==j-1 .AND. ypos<=eps)) nlimit(2) = nlimit(2) + 1
                IF ((kpos==k+1 .AND. zpos>=1.0-eps) .OR. (kpos==k-1 .AND. zpos<=eps)) nlimit(3) = nlimit(3) + 1

                IF (prev == NULL) THEN
                   particle_ptr(i,j,k,cat) = next
                ELSE
                   p_next(prev) = next
                END IF
                IF (next /= NULL) p_prev(next) = prev

                particle_cnt(i,j,k,cat) = particle_cnt(i,j,k,cat) - 1

                IF (moved_ptr(tid,ipos,jpos,kpos) == NULL) THEN
                   p_next(ptr) = ptr
                   p_prev(ptr) = ptr
                   moved_ptr(tid,ipos,jpos,kpos) = ptr
                ELSE
!REQUIRE RENEWED IMPLEMENTATION!!!! 2020/11/16
                   p_next(ptr) = p_next(moved_ptr(tid,ipos,jpos,kpos))
                   p_prev(ptr) = moved_ptr(tid,ipos,jpos,kpos)
                   p_next(moved_ptr(tid,ipos,jpos,kpos)) = ptr
                   p_prev(p_next(ptr)) = ptr
                END IF

                moved_cnt(tid,ipos,jpos,kpos) =  moved_cnt(tid,ipos,jpos,kpos) + 1
             END IF

             IF (p_track_index(ptr) > 0) CALL track_particle(ptr)

             ptr = next
          END DO
       END DO
!$OMP END PARALLEL

#if defined(PARALLEL_MPI) && defined(PROFILE)
       CALL mpi_barrier(comm, ierr)
#endif

PROFILE_END('padv1')

PROFILE_BEGIN('padv2')
!$OMP PARALLEL DO SCHEDULE(dynamic, ompchunk*4) PRIVATE(i, j, k, t)
       DO ijk=0, (isize+2)*(jsize+2)*(ksize+2)-1
          k = ijk / ((isize+2)*(jsize+2))
          j = (ijk - k*(isize+2)*(jsize+2))/(isize+2)
          i = ijk - k*(isize+2)*(jsize+2) - j*(isize+2)

          DO t=0, nthreads-1
             IF (moved_ptr(t,i,j,k) == NULL) CYCLE

             IF (particle_ptr(i,j,k,cat) /= NULL) p_prev(particle_ptr(i,j,k,cat)) = p_prev(moved_ptr(t,i,j,k))

             p_next(p_prev(moved_ptr(t,i,j,k))) = particle_ptr(i,j,k,cat)
             particle_ptr(i,j,k,cat) = moved_ptr(t,i,j,k)
             particle_cnt(i,j,k,cat) = particle_cnt(i,j,k,cat) + moved_cnt(t,i,j,k)
             p_prev(particle_ptr(i,j,k,cat)) = NULL
          END DO
       END DO
PROFILE_END('padv2')

!-----
       IF (cycle_x) THEN
          propindex_cx = property_index(cat, 'CYCLECOUNT_X')
          IF (propindex_cx /= 0) THEN
             IF (icoord==ipes-1) THEN
!$OMP PARALLEL DO SCHEDULE(dynamic, jsize+2) COLLAPSE(2) PRIVATE(ptr)
                DO k=0, ksize+1
                DO j=0, jsize+1
                   ptr = particle_ptr(isize+1,j,k,cat)
                   DO WHILE (ptr/=NULL)
                      p_property(propindex_cx,ptr) = p_property(propindex_cx,ptr) + 1.0
                      ptr = p_next(ptr)
                   END DO
                END DO
                END DO
             END IF
             IF (icoord==0) THEN
!$OMP PARALLEL DO SCHEDULE(dynamic, jsize+2) COLLAPSE(2) PRIVATE(ptr)
                DO k=0, ksize+1
                DO j=0, jsize+1
                   ptr = particle_ptr(0,j,k,cat)
                   DO WHILE (ptr/=NULL)
                      p_property(propindex_cx,ptr) = p_property(propindex_cx,ptr) - 1.0
                      ptr = p_next(ptr)
                   END DO
                END DO
                END DO
             END IF
          END IF
       END IF

       IF (cycle_y) THEN
          propindex_cy = property_index(cat, 'CYCLECOUNT_Y')
          IF (propindex_cy /= 0) THEN
             IF (jcoord==jpes-1) THEN
!$OMP PARALLEL DO SCHEDULE(dynamic, isize+2) COLLAPSE(2) PRIVATE(ptr)
                DO k=0, ksize+1
                DO i=0, isize+1
                   ptr = particle_ptr(i,jsize+1,k,cat)
                   DO WHILE (ptr/=NULL)
                      p_property(propindex_cy,ptr) = p_property(propindex_cy,ptr) + 1.0
                      ptr = p_next(ptr)
                   END DO
                END DO
                END DO
             END IF
             IF (jcoord==0) THEN
!$OMP PARALLEL DO SCHEDULE(dynamic, isize+2) COLLAPSE(2) PRIVATE(ptr)
                DO k=0, ksize+1
                DO i=0, isize+1
                   ptr = particle_ptr(i,0,k,cat)
                   DO WHILE (ptr/=NULL)
                      p_property(propindex_cy,ptr) = p_property(propindex_cy,ptr) - 1.0
                      ptr = p_next(ptr)
                   END DO
                END DO
                END DO
             END IF
          END IF
       END IF

       IF (cycle_z) THEN
          propindex_cz = property_index(cat, 'CYCLECOUNT_Z')
          IF (propindex_cz /= 0) THEN
             IF (kcoord==kpes-1) THEN
!$OMP PARALLEL DO SCHEDULE(dynamic, isize+2) COLLAPSE(2) PRIVATE(ptr)
                DO j=0, jsize+1
                DO i=0, isize+1
                   ptr = particle_ptr(i,j,ksize+1,cat)
                   DO WHILE (ptr/=NULL)
                      p_property(propindex_cz,ptr) = p_property(propindex_cz,ptr) + 1.0
                      ptr = p_next(ptr)
                   END DO
                END DO
                END DO
             END IF
             IF (kcoord==0) THEN
!$OMP PARALLEL DO SCHEDULE(dynamic, isize+2) COLLAPSE(2) PRIVATE(ptr)
                DO j=0, jsize+1
                DO i=0, isize+1
                   ptr = particle_ptr(i,j,0,cat)
                   DO WHILE (ptr/=NULL)
                      p_property(propindex_cz,ptr) = p_property(propindex_cz,ptr) - 1.0
                      ptr = p_next(ptr)
                   END DO
                END DO
                END DO
             END IF
          END IF
       END IF

       IF (tripolar) THEN
          propindex_cy = property_index(cat, 'CYCLECOUNT_Y')
          IF (propindex_cy /= 0) THEN
             IF (jcoord==jpes-1 .AND. icoord<ipes/2) THEN
!$OMP PARALLEL DO SCHEDULE(dynamic, isize+2) COLLAPSE(2) PRIVATE(ptr)
                DO k=0, ksize+1
                DO i=0, isize+1
                   ptr = particle_ptr(i,jsize+1,k,cat)
                   DO WHILE (ptr/=NULL)
                      p_property(propindex_cy,ptr) = p_property(propindex_cy,ptr) + 1.0
                      ptr = p_next(ptr)
                   END DO
                END DO
                END DO
             END IF
             IF (jcoord==jpes-1 .AND. icoord>=ipes/2) THEN
!$OMP PARALLEL DO SCHEDULE(dynamic, isize+2) COLLAPSE(2) PRIVATE(ptr)
                DO k=0, ksize+1
                DO i=0, isize+1
                   ptr = particle_ptr(i,jsize+1,k,cat)
                   DO WHILE (ptr/=NULL)
                      p_property(propindex_cy,ptr) = p_property(propindex_cy,ptr) - 1.0
                      ptr = p_next(ptr)
                   END DO
                END DO
                END DO
             END IF
          END IF
       END IF

    END DO

  CONTAINS
    REAL(4) PURE FUNCTION interp1d(a, x)
      REAL(4), INTENT(IN) :: a(0:1)
      REAL(4), INTENT(IN) :: x

      interp1d = a(0)*(1.0-x) + a(1)*x
    END FUNCTION interp1d

    REAL(4) PURE FUNCTION interp2d(a, x, y)
      REAL(4), INTENT(IN) :: a(0:1, 0:1)
      REAL(4), INTENT(IN) :: x, y

      interp2d = (a(0,0)*(1.0-x) + a(1,0)*x)*(1.0-y) + (a(0,1)*(1.0-x) + a(1,1)*x)*y
    END FUNCTION interp2d

    REAL(4) PURE FUNCTION interp3d(a, x, y, z)
      REAL(4), INTENT(IN) :: a(0:1, 0:1, 0:1)
      REAL(4), INTENT(IN) :: x, y, z

      interp3d = ((a(0,0,0)*(1.0-x) + a(1,0,0)*x)*(1.0-y) + (a(0,1,0)*(1.0-x) + a(1,1,0)*x)*y)*(1.0-z) &
              +  ((a(0,0,1)*(1.0-x) + a(1,0,1)*x)*(1.0-y) + (a(0,1,1)*(1.0-x) + a(1,1,1)*x)*y)*z
    END FUNCTION interp3d

    REAL(4) PURE FUNCTION interp_l(u, x, y, z) ! bi-linear
      REAL(4), INTENT(IN) :: u(0:1,0:4)
      REAL(4), INTENT(IN) :: x, y, z

      IF (z < 0.5) THEN
         IF (y < 0.5) THEN
            interp_l = (u(0,0) + (u(0,1)-u(0,0))*(0.5-y) + (u(0,3)-u(0,0))*(0.5-z))*(1.0-x) &
                     + (u(1,0) + (u(1,1)-u(1,0))*(0.5-y) + (u(1,3)-u(1,0))*(0.5-z))*x
         ELSE
            interp_l = (u(0,0) + (u(0,2)-u(0,0))*(y-0.5) + (u(0,3)-u(0,0))*(0.5-z))*(1.0-x) &
                     + (u(1,0) + (u(1,2)-u(1,0))*(y-0.5) + (u(1,3)-u(1,0))*(0.5-z))*x
         END IF
      ELSE
         IF (y < 0.5) THEN
            interp_l = (u(0,0) + (u(0,1)-u(0,0))*(0.5-y) + (u(0,4)-u(0,0))*(z-0.5))*(1.0-x) &
                     + (u(1,0) + (u(1,1)-u(1,0))*(0.5-y) + (u(1,4)-u(1,0))*(z-0.5))*x
         ELSE
            interp_l = (u(0,0) + (u(0,2)-u(0,0))*(y-0.5) + (u(0,4)-u(0,0))*(z-0.5))*(1.0-x) &
                     + (u(1,0) + (u(1,2)-u(1,0))*(y-0.5) + (u(1,4)-u(1,0))*(z-0.5))*x
         END IF
      END IF

    END FUNCTION interp_l

    REAL(4) PURE FUNCTION interp_q(u, x, y, z) ! bi-quadratic for y-z plane, linear for x-axis
      REAL(4), INTENT(IN) :: u(0:1,0:4)
      REAL(4), INTENT(IN) :: x, y, z

      interp_q = (u(0,0) + ((0.5*(u(0,1)+u(0,2))-u(0,0))*(y-0.5) + 0.5*(u(0,2)-u(0,1)))*(y-0.5) &
                         + ((0.5*(u(0,3)+u(0,4))-u(0,0))*(z-0.5) + 0.5*(u(0,4)-u(0,3)))*(z-0.5)) * (1.0-x) &
               + (u(1,0) + ((0.5*(u(1,1)+u(1,2))-u(1,0))*(y-0.5) + 0.5*(u(1,2)-u(1,1)))*(y-0.5) &
                         + ((0.5*(u(1,3)+u(1,4))-u(1,0))*(z-0.5) + 0.5*(u(1,4)-u(1,3)))*(z-0.5)) * x
    END FUNCTION interp_q

!------------------------------

    SUBROUTINE cross_x(ipos, jpos, kpos, xpos)
      INTEGER, INTENT(INOUT) :: ipos
      INTEGER, INTENT(IN)    :: jpos
      INTEGER, INTENT(IN)    :: kpos
      REAL(4), INTENT(INOUT) :: xpos

      IF (xpos >= 1.0 .AND. lmask3d(ipos+1,jpos,kpos)) THEN
         xpos = min((xpos-1.0)*dx(ipos,jpos)*idx0(ipos+1,jpos), 1.0-eps)
         ipos = ipos+1
      ELSE IF (xpos <= 0.0 .AND. lmask3d(ipos-1,jpos,kpos)) THEN
         xpos = max(1.0 + xpos*dx(ipos,jpos)*idx0(ipos-1,jpos), eps)
         ipos = ipos-1
      END IF

      IF (.NOT. lmask3d(ipos-1,jpos,kpos)) xpos = max(xpos, repul)
      IF (.NOT. lmask3d(ipos+1,jpos,kpos)) xpos = min(xpos, 1.0-repul)

    END SUBROUTINE cross_x

!------------------------------

    SUBROUTINE cross_y(ipos, jpos, kpos, ypos)
      INTEGER, INTENT(IN)    :: ipos
      INTEGER, INTENT(INOUT) :: jpos
      INTEGER, INTENT(IN)    :: kpos
      REAL(4), INTENT(INOUT) :: ypos

      IF (ypos >= 1.0 .AND. lmask3d(ipos,jpos+1,kpos)) THEN
         ypos = min((ypos-1.0)*dy(ipos,jpos)*idy0(ipos,jpos+1), 1.0-eps)
         jpos = jpos+1
      ELSE IF (ypos <= 0.0 .AND. lmask3d(ipos,jpos-1,kpos)) THEN
         ypos = max(1.0 + ypos*dy(ipos,jpos)*idy0(ipos-1,jpos), eps)
         jpos = jpos-1
      END IF

      IF (.NOT. lmask3d(ipos,jpos-1,kpos)) ypos = max(ypos, repul)
      IF (.NOT. lmask3d(ipos,jpos+1,kpos)) ypos = min(ypos, 1.0-repul)

    END SUBROUTINE cross_y

!------------------------------

    SUBROUTINE cross_z(ipos, jpos, kpos, zpos)
      INTEGER, INTENT(IN)    :: ipos
      INTEGER, INTENT(IN)    :: jpos
      INTEGER, INTENT(INOUT) :: kpos
      REAL(4), INTENT(INOUT) :: zpos

      IF (zpos >= 1.0 .AND. lmask3d(ipos,jpos,kpos+1)) THEN
         zpos = min((zpos-1.0)*dz(kpos)*idz0(kpos+1), 1.0-eps)
         kpos = kpos+1
      ELSE IF (zpos <= 0.0 .AND. lmask3d(ipos,jpos,kpos-1)) THEN
         zpos = max(1.0 + zpos*dz(kpos)*idz0(kpos-1), eps)
         kpos = kpos-1
      END IF

      IF (.NOT. lmask3d(ipos,jpos,kpos-1)) zpos = max(zpos, (dz(kpos) - dz_ref(ipos,jpos,kpos))*idz0(kpos) + repul_btm, 0.0)
      IF (.NOT. lmask3d(ipos,jpos,kpos+1)) zpos = min(zpos, 1.0-repul_sfc)

    END SUBROUTINE cross_z

!------------------------------

    SUBROUTINE track_particle(ptr)
      INTEGER, INTENT(IN) :: ptr
      INTEGER(8) :: tmp

      tmp = p_id(ptr)
      p_id(ptr) = INT(t_current + dtime, 8)
      CALL append_record(get_particle(ptr), track_record(p_track_index(ptr)), flush=.TRUE.)
      p_id(ptr) = tmp
    END SUBROUTINE track_particle

!--------------------------

    REAL(4) PURE FUNCTION w_sigma(i, j, k, zpos, tgt_sigma0)
      USE state, ONLY: sigma0
      INTEGER, INTENT(IN) :: i, j, k
      REAL(4), INTENT(IN) :: zpos
      REAL(4), INTENT(IN) :: tgt_sigma0

      REAL(4) :: cur_sigma0
      REAL(4) :: dsigmadz

      dsigmadz = (sigma0(i,j,k+1) - sigma0(i,j,k-1)) * idz2(k)
      dsigmadz = min(dsigmadz, -1.0E-6)

      cur_sigma0 = sigma0(i,j,k) + dsigmadz*(zpos - 0.5)*dz(k)

      w_sigma = REAL(min(max((tgt_sigma0 - cur_sigma0) / dsigmadz, -dz(k)), dz(k)) * idtime, 4)

    END FUNCTION w_sigma

!--------------------------

    SUBROUTINE update_npzd(ptr)
      USE npzd
      INTEGER, INTENT(IN) :: ptr
      INTEGER :: i, j, k
      REAL(4) :: state, count, w
      INTEGER :: tmp, seed
      REAL(4) :: rand

      i = p_ijk(1,ptr)
      j = p_ijk(2,ptr)
      k = p_ijk(3,ptr)

      state = p_property(1, ptr)
      count = p_property(2, ptr)
      w     = p_property(3, ptr)

      seed = xorshift(p_seed(ptr))
      tmp = abs(seed)
      seed = xorshift(seed)
      tmp = max(tmp, abs(seed))
      seed = xorshift(seed)
      tmp = max(tmp, abs(seed))
      seed = xorshift(seed)
      tmp = max(tmp, abs(seed))
      p_seed(ptr) = seed
      rand = tmp * 0.5E0**31

      SELECT CASE(INT(state))
      CASE (0) ! N
         IF (rand < sqrt(sqrt(convrate_n2p(i,j,k)))) THEN
            state = 1.0
            count = count + 1.0
            conv_n2p(i,j,k) = conv_n2p(i,j,k)  + 1
         END IF
      CASE (1) ! P
         IF (rand < sqrt(sqrt(convrate_p2n(i,j,k)))) THEN
            state = 0.0
            conv_p2n(i,j,k) = conv_p2n(i,j,k)  + 1
         ELSE IF (rand < sqrt(sqrt(convrate_p2n(i,j,k)+convrate_p2z(i,j,k)))) THEN
            state = 2.0
            conv_p2z(i,j,k) = conv_p2z(i,j,k)  + 1
         ELSE IF (rand < sqrt(sqrt(convrate_p2n(i,j,k)+convrate_p2z(i,j,k)+convrate_p2d(i,j,k)))) THEN
            state = 3.0
            w = -PON_SVn*(1.0 + npzd_svsigma*nrand(p_seed(ptr)))
            conv_p2d(i,j,k) = conv_p2d(i,j,k)  + 1
         END IF
      CASE (2) ! Z
         IF (rand < sqrt(sqrt(convrate_z2n(i,j,k)))) THEN
            state = 0.0
            conv_z2n(i,j,k) = conv_z2n(i,j,k)  + 1
         ELSE IF (rand < sqrt(sqrt(convrate_z2n(i,j,k)+convrate_z2d(i,j,k)))) THEN
            state = 3.0
            w = -PON_SVn*(1.0 + npzd_svsigma*nrand(p_seed(ptr)))
            conv_z2d(i,j,k) = conv_z2d(i,j,k)  + 1
         END IF
      CASE (3) ! D
         IF (rand < sqrt(sqrt(convrate_d2n(i,j,k)))) THEN
            state = 0.0
            w = 0.0
            p_xyz(3, ptr) = abs(urand(p_seed(ptr)))
            conv_d2n(i,j,k) = conv_d2n(i,j,k)  + 1
         END IF
      END SELECT

      p_property(1, ptr) = state
      p_property(2, ptr) = count
      p_property(3, ptr) = w
    END SUBROUTINE update_npzd

    SUBROUTINE update_fish_kinesis(ptr, t_lower, t_upper, swim)
      ! Kinesis argorithm for fish migration (Humston et al. 2000, fish. orceanogr.)
      INTEGER, INTENT(IN) :: ptr
      REAL(4), INTENT(IN) :: t_lower, t_upper, swim
      REAL(4) :: t, hi, rad
      REAL(4), PARAMETER :: h1 = 0.75, h2 = 0.90

      t = tracer(p_ijk(1,ptr), p_ijk(2,ptr), p_ijk(3,ptr), tracer_index_t) !in-situ t. but substituted by pot.t.

      !stepwise habitat index
      !IF ((t > t_lower .OR. t_lower==UNDEF) .AND. (t < t_upper .OR. t_upper==UNDEF)) THEN
      !   hi = 1.0
      !ELSE
      !   hi = 0.0
      !END IF

      !Gaussian-shape habitat index with sigma=(t_upper-t_lower)/2
      hi = exp(-0.5*((t-(t_lower+t_upper)/2)/((t_upper-t_lower)/2))**2)

      rad = urand(p_seed(ptr)) * 2*PI

      p_property(propindex_u, ptr) = p_property(propindex_u, ptr) * h1 * hi + swim * cos(rad) * (1.0 - h2*hi)
      p_property(propindex_v, ptr) = p_property(propindex_v, ptr) * h1 * hi + swim * sin(rad) * (1.0 - h2*hi)

    END SUBROUTINE update_fish_kinesis

!--------------------------

  END SUBROUTINE main

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE probe(info, i, j, k, x, y, z, prop)
    USE velocity,   ONLY: u, v, w, u_old, v_old, w_old, ke, dudx, dudy, dvdx, dvdy
    USE tracers,    ONLY: tracer
    USE state,      ONLY: sigma0, rho

    TYPE(category_info_struct), INTENT(IN) :: info
    INTEGER, INTENT(IN)  :: i, j, k
    REAL(4), INTENT(IN)  :: x, y, z
    REAL(4), INTENT(INOUT) :: prop(n_prop)

    INTEGER :: m
    REAL(4) :: t
    m = 1
    DO WHILE(info%propindex_probe_tracer(1,m)/=0)
       t = tracer(i,j,k,info%propindex_probe_tracer(2,m))
       IF (t /= UNDEF) prop(info%propindex_probe_tracer(1,m)) = t
       m = m+1
    END DO

   !probe_vars = (/"U", "V", "W", "KE", "SIGMA0", "RHO", "DIVH", "ZETA"/)

    IF (info%propindex_probe(1) /= 0) prop(info%propindex_probe(1)) = 0.5*((u(i-1,j,k)+u_old(i-1,j,k))*(1.0-x)+(u(i,j,k)+u_old(i,j,k))*x)
    IF (info%propindex_probe(2) /= 0) prop(info%propindex_probe(2)) = 0.5*((v(i,j-1,k)+v_old(i,j-1,k))*(1.0-y)+(v(i,j,k)+v_old(i,j,k))*y)
    IF (info%propindex_probe(3) /= 0) prop(info%propindex_probe(3)) = 0.5*((w(i,j,k-1)+w_old(i,j,k-1))*(1.0-z)+(w(i,j,k)+w_old(i,j,k))*z)
    IF (info%propindex_probe(4) /= 0) prop(info%propindex_probe(4)) = ke(i,j,k)
    IF (info%propindex_probe(5) /= 0) prop(info%propindex_probe(5)) = sigma0(i,j,k)
    IF (info%propindex_probe(6) /= 0) prop(info%propindex_probe(6)) = rho(i,j,k)
    IF (info%propindex_probe(7) /= 0) prop(info%propindex_probe(7)) = dudx(i,j,k) + dvdy(i,j,k)
    IF (info%propindex_probe(8) /= 0) prop(info%propindex_probe(8)) = 0.25*(sum(dvdx(i-1:i,j-1:j,k))-sum(dudy(i-1:i,j-1:j,k)))

  END SUBROUTINE probe

  SUBROUTINE cumul(info, i, j, k, x, y, z, prop)
    USE tracers,    ONLY: tracer

    TYPE(category_info_struct), INTENT(IN) :: info
    INTEGER, INTENT(IN)  :: i, j, k
    REAL(4), INTENT(IN)  :: x, y, z
    REAL(4), INTENT(OUT) :: prop(n_prop)

    INTEGER :: m, n
    REAL(4) :: tmp1, tmp2

    m = 1
    DO WHILE(info%propindex_cumul_tracer(1,m)/=0)
       n = info%propindex_cumul_tracer(1,m)

       IF (prop(n)   == UNDEF) prop(n)   = 0.0
       IF (prop(n+1) == UNDEF) prop(n+1) = 0.0

       tmp1 = tracer(i,j,k,info%propindex_cumul_tracer(2,m))

       IF (tmp1==UNDEF) CYCLE

       tmp1 = tmp1*dtime - prop(n+1)
       tmp2 = prop(n) + tmp1
       prop(n+1) = (tmp2 - prop(n)) - tmp1
       prop(n)   = tmp2

       m = m+1
    END DO
  END SUBROUTINE cumul

  SUBROUTINE set_passflag(label, flag)
    REAL(4), INTENT(IN)    :: label
    REAL(4), INTENT(INOUT) :: flag

    INTEGER :: ilabel, iflag

    ilabel = INT(label, 4)
    iflag  = INT(flag,  4)

    iflag = IOR(iflag, ilabel)

    flag = REAL(iflag, 4)

  END SUBROUTINE set_passflag


  SUBROUTINE transport
    INTEGER :: ptr, next

    TYPE(particle_struct) :: tmp_particle

#ifdef PARALLEL_MPI
    INTEGER :: nsend(n_sendrecv)
    INTEGER :: nrecv(n_sendrecv)

    INTEGER :: sendreq(n_sendrecv)
    INTEGER :: recvreq(n_sendrecv)

    INTEGER :: offset(n_sendrecv)

    INTEGER :: ierr
#endif
    INTEGER :: i, j, k, n, cat

!$OMP PARALLEL PRIVATE(i, j, k, cat)

    IF (open_w .AND. icoord==0) THEN
!$OMP DO SCHEDULE(dynamic, jsize+2) COLLAPSE(3)
       DO cat=1, n_categories
       DO k=0, ksize+1
       DO j=0, jsize+1
          IF (particle_cnt(0,j,k,cat) > 0) THEN
!$OMP CRITICAL
             CALL delete_particles(0, j, k, cat)
!$OMP END CRITICAL
          END IF
       END DO
       END DO
       END DO
    END IF

    IF (open_e .AND. icoord==ipes-1) THEN
!$OMP DO SCHEDULE(dynamic, jsize+2) COLLAPSE(3)
       DO cat=1, n_categories
       DO k=0, ksize+1
       DO j=0, jsize+1
          IF (particle_cnt(isize+1,j,k,cat) > 0) THEN
!$OMP CRITICAL
             CALL delete_particles(isize+1, j, k, cat)
!$OMP END CRITICAL
          END IF
       END DO
       END DO
       END DO
    END IF

    IF (open_s .AND. jcoord==0) THEN
!$OMP DO SCHEDULE(dynamic, isize+2) COLLAPSE(3)
       DO cat=1, n_categories
       DO k=0, ksize+1
       DO i=0, isize+1
          IF (particle_cnt(i,0,k,cat) > 0) THEN
!$OMP CRITICAL
             CALL delete_particles(i, 0, k, cat)
!$OMP END CRITICAL
          END IF
       END DO
       END DO
       END DO
    END IF

    IF (open_n .AND. jcoord==jpes-1) THEN
!$OMP DO SCHEDULE(dynamic, isize+2) COLLAPSE(3)
       DO cat=1, n_categories
       DO k=0, ksize+1
       DO i=0, isize+1
          IF (particle_cnt(i,jsize+1,k,cat) > 0) THEN
!$OMP CRITICAL
             CALL delete_particles(i, jsize+1, k, cat)
!$OMP END CRITICAL
          END IF
       END DO
       END DO
       END DO
    END IF
!$OMP END PARALLEL

#ifdef PARALLEL_MPI
!$OMP PARALLEL WORKSHARE
    nsend(tag_e)   = sum(particle_cnt(isize+1,1:jsize,1:ksize,1:n_categories)) !E
    nsend(tag_w)   = sum(particle_cnt(      0,1:jsize,1:ksize,1:n_categories)) !W
    nsend(tag_n)   = sum(particle_cnt(1:isize,jsize+1,1:ksize,1:n_categories)) !N
    nsend(tag_s)   = sum(particle_cnt(1:isize,      0,1:ksize,1:n_categories)) !S

    nsend(tag_ne)  = sum(particle_cnt(isize+1,jsize+1,1:ksize,1:n_categories)) !NE
    nsend(tag_nw)  = sum(particle_cnt(      0,jsize+1,1:ksize,1:n_categories)) !NW
    nsend(tag_se)  = sum(particle_cnt(isize+1,      0,1:ksize,1:n_categories)) !SE
    nsend(tag_sw)  = sum(particle_cnt(      0,      0,1:ksize,1:n_categories)) !SW

#ifdef PARALLEL3D
    nsend(tag_u)   = sum(particle_cnt(1:isize,1:jsize,ksize+1,1:n_categories)) !U
    nsend(tag_l)   = sum(particle_cnt(1:isize,1:jsize,      0,1:n_categories)) !L

    nsend(tag_lw)  = sum(particle_cnt(      0,1:jsize,      0,1:n_categories)) !LW
    nsend(tag_le)  = sum(particle_cnt(isize+1,1:jsize,      0,1:n_categories)) !LE
    nsend(tag_ls)  = sum(particle_cnt(1:isize,      0,      0,1:n_categories)) !LS
    nsend(tag_ln)  = sum(particle_cnt(1:isize,jsize+1,      0,1:n_categories)) !LN
    nsend(tag_ue)  = sum(particle_cnt(isize+1,1:jsize,ksize+1,1:n_categories)) !UE
    nsend(tag_uw)  = sum(particle_cnt(      0,1:jsize,ksize+1,1:n_categories)) !UW
    nsend(tag_un)  = sum(particle_cnt(1:isize,jsize+1,ksize+1,1:n_categories)) !UN
    nsend(tag_us)  = sum(particle_cnt(1:isize,      0,ksize+1,1:n_categories)) !US

    nsend(tag_lne) = sum(particle_cnt(isize+1,jsize+1,      0,1:n_categories)) !LNE
    nsend(tag_lnw) = sum(particle_cnt(      0,jsize+1,      0,1:n_categories)) !LNW
    nsend(tag_lse) = sum(particle_cnt(isize+1,      0,      0,1:n_categories)) !LSE
    nsend(tag_lsw) = sum(particle_cnt(      0,      0,      0,1:n_categories)) !LSW

    nsend(tag_une) = sum(particle_cnt(isize+1,jsize+1,ksize+1,1:n_categories)) !UNE
    nsend(tag_unw) = sum(particle_cnt(      0,jsize+1,ksize+1,1:n_categories)) !UNW
    nsend(tag_use) = sum(particle_cnt(isize+1,      0,ksize+1,1:n_categories)) !USE
    nsend(tag_usw) = sum(particle_cnt(      0,      0,ksize+1,1:n_categories)) !USW
#endif
!$OMP END PARALLEL WORKSHARE

    nrecv(:) = 0

    DO n=1, n_sendrecv
       CALL mpi_isend(nsend(n), 1, MPI_INTEGER, sendrank(n), n, comm, sendreq(n), ierr)
       CALL mpi_irecv(nrecv(n), 1, MPI_INTEGER, recvrank(n), n, comm, recvreq(n), ierr)
    END DO

    CALL mpi_waitall(n_sendrecv, sendreq, MPI_STATUSES_IGNORE, ierr)
    CALL mpi_waitall(n_sendrecv, recvreq, MPI_STATUSES_IGNORE, ierr)

    IF (buffersize < max(sum(nsend), sum(nrecv))) THEN
       buffersize = 2*max(sum(nsend), sum(nrecv))
       DEALLOCATE(sendbuf)
       DEALLOCATE(recvbuf)
       ALLOCATE(sendbuf(size_of_particle,buffersize))
       ALLOCATE(recvbuf(size_of_particle,buffersize))
    END IF

    offset(1) = 1
    DO n=2, n_sendrecv
       offset(n) = offset(n-1) + nrecv(n-1)
    END DO

    recvreq(:) = MPI_REQUEST_NULL

    DO n=1, n_sendrecv
       IF (nrecv(n) > 0) CALL mpi_irecv(recvbuf(1,offset(n)), size_of_particle*nrecv(n), MPI_BYTE, recvrank(n), n, comm, recvreq(n), ierr)
    END DO

    offset(1) = 1
    DO n=2, n_sendrecv
       offset(n) = offset(n-1) + nsend(n-1)
    END DO

    sendreq(:) = MPI_REQUEST_NULL

!$OMP PARALLEL PRIVATE(n, i, j, k, cat, ptr, ierr)
!$OMP SECTIONS
!$OMP SECTION
    IF (nsend(tag_e) > 0) THEN
       n = offset(tag_e)
       DO cat=1, n_categories
       DO k=1, ksize
       DO j=1, jsize
          ptr = particle_ptr(isize+1,j,k,cat)
          DO WHILE (ptr /= NULL)
             CALL serialize(get_particle(ptr), sendbuf(:,n), convert_gridpos=.TRUE.)
             ptr = p_next(ptr)
             n = n + 1
          END DO
       END DO
       END DO
       END DO
!$OMP CRITICAL
       CALL mpi_isend(sendbuf(1,offset(tag_e)), size_of_particle*nsend(tag_e), MPI_BYTE, rank_e, tag_e, comm, sendreq(tag_e), ierr)
!$OMP END CRITICAL
    END IF

!$OMP SECTION
    IF (nsend(tag_w) > 0) THEN
       n = offset(tag_w)
       DO cat=1, n_categories
       DO k=1, ksize
       DO j=1, jsize
          ptr = particle_ptr(0,j,k,cat)
          DO WHILE (ptr /= NULL)
             CALL serialize(get_particle(ptr), sendbuf(:,n), convert_gridpos=.TRUE.)
             ptr = p_next(ptr)
             n = n + 1
          END DO
       END DO
       END DO
       END DO
!$OMP CRITICAL
       CALL mpi_isend(sendbuf(1,offset(tag_w)), size_of_particle*nsend(tag_w), MPI_BYTE, rank_w, tag_w, comm, sendreq(tag_w), ierr)
!$OMP END CRITICAL
    END IF

!$OMP SECTION
    IF (nsend(tag_n) > 0) THEN
       n = offset(tag_n)
       DO cat=1, n_categories
       DO k=1, ksize
       DO i=1, isize
          ptr = particle_ptr(i,jsize+1,k,cat)
          DO WHILE (ptr /= NULL)
             CALL serialize(get_particle(ptr), sendbuf(:,n), convert_gridpos=.TRUE.)
             ptr = p_next(ptr)
             n = n + 1
          END DO
       END DO
       END DO
       END DO
!$OMP CRITICAL
       CALL mpi_isend(sendbuf(1,offset(tag_n)), size_of_particle*nsend(tag_n), MPI_BYTE, rank_n, tag_n, comm, sendreq(tag_n), ierr)
!$OMP END CRITICAL
    END IF

!$OMP SECTION
    IF (nsend(tag_s) > 0) THEN
       n = offset(tag_s)
       DO cat=1, n_categories
       DO k=1, ksize
       DO i=1, isize
          ptr = particle_ptr(i,0,k,cat)
          DO WHILE (ptr /= NULL)
             CALL serialize(get_particle(ptr), sendbuf(:,n), convert_gridpos=.TRUE.)
             ptr = p_next(ptr)
             n = n + 1
          END DO
       END DO
       END DO
       END DO
!$OMP CRITICAL
       CALL mpi_isend(sendbuf(1,offset(tag_s)), size_of_particle*nsend(tag_s), MPI_BYTE, rank_s, tag_s, comm, sendreq(tag_s), ierr)
!$OMP END CRITICAL
    END IF

!$OMP SECTION
    IF (nsend(tag_ne) > 0) THEN
       n = offset(tag_ne)
       DO cat=1, n_categories
       DO k=1, ksize
          ptr = particle_ptr(isize+1,jsize+1,k,cat)
          DO WHILE (ptr /= NULL)
             CALL serialize(get_particle(ptr), sendbuf(:,n), convert_gridpos=.TRUE.)
             ptr = p_next(ptr)
             n = n + 1
          END DO
       END DO
       END DO
!$OMP CRITICAL
       CALL mpi_isend(sendbuf(1,offset(tag_ne)), size_of_particle*nsend(tag_ne), MPI_BYTE, rank_ne, tag_ne, comm, sendreq(tag_ne), ierr)
!$OMP END CRITICAL
    END IF

!$OMP SECTION
    IF (nsend(tag_nw) > 0) THEN
       n = offset(tag_nw)
       DO cat=1, n_categories
       DO k=1, ksize
          ptr = particle_ptr(0,jsize+1,k,cat)
          DO WHILE (ptr /= NULL)
             CALL serialize(get_particle(ptr), sendbuf(:,n), convert_gridpos=.TRUE.)
             ptr = p_next(ptr)
             n = n + 1
          END DO
       END DO
       END DO
!$OMP CRITICAL
       CALL mpi_isend(sendbuf(1,offset(tag_nw)), size_of_particle*nsend(tag_nw), MPI_BYTE, rank_nw, tag_nw, comm, sendreq(tag_nw), ierr)
!$OMP END CRITICAL
    END IF

!$OMP SECTION
    IF (nsend(tag_se) > 0) THEN
       n = offset(tag_se)
       DO cat=1, n_categories
       DO k=1, ksize
          ptr = particle_ptr(isize+1,0,k,cat)
          DO WHILE (ptr /= NULL)
             CALL serialize(get_particle(ptr), sendbuf(:,n), convert_gridpos=.TRUE.)
             ptr = p_next(ptr)
             n = n + 1
         END DO
       END DO
       END DO
!$OMP CRITICAL
       CALL mpi_isend(sendbuf(1,offset(tag_se)), size_of_particle*nsend(tag_se), MPI_BYTE, rank_se, tag_se, comm, sendreq(tag_se), ierr)
!$OMP END CRITICAL
    END IF

!$OMP SECTION
    IF (nsend(tag_sw) > 0) THEN
       n = offset(tag_sw)
       DO cat=1, n_categories
       DO k=1, ksize
          ptr = particle_ptr(0,0,k,cat)
          DO WHILE (ptr /= NULL)
             CALL serialize(get_particle(ptr), sendbuf(:,n), convert_gridpos=.TRUE.)
             ptr = p_next(ptr)
             n = n + 1
          END DO
       END DO
       END DO
!$OMP CRITICAL
       CALL mpi_isend(sendbuf(1,offset(tag_sw)), size_of_particle*nsend(tag_sw), MPI_BYTE, rank_sw, tag_sw, comm, sendreq(tag_sw), ierr)
!$OMP END CRITICAL
    END IF

#ifdef PARALLEL3D
!$OMP SECTION
    IF (nsend(tag_u) > 0) THEN
       n = offset(tag_u)
       DO cat=1, n_categories
       DO j=1, jsize
       DO i=1, isize
          ptr = particle_ptr(i,j,ksize+1,cat)
          DO WHILE (ptr /= NULL)
             CALL serialize(get_particle(ptr), sendbuf(:,n), convert_gridpos=.TRUE.)
             ptr = p_next(ptr)
             n = n + 1
          END DO
       END DO
       END DO
       END DO
!$OMP CRITICAL
       CALL mpi_isend(sendbuf(1,offset(tag_u)), size_of_particle*nsend(tag_u), MPI_BYTE, rank_u, tag_u, comm, sendreq(tag_u), ierr)
!$OMP END CRITICAL
    END IF

!$OMP SECTION
    IF (nsend(tag_l) > 0) THEN
       n = offset(tag_l)
       DO cat=1, n_categories
       DO j=1, jsize
       DO i=1, isize
          ptr = particle_ptr(i,j,0,cat)
          DO WHILE (ptr /= NULL)
             CALL serialize(get_particle(ptr), sendbuf(:,n), convert_gridpos=.TRUE.)
             ptr = p_next(ptr)
             n = n + 1
          END DO
       END DO
       END DO
       END DO
!$OMP CRITICAL
       CALL mpi_isend(sendbuf(1,offset(tag_l)), size_of_particle*nsend(tag_l), MPI_BYTE, rank_l, tag_l, comm, sendreq(tag_l), ierr)
!$OMP END CRITICAL
    END IF

!$OMP SECTION
    IF (nsend(tag_ue) > 0) THEN
       n = offset(tag_ue)
       DO cat=1, n_categories
       DO j=1, jsize
          ptr = particle_ptr(isize+1,j,ksize+1,cat)
          DO WHILE (ptr /= NULL)
             CALL serialize(get_particle(ptr), sendbuf(:,n), convert_gridpos=.TRUE.)
             ptr = p_next(ptr)
             n = n + 1
          END DO
       END DO
       END DO
!$OMP CRITICAL
       CALL mpi_isend(sendbuf(1,offset(tag_ue)), size_of_particle*nsend(tag_ue), MPI_BYTE, rank_ue, tag_ue, comm, sendreq(tag_ue), ierr)
!$OMP END CRITICAL
    END IF

!$OMP SECTION
    IF (nsend(tag_uw) > 0) THEN
       n = offset(tag_uw)
       DO cat=1, n_categories
       DO j=1, jsize
          ptr = particle_ptr(0,j,ksize+1,cat)
          DO WHILE (ptr /= NULL)
             CALL serialize(get_particle(ptr), sendbuf(:,n), convert_gridpos=.TRUE.)
             ptr = p_next(ptr)
             n = n + 1
          END DO
       END DO
       END DO
!$OMP CRITICAL
       CALL mpi_isend(sendbuf(1,offset(tag_uw)), size_of_particle*nsend(tag_uw), MPI_BYTE, rank_uw, tag_uw, comm, sendreq(tag_uw), ierr)
!$OMP END CRITICAL
    END IF

!$OMP SECTION
    IF (nsend(tag_un) > 0) THEN
       n = offset(tag_un)
       DO cat=1, n_categories
       DO i=1, isize
          ptr = particle_ptr(i,jsize+1,ksize+1,cat)
          DO WHILE (ptr /= NULL)
             CALL serialize(get_particle(ptr), sendbuf(:,n), convert_gridpos=.TRUE.)
             ptr = p_next(ptr)
             n = n + 1
          END DO
       END DO
       END DO
!$OMP CRITICAL
       CALL mpi_isend(sendbuf(1,offset(tag_un)), size_of_particle*nsend(tag_un), MPI_BYTE, rank_un, tag_un, comm, sendreq(tag_un), ierr)
!$OMP END CRITICAL
    END IF

!$OMP SECTION
    IF (nsend(tag_us) > 0) THEN
       n = offset(tag_us)
       DO cat=1, n_categories
       DO i=1, isize
          ptr = particle_ptr(i,0,ksize+1,cat)
          DO WHILE (ptr /= NULL)
             CALL serialize(get_particle(ptr), sendbuf(:,n), convert_gridpos=.TRUE.)
             ptr = p_next(ptr)
             n = n + 1
          END DO
       END DO
       END DO
!$OMP CRITICAL
       CALL mpi_isend(sendbuf(1,offset(tag_us)), size_of_particle*nsend(tag_us), MPI_BYTE, rank_us, tag_us, comm, sendreq(tag_us), ierr)
!$OMP END CRITICAL
    END IF

!$OMP SECTION
    IF (nsend(tag_le) > 0) THEN
       n = offset(tag_le)
       DO cat=1, n_categories
       DO j=1, jsize
          ptr = particle_ptr(isize+1,j,0,cat)
          DO WHILE (ptr /= NULL)
             CALL serialize(get_particle(ptr), sendbuf(:,n), convert_gridpos=.TRUE.)
             ptr = p_next(ptr)
             n = n + 1
          END DO
       END DO
       END DO
!$OMP CRITICAL
       CALL mpi_isend(sendbuf(1,offset(tag_le)), size_of_particle*nsend(tag_le), MPI_BYTE, rank_le, tag_le, comm, sendreq(tag_le), ierr)
!$OMP END CRITICAL
    END IF

!$OMP SECTION
    IF (nsend(tag_lw) > 0) THEN
       n = offset(tag_lw)
       DO cat=1, n_categories
       DO j=1, jsize
          ptr = particle_ptr(0,j,0,cat)
          DO WHILE (ptr /= NULL)
             CALL serialize(get_particle(ptr), sendbuf(:,n), convert_gridpos=.TRUE.)
             ptr = p_next(ptr)
             n = n + 1
          END DO
       END DO
       END DO
!$OMP CRITICAL
       CALL mpi_isend(sendbuf(1,offset(tag_lw)), size_of_particle*nsend(tag_lw), MPI_BYTE, rank_lw, tag_lw, comm, sendreq(tag_lw), ierr)
!$OMP END CRITICAL
    END IF

!$OMP SECTION
    IF (nsend(tag_ln) > 0) THEN
       n = offset(tag_ln)
       DO cat=1, n_categories
       DO i=1, isize
          ptr = particle_ptr(i,jsize+1,0,cat)
          DO WHILE (ptr /= NULL)
             CALL serialize(get_particle(ptr), sendbuf(:,n), convert_gridpos=.TRUE.)
             ptr = p_next(ptr)
             n = n + 1
          END DO
       END DO
       END DO
!$OMP CRITICAL
       CALL mpi_isend(sendbuf(1,offset(tag_ln)), size_of_particle*nsend(tag_ln), MPI_BYTE, rank_ln, tag_ln, comm, sendreq(tag_ln), ierr)
!$OMP END CRITICAL
    END IF

!$OMP SECTION
    IF (nsend(tag_ls) > 0) THEN
       n = offset(tag_ls)
       DO cat=1, n_categories
       DO i=1, isize
          ptr = particle_ptr(i,0,0,cat)
          DO WHILE (ptr /= NULL)
             CALL serialize(get_particle(ptr), sendbuf(:,n), convert_gridpos=.TRUE.)
             ptr = p_next(ptr)
             n = n + 1
          END DO
       END DO
       END DO
!$OMP CRITICAL
       CALL mpi_isend(sendbuf(1,offset(tag_ls)), size_of_particle*nsend(tag_ls), MPI_BYTE, rank_ls, tag_ls, comm, sendreq(tag_ls), ierr)
!$OMP END CRITICAL
    END IF

!$OMP SECTION
    IF (nsend(tag_une) > 0) THEN
       n = offset(tag_une)
       DO cat=1, n_categories
          ptr = particle_ptr(isize+1,jsize+1,ksize+1,cat)
          DO WHILE (ptr /= NULL)
             CALL serialize(get_particle(ptr), sendbuf(:,n), convert_gridpos=.TRUE.)
             ptr = p_next(ptr)
             n = n + 1
          END DO
       END DO
!$OMP CRITICAL
       CALL mpi_isend(sendbuf(1,offset(tag_une)), size_of_particle*nsend(tag_une), MPI_BYTE, rank_une, tag_une, comm, sendreq(tag_une), ierr)
!$OMP END CRITICAL
    END IF

!$OMP SECTION
    IF (nsend(tag_unw) > 0) THEN
       n = offset(tag_unw)
       DO cat=1, n_categories
          ptr = particle_ptr(0,jsize+1,ksize+1,cat)
          DO WHILE (ptr /= NULL)
             CALL serialize(get_particle(ptr), sendbuf(:,n), convert_gridpos=.TRUE.)
             ptr = p_next(ptr)
             n = n + 1
          END DO
       END DO
!$OMP CRITICAL
       CALL mpi_isend(sendbuf(1,offset(tag_unw)), size_of_particle*nsend(tag_unw), MPI_BYTE, rank_unw, tag_unw, comm, sendreq(tag_unw), ierr)
!$OMP END CRITICAL
    END IF

!$OMP SECTION
    IF (nsend(tag_use) > 0) THEN
       n = offset(tag_use)
       DO cat=1, n_categories
          ptr = particle_ptr(isize+1,0,ksize+1,cat)
          DO WHILE (ptr /= NULL)
             CALL serialize(get_particle(ptr), sendbuf(:,n), convert_gridpos=.TRUE.)
             ptr = p_next(ptr)
             n = n + 1
          END DO
       END DO
!$OMP CRITICAL
       CALL mpi_isend(sendbuf(1,offset(tag_use)), size_of_particle*nsend(tag_use), MPI_BYTE, rank_use, tag_use, comm, sendreq(tag_use), ierr)
!$OMP END CRITICAL
    END IF

!$OMP SECTION
    IF (nsend(tag_usw) > 0) THEN
       n = offset(tag_usw)
       DO cat=1, n_categories
          ptr = particle_ptr(0,0,ksize+1,cat)
          DO WHILE (ptr /= NULL)
             CALL serialize(get_particle(ptr), sendbuf(:,n), convert_gridpos=.TRUE.)
             ptr = p_next(ptr)
             n = n + 1
          END DO
       END DO
!$OMP CRITICAL
       CALL mpi_isend(sendbuf(1,offset(tag_usw)), size_of_particle*nsend(tag_usw), MPI_BYTE, rank_usw, tag_usw, comm, sendreq(tag_usw), ierr)
!$OMP END CRITICAL
    END IF

!$OMP SECTION
    IF (nsend(tag_lne) > 0) THEN
       n = offset(tag_lne)
       DO cat=1, n_categories
          ptr = particle_ptr(isize+1,jsize+1,0,cat)
          DO WHILE (ptr /= NULL)
             CALL serialize(get_particle(ptr), sendbuf(:,n), convert_gridpos=.TRUE.)
             ptr = p_next(ptr)
             n = n + 1
          END DO
       END DO
!$OMP CRITICAL
       CALL mpi_isend(sendbuf(1,offset(tag_lne)), size_of_particle*nsend(tag_lne), MPI_BYTE, rank_lne, tag_lne, comm, sendreq(tag_lne), ierr)
!$OMP END CRITICAL
    END IF

!$OMP SECTION
    IF (nsend(tag_lnw) > 0) THEN
       n = offset(tag_lnw)
       DO cat=1, n_categories
          ptr = particle_ptr(0,jsize+1,0,cat)
          DO WHILE (ptr /= NULL)
             CALL serialize(get_particle(ptr), sendbuf(:,n), convert_gridpos=.TRUE.)
             ptr = p_next(ptr)
             n = n + 1
          END DO
       END DO
!$OMP CRITICAL
       CALL mpi_isend(sendbuf(1,offset(tag_lnw)), size_of_particle*nsend(tag_lnw), MPI_BYTE, rank_lnw, tag_lnw, comm, sendreq(tag_lnw), ierr)
!$OMP END CRITICAL
    END IF

!$OMP SECTION
    IF (nsend(tag_lse) > 0) THEN
       n = offset(tag_lse)
       DO cat=1, n_categories
          ptr = particle_ptr(isize+1,0,0,cat)
          DO WHILE (ptr /= NULL)
             CALL serialize(get_particle(ptr), sendbuf(:,n), convert_gridpos=.TRUE.)
             ptr = p_next(ptr)
             n = n + 1
          END DO
       END DO
!$OMP CRITICAL
       CALL mpi_isend(sendbuf(1,offset(tag_lse)), size_of_particle*nsend(tag_lse), MPI_BYTE, rank_lse, tag_lse, comm, sendreq(tag_lse), ierr)
!$OMP END CRITICAL
    END IF

!$OMP SECTION
    IF (nsend(tag_lsw) > 0) THEN
       n = offset(tag_lsw)
       DO cat=1, n_categories
          ptr = particle_ptr(0,0,0,cat)
          DO WHILE (ptr /= NULL)
             CALL serialize(get_particle(ptr), sendbuf(:,n), convert_gridpos=.TRUE.)
             ptr = p_next(ptr)
             n = n + 1
          END DO
       END DO
!$OMP CRITICAL
       CALL mpi_isend(sendbuf(1,offset(tag_lsw)), size_of_particle*nsend(tag_lsw), MPI_BYTE, rank_lsw, tag_lsw, comm, sendreq(tag_lsw), ierr)
!$OMP END CRITICAL
    END IF
!$OMP END SECTIONS
!$OMP END PARALLEL
#endif


!$OMP PARALLEL PRIVATE(i,j,k,cat)
    DO cat=1, n_categories
!$OMP DO SCHEDULE(dynamic, jsize+2) COLLAPSE(2)
       DO k=1, ksize
       DO j=0, jsize+1
          IF (particle_cnt(0,j,k,cat) > 0) THEN
!$OMP CRITICAL
             CALL delete_particles(0, j, k, cat, record=.FALSE.)
!$OMP END CRITICAL
          END IF
          IF (particle_cnt(isize+1,j,k,cat) > 0) THEN
!$OMP CRITICAL
             CALL delete_particles(isize+1, j, k, cat, record=.FALSE.)
!$OMP END CRITICAL
          END IF
       END DO
       END DO

!$OMP DO SCHEDULE(dynamic, isize+2) COLLAPSE(2)
       DO k=1, ksize
       DO i=0, isize+1
          IF (particle_cnt(i,0,k,cat) > 0) THEN
!$OMP CRITICAL
             CALL delete_particles(i, 0, k, cat, record=.FALSE.)
!$OMP END CRITICAL
          END IF
          IF (particle_cnt(i,jsize+1,k,cat) > 0) THEN
!$OMP CRITICAL
             CALL delete_particles(i, jsize+1, k, cat, record=.FALSE.)
!$OMP END CRITICAL
          END IF
       END DO
       END DO

#ifdef PARALLEL3D
!$OMP DO SCHEDULE(dynamic, isize+2) COLLAPSE(2)
       DO j=0, jsize+1
       DO i=0, isize+1
          IF (particle_cnt(i,j,0,cat) > 0) THEN
!$OMP CRITICAL
             CALL delete_particles(i, j, 0, cat, record=.FALSE.)
!$OMP END CRITICAL
          END IF
          IF (particle_cnt(i,j,ksize+1,cat) > 0) THEN
!$OMP CRITICAL
             CALL delete_particles(i, j, ksize+1, cat, record=.FALSE.)
!$OMP END CRITICAL
          END IF
       END DO
       END DO
#endif
    END DO
!$OMP END PARALLEL

    CALL mpi_waitall(n_sendrecv, recvreq, MPI_STATUSES_IGNORE, ierr)

    DO n=1, sum(nrecv)
       CALL deserialize(recvbuf(:,n), tmp_particle, convert_gridpos=.TRUE.)
       CALL create_particle(tmp_particle, record=.FALSE.)
    END DO

    CALL mpi_waitall(n_sendrecv, sendreq, MPI_STATUSES_IGNORE, ierr)

#else
    IF (cycle_x) THEN
       DO cat=1, n_categories
       DO k=0, ksize+1
       DO j=0, jsize+1
          IF (particle_cnt(isize+1,j,k,cat) /= 0) THEN
             ptr = particle_ptr(isize+1,j,k,cat)
             DO
                p_ijk(1,ptr) = 1
                next = p_next(ptr)
                IF (next==NULL) THEN
                   IF (particle_ptr(1,j,k,cat)/=NULL) p_prev(particle_ptr(1,j,k,cat)) = ptr
                   p_next(ptr)  = particle_ptr(1,j,k,cat)
                   particle_ptr(1,j,k,cat) = particle_ptr(isize+1,j,k,cat)
                   EXIT
                END IF
                ptr = next
             END DO
             particle_ptr(isize+1,j,k,cat) = NULL
             particle_cnt(1,j,k,cat) = particle_cnt(1,j,k,cat) + particle_cnt(isize+1,j,k,cat)
             particle_cnt(isize+1,j,k,cat) = 0
          END IF

          IF (particle_cnt(0,j,k,cat) /= 0) THEN
             ptr = particle_ptr(0,j,k,cat)
             DO
                p_ijk(1,ptr) = isize
                next = p_next(ptr)
                IF (next==NULL) THEN
                   IF (particle_ptr(isize,j,k,cat)/=NULL) p_prev(particle_ptr(isize,j,k,cat)) = ptr
                   p_next(ptr)  = particle_ptr(isize,j,k,cat)
                   particle_ptr(isize,j,k,cat) = particle_ptr(0,j,k,cat)
                   EXIT
                END IF
                ptr = next
             END DO
             particle_ptr(0,j,k,cat) = NULL
             particle_cnt(isize,j,k,cat) = particle_cnt(isize,j,k,cat) + particle_cnt(0,j,k,cat)
             particle_cnt(0,j,k,cat) = 0
          END IF
       END DO
       END DO
       END DO
    END IF

    IF (cycle_y) THEN
       DO cat=1, n_categories
       DO k=0, ksize+1
       DO i=0, isize+1
          IF (particle_cnt(i,jsize+1,k,cat) /= 0) THEN
             ptr = particle_ptr(i,jsize+1,k,cat)
             DO
                p_ijk(2,ptr) = 1
                next = p_next(ptr)
                IF (next==NULL) THEN
                   IF (particle_ptr(i,1,k,cat)/=NULL) p_prev(particle_ptr(i,1,k,cat)) = ptr
                   p_next(ptr)  = particle_ptr(i,1,k,cat)
                   particle_ptr(i,1,k,cat) = particle_ptr(i,jsize+1,k,cat)
                   EXIT
                END IF
                ptr = next
             END DO
             particle_ptr(i,jsize+1,k,cat) = NULL
             particle_cnt(i,1,k,cat) = particle_cnt(i,1,k,cat) + particle_cnt(i,jsize+1,k,cat)
             particle_cnt(i,jsize+1,k,cat) = 0
          END IF

          IF (particle_cnt(i,0,k,cat) /= 0) THEN
             ptr = particle_ptr(i,0,k,cat)
             DO
                p_ijk(2,ptr) = jsize
                next = p_next(ptr)
                IF (next==NULL) THEN
                   IF (particle_ptr(i,jsize,k,cat)/=NULL) p_prev(particle_ptr(i,jsize,k,cat)) = ptr
                   p_next(ptr)      = particle_ptr(i,jsize,k,cat)
                   particle_ptr(i,jsize,k,cat) = particle_ptr(i,0,k,cat)
                   EXIT
                END IF
                ptr = next
             END DO
             particle_ptr(i,0,k,cat) = NULL
             particle_cnt(i,jsize,k,cat) = particle_cnt(i,jsize,k,cat) + particle_cnt(i,0,k,cat)
             particle_cnt(i,0,k,cat) = 0
          END IF
       END DO
       END DO
       END DO
    END IF

    IF (cycle_z) THEN
       DO cat=1, n_categories
       DO j=0, jsize+1
       DO i=0, isize+1
          IF (particle_cnt(i,j,ksize+1,cat) /= 0) THEN
             ptr = particle_ptr(i,j,ksize+1,cat)
             DO
                p_ijk(3,ptr) = 1
                next = p_next(ptr)
                IF (next==NULL) THEN
                   IF (particle_ptr(i,j,1,cat)/=NULL) p_prev(particle_ptr(i,j,1,cat)) = ptr
                   p_next(ptr)  = particle_ptr(i,j,1,cat)
                   particle_ptr(i,j,1,cat) = particle_ptr(i,j,ksize+1,cat)
                   EXIT
                END IF
                ptr = next
             END DO
             particle_ptr(i,j,ksize+1,cat) = NULL
             particle_cnt(i,j,1,cat) = particle_cnt(i,j,1,cat) + particle_cnt(i,j,ksize+1,cat)
             particle_cnt(i,j,ksize+1,cat) = 0
          END IF

          IF (particle_cnt(i,j,0,cat) /= 0) THEN
             ptr = particle_ptr(i,j,0,cat)
             DO
                p_ijk(3,ptr) = ksize
                next = p_next(ptr)
                IF (next==NULL) THEN
                   IF (particle_ptr(i,j,ksize,cat)/=NULL) p_prev(particle_ptr(i,j,ksize,cat)) = ptr
                   p_next(ptr)      = particle_ptr(i,j,ksize,cat)
                   particle_ptr(i,j,ksize,cat) = particle_ptr(i,j,0,cat)
                   EXIT
                END IF
                ptr = next
             END DO
             particle_ptr(i,j,0,cat) = NULL
             particle_cnt(i,j,ksize,cat) = particle_cnt(i,j,ksize,cat) + particle_cnt(i,j,0,cat)
             particle_cnt(i,j,0,cat) = 0
          END IF
       END DO
       END DO
       END DO
    END IF

    IF (tripolar) THEN
       DO cat=1, n_categories
       DO k=1, dimz
       DO i=1, dimx
          IF (particle_cnt(i,dimy+1,k,cat) /= 0) THEN

             ptr = particle_ptr(i,dimy+1,k,cat)
             DO
                p_ijk(2,ptr) = dimy
                p_ijk(1,ptr) = dimx - i + 1
                p_xyz(1,ptr) = 1.0 - p_xyz(1,ptr)
                p_xyz(2,ptr) = 1.0 - p_xyz(2,ptr)

                next = p_next(ptr)
                IF (next==NULL) THEN
                   IF (particle_ptr(dimx-i+1,dimy,k,cat)/=NULL) p_prev(particle_ptr(dimx-i+1,dimy,k,cat)) = ptr
                   p_next(ptr)  = particle_ptr(dimx-i+1,dimy,k,cat)
                   particle_ptr(dimx-i+1,dimy,k,cat) = particle_ptr(i,dimy+1,k,cat)
                   EXIT
                END IF
                ptr = next
             END DO
             particle_ptr(i,dimy+1,k,cat) = NULL
             particle_cnt(dimx-i+1,dimy,k,cat) = particle_cnt(dimx-i+1,dimy,k,cat) + particle_cnt(i,dimy+1,k,cat)
             particle_cnt(i,dimy+1,k,cat) = 0
          END IF
       END DO
       END DO
       END DO
    END IF

#endif
  END SUBROUTINE transport

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE step_particles
    INTEGER :: n, m
    INTEGER :: hist_result(0:32)
PROFILE_BEGIN('pinput')
    DO n=1, n_input
       CALL input(input_info(n))
    END DO
PROFILE_END('pinput')

    CALL spawn
    CALL main
    CALL transport

    CALL checkout_particles

    DO n=1, n_categories
       IF (category_info(n)%hist_prop==0) CYCLE
       IF (t_current - INT(t_current/category_info(n)%hist_interval)*category_info(n)%hist_interval < dtime) THEN
          CALL particle_histogram(n, category_info(n)%hist_prop, category_info(n)%hist_bins, hist_result, all=.FALSE.)

          IF (rank==0) THEN
             WRITE(REPORT_UNIT, '(A, ": particle histogram for ",A,"_",A)', ADVANCE='NO') &
                  trim(current_datetime), trim(category_info(n)%name), trim(category_info(n)%property_name(category_info(n)%hist_prop))
             WRITE(REPORT_UNIT, '(" ",I0)', ADVANCE='NO') hist_result(0)
             DO m=1, 32
                IF (category_info(n)%hist_bins(m)==UNDEF) EXIT
                WRITE(REPORT_UNIT, '(", ",I0)', ADVANCE='NO') hist_result(m)
             END DO
             WRITE(REPORT_UNIT, *)
          END IF
       END IF
    END DO

  END SUBROUTINE step_particles

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE spawn
    USE velocity, ONLY: u, v, w

    TYPE(particle_struct) :: tmp_particle

    INTEGER :: i, j, k, n, nn, ijk, cat
    LOGICAL :: stat(9)

    REAL(8) :: source(0:isize,0:jsize,0:ksize)
    REAL(8) :: source_yz(0:jsize,0:ksize)
    REAL(8) :: source_xz(0:isize,0:ksize)

    INTEGER :: propindex_spawn(3)

    DO cat = 1, n_categories
       stat(1) = require_checkin('SPAWN_'      //category_info(cat)%name)
       stat(2) = require_checkin('SPAWN_ABSV_' //category_info(cat)%name)
       stat(3) = require_checkin('SPAWN_FX_'   //category_info(cat)%name)
       stat(4) = require_checkin('SPAWN_FY_'   //category_info(cat)%name)
       stat(5) = require_checkin('SPAWN_FZ_'   //category_info(cat)%name)
       stat(6) = require_checkin('SPAWN_WEST_' //category_info(cat)%name)
       stat(7) = require_checkin('SPAWN_EAST_' //category_info(cat)%name)
       stat(8) = require_checkin('SPAWN_SOUTH_'//category_info(cat)%name)
       stat(9) = require_checkin('SPAWN_NORTH_'//category_info(cat)%name)

       IF (.NOT. any(stat)) CYCLE

       IF (.NOT. ASSOCIATED(category_info(cat)%aspawn)) THEN
          ALLOCATE(category_info(cat)%aspawn(isize,jsize,ksize))
          category_info(cat)%aspawn(:,:,:) = 0.0
          CALL initial_data('ASPAWN_'//category_info(cat)%name, category_info(cat)%aspawn, default=perfect_restart)
       END IF

       IF (stat(1)) THEN
          CALL checkin('SPAWN_'//category_info(cat)%name, source)
!$OMP PARALLEL DO
          DO k=1, ksize
          DO j=1, jsize
          DO i=1, isize
             IF (.NOT. lmask3d(i,j,k)) CYCLE
             IF (source(i,j,k)==0.0)   CYCLE
             category_info(cat)%aspawn(i,j,k) = category_info(cat)%aspawn(i,j,k) + source(i,j,k)*dtime
          END DO
          END DO
          END DO
       END IF

       IF (stat(2)) THEN
          CALL checkin('SPAWN_ABSV_'//category_info(cat)%name, source)
!$OMP PARALLEL DO
          DO k=1, ksize
          DO j=1, jsize
          DO i=1, isize
             IF (.NOT. lmask3d(i,j,k)) CYCLE
             IF (source(i,j,k)==0.0)   CYCLE
             category_info(cat)%aspawn(i,j,k) = category_info(cat)%aspawn(i,j,k) &
                  + source(i,j,k) * 0.5*sqrt((u(i,j,k)+u(i-1,j,k))**2+(v(i,j,k)+v(i,j-1,k))**2+(w(i,j,k)+w(i,j,k-1))**2) * dtime
          END DO
          END DO
          END DO
       END IF

       IF (stat(3)) THEN
          CALL checkin('SPAWN_FX_'//category_info(cat)%name, source)
!$OMP PARALLEL DO
          DO k=1, ksize
          DO j=1, jsize
          DO i=1, isize
             IF (.NOT. lmask3d(i,j,k)) CYCLE
             IF (source(i-1,j,k) <= 0.0 .AND. source(i,j,k) >= 0.0) CYCLE
             category_info(cat)%aspawn(i,j,k) = category_info(cat)%aspawn(i,j,k) &
                  + (max(source(i-1,j,k), 0.0D0)*max(u(i-1,j,k), 0.0D0)*dsx(i-1,j,k) &
                   + min(source(i,  j,k), 0.0D0)*min(u(i,  j,k), 0.0D0)*dsx(i,  j,k)) * dtime
          END DO
          END DO
          END DO
       END IF

       IF (stat(4)) THEN
          CALL checkin('SPAWN_FY_'//category_info(cat)%name, source)
!$OMP PARALLEL DO
          DO k=1, ksize
          DO j=1, jsize
          DO i=1, isize
             IF (.NOT. lmask3d(i,j,k)) CYCLE
             IF (source(i,j-1,k) <= 0.0 .AND. source(i,j,k) >= 0.0) CYCLE
             category_info(cat)%aspawn(i,j,k) = category_info(cat)%aspawn(i,j,k) &
                  + (max(source(i,j-1,k), 0.0D0)*max(v(i,j-1,k), 0.0D0)*dsy(i,j-1,k) &
                   + min(source(i,j,  k), 0.0D0)*min(v(i,j,  k), 0.0D0)*dsy(i,j,  k)) * dtime
          END DO
          END DO
          END DO
       END IF

       IF (stat(5)) THEN
          CALL checkin('SPAWN_FZ_'//category_info(cat)%name, source)
!$OMP PARALLEL DO
          DO k=1, ksize
          DO j=1, jsize
          DO i=1, isize
             IF (.NOT. lmask3d(i,j,k)) CYCLE
             IF (source(i,j,k-1) <= 0.0 .AND. source(i,j,k) >= 0.0) CYCLE
             category_info(cat)%aspawn(i,j,k) = category_info(cat)%aspawn(i,j,k) &
                  + (max(source(i,j,k-1), 0.0D0)*max(w(i,j,k-1), 0.0D0)*dsz(i,j) &
                   + min(source(i,j,k),   0.0D0)*min(w(i,j,k),   0.0D0)*dsz(i,j)) * dtime
          END DO
          END DO
          END DO
       END IF

       IF (stat(6)) THEN
          CALL checkin('SPAWN_WEST_'//category_info(cat)%name, source_yz, section='YZ')

          IF (open_w .AND. icoord==0) THEN
             i = 11 !10-grid margin from the open boundary
!$OMP PARALLEL DO
             DO k=1, ksize
             DO j=1, jsize
                IF (.NOT. lmask3d(i,j,k)) CYCLE
                IF (source_yz(j,k)<=0.0)  CYCLE
                category_info(cat)%aspawn(i,j,k) = category_info(cat)%aspawn(i,j,k) &
                     + source_yz(j,k)*max(u(i-1,j,k), 0.0D0)*dsx(i-1,j,k) * dtime
             END DO
             END DO
          END IF
       END IF

       IF (stat(7)) THEN
          CALL checkin('SPAWN_EAST_'//category_info(cat)%name, source_yz, section='YZ')

          IF (open_e .AND. icoord==ipes-1) THEN
             i = isize-10
!$OMP PARALLEL DO
             DO k=1, ksize
             DO j=1, jsize
                IF (.NOT. lmask3d(i,j,k)) CYCLE
                IF (source_yz(j,k)<=0.0)  CYCLE
                category_info(cat)%aspawn(i,j,k) = category_info(cat)%aspawn(i,j,k) &
                     - source_yz(j,k)*min(u(i,j,k), 0.0D0)*dsx(i,j,k) * dtime
             END DO
             END DO
          END IF
       END IF

       IF (stat(8)) THEN
          CALL checkin('SPAWN_SOUTH_'//category_info(cat)%name, source_xz, section='XZ')

          IF (open_s .AND. jcoord==0) THEN
             j = 11
!$OMP PARALLEL DO
             DO k=1, ksize
             DO i=1, isize
                IF (.NOT. lmask3d(i,j,k)) CYCLE
                IF (source_xz(i,k)<=0.0)  CYCLE
                category_info(cat)%aspawn(i,j,k) = category_info(cat)%aspawn(i,j,k) &
                     + source_xz(i,k)*max(v(i,j-1,k), 0.0D0)*dsy(i,j-1,k) * dtime
             END DO
             END DO
          END IF
       END IF

       IF (stat(9)) THEN
          CALL checkin('SPAWN_NORTH_'//category_info(cat)%name, source_xz, section='XZ')

          IF (open_n .AND. jcoord==jpes-1) THEN
             j = jsize-10
!$OMP PARALLEL DO
             DO k=1, ksize
             DO i=1, isize
                IF (.NOT. lmask3d(i,j,k)) CYCLE
                IF (source_xz(i,k)<=0.0)  CYCLE
                category_info(cat)%aspawn(i,j,k) = category_info(cat)%aspawn(i,j,k) &
                     - source_xz(i,k)*min(v(i,j,k), 0.0D0)*dsy(i,j,k) * dtime
             END DO
             END DO
          END IF
       END IF

!$OMP PARALLEL PRIVATE(i, j, k, n, nn, tmp_particle)
       tmp_particle%id   = 0_8
       tmp_particle%category = cat
       tmp_particle%xpos = UNDEF
       tmp_particle%ypos = UNDEF
       tmp_particle%zpos = UNDEF
       tmp_particle%seed = 0
       tmp_particle%property(:) = UNDEF
       tmp_particle%track_index = -1

!$OMP DO SCHEDULE(dynamic, ompchunk)
       DO ijk=0, isize*jsize*ksize-1
          k = ijk/(isize*jsize) + 1
          j = (ijk - (k-1)*isize*jsize)/isize + 1
          i = ijk - (k-1)*isize*jsize - (j-1)*isize + 1

          IF (.NOT. lmask3d(i,j,k) .OR. category_info(cat)%aspawn(i,j,k) < 1.0) CYCLE

          tmp_particle%ipos = i
          tmp_particle%jpos = j
          tmp_particle%kpos = k

          nn = int(category_info(cat)%aspawn(i,j,k))
!$OMP CRITICAL
          DO n=1, nn
             CALL create_particle(tmp_particle)
          END DO
!$OMP END CRITICAL
          category_info(cat)%aspawn(i,j,k) = category_info(cat)%aspawn(i,j,k) - nn
       END DO
!$OMP END PARALLEL

       CALL checkout('ASPAWN_'//category_info(cat)%name, category_info(cat)%aspawn)
    END DO

    CALL deploy_npzd

  END SUBROUTINE spawn

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE checkout_particles
    REAL(8)       :: tmp(1:isize, 1:jsize, 1:ksize)
    CHARACTER(64) :: varname

    LOGICAL :: watchmask(isize,jsize,ksize)
    LOGICAL :: stat

    INTEGER :: i, j, k, ijk, n, cat, p

    stat = watch_record%initialized .AND. t_current >= watch_record%start .AND. t_current < watch_record%end
    IF (stat) CALL checkin('WATCHMASK', watchmask, stat=stat)

    DO cat=1, n_categories
       IF (stat .AND. category_info(cat)%record_watch) THEN
!$OMP PARALLEL DO SCHEDULE(dynamic, ompchunk) PRIVATE(i, j, k, p)
          DO ijk=0, isize*jsize*ksize-1
             k = ijk/(isize*jsize) + 1
             j = (ijk - (k-1)*isize*jsize)/isize + 1
             i = ijk - (k-1)*isize*jsize - (j-1)*isize + 1

             IF (.NOT. watchmask(i,j,k) .OR. particle_cnt(i,j,k,cat)==0) CYCLE

!$OMP CRITICAL
             p = particle_handle(i,j,k,cat)
             DO WHILE (p /= NULL)
                CALL append_record(get_particle(p, seed=n_timestep), watch_record)
                p = p_next(p)
             END DO
!$OMP END CRITICAL
          END DO
       END IF

       varname = trim(category_info(cat)%name)

       IF (require_checkout(trim(varname))) THEN
!$OMP PARALLEL DO
          DO k=1, ksize
          DO j=1, jsize
          DO i=1, isize
             tmp(i,j,k) = particle_cnt(i,j,k,cat) / dvol(i,j,k)
          END DO
          END DO
          END DO
          CALL checkout(trim(varname), tmp)
       END IF

       DO n=1, n_prop
          varname = trim(category_info(cat)%name) // '_' // trim(category_info(cat)%property_name(n))
          IF (require_checkout(trim(varname))) THEN
!$OMP PARALLEL DO SCHEDULE(dynamic, ompchunk) PRIVATE(i, j, k, p)
             DO ijk=0, isize*jsize*ksize-1
                k = ijk/(isize*jsize) + 1
                j = (ijk - (k-1)*isize*jsize)/isize + 1
                i = ijk - (k-1)*isize*jsize - (j-1)*isize + 1

                tmp(i,j,k) = 0.0
                p = particle_handle(i,j,k,cat)
                DO WHILE (p /= NULL)
                   tmp(i,j,k) = tmp(i,j,k) + p_property(n,p)
                   p = p_next(p)
                END DO
                tmp(i,j,k) = tmp(i,j,k) / dvol(i,j,k)
             END DO
             CALL checkout(trim(varname), tmp)
          END IF
       END DO
    END DO

    IF (require_checkout('PARTICLE_NCREATE')) THEN
       CALL assert(maxval(ncreate) < 2**id_nshift, "PARTICLE_NCREATE exceeds the limit")
!$OMP PARALLEL WORKSHARE
        tmp(:,:,:) = INT(ncreate(:,:,:),8)
!$OMP END PARALLEL WORKSHARE
       CALL checkout('PARTICLE_NCREATE', tmp)
    END IF


  END SUBROUTINE checkout_particles

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE particle_density(result, category_name, property_name)
    REAL(8), INTENT(OUT) :: result(isize, jsize, ksize)
    CHARACTER(*), INTENT(IN)           :: category_name
    CHARACTER(*), INTENT(IN), OPTIONAL :: property_name

    INTEGER :: cat, prop
    INTEGER :: i, j, k, ijk, p

    cat = category_index(category_name)
    CALL assert(cat /= 0, "category_name '" // category_name // "' is not defined")

    IF (PRESENT(property_name)) THEN
       prop = property_index(cat, property_name)
       CALL assert(prop /= 0,  "propety_name '"// property_name // "' is not defined for category '" // category_name // "'")

!$OMP PARALLEL DO SCHEDULE(dynamic, ompchunk) PRIVATE(i, j, k, p)
       DO ijk=0, isize*jsize*ksize-1
          k = ijk/(isize*jsize) + 1
          j = (ijk - (k-1)*isize*jsize)/isize + 1
          i = ijk - (k-1)*isize*jsize - (j-1)*isize + 1

          result(i,j,k) = 0.0
          p = particle_handle(i,j,k,cat)
          DO WHILE (p /= NULL)
             result(i,j,k) = result(i,j,k) + p_property(prop,p)
             p = p_next(p)
          END DO
          result(i,j,k) = result(i,j,k) / dvol(i,j,k)
       END DO
    ELSE
!$OMP PARALLEL DO
       DO k=1, ksize
       DO j=1, jsize
       DO i=1, isize
          result(i,j,k) = particle_cnt(i,j,k,cat) / dvol(i,j,k)
       END DO
       END DO
       END DO
    END IF

  END SUBROUTINE particle_density


!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE flush_particles
    IF (grave_record%initialized) THEN
       IF (check_time(t_current, grave_record%start, grave_record%end, grave_record%interval, grave_record%lastsync)) &
          CALL flush_record(grave_record, all=.TRUE., timestamp=.TRUE.)
    END IF

    IF (birth_record%initialized) THEN
       IF (check_time(t_current, birth_record%start, birth_record%end, birth_record%interval, birth_record%lastsync)) &
          CALL flush_record(birth_record, all=.TRUE., timestamp=.TRUE.)
    END IF

    IF (watch_record%initialized) THEN
       IF (check_time(t_current, watch_record%start, watch_record%end, watch_record%interval, watch_record%lastsync)) &
          CALL flush_record(watch_record, all=.TRUE., timestamp=.TRUE.)
    END IF

    IF (.NOT. check_time(t_current, output_info%start, output_info%end, output_info%interval, output_info%lastwrite)) RETURN

    CALL check_particles

    output_info%lastwrite = t_current - mod(t_current - output_info%start, output_info%interval)

    CALL write_particles(trim(output_info%filepath) // '.' // trim(format_datetime(output_info%lastwrite, omit_time=output_info%omit_time)))

  END SUBROUTINE flush_particles

!-----------------------------------------------------------------------------------------------------------------------

  LOGICAL FUNCTION in_heap(p)
    INTEGER, INTENT(IN) :: p
    INTEGER :: ptr

    in_heap = .FALSE.

    ptr = heap_ptr

    DO WHILE(ptr /= NULL)
       IF (ptr==p) THEN
          in_heap = .TRUE.
          RETURN
       END IF

       ptr = p_next(ptr)
    END DO

  END FUNCTION in_heap

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE check_heap
    INTEGER :: ptr
    INTEGER(8) :: n

    ptr = heap_ptr
    n = 0

    DO WHILE (ptr /= NULL)
       CALL assert(p_id(ptr) == 0, "particle-ID in the heap should be 0")

       IF (p_next(ptr) /= NULL) THEN
          CALL assert(p_prev(p_next(ptr))==ptr, "p_prev and p_next are inconsistent")
       END IF

       ptr = p_next(ptr)
       n = n + 1
    END DO
    CALL assert(heap_cnt==n, "heap_cnt incorrect")

  END SUBROUTINE check_heap

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE check_particles
    INTEGER :: i, j, k, n, cat
    INTEGER(8) :: num(n_categories)
    INTEGER :: ptr

#ifdef PARALLEL_MPI
    INTEGER :: ierr
#endif

!    CALL check_heap
    num(:) = 0

    DO cat=1, n_categories
       DO k=1, ksize
       DO j=1, jsize
       DO i=1, isize
          ptr = particle_ptr(i,j,k,cat)

        ! IF (ptr /= NULL) CALL assert(.NOT. in_heap(ptr), "particle_ptr(i,j,k) is in heap !!!!")
          IF (ptr /= NULL) CALL assert(p_prev(ptr)==NULL, "p_prev of head-of-list particle is not NULL")

          n = 0
          DO WHILE (ptr /= NULL)
             !  CALL assert(.NOT. in_heap(ptr), "an active particle ptr in heap")

             CALL assert(p_id(ptr) > 0, "an invalid particle contaminates")

             CALL assert(p_ijk(1,ptr)==i, "particle%ipos incorrect")
             CALL assert(p_ijk(2,ptr)==j, "particle%jpos incorrect")
             CALL assert(p_ijk(3,ptr)==k, "particle%kpos incorrect")

             CALL assert(p_xyz(1,ptr) >= 0.0 .AND. p_xyz(1,ptr) <= 1.0, "particle%xpos is not in the valid range")
             CALL assert(p_xyz(2,ptr) >= 0.0 .AND. p_xyz(2,ptr) <= 1.0, "particle%ypos is not in the valid range")
             CALL assert(p_xyz(3,ptr) >= 0.0 .AND. p_xyz(3,ptr) <= 1.0, "particle%zpos is not in the valid range")

             IF (p_next(ptr) /= NULL) THEN
                CALL assert(p_prev(p_next(ptr))==ptr, "p_prev is invalid")
             END IF

             ptr = p_next(ptr)
             n = n + 1
          END DO

          CALL assert(particle_cnt(i,j,k,cat) == n, "particle_cnt inconsistent")
          CALL assert(logical(lmask3d(i,j,k)) .OR. n==0,   "particles in masked cell")

          num(cat) = num(cat) + n
       END DO
       END DO
       END DO
    END DO

    CALL assert(sum(num) + heap_cnt >= max_particles, "some particles are lost from link-list")
    CALL assert(sum(num) + heap_cnt <= max_particles, "some particles are double-assigned")

    CALL gsum(num)

    IF (rank==0) THEN
       WRITE(REPORT_UNIT, '(A,I0,A)', ADVANCE='NO') "total number of particles: ", sum(num), " ("
       DO cat=1, n_categories
          WRITE(REPORT_UNIT, '(A,I0)', ADVANCE='NO') trim(category_info(cat)%name)//":", num(cat)
          IF (cat==n_categories) THEN
             WRITE(REPORT_UNIT, '(A)') ")"
          ELSE
             WRITE(REPORT_UNIT, '(A)', ADVANCE='NO') ", "
          END IF
       END DO
    END IF

    CALL gsum(nlimit)
    IF (rank==0 .AND. sum(nlimit)>0) WRITE(REPORT_UNIT, '("CFL limiter trigered (X:", I0, ", Y:", I0, ", Z:",I0, ")" )') nlimit
    nlimit(:) = 0

  END SUBROUTINE check_particles

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE init_ncreate
    INTEGER(8) :: imax
    INTEGER    :: i, j, k, ijk, cat, ptr

    REAL(8)    :: tmp(isize,jsize,ksize)

    ALLOCATE(ncreate(isize,jsize,ksize))
!$OMP PARALLEL WORKSHARE
    tmp(:,:,:) = 0.0
!$OMP END PARALLEL WORKSHARE
    CALL initial_data("PARTICLE_NCREATE", tmp, default=perfect_restart)
!$OMP PARALLEL WORKSHARE
    ncreate(:,:,:) = INT(tmp)
!$OMP END PARALLEL WORKSHARE

  END SUBROUTINE init_ncreate

!-----------------------------------------------------------------------------------------------------------------------

   SUBROUTINE register_track(id, filepath, csv, append, start, end)
    INTEGER(8),   INTENT(IN) :: id
    CHARACTER(*), INTENT(IN) :: filepath
    LOGICAL,      INTENT(IN) :: csv
    LOGICAL,      INTENT(IN) :: append
    CHARACTER(*), INTENT(IN) :: start
    CHARACTER(*), INTENT(IN) :: end

    INTEGER :: ptr
    INTEGER :: n

#ifdef PARALLEL_MPI
    INTEGER :: ierr
#endif

    DO n = 1, n_tracks
       CALL assert(id /= track_id(n), "particle ID=" // trim(format(id)) // " is already registered")
       CALL assert(trim(filepath)/=trim(track_record(n)%filepath), "the same filename is assigned for different particles")
    END DO

    n_tracks = n_tracks + 1

    CALL assert(n_tracks <= max_tracks, "number of particle-track exceeds MAX_TRACKS")

    IF (n_tracks == 1) THEN
       max_track_id = id
       min_track_id = id
    ELSE
       max_track_id = max(max_track_id, id)
       min_track_id = min(min_track_id, id)
    END IF

    track_id(n_tracks) = id

    CALL init_record(trim(filepath), csv, append, start, end, "", track_record(n_tracks))

  END SUBROUTINE register_track

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE init_record(filepath, csv, append, start, end, interval, record)
    CHARACTER(*), INTENT(IN) :: filepath
    LOGICAL,      INTENT(IN) :: csv
    LOGICAL,      INTENT(IN) :: append
    CHARACTER(*), INTENT(IN) :: start
    CHARACTER(*), INTENT(IN) :: end
    CHARACTER(*), INTENT(IN) :: interval
    TYPE(record_struct), INTENT(INOUT) :: record

    CHARACTER(160) :: line
    INTEGER :: iostat, ierr
    TYPE(particle_struct) :: ts
    INTEGER(1) :: buf(size_of_particle)

    record%initialized = .TRUE.

    record%filepath = trim(filepath)
    record%count    = 0
    record%n_write  = 0
    record%bufsize  = chunksize
    ALLOCATE(record%buffer(size_of_particle, record%bufsize))

    IF (rank==0 .AND. .NOT. append) THEN
       OPEN(TMP_UNIT, FILE=record%filepath, FORM='UNFORMATTED', STATUS='REPLACE', ACTION='WRITE', IOSTAT=iostat)
       CALL assert(iostat==0, "failed to create file '"//trim(record%filepath)//"'")
       CLOSE(TMP_UNIT)
    END IF

    record%csv = (csv .OR. is_csvfile(filepath))

#ifdef MPIIO
    IF (mpiio .AND. .NOT. record%csv) THEN
       CALL mpi_barrier(comm, ierr)
       CALL mpi_file_open(comm, trim(record%filepath), MPI_MODE_CREATE+MPI_MODE_APPEND+MPI_MODE_WRONLY, MPI_INFO_NULL, record%filehandle, ierr)
       CALL mpi_file_set_atomicity(record%filehandle, .TRUE., ierr)
    END IF
#endif

    IF (rank==0) THEN
       WRITE(REPORT_UNIT, '(A)', ADVANCE='NO') "'"//trim(record%filepath)//"'" // ' is activated'
       IF (record%csv) THEN
          WRITE(REPORT_UNIT, '(A)') ' (csv-text)'
       ELSE
          WRITE(REPORT_UNIT, '(A)') ' (binary)'
       END IF
    END IF

    record%start = datetime_seconds(start)
    IF (end(1:1)=='+') THEN
       record%end = record%start + interval_seconds(end(2:)//' ')
    ELSE
       record%end = datetime_seconds(end)
    END IF

    IF (interval /= "") record%interval = interval_seconds(interval)

    record%lastsync = record%start

  END SUBROUTINE init_record

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE close_record(record)
    TYPE(record_struct), INTENT(INOUT) :: record

#ifdef PARALLEL_MPI
    INTEGER :: ierr
#endif

    IF (.NOT. record%initialized) RETURN

    CALL flush_record(record, all=.TRUE.)

#ifdef MPIIO
    IF (mpiio) CALL mpi_file_close(record%filehandle, ierr)
#endif

#ifdef NO_ALLOCATABLE_IN_TYPE
    IF (ASSOCIATED(record%buffer)) DEALLOCATE(record%buffer)
#else
    IF (ALLOCATED(record%buffer))  DEALLOCATE(record%buffer)
#endif

    record%initialized = .FALSE.

  END SUBROUTINE close_record

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE flush_record(record, all, timestamp)
    TYPE(record_struct), INTENT(INOUT) :: record
    LOGICAL,   OPTIONAL, INTENT(IN)    :: all
    LOGICAL,   OPTIONAL, INTENT(IN)    :: timestamp

    TYPE(particle_struct) :: tmp_particle
    CHARACTER(160) :: line
    LOGICAL :: all_, ts_
    INTEGER :: l, n
    INTEGER :: iostat
#ifdef MPIIO
    INTEGER :: ierr
    INTEGER(MPI_OFFSET_KIND) :: offset
#endif

    IF (.NOT. record%initialized) RETURN

    all_ = .FALSE.
    IF (present(all)) all_ = all

    ts_ = .FALSE.
    IF (present(timestamp)) ts_ = timestamp
    IF (all_) ts_ = ts_ .AND. rank==npes-1

    IF (ts_) CALL append_record(get_timestamp(t_current), record)

    IF (record%csv) THEN
       IF (all_) CALL serial_begin

       IF (record%n_write > 0) THEN
          OPEN(UNIT     = TMP_UNIT,        &
               FILE     = record%filepath, &
               FORM     = 'FORMATTED',     &
               ACCESS   = 'SEQUENTIAL',    &
               STATUS   = 'UNKNOWN',       &
               POSITION = 'APPEND',        &
               ACTION   = 'WRITE',         &
               IOSTAT   = iostat)

          CALL assert(iostat==0, "failed to open '"//trim(record%filepath)//"'")

          DO n=1, record%n_write
             line = repeat(" ", 160)
             CALL deserialize(record%buffer(:,n), tmp_particle, convert_bigendian=.TRUE.)
             IF (tmp_particle%category==cat_timestamp) THEN
                line = format_timestamp(tmp_particle)
                WRITE(TMP_UNIT, '(A)') trim(line)
             END IF

             line = format_particle(tmp_particle)
             WRITE(TMP_UNIT, '(A)') trim(line)
          END DO

          CLOSE(TMP_UNIT)
       END IF

       IF (all_) CALL serial_end
    ELSE
#ifdef MPIIO
       IF (mpiio) THEN
       IF (.NOT. all_ .OR. rank==0) CALL mpi_file_get_size(record%filehandle, offset, ierr)

       IF (all_) THEN
          IF (rank > 0)      CALL mpi_recv(offset,                                 1, MPI_INTEGER8, rank-1, 0, comm, MPI_STATUS_IGNORE, ierr)
          IF (rank < npes-1) CALL mpi_send(offset+record%n_write*size_of_particle, 1, MPI_INTEGER8, rank+1, 0, comm, ierr)
       END IF

       IF (record%n_write > 0) CALL mpi_file_write_at(record%filehandle, offset, record%buffer, record%n_write*size_of_particle, MPI_BYTE, MPI_STATUS_IGNORE, ierr)
       ELSE
#endif
       IF (all_) CALL serial_begin

       IF (record%n_write > 0) THEN
          OPEN(UNIT     = TMP_UNIT,        &
               FILE     = record%filepath, &
               FORM     = 'UNFORMATTED',   &
               ACCESS   = 'STREAM',        &
               STATUS   = 'UNKNOWN',       &
               POSITION = 'APPEND',        &
               ACTION   = 'WRITE',         &
               IOSTAT   = iostat)

          CALL assert(iostat==0, "failed to open '"//trim(record%filepath)//"'")
          WRITE(TMP_UNIT) record%buffer(:,1:record%n_write)
          CLOSE(TMP_UNIT)
       END IF

       IF (all_) CALL serial_end
#ifdef MPIIO
       END IF
#endif
    END IF

    record%count   = record%count + record%n_write
    record%n_write = 0

  END SUBROUTINE flush_record

!-----------------------------------------------------------------------------------------------------------------------

  INTEGER PURE FUNCTION get_track_index(id)
    INTEGER(8), INTENT(IN) :: id

    INTEGER :: i

    get_track_index = 0

    IF (id < min_track_id .OR. id > max_track_id) RETURN

    DO i=1, n_tracks
       IF (track_id(i) == id) THEN
          get_track_index = i
          EXIT
       END IF
    END DO

  END FUNCTION get_track_index

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE write_particles(filename)
    CHARACTER(*), INTENT(IN) :: filename

    INTEGER(1), ALLOCATABLE :: iobuffer(:,:)

    INTEGER :: i, n
    INTEGER :: cat, ptr
    INTEGER :: iostat, ierr

#ifdef MPIIO
    INTEGER    :: fh
    INTEGER(8) :: offset

    INTEGER(1) :: tmpbuf(size_of_particle)
    LOGICAL    :: check
#endif

    n = sum(particle_cnt(1:isize,1:jsize,1:ksize,1:n_categories))
    IF (n > 0) THEN
       ALLOCATE(iobuffer(size_of_particle,n))

       n = 0
       DO i=1, max_particles
          IF (p_category(i) /= 0) THEN
             n = n+1
             CALL serialize(get_particle(i), iobuffer(:,n), convert_gridpos=.TRUE., convert_categorycode=.TRUE., convert_bigendian=.TRUE.)
          END IF
       END DO
       CALL assert(n==size(iobuffer,2), "particle-list has been corrupted")
    END IF

#ifdef MPIIO
    IF (mpiio) THEN
    IF (rank==0) THEN
       ! delete garbage if file exists
       OPEN(TMP_UNIT, FILE=trim(filename), FORM='UNFORMATTED', ACCESS='STREAM', STATUS='REPLACE', ACTION='WRITE', IOSTAT=iostat)
       CALL assert(iostat==0, "failed to create file '"//trim(filename)//"'")
       CLOSE(TMP_UNIT)
    END IF
    CALL mpi_barrier(comm, ierr)

    CALL mpi_file_open(comm, trim(filename), MPI_MODE_CREATE + MPI_MODE_WRONLY, MPI_INFO_NULL, fh, ierr)

    offset = 0
    IF (rank > 0)      CALL mpi_recv(offset,   1, MPI_INTEGER8, rank-1, 0, comm, MPI_STATUS_IGNORE, ierr)
    IF (rank < npes-1) CALL mpi_send(offset+n, 1, MPI_INTEGER8, rank+1, 0, comm, ierr)

    IF (n > 0) CALL mpi_file_write_at(fh, INT(offset*size_of_particle, MPI_OFFSET_KIND), iobuffer, n*size_of_particle, MPI_BYTE, MPI_STATUS_IGNORE, ierr)

    CALL mpi_barrier(comm, ierr)

    CALL mpi_file_close(fh, ierr)

#ifdef DEBUG
#ifdef A64FX
    !******* CHECK CODE for particle lost on the Wisteria system (2021/12/7) **********

    CALL mpi_file_open(comm, trim(filename), MPI_MODE_RDONLY, MPI_INFO_NULL, fh, ierr)
    IF (n > 0) THEN
       CALL mpi_file_read_at(fh, INT(offset*size_of_particle, MPI_OFFSET_KIND), tmpbuf, size_of_particle, MPI_BYTE, MPI_STATUS_IGNORE, ierr)
       check = .TRUE.
       DO i=1, size_of_particle
          check = check .AND. (tmpbuf(i) == iobuffer(i,1))
       END DO
       CALL assert(check, "**** WRITE_PARTICLES (MPI-IO) FAILED!!! ****")
    END IF
    CALL mpi_file_close(fh, ierr)

    !**********************************************************************************
#endif
#endif

    ELSE
#endif

    CALL serial_begin

    IF (rank==0) THEN
       OPEN(UNIT     = TMP_UNIT,      &
            FILE     = filename,      &
            FORM     = 'UNFORMATTED', &
            ACCESS   = 'STREAM',      &
            STATUS   = 'REPLACE',     &
            POSITION = 'REWIND',      &
            ACTION   = 'WRITE',       &
            IOSTAT   = iostat)
    ELSE IF (n > 0) THEN
       OPEN(UNIT     = TMP_UNIT,      &
            FILE     = filename,      &
            FORM     = 'UNFORMATTED', &
            ACCESS   = 'STREAM',      &
            STATUS   = 'OLD',         &
            POSITION = 'APPEND',      &
            ACTION   = 'WRITE',       &
            IOSTAT   = iostat)
    ELSE
       iostat = 0
    END IF

    CALL assert(iostat==0, "failed to create file '"//trim(filename)//"'")

    IF (n > 0) WRITE(TMP_UNIT) iobuffer

    IF (rank==0 .OR. n > 0) CLOSE(TMP_UNIT)

    CALL serial_end

#ifdef MPIIO
    END IF
#endif

    IF (allocated(iobuffer)) DEALLOCATE(iobuffer)

    IF (rank==0) WRITE(REPORT_UNIT, *) 'write particles to ' // trim(filename)

  END SUBROUTINE write_particles

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE read_particles(filename, csv, default_category, offset_id)
    CHARACTER(*), INTENT(IN) :: filename
    LOGICAL,      INTENT(IN), OPTIONAL :: csv
    INTEGER,      INTENT(IN), OPTIONAL :: default_category
    INTEGER(8),   INTENT(IN), OPTIONAL :: offset_id

    TYPE(particle_struct) :: tmp_particle

    INTEGER(1) :: iobuffer(size_of_particle, chunksize)
    INTEGER :: nread
    INTEGER :: iostat
    INTEGER :: n, m
    INTEGER :: ptr

    LOGICAL :: csv_

    INTEGER :: default_category_
    INTEGER(8) :: offset_id_

#ifdef PARALLEL_MPI
    INTEGER :: fh

    INTEGER :: mpistat(MPI_STATUS_SIZE)
    INTEGER :: ierr
#endif

    IF (trim(basename(filename)) == '$NONE' .OR. trim(basename(filename))=='$none') RETURN

    csv_ = is_csvfile(filename)
    IF (present(csv)) csv_ = csv_ .OR. csv

    default_category_ = 0
    IF (present(default_category)) default_category_ = default_category

    offset_id_ = 0_8
    IF (present(offset_id)) offset_id_ = offset_id

    iobuffer(:,:) = 0_1

    IF (csv_) THEN
       IF (rank==0) THEN
          OPEN(UNIT   = TMP_UNIT,         &
               FORM   = 'UNFORMATTED',    &
               ACCESS = 'DIRECT',         &
               STATUS = 'SCRATCH',        &
               ACTION = 'READWRITE',      &
               RECL   = size_of_particle)

          OPEN(UNIT   = TMP_UNIT+1,       &
               FILE   = filename,         &
               FORM   = 'FORMATTED',      &
               ACCESS = 'SEQUENTIAL',     &
               STATUS = 'OLD',            &
               ACTION = 'READ',           &
               IOSTAT = iostat)

          CALL assert(iostat==0, "failed to open '"//trim(filename)//"'")

          CALL encode_particles(TMP_UNIT+1, TMP_UNIT)
          CLOSE(TMP_UNIT+1)
       END IF

       m = 0
       DO
          IF (rank==0) THEN
             nread = 0
             DO WHILE (nread < chunksize)
                READ(TMP_UNIT, REC=m*chunksize + nread+1, IOSTAT=iostat) iobuffer(:,nread+1)
                IF (iostat /= 0) EXIT
                nread = nread+1
             END DO
          END IF

          CALL bcast(nread)
          CALL bcast(iobuffer(:,1:chunksize))

          DO n=1, nread
             CALL restore_particle(iobuffer(:,n), tmp_particle, default_category=default_category_)

             IF (tmp_particle%id > 0) tmp_particle%id = tmp_particle%id + offset_id_

             IF (tmp_particle%id < 0)      CYCLE
             IF (tmp_particle%category==0) CYCLE

             CALL create_particle(tmp_particle, record=.FALSE., ptr=ptr)

             IF (ptr /= NULL) max_id = max(max_id, p_id(ptr))
          END DO

          IF (nread /= chunksize) EXIT
          m = m+1
       END DO
       IF (rank==0) CLOSE(TMP_UNIT)
    ELSE
#ifdef MPIIO
       IF (mpiio) THEN
       CALL mpi_file_open(comm, trim(filename), MPI_MODE_RDONLY, MPI_INFO_NULL, fh, ierr)
       CALL assert(ierr==MPI_SUCCESS, "failed to open '"//trim(filename)//"'")

       m = 0
       DO
          CALL mpi_file_read(fh, iobuffer, chunksize*size_of_particle, MPI_BYTE, mpistat, ierr)
          CALL mpi_get_count(mpistat, MPI_BYTE, nread, ierr)

          nread = nread / size_of_particle

          DO n=1, nread
             CALL restore_particle(iobuffer(:,n), tmp_particle, default_category=default_category_)

             IF (tmp_particle%id > 0) tmp_particle%id = tmp_particle%id + offset_id_

             IF (tmp_particle%id < 0)      CYCLE
             IF (tmp_particle%category==0) CYCLE

             CALL create_particle(tmp_particle, record=.FALSE., ptr=ptr)

             IF (ptr /= NULL) max_id = max(max_id, p_id(ptr))
          END DO
          IF (nread < chunksize) EXIT
          m = m + 1
       END DO

       CALL mpi_file_close(fh, ierr)
       ELSE
#endif
       OPEN(UNIT   = TMP_UNIT,         &
            FILE   = filename,         &
            FORM   = 'UNFORMATTED',    &
            ACCESS = 'DIRECT',         &
            STATUS = 'OLD',            &
            ACTION = 'READ',           &
            RECL   = size_of_particle, &
            IOSTAT = iostat)

       CALL assert(iostat==0, "failed to open '"//trim(filename)//"'")

       n = 0
       DO
          n = n+1
          READ(TMP_UNIT, REC=n, IOSTAT=iostat) iobuffer(:,1)
          IF (iostat /= 0) EXIT

          CALL restore_particle(iobuffer(:,1), tmp_particle, default_category=default_category_)

          IF (tmp_particle%id > 0) tmp_particle%id = tmp_particle%id + offset_id_

          IF (tmp_particle%id < 0)      CYCLE
          IF (tmp_particle%category==0) CYCLE

          CALL create_particle(tmp_particle, record=.FALSE., ptr=ptr)

          IF (ptr /= NULL) max_id = max(max_id, p_id(ptr))
       END DO

       CLOSE(TMP_UNIT)
#ifdef MPIIO
       END IF
#endif
    END IF

    str_format = format(m*chunksize+nread)
    IF (rank==0) WRITE(REPORT_UNIT, *) 'read '// trim(str_format) // ' particles from ' // "'"//trim(filename)//"'"

    CALL gmax(max_id, all=.TRUE.)
  END SUBROUTINE read_particles

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE restore_particle(buffer, particle, default_category)
    INTEGER(1),            INTENT(IN)    :: buffer(size_of_particle)
    TYPE(particle_struct), INTENT(INOUT) :: particle
    INTEGER,     OPTIONAL, INTENT(IN)    :: default_category

    INTEGER :: n

    CALL deserialize(buffer, particle, convert_gridpos=.TRUE., convert_categorycode=.TRUE., convert_bigendian=.TRUE.)

    IF (particle%category < 0) RETURN

    IF (particle%id < 0) RETURN

    IF (present(default_category)) THEN
       IF (particle%category == 0) particle%category = category_index_bycode(default_category)
    END IF

    IF (particle%category == 0) RETURN

    DO n=1, n_prop
       IF (particle%property(n) == UNDEF) particle%property(n) = category_info(particle%category)%property_default(n)
    END DO

  END SUBROUTINE restore_particle

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE encode_particles(csvunit, binunit)
    INTEGER, INTENT(IN) :: csvunit
    INTEGER, INTENT(IN) :: binunit

    TYPE(particle_struct) :: tmp_particle

    INTEGER(1) :: iobuffer(size_of_particle)

    INTEGER :: l, m, n
    REAL(8) :: x, y, z
    INTEGER :: iostat
    INTEGER :: np
    CHARACTER(640) :: line

    native_bigendian = check_endian()

    tmp_particle%track_index = -1

    l = 0
    n = 1
    DO
       l = l+1

       READ(csvunit, '(A)', IOSTAT=iostat) line
       IF (iostat < 0) EXIT
       CALL assert(iostat == 0, "line "//trim(format(l))//": too long")

       line = adjustl(trim(line))
       IF (line(1:1)=='#')    CYCLE
       IF (len_trim(line)==0) CYCLE

       tmp_particle%property(:) = UNDEF
       tmp_particle%category    = 0
       tmp_particle%seed        = 0

       np = 0
       DO m=1, len(line)
          IF (line(m:m)==",") np = np + 1
       END DO
       np = np - 4

       CALL assert(np >= 0 .AND. np <= n_prop+1, "line "//trim(format(l))//": invalid format")

       IF (np == n_prop+1) THEN
          READ(line, *) tmp_particle%id, x, y, z, tmp_particle%category, tmp_particle%property(1:n_prop), tmp_particle%seed
       ELSE IF (np > 0) THEN
          READ(line, *) tmp_particle%id, x, y, z, tmp_particle%category, tmp_particle%property(1:np)
       ELSE
          READ(line, *) tmp_particle%id, x, y, z, tmp_particle%category
       END IF

       CALL assert(tmp_particle%category >= -1, "negadive category_code is not allowed")

       tmp_particle%ipos = FLOOR(x) + 1
       tmp_particle%jpos = FLOOR(y) + 1
       tmp_particle%kpos = FLOOR(z) + 1
       tmp_particle%xpos = REAL(x - FLOOR(x))
       tmp_particle%ypos = REAL(y - FLOOR(y))
       tmp_particle%zpos = REAL(z - FLOOR(z))

       IF (tmp_particle%xpos < 0) THEN
          tmp_particle%ipos = tmp_particle%ipos - 1
          tmp_particle%xpos = tmp_particle%xpos + 1.0
       END IF

       IF (tmp_particle%ypos < 0) THEN
          tmp_particle%jpos = tmp_particle%jpos - 1
          tmp_particle%ypos = tmp_particle%ypos + 1.0
       END IF

       IF (tmp_particle%zpos < 0) THEN
          tmp_particle%kpos = tmp_particle%kpos - 1
          tmp_particle%zpos = tmp_particle%zpos + 1.0
       END IF

       CALL serialize(tmp_particle, iobuffer, convert_bigendian=.TRUE.)
       WRITE(binunit, REC=n) iobuffer
       n = n+1
    END DO

  END SUBROUTINE encode_particles

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE decode_particles(binunit, csvunit, category, propfmt, calfmt)
    INTEGER, INTENT(IN) :: binunit
    INTEGER, INTENT(IN) :: csvunit
    INTEGER,      INTENT(IN), OPTIONAL :: category
    CHARACTER(*), INTENT(IN), OPTIONAL :: propfmt
    INTEGER,      INTENT(IN), OPTIONAL :: calfmt

    INTEGER(1) :: iobuffer(size_of_particle)
    TYPE(particle_struct) :: tmp_particle

    INTEGER :: n
    INTEGER :: iostat

    CHARACTER(160) :: line

    native_bigendian = check_endian()

    WRITE(csvunit, '(A)') "#ID, x-position, y-position, z-position, category, properties..., random_seed"

    n = 0
    DO
       n = n+1
       READ(binunit, REC=n, IOSTAT=iostat) iobuffer
       IF (iostat /= 0) EXIT

       CALL deserialize(iobuffer, tmp_particle, convert_bigendian=.TRUE.)

       IF (tmp_particle%category==cat_timestamp) THEN
          line = format_timestamp(tmp_particle, calfmt)
          WRITE(csvunit, '(A)') trim(line)
       END IF

       IF (present(category)) THEN
          IF (tmp_particle%category /= category) CYCLE
       END IF

       line = format_particle(tmp_particle, propfmt)
       WRITE(csvunit, '(A)') trim(line)
    END DO

  END SUBROUTINE decode_particles

!-----------------------------------------------------------------------------------------------------------------------

  CHARACTER(160) FUNCTION format_particle(particle, propfmt, convert_gridpos, convert_categorycode)
    TYPE(particle_struct),  INTENT(IN) :: particle
    CHARACTER(*), OPTIONAL, INTENT(IN) :: propfmt
    LOGICAL,      OPTIONAL, INTENT(IN) :: convert_gridpos
    LOGICAL,      OPTIONAL, INTENT(IN) :: convert_categorycode

    INTEGER :: cat, i, j, k
    CHARACTER(16) :: fmt
    CHARACTER(1)  :: np

    IF (present(propfmt)) THEN
       fmt = trim(propfmt)
    ELSE
       fmt = 'ES13.6'
    END IF

    cat = particle%category
    IF (present(convert_categorycode)) THEN
       IF (convert_categorycode) cat = category_info(cat)%code
    END IF

    i = particle%ipos
    j = particle%jpos
    k = particle%kpos
    IF (present(convert_gridpos)) THEN
       IF (convert_gridpos) THEN
          i = icoord*isize + i
          j = jcoord*jsize + j
          k = kcoord*ksize + k
       END IF
    END IF


    np = trim(format(n_prop,'I1'))

    WRITE(format_particle, '(I20,3(",",F12.6),",",I3,'//np//'(","'//trim(fmt)//'),",",I12)') &
            particle%id, REAL(i - 1 + particle%xpos), &
                         REAL(j - 1 + particle%ypos), &
                         REAL(k - 1 + particle%zpos), &
                         cat, particle%property(:), particle%seed

  END function format_particle

!-----------------------------------------------------------------------------------------------------------------------

  CHARACTER(20) FUNCTION format_timestamp(ts, fmt)
    TYPE(particle_struct), INTENT(IN) :: ts
    INTEGER,    INTENT(IN), OPTIONAL :: fmt
    CALL assert(ts%category == cat_timestamp, "FORMAT_TIMESTAMP only accepts TIMESTAMP particle (category=-1)")
    format_timestamp = "#TS " // format_datetime(real(ts%id,8), fmt)

  END FUNCTION format_timestamp

!-----------------------------------------------------------------------------------------------------------------------

  INTEGER PURE FUNCTION particle_handle(i, j, k, cat)
    INTEGER, INTENT(IN) :: i
    INTEGER, INTENT(IN) :: j
    INTEGER, INTENT(IN) :: k
    INTEGER, INTENT(IN) :: cat

    particle_handle = particle_ptr(i, j, k, cat)

  END FUNCTION particle_handle

!-----------------------------------------------------------------------------------------------------------------------

  INTEGER PURE FUNCTION particle_handle_next(p)
    INTEGER, INTENT(IN) :: p

    IF (p == NULL) THEN
       particle_handle_next = NULL
       RETURN
    END IF

    particle_handle_next = p_next(p)

  END FUNCTION particle_handle_next

!-----------------------------------------------------------------------------------------------------------------------

  TYPE(particle_struct) PURE FUNCTION get_particle(p, seed)
    INTEGER, INTENT(IN) :: p
    INTEGER, INTENT(IN), OPTIONAL :: seed

    get_particle%id   = p_id(p)
    get_particle%ipos = p_ijk(1,p)
    get_particle%jpos = p_ijk(2,p)
    get_particle%kpos = p_ijk(3,p)
    get_particle%xpos = p_xyz(1,p)
    get_particle%ypos = p_xyz(2,p)
    get_particle%zpos = p_xyz(3,p)
    get_particle%property(:)  = p_property(:,p)
    get_particle%category     = p_category(p)

    IF (present(seed)) THEN
       get_particle%seed = seed
    ELSE
       get_particle%seed = p_seed(p)
    END IF

    get_particle%track_index  = p_track_index(p)

  END FUNCTION get_particle

!-----------------------------------------------------------------------------------------------------------------------

  INTEGER(8) FUNCTION get_particle_id(p)
    INTEGER, INTENT(IN) :: p

    get_particle_id = p_id(p)
  END FUNCTION get_particle_id

!-----------------------------------------------------------------------------------------------------------------------

  REAL(4) FUNCTION get_particle_xpos(p)
    INTEGER, INTENT(IN) :: p

    get_particle_xpos = p_xyz(1,p)
  END FUNCTION get_particle_xpos

!-----------------------------------------------------------------------------------------------------------------------

  REAL(4) FUNCTION get_particle_ypos(p)
    INTEGER, INTENT(IN) :: p

    get_particle_ypos = p_xyz(2,p)
  END FUNCTION get_particle_ypos

!-----------------------------------------------------------------------------------------------------------------------

  REAL(4) FUNCTION get_particle_zpos(p)
    INTEGER, INTENT(IN) :: p

    get_particle_zpos = p_xyz(3,p)
  END FUNCTION get_particle_zpos

!-----------------------------------------------------------------------------------------------------------------------

  REAL(PROP_KIND) PURE FUNCTION get_particle_property(p, prop)
    INTEGER, INTENT(IN) :: p
    INTEGER, INTENT(IN) :: prop

    get_particle_property = p_property(prop,p)
  END FUNCTION get_particle_property

!-----------------------------------------------------------------------------------------------------------------------

  REAL(PROP_KIND) PURE FUNCTION get_default_property(cat, prop)
    INTEGER, INTENT(IN) :: cat
    INTEGER, INTENT(IN) :: prop

    get_default_property = category_info(cat)%property_default(prop)
  END FUNCTION get_default_property

!-----------------------------------------------------------------------------------------------------------------------

  INTEGER PURE FUNCTION get_particle_category(p)
    INTEGER, INTENT(IN) :: p

    get_particle_category = p_category(p)
  END FUNCTION get_particle_category

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE set_particle_property(p, prop, value)
    INTEGER, INTENT(IN) :: p
    INTEGER, INTENT(IN) :: prop
    REAL(PROP_KIND), INTENT(IN) :: value

    p_property(prop,p) = value

  END SUBROUTINE set_particle_property

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE set_particle_pos_in_grid(p, xpos, ypos, zpos)
    INTEGER, INTENT(IN) :: p
    REAL(4), INTENT(IN) :: xpos
    REAL(4), INTENT(IN) :: ypos
    REAL(4), INTENT(IN) :: zpos

#ifdef DEBUG__
    CALL assert(xpos >= 0.0 .AND. xpos <= 1.0)
    CALL assert(ypos >= 0.0 .AND. ypos <= 1.0)
    CALL assert(zpos >= 0.0 .AND. zpos <= 1.0)
#endif

    p_xyz(1,p) = xpos
    p_xyz(2,p) = ypos
    p_xyz(3,p) = zpos

  END SUBROUTINE set_particle_pos_in_grid

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE set_particle_xpos(p, xpos)
    INTEGER, INTENT(IN) :: p
    REAL(4), INTENT(IN) :: xpos

#ifdef DEBUG__
    CALL assert(xpos >= 0.0 .AND. xpos <= 1.0)
#endif

    p_xyz(1,p) = xpos

  END SUBROUTINE set_particle_xpos

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE set_particle_ypos(p, ypos)
    INTEGER, INTENT(IN) :: p
    REAL(4), INTENT(IN) :: ypos

#ifdef DEBUG__
    CALL assert(ypos >= 0.0 .AND. ypos <= 1.0)
#endif

    p_xyz(2,p) = ypos

  END SUBROUTINE set_particle_ypos

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE set_particle_zpos(p, zpos)
    INTEGER, INTENT(IN) :: p
    REAL(4), INTENT(IN) :: zpos

#ifdef DEBUG__
    CALL assert(zpos >= 0.0 .AND. zpos <= 1.0)
#endif

    p_xyz(3,p) = zpos

  END SUBROUTINE set_particle_zpos

!-----------------------------------------------------------------------------------------------------------------------

  INTEGER PURE FUNCTION category_index_byname(name)
    CHARACTER(*), INTENT(IN) :: name

    INTEGER :: n

    category_index_byname = 0
    IF (name == '') RETURN

    IF (name == 'TIMESTAMP') THEN
       category_index_byname = cat_timestamp
       RETURN
    END IF

    DO n=1, n_categories
       IF (trim(name) == trim(category_info(n)%name)) THEN
          category_index_byname = n
          RETURN
       END IF
    END DO

  END FUNCTION category_index_byname

!-----------------------------------------------------------------------------------------------------------------------

  INTEGER PURE FUNCTION category_index_bycode(code)
    INTEGER, INTENT(IN) :: code

    INTEGER :: n

    IF (code < 0) THEN
       category_index_bycode = code
       RETURN
    END IF

    category_index_bycode = 0

    DO n=1, n_categories
       IF (code == category_info(n)%code) THEN
          category_index_bycode = n
          RETURN
       END IF
    END DO

  END FUNCTION category_index_bycode

!-----------------------------------------------------------------------------------------------------------------------

  INTEGER PURE FUNCTION property_index(category, property_name)
    INTEGER,      INTENT(IN) :: category
    CHARACTER(*), INTENT(IN) :: property_name

    INTEGER :: n

    property_index = 0

    IF (property_name == '') RETURN

    DO n=1, n_prop
       IF (trim(property_name) == trim(category_info(category)%property_name(n))) THEN
          property_index = n
          RETURN
       END IF
    END DO

  END FUNCTION property_index

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE particle_driven_forcing(gx, gy, gz)
    REAL(8), INTENT(INOUT) :: gx(0:isize, 1:jsize, 1:ksize)
    REAL(8), INTENT(INOUT) :: gy(1:isize, 0:jsize, 1:ksize)
    REAL(8), INTENT(INOUT) :: gz(1:isize, 1:jsize, 0:ksize)

    CALL particles_buoyancy(gz)

  END SUBROUTINE particle_driven_forcing

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE input(info)
    TYPE(input_info_struct), INTENT(INOUT) :: info
    INTEGER :: m, n
    INTEGER :: ptr

    TYPE(particle_struct) :: tmp_particle

    INTEGER :: count, count_total

    INTEGER :: itmp
    REAL(8) :: rtmp

    IF (t_current <  info%start) RETURN

    IF (t_current >= info%end) THEN
       IF (info%initialized) CALL finalize_particle_input(info)
       RETURN
    END IF

    count = 0

    IF (info%mode == 10) THEN
       ! HISTORICAL mode specific procedure
       IF (info%initialized) THEN
          IF (t_current - info%lastread >= info%interval) THEN
             info%lastread = info%lastread + info%interval
             CALL finalize_particle_input(info)

             itmp = index(info%filepath, '.', back=.TRUE.)
             info%filepath = info%filepath(1:itmp) // format_datetime(info%lastread)
          END IF
       ELSE
          rtmp = t_current - info%start
          itmp = int(rtmp/info%interval)

          info%lastread = info%start + info%interval * itmp

          info%filepath = trim(info%filepath) // '.' // format_datetime(info%lastread)
       END IF
    END IF

    IF (.NOT. info%initialized) CALL init_particle_input(info)

    IF (.NOT. info%cyclic .AND. info%count >= info%length) RETURN
    IF (mod(t_current - info%start,info%interval) < dtime) info%residual = 0.0

    IF (info%mode == 20) THEN
       DO m=info%count+1, info%length
          IF (info%particles(m)%category /= cat_timestamp) CYCLE
          IF (info%particles(m)%id >= t_current + dtime) EXIT
       END DO

       info%residual = m - 1 - info%count
    ELSE
       info%residual = info%residual + get_input_count(info)
    END IF

    DO WHILE (info%residual >= 1.0)
       IF (info%count >= info%length .AND. .NOT. info%cyclic) EXIT

       info%count    = info%count    + 1
       info%residual = info%residual - 1.0

       m = mod(info%count-1, info%length)+1

       IF (info%particles(m)%category <= 0) CYCLE
       IF (info%particles(m)%id       <  0) CYCLE

       tmp_particle = info%particles(m)

       IF (info%cyclic .AND. info%count > info%length .AND. tmp_particle%id /= 0) THEN
          IF (info%cyclic_id==0) THEN
             tmp_particle%id = 0
          ELSE IF (info%cyclic_id < 0) THEN
             tmp_particle%id = tmp_particle%id + INT(info%count/info%length,KIND=8)*(-info%cyclic_id*info%length)
          ELSE
             tmp_particle%id = tmp_particle%id + INT(info%count/info%length,KIND=8)*info%cyclic_id
          END IF
       END IF

       CALL create_particle(tmp_particle, ptr=ptr)

       IF (ptr /= NULL) count = count + 1

#ifdef DEBUG
       IF (ptr /= NULL) THEN
          WRITE(REPORT_UNIT, '(A)', ADVANCE='NO') ' input particle:' // trim(format_particle(tmp_particle, convert_categorycode=.TRUE.))
       END IF
#endif
    END DO

    IF (info%report) THEN
       count_total = count
       CALL barrier !OFP requires this barrier (3/13/2018, Y. Matsumura)
       CALL gsum(count_total)

       str_format = format(count_total)
       IF (rank==0) write(REPORT_UNIT,*) 'particle_input: total ' // trim(str_format) // &
            ' particles from ' // "'"//trim(info%filepath)//"'"
    END IF

  END SUBROUTINE input

!-----------------------------------------------------------------------------------------------------------------------

  REAL(8) PURE FUNCTION get_input_count(info, time)
    TYPE(input_info_struct), INTENT(IN) :: info
    REAL(8),                 INTENT(IN), OPTIONAL :: time

    REAL(8) :: t_
    REAL(4) :: alpha

    INTEGER :: itmp
    REAL(8) :: rtmp

    t_ = t_current
    IF (PRESENT(time)) t_ = time

    get_input_count = 0.0
    IF (t_ < info%start .OR. t_ >= info%end) RETURN

    rtmp = t_ - info%start
    itmp = int(rtmp/info%interval)

    SELECT CASE (info%mode)
    CASE(2)
       get_input_count = info%rate(mod(itmp,info%periods)) * dtime
    CASE(3)
       alpha = REAL(mod(rtmp, info%interval)/info%interval)

       get_input_count = ((1.0-alpha)*info%rate(mod(itmp,  info%periods)) &
                              +alpha *info%rate(mod(itmp+1,info%periods))) * dtime
    CASE(4)
       IF (rtmp - itmp*info%interval < dtime) THEN
          get_input_count = info%rate(mod(itmp,info%periods))
       ELSE
          get_input_count = 0.0
       END IF

    CASE(5)
       get_input_count = max(0.0, info%rate(0) + info%rate(1)*REAL(sin(2.0*pi*rtmp/info%interval))) * dtime

    CASE(10)
       get_input_count = info%rate(0)*dtime

    END SELECT

    IF (info%rate_rel) THEN
       get_input_count = get_input_count * info%length
    END IF

  END FUNCTION get_input_count

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE particles_buoyancy(gz)
    REAL(8), INTENT(INOUT) :: gz(1:isize, 1:jsize, 0:ksize)

    REAL(8) :: b_l(isize, jsize, 1:ksize+1)
    REAL(8) :: b_u(isize, jsize, 0:ksize)
    REAL(8) :: bfactor(isize, jsize, ksize)
    INTEGER :: p

    REAL(8) :: zpos, btmp

    REAL(8) :: gp

    INTEGER :: i, j, k, ijk, cat
    INTEGER :: ierr

    LOGICAL :: stat

!$OMP PARALLEL WORKSHARE
    b_l(:,:,:) = 0.0
    b_u(:,:,:) = 0.0
!$OMP END PARALLEL WORKSHARE

    DO cat=1, n_categories
       IF (.NOT. category_info(cat)%buoyancy)    CYCLE
       IF (category_info(cat)%propindex_mass==0) CYCLE

       CALL checkin('BFACTOR_' // trim(category_info(cat)%name), bfactor, stat)
       IF (.NOT. stat) bfactor(:,:,:) = 1.0

       gp = gravity * (rho_0 - category_info(cat)%rho) / category_info(cat)%rho

!$OMP PARALLEL DO SCHEDULE(dynamic, ompchunk) PRIVATE(i, j, k, p, zpos, btmp)
       DO ijk=0, isize*jsize*ksize-1
          k = ijk/(isize*jsize) + 1
          j = (ijk - (k-1)*isize*jsize)/isize + 1
          i = ijk - (k-1)*isize*jsize - (j-1)*isize + 1

          p = particle_handle(i,j,k,cat)
          DO WHILE (p /= NULL)
             zpos = get_particle_zpos(p)
             btmp = get_particle_property(p, category_info(cat)%propindex_mass) * gp * bfactor(i,j,k)

             b_l(i,j,k) = b_l(i,j,k) + btmp*(1.0-zpos)
             b_u(i,j,k) = b_u(i,j,k) + btmp*zpos

             p = particle_handle_next(p)
          END DO
       END DO
    END DO

#ifdef PARALLEL3D
    CALL mpi_sendrecv(b_l(1,1,1),       isize*jsize, MPI_REAL8, rank_l, 0, &
                      b_l(1,1,ksize+1), isize*jsize, MPI_REAL8, rank_u, 0, comm, MPI_STATUS_IGNORE, ierr)

    CALL mpi_sendrecv(b_u(1,1,ksize),   isize*jsize, MPI_REAL8, rank_u, 1, &
                      b_u(1,1,0),       isize*jsize, MPI_REAL8, rank_l, 1, comm, MPI_STATUS_IGNORE, ierr)
#endif

!$OMP PARALLEL DO
    DO k=0, ksize
    DO j=1, jsize
    DO i=1, isize
       gz(i,j,k) = gz(i,j,k) + imask3d(i,j,k)*imask3d(i,j,k+1)*REAL(b_l(i,j,k+1)+b_u(i,j,k),4)*2.0/((dvol(i,j,k)+dvol(i,j,k+1))*rho_0)
    END DO
    END DO
    END DO

  END SUBROUTINE particles_buoyancy

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE deploy_npzd
    USE velocity,   ONLY: u, v
    USE tracers
    USE npzd

    INTEGER :: i, j, k, ijk
    INTEGER :: cat
    INTEGER :: n, nn, s
    REAL(4) :: r
    INTEGER :: seed
    LOGICAL :: stat
    INTEGER :: num(0:3)
    INTEGER :: ptr

    TYPE(particle_struct) :: tmp_particle

    INTEGER :: tracer_index_npzd(0:3)

    INTEGER(1) :: flag(isize,jsize,ksize)

    INTEGER(1), PARAMETER :: fmask(0:3) = (/1_1, 2_1, 4_1, 8_1/)

    cat = category_index('NPZD')
    IF (cat==0) RETURN

    CALL checkin("PNPZD_FLAG", flag, stat=stat)
    IF (.NOT. stat) RETURN

    tracer_index_npzd(0) = tracer_index_nut
    tracer_index_npzd(1) = tracer_index_phy
    tracer_index_npzd(2) = tracer_index_zoo
    tracer_index_npzd(3) = tracer_index_det

    tmp_particle%id   = 0_8
    tmp_particle%category = cat
    tmp_particle%xpos = UNDEF
    tmp_particle%ypos = UNDEF
    tmp_particle%zpos = UNDEF
    tmp_particle%seed = 0
    tmp_particle%property(2) = 0.0
    tmp_particle%property(3) = 0.0
    tmp_particle%property(4:) = UNDEF
    tmp_particle%track_index = -1

    seed = mkseed(i4=kcoord*ipes*jpes+jcoord*ipes+icoord, r8=t_current)

    IF (stat) THEN
       DO k=1, ksize
       DO j=1, jsize
       DO i=1, isize
          IF (flag(i,j,k)==0 .OR. .NOT. lmask3d(i,j,k)) CYCLE

          num(:) = 0
          IF (particle_cnt(i,j,k,cat) > 0) THEN
             ptr = particle_ptr(i,j,k,cat)
             DO WHILE (ptr /= NULL)
                s = INT(p_property(1,ptr))
                num(s) = num(s)+1
                ptr = p_next(ptr)
             END DO
          END IF

          tmp_particle%ipos = i
          tmp_particle%jpos = j
          tmp_particle%kpos = k

          tmp_particle%property(3) = 0.0

          DO s=0, 3
             IF (iand(flag(i,j,k), fmask(s)) == 0_1) CYCLE

             tmp_particle%property(1) = REAL(s)
             IF (s==3) tmp_particle%property(3) = -PON_SVn*(1.0 + npzd_svsigma*nrand(seed))

             r = max(tracer(i,j,k,tracer_index_npzd(s))*dvol(i,j,k)/npzd_permoln - num(s), 0.0) * npzd_dfactor
             nn = INT(r)
             r = r - nn
             IF (abs(urand(seed)) < r) nn = nn + 1
             DO n=1, nn
                CALL create_particle(tmp_particle)
             END DO
          END DO
       END DO
       END DO
       END DO
    END IF

  END SUBROUTINE deploy_npzd

  SUBROUTINE deploy_npzd_boundary
    USE velocity,   ONLY: u, v
    USE tracers
    USE npzd

    INTEGER :: i, j, k, ijk
    INTEGER :: p, s
    INTEGER :: cat
    INTEGER :: n, nn
    REAL(4) :: r
    INTEGER :: seed

    TYPE(particle_struct) :: tmp_particle

    INTEGER :: tracer_index_npzd(0:3)

    cat = category_index('NPZD')
    IF (cat==0) RETURN

    tracer_index_npzd(0) = tracer_index_nut
    tracer_index_npzd(1) = tracer_index_phy
    tracer_index_npzd(2) = tracer_index_zoo
    tracer_index_npzd(3) = tracer_index_det

    tmp_particle%id   = 0_8
    tmp_particle%category = cat
    tmp_particle%xpos = UNDEF
    tmp_particle%ypos = UNDEF
    tmp_particle%zpos = UNDEF
    tmp_particle%seed = 0
    tmp_particle%property(2) = 0.0
    tmp_particle%property(3) = 0.0
    tmp_particle%property(4:) = UNDEF
    tmp_particle%track_index = -1

    seed = mkseed(i4=kcoord*ipes*jpes+jcoord*ipes+icoord, r8=t_current)

    IF (open_w .AND. icoord==0) THEN
       DO ijk=0, jsize*ksize-1
          k = ijk/jsize + 1
          j = ijk - (k-1)*isize + 1

          IF (u(1,j,k) <= 0.0)  CYCLE
          IF (open_s .AND. jcoord==0      .AND. j==1)     CYCLE
          IF (open_n .AND. jcoord==jpes-1 .AND. j==jsize) CYCLE

          tmp_particle%ipos = 2
          tmp_particle%jpos = j
          tmp_particle%kpos = k

          CALL deploy(1,j,k,u(1,j,k)*dsx(1,j,k)*dtime)
       END DO
    END IF

    IF (open_e .AND. icoord==ipes-1) THEN
       DO ijk=0, jsize*ksize-1
          k = ijk/jsize + 1
          j = ijk - (k-1)*isize + 1

          IF (u(isize-1,j,k) >= 0.0) CYCLE
          IF (open_s .AND. jcoord==0      .AND. j==1)     CYCLE
          IF (open_n .AND. jcoord==jpes-1 .AND. j==jsize) CYCLE

          tmp_particle%ipos = isize-1
          tmp_particle%jpos = j
          tmp_particle%kpos = k

          CALL deploy(isize,j,k,-u(isize-1,j,k)*dsx(isize-1,j,k)*dtime)
       END DO
    END IF

    IF (open_s .AND. jcoord==0) THEN
       DO ijk=0, isize*ksize-1
          k = ijk/isize + 1
          i = ijk - (k-1)*isize + 1

          IF (v(i,1,k) <= 0.0) CYCLE
          IF (open_w .AND. icoord==0      .AND. i==1)     CYCLE
          IF (open_e .AND. icoord==ipes-1 .AND. i==isize) CYCLE

          tmp_particle%ipos = i
          tmp_particle%jpos = 2
          tmp_particle%kpos = k

          CALL deploy(i,1,k,v(i,1,k)*dsy(i,1,k)*dtime)
       END DO
    END IF

    IF (open_n .AND. jcoord==jpes-1) THEN
       DO ijk=0, isize*ksize-1
          k = ijk/isize + 1
          i = ijk - (k-1)*isize + 1

          IF (v(i,jsize-1,k) >= 0.0)  CYCLE
          IF (open_w .AND. icoord==0      .AND. i==1)     CYCLE
          IF (open_e .AND. icoord==ipes-1 .AND. i==isize) CYCLE

          tmp_particle%ipos = i
          tmp_particle%jpos = jsize-1
          tmp_particle%kpos = k

          CALL deploy(i,jsize,k,-v(i,jsize-1,k)*dsy(i,jsize-1,k)*dtime)
       END DO
    END IF

  CONTAINS
    SUBROUTINE deploy(i, j, k, dvol)
      INTEGER, INTENT(IN) :: i, j, k
      REAL(8), INTENT(IN) :: dvol

      INTEGER :: s, n, nn
      REAL(4) :: r

      DO s=0, 3
         r = REAL(tracer(i,j,k,tracer_index_npzd(s))*dvol / npzd_permoln)
         nn = INT(r)
         r = r - nn
         IF (abs(urand(seed)) < r) nn = nn + 1
         tmp_particle%property(1) = REAL(s)
         IF (s==3) THEN
            tmp_particle%property(3) = -PON_SVn*(1.0 + npzd_svsigma*nrand(seed))
         ELSE
            tmp_particle%property(3) = 0.0
         END IF
         DO n=1, nn
            CALL create_particle(tmp_particle)
         END DO
      END DO
    END SUBROUTINE deploy

  END SUBROUTINE deploy_npzd_boundary

!----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE init_npzd_particles
    USE tracers
    USE npzd

    INTEGER :: i, j, k, ijk
    INTEGER :: p, s
    INTEGER :: cat
    INTEGER :: n, nn
    REAL(4) :: r

    INTEGER :: seed

    TYPE(particle_struct) :: tmp_particle

    INTEGER :: tracer_index_npzd(0:3)

    INTEGER :: ncreate(0:3)

    cat = category_index('NPZD')
    IF (cat==0) RETURN

    tracer_index_npzd(0) = tracer_index_nut
    tracer_index_npzd(1) = tracer_index_phy
    tracer_index_npzd(2) = tracer_index_zoo
    tracer_index_npzd(3) = tracer_index_det

!$OMP PARALLEL PRIVATE(ijk, i, j, k, n, nn, s, r, p, seed, tmp_particle)
    tmp_particle%id   = 0_8
    tmp_particle%category = cat
    tmp_particle%xpos = UNDEF
    tmp_particle%ypos = UNDEF
    tmp_particle%zpos = UNDEF
    tmp_particle%seed = 0
    tmp_particle%property(2) = 0.0
    tmp_particle%property(3) = 0.0
    tmp_particle%property(4:) = UNDEF
    tmp_particle%track_index = -1

    ncreate(:) = 0
!$OMP DO SCHEDULE(dynamic, ompchunk)
    DO ijk=0, isize*jsize*ksize-1
       k = ijk/(isize*jsize) + 1
       j = (ijk - (k-1)*(isize*jsize))/isize + 1
       i = ijk - (k-1)*(isize*jsize) - (j-1)*isize + 1

       IF (.NOT. lmask3d(i,j,k)) CYCLE

       seed = mkseed(i8=1_8*(ksize*kcoord+k)*dimx*dimy + (jsize*jcoord+j)*dimx + (isize*icoord+i), r8=t_current)
       tmp_particle%ipos = i
       tmp_particle%jpos = j
       tmp_particle%kpos = k

       DO s=0, 3
          r = max(tracer(i,j,k,tracer_index_npzd(s))*dvol(i,j,k) / npzd_permoln, 0.0)
          nn = INT(r)
          r = r - nn
          IF (abs(urand(seed)) < r) nn = nn + 1
          tmp_particle%property(1) = REAL(s)
          IF (s==3) THEN
             tmp_particle%property(3) = -PON_SVn*(1.0 + npzd_svsigma*nrand(seed))
          ELSE
             tmp_particle%property(3) = 0.0
          END IF

          DO n=1, nn
!$OMP CRITICAL
             CALL create_particle(tmp_particle)
             ncreate(s) = ncreate(s) + 1
!$OMP END CRITICAL
          END DO
       END DO
    END DO
!$OMP END PARALLEL

    CALL gsum(ncreate)

    IF (rank==0) THEN
       WRITE(REPORT_UNIT,'(A)')     "Initialize NPZD particles from corresponding tracers"
       WRITE(REPORT_UNIT,'(A,I10)') "   N: ", ncreate(0)
       WRITE(REPORT_UNIT,'(A,I10)') "   P: ", ncreate(1)
       WRITE(REPORT_UNIT,'(A,I10)') "   Z: ", ncreate(2)
       WRITE(REPORT_UNIT,'(A,I10)') "   D: ", ncreate(3)
    END IF

  END SUBROUTINE init_npzd_particles

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE particle_histogram(cat, prop, bins, result, all)
    INTEGER, INTENT(IN) :: cat
    INTEGER, INTENT(IN) :: prop
    REAL(4), INTENT(IN) :: bins(:)
    INTEGER, INTENT(OUT) :: result(0:)
    LOGICAL, INTENT(IN), OPTIONAL :: all

    INTEGER :: nbin
    INTEGER :: ptr
    INTEGER :: ijk, i, j, k, n
    LOGICAL :: all_

    all_ = .TRUE.
    IF (present(all)) all_ = all

    DO nbin=0, size(bins)-1
       IF (bins(nbin+1)==UNDEF) EXIT
    END DO

    result(0:nbin) = 0

!$OMP PARALLEL DO SCHEDULE(dynamic, ompchunk) REDUCTION(+:result) PRIVATE(i, j, k, n, ptr)
    DO ijk=0, isize*jsize*ksize-1
       k = ijk/(isize*jsize) + 1
       j = (ijk - (k-1)*(isize*jsize))/isize + 1
       i = ijk - (k-1)*(isize*jsize) - (j-1)*isize + 1

       IF (particle_cnt(i,j,k,cat) == 0) CYCLE

       ptr = particle_ptr(i,j,k,cat)
       DO WHILE (ptr/=NULL)
          IF (p_property(prop,ptr) == UNDEF) CYCLE
          DO n=1, nbin
             IF (p_property(prop,ptr) < bins(n)) THEN
                result(n-1) = result(n-1) + 1
                EXIT
             END IF
          END DO
          IF (n==nbin+1) result(nbin) = result(nbin) + 1

          ptr = p_next(ptr)
       END DO
    END DO

    CALL gsum(result, all=all_)

  END SUBROUTINE particle_histogram

!-----------------------------------------------------------------------------------------------------------------------

! terminal buoyant rise/fall velocity (relative to the embient vertical velocity) of a sphere of diameter 'd'

  REAL(4) PURE FUNCTION terminal_velocity(d, rho, scheme)
    REAL(4), INTENT(IN) :: d
    REAL(4), INTENT(IN) :: rho
    INTEGER, INTENT(IN), OPTIONAL :: scheme

    INTEGER :: s
    REAL(4) :: w

    s = 1
    IF (present(scheme)) s = scheme

    SELECT CASE (s)
    CASE (1)
       IF (rho > rho_0) THEN
          terminal_velocity = max(terminal_velocity_stokes(d, rho), &
                                  terminal_velocity_allen( d, rho), &
                                  terminal_velocity_newton(d, rho))
       ELSE
          terminal_velocity = min(terminal_velocity_stokes(d, rho), &
                                  terminal_velocity_allen( d, rho), &
                                  terminal_velocity_newton(d, rho))
       END IF

    CASE (2)
       terminal_velocity = terminal_velocity_rubey(d, rho)
    CASE (3)
       terminal_velocity = terminal_velocity_stokes(d, rho)
    CASE (4)
       terminal_velocity = terminal_velocity_allen(d, rho)
    CASE (5)
       terminal_velocity = terminal_velocity_newton(d, rho)
    CASE (6)
       terminal_velocity = terminal_velocity_kaiser(d, rho)
    CASE (7) ! for microplastic of wide range
       IF (rho > rho_0) THEN
          terminal_velocity = max(terminal_velocity_kaiser(d, rho), &
                                  terminal_velocity_allen(  d, rho), &
                                  terminal_velocity_newton(d, rho))
       ELSE
          terminal_velocity = min(terminal_velocity_kaiser(d, rho), &
                                  terminal_velocity_allen(  d, rho), &
                                  terminal_velocity_newton(d, rho))
       END IF
    CASE (8)
       terminal_velocity = terminal_velocity_jsce(d, rho)
    CASE (9)
       terminal_velocity = terminal_velocity_kaiser_mod(d, rho)

    CASE DEFAULT
       terminal_velocity = 0.0
    END SELECT

  END FUNCTION terminal_velocity

  REAL(4) PURE FUNCTION terminal_velocity_default(d, rho)
    ! mix of Stokes', Allen's and Newton's terminal velocity based on the Reynolds number
    !  (currently unused)
    REAL(4), INTENT(IN) :: d
    REAL(4), INTENT(IN) :: rho

    REAL(4) :: gp, w

    REAL(4) :: re1
    REAL(4) :: re2

    re1 = 2
    re2 = 500

    gp = abs((rho - rho_0)/rho_0) * gravity

    w = d**2 * gp/(18*nu_mol) ! Stokes's velocity

    IF (w * d > re1 * nu_mol) THEN ! Re > 2
       w = (4*gp**2 / (225*nu_mol))**(0.3333) * d ! Allen's velocity
       IF (w * d > re2 * nu_mol) THEN ! Re > 500
          w = sqrt(4 * gp * d / (3*0.44)) ! Newton's velocity
       END IF
    END IF

    terminal_velocity_default = sign(w, rho_0 - rho)

  END FUNCTION terminal_velocity_default

  REAL(4) PURE FUNCTION terminal_velocity_jsce(d, rho)
    ! mix of Stokes', Allen's and Newton's terminal velocity based on the Reynolds number
    !  following Japan Socierty of Civil Engeener (1999), used in Isobe et al. (2014)
    REAL(4), INTENT(IN) :: d
    REAL(4), INTENT(IN) :: rho

    REAL(4) :: gp, w

    REAL(4) :: re1
    REAL(4) :: re2

    re1 = 1
    re2 = 100

    gp = abs((rho - rho_0)/rho_0) * gravity

    w = d**2 * gp/(18*nu_mol) ! Stokes's velocity

    IF (w * d > re1 * nu_mol) THEN ! Re > 1.0
       w = 0.223 * (gp**2 /nu_mol)**(0.3333) * d ! Allen's velocity
       IF (w * d > re2 * nu_mol) THEN ! Re > 100
          w = 1.82*sqrt(gp * d) ! Newton's velocity
       END IF
    END IF

    terminal_velocity_jsce = sign(w, rho_0 - rho)

  END FUNCTION terminal_velocity_jsce

  REAL(4) PURE FUNCTION terminal_velocity_stokes(d, rho)
    ! Stokes' terminal velocity proportional to d^2
    REAL(4), INTENT(IN) :: d
    REAL(4), INTENT(IN) :: rho

    REAL(4) :: gp, w

    gp = abs((rho - rho_0)/rho_0) * gravity
    w = d**2 * gp/(18*nu_mol) ! Stokes velocity

    terminal_velocity_stokes = sign(w, rho_0 - rho)

  END FUNCTION terminal_velocity_stokes

  REAL(4) PURE FUNCTION terminal_velocity_allen(d, rho)
    ! Allen's terminal velocity proportional to d
    REAL(4), INTENT(IN) :: d
    REAL(4), INTENT(IN) :: rho

    REAL(4) :: gp, w

    gp = abs((rho - rho_0)/rho_0) * gravity
    w = (4*gp**2 / (225*nu_mol))**(0.3333) * d

    terminal_velocity_allen = sign(w, rho_0 - rho)

  END FUNCTION terminal_velocity_allen

  REAL(4) PURE FUNCTION terminal_velocity_newton(d, rho)
    ! Newton's terminal velocity proportional to sqrt(d)
    REAL(4), INTENT(IN) :: d
    REAL(4), INTENT(IN) :: rho

    REAL(4) :: gp, w

    gp = abs((rho - rho_0)/rho_0) * gravity
    w = sqrt(4 * gp * d / (3*0.44))

    terminal_velocity_newton = sign(w, rho_0 - rho)

  END FUNCTION terminal_velocity_newton

  REAL(4) PURE FUNCTION terminal_velocity_rubey(d, rho)
    ! following Rubey 1933 formula
    REAL(4), INTENT(IN) :: d
    REAL(4), INTENT(IN) :: rho

    REAL(8) :: gp, dp
    REAL(4) :: w

    gp = abs((rho - rho_0)/rho_0) * gravity

    dp = gp * (d**3) / (nu_mol**2)
    w = sqrt(gp*d)*(sqrt(0.66 + 36/dp) - sqrt(36/dp))

    terminal_velocity_rubey = sign(w, rho_0 - rho)

  END FUNCTION terminal_velocity_rubey

  REAL(4) PURE FUNCTION terminal_velocity_kaiser(d, rho)
    ! following Keiser 2019 quadric regression for microplastic particle
    !  (gravity is fixed as 9.8 [m/s^2]), d is nominal diameter (ESD) for
    REAL(4), INTENT(IN) :: d
    REAL(4), INTENT(IN) :: rho

    REAL(4) :: drho, esd, w

    esd  = d*1.0E6           ![um]
    drho = abs(rho - rho_0)  ![kg/m^3]

    w = 11.68 + 0.1991*esd + 0.0004*(esd**2) - 0.0993*drho + 0.0002*(drho**2)
    w = w / 86400 !unit conversion [m/day] -> [m/s]

    terminal_velocity_kaiser = sign(w, rho_0 - rho)

  END FUNCTION terminal_velocity_kaiser

  REAL(4) PURE FUNCTION terminal_velocity_kaiser_mod(d, rho)
    ! following Keiser 2019 quadric regression for microplastic particle
    !  (gravity is fixed as 9.8 [m/s^2]), d is nominal diameter (ESD) for
    REAL(4), INTENT(IN) :: d
    REAL(4), INTENT(IN) :: rho

    REAL(4) :: drho, esd, w

    esd  = d*1.0E6           ![um]
    drho = abs(rho - rho_0)  ![kg/m^3]

    w = 11.68 + 0.01991*esd + 0.004*(esd**2) - 0.0993*drho + 0.0002*(drho**2)
     ! coefficients for 2nd and 3rd terms in the original paper might be typo and corrected

    w = w / 86400 !unit conversion [m/day] -> [m/s]

    terminal_velocity_kaiser_mod = sign(w, rho_0 - rho)

  END FUNCTION terminal_velocity_kaiser_mod

!------------------------------------------------------------------------------------------------------------

  LOGICAL PURE FUNCTION is_csvfile(filename)
    CHARACTER(*), INTENT(IN) :: filename

    is_csvfile = (len_trim(filename) >= 4 .AND. &
                 (index(trim(filename), '.csv', back=.TRUE.)==len_trim(filename)-3 .OR.  &
                  index(trim(filename), '.CSV', back=.TRUE.)==len_trim(filename)-3))
  END FUNCTION is_csvfile

!------------------------------------------------------------------------------------------------------------

  INTEGER(4) PURE FUNCTION xorshift(x)
    INTEGER(4), INTENT(IN) :: x

    xorshift = x

    IF (xorshift > 0) THEN
       xorshift = IEOR(ISHFT(xorshift,  13), xorshift)
       xorshift = IEOR(ISHFT(xorshift, -17), xorshift)
       xorshift = IEOR(ISHFT(xorshift,   5), xorshift)
    ELSE
       xorshift = IEOR(ISHFT(xorshift,  -3), xorshift)
       xorshift = IEOR(ISHFT(xorshift,   7), xorshift)
       xorshift = IEOR(ISHFT(xorshift, -29), xorshift)
    END IF
  END FUNCTION xorshift

  REAL(4) FUNCTION urand(seed)
    INTEGER, INTENT(INOUT) :: seed
    REAL(4), PARAMETER :: base = 0.5E0**31
    !return [-1:1] uniform distribution random number (variance=1/3)

    seed = xorshift(seed)
    urand = seed * base

  END FUNCTION urand

  REAL(4) FUNCTION nrand(seed)
    INTEGER, INTENT(INOUT) :: seed

    nrand = urand(seed)
    nrand = nrand + urand(seed)
    nrand = nrand + urand(seed)
  END FUNCTION nrand

END MODULE particles
