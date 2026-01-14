#include "macro.h"

MODULE tracers
  USE misc
  USE velocity
  USE subgrid
  USE geometry
  USE io

  IMPLICIT NONE
  SAVE
  PRIVATE

  PUBLIC tracer
  PUBLIC diff
  PUBLIC init_tracers, finalize_tracers
  PUBLIC add_tracer
  PUBLIC step_tracers
  PUBLIC checkout_tracers
  PUBLIC tracer_name,  n_tracer
  PUBLIC tracer_index, tracer_index_t, tracer_index_s
  PUBLIC flux_quickest, flux_upwind
  PUBLIC report_tracers
  PUBLIC update_tracer_boundary
  PUBLIC max_tracer
  PUBLIC tracer_info
  PUBLIC TRC_KIND

#ifdef TRACER_REAL4
  INTEGER, PARAMETER :: TRC_KIND = 4
#else
  INTEGER, PARAMETER :: TRC_KIND = 8
#endif

  INTEGER, PARAMETER :: max_tracer = 64

  INTEGER :: tracer_index_t = 0
  INTEGER :: tracer_index_s = 0

  INTEGER :: n_tracer

  REAL(TRC_KIND), ALLOCATABLE :: tracer(:,:,:,:)

  REAL(TRC_KIND), ALLOCATABLE :: diff(:,:,:,:)

  TYPE tracer_info_struct
     CHARACTER(32) :: name
     INTEGER       :: advscheme
     INTEGER       :: diffscheme
     LOGICAL       :: correctdiv
     LOGICAL       :: maskflux
     REAL(4)       :: diffh
     REAL(4)       :: diffv
     REAL(4)       :: csmag
     REAL(4)       :: cgls
     REAL(4)       :: cles
     REAL(4)       :: u
     REAL(4)       :: v
     REAL(4)       :: w
     REAL(TRC_KIND):: offsetx
     REAL(TRC_KIND):: offsety
     REAL(TRC_KIND):: offsetz
     LOGICAL       :: rmvsfc
     LOGICAL       :: rmvbtm
     LOGICAL       :: fixed
     CHARACTER(32) :: convtgt
     INTEGER       :: n_convtgt
     INTEGER(1), POINTER :: convflag(:,:,:)
     REAL(8),    POINTER :: convrate(:,:,:)
  END type tracer_info_struct

  TYPE(tracer_info_struct) :: tracer_info(max_tracer)
  TYPE(tracer_info_struct) :: default_info


CONTAINS
  SUBROUTINE init_tracers
    USE parameters, ONLY: tracer_advection_scheme, &
                          tracer_diffusion_scheme, &
                          tracer_correct_div,      &
                          tracer_mask_sfcflux,     &
                          default_diffh => diffh,  &
                          default_diffv => diffv,  &
                          c_smagorinsky_diff
    INTEGER :: n

    default_info%name = ''
    default_info%advscheme  = tracer_advection_scheme
    default_info%diffscheme = tracer_diffusion_scheme
    default_info%correctdiv = tracer_correct_div
    default_info%maskflux   = tracer_mask_sfcflux
    default_info%diffh      = default_diffh
    default_info%diffv      = default_diffv
    default_info%csmag      = c_smagorinsky_diff
    default_info%cgls       = 1.0
    default_info%cles       = 1.0/prandtl_les
    default_info%u          = 0.0
    default_info%v          = 0.0
    default_info%w          = 0.0
    default_info%offsetx    = 0.0
    default_info%offsety    = 0.0
    default_info%offsetz    = 0.0
    default_info%rmvsfc     = .FALSE.
    default_info%rmvbtm     = .FALSE.
    default_info%convtgt    = ""

    DO n=1, max_tracer
       tracer_info(n) = default_info
    END DO

    CALL read_namelist

    IF (n_tracer > 0) THEN
       ALLOCATE(tracer(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv, 1:n_tracer))
       tracer(:,:,:,:) = UNDEF
    END IF

    DO n=1, n_tracer
       CALL initialize(n)
    END DO

    tracer_index_t = tracer_index('T')
    tracer_index_s = tracer_index('S')

    ALLOCATE(diff(0:isize+1,0:jsize+1,0:ksize+1,3))

!$OMP PARALLEL WORKSHARE
    diff(:,:,:,:) = 0.0D0
!$OMP END PARALLEL WORKSHARE

  CONTAINS
    SUBROUTINE read_namelist
      CHARACTER(32) :: varname
      INTEGER       :: advection_scheme
      INTEGER       :: diffusion_scheme
      REAL(4)       :: diffh, diff_h
      REAL(4)       :: diffv, diff_v
      REAL(4)       :: c_smagorinsky
      REAL(4)       :: c_gls
      REAL(4)       :: c_les
      LOGICAL       :: correct_div
      LOGICAL       :: mask_sfcflux
      REAL(4)       :: u
      REAL(4)       :: v
      REAL(4)       :: w
      REAL(4)       :: addin_u
      REAL(4)       :: addin_v
      REAL(4)       :: addin_w
      LOGICAL       :: remove_surface
      LOGICAL       :: remove_bottom
      LOGICAL       :: fixed
      INTEGER       :: report_interval
      REAL(TRC_KIND):: cycleoffset_x
      REAL(TRC_KIND):: cycleoffset_y
      REAL(TRC_KIND):: cycleoffset_z
      REAL(TRC_KIND):: x_offset
      REAL(TRC_KIND):: y_offset
      REAL(TRC_KIND):: z_offset
      CHARACTER(32) :: convert_target

      INTEGER :: n
      INTEGER :: iostat
      CHARACTER(256) :: iomsg

      NAMELIST / tracer /    &
           varname,          &
           advection_scheme, &
           diffusion_scheme, &
           correct_div,      &
           diffh, diff_h,    &
           diffv, diff_v,    &
           c_smagorinsky,    &
           c_gls,            &
           c_les,            &
           u,                &
           v,                &
           w,                &
           addin_u,          &
           addin_v,          &
           addin_w,          &
           cycleoffset_x,    &
           cycleoffset_y,    &
           cycleoffset_z,    &
           x_offset,         &
           y_offset,         &
           z_offset,         &
           remove_surface,   &
           remove_bottom,    &
           convert_target,   &
           fixed,            &
           report_interval

      IF (rank==0) THEN
         n_tracer = 0

         REWIND(CONFIG_UNIT)
         DO
            varname          = ''
            advection_scheme = default_info%advscheme
            diffusion_scheme = default_info%diffscheme
            correct_div      = default_info%correctdiv
            mask_sfcflux     = default_info%maskflux
            diffh            = default_info%diffh
            diffv            = default_info%diffv
            diff_h           = UNDEF
            diff_v           = UNDEF
            c_smagorinsky    = default_info%csmag
            c_gls            = default_info%cgls
            c_les            = default_info%cles
            u                = default_info%u
            v                = default_info%v
            w                = default_info%w
            addin_u          = UNDEF
            addin_v          = UNDEF
            addin_w          = UNDEF
            remove_surface   = default_info%rmvsfc
            remove_bottom    = default_info%rmvbtm
            cycleoffset_x    = default_info%offsetx
            cycleoffset_y    = default_info%offsety
            cycleoffset_z    = default_info%offsetz
            x_offset         = UNDEF
            y_offset         = UNDEF
            z_offset         = UNDEF
            convert_target   = ""
            fixed            = .FALSE.

            READ(CONFIG_UNIT, NML=tracer, IOSTAT=iostat, IOMSG=iomsg)

            IF (iostat < 0) EXIT

            CALL assert(iostat == 0, "failed to read TRACER namelist for varname = '"//trim(varname)// "'", iomsg)
            CALL assert(trim(varname)/='',   "varname is not specified in TRACER namelist.")
            CALL assert(valid_name(varname), "invalid varname '"//trim(varname)//"' in TRACER namelist.")

            n = tracer_index(trim(varname))

            IF (n==0) THEN
               n_tracer = n_tracer + 1
               n = n_tracer
               tracer_info(n)%name = trim(varname)
            END IF
            !for backward compatibility
            IF (addin_u  /= UNDEF) u = addin_u
            IF (addin_v  /= UNDEF) v = addin_v
            IF (addin_w  /= UNDEF) w = addin_w
            IF (diff_h   /= UNDEF) diffh = diff_h
            IF (diff_v   /= UNDEF) diffv = diff_v
            IF (x_offset /= UNDEF) cycleoffset_x = x_offset
            IF (y_offset /= UNDEF) cycleoffset_y = y_offset
            IF (z_offset /= UNDEF) cycleoffset_z = z_offset
            tracer_info(n)%advscheme  = advection_scheme
            tracer_info(n)%diffscheme = diffusion_scheme
            tracer_info(n)%correctdiv = correct_div .AND. (advection_scheme /= 0)
            tracer_info(n)%maskflux   = mask_sfcflux
            tracer_info(n)%diffh      = diffh
            tracer_info(n)%diffv      = diffv
            tracer_info(n)%csmag      = c_smagorinsky
            tracer_info(n)%cgls       = c_gls
            tracer_info(n)%cles       = c_les
            tracer_info(n)%u          = u
            tracer_info(n)%v          = v
            tracer_info(n)%w          = w
            tracer_info(n)%rmvsfc     = remove_surface
            tracer_info(n)%rmvbtm     = remove_bottom
            tracer_info(n)%offsetx    = cycleoffset_x
            tracer_info(n)%offsety    = cycleoffset_y
            tracer_info(n)%offsetz    = cycleoffset_z
            tracer_info(n)%convtgt    = convert_target
            tracer_info(n)%fixed      = fixed
         END DO

         CALL assert(n_tracer <= max_tracer, "n_tracer should be less than MAX_TRACER!")

#ifdef DEBUG
         WRITE(REPORT_UNIT, *) '--tracer varname  advection_scheme  diffusion_scheme  correct_div --'
         DO n=1, n_tracer
            WRITE(REPORT_UNIT, '(2X,A10)',    ADVANCE='NO') "'"//trim(tracer_info(n)%name)//"'"
            WRITE(REPORT_UNIT, '(2X,I6,4X)',  ADVANCE='NO') tracer_info(n)%advscheme
            WRITE(REPORT_UNIT, '(2X,I6,4X)',  ADVANCE='NO') tracer_info(n)%diffscheme
            WRITE(REPORT_UNIT, '(2X,L6,4X)',  ADVANCE='NO') tracer_info(n)%correctdiv
            WRITE(REPORT_UNIT, *)
         END DO
         WRITE(REPORT_UNIT,*)
#endif
      END IF

      CALL bcast(n_tracer)
      DO n=1, n_tracer
         CALL bcast(tracer_info(n)%name)
         CALL bcast(tracer_info(n)%advscheme)
         CALL bcast(tracer_info(n)%diffscheme)
         CALL bcast(tracer_info(n)%correctdiv)
         CALL bcast(tracer_info(n)%diffh)
         CALL bcast(tracer_info(n)%diffv)
         CALL bcast(tracer_info(n)%csmag)
         CALL bcast(tracer_info(n)%cgls)
         CALL bcast(tracer_info(n)%u)
         CALL bcast(tracer_info(n)%v)
         CALL bcast(tracer_info(n)%w)
         CALL bcast(tracer_info(n)%offsetx)
         CALL bcast(tracer_info(n)%offsety)
         CALL bcast(tracer_info(n)%offsetz)
         CALL bcast(tracer_info(n)%rmvsfc)
         CALL bcast(tracer_info(n)%rmvbtm)
         CALL bcast(tracer_info(n)%convtgt)
         CALL bcast(tracer_info(n)%fixed)
      END DO

    END SUBROUTINE read_namelist

  END SUBROUTINE init_tracers

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE add_tracer(name, index)
    CHARACTER(*), INTENT(IN) :: name
    INTEGER,      INTENT(OUT), OPTIONAL :: index

    REAL(TRC_KIND), ALLOCATABLE :: tmp(:,:,:,:)

    INTEGER :: i, j, k

    IF (tracer_index(name) /= 0) THEN
       IF (present(index)) index = tracer_index(name)
       RETURN
    END IF

    IF (n_tracer==0) THEN
       ALLOCATE(tracer(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv, 1:1))
       tracer(:,:,:,1) = UNDEF
    ELSE
       ALLOCATE(tmp(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv, 1:n_tracer))
!$OMP PARALLEL WORKSHARE
       tmp(:,:,:,:) = tracer(:,:,:,:)
!$OMP END PARALLEL WORKSHARE

       DEALLOCATE(tracer)
       ALLOCATE(tracer(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv, 1:n_tracer+1))

!$OMP PARALLEL WORKSHARE
       tracer(:,:,:,1:n_tracer) = tmp(:,:,:,:)
       tracer(:,:,:,n_tracer+1) = UNDEF
!$OMP END PARALLEL WORKSHARE
       DEALLOCATE(tmp)
    END IF

    n_tracer = n_tracer + 1

    tracer_info(n_tracer) = default_info
    tracer_info(n_tracer)%name = trim(name)

    CALL initialize(n_tracer)

    IF (present(index)) index = n_tracer
  END SUBROUTINE add_tracer

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE initialize(index)
    INTEGER, INTENT(IN) :: index

    REAL(TRC_KIND) :: tmp_xz(isize,ksize)
    REAL(TRC_KIND) :: tmp_yz(jsize,ksize)
    CHARACTER(32) :: name
    INTEGER :: i, j, k
    LOGICAL :: stat

    name = tracer_name(index)

    CALL initial_data(name, tracer(:,:,:,index), default=.TRUE.)

    IF (open_e) THEN
       tmp_yz(:,:) = UNDEF
       IF (radi_tracers .AND. radi_e) THEN
          CALL initial_data(trim(name)//'_EAST', tmp_yz, section='YZ', default=perfect_restart)
       ELSE
          CALL checkin(trim(name)//'_EAST', tmp_yz, section='YZ', time=t_start)
       END IF
       IF (icoord==ipes-1) THEN
!$OMP PARALLEL DO
          DO k=1, ksize
          DO j=1, jsize
             IF (tmp_yz(j,k)==UNDEF) THEN
!!!!!!!!!!!YJKIM               
!                tracer(isize,j,k,index) = tracer(isize-1,j,k,index)
                tracer(isize+1,j,k,index) = tracer(isize,j,k,index)
             ELSE
                tracer(isize+1,j,k,index) = tmp_yz(j,k)
             END IF
          END DO
          END DO
       END IF
    END IF
    IF (open_w) THEN
       tmp_yz(:,:) = UNDEF
       IF (radi_tracers .AND. radi_w) THEN
          CALL initial_data(trim(name)//'_WEST', tmp_yz, section='YZ', default=perfect_restart)
       ELSE
          CALL checkin(trim(name)//'_WEST', tmp_yz, section='YZ', time=t_start)
       END IF
       IF (icoord==0) THEN
!$OMP PARALLEL DO
          DO k=1, ksize
          DO j=1, jsize
             IF (tmp_yz(j,k)==UNDEF) THEN
!!!!!!!!!!!YJKIM
!                tracer(1,j,k,index) = tracer(2,j,k,index)
                tracer(0,j,k,index) = tracer(1,j,k,index)
             ELSE
                tracer(0,j,k,index) = tmp_yz(j,k)
             END IF
          END DO
          END DO
       END IF
    END IF
    IF (open_n) THEN
       tmp_xz(:,:) = UNDEF
       IF (radi_tracers .AND. radi_n) THEN
          CALL initial_data(trim(name)//'_NORTH', tmp_xz, section='XZ', default=perfect_restart)
       ELSE
          CALL checkin(trim(name)//'_NORTH', tmp_xz, section='XZ', time=t_start)
       END IF
       IF (jcoord==jpes-1) THEN
!$OMP PARALLEL DO
          DO k=1, ksize
          DO i=1, isize
             IF (tmp_xz(i,k)==UNDEF) THEN
!!!!!!!!!!!YJKIM               
!                tracer(i,jsize,k,index) = tracer(i,jsize-1,k,index)
                tracer(i,jsize+1,k,index) = tracer(i,jsize,k,index)
             ELSE
                tracer(i,jsize+1,k,index) = tmp_xz(i,k)
             END IF
          END DO
          END DO
       END IF
    END IF
    IF (open_s) THEN
       tmp_xz(:,:) = UNDEF
       IF (radi_tracers .AND. radi_s) THEN
          CALL initial_data(trim(name)//'_SOUTH', tmp_xz, section='XZ', default=perfect_restart)
       ELSE
          CALL checkin(trim(name)//'_SOUTH', tmp_xz, section='XZ', time=t_start)
       END IF
       IF (jcoord==0) THEN
!$OMP PARALLEL DO
          DO k=1, ksize
          DO i=1, isize
             IF (tmp_xz(i,k)==UNDEF) THEN
!!!!!!!!!!!YJKIM
!                tracer(i,1,k,index) = tracer(i,2,k,index)
                tracer(i,0,k,index) = tracer(i,1,k,index)
             ELSE
                tracer(i,0,k,index) = tmp_xz(i,k)
             END IF
          END DO
          END DO
       END IF
    END IF

    CALL update_tracer_boundary(index)

  END SUBROUTINE initialize

!-----------------------------------------------------------------------------------------------------------------------


  SUBROUTINE update_tracer_boundary(index)
    INTEGER, INTENT(IN) :: index
    REAL(TRC_KIND) :: tmp_xz(1-slv:isize+slv,1-slv:ksize+slv)
    INTEGER :: k

    CALL update_boundary(tracer(:,:,:,index))

    IF (open_n .AND. jcoord==jpes-1) CALL update_boundary_xz(tracer(:,jsize+1,:,index))
    IF (open_s .AND. jcoord==0)      CALL update_boundary_xz(tracer(:,0,      :,index))
!!!!!!!!!!!!!!!YJKIM
    IF (open_e .AND. icoord==ipes-1) THEN
!$OMP PARALLEL WORKSHARE
      !  tracer(isize-4,:,:,index) = tracer(isize-5,:,:,index)
      !  tracer(isize-3,:,:,index) = tracer(isize-4,:,:,index)
      !  tracer(isize-2,:,:,index) = tracer(isize-3,:,:,index)
      !  tracer(isize-1,:,:,index) = tracer(isize-2,:,:,index)
      !  tracer(isize,:,:,index) = tracer(isize-1,:,:,index)
      !  tracer(isize+1,:,:,index) = tracer(isize,:,:,index)
       tracer(isize+2,:,:,index) = tracer(isize+1,:,:,index)
!$OMP END PARALLEL WORKSHARE
    END IF
    IF (open_w .AND. icoord==0) THEN
!$OMP PARALLEL WORKSHARE
      !  tracer(5,:,:,index) = tracer(6,:,:,index)
      !  tracer(4,:,:,index) = tracer(5,:,:,index)
      !  tracer(3,:,:,index) = tracer(4,:,:,index)
      !  tracer(2,:,:,index) = tracer(3,:,:,index)
      !  tracer(1,:,:,index) = tracer(2,:,:,index)
      !  tracer(0,:,:,index) = tracer(1,:,:,index)
       tracer(-1,:,:,index) = tracer(0,:,:,index)
!$OMP END PARALLEL WORKSHARE
    END IF
    IF (open_n .AND. jcoord==jpes-1) THEN
!$OMP PARALLEL WORKSHARE
      !  tracer(:,jsize-4,:,index) = tracer(:,jsize-5,:,index)
      !  tracer(:,jsize-3,:,index) = tracer(:,jsize-4,:,index)
      !  tracer(:,jsize-2,:,index) = tracer(:,jsize-3,:,index)
      !  tracer(:,jsize-1,:,index) = tracer(:,jsize-2,:,index)
      !  tracer(:,jsize,:,index) = tracer(:,jsize-1,:,index)
      !  tracer(:,jsize+1,:,index) = tracer(:,jsize,:,index)
       tracer(:,jsize+2,:,index) = tracer(:,jsize+1,:,index)
!$OMP END PARALLEL WORKSHARE
    END IF
    IF (open_s .AND. jcoord==0) THEN
!$OMP PARALLEL WORKSHARE
      !  tracer(:,5,:,index) = tracer(:,6,:,index)
      !  tracer(:,4,:,index) = tracer(:,5,:,index)
      !  tracer(:,3,:,index) = tracer(:,4,:,index)
      !  tracer(:,2,:,index) = tracer(:,3,:,index)
      !  tracer(:,1,:,index) = tracer(:,2,:,index)
      !  tracer(:,0,:,index) = tracer(:,1,:,index)
       tracer(:,-1,:,index) = tracer(:,0,:,index)
!$OMP END PARALLEL WORKSHARE
    END IF

    IF (open_e .AND. open_n .AND. icoord==ipes-1 .AND. jcoord==jpes-1) THEN
!$OMP PARALLEL DO
       DO k=1-slv, ksize+slv
          tracer(isize+1:isize+2,jsize+1:jsize+2,k,index) = 0.5*(tracer(isize+1,jsize,k,index)+tracer(isize,jsize+1,k,index))
       END DO
    END IF
    IF (open_e .AND. open_s .AND. icoord==ipes-1 .AND. jcoord==0) THEN
!$OMP PARALLEL DO
       DO k=1-slv, ksize+slv
          tracer(isize+1:isize+2,-1:0,k,index) = 0.5*(tracer(isize+1,1,k,index)+tracer(isize,0,k,index))
       END DO
    END IF
    IF (open_w .AND. open_n .AND. icoord==0 .AND. jcoord==jpes-1) THEN
!$OMP PARALLEL DO
       DO k=1-slv, ksize+slv
          tracer(-1:0,jsize+1:jsize+2,k,index) = 0.5*(tracer(0,jsize,k,index)+tracer(1,jsize+1,k,index))
       END DO
    END IF
    IF (open_w .AND. open_s .AND. icoord==0 .AND. jcoord==0) THEN
!$OMP PARALLEL DO
       DO k=1-slv, ksize+slv
          tracer(-1:0,-1:0,k,index) = 0.5*(tracer(0,1,k,index)+tracer(1,0,k,index))
       END DO
    END IF

    CALL cycleoffset(index)

  END SUBROUTINE update_tracer_boundary

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE cycleoffset(index)
    INTEGER, INTENT(IN) :: index

!$OMP PARALLEL
    IF (cycle_x .AND. tracer_info(index)%offsetx /= 0.0) THEN
       IF (icoord==0) THEN
!$OMP WORKSHARE
          tracer(1-slv:0,:,:,index) = tracer(1-slv:0,:,:,index) - tracer_info(index)%offsetx
!$OMP END WORKSHARE
       END IF
       IF (icoord==ipes-1) THEN
!$OMP WORKSHARE
          tracer(isize+1:isize+slv,:,:,index) = tracer(isize+1:isize+slv,:,:,index) + tracer_info(index)%offsetx
!$OMP END WORKSHARE
       END IF
    END IF

    IF (cycle_y .AND. tracer_info(index)%offsety /= 0.0) THEN
       IF (jcoord==0) THEN
!$OMP WORKSHARE
          tracer(:,1-slv:0,:,index) = tracer(:,1-slv:0,:,index) - tracer_info(index)%offsety
!$OMP END WORKSHARE
       END IF
       IF (jcoord==jpes-1) THEN
!$OMP WORKSHARE
          tracer(:,jsize+1:jsize+slv,:,index) = tracer(:,jsize+1:jsize+slv,:,index) + tracer_info(index)%offsety
!$OMP END WORKSHARE
       END IF
    END IF

    IF (cycle_z .AND. tracer_info(index)%offsetz /= 0.0) THEN
       IF (kcoord==0) THEN
!$OMP WORKSHARE
          tracer(:,:,1-slv:0,index) = tracer(:,:,1-slv:0,index) - tracer_info(index)%offsetz
!$OMP END WORKSHARE
       END IF
       IF (kcoord==kpes-1) THEN
!$OMP WORKSHARE
          tracer(:,:,ksize+1:ksize+slv,index) = tracer(:,:,ksize+1:ksize+slv,index) + tracer_info(index)%offsetz
!$OMP END WORKSHARE
       END IF
    END IF
!$OMP END PARALLEL

  END SUBROUTINE cycleoffset

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE finalize_tracers
    IF (ALLOCATED(tracer)) DEALLOCATE(tracer)
  END SUBROUTINE finalize_tracers

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE advection_upwind(a, u, v, w, fx, fy, fz)
    REAL(TRC_KIND), INTENT(IN)  :: a(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)

    REAL(TRC_KIND), INTENT(IN)  :: u( -slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
    REAL(TRC_KIND), INTENT(IN)  :: v(1-slv:isize+slv,  -slv:jsize+slv, 1-slv:ksize+slv)
    REAL(TRC_KIND), INTENT(IN)  :: w(1-slv:isize+slv, 1-slv:jsize+slv,  -slv:ksize+slv)

    REAL(TRC_KIND), INTENT(OUT) :: fx(0:isize, 1:jsize, 1:ksize)
    REAL(TRC_KIND), INTENT(OUT) :: fy(1:isize, 0:jsize, 1:ksize)
    REAL(TRC_KIND), INTENT(OUT) :: fz(1:isize, 1:jsize, 0:ksize)

    INTEGER :: i,j,k

!$OMP PARALLEL DO
    DO k=0, ksize
       DO j=1, jsize
       DO i=1, isize
          fz(i,j,k) = flux_upwind(w(i,j,k), a(i,j,k), a(i,j,k+1)) * imask3d(i,j,k)*imask3d(i,j,k+1)
       END DO
       END DO

       IF (k < 1) CYCLE

       DO j=1, jsize
       DO i=0, isize
          fx(i,j,k) = flux_upwind(u(i,j,k), a(i,j,k), a(i+1,j,k)) * imask3d(i,j,k)*imask3d(i+1,j,k)
       END DO
       END DO

       DO j=0, jsize
       DO i=1, isize
          fy(i,j,k) = flux_upwind(v(i,j,k), a(i,j,k), a(i,j+1,k)) * imask3d(i,j,k)*imask3d(i,j+1,k)
       END DO
       END DO
    END DO

  END SUBROUTINE advection_upwind

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE advection_quickest(a, u, v, w, fx, fy, fz)
    REAL(TRC_KIND), INTENT(IN)  :: a(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)

    REAL(TRC_KIND), INTENT(IN)  :: u( -slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
    REAL(TRC_KIND), INTENT(IN)  :: v(1-slv:isize+slv,  -slv:jsize+slv, 1-slv:ksize+slv)
    REAL(TRC_KIND), INTENT(IN)  :: w(1-slv:isize+slv, 1-slv:jsize+slv,  -slv:ksize+slv)

    REAL(TRC_KIND), INTENT(OUT) :: fx(0:isize, 1:jsize, 1:ksize)
    REAL(TRC_KIND), INTENT(OUT) :: fy(1:isize, 0:jsize, 1:ksize)
    REAL(TRC_KIND), INTENT(OUT) :: fz(1:isize, 1:jsize, 0:ksize)

    INTEGER :: i,j,k

!$OMP PARALLEL DO
    DO k=0, ksize
       DO j=1, jsize
       DO i=1, isize
          fz(i,j,k) = flux_quickest(w(i,j,k), a(i,j,k-1),       a(i,j,k),       a(i,j,k+1),       a(i,j,k+2), &
                                             dz(    k-1),      dz(    k),      dz(    k+1),      dz(    k+2), &
                                        imask3d(i,j,k-1), imask3d(i,j,k), imask3d(i,j,k+1), imask3d(i,j,k+2), dtime)
       END DO
       END DO

       IF (k < 1) CYCLE

       DO j=1, jsize
       DO i=0, isize
          fx(i,j,k) = flux_quickest(u(i,j,k), a(i-1,j,k),       a(i,j,k),       a(i+1,j,k),       a(i+2,j,k), &
                                             dx(i-1,j),        dx(i,j),        dx(i+1,j),        dx(i+2,j),   &
                                        imask3d(i-1,j,k), imask3d(i,j,k), imask3d(i+1,j,k), imask3d(i+2,j,k), dtime)
       END DO
       END DO

       DO j=0, jsize
       DO i=1, isize
          fy(i,j,k) = flux_quickest(v(i,j,k), a(i,j-1,k),       a(i,j,k),       a(i,j+1,k),       a(i,j+2,k), &
                                             dy(i,j-1),        dy(i,j),        dy(i,j+1),        dy(i,j+2),   &
                                        imask3d(i,j-1,k), imask3d(i,j,k), imask3d(i,j+1,k), imask3d(i,j+2,k), dtime)
       END DO
       END DO
    END DO

  END SUBROUTINE advection_quickest

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE advection_cosmic_nottuned(a, u, v, w, fx, fy, fz)
    REAL(TRC_KIND), INTENT(IN)  :: a(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)

    REAL(TRC_KIND), INTENT(IN)  :: u( -slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
    REAL(TRC_KIND), INTENT(IN)  :: v(1-slv:isize+slv,  -slv:jsize+slv, 1-slv:ksize+slv)
    REAL(TRC_KIND), INTENT(IN)  :: w(1-slv:isize+slv, 1-slv:jsize+slv,  -slv:ksize+slv)

    REAL(TRC_KIND), INTENT(OUT) :: fx(0:isize, 1:jsize, 1:ksize)
    REAL(TRC_KIND), INTENT(OUT) :: fy(1:isize, 0:jsize, 1:ksize)
    REAL(TRC_KIND), INTENT(OUT) :: fz(1:isize, 1:jsize, 0:ksize)

    REAL(TRC_KIND) :: a_x( 0:isize+1, -1:jsize+2, -1:ksize+2)
    REAL(TRC_KIND) :: a_y(-1:isize+2,  0:jsize+1, -1:ksize+2)
    REAL(TRC_KIND) :: a_z(-1:isize+2, -1:jsize+2,  1:ksize)

    REAL(TRC_KIND) :: a_yz(-1:isize+2, 1:jsize,   1:ksize)
    REAL(TRC_KIND) :: a_zx( 1:isize,  -1:jsize+2, 1:ksize)
    REAL(TRC_KIND) :: a_xy( 1:isize,   1:jsize,  -1:ksize+2)

    REAL(TRC_KIND) :: udt( 0:isize+1, -1:jsize+2, -1:ksize+2)
    REAL(TRC_KIND) :: vdt(-1:isize+2,  0:jsize+1, -1:ksize+2)
    REAL(TRC_KIND) :: wdt(-1:isize+2, -1:jsize+2,  1:ksize)

    INTEGER :: i,j,k

!$OMP PARALLEL
!$OMP DO
    DO k=-1, ksize+2
       DO j=-1, jsize+2
       DO i= 0, isize+1
          udt(i,j,k) = imask3d(i-1,j,k)*imask3d(i,j,k)*imask3d(i+1,j,k) * 0.5 * (u(i-1,j,k)+u(i,j,k)) * dtime * idx2(i,j)
          a_x(i,j,k) = a(i,j,k) - udt(i,j,k) * (a(i+1,j,k)-a(i-1,j,k))
       END DO
       END DO

       DO j= 0, jsize+1
       DO i=-1, isize+2
          vdt(i,j,k) = imask3d(i,j-1,k)*imask3d(i,j,k)*imask3d(i,j+1,k) * 0.5 * (v(i,j-1,k)+v(i,j,k)) * dtime * idy2(i,j)
          a_y(i,j,k) = a(i,j,k) - vdt(i,j,k) * (a(i,j+1,k)-a(i,j-1,k))
       END DO
       END DO

       IF (k < 1 .OR. k > ksize) CYCLE

       DO j=-1, jsize+2
       DO i=-1, isize+2
          wdt(i,j,k) = imask3d(i,j,k-1)*imask3d(i,j,k)*imask3d(i,j,k+1) * 0.5 * (w(i,j,k-1)+w(i,j,k)) * dtime * idz2(k)
          a_z(i,j,k) = a(i,j,k) - wdt(i,j,k) * (a(i,j,k+1)-a(i,j,k-1))
       END DO
       END DO
    END DO

!$OMP DO
    DO k=-1, ksize+2
       DO j= 1, jsize
       DO i= 1, isize
          a_xy(i,j,k) = (a(i,j,k) + a_x(i,j,k) + a_y(i,j,k) - 0.5*vdt(i,j,k)*(a_x(i,j+1,k)-a_x(i,j-1,k)) &
                                                            - 0.5*udt(i,j,k)*(a_y(i+1,j,k)-a_y(i-1,j,k))) / 3.0
       END DO
       END DO

       IF (k < 1 .OR. k > ksize) CYCLE

       DO j= 1, jsize
       DO i=-1, isize+2
          a_yz(i,j,k) = (a(i,j,k) + a_y(i,j,k) + a_z(i,j,k) - 0.5*wdt(i,j,k)*(a_y(i,j,k+1)-a_y(i,j,k-1)) &
                                                            - 0.5*vdt(i,j,k)*(a_z(i,j+1,k)-a_z(i,j-1,k))) / 3.0
       END DO
       END DO

       DO j=-1, jsize+2
       DO i= 1, isize
          a_zx(i,j,k) = (a(i,j,k) + a_z(i,j,k) + a_x(i,j,k) - 0.5*udt(i,j,k)*(a_z(i+1,j,k)-a_z(i-1,j,k)) &
                                                            - 0.5*wdt(i,j,k)*(a_x(i,j,k+1)-a_x(i,j,k-1))) / 3.0
       END DO
       END DO
    END DO

!$OMP DO
    DO k=0, ksize
       DO j=1, jsize
       DO i=1, isize
          fz(i,j,k) = flux_quickest(w(i,j,k), a_xy(i,j,k-1),    a_xy(i,j,k),    a_xy(i,j,k+1),    a_xy(i,j,k+2), &
                                                dz(    k-1),      dz(    k),      dz(    k+1),      dz(    k+2), &
                                           imask3d(i,j,k-1), imask3d(i,j,k), imask3d(i,j,k+1), imask3d(i,j,k+2), dtime)
       END DO
       END DO

       IF (k < 1) CYCLE

       DO j=1, jsize
       DO i=0, isize
          fx(i,j,k) = flux_quickest(u(i,j,k), a_yz(i-1,j,k),    a_yz(i,j,k),    a_yz(i+1,j,k),    a_yz(i+2,j,k), &
                                                dx(i-1,j),        dx(i,j),        dx(i+1,j),        dx(i+2,j),   &
                                           imask3d(i-1,j,k), imask3d(i,j,k), imask3d(i+1,j,k), imask3d(i+2,j,k), dtime)

       END DO
       END DO

       DO j=0, jsize
       DO i=1, isize
          fy(i,j,k) = flux_quickest(v(i,j,k), a_zx(i,j-1,k),    a_zx(i,j,k),    a_zx(i,j+1,k),    a_zx(i,j+2,k), &
                                                dy(i,j-1),        dy(i,j),        dy(i,j+1),        dy(i,j+2),   &
                                           imask3d(i,j-1,k), imask3d(i,j,k), imask3d(i,j+1,k), imask3d(i,j+2,k), dtime)
       END DO
       END DO
    END DO
!$OMP END PARALLEL

  END SUBROUTINE advection_cosmic_nottuned

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE advection_cosmic(a, u, v, w, fx, fy, fz)
    REAL(TRC_KIND), INTENT(IN)  :: a(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)

    REAL(TRC_KIND), INTENT(IN)  :: u( -slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
    REAL(TRC_KIND), INTENT(IN)  :: v(1-slv:isize+slv,  -slv:jsize+slv, 1-slv:ksize+slv)
    REAL(TRC_KIND), INTENT(IN)  :: w(1-slv:isize+slv, 1-slv:jsize+slv,  -slv:ksize+slv)

    REAL(TRC_KIND), INTENT(OUT) :: fx(0:isize, 1:jsize, 1:ksize)
    REAL(TRC_KIND), INTENT(OUT) :: fy(1:isize, 0:jsize, 1:ksize)
    REAL(TRC_KIND), INTENT(OUT) :: fz(1:isize, 1:jsize, 0:ksize)

    REAL(TRC_KIND) :: a_x(0:3,  0:isize+1, -1:jsize+2)
    REAL(TRC_KIND) :: a_y(0:3, -1:isize+2,  0:jsize+1)
    REAL(TRC_KIND) :: a_z(     -1:isize+2, -1:jsize+2)

    REAL(TRC_KIND) :: a_yz(    -1:isize+2, 1:jsize)
    REAL(TRC_KIND) :: a_zx(     1:isize,  -1:jsize+2)
    REAL(TRC_KIND) :: a_xy(0:3, 1:isize,   1:jsize)

    REAL(TRC_KIND) :: udt( 0:isize+1, -1:jsize+2, 0:3)
    REAL(TRC_KIND) :: vdt(-1:isize+2,  0:jsize+1, 0:3)
    REAL(TRC_KIND) :: wdt(-1:isize+2, -1:jsize+2)

    INTEGER :: i, j, k
    INTEGER :: k0, k1, k2, k3
    INTEGER :: km0, km1, km2, km3
    INTEGER :: tmp
    INTEGER :: kstr, kend
    INTEGER :: tid

    tid  = 0
    kstr = 1
    kend = ksize

!$OMP PARALLEL PRIVATE(i, j, k, tid) &
!$OMP          PRIVATE(kstr, kend) &
!$OMP          PRIVATE(a_x, a_y, a_z, a_xy, a_yz, a_zx) &
!$OMP          PRIVATE(udt, vdt, wdt) &
!$OMP          PRIVATE(k0, k1, k2, k3, tmp) &
!$OMP          PRIVATE(km0, km1, km2, km3) &
!$OMP          shared(a, u, v, w)

!$  tid = omp_get_thread_num()
!$  kstr = kstr_t(tid)
!$  kend = kend_t(tid)

    k0=0
    k1=1
    k2=2
    k3=3

    DO k=kstr-4, kend
       tmp = k0
       k0 = k1
       k1 = k2
       k2 = k3
       k3 = tmp

       km3 = maskindex(k+2)

       DO j=-1, jsize+2
       DO i= 0, isize+1
          udt(i,j,k3) = imask3d(i-1,j,km3)*imask3d(i,j,km3)*imask3d(i+1,j,km3) * (u(i-1,j,k+2)+u(i,j,k+2))*0.5 * dtime * idx2(i,j)
          a_x(k3,i,j) = a(i,j,k+2) - udt(i,j,k3) * (a(i+1,j,k+2)-a(i-1,j,k+2))
       END DO
       END DO

       DO j= 0, jsize+1
       DO i=-1, isize+2
          vdt(i,j,k3) = imask3d(i,j-1,km3)*imask3d(i,j,km3)*imask3d(i,j+1,km3) * (v(i,j-1,k+2)+v(i,j,k+2))*0.5 * dtime * idy2(i,j)
          a_y(k3,i,j) = a(i,j,k+2) - vdt(i,j,k3) * (a(i,j+1,k+2)-a(i,j-1,k+2))
       END DO
       END DO

       DO j=1, jsize
       DO i=1, isize
          a_xy(k3,i,j) = (a(i,j,k+2) + a_x(k3,i,j) + a_y(k3,i,j) - 0.5*vdt(i,j,k3)*(a_x(k3,i,j+1)-a_x(k3,i,j-1)) &
                                                                 - 0.5*udt(i,j,k3)*(a_y(k3,i+1,j)-a_y(k3,i-1,j))) / 3.0
       END DO
       END DO

       IF (k <= kstr-2) CYCLE
!$   IF (tid/=0 .AND. k==kstr-1) CYCLE

       km0 = maskindex(k-1)
       km1 = maskindex(k)
       km2 = maskindex(k+1)

       DO j=1, jsize
       DO i=1, isize
          fz(i,j,k) = flux_quickest(w(i,j,k), a_xy(k0,i,j),     a_xy(k1,i,j),     a_xy(k2,i,j),     a_xy(k3,i,j), &
                                                   dz(k-1),            dz(k),          dz(k+1),          dz(k+2), &
                                          imask3d(i,j,km0), imask3d(i,j,km1), imask3d(i,j,km2), imask3d(i,j,km3), dtime)
       END DO
       END DO

       IF (k==kstr-1) CYCLE

       DO j=-1, jsize+2
       DO i=-1, isize+2
          wdt(i,j) = imask3d(i,j,km0)*imask3d(i,j,km1)*imask3d(i,j,km2) * (w(i,j,k-1)+w(i,j,k))*0.5 * dtime * idz2(k)
          a_z(i,j) = a(i,j,k) - wdt(i,j) * (a(i,j,k+1)-a(i,j,k-1))
       END DO
       END DO

       DO j= 1, jsize
       DO i=-1, isize+2
          a_yz(i,j) = (a(i,j,k) + a_y(k1,i,j) + a_z(i,j) - 0.5*wdt(i,j)   *(a_y(k2,i,j)-a_y(k0,i,j)) &
                                                         - 0.5*vdt(i,j,k1)*(a_z(i,j+1) -a_z(i,j-1))) / 3.0
       END DO
       END DO

       DO j=1, jsize
       DO i=0, isize
          fx(i,j,k) = flux_quickest(u(i,j,k), a_yz(i-1,j),        a_yz(i,j),        a_yz(i+1,j),        a_yz(i+2,j), &
                                                dx(i-1,j),          dx(i,j),          dx(i+1,j),          dx(i+2,j), &
                                           imask3d(i-1,j,km1), imask3d(i,j,km1), imask3d(i+1,j,km1), imask3d(i+2,j,km1), dtime)
       END DO
       END DO

       DO j=-1, jsize+2
       DO i= 1, isize
             a_zx(i,j) = (a(i,j,k) + a_z(i,j) + a_x(k1,i,j) - 0.5*udt(i,j,k1)*(a_z(i+1,j) -a_z(i-1,j)) &
                                                            - 0.5*wdt(i,j)   *(a_x(k2,i,j)-a_x(k0,i,j))) / 3.0
       END DO
       END DO

       DO j=0, jsize
       DO i=1, isize
          fy(i,j,k) = flux_quickest(v(i,j,k), a_zx(i,j-1),        a_zx(i,j),        a_zx(i,j+1),        a_zx(i,j+2), &
                                                dy(i,j-1),          dy(i,j),          dy(i,j+1),          dy(i,j+2), &
                                           imask3d(i,j-1,km1), imask3d(i,j,km1), imask3d(i,j+1,km1), imask3d(i,j+2,km1), dtime)
       END DO
       END DO
    END DO
!$OMP END PARALLEL

  END SUBROUTINE advection_cosmic

!-----------------------------------------------------------------------------------------------------------------------

  !---shoud be inlined---
  REAL(TRC_KIND) PURE FUNCTION flux_upwind(v, a0, a1)
    REAL(TRC_KIND), INTENT(IN) :: v
    REAL(TRC_KIND), INTENT(IN) :: a0, a1

    IF (v > 0.0) THEN
       flux_upwind = v*a0
    ELSE
       flux_upwind = v*a1
    END IF

  END FUNCTION flux_upwind

!-----------------------------------------------------------------------------------------------------------------------

  !---shoud be inlined---
  REAL(TRC_KIND) PURE FUNCTION flux_quickest(v, a0, a1, a2, a3, d0, d1, d2, d3, m0, m1, m2, m3, dtime)
    REAL(TRC_KIND), INTENT(IN) :: v
    REAL(TRC_KIND), INTENT(IN) :: a0, a1, a2, a3
    REAL(8),        INTENT(IN) :: d0, d1, d2, d3
    INTEGER(1),     INTENT(IN) :: m0, m1, m2, m3
    REAL(4),        INTENT(IN) :: dtime

    REAL(TRC_KIND) :: au, ac, ad, a_eff
    REAL(TRC_KIND) :: v_abs
    REAL(8)        :: du, dc, dd
    REAL(TRC_KIND) :: c0, c1, c2

    IF (m1*m2==0) THEN
       flux_quickest = 0.0
       RETURN
    END IF

    IF (v > 0.0) THEN
       IF (m0==0) THEN
          flux_quickest = v * a1
          RETURN
       END IF
       au = a0
       ac = a1
       ad = a2
       du = d0*0.5 + d1
       dc = d1*0.5
       dd = d2*0.5
       v_abs = v
    ELSE
       IF (m3==0) THEN
          flux_quickest = v * a2
          RETURN
       END IF

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

    flux_quickest = v*limiter(a_eff, au, ac, ad, REAL(v_abs*dtime/(dc+dd), TRC_KIND))

  END FUNCTION flux_quickest

!-----------------------------------------------------------------------------------------------------------------------

  !---shoud be inlined---
  REAL(TRC_KIND) PURE FUNCTION limiter(a, au, ac, ad, alpha)
    REAL(TRC_KIND), INTENT(IN) :: a, au, ac, ad, alpha

    REAL(TRC_KIND) :: ar

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

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE div_correction(a)
    REAL(TRC_KIND), INTENT(INOUT) :: a(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)

    REAL(8) :: err(1:isize, 1:jsize, 1:ksize)
    REAL(8) :: err_sum

    INTEGER :: i, j, k, kl

!$OMP PARALLEL DO PRIVATE(kl)
    DO k=1, ksize
       kl=dzindex(k)
       DO j=1, jsize
       DO i=1, isize
          err(i,j,k) =  imask3d(i,j,k) * a(i,j,k) * (dvol(i,j,kl) - dvol_old(i,j,kl)  &
                      +  ( u(i,j,k)*dsx_old(i,j,kl) - u(i-1,j,k)*dsx_old(i-1,j,kl) &
                         + v(i,j,k)*dsy_old(i,j,kl) - v(i,j-1,k)*dsy_old(i,j-1,kl) &
                         + w(i,j,k)*dsz(i,j)        - w(i,j,k-1)*dsz(i,j)) * a_ebcn * dtime &
                      +  ( u_old(i,j,k)*dsx_old(i,j,kl) - u_old(i-1,j,k)*dsx_old(i-1,j,kl) &
                         + v_old(i,j,k)*dsy_old(i,j,kl) - v_old(i,j-1,k)*dsy_old(i,j-1,kl) &
                         + w_old(i,j,k)*dsz(i,j)        - w_old(i,j,k-1)*dsz(i,j)) * (1.0-a_ebcn) * dtime)
       END DO
       END DO
    END DO

    IF (use_landwater .AND. vrank==0) THEN
       DO j=1, jsize
       DO i=1, isize
          IF (lwdried(i,j)) err(i,j,ksize) =  0.0D0
       END DO
       END DO
    END IF

    err_sum = sum(err(1:isize, 1:jsize, 1:ksize))
    CALL gsum(err_sum, all=.TRUE.)

!$OMP PARALLEL DO PRIVATE(kl)
    DO k=1, ksize
       kl=dzindex(k)
       DO j=1, jsize
       DO i=1, isize
        ! a(i,j,k) = a(i,j,k) + imask3d(i,j,k) * (err(i,j,k)/dvol(i,j,kl) - err_sum/total_vol)
          a(i,j,k) = a(i,j,k) + imask3d(i,j,k) * (err(i,j,k)/dvol(i,j,kl) - REAL(err_sum/total_vol, KIND=4))
        ! down-convert err_sum/total_vol to REAL4 to keep binary-level consistency for different parallel-configuration
       END DO
       END DO
    END DO

  END SUBROUTINE div_correction

!-----------------------------------------------------------------------------------------------------------------------

  CHARACTER(32) PURE FUNCTION tracer_name(n)
    INTEGER, INTENT(IN) :: n

    IF (n <= 0 .OR. n > n_tracer) THEN
       tracer_name = '--ERROR--'
    ELSE
       tracer_name = tracer_info(n)%name
    END IF
  END FUNCTION tracer_name

!-----------------------------------------------------------------------------------------------------------------------

  INTEGER PURE FUNCTION tracer_index(varname)
    CHARACTER(*), INTENT(IN) :: varname

    INTEGER :: n

    tracer_index = 0

    DO n=1, n_tracer
       IF (trim(varname) == trim(tracer_info(n)%name)) THEN
          tracer_index = n
          RETURN
       END IF
    END DO
  END FUNCTION tracer_index

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE step_tracers
    REAL(TRC_KIND) :: fx(0:isize, 1:jsize, 1:ksize)
    REAL(TRC_KIND) :: fy(1:isize, 0:jsize, 1:ksize)
    REAL(TRC_KIND) :: fz(1:isize, 1:jsize, 0:ksize)

    REAL(TRC_KIND) :: fx_a(0:isize, 1:jsize, 1:ksize)
    REAL(TRC_KIND) :: fy_a(1:isize, 0:jsize, 1:ksize)
    REAL(TRC_KIND) :: fz_a(1:isize, 1:jsize, 0:ksize)

    REAL(TRC_KIND) :: fx_d(0:isize, 1:jsize, 1:ksize)
    REAL(TRC_KIND) :: fy_d(1:isize, 0:jsize, 1:ksize)
    REAL(TRC_KIND) :: fz_d(1:isize, 1:jsize, 0:ksize)

    REAL(TRC_KIND) :: fz_i(1:isize, 1:jsize, 0:ksize)

    REAL(TRC_KIND) :: uadv( -slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
    REAL(TRC_KIND) :: vadv(1-slv:isize+slv,  -slv:jsize+slv, 1-slv:ksize+slv)
    REAL(TRC_KIND) :: wadv(1-slv:isize+slv, 1-slv:jsize+slv,  -slv:ksize+slv)

    REAL(TRC_KIND) :: rmv(1:isize, 1:jsize, 0:1)

    REAL(TRC_KIND) :: input(  1:isize, 1:jsize, 1:ksize)
    REAL(TRC_KIND) :: resrate(1:isize, 1:jsize, 1:ksize)

    REAL(TRC_KIND) :: old_e(isize-1:isize+1, 0:jsize+1, 1:ksize, 1:n_tracer)
    REAL(TRC_KIND) :: old_w(0:2,             0:jsize+1, 1:ksize, 1:n_tracer)
    REAL(TRC_KIND) :: old_n(0:isize+1, jsize-1:jsize+1, 1:ksize, 1:n_tracer)
    REAL(TRC_KIND) :: old_s(0:isize+1, 0:2,             1:ksize, 1:n_tracer)

    REAL(8) :: u_st(-slv:ksize+slv)
    REAL(8) :: v_st(-slv:ksize+slv)

    INTEGER :: i, j, k, kl, km, n, m
    LOGICAL :: stat_uadv, stat_vadv, stat_wadv
    LOGICAL :: stat_diff
    LOGICAL :: stat_tconv

    IF (radi_tracers) THEN
!$OMP PARALLEL
       IF (radi_e .AND. icoord==ipes-1) THEN
!$OMP WORKSHARE
          old_e(:,:,:,:) = tracer(isize-1:isize+1,0:jsize+1,1:ksize,1:n_tracer)
!$OMP END WORKSHARE
       END IF
       IF (radi_w .AND. icoord==0) THEN
!$OMP WORKSHARE
          old_w(:,:,:,:) = tracer(0:2,0:jsize+1,1:ksize,1:n_tracer)
!$OMP END WORKSHARE
       END IF
       IF (radi_n .AND. jcoord==jpes-1) THEN
!$OMP WORKSHARE
          old_n(:,:,:,:) = tracer(0:isize+1,jsize-1:jsize+1,1:ksize,1:n_tracer)
!$OMP END WORKSHARE
       END IF
       IF (radi_s .AND. jcoord==0) THEN
!$OMP WORKSHARE
          old_s(:,:,:,:) = tracer(0:isize+1,0:2,1:ksize,1:n_tracer)
!$OMP END WORKSHARE
       END IF
!$OMP END PARALLEL
    END IF

    input(:,:,:)   = UNDEF
    resrate(:,:,:) = UNDEF

    u_st(:) = 0.0
    v_st(:) = 0.0
    CALL checkin('U_STOKES', u_st, axis='Z')
    CALL checkin('V_STOKES', v_st, axis='Z')

    stat_uadv = .TRUE.
    stat_vadv = .TRUE.
    stat_wadv = .TRUE.
    stat_diff = .TRUE.

    DO n=1, n_tracer
       IF (tracer_info(n)%fixed) THEN
          CALL checkin('FIX_'//trim(tracer_name(n)), tracer(:,:,:,n))
          CYCLE
       END IF

       IF (.NOT. stat_uadv) stat_uadv = tracer_info(n)%u /= tracer_info(n-1)%u
       IF (.NOT. stat_vadv) stat_vadv = tracer_info(n)%v /= tracer_info(n-1)%v
       IF (.NOT. stat_wadv) stat_wadv = tracer_info(n)%w /= tracer_info(n-1)%w

       IF (stat_uadv) THEN
!$OMP PARALLEL DO
          DO k=1-slv, ksize+slv
          DO j=1-slv, jsize+slv
          DO i= -slv, isize+slv
             uadv(i,j,k) = a_ebcn*u(i,j,k)+(1.0-a_ebcn)*u_old(i,j,k) + tracer_info(n)%u + 0.5*(u_st(k)+u_st(k-1))
          END DO
          END DO
          END DO

          stat_uadv = .FALSE.
       END IF

       IF (stat_vadv) THEN
!$OMP PARALLEL DO
          DO k=1-slv, ksize+slv
          DO j= -slv, jsize+slv
          DO i=1-slv, isize+slv
             vadv(i,j,k) = a_ebcn*v(i,j,k)+(1.0-a_ebcn)*v_old(i,j,k) + tracer_info(n)%v + 0.5*(v_st(k)+v_st(k-1))
          END DO
          END DO
          END DO

          stat_vadv = .FALSE.
       END IF

       IF (stat_wadv) THEN
!$OMP PARALLEL DO
          DO k= -slv, ksize+slv
          DO j=1-slv, jsize+slv
          DO i=1-slv, isize+slv
             wadv(i,j,k) = a_ebcn*w(i,j,k)+(1.0-a_ebcn)*w_old(i,j,k) + tracer_info(n)%w
          END DO
          END DO
          END DO

          stat_wadv = .FALSE.
       END IF

       CALL checkin('U_'//trim(tracer_name(n)), uadv, add=.TRUE., stat=stat_uadv)
       CALL checkin('V_'//trim(tracer_name(n)), vadv, add=.TRUE., stat=stat_vadv)
       CALL checkin('W_'//trim(tracer_name(n)), wadv, add=.TRUE., stat=stat_wadv)

PROFILE_BEGIN('tadv')
       SELECT CASE (tracer_info(n)%advscheme)
       CASE (0)
!$OMP PARALLEL WORKSHARE
          fx_a(:,:,:) = 0.0
          fy_a(:,:,:) = 0.0
          fz_a(:,:,:) = 0.0
!$OMP END PARALLEL WORKSHARE
       CASE (1)
          CALL advection_upwind(tracer(:,:,:,n), uadv, vadv, wadv, fx_a, fy_a, fz_a)
       CASE (2)
          CALL advection_quickest(tracer(:,:,:,n), uadv, vadv, wadv, fx_a, fy_a, fz_a)
       CASE (3)
          CALL advection_cosmic(tracer(:,:,:,n), uadv, vadv, wadv, fx_a, fy_a, fz_a)
       CASE (-3)
          !for gfortran segmentation fault issue of OpenMP (Y. Matsumura, 2022/8/23)
          CALL advection_cosmic_nottuned(tracer(:,:,:,n), uadv, vadv, wadv, fx_a, fy_a, fz_a)
       CASE DEFAULT
          CALL assert(.FALSE., "invalid tracer advection scheme for '"//trim(tracer_name(n))//"'")
       END SELECT

       CALL advection_open(trim(tracer_name(n)), tracer(:,:,:,n), uadv, vadv, wadv, fx_a, fy_a, fz_a)

       CALL checkout_flux(trim(tracer_name(n))//'_ADVC', fx_a, fy_a, fz_a)

PROFILE_END('tadv')

PROFILE_BEGIN('diff')
       IF (.NOT. stat_diff) THEN
          stat_diff = tracer_info(n)%diffh /= tracer_info(n-1)%diffh .OR. &
                      tracer_info(n)%diffv /= tracer_info(n-1)%diffv .OR. &
                      tracer_info(n)%csmag /= tracer_info(n-1)%csmag .OR. &
                      (les_scheme /= 0 .AND. tracer_info(n)%cles  /= tracer_info(n-1)%cles)
       END IF

       IF (stat_diff) THEN
!$OMP PARALLEL DO
          DO k=0, ksize+1
          DO j=0, jsize+1
          DO i=0, isize+1
             diff(i,j,k,1) = tracer_info(n)%diffh
             diff(i,j,k,2) = tracer_info(n)%diffv
          END DO
          END DO
          END DO

          CALL checkin('DIFF_H', diff(:,:,:,1), add=.TRUE.)
          CALL checkin('DIFF_V', diff(:,:,:,2), add=.TRUE.)

          CALL checkin('DIFFH', diff(:,:,:,1), add=.TRUE.)
          CALL checkin('DIFFV', diff(:,:,:,2), add=.TRUE.)

          IF (tracer_info(n)%csmag /= 0.0) THEN
!$OMP PARALLEL DO
             DO k=0, ksize+1
             DO j=0, jsize+1
             DO i=0, isize+1
                diff(i,j,k,1) = diff(i,j,k,1) + tracer_info(n)%csmag**2 * delta2(i,j,k) * sqrt2 * srate(i,j,k)
             END DO
             END DO
             END DO
          END IF

          IF (les_scheme /= 0) THEN
!$OMP PARALLEL DO
             DO k=0, ksize+1
             DO j=0, jsize+1
             DO i=0, isize+1
                diff(i,j,k,1) = diff(i,j,k,1) + visc(i,j,k,3)*tracer_info(n)%cles
                diff(i,j,k,2) = diff(i,j,k,2) + visc(i,j,k,3)*tracer_info(n)%cles
             END DO
             END DO
             END DO
          END IF

          stat_diff = .FALSE.
       END IF

       CALL checkin('DIFFH_'//trim(tracer_name(n)), diff(:,:,:,1), add=.TRUE., stat=stat_diff)
       CALL checkin('DIFFV_'//trim(tracer_name(n)), diff(:,:,:,2), add=.TRUE., stat=stat_diff)

       CALL checkout('DIFFH_'//trim(tracer_name(n)), diff(:,:,:,1)) ! horizontal diffusivity coefficient for tracer [m^2/s]
       CALL checkout('DIFFV_'//trim(tracer_name(n)), diff(:,:,:,2)) ! vertical diffusivity coefficient for tracer [m^2/s]

       SELECT CASE (tracer_info(n)%diffscheme)
       CASE (0)
!$OMP PARALLEL WORKSHARE
          fx_d(:,:,:) = 0.0
          fy_d(:,:,:) = 0.0
          fz_d(:,:,:) = 0.0
!$OMP END PARALLEL WORKSHARE

       CASE (1)
          CALL diffusion_iso(tracer(:,:,:,n), diff(:,:,:,1), fx_d, fy_d, fz_d)
       CASE (2)
          CALL diffusion_h(tracer(:,:,:,n), diff(:,:,:,1), fx_d, fy_d)
          CALL diffusion_v(tracer(:,:,:,n), diff(:,:,:,2), fz_d)
    !  CASE (3)
    !     CALL diffusion_idp(tracer(:,:,:,n), diff, n_x, n_y, n_z, fx_d, fy_d, fz_d)
       CASE (4)
          CALL diffusion_h(   tracer(:,:,:,n), diff(:,:,:,1), fx_d, fy_d)
          CALL diffusion_vimp(tracer(:,:,:,n), diff(:,:,:,2), fz_d)

       CASE DEFAULT
          CALL assert(.FALSE., "invalid tracer diffusion scheme for '"//trim(tracer_name(n))//"'")
       END SELECT

       IF (use_gls) THEN
          IF (tracer_info(n)%cgls < 0) THEN !for TKE and GLS, cgls is set negative and use KM instead of KH
             CALL diffusion_vimp(tracer(:,:,:,n), (-tracer_info(n)%cgls)*visc(:,:,:,3), fz_i)
          ELSE
             CALL diffusion_vimp(tracer(:,:,:,n), diff(:,:,:,3), fz_i, factor=REAL(tracer_info(n)%cgls, TRC_KIND))
          END IF
!$OMP PARALLEL WORKSHARE
          fz_d(:,:,:) = fz_d(:,:,:) + fz_i(:,:,:)
!$OMP END PARALLEL WORKSHARE
       END IF
PROFILE_END('diff')

       CALL checkout_flux(trim(tracer_name(n))//'_DIFF', fx_d, fy_d, fz_d)

!$OMP PARALLEL WORKSHARE
       fx(:,:,:) = fx_a(:,:,:) + fx_d(:,:,:)
       fy(:,:,:) = fy_a(:,:,:) + fy_d(:,:,:)
       fz(:,:,:) = fz_a(:,:,:) + fz_d(:,:,:)
!$OMP END PARALLEL WORKSHARE

       CALL checkout_flux(tracer_name(n), fx, fy, fz)

       IF (tracer_info(n)%rmvsfc) THEN
          rmv(:,:,1) = 0.0

          DO j=1, jsize
          DO i=1, isize
             IF (surface_flag(i,j)) THEN
                k = surface_k(i,j)
                rmv(i,j,1) = max(wadv(i,j,k), 0.0) * tracer(i,j,k,n)
             END IF
          END DO
          END DO

          CALL vsum(rmv(:,:,1))
          CALL checkout('SFCRMV_'//trim(tracer_name(n)), rmv(:,:,1)) ! removed tracer at the surface per unit area and unit time [(TR) m/s]
       END IF

       IF (tracer_info(n)%rmvbtm) THEN
          rmv(:,:,0) = 0.0
          DO j=1, jsize
          DO i=1, isize
             IF (bottom_flag(i,j)) THEN
                k = bottom_k(i,j)
                rmv(i,j,0) = -min(wadv(i,j,k), 0.0) * tracer(i,j,k,n)
             END IF
          END DO
          END DO

          CALL vsum(rmv(:,:,0))
          CALL checkout('BTMRMV_'//trim(tracer_name(n)), rmv(:,:,0)) ! removed tracer at the bottom per unit area and unit time [(TR) m/s]
       END IF

!$OMP PARALLEL DO PRIVATE(kl)
       DO k=1, ksize
          kl=dzindex(k)
          DO j=1, jsize
          DO i=1, isize
             tracer(i,j,k,n) = (tracer(i,j,k,n)*dvol_old(i,j,kl) &
                  + ( fx(i-1,j,k)*dsx_old(i-1,j,kl) - fx(i,j,k)*dsx_old(i,j,kl) &
                    + fy(i,j-1,k)*dsy_old(i,j-1,kl) - fy(i,j,k)*dsy_old(i,j,kl) &
                    + fz(i,j,k-1)*dsz(i,j)          - fz(i,j,k)*dsz(i,j)) * dtime)/dvol(i,j,kl)
          END DO
          END DO
       END DO

       IF (tracer_info(n)%convtgt /= "") THEN
          IF (.NOT. ASSOCIATED(tracer_info(n)%convflag)) THEN
             ALLOCATE(tracer_info(n)%convflag(0:isize+1, 0:jsize+1, 0:ksize+1))
             ALLOCATE(tracer_info(n)%convrate(1:isize,   1:jsize,   1:ksize))

             tracer_info(n)%n_convtgt = tracer_index(tracer_info(n)%convtgt)
             CALL assert(tracer_info(n)%n_convtgt /= 0, "undefined tracer convert target '" // trim(tracer_info(n)%convtgt) // "'")
          END IF

!$OMP PARALLEL WORKSHARE
          tracer_info(n)%convrate(:,:,:) = 0.0
!$OMP END PARALLEL WORKSHARE

          CALL checkin("TCONVFLAG_"//tracer_name(n), tracer_info(n)%convflag, stat=stat_tconv)
          IF (.NOT. stat_tconv) CALL checkin("TCONVFLAG", tracer_info(n)%convflag, stat=stat_tconv)

          IF (stat_tconv) THEN
!$OMP PARALLEL DO
             DO k=1, ksize
             DO j=1, jsize
             DO i=1, isize
                IF (tracer_info(n)%convflag(i,j,k)==1 .OR. tracer_info(n)%convflag(i,j,k)==3) THEN
                   IF (tracer_info(n)%convflag(i-1,j,k)==0 .AND. fx_a(i-1,j,k) > 0.0) tracer_info(n)%convrate(i,j,k) = tracer_info(n)%convrate(i,j,k) + fx_a(i-1,j,k)*dsx_old(i-1,j,k)
                   IF (tracer_info(n)%convflag(i+1,j,k)==0 .AND. fx_a(i,  j,k) < 0.0) tracer_info(n)%convrate(i,j,k) = tracer_info(n)%convrate(i,j,k) - fx_a(i,  j,k)*dsx_old(i,  j,k)
                   IF (tracer_info(n)%convflag(i,j-1,k)==0 .AND. fy_a(i,j-1,k) > 0.0) tracer_info(n)%convrate(i,j,k) = tracer_info(n)%convrate(i,j,k) + fy_a(i,j-1,k)*dsy_old(i,j-1,k)
                   IF (tracer_info(n)%convflag(i,j+1,k)==0 .AND. fy_a(i,j,  k) < 0.0) tracer_info(n)%convrate(i,j,k) = tracer_info(n)%convrate(i,j,k) - fy_a(i,j+1,k)*dsy_old(i,j,  k)
                   IF (tracer_info(n)%convflag(i,j,k-1)==0 .AND. fz_a(i,j,k-1) > 0.0) tracer_info(n)%convrate(i,j,k) = tracer_info(n)%convrate(i,j,k) + fz_a(i,j,k-1)*dsz(i,j)
                   IF (tracer_info(n)%convflag(i,j,k+1)==0 .AND. fz_a(i,j,k)   < 0.0) tracer_info(n)%convrate(i,j,k) = tracer_info(n)%convrate(i,j,k) - fz_a(i,j,k+1)*dsz(i,j)
                END IF
                IF (tracer_info(n)%convflag(i,j,k)==2 .OR. tracer_info(n)%convflag(i,j,k)==3) THEN
                   IF (tracer_info(n)%convflag(i-1,j,k)==0 .AND. fx_d(i-1,j,k) > 0.0) tracer_info(n)%convrate(i,j,k) = tracer_info(n)%convrate(i,j,k) + fx_d(i-1,j,k)*dsx_old(i-1,j,k)
                   IF (tracer_info(n)%convflag(i+1,j,k)==0 .AND. fx_d(i,  j,k) < 0.0) tracer_info(n)%convrate(i,j,k) = tracer_info(n)%convrate(i,j,k) - fx_d(i,  j,k)*dsx_old(i,  j,k)
                   IF (tracer_info(n)%convflag(i,j-1,k)==0 .AND. fy_d(i,j-1,k) > 0.0) tracer_info(n)%convrate(i,j,k) = tracer_info(n)%convrate(i,j,k) + fy_d(i,j-1,k)*dsy_old(i,j-1,k)
                   IF (tracer_info(n)%convflag(i,j+1,k)==0 .AND. fy_d(i,j,  k) < 0.0) tracer_info(n)%convrate(i,j,k) = tracer_info(n)%convrate(i,j,k) - fy_d(i,j+1,k)*dsy_old(i,j,  k)
                   IF (tracer_info(n)%convflag(i,j,k-1)==0 .AND. fz_d(i,j,k-1) > 0.0) tracer_info(n)%convrate(i,j,k) = tracer_info(n)%convrate(i,j,k) + fz_d(i,j,k-1)*dsz(i,j)
                   IF (tracer_info(n)%convflag(i,j,k+1)==0 .AND. fz_d(i,j,k)   < 0.0) tracer_info(n)%convrate(i,j,k) = tracer_info(n)%convrate(i,j,k) - fz_d(i,j,k+1)*dsz(i,j)
                END IF
                tracer_info(n)%convrate(i,j,k) = tracer_info(n)%convrate(i,j,k) / dvol(i,j,k)
             END DO
             END DO
             END DO
          END IF
       ELSE
          tracer_info(n)%n_convtgt = 0
       END IF

       IF (tracer_info(n)%correctdiv) CALL div_correction(tracer(:,:,:,n))

       IF (tracer_info(n)%rmvsfc) THEN
          DO j=1, jsize
          DO i=1, isize
             IF (use_landwater .AND. vrank==0) THEN
                IF (lwdried(i,j)) CYCLE
             END IF

             IF (surface_flag(i,j)) THEN
                k = surface_k(i,j)
                tracer(i,j,k,n) = tracer(i,j,k,n) - rmv(i,j,1)*dsz(i,j)*dtime/dvol(i,j,k)
             END IF
          END DO
          END DO
       END IF

       IF (tracer_info(n)%rmvbtm) THEN
          DO j=1, jsize
          DO i=1, isize
             IF (use_landwater .AND. vrank==0) THEN
                IF (lwdried(i,j)) CYCLE
             END IF

             IF (bottom_flag(i,j)) THEN
                k = bottom_k(i,j)
                tracer(i,j,k,n) = tracer(i,j,k,n) - rmv(i,j,0)*dsz(i,j)*dtime/dvol(i,j,k)
             END IF
          END DO
          END DO
       END IF

       CALL forcing(n)

       CALL openboundary(n)
PROFILE_BEGIN('tb3d')
       CALL update_tracer_boundary(n)
PROFILE_END('tb3d')
    END DO

    DO n=1, n_tracer
       m = tracer_info(n)%n_convtgt
       IF (m==0) CYCLE

!$OMP PARALLEL DO
       DO k=1, ksize
       DO j=1, jsize
       DO i=1, isize
          tracer(i,j,k, n) = tracer(i,j,k, n) - tracer_info(n)%convrate(i,j,k)*dtime
          tracer(i,j,k, m) = tracer(i,j,k, m) + tracer_info(n)%convrate(i,j,k)*dtime
       END DO
       END DO
       END DO
       CALL update_tracer_boundary(n)
       CALL update_tracer_boundary(m)
       CALL checkout("TCONVRATE_"//trim(tracer_name(n))//"_"//trim(tracer_name(m)), tracer_info(n)%convrate)
    END DO

!$OMP PARALLEL WORKSHARE
    diff(:,:,:,:) = 0.0D0
!$OMP END PARALLEL WORKSHARE

  CONTAINS
    SUBROUTINE forcing(index)
      USE parameters, ONLY: rho_0, cp
      INTEGER, INTENT(IN) :: index

      LOGICAL :: stat
      INTEGER :: i, j, k

      REAL(8) :: incr(isize, jsize, ksize)
      REAL(8) :: tmp
      LOGICAL :: req


      req = require_checkout('FLXINCR_'//trim(tracer_name(index)))
      IF (req) THEN
!$OMP PARALLEL WORKSHARE
         incr(:,:,:) = 0.0
!$OMP END PARALLEL WORKSHARE
      END IF

      ! surface flux, upward positive
      stat = .FALSE.
      IF (index == tracer_index_t) THEN
         CALL checkin('SFCFLUX_HEAT', input(:,:,1), stat)       ! [w/m^2]
         IF (stat) input(:,:,1) = input(:,:,1) / (rho_0 * cp)   ! [K m/s]
      ELSE IF (index == tracer_index_s) THEN
         CALL checkin('SFCFLUX_SALT', input(:,:,1), stat)       ! [g/(m^2 s)]
         IF (stat) input(:,:,1) = input(:,:,1) / rho_0          ! [psu m/s]
      END IF
      CALL checkin('SFCFLUX_'//trim(tracer_name(index)), input(:,:,1), stat, add=stat.EQV..TRUE.)

      IF (stat) THEN
!$OMP PARALLEL DO PRIVATE(tmp, k)
         DO j=1, jsize
         DO i=1, isize
            IF (use_landwater .AND. vrank==0) THEN
               IF (lwdried(i,j)) CYCLE
            END IF

            IF (surface_flag(i,j)) THEN
               k = surface_k(i,j)
               IF (tracer_info(index)%maskflux .AND. .NOT. (vrank==0 .AND. k==ksize)) CYCLE
               tmp = imask3d(i,j,k)*input(i,j,1)/dz_star(i,j,k)
               tracer(i,j,k,index) = tracer(i,j,k,index) - tmp*dtime
               IF (req) incr(i,j,k) = incr(i,j,k) - tmp
            END IF
         END DO
         END DO
      END IF

      ! bottom flux,  upward positive
      stat = .FALSE.
      IF (index == tracer_index_t) THEN
         CALL checkin('BTMFLUX_HEAT', input(:,:,1), stat)       ! [w/m^2]
         IF (stat) input(:,:,1) = input(:,:,1) / (rho_0 * cp)   ! [K m/s]
      ELSE IF (index == tracer_index_s) THEN
         CALL checkin('BTMFLUX_SALT', input(:,:,1), stat)       ! [g/(m^2 s)]
         IF (stat) input(:,:,1) = input(:,:,1) / rho_0          ! [psu m/s]
      END IF
      CALL checkin('BTMFLUX_'//trim(tracer_name(index)), input(:,:,1), stat, add=stat.EQV..TRUE.)

      IF (stat) THEN
!$OMP PARALLEL DO PRIVATE(tmp, k)
         DO j=1, jsize
         DO i=1, isize
            IF (bottom_flag(i,j)) THEN
               k = bottom_k(i,j)
               tmp = imask3d(i,j,k)*input(i,j,1)/dz_star(i,j,k)
               tracer(i,j,k,index) = tracer(i,j,k,index) + tmp*dtime
               IF (req) incr(i,j,k) = incr(i,j,k) + tmp
            END IF
         END DO
         END DO
      END IF

      IF (req) CALL checkout('FLXINCR_'//trim(tracer_name(index)), incr)

      ! source/sink
      stat = .FALSE.
      IF (index == tracer_index_t) THEN
         CALL checkin('SOURCE_HEAT', input, stat)              ! [W/m^3]
         IF (stat) THEN
!$OMP PARALLEL WORKSHARE
            input(:,:,:) = input(:,:,:) / (rho_0 * cp)         ! [K/s]
!$OMP END PARALLEL WORKSHARE
         END IF
      ELSE IF (index == tracer_index_s) THEN
         CALL checkin('SOURCE_SALT', input, stat)              ! [g/(m^3 s)]
         IF (stat) THEN
!$OMP PARALLEL WORKSHARE
            input(:,:,:) = input(:,:,:) / rho_0                ! [psu/s]
!$OMP END PARALLEL WORKSHARE
         END IF
      END IF
      CALL checkin('SOURCE_'//trim(tracer_name(index)), input, stat, add=stat.EQV..TRUE.)

      IF (stat) THEN
!$OMP PARALLEL DO
         DO k=1, ksize
         DO j=1, jsize
         DO i=1, isize
            IF (input(i,j,k)==UNDEF) CYCLE
            tracer(i,j,k,index) = tracer(i,j,k,index) + imask3d(i,j,k)*input(i,j,k)*dtime
         END DO
         END DO
         END DO
      END IF

      ! restore
      req = require_checkout('RESINCR_'//trim(tracer_name(index)))
      IF (req) THEN
!$OMP PARALLEL WORKSHARE
         incr(:,:,:) = 0.0
!$OMP END PARALLEL WORKSHARE
      END IF

      CALL checkin('RESTORE_'//trim(tracer_name(index)), input, stat)
      IF (stat) THEN
         CALL checkin('RESRATE_'//trim(tracer_name(index)), resrate, stat)
         IF (.NOT. stat) CALL checkin('RESRATE_TRACERS', resrate, stat)
         IF (.NOT. stat) CALL checkin('RESRATE',         resrate, stat)
         CALL assert(stat, "restoring rate (RESRATE) for '"//trim(tracer_name(index))//"' is not specified.")

!$OMP PARALLEL DO PRIVATE(tmp)
         DO k=1, ksize
         DO j=1, jsize
         DO i=1, isize
            IF (input(i,j,k)==UNDEF) CYCLE

            tmp = imask3d(i,j,k)*(input(i,j,k)-tracer(i,j,k,index))*resrate(i,j,k)
            tracer(i,j,k,index) = tracer(i,j,k,index) + tmp*dtime
            IF (req) incr(i,j,k) = incr(i,j,k) + tmp
         END DO
         END DO
         END DO
      END IF

      CALL checkin('RESTORE_SFC_'//trim(tracer_name(index)), input(:,:,1), stat)
      IF (stat) THEN
         CALL checkin('RESRATE_SFC_'//trim(tracer_name(index)), resrate(:,:,1), stat)
         IF (.NOT. stat) CALL checkin('RESRATE_SURFACE', resrate(:,:,1), stat)
         CALL assert(stat, "restoring rate (RESRATE) for 'SFC_"//trim(tracer_name(index))//"' is not specified.")

!$OMP PARALLEL DO PRIVATE(tmp, k)
         DO j=1, jsize
         DO i=1, isize
            IF (input(i,j,1)==UNDEF) CYCLE
            IF (surface_flag(i,j)) THEN
               k = surface_k(i,j)
               tmp = imask3d(i,j,k)*(input(i,j,1)-tracer(i,j,k,index))*resrate(i,j,1)
               tracer(i,j,k,index) = tracer(i,j,k,index) + tmp*dtime
               IF (req) incr(i,j,k) = incr(i,j,k) + tmp
            END IF
         END DO
         END DO
      END IF

      CALL checkin('RESTORE_BTM_'//trim(tracer_name(index)), input(:,:,1), stat)
      IF (stat) THEN
         CALL checkin('RESRATE_BTM_'//trim(tracer_name(index)), resrate(:,:,1), stat)
         IF (.NOT. stat) CALL checkin('RESRATE_BOTTOM', resrate(:,:,1), stat)
         CALL assert(stat, "restoring rate (RESRATE) for 'BTM_"//trim(tracer_name(index))//"' is not specified.")

!$OMP PARALLEL DO PRIVATE(tmp, k)
         DO j=1, jsize
         DO i=1, isize
            IF (input(i,j,1)==UNDEF) CYCLE
            IF (bottom_flag(i,j)) THEN
               k = bottom_k(i,j)
               tmp = imask3d(i,j,k)*(input(i,j,1)-tracer(i,j,k,index))*resrate(i,j,1)
               tracer(i,j,k,index) = tracer(i,j,k,index) + tmp*dtime
               IF (req) incr(i,j,k) = incr(i,j,k) + tmp
            END IF
         END DO
         END DO
      END IF

      IF (req) CALL checkout('RESINCR_'//trim(tracer_name(index)), incr)

      ! fix
      CALL checkin('FIX_'//trim(tracer_name(index)), input, stat)
      IF (stat) THEN
!$OMP PARALLEL DO
         DO k=1, ksize
         DO j=1, jsize
         DO i=1, isize
            IF (input(i,j,k) /= UNDEF) tracer(i,j,k,index) = input(i,j,k)
         END DO
         END DO
         END DO
      END IF

      CALL checkin('FIX_SFC_'//trim(tracer_name(index)), input(:,:,1), stat)
      IF (stat) THEN
!$OMP PARALLEL DO PRIVATE(k)
         DO j=1, jsize
         DO i=1, isize
            IF (surface_flag(i,j)) THEN
               k = surface_k(i,j)
               IF (input(i,j,1) /= UNDEF) tracer(i,j,k,index) = input(i,j,1)
            END IF
         END DO
         END DO
      END IF

      CALL checkin('FIX_BTM_'//trim(tracer_name(index)), input(:,:,1), stat)
      IF (stat) THEN
!$OMP PARALLEL DO PRIVATE(k)
         DO j=1, jsize
         DO i=1, isize
            IF (bottom_flag(i,j)) THEN
               k = bottom_k(i,j)
               IF (input(i,j,1) /= UNDEF) tracer(i,j,k,index) = input(i,j,1)
            END IF
         END DO
         END DO
      END IF

    END SUBROUTINE forcing

!-----------------------------------------------------------------------------------------------------------------------

    SUBROUTINE openboundary(index)
      INTEGER, INTENT(IN) :: index

      CHARACTER(32) :: name
      REAL(8) :: tmp_xz(isize,ksize)
      REAL(8) :: tmp_yz(jsize,ksize)
      REAL(8) :: a, bx, by, cx, cy, d
      INTEGER :: i, j, k
      LOGICAL :: stat

      name = tracer_name(index)

      IF (open_e) THEN
         tmp_yz(:,:) = UNDEF
         CALL checkin(trim(name)//'_EAST', tmp_yz, section='YZ', time=t_current+dtime)

         IF (icoord==ipes-1) THEN
            IF (radi_tracers .AND. radi_e) THEN
!$OMP PARALLEL DO PRIVATE(a, bx, by, cx, cy, d)
               DO k=1, ksize
               DO j=1, jsize
                  IF (.NOT. lmask3d(isize,j,k)) CYCLE

                  ! a = tracer(isize,j,k,index) - old_e(isize,j,k,index)

                  ! bx = old_e(isize-1,j,k,index) - old_e(isize,j,k,index)

                  ! IF (radi_oblique) THEN
                  !    by = 0.5*(old_e(isize, j-1, k, index) - old_e(isize,j+1,k,index))
                  ! ELSE
                  !    by = 0.0
                  ! END IF

                  ! IF (a*bx > 0.0) THEN
                  !    cx = a*bx / (bx**2 + by**2)
                  !    cy = a*by / (bx**2 + by**2)
                  !    d = radi_nudge_e(1)
                  ! ELSE
                  !    cx = 0.0
                  !    cy = 0.0
                  !    d = radi_nudge_e(2)
                  ! END IF

                  ! tracer(isize+1,j,k,index) = old_e(isize+1,j,k,index) + cx*tracer(isize,j,k,index)

                  ! IF (cy > 0.0) THEN
                  !    tracer(isize+1,j,k,index) = tracer(isize+1,j,k,index) + cy*old_e(isize+1,j-1,k,index)
                  ! ELSE
                  !    tracer(isize+1,j,k,index) = tracer(isize+1,j,k,index) - cy*old_e(isize+1,j+1,k,index)
                  ! END IF

                  ! tracer(isize+1,j,k,index) = tracer(isize+1,j,k,index) / (1.0 + cx + abs(cy))

                  ! IF (tmp_yz(j,k) /= UNDEF) tracer(isize+1,j,k,index) = (1.0-d)*tracer(isize+1,j,k,index) + d*tmp_yz(j,k)

                  ! tmp_yz(j,k) = tracer(isize+1,j,k,index)
                  tracer(isize,j,k,index) = tracer(isize-1,j,k,index)
                  tracer(isize+1,j,k,index) = tracer(isize,j,k,index)
               END DO
               END DO
            ELSE
!$OMP PARALLEL DO
!!!!!!!!!!!!!!!YJKIM
               DO k=1, ksize
               DO j=1, jsize
                  IF (tmp_yz(j,k) == UNDEF) tmp_yz(j,k) = tracer(isize,j,k,index)
!                  tracer(isize,j,k,index) = tmp_yz(j,k)
                  tracer(isize+1,j,k,index) = tmp_yz(j,k)
!                  tracer(isize+2,j,k,index) = tmp_yz(j,k)
               END DO
               END DO
           END IF
         ELSE
            tmp_yz(:,:) = 0.0
         END IF

         CALL checkout(trim(name)//'_EAST', tmp_yz, section='YZ')
      END IF

      IF (open_w) THEN
         tmp_yz(:,:) = UNDEF
         CALL checkin(trim(name)//'_WEST', tmp_yz, section='YZ', time=t_current+dtime)

         IF (icoord==0) THEN
            IF (radi_tracers .AND. radi_w) THEN
!$OMP PARALLEL DO PRIVATE(a, bx, by, cx, cy, d)
               DO k=1, ksize
               DO j=1, jsize
                  IF (.NOT. lmask3d(1,j,k)) CYCLE

                  a = tracer(1,j,k,index) - old_w(1,j,k,index)

                  bx = old_w(2,j,k,index) - old_w(1,j,k,index)

                  IF (radi_oblique) THEN
                     by = 0.5*(old_w(1,j-1,k,index) - old_w(1,j+1,k,index))
                  ELSE
                     by = 0.0
                  END IF

                  IF (a*bx > 0.0) THEN
                     cx = a*bx / (bx**2 + by**2)
                     cy = a*by / (bx**2 + by**2)
                     d = radi_nudge_w(1)
                  ELSE
                     cx = 0.0
                     cy = 0.0
                     d = radi_nudge_w(2)
                  END IF

                  tracer(0,j,k,index) = old_w(0,j,k,index) + cx*tracer(1,j,k,index)

                  IF (cy > 0.0) THEN
                     tracer(0,j,k,index) = tracer(0,j,k,index) + cy*old_w(0,j-1,k,index)
                  ELSE
                     tracer(0,j,k,index) = tracer(0,j,k,index) - cy*old_w(0,j+1,k,index)
                  END IF

                  tracer(0,j,k,index) = tracer(0,j,k,index) / (1.0 + cx + abs(cy))

                  IF (tmp_yz(j,k) /= UNDEF) tracer(0,j,k,index) = (1.0-d)*tracer(0,j,k,index) + d*tmp_yz(j,k)

                  tmp_yz(j,k) = tracer(0,j,k,index)
               END DO
               END DO
            ELSE
!$OMP PARALLEL DO
!!!!!!!!!!!!!!!YJKIM
               DO k=1, ksize
               DO j=1, jsize
                  IF (tmp_yz(j,k) == UNDEF) tmp_yz(j,k) = tracer(1,j,k,index)
!                  tracer(1,j,k,index) = tmp_yz(j,k)
                  tracer(0,j,k,index) = tmp_yz(j,k)
!                  tracer(-1,j,k,index) = tmp_yz(j,k)
               END DO
               END DO
            END IF
         ELSE
            tmp_yz(:,:) = 0.0
         END IF

         CALL checkout(trim(name)//'_WEST', tmp_yz, section='YZ')
      END IF


      IF (open_n) THEN
         tmp_xz(:,:) = UNDEF
         CALL checkin(trim(name)//'_NORTH', tmp_xz, section='XZ', time=t_current+dtime)

         IF (jcoord==jpes-1) THEN
            IF (radi_tracers .AND. radi_n) THEN
!$OMP PARALLEL DO PRIVATE(a, bx, by, cx, cy, d)
               DO k=1, ksize
               DO i=1, isize
                  IF (.NOT. lmask3d(i,jsize,k)) CYCLE

                  ! a = tracer(i,jsize,k,index) - old_n(i,jsize,k,index)

                  ! by = old_n(i,jsize-1,k,index) - old_n(i,jsize,k,index)

                  ! IF (radi_oblique) THEN
                  !    bx = 0.5*(old_n(i-1,jsize,k,index) - old_n(i+1,jsize,k,index))
                  ! ELSE
                  !    bx = 0.0
                  ! END IF

                  ! IF (a*by > 0.0) THEN
                  !    cx = a*bx / (bx**2 + by**2)
                  !    cy = a*by / (bx**2 + by**2)
                  !    d = radi_nudge_n(1)
                  ! ELSE
                  !    cx = 0.0
                  !    cy = 0.0
                  !    d = radi_nudge_n(2)
                  ! END IF

                  ! tracer(i,jsize+1,k,index) = old_n(i,jsize+1,k,index) + cy*tracer(i,jsize,k,index)

                  ! IF (cx > 0.0) THEN
                  !    tracer(i,jsize+1,k,index) = tracer(i,jsize+1,k,index) + cx*old_n(i-1,jsize+1,k,index)
                  ! ELSE
                  !    tracer(i,jsize+1,k,index) = tracer(i,jsize+1,k,index) - cx*old_n(i+1,jsize+1,k,index)
                  ! END IF

                  ! tracer(i,jsize+1,k,index) = tracer(i,jsize+1,k,index) / (1.0 + cy + abs(cx))

                  ! IF (tmp_xz(i,k) /= UNDEF) tracer(i,jsize+1,k,index) = (1.0-d)*tracer(i,jsize+1,k,index) + d*tmp_xz(i,k)

                  ! tmp_xz(i,k) = tracer(i,jsize+1,k,index)
                  tracer(i,jsize,k,index)=tracer(i,jsize-1,k,index)
                  tracer(i,jsize+1,k,index)=tracer(i,jsize,k,index)
               END DO
               END DO
            ELSE
!$OMP PARALLEL DO
!!!!!!!!!!!!!!!YJKIM
               DO k=1, ksize
               DO i=1, isize
                  IF (tmp_xz(i,k) == UNDEF) tmp_xz(i,k) = tracer(i,jsize,k,index)
!                  tracer(i,jsize,k,index) = tmp_xz(i,k)
                  tracer(i,jsize+1,k,index) = tmp_xz(i,k)
!                  tracer(i,jsize+2,k,index) = tmp_xz(i,k)
               END DO
               END DO
            END IF
         ELSE
            tmp_xz(:,:) = 0.0
         END IF

         CALL checkout(trim(name)//'_NORTH', tmp_xz, section='XZ')
      END IF

      IF (open_s) THEN
         tmp_xz(:,:) = UNDEF
         CALL checkin(trim(name)//'_SOUTH', tmp_xz, section='XZ', time=t_current+dtime)

         IF (jcoord==0) THEN
            IF (radi_tracers .AND. radi_s) THEN
!$OMP PARALLEL DO PRIVATE(a, bx, by, cx, cy, d)
               DO k=1, ksize
               DO i=1, isize
                  IF (.NOT. lmask3d(i,1,k)) CYCLE

                  ! a = tracer(i,1,k,index) - old_s(i,1,k,index)

                  ! by = old_s(i,2,k,index) - old_s(i,1,k,index)

                  ! IF (radi_oblique) THEN
                  !    bx = 0.5*(old_s(i-1,1,k,index) - old_s(i+1,1,k,index))
                  ! ELSE
                  !    bx = 0.0
                  ! END IF

                  ! IF (a*by > 0.0) THEN
                  !    cx = a*bx / (bx**2 + by**2)
                  !    cy = a*by / (bx**2 + by**2)
                  !    d = radi_nudge_s(1)
                  ! ELSE
                  !    cx = 0.0
                  !    cy = 0.0
                  !    d = radi_nudge_s(2)
                  ! END IF

                  ! tracer(i,0,k,index) = old_s(i,0,k,index) + cy*tracer(i,1,k,index)

                  ! IF (cx > 0.0) THEN
                  !    tracer(i,0,k,index) = tracer(i,0,k,index) + cx*old_s(i-1,0,k,index)
                  ! ELSE
                  !    tracer(i,0,k,index) = tracer(i,0,k,index) - cx*old_s(i+1,0,k,index)
                  ! END IF

                  ! tracer(i,0,k,index) = tracer(i,0,k,index) / (1.0 + cy + abs(cx))

                  ! IF (tmp_xz(i,k) /= UNDEF) tracer(i,0,k,index) = (1.0-d)*tracer(i,0,k,index) + d*tmp_xz(i,k)

                  ! tmp_xz(i,k) = tracer(i,0,k,index)
                  tracer(i,2,k,index) = tracer(i,3,k,index)
                  tracer(i,1,k,index) = tracer(i,2,k,index)
                  tracer(i,0,k,index) = tracer(i,1,k,index)
               END DO
               END DO
            ELSE
!$OMP PARALLEL DO
!!!!!!!!!!!!!!!YJKIM               
               DO k=1, ksize
               DO i=1, isize
                  IF (tmp_xz(i,k) == UNDEF) tmp_xz(i,k) = tracer(i,1,k,index)
!                  tracer(i,1,k,index) = tmp_xz(i,k)
                  tracer(i,0,k,index) = tmp_xz(i,k)
!                  tracer(i,-1,k,index) = tmp_xz(i,k)
               END DO
               END DO
            END IF
         ELSE
            tmp_xz(:,:) = 0.0
         END IF

         CALL checkout(trim(name)//'_SOUTH', tmp_xz, section='XZ')
      END IF

    END SUBROUTINE openboundary

  END SUBROUTINE step_tracers

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE checkout_tracers
    LOGICAL :: stat
    INTEGER :: i, j, k, n

    REAL(TRC_KIND) :: tmp(1:isize,1:jsize,1:ksize)

    DO n=1, n_tracer
       CALL checkout(trim(tracer_name(n)), tracer(:,:,:,n)) ! tracer concentration

       IF (require_checkout(trim(tracer_name(n))//'_VOL')) THEN
!$OMP PARALLEL DO
          DO k=1, ksize
          DO j=1, jsize
          DO i=1, isize
             tmp(i,j,k) = imask3d(i,j,k)*tracer(i,j,k,n)*dvol(i,j,k)
          END DO
          END DO
          END DO
        END IF

        CALL checkout(trim(tracer_name(n))//'_VOL', tmp) !volume integral of tracer
     END DO

  END SUBROUTINE checkout_tracers

  SUBROUTINE checkout_flux(varname, fx, fy, fz)
    CHARACTER(*),   INTENT(IN) :: varname
    REAL(TRC_KIND), INTENT(IN) :: fx(0:isize, 1:jsize, 1:ksize)
    REAL(TRC_KIND), INTENT(IN) :: fy(1:isize, 0:jsize, 1:ksize)
    REAL(TRC_KIND), INTENT(IN) :: fz(1:isize, 1:jsize, 0:ksize)

    CALL checkout('FX_'//trim(varname), fx)
    CALL checkout('FY_'//trim(varname), fy)
    CALL checkout('FZ_'//trim(varname), fz)

    CALL checkout_div(trim(varname), fx, fy, fz)

    ! tracer transport (flux * sectin area) [tracer * m^3 / s]
    IF (require_checkout('TX_'//trim(varname))) CALL checkout('TX_'//trim(varname), fx*dsx_old(0:isize,1:jsize,1:ksize))
    IF (require_checkout('TY_'//trim(varname))) CALL checkout('TY_'//trim(varname), fy*dsy_old(1:isize,0:jsize,1:ksize))

  END SUBROUTINE checkout_flux

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE diffusion_iso(a, kappa, fx, fy, fz)
    REAL(TRC_KIND), INTENT(IN)  :: a(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
    REAL(TRC_KIND), INTENT(IN)  :: kappa(0:isize+1, 0:jsize+1, 0:ksize+1)
    REAL(TRC_KIND), INTENT(OUT) :: fx(0:isize, 1:jsize, 1:ksize)
    REAL(TRC_KIND), INTENT(OUT) :: fy(1:isize, 0:jsize, 1:ksize)
    REAL(TRC_KIND), INTENT(OUT) :: fz(1:isize, 1:jsize, 0:ksize)

    INTEGER :: i, j, k
    INTEGER :: km1, km2

!$OMP PARALLEL DO PRIVATE(km1, km2)
    DO k=0, ksize
       km1 = maskindex(k)
       km2 = maskindex(k+1)

       DO j=1, jsize
       DO i=1, isize
          fz(i,j,k) = - imask3d(i,j,km1)*imask3d(i,j,km2) * (kappa(i,j,k)+kappa(i,j,k+1))*0.5 * (a(i,j,k+1)-a(i,j,k))*idz1(k)
       END DO
       END DO

       IF (k==0) CYCLE

       DO j=1, jsize
       DO i=0, isize
          fx(i,j,k) = - imask3d(i,j,km1)*imask3d(i+1,j,km1) * (kappa(i,j,k)+kappa(i+1,j,k))*0.5 * (a(i+1,j,k)-a(i,j,k))*idx1(i,j)
       END DO
       END DO

       DO j=0, jsize
       DO i=1, isize
          fy(i,j,k) = - imask3d(i,j,km1)*imask3d(i,j+1,km1) * (kappa(i,j,k)+kappa(i,j+1,k))*0.5 * (a(i,j+1,k)-a(i,j,k))*idy1(i,j)
       END DO
       END DO
    END DO

  END SUBROUTINE diffusion_iso

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE diffusion_h(a, kappa, fx, fy)
    REAL(TRC_KIND), INTENT(IN)  :: a(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)

    REAL(TRC_KIND), INTENT(IN)  :: kappa(0:isize+1, 0:jsize+1, 0:ksize+1)

    REAL(TRC_KIND), INTENT(OUT) :: fx(0:isize, 1:jsize, 1:ksize)
    REAL(TRC_KIND), INTENT(OUT) :: fy(1:isize, 0:jsize, 1:ksize)

    INTEGER :: i, j, k
    INTEGER :: km

!$OMP PARALLEL DO PRIVATE(km)
    DO k=1, ksize
       km = maskindex(k)

       DO j=1, jsize
       DO i=0, isize
          fx(i,j,k) = - imask3d(i,j,km)*imask3d(i+1,j,km) * (kappa(i,j,k)+kappa(i+1,j,k))*0.5 * (a(i+1,j,k)-a(i,j,k))*idx1(i,j)
       END DO
       END DO

       DO j=0, jsize
       DO i=1, isize
          fy(i,j,k) = - imask3d(i,j,km)*imask3d(i,j+1,km) * (kappa(i,j,k)+kappa(i,j+1,k))*0.5 * (a(i,j+1,k)-a(i,j,k))*idy1(i,j)
       END DO
       END DO
    END DO

  END SUBROUTINE diffusion_h

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE diffusion_v(a, kappa, fz)
    REAL(TRC_KIND), INTENT(IN)  :: a(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)

    REAL(TRC_KIND), INTENT(IN)  :: kappa(0:isize+1, 0:jsize+1, 0:ksize+1)

    REAL(TRC_KIND), INTENT(OUT) :: fz(1:isize, 1:jsize, 0:ksize)

    INTEGER :: i, j, k
    INTEGER :: km1, km2

!$OMP PARALLEL DO PRIVATE(km1, km2)
    DO k=0, ksize
       km1 = maskindex(k)
       km2 = maskindex(k+1)

       DO j=1, jsize
       DO i=1, isize
          fz(i,j,k) = - imask3d(i,j,km1)*imask3d(i,j,km2) * (kappa(i,j,k)+kappa(i,j,k+1))*0.5 * (a(i,j,k+1)-a(i,j,k))*idz1(k)
       END DO
       END DO
    END DO

  END SUBROUTINE diffusion_v

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE diffusion_idp(a, kappa1, kappa2, n_x, n_y, n_z, fx, fy, fz)
    REAL(TRC_KIND), INTENT(IN)  :: a(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)

    REAL(TRC_KIND), INTENT(IN)  :: kappa1(0:isize+1, 0:jsize+1, 0:ksize+1)
    REAL(TRC_KIND), INTENT(IN)  :: kappa2(0:isize+1, 0:jsize+1, 0:ksize+1)

    REAL(TRC_KIND), INTENT(IN)  :: n_x( -slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
    REAL(TRC_KIND), INTENT(IN)  :: n_y(1-slv:isize+slv,  -slv:jsize+slv, 1-slv:ksize+slv)
    REAL(TRC_KIND), INTENT(IN)  :: n_z(1-slv:isize+slv, 1-slv:jsize+slv,  -slv:ksize+slv)

    REAL(TRC_KIND), INTENT(OUT) :: fx(0:isize, 1:jsize, 1:ksize)
    REAL(TRC_KIND), INTENT(OUT) :: fy(1:isize, 0:jsize, 1:ksize)
    REAL(TRC_KIND), INTENT(OUT) :: fz(1:isize, 1:jsize, 0:ksize)

    REAL(TRC_KIND) :: dnx(0:isize+1, 0:jsize+1, 1:2)
    REAL(TRC_KIND) :: dny(0:isize+1, 0:jsize+1, 1:2)
    REAL(TRC_KIND) :: dnz(0:isize+1, 0:jsize+1, 1:2)

    INTEGER :: i, j, k
    INTEGER :: kstr, kend
    INTEGER :: ka, kb, kt
    INTEGER :: km1, km2, km3

    INTEGER :: tid

    tid  = 0
    kstr = 1
    kend = ksize

!$OMP PARALLEL PRIVATE(i, j, k, tid)  &
!$OMP          PRIVATE(dnx, dny, dnz) &
!$OMP          PRIVATE(kstr, kend) &
!$OMP          PRIVATE(ka, kb, kt) &
!$OMP          PRIVATE(km1, km2, km3)

    ka = 1
    kb = 2

!$  tid = omp_get_thread_num()
!$  kstr=kstr_t(tid)
!$  kend=kend_t(tid)

    DO k=kstr-2, kend
       kt=ka
       ka=kb
       kb=kt

       km1 = maskindex(k)
       km2 = maskindex(k+1)
       km3 = maskindex(k+2)

       DO j=0, jsize+1
       DO i=0, isize+1
          dnx(i,j,kb) = imask3d(i,j,km2)*imask3d(i+1,j,km2) * n_x(i,j,k+1) * (a(i+1,j,k+1)-a(i,j,k+1))*idx1(i,j)
          dny(i,j,kb) = imask3d(i,j,km2)*imask3d(i,j+1,km2) * n_y(i,j,k+1) * (a(i,j+1,k+1)-a(i,j,k+1))*idy1(i,j)
          dnz(i,j,kb) = imask3d(i,j,km1)*imask3d(i,j,  km2) * n_z(i,j,k)   * (a(i,j,  k+1)-a(i,j,k))*idz1(k)
       END DO
       END DO

       IF (k==kstr-2) CYCLE
!$     IF (tid/=0 .AND. k==kstr-1) CYCLE

       DO j=1, jsize
       DO i=1, isize
          fz(i,j,k) = - imask3d(i,j,km1)*imask3d(i,j,km2) &
                  * ( (kappa1(i,j,k)+kappa1(i,j,k+1))*0.5 * (a(i,j,k+1)-a(i,j,k))*idz1(k)            &
                    - (kappa1(i,j,k)+kappa1(i,j,k+1)-kappa2(i,j,k)-kappa2(i,j,k+1))*0.5 * n_z(i,j,k) &
                       * (dnz(i,j,kb) + (dnx(i-1,j,ka)+dnx(i,j,ka)+dnx(i-1,j,kb)+dnx(i,j,kb))*0.25   &
                                      + (dny(i,j-1,ka)+dny(i,j,ka)+dny(i,j-1,kb)+dny(i,j,kb))*0.25))
       END DO
       END DO

       IF (k==kstr-1) CYCLE

       DO j=1, jsize
       DO i=0, isize
          fx(i,j,k) = - imask3d(i,j,km1)*imask3d(i+1,j,km1) &
               * ( (kappa1(i,j,k)+kappa1(i+1,j,k))*0.5 * (a(i+1,j,k)-a(i,j,k))*idx1(i,j)             &
                 - (kappa1(i,j,k)+kappa1(i+1,j,k)-kappa2(i,j,k)-kappa2(i+1,j,k))*0.5 * n_x(i,j,k)    &
                    * (dnx(i,j,ka) + (dny(i,j-1,ka)+dny(i,j,ka)+dny(i+1,j-1,ka)+dny(i+1,j,ka))*0.25  &
                                   + (dnz(i,j,  ka)+dnz(i,j,kb)+dnz(i+1,j,  ka)+dnz(i+1,j,kb))*0.25))
          END DO
       END DO

       DO j=0, jsize
       DO i=1, isize
          fy(i,j,k) = - imask3d(i,j,km1)*imask3d(i,j+1,km1) &
               * ( (kappa1(i,j,k)+kappa1(i,j+1,k))*0.5 * (a(i,j+1,k)-a(i,j,k))*idy1(i,j)             &
                 - (kappa1(i,j,k)+kappa1(i,j+1,k)-kappa2(i,j,k)-kappa2(i,j+1,k))*0.5 * n_y(i,j,k)    &
                    * (dny(i,j,ka) + (dnz(i,  j,ka)+dnz(i,j,kb)+dnz(i,  j+1,ka)+dnz(i,j+1,kb))*0.25  &
                                   + (dnx(i-1,j,ka)+dnx(i,j,ka)+dnx(i-1,j+1,ka)+dnx(i,j+1,ka))*0.25))
       END DO
       END DO
    END DO
!$OMP END PARALLEL

  END SUBROUTINE diffusion_idp

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE diffusion_nottuned(a, kappa1, kappa2, n_x, n_y, n_z, fx, fy, fz)
    REAL(TRC_KIND), INTENT(IN)  :: a(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)

    REAL(TRC_KIND), INTENT(IN)  :: kappa1(0:isize+1, 0:jsize+1, 0:ksize+1)
    REAL(TRC_KIND), INTENT(IN)  :: kappa2(0:isize+1, 0:jsize+1, 0:ksize+1)

    REAL(TRC_KIND), INTENT(IN)  :: n_x( -slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
    REAL(TRC_KIND), INTENT(IN)  :: n_y(1-slv:isize+slv,  -slv:jsize+slv, 1-slv:ksize+slv)
    REAL(TRC_KIND), INTENT(IN)  :: n_z(1-slv:isize+slv, 1-slv:jsize+slv,  -slv:ksize+slv)

    REAL(TRC_KIND), INTENT(OUT) :: fx(0:isize, 1:jsize, 1:ksize)
    REAL(TRC_KIND), INTENT(OUT) :: fy(1:isize, 0:jsize, 1:ksize)
    REAL(TRC_KIND), INTENT(OUT) :: fz(1:isize, 1:jsize, 0:ksize)

    REAL(TRC_KIND) :: dnx(0:isize+1, 0:jsize+1, 0:ksize+1)
    REAL(TRC_KIND) :: dny(0:isize+1, 0:jsize+1, 0:ksize+1)
    REAL(TRC_KIND) :: dnz(0:isize+1, 0:jsize+1, 0:ksize+1)

    INTEGER :: i, j, k

    DO k=0, ksize+1
    DO j=0, jsize+1
    DO i=0, isize+1
       dnx(i,j,k) = imask3d(i,j,k)*imask3d(i+1,j,k) * n_x(i,j,k) * (a(i+1,j,k)-a(i,j,k))*idx1(i,j)
       dny(i,j,k) = imask3d(i,j,k)*imask3d(i,j+1,k) * n_y(i,j,k) * (a(i,j+1,k)-a(i,j,k))*idy1(i,j)
       dnz(i,j,k) = imask3d(i,j,k)*imask3d(i,j,k+1) * n_z(i,j,k) * (a(i,j,k+1)-a(i,j,k))*idz1(k)
    END DO
    END DO
    END DO

    DO k=1, ksize
    DO j=1, jsize
    DO i=0, isize
       fx(i,j,k) = - imask3d(i,j,k)*imask3d(i+1,j,k) &
            * ( (kappa1(i,j,k)+kappa1(i+1,j,k))*0.5 * (a(i+1,j,k)-a(i,j,k))*idx1(i,j)        &
              + (kappa2(i,j,k)+kappa2(i+1,j,k))*0.5 * n_x(i,j,k)                             &
                 * (dnx(i,j,k) + (dny(i,j-1,k)+dny(i,j,k)+dny(i+1,j-1,k)+dny(i+1,j,k))*0.25  &
                               + (dnz(i,j,k-1)+dnz(i,j,k)+dnz(i+1,j,k-1)+dnz(i+1,j,k))*0.25))
    END DO
    END DO
    END DO

    DO k=1, ksize
    DO j=0, jsize
    DO i=1, isize
       fy(i,j,k) = - imask3d(i,j,k)*imask3d(i,j+1,k) &
            * ( (kappa1(i,j,k)+kappa1(i,j+1,k))*0.5 * (a(i,j+1,k)-a(i,j,k))*idy1(i,j)        &
              + (kappa2(i,j,k)+kappa2(i,j+1,k))*0.5 * n_y(i,j,k)                             &
                 * (dny(i,j,k) + (dnz(i,j,k-1)+dnz(i,j,k)+dnz(i,j+1,k-1)+dnz(i,j+1,k))*0.25  &
                               + (dnx(i-1,j,k)+dnx(i,j,k)+dnx(i-1,j+1,k)+dnx(i,j+1,k))*0.25))
    END DO
    END DO
    END DO

    DO k=0, ksize
    DO j=1, jsize
    DO i=1, isize
       fz(i,j,k) = - imask3d(i,j,k)*imask3d(i,j,k+1) &
            * ( (kappa1(i,j,k)+kappa1(i,j,k+1))*0.5 * (a(i,j,k+1)-a(i,j,k))*idz1(k)          &
              + (kappa2(i,j,k)+kappa2(i,j,k+1))*0.5 * n_z(i,j,k)                             &
                 * (dnz(i,j,k) + (dnx(i-1,j,k)+dnx(i,j,k)+dnx(i-1,j,k+1)+dnx(i,j,k+1))*0.25  &
                               + (dny(i,j-1,k)+dny(i,j,k)+dny(i,j-1,k+1)+dny(i,j,k+1))*0.25))
    END DO
    END DO
    END DO

  END SUBROUTINE diffusion_nottuned

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE diffusion_vimp(a, kappa, fz, factor)
    REAL(TRC_KIND), INTENT(IN)  :: a(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
    REAL(TRC_KIND), INTENT(IN)  :: kappa(0:isize+1, 0:jsize+1, 0:ksize+1)
    REAL(TRC_KIND), INTENT(OUT) :: fz(1:isize, 1:jsize, 0:ksize)
    REAL(TRC_KIND), INTENT(IN), OPTIONAL  :: factor

    REAL(TRC_KIND) :: x(isize,jsize,0:ksize+1)
    REAL(TRC_KIND) :: c(isize,jsize,0:ksize)
    REAL(TRC_KIND) :: w(isize,jsize,ksize+1)
    REAL(TRC_KIND) :: t(isize,jsize)

    INTEGER :: i, j, k
    INTEGER :: km0, km1

    REAL(TRC_KIND) :: f

    CALL assert(.NOT. cycle_z, "implicit vertical diffusivity is not supported for CYCLE_Z")

    f = 1.0
    IF (PRESENT(factor)) f = factor

!$OMP PARALLEL DO
    DO k=0, ksize
    DO j=1, jsize
    DO i=1, isize
       c(i,j,k) = 0.5*f*(kappa(i,j,k)+kappa(i,j,k+1))*idz1(k)*imask3d(i,j,k)*imask3d(i,j,k+1)*dtime
    END DO
    END DO
    END DO

!$OMP PARALLEL DO
    DO j=1, jsize
    DO i=1, isize
       t(i,j)   = 1.0
       x(i,j,0) = 0.0
    END DO
    END DO

    CALL lrecv(t,        tag=1)
    CALL lrecv(x(:,:,0), tag=2)

    DO k=1, ksize
!$OMP PARALLEL DO
       DO j=1, jsize
       DO i=1, isize
          w(i,j,k) = -c(i,j,k-1) / t(i,j)
          t(i,j)   = dz_star(i,j,k) + c(i,j,k-1) + c(i,j,k) + c(i,j,k-1)*w(i,j,k)
          x(i,j,k) = (dz_star(i,j,k)*a(i,j,k) + c(i,j,k-1)*x(i,j,k-1)) / t(i,j)
       END DO
       END DO
    END DO

!$OMP PARALLEL DO PRIVATE(km0,km1)
    DO j=1, jsize
    DO i=1, isize
       w(i,j,ksize+1) = -c(i,j,ksize) / t(i,j)
    END DO
    END DO

    CALL usend(t,              tag=1)
    CALL usend(x(:,:,ksize),   tag=2)

    CALL urecv(x(:,:,ksize+1), tag=3)

    DO k=ksize, 0, -1
       km0 = maskindex(k)
       km1 = maskindex(k+1)

       IF (k==ksize .AND. vrank==0) THEN
          fz(:,:,k) = 0.0
          CYCLE
       END IF
!$OMP PARALLEL DO COLLAPSE(2)
       DO j=1, jsize
       DO i=1, isize
          x(i,j,k)  = x(i,j,k) - w(i,j,k+1)*x(i,j,k+1)
          fz(i,j,k) = - (x(i,j,k+1) - x(i,j,k))*c(i,j,k) * imask3d(i,j,km0)*imask3d(i,j,km1) * idtime
       END DO
       END DO
    END DO

    CALL lsend(x(:,:,1), tag=3)

  END SUBROUTINE diffusion_vimp

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE advection_open(name, a, u, v, w, fx, fy, fz)
    CHARACTER(*), INTENT(IN) :: name

    REAL(TRC_KIND), INTENT(IN) :: a(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
    REAL(TRC_KIND), INTENT(IN) :: u( -slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
    REAL(TRC_KIND), INTENT(IN) :: v(1-slv:isize+slv,  -slv:jsize+slv, 1-slv:ksize+slv)
    REAL(TRC_KIND), INTENT(IN) :: w(1-slv:isize+slv, 1-slv:jsize+slv,  -slv:ksize+slv)

    REAL(TRC_KIND), INTENT(INOUT) :: fx(0:isize, 1:jsize, 1:ksize)
    REAL(TRC_KIND), INTENT(INOUT) :: fy(1:isize, 0:jsize, 1:ksize)
    REAL(TRC_KIND), INTENT(INOUT) :: fz(1:isize, 1:jsize, 0:ksize)

    REAL(TRC_KIND) :: tmp_xy(isize,jsize)
    REAL(TRC_KIND) :: tmp_xz(isize,ksize)
    REAL(TRC_KIND) :: tmp_yz(jsize,ksize)

    INTEGER :: i, j, k

    IF (open_w) THEN

       IF (icoord==0) THEN
!$OMP PARALLEL DO
          DO k=1, ksize
          DO j=1, jsize
             fx(0,j,k) = flux_upwind(u(0,j,k), a(0,j,k), a(1,j,k))
          END DO
          END DO
       END IF

       IF (require_checkout('TIN_'//trim(name)//'_WEST')) THEN
          tmp_yz(:,:) = 0.0

          IF (icoord==0) THEN
!$OMP PARALLEL DO
             DO k=1, ksize
             DO j=1, jsize
                IF (u(0,j,k) > 0.0) tmp_yz(j,k) = fx(0,j,k)*dsx_old(0,j,k)
             END DO
             END DO
          END IF

          CALL checkout('TIN_'//trim(name)//'_WEST', tmp_yz, section='YZ') !tracer inflow transport (volume integral, per second)
       END IF

       IF (require_checkout('TOUT_'//trim(name)//'_WEST')) THEN
          tmp_yz(:,:) = 0.0

          IF (icoord==0) THEN
!$OMP PARALLEL DO
             DO k=1, ksize
             DO j=1, jsize
                IF (u(0,j,k) < 0.0) tmp_yz(j,k) = -fx(0,j,k)*dsx_old(0,j,k)
             END DO
             END DO
          END IF

          CALL checkout('TOUT_'//trim(name)//'_WEST', tmp_yz, section='YZ') !tracer outflow transport (volume integral, per second)
       END IF
    END IF

    IF (open_e) THEN
       IF (icoord==ipes-1) THEN
!$OMP PARALLEL DO
          DO k=1, ksize
          DO j=1, jsize
             fx(isize,j,k) = flux_upwind(u(isize,j,k), a(isize,j,k), a(isize+1,j,k))
          END DO
          END DO
       END IF

       IF (require_checkout('TIN_'//trim(name)//'_EAST')) THEN
          tmp_yz(:,:) = 0.0

          IF (icoord==ipes-1) THEN
!$OMP PARALLEL DO
             DO k=1, ksize
             DO j=1, jsize
                IF (u(isize,j,k) < 0.0) tmp_yz(j,k) = -fx(isize,j,k)*dsx_old(isize,j,k)
             END DO
             END DO
          END IF

          CALL checkout('TIN_'//trim(name)//'_EAST', tmp_yz, section='YZ') !tracer inflow transport (volume integral, per second)
       END IF

       IF (require_checkout('TOUT_'//trim(name)//'_EAST')) THEN
          tmp_yz(:,:) = 0.0

          IF (icoord==ipes-1) THEN
!$OMP PARALLEL DO
             DO k=1, ksize
             DO j=1, jsize
                IF (u(isize,j,k) > 0.0) tmp_yz(j,k) = fx(isize,j,k)*dsx_old(isize,j,k)
             END DO
             END DO
          END IF

          CALL checkout('TOUT_'//trim(name)//'_EAST', tmp_yz, section='YZ') !tracer outflow transport (volume integral, per second)
       END IF
    END IF

    IF (open_s) THEN
       IF (jcoord==0) THEN
!$OMP PARALLEL DO
          DO k=1, ksize
          DO i=1, isize
             fy(i,0,k) = flux_upwind(v(i,0,k), a(i,0,k), a(i,1,k))
          END DO
          END DO
       END IF

       IF (require_checkout('TIN_'//trim(name)//'_SOUTH')) THEN
          tmp_xz(:,:) = 0.0

          IF (jcoord==0) THEN
!$OMP PARALLEL DO
             DO k=1, ksize
             DO i=1, isize
                IF (v(i,0,k) > 0.0) tmp_xz(i,k) = fy(i,0,k)*dsy_old(i,0,k)
             END DO
             END DO
          END IF

          CALL checkout('TIN_'//trim(name)//'_SOUTH', tmp_xz, section='XZ') !tracer inflow transport (volume integral, per second)
       END IF

       IF (require_checkout('TOUT_'//trim(name)//'_SOUTH')) THEN
          tmp_xz(:,:) = 0.0

          IF (jcoord==0) THEN
!$OMP PARALLEL DO
             DO k=1, ksize
             DO i=1, isize
                IF (v(i,0,k) < 0.0) tmp_xz(i,k) = -fy(i,0,k)*dsy_old(i,0,k)
             END DO
             END DO
          END IF

          CALL checkout('TOUT_'//trim(name)//'_SOUTH', tmp_xz, section='XZ') !tracer outflow transport (volume integral, per second)
       END IF
    END IF

    IF (open_n) THEN
       IF (jcoord==jpes-1) THEN
!$OMP PARALLEL DO
          DO k=1, ksize
          DO i=1, isize
             fy(i,jsize,k) = flux_upwind(v(i,jsize,k), a(i,jsize,k), a(i,jsize+1,k))
          END DO
          END DO
       END IF

       IF (require_checkout('TIN_'//trim(name)//'_NORTH')) THEN
          tmp_xz(:,:) = 0.0

          IF (jcoord==jpes-1) THEN
!$OMP PARALLEL DO
             DO k=1, ksize
             DO i=1, isize
                IF (v(i,jsize,k) < 0.0) tmp_xz(i,k) = -fy(i,jsize,k)*dsy_old(i,jsize,k)
             END DO
             END DO
          END IF

          CALL checkout('TIN_'//trim(name)//'_NORTH', tmp_xz, section='XZ') !tracer inflow transport (volume integral, per second)
       END IF

       IF (require_checkout('TOUT_'//trim(name)//'_NORTH')) THEN
          tmp_xz(:,:) = 0.0

          IF (jcoord==jpes-1) THEN
!$OMP PARALLEL DO
             DO k=1, ksize
             DO i=1, isize
                IF (v(i,jsize,k) > 0.0) tmp_xz(i,k) = fy(i,jsize,k)*dsy_old(i,jsize,k)
             END DO
             END DO
          END IF

          CALL checkout('TOUT_'//trim(name)//'_NORTH', tmp_xz, section='XZ') !tracer outflow transport (volume integral, per second)
       END IF
    END IF

    IF (input_registered('W_BOTTOM')) THEN
       tmp_xy(:,:) = UNDEF

       CALL checkin(trim(name)//'_BOTTOM', tmp_xy)

       DO j=1, jsize
       DO i=1, isize
          IF (.NOT. bottom_flag(i,j)) CYCLE
          k = bottom_k(i,j)

          IF (tmp_xy(i,j) == UNDEF .OR. w(i,j,k-1) < 0.0) THEN
             fz(i,j,k-1) = w(i,j,k-1)*a(i,j,k)
          ELSE
             fz(i,j,k-1) = w(i,j,k-1)*tmp_xy(i,j)
          END IF
       END DO
       END DO

    END IF

  END SUBROUTINE advection_open

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE report_tracers(force)
    USE calendar
    LOGICAL, INTENT(IN), OPTIONAL :: force
    REAL(8) :: max, min, total, tmp
    INTEGER :: i, j, k, n
    LOGICAL :: flag

    flag = .FALSE.

    IF (tracer_report_interval > 0) flag = mod(n_timestep, tracer_report_interval) == 0
    IF (present(force)) flag = flag .OR. force

    IF (.NOT. flag) RETURN

    DO n=1, n_tracer
       IF (.NOT. present(force)) THEN
          IF (tracer_report_interval <= 0) CYCLE
          IF (mod(n_timestep, tracer_report_interval)/=0) CYCLE
       END IF

       max = maxval(tracer(1:isize,1:jsize,1:ksize,n), mask=lmask3d(1:isize,1:jsize,1:ksize))
       min = minval(tracer(1:isize,1:jsize,1:ksize,n), mask=lmask3d(1:isize,1:jsize,1:ksize))
       total = sum(tracer(1:isize,1:jsize,1:ksize,n)*imask3d(1:isize,1:jsize,1:ksize)*dvol(1:isize,1:jsize,1:ksize))

       CALL gmax(max)
       CALL gmin(min)
       CALL gsum(total)

       IF (rank==0) THEN
          WRITE(REPORT_UNIT, '(A, ": ",  A8, " max=", ES10.3, ", min=", ES10.3, ", ave=", ES10.3, ", total=", ES13.6)') &
               trim(current_datetime), trim(tracer_info(n)%name), max, min, total/total_vol, total
       END IF
    END DO


  END SUBROUTINE report_tracers


END MODULE tracers
