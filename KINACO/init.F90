#include "macro.h"

SUBROUTINE init
  USE misc
  USE parameters
  USE geometry
  USE velocity
  USE subgrid
  USE tracers
  USE state
  USE calendar
  USE io
  USE particles
  USE seaice
  USE npzd
  USE ecosystem
  USE gls
#if defined(DLES_VISCOSITY) || defined(DLES_DIFFUSIVITY)
  USE dles
#endif
#ifdef PROFILE
  USE profile
#endif
  IMPLICIT NONE

  CHARACTER(16) :: start, start_date
  CHARACTER(16) :: end,   end_date

  CHARACTER(4096) :: line
  INTEGER :: iosize

  CHARACTER(512) :: reportfile = ''
  CHARACTER(512) :: inputdir   = '.'
  CHARACTER(512) :: outputdir  = '.'
  CHARACTER(512) :: initialdir  = ''

  CHARACTER(512) :: tmp = ""

  LOGICAL :: reportfile_allproc  = .FALSE.
  LOGICAL :: reportfile_append   = .TRUE.

  INTEGER :: volume_report_interval = 0 ! for backward-comaptibility, use ssh_report_interval

  INTEGER        :: iostat
  CHARACTER(256) :: iomsg

  INTEGER :: i, j, k, n

  NAMELIST / runinfo / &
       runname,       &
       hydrostatic,   &
       rigid_lid,     &
       fix_meanssh,   &
       iceshelf_melt, &
       vcoord_zstar,  &
       start_date,    &
       start,         &
       end_date,      &
       end,           &
       dtime,         &
       calendar_type,   &
       calendar_format, &
       perfect_restart,&
       inputdir,       &
       outputdir,      &
       initialdir,     &
       cachedir,       &
       workdir,        &
       reportfile,     &
       reportfile_allproc,     &
       reportfile_append,      &
       cfl_report_interval,    &
       tracer_report_interval, &
       ssh_report_interval,    &
       volume_report_interval, &
       open_report_interval,   &
       momentum_report_interval, &
       energy_report_interval, &
       stress_report_interval, &
       les_report_interval,    &
       offline, &
       use_landwater,  &
       use_gls, &
       use_particles,  &
       use_ecosystem,  &
       use_npzd,       &
       report_qsum, &
       assert_nan,  &
       bypass_solver

  CALL init_parallel

#ifdef DEBUG
  report_qsum = .TRUE.
  assert_nan  = .TRUE.

  cfl_report_interval    = 1
  tracer_report_interval = 1
  volume_report_interval = 1
  stress_report_interval = 1
  les_report_interval    = 1
#endif

  IF (rank==0) THEN
     calendar_type   = 0
     calendar_format = 0
     start      = '0000000_000000'
     end        = ''
     start_date = ''
     end_date   = ''

     OPEN(CONFIG_UNIT, &
          STATUS = 'SCRATCH',    &
          ACTION = 'READWRITE',  &
          ACCESS = 'SEQUENTIAL', &
          IOSTAT = iostat)
     CALL assert(iostat==0, "failed to create a temporary file for the configuration namelist")

     DO WHILE(iostat/=-1)
        line=''
        READ(*, '(A)', IOSTAT=iostat, SIZE=iosize, ADVANCE='NO') line
        CALL assert(iosize/=len(line), "line length in the configuration file is too long")
        IF (index(adjustl(line),'#') == 1) CYCLE
        WRITE(CONFIG_UNIT, '(A)') trim(line)
     END DO

     REWIND(CONFIG_UNIT)
     READ(CONFIG_UNIT, NML=runinfo, IOSTAT=iostat, IOMSG=iomsg)
     CALL assert(iostat==0, "failed to read RUNINFO namelist", iomsg)

     IF (start_date /='') start = start_date
     IF (end_date   /='') end   = end_date

     !for backward-compatibility
     IF (ssh_report_interval==0 .AND. volume_report_interval/=0) ssh_report_interval = volume_report_interval

#ifdef F2003
     IF (workdir(1:1)=="$") THEN
        CALL get_environment_variable(name=workdir(2:), value=tmp, status=iostat)
        CALL assert(iostat==0, "failed to get environment variable "//trim(workdir)//" for WORKDIR")
        workdir = tmp
     END IF

     IF (cachedir(1:1)=="$") THEN
        CALL get_environment_variable(name=cachedir(2:), value=tmp, status=iostat)
        CALL assert(iostat==0, "failed to get environment variable "//trim(cachedir)//" for CACHEDIR")
        cachedir = tmp
     END IF
#endif
  END IF

  CALL bcast(runname)

  CALL bcast(hydrostatic)
  CALL bcast(rigid_lid)
  CALL bcast(fix_meanssh)
  CALL bcast(iceshelf_melt)
  CALL bcast(vcoord_zstar)

  CALL bcast(offline)

  CALL bcast(dtime)
  CALL bcast(start)
  CALL bcast(end)

  CALL assert(dtime > 0.0, "DTIME is mandatory in RUNINFO")
  CALL assert(end/='',     "END is mandatory in RUNINFO")

  CALL bcast(perfect_restart)

  CALL bcast(calendar_type)
  CALL bcast(calendar_format)

  CALL assert(calendar_type   >= 0 .AND. calendar_type   <= 2, "unsupported CALENDAR_TYPE")
  CALL assert(calendar_format >= 0 .AND. calendar_format <= 2, "unsupported CALENDAR_FORMAT")

  t_start = datetime_seconds(start)
  IF (end(1:1)=='+') THEN
     t_end = t_start + interval_seconds(end(2:)//' ')
  ELSE
     t_end = datetime_seconds(end)
  END IF

  idtime = 1.0 / dtime
  idtime_r4 = REAL(idtime,4)

  t_current = t_start

  current_datetime = format_datetime(t_current)
  start_datetime   = format_datetime(t_start)
  end_datetime     = format_datetime(t_end)

  CALL bcast(inputdir)
  CALL bcast(outputdir)
  CALL bcast(initialdir)

  CALL bcast(workdir)
  CALL bcast(cachedir)

  CALL bcast(reportfile)
  CALL bcast(reportfile_allproc)
  CALL bcast(reportfile_append)

  CALL replace_vars(inputdir)
  CALL replace_vars(outputdir)
  CALL replace_vars(initialdir)
  CALL replace_vars(reportfile)

  CALL bcast(use_landwater)
  CALL bcast(use_gls)
  CALL bcast(use_particles)
  CALL bcast(use_ecosystem)
  CALL bcast(use_npzd)

  CALL init_reportfile

  IF (rank==0) THEN
     CALL timestamp(REPORT_UNIT)

     WRITE(REPORT_UNIT, *) "-- initialization start --"

     IF (check_endian()) THEN
        WRITE(REPORT_UNIT, *) "Running on a Big-Endian system"
     ELSE
        WRITE(REPORT_UNIT, *) "Running on a Little-Endian system"
     END IF

     WRITE(REPORT_UNIT, *) "RUNNAME='"//trim(runname)//"'"

     WRITE(REPORT_UNIT, *) "START='"//trim(start)//"' , END='"//trim(end)//"', dtime=", dtime, "[s]"
     IF (offline) WRITE(REPORT_UNIT, *) "OFFLINE mode"
  END IF

  CALL bcast(     cfl_report_interval)
  CALL bcast(  tracer_report_interval)
  CALL bcast(     ssh_report_interval)
  CALL bcast(    open_report_interval)
  CALL bcast(  energy_report_interval)
  CALL bcast(momentum_report_interval)
  CALL bcast(  stress_report_interval)
  CALL bcast(     les_report_interval)

  CALL bcast(report_qsum)
  CALL bcast(assert_nan)

  CALL bcast(bypass_solver)
  CALL assert(.NOT. rigid_lid .OR. .NOT. use_landwater, "RIGID_LID and USE_LANDWATER cannot be enabled simultaneously")

  CALL init_parameters

  CALL init_geometry(inputdir, outputdir)

#ifdef PARALLEL_MPI
  IF (comm == MPI_COMM_NULL) RETURN
#endif

  IF (.NOT. bypass_solver .AND. .NOT. offline) CALL init_solver

  CALL init_io(inputdir, outputdir, initialdir)

  IF (.NOT. rigid_lid) THEN
     CALL initial_data('SSH', ssh, default=.TRUE.)
     ssh_old(:,:) = ssh(:,:)

     CALL update_geometry
  END IF

  CALL init_velocity
  CALL init_subgrid

  CALL init_tracers

  CALL init_state

  CALL init_gls

  CALL init_npzd
  CALL init_ecosystem

  CALL assert(.NOT. (use_npzd .AND. use_ecosystem), "NPZD and ECOSYSTEM cannot be enabled simultaneously")

#if defined(DLES_VISCOSITY) || defined(DLES_DIFFUSIVITY)
  CALL init_dles
#endif

  CALL init_particles

  CALL init_seaice
#ifdef PARALLEL_MPI
  cputime1 = mpi_wtime()
#else
  CALL cpu_time(cputime1)
#endif

  CALL barrier

  IF (rank==0) THEN
     str_format = format(n_tracer)
     WRITE(REPORT_UNIT, '(A)', ADVANCE='NO') "using " // trim(str_format) // " tracers: "
     DO n=1, n_tracer
        WRITE(REPORT_UNIT, '(A)', ADVANCE='NO') trim(tracer_name(n))
        IF (n/=n_tracer)  WRITE(REPORT_UNIT, '(A)', ADVANCE='NO') ", "
     END DO
     IF (TRC_KIND==4) WRITE(REPORT_UNIT, '(A)') " (single-precision tracer computation)"
     WRITE(REPORT_UNIT, *)
  END IF

  CALL report_ssh(force=.TRUE.)
  CALL report_tracers(force=.TRUE.)
  CALL report_momentum(force=.TRUE.)
  CALL report_energy(force=.TRUE.)
  CALL report_cfl(force=.TRUE.)

  IF (rank==0) WRITE(REPORT_UNIT, *) "-- initialization done! --"

CONTAINS

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE init_reportfile
  INTEGER       :: iostat
  CHARACTER(256):: file


  IF (rank==0 .OR. reportfile_allproc) THEN
     IF (reportfile /= '') THEN
        file = trim(path(outputdir, reportfile))
        IF (reportfile_allproc) THEN
           file = trim(file) // '.' // strrank
        END IF

        IF (reportfile_append) THEN
           OPEN(REPORT_UNIT, FILE=trim(file), STATUS='UNKNOWN', ACTION='WRITE', POSITION='APPEND', IOSTAT=iostat)
        ELSE
           OPEN(REPORT_UNIT, FILE=trim(file), STATUS='REPLACE', ACTION='WRITE',                    IOSTAT=iostat)
        END IF
        CALL assert(iostat==0, "fail to create reportfile '"//trim(file)//"'")
     ELSE
        CALL ASSERT(.NOT. reportfile_allproc, "REPORTFILE_ALLPROC requires to specify REPORTFILE")
        REPORT_UNIT = STDOUT_UNIT
     END IF
  ELSE
     REPORT_UNIT = -1
  END IF

END SUBROUTINE init_reportfile

SUBROUTINE init_solver
  USE solver2d
  USE solver3d

  CHARACTER(8) :: method
  INTEGER :: cmode ! 0 for V-cycle, 1 for W-cycle, 2 for F-cycle
  REAL(8) :: eps
  INTEGER :: itmax
  LOGICAL :: report

  INTEGER        :: iostat
  CHARACTER(256) :: iomsg

  INTEGER :: i, j, k

  REAL(8), ALLOCATABLE :: tmp(:,:,:)

  NAMELIST / solver /   &
       method, cmode, eps, itmax, report

  IF (rank==0) THEN
     method = 'MGCG'
     cmode  = 0
     eps    = 1.0D-8
     itmax  = 100
     report = .FALSE.
     REWIND(CONFIG_UNIT)
     READ(CONFIG_UNIT, NML=solver, IOSTAT=iostat, IOMSG=iomsg)
     CALL assert(iostat <= 0, "failed to read SOLVER namelist", iomsg)
  END IF

  CALL bcast(method)
  CALL bcast(cmode)
  CALL bcast(eps)
  CALL bcast(itmax)
  CALL bcast(report)

  IF (hydrostatic) THEN
     ALLOCATE(da_solver(0:isize+1, 0:jsize+1, 1:1))

     IF (rigid_lid) THEN
        da_solver(:,:,:) = 0.0D0
     ELSE
        DO j=0, jsize+1
        DO i=0, isize+1
           da_solver(i,j,1) = 2.0*dsz(i,j)/(gravity*dtime**2) * b_ebcn
        END DO
        END DO
     END IF

     IF (vrank==0)  THEN
        CALL init_solver2d(isize, jsize,                                             &
                           dx, dy, dsx2d, dsy2d, da_solver(:,:,1), lmask2d,          &
                           (/cycle_x, cycle_y/), (/open_e, open_w, open_n, open_s/), &
                           method, eps, cmode, itmax, report)
     END IF
  ELSE ! non-hydrostatic
     ALLOCATE(da_solver(0:isize+1, 0:jsize+1, 0:ksize+1))
     da_solver(:,:,:) = 0.0D0
     IF (.NOT. rigid_lid) THEN
        DO k=0, ksize+1
        DO j=0, jsize+1
        DO i=0, isize+1
           da_solver(i,j,k) = 2.0*dsz(i,j)/(gravity*dtime**2) * b_ebcn * imask3d(i,j,k)*cz_star(i,j,k)
        END DO
        END DO
        END DO
     END IF

     SELECT CASE (pgradient_scheme)
     CASE (1)
        CALL init_solver3d(isize, jsize, ksize,                      &
                           dx_h=dx, dy_h=dy, dz_v=dz,                &
                           dsx=dsx, dsy=dsy, dsz_h=dsz,              &
                           dvol=dvol, alpha=da_solver, mask=lmask3d, &
                           periods=(/cycle_x, cycle_y, cycle_z/),    &
                           opens=(/open_e, open_w, open_n, open_s/), &
                           method=method, epsilon=eps, c_mode=cmode, &
                           it_max=itmax, report=report, scheme=pgradient_scheme)
     CASE (2,3)
        CALL init_solver3d(isize, jsize, ksize,                      &
                           dx_h=dx, dy_h=dy, dz_v=dz, dz=dz_ref,     &
                           dsx=dsx, dsy=dsy, dsz_h=dsz,              &
                           dvol=dvol, alpha=da_solver, mask=lmask3d, &
                           periods=(/cycle_x, cycle_y, cycle_z/),    &
                           opens=(/open_e, open_w, open_n, open_s/), &
                           method=method, epsilon=eps, c_mode=cmode, &
                           it_max=itmax, report=report, scheme=pgradient_scheme)
     END SELECT

#ifdef __DEBUG__
     ALLOCATE(tmp(1:isize, 1:jsize, 1:ksize))
     CALL solver3d_coefficient(0, tmp)
     CALL write_data_3d(tmp, 'L3D')
     DEALLOCATE(tmp)
#endif
  END IF

END SUBROUTINE init_solver

SUBROUTINE timestamp(unit)
  INTEGER, INTENT(IN) :: unit

  CHARACTER(8)  :: date
  CHARACTER(10) :: time
  CHARACTER(5)  :: zone

  CALL date_and_time(date, time, zone)

  WRITE(unit, *) "******* " // date // "-" // time // zone // " *******"
END SUBROUTINE timestamp

END SUBROUTINE init

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE finalize
  USE misc
  USE geometry
  USE velocity
  USE subgrid
  USE tracers
  USE io
  USE solver2d
  USE solver3d
#if defined(DLES_VISCOSITY) || defined(DLES_DIFFUSIVITY)
  USE dles
#endif
  USE particles
#ifdef PROFILE
  USE profile
#endif
  IMPLICIT NONE

  INTEGER :: i, n

  IF (rank==0) THEN
#ifdef PARALLEL_MPI
     cputime2 = mpi_wtime()
#else
     CALL cpu_time(cputime2)
#endif
     WRITE(REPORT_UNIT, '("cpu time:",EN12.3,"[sec]")') cputime2 - cputime1

#ifdef PROFILE
     WRITE(REPORT_UNIT, '("main      :",EN12.3,"[sec]")') profile_time('main')
     WRITE(REPORT_UNIT, '("init      :",EN12.3,"[sec]")') profile_time('init')
     IF (hydrostatic) THEN
        WRITE(REPORT_UNIT, '("2D-solver  :",EN12.3,"[sec]")') profile_time('solver2d')
     ELSE
        WRITE(REPORT_UNIT, '("3D-solver  :",EN12.3,"[sec]")') profile_time('solver3d')
     END IF
     IF (.NOT. rigid_lid) THEN
        WRITE(REPORT_UNIT, '("solver-reset:",EN12.3,"[sec]")') profile_time('reset_solver')
     END IF
     WRITE(REPORT_UNIT, '("nonlinear  :",EN12.3,"[sec]")') profile_time('vadv')
     WRITE(REPORT_UNIT, '("viscosity  :",EN12.3,"[sec]")') profile_time('visc')
     WRITE(REPORT_UNIT, '("advection  :",EN12.3,"[sec]")') profile_time('tadv')
     WRITE(REPORT_UNIT, '("diffusion  :",EN12.3,"[sec]")') profile_time('diff')
     IF (use_particles) WRITE(REPORT_UNIT, '("particles  :",EN12.3,"[sec]")') profile_time('particle')
     WRITE(REPORT_UNIT, '("IO read    :",EN12.3,"[sec]")') profile_time('ioread')
     WRITE(REPORT_UNIT, '("IO write   :",EN12.3,"[sec]")') profile_time('iowrite')
     WRITE(REPORT_UNIT, *)
#endif
  END IF

  CALL finalize_solver2d

  IF (.NOT. hydrostatic) THEN
     CALL finalize_solver3d
  END IF

#if defined(DLES_VISCOSITY) || defined(DLES_DFFUSIVITY)
  CALL finalize_dles
#endif

  IF (use_particles) CALL finalize_particles

  CALL finalize_tracers
  CALL finalize_subgrid
  CALL finalize_velocity

  CALL finalize_io

  CALL finalize_geometry

  IF (REPORT_UNIT > 0 .AND. REPORT_UNIT /= STDOUT_UNIT) CLOSE(REPORT_UNIT)
  CLOSE (CONFIG_UNIT)

  CALL finalize_parallel

END SUBROUTINE finalize

!-----------------------------------------------------------------------------------------------------------------------
