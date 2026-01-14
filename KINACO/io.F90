#include "macro.h"

MODULE io
  USE misc
  USE geometry
  USE calendar
  IMPLICIT NONE
  PRIVATE
  PUBLIC init_io, finalize_io, flush_io
  PUBLIC checkin, checkout, require_checkin, require_checkout, checkout_zeros, check_time
  PUBLIC input_registered, output_registered
  PUBLIC default_input, default_output
  PUBLIC count_input, count_output
  PUBLIC initial_data, has_initial
  PUBLIC checkout_div
  PUBLIC global_initialdir, global_outputdir, global_inputdir
  PUBLIC workdir, cachedir
  PUBLIC replace_vars, replace_datetime
  PUBLIC dump_start, dump_end, dump_interval, dumpdir, dump_data

  CHARACTER(512) :: global_inputdir  = ''
  CHARACTER(512) :: global_outputdir = ''
  CHARACTER(512) :: global_initialdir = ''

  CHARACTER(16)  :: restart_start
  CHARACTER(16)  :: restart_end
  CHARACTER(16)  :: restart_interval
  CHARACTER(512) :: restart_outputdir

  INTERFACE checkin
     MODULE PROCEDURE checkin_1d
     MODULE PROCEDURE checkin_2d
     MODULE PROCEDURE checkin_3d
     MODULE PROCEDURE checkin_1d_r4
     MODULE PROCEDURE checkin_2d_r4
     MODULE PROCEDURE checkin_3d_r4
     MODULE PROCEDURE checkin_1d_i1
     MODULE PROCEDURE checkin_2d_i1
     MODULE PROCEDURE checkin_3d_i1
     MODULE PROCEDURE checkin_1d_logical
     MODULE PROCEDURE checkin_2d_logical
     MODULE PROCEDURE checkin_3d_logical
  END INTERFACE

  INTERFACE checkout
     MODULE PROCEDURE checkout_1d
     MODULE PROCEDURE checkout_2d
     MODULE PROCEDURE checkout_3d
     MODULE PROCEDURE checkout_1d_r4
     MODULE PROCEDURE checkout_2d_r4
     MODULE PROCEDURE checkout_3d_r4
     MODULE PROCEDURE checkout_1d_i1
     MODULE PROCEDURE checkout_2d_i1
     MODULE PROCEDURE checkout_3d_i1
  END INTERFACE

  INTERFACE checkin_report
     MODULE PROCEDURE checkin_report_1d
     MODULE PROCEDURE checkin_report_2d
     MODULE PROCEDURE checkin_report_3d
  END INTERFACE

  INTERFACE require_checkout
     MODULE PROCEDURE require_checkout_by_name
     MODULE PROCEDURE require_checkout_by_id
  END INTERFACE

  INTERFACE require_checkin
     MODULE PROCEDURE require_checkin_by_name
     MODULE PROCEDURE require_checkin_by_id
  END INTERFACE

  INTERFACE initial_data
     MODULE PROCEDURE initial_data_1d_r4
     MODULE PROCEDURE initial_data_1d_r8
     MODULE PROCEDURE initial_data_2d_r4
     MODULE PROCEDURE initial_data_2d_r8
     MODULE PROCEDURE initial_data_3d_r4
     MODULE PROCEDURE initial_data_3d_r8
  END INTERFACE initial_data

  INTERFACE dump_data
     MODULE PROCEDURE dump_data_1d
     MODULE PROCEDURE dump_data_2d
     MODULE PROCEDURE dump_data_3d
     MODULE PROCEDURE dump_data_4d
  END INTERFACE dump_data

  INTEGER, PARAMETER :: MAX_INPUT  = 1024
  INTEGER, PARAMETER :: MAX_OUTPUT = 1024

  TYPE :: output_registry_entry
#ifdef F2008
     REAL(8), POINTER, CONTIGUOUS :: buf(:,:,:)
#else
     REAL(8), POINTER :: buf(:,:,:)
#endif
     CHARACTER(32)  :: varname
     CHARACTER(128) :: basename
     CHARACTER(512) :: outputdir
     INTEGER :: kind
     INTEGER :: mode
     INTEGER :: region
     INTEGER :: slice
     INTEGER :: intspan(0:1)
     LOGICAL :: mean
     INTEGER :: dimcode
     REAL(8) :: start
     REAL(8) :: end
     REAL(8) :: interval
     REAL(8) :: lastwrite
     REAL(8) :: meanspan
     INTEGER :: meancount
     LOGICAL :: checked
     LOGICAL :: initialized
     LOGICAL :: omit_time
! dimcode
! 0: 3d
! 1: 2d (xy-plane)
! 2: xz plane
! 3: yz plane

! mode
! 0: full-3d
! 1: yz-slice at i=slice_x
! 2: zx-slice at j=slice_y
! 3: xy-slice at k=slice_z
! 4: surface
! 5: bottom
! 6: x-integrated
! 7: y-integrated
! 8: z-integrated
  END TYPE output_registry_entry

  TYPE :: input_registry_entry
#ifdef F2008
     REAL(8), POINTER, CONTIGUOUS :: buf(:,:,:,:)
#else
     REAL(8), POINTER :: buf(:,:,:,:)
#endif
     CHARACTER(128), POINTER :: filename(:)
     CHARACTER(32)  :: varname
     CHARACTER(512) :: inputdir
     CHARACTER(32)  :: suffix
     INTEGER :: kind
     INTEGER :: mode
     INTEGER :: periods
     INTEGER :: region
     INTEGER :: view
     LOGICAL :: descend
     REAL(4) :: scale
     REAL(4) :: offset
     REAL(4) :: minlimit
     REAL(4) :: maxlimit
     REAL(4) :: missing(2)
     REAL(4) :: noise
     REAL(8) :: start
     REAL(8) :: end
     REAL(8) :: interval
     REAL(8) :: tshift
     REAL(8) :: tcycle
     REAL(8) :: lastread
     REAL(8) :: offset_cx
     REAL(8) :: offset_cy
     REAL(8) :: offset_cz
     INTEGER :: current
     LOGICAL :: initialized
     LOGICAL :: report
     LOGICAL :: omit_time

! mode
! 0: constant
! 1: linear  interporation
! 2: sigmoid interporation
!10: cyclic with constant value
!11: cyclic with linear  interpolation
!12: cyclic with sigmoid interpolation
!20: historical with constant value
!21: historical with linear  interpolation
!22: historical with sigmoid interpolation
!30: sine-curve
!40: gaussian
  END TYPE input_registry_entry

  TYPE :: initial_params
     CHARACTER(1024) :: filepath
     INTEGER :: kind
     REAL(4) :: scale
     REAL(4) :: offset
     REAL(4) :: minlimit
     REAL(4) :: maxlimit
     REAL(4) :: missing(2)
     REAL(4) :: noise
     INTEGER :: fileview
     LOGICAL :: descend
  END TYPE initial_params

  TYPE(input_registry_entry),  SAVE :: input_registry(MAX_INPUT)
  TYPE(output_registry_entry), SAVE :: output_registry(MAX_OUTPUT)

  INTEGER, SAVE :: input_registry_index( 0:128)
  INTEGER, SAVE :: output_registry_index(0:128)

  INTEGER, SAVE :: n_output = 0
  INTEGER, SAVE :: n_input = 0

  CHARACTER(512) :: workdir  = ""
  CHARACTER(512) :: cachedir = ""

  REAL(8) :: dump_start = UNDEF
  REAL(8) :: dump_end   = UNDEF
  REAL(8) :: dump_interval = 0.0
  CHARACTER(256) :: dumpdir = ""
  REAL(8) :: lastdump = UNDEF

CONTAINS
  SUBROUTINE init_io(inputdir, outputdir, initialdir)
    CHARACTER(512), INTENT(IN) :: inputdir
    CHARACTER(512), INTENT(IN) :: outputdir
    CHARACTER(512), INTENT(IN) :: initialdir

    global_inputdir   = inputdir
    global_outputdir  = outputdir
    IF (trim(initialdir) == "") THEN
       global_initialdir = inputdir
     ELSE
       global_initialdir = initialdir
    END IF

    CALL read_input_namelist
    CALL read_output_namelist

    CALL read_dump_namelist

  END SUBROUTINE init_io

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE read_input_namelist
    INTEGER, PARAMETER :: nmax = 4096

    CHARACTER(32)  :: varname
    CHARACTER(128) :: filename(0:nmax-1)
    CHARACTER(8)   :: precision
    LOGICAL        :: descend
    CHARACTER(20)  :: mode
    INTEGER        :: periods
    REAL(4)        :: scale
    REAL(4)        :: offset
    REAL(4)        :: minlimit
    REAL(4)        :: maxlimit
    REAL(4)        :: missing(2)
    REAL(4)        :: noise
    REAL(8)        :: cycleoffset_x
    REAL(8)        :: cycleoffset_y
    REAL(8)        :: cycleoffset_z
    INTEGER        :: subregion
    INTEGER        :: fileview
    CHARACTER(16)  :: interval
    CHARACTER(16)  :: start, start_date
    CHARACTER(16)  :: end,   end_date
    CHARACTER(16)  :: histshift
    CHARACTER(16)  :: histcycle
    CHARACTER(512) :: inputdir
    CHARACTER(32)  :: suffix
    LOGICAL        :: omit_time
    LOGICAL        :: inquire_exist
    LOGICAL        :: report

    CHARACTER(512) :: default_inputdir
    CHARACTER(8)   :: default_precision
    REAL(4)        :: default_missing(2)
    LOGICAL        :: default_descend
    INTEGER        :: default_fileview
    CHARACTER(20)  :: default_mode
    INTEGER        :: default_periods
    CHARACTER(16)  :: default_interval
    CHARACTER(16)  :: default_start
    CHARACTER(16)  :: default_end
    CHARACTER(16)  :: default_histshift
    CHARACTER(16)  :: default_histcycle
    LOGICAL        :: default_omit_time
    LOGICAL        :: default_inquire

    INTEGER        :: iostat
    CHARACTER(256) :: iomsg

    INTEGER :: imode
    INTEGER :: n, i

    NAMELIST / input /  &
         varname,       &
         filename,      &
         inputdir,      &
         suffix,        &
         precision,     &
         descend,       &
         mode,          &
         periods,       &
         scale,         &
         offset,        &
         minlimit,      &
         maxlimit,      &
         missing,       &
         noise,         &
         cycleoffset_x, &
         cycleoffset_y, &
         cycleoffset_z, &
         interval,      &
         start,         &
         start_date,    &
         end,           &
         end_date,      &
         histshift,     &
         histcycle,     &
         subregion,     &
         fileview,      &
         omit_time,     &
         inquire_exist, &
         report

    CALL default_input(default_inputdir  = default_inputdir,  &
                       default_precision = default_precision, &
                       default_missing   = default_missing,   &
                       default_descend   = default_descend,   &
                       default_fileview  = default_fileview,  &
                       default_mode      = default_mode,      &
                       default_periods   = default_periods,   &
                       default_interval  = default_interval,  &
                       default_histshift = default_histshift, &
                       default_histcycle = default_histcycle, &
                       default_start     = default_start,     &
                       default_end       = default_end,       &
                       default_omit_time   = default_omit_time,   &
                       default_inquire   = default_inquire)

    IF (rank==0) REWIND(CONFIG_UNIT)
    DO
       IF (rank==0) THEN
          varname     = ''
          filename(:) = ''
          inputdir    = default_inputdir
          precision   = default_precision
          mode        = default_mode
          periods     = default_periods
          descend     = default_descend
          fileview    = default_fileview
          suffix      = ""
          scale       = 1.0
          offset      = 0.0
          minlimit    = UNDEF
          maxlimit    = UNDEF
          missing     = default_missing
          noise       = 0.0
          interval    = default_interval
          start       = default_start
          end         = default_end
          histshift   = default_histshift
          histcycle   = default_histcycle
          start_date  = ''
          end_date    = ''
          subregion   = 0
          cycleoffset_x = 0.0
          cycleoffset_y = 0.0
          cycleoffset_z = 0.0
          omit_time     = default_omit_time
          inquire_exist = default_inquire
#ifdef DEBUG
          report      = .TRUE.
#else
          report      = .FALSE.
#endif

          READ(CONFIG_UNIT, NML=input, IOSTAT=iostat, IOMSG=iomsg)
          CALL assert(iostat <= 0, "failed to read INPUT namelist for VARNAME='"//trim(varname)//"'", iomsg)
       END IF

       CALL bcast(iostat)
       IF (iostat < 0) EXIT

       IF (rank==0) THEN
          CALL assert(trim(varname) /= '', "VARNAME is mandatory for INPUT namelist")
          CALL assert(valid_name(varname), "invalid varname '"//trim(varname)//"' in INPUT namelist")
          CALL assert(periods > 0,         "PERIODS in INPUT namelist should be a positive integer")
          CALL assert(periods <= nmax,     "PERIODS in INPUT namelist exceeds the limit")

          IF (start_date /= '') start = start_date
          IF (end_date   /= '') end   = end_date

          CALL replace_vars(inputdir, default=default_inputdir)

          SELECT CASE (trim(mode))
          CASE ('CONST', 'CONSTANT', 'const', 'constant')
             IF (periods == 1) THEN
                imode = 0
                IF (trim(filename(0)) == '') filename(0) = trim(varname)
             ELSE
                imode = 10
                DO n=0, periods-1
                   IF (trim(filename(n)) == '') filename(n) = trim(varname) // '.' // trim(format(n, 'I0'))
                END DO
             END IF

          CASE ('LINEAR', 'LIN', 'linear', 'lin')
             IF (periods == 1) THEN
                imode = 1

                IF (trim(filename(0)) == '') filename(0) = trim(varname) // '.0'
                IF (trim(filename(1)) == '') filename(1) = trim(varname) // '.1'
             ELSE
                imode = 11
                DO n=0, periods-1
                   IF (trim(filename(n)) == '') filename(n) = trim(varname) // '.' // trim(format(n, 'I0'))
                END DO
             END IF

          CASE ('SIGMOID', 'SIG', 'sigmoid', 'sig')
             IF (periods == 1) THEN
                imode = 2

                IF (trim(filename(0)) == '') filename(0) = trim(varname) // '.0'
                IF (trim(filename(1)) == '') filename(1) = trim(varname) // '.1'
             ELSE
                imode = 12
                DO n=0, periods-1
                   IF (trim(filename(n)) == '') filename(n) = trim(varname) // '.' // trim(format(n, 'I0'))
                END DO
             END IF

          CASE ('HIST', 'HISTORICAL', 'hist', 'historical')
             imode = 20
             periods = 1
             IF (trim(filename(0)) == '') filename(0) = trim(varname)

          CASE ('HISTLIN', 'HISTORICAL-LINEAR', 'histlin', 'historical-linear')
             imode = 21
             periods = 1
             IF (trim(filename(0)) == '') filename(0) = trim(varname)

          CASE ('HISTSIG', 'HISTORICAL-SIGMOID', 'histsig', 'historical-sigmoid')
             imode = 22
             periods = 1
             IF (trim(filename(0)) == '') filename(0) = trim(varname)

          CASE ('SIN', 'SINE', 'SINE-CURVE', 'sin', 'sine', 'sine-curve')
             imode = 30
             periods = 2
             IF (trim(filename(0)) == '') filename(0) = trim(varname)

          CASE ('GAUSSIAN', 'GAUSS', 'gaussian', 'gauss')
             imode = 40
             periods = 1
             IF (trim(filename(0)) == '') filename(0) = trim(varname)

          CASE DEFAULT
             CALL assert(.FALSE., "input mode '"//trim(mode)//"' is not supported")
          END SELECT

          DO n=0, nmax-1
             CALL replace('$RUNNAME', trim(runname), filename(n))
          END DO

          IF (subregion /= 0) THEN
             CALL assert(subregion > 0 .AND. subregion <= MAX_SUBREGION, "invalid SUBREGION")
             CALL assert(regions(subregion)%defined, "SUBREGION ID="//trim(format(subregion))//" is not defined")
          END IF

          IF (fileview /= 0) THEN
#ifdef MPIIO
             CALL assert(fileview > 0 .AND. fileview <= MAX_FILEVIEW, "invalid FILEVIEW")
#else
             CALL assert(.FALSE., "FILEVIEW of INPUT is supported only on MPI-IO environment")
#endif
          END IF

          IF (cycleoffset_x/=0.0 .OR. cycleoffset_y/=0.0 .OR. cycleoffset_z/=0.0) THEN
             CALL assert(imode < 30,   "CYCLEOFFSET is not allowed for SINE/GAUSSIAN mode input")
             CALL assert(subregion==0, "CYCLEOFFSET is not allowed for subregion input")
          END IF
       END IF

       CALL bcast(varname)
       CALL bcast(imode)
       CALL bcast(periods)
       CALL bcast(descend)
       CALL bcast(scale)
       CALL bcast(offset)
       CALL bcast(minlimit)
       CALL bcast(maxlimit)
       CALL bcast(missing)
       CALL bcast(noise)
       CALL bcast(cycleoffset_x)
       CALL bcast(cycleoffset_y)
       CALL bcast(cycleoffset_z)
       CALL bcast(interval)
       CALL bcast(start)
       CALL bcast(end)
       CALL bcast(histshift)
       CALL bcast(histcycle)
       CALL bcast(inputdir)
       CALL bcast(precision)
       CALL bcast(subregion)
       CALL bcast(fileview)
       CALL bcast(report)
       CALL bcast(omit_time)
       CALL bcast(inquire_exist)

       DO n=0, nmax-1
          CALL bcast(filename(n))
       END DO
       CALL bcast(suffix)

       CALL assert(n_input < MAX_INPUT, "number of INPUT exceeds the MAX_INPUT")

       n_input = n_input + 1

       input_registry(n_input)%mode      = imode
       input_registry(n_input)%varname   = trim(varname)
       input_registry(n_input)%inputdir  = trim(inputdir)
       input_registry(n_input)%suffix    = trim(suffix)
       input_registry(n_input)%interval  = interval_seconds(interval)
       input_registry(n_input)%start     = datetime_seconds(start)

       IF (end(1:1)=='+') THEN
          input_registry(n_input)%end = input_registry(n_input)%start + interval_seconds(end(2:))
       ELSE
          input_registry(n_input)%end = datetime_seconds(end)
       END IF
       IF (trim(histshift) == "") THEN
          input_registry(n_input)%tshift = 0.0
       ELSE
          SELECT CASE (histshift(1:1))
          CASE ('+')
             input_registry(n_input)%tshift = interval_seconds(histshift(2:))
          CASE ('-')
             input_registry(n_input)%tshift =-interval_seconds(histshift(2:))
          CASE DEFAULT
             input_registry(n_input)%tshift = interval_seconds(histshift)
          END SELECT
       END IF
       IF (trim(histcycle) == "") THEN
          input_registry(n_input)%tcycle = 0.0
       ELSE
          input_registry(n_input)%tcycle = interval_seconds(histcycle)
       END IF
       CALL assert(input_registry(n_input)%tcycle >= 0.0, "invalid HISTCYCLE")


       input_registry(n_input)%periods   = periods
       input_registry(n_input)%descend   = descend
       input_registry(n_input)%scale     = scale
       input_registry(n_input)%offset    = offset
       input_registry(n_input)%minlimit  = minlimit
       input_registry(n_input)%maxlimit  = maxlimit
       input_registry(n_input)%missing   = missing
       input_registry(n_input)%noise     = noise
       input_registry(n_input)%offset_cx = cycleoffset_x
       input_registry(n_input)%offset_cy = cycleoffset_y
       input_registry(n_input)%offset_cz = cycleoffset_z
       input_registry(n_input)%region    = subregion
       input_registry(n_input)%view      = fileview
       input_registry(n_input)%report    = report
       input_registry(n_input)%omit_time = omit_time

       SELECT CASE (trim(precision))
          CASE ('BYTE',   'byte',   'INT1',  'int1',  'i1', 'I1', '1')
             input_registry(n_input)%kind = 1
          CASE ('SINGLE', 'single', 'REAL4', 'real4', 'r4', 'R4', '4')
             input_registry(n_input)%kind = 4
          CASE ('DOUBLE', 'double', 'REAL8', 'real8', 'r8', 'R8', '8')
             input_registry(n_input)%kind = 8
          CASE DEFAULT
             CALL assert(.FALSE., "unsupported PRECISION '"//trim(precision)//"'")
       END SELECT

       ALLOCATE(input_registry(n_input)%filename(0:max(1,periods-1)))
       input_registry(n_input)%filename(0:max(1,periods-1)) = filename(0:max(1,periods-1))

       IF (rank==0) THEN
          SELECT CASE (imode)
          CASE (0)
             WRITE(REPORT_UNIT, *) "register input for '"//trim(varname)//"', constant mode"
          CASE (1)
             WRITE(REPORT_UNIT, *) "register input for '"//trim(varname)//"', linear-interpolation mode"
          CASE (2)
             WRITE(REPORT_UNIT, *) "register input for '"//trim(varname)//"', sigmoid-interpolation mode"
          CASE (10)
             WRITE(REPORT_UNIT, *) "register input for '"//trim(varname)//"', cyclic mode with constant value for each period"
          CASE (11)
             WRITE(REPORT_UNIT, *) "register input for '"//trim(varname)//"', cyclic mode with linear interpolation"
          CASE (12)
             WRITE(REPORT_UNIT, *) "register input for '"//trim(varname)//"', cyclic mode with sigmoid interpolation"
          CASE (20)
             WRITE(REPORT_UNIT, *) "register input for '"//trim(varname)//"', historical mode with constant value for each period"
          CASE (21)
             WRITE(REPORT_UNIT, *) "register input for '"//trim(varname)//"', historical mode with linear interpolation"
          CASE (22)
             WRITE(REPORT_UNIT, *) "register input for '"//trim(varname)//"', historical mode with sigmoid interpolation"
          CASE (30)
             WRITE(REPORT_UNIT, *) "register input for '"//trim(varname)//"', sine-curve mode"
          CASE (40)
             WRITE(REPORT_UNIT, *) "register input for '"//trim(varname)//"', gaussian mode"
          CASE DEFAULT
             CALL assert(.FALSE., "invalid input mode")
          END SELECT
       END IF

       IF (inquire_exist) CALL inquire_input_files(input_registry(n_input))
    END DO

    input_registry_index(:) = 0
    DO n=1, n_input
       i = iachar(input_registry(n)%varname(1:1))
       input_registry_index(i) = input_registry_index(i) + 1
    END DO

    input_registry_index(128) = sum(input_registry_index(1:127)) + 1
    DO i=127, 0, -1
       input_registry_index(i) = input_registry_index(i+1) - input_registry_index(i)
    END DO

    IF (n_input > 0) CALL merge_sort(input_registry(1:n_input))

#ifdef DEBUG
    IF (rank==0) THEN
       WRITE(STDERR_UNIT,*) '-----INPUT LIST------'
       DO n=1, n_input
          WRITE(STDERR_UNIT,*) n, ': ', input_registry(n)%varname
       END DO
       WRITE(STDERR_UNIT,*) '---------------------'
    END IF
#endif

  CONTAINS
    SUBROUTINE inquire_input_files(reg)
      TYPE(input_registry_entry), INTENT(INOUT) :: reg
      INTEGER :: n
      CHARACTER(1024) :: file
      LOGICAL :: exist
      REAL(8) :: t

      IF (rank/=0) RETURN

      IF (reg%mode >= 20 .AND. reg%mode < 30) THEN !historical mode

         t = reg%start - reg%interval

         DO WHILE (t < min(reg%end, t_end))
            t = t + reg%interval

            IF (t <= t_start - reg%interval) CYCLE
            IF (reg%mode==20 .AND. t >= min(reg%end, t_end)) EXIT

            IF (reg%tcycle /= 0) t = mod(t, reg%tcycle)

            file = truncate_atmark(path(reg%inputdir, trim(reg%filename(0))//'.'//format_datetime(t+reg%tshift, omit_time=reg%omit_time), reg%suffix))
            CALL replace_datetime(file, t+reg%tshift)
            INQUIRE (FILE=file, EXIST=exist)
            CALL assert(exist, "file '"//trim(file)//"' required for '"//trim(reg%varname)//"' input is not found")
         END DO
      ELSE
         DO n=0, size(reg%filename(:))-1
            IF (trim(reg%filename(n))=='') CYCLE
            IF (check_dollar(reg%filename(n))) CYCLE

            file = truncate_atmark(path(reg%inputdir, reg%filename(n), reg%suffix))
            INQUIRE (FILE=file, EXIST=exist)
            CALL assert(exist, "file '"//trim(file)//"' required for '"//trim(reg%varname)//"' input is not found")
         END DO
      END IF

    END SUBROUTINE inquire_input_files

    RECURSIVE SUBROUTINE merge_sort(reg, work)
      TYPE(input_registry_entry), INTENT(INOUT) :: reg(:)
      TYPE(input_registry_entry), INTENT(INOUT), OPTIONAL, TARGET :: work(:)

      TYPE(input_registry_entry), POINTER :: ptr(:)
      TYPE(input_registry_entry) :: tmp

      INTEGER :: n, m
      INTEGER :: i, j, k

      n = size(reg)

      IF (n==1) RETURN

      IF (n==2) THEN
         IF (lgt(reg(1)%varname, reg(2)%varname)) THEN
            tmp    = reg(2)
            reg(2) = reg(1)
            reg(1) = tmp
         END IF
         RETURN
      END IF

      IF (.NOT. present(work)) THEN
         ALLOCATE(ptr(n))
      ELSE
#ifdef DEBUG
         CALL assert(n == size(work), "invalid work array in READ_INPUT_NAMELIST%MERGE_SORT")
#endif
         ptr => work
      END IF

      m = INT(n / 2)
      CALL merge_sort(reg(1:m),   ptr(1:m))
      CALL merge_sort(reg(m+1:n), ptr(m+1:n))

      ptr(1:m)   = reg(1:m)
      ptr(m+1:n) = reg(m+1:n)

      j=1
      k=1
      DO i=1, n
         IF (m+k > n) THEN
            reg(i:n) = ptr(j:m)
            EXIT
         ELSE IF (j > m) THEN
            reg(i:n) = ptr(m+k:n)
            EXIT
         ELSE IF (lgt(ptr(j)%varname, ptr(m+k)%varname)) THEN
            reg(i) = ptr(m+k)
            k = k+1
         ELSE
            reg(i) = ptr(j)
            j = j+1
         END IF
      END DO

      IF (.NOT. present(work)) DEALLOCATE(ptr)
    END SUBROUTINE merge_sort

  END SUBROUTINE read_input_namelist

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE read_output_namelist
    CHARACTER(32)  :: varname
    CHARACTER(512) :: outputdir
    CHARACTER(8)   :: precision
    CHARACTER(16)  :: start, start_date
    CHARACTER(16)  :: end,   end_date
    CHARACTER(16)  :: interval
    INTEGER        :: subregion
    INTEGER        :: slice_x
    INTEGER        :: slice_y
    INTEGER        :: slice_z
    LOGICAL        :: surface
    LOGICAL        :: bottom
    LOGICAL        :: int_x
    LOGICAL        :: int_y
    LOGICAL        :: int_z
    INTEGER        :: int_span(0:1)
    LOGICAL        :: dif_x
    LOGICAL        :: dif_y
    LOGICAL        :: dif_z
    LOGICAL        :: mean
    CHARACTER(16)  :: mean_span
    LOGICAL        :: omit_time

    CHARACTER(512) :: default_outputdir
    CHARACTER(8)   :: default_precision
    CHARACTER(16)  :: default_start
    CHARACTER(16)  :: default_end
    CHARACTER(16)  :: default_interval
    LOGICAL        :: default_omit_time

    INTEGER        :: iostat
    CHARACTER(256) :: iomsg

    INTEGER :: n

    NAMELIST / output / &
         varname,       &
         outputdir,     &
         precision,     &
         start,         &
         start_date,    &
         end,           &
         end_date,      &
         interval,      &
         mean,          &
         mean_span,     &
         subregion,&
         slice_x,       &
         slice_y,       &
         slice_z,       &
         surface,       &
         bottom,        &
         int_x,         &
         int_y,         &
         int_z,         &
         int_span,      &
         dif_x,         &
         dif_y,         &
         dif_z,         &
         omit_time

    NAMELIST / output_restart / &
         start, end, interval, outputdir

    CALL default_output(default_outputdir = default_outputdir, &
                        default_precision = default_precision, &
                        default_start     = default_start,     &
                        default_end       = default_end,       &
                        default_interval  = default_interval,  &
                        default_omit_time = default_omit_time)

    IF (rank==0) REWIND(CONFIG_UNIT)
    DO
       IF (rank==0) THEN
          varname   = ''
          outputdir = default_outputdir
          precision = default_precision
          start     = default_start
          end       = default_end
          start_date= ''
          end_date  = ''
          interval  = default_interval
          mean      = .FALSE.
          mean_span = ''
          subregion = 0
          slice_x   = -1
          slice_y   = -1
          slice_z   = -1
          surface   = .FALSE.
          bottom    = .FALSE.
          int_x     = .FALSE.
          int_y     = .FALSE.
          int_z     = .FALSE.
          int_span  = -1
          dif_x     = .FALSE.
          dif_y     = .FALSE.
          dif_z     = .FALSE.
          omit_time = default_omit_time

          READ(CONFIG_UNIT, NML=output, IOSTAT=iostat, IOMSG=iomsg)
          CALL assert(iostat <= 0, "failed to read OUTPUT namelist for VARNAME='"//trim(varname)//"'", iomsg)
       END IF

       CALL bcast(iostat)
       IF (iostat < 0) EXIT

       IF (rank==0) THEN
          CALL assert(trim(varname) /= '', "VARNAME is mandatory for OUTPUT namelist")
          CALL assert(valid_name(varname), "invalid VARNAME '"//trim(varname)//"' in OUTPUT namelist")

          !for backward compatibility
          IF (start_date /= '') start = start_date
          IF (end_date   /= '') end   = end_date

          CALL replace_vars(outputdir, default=default_outputdir)

          CALL assert(valid_dir(outputdir), "OUTPUTDIR '"//trim(outputdir)//"' dose not exist")
       END IF

       CALL bcast(varname)
       CALL bcast(outputdir)
       CALL bcast(precision)
       CALL bcast(start)
       CALL bcast(end)
       CALL bcast(interval)
       CALL bcast(mean)
       CALL bcast(mean_span)
       CALL bcast(subregion)
       CALL bcast(slice_x)
       CALL bcast(slice_y)
       CALL bcast(slice_z)
       CALL bcast(surface)
       CALL bcast(bottom)
       CALL bcast(int_x)
       CALL bcast(int_y)
       CALL bcast(int_z)
       CALL bcast(int_span)
       CALL bcast(dif_x)
       CALL bcast(dif_y)
       CALL bcast(dif_z)
       CALL bcast(omit_time)

       CALL assert(n_output < MAX_OUTPUT, "number of OUTPUT exceeds the MAX_OUTPUT")

       n = 0
       IF (slice_x >= 0)       n = n+1
       IF (slice_y >= 0)       n = n+1
       IF (slice_z >= 0)       n = n+1
       IF (surface)            n = n+1
       IF (bottom)             n = n+1
       IF (int_x)              n = n+1
       IF (int_y)              n = n+1
       IF (int_z)              n = n+1

       CALL assert(n <= 1, "SLICE_X/Y/Z, INT_X/Y/Z and SURFACE/BOTTOM flags cannot be specified simultaneously")

       IF (subregion /= 0) THEN
          CALL assert(subregion > 0 .AND. subregion <= MAX_SUBREGION, "invalid SUBREGION")
          CALL assert(regions(subregion)%defined, "SUBREGION="//trim(format(subregion))//" is not defined")
       END IF

       CALL assert(slice_x >= -1 .AND. slice_x <= dimx, "SLICE_X is out of range")
       CALL assert(slice_y >= -1 .AND. slice_y <= dimy, "SLICE_Y is out of range")
       CALL assert(slice_z >= -1 .AND. slice_z <= dimz, "SLICE_Z is out of range")

       IF (int_x) THEN
          IF (int_span(0)<=0) int_span(0) = 1
          IF (int_span(1)<=0) int_span(1) = dimx
          int_span(1) = min(int_span(1), dimx)
       END IF
       IF (int_y) THEN
          IF (int_span(0)<=0) int_span(0) = 1
          IF (int_span(1)<=0) int_span(1) = dimy
          int_span(1) = min(int_span(1), dimy)
       END IF
       IF (int_z) THEN
          IF (int_span(0)<=0) int_span(0) = 1
          IF (int_span(1)<=0) int_span(1) = dimz
          int_span(1) = min(int_span(1), dimz)
       END IF
       CALL assert(int_span(0) <= int_span(1), "INT_SPAN is out of range")

       n_output = n_output + 1

       output_registry(n_output)%varname    = trim(varname)
       output_registry(n_output)%outputdir  = trim(outputdir)
       output_registry(n_output)%start      = datetime_seconds(start)
       IF (end(1:1)=='+') THEN
          output_registry(n_output)%end = output_registry(n_output)%start + interval_seconds(end(2:))
       ELSE
          output_registry(n_output)%end = datetime_seconds(end)
       END IF
       output_registry(n_output)%interval   = interval_seconds(interval)
       IF (omit_time) CALL assert(interval_seconds(interval) >= 86400.0, "Cannot specify True for OMIT_TIME when output interval is less than one day")
       output_registry(n_output)%omit_time  = omit_time
       output_registry(n_output)%mean       = mean
       output_registry(n_output)%meancount  = 0
       IF (mean_span=='') THEN
          output_registry(n_output)%meanspan = interval_seconds(interval)
       ELSE
          output_registry(n_output)%meanspan = interval_seconds(mean_span)
       END IF
       output_registry(n_output)%lastwrite   = datetime_seconds(start)

       SELECT CASE (trim(precision))
          CASE ('BYTE',   'byte',   'INT1',  'int1',  'i1', 'I1', '1')
             output_registry(n_output)%kind = 1
          CASE ('SINGLE', 'single', 'REAL4', 'real4', 'r4', 'R4', '4')
             output_registry(n_output)%kind = 4
          CASE ('DOUBLE', 'double', 'REAL8', 'real8', 'r8', 'R8', '8')
             output_registry(n_output)%kind = 8
          CASE DEFAULT
             CALL assert(.FALSE., "unsupported PRECISION '"//trim(precision)//"'")
       END SELECT

       IF (slice_x >= 0) THEN
          output_registry(n_output)%mode = 1
          output_registry(n_output)%slice = slice_x
          output_registry(n_output)%basename = trim(varname) // '_X' // trim(format(slice_x, 'I0'))
       ELSE IF (slice_y >= 0) THEN
          output_registry(n_output)%mode = 2
          output_registry(n_output)%slice = slice_y
          output_registry(n_output)%basename = trim(varname) // '_Y' // trim(format(slice_y, 'I0'))
       ELSE IF (slice_z >= 0) THEN
          output_registry(n_output)%mode = 3
          output_registry(n_output)%slice = slice_z
          output_registry(n_output)%basename = trim(varname) // '_Z' // trim(format(slice_z, 'I0'))
       ELSE IF (surface) THEN
          output_registry(n_output)%mode = 4
          output_registry(n_output)%basename = trim(varname) // '_SFC'
       ELSE IF (bottom) THEN
          output_registry(n_output)%mode = 5
          output_registry(n_output)%basename = trim(varname) // '_BTM'
       ELSE IF (int_x) THEN
          output_registry(n_output)%mode = 6
          output_registry(n_output)%intspan(:) = int_span(:)
          output_registry(n_output)%basename = trim(varname) // '_XINT'
          IF (int_span(0) /= 1 .OR. int_span(1) /= dimx) THEN
             output_registry(n_output)%basename = trim(output_registry(n_output)%basename) &
                  // trim(format(int_span(0), 'I0')) // '-' // trim(format(int_span(1), 'I0'))
          END IF
       ELSE IF (int_y) THEN
          output_registry(n_output)%mode = 7
          output_registry(n_output)%intspan(:) = int_span(:)
          output_registry(n_output)%basename = trim(varname) // '_YINT'
          IF (int_span(0) /= 1 .OR.  int_span(1) /= dimy) THEN
             output_registry(n_output)%basename = trim(output_registry(n_output)%basename) &
                  // trim(format(int_span(0), 'I0')) // '-' // trim(format(int_span(1), 'I0'))
          END IF
       ELSE IF (int_z) THEN
          output_registry(n_output)%mode = 8
          output_registry(n_output)%intspan(:) = int_span(:)
          output_registry(n_output)%basename = trim(varname) // '_ZINT'
          IF (int_span(0) /= 1 .OR.  int_span(1) /= dimz) THEN
             output_registry(n_output)%basename = trim(output_registry(n_output)%basename) &
                  // trim(format(int_span(0), 'I0')) // '-' // trim(format(int_span(1), 'I0'))
          END IF
       ELSE
          output_registry(n_output)%mode = 0
          output_registry(n_output)%basename = trim(varname)
       END IF

       output_registry(n_output)%region = subregion
       IF (subregion /= 0) THEN
          output_registry(n_output)%basename = trim(output_registry(n_output)%basename) // '_RGN' // trim(format(subregion, 'I0'))
       END IF

       IF (rank==0) THEN
          WRITE(REPORT_UNIT, '(A)', ADVANCE='NO') "register output for '"//trim(output_registry(n_output)%basename)//"'"
          IF (output_registry(n_output)%mean) THEN
             WRITE(REPORT_UNIT, '(A)', ADVANCE='NO') ", time-mean "
          ELSE
             WRITE(REPORT_UNIT, '(A)', ADVANCE='NO') ", snap-shot "
          END IF
          WRITE(REPORT_UNIT, '(A)', ADVANCE='NO') "with interval '"//trim(interval)//"'"
          SELECT CASE(output_registry(n_output)%kind)
          CASE (8)
             WRITE(REPORT_UNIT, '(A)', ADVANCE='NO') ", double-precision"
          CASE (4)
             WRITE(REPORT_UNIT, '(A)', ADVANCE='NO') ", single-precision"
          CASE (1)
             WRITE(REPORT_UNIT, '(A)', ADVANCE='NO') ", 8bit-integer"
          END SELECT
          WRITE(REPORT_UNIT, *)
       END IF

       output_registry(n_output)%dimcode     = 0
       output_registry(n_output)%checked     = .FALSE.
       output_registry(n_output)%initialized = .FALSE.
    END DO

    CALL organize_output_registry

    IF (rank==0) THEN
       outputdir = ''
       start     = default_start
       end       = default_end
       interval  = ''

       REWIND(CONFIG_UNIT)
       READ(CONFIG_UNIT, NML=output_restart, IOSTAT=iostat, IOMSG=iomsg)
       CALL assert(iostat <= 0, "failed to read OUTPUT_RESTART namelist", iomsg)

       IF (iostat==0) THEN
          CALL assert(outputdir/='', "OUTPUTDIR is mandatory for OUTPUT_RESTART namelist")
          CALL assert(interval /='',  "INTERVAL is mandatory for OUTPUT_RESTART namelist")

          CALL replace_vars(outputdir, default=default_outputdir)

          CALL assert(valid_dir(outputdir), "OUTPTUDIR '"//trim(outputdir)//"' does not exist")
       END IF

       restart_outputdir = outputdir
       restart_start     = start
       restart_end       = end
       restart_interval  = interval
    END IF
    CALL bcast(restart_outputdir)
    CALL bcast(restart_start)
    CALL bcast(restart_end)
    CALL bcast(restart_interval)

#ifdef DEBUG
    IF (rank==0) THEN
       WRITE(STDERR_UNIT,*) '------OUTPUT LIST-------'
       DO n=1, n_output
          WRITE(STDERR_UNIT,*) n, ': ', output_registry(n)%varname
       END DO
       WRITE(STDERR_UNIT,*) '------------------------'
    END IF
#endif

  END SUBROUTINE read_output_namelist

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE organize_output_registry
    INTEGER :: i, n, m

    IF (n_output == 0) RETURN

    IF (rank==0) THEN
       DO m=2, n_output
          DO n=1, m-1
            CALL assert((output_registry(n)%basename   /= output_registry(m)%basename)   .OR. &
                        (output_registry(n)%outputdir  /= output_registry(m)%outputdir)  .OR. &
                        (output_registry(n)%mean   .NEQV. output_registry(m)%mean)       .OR. &
                        (output_registry(n)%end        <= output_registry(m)%start)      .OR. &
                        (output_registry(n)%start      >= output_registry(m)%end),            &
                        "output for '"//trim(output_registry(n)%basename)//"' is overlaped."  &
                        // " Please specify different OUTPUTDIR to avoid filename collision.")
          END DO
       END DO
    END IF

    output_registry_index(:) = 0
    DO n=1, n_output
       i = iachar(output_registry(n)%varname(1:1))
       output_registry_index(i) = output_registry_index(i) + 1
    END DO

    output_registry_index(128) = sum(output_registry_index(1:127)) + 1
    DO i=127, 0, -1
       output_registry_index(i) = output_registry_index(i+1) - output_registry_index(i)
    END DO

    CALL merge_sort(output_registry(1:n_output))

  CONTAINS
    RECURSIVE SUBROUTINE merge_sort(reg, work)
      TYPE(output_registry_entry), INTENT(INOUT) :: reg(:)
      TYPE(output_registry_entry), INTENT(INOUT), OPTIONAL, TARGET :: work(:)

      TYPE(output_registry_entry), POINTER :: ptr(:)
      TYPE(output_registry_entry) :: tmp

      INTEGER :: n, m
      INTEGER :: i, j, k

      n = size(reg)

      IF (n==1) RETURN

      IF (n==2) THEN
         IF (lgt(reg(1)%varname, reg(2)%varname)) THEN
            tmp    = reg(2)
            reg(2) = reg(1)
            reg(1) = tmp
         END IF
         RETURN
      END IF

      IF (.NOT. present(work)) THEN
         ALLOCATE(ptr(n))
      ELSE
#ifdef DEBUG
         CALL assert(n == size(work), "invalid work array in READ_OUTPUT_NAMELIST%MERGE_SORT")
#endif
         ptr => work
      END IF

      m = INT(n / 2)
      CALL merge_sort(reg(1:m),   ptr(1:m))
      CALL merge_sort(reg(m+1:n), ptr(m+1:n))

      ptr(1:m)   = reg(1:m)
      ptr(m+1:n) = reg(m+1:n)

      j=1
      k=1
      DO i=1, n
         IF (m+k > n) THEN
            reg(i:n) = ptr(j:m)
            EXIT
         ELSE IF (j > m) THEN
            reg(i:n) = ptr(m+k:n)
            EXIT
         ELSE IF (lgt(ptr(j)%varname, ptr(m+k)%varname)) THEN
            reg(i) = ptr(m+k)
            k = k+1
         ELSE
            reg(i) = ptr(j)
            j = j+1
         END IF
      END DO

      IF (.NOT. present(work)) DEALLOCATE(ptr)
    END SUBROUTINE merge_sort

  END SUBROUTINE organize_output_registry

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE add_restart_var(varname, kind)
    CHARACTER(*), INTENT(IN) :: varname
    INTEGER,      INTENT(IN), OPTIONAL :: kind

    INTEGER :: kind_

    IF (restart_interval == '') RETURN

    kind_ = 8
    IF (present(kind)) kind_ = kind

    n_output = n_output + 1
    output_registry(n_output)%varname   = trim(varname)
    output_registry(n_output)%basename  = trim(varname)
    output_registry(n_output)%outputdir = restart_outputdir
    output_registry(n_output)%kind      = kind_
    output_registry(n_output)%mode      = 0
    output_registry(n_output)%region    = 0
    output_registry(n_output)%mean      = .FALSE.
    output_registry(n_output)%start     = datetime_seconds(restart_start)
    IF (restart_end(1:1)=='+') THEN
       output_registry(n_output)%end = output_registry(n_output)%start + interval_seconds(restart_end(2:))
    ELSE
       output_registry(n_output)%end = datetime_seconds(restart_end)
    END IF
    output_registry(n_output)%interval    = interval_seconds(restart_interval)
    output_registry(n_output)%lastwrite   = datetime_seconds(restart_start)
    output_registry(n_output)%meanspan    = interval_seconds(restart_interval)
    output_registry(n_output)%meancount   = 0
    output_registry(n_output)%checked     = .FALSE.
    output_registry(n_output)%initialized = .FALSE.

    CALL organize_output_registry

  END SUBROUTINE add_restart_var

!-----------------------------------------------------------------------------------------------------------------------

  LOGICAL PURE FUNCTION require_checkin_by_name(varname, time)
    CHARACTER(*), INTENT(IN) :: varname
    REAL(8),      INTENT(IN), OPTIONAL :: time

    INTEGER :: id

    id = lookup_input_id(varname)

    DO WHILE (id /= 0)
       IF (require_checkin_by_id(id, time)) THEN
          require_checkin_by_name = .TRUE.
          RETURN
       END IF

       id = lookup_input_id(varname, skip=id)
    END DO

    require_checkin_by_name = .FALSE.

  END FUNCTION require_checkin_by_name

!-----------------------------------------------------------------------------------------------------------------------

  LOGICAL PURE FUNCTION require_checkin_by_id(id, time)
    INTEGER, INTENT(IN) :: id
    REAL(8), INTENT(IN), OPTIONAL :: time

    REAL(8) :: t

    IF (present(time)) THEN
       t = time
    ELSE
       t = t_current
    END IF

    require_checkin_by_id = (t >= input_registry(id)%start .AND. t <  input_registry(id)%end)

  END FUNCTION require_checkin_by_id

!-----------------------------------------------------------------------------------------------------------------------

  LOGICAL PURE FUNCTION input_registered(varname)
    CHARACTER(*), INTENT(IN) :: varname

    input_registered = (lookup_input_id(varname) /= 0)
  END FUNCTION input_registered

!-----------------------------------------------------------------------------------------------------------------------

  LOGICAL PURE FUNCTION output_registered(varname)
    CHARACTER(*), INTENT(IN) :: varname

    output_registered = (lookup_output_id(varname) /= 0)
  END FUNCTION output_registered

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE checkin_1d_r4(varname, data, stat, add, time, axis)
    CHARACTER(*), INTENT(IN)    :: varname
    REAL(4),      INTENT(INOUT) :: data(:)
    LOGICAL,      INTENT(INOUT), OPTIONAL :: stat
    LOGICAL,      INTENT(IN),    OPTIONAL :: add
    REAL(8),      INTENT(IN),    OPTIONAL :: time
    CHARACTER(1), INTENT(IN),    OPTIONAL :: axis
#ifdef F2008
    CONTIGUOUS data
#endif

    REAL(8) :: tmp(size(data))

    IF (.NOT. require_checkin(varname, time)) THEN
       IF (present(stat)) THEN
          IF (present(add)) THEN
             stat = add .AND. stat
          ELSE
             stat = .FALSE.
          END IF
       END IF
       RETURN
    END IF

!$OMP PARALLEL WORKSHARE
    tmp(:) = REAL(data(:), KIND=8)
!$OMP END PARALLEL WORKSHARE

    CALL checkin_1d(varname, tmp, stat, add, time, axis)

!$OMP PARALLEL WORKSHARE
    data(:) = REAL(tmp(:), KIND=4)
!$OMP END PARALLEL WORKSHARE

  END SUBROUTINE checkin_1d_r4

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE checkin_1d_logical(varname, data, stat, time, axis)
    CHARACTER(*), INTENT(IN)    :: varname
    LOGICAL,      INTENT(INOUT) :: data(:)
    LOGICAL,      INTENT(INOUT), OPTIONAL :: stat
    REAL(8),      INTENT(IN),    OPTIONAL :: time
    CHARACTER(1), INTENT(IN),    OPTIONAL :: axis
#ifdef F2008
    CONTIGUOUS data
#endif

    REAL(8) :: tmp(size(data))

    IF (.NOT. require_checkin(varname, time)) THEN
       IF (present(stat)) stat = .FALSE.
       RETURN
    END IF

!$OMP PARALLEL WORKSHARE
    tmp(:) = 0.0
!$OMP END PARALLEL WORKSHARE

    CALL checkin_1d(varname, tmp, stat, add=.FALSE., time=time, axis=axis)

!$OMP PARALLEL WORKSHARE
    data(:) = tmp(:) /= 0.0
!$OMP END PARALLEL WORKSHARE

  END SUBROUTINE checkin_1d_logical

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE checkin_1d_i1(varname, data, stat, add, time, axis)
    CHARACTER(*), INTENT(IN)    :: varname
    INTEGER(1),   INTENT(INOUT) :: data(:)
    LOGICAL,      INTENT(INOUT), OPTIONAL :: stat
    LOGICAL,      INTENT(IN),    OPTIONAL :: add
    REAL(8),      INTENT(IN),    OPTIONAL :: time
    CHARACTER(1), INTENT(IN),    OPTIONAL :: axis
#ifdef F2008
    CONTIGUOUS data
#endif

    REAL(8) :: tmp(size(data))

    IF (.NOT. require_checkin(varname, time)) THEN
       IF (present(stat)) THEN
          IF (present(add)) THEN
             stat = add .AND. stat
          ELSE
             stat = .FALSE.
          END IF
       END IF
       RETURN
    END IF

!$OMP PARALLEL WORKSHARE
    tmp(:) = REAL(data(:), KIND=8)
!$OMP END PARALLEL WORKSHARE

    CALL checkin_1d(varname, tmp, stat, add, time, axis)

!$OMP PARALLEL WORKSHARE
    data(:) = INT(tmp(:), KIND=1)
!$OMP END PARALLEL WORKSHARE

  END SUBROUTINE checkin_1d_i1

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE checkin_1d(varname, data, stat, add, time, axis)
    CHARACTER(*), INTENT(IN)    :: varname
    REAL(8),      INTENT(INOUT) :: data(:)
    LOGICAL,      INTENT(INOUT), OPTIONAL :: stat
    LOGICAL,      INTENT(IN),    OPTIONAL :: add
    REAL(8),      INTENT(IN),    OPTIONAL :: time
    CHARACTER(1), INTENT(IN),    OPTIONAL :: axis
#ifdef F2008
    CONTIGUOUS data
#endif

    REAL(8) :: tmp(lbound(data,1):ubound(data,1), 1, 1)

    INTEGER :: id
    INTEGER :: dimcode
    REAL(8) :: time_
    LOGICAL :: add_, stat_, report_

    add_ = .FALSE.
    IF (present(add)) add_ = add

    IF (.NOT. require_checkin(varname, time)) THEN
       IF (present(stat)) stat = stat .AND. add_
       RETURN
    END IF

    stat_   = .FALSE.
    report_ = .FALSE.

    IF (present(time)) THEN
       time_ = time
    ELSE
       time_ = t_current
    END IF

    dimcode = -3
    IF (present(axis)) THEN
       SELECT CASE(axis)
       CASE ('X', 'x')
          dimcode = -1
       CASE ('Y', 'y')
          dimcode = -2
       CASE ('Z', 'z')
          dimcode = -3
       CASE DEFAULT
          CALL assert(.FALSE., "invalid AXIS in CHECKIN_1D")
       END SELECT
    END IF

    id = lookup_input_id(varname)

!$OMP PARALLEL WORKSHARE
    tmp(:,1,1) = data(:)
!$OMP END PARALLEL WORKSHARE

    DO WHILE (id /= 0)
       IF (require_checkin(id, time_)) CALL checkin_private(input_registry(id), tmp, stat_, add_, time_, dimcode)
       add_    = add_    .OR. stat_
       report_ = report_ .OR. input_registry(id)%report

       id = lookup_input_id(varname, skip=id)
    END DO

    IF (stat_) THEN
!$OMP PARALLEL WORKSHARE
       data(:) = tmp(:,1,1)
!$OMP END PARALLEL WORKSHARE
    END IF

    IF (stat_ .AND. report_) CALL checkin_report(varname, data, dimcode)

    IF (present(stat)) stat = stat_ .OR. (add_ .AND. stat)

  END SUBROUTINE checkin_1d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE checkin_2d_r4(varname, data, stat, add, time, section)
    CHARACTER(*), INTENT(IN)    :: varname
    REAL(4),      INTENT(INOUT) :: data(:,:)
    LOGICAL,      INTENT(INOUT), OPTIONAL :: stat
    LOGICAL,      INTENT(IN),    OPTIONAL :: add
    REAL(8),      INTENT(IN),    OPTIONAL :: time
    CHARACTER(2), INTENT(IN),    OPTIONAL :: section
#ifdef F2008
    CONTIGUOUS data
#endif

    REAL(8) :: tmp(size(data,1),size(data,2))

    IF (.NOT. require_checkin(varname, time)) THEN
       IF (present(stat)) THEN
          IF (present(add)) THEN
             stat = add .AND. stat
          ELSE
             stat = .FALSE.
          END IF
       END IF
       RETURN
    END IF

!$OMP PARALLEL WORKSHARE
    tmp(:,:) = REAL(data(:,:), KIND=8)
!$OMP END PARALLEL WORKSHARE

    CALL checkin_2d(varname, tmp, stat, add, time, section)

!$OMP PARALLEL WORKSHARE
    data(:,:) = REAL(tmp(:,:), KIND=4)
!$OMP END PARALLEL WORKSHARE

  END SUBROUTINE checkin_2d_r4

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE checkin_2d_i1(varname, data, stat, add, time, section)
    CHARACTER(*), INTENT(IN)    :: varname
    INTEGER(1),   INTENT(INOUT) :: data(:,:)
    LOGICAL,      INTENT(INOUT), OPTIONAL :: stat
    LOGICAL,      INTENT(IN),    OPTIONAL :: add
    REAL(8),      INTENT(IN),    OPTIONAL :: time
    CHARACTER(2), INTENT(IN),    OPTIONAL :: section
#ifdef F2008
    CONTIGUOUS data
#endif

    REAL(8) :: tmp(size(data,1),size(data,2))

    IF (.NOT. require_checkin(varname, time)) THEN
       IF (present(stat)) THEN
          IF (present(add)) THEN
             stat = add .AND. stat
          ELSE
             stat = .FALSE.
          END IF
       END IF
       RETURN
    END IF

!$OMP PARALLEL WORKSHARE
    tmp(:,:) = REAL(data(:,:), KIND=8)
!$OMP END PARALLEL WORKSHARE

    CALL checkin_2d(varname, tmp, stat, add, time, section)

!$OMP PARALLEL WORKSHARE
    data(:,:) = INT(tmp(:,:), KIND=1)
!$OMP END PARALLEL WORKSHARE

  END SUBROUTINE checkin_2d_i1

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE checkin_2d_logical(varname, data, stat, time, section)
    CHARACTER(*), INTENT(IN)    :: varname
    LOGICAL,      INTENT(INOUT) :: data(:,:)
    LOGICAL,      INTENT(INOUT), OPTIONAL :: stat
    REAL(8),      INTENT(IN),    OPTIONAL :: time
    CHARACTER(2), INTENT(IN),    OPTIONAL :: section
#ifdef F2008
    CONTIGUOUS data
#endif

    REAL(8) :: tmp(size(data,1),size(data,2))

    IF (.NOT. require_checkin(varname, time)) THEN
       IF (present(stat)) stat = .FALSE.
       RETURN
    END IF

!$OMP PARALLEL WORKSHARE
    tmp(:,:) = 0.0
!$OMP END PARALLEL WORKSHARE

    CALL checkin_2d(varname, tmp, stat, add=.FALSE., time=time, section=section)

!$OMP PARALLEL WORKSHARE
    data(:,:) = tmp(:,:) /= 0.0
!$OMP END PARALLEL WORKSHARE

  END SUBROUTINE checkin_2d_logical

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE checkin_2d(varname, data, stat, add, time, section)
    CHARACTER(*), INTENT(IN)    :: varname
    REAL(8),      INTENT(INOUT) :: data(:,:)
    LOGICAL,      INTENT(INOUT), OPTIONAL :: stat
    LOGICAL,      INTENT(IN),    OPTIONAL :: add
    REAL(8),      INTENT(IN),    OPTIONAL :: time
    CHARACTER(2), INTENT(IN),    OPTIONAL :: section

    REAL(8) :: tmp(lbound(data,1):ubound(data,1), lbound(data,2):ubound(data,2), 1)

    INTEGER :: id
    INTEGER :: dimcode
    REAL(8) :: time_
    LOGICAL :: add_, stat_, report_

    add_ = .FALSE.
    IF (present(add)) add_ = add

    IF (.NOT. require_checkin(varname, time)) THEN
       IF (present(stat)) stat = stat.AND. add_
       RETURN
    END IF

    stat_   = .FALSE.
    report_ = .FALSE.

    IF (present(time)) THEN
       time_ = time
    ELSE
       time_ = t_current
    END IF

    dimcode = 1
    IF (present(section)) THEN
       SELECT CASE(section)
       CASE ('XY', 'xy')
          dimcode = 1
       CASE ('XZ', 'xz')
          dimcode = 2
       CASE ('YZ', 'yz')
          dimcode = 3
       CASE DEFAULT
          CALL assert(.FALSE., "invalid SECTION in CHECKIN_2D")
       END SELECT
    END IF

    id = lookup_input_id(varname)

!$OMP PARALLEL WORKSHARE
    tmp(:,:,1) = data(:,:)
!$OMP END PARALLEL WORKSHARE

    DO WHILE (id /= 0)
       IF (require_checkin(id, time_)) CALL checkin_private(input_registry(id), tmp, stat_, add_, time_, dimcode)
       add_    = add_    .OR. stat_
       report_ = report_ .OR. input_registry(id)%report

       id = lookup_input_id(varname, skip=id)
    END DO

    IF (stat_) THEN
!$OMP PARALLEL WORKSHARE
       data(:,:) = tmp(:,:,1)
!$OMP END PARALLEL WORKSHARE
    END IF

    IF (stat_ .AND. report_) CALL checkin_report(varname, data, dimcode)

    IF (present(stat)) stat = stat_ .OR. (add_ .AND. stat)

  END SUBROUTINE checkin_2d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE checkin_private(reg, data, stat, add, time, dimcode)
    TYPE(input_registry_entry), INTENT(INOUT) :: reg
    REAL(8), INTENT(INOUT) :: data(:,:,:)
    LOGICAL, INTENT(INOUT) :: stat
    LOGICAL, INTENT(IN)    :: add
    REAL(8), INTENT(IN)    :: time
    INTEGER, INTENT(IN)    :: dimcode
#ifdef F2008
    CONTIGUOUS data
#endif

    REAL(8) :: tmp(size(data,1), size(data,2), size(data,3))

    REAL(8) :: alpha, beta
    REAL(8) :: t

    INTEGER :: n1, n2, n3
    INTEGER :: slv1, slv2, slv3

    CHARACTER(1024) :: filepath

    INTEGER :: i, j, k

    n1 = size(data,1)
    n2 = size(data,2)
    n3 = size(data,3)

    IF (add) THEN
!$OMP PARALLEL DO COLLAPSE(2)
       DO k=1, n3
       DO j=1, n2
       DO i=1, n1
          tmp(i,j,k) = data(i,j,k)
       END DO
       END DO
       END DO
    END IF

    SELECT CASE (dimcode)
    CASE (0) !3d
       CALL assert(n1-isize >= 0 .AND. n1-isize <= 5, "unsupported dimension for CHECKIN with DIMCODE=0 (3D)")
       CALL assert(n2-jsize >= 0 .AND. n2-jsize <= 5, "unsupported dimension for CHECKIN with DIMCODE=0 (3D)")
       CALL assert(n3-ksize >= 0 .AND. n3-ksize <= 5, "unsupported dimension for CHECKIN with DIMCODE=0 (3D)")

       slv1 = (n1-isize+1)/2
       slv2 = (n2-jsize+1)/2
       slv3 = (n3-ksize+1)/2

    CASE (1) !2d (xy-plane)
       CALL assert(n1-isize >= 0 .AND. n1-isize <= 5, "unsupported dimension for CHECKIN with DIMCODE=1 (2D)")
       CALL assert(n2-jsize >= 0 .AND. n2-jsize <= 5, "unsupported dimension for CHECKIN with DIMCODE=1 (2D)")
       CALL assert(n3==1,                             "dimension conflicts for CHECKIN with DIMCODE=1 (2D)")

       slv1 = (n1-isize+1)/2
       slv2 = (n2-jsize+1)/2
       slv3 = 0

    CASE (2) !xz-plane
       CALL assert(n1-isize >= 0 .AND. n1-isize <= 5, "unsupported dimension for CHECKIN with DIMCODE=2 (2D, XZ-section)")
       CALL assert(n2-ksize >= 0 .AND. n2-ksize <= 5, "unsupported dimension for CHECKIN with DIMCODE=2 (2D, XZ-section)")
       CALL assert(n3==1,                             "dimension conflicts for CHECKIN with DIMCODE=2 (2D, XZ-section)")

       slv1 = (n1-isize+1)/2
       slv2 = (n2-ksize+1)/2
       slv3 = 0

    CASE (3) !yz-plane
       CALL assert(n1-jsize >= 0 .AND. n1-jsize <= 5, "unsupported dimension for CHECKIN with DIMCODE=3 (2D, YZ-section)")
       CALL assert(n2-ksize >= 0 .AND. n2-ksize <= 5, "unsupported dimension for CHECKIN with DIMCODE=3 (2D, YZ-section)")
       CALL assert(n3==1,                             "dimension conflicts for CHECKIN with DIMCODE=3 (2D, YZ-section)")

       slv1 = (n1-jsize+1)/2
       slv2 = (n2-ksize+1)/2
       slv3 = 0

    CASE (-1) !1d x-axis
       CALL assert(n1-isize >= 0 .AND. n1-isize <= 5, "unsupported dimension for CHECKIN with DIMCODE=-1 (1D, X-axis)")
       CALL assert(n2==1,                             "dimension conflicts for CHECKIN with DIMCODE=-1 (1D, X-axis)")
       CALL assert(n3==1,                             "dimension conflicts for CHECKIN with DIMCODE=-1 (1D, X-axis)")

       slv1 = (n1-isize+1)/2
       slv2 = 0
       slv3 = 0

    CASE (-2) !1d y-axis
       CALL assert(n1-jsize >= 0 .AND. n1-jsize <= 5, "unsupported dimension for CHECKIN with DIMCODE=-2 (1D, Y-axis)")
       CALL assert(n2==1,                             "dimension conflicts for CHECKIN with DIMCODE=-2 (1D, Y-axis)")
       CALL assert(n3==1,                             "dimension conflicts for CHECKIN with DIMCODE=-2 (1D, Y-axis)")

       slv1 = (n1-jsize+1)/2
       slv2 = 0
       slv3 = 0

    CASE (-3) !1d z-axis
       CALL assert(n1-ksize >= 0 .AND. n1-ksize <= 5, "unsupported dimension for CHECKIN with DIMCODE=-3 (1D, Z-axis)")
       CALL assert(n2==1,                             "dimension conflicts for CHECKIN with DIMCODE=-3 (1D, Z-axis)")
       CALL assert(n3==1,                             "dimension conflicts for CHECKIN with DIMCODE=-3 (1D, Z-axis)")

       slv1 = (n1-ksize+1)/2
       slv2 = 0
       slv3 = 0

    CASE DEFAULT
       CALL assert(.FALSE., "invalid DIMCODE in CHECKIN_PRIVATE")

    END SELECT

    IF (.NOT. reg%initialized) CALL init_input

    SELECT CASE (reg%mode)
    CASE (0) ! constant mode
!$OMP PARALLEL DO COLLAPSE(2)
       DO k=1, n3
       DO j=1, n2
       DO i=1, n1
          data(i,j,k) = reg%buf(i-slv1,j-slv2,k-slv3,0)
       END DO
       END DO
       END DO

    CASE (1,2) ! interpolation mode
       alpha = (time - reg%start)/(reg%end - reg%start)  !for linear  interpolation
       IF (reg%mode==2) alpha = logistic(12*(alpha-0.5)) !for sigmoid interpolation

!$OMP PARALLEL DO COLLAPSE(2)
       DO k=1, n3
       DO j=1, n2
       DO i=1, n1
          IF (reg%buf(i-slv1,j-slv2,k-slv3,0)==UNDEF .OR. reg%buf(i-slv1,j-slv2,k-slv3,1)==UNDEF) THEN
             data(i,j,k) = UNDEF
          ELSE
             data(i,j,k) = reg%buf(i-slv1,j-slv2,k-slv3,0) * (1.0-alpha) &
                         + reg%buf(i-slv1,j-slv2,k-slv3,1) * alpha
          END IF
       END DO
       END DO
       END DO

       CALL limitminmax(data, reg%minlimit, reg%maxlimit)

    CASE (10) ! cyclic-constant mode
       IF (time - reg%lastread >= reg%interval) THEN
          reg%current = mod(reg%current+1, reg%periods)
          reg%lastread = reg%lastread + reg%interval

          CALL read_buf(0, path(reg%inputdir, reg%filename(reg%current), reg%suffix))

          CALL substmissing(reg%buf(:,:,:,0), reg%missing(1), reg%missing(2))
          CALL scaleoffset(reg%buf(:,:,:,0), reg%scale, reg%offset)
          CALL limitminmax(reg%buf(:,:,:,0), reg%minlimit, reg%maxlimit)
       END IF

!$OMP PARALLEL DO COLLAPSE(2)
       DO k=1, n3
       DO j=1, n2
       DO i=1, n1
          data(i,j,k) = reg%buf(i-slv1,j-slv2,k-slv3,0)
       END DO
       END DO
       END DO

    CASE (11,12) ! cyclic with interpolation
       IF (time - reg%lastread >= reg%interval) THEN
          reg%current = mod(reg%current+1, reg%periods)
          reg%lastread = reg%lastread + reg%interval

!$OMP PARALLEL WORKSHARE
          reg%buf(:,:,:,0) = reg%buf(:,:,:,1)
!$OMP END PARALLEL WORKSHARE

          CALL read_buf(1, path(reg%inputdir, reg%filename(mod(reg%current+1, reg%periods)), reg%suffix))

          CALL substmissing(reg%buf(:,:,:,1), reg%missing(1), reg%missing(2))
          CALL scaleoffset(reg%buf(:,:,:,1), reg%scale, reg%offset)
       END IF

       alpha = (time - reg%lastread)/(reg%interval)      !for linear  interpolation
       IF (reg%mode==12) alpha = logistic(12*(alpha-0.5)) !for sigmoid interpolation

!$OMP PARALLEL DO COLLAPSE(2)
       DO k=1, n3
       DO j=1, n2
       DO i=1, n1
          IF (reg%buf(i-slv1,j-slv2,k-slv3,0)==UNDEF .OR. reg%buf(i-slv1,j-slv2,k-slv3,1)==UNDEF) THEN
             data(i,j,k) = UNDEF
          ELSE
             data(i,j,k) = reg%buf(i-slv1,j-slv2,k-slv3,0) * (1.0-alpha) &
                         + reg%buf(i-slv1,j-slv2,k-slv3,1) * alpha
          END IF
       END DO
       END DO
       END DO

       CALL limitminmax(data, reg%minlimit, reg%maxlimit)

    CASE (20) ! historical mode with constant value
       IF (time - reg%lastread >= reg%interval) THEN
          reg%lastread = reg%lastread + reg%interval

          t = reg%lastread
          IF (reg%tcycle > 0.0) t = mod(t, reg%tcycle)
          t = t + reg%tshift
          filepath = path(reg%inputdir, trim(reg%filename(0))//'.'//format_datetime(t, omit_time=reg%omit_time), reg%suffix)
          CALL replace_datetime(filepath, t)
          CALL read_buf(0, filepath)

          CALL substmissing(reg%buf(:,:,:,0), reg%missing(1), reg%missing(2))
          CALL scaleoffset(reg%buf(:,:,:,0), reg%scale, reg%offset)
          CALL limitminmax(reg%buf(:,:,:,0), reg%minlimit, reg%maxlimit)
       END IF

!$OMP PARALLEL DO COLLAPSE(2)
       DO k=1, n3
       DO j=1, n2
       DO i=1, n1
          data(i,j,k) = reg%buf(i-slv1,j-slv2,k-slv3,0)
       END DO
       END DO
       END DO

    CASE (21,22) ! historical mode with interpolation
       IF (time - reg%lastread >= reg%interval) THEN
          reg%lastread = reg%lastread + reg%interval

!$OMP PARALLEL WORKSHARE
          reg%buf(:,:,:,0) = reg%buf(:,:,:,1)
!$OMP END PARALLEL WORKSHARE

          t = reg%lastread + reg%interval
          IF (reg%tcycle > 0.0) t = mod(t, reg%tcycle)
          t = t + reg%tshift
          filepath = path(reg%inputdir, trim(reg%filename(0))//'.'//format_datetime(t, omit_time=reg%omit_time), reg%suffix)
          CALL replace_datetime(filepath, t)

          CALL read_buf(1, filepath)

          CALL substmissing(reg%buf(:,:,:,1), reg%missing(1), reg%missing(2))
          CALL scaleoffset(reg%buf(:,:,:,1), reg%scale, reg%offset)
       END IF

       alpha = (time - reg%lastread)/(reg%interval)      !for linear  interpolation
       IF (reg%mode==22) alpha = logistic(12*(alpha-0.5)) !for sigmoid interpolation

!$OMP PARALLEL DO COLLAPSE(2)
       DO k=1, n3
       DO j=1, n2
       DO i=1, n1
          IF (reg%buf(i-slv1,j-slv2,k-slv3,0)==UNDEF .OR. reg%buf(i-slv1,j-slv2,k-slv3,1)==UNDEF) THEN
             data(i,j,k) = UNDEF
          ELSE
             data(i,j,k) = reg%buf(i-slv1,j-slv2,k-slv3,0) * (1.0-alpha) &
                         + reg%buf(i-slv1,j-slv2,k-slv3,1) * alpha
          END IF
       END DO
       END DO
       END DO

       CALL limitminmax(data, reg%minlimit, reg%maxlimit)

    CASE (30) ! sin-curve mode
       alpha = sin(2*pi*(time - reg%start) / (reg%interval))
       beta  = cos(2*pi*(time - reg%start) / (reg%interval))

!$OMP PARALLEL DO COLLAPSE(2)
       DO k=1, n3
       DO j=1, n2
       DO i=1, n1
          IF (reg%buf(i-slv1,j-slv2,k-slv3,0)==UNDEF) THEN
             data(i,j,k) = UNDEF
          ELSE
             data(i,j,k) = reg%offset + reg%scale*reg%buf(i-slv1,j-slv2,k-slv3,0) * (reg%buf(i-slv1,j-slv2,k-slv3,1)*beta &
                                                                                   + reg%buf(i-slv1,j-slv2,k-slv3,2)*alpha)
          END IF
       END DO
       END DO
       END DO

       CALL limitminmax(data, reg%minlimit, reg%maxlimit)

    CASE (40) ! gaussian mode
       alpha = exp(-pi*(3*((time - reg%start)/(reg%end - reg%start)-0.5))**2)

!$OMP PARALLEL DO COLLAPSE(2)
       DO k=1, n3
       DO j=1, n2
       DO i=1, n1
          IF (reg%buf(i-slv1,j-slv2,k-slv3,0)==UNDEF) THEN
             data(i,j,k) = UNDEF
          ELSE
             data(i,j,k) = reg%offset + reg%scale*reg%buf(i-slv1,j-slv2,k-slv3,0)*alpha
          END IF
       END DO
       END DO
       END DO

       CALL limitminmax(data, reg%minlimit, reg%maxlimit)

    END SELECT


    SELECT CASE (dimcode)
    CASE (0) !3d
       IF (reg%noise /= 0.0) THEN
          CALL whitenoise(data(slv1+1:slv1+isize,slv2+1:slv2+jsize,slv3+1:slv3+ksize), reg%noise, &
                          mkseed(r8=t_current, i4=kcoord*ipes*jpes+jcoord*ipes+icoord, str=reg%varname))
          CALL update_boundary_3d(data)
       END IF

       IF (cycle_x .AND. reg%offset_cx /= 0.0) THEN
          IF (icoord==0 .AND. slv1 > 0) THEN
!$OMP PARALLEL DO
             DO k=1, n3
             DO j=1, n2
             DO i=1, slv1
                IF (data(i,j,k) == UNDEF) CYCLE
                data(i,j,k) = data(i,j,k) - reg%offset_cx
             END DO
             END DO
             END DO
          END IF

          IF (icoord==ipes-1 .AND. n1 > slv1+isize) THEN
!$OMP PARALLEL DO
             DO k=1, n3
             DO j=1, n2
             DO i=slv1+isize+1, n1
                IF (data(i,j,k) == UNDEF) CYCLE
                data(i,j,k) = data(i,j,k) + reg%offset_cx
             END DO
             END DO
             END DO
          END IF
       END IF

       IF (cycle_y .AND. reg%offset_cy /= 0.0) THEN
          IF (jcoord==0 .AND. slv2 > 0) THEN
!$OMP PARALLEL DO
             DO k=1, n3
             DO j=1, slv2
             DO i=1, n1
                IF (data(i,j,k) == UNDEF) CYCLE
                data(i,j,k) = data(i,j,k) - reg%offset_cy
             END DO
             END DO
             END DO
          END IF

          IF (jcoord==jpes-1 .AND. n2 > slv2+jsize) THEN
!$OMP PARALLEL DO
             DO k=1, n3
             DO j=slv2+jsize+1, n2
             DO i=1, n1
                IF (data(i,j,k) == UNDEF) CYCLE
                data(i,j,k) = data(i,j,k) + reg%offset_cy
             END DO
             END DO
             END DO
          END IF
       END IF

       IF (cycle_z .AND. reg%offset_cz /= 0.0) THEN
          IF (kcoord==0 .AND. slv3 > 0) THEN
!$OMP PARALLEL DO
             DO k=1, slv3
             DO j=1, n2
             DO i=1, n1
                IF (data(i,j,k) == UNDEF) CYCLE
                data(i,j,k) = data(i,j,k) - reg%offset_cz
             END DO
             END DO
             END DO
          END IF

          IF (kcoord==kpes-1 .AND. n3 > +slv3+ksize) THEN
!$OMP PARALLEL DO
             DO k=slv3+ksize+1, n3
             DO j=1, n2
             DO i=1, n1
                IF (data(i,j,k) == UNDEF) CYCLE
                data(i,j,k) = data(i,j,k) + reg%offset_cz
             END DO
             END DO
             END DO
          END IF
       END IF

    CASE (1) !xy-plane
       IF (reg%noise /= 0.0) THEN
         CALL whitenoise(data(slv1+1:slv1+isize,slv2+1:slv2+jsize,1), reg%noise, &
                         mkseed(r8=t_current, i4=jcoord*ipes+icoord, str=reg%varname))

          CALL update_boundary_2d(data(:,:,1))
       END IF

       IF (cycle_x .AND. reg%offset_cx /= 0.0) THEN
          IF (icoord==0 .AND. slv1 > 0) THEN
!$OMP PARALLEL DO COLLAPSE(2)
             DO j=1, n2
             DO i=1, slv1
                IF (data(i,j,1) == UNDEF) CYCLE
                data(i,j,1) = data(i,j,1) - reg%offset_cx
             END DO
             END DO
          END IF

          IF (icoord==ipes-1 .AND. n1 > slv1+isize) THEN
!$OMP PARALLEL DO COLLAPSE(2)
             DO j=1, n2
             DO i=slv1+isize+1, n1
                IF (data(i,j,1) == UNDEF) CYCLE
                data(i,j,1) = data(i,j,1) + reg%offset_cx
             END DO
             END DO
          END IF
       END IF

       IF (cycle_y .AND. reg%offset_cy /= 0.0) THEN
          IF (jcoord==0 .AND. slv2 > 0) THEN
!$OMP PARALLEL DO COLLAPSE(2)
             DO j=1, slv2
             DO i=1, n1
                IF (data(i,j,1) == UNDEF) CYCLE
                data(i,j,1) = data(i,j,1) - reg%offset_cy
             END DO
             END DO
          END IF

          IF (jcoord==jpes-1 .AND. n2 > slv2+jsize) THEN
!$OMP PARALLEL DO COLLAPSE(2)
             DO j=slv2+jsize+1, n2
             DO i=1, n1
                IF (data(i,j,1) == UNDEF) CYCLE
                data(i,j,1) = data(i,j,1) + reg%offset_cy
             END DO
             END DO
          END IF
       END IF

    CASE (2) !xz-plane
       IF (reg%noise /= 0.0) THEN
          CALL whitenoise(data(slv1+1:slv1+isize,slv2+1:slv2+jsize,1), reg%noise, &
                          mkseed(r8=t_current, i4=kcoord*ipes+icoord, str=reg%varname))
          CALL update_boundary_xz(data(:,:,1))
       END IF

       IF (cycle_x .AND. reg%offset_cx /= 0.0) THEN
          IF (icoord==0 .AND. slv1 > 0) THEN
!$OMP PARALLEL DO COLLAPSE(2)
             DO k=1, n2
             DO i=1, slv1
                IF (data(i,k,1) == UNDEF) CYCLE
                data(i,k,1) = data(i,k,1) - reg%offset_cx
             END DO
             END DO
          END IF

          IF (icoord==ipes-1 .AND. n1 > slv1+isize) THEN
!$OMP PARALLEL DO COLLAPSE(2)
             DO k=1, n2
             DO i=slv1+isize+1, n1
                IF (data(i,k,1) == UNDEF) CYCLE
                data(i,k,1) = data(i,k,1) + reg%offset_cx
             END DO
             END DO
          END IF
       END IF

       IF (cycle_z .AND. reg%offset_cz /= 0.0) THEN
          IF (kcoord==0 .AND. slv2 > 0) THEN
!$OMP PARALLEL DO COLLAPSE(2)
             DO k=1, slv2
             DO i=1, n1
                IF (data(i,k,1) == UNDEF) CYCLE
                data(i,k,1) = data(i,k,1) - reg%offset_cz
             END DO
             END DO
          END IF

          IF (kcoord==kpes-1 .AND. n2 > slv2+ksize) THEN
!$OMP PARALLEL DO COLLAPSE(2)
             DO k=slv2+ksize+1, n2
             DO i=1, n1
                IF (data(i,k,1) == UNDEF) CYCLE
                data(i,k,1) = data(i,k,1) + reg%offset_cz
             END DO
             END DO
          END IF
       END IF

    CASE (3) !yz-plane
       IF (reg%noise /= 0.0) THEN
          CALL whitenoise(data(slv1+1:slv1+isize,slv2+1:slv2+jsize,1), reg%noise, &
                          mkseed(r8=t_current, i4=kcoord*jpes+jcoord, str=reg%varname))
          CALL update_boundary_yz(data(:,:,1))
       END IF

       IF (cycle_y .AND. reg%offset_cy /= 0.0) THEN
          IF (jcoord==0 .AND. slv1 > 0) THEN
!$OMP PARALLEL DO COLLAPSE(2)
             DO k=1, n2
             DO j=1, slv1
                IF (data(j,k,1) == UNDEF) CYCLE
                data(j,k,1) = data(j,k,1) - reg%offset_cy
             END DO
             END DO
          END IF

          IF (jcoord==jpes-1 .AND. n1 > slv1+jsize) THEN
!$OMP PARALLEL DO COLLAPSE(2)
             DO k=1, n2
             DO j=slv1+jsize+1, n1
                IF (data(j,k,1) == UNDEF) CYCLE
                data(j,k,1) = data(j,k,1) + reg%offset_cy
             END DO
             END DO
          END IF
       END IF

       IF (cycle_z .AND. reg%offset_cz /= 0.0) THEN
          IF (kcoord==0 .AND. slv2 > 0) THEN
!$OMP PARALLEL DO COLLAPSE(2)
             DO k=1, slv2
             DO j=1, n1
                IF (data(j,k,1) == UNDEF) CYCLE
                data(j,k,1) = data(j,k,1) - reg%offset_cz
             END DO
             END DO
          END IF

          IF (kcoord==kpes-1 .AND. n2 > slv2+ksize) THEN
!$OMP PARALLEL DO COLLAPSE(2)
             DO k=slv2+ksize+1, n2
             DO j=1, n1
                IF (data(j,k,1) == UNDEF) CYCLE
                data(j,k,1) = data(j,k,1) + reg%offset_cz
             END DO
             END DO
          END IF
       END IF

    CASE (-1) !x-axis
       IF (reg%noise /= 0.0) THEN
          CALL whitenoise(data(slv1+1:slv1+isize,slv2+1:slv2+jsize,1), reg%noise, &
                          mkseed(r8=t_current, i4=icoord, str=reg%varname))
          CALL update_boundary_x(data(:,1,1))
       END IF

       IF (cycle_x .AND. reg%offset_cx /= 0.0) THEN
          IF (icoord==0 .AND. slv1 > 0) THEN
             DO i=1, slv1
                IF (data(i,1,1) == UNDEF) CYCLE
                data(i,1,1) = data(i,1,1) - reg%offset_cx
             END DO
          END IF

          IF (icoord==ipes-1 .AND. n1 > slv1+isize) THEN
             DO i=slv1+isize+1, n1
                IF (data(i,1,1) == UNDEF) CYCLE
                data(i,1,1) = data(i,1,1) + reg%offset_cx
             END DO
          END IF
       END IF

    CASE (-2) !y-axis
       IF (reg%noise /= 0.0) THEN
          CALL whitenoise(data(slv1+1:slv1+isize,slv2+1:slv2+jsize,1), reg%noise, &
                          mkseed(r8=t_current, i4=jcoord, str=reg%varname))
          CALL update_boundary_y(data(:,1,1))
       END IF

       IF (cycle_y .AND. reg%offset_cy /= 0.0) THEN
          IF (jcoord==0 .AND. slv1 > 0) THEN
             DO j=1, slv1
                IF (data(j,1,1) == UNDEF) CYCLE
                data(j,1,1) = data(j,1,1) - reg%offset_cy
             END DO
          END IF

          IF (jcoord==jpes-1 .AND. n1 > slv1+jsize) THEN
             DO j=slv1+jsize+1, n1
                IF (data(j,1,1) == UNDEF) CYCLE
                data(j,1,1) = data(j,1,1) + reg%offset_cy
             END DO
          END IF
       END IF

    CASE (-3) !z-axis
       IF (reg%noise /= 0.0) THEN
          CALL whitenoise(data(slv1+1:slv1+isize,slv2+1:slv2+jsize,1), reg%noise, &
                          mkseed(r8=t_current, i4=kcoord, str=reg%varname))
          CALL update_boundary_z(data(:,1,1))
       END IF

       IF (cycle_z .AND. reg%offset_cz /= 0.0) THEN
          IF (kcoord==0 .AND. slv1 > 0) THEN
             DO k=1, slv1
                IF (data(k,1,1) == UNDEF) CYCLE
                data(k,1,1) = data(k,1,1) - reg%offset_cz
             END DO
          END IF

          IF (kcoord==kpes-1 .AND. n1 > slv1+ksize) THEN
             DO k=slv1+ksize+1, n1
                IF (data(k,1,1) == UNDEF) CYCLE
                data(k,1,1) = data(k,1,1) + reg%offset_cz
             END DO
          END IF
       END IF

    END SELECT

    IF (add) THEN
!$OMP PARALLEL DO COLLAPSE(2)
       DO k=1, n3
       DO j=1, n2
       DO i=1, n1
          IF ( tmp(i,j,k) == UNDEF) CYCLE
          IF (data(i,j,k) == UNDEF) THEN
             data(i,j,k) = tmp(i,j,k)
          ELSE
             data(i,j,k) = tmp(i,j,k) + data(i,j,k)
          END IF
       END DO
       END DO
       END DO
    END IF

    stat = .TRUE.

  CONTAINS
    SUBROUTINE read_buf(n, filepath)
      INTEGER,      INTENT(IN) :: n
      CHARACTER(*), INTENT(IN) :: filepath

      IF (trim(filepath) == '') RETURN

#ifdef MPIIO
      CALL assert(reg%view==0 .OR. dimcode==0 .OR. dimcode==1, "FILEVIEW is not supported for 1D or virtical slice data")
#else
      CALL assert(reg%view==0, "FILEVIEW is suppored only on MPI-IO environment")
#endif

      SELECT CASE(dimcode)
      CASE (0)
         CALL read_data_3d(reg%buf(1:isize,1:jsize,1:ksize,n), trim(filepath), reg%kind, reg%region, reg%view, reg%descend)
         CALL update_boundary_3d(reg%buf(:,:,:,n), fill=0.0D0)

      CASE (1)
         CALL read_data_2d(reg%buf(1:isize,1:jsize,1,n), trim(filepath), reg%kind, reg%region, reg%view)
         CALL update_boundary_2d(reg%buf(:,:,1,n), fill=0.0D0)

      CASE (2)
         CALL read_data_xz(reg%buf(1:isize,1:ksize,1,n), trim(filepath), reg%kind, reg%region, reg%descend)
         CALL update_boundary_xz(reg%buf(:,:,1,n), fill=0.0D0)

      CASE (3)
         CALL read_data_yz(reg%buf(1:jsize,1:ksize,1,n), trim(filepath), reg%kind, reg%region, reg%descend)
         CALL update_boundary_yz(reg%buf(:,:,1,n), fill=0.0D0)

      CASE (-1)
         CALL read_data_x(reg%buf(1:isize,1,1,n), trim(filepath), reg%kind, reg%region)
         CALL update_boundary_x(reg%buf(:,1,1,n), fill=0.0D0)

      CASE (-2)
         CALL read_data_y(reg%buf(1:jsize,1,1,n), trim(filepath), reg%kind, reg%region)
         CALL update_boundary_y(reg%buf(:,1,1,n), fill=0.0D0)

      CASE (-3)
         CALL read_data_z(reg%buf(1:ksize,1,1,n), trim(filepath), reg%kind, reg%region, reg%descend)
         CALL update_boundary_z(reg%buf(:,1,1,n), fill=0.0D0)

      END SELECT

    END SUBROUTINE read_buf

    SUBROUTINE init_buf(num)
      INTEGER, INTENT(IN) :: num

      SELECT CASE (dimcode)
      CASE (0) !3d
         ALLOCATE(reg%buf(-slv:isize+slv, -slv:jsize+slv, -slv:ksize+slv, 0:num-1))
      CASE (1) !xy
         ALLOCATE(reg%buf(-slv:isize+slv, -slv:jsize+slv, 1, 0:num-1))
      CASE (2) !xz
         ALLOCATE(reg%buf(-slv:isize+slv, -slv:ksize+slv, 1, 0:num-1))
      CASE (3) !yz
         ALLOCATE(reg%buf(-slv:jsize+slv, -slv:ksize+slv, 1, 0:num-1))
      CASE (-1) !x
         ALLOCATE(reg%buf(-slv:isize+slv, 1, 1, 0:num-1))
      CASE (-2) !y
         ALLOCATE(reg%buf(-slv:jsize+slv, 1, 1, 0:num-1))
      CASE (-3) !z
         ALLOCATE(reg%buf(-slv:ksize+slv, 1, 1, 0:num-1))
      END SELECT

!$OMP PARALLEL WORKSHARE
      reg%buf(:,:,:,:) = UNDEF
!$OMP END PARALLEL WORKSHARE
    END SUBROUTINE init_buf

    SUBROUTINE init_input
      INTEGER :: itmp
      REAL(8) :: rtmp

      CHARACTER(1024) :: filepath

      INTEGER :: i, j, k

      IF (reg%initialized) RETURN

      CALL assert(dimcode>=-3 .AND. dimcode<=3, "invalid DIMCODE in INIT_INPUT")

      SELECT CASE (reg%mode)
      CASE (0) ! constant mode
         CALL init_buf(1)

         CALL read_buf(0, path(reg%inputdir, reg%filename(0), reg%suffix))

         CALL substmissing(reg%buf(:,:,:,0), reg%missing(1), reg%missing(2))
         CALL scaleoffset(reg%buf(:,:,:,0), reg%scale, reg%offset)
         CALL limitminmax(reg%buf(:,:,:,0), reg%minlimit, reg%maxlimit)

      CASE (1,2) ! interpolation mode
         CALL init_buf(2)

         CALL read_buf(0, path(reg%inputdir, reg%filename(0), reg%suffix))
         CALL read_buf(1, path(reg%inputdir, reg%filename(1), reg%suffix))

         CALL substmissing(reg%buf(:,:,:,0), reg%missing(1), reg%missing(2))
         CALL substmissing(reg%buf(:,:,:,1), reg%missing(1), reg%missing(2))
         CALL scaleoffset(reg%buf(:,:,:,0), reg%scale, reg%offset)
         CALL scaleoffset(reg%buf(:,:,:,1), reg%scale, reg%offset)

      CASE (10) ! cyclic mode with constant value
         CALL init_buf(1)

         rtmp = time - reg%start
         itmp = int(rtmp/reg%interval)
         reg%current = mod(itmp, reg%periods)
         reg%lastread = reg%start + reg%interval * itmp

         CALL read_buf(0, path(reg%inputdir, reg%filename(reg%current), reg%suffix))

         CALL substmissing(reg%buf(:,:,:,0), reg%missing(1), reg%missing(2))
         CALL scaleoffset(reg%buf(:,:,:,0), reg%scale, reg%offset)
         CALL limitminmax(reg%buf(:,:,:,0), reg%minlimit, reg%maxlimit)

      CASE (11,12) ! cyclic mode with interpolation
         CALL init_buf(2)

         rtmp = time - reg%start
         itmp = int(rtmp/reg%interval)
         reg%current = mod(itmp, reg%periods)
         reg%lastread = reg%start + reg%interval * itmp

         CALL read_buf(0, path(reg%inputdir, reg%filename(reg%current), reg%suffix))
         CALL read_buf(1, path(reg%inputdir, reg%filename(mod(reg%current+1, reg%periods)), reg%suffix))

         CALL substmissing(reg%buf(:,:,:,0), reg%missing(1), reg%missing(2))
         CALL substmissing(reg%buf(:,:,:,1), reg%missing(1), reg%missing(2))
         CALL scaleoffset(reg%buf(:,:,:,0), reg%scale, reg%offset)
         CALL scaleoffset(reg%buf(:,:,:,1), reg%scale, reg%offset)

      CASE (20) ! historical mode with constant value
         CALL init_buf(1)

         rtmp = time - reg%start
         itmp = int(rtmp/reg%interval)
         reg%lastread = reg%start + reg%interval * itmp

         t = reg%lastread
         IF (reg%tcycle > 0.0) t = mod(t, reg%tcycle)
         t = t + reg%tshift
         filepath = path(reg%inputdir, trim(reg%filename(0))//'.'//format_datetime(t, omit_time=reg%omit_time), reg%suffix)
         CALL replace_datetime(filepath, t)
         CALL read_buf(0, filepath)

         CALL substmissing(reg%buf(:,:,:,0), reg%missing(1), reg%missing(2))
         CALL scaleoffset(reg%buf(:,:,:,0), reg%scale, reg%offset)
         CALL limitminmax(reg%buf(:,:,:,0), reg%minlimit, reg%maxlimit)

      CASE (21,22) ! historical mode with interpolation
         CALL init_buf(2)

         rtmp = time - reg%start
         itmp = int(rtmp/reg%interval)
         reg%lastread = reg%start + reg%interval * itmp

         t = reg%lastread
         IF (reg%tcycle > 0.0) t = mod(t, reg%tcycle)
         t = t + reg%tshift
         filepath = path(reg%inputdir, trim(reg%filename(0))//'.'//format_datetime(t, omit_time=reg%omit_time), reg%suffix)
         CALL replace_datetime(filepath, t)
         CALL read_buf(0, filepath)

         t = reg%lastread + reg%interval
         IF (reg%tcycle > 0.0) t = mod(t, reg%tcycle)
         t = t + reg%tshift
         filepath = path(reg%inputdir, trim(reg%filename(0))//'.'//format_datetime(t, omit_time=reg%omit_time), reg%suffix)
         CALL replace_datetime(filepath, t)
         CALL read_buf(1, filepath)

         CALL substmissing(reg%buf(:,:,:,0), reg%missing(1), reg%missing(2))
         CALL substmissing(reg%buf(:,:,:,1), reg%missing(1), reg%missing(2))
         CALL scaleoffset(reg%buf(:,:,:,0), reg%scale, reg%offset)
         CALL scaleoffset(reg%buf(:,:,:,1), reg%scale, reg%offset)

      CASE (30) ! sin-curve mode
         CALL init_buf(3)
         CALL read_buf(0, path(reg%inputdir, reg%filename(0), reg%suffix))

         CALL substmissing(reg%buf(:,:,:,0), reg%missing(1), reg%missing(2))

         IF (reg%filename(1) /= '') THEN
            CALL read_buf(1, path(reg%inputdir, reg%filename(1), reg%suffix))
            CALL substmissing(reg%buf(:,:,:,1), UNDEF, 0.0)
         ELSE
!$OMP PARALLEL WORKSHARE
            reg%buf(:,:,:,1) = 0.0
!$OMP END PARALLEL WORKSHARE
         END IF

!$OMP PARALLEL DO PRIVATE(rtmp)
         DO k=lbound(reg%buf, 3), ubound(reg%buf, 3)
         DO j=lbound(reg%buf, 2), ubound(reg%buf, 2)
         DO i=lbound(reg%buf, 1), ubound(reg%buf, 1)
            rtmp = 2*pi*reg%buf(i,j,k,1)
            reg%buf(i,j,k,1) = sin(rtmp)
            reg%buf(i,j,k,2) = cos(rtmp)
         END DO
         END DO
         END DO

      CASE (40) ! gaussian mode
         CALL init_buf(1)
         CALL read_buf(0, path(reg%inputdir, reg%filename(0), reg%suffix))

         CALL substmissing(reg%buf(:,:,:,0), reg%missing(1), reg%missing(2))
      END SELECT

      reg%initialized = .TRUE.

    END SUBROUTINE init_input

  END SUBROUTINE checkin_private

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE checkin_3d_r4(varname, data, stat, add, time)
    CHARACTER(*), INTENT(IN)    :: varname
    REAL(4),      INTENT(INOUT) :: data(:,:,:)
    LOGICAL,      INTENT(INOUT), OPTIONAL :: stat
    LOGICAL,      INTENT(IN),    OPTIONAL :: add
    REAL(8),      INTENT(IN),    OPTIONAL :: time
#ifdef F2008
    CONTIGUOUS data
#endif

    REAL(8) :: tmp(size(data,1),size(data,2),size(data,3))
    INTEGER :: n1, n2, n3

    IF (.NOT. require_checkin(varname, time)) THEN
       IF (present(stat)) THEN
          IF (present(add)) THEN
             stat = add .AND. stat
          ELSE
             stat = .FALSE.
          END IF
       END IF
       RETURN
    END IF

!$OMP PARALLEL WORKSHARE
    tmp(:,:,:) = REAL(data(:,:,:), KIND=8)
!$OMP END PARALLEL WORKSHARE

    CALL checkin_3d(varname, tmp, stat, add, time)

!$OMP PARALLEL WORKSHARE
    data(:,:,:) = REAL(tmp(:,:,:), KIND=4)
!$OMP END PARALLEL WORKSHARE

  END SUBROUTINE checkin_3d_r4

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE checkin_3d_i1(varname, data, stat, add, time)
    CHARACTER(*), INTENT(IN)    :: varname
    INTEGER(1),   INTENT(INOUT) :: data(:,:,:)
    LOGICAL,      INTENT(INOUT), OPTIONAL :: stat
    LOGICAL,      INTENT(IN),    OPTIONAL :: add
    REAL(8),      INTENT(IN),    OPTIONAL :: time
#ifdef F2008
    CONTIGUOUS data
#endif

    REAL(8) :: tmp(size(data,1),size(data,2),size(data,3))
    INTEGER :: n1, n2, n3

    IF (.NOT. require_checkin(varname, time)) THEN
       IF (present(stat)) THEN
          IF (present(add)) THEN
             stat = add .AND. stat
          ELSE
             stat = .FALSE.
          END IF
       END IF
       RETURN
    END IF

!$OMP PARALLEL WORKSHARE
    tmp(:,:,:) = REAL(data(:,:,:), KIND=8)
!$OMP END PARALLEL WORKSHARE

    CALL checkin_3d(varname, tmp, stat, add, time)

!$OMP PARALLEL WORKSHARE
    data(:,:,:) = INT(tmp(:,:,:), KIND=1)
!$OMP END PARALLEL WORKSHARE

  END SUBROUTINE checkin_3d_i1

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE checkin_3d_logical(varname, data, stat, time)
    CHARACTER(*), INTENT(IN)    :: varname
    LOGICAL,      INTENT(INOUT) :: data(:,:,:)
    LOGICAL,      INTENT(INOUT), OPTIONAL :: stat
    REAL(8),      INTENT(IN),    OPTIONAL :: time
#ifdef F2008
    CONTIGUOUS data
#endif

    REAL(8) :: tmp(size(data,1),size(data,2),size(data,3))

    IF (.NOT. require_checkin(varname, time)) THEN
       IF (present(stat)) stat = .FALSE.
       RETURN
    END IF

!$OMP PARALLEL WORKSHARE
    tmp(:,:,:) = 0.0
!$OMP END PARALLEL WORKSHARE

    CALL checkin_3d(varname, tmp, stat, add=.FALSE., time=time)

!$OMP PARALLEL WORKSHARE
    data(:,:,:) = tmp(:,:,:) /= 0.0
!$OMP END PARALLEL WORKSHARE

  END SUBROUTINE checkin_3d_logical

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE checkin_3d(varname, data, stat, add, time)
    CHARACTER(*), INTENT(IN) :: varname
    REAL(8), INTENT(INOUT)   :: data(:,:,:)
    LOGICAL, INTENT(INOUT), OPTIONAL :: stat
    LOGICAL, INTENT(IN),    OPTIONAL :: add
    REAL(8), INTENT(IN),    OPTIONAL :: time
#ifdef F2008
    CONTIGUOUS data
#endif

    INTEGER :: id
    REAL(8) :: time_
    LOGICAL :: add_, stat_, report_

    add_ = .FALSE.
    IF (present(add)) add_ = add

    IF (.NOT. require_checkin(varname, time)) THEN
       IF (present(stat)) stat = stat .AND. add_
       RETURN
    END IF

    stat_   = .FALSE.
    report_ = .FALSE.

    id = lookup_input_id(varname)

    IF (present(time)) THEN
       time_ = time
    ELSE
       time_ = t_current
    END IF

    DO WHILE (id /= 0)
       IF (require_checkin(id, time_)) CALL checkin_private(input_registry(id), data, stat_, add_, time_, 0)
       add_    = add_    .OR. stat_
       report_ = report_ .OR. input_registry(id)%report

       id = lookup_input_id(varname, skip=id)
    END DO

    IF (stat_ .AND. report_) CALL checkin_report(varname, data)

    IF (present(stat)) stat = stat_ .OR. (add_ .AND. stat)
  END SUBROUTINE checkin_3d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE checkin_report_1d(varname, data, dimcode)
    CHARACTER(*), INTENT(IN) :: varname
    REAL(8),      INTENT(IN) :: data(:)
    INTEGER,      INTENT(IN) :: dimcode

    REAL(8) :: max, min

    INTEGER :: s, n

    SELECT CASE (dimcode)
    CASE (-1)
       s = (size(data) - isize + 1)/2
       n = isize
    CASE (-2)
       s = (size(data) - jsize + 1)/2
       n = jsize
    CASE (-3)
       s = (size(data) - ksize + 1)/2
       n = ksize
    END SELECT

    max = maxval(data(s+1:s+n))
    min = minval(data(s+1:s+n))

    CALL gmax(max)
    CALL gmin(min)

    IF (rank==0) WRITE(REPORT_UNIT, '(A,ES10.3,X,ES10.3)') "input '"//trim(varname)//"':  ", max, min

  END SUBROUTINE checkin_report_1d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE checkin_report_2d(varname, data, dimcode)
    CHARACTER(*), INTENT(IN) :: varname
    REAL(8),      INTENT(IN)  :: data(:,:)
    INTEGER,      INTENT(IN) :: dimcode

    REAL(8) :: max, min

    INTEGER :: s(2), n(2)

    SELECT CASE (dimcode)
    CASE (1)
       s(1) = (size(data,1) - isize + 1)/2
       s(2) = (size(data,2) - jsize + 1)/2
       n(1) = isize
       n(2) = jsize
    CASE (2)
       s(1) = (size(data,1) - isize + 1)/2
       s(2) = (size(data,2) - ksize + 1)/2
       n(1) = isize
       n(2) = ksize
    CASE (3)
       s(1) = (size(data,1) - jsize + 1)/2
       s(2) = (size(data,2) - ksize + 1)/2
       n(1) = jsize
       n(2) = ksize
    END SELECT

    max = maxval(data(s(1)+1:s(1)+n(1), s(2)+1:s(2)+n(2)))
    min = minval(data(s(1)+1:s(1)+n(1), s(2)+1:s(2)+n(2)))

    CALL gmax(max)
    CALL gmin(min)

    IF (rank==0) WRITE(REPORT_UNIT, '(A,ES10.3,X,ES10.3)') "input '"//trim(varname)//"':  ", max, min

  END SUBROUTINE checkin_report_2d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE checkin_report_3d(varname, data)
    CHARACTER(*), INTENT(IN) :: varname
    REAL(8),      INTENT(IN) :: data(:,:,:)

    REAL(8) :: max, min

    INTEGER :: s(3), n(3)

    s(1) = (size(data,1) - isize + 1)/2
    s(2) = (size(data,2) - jsize + 1)/2
    s(3) = (size(data,3) - ksize + 1)/2
    n(1) = isize
    n(2) = jsize
    n(3) = ksize

    max = maxval(data(s(1)+1:s(1)+n(1), s(2)+1:s(2)+n(2), s(3)+1:s(3)+n(3)))
    min = minval(data(s(1)+1:s(1)+n(1), s(2)+1:s(2)+n(2), s(3)+1:s(3)+n(3)))

    CALL gmax(max)
    CALL gmin(min)

    IF (rank==0) WRITE(REPORT_UNIT, '(A,ES10.3,X,ES10.3)') "input "//trim(varname)//"':  ", max, min

  END SUBROUTINE checkin_report_3d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE checkout_3d_r4(varname, data)
    CHARACTER(*), INTENT(IN) :: varname
    REAL(4),      INTENT(IN) :: data(:,:,:)
#ifdef F2008
    CONTIGUOUS data
#endif

    REAL(8) :: tmp(size(data,1),size(data,2),size(data,3))

    IF (.NOT. require_checkout(varname)) RETURN

!$OMP PARALLEL WORKSHARE
    tmp(:,:,:) = REAL(data(:,:,:), KIND=8)
!$OMP END PARALLEL WORKSHARE

    CALL checkout_3d(varname, tmp)

  END SUBROUTINE checkout_3d_r4

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE checkout_3d_i1(varname, data)
    CHARACTER(*), INTENT(IN) :: varname
    INTEGER(1)  , INTENT(IN) :: data(:,:,:)
#ifdef F2008
    CONTIGUOUS data
#endif

    REAL(8) :: tmp(size(data,1),size(data,2),size(data,3))

    IF (.NOT. require_checkout(varname)) RETURN

!$OMP PARALLEL WORKSHARE
    tmp(:,:,:) = REAL(data(:,:,:), KIND=8)
!$OMP END PARALLEL WORKSHARE

    CALL checkout_3d(varname, tmp)

  END SUBROUTINE checkout_3d_i1

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE checkout_3d(varname, data)
    CHARACTER(*), INTENT(IN) :: varname
    REAL(8),      INTENT(IN) :: data(:,:,:)
#ifdef F2008
    CONTIGUOUS data
#endif

    INTEGER :: id

    INTEGER :: nx, ny, nz
    INTEGER :: i0, j0, k0

    LOGICAL :: gridx, gridy, gridz

    INTEGER :: start, end


    INTEGER :: i, j, k
    INTEGER :: m

    id = lookup_output_id(varname)
    IF (id == 0) RETURN

    nx = size(data,1) - isize
    ny = size(data,2) - jsize
    nz = size(data,3) - ksize

    CALL assert(nx >= 0 .AND. nx <= 5, "unsupported dimension in CHECKOUT_3D for '"//trim(varname)//"'")
    CALL assert(ny >= 0 .AND. ny <= 5, "unsupported dimension in CHECKOUT_3D for '"//trim(varname)//"'")
    CALL assert(nz >= 0 .AND. nz <= 5, "unsupported dimension in CHECKOUT_3D for '"//trim(varname)//"'")

    i0 = (nx+1)/2
    j0 = (ny+1)/2
    k0 = (nz+1)/2

    gridx = (mod(nx,2)==0)
    gridy = (mod(ny,2)==0)
    gridz = (mod(nz,2)==0)

    DO WHILE (id /= 0)
       CALL assert(.NOT. output_registry(id)%checked, "'"//trim(output_registry(id)%varname)//"' has been checked-out multiple time")

       IF (require_checkout_by_id(id)) THEN
          IF (.NOT. output_registry(id)%initialized) THEN
             SELECT CASE (output_registry(id)%mode)
             CASE (0)
                ALLOCATE(output_registry(id)%buf(isize,jsize,ksize))
             CASE (1,6)
                ALLOCATE(output_registry(id)%buf(jsize,ksize,1))
             CASE (2,7)
                ALLOCATE(output_registry(id)%buf(isize,ksize,1))
             CASE (3,4,5,8)
                ALLOCATE(output_registry(id)%buf(isize,jsize,1))
             END SELECT

             output_registry(id)%buf(:,:,:) = 0.0D0
             output_registry(id)%meancount  = 0
             output_registry(id)%dimcode    = 0
             output_registry(id)%initialized  = .TRUE.
          END IF

          CALL assert(output_registry(id)%dimcode==0, "DIMCODE comflict in CHECKOUT_3D")

          IF (output_registry(id)%mean) THEN
             output_registry(id)%meancount = output_registry(id)%meancount + 1
             m = 1
          ELSE
             m = 0
          END IF

!$OMP PARALLEL PRIVATE(i,j,k)
          SELECT CASE (output_registry(id)%mode)
          CASE (0) ! full-3D
!$OMP DO
             DO k=1, ksize
             DO j=1, jsize
             DO i=1, isize
                output_registry(id)%buf(i,j,k) = output_registry(id)%buf(i,j,k)*m + data(i0+i,j0+j,k0+k)
             END DO
             END DO
             END DO

          CASE (1) ! slice_x
             i = output_registry(id)%slice - isize*icoord
             IF ((i > 0 .AND. i  <= isize) .OR. (i == 0 .AND. icoord == 0 .AND. i0 > 0)) THEN
!$OMP DO
                DO k=1, ksize
                DO j=1, jsize
                   output_registry(id)%buf(j,k,1) = output_registry(id)%buf(j,k,1)*m + data(i0+i,j0+j,k0+k)
                END DO
                END DO
             ELSE
!$OMP WORKSHARE
                output_registry(id)%buf(:,:,1) = 0.0D0
!$OMP END WORKSHARE
             END IF

          CASE (2) ! slice_y
             j = output_registry(id)%slice - jsize*jcoord
             IF ((j > 0 .AND. j  <= jsize) .OR. (j == 0 .AND. jcoord == 0 .AND. j0 > 0)) THEN
!$OMP DO
                DO k=1, ksize
                DO i=1, isize
                   output_registry(id)%buf(i,k,1) = output_registry(id)%buf(i,k,1)*m + data(i0+i,j0+j,k0+k)
                END DO
                END DO
             ELSE
!$OMP WORKSHARE
                output_registry(id)%buf(:,:,1) = 0.0D0
!$OMP END WORKSHARE
             END IF

          CASE (3) ! slice_z
             k = output_registry(id)%slice - ksize*kcoord
             IF ((k > 0 .AND. k  <= ksize) .OR. (k == 0 .AND. kcoord == 0 .AND. k0 > 0)) THEN
!$OMP DO
                DO j=1, jsize
                DO i=1, isize
                   output_registry(id)%buf(i,j,1) = output_registry(id)%buf(i,j,1)*m + data(i0+i,j0+j,k0+k)
                END DO
                END DO
             ELSE
!$OMP WORKSHARE
                output_registry(id)%buf(:,:,1) = 0.0D0
!$OMP END WORKSHARE
             END IF

          CASE (4) ! surface

             IF (gridx .AND. gridy .AND. gridz) THEN
!$OMP DO
                DO j=1, jsize
                DO i=1, isize
                   k = surface_k(i,j)
                   IF (k >= 1 .AND. k <= ksize) THEN
                      output_registry(id)%buf(i,j,1) = output_registry(id)%buf(i,j,1)*m + data(i0+i,j0+j,k0+k)
                   ELSE
                      output_registry(id)%buf(i,j,1) = 0.0D0
                   END IF
                END DO
                END DO

             ELSE IF (.NOT. gridx .AND. gridy .AND. gridz) THEN ! x-surface
!$OMP DO
                DO j=1, jsize
                DO i=1, isize
                   k = min(surface_k(i,j), surface_k(i+1,j))
                   IF (k >= 1 .AND. k <= ksize) THEN
                      output_registry(id)%buf(i,j,1) = output_registry(id)%buf(i,j,1)*m + data(i0+i,j0+j,k0+k)
                   ELSE
                      output_registry(id)%buf(i,j,1) = 0.0D0
                   END IF
                END DO
                END DO


             ELSE IF (.NOT. gridy .AND. gridx .AND. gridz) THEN ! y-surfae
!$OMP DO
                DO j=1, jsize
                DO i=1, isize
                   k = min(surface_k(i,j), surface_k(i,j+1))
                   IF (k >= 1 .AND. k <= ksize) THEN
                      output_registry(id)%buf(i,j,1) = output_registry(id)%buf(i,j,1)*m + data(i0+i,j0+j,k0+k)
                   ELSE
                      output_registry(id)%buf(i,j,1) = 0.0D0
                   END IF
                END DO
                END DO

             ELSE IF (.NOT. gridz .AND. gridx .AND. gridy) THEN ! z-surface
!$OMP DO
                DO j=1, jsize
                DO i=1, isize
                   k = surface_k(i,j)
                   IF (k >= 1 .AND. k <= ksize) THEN
                      output_registry(id)%buf(i,j,1) = output_registry(id)%buf(i,j,1)*m + data(i0+i,j0+j,k0+k)
                   ELSE
                      output_registry(id)%buf(i,j,1) = 0.0D0
                   END IF
                END DO
                END DO

             ELSE
                CALL assert(.FALSE., "surface output is supported only for cell-center or cell-surface value")
             END IF

          CASE (5) ! bottom
             IF (gridx .AND. gridy .AND. gridz) THEN
!$OMP DO
                DO j=1, jsize
                DO i=1, isize
                   k = bottom_k(i,j)
                   IF (k >= 1 .AND. k <= ksize) THEN
                      output_registry(id)%buf(i,j,1) = output_registry(id)%buf(i,j,1)*m + data(i0+i,j0+j,k0+k)
                   ELSE
                      output_registry(id)%buf(i,j,1) = 0.0D0
                   END IF
                END DO
                END DO

             ELSE IF (.NOT. gridx .AND. gridy .AND. gridz) THEN ! x-surface
!$OMP DO
                DO j=1, jsize
                DO i=1, isize
                   k = max(bottom_k(i,j), bottom_k(i+1,j))
                   IF (k >= 1 .AND. k <= ksize) THEN
                      output_registry(id)%buf(i,j,1) = output_registry(id)%buf(i,j,1)*m + data(i0+i,j0+j,k0+k)
                   ELSE
                      output_registry(id)%buf(i,j,1) = 0.0D0
                   END IF
                END DO
                END DO


             ELSE IF (.NOT. gridy .AND. gridx .AND. gridz) THEN ! y-surfae
!$OMP DO
                DO j=1, jsize
                DO i=1, isize
                   k = max(bottom_k(i,j), bottom_k(i,j+1))
                   IF (k >= 1 .AND. k <= ksize) THEN
                      output_registry(id)%buf(i,j,1) = output_registry(id)%buf(i,j,1)*m + data(i0+i,j0+j,k0+k)
                   ELSE
                      output_registry(id)%buf(i,j,1) = 0.0D0
                   END IF
                END DO
                END DO

             ELSE IF (.NOT. gridz .AND. gridx .AND. gridy) THEN ! z-surface
!$OMP DO
                DO j=1, jsize
                DO i=1, isize
                   k = bottom_k(i,j)-1
                   IF (k >= 0 .AND. k <= ksize-1) THEN
                      output_registry(id)%buf(i,j,1) = output_registry(id)%buf(i,j,1)*m + data(i0+i,j0+j,k0+k)
                   ELSE
                      output_registry(id)%buf(i,j,1) = 0.0D0
                   END IF
                END DO
                END DO

             ELSE
                CALL assert(.FALSE., "bottom output is supported only for cell-center or cell-surface value")
             END IF

          CASE (6) ! int_x
             start = max(output_registry(id)%intspan(0) - isize*icoord, 1)
             end   = min(output_registry(id)%intspan(1) - isize*icoord, isize)
!$OMP WORKSHARE
             output_registry(id)%buf(:,:,1) = output_registry(id)%buf(:,:,1)*m
!$OMP END WORKSHARE

             IF (gridx .AND. gridy .AND. gridz) THEN ! cell-center
!$OMP DO
                DO k=1, ksize
                DO j=1, jsize
                   DO i=start, end
                      IF (.NOT. lmask3d(i,j,k)) CYCLE
                      output_registry(id)%buf(j,k,1) = output_registry(id)%buf(j,k,1) &
                           + data(i0+i,j0+j,k0+k)*dx(i,j)
                   END DO
                END DO
                END DO
             ELSE IF (.NOT. gridx .AND. gridy .AND. gridz) THEN ! x-surface
!$OMP DO
                DO k=1, ksize
                DO j=1, jsize
                   DO i=start, end
                      IF (.NOT. (lmask3d(i,j,k) .OR. lmask3d(i+1,j,k))) CYCLE
                      output_registry(id)%buf(j,k,1) = output_registry(id)%buf(j,k,1) &
                           + data(i0+i,j0+j,k0+k) * 0.5D0*(dx(i,j)+dx(i+1,j))
                   END DO
                END DO
                END DO
             ELSE IF (.NOT. gridy .AND. gridx .AND. gridz) THEN ! y-surfae
!$OMP DO
                DO k=1, ksize
                DO j=1, jsize
                   DO i=start, end
                      IF (.NOT. (lmask3d(i,j,k) .OR. lmask3d(i,j+1,k))) CYCLE
                      output_registry(id)%buf(j,k,1) = output_registry(id)%buf(j,k,1) &
                           + data(i0+i,j0+j,k0+k) * 0.5D0*(dx(i,j)+dx(i,j+1))
                   END DO
                END DO
                END DO
             ELSE IF (.NOT. gridz .AND. gridx .AND. gridy) THEN ! z-surfae
!$OMP DO
                DO k=1, ksize
                DO j=1, jsize
                   DO i=start, end
                      IF (.NOT. (lmask3d(i,j,k) .OR. lmask3d(i,j,k+1))) CYCLE
                      output_registry(id)%buf(j,k,1) = output_registry(id)%buf(j,k,1) &
                           + data(i0+i,j0+j,k0+k) * dx(i,j)
                   END DO
                END DO
                END DO
             ELSE
                CALL assert(.FALSE., "INT_X output is supported only for cell-center or cell-surface value")
             END IF

          CASE (7) ! int_y
             start = max(output_registry(id)%intspan(0) - jsize*jcoord, 1)
             end   = min(output_registry(id)%intspan(1) - jsize*jcoord, jsize)
!$OMP WORKSHARE
             output_registry(id)%buf(:,:,1) = output_registry(id)%buf(:,:,1)*m
!$OMP END WORKSHARE

             IF (gridx .AND. gridy .AND. gridz) THEN ! cell-center
!$OMP DO
                DO k=1, ksize
                DO i=1, isize
                   DO j=start, end
                      IF (.NOT. lmask3d(i,j,k)) CYCLE
                      output_registry(id)%buf(i,k,1) = output_registry(id)%buf(i,k,1) &
                           + data(i0+i,j0+j,k0+k)*dy(i,j)
                   END DO
                END DO
                END DO
             ELSE IF (.NOT. gridx .AND. gridy .AND. gridz) THEN ! x-surface
!$OMP DO
                DO k=1, ksize
                DO i=1, isize
                   DO j=start, end
                      IF (.NOT. (lmask3d(i,j,k) .OR. lmask3d(i+1,j,k))) CYCLE
                      output_registry(id)%buf(i,k,1) = output_registry(id)%buf(i,k,1) &
                           + data(i0+i,j0+j,k0+k) * 0.5D0*(dy(i,j)+dy(i+1,j))
                   END DO
                END DO
                END DO
             ELSE IF (.NOT. gridy .AND. gridx .AND. gridz) THEN ! y-surfae
!$OMP DO
                DO k=1, ksize
                DO i=1, isize
                   DO j=start, end
                      IF (.NOT. (lmask3d(i,j,k) .OR. lmask3d(i,j+1,k))) CYCLE
                      output_registry(id)%buf(i,k,1) = output_registry(id)%buf(i,k,1) &
                           + data(i0+i,j0+j,k0+k) * 0.5D0*(dy(i,j)+dy(i,j+1))
                   END DO
                END DO
                END DO
             ELSE IF (.NOT. gridz .AND. gridx .AND. gridy) THEN ! z-surfae
!$OMP DO
                DO k=1, ksize
                DO i=1, isize
                   DO j=start, end
                      IF (.NOT. (lmask3d(i,j,k) .OR. lmask3d(i,j,k+1))) CYCLE
                      output_registry(id)%buf(i,k,1) = output_registry(id)%buf(i,k,1) &
                           + data(i0+i,j0+j,k0+k) * dy(i,j)
                   END DO
                END DO
                END DO
             ELSE
                CALL assert(.FALSE., "INT_Y output is supported only for cell-center or cell-surface value")
             END IF

          CASE (8) ! int_z
             start = max(output_registry(id)%intspan(0) - ksize*kcoord, 1)
             end   = min(output_registry(id)%intspan(1) - ksize*kcoord, ksize)

!$OMP WORKSHARE
             output_registry(id)%buf(:,:,1) = output_registry(id)%buf(:,:,1)*m
!$OMP END WORKSHARE

             IF (gridx .AND. gridy .AND. gridz) THEN
                DO k=start, end
!$OMP DO
                   DO j=1, jsize
                   DO i=1, isize
                      IF (.NOT. lmask3d(i,j,k)) CYCLE
                      output_registry(id)%buf(i,j,1) = output_registry(id)%buf(i,j,1) &
                           + data(i0+i,j0+j,k0+k)*dz_ref(i,j,k)
                   END DO
                   END DO
                END DO
             ELSE IF (.NOT. gridx .AND. gridy .AND. gridz) THEN ! x-surface
                DO k=start, end
!$OMP DO
                   DO j=1, jsize
                   DO i=1, isize
                      IF (.NOT. (lmask3d(i,j,k) .OR. lmask3d(i+1,j,k))) CYCLE
                      output_registry(id)%buf(i,j,1) = output_registry(id)%buf(i,j,1) &
                           + data(i0+i,j0+j,k0+k) * (dvol(i,j,k)+dvol(i+1,j,k))/(dsz(i,j)+dsz(i+1,j))
                   END DO
                   END DO
                END DO
             ELSE IF (.NOT. gridy .AND. gridx .AND. gridz) THEN ! y-surfae
                DO k=start, end
!$OMP DO
                   DO j=1, jsize
                   DO i=1, isize
                      IF (.NOT. (lmask3d(i,j,k) .OR. lmask3d(i,j+1,k))) CYCLE
                      output_registry(id)%buf(i,j,1) = output_registry(id)%buf(i,j,1) &
                           + data(i0+i,j0+j,k0+k) * (dvol(i,j,k)+dvol(i,j+1,k))/(dsz(i,j)+dsz(i,j+1))
                   END DO
                   END DO
                END DO
             ELSE
                CALL assert(.FALSE., "INT_Z output is supported only for cell-center or cell-surface value")
             END IF

          END SELECT
!$OMP END PARALLEL

          output_registry(id)%checked = .TRUE.
       END IF

       id = lookup_output_id(varname, skip=id)
    END DO

  END SUBROUTINE checkout_3d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE checkout_2d_r4(varname, data, section)
    CHARACTER(*), INTENT(IN) :: varname
    REAL(4),      INTENT(IN) :: data(:,:)
    CHARACTER(2), INTENT(IN), OPTIONAL :: section
#ifdef F2008
    CONTIGUOUS data
#endif

    REAL(8) :: tmp(size(data,1),size(data,2))

    IF (.NOT. require_checkout(varname)) RETURN

!$OMP PARALLEL WORKSHARE
    tmp(:,:) = REAL(data(:,:), KIND=8)
!$OMP END PARALLEL WORKSHARE

    CALL checkout_2d(varname, tmp, section)

  END SUBROUTINE checkout_2d_r4

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE checkout_2d_i1(varname, data, section)
    CHARACTER(*), INTENT(IN) :: varname
    INTEGER(1),   INTENT(IN) :: data(:,:)
    CHARACTER(2), INTENT(IN), OPTIONAL :: section
#ifdef F2008
    CONTIGUOUS data
#endif

    REAL(8) :: tmp(size(data,1),size(data,2))

    IF (.NOT. require_checkout(varname)) RETURN

!$OMP PARALLEL WORKSHARE
    tmp(:,:) = REAL(data(:,:), KIND=8)
!$OMP END PARALLEL WORKSHARE

    CALL checkout_2d(varname, tmp, section)

  END SUBROUTINE checkout_2d_i1

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE checkout_2d(varname, data, section)
    CHARACTER(*), INTENT(IN) :: varname
    REAL(8),      INTENT(IN) :: data(:,:)
    CHARACTER(2), INTENT(IN), OPTIONAL :: section

    INTEGER :: id

    INTEGER :: nx, ny
    INTEGER :: i0, j0

    INTEGER :: i, j
    INTEGER :: m

    id = lookup_output_id(varname)
    IF (id == 0) RETURN

    IF (present(section)) THEN
       SELECT CASE(section)
       CASE ('XZ', 'xz')
          CALL checkout_xz(varname, data)
          RETURN
       CASE ('YZ', 'yz')
          CALL checkout_yz(varname, data)
          RETURN
       CASE ('XY', 'xy')
       CASE DEFAULT
          CALL assert(.FALSE., "invalid SECTION '"//section//"' in CHECKOUT_2D")
       END SELECT
    END IF

    nx = size(data,1) - isize
    ny = size(data,2) - jsize

    CALL assert(nx >= 0 .AND. nx <= 5, "unsupported dimension in CHECKOUT_2D for '"//trim(varname)//"'")
    CALL assert(ny >= 0 .AND. ny <= 5, "unsupported dimension in CHECKOUT_2D for '"//trim(varname)//"'")

    i0 = (nx+1)/2
    j0 = (ny+1)/2

    DO WHILE (id /= 0)
       CALL assert(.NOT. output_registry(id)%checked, "'"//trim(output_registry(id)%varname)//"' has been checked-out multiple time")
       CALL assert(output_registry(id)%mode==0, "SLICE_X/Y/Z, INT_X/Y/Z and SURFACE/BOTTOM flags are not allowd for output of 2D-variable" &
                                                 // "'"//trim(output_registry(id)%varname)//"'")

       IF (require_checkout_by_id(id)) THEN
          IF (.NOT. output_registry(id)%initialized) THEN
             ALLOCATE(output_registry(id)%buf(isize,jsize,1))

             output_registry(id)%buf(:,:,:) = 0.0D0
             output_registry(id)%meancount  = 0
             output_registry(id)%dimcode    = 1
             output_registry(id)%initialized  = .TRUE.
          END IF

          CALL assert(output_registry(id)%dimcode==1, "DIMCODE conflict in CHECKOUT_2D")

          IF (output_registry(id)%mean) THEN
             output_registry(id)%meancount = output_registry(id)%meancount + 1
             m = 1
          ELSE
             m = 0
          END IF

!$OMP PARALLEL DO
          DO j=1, jsize
          DO i=1, isize
             output_registry(id)%buf(i,j,1) = output_registry(id)%buf(i,j,1)*m + data(i0+i,j0+j)
          END DO
          END DO

          output_registry(id)%checked = .TRUE.
       END IF

       id = lookup_output_id(varname, skip=id)
    END DO

  END SUBROUTINE checkout_2d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE checkout_xz(varname, data)
    CHARACTER(*), INTENT(IN) :: varname
    REAL(8), INTENT(IN) :: data(:,:)
#ifdef F2008
    CONTIGUOUS data
#endif

    INTEGER :: id

    INTEGER :: nx, nz
    INTEGER :: i0, k0

    INTEGER :: i, k
    INTEGER :: m

    id = lookup_output_id(varname)
    IF (id == 0) RETURN

    nx = size(data,1) - isize
    nz = size(data,2) - ksize

    CALL assert(nx >= 0 .AND. nx <= 5, "unsupported dimension in CHECKOUT_XZ for '"//trim(varname)//"'")
    CALL assert(nz >= 0 .AND. nz <= 5, "unsupported dimension in CHECKOUT_XZ for '"//trim(varname)//"'")

    i0 = (nx+1)/2
    k0 = (nz+1)/2

    DO WHILE (id /= 0)
       CALL assert(.NOT. output_registry(id)%checked, "'"//trim(output_registry(id)%varname)//"' has been checked-out multiple time")
       CALL assert(output_registry(id)%mode==0, "SLICE_X/Y/Z, INT_X/Y/Z and SURFACE/BOTTOM flags are not allowd for output of 2D-variable" &
                                                 // "'"//trim(output_registry(id)%varname)//"'")

       IF (require_checkout_by_id(id)) THEN
          IF (.NOT. output_registry(id)%initialized) THEN
             ALLOCATE(output_registry(id)%buf(isize,ksize,1))

             output_registry(id)%buf(:,:,:) = 0.0D0
             output_registry(id)%meancount  = 0
             output_registry(id)%dimcode    = 2
             output_registry(id)%initialized  = .TRUE.
          END IF

          CALL assert(output_registry(id)%dimcode==2, "DIMCODE conflict in CHECKOUT_XZ")

          IF (output_registry(id)%mean) THEN
             output_registry(id)%meancount = output_registry(id)%meancount + 1
             m = 1
          ELSE
             m = 0
          END IF

!$OMP PARALLEL DO
          DO k=1, ksize
          DO i=1, isize
             output_registry(id)%buf(i,k,1) = output_registry(id)%buf(i,k,1)*m + data(i0+i,k0+k)
          END DO
          END DO

          output_registry(id)%checked = .TRUE.
       END IF

       id = lookup_output_id(varname, skip=id)
    END DO

  END SUBROUTINE checkout_xz

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE checkout_yz(varname, data)
    CHARACTER(*), INTENT(IN) :: varname
    REAL(8), INTENT(IN) :: data(:,:)
#ifdef F2008
    CONTIGUOUS data
#endif

    INTEGER :: id

    INTEGER :: ny, nz
    INTEGER :: j0, k0

    INTEGER :: j, k
    INTEGER :: m

    id = lookup_output_id(varname)
    IF (id == 0) RETURN

    ny = size(data,1) - jsize
    nz = size(data,2) - ksize

    CALL assert(ny >= 0 .AND. ny <= 5, "unsupported dimension in CHECKOUT_YZ for '"//trim(varname)//"'")
    CALL assert(nz >= 0 .AND. nz <= 5, "unsupported dimension in CHECKOUT_YZ for '"//trim(varname)//"'")

    j0 = (ny+1)/2
    k0 = (nz+1)/2

    DO WHILE (id /= 0)
       CALL assert(.NOT. output_registry(id)%checked, "'"//trim(output_registry(id)%varname)//"' has been checked-out multiple time")
       CALL assert(output_registry(id)%mode==0, "SLICE_X/Y/Z, INT_X/Y/Z and SURFACE/BOTTOM flags are not allowd for output of 2D-variable" &
                                                 // "'"//trim(output_registry(id)%varname)//"'")

       IF (require_checkout_by_id(id)) THEN
          IF (.NOT. output_registry(id)%initialized) THEN
             ALLOCATE(output_registry(id)%buf(jsize,ksize,1))

             output_registry(id)%buf(:,:,:) = 0.0D0
             output_registry(id)%meancount  = 0
             output_registry(id)%dimcode    = 3
             output_registry(id)%initialized  = .TRUE.
          END IF

          CALL assert(output_registry(id)%dimcode==3, "DIMCODE conflict in CHECKOUT_YZ")

          IF (output_registry(id)%mean) THEN
             output_registry(id)%meancount = output_registry(id)%meancount + 1
             m = 1
          ELSE
             m = 0
          END IF

!$OMP PARALLEL DO
          DO k=1, ksize
          DO j=1, jsize
             output_registry(id)%buf(j,k,1) = output_registry(id)%buf(j,k,1)*m + data(j0+j,k0+k)
          END DO
          END DO

          output_registry(id)%checked = .TRUE.
       END IF

       id = lookup_output_id(varname, skip=id)
    END DO

  END SUBROUTINE checkout_yz

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE checkout_1d_r4(varname, data, axis)
    CHARACTER(*), INTENT(IN) :: varname
    REAL(4),      INTENT(IN) :: data(:)
    CHARACTER(1), INTENT(IN) :: axis
#ifdef F2008
    CONTIGUOUS data
#endif

    REAL(8) :: tmp(size(data))

    IF (.NOT. require_checkout(varname)) RETURN

    tmp(:) = REAL(data(:), KIND=8)

    CALL checkout_1d(varname, tmp, axis)

  END SUBROUTINE checkout_1d_r4

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE checkout_1d_i1(varname, data, axis)
    CHARACTER(*), INTENT(IN) :: varname
    INTEGER(1),   INTENT(IN) :: data(:)
    CHARACTER(1), INTENT(IN) :: axis
#ifdef F2008
    CONTIGUOUS data
#endif

    REAL(8) :: tmp(size(data))

    IF (.NOT. require_checkout(varname)) RETURN

    tmp(:) = REAL(data(:), KIND=8)

    CALL checkout_1d(varname, tmp, axis)

  END SUBROUTINE checkout_1d_i1

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE checkout_1d(varname, data, axis)
    CHARACTER(*), INTENT(IN) :: varname
    REAL(8),      INTENT(IN) :: data(:)
    CHARACTER(1), INTENT(IN) :: axis
#ifdef F2008
    CONTIGUOUS data
#endif

    INTEGER :: id

    INTEGER :: n
    INTEGER :: slv

    INTEGER :: i
    INTEGER :: m

    INTEGER :: size_
    INTEGER :: dimcode

    id = lookup_output_id(varname)
    IF (id == 0) RETURN

    SELECT CASE(axis)
    CASE ('X', 'x')
       dimcode = -1
       size_ = isize
    CASE ('Y', 'y')
       dimcode = -2
       size_ = jsize
    CASE ('Z', 'z')
       dimcode = -3
       size_ = ksize
    CASE DEFAULT
       CALL assert(.FALSE., "invalid AXIS '"//axis//"' in CHECKOUT_1D")
    END SELECT

    n = size(data) - size_

    CALL assert(n >= 0 .AND. n <= 5, "unsupported dimension in CHECKOUT_1D for '"//trim(varname)//"'")

    slv = (n+1)/2

    DO WHILE (id /= 0)
       CALL assert(.NOT. output_registry(id)%checked, "'"//trim(output_registry(id)%varname)//"' has been checked-out multiple time")
       CALL assert(output_registry(id)%mode==0, "SLICE_X/Y/Z, INT_X/Y/Z and SURFACE/BOTTOM flags are not allowd for output of 1D-variable " &
                                                 // "'"//trim(output_registry(id)%varname)//"'")

       IF (require_checkout_by_id(id)) THEN
          IF (.NOT. output_registry(id)%initialized) THEN
             ALLOCATE(output_registry(id)%buf(size_,1,1))

             output_registry(id)%buf(:,:,:) = 0.0D0
             output_registry(id)%meancount  = 0
             output_registry(id)%dimcode    = dimcode
             output_registry(id)%initialized  = .TRUE.
          END IF

          CALL assert(output_registry(id)%dimcode==dimcode, "DIMCODE conflict in CHECKOUT_1D")

          IF (output_registry(id)%mean) THEN
             output_registry(id)%meancount = output_registry(id)%meancount + 1
             m = 1
          ELSE
             m = 0
          END IF

          DO i=1, size_
             output_registry(id)%buf(i,1,1) = output_registry(id)%buf(i,1,1)*m + data(slv+i)
          END DO

          output_registry(id)%checked = .TRUE.
       END IF

       id = lookup_output_id(varname, skip=id)
    END DO

  END SUBROUTINE checkout_1d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE checkout_zeros(varname, axis, section)
    CHARACTER(*), INTENT(IN) :: varname
    CHARACTER(1), INTENT(IN), OPTIONAL :: axis
    CHARACTER(2), INTENT(IN), OPTIONAL :: section
    INTEGER :: id

    INTEGER :: dimcode
    INTEGER :: size_(3)

    id = lookup_output_id(varname)
    IF (id == 0) RETURN

    CALL assert(.NOT. (present(axis) .AND. present(section)), "argument AXIS and SECTION is exclusive")

    IF (present(axis)) THEN
       SELECT CASE(axis)
       CASE ('X', 'x')
          dimcode = -1
          size_ = (/isize,1,1/)
       CASE ('Y', 'y')
          dimcode = -2
          size_ = (/jsize,1,1/)
       CASE ('Z', 'z')
          dimcode = -3
          size_ = (/ksize,1,1/)
       CASE DEFAULT
          CALL assert(.FALSE., "invalid AXIS '"//axis//"' in CHECKOUT_ZEROS")
       END SELECT
    ELSE IF (present(section)) THEN
       SELECT CASE(section)
       CASE ('XY', 'xy')
          dimcode = 1
          size_ = (/isize,jsize,1/)
       CASE ('XZ', 'xz')
          dimcode = 2
          size_ = (/isize,ksize,1/)
       CASE ('YZ', 'yz')
          dimcode = 3
          size_ = (/jsize,ksize,1/)
       CASE DEFAULT
          CALL assert(.FALSE., "invalid SECTION '"//section//"' in CHECKOUT_ZEROS")
       END SELECT
    ELSE
       dimcode = 0
       size_ = (/isize,jsize,ksize/)
    END IF

    DO WHILE (id /= 0)
       CALL assert(.NOT. output_registry(id)%checked, "'"//trim(output_registry(id)%varname)//"' has been checked-out multiple time")

       IF (require_checkout_by_id(id)) THEN
          IF (.NOT. output_registry(id)%initialized) THEN
             ALLOCATE(output_registry(id)%buf(size_(1),size_(2),size_(3)))

             output_registry(id)%buf(:,:,:) = 0.0D0
             output_registry(id)%meancount  = 0
             output_registry(id)%dimcode    = dimcode
             output_registry(id)%initialized  = .TRUE.
          END IF

          CALL assert(output_registry(id)%dimcode==dimcode, "DIMCODE conflict in CHECKOUT_ZEROS")

          IF (output_registry(id)%mean) THEN
             output_registry(id)%meancount = output_registry(id)%meancount + 1
          END IF

          output_registry(id)%checked = .TRUE.
       END IF

       id = lookup_output_id(varname, skip=id)
    END DO

  END SUBROUTINE checkout_zeros

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE checkout_div(varname, u, v, w)
    CHARACTER(*), INTENT(IN) :: varname
    REAL(8), INTENT(IN) :: u(:,:,:)
    REAL(8), INTENT(IN) :: v(:,:,:)
    REAL(8), INTENT(IN) :: w(:,:,:)
#ifdef F2008
    CONTIGUOUS u, v, w
#endif

    REAL(8) :: tmp(isize,jsize,ksize)

    INTEGER :: u_nx, u_ny, u_nz
    INTEGER :: v_nx, v_ny, v_nz
    INTEGER :: w_nx, w_ny, w_nz

    INTEGER :: u_i0, u_j0, u_k0
    INTEGER :: v_i0, v_j0, v_k0
    INTEGER :: w_i0, w_j0, w_k0

    INTEGER :: i, j, k

    u_nx = size(u,1) - isize
    u_ny = size(u,2) - jsize
    u_nz = size(u,3) - ksize
    CALL assert(mod(u_nx,2)==1 .AND. mod(u_ny,2)==0 .AND. mod(u_nz,2)==0, "the 2nd argument in CHECKOUT_DIV_3D is not x-component of a vector")
    u_i0 = (u_nx+1)/2
    u_j0 = (u_ny+1)/2
    u_k0 = (u_nz+1)/2

    v_nx = size(v,1) - isize
    v_ny = size(v,2) - jsize
    v_nz = size(v,3) - ksize
    CALL assert(mod(v_nx,2)==0 .AND. mod(v_ny,2)==1 .AND. mod(v_nz,2)==0, "the 3rd argument in CHECKOUT_DIV_3D is not y-component of a vector")
    v_i0 = (v_nx+1)/2
    v_j0 = (v_ny+1)/2
    v_k0 = (v_nz+1)/2

    w_nx = size(w,1) - isize
    w_ny = size(w,2) - jsize
    w_nz = size(w,3) - ksize
    CALL assert(mod(w_nx,2)==0 .AND. mod(w_ny,2)==0 .AND. mod(w_nz,2)==1, "the 4th argument in CHECKOUT_DIV_3D is not z-component of a vector")
    w_i0 = (w_nx+1)/2
    w_j0 = (w_ny+1)/2
    w_k0 = (w_nz+1)/2

    IF (require_checkout('DIV_'//trim(varname))) THEN
!$OMP PARALLEL DO
       DO k=1, ksize
       DO j=1, jsize
       DO i=1, isize
          tmp(i,j,k) = imask3d(i,j,k)*(  u(u_i0+i,  u_j0+j,  u_k0+k  )*dsx_old(i,  j,  k) &
                                       - u(u_i0+i-1,u_j0+j,  u_k0+k  )*dsx_old(i-1,j,  k) &
                                       + v(v_i0+i,  v_j0+j,  v_k0+k  )*dsy_old(i,  j,  k) &
                                       - v(v_i0+i,  v_j0+j-1,v_k0+k  )*dsy_old(i,  j-1,k) &
                                       + w(w_i0+i,  w_j0+j,  w_k0+k  )*dsz(i,j)           &
                                       - w(w_i0+i,  w_j0+j,  w_k0+k-1)*dsz(i,j)) / dvol(i,j,k)
       END DO
       END DO
       END DO
       CALL checkout('DIV_'//trim(varname), tmp)
    END IF

    IF (require_checkout('DIVH_'//trim(varname))) THEN
!$OMP PARALLEL DO
       DO k=1, ksize
       DO j=1, jsize
       DO i=1, isize
          tmp(i,j,k) = imask3d(i,j,k)*(  u(u_i0+i,  u_j0+j,  u_k0+k)*dsx_old(i,  j,  k) &
                                       - u(u_i0+i-1,u_j0+j,  u_k0+k)*dsx_old(i-1,j,  k) &
                                       + v(v_i0+i,  v_j0+j,  v_k0+k)*dsy_old(i,  j,  k) &
                                       - v(v_i0+i,  v_j0+j-1,v_k0+k)*dsy_old(i,  j-1,k)) / dvol(i,j,k)
       END DO
       END DO
       END DO
       CALL checkout('DIVH_'//trim(varname), tmp)
    END IF

    IF (require_checkout('DIVX_'//trim(varname))) THEN
!$OMP PARALLEL DO
       DO k=1, ksize
       DO j=1, jsize
       DO i=1, isize
          tmp(i,j,k) = imask3d(i,j,k)*(  u(u_i0+i,  u_j0+j,  u_k0+k)*dsx_old(i,  j,  k) &
                                       - u(u_i0+i-1,u_j0+j,  u_k0+k)*dsx_old(i-1,j,  k)) / dvol(i,j,k)
       END DO
       END DO
       END DO
       CALL checkout('DIVX_'//trim(varname), tmp)
    END IF

    IF (require_checkout('DIVY_'//trim(varname))) THEN
!$OMP PARALLEL DO
       DO k=1, ksize
       DO j=1, jsize
       DO i=1, isize
          tmp(i,j,k) = imask3d(i,j,k)*(  v(v_i0+i,  v_j0+j,  v_k0+k)*dsy_old(i,  j,  k) &
                                       - v(v_i0+i,  v_j0+j-1,v_k0+k)*dsy_old(i,  j-1,k)) / dvol(i,j,k)
       END DO
       END DO
       END DO
       CALL checkout('DIVY_'//trim(varname), tmp)
    END IF

    IF (require_checkout('DIVZ_'//trim(varname))) THEN
!$OMP PARALLEL DO
       DO k=1, ksize
       DO j=1, jsize
       DO i=1, isize
          tmp(i,j,k) = imask3d(i,j,k)*(  w(w_i0+i,  w_j0+j,  w_k0+k  )*dsz(i,j) &
                                       - w(w_i0+i,  w_j0+j,  w_k0+k-1)*dsz(i,j)) / dvol(i,j,k)
       END DO
       END DO
       END DO
       CALL checkout('DIVZ_'//trim(varname), tmp)
    END IF

  END SUBROUTINE checkout_div

!-----------------------------------------------------------------------------------------------------------------------

  LOGICAL PURE FUNCTION require_checkout_by_name(varname, time)
    CHARACTER(*), INTENT(IN) :: varname
    REAL(8),      INTENT(IN), OPTIONAL :: time

    INTEGER :: id

    id = lookup_output_id(varname)

    DO WHILE (id /= 0)
       IF (require_checkout_by_id(id, time)) THEN
          require_checkout_by_name = .TRUE.
          RETURN
       END IF

       id = lookup_output_id(varname, skip=id)
    END DO

    require_checkout_by_name = .FALSE.

  END FUNCTION require_checkout_by_name

!-----------------------------------------------------------------------------------------------------------------------

  LOGICAL PURE FUNCTION require_checkout_by_id(id, time)
    INTEGER, INTENT(IN) :: id
    REAL(8), INTENT(IN), OPTIONAL :: time

    REAL(8) :: t

    IF (present(time)) THEN
       t = time + dtime
    ELSE
       t = t_current + dtime
    END IF

    require_checkout_by_id = .FALSE.

    IF (id < 1 .OR. id > n_output) RETURN

    IF (t <= output_registry(id)%start .OR. t > output_registry(id)%end) RETURN

    IF (output_registry(id)%mean) THEN
       require_checkout_by_id = t > output_registry(id)%lastwrite + output_registry(id)%interval - output_registry(id)%meanspan
    ELSE
       require_checkout_by_id = t >=output_registry(id)%lastwrite + output_registry(id)%interval .OR. (t <  output_registry(id)%start + dtime)
    END IF
  END FUNCTION require_checkout_by_id

!-----------------------------------------------------------------------------------------------------------------------

  LOGICAL PURE FUNCTION check_time(time, start, end, interval, last)
    REAL(8), INTENT(IN) :: time, start, end, interval, last
    check_time = (time >= start) .AND. (time <= end) .AND. (time-last >= interval .OR. time < start+dtime)
  END FUNCTION check_time

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE flush_io

    INTEGER :: n, i, j, k

    CHARACTER(1024) :: filepath

    REAL(8) :: t

    DO n=1, n_output
       output_registry(n)%checked = .FALSE.

       IF (.NOT. check_time(t_current, output_registry(n)%start,    output_registry(n)%end, &
                                       output_registry(n)%interval, output_registry(n)%lastwrite)) CYCLE

       IF (.NOT. output_registry(n)%initialized) THEN
          IF (rank==0) WRITE(REPORT_UNIT, *) "'"//trim(output_registry(n)%varname)//"' has not been checked-out"
          CYCLE
       END IF

       t = t_current - mod(t_current - output_registry(n)%start, output_registry(n)%interval)
       output_registry(n)%lastwrite = t

       filepath = path(output_registry(n)%outputdir, output_registry(n)%basename, '.'//format_datetime(t, omit_time=output_registry(n)%omit_time))
       CALL replace_datetime(filepath, t)

       IF (output_registry(n)%mean) THEN
          filepath = trim(filepath) // '_MEAN'
          CALL scal(1.0D0 / output_registry(n)%meancount, output_registry(n)%buf(:,:,:))
       END IF

       IF (output_registry(n)%dimcode==1) THEN
          CALL write_data_2d(output_registry(n)%buf(:,:,1), filepath, kind=output_registry(n)%kind, region=output_registry(n)%region)
       ELSE IF (output_registry(n)%dimcode==2) THEN
          CALL write_data_xz(output_registry(n)%buf(:,:,1), filepath, kind=output_registry(n)%kind, region=output_registry(n)%region)
       ELSE IF (output_registry(n)%dimcode==3) THEN
          CALL write_data_yz(output_registry(n)%buf(:,:,1), filepath, kind=output_registry(n)%kind, region=output_registry(n)%region)
       ELSE IF (output_registry(n)%dimcode==-1) THEN
          CALL write_data_x(output_registry(n)%buf(:,1,1), filepath, kind=output_registry(n)%kind, region=output_registry(n)%region)
       ELSE IF (output_registry(n)%dimcode==-2) THEN
          CALL write_data_y(output_registry(n)%buf(:,1,1), filepath, kind=output_registry(n)%kind, region=output_registry(n)%region)
       ELSE IF (output_registry(n)%dimcode==-3) THEN
          CALL write_data_z(output_registry(n)%buf(:,1,1), filepath, kind=output_registry(n)%kind, region=output_registry(n)%region)
       ELSE
          SELECT CASE (output_registry(n)%mode)
          CASE (0)
             CALL write_data_3d(output_registry(n)%buf(:,:,:), filepath, kind=output_registry(n)%kind, region=output_registry(n)%region)

          CASE (1,6)
             CALL write_data_yz(output_registry(n)%buf(:,:,1), filepath, kind=output_registry(n)%kind, region=output_registry(n)%region)

          CASE (2,7)
             CALL write_data_xz(output_registry(n)%buf(:,:,1), filepath, kind=output_registry(n)%kind, region=output_registry(n)%region)

          CASE (3,4,5,8)
             CALL vsum(output_registry(n)%buf(:,:,1))
             CALL write_data_2d(output_registry(n)%buf(:,:,1), filepath, kind=output_registry(n)%kind, region=output_registry(n)%region)
          END SELECT
       END IF

       IF (output_registry(n)%mean) THEN
!$OMP PARALLEL WORKSHARE
          output_registry(n)%buf(:,:,:) = 0.0D0
!$OMP END PARALLEL WORKSHARE
          output_registry(n)%meancount = 0
       END IF
    END DO

    CALL barrier

  END SUBROUTINE flush_io

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE finalize_io
    INTEGER :: n

    DO n=1, n_output
       IF (associated(output_registry(n)%buf)) DEALLOCATE(output_registry(n)%buf)
    END DO

    DO n=1, n_input
       IF (input_registry(n)%start <= t_end .AND. input_registry(n)%end > t_start) THEN
          IF (rank==0) CALL warning(input_registry(n)%initialized, &
                                    "'" // trim(input_registry(n)%varname) // "'" // ' has not been checked-in')
       END IF

       IF (associated(input_registry(n)%buf))      DEALLOCATE(input_registry(n)%buf)
       IF (associated(input_registry(n)%filename)) DEALLOCATE(input_registry(n)%filename)
    END DO

  END SUBROUTINE finalize_io

!-----------------------------------------------------------------------------------------------------------------------

  INTEGER PURE FUNCTION lookup_input_id(varname, skip)
    CHARACTER(*), INTENT(IN) :: varname
    INTEGER,      INTENT(IN), OPTIONAL :: skip

    INTEGER :: n, s, i

    lookup_input_id = 0

    IF (len_trim(varname) == 0) RETURN

    i = iachar(varname(1:1))

    s = input_registry_index(i)

    IF (present(skip)) s = max(s, skip+1)

    DO n=s, input_registry_index(i+1)-1
       IF (trim(input_registry(n)%varname) == trim(varname)) THEN
          lookup_input_id = n
          RETURN
       END IF
    END DO

  END FUNCTION lookup_input_id

!-----------------------------------------------------------------------------------------------------------------------

  INTEGER PURE FUNCTION lookup_output_id(varname, skip)
    CHARACTER(*), INTENT(IN) :: varname
    INTEGER,      INTENT(IN), OPTIONAL :: skip

    INTEGER :: n, s, i

    lookup_output_id = 0

    IF (len_trim(varname) == 0) RETURN

    i = iachar(varname(1:1))

    s = output_registry_index(i)

    IF (present(skip)) s = max(s, skip+1)

    DO n=s, output_registry_index(i+1)-1
       IF (trim(output_registry(n)%varname) == trim(varname)) THEN
          lookup_output_id = n
          RETURN
       END IF
    END DO

  END FUNCTION lookup_output_id

!-----------------------------------------------------------------------------------------------------------------------

  INTEGER PURE FUNCTION count_input(varname)
    CHARACTER(*), INTENT(IN) :: varname

    INTEGER :: n

    count_input = 0
    DO n=1, n_input
       IF (trim(input_registry(n)%varname) == trim(varname)) count_input = count_input + 1
    END DO

  END FUNCTION count_input

!-----------------------------------------------------------------------------------------------------------------------

  INTEGER PURE FUNCTION count_output(varname)
    CHARACTER(*), INTENT(IN) :: varname

    INTEGER :: n

    count_output = 0
    DO n=1, n_output
       IF (trim(output_registry(n)%varname) == trim(varname)) count_output = count_output + 1
    END DO

  END FUNCTION count_output

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE default_output(default_outputdir, default_precision, &
                            default_start, default_end, default_interval, &
                            default_omit_time)
    CHARACTER(512), INTENT(OUT), OPTIONAL :: default_outputdir
    CHARACTER(8),   INTENT(OUT), OPTIONAL :: default_precision
    CHARACTER(16),  INTENT(OUT), OPTIONAL :: default_start
    CHARACTER(16),  INTENT(OUT), OPTIONAL :: default_end
    CHARACTER(16),  INTENT(OUT), OPTIONAL :: default_interval
    LOGICAL,        INTENT(OUT), OPTIONAL :: default_omit_time

    CHARACTER(512) :: outputdir
    CHARACTER(8)   :: precision
    CHARACTER(16)  :: start, start_date
    CHARACTER(16)  :: end,   end_date
    CHARACTER(16)  :: interval
    LOGICAL        :: omit_time

    INTEGER        :: iostat
    CHARACTER(256) :: iomsg

    NAMELIST / output_default / &
         outputdir, &
         precision, &
         start,     &
         start_date,&
         end,       &
         end_date,  &
         interval,  &
         omit_time

    IF (rank==0) THEN
       outputdir = global_outputdir
       precision = 'REAL8'
       start     = format_datetime(t_start)
       end       = format_datetime(t_end)
       start_date= ''
       end_date  = ''
       interval  = '0000001_000000'
       omit_time = .FALSE.

       REWIND(CONFIG_UNIT)
       READ(CONFIG_UNIT, NML=output_default, IOSTAT=iostat, IOMSG=iomsg)

       CALL assert(iostat <= 0, "failed to read OUTPUT_DEFAULT namelist", iomsg)

       IF (start_date/= '') start = start_date
       IF (end_date  /= '') end   = end_date
    END IF

    CALL bcast(outputdir)
    CALL bcast(precision)
    CALL bcast(start)
    CALL bcast(end)
    CALL bcast(interval)
    CALL bcast(omit_time)

    IF (present(default_outputdir))    default_outputdir    = outputdir
    IF (present(default_precision))    default_precision    = precision
    IF (present(default_start))        default_start        = start
    IF (present(default_end))          default_end          = end
    IF (present(default_interval))     default_interval     = interval
    IF (present(default_omit_time))    default_omit_time  = omit_time

  END SUBROUTINE default_output

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE default_input(default_inputdir, default_precision, default_missing, default_descend, default_fileview, &
                           default_mode, default_periods,       &
                           default_start, default_end, default_interval, &
                           default_histshift, default_histcycle, &
                           default_omit_time, default_inquire)
    CHARACTER(512), INTENT(OUT), OPTIONAL :: default_inputdir
    CHARACTER(8),   INTENT(OUT), OPTIONAL :: default_precision
    REAL(4),        INTENT(OUT), OPTIONAL :: default_missing(2)
    LOGICAL,        INTENT(OUT), OPTIONAL :: default_descend
    INTEGER,        INTENT(OUT), OPTIONAL :: default_fileview
    CHARACTER(20),  INTENT(OUT), OPTIONAL :: default_mode
    INTEGER,        INTENT(OUT), OPTIONAL :: default_periods
    CHARACTER(16),  INTENT(OUT), OPTIONAL :: default_start
    CHARACTER(16),  INTENT(OUT), OPTIONAL :: default_end
    CHARACTER(16),  INTENT(OUT), OPTIONAL :: default_interval
    CHARACTER(16),  INTENT(OUT), OPTIONAL :: default_histshift
    CHARACTER(16),  INTENT(OUT), OPTIONAL :: default_histcycle
    LOGICAL,        INTENT(OUT), OPTIONAL :: default_omit_time
    LOGICAL,        INTENT(OUT), OPTIONAL :: default_inquire

    CHARACTER(512) :: inputdir
    CHARACTER(8)   :: precision
    REAL(4)        :: missing(2)
    LOGICAL        :: descend
    INTEGER        :: fileview
    CHARACTER(20)  :: mode
    INTEGER        :: periods
    CHARACTER(16)  :: interval
    CHARACTER(16)  :: start, start_date
    CHARACTER(16)  :: end,   end_date
    CHARACTER(16)  :: histshift, histcycle
    LOGICAL        :: omit_time
    LOGICAL        :: inquire_exist

    INTEGER        :: iostat
    CHARACTER(256) :: iomsg

    NAMELIST / input_default / &
         inputdir,  &
         precision, &
         missing,   &
         descend,   &
         fileview,  &
         mode,      &
         periods,   &
         interval,  &
         start,     &
         start_date,&
         end,       &
         end_date,  &
         histshift, &
         histcycle, &
         omit_time, &
         inquire_exist

    IF (rank==0) THEN
       inputdir  = global_inputdir
       precision = 'REAL8'
       missing   = (/UNDEF, UNDEF/)
       descend   = .FALSE.
       fileview  = 0
       mode      = 'CONST'
       periods   = 1
       interval  = '0000001_000000'
       start     = format_datetime(t_start)
       end       = format_datetime(t_end)
       start_date= ''
       end_date  = ''
       histshift = ''
       histcycle = ''
       omit_time = .FALSE.
       inquire_exist = .TRUE.

       REWIND(CONFIG_UNIT)
       READ(CONFIG_UNIT, NML=input_default, IOSTAT=iostat, IOMSG=iomsg)

       CALL assert(iostat <= 0, "failed to read INPUT_DEFAULT namelist", iomsg)

       IF (start_date /= '') start = start_date
       IF (end_date   /= '') end   = end_date
    END IF

    CALL bcast(inputdir)
    CALL bcast(precision)
    CALL bcast(missing)
    CALL bcast(descend)
    CALL bcast(fileview)
    CALL bcast(mode)
    CALL bcast(periods)
    CALL bcast(interval)
    CALL bcast(start)
    CALL bcast(end)
    CALL bcast(histshift)
    CALL bcast(histcycle)
    CALL bcast(omit_time)
    CALL bcast(inquire_exist)

    IF (present(default_inputdir))     default_inputdir     = inputdir
    IF (present(default_precision))    default_precision    = precision
    IF (present(default_precision))    default_missing      = missing
    IF (present(default_descend))      default_descend      = descend
    IF (present(default_fileview))     default_fileview     = fileview
    IF (present(default_mode))         default_mode         = mode
    IF (present(default_periods))      default_periods      = periods
    IF (present(default_start))        default_start        = start
    IF (present(default_end))          default_end          = end
    IF (present(default_interval))     default_interval     = interval
    IF (present(default_omit_time))    default_omit_time    = omit_time
    IF (present(default_inquire))      default_inquire      = inquire_exist
    IF (present(default_histshift))    default_histshift    = histshift
    IF (present(default_histcycle))    default_histcycle    = histcycle

  END SUBROUTINE default_input

!-----------------------------------------------------------------------------------------------------------------------

  LOGICAL FUNCTION has_initial(varname)
    CHARACTER(*), INTENT(IN)  :: varname
    TYPE(initial_params) :: params

    CALL read_initial_namelist(varname, params, stat=has_initial)

  END FUNCTION has_initial

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE initial_data_2d_r4(varname, data, section, default, stat)
    CHARACTER(*), INTENT(IN)    :: varname
    REAL(4),      INTENT(INOUT) :: data(:,:)
    CHARACTER(2), INTENT(IN),  OPTIONAL :: section
    LOGICAL,      INTENT(IN),  OPTIONAL :: default
    LOGICAL,      INTENT(OUT), OPTIONAL :: stat
#ifdef F2008
    CONTIGUOUS data
#endif

    INTEGER :: nx, ny
    INTEGER :: i0, j0

    TYPE(initial_params) :: params

    REAL(8) :: tmp(isize, jsize)
    LOGICAL :: stat_

    IF (present(section)) THEN
       SELECT CASE(section)
       CASE ('XZ', 'xz')
          CALL initial_data_xz_r4(varname, data, default, stat)
          RETURN
       CASE ('YZ', 'yz')
          CALL initial_data_yz_r4(varname, data, default, stat)
          RETURN
       CASE ('XY', 'xy')
       CASE DEFAULT
          CALL assert(.FALSE., "invalid SECTION '"//section//"' in INITIAL_DATA_2D")
       END SELECT
    END IF

    CALL add_restart_var(varname, kind=4)

    IF (present(stat)) stat = .FALSE.

    nx = size(data,1) - isize
    ny = size(data,2) - jsize

    CALL assert(nx >= 0 .AND. nx <= 5, "unsupported dimension in INITIAL_DATA_2D for '"//trim(varname)//"'")
    CALL assert(ny >= 0 .AND. ny <= 5, "unsupported dimension in INITIAL_DATA_2D for '"//trim(varname)//"'")

    i0 = (nx+1)/2
    j0 = (ny+1)/2

    CALL read_initial_namelist(varname, params, default, stat=stat_)
    IF (.NOT. stat_) RETURN

!$OMP PARALLEL WORKSHARE
    data(:,:) = 0.0
!$OMP END PARALLEL WORKSHARE

    DO WHILE (stat_)
       CALL read_data_2d(tmp, params%filepath, kind=params%kind, view=params%fileview)
       CALL substmissing(tmp, params%missing(1), params%missing(2))
       CALL scaleoffset( tmp, params%scale,      params%offset)
       CALL limitminmax( tmp, params%minlimit,   params%maxlimit)
       CALL whitenoise(  tmp, params%noise, mkseed(r8=t_start, i4=jcoord*ipes+icoord, str=varname))

!$OMP PARALLEL WORKSHARE
       data(i0+1:i0+isize, j0+1:j0+jsize) = data(i0+1:i0+isize, j0+1:j0+jsize) + REAL(tmp(:,:), KIND=4)
!$OMP END PARALLEL WORKSHARE

       CALL read_initial_namelist(varname, params, rewind=.FALSE., stat=stat_)
    END DO

    CALL update_boundary_2d(data)

    IF (present(stat)) stat = .TRUE.

  END SUBROUTINE initial_data_2d_r4

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE initial_data_2d_r8(varname, data, section, default, stat)
    CHARACTER(*), INTENT(IN)    :: varname
    REAL(8),      INTENT(INOUT) :: data(:,:)
    CHARACTER(2), INTENT(IN),  OPTIONAL :: section
    LOGICAL,      INTENT(IN),  OPTIONAL :: default
    LOGICAL,      INTENT(OUT), OPTIONAL :: stat
#ifdef F2008
    CONTIGUOUS data
#endif

    INTEGER :: nx, ny
    INTEGER :: i0, j0

    TYPE(initial_params) :: params

    REAL(8) :: tmp(isize, jsize)
    LOGICAL :: stat_

    IF (present(section)) THEN
       SELECT CASE(section)
       CASE ('XZ', 'xz')
          CALL initial_data_xz_r8(varname, data, default, stat)
          RETURN
       CASE ('YZ', 'yz')
          CALL initial_data_yz_r8(varname, data, default, stat)
          RETURN
       CASE ('XY', 'xy')
       CASE DEFAULT
          CALL assert(.FALSE., "invalid SECTION '"//section//"' in INITIAL_DATA_2D")
       END SELECT
    END IF

    CALL add_restart_var(varname)

    IF (present(stat)) stat = .FALSE.

    nx = size(data,1) - isize
    ny = size(data,2) - jsize

    CALL assert(nx >= 0 .AND. nx <= 5, "unsupported dimension in INITIAL_DATA_2D for '"//trim(varname)//"'")
    CALL assert(ny >= 0 .AND. ny <= 5, "unsupported dimension in INITIAL_DATA_2D for '"//trim(varname)//"'")

    i0 = (nx+1)/2
    j0 = (ny+1)/2

    CALL read_initial_namelist(varname, params, default, stat=stat_)
    IF (.NOT. stat_) RETURN

    CALL assert(params%kind==8 .OR. .NOT. perfect_restart, "PERFECT_RESTART requires REAL8 (double-precision) initial data")

!$OMP PARALLEL WORKSHARE
    data(:,:) = 0.0
!$OMP END PARALLEL WORKSHARE

    DO WHILE (stat_)
       CALL read_data_2d(tmp, params%filepath, kind=params%kind, view=params%fileview)
       CALL substmissing(tmp, params%missing(1), params%missing(2))
       CALL scaleoffset( tmp, params%scale,      params%offset)
       CALL limitminmax( tmp, params%minlimit,   params%maxlimit)
       CALL whitenoise(  tmp, params%noise, mkseed(r8=t_start, i4=jcoord*ipes+icoord, str=varname))

!$OMP PARALLEL WORKSHARE
       data(i0+1:i0+isize, j0+1:j0+jsize) = data(i0+1:i0+isize, j0+1:j0+jsize) + tmp(:,:)
!$OMP END PARALLEL WORKSHARE

       CALL read_initial_namelist(varname, params, rewind=.FALSE., stat=stat_)
    END DO

    CALL update_boundary_2d(data)

    IF (present(stat)) stat = .TRUE.

  END SUBROUTINE initial_data_2d_r8

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE initial_data_xz_r4(varname, data, default, stat)
    CHARACTER(*), INTENT(IN)    :: varname
    REAL(4),      INTENT(INOUT) :: data(:,:)
    LOGICAL,      INTENT(IN),  OPTIONAL :: default
    LOGICAL,      INTENT(OUT), OPTIONAL :: stat
#ifdef F2008
    CONTIGUOUS data
#endif

    INTEGER :: nx, nz
    INTEGER :: i0, k0

    TYPE(initial_params) :: params

    REAL(8) :: tmp(isize, ksize)
    LOGICAL :: stat_

    CALL add_restart_var(varname, kind=4)

    IF (present(stat)) stat = .FALSE.

    nx = size(data,1) - isize
    nz = size(data,2) - ksize

    CALL assert(nx >= 0 .AND. nx <= 5, "unsupported dimension in INITIAL_DATA_XZ for '"//trim(varname)//"'")
    CALL assert(nz >= 0 .AND. nz <= 5, "unsupported dimension in INITIAL_DATA_XZ for '"//trim(varname)//"'")

    i0 = (nx+1)/2
    k0 = (nz+1)/2

    CALL read_initial_namelist(varname, params, default, stat=stat_)
    IF (.NOT. stat_) RETURN

!$OMP PARALLEL WORKSHARE
    data(:,:) = 0.0
!$OMP END PARALLEL WORKSHARE

    DO WHILE (stat_)
       CALL read_data_xz(tmp, params%filepath, kind=params%kind, descend=params%descend)
       CALL substmissing(tmp, params%missing(1), params%missing(2))
       CALL scaleoffset( tmp, params%scale,      params%offset)
       CALL limitminmax( tmp, params%minlimit,   params%maxlimit)
       CALL whitenoise(  tmp, params%noise, mkseed(r8=t_start, i4=kcoord*ipes+icoord, str=varname))

!$OMP PARALLEL WORKSHARE
       data(i0+1:i0+isize, k0+1:k0+ksize) = data(i0+1:i0+isize, k0+1:k0+ksize) + REAL(tmp(:,:), KIND=4)
!$OMP END PARALLEL WORKSHARE

       CALL read_initial_namelist(varname, params, rewind=.FALSE., stat=stat_)
    END DO

    CALL update_boundary_xz(data)

    IF (present(stat)) stat = .TRUE.

  END SUBROUTINE initial_data_xz_r4

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE initial_data_xz_r8(varname, data, default, stat)
    CHARACTER(*), INTENT(IN)    :: varname
    REAL(8),      INTENT(INOUT) :: data(:,:)
    LOGICAL,      INTENT(IN),  OPTIONAL :: default
    LOGICAL,      INTENT(OUT), OPTIONAL :: stat
#ifdef F2008
    CONTIGUOUS data
#endif

    INTEGER :: nx, nz
    INTEGER :: i0, k0

    TYPE(initial_params) :: params

    REAL(8) :: tmp(isize, ksize)
    LOGICAL :: stat_

    CALL add_restart_var(varname)

    IF (present(stat)) stat = .FALSE.

    nx = size(data,1) - isize
    nz = size(data,2) - ksize

    CALL assert(nx >= 0 .AND. nx <= 5, "unsupported dimension in INITIAL_DATA_XZ for '"//trim(varname)//"'")
    CALL assert(nz >= 0 .AND. nz <= 5, "unsupported dimension in INITIAL_DATA_XZ for '"//trim(varname)//"'")

    i0 = (nx+1)/2
    k0 = (nz+1)/2

    CALL read_initial_namelist(varname, params, default, stat=stat_)
    IF (.NOT. stat_) RETURN

    CALL assert(params%kind==8 .OR. .NOT. perfect_restart, "PERFECT_RESTART requires REAL8 (double-precision) initial data")

!$OMP PARALLEL WORKSHARE
    data(:,:) = 0.0
!$OMP END PARALLEL WORKSHARE

    DO WHILE (stat_)
       CALL read_data_xz(tmp, params%filepath, kind=params%kind, descend=params%descend)
       CALL substmissing(tmp, params%missing(1), params%missing(2))
       CALL scaleoffset( tmp, params%scale,      params%offset)
       CALL limitminmax( tmp, params%minlimit,   params%maxlimit)
       CALL whitenoise(  tmp, params%noise, mkseed(r8=t_start, i4=kcoord*ipes+icoord, str=varname))

!$OMP PARALLEL WORKSHARE
       data(i0+1:i0+isize, k0+1:k0+ksize) = data(i0+1:i0+isize, k0+1:k0+ksize) + tmp(:,:)
!$OMP END PARALLEL WORKSHARE

       CALL read_initial_namelist(varname, params, rewind=.FALSE., stat=stat_)
    END DO

    CALL update_boundary_xz(data)

    IF (present(stat)) stat = .TRUE.

  END SUBROUTINE initial_data_xz_r8

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE initial_data_yz_r4(varname, data, default, stat)
    CHARACTER(*), INTENT(IN)    :: varname
    REAL(4),      INTENT(INOUT) :: data(:,:)
    LOGICAL,      INTENT(IN),  OPTIONAL :: default
    LOGICAL,      INTENT(OUT), OPTIONAL :: stat
#ifdef F2008
    CONTIGUOUS data
#endif

    INTEGER :: ny, nz
    INTEGER :: j0, k0

    TYPE(initial_params) :: params

    REAL(8) :: tmp(jsize, ksize)
    LOGICAL :: stat_

    CALL add_restart_var(varname, kind=4)

    IF (present(stat)) stat = .FALSE.

    ny = size(data,1) - jsize
    nz = size(data,2) - ksize

    CALL assert(ny >= 0 .AND. ny <= 5, "unsupported dimension in INITIAL_DATA_YZ for '"//trim(varname)//"'")
    CALL assert(nz >= 0 .AND. nz <= 5, "unsupported dimension in INITIAL_DATA_YZ for '"//trim(varname)//"'")

    j0 = (ny+1)/2
    k0 = (nz+1)/2

    CALL read_initial_namelist(varname, params, default, stat=stat_)
    IF (.NOT. stat_) RETURN

!$OMP PARALLEL WORKSHARE
    data(:,:) = 0.0
!$OMP END PARALLEL WORKSHARE

    DO WHILE (stat_)
       CALL read_data_yz(tmp, params%filepath, kind=params%kind, descend=params%descend)
       CALL substmissing(tmp, params%missing(1), params%missing(2))
       CALL scaleoffset( tmp, params%scale,      params%offset)
       CALL limitminmax( tmp, params%minlimit,   params%maxlimit)
       CALL whitenoise(  tmp, params%noise, mkseed(r8=t_start, i4=kcoord*jpes+jcoord, str=varname))

!$OMP PARALLEL WORKSHARE
       data(j0+1:j0+jsize, k0+1:k0+ksize) = data(j0+1:j0+jsize, k0+1:k0+ksize) + REAL(tmp(:,:), KIND=4)
!$OMP END PARALLEL WORKSHARE

       CALL read_initial_namelist(varname, params, rewind=.FALSE., stat=stat_)
    END DO

    CALL update_boundary_yz(data)

    IF (present(stat)) stat = .TRUE.

  END SUBROUTINE initial_data_yz_r4

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE initial_data_yz_r8(varname, data, default, stat)
    CHARACTER(*), INTENT(IN)    :: varname
    REAL(8),      INTENT(INOUT) :: data(:,:)
    LOGICAL,      INTENT(IN),  OPTIONAL :: default
    LOGICAL,      INTENT(OUT), OPTIONAL :: stat
#ifdef F2008
    CONTIGUOUS data
#endif

    INTEGER :: ny, nz
    INTEGER :: j0, k0

    TYPE(initial_params) :: params

    REAL(8) :: tmp(jsize, ksize)
    LOGICAL :: stat_

    CALL add_restart_var(varname)

    IF (present(stat)) stat = .FALSE.

    ny = size(data,1) - jsize
    nz = size(data,2) - ksize

    CALL assert(ny >= 0 .AND. ny <= 5, "unsupported dimension in INITIAL_DATA_YZ for '"//trim(varname)//"'")
    CALL assert(nz >= 0 .AND. nz <= 5, "unsupported dimension in INITIAL_DATA_YZ for '"//trim(varname)//"'")

    j0 = (ny+1)/2
    k0 = (nz+1)/2

    CALL read_initial_namelist(varname, params, default, stat=stat_)
    IF (.NOT. stat_) RETURN

    CALL assert(params%kind==8 .OR. .NOT. perfect_restart, "PERFECT_RESTART requires REAL8 (double-precision) initial data")

!$OMP PARALLEL WORKSHARE
    data(:,:) = 0.0
!$OMP END PARALLEL WORKSHARE

    DO WHILE (stat_)
       CALL read_data_yz(tmp, params%filepath, kind=params%kind, descend=params%descend)
       CALL substmissing(tmp, params%missing(1), params%missing(2))
       CALL scaleoffset( tmp, params%scale,      params%offset)
       CALL limitminmax( tmp, params%minlimit,   params%maxlimit)
       CALL whitenoise(  tmp, params%noise, mkseed(r8=t_start, i4=kcoord*jpes+jcoord, str=varname))

!$OMP PARALLEL WORKSHARE
       data(j0+1:j0+jsize, k0+1:k0+ksize) = data(j0+1:j0+jsize, k0+1:k0+ksize) + tmp(:,:)
!$OMP END PARALLEL WORKSHARE

       CALL read_initial_namelist(varname, params, rewind=.FALSE., stat=stat_)
    END DO

    CALL update_boundary_yz(data)

    IF (present(stat)) stat = .TRUE.

  END SUBROUTINE initial_data_yz_r8

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE initial_data_3d_r4(varname, data, default, stat)
    CHARACTER(*), INTENT(IN)    :: varname
    REAL(4),      INTENT(INOUT) :: data(:,:,:)
    LOGICAL,      INTENT(IN),  OPTIONAL :: default
    LOGICAL,      INTENT(OUT), OPTIONAL :: stat
#ifdef F2008
    CONTIGUOUS data
#endif

    INTEGER :: nx, ny, nz
    INTEGER :: i0, j0, k0

    TYPE(initial_params) :: params

    REAL(8) :: tmp(isize, jsize, ksize)
    LOGICAL :: stat_

    CALL add_restart_var(varname, kind=4)

    IF (present(stat)) stat = .FALSE.

    nx = size(data,1) - isize
    ny = size(data,2) - jsize
    nz = size(data,3) - ksize

    CALL assert(nx >= 0 .AND. nx <= 5, "unsupported dimension in INITIAL_DATA_3D for '"//trim(varname)//"'")
    CALL assert(ny >= 0 .AND. ny <= 5, "unsupported dimension in INITIAL_DATA_3D for '"//trim(varname)//"'")
    CALL assert(nz >= 0 .AND. nz <= 5, "unsupported dimension in INITIAL_DATA_3D for '"//trim(varname)//"'")

    i0 = (nx+1)/2
    j0 = (ny+1)/2
    k0 = (nz+1)/2

    CALL read_initial_namelist(varname, params, default, stat=stat_)
    IF (.NOT. stat_) RETURN

!$OMP PARALLEL WORKSHARE
       data(:,:,:) = 0.0
!$OMP END PARALLEL WORKSHARE

    DO WHILE (stat_)
       CALL read_data_3d(tmp, params%filepath, kind=params%kind, view=params%fileview, descend=params%descend)
       CALL substmissing(tmp, params%missing(1), params%missing(2))
       CALL scaleoffset( tmp, params%scale,      params%offset)
       CALL limitminmax( tmp, params%minlimit,   params%maxlimit)
       CALL whitenoise(  tmp, params%noise, mkseed(r8=t_start, i4=kcoord*jpes*ipes+jcoord*ipes+icoord, str=varname))

!$OMP PARALLEL WORKSHARE
       data(i0+1:i0+isize, j0+1:j0+jsize, k0+1:k0+ksize) = data(i0+1:i0+isize, j0+1:j0+jsize, k0+1:k0+ksize) + REAL(tmp(:,:,:), KIND=4)
!$OMP END PARALLEL WORKSHARE

       CALL read_initial_namelist(varname, params, rewind=.FALSE., stat=stat_)
    END DO

    CALL update_boundary_3d(data)

    IF (present(stat)) stat = .TRUE.

  END SUBROUTINE initial_data_3d_r4

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE initial_data_3d_r8(varname, data, default, stat)
    CHARACTER(*), INTENT(IN)    :: varname
    REAL(8),      INTENT(INOUT) :: data(:,:,:)
    LOGICAL,      INTENT(IN),  OPTIONAL :: default
    LOGICAL,      INTENT(OUT), OPTIONAL :: stat
#ifdef F2008
    CONTIGUOUS data
#endif

    INTEGER :: nx, ny, nz
    INTEGER :: i0, j0, k0

    TYPE(initial_params) :: params

    REAL(8) :: tmp(isize, jsize, ksize)
    LOGICAL :: stat_

    CALL add_restart_var(varname)

    IF (present(stat)) stat = .FALSE.

    nx = size(data,1) - isize
    ny = size(data,2) - jsize
    nz = size(data,3) - ksize

    CALL assert(nx >= 0 .AND. nx <= 5, "unsupported dimension in INITIAL_DATA_3D for '"//trim(varname)//"'")
    CALL assert(ny >= 0 .AND. ny <= 5, "unsupported dimension in INITIAL_DATA_3D for '"//trim(varname)//"'")
    CALL assert(nz >= 0 .AND. nz <= 5, "unsupported dimension in INITIAL_DATA_3D for '"//trim(varname)//"'")

    i0 = (nx+1)/2
    j0 = (ny+1)/2
    k0 = (nz+1)/2

    CALL read_initial_namelist(varname, params, default, stat=stat_)
    IF (.NOT. stat_) RETURN

    CALL assert(params%kind==8 .OR. .NOT. perfect_restart, "PERFECT_RESTART requires REAL8 (double-precision) initial data")

!$OMP PARALLEL WORKSHARE
       data(:,:,:) = 0.0
!$OMP END PARALLEL WORKSHARE

    DO WHILE (stat_)
       CALL read_data_3d(tmp, params%filepath, kind=params%kind, view=params%fileview, descend=params%descend)
       CALL substmissing(tmp, params%missing(1), params%missing(2))
       CALL scaleoffset( tmp, params%scale,      params%offset)
       CALL limitminmax( tmp, params%minlimit,   params%maxlimit)
       CALL whitenoise(  tmp, params%noise, mkseed(r8=t_start, i4=kcoord*jpes*ipes+jcoord*ipes+icoord, str=varname))

!$OMP PARALLEL WORKSHARE
       data(i0+1:i0+isize, j0+1:j0+jsize, k0+1:k0+ksize) = data(i0+1:i0+isize, j0+1:j0+jsize, k0+1:k0+ksize) + tmp(:,:,:)
!$OMP END PARALLEL WORKSHARE

       CALL read_initial_namelist(varname, params, rewind=.FALSE., stat=stat_)
    END DO

    CALL update_boundary_3d(data)

    IF (present(stat)) stat = .TRUE.

  END SUBROUTINE initial_data_3d_r8

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE initial_data_1d_r4(varname, data, axis, default, stat)
    CHARACTER(*), INTENT(IN)    :: varname
    REAL(4),      INTENT(INOUT) :: data(:)
    CHARACTER(1), INTENT(IN)    :: axis
    LOGICAL,      INTENT(IN),  OPTIONAL :: default
    LOGICAL,      INTENT(OUT), OPTIONAL :: stat
#ifdef F2008
    CONTIGUOUS data
#endif

    SELECT CASE(axis)
    CASE ('X', 'x')
       CALL initial_data_x_r4(varname, data, default, stat)
    CASE ('Y', 'y')
       CALL initial_data_y_r4(varname, data, default, stat)
    CASE ('Z', 'z')
       CALL initial_data_z_r4(varname, data, default, stat)
    CASE DEFAULT
       CALL assert(.FALSE., "invalid AXIS '"//axis//"' in INITIAL_DATA_1D")
    END SELECT

  END SUBROUTINE initial_data_1d_r4

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE initial_data_1d_r8(varname, data, axis, default, stat)
    CHARACTER(*), INTENT(IN)    :: varname
    REAL(8),      INTENT(INOUT) :: data(:)
    CHARACTER(1), INTENT(IN)    :: axis
    LOGICAL,      INTENT(IN),  OPTIONAL :: default
    LOGICAL,      INTENT(OUT), OPTIONAL :: stat
#ifdef F2008
    CONTIGUOUS data
#endif

    SELECT CASE(axis)
    CASE ('X', 'x')
       CALL initial_data_x_r8(varname, data, default, stat)
    CASE ('Y', 'y')
       CALL initial_data_y_r8(varname, data, default, stat)
    CASE ('Z', 'z')
       CALL initial_data_z_r8(varname, data, default, stat)
    CASE DEFAULT
       CALL assert(.FALSE., "invalid AXIS '"//axis//"' in INITIAL_DATA_1D")
    END SELECT

  END SUBROUTINE initial_data_1d_r8

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE initial_data_x_r4(varname, data, default, stat)
    CHARACTER(*), INTENT(IN)    :: varname
    REAL(4),      INTENT(INOUT) :: data(:)
    LOGICAL,      INTENT(IN),  OPTIONAL :: default
    LOGICAL,      INTENT(OUT), OPTIONAL :: stat
#ifdef F2008
    CONTIGUOUS data
#endif

    INTEGER :: nx
    INTEGER :: i0

    TYPE(initial_params) :: params

    REAL(8) :: tmp(isize)
    LOGICAL :: stat_

    CALL add_restart_var(varname, kind=4)

    IF (present(stat)) stat = .FALSE.

    nx = size(data) - isize

    CALL assert(nx >= 0 .AND. nx <= 5, "unsupported dimension in INITIAL_DATA_X for '"//trim(varname)//"'")

    i0 = (nx+1)/2

    CALL read_initial_namelist(varname, params, default, stat=stat_)
    IF (.NOT. stat_) RETURN

    data(:) = 0.0

    CALL assert(params%fileview==0, "FILEVIEW is not supported for 1D data")

    DO WHILE (stat_)
       CALL read_data_x(tmp, params%filepath, kind=params%kind)
       CALL substmissing(tmp, params%missing(1), params%missing(2))
       CALL scaleoffset( tmp, params%scale,      params%offset)
       CALL limitminmax( tmp, params%minlimit,   params%maxlimit)
       CALL whitenoise(  tmp, params%noise, mkseed(r8=t_start, i4=icoord, str=varname))

       data(i0+1:i0+isize) = data(i0+1:i0+isize) + REAL(tmp(:), KIND=4)

       CALL read_initial_namelist(varname, params, rewind=.FALSE., stat=stat_)
    END DO

    CALL update_boundary_x(data)

    IF (present(stat)) stat = .TRUE.

  END SUBROUTINE initial_data_x_r4

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE initial_data_x_r8(varname, data, default, stat)
    CHARACTER(*), INTENT(IN)    :: varname
    REAL(8),      INTENT(INOUT) :: data(:)
    LOGICAL,      INTENT(IN),  OPTIONAL :: default
    LOGICAL,      INTENT(OUT), OPTIONAL :: stat
#ifdef F2008
    CONTIGUOUS data
#endif

    INTEGER :: nx
    INTEGER :: i0

    TYPE(initial_params) :: params

    REAL(8) :: tmp(isize)
    LOGICAL :: stat_

    CALL add_restart_var(varname)

    IF (present(stat)) stat = .FALSE.

    nx = size(data) - isize

    CALL assert(nx >= 0 .AND. nx <= 5, "unsupported dimension in INITIAL_DATA_X for '"//trim(varname)//"'")

    i0 = (nx+1)/2

    CALL read_initial_namelist(varname, params, default, stat=stat_)
    IF (.NOT. stat_) RETURN

    data(:) = 0.0

    CALL assert(params%fileview==0, "FILEVIEW is not supported for 1D data")
    CALL assert(params%kind==8 .OR. .NOT. perfect_restart, "PERFECT_RESTART requires REAL8 (double-precision) initial data")

    DO WHILE (stat_)
       CALL read_data_x(tmp, params%filepath, kind=params%kind)
       CALL substmissing(tmp, params%missing(1), params%missing(2))
       CALL scaleoffset( tmp, params%scale,      params%offset)
       CALL limitminmax( tmp, params%minlimit,   params%maxlimit)
       CALL whitenoise(  tmp, params%noise, mkseed(r8=t_start, i4=icoord, str=varname))

       data(i0+1:i0+isize) = data(i0+1:i0+isize) + tmp(:)

       CALL read_initial_namelist(varname, params, rewind=.FALSE., stat=stat_)
    END DO

    CALL update_boundary_x(data)

    IF (present(stat)) stat = .TRUE.

  END SUBROUTINE initial_data_x_r8

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE initial_data_y_r4(varname, data, default, stat)
    CHARACTER(*), INTENT(IN)    :: varname
    REAL(4),      INTENT(INOUT) :: data(:)
    LOGICAL,      INTENT(IN),  OPTIONAL :: default
    LOGICAL,      INTENT(OUT), OPTIONAL :: stat
#ifdef F2008
    CONTIGUOUS data
#endif

    INTEGER :: ny
    INTEGER :: j0

    TYPE(initial_params) :: params

    REAL(8) :: tmp(jsize)
    LOGICAL :: stat_

    CALL add_restart_var(varname, kind=4)

    IF (present(stat)) stat = .FALSE.

    ny = size(data) - jsize

    CALL assert(ny >= 0 .AND. ny <= 5, "unsupported dimension in INITIAL_DATA_Y for '"//trim(varname)//"'")

    j0 = (ny+1)/2

    CALL read_initial_namelist(varname, params, default, stat=stat_)
    IF (.NOT. stat_) RETURN

    CALL assert(params%fileview==0, "FILEVIEW is not supported for 1D data")

    data(:) = 0.0

    DO WHILE (stat_)
       CALL read_data_y(tmp, params%filepath, kind=params%kind)
       CALL substmissing(tmp, params%missing(1), params%missing(2))
       CALL scaleoffset( tmp, params%scale,      params%offset)
       CALL limitminmax( tmp, params%minlimit,   params%maxlimit)
       CALL whitenoise(  tmp, params%noise, mkseed(r8=t_start, i4=jcoord, str=varname))

       data(j0+1:j0+jsize) = REAL(tmp(:), KIND=4)

       CALL read_initial_namelist(varname, params, rewind=.FALSE., stat=stat_)
    END DO

    CALL update_boundary_y(data)

    IF (present(stat)) stat = .TRUE.

  END SUBROUTINE initial_data_y_r4

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE initial_data_y_r8(varname, data, default, stat)
    CHARACTER(*), INTENT(IN)    :: varname
    REAL(8),      INTENT(INOUT) :: data(:)
    LOGICAL,      INTENT(IN),  OPTIONAL :: default
    LOGICAL,      INTENT(OUT), OPTIONAL :: stat
#ifdef F2008
    CONTIGUOUS data
#endif

    INTEGER :: ny
    INTEGER :: j0

    TYPE(initial_params) :: params

    REAL(8) :: tmp(jsize)
    LOGICAL :: stat_

    CALL add_restart_var(varname)

    IF (present(stat)) stat = .FALSE.

    ny = size(data) - jsize

    CALL assert(ny >= 0 .AND. ny <= 5, "unsupported dimension in INITIAL_DATA_Y for '"//trim(varname)//"'")

    j0 = (ny+1)/2

    CALL read_initial_namelist(varname, params, default, stat=stat_)
    IF (.NOT. stat_) RETURN

    CALL assert(params%fileview==0, "FILEVIEW is not supported for 1D data")
    CALL assert(params%kind==8 .OR. .NOT. perfect_restart, "PERFECT_RESTART requires REAL8 (double-precision) initial data")

    data(:) = 0.0

    DO WHILE (stat_)
       CALL read_data_y(tmp, params%filepath, kind=params%kind)
       CALL substmissing(tmp, params%missing(1), params%missing(2))
       CALL scaleoffset( tmp, params%scale,      params%offset)
       CALL limitminmax( tmp, params%minlimit,   params%maxlimit)
       CALL whitenoise(  tmp, params%noise, mkseed(r8=t_start, i4=jcoord, str=varname))

       data(j0+1:j0+jsize) = tmp(:)

       CALL read_initial_namelist(varname, params, rewind=.FALSE., stat=stat_)
    END DO

    CALL update_boundary_y(data)

    IF (present(stat)) stat = .TRUE.

  END SUBROUTINE initial_data_y_r8

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE initial_data_z_r4(varname, data, default, stat)
    CHARACTER(*), INTENT(IN)    :: varname
    REAL(4),      INTENT(INOUT) :: data(:)
    LOGICAL,      INTENT(IN),  OPTIONAL :: default
    LOGICAL,      INTENT(OUT), OPTIONAL :: stat
#ifdef F2008
    CONTIGUOUS data
#endif

    INTEGER :: nz
    INTEGER :: k0

    TYPE(initial_params) :: params

    REAL(8) :: tmp(ksize)
    LOGICAL :: stat_

    CALL add_restart_var(varname, kind=4)

    IF (present(stat)) stat = .FALSE.

    nz = size(data) - ksize

    CALL assert(nz >= 0 .AND. nz <= 5, "unsupported dimension in INITIAL_DATA_Z for '"//trim(varname)//"'")

    k0 = (nz+1)/2

    CALL read_initial_namelist(varname, params, default, stat=stat_)
    IF (.NOT. stat_) RETURN

    CALL assert(params%fileview==0, "FILEVIEW is not supported for 1D data")

    data(:) = 0.0

    DO WHILE (stat_)
       CALL read_data_z(tmp, params%filepath, kind=params%kind)
       CALL substmissing(tmp, params%missing(1), params%missing(2))
       CALL scaleoffset( tmp, params%scale,      params%offset)
       CALL limitminmax( tmp, params%minlimit,   params%maxlimit)
       CALL whitenoise(  tmp, params%noise, mkseed(r8=t_start, i4=kcoord, str=varname))

       data(k0+1:k0+ksize) = REAL(tmp(:), KIND=4)

       CALL read_initial_namelist(varname, params, rewind=.FALSE., stat=stat_)
    END DO

    CALL update_boundary_z(data)

    IF (present(stat)) stat = .TRUE.

  END SUBROUTINE initial_data_z_r4

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE initial_data_z_r8(varname, data, default, stat)
    CHARACTER(*), INTENT(IN)    :: varname
    REAL(8),      INTENT(INOUT) :: data(:)
    LOGICAL,      INTENT(IN),  OPTIONAL :: default
    LOGICAL,      INTENT(OUT), OPTIONAL :: stat
#ifdef F2008
    CONTIGUOUS data
#endif

    INTEGER :: nz
    INTEGER :: k0

    TYPE(initial_params) :: params

    REAL(8) :: tmp(ksize)
    LOGICAL :: stat_

    CALL add_restart_var(varname)

    IF (present(stat)) stat = .FALSE.

    nz = size(data) - ksize

    CALL assert(nz >= 0 .AND. nz <= 5, "unsupported dimension in INITIAL_DATA_Z for '"//trim(varname)//"'")

    k0 = (nz+1)/2

    CALL read_initial_namelist(varname, params, default, stat=stat_)
    IF (.NOT. stat_) RETURN

    data(:) = 0.0

    CALL assert(params%fileview==0, "FILEVIEW is not supported for 1D data")
    CALL assert(params%kind==8 .OR. .NOT. perfect_restart, "PERFECT_RESTART requires REAL8 (double-precision) initial data")

    DO WHILE (stat_)
       CALL read_data_z(tmp, params%filepath, kind=params%kind, descend=params%descend)
       CALL substmissing(tmp, params%missing(1), params%missing(2))
       CALL scaleoffset( tmp, params%scale,      params%offset)
       CALL limitminmax( tmp, params%minlimit,   params%maxlimit)
       CALL whitenoise(  tmp, params%noise, mkseed(r8=t_start, i4=kcoord, str=varname))

       data(k0+1:k0+ksize) = data(k0+1:k0+ksize) + tmp(:)

       CALL read_initial_namelist(varname, params, rewind=.FALSE., stat=stat_)
    END DO

    CALL update_boundary_z(data)

    IF (present(stat)) stat = .TRUE.

  END SUBROUTINE initial_data_z_r8

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE read_initial_namelist(var, params, default, rewind, stat)
    CHARACTER(*), INTENT(IN) :: var
    TYPE(initial_params), INTENT(OUT) :: params
    LOGICAL, INTENT(IN),  OPTIONAL :: default
    LOGICAL, INTENT(IN),  OPTIONAL :: rewind
    LOGICAL, INTENT(OUT), OPTIONAL :: stat

    CHARACTER(32)  :: varname
    CHARACTER(512) :: initialdir
    CHARACTER(512) :: inputdir
    CHARACTER(8)   :: precision
    CHARACTER(128) :: filename
    REAL(4)        :: scale
    REAL(4)        :: minlimit
    REAL(4)        :: maxlimit
    REAL(4)        :: offset
    REAL(4)        :: missing(2)
    REAL(4)        :: noise
    INTEGER        :: fileview
    LOGICAL        :: descend

    INTEGER        :: iostat
    CHARACTER(256) :: iomsg
    LOGICAL        :: rewind_
    LOGICAL        :: stat_

    INTEGER :: i

    NAMELIST / initial / &
         varname,   &
         initialdir,&
         inputdir,  &
         filename,  &
         scale,     &
         offset,    &
         minlimit,  &
         maxlimit,  &
         missing,   &
         noise,     &
         fileview,  &
         descend,   &
         precision

    rewind_ = .TRUE.
    IF (present(rewind)) rewind_ = rewind

    stat_ = .FALSE.

    IF (rank==0) THEN
       IF (rewind_) REWIND(CONFIG_UNIT)
       DO
          varname  = ''
          initialdir = global_initialdir
          inputdir = ''
          filename = ''
          scale    = 1.0
          offset   = 0.0
          minlimit = UNDEF
          maxlimit = UNDEF
          missing  = (/UNDEF, UNDEF/)
          noise    = 0.0
          fileview = 0
          descend  = .FALSE.
          precision = 'REAL8'
          read(CONFIG_UNIT, NML=initial, IOSTAT=iostat, IOMSG=iomsg)

          IF (iostat < 0) EXIT

          CALL assert(iostat == 0, "failed to read INITIAL namelist for VARNAME='"//trim(varname)//"'", iomsg)
          CALL assert(trim(varname)/='', "VARNAME is mandatory for INITIAL namelist")
          IF (trim(varname)==trim(var)) THEN
             stat_ = .TRUE.
             IF (inputdir /= '') initialdir = inputdir !for backward compatibility
             EXIT
          END IF
       END DO

       IF (present(default)) stat_ = stat_ .OR. default

       IF (stat_ .AND. filename == "")  filename = trim(var) // "." // format_datetime(t_current)

       CALL replace_vars(filename)
       CALL replace_vars(initialdir, default=global_initialdir)
    END IF

    CALL bcast(stat_)

    CALL bcast(filename)
    CALL bcast(initialdir)
    CALL bcast(scale)
    CALL bcast(offset)
    CALL bcast(minlimit)
    CALL bcast(maxlimit)
    CALL bcast(missing)
    CALL bcast(noise)
    CALL bcast(fileview)
    CALL bcast(descend)
    CALL bcast(precision)

    params%filepath   = path(initialdir, filename)
    params%scale      = scale
    params%offset     = offset
    params%minlimit   = minlimit
    params%maxlimit   = maxlimit
    params%missing    = missing
    params%noise      = noise
    params%fileview   = fileview
    params%descend    = descend

    SELECT CASE (trim(precision))
    CASE ('BYTE',   'byte',   'INT1',  'int1',  'I1', 'i1', '1')
       params%kind = 1
    CASE ('SINGLE', 'single', 'REAL4', 'real4', 'R4', 'r4', '4')
       params%kind = 4
    CASE ('DOUBLE', 'double', 'REAL8', 'real8', 'R8', 'r8', '8')
       params%kind = 8
    CASE DEFAULT
       CALL assert(.FALSE., "unsupported PRECISION '"//trim(precision)//"'")
    END SELECT

    IF (present(stat)) stat = stat_

  END SUBROUTINE read_initial_namelist

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE replace_vars(str, default)
    CHARACTER(*), INTENT(INOUT) :: str
    CHARACTER(*), INTENT(IN), OPTIONAL :: default

    INTEGER :: stat

    CALL replace('$RUNNAME',      trim(runname),                    str)
    CALL replace('$START_DATE',   trim(strip_date(start_datetime)), str)
    CALL replace('$START_TIME',   trim(strip_time(start_datetime)), str)
    CALL replace('$START',        trim(start_datetime),             str)

    IF (present(default)) CALL replace('$DEFAULT', trim(default), str)

    CALL replace('$WORKDIR',  trim(workdir),  str)
    CALL replace('$CACHEDIR', trim(cachedir), str)

  END SUBROUTINE replace_vars

  SUBROUTINE replace_datetime(str, t)
    CHARACTER(*), INTENT(INOUT) :: str
    REAL(8),      INTENT(IN)    :: t

    CHARACTER(16) :: s

    s = format_datetime(t, fmt=0)

    CALL replace('$YEAR',  s(1:4),  str)
    CALL replace('$YEAR2', s(3:4),  str)
    CALL replace('$DOY',   s(5:7),  str)

    CALL replace('$HH',    s( 9:10), str)
    CALL replace('$MM',    s(11:12), str)
    CALL replace('$SS',    s(13:14), str)

    s = format_datetime(t, fmt=1, omit_time=.TRUE.)

    CALL replace('$MON',   s(5:6),  str)
    CALL replace('$DOM',   s(7:8),  str)

  END SUBROUTINE replace_datetime

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE read_dump_namelist
    CHARACTER(16) :: start    = ""
    CHARACTER(16) :: end      = ""
    CHARACTER(16) :: interval = ""

    INTEGER :: iostat
    CHARACTER(256) :: iomsg

    NAMELIST / dump /   &
         start,         &
         end,           &
         interval,      &
         dumpdir

    IF (rank==0) THEN
       REWIND(CONFIG_UNIT)
       READ(CONFIG_UNIT, NML=dump, IOSTAT=iostat, IOMSG=iomsg)
       CALL assert(iostat <= 0, "failed to read DUMP namelist", iomsg)
    END IF
    CALL bcast(iostat)
    IF (iostat < 0) RETURN

    CALL bcast(start)
    CALL bcast(end)
    CALL bcast(interval)
    CALL bcast(dumpdir)

    IF (start == '') THEN
       dump_start = t_start
    ELSE
       dump_start = datetime_seconds(start)
    END IF

    IF (end == '') THEN
       dump_end = t_end
    ELSE IF (end(1:1)=='+') THEN
       dump_end = dump_start + interval_seconds(end(2:))
    ELSE
       dump_end = datetime_seconds(end)
    END IF

    IF (interval == '') THEN
       dump_interval  = dtime
    ELSE
       dump_interval  = interval_seconds(interval)
    END IF

    CALL assert(dumpdir/='', "DUMPDIR is mandatory for DUMP namelist")
    CALL replace_vars(dumpdir)
    IF (rank==0) WRITE(REPORT_UNIT, *) "dumpdir = ", dumpdir

    lastdump = dump_start - dump_interval


  END SUBROUTINE read_dump_namelist

  SUBROUTINE dump_data_4d(x)
    REAL(8), INTENT(IN) :: x(:,:,:,:)

    CHARACTER(512) :: dumpfile
    CHARACTER(8)   :: status
    INTEGER        :: iostat
    CHARACTER(512) :: iomsg

    IF (dump_interval < 0) RETURN
    IF (t_current < dump_start .OR. t_current > dump_end) RETURN

    IF (t_current >= lastdump + dump_interval) THEN
       lastdump = t_current
       status   = 'REPLACE'
    ELSE IF (t_current == lastdump) THEN
       status = 'OLD'
    ELSE
       RETURN
    END IF

    dumpfile = trim(dumpdir) // '/DUMP.' // format_datetime(lastdump)
#ifdef PARALLEL_MPI
    dumpfile = trim(dumpfile) // '.' // strrank
#endif

    OPEN(UNIT     = DUMP_UNIT,       &
         FILE     = trim(dumpfile),  &
         FORM     = 'UNFORMATTED',   &
         ACCESS   = 'STREAM',        &
         STATUS   = status,          &
         POSITION = 'APPEND',        &
         ACTION   = 'WRITE',         &
         IOMSG    = iomsg,           &
         IOSTAT   = iostat)
    CALL assert(iostat==0, iomsg)

    WRITE(DUMP_UNIT, IOSTAT=iostat, IOMSG=iomsg) x
    CALL assert(iostat==0, iomsg)

    CLOSE(DUMP_UNIT, IOSTAT=iostat, IOMSG=iomsg)
    CALL assert(iostat==0, iomsg)

  END SUBROUTINE dump_data_4d

  SUBROUTINE dump_data_3d(x)
    REAL(8), INTENT(IN) :: x(:,:,:)

    CALL dump_data_4d(reshape(x, (/size(x,1), size(x,2), size(x,3), 1/)))

  END SUBROUTINE dump_data_3d

  SUBROUTINE dump_data_2d(x)
    REAL(8), INTENT(IN) :: x(:,:)

    CALL dump_data_4d(reshape(x, (/size(x,1), size(x,2), 1, 1/)))

  END SUBROUTINE dump_data_2d

  SUBROUTINE dump_data_1d(x)
    REAL(8), INTENT(IN) :: x(:)

    CALL dump_data_4d(reshape(x, (/size(x), 1, 1, 1/)))

  END SUBROUTINE dump_data_1d

END MODULE io
