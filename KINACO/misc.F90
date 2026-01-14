#include "macro.h"

MODULE misc
  IMPLICIT NONE
  INTEGER, PARAMETER :: NULL  = 0
  REAL(8), PARAMETER :: PI    = 3.141592653589793D0

  REAL(4), PARAMETER :: UNDEF = 999.0E9

  REAL(8), PARAMETER :: SQRT2 = 1.414213562373095D0
  REAL(8), PARAMETER :: SQRT3 = 1.732050807568877D0

  INTEGER, PARAMETER :: STDIN_UNIT  = 5
  INTEGER, PARAMETER :: STDOUT_UNIT = 6
  INTEGER, PARAMETER :: STDERR_UNIT = 0

  INTEGER, PARAMETER :: CONFIG_UNIT = 99

  INTEGER, PARAMETER :: TMP_UNIT = 90

  CHARACTER(32), SAVE :: runname = 'UNTITLED'

  INTEGER, SAVE :: REPORT_UNIT = 98
  INTEGER, SAVE :: DUMP_UNIT = 97

  CHARACTER(32) :: str_format
#ifdef PARALLEL_MPI
  REAL(8) :: cputime1, cputime2
#else
  REAL    :: cputime1, cputime2
#endif

  INTERFACE scaleoffset
     MODULE PROCEDURE scaleoffset_3d_r4
     MODULE PROCEDURE scaleoffset_3d_r8
     MODULE PROCEDURE scaleoffset_2d_r4
     MODULE PROCEDURE scaleoffset_2d_r8
     MODULE PROCEDURE scaleoffset_1d_r4
     MODULE PROCEDURE scaleoffset_1d_r8
  END INTERFACE scaleoffset

  INTERFACE limitminmax
     MODULE PROCEDURE limitminmax_3d_r4
     MODULE PROCEDURE limitminmax_3d_r8
     MODULE PROCEDURE limitminmax_2d_r4
     MODULE PROCEDURE limitminmax_2d_r8
     MODULE PROCEDURE limitminmax_1d_r4
     MODULE PROCEDURE limitminmax_1d_r8
  END INTERFACE limitminmax

  INTERFACE substmissing
     MODULE PROCEDURE substmissing_3d_r4
     MODULE PROCEDURE substmissing_3d_r8
     MODULE PROCEDURE substmissing_2d_r4
     MODULE PROCEDURE substmissing_2d_r8
     MODULE PROCEDURE substmissing_1d_r4
     MODULE PROCEDURE substmissing_1d_r8
  END INTERFACE substmissing

  INTERFACE whitenoise
     MODULE PROCEDURE whitenoise_3d_r4
     MODULE PROCEDURE whitenoise_3d_r8
     MODULE PROCEDURE whitenoise_2d_r4
     MODULE PROCEDURE whitenoise_2d_r8
     MODULE PROCEDURE whitenoise_1d_r4
     MODULE PROCEDURE whitenoise_1d_r8
  END INTERFACE whitenoise

  INTERFACE write_file
     MODULE PROCEDURE write_file_3d_i1
     MODULE PROCEDURE write_file_3d_r4
     MODULE PROCEDURE write_file_3d_r8
     MODULE PROCEDURE write_file_2d_i1
     MODULE PROCEDURE write_file_2d_r4
     MODULE PROCEDURE write_file_2d_r8
     MODULE PROCEDURE write_file_1d_i1
     MODULE PROCEDURE write_file_1d_r4
     MODULE PROCEDURE write_file_1d_r8
  END INTERFACE

  INTERFACE read_file
     MODULE PROCEDURE read_file_3d_i1
     MODULE PROCEDURE read_file_3d_r4
     MODULE PROCEDURE read_file_3d_r8
     MODULE PROCEDURE read_file_2d_i1
     MODULE PROCEDURE read_file_2d_r4
     MODULE PROCEDURE read_file_2d_r8
     MODULE PROCEDURE read_file_1d_i1
     MODULE PROCEDURE read_file_1d_r4
     MODULE PROCEDURE read_file_1d_r8
  END INTERFACE


  INTERFACE kreverse
     MODULE PROCEDURE kreverse_1d_r4
     MODULE PROCEDURE kreverse_1d_r8
     MODULE PROCEDURE kreverse_2d_r4
     MODULE PROCEDURE kreverse_2d_r8
     MODULE PROCEDURE kreverse_3d_r4
     MODULE PROCEDURE kreverse_3d_r8
  END INTERFACE


  INTERFACE format
     MODULE PROCEDURE format_int4
     MODULE PROCEDURE format_int8
     MODULE PROCEDURE format_real4
     MODULE PROCEDURE format_real8
  END INTERFACE

  INTERFACE nan_to_undef
     MODULE PROCEDURE nan_to_undef_1d_r4
     MODULE PROCEDURE nan_to_undef_2d_r4
     MODULE PROCEDURE nan_to_undef_3d_r4
     MODULE PROCEDURE nan_to_undef_1d_r8
     MODULE PROCEDURE nan_to_undef_2d_r8
     MODULE PROCEDURE nan_to_undef_3d_r8
  END INTERFACE

  INTERFACE check_nan
     MODULE PROCEDURE check_nan_1d_r4
     MODULE PROCEDURE check_nan_2d_r4
     MODULE PROCEDURE check_nan_3d_r4
     MODULE PROCEDURE check_nan_1d_r8
     MODULE PROCEDURE check_nan_2d_r8
     MODULE PROCEDURE check_nan_3d_r8
  END INTERFACE

  INTERFACE copy
     MODULE PROCEDURE copy_1d_r4
     MODULE PROCEDURE copy_2d_r4
     MODULE PROCEDURE copy_3d_r4
     MODULE PROCEDURE copy_1d_r8
     MODULE PROCEDURE copy_2d_r8
     MODULE PROCEDURE copy_3d_r8
  END INTERFACE copy

  INTERFACE scal
     MODULE PROCEDURE scal_1d_r4
     MODULE PROCEDURE scal_2d_r4
     MODULE PROCEDURE scal_3d_r4
     MODULE PROCEDURE scal_1d_r8
     MODULE PROCEDURE scal_2d_r8
     MODULE PROCEDURE scal_3d_r8
  END INTERFACE scal

  INTERFACE axpy
     MODULE PROCEDURE axpy_1d_r4
     MODULE PROCEDURE axpy_2d_r4
     MODULE PROCEDURE axpy_3d_r4
     MODULE PROCEDURE axpy_1d_r8
     MODULE PROCEDURE axpy_2d_r8
     MODULE PROCEDURE axpy_3d_r8
  END INTERFACE axpy

  INTERFACE increase_buffer
     MODULE PROCEDURE increase_buffer_i4
     MODULE PROCEDURE increase_buffer_i8
     MODULE PROCEDURE increase_buffer_r4
     MODULE PROCEDURE increase_buffer_r8
     MODULE PROCEDURE increase_buffer_i4_2d
     MODULE PROCEDURE increase_buffer_i8_2d
     MODULE PROCEDURE increase_buffer_r4_2d
     MODULE PROCEDURE increase_buffer_r8_2d
  END INTERFACE increase_buffer

CONTAINS
!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE assert(exp, msg, msg2)
    LOGICAL, INTENT(IN) :: exp
    CHARACTER(*), INTENT(IN) :: msg
    CHARACTER(*), INTENT(IN), OPTIONAL :: msg2

    IF (.NOT. exp) THEN
       IF (REPORT_UNIT > 0 .AND. REPORT_UNIT /= STDOUT_UNIT) THEN
          WRITE(REPORT_UNIT, *) "ERROR: " // msg
          IF (PRESENT(msg2)) WRITE(REPORT_UNIT, *) trim(msg2)
       END IF

       WRITE(STDERR_UNIT, *) "ERROR: " // trim(msg)
       IF (PRESENT(msg2)) WRITE(STDERR_UNIT, *) trim(msg2)

       STOP 1
    END IF

  END SUBROUTINE assert

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE warning(exp, msg)
    LOGICAL, INTENT(IN) :: exp
    CHARACTER(*), INTENT(IN) :: msg

    IF (.NOT. exp) THEN
       WRITE(STDERR_UNIT, *) msg
    END IF

  END SUBROUTINE warning

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE scaleoffset_3d_r4(data, scale, offset)
    REAL(4), INTENT(INOUT) :: data(:,:,:)
    REAL(4), INTENT(IN)    :: scale
    REAL(4), INTENT(IN)    :: offset

    INTEGER :: i, j, k

    IF (scale==1.0 .AND. offset==0.0) RETURN

!$OMP PARALLEL DO
    DO k=1, size(data,3)
    DO j=1, size(data,2)
    DO i=1, size(data,1)
       IF (data(i,j,k) == UNDEF) CYCLE
       data(i,j,k) = data(i,j,k)*scale + offset
    END DO
    END DO
    END DO
  END SUBROUTINE scaleoffset_3d_r4

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE scaleoffset_3d_r8(data, scale, offset)
    REAL(8), INTENT(INOUT) :: data(:,:,:)
    REAL(4), INTENT(IN)    :: scale
    REAL(4), INTENT(IN)    :: offset

    INTEGER :: i, j, k

    IF (scale==1.0 .AND. offset==0.0) RETURN

!$OMP PARALLEL DO
    DO k=1, size(data,3)
    DO j=1, size(data,2)
    DO i=1, size(data,1)
       IF (data(i,j,k) == UNDEF) CYCLE
       data(i,j,k) = data(i,j,k)*scale + offset
    END DO
    END DO
    END DO
  END SUBROUTINE scaleoffset_3d_r8

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE scaleoffset_2d_r4(data, scale, offset)
    REAL(4), INTENT(INOUT) :: data(:,:)
    REAL(4), INTENT(IN)    :: scale
    REAL(4), INTENT(IN)    :: offset

    INTEGER :: i, j

    IF (scale==1.0 .AND. offset==0.0) RETURN

!$OMP PARALLEL DO
    DO j=1, size(data,2)
    DO i=1, size(data,1)
       IF (data(i,j) == UNDEF) CYCLE
       data(i,j) = data(i,j)*scale + offset
    END DO
    END DO
  END SUBROUTINE scaleoffset_2d_r4

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE scaleoffset_2d_r8(data, scale, offset)
    REAL(8), INTENT(INOUT) :: data(:,:)
    REAL(4), INTENT(IN)    :: scale
    REAL(4), INTENT(IN)    :: offset

    INTEGER :: i, j

    IF (scale==1.0 .AND. offset==0.0) RETURN

!$OMP PARALLEL DO
    DO j=1, size(data,2)
    DO i=1, size(data,1)
       IF (data(i,j) == UNDEF) CYCLE
       data(i,j) = data(i,j)*scale + offset
    END DO
    END DO
  END SUBROUTINE scaleoffset_2d_r8

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE scaleoffset_1d_r4(data, scale, offset)
    REAL(4), INTENT(INOUT) :: data(:)
    REAL(4), INTENT(IN)    :: scale
    REAL(4), INTENT(IN)    :: offset

    INTEGER :: i

    IF (scale==1.0 .AND. offset==0.0) RETURN

!$OMP PARALLEL DO
    DO i=1, size(data)
       IF (data(i) == UNDEF) CYCLE
       data(i) = data(i)*scale + offset
    END DO
  END SUBROUTINE scaleoffset_1d_r4

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE scaleoffset_1d_r8(data, scale, offset)
    REAL(8), INTENT(INOUT) :: data(:)
    REAL(4), INTENT(IN)    :: scale
    REAL(4), INTENT(IN)    :: offset

    INTEGER :: i

    IF (scale==1.0 .AND. offset==0.0) RETURN

!$OMP PARALLEL DO
    DO i=1, size(data)
       IF (data(i) == UNDEF) CYCLE
       data(i) = data(i)*scale + offset
    END DO
  END SUBROUTINE scaleoffset_1d_r8

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE limitminmax_3d_r4(data, minlimit, maxlimit)
    REAL(4), INTENT(INOUT) :: data(:,:,:)
    REAL(4), INTENT(IN)    :: minlimit
    REAL(4), INTENT(IN)    :: maxlimit

    INTEGER :: i, j, k

    IF (minlimit /= UNDEF) THEN
!$OMP PARALLEL DO
       DO k=1, size(data,3)
       DO j=1, size(data,2)
       DO i=1, size(data,1)
          IF (data(i,j,k) == UNDEF) CYCLE
          data(i,j,k) = max(data(i,j,k), minlimit)
       END DO
       END DO
       END DO
    END IF

    IF (maxlimit /= UNDEF) THEN
!$OMP PARALLEL DO
       DO k=1, size(data,3)
       DO j=1, size(data,2)
       DO i=1, size(data,1)
          IF (data(i,j,k) == UNDEF) CYCLE
          data(i,j,k) = min(data(i,j,k), maxlimit)
       END DO
       END DO
       END DO
    END IF
  END SUBROUTINE limitminmax_3d_r4

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE limitminmax_3d_r8(data, minlimit, maxlimit)
    REAL(8), INTENT(INOUT) :: data(:,:,:)
    REAL(4), INTENT(IN)    :: minlimit
    REAL(4), INTENT(IN)    :: maxlimit

    INTEGER :: i, j, k

    IF (minlimit /= UNDEF) THEN
!$OMP PARALLEL DO
       DO k=1, size(data,3)
       DO j=1, size(data,2)
       DO i=1, size(data,1)
          IF (data(i,j,k) == UNDEF) CYCLE
          data(i,j,k) = max(data(i,j,k), REAL(minlimit,8))
       END DO
       END DO
       END DO
    END IF

    IF (maxlimit /= UNDEF) THEN
!$OMP PARALLEL DO
       DO k=1, size(data,3)
       DO j=1, size(data,2)
       DO i=1, size(data,1)
          IF (data(i,j,k) == UNDEF) CYCLE
          data(i,j,k) = min(data(i,j,k), REAL(maxlimit,8))
       END DO
       END DO
       END DO
    END IF
  END SUBROUTINE limitminmax_3d_r8

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE limitminmax_2d_r4(data, minlimit, maxlimit)
    REAL(4), INTENT(INOUT) :: data(:,:)
    REAL(4), INTENT(IN)    :: minlimit
    REAL(4), INTENT(IN)    :: maxlimit

    INTEGER :: i, j

    IF (minlimit /= UNDEF) THEN
!$OMP PARALLEL DO
       DO j=1, size(data,2)
       DO i=1, size(data,1)
          IF (data(i,j) == UNDEF) CYCLE
          data(i,j) = max(data(i,j), minlimit)
       END DO
       END DO
    END IF

    IF (maxlimit /= UNDEF) THEN
!$OMP PARALLEL DO
       DO j=1, size(data,2)
       DO i=1, size(data,1)
          IF (data(i,j) == UNDEF) CYCLE
          data(i,j) = min(data(i,j), maxlimit)
       END DO
       END DO
    END IF
  END SUBROUTINE limitminmax_2d_r4

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE limitminmax_2d_r8(data, minlimit, maxlimit)
    REAL(8), INTENT(INOUT) :: data(:,:)
    REAL(4), INTENT(IN)    :: minlimit
    REAL(4), INTENT(IN)    :: maxlimit

    INTEGER :: i, j

    IF (minlimit /= UNDEF) THEN
!$OMP PARALLEL DO
       DO j=1, size(data,2)
       DO i=1, size(data,1)
          IF (data(i,j) == UNDEF) CYCLE
          data(i,j) = max(data(i,j), REAL(minlimit,8))
       END DO
       END DO
    END IF

    IF (maxlimit /= UNDEF) THEN
!$OMP PARALLEL DO
       DO j=1, size(data,2)
       DO i=1, size(data,1)
          IF (data(i,j) == UNDEF) CYCLE
          data(i,j) = min(data(i,j), REAL(maxlimit,8))
       END DO
       END DO
    END IF
  END SUBROUTINE limitminmax_2d_r8

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE limitminmax_1d_r4(data, minlimit, maxlimit)
    REAL(4), INTENT(INOUT) :: data(:)
    REAL(4), INTENT(IN)    :: minlimit
    REAL(4), INTENT(IN)    :: maxlimit

    INTEGER :: i

    IF (minlimit /= UNDEF) THEN
!$OMP PARALLEL DO
       DO i=1, size(data,1)
          IF (data(i) == UNDEF) CYCLE
          data(i) = max(data(i), minlimit)
       END DO
    END IF

    IF (maxlimit /= UNDEF) THEN
!$OMP PARALLEL DO
       DO i=1, size(data,1)
          IF (data(i) == UNDEF) CYCLE
          data(i) = min(data(i), maxlimit)
       END DO
    END IF
  END SUBROUTINE limitminmax_1d_r4

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE limitminmax_1d_r8(data, minlimit, maxlimit)
    REAL(8), INTENT(INOUT) :: data(:)
    REAL(4), INTENT(IN)    :: minlimit
    REAL(4), INTENT(IN)    :: maxlimit

    INTEGER :: i

    IF (minlimit /= UNDEF) THEN
!$OMP PARALLEL DO
       DO i=1, size(data,1)
          IF (data(i) == UNDEF) CYCLE
          data(i) = max(data(i), REAL(minlimit,8))
       END DO
    END IF

    IF (maxlimit /= UNDEF) THEN
!$OMP PARALLEL DO
       DO i=1, size(data,1)
          IF (data(i) == UNDEF) CYCLE
          data(i) = min(data(i), REAL(maxlimit,8))
       END DO
    END IF
  END SUBROUTINE limitminmax_1d_r8

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE substmissing_3d_r4(data, missing, subst)
    REAL(4), INTENT(INOUT) :: data(:,:,:)
    REAL(4), INTENT(IN) :: missing
    REAL(4), INTENT(IN) :: subst

    INTEGER :: i, j, k

    IF (missing==subst) RETURN

!$OMP PARALLEL DO
    DO k=1, size(data,3)
    DO j=1, size(data,2)
    DO i=1, size(data,1)
       IF (data(i,j,k) == missing) data(i,j,k) = subst
    END DO
    END DO
    END DO
  END SUBROUTINE substmissing_3d_r4

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE substmissing_3d_r8(data, missing, subst)
    REAL(8), INTENT(INOUT) :: data(:,:,:)
    REAL(4), INTENT(IN) :: missing
    REAL(4), INTENT(IN) :: subst

    INTEGER :: i, j, k

    IF (missing==subst) RETURN

!$OMP PARALLEL DO
    DO k=1, size(data,3)
    DO j=1, size(data,2)
    DO i=1, size(data,1)
       IF (data(i,j,k) == missing) data(i,j,k) = subst
    END DO
    END DO
    END DO
  END SUBROUTINE substmissing_3d_r8

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE substmissing_2d_r4(data, missing, subst)
    REAL(4), INTENT(INOUT) :: data(:,:)
    REAL(4), INTENT(IN) :: missing
    REAL(4), INTENT(IN) :: subst

    INTEGER :: i, j

    IF (missing==subst) RETURN

!$OMP PARALLEL DO
    DO j=1, size(data,2)
    DO i=1, size(data,1)
       IF (data(i,j) == missing) data(i,j) = subst
    END DO
    END DO
  END SUBROUTINE substmissing_2d_r4

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE substmissing_2d_r8(data, missing, subst)
    REAL(8), INTENT(INOUT) :: data(:,:)
    REAL(4), INTENT(IN) :: missing
    REAL(4), INTENT(IN) :: subst

    INTEGER :: i, j

    IF (missing==subst) RETURN

!$OMP PARALLEL DO
    DO j=1, size(data,2)
    DO i=1, size(data,1)
       IF (data(i,j) == missing) data(i,j) = subst
    END DO
    END DO
  END SUBROUTINE substmissing_2d_r8

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE substmissing_1d_r4(data, missing, subst)
    REAL(4), INTENT(INOUT) :: data(:)
    REAL(4), INTENT(IN) :: missing
    REAL(4), INTENT(IN) :: subst

    INTEGER :: i

    IF (missing==subst) RETURN

!$OMP PARALLEL DO
    DO i=1, size(data)
       IF (data(i) == missing) data(i) = subst
    END DO
  END SUBROUTINE substmissing_1d_r4

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE substmissing_1d_r8(data, missing, subst)
    REAL(8), INTENT(INOUT) :: data(:)
    REAL(4), INTENT(IN) :: missing
    REAL(4), INTENT(IN) :: subst

    INTEGER :: i

    IF (missing==subst) RETURN

!$OMP PARALLEL DO
    DO i=1, size(data)
       IF (data(i) == missing) data(i) = subst
    END DO

  END SUBROUTINE substmissing_1d_r8

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE whitenoise_3d_r4(data, amp, seed)
    REAL(4), INTENT(INOUT) :: data(:,:,:)
    REAL(4), INTENT(IN) :: amp
    INTEGER, INTENT(IN) :: seed

    INTEGER :: i, j, k, s
    REAL(4), PARAMETER :: base = 0.5**31

    IF (amp==0.0) RETURN

    s = seed
    DO k=1, size(data,3)
    DO j=1, size(data,2)
    DO i=1, size(data,1)
       s = xorshift32(s)
       IF (data(i,j,k) /= UNDEF) data(i,j,k) = data(i,j,k) + amp*(s*base)
    END DO
    END DO
    END DO

  END SUBROUTINE whitenoise_3d_r4

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE whitenoise_3d_r8(data, amp, seed)
    REAL(8), INTENT(INOUT) :: data(:,:,:)
    REAL(4), INTENT(IN) :: amp
    INTEGER, INTENT(IN) :: seed

    INTEGER :: i, j, k, s
    REAL(4), PARAMETER :: base = 0.5**31

    IF (amp==0.0) RETURN

    s = seed
    DO k=1, size(data,3)
    DO j=1, size(data,2)
    DO i=1, size(data,1)
       s = xorshift32(s)
       IF (data(i,j,k) /= UNDEF) data(i,j,k) = data(i,j,k) + amp*(s*base)
    END DO
    END DO
    END DO

  END SUBROUTINE whitenoise_3d_r8

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE whitenoise_2d_r4(data, amp, seed)
    REAL(4), INTENT(INOUT) :: data(:,:)
    REAL(4), INTENT(IN) :: amp
    INTEGER, INTENT(IN) :: seed

    INTEGER :: i, j, s
    REAL(4), PARAMETER :: base = 0.5**31

    IF (amp==0.0) RETURN

    s = seed
    DO j=1, size(data,2)
    DO i=1, size(data,1)
       s = xorshift32(s)
       IF (data(i,j) /= UNDEF) data(i,j) = data(i,j) + amp*(s*base)
    END DO
    END DO

  END SUBROUTINE whitenoise_2d_r4

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE whitenoise_2d_r8(data, amp, seed)
    REAL(8), INTENT(INOUT) :: data(:,:)
    REAL(4), INTENT(IN) :: amp
    INTEGER, INTENT(IN) :: seed

    INTEGER :: i, j, s
    REAL(4), PARAMETER :: base = 0.5**31

    IF (amp==0.0) RETURN

    s = seed
    DO j=1, size(data,2)
    DO i=1, size(data,1)
       s = xorshift32(s)
       IF (data(i,j) /= UNDEF) data(i,j) = data(i,j) + amp*(s*base)
    END DO
    END DO

  END SUBROUTINE whitenoise_2d_r8

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE whitenoise_1d_r4(data, amp, seed)
    REAL(4), INTENT(INOUT) :: data(:)
    REAL(4), INTENT(IN) :: amp
    INTEGER, INTENT(IN) :: seed

    INTEGER :: i, s
    REAL(4), PARAMETER :: base = 0.5**31

    IF (amp==0.0) RETURN

    s = seed
    DO i=1, size(data)
       s = xorshift32(s)
       IF (data(i) /= UNDEF) data(i) = data(i) + amp*(s*base)
    END DO

  END SUBROUTINE whitenoise_1d_r4

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE whitenoise_1d_r8(data, amp, seed)
    REAL(8), INTENT(INOUT) :: data(:)
    REAL(4), INTENT(IN) :: amp
    INTEGER, INTENT(IN) :: seed

    INTEGER :: i, s
    REAL(4), PARAMETER :: base = 0.5**31

    IF (amp==0.0) RETURN

    s = seed
    DO i=1, size(data)
       s = xorshift32(s)
       IF (data(i) /= UNDEF) data(i) = data(i) + amp*(s*base)
    END DO

  END SUBROUTINE whitenoise_1d_r8

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE write_file_3d_i1(data, filename, kind)
    INTEGER(1),   INTENT(IN) :: data(:,:,:)
    CHARACTER(*), INTENT(IN) :: filename
    INTEGER, INTENT(IN), OPTIONAL :: kind

    INTEGER :: dim1
    INTEGER :: dim2
    INTEGER :: dim3

    INTEGER :: kind_
    INTEGER :: iostat

    INTEGER :: k

    dim1=size(data,1)
    dim2=size(data,2)
    dim3=size(data,3)

    kind_ = 1
    IF (present(kind)) kind_ = kind

    CALL assert(kind_==1 .OR. kind_==4 .OR. kind_==8, "unsupported KIND in WRITE_FILE_3D for FILENAME='"//trim(filename)//"'")

    OPEN(UNIT   = TMP_UNIT,        &
         FILE   = trim(filename),  &
         FORM   = 'UNFORMATTED',   &
         ACCESS = 'DIRECT',        &
         STATUS = 'REPLACE',       &
         ACTION = 'WRITE',         &
         RECL   = dim1*dim2*kind_, &
         IOSTAT = iostat)

    CALL assert(iostat==0, "failed to create '"//trim(filename)//"'")

    DO k=1, dim3
       SELECT CASE (kind_)
       CASE (1)
          WRITE(TMP_UNIT, REC=k) INT( data(:,:,k), KIND=1)
       CASE (4)
          WRITE(TMP_UNIT, REC=k) REAL(data(:,:,k), KIND=4)
       CASE (8)
          WRITE(TMP_UNIT, REC=k) REAL(data(:,:,k), KIND=8)
       END SELECT
    END DO

    CLOSE(TMP_UNIT)

  END SUBROUTINE write_file_3d_i1

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE write_file_3d_r4(data, filename, kind)
    REAL(4),      INTENT(IN) :: data(:,:,:)
    CHARACTER(*), INTENT(IN) :: filename
    INTEGER, INTENT(IN), OPTIONAL :: kind

    INTEGER :: dim1
    INTEGER :: dim2
    INTEGER :: dim3

    INTEGER :: kind_
    INTEGER :: iostat

    INTEGER :: k

    dim1=size(data,1)
    dim2=size(data,2)
    dim3=size(data,3)

    kind_ = 4
    IF (present(kind)) kind_ = kind

    CALL assert(kind_==1 .OR. kind_==4 .OR. kind_==8, "unsupported KIND in WRITE_FILE_3D for FILENAME='"//trim(filename)//"'")

    OPEN(UNIT   = TMP_UNIT,        &
         FILE   = trim(filename),  &
         FORM   = 'UNFORMATTED',   &
         ACCESS = 'DIRECT',        &
         STATUS = 'REPLACE',       &
         ACTION = 'WRITE',         &
         RECL   = dim1*dim2*kind_, &
         IOSTAT = iostat)

    CALL assert(iostat==0, "failed to create '"//trim(filename)//"'")

    DO k=1, dim3
       SELECT CASE (kind_)
       CASE (1)
          WRITE(TMP_UNIT, REC=k) INT( data(:,:,k), KIND=1)
       CASE (4)
          WRITE(TMP_UNIT, REC=k) REAL(data(:,:,k), KIND=4)
       CASE (8)
          WRITE(TMP_UNIT, REC=k) REAL(data(:,:,k), KIND=8)
       END SELECT
    END DO

    CLOSE(TMP_UNIT)

  END SUBROUTINE write_file_3d_r4

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE write_file_3d_r8(data, filename, kind)
    REAL(8),      INTENT(IN) :: data(:,:,:)
    CHARACTER(*), INTENT(IN) :: filename
    INTEGER, INTENT(IN), OPTIONAL :: kind

    INTEGER :: dim1
    INTEGER :: dim2
    INTEGER :: dim3

    INTEGER :: kind_
    INTEGER :: iostat

    INTEGER :: k

    dim1=size(data,1)
    dim2=size(data,2)
    dim3=size(data,3)

    kind_ = 8
    IF (present(kind)) kind_ = kind

    CALL assert(kind_==1 .OR. kind_==4 .OR. kind_==8, "unsupported KIND in WRITE_FILE_3D for FILENAME='"//trim(filename)//"'")

    OPEN(UNIT   = TMP_UNIT,        &
         FILE   = trim(filename),  &
         FORM   = 'UNFORMATTED',   &
         ACCESS = 'DIRECT',        &
         STATUS = 'REPLACE',       &
         ACTION = 'WRITE',         &
         RECL   = dim1*dim2*kind_, &
         IOSTAT = iostat)

    CALL assert(iostat==0, "failed to create '"//trim(filename)//"'")

    DO k=1, dim3
       SELECT CASE (kind_)
       CASE (1)
          WRITE(TMP_UNIT, REC=k) INT( data(:,:,k), KIND=1)
       CASE (4)
          WRITE(TMP_UNIT, REC=k) REAL(data(:,:,k), KIND=4)
       CASE (8)
          WRITE(TMP_UNIT, REC=k) REAL(data(:,:,k), KIND=8)
       END SELECT
    END DO

    CLOSE(TMP_UNIT)

  END SUBROUTINE write_file_3d_r8

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE write_file_2d_i1(data, filename, kind)
    INTEGER(1),   INTENT(IN) :: data(:,:)
    CHARACTER(*), INTENT(IN) :: filename
    INTEGER, INTENT(IN), OPTIONAL :: kind

    INTEGER :: dim1
    INTEGER :: dim2

    INTEGER :: kind_
    INTEGER :: iostat

    dim1=size(data,1)
    dim2=size(data,2)

    kind_ = 1
    IF (present(kind)) kind_ = kind

    CALL assert(kind_==1 .OR. kind_==4 .OR. kind_==8, "unsupported KIND in WRITE_FILE_2D for FILENAME='"//trim(filename)//"'")

    OPEN(UNIT   = TMP_UNIT,        &
         FILE   = trim(filename),  &
         FORM   = 'UNFORMATTED',   &
         ACCESS = 'DIRECT',        &
         STATUS = 'REPLACE',       &
         ACTION = 'WRITE',         &
         RECL   = dim1*dim2*kind_, &
         IOSTAT = iostat)

    CALL assert(iostat==0, "failed to create '"//trim(filename)//"'")

    SELECT CASE (kind_)
    CASE (1)
       WRITE(TMP_UNIT, REC=1) INT( data(:,:), KIND=1)
    CASE (4)
       WRITE(TMP_UNIT, REC=1) REAL(data(:,:), KIND=4)
    CASE (8)
       WRITE(TMP_UNIT, REC=1) REAL(data(:,:), KIND=8)
    END SELECT

    CLOSE(TMP_UNIT)

  END SUBROUTINE write_file_2d_i1

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE write_file_2d_r4(data, filename, kind)
    REAL(4),      INTENT(IN) :: data(:,:)
    CHARACTER(*), INTENT(IN) :: filename
    INTEGER, INTENT(IN), OPTIONAL :: kind

    INTEGER :: dim1
    INTEGER :: dim2

    INTEGER :: kind_
    INTEGER :: iostat

    dim1=size(data,1)
    dim2=size(data,2)

    kind_ = 4
    IF (present(kind)) kind_ = kind

    CALL assert(kind_==1 .OR. kind_==4 .OR. kind_==8, "unsupported KIND in WRITE_FILE_2D for FILENAME='"//trim(filename)//"'")

    OPEN(UNIT   = TMP_UNIT,        &
         FILE   = trim(filename),  &
         FORM   = 'UNFORMATTED',   &
         ACCESS = 'DIRECT',        &
         STATUS = 'REPLACE',       &
         ACTION = 'WRITE',         &
         RECL   = dim1*dim2*kind_, &
         IOSTAT = iostat)

    CALL assert(iostat==0, "failed to create '"//trim(filename)//"'")

    SELECT CASE (kind_)
    CASE (1)
       WRITE(TMP_UNIT, REC=1) INT( data(:,:), KIND=1)
    CASE (4)
       WRITE(TMP_UNIT, REC=1) REAL(data(:,:), KIND=4)
    CASE (8)
       WRITE(TMP_UNIT, REC=1) REAL(data(:,:), KIND=8)
    END SELECT

    CLOSE(TMP_UNIT)

  END SUBROUTINE write_file_2d_r4

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE write_file_2d_r8(data, filename, kind)
    REAL(8),      INTENT(IN) :: data(:,:)
    CHARACTER(*), INTENT(IN) :: filename
    INTEGER, INTENT(IN), OPTIONAL :: kind

    INTEGER :: dim1
    INTEGER :: dim2

    INTEGER :: kind_
    INTEGER :: iostat

    dim1=size(data,1)
    dim2=size(data,2)

    kind_ = 8
    IF (present(kind)) kind_ = kind

    CALL assert(kind_==1 .OR. kind_==4 .OR. kind_==8, "unsupported KIND in WRITE_FILE_2D for FILENAME='"//trim(filename)//"'")

    OPEN(UNIT   = TMP_UNIT,        &
         FILE   = trim(filename),  &
         FORM   = 'UNFORMATTED',   &
         ACCESS = 'DIRECT',        &
         STATUS = 'REPLACE',       &
         ACTION = 'WRITE',         &
         RECL   = dim1*dim2*kind_, &
         IOSTAT = iostat)

    CALL assert(iostat==0, "failed to create '"//trim(filename)//"'")

    SELECT CASE (kind_)
    CASE (1)
       WRITE(TMP_UNIT, REC=1) INT( data(:,:), KIND=1)
    CASE (4)
       WRITE(TMP_UNIT, REC=1) REAL(data(:,:), KIND=4)
    CASE (8)
       WRITE(TMP_UNIT, REC=1) REAL(data(:,:), KIND=8)
    END SELECT

    CLOSE(TMP_UNIT)

  END SUBROUTINE write_file_2d_r8

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE write_file_1d_i1(data, filename, kind)
    INTEGER(1),   INTENT(IN) :: data(:)
    CHARACTER(*), INTENT(IN) :: filename
    INTEGER, INTENT(IN), OPTIONAL :: kind

    INTEGER :: dim

    INTEGER :: kind_
    INTEGER :: iostat

    dim=size(data,1)

    kind_ = 1
    IF (present(kind)) kind_ = kind

    CALL assert(kind_==1 .OR. kind_==4 .OR. kind_==8, "unsupported KIND in WRITE_FILE_1D for FILENAME='"//trim(filename)//"'")

    OPEN(UNIT   = TMP_UNIT,      &
         FILE   = trim(filename),&
         FORM   = 'UNFORMATTED', &
         ACCESS = 'DIRECT',      &
         STATUS = 'REPLACE',     &
         ACTION = 'WRITE',       &
         RECL   = dim*kind_,     &
         IOSTAT = iostat)

    CALL assert(iostat==0, "failed to create '"//trim(filename)//"'")

    SELECT CASE (kind_)
    CASE (1)
       WRITE(TMP_UNIT, REC=1) INT( data(:), KIND=1)
    CASE (4)
       WRITE(TMP_UNIT, REC=1) REAL(data(:), KIND=4)
    CASE (8)
       WRITE(TMP_UNIT, REC=1) REAL(data(:), KIND=8)
    END SELECT

    CLOSE(TMP_UNIT)

  END SUBROUTINE write_file_1d_i1

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE write_file_1d_r4(data, filename, kind)
    REAL(4),      INTENT(IN) :: data(:)
    CHARACTER(*), INTENT(IN) :: filename
    INTEGER, INTENT(IN), OPTIONAL :: kind

    INTEGER :: dim

    INTEGER :: kind_
    INTEGER :: iostat

    dim=size(data,1)

    kind_ = 4
    IF (present(kind)) kind_ = kind

    CALL assert(kind_==4 .OR. kind_==8, "unsupported KIND in WRITE_FILE_1D for FILENAME='"//trim(filename)//"'")

    OPEN(UNIT   = TMP_UNIT,      &
         FILE   = trim(filename),&
         FORM   = 'UNFORMATTED', &
         ACCESS = 'DIRECT',      &
         STATUS = 'REPLACE',     &
         ACTION = 'WRITE',       &
         RECL   = dim*kind_,     &
         IOSTAT = iostat)

    CALL assert(iostat==0, "failed to create '"//trim(filename)//"'")

    SELECT CASE (kind_)
    CASE (1)
       WRITE(TMP_UNIT, REC=1) INT( data(:), KIND=1)
    CASE (4)
       WRITE(TMP_UNIT, REC=1) REAL(data(:), KIND=4)
    CASE (8)
       WRITE(TMP_UNIT, REC=1) REAL(data(:), KIND=8)
    END SELECT

    CLOSE(TMP_UNIT)

  END SUBROUTINE write_file_1d_r4

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE write_file_1d_r8(data, filename, kind)
    REAL(8),      INTENT(IN) :: data(:)
    CHARACTER(*), INTENT(IN) :: filename
    INTEGER, INTENT(IN), OPTIONAL :: kind

    INTEGER :: dim

    INTEGER :: kind_
    INTEGER :: iostat

    dim=size(data,1)

    kind_ = 8
    IF (present(kind)) kind_ = kind

    CALL assert(kind_==1 .OR. kind_==4 .OR. kind_==8, "unsupported KIND in WRITE_FILE_1D for FILENAME='"//trim(filename)//"'")

    OPEN(UNIT   = TMP_UNIT,      &
         FILE   = trim(filename),&
         FORM   = 'UNFORMATTED', &
         ACCESS = 'DIRECT',      &
         STATUS = 'REPLACE',     &
         ACTION = 'WRITE',       &
         RECL   = dim*kind_,     &
         IOSTAT = iostat)

    CALL assert(iostat==0, "failed to create '"//trim(filename)//"'")

    SELECT CASE (kind_)
    CASE (1)
       WRITE(TMP_UNIT, REC=1) INT( data(:), KIND=1)
    CASE (4)
       WRITE(TMP_UNIT, REC=1) REAL(data(:), KIND=4)
    CASE (8)
       WRITE(TMP_UNIT, REC=1) REAL(data(:), KIND=8)
    END SELECT

    CLOSE(TMP_UNIT)

  END SUBROUTINE write_file_1d_r8

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE read_file_3d_i1(data, filename, kind)
    INTEGER(1),   INTENT(OUT) :: data(:,:,:)
    CHARACTER(*), INTENT(IN)  :: filename
    INTEGER, INTENT(IN), OPTIONAL :: kind

    INTEGER :: dim1
    INTEGER :: dim2
    INTEGER :: dim3

    INTEGER :: kind_
    INTEGER :: iostat

    INTEGER(1) :: tmp1(size(data,1),size(data,2))
    REAL(4)    :: tmp4(size(data,1),size(data,2))
    REAL(8)    :: tmp8(size(data,1),size(data,2))

    INTEGER :: k

    dim1=size(data,1)
    dim2=size(data,2)
    dim3=size(data,3)

    kind_ = 1
    IF (present(kind)) kind_ = kind

    CALL assert(kind_==1 .OR. kind_==4 .OR. kind_==8, "unsupported KIND in READ_FILE_3D for FILENAME='"//trim(filename)//"'")

    OPEN(UNIT   = TMP_UNIT,        &
         FILE   = trim(filename),  &
         FORM   = 'UNFORMATTED',   &
         ACCESS = 'DIRECT',        &
         STATUS = 'OLD',           &
         ACTION = 'READ',          &
         RECL   = dim1*dim2*kind_, &
         IOSTAT = iostat)

    CALL assert(iostat==0, "failed to open '"//trim(filename)//"'")

    SELECT CASE (kind_)
    CASE (1)
       DO k=1, dim3
          READ(TMP_UNIT, REC=k) tmp1(:,:)
          data(:,:,k) = INT(tmp1(:,:), KIND=1)
       END DO
    CASE (4)
       DO k=1, dim3
          READ(TMP_UNIT, REC=k) tmp4(:,:)
          data(:,:,k) = INT(tmp4(:,:), KIND=1)
       END DO
    CASE (8)
       DO k=1, dim3
          READ(TMP_UNIT, REC=k) tmp8(:,:)
          data(:,:,k) = INT(tmp8(:,:), KIND=1)
       END DO
    END SELECT

    CLOSE(TMP_UNIT)

  END SUBROUTINE read_file_3d_i1

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE read_file_3d_r4(data, filename, kind)
    REAL(4),      INTENT(OUT) :: data(:,:,:)
    CHARACTER(*), INTENT(IN)  :: filename
    INTEGER, INTENT(IN), OPTIONAL :: kind

    INTEGER :: dim1
    INTEGER :: dim2
    INTEGER :: dim3

    INTEGER :: kind_
    INTEGER :: iostat

    INTEGER(1) :: tmp1(size(data,1),size(data,2))
    REAL(4)    :: tmp4(size(data,1),size(data,2))
    REAL(8)    :: tmp8(size(data,1),size(data,2))

    INTEGER :: k

    dim1=size(data,1)
    dim2=size(data,2)
    dim3=size(data,3)

    kind_ = 4
    IF (present(kind)) kind_ = kind

    CALL assert(kind_==1 .OR. kind_==4 .OR. kind_==8, "unsupported KIND in READ_FILE_3D for FILENAME='"//trim(filename)//"'")

    OPEN(UNIT   = TMP_UNIT,        &
         FILE   = trim(filename),  &
         FORM   = 'UNFORMATTED',   &
         ACCESS = 'DIRECT',        &
         STATUS = 'OLD',           &
         ACTION = 'READ',          &
         RECL   = dim1*dim2*kind_, &
         IOSTAT = iostat)

    CALL assert(iostat==0, "failed to open '"//trim(filename)//"'")

    SELECT CASE (kind_)
    CASE (1)
       DO k=1, dim3
          READ(TMP_UNIT, REC=k) tmp1(:,:)
          data(:,:,k) = REAL(tmp1(:,:), KIND=4)
       END DO
    CASE (4)
       DO k=1, dim3
          READ(TMP_UNIT, REC=k) tmp4(:,:)
          data(:,:,k) = REAL(tmp4(:,:), KIND=4)
       END DO
    CASE (8)
       DO k=1, dim3
          READ(TMP_UNIT, REC=k) tmp8(:,:)
          data(:,:,k) = REAL(tmp8(:,:), KIND=4)
       END DO
    END SELECT

    CLOSE(TMP_UNIT)

  END SUBROUTINE read_file_3d_r4

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE read_file_3d_r8(data, filename, kind)
    REAL(8),      INTENT(OUT) :: data(:,:,:)
    CHARACTER(*), INTENT(IN)  :: filename
    INTEGER, INTENT(IN), OPTIONAL :: kind

    INTEGER :: dim1
    INTEGER :: dim2
    INTEGER :: dim3

    INTEGER :: kind_
    INTEGER :: iostat

    INTEGER(1) :: tmp1(size(data,1),size(data,2))
    REAL(4)    :: tmp4(size(data,1),size(data,2))
    REAL(8)    :: tmp8(size(data,1),size(data,2))

    INTEGER :: k

    dim1=size(data,1)
    dim2=size(data,2)
    dim3=size(data,3)

    kind_ = 8
    IF (present(kind)) kind_ = kind

    CALL assert(kind_==1 .OR. kind_==4 .OR. kind_==8, "unsupported KIND in READ_FILE_3D for FILENAME='"//trim(filename)//"'")

    OPEN(UNIT   = TMP_UNIT,        &
         FILE   = trim(filename),  &
         FORM   = 'UNFORMATTED',   &
         ACCESS = 'DIRECT',        &
         STATUS = 'OLD',           &
         ACTION = 'READ',          &
         RECL   = dim1*dim2*kind_, &
         IOSTAT = iostat)

    CALL assert(iostat==0, "failed to open '"//trim(filename)//"'")

    SELECT CASE (kind_)
    CASE (1)
       DO k=1, dim3
          READ(TMP_UNIT, REC=k) tmp1(:,:)
          data(:,:,k) = REAL(tmp1(:,:), KIND=8)
       END DO
    CASE (4)
       DO k=1, dim3
          READ(TMP_UNIT, REC=k) tmp4(:,:)
          data(:,:,k) = REAL(tmp4(:,:), KIND=8)
       END DO
    CASE (8)
       DO k=1, dim3
          READ(TMP_UNIT, REC=k) tmp8(:,:)
          data(:,:,k) = REAL(tmp8(:,:), KIND=8)
       END DO
    END SELECT

    CLOSE(TMP_UNIT)

  END SUBROUTINE read_file_3d_r8

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE read_file_2d_i1(data, filename, kind)
    INTEGER(1),   INTENT(OUT) :: data(:,:)
    CHARACTER(*), INTENT(IN)  :: filename
    INTEGER, INTENT(IN), OPTIONAL :: kind

    INTEGER :: dim1
    INTEGER :: dim2

    INTEGER :: kind_
    INTEGER :: iostat

    INTEGER(1) :: tmp1(size(data,1),size(data,2))
    REAL(4)    :: tmp4(size(data,1),size(data,2))
    REAL(8)    :: tmp8(size(data,1),size(data,2))

    dim1=size(data,1)
    dim2=size(data,2)

    kind_ = 1
    IF (present(kind)) kind_ = kind

    CALL assert(kind_==1 .OR. kind_==4 .OR. kind_==8, "unsupported KIND in READ_FILE_2D for FILENAME='"//trim(filename)//"'")

    OPEN(UNIT   = TMP_UNIT,        &
         FILE   = trim(filename),  &
         FORM   = 'UNFORMATTED',   &
         ACCESS = 'DIRECT',        &
         STATUS = 'OLD',           &
         ACTION = 'READ',          &
         RECL   = dim1*dim2*kind_, &
         IOSTAT = iostat)

    CALL assert(iostat==0, "failed to open '"//trim(filename)//"'")

    SELECT CASE (kind_)
    CASE (1)
       READ(TMP_UNIT, REC=1) tmp1
       data(:,:) = INT(tmp1(:,:), KIND=1)
    CASE (4)
       READ(TMP_UNIT, REC=1) tmp4
       data(:,:) = INT(tmp4(:,:), KIND=1)
    CASE (8)
       READ(TMP_UNIT, REC=1) tmp8
       data(:,:) = INT(tmp8(:,:), KIND=1)
    END SELECT

    CLOSE(TMP_UNIT)

  END SUBROUTINE read_file_2d_i1

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE read_file_2d_r4(data, filename, kind)
    REAL(4),      INTENT(OUT) :: data(:,:)
    CHARACTER(*), INTENT(IN)  :: filename
    INTEGER, INTENT(IN), OPTIONAL :: kind

    INTEGER :: dim1
    INTEGER :: dim2

    INTEGER :: kind_
    INTEGER :: iostat

    INTEGER(1) :: tmp1(size(data,1),size(data,2))
    REAL(4)    :: tmp4(size(data,1),size(data,2))
    REAL(8)    :: tmp8(size(data,1),size(data,2))

    dim1=size(data,1)
    dim2=size(data,2)

    kind_ = 4
    IF (present(kind)) kind_ = kind

    CALL assert(kind_==1 .OR. kind_==4 .OR. kind_==8, "unsupported KIND in READ_FILE_2D for FILENAME='"//trim(filename)//"'")

    OPEN(UNIT   = TMP_UNIT,        &
         FILE   = trim(filename),  &
         FORM   = 'UNFORMATTED',   &
         ACCESS = 'DIRECT',        &
         STATUS = 'OLD',           &
         ACTION = 'READ',          &
         RECL   = dim1*dim2*kind_, &
         IOSTAT = iostat)

    CALL assert(iostat==0, "failed to open '"//trim(filename)//"'")

    SELECT CASE (kind_)
    CASE (1)
       READ(TMP_UNIT, REC=1) tmp1
       data(:,:) = REAL(tmp1(:,:), KIND=4)
    CASE (4)
       READ(TMP_UNIT, REC=1) tmp4
       data(:,:) = REAL(tmp4(:,:), KIND=4)
    CASE (8)
       READ(TMP_UNIT, REC=1) tmp8
       data(:,:) = REAL(tmp8(:,:), KIND=4)
    END SELECT

    CLOSE(TMP_UNIT)

  END SUBROUTINE read_file_2d_r4

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE read_file_2d_r8(data, filename, kind)
    REAL(8),      INTENT(OUT) :: data(:,:)
    CHARACTER(*), INTENT(IN)  :: filename
    INTEGER, INTENT(IN), OPTIONAL :: kind

    INTEGER :: dim1
    INTEGER :: dim2

    INTEGER :: kind_
    INTEGER :: iostat

    INTEGER(1) :: tmp1(size(data,1),size(data,2))
    REAL(4)    :: tmp4(size(data,1),size(data,2))
    REAL(8)    :: tmp8(size(data,1),size(data,2))

    dim1=size(data,1)
    dim2=size(data,2)

    kind_ = 8
    IF (present(kind)) kind_ = kind

    CALL assert(kind_==1 .OR. kind_==4 .OR. kind_==8, "unsupported KIND in READ_FILE_2D for FILENAME='"//trim(filename)//"'")

    OPEN(UNIT   = TMP_UNIT,        &
         FILE   = trim(filename),  &
         FORM   = 'UNFORMATTED',   &
         ACCESS = 'DIRECT',        &
         STATUS = 'OLD',           &
         ACTION = 'READ',          &
         RECL   = dim1*dim2*kind_, &
         IOSTAT = iostat)

    CALL assert(iostat==0, "failed to open '"//trim(filename)//"'")

    SELECT CASE (kind_)
    CASE (1)
       READ(TMP_UNIT, REC=1) tmp1
       data(:,:) = REAL(tmp1(:,:), KIND=8)
    CASE (4)
       READ(TMP_UNIT, REC=1) tmp4
       data(:,:) = REAL(tmp4(:,:), KIND=8)
    CASE (8)
       READ(TMP_UNIT, REC=1) tmp8
       data(:,:) = REAL(tmp8(:,:), KIND=8)
    END SELECT

    CLOSE(TMP_UNIT)

  END SUBROUTINE read_file_2d_r8

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE read_file_1d_i1(data, filename, kind)
    INTEGER(1),   INTENT(OUT) :: data(:)
    CHARACTER(*), INTENT(IN)  :: filename
    INTEGER, INTENT(IN), OPTIONAL :: kind

    INTEGER :: dim

    INTEGER :: kind_
    INTEGER :: iostat

    INTEGER(1) :: tmp1(size(data))
    REAL(4)    :: tmp4(size(data))
    REAL(8)    :: tmp8(size(data))

    dim=size(data)

    kind_ = 4
    IF (present(kind)) kind_ = kind

    CALL assert(kind_==1 .OR. kind_==4 .OR. kind_==8, "unsupported KIND in READ_FILE_1D for FILENAME='"//trim(filename)//"'")

    OPEN(UNIT   = TMP_UNIT,      &
         FILE   = trim(filename),&
         FORM   = 'UNFORMATTED', &
         ACCESS = 'DIRECT',      &
         STATUS = 'OLD',         &
         ACTION = 'READ',        &
         RECL   = dim*kind_,     &
         IOSTAT = iostat)

    CALL assert(iostat==0, "failed to open '"//trim(filename)//"'")

    SELECT CASE (kind_)
    CASE (1)
       READ(TMP_UNIT, REC=1) tmp1
       data(:) = INT(tmp1(:), KIND=1)
    CASE (4)
       READ(TMP_UNIT, REC=1) tmp4
       data(:) = INT(tmp4(:), KIND=1)
    CASE (8)
       READ(TMP_UNIT, REC=1) tmp8
       data(:) = INT(tmp8(:), KIND=1)
    END SELECT

    CLOSE(TMP_UNIT)

  END SUBROUTINE read_file_1d_i1

!-----------------------------------------------------------------------------------------------------------------------
  SUBROUTINE read_file_1d_r4(data, filename, kind)
    REAL(4),      INTENT(OUT) :: data(:)
    CHARACTER(*), INTENT(IN)  :: filename
    INTEGER, INTENT(IN), OPTIONAL :: kind

    INTEGER :: dim

    INTEGER :: kind_
    INTEGER :: iostat

    INTEGER(1) :: tmp1(size(data))
    REAL(4)    :: tmp4(size(data))
    REAL(8)    :: tmp8(size(data))

    dim=size(data)

    kind_ = 4
    IF (present(kind)) kind_ = kind

    CALL assert(kind_==1 .OR. kind_==4 .OR. kind_==8, "unsupported KIND in READ_FILE_1D for FILENAME='"//trim(filename)//"'")

    OPEN(UNIT   = TMP_UNIT,      &
         FILE   = trim(filename),&
         FORM   = 'UNFORMATTED', &
         ACCESS = 'DIRECT',      &
         STATUS = 'OLD',         &
         ACTION = 'READ',        &
         RECL   = dim*kind_,     &
         IOSTAT = iostat)

    CALL assert(iostat==0, "failed to open '"//trim(filename)//"'")

    SELECT CASE (kind_)
    CASE (1)
       READ(TMP_UNIT, REC=1) tmp1
       data(:) = REAL(tmp1(:), KIND=4)
    CASE (4)
       READ(TMP_UNIT, REC=1) tmp4
       data(:) = REAL(tmp4(:), KIND=4)
    CASE (8)
       READ(TMP_UNIT, REC=1) tmp8
       data(:) = REAL(tmp8(:), KIND=4)
    END SELECT

    CLOSE(TMP_UNIT)

  END SUBROUTINE read_file_1d_r4

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE read_file_1d_r8(data, filename, kind)
    REAL(8),      INTENT(OUT) :: data(:)
    CHARACTER(*), INTENT(IN)  :: filename
    INTEGER, INTENT(IN), OPTIONAL :: kind

    INTEGER :: dim

    INTEGER :: kind_
    INTEGER :: iostat

    INTEGER(1) :: tmp1(size(data))
    REAL(4)    :: tmp4(size(data))
    REAL(8)    :: tmp8(size(data))

    dim=size(data)

    kind_ = 8
    IF (present(kind)) kind_ = kind

    CALL assert(kind_==1 .OR. kind_==4 .OR. kind_==8, "unsupported KIND in READ_FILE_1D for FILENAME='"//trim(filename)//"'")

    OPEN(UNIT   = TMP_UNIT,      &
         FILE   = trim(filename),&
         FORM   = 'UNFORMATTED', &
         ACCESS = 'DIRECT',      &
         STATUS = 'OLD',         &
         ACTION = 'READ',        &
         RECL   = dim*kind_,     &
         IOSTAT = iostat)

    CALL assert(iostat==0, "failed to open '"//trim(filename)//"'")

    SELECT CASE (kind_)
    CASE (1)
       READ(TMP_UNIT, REC=1) tmp1
       data(:) = REAL(tmp1(:), KIND=8)
    CASE (4)
       READ(TMP_UNIT, REC=1) tmp4
       data(:) = REAL(tmp4(:), KIND=8)
    CASE (8)
       READ(TMP_UNIT, REC=1) tmp8
       data(:) = REAL(tmp8(:), KIND=8)
    END SELECT

    CLOSE(TMP_UNIT)

  END SUBROUTINE read_file_1d_r8

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE dump_file(x, filename, size, stride, kind)
    REAL(4),      INTENT(IN) :: x
    CHARACTER(*), INTENT(IN) :: filename
    INTEGER,      INTENT(IN) :: size
    INTEGER,      INTENT(IN), OPTIONAL :: stride
    INTEGER,      INTENT(IN), OPTIONAL :: kind

    REAL(4), ALLOCATABLE :: buf(:)

    INTEGER :: iostat
    INTEGER :: s, n, i
    INTEGER :: kind_

    kind_ = 4
    IF (present(kind)) kind_ = kind
    CALL assert(kind_==1 .OR. kind_==4 .OR. kind_==8, "usupported KIND in DUMP_FILE")

    IF (present(stride)) THEN
       CALL assert(mod(size,stride)==0, "STRIDE shoud be a divisor of SIZE")
       n = size/stride
       s = stride
    ELSE
       n = 1
       s = size
    END IF

    ALLOCATE(buf(s))

    OPEN(UNIT   = TMP_UNIT,        &
         FILE   = trim(filename),  &
         FORM   = 'UNFORMATTED',   &
         ACCESS = 'DIRECT',        &
         STATUS = 'REPLACE',       &
         ACTION = 'WRITE',         &
         RECL   = s*kind_,         &
         IOSTAT = iostat)

    CALL assert(iostat==0, "failed to create '"//trim(filename)//"'")

    buf(:) = x

    DO i=1, n
       SELECT CASE (kind_)
       CASE (1)
          WRITE(TMP_UNIT, REC=i) INT( buf, KIND=1)
       CASE (4)
          WRITE(TMP_UNIT, REC=i) REAL(buf, KIND=4)
       CASE (8)
          WRITE(TMP_UNIT, REC=i) REAL(buf, KIND=8)
       END SELECT
    END DO

    DEALLOCATE (buf)
    CLOSE(TMP_UNIT)

  END SUBROUTINE dump_file

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE kreverse_1d_r4(data)
    REAL(4), INTENT(INOUT) :: data(:)

    REAL(4) :: tmp
    INTEGER :: k, ksize

    ksize = size(data)
    DO k=1, ksize/2
       tmp = data(k)
       data(k) = data(ksize+1-k)
       data(ksize+1-k) = tmp
    END DO

  END SUBROUTINE kreverse_1d_r4

  SUBROUTINE kreverse_1d_r8(data)
    REAL(8), INTENT(INOUT) :: data(:)

    REAL(8) :: tmp
    INTEGER :: k, ksize

    ksize = size(data)
    DO k=1, ksize/2
       tmp = data(k)
       data(k) = data(ksize+1-k)
       data(ksize+1-k) = tmp
    END DO

  END SUBROUTINE kreverse_1d_r8

  SUBROUTINE kreverse_2d_r4(data)
    REAL(4), INTENT(INOUT) :: data(:,:)

    REAL(4) :: tmp(size(data,1))
    INTEGER :: k, ksize

    ksize = size(data,2)
    DO k=1, ksize/2
       tmp(:) = data(:,k)
       data(:,k) = data(:,ksize+1-k)
       data(:,ksize+1-k) = tmp(:)
    END DO

  END SUBROUTINE kreverse_2d_r4

  SUBROUTINE kreverse_2d_r8(data)
    REAL(8), INTENT(INOUT) :: data(:,:)

    REAL(8) :: tmp(size(data,1))
    INTEGER :: k, ksize

    ksize = size(data, 2)
    DO k=1, ksize/2
       tmp(:) = data(:,k)
       data(:,k) = data(:,ksize+1-k)
       data(:,ksize+1-k) = tmp(:)
    END DO

  END SUBROUTINE kreverse_2d_r8

  SUBROUTINE kreverse_3d_r4(data)
    REAL(4), INTENT(INOUT) :: data(:,:,:)

    REAL(4) :: tmp(size(data,1),size(data,2))
    INTEGER :: k, ksize

    ksize = size(data, 3)
    DO k=1, ksize/2
       tmp(:,:) = data(:,:,k)
       data(:,:,k) = data(:,:,ksize+1-k)
       data(:,:,ksize+1-k) = tmp(:,:)
    END DO

  END SUBROUTINE kreverse_3d_r4

  SUBROUTINE kreverse_3d_r8(data)
    REAL(8), INTENT(INOUT) :: data(:,:,:)

    REAL(8) :: tmp(size(data,1),size(data,2))
    INTEGER :: k, ksize

    ksize = size(data, 3)
    DO k=1, ksize/2
       tmp(:,:) = data(:,:,k)
       data(:,:,k) = data(:,:,ksize+1-k)
       data(:,:,ksize+1-k) = tmp(:,:)
    END DO

  END SUBROUTINE kreverse_3d_r8

!-----------------------------------------------------------------------------------------------------------------------

  RECURSIVE SUBROUTINE replace(a, b, str)
    CHARACTER(*), INTENT(IN) :: a
    CHARACTER(*), INTENT(IN) :: b
    CHARACTER(*), INTENT(INOUT) :: str

    INTEGER :: i, l, la, lb, ls
    la = len(a)
    lb = len(b)
    ls = len_trim(str)
    i = index(str, a)

    IF (i == 0) RETURN

    CALL assert(len(str) >= ls - la + lb, "str length is not enough to replace")

    str(i+lb:i+lb+len_trim(str(i+la:))-1) = trim(str(i+la:))
    str(i:i+lb-1) = b

    IF (la > lb) THEN
       DO i=ls+lb-la+1, ls
          str(i:i) = ' '
       END DO
    END IF

    CALL replace(a, b, str)

  END SUBROUTINE replace

!-----------------------------------------------------------------------------------------------------------------------

  CHARACTER(1024) FUNCTION basename(filename)
    CHARACTER(*), INTENT(IN) :: filename
    basename = trim(filename(index(filename, '/', back=.TRUE.)+1:))
  END FUNCTION basename

!-----------------------------------------------------------------------------------------------------------------------

  CHARACTER(1024) PURE FUNCTION path(dir, basename, suffix)
    CHARACTER(*), INTENT(IN) :: dir
    CHARACTER(*), INTENT(IN) :: basename
    CHARACTER(*), INTENT(IN), OPTIONAL :: suffix

    INTEGER :: i

    path = trim(basename)

    IF (path == '') RETURN

    IF (present(suffix)) THEN
       path = trim(path) // trim(suffix)
    END IF

    IF (path(1:1) /= '/' .AND. trim(dir) /= '') THEN
       i = len_trim(dir)
       IF (dir(i:i) /= '/') path = '/' // trim(path)
       path = trim(dir) // trim(path)
    END IF

  END FUNCTION path

!-----------------------------------------------------------------------------------------------------------------------

  CHARACTER(32) PURE FUNCTION format_int4(n, fmt)
    INTEGER(4),   INTENT(IN) :: n
    CHARACTER(*), INTENT(IN), OPTIONAL :: fmt

    IF (present(fmt)) THEN
       WRITE(format_int4, FMT='('//trim(fmt)//')') n
    ELSE
       WRITE(format_int4, FMT=*) n
    END IF
    format_int4 = adjustl(format_int4)
  END FUNCTION format_int4

!-----------------------------------------------------------------------------------------------------------------------

  CHARACTER(32) PURE FUNCTION format_int8(n, fmt)
    INTEGER(8),   INTENT(IN) :: n
    CHARACTER(*), INTENT(IN), OPTIONAL  :: fmt

    IF (present(fmt)) THEN
       WRITE(format_int8, FMT='('//trim(fmt)//')') n
    ELSE
       WRITE(format_int8, FMT=*) n
    END IF
    format_int8 = adjustl(format_int8)
  END FUNCTION format_int8

!-----------------------------------------------------------------------------------------------------------------------

  CHARACTER(32) PURE FUNCTION format_real4(x, fmt)
    REAL(4),      INTENT(IN) :: x
    CHARACTER(*), INTENT(IN), OPTIONAL  :: fmt

    IF (present(fmt)) THEN
       WRITE(format_real4, FMT='('//trim(fmt)//')') x
    ELSE
       WRITE(format_real4, FMT=*) x
    END IF
    format_real4 = adjustl(format_real4)
  END FUNCTION format_real4

!-----------------------------------------------------------------------------------------------------------------------

  CHARACTER(32) PURE FUNCTION format_real8(x, fmt)
    REAL(8),      INTENT(IN) :: x
    CHARACTER(*), INTENT(IN), OPTIONAL  :: fmt

    IF (present(fmt)) THEN
       WRITE(format_real8, FMT='('//trim(fmt)//')') x
    ELSE
       WRITE(format_real8, FMT=*) x
    END IF
    format_real8 = adjustl(format_real8)
  END FUNCTION format_real8

!-----------------------------------------------------------------------------------------------------------------------

  REAL(8) FUNCTION logistic(x)
    REAL(8), INTENT(IN) :: x

    logistic = 1.0D0 / (1.0D0 + exp(-x))
  END FUNCTION logistic

!-----------------------------------------------------------------------------------------------------------------------

  LOGICAL PURE FUNCTION check_endian()
    INTEGER(4) :: x
    INTEGER(1) :: y(4)

    x = 1_4
    y = transfer(x, y)

    check_endian = (y(4) == 1_1)
    ! .TRUE.  if big-endian
    ! .FALSE. if little-endian
  END FUNCTION check_endian

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE convert_endian(buffer, size, kind)
    INTEGER, INTENT(IN) :: size
    INTEGER, INTENT(IN) :: kind
    INTEGER(1), INTENT(INOUT) :: buffer(size*kind)

    INTEGER(1) :: tmp
    INTEGER :: i, j, k

    k = INT(kind/2)

!$OMP PARALLEL DO PRIVATE(tmp)
    DO i=1, size
       DO j=1, k
          tmp = buffer((i-1)*kind+j)
          buffer((i-1)*kind+j) = buffer(i*kind+1-j)
          buffer(i*kind+1-j) = tmp
       END DO
    END DO

  END SUBROUTINE convert_endian

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE nan_to_undef_1d_r8(data)
    REAL(8) :: data(:)
    INTEGER :: i

    DO i=1, size(data)
       IF (.NOT. data(i)==data(i)) data(i) = UNDEF
    END DO
  END SUBROUTINE nan_to_undef_1d_r8

  SUBROUTINE nan_to_undef_2d_r8(data)
    REAL(8) :: data(:,:)
    INTEGER :: i, j

    DO j=1, size(data,2)
    DO i=1, size(data,1)
       IF (.NOT. data(i,j)==data(i,j)) data(i,j) = UNDEF
    END DO
    END DO
  END SUBROUTINE nan_to_undef_2d_r8

  SUBROUTINE nan_to_undef_3d_r8(data)
    REAL(8) :: data(:,:,:)
    INTEGER :: i, j, k

    DO k=1, size(data,3)
    DO j=1, size(data,2)
    DO i=1, size(data,1)
       IF (.NOT. data(i,j,k)==data(i,j,k)) data(i,j,k) = UNDEF
    END DO
    END DO
    END DO
  END SUBROUTINE nan_to_undef_3d_r8

  SUBROUTINE nan_to_undef_1d_r4(data)
    REAL(4) :: data(:)
    INTEGER :: i

    DO i=1, size(data)
       IF (.NOT. data(i)==data(i)) data(i) = UNDEF
    END DO
  END SUBROUTINE nan_to_undef_1d_r4

  SUBROUTINE nan_to_undef_2d_r4(data)
    REAL(4) :: data(:,:)
    INTEGER :: i, j

    DO j=1, size(data,2)
    DO i=1, size(data,1)
       IF (.NOT. data(i,j)==data(i,j)) data(i,j) = UNDEF
    END DO
    END DO
  END SUBROUTINE nan_to_undef_2d_r4

  SUBROUTINE nan_to_undef_3d_r4(data)
    REAL(4) :: data(:,:,:)
    INTEGER :: i, j, k

    DO k=1, size(data,3)
    DO j=1, size(data,2)
    DO i=1, size(data,1)
       IF (.NOT. data(i,j,k)==data(i,j,k)) data(i,j,k) = UNDEF
    END DO
    END DO
    END DO
  END SUBROUTINE nan_to_undef_3d_r4

!-----------------------------------------------------------------------------------------------------------------------

  FUNCTION check_nan_1d_r8(data) RESULT(result)
#ifdef F2003
    USE, INTRINSIC :: ieee_arithmetic
#endif
    REAL(8), INTENT(IN) :: data(:)
    LOGICAL :: result
    INTEGER :: i

    result = .TRUE.
!$OMP PARALLEL DO REDUCTION (.AND. : result)
    DO i=1, size(data)
#ifdef F2003
       IF (ieee_is_nan(data(i)))     result = result .AND. .FALSE.
#else
       IF (.NOT. (data(i)==data(i))) result = result .AND. .FALSE.
#endif
    END DO
  END FUNCTION check_nan_1d_r8

  FUNCTION check_nan_1d_r4(data) RESULT(result)
#ifdef F2003
    USE, INTRINSIC :: ieee_arithmetic
#endif
    REAL(4), INTENT(IN) :: data(:)
    LOGICAL :: result
    INTEGER :: i

    result = .TRUE.
!$OMP PARALLEL DO REDUCTION (.AND. : result)
    DO i=1, size(data)
#ifdef F2003
       IF (ieee_is_nan(data(i)))     result = result .AND. .FALSE.
#else
       IF (.NOT. (data(i)==data(i))) result = result .AND. .FALSE.
#endif
    END DO
  END FUNCTION check_nan_1d_r4

  FUNCTION check_nan_2d_r8(data) RESULT(result)
#ifdef F2003
    USE, INTRINSIC :: ieee_arithmetic
#endif
    REAL(8), INTENT(IN) :: data(:,:)
    LOGICAL :: result
    INTEGER :: i, j

    result = .TRUE.
!$OMP PARALLEL DO REDUCTION (.AND. : result) COLLAPSE(2)
    DO j=1, size(data,2)
    DO i=1, size(data,1)
#ifdef F2003
       IF (ieee_is_nan(data(i,j)))       result = result .AND. .FALSE.
#else
       IF (.NOT. (data(i,j)==data(i,j))) result = result .AND. .FALSE.
#endif
    END DO
    END DO
  END FUNCTION check_nan_2d_r8

  FUNCTION check_nan_2d_r4(data) RESULT(result)
#ifdef F2003
    USE, INTRINSIC :: ieee_arithmetic
#endif
    REAL(4), INTENT(IN) :: data(:,:)
    LOGICAL :: result
    INTEGER :: i, j

    result = .TRUE.
!$OMP PARALLEL DO REDUCTION (.AND. : result) COLLAPSE(2)
    DO j=1, size(data,2)
    DO i=1, size(data,1)
#ifdef F2003
       IF (ieee_is_nan(data(i,j)))       result = result .AND. .FALSE.
#else
       IF (.NOT. (data(i,j)==data(i,j))) result = result .AND. .FALSE.
#endif
    END DO
    END DO
  END FUNCTION check_nan_2d_r4

  FUNCTION check_nan_3d_r8(data) RESULT(result)
#ifdef F2003
    USE, INTRINSIC :: ieee_arithmetic
#endif
    REAL(8), INTENT(IN) :: data(:,:,:)
    LOGICAL :: result
    INTEGER :: i, j, k

    result = .TRUE.
!$OMP PARALLEL DO REDUCTION (.AND. : result)
    DO k=1, size(data,3)
    DO j=1, size(data,2)
    DO i=1, size(data,1)
#ifdef F2003
       IF (ieee_is_nan(data(i,j,k)))         result = result .AND. .FALSE.
#else
       IF (.NOT. (data(i,j,k)==data(i,j,k))) result = result .AND. .FALSE.
#endif
    END DO
    END DO
    END DO
  END FUNCTION check_nan_3d_r8

  FUNCTION check_nan_3d_r4(data) RESULT(result)
#ifdef F2003
    USE, INTRINSIC :: ieee_arithmetic
#endif
    REAL(4), INTENT(IN) :: data(:,:,:)
    LOGICAL :: result
    INTEGER :: i, j, k

    result = .TRUE.
!$OMP PARALLEL DO REDUCTION (.AND. : result)
    DO k=1, size(data,3)
    DO j=1, size(data,2)
    DO i=1, size(data,1)
#ifdef F2003
       IF (ieee_is_nan(data(i,j,k)))         result = result .AND. .FALSE.
#else
       IF (.NOT. (data(i,j,k)==data(i,j,k))) result = result .AND. .FALSE.
#endif
    END DO
    END DO
    END DO
  END FUNCTION check_nan_3d_r4

!-----------------------------------------------------------------------------------------------------------------------

  LOGICAL PURE FUNCTION valid_char(a, allow)
    CHARACTER(1), INTENT(IN) :: a
    CHARACTER(*), INTENT(IN), OPTIONAL :: allow
    INTEGER :: i

    valid_char = (iachar(a) >= iachar('0') .AND. iachar(a) <= iachar('9')) .OR. &
                 (iachar(a) >= iachar('a') .AND. iachar(a) <= iachar('z')) .OR. &
                 (iachar(a) >= iachar('A') .AND. iachar(a) <= iachar('Z')) .OR. &
                 a == '+' .OR. a == '-' .OR. a == '_'

    IF (present(allow)) THEN
       DO i=1, len(allow)
          valid_char = valid_char .OR. a == allow(i:i)
       END DO
    END IF
  END FUNCTION valid_char

  LOGICAL PURE FUNCTION valid_name(a, allow)
    CHARACTER(*), INTENT(IN) :: a
    CHARACTER(*), INTENT(IN), OPTIONAL :: allow

    INTEGER :: i

    valid_name = .FALSE.

    DO i=1, len_trim(a)
       valid_name = valid_char(a(i:i), allow)

       IF (.NOT. valid_name) RETURN
    END DO

  END FUNCTION valid_name

!-----------------------------------------------------------------------------------------------------------------------

  LOGICAL FUNCTION valid_dir(dir)
    CHARACTER(*), INTENT(IN) :: dir
    LOGICAL :: tmp
    INTEGER :: stat

!    INQUIRE(file=trim(dir) // '/.', exist=tmp)
    OPEN(UNIT=TMP_UNIT, FILE=trim(dir) // '/.KINACO_DIRCHK', IOSTAT=stat)
    IF (stat==0) THEN
       valid_dir = .TRUE.
       CLOSE(TMP_UNIT)
    ELSE
       valid_dir = .FALSE.
    END IF

  END FUNCTION valid_dir

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE copy_1d_r4(x, y)
    REAL(4), INTENT(IN)  :: x(:)
    REAL(4), INTENT(OUT) :: y(:)
#ifdef F2008
    CONTIGUOUS x, y
#endif
!$OMP PARALLEL WORKSHARE
    y(:) = x(:)
!$OMP END PARALLEL WORKSHARE
  END SUBROUTINE copy_1d_r4

  SUBROUTINE copy_1d_r8(x, y)
    REAL(8), INTENT(IN)  :: x(:)
    REAL(8), INTENT(OUT) :: y(:)
#ifdef F2008
    CONTIGUOUS x, y
#endif
!$OMP PARALLEL WORKSHARE
    y(:) = x(:)
!$OMP END PARALLEL WORKSHARE
  END SUBROUTINE copy_1d_r8

  SUBROUTINE copy_2d_r4(x, y)
    REAL(4), INTENT(IN)  :: x(:,:)
    REAL(4), INTENT(OUT) :: y(:,:)
#ifdef F2008
    CONTIGUOUS x, y
#endif
!$OMP PARALLEL WORKSHARE
    y(:,:) = x(:,:)
!$OMP END PARALLEL WORKSHARE
  END SUBROUTINE copy_2d_r4

  SUBROUTINE copy_2d_r8(x, y)
    REAL(8), INTENT(IN)  :: x(:,:)
    REAL(8), INTENT(OUT) :: y(:,:)
#ifdef F2008
    CONTIGUOUS x, y
#endif
!$OMP PARALLEL WORKSHARE
    y(:,:) = x(:,:)
!$OMP END PARALLEL WORKSHARE
  END SUBROUTINE copy_2d_r8

  SUBROUTINE copy_3d_r4(x, y)
    REAL(4), INTENT(IN)  :: x(:,:,:)
    REAL(4), INTENT(OUT) :: y(:,:,:)
#ifdef F2008
    CONTIGUOUS x, y
#endif
!$OMP PARALLEL WORKSHARE
    y(:,:,:) = x(:,:,:)
!$OMP END PARALLEL WORKSHARE
  END SUBROUTINE copy_3d_r4

  SUBROUTINE copy_3d_r8(x, y)
    REAL(8), INTENT(IN)  :: x(:,:,:)
    REAL(8), INTENT(OUT) :: y(:,:,:)
#ifdef F2008
    CONTIGUOUS x, y
#endif
!$OMP PARALLEL WORKSHARE
    y(:,:,:) = x(:,:,:)
!$OMP END PARALLEL WORKSHARE
  END SUBROUTINE copy_3d_r8

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE scal_1d_r4(a, x, y)
    REAL(4), INTENT(IN)    :: a
    REAL(4), INTENT(INOUT) :: x(:)
    REAL(4), INTENT(OUT), OPTIONAL :: y(:)
#ifdef F2008
    CONTIGUOUS x, y
#endif
    IF (PRESENT(y)) THEN
!$OMP PARALLEL WORKSHARE
       y(:) = a*x(:)
!$OMP END PARALLEL WORKSHARE
    ELSE
!$OMP PARALLEL WORKSHARE
       x(:) = a*x(:)
!$OMP END PARALLEL WORKSHARE
    END IF
  END SUBROUTINE scal_1d_r4

  SUBROUTINE scal_1d_r8(a, x, y)
    REAL(8), INTENT(IN)    :: a
    REAL(8), INTENT(INOUT) :: x(:)
    REAL(8), INTENT(OUT), OPTIONAL :: y(:)
#ifdef F2008
    CONTIGUOUS x, y
#endif
    IF (PRESENT(y)) THEN
!$OMP PARALLEL WORKSHARE
       y(:) = a*x(:)
!$OMP END PARALLEL WORKSHARE
    ELSE
!$OMP PARALLEL WORKSHARE
       x(:) = a*x(:)
!$OMP END PARALLEL WORKSHARE
    END IF
  END SUBROUTINE scal_1d_r8

  SUBROUTINE scal_2d_r4(a, x, y)
    REAL(4), INTENT(IN)    :: a
    REAL(4), INTENT(INOUT) :: x(:,:)
    REAL(4), INTENT(OUT), OPTIONAL :: y(:,:)
#ifdef F2008
    CONTIGUOUS x, y
#endif
    IF (PRESENT(y)) THEN
!$OMP PARALLEL WORKSHARE
       y(:,:) = a*x(:,:)
!$OMP END PARALLEL WORKSHARE
    ELSE
!$OMP PARALLEL WORKSHARE
       x(:,:) = a*x(:,:)
!$OMP END PARALLEL WORKSHARE
    END IF
  END SUBROUTINE scal_2d_r4

  SUBROUTINE scal_2d_r8(a, x, y)
    REAL(8), INTENT(IN)    :: a
    REAL(8), INTENT(INOUT) :: x(:,:)
    REAL(8), INTENT(OUT), OPTIONAL :: y(:,:)
#ifdef F2008
    CONTIGUOUS x, y
#endif
    IF (PRESENT(y)) THEN
!$OMP PARALLEL WORKSHARE
       y(:,:) = a*x(:,:)
!$OMP END PARALLEL WORKSHARE
    ELSE
!$OMP PARALLEL WORKSHARE
       x(:,:) = a*x(:,:)
!$OMP END PARALLEL WORKSHARE
    END IF
  END SUBROUTINE scal_2d_r8

  SUBROUTINE scal_3d_r4(a, x, y)
    REAL(4), INTENT(IN)    :: a
    REAL(4), INTENT(INOUT) :: x(:,:,:)
    REAL(4), INTENT(OUT), OPTIONAL :: y(:,:,:)
#ifdef F2008
    CONTIGUOUS x, y
#endif
    IF (PRESENT(y)) THEN
!$OMP PARALLEL WORKSHARE
       y(:,:,:) = a*x(:,:,:)
!$OMP END PARALLEL WORKSHARE
    ELSE
!$OMP PARALLEL WORKSHARE
       x(:,:,:) = a*x(:,:,:)
!$OMP END PARALLEL WORKSHARE
    END IF
  END SUBROUTINE scal_3d_r4

  SUBROUTINE scal_3d_r8(a, x, y)
    REAL(8), INTENT(IN)    :: a
    REAL(8), INTENT(INOUT) :: x(:,:,:)
    REAL(8), INTENT(OUT), OPTIONAL :: y(:,:,:)
#ifdef F2008
    CONTIGUOUS x, y
#endif
    IF (PRESENT(y)) THEN
!$OMP PARALLEL WORKSHARE
       y(:,:,:) = a*x(:,:,:)
!$OMP END PARALLEL WORKSHARE
    ELSE
!$OMP PARALLEL WORKSHARE
       x(:,:,:) = a*x(:,:,:)
!$OMP END PARALLEL WORKSHARE
    END IF
  END SUBROUTINE scal_3d_r8

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE axpy_1d_r4(a, x, y, z)
    REAL(4), INTENT(IN)    :: a
    REAL(4), INTENT(IN)    :: x(:)
    REAL(4), INTENT(INOUT) :: y(:)
    REAL(4), INTENT(OUT), OPTIONAL :: z(:)
#ifdef F2008
    CONTIGUOUS x, y, z
#endif
    IF (PRESENT(z)) THEN
!$OMP PARALLEL WORKSHARE
       z(:) = a*x(:) + y(:)
!$OMP END PARALLEL WORKSHARE
    ELSE
!$OMP PARALLEL WORKSHARE
       y(:) = a*x(:) + y(:)
!$OMP END PARALLEL WORKSHARE
    END IF
  END SUBROUTINE axpy_1d_r4

  SUBROUTINE axpy_1d_r8(a, x, y, z)
    REAL(8), INTENT(IN)    :: a
    REAL(8), INTENT(IN)    :: x(:)
    REAL(8), INTENT(INOUT) :: y(:)
    REAL(8), INTENT(OUT), OPTIONAL :: z(:)
#ifdef F2008
    CONTIGUOUS x, y, z
#endif
    IF (PRESENT(z)) THEN
!$OMP PARALLEL WORKSHARE
       z(:) = a*x(:) + y(:)
!$OMP END PARALLEL WORKSHARE
    ELSE
!$OMP PARALLEL WORKSHARE
       y(:) = a*x(:) + y(:)
!$OMP END PARALLEL WORKSHARE
    END IF
  END SUBROUTINE axpy_1d_r8

  SUBROUTINE axpy_2d_r4(a, x, y, z)
    REAL(4), INTENT(IN)    :: a
    REAL(4), INTENT(IN)    :: x(:,:)
    REAL(4), INTENT(INOUT) :: y(:,:)
    REAL(4), INTENT(OUT), OPTIONAL :: z(:,:)
#ifdef F2008
    CONTIGUOUS x, y, z
#endif
    IF (PRESENT(z)) THEN
!$OMP PARALLEL WORKSHARE
       z(:,:) = a*x(:,:) + y(:,:)
!$OMP END PARALLEL WORKSHARE
    ELSE
!$OMP PARALLEL WORKSHARE
       y(:,:) = a*x(:,:) + y(:,:)
!$OMP END PARALLEL WORKSHARE
    END IF
  END SUBROUTINE axpy_2d_r4

  SUBROUTINE axpy_2d_r8(a, x, y, z)
    REAL(8), INTENT(IN)    :: a
    REAL(8), INTENT(IN)    :: x(:,:)
    REAL(8), INTENT(INOUT) :: y(:,:)
    REAL(8), INTENT(OUT), OPTIONAL :: z(:,:)
#ifdef F2008
    CONTIGUOUS x, y, z
#endif
    IF (PRESENT(z)) THEN
!$OMP PARALLEL WORKSHARE
       z(:,:) = a*x(:,:) + y(:,:)
!$OMP END PARALLEL WORKSHARE
    ELSE
!$OMP PARALLEL WORKSHARE
       y(:,:) = a*x(:,:) + y(:,:)
!$OMP END PARALLEL WORKSHARE
    END IF
  END SUBROUTINE axpy_2d_r8

  SUBROUTINE axpy_3d_r4(a, x, y, z)
    REAL(4), INTENT(IN)    :: a
    REAL(4), INTENT(IN)    :: x(:,:,:)
    REAL(4), INTENT(INOUT) :: y(:,:,:)
    REAL(4), INTENT(OUT), OPTIONAL :: z(:,:,:)
#ifdef F2008
    CONTIGUOUS x, y, z
#endif
    IF (PRESENT(z)) THEN
!$OMP PARALLEL WORKSHARE
       z(:,:,:) = a*x(:,:,:) + y(:,:,:)
!$OMP END PARALLEL WORKSHARE
    ELSE
!$OMP PARALLEL WORKSHARE
       y(:,:,:) = a*x(:,:,:) + y(:,:,:)
!$OMP END PARALLEL WORKSHARE
    END IF
  END SUBROUTINE axpy_3d_r4

  SUBROUTINE axpy_3d_r8(a, x, y, z)
    REAL(8), INTENT(IN)    :: a
    REAL(8), INTENT(IN)    :: x(:,:,:)
    REAL(8), INTENT(INOUT) :: y(:,:,:)
    REAL(8), INTENT(OUT), OPTIONAL :: z(:,:,:)
#ifdef F2008
    CONTIGUOUS x, y, z
#endif
    IF (PRESENT(z)) THEN
!$OMP PARALLEL WORKSHARE
       z(:,:,:) = a*x(:,:,:) + y(:,:,:)
!$OMP END PARALLEL WORKSHARE
    ELSE
!$OMP PARALLEL WORKSHARE
       y(:,:,:) = a*x(:,:,:) + y(:,:,:)
!$OMP END PARALLEL WORKSHARE
    END IF
  END SUBROUTINE axpy_3d_r8

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE increase_buffer_i4(x, n)
    INTEGER, ALLOCATABLE, INTENT(INOUT) :: x(:)
    INTEGER, INTENT(IN) :: n
    INTEGER, ALLOCATABLE :: tmp(:)

    INTEGER :: m

    m = 0
    IF (allocated(x)) THEN
       m = size(x)
       ALLOCATE(tmp(m))
       tmp(:) = x(:)
       DEALLOCATE(x)
    END IF
    ALLOCATE(x(m+n))
    IF (allocated(tmp)) THEN
       x(1:m) = tmp(:)
       DEALLOCATE(tmp)
    END IF
  END SUBROUTINE increase_buffer_i4


  SUBROUTINE increase_buffer_i8(x, n)
    INTEGER(8), ALLOCATABLE, INTENT(INOUT) :: x(:)
    INTEGER,    INTENT(IN) :: n
    INTEGER(8), ALLOCATABLE :: tmp(:)

    INTEGER :: m

    m = 0
    IF (allocated(x)) THEN
       m = size(x)
       ALLOCATE(tmp(m))
       tmp(:) = x(:)
       DEALLOCATE(x)
    END IF
    ALLOCATE(x(m+n))
    IF (allocated(tmp)) THEN
       x(1:m) = tmp(:)
       DEALLOCATE(tmp)
    END IF
  END SUBROUTINE increase_buffer_i8

  SUBROUTINE increase_buffer_r4(x, n)
    REAL(4), ALLOCATABLE, INTENT(INOUT) :: x(:)
    INTEGER, INTENT(IN) :: n
    REAL(4), ALLOCATABLE :: tmp(:)

    INTEGER :: m

    m = 0
    IF (allocated(x)) THEN
       m = size(x)
       ALLOCATE(tmp(m))
       tmp(:) = x(:)
       DEALLOCATE(x)
    END IF
    ALLOCATE(x(m+n))
    IF (allocated(tmp)) THEN
       x(1:m) = tmp(:)
       DEALLOCATE(tmp)
    END IF
  END SUBROUTINE increase_buffer_r4

  SUBROUTINE increase_buffer_r8(x, n)
    REAL(8), ALLOCATABLE, INTENT(INOUT) :: x(:)
    INTEGER, INTENT(IN) :: n
    REAL(8), ALLOCATABLE :: tmp(:)

    INTEGER :: m

    m = 0
    IF (allocated(x)) THEN
       m = size(x)
       ALLOCATE(tmp(m))
       tmp(:) = x(:)
       DEALLOCATE(x)
    END IF
    ALLOCATE(x(m+n))
    IF (allocated(tmp)) THEN
       x(1:m) = tmp(:)
       DEALLOCATE(tmp)
    END IF
  END SUBROUTINE increase_buffer_r8

  SUBROUTINE increase_buffer_i4_2d(x, n, l)
    INTEGER(4), ALLOCATABLE, INTENT(INOUT) :: x(:,:)
    INTEGER, INTENT(IN) :: n
    INTEGER, INTENT(IN) :: l
    INTEGER(4), ALLOCATABLE :: tmp(:,:)

    INTEGER :: m

    m = 0
    IF (allocated(x)) THEN
       CALL assert(l==size(x,1), "dimension mismatch in INCREASE_BUFFER_R4_2D")
       m = size(x,2)
       ALLOCATE(tmp(l,m))
       tmp(:,:) = x(:,:)
       DEALLOCATE(x)
    END IF
    ALLOCATE(x(l,m+n))
    IF (allocated(tmp)) THEN
       x(:,1:m) = tmp(:,:)
       DEALLOCATE(tmp)
    END IF
  END SUBROUTINE increase_buffer_i4_2d

  SUBROUTINE increase_buffer_i8_2d(x, n, l)
    INTEGER(8), ALLOCATABLE, INTENT(INOUT) :: x(:,:)
    INTEGER, INTENT(IN) :: n
    INTEGER, INTENT(IN) :: l
    INTEGER(8), ALLOCATABLE :: tmp(:,:)

    INTEGER :: m

    m = 0
    IF (allocated(x)) THEN
       CALL assert(l==size(x,1), "dimension mismatch in INCREASE_BUFFER_R4_2D")
       m = size(x,2)
       ALLOCATE(tmp(l,m))
       tmp(:,:) = x(:,:)
       DEALLOCATE(x)
    END IF
    ALLOCATE(x(l,m+n))
    IF (allocated(tmp)) THEN
       x(:,1:m) = tmp(:,:)
       DEALLOCATE(tmp)
    END IF
  END SUBROUTINE increase_buffer_i8_2d

  SUBROUTINE increase_buffer_r4_2d(x, n, l)
    REAL(4), ALLOCATABLE, INTENT(INOUT) :: x(:,:)
    INTEGER, INTENT(IN) :: n
    INTEGER, INTENT(IN) :: l
    REAL(4), ALLOCATABLE :: tmp(:,:)

    INTEGER :: m

    m = 0
    IF (allocated(x)) THEN
       CALL assert(l==size(x,1), "dimension mismatch in INCREASE_BUFFER_R4_2D")
       m = size(x,2)
       ALLOCATE(tmp(l,m))
       tmp(:,:) = x(:,:)
       DEALLOCATE(x)
    END IF
    ALLOCATE(x(l,m+n))
    IF (allocated(tmp)) THEN
       x(:,1:m) = tmp(:,:)
       DEALLOCATE(tmp)
    END IF
  END SUBROUTINE increase_buffer_r4_2d

  SUBROUTINE increase_buffer_r8_2d(x, n, l)
    REAL(8), ALLOCATABLE, INTENT(INOUT) :: x(:,:)
    INTEGER, INTENT(IN) :: n
    INTEGER, INTENT(IN) :: l
    REAL(8), ALLOCATABLE :: tmp(:,:)

    INTEGER :: m

    m = 0
    IF (allocated(x)) THEN
       CALL assert(l==size(x,1), "dimension mismatch in INCREASE_BUFFER_R8_2D")
       m = size(x,2)
       ALLOCATE(tmp(l,m))
       tmp(:,:) = x(:,:)
       DEALLOCATE(x)
    END IF
    ALLOCATE(x(l,m+n))
    IF (allocated(tmp)) THEN
       x(:,1:m) = tmp(:,:)
       DEALLOCATE(tmp)
    END IF
  END SUBROUTINE increase_buffer_r8_2d

!-----------------------------------------------------------------------------------------------------------------------

  LOGICAL PURE FUNCTION is_number(str)
    CHARACTER(*), INTENT(IN) :: str

    INTEGER :: i
    LOGICAL :: s

    is_number = .FALSE.
    DO i=1, len(str)
       IF (iachar(str(i:i)) < iachar('0') .OR. iachar(str(i:i)) > iachar('9')) RETURN
    END DO
    is_number = .TRUE.

  END FUNCTION is_number

!-----------------------------------------------------------------------------------------------------------------------

  INTEGER(4) PURE FUNCTION xorshift32(x)
    INTEGER(4), INTENT(IN) :: x

    xorshift32 = x
    IF (xorshift32 > 0) THEN
       xorshift32 = IEOR(ISHFT(xorshift32,  13), xorshift32)
       xorshift32 = IEOR(ISHFT(xorshift32, -17), xorshift32)
       xorshift32 = IEOR(ISHFT(xorshift32,   5), xorshift32)
    ELSE
       xorshift32 = IEOR(ISHFT(xorshift32,  -3), xorshift32)
       xorshift32 = IEOR(ISHFT(xorshift32,   7), xorshift32)
       xorshift32 = IEOR(ISHFT(xorshift32, -29), xorshift32)
    END IF
  END FUNCTION xorshift32

  INTEGER(8) PURE FUNCTION xorshift64(x)
    INTEGER(8), INTENT(IN) :: x

    xorshift64 = x
    IF (xorshift64 > 0) THEN
       xorshift64 = IEOR(ISHFT(xorshift64,  13), xorshift64)
       xorshift64 = IEOR(ISHFT(xorshift64,  -7), xorshift64)
       xorshift64 = IEOR(ISHFT(xorshift64,  17), xorshift64)
    ELSE
       xorshift64 = IEOR(ISHFT(xorshift64,  -3), xorshift64)
       xorshift64 = IEOR(ISHFT(xorshift64,  29), xorshift64)
       xorshift64 = IEOR(ISHFT(xorshift64, -47), xorshift64)
    END IF
  END FUNCTION xorshift64

!-----------------------------------------------------------------------------------------------------------------------

  INTEGER(4) PURE FUNCTION mkseed(i4, i8, r4, r8, str)
    INTEGER(4),   INTENT(IN), OPTIONAL :: i4
    INTEGER(8),   INTENT(IN), OPTIONAL :: i8
    REAL(4),      INTENT(IN), OPTIONAL :: r4
    REAL(8),      INTENT(IN), OPTIONAL :: r8
    CHARACTER(*), INTENT(IN), OPTIONAL :: str
    INTEGER :: i

    mkseed = -1831433054 !magic number from Marsagalia (2003) but transfeed to signed int32)

    IF (PRESENT(r4)) mkseed = IEOR(xorshift32(mkseed), TRANSFER(r4, 0_4))
    IF (PRESENT(r8)) mkseed = IEOR(xorshift32(mkseed), IEOR(TRANSFER(r8, 0_4), INT(ISHFT(TRANSFER(r8, 0_8), -32), 4)))
    IF (PRESENT(i4)) mkseed = IEOR(xorshift32(mkseed), i4)
    IF (PRESENT(i8)) mkseed = IEOR(xorshift32(mkseed), IEOR(TRANSFER(i8, 0_4), INT(ISHFT(i8, -32), 4)))

    IF (PRESENT(str)) THEN
       DO i=1, len(str), 4
          mkseed = xorshift32(mkseed)
          mkseed = IEOR(mkseed, TRANSFER(str(i:min(i+3,len(str))), 0_4))
       END DO
    END IF

    mkseed = xorshift32(xorshift32(mkseed))
  END FUNCTION mkseed

!-----------------------------------------------------------------------------------------------------------------------

END MODULE misc
