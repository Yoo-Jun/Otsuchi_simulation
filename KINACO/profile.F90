#include "macro.h"

MODULE profile
#ifdef PARALLEL_MPI
  USE mpi
#endif
  IMPLICIT NONE
  PRIVATE
  PUBLIC profile_start, profile_stop, profile_time

  INTEGER, PARAMETER :: max_name  = 256
  INTEGER, PARAMETER :: max_level = 32

  CHARACTER(32), SAVE :: names(max_name) = ''
  REAL(8),       SAVE :: time(0:max_name,0:max_level) = 0.0D0

  INTEGER, SAVE :: count = 0

CONTAINS

  SUBROUTINE profile_start(name, level)
    CHARACTER(*), INTENT(IN) :: name
    INTEGER,      INTENT(IN) :: level

    INTEGER :: n
    REAL    :: cputime

    n = lookup(name)

    IF (n == 0) THEN
       count = count+1
       n = count
       names(n) = name
    END IF

#ifdef PARALLEL_MPI
    time(n,level) = time(n,level) - mpi_wtime()
#else
    CALL cpu_time(cputime)
    time(n,level) = time(n,level) - cputime
#endif
  END SUBROUTINE profile_start

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE profile_stop(name, level)
    CHARACTER(*), INTENT(IN) :: name
    INTEGER,      INTENT(IN) :: level

    INTEGER :: n
    REAL    :: cputime

    n = lookup(name)

#ifdef PARALLEL_MPI
    time(n,level) = time(n,level) + mpi_wtime()
#else
    CALL cpu_time(cputime)
    time(n,level) = time(n,level) + cputime
#endif
  END SUBROUTINE profile_stop

!-----------------------------------------------------------------------------------------------------------------------

  REAL(8) PURE FUNCTION profile_time(name, level)
    CHARACTER(*), INTENT(IN) :: name
    INTEGER,      INTENT(IN), OPTIONAL :: level

    INTEGER :: n

    n = lookup(name)

    IF (present(level)) THEN
       profile_time = time(n,level)
    ELSE
       profile_time = sum(time(n,0:max_level))
    END IF

  END FUNCTION profile_time

!-----------------------------------------------------------------------------------------------------------------------

  INTEGER PURE FUNCTION lookup(name)
    CHARACTER(*), INTENT(IN) :: name

    INTEGER :: n
    lookup = 0
    DO n=1, max_name
       IF (trim(name) == trim(names(n))) THEN
          lookup = n
          EXIT
       END IF
    END DO
  END FUNCTION lookup

END MODULE profile
