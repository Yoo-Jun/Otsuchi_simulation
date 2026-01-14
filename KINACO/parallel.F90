#include "macro.h"

MODULE parallel
  USE misc
#ifdef PARALLEL_MPI
#ifndef NO_MPI_MOD
  USE mpi
#endif
#endif
#ifndef NO_OMPLIB_MOD
!$ USE omp_lib
#endif
  IMPLICIT NONE
#ifdef NO_MPI_MOD
!$ include "mpif.h"
#endif
#ifdef NO_OMPLIB_MOD
!$ include "omp_lib.h"
#endif

  INTEGER, SAVE :: npes = 1
  INTEGER, SAVE :: rank = 0

  INTEGER, SAVE :: hpes  = 1
  INTEGER, SAVE :: hrank = 0

  INTEGER, SAVE :: vpes  = 1
  INTEGER, SAVE :: vrank = 0

  INTEGER, SAVE :: ipes = 1
  INTEGER, SAVE :: jpes = 1
  INTEGER, SAVE :: kpes = 1

#ifdef PARALLEL_MPI
  INTEGER, SAVE :: comm   = MPI_COMM_WORLD
  INTEGER, SAVE :: hcomm  = MPI_COMM_WORLD
  INTEGER, SAVE :: vcomm  = MPI_COMM_WORLD
  INTEGER, SAVE :: scomm  = MPI_COMM_NULL
  INTEGER, SAVE :: gcomm  = MPI_COMM_NULL

  INTEGER, SAVE :: node_dim(3) = (/1, 1, 1/)

  INTEGER, SAVE :: spes  = 1
  INTEGER, SAVE :: srank = 0

  INTEGER, SAVE :: gpes  = 1
  INTEGER, SAVE :: grank = 0

  INTEGER, SAVE :: rank_e
  INTEGER, SAVE :: rank_w
  INTEGER, SAVE :: rank_n
  INTEGER, SAVE :: rank_s
  INTEGER, SAVE :: rank_u
  INTEGER, SAVE :: rank_l

  INTEGER, SAVE :: rank_ne
  INTEGER, SAVE :: rank_nw
  INTEGER, SAVE :: rank_se
  INTEGER, SAVE :: rank_sw
  INTEGER, SAVE :: rank_ue
  INTEGER, SAVE :: rank_uw
  INTEGER, SAVE :: rank_un
  INTEGER, SAVE :: rank_us
  INTEGER, SAVE :: rank_le
  INTEGER, SAVE :: rank_lw
  INTEGER, SAVE :: rank_ln
  INTEGER, SAVE :: rank_ls

  INTEGER, SAVE :: rank_une
  INTEGER, SAVE :: rank_unw
  INTEGER, SAVE :: rank_use
  INTEGER, SAVE :: rank_usw
  INTEGER, SAVE :: rank_lne
  INTEGER, SAVE :: rank_lnw
  INTEGER, SAVE :: rank_lse
  INTEGER, SAVE :: rank_lsw

  INTEGER, SAVE :: rank_tp(-1:1, -1:1)

  INTEGER, PARAMETER :: tag_e   =  1
  INTEGER, PARAMETER :: tag_w   =  2
  INTEGER, PARAMETER :: tag_n   =  3
  INTEGER, PARAMETER :: tag_s   =  4
  INTEGER, PARAMETER :: tag_ne  =  5
  INTEGER, PARAMETER :: tag_nw  =  6
  INTEGER, PARAMETER :: tag_se  =  7
  INTEGER, PARAMETER :: tag_sw  =  8
#ifdef PARALLEL3D
  INTEGER, PARAMETER :: tag_u   =  9
  INTEGER, PARAMETER :: tag_l   = 10
  INTEGER, PARAMETER :: tag_ue  = 11
  INTEGER, PARAMETER :: tag_uw  = 12
  INTEGER, PARAMETER :: tag_un  = 13
  INTEGER, PARAMETER :: tag_us  = 14
  INTEGER, PARAMETER :: tag_le  = 15
  INTEGER, PARAMETER :: tag_lw  = 16
  INTEGER, PARAMETER :: tag_ln  = 17
  INTEGER, PARAMETER :: tag_ls  = 18
  INTEGER, PARAMETER :: tag_une = 19
  INTEGER, PARAMETER :: tag_unw = 20
  INTEGER, PARAMETER :: tag_use = 21
  INTEGER, PARAMETER :: tag_usw = 22
  INTEGER, PARAMETER :: tag_lne = 23
  INTEGER, PARAMETER :: tag_lnw = 24
  INTEGER, PARAMETER :: tag_lse = 25
  INTEGER, PARAMETER :: tag_lsw = 26

  INTEGER, PARAMETER :: n_sendrecv = 26
#else
  INTEGER, PARAMETER :: n_sendrecv = 8
#endif

  INTEGER :: sendrank(n_sendrecv)
  INTEGER :: recvrank(n_sendrecv)

  INTEGER, ALLOCATABLE :: icoords(:)
  INTEGER, ALLOCATABLE :: jcoords(:)
  INTEGER, ALLOCATABLE :: kcoords(:)

  INTEGER, ALLOCATABLE :: icoords_h(:)
  INTEGER, ALLOCATABLE :: jcoords_h(:)
  INTEGER, ALLOCATABLE :: kcoords_v(:)

  INTEGER, ALLOCATABLE :: ranks(:,:,:)
  INTEGER, ALLOCATABLE :: hranks(:,:)
  INTEGER, ALLOCATABLE :: vranks(:)

  LOGICAL, SAVE :: remove_masked_pe = .FALSE.
#endif
  INTEGER, SAVE :: icoord = 0
  INTEGER, SAVE :: jcoord = 0
  INTEGER, SAVE :: kcoord = 0

  INTEGER, PARAMETER :: max_threads = 16
  INTEGER, SAVE      :: nthreads = 1

  CHARACTER(8), SAVE :: strrank = '0'

  INTERFACE bcast
     MODULE PROCEDURE bcast_real4
     MODULE PROCEDURE bcast_real4_array
     MODULE PROCEDURE bcast_real8
     MODULE PROCEDURE bcast_real8_array
     MODULE PROCEDURE bcast_integer4
     MODULE PROCEDURE bcast_integer4_array
     MODULE PROCEDURE bcast_integer8
     MODULE PROCEDURE bcast_integer8_array
     MODULE PROCEDURE bcast_byte
     MODULE PROCEDURE bcast_byte_array
     MODULE PROCEDURE bcast_byte_array2d
     MODULE PROCEDURE bcast_logical
     MODULE PROCEDURE bcast_character
  END INTERFACE

  INTERFACE vsum
     MODULE PROCEDURE vsum_r4_3d
     MODULE PROCEDURE vsum_r8_3d
     MODULE PROCEDURE vsum_r4_2d
     MODULE PROCEDURE vsum_r8_2d
     MODULE PROCEDURE vsum_r4_1d
     MODULE PROCEDURE vsum_r8_1d
     MODULE PROCEDURE vsum_r4
     MODULE PROCEDURE vsum_r8
  END INTERFACE

  INTERFACE vcast
     MODULE PROCEDURE vcast_r4_3d
     MODULE PROCEDURE vcast_r8_3d
     MODULE PROCEDURE vcast_r4_2d
     MODULE PROCEDURE vcast_r8_2d
  END INTERFACE

  INTERFACE vmax
     MODULE PROCEDURE vmax_r8_2d
     MODULE PROCEDURE vmax_r4_2d
     MODULE PROCEDURE vmax_i8_2d
     MODULE PROCEDURE vmax_i4_2d
     MODULE PROCEDURE vmax_r8
     MODULE PROCEDURE vmax_r4
  END INTERFACE

  INTERFACE vmin
     MODULE PROCEDURE vmin_r8_2d
     MODULE PROCEDURE vmin_r4_2d
     MODULE PROCEDURE vmin_i8_2d
     MODULE PROCEDURE vmin_i4_2d
     MODULE PROCEDURE vmin_r8
     MODULE PROCEDURE vmin_r4
  END INTERFACE

  INTERFACE hsum
     MODULE PROCEDURE hsum_r4_3d
     MODULE PROCEDURE hsum_r8_3d
     MODULE PROCEDURE hsum_r4_2d
     MODULE PROCEDURE hsum_r8_2d
     MODULE PROCEDURE hsum_r4_1d
     MODULE PROCEDURE hsum_r8_1d
  END INTERFACE

  INTERFACE hcast
     MODULE PROCEDURE hcast_r4_3d
     MODULE PROCEDURE hcast_r8_3d
     MODULE PROCEDURE hcast_r4_2d
     MODULE PROCEDURE hcast_r8_2d
  END INTERFACE

  INTERFACE hmax
     MODULE PROCEDURE hmax_r8
     MODULE PROCEDURE hmax_r4
  END INTERFACE

  INTERFACE hmin
     MODULE PROCEDURE hmin_r8
     MODULE PROCEDURE hmin_r4
  END INTERFACE


  INTERFACE gsum
     MODULE PROCEDURE gsum_r4
     MODULE PROCEDURE gsum_r8
     MODULE PROCEDURE gsum_i4
     MODULE PROCEDURE gsum_i8
     MODULE PROCEDURE gsum_r4_1d
     MODULE PROCEDURE gsum_r8_1d
     MODULE PROCEDURE gsum_i4_1d
     MODULE PROCEDURE gsum_i8_1d
     MODULE PROCEDURE gsum_r4_2d
     MODULE PROCEDURE gsum_r8_2d
     MODULE PROCEDURE gsum_i4_2d
     MODULE PROCEDURE gsum_i8_2d
  END INTERFACE

  INTERFACE gmax
     MODULE PROCEDURE gmax_r4
     MODULE PROCEDURE gmax_r8
     MODULE PROCEDURE gmax_i4
     MODULE PROCEDURE gmax_i8
     MODULE PROCEDURE gmax_r4_1d
     MODULE PROCEDURE gmax_r8_1d
     MODULE PROCEDURE gmax_i4_1d
     MODULE PROCEDURE gmax_i8_1d
     MODULE PROCEDURE gmax_r4_2d
     MODULE PROCEDURE gmax_r8_2d
     MODULE PROCEDURE gmax_i4_2d
     MODULE PROCEDURE gmax_i8_2d
  END INTERFACE

  INTERFACE gmin
     MODULE PROCEDURE gmin_r4
     MODULE PROCEDURE gmin_r8
     MODULE PROCEDURE gmin_i4
     MODULE PROCEDURE gmin_i8
     MODULE PROCEDURE gmin_r4_1d
     MODULE PROCEDURE gmin_r8_1d
     MODULE PROCEDURE gmin_i4_1d
     MODULE PROCEDURE gmin_i8_1d
     MODULE PROCEDURE gmin_r4_2d
     MODULE PROCEDURE gmin_r8_2d
     MODULE PROCEDURE gmin_i4_2d
     MODULE PROCEDURE gmin_i8_2d
  END INTERFACE

  INTERFACE usend
     MODULE PROCEDURE usend_r4_3d
     MODULE PROCEDURE usend_r8_3d
     MODULE PROCEDURE usend_r4_2d
     MODULE PROCEDURE usend_r8_2d
     MODULE PROCEDURE usend_r4_1d
     MODULE PROCEDURE usend_r8_1d
     MODULE PROCEDURE usend_r4
     MODULE PROCEDURE usend_r8
  END INTERFACE usend

  INTERFACE urecv
     MODULE PROCEDURE urecv_r4_3d
     MODULE PROCEDURE urecv_r8_3d
     MODULE PROCEDURE urecv_r4_2d
     MODULE PROCEDURE urecv_r8_2d
     MODULE PROCEDURE urecv_r4_1d
     MODULE PROCEDURE urecv_r8_1d
     MODULE PROCEDURE urecv_r4
     MODULE PROCEDURE urecv_r8
  END INTERFACE urecv

  INTERFACE lsend
     MODULE PROCEDURE lsend_r4_3d
     MODULE PROCEDURE lsend_r8_3d
     MODULE PROCEDURE lsend_r4_2d
     MODULE PROCEDURE lsend_r8_2d
     MODULE PROCEDURE lsend_r4_1d
     MODULE PROCEDURE lsend_r8_1d
     MODULE PROCEDURE lsend_r4
     MODULE PROCEDURE lsend_r8
  END INTERFACE lsend

  INTERFACE lrecv
     MODULE PROCEDURE lrecv_r4_3d
     MODULE PROCEDURE lrecv_r8_3d
     MODULE PROCEDURE lrecv_r4_2d
     MODULE PROCEDURE lrecv_r8_2d
     MODULE PROCEDURE lrecv_r4_1d
     MODULE PROCEDURE lrecv_r8_1d
     MODULE PROCEDURE lrecv_r4
     MODULE PROCEDURE lrecv_r8
  END INTERFACE lrecv

CONTAINS
  SUBROUTINE init_parallel
    INTEGER :: ierr
    CHARACTER(16) :: env

#ifdef PARALLEL_MPI
    CALL mpi_init(ierr)
!   CALL mpi_pcontrol(0, ierr)
    comm = MPI_COMM_WORLD
    CALL mpi_comm_rank(comm, rank, ierr)

    strrank = adjustl(format(rank))
#endif

#ifdef F2003
!$  CALL get_environment_variable('OMP_NUM_THREADS', env, status=ierr)
!$  IF (ierr/=0) THEN
!$     IF (rank==0) WRITE(0, *) "OMP_NUM_THREADS is not set!"
!$     CALL omp_set_num_threads(1)
!$  END IF
#endif
!$  nthreads = omp_get_max_threads()

  END SUBROUTINE init_parallel

  SUBROUTINE finalize_parallel
    INTEGER :: ierr

#ifdef PARALLEL_MPI
  CALL mpi_barrier(comm, ierr)
  CALL mpi_finalize(ierr)
#endif
  END SUBROUTINE finalize_parallel

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE bcast_real4(x)
    REAL(4), INTENT(INOUT) :: x
    INTEGER :: ierr
#ifdef PARALLEL_MPI
    CALL mpi_bcast(x, 1, MPI_REAL4, 0, comm, ierr)
#endif
  END SUBROUTINE bcast_real4

  SUBROUTINE bcast_real4_array(x)
    REAL(4), INTENT(INOUT) :: x(:)
    INTEGER :: ierr
#ifdef F2008
    CONTIGUOUS x
#endif
#ifdef PARALLEL_MPI
    CALL mpi_bcast(x, size(x), MPI_REAL4, 0, comm, ierr)
#endif
  END SUBROUTINE bcast_real4_array

  SUBROUTINE bcast_real8(x)
    REAL(8), INTENT(INOUT) :: x
    INTEGER :: ierr
#ifdef PARALLEL_MPI
    CALL mpi_bcast(x, 1, MPI_REAL8, 0, comm, ierr)
#endif
  END SUBROUTINE bcast_real8

  SUBROUTINE bcast_real8_array(x)
    REAL(8), INTENT(INOUT) :: x(:)
    INTEGER :: ierr
#ifdef F2008
    CONTIGUOUS x
#endif
#ifdef PARALLEL_MPI
    CALL mpi_bcast(x, size(x), MPI_REAL8, 0, comm, ierr)
#endif
  END SUBROUTINE bcast_real8_array

  SUBROUTINE bcast_integer4(x)
    INTEGER(4), INTENT(INOUT) :: x
    INTEGER :: ierr
#ifdef PARALLEL_MPI
    CALL mpi_bcast(x, 1, MPI_INTEGER4, 0, comm, ierr)
#endif
  END SUBROUTINE bcast_integer4

  SUBROUTINE bcast_integer4_array(x)
    INTEGER(4), INTENT(INOUT) :: x(:)
    INTEGER :: ierr
#ifdef F2008
    CONTIGUOUS x
#endif
#ifdef PARALLEL_MPI
    CALL mpi_bcast(x, size(x), MPI_INTEGER4, 0, comm, ierr)
#endif
  END SUBROUTINE bcast_integer4_array

  SUBROUTINE bcast_integer8(x)
    INTEGER(8), INTENT(INOUT) :: x
    INTEGER :: ierr
#ifdef PARALLEL_MPI
    CALL mpi_bcast(x, 1, MPI_INTEGER8, 0, comm, ierr)
#endif
  END SUBROUTINE bcast_integer8

  SUBROUTINE bcast_integer8_array(x)
    INTEGER(8), INTENT(INOUT) :: x(:)
    INTEGER :: ierr
#ifdef F2008
    CONTIGUOUS x
#endif
#ifdef PARALLEL_MPI
    CALL mpi_bcast(x, size(x), MPI_INTEGER8, 0, comm, ierr)
#endif
  END SUBROUTINE bcast_integer8_array

  SUBROUTINE bcast_byte(x)
    INTEGER(1), INTENT(INOUT) :: x
    INTEGER :: ierr
#ifdef PARALLEL_MPI
    CALL mpi_bcast(x, 1, MPI_BYTE, 0, comm, ierr)
#endif
  END SUBROUTINE bcast_byte

  SUBROUTINE bcast_byte_array(x)
    INTEGER(1), INTENT(INOUT) :: x(:)
    INTEGER :: ierr
#ifdef F2008
    CONTIGUOUS x
#endif
#ifdef PARALLEL_MPI
    CALL mpi_bcast(x, size(x), MPI_BYTE, 0, comm, ierr)
#endif
  END SUBROUTINE bcast_byte_array

  SUBROUTINE bcast_byte_array2d(x)
    INTEGER(1), INTENT(INOUT) :: x(:,:)
    INTEGER :: ierr
#ifdef F2008
    CONTIGUOUS x
#endif
#ifdef PARALLEL_MPI
    CALL mpi_bcast(x, size(x), MPI_BYTE, 0, comm, ierr)
#endif
  END SUBROUTINE bcast_byte_array2d

  SUBROUTINE bcast_logical(x)
    LOGICAL, INTENT(INOUT) :: x
    INTEGER :: ierr
#ifdef PARALLEL_MPI
    CALL mpi_bcast(x, 1, MPI_LOGICAL, 0, comm, ierr)
#endif
  END SUBROUTINE bcast_logical

  SUBROUTINE bcast_character(str)
    CHARACTER(*), INTENT(INOUT) :: str
    INTEGER :: ierr
#ifdef PARALLEL_MPI
    CALL mpi_bcast(str, len(str), MPI_CHARACTER, 0, comm, ierr)
#endif
  END SUBROUTINE bcast_character

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE vsum_r8_3d(x, all)
    REAL(8), INTENT(INOUT) :: x(:,:,:)
    LOGICAL, INTENT(IN), OPTIONAL :: all
#ifdef F2008
    CONTIGUOUS x
#endif
    LOGICAL :: all_
    INTEGER :: dummy, ierr

#ifdef PARALLEL3D
    IF (vpes == 1) RETURN

    all_ = .FALSE.
    IF (present(all)) all_ = all

    IF (all_) THEN
       CALL mpi_allreduce(MPI_IN_PLACE, x, size(x), MPI_REAL8, MPI_SUM, vcomm, ierr)
    ELSE
       IF (vrank==0) THEN
          CALL mpi_reduce(MPI_IN_PLACE, x, size(x), MPI_REAL8, MPI_SUM, 0, vcomm, ierr)
       ELSE
          CALL mpi_reduce(x,        dummy, size(x), MPI_REAL8, MPI_SUM, 0, vcomm, ierr)
       END IF
    END IF
#ifdef REDUCE_CONSISTENT
    x(:,:,:) = REAL(x, KIND=4)
#endif
#endif

  END SUBROUTINE vsum_r8_3d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE vsum_r4_3d(x, all)
    REAL(4), INTENT(INOUT) :: x(:,:,:)
    LOGICAL, INTENT(IN), OPTIONAL :: all
#ifdef F2008
    CONTIGUOUS x
#endif
    LOGICAL :: all_
    INTEGER :: dummy, ierr

#ifdef PARALLEL3D
    IF (vpes == 1) RETURN

    all_ = .FALSE.
    IF (present(all)) all_ = all

    IF (all_) THEN
       CALL mpi_allreduce(MPI_IN_PLACE, x, size(x), MPI_REAL4, MPI_SUM, vcomm, ierr)
    ELSE
       IF (vrank==0) THEN
          CALL mpi_reduce(MPI_IN_PLACE, x, size(x), MPI_REAL4, MPI_SUM, 0, vcomm, ierr)
       ELSE
          CALL mpi_reduce(x,        dummy, size(x), MPI_REAL4, MPI_SUM, 0, vcomm, ierr)
       END IF
    END IF
#endif

  END SUBROUTINE vsum_r4_3d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE vsum_r8_2d(x, all)
    REAL(8), INTENT(INOUT) :: x(:,:)
    LOGICAL, INTENT(IN), OPTIONAL :: all
#ifdef F2008
    CONTIGUOUS x
#endif
    LOGICAL :: all_
    INTEGER :: dummy, ierr

#ifdef PARALLEL3D
    IF (vpes == 1) RETURN

    all_ = .FALSE.
    IF (present(all)) all_ = all

    IF (all_) THEN
       CALL mpi_allreduce(MPI_IN_PLACE, x, size(x), MPI_REAL8, MPI_SUM, vcomm, ierr)
    ELSE
       IF (vrank==0) THEN
          CALL mpi_reduce(MPI_IN_PLACE, x, size(x), MPI_REAL8, MPI_SUM, 0, vcomm, ierr)
       ELSE
          CALL mpi_reduce(x,        dummy, size(x), MPI_REAL8, MPI_SUM, 0, vcomm, ierr)
       END IF
    END IF
#ifdef REDUCE_CONSISTENT
    x(:,:) = REAL(x, KIND=4)
#endif
#endif

  END SUBROUTINE vsum_r8_2d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE vsum_r4_2d(x, all)
    REAL(4), INTENT(INOUT) :: x(:,:)
    LOGICAL, INTENT(IN), OPTIONAL :: all
#ifdef F2008
    CONTIGUOUS x
#endif
    LOGICAL :: all_
    INTEGER :: dummy, ierr

#ifdef PARALLEL3D
    IF (vpes == 1) RETURN

    all_ = .FALSE.
    IF (present(all)) all_ = all

    IF (all_) THEN
       CALL mpi_allreduce(MPI_IN_PLACE, x, size(x), MPI_REAL4, MPI_SUM, vcomm, ierr)
    ELSE
       IF (vrank==0) THEN
          CALL mpi_reduce(MPI_IN_PLACE, x, size(x), MPI_REAL4, MPI_SUM, 0, vcomm, ierr)
       ELSE
          CALL mpi_reduce(x,        dummy, size(x), MPI_REAL4, MPI_SUM, 0, vcomm, ierr)
       END IF
    END IF
#endif

  END SUBROUTINE vsum_r4_2d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE vsum_r8_1d(x, all)
    REAL(8), INTENT(INOUT) :: x(:)
    LOGICAL, INTENT(IN), OPTIONAL :: all
#ifdef F2008
    CONTIGUOUS x
#endif
    LOGICAL :: all_
    INTEGER :: dummy, ierr

#ifdef PARALLEL3D
    IF (vpes == 1) RETURN

    all_ = .FALSE.
    IF (present(all)) all_ = all

    IF (all_) THEN
       CALL mpi_allreduce(MPI_IN_PLACE, x, size(x), MPI_REAL8, MPI_SUM, vcomm, ierr)
    ELSE
       IF (vrank==0) THEN
          CALL mpi_reduce(MPI_IN_PLACE, x, size(x), MPI_REAL8, MPI_SUM, 0, vcomm, ierr)
       ELSE
          CALL mpi_reduce(x,        dummy, size(x), MPI_REAL8, MPI_SUM, 0, vcomm, ierr)
       END IF
    END IF
#ifdef REDUCE_CONSISTENT
    x(:) = REAL(x, KIND=4)
#endif
#endif

  END SUBROUTINE vsum_r8_1d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE vsum_r4_1d(x, all)
    REAL(4), INTENT(INOUT) :: x(:)
    LOGICAL, INTENT(IN), OPTIONAL :: all
#ifdef F2008
    CONTIGUOUS x
#endif
    LOGICAL :: all_
    INTEGER :: dummy, ierr

#ifdef PARALLEL3D
    IF (vpes == 1) RETURN

    all_ = .FALSE.
    IF (present(all)) all_ = all

    IF (all_) THEN
       CALL mpi_allreduce(MPI_IN_PLACE, x, size(x), MPI_REAL4, MPI_SUM, vcomm, ierr)
    ELSE
       IF (vrank==0) THEN
          CALL mpi_reduce(MPI_IN_PLACE, x, size(x), MPI_REAL4, MPI_SUM, 0, vcomm, ierr)
       ELSE
          CALL mpi_reduce(x,        dummy, size(x), MPI_REAL4, MPI_SUM, 0, vcomm, ierr)
       END IF
    END IF
#endif

  END SUBROUTINE vsum_r4_1d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE vsum_r8(x, all)
    REAL(8), INTENT(INOUT) :: x
    LOGICAL, INTENT(IN), OPTIONAL :: all
    LOGICAL :: all_
    INTEGER :: dummy, ierr

#ifdef PARALLEL3D
    IF (vpes == 1) RETURN

    all_ = .FALSE.
    IF (present(all)) all_ = all

    IF (all_) THEN
       CALL mpi_allreduce(MPI_IN_PLACE, x, 1, MPI_REAL8, MPI_SUM, vcomm, ierr)
    ELSE
       IF (vrank==0) THEN
          CALL mpi_reduce(MPI_IN_PLACE, x, 1, MPI_REAL8, MPI_SUM, 0, vcomm, ierr)
       ELSE
          CALL mpi_reduce(x,        dummy, 1, MPI_REAL8, MPI_SUM, 0, vcomm, ierr)
       END IF
    END IF
#ifdef REDUCE_CONSISTENT
    x = REAL(x, KIND=4)
#endif
#endif

  END SUBROUTINE vsum_r8

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE vsum_r4(x, all)
    REAL(4), INTENT(INOUT) :: x
    LOGICAL, INTENT(IN), OPTIONAL :: all
    LOGICAL :: all_
    INTEGER :: dummy, ierr

#ifdef PARALLEL3D
    IF (vpes == 1) RETURN

    all_ = .FALSE.
    IF (present(all)) all_ = all

    IF (all_) THEN
       CALL mpi_allreduce(MPI_IN_PLACE, x, 1, MPI_REAL4, MPI_SUM, vcomm, ierr)
    ELSE
       IF (vrank==0) THEN
          CALL mpi_reduce(MPI_IN_PLACE, x, 1, MPI_REAL4, MPI_SUM, 0, vcomm, ierr)
       ELSE
          CALL mpi_reduce(x,        dummy, 1, MPI_REAL4, MPI_SUM, 0, vcomm, ierr)
       END IF
    END IF
#endif

  END SUBROUTINE vsum_r4

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE vcast_r8_2d(x)
    REAL(8), INTENT(INOUT) :: x(:,:)
#ifdef F2008
    CONTIGUOUS x
#endif
    INTEGER :: ierr

#ifdef PARALLEL3D
    IF (vpes == 1) RETURN

    CALL mpi_bcast(x, size(x), MPI_REAL8, 0, vcomm, ierr)
#endif

  END SUBROUTINE vcast_r8_2d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE vcast_r4_2d(x)
    REAL(4), INTENT(INOUT) :: x(:,:)
#ifdef F2008
    CONTIGUOUS x
#endif
    INTEGER :: ierr

#ifdef PARALLEL3D
    IF (vpes == 1) RETURN

    CALL mpi_bcast(x, size(x), MPI_REAL4, 0, vcomm, ierr)
#endif

  END SUBROUTINE vcast_r4_2d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE vcast_r8_3d(x)
    REAL(8), INTENT(INOUT) :: x(:,:,:)
#ifdef F2008
    CONTIGUOUS x
#endif
    INTEGER :: ierr

#ifdef PARALLEL3D
    IF (vpes == 1) RETURN

    CALL mpi_bcast(x, size(x), MPI_REAL8, 0, vcomm, ierr)
#endif

  END SUBROUTINE vcast_r8_3d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE vcast_r4_3d(x)
    REAL(4), INTENT(INOUT) :: x(:,:,:)
#ifdef F2008
    CONTIGUOUS x
#endif
    INTEGER :: ierr

#ifdef PARALLEL3D
    IF (vpes == 1) RETURN

    CALL mpi_bcast(x, size(x), MPI_REAL4, 0, vcomm, ierr)
#endif

  END SUBROUTINE vcast_r4_3d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE vmax_r4_2d(x, all)
    REAL(4), INTENT(INOUT) :: x(:,:)
    LOGICAL, INTENT(IN), OPTIONAL :: all
#ifdef F2008
    CONTIGUOUS x
#endif
    LOGICAL :: all_
    INTEGER :: dummy, ierr

#ifdef PARALLEL3D
    IF (vpes == 1) RETURN

    all_ = .FALSE.
    IF (present(all)) all_ = all

    IF (all_) THEN
       CALL mpi_allreduce(MPI_IN_PLACE, x, size(x), MPI_REAL4, MPI_MAX, vcomm, ierr)
    ELSE
       IF (vrank==0) THEN
          CALL mpi_reduce(MPI_IN_PLACE, x, size(x), MPI_REAL4, MPI_MAX, 0, vcomm, ierr)
       ELSE
          CALL mpi_reduce(x,        dummy, size(x), MPI_REAL4, MPI_MAX, 0, vcomm, ierr)
       END IF
    END IF
#endif

  END SUBROUTINE vmax_r4_2d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE vmax_r8_2d(x, all)
    REAL(8), INTENT(INOUT) :: x(:,:)
    LOGICAL, INTENT(IN), OPTIONAL :: all
#ifdef F2008
    CONTIGUOUS x
#endif
    LOGICAL :: all_
    INTEGER :: dummy, ierr

#ifdef PARALLEL3D
    IF (vpes == 1) RETURN

    all_ = .FALSE.
    IF (present(all)) all_ = all

    IF (all_) THEN
       CALL mpi_allreduce(MPI_IN_PLACE, x, size(x), MPI_REAL8, MPI_MAX, vcomm, ierr)
    ELSE
       IF (vrank==0) THEN
          CALL mpi_reduce(MPI_IN_PLACE, x, size(x), MPI_REAL8, MPI_MAX, 0, vcomm, ierr)
       ELSE
          CALL mpi_reduce(x,        dummy, size(x), MPI_REAL8, MPI_MAX, 0, vcomm, ierr)
       END IF
    END IF
#endif

  END SUBROUTINE vmax_r8_2d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE vmax_i4_2d(x, all)
    INTEGER(4), INTENT(INOUT) :: x(:,:)
    LOGICAL,    INTENT(IN), OPTIONAL :: all
#ifdef F2008
    CONTIGUOUS x
#endif
    LOGICAL :: all_
    INTEGER :: dummy, ierr

#ifdef PARALLEL3D
    IF (vpes == 1) RETURN

    all_ = .FALSE.
    IF (present(all)) all_ = all

    IF (all_) THEN
       CALL mpi_allreduce(MPI_IN_PLACE, x, size(x), MPI_INTEGER4, MPI_MAX, vcomm, ierr)
    ELSE
       IF (vrank==0) THEN
          CALL mpi_reduce(MPI_IN_PLACE, x, size(x), MPI_INTEGER4, MPI_MAX, 0, vcomm, ierr)
       ELSE
          CALL mpi_reduce(x,        dummy, size(x), MPI_INTEGER4, MPI_MAX, 0, vcomm, ierr)
       END IF
    END IF
#endif

  END SUBROUTINE vmax_i4_2d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE vmax_i8_2d(x, all)
    INTEGER(8), INTENT(INOUT) :: x(:,:)
    LOGICAL,    INTENT(IN), OPTIONAL :: all
#ifdef F2008
    CONTIGUOUS x
#endif
    LOGICAL :: all_
    INTEGER :: dummy, ierr

#ifdef PARALLEL3D
    IF (vpes == 1) RETURN

    all_ = .FALSE.
    IF (present(all)) all_ = all

    IF (all_) THEN
       CALL mpi_allreduce(MPI_IN_PLACE, x, size(x), MPI_INTEGER8, MPI_MAX, vcomm, ierr)
    ELSE
       IF (vrank==0) THEN
          CALL mpi_reduce(MPI_IN_PLACE, x, size(x), MPI_INTEGER8, MPI_MAX, 0, vcomm, ierr)
       ELSE
          CALL mpi_reduce(x,        dummy, size(x), MPI_INTEGER8, MPI_MAX, 0, vcomm, ierr)
       END IF
    END IF
#endif

  END SUBROUTINE vmax_i8_2d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE vmax_r4(x, all)
    REAL(4), INTENT(INOUT) :: x
    LOGICAL, INTENT(IN), OPTIONAL :: all
    LOGICAL :: all_
    INTEGER :: dummy, ierr

#ifdef PARALLEL3D
    IF (vpes == 1) RETURN

    all_ = .FALSE.
    IF (present(all)) all_ = all

    IF (all_) THEN
       CALL mpi_allreduce(MPI_IN_PLACE, x, 1, MPI_REAL4, MPI_MAX, vcomm, ierr)
    ELSE
       IF (vrank==0) THEN
          CALL mpi_reduce(MPI_IN_PLACE, x, 1, MPI_REAL4, MPI_MAX, 0, vcomm, ierr)
       ELSE
          CALL mpi_reduce(x,        dummy, 1, MPI_REAL4, MPI_MAX, 0, vcomm, ierr)
       END IF
    END IF
#endif

  END SUBROUTINE vmax_r4

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE vmax_r8(x, all)
    REAL(8), INTENT(INOUT) :: x
    LOGICAL, INTENT(IN), OPTIONAL :: all
    LOGICAL :: all_
    INTEGER :: dummy, ierr

#ifdef PARALLEL3D
    IF (vpes == 1) RETURN

    all_ = .FALSE.
    IF (present(all)) all_ = all

    IF (all_) THEN
       CALL mpi_allreduce(MPI_IN_PLACE, x, 1, MPI_REAL8, MPI_MAX, vcomm, ierr)
    ELSE
       IF (vrank==0) THEN
          CALL mpi_reduce(MPI_IN_PLACE, x, 1, MPI_REAL8, MPI_MAX, 0, vcomm, ierr)
       ELSE
          CALL mpi_reduce(x,        dummy, 1, MPI_REAL8, MPI_MAX, 0, vcomm, ierr)
       END IF
    END IF
#endif

  END SUBROUTINE vmax_r8

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE vmin_r4_2d(x, all)
    REAL(4), INTENT(INOUT) :: x(:,:)
    LOGICAL, INTENT(IN), OPTIONAL :: all
#ifdef F2008
    CONTIGUOUS x
#endif
    LOGICAL :: all_
    INTEGER :: dummy, ierr

#ifdef PARALLEL3D
    IF (vpes == 1) RETURN

    all_ = .FALSE.
    IF (present(all)) all_ = all

    IF (all_) THEN
       CALL mpi_allreduce(MPI_IN_PLACE, x, size(x), MPI_REAL4, MPI_MIN, vcomm, ierr)
    ELSE
       IF (vrank==0) THEN
          CALL mpi_reduce(MPI_IN_PLACE, x, size(x), MPI_REAL4, MPI_MIN, 0, vcomm, ierr)
       ELSE
          CALL mpi_reduce(x,        dummy, size(x), MPI_REAL4, MPI_MIN, 0, vcomm, ierr)
       END IF
    END IF
#endif

  END SUBROUTINE vmin_r4_2d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE vmin_r8_2d(x, all)
    REAL(8), INTENT(INOUT) :: x(:,:)
    LOGICAL, INTENT(IN), OPTIONAL :: all
#ifdef F2008
    CONTIGUOUS x
#endif
    LOGICAL :: all_
    INTEGER :: dummy, ierr

#ifdef PARALLEL3D
    IF (vpes == 1) RETURN

    all_ = .FALSE.
    IF (present(all)) all_ = all

    IF (all_) THEN
       CALL mpi_allreduce(MPI_IN_PLACE, x, size(x), MPI_REAL8, MPI_MIN, vcomm, ierr)
    ELSE
       IF (vrank==0) THEN
          CALL mpi_reduce(MPI_IN_PLACE, x, size(x), MPI_REAL8, MPI_MIN, 0, vcomm, ierr)
       ELSE
          CALL mpi_reduce(x,        dummy, size(x), MPI_REAL8, MPI_MIN, 0, vcomm, ierr)
       END IF
    END IF
#endif

  END SUBROUTINE vmin_r8_2d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE vmin_i4_2d(x, all)
    INTEGER(4), INTENT(INOUT) :: x(:,:)
    LOGICAL,    INTENT(IN), OPTIONAL :: all
#ifdef F2008
    CONTIGUOUS x
#endif
    LOGICAL :: all_
    INTEGER :: dummy, ierr

#ifdef PARALLEL3D
    IF (vpes == 1) RETURN

    all_ = .FALSE.
    IF (present(all)) all_ = all

    IF (all_) THEN
       CALL mpi_allreduce(MPI_IN_PLACE, x, size(x), MPI_INTEGER4, MPI_MIN, vcomm, ierr)
    ELSE
       IF (vrank==0) THEN
          CALL mpi_reduce(MPI_IN_PLACE, x, size(x), MPI_INTEGER4, MPI_MIN, 0, vcomm, ierr)
       ELSE
          CALL mpi_reduce(x,        dummy, size(x), MPI_INTEGER4, MPI_MIN, 0, vcomm, ierr)
       END IF
    END IF
#endif

  END SUBROUTINE vmin_i4_2d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE vmin_i8_2d(x, all)
    INTEGER(8), INTENT(INOUT) :: x(:,:)
    LOGICAL,    INTENT(IN), OPTIONAL :: all
#ifdef F2008
    CONTIGUOUS x
#endif
    LOGICAL :: all_
    INTEGER :: dummy, ierr

#ifdef PARALLEL3D
    IF (vpes == 1) RETURN

    all_ = .FALSE.
    IF (present(all)) all_ = all

    IF (all_) THEN
       CALL mpi_allreduce(MPI_IN_PLACE, x, size(x), MPI_INTEGER8, MPI_MIN, vcomm, ierr)
    ELSE
       IF (vrank==0) THEN
          CALL mpi_reduce(MPI_IN_PLACE, x, size(x), MPI_INTEGER8, MPI_MIN, 0, vcomm, ierr)
       ELSE
          CALL mpi_reduce(x,        dummy, size(x), MPI_INTEGER8, MPI_MIN, 0, vcomm, ierr)
       END IF
    END IF
#endif

  END SUBROUTINE vmin_i8_2d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE vmin_r4(x, all)
    REAL(4), INTENT(INOUT) :: x
    LOGICAL, INTENT(IN), OPTIONAL :: all
    LOGICAL :: all_
    INTEGER :: dummy, ierr

#ifdef PARALLEL3D
    IF (vpes == 1) RETURN

    all_ = .FALSE.
    IF (present(all)) all_ = all

    IF (all_) THEN
       CALL mpi_allreduce(MPI_IN_PLACE, x, 1, MPI_REAL4, MPI_MIN, vcomm, ierr)
    ELSE
       IF (vrank==0) THEN
          CALL mpi_reduce(MPI_IN_PLACE, x, 1, MPI_REAL4, MPI_MIN, 0, vcomm, ierr)
       ELSE
          CALL mpi_reduce(x,        dummy, 1, MPI_REAL4, MPI_MIN, 0, vcomm, ierr)
       END IF
    END IF
#endif

  END SUBROUTINE vmin_r4

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE vmin_r8(x, all)
    REAL(8), INTENT(INOUT) :: x
    LOGICAL, INTENT(IN), OPTIONAL :: all
    LOGICAL :: all_
    INTEGER :: dummy, ierr

#ifdef PARALLEL3D
    IF (vpes == 1) RETURN

    all_ = .FALSE.
    IF (present(all)) all_ = all

    IF (all_) THEN
       CALL mpi_allreduce(MPI_IN_PLACE, x, 1, MPI_REAL8, MPI_MIN, vcomm, ierr)
    ELSE
       IF (vrank==0) THEN
          CALL mpi_reduce(MPI_IN_PLACE, x, 1, MPI_REAL8, MPI_MIN, 0, vcomm, ierr)
       ELSE
          CALL mpi_reduce(x,        dummy, 1, MPI_REAL8, MPI_MIN, 0, vcomm, ierr)
       END IF
    END IF
#endif

  END SUBROUTINE vmin_r8

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE hsum_r8_3d(x, all)
    REAL(8), INTENT(INOUT) :: x(:,:,:)
    LOGICAL, INTENT(IN), OPTIONAL :: all
#ifdef F2008
    CONTIGUOUS x
#endif
    LOGICAL :: all_
    INTEGER :: dummy, ierr

#ifdef PARALLEL_MPI
    all_ = .FALSE.
    IF (present(all)) all_ = all

    IF (all_) THEN
       CALL mpi_allreduce(MPI_IN_PLACE, x, size(x), MPI_REAL8, MPI_SUM, hcomm, ierr)
    ELSE
       IF (hrank==0) THEN
          CALL mpi_reduce(MPI_IN_PLACE, x, size(x), MPI_REAL8, MPI_SUM, 0, hcomm, ierr)
       ELSE
          CALL mpi_reduce(x,        dummy, size(x), MPI_REAL8, MPI_SUM, 0, hcomm, ierr)
       END IF
    END IF
#ifdef REDUCE_CONSISTENT
    x(:,:,:) = REAL(x, KIND=4)
#endif
#endif

  END SUBROUTINE hsum_r8_3d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE hsum_r4_3d(x, all)
    REAL(4), INTENT(INOUT) :: x(:,:,:)
    LOGICAL, INTENT(IN), OPTIONAL :: all
#ifdef F2008
    CONTIGUOUS x
#endif
    LOGICAL :: all_
    INTEGER :: dummy, ierr

#ifdef PARALLEL_MPI
    all_ = .FALSE.
    IF (present(all)) all_ = all

    IF (all_) THEN
       CALL mpi_allreduce(MPI_IN_PLACE, x, size(x), MPI_REAL4, MPI_SUM, hcomm, ierr)
    ELSE
       IF (hrank==0) THEN
          CALL mpi_reduce(MPI_IN_PLACE, x, size(x), MPI_REAL4, MPI_SUM, 0, hcomm, ierr)
       ELSE
          CALL mpi_reduce(x,        dummy, size(x), MPI_REAL4, MPI_SUM, 0, hcomm, ierr)
       END IF
    END IF
#endif

  END SUBROUTINE hsum_r4_3d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE hsum_r8_2d(x, all)
    REAL(8), INTENT(INOUT) :: x(:,:)
    LOGICAL, INTENT(IN), OPTIONAL :: all
#ifdef F2008
    CONTIGUOUS x
#endif
    LOGICAL :: all_
    INTEGER :: dummy, ierr

#ifdef PARALLEL_MPI
    all_ = .FALSE.
    IF (present(all)) all_ = all

    IF (all_) THEN
       CALL mpi_allreduce(MPI_IN_PLACE, x, size(x), MPI_REAL8, MPI_SUM, hcomm, ierr)
    ELSE
       IF (hrank==0) THEN
          CALL mpi_reduce(MPI_IN_PLACE, x, size(x), MPI_REAL8, MPI_SUM, 0, hcomm, ierr)
       ELSE
          CALL mpi_reduce(x,        dummy, size(x), MPI_REAL8, MPI_SUM, 0, hcomm, ierr)
       END IF
    END IF
#ifdef REDUCE_CONSISTENT
    x(:,:) = REAL(x, KIND=4)
#endif
#endif

  END SUBROUTINE hsum_r8_2d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE hsum_r4_2d(x, all)
    REAL(4), INTENT(INOUT) :: x(:,:)
    LOGICAL, INTENT(IN), OPTIONAL :: all
#ifdef F2008
    CONTIGUOUS x
#endif
    LOGICAL :: all_
    INTEGER :: dummy, ierr

#ifdef PARALLEL_MPI
    all_ = .FALSE.
    IF (present(all)) all_ = all

    IF (all_) THEN
       CALL mpi_allreduce(MPI_IN_PLACE, x, size(x), MPI_REAL4, MPI_SUM, hcomm, ierr)
    ELSE
       IF (hrank==0) THEN
          CALL mpi_reduce(MPI_IN_PLACE, x, size(x), MPI_REAL4, MPI_SUM, 0, hcomm, ierr)
       ELSE
          CALL mpi_reduce(x,        dummy, size(x), MPI_REAL4, MPI_SUM, 0, hcomm, ierr)
       END IF
    END IF
#endif

  END SUBROUTINE hsum_r4_2d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE hsum_r8_1d(x, all)
    REAL(8), INTENT(INOUT) :: x(:)
    LOGICAL, INTENT(IN), OPTIONAL :: all
#ifdef F2008
    CONTIGUOUS x
#endif
    LOGICAL :: all_
    INTEGER :: dummy, ierr

#ifdef PARALLEL_MPI
    all_ = .FALSE.
    IF (present(all)) all_ = all

    IF (all_) THEN
       CALL mpi_allreduce(MPI_IN_PLACE, x, size(x), MPI_REAL8, MPI_SUM, hcomm, ierr)
    ELSE
       IF (hrank==0) THEN
          CALL mpi_reduce(MPI_IN_PLACE, x, size(x), MPI_REAL8, MPI_SUM, 0, hcomm, ierr)
       ELSE
          CALL mpi_reduce(x,        dummy, size(x), MPI_REAL8, MPI_SUM, 0, hcomm, ierr)
       END IF
    END IF
#ifdef REDUCE_CONSISTENT
    x(:) = REAL(x, KIND=4)
#endif
#endif

  END SUBROUTINE hsum_r8_1d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE hsum_r4_1d(x, all)
    REAL(4), INTENT(INOUT) :: x(:)
    LOGICAL, INTENT(IN), OPTIONAL :: all
#ifdef F2008
    CONTIGUOUS x
#endif
    LOGICAL :: all_
    INTEGER :: dummy, ierr

#ifdef PARALLEL_MPI
    all_ = .FALSE.
    IF (present(all)) all_ = all

    IF (all_) THEN
       CALL mpi_allreduce(MPI_IN_PLACE, x, size(x), MPI_REAL4, MPI_SUM, hcomm, ierr)
    ELSE
       IF (hrank==0) THEN
          CALL mpi_reduce(MPI_IN_PLACE, x, size(x), MPI_REAL4, MPI_SUM, 0, hcomm, ierr)
       ELSE
          CALL mpi_reduce(x,        dummy, size(x), MPI_REAL4, MPI_SUM, 0, hcomm, ierr)
       END IF
    END IF
#endif

  END SUBROUTINE hsum_r4_1d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE hcast_r8_2d(x)
    REAL(8), INTENT(INOUT) :: x(:,:)
#ifdef F2008
    CONTIGUOUS x
#endif
    INTEGER :: ierr

#ifdef PARALLEL_MPI
    CALL mpi_bcast(x, size(x), MPI_REAL8, 0, hcomm, ierr)
#endif
  END SUBROUTINE hcast_r8_2d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE hcast_r4_2d(x)
    REAL(4), INTENT(INOUT) :: x(:,:)
#ifdef F2008
    CONTIGUOUS x
#endif
    INTEGER :: ierr

#ifdef PARALLEL_MPI
    CALL mpi_bcast(x, size(x), MPI_REAL4, 0, hcomm, ierr)
#endif
  END SUBROUTINE hcast_r4_2d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE hcast_r8_3d(x)
    REAL(8), INTENT(INOUT) :: x(:,:,:)
#ifdef F2008
    CONTIGUOUS x
#endif
    INTEGER :: ierr

#ifdef PARALLEL_MPI
    CALL mpi_bcast(x, size(x), MPI_REAL8, 0, hcomm, ierr)
#endif
  END SUBROUTINE hcast_r8_3d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE hcast_r4_3d(x)
    REAL(4), INTENT(INOUT) :: x(:,:,:)
#ifdef F2008
    CONTIGUOUS x
#endif
    INTEGER :: ierr

#ifdef PARALLEL_MPI
    CALL mpi_bcast(x, size(x), MPI_REAL4, 0, hcomm, ierr)
#endif
  END SUBROUTINE hcast_r4_3d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE hmax_r4(x, all)
    REAL(4), INTENT(INOUT) :: x
    LOGICAL, INTENT(IN), OPTIONAL :: all
    LOGICAL :: all_
    INTEGER :: dummy, ierr

#ifdef PARALLEL_MPI
    all_ = .FALSE.
    IF (present(all)) all_ = all

    IF (all_) THEN
       CALL mpi_allreduce(MPI_IN_PLACE, x, 1, MPI_REAL4, MPI_MAX, hcomm, ierr)
    ELSE
       IF (rank==0) THEN
          CALL mpi_reduce(MPI_IN_PLACE, x, 1, MPI_REAL4, MPI_MAX, 0, hcomm, ierr)
       ELSE
          CALL mpi_reduce(x,        dummy, 1, MPI_REAL4, MPI_MAX, 0, hcomm, ierr)
       END IF
    END IF
#endif

  END SUBROUTINE hmax_r4

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE hmax_r8(x, all)
    REAL(8), INTENT(INOUT) :: x
    LOGICAL, INTENT(IN), OPTIONAL :: all
    LOGICAL :: all_
    INTEGER :: dummy, ierr

#ifdef PARALLEL_MPI
    all_ = .FALSE.
    IF (present(all)) all_ = all

    IF (all_) THEN
       CALL mpi_allreduce(MPI_IN_PLACE, x, 1, MPI_REAL8, MPI_MAX, hcomm, ierr)
    ELSE
       IF (rank==0) THEN
          CALL mpi_reduce(MPI_IN_PLACE, x, 1, MPI_REAL8, MPI_MAX, 0, hcomm, ierr)
       ELSE
          CALL mpi_reduce(x,        dummy, 1, MPI_REAL8, MPI_MAX, 0, hcomm, ierr)
       END IF
    END IF
#endif

  END SUBROUTINE hmax_r8

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE hmin_r4(x, all)
    REAL(4), INTENT(INOUT) :: x
    LOGICAL, INTENT(IN), OPTIONAL :: all
    LOGICAL :: all_
    INTEGER :: dummy, ierr

#ifdef PARALLEL_MPI
    all_ = .FALSE.
    IF (present(all)) all_ = all

    IF (all_) THEN
       CALL mpi_allreduce(MPI_IN_PLACE, x, 1, MPI_REAL4, MPI_MIN, hcomm, ierr)
    ELSE
       IF (rank==0) THEN
          CALL mpi_reduce(MPI_IN_PLACE, x, 1, MPI_REAL4, MPI_MIN, 0, hcomm, ierr)
       ELSE
          CALL mpi_reduce(x,        dummy, 1, MPI_REAL4, MPI_MIN, 0, hcomm, ierr)
       END IF
    END IF
#endif

  END SUBROUTINE hmin_r4

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE hmin_r8(x, all)
    REAL(8), INTENT(INOUT) :: x
    LOGICAL, INTENT(IN), OPTIONAL :: all
    LOGICAL :: all_
    INTEGER :: dummy, ierr

#ifdef PARALLEL_MPI
    all_ = .FALSE.
    IF (present(all)) all_ = all

    IF (all_) THEN
       CALL mpi_allreduce(MPI_IN_PLACE, x, 1, MPI_REAL8, MPI_MIN, hcomm, ierr)
    ELSE
       IF (rank==0) THEN
          CALL mpi_reduce(MPI_IN_PLACE, x, 1, MPI_REAL8, MPI_MIN, 0, hcomm, ierr)
       ELSE
          CALL mpi_reduce(x,        dummy, 1, MPI_REAL8, MPI_MIN, 0, hcomm, ierr)
       END IF
    END IF
#endif

  END SUBROUTINE hmin_r8

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE gsum_i8(x, all)
    INTEGER(8), INTENT(INOUT) :: x
    LOGICAL,    INTENT(IN), OPTIONAL :: all
    LOGICAL    :: all_
    INTEGER :: dummy, ierr

#ifdef PARALLEL_MPI
    all_ = .FALSE.
    IF (present(all)) all_ = all

    IF (all_) THEN
       CALL mpi_allreduce(MPI_IN_PLACE, x, 1, MPI_INTEGER8, MPI_SUM, comm, ierr)
    ELSE
       IF (rank==0) THEN
          CALL mpi_reduce(MPI_IN_PLACE, x, 1, MPI_INTEGER8, MPI_SUM, 0, comm, ierr)
       ELSE
          CALL mpi_reduce(x,        dummy, 1, MPI_INTEGER8, MPI_SUM, 0, comm, ierr)
       END IF
    END IF
#endif

  END SUBROUTINE gsum_i8

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE gsum_i4(x, all)
    INTEGER(4), INTENT(INOUT) :: x
    LOGICAL,    INTENT(IN), OPTIONAL :: all
    LOGICAL :: all_
    INTEGER :: dummy, ierr

#ifdef PARALLEL_MPI
    all_ = .FALSE.
    IF (present(all)) all_ = all

    IF (all_) THEN
       CALL mpi_allreduce(MPI_IN_PLACE, x, 1, MPI_INTEGER4, MPI_SUM, comm, ierr)
    ELSE
       IF (rank==0) THEN
          CALL mpi_reduce(MPI_IN_PLACE, x, 1, MPI_INTEGER4, MPI_SUM, 0, comm, ierr)
       ELSE
          CALL mpi_reduce(x,        dummy, 1, MPI_INTEGER4, MPI_SUM, 0, comm, ierr)
       END IF
    END IF
#endif

  END SUBROUTINE gsum_i4

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE gsum_r8(x, all)
    REAL(8), INTENT(INOUT) :: x
    LOGICAL, INTENT(IN), OPTIONAL :: all
    LOGICAL :: all_
    INTEGER :: dummy, ierr

#ifdef PARALLEL_MPI
    all_ = .FALSE.
    IF (present(all)) all_ = all

    IF (all_) THEN
       CALL mpi_allreduce(MPI_IN_PLACE, x, 1, MPI_REAL8, MPI_SUM, comm, ierr)
    ELSE
       IF (rank==0) THEN
          CALL mpi_reduce(MPI_IN_PLACE, x, 1, MPI_REAL8, MPI_SUM, 0, comm, ierr)
       ELSE
          CALL mpi_reduce(x,        dummy, 1, MPI_REAL8, MPI_SUM, 0, comm, ierr)
       END IF
    END IF
#ifdef REDUCE_CONSISTENT
    x = REAL(x, KIND=4)
#endif
#endif

  END SUBROUTINE gsum_r8

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE gsum_r4(x, all)
    REAL(4), INTENT(INOUT) :: x
    LOGICAL, INTENT(IN), OPTIONAL :: all
    LOGICAL :: all_
    INTEGER :: dummy, ierr

#ifdef PARALLEL_MPI
    all_ = .FALSE.
    IF (present(all)) all_ = all

    IF (all_) THEN
       CALL mpi_allreduce(MPI_IN_PLACE, x, 1, MPI_REAL4, MPI_SUM, comm, ierr)
    ELSE
       IF (rank==0) THEN
          CALL mpi_reduce(MPI_IN_PLACE, x, 1, MPI_REAL4, MPI_SUM, 0, comm, ierr)
       ELSE
          CALL mpi_reduce(x,        dummy, 1, MPI_REAL4, MPI_SUM, 0, comm, ierr)
       END IF
    END IF
#endif

  END SUBROUTINE gsum_r4

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE gsum_r8_1d(x, all)
    REAL(8), INTENT(INOUT) :: x(:)
    LOGICAL, INTENT(IN), OPTIONAL :: all
#ifdef F2008
    CONTIGUOUS x
#endif
    LOGICAL :: all_
    INTEGER :: dummy, ierr

#ifdef PARALLEL_MPI
    all_ = .FALSE.
    IF (present(all)) all_ = all

    IF (all_) THEN
       CALL mpi_allreduce(MPI_IN_PLACE, x, size(x), MPI_REAL8, MPI_SUM, comm, ierr)
    ELSE
       IF (rank==0) THEN
          CALL mpi_reduce(MPI_IN_PLACE, x, size(x), MPI_REAL8, MPI_SUM, 0, comm, ierr)
       ELSE
          CALL mpi_reduce(x,        dummy, size(x), MPI_REAL8, MPI_SUM, 0, comm, ierr)
       END IF
    END IF
#ifdef REDUCE_CONSISTENT
    x(:) = REAL(x, KIND=4)
#endif
#endif

  END SUBROUTINE gsum_r8_1d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE gsum_r4_1d(x, all)
    REAL(4), INTENT(INOUT) :: x(:)
    LOGICAL, INTENT(IN), OPTIONAL :: all
#ifdef F2008
    CONTIGUOUS x
#endif
    LOGICAL :: all_
    INTEGER :: dummy, ierr

#ifdef PARALLEL_MPI
    all_ = .FALSE.
    IF (present(all)) all_ = all

    IF (all_) THEN
       CALL mpi_allreduce(MPI_IN_PLACE, x, size(x), MPI_REAL4, MPI_SUM, comm, ierr)
    ELSE
       IF (rank==0) THEN
          CALL mpi_reduce(MPI_IN_PLACE, x, size(x), MPI_REAL4, MPI_SUM, 0, comm, ierr)
       ELSE
          CALL mpi_reduce(x,        dummy, size(x), MPI_REAL4, MPI_SUM, 0, comm, ierr)
       END IF
    END IF
#endif

  END SUBROUTINE gsum_r4_1d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE gsum_i8_1d(x, all)
    INTEGER(8), INTENT(INOUT) :: x(:)
    LOGICAL,    INTENT(IN), OPTIONAL :: all
#ifdef F2008
    CONTIGUOUS x
#endif
    LOGICAL :: all_
    INTEGER :: dummy, ierr

#ifdef PARALLEL_MPI
    all_ = .FALSE.
    IF (present(all)) all_ = all

    IF (all_) THEN
       CALL mpi_allreduce(MPI_IN_PLACE, x, size(x), MPI_INTEGER8, MPI_SUM, comm, ierr)
    ELSE
       IF (rank==0) THEN
          CALL mpi_reduce(MPI_IN_PLACE, x, size(x), MPI_INTEGER8, MPI_SUM, 0, comm, ierr)
       ELSE
          CALL mpi_reduce(x,        dummy, size(x), MPI_INTEGER8, MPI_SUM, 0, comm, ierr)
       END IF
    END IF
#endif

  END SUBROUTINE gsum_i8_1d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE gsum_i4_1d(x, all)
    INTEGER(4), INTENT(INOUT) :: x(:)
    LOGICAL,    INTENT(IN), OPTIONAL :: all
#ifdef F2008
    CONTIGUOUS x
#endif
    LOGICAL :: all_
    INTEGER :: dummy, ierr

#ifdef PARALLEL_MPI
    all_ = .FALSE.
    IF (present(all)) all_ = all

    IF (all_) THEN
       CALL mpi_allreduce(MPI_IN_PLACE, x, size(x), MPI_INTEGER4, MPI_SUM, comm, ierr)
    ELSE
       IF (rank==0) THEN
          CALL mpi_reduce(MPI_IN_PLACE, x, size(x), MPI_INTEGER4, MPI_SUM, 0, comm, ierr)
       ELSE
          CALL mpi_reduce(x,        dummy, size(x), MPI_INTEGER4, MPI_SUM, 0, comm, ierr)
       END IF
    END IF
#endif

  END SUBROUTINE gsum_i4_1d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE gsum_r8_2d(x, all)
    REAL(8), INTENT(INOUT) :: x(:,:)
    LOGICAL, INTENT(IN), OPTIONAL :: all
#ifdef F2008
    CONTIGUOUS x
#endif
    LOGICAL :: all_
    INTEGER :: dummy, ierr

#ifdef PARALLEL_MPI
    all_ = .FALSE.
    IF (present(all)) all_ = all

    IF (all_) THEN
       CALL mpi_allreduce(MPI_IN_PLACE, x, size(x), MPI_REAL8, MPI_SUM, comm, ierr)
    ELSE
       IF (rank==0) THEN
          CALL mpi_reduce(MPI_IN_PLACE, x, size(x), MPI_REAL8, MPI_SUM, 0, comm, ierr)
       ELSE
          CALL mpi_reduce(x,        dummy, size(x), MPI_REAL8, MPI_SUM, 0, comm, ierr)
       END IF
    END IF
#ifdef REDUCE_CONSISTENT
    x(:,:) = REAL(x, KIND=4)
#endif
#endif

  END SUBROUTINE gsum_r8_2d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE gsum_r4_2d(x, all)
    REAL(4), INTENT(INOUT) :: x(:,:)
    LOGICAL, INTENT(IN), OPTIONAL :: all
#ifdef F2008
    CONTIGUOUS x
#endif
    LOGICAL :: all_
    INTEGER :: dummy, ierr

#ifdef PARALLEL_MPI
    all_ = .FALSE.
    IF (present(all)) all_ = all

    IF (all_) THEN
       CALL mpi_allreduce(MPI_IN_PLACE, x, size(x), MPI_REAL4, MPI_SUM, comm, ierr)
    ELSE
       IF (rank==0) THEN
          CALL mpi_reduce(MPI_IN_PLACE, x, size(x), MPI_REAL4, MPI_SUM, 0, comm, ierr)
       ELSE
          CALL mpi_reduce(x,        dummy, size(x), MPI_REAL4, MPI_SUM, 0, comm, ierr)
       END IF
    END IF
#endif

  END SUBROUTINE gsum_r4_2d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE gsum_i8_2d(x, all)
    INTEGER(8), INTENT(INOUT) :: x(:,:)
    LOGICAL,    INTENT(IN), OPTIONAL :: all
#ifdef F2008
    CONTIGUOUS x
#endif
    LOGICAL :: all_
    INTEGER :: dummy, ierr

#ifdef PARALLEL_MPI
    all_ = .FALSE.
    IF (present(all)) all_ = all

    IF (all_) THEN
       CALL mpi_allreduce(MPI_IN_PLACE, x, size(x), MPI_INTEGER8, MPI_SUM, comm, ierr)
    ELSE
       IF (rank==0) THEN
          CALL mpi_reduce(MPI_IN_PLACE, x, size(x), MPI_INTEGER8, MPI_SUM, 0, comm, ierr)
       ELSE
          CALL mpi_reduce(x,        dummy, size(x), MPI_INTEGER8, MPI_SUM, 0, comm, ierr)
       END IF
    END IF
#endif

  END SUBROUTINE gsum_i8_2d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE gsum_i4_2d(x, all)
    INTEGER(4), INTENT(INOUT) :: x(:,:)
    LOGICAL,    INTENT(IN), OPTIONAL :: all
#ifdef F2008
    CONTIGUOUS x
#endif
    LOGICAL :: all_
    INTEGER :: dummy, ierr

#ifdef PARALLEL_MPI
    all_ = .FALSE.
    IF (present(all)) all_ = all

    IF (all_) THEN
       CALL mpi_allreduce(MPI_IN_PLACE, x, size(x), MPI_INTEGER4, MPI_SUM, comm, ierr)
    ELSE
       IF (rank==0) THEN
          CALL mpi_reduce(MPI_IN_PLACE, x, size(x), MPI_INTEGER4, MPI_SUM, 0, comm, ierr)
       ELSE
          CALL mpi_reduce(x,        dummy, size(x), MPI_INTEGER4, MPI_SUM, 0, comm, ierr)
       END IF
    END IF
#endif

  END SUBROUTINE gsum_i4_2d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE gmax_r4(x, all)
    REAL(4), INTENT(INOUT) :: x
    LOGICAL, INTENT(IN), OPTIONAL :: all
    LOGICAL :: all_
    INTEGER :: dummy, ierr

#ifdef PARALLEL_MPI
    all_ = .FALSE.
    IF (present(all)) all_ = all

    IF (all_) THEN
       CALL mpi_allreduce(MPI_IN_PLACE, x, 1, MPI_REAL4, MPI_MAX, comm, ierr)
    ELSE
       IF (rank==0) THEN
          CALL mpi_reduce(MPI_IN_PLACE, x, 1, MPI_REAL4, MPI_MAX, 0, comm, ierr)
       ELSE
          CALL mpi_reduce(x,        dummy, 1, MPI_REAL4, MPI_MAX, 0, comm, ierr)
       END IF
    END IF
#endif

  END SUBROUTINE gmax_r4

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE gmax_r8(x, all)
    REAL(8), INTENT(INOUT) :: x
    LOGICAL, INTENT(IN), OPTIONAL :: all
    LOGICAL :: all_
    INTEGER :: dummy, ierr

#ifdef PARALLEL_MPI
    all_ = .FALSE.
    IF (present(all)) all_ = all

    IF (all_) THEN
       CALL mpi_allreduce(MPI_IN_PLACE, x, 1, MPI_REAL8, MPI_MAX, comm, ierr)
    ELSE
       IF (rank==0) THEN
          CALL mpi_reduce(MPI_IN_PLACE, x, 1, MPI_REAL8, MPI_MAX, 0, comm, ierr)
       ELSE
          CALL mpi_reduce(x,        dummy, 1, MPI_REAL8, MPI_MAX, 0, comm, ierr)
       END IF
    END IF
#endif

  END SUBROUTINE gmax_r8

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE gmax_i4(x, all)
    INTEGER(4), INTENT(INOUT) :: x
    LOGICAL,    INTENT(IN), OPTIONAL :: all
    LOGICAL :: all_
    INTEGER :: dummy, ierr

#ifdef PARALLEL_MPI
    all_ = .FALSE.
    IF (present(all)) all_ = all

    IF (all_) THEN
       CALL mpi_allreduce(MPI_IN_PLACE, x, 1, MPI_INTEGER4, MPI_MAX, comm, ierr)
    ELSE
       IF (rank==0) THEN
          CALL mpi_reduce(MPI_IN_PLACE, x, 1, MPI_INTEGER4, MPI_MAX, 0, comm, ierr)
       ELSE
          CALL mpi_reduce(x,        dummy, 1, MPI_INTEGER4, MPI_MAX, 0, comm, ierr)
       END IF
    END IF
#endif

  END SUBROUTINE gmax_i4

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE gmax_i8(x, all)
    INTEGER(8), INTENT(INOUT) :: x
    LOGICAL,    INTENT(IN), OPTIONAL :: all
    LOGICAL :: all_
    INTEGER :: dummy, ierr

#ifdef PARALLEL_MPI
    all_ = .FALSE.
    IF (present(all)) all_ = all

    IF (all_) THEN
       CALL mpi_allreduce(MPI_IN_PLACE, x, 1, MPI_INTEGER8, MPI_MAX, comm, ierr)
    ELSE
       IF (rank==0) THEN
          CALL mpi_reduce(MPI_IN_PLACE, x, 1, MPI_INTEGER8, MPI_MAX, 0, comm, ierr)
       ELSE
          CALL mpi_reduce(x,        dummy, 1, MPI_INTEGER8, MPI_MAX, 0, comm, ierr)
       END IF
    END IF
#endif

  END SUBROUTINE gmax_i8

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE gmax_r4_1d(x, all)
    REAL(4), INTENT(INOUT) :: x(:)
    LOGICAL, INTENT(IN), OPTIONAL :: all
#ifdef F2008
    CONTIGUOUS x
#endif
    LOGICAL :: all_
    INTEGER :: dummy, ierr

#ifdef PARALLEL_MPI
    all_ = .FALSE.
    IF (present(all)) all_ = all

    IF (all_) THEN
       CALL mpi_allreduce(MPI_IN_PLACE, x, size(x), MPI_REAL4, MPI_MAX, comm, ierr)
    ELSE
       IF (rank==0) THEN
          CALL mpi_reduce(MPI_IN_PLACE, x, size(x), MPI_REAL4, MPI_MAX, 0, comm, ierr)
       ELSE
          CALL mpi_reduce(x,        dummy, size(x), MPI_REAL4, MPI_MAX, 0, comm, ierr)
       END IF
    END IF
#endif

  END SUBROUTINE gmax_r4_1d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE gmax_r8_1d(x, all)
    REAL(8), INTENT(INOUT) :: x(:)
    LOGICAL, INTENT(IN), OPTIONAL :: all
#ifdef F2008
    CONTIGUOUS x
#endif
    LOGICAL :: all_
    INTEGER :: dummy, ierr

#ifdef PARALLEL_MPI
    all_ = .FALSE.
    IF (present(all)) all_ = all

    IF (all_) THEN
       CALL mpi_allreduce(MPI_IN_PLACE, x, size(x), MPI_REAL8, MPI_MAX, comm, ierr)
    ELSE
       IF (rank==0) THEN
          CALL mpi_reduce(MPI_IN_PLACE, x, size(x), MPI_REAL8, MPI_MAX, 0, comm, ierr)
       ELSE
          CALL mpi_reduce(x,        dummy, size(x), MPI_REAL8, MPI_MAX, 0, comm, ierr)
       END IF
    END IF
#endif

  END SUBROUTINE gmax_r8_1d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE gmax_i4_1d(x, all)
    INTEGER(4), INTENT(INOUT) :: x(:)
    LOGICAL,    INTENT(IN), OPTIONAL :: all
    LOGICAL :: all_
    INTEGER :: dummy, ierr

#ifdef PARALLEL_MPI
    all_ = .FALSE.
    IF (present(all)) all_ = all

    IF (all_) THEN
       CALL mpi_allreduce(MPI_IN_PLACE, x, size(x), MPI_INTEGER4, MPI_MAX, comm, ierr)
    ELSE
       IF (rank==0) THEN
          CALL mpi_reduce(MPI_IN_PLACE, x, size(x), MPI_INTEGER4, MPI_MAX, 0, comm, ierr)
       ELSE
          CALL mpi_reduce(x,        dummy, size(x), MPI_INTEGER4, MPI_MAX, 0, comm, ierr)
       END IF
    END IF
#endif

  END SUBROUTINE gmax_i4_1d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE gmax_i8_1d(x, all)
    INTEGER(8), INTENT(INOUT) :: x(:)
    LOGICAL,    INTENT(IN), OPTIONAL :: all
    LOGICAL :: all_
    INTEGER :: dummy, ierr

#ifdef PARALLEL_MPI
    all_ = .FALSE.
    IF (present(all)) all_ = all

    IF (all_) THEN
       CALL mpi_allreduce(MPI_IN_PLACE, x, size(x), MPI_INTEGER8, MPI_MAX, comm, ierr)
    ELSE
       IF (rank==0) THEN
          CALL mpi_reduce(MPI_IN_PLACE, x, size(x), MPI_INTEGER8, MPI_MAX, 0, comm, ierr)
       ELSE
          CALL mpi_reduce(x,        dummy, size(x), MPI_INTEGER8, MPI_MAX, 0, comm, ierr)
       END IF
    END IF
#endif

  END SUBROUTINE gmax_i8_1d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE gmax_r4_2d(x, all)
    REAL(4), INTENT(INOUT) :: x(:,:)
    LOGICAL, INTENT(IN), OPTIONAL :: all
#ifdef F2008
    CONTIGUOUS x
#endif
    LOGICAL :: all_
    INTEGER :: dummy, ierr

#ifdef PARALLEL_MPI
    all_ = .FALSE.
    IF (present(all)) all_ = all

    IF (all_) THEN
       CALL mpi_allreduce(MPI_IN_PLACE, x, size(x), MPI_REAL4, MPI_MAX, comm, ierr)
    ELSE
       IF (rank==0) THEN
          CALL mpi_reduce(MPI_IN_PLACE, x, size(x), MPI_REAL4, MPI_MAX, 0, comm, ierr)
       ELSE
          CALL mpi_reduce(x,        dummy, size(x), MPI_REAL4, MPI_MAX, 0, comm, ierr)
       END IF
    END IF
#endif

  END SUBROUTINE gmax_r4_2d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE gmax_r8_2d(x, all)
    REAL(8), INTENT(INOUT) :: x(:,:)
    LOGICAL, INTENT(IN), OPTIONAL :: all
#ifdef F2008
    CONTIGUOUS x
#endif
    LOGICAL :: all_
    INTEGER :: dummy, ierr

#ifdef PARALLEL_MPI
    all_ = .FALSE.
    IF (present(all)) all_ = all

    IF (all_) THEN
       CALL mpi_allreduce(MPI_IN_PLACE, x, size(x), MPI_REAL8, MPI_MAX, comm, ierr)
    ELSE
       IF (rank==0) THEN
          CALL mpi_reduce(MPI_IN_PLACE, x, size(x), MPI_REAL8, MPI_MAX, 0, comm, ierr)
       ELSE
          CALL mpi_reduce(x,        dummy, size(x), MPI_REAL8, MPI_MAX, 0, comm, ierr)
       END IF
    END IF
#endif

  END SUBROUTINE gmax_r8_2d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE gmax_i4_2d(x, all)
    INTEGER(4), INTENT(INOUT) :: x(:,:)
    LOGICAL,    INTENT(IN), OPTIONAL :: all
    LOGICAL :: all_
    INTEGER :: dummy, ierr

#ifdef PARALLEL_MPI
    all_ = .FALSE.
    IF (present(all)) all_ = all

    IF (all_) THEN
       CALL mpi_allreduce(MPI_IN_PLACE, x, size(x), MPI_INTEGER4, MPI_MAX, comm, ierr)
    ELSE
       IF (rank==0) THEN
          CALL mpi_reduce(MPI_IN_PLACE, x, size(x), MPI_INTEGER4, MPI_MAX, 0, comm, ierr)
       ELSE
          CALL mpi_reduce(x,        dummy, size(x), MPI_INTEGER4, MPI_MAX, 0, comm, ierr)
       END IF
    END IF
#endif

  END SUBROUTINE gmax_i4_2d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE gmax_i8_2d(x, all)
    INTEGER(8), INTENT(INOUT) :: x(:,:)
    LOGICAL,    INTENT(IN), OPTIONAL :: all
    LOGICAL :: all_
    INTEGER :: dummy, ierr

#ifdef PARALLEL_MPI
    all_ = .FALSE.
    IF (present(all)) all_ = all

    IF (all_) THEN
       CALL mpi_allreduce(MPI_IN_PLACE, x, size(x), MPI_INTEGER8, MPI_MAX, comm, ierr)
    ELSE
       IF (rank==0) THEN
          CALL mpi_reduce(MPI_IN_PLACE, x, size(x), MPI_INTEGER8, MPI_MAX, 0, comm, ierr)
       ELSE
          CALL mpi_reduce(x,        dummy, size(x), MPI_INTEGER8, MPI_MAX, 0, comm, ierr)
       END IF
    END IF
#endif

  END SUBROUTINE gmax_i8_2d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE gmin_r4(x, all)
    REAL(4), INTENT(INOUT) :: x
    LOGICAL, INTENT(IN), OPTIONAL :: all
    LOGICAL :: all_
    INTEGER :: dummy, ierr

#ifdef PARALLEL_MPI
    all_ = .FALSE.
    IF (present(all)) all_ = all

    IF (all_) THEN
       CALL mpi_allreduce(MPI_IN_PLACE, x, 1, MPI_REAL4, MPI_MIN, comm, ierr)
    ELSE
       IF (rank==0) THEN
          CALL mpi_reduce(MPI_IN_PLACE, x, 1, MPI_REAL4, MPI_MIN, 0, comm, ierr)
       ELSE
          CALL mpi_reduce(x,        dummy, 1, MPI_REAL4, MPI_MIN, 0, comm, ierr)
       END IF
    END IF
#endif

  END SUBROUTINE gmin_r4

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE gmin_r8(x, all)
    REAL(8), INTENT(INOUT) :: x
    LOGICAL, INTENT(IN), OPTIONAL :: all
    LOGICAL :: all_
    INTEGER :: dummy, ierr

#ifdef PARALLEL_MPI
    all_ = .FALSE.
    IF (present(all)) all_ = all

    IF (all_) THEN
       CALL mpi_allreduce(MPI_IN_PLACE, x, 1, MPI_REAL8, MPI_MIN, comm, ierr)
    ELSE
       IF (rank==0) THEN
          CALL mpi_reduce(MPI_IN_PLACE, x, 1, MPI_REAL8, MPI_MIN, 0, comm, ierr)
       ELSE
          CALL mpi_reduce(x,        dummy, 1, MPI_REAL8, MPI_MIN, 0, comm, ierr)
       END IF
    END IF
#endif

  END SUBROUTINE gmin_r8

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE gmin_i4(x, all)
    INTEGER(4), INTENT(INOUT) :: x
    LOGICAL,    INTENT(IN), OPTIONAL :: all
    LOGICAL :: all_
    INTEGER :: dummy, ierr

#ifdef PARALLEL_MPI
    all_ = .FALSE.
    IF (present(all)) all_ = all

    IF (all_) THEN
       CALL mpi_allreduce(MPI_IN_PLACE, x, 1, MPI_INTEGER4, MPI_MIN, comm, ierr)
    ELSE
       IF (rank==0) THEN
          CALL mpi_reduce(MPI_IN_PLACE, x, 1, MPI_INTEGER4, MPI_MIN, 0, comm, ierr)
       ELSE
          CALL mpi_reduce(x,        dummy, 1, MPI_INTEGER4, MPI_MIN, 0, comm, ierr)
       END IF
    END IF
#endif

  END SUBROUTINE gmin_i4

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE gmin_i8(x, all)
    INTEGER(8), INTENT(INOUT) :: x
    LOGICAL,    INTENT(IN), OPTIONAL :: all
    LOGICAL :: all_
    INTEGER :: dummy, ierr

#ifdef PARALLEL_MPI
    all_ = .FALSE.
    IF (present(all)) all_ = all

    IF (all_) THEN
       CALL mpi_allreduce(MPI_IN_PLACE, x, 1, MPI_INTEGER8, MPI_MIN, comm, ierr)
    ELSE
       IF (rank==0) THEN
          CALL mpi_reduce(MPI_IN_PLACE, x, 1, MPI_INTEGER8, MPI_MIN, 0, comm, ierr)
       ELSE
          CALL mpi_reduce(x,        dummy, 1, MPI_INTEGER8, MPI_MIN, 0, comm, ierr)
       END IF
    END IF
#endif

  END SUBROUTINE gmin_i8

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE gmin_r4_1d(x, all)
    REAL(4), INTENT(INOUT) :: x(:)
    LOGICAL, INTENT(IN), OPTIONAL :: all
#ifdef F2008
    CONTIGUOUS x
#endif
    LOGICAL :: all_
    INTEGER :: dummy, ierr

#ifdef PARALLEL_MPI
    all_ = .FALSE.
    IF (present(all)) all_ = all

    IF (all_) THEN
       CALL mpi_allreduce(MPI_IN_PLACE, x, size(x), MPI_REAL4, MPI_MIN, comm, ierr)
    ELSE
       IF (rank==0) THEN
          CALL mpi_reduce(MPI_IN_PLACE, x, size(x), MPI_REAL4, MPI_MIN, 0, comm, ierr)
       ELSE
          CALL mpi_reduce(x,        dummy, size(x), MPI_REAL4, MPI_MIN, 0, comm, ierr)
       END IF
    END IF
#endif

  END SUBROUTINE gmin_r4_1d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE gmin_r8_1d(x, all)
    REAL(8), INTENT(INOUT) :: x(:)
    LOGICAL, INTENT(IN), OPTIONAL :: all
#ifdef F2008
    CONTIGUOUS x
#endif
    LOGICAL :: all_
    INTEGER :: dummy, ierr

#ifdef PARALLEL_MPI
    all_ = .FALSE.
    IF (present(all)) all_ = all

    IF (all_) THEN
       CALL mpi_allreduce(MPI_IN_PLACE, x, size(x), MPI_REAL8, MPI_MIN, comm, ierr)
    ELSE
       IF (rank==0) THEN
          CALL mpi_reduce(MPI_IN_PLACE, x, size(x), MPI_REAL8, MPI_MIN, 0, comm, ierr)
       ELSE
          CALL mpi_reduce(x,        dummy, size(x), MPI_REAL8, MPI_MIN, 0, comm, ierr)
       END IF
    END IF
#endif

  END SUBROUTINE gmin_r8_1d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE gmin_i4_1d(x, all)
    INTEGER(4), INTENT(INOUT) :: x(:)
    LOGICAL,    INTENT(IN), OPTIONAL :: all
    LOGICAL :: all_
    INTEGER :: dummy, ierr

#ifdef PARALLEL_MPI
    all_ = .FALSE.
    IF (present(all)) all_ = all

    IF (all_) THEN
       CALL mpi_allreduce(MPI_IN_PLACE, x, size(x), MPI_INTEGER4, MPI_MIN, comm, ierr)
    ELSE
       IF (rank==0) THEN
          CALL mpi_reduce(MPI_IN_PLACE, x, size(x), MPI_INTEGER4, MPI_MIN, 0, comm, ierr)
       ELSE
          CALL mpi_reduce(x,        dummy, size(x), MPI_INTEGER4, MPI_MIN, 0, comm, ierr)
       END IF
    END IF
#endif

  END SUBROUTINE gmin_i4_1d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE gmin_i8_1d(x, all)
    INTEGER(8), INTENT(INOUT) :: x(:)
    LOGICAL,    INTENT(IN), OPTIONAL :: all
    LOGICAL :: all_
    INTEGER :: dummy, ierr

#ifdef PARALLEL_MPI
    all_ = .FALSE.
    IF (present(all)) all_ = all

    IF (all_) THEN
       CALL mpi_allreduce(MPI_IN_PLACE, x, size(x), MPI_INTEGER8, MPI_MIN, comm, ierr)
    ELSE
       IF (rank==0) THEN
          CALL mpi_reduce(MPI_IN_PLACE, x, size(x), MPI_INTEGER8, MPI_MIN, 0, comm, ierr)
       ELSE
          CALL mpi_reduce(x,        dummy, size(x), MPI_INTEGER8, MPI_MIN, 0, comm, ierr)
       END IF
    END IF
#endif

  END SUBROUTINE gmin_i8_1d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE gmin_r4_2d(x, all)
    REAL(4), INTENT(INOUT) :: x(:,:)
    LOGICAL, INTENT(IN), OPTIONAL :: all
#ifdef F2008
    CONTIGUOUS x
#endif
    LOGICAL :: all_
    INTEGER :: dummy, ierr

#ifdef PARALLEL_MPI
    all_ = .FALSE.
    IF (present(all)) all_ = all

    IF (all_) THEN
       CALL mpi_allreduce(MPI_IN_PLACE, x, size(x), MPI_REAL4, MPI_MIN, comm, ierr)
    ELSE
       IF (rank==0) THEN
          CALL mpi_reduce(MPI_IN_PLACE, x, size(x), MPI_REAL4, MPI_MIN, 0, comm, ierr)
       ELSE
          CALL mpi_reduce(x,        dummy, size(x), MPI_REAL4, MPI_MIN, 0, comm, ierr)
       END IF
    END IF
#endif

  END SUBROUTINE gmin_r4_2d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE gmin_r8_2d(x, all)
    REAL(8), INTENT(INOUT) :: x(:,:)
    LOGICAL, INTENT(IN), OPTIONAL :: all
#ifdef F2008
    CONTIGUOUS x
#endif
    LOGICAL :: all_
    INTEGER :: dummy, ierr

#ifdef PARALLEL_MPI
    all_ = .FALSE.
    IF (present(all)) all_ = all

    IF (all_) THEN
       CALL mpi_allreduce(MPI_IN_PLACE, x, size(x), MPI_REAL8, MPI_MIN, comm, ierr)
    ELSE
       IF (rank==0) THEN
          CALL mpi_reduce(MPI_IN_PLACE, x, size(x), MPI_REAL8, MPI_MIN, 0, comm, ierr)
       ELSE
          CALL mpi_reduce(x,        dummy, size(x), MPI_REAL8, MPI_MIN, 0, comm, ierr)
       END IF
    END IF
#endif

  END SUBROUTINE gmin_r8_2d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE gmin_i4_2d(x, all)
    INTEGER(4), INTENT(INOUT) :: x(:,:)
    LOGICAL,    INTENT(IN), OPTIONAL :: all
    LOGICAL :: all_
    INTEGER :: dummy, ierr

#ifdef PARALLEL_MPI
    all_ = .FALSE.
    IF (present(all)) all_ = all

    IF (all_) THEN
       CALL mpi_allreduce(MPI_IN_PLACE, x, size(x), MPI_INTEGER4, MPI_MIN, comm, ierr)
    ELSE
       IF (rank==0) THEN
          CALL mpi_reduce(MPI_IN_PLACE, x, size(x), MPI_INTEGER4, MPI_MIN, 0, comm, ierr)
       ELSE
          CALL mpi_reduce(x,        dummy, size(x), MPI_INTEGER4, MPI_MIN, 0, comm, ierr)
       END IF
    END IF
#endif

  END SUBROUTINE gmin_i4_2d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE gmin_i8_2d(x, all)
    INTEGER(8), INTENT(INOUT) :: x(:,:)
    LOGICAL,    INTENT(IN), OPTIONAL :: all
    LOGICAL :: all_
    INTEGER :: dummy, ierr

#ifdef PARALLEL_MPI
    all_ = .FALSE.
    IF (present(all)) all_ = all

    IF (all_) THEN
       CALL mpi_allreduce(MPI_IN_PLACE, x, size(x), MPI_INTEGER8, MPI_MIN, comm, ierr)
    ELSE
       IF (rank==0) THEN
          CALL mpi_reduce(MPI_IN_PLACE, x, size(x), MPI_INTEGER8, MPI_MIN, 0, comm, ierr)
       ELSE
          CALL mpi_reduce(x,        dummy, size(x), MPI_INTEGER8, MPI_MIN, 0, comm, ierr)
       END IF
    END IF
#endif

  END SUBROUTINE gmin_i8_2d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE usend_r4_3d(x, tag)
    REAL(4), INTENT(INOUT) :: x(:,:,:)
    INTEGER, INTENT(IN), OPTIONAL :: tag
#ifdef F2008
    CONTIGUOUS x
#endif
    INTEGER :: itag, ierr

#ifdef PARALLEL3D
    IF (rank_u == MPI_PROC_NULL) RETURN

    itag = 0
    IF (present(tag)) itag = tag

    CALL mpi_send(x, size(x), MPI_REAL4, rank_u, itag, comm, ierr)
#endif

  END SUBROUTINE usend_r4_3d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE usend_r8_3d(x, tag)
    REAL(8), INTENT(INOUT) :: x(:,:,:)
    INTEGER, INTENT(IN), OPTIONAL :: tag
#ifdef F2008
    CONTIGUOUS x
#endif
    INTEGER :: itag, ierr

#ifdef PARALLEL3D
    IF (rank_u == MPI_PROC_NULL) RETURN

    itag = 0
    IF (present(tag)) itag = tag

    CALL mpi_send(x, size(x), MPI_REAL8, rank_u, itag, comm, ierr)
#endif

  END SUBROUTINE usend_r8_3d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE usend_r4_2d(x, tag)
    REAL(4), INTENT(INOUT) :: x(:,:)
    INTEGER, INTENT(IN), OPTIONAL :: tag
#ifdef F2008
    CONTIGUOUS x
#endif
    INTEGER :: itag, ierr

#ifdef PARALLEL3D
    IF (rank_u == MPI_PROC_NULL) RETURN

    itag = 0
    IF (present(tag)) itag = tag

    CALL mpi_send(x, size(x), MPI_REAL4, rank_u, itag, comm, ierr)
#endif

  END SUBROUTINE usend_r4_2d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE usend_r8_2d(x, tag)
    REAL(8), INTENT(INOUT) :: x(:,:)
    INTEGER, INTENT(IN), OPTIONAL :: tag
#ifdef F2008
    CONTIGUOUS x
#endif
    INTEGER :: itag, ierr

#ifdef PARALLEL3D
    IF (rank_u == MPI_PROC_NULL) RETURN

    itag = 0
    IF (present(tag)) itag = tag

    CALL mpi_send(x, size(x), MPI_REAL8, rank_u, itag, comm, ierr)
#endif

  END SUBROUTINE usend_r8_2d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE usend_r4_1d(x, tag)
    REAL(4), INTENT(INOUT) :: x(:)
    INTEGER, INTENT(IN), OPTIONAL :: tag
#ifdef F2008
    CONTIGUOUS x
#endif
    INTEGER :: itag, ierr

#ifdef PARALLEL3D
    IF (rank_u == MPI_PROC_NULL) RETURN

    itag = 0
    IF (present(tag)) itag = tag

    CALL mpi_send(x, size(x), MPI_REAL4, rank_u, itag, comm, ierr)
#endif

  END SUBROUTINE usend_r4_1d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE usend_r8_1d(x, tag)
    REAL(8), INTENT(INOUT) :: x(:)
    INTEGER, INTENT(IN), OPTIONAL :: tag
#ifdef F2008
    CONTIGUOUS x
#endif
    INTEGER :: itag, ierr

#ifdef PARALLEL3D
    IF (rank_u == MPI_PROC_NULL) RETURN

    itag = 0
    IF (present(tag)) itag = tag

    CALL mpi_send(x, size(x), MPI_REAL8, rank_u, itag, comm, ierr)
#endif

  END SUBROUTINE usend_r8_1d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE usend_r4(x, tag)
    REAL(4), INTENT(INOUT) :: x
    INTEGER, INTENT(IN), OPTIONAL :: tag
    INTEGER :: itag, ierr

#ifdef PARALLEL3D
    IF (rank_u == MPI_PROC_NULL) RETURN

    itag = 0
    IF (present(tag)) itag = tag

    CALL mpi_send(x, 1, MPI_REAL4, rank_u, itag, comm, ierr)
#endif

  END SUBROUTINE usend_r4

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE usend_r8(x, tag)
    REAL(8), INTENT(INOUT) :: x
    INTEGER, INTENT(IN), OPTIONAL :: tag
    INTEGER :: itag, ierr

#ifdef PARALLEL3D
    IF (rank_u == MPI_PROC_NULL) RETURN

    itag = 0
    IF (present(tag)) itag = tag

    CALL mpi_send(x, 1, MPI_REAL8, rank_u, itag, comm, ierr)
#endif

  END SUBROUTINE usend_r8

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE lsend_r4_3d(x, tag)
    REAL(4), INTENT(INOUT) :: x(:,:,:)
    INTEGER, INTENT(IN), OPTIONAL :: tag
#ifdef F2008
    CONTIGUOUS x
#endif
    INTEGER :: itag, ierr

#ifdef PARALLEL3D
    IF (rank_l == MPI_PROC_NULL) RETURN

    itag = 0
    IF (present(tag)) itag = tag

    CALL mpi_send(x, size(x), MPI_REAL4, rank_l, itag, comm, ierr)
#endif

  END SUBROUTINE lsend_r4_3d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE lsend_r8_3d(x, tag)
    REAL(8), INTENT(INOUT) :: x(:,:,:)
    INTEGER, INTENT(IN), OPTIONAL :: tag
#ifdef F2008
    CONTIGUOUS x
#endif
    INTEGER :: itag, ierr

#ifdef PARALLEL3D
    IF (rank_l == MPI_PROC_NULL) RETURN

    itag = 0
    IF (present(tag)) itag = tag

    CALL mpi_send(x, size(x), MPI_REAL8, rank_l, itag, comm, ierr)
#endif

  END SUBROUTINE lsend_r8_3d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE lsend_r4_2d(x, tag)
    REAL(4), INTENT(INOUT) :: x(:,:)
    INTEGER, INTENT(IN), OPTIONAL :: tag
#ifdef F2008
    CONTIGUOUS x
#endif
    INTEGER :: itag, ierr

#ifdef PARALLEL3D
    IF (rank_l == MPI_PROC_NULL) RETURN

    itag = 0
    IF (present(tag)) itag = tag

    CALL mpi_send(x, size(x), MPI_REAL4, rank_l, itag, comm, ierr)
#endif

  END SUBROUTINE lsend_r4_2d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE lsend_r8_2d(x, tag)
    REAL(8), INTENT(INOUT) :: x(:,:)
    INTEGER, INTENT(IN), OPTIONAL :: tag
#ifdef F2008
    CONTIGUOUS x
#endif
    INTEGER :: itag, ierr

#ifdef PARALLEL3D
    IF (rank_l == MPI_PROC_NULL) RETURN

    itag = 0
    IF (present(tag)) itag = tag

    CALL mpi_send(x, size(x), MPI_REAL8, rank_l, itag, comm, ierr)
#endif

  END SUBROUTINE lsend_r8_2d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE lsend_r4_1d(x, tag)
    REAL(4), INTENT(INOUT) :: x(:)
    INTEGER, INTENT(IN), OPTIONAL :: tag
#ifdef F2008
    CONTIGUOUS x
#endif
    INTEGER :: itag, ierr

#ifdef PARALLEL3D
    IF (rank_l == MPI_PROC_NULL) RETURN

    itag = 0
    IF (present(tag)) itag = tag

    CALL mpi_send(x, size(x), MPI_REAL4, rank_l, itag, comm, ierr)
#endif

  END SUBROUTINE lsend_r4_1d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE lsend_r8_1d(x, tag)
    REAL(8), INTENT(INOUT) :: x(:)
    INTEGER, INTENT(IN), OPTIONAL :: tag
#ifdef F2008
    CONTIGUOUS x
#endif
    INTEGER :: itag, ierr

#ifdef PARALLEL3D
    IF (rank_l == MPI_PROC_NULL) RETURN

    itag = 0
    IF (present(tag)) itag = tag

    CALL mpi_send(x, size(x), MPI_REAL8, rank_l, itag, comm, ierr)
#endif

  END SUBROUTINE lsend_r8_1d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE lsend_r4(x, tag)
    REAL(4), INTENT(INOUT) :: x
    INTEGER, INTENT(IN), OPTIONAL :: tag
    INTEGER :: itag, ierr

#ifdef PARALLEL3D
    IF (rank_l == MPI_PROC_NULL) RETURN

    itag = 0
    IF (present(tag)) itag = tag

    CALL mpi_send(x, 1, MPI_REAL4, rank_l, itag, comm, ierr)
#endif

  END SUBROUTINE lsend_r4

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE lsend_r8(x, tag)
    REAL(8), INTENT(INOUT) :: x
    INTEGER, INTENT(IN), OPTIONAL :: tag
    INTEGER :: itag, ierr

#ifdef PARALLEL3D
    IF (rank_l == MPI_PROC_NULL) RETURN

    itag = 0
    IF (present(tag)) itag = tag

    CALL mpi_send(x, 1, MPI_REAL8, rank_l, itag, comm, ierr)
#endif

  END SUBROUTINE lsend_r8

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE urecv_r4_3d(x, tag)
    REAL(4), INTENT(INOUT) :: x(:,:,:)
    INTEGER, INTENT(IN), OPTIONAL :: tag
#ifdef F2008
    CONTIGUOUS x
#endif
    INTEGER :: itag, ierr

#ifdef PARALLEL3D
    IF (rank_u == MPI_PROC_NULL) RETURN

    itag = 0
    IF (present(tag)) itag = tag

    CALL mpi_recv(x, size(x), MPI_REAL4, rank_u, itag, comm, MPI_STATUS_IGNORE, ierr)
#endif

  END SUBROUTINE urecv_r4_3d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE urecv_r8_3d(x, tag)
    REAL(8), INTENT(INOUT) :: x(:,:,:)
    INTEGER, INTENT(IN), OPTIONAL :: tag
#ifdef F2008
    CONTIGUOUS x
#endif
    INTEGER :: itag, ierr

#ifdef PARALLEL3D
    IF (rank_u == MPI_PROC_NULL) RETURN

    itag = 0
    IF (present(tag)) itag = tag

    CALL mpi_recv(x, size(x), MPI_REAL8, rank_u, itag, comm, MPI_STATUS_IGNORE, ierr)
#endif

  END SUBROUTINE urecv_r8_3d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE urecv_r4_2d(x, tag)
    REAL(4), INTENT(INOUT) :: x(:,:)
    INTEGER, INTENT(IN), OPTIONAL :: tag
#ifdef F2008
    CONTIGUOUS x
#endif
    INTEGER :: itag, ierr

#ifdef PARALLEL3D
    IF (rank_u == MPI_PROC_NULL) RETURN

    itag = 0
    IF (present(tag)) itag = tag

    CALL mpi_recv(x, size(x), MPI_REAL4, rank_u, itag, comm, MPI_STATUS_IGNORE, ierr)
#endif

  END SUBROUTINE urecv_r4_2d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE urecv_r8_2d(x, tag)
    REAL(8), INTENT(INOUT) :: x(:,:)
    INTEGER, INTENT(IN), OPTIONAL :: tag
#ifdef F2008
    CONTIGUOUS x
#endif
    INTEGER :: itag, ierr

#ifdef PARALLEL3D
    IF (rank_u == MPI_PROC_NULL) RETURN

    itag = 0
    IF (present(tag)) itag = tag

    CALL mpi_recv(x, size(x), MPI_REAL8, rank_u, itag, comm, MPI_STATUS_IGNORE, ierr)
#endif

  END SUBROUTINE urecv_r8_2d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE urecv_r4_1d(x, tag)
    REAL(4), INTENT(INOUT) :: x(:)
    INTEGER, INTENT(IN), OPTIONAL :: tag
#ifdef F2008
    CONTIGUOUS x
#endif
    INTEGER :: itag, ierr

#ifdef PARALLEL3D
    IF (rank_u == MPI_PROC_NULL) RETURN

    itag = 0
    IF (present(tag)) itag = tag

    CALL mpi_recv(x, size(x), MPI_REAL4, rank_u, itag, comm, MPI_STATUS_IGNORE, ierr)
#endif

  END SUBROUTINE urecv_r4_1d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE urecv_r8_1d(x, tag)
    REAL(8), INTENT(INOUT) :: x(:)
    INTEGER, INTENT(IN), OPTIONAL :: tag
#ifdef F2008
    CONTIGUOUS x
#endif
    INTEGER :: itag, ierr

#ifdef PARALLEL3D
    IF (rank_u == MPI_PROC_NULL) RETURN

    itag = 0
    IF (present(tag)) itag = tag

    CALL mpi_recv(x, size(x), MPI_REAL8, rank_u, itag, comm, MPI_STATUS_IGNORE, ierr)
#endif

  END SUBROUTINE urecv_r8_1d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE urecv_r4(x, tag)
    REAL(4), INTENT(INOUT) :: x
    INTEGER, INTENT(IN), OPTIONAL :: tag
    INTEGER :: itag, ierr

#ifdef PARALLEL3D
    IF (rank_u == MPI_PROC_NULL) RETURN

    itag = 0
    IF (present(tag)) itag = tag

    CALL mpi_recv(x, 1, MPI_REAL4, rank_u, itag, comm, MPI_STATUS_IGNORE, ierr)
#endif

  END SUBROUTINE urecv_r4

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE urecv_r8(x, tag)
    REAL(8), INTENT(INOUT) :: x
    INTEGER, INTENT(IN), OPTIONAL :: tag
    INTEGER :: itag, ierr

#ifdef PARALLEL3D
    IF (rank_u == MPI_PROC_NULL) RETURN

    itag = 0
    IF (present(tag)) itag = tag

    CALL mpi_recv(x, 1, MPI_REAL8, rank_u, itag, comm, MPI_STATUS_IGNORE, ierr)
#endif

  END SUBROUTINE urecv_r8

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE lrecv_r4_3d(x, tag)
    REAL(4), INTENT(INOUT) :: x(:,:,:)
    INTEGER, INTENT(IN), OPTIONAL :: tag
#ifdef F2008
    CONTIGUOUS x
#endif
    INTEGER :: itag, ierr

#ifdef PARALLEL3D
    IF (rank_l == MPI_PROC_NULL) RETURN

    itag = 0
    IF (present(tag)) itag = tag

    CALL mpi_recv(x, size(x), MPI_REAL4, rank_l, itag, comm, MPI_STATUS_IGNORE, ierr)
#endif

  END SUBROUTINE lrecv_r4_3d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE lrecv_r8_3d(x, tag)
    REAL(8), INTENT(INOUT) :: x(:,:,:)
    INTEGER, INTENT(IN), OPTIONAL :: tag
#ifdef F2008
    CONTIGUOUS x
#endif
    INTEGER :: itag, ierr

#ifdef PARALLEL3D
    IF (rank_l == MPI_PROC_NULL) RETURN

    itag = 0
    IF (present(tag)) itag = tag

    CALL mpi_recv(x, size(x), MPI_REAL8, rank_l, itag, comm, MPI_STATUS_IGNORE, ierr)
#endif

  END SUBROUTINE lrecv_r8_3d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE lrecv_r4_2d(x, tag)
    REAL(4), INTENT(INOUT) :: x(:,:)
    INTEGER, INTENT(IN), OPTIONAL :: tag
#ifdef F2008
    CONTIGUOUS x
#endif
    INTEGER :: itag, ierr

#ifdef PARALLEL3D
    IF (rank_l == MPI_PROC_NULL) RETURN

    itag = 0
    IF (present(tag)) itag = tag

    CALL mpi_recv(x, size(x), MPI_REAL4, rank_l, itag, comm, MPI_STATUS_IGNORE, ierr)
#endif

  END SUBROUTINE lrecv_r4_2d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE lrecv_r8_2d(x, tag)
    REAL(8), INTENT(INOUT) :: x(:,:)
    INTEGER, INTENT(IN), OPTIONAL :: tag
#ifdef F2008
    CONTIGUOUS x
#endif
    INTEGER :: itag, ierr

#ifdef PARALLEL3D
    IF (rank_l == MPI_PROC_NULL) RETURN

    itag = 0
    IF (present(tag)) itag = tag

    CALL mpi_recv(x, size(x), MPI_REAL8, rank_l, itag, comm, MPI_STATUS_IGNORE, ierr)
#endif

  END SUBROUTINE lrecv_r8_2d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE lrecv_r4_1d(x, tag)
    REAL(4), INTENT(INOUT) :: x(:)
    INTEGER, INTENT(IN), OPTIONAL :: tag
#ifdef F2008
    CONTIGUOUS x
#endif
    INTEGER :: itag, ierr

#ifdef PARALLEL3D
    IF (rank_l == MPI_PROC_NULL) RETURN

    itag = 0
    IF (present(tag)) itag = tag

    CALL mpi_recv(x, size(x), MPI_REAL4, rank_l, itag, comm, MPI_STATUS_IGNORE, ierr)
#endif

  END SUBROUTINE lrecv_r4_1d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE lrecv_r8_1d(x, tag)
    REAL(8), INTENT(INOUT) :: x(:)
    INTEGER, INTENT(IN), OPTIONAL :: tag
    INTEGER :: itag, ierr

#ifdef PARALLEL3D
    IF (rank_l == MPI_PROC_NULL) RETURN

    itag = 0
    IF (present(tag)) itag = tag

    CALL mpi_recv(x, size(x), MPI_REAL8, rank_l, itag, comm, MPI_STATUS_IGNORE, ierr)
#endif

  END SUBROUTINE lrecv_r8_1d

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE lrecv_r4(x, tag)
    REAL(4), INTENT(INOUT) :: x
    INTEGER, INTENT(IN), OPTIONAL :: tag
    INTEGER :: itag, ierr

#ifdef PARALLEL3D
    IF (rank_l == MPI_PROC_NULL) RETURN

    itag = 0
    IF (present(tag)) itag = tag

    CALL mpi_recv(x, 1, MPI_REAL4, rank_l, itag, comm, MPI_STATUS_IGNORE, ierr)
#endif

  END SUBROUTINE lrecv_r4

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE lrecv_r8(x, tag)
    REAL(8), INTENT(INOUT) :: x
    INTEGER, INTENT(IN), OPTIONAL :: tag
    INTEGER :: itag, ierr

#ifdef PARALLEL3D
    IF (rank_l == MPI_PROC_NULL) RETURN

    itag = 0
    IF (present(tag)) itag = tag

    CALL mpi_recv(x, 1, MPI_REAL8, rank_l, itag, comm, MPI_STATUS_IGNORE, ierr)
#endif

  END SUBROUTINE lrecv_r8

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE barrier
#ifdef PARALLEL_MPI
    INTEGER :: ierr
    CALL mpi_barrier(comm, ierr)
#endif
  END SUBROUTINE barrier

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE serial_begin
#ifdef PARALLEL_MPI
    INTEGER :: dummy, ierr
    dummy = 0
    IF (rank > 0) CALL mpi_recv(dummy, 1, MPI_INTEGER, rank-1, 99, comm, MPI_STATUS_IGNORE, ierr)
#endif
  END SUBROUTINE serial_begin

  SUBROUTINE serial_end
#ifdef PARALLEL_MPI
    INTEGER :: dummy, ierr
    dummy = 0
    IF (rank < npes-1) CALL mpi_send(dummy, 1, MPI_INTEGER, rank+1, 99, comm, ierr)
    CALL mpi_barrier(comm, ierr)
#endif
  END SUBROUTINE serial_end

!-----------------------------------------------------------------------------------------------------------------------

END MODULE parallel
