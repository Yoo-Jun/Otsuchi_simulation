#include "macro.h"

MODULE solver_util
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

  INTEGER, PARAMETER :: max_level = 32
  INTEGER, PARAMETER :: min_dim   = 2

#ifdef MG_REAL4
  INTEGER, PARAMETER :: MG_KIND = 4
#else
  INTEGER, PARAMETER :: MG_KIND = 8
#endif

  INTERFACE fine2d
     MODULE PROCEDURE fine2d_real4
     MODULE PROCEDURE fine2d_real8
     MODULE PROCEDURE fine2d_mix
     MODULE PROCEDURE fine2d_blocked_real4
     MODULE PROCEDURE fine2d_blocked_real8
     MODULE PROCEDURE fine2d_blocked_mix
  END INTERFACE fine2d

  INTERFACE fine3d
     MODULE PROCEDURE fine3d_real4
     MODULE PROCEDURE fine3d_real8
     MODULE PROCEDURE fine3d_mix
     MODULE PROCEDURE fine3d_blocked_real4
     MODULE PROCEDURE fine3d_blocked_real8
     MODULE PROCEDURE fine3d_blocked_mix
  END INTERFACE fine3d

  INTERFACE fine2d3d
     MODULE PROCEDURE fine2d3d_real4
     MODULE PROCEDURE fine2d3d_real8
     MODULE PROCEDURE fine2d3d_mix
     MODULE PROCEDURE fine2d3d_blocked_real4
     MODULE PROCEDURE fine2d3d_blocked_real8
     MODULE PROCEDURE fine2d3d_blocked_mix
  END INTERFACE fine2d3d

  INTERFACE coarse2d
     MODULE PROCEDURE coarse2d_real4
     MODULE PROCEDURE coarse2d_real8
     MODULE PROCEDURE coarse2d_mix
     MODULE PROCEDURE coarse2d_blocked_real4
     MODULE PROCEDURE coarse2d_blocked_real8
     MODULE PROCEDURE coarse2d_blocked_mix
  END INTERFACE coarse2d

  INTERFACE coarse3d
     MODULE PROCEDURE coarse3d_real4
     MODULE PROCEDURE coarse3d_real8
     MODULE PROCEDURE coarse3d_mix
     MODULE PROCEDURE coarse3d_blocked_real4
     MODULE PROCEDURE coarse3d_blocked_real8
     MODULE PROCEDURE coarse3d_blocked_mix
  END INTERFACE coarse3d

  INTERFACE coarse3d2d
     MODULE PROCEDURE coarse3d2d_real4
     MODULE PROCEDURE coarse3d2d_real8
     MODULE PROCEDURE coarse3d2d_mix
     MODULE PROCEDURE coarse3d2d_blocked_real4
     MODULE PROCEDURE coarse3d2d_blocked_real8
     MODULE PROCEDURE coarse3d2d_blocked_mix
  END INTERFACE coarse3d2d

  INTERFACE apply2d1
     MODULE PROCEDURE apply2d1_real4
     MODULE PROCEDURE apply2d1_real8
     MODULE PROCEDURE apply2d1_mix
     MODULE PROCEDURE apply2d1_blocked_real4
     MODULE PROCEDURE apply2d1_blocked_real8
     MODULE PROCEDURE apply2d1_blocked_mix
  END INTERFACE apply2d1

  INTERFACE apply2d2
     MODULE PROCEDURE apply2d2_real4
     MODULE PROCEDURE apply2d2_real8
     MODULE PROCEDURE apply2d2_mix
     MODULE PROCEDURE apply2d2_blocked_real4
     MODULE PROCEDURE apply2d2_blocked_real8
     MODULE PROCEDURE apply2d2_blocked_mix
  END INTERFACE apply2d2

  INTERFACE apply2d3
     MODULE PROCEDURE apply2d3_real4
     MODULE PROCEDURE apply2d3_real8
     MODULE PROCEDURE apply2d3_mix
     MODULE PROCEDURE apply2d3_blocked_real4
     MODULE PROCEDURE apply2d3_blocked_real8
     MODULE PROCEDURE apply2d3_blocked_mix
  END INTERFACE apply2d3

  INTERFACE apply3d1
     MODULE PROCEDURE apply3d1_real4
     MODULE PROCEDURE apply3d1_real8
     MODULE PROCEDURE apply3d1_mix
     MODULE PROCEDURE apply3d1_blocked_real4
     MODULE PROCEDURE apply3d1_blocked_real8
     MODULE PROCEDURE apply3d1_blocked_mix
  END INTERFACE apply3d1

  INTERFACE apply3d2
     MODULE PROCEDURE apply3d2_real4
     MODULE PROCEDURE apply3d2_real8
     MODULE PROCEDURE apply3d2_mix
     MODULE PROCEDURE apply3d2_blocked_real4
     MODULE PROCEDURE apply3d2_blocked_real8
     MODULE PROCEDURE apply3d2_blocked_mix
  END INTERFACE apply3d2

  INTERFACE apply3d3
     MODULE PROCEDURE apply3d3_real4
     MODULE PROCEDURE apply3d3_real8
     MODULE PROCEDURE apply3d3_mix
     MODULE PROCEDURE apply3d3_blocked_real4
     MODULE PROCEDURE apply3d3_blocked_real8
     MODULE PROCEDURE apply3d3_blocked_mix
  END INTERFACE apply3d3

  INTERFACE invmat
     MODULE PROCEDURE invmat_real4
     MODULE PROCEDURE invmat_real8
  END INTERFACE invmat

CONTAINS

!#define FINE_P
!#define COARSE_P
! enable COARSE_P for consistent threading (binary-level matching of the results with diferent thread_num)

#ifdef OMP_CONSISTENT
#define COARSE_P
#endif

SUBROUTINE fine2d_real4(n1, n2, m1, m2, p1, p2, q1, q2, in, out)
  INTEGER, INTENT(IN)  :: n1
  INTEGER, INTENT(IN)  :: n2
  INTEGER, INTENT(IN)  :: m1
  INTEGER, INTENT(IN)  :: m2
  INTEGER, INTENT(IN)  :: p1(0:m1)
  INTEGER, INTENT(IN)  :: p2(0:m2)
  INTEGER, INTENT(IN)  :: q1(1:n1)
  INTEGER, INTENT(IN)  :: q2(1:n2)
  REAL(4), INTENT(IN)    :: in( 0:m1+1, 0:m2+1)
  REAL(4), INTENT(INOUT) :: out(0:n1+1, 0:n2+1)

  INTEGER :: i, j
  INTEGER :: l, m

#ifdef FINE_P
!$OMP PARALLEL DO PRIVATE(l, m)
  DO j=1, m2
  DO m=p2(j-1)+1, p2(j)
     DO i=1, m1
     DO l=p1(i-1)+1, p1(i)
        out(l,m) = out(l,m) + in(i,j)
     END DO
     END DO
  END DO
  END DO
! out(p1(i-1)+1:p1(i),p2(j-1)+1:p2(j)) = out(p1(i-1)+1:p1(i),p2(j-1)+1:p2(j)) + in(i,j)
#else
!$OMP PARALLEL DO PRIVATE(l, m)
  DO j=1, n2
     m = q2(j)
     DO i=1, n1
        l = q1(i)
        out(i,j) = out(i,j) + in(l,m)
     END DO
  END DO
#endif

END SUBROUTINE fine2d_real4

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE fine2d_real8(n1, n2, m1, m2, p1, p2, q1, q2, in, out)
  INTEGER, INTENT(IN)  :: n1
  INTEGER, INTENT(IN)  :: n2
  INTEGER, INTENT(IN)  :: m1
  INTEGER, INTENT(IN)  :: m2
  INTEGER, INTENT(IN)  :: p1(0:m1)
  INTEGER, INTENT(IN)  :: p2(0:m2)
  INTEGER, INTENT(IN)  :: q1(1:n1)
  INTEGER, INTENT(IN)  :: q2(1:n2)
  REAL(8), INTENT(IN)    :: in( 0:m1+1, 0:m2+1)
  REAL(8), INTENT(INOUT) :: out(0:n1+1, 0:n2+1)

  INTEGER :: i, j
  INTEGER :: l, m

#ifdef FINE_P
!$OMP PARALLEL DO PRIVATE(l, m)
  DO j=1, m2
  DO m=p2(j-1)+1, p2(j)
     DO i=1, m1
     DO l=p1(i-1)+1, p1(i)
        out(l,m) = out(l,m) + in(i,j)
     END DO
     END DO
  END DO
  END DO
! out(p1(i-1)+1:p1(i),p2(j-1)+1:p2(j)) = out(p1(i-1)+1:p1(i),p2(j-1)+1:p2(j)) + in(i,j)
#else
!$OMP PARALLEL DO PRIVATE(l, m)
  DO j=1, n2
     m = q2(j)
     DO i=1, n1
        l = q1(i)
        out(i,j) = out(i,j) + in(l,m)
     END DO
  END DO
#endif

END SUBROUTINE fine2d_real8

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE fine2d_mix(n1, n2, m1, m2, p1, p2, q1, q2, in, out)
  INTEGER, INTENT(IN)  :: n1
  INTEGER, INTENT(IN)  :: n2
  INTEGER, INTENT(IN)  :: m1
  INTEGER, INTENT(IN)  :: m2
  INTEGER, INTENT(IN)  :: p1(0:m1)
  INTEGER, INTENT(IN)  :: p2(0:m2)
  INTEGER, INTENT(IN)  :: q1(1:n1)
  INTEGER, INTENT(IN)  :: q2(1:n2)
  REAL(4), INTENT(IN)    :: in( 0:m1+1, 0:m2+1)
  REAL(8), INTENT(INOUT) :: out(0:n1+1, 0:n2+1)

  INTEGER :: i, j
  INTEGER :: l, m

#ifdef FINE_P
!$OMP PARALLEL DO PRIVATE(l, m)
  DO j=1, m2
  DO m=p2(j-1)+1, p2(j)
     DO i=1, m1
     DO l=p1(i-1)+1, p1(i)
        out(l,m) = out(l,m) + in(i,j)
     END DO
     END DO
  END DO
  END DO
! out(p1(i-1)+1:p1(i),p2(j-1)+1:p2(j)) = out(p1(i-1)+1:p1(i),p2(j-1)+1:p2(j)) + in(i,j)
#else
!$OMP PARALLEL DO PRIVATE(l, m)
  DO j=1, n2
     m = q2(j)
     DO i=1, n1
        l = q1(i)
        out(i,j) = out(i,j) + in(l,m)
     END DO
  END DO
#endif

END SUBROUTINE fine2d_mix

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE fine2d_blocked_real4(n1, n2, m1, m2, p1, p2, q1, q2, in, out, jstart, jend)
  INTEGER, INTENT(IN)  :: n1
  INTEGER, INTENT(IN)  :: n2
  INTEGER, INTENT(IN)  :: m1
  INTEGER, INTENT(IN)  :: m2
  INTEGER, INTENT(IN)  :: p1(0:m1)
  INTEGER, INTENT(IN)  :: p2(0:m2)
  INTEGER, INTENT(IN)  :: q1(1:n1)
  INTEGER, INTENT(IN)  :: q2(1:n2)
  REAL(4), INTENT(IN)    :: in( 0:m1+1, 0:m2+1)
  REAL(4), INTENT(INOUT) :: out(0:n1+1, 0:n2+1)
  INTEGER, INTENT(IN)  :: jstart
  INTEGER, INTENT(IN)  :: jend

  INTEGER :: i, j
  INTEGER :: l, m

  DO j=jstart, jend
     m = q2(j)
     DO i=1, n1
        l = q1(i)
        out(i,j) = out(i,j) + in(l,m)
     END DO
  END DO

END SUBROUTINE fine2d_blocked_real4

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE fine2d_blocked_real8(n1, n2, m1, m2, p1, p2, q1, q2, in, out, jstart, jend)
  INTEGER, INTENT(IN)  :: n1
  INTEGER, INTENT(IN)  :: n2
  INTEGER, INTENT(IN)  :: m1
  INTEGER, INTENT(IN)  :: m2
  INTEGER, INTENT(IN)  :: p1(0:m1)
  INTEGER, INTENT(IN)  :: p2(0:m2)
  INTEGER, INTENT(IN)  :: q1(1:n1)
  INTEGER, INTENT(IN)  :: q2(1:n2)
  REAL(8), INTENT(IN)    :: in( 0:m1+1, 0:m2+1)
  REAL(8), INTENT(INOUT) :: out(0:n1+1, 0:n2+1)
  INTEGER, INTENT(IN)  :: jstart
  INTEGER, INTENT(IN)  :: jend

  INTEGER :: i, j
  INTEGER :: l, m

  DO j=jstart, jend
     m = q2(j)
     DO i=1, n1
        l = q1(i)
        out(i,j) = out(i,j) + in(l,m)
     END DO
  END DO

END SUBROUTINE fine2d_blocked_real8

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE fine2d_blocked_mix(n1, n2, m1, m2, p1, p2, q1, q2, in, out, jstart, jend)
  INTEGER, INTENT(IN)  :: n1
  INTEGER, INTENT(IN)  :: n2
  INTEGER, INTENT(IN)  :: m1
  INTEGER, INTENT(IN)  :: m2
  INTEGER, INTENT(IN)  :: p1(0:m1)
  INTEGER, INTENT(IN)  :: p2(0:m2)
  INTEGER, INTENT(IN)  :: q1(1:n1)
  INTEGER, INTENT(IN)  :: q2(1:n2)
  REAL(4), INTENT(IN)    :: in( 0:m1+1, 0:m2+1)
  REAL(8), INTENT(INOUT) :: out(0:n1+1, 0:n2+1)
  INTEGER, INTENT(IN)  :: jstart
  INTEGER, INTENT(IN)  :: jend

  INTEGER :: i, j
  INTEGER :: l, m

  DO j=jstart, jend
     m = q2(j)
     DO i=1, n1
        l = q1(i)
        out(i,j) = out(i,j) + in(l,m)
     END DO
  END DO

END SUBROUTINE fine2d_blocked_mix

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE fine3d_real4(n1, n2, n3, m1, m2, m3, p1, p2, p3, q1, q2, q3, in, out)
  INTEGER, INTENT(IN)  :: n1
  INTEGER, INTENT(IN)  :: n2
  INTEGER, INTENT(IN)  :: n3
  INTEGER, INTENT(IN)  :: m1
  INTEGER, INTENT(IN)  :: m2
  INTEGER, INTENT(IN)  :: m3
  INTEGER, INTENT(IN)  :: p1(0:m1)
  INTEGER, INTENT(IN)  :: p2(0:m2)
  INTEGER, INTENT(IN)  :: p3(0:m3)
  INTEGER, INTENT(IN)  :: q1(1:n1)
  INTEGER, INTENT(IN)  :: q2(1:n2)
  INTEGER, INTENT(IN)  :: q3(1:n3)
  REAL(4), INTENT(IN)    :: in( 0:m1+1, 0:m2+1, 0:m3+1)
  REAL(4), INTENT(INOUT) :: out(0:n1+1, 0:n2+1, 0:n3+1)

  INTEGER :: i, j, k
  INTEGER :: l, m, n

#ifdef FINE_P
!$OMP PARALLEL DO PRIVATE(l, m, n)
  DO k=1, m3
  DO n=p3(k-1)+1, p3(k)
     DO j=1, m2
     DO m=p2(j-1)+1, p2(j)
        DO i=1, m1
        DO l=p1(i-1)+1, p1(i)
           out(l,m,n) = out(l,m,n) + in(i,j,k)
        END DO
        END DO
     END DO
     END DO
  END DO
  END DO
!     out(p1(i-1)+1:p1(i),p2(j-1)+1:p2(j),p3(k-1)+1:p3(k)) = out(p1(i-1)+1:p1(i),p2(j-1)+1:p2(j),p3(k-1)+1:p3(k)) + in(i,j,k)
#else
!$OMP PARALLEL DO PRIVATE(l, m, n)
  DO k=1, n3
     n = q3(k)
     DO j=1, n2
        m = q2(j)
        DO i=1, n1
           l = q1(i)
           out(i,j,k) = out(i,j,k) + in(l,m,n)
        END DO
     END DO
  END DO
#endif

END SUBROUTINE fine3d_real4

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE fine3d_real8(n1, n2, n3, m1, m2, m3, p1, p2, p3, q1, q2, q3, in, out)
  INTEGER, INTENT(IN)  :: n1
  INTEGER, INTENT(IN)  :: n2
  INTEGER, INTENT(IN)  :: n3
  INTEGER, INTENT(IN)  :: m1
  INTEGER, INTENT(IN)  :: m2
  INTEGER, INTENT(IN)  :: m3
  INTEGER, INTENT(IN)  :: p1(0:m1)
  INTEGER, INTENT(IN)  :: p2(0:m2)
  INTEGER, INTENT(IN)  :: p3(0:m3)
  INTEGER, INTENT(IN)  :: q1(1:n1)
  INTEGER, INTENT(IN)  :: q2(1:n2)
  INTEGER, INTENT(IN)  :: q3(1:n3)
  REAL(8), INTENT(IN)    :: in( 0:m1+1, 0:m2+1, 0:m3+1)
  REAL(8), INTENT(INOUT) :: out(0:n1+1, 0:n2+1, 0:n3+1)

  INTEGER :: i, j, k
  INTEGER :: l, m, n

#ifdef FINE_P
!$OMP PARALLEL DO PRIVATE(l, m, n)
  DO k=1, m3
  DO n=p3(k-1)+1, p3(k)
     DO j=1, m2
     DO m=p2(j-1)+1, p2(j)
        DO i=1, m1
        DO l=p1(i-1)+1, p1(i)
           out(l,m,n) = out(l,m,n) + in(i,j,k)
        END DO
        END DO
     END DO
     END DO
  END DO
  END DO
!     out(p1(i-1)+1:p1(i),p2(j-1)+1:p2(j),p3(k-1)+1:p3(k)) = out(p1(i-1)+1:p1(i),p2(j-1)+1:p2(j),p3(k-1)+1:p3(k)) + in(i,j,k)
#else
!$OMP PARALLEL DO PRIVATE(l, m, n)
  DO k=1, n3
     n = q3(k)
     DO j=1, n2
        m = q2(j)
        DO i=1, n1
           l = q1(i)
           out(i,j,k) = out(i,j,k) + in(l,m,n)
        END DO
     END DO
  END DO
#endif

END SUBROUTINE fine3d_real8

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE fine3d_mix(n1, n2, n3, m1, m2, m3, p1, p2, p3, q1, q2, q3, in, out)
  INTEGER, INTENT(IN)  :: n1
  INTEGER, INTENT(IN)  :: n2
  INTEGER, INTENT(IN)  :: n3
  INTEGER, INTENT(IN)  :: m1
  INTEGER, INTENT(IN)  :: m2
  INTEGER, INTENT(IN)  :: m3
  INTEGER, INTENT(IN)  :: p1(0:m1)
  INTEGER, INTENT(IN)  :: p2(0:m2)
  INTEGER, INTENT(IN)  :: p3(0:m3)
  INTEGER, INTENT(IN)  :: q1(1:n1)
  INTEGER, INTENT(IN)  :: q2(1:n2)
  INTEGER, INTENT(IN)  :: q3(1:n3)
  REAL(4), INTENT(IN)    :: in( 0:m1+1, 0:m2+1, 0:m3+1)
  REAL(8), INTENT(INOUT) :: out(0:n1+1, 0:n2+1, 0:n3+1)

  INTEGER :: i, j, k
  INTEGER :: l, m, n

#ifdef FINE_P
!$OMP PARALLEL DO PRIVATE(l, m, n)
  DO k=1, m3
  DO n=p3(k-1)+1, p3(k)
     DO j=1, m2
     DO m=p2(j-1)+1, p2(j)
        DO i=1, m1
        DO l=p1(i-1)+1, p1(i)
           out(l,m,n) = out(l,m,n) + in(i,j,k)
        END DO
        END DO
     END DO
     END DO
  END DO
  END DO
!     out(p1(i-1)+1:p1(i),p2(j-1)+1:p2(j),p3(k-1)+1:p3(k)) = out(p1(i-1)+1:p1(i),p2(j-1)+1:p2(j),p3(k-1)+1:p3(k)) + in(i,j,k)
#else
!$OMP PARALLEL DO PRIVATE(l, m, n)
  DO k=1, n3
     n = q3(k)
     DO j=1, n2
        m = q2(j)
        DO i=1, n1
           l = q1(i)
           out(i,j,k) = out(i,j,k) + in(l,m,n)
        END DO
     END DO
  END DO
#endif

END SUBROUTINE fine3d_mix

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE fine3d_blocked_real4(n1, n2, n3, m1, m2, m3, p1, p2, p3, q1, q2, q3, in, out, kstart, kend)
  INTEGER, INTENT(IN)  :: n1
  INTEGER, INTENT(IN)  :: n2
  INTEGER, INTENT(IN)  :: n3
  INTEGER, INTENT(IN)  :: m1
  INTEGER, INTENT(IN)  :: m2
  INTEGER, INTENT(IN)  :: m3
  INTEGER, INTENT(IN)  :: p1(0:m1)
  INTEGER, INTENT(IN)  :: p2(0:m2)
  INTEGER, INTENT(IN)  :: p3(0:m3)
  INTEGER, INTENT(IN)  :: q1(1:n1)
  INTEGER, INTENT(IN)  :: q2(1:n2)
  INTEGER, INTENT(IN)  :: q3(1:n3)
  REAL(4), INTENT(IN)    :: in( 0:m1+1, 0:m2+1, 0:m3+1)
  REAL(4), INTENT(INOUT) :: out(0:n1+1, 0:n2+1, 0:n3+1)
  INTEGER, INTENT(IN)  :: kstart
  INTEGER, INTENT(IN)  :: kend

  INTEGER :: i, j, k
  INTEGER :: l, m, n

  DO k=kstart, kend
     n = q3(k)
     DO j=1, n2
        m = q2(j)
        DO i=1, n1
           l = q1(i)
           out(i,j,k) = out(i,j,k) + in(l,m,n)
        END DO
     END DO
  END DO

END SUBROUTINE fine3d_blocked_real4

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE fine3d_blocked_real8(n1, n2, n3, m1, m2, m3, p1, p2, p3, q1, q2, q3, in, out, kstart, kend)
  INTEGER, INTENT(IN)  :: n1
  INTEGER, INTENT(IN)  :: n2
  INTEGER, INTENT(IN)  :: n3
  INTEGER, INTENT(IN)  :: m1
  INTEGER, INTENT(IN)  :: m2
  INTEGER, INTENT(IN)  :: m3
  INTEGER, INTENT(IN)  :: p1(0:m1)
  INTEGER, INTENT(IN)  :: p2(0:m2)
  INTEGER, INTENT(IN)  :: p3(0:m3)
  INTEGER, INTENT(IN)  :: q1(1:n1)
  INTEGER, INTENT(IN)  :: q2(1:n2)
  INTEGER, INTENT(IN)  :: q3(1:n3)
  REAL(8), INTENT(IN)    :: in( 0:m1+1, 0:m2+1, 0:m3+1)
  REAL(8), INTENT(INOUT) :: out(0:n1+1, 0:n2+1, 0:n3+1)
  INTEGER, INTENT(IN)  :: kstart
  INTEGER, INTENT(IN)  :: kend

  INTEGER :: i, j, k
  INTEGER :: l, m, n

  DO k=kstart, kend
     n = q3(k)
     DO j=1, n2
        m = q2(j)
        DO i=1, n1
           l = q1(i)
           out(i,j,k) = out(i,j,k) + in(l,m,n)
        END DO
     END DO
  END DO

END SUBROUTINE fine3d_blocked_real8

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE fine3d_blocked_mix(n1, n2, n3, m1, m2, m3, p1, p2, p3, q1, q2, q3, in, out, kstart, kend)
  INTEGER, INTENT(IN)  :: n1
  INTEGER, INTENT(IN)  :: n2
  INTEGER, INTENT(IN)  :: n3
  INTEGER, INTENT(IN)  :: m1
  INTEGER, INTENT(IN)  :: m2
  INTEGER, INTENT(IN)  :: m3
  INTEGER, INTENT(IN)  :: p1(0:m1)
  INTEGER, INTENT(IN)  :: p2(0:m2)
  INTEGER, INTENT(IN)  :: p3(0:m3)
  INTEGER, INTENT(IN)  :: q1(1:n1)
  INTEGER, INTENT(IN)  :: q2(1:n2)
  INTEGER, INTENT(IN)  :: q3(1:n3)
  REAL(4), INTENT(IN)    :: in( 0:m1+1, 0:m2+1, 0:m3+1)
  REAL(8), INTENT(INOUT) :: out(0:n1+1, 0:n2+1, 0:n3+1)
  INTEGER, INTENT(IN)  :: kstart
  INTEGER, INTENT(IN)  :: kend

  INTEGER :: i, j, k
  INTEGER :: l, m, n

  DO k=kstart, kend
     n = q3(k)
     DO j=1, n2
        m = q2(j)
        DO i=1, n1
           l = q1(i)
           out(i,j,k) = out(i,j,k) + in(l,m,n)
        END DO
     END DO
  END DO

END SUBROUTINE fine3d_blocked_mix

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE fine2d3d_real4(n1, n2, n3, m1, m2, p1, p2, q1, q2, in, out)
  INTEGER, INTENT(IN)  :: n1
  INTEGER, INTENT(IN)  :: n2
  INTEGER, INTENT(IN)  :: n3
  INTEGER, INTENT(IN)  :: m1
  INTEGER, INTENT(IN)  :: m2
  INTEGER, INTENT(IN)  :: p1(0:m1)
  INTEGER, INTENT(IN)  :: p2(0:m2)
  INTEGER, INTENT(IN)  :: q1(1:n1)
  INTEGER, INTENT(IN)  :: q2(1:n2)
  REAL(4), INTENT(IN)    :: in( 0:m1+1, 0:m2+1)
  REAL(4), INTENT(INOUT) :: out(0:n1+1, 0:n2+1, 0:n3+1)

  INTEGER :: i, j, k
  INTEGER :: l, m, n

#ifdef FINE_P
!$OMP PARALLEL DO PRIVATE(l, m, n)
  DO n=1, n3
     DO j=1, m2
     DO m=p2(j-1)+1, p2(j)
        DO i=1, m1
        DO l=p1(i-1)+1, p1(i)
           out(l,m,n) = out(l,m,n) + in(i,j)
        END DO
        END DO
     END DO
     END DO
  END DO
!     out(p1(i-1)+1:p1(i),p2(j-1)+1:p2(j),1:n3) = out(p1(i-1)+1:p1(i),p2(j-1)+1:p2(j),1:n3) + in(i,j)
#else
!$OMP PARALLEL DO PRIVATE(l, m, n)
  DO k=1, n3
     DO j=1, n2
        m = q2(j)
        DO i=1, n1
           l = q1(i)
           out(i,j,k) = out(i,j,k) + in(l,m)
        END DO
     END DO
  END DO
#endif

END SUBROUTINE fine2d3d_real4

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE fine2d3d_real8(n1, n2, n3, m1, m2, p1, p2, q1, q2, in, out)
  INTEGER, INTENT(IN)  :: n1
  INTEGER, INTENT(IN)  :: n2
  INTEGER, INTENT(IN)  :: n3
  INTEGER, INTENT(IN)  :: m1
  INTEGER, INTENT(IN)  :: m2
  INTEGER, INTENT(IN)  :: p1(0:m1)
  INTEGER, INTENT(IN)  :: p2(0:m2)
  INTEGER, INTENT(IN)  :: q1(1:n1)
  INTEGER, INTENT(IN)  :: q2(1:n2)
  REAL(8), INTENT(IN)    :: in( 0:m1+1, 0:m2+1)
  REAL(8), INTENT(INOUT) :: out(0:n1+1, 0:n2+1, 0:n3+1)

  INTEGER :: i, j, k
  INTEGER :: l, m, n

#ifdef FINE_P
!$OMP PARALLEL DO PRIVATE(l, m, n)
  DO n=1, n3
     DO j=1, m2
     DO m=p2(j-1)+1, p2(j)
        DO i=1, m1
        DO l=p1(i-1)+1, p1(i)
           out(l,m,n) = out(l,m,n) + in(i,j)
        END DO
        END DO
     END DO
     END DO
  END DO
!     out(p1(i-1)+1:p1(i),p2(j-1)+1:p2(j),1:n3) = out(p1(i-1)+1:p1(i),p2(j-1)+1:p2(j),1:n3) + in(i,j)
#else
!$OMP PARALLEL DO PRIVATE(l, m, n)
  DO k=1, n3
     DO j=1, n2
        m = q2(j)
        DO i=1, n1
           l = q1(i)
           out(i,j,k) = out(i,j,k) + in(l,m)
        END DO
     END DO
  END DO
#endif

END SUBROUTINE fine2d3d_real8

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE fine2d3d_mix(n1, n2, n3, m1, m2, p1, p2, q1, q2, in, out)
  INTEGER, INTENT(IN)  :: n1
  INTEGER, INTENT(IN)  :: n2
  INTEGER, INTENT(IN)  :: n3
  INTEGER, INTENT(IN)  :: m1
  INTEGER, INTENT(IN)  :: m2
  INTEGER, INTENT(IN)  :: p1(0:m1)
  INTEGER, INTENT(IN)  :: p2(0:m2)
  INTEGER, INTENT(IN)  :: q1(1:n1)
  INTEGER, INTENT(IN)  :: q2(1:n2)
  REAL(4), INTENT(IN)    :: in( 0:m1+1, 0:m2+1)
  REAL(8), INTENT(INOUT) :: out(0:n1+1, 0:n2+1, 0:n3+1)

  INTEGER :: i, j, k
  INTEGER :: l, m, n

#ifdef FINE_P
!$OMP PARALLEL DO PRIVATE(l, m, n)
  DO n=1, n3
     DO j=1, m2
     DO m=p2(j-1)+1, p2(j)
        DO i=1, m1
        DO l=p1(i-1)+1, p1(i)
           out(l,m,n) = out(l,m,n) + in(i,j)
        END DO
        END DO
     END DO
     END DO
  END DO
!     out(p1(i-1)+1:p1(i),p2(j-1)+1:p2(j),1:n3) = out(p1(i-1)+1:p1(i),p2(j-1)+1:p2(j),1:n3) + in(i,j)
#else
!$OMP PARALLEL DO PRIVATE(l, m, n)
  DO k=1, n3
     DO j=1, n2
        m = q2(j)
        DO i=1, n1
           l = q1(i)
           out(i,j,k) = out(i,j,k) + in(l,m)
        END DO
     END DO
  END DO
#endif

END SUBROUTINE fine2d3d_mix

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE fine2d3d_blocked_real4(n1, n2, n3, m1, m2, p1, p2, q1, q2, in, out, kstart, kend)
  INTEGER, INTENT(IN)  :: n1
  INTEGER, INTENT(IN)  :: n2
  INTEGER, INTENT(IN)  :: n3
  INTEGER, INTENT(IN)  :: m1
  INTEGER, INTENT(IN)  :: m2
  INTEGER, INTENT(IN)  :: p1(0:m1)
  INTEGER, INTENT(IN)  :: p2(0:m2)
  INTEGER, INTENT(IN)  :: q1(1:n1)
  INTEGER, INTENT(IN)  :: q2(1:n2)
  REAL(4), INTENT(IN)    :: in( 0:m1+1, 0:m2+1)
  REAL(4), INTENT(INOUT) :: out(0:n1+1, 0:n2+1, 0:n3+1)
  INTEGER, INTENT(IN)  :: kstart
  INTEGER, INTENT(IN)  :: kend

  INTEGER :: i, j, k
  INTEGER :: l, m

  DO k=kstart, kend
     DO j=1, n2
        m = q2(j)
        DO i=1, n1
           l = q1(i)
           out(i,j,k) = out(i,j,k) + in(l,m)
        END DO
     END DO
  END DO

END SUBROUTINE fine2d3d_blocked_real4

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE fine2d3d_blocked_real8(n1, n2, n3, m1, m2, p1, p2, q1, q2, in, out, kstart, kend)
  INTEGER, INTENT(IN)  :: n1
  INTEGER, INTENT(IN)  :: n2
  INTEGER, INTENT(IN)  :: n3
  INTEGER, INTENT(IN)  :: m1
  INTEGER, INTENT(IN)  :: m2
  INTEGER, INTENT(IN)  :: p1(0:m1)
  INTEGER, INTENT(IN)  :: p2(0:m2)
  INTEGER, INTENT(IN)  :: q1(1:n1)
  INTEGER, INTENT(IN)  :: q2(1:n2)
  REAL(8), INTENT(IN)    :: in( 0:m1+1, 0:m2+1)
  REAL(8), INTENT(INOUT) :: out(0:n1+1, 0:n2+1, 0:n3+1)
  INTEGER, INTENT(IN)  :: kstart
  INTEGER, INTENT(IN)  :: kend

  INTEGER :: i, j, k
  INTEGER :: l, m

  DO k=kstart, kend
     DO j=1, n2
        m = q2(j)
        DO i=1, n1
           l = q1(i)
           out(i,j,k) = out(i,j,k) + in(l,m)
        END DO
     END DO
  END DO

END SUBROUTINE fine2d3d_blocked_real8

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE fine2d3d_blocked_mix(n1, n2, n3, m1, m2, p1, p2, q1, q2, in, out, kstart, kend)
  INTEGER, INTENT(IN)  :: n1
  INTEGER, INTENT(IN)  :: n2
  INTEGER, INTENT(IN)  :: n3
  INTEGER, INTENT(IN)  :: m1
  INTEGER, INTENT(IN)  :: m2
  INTEGER, INTENT(IN)  :: p1(0:m1)
  INTEGER, INTENT(IN)  :: p2(0:m2)
  INTEGER, INTENT(IN)  :: q1(1:n1)
  INTEGER, INTENT(IN)  :: q2(1:n2)
  REAL(4), INTENT(IN)    :: in( 0:m1+1, 0:m2+1)
  REAL(8), INTENT(INOUT) :: out(0:n1+1, 0:n2+1, 0:n3+1)
  INTEGER, INTENT(IN)  :: kstart
  INTEGER, INTENT(IN)  :: kend

  INTEGER :: i, j, k
  INTEGER :: l, m

  DO k=kstart, kend
     DO j=1, n2
        m = q2(j)
        DO i=1, n1
           l = q1(i)
           out(i,j,k) = out(i,j,k) + in(l,m)
        END DO
     END DO
  END DO

END SUBROUTINE fine2d3d_blocked_mix

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE coarse2d_real4(n1, n2, m1, m2, p1, p2, q1, q2, in, out)
  INTEGER, INTENT(IN)  :: n1
  INTEGER, INTENT(IN)  :: n2
  INTEGER, INTENT(IN)  :: m1
  INTEGER, INTENT(IN)  :: m2
  INTEGER, INTENT(IN)  :: p1(0:m1)
  INTEGER, INTENT(IN)  :: p2(0:m2)
  INTEGER, INTENT(IN)  :: q1(1:n1)
  INTEGER, INTENT(IN)  :: q2(1:n2)
  REAL(4), INTENT(IN)  :: in( 0:n1+1, 0:n2+1)
  REAL(4), INTENT(OUT) :: out(0:m1+1, 0:m2+1)

  INTEGER :: i, j
  INTEGER :: l, m

#ifdef COARSE_P
!$OMP PARALLEL DO
  DO j=1, m2
  DO i=1, m1
     out(i,j) = sum(in(p1(i-1)+1:p1(i), p2(j-1)+1:p2(j)))
  END DO
  END DO
#else
  out(:,:) = 0.0

!POPTION INDEP(out)
!$OMP PARALLEL DO PRIVATE(l, m)
  DO j=1, n2
     m = q2(j)
     DO i=1, n1
        l = q1(i)
        out(l,m)  = out(l,m) + in(i,j)
     END DO
  END DO
#endif

END SUBROUTINE coarse2d_real4

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE coarse2d_real8(n1, n2, m1, m2, p1, p2, q1, q2, in, out)
  INTEGER, INTENT(IN)  :: n1
  INTEGER, INTENT(IN)  :: n2
  INTEGER, INTENT(IN)  :: m1
  INTEGER, INTENT(IN)  :: m2
  INTEGER, INTENT(IN)  :: p1(0:m1)
  INTEGER, INTENT(IN)  :: p2(0:m2)
  INTEGER, INTENT(IN)  :: q1(1:n1)
  INTEGER, INTENT(IN)  :: q2(1:n2)
  REAL(8), INTENT(IN)  :: in( 0:n1+1, 0:n2+1)
  REAL(8), INTENT(OUT) :: out(0:m1+1, 0:m2+1)

  INTEGER :: i, j
  INTEGER :: l, m

#ifdef COARSE_P
!$OMP PARALLEL DO
  DO j=1, m2
  DO i=1, m1
     out(i,j) = sum(in(p1(i-1)+1:p1(i), p2(j-1)+1:p2(j)))
  END DO
  END DO
#else
  out(:,:) = 0.0D0

!POPTION INDEP(out)
!$OMP PARALLEL DO PRIVATE(l, m)
  DO j=1, n2
     m = q2(j)
     DO i=1, n1
        l = q1(i)
        out(l,m)  = out(l,m) + in(i,j)
     END DO
  END DO
#endif

END SUBROUTINE coarse2d_real8

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE coarse2d_mix(n1, n2, m1, m2, p1, p2, q1, q2, in, out)
  INTEGER, INTENT(IN)  :: n1
  INTEGER, INTENT(IN)  :: n2
  INTEGER, INTENT(IN)  :: m1
  INTEGER, INTENT(IN)  :: m2
  INTEGER, INTENT(IN)  :: p1(0:m1)
  INTEGER, INTENT(IN)  :: p2(0:m2)
  INTEGER, INTENT(IN)  :: q1(1:n1)
  INTEGER, INTENT(IN)  :: q2(1:n2)
  REAL(8), INTENT(IN)  :: in( 0:n1+1, 0:n2+1)
  REAL(4), INTENT(OUT) :: out(0:m1+1, 0:m2+1)

  INTEGER :: i, j
  INTEGER :: l, m

#ifdef COARSE_P
!$OMP PARALLEL DO
  DO j=1, m2
  DO i=1, m1
     out(i,j) = real(sum(in(p1(i-1)+1:p1(i), p2(j-1)+1:p2(j))))
  END DO
  END DO
#else
  out(:,:) = 0.0

!POPTION INDEP(out)
!$OMP PARALLEL DO PRIVATE(l, m)
  DO j=1, n2
     m = q2(j)
     DO i=1, n1
        l = q1(i)
        out(l,m)  = out(l,m) + real(in(i,j))
     END DO
  END DO
#endif

END SUBROUTINE coarse2d_mix

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE coarse2d_blocked_real4(n1, n2, m1, m2, p1, p2, q1, q2, in, out, jstart, jend)
  INTEGER, INTENT(IN)  :: n1
  INTEGER, INTENT(IN)  :: n2
  INTEGER, INTENT(IN)  :: m1
  INTEGER, INTENT(IN)  :: m2
  INTEGER, INTENT(IN)  :: p1(0:m1)
  INTEGER, INTENT(IN)  :: p2(0:m2)
  INTEGER, INTENT(IN)  :: q1(1:n1)
  INTEGER, INTENT(IN)  :: q2(1:n2)
  REAL(4), INTENT(IN)  :: in( 0:n1+1, 0:n2+1)
  REAL(4), INTENT(OUT) :: out(0:m1+1, 0:m2+1)
  INTEGER, INTENT(IN)  :: jstart
  INTEGER, INTENT(IN)  :: jend

  INTEGER :: i, j

  DO j=jstart, jend
  DO i=1, m1
     out(i,j) = sum(in(p1(i-1)+1:p1(i), p2(j-1)+1:p2(j)))
  END DO
  END DO

END SUBROUTINE coarse2d_blocked_real4

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE coarse2d_blocked_real8(n1, n2, m1, m2, p1, p2, q1, q2, in, out, jstart, jend)
  INTEGER, INTENT(IN)  :: n1
  INTEGER, INTENT(IN)  :: n2
  INTEGER, INTENT(IN)  :: m1
  INTEGER, INTENT(IN)  :: m2
  INTEGER, INTENT(IN)  :: p1(0:m1)
  INTEGER, INTENT(IN)  :: p2(0:m2)
  INTEGER, INTENT(IN)  :: q1(1:n1)
  INTEGER, INTENT(IN)  :: q2(1:n2)
  REAL(8), INTENT(IN)  :: in( 0:n1+1, 0:n2+1)
  REAL(8), INTENT(OUT) :: out(0:m1+1, 0:m2+1)
  INTEGER, INTENT(IN)  :: jstart
  INTEGER, INTENT(IN)  :: jend

  INTEGER :: i, j

  DO j=jstart, jend
  DO i=1, m1
     out(i,j) = sum(in(p1(i-1)+1:p1(i), p2(j-1)+1:p2(j)))
  END DO
  END DO

END SUBROUTINE coarse2d_blocked_real8

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE coarse2d_blocked_mix(n1, n2, m1, m2, p1, p2, q1, q2, in, out, jstart, jend)
  INTEGER, INTENT(IN)  :: n1
  INTEGER, INTENT(IN)  :: n2
  INTEGER, INTENT(IN)  :: m1
  INTEGER, INTENT(IN)  :: m2
  INTEGER, INTENT(IN)  :: p1(0:m1)
  INTEGER, INTENT(IN)  :: p2(0:m2)
  INTEGER, INTENT(IN)  :: q1(1:n1)
  INTEGER, INTENT(IN)  :: q2(1:n2)
  REAL(8), INTENT(IN)  :: in( 0:n1+1, 0:n2+1)
  REAL(4), INTENT(OUT) :: out(0:m1+1, 0:m2+1)
  INTEGER, INTENT(IN)  :: jstart
  INTEGER, INTENT(IN)  :: jend

  INTEGER :: i, j

  DO j=jstart, jend
  DO i=1, m1
     out(i,j) = real(sum(in(p1(i-1)+1:p1(i), p2(j-1)+1:p2(j))))
  END DO
  END DO

END SUBROUTINE coarse2d_blocked_mix

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE coarse3d_real4(n1, n2, n3, m1, m2, m3, p1, p2, p3, q1, q2, q3, in, out)
  INTEGER, INTENT(IN)  :: n1
  INTEGER, INTENT(IN)  :: n2
  INTEGER, INTENT(IN)  :: n3
  INTEGER, INTENT(IN)  :: m1
  INTEGER, INTENT(IN)  :: m2
  INTEGER, INTENT(IN)  :: m3
  INTEGER, INTENT(IN)  :: p1(0:m1)
  INTEGER, INTENT(IN)  :: p2(0:m2)
  INTEGER, INTENT(IN)  :: p3(0:m3)
  INTEGER, INTENT(IN)  :: q1(1:n1)
  INTEGER, INTENT(IN)  :: q2(1:n2)
  INTEGER, INTENT(IN)  :: q3(1:n3)
  REAL(4), INTENT(IN)  :: in( 0:n1+1, 0:n2+1, 0:n3+1)
  REAL(4), INTENT(OUT) :: out(0:m1+1, 0:m2+1, 0:m3+1)

  INTEGER :: i, j, k
  INTEGER :: l, m, n

#ifdef COARSE_P
!$OMP PARALLEL DO
  DO k=1, m3
  DO j=1, m2
  DO i=1, m1
     out(i,j,k) = sum(in(p1(i-1)+1:p1(i), p2(j-1)+1:p2(j), p3(k-1)+1:p3(k)))
  END DO
  END DO
  END DO
#else
  out(:,:,:) = 0.0

!POPTION INDEP(out)
!$OMP PARALLEL DO PRIVATE(l, m, n)
  DO k=1, n3
     n = q3(k)
     DO j=1, n2
        m = q2(j)
        DO i=1, n1
           l = q1(i)
           out(l,m,n)  = out(l,m,n) + in(i,j,k)
        END DO
     END DO
  END DO
#endif

END SUBROUTINE coarse3d_real4

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE coarse3d_real8(n1, n2, n3, m1, m2, m3, p1, p2, p3, q1, q2, q3, in, out)
  INTEGER, INTENT(IN)  :: n1
  INTEGER, INTENT(IN)  :: n2
  INTEGER, INTENT(IN)  :: n3
  INTEGER, INTENT(IN)  :: m1
  INTEGER, INTENT(IN)  :: m2
  INTEGER, INTENT(IN)  :: m3
  INTEGER, INTENT(IN)  :: p1(0:m1)
  INTEGER, INTENT(IN)  :: p2(0:m2)
  INTEGER, INTENT(IN)  :: p3(0:m3)
  INTEGER, INTENT(IN)  :: q1(1:n1)
  INTEGER, INTENT(IN)  :: q2(1:n2)
  INTEGER, INTENT(IN)  :: q3(1:n3)
  REAL(8), INTENT(IN)  :: in( 0:n1+1, 0:n2+1, 0:n3+1)
  REAL(8), INTENT(OUT) :: out(0:m1+1, 0:m2+1, 0:m3+1)

  INTEGER :: i, j, k
  INTEGER :: l, m, n

#ifdef COARSE_P
!$OMP PARALLEL DO
  DO k=1, m3
  DO j=1, m2
  DO i=1, m1
     out(i,j,k) = sum(in(p1(i-1)+1:p1(i), p2(j-1)+1:p2(j), p3(k-1)+1:p3(k)))
  END DO
  END DO
  END DO
#else
  out(:,:,:) = 0.0D0

!POPTION INDEP(out)
!$OMP PARALLEL DO PRIVATE(l, m, n)
  DO k=1, n3
     n = q3(k)
     DO j=1, n2
        m = q2(j)
        DO i=1, n1
           l = q1(i)
           out(l,m,n)  = out(l,m,n) + in(i,j,k)
        END DO
     END DO
  END DO
#endif

END SUBROUTINE coarse3d_real8

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE coarse3d_mix(n1, n2, n3, m1, m2, m3, p1, p2, p3, q1, q2, q3, in, out)
  INTEGER, INTENT(IN)  :: n1
  INTEGER, INTENT(IN)  :: n2
  INTEGER, INTENT(IN)  :: n3
  INTEGER, INTENT(IN)  :: m1
  INTEGER, INTENT(IN)  :: m2
  INTEGER, INTENT(IN)  :: m3
  INTEGER, INTENT(IN)  :: p1(0:m1)
  INTEGER, INTENT(IN)  :: p2(0:m2)
  INTEGER, INTENT(IN)  :: p3(0:m3)
  INTEGER, INTENT(IN)  :: q1(1:n1)
  INTEGER, INTENT(IN)  :: q2(1:n2)
  INTEGER, INTENT(IN)  :: q3(1:n3)
  REAL(8), INTENT(IN)  :: in( 0:n1+1, 0:n2+1, 0:n3+1)
  REAL(4), INTENT(OUT) :: out(0:m1+1, 0:m2+1, 0:m3+1)

  INTEGER :: i, j, k
  INTEGER :: l, m, n

#ifdef COARSE_P
!$OMP PARALLEL DO
  DO k=1, m3
  DO j=1, m2
  DO i=1, m1
     out(i,j,k) = real(sum(in(p1(i-1)+1:p1(i), p2(j-1)+1:p2(j), p3(k-1)+1:p3(k))))
  END DO
  END DO
  END DO
#else
  out(:,:,:) = 0.0

!POPTION INDEP(out)
!$OMP PARALLEL DO PRIVATE(l, m, n)
  DO k=1, n3
     n = q3(k)
     DO j=1, n2
        m = q2(j)
        DO i=1, n1
           l = q1(i)
           out(l,m,n)  = out(l,m,n) + real(in(i,j,k))
        END DO
     END DO
  END DO
#endif

END SUBROUTINE coarse3d_mix

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE coarse3d_blocked_real4(n1, n2, n3, m1, m2, m3, p1, p2, p3, q1, q2, q3, in, out, kstart, kend)
  INTEGER, INTENT(IN)  :: n1
  INTEGER, INTENT(IN)  :: n2
  INTEGER, INTENT(IN)  :: n3
  INTEGER, INTENT(IN)  :: m1
  INTEGER, INTENT(IN)  :: m2
  INTEGER, INTENT(IN)  :: m3
  INTEGER, INTENT(IN)  :: p1(0:m1)
  INTEGER, INTENT(IN)  :: p2(0:m2)
  INTEGER, INTENT(IN)  :: p3(0:m3)
  INTEGER, INTENT(IN)  :: q1(1:n1)
  INTEGER, INTENT(IN)  :: q2(1:n2)
  INTEGER, INTENT(IN)  :: q3(1:n3)
  REAL(4), INTENT(IN)  :: in( 0:n1+1, 0:n2+1, 0:n3+1)
  REAL(4), INTENT(OUT) :: out(0:m1+1, 0:m2+1, 0:m3+1)
  INTEGER, INTENT(IN)  :: kstart
  INTEGER, INTENT(IN)  :: kend

  INTEGER :: i, j, k

  DO k=kstart, kend
     DO j=1, m2
     DO i=1, m1
        out(i,j,k) = sum(in(p1(i-1)+1:p1(i), p2(j-1)+1:p2(j), p3(k-1)+1:p3(k)))
     END DO
     END DO
  END DO

END SUBROUTINE coarse3d_blocked_real4

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE coarse3d_blocked_real8(n1, n2, n3, m1, m2, m3, p1, p2, p3, q1, q2, q3, in, out, kstart, kend)
  INTEGER, INTENT(IN)  :: n1
  INTEGER, INTENT(IN)  :: n2
  INTEGER, INTENT(IN)  :: n3
  INTEGER, INTENT(IN)  :: m1
  INTEGER, INTENT(IN)  :: m2
  INTEGER, INTENT(IN)  :: m3
  INTEGER, INTENT(IN)  :: p1(0:m1)
  INTEGER, INTENT(IN)  :: p2(0:m2)
  INTEGER, INTENT(IN)  :: p3(0:m3)
  INTEGER, INTENT(IN)  :: q1(1:n1)
  INTEGER, INTENT(IN)  :: q2(1:n2)
  INTEGER, INTENT(IN)  :: q3(1:n3)
  REAL(8), INTENT(IN)  :: in( 0:n1+1, 0:n2+1, 0:n3+1)
  REAL(8), INTENT(OUT) :: out(0:m1+1, 0:m2+1, 0:m3+1)
  INTEGER, INTENT(IN)  :: kstart
  INTEGER, INTENT(IN)  :: kend

  INTEGER :: i, j, k

  DO k=kstart, kend
     DO j=1, m2
     DO i=1, m1
        out(i,j,k) = sum(in(p1(i-1)+1:p1(i), p2(j-1)+1:p2(j), p3(k-1)+1:p3(k)))
     END DO
     END DO
  END DO

END SUBROUTINE coarse3d_blocked_real8

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE coarse3d_blocked_mix(n1, n2, n3, m1, m2, m3, p1, p2, p3, q1, q2, q3, in, out, kstart, kend)
  INTEGER, INTENT(IN)  :: n1
  INTEGER, INTENT(IN)  :: n2
  INTEGER, INTENT(IN)  :: n3
  INTEGER, INTENT(IN)  :: m1
  INTEGER, INTENT(IN)  :: m2
  INTEGER, INTENT(IN)  :: m3
  INTEGER, INTENT(IN)  :: p1(0:m1)
  INTEGER, INTENT(IN)  :: p2(0:m2)
  INTEGER, INTENT(IN)  :: p3(0:m3)
  INTEGER, INTENT(IN)  :: q1(1:n1)
  INTEGER, INTENT(IN)  :: q2(1:n2)
  INTEGER, INTENT(IN)  :: q3(1:n3)
  REAL(8), INTENT(IN)  :: in( 0:n1+1, 0:n2+1, 0:n3+1)
  REAL(4), INTENT(OUT) :: out(0:m1+1, 0:m2+1, 0:m3+1)
  INTEGER, INTENT(IN)  :: kstart
  INTEGER, INTENT(IN)  :: kend

  INTEGER :: i, j, k

  DO k=kstart, kend
     DO j=1, m2
     DO i=1, m1
        out(i,j,k) = real(sum(in(p1(i-1)+1:p1(i), p2(j-1)+1:p2(j), p3(k-1)+1:p3(k))))
     END DO
     END DO
  END DO
END SUBROUTINE coarse3d_blocked_mix

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE coarse3d2d_real4(n1, n2, n3, m1, m2, p1, p2, q1, q2, in, out)
  INTEGER, INTENT(IN)  :: n1
  INTEGER, INTENT(IN)  :: n2
  INTEGER, INTENT(IN)  :: n3
  INTEGER, INTENT(IN)  :: m1
  INTEGER, INTENT(IN)  :: m2
  INTEGER, INTENT(IN)  :: p1(0:m1)
  INTEGER, INTENT(IN)  :: p2(0:m2)
  INTEGER, INTENT(IN)  :: q1(1:n1)
  INTEGER, INTENT(IN)  :: q2(1:n2)
  REAL(4), INTENT(IN)  :: in( 0:n1+1, 0:n2+1, 0:n3+1)
  REAL(4), INTENT(OUT) :: out(0:m1+1, 0:m2+1)

  INTEGER :: i, j, k
  INTEGER :: l, m

#ifdef COARSE_P
!$OMP PARALLEL DO
  DO j=1, m2
  DO i=1, m1
     out(i,j) = sum(in(p1(i-1)+1:p1(i), p2(j-1)+1:p2(j), 1:n3))
  END DO
  END DO
#else
  out(:,:) = 0.0

!POPTION INDEP(out)
!$OMP PARALLEL DO PRIVATE(l, m)
  DO k=1, n3
     DO j=1, n2
        m = q2(j)
        DO i=1, n1
           l = q1(i)
           out(l,m)  = out(l,m) + in(i,j,k)
        END DO
     END DO
  END DO
#endif

END SUBROUTINE coarse3d2d_real4

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE coarse3d2d_real8(n1, n2, n3, m1, m2, p1, p2, q1, q2, in, out)
  INTEGER, INTENT(IN)  :: n1
  INTEGER, INTENT(IN)  :: n2
  INTEGER, INTENT(IN)  :: n3
  INTEGER, INTENT(IN)  :: m1
  INTEGER, INTENT(IN)  :: m2
  INTEGER, INTENT(IN)  :: p1(0:m1)
  INTEGER, INTENT(IN)  :: p2(0:m2)
  INTEGER, INTENT(IN)  :: q1(1:n1)
  INTEGER, INTENT(IN)  :: q2(1:n2)
  REAL(8), INTENT(IN)  :: in( 0:n1+1, 0:n2+1, 0:n3+1)
  REAL(8), INTENT(OUT) :: out(0:m1+1, 0:m2+1)

  INTEGER :: i, j, k
  INTEGER :: l, m

#ifdef COARSE_P
!$OMP PARALLEL DO
  DO j=1, m2
  DO i=1, m1
     out(i,j) = sum(in(p1(i-1)+1:p1(i), p2(j-1)+1:p2(j), 1:n3))
  END DO
  END DO
#else
  out(:,:) = 0.0D0

!POPTION INDEP(out)
!$OMP PARALLEL DO PRIVATE(l, m)
  DO k=1, n3
     DO j=1, n2
        m = q2(j)
        DO i=1, n1
           l = q1(i)
           out(l,m)  = out(l,m) + in(i,j,k)
        END DO
     END DO
  END DO
#endif

END SUBROUTINE coarse3d2d_real8

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE coarse3d2d_mix(n1, n2, n3, m1, m2, p1, p2, q1, q2, in, out)
  INTEGER, INTENT(IN)  :: n1
  INTEGER, INTENT(IN)  :: n2
  INTEGER, INTENT(IN)  :: n3
  INTEGER, INTENT(IN)  :: m1
  INTEGER, INTENT(IN)  :: m2
  INTEGER, INTENT(IN)  :: p1(0:m1)
  INTEGER, INTENT(IN)  :: p2(0:m2)
  INTEGER, INTENT(IN)  :: q1(1:n1)
  INTEGER, INTENT(IN)  :: q2(1:n2)
  REAL(8), INTENT(IN)  :: in( 0:n1+1, 0:n2+1, 0:n3+1)
  REAL(4), INTENT(OUT) :: out(0:m1+1, 0:m2+1)

  INTEGER :: i, j, k
  INTEGER :: l, m

#ifdef COARSE_P
!$OMP PARALLEL DO
  DO j=1, m2
  DO i=1, m1
     out(i,j) = real(sum(in(p1(i-1)+1:p1(i), p2(j-1)+1:p2(j), 1:n3)))
  END DO
  END DO
#else
  out(:,:) = 0.0

!POPTION INDEP(out)
!$OMP PARALLEL DO PRIVATE(l, m)
  DO k=1, n3
     DO j=1, n2
        m = q2(j)
        DO i=1, n1
           l = q1(i)
           out(l,m)  = out(l,m) + real(in(i,j,k))
        END DO
     END DO
  END DO
#endif

END SUBROUTINE coarse3d2d_mix

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE coarse3d2d_blocked_real4(n1, n2, n3, m1, m2, p1, p2, q1, q2, in, out, jstart, jend)
  INTEGER, INTENT(IN)  :: n1
  INTEGER, INTENT(IN)  :: n2
  INTEGER, INTENT(IN)  :: n3
  INTEGER, INTENT(IN)  :: m1
  INTEGER, INTENT(IN)  :: m2
  INTEGER, INTENT(IN)  :: p1(0:m1)
  INTEGER, INTENT(IN)  :: p2(0:m2)
  INTEGER, INTENT(IN)  :: q1(1:n1)
  INTEGER, INTENT(IN)  :: q2(1:n2)
  REAL(4), INTENT(IN)  :: in( 0:n1+1, 0:n2+1, 0:n3+1)
  REAL(4), INTENT(OUT) :: out(0:m1+1, 0:m2+1)
  INTEGER, INTENT(IN)  :: jstart
  INTEGER, INTENT(IN)  :: jend

  INTEGER :: i, j

  DO j=jstart, jend
  DO i=1, m1
     out(i,j) = sum(in(p1(i-1)+1:p1(i), p2(j-1)+1:p2(j), 1:n3))
  END DO
  END DO

END SUBROUTINE coarse3d2d_blocked_real4

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE coarse3d2d_blocked_real8(n1, n2, n3, m1, m2, p1, p2, q1, q2, in, out, jstart, jend)
  INTEGER, INTENT(IN)  :: n1
  INTEGER, INTENT(IN)  :: n2
  INTEGER, INTENT(IN)  :: n3
  INTEGER, INTENT(IN)  :: m1
  INTEGER, INTENT(IN)  :: m2
  INTEGER, INTENT(IN)  :: p1(0:m1)
  INTEGER, INTENT(IN)  :: p2(0:m2)
  INTEGER, INTENT(IN)  :: q1(1:n1)
  INTEGER, INTENT(IN)  :: q2(1:n2)
  REAL(8), INTENT(IN)  :: in( 0:n1+1, 0:n2+1, 0:n3+1)
  REAL(8), INTENT(OUT) :: out(0:m1+1, 0:m2+1)
  INTEGER, INTENT(IN)  :: jstart
  INTEGER, INTENT(IN)  :: jend

  INTEGER :: i, j

  DO j=jstart, jend
  DO i=1, m1
     out(i,j) = sum(in(p1(i-1)+1:p1(i), p2(j-1)+1:p2(j), 1:n3))
  END DO
  END DO

END SUBROUTINE coarse3d2d_blocked_real8

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE coarse3d2d_blocked_mix(n1, n2, n3, m1, m2, p1, p2, q1, q2, in, out, jstart, jend)
  INTEGER, INTENT(IN)  :: n1
  INTEGER, INTENT(IN)  :: n2
  INTEGER, INTENT(IN)  :: n3
  INTEGER, INTENT(IN)  :: m1
  INTEGER, INTENT(IN)  :: m2
  INTEGER, INTENT(IN)  :: p1(0:m1)
  INTEGER, INTENT(IN)  :: p2(0:m2)
  INTEGER, INTENT(IN)  :: q1(1:n1)
  INTEGER, INTENT(IN)  :: q2(1:n2)
  REAL(8), INTENT(IN)  :: in( 0:n1+1, 0:n2+1, 0:n3+1)
  REAL(4), INTENT(OUT) :: out(0:m1+1, 0:m2+1)
  INTEGER, INTENT(IN)  :: jstart
  INTEGER, INTENT(IN)  :: jend

  INTEGER :: i, j

  DO j=jstart, jend
  DO i=1, m1
     out(i,j) = real(sum(in(p1(i-1)+1:p1(i), p2(j-1)+1:p2(j), 1:n3)))
  END DO
  END DO

END SUBROUTINE coarse3d2d_blocked_mix

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE apply2d1_real4(n1, n2, a, x, out)
  INTEGER, INTENT(IN)  :: n1
  INTEGER, INTENT(IN)  :: n2
  REAL(4), INTENT(IN)  :: a(-2:2, 1:n1, 1:n2)
  REAL(4), INTENT(IN)  :: x(  0:n1+1, 0:n2+1)
  REAL(4), INTENT(OUT) :: out(0:n1+1, 0:n2+1)

  INTEGER :: i, j

!$OMP PARALLEL DO
  DO j=1, n2
  DO i=1, n1
     out(i,j) = a(-2,i,j) * x(i,  j-1) &
              + a(-1,i,j) * x(i-1,j)   &
              + a( 0,i,j) * x(i,  j)   &
              + a( 1,i,j) * x(i+1,j)   &
              + a( 2,i,j) * x(i,  j+1)
  END DO
  END DO

END SUBROUTINE apply2d1_real4

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE apply2d1_real8(n1, n2, a, x, out)
  INTEGER, INTENT(IN)  :: n1
  INTEGER, INTENT(IN)  :: n2
  REAL(8), INTENT(IN)  :: a(-2:2, 1:n1, 1:n2)
  REAL(8), INTENT(IN)  :: x(  0:n1+1, 0:n2+1)
  REAL(8), INTENT(OUT) :: out(0:n1+1, 0:n2+1)

  INTEGER :: i, j

!$OMP PARALLEL DO
  DO j=1, n2
  DO i=1, n1
     out(i,j) = a(-2,i,j) * x(i,  j-1) &
              + a(-1,i,j) * x(i-1,j)   &
              + a( 0,i,j) * x(i,  j)   &
              + a( 1,i,j) * x(i+1,j)   &
              + a( 2,i,j) * x(i,  j+1)
  END DO
  END DO

END SUBROUTINE apply2d1_real8

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE apply2d1_mix(n1, n2, a, x, out)
  INTEGER, INTENT(IN)  :: n1
  INTEGER, INTENT(IN)  :: n2
  REAL(4), INTENT(IN)  :: a(-2:2, 1:n1, 1:n2)
  REAL(8), INTENT(IN)  :: x(  0:n1+1, 0:n2+1)
  REAL(8), INTENT(OUT) :: out(0:n1+1, 0:n2+1)

  INTEGER :: i, j

!$OMP PARALLEL DO
  DO j=1, n2
  DO i=1, n1
     out(i,j) = a(-2,i,j) * x(i,  j-1) &
              + a(-1,i,j) * x(i-1,j)   &
              + a( 0,i,j) * x(i,  j)   &
              + a( 1,i,j) * x(i+1,j)   &
              + a( 2,i,j) * x(i,  j+1)
  END DO
  END DO

END SUBROUTINE apply2d1_mix

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE apply2d1_blocked_real4(n1, n2, a, x, out, jstart, jend)
  INTEGER, INTENT(IN)  :: n1
  INTEGER, INTENT(IN)  :: n2
  REAL(4), INTENT(IN)  :: a(-2:2, 1:n1, 1:n2)
  REAL(4), INTENT(IN)  :: x(  0:n1+1, 0:n2+1)
  REAL(4), INTENT(OUT) :: out(0:n1+1, 0:n2+1)
  INTEGER, INTENT(IN)  :: jstart
  INTEGER, INTENT(IN)  :: jend

  INTEGER :: i, j

  DO j=jstart, jend
  DO i=1, n1
     out(i,j) = a(-2,i,j) * x(i,  j-1) &
              + a(-1,i,j) * x(i-1,j)   &
              + a( 0,i,j) * x(i,  j)   &
              + a( 1,i,j) * x(i+1,j)   &
              + a( 2,i,j) * x(i,  j+1)
  END DO
  END DO

END SUBROUTINE apply2d1_blocked_real4

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE apply2d1_blocked_real8(n1, n2, a, x, out, jstart, jend)
  INTEGER, INTENT(IN)  :: n1
  INTEGER, INTENT(IN)  :: n2
  REAL(8), INTENT(IN)  :: a(-2:2, 1:n1, 1:n2)
  REAL(8), INTENT(IN)  :: x(  0:n1+1, 0:n2+1)
  REAL(8), INTENT(OUT) :: out(0:n1+1, 0:n2+1)
  INTEGER, INTENT(IN)  :: jstart
  INTEGER, INTENT(IN)  :: jend

  INTEGER :: i, j

  DO j=jstart, jend
  DO i=1, n1
     out(i,j) = a(-2,i,j) * x(i,  j-1) &
              + a(-1,i,j) * x(i-1,j)   &
              + a( 0,i,j) * x(i,  j)   &
              + a( 1,i,j) * x(i+1,j)   &
              + a( 2,i,j) * x(i,  j+1)
  END DO
  END DO

END SUBROUTINE apply2d1_blocked_real8

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE apply2d1_blocked_mix(n1, n2, a, x, out, jstart, jend)
  INTEGER, INTENT(IN)  :: n1
  INTEGER, INTENT(IN)  :: n2
  REAL(4), INTENT(IN)  :: a(-2:2, 1:n1, 1:n2)
  REAL(8), INTENT(IN)  :: x(  0:n1+1, 0:n2+1)
  REAL(8), INTENT(OUT) :: out(0:n1+1, 0:n2+1)
  INTEGER, INTENT(IN)  :: jstart
  INTEGER, INTENT(IN)  :: jend

  INTEGER :: i, j

  DO j=jstart, jend
  DO i=1, n1
     out(i,j) = a(-2,i,j) * x(i,  j-1) &
               + a(-1,i,j) * x(i-1,j)   &
               + a( 0,i,j) * x(i,  j)   &
               + a( 1,i,j) * x(i+1,j)   &
               + a( 2,i,j) * x(i,  j+1)
  END DO
  END DO

END SUBROUTINE apply2d1_blocked_mix

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE apply2d2_real4(n1, n2, a, x, y, out)
  INTEGER, INTENT(IN)  :: n1
  INTEGER, INTENT(IN)  :: n2
  REAL(4), INTENT(IN)  :: a(-2:2, 1:n1, 1:n2)
  REAL(4), INTENT(IN)  :: x(  0:n1+1, 0:n2+1)
  REAL(4), INTENT(IN)  :: y(  0:n1+1, 0:n2+1)
  REAL(4), INTENT(OUT) :: out(0:n1+1, 0:n2+1)

  INTEGER :: i, j

!$OMP PARALLEL DO
  DO j=1, n2
  DO i=1, n1
     out(i,j) = x(i,j)             &
          - a(-2,i,j) * y(i,  j-1) &
          - a(-1,i,j) * y(i-1,j)   &
          - a( 0,i,j) * y(i,  j)   &
          - a( 1,i,j) * y(i+1,j)   &
          - a( 2,i,j) * y(i,  j+1)
  END DO
  END DO

END SUBROUTINE apply2d2_real4

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE apply2d2_real8(n1, n2, a, x, y, out)
  INTEGER, INTENT(IN)  :: n1
  INTEGER, INTENT(IN)  :: n2
  REAL(8), INTENT(IN)  :: a(-2:2, 1:n1, 1:n2)
  REAL(8), INTENT(IN)  :: x(  0:n1+1, 0:n2+1)
  REAL(8), INTENT(IN)  :: y(  0:n1+1, 0:n2+1)
  REAL(8), INTENT(OUT) :: out(0:n1+1, 0:n2+1)

  INTEGER :: i, j

!$OMP PARALLEL DO
  DO j=1, n2
  DO i=1, n1
     out(i,j) = x(i,j)             &
          - a(-2,i,j) * y(i,  j-1) &
          - a(-1,i,j) * y(i-1,j)   &
          - a( 0,i,j) * y(i,  j)   &
          - a( 1,i,j) * y(i+1,j)   &
          - a( 2,i,j) * y(i,  j+1)
  END DO
  END DO

END SUBROUTINE apply2d2_real8

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE apply2d2_mix(n1, n2, a, x, y, out)
  INTEGER, INTENT(IN)  :: n1
  INTEGER, INTENT(IN)  :: n2
  REAL(4), INTENT(IN)  :: a(-2:2, 1:n1, 1:n2)
  REAL(8), INTENT(IN)  :: x(  0:n1+1, 0:n2+1)
  REAL(8), INTENT(IN)  :: y(  0:n1+1, 0:n2+1)
  REAL(8), INTENT(OUT) :: out(0:n1+1, 0:n2+1)

  INTEGER :: i, j

!$OMP PARALLEL DO
  DO j=1, n2
  DO i=1, n1
     out(i,j) = x(i,j)             &
          - a(-2,i,j) * y(i,  j-1) &
          - a(-1,i,j) * y(i-1,j)   &
          - a( 0,i,j) * y(i,  j)   &
          - a( 1,i,j) * y(i+1,j)   &
          - a( 2,i,j) * y(i,  j+1)
  END DO
  END DO

END SUBROUTINE apply2d2_mix

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE apply2d2_blocked_real4(n1, n2, a, x, y, out, jstart, jend)
  INTEGER, INTENT(IN)  :: n1
  INTEGER, INTENT(IN)  :: n2
  REAL(4), INTENT(IN)  :: a(-2:2, 1:n1, 1:n2)
  REAL(4), INTENT(IN)  :: x(  0:n1+1, 0:n2+1)
  REAL(4), INTENT(IN)  :: y(  0:n1+1, 0:n2+1)
  REAL(4), INTENT(OUT) :: out(0:n1+1, 0:n2+1)
  INTEGER, INTENT(IN)  :: jstart
  INTEGER, INTENT(IN)  :: jend

  INTEGER :: i, j

  DO j=jstart, jend
  DO i=1, n1
     out(i,j) = x(i,j)             &
          - a(-2,i,j) * y(i,  j-1) &
          - a(-1,i,j) * y(i-1,j)   &
          - a( 0,i,j) * y(i,  j)   &
          - a( 1,i,j) * y(i+1,j)   &
          - a( 2,i,j) * y(i,  j+1)
  END DO
  END DO

END SUBROUTINE apply2d2_blocked_real4

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE apply2d2_blocked_real8(n1, n2, a, x, y, out, jstart, jend)
  INTEGER, INTENT(IN)  :: n1
  INTEGER, INTENT(IN)  :: n2
  REAL(8), INTENT(IN)  :: a(-2:2, 1:n1, 1:n2)
  REAL(8), INTENT(IN)  :: x(  0:n1+1, 0:n2+1)
  REAL(8), INTENT(IN)  :: y(  0:n1+1, 0:n2+1)
  REAL(8), INTENT(OUT) :: out(0:n1+1, 0:n2+1)
  INTEGER, INTENT(IN)  :: jstart
  INTEGER, INTENT(IN)  :: jend

  INTEGER :: i, j

  DO j=jstart, jend
  DO i=1, n1
     out(i,j) = x(i,j)             &
          - a(-2,i,j) * y(i,  j-1) &
          - a(-1,i,j) * y(i-1,j)   &
          - a( 0,i,j) * y(i,  j)   &
          - a( 1,i,j) * y(i+1,j)   &
          - a( 2,i,j) * y(i,  j+1)
  END DO
  END DO

END SUBROUTINE apply2d2_blocked_real8

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE apply2d2_blocked_mix(n1, n2, a, x, y, out, jstart, jend)
  INTEGER, INTENT(IN)  :: n1
  INTEGER, INTENT(IN)  :: n2
  REAL(4), INTENT(IN)  :: a(-2:2, 1:n1, 1:n2)
  REAL(8), INTENT(IN)  :: x(  0:n1+1, 0:n2+1)
  REAL(8), INTENT(IN)  :: y(  0:n1+1, 0:n2+1)
  REAL(8), INTENT(OUT) :: out(0:n1+1, 0:n2+1)
  INTEGER, INTENT(IN)  :: jstart
  INTEGER, INTENT(IN)  :: jend

  INTEGER :: i, j

  DO j=jstart, jend
  DO i=1, n1
     out(i,j) = x(i,j)             &
          - a(-2,i,j) * y(i,  j-1) &
          - a(-1,i,j) * y(i-1,j)   &
          - a( 0,i,j) * y(i,  j)   &
          - a( 1,i,j) * y(i+1,j)   &
          - a( 2,i,j) * y(i,  j+1)
  END DO
  END DO

END SUBROUTINE apply2d2_blocked_mix

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE apply2d3_real4(n1, n2, a, x, out)
  INTEGER, INTENT(IN)    :: n1
  INTEGER, INTENT(IN)    :: n2
  REAL(4), INTENT(IN)    :: a(-2:2, 1:n1, 1:n2)
  REAL(4), INTENT(IN)    :: x(  0:n1+1, 0:n2+1)
  REAL(4), INTENT(INOUT) :: out(0:n1+1, 0:n2+1)

  INTEGER :: i, j

!$OMP PARALLEL DO
  DO j=1, n2
  DO i=1, n1
     out(i,j) = out(i,j)           &
          + a(-2,i,j) * x(i,  j-1) &
          + a(-1,i,j) * x(i-1,j)   &
          + a( 0,i,j) * x(i,  j)   &
          + a( 1,i,j) * x(i+1,j)   &
          + a( 2,i,j) * x(i,  j+1)
  END DO
  END DO

END SUBROUTINE apply2d3_real4

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE apply2d3_real8(n1, n2, a, x, out)
  INTEGER, INTENT(IN)    :: n1
  INTEGER, INTENT(IN)    :: n2
  REAL(8), INTENT(IN)    :: a(-2:2, 1:n1, 1:n2)
  REAL(8), INTENT(IN)    :: x(  0:n1+1, 0:n2+1)
  REAL(8), INTENT(INOUT) :: out(0:n1+1, 0:n2+1)

  INTEGER :: i, j

!$OMP PARALLEL DO
  DO j=1, n2
  DO i=1, n1
     out(i,j) = out(i,j)           &
          + a(-2,i,j) * x(i,  j-1) &
          + a(-1,i,j) * x(i-1,j)   &
          + a( 0,i,j) * x(i,  j)   &
          + a( 1,i,j) * x(i+1,j)   &
          + a( 2,i,j) * x(i,  j+1)
  END DO
  END DO

END SUBROUTINE apply2d3_real8

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE apply2d3_mix(n1, n2, a, x, out)
  INTEGER, INTENT(IN)    :: n1
  INTEGER, INTENT(IN)    :: n2
  REAL(4), INTENT(IN)    :: a(-2:2, 1:n1, 1:n2)
  REAL(8), INTENT(IN)    :: x(  0:n1+1, 0:n2+1)
  REAL(8), INTENT(INOUT) :: out(0:n1+1, 0:n2+1)

  INTEGER :: i, j

!$OMP PARALLEL DO
  DO j=1, n2
  DO i=1, n1
     out(i,j) = out(i,j)           &
          + a(-2,i,j) * x(i,  j-1) &
          + a(-1,i,j) * x(i-1,j)   &
          + a( 0,i,j) * x(i,  j)   &
          + a( 1,i,j) * x(i+1,j)   &
          + a( 2,i,j) * x(i,  j+1)
  END DO
  END DO

END SUBROUTINE apply2d3_mix

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE apply2d3_blocked_real4(n1, n2, a, x, out, jstart, jend)
  INTEGER, INTENT(IN)    :: n1
  INTEGER, INTENT(IN)    :: n2
  REAL(4), INTENT(IN)    :: a(-2:2, 1:n1, 1:n2)
  REAL(4), INTENT(IN)    :: x(  0:n1+1, 0:n2+1)
  REAL(4), INTENT(INOUT) :: out(0:n1+1, 0:n2+1)
  INTEGER, INTENT(IN)    :: jstart
  INTEGER, INTENT(IN)    :: jend

  INTEGER :: i, j

  DO j=jstart, jend
  DO i=1, n1
     out(i,j) = out(i,j)           &
          + a(-2,i,j) * x(i,  j-1) &
          + a(-1,i,j) * x(i-1,j)   &
          + a( 0,i,j) * x(i,  j)   &
          + a( 1,i,j) * x(i+1,j)   &
          + a( 2,i,j) * x(i,  j+1)
  END DO
  END DO

END SUBROUTINE apply2d3_blocked_real4

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE apply2d3_blocked_real8(n1, n2, a, x, out, jstart, jend)
  INTEGER, INTENT(IN)    :: n1
  INTEGER, INTENT(IN)    :: n2
  REAL(8), INTENT(IN)    :: a(-2:2, 1:n1, 1:n2)
  REAL(8), INTENT(IN)    :: x(  0:n1+1, 0:n2+1)
  REAL(8), INTENT(INOUT) :: out(0:n1+1, 0:n2+1)
  INTEGER, INTENT(IN)    :: jstart
  INTEGER, INTENT(IN)    :: jend

  INTEGER :: i, j

  DO j=jstart, jend
  DO i=1, n1
     out(i,j) = out(i,j)           &
          + a(-2,i,j) * x(i,  j-1) &
          + a(-1,i,j) * x(i-1,j)   &
          + a( 0,i,j) * x(i,  j)   &
          + a( 1,i,j) * x(i+1,j)   &
          + a( 2,i,j) * x(i,  j+1)
  END DO
  END DO

END SUBROUTINE apply2d3_blocked_real8

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE apply2d3_blocked_mix(n1, n2, a, x, out, jstart, jend)
  INTEGER, INTENT(IN)    :: n1
  INTEGER, INTENT(IN)    :: n2
  REAL(4), INTENT(IN)    :: a(-2:2, 1:n1, 1:n2)
  REAL(8), INTENT(IN)    :: x(  0:n1+1, 0:n2+1)
  REAL(8), INTENT(INOUT) :: out(0:n1+1, 0:n2+1)
  INTEGER, INTENT(IN)    :: jstart
  INTEGER, INTENT(IN)    :: jend

  INTEGER :: i, j

  DO j=jstart, jend
  DO i=1, n1
     out(i,j) = out(i,j)           &
          + a(-2,i,j) * x(i,  j-1) &
          + a(-1,i,j) * x(i-1,j)   &
          + a( 0,i,j) * x(i,  j)   &
          + a( 1,i,j) * x(i+1,j)   &
          + a( 2,i,j) * x(i,  j+1)
  END DO
  END DO

END SUBROUTINE apply2d3_blocked_mix

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE apply3d1_real4(n1, n2, n3, a, x, out, kindex)
  INTEGER, INTENT(IN)  :: n1
  INTEGER, INTENT(IN)  :: n2
  INTEGER, INTENT(IN)  :: n3
  REAL(4), INTENT(IN)  :: a(-3:3, 1:n1, 1:n2, 1:n3)
  REAL(4), INTENT(IN)  :: x(  0:n1+1, 0:n2+1, 0:n3+1)
  REAL(4), INTENT(OUT) :: out(0:n1+1, 0:n2+1, 0:n3+1)
  INTEGER, INTENT(IN), OPTIONAL :: kindex(1:n3)

  INTEGER :: i, j, k, kk

!$OMP PARALLEL DO PRIVATE(kk)
  DO k=1, n3
     IF (present(kindex)) THEN
        kk = kindex(k)
     ELSE
        kk = k
     END IF

     DO j=1, n2
     DO i=1, n1
         out(i,j,k) = a(-3,i,j,kk) * x(i,  j,  k-1) &
                    + a(-2,i,j,kk) * x(i,  j-1,k)   &
                    + a(-1,i,j,kk) * x(i-1,j,  k)   &
                    + a( 0,i,j,kk) * x(i,  j,  k)   &
                    + a( 1,i,j,kk) * x(i+1,j,  k)   &
                    + a( 2,i,j,kk) * x(i,  j+1,k)   &
                    + a( 3,i,j,kk) * x(i,  j,  k+1)
     END DO
     END DO
  END DO

END SUBROUTINE apply3d1_real4

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE apply3d1_real8(n1, n2, n3, a, x, out, kindex)
  INTEGER, INTENT(IN)  :: n1
  INTEGER, INTENT(IN)  :: n2
  INTEGER, INTENT(IN)  :: n3
  REAL(8), INTENT(IN)  :: a(-3:3, 1:n1, 1:n2, 1:n3)
  REAL(8), INTENT(IN)  :: x(  0:n1+1, 0:n2+1, 0:n3+1)
  REAL(8), INTENT(OUT) :: out(0:n1+1, 0:n2+1, 0:n3+1)
  INTEGER, INTENT(IN), OPTIONAL :: kindex(1:n3)

  INTEGER :: i, j, k, kk

!$OMP PARALLEL DO PRIVATE(kk)
  DO k=1, n3
     IF (present(kindex)) THEN
        kk = kindex(k)
     ELSE
        kk = k
     END IF

     DO j=1, n2
     DO i=1, n1
         out(i,j,k) = a(-3,i,j,kk) * x(i,  j,  k-1) &
                    + a(-2,i,j,kk) * x(i,  j-1,k)   &
                    + a(-1,i,j,kk) * x(i-1,j,  k)   &
                    + a( 0,i,j,kk) * x(i,  j,  k)   &
                    + a( 1,i,j,kk) * x(i+1,j,  k)   &
                    + a( 2,i,j,kk) * x(i,  j+1,k)   &
                    + a( 3,i,j,kk) * x(i,  j,  k+1)
     END DO
     END DO
  END DO

END SUBROUTINE apply3d1_real8

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE apply3d1_mix(n1, n2, n3, a, x, out, kindex)
  INTEGER, INTENT(IN)  :: n1
  INTEGER, INTENT(IN)  :: n2
  INTEGER, INTENT(IN)  :: n3
  REAL(4), INTENT(IN)  :: a(-3:3, 1:n1, 1:n2, 1:n3)
  REAL(8), INTENT(IN)  :: x(  0:n1+1, 0:n2+1, 0:n3+1)
  REAL(8), INTENT(OUT) :: out(0:n1+1, 0:n2+1, 0:n3+1)
  INTEGER, INTENT(IN), OPTIONAL :: kindex(1:n3)

  INTEGER :: i, j, k, kk

!$OMP PARALLEL DO PRIVATE(kk)
  DO k=1, n3
     IF (present(kindex)) THEN
        kk = kindex(k)
     ELSE
        kk = k
     END IF

     DO j=1, n2
     DO i=1, n1
         out(i,j,k) = a(-3,i,j,kk) * x(i,  j,  k-1) &
                    + a(-2,i,j,kk) * x(i,  j-1,k)   &
                    + a(-1,i,j,kk) * x(i-1,j,  k)   &
                    + a( 0,i,j,kk) * x(i,  j,  k)   &
                    + a( 1,i,j,kk) * x(i+1,j,  k)   &
                    + a( 2,i,j,kk) * x(i,  j+1,k)   &
                    + a( 3,i,j,kk) * x(i,  j,  k+1)
     END DO
     END DO
  END DO

END SUBROUTINE apply3d1_mix

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE apply3d1_blocked_real4(n1, n2, n3, a, x, out, kstart, kend, kindex)
  INTEGER, INTENT(IN)  :: n1
  INTEGER, INTENT(IN)  :: n2
  INTEGER, INTENT(IN)  :: n3
  REAL(4), INTENT(IN)  :: a(-3:3, 1:n1, 1:n2, 1:n3)
  REAL(4), INTENT(IN)  :: x(  0:n1+1, 0:n2+1, 0:n3+1)
  REAL(4), INTENT(OUT) :: out(0:n1+1, 0:n2+1, 0:n3+1)
  INTEGER, INTENT(IN)   :: kstart
  INTEGER, INTENT(IN)   :: kend
  INTEGER, INTENT(IN), OPTIONAL :: kindex(1:n3)

  INTEGER :: i, j, k, kk

  DO k=kstart, kend
     IF (present(kindex)) THEN
        kk = kindex(k)
     ELSE
        kk = k
     END IF

     DO j=1, n2
     DO i=1, n1
         out(i,j,k) = a(-3,i,j,kk) * x(i,  j,  k-1) &
                    + a(-2,i,j,kk) * x(i,  j-1,k)   &
                    + a(-1,i,j,kk) * x(i-1,j,  k)   &
                    + a( 0,i,j,kk) * x(i,  j,  k)   &
                    + a( 1,i,j,kk) * x(i+1,j,  k)   &
                    + a( 2,i,j,kk) * x(i,  j+1,k)   &
                    + a( 3,i,j,kk) * x(i,  j,  k+1)
     END DO
     END DO
  END DO

END SUBROUTINE apply3d1_blocked_real4

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE apply3d1_blocked_real8(n1, n2, n3, a, x, out, kstart, kend, kindex)
  INTEGER, INTENT(IN)  :: n1
  INTEGER, INTENT(IN)  :: n2
  INTEGER, INTENT(IN)  :: n3
  REAL(8), INTENT(IN)  :: a(-3:3, 1:n1, 1:n2, 1:n3)
  REAL(8), INTENT(IN)  :: x(  0:n1+1, 0:n2+1, 0:n3+1)
  REAL(8), INTENT(OUT) :: out(0:n1+1, 0:n2+1, 0:n3+1)
  INTEGER, INTENT(IN)  :: kstart
  INTEGER, INTENT(IN)  :: kend
  INTEGER, INTENT(IN), OPTIONAL :: kindex(1:n3)

  INTEGER :: i, j, k, kk

  DO k=kstart, kend
     IF (present(kindex)) THEN
        kk = kindex(k)
     ELSE
        kk = k
     END IF

     DO j=1, n2
     DO i=1, n1
         out(i,j,k) = a(-3,i,j,kk) * x(i,  j,  k-1) &
                    + a(-2,i,j,kk) * x(i,  j-1,k)   &
                    + a(-1,i,j,kk) * x(i-1,j,  k)   &
                    + a( 0,i,j,kk) * x(i,  j,  k)   &
                    + a( 1,i,j,kk) * x(i+1,j,  k)   &
                    + a( 2,i,j,kk) * x(i,  j+1,k)   &
                    + a( 3,i,j,kk) * x(i,  j,  k+1)
     END DO
     END DO
  END DO

END SUBROUTINE apply3d1_blocked_real8

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE apply3d1_blocked_mix(n1, n2, n3, a, x, out, kstart, kend, kindex)
  INTEGER, INTENT(IN)  :: n1
  INTEGER, INTENT(IN)  :: n2
  INTEGER, INTENT(IN)  :: n3
  REAL(4), INTENT(IN)  :: a(-3:3, 1:n1, 1:n2, 1:n3)
  REAL(8), INTENT(IN)  :: x(  0:n1+1, 0:n2+1, 0:n3+1)
  REAL(8), INTENT(OUT) :: out(0:n1+1, 0:n2+1, 0:n3+1)
  INTEGER, INTENT(IN)  :: kstart
  INTEGER, INTENT(IN)  :: kend
  INTEGER, INTENT(IN), OPTIONAL :: kindex(1:n3)

  INTEGER :: i, j, k, kk

  DO k=kstart, kend
    IF (present(kindex)) THEN
        kk = kindex(k)
     ELSE
        kk = k
     END IF

     DO j=1, n2
     DO i=1, n1
         out(i,j,k) = a(-3,i,j,kk) * x(i,  j,  k-1) &
                    + a(-2,i,j,kk) * x(i,  j-1,k)   &
                    + a(-1,i,j,kk) * x(i-1,j,  k)   &
                    + a( 0,i,j,kk) * x(i,  j,  k)   &
                    + a( 1,i,j,kk) * x(i+1,j,  k)   &
                    + a( 2,i,j,kk) * x(i,  j+1,k)   &
                    + a( 3,i,j,kk) * x(i,  j,  k+1)
     END DO
     END DO
  END DO

END SUBROUTINE apply3d1_blocked_mix

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE apply3d2_real4(n1, n2, n3, a, x, y, out, kindex)
  INTEGER, INTENT(IN)  :: n1
  INTEGER, INTENT(IN)  :: n2
  INTEGER, INTENT(IN)  :: n3
  REAL(4), INTENT(IN)  :: a(-3:3, 1:n1, 1:n2, 1:n3)
  REAL(4), INTENT(IN)  :: x(  0:n1+1, 0:n2+1, 0:n3+1)
  REAL(4), INTENT(IN)  :: y(  0:n1+1, 0:n2+1, 0:n3+1)
  REAL(4), INTENT(OUT) :: out(0:n1+1, 0:n2+1, 0:n3+1)
  INTEGER, INTENT(IN), OPTIONAL :: kindex(1:n3)

  INTEGER :: i, j, k, kk

!$OMP PARALLEL DO PRIVATE(kk)
  DO k=1, n3
    IF (present(kindex)) THEN
        kk = kindex(k)
     ELSE
        kk = k
     END IF

     DO j=1, n2
     DO i=1, n1
        out(i,j,k) = x(i,j,k)                &
             - a(-3,i,j,kk) * y(i,  j,  k-1) &
             - a(-2,i,j,kk) * y(i,  j-1,k)   &
             - a(-1,i,j,kk) * y(i-1,j,  k)   &
             - a( 0,i,j,kk) * y(i,  j,  k)   &
             - a( 1,i,j,kk) * y(i+1,j,  k)   &
             - a( 2,i,j,kk) * y(i,  j+1,k)   &
             - a( 3,i,j,kk) * y(i,  j,  k+1)
     END DO
     END DO
  END DO

END SUBROUTINE apply3d2_real4

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE apply3d2_real8(n1, n2, n3, a, x, y, out, kindex)
  INTEGER, INTENT(IN)  :: n1
  INTEGER, INTENT(IN)  :: n2
  INTEGER, INTENT(IN)  :: n3
  REAL(8), INTENT(IN)  :: a(-3:3, 1:n1, 1:n2, 1:n3)
  REAL(8), INTENT(IN)  :: x(  0:n1+1, 0:n2+1, 0:n3+1)
  REAL(8), INTENT(IN)  :: y(  0:n1+1, 0:n2+1, 0:n3+1)
  REAL(8), INTENT(OUT) :: out(0:n1+1, 0:n2+1, 0:n3+1)
  INTEGER, INTENT(IN), OPTIONAL :: kindex(1:n3)

  INTEGER :: i, j, k, kk

!$OMP PARALLEL DO PRIVATE(kk)
  DO k=1, n3
    IF (present(kindex)) THEN
        kk = kindex(k)
     ELSE
        kk = k
     END IF

     DO j=1, n2
     DO i=1, n1
        out(i,j,k) = x(i,j,k)                &
             - a(-3,i,j,kk) * y(i,  j,  k-1) &
             - a(-2,i,j,kk) * y(i,  j-1,k)   &
             - a(-1,i,j,kk) * y(i-1,j,  k)   &
             - a( 0,i,j,kk) * y(i,  j,  k)   &
             - a( 1,i,j,kk) * y(i+1,j,  k)   &
             - a( 2,i,j,kk) * y(i,  j+1,k)   &
             - a( 3,i,j,kk) * y(i,  j,  k+1)
     END DO
     END DO
  END DO

END SUBROUTINE apply3d2_real8

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE apply3d2_mix(n1, n2, n3, a, x, y, out, kindex)
  INTEGER, INTENT(IN)  :: n1
  INTEGER, INTENT(IN)  :: n2
  INTEGER, INTENT(IN)  :: n3
  REAL(4), INTENT(IN)  :: a(-3:3, 1:n1, 1:n2, 1:n3)
  REAL(8), INTENT(IN)  :: x(  0:n1+1, 0:n2+1, 0:n3+1)
  REAL(8), INTENT(IN)  :: y(  0:n1+1, 0:n2+1, 0:n3+1)
  REAL(8), INTENT(OUT) :: out(0:n1+1, 0:n2+1, 0:n3+1)
  INTEGER, INTENT(IN), OPTIONAL :: kindex(1:n3)

  INTEGER :: i, j, k, kk

!$OMP PARALLEL DO PRIVATE(kk)
  DO k=1, n3
    IF (present(kindex)) THEN
        kk = kindex(k)
     ELSE
        kk = k
     END IF

     DO j=1, n2
     DO i=1, n1
        out(i,j,k) = x(i,j,k)                &
             - a(-3,i,j,kk) * y(i,  j,  k-1) &
             - a(-2,i,j,kk) * y(i,  j-1,k)   &
             - a(-1,i,j,kk) * y(i-1,j,  k)   &
             - a( 0,i,j,kk) * y(i,  j,  k)   &
             - a( 1,i,j,kk) * y(i+1,j,  k)   &
             - a( 2,i,j,kk) * y(i,  j+1,k)   &
             - a( 3,i,j,kk) * y(i,  j,  k+1)
     END DO
     END DO
  END DO

END SUBROUTINE apply3d2_mix

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE apply3d2_blocked_real4(n1, n2, n3, a, x, y, out, kstart, kend, kindex)
  INTEGER, INTENT(IN)  :: n1
  INTEGER, INTENT(IN)  :: n2
  INTEGER, INTENT(IN)  :: n3
  REAL(4), INTENT(IN)  :: a(-3:3, 1:n1, 1:n2, 1:n3)
  REAL(4), INTENT(IN)  :: x(  0:n1+1, 0:n2+1, 0:n3+1)
  REAL(4), INTENT(IN)  :: y(  0:n1+1, 0:n2+1, 0:n3+1)
  REAL(4), INTENT(OUT) :: out(0:n1+1, 0:n2+1, 0:n3+1)
  INTEGER, INTENT(IN)  :: kstart
  INTEGER, INTENT(IN)  :: kend
  INTEGER, INTENT(IN), OPTIONAL :: kindex(1:n3)

  INTEGER :: i, j, k, kk

  DO k=kstart, kend
    IF (present(kindex)) THEN
        kk = kindex(k)
     ELSE
        kk = k
     END IF

     DO j=1, n2
     DO i=1, n1
        out(i,j,k) = x(i,j,k)                &
             - a(-3,i,j,kk) * y(i,  j,  k-1) &
             - a(-2,i,j,kk) * y(i,  j-1,k)   &
             - a(-1,i,j,kk) * y(i-1,j,  k)   &
             - a( 0,i,j,kk) * y(i,  j,  k)   &
             - a( 1,i,j,kk) * y(i+1,j,  k)   &
             - a( 2,i,j,kk) * y(i,  j+1,k)   &
             - a( 3,i,j,kk) * y(i,  j,  k+1)
     END DO
     END DO
  END DO

END SUBROUTINE apply3d2_blocked_real4

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE apply3d2_blocked_real8(n1, n2, n3, a, x, y, out, kstart, kend, kindex)
  INTEGER, INTENT(IN)  :: n1
  INTEGER, INTENT(IN)  :: n2
  INTEGER, INTENT(IN)  :: n3
  REAL(8), INTENT(IN)  :: a(-3:3, 1:n1, 1:n2, 1:n3)
  REAL(8), INTENT(IN)  :: x(  0:n1+1, 0:n2+1, 0:n3+1)
  REAL(8), INTENT(IN)  :: y(  0:n1+1, 0:n2+1, 0:n3+1)
  REAL(8), INTENT(OUT) :: out(0:n1+1, 0:n2+1, 0:n3+1)
  INTEGER, INTENT(IN)  :: kstart
  INTEGER, INTENT(IN)  :: kend
  INTEGER, INTENT(IN), OPTIONAL :: kindex(1:n3)

  INTEGER :: i, j, k, kk

  DO k=kstart, kend
    IF (present(kindex)) THEN
        kk = kindex(k)
     ELSE
        kk = k
     END IF

     DO j=1, n2
     DO i=1, n1
        out(i,j,k) = x(i,j,k)                &
             - a(-3,i,j,kk) * y(i,  j,  k-1) &
             - a(-2,i,j,kk) * y(i,  j-1,k)   &
             - a(-1,i,j,kk) * y(i-1,j,  k)   &
             - a( 0,i,j,kk) * y(i,  j,  k)   &
             - a( 1,i,j,kk) * y(i+1,j,  k)   &
             - a( 2,i,j,kk) * y(i,  j+1,k)   &
             - a( 3,i,j,kk) * y(i,  j,  k+1)
     END DO
     END DO
  END DO

END SUBROUTINE apply3d2_blocked_real8

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE apply3d2_blocked_mix(n1, n2, n3, a, x, y, out, kstart, kend, kindex)
  INTEGER, INTENT(IN)  :: n1
  INTEGER, INTENT(IN)  :: n2
  INTEGER, INTENT(IN)  :: n3
  REAL(4), INTENT(IN)  :: a(-3:3, 1:n1, 1:n2, 1:n3)
  REAL(8), INTENT(IN)  :: x(  0:n1+1, 0:n2+1, 0:n3+1)
  REAL(8), INTENT(IN)  :: y(  0:n1+1, 0:n2+1, 0:n3+1)
  REAL(8), INTENT(OUT) :: out(0:n1+1, 0:n2+1, 0:n3+1)
  INTEGER, INTENT(IN)  :: kstart
  INTEGER, INTENT(IN)  :: kend
  INTEGER, INTENT(IN), OPTIONAL :: kindex(1:n3)

  INTEGER :: i, j, k, kk

  DO k=kstart, kend
    IF (present(kindex)) THEN
        kk = kindex(k)
     ELSE
        kk = k
     END IF

     DO j=1, n2
     DO i=1, n1
        out(i,j,k) = x(i,j,k)                &
             - a(-3,i,j,kk) * y(i,  j,  k-1) &
             - a(-2,i,j,kk) * y(i,  j-1,k)   &
             - a(-1,i,j,kk) * y(i-1,j,  k)   &
             - a( 0,i,j,kk) * y(i,  j,  k)   &
             - a( 1,i,j,kk) * y(i+1,j,  k)   &
             - a( 2,i,j,kk) * y(i,  j+1,k)   &
             - a( 3,i,j,kk) * y(i,  j,  k+1)
     END DO
     END DO
  END DO

END SUBROUTINE apply3d2_blocked_mix

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE apply3d3_real4(n1, n2, n3, a, x, out, kindex)
  INTEGER, INTENT(IN)    :: n1
  INTEGER, INTENT(IN)    :: n2
  INTEGER, INTENT(IN)    :: n3
  REAL(4), INTENT(IN)    :: a(-3:3, 1:n1, 1:n2, 1:n3)
  REAL(4), INTENT(IN)    :: x(  0:n1+1, 0:n2+1, 0:n3+1)
  REAL(4), INTENT(INOUT) :: out(0:n1+1, 0:n2+1, 0:n3+1)
  INTEGER, INTENT(IN), OPTIONAL :: kindex(1:n3)

  INTEGER :: i, j, k, kk

!$OMP PARALLEL DO PRIVATE(kk)
  DO k=1, n3
    IF (present(kindex)) THEN
        kk = kindex(k)
     ELSE
        kk = k
     END IF

     DO j=1, n2
     DO i=1, n1
        out(i,j,k) = out(i,j,k)                  &
             + a(-3,i,j,kk) * x(i,  j,  k-1) &
             + a(-2,i,j,kk) * x(i,  j-1,k)   &
             + a(-1,i,j,kk) * x(i-1,j,  k)   &
             + a( 0,i,j,kk) * x(i,  j,  k)   &
             + a( 1,i,j,kk) * x(i+1,j,  k)   &
             + a( 2,i,j,kk) * x(i,  j+1,k)   &
             + a( 3,i,j,kk) * x(i,  j,  k+1)
     END DO
     END DO
  END DO

END SUBROUTINE apply3d3_real4

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE apply3d3_real8(n1, n2, n3, a, x, out, kindex)
  INTEGER, INTENT(IN)    :: n1
  INTEGER, INTENT(IN)    :: n2
  INTEGER, INTENT(IN)    :: n3
  REAL(8), INTENT(IN)    :: a(-3:3, 1:n1, 1:n2, 1:n3)
  REAL(8), INTENT(IN)    :: x(  0:n1+1, 0:n2+1, 0:n3+1)
  REAL(8), INTENT(INOUT) :: out(0:n1+1, 0:n2+1, 0:n3+1)
  INTEGER, INTENT(IN), OPTIONAL :: kindex(1:n3)

  INTEGER :: i, j, k, kk

!$OMP PARALLEL DO PRIVATE(kk)
  DO k=1, n3
    IF (present(kindex)) THEN
        kk = kindex(k)
     ELSE
        kk = k
     END IF

     DO j=1, n2
     DO i=1, n1
        out(i,j,k) = out(i,j,k)                  &
             + a(-3,i,j,kk) * x(i,  j,  k-1) &
             + a(-2,i,j,kk) * x(i,  j-1,k)   &
             + a(-1,i,j,kk) * x(i-1,j,  k)   &
             + a( 0,i,j,kk) * x(i,  j,  k)   &
             + a( 1,i,j,kk) * x(i+1,j,  k)   &
             + a( 2,i,j,kk) * x(i,  j+1,k)   &
             + a( 3,i,j,kk) * x(i,  j,  k+1)
     END DO
     END DO
  END DO

END SUBROUTINE apply3d3_real8

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE apply3d3_mix(n1, n2, n3, a, x, out, kindex)
  INTEGER, INTENT(IN)    :: n1
  INTEGER, INTENT(IN)    :: n2
  INTEGER, INTENT(IN)    :: n3
  REAL(4), INTENT(IN)    :: a(-3:3, 1:n1, 1:n2, 1:n3)
  REAL(8), INTENT(IN)    :: x(  0:n1+1, 0:n2+1, 0:n3+1)
  REAL(8), INTENT(INOUT) :: out(0:n1+1, 0:n2+1, 0:n3+1)
  INTEGER, INTENT(IN), OPTIONAL :: kindex(1:n3)

  INTEGER :: i, j, k, kk

!$OMP PARALLEL DO PRIVATE(kk)
  DO k=1, n3
    IF (present(kindex)) THEN
        kk = kindex(k)
     ELSE
        kk = k
     END IF

     DO j=1, n2
     DO i=1, n1
        out(i,j,k) = out(i,j,k)                  &
             + a(-3,i,j,kk) * x(i,  j,  k-1) &
             + a(-2,i,j,kk) * x(i,  j-1,k)   &
             + a(-1,i,j,kk) * x(i-1,j,  k)   &
             + a( 0,i,j,kk) * x(i,  j,  k)   &
             + a( 1,i,j,kk) * x(i+1,j,  k)   &
             + a( 2,i,j,kk) * x(i,  j+1,k)   &
             + a( 3,i,j,kk) * x(i,  j,  k+1)
     END DO
     END DO
  END DO

END SUBROUTINE apply3d3_mix

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE apply3d3_blocked_real4(n1, n2, n3, a, x, out, kstart, kend, kindex)
  INTEGER, INTENT(IN)    :: n1
  INTEGER, INTENT(IN)    :: n2
  INTEGER, INTENT(IN)    :: n3
  REAL(4), INTENT(IN)    :: a(-3:3, 1:n1, 1:n2, 1:n3)
  REAL(4), INTENT(IN)    :: x(  0:n1+1, 0:n2+1, 0:n3+1)
  REAL(4), INTENT(INOUT) :: out(0:n1+1, 0:n2+1, 0:n3+1)
  INTEGER, INTENT(IN)    :: kstart
  INTEGER, INTENT(IN)    :: kend
  INTEGER, INTENT(IN), OPTIONAL :: kindex(1:n3)

  INTEGER :: i, j, k, kk

  DO k=kstart, kend
    IF (present(kindex)) THEN
        kk = kindex(k)
     ELSE
        kk = k
     END IF

     DO j=1, n2
     DO i=1, n1
        out(i,j,k) = out(i,j,k)                  &
             + a(-3,i,j,kk) * x(i,  j,  k-1) &
             + a(-2,i,j,kk) * x(i,  j-1,k)   &
             + a(-1,i,j,kk) * x(i-1,j,  k)   &
             + a( 0,i,j,kk) * x(i,  j,  k)   &
             + a( 1,i,j,kk) * x(i+1,j,  k)   &
             + a( 2,i,j,kk) * x(i,  j+1,k)   &
             + a( 3,i,j,kk) * x(i,  j,  k+1)
     END DO
     END DO
  END DO

END SUBROUTINE apply3d3_blocked_real4

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE apply3d3_blocked_real8(n1, n2, n3, a, x, out, kstart, kend, kindex)
  INTEGER, INTENT(IN)    :: n1
  INTEGER, INTENT(IN)    :: n2
  INTEGER, INTENT(IN)    :: n3
  REAL(8), INTENT(IN)    :: a(-3:3, 1:n1, 1:n2, 1:n3)
  REAL(8), INTENT(IN)    :: x(  0:n1+1, 0:n2+1, 0:n3+1)
  REAL(8), INTENT(INOUT) :: out(0:n1+1, 0:n2+1, 0:n3+1)
  INTEGER, INTENT(IN)    :: kstart
  INTEGER, INTENT(IN)    :: kend
  INTEGER, INTENT(IN), OPTIONAL :: kindex(1:n3)

  INTEGER :: i, j, k, kk

  DO k=kstart, kend
    IF (present(kindex)) THEN
        kk = kindex(k)
     ELSE
        kk = k
     END IF

     DO j=1, n2
     DO i=1, n1
        out(i,j,k) = out(i,j,k)                  &
             + a(-3,i,j,kk) * x(i,  j,  k-1) &
             + a(-2,i,j,kk) * x(i,  j-1,k)   &
             + a(-1,i,j,kk) * x(i-1,j,  k)   &
             + a( 0,i,j,kk) * x(i,  j,  k)   &
             + a( 1,i,j,kk) * x(i+1,j,  k)   &
             + a( 2,i,j,kk) * x(i,  j+1,k)   &
             + a( 3,i,j,kk) * x(i,  j,  k+1)
     END DO
     END DO
  END DO

END SUBROUTINE apply3d3_blocked_real8

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE apply3d3_blocked_mix(n1, n2, n3, a, x, out, kstart, kend, kindex)
  INTEGER, INTENT(IN)    :: n1
  INTEGER, INTENT(IN)    :: n2
  INTEGER, INTENT(IN)    :: n3
  REAL(4), INTENT(IN)    :: a(-3:3, 1:n1, 1:n2, 1:n3)
  REAL(8), INTENT(IN)    :: x(  0:n1+1, 0:n2+1, 0:n3+1)
  REAL(8), INTENT(INOUT) :: out(0:n1+1, 0:n2+1, 0:n3+1)
  INTEGER, INTENT(IN)    :: kstart
  INTEGER, INTENT(IN)    :: kend
  INTEGER, INTENT(IN), OPTIONAL :: kindex(1:n3)

  INTEGER :: i, j, k, kk

  DO k=kstart, kend
    IF (present(kindex)) THEN
        kk = kindex(k)
     ELSE
        kk = k
     END IF

     DO j=1, n2
     DO i=1, n1
        out(i,j,k) = out(i,j,k)                  &
             + a(-3,i,j,kk) * x(i,  j,  k-1) &
             + a(-2,i,j,kk) * x(i,  j-1,k)   &
             + a(-1,i,j,kk) * x(i-1,j,  k)   &
             + a( 0,i,j,kk) * x(i,  j,  k)   &
             + a( 1,i,j,kk) * x(i+1,j,  k)   &
             + a( 2,i,j,kk) * x(i,  j+1,k)   &
             + a( 3,i,j,kk) * x(i,  j,  k+1)
     END DO
     END DO
  END DO

END SUBROUTINE apply3d3_blocked_mix

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE axpy2d(n1, n2, alpha, x, y)
  INTEGER, INTENT(IN)    :: n1
  INTEGER, INTENT(IN)    :: n2
  REAL(8), INTENT(IN)    :: alpha
  REAL(8), INTENT(IN)    :: x(0:n1+1,0:n2+1)
  REAL(8), INTENT(INOUT) :: y(0:n1+1,0:n2+1)

  INTEGER :: i, j

#ifdef BLAS
  CALL daxpy((n1+2)*(n2+2), alpha, x(0,0), 1, y(0,0), 1)
#else
!$OMP PARALLEL DO
  DO j=0, n2+1
  DO i=0, n1+1
     y(i,j) = y(i,j) + alpha*x(i,j)
  END DO
  END DO
#endif

END SUBROUTINE axpy2d

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE xpby2d(n1, n2, x, beta, y)
  INTEGER, INTENT(IN)    :: n1
  INTEGER, INTENT(IN)    :: n2
  REAL(8), INTENT(IN)    :: x(0:n1+1,0:n2+1)
  REAL(8), INTENT(IN)    :: beta
  REAL(8), INTENT(INOUT) :: y(0:n1+1,0:n2+1)

  INTEGER :: i, j

!$OMP PARALLEL DO
  DO j=0, n2+1
  DO i=0, n1+1
     y(i,j) = beta*y(i,j) + x(i,j)
  END DO
  END DO
END SUBROUTINE xpby2d

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE axpy3d(n1, n2, n3, alpha, x, y)
  INTEGER, INTENT(IN)    :: n1
  INTEGER, INTENT(IN)    :: n2
  INTEGER, INTENT(IN)    :: n3
  REAL(8), INTENT(IN)    :: alpha
  REAL(8), INTENT(IN)    :: x(0:n1+1,0:n2+1,0:n3+1)
  REAL(8), INTENT(INOUT) :: y(0:n1+1,0:n2+1,0:n3+1)

  INTEGER :: i, j, k

#ifdef BLAS
  CALL daxpy((n1+2)*(n2+2)*(n3+2), alpha, x(0,0,0), 1, y(0,0,0), 1)
#else
!$OMP PARALLEL DO
  DO k=0, n3+1
  DO j=0, n2+1
  DO i=0, n1+1
     y(i,j,k) = y(i,j,k) + alpha*x(i,j,k)
  END DO
  END DO
  END DO
#endif

END SUBROUTINE axpy3d

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE xpby3d(n1, n2, n3, x, beta, y)
  INTEGER, INTENT(IN)    :: n1
  INTEGER, INTENT(IN)    :: n2
  INTEGER, INTENT(IN)    :: n3
  REAL(8), INTENT(IN)    :: x(0:n1+1,0:n2+1,0:n3+1)
  REAL(8), INTENT(IN)    :: beta
  REAL(8), INTENT(INOUT) :: y(0:n1+1,0:n2+1,0:n3+1)

  INTEGER :: i, j, k

!$OMP PARALLEL DO
  DO k=0, n3+1
  DO j=0, n2+1
  DO i=0, n1+1
     y(i,j,k) = beta*y(i,j,k) + x(i,j,k)
  END DO
  END DO
  END DO

END SUBROUTINE xpby3d

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE norm2d(n1, n2, x, result, comm)
  INTEGER, INTENT(IN)  :: n1
  INTEGER, INTENT(IN)  :: n2
  REAL(8), INTENT(IN)  :: x(0:n1+1,0:n2+1)
  REAL(8), INTENT(OUT) :: result
  INTEGER, INTENT(IN), OPTIONAL :: comm

  INTEGER :: i, j

#ifdef PARALLEL_MPI
  INTEGER :: ierr
#endif

#ifdef OMP_CONSISTENT
  REAL(8) :: tmp(n2)

!$OMP PARALLEL DO shared(tmp)
  DO j=1, n2
     tmp(j) = 0.0D0
     DO i=1, n1
        tmp(j) = tmp(j) + x(i,j)**2
     END DO
  END DO
  result = sum(tmp(:))
#else
  result = 0.0D0
!$OMP PARALLEL DO REDUCTION(+:result)
  DO j=1, n2
  DO i=1, n1
     result = result + x(i,j)**2
  END DO
  END DO
#endif

#ifdef PARALLEL_MPI
  IF (present(comm)) THEN
     CALL mpi_allreduce(MPI_IN_PLACE, result, 1, MPI_REAL8, MPI_SUM, comm, ierr)
#ifdef REDUCE_CONSISTENT
     result = REAL(result, KIND=4)
#endif
  END IF
#endif
END SUBROUTINE norm2d

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE norm3d(n1, n2, n3, x, result, comm)
  INTEGER, INTENT(IN)  :: n1
  INTEGER, INTENT(IN)  :: n2
  INTEGER, INTENT(IN)  :: n3
  REAL(8), INTENT(IN)  :: x(0:n1+1,0:n2+1,0:n3+1)
  REAL(8), INTENT(OUT) :: result
  INTEGER, INTENT(IN), OPTIONAL :: comm

  INTEGER :: i, j, k

#ifdef PARALLEL_MPI
  INTEGER :: ierr
#endif

#ifdef OMP_CONSISTENT
  REAL(8) :: tmp(n3)

!$OMP PARALLEL DO shared(tmp)
  DO k=1, n3
     tmp(k) = 0.0D0
     DO j=1, n2
     DO i=1, n1
        tmp(k) = tmp(k) + x(i,j,k)**2
     END DO
     END DO
  END DO
  result = sum(tmp(:))
#else
  result = 0.0D0
!$OMP PARALLEL DO REDUCTION(+:result)
  DO k=1, n3
  DO j=1, n2
  DO i=1, n1
     result = result + x(i,j,k)**2
  END DO
  END DO
  END DO
#endif

#ifdef PARALLEL_MPI
  IF (present(comm)) THEN
     CALL mpi_allreduce(MPI_IN_PLACE, result, 1, MPI_REAL8, MPI_SUM, comm, ierr)
#ifdef REDUCE_CONSISTENT
     result = REAL(result, KIND=4)
#endif
  END IF
#endif
END SUBROUTINE norm3d

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE dot2d(n1, n2, x, y, result, comm)
  INTEGER, INTENT(IN)  :: n1
  INTEGER, INTENT(IN)  :: n2
  REAL(8), INTENT(IN)  :: x(0:n1+1,0:n2+1)
  REAL(8), INTENT(IN)  :: y(0:n1+1,0:n2+1)
  REAL(8), INTENT(OUT) :: result
  INTEGER, INTENT(IN), OPTIONAL :: comm

  INTEGER :: i, j

#ifdef PARALLEL_MPI
  INTEGER :: ierr
#endif

#ifdef OMP_CONSISTENT
  REAL(8) :: tmp(n2)

!$OMP PARALLEL DO shared(tmp)
  DO j=1, n2
     tmp(j) = 0.0D0
     DO i=1, n1
        tmp(j) = tmp(j) + x(i,j)*y(i,j)
     END DO
  END DO
  result = sum(tmp(:))
#else
  result = 0.0D0
!$OMP PARALLEL DO REDUCTION(+:result)
  DO j=1, n2
  DO i=1, n1
     result = result + x(i,j)*y(i,j)
  END DO
  END DO
#endif

#ifdef PARALLEL_MPI
  IF (present(comm)) THEN
     CALL mpi_allreduce(MPI_IN_PLACE, result, 1, MPI_REAL8, MPI_SUM, comm, ierr)
#ifdef REDUCE_CONSISTENT
     result = REAL(result, KIND=4)
#endif
  END IF
#endif
END SUBROUTINE dot2d

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE dot2d2(n1, n2, x, y1, y2, result1, result2, comm)
  INTEGER, INTENT(IN)  :: n1
  INTEGER, INTENT(IN)  :: n2
  REAL(8), INTENT(IN)  :: x( 0:n1+1,0:n2+1)
  REAL(8), INTENT(IN)  :: y1(0:n1+1,0:n2+1)
  REAL(8), INTENT(IN)  :: y2(0:n1+1,0:n2+1)
  REAL(8), INTENT(OUT) :: result1
  REAL(8), INTENT(OUT) :: result2
  INTEGER, INTENT(IN), OPTIONAL :: comm

  INTEGER :: i, j

#ifdef PARALLEL_MPI
  REAL(8) :: buf(2)
  INTEGER :: ierr
#endif

#ifdef OMP_CONSISTENT
  REAL(8) :: tmp1(n2)
  REAL(8) :: tmp2(n2)

!$OMP PARALLEL DO shared(tmp1, tmp2)
  DO j=1, n2
     tmp1(j) = 0.0D0
     tmp2(j) = 0.0D0
     DO i=1, n1
        tmp1(j) = tmp1(j) + x(i,j)*y1(i,j)
        tmp2(j) = tmp2(j) + x(i,j)*y2(i,j)
     END DO
  END DO
  result1 = sum(tmp1(:))
  result2 = sum(tmp2(:))
#else
  result1 = 0.0D0
  result2 = 0.0D0
!$OMP PARALLEL DO REDUCTION(+:result1, result2)
  DO j=1, n2
  DO i=1, n1
     result1 = result1 + x(i,j)*y1(i,j)
     result2 = result2 + x(i,j)*y2(i,j)
  END DO
  END DO
#endif

#ifdef PARALLEL_MPI
  IF (present(comm)) THEN
     buf(1) = result1
     buf(2) = result2
     CALL mpi_allreduce(MPI_IN_PLACE, buf, 2, MPI_REAL8, MPI_SUM, comm, ierr)
     result1 = buf(1)
     result2 = buf(2)
#ifdef REDUCE_CONSISTENT
     result1 = REAL(result1, KIND=4)
     result2 = REAL(result2, KIND=4)
#endif
  END IF
#endif

END SUBROUTINE dot2d2

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE dot3d(n1, n2, n3, x, y, result, comm)
  INTEGER, INTENT(IN)  :: n1
  INTEGER, INTENT(IN)  :: n2
  INTEGER, INTENT(IN)  :: n3
  REAL(8), INTENT(IN)  :: x(0:n1+1,0:n2+1,0:n3+1)
  REAL(8), INTENT(IN)  :: y(0:n1+1,0:n2+1,0:n3+1)
  REAL(8), INTENT(OUT) :: result
  INTEGER, INTENT(IN), OPTIONAL :: comm

  INTEGER :: i, j, k

#ifdef PARALLEL_MPI
  INTEGER :: ierr
#endif

#ifdef OMP_CONSISTENT
  REAL(8) :: tmp(n3)

!$OMP PARALLEL DO shared(tmp)
  DO k=1, n3
     tmp(k) = 0.0D0
     DO j=1, n2
     DO i=1, n1
        tmp(k) = tmp(k) + x(i,j,k)*y(i,j,k)
     END DO
     END DO
  END DO
  result = sum(tmp(:))
#else
  result = 0.0D0
!$OMP PARALLEL DO REDUCTION(+:result)
  DO k=1, n3
  DO j=1, n2
  DO i=1, n1
     result = result + x(i,j,k)*y(i,j,k)
  END DO
  END DO
  END DO
#endif

#ifdef PARALLEL_MPI
  IF (present(comm)) THEN
     CALL mpi_allreduce(MPI_IN_PLACE, result, 1, MPI_REAL8, MPI_SUM, comm, ierr)
#ifdef REDUCE_CONSISTENT
     result = REAL(result, KIND=4)
#endif
  END IF
#endif
END SUBROUTINE dot3d

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE dot3d2(n1, n2, n3, x, y1, y2, result1, result2, comm)
  INTEGER, INTENT(IN)  :: n1
  INTEGER, INTENT(IN)  :: n2
  INTEGER, INTENT(IN)  :: n3
  REAL(8), INTENT(IN)  :: x( 0:n1+1,0:n2+1,0:n3+1)
  REAL(8), INTENT(IN)  :: y1(0:n1+1,0:n2+1,0:n3+1)
  REAL(8), INTENT(IN)  :: y2(0:n1+1,0:n2+1,0:n3+1)
  REAL(8), INTENT(OUT) :: result1
  REAL(8), INTENT(OUT) :: result2
  INTEGER, INTENT(IN), OPTIONAL :: comm

  INTEGER :: i, j, k

#ifdef PARALLEL_MPI
  REAL(8) :: buf(2)
  INTEGER :: ierr
#endif

#ifdef OMP_CONSISTENT
  REAL(8) :: tmp1(n3)
  REAL(8) :: tmp2(n3)

!$OMP PARALLEL DO shared(tmp1, tmp2)
  DO k=1, n3
     tmp1(k) = 0.0D0
     tmp2(k) = 0.0D0
     DO j=1, n2
     DO i=1, n1
        tmp1(k) = tmp1(k) + x(i,j,k)*y1(i,j,k)
        tmp2(k) = tmp2(k) + x(i,j,k)*y2(i,j,k)
     END DO
     END DO
  END DO
  result1 = sum(tmp1(:))
  result2 = sum(tmp2(:))
#else
  result1 = 0.0D0
  result2 = 0.0D0
!$OMP PARALLEL DO REDUCTION(+:result1, result2)
  DO k=1, n3
  DO j=1, n2
  DO i=1, n1
     result1 = result1 + x(i,j,k)*y1(i,j,k)
     result2 = result2 + x(i,j,k)*y2(i,j,k)
  END DO
  END DO
  END DO
#endif

#ifdef PARALLEL_MPI
  IF (present(comm)) THEN
     buf(1) = result1
     buf(2) = result2
     CALL mpi_allreduce(MPI_IN_PLACE, buf, 2, MPI_REAL8, MPI_SUM, comm, ierr)
     result1 = buf(1)
     result2 = buf(2)
#ifdef REDUCE_CONSISTENT
     result1 = REAL(result1, KIND=4)
     result2 = REAL(result2, KIND=4)
#endif
  END IF
#endif
END SUBROUTINE dot3d2

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE cholesky(n, a, x)
  INTEGER, INTENT(IN) :: n
  REAL(8), INTENT(INOUT) :: a(n, n)
  REAL(8), INTENT(INOUT) :: x(n)

  REAL(8) :: rho
  INTEGER :: i, j, k

  a(1,1) = 1.0D0/sqrt(a(1,1))
  DO j=2, n
     a(j,1) = a(1,j)*a(1,1)
  END DO
  DO i=2, n-1
     rho = a(i,i)
     DO k=1, i-1
        rho = rho-a(i,k)**2
     END DO
     a(i,i) = 1.0D0/sqrt(rho)

     DO j=i+1, n
        rho = a(i,j)
        DO k=1, i-1
           rho = rho-a(i,k)*a(j,k)
        END DO
        a(j,i) = rho*a(i,i)
     END DO
  END DO
  rho = a(n,n)
  DO k=1, n-1
     rho = rho-a(n,k)**2
  END DO
  a(n,n) = 1.0D0/sqrt(rho)

  x(1) = x(1)*a(1,1)
  DO i=2, n
     rho = x(i)
     DO k=1, i-1
        rho = rho-a(i,k)*x(k)
     END DO
     x(i) = rho*a(i,i)
  END DO

  x(n) = x(n)*a(n,n)
  DO i=n-1, 1, -1
     rho = x(i)
     DO k = i+1, n
        rho = rho-a(k,i)*x(k)
     END DO
     x(i) = rho*a(i,i)
  END DO
END SUBROUTINE cholesky

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE invmat_real8(n, a)
  INTEGER, INTENT(IN) :: n
  REAL(8), INTENT(INOUT) :: a(n, n)

  INTEGER   :: i, j
  REAL(8) :: p, q

  DO j=1, n
     p = a(j,j)
     a(j,j) = 1.0D0
     a(j,:) = a(j,:)/p
     DO i=1, n
        IF (i == j) CYCLE
        q = a(i,j)
        a(i,j) = 0.0D0
        a(i,:) = a(i,:) - q*a(j,:)
     END DO
  END DO

END SUBROUTINE invmat_real8

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE invmat_real4(n, a)
  INTEGER, INTENT(IN) :: n
  REAL(4), INTENT(INOUT) :: a(n, n)

  INTEGER :: i, j
  REAL(4) :: p, q

  DO j=1, n
     p = a(j,j)
     a(j,j) = 1.0
     a(j,:) = a(j,:)/p
     DO i=1, n
        IF (i == j) CYCLE
        q = a(i,j)
        a(i,j) = 0.0
        a(i,:) = a(i,:) - q*a(j,:)
     END DO
  END DO

END SUBROUTINE invmat_real4

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE indicator(lev, char)
  INTEGER, INTENT(IN) :: lev
  CHARACTER(1), INTENT(IN) :: char

  INTEGER :: n

  DO n=1, lev
     WRITE(STDERR_UNIT, '(A)', ADVANCE="NO")  char
  END DO
  WRITE(STDERR_UNIT, *)
END SUBROUTINE indicator

!-----------------------------------------------------------------------------------------------------------------------

REAL(8) FUNCTION mean2d(x, comm)
  REAL(8), INTENT(IN) :: x(:,:)
  INTEGER, INTENT(IN) :: comm

  REAL(8)    :: rtmp
  INTEGER(8) :: itmp

#ifdef PARALLEL_MPI
  INTEGER :: ierr
#endif

  rtmp = sum(x)
  itmp = size(x)

#ifdef PARALLEL_MPI
  CALL mpi_allreduce(MPI_IN_PLACE, rtmp, 1, MPI_REAL8,    MPI_SUM, comm, ierr)
  CALL mpi_allreduce(MPI_IN_PLACE, itmp, 1, MPI_INTEGER8, MPI_SUM, comm, ierr)
#ifdef REDUCE_CONSISTENT
  rtmp = REAL(rtmp, KIND=4)
#endif
#endif

  mean2d = rtmp/itmp
END FUNCTION mean2d

!-----------------------------------------------------------------------------------------------------------------------

REAL(8) FUNCTION mean3d(x, comm)
  REAL(8), INTENT(IN) :: x(:,:,:)
  INTEGER, INTENT(IN) :: comm

  REAL(8)    :: rtmp
  INTEGER(8) :: itmp

#ifdef PARALLEL_MPI
  INTEGER :: ierr
#endif

  rtmp = sum(x)
  itmp = size(x)

#ifdef PARALLEL_MPI
  CALL mpi_allreduce(MPI_IN_PLACE, rtmp, 1, MPI_REAL8,    MPI_SUM, comm, ierr)
  CALL mpi_allreduce(MPI_IN_PLACE, itmp, 1, MPI_INTEGER8, MPI_SUM, comm, ierr)
#ifdef REDUCE_CONSISTENT
  rtmp = REAL(rtmp, KIND=4)
#endif
#endif

  mean3d = rtmp/itmp
END FUNCTION mean3d

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE allreduce_sum(in, out, comm)
  REAL(8), INTENT(IN)  :: in
  REAL(8), INTENT(OUT) :: out
  INTEGER, INTENT(IN)  :: comm

#ifdef PARALLEL_MPI
  INTEGER :: err

  CALL mpi_allreduce(in, out, 1, MPI_REAL8, MPI_SUM, comm, err)
#ifdef REDUCE_CONSISTENT
     out = REAL(out, KIND=4)
#endif
#else
  out = in
#endif
END SUBROUTINE allreduce_sum

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE setparam3d(n1, n2, n3, param, a, b, c, d)
  INTEGER, INTENT(IN)  :: n1, n2, n3
  REAL(8), INTENT(OUT) :: param(0:,0:,0:)
  REAL(8), INTENT(IN), OPTIONAL :: a(:,:,:)
  REAL(8), INTENT(IN), OPTIONAL :: b(:,:)
  REAL(8), INTENT(IN), OPTIONAL :: c(:)
  REAL(8), INTENT(IN), OPTIONAL :: d

  INTEGER :: s1, s2, s3
  INTEGER :: h1, h2, h3
  INTEGER :: i, j, k

  IF (present(a)) THEN
     s1 = size(a,1)
     s2 = size(a,2)
     s3 = size(a,3)
     h1 = (s1 - n1 + 1)/2
     h2 = (s2 - n2 + 1)/2
     h3 = (s3 - n3 + 1)/2

     DO k=0, n3+1
     DO j=0, n2+1
     DO i=0, n1+1
        param(i,j,k) = a(min(max(i+h1,1),s1), min(max(j+h2,1),s2), min(max(k+h3,1),s3))
     END DO
     END DO
     END DO
  ELSE IF (present(b)) THEN
     s1 = size(b,1)
     s2 = size(b,2)
     h1 = (s1 - n1 + 1)/2
     h2 = (s2 - n2 + 1)/2

     DO k=0, n3+1
     DO j=0, n2+1
     DO i=0, n1+1
        param(i,j,k) = b(min(max(i+h1,1),s1), min(max(j+h2,1),s2))
     END DO
     END DO
     END DO

  ELSE IF (present(c)) THEN
     s3 = size(c)
     h3 = (s3 - n3 + 1)/2

     DO k=0, n3+1
     DO j=0, n2+1
     DO i=0, n1+1
        param(i,j,k) = c(min(max(k+(size(c)-n3+1)/2,1),size(c)))
     END DO
     END DO
     END DO

  ELSE IF (present(d)) THEN
     DO k=0, n3+1
     DO j=0, n2+1
     DO i=0, n1+1
        param(i,j,k) = d
     END DO
     END DO
     END DO
  END IF

END SUBROUTINE setparam3d

SUBROUTINE setparam3d_logical(n1, n2, n3, param, a, b, c, d)
  INTEGER, INTENT(IN)  :: n1, n2, n3
  REAL(8), INTENT(OUT) :: param(0:,0:,0:)
  LOGICAL(1), INTENT(IN), OPTIONAL :: a(:,:,:)
  LOGICAL(1), INTENT(IN), OPTIONAL :: b(:,:)
  LOGICAL(1), INTENT(IN), OPTIONAL :: c(:)
  LOGICAL(1), INTENT(IN), OPTIONAL :: d

  INTEGER :: s1, s2, s3
  INTEGER :: h1, h2, h3
  INTEGER :: i, j, k

  IF (present(a)) THEN
     s1 = size(a,1)
     s2 = size(a,2)
     s3 = size(a,3)
     h1 = (s1 - n1 + 1)/2
     h2 = (s2 - n2 + 1)/2
     h3 = (s3 - n3 + 1)/2

     DO k=0, n3+1
     DO j=0, n2+1
     DO i=0, n1+1
        IF (a(min(max(i+h1,1),s1), min(max(j+h2,1),s2), min(max(k+h3,1),s3))) THEN
           param(i,j,k) = 1.0
        ELSE
           param(i,j,k) = 0.0
        END IF
     END DO
     END DO
     END DO

  ELSE IF (present(b)) THEN
     s1 = size(b,1)
     s2 = size(b,2)
     h1 = (s1 - n1 + 1)/2
     h2 = (s2 - n2 + 1)/2

     DO k=0, n3+1
     DO j=0, n2+1
     DO i=0, n1+1
        IF (b(min(max(i+h1,1),s1), min(max(j+h2,1),s2))) THEN
           param(i,j,k) = 1.0
        ELSE
           param(i,j,k) = 0.0
        END IF
     END DO
     END DO
     END DO

  ELSE IF (present(c)) THEN
     s3 = size(c)
     h3 = (s3 - n3 + 1)/2

     DO k=0, n3+1
     DO j=0, n2+1
     DO i=0, n1+1
        IF (c(min(max(k+(size(c)-n3+1)/2,1),size(c)))) THEN
           param(i,j,k) = 1.0
        ELSE
           param(i,j,k) = 0.0
        END IF
     END DO
     END DO
     END DO
  ELSE IF (present(d)) THEN
     DO k=0, n3+1
     DO j=0, n2+1
     DO i=0, n1+1
        IF (d) THEN
           param(i,j,k) = 1.0
        ELSE
           param(i,j,k) = 0.0
        END IF
     END DO
     END DO
     END DO
  END IF

END SUBROUTINE setparam3d_logical

SUBROUTINE setparam2d(n1, n2, param, a, d)
  INTEGER, INTENT(IN)  :: n1, n2
  REAL(8), INTENT(OUT) :: param(0:,0:)
  REAL(8), INTENT(IN), OPTIONAL :: a(:,:)
  REAL(8), INTENT(IN), OPTIONAL :: d

  INTEGER :: s1, s2
  INTEGER :: h1, h2
  INTEGER :: i, j

  IF (present(a)) THEN
     s1 = size(a,1)
     s2 = size(a,2)
     h1 = (s1 - n1 + 1)/2
     h2 = (s2 - n2 + 1)/2

     DO j=0, n2+1
     DO i=0, n1+1
        param(i,j) = a(min(max(i+h1,1),s1), min(max(j+h2,1),s2))
     END DO
     END DO
  ELSE IF (present(d)) THEN
     DO j=0, n2+1
     DO i=0, n1+1
        param(i,j) = d
     END DO
     END DO
  END IF

END SUBROUTINE setparam2d

SUBROUTINE setparam2d_logical(n1, n2, param, a, d)
  INTEGER,    INTENT(IN)  :: n1, n2
  REAL(8),    INTENT(OUT) :: param(0:,0:)
  LOGICAL(1), INTENT(IN), OPTIONAL :: a(:,:)
  LOGICAL(1), INTENT(IN), OPTIONAL :: d

  INTEGER :: s1, s2
  INTEGER :: h1, h2
  INTEGER :: i, j

  IF (present(a)) THEN
     s1 = size(a,1)
     s2 = size(a,2)
     h1 = (s1 - n1 + 1)/2
     h2 = (s2 - n2 + 1)/2

     DO j=0, n2+1
     DO i=0, n1+1
        IF (a(min(max(i+h1,1),s1), min(max(j+h2,1),s2))) THEN
           param(i,j) = 1.0
        ELSE
           param(i,j) = 0.0
        END IF
     END DO
     END DO
  ELSE IF (present(d)) THEN
     DO j=0, n2+1
     DO i=0, n1+1
        IF (d) THEN
           param(i,j) = 1.0
        ELSE
           param(i,j) = 0.0
        END IF
     END DO
     END DO
  END IF

END SUBROUTINE setparam2d_logical

END MODULE solver_util
