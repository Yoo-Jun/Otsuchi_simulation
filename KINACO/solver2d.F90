#include "macro.h"

MODULE solver2d
  USE misc
  USE parallel
  USE solver_util
  IMPLICIT NONE
  SAVE
  PRIVATE
  PUBLIC init_solver2d, solve2d, reset_solver2d, finalize_solver2d

  LOGICAL :: initialized = .FALSE.

  INTEGER :: solver_method = 2 ! 0 : pure CG,
                               ! 1 : Stand-alone MG
                               ! 2 : MGCG
                               ! 3 : MGCGS

  INTEGER :: n1
  INTEGER :: n2

  REAL(8), ALLOCATABLE :: l2d(:,:,:)

  TYPE :: mg_struct
     INTEGER :: n1 = 1
     INTEGER :: n2 = 1

#ifdef NO_ALLOCATABLE_IN_TYPE
     REAL(MG_KIND), POINTER :: r(:,:)
     REAL(MG_KIND), POINTER :: e(:,:)
     REAL(MG_KIND), POINTER :: t(:,:)

     REAL(MG_KIND), POINTER :: l(:,:,:)
     REAL(MG_KIND), POINTER :: m(:,:,:)

     INTEGER, POINTER :: p1(:)
     INTEGER, POINTER :: p2(:)

     INTEGER, POINTER :: q1(:)
     INTEGER, POINTER :: q2(:)

#else
     REAL(MG_KIND), ALLOCATABLE :: r(:,:)
     REAL(MG_KIND), ALLOCATABLE :: e(:,:)
     REAL(MG_KIND), ALLOCATABLE :: t(:,:)

     REAL(MG_KIND), ALLOCATABLE :: l(:,:,:)
     REAL(MG_KIND), ALLOCATABLE :: m(:,:,:)

     INTEGER, ALLOCATABLE :: p1(:)
     INTEGER, ALLOCATABLE :: p2(:)

     INTEGER, ALLOCATABLE :: q1(:)
     INTEGER, ALLOCATABLE :: q2(:)
#endif

#ifdef PARALLEL_MPI
     INTEGER :: rank_e
     INTEGER :: rank_w
     INTEGER :: rank_n
     INTEGER :: rank_s
#endif
  END TYPE mg_struct



  TYPE(mg_struct) :: mg(0:max_level)

  REAL(8) :: eps      = 1.0D-8
  INTEGER :: cmode    = 0
  INTEGER :: itmax    = 1000
  LOGICAL :: itreport = .FALSE.
  LOGICAL :: fallback = .FALSE.
  REAL(8) :: fallback_eps = 1.0D-16

  LOGICAL :: open_e = .FALSE.
  LOGICAL :: open_w = .FALSE.
  LOGICAL :: open_n = .FALSE.
  LOGICAL :: open_s = .FALSE.

  LOGICAL :: period1 = .FALSE.
  LOGICAL :: period2 = .FALSE.

  LOGICAL :: use_il  = .FALSE.
  LOGICAL :: singular

  INTEGER :: nlev = max_level
  INTEGER :: glev = max_level

#ifdef PARALLEL_MPI
  REAL(4), ALLOCATABLE :: gsbuf_r4(:)
  REAL(8), ALLOCATABLE :: gsbuf_r8(:)

  REAL(4), ALLOCATABLE :: sendbuf_e_r4(:)
  REAL(4), ALLOCATABLE :: sendbuf_w_r4(:)

  REAL(4), ALLOCATABLE :: sendbuf_n_r4(:)
  REAL(4), ALLOCATABLE :: sendbuf_s_r4(:)

  REAL(8), ALLOCATABLE :: sendbuf_e_r8(:)
  REAL(8), ALLOCATABLE :: sendbuf_w_r8(:)
  REAL(8), ALLOCATABLE :: sendbuf_n_r8(:)
  REAL(8), ALLOCATABLE :: sendbuf_s_r8(:)

  REAL(4), ALLOCATABLE :: recvbuf_e_r4(:)
  REAL(4), ALLOCATABLE :: recvbuf_w_r4(:)
  REAL(4), ALLOCATABLE :: recvbuf_n_r4(:)
  REAL(4), ALLOCATABLE :: recvbuf_s_r4(:)

  REAL(8), ALLOCATABLE :: recvbuf_e_r8(:)
  REAL(8), ALLOCATABLE :: recvbuf_w_r8(:)
  REAL(8), ALLOCATABLE :: recvbuf_n_r8(:)
  REAL(8), ALLOCATABLE :: recvbuf_s_r8(:)

  INTEGER :: req2d_r4(0:7, 0:max_level)
  INTEGER :: req2d_r8(0:7, 0:max_level)
#else
  INTEGER :: comm ! dummy
#endif

  REAL(MG_KIND), ALLOCATABLE :: il(:,:)
  REAL(8),       ALLOCATABLE :: work(:,:)

  INTEGER :: n_call = 0
  INTEGER :: n_iter = 0

  INTEGER :: num_threads = 1

  INTEGER :: n2_start(0:max_threads-1, 0:max_level) = 0
  INTEGER :: n2_span( 0:max_threads-1, 0:max_level) = 0

  INTERFACE sync2d
     MODULE PROCEDURE sync2d_real4
     MODULE PROCEDURE sync2d_real8
  END INTERFACE sync2d

  INTERFACE sync2d_begin
     MODULE PROCEDURE sync2d_begin_real4
     MODULE PROCEDURE sync2d_begin_real8
  END INTERFACE sync2d_begin

  INTERFACE sync2d_end
     MODULE PROCEDURE sync2d_end_real4
     MODULE PROCEDURE sync2d_end_real8
  END INTERFACE sync2d_end

#ifdef PARALLEL_MPI
  INTERFACE gather
     MODULE PROCEDURE gather_real4
     MODULE PROCEDURE gather_real8
  END INTERFACE gather

  INTERFACE scatter
     MODULE PROCEDURE scatter_real4
     MODULE PROCEDURE scatter_real8
  END INTERFACE scatter
#endif

CONTAINS

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE init_common(periods, opens)
  LOGICAL, INTENT(IN), OPTIONAL :: periods(2)
  LOGICAL, INTENT(IN), OPTIONAL :: opens(4)

  INTEGER :: n
  LOGICAL :: tmp

#ifdef PARALLEL_MPI
  INTEGER :: ierr
#endif

  CALL assert(.NOT. initialized, "SOLER2D is already initialized!")

  initialized = .TRUE.

  IF (present(periods)) THEN
     period1 = periods(1)
     period2 = periods(2)
  END IF

#ifdef PARALLEL_MPI
  ALLOCATE(gsbuf_r4(hpes*MG_GDIM**2))
  ALLOCATE(gsbuf_r8(hpes*MG_GDIM**2))

  gsbuf_r4(:) = UNDEF
  gsbuf_r8(:) = UNDEF
#endif

  IF (hrank==0) THEN
     WRITE(REPORT_UNIT, '(A)') "initialize solver2d"
  END IF

  IF (PRESENT(opens)) THEN
     open_e = opens(1)
     open_w = opens(2)
     open_n = opens(3)
     open_s = opens(4)
  END IF

  n_call = 0
  n_iter = 0

END SUBROUTINE init_common

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE init_solver2d(nx, ny, dx, dy, dsx, dsy, alpha, mask, &
                         periods, opens, method, epsilon, c_mode, it_max, report)
  INTEGER, INTENT(IN) :: nx
  INTEGER, INTENT(IN) :: ny

  REAL(8), INTENT(IN) :: dx(:,:)
  REAL(8), INTENT(IN) :: dy(:,:)

  REAL(8), INTENT(IN) :: dsx(:,:)
  REAL(8), INTENT(IN) :: dsy(:,:)

  REAL(8), INTENT(IN) :: alpha(:,:)

  LOGICAL(1), INTENT(IN) :: mask(:,:)

  LOGICAL, INTENT(IN), OPTIONAL :: periods(2)
  LOGICAL, INTENT(IN), OPTIONAL :: opens(4)

  CHARACTER(*), INTENT(IN), OPTIONAL :: method  ! default: 'MGCG'
  REAL(8),      INTENT(IN), OPTIONAL :: epsilon ! default: 1.0D-8
  INTEGER,      INTENT(IN), OPTIONAL :: c_mode  ! 0 for V-cycle, 1 for W-cycle, 2 for F-cycle
  INTEGER,      INTENT(IN), OPTIONAL :: it_max  ! default: 1000
  LOGICAL,      INTENT(IN), OPTIONAL :: report  ! default: .FALSE.

  REAL(8), ALLOCATABLE :: d1(:,:)
  REAL(8), ALLOCATABLE :: d2(:,:)

  REAL(8), ALLOCATABLE :: ds1(:,:)
  REAL(8), ALLOCATABLE :: ds2(:,:)
  REAL(8), ALLOCATABLE :: da(:,:)
  REAL(8), ALLOCATABLE :: dm(:,:)

  REAL(8), ALLOCATABLE :: dtmp(:,:)
  REAL(8), ALLOCATABLE :: ltmp(:,:,:)
  REAL(8), ALLOCATABLE :: mtmp(:,:,:)

  REAL(8) :: d1_, d2_

  INTEGER :: lev

  INTEGER :: i, j, n

  INTEGER :: tmp

#ifdef PARALLEL_MPI
  INTEGER :: ierr
#endif

  CALL init_common(periods, opens)

  n1 = nx
  n2 = ny

  mg(0)%n1 = n1
  mg(0)%n2 = n2

#ifdef PARALLEL_MPI
  ALLOCATE(d1(0:max(n1,ipes*MG_GDIM)+1, 0:max(n2,jpes*MG_GDIM)+1))
  ALLOCATE(d2(0:max(n1,ipes*MG_GDIM)+1, 0:max(n2,jpes*MG_GDIM)+1))

  ALLOCATE(ds1(0:max(n1,ipes*MG_GDIM)+1, 0:max(n2,jpes*MG_GDIM)+1))
  ALLOCATE(ds2(0:max(n1,ipes*MG_GDIM)+1, 0:max(n2,jpes*MG_GDIM)+1))

  ALLOCATE(da(0:max(n1,ipes*MG_GDIM)+1, 0:max(n2,jpes*MG_GDIM)+1))

  ALLOCATE(dm(0:max(n1,ipes*MG_GDIM)+1, 0:max(n2,jpes*MG_GDIM)+1))
#else
  ALLOCATE(d1(0:n1+1, 0:n2+1))
  ALLOCATE(d2(0:n1+1, 0:n2+1))

  ALLOCATE(ds1(0:n1+1, 0:n2+1))
  ALLOCATE(ds2(0:n1+1, 0:n2+1))

  ALLOCATE(da(0:n1+1, 0:n2+1))

  ALLOCATE(dm(0:n1+1, 0:n2+1))
#endif

  d1(:,:) = UNDEF
  d2(:,:) = UNDEF

#ifdef PARALLEL_MPI
  ALLOCATE(sendbuf_e_r4(n2))
  ALLOCATE(sendbuf_w_r4(n2))
  ALLOCATE(sendbuf_n_r4(n1))
  ALLOCATE(sendbuf_s_r4(n1))

  ALLOCATE(sendbuf_e_r8(n2))
  ALLOCATE(sendbuf_w_r8(n2))
  ALLOCATE(sendbuf_n_r8(n1))
  ALLOCATE(sendbuf_s_r8(n1))

  sendbuf_e_r4 = 0.0
  sendbuf_w_r4 = 0.0
  sendbuf_n_r4 = 0.0
  sendbuf_s_r4 = 0.0

  sendbuf_e_r8 = 0.0D0
  sendbuf_w_r8 = 0.0D0
  sendbuf_n_r8 = 0.0D0
  sendbuf_s_r8 = 0.0D0

  ALLOCATE(recvbuf_e_r4(n2))
  ALLOCATE(recvbuf_w_r4(n2))
  ALLOCATE(recvbuf_n_r4(n1))
  ALLOCATE(recvbuf_s_r4(n1))

  ALLOCATE(recvbuf_e_r8(n2))
  ALLOCATE(recvbuf_w_r8(n2))
  ALLOCATE(recvbuf_n_r8(n1))
  ALLOCATE(recvbuf_s_r8(n1))

  recvbuf_e_r4 = 0.0
  recvbuf_w_r4 = 0.0
  recvbuf_n_r4 = 0.0
  recvbuf_s_r4 = 0.0

  recvbuf_e_r8 = 0.0D0
  recvbuf_w_r8 = 0.0D0
  recvbuf_n_r8 = 0.0D0
  recvbuf_s_r8 = 0.0D0
#endif

  ALLOCATE(work(0:n1+1, 0:n2+1))
  work(:,:) = 0.0D0

  DO lev=0, max_level
     IF (lev==0) THEN
        CALL setparam2d(n1, n2, d1, dx)
        CALL setparam2d(n1, n2, d2, dy)
        CALL setparam2d(n1, n2, ds1, dsx)
        CALL setparam2d(n1, n2, ds2, dsy)
        CALL setparam2d(n1, n2, da,  alpha, 0.0D0)
        CALL setparam2d_logical(n1, n2, dm, mask, LOGICAL(.TRUE.,1))

#ifdef PARALLEL_MPI
     ELSE IF (lev == glev+1) THEN
        CALL assert(mg(lev-1)%n1==1, "DIM1 inconsistentcy for GLEV")
        CALL assert(mg(lev-1)%n2==1, "DIM2 inconsistentcy for GLEV")

        mg(lev)%n1 = mg(lev-1)%n1*ipes
        mg(lev)%n2 = mg(lev-1)%n2*jpes

        CALL gather( d1(0:mg(lev-1)%n1+1,0:mg(lev-1)%n2+1),  d1(0:mg(lev)%n1+1,0:mg(lev)%n2+1), default=1.0D0)
        CALL gather( d2(0:mg(lev-1)%n1+1,0:mg(lev-1)%n2+1),  d2(0:mg(lev)%n1+1,0:mg(lev)%n2+1), default=1.0D0)

        CALL gather(ds1(0:mg(lev-1)%n1+1,0:mg(lev-1)%n2+1), ds1(0:mg(lev)%n1+1,0:mg(lev)%n2+1), default=1.0D0)
        CALL gather(ds2(0:mg(lev-1)%n1+1,0:mg(lev-1)%n2+1), ds2(0:mg(lev)%n1+1,0:mg(lev)%n2+1), default=1.0D0)

        CALL gather( da(0:mg(lev-1)%n1+1,0:mg(lev-1)%n2+1),  da(0:mg(lev)%n1+1,0:mg(lev)%n2+1))

        CALL gather( dm(0:mg(lev-1)%n1+1,0:mg(lev-1)%n2+1),  dm(0:mg(lev)%n1+1,0:mg(lev)%n2+1))

        IF (hrank==0) THEN
           ds1(0,1:mg(lev)%n2) = 0.0D0
           ds2(1:mg(lev)%n1,0) = 0.0D0

           dm(0,           1:mg(lev)%n2) = 0.0D0
           dm(mg(lev)%n1+1,1:mg(lev)%n2) = 0.0D0
           dm(1:mg(lev)%n1,0)            = 0.0D0
           dm(1:mg(lev)%n1,mg(lev)%n2+1) = 0.0D0
        END IF

#endif
     ELSE
#ifdef PARALLEL_MPI
        IF (lev > glev) THEN
           IF (hrank==0) THEN
              d1_ = sum(d1(1:mg(lev-1)%n1,1:mg(lev-1)%n2))/(mg(lev-1)%n1*mg(lev-1)%n2)
              d2_ = sum(d2(1:mg(lev-1)%n1,1:mg(lev-1)%n2))/(mg(lev-1)%n1*mg(lev-1)%n2)
           END IF
        ELSE
           d1_ = mean2d(d1(1:mg(lev-1)%n1,1:mg(lev-1)%n2), hcomm)
           d2_ = mean2d(d2(1:mg(lev-1)%n1,1:mg(lev-1)%n2), hcomm)
        END IF
#else
        d1_ = sum(d1(1:mg(lev-1)%n1,1:mg(lev-1)%n2))/(mg(lev-1)%n1*mg(lev-1)%n2)
        d2_ = sum(d2(1:mg(lev-1)%n1,1:mg(lev-1)%n2))/(mg(lev-1)%n1*mg(lev-1)%n2)
#endif
        mg(lev-1)%p1(:) = 0
        mg(lev-1)%p2(:) = 0

        mg(lev-1)%q1(:) = 0
        mg(lev-1)%q2(:) = 0

        IF (mg(lev-1)%n1 == 1 .OR. (lev < glev      .AND. mg(lev-1)%n1 <= MG_GDIM) &
                              .OR. (d1_ > d2_*1.5D0 .AND. mg(lev-1)%n2 > 1)) THEN
           mg(lev)%n1 = mg(lev-1)%n1
           DO i=1, mg(lev)%n1
              mg(lev-1)%p1(i) = i
              mg(lev-1)%q1(i) = i
           END DO
        ELSE
           mg(lev)%n1 = mg(lev-1)%n1/2
           DO i=1, mg(lev)%n1
              mg(lev-1)%p1(i)         = 2*i
              mg(lev-1)%q1(2*i-1:2*i) =   i
           END DO

           IF (mg(lev-1)%n1 /= mg(lev)%n1*2) THEN
              mg(lev)%n1 = mg(lev)%n1+1
              mg(lev-1)%p1(mg(lev  )%n1) = mg(lev-1)%n1
              mg(lev-1)%q1(mg(lev-1)%n1) = mg(lev  )%n1
           END IF
        END IF

        IF (mg(lev-1)%n2 == 1 .OR. (lev < glev      .AND. mg(lev-1)%n2 <= MG_GDIM) &
                              .OR. (d2_ > d1_*1.5D0 .AND. mg(lev-1)%n1 > 1)) THEN
           mg(lev)%n2 = mg(lev-1)%n2
           DO j=1, mg(lev)%n2
              mg(lev-1)%p2(j) = j
              mg(lev-1)%q2(j) = j
           END DO
        ELSE
           mg(lev)%n2 = mg(lev-1)%n2/2
           DO j=1, mg(lev)%n2
              mg(lev-1)%p2(j)         = 2*j
              mg(lev-1)%q2(2*j-1:2*j) =   j
           END DO

           IF (mg(lev-1)%n2 /= mg(lev)%n2*2) THEN
              mg(lev)%n2 = mg(lev)%n2+1
              mg(lev-1)%p2(mg(lev  )%n2) = mg(lev-1)%n2
              mg(lev-1)%q2(mg(lev-1)%n2) = mg(lev  )%n2
           END IF
        END IF

        DO j=1, mg(lev)%n2
        DO i=1, mg(lev)%n1
           d1(i,j) = sum(d1(mg(lev-1)%p1(i-1)+1:mg(lev-1)%p1(i),  &
                            mg(lev-1)%p2(j-1)+1:mg(lev-1)%p2(j))) &
                    /(mg(lev-1)%p2(j)-mg(lev-1)%p2(j-1))
        END DO
        END DO

        DO j=1, mg(lev)%n2
        DO i=1, mg(lev)%n1
           d2(i,j) = sum(d2(mg(lev-1)%p1(i-1)+1:mg(lev-1)%p1(i),  &
                            mg(lev-1)%p2(j-1)+1:mg(lev-1)%p2(j))) &
                    /(mg(lev-1)%p1(i)-mg(lev-1)%p1(i-1))
        END DO
        END DO

        DO j=1, mg(lev)%n2
        DO i=0, mg(lev)%n1
           ds1(i,j) = sum(ds1(mg(lev-1)%p1(i), mg(lev-1)%p2(j-1)+1:mg(lev-1)%p2(j)))
        END DO
        END DO

        DO j=0, mg(lev)%n2
        DO i=1, mg(lev)%n1
           ds2(i,j) = sum(ds2(mg(lev-1)%p1(i-1)+1:mg(lev-1)%p1(i), mg(lev-1)%p2(j)))
        END DO
        END DO

        DO j=1, mg(lev)%n2
        DO i=1, mg(lev)%n1
           da(i,j) = sum(da(mg(lev-1)%p1(i-1)+1:mg(lev-1)%p1(i), mg(lev-1)%p2(j-1)+1:mg(lev-1)%p2(j)))
        END DO
        END DO

        dm(0,0:mg(lev)%n2+1) = 0.0D0
        dm(0:mg(lev)%n1+1,0) = 0.0D0
        DO j=1, mg(lev)%n2
        DO i=1, mg(lev)%n1
           dm(i,j) = 1.0D0 - product(1.0D0-dm(mg(lev-1)%p1(i-1)+1:mg(lev-1)%p1(i), &
                                              mg(lev-1)%p2(j-1)+1:mg(lev-1)%p2(j)))
        END DO
        END DO
        dm(  mg(lev)%n1+1,0:mg(lev)%n2+1) = 0.0D0
        dm(0:mg(lev)%n1+1,  mg(lev)%n2+1) = 0.0D0
     END IF

     IF (hrank==0) THEN
#ifdef PARALLEL_MPI
        IF (lev > glev) THEN
           WRITE(REPORT_UNIT, '("  lev",I3,": dim = (",I4,"x",I3,",",I4,"x",I3,")")') &
                lev, mg(lev)%n1, 1, mg(lev)%n2, 1
        ELSE
           WRITE(REPORT_UNIT, '("  lev",I3,": dim = (",I4,"x",I3,",",I4,"x",I3,")")') &
                lev, mg(lev)%n1, ipes, mg(lev)%n2, jpes
        END IF
#else
        WRITE(REPORT_UNIT, '("  lev",I3,": dim = (",I4,",",I4,")")') lev, mg(lev)%n1, mg(lev)%n2
#endif
     END IF

#ifdef PARALLEL_MPI
     IF (lev <= glev) THEN
        mg(lev)%rank_e = hranks(icoord+1, jcoord)
        mg(lev)%rank_w = hranks(icoord-1, jcoord)
        mg(lev)%rank_n = hranks(icoord, jcoord+1)
        mg(lev)%rank_s = hranks(icoord, jcoord-1)

        CALL mpi_recv_init(recvbuf_w_r4, mg(lev)%n2, MPI_REAL4, mg(lev)%rank_w, 1, hcomm, req2d_r4(0,lev), ierr)
        CALL mpi_recv_init(recvbuf_w_r8, mg(lev)%n2, MPI_REAL8, mg(lev)%rank_w, 1, hcomm, req2d_r8(0,lev), ierr)
        CALL mpi_recv_init(recvbuf_e_r4, mg(lev)%n2, MPI_REAL4, mg(lev)%rank_e, 2, hcomm, req2d_r4(1,lev), ierr)
        CALL mpi_recv_init(recvbuf_e_r8, mg(lev)%n2, MPI_REAL8, mg(lev)%rank_e, 2, hcomm, req2d_r8(1,lev), ierr)
        CALL mpi_recv_init(recvbuf_s_r4, mg(lev)%n1, MPI_REAL4, mg(lev)%rank_s, 3, hcomm, req2d_r4(2,lev), ierr)
        CALL mpi_recv_init(recvbuf_s_r8, mg(lev)%n1, MPI_REAL8, mg(lev)%rank_s, 3, hcomm, req2d_r8(2,lev), ierr)
        CALL mpi_recv_init(recvbuf_n_r4, mg(lev)%n1, MPI_REAL4, mg(lev)%rank_n, 4, hcomm, req2d_r4(3,lev), ierr)
        CALL mpi_recv_init(recvbuf_n_r8, mg(lev)%n1, MPI_REAL8, mg(lev)%rank_n, 4, hcomm, req2d_r8(3,lev), ierr)

        CALL mpi_send_init(sendbuf_e_r4, mg(lev)%n2, MPI_REAL4, mg(lev)%rank_e, 1, hcomm, req2d_r4(4,lev), ierr)
        CALL mpi_send_init(sendbuf_e_r8, mg(lev)%n2, MPI_REAL8, mg(lev)%rank_e, 1, hcomm, req2d_r8(4,lev), ierr)
        CALL mpi_send_init(sendbuf_w_r4, mg(lev)%n2, MPI_REAL4, mg(lev)%rank_w, 2, hcomm, req2d_r4(5,lev), ierr)
        CALL mpi_send_init(sendbuf_w_r8, mg(lev)%n2, MPI_REAL8, mg(lev)%rank_w, 2, hcomm, req2d_r8(5,lev), ierr)
        CALL mpi_send_init(sendbuf_n_r4, mg(lev)%n1, MPI_REAL4, mg(lev)%rank_n, 3, hcomm, req2d_r4(6,lev), ierr)
        CALL mpi_send_init(sendbuf_n_r8, mg(lev)%n1, MPI_REAL8, mg(lev)%rank_n, 3, hcomm, req2d_r8(6,lev), ierr)
        CALL mpi_send_init(sendbuf_s_r4, mg(lev)%n1, MPI_REAL4, mg(lev)%rank_s, 4, hcomm, req2d_r4(7,lev), ierr)
        CALL mpi_send_init(sendbuf_s_r8, mg(lev)%n1, MPI_REAL8, mg(lev)%rank_s, 4, hcomm, req2d_r8(7,lev), ierr)
     ELSE
        IF (period1) THEN
           mg(lev)%rank_e = hrank
           mg(lev)%rank_w = hrank
        ELSE
           mg(lev)%rank_e = MPI_PROC_NULL
           mg(lev)%rank_w = MPI_PROC_NULL
        END IF

        IF (period2) THEN
           mg(lev)%rank_n = hrank
           mg(lev)%rank_s = hrank
        ELSE
           mg(lev)%rank_n = MPI_PROC_NULL
           mg(lev)%rank_s = MPI_PROC_NULL
        END IF
     END IF
#endif

     IF (open_e) THEN
        DO j=0, mg(lev)%n2+1
           d1(mg(lev)%n1+1,j) = d1(mg(lev)%n1,j)
           d2(mg(lev)%n1+1,j) = d2(mg(lev)%n1,j)
           !dm(mg(lev)%n1+1,j) = dm(mg(lev)%n1,j)
           dm(mg(lev)%n1+1,j) = 0.0D0
        END DO
     END IF

     IF (open_w) THEN
        DO j=0, mg(lev)%n2+1
           d1(0,j) = d1(1,j)
           d2(0,j) = d2(1,j)
           !dm(0,j) = dm(1,j)
           dm(0,j) = 0.0D0
        END DO
     END IF

     IF (open_n) THEN
        DO i=0, mg(lev)%n1+1
           d1(i,mg(lev)%n2+1) = d1(i,mg(lev)%n2)
           d2(i,mg(lev)%n2+1) = d2(i,mg(lev)%n2)
           !dm(i,mg(lev)%n2+1) = dm(i,mg(lev)%n2)
           dm(i,mg(lev)%n2+1) = 0.0D0
        END DO
     END IF

     IF (open_s) THEN
        DO i=0, mg(lev)%n1+1
           d1(i,0) = d1(i,1)
           d2(i,0) = d2(i,1)
           !dm(i,0) = dm(i,1)
           dm(i,0) = 0.0D0
        END DO
     END IF

     ALLOCATE(mg(lev)%l(-2:2,1:mg(lev)%n1,1:mg(lev)%n2))
     ALLOCATE(mg(lev)%m(-2:2,1:mg(lev)%n1,1:mg(lev)%n2))

#ifdef PARALLEL_MPI
     IF (lev > glev .AND. hrank /= 0) THEN
#else
     IF (hrank /= 0) THEN
#endif
        mg(lev)%l(:,:,:) = 0.0D0
        mg(lev)%m(:,:,:) = 0.0D0
     ELSE
        ALLOCATE(dtmp(0:mg(lev)%n1+1,0:mg(lev)%n2+1))

        dtmp(:,:) = d1(0:mg(lev)%n1+1,0:mg(lev)%n2+1)
        CALL sync2d(lev, dtmp)
        d1(0:mg(lev)%n1+1,0:mg(lev)%n2+1) = dtmp(:,:)

        dtmp(:,:) = d2(0:mg(lev)%n1+1,0:mg(lev)%n2+1)
        CALL sync2d(lev, dtmp)
        d2(0:mg(lev)%n1+1,0:mg(lev)%n2+1) = dtmp(:,:)

        dtmp(:,:) = ds1(0:mg(lev)%n1+1,0:mg(lev)%n2+1)
        CALL sync2d(lev, dtmp)
        ds1(0:mg(lev)%n1+1,0:mg(lev)%n2+1) = dtmp(:,:)

        dtmp(:,:) = ds2(0:mg(lev)%n1+1,0:mg(lev)%n2+1)
        CALL sync2d(lev, dtmp)
        ds2(0:mg(lev)%n1+1,0:mg(lev)%n2+1) = dtmp(:,:)

        dtmp(:,:) = dm(0:mg(lev)%n1+1,0:mg(lev)%n2+1)
        CALL sync2d(lev, dtmp)
        dm(0:mg(lev)%n1+1,0:mg(lev)%n2+1) = dtmp(:,:)

        DEALLOCATE(dtmp)

        ALLOCATE(ltmp(-2:2,1:mg(lev)%n1,1:mg(lev)%n2))
        ALLOCATE(mtmp(-2:2,1:mg(lev)%n1,1:mg(lev)%n2))

#ifdef PARALLEL_MPI
        IF (mg(lev)%n1 == 1 .AND. mg(lev)%rank_e == hrank .AND. mg(lev)%rank_w == hrank) THEN
#else
        IF (mg(lev)%n1 == 1) THEN
#endif
           ltmp(-1,1,1:mg(lev)%n2) = 0.0D0
           ltmp( 1,1,1:mg(lev)%n2) = 0.0D0
        ELSE
           DO j=1, mg(lev)%n2
           DO i=1, mg(lev)%n1
              ltmp(-1,i,j) = ds1(i-1,j) * 2.0D0/(d1(i,j)+d1(i-1,j)) * dm(i,j)*dm(i-1,j)
              ltmp( 1,i,j) = ds1(i,  j) * 2.0D0/(d1(i,j)+d1(i+1,j)) * dm(i,j)*dm(i+1,j)
           END DO
           END DO
        END IF

#ifdef PARALLEL_MPI
        IF (mg(lev)%n2 == 1 .AND. mg(lev)%rank_n == hrank .AND. mg(lev)%rank_s == hrank) THEN
#else
        IF (mg(lev)%n2 == 1) THEN
#endif
           ltmp(-2,1:mg(lev)%n1,1) = 0.0D0
           ltmp( 2,1:mg(lev)%n1,1) = 0.0D0
        ELSE
           DO j=1, mg(lev)%n2
           DO i=1, mg(lev)%n1
              ltmp(-2,i,j) = ds2(i,j-1) * 2.0D0/(d2(i,j)+d2(i,j-1)) * dm(i,j)*dm(i,j-1)
              ltmp( 2,i,j) = ds2(i,j  ) * 2.0D0/(d2(i,j)+d2(i,j+1)) * dm(i,j)*dm(i,j+1)
           END DO
           END DO
        END IF

        DO j=1, mg(lev)%n2
        DO i=1, mg(lev)%n1
           ltmp( 0,i,j) = -da(i,j) - ltmp(-1,i,j) - ltmp(1,i,j) &
                                   - ltmp(-2,i,j) - ltmp(2,i,j)
        END DO
        END DO

        CALL sai2d(lev, ltmp, mtmp)

        mg(lev)%l(:,:,:) = real(ltmp(:,:,:), KIND=MG_KIND)
        mg(lev)%m(:,:,:) = real(mtmp(:,:,:), KIND=MG_KIND)
     END IF

     IF (lev==0) THEN
        ALLOCATE(l2d(-2:2,1:n1,1:n2))
        l2d(:,:,:) = ltmp(:,:,:)
     END IF

     IF (ALLOCATED(ltmp)) DEALLOCATE(ltmp)
     IF (ALLOCATED(mtmp)) DEALLOCATE(mtmp)

     ALLOCATE(mg(lev)%r(0:mg(lev)%n1+1, 0:mg(lev)%n2+1))
     ALLOCATE(mg(lev)%e(0:mg(lev)%n1+1, 0:mg(lev)%n2+1))
     ALLOCATE(mg(lev)%t(0:mg(lev)%n1+1, 0:mg(lev)%n2+1))

     mg(lev)%r(:,:) = 0.0D0
     mg(lev)%e(:,:) = 0.0D0
     mg(lev)%t(:,:) = 0.0D0

     ALLOCATE(mg(lev)%p1(0:mg(lev)%n1))
     ALLOCATE(mg(lev)%p2(0:mg(lev)%n2))

     ALLOCATE(mg(lev)%q1(1:mg(lev)%n1))
     ALLOCATE(mg(lev)%q2(1:mg(lev)%n2))

#ifdef PARALLEL_MPI
     IF (glev > lev) THEN
        IF (mg(lev)%n1 <= MG_GDIM  .AND. mg(lev)%n2 <= MG_GDIM) THEN
           glev = lev
           CALL assert(mg(glev)%n1 == MG_GDIM .AND. mg(glev)%n2 == MG_GDIM, "MG_GDIM does not fit")
        ENDIF
     ELSE
        IF (mg(lev)%n1 <= min_dim .AND. mg(lev)%n2 <= min_dim) THEN
           nlev = lev
           EXIT
        END IF
     END IF

     IF (lev > glev .AND. hrank /= 0) EXIT

#else
     IF (mg(lev)%n1 <= min_dim .AND. mg(lev)%n2 <= min_dim) THEN
        nlev = lev
        EXIT
     END IF
#endif

  END DO

  IF (hrank==0) THEN
     singular=(maxval(da(1:mg(nlev)%n1,1:mg(nlev)%n2))==0.0D0 .AND. &
               minval(da(1:mg(nlev)%n1,1:mg(nlev)%n2))==0.0D0)
  END IF

!$OMP PARALLEL
!$ num_threads = omp_get_num_threads()
!$OMP END PARALLEL
!$ CALL assert(num_threads <= max_threads, "OMP_NUM_THREADS exceeds MAX_THREADS")

  IF (present(method)) THEN
     SELECT CASE (trim(method))
     CASE ('CG')
        solver_method = 0
     CASE ('MG')
        solver_method = 1
     CASE ('MGCG')
        solver_method = 2
     CASE ('MGCGS')
        solver_method = 3
     CASE DEFAULT
        CALL assert(.FALSE., "invalid solver method '"//trim(method)//"'")
     END SELECT
  END IF

  IF (num_threads > 1 .AND. (solver_method == 1 .OR. solver_method == 2)) solver_method = -solver_method

  IF (present(epsilon)) eps      = epsilon
  IF (present(c_mode))  cmode    = c_mode
  IF (present(it_max))  itmax    = it_max
  IF (present(it_max))  itreport = report

  IF (hrank==0) THEN
     WRITE(REPORT_UNIT, '(A)', ADVANCE='NO') "solver2d using "
     SELECT CASE (solver_method)
     CASE (0)
        WRITE(REPORT_UNIT, '(A)', ADVANCE='NO') "Pure CG method"
     CASE (1)
        WRITE(REPORT_UNIT, '(A)', ADVANCE='NO') "Stand-alone MG method"
     CASE (-1)
        WRITE(REPORT_UNIT, '(A)', ADVANCE='NO') "Stand-alone MG method (OpenMP version)"
     CASE (2)
        WRITE(REPORT_UNIT, '(A)', ADVANCE='NO') "MG-CG method"
     CASE (-2)
        WRITE(REPORT_UNIT, '(A)', ADVANCE='NO') "MG-CG method (OpenMP version)"
     CASE (3)
        WRITE(REPORT_UNIT, '(A)', ADVANCE='NO') "MG-CGS method"
     END SELECT
     IF (solver_method /=0) THEN
        SELECT CASE (cmode)
        CASE (0)
           WRITE(REPORT_UNIT, '(A)', ADVANCE='NO') " with V-cycle"
        CASE (1)
           WRITE(REPORT_UNIT, '(A)', ADVANCE='NO') " with W-cycle"
        CASE (2)
           WRITE(REPORT_UNIT, '(A)', ADVANCE='NO') " with F-cycle"
        END SELECT
     END IF
     WRITE(REPORT_UNIT, '(A, ES10.3, A, I0)') ", eps =", eps, ", itmax=", itmax
  END IF

  DEALLOCATE(d1)
  DEALLOCATE(d2)

  DEALLOCATE(ds1)
  DEALLOCATE(ds2)

  DEALLOCATE(da)

  DEALLOCATE(dm)

  IF (mg(nlev)%n1 > 1 .AND. mg(nlev)%n2 > 1) THEN
!    use_il = .TRUE.
     IF (hrank==0) CALL init_il
  END IF

#ifdef PARALLEL_MPI
  CALL mpi_allreduce(MPI_IN_PLACE, nlev, 1, MPI_INTEGER, MPI_MIN, hcomm, ierr)
#endif

  DO lev=0, nlev
#ifdef PARALLEL_MPI
     CALL mpi_barrier(hcomm, ierr)
#endif

     IF (num_threads > mg(lev)%n2) THEN
        DO n=0, mg(lev)%n2-1
           n2_start(n,lev) = n+1
           n2_span(n,lev)  = 1
        END DO
        n2_start(mg(lev)%n2:num_threads-1,lev) = 0
        n2_span( mg(lev)%n2:num_threads-1,lev) = 0
     ELSE
        n2_span(0:num_threads-1,lev) = mg(lev)%n2/num_threads

        DO n=1, mg(lev)%n2 - n2_span(0,lev)*num_threads
           n2_span(num_threads-n,lev) = n2_span(num_threads-n,lev) + 1
        END DO
        DO n=0, num_threads-1
           n2_start(n,lev) = mg(lev)%n2 - sum(n2_span(n:num_threads-1,lev)) + 1
        END DO
        CALL assert(sum(n2_span(0:num_threads-1,lev))==mg(lev)%n2, "bad thread-assignment in solver2d")
     END IF

#ifdef DEBUG
!$   IF (hrank==0) THEN
!$      WRITE(REPORT_UNIT, *) "OMP_NUM_THREADS=", num_threads
!$      WRITE(REPORT_UNIT, '("thread-assignment for lev=", I2)') lev
!$      DO n=0, num_threads-1
!$         WRITE(REPORT_UNIT, '("#", I2, "  j=", I3, ":", I3)') &
!$            n, n2_start(n,lev), n2_start(n,lev)+n2_span(n,lev)-1
!$      END DO
!$   END IF
#endif
  END DO

END SUBROUTINE init_solver2d

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE reset_solver2d(dx, dy, dsx, dsy, alpha, mask)
  REAL(8), INTENT(IN) :: dx(:,:)
  REAL(8), INTENT(IN) :: dy(:,:)

  REAL(8), INTENT(IN) :: dsx(:,:)
  REAL(8), INTENT(IN) :: dsy(:,:)

  REAL(8), INTENT(IN) :: alpha(:,:)

  LOGICAL(1), INTENT(IN) :: mask(:,:)

  REAL(8) :: d1( 0:n1+1,0:n2+1)
  REAL(8) :: d2( 0:n1+1,0:n2+1)
  REAL(8) :: ds1(0:n1+1,0:n2+1)
  REAL(8) :: ds2(0:n1+1,0:n2+1)
  REAL(8) :: da( 0:n1+1,0:n2+1)
  REAL(8) :: dm(0:n1+1,0:n2+1)

  INTEGER :: i, j

#ifdef PARALLEL_MPI
  IF (comm == MPI_COMM_NULL) RETURN
#endif

  CALL assert(initialized, "solver2d is not initialized.")

  CALL setparam2d(n1, n2, d1,  dx)
  CALL setparam2d(n1, n2, d2,  dy)
  CALL setparam2d(n1, n2, ds1, dsx)
  CALL setparam2d(n1, n2, ds2, dsy)
  CALL setparam2d(n1, n2, da,  alpha, 0.0D0)
  CALL setparam2d_logical(n1, n2, dm, mask, LOGICAL(.TRUE.,1))

#ifdef PARALLEL_MPI
  IF (n1 == 1 .AND. mg(0)%rank_e == hrank .AND. mg(0)%rank_w == hrank) THEN
#else
  IF (n1 == 1) THEN
#endif
     l2d(-1,1,1:n2) = 0.0D0
     l2d( 1,1,1:n2) = 0.0D0
  ELSE
     DO j=1, n2
     DO i=1, n1
        l2d(-1,i,j) = ds1(i-1,j) * 2.0/(d1(i,j)+d1(i-1,j)) * dm(i,j)*dm(i-1,j)
        l2d( 1,i,j) = ds1(i,  j) * 2.0/(d1(i,j)+d1(i+1,j)) * dm(i,j)*dm(i+1,j)
     END DO
     END DO
  END IF

#ifdef PARALLEL_MPI
  IF (n2 == 1 .AND. mg(0)%rank_n == hrank .AND. mg(0)%rank_s == hrank) THEN
#else
  IF (n2 == 1) THEN
#endif
     l2d(-2,1:n1,1) = 0.0D0
     l2d( 2,1:n1,1) = 0.0D0
  ELSE
     DO j=1, n2
     DO i=1, n1
        l2d(-2,i,j) = ds2(i,j-1) * 2.0/(d2(i,j)+d2(i,j-1)) * dm(i,j)*dm(i,j-1)
        l2d( 2,i,j) = ds2(i,j)   * 2.0/(d2(i,j)+d2(i,j+1)) * dm(i,j)*dm(i,j+1)
     END DO
     END DO
  END IF

  DO j=1, n2
  DO i=1, n1
     l2d( 0,i,j) = -da(i,j) - l2d(-1,i,j) - l2d(1,i,j) &
                            - l2d(-2,i,j) - l2d(2,i,j)

     mg(0)%l(:,i,j) = real(l2d(:,i,j), KIND=MG_KIND)
  END DO
  END DO

CONTAINS
END SUBROUTINE reset_solver2d

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE sai2d(lev, l, m)
  INTEGER, INTENT(IN)  :: lev
  REAL(8), INTENT(IN)  :: l(-2:2, 1:mg(lev)%n1, 1:mg(lev)%n2)
  REAL(8), INTENT(OUT) :: m(-2:2, 1:mg(lev)%n1, 1:mg(lev)%n2)

  REAL(8) :: a(-6:6, -2:2), aa(-2:2, -2:2)

  REAL(8) :: tmp(-2:2,0:mg(lev)%n1+1,0:mg(lev)%n2+1)

#ifdef PARALLEL_MPI
  REAL(8) :: sendbuf1(-2:2,0:mg(lev)%n2+1)
  REAL(8) :: recvbuf1(-2:2,0:mg(lev)%n2+1)
  REAL(8) :: sendbuf2(-2:2,0:mg(lev)%n1+1)
  REAL(8) :: recvbuf2(-2:2,0:mg(lev)%n1+1)

  INTEGER :: ierr
#endif

  INTEGER :: i, j
  INTEGER :: p, q

  INTEGER :: invmat_i
  INTEGER :: invmat_j
  REAL(8) :: invmat_p
  REAL(8) :: invmat_q

  INTEGER :: n1
  INTEGER :: n2

  n1 = mg(lev)%n1
  n2 = mg(lev)%n2

  tmp(:,:,:) = 0.0D0
  tmp(-2:2,1:n1,1:n2) = l(-2:2,1:n1,1:n2)

#ifdef PARALLEL_MPI
  IF (lev > glev) THEN
     IF (hrank /= 0) RETURN

     IF (period1) THEN
        tmp(:,0,:)    = tmp(:,n1,:)
        tmp(:,n1+1,:) = tmp(:,1,:)
     END IF

     IF (period2) THEN
        tmp(:,:,0)    = tmp(:,:,n2)
        tmp(:,:,n2+1) = tmp(:,:,1)
     END IF
  ELSE
     sendbuf1(:,:) = tmp(:,n1,:)
     recvbuf1(:,:) = 0.0D0
     CALL mpi_sendrecv(sendbuf1(-2,0), (n2+2)*5, MPI_REAL8, mg(lev)%rank_e, 1, &
                       recvbuf1(-2,0), (n2+2)*5, MPI_REAL8, mg(lev)%rank_w, 1, &
                       hcomm, MPI_STATUS_IGNORE, ierr)
     tmp(:,0,:) = recvbuf1(:,:)

     sendbuf1(:,:) = tmp(:,1,:)
     recvbuf1(:,:) = 0.0D0
     CALL mpi_sendrecv(sendbuf1(-2,0), (n2+2)*5, MPI_REAL8, mg(lev)%rank_w, 2, &
                       recvbuf1(-2,0), (n2+2)*5, MPI_REAL8, mg(lev)%rank_e, 2, &
                       hcomm, MPI_STATUS_IGNORE, ierr)
     tmp(:,n1+1,:) = recvbuf1(:,:)

     sendbuf2(:,:) = tmp(:,:,n2)
     recvbuf2(:,:) = 0.0D0
     CALL mpi_sendrecv(sendbuf2(-2,0), (n1+2)*5, MPI_REAL8, mg(lev)%rank_n, 3, &
                       recvbuf2(-2,0), (n1+2)*5, MPI_REAL8, mg(lev)%rank_s, 3, &
                       hcomm, MPI_STATUS_IGNORE, ierr)
     tmp(:,:,0) = recvbuf2(:,:)

     sendbuf2(:,:) = tmp(:,:,1)
     recvbuf2(:,:) = 0.0D0
     CALL mpi_sendrecv(sendbuf2(-2,0), (n1+2)*5, MPI_REAL8, mg(lev)%rank_s, 4, &
                       recvbuf2(-2,0), (n1+2)*5, MPI_REAL8, mg(lev)%rank_n, 4, &
                       hcomm, MPI_STATUS_IGNORE, ierr)
     tmp(:,:,n2+1) = recvbuf2(:,:)
END IF

#else
  IF (period1) THEN
     tmp(:,0,:)              = tmp(:,n1,:)
     tmp(:,n1+1,:) = tmp(:,1,:)
  END IF

  IF (period2) THEN
     tmp(:,:,0)              = tmp(:,:,n2)
     tmp(:,:,n2+1) = tmp(:,:,1)
  END IF
#endif

  DO j=1, n2
     DO i=1, n1

        IF ( (.NOT. period1 .AND. n1==1) .OR. &
             (.NOT. period2 .AND. n2==1)) THEN
           m(:,i,j) = 0.0D0
           m(0,i,j) = sum(l(-2:2,i,j)**2)
           IF (m(0,i,j) > 0.0D0) THEN
              m(0,i,j) = l(0,i,j)/m(0,i,j)
           END IF
           CYCLE
        END IF

        a(:,:)  = 0.0D0

        a(-6, -2) = tmp(-2,i,j-1)
        a(-5, -2) = tmp(-1,i,j-1)
        a(-4, -2) = tmp( 0,i,j-1)
        a(-3, -2) = tmp( 1,i,j-1)
        a( 0, -2) = tmp( 2,i,j-1)

        a(-5, -1) = tmp(-2,i-1,j)
        a(-2, -1) = tmp(-1,i-1,j)
        a(-1, -1) = tmp( 0,i-1,j)
        a( 0, -1) = tmp( 1,i-1,j)
        a( 3, -1) = tmp( 2,i-1,j)

        a(-4,  0) = tmp(-2,i,j)
        a(-1,  0) = tmp(-1,i,j)
        a( 0,  0) = tmp( 0,i,j)
        a( 1,  0) = tmp( 1,i,j)
        a( 4,  0) = tmp( 2,i,j)

        a(-3,  1) = tmp(-2,i+1,j)
        a( 0,  1) = tmp(-1,i+1,j)
        a( 1,  1) = tmp( 0,i+1,j)
        a( 2,  1) = tmp( 1,i+1,j)
        a( 5,  1) = tmp( 2,i+1,j)

        a( 0,  2) = tmp(-2,i,j+1)
        a( 3,  2) = tmp(-1,i,j+1)
        a( 4,  2) = tmp( 0,i,j+1)
        a( 5,  2) = tmp( 1,i,j+1)
        a( 6,  2) = tmp( 2,i,j+1)

#ifdef CHOLESKY
        DO q=-2, 2
           DO p=-2, q
              aa(p,q) = sum(a(:,p)*a(:,q))
           END DO
           IF (aa(q,q) == 0.0D0) aa(q,q) = 1.0D0
        END DO

        m(-2,i,j) = tmp( 2,i,  j-1)
        m(-1,i,j) = tmp( 1,i-1,j)
        m( 0,i,j) = tmp( 0,i,  j)
        m( 1,i,j) = tmp(-1,i+1,j)
        m( 2,i,j) = tmp(-2,i,  j+1)

        CALL cholesky(5, aa, m(:,i,j))
#else
        aa =  matmul(transpose(a),a)
        DO q=-2, 2
           IF (aa(q,q) == 0.0D0) aa(q,q) = 1.0D0
        END DO

!        CALL invmat(5, aa)
        DO invmat_j=-2, 2
           invmat_p = aa(invmat_j,invmat_j)
           aa(invmat_j,invmat_j) = 1.0D0
           aa(invmat_j,:)        = aa(invmat_j,:)/invmat_p
           DO invmat_i=-2, 2
              IF (invmat_i == invmat_j) CYCLE
              invmat_q = aa(invmat_i,invmat_j)
              aa(invmat_i,invmat_j) = 0.0D0
              aa(invmat_i,:)        = aa(invmat_i,:) - invmat_q*aa(invmat_j,:)
           END DO
        END DO
!  invmat is not expanded in iniline form on SR16000, I don"t know why.

        DO q=-2, 2
           m(q,i,j) = tmp(0,i,j)*aa(q,0) &
                    + tmp(-1,i+1,j)*aa(q,1) + tmp(1,i-1,j)*aa(q,-1) &
                    + tmp(-2,i,j+1)*aa(q,2) + tmp(2,i,j-1)*aa(q,-2)
        END DO
#endif
     END DO
  END DO
END SUBROUTINE sai2d

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE init_il
  INTEGER :: col(-2:2)
  INTEGER :: n1, n2
  INTEGER :: i, j, l, m

  n1 = mg(nlev)%n1
  n2 = mg(nlev)%n2

  ALLOCATE(il(n1*n2,n1*n2))

  il(:,:) = 0.0_MG_KIND

  DO j=1, n2
     DO i=1, n1
        col(0) = (j-1)*n1 + i

        IF (i==1) THEN
           col(-1) = j*n1
        ELSE
           col(-1) = (j-1)*n1 + i-1
        END IF

        IF (i==n1) THEN
           col( 1) = (j-1)*n1 + 1
        ELSE
           col( 1) = (j-1)*n1 + i+1
        END IF

        IF (j==1) THEN
           col(-2) = (n2-1)*n1 + i
        ELSE
           col(-2) = (j-2)*n1 + i
        END IF

        IF (j==n2) THEN
           col( 2) = i
        ELSE
           col( 2) = j*n1 + i
        END IF

        IF (mg(nlev)%l(0,i,j) /= 0.0D0) THEN
           il(col(0),col(-2)) = mg(nlev)%l(-2,i,j)
           il(col(0),col(-1)) = mg(nlev)%l(-1,i,j)
           il(col(0),col( 0)) = mg(nlev)%l( 0,i,j)
           il(col(0),col( 1)) = mg(nlev)%l( 1,i,j)
           il(col(0),col( 2)) = mg(nlev)%l( 2,i,j)
        ELSE
           il(col(0),col( 0)) = 1.0_MG_KIND
        END IF
     END DO
  END DO
  IF (singular) il(1,:) = 1.0_MG_KIND

  CALL invmat(n1*n2, il)

END SUBROUTINE init_il

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE apply_il(n1, n2, il, in, out)
  INTEGER,       INTENT(IN)  :: n1, n2
  REAL(MG_KIND), INTENT(IN)  :: il(n1*n2, n1*n2)
  REAL(MG_KIND), INTENT(IN)  :: in( 0:n1+1, 0:n2+1)
  REAL(MG_KIND), INTENT(OUT) :: out(0:n1+1, 0:n2+1)

  REAL(MG_KIND) :: tmp(n1*n2)

  INTEGER :: i, j

  tmp(:) = reshape(in(1:n1, 1:n2), (/n1*n2/))
  IF (singular) tmp(1) = sum(out(1:n1,1:n2))

  out(1:n1,1:n2) = reshape(matmul(il,tmp), (/n1, n2/))

  IF (period1) THEN
     out(0,   1:n2) = out(n1,1:n2)
     out(n1+1,1:n2) = out(1, 1:n2)
  END IF

  IF (period2) THEN
     out(1:n1,0)    = out(1:n1,n2)
     out(1:n1,n2+1) = out(1:n1,1)
  END IF

END SUBROUTINE apply_il

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE finalize_solver2d
  INTEGER :: lev
  REAL(8) :: time
  REAL(8) :: tmp

  INTEGER :: i

  IF (.NOT. initialized) RETURN

  DEALLOCATE(work)

  DEALLOCATE(l2d)

  DO lev=0, nlev
#ifdef NO_ALLOCATABLE_IN_TYPE
     IF (ASSOCIATED(mg(lev)%r)) DEALLOCATE(mg(lev)%r)
     IF (ASSOCIATED(mg(lev)%e)) DEALLOCATE(mg(lev)%e)
     IF (ASSOCIATED(mg(lev)%t)) DEALLOCATE(mg(lev)%t)

     IF (ASSOCIATED(mg(lev)%l)) DEALLOCATE(mg(lev)%l)
     IF (ASSOCIATED(mg(lev)%m)) DEALLOCATE(mg(lev)%m)

     IF (ASSOCIATED(mg(lev)%p1)) DEALLOCATE(mg(lev)%p1)
     IF (ASSOCIATED(mg(lev)%p2)) DEALLOCATE(mg(lev)%p2)

     IF (ASSOCIATED(mg(lev)%q1)) DEALLOCATE(mg(lev)%q1)
     IF (ASSOCIATED(mg(lev)%q2)) DEALLOCATE(mg(lev)%q2)
#else
     IF (ALLOCATED(mg(lev)%r)) DEALLOCATE(mg(lev)%r)
     IF (ALLOCATED(mg(lev)%e)) DEALLOCATE(mg(lev)%e)
     IF (ALLOCATED(mg(lev)%t)) DEALLOCATE(mg(lev)%t)

     IF (ALLOCATED(mg(lev)%l)) DEALLOCATE(mg(lev)%l)
     IF (ALLOCATED(mg(lev)%m)) DEALLOCATE(mg(lev)%m)

     IF (ALLOCATED(mg(lev)%p1)) DEALLOCATE(mg(lev)%p1)
     IF (ALLOCATED(mg(lev)%p2)) DEALLOCATE(mg(lev)%p2)

     IF (ALLOCATED(mg(lev)%q1)) DEALLOCATE(mg(lev)%q1)
     IF (ALLOCATED(mg(lev)%q2)) DEALLOCATE(mg(lev)%q2)
#endif
  END DO

  IF (ALLOCATED(il)) DEALLOCATE(il)

#ifdef PARALLEL_MPI
  IF (ALLOCATED(gsbuf_r4)) DEALLOCATE(gsbuf_r4)
  IF (ALLOCATED(gsbuf_r8)) DEALLOCATE(gsbuf_r8)

  IF (ALLOCATED(sendbuf_e_r4)) DEALLOCATE(sendbuf_e_r4)
  IF (ALLOCATED(sendbuf_e_r8)) DEALLOCATE(sendbuf_e_r8)
  IF (ALLOCATED(sendbuf_w_r4)) DEALLOCATE(sendbuf_w_r4)
  IF (ALLOCATED(sendbuf_w_r8)) DEALLOCATE(sendbuf_w_r8)
  IF (ALLOCATED(sendbuf_n_r4)) DEALLOCATE(sendbuf_n_r4)
  IF (ALLOCATED(sendbuf_n_r8)) DEALLOCATE(sendbuf_n_r8)
  IF (ALLOCATED(sendbuf_s_r4)) DEALLOCATE(sendbuf_s_r4)
  IF (ALLOCATED(sendbuf_s_r8)) DEALLOCATE(sendbuf_s_r8)

  IF (ALLOCATED(recvbuf_e_r4)) DEALLOCATE(recvbuf_e_r4)
  IF (ALLOCATED(recvbuf_e_r8)) DEALLOCATE(recvbuf_e_r8)
  IF (ALLOCATED(recvbuf_w_r4)) DEALLOCATE(recvbuf_w_r4)
  IF (ALLOCATED(recvbuf_w_r8)) DEALLOCATE(recvbuf_w_r8)
  IF (ALLOCATED(recvbuf_n_r4)) DEALLOCATE(recvbuf_n_r4)
  IF (ALLOCATED(recvbuf_n_r8)) DEALLOCATE(recvbuf_n_r8)
  IF (ALLOCATED(recvbuf_s_r4)) DEALLOCATE(recvbuf_s_r4)
  IF (ALLOCATED(recvbuf_s_r8)) DEALLOCATE(recvbuf_s_r8)
#endif

  IF (hrank==0) THEN
     IF (n_call >= 1) THEN
        WRITE(REPORT_UNIT, '(A,ES10.3)') "finalize solver2d: average iteration =", DBLE(n_iter)/DBLE(n_call)
     END IF
  END IF

  initialized = .FALSE.
  eps = 1.0D-8
  nlev = max_level

END SUBROUTINE finalize_solver2d

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE solve2d(q, p)
  REAL(8), INTENT(IN)    :: q(0:n1+1, 0:n2+1)
  REAL(8), INTENT(INOUT) :: p(0:n1+1, 0:n2+1)

#ifdef DEBUG__
  REAL(8) :: sumq
  INTEGER :: ierr
#endif

#ifdef PARALLEL_MPI
  IF (comm == MPI_COMM_NULL) RETURN
#endif

#ifdef DEBUG__
  CALL assert(initialized, "solver2d is not initialized.")

  CALL mpi_reduce(sum(q(1:n1, 1:n2)), sumq, 1, MPI_REAL8, MPI_SUM, 0, hcomm, ierr)
  IF (hrank==0) WRITE(REPORT_UNIT, '("SUM(Q)=",ES10.3)') sumq
#endif

  SELECT CASE(solver_method)
     CASE(0)
        CALL solver_cg2d(q, p)
     CASE(1)
        CALL solver_mg2d(q, p)
     CASE(2)
        CALL solver_mgcg2d(q, p)
     CASE(3)
        CALL solver_mgcgs2d(q, p)
     CASE(-1)
        CALL solver_mg2d_omp(q, p)
     CASE(-2)
        CALL solver_mgcg2d_omp(q, p)
  END SELECT

END SUBROUTINE solve2d

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE solver_mgcg2d(q, p)
  REAL(8), INTENT(IN)    :: q(0:n1+1, 0:n2+1)
  REAL(8), INTENT(INOUT) :: p(0:n1+1, 0:n2+1)

  REAL(8) :: res, res_c

  REAL(8) :: alpha, beta, rho

  REAL(8) :: u(0:n1+1, 0:n2+1)
  REAL(8) :: v(0:n1+1, 0:n2+1)
  REAL(8) :: r(0:n1+1, 0:n2+1)
  REAL(8) :: e(0:n1+1, 0:n2+1)

  INTEGER :: n


!$OMP PARALLEL WORKSHARE
  v(:,:) = 0.0D0
  r(:,:) = 0.0D0
  e(:,:) = 0.0D0
!$OMP END PARALLEL WORKSHARE

  CALL norm2d(n1, n2, q, res_c, comm)
  res_c = res_c *eps*eps

  CALL apply2d2(n1, n2, l2d, q, p, r)

  CALL sync2d_begin(0, r)

  CALL norm2d(n1, n2, r, res, comm)

  CALL sync2d_end(0, r)

  IF (res <= res_c) THEN
     CALL checkout(0, res, res_c)
     RETURN
  END IF

  CALL multigrid2d(r, e)

  CALL sync2d(0, e)

!$OMP PARALLEL WORKSHARE
  u(:,:) = e(:,:)
!$OMP END PARALLEL WORKSHARE

  DO n=1, itmax
     CALL apply2d1(n1, n2, l2d, u, v)

     CALL sync2d_begin(0, v)

     CALL dot2d2(n1, n2, u, v, r, rho, alpha, comm)

     CALL sync2d_end(0, v)

     IF (fallback .AND. abs(rho) < fallback_eps) THEN
        IF (itreport .AND. hrank==0) WRITE(REPORT_UNIT, '(A, I3)') "MGCG unstable, fallback to stand-alone multigrid: n=", n

        CALL solver_mg2d(q, p)
        RETURN
     END IF

     alpha = alpha / rho

     CALL axpy2d(n1, n2,  alpha, u, p)
     CALL axpy2d(n1, n2, -alpha, v, r)

     CALL norm2d(n1, n2, r, res, comm)

     IF (res <= res_c) THEN
        CALL checkout(n, res, res_c)
        RETURN
     END IF

     CALL multigrid2d(r, e)

     CALL sync2d_begin(0, e)

     CALL dot2d(n1, n2, v, e, beta, comm)
     beta = -beta / rho

     CALL sync2d_end(0, e)

     CALL xpby2d(n1, n2, e, beta, u)
  END DO

  CALL quit_message(res, res_c)

END SUBROUTINE solver_mgcg2d

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE solver_mgcg2d_omp(q, p)
  REAL(8), INTENT(IN)    :: q(0:n1+1, 0:n2+1)
  REAL(8), INTENT(INOUT) :: p(0:n1+1, 0:n2+1)

  REAL(8) :: res, res_c

  REAL(8) :: alpha, beta, rho

  REAL(8) :: u(0:n1+1, 0:n2+1)
  REAL(8) :: v(0:n1+1, 0:n2+1)
  REAL(8) :: r(0:n1+1, 0:n2+1)
  REAL(8) :: e(0:n1+1, 0:n2+1)

  INTEGER :: n
  INTEGER :: i, j

  INTEGER, SAVE :: tid = 0
  INTEGER, SAVE :: jstart, jend, start_offset, end_offset
!$OMP threadprivate(tid, jstart, jend, start_offset, end_offset)

  REAL(8) :: tmp(n2,2)

!$OMP PARALLEL
!$ tid = omp_get_thread_num()
  jstart = n2_start(tid,0)
  jend   = n2_start(tid,0) + n2_span(tid,0) - 1

  IF (jstart==1) THEN
     start_offset = 1
  ELSE
     start_offset = 0
  END IF

  IF (jend==n2) THEN
     end_offset = 1
  ELSE
     end_offset = 0
  END IF

  v(:,jstart-start_offset:jend+end_offset) = 0.0D0
  r(:,jstart-start_offset:jend+end_offset) = 0.0D0
  e(:,jstart-start_offset:jend+end_offset) = 0.0D0

! CALL norm2d(n1, n2, q, res_c, comm)
  DO j=jstart, jend
     tmp(j,1) = 0.0D0
     DO i=1, n1
        tmp(j,1) = tmp(j,1) + q(i,j)**2
     END DO
  END DO
!$OMP BARRIER
!$OMP SINGLE
  CALL allreduce_sum(sum(tmp(1:n2,1)), res_c, comm)
  res_c = res_c *eps*eps
!$OMP END SINGLE nowait

  CALL apply2d2(n1, n2, l2d, q, p, r, jstart, jend)

  CALL sync2d_begin(0, r, tid)

! CALL norm2d(n1, n2, r, res, comm)
  DO j=jstart, jend
     tmp(j,1) = 0.0D0
     DO i=1, n1
        tmp(j,1) = tmp(j,1) + r(i,j)**2
     END DO
  END DO
!$OMP BARRIER
!$OMP SINGLE
  CALL allreduce_sum(sum(tmp(1:n2,1)), res, comm)
!$OMP END SINGLE

  CALL sync2d_end(0, r, tid)
!$OMP END PARALLEL

  IF (res <= res_c) THEN
     CALL checkout(0, res, res_c)
     RETURN
  END IF

!$OMP PARALLEL
  CALL multigrid2d_omp(r, e, tid)

  CALL sync2d(0, e, tid)

  u(:,jstart-start_offset:jend+end_offset) = e(:,jstart-start_offset:jend+end_offset)
!$OMP END PARALLEL

  DO n=1, itmax
!$OMP PARALLEL
     CALL apply2d1(n1, n2, l2d, u, v, jstart, jend)

     CALL sync2d_begin(0, v, tid)

!    CALL dot2d2(n1, n2, u, v, r, rho, alpha, comm)
     DO j=jstart, jend
        tmp(j,1) = 0.0D0
        tmp(j,2) = 0.0D0
        DO i=1, n1
           tmp(j,1) = tmp(j,1) + u(i,j)*v(i,j)
           tmp(j,2) = tmp(j,2) + u(i,j)*r(i,j)
        END DO
     END DO
!$OMP BARRIER
!$OMP SINGLE
     CALL allreduce_sum(sum(tmp(1:n2,1)), rho, comm)
     CALL allreduce_sum(sum(tmp(1:n2,2)), alpha, comm)
!$OMP END SINGLE

     CALL sync2d_end(0, v, tid)
!$OMP END PARALLEL
     IF (fallback .AND. abs(rho) < fallback_eps) THEN
        IF (itreport .AND. hrank==0) WRITE(REPORT_UNIT, '(A, I3)') "MGCG unstable, fallback to stand-alone multigrid: n=", n

        CALL solver_mg2d_omp(q, p)
        RETURN
     END IF

     alpha = alpha / rho

!$OMP PARALLEL
!    CALL axpy2d(n1, n2,  alpha, u, p, jstart-start_offset, jend+end_offset)
!    CALL axpy2d(n1, n2, -alpha, v, r, jstart-start_offset, jend+end_offset)
     p(:,jstart-start_offset:jend+end_offset) = p(:,jstart-start_offset:jend+end_offset) + alpha*u(:,jstart-start_offset:jend+end_offset)
     r(:,jstart-start_offset:jend+end_offset) = r(:,jstart-start_offset:jend+end_offset) - alpha*v(:,jstart-start_offset:jend+end_offset)

!    CALL norm2d(n1, n2, r, res, comm)
     DO j=jstart, jend
        tmp(j,1) = 0.0D0
        DO i=1, n1
           tmp(j,1) = tmp(j,1) + r(i,j)**2
        END DO
     END DO
!$OMP BARRIER
!$OMP SINGLE
     CALL allreduce_sum(sum(tmp(1:n2,1)), res, comm)
!$OMP END SINGLE
!$OMP END PARALLEL
     IF (res <= res_c) THEN
        CALL checkout(n, res, res_c)
        RETURN
     END IF

!$OMP PARALLEL
     CALL multigrid2d_omp(r, e, tid)

     CALL sync2d_begin(0, e, tid)

!    CALL dot2d(n1, n2, v, e, beta, comm)
     DO j=jstart, jend
        tmp(j,1) = 0.0D0
        DO i=1, n1
           tmp(j,1) = tmp(j,1) + v(i,j)*e(i,j)
        END DO
     END DO
!$OMP BARRIER
!$OMP SINGLE
     CALL allreduce_sum(sum(tmp(1:n2,1)), beta, comm)
     beta = -beta / rho
!$OMP END SINGLE
     CALL sync2d_end(0, e, tid)

!    CALL xpby2d(n1, n2, e, beta, u)
     u(:,jstart-start_offset:jend+end_offset) = beta*u(:,jstart-start_offset:jend+end_offset) + e(:,jstart-start_offset:jend+end_offset)
!$OMP END PARALLEL
  END DO

  CALL quit_message(res, res_c)

END SUBROUTINE solver_mgcg2d_omp

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE solver_mgcgs2d(q, p)
  REAL(8), INTENT(IN)    :: q(0:n1+1, 0:n2+1)
  REAL(8), INTENT(INOUT) :: p(0:n1+1, 0:n2+1)

  REAL(8) :: res, res_c

  REAL(8) :: alpha, beta, rho

  REAL(8) :: r(0:n1+1, 0:n2+1)
  REAL(8) :: e(0:n1+1, 0:n2+1)

  REAL(8) :: s(0:n1+1, 0:n2+1)
  REAL(8) :: t(0:n1+1, 0:n2+1)
  REAL(8) :: u(0:n1+1, 0:n2+1)
  REAL(8) :: v(0:n1+1, 0:n2+1)

  REAL(8) :: r0(0:n1+1, 0:n2+1)

  INTEGER :: n

!$OMP PARALLEL WORKSHARE
  r0(:,:) = 0.0D0

  s(:,:) = 0.0D0
  t(:,:) = 0.0D0
  e(:,:) = 0.0D0
!$OMP END PARALLEL WORKSHARE

  CALL norm2d(n1, n2, q, res_c, comm)
  res_c = res_c *eps*eps

  CALL apply2d2(n1, n2, l2d, q, p, r0)

  CALL sync2d_begin(0, r0)

  CALL norm2d(n1, n2, r0, res, comm)

  CALL sync2d_end(0, r0)

  IF (res <= res_c) THEN
     CALL checkout(0, res, res_c)
     RETURN
  END IF

!$OMP PARALLEL WORKSHARE
  u(:,:) = r0(:,:)
  v(:,:) = r0(:,:)
  r(:,:) = r0(:,:)
!$OMP END PARALLEL WORKSHARE

  rho = res

  DO n=1, itmax
     IF (fallback .AND. abs(rho) < fallback_eps) THEN
        IF (itreport .AND. hrank==0) WRITE(REPORT_UNIT, '(A,I3)') "MGCGS unstable, fallback to stand-alone multigrid: n=", n

        CALL solver_mg2d(q, p)
        RETURN
     END IF

     CALL multigrid2d(u, e)

     CALL apply2d1(n1, n2, l2d, e, s)

     CALL sync2d_begin(0, s)

     CALL dot2d(n1, n2, r0, s, alpha, comm)

     alpha = rho / alpha

     CALL sync2d_end(0, s)

!$OMP PARALLEL WORKSHARE
     t(:,:) = v(:,:) - alpha * s(:,:)
!$OMP END PARALLEL WORKSHARE

     CALL multigrid2d(v + t, e)

     CALL axpy2d(n1, n2, alpha, e, p)

     CALL apply2d1(n1, n2, l2d, e, s)

     CALL sync2d(0, s)

     CALL axpy2d(n1, n2, -alpha, s, r)

     CALL norm2d(n1, n2, r, res, comm)

     IF (res <= res_c) THEN
        CALL checkout(n, res, res_c)
        RETURN
     END IF

     beta = rho
     CALL dot2d(n1, n2, r, r0, rho, comm)
     beta = rho/beta

!$OMP PARALLEL WORKSHARE
     v(:,:) = r(:,:) + beta * t(:,:)
     u(:,:) = v(:,:) + beta * (t(:,:) + beta * u(:,:))
!$OMP END PARALLEL WORKSHARE
  END DO

  CALL quit_message(res, res_c)

END SUBROUTINE solver_mgcgs2d

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE solver_cg2d(q, p)
  REAL(8), INTENT(IN)    :: q(0:n1+1, 0:n2+1)
  REAL(8), INTENT(INOUT) :: p(0:n1+1, 0:n2+1)

  REAL(8) :: res, res_c

  REAL(8) :: alpha, beta, rho

  REAL(8) :: u(0:n1+1, 0:n2+1)
  REAL(8) :: v(0:n1+1, 0:n2+1)

  REAL(8) :: r(0:n1+1, 0:n2+1)
  REAL(8) :: e(0:n1+1, 0:n2+1)

  INTEGER :: n

!$OMP PARALLEL WORKSHARE
  v(:,:) = 0.0D0
  r(:,:) = 0.0D0
!$OMP END PARALLEL WORKSHARE

  CALL norm2d(n1, n2, q, res_c, comm)
  res_c = res_c *eps*eps

  CALL apply2d2(n1, n2, l2d, q, p, r)

  CALL sync2d_begin(0, r)

  CALL norm2d(n1, n2, r, res, comm)

  CALL sync2d_end(0, r)

  IF (res <= res_c) THEN
     CALL checkout(0, res, res_c)
     RETURN
  END IF

  CALL sync2d(0, r)

!$OMP PARALLEL WORKSHARE
  u(:,:) = r(:,:)
!$OMP END PARALLEL WORKSHARE

  DO n=1, itmax
     CALL apply2d1(n1, n2, l2d, u, v)

     CALL sync2d_begin(0, v)

     CALL dot2d2(n1, n2, u, v, r, rho, alpha, comm)
     alpha = alpha / rho

     CALL sync2d_end(0, v)

     CALL axpy2d(n1, n2,  alpha, u, p)
     CALL axpy2d(n1, n2, -alpha, v, r)

     CALL norm2d(n1, n2, r, res, comm)

     IF (res <= res_c) THEN
        CALL checkout(n, res, res_c)
        RETURN
     END IF

     CALL sync2d_begin(0, r)

     CALL dot2d(n1, n2, v, r, beta, comm)
     beta = -beta / rho

     CALL sync2d_end(0, r)

     CALL xpby2d(n1, n2, r, beta, u)
  END DO

  CALL quit_message(res, res_c)

END SUBROUTINE solver_cg2d

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE solver_mg2d(q, p)
  REAL(8), INTENT(IN)    :: q(0:n1+1, 0:n2+1)
  REAL(8), INTENT(INOUT) :: p(0:n1+1, 0:n2+1)

  REAL(8) :: res, res_c

  REAL(8) :: r(0:n1+1, 0:n2+1)
  REAL(8) :: e(0:n1+1, 0:n2+1)

  INTEGER :: n

!$OMP PARALLEL WORKSHARE
  r(:,:) = 0.0D0
  e(:,:) = 0.0D0
!$OMP END PARALLEL WORKSHARE

  CALL norm2d(n1, n2, q, res_c, comm)
  res_c = res_c *eps*eps

  DO n=1, itmax
     CALL apply2d2(n1, n2, l2d, q, p, r)

     CALL sync2d_begin(0, r)

     CALL norm2d(n1, n2, r, res, comm)

     CALL sync2d_end(0, r)

     IF (res <= res_c) THEN
        CALL checkout(n, res, res_c)
        RETURN
     END IF

     CALL multigrid2d(r, e)

     CALL sync2d(0, e)

!$OMP PARALLEL WORKSHARE
     p(:,:) = p(:,:) + e(:,:)
!$OMP END PARALLEL WORKSHARE
  END DO

  CALL quit_message(res, res_c)

END SUBROUTINE solver_mg2d

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE solver_mg2d_omp(q, p)
  REAL(8), INTENT(IN)    :: q(0:n1+1, 0:n2+1)
  REAL(8), INTENT(INOUT) :: p(0:n1+1, 0:n2+1)

  REAL(8) :: res, res_c

  REAL(8) :: r(0:n1+1, 0:n2+1)
  REAL(8) :: e(0:n1+1, 0:n2+1)

  INTEGER :: n
  INTEGER :: i, j

  INTEGER, SAVE :: tid = 0
  INTEGER, SAVE :: jstart, jend, start_offset, end_offset
!$OMP threadprivate(tid, jstart, jend, start_offset, end_offset)

  REAL(8) :: tmp(n2)

!$OMP PARALLEL
!$ tid = omp_get_thread_num()
  jstart = n2_start(tid,0)
  jend   = n2_start(tid,0) + n2_span(tid,0) - 1

  IF (jstart==1) THEN
     start_offset = 1
  ELSE
     start_offset = 0
  END IF

  IF (jend==n2) THEN
     end_offset = 1
  ELSE
     end_offset = 0
  END IF

  r(:,jstart-start_offset:jend+end_offset) = 0.0D0
  e(:,jstart-start_offset:jend+end_offset) = 0.0D0

! CALL norm2d(n1, n2, q, res_c, comm)
  DO j=jstart, jend
     tmp(j) = 0.0D0
     DO i=1, n1
        tmp(j) = tmp(j) + q(i,j)**2
     END DO
  END DO

!$OMP BARRIER
!$OMP SINGLE
  CALL allreduce_sum(sum(tmp(1:n2)), res_c, comm)
  res_c = res_c *eps*eps
!$OMP END SINGLE
!$OMP END PARALLEL

  DO n=1, itmax
!$OMP PARALLEL
     CALL apply2d2(n1, n2, l2d, q, p, r, jstart, jend)

     CALL sync2d_begin(0, r, tid)

!    CALL norm2d(n1, n2, r, res, comm)
     DO j=jstart, jend
        tmp(j) = 0.0D0
        DO i=1, n1
           tmp(j) = tmp(j) + r(i,j)**2
        END DO
     END DO
!$OMP BARRIER
!$OMP SINGLE
     CALL allreduce_sum(sum(tmp(1:n2)), res, comm)
!$OMP END SINGLE

     CALL sync2d_end(0, r, tid)
!$OMP END PARALLEL
     IF (res <= res_c) THEN
        CALL checkout(n, res, res_c)
        RETURN
     END IF

!$OMP PARALLEL
     CALL multigrid2d_omp(r, e, tid)

     CALL sync2d(0, e, tid)

!    CALL axpy2d(n1, n2, 1.0D0, e, p, jstart-start_offset, jend+end_offset)
     p(:,jstart-start_offset:jend+end_offset) = p(:,jstart-start_offset:jend+end_offset) + e(:,jstart-start_offset:jend+end_offset)

!$OMP END PARALLEL
  END DO

  CALL quit_message(res, res_c)

END SUBROUTINE solver_mg2d_omp

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE checkout(n, res, res_c)
  INTEGER, INTENT(IN) :: n
  REAL(8), INTENT(IN) :: res
  REAL(8), INTENT(IN) :: res_c

  IF (itreport .AND. res > 0 .AND. hrank==0)  THEN
     WRITE(REPORT_UNIT, '(A,I5,A,ES10.3)') "solver2d: iteration count=", n, &
                                           ", relative residual=", sqrt(res/res_c)*eps
  END IF

  IF (n > 0) THEN
     n_call = n_call + 1
     n_iter = n_iter + n
  END IF
END SUBROUTINE checkout

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE quit_message(res, res_c)
  REAL(8), INTENT(IN) :: res
  REAL(8), INTENT(IN) :: res_c

  IF (itreport .AND. hrank==0)  THEN
     WRITE(REPORT_UNIT, '(A,I5,A,ES10.3)') "solver2d: iteration count exceeds ", itmax, &
                                           ", relative residual=", sqrt(res/res_c)*eps
  END IF

  n_call = n_call + 1
  n_iter = n_iter + itmax

END SUBROUTINE quit_message

!-----------------------------------------------------------------------------------------------------------------------

RECURSIVE SUBROUTINE mgsmooth2d(lev, cyc)
  INTEGER, INTENT(IN) :: lev
  INTEGER, INTENT(IN) :: cyc

  INTEGER :: i, n

#ifdef PARALLEL_MPI
  INTEGER :: ierr
#endif

!  IF (hrank==0) CALL indicator(lev,'-')

  IF (lev==nlev) THEN
     IF (use_il) THEN
        CALL apply_il(mg(lev)%n1, mg(lev)%n2, il, mg(lev)%r, mg(lev)%e)
     ELSE
        CALL apply2d1(mg(lev)%n1, mg(lev)%n2, mg(lev)%m, mg(lev)%r, mg(lev)%e)
        CALL sync2d(lev, mg(lev)%e)
     END IF
     RETURN
  END IF

#ifdef PARALLEL_MPI
  IF (lev == glev) THEN
!     CALL apply2d1(mg(lev)%n1, mg(lev)%n2, mg(lev)%m, mg(lev)%r, mg(lev)%e)
!     CALL sync2d(lev, mg(lev)%e)
     CALL gather(mg(lev)%r, mg(lev+1)%r)
     IF (hrank==0) CALL mgsmooth2d(lev+1, cyc)
     CALL mpi_barrier(hcomm, ierr)
     CALL scatter(mg(lev+1)%e, mg(lev)%e)
     RETURN
  END IF
#endif

  IF (cmode == 1  .OR. (cmode==2 .AND. cyc==0)) THEN
     n=1
  ELSE
     n=0
  END IF

  CALL apply2d1(mg(lev)%n1, mg(lev)%n2, mg(lev)%m, mg(lev)%r, mg(lev)%e)

  DO i=0, n
     CALL sync2d(lev, mg(lev)%e)

     CALL apply2d2(mg(lev)%n1, mg(lev)%n2, mg(lev)%l, mg(lev)%r, mg(lev)%e, mg(lev)%t)

     CALL coarse2d(mg(lev)%n1, mg(lev)%n2, mg(lev+1)%n1, mg(lev+1)%n2, mg(lev)%p1, mg(lev)%p2, mg(lev)%q1, mg(lev)%q2, mg(lev)%t, mg(lev+1)%r)

     CALL sync2d(lev+1, mg(lev+1)%r)

     CALL mgsmooth2d(lev+1, cyc+i)

     CALL fine2d(mg(lev)%n1, mg(lev)%n2, mg(lev+1)%n1, mg(lev+1)%n2, mg(lev)%p1, mg(lev)%p2, mg(lev)%q1, mg(lev)%q2, mg(lev+1)%e, mg(lev)%e)

     CALL sync2d(lev, mg(lev)%e)

     CALL apply2d2(mg(lev)%n1, mg(lev)%n2, mg(lev)%l, mg(lev)%r, mg(lev)%e, mg(lev)%t)

     CALL sync2d(lev, mg(lev)%t)

     CALL apply2d3(mg(lev)%n1, mg(lev)%n2, mg(lev)%m, mg(lev)%t, mg(lev)%e)
  END DO

!  IF (hrank==0) CALL indicator(lev,'-')

END SUBROUTINE mgsmooth2d

!-----------------------------------------------------------------------------------------------------------------------

RECURSIVE SUBROUTINE mgsmooth2d_omp(lev, cyc, tid)
  INTEGER, INTENT(IN) :: lev
  INTEGER, INTENT(IN) :: cyc
  INTEGER, INTENT(IN) :: tid

  INTEGER :: i, n
  INTEGER :: jstart, jend

#ifdef PARALLEL_MPI
  INTEGER :: ierr
#endif

  jstart = n2_start(tid,lev)
  jend   = n2_start(tid,lev) + n2_span(tid,lev) - 1

!  IF (hrank==0) CALL indicator(lev,'-')

  IF (lev==nlev) THEN
     IF (use_il) THEN
!$OMP SINGLE
        CALL apply_il(mg(lev)%n1, mg(lev)%n2, il, mg(lev)%r, mg(lev)%e)
!$OMP END SINGLE
     ELSE
!$OMP SINGLE
        CALL apply2d1(mg(lev)%n1, mg(lev)%n2, mg(lev)%m, mg(lev)%r, mg(lev)%e, 1, mg(lev)%n2)
!$OMP END SINGLE
        CALL sync2d(lev, mg(lev)%e, tid)
     END IF
     RETURN
  END IF

#ifdef PARALLEL_MPI
  IF (lev == glev) THEN
!     CALL apply2d1(mg(lev)%n1, mg(lev)%n2, mg(lev)%m, mg(lev)%r, mg(lev)%e)
!     CALL sync2d(lev, mg(lev)%e)
!$OMP SINGLE
     CALL gather(mg(lev)%r, mg(lev+1)%r)
!$OMP END SINGLE
     IF (hrank==0) CALL mgsmooth2d_omp(lev+1, cyc, tid)
!$OMP SINGLE
     CALL mpi_barrier(hcomm, ierr)
     CALL scatter(mg(lev+1)%e, mg(lev)%e)
!$OMP END SINGLE
     RETURN
  END IF
#endif

  IF (cmode == 1  .OR. (cmode==2 .AND. cyc==0)) THEN
     n=1
  ELSE
     n=0
  END IF

  CALL apply2d1(mg(lev)%n1, mg(lev)%n2, mg(lev)%m, mg(lev)%r, mg(lev)%e, jstart, jend)

  DO i=0, n
     CALL sync2d(lev, mg(lev)%e, tid)

     CALL apply2d2(mg(lev)%n1, mg(lev)%n2, mg(lev)%l, mg(lev)%r, mg(lev)%e, mg(lev)%t, jstart, jend)

!$OMP BARRIER
     CALL coarse2d(mg(lev)%n1, mg(lev)%n2, mg(lev+1)%n1, mg(lev+1)%n2, mg(lev)%p1, mg(lev)%p2, &
                   mg(lev)%q1, mg(lev)%q2, mg(lev)%t, mg(lev+1)%r, &
                   n2_start(tid,lev+1), n2_start(tid,lev+1)+n2_span(tid,lev+1)-1)

     CALL sync2d(lev+1, mg(lev+1)%r, tid)

     CALL mgsmooth2d_omp(lev+1, cyc+i, tid)

     CALL fine2d(mg(lev)%n1, mg(lev)%n2, mg(lev+1)%n1, mg(lev+1)%n2, mg(lev)%p1, mg(lev)%p2, mg(lev)%q1, mg(lev)%q2, mg(lev+1)%e, mg(lev)%e, jstart, jend)

     CALL sync2d(lev, mg(lev)%e, tid)

     CALL apply2d2(mg(lev)%n1, mg(lev)%n2, mg(lev)%l, mg(lev)%r, mg(lev)%e, mg(lev)%t, jstart, jend)

     CALL sync2d(lev, mg(lev)%t, tid)

     CALL apply2d3(mg(lev)%n1, mg(lev)%n2, mg(lev)%m, mg(lev)%t, mg(lev)%e, jstart, jend)
  END DO
!$OMP BARRIER

!  IF (hrank==0) CALL indicator(lev,'-')

END SUBROUTINE mgsmooth2d_omp

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE multigrid2d(r, e)
  REAL(8), INTENT(IN)    :: r(0:n1+1, 0:n2+1)
  REAL(8), INTENT(INOUT) :: e(0:n1+1, 0:n2+1)

!  IF (hrank==0) CALL indicator(0,'-')

  CALL apply2d1(n1, n2, mg(0)%m, r, e)

  CALL sync2d(0, e)

  CALL apply2d2(n1, n2, mg(0)%l, r, e, work)

  CALL coarse2d(mg(0)%n1, mg(0)%n2, mg(1)%n1, mg(1)%n2, mg(0)%p1, mg(0)%p2, mg(0)%q1, mg(0)%q2, work, mg(1)%r)

  CALL sync2d(1, mg(1)%r)

  CALL mgsmooth2d(1, 0)

  CALL fine2d(mg(0)%n1, mg(0)%n2, mg(1)%n1, mg(1)%n2, mg(0)%p1, mg(0)%p2, mg(0)%q1, mg(0)%q2, mg(1)%e, e)

  CALL sync2d(0, e)

  CALL apply2d2(n1, n2, mg(0)%l, r, e, work)

  CALL sync2d(0, work)

  CALL apply2d3(n1, n2, mg(0)%m, work, e)

!  IF (hrank==0) CALL indicator(0,'-')

END SUBROUTINE multigrid2d

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE multigrid2d_omp(r, e, tid)
  REAL(8), INTENT(IN)    :: r(0:n1+1, 0:n2+1)
  REAL(8), INTENT(INOUT) :: e(0:n1+1, 0:n2+1)
  INTEGER, INTENT(IN)    :: tid

  INTEGER :: jstart, jend

  jstart = n2_start(tid,0)
  jend   = n2_start(tid,0) + n2_span(tid,0) - 1
!  IF (hrank==0) CALL indicator(0,'-')

  CALL apply2d1(n1, n2, mg(0)%m, r, e, jstart, jend)

  CALL sync2d(0, e, tid)

  CALL apply2d2(n1, n2, mg(0)%l, r, e, work, jstart, jend)
!$OMP BARRIER

  CALL coarse2d(mg(0)%n1, mg(0)%n2, mg(1)%n1, mg(1)%n2, mg(0)%p1, mg(0)%p2, mg(0)%q1, mg(0)%q2, work, mg(1)%r, &
                n2_start(tid,1), n2_start(tid,1)+n2_span(tid,1)-1)

  CALL sync2d(1, mg(1)%r, tid)

  CALL mgsmooth2d_omp(1, 0, tid)

  CALL fine2d(mg(0)%n1, mg(0)%n2, mg(1)%n1, mg(1)%n2, mg(0)%p1, mg(0)%p2, mg(0)%q1, mg(0)%q2, mg(1)%e, e, jstart, jend)

  CALL sync2d(0, e, tid)

  CALL apply2d2(n1, n2, mg(0)%l, r, e, work, jstart, jend)

  CALL sync2d(0, work, tid)

  CALL apply2d3(n1, n2, mg(0)%m, work, e, jstart, jend)

!  IF (hrank==0) CALL indicator(0,'-')

END SUBROUTINE multigrid2d_omp

!-----------------------------------------------------------------------------------------------------------------------

#ifdef PARALLEL_MPI

SUBROUTINE gather_real4(in, out, default)
  REAL(4), INTENT(IN)  :: in( 0:mg(glev  )%n1+1,0:mg(glev  )%n2+1)
  REAL(4), INTENT(OUT) :: out(0:mg(glev+1)%n1+1,0:mg(glev+1)%n2+1)
  REAL(4), INTENT(IN), OPTIONAL :: default

  REAL(4) :: default_
  REAL(4) :: tmp(MG_GDIM,MG_GDIM)
  INTEGER :: i, j, l, m, n
  INTEGER :: ierr

  default_ = 0.0
  IF (present(default)) default_ = default

  tmp(:,:) = in(1:MG_GDIM,1:MG_GDIM)
  CALL mpi_gather(tmp, MG_GDIM**2, MPI_REAL4, gsbuf_r4, MG_GDIM**2, MPI_REAL4, 0, hcomm, ierr)

  IF (hrank /= 0) RETURN

!$OMP PARALLEL PRIVATE(l,m,n)
!$OMP DO COLLAPSE(2)
  DO j=0, jpes-1
  DO i=0, ipes-1
     l = hranks(i,j)
     IF (l == MPI_PROC_NULL) THEN
        out(i*MG_GDIM+1:(i+1)*MG_GDIM,j*MG_GDIM+1:(j+1)*MG_GDIM) = default_
     ELSE
        l = l*MG_GDIM**2
        DO n=j*MG_GDIM+1, (j+1)*MG_GDIM
        DO m=i*MG_GDIM+1, (i+1)*MG_GDIM
           l = l+1
           out(m,n) = gsbuf_r4(l)
        END DO
        END DO
     END IF
  END DO
  END DO

  IF (period1) THEN
!$OMP DO
     DO j=1, mg(glev+1)%n2
        out(0,              j) = out(mg(glev+1)%n1,j)
        out(mg(glev+1)%n1+1,j) = out(1,            j)
     END DO
  END IF

  IF (period2) THEN
!$OMP DO
     DO i=1, mg(glev+1)%n1
        out(i,0)               = out(i,mg(glev+1)%n2)
        out(i,mg(glev+1)%n2+1) = out(i,1)
     END DO
  END IF
!$OMP END PARALLEL

END SUBROUTINE gather_real4

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE gather_real8(in, out, default)
  REAL(8), INTENT(IN)  :: in( 0:mg(glev  )%n1+1,0:mg(glev  )%n2+1)
  REAL(8), INTENT(OUT) :: out(0:mg(glev+1)%n1+1,0:mg(glev+1)%n2+1)
  REAL(8), INTENT(IN), OPTIONAL :: default

  REAL(8) :: default_
  REAL(8) :: tmp(MG_GDIM, MG_GDIM)
  INTEGER :: i, j, l, m, n
  INTEGER :: ierr

  default_ = 0.0
  IF (present(default)) default_ = default

  tmp(:,:) = in(1:MG_GDIM,1:MG_GDIM)
  CALL mpi_gather(tmp, MG_GDIM**2, MPI_REAL8, gsbuf_r8, MG_GDIM**2, MPI_REAL8, 0, hcomm, ierr)

  IF (hrank /= 0) RETURN

!$OMP PARALLEL PRIVATE(l,m,n)
!$OMP DO COLLAPSE(2)
  DO j=0, jpes-1
  DO i=0, ipes-1
     l = hranks(i,j)
     IF (l == MPI_PROC_NULL) THEN
        out(i*MG_GDIM+1:(i+1)*MG_GDIM,j*MG_GDIM+1:(j+1)*MG_GDIM) = default_
     ELSE
        l = l*MG_GDIM**2
        DO n=j*MG_GDIM+1, (j+1)*MG_GDIM
        DO m=i*MG_GDIM+1, (i+1)*MG_GDIM
           l = l+1
           out(m,n) = gsbuf_r8(l)
        END DO
        END DO
     END IF
  END DO
  END DO

  IF (period1) THEN
!$OMP DO
     DO j=1, mg(glev+1)%n2
        out(0,              j) = out(mg(glev+1)%n1,j)
        out(mg(glev+1)%n1+1,j) = out(1,            j)
     END DO
  END IF

  IF (period2) THEN
!$OMP DO
     DO i=1, mg(glev+1)%n1
        out(i,0)               = out(i,mg(glev+1)%n2)
        out(i,mg(glev+1)%n2+1) = out(i,1)
     END DO
  END IF
!$OMP END PARALLEL

END SUBROUTINE gather_real8

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE scatter_real4(in, out)
  REAL(4), INTENT(IN)    :: in( 0:mg(glev+1)%n1+1,0:mg(glev+1)%n2+1)
  REAL(4), INTENT(INOUT) :: out(0:mg(glev  )%n1+1,0:mg(glev  )%n2+1)

  REAL(4) :: tmp(MG_GDIM,MG_GDIM)
  INTEGER :: i, j, l, m, n
  INTEGER :: ierr

  IF (hrank==0) THEN
!$OMP PARALLEL DO COLLAPSE(2) PRIVATE(l,m,n)
     DO j=0, jpes-1
     DO i=0, ipes-1
        l = hranks(i,j)
        IF (l == MPI_PROC_NULL) CYCLE
        l = l*MG_GDIM**2
        DO n=j*MG_GDIM+1, (j+1)*MG_GDIM
        DO m=i*MG_GDIM+1, (i+1)*MG_GDIM
           l = l+1
           gsbuf_r4(l) = in(m,n)
        END DO
        END DO
     END DO
     END DO
  END IF

  CALL mpi_scatter(gsbuf_r4, MG_GDIM**2, MPI_REAL4, tmp, MG_GDIM**2, MPI_REAL4, 0, hcomm, ierr)
  out(1:MG_GDIM,1:MG_GDIM) = tmp(:,:)

END SUBROUTINE scatter_real4

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE scatter_real8(in, out)
  REAL(8), INTENT(IN)    :: in( 0:mg(glev+1)%n1+1,0:mg(glev+1)%n2+1)
  REAL(8), INTENT(INOUT) :: out(0:mg(glev  )%n1+1,0:mg(glev  )%n2+1)

  REAL(8) :: tmp(MG_GDIM,MG_GDIM)
  INTEGER :: i, j, l, m, n
  INTEGER :: ierr

  IF (hrank==0) THEN
!$OMP PARALLEL DO COLLAPSE(2) PRIVATE(l,m,n)
     DO j=0, jpes-1
     DO i=0, ipes-1
        l = hranks(i,j)
        IF (l == MPI_PROC_NULL) CYCLE
        l = l*MG_GDIM**2
        DO n=j*MG_GDIM+1, (j+1)*MG_GDIM
        DO m=i*MG_GDIM+1, (i+1)*MG_GDIM
           l = l+1
           gsbuf_r8(l) = in(m,n)
        END DO
        END DO
     END DO
     END DO
  END IF

  CALL mpi_scatter(gsbuf_r8, MG_GDIM**2, MPI_REAL8, tmp, MG_GDIM**2, MPI_REAL8, 0, hcomm, ierr)
  out(1:MG_GDIM,1:MG_GDIM) = tmp(:,:)

END SUBROUTINE scatter_real8

!-----------------------------------------------------------------------------------------------------------------------
#endif

SUBROUTINE sync2d_begin_real4(lev, buf, tid)
  INTEGER, INTENT(IN)    :: lev
  REAL(4), INTENT(INOUT) :: buf(0:mg(lev)%n1+1, 0:mg(lev)%n2+1)
  INTEGER, INTENT(IN), OPTIONAL :: tid

  INTEGER :: n1, n2
  INTEGER :: j, jstart, jend

#ifdef PARALLEL_MPI
  INTEGER :: ierr
#endif

  n1 = mg(lev)%n1
  n2 = mg(lev)%n2

  IF (present(tid)) THEN
     jstart = n2_start(tid,lev)
     jend   = n2_start(tid,lev) + n2_span(tid,lev)-1
  ELSE
     jstart = 1
     jend   = n2
  END IF
!$OMP BARRIER

#ifdef PARALLEL_MPI
  IF (lev > glev) RETURN

!$OMP SINGLE
  CALL mpi_startall(4, req2d_r4(0:3,lev), ierr)
!$OMP END SINGLE nowait

  DO j=jstart, jend
     sendbuf_e_r4(j) = buf(n1,j)
     sendbuf_w_r4(j) = buf(1, j)
  END DO

!$OMP SECTIONS
  sendbuf_n_r4(1:n1) = buf(1:n1,n2)
!$OMP SECTION
  sendbuf_s_r4(1:n1) = buf(1:n1,1)
!$OMP END SECTIONS

!$OMP SINGLE
  CALL mpi_startall(4, req2d_r4(4:7,lev), ierr)
!$OMP END SINGLE
#endif

END SUBROUTINE sync2d_begin_real4

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE sync2d_begin_real8(lev, buf, tid)
  INTEGER, INTENT(IN)    :: lev
  REAL(8), INTENT(INOUT) :: buf(0:mg(lev)%n1+1, 0:mg(lev)%n2+1)
  INTEGER, INTENT(IN), OPTIONAL :: tid

  INTEGER :: n1, n2
  INTEGER :: j, jstart, jend

#ifdef PARALLEL_MPI
  INTEGER :: ierr
#endif

  n1 = mg(lev)%n1
  n2 = mg(lev)%n2

  IF (present(tid)) THEN
     jstart = n2_start(tid,lev)
     jend   = n2_start(tid,lev) + n2_span(tid,lev)-1
  ELSE
     jstart = 1
     jend   = n2
  END IF
!$OMP BARRIER

#ifdef PARALLEL_MPI
  IF (lev > glev) RETURN

!$OMP SINGLE
  CALL mpi_startall(4, req2d_r8(0:3,lev), ierr)
!$OMP END SINGLE nowait

  DO j=jstart, jend
     sendbuf_e_r8(j) = buf(n1,j)
     sendbuf_w_r8(j) = buf(1, j)
  END DO

!$OMP SECTIONS
  sendbuf_n_r8(1:n1) = buf(1:n1,n2)
!$OMP SECTION
  sendbuf_s_r8(1:n1) = buf(1:n1,1)
!$OMP END SECTIONS

!$OMP SINGLE
  CALL mpi_startall(4, req2d_r8(4:7,lev), ierr)
!$OMP END SINGLE
#endif

END SUBROUTINE sync2d_begin_real8

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE sync2d_end_real4(lev, buf, tid)
  INTEGER, INTENT(IN)    :: lev
  REAL(4), INTENT(INOUT) :: buf(0:mg(lev)%n1+1, 0:mg(lev)%n2+1)
  INTEGER, INTENT(IN), OPTIONAL :: tid

  INTEGER :: n1, n2
  INTEGER :: i, j
  INTEGER :: jstart, jend

#ifdef PARALLEL_MPI
  INTEGER :: ierr
#endif

  n1 = mg(lev)%n1
  n2 = mg(lev)%n2

  IF (present(tid)) THEN
     jstart = n2_start(tid,lev)
     jend   = n2_start(tid,lev) + n2_span(tid,lev)-1
  ELSE
     jstart = 1
     jend   = n2
  END IF

#ifdef PARALLEL_MPI
  IF (lev > glev) THEN
     IF (period1) THEN
        DO j=jstart, jend
           buf(0,j)    = buf(n1,j)
           buf(n1+1,j) = buf(1,j)
        END DO
     ELSE
        DO j=jstart, jend
           buf(0,j)    = 0.0
           buf(n1+1,j) = 0.0
        END DO
     END IF

     IF (period2) THEN
!$OMP SECTIONS
        buf(1:n1,0)    = buf(1:n1,n2)
!$OMP SECTION
        buf(1:n1,n2+1) = buf(1:n1,1)
!$OMP END SECTIONS
     ELSE
!$OMP SECTIONS
        buf(1:n1,0)    = 0.0
!$OMP SECTION
        buf(1:n1,n2+1) = 0.0
!$OMP END SECTIONS
     END IF
  ELSE

!$OMP SINGLE
     CALL mpi_waitall(8, req2d_r4(0:7,lev), MPI_STATUSES_IGNORE, ierr)
!$OMP END SINGLE

!$OMP SECTIONS
     IF (mg(lev)%rank_e /= MPI_PROC_NULL) buf(n1+1,1:n2) = recvbuf_e_r4(1:n2)
!$OMP SECTION
     IF (mg(lev)%rank_w /= MPI_PROC_NULL) buf(0,   1:n2) = recvbuf_w_r4(1:n2)
!$OMP SECTION
     IF (mg(lev)%rank_n /= MPI_PROC_NULL) buf(1:n1,n2+1) = recvbuf_n_r4(1:n1)
!$OMP SECTION
     IF (mg(lev)%rank_s /= MPI_PROC_NULL) buf(1:n1,0)    = recvbuf_s_r4(1:n1)
!$OMP END SECTIONS
  END IF
#else
  IF (period1) THEN
     DO j=jstart, jend
        buf(0,   j) = buf(n1,j)
        buf(n1+1,j) = buf(1, j)
     END DO
  END IF

  IF (period2) THEN
!$OMP SECTIONS
     buf(1:n1,0)    = buf(1:n1,n2)
!$OMP SECTION
     buf(1:n1,n2+1) = buf(1:n1,1)
!$OMP END SECTIONS
  END IF
#endif
!$OMP BARRIER
END SUBROUTINE sync2d_end_real4

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE sync2d_end_real8(lev, buf, tid)
  INTEGER, INTENT(IN)    :: lev
  REAL(8), INTENT(INOUT) :: buf(0:mg(lev)%n1+1, 0:mg(lev)%n2+1)
  INTEGER, INTENT(IN), OPTIONAL :: tid

  INTEGER :: n1, n2
  INTEGER :: i, j
  INTEGER :: jstart, jend

#ifdef PARALLEL_MPI
  INTEGER :: ierr
#endif

  n1 = mg(lev)%n1
  n2 = mg(lev)%n2

  IF (present(tid)) THEN
     jstart = n2_start(tid,lev)
     jend   = n2_start(tid,lev) + n2_span(tid,lev)-1
  ELSE
     jstart = 1
     jend   = n2
  END IF

#ifdef PARALLEL_MPI
  IF (lev > glev) THEN
     IF (period1) THEN
        DO j=jstart, jend
           buf(0,   j) = buf(n1,j)
           buf(n1+1,j) = buf(1,j)
        END DO
     ELSE
        DO j=jstart, jend
           buf(0,j)            = 0.0D0
           buf(n1+1,j) = 0.0D0
        END DO
     END IF

     IF (period2) THEN
!$OMP SECTIONS
        buf(1:n1,0)    = buf(1:n1,n2)
!$OMP SECTION
        buf(1:n1,n2+1) = buf(1:n1,1)
!$OMP END SECTIONS
     ELSE
!$OMP SECTIONS
        buf(1:n1,0)    = 0.0D0
!$OMP SECTION
        buf(1:n1,n2+1) = 0.0D0
!$OMP END SECTIONS
     END IF
  ELSE

!$OMP SINGLE
     CALL mpi_waitall(8, req2d_r8(0:7,lev), MPI_STATUSES_IGNORE, ierr)
!$OMP END SINGLE

!$OMP SECTIONS
     IF (mg(lev)%rank_e /= MPI_PROC_NULL) buf(n1+1,1:n2) = recvbuf_e_r8(1:n2)
!$OMP SECTION
     IF (mg(lev)%rank_w /= MPI_PROC_NULL) buf(0,   1:n2) = recvbuf_w_r8(1:n2)
!$OMP SECTION
     IF (mg(lev)%rank_n /= MPI_PROC_NULL) buf(1:n1,n2+1) = recvbuf_n_r8(1:n1)
!$OMP SECTION
     IF (mg(lev)%rank_s /= MPI_PROC_NULL) buf(1:n1,0)    = recvbuf_s_r8(1:n1)
!$OMP END SECTIONS
  END IF
#else
  IF (period1) THEN
     DO j=jstart, jend
        buf(0,j)    = buf(n1,j)
        buf(n1+1,j) = buf(1, j)
     END DO
  END IF

  IF (period2) THEN
!$OMP SECTIONS
     buf(1:n1,0)    = buf(1:n1,n2)
!$OMP SECTION
     buf(1:n1,n2+1) = buf(1:n1,1)
!$OMP END SECTIONS
  END IF
#endif
!$OMP BARRIER
END SUBROUTINE sync2d_end_real8

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE sync2d_real4(lev, buf, tid)
  INTEGER, INTENT(IN)    :: lev
  REAL(4), INTENT(INOUT) :: buf(0:mg(lev)%n1+1, 0:mg(lev)%n2+1)
  INTEGER, INTENT(IN), OPTIONAL :: tid

  IF (present(tid)) THEN
     CALL sync2d_begin(lev, buf, tid)
     CALL sync2d_end(lev, buf, tid)
  ELSE
     CALL sync2d_begin(lev, buf)
     CALL sync2d_end(lev, buf)
  END IF

END SUBROUTINE sync2d_real4

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE sync2d_real8(lev, buf, tid)
  INTEGER, INTENT(IN)    :: lev
  REAL(8), INTENT(INOUT) :: buf(0:mg(lev)%n1+1, 0:mg(lev)%n2+1)
  INTEGER, INTENT(IN), OPTIONAL :: tid

  IF (present(tid)) THEN
     CALL sync2d_begin(lev, buf, tid)
     CALL sync2d_end(lev, buf, tid)
  ELSE
     CALL sync2d_begin(lev, buf)
     CALL sync2d_end(lev, buf)
  END IF

END SUBROUTINE sync2d_real8

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE get_coefficient(row, data)
  INTEGER, INTENT(IN) :: row
  REAL(8), INTENT(OUT) :: data(:,:)

  INTEGER :: i, j, k

  CALL assert(initialized, "solver2d is not initialized.")
  CALL assert(row >= -3 .AND. row <= 3, "must be -3 <= row <= 3")

  DO j=1, n2
  DO i=1, n1
     data(i,j) = l2d(row,i,j)
  END DO
  END DO

END SUBROUTINE get_coefficient

!-----------------------------------------------------------------------------------------------------------------------

END MODULE solver2d
