#include "macro.h"

MODULE solver3d
  USE misc
  USE parallel
  USE solver_util
#ifdef PROFILE
  USE profile
#endif
  IMPLICIT NONE
  SAVE
  PRIVATE
  PUBLIC init_solver3d, solve3d, reset_solver3d, finalize_solver3d

  LOGICAL :: initialized = .FALSE.

  INTEGER :: solver_method = 2 ! 0 : pure CG,
                               ! 1 : Stand-alone MG
                               ! 2 : MGCG
                               ! 3 : MGCGS

  INTEGER :: n1
  INTEGER :: n2
  INTEGER :: n3

  REAL(8), ALLOCATABLE :: l3d(:,:,:,:)

  TYPE :: mg_struct
     INTEGER :: n1 = 1
     INTEGER :: n2 = 1
     INTEGER :: n3 = 1

#ifdef NO_ALLOCATABLE_IN_TYPE
     REAL(MG_KIND), POINTER :: r(:,:,:)
     REAL(MG_KIND), POINTER :: e(:,:,:)
     REAL(MG_KIND), POINTER :: t(:,:,:)

     REAL(MG_KIND), POINTER :: l(:,:,:,:)
     REAL(MG_KIND), POINTER :: m(:,:,:,:)

     INTEGER, POINTER :: p1(:)
     INTEGER, POINTER :: p2(:)
     INTEGER, POINTER :: p3(:)

     INTEGER, POINTER :: q1(:)
     INTEGER, POINTER :: q2(:)
     INTEGER, POINTER :: q3(:)

     INTEGER, POINTER :: lindx(:)
     INTEGER, POINTER :: mindx(:)
#else
     REAL(MG_KIND), ALLOCATABLE :: r(:,:,:)
     REAL(MG_KIND), ALLOCATABLE :: e(:,:,:)
     REAL(MG_KIND), ALLOCATABLE :: t(:,:,:)

     REAL(MG_KIND), ALLOCATABLE :: l(:,:,:,:)
     REAL(MG_KIND), ALLOCATABLE :: m(:,:,:,:)

     INTEGER, ALLOCATABLE :: p1(:)
     INTEGER, ALLOCATABLE :: p2(:)
     INTEGER, ALLOCATABLE :: p3(:)

     INTEGER, ALLOCATABLE :: q1(:)
     INTEGER, ALLOCATABLE :: q2(:)
     INTEGER, ALLOCATABLE :: q3(:)

     INTEGER, ALLOCATABLE :: lindx(:)
     INTEGER, ALLOCATABLE :: mindx(:)
#endif

#ifdef PARALLEL_MPI
     INTEGER :: rank_e
     INTEGER :: rank_w
     INTEGER :: rank_n
     INTEGER :: rank_s
     INTEGER :: rank_u
     INTEGER :: rank_l
#endif
  END TYPE mg_struct

  TYPE(mg_struct) :: mg(0:max_level)

  REAL(8) :: eps      = 1.0D-8
  INTEGER :: cmode    = 0
  INTEGER :: itmax    = 1000
  LOGICAL :: itreport = .FALSE.
  LOGICAL :: fallback = .FALSE.
  REAL(8) :: fallback_eps = 1.0D-16
  INTEGER :: coeffscheme  = 1

  LOGICAL :: open_e = .FALSE.
  LOGICAL :: open_w = .FALSE.
  LOGICAL :: open_n = .FALSE.
  LOGICAL :: open_s = .FALSE.

  LOGICAL :: period1 = .FALSE.
  LOGICAL :: period2 = .FALSE.
  LOGICAL :: period3 = .FALSE.

  LOGICAL :: use_il = .FALSE.
  LOGICAL :: singular

  INTEGER :: nlev = max_level

  INTEGER :: vlev = max_level
  INTEGER :: glev = max_level

#ifdef PARALLEL_MPI
  REAL(4), ALLOCATABLE :: gsbuf_r4(:)
  REAL(8), ALLOCATABLE :: gsbuf_r8(:)

  REAL(4), ALLOCATABLE :: sendbuf_e_r4(:)
  REAL(4), ALLOCATABLE :: sendbuf_w_r4(:)
  REAL(4), ALLOCATABLE :: sendbuf_n_r4(:)
  REAL(4), ALLOCATABLE :: sendbuf_s_r4(:)
  REAL(4), ALLOCATABLE :: sendbuf_u_r4(:)
  REAL(4), ALLOCATABLE :: sendbuf_l_r4(:)

  REAL(8), ALLOCATABLE :: sendbuf_e_r8(:)
  REAL(8), ALLOCATABLE :: sendbuf_w_r8(:)
  REAL(8), ALLOCATABLE :: sendbuf_n_r8(:)
  REAL(8), ALLOCATABLE :: sendbuf_s_r8(:)
  REAL(8), ALLOCATABLE :: sendbuf_u_r8(:)
  REAL(8), ALLOCATABLE :: sendbuf_l_r8(:)

  REAL(4), ALLOCATABLE :: recvbuf_e_r4(:)
  REAL(4), ALLOCATABLE :: recvbuf_w_r4(:)
  REAL(4), ALLOCATABLE :: recvbuf_n_r4(:)
  REAL(4), ALLOCATABLE :: recvbuf_s_r4(:)
  REAL(4), ALLOCATABLE :: recvbuf_u_r4(:)
  REAL(4), ALLOCATABLE :: recvbuf_l_r4(:)

  REAL(8), ALLOCATABLE :: recvbuf_e_r8(:)
  REAL(8), ALLOCATABLE :: recvbuf_w_r8(:)
  REAL(8), ALLOCATABLE :: recvbuf_n_r8(:)
  REAL(8), ALLOCATABLE :: recvbuf_s_r8(:)
  REAL(8), ALLOCATABLE :: recvbuf_u_r8(:)
  REAL(8), ALLOCATABLE :: recvbuf_l_r8(:)

  INTEGER :: req2d_r4(0:7,  0:max_level)
  INTEGER :: req2d_r8(0:7,  0:max_level)
  INTEGER :: req3d_r4(0:11, 0:max_level)
  INTEGER :: req3d_r8(0:11, 0:max_level)
#else
  INTEGER :: comm  = 0 !dummy
#endif

  REAL(MG_KIND), ALLOCATABLE :: il(:,:)
  REAL(8),       ALLOCATABLE :: work(:,:,:)


  INTEGER :: n_call = 0
  INTEGER :: n_iter = 0

  INTEGER :: num_threads = 1

  INTEGER :: n3_start(0:max_threads-1, 0:max_level) = 0
  INTEGER :: n3_span( 0:max_threads-1, 0:max_level) = 0

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

  INTERFACE sync3d
     MODULE PROCEDURE sync3d_real4
     MODULE PROCEDURE sync3d_real8
  END INTERFACE sync3d

  INTERFACE sync3d_begin
     MODULE PROCEDURE sync3d_begin_real4
     MODULE PROCEDURE sync3d_begin_real8
  END INTERFACE sync3d_begin

  INTERFACE sync3d_end
     MODULE PROCEDURE sync3d_end_real4
     MODULE PROCEDURE sync3d_end_real8
  END INTERFACE sync3d_end

#ifdef PARALLEL_MPI
  INTERFACE gather_h
     MODULE PROCEDURE gather_h_real4
     MODULE PROCEDURE gather_h_real8
  END INTERFACE gather_h

  INTERFACE gather_v
     MODULE PROCEDURE gather_v_real4
     MODULE PROCEDURE gather_v_real8
  END INTERFACE gather_v

  INTERFACE scatter_h
     MODULE PROCEDURE scatter_h_real4
     MODULE PROCEDURE scatter_h_real8
  END INTERFACE scatter_h

  INTERFACE scatter_v
     MODULE PROCEDURE scatter_v_real4
     MODULE PROCEDURE scatter_v_real8
  END INTERFACE scatter_v
#endif
CONTAINS

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE init_common(periods, opens)
  LOGICAL, INTENT(IN), OPTIONAL :: periods(3)
  LOGICAL, INTENT(IN), OPTIONAL :: opens(4)

  INTEGER :: n
  LOGICAL :: tmp

#ifdef PARALLEL_MPI
  INTEGER :: ierr
#endif

  CALL assert(.NOT. initialized, "SOLVER3D is already initialized!")

  initialized = .TRUE.

  IF (present(periods)) THEN
     period1 = periods(1)
     period2 = periods(2)
     period3 = periods(3)
  END IF

#ifdef PARALLEL_MPI
  ALLOCATE(gsbuf_r4(npes*MG_GDIM**2))
  ALLOCATE(gsbuf_r8(npes*MG_GDIM**2))

  gsbuf_r4(:) = UNDEF
  gsbuf_r8(:) = UNDEF
#endif

  IF (rank==0) WRITE(REPORT_UNIT, *) 'initialize solver3d'

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

SUBROUTINE init_solver3d(nx, ny, nz,  &
                         dx,   dy,   dz,   dsx,   dsy,   dsz,   dvol,   alpha,   mask,   &
                         dx_h, dy_h, dz_h, dsx_h, dsy_h, dsz_h, dvol_h, alpha_h, mask_h, &
                         dx_v, dy_v, dz_v, dsx_v, dsy_v, dsz_v, dvol_v, &
                         periods, opens, method, epsilon, c_mode, it_max, report, scheme)
  INTEGER, INTENT(IN) :: nx
  INTEGER, INTENT(IN) :: ny
  INTEGER, INTENT(IN) :: nz

  REAL(8), INTENT(IN), OPTIONAL  :: dx(:,:,:), dx_h(:,:), dx_v(:)
  REAL(8), INTENT(IN), OPTIONAL  :: dy(:,:,:), dy_h(:,:), dy_v(:)
  REAL(8), INTENT(IN), OPTIONAL  :: dz(:,:,:), dz_h(:,:), dz_v(:)

  REAL(8), INTENT(IN), OPTIONAL  :: dsx(:,:,:), dsx_h(:,:), dsx_v(:)
  REAL(8), INTENT(IN), OPTIONAL  :: dsy(:,:,:), dsy_h(:,:), dsy_v(:)
  REAL(8), INTENT(IN), OPTIONAL  :: dsz(:,:,:), dsz_h(:,:), dsz_v(:)

  REAL(8), INTENT(IN), OPTIONAL  :: dvol(:,:,:), dvol_h(:,:), dvol_v(:)

  REAL(8), INTENT(IN), OPTIONAL  :: alpha(:,:,:), alpha_h(:,:)

  LOGICAL(1), INTENT(IN), OPTIONAL  :: mask(:,:,:), mask_h(:,:)

  LOGICAL, INTENT(IN), OPTIONAL :: periods(3)

  LOGICAL, INTENT(IN), OPTIONAL :: opens(4)

  CHARACTER(*), INTENT(IN), OPTIONAL :: method  ! default: 'MGCG'
  REAL(8),      INTENT(IN), OPTIONAL :: epsilon ! default: 1.0D-8
  INTEGER,      INTENT(IN), OPTIONAL :: c_mode  ! 0 for V-cycle, 1 for W-cycle, 2 for F-cycle
  INTEGER,      INTENT(IN), OPTIONAL :: it_max  ! default: 1000
  LOGICAL,      INTENT(IN), OPTIONAL :: report  ! default: .FALSE.
  INTEGER,      INTENT(IN), OPTIONAL :: scheme  ! default: 1

  REAL(8), ALLOCATABLE :: d1(:,:,:)
  REAL(8), ALLOCATABLE :: d2(:,:,:)
  REAL(8), ALLOCATABLE :: d3(:,:,:)

  REAL(8), ALLOCATABLE :: ds1(:,:,:)
  REAL(8), ALLOCATABLE :: ds2(:,:,:)
  REAL(8), ALLOCATABLE :: ds3(:,:,:)

  REAL(8), ALLOCATABLE :: dv(:,:,:)

  REAL(8), ALLOCATABLE :: da(:,:,:)

  REAL(8), ALLOCATABLE :: dm(:,:,:)

  REAL(8), ALLOCATABLE :: dtmp(:,:,:)
  REAL(8), ALLOCATABLE :: ltmp(:,:,:,:)
  REAL(8), ALLOCATABLE :: mtmp(:,:,:,:)

  LOGICAL, ALLOCATABLE :: lflag(:)
  LOGICAL, ALLOCATABLE :: mflag(:)


  REAL(8) :: d1_, d2_, d3_
  REAL(8) :: dmax

  REAL(8) :: d3_mean(1:nz+1)
  REAL(8) :: d3_tmp

  INTEGER :: lev

  INTEGER :: i, j, k, l, n

  INTEGER :: n3_max

#ifdef PARALLEL_MPI
  INTEGER :: ierr
#endif

  REAL(8) :: tmp

  CALL init_common(periods, opens)

  n1 = nx
  n2 = ny
  n3 = nz

  mg(0)%n1 = n1
  mg(0)%n2 = n2
  mg(0)%n3 = n3


#ifdef PARALLEL_MPI
  ALLOCATE(d1(0:max(n1,ipes*MG_GDIM)+1, 0:max(n2,jpes*MG_GDIM)+1, 0:max(n3,kpes)+1))
  ALLOCATE(d2(0:max(n1,ipes*MG_GDIM)+1, 0:max(n2,jpes*MG_GDIM)+1, 0:max(n3,kpes)+1))
  ALLOCATE(d3(0:max(n1,ipes*MG_GDIM)+1, 0:max(n2,jpes*MG_GDIM)+1, 0:max(n3,kpes)+1))

  ALLOCATE(ds1(0:max(n1,ipes*MG_GDIM)+1, 0:max(n2,jpes*MG_GDIM)+1, 0:max(n3,kpes)+1))
  ALLOCATE(ds2(0:max(n1,ipes*MG_GDIM)+1, 0:max(n2,jpes*MG_GDIM)+1, 0:max(n3,kpes)+1))
  ALLOCATE(ds3(0:max(n1,ipes*MG_GDIM)+1, 0:max(n2,jpes*MG_GDIM)+1, 0:max(n3,kpes)+1))

  ALLOCATE(dv(0:max(n1,ipes*MG_GDIM)+1, 0:max(n2,jpes*MG_GDIM)+1, 0:max(n3,kpes)+1))
  ALLOCATE(da(0:max(n1,ipes*MG_GDIM)+1, 0:max(n2,jpes*MG_GDIM)+1, 0:max(n3,kpes)+1))

  ALLOCATE(dm(0:max(n1,ipes*MG_GDIM)+1, 0:max(n2,jpes*MG_GDIM)+1, 0:max(n3,kpes)+1))
#else
  ALLOCATE(d1(0:n1+1, 0:n2+1, 0:n3+1))
  ALLOCATE(d2(0:n1+1, 0:n2+1, 0:n3+1))
  ALLOCATE(d3(0:n1+1, 0:n2+1, 0:n3+1))

  ALLOCATE(ds1(0:n1+1, 0:n2+1, 0:n3+1))
  ALLOCATE(ds2(0:n1+1, 0:n2+1, 0:n3+1))
  ALLOCATE(ds3(0:n1+1, 0:n2+1, 0:n3+1))

  ALLOCATE(dv(0:n1+1, 0:n2+1, 0:n3+1))
  ALLOCATE(da(0:n1+1, 0:n2+1, 0:n3+1))

  ALLOCATE(dm(0:n1+1, 0:n2+1, 0:n3+1))
#endif

  ALLOCATE(lflag(1:n3))
  ALLOCATE(mflag(1:n3))

  d1(:,:,:) = UNDEF
  d2(:,:,:) = UNDEF
  d3(:,:,:) = UNDEF

#ifdef PARALLEL_MPI
  ALLOCATE(sendbuf_e_r4(n2*n3))
  ALLOCATE(sendbuf_w_r4(n2*n3))
  ALLOCATE(sendbuf_n_r4(n1*n3))
  ALLOCATE(sendbuf_s_r4(n1*n3))
  ALLOCATE(sendbuf_u_r4(n1*n2))
  ALLOCATE(sendbuf_l_r4(n1*n2))

  ALLOCATE(sendbuf_e_r8(n2*n3))
  ALLOCATE(sendbuf_w_r8(n2*n3))
  ALLOCATE(sendbuf_n_r8(n1*n3))
  ALLOCATE(sendbuf_s_r8(n1*n3))
  ALLOCATE(sendbuf_u_r8(n1*n2))
  ALLOCATE(sendbuf_l_r8(n1*n2))

  sendbuf_e_r4 = 0.0
  sendbuf_w_r4 = 0.0
  sendbuf_n_r4 = 0.0
  sendbuf_s_r4 = 0.0
  sendbuf_u_r4 = 0.0
  sendbuf_l_r4 = 0.0

  sendbuf_e_r8 = 0.0D0
  sendbuf_w_r8 = 0.0D0
  sendbuf_n_r8 = 0.0D0
  sendbuf_s_r8 = 0.0D0
  sendbuf_u_r8 = 0.0D0
  sendbuf_l_r8 = 0.0D0

  ALLOCATE(recvbuf_e_r4(n2*n3))
  ALLOCATE(recvbuf_w_r4(n2*n3))
  ALLOCATE(recvbuf_n_r4(n1*n3))
  ALLOCATE(recvbuf_s_r4(n1*n3))
  ALLOCATE(recvbuf_u_r4(n1*n2))
  ALLOCATE(recvbuf_l_r4(n1*n2))

  ALLOCATE(recvbuf_e_r8(n2*n3))
  ALLOCATE(recvbuf_w_r8(n2*n3))
  ALLOCATE(recvbuf_n_r8(n1*n3))
  ALLOCATE(recvbuf_s_r8(n1*n3))
  ALLOCATE(recvbuf_u_r8(n1*n2))
  ALLOCATE(recvbuf_l_r8(n1*n2))

  recvbuf_e_r4 = 0.0
  recvbuf_w_r4 = 0.0
  recvbuf_n_r4 = 0.0
  recvbuf_s_r4 = 0.0
  recvbuf_u_r4 = 0.0
  recvbuf_l_r4 = 0.0

  recvbuf_e_r8 = 0.0D0
  recvbuf_w_r8 = 0.0D0
  recvbuf_n_r8 = 0.0D0
  recvbuf_s_r8 = 0.0D0
  recvbuf_u_r8 = 0.0D0
  recvbuf_l_r8 = 0.0D0
#endif

  ALLOCATE(work(0:n1+1, 0:n2+1, 0:n3+1))
  work = 0.0D0

  DO lev=0, max_level

     IF (lev==0) THEN
        CALL setparam3d(n1, n2, n3, d1, dx, dx_h, dx_v)
        CALL setparam3d(n1, n2, n3, d2, dy, dy_h, dy_v)
        CALL setparam3d(n1, n2, n3, d3, dz, dz_h, dz_v)

        CALL setparam3d(n1, n2, n3, ds1, dsx, dsx_h, dsx_v)
        CALL setparam3d(n1, n2, n3, ds2, dsy, dsy_h, dsy_v)
        CALL setparam3d(n1, n2, n3, ds3, dsz, dsz_h, dsz_v)

        CALL setparam3d(n1, n2, n3, dv, dvol, dvol_h, dvol_v)
        CALL setparam3d(n1, n2, n3, da, alpha, alpha_h, d=0.0D0)

        CALL setparam3d_logical(n1, n2, n3, dm, mask, mask_h, d=LOGICAL(.TRUE.,1))

#ifdef PARALLEL_MPI
     ELSE IF (lev == vlev+1) THEN
        mg(lev)%n1 = mg(lev-1)%n1
        mg(lev)%n2 = mg(lev-1)%n2
        mg(lev)%n3 = kpes

        CALL gather_v( d1(0:mg(lev)%n1+1,0:mg(lev)%n2+1,1),  d1(0:mg(lev)%n1+1,0:mg(lev)%n2+1,0:mg(lev)%n3+1), default=1.0D0)
        CALL gather_v( d2(0:mg(lev)%n1+1,0:mg(lev)%n2+1,1),  d2(0:mg(lev)%n1+1,0:mg(lev)%n2+1,0:mg(lev)%n3+1), default=1.0D0)
        CALL gather_v( d3(0:mg(lev)%n1+1,0:mg(lev)%n2+1,1),  d3(0:mg(lev)%n1+1,0:mg(lev)%n2+1,0:mg(lev)%n3+1), default=1.0D0)

        CALL gather_v(ds1(0:mg(lev)%n1+1,0:mg(lev)%n2+1,1), ds1(0:mg(lev)%n1+1,0:mg(lev)%n2+1,0:mg(lev)%n3+1), default=1.0D0)
        CALL gather_v(ds2(0:mg(lev)%n1+1,0:mg(lev)%n2+1,1), ds2(0:mg(lev)%n1+1,0:mg(lev)%n2+1,0:mg(lev)%n3+1), default=1.0D0)
        CALL gather_v(ds3(0:mg(lev)%n1+1,0:mg(lev)%n2+1,1), ds3(0:mg(lev)%n1+1,0:mg(lev)%n2+1,0:mg(lev)%n3+1), default=1.0D0)

        CALL gather_v( dv(0:mg(lev)%n1+1,0:mg(lev)%n2+1,1),  dv(0:mg(lev)%n1+1,0:mg(lev)%n2+1,0:mg(lev)%n3+1))
        CALL gather_v( da(0:mg(lev)%n1+1,0:mg(lev)%n2+1,1),  da(0:mg(lev)%n1+1,0:mg(lev)%n2+1,0:mg(lev)%n3+1))

        CALL gather_v( dm(0:mg(lev)%n1+1,0:mg(lev)%n2+1,1),  dm(0:mg(lev)%n1+1,0:mg(lev)%n2+1,0:mg(lev)%n3+1))

        IF (vrank==0) THEN
           ds1(0,1:mg(lev)%n2,1:mg(lev)%n3) = 0.0D0
           ds2(1:mg(lev)%n1,0,1:mg(lev)%n3) = 0.0D0
           ds3(1:mg(lev)%n1,1:mg(lev)%n2,0) = 0.0D0

           dm(0,           1:mg(lev)%n2,1:mg(lev)%n3) = 0.0D0
           dm(mg(lev)%n1+1,1:mg(lev)%n2,1:mg(lev)%n3) = 0.0D0

           dm(1:mg(lev)%n1,0,           1:mg(lev)%n3) = 0.0D0
           dm(1:mg(lev)%n1,mg(lev)%n2+1,1:mg(lev)%n3) = 0.0D0

           dm(1:mg(lev)%n1,1:mg(lev)%n2,0)            = 0.0D0
           dm(1:mg(lev)%n1,1:mg(lev)%n2,mg(lev)%n3+1) = 0.0D0
        END IF

!     ELSE IF (lev > vlev+1 .AND. lev == glev + 1) THEN
     ELSE IF ((lev > vlev+1 .OR. kpes==1) .AND. lev == glev + 1) THEN
        CALL assert(mg(lev-1)%n3 == 1, "INIT_SOLVER3D: MG(GLEV)%DIME3 /= 1")

        mg(lev)%n1 = mg(lev-1)%n1*ipes
        mg(lev)%n2 = mg(lev-1)%n2*jpes
        mg(lev)%n3 = 1

        CALL gather_h( d1(0:mg(lev-1)%n1+1,0:mg(lev-1)%n2+1,1),  d1(0:mg(lev)%n1+1,0:mg(lev)%n2+1,1), default=1.0D0)
        CALL gather_h( d2(0:mg(lev-1)%n1+1,0:mg(lev-1)%n2+1,1),  d2(0:mg(lev)%n1+1,0:mg(lev)%n2+1,1), default=1.0D0)
        CALL gather_h( d3(0:mg(lev-1)%n1+1,0:mg(lev-1)%n2+1,1),  d3(0:mg(lev)%n1+1,0:mg(lev)%n2+1,1), default=1.0D0)

        CALL gather_h(ds1(0:mg(lev-1)%n1+1,0:mg(lev-1)%n2+1,1), ds1(0:mg(lev)%n1+1,0:mg(lev)%n2+1,1), default=1.0D0)
        CALL gather_h(ds2(0:mg(lev-1)%n1+1,0:mg(lev-1)%n2+1,1), ds2(0:mg(lev)%n1+1,0:mg(lev)%n2+1,1), default=1.0D0)
        CALL gather_h(ds3(0:mg(lev-1)%n1+1,0:mg(lev-1)%n2+1,1), ds3(0:mg(lev)%n1+1,0:mg(lev)%n2+1,1), default=1.0D0)

        CALL gather_h( dv(0:mg(lev-1)%n1+1,0:mg(lev-1)%n2+1,1),  dv(0:mg(lev)%n1+1,0:mg(lev)%n2+1,1))
        CALL gather_h( da(0:mg(lev-1)%n1+1,0:mg(lev-1)%n2+1,1),  da(0:mg(lev)%n1+1,0:mg(lev)%n2+1,1))

        CALL gather_h( dm(0:mg(lev-1)%n1+1,0:mg(lev-1)%n2+1,1),  dm(0:mg(lev)%n1+1,0:mg(lev)%n2+1,1))

        IF (hrank == 0) THEN
           ds1(0,1:mg(lev)%n2,1) = 0.0D0
           ds2(1:mg(lev)%n1,0,1) = 0.0D0
           ds3(1:mg(lev)%n1,1:mg(lev)%n2,0) = 0.0D0

           dm(0,           1:mg(lev)%n2,1) = 0.0D0
           dm(mg(lev)%n1+1,1:mg(lev)%n2,1) = 0.0D0

           dm(1:mg(lev)%n1,0,           1) = 0.0D0
           dm(1:mg(lev)%n1,mg(lev)%n2+1,1) = 0.0D0

           dm(1:mg(lev)%n1,1:mg(lev)%n2,0) = 0.0D0
           dm(1:mg(lev)%n1,1:mg(lev)%n2,2) = 0.0D0
        END IF
    !cyclic-boundary is NOT cared?(2018/8/22)

     ELSE IF (lev == glev+1) THEN
        CALL assert(.FALSE., "GLEV and VLEV inconsistent")
#endif
     ELSE
#ifdef PARALLEL_MPI
        IF (lev > glev) THEN
           d1_ = sum(d1(1:mg(lev-1)%n1,1:mg(lev-1)%n2,1:mg(lev-1)%n3)) &
                / (mg(lev-1)%n1*mg(lev-1)%n2*mg(lev-1)%n3)
           d2_ = sum(d2(1:mg(lev-1)%n1,1:mg(lev-1)%n2,1:mg(lev-1)%n3)) &
                / (mg(lev-1)%n1*mg(lev-1)%n2*mg(lev-1)%n3)
           d3_ = sum(d3(1:mg(lev-1)%n1,1:mg(lev-1)%n2,1:mg(lev-1)%n3)) &
                / (mg(lev-1)%n1*mg(lev-1)%n2*mg(lev-1)%n3)
        ELSE IF (lev > vlev) THEN
           d1_ = mean3d(d1(1:mg(lev-1)%n1,1:mg(lev-1)%n2,1:mg(lev-1)%n3), hcomm)
           d2_ = mean3d(d2(1:mg(lev-1)%n1,1:mg(lev-1)%n2,1:mg(lev-1)%n3), hcomm)
           d3_ = mean3d(d3(1:mg(lev-1)%n1,1:mg(lev-1)%n2,1:mg(lev-1)%n3), hcomm)
        ELSE
           d1_ = mean3d(d1(1:mg(lev-1)%n1,1:mg(lev-1)%n2,1:mg(lev-1)%n3), comm)
           d2_ = mean3d(d2(1:mg(lev-1)%n1,1:mg(lev-1)%n2,1:mg(lev-1)%n3), comm)
           d3_ = mean3d(d3(1:mg(lev-1)%n1,1:mg(lev-1)%n2,1:mg(lev-1)%n3), comm)
        END IF
#else
        d1_ = sum(d1(1:mg(lev-1)%n1,1:mg(lev-1)%n2,1:mg(lev-1)%n3)) &
             / (mg(lev-1)%n1*mg(lev-1)%n2*mg(lev-1)%n3)
        d2_ = sum(d2(1:mg(lev-1)%n1,1:mg(lev-1)%n2,1:mg(lev-1)%n3)) &
             / (mg(lev-1)%n1*mg(lev-1)%n2*mg(lev-1)%n3)
        d3_ = sum(d3(1:mg(lev-1)%n1,1:mg(lev-1)%n2,1:mg(lev-1)%n3)) &
             / (mg(lev-1)%n1*mg(lev-1)%n2*mg(lev-1)%n3)
#endif
        mg(lev-1)%q1(:) = 0
        mg(lev-1)%q2(:) = 0
        mg(lev-1)%q3(:) = 0

        mg(lev-1)%p1(:) = 0
        mg(lev-1)%p2(:) = 0
        mg(lev-1)%p3(:) = 0


        IF (mg(lev-1)%n1 == 1 .OR. (lev < glev .AND. mg(lev-1)%n1 <= MG_GDIM) &
                              .OR. (d1_ > d2_*1.5 .AND. mg(lev-1)%n2 > 1)     &
                              .OR. (d1_ > d3_*1.5 .AND. mg(lev-1)%n3 > 1)) THEN
           mg(lev)%n1 = mg(lev-1)%n1
           DO i=1, mg(lev)%n1
              mg(lev-1)%p1(i) = i
              mg(lev-1)%q1(i) = i
           END DO
        ELSE
           mg(lev)%n1 = mg(lev-1)%n1/2
           DO i=1, mg(lev)%n1
              mg(lev-1)%p1(i) = 2*i
              mg(lev-1)%q1(2*i-1:2*i) = i
           END DO

           IF (mg(lev-1)%n1 /= mg(lev)%n1*2) THEN
              mg(lev)%n1 = mg(lev)%n1+1
              mg(lev-1)%p1(mg(lev  )%n1) = mg(lev-1)%n1
              mg(lev-1)%q1(mg(lev-1)%n1) = mg(lev  )%n1
           END IF
        END IF

        IF (mg(lev-1)%n2 == 1 .OR. (lev < glev .AND. mg(lev-1)%n2 <= MG_GDIM) &
                              .OR. (d2_ > d1_*1.5 .AND. mg(lev-1)%n1 > 1)     &
                              .OR. (d2_ > d3_*1.5 .AND. mg(lev-1)%n3 > 1)) THEN
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

        DO k=1, mg(lev-1)%n3
#ifdef PARALLEL_MPI
           IF (lev > glev) THEN
              d3_mean(k) = sum(d3(1:mg(lev-1)%n1,1:mg(lev-1)%n2,k))/(mg(lev-1)%n1*mg(lev-1)%n2)
           ELSE IF (lev > vlev) THEN
              d3_mean(k) = mean2d(d3(1:mg(lev-1)%n1,1:mg(lev-1)%n2,k), hcomm)
           ELSE
              d3_mean(k) = mean2d(d3(1:mg(lev-1)%n1,1:mg(lev-1)%n2,k), comm)
           END IF
#else
           d3_mean(k) = sum(d3(1:mg(lev-1)%n1,1:mg(lev-1)%n2,k))/(mg(lev-1)%n1*mg(lev-1)%n2)
#endif
        END DO
        d3_mean(mg(lev-1)%n3+1) = 0.0D0

        dmax = max(d3_*4.0D0, maxval(d3_mean(1:mg(lev-1)%n3)))
!        dmax = d3_*2.0D0
        IF (mg(lev)%n1 > 1) dmax = min(dmax,d1_*2.0D0)
        IF (mg(lev)%n2 > 1) dmax = min(dmax,d2_*2.0D0)

        mg(lev-1)%p3(0) = 0
        DO k=1, mg(lev-1)%n3
           d3_tmp = 0.0D0

           l = mg(lev-1)%p3(k-1)
           DO WHILE(l < mg(lev-1)%n3)
              l = l+1
              d3_tmp = d3_tmp + d3_mean(l)
              mg(lev-1)%q3(l) = k
              mg(lev-1)%p3(k) = l
              IF (d3_tmp+d3_mean(l+1) > dmax) EXIT
           END DO
           mg(lev)%n3 = k
           IF (l >= mg(lev-1)%n3) EXIT
        END DO

        IF (mg(lev)%n3 > 1 .AND. mg(lev)%n3==mg(lev-1)%n3) THEN
              mg(lev)%n3 = mg(lev)%n3 - 1
              mg(lev-1)%p3(mg(lev  )%n3) = mg(lev-1)%n3
              mg(lev-1)%q3(mg(lev-1)%n3) = mg(lev  )%n3
        END IF

        DO k=1, mg(lev)%n3
        DO j=1, mg(lev)%n2
        DO i=1, mg(lev)%n1
           d1(i,j,k) = sum(d1(mg(lev-1)%p1(i-1)+1:mg(lev-1)%p1(i),  &
                              mg(lev-1)%p2(j-1)+1:mg(lev-1)%p2(j),  &
                              mg(lev-1)%p3(k-1)+1:mg(lev-1)%p3(k))) &
                      / ((mg(lev-1)%p2(j)-mg(lev-1)%p2(j-1))*(mg(lev-1)%p3(k)-mg(lev-1)%p3(k-1)))
        END DO
        END DO
        END DO

        DO k=1, mg(lev)%n3
        DO j=1, mg(lev)%n2
        DO i=1, mg(lev)%n1
           d2(i,j,k) = sum(d2(mg(lev-1)%p1(i-1)+1:mg(lev-1)%p1(i),  &
                              mg(lev-1)%p2(j-1)+1:mg(lev-1)%p2(j),  &
                              mg(lev-1)%p3(k-1)+1:mg(lev-1)%p3(k))) &
                      / ((mg(lev-1)%p1(i)-mg(lev-1)%p1(i-1))*(mg(lev-1)%p3(k)-mg(lev-1)%p3(k-1)))
        END DO
        END DO
        END DO

        DO k=1, mg(lev)%n3
        DO j=1, mg(lev)%n2
        DO i=1, mg(lev)%n1
           d3(i,j,k) = sum(d3(mg(lev-1)%p1(i-1)+1:mg(lev-1)%p1(i),  &
                              mg(lev-1)%p2(j-1)+1:mg(lev-1)%p2(j),  &
                              mg(lev-1)%p3(k-1)+1:mg(lev-1)%p3(k))) &
                      / ((mg(lev-1)%p1(i)-mg(lev-1)%p1(i-1))*(mg(lev-1)%p2(j)-mg(lev-1)%p2(j-1)))
        END DO
        END DO
        END DO

        DO k=1, mg(lev)%n3
        DO j=1, mg(lev)%n2
        DO i=0, mg(lev)%n1
           ds1(i,j,k) = sum(ds1(mg(lev-1)%p1(i),                    &
                                mg(lev-1)%p2(j-1)+1:mg(lev-1)%p2(j), &
                                mg(lev-1)%p3(k-1)+1:mg(lev-1)%p3(k)))
        END DO
        END DO
        END DO

        DO k=1, mg(lev)%n3
        DO j=0, mg(lev)%n2
        DO i=1, mg(lev)%n1
           ds2(i,j,k) = sum(ds2(mg(lev-1)%p1(i-1)+1:mg(lev-1)%p1(i), &
                                mg(lev-1)%p2(j),                    &
                                mg(lev-1)%p3(k-1)+1:mg(lev-1)%p3(k)))
        END DO
        END DO
        END DO

        DO k=0, mg(lev)%n3
        DO j=1, mg(lev)%n2
        DO i=1, mg(lev)%n1
           ds3(i,j,k) = sum(ds3(mg(lev-1)%p1(i-1)+1:mg(lev-1)%p1(i), &
                                mg(lev-1)%p2(j-1)+1:mg(lev-1)%p2(j), &
                                mg(lev-1)%p3(k)))
        END DO
        END DO
        END DO

        DO k=1, mg(lev)%n3
        DO j=1, mg(lev)%n2
        DO i=1, mg(lev)%n1
           dv(i,j,k) = sum(dv(mg(lev-1)%p1(i-1)+1:mg(lev-1)%p1(i), &
                              mg(lev-1)%p2(j-1)+1:mg(lev-1)%p2(j), &
                              mg(lev-1)%p3(k-1)+1:mg(lev-1)%p3(k)))
        END DO
        END DO
        END DO

        DO k=1, mg(lev)%n3
        DO j=1, mg(lev)%n2
        DO i=1, mg(lev)%n1
           da(i,j,k) = sum(da(mg(lev-1)%p1(i-1)+1:mg(lev-1)%p1(i), &
                              mg(lev-1)%p2(j-1)+1:mg(lev-1)%p2(j), &
                              mg(lev-1)%p3(k-1)+1:mg(lev-1)%p3(k)))
        END DO
        END DO
        END DO

        dm(0,0:mg(lev)%n2+1,0:mg(lev)%n3+1) = 0.0D0
        dm(0:mg(lev)%n1+1,0,0:mg(lev)%n3+1) = 0.0D0
        dm(0:mg(lev)%n1+1,0:mg(lev)%n2+1,0) = 0.0D0
        DO k=1, mg(lev)%n3
        DO j=1, mg(lev)%n2
        DO i=1, mg(lev)%n1
           dm(i,j,k) = 1.0D0 - product(1.0D0-dm(mg(lev-1)%p1(i-1)+1:mg(lev-1)%p1(i), &
                                                mg(lev-1)%p2(j-1)+1:mg(lev-1)%p2(j), &
                                                mg(lev-1)%p3(k-1)+1:mg(lev-1)%p3(k)))
        END DO
        END DO
        END DO
        dm(  mg(lev)%n1+1,0:mg(lev)%n2+1,0:mg(lev)%n3+1) = 0.0D0
        dm(0:mg(lev)%n1+1,  mg(lev)%n2+1,0:mg(lev)%n3+1) = 0.0D0
        dm(0:mg(lev)%n1+1,0:mg(lev)%n2+1,  mg(lev)%n3+1) = 0.0D0
     END IF

     IF (rank == 0) THEN
#ifdef PARALLEL_MPI
        IF (lev > glev) THEN
           WRITE(REPORT_UNIT, '("  lev",I3,": dim = (",I4,"x",I3,",",I4,"x",I3,",",I4,"x",I3,")")') &
                lev, mg(lev)%n1, 1, mg(lev)%n2, 1, mg(lev)%n3, 1
        ELSE IF (lev > vlev) THEN
           WRITE(REPORT_UNIT, '("  lev",I3,": dim = (",I4,"x",I3,",",I4,"x",I3,",",I4,"x",I3,")", ES10.3)') &
                lev, mg(lev)%n1, ipes, mg(lev)%n2, jpes, mg(lev)%n3, 1
        ELSE
           WRITE(REPORT_UNIT, '("  lev",I3,": dim = (",I4,"x",I3,",",I4,"x",I3,",",I4,"x",I3,")", ES10.3)') &
                lev, mg(lev)%n1, ipes, mg(lev)%n2, jpes, mg(lev)%n3, kpes
        END IF
#else
        WRITE(REPORT_UNIT, '("  lev",I3,": dim = (",I4,",",I4,",",I4,")")') lev, mg(lev)%n1, mg(lev)%n2, mg(lev)%n3
#endif
     END IF

#ifdef PARALLEL_MPI
     IF (lev <= glev) THEN
        mg(lev)%rank_e = ranks(icoord+1, jcoord,   kcoord)
        mg(lev)%rank_w = ranks(icoord-1, jcoord,   kcoord)
        mg(lev)%rank_n = ranks(icoord,   jcoord+1, kcoord)
        mg(lev)%rank_s = ranks(icoord,   jcoord-1, kcoord)
        mg(lev)%rank_u = ranks(icoord,   jcoord,   kcoord+1)
        mg(lev)%rank_l = ranks(icoord,   jcoord,   kcoord-1)

#ifdef PARALLEL3D
        IF (lev > vlev) THEN
           IF (vrank==0) THEN
              mg(lev)%rank_e = ranks(icoord+1, jcoord,   kcoord)
              mg(lev)%rank_w = ranks(icoord-1, jcoord,   kcoord)
              mg(lev)%rank_n = ranks(icoord,   jcoord+1, kcoord)
              mg(lev)%rank_s = ranks(icoord,   jcoord-1, kcoord)
           ELSE
              mg(lev)%rank_e = MPI_PROC_NULL
              mg(lev)%rank_w = MPI_PROC_NULL
              mg(lev)%rank_n = MPI_PROC_NULL
              mg(lev)%rank_s = MPI_PROC_NULL
           END IF
           mg(lev)%rank_u = MPI_PROC_NULL
           mg(lev)%rank_l = MPI_PROC_NULL
        END IF
#endif
        CALL mpi_recv_init(recvbuf_w_r4, mg(lev)%n2, MPI_REAL4, mg(lev)%rank_w, 1, comm, req2d_r4(0,lev), ierr)
        CALL mpi_recv_init(recvbuf_w_r8, mg(lev)%n2, MPI_REAL8, mg(lev)%rank_w, 1, comm, req2d_r8(0,lev), ierr)
        CALL mpi_recv_init(recvbuf_e_r4, mg(lev)%n2, MPI_REAL4, mg(lev)%rank_e, 2, comm, req2d_r4(1,lev), ierr)
        CALL mpi_recv_init(recvbuf_e_r8, mg(lev)%n2, MPI_REAL8, mg(lev)%rank_e, 2, comm, req2d_r8(1,lev), ierr)
        CALL mpi_recv_init(recvbuf_s_r4, mg(lev)%n1, MPI_REAL4, mg(lev)%rank_s, 3, comm, req2d_r4(2,lev), ierr)
        CALL mpi_recv_init(recvbuf_s_r8, mg(lev)%n1, MPI_REAL8, mg(lev)%rank_s, 3, comm, req2d_r8(2,lev), ierr)
        CALL mpi_recv_init(recvbuf_n_r4, mg(lev)%n1, MPI_REAL4, mg(lev)%rank_n, 4, comm, req2d_r4(3,lev), ierr)
        CALL mpi_recv_init(recvbuf_n_r8, mg(lev)%n1, MPI_REAL8, mg(lev)%rank_n, 4, comm, req2d_r8(3,lev), ierr)

        CALL mpi_send_init(sendbuf_e_r4, mg(lev)%n2, MPI_REAL4, mg(lev)%rank_e, 1, comm, req2d_r4(4,lev), ierr)
        CALL mpi_send_init(sendbuf_e_r8, mg(lev)%n2, MPI_REAL8, mg(lev)%rank_e, 1, comm, req2d_r8(4,lev), ierr)
        CALL mpi_send_init(sendbuf_w_r4, mg(lev)%n2, MPI_REAL4, mg(lev)%rank_w, 2, comm, req2d_r4(5,lev), ierr)
        CALL mpi_send_init(sendbuf_w_r8, mg(lev)%n2, MPI_REAL8, mg(lev)%rank_w, 2, comm, req2d_r8(5,lev), ierr)
        CALL mpi_send_init(sendbuf_n_r4, mg(lev)%n1, MPI_REAL4, mg(lev)%rank_n, 3, comm, req2d_r4(6,lev), ierr)
        CALL mpi_send_init(sendbuf_n_r8, mg(lev)%n1, MPI_REAL8, mg(lev)%rank_n, 3, comm, req2d_r8(6,lev), ierr)
        CALL mpi_send_init(sendbuf_s_r4, mg(lev)%n1, MPI_REAL4, mg(lev)%rank_s, 4, comm, req2d_r4(7,lev), ierr)
        CALL mpi_send_init(sendbuf_s_r8, mg(lev)%n1, MPI_REAL8, mg(lev)%rank_s, 4, comm, req2d_r8(7,lev), ierr)

        CALL mpi_recv_init(recvbuf_w_r4, mg(lev)%n2*mg(lev)%n3, MPI_REAL4, mg(lev)%rank_w, 1, comm, req3d_r4(0,lev), ierr)
        CALL mpi_recv_init(recvbuf_w_r8, mg(lev)%n2*mg(lev)%n3, MPI_REAL8, mg(lev)%rank_w, 1, comm, req3d_r8(0,lev), ierr)
        CALL mpi_recv_init(recvbuf_e_r4, mg(lev)%n2*mg(lev)%n3, MPI_REAL4, mg(lev)%rank_e, 2, comm, req3d_r4(1,lev), ierr)
        CALL mpi_recv_init(recvbuf_e_r8, mg(lev)%n2*mg(lev)%n3, MPI_REAL8, mg(lev)%rank_e, 2, comm, req3d_r8(1,lev), ierr)
        CALL mpi_recv_init(recvbuf_s_r4, mg(lev)%n1*mg(lev)%n3, MPI_REAL4, mg(lev)%rank_s, 3, comm, req3d_r4(2,lev), ierr)
        CALL mpi_recv_init(recvbuf_s_r8, mg(lev)%n1*mg(lev)%n3, MPI_REAL8, mg(lev)%rank_s, 3, comm, req3d_r8(2,lev), ierr)
        CALL mpi_recv_init(recvbuf_n_r4, mg(lev)%n1*mg(lev)%n3, MPI_REAL4, mg(lev)%rank_n, 4, comm, req3d_r4(3,lev), ierr)
        CALL mpi_recv_init(recvbuf_n_r8, mg(lev)%n1*mg(lev)%n3, MPI_REAL8, mg(lev)%rank_n, 4, comm, req3d_r8(3,lev), ierr)

        CALL mpi_send_init(sendbuf_e_r4, mg(lev)%n2*mg(lev)%n3, MPI_REAL4, mg(lev)%rank_e, 1, comm, req3d_r4(6,lev), ierr)
        CALL mpi_send_init(sendbuf_e_r8, mg(lev)%n2*mg(lev)%n3, MPI_REAL8, mg(lev)%rank_e, 1, comm, req3d_r8(6,lev), ierr)
        CALL mpi_send_init(sendbuf_w_r4, mg(lev)%n2*mg(lev)%n3, MPI_REAL4, mg(lev)%rank_w, 2, comm, req3d_r4(7,lev), ierr)
        CALL mpi_send_init(sendbuf_w_r8, mg(lev)%n2*mg(lev)%n3, MPI_REAL8, mg(lev)%rank_w, 2, comm, req3d_r8(7,lev), ierr)
        CALL mpi_send_init(sendbuf_n_r4, mg(lev)%n1*mg(lev)%n3, MPI_REAL4, mg(lev)%rank_n, 3, comm, req3d_r4(8,lev), ierr)
        CALL mpi_send_init(sendbuf_n_r8, mg(lev)%n1*mg(lev)%n3, MPI_REAL8, mg(lev)%rank_n, 3, comm, req3d_r8(8,lev), ierr)
        CALL mpi_send_init(sendbuf_s_r4, mg(lev)%n1*mg(lev)%n3, MPI_REAL4, mg(lev)%rank_s, 4, comm, req3d_r4(9,lev), ierr)
        CALL mpi_send_init(sendbuf_s_r8, mg(lev)%n1*mg(lev)%n3, MPI_REAL8, mg(lev)%rank_s, 4, comm, req3d_r8(9,lev), ierr)

        req3d_r4( 4,lev) = MPI_REQUEST_NULL
        req3d_r4( 5,lev) = MPI_REQUEST_NULL
        req3d_r4(10,lev) = MPI_REQUEST_NULL
        req3d_r4(11,lev) = MPI_REQUEST_NULL

        req3d_r8( 4,lev) = MPI_REQUEST_NULL
        req3d_r8( 5,lev) = MPI_REQUEST_NULL
        req3d_r8(10,lev) = MPI_REQUEST_NULL
        req3d_r8(11,lev) = MPI_REQUEST_NULL

#ifdef PARALLEL3D
        IF (lev < vlev) THEN
           CALL mpi_recv_init(recvbuf_l_r4, mg(lev)%n1*mg(lev)%n2, MPI_REAL4, mg(lev)%rank_l, 5, comm, req3d_r4( 4,lev), ierr)
           CALL mpi_recv_init(recvbuf_l_r8, mg(lev)%n1*mg(lev)%n2, MPI_REAL8, mg(lev)%rank_l, 5, comm, req3d_r8( 4,lev), ierr)
           CALL mpi_recv_init(recvbuf_u_r4, mg(lev)%n1*mg(lev)%n2, MPI_REAL4, mg(lev)%rank_u, 6, comm, req3d_r4( 5,lev), ierr)
           CALL mpi_recv_init(recvbuf_u_r8, mg(lev)%n1*mg(lev)%n2, MPI_REAL8, mg(lev)%rank_u, 6, comm, req3d_r8( 5,lev), ierr)

           CALL mpi_send_init(sendbuf_u_r4, mg(lev)%n1*mg(lev)%n2, MPI_REAL4, mg(lev)%rank_u, 5, comm, req3d_r4(10,lev), ierr)
           CALL mpi_send_init(sendbuf_u_r8, mg(lev)%n1*mg(lev)%n2, MPI_REAL8, mg(lev)%rank_u, 5, comm, req3d_r8(10,lev), ierr)
           CALL mpi_send_init(sendbuf_l_r4, mg(lev)%n1*mg(lev)%n2, MPI_REAL4, mg(lev)%rank_l, 6, comm, req3d_r4(11,lev), ierr)
           CALL mpi_send_init(sendbuf_l_r8, mg(lev)%n1*mg(lev)%n2, MPI_REAL8, mg(lev)%rank_l, 6, comm, req3d_r8(11,lev), ierr)
        ELSE
           IF (vrank/=0) THEN
              req3d_r4(:,lev) = MPI_REQUEST_NULL
              req3d_r8(:,lev) = MPI_REQUEST_NULL
           END IF
        END IF
#endif
     END IF
#endif

     IF (open_e) THEN
        DO k=0, mg(lev)%n3+1
        DO j=0, mg(lev)%n2+1
           d1(mg(lev)%n1+1,j,k) = d1(mg(lev)%n1,j,k)
           d2(mg(lev)%n1+1,j,k) = d2(mg(lev)%n1,j,k)
           d3(mg(lev)%n1+1,j,k) = d3(mg(lev)%n1,j,k)
           !dm(mg(lev)%n1+1,j,k) = dm(mg(lev)%n1,j,k)
           dm(mg(lev)%n1+1,j,k) = 0.0D0
        END DO
        END DO
     END IF

     IF (open_w) THEN
        DO k=0, mg(lev)%n3+1
        DO j=0, mg(lev)%n2+1
           d1(0,j,k) = d1(1,j,k)
           d2(0,j,k) = d2(1,j,k)
           d3(0,j,k) = d3(1,j,k)
           !dm(0,j,k) = dm(1,j,k)
           dm(0,j,k) = 0.0D0
        END DO
        END DO
     END IF

     IF (open_n) THEN
        DO k=0, mg(lev)%n3+1
        DO i=0, mg(lev)%n1+1
           d1(i,mg(lev)%n2+1,k) = d1(i,mg(lev)%n2,k)
           d2(i,mg(lev)%n2+1,k) = d2(i,mg(lev)%n2,k)
           d3(i,mg(lev)%n2+1,k) = d3(i,mg(lev)%n2,k)
           !dm(i,mg(lev)%n2+1,k) = dm(i,mg(lev)%n2,k)
           dm(i,mg(lev)%n2+1,k) = 0.0D0
        END DO
        END DO
     END IF

     IF (open_s) THEN
        DO k=0, mg(lev)%n3+1
        DO i=0, mg(lev)%n1+1
           d1(i,0,k) = d1(i,1,k)
           d2(i,0,k) = d2(i,1,k)
           d3(i,0,k) = d3(i,1,k)
           !dm(i,0,k) = dm(i,1,k)
           dm(i,0,k) = 0.0D0
        END DO
        END DO
     END IF

     IF (mg(lev)%n3 /= 1) THEN
        DO j=0, mg(lev)%n2+1
        DO i=0, mg(lev)%n1+1
           d1(i,j,0) = d1(i,j,1)
           d2(i,j,0) = d2(i,j,1)
           d3(i,j,0) = d3(i,j,1)

           d1(i,j,mg(lev)%n3+1) = d1(i,j,mg(lev)%n3)
           d2(i,j,mg(lev)%n3+1) = d2(i,j,mg(lev)%n3)
           d3(i,j,mg(lev)%n3+1) = d3(i,j,mg(lev)%n3)

  !        dm(i,j,0)            = dm(i,j,1)
  !        dm(i,j,mg(lev)%n3+1) = dm(i,j,mg(lev)%n3)
        END DO
        END DO
     END IF

#ifdef PARALLEL_MPI
     IF ((lev > vlev .OR. kpes==1).AND. mg(lev)%n3==1) THEN
#else
     IF (mg(lev)%n3==1) THEN
#endif
        ALLOCATE(mg(lev)%l(-2:2,1:mg(lev)%n1,1:mg(lev)%n2,1))
        ALLOCATE(mg(lev)%m(-2:2,1:mg(lev)%n1,1:mg(lev)%n2,1))

#ifdef PARALLEL_MPI
        IF (lev > glev .AND. rank/=0) THEN
#else
        IF (rank/=0) THEN
#endif
           mg(lev)%l(:,:,:,:) = 0.0D0
           mg(lev)%m(:,:,:,:) = 0.0D0
        ELSE
           ALLOCATE(dtmp(0:mg(lev)%n1+1,0:mg(lev)%n2+1,1:1))

           dtmp(:,:,1) = d1(0:mg(lev)%n1+1,0:mg(lev)%n2+1,1)
           CALL sync2d(lev, dtmp(:,:,1))
           d1(0:mg(lev)%n1+1,0:mg(lev)%n2+1,1) = dtmp(:,:,1)

           dtmp(:,:,1) = d2(0:mg(lev)%n1+1,0:mg(lev)%n2+1,1)
           CALL sync2d(lev, dtmp(:,:,1))
           d2(0:mg(lev)%n1+1,0:mg(lev)%n2+1,1) = dtmp(:,:,1)

           dtmp(:,:,1) = d3(0:mg(lev)%n1+1,0:mg(lev)%n2+1,1)
           CALL sync2d(lev, dtmp(:,:,1))
           d3(0:mg(lev)%n1+1,0:mg(lev)%n2+1,1) = dtmp(:,:,1)

           dtmp(:,:,1) = ds1(0:mg(lev)%n1+1,0:mg(lev)%n2+1,1)
           CALL sync2d(lev, dtmp(:,:,1))
           ds1(0:mg(lev)%n1+1,0:mg(lev)%n2+1,1) = dtmp(:,:,1)

           dtmp(:,:,1) = ds2(0:mg(lev)%n1+1,0:mg(lev)%n2+1,1)
           CALL sync2d(lev, dtmp(:,:,1))
           ds2(0:mg(lev)%n1+1,0:mg(lev)%n2+1,1) = dtmp(:,:,1)

           dtmp(:,:,1) = ds3(0:mg(lev)%n1+1,0:mg(lev)%n2+1,1)
           CALL sync2d(lev, dtmp(:,:,1))
           ds3(0:mg(lev)%n1+1,0:mg(lev)%n2+1,1) = dtmp(:,:,1)

           dtmp(:,:,1) = dv(0:mg(lev)%n1+1,0:mg(lev)%n2+1,1)
           CALL sync2d(lev, dtmp(:,:,1))
           dv(0:mg(lev)%n1+1,0:mg(lev)%n2+1,1) = dtmp(:,:,1)

           dtmp(:,:,1) = da(0:mg(lev)%n1+1,0:mg(lev)%n2+1,1)
           CALL sync2d(lev, dtmp(:,:,1))
           da(0:mg(lev)%n1+1,0:mg(lev)%n2+1,1) = dtmp(:,:,1)

           dtmp(:,:,1) = dm(0:mg(lev)%n1+1,0:mg(lev)%n2+1,1)
           CALL sync2d(lev, dtmp(:,:,1))
           dm(0:mg(lev)%n1+1,0:mg(lev)%n2+1,1) = dtmp(:,:,1)

           DEALLOCATE(dtmp)

           ALLOCATE(ltmp(-2:2,1:mg(lev)%n1,1:mg(lev)%n2,1))
           ALLOCATE(mtmp(-2:2,1:mg(lev)%n1,1:mg(lev)%n2,1))

#ifdef PARALLEL_MPI
           IF (mg(lev)%n1 == 1 .AND. mg(lev)%rank_e == rank) THEN
#else
           IF (mg(lev)%n1 == 1) THEN
#endif
              DO j=1, mg(lev)%n2
                 ltmp(-1,1,j,1) = 0.0D0
                 ltmp( 1,1,j,1) = 0.0D0
              END DO
           ELSE
              SELECT CASE (coeffscheme)
              CASE (1)
                 DO j=1, mg(lev)%n2
                 DO i=1, mg(lev)%n1
                    ltmp(-1,i,j,1) = ds1(i-1,j,1) * 2.0/(d1(i,j,1)+d1(i-1,j,1)) * dm(i,j,1)*dm(i-1,j,1)
                    ltmp( 1,i,j,1) = ds1(i,  j,1) * 2.0/(d1(i,j,1)+d1(i+1,j,1)) * dm(i,j,1)*dm(i+1,j,1)
                 END DO
                 END DO
              CASE (2)
                 DO j=1, mg(lev)%n2
                 DO i=1, mg(lev)%n1
                    ltmp(-1,i,j,1) = ds1(i-1,j,1) * d2(i-1,j,1)*d3(i-1,j,1) * 2.0/(dv(i,j,k)+dv(i-1,j,k)) * dm(i,j,1)*dm(i-1,j,1)
                    ltmp( 1,i,j,1) = ds1(i,  j,1) * d2(i+1,j,1)*d3(i+1,j,1) * 2.0/(dv(i,j,k)+dv(i+1,j,k)) * dm(i,j,1)*dm(i+1,j,1)
                 END DO
                 END DO
              CASE (3)
                 DO j=1, mg(lev)%n2
                 DO i=1, mg(lev)%n1
                    ltmp(-1,i,j,1) = ds1(i-1,j,1) / d1(i-1,j,1) * 2.0/(dv(i,j,k)+dv(i-1,j,k)) * dm(i,j,1)*dm(i-1,j,1)
                    ltmp( 1,i,j,1) = ds1(i,  j,1) / d1(i+1,j,1) * 2.0/(dv(i,j,k)+dv(i+1,j,k)) * dm(i,j,1)*dm(i+1,j,1)
                 END DO
                 END DO
              END SELECT
           END IF

#ifdef PARALLEL_MPI
           IF (mg(lev)%n2 == 1 .AND. mg(lev)%rank_n == rank) THEN
#else
           IF (mg(lev)%n2 == 1) THEN
#endif
              DO i=1, mg(lev)%n1
                 ltmp(-2,i,1,1) = 0.0D0
                 ltmp( 2,i,1,1) = 0.0D0
              END DO
           ELSE
              SELECT CASE (coeffscheme)
              CASE (1)
                 DO j=1, mg(lev)%n2
                 DO i=1, mg(lev)%n1
                    ltmp(-2,i,j,1) = ds2(i,j-1,1) * 2.0/(d2(i,j,1)+d2(i,j-1,1)) * dm(i,j,1)*dm(i,j-1,1)
                    ltmp( 2,i,j,1) = ds2(i,j,  1) * 2.0/(d2(i,j,1)+d2(i,j+1,1)) * dm(i,j,1)*dm(i,j+1,1)
                 END DO
                 END DO
              CASE (2)
                 DO j=1, mg(lev)%n2
                 DO i=1, mg(lev)%n1
                    ltmp(-2,i,j,1) = ds2(i,j-1,1) * d1(i,j-1,1)*d3(i,j-1,1) * 2.0/(dv(i,j,1)+dv(i,j-1,1)) * dm(i,j,1)*dm(i,j-1,1)
                    ltmp( 2,i,j,1) = ds2(i,j,  1) * d1(i,j+1,1)*d3(i,j+1,1) * 2.0/(dv(i,j,1)+dv(i,j+1,1)) * dm(i,j,1)*dm(i,j+1,1)
                 END DO
                 END DO
              CASE (3)
                 DO j=1, mg(lev)%n2
                 DO i=1, mg(lev)%n1
                    ltmp(-2,i,j,1) = ds2(i,j-1,1) / d2(i,j-1,1) * 2.0/(dv(i,j,1)+dv(i,j-1,1)) * dm(i,j,1)*dm(i,j-1,1)
                    ltmp( 2,i,j,1) = ds2(i,j,  1) / d2(i,j+1,1) * 2.0/(dv(i,j,1)+dv(i,j+1,1)) * dm(i,j,1)*dm(i,j+1,1)
                 END DO
                 END DO
              END SELECT
           END IF

           SELECT CASE (coeffscheme)
           CASE (1)
              DO j=1, mg(lev)%n2
              DO i=1, mg(lev)%n1
                 ltmp( 0,i,j,1) = -da(i,j,1) - ltmp(-1,i,j,1) - ltmp(1,i,j,1) &
                                             - ltmp(-2,i,j,1) - ltmp(2,i,j,1)
              END DO
              END DO
           CASE (2)
              DO j=1, mg(lev)%n2
              DO i=1, mg(lev)%n1
                 ltmp( 0,i,j,1) = -da(i,j,1) - ds1(i-1,j,1) * d2(i,j,1)*d3(i,j,1) * 2.0/(dv(i,j,k)+dv(i-1,j,k)) * dm(i,j,1)*dm(i-1,j,1) &
                                             - ds1(i,  j,1) * d2(i,j,1)*d3(i,j,1) * 2.0/(dv(i,j,k)+dv(i+1,j,k)) * dm(i,j,1)*dm(i+1,j,1) &
                                             - ds2(i,j-1,1) * d1(i,j,1)*d3(i,j,1) * 2.0/(dv(i,j,k)+dv(i,j-1,k)) * dm(i,j,1)*dm(i,j-1,1) &
                                             - ds2(i,j,  1) * d1(i,j,1)*d3(i,j,1) * 2.0/(dv(i,j,k)+dv(i,j+1,k)) * dm(i,j,1)*dm(i,j+1,1)
              END DO
              END DO
           CASE (3)
              DO j=1, mg(lev)%n2
              DO i=1, mg(lev)%n1
                 ltmp( 0,i,j,1) = -da(i,j,1) - ds1(i-1,j,1) / d1(i,j,1) * 2.0/(dv(i,j,k)+dv(i-1,j,k)) * dm(i,j,1)*dm(i-1,j,1) &
                                             - ds1(i,  j,1) / d1(i,j,1) * 2.0/(dv(i,j,k)+dv(i+1,j,k)) * dm(i,j,1)*dm(i+1,j,1) &
                                             - ds2(i,j-1,1) / d2(i,j,1) * 2.0/(dv(i,j,k)+dv(i,j-1,k)) * dm(i,j,1)*dm(i,j-1,1) &
                                             - ds2(i,j,  1) / d2(i,j,1) * 2.0/(dv(i,j,k)+dv(i,j+1,k)) * dm(i,j,1)*dm(i,j+1,1)
              END DO
              END DO
           END SELECT

           CALL sai2d(lev, ltmp, mtmp)

           mg(lev)%l(:,:,:,:) = real(ltmp(:,:,:,:), KIND=MG_KIND)
           mg(lev)%m(:,:,:,:) = real(mtmp(:,:,:,:), KIND=MG_KIND)
        END IF

        ALLOCATE(mg(lev)%r(0:mg(lev)%n1+1, 0:mg(lev)%n2+1, 1))
        ALLOCATE(mg(lev)%e(0:mg(lev)%n1+1, 0:mg(lev)%n2+1, 1))
        ALLOCATE(mg(lev)%t(0:mg(lev)%n1+1, 0:mg(lev)%n2+1, 1))
     ELSE
        ALLOCATE(dtmp(0:mg(lev)%n1+1, 0:mg(lev)%n2+1, 0:mg(lev)%n3+1))

        dtmp(:,:,:) = d1(0:mg(lev)%n1+1, 0:mg(lev)%n2+1, 0:mg(lev)%n3+1)
        CALL sync3d(lev, dtmp)
        d1(0:mg(lev)%n1+1, 0:mg(lev)%n2+1, 0:mg(lev)%n3+1) = dtmp(:,:,:)

        dtmp(:,:,:) = d2(0:mg(lev)%n1+1,0:mg(lev)%n2+1,0:mg(lev)%n3+1)
        CALL sync3d(lev, dtmp)
        d2(0:mg(lev)%n1+1,0:mg(lev)%n2+1,0:mg(lev)%n3+1) = dtmp(:,:,:)

        dtmp(:,:,:) = d3(0:mg(lev)%n1+1,0:mg(lev)%n2+1,0:mg(lev)%n3+1)
        CALL sync3d(lev, dtmp)
        d3(0:mg(lev)%n1+1,0:mg(lev)%n2+1,0:mg(lev)%n3+1) = dtmp(:,:,:)

        dtmp(:,:,:) = ds1(0:mg(lev)%n1+1,0:mg(lev)%n2+1,0:mg(lev)%n3+1)
        CALL sync3d(lev, dtmp)
        ds1(0:mg(lev)%n1+1,0:mg(lev)%n2+1,0:mg(lev)%n3+1) = dtmp(:,:,:)

        dtmp(:,:,:) = ds2(0:mg(lev)%n1+1,0:mg(lev)%n2+1,0:mg(lev)%n3+1)
        CALL sync3d(lev, dtmp)
        ds2(0:mg(lev)%n1+1,0:mg(lev)%n2+1,0:mg(lev)%n3+1) = dtmp(:,:,:)

        dtmp(:,:,:) = ds3(0:mg(lev)%n1+1,0:mg(lev)%n2+1,0:mg(lev)%n3+1)
        CALL sync3d(lev, dtmp)
        ds3(0:mg(lev)%n1+1,0:mg(lev)%n2+1,0:mg(lev)%n3+1) = dtmp(:,:,:)

        dtmp(:,:,:) = dv(0:mg(lev)%n1+1,0:mg(lev)%n2+1,0:mg(lev)%n3+1)
        CALL sync3d(lev, dtmp)
        dv(0:mg(lev)%n1+1,0:mg(lev)%n2+1,0:mg(lev)%n3+1) = dtmp(:,:,:)

        dtmp(:,:,:) = da(0:mg(lev)%n1+1,0:mg(lev)%n2+1,0:mg(lev)%n3+1)
        CALL sync3d(lev, dtmp)
        da(0:mg(lev)%n1+1,0:mg(lev)%n2+1,0:mg(lev)%n3+1) = dtmp(:,:,:)

        dtmp(:,:,:) = dm(0:mg(lev)%n1+1,0:mg(lev)%n2+1,0:mg(lev)%n3+1)
        CALL sync3d(lev, dtmp)
        dm(0:mg(lev)%n1+1,0:mg(lev)%n2+1,0:mg(lev)%n3+1) = dtmp(:,:,:)

        DEALLOCATE(dtmp)

        ALLOCATE(ltmp(-3:3,1:mg(lev)%n1,1:mg(lev)%n2,1:mg(lev)%n3))
        ALLOCATE(mtmp(-3:3,1:mg(lev)%n1,1:mg(lev)%n2,1:mg(lev)%n3))

#ifdef PARALLEL_MPI
        IF (mg(lev)%n1 == 1 .AND. mg(lev)%rank_e == rank .AND. mg(lev)%rank_w == rank) THEN
#else
        IF (mg(lev)%n1 == 1) THEN
#endif
           DO k=1, mg(lev)%n3
           DO j=1, mg(lev)%n2
              ltmp(-1,1,j,k) = 0.0D0
              ltmp( 1,1,j,k) = 0.0D0
           END DO
           END DO
        ELSE
           SELECT CASE (coeffscheme)
           CASE (1)
              DO k=1, mg(lev)%n3
              DO j=1, mg(lev)%n2
              DO i=1, mg(lev)%n1
                 ltmp(-1,i,j,k) = ds1(i-1,j,k) * 2.0/(d1(i,j,k)+d1(i-1,j,k)) * dm(i,j,k)*dm(i-1,j,k)
                 ltmp( 1,i,j,k) = ds1(i,  j,k) * 2.0/(d1(i,j,k)+d1(i+1,j,k)) * dm(i,j,k)*dm(i+1,j,k)
              END DO
              END DO
              END DO
           CASE (2)
              DO k=1, mg(lev)%n3
              DO j=1, mg(lev)%n2
              DO i=1, mg(lev)%n1
                 ltmp(-1,i,j,k) = ds1(i-1,j,k) * d2(i-1,j,k)*d3(i-1,j,k) * 2.0/(dv(i,j,k)+dv(i-1,j,k)) * dm(i,j,k)*dm(i-1,j,k)
                 ltmp( 1,i,j,k) = ds1(i,  j,k) * d2(i+1,j,k)*d3(i+1,j,k) * 2.0/(dv(i,j,k)+dv(i+1,j,k)) * dm(i,j,k)*dm(i+1,j,k)
              END DO
              END DO
              END DO
           CASE (3)
              DO k=1, mg(lev)%n3
              DO j=1, mg(lev)%n2
              DO i=1, mg(lev)%n1
                 ltmp(-1,i,j,k) = ds1(i-1,j,k) / d1(i-1,j,k) * 2.0/(dv(i,j,k)+dv(i-1,j,k)) * dm(i,j,k)*dm(i-1,j,k)
                 ltmp( 1,i,j,k) = ds1(i,  j,k) / d1(i+1,j,k) * 2.0/(dv(i,j,k)+dv(i+1,j,k)) * dm(i,j,k)*dm(i+1,j,k)
              END DO
              END DO
              END DO
           END SELECT
        END IF

#ifdef PARALLEL_MPI
        IF (mg(lev)%n2 == 1 .AND. mg(lev)%rank_n == rank .AND. mg(lev)%rank_s == rank) THEN
#else
        IF (mg(lev)%n2 == 1) THEN
#endif
           DO k=1, mg(lev)%n3
           DO i=1, mg(lev)%n1
              ltmp(-2,i,1,k) = 0.0D0
              ltmp( 2,i,1,k) = 0.0D0
           END DO
           END DO
        ELSE
           SELECT CASE (coeffscheme)
           CASE (1)
              DO k=1, mg(lev)%n3
              DO j=1, mg(lev)%n2
              DO i=1, mg(lev)%n1
                 ltmp(-2,i,j,k) = ds2(i,j-1,k) * 2.0/(d2(i,j,k)+d2(i,j-1,k)) * dm(i,j,k)*dm(i,j-1,k)
                 ltmp( 2,i,j,k) = ds2(i,j,  k) * 2.0/(d2(i,j,k)+d2(i,j+1,k)) * dm(i,j,k)*dm(i,j+1,k)
              END DO
              END DO
              END DO
           CASE (2)
              DO k=1, mg(lev)%n3
              DO j=1, mg(lev)%n2
              DO i=1, mg(lev)%n1
                 ltmp(-2,i,j,k) = ds2(i,j-1,k) * d1(i,j-1,k)*d3(i,j-1,k) * 2.0/(dv(i,j,k)+dv(i,j-1,k)) * dm(i,j,k)*dm(i,j-1,k)
                 ltmp( 2,i,j,k) = ds2(i,j,  k) * d1(i,j+1,k)*d3(i,j+1,k) * 2.0/(dv(i,j,k)+dv(i,j+1,k)) * dm(i,j,k)*dm(i,j+1,k)
              END DO
              END DO
              END DO
           CASE (3)
              DO k=1, mg(lev)%n3
              DO j=1, mg(lev)%n2
              DO i=1, mg(lev)%n1
                 ltmp(-2,i,j,k) = ds2(i,j-1,k) / d2(i,j-1,k) * 2.0/(dv(i,j,k)+dv(i,j-1,k)) * dm(i,j,k)*dm(i,j-1,k)
                 ltmp( 2,i,j,k) = ds2(i,j,  k) / d2(i,j+1,k) * 2.0/(dv(i,j,k)+dv(i,j+1,k)) * dm(i,j,k)*dm(i,j+1,k)
              END DO
              END DO
              END DO
           END SELECT
        END IF

        SELECT CASE (coeffscheme)
        CASE (1)
           DO k=1, mg(lev)%n3
           DO j=1, mg(lev)%n2
           DO i=1, mg(lev)%n1
              ltmp(-3,i,j,k) = ds3(i,j,k-1) * 2.0/(d3(i,j,k)+d3(i,j,k-1)) * dm(i,j,k)*dm(i,j,k-1)
              ltmp( 3,i,j,k) = ds3(i,j,k  ) * 2.0/(d3(i,j,k)+d3(i,j,k+1)) * dm(i,j,k)*dm(i,j,k+1)
           END DO
           END DO
           END DO
        CASE (2)
           DO k=1, mg(lev)%n3
           DO j=1, mg(lev)%n2
           DO i=1, mg(lev)%n1
              ltmp(-3,i,j,k) = ds3(i,j,k-1) * d1(i,j,k-1)*d2(i,j,k-1) * 2.0/(dv(i,j,k)+dv(i,j,k-1)) * dm(i,j,k)*dm(i,j,k-1)
              ltmp( 3,i,j,k) = ds3(i,j,k  ) * d1(i,j,k+1)*d2(i,j,k+1) * 2.0/(dv(i,j,k)+dv(i,j,k+1)) * dm(i,j,k)*dm(i,j,k+1)
           END DO
           END DO
           END DO
        CASE (3)
           DO k=1, mg(lev)%n3
           DO j=1, mg(lev)%n2
           DO i=1, mg(lev)%n1
              ltmp(-3,i,j,k) = ds3(i,j,k-1) / d3(i,j,k-1) * 2.0/(dv(i,j,k)+dv(i,j,k-1)) * dm(i,j,k)*dm(i,j,k-1)
              ltmp( 3,i,j,k) = ds3(i,j,k  ) / d3(i,j,k+1) * 2.0/(dv(i,j,k)+dv(i,j,k+1)) * dm(i,j,k)*dm(i,j,k+1)
           END DO
           END DO
           END DO
        END SELECT

        SELECT CASE (coeffscheme)
        CASE (1)
           DO k=1, mg(lev)%n3
           DO j=1, mg(lev)%n2
           DO i=1, mg(lev)%n1
              ltmp( 0,i,j,k) = -da(i,j,k) - ltmp(-1,i,j,k) - ltmp(1,i,j,k) &
                                          - ltmp(-2,i,j,k) - ltmp(2,i,j,k) &
                                          - ltmp(-3,i,j,k) - ltmp(3,i,j,k)
           END DO
           END DO
           END DO
        CASE (2)
           DO k=1, mg(lev)%n3
           DO j=1, mg(lev)%n2
           DO i=1, mg(lev)%n1
              ltmp( 0,i,j,k) = -da(i,j,k) - ds1(i-1,j,k) * d2(i,j,k)*d3(i,j,k) * 2.0/(dv(i,j,k)+dv(i-1,j,k)) * dm(i,j,k)*dm(i-1,j,k) &
                                          - ds1(i,  j,k) * d2(i,j,k)*d3(i,j,k) * 2.0/(dv(i,j,k)+dv(i+1,j,k)) * dm(i,j,k)*dm(i+1,j,k) &
                                          - ds2(i,j-1,k) * d1(i,j,k)*d3(i,j,k) * 2.0/(dv(i,j,k)+dv(i,j-1,k)) * dm(i,j,k)*dm(i,j-1,k) &
                                          - ds2(i,j,  1) * d1(i,j,k)*d3(i,j,k) * 2.0/(dv(i,j,k)+dv(i,j+1,k)) * dm(i,j,k)*dm(i,j+1,k) &
                                          - ds3(i,j,k-1) * d1(i,j,k)*d2(i,j,k) * 2.0/(dv(i,j,k)+dv(i,j,k-1)) * dm(i,j,k)*dm(i,j,k-1) &
                                          - ds3(i,j,k+1) * d1(i,j,k)*d2(i,j,k) * 2.0/(dv(i,j,k)+dv(i,j,k+1)) * dm(i,j,k)*dm(i,j,k+1)
           END DO
           END DO
           END DO
        CASE (3)
           DO k=1, mg(lev)%n3
           DO j=1, mg(lev)%n2
           DO i=1, mg(lev)%n1
              ltmp( 0,i,j,k) = -da(i,j,k) - ds1(i-1,j,k) / d1(i,j,k) * 2.0/(dv(i,j,k)+dv(i-1,j,k)) * dm(i,j,k)*dm(i-1,j,k) &
                                          - ds1(i,  j,k) / d1(i,j,k) * 2.0/(dv(i,j,k)+dv(i+1,j,k)) * dm(i,j,k)*dm(i+1,j,k) &
                                          - ds2(i,j-1,k) / d2(i,j,k) * 2.0/(dv(i,j,k)+dv(i,j-1,k)) * dm(i,j,k)*dm(i,j-1,k) &
                                          - ds2(i,j,  1) / d2(i,j,k) * 2.0/(dv(i,j,k)+dv(i,j+1,k)) * dm(i,j,k)*dm(i,j+1,k) &
                                          - ds3(i,j,k-1) / d3(i,j,k) * 2.0/(dv(i,j,k)+dv(i,j,k-1)) * dm(i,j,k)*dm(i,j,k-1) &
                                          - ds3(i,j,k+1) / d3(i,j,k) * 2.0/(dv(i,j,k)+dv(i,j,k+1)) * dm(i,j,k)*dm(i,j,k+1)
           END DO
           END DO
           END DO
        END SELECT

        CALL sai3d(lev, ltmp, mtmp)

        ALLOCATE(mg(lev)%l(-3:3,1:mg(lev)%n1,1:mg(lev)%n2,1:mg(lev)%n3))
        ALLOCATE(mg(lev)%m(-3:3,1:mg(lev)%n1,1:mg(lev)%n2,1:mg(lev)%n3))

        mg(lev)%l(:,:,:,:) = real(ltmp(:,:,:,:), KIND=MG_KIND)
        mg(lev)%m(:,:,:,:) = real(mtmp(:,:,:,:), KIND=MG_KIND)

        lflag(1:mg(lev)%n3) = .FALSE.
        mflag(1:mg(lev)%n3) = .FALSE.

        DO k=2, mg(lev)%n3
           lflag(k) = .TRUE.
           mflag(k) = .TRUE.
           DO j=1, mg(lev)%n2
           DO i=1, mg(lev)%n1
              lflag(k) = lflag(k) &
                   .AND. ltmp(-3,i,j,k) == ltmp(-3,i,j,k-1) &
                   .AND. ltmp(-2,i,j,k) == ltmp(-2,i,j,k-1) &
                   .AND. ltmp(-1,i,j,k) == ltmp(-1,i,j,k-1) &
                   .AND. ltmp( 0,i,j,k) == ltmp( 0,i,j,k-1) &
                   .AND. ltmp( 1,i,j,k) == ltmp( 1,i,j,k-1) &
                   .AND. ltmp( 2,i,j,k) == ltmp( 2,i,j,k-1) &
                   .AND. ltmp( 3,i,j,k) == ltmp( 3,i,j,k-1)

              mflag(k) = mflag(k) &
                   .AND. mtmp(-3,i,j,k) == mtmp(-3,i,j,k-1) &
                   .AND. mtmp(-2,i,j,k) == mtmp(-2,i,j,k-1) &
                   .AND. mtmp(-1,i,j,k) == mtmp(-1,i,j,k-1) &
                   .AND. mtmp( 0,i,j,k) == mtmp( 0,i,j,k-1) &
                   .AND. mtmp( 1,i,j,k) == mtmp( 1,i,j,k-1) &
                   .AND. mtmp( 2,i,j,k) == mtmp( 2,i,j,k-1) &
                   .AND. mtmp( 3,i,j,k) == mtmp( 3,i,j,k-1)
           END DO
           END DO
        END DO

        ALLOCATE(mg(lev)%lindx(1:mg(lev)%n3))
        ALLOCATE(mg(lev)%mindx(1:mg(lev)%n3))

        mg(lev)%lindx(1) = 1
        mg(lev)%mindx(1) = 1
        DO k=2, mg(lev)%n3
           IF (lflag(k)) THEN
              mg(lev)%lindx(k) = mg(lev)%lindx(k-1)
           ELSE
              mg(lev)%lindx(k) = k
           END IF

           IF (mflag(k)) THEN
              mg(lev)%mindx(k) = mg(lev)%mindx(k-1)
           ELSE
              mg(lev)%mindx(k) = k
           END IF
        END DO

        ALLOCATE(mg(lev)%r(0:mg(lev)%n1+1, 0:mg(lev)%n2+1, 0:mg(lev)%n3+1))
        ALLOCATE(mg(lev)%e(0:mg(lev)%n1+1, 0:mg(lev)%n2+1, 0:mg(lev)%n3+1))
        ALLOCATE(mg(lev)%t(0:mg(lev)%n1+1, 0:mg(lev)%n2+1, 0:mg(lev)%n3+1))
     END IF

     IF (lev==0) THEN
        ALLOCATE(l3d(-3:3,1:n1,1:n2,1:n3))

        l3d(:,:,:,:) = ltmp(:,:,:,:)
     END IF

     IF (ALLOCATED(ltmp)) DEALLOCATE(ltmp)
     IF (ALLOCATED(mtmp)) DEALLOCATE(mtmp)

     mg(lev)%r(:,:,:) = 0.0D0
     mg(lev)%e(:,:,:) = 0.0D0
     mg(lev)%t(:,:,:) = 0.0D0

     ALLOCATE(mg(lev)%p1(0:mg(lev)%n1))
     ALLOCATE(mg(lev)%p2(0:mg(lev)%n2))
     ALLOCATE(mg(lev)%p3(0:mg(lev)%n3))

     ALLOCATE(mg(lev)%q1(1:mg(lev)%n1))
     ALLOCATE(mg(lev)%q2(1:mg(lev)%n2))
     ALLOCATE(mg(lev)%q3(1:mg(lev)%n3))

#ifdef PARALLEL_MPI
     IF (kpes > 1 .AND. vlev > lev) THEN
        CALL mpi_allreduce(mg(lev)%n3, n3_max, 1, MPI_INTEGER, MPI_MAX, comm, ierr)
        IF (n3_max == 1) THEN
           vlev = lev
        ENDIF
     ELSE IF (ipes*jpes > 1 .AND. glev > lev) THEN
        IF (mg(lev)%n1 <= MG_GDIM .AND. mg(lev)%n2 <= MG_GDIM .AND. mg(lev)%n3 == 1) THEN
           glev = lev
           CALL assert(mg(glev)%n1 == MG_GDIM .AND. mg(glev)%n2 == MG_GDIM, "MG_GDIM does not fit")
        ENDIF
!     IF (glev > lev) THEN
!        CALL mpi_allreduce(mg(lev)%n3, n3_max, 1, MPI_INTEGER, MPI_MAX, comm, ierr)
!        IF (mg(lev)%n1 == 1 .AND. mg(lev)%n2 == 1 .AND. n3_max == 1) THEN
!           glev = lev
!        ENDIF
     ELSE IF (mg(lev)%n1 <= min_dim .AND. mg(lev)%n2 <= min_dim .AND. mg(lev)%n3==1) THEN
           nlev = lev
           EXIT
     END IF

     IF (lev > vlev .AND. vrank /= 0) EXIT
     IF (lev > glev .AND. hrank /= 0) EXIT
#else
     IF (mg(lev)%n1 <= min_dim .AND. mg(lev)%n2 <= min_dim .AND. mg(lev)%n3==1) THEN
        nlev = lev
        EXIT
     END IF
#endif
  END DO

#ifdef PARALLEL_MPI
  CALL mpi_barrier(comm, ierr)
#endif

  IF (rank==0) THEN
     singular=(maxval(da(1:mg(nlev)%n1,1:mg(nlev)%n2,1))==0.0D0 .AND. &
               minval(da(1:mg(nlev)%n1,1:mg(nlev)%n2,1))==0.0D0)
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
        CALL assert(.FALSE., "invalid solver method '"//trim(method)//".")
     END SELECT
  END IF

  IF (num_threads > 1 .AND. (solver_method == 1 .OR. solver_method == 2)) solver_method = -solver_method

  IF (present(epsilon)) eps   = epsilon
  IF (present(c_mode))  cmode = c_mode
  IF (present(it_max))  itmax = it_max
  IF (present(report))  itreport = report
  IF (present(scheme))  coeffscheme = scheme

  CALL assert(coeffscheme==1 .OR. coeffscheme==2 .OR. coeffscheme==3, "unsupported coeffscheme")

  IF (rank == 0) THEN
     WRITE(REPORT_UNIT, '(A)', ADVANCE='NO') "solver3d using "
     SELECT CASE (solver_method)
     CASE (0)
        WRITE(REPORT_UNIT, '(A)', ADVANCE='NO') "Pure CG method"
     CASE (1)
        WRITE(REPORT_UNIT, '(A)', ADVANCE='NO') "Stand-alone MG method"
     CASE (-1)
        WRITE(REPORT_UNIT, '(A)', ADVANCE='NO') "Stand-alone MG method (OpenMP)"
     CASE (2)
        WRITE(REPORT_UNIT, '(A)', ADVANCE='NO') "MG-CG method"
     CASE (-2)
        WRITE(REPORT_UNIT, '(A)', ADVANCE='NO') "MG-CG method (OpenMP)"
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
     WRITE(REPORT_UNIT, '(A,ES10.3,A,I0)') ", eps =", eps, ", itmax=", itmax
  END IF

  DEALLOCATE(d1)
  DEALLOCATE(d2)
  DEALLOCATE(d3)

  DEALLOCATE(ds1)
  DEALLOCATE(ds2)
  DEALLOCATE(ds3)

  DEALLOCATE(dv)
  DEALLOCATE(da)

  DEALLOCATE(dm)

  DEALLOCATE(lflag)
  DEALLOCATE(mflag)

!  use_il = .TRUE.
  IF (rank==0) CALL init_il

#ifdef PARALLEL_MPI
  CALL mpi_allreduce(MPI_IN_PLACE, glev, 1, MPI_INTEGER, MPI_MIN, comm, ierr)
  CALL mpi_allreduce(MPI_IN_PLACE, nlev, 1, MPI_INTEGER, MPI_MIN, comm, ierr)
#endif

  DO lev=0, nlev
#ifdef PARALLEL_MPI
     CALL mpi_barrier(comm, ierr)
#endif

     IF (num_threads > mg(lev)%n3) THEN
        DO n=0, mg(lev)%n3-1
           n3_start(n,lev) = n+1
           n3_span( n,lev) = 1
        END DO
        n3_start(mg(lev)%n3:num_threads-1,lev) = 0
        n3_span( mg(lev)%n3:num_threads-1,lev) = 0
     ELSE
        n3_span(0:num_threads-1,lev) = mg(lev)%n3/num_threads

        DO n=1, mg(lev)%n3 - n3_span(0,lev)*num_threads
           n3_span(num_threads-n,lev) = n3_span(num_threads-n,lev) + 1
        END DO
        DO n=0, num_threads-1
           n3_start(n,lev) = mg(lev)%n3 - sum(n3_span(n:num_threads-1,lev)) + 1
        END DO
        CALL assert(sum(n3_span(0:num_threads-1,lev))==mg(lev)%n3, "bad thread-assignment in solver3d")
     END IF

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
        CALL assert(sum(n2_span(0:num_threads-1,lev))==mg(lev)%n2, "bad thread-assignment in solver3d")
     END IF

#ifdef DEBUG
!$   IF (rank==0) THEN
!$      WRITE(REPORT_UNIT, *) "OMP_NUM_THREADS=", num_threads
!$      WRITE(REPORT_UNIT, '("thread-assignment for lev=", I2)') lev
!$      DO n=0, num_threads-1
!$         WRITE(REPORT_UNIT,'("#", I2,"   k=", I3, ":", I3, "  j=", I3, ":", I3)') &
!$            n, n3_start(n,lev), n3_start(n,lev)+n3_span(n,lev)-1,  &
!$               n2_start(n,lev), n2_start(n,lev)+n2_span(n,lev)-1
!$      END DO
!$   END IF
#endif
  END DO

  CALL reset_solver3d(dx,   dy,   dz,   dsx,   dsy,   dsz,   dvol,   alpha,   mask,   &
                      dx_h, dy_h, dz_h, dsx_h, dsy_h, dsz_h, dvol_h, alpha_h, mask_h, &
                      dx_v, dy_v, dz_v, dsx_v, dsy_v, dsz_v, dvol_v)

END SUBROUTINE init_solver3d

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE reset_solver3d(dx,   dy,   dz,   dsx,   dsy,   dsz,   dvol,   alpha,   mask,   &
                          dx_h, dy_h, dz_h, dsx_h, dsy_h, dsz_h, dvol_h, alpha_h, mask_h, &
                          dx_v, dy_v, dz_v, dsx_v, dsy_v, dsz_v, dvol_v)
  REAL(8), INTENT(IN), OPTIONAL :: dx(:,:,:), dx_h(:,:), dx_v(:)
  REAL(8), INTENT(IN), OPTIONAL :: dy(:,:,:), dy_h(:,:), dy_v(:)
  REAL(8), INTENT(IN), OPTIONAL :: dz(:,:,:), dz_h(:,:), dz_v(:)

  REAL(8), INTENT(IN), OPTIONAL :: dsx(:,:,:), dsx_h(:,:), dsx_v(:)
  REAL(8), INTENT(IN), OPTIONAL :: dsy(:,:,:), dsy_h(:,:), dsy_v(:)
  REAL(8), INTENT(IN), OPTIONAL :: dsz(:,:,:), dsz_h(:,:), dsz_v(:)

  REAL(8), INTENT(IN), OPTIONAL :: dvol(:,:,:), dvol_h(:,:), dvol_v(:)

  REAL(8), INTENT(IN), OPTIONAL :: alpha(:,:,:), alpha_h(:,:)

  LOGICAL(1), INTENT(IN), OPTIONAL :: mask(:,:,:), mask_h(:,:)

  REAL(8) :: d1(0:n1+1,0:n2+1,0:n3+1)
  REAL(8) :: d2(0:n1+1,0:n2+1,0:n3+1)
  REAL(8) :: d3(0:n1+1,0:n2+1,0:n3+1)

  REAL(8) :: ds1(0:n1+1,0:n2+1,0:n3+1)
  REAL(8) :: ds2(0:n1+1,0:n2+1,0:n3+1)
  REAL(8) :: ds3(0:n1+1,0:n2+1,0:n3+1)

  REAL(8) :: dv(0:n1+1,0:n2+1,0:n3+1)
  REAL(8) :: da(0:n1+1,0:n2+1,0:n3+1)

  REAL(8) :: dm(0:n1+1,0:n2+1,0:n3+1)

  INTEGER :: i, j, k

#ifdef PARALLEL_MPI
  INTEGER :: ierr

  IF (comm == MPI_COMM_NULL) RETURN
#endif

  CALL assert(initialized, "solver3d is not initialized.")
  CALL assert(n3>1, "solver dimension missmatch")

  CALL setparam3d(n1, n2, n3, d1, dx, dx_h, dx_v)
  CALL setparam3d(n1, n2, n3, d2, dy, dy_h, dy_v)
  CALL setparam3d(n1, n2, n3, d3, dz, dz_h, dz_v)

  CALL setparam3d(n1, n2, n3, ds1, dsx, dsx_h, dsx_v)
  CALL setparam3d(n1, n2, n3, ds2, dsy, dsy_h, dsy_v)
  CALL setparam3d(n1, n2, n3, ds3, dsz, dsz_h, dsz_v)

  CALL setparam3d(n1, n2, n3, dv, dvol, dvol_h, dvol_v)

  CALL setparam3d(n1, n2, n3, da, alpha, alpha_h, d=0.0D0)

  CALL setparam3d_logical(n1, n2, n3, dm, mask, mask_h, d=LOGICAL(.TRUE.,1))

  DO k=1, n3
!$OMP PARALLEL PRIVATE(i, j)

#ifdef PARALLEL_MPI
     IF (n1 == 1 .AND. mg(0)%rank_e == rank .AND. mg(0)%rank_w == rank) THEN
#else
     IF (n1 == 1) THEN
#endif
        l3d(-1,1,1:n2,k) = 0.0D0
        l3d( 1,1,1:n2,k) = 0.0D0
     ELSE
        SELECT CASE (coeffscheme)
        CASE (1)
!$OMP DO
           DO j=1, n2
           DO i=1, n1
              l3d(-1,i,j,k) = ds1(i-1,j,k) * 2.0/(d1(i,j,k)+d1(i-1,j,k)) * dm(i,j,k)*dm(i-1,j,k)
              l3d( 1,i,j,k) = ds1(i,  j,k) * 2.0/(d1(i,j,k)+d1(i+1,j,k)) * dm(i,j,k)*dm(i+1,j,k)
           END DO
           END DO

        CASE (2)
!$OMP DO
           DO j=1, n2
           DO i=1, n1
              l3d(-1,i,j,k) = dm(i,j,k)*dm(i-1,j,k) * ds1(i-1,j,k) * d2(i-1,j,k)*d3(i-1,j,k) * 2.0/(dv(i,j,k)+dv(i-1,j,k))
              l3d( 1,i,j,k) = dm(i,j,k)*dm(i+1,j,k) * ds1(i,  j,k) * d2(i+1,j,k)*d3(i+1,j,k) * 2.0/(dv(i,j,k)+dv(i+1,j,k))
           END DO
           END DO

        CASE (3)
!$OMP DO
           DO j=1, n2
           DO i=1, n1
              l3d(-1,i,j,k) = dm(i,j,k)*dm(i-1,j,k) * ds1(i-1,j,k) / d1(i-1,j,k) * 2.0/(dv(i,j,k)+dv(i-1,j,k))
              l3d( 1,i,j,k) = dm(i,j,k)*dm(i+1,j,k) * ds1(i,  j,k) / d1(i+1,j,k) * 2.0/(dv(i,j,k)+dv(i+1,j,k))
           END DO
           END DO
        END SELECT
     END IF

#ifdef PARALLEL_MPI
     IF (n2 == 1 .AND. mg(0)%rank_n == rank .AND. mg(0)%rank_s == rank) THEN
#else
     IF (n2 == 1) THEN
#endif
        l3d(-2,1:n1,1,k) = 0.0D0
        l3d( 2,1:n1,1,k) = 0.0D0
     ELSE
        SELECT CASE (coeffscheme)
        CASE (1)
!$OMP DO
           DO j=1, n2
           DO i=1, n1
              l3d(-2,i,j,k) = ds2(i,j-1,k) * 2.0/(d2(i,j,k)+d2(i,j-1,k)) * dm(i,j,k)*dm(i,j-1,k)
              l3d( 2,i,j,k) = ds2(i,j,  k) * 2.0/(d2(i,j,k)+d2(i,j+1,k)) * dm(i,j,k)*dm(i,j+1,k)
           END DO
           END DO
        CASE (2)
!$OMP DO
           DO j=1, n2
           DO i=1, n1
              l3d(-2,i,j,k) = dm(i,j,k)*dm(i,j-1,k) * ds2(i,j-1,k) * d1(i,j-1,k)*d3(i,j-1,k) * 2.0/(dv(i,j,k)+dv(i,j-1,k))
              l3d( 2,i,j,k) = dm(i,j,k)*dm(i,j+1,k) * ds2(i,j,  k) * d1(i,j+1,k)*d3(i,j+1,k) * 2.0/(dv(i,j,k)+dv(i,j+1,k))
           END DO
           END DO
        CASE (3)
!$OMP DO
           DO j=1, n2
           DO i=1, n1
              l3d(-2,i,j,k) = dm(i,j,k)*dm(i,j-1,k) * ds2(i,j-1,k) / d2(i,j-1,k) * 2.0/(dv(i,j,k)+dv(i,j-1,k))
              l3d( 2,i,j,k) = dm(i,j,k)*dm(i,j+1,k) * ds2(i,j,  k) / d2(i,j+1,k) * 2.0/(dv(i,j,k)+dv(i,j+1,k))
           END DO
           END DO
        END SELECT
     END IF

     SELECT CASE (coeffscheme)
     CASE (1)
!$OMP DO
        DO j=1, n2
        DO i=1, n1
           l3d(-3,i,j,k) = ds3(i,j,k-1) * 2.0/(d3(i,j,k)+d3(i,j,k-1)) * dm(i,j,k)*dm(i,j,k-1)
           l3d( 3,i,j,k) = ds3(i,j,k)   * 2.0/(d3(i,j,k)+d3(i,j,k+1)) * dm(i,j,k)*dm(i,j,k+1)
        END DO
        END DO
     CASE (2)
!$OMP DO
        DO j=1, n2
        DO i=1, n1
           l3d(-3,i,j,k) = dm(i,j,k)*dm(i,j,k-1) * ds3(i,j,k-1) * d1(i,j,k-1)*d2(i,j,k-1) * 2.0/(dv(i,j,k)+dv(i,j,k-1))
           l3d( 3,i,j,k) = dm(i,j,k)*dm(i,j,k+1) * ds3(i,j,k)   * d1(i,j,k+1)*d2(i,j,k+1) * 2.0/(dv(i,j,k)+dv(i,j,k+1))
        END DO
        END DO
     CASE (3)
!$OMP DO
        DO j=1, n2
        DO i=1, n1
           l3d(-3,i,j,k) = dm(i,j,k)*dm(i,j,k-1) * ds3(i,j,k-1) / d3(i,j,k-1) * 2.0/(dv(i,j,k)+dv(i,j,k-1))
           l3d( 3,i,j,k) = dm(i,j,k)*dm(i,j,k+1) * ds3(i,j,k)   / d3(i,j,k+1) * 2.0/(dv(i,j,k)+dv(i,j,k+1))
        END DO
        END DO
     END SELECT

     SELECT CASE (coeffscheme)
     CASE (1)
!$OMP DO
        DO j=1, n2
        DO i=1, n1
           l3d(0,i,j,k) = -da(i,j,k) - l3d(-1,i,j,k) - l3d(1,i,j,k) &
                                     - l3d(-2,i,j,k) - l3d(2,i,j,k) &
                                     - l3d(-3,i,j,k) - l3d(3,i,j,k)
        END DO
        END DO
     CASE (2)
!$OMP DO
        DO j=1, n2
        DO i=1, n1
           l3d(0,i,j,k) = -da(i,j,k) - dm(i,j,k)*dm(i-1,j,k) * ds1(i-1,j,k) * d2(i,j,k)*d3(i,j,k) * 2.0/(dv(i,j,k)+dv(i-1,j,k)) &
                                     - dm(i,j,k)*dm(i+1,j,k) * ds1(i,  j,k) * d2(i,j,k)*d3(i,j,k) * 2.0/(dv(i,j,k)+dv(i+1,j,k)) &
                                     - dm(i,j,k)*dm(i,j-1,k) * ds2(i,j-1,k) * d1(i,j,k)*d3(i,j,k) * 2.0/(dv(i,j,k)+dv(i,j-1,k)) &
                                     - dm(i,j,k)*dm(i,j+1,k) * ds2(i,j,  k) * d1(i,j,k)*d3(i,j,k) * 2.0/(dv(i,j,k)+dv(i,j+1,k)) &
                                     - dm(i,j,k)*dm(i,j,k-1) * ds3(i,j,k-1) * d1(i,j,k)*d2(i,j,k) * 2.0/(dv(i,j,k)+dv(i,j,k-1)) &
                                     - dm(i,j,k)*dm(i,j,k+1) * ds3(i,j,k)   * d1(i,j,k)*d2(i,j,k) * 2.0/(dv(i,j,k)+dv(i,j,k+1))
        END DO
        END DO
     CASE (3)
!$OMP DO
        DO j=1, n2
        DO i=1, n1
           l3d( 0,i,j,k) = -da(i,j,k) - dm(i,j,k)*dm(i-1,j,k) * ds1(i-1,j,k) / d1(i,j,k) * 2.0/(dv(i,j,k)+dv(i-1,j,k)) &
                                      - dm(i,j,k)*dm(i+1,j,k) * ds1(i,  j,k) / d1(i,j,k) * 2.0/(dv(i,j,k)+dv(i+1,j,k)) &
                                      - dm(i,j,k)*dm(i,j-1,k) * ds2(i,j-1,k) / d2(i,j,k) * 2.0/(dv(i,j,k)+dv(i,j-1,k)) &
                                      - dm(i,j,k)*dm(i,j+1,k) * ds2(i,j,  k) / d2(i,j,k) * 2.0/(dv(i,j,k)+dv(i,j+1,k)) &
                                      - dm(i,j,k)*dm(i,j,k-1) * ds3(i,j,k-1) / d3(i,j,k) * 2.0/(dv(i,j,k)+dv(i,j,k-1)) &
                                      - dm(i,j,k)*dm(i,j,k+1) * ds3(i,j,k)   / d3(i,j,k) * 2.0/(dv(i,j,k)+dv(i,j,k+1))
        END DO
        END DO
     END SELECT

!$OMP DO
     DO j=1, n2
     DO i=1, n1
        mg(0)%l(:,i,j,k) = real(l3d(:,i,j,k), KIND=MG_KIND)
     END DO
     END DO
!$OMP END PARALLEL

     IF (k>1) THEN
        DO j=1, n2
           IF (mg(0)%lindx(k)==k) EXIT
           DO i=1, n1
              IF (mg(0)%l(0,i,j,k) /= mg(0)%l(0,i,j,k-1)) mg(0)%lindx(k) = k
              IF (mg(0)%lindx(k)==k) EXIT
           END DO
        END DO
     END IF
  END DO

END SUBROUTINE reset_solver3d

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

  INTEGER :: n1, n2

  n1 = mg(lev)%n1
  n2 = mg(lev)%n2

  tmp(:,:,:) = 0.0D0
  tmp(:,1:n1,1:n2) = l(:,1:n1,1:n2)

#ifdef PARALLEL_MPI
  IF (lev > glev) THEN
     IF (hrank /= 0) RETURN

     IF (period1) THEN
        tmp(:,0,   :) = tmp(:,n1,:)
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
                       comm, MPI_STATUS_IGNORE, ierr)
     tmp(:,0,:) = recvbuf1(:,:)

     sendbuf1(:,:) = tmp(:,1,:)
     recvbuf1(:,:) = 0.0D0
     CALL mpi_sendrecv(sendbuf1(-2,0), (n2+2)*5, MPI_REAL8, mg(lev)%rank_w, 2, &
                       recvbuf1(-2,0), (n2+2)*5, MPI_REAL8, mg(lev)%rank_e, 2, &
                       comm, MPI_STATUS_IGNORE, ierr)
     tmp(:,n1+1,:) = recvbuf1(:,:)

     sendbuf2(:,:) = tmp(:,:,n2)
     recvbuf2(:,:) = 0.0D0
     CALL mpi_sendrecv(sendbuf2(-2,0), (n1+2)*5, MPI_REAL8, mg(lev)%rank_n, 3, &
                       recvbuf2(-2,0), (n1+2)*5, MPI_REAL8, mg(lev)%rank_s, 3, &
                       comm, MPI_STATUS_IGNORE, ierr)
     tmp(:,:,0) = recvbuf2(:,:)

     sendbuf2(:,:) = tmp(:,:,1)
     recvbuf2(:,:) = 0.0D0
     CALL mpi_sendrecv(sendbuf2(-2,0), (n1+2)*5, MPI_REAL8, mg(lev)%rank_s, 4, &
                       recvbuf2(-2,0), (n1+2)*5, MPI_REAL8, mg(lev)%rank_n, 4, &
                       comm, MPI_STATUS_IGNORE, ierr)
     tmp(:,:,n2+1) = recvbuf2(:,:)
END IF

#else
  IF (period1) THEN
     tmp(:,0,:)    = tmp(:,n1,:)
     tmp(:,n1+1,:) = tmp(:,1,:)
  END IF

  IF (period2) THEN
     tmp(:,:,0)    = tmp(:,:,n2)
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
        DO p=-2, 2
           IF (aa(p,p) == 0.0D0) aa(p,p) = 1.0D0
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

        DO p=-2, 2
           m(p,i,j) = tmp(0,i,j)*aa(p,0) &
                    + tmp(-1,i+1,j)*aa(p,1) + tmp(1,i-1,j)*aa(p,-1) &
                    + tmp(-2,i,j+1)*aa(p,2) + tmp(2,i,j-1)*aa(p,-2)
        END DO
#endif
     END DO
  END DO
END SUBROUTINE sai2d

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE sai3d(lev, l, m)
  INTEGER, INTENT(IN)  :: lev
  REAL(8), INTENT(IN)  :: l(-3:3, 1:mg(lev)%n1, 1:mg(lev)%n2, 1:mg(lev)%n3)
  REAL(8), INTENT(OUT) :: m(-3:3, 1:mg(lev)%n1, 1:mg(lev)%n2, 1:mg(lev)%n3)

  REAL(8) :: a( -12:12, -3:3)
  REAL(8) :: aa(-12:12, -3:3)

  REAL(8) :: tmp(-3:3,0:mg(lev)%n1+1,0:mg(lev)%n2+1,0:mg(lev)%n3+1)

#ifdef PARALLEL_MPI
  REAL(8) :: sendbuf1(-3:3,0:mg(lev)%n2+1,0:mg(lev)%n3+1)
  REAL(8) :: recvbuf1(-3:3,0:mg(lev)%n2+1,0:mg(lev)%n3+1)
  REAL(8) :: sendbuf2(-3:3,0:mg(lev)%n1+1,0:mg(lev)%n3+1)
  REAL(8) :: recvbuf2(-3:3,0:mg(lev)%n1+1,0:mg(lev)%n3+1)
#ifdef PARALLEL3D
  REAL(8) :: sendbuf3(-3:3,0:mg(lev)%n1+1,0:mg(lev)%n2+1)
  REAL(8) :: recvbuf3(-3:3,0:mg(lev)%n1+1,0:mg(lev)%n2+1)
#endif

  INTEGER :: ierr
#endif

  INTEGER :: i, j, k
  INTEGER :: p, q

  INTEGER :: invmat_i
  INTEGER :: invmat_j
  REAL(8) :: invmat_p
  REAL(8) :: invmat_q

  INTEGER :: n1, n2, n3

  n1 = mg(lev)%n1
  n2 = mg(lev)%n2
  n3 = mg(lev)%n3

  tmp(:,:,:,:) = 0.0D0
  tmp(-3:3,1:n1,1:n2,1:n3) = l(-3:3,1:n1,1:n2,1:n3)

#ifdef PARALLEL_MPI
  IF (lev > glev) THEN
     IF (period1) THEN
        tmp(:,0,:,:)    = tmp(:,n1,:,:)
        tmp(:,n1+1,:,:) = tmp(:,1,:,:)
     END IF

     IF (period2) THEN
        tmp(:,:,0,:)    = tmp(:,:,n2,:)
        tmp(:,:,n2+1,:) = tmp(:,:,1,:)
     END IF

     IF (period3) THEN
        tmp(:,:,:,0)    = tmp(:,:,:,n3)
        tmp(:,:,:,n3+1) = tmp(:,:,:,1)
     END IF
  ELSE
     sendbuf1(:,:,:) = tmp(-3:3,n1,:,:)
     recvbuf1(:,:,:) = 0.0D0

     CALL mpi_sendrecv(sendbuf1(-3,0,0), (n2+2)*(n3+2)*7, MPI_REAL8, mg(lev)%rank_e, 1, &
                       recvbuf1(-3,0,0), (n2+2)*(n3+2)*7, MPI_REAL8, mg(lev)%rank_w, 1, &
                       comm, MPI_STATUS_IGNORE, ierr)

     tmp(:,0,:,:) = recvbuf1(:,:,:)

     sendbuf1(:,:,:) = tmp(:,1,:,:)
     recvbuf1(:,:,:) = 0.0D0
     CALL mpi_sendrecv(sendbuf1(-3,0,0), (n2+2)*(n3+2)*7, MPI_REAL8, mg(lev)%rank_w, 2, &
                       recvbuf1(-3,0,0), (n2+2)*(n3+2)*7, MPI_REAL8, mg(lev)%rank_e, 2, &
                       comm, MPI_STATUS_IGNORE, ierr)

     tmp(:,n1+1,:,:) = recvbuf1(:,:,:)

     sendbuf2(:,:,:) = tmp(:,:,n2,:)
     recvbuf2(:,:,:) = 0.0D0
     CALL mpi_sendrecv(sendbuf2(-3,0,0), (n1+2)*(n3+2)*7, MPI_REAL8, mg(lev)%rank_n, 3, &
                       recvbuf2(-3,0,0), (n1+2)*(n3+2)*7, MPI_REAL8, mg(lev)%rank_s, 3, &
                       comm, MPI_STATUS_IGNORE, ierr)

     tmp(:,:,0,:) = recvbuf2(:,:,:)

     sendbuf2(:,:,:) = tmp(:,:,1,:)
     recvbuf2(:,:,:) = 0.0D0
     CALL mpi_sendrecv(sendbuf2(-3,0,0), (n1+2)*(n3+2)*7, MPI_REAL8, mg(lev)%rank_s, 4, &
                       recvbuf2(-3,0,0), (n1+2)*(n3+2)*7, MPI_REAL8, mg(lev)%rank_n, 4, &
                       comm, MPI_STATUS_IGNORE, ierr)

     tmp(:,:,n2+1,:) = recvbuf2(:,:,:)

#ifdef PARALLEL3D
     IF (lev <= vlev) THEN
        sendbuf3(:,:,:) = tmp(:,:,:,n3)
        recvbuf3(:,:,:) = 0.0D0
        CALL mpi_sendrecv(sendbuf3(-3,0,0), (n1+2)*(n2+2)*7, MPI_REAL8, mg(lev)%rank_u, 5, &
                          recvbuf3(-3,0,0), (n1+2)*(n2+2)*7, MPI_REAL8, mg(lev)%rank_l, 5, &
                          comm, MPI_STATUS_IGNORE, ierr)

        tmp(:,:,:,0) = recvbuf3(:,:,:)

        sendbuf3(:,:,:) = tmp(:,:,:,1)
        recvbuf3(:,:,:) = 0.0D0
        CALL mpi_sendrecv(sendbuf3(-3,0,0), (n1+2)*(n2+2)*7, MPI_REAL8, mg(lev)%rank_l, 6, &
                          recvbuf3(-3,0,0), (n1+2)*(n2+2)*7, MPI_REAL8, mg(lev)%rank_u, 6, &
                          comm, MPI_STATUS_IGNORE, ierr)

        tmp(:,:,:,n3+1) = recvbuf3(:,:,:)
     ELSE
        IF (period3) THEN
           tmp(:,:,:,0)    = tmp(:,:,:,n3)
           tmp(:,:,:,n3+1) = tmp(:,:,:,1)
        END IF
     END IF
#else
     IF (period3) THEN
        tmp(:,:,:,0)    = tmp(:,:,:,n3)
        tmp(:,:,:,n3+1) = tmp(:,:,:,1)
     END IF
#endif
  END IF

#else
  IF (period1) THEN
     tmp(:,0,:,:)    = tmp(:,n1,:,:)
     tmp(:,n1+1,:,:) = tmp(:,1,:,:)
  END IF

  IF (period2) THEN
     tmp(:,:,0,:)    = tmp(:,:,n2,:)
     tmp(:,:,n2+1,:) = tmp(:,:,1,:)
  END IF

  IF (period3) THEN
     tmp(:,:,:,0)    = tmp(:,:,:,n3)
     tmp(:,:,:,n3+1) = tmp(:,:,:,1)
  END IF

#endif

  DO k=1, n3
     DO j=1, n2
        DO i=1, n1

           IF ( (.NOT. period1 .AND. n1==1) .OR. &
                (.NOT. period2 .AND. n2==1)) THEN
              m(:,i,j,k) = 0.0D0
              m(0,i,j,k) = sum(l(-3:3,i,j,k)**2)
              IF (m(0,i,j,k) > 0.0D0) THEN
                 m(0,i,j,k) = l(0,i,j,k)/m(0,i,j,k)
              END IF
              CYCLE
           END IF

           a(:,:) = 0.0D0

           a(-12,-3) = tmp(-3,i,j,k-1)
           a(-11,-3) = tmp(-2,i,j,k-1)
           a(-10,-3) = tmp(-1,i,j,k-1)
           a(-9, -3) = tmp( 0,i,j,k-1)
           a(-8, -3) = tmp( 1,i,j,k-1)
           a(-7 ,-3) = tmp( 2,i,j,k-1)
           a( 0, -3) = tmp( 3,i,j,k-1)

           a(-11,-2) = tmp(-3,i,j-1,k)
           a(-6, -2) = tmp(-2,i,j-1,k)
           a(-5, -2) = tmp(-1,i,j-1,k)
           a(-4, -2) = tmp( 0,i,j-1,k)
           a(-3, -2) = tmp( 1,i,j-1,k)
           a( 0, -2) = tmp( 2,i,j-1,k)
           a( 7, -2) = tmp( 3,i,j-1,k)

           a(-10,-1) = tmp(-3,i-1,j,k)
           a(-5, -1) = tmp(-2,i-1,j,k)
           a(-2, -1) = tmp(-1,i-1,j,k)
           a(-1, -1) = tmp( 0,i-1,j,k)
           a( 0, -1) = tmp( 1,i-1,j,k)
           a( 3, -1) = tmp( 2,i-1,j,k)
           a( 8, -1) = tmp( 3,i-1,j,k)

           a(-9,  0) = tmp(-3,i,j,k)
           a(-4,  0) = tmp(-2,i,j,k)
           a(-1,  0) = tmp(-1,i,j,k)
           a( 0,  0) = tmp( 0,i,j,k)
           a( 1,  0) = tmp( 1,i,j,k)
           a( 4,  0) = tmp( 2,i,j,k)
           a( 9,  0) = tmp( 3,i,j,k)

           a(-8,  1) = tmp(-3,i+1,j,k)
           a(-3,  1) = tmp(-2,i+1,j,k)
           a( 0,  1) = tmp(-1,i+1,j,k)
           a( 1,  1) = tmp( 0,i+1,j,k)
           a( 2,  1) = tmp( 1,i+1,j,k)
           a( 5,  1) = tmp( 2,i+1,j,k)
           a(10,  1) = tmp( 3,i+1,j,k)

           a(-7,  2) = tmp(-3,i,j+1,k)
           a( 0,  2) = tmp(-2,i,j+1,k)
           a( 3,  2) = tmp(-1,i,j+1,k)
           a( 4,  2) = tmp( 0,i,j+1,k)
           a( 5,  2) = tmp( 1,i,j+1,k)
           a( 6,  2) = tmp( 2,i,j+1,k)
           a(11,  2) = tmp( 3,i,j+1,k)

           a(0,   3) = tmp(-3,i,j,k+1)
           a(7,   3) = tmp(-2,i,j,k+1)
           a(8,   3) = tmp(-1,i,j,k+1)
           a(9,   3) = tmp( 0,i,j,k+1)
           a(10,  3) = tmp( 1,i,j,k+1)
           a(11,  3) = tmp( 2,i,j,k+1)
           a(12,  3) = tmp( 3,i,j,k+1)

#ifdef CHOLESKY
           DO q=-3, 3
              DO p=-3, q
                 aa(p,q) = sum(a(:,p)*a(:,q))
              END DO
              IF (aa(q,q) == 0.0D0) aa(q,q) = 1.0D0
           END DO

           m(-3,i,j,k) = tmp( 3,i,  j,  k-1)
           m(-2,i,j,k) = tmp( 2,i,  j-1,k)
           m(-1,i,j,k) = tmp( 1,i-1,j,  k)
           m( 0,i,j,k) = tmp( 0,i,  j,  k)
           m( 1,i,j,k) = tmp(-1,i+1,j,  k)
           m( 2,i,j,k) = tmp(-2,i,  j+1,k)
           m( 3,i,j,k) = tmp(-3,i,  j,  k+1)

           CALL cholesky(7, aa(:,:), m(:,i,j,k))
#else
!           aa(:,:) =  matmul(transpose(a(:,:)),a(:,:))
!           DO n=-3, 3
!              IF (aa(n,n) == 0.0D0) aa(n,n) = 1.0D0
!           END DO
           DO q=-3, 3
              DO p=-3, 3
                 aa(p,q) = sum(a(:,p)*a(:,q))
              END DO
              IF (aa(q,q) == 0.0D0) aa(q,q) = 1.0D0
           END DO

!           CALL invmat(7, aa(:,:))
           DO invmat_j=-3, 3
              invmat_p = aa(invmat_j,invmat_j)
              aa(invmat_j,invmat_j) = 1.0D0
              aa(invmat_j,:)        = aa(invmat_j,:)/invmat_p
              DO invmat_i=-3, 3
                 IF (invmat_i == invmat_j) CYCLE
                 invmat_q = aa(invmat_i,invmat_j)
                 aa(invmat_i,invmat_j) = 0.0D0
                 aa(invmat_i,:)        = aa(invmat_i,:) - invmat_q*aa(invmat_j,:)
              END DO
           END DO
!  invmat is not expanded in iniline form on SR16000, I don"t know why.

           DO q=-3, 3
              m(q,i,j,k) = tmp(0,i,j,k)*aa(q,0)                                  &
                         + tmp(-1,i+1,j,k)*aa(q,1) + tmp(1,i-1,j,k)*aa(q,-1) &
                         + tmp(-2,i,j+1,k)*aa(q,2) + tmp(2,i,j-1,k)*aa(q,-2) &
                         + tmp(-3,i,j,k+1)*aa(q,3) + tmp(3,i,j,k-1)*aa(q,-3)
           END DO
#endif
        END DO
     END DO
  END DO
END SUBROUTINE sai3d

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE init_il
  INTEGER :: col(-2:2)
  INTEGER :: n1, n2
  INTEGER :: i, j, l, m

  CALL assert(mg(nlev)%n3==1, "error in solver3d%init_il")

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

        IF (mg(nlev)%l(0,i,j,1) /= 0.0D0) THEN
           il(col(0),col(-2)) = mg(nlev)%l(-2,i,j,1)
           il(col(0),col(-1)) = mg(nlev)%l(-1,i,j,1)
           il(col(0),col( 0)) = mg(nlev)%l( 0,i,j,1)
           il(col(0),col( 1)) = mg(nlev)%l( 1,i,j,1)
           il(col(0),col( 2)) = mg(nlev)%l( 2,i,j,1)
        ELSE
           il(col(0),col( 0)) = 1.0_MG_KIND
        END IF
     END DO
  END DO
  IF (singular) il(1,:) = 1.0_MG_KIND

  CALL invmat(n1*n2, il(:,:))

END SUBROUTINE init_il

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE apply_il(n1, n2, il, in, out)
  INTEGER,       INTENT(IN)  :: n1, n2
  REAL(MG_KIND), INTENT(IN)  :: il(n1*n2,n1*n2)
  REAL(MG_KIND), INTENT(IN)  :: in( 0:n1+1, 0:n2+1)
  REAL(MG_KIND), INTENT(OUT) :: out(0:n1+1, 0:n2+1)

  REAL(MG_KIND) :: tmp(n1*n2)

  INTEGER :: i, j

  tmp(:) = reshape(in(1:n1, 1:n2), (/n1*n2/))
  IF (singular) tmp(1) = sum(out(1:n1,1:n2))

  out(1:n1,1:n2) = reshape(matmul(il, tmp), (/n1, n2/))

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

SUBROUTINE finalize_solver3d
  INTEGER :: lev
  REAL(8) :: time
  REAL(8) :: tmp

  INTEGER :: i

  IF (.NOT. initialized) RETURN

  DEALLOCATE(work)

  DEALLOCATE(l3d)

  DO lev=0, nlev
#ifdef NO_ALLOCATABLE_IN_TYPE
     IF (ASSOCIATED(mg(lev)%r)) DEALLOCATE(mg(lev)%r)
     IF (ASSOCIATED(mg(lev)%e)) DEALLOCATE(mg(lev)%e)
     IF (ASSOCIATED(mg(lev)%t)) DEALLOCATE(mg(lev)%t)

     IF (ASSOCIATED(mg(lev)%l)) DEALLOCATE(mg(lev)%l)
     IF (ASSOCIATED(mg(lev)%m)) DEALLOCATE(mg(lev)%m)

     IF (ASSOCIATED(mg(lev)%p1)) DEALLOCATE(mg(lev)%p1)
     IF (ASSOCIATED(mg(lev)%p2)) DEALLOCATE(mg(lev)%p2)
     IF (ASSOCIATED(mg(lev)%p3)) DEALLOCATE(mg(lev)%p3)

     IF (ASSOCIATED(mg(lev)%q1)) DEALLOCATE(mg(lev)%q1)
     IF (ASSOCIATED(mg(lev)%q2)) DEALLOCATE(mg(lev)%q2)
     IF (ASSOCIATED(mg(lev)%q3)) DEALLOCATE(mg(lev)%q3)

     IF (ASSOCIATED(mg(lev)%lindx)) DEALLOCATE(mg(lev)%lindx)
     IF (ASSOCIATED(mg(lev)%mindx)) DEALLOCATE(mg(lev)%mindx)

#else
     IF (ALLOCATED(mg(lev)%r)) DEALLOCATE(mg(lev)%r)
     IF (ALLOCATED(mg(lev)%e)) DEALLOCATE(mg(lev)%e)
     IF (ALLOCATED(mg(lev)%t)) DEALLOCATE(mg(lev)%t)

     IF (ALLOCATED(mg(lev)%l)) DEALLOCATE(mg(lev)%l)
     IF (ALLOCATED(mg(lev)%m)) DEALLOCATE(mg(lev)%m)

     IF (ALLOCATED(mg(lev)%p1)) DEALLOCATE(mg(lev)%p1)
     IF (ALLOCATED(mg(lev)%p2)) DEALLOCATE(mg(lev)%p2)
     IF (ALLOCATED(mg(lev)%p3)) DEALLOCATE(mg(lev)%p3)

     IF (ALLOCATED(mg(lev)%q1)) DEALLOCATE(mg(lev)%q1)
     IF (ALLOCATED(mg(lev)%q2)) DEALLOCATE(mg(lev)%q2)
     IF (ALLOCATED(mg(lev)%q3)) DEALLOCATE(mg(lev)%q3)

     IF (ALLOCATED(mg(lev)%lindx)) DEALLOCATE(mg(lev)%lindx)
     IF (ALLOCATED(mg(lev)%mindx)) DEALLOCATE(mg(lev)%mindx)
#endif

  END DO

  IF (ALLOCATED(il)) DEALLOCATE(il)

#ifdef PARALLEL_MPI
  IF (ALLOCATED(ranks))  DEALLOCATE(ranks)
  IF (ALLOCATED(hranks)) DEALLOCATE(hranks)

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
  IF (ALLOCATED(sendbuf_u_r4)) DEALLOCATE(sendbuf_u_r4)
  IF (ALLOCATED(sendbuf_u_r8)) DEALLOCATE(sendbuf_u_r8)
  IF (ALLOCATED(sendbuf_l_r4)) DEALLOCATE(sendbuf_l_r4)
  IF (ALLOCATED(sendbuf_l_r8)) DEALLOCATE(sendbuf_l_r8)

  IF (ALLOCATED(recvbuf_e_r4)) DEALLOCATE(recvbuf_e_r4)
  IF (ALLOCATED(recvbuf_e_r8)) DEALLOCATE(recvbuf_e_r8)
  IF (ALLOCATED(recvbuf_w_r4)) DEALLOCATE(recvbuf_w_r4)
  IF (ALLOCATED(recvbuf_w_r8)) DEALLOCATE(recvbuf_w_r8)
  IF (ALLOCATED(recvbuf_n_r4)) DEALLOCATE(recvbuf_n_r4)
  IF (ALLOCATED(recvbuf_n_r8)) DEALLOCATE(recvbuf_n_r8)
  IF (ALLOCATED(recvbuf_s_r4)) DEALLOCATE(recvbuf_s_r4)
  IF (ALLOCATED(recvbuf_s_r8)) DEALLOCATE(recvbuf_s_r8)
  IF (ALLOCATED(recvbuf_u_r4)) DEALLOCATE(recvbuf_u_r4)
  IF (ALLOCATED(recvbuf_u_r8)) DEALLOCATE(recvbuf_u_r8)
  IF (ALLOCATED(recvbuf_l_r4)) DEALLOCATE(recvbuf_l_r4)
  IF (ALLOCATED(recvbuf_l_r8)) DEALLOCATE(recvbuf_l_r8)
#endif

  IF (rank==0) THEN
     IF (n_call >= 1) THEN
        WRITE(REPORT_UNIT, '(A,EN12.3)') "finalize solver3d: average iteration =", DBLE(n_iter)/DBLE(n_call)
#ifdef PROFILE
        WRITE(REPORT_UNIT, '("profile of solver3d")')
        WRITE(REPORT_UNIT, '("  normtime: ", EN12.3, " sec.")') profile_time('norm')
        WRITE(REPORT_UNIT, '("  dottime:  ", EN12.3, ", ", EN12.3, " sec.")') profile_time('dot'), profile_time('dot2')
        WRITE(REPORT_UNIT, '("  axpytime:  ", EN12.3, ", ", EN12.3, " sec.")') profile_time('axpy'), profile_time('xpby')

        WRITE(REPORT_UNIT, '("  finetime: ", EN12.3, " sec.")') profile_time('fine')
        DO lev=0, nlev
           WRITE(REPORT_UNIT, '("    level=", I2, ": ", EN12.3, " sec.")') lev, profile_time('fine', lev)
        END DO

        WRITE(REPORT_UNIT, '("  coarsetime: ", EN12.3, " sec.")') profile_time('coarse')
        DO lev=0, nlev
           WRITE(REPORT_UNIT, '("    level=", I2, ": ", EN12.3, " sec.")') lev, profile_time('coarse', lev)
        END DO

        WRITE(REPORT_UNIT, '("  applytime: ", EN12.3, " sec.")') profile_time('apply1') + profile_time('apply2') + profile_time('apply3')
        DO lev=0, nlev
           WRITE(REPORT_UNIT, '("    level=", I2, ": ", EN12.3, " ", EN12.3, " ", EN12.3, " sec.")') &
                lev, profile_time('apply1', lev), profile_time('apply2', lev), profile_time('apply3', lev)
        END DO

        WRITE(REPORT_UNIT, '("  synctime: ", EN12.3, " sec.")') profile_time('sync')
        DO lev=0, nlev
           WRITE(REPORT_UNIT, '("    level=", I2, ": ", EN12.3, " sec.")') lev, profile_time('sync', lev)
        END DO

        WRITE(REPORT_UNIT, '("  gathtime: ", EN12.3, " sec.")') profile_time('gath')
        WRITE(REPORT_UNIT, '("  scattime: ", EN12.3, " sec.")') profile_time('scat')
#endif
     END IF
  END IF

  initialized = .FALSE.
  eps = 1.0D-8
  nlev = max_level

END SUBROUTINE finalize_solver3d

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE solve3d(q, p)
  REAL(8), INTENT(IN)    :: q(0:n1+1, 0:n2+1, 0:n3+1)
  REAL(8), INTENT(INOUT) :: p(0:n1+1, 0:n2+1, 0:n3+1)

#ifdef DEBUG__
  REAL(8) :: sumq, tmp
  REAL(8) :: sump
  INTEGER :: ierr
#endif

#ifdef PARALLEL_MPI
  IF (comm == MPI_COMM_NULL) RETURN
#endif

#ifdef DEBUG__
  CALL assert(initialized, "solver3d is not initialized.")
  CALL assert(n3>1, "solver dimension  missmatch")

  sumq = sum(q(1:n1, 1:n2, 1:n3))
#ifdef PARALLEL_MPI
  tmp = sumq
  CALL mpi_reduce(tmp, sumq, 1, MPI_REAL8, MPI_SUM, 0, comm, ierr)
#endif
  IF (rank==0) WRITE(REPORT_UNIT, '("SUM(Q)=",EN12.3)') sumq

!  CALL mpi_reduce(sum(p(1:n1, 1:n2, 1:n3)), sump, 1, MPI_REAL8, MPI_SUM, 0, comm, ierr)
!  IF (rank==0) WRITE(REPORT_UNIT,'("SUM(P)=",EN12.3)') sump
#endif

!2013/6/16: for iceshelf (MGCG is unstable without zero-fill, I don"t know why)
!  p(:,:,:)   = 0.0D0

  SELECT CASE(solver_method)
     CASE(0)
        CALL solver_cg3d(q, p)
     CASE(1)
        CALL solver_mg3d(q, p)
     CASE(2)
        CALL solver_mgcg3d(q, p)
     CASE(3)
        CALL solver_mgcgs3d(q, p)
     CASE(-1)
        CALL solver_mg3d_omp(q, p)
     CASE(-2)
        CALL solver_mgcg3d_omp(q, p)
  END SELECT

END SUBROUTINE solve3d

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE solver_mgcg3d(q, p)
  REAL(8), INTENT(IN)    :: q(0:n1+1, 0:n2+1, 0:n3+1)
  REAL(8), INTENT(INOUT) :: p(0:n1+1, 0:n2+1, 0:n3+1)

  REAL(8) :: res, res_c

  REAL(8) :: alpha, beta, rho

  REAL(8) :: u(0:n1+1, 0:n2+1, 0:n3+1)
  REAL(8) :: v(0:n1+1, 0:n2+1, 0:n3+1)
  REAL(8) :: r(0:n1+1, 0:n2+1, 0:n3+1)
  REAL(8) :: e(0:n1+1, 0:n2+1, 0:n3+1)

  INTEGER :: n

!$OMP PARALLEL WORKSHARE
  v(:,:,:) = 0.0D0
  r(:,:,:) = 0.0D0
  e(:,:,:) = 0.0D0
!$OMP END PARALLEL WORKSHARE

MGPROFILE_BEGIN('norm', 0)
  CALL norm3d(n1, n2, n3, q, res_c, comm)
MGPROFILE_END('norm', 0)
  res_c = res_c * eps*eps

MGPROFILE_BEGIN('apply2', 0)
  CALL apply3d2(n1, n2, n3, l3d, q, p, r, kindex=mg(0)%lindx)
MGPROFILE_END('apply2', 0)

  CALL sync3d_begin(0, r)

MGPROFILE_BEGIN('norm', 0)
  CALL norm3d(n1, n2, n3, r, res, comm)
MGPROFILE_END('norm', 0)

  CALL sync3d_end(0, r)

  IF (res <= res_c) THEN
     CALL checkout(0, res, res_c)
     RETURN
  END IF

  CALL multigrid3d(r, e)

  CALL sync3d(0, e)

MGPROFILE_BEGIN('copy', 0)
!$OMP PARALLEL WORKSHARE
  u(:,:,:) = e(:,:,:)
!$OMP END PARALLEL WORKSHARE
MGPROFILE_END('copy', 0)

  DO n=1, itmax
MGPROFILE_BEGIN('apply1', 0)
     CALL apply3d1(n1, n2, n3, l3d, u, v, kindex=mg(0)%lindx)
MGPROFILE_END('apply1', 0)

     CALL sync3d_begin(0, v)

MGPROFILE_BEGIN('dot2', 0)
     CALL dot3d2(n1, n2, n3, u, v, r, rho, alpha, comm)
MGPROFILE_END('dot2', 0)

     CALL sync3d_end(0, v)

     IF (fallback .AND. abs(rho) < fallback_eps) THEN
        IF (itreport .AND. rank==0) WRITE(REPORT_UNIT, '(A,I3)') 'MGCG unstable, fallback to stand-alone multigrid: n=', n

        CALL solver_mg3d(q, p)
        RETURN
     END IF

     alpha = alpha / rho

MGPROFILE_BEGIN('axpy', 0)
     CALL axpy3d(n1, n2, n3,  alpha, u, p)
     CALL axpy3d(n1, n2, n3, -alpha, v, r)
MGPROFILE_END('axpy', 0)

MGPROFILE_BEGIN('norm', 0)
     CALL norm3d(n1, n2, n3, r, res, comm)
MGPROFILE_END('norm', 0)

     IF (res <= res_c) THEN
        CALL checkout(n, res, res_c)
        RETURN
     END IF

     CALL multigrid3d(r, e)

     CALL sync3d_begin(0, e)

MGPROFILE_BEGIN('dot', 0)
     CALL dot3d(n1, n2, n3, v, e, beta, comm)
MGPROFILE_END('dot', 0)
     beta = -beta / rho

     CALL sync3d_end(0, e)

MGPROFILE_BEGIN('xpby', 0)
     CALL xpby3d(n1, n2, n3, e, beta, u)
MGPROFILE_END('xpby', 0)
  END DO

  IF (rank==0) WRITE(REPORT_UNIT,'(A,I3)') 'MGCG unstable, fallback to stand-alone multigrid: n=', n

!$OMP PARALLEL WORKSHARE
  p(:,:,:) = 0.0D0
!$OMP END PARALLEL WORKSHARE

  CALL solver_mg3d(q, p)
  RETURN

  CALL quit_message(res, res_c)

END SUBROUTINE solver_mgcg3d

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE solver_mgcg3d_omp(q, p)
  REAL(8), INTENT(IN)    :: q(0:n1+1, 0:n2+1, 0:n3+1)
  REAL(8), INTENT(INOUT) :: p(0:n1+1, 0:n2+1, 0:n3+1)

  REAL(8) :: res, res_c

  REAL(8) :: alpha, beta, rho

  REAL(8) :: u(0:n1+1, 0:n2+1, 0:n3+1)
  REAL(8) :: v(0:n1+1, 0:n2+1, 0:n3+1)
  REAL(8) :: r(0:n1+1, 0:n2+1, 0:n3+1)
  REAL(8) :: e(0:n1+1, 0:n2+1, 0:n3+1)

  INTEGER :: n
  INTEGER :: i, j, k

  INTEGER, SAVE :: tid = 0
  INTEGER, SAVE :: kstart, kend, start_offset, end_offset
!$OMP threadprivate(tid, kstart, kend, start_offset, end_offset)

  REAL(8) :: tmp(n3,2)

!$OMP PARALLEL
!$ tid = omp_get_thread_num()
  kstart = n3_start(tid,0)
  kend   = n3_start(tid,0) + n3_span(tid,0) - 1

  IF (kstart==1) THEN
     start_offset = 1
  ELSE
     start_offset = 0
  END IF

  IF (kend==n3) THEN
     end_offset = 1
  ELSE
     end_offset = 0
  END IF

  r(:,:,kstart-start_offset:kend+end_offset) = 0.0D0
  v(:,:,kstart-start_offset:kend+end_offset) = 0.0D0
  e(:,:,kstart-start_offset:kend+end_offset) = 0.0D0

MGPROFILE_BEGIN_THREAD('norm', 0, tid)
! CALL norm3d(n1, n2, n3, q, res_c, comm)
  DO k=kstart, kend
     tmp(k,1) = 0.0D0
     DO j=1, n2
     DO i=1, n1
        tmp(k,1) = tmp(k,1) + q(i,j,k)**2
     END DO
     END DO
  END DO
!$OMP BARRIER
!$OMP SINGLE
  CALL allreduce_sum(sum(tmp(1:n3,1)), res_c, comm)
  res_c = res_c * eps*eps
!$OMP END SINGLE
MGPROFILE_END_THREAD('norm', 0, tid)

MGPROFILE_BEGIN_THREAD('apply2', 0, tid)
  CALL apply3d2(n1, n2, n3, l3d, q, p, r, kstart, kend, kindex=mg(0)%lindx)
MGPROFILE_END_THREAD('apply2', 0, tid)

  CALL sync3d_begin(0, r, tid)

MGPROFILE_BEGIN_THREAD('norm', 0, tid)
! CALL norm3d(n1, n2, n3, r, res, comm)
  DO k=kstart, kend
     tmp(k,1) = 0.0D0
     DO j=1, n2
     DO i=1, n1
        tmp(k,1) = tmp(k,1) + r(i,j,k)**2
     END DO
     END DO
  END DO
!$OMP BARRIER
!$OMP SINGLE
  CALL allreduce_sum(sum(tmp(1:n3,1)), res, comm)
!$OMP END SINGLE
MGPROFILE_END_THREAD('norm', 0, tid)

  CALL sync3d_end(0, r, tid)
!$OMP END PARALLEL

  IF (res <= res_c) THEN
     CALL checkout(0, res, res_c)
     RETURN
  END IF

!$OMP PARALLEL
  CALL multigrid3d_omp(r, e, tid)

  CALL sync3d(0, e, tid)

MGPROFILE_BEGIN_THREAD('copy', 0, tid)
  u(:,:,kstart-start_offset:kend+end_offset) = e(:,:,kstart-start_offset:kend+end_offset)
MGPROFILE_END_THREAD('copy', 0, tid)
!$OMP END PARALLEL

  DO n=1, itmax
!$OMP PARALLEL
MGPROFILE_BEGIN_THREAD('apply1', 0, tid)
     CALL apply3d1(n1, n2, n3, l3d, u, v, kstart, kend, kindex=mg(0)%lindx)
MGPROFILE_END_THREAD('apply1', 0, tid)

     CALL sync3d_begin(0, v, tid)

MGPROFILE_BEGIN_THREAD('dot2', 0, tid)
!    CALL dot3d2(n1, n2, n3, u, v, r, rho, alpha, comm)
     DO k=kstart, kend
        tmp(k,1) = 0.0D0
        tmp(k,2) = 0.0D0
        DO j=1, n2
        DO i=1, n1
           tmp(k,1) = tmp(k,1) + u(i,j,k)*v(i,j,k)
           tmp(k,2) = tmp(k,2) + u(i,j,k)*r(i,j,k)
        END DO
     END DO
     END DO
!$OMP BARRIER
!$OMP SINGLE
     CALL allreduce_sum(sum(tmp(1:n3,1)), rho,   comm)
     CALL allreduce_sum(sum(tmp(1:n3,2)), alpha, comm)
!$OMP END SINGLE
MGPROFILE_END_THREAD('dot2', 0, tid)

     CALL sync3d_end(0, v, tid)
!$OMP END PARALLEL
     IF (fallback .AND. abs(rho) < fallback_eps) THEN
        IF (itreport .AND. rank==0) WRITE(REPORT_UNIT,'(A,I3)') 'MGCG unstable, fallback to stand-alone multigrid: n=', n

        CALL solver_mg3d_omp(q, p)
        RETURN
     END IF

     alpha = alpha / rho

!$OMP PARALLEL
MGPROFILE_BEGIN_THREAD('axpy', 0, tid)
!    CALL axpy3d(n1, n2, n3,  alpha, u, p)
!    CALL axpy3d(n1, n2, n3, -alpha, v, r)
     p(:,:,kstart-start_offset:kend+end_offset) = p(:,:,kstart-start_offset:kend+end_offset) + alpha*u(:,:,kstart-start_offset:kend+end_offset)
     r(:,:,kstart-start_offset:kend+end_offset) = r(:,:,kstart-start_offset:kend+end_offset) - alpha*v(:,:,kstart-start_offset:kend+end_offset)
MGPROFILE_END_THREAD('axpy', 0, tid)

MGPROFILE_BEGIN_THREAD('norm', 0, tid)
!    CALL norm3d(n1, n2, n3, r, res, comm)
     DO k=kstart, kend
        tmp(k,1) = 0.0D0
        DO j=1, n2
        DO i=1, n1
           tmp(k,1) = tmp(k,1) + r(i,j,k)**2
        END DO
     END DO
     END DO
!$OMP BARRIER
!$OMP SINGLE
     CALL allreduce_sum(sum(tmp(1:n3,1)), res, comm)
!$OMP END SINGLE
MGPROFILE_END_THREAD('norm', 0, tid)
!$OMP END PARALLEL

     IF (res <= res_c) THEN
        CALL checkout(n, res, res_c)
        RETURN
     END IF

!$OMP PARALLEL
     CALL multigrid3d_omp(r, e, tid)

     CALL sync3d_begin(0, e, tid)

MGPROFILE_BEGIN_THREAD('dot', 0, tid)
!    CALL dot3d(n1, n2, n3, v, e, beta, comm)
     DO k=kstart, kend
        tmp(k,1) = 0.0D0
        DO j=1, n2
        DO i=1, n1
           tmp(k,1) = tmp(k,1) + v(i,j,k)*e(i,j,k)
        END DO
     END DO
     END DO
!$OMP BARRIER
!$OMP SINGLE
     CALL allreduce_sum(sum(tmp(1:n3,1)), beta, comm)
!$OMP END SINGLE
MGPROFILE_END_THREAD('dot', 0, tid)

!$OMP SINGLE
     beta = -beta / rho
!$OMP END SINGLE

     CALL sync3d_end(0, e, tid)

MGPROFILE_BEGIN_THREAD('xpby', 0, tid)
!    CALL xpby3d(n1, n2, n3, e, beta, u)
     u(:,:,kstart-start_offset:kend+end_offset) = beta*u(:,:,kstart-start_offset:kend+end_offset) + e(:,:,kstart-start_offset:kend+end_offset)
MGPROFILE_END_THREAD('xpby', 0, tid)
!$OMP END PARALLEL
  END DO

  IF (rank==0) WRITE(REPORT_UNIT, '(A,I3)') 'MGCG unstable, fallback to stand-alone multigrid: n=', n
  p(:,:,:) = 0.0D0
  CALL solver_mg3d_omp(q, p)
  RETURN

  CALL quit_message(res, res_c)

END SUBROUTINE solver_mgcg3d_omp

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE solver_mgcgs3d(q, p)
  REAL(8), INTENT(IN)    :: q(0:n1+1, 0:n2+1, 0:n3+1)
  REAL(8), INTENT(INOUT) :: p(0:n1+1, 0:n2+1, 0:n3+1)

  REAL(8) :: res, res_c

  REAL(8) :: alpha, beta, rho

  REAL(8) :: r(0:n1+1, 0:n2+1, 0:n3+1)
  REAL(8) :: e(0:n1+1, 0:n2+1, 0:n3+1)

  REAL(8) :: s(0:n1+1, 0:n2+1, 0:n3+1)
  REAL(8) :: t(0:n1+1, 0:n2+1, 0:n3+1)
  REAL(8) :: u(0:n1+1, 0:n2+1, 0:n3+1)
  REAL(8) :: v(0:n1+1, 0:n2+1, 0:n3+1)

  REAL(8) :: r0(0:n1+1, 0:n2+1, 0:n3+1)

  INTEGER :: n

!$OMP PARALLEL WORKSHARE
  r0(:,:,:) = 0.0D0

  s(:,:,:)  = 0.0D0
  t(:,:,:)  = 0.0D0
  e(:,:,:)  = 0.0D0
!$OMP END PARALLEL WORKSHARE

  CALL norm3d(n1, n2, n3, q, res_c, comm)
  res_c = res_c *eps*eps

MGPROFILE_BEGIN('apply2', 0)
  CALL apply3d2(n1, n2, n3, l3d, q, p, r0, kindex=mg(0)%lindx)
MGPROFILE_END('apply2', 0)

  CALL sync3d_begin(0, r0)

  CALL norm3d(n1, n2, n3, r0, res, comm)

  CALL sync3d_end(0, r0)

  IF (res <= res_c) THEN
     CALL checkout(0, res, res_c)
     RETURN
  END IF

!$OMP PARALLEL WORKSHARE
  u(:,:,:) = r0(:,:,:)
  v(:,:,:) = r0(:,:,:)
  r(:,:,:) = r0(:,:,:)
!$OMP END PARALLEL WORKSHARE
  rho = res

  DO n=1, itmax
     IF (fallback .AND. abs(rho) < fallback_eps) THEN
        IF (itreport .AND. rank==0) WRITE(REPORT_UNIT, '(A,I3)') 'MGCGS unstable, fallback to stand-alone multigrid: n=', n

        CALL solver_mg3d(q, p)
        RETURN
     END IF

     CALL multigrid3d(u, e)

MGPROFILE_BEGIN('apply1', 0)
     CALL apply3d1(n1, n2, n3, l3d, e, s, kindex=mg(0)%lindx)
MGPROFILE_END('apply1', 0)

     CALL sync3d_begin(0, s)

     CALL dot3d(n1, n2, n3, r0, s, alpha, comm)
     alpha = rho / alpha

     CALL sync3d_end(0, s)


!$OMP PARALLEL WORKSHARE
     t(:,:,:) = v(:,:,:) - alpha * s(:,:,:)
!$OMP END PARALLEL WORKSHARE

     CALL multigrid3d(v + t, e)

     CALL axpy3d(n1, n2, n3, alpha, e, p)

MGPROFILE_BEGIN('apply1', 0)
     CALL apply3d1(n1, n2, n3, l3d, e, s, kindex=mg(0)%lindx)
MGPROFILE_END('apply1', 0)

     CALL sync3d(0, s)

     CALL axpy3d(n1, n2, n3, -alpha, s, r)

     CALL norm3d(n1, n2, n3, r, res, comm)

     IF (res <= res_c) THEN
        CALL checkout(n, res, res_c)
        RETURN
     END IF

     beta = rho
     CALL dot3d(n1, n2, n3, r0, r, rho, comm)
     beta = rho/beta

!$OMP PARALLEL WORKSHARE
     v(:,:,:) = r(:,:,:) + beta * t(:,:,:)
     u(:,:,:) = v(:,:,:) + beta * (t(:,:,:) + beta * u(:,:,:))
!$OMP END PARALLEL WORKSHARE
  END DO

  CALL quit_message(res, res_c)

END SUBROUTINE solver_mgcgs3d

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE solver_cg3d(q, p)
  REAL(8), INTENT(IN)    :: q(0:n1+1, 0:n2+1, 0:n3+1)
  REAL(8), INTENT(INOUT) :: p(0:n1+1, 0:n2+1, 0:n3+1)

  REAL(8) :: res, res_c

  REAL(8) :: alpha, beta, rho

  REAL(8) :: u(0:n1+1, 0:n2+1, 0:n3+1)
  REAL(8) :: v(0:n1+1, 0:n2+1, 0:n3+1)
  REAL(8) :: r(0:n1+1, 0:n2+1, 0:n3+1)

  INTEGER :: n
  INTEGER :: i, j, k

!$OMP PARALLEL WORKSHARE
  v(:,:,:) = 0.0D0
  r(:,:,:) = 0.0D0
!$OMP END PARALLEL WORKSHARE

  CALL norm3d(n1, n2, n3, q, res_c, comm)
  res_c = res_c * eps*eps

MGPROFILE_BEGIN('apply2', 0)
  CALL apply3d2(n1, n2, n3, l3d, q, p, r, kindex=mg(0)%lindx)
MGPROFILE_END('apply2', 0)

  CALL sync3d_begin(0, r)

  CALL norm3d(n1, n2, n3, r, res, comm)

  CALL sync3d_end(0, r)

  IF (res <= res_c) THEN
     CALL checkout(0, res, res_c)
     RETURN
  END IF

  CALL sync3d(0, r)

!$OMP PARALLEL WORKSHARE
  u(:,:,:) = r(:,:,:)
!$OMP END PARALLEL WORKSHARE

  DO n=1, itmax
MGPROFILE_BEGIN('apply1', 0)
     CALL apply3d1(n1, n2, n3, l3d, u, v, kindex=mg(0)%lindx)
MGPROFILE_END('apply1', 0)

     CALL sync3d_begin(0, v)

     CALL dot3d2(n1, n2, n3, u, v, r, rho, alpha, comm)

     CALL sync3d_end(0, v)

     alpha = alpha / rho

     CALL axpy3d(n1, n2, n3,  alpha, u, p)
     CALL axpy3d(n1, n2, n3, -alpha, v, r)

     CALL norm3d(n1, n2, n3, r, res, comm)

     IF (res <= res_c) THEN
        CALL checkout(n, res, res_c)
        RETURN
     END IF

     CALL sync3d_begin(0, r)

     CALL dot3d(n1, n2, n3, v, r, beta, comm)

     beta = -beta / rho

     CALL sync3d_end(0, r)

     CALL xpby3d(n1, n2, n3, r, beta, u)
  END DO

  CALL quit_message(res, res_c)

END SUBROUTINE solver_cg3d

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE solver_mg3d(q, p)
  REAL(8), INTENT(IN)    :: q(0:n1+1, 0:n2+1, 0:n3+1)
  REAL(8), INTENT(INOUT) :: p(0:n1+1, 0:n2+1, 0:n3+1)

  REAL(8) :: res, res_c

  REAL(8) :: r(0:n1+1, 0:n2+1, 0:n3+1)
  REAL(8) :: e(0:n1+1, 0:n2+1, 0:n3+1)

  INTEGER :: n

!$OMP PARALLEL WORKSHARE
  r(:,:,:) = 0.0D0
  e(:,:,:) = 0.0D0
!$OMP END PARALLEL WORKSHARE

  CALL norm3d(n1, n2, n3, q, res_c, comm)
  res_c = res_c *eps*eps

  DO n=1, itmax
MGPROFILE_BEGIN('apply2', 0)
     CALL apply3d2(n1, n2, n3, l3d, q, p, r, kindex=mg(0)%lindx)
MGPROFILE_END('apply2', 0)

     CALL sync3d_begin(0, r)

     CALL norm3d(n1, n2, n3, r, res, comm)

     CALL sync3d_end(0, r)

     IF (res <= res_c) THEN
        CALL checkout(n, res, res_c)
        RETURN
     END IF

     CALL multigrid3d(r, e)

     CALL sync3d(0, e)

!$OMP PARALLEL WORKSHARE
     p(:,:,:) = p(:,:,:) + e(:,:,:)
!$OMP END PARALLEL WORKSHARE
  END DO

  CALL quit_message(res, res_c)

END SUBROUTINE solver_mg3d

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE solver_mg3d_omp(q, p)
  REAL(8), INTENT(IN)    :: q(0:n1+1, 0:n2+1, 0:n3+1)
  REAL(8), INTENT(INOUT) :: p(0:n1+1, 0:n2+1, 0:n3+1)

  REAL(8) :: res, res_c

  REAL(8) :: r(0:n1+1, 0:n2+1, 0:n3+1)
  REAL(8) :: e(0:n1+1, 0:n2+1, 0:n3+1)

  INTEGER :: n
  INTEGER :: i, j, k

  INTEGER, SAVE :: tid = 0
  INTEGER, SAVE :: kstart, kend, start_offset, end_offset
!$OMP threadprivate(tid, kstart, kend, start_offset, end_offset)

  REAL(8) :: tmp(n3)

!$OMP PARALLEL
!$ tid = omp_get_thread_num()
  kstart = n3_start(tid,0)
  kend   = n3_start(tid,0) + n3_span(tid,0) - 1

  IF (kstart==1) THEN
     start_offset = 1
  ELSE
     start_offset = 0
  END IF

  IF (kend==n3) THEN
     end_offset = 1
  ELSE
     end_offset = 0
  END IF

  r(:,:,kstart-start_offset:kend+end_offset) = 0.0D0
  e(:,:,kstart-start_offset:kend+end_offset) = 0.0D0

! CALL norm3d(n1, n2, n3, q, res_c, comm)
  DO k=kstart, kend
     tmp(k) = 0.0D0
     DO j=1, n2
     DO i=1, n1
        tmp(k) = tmp(k) + q(i,j,k)**2
     END DO
     END DO
  END DO
!$OMP BARRIER
!$OMP SINGLE
  CALL allreduce_sum(sum(tmp(1:n3)), res_c, comm)
  res_c = res_c *eps*eps
!$OMP END SINGLE
!$OMP END PARALLEL

  DO n=1, itmax
!$OMP PARALLEL
MGPROFILE_BEGIN_THREAD('apply2', 0, tid)
     CALL apply3d2(n1, n2, n3, l3d, q, p, r, kstart, kend, kindex=mg(0)%lindx)
MGPROFILE_END_THREAD('apply2', 0, tid)

     CALL sync3d_begin(0, r, tid)

!    CALL norm3d(n1, n2, n3, r, res, comm)
     DO k=kstart, kend
        tmp(k) = 0.0D0
        DO j=1, n2
        DO i=1, n1
           tmp(k) = tmp(k) + r(i,j,k)**2
        END DO
        END DO
     END DO
!$OMP BARRIER
!$OMP SINGLE
     CALL allreduce_sum(sum(tmp(1:n3)), res, comm)
!$OMP END SINGLE

     CALL sync3d_end(0, r, tid)

!$OMP END PARALLEL
     IF (res <= res_c) THEN
        CALL checkout(n, res, res_c)
        RETURN
     END IF

!$OMP PARALLEL
     CALL multigrid3d_omp(r, e, tid)

     CALL sync3d(0, e, tid)

MGPROFILE_BEGIN_THREAD('axpy', 0, tid)
!    CALL axpy3d(n1, n2, n3, 1.0D0, e, p, kstart-start_offset, kend+end_offset)
     p(:,:,kstart-start_offset:kend+end_offset) = p(:,:,kstart-start_offset:kend+end_offset) + e(:,:,kstart-start_offset:kend+end_offset)
MGPROFILE_END_THREAD('axpy', 0, tid)
!$OMP END PARALLEL
  END DO

  CALL quit_message(res, res_c)

END SUBROUTINE solver_mg3d_omp

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE checkout(n, res, res_c)
  INTEGER, INTENT(IN) :: n
  REAL(8), INTENT(IN) :: res
  REAL(8), INTENT(IN) :: res_c

  IF (itreport .AND. res > 0 .AND. rank==0)  THEN
     WRITE(REPORT_UNIT, '(A,I5,A,ES10.3)') 'solver3d: iteration count=', n, &
                                           ', relative residual=', sqrt(res/res_c)*eps
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

  IF (itreport .AND. rank==0)  THEN
     WRITE(REPORT_UNIT,'(A,I5,A,ES10.3)') 'solver3d: iteration count exceeds ', itmax, &
                                          ', relative residual=', sqrt(res/res_c)*eps
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

  IF (lev==nlev) THEN
     IF (use_il) THEN
        CALL apply_il(mg(lev)%n1, mg(lev)%n2, il(:,:), mg(lev)%r(:,:,1), mg(lev)%e(:,:,1))
     ELSE
        CALL apply2d1(mg(lev)%n1, mg(lev)%n2, mg(lev)%m(:,:,:,1), mg(lev)%r(:,:,1), mg(lev)%e(:,:,1))
        CALL sync2d(lev, mg(lev)%e(:,:,1))
     END IF
     RETURN
  END IF

#ifdef PARALLEL_MPI
  IF (lev == glev) THEN
!     CALL apply2d1(mg(lev)%n1, mg(lev)%n2, mg(lev)%m(:,:,:,1), mg(lev)%r(:,:,1), mg(lev)%e(:,:,1))
!     CALL sync2d(lev, mg(lev)%e(:,:,1))
     CALL gather_h(mg(lev)%r(:,:,1), mg(lev+1)%r(:,:,1))
     IF (hrank == 0) CALL mgsmooth2d(lev+1, cyc)
     CALL mpi_barrier(hcomm, ierr)
     CALL scatter_h(mg(lev+1)%e(:,:,1), mg(lev)%e(:,:,1))
     RETURN
  END IF
#endif

  IF (cmode == 1  .OR. (cmode==2 .AND. cyc==0)) THEN
     n = 1
  ELSE
     n = 0
  END IF

  CALL apply2d1(mg(lev)%n1, mg(lev)%n2, mg(lev)%m(:,:,:,1), mg(lev)%r(:,:,1), mg(lev)%e(:,:,1))

  DO i=0, n
     CALL sync2d(lev, mg(lev)%e(:,:,1))

     CALL apply2d2(mg(lev)%n1, mg(lev)%n2, mg(lev)%l(:,:,:,1), mg(lev)%r(:,:,1), mg(lev)%e(:,:,1), mg(lev)%t(:,:,1))

     CALL coarse2d(mg(lev)%n1, mg(lev)%n2, mg(lev+1)%n1, mg(lev+1)%n2, mg(lev)%p1, mg(lev)%p2, mg(lev)%q1, mg(lev)%q2, mg(lev)%t(:,:,1), mg(lev+1)%r(:,:,1))

     CALL sync2d(lev+1, mg(lev+1)%r(:,:,1))

     CALL mgsmooth2d(lev+1, cyc+i)

     CALL fine2d(mg(lev)%n1, mg(lev)%n2, mg(lev+1)%n1, mg(lev+1)%n2, mg(lev)%p1, mg(lev)%p2, mg(lev)%q1, mg(lev)%q2, mg(lev+1)%e(:,:,1), mg(lev)%e(:,:,1))

     CALL sync2d(lev, mg(lev)%e(:,:,1))

     CALL apply2d2(mg(lev)%n1, mg(lev)%n2, mg(lev)%l(:,:,:,1), mg(lev)%r(:,:,1), mg(lev)%e(:,:,1), mg(lev)%t(:,:,1))

     CALL sync2d(lev, mg(lev)%t(:,:,1))

     CALL apply2d3(mg(lev)%n1, mg(lev)%n2, mg(lev)%m(:,:,:,1), mg(lev)%t(:,:,1), mg(lev)%e(:,:,1))
  END DO

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

  IF (lev==nlev) THEN
     IF (use_il) THEN
!$OMP SINGLE
        CALL apply_il(mg(lev)%n1, mg(lev)%n2, il(:,:), mg(lev)%r(:,:,1), mg(lev)%e(:,:,1))
!$OMP END SINGLE
     ELSE
!$OMP SINGLE
        CALL apply2d1(mg(lev)%n1, mg(lev)%n2, mg(lev)%m(:,:,:,1), mg(lev)%r(:,:,1), mg(lev)%e(:,:,1), 1, mg(lev)%n2)
!$OMP END SINGLE
        CALL sync2d(lev, mg(lev)%e(:,:,1), tid)
     END IF
     RETURN
  END IF

#ifdef PARALLEL_MPI
  IF (lev == glev) THEN
!     CALL apply2d1(mg(lev)%n1, mg(lev)%n2, mg(lev)%m(:,:,:,1), mg(lev)%r(:,:,1), mg(lev)%e(:,:,1))
!     CALL sync2d(lev, mg(lev)%e(:,:,1))
!$OMP SINGLE
     CALL gather_h(mg(lev)%r(:,:,1), mg(lev+1)%r(:,:,1))
!$OMP END SINGLE
     IF (hrank == 0) CALL mgsmooth2d_omp(lev+1, cyc, tid)
!$OMP SINGLE
     CALL mpi_barrier(hcomm, ierr)
     CALL scatter_h(mg(lev+1)%e(:,:,1), mg(lev)%e(:,:,1))
!$OMP END SINGLE
     RETURN
  END IF
#endif

  IF (cmode == 1  .OR. (cmode==2 .AND. cyc==0)) THEN
     n = 1
  ELSE
     n = 0
  END IF

  CALL apply2d1(mg(lev)%n1, mg(lev)%n2, mg(lev)%m(:,:,:,1), mg(lev)%r(:,:,1), mg(lev)%e(:,:,1), jstart, jend)

  DO i=0, n
     CALL sync2d(lev, mg(lev)%e(:,:,1), tid)

     CALL apply2d2(mg(lev)%n1, mg(lev)%n2, mg(lev)%l(:,:,:,1), mg(lev)%r(:,:,1), mg(lev)%e(:,:,1), mg(lev)%t(:,:,1), jstart, jend)

!$OMP BARRIER
     CALL coarse2d(mg(lev)%n1, mg(lev)%n2, mg(lev+1)%n1, mg(lev+1)%n2, mg(lev)%p1, mg(lev)%p2, &
                   mg(lev)%q1, mg(lev)%q2, mg(lev)%t(:,:,1), mg(lev+1)%r(:,:,1), &
                   n2_start(tid,lev+1), n2_start(tid,lev+1)+n2_span(tid,lev+1)-1)

     CALL sync2d(lev+1, mg(lev+1)%r(:,:,1), tid)

     CALL mgsmooth2d_omp(lev+1, cyc+i, tid)

     CALL fine2d(mg(lev)%n1, mg(lev)%n2, mg(lev+1)%n1, mg(lev+1)%n2, mg(lev)%p1, mg(lev)%p2, &
                 mg(lev)%q1, mg(lev)%q2, mg(lev+1)%e(:,:,1), mg(lev)%e(:,:,1), jstart, jend)

     CALL sync2d(lev, mg(lev)%e(:,:,1), tid)

     CALL apply2d2(mg(lev)%n1, mg(lev)%n2, mg(lev)%l(:,:,:,1), mg(lev)%r(:,:,1), mg(lev)%e(:,:,1), mg(lev)%t(:,:,1), jstart, jend)

     CALL sync2d(lev, mg(lev)%t(:,:,1), tid)

     CALL apply2d3(mg(lev)%n1, mg(lev)%n2, mg(lev)%m(:,:,:,1), mg(lev)%t(:,:,1), mg(lev)%e(:,:,1), jstart, jend)
  END DO
!$OMP BARRIER

END SUBROUTINE mgsmooth2d_omp

!-----------------------------------------------------------------------------------------------------------------------

RECURSIVE SUBROUTINE mgsmooth3d(lev, cyc)
  INTEGER, INTENT(IN) :: lev
  INTEGER, INTENT(IN) :: cyc

  INTEGER :: i, n

#ifdef PARALLEL_MPI
  INTEGER :: ierr
#endif

#ifdef PARALLEL_MPI
  IF (lev == glev) THEN
     CALL gather_h(mg(lev)%r(:,:,1), mg(lev+1)%r(:,:,1))
     IF (rank==0) CALL mgsmooth3d(lev+1, cyc)
     CALL mpi_barrier(comm, ierr) ! really required?
     CALL scatter_h(mg(lev+1)%e(:,:,1), mg(lev)%e(:,:,1))

     RETURN
  ELSE IF (lev == vlev) THEN
     CALL gather_v(mg(lev)%r(:,:,1), mg(lev+1)%r)
     IF (vrank==0) CALL mgsmooth3d(lev+1, cyc)
     CALL mpi_barrier(comm, ierr) ! really required?
     CALL scatter_v(mg(lev+1)%e, mg(lev)%e(:,:,1))

     RETURN
  END IF
#endif

  IF (cmode == 1  .OR. (cmode==2 .AND. cyc==0)) THEN
     n = 1
  ELSE
     n = 0
  END IF

MGPROFILE_BEGIN('apply1', lev)
  CALL apply3d1(mg(lev)%n1, mg(lev)%n2, mg(lev)%n3, mg(lev)%m, mg(lev)%r, mg(lev)%e, kindex=mg(lev)%mindx)
MGPROFILE_END('apply1', lev)

  DO i=0, n
     CALL sync3d(lev, mg(lev)%e)

MGPROFILE_BEGIN('apply2', lev)
     CALL apply3d2(mg(lev)%n1, mg(lev)%n2, mg(lev)%n3, mg(lev)%l, mg(lev)%r, mg(lev)%e, mg(lev)%t, kindex=mg(lev)%lindx)
MGPROFILE_END('apply2', lev)

#ifdef PARALLEL_MPI
     IF ((lev > vlev .OR. kpes==1) .AND. mg(lev+1)%n3==1) THEN
#else
     IF (mg(lev+1)%n3==1) THEN
#endif

MGPROFILE_BEGIN('coarse', lev)
        CALL coarse3d2d(mg(lev)%n1, mg(lev)%n2, mg(lev)%n3, mg(lev+1)%n1, mg(lev+1)%n2, &
                        mg(lev)%p1, mg(lev)%p2, mg(lev)%q1, mg(lev)%q2, mg(lev)%t, mg(lev+1)%r(:,:,1))
MGPROFILE_END('coarse', lev)

        CALL sync2d(lev+1, mg(lev+1)%r(:,:,1))

        CALL mgsmooth2d(lev+1, cyc+i)

MGPROFILE_BEGIN('fine', lev)
        CALL fine2d3d(mg(lev)%n1, mg(lev)%n2, mg(lev)%n3, mg(lev+1)%n1, mg(lev+1)%n2, &
                      mg(lev)%p1, mg(lev)%p2, mg(lev)%q1, mg(lev)%q2, mg(lev+1)%e(:,:,1), mg(lev)%e)
MGPROFILE_END('fine', lev)

     ELSE
MGPROFILE_BEGIN('coarse', lev)
        CALL coarse3d(mg(lev)%n1, mg(lev)%n2, mg(lev)%n3, mg(lev+1)%n1, mg(lev+1)%n2, mg(lev+1)%n3, &
                      mg(lev)%p1, mg(lev)%p2, mg(lev)%p3, mg(lev)%q1, mg(lev)%q2, mg(lev)%q3, mg(lev)%t, mg(lev+1)%r)
MGPROFILE_END('coarse', lev)

        CALL sync3d(lev+1, mg(lev+1)%r)

        CALL mgsmooth3d(lev+1, cyc+i)

MGPROFILE_BEGIN('fine', lev)
        CALL fine3d(mg(lev)%n1, mg(lev)%n2, mg(lev)%n3, mg(lev+1)%n1, mg(lev+1)%n2, mg(lev+1)%n3, &
                    mg(lev)%p1, mg(lev)%p2, mg(lev)%p3, mg(lev)%q1, mg(lev)%q2, mg(lev)%q3, mg(lev+1)%e, mg(lev)%e)
MGPROFILE_END('fine', lev)
     END IF

     CALL sync3d(lev, mg(lev)%e)

MGPROFILE_BEGIN('apply2', lev)
     CALL apply3d2(mg(lev)%n1, mg(lev)%n2, mg(lev)%n3, mg(lev)%l, mg(lev)%r, mg(lev)%e, mg(lev)%t, kindex=mg(lev)%lindx)
MGPROFILE_END('apply2', lev)

     CALL sync3d(lev, mg(lev)%t)

MGPROFILE_BEGIN('apply3', lev)
     CALL apply3d3(mg(lev)%n1, mg(lev)%n2, mg(lev)%n3, mg(lev)%m, mg(lev)%t, mg(lev)%e, kindex=mg(lev)%mindx)
MGPROFILE_END('apply3', lev)
  END DO

END SUBROUTINE mgsmooth3d

!-----------------------------------------------------------------------------------------------------------------------

RECURSIVE SUBROUTINE mgsmooth3d_omp(lev, cyc, tid)
  INTEGER, INTENT(IN) :: lev
  INTEGER, INTENT(IN) :: cyc
  INTEGER, INTENT(IN) :: tid

  INTEGER :: i, n

  INTEGER :: kstart, kend

#ifdef PARALLEL_MPI
  INTEGER :: ierr
#endif

  kstart = n3_start(tid,lev)
  kend   = n3_start(tid,lev) + n3_span(tid,lev) - 1

#ifdef PARALLEL_MPI
  IF (lev == glev) THEN
!$OMP SINGLE
     CALL gather_h(mg(lev)%r(:,:,1), mg(lev+1)%r(:,:,1))
!$OMP END SINGLE
     IF (rank==0) CALL mgsmooth3d_omp(lev+1, cyc, tid)
!$OMP SINGLE
     CALL mpi_barrier(comm, ierr)
     CALL scatter_h(mg(lev+1)%e(:,:,1), mg(lev)%e(:,:,1))
!$OMP END SINGLE
     RETURN
  ELSE IF (lev == vlev) THEN
!$OMP SINGLE
     CALL gather_v(mg(lev)%r(:,:,1), mg(lev+1)%r)
!$OMP END SINGLE
     IF (vrank==0) CALL mgsmooth3d_omp(lev+1, cyc, tid)
!$OMP SINGLE
     CALL mpi_barrier(comm, ierr)
     CALL scatter_v(mg(lev+1)%e, mg(lev)%e(:,:,1))
!$OMP END SINGLE
     RETURN
  END IF
#endif

  IF (cmode == 1  .OR. (cmode==2 .AND. cyc==0)) THEN
     n = 1
  ELSE
     n = 0
  END IF

MGPROFILE_BEGIN_THREAD('apply1', lev, tid)
  CALL apply3d1(mg(lev)%n1, mg(lev)%n2, mg(lev)%n3, mg(lev)%m, mg(lev)%r, mg(lev)%e, kstart, kend, kindex=mg(lev)%mindx)
MGPROFILE_END_THREAD('apply1', lev, tid)

  DO i=0, n
     CALL sync3d(lev, mg(lev)%e, tid)

MGPROFILE_BEGIN_THREAD('apply2', lev, tid)
     CALL apply3d2(mg(lev)%n1, mg(lev)%n2, mg(lev)%n3, mg(lev)%l, mg(lev)%r, mg(lev)%e, mg(lev)%t, kstart, kend, kindex=mg(lev)%lindx)
MGPROFILE_END_THREAD('apply2', lev, tid)
!$OMP BARRIER

#ifdef PARALLEL_MPI
     IF ((lev > vlev .OR. kpes==1) .AND. mg(lev+1)%n3==1) THEN
#else
     IF (mg(lev+1)%n3==1) THEN
#endif
MGPROFILE_BEGIN('coarse', lev)
        CALL coarse3d2d(mg(lev)%n1, mg(lev)%n2, mg(lev)%n3, mg(lev+1)%n1, mg(lev+1)%n2, &
                        mg(lev)%p1, mg(lev)%p2, mg(lev)%q1, mg(lev)%q2, mg(lev)%t, mg(lev+1)%r(:,:,1), &
                        n2_start(tid,lev+1), n2_start(tid,lev+1)+n2_span(tid,lev+1)-1)
MGPROFILE_END('coarse', lev)

        CALL sync2d(lev+1, mg(lev+1)%r(:,:,1), tid)

        CALL mgsmooth2d_omp(lev+1, cyc+i, tid)

MGPROFILE_BEGIN('fine', lev)
        CALL fine2d3d(mg(lev)%n1, mg(lev)%n2, mg(lev)%n3, mg(lev+1)%n1, mg(lev+1)%n2, &
                      mg(lev)%p1, mg(lev)%p2, mg(lev)%q1, mg(lev)%q2, mg(lev+1)%e(:,:,1), mg(lev)%e, kstart, kend)
MGPROFILE_END('fine', lev)
     ELSE
MGPROFILE_BEGIN_THREAD('coarse', lev, tid)
        CALL coarse3d(mg(lev)%n1, mg(lev)%n2, mg(lev)%n3, mg(lev+1)%n1, mg(lev+1)%n2, mg(lev+1)%n3, &
                      mg(lev)%p1, mg(lev)%p2, mg(lev)%p3, mg(lev)%q1, mg(lev)%q2, mg(lev)%q3, mg(lev)%t, mg(lev+1)%r, &
                      n3_start(tid,lev+1), n3_start(tid,lev+1)+n3_span(tid,lev+1)-1)
MGPROFILE_END_THREAD('coarse', lev, tid)

        CALL sync3d(lev+1, mg(lev+1)%r, tid)

        CALL mgsmooth3d_omp(lev+1, cyc+i, tid)

MGPROFILE_BEGIN_THREAD('fine', lev, tid)
        CALL fine3d(mg(lev)%n1, mg(lev)%n2, mg(lev)%n3, mg(lev+1)%n1, mg(lev+1)%n2, mg(lev+1)%n3, &
                    mg(lev)%p1, mg(lev)%p2, mg(lev)%p3, mg(lev)%q1, mg(lev)%q2, mg(lev)%q3, mg(lev+1)%e, mg(lev)%e, kstart, kend)
MGPROFILE_END_THREAD('fine', lev, tid)
     END IF

     CALL sync3d(lev, mg(lev)%e, tid)

MGPROFILE_BEGIN_THREAD('apply2', lev, tid)
     CALL apply3d2(mg(lev)%n1, mg(lev)%n2, mg(lev)%n3, mg(lev)%l, mg(lev)%r, mg(lev)%e, mg(lev)%t, kstart, kend, kindex=mg(lev)%lindx)
MGPROFILE_END_THREAD('apply2', lev, tid)
     CALL sync3d(lev, mg(lev)%t, tid)

MGPROFILE_BEGIN_THREAD('apply3', lev, tid)
     CALL apply3d3(mg(lev)%n1, mg(lev)%n2, mg(lev)%n3, mg(lev)%m, mg(lev)%t, mg(lev)%e, kstart, kend, kindex=mg(lev)%mindx)
MGPROFILE_END_THREAD('apply3', lev, tid)
  END DO
!$OMP BARRIER

END SUBROUTINE mgsmooth3d_omp

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE multigrid3d(r, e)
  REAL(8), INTENT(IN)    :: r(0:n1+1, 0:n2+1, 0:n3+1)
  REAL(8), INTENT(INOUT) :: e(0:n1+1, 0:n2+1, 0:n3+1)

MGPROFILE_BEGIN('apply1', 0)
  CALL apply3d1(n1, n2, n3, mg(0)%m, r, e, kindex=mg(0)%mindx)
MGPROFILE_END('apply1', 0)

  CALL sync3d(0, e)

MGPROFILE_BEGIN('apply2', 0)
  CALL apply3d2(n1, n2, n3, mg(0)%l, r, e, work, kindex=mg(0)%lindx)
MGPROFILE_END('apply2', 0)

#ifdef PARALLEL_MPI
  IF (kpes==1 .AND. mg(1)%n3==1) THEN
#else
  IF (mg(1)%n3==1) THEN
#endif

MGPROFILE_BEGIN('coarse', 0)
     CALL coarse3d2d(mg(0)%n1, mg(0)%n2, mg(0)%n3, mg(1)%n1, mg(1)%n2, &
                     mg(0)%p1, mg(0)%p2, mg(0)%q1, mg(0)%q2, work, mg(1)%r(:,:,1))
MGPROFILE_END('coarse', 0)

     CALL sync2d(1, mg(1)%r(:,:,1))

     CALL mgsmooth2d(1, 0)

MGPROFILE_BEGIN('fine', 0)
     CALL fine2d3d(mg(0)%n1, mg(0)%n2, mg(0)%n3, mg(1)%n1, mg(1)%n2, &
                   mg(0)%p1, mg(0)%p2, mg(0)%q1, mg(0)%q2, mg(1)%e(:,:,1), e)
MGPROFILE_END('fine', 0)
  ELSE
MGPROFILE_BEGIN('coarse', 0)
     CALL coarse3d(mg(0)%n1, mg(0)%n2, mg(0)%n3, mg(1)%n1, mg(1)%n2, mg(1)%n3, &
                   mg(0)%p1, mg(0)%p2, mg(0)%p3, mg(0)%q1, mg(0)%q2, mg(0)%q3, work, mg(1)%r)
MGPROFILE_END('coarse', 0)

     CALL sync3d(1, mg(1)%r)

     CALL mgsmooth3d(1, 0)

MGPROFILE_BEGIN('fine', 0)
     CALL fine3d(mg(0)%n1, mg(0)%n2, mg(0)%n3, mg(1)%n1, mg(1)%n2, mg(1)%n3, &
                 mg(0)%p1, mg(0)%p2, mg(0)%p3, mg(0)%q1, mg(0)%q2, mg(0)%q3, mg(1)%e, e)
MGPROFILE_END('fine', 0)
  END IF

  CALL sync3d(0, e)

MGPROFILE_BEGIN('apply2', 0)
  CALL apply3d2(n1, n2, n3, mg(0)%l, r, e, work, kindex=mg(0)%lindx)
MGPROFILE_END('apply2', 0)

  CALL sync3d(0, work)

MGPROFILE_BEGIN('apply3', 0)
  CALL apply3d3(n1, n2, n3, mg(0)%m, work, e, kindex=mg(0)%mindx)
MGPROFILE_END('apply3', 0)

END SUBROUTINE multigrid3d

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE multigrid3d_omp(r, e, tid)
  REAL(8), INTENT(IN)    :: r(0:n1+1, 0:n2+1, 0:n3+1)
  REAL(8), INTENT(INOUT) :: e(0:n1+1, 0:n2+1, 0:n3+1)
  INTEGER, INTENT(IN)    :: tid

  INTEGER :: kstart, kend

  kstart = n3_start(tid,0)
  kend   = n3_start(tid,0) + n3_span(tid,0) - 1

MGPROFILE_BEGIN_THREAD('apply1', 0, tid)
  CALL apply3d1(n1, n2, n3, mg(0)%m, r, e, kstart, kend, kindex=mg(0)%mindx)
MGPROFILE_END_THREAD('apply1', 0, tid)

  CALL sync3d(0, e, tid)

MGPROFILE_BEGIN_THREAD('apply2', 0, tid)
  CALL apply3d2(n1, n2, n3, mg(0)%l, r, e, work, kstart, kend, kindex=mg(0)%lindx)
MGPROFILE_END_THREAD('apply2', 0, tid)
!$OMP BARRIER

#ifdef PARALLEL_MPI
  IF (kpes==1 .AND. mg(1)%n3==1) THEN
#else
  IF (mg(1)%n3==1) THEN
#endif

MGPROFILE_BEGIN_THREAD('coarse', 0, tid)
     CALL coarse3d2d(mg(0)%n1, mg(0)%n2, mg(0)%n3, mg(1)%n1, mg(1)%n2, &
                     mg(0)%p1, mg(0)%p2, mg(0)%q1, mg(0)%q2, work, mg(1)%r(:,:,1), &
                     n2_start(tid,1), n2_start(tid,1)+n2_span(tid,1)-1)
MGPROFILE_END_THREAD('coarse', 0, tid)

     CALL sync2d(1, mg(1)%r(:,:,1), tid)

     CALL mgsmooth2d_omp(1, 0, tid)

MGPROFILE_BEGIN_THREAD('fine', 0, tid)
     CALL fine2d3d(mg(0)%n1, mg(0)%n2, mg(0)%n3, mg(1)%n1, mg(1)%n2, &
                   mg(0)%p1, mg(0)%p2, mg(0)%q1, mg(0)%q2, mg(1)%e(:,:,1), e, kstart, kend)
MGPROFILE_END_THREAD('fine', 0, tid)

  ELSE
MGPROFILE_BEGIN_THREAD('coarse', 0, tid)
     CALL coarse3d(mg(0)%n1, mg(0)%n2, mg(0)%n3, mg(1)%n1, mg(1)%n2, mg(1)%n3, &
                   mg(0)%p1, mg(0)%p2, mg(0)%p3, mg(0)%q1, mg(0)%q2, mg(0)%q3, work, mg(1)%r, &
                   n3_start(tid,1), n3_start(tid,1)+n3_span(tid,1)-1)

MGPROFILE_END_THREAD('coarse', 0, tid)

     CALL sync3d(1, mg(1)%r, tid)

     CALL mgsmooth3d_omp(1, 0, tid)

MGPROFILE_BEGIN_THREAD('fine', 0, tid)
     CALL fine3d(mg(0)%n1, mg(0)%n2, mg(0)%n3, mg(1)%n1, mg(1)%n2, mg(1)%n3, &
                 mg(0)%p1, mg(0)%p2, mg(0)%p3, mg(0)%q1, mg(0)%q2, mg(0)%q3, mg(1)%e, e, kstart, kend)
MGPROFILE_END_THREAD('fine', 0, tid)
  END IF

  CALL sync3d(0, e, tid)

MGPROFILE_BEGIN_THREAD('apply2', 0, tid)
  CALL apply3d2(n1, n2, n3, mg(0)%l, r, e, work, kstart, kend, kindex=mg(0)%lindx)
MGPROFILE_END_THREAD('apply2', 0, tid)

  CALL sync3d(0, work, tid)

MGPROFILE_BEGIN_THREAD('apply3', 0, tid)
  CALL apply3d3(n1, n2, n3, mg(0)%m, work, e, kstart, kend, kindex=mg(0)%mindx)
MGPROFILE_END_THREAD('apply3', 0, tid)

END SUBROUTINE multigrid3d_omp

!-----------------------------------------------------------------------------------------------------------------------

#ifdef PARALLEL_MPI

SUBROUTINE gather_h_real4(in, out, default)
  REAL(4), INTENT(IN)  :: in( 0:mg(glev  )%n1+1,0:mg(glev  )%n2+1)
  REAL(4), INTENT(OUT) :: out(0:mg(glev+1)%n1+1,0:mg(glev+1)%n2+1)
  REAL(4), INTENT(IN), OPTIONAL :: default

  REAL(4) :: default_
  REAL(4) :: tmp(MG_GDIM,MG_GDIM)
  INTEGER :: i, j, l, m, n
  INTEGER :: ierr

MGPROFILE_BEGIN('gath',0)
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

MGPROFILE_END('gath',0)

END SUBROUTINE gather_h_real4

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE gather_h_real8(in, out, default)
  REAL(8), INTENT(IN)  :: in( 0:mg(glev  )%n1+1,0:mg(glev  )%n2+1)
  REAL(8), INTENT(OUT) :: out(0:mg(glev+1)%n1+1,0:mg(glev+1)%n2+1)
  REAL(8), INTENT(IN), OPTIONAL :: default

  REAL(8) :: default_
  REAL(8) :: tmp(MG_GDIM,MG_GDIM)
  INTEGER :: i, j, l, m, n
  INTEGER :: ierr

MGPROFILE_BEGIN('gath',0)
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

MGPROFILE_END('gath',0)

END SUBROUTINE gather_h_real8

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE gather_v_real4(in, out, default)
  REAL(4), INTENT(IN)  :: in( 0:mg(vlev)%n1+1, 0:mg(vlev)%n2+1)
  REAL(4), INTENT(OUT) :: out(0:mg(vlev)%n1+1, 0:mg(vlev)%n2+1, 0:kpes+1)
  REAL(4), INTENT(IN), OPTIONAL :: default

  REAL(4) :: default_
  INTEGER :: k
  INTEGER :: req(0:kpes-1)
  INTEGER :: ierr

  default_ = 0.0
  IF (present(default)) default_ = default

  IF (vrank==0) THEN
     DO k=0, kpes-1
        IF (ranks(icoord, jcoord, k)==MPI_PROC_NULL) THEN
           out(:,:,k+1) = default_
           req(k) = MPI_REQUEST_NULL
        ELSE
           CALL mpi_irecv(out(0,0,k+1), (mg(vlev)%n1+2)*(mg(vlev)%n2+2), MPI_REAL4, vranks(k), 0, vcomm, req(k), ierr)
        END IF
     END DO
  END IF

  CALL mpi_send(in(0,0), (mg(vlev)%n1+2)*(mg(vlev)%n2+2), MPI_REAL4, 0, 0, vcomm, ierr)

  IF (vrank==0) CALl mpi_waitall(kpes, req, MPI_STATUSES_IGNORE, ierr)
END SUBROUTINE gather_v_real4

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE gather_v_real8(in, out, default)
  REAL(8), INTENT(IN)  :: in( 0:mg(vlev)%n1+1, 0:mg(vlev)%n2+1)
  REAL(8), INTENT(OUT) :: out(0:mg(vlev)%n1+1, 0:mg(vlev)%n2+1, 0:kpes+1)
  REAL(8), INTENT(IN), OPTIONAL :: default

  REAL(8) :: default_
  INTEGER :: k
  INTEGER :: req(0:kpes-1)
  INTEGER :: ierr

  default_ = 0.0
  IF (present(default)) default_ = default

  IF (vrank==0) THEN
     DO k=0, kpes-1
        IF (ranks(icoord, jcoord, k)==MPI_PROC_NULL) THEN
           out(:,:,k+1) = default_
           req(k) = MPI_REQUEST_NULL
        ELSE
           CALL mpi_irecv(out(0,0,k+1), (mg(vlev)%n1+2)*(mg(vlev)%n2+2), MPI_REAL8, vranks(k), 0, vcomm, req(k), ierr)
        END IF
     END DO
  END IF

  CALL mpi_send(in(0,0), (mg(vlev)%n1+2)*(mg(vlev)%n2+2), MPI_REAL8, 0, 0, vcomm, ierr)

  IF (vrank==0) CALl mpi_waitall(kpes, req, MPI_STATUSES_IGNORE, ierr)

END SUBROUTINE gather_v_real8

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE scatter_h_real4(in, out)
  REAL(4), INTENT(IN)    :: in( 0:mg(glev+1)%n1+1, 0:mg(glev+1)%n2+1)
  REAL(4), INTENT(INOUT) :: out(0:mg(glev  )%n1+1, 0:mg(glev  )%n2+1)

  REAL(4) :: tmp(MG_GDIM,MG_GDIM)
  INTEGER :: i, j, l, m, n
  INTEGER :: ierr

MGPROFILE_BEGIN('scat',0)

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

MGPROFILE_END('scat',0)

END SUBROUTINE scatter_h_real4

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE scatter_h_real8(in, out)
  REAL(8), INTENT(IN)    :: in( 0:mg(glev+1)%n1+1,0:mg(glev+1)%n2+1)
  REAL(8), INTENT(INOUT) :: out(0:mg(glev  )%n1+1,0:mg(glev  )%n2+1)

  REAL(8) :: tmp(MG_GDIM,MG_GDIM)
  INTEGER :: i, j, l, m, n
  INTEGER :: ierr

MGPROFILE_BEGIN('scat',0)

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

MGPROFILE_END('scat',0)

END SUBROUTINE scatter_h_real8

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE scatter_v_real4(in, out)
  REAL(4), INTENT(IN)  :: in( 0:mg(vlev)%n1+1, 0:mg(vlev)%n2+1, 0:kpes+1)
  REAL(4), INTENT(OUT) :: out(0:mg(vlev)%n1+1, 0:mg(vlev)%n2+1)

  INTEGER :: k
  INTEGER :: req(0:kpes-1)
  INTEGER :: ierr

  IF (vrank==0) THEN
     DO k=0, kpes-1
        CALL mpi_isend(in(0,0,k+1), (mg(vlev)%n1+2)*(mg(vlev)%n2+2), MPI_REAL4, vranks(k), 0, vcomm, req(k), ierr)
     END DO
  END IF

  CALL mpi_recv(out(0,0), (mg(vlev)%n1+2)*(mg(vlev)%n2+2), MPI_REAL4, 0, 0, vcomm, MPI_STATUS_IGNORE, ierr)

  IF (vrank==0) CALl mpi_waitall(kpes, req, MPI_STATUSES_IGNORE, ierr)
END SUBROUTINE scatter_v_real4

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE scatter_v_real8(in, out)
  REAL(8), INTENT(IN)  :: in( 0:mg(vlev)%n1+1, 0:mg(vlev)%n2+1, 0:kpes+1)
  REAL(8), INTENT(OUT) :: out(0:mg(vlev)%n1+1, 0:mg(vlev)%n2+1)

  INTEGER :: k
  INTEGER :: req(0:kpes-1)
  INTEGER :: ierr

  IF (vrank==0) THEN
     DO k=0, kpes-1
        CALL mpi_isend(in(0,0,k+1), (mg(vlev)%n1+2)*(mg(vlev)%n2+2), MPI_REAL8, vranks(k), 0, vcomm, req(k), ierr)
     END DO
  END IF

  CALL mpi_recv(out(0,0), (mg(vlev)%n1+2)*(mg(vlev)%n2+2), MPI_REAL8, 0, 0, vcomm, MPI_STATUS_IGNORE, ierr)


  IF (vrank==0) CALl mpi_waitall(kpes, req, MPI_STATUSES_IGNORE, ierr)
END SUBROUTINE scatter_v_real8

#endif

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE sync2d_begin_real4(lev, buf, tid)
  INTEGER, INTENT(IN)    :: lev
  REAL(4), INTENT(INOUT) :: buf(0:mg(lev)%n1+1, 0:mg(lev)%n2+1)
  INTEGER, INTENT(IN), OPTIONAL :: tid

  INTEGER :: j, jstart, jend

#ifdef PARALLEL_MPI
  INTEGER :: ierr
#endif

  IF (present(tid)) THEN
     jstart = n2_start(tid,lev)
     jend   = n2_start(tid,lev) + n2_span(tid,lev)-1
  ELSE
     jstart = 1
     jend   = mg(lev)%n2
  END IF
!$OMP BARRIER

#ifdef PARALLEL_MPI
  IF (lev > glev) RETURN

!$OMP SINGLE
  CALL mpi_startall(4, req2d_r4(0:3,lev), ierr)
!$OMP END SINGLE nowait

  DO j=jstart, jend
     sendbuf_e_r4(j) = buf(mg(lev)%n1,j)
     sendbuf_w_r4(j) = buf(1,j)
  END DO

!$OMP SECTIONS
  sendbuf_n_r4(1:mg(lev)%n1) = buf(1:mg(lev)%n1,mg(lev)%n2)
!$OMP SECTION
  sendbuf_s_r4(1:mg(lev)%n1) = buf(1:mg(lev)%n1,1)
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
     sendbuf_w_r8(j) = buf(1,j)
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
           buf(0,   j) = buf(n1,j)
           buf(n1+1,j) = buf(1, j)
        END DO
     ELSE
        DO j=jstart, jend
           buf(0,   j) = 0.0
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
MGPROFILE_BEGIN('sync',lev)
     CALL mpi_waitall(8, req2d_r4(0:7,lev), MPI_STATUSES_IGNORE, ierr)
MGPROFILE_END('sync',lev)
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
        buf(0,j)    = buf(n1,j)
        buf(n1+1,j) = buf(1,j)
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
           buf(0,   j) = 0.0D0
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
MGPROFILE_BEGIN('sync',lev)
     CALL mpi_waitall(8, req2d_r8(0:7,lev), MPI_STATUSES_IGNORE, ierr)
MGPROFILE_END('sync',lev)
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
        buf(0,   j) = buf(n1,j)
        buf(n1+1,j) = buf(1,j)
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

SUBROUTINE sync3d_begin_real4(lev, buf, tid)
  INTEGER, INTENT(IN)    :: lev
  REAL(4), INTENT(INOUT) :: buf(0:mg(lev)%n1+1, 0:mg(lev)%n2+1, 0:mg(lev)%n3+1)
  INTEGER, INTENT(IN), OPTIONAL :: tid

  INTEGER :: n1, n2, n3
  INTEGER :: i, j, k
  INTEGER :: kstart, kend, jstart, jend

#ifdef PARALLEL_MPI
  INTEGER :: ierr
#endif

  n1 = mg(lev)%n1
  n2 = mg(lev)%n2
  n3 = mg(lev)%n3

  IF (present(tid)) THEN
     kstart = n3_start(tid,lev)
     kend   = n3_start(tid,lev) + n3_span(tid,lev)-1

     jstart = n2_start(tid,lev)
     jend   = n2_start(tid,lev) + n2_span(tid,lev)-1
  ELSE
     kstart = 1
     kend   = n3

     jstart = 1
     jend   = n2
  END IF
!$OMP BARRIER

#ifdef PARALLEL_MPI
  IF (lev > glev) RETURN

!$OMP SINGLE
#ifdef PARALLEL3D
  IF (lev <= vlev) THEN
     CALL mpi_startall(6, req3d_r4(0:5,lev), ierr)
  ELSE
     IF (vrank==0) CALL mpi_startall(4, req3d_r4(0:3,lev), ierr)
  END IF
#else
  CALL mpi_startall(4, req3d_r4(0:3,lev), ierr)
#endif
!$OMP END SINGLE nowait

#ifdef PARALLEL3D
  DO j=jstart, jend
     DO i=1, n1
        sendbuf_l_r4(n1*(j-1) + i) = buf(i,j,1)
     END DO
  END DO
#endif

  DO k=kstart, kend
     DO j=1, n2
        sendbuf_e_r4(n2*(k-1) + j) = buf(n1,j,k)
        sendbuf_w_r4(n2*(k-1) + j) = buf(1, j,k)
     END DO
     DO i=1, n1
        sendbuf_n_r4(n1*(k-1) + i) = buf(i,n2,k)
        sendbuf_s_r4(n1*(k-1) + i) = buf(i,1, k)
     END DO
  END DO

#ifdef PARALLEL3D
  DO j=jstart, jend
     DO i=1, n1
        sendbuf_u_r4(n1*(j-1) + i) = buf(i,j,n3)
     END DO
  END DO
#endif

!$OMP BARRIER
!$OMP SINGLE
#ifdef PARALLEL3D
  IF (lev <= vlev) THEN
     CALL mpi_startall(6, req3d_r4(6:11,lev), ierr)
  ELSE
     IF (vrank==0) CALL mpi_startall(4, req3d_r4(6:9, lev), ierr)
  END IF
#else
  CALL mpi_startall(4, req3d_r4(6:9,lev), ierr)
#endif
!$OMP END SINGLE
#endif

END SUBROUTINE sync3d_begin_real4

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE sync3d_begin_real8(lev, buf, tid)
  INTEGER, INTENT(IN)    :: lev
  REAL(8), INTENT(INOUT) :: buf(0:mg(lev)%n1+1, 0:mg(lev)%n2+1, 0:mg(lev)%n3+1)
  INTEGER, INTENT(IN), OPTIONAL :: tid

  INTEGER :: n1, n2, n3
  INTEGER :: i, j, k
  INTEGER :: kstart, kend, jstart, jend

#ifdef PARALLEL_MPI
  INTEGER :: ierr
#endif

  n1 = mg(lev)%n1
  n2 = mg(lev)%n2
  n3 = mg(lev)%n3

  IF (present(tid)) THEN
     kstart = n3_start(tid,lev)
     kend   = n3_start(tid,lev) + n3_span(tid,lev)-1

     jstart = n2_start(tid,lev)
     jend   = n2_start(tid,lev) + n2_span(tid,lev)-1
  ELSE
     kstart = 1
     kend   = n3

     jstart = 1
     jend   = n2
  END IF
!$OMP BARRIER

#ifdef PARALLEL_MPI
  IF (lev > glev) RETURN

!$OMP SINGLE
#ifdef PARALLEL3D
  IF (lev <= vlev) THEN
     CALL mpi_startall(6, req3d_r8(0:5,lev), ierr)
  ELSE
     IF (vrank==0) CALL mpi_startall(4, req3d_r8(0:3,lev), ierr)
  END IF
#else
  CALL mpi_startall(4, req3d_r8(0:3,lev), ierr)
#endif
!$OMP END SINGLE nowait

#ifdef PARALLEL3D
  DO j=jstart, jend
     DO i=1, n1
        sendbuf_l_r8(n1*(j-1) + i) = buf(i,j,1)
     END DO
  END DO
#endif

  DO k=kstart, kend
     DO j=1, n2
        sendbuf_e_r8(n2*(k-1) + j) = buf(n1,j,k)
        sendbuf_w_r8(n2*(k-1) + j) = buf(1, j,k)
     END DO
     DO i=1, n1
        sendbuf_n_r8(n1*(k-1) + i) = buf(i,n2,k)
        sendbuf_s_r8(n1*(k-1) + i) = buf(i,1, k)
     END DO
  END DO

#ifdef PARALLEL3D
  DO j=jstart, jend
     DO i=1, n1
        sendbuf_u_r8(n1*(j-1) + i) = buf(i,j,n3)
     END DO
  END DO
#endif

!$OMP BARRIER
!$OMP SINGLE
#ifdef PARALLEL3D
  IF (lev <= vlev) THEN
     CALL mpi_startall(6, req3d_r8(6:11,lev), ierr)
  ELSE
     IF (vrank==0) CALL mpi_startall(4, req3d_r8(6:9, lev), ierr)
  END IF
#else
  CALL mpi_startall(4, req3d_r8(6:9,lev), ierr)
#endif
!$OMP END SINGLE
#endif

END SUBROUTINE sync3d_begin_real8

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE sync3d_end_real4(lev, buf, tid)
  INTEGER, INTENT(IN)    :: lev
  REAL(4), INTENT(INOUT) :: buf(0:mg(lev)%n1+1, 0:mg(lev)%n2+1, 0:mg(lev)%n3+1)
  INTEGER, INTENT(IN), OPTIONAL :: tid

  INTEGER :: n1, n2, n3
  INTEGER :: i, j, k
  INTEGER :: kstart, kend, jstart, jend

#ifdef PARALLEL_MPI
  INTEGER :: ierr
#endif

  n1 = mg(lev)%n1
  n2 = mg(lev)%n2
  n3 = mg(lev)%n3

  IF (present(tid)) THEN
     kstart = n3_start(tid,lev)
     kend   = n3_start(tid,lev) + n3_span(tid,lev) - 1

     jstart = n2_start(tid,lev)
     jend   = n2_start(tid,lev) + n2_span(tid,lev) - 1
  ELSE
     kstart = 1
     kend   = n3

     jstart = 1
     jend   = n2
  END IF

#ifdef PARALLEL_MPI
  IF (lev > glev) THEN
     IF (period1) THEN
        buf(0,   1:n2,1:n3) = buf(n1,1:n2,1:n3)
        buf(n1+1,1:n2,1:n3) = buf(1, 1:n2,1:n3)
     ELSE
        buf(0,   1:n2,1:n3) = 0.0
        buf(n1+1,1:n2,1:n3) = 0.0
     END IF

     IF (period2) THEN
        buf(1:n1,0,   1:n3) = buf(1:n1,n2,1:n3)
        buf(1:n1,n2+1,1:n3) = buf(1:n1,1, 1:n3)
     ELSE
        buf(1:n1,0,   1:n3) = 0.0
        buf(1:n1,n2+1,1:n3) = 0.0
     END IF

     IF (period3) THEN
        buf(1:n1,1:n2,0)    = buf(1:n1,1:n2,n3+1)
        buf(1:n1,1:n2,n3+1) = buf(1:n1,1:n2,1)
     END IF
  ELSE
!$OMP SINGLE
MGPROFILE_BEGIN('sync',lev)
     CALL mpi_waitall(12, req3d_r4(0:11,lev), MPI_STATUSES_IGNORE, ierr)
MGPROFILE_END('sync',lev)
!$OMP END SINGLE

     IF (mg(lev)%rank_w /= MPI_PROC_NULL) THEN
        DO k=kstart, kend
           DO j=1, n2
              buf(0,j,k) = recvbuf_w_r4(n2*(k-1) + j)
           END DO
        END DO
     END IF

     IF (mg(lev)%rank_e /= MPI_PROC_NULL) THEN
        DO k=kstart, kend
           DO j=1, n2
              buf(n1+1,j,k) = recvbuf_e_r4(n2*(k-1) + j)
           END DO
        END DO
     END IF

     IF (mg(lev)%rank_s /= MPI_PROC_NULL) THEN
        DO k=kstart, kend
           DO i=1, n1
              buf(i,0,k) = recvbuf_s_r4(n1*(k-1) + i)
           END DO
        END DO
     END IF

     IF (mg(lev)%rank_n /= MPI_PROC_NULL) THEN
        DO k=kstart, kend
           DO i=1, n1
              buf(i,n2+1,k) = recvbuf_n_r4(n1*(k-1) + i)
           END DO
        END DO
     END IF

#ifdef PARALLEL3D
     IF (lev <= vlev) THEN
        IF (mg(lev)%rank_l /= MPI_PROC_NULL) THEN
           DO j=jstart, jend
              DO i=1, n1
                 buf(i,j,0) = recvbuf_l_r4(n1*(j-1) + i)
              END DO
           END DO
        END IF

        IF (mg(lev)%rank_u /= MPI_PROC_NULL) THEN
           DO j=jstart, jend
              DO i=1, n1
                 buf(i,j,n3+1) = recvbuf_u_r4(n1*(j-1) + i)
              END DO
           END DO
        END IF
     ELSE
        IF (period3) THEN
           buf(1:n1,1:n2,0)    = buf(1:n1,1:n2,n3+1)
           buf(1:n1,1:n2,n3+1) = buf(1:n1,1:n2,1)
        END IF
     END IF
#else
     IF (period3) THEN
        buf(1:n1,1:n2,0)    = buf(1:n1,1:n2,n3+1)
        buf(1:n1,1:n2,n3+1) = buf(1:n1,1:n2,1)
     END IF
#endif
  END IF
#else
  IF (period1) THEN
     buf(0,   1:n2,kstart:kend) = buf(n1,1:n2,kstart:kend)
     buf(n1+1,1:n2,kstart:kend) = buf(1, 1:n2,kstart:kend)
  END IF

  IF (period2) THEN
     buf(1:n1,0,   kstart:kend) = buf(1:n1,n2,kstart:kend)
     buf(1:n1,n2+1,kstart:kend) = buf(1:n1,1, kstart:kend)
  END IF

  IF (period3) THEN
     buf(1:n1,jstart:jend,0)    = buf(1:n1,jstart:jend, n3)
     buf(1:n1,jstart:jend,n3+1) = buf(1:n1,jstart:jend, 1)
  END IF
#endif
!$OMP BARRIER
END SUBROUTINE sync3d_end_real4

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE sync3d_end_real8(lev, buf, tid)
  INTEGER, INTENT(IN)    :: lev
  REAL(8), INTENT(INOUT) :: buf(0:mg(lev)%n1+1, 0:mg(lev)%n2+1, 0:mg(lev)%n3+1)
  INTEGER, INTENT(IN), OPTIONAL :: tid

  INTEGER :: n1, n2, n3
  INTEGER :: i, j, k
  INTEGER :: kstart, kend, jstart, jend

#ifdef PARALLEL_MPI
  INTEGER :: ierr
#endif

  n1 = mg(lev)%n1
  n2 = mg(lev)%n2
  n3 = mg(lev)%n3

  IF (present(tid)) THEN
     kstart = n3_start(tid,lev)
     kend   = n3_start(tid,lev) + n3_span(tid,lev) - 1

     jstart = n2_start(tid,lev)
     jend   = n2_start(tid,lev) + n2_span(tid,lev)-1
  ELSE
     kstart = 1
     kend   = n3

     jstart = 1
     jend   = n2
  END IF

#ifdef PARALLEL_MPI
  IF (lev > glev) THEN
     IF (period1) THEN
        buf(0,   1:n2,1:n3) = buf(n1,1:n2,1:n3)
        buf(n1+1,1:n2,1:n3) = buf(1, 1:n2,1:n3)
     ELSE
        buf(0,   1:n2,1:n3) = 0.0D0
        buf(n1+1,1:n2,1:n3) = 0.0D0
     END IF

     IF (period2) THEN
        buf(1:n1,0,   1:n3) = buf(1:n1,n2,1:n3)
        buf(1:n1,n2+1,1:n3) = buf(1:n1,1, 1:n3)
     ELSE
        buf(1:n1,0,   1:n3) = 0.0D0
        buf(1:n1,n2+1,1:n3) = 0.0D0
     END IF

     IF (period3) THEN
        buf(1:n1,1:n2,0)    = buf(1:n1,1:n2,n3+1)
        buf(1:n1,1:n2,n3+1) = buf(1:n1,1:n2,1)
     END IF
  ELSE
!$OMP SINGLE
MGPROFILE_BEGIN('sync',lev)
     CALL mpi_waitall(12, req3d_r8(0:11,lev), MPI_STATUSES_IGNORE, ierr)
MGPROFILE_END('sync',lev)
!$OMP END SINGLE

     IF (mg(lev)%rank_w /= MPI_PROC_NULL) THEN
        DO k=kstart, kend
           DO j=1, n2
              buf(0,j,k) = recvbuf_w_r8(n2*(k-1) + j)
           END DO
        END DO
     END IF

     IF (mg(lev)%rank_e /= MPI_PROC_NULL) THEN
        DO k=kstart, kend
           DO j=1, n2
              buf(n1+1,j,k) = recvbuf_e_r8(n2*(k-1) + j)
           END DO
        END DO
     END IF

     IF (mg(lev)%rank_s /= MPI_PROC_NULL) THEN
        DO k=kstart, kend
           DO i=1, n1
              buf(i,0,k) = recvbuf_s_r8(n1*(k-1) + i)
           END DO
        END DO
     END IF

     IF (mg(lev)%rank_n /= MPI_PROC_NULL) THEN
        DO k=kstart, kend
           DO i=1, n1
              buf(i,n2+1,k) = recvbuf_n_r8(n1*(k-1) + i)
           END DO
        END DO
     END IF

#ifdef PARALLEL3D
     IF (lev <= vlev) THEN
        IF (mg(lev)%rank_l /= MPI_PROC_NULL) THEN
           DO j=jstart, jend
              DO i=1, n1
                 buf(i,j,0) = recvbuf_l_r8(n1*(j-1) + i)
              END DO
           END DO
        END IF

        IF (mg(lev)%rank_u /= MPI_PROC_NULL) THEN
           DO j=jstart, jend
              DO i=1, n1
                 buf(i,j,n3+1) = recvbuf_u_r8(n1*(j-1) + i)
              END DO
           END DO
        END IF
     ELSE
        IF (period3) THEN
           buf(1:n1,1:n2,0)    = buf(1:n1,1:n2,n3+1)
           buf(1:n1,1:n2,n3+1) = buf(1:n1,1:n2,1)
        END IF
     END IF
#else
     IF (period3) THEN
        buf(1:n1,1:n2,0)    = buf(1:n1,1:n2,n3+1)
        buf(1:n1,1:n2,n3+1) = buf(1:n1,1:n2,1)
     END IF
#endif
  END IF
#else
  IF (period1) THEN
     buf(0,   1:n2,kstart:kend) = buf(n1,1:n2,kstart:kend)
     buf(n1+1,1:n2,kstart:kend) = buf(1, 1:n2,kstart:kend)
  END IF

  IF (period2) THEN
     buf(1:n1,0,   kstart:kend) = buf(1:n1,n2,kstart:kend)
     buf(1:n1,n2+1,kstart:kend) = buf(1:n1,1, kstart:kend)
  END IF

  IF (period3) THEN
     buf(1:n1,jstart:jend,0)    = buf(1:n1,jstart:jend, n3)
     buf(1:n1,jstart:jend,n3+1) = buf(1:n1,jstart:jend, 1)
  END IF
#endif
!$OMP BARRIER
END SUBROUTINE sync3d_end_real8

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE sync3d_real4(lev, buf, tid)
  INTEGER, INTENT(IN)    :: lev
  REAL(4), INTENT(INOUT) :: buf(0:mg(lev)%n1+1, 0:mg(lev)%n2+1, 0:mg(lev)%n3+1)
  INTEGER, INTENT(IN), OPTIONAL :: tid

  IF (present(tid)) THEN
     CALL sync3d_begin(lev, buf, tid)
     CALL sync3d_end(lev, buf, tid)
  ELSE
     CALL sync3d_begin(lev, buf)
     CALL sync3d_end(lev, buf)
  END IF

END SUBROUTINE sync3d_real4

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE sync3d_real8(lev, buf, tid)
  INTEGER, INTENT(IN)    :: lev
  REAL(8), INTENT(INOUT) :: buf(0:mg(lev)%n1+1, 0:mg(lev)%n2+1, 0:mg(lev)%n3+1)
  INTEGER, INTENT(IN), OPTIONAL :: tid

  IF (present(tid)) THEN
     CALL sync3d_begin(lev, buf, tid)
     CALL sync3d_end(lev, buf, tid)
  ELSE
     CALL sync3d_begin(lev, buf)
     CALL sync3d_end(lev, buf)
  END IF

END SUBROUTINE sync3d_real8

!-----------------------------------------------------------------------------------------------------------------------

SUBROUTINE get_coefficient(row, data)
  INTEGER, INTENT(IN) :: row
  REAL(8), INTENT(OUT) :: data(:,:,:)

  INTEGER :: i, j, k

  CALL assert(initialized, "solver3d is not initialized")
  CALL assert(row >= -3 .AND. row <= 3, "invalid row")

  DO k=1, n3
  DO j=1, n2
  DO i=1, n1
     data(i,j,k) = l3d(row,i,j,k)
  END DO
  END DO
  END DO

END SUBROUTINE get_coefficient

!-----------------------------------------------------------------------------------------------------------------------

END MODULE solver3d
