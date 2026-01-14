#include "macro.h"

MODULE velocity
  USE misc
  USE geometry
  USE io
  IMPLICIT NONE

  REAL(8), ALLOCATABLE :: u(:,:,:)
  REAL(8), ALLOCATABLE :: v(:,:,:)
  REAL(8), ALLOCATABLE :: w(:,:,:)

  REAL(8), ALLOCATABLE :: u_old(:,:,:)
  REAL(8), ALLOCATABLE :: v_old(:,:,:)
  REAL(8), ALLOCATABLE :: w_old(:,:,:)

  REAL(8), ALLOCATABLE :: gx_ab(:,:,:)
  REAL(8), ALLOCATABLE :: gy_ab(:,:,:)
  REAL(8), ALLOCATABLE :: gz_ab(:,:,:)

  REAL(8), ALLOCATABLE :: gx(:,:,:)
  REAL(8), ALLOCATABLE :: gy(:,:,:)
  REAL(8), ALLOCATABLE :: gz(:,:,:)

  REAL(8), ALLOCATABLE :: p2d(:,:)
  REAL(8), ALLOCATABLE :: p3d(:,:,:)
  REAL(8), ALLOCATABLE :: q2d(:,:)
  REAL(8), ALLOCATABLE :: q3d(:,:,:)

  REAL(8), ALLOCATABLE :: dudx(:,:,:)
  REAL(8), ALLOCATABLE :: dudy(:,:,:)
  REAL(8), ALLOCATABLE :: dudz(:,:,:)

  REAL(8), ALLOCATABLE :: dvdx(:,:,:)
  REAL(8), ALLOCATABLE :: dvdy(:,:,:)
  REAL(8), ALLOCATABLE :: dvdz(:,:,:)

  REAL(8), ALLOCATABLE :: dwdx(:,:,:)
  REAL(8), ALLOCATABLE :: dwdy(:,:,:)
  REAL(8), ALLOCATABLE :: dwdz(:,:,:)


  REAL(8), ALLOCATABLE :: ke(:,:,:)
  REAL(8), ALLOCATABLE :: srate(:,:,:)  ! strain rate [s^-1]
  REAL(8), ALLOCATABLE :: sfreq2(:,:,:) ! square of shear frequency [s^-2]

  REAL(8), ALLOCATABLE :: taux_sfc(:,:)
  REAL(8), ALLOCATABLE :: tauy_sfc(:,:)

  REAL(8), ALLOCATABLE :: taux_btm(:,:)
  REAL(8), ALLOCATABLE :: tauy_btm(:,:)

CONTAINS
  SUBROUTINE init_velocity
    USE parameters
    INTEGER :: i, j, k
    LOGICAL :: stat, stat_w

    REAL(8) :: tmp_xz(isize,ksize)
    REAL(8) :: tmp_yz(jsize,ksize)

    ALLOCATE(u( -slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv))
    ALLOCATE(v(1-slv:isize+slv,  -slv:jsize+slv, 1-slv:ksize+slv))
    ALLOCATE(w(1-slv:isize+slv, 1-slv:jsize+slv,  -slv:ksize+slv))

    ALLOCATE(u_old( -slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv))
    ALLOCATE(v_old(1-slv:isize+slv,  -slv:jsize+slv, 1-slv:ksize+slv))
    ALLOCATE(w_old(1-slv:isize+slv, 1-slv:jsize+slv,  -slv:ksize+slv))

!$OMP PARALLEL WORKSHARE
    u(:,:,:) = 0.0
    v(:,:,:) = 0.0
    w(:,:,:) = 0.0

    u_old(:,:,:) = UNDEF
    v_old(:,:,:) = UNDEF
    w_old(:,:,:) = UNDEF
!$OMP END PARALLEL WORKSHARE

    IF (.NOT. perfect_restart .AND. has_initial('PSI')) THEN
       CALL initial_stream
    ELSE
       CALL initial_data('U', u, default=.TRUE.)
       CALL initial_data('V', v, default=.TRUE.)
    END IF

!$OMP PARALLEL DO
    DO k=1, ksize
       DO j=1, jsize
       DO i=0, isize
          u(i,j,k) = u(i,j,k) * imask3d(i,j,k)*imask3d(i+1,j,k)
       END DO
       END DO

       DO j=0, jsize
       DO i=1, isize
          v(i,j,k) = v(i,j,k) * imask3d(i,j,k)*imask3d(i,j+1,k)
       END DO
       END DO
    END DO

    CALL initial_data('W', w, default=perfect_restart, stat=stat_w)
    IF (stat_w) THEN
!$OMP PARALLEL DO
       DO k=0, ksize
          DO j=1, jsize
          DO i=1, isize
             w(i,j,k) = w(i,j,k) * imask3d(i,j,k)*imask3d(i,j,k+1)
          END DO
          END DO
       END DO
    END IF

    IF (open_e) CALL ob_east
    IF (open_w) CALL ob_west
    IF (open_n) CALL ob_north
    IF (open_s) CALL ob_south

    IF (.NOT. stat_w)  CALL update_w(u, v, w, topdown=w_topdown)

    CALL update_velocity_boundary(u, v, w)
    CALL update_velocity_external

    ALLOCATE(p2d(0:isize+1, 0:jsize+1))
    p2d(:,:) = 0.0

    IF (hydrostatic) THEN
       ALLOCATE(q2d(0:isize+1, 0:jsize+1))
       q2d(:,:) = 0.0

       IF (.NOT. bypass_solver) CALL initial_data('P2D', p2d, default=perfect_restart)
    ELSE
       ALLOCATE(p3d(0:isize+1, 0:jsize+1, 0:ksize+1))
       ALLOCATE(q3d(0:isize+1, 0:jsize+1, 0:ksize+1))
       p3d(:,:,:) = 0.0
       q3d(:,:,:) = 0.0

       IF (.NOT. bypass_solver) CALL initial_data('P3D', p3d, default=perfect_restart)
    END IF

    ALLOCATE(gx(0:isize, 1:jsize, 1:ksize))
    ALLOCATE(gy(1:isize, 0:jsize, 1:ksize))
    ALLOCATE(gz(1:isize, 1:jsize, 0:ksize))

    ALLOCATE(gx_ab(0:isize, 1:jsize, 1:ksize))
    ALLOCATE(gy_ab(1:isize, 0:jsize, 1:ksize))
    ALLOCATE(gz_ab(1:isize, 1:jsize, 0:ksize))

!$OMP PARALLEL WORKSHARE
    gx(:,:,:) = 0.0
    gy(:,:,:) = 0.0
    gz(:,:,:) = 0.0

    gx_ab(:,:,:) = 0.0
    gy_ab(:,:,:) = 0.0
    gz_ab(:,:,:) = 0.0
!$OMP END PARALLEL WORKSHARE

    ALLOCATE(dudx( 0:isize+1, 0:jsize+1, 0:ksize+1))
    ALLOCATE(dudy(-1:isize+1,-1:jsize+1, 0:ksize+1))
    ALLOCATE(dudz(-1:isize+1, 0:jsize+1,-1:ksize+1))

    ALLOCATE(dvdx(-1:isize+1,-1:jsize+1, 0:ksize+1))
    ALLOCATE(dvdy( 0:isize+1, 0:jsize+1, 0:ksize+1))
    ALLOCATE(dvdz( 0:isize+1,-1:jsize+1,-1:ksize+1))

    ALLOCATE(dwdx(-1:isize+1, 0:jsize+1,-1:ksize+1))
    ALLOCATE(dwdy( 0:isize+1,-1:jsize+1,-1:ksize+1))
    ALLOCATE(dwdz( 0:isize+1, 0:jsize+1, 0:ksize+1))

    ALLOCATE(ke(0:isize+1, 0:jsize+1, 0:ksize+1))

    ALLOCATE(srate(0:isize+1, 0:jsize+1, 0:ksize+1))

    ALLOCATE(sfreq2(0:isize+1, 0:jsize+1, -1:ksize+1))

    ALLOCATE(taux_sfc(0:isize, 1:jsize))
    ALLOCATE(tauy_sfc(1:isize, 0:jsize))
    ALLOCATE(taux_btm(0:isize, 1:jsize))
    ALLOCATE(tauy_btm(1:isize, 0:jsize))

!$OMP PARALLEL WORKSHARE
    dudx(:,:,:) = 0.0
    dudy(:,:,:) = 0.0
    dudz(:,:,:) = 0.0

    dvdx(:,:,:) = 0.0
    dvdy(:,:,:) = 0.0
    dvdz(:,:,:) = 0.0

    dwdx(:,:,:) = 0.0
    dwdy(:,:,:) = 0.0
    dwdz(:,:,:) = 0.0

    ke(:,:,:)   = 0.0

    srate(:,:,:)  = 0.0

    sfreq2(:,:,:)  = 0.0

    taux_sfc(:,:) = 0.0
    tauy_sfc(:,:) = 0.0
    taux_btm(:,:) = 0.0
    tauy_btm(:,:) = 0.0
!$OMP END PARALLEL WORKSHARE

    CALL update_velocity_diagnostic(no_checkout=.TRUE.)

    CALL nonlinear(gx_ab, gy_ab, gz_ab)

    CALL coriolis(gx_ab, gy_ab, gz_ab)

    CALL initial_data('GX', gx_ab, default=perfect_restart)
    CALL initial_data('GY', gy_ab, default=perfect_restart)
    CALL initial_data('GZ', gz_ab, default=perfect_restart)

  CONTAINS
     SUBROUTINE ob_east
       REAL(8) :: u_ext(1:jsize,1:ksize)
       REAL(8) :: v_ext(0:jsize,1:ksize)
       REAL(8) :: w_ext(1:jsize,0:ksize)
       LOGICAL :: stat
       INTEGER :: j, k

!$OMP PARALLEL WORKSHARE
       u_ext(:,:) = UNDEF
       v_ext(:,:) = UNDEF
       w_ext(:,:) = UNDEF
!$OMP END PARALLEl WORKSHARE

       stat = .FALSE.
       CALL initial_data('U_EAST', u_ext, section='YZ', default=perfect_restart, stat=stat)
       IF (.not. stat) CALL checkin('U_EAST', u_ext, section='YZ', time=t_start)

       stat = .FALSE.
       IF (radi_e .AND. radi_tangential) CALL initial_data('V_EAST', v_ext, section='YZ', default=perfect_restart, stat=stat)
       IF (.not. stat) CALL checkin('V_EAST', v_ext, section='YZ', time=t_start)

       stat = .FALSE.
       IF (radi_e .AND. radi_tangential) CALL initial_data('W_EAST', w_ext, section='YZ', default=perfect_restart, stat=stat)
       IF (.not. stat) CALL checkin('W_EAST', w_ext, section='YZ', time=t_start)

       IF (icoord==ipes-1) THEN
!$OMP PARALLEL
!$OMP DO
          DO k=1, ksize
          DO j=1, jsize
             IF (u_ext(j,k)==UNDEF) u_ext(j,k) = u(isize,j,k)
             u(isize,j,k) = u_ext(j,k)*imask3d(isize,j,k)
          END DO
          END DO

!$OMP DO
          DO k=1, ksize
          DO j=0, jsize
             IF (v_ext(j,k)==UNDEF) v_ext(j,k) = v(isize,j,k)
             v(isize+1,j,k) = v_ext(j,k)*imask3d(isize,j,k)*imask3d(isize,j+1,k)
          END DO
          END DO

!$OMP DO
          DO k=0, ksize
          DO j=1, jsize
             IF (w_ext(j,k)==UNDEF) w_ext(j,k) = w(isize,j,k)
             w(isize+1,j,k) = w_ext(j,k)*imask3d(isize,j,k)*imask3d(isize,j,k+1)
          END DO
          END DO
!$OMP END PARALLEL
       END IF
     END SUBROUTINE ob_east

     SUBROUTINE ob_west
       REAL(8) :: u_ext(1:jsize,1:ksize)
       REAL(8) :: v_ext(0:jsize,1:ksize)
       REAL(8) :: w_ext(1:jsize,0:ksize)
       LOGICAL :: stat
       INTEGER :: j, k

!$OMP PARALLEL WORKSHARE
       u_ext(:,:) = UNDEF
       v_ext(:,:) = UNDEF
       w_ext(:,:) = UNDEF
!$OMP END PARALLEl WORKSHARE

       stat = .FALSE.
       CALL initial_data('U_WEST', u_ext, section='YZ', default=perfect_restart, stat=stat)
       IF (.not. stat) CALL checkin('U_WEST', u_ext, section='YZ', time=t_start)

       stat = .FALSE.
       IF (radi_w .AND. radi_tangential) CALL initial_data('V_WEST', v_ext, section='YZ', default=perfect_restart, stat=stat)
       IF (.not. stat) CALL checkin('V_WEST', v_ext, section='YZ', time=t_start)

       stat = .FALSE.
       IF (radi_w .AND. radi_tangential) CALL initial_data('W_WEST', w_ext, section='YZ', default=perfect_restart, stat=stat)
       IF (.not. stat) CALL checkin('W_WEST', w_ext, section='YZ', time=t_start)

       IF (icoord==0) THEN
!$OMP PARALLEL
!$OMP DO
          DO k=1, ksize
          DO j=1, jsize
             IF (u_ext(j,k)==UNDEF) u_ext(j,k) = u(1,j,k)
             u(0,j,k) = u_ext(j,k)*imask3d(1,j,k)
          END DO
          END DO

!$OMP DO
          DO k=1, ksize
          DO j=0, jsize
             IF (v_ext(j,k)==UNDEF) v_ext(j,k) = v(1,j,k)
             v(0,j,k) = v_ext(j,k)*imask3d(1,j,k)*imask3d(1,j+1,k)
          END DO
          END DO

!$OMP DO
          DO k=0, ksize
          DO j=1, jsize
             IF (w_ext(j,k)==UNDEF) w_ext(j,k) = w(1,j,k)
             w(0,j,k) = w_ext(j,k)*imask3d(1,j,k)*imask3d(1,j,k+1)
          END DO
          END DO
!$OMP END PARALLEL
       END IF
     END SUBROUTINE ob_west

     SUBROUTINE ob_north
       REAL(8) :: u_ext(0:isize,1:ksize)
       REAL(8) :: v_ext(1:isize,1:ksize)
       REAL(8) :: w_ext(1:isize,0:ksize)
       LOGICAL :: stat
       INTEGER :: i, k

!$OMP PARALLEL WORKSHARE
       u_ext(:,:) = UNDEF
       v_ext(:,:) = UNDEF
       w_ext(:,:) = UNDEF
!$OMP END PARALLEl WORKSHARE

       stat = .FALSE.
       CALL initial_data('V_NORTH', v_ext, section='XZ', default=perfect_restart, stat=stat)
       IF (.not. stat) CALL checkin('V_NORTH', v_ext, section='XZ', time=t_start)

       stat = .FALSE.
       IF (radi_n .AND. radi_tangential) CALL initial_data('U_NORTH', u_ext, section='XZ', default=perfect_restart, stat=stat)
       IF (.not. stat) CALL checkin('V_NORTH', u_ext, section='XZ', time=t_start)

       stat = .FALSE.
       IF (radi_n .AND. radi_tangential) CALL initial_data('W_NORTH', w_ext, section='XZ', default=perfect_restart, stat=stat)
       IF (.not. stat) CALL checkin('W_NORTH', w_ext, section='XZ', time=t_start)

       IF (jcoord==jpes-1) THEN
!$OMP PARALLEL
!$OMP DO
          DO k=1, ksize
          DO i=1, isize
             IF (v_ext(i,k)==UNDEF) v_ext(i,k) = v(i,jsize,k)
             v(i,jsize,k) = v_ext(i,k)*imask3d(i,jsize,k)
          END DO
          END DO

!$OMP DO
          DO k=1, ksize
          DO i=0, isize
             IF (u_ext(i,k)==UNDEF) u_ext(i,k) = u(i,jsize,k)
             u(i,jsize+1,k) = u_ext(i,k)*imask3d(i,jsize,k)*imask3d(i+1,jsize,k)
          END DO
          END DO

!$OMP DO
          DO k=0, ksize
          DO i=1, isize
             IF (w_ext(i,k)==UNDEF) w_ext(i,k) = w(i,jsize,k)
             w(i,jsize+1,k) = w_ext(i,k)*imask3d(i,jsize,k)*imask3d(i,jsize,k+1)
          END DO
          END DO
!$OMP END PARALLEL
       END IF
     END SUBROUTINE ob_north

     SUBROUTINE ob_south
       REAL(8) :: u_ext(0:isize,1:ksize)
       REAL(8) :: v_ext(1:isize,1:ksize)
       REAL(8) :: w_ext(1:isize,0:ksize)
       LOGICAL :: stat
       INTEGER :: i, k

!$OMP PARALLEL WORKSHARE
       u_ext(:,:) = UNDEF
       v_ext(:,:) = UNDEF
       w_ext(:,:) = UNDEF
!$OMP END PARALLEl WORKSHARE

       stat = .FALSE.
       CALL initial_data('V_SOUTH', v_ext, section='XZ', default=perfect_restart, stat=stat)
       IF (.not. stat) CALL checkin('V_SOUTH', v_ext, section='XZ', time=t_start)

       stat = .FALSE.
       IF (radi_s .AND. radi_tangential) CALL initial_data('U_SOUTH', u_ext, section='XZ', default=perfect_restart, stat=stat)
       IF (.not. stat) CALL checkin('V_SOUTH', u_ext, section='XZ', time=t_start)

       stat = .FALSE.
       IF (radi_s .AND. radi_tangential) CALL initial_data('W_SOUTH', w_ext, section='XZ', default=perfect_restart, stat=stat)
       IF (.not. stat) CALL checkin('W_SOUTH', w_ext, section='XZ', time=t_start)

       IF (jcoord==0) THEN
!$OMP PARALLEL
!$OMP DO
          DO k=1, ksize
          DO i=1, isize
             IF (v_ext(i,k)==UNDEF) v_ext(i,k) = v(i,1,k)
             v(i,0,k) = v_ext(i,k)*imask3d(i,1,k)
          END DO
          END DO

!$OMP DO
          DO k=1, ksize
          DO i=0, isize
             IF (u_ext(i,k)==UNDEF) u_ext(i,k) = u(i,1,k)
             u(i,0,k) = u_ext(i,k)*imask3d(i,1,k)*imask3d(i+1,1,k)
          END DO
          END DO

!$OMP DO
          DO k=0, ksize
          DO i=1, isize
             IF (w_ext(i,k)==UNDEF) w_ext(i,k) = w(i,1,k)
             w(i,0,k) = w_ext(i,k)*imask3d(i,1,k)*imask3d(i,1,k+1)
          END DO
          END DO
!$OMP END PARALLEL
       END IF
     END SUBROUTINE ob_south


     SUBROUTINE initial_stream
       REAL(8) :: tmp3d(0:isize,0:jsize,1:ksize)
       tmp3d(:,:,:) = 0.0
       CALL initial_data('PSI', tmp3d)
       CALL stream_uv(tmp3d, u, v)
       CALL update_velocity_boundary(u, v)
     END SUBROUTINE initial_stream

  END SUBROUTINE init_velocity

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE finalize_velocity
    DEALLOCATE(u, v, w)

    DEALLOCATE(gx, gx_ab)
    DEALLOCATE(gy, gy_ab)
    DEALLOCATE(gz, gz_ab)

    IF (ALLOCATED(p2d)) DEALLOCATE(p2d)
    IF (ALLOCATED(p3d)) DEALLOCATE(p3d)
    IF (ALLOCATED(q2d)) DEALLOCATE(q2d)
    IF (ALLOCATED(q3d)) DEALLOCATE(q3d)

    DEALLOCATE(dudx, dudy, dudz)
    DEALLOCATE(dvdx, dvdy, dvdz)
    DEALLOCATE(dwdx, dwdy, dwdz)
    DEALLOCATE(ke, srate)

  END SUBROUTINE finalize_velocity

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE openboundary_velocity
    IF (open_e) CALL ob_east
    IF (open_w) CALL ob_west
    IF (open_n) CALL ob_north
    IF (open_s) CALL ob_south

!!    CALL tidal

    IF (input_registered('W_BOTTOM')) CALL ob_bottom

    CALL update_velocity_external

  CONTAINS
     SUBROUTINE ob_east
       REAL(8) :: u_ext(1:jsize,1:ksize)
       REAL(8) :: v_ext(0:jsize,1:ksize)
       REAL(8) :: w_ext(1:jsize,0:ksize)
!!!!!!!!YJKIM
       REAL(8) :: u_barocli(1:jsize,1:ksize)
       REAL(8) :: u_barotro(1:jsize) 
       REAL(8) :: u_barotro_e(1:jsize)       
       REAL(8) :: ext_ssh_e(1:jsize)
       REAL(8) :: ext_tide_e(1:jsize)
       LOGICAL :: stat1
       LOGICAL :: stat2
!!!!!!!!!!!!!

       REAL(8) :: du(1:jsize)
       REAL(8) :: a, bx, by, cx, cy, d
       REAL(8) :: relax(1:10)
       INTEGER :: j, k

!$OMP PARALLEL WORKSHARE
       u_ext(:,:) = UNDEF
       v_ext(:,:) = 0.0
       w_ext(:,:) = 0.0
!!!!!!!!!YJKIM     
       u_barotro_e(:) = 0.0
!!!!!!!!!!!!!       
!$OMP END PARALLEL WORKSHARE

       CALL checkin('U_EAST', u_ext, time=t_current+dtime, section='YZ')
       CALL checkin('V_EAST', v_ext, time=t_current+dtime, section='YZ')
       CALL checkin('W_EAST', w_ext, time=t_current+dtime, section='YZ')
!!!!!!!!!YJKIM       
       CALL checkin('SSH_EAST', ext_ssh_e, time=t_current+dtime, axis='Y',stat=stat1)
       CALL checkin('UTIDE_EAST', ext_tide_e, time=t_current+dtime, axis='Y',stat=stat2)
!!!!!!!!!!!!!
       IF (icoord/=ipes-1) RETURN
!!!!!!!!!!!!! YJKIM 1. barotro Flather / barocli relaxation (Carter and Merrifield, 2007)
!!!!!!!!!!!!!       2. barotro Flather / barocli nogradient <- current version       
!!!!!!!!!!!!! to turn off Tidal forcing, put SSH / Vel to zero by boundary condition --> passive Flather        
       IF (stat1 .AND. stat2 .AND. icoord==ipes-1) THEN
            DO j=1, jsize
               IF (lmask2d(isize,j) .AND. ext_tide_e(j) /= UNDEF .AND. ext_ssh_e(j) /= UNDEF) THEN
                   u_barotro_e(j) = ext_tide_e(j) + (cext(isize,j)/(wct_ref(isize,j)+ssh(isize,j)))*(ssh(isize,j)-ext_ssh_e(j))   
               ELSE
                   u_barotro_e(j) = 0.0
               END IF
            END DO

         u_barocli(:,:)=0
         u_barotro(:)=0
          DO k=1, ksize
          DO j=1, jsize
                u_barotro(j) = u_barotro(j) + (u(isize-2,j,k)*dsx(isize-1,j,k))*imask3d(isize-1,j,k)
          END DO
          END DO

          CALL vsum(u_barotro, all=.TRUE.)

          DO j=1, jsize
             IF (dsx2d(isize-1,j) > 0.0) THEN
                u_barotro(j) = u_barotro(j) / dsx2d(isize-1,j)
             ELSE
                u_barotro(j) = 0.0
             END IF
          END DO

         DO k=1, ksize
         DO j=1, jsize
            u_barocli(j,k)=u(isize-2,j,k)-u_barotro(j)*imask3d(isize-1,j,k)
         END DO
         END DO

         DO k=1, ksize
         DO j=1, jsize
            u(isize-1,j,k)=u_barocli(j,k)*imask3d(isize-1,j,k)+u_barotro_e(j)*imask3d(isize-1,j,k)
            u(isize,j,k) = u(isize-1,j,k)
            v(isize,j,k) = v(isize-1,j,k)
            w(isize,j,k) = w(isize-1,j,k)
         END DO
         END DO   



!             DO j=1,10
!                a=(j-1)/2
!                relax(j)=DTANH(a)
!                ! relax(j)=0.1*(j-1)
!                ! relax(j)=((10-j+1)/10)**2
!             END DO
!             DO k=1, ksize
!             DO j=1, jsize
!              u(isize-9,j,k)= relax(10)*u(isize-9,j,k)+ (1-relax(10))*u_barotro_e(j)*imask3d(isize-9,j,k) 
!              u(isize-8,j,k) = relax(9)*u(isize-8,j,k) + (1-relax(9))*u_barotro_e(j)*imask3d(isize-8,j,k) 
!              u(isize-7,j,k) = relax(8)*u(isize-7,j,k) + (1-relax(8))*u_barotro_e(j)*imask3d(isize-7,j,k) 
!              u(isize-6,j,k) = relax(7)*u(isize-6,j,k) + (1-relax(7))*u_barotro_e(j)*imask3d(isize-6,j,k) 
!              u(isize-5,j,k) = relax(6)*u(isize-5,j,k) + (1-relax(6))*u_barotro_e(j)*imask3d(isize-5,j,k) 
!              u(isize-4,j,k) = relax(5)*u(isize-4,j,k) + (1-relax(5))*u_barotro_e(j)*imask3d(isize-4,j,k) 
!              u(isize-3,j,k) = relax(4)*u(isize-3,j,k) + (1-relax(4))*u_barotro_e(j)*imask3d(isize-3,j,k) 
!              u(isize-2,j,k) = relax(3)*u(isize-2,j,k) + (1-relax(3))*u_barotro_e(j)*imask3d(isize-2,j,k) 
!              u(isize-1,j,k) = relax(2)*u(isize-1,j,k) + (1-relax(2))*u_barotro_e(j)*imask3d(isize-1,j,k)
              
!              u(isize,j,k)=  u_barotro_e(j)*imask3d(isize,j,k)
! !             u(isize,j,k) = u(isize-1,j,k)
!              v(isize,j,k) = v(isize-1,j,k)

!              w(isize-10,j,k)=relax(10)*w(isize-10,j,k)
!              w(isize-9,j,k)=relax(9)*w(isize-9,j,k)
!              w(isize-8,j,k)=relax(8)*w(isize-8,j,k)
!              w(isize-7,j,k)=relax(7)*w(isize-7,j,k)
!              w(isize-6,j,k)=relax(6)*w(isize-6,j,k)
!              w(isize-5,j,k)=relax(5)*w(isize-5,j,k)
!              w(isize-4,j,k)=relax(4)*w(isize-4,j,k)
!              w(isize-3,j,k)=relax(3)*w(isize-3,j,k)
!              w(isize-2,j,k)=relax(2)*w(isize-2,j,k) 
!              w(isize-1,j,k)=0
!              w(isize,j,k) = w(isize-1,j,k)
!             ENDDO
!             ENDDO    
      END IF   


! !$OMP PARALLEL DO PRIVATE(a, bx, by, cx, cy, d)
!           DO k=1, ksize
!           DO j=1, jsize
!             IF (.NOT. lmask3d(isize,j,k)) CYCLE
!             ! a = u(isize-1,j,k) - u_old(isize-1,j,k)

!             ! bx = u_old(isize-2,j,k) - u_old(isize-1,j,k)

!             ! IF (radi_oblique) THEN
!             !    by = 0.5*(u_old(isize-1,j-1,k) - u_old(isize-1,j+1,k))
!             ! ELSE
!             !    by = 0.0
!             ! END IF

!             ! IF (a*bx > 0.0) THEN
!             !    cx = a*bx/(bx**2 + by**2)
!             !    cy = a*by/(bx**2 + by**2)
!             !    d = radi_nudge_e(1)
!             ! ELSE
!             !    cx = 0.0
!             !    cy = 0.0
!             !    d = radi_nudge_e(2)
!             ! END IF

!             ! u(isize,j,k) = u_old(isize,j,k) + cx*u(isize-1,j,k)

!             ! IF (cy > 0.0) THEN
!             !    u(isize,j,k) = u(isize,j,k) + cy*u_old(isize,j-1,k)
!             ! ELSE
!             !    u(isize,j,k) = u(isize,j,k) - cy*u_old(isize,j+1,k)
!             ! END IF

!             ! u(isize,j,k) = u(isize,j,k) / (1.0 + cx + abs(cy))

! !            IF (u_ext(j,k) /= UNDEF) u(isize,j,k) = (1.0-d)*u(isize,j,k) + d*u_ext(j,k)
!              u(isize,j,k) = u(isize-1,j,k)
!              v(isize,j,k) = v(isize-1,j,k)              
!              w(isize,j,k) = w(isize-1,j,k)
!          END DO
!          END DO
!        ELSE
!           du(:) = 0.0
!           DO k=1, ksize
!              DO j=1, jsize
!                 IF (u_ext(j,k) == UNDEF) u_ext(j,k) = 0.0
!                 du(j) = du(j) + (u_ext(j,k)*dsx_ref(isize,j,k) - (u(isize-1,j,k)*open_alpha_e + u_ext(j,k)*(1.0-open_alpha_e))*dsx(isize,j,k))*imask3d(isize,j,k)
!              END DO
!           END DO

!           CALL vsum(du, all=.TRUE.)

!           DO j=1, jsize
!              IF (dsx2d(isize,j) > 0.0) THEN
!                 du(j) = du(j) / dsx2d(isize,j)
!              ELSE
!                 du(j) = 0.0
!              END IF
!           END DO

! !$OMP PARALLEL DO
!           DO k=1, ksize
!           DO j=1, jsize
!              u(isize,j,k) = (u(isize-1,j,k)*open_alpha_e + u_ext(j,k)*(1.0-open_alpha_e) + du(j))*imask3d(isize,j,k)
!           END DO
!           END DO
!        END IF

!        IF (radi_e .AND. radi_tangential) THEN
! !$OMP PARALLEL DO PRIVATE(a, bx, by, cx, cy, d)
!           DO k=1, ksize
!           DO j=0, jsize
!              IF (.NOT. (lmask3d(isize,j,k) .AND. lmask3d(isize,j+1,k))) CYCLE

!              a = v(isize,  j,k) - v_old(isize,j,k)

!              bx = v_old(isize-1,j,k) - v_old(isize,j,k)

!              IF (radi_oblique) THEN
!                 by = 0.5*(v_old(isize,j-1,k) - v_old(isize,j+1,k))
!              ELSE
!                 by = 0.0
!              END IF

!              IF (a*bx > 0.0) THEN
!                 cx = a*bx / (bx**2 + by**2)
!                 cy = a*by / (bx**2 + by**2)
!                 d = radi_nudge_e(1)
!              ELSE
!                 cx = 0.0
!                 cy = 0.0
!                 d = radi_nudge_e(2)
!              END IF

!              v(isize+1,j,k) = v_old(isize+1,j,k) + cx*v(isize,j,k)

!              IF (cy > 0) THEN
!                 v(isize+1,j,k) = v(isize+1,j,k) + cy*v_old(isize+1,j-1,k)
!              ELSE
!                 v(isize+1,j,k) = v(isize+1,j,k) - cy*v_old(isize+1,j+1,k)
!              END IF

!              v(isize+1,j,k) = v(isize+1,j,k) / (1.0 + cx + abs(cy))

! !             IF (v_ext(j,k) /= UNDEF) v(isize+1,j,k) = (1.0-d)*v(isize+1,j,k) + d*v_ext(j,k)
!           END DO
!           END DO

! !$OMP PARALLEL DO PRIVATE(a, bx, by, cx, cy, d)
!           DO k=0, ksize
!           DO j=1, jsize
!              IF (.NOT. (lmask3d(isize,j,k) .AND. lmask3d(isize,j,k+1))) CYCLE

!              a = w(isize,j,k) - w_old(isize,j,k)

!              bx = w_old(isize-1,j,k) - w_old(isize,j,k)

!              IF (radi_oblique) THEN
!                 by = 0.5*(w_old(isize,j-1,k) - w_old(isize,j+1,k))
!              ELSE
!                 by = 0.0
!              END IF

!              IF (a*bx > 0.0) THEN
!                 cx = a*bx / (bx**2 + by**2)
!                 cy = a*by / (bx**2 + by**2)
!                 d = radi_nudge_e(1)
!              ELSE
!                 cx = 0.0
!                 cy = 0.0
!                 d = radi_nudge_e(2)
!              END IF

!              w(isize+1,j,k) = w_old(isize+1,j,k) + cx*w(isize,j,k)

!              IF (cy > 0.0) THEN
!                 w(isize+1,j,k) = w(isize+1,j,k) + cy*w_old(isize+1,j-1,k)
!              ELSE
!                 w(isize+1,j,k) = w(isize+1,j,k) - cy*w_old(isize+1,j+1,k)
!              END IF

!              w(isize+1,j,k) = w(isize+1,j,k) / (1.0 + cx + abs(cy))

! !             IF (w_ext(j,k) /= UNDEF) w(isize+1,j,k) = (1.0-d)*w(isize+1,j,k) + d*w_ext(j,k)
!           END DO
!           END DO
!        ELSE
! !$OMP PARALLEL DO
!           DO k=1, ksize
!           DO j=0, jsize
!              IF (v_ext(j,k) == UNDEF) THEN
!                 v(isize+1,j,k) = v(isize,j,k)
!              ELSE
!                 v(isize+1,j,k) = v_ext(j,k)*imask3d(isize,j,k)*imask3d(isize,j+1,k)
!              END IF
!           END DO
!           END DO

! !$OMP PARALLEL DO
!           DO k=0, ksize
!           DO j=1, jsize
!              IF (w_ext(j,k) == UNDEF) THEN
!                 w(isize+1,j,k) = w(isize,j,k)
!              ELSE
!                 w(isize+1,j,k) = w_ext(j,k)*imask3d(isize,j,k)*imask3d(isize,j,k+1)
!              END IF
!           END DO
!           END DO
!        END IF

     END SUBROUTINE ob_east

     SUBROUTINE ob_west
       REAL(8) :: u_ext(1:jsize,1:ksize)
       REAL(8) :: v_ext(0:jsize,1:ksize)
       REAL(8) :: w_ext(1:jsize,0:ksize)
       REAL(8) :: du(1:jsize)
       REAL(8) :: a, bx, by, cx, cy, d
       INTEGER :: j, k

!$OMP PARALLEL WORKSHARE
       u_ext(:,:) = UNDEF
       v_ext(:,:) = 0.0
       w_ext(:,:) = 0.0
!$OMP END PARALLEL WORKSHARE

       CALL checkin('U_WEST', u_ext, time=t_current+dtime, section='YZ')
       CALL checkin('V_WEST', v_ext, time=t_current+dtime, section='YZ')
       CALL checkin('W_WEST', w_ext, time=t_current+dtime, section='YZ')

       IF (icoord/=0) RETURN

       IF (radi_w) THEN
!$OMP PARALLEL DO PRIVATE(a, bx, by, cx, cy, d)
          DO k=1, ksize
          DO j=1, jsize
             IF (.NOT. lmask3d(1,j,k)) CYCLE

             a = u(1,j,k) - u_old(1,j,k)

             bx = u_old(2,j,k) - u_old(1,j,k)

             IF (radi_oblique) THEN
                by = 0.5*(u_old(1,j-1,k) - u_old(1,j+1,k))
             ELSE
                by = 0.0
             END IF

             IF (a*bx > 0.0) THEN
                cx = a*bx/(bx**2 + by**2)
                cy = a*by/(bx**2 + by**2)
                d = radi_nudge_w(1)
             ELSE
                cx = 0.0
                cy = 0.0
                d = radi_nudge_w(2)
             END IF

             u(0,j,k) = u_old(0,j,k) + cx*u(1,j,k)
             IF (cy > 0.0) THEN
                u(0,j,k) = u(0,j,k) + cy*u_old(0,j-1,k)
             ELSE
                u(0,j,k) = u(0,j,k) - cy*u_old(0,j+1,k)
             END IF
             u(0,j,k) = u(0,j,k) / (1.0 + cx + abs(cy))

!             IF (u_ext(j,k) /= UNDEF) u(0,j,k) = (1.0-d)*u(0,j,k) + d*u_ext(j,k)

          END DO
          END DO
       ELSE
          du(:) = 0.0
          DO k=1, ksize
             DO j=1, jsize
                IF (u_ext(j,k) == UNDEF) u_ext(j,k) = 0.0
                du(j) = du(j) + (u_ext(j,k)*dsx_ref(0,j,k) - (u(1,j,k)*open_alpha_w + u_ext(j,k)*(1.0-open_alpha_w))*dsx(0,j,k))*imask3d(1,j,k)
             END DO
          END DO

          CALL vsum(du, all=.TRUE.)

          DO j=1, jsize
             IF (dsx2d(0,j) > 0.0) THEN
                du(j) = du(j) / dsx2d(0,j)
             ELSE
                du(j) = 0.0
             END IF
          END DO

!$OMP PARALLEL DO
          DO k=1, ksize
          DO j=1, jsize
             u(0,j,k) = (u(0,j,k)*open_alpha_e + u_ext(j,k)*(1.0-open_alpha_e) + du(j))*imask3d(1,j,k)
          END DO
          END DO
       END IF

       IF (radi_w .AND. radi_tangential) THEN
!$OMP PARALLEL DO PRIVATE(a, bx, by, cx, cy, d)
          DO k=1, ksize
          DO j=0, jsize
             IF (.NOT. (lmask3d(1,j,k) .AND. lmask3d(1,j+1,k))) CYCLE

             a = v(1,j,k) - v_old(1,j,k)

             bx = v_old(2,j,k) - v_old(1,j,k)

             IF (radi_oblique) THEN
                by = 0.5*(v_old(1,j-1,k) - v_old(1,j+1,k))
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

             v(0,j,k) = v_old(0,j,k) + cx*v(1,j,k)

             IF (cy > 0.0) THEN
                v(0,j,k) = v(0,j,k) + cy*v_old(0,j-1,k)
             ELSE
                v(0,j,k) = v(0,j,k) - cy*v_old(0,j+1,k)
             END IF

             v(0,j,k) = v(0,j,k) / (1.0 + cx + abs(cy))

!             IF (v_ext(j,k) /= UNDEF) v(0,j,k) = (1.0-d)*v(0,j,k) + d*v_ext(j,k)
          END DO
          END DO

!$OMP PARALLEL DO PRIVATE(a, bx, by, cx, cy, d)
          DO k=0, ksize
          DO j=1, jsize
             IF (.NOT. (lmask3d(1,j,k) .AND. lmask3d(1,j,k+1))) CYCLE

             a = w(1,j,k) - w_old(1,j,k)

             bx = w_old(2,j,k) - w_old(1,j,k)

             IF (radi_oblique) THEN
                by = 0.5*(w_old(1,j-1,k) - w_old(1,j+1,k))
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

             w(0,j,k) = w_old(0,j,k) + cx*w(1,j,k)

             IF (cy > 0.0) THEN
                w(0,j,k) = w(0,j,k) + cy*w_old(0,j-1,k)
             ELSE
                w(0,j,k) = w(0,j,k) - cy*w_old(0,j+1,k)
             END IF

             w(0,j,k) = w(0,j,k) / (1.0 + cx + abs(cy))

!             IF (w_ext(j,k) /= UNDEF) w(0,j,k) = (1.0-d)*w(0,j,k) + d*w_ext(j,k)
          END DO
          END DO
       ELSE
!$OMP PARALLEL DO
          DO k=1, ksize
          DO j=0, jsize
             IF (v_ext(j,k) == UNDEF) THEN
                v(0,j,k) = v(1,j,k)
             ELSE
                v(0,j,k) = v_ext(j,k)*imask3d(1,j,k)*imask3d(1,j+1,k)
             END IF
          END DO
          END DO

!$OMP PARALLEL DO
          DO k=0, ksize
          DO j=1, jsize
             IF (w_ext(j,k) == UNDEF) THEN
                w(0,j,k) = w(1,j,k)
             ELSE
                w(0,j,k) = w_ext(j,k)*imask3d(1,j,k)*imask3d(1,j,k+1)
             END IF
          END DO
          END DO
       END IF

     END SUBROUTINE ob_west

     SUBROUTINE ob_north
       REAL(8) :: u_ext(0:isize,1:ksize)
       REAL(8) :: v_ext(1:isize,1:ksize)
       REAL(8) :: w_ext(1:isize,0:ksize)
!!!!!!!!YJKIM
       REAL(8) :: v_barocli(1:isize,1:ksize)
       REAL(8) :: v_barotro(1:isize)
       REAL(8) :: v_barotro_n(1:isize)       
       REAL(8) :: ext_ssh_n(1:isize)
       REAL(8) :: ext_tide_n(1:isize)
       LOGICAL :: stat1
       LOGICAL :: stat2
!!!!!!!!!!!!!
       REAL(8) :: dv(1:isize)
       REAL(8) :: relax(1:10)
       REAL(8) :: a, bx, by, cx, cy, d
       INTEGER :: i, k

!$OMP PARALLEL WORKSHARE
       u_ext(:,:) = 0.0
       v_ext(:,:) = UNDEF
       w_ext(:,:) = 0.0
!!!!!!!!!YJKIM       
       v_barotro_n(:) = 0.0
!!!!!!!!!!!!!              
!$OMP END PARALLEL WORKSHARE

       CALL checkin('U_NORTH', u_ext, time=t_current+dtime, section='XZ')
       CALL checkin('V_NORTH', v_ext, time=t_current+dtime, section='XZ')
       CALL checkin('W_NORTH', w_ext, time=t_current+dtime, section='XZ')
!!!!!!!!!YJKIM       
       CALL checkin('SSH_NORTH', ext_ssh_n, time=t_current+dtime, axis='X', stat=stat1)
       CALL checkin('VTIDE_NORTH', ext_tide_n, time=t_current+dtime, axis='X', stat=stat2)
!!!!!!!!!!!!!
       IF (jcoord/=jpes-1) RETURN

!!!!!!!!!!!!! YJKIM 1. barotro Flather / barocli relaxation (Carter and Merrifield, 2007)
!!!!!!!!!!!!!       2. barotro Flather / barocli nogradient <- current version        
!!!!!!!!!!!!! to turn off Tidal forcing, put SSH / Vel to zero by boundary condition --> passive Flather        
       IF (stat1 .AND. stat2 .AND. jcoord==jpes-1) THEN
         DO i=1, isize
            IF (lmask2d(i,jsize) .AND. ext_tide_n(i) /= UNDEF .AND. ext_ssh_n(i) /= UNDEF) THEN
                v_barotro_n(i) = ext_tide_n(i) + (cext(i,jsize)/(wct_ref(i,jsize)+ssh(i,jsize)))*(ssh(i,jsize)-ext_ssh_n(i)) 
            ELSE
                v_barotro_n(i) = 0.0
            END IF
         END DO

         v_barocli(:,:)=0
         v_barotro(:)=0

          DO k=1, ksize
          DO i=1, isize
               v_barotro(i) =v_barotro(i) + (v(i,jsize-2,k)*dsy(i,jsize-1,k))*imask3d(i,jsize-1,k)
          END DO
          END DO

          CALL vsum(v_barotro, all=.TRUE.)

          DO i=1, isize
             IF (dsy2d(i,jsize-1) > 0.0) THEN
               v_barotro(i) =v_barotro(i) / dsy2d(i,jsize-1)
             ELSE
               v_barotro(i) = 0.0
             END IF
          END DO

         DO k=1, ksize
         DO i=1, isize
            v_barocli(i,k)=v(i,jsize-2,k)-v_barotro(i)*imask3d(i,jsize-1,k)
         END DO
         END DO

         DO k=1, ksize
         DO i=1, isize
            u(i,jsize,k) = u(i,jsize-1,k)
            v(i,jsize-1,k)=v_barocli(i,k)*imask3d(i,jsize-1,k)+v_barotro_n(i)*imask3d(i,jsize-1,k)
            v(i,jsize,k) = v(i,jsize-1,k)
            w(i,jsize,k) = w(i,jsize-1,k)
         END DO
         END DO   
         ! DO i=1,10
         !   a=(i-1)/2
         !   relax(i)=DTANH(a)
         !    ! relax(i)=0.1*(i-1)
         !    ! relax(i)=((10-i+1)/10)**2
         ! END DO
         ! DO k=1, ksize
         ! DO i=1, isize
         !  v(i,jsize-9,k)= relax(10)*v(i,jsize-9,k)+ (1-relax(10))*v_barotro_n(i)*imask3d(i,jsize-9,k) 
         !  v(i,jsize-8,k) = relax(9)*v(i,jsize-8,k) + (1-relax(9))*v_barotro_n(i)*imask3d(i,jsize-8,k) 
         !  v(i,jsize-7,k) = relax(8)*v(i,jsize-7,k) + (1-relax(8))*v_barotro_n(i)*imask3d(i,jsize-7,k) 
         !  v(i,jsize-6,k) = relax(7)*v(i,jsize-6,k) + (1-relax(7))*v_barotro_n(i)*imask3d(i,jsize-6,k) 
         !  v(i,jsize-5,k) = relax(6)*v(i,jsize-5,k) + (1-relax(6))*v_barotro_n(i)*imask3d(i,jsize-5,k) 
         !  v(i,jsize-4,k) = relax(5)*v(i,jsize-4,k) + (1-relax(5))*v_barotro_n(i)*imask3d(i,jsize-4,k) 
         !  v(i,jsize-3,k) = relax(4)*v(i,jsize-3,k) + (1-relax(4))*v_barotro_n(i)*imask3d(i,jsize-3,k) 
         !  v(i,jsize-2,k) = relax(3)*v(i,jsize-2,k) + (1-relax(3))*v_barotro_n(i)*imask3d(i,jsize-2,k) 
         !  v(i,jsize-1,k) = relax(2)*v(i,jsize-1,k) + (1-relax(2))*v_barotro_n(i)*imask3d(i,jsize-1,k) 
         !  v(i,jsize,k)=  v_barotro_n(i)*imask3d(i,jsize,k)
         ! !  v(i,jsize,k) = v(i,jsize-1,k)
         !  u(i,jsize,k) = u(i,jsize-1,k)
         
         !  w(i,jsize-10,k)= relax(10)*w(i,jsize-10,k)
         !  w(i,jsize-9,k)= relax(9)*w(i,jsize-9,k)
         !  w(i,jsize-8,k)= relax(8)*w(i,jsize-8,k)
         !  w(i,jsize-7,k)= relax(7)*w(i,jsize-7,k)
         !  w(i,jsize-6,k)= relax(6)*w(i,jsize-6,k)
         !  w(i,jsize-5,k)= relax(5)*w(i,jsize-5,k)
         !  w(i,jsize-4,k)= relax(4)*w(i,jsize-4,k)
         !  w(i,jsize-3,k)= relax(3)*w(i,jsize-3,k)
         !  w(i,jsize-2,k)= relax(2)*w(i,jsize-2,k)
         !  w(i,jsize-1,k)=0
         !  w(i,jsize,k) = w(i,jsize-1,k)
         ! ENDDO
         ! ENDDO

       END IF          
  
! !$OMP PARALLEL DO PRIVATE(a, bx, by, cx, cy, d)
!           DO k=1, ksize
!           DO i=1, isize
!              IF (.NOT. lmask3d(i,jsize,k)) CYCLE

!             !  a = v(i,jsize-1,k) - v_old(i,jsize-1,k)

!             !  by = v_old(i,jsize-2,k) - v_old(i,jsize-1,k)

!             !  IF (radi_oblique) THEN
!             !     bx = 0.5*(v_old(i-1,jsize-1,k) - v_old(i+1,jsize-1,k))
!             !  ELSE
!             !     bx = 0.0
!             !  END IF

!             !  IF (a*by > 0.0) THEN
!             !     cx = a*bx/(bx**2 + by**2)
!             !     cy = a*by/(bx**2 + by**2)
!             !     d = radi_nudge_n(1)
!             !  ELSE
!             !     cx = 0.0
!             !     cy = 0.0
!             !     d = radi_nudge_n(2)
!             !  END IF

!             !  v(i,jsize,k) = v_old(i,jsize,k) + cy*v(i,jsize-1,k)

!             !  IF (cx > 0.0) THEN
!             !     v(i,jsize,k) = v(i,jsize,k) + cx*v_old(i-1,jsize,k)
!             !  ELSE
!             !     v(i,jsize,k) = v(i,jsize,k) - cx*v_old(i+1,jsize,k)
!             !  END IF

!             !  v(i,jsize,k) = v(i,jsize,k) / (1.0 + cy + abs(cx))

! !             IF (v_ext(i,k) /= UNDEF) v(i,jsize,k) = (1.0-d)*v(i,jsize,k) + d*v_ext(i,k)
!               v(i,jsize,k) = v(i,jsize-1,k)
              
!               u(i,jsize,k) = u(i,jsize-1,k)
!               w(i,jsize,k) = w(i,jsize-1,k)
!           END DO
!           END DO
!        ELSE
!           dv(:) = 0.0
!           DO k=1, ksize
!              DO i=1, isize
!                 IF (v_ext(i,k) == UNDEF) v_ext(i,k) = 0.0
!                 dv(i) = dv(i) + (v_ext(i,k)*dsy_ref(i,jsize,k) - (v(i,jsize-1,k)*open_alpha_n + v_ext(i,k)*(1.0-open_alpha_n))*dsy(i,jsize,k))*imask3d(i,jsize,k)
!              END DO
!           END DO

!           CALL vsum(dv, all=.TRUE.)

!           DO i=1, isize
!              IF (dsy2d(i,jsize) > 0.0) THEN
!                 dv(i) = dv(i) / dsy2d(i,jsize)
!              ELSE
!                 dv(i) = 0.0
!              END IF
!           END DO

! !$OMP PARALLEL DO
!           DO k=1, ksize
!           DO i=1, isize
!              v(i,jsize,  k) = (v(i,jsize-1,k)*open_alpha_n + v_ext(i,k)*(1.0-open_alpha_n) + dv(i))*imask3d(i,jsize,k)
!           END DO
!           END DO
!        END IF

!        IF (radi_n .AND. radi_tangential) THEN
! !$OMP PARALLEL DO PRIVATE(a, bx, by, cx, cy, d)
!           DO k=1, ksize
!           DO i=0, isize
!              IF (.NOT. (lmask3d(i,jsize,k) .AND. lmask3d(i+1,jsize,k))) CYCLE

!              a = u(i,jsize,k) - u_old(i,jsize,k)

!              by = u_old(i,jsize-1,k) - u_old(i,jsize,k)

!              IF (radi_oblique) THEN
!                 bx = 0.5*(u_old(i-1,jsize,k) - u_old(i+1,jsize,k))
!              ELSE
!                 bx = 0.0
!              END IF

!              IF (a*by > 0.0) THEN
!                 cx = a*bx / (bx**2 + by**2)
!                 cy = a*by / (bx**2 + by**2)
!                 d = radi_nudge_n(1)
!              ELSE
!                 cx = 0.0
!                 cy = 0.0
!                 d = radi_nudge_n(2)
!              END IF

!              u(i,jsize+1,k) = u_old(i,jsize+1,k) + cy*u(i,jsize,k)

!              IF (cx > 0.0) THEN
!                 u(i,jsize+1,k) = u(i,jsize+1,k) + cx*u_old(i-1,jsize+1,k)
!              ELSE
!                 u(i,jsize+1,k) = u(i,jsize+1,k) - cx*u_old(i+1,jsize+1,k)
!              END IF

!              u(i,jsize+1,k) = u(i,jsize+1,k) / (1.0 + cy + abs(cx))

! !             IF (u_ext(i,k) /= UNDEF) u(i,jsize+1,k) = (1.0-d)*u(i,jsize+1,k) + d*u_ext(i,k)
!           END DO
!           END DO

! !$OMP PARALLEL DO PRIVATE(a, bx, by, cx, cy, d)
!           DO k=0, ksize
!           DO i=1, isize
!              IF (.NOT. (lmask3d(i,jsize,k) .AND. lmask3d(i,jsize,k+1))) CYCLE

!              a = w(i,jsize,k) - w_old(i,jsize,k)

!              by = w_old(i,jsize-1,k) - w_old(i,jsize,k)

!              IF (radi_oblique) THEN
!                 bx = 0.5*(w_old(i-1,jsize,k) - w_old(i+1,jsize,k))
!              ELSE
!                 bx = 0.0
!              END IF

!              IF (a*by > 0.0) THEN
!                 cx = a*bx / (bx**2 + by**2)
!                 cy = a*by / (bx**2 + by**2)
!                 d = radi_nudge_n(1)
!              ELSE
!                 cx = 0.0
!                 cy = 0.0
!                 d = radi_nudge_n(2)
!              END IF

!              w(i,jsize+1,k) = w_old(i,jsize+1,k) + cy*w(i,jsize,k)

!              IF (cx > 0.0) THEN
!                 w(i,jsize+1,k) = w(i,jsize+1,k) + cx*w_old(i-1,jsize+1,k)
!              ELSE
!                 w(i,jsize+1,k) = w(i,jsize+1,k) - cx*w_old(i+1,jsize+1,k)
!              END IF

!              w(i,jsize+1,k) = w(i,jsize+1,k) / (1.0 + cy + abs(cx))

! !             IF (w_ext(i,k) /= UNDEF) w(i,jsize+1,k) = (1.0-d)*w(i,jsize+1,k) + d*w_ext(i,k)
!           END DO
!           END DO
!        ELSE
! !$OMP PARALLEL DO
!           DO k=1, ksize
!           DO i=0, isize
!              IF (u_ext(i,k) == UNDEF) THEN
!                 u(i,jsize+1,k) = u(i,jsize,k)
!              ELSE
!                 u(i,jsize+1,k) = u_ext(i,k)*imask3d(i,jsize,k)*imask3d(i+1,jsize,k)
!              END IF
!           END DO
!           END DO

! !$OMP PARALLEL DO
!           DO k=0, ksize
!           DO i=1, isize
!              IF (w_ext(i,k) == UNDEF) THEN
!                 w(i,jsize+1,k) = w(i,jsize,k)
!              ELSE
!                 w(i,jsize+1,k) = w_ext(i,k)*imask3d(i,jsize,k)*imask3d(i,jsize,k+1)
!              END IF
!           END DO
!           END DO
!        END IF

     END SUBROUTINE ob_north

     SUBROUTINE ob_south
       REAL(8) :: u_ext(0:isize,1:ksize)
       REAL(8) :: v_ext(1:isize,1:ksize)
       REAL(8) :: w_ext(1:isize,0:ksize)
!!!!!!!!YJKIM
       REAL(8) :: v_barocli(1:isize,1:ksize)
       REAL(8) :: v_barotro(1:isize)
       REAL(8) :: v_barotro_s(1:isize)       
       REAL(8) :: ext_ssh_s(1:isize)
       REAL(8) :: ext_tide_s(1:isize)
       REAL(8) :: relax(1:10)
       LOGICAL :: stat1
       LOGICAL :: stat2
!!!!!!!!!!!!!
       REAL(8) :: dv(1:isize)
       REAL(8) :: a, bx, by, cx, cy, d
       INTEGER :: i, k

!$OMP PARALLEL WORKSHARE
       u_ext(:,:) = 0.0
       v_ext(:,:) = UNDEF
       w_ext(:,:) = 0.0
!!!!!!!!!YJKIM       
       v_barotro_s(:) = 0.0
!!!!!!!!!!!!!         
!$OMP END PARALLEL WORKSHARE

       CALL checkin('U_SOUTH', u_ext, time=t_current+dtime, section='XZ')
       CALL checkin('V_SOUTH', v_ext, time=t_current+dtime, section='XZ')
       CALL checkin('W_SOUTH', w_ext, time=t_current+dtime, section='XZ')
!!!!!!!!!YJKIM       
       CALL checkin('SSH_SOUTH', ext_ssh_s, time=t_current+dtime, axis='X', stat=stat1)
       CALL checkin('VTIDE_SOUTH', ext_tide_s, time=t_current+dtime, axis='X', stat=stat2)
!!!!!!!!!!!!!
       IF (jcoord/=0) RETURN
!!!!!!!!!!!!! YJKIM 1. barotro Flather / barocli relaxation (Carter and Merrifield, 2007)
!!!!!!!!!!!!!       2. barotro Flather / barocli nogradient <- current version       
!!!!!!!!!!!!! to turn off Tidal forcing, put SSH / Vel to zero by boundary condition --> passive Flather            
       IF (stat1 .AND. stat2 .AND. jcoord==0) THEN
         DO i=1, isize
            IF (lmask2d(i,2) .AND. ext_tide_s(i) /= UNDEF .AND. ext_ssh_s(i) /= UNDEF) THEN
                v_barotro_s(i) = ext_tide_s(i) - (cext(i,2)/(wct_ref(i,2)+ssh(i,2)))*(ssh(i,2)-ext_ssh_s(i)) 
            ELSE
                v_barotro_s(i) = 0.0
            END IF
         END DO

         v_barocli(:,:)=0
         v_barotro(:)=0

          DO k=1, ksize
          DO i=1, isize
             v_barotro(i) = v_barotro(i) + (v(i,3,k)*dsy(i,2,k))*imask3d(i,3,k)
          END DO
          END DO

          CALL vsum(v_barotro, all=.TRUE.)

          DO i=1, isize
             IF (dsy2d(i,2) > 0.0) THEN
                v_barotro(i) = v_barotro(i) / dsy2d(i,2)
             ELSE
                v_barotro(i) = 0.0
             END IF
          END DO

         DO k=1, ksize
         DO i=1, isize
            v_barocli(i,k)=v(i,3,k)-v_barotro(i)*imask3d(i,3,k)
         END DO
         END DO

         DO k=1, ksize
         DO i=1, isize
          u(i,2,k) = u(i,3,k)
          u(i,1,k) = u(i,2,k)
          v(i,2,k) = v_barocli(i,k)*imask3d(i,2,k)+v_barotro_s(i)*imask3d(i,2,k)
          v(i,1,k) = v(i,2,k)
          w(i,2,k) = w(i,3,k)
          w(i,1,k) = w(i,2,k)
         ENDDO
         ENDDO

         ! DO i=1,10
         !   a=(i-1)/2
         !   relax(i)=DTANH(a)
         !    ! relax(i)=0.1*(i-1)
         !    ! relax(i)=((10-i+1)/10)**2
         ! END DO
         ! DO k=1, ksize
         ! DO i=1, isize
         !  v(i,10,k)= relax(10)*v(i,10,k)+ (1-relax(10))*v_barotro_s(i)*imask3d(i,10,k) 
         !  v(i,9,k) = relax(9)*v(i,9,k) + (1-relax(9))*v_barotro_s(i)*imask3d(i,9,k) 
         !  v(i,8,k) = relax(8)*v(i,8,k) + (1-relax(8))*v_barotro_s(i)*imask3d(i,8,k) 
         !  v(i,7,k) = relax(7)*v(i,7,k) + (1-relax(7))*v_barotro_s(i)*imask3d(i,7,k) 
         !  v(i,6,k) = relax(6)*v(i,6,k) + (1-relax(6))*v_barotro_s(i)*imask3d(i,6,k) 
         !  v(i,5,k) = relax(5)*v(i,5,k) + (1-relax(5))*v_barotro_s(i)*imask3d(i,5,k) 
         !  v(i,4,k) = relax(4)*v(i,4,k) + (1-relax(4))*v_barotro_s(i)*imask3d(i,4,k) 
         !  v(i,3,k) = relax(3)*v(i,3,k) + (1-relax(3))*v_barotro_s(i)*imask3d(i,3,k) 
         !  v(i,2,k) = relax(2)*v(i,2,k) + (1-relax(2))*v_barotro_s(i)*imask3d(i,2,k) 
         !  v(i,1,k)=  v_barotro_s(i)*imask3d(i,1,k)
         !  u(i,2,k) = u(i,3,k)
         !  u(i,1,k) = u(i,2,k)
         !  w(i,12,k)= relax(10)*w(i,12,k)
         !  w(i,11,k)= relax(9)*w(i,11,k)
         !  w(i,10,k)= relax(8)*w(i,10,k)
         !  w(i,9,k)= relax(7)*w(i,9,k)
         !  w(i,8,k)= relax(6)*w(i,8,k)
         !  w(i,7,k)= relax(5)*w(i,7,k)
         !  w(i,6,k)= relax(4)*w(i,6,k)
         !  w(i,5,k)= relax(3)*w(i,5,k)
         !  w(i,4,k)= relax(2)*w(i,4,k)
         !  w(i,3,k) = 0
         !  w(i,2,k) = w(i,3,k)
         !  w(i,1,k) = w(i,2,k)
         ! ENDDO
         ! ENDDO
       END IF     
! !$OMP PARALLEL DO PRIVATE(a, bx, by, cx, cy, d)
!           DO k=1, ksize
!           DO i=1, isize
!              IF (.NOT. lmask3d(i,1,k)) CYCLE

!             !  a = v(i,1,k) - v_old(i,1,k)

!             !  by = v_old(i,2,k) - v_old(i,1,k)

!             !  IF (radi_oblique) THEN
!             !     bx = 0.5*(v_old(i-1,1,k) - v_old(i+1,1,k))
!             !  ELSE
!             !     bx = 0.0
!             !  END IF

!             !  IF (a*by > 0.0) THEN
!             !     cx = a*bx/(bx**2 + by**2)
!             !     cy = a*by/(bx**2 + by**2)
!             !     d = radi_nudge_s(1)
!             !  ELSE
!             !     cx = 0.0
!             !     cy = 0.0
!             !     d = radi_nudge_s(2)
!             !  END IF

!             !  v(i,0,k) = v_old(i,0,k) + cy*v(i,1,k)

!             !  IF (cx > 0.0) THEN
!             !     v(i,0,k) = v(i,0,k) + cx*v_old(i-1,0,k)
!             !  ELSE
!             !     v(i,0,k) = v(i,0,k) - cx*v_old(i+1,0,k)
!             !  END IF

!             !  v(i,0,k) = v(i,0,k) / (1.0 + cy + abs(cx))

! !             IF (v_ext(i,k) /= UNDEF) v(i,0,k) = (1.0-d)*v(i,0,k) + d*v_ext(i,k)
!              v(i,1,k) = v(i,2,k) 
!              v(i,0,k) = v(i,1,k)

!              u(i,1,k) = u(i,2,k) 
!              u(i,0,k) = u(i,1,k)
!              w(i,1,k) = w(i,2,k) 
!              w(i,0,k) = w(i,1,k)
!           END DO
!           END DO
!        ELSE
!           dv(:) = 0.0
!           DO k=1, ksize
!              DO i=1, isize
!                 IF (v_ext(i,k) == UNDEF) v_ext(i,k) = 0.0
!                 dv(i) = dv(i) + (v_ext(i,k)*dsy_ref(i,0,k) - (v(i,1,k)*open_alpha_s + v_ext(i,k)*(1.0-open_alpha_s))*dsy(i,0,k))*imask3d(i,1,k)
!              END DO
!           END DO

!           CALL vsum(dv, all=.TRUE.)

!           DO i=1, isize
!              IF (dsy2d(i,0) > 0.0) THEN
!                 dv(i) = dv(i) / dsy2d(i,0)
!              ELSE
!                 dv(i) = 0.0
!              END IF
!           END DO

! !$OMP PARALLEL DO
!           DO k=1, ksize
!           DO i=1, isize
!              v(i, 0,k) = (v(i,1,k)*open_alpha_s + v_ext(i,k)*(1.0-open_alpha_s) + dv(i))*imask3d(i,1,k)
!           END DO
!           END DO
!        END IF

!        IF (radi_s .AND. radi_tangential) THEN
! !$OMP PARALLEL DO PRIVATE(a, bx, by, cx, cy, d)
!           DO k=1, ksize
!           DO i=0, isize
!              IF (.NOT. (lmask3d(i,1,k) .AND. lmask3d(i+1,1,k))) CYCLE

!              a = u(i,1,k) - u_old(i,1,k)

!              by = u_old(i,2,k) - u_old(i,1,k)

!              IF (radi_oblique) THEN
!                 bx = 0.5*(u_old(i-1,1,k) - u_old(i+1,1,k))
!              ELSE
!                 bx = 0.0
!              END IF

!              IF (a*by > 0.0) THEN
!                 cx = a*bx / (bx**2 + by**2)
!                 cy = a*by / (bx**2 + by**2)
!                 d = radi_nudge_s(1)
!              ELSE
!                 cx = 0.0
!                 cy = 0.0
!                 d = radi_nudge_s(2)
!              END IF

!              u(i,0,k) = u_old(i,0,k) + cy*u(i,1,k)

!              IF (cx > 0.0) THEN
!                 u(i,0,k) = u(i,0,k) + cx*u_old(i-1,0,k)
!              ELSE
!                 u(i,0,k) = u(i,0,k) - cx*u_old(i+1,0,k)
!              END IF

!              u(i,0,k) = u(i,0,k) / (1.0 + cy + abs(cx))

! !             IF (u_ext(i,k) /= UNDEF) u(i,0,k) = (1.0-d)*u(i,0,k) + d*u_ext(i,k)
!           END DO
!           END DO

! !$OMP PARALLEL DO PRIVATE(a, bx, by, cx, cy, d)
!           DO k=0, ksize
!           DO i=1, isize
!              IF (.NOT. (lmask3d(i,1,k) .AND. lmask3d(i,1,k+1))) CYCLE

!              a = w(i,1,k) - w_old(i,1,k)

!              by = w_old(i,2,k) - w_old(i,1,k)

!              IF (radi_oblique) THEN
!                 bx = 0.5*(w_old(i-1,1,k) - w_old(i+1,1,k))
!              ELSE
!                 bx = 0.0
!              END IF

!              IF (a*by > 0.0) THEN
!                 cx = a*bx / (bx**2 + by**2)
!                 cy = a*by / (bx**2 + by**2)
!                 d = radi_nudge_s(1)
!              ELSE
!                 cx = 0.0
!                 cy = 0.0
!                 d = radi_nudge_s(2)
!              END IF

!              w(i,0,k) = w_old(i,0,k) + cy*w(i,1,k)

!              IF (cx > 0.0) THEN
!                 w(i,0,k) = w(i,0,k) + cx*w_old(i-1,0,k)
!              ELSE
!                 w(i,0,k) = w(i,0,k) - cx*w_old(i+1,0,k)
!              END IF

!              w(i,0,k) = w(i,0,k) / (1.0 + cy + abs(cx))

! !             IF (w_ext(i,k) /= UNDEF) w(i,0,k) = (1.0-d)*w(i,0,k) + d*w_ext(i,k)
!           END DO
!           END DO
!        ELSE
! !$OMP PARALLEL DO
!           DO k=1, ksize
!           DO i=0, isize
!              IF (u_ext(i,k) == UNDEF) THEN
!                 u(i,0,k) = u(i,1,k)
!              ELSE
!                 u(i,0,k) = u_ext(i,k)*imask3d(i,1,k)*imask3d(i+1,1,k)
!              END IF
!           END DO
!           END DO

! !$OMP PARALLEL DO
!           DO k=0, ksize
!           DO i=1, isize
!              IF (w_ext(i,k) == UNDEF) THEN
!                 w(i,0,k) = w(i,1,k)
!              ELSE
!                 w(i,0,k) = w_ext(i,k)*imask3d(i,1,k)*imask3d(i,1,k+1)
!              END IF
!           END DO
!           END DO
!        END IF

     END SUBROUTINE ob_south

     SUBROUTINE ob_bottom
       REAL(8) :: w_ext(1:isize, 1:jsize)

       INTEGER :: i, j

       w_ext(:,:) = 0.0
       CALL checkin('W_BOTTOM', w_ext(:,:), time=t_current+dtime)

       DO j=1, jsize
       DO i=1, isize
          IF (.NOT. bottom_flag(i,j)) CYCLE
          w(i,j,bottom_k(i,j)-1) = w_ext(i,j)

       END DO
       END DO
     END SUBROUTINE ob_bottom
!!!flather
     SUBROUTINE tidal
       REAL(8) :: tmp_x(1:isize)
       REAL(8) :: tmp_x_vel(1:isize)
       REAL(8) :: tmp_y(1:jsize)
       REAL(8) :: tmp_y_vel(1:jsize)
       REAL(8) :: d_e(1:jsize)
       REAL(8) :: d_w(1:jsize)
       REAL(8) :: d_n(1:isize)
       REAL(8) :: d_s(1:isize)

       INTEGER :: i, j, k
       LOGICAL :: stat


       IF (open_e) THEN
         tmp_y(:) = 0.0
         tmp_y_vel(:) = 0.0
          CALL checkin('SSH_EAST', tmp_y, time=t_current+dtime, axis='Y', stat=stat)
          CALL checkin('UTIDE_EAST', tmp_y_vel, time=t_current+dtime, axis='Y', stat=stat)
          d_e(:) = 0.0
          IF (stat .AND. icoord==ipes-1) THEN
!          IF (icoord==ipes-1) THEN
             DO k=1, ksize
                DO j=1, jsize
                   d_e(j) = d_e(j) + u(isize-1,j,k)*dsx(isize-1,j,k)-u(isize,j,k)*dsx(isize,j,k) &
                                   + v(isize,j-1,k)*dsy(isize,j-1,k)-v(isize,j,k)*dsy(isize,j,k)
                END DO
             END DO
             CALL vsum(d_e, all=.TRUE.)
             DO j=1, jsize
                d_e(j) = d_e(j)*dtime / dsz(isize,j)
             END DO

             DO j=1, jsize
                IF (lmask2d(isize,j) .AND. tmp_y(j) /= UNDEF) THEN                  
                   d_e(j) = (ssh(isize,j) + d_e(j) - tmp_y(j))*cext(isize,j)*dsz(isize,j) / (dsx2d(isize,j)*(dx(isize,j) + cext(isize,j)*dtime))
                  !  d_e(j) = tmp_y_vel(j) + (ssh(isize,j) - tmp_y(j))*cext(isize,j)*dsz(isize,j) / (dsx2d(isize,j)*(dx(isize,j) + cext(isize,j)*dtime))
!                   d_e(j) = tmp_y_vel(j) + (cext(isize,j)/(wct_ref(isize,j)+ssh(isize,j)))*(ssh(isize,j)-tmp_y(j))
                ELSE
                   d_e(j) = 0.0
                END IF
             END DO

             IF (open_n .AND. jcoord==jpes-1) d_e(jsize) = 0.5*d_e(jsize)
             IF (open_s .AND. jcoord==0)      d_e(1)     = 0.5*d_e(1)
          END IF
       END IF
!!!!flather
       IF (open_w) THEN
          CALL checkin('SSH_WEST', tmp_y, time=t_current+dtime, axis='Y', stat=stat)

          d_w(:) = 0.0
          IF (stat .AND. icoord==0) THEN
             DO k=1, ksize
                DO j=1, jsize
                   d_w(j) = d_w(j) + u(0,j,  k)*dsx(0,j,  k)-u(1,j,k)*dsx(1,j,k) &
                                   + v(1,j-1,k)*dsy(1,j-1,k)-v(1,j,k)*dsy(1,j,k)
                END DO
             END DO
             CALL vsum(d_w, all=.TRUE.)
             DO j=1, jsize
                d_w(j) = d_w(j)*dtime / dsz(1,j)
             END DO

             DO j=1, jsize
                IF (lmask2d(isize,j) .AND. tmp_y(j) /= UNDEF) THEN
                   d_w(j) = -(ssh(1,j) + d_w(j) - tmp_y(j))*cext(1,j)*dsz(1,j) / (dsx2d(0,j)*(dx(1,j) + cext(1,j)*dtime))
                ELSE
                   d_w(j) = 0.0
                END IF
             END DO

             IF (open_n .AND. jcoord==jpes-1) d_w(jsize) = 0.5*d_w(jsize)
             IF (open_s .AND. jcoord==0)      d_w(1)     = 0.5*d_w(1)
          END IF
       END IF

       IF (open_n) THEN
         tmp_x(:) = 0.0
         tmp_x_vel(:) = 0.0
          CALL checkin('SSH_NORTH', tmp_x, time=t_current+dtime, axis='X', stat=stat)
          CALL checkin('VTIDE_NORTH', tmp_x_vel, time=t_current+dtime, axis='X', stat=stat)

          d_n(:) = 0.0
          IF (stat .AND. jcoord==jpes-1) THEN
             DO k=1, ksize
                DO i=1, isize
                   d_n(i) = d_n(i) + u(i-1,jsize,k)*dsx(i-1,jsize,k)-u(i,jsize,k)*dsx(i,jsize,k) &
                                   + v(i,jsize-1,k)*dsy(i,jsize-1,k)-v(i,jsize,k)*dsy(i,jsize,k)
                END DO
             END DO
             CALL vsum(d_n, all=.TRUE.)
             DO i=1, isize
                d_n(i) = d_n(i)*dtime / dsz(i,jsize)
             END DO

             DO i=1, isize
                IF (lmask2d(i,jsize).AND. tmp_x(i) /= UNDEF) THEN
                    d_n(i) = (ssh(i,jsize)+ d_n(i) - tmp_x(i))*cext(i,jsize)*dsz(i,jsize) / (dsy2d(i,jsize)*(dy(i,jsize) + cext(i,jsize)*dtime))
                  !  d_n(i) = tmp_x_vel(i) + (ssh(i,jsize) - tmp_x(i))*cext(i,jsize)*dsz(i,jsize) / (dsy2d(i,jsize)*(dy(i,jsize) + cext(i,jsize)*dtime))
!                   d_n(i) = tmp_x_vel(i) + (cext(i,jsize)/(wct_ref(i,jsize)+ssh(i,jsize)))*(ssh(i,jsize)-tmp_x(i))
                ELSE
                   d_n(i) = 0.0
                END IF
             END DO

             IF (open_e .AND. icoord==ipes-1) d_n(isize) = 0.5*d_n(isize)
             IF (open_w .AND. icoord==0)      d_n(1)     = 0.5*d_n(1)
          END IF
       END IF

       IF (open_s) THEN
          tmp_x(:) = 0.0
          tmp_x_vel(:) = 0.0
          CALL checkin('SSH_SOUTH', tmp_x, time=t_current+dtime, axis='X', stat=stat)
          CALL checkin('VTIDE_SOUTH', tmp_x, time=t_current+dtime, axis='X', stat=stat)
          d_s(:) = 0.0
          IF (stat .AND. jcoord==0) THEN
             DO k=1, ksize
                DO i=1, isize
                   d_s(i) = d_s(i) + u(i-1,1,k)*dsx(i-1,1,k)-u(i,1,k)*dsx(i,1,k) &
                                   + v(i,  0,k)*dsy(i,  0,k)-v(i,1,k)*dsy(i,1,k)
                END DO
             END DO
             CALL vsum(d_s, all=.TRUE.)
             DO i=1, isize
                d_s(i) = d_s(i)*dtime / dsz(i,1)
             END DO

             DO i=1, isize
                IF (lmask2d(i,1) .AND. tmp_x(i) /= UNDEF) THEN
                    d_s(i) = -(ssh(i,1) + d_s(i) - tmp_x(i))*cext(i,1)*dsz(i,1) / (dsy2d(i,0)*(dy(i,1) + cext(i,1)*dtime))
                  !  d_s(i) = tmp_x_vel(i)-(ssh(i,1) - tmp_x(i))*cext(i,1)*dsz(i,1) / (dsy2d(i,0)*(dy(i,1) + cext(i,1)*dtime))
!                   d_s(i) = tmp_x_vel(i) - (cext(i,1)/(wct_ref(i,1)+ssh(i,1)))*(ssh(i,1)-tmp_x(i))
                  ELSE
                   d_s(i) = 0.0
                END IF
             END DO

             IF (open_e .AND. icoord==ipes-1) d_s(isize) = 0.5*d_s(isize)
             IF (open_w .AND. icoord==0)      d_s(1)     = 0.5*d_s(1)
          END IF
       END IF

       IF (open_e .AND. icoord==ipes-1) THEN
!$OMP PARALLEL DO
          DO k=1, ksize
          DO j=1, jsize
             u(isize,j,k) = u(isize,j,k) + d_e(j)*imask3d(isize,j,k)
          END DO
          END DO
       END IF

       IF (open_w .AND. icoord==0) THEN
!$OMP PARALLEL DO
          DO k=1, ksize
          DO j=1, jsize
             u(0,j,k) = u(0,j,k) + d_w(j)*imask3d(1,j,k)
          END DO
          END DO
       END IF

       IF (open_n .AND. jcoord==jpes-1) THEN
!$OMP PARALLEL DO
          DO k=1, ksize
          DO i=1, isize
            v(i,jsize,k) = v(i,jsize,k) + d_n(i)*imask3d(i,jsize,k)
          END DO
          END DO
       END IF

       IF (open_s .AND. jcoord==0) THEN
!$OMP PARALLEL DO
          DO k=1, ksize
          DO i=1, isize
             v(i,0,k) = v(i,0,k) + d_s(i)*imask3d(i,1,k)
          END DO
          END DO
       END IF

     END SUBROUTINE tidal

   END SUBROUTINE openboundary_velocity

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE update_velocity_external
    INTEGER :: k

    IF (open_e .AND. icoord==ipes-1) THEN
       CALL update_boundary_yz(u(isize,  :,:))
       CALL update_boundary_yz(v(isize+1,:,:))
       CALL update_boundary_yz(w(isize+1,:,:))
    END IF

    IF (open_w .AND. icoord==0) THEN
       CALL update_boundary_yz(u(0,:,:))
       CALL update_boundary_yz(v(0,:,:))
       CALL update_boundary_yz(w(0,:,:))
    END IF

    IF (open_n .AND. jcoord==jpes-1) THEN
       CALL update_boundary_xz(u(:,jsize+1,:))
       CALL update_boundary_xz(v(:,jsize,  :))
       CALL update_boundary_xz(w(:,jsize+1,:))
    END IF

    IF (open_s .AND. jcoord==0) THEN
       CALL update_boundary_xz(u(:,0,:))
       CALL update_boundary_xz(v(:,0,:))
       CALL update_boundary_xz(w(:,0,:))
    END IF

    IF (open_e .AND. icoord==ipes-1) THEN
!$OMP PARALLEL WORKSHARE
!!!!!!!!YJKIM
!       u(isize,:,:)   = u(isize-1,:,:)     
       u(isize+1,:,:) = u(isize,  :,:)
       u(isize+2,:,:) = u(isize,  :,:)

!       v(isize,:,:)   = v(isize-1,:,:)     
       v(isize+1,:,:) = v(isize,  :,:)
       v(isize+2,:,:) = v(isize+1,:,:)

!       w(isize,:,:)   = w(isize-1,:,:)     
       w(isize+1,:,:) = w(isize,  :,:)
       w(isize+2,:,:) = w(isize+1,:,:)
!$OMP END PARALLEL WORKSHARE
    END IF

    IF (open_w .AND. icoord==0) THEN
!$OMP PARALLEL WORKSHARE
       u(-1,:,:) = u(0,:,:)
       u(-2,:,:) = u(0,:,:)

       v(-1,:,:) = v(0,:,:)
       w(-1,:,:) = w(0,:,:)
!$OMP END PARALLEL WORKSHARE
    END IF

    IF (open_n .AND. jcoord==jpes-1) THEN
!$OMP PARALLEL WORKSHARE
!!!!!!!!YJKIM
!       v(:,jsize,:)   = v(:,jsize-1,:)       
       v(:,jsize+1,:) = v(:,jsize,:)
       v(:,jsize+2,:) = v(:,jsize,:)

!       u(:,jsize,:)   = u(:,jsize-1,:)       
       u(:,jsize+1,:) = u(:,jsize,:)
       u(:,jsize+2,:) = u(:,jsize+1,:)

!       w(:,jsize,:)   = w(:,jsize-1,:)       
       w(:,jsize+1,:) = w(:,jsize,:)
       w(:,jsize+2,:) = w(:,jsize+1,:)
!$OMP END PARALLEL WORKSHARE
    END IF

    IF (open_s .AND. jcoord==0) THEN
!$OMP PARALLEL WORKSHARE
!!!!!!!YJKIM
!       v(:,1,:) = v(:,2,:)
       v(:,0,:) = v(:,1,:)      
       v(:,-1,:) = v(:,0,:)
       v(:,-2,:) = v(:,0,:)

!       u(:,1,:) = u(:,2,:)
       u(:,0,:) = u(:,1,:)      
       u(:,-1,:) = u(:,0,:)
!       u(:,-2,:) = u(:,0,:)

!       w(:,1,:) = w(:,2,:) 
       w(:,0,:) = w(:,1,:)      
       w(:,-1,:) = w(:,0,:)
!       w(:,-2,:) = w(:,0,:)
!$OMP END PARALLEl WORKSHARE
    END IF

    IF (open_e .AND. open_n .AND. icoord==ipes-1 .AND. jcoord==jpes-1) THEN
!$OMP PARALLEL DO
       DO k=0, ksize
          w(isize+1:isize+2,jsize+1:jsize+2,k) = 0.5*(w(isize+1,jsize,k) + w(isize,  jsize+1,k))

          IF (k==0) CYCLE

          u(isize  :isize+2,jsize+1:jsize+2,k) = 0.5*(u(isize+1,jsize,k) + u(isize-1,jsize+1,k))
          v(isize+1:isize+2,jsize  :jsize+2,k) = 0.5*(v(isize,jsize+1,k) + v(isize+1,jsize-1,k))
       END DO
    END IF

    IF (open_w .AND. open_n .AND. icoord==0 .AND. jcoord==jpes-1) THEN
!$OMP PARALLEL DO
       DO k=0, ksize
          w(-1:0,jsize+1:jsize+2,k) = 0.5*(w( 0,jsize,  k) + w(1,jsize+1,k))

          IF (k==0) CYCLE

          u(-2:0,jsize+1:jsize+2,k) = 0.5*(u(-1,jsize,  k) + u(1,jsize+1,k))
          v(-1:0,jsize  :jsize+2,k) = 0.5*(v( 0,jsize-1,k) + v(1,jsize+1,k))
       END DO
    END IF

    IF (open_e .AND. open_s .AND. icoord==ipes-1 .AND. jcoord==0) THEN
!$OMP PARALLEL DO
       DO k=0, ksize
          w(isize+1:isize+2,-1:0,k) = 0.5*(w(isize+1,1,k) + w(isize,  0,k))

          IF (k==0) CYCLE

          u(isize  :isize+2,-1:0,k) = 0.5*(u(isize+1,1,k) + u(isize-1,0,k))
          v(isize+1:isize+2,-2:0,k) = 0.5*(v(isize+1,1,k) + v(isize, -1,k))
       END DO
    END IF

    IF (open_w .AND. open_s .AND. icoord==0 .AND. jcoord==0) THEN
!$OMP PARALLEL DO
       DO k=0, ksize
          w(-1:0,-1:0,k) = 0.5*(w(0, 1,k) + w(1, 0,k))

          IF (k==0) CYCLE

          u(-2:0,-1:0,k) = 0.5*(u(-1,1,k) + u(1, 0,k))
          v(-1:0,-2:0,k) = 0.5*(v( 0,1,k) + v(1,-1,k))
       END DO
    END IF

  END SUBROUTINE update_velocity_external

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE update_velocity_diagnostic(no_checkout)
    LOGICAL, INTENT(IN), OPTIONAL :: no_checkout

    REAL(8) :: tmp(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv, 4)
    INTEGER :: i, j, k
    INTEGER :: kl, km

    CALL velocity_gradient(u, v, w, dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz)
    CALL kinetic_energy(u, v, w, ke)
    CALL strain_rate(dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz, srate)
    CALL shear_freq2(dudz, dvdz, sfreq2)

    IF (present(no_checkout)) THEN
       IF (no_checkout) RETURN
    END IF

    CALL checkout('KE',    ke)    ! kinetic energy (grid-resolved-scale) per unit mass [m^2/s^2]
    CALL checkout('SRATE', srate) ! strain rate, sqrt(|S_ij|^2) [1/s]

    IF (require_checkout('DIV')) THEN
!$OMP PARALLEL DO
       DO k=1, ksize
          kl = dzindex(k)
          km = maskindex(k)

          DO j=1, jsize
          DO i=1, isize
             tmp(i,j,k,1) = imask3d(i,j,km)* &
                  ( ( u(i,j,k)*dsx_old(i,j,kl) - u(i-1,j,k)*dsx_old(i-1,j,kl) &
                    + v(i,j,k)*dsy_old(i,j,kl) - v(i,j-1,k)*dsy_old(i,j-1,kl) &
                    + w(i,j,k)*dsz(i,j)        - w(i,j,k-1)*dsz(i,j)) * a_ebcn &
                  + ( u_old(i,j,k)*dsx_old(i,j,kl) - u_old(i-1,j,k)*dsx_old(i-1,j,kl) &
                    + v_old(i,j,k)*dsy_old(i,j,kl) - v_old(i,j-1,k)*dsy_old(i,j-1,kl) &
                    + w_old(i,j,k)*dsz(i,j)        - w_old(i,j,k-1)*dsz(i,j)) * (1.0-a_ebcn) &
                  + (dvol(i,j,kl) - dvol_old(i,j,kl))/dtime) / dvol(i,j,kl)
          END DO
          END DO
       END DO

       IF (use_landwater .AND. vrank==0) THEN
          DO j=1, jsize
          DO i=1, isize
             IF (lwdried(i,j)) tmp(i,j,ksize,1) = 0.0
          END DO
          END DO
       END IF

       CALL checkout('DIV', tmp(:,:,:,1)) ! fake divergence (due to numerical/truncate error) [1/s]
    END IF

    IF (require_checkout('UV') .OR. require_checkout('UW') .OR. require_checkout('VW')) THEN
!$OMP PARALLEL DO
       DO k=1, ksize
       DO j=1, jsize
       DO i=1, isize
          tmp(i,j,k,1) = 0.25*(u(i,j,k)+u(i,j+1,k))*(v(i,j,k)+v(i+1,j,k))
          tmp(i,j,k,2) = 0.25*(u(i,j,k)+u(i,j,k+1))*(w(i,j,k)+w(i+1,j,k))
          tmp(i,j,k,3) = 0.25*(v(i,j,k)+v(i,j,k+1))*(w(i,j,k)+w(i,j+1,k))
       END DO
       END DO
       END DO
       CALL checkout('UV', tmp(:,:,:,1))
       CALL checkout('UW', tmp(:,:,:,2))
       CALL checkout('VW', tmp(:,:,:,3))
    END IF

    IF (require_checkout('ENSTROPHY')) THEN
       CALL enstrophy(tmp(:,:,:,1))
       CALL checkout('ENSTROPHY', tmp(:,:,:,1))
    END IF

    IF (require_checkout('PRPLS')) THEN
       CALL principal_strain_rate(tmp(:,:,:,1), tmp(:,:,:,2), tmp(:,:,:,3), tmp(:,:,:,4))
       CALL checkout('PRPLS', tmp(:,:,:,1))   ! principal strain rate, [1/s]
       CALL checkout('PRPLV_X', tmp(:,:,:,2)) ! x-compornent of principal axis for strain rate
       CALL checkout('PRPLV_Y', tmp(:,:,:,3)) ! y-compornent of principal axis for strain rate
       CALL checkout('PRPLV_Z', tmp(:,:,:,4)) ! z-compornent of principal axis for strain rate
    END IF

    IF (require_checkout('LAMBDA2')) THEN
       CALL lambda2(tmp(:,:,:,1))
       CALL checkout('LAMBDA2', tmp(:,:,:,1)) ! labmda-2, the second eigenvalue for S^2 + \Omega^2, see Jeong and Hussain 1995 JFM, [1/s^2]
    END IF

    IF (require_checkout('OWPARAM')) THEN
       CALL owparam(tmp(:,:,:,1))
       CALL checkout('OWPARAM', tmp(:,:,:,1)) ! The Okubo-Weiss parameter, see Okubo 1970 DSR, [1/s^2]
    END IF
    IF (require_checkout('OWPARAM_XZ')) THEN
       CALL owparam_xz(tmp(:,:,:,2))
       CALL checkout('OWPARAM_XZ', tmp(:,:,:,2))
    END IF
    IF (require_checkout('OWPARAM_YZ')) THEN
       CALL owparam_yz(tmp(:,:,:,3))
       CALL checkout('OWPARAM_YZ', tmp(:,:,:,3))
    END IF


  END SUBROUTINE update_velocity_diagnostic

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE velocity_gradient(u, v, w, dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz)
    REAL(8), INTENT(IN)  :: u( -slv:isize+slv,1-slv:jsize+slv,1-slv:ksize+slv)
    REAL(8), INTENT(IN)  :: v(1-slv:isize+slv, -slv:jsize+slv,1-slv:ksize+slv)
    REAL(8), INTENT(IN)  :: w(1-slv:isize+slv,1-slv:jsize+slv, -slv:ksize+slv)

    REAL(8), INTENT(OUT) :: dudx( 0:isize+1, 0:jsize+1, 0:ksize+1)
    REAL(8), INTENT(OUT) :: dudy(-1:isize+1,-1:jsize+1, 0:ksize+1)
    REAL(8), INTENT(OUT) :: dudz(-1:isize+1, 0:jsize+1,-1:ksize+1)

    REAL(8), INTENT(OUT) :: dvdx(-1:isize+1,-1:jsize+1, 0:ksize+1)
    REAL(8), INTENT(OUT) :: dvdy( 0:isize+1, 0:jsize+1, 0:ksize+1)
    REAL(8), INTENT(OUT) :: dvdz( 0:isize+1,-1:jsize+1,-1:ksize+1)

    REAL(8), INTENT(OUT) :: dwdx(-1:isize+1, 0:jsize+1,-1:ksize+1)
    REAL(8), INTENT(OUT) :: dwdy( 0:isize+1,-1:jsize+1,-1:ksize+1)
    REAL(8), INTENT(OUT) :: dwdz( 0:isize+1, 0:jsize+1, 0:ksize+1)

    INTEGER :: i, j, k
    INTEGER :: km0, km1

    REAL(4) :: sf_u, sf_l, sf_s

    sf_l = slipfactor_bottom
    sf_s = slipfactor_side

!$OMP PARALLEL private (i, j, k, km0, km1, sf_u)
!$OMP DO
    DO k=-1, ksize+1
       km0=maskindex(k)
       km1=maskindex(k+1)

       IF (kcoord==kpes-1 .AND. k==ksize) THEN
          sf_u = slipfactor_surface
       ELSE
          sf_u = slipfactor_top
       END IF

       DO j= 0, jsize+1
       DO i=-1, isize+1
          dudz(i,j,k) = masked_shear(u(i,j,k+1), u(i,j,k), imask3d(i,  j,km1)*imask3d(i+1,j,km1), &
                                                           imask3d(i,  j,km0)*imask3d(i+1,j,km0), sf_u, sf_l)*idz1(k)
          dwdx(i,j,k) = masked_shear(w(i+1,j,k), w(i,j,k), imask3d(i+1,j,km0)*imask3d(i+1,j,km1), &
                                                           imask3d(i,  j,km0)*imask3d(i,  j,km1), sf_s, sf_s)*idx1(i,j)
       END DO
       END DO

       DO j=-1, jsize+1
       DO i= 0, isize+1
          dvdz(i,j,k) = masked_shear(v(i,j,k+1), v(i,j,k), imask3d(i,j,  km1)*imask3d(i,j+1,km1), &
                                                           imask3d(i,j,  km0)*imask3d(i,j+1,km0), sf_u, sf_l)*idz1(k)
          dwdy(i,j,k) = masked_shear(w(i,j+1,k), w(i,j,k), imask3d(i,j+1,km0)*imask3d(i,j+1,km1), &
                                                           imask3d(i,j,  km0)*imask3d(i,j,  km1), sf_s, sf_s)*idy1(i,j)
       END DO
       END DO

       IF (k==-1) CYCLE

       DO j=-1, jsize+1
       DO i=-1, isize+1
          dudy(i,j,k) = masked_shear(u(i,j+1,k), u(i,j,k), imask3d(i,j+1,km0)*imask3d(i+1,j+1,km0), &
                                                           imask3d(i,j,  km0)*imask3d(i+1,j,  km0), sf_s, sf_s)*0.5*(idy1(i,j)+idy1(i+1,j)) &
                        - 0.5*(v(i,j,k)+v(i+1,j,k))*metyx(i,j)

          dvdx(i,j,k) = masked_shear(v(i+1,j,k), v(i,j,k), imask3d(i+1,j,km0)*imask3d(i+1,j+1,km0), &
                                                           imask3d(i,  j,km0)*imask3d(i,  j+1,km0), sf_s, sf_s)*0.5*(idx1(i,j)+idx1(i,j+1)) &
                        - 0.5*(u(i,j,k)+u(i,j+1,k))*metxy(i,j)
       END DO
       END DO

       DO j=0, jsize+1
       DO i=0, isize+1
          dudx(i,j,k) = (u(i,j,k)-u(i-1,j,k))*idx0(i,j) + 0.25*(v(i,j-1,k)*(metxy(i-1,j-1)+metxy(i,j-1)) + v(i,j,k)*(metxy(i-1,j)+metxy(i,j)))
          dvdy(i,j,k) = (v(i,j,k)-v(i,j-1,k))*idy0(i,j) + 0.25*(u(i-1,j,k)*(metyx(i-1,j-1)+metyx(i-1,j)) + u(i,j,k)*(metyx(i,j-1)+metyx(i,j)))
          dwdz(i,j,k) = (w(i,j,k)-w(i,j,k-1))*idz0(k)
       END DO
       END DO
    END DO

    IF (open_e .AND. icoord==ipes-1) THEN
!$OMP DO
       DO k=-1, ksize+1
          km0=maskindex(k)
          km1=maskindex(k+1)

          IF (kcoord==kpes-1 .AND. k==ksize) THEN
             sf_u = slipfactor_surface
          ELSE
             sf_u = slipfactor_top
          END IF

          DO j= 0, jsize+1
             dudz(isize,j,k) = masked_shear(u(isize,j,k+1), u(isize,j,k), imask3d(isize,j,km1), imask3d(isize,j,km0), sf_u, sf_l)*idz1(k)
             dwdx(isize,j,k) = (w(isize+1,j,k) - w(isize,j,k)) * imask3d(isize,j,km0)*imask3d(isize,j,km1) * idx1(isize,j)

             dudz(isize+1,j,k) = 0.0
             dwdx(isize+1,j,k) = 0.0
          END DO

          DO j=-1, jsize+1
             dvdz(isize+1,j,k) = masked_shear(v(isize+1,j,k+1), v(isize+1,j,k), imask3d(isize,j,  km1)*imask3d(isize,j+1,km1), &
                                                                                imask3d(isize,j,  km0)*imask3d(isize,j+1,km0), sf_u, sf_l)*idz1(k)
             dwdy(isize+1,j,k) = masked_shear(w(isize+1,j+1,k), w(isize+1,j,k), imask3d(isize,j+1,km0)*imask3d(isize,j+1,km1), &
                                                                                imask3d(isize,j,  km0)*imask3d(isize,j,  km1), sf_s, sf_s)*idy1(isize+1,j)
          END DO

          IF (k==-1) CYCLE

          DO j=-1, jsize+1
             dudy(isize,j,k) = masked_shear(u(isize,j+1,k), u(isize,j,k), imask3d(isize,j+1,km0), &
                                                                          imask3d(isize,j,  km0), sf_s, sf_s) * 0.5*(idy1(isize,j)+idy1(isize+1,j)) &
                               - 0.5*(v(isize,j,k)+v(isize+1,j,k))*metyx(isize,j)

             dvdx(isize,j,k) = (v(isize+1,j,k) - v(isize,j,k)) * imask3d(isize,j,km0)*imask3d(isize,j+1,km0)  * 0.5*(idx1(isize,j)+idx1(isize,j+1)) &
                               - 0.5*(u(isize,j,k)+u(isize,j+1,k))*metxy(isize,j)

             dudy(isize+1,j,k) = 0.0
             dvdx(isize+1,j,k) = 0.0
          END DO

          DO j=0, jsize+1
             dudx(isize+1,j,k) = 0.0
             dvdy(isize+1,j,k) = 0.0
             dwdz(isize+1,j,k) = 0.0
          END DO
       END DO
    END IF

    IF (open_w .AND. icoord==0) THEN
!$OMP DO
       DO k=-1, ksize+1
          km0=maskindex(k)
          km1=maskindex(k+1)

          IF (kcoord==kpes-1 .AND. k==ksize) THEN
             sf_u = slipfactor_surface
          ELSE
             sf_u = slipfactor_top
          END IF

          DO j= 0, jsize+1
             dudz(0,j,k) = masked_shear(u(0,j,k+1), u(0,j,k), imask3d(1,j,km1), imask3d(1,  j,km0), sf_u, sf_l)*idz1(k)
             dwdx(0,j,k) = (w(1,j,k) - w(0,j,k)) * imask3d(1,j,km0)*imask3d(1,j,km1) * idx1(0,j)

             dudz(-1,j,k) = 0.0
             dwdx(-1,j,k) = 0.0
          END DO

          DO j=-1, jsize+1
             dvdz(0,j,k) = masked_shear(v(0,j,k+1), v(0,j,k), imask3d(1,j,  km1)*imask3d(1,j+1,km1), &
                                                              imask3d(1,j,  km0)*imask3d(1,j+1,km0), sf_u, sf_l)*idz1(k)
             dwdy(0,j,k) = masked_shear(w(0,j+1,k), w(0,j,k), imask3d(1,j+1,km0)*imask3d(1,j+1,km1), &
                                                              imask3d(1,j,  km0)*imask3d(1,j,  km1), sf_s, sf_s)*idy1(0,j)
          END DO


          IF (k==-1) CYCLE

          DO j=-1, jsize+1
             dudy(0,j,k) = masked_shear(u(0,j+1,k), u(0,j,k), imask3d(0,j+1,km0), &
                                                              imask3d(0,j,  km0), sf_s, sf_s) * 0.5*(idy1(0,j)+idy1(1,j))   &
                               - 0.5*(v(0,j,k)+v(1,j,k))*metyx(0,j)

             dvdx(0,j,k) = (v(1,j,k) - v(0,j,k)) * imask3d(1,j,km0)*imask3d(1,j+1,km0)        * 0.5*(idx1(0,j)+idx1(0,j+1)) &
                               - 0.5*(u(0,j,k)+u(0,j+1,k))*metxy(0,j)

             dudy(-1,j,k) = 0.0
             dvdx(-1,j,k) = 0.0
          END DO

          DO j=0, jsize+1
             dudx(0,j,k) = 0.0
             dvdy(0,j,k) = 0.0
             dwdz(0,j,k) = 0.0
          END DO
       END DO
    END IF

    IF (open_n .AND. jcoord==jpes-1) THEN
!$OMP DO
       DO k=-1, ksize+1
          km0=maskindex(k)
          km1=maskindex(k+1)

          IF (kcoord==kpes-1 .AND. k==ksize) THEN
             sf_u = slipfactor_surface
          ELSE
             sf_u = slipfactor_top
          END IF

          DO i=-1, isize+1
             dudz(i,jsize+1,k) = masked_shear(u(i,jsize+1,k+1), u(i,jsize+1,k), imask3d(i,  jsize,km1)*imask3d(i+1,jsize,km1), &
                                                                                imask3d(i,  jsize,km0)*imask3d(i+1,jsize,km0), sf_u, sf_l)*idz1(k)
             dwdx(i,jsize+1,k) = masked_shear(w(i+1,jsize+1,k), w(i,jsize+1,k), imask3d(i+1,jsize,km0)*imask3d(i+1,jsize,km1), &
                                                                                imask3d(i,  jsize,km0)*imask3d(i,  jsize,km1), sf_s, sf_s)*idx1(i,jsize+1)
          END DO

          DO i= 0, isize+1
             dvdz(i,jsize,k) = masked_shear(v(i,jsize,k+1), v(i,jsize,k), imask3d(i,jsize,km1), imask3d(i,jsize,km0), sf_u, sf_l)*idz1(k)
             dwdy(i,jsize,k) =(w(i,jsize+1,k) - w(i,jsize,k)) * imask3d(i,jsize,km0)*imask3d(i,jsize,km1) * idy1(i,jsize)

             dvdz(i,jsize+1,k) = 0.0
             dwdy(i,jsize+1,k) = 0.0
          END DO

          IF (k==-1) CYCLE

          DO i=-1, isize+1
             dudy(i,jsize,k) = (u(i,jsize+1,k)-u(i,jsize,k)) * imask3d(i,jsize,km0)*imask3d(i+1,jsize,km0)    * 0.5*(idy1(i,jsize)+idy1(i+1,jsize)) &
                               - 0.5*(v(i,jsize,k)+v(i+1,jsize,k))*metyx(i,jsize)

             dvdx(i,jsize,k) = masked_shear(v(i+1,jsize,k), v(i,jsize,k), imask3d(i+1,jsize,km0), &
                                                                          imask3d(i,  jsize,km0), sf_s, sf_s) * 0.5*(idx1(i,jsize)+idx1(i,jsize+1)) &
                               - 0.5*(u(i,jsize,k)+u(i,jsize+1,k))*metxy(i,jsize)


             dudy(i,jsize+1,k) = 0.0
             dvdx(i,jsize+1,k) = 0.0
          END DO

          DO i=0, isize+1
             dudx(i,jsize+1,k) = 0.0
             dvdy(i,jsize+1,k) = 0.0
             dwdz(i,jsize+1,k) = 0.0
          END DO
       END DO
    END IF

    IF (open_s .AND. jcoord==0) THEN
!$OMP DO
       DO k=-1, ksize+1
          km0=maskindex(k)
          km1=maskindex(k+1)

          IF (kcoord==kpes-1 .AND. k==ksize) THEN
             sf_u = slipfactor_surface
          ELSE
             sf_u = slipfactor_top
          END IF

          DO i=-1, isize+1
             dudz(i,0,k) = masked_shear(u(i,0,k+1), u(i,0,k), imask3d(i,  1,km1)*imask3d(i+1,1,km1), &
                                                              imask3d(i,  1,km0)*imask3d(i+1,1,km0), sf_u, sf_l)*idz1(k)
             dwdx(i,0,k) = masked_shear(w(i+1,0,k), w(i,0,k), imask3d(i+1,1,km0)*imask3d(i+1,1,km1), &
                                                              imask3d(i,  1,km0)*imask3d(i,  1,km1), sf_s, sf_s)*idx1(i,0)
          END DO

          DO i= 0, isize+1
             dvdz(i,0,k) = masked_shear(v(i,0,k+1), v(i,0,k), imask3d(i,1,km1), imask3d(i,1,km0), sf_u, sf_l)*idz1(k)
             dwdy(i,0,k) =(w(i,1,k) - w(i,0,k)) * imask3d(i,1,km0)*imask3d(i,1,km1) * idy1(i,0)

             dvdz(i,-1,k) = 0.0
             dwdy(i,-1,k) = 0.0
          END DO

          IF (k==-1) CYCLE

          DO i=-1, isize+1
             dudy(i,0,k) = (u(i,1,k)-u(i,0,k)) * imask3d(i,1,km0)*imask3d(i+1,1,km0)          * 0.5*(idy1(i,0)+idy1(i+1,0)) &
                           - 0.5*(v(i,0,k)+v(i+1,0,k))*metyx(i,0)

             dvdx(i,0,k) = masked_shear(v(i+1,0,k), v(i,0,k), imask3d(i+1,1,km0), &
                                                              imask3d(i,  1,km0), sf_s, sf_s) * 0.5*(idx1(i,0)+idx1(i,1)) &
                           - 0.5*(u(i,0,k)+u(i,0+1,k))*metxy(i,0)

             dudy(i,-1,k) = 0.0
             dvdx(i,-1,k) = 0.0
          END DO

          DO i=0, isize+1
             dudx(i,0,k) = 0.0
             dvdy(i,0,k) = 0.0
             dwdz(i,0,k) = 0.0
          END DO
       END DO
    END IF
!$OMP END PARALLEL

  CONTAINS
!---shoud be inlined---
    REAL(8) PURE FUNCTION masked_shear(v1, v0, m1, m0, s1, s0)
      REAL(8),    INTENT(IN) :: v1, v0
      INTEGER(1), INTENT(IN) :: m1, m0
      REAL(4),    INTENT(IN) :: s1, s0

      masked_shear = m0*m1*(v1-v0) + m1*(1-m0)*v1*(1.0-s0) + m0*(1-m1)*v0*(s1-1.0)

    END FUNCTION masked_shear
  END SUBROUTINE velocity_gradient

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE nonlinear(gx, gy, gz)
    USE parameters, ONLY: nonlinear_scheme

    REAL(8), INTENT(INOUT) :: gx(0:isize, 1:jsize, 1:ksize)
    REAL(8), INTENT(INOUT) :: gy(1:isize, 0:jsize, 1:ksize)
    REAL(8), INTENT(INOUT) :: gz(1:isize, 1:jsize, 0:ksize)

    SELECT CASE(nonlinear_scheme)
    CASE (0)
       ! no nonlinear term
    CASE (1)
       CALL nonlinear_o2(gx, gy, gz)
    CASE (2)
       CALL nonlinear_o4(gx, gy, gz)
    CASE (3)
       CALL nonlinear_vi(gx, gy, gz)
    END SELECT

  END SUBROUTINE nonlinear

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE nonlinear_o2(gx, gy, gz)
    REAL(8), INTENT(INOUT) :: gx(0:isize, 1:jsize, 1:ksize)
    REAL(8), INTENT(INOUT) :: gy(1:isize, 0:jsize, 1:ksize)
    REAL(8), INTENT(INOUT) :: gz(1:isize, 1:jsize, 0:ksize)

    REAL(8) :: fx(0:isize+1, 0:jsize,   0:ksize)
    REAL(8) :: fy(0:isize,   0:jsize+1, 0:ksize)
    REAL(8) :: fz(0:isize,   0:jsize,   0:ksize+1)

    INTEGER :: i,j,k

    DO k=1, ksize
    DO j=1, jsize
    DO i=0, isize+1
       fx(i,j,k) = imask3d(i-1,j,k)*imask3d(i,j,k)*imask3d(i+1,j,k) * ((u(i-1,j,k)+u(i,j,k))*0.5)**2
    END DO
    END DO
    END DO

    DO k=1, ksize
    DO j=0, jsize
    DO i=0, isize
       fy(i,j,k) = imask3d(i,j,k)*imask3d(i+1,j,k)*imask3d(i,j+1,k)*imask3d(i+1,j+1,k) * (v(i,j,k)+v(i+1,j,k))*0.5 &
            * (u(i,j,k)*(dy(i,j+1)+dy(i+1,j+1)) + u(i,j+1,k)*(dy(i,j)+dy(i+1,j))) &
            / (dy(i,j)+dy(i+1,j)+dy(i,j+1)+dy(i+1,j+1))
    END DO
    END DO
    END DO

    DO k=0, ksize
    DO j=1, jsize
    DO i=0, isize
       fz(i,j,k) = imask3d(i,j,k)*imask3d(i+1,j,k)*imask3d(i,j,k+1)*imask3d(i+1,j,k+1) * (w(i,j,k)+w(i+1,j,k))*0.5 &
            * (u(i,j,k)*dz(k+1) + u(i,j,k+1)*dz(k)) / (dz(k)+dz(k+1))
    END DO
    END DO
    END DO

    DO k=1, ksize
    DO j=1, jsize
    DO i=0, isize
       gx(i,j,k) = gx(i,j,k) + imask3d(i,j,k)*imask3d(i+1,j,k) &
            * ( fx(i,  j,k)*dsx(i,  j,k) &
              - fx(i+1,j,k)*dsx(i+1,j,k) &
              + fy(i,j-1,k)*(dsy(i,j-1,k)+dsy(i+1,j-1,k))*0.5  &
              - fy(i,j,  k)*(dsy(i,j,  k)+dsy(i+1,j,  k))*0.5  &
              + fz(i,j,k-1)*(dsz(i,j)    +dsz(i+1,j)    )*0.5  &
              - fz(i,j,k)  *(dsz(i,j)    +dsz(i+1,j)    )*0.5) &
            * 2.0/(dvol(i,j,k)+dvol(i+1,j,k))
    END DO
    END DO
    END DO

    DO k=1, ksize
    DO j=0, jsize
    DO i=0, isize
       fx(i,j,k) = imask3d(i,j,k)*imask3d(i,j+1,k)*imask3d(i+1,j,k)*imask3d(i+1,j+1,k) * (u(i,j,k)+u(i,j+1,k))*0.5 &
            * (v(i,j,k)*(dx(i+1,j)+dx(i+1,j+1)) + v(i+1,j,k)*(dx(i,j)+dx(i,j+1))) &
            / (dx(i,j)+dx(i+1,j)+dx(i,j+1)+dx(i+1,j+1))
    END DO
    END DO
    END DO

    DO k=1, ksize
    DO j=0, jsize+1
    DO i=1, isize
       fy(i,j,k) = imask3d(i,j-1,k)*imask3d(i,j,k)*imask3d(i,j+1,k) * ((v(i,j-1,k)+v(i,j,k))*0.5)**2
    END DO
    END DO
    END DO

    DO k=0, ksize
    DO j=0, jsize
    DO i=1, isize
       fz(i,j,k) = imask3d(i,j,k)*imask3d(i,j+1,k)*imask3d(i,j,k+1)*imask3d(i,j+1,k+1) * (w(i,j,k)+w(i,j+1,k))*0.5 &
            * (v(i,j,k)*dz(k+1) + v(i,j,k+1)*dz(k)) / (dz(k)+dz(k+1))
    END DO
    END DO
    END DO

    DO k=1, ksize
    DO j=0, jsize
    DO i=1, isize
       gy(i,j,k) = gy(i,j,k) + imask3d(i,j,k)*imask3d(i,j+1,k) &
            * ( fy(i,j,  k)*dsy(i,j,  k) &
              - fy(i,j+1,k)*dsy(i,j+1,k) &
              + fz(i,j,k-1)*(dsz(i,j)    +dsz(i,j+1)    )*0.5  &
              - fz(i,j,k  )*(dsz(i,j)    +dsz(i,j+1)    )*0.5  &
              + fx(i-1,j,k)*(dsx(i-1,j,k)+dsx(i-1,j+1,k))*0.5  &
              - fx(i,  j,k)*(dsx(i,  j,k)+dsx(i,  j+1,k))*0.5) &
            * 2.0/(dvol(i,j,k)+dvol(i,j+1,k))
    END DO
    END DO
    END DO

    DO k=0, ksize
    DO j=1, jsize
    DO i=0, isize
       fx(i,j,k) = imask3d(i,j,k)*imask3d(i,j,k+1)*imask3d(i+1,j,k)*imask3d(i+1,j,k+1) * (u(i,j,k)+u(i,j,k+1))*0.5 &
            * (w(i,j,k)*dx(i+1,j) + w(i+1,j,k)*dx(i,j)) / (dx(i,j)+dx(i+1,j))
    END DO
    END DO
    END DO

    DO k=0, ksize
    DO j=0, jsize
    DO i=1, isize
       fy(i,j,k) = imask3d(i,j,k)*imask3d(i,j,k+1)*imask3d(i,j+1,k)*imask3d(i,j+1,k+1) * (v(i,j,k)+v(i,j,k+1))*0.5 &
            * (w(i,j,k)*dy(i,j+1) + w(i,j+1,k)*dy(i,j)) / (dy(i,j)+dy(i,j+1))
    END DO
    END DO
    END DO

    DO k=0, ksize+1
    DO j=1, jsize
    DO i=1, isize
       fz(i,j,k) = imask3d(i,j,k-1)*imask3d(i,j,k)*imask3d(i,j,k+1) * ((w(i,j,k-1)+w(i,j,k))*0.5)**2
    END DO
    END DO
    END DO

    DO k=0, ksize
    DO j=1, jsize
    DO i=1, isize
       gz(i,j,k) = gz(i,j,k) + imask3d(i,j,k)*imask3d(i,j,k+1) &
            * ( fz(i,j,k  )*dsz(i,j) &
              - fz(i,j,k+1)*dsz(i,j) &
              + fx(i-1,j,k)*(dsx(i-1,j,k)+dsx(i-1,j,k+1))*0.5  &
              - fx(i,  j,k)*(dsx(i,  j,k)+dsx(i,  j,k+1))*0.5  &
              + fy(i,j-1,k)*(dsy(i,j-1,k)+dsy(i,j-1,k+1))*0.5  &
              - fy(i,j,  k)*(dsy(i,j,  k)+dsy(i,j,  k+1))*0.5) &
            * 2.0/(dvol(i,j,k)+dvol(i,j,k+1))
    END DO
    END DO
    END DO

    CALL metric_term(gx, gy, gz)

  END SUBROUTINE nonlinear_o2

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE nonlinear_o4_nottuned(gx, gy, gz)
    REAL(8), INTENT(INOUT) :: gx(0:isize, 1:jsize, 1:dimz)
    REAL(8), INTENT(INOUT) :: gy(1:isize, 0:jsize, 1:dimz)
    REAL(8), INTENT(INOUT) :: gz(1:isize, 1:jsize, 0:dimz)

    REAL(8) :: fx(0:isize+1, 0:jsize,   0:dimz)
    REAL(8) :: fy(0:isize,   0:jsize+1, 0:dimz)
    REAL(8) :: fz(0:isize,   0:jsize,   0:dimz+1)

    REAL(8) :: d1, d2, d3, d4
    REAL(8) :: a1, a2, a3, a4

    INTEGER :: i,j,k

    DO k=1, ksize
    DO j=1, jsize
    DO i=0, isize+1
       d1 =  dx(i,j)*0.5 + dx(i+1,j)
       d2 =  dx(i,j)*0.5
       d3 = -dx(i,j)*0.5
       d4 = -dx(i,j)*0.5 - dx(i-1,j)

       a1 =  d2*d3*d4*(d2*d3*(d3-d2) + d3*d4*(d4-d3) + d4*d2*(d2-d4))
       a2 = -d1*d3*d4*(d1*d3*(d3-d1) + d3*d4*(d4-d3) + d4*d1*(d1-d4))
       a3 =  d1*d2*d4*(d1*d2*(d2-d1) + d2*d4*(d4-d2) + d4*d1*(d1-d4))
       a4 = -d1*d2*d3*(d1*d2*(d2-d1) + d2*d3*(d3-d2) + d3*d1*(d1-d3))

       fx(i,j,k) = imask3d(i-1,j,k)*imask3d(i,j,k)*imask3d(i+1,j,k) * (u(i-1,j,k)+u(i,j,k))*0.5 &
            * (a1*u(i+1,j,k) + a2*u(i,j,k) + a3*u(i-1,j,k) + a4*u(i-2,j,k)) / (a1+a2+a3+a4)
    END DO
    END DO
    END DO

    DO k=1, ksize
    DO j=0, jsize
    DO i=0, isize
       d1 =  (dy(i,j+2)+dy(i+1,j+2))*0.25 + (dy(i,j+1)+dy(i+1,j+1))*0.5
       d2 =  (dy(i,j+1)+dy(i+1,j+1))*0.25
       d3 = -(dy(i,j)  +dy(i+1,j  ))*0.25
       d4 = -(dy(i,j-1)+dy(i+1,j-1))*0.25 - (dy(i,j)  +dy(i+1,j)  )*0.5

       a1 =  d2*d3*d4*(d2*d3*(d3-d2) + d3*d4*(d4-d3) + d4*d2*(d2-d4))
       a2 = -d1*d3*d4*(d1*d3*(d3-d1) + d3*d4*(d4-d3) + d4*d1*(d1-d4))
       a3 =  d1*d2*d4*(d1*d2*(d2-d1) + d2*d4*(d4-d2) + d4*d1*(d1-d4))
       a4 = -d1*d2*d3*(d1*d2*(d2-d1) + d2*d3*(d3-d2) + d3*d1*(d1-d3))

       fy(i,j,k) = imask3d(i,j+1,k)*imask3d(i+1,j+1,k)*imask3d(i,j,k)*imask3d(i+1,j,k) * (v(i,j,k)+v(i+1,j,k))*0.5 &
            * (a1*u(i,j+2,k) + a2*u(i,j+1,k) + a3*u(i,j,k) + a4*u(i,j-1,k)) / (a1+a2+a3+a4)
    END DO
    END DO
    END DO

    DO k=0, ksize
    DO j=1, jsize
    DO i=0, isize
       d1 =  dz(k+2)*0.5 + dz(k+1)
       d2 =  dz(k+1)*0.5
       d3 = -dz(k  )*0.5
       d4 = -dz(k-1)*0.5 - dz(k)

       a1 =  d2*d3*d4*(d2*d3*(d3-d2) + d3*d4*(d4-d3) + d4*d2*(d2-d4))
       a2 = -d1*d3*d4*(d1*d3*(d3-d1) + d3*d4*(d4-d3) + d4*d1*(d1-d4))
       a3 =  d1*d2*d4*(d1*d2*(d2-d1) + d2*d4*(d4-d2) + d4*d1*(d1-d4))
       a4 = -d1*d2*d3*(d1*d2*(d2-d1) + d2*d3*(d3-d2) + d3*d1*(d1-d3))

       fz(i,j,k) = imask3d(i,j,k+1)*imask3d(i+1,j,k+1)*imask3d(i,j,k)*imask3d(i+1,j,k) * (w(i,j,k)+w(i+1,j,k))*0.5 &
            * (a1*u(i,j,k+2) + a2*u(i,j,k+1) + a3*u(i,j,k) + a4*u(i,j,k-1)) / (a1+a2+a3+a4)
    END DO
    END DO
    END DO

    DO k=1, ksize
    DO j=1, jsize
    DO i=0, isize
       gx(i,j,k) = gx(i,j,k) + imask3d(i,j,k)*imask3d(i+1,j,k) &
            * ( fx(i,  j,k)*dsx(i,  j,k) &
              - fx(i+1,j,k)*dsx(i+1,j,k) &
              + fy(i,j-1,k)*(dsy(i,j-1,k)+dsy(i+1,j-1,k))*0.5  &
              - fy(i,j,  k)*(dsy(i,j,  k)+dsy(i+1,j,  k))*0.5  &
              + fz(i,j,k-1)*(dsz(i,j)    +dsz(i+1,j)    )*0.5  &
              - fz(i,j,k)  *(dsz(i,j)    +dsz(i+1,j)    )*0.5) &
            * 2.0/(dvol(i,j,k)+dvol(i+1,j,k))
    END DO
    END DO
    END DO

    DO k=1, ksize
    DO j=0, jsize
    DO i=0, isize
       d1 =  (dx(i+2,j)+dx(i+2,j+1))*0.25 + (dx(i+1,j)+dx(i+1,j+1))*0.5
       d2 =  (dx(i+1,j)+dx(i+1,j+1))*0.25
       d3 = -(dx(i,  j)+dx(i,  j+1))*0.25
       d4 = -(dx(i-1,j)+dx(i-1,j+1))*0.25 - (dx(i,  j)+dx(i,  j+1))*0.5

       a1 =  d2*d3*d4*(d2*d3*(d3-d2) + d3*d4*(d4-d3) + d4*d2*(d2-d4))
       a2 = -d1*d3*d4*(d1*d3*(d3-d1) + d3*d4*(d4-d3) + d4*d1*(d1-d4))
       a3 =  d1*d2*d4*(d1*d2*(d2-d1) + d2*d4*(d4-d2) + d4*d1*(d1-d4))
       a4 = -d1*d2*d3*(d1*d2*(d2-d1) + d2*d3*(d3-d2) + d3*d1*(d1-d3))

       fx(i,j,k) = imask3d(i,j,k)*imask3d(i,j+1,k)*imask3d(i+1,j,k)*imask3d(i+1,j+1,k) * (u(i,j,k)+u(i,j+1,k))*0.5 &
            * (a1*v(i+2,j,k) + a2*v(i+1,j,k) + a3*v(i,j,k) + a4*v(i-1,j,k)) / (a1+a2+a3+a4)
    END DO
    END DO
    END DO

    DO k=1, ksize
    DO j=0, jsize+1
    DO i=1, isize
       d1 =  dy(i,j)*0.5 + dy(i,j+1)
       d2 =  dy(i,j)*0.5
       d3 = -dy(i,j)*0.5
       d4 = -dy(i,j)*0.5 - dy(i,j-1)

       a1 =  d2*d3*d4*(d2*d3*(d3-d2) + d3*d4*(d4-d3) + d4*d2*(d2-d4))
       a2 = -d1*d3*d4*(d1*d3*(d3-d1) + d3*d4*(d4-d3) + d4*d1*(d1-d4))
       a3 =  d1*d2*d4*(d1*d2*(d2-d1) + d2*d4*(d4-d2) + d4*d1*(d1-d4))
       a4 = -d1*d2*d3*(d1*d2*(d2-d1) + d2*d3*(d3-d2) + d3*d1*(d1-d3))

       fy(i,j,k) = imask3d(i,j-1,k)*imask3d(i,j,k)*imask3d(i,j+1,k) * (v(i,j-1,k)+v(i,j,k))*0.5 &
            * (a1*v(i,j+1,k) + a2*v(i,j,k) + a3*v(i,j-1,k) + a4*v(i,j-2,k)) / (a1+a2+a3+a4)
    END DO
    END DO
    END DO

    DO k=0, ksize
    DO j=0, jsize
    DO i=1, isize
       d1 =  dz(k+2)*0.5 + dz(k+1)
       d2 =  dz(k+1)*0.5
       d3 = -dz(k  )*0.5
       d4 = -dz(k-1)*0.5 - dz(k)

       a1 =  d2*d3*d4*(d2*d3*(d3-d2) + d3*d4*(d4-d3) + d4*d2*(d2-d4))
       a2 = -d1*d3*d4*(d1*d3*(d3-d1) + d3*d4*(d4-d3) + d4*d1*(d1-d4))
       a3 =  d1*d2*d4*(d1*d2*(d2-d1) + d2*d4*(d4-d2) + d4*d1*(d1-d4))
       a4 = -d1*d2*d3*(d1*d2*(d2-d1) + d2*d3*(d3-d2) + d3*d1*(d1-d3))

       fz(i,j,k) = imask3d(i,j,k)*imask3d(i,j+1,k)*imask3d(i,j,k+1)*imask3d(i,j+1,k+1) * (w(i,j,k)+w(i,j+1,k))*0.5 &
            * (a1*v(i,j,k+2) + a2*v(i,j,k+1) + a3*v(i,j,k) + a4*v(i,j,k-1)) / (a1+a2+a3+a4)
    END DO
    END DO
    END DO

    DO k=1, ksize
    DO j=0, jsize
    DO i=1, isize
       gy(i,j,k) = gy(i,j,k) + imask3d(i,j,k)*imask3d(i,j+1,k) &
            * ( fy(i,j,  k)*dsy(i,j,  k) &
              - fy(i,j+1,k)*dsy(i,j+1,k) &
              + fz(i,j,k-1)*(dsz(i,j)    +dsz(i,j+1)    )*0.5  &
              - fz(i,j,k  )*(dsz(i,j)    +dsz(i,j+1)    )*0.5  &
              + fx(i-1,j,k)*(dsx(i-1,j,k)+dsx(i-1,j+1,k))*0.5  &
              - fx(i,  j,k)*(dsx(i,  j,k)+dsx(i,  j+1,k))*0.5) &
            * 2.0/(dvol(i,j,k)+dvol(i,j+1,k))
    END DO
    END DO
    END DO

    DO k=0, ksize
    DO j=1, jsize
    DO i=0, isize
       d1 =  dx(i+2,j)*0.5 + dx(i+1,j)
       d2 =  dx(i+1,j)*0.5
       d3 = -dx(i,  j)*0.5
       d4 = -dx(i-1,j)*0.5 - dx(i,  j)

       a1 =  d2*d3*d4*(d2*d3*(d3-d2) + d3*d4*(d4-d3) + d4*d2*(d2-d4))
       a2 = -d1*d3*d4*(d1*d3*(d3-d1) + d3*d4*(d4-d3) + d4*d1*(d1-d4))
       a3 =  d1*d2*d4*(d1*d2*(d2-d1) + d2*d4*(d4-d2) + d4*d1*(d1-d4))
       a4 = -d1*d2*d3*(d1*d2*(d2-d1) + d2*d3*(d3-d2) + d3*d1*(d1-d3))

       fx(i,j,k) = imask3d(i,j,k)*imask3d(i,j,k+1)*imask3d(i+1,j,k)*imask3d(i+1,j,k+1) * (u(i,j,k)+u(i,j,k+1))*0.5 &
            * (a1*w(i+2,j,k) + a2*w(i+1,j,k) + a3*w(i,j,k) + a4*w(i-1,j,k)) / (a1+a2+a3+a4)
    END DO
    END DO
    END DO

    DO k=0, ksize
    DO j=0, jsize
    DO i=1, isize
       d1 =  dy(i,j+2)*0.5 + dy(i,j+1)
       d2 =  dy(i,j+1)*0.5
       d3 = -dy(i,j)  *0.5
       d4 = -dy(i,j-1)*0.5 - dy(i,j)

       a1 =  d2*d3*d4*(d2*d3*(d3-d2) + d3*d4*(d4-d3) + d4*d2*(d2-d4))
       a2 = -d1*d3*d4*(d1*d3*(d3-d1) + d3*d4*(d4-d3) + d4*d1*(d1-d4))
       a3 =  d1*d2*d4*(d1*d2*(d2-d1) + d2*d4*(d4-d2) + d4*d1*(d1-d4))
       a4 = -d1*d2*d3*(d1*d2*(d2-d1) + d2*d3*(d3-d2) + d3*d1*(d1-d3))

       fy(i,j,k) = imask3d(i,j,k)*imask3d(i,j,k+1)*imask3d(i,j+1,k)*imask3d(i,j+1,k+1) * (v(i,j,k)+v(i,j,k+1))*0.5 &
            * (a1*w(i,j+2,k) + a2*w(i,j+1,k) + a3*w(i,j,k) + a4*w(i,j-1,k)) / (a1+a2+a3+a4)
    END DO
    END DO
    END DO

    DO k=0, ksize+1
    DO j=1, jsize
    DO i=1, isize
       d1 =  dz(k)*0.5 + dz(k+1)
       d2 =  dz(k)*0.5
       d3 = -dz(k)*0.5
       d4 = -dz(k)*0.5 - dz(k-1)

       a1 =  d2*d3*d4*(d2*d3*(d3-d2) + d3*d4*(d4-d3) + d4*d2*(d2-d4))
       a2 = -d1*d3*d4*(d1*d3*(d3-d1) + d3*d4*(d4-d3) + d4*d1*(d1-d4))
       a3 =  d1*d2*d4*(d1*d2*(d2-d1) + d2*d4*(d4-d2) + d4*d1*(d1-d4))
       a4 = -d1*d2*d3*(d1*d2*(d2-d1) + d2*d3*(d3-d2) + d3*d1*(d1-d3))

       fz(i,j,k) = imask3d(i,j,k-1)*imask3d(i,j,k)*imask3d(i,j,k+1) * (w(i,j,k-1)+w(i,j,k))*0.5 &
            * (a1*w(i,j,k+1) + a2*w(i,j,k) + a3*w(i,j,k-1) + a4*w(i,j,k-2)) / (a1+a2+a3+a4)
    END DO
    END DO
    END DO

    DO k=0, ksize
    DO j=1, jsize
    DO i=1, isize
       gz(i,j,k) = gz(i,j,k) + imask3d(i,j,k)*imask3d(i,j,k+1) &
            * ( fz(i,j,k  )*dsz(i,j) &
              - fz(i,j,k+1)*dsz(i,j) &
              + fx(i-1,j,k)*(dsx(i-1,j,k)+dsx(i-1,j,k+1))*0.5  &
              - fx(i,  j,k)*(dsx(i,  j,k)+dsx(i,  j,k+1))*0.5  &
              + fy(i,j-1,k)*(dsy(i,j-1,k)+dsy(i,j-1,k+1))*0.5  &
              - fy(i,j,  k)*(dsy(i,j,  k)+dsy(i,j,  k+1))*0.5) &
            * 2.0/(dvol(i,j,k)+dvol(i,j,k+1))
    END DO
    END DO
    END DO

    CALL metric_term(gx, gy, gz)

  END SUBROUTINE nonlinear_o4_nottuned

!-----------------------------------------------------------------------------------------------------------------------

#define COEFFICIENTS_INIT_ONCE
  SUBROUTINE nonlinear_o4(gx, gy, gz)
    REAL(8), INTENT(INOUT) :: gx(0:isize, 1:jsize, 1:ksize)
    REAL(8), INTENT(INOUT) :: gy(1:isize, 0:jsize, 1:ksize)
    REAL(8), INTENT(INOUT) :: gz(1:isize, 1:jsize, 0:ksize)

    REAL(8) :: fx(0:isize+1, 0:jsize)
    REAL(8) :: fy(0:isize,   0:jsize+1)
    REAL(8) :: fz(0:1, 0:isize, 0:jsize)

#ifdef COEFFICIENTS_INIT_ONCE
    LOGICAL, SAVE :: initialized = .FALSE.

    REAL(8), ALLOCATABLE, SAVE :: axx(:,:,:)
    REAL(8), ALLOCATABLE, SAVE :: axy(:,:,:)
    REAL(8), ALLOCATABLE, SAVE :: axz(:,:)

    REAL(8), ALLOCATABLE, SAVE :: ayx(:,:,:)
    REAL(8), ALLOCATABLE, SAVE :: ayy(:,:,:)
    REAL(8), ALLOCATABLE, SAVE :: ayz(:,:)

    REAL(8), ALLOCATABLE, SAVE :: azx(:,:,:)
    REAL(8), ALLOCATABLE, SAVE :: azy(:,:,:)
    REAL(8), ALLOCATABLE, SAVE :: azz(:,:)
#else
    REAL(8) :: axx(0:3, 0:isize+1, 1:jsize)
    REAL(8) :: axy(0:3, 0:isize,   0:jsize)
    REAL(8) :: axz(0:3, 0:ksize)

    REAL(8) :: ayx(0:3, 0:isize, 0:jsize)
    REAL(8) :: ayy(0:3, 1:isize, 0:jsize+1)
    REAL(8) :: ayz(0:3, 0:ksize)

    REAL(8) :: azx(0:3, 0:isize, 1:jsize)
    REAL(8) :: azy(0:3, 1:isize, 0:jsize)
    REAL(8) :: azz(0:3, 0:ksize+1)
#endif

    INTEGER :: i, j, k
    INTEGER :: kstr, kend
    INTEGER :: k0, k1, kt
    INTEGER :: kl0, kl1
    INTEGER :: km0, km1, km2
    INTEGER :: tid

#ifdef COEFFICIENTS_INIT_ONCE
    IF (.NOT. initialized) THEN
       initialized = .TRUE.

       ALLOCATE(axx(0:3, 0:isize+1, 1:jsize))
       ALLOCATE(axy(0:3, 0:isize,   0:jsize))
       ALLOCATE(axz(0:3, 0:ksize))

       ALLOCATE(ayx(0:3, 0:isize, 0:jsize))
       ALLOCATE(ayy(0:3, 1:isize, 0:jsize+1))
       ALLOCATE(ayz(0:3, 0:ksize))

       ALLOCATE(azx(0:3, 0:isize, 1:jsize))
       ALLOCATE(azy(0:3, 1:isize, 0:jsize))
       ALLOCATE(azz(0:3, 0:ksize+1))
#endif

       DO j=1, jsize
       DO i=0, isize+1
          CALL o4_coefficients(-dx(i,j)*0.5 - dx(i-1,j), &
                               -dx(i,j)*0.5,             &
                                dx(i,j)*0.5,             &
                                dx(i,j)*0.5 + dx(i+1,j), &
                               axx(0,i,j), axx(1,i,j), axx(2,i,j), axx(3,i,j))
       END DO
       END DO

       DO j=0, jsize
       DO i=0, isize
          CALL o4_coefficients(-(dy(i,j-1)+dy(i+1,j-1))*0.25 - (dy(i,j)  +dy(i+1,j)  )*0.5, &
                               -(dy(i,j)  +dy(i+1,j)  )*0.25, &
                                (dy(i,j+1)+dy(i+1,j+1))*0.25, &
                                (dy(i,j+2)+dy(i+1,j+2))*0.25 + (dy(i,j+1)+dy(i+1,j+1))*0.5, &
                               axy(0,i,j), axy(1,i,j), axy(2,i,j), axy(3,i,j))
       END DO
       END DO

       DO k=0, ksize
          CALL o4_coefficients(-dz(k-1)*0.5 - dz(k),   &
                               -dz(k)  *0.5,           &
                                dz(k+1)*0.5,           &
                                dz(k+2)*0.5 + dz(k+1), &
                               axz(0,k), axz(1,k), axz(2,k), axz(3,k))
       END DO

       DO j=0, jsize
       DO i=0, isize
          CALL o4_coefficients(-(dx(i-1,j)+dx(i-1,j+1))*0.25 - (dx(i,  j)+dx(i,  j+1))*0.5, &
                               -(dx(i,  j)+dx(i,  j+1))*0.25,                                 &
                                (dx(i+1,j)+dx(i+1,j+1))*0.25,                                 &
                                (dx(i+2,j)+dx(i+2,j+1))*0.25 + (dx(i+1,j)+dx(i+1,j+1))*0.5, &
                               ayx(0,i,j), ayx(1,i,j), ayx(2,i,j), ayx(3,i,j))
       END DO
       END DO

       DO j=0, jsize+1
       DO i=1, isize
          CALL o4_coefficients(-dy(i,j)*0.5 - dy(i,j-1), &
                               -dy(i,j)*0.5,             &
                                dy(i,j)*0.5,             &
                                dy(i,j)*0.5 + dy(i,j+1), &
                               ayy(0,i,j), ayy(1,i,j), ayy(2,i,j), ayy(3,i,j))
       END DO
       END DO

       DO k=0, ksize
          CALL o4_coefficients(-dz(k-1)*0.5 - dz(k),   &
                               -dz(k)  *0.5,           &
                                dz(k+1)*0.5,           &
                                dz(k+2)*0.5 + dz(k+1), &
                               ayz(0,k), ayz(1,k), ayz(2,k), ayz(3,k))
       END DO

       DO j=1, jsize
       DO i=0, isize
          CALL o4_coefficients(-dx(i-1,j)*0.5 - dx(i,  j), &
                               -dx(i,  j)*0.5,             &
                                dx(i+1,j)*0.5,             &
                                dx(i+2,j)*0.5 + dx(i+1,j), &
                               azx(0,i,j), azx(1,i,j), azx(2,i,j), azx(3,i,j))
       END DO
       END DO

       DO j=0, jsize
       DO i=1, isize
          CALL o4_coefficients(-dy(i,j-1)*0.5 - dy(i,j),   &
                               -dy(i,j)  *0.5,             &
                                dy(i,j+1)*0.5,             &
                                dy(i,j+2)*0.5 + dy(i,j+1), &
                               azy(0,i,j), azy(1,i,j), azy(2,i,j), azy(3,i,j))
       END DO
       END DO

       DO k=0, ksize+1
          CALL o4_coefficients(-dz(k)*0.5 - dz(k-1), &
                               -dz(k)*0.5,           &
                                dz(k)*0.5,           &
                                dz(k)*0.5 + dz(k+1), &
                               azz(0,k), azz(1,k), azz(2,k), azz(3,k))
       END DO

#ifdef COEFFICIENTS_INIT_ONCE
    END IF
#endif

    tid  = 0
    kstr = 1
    kend = ksize

!$OMP PARALLEL PRIVATE(i, j, k, tid) &
!$OMP          PRIVATE(fx, fy, fz) &
!$OMP          PRIVATE(kstr, kend) &
!$OMP          PRIVATE(k0, k1, kt) &
!$OMP          PRIVATE(kl0, kl1, km0, km1, km2)

!$  tid  = omp_get_thread_num()
!$  kstr = kstr_t(tid)
!$  kend = kend_t(tid)

    k0=0
    k1=1

    DO k=kstr-1, kend
       kt = k0
       k0 = k1
       k1 = kt

       kl0=dzindex(k)

       km0=maskindex(k)
       km1=maskindex(k+1)

       DO j=1, jsize
       DO i=0, isize
          fz(k1,i,j) = imask3d(i,j,km0)*imask3d(i+1,j,km0)*imask3d(i,j,km1)*imask3d(i+1,j,km1) * (w(i,j,k)+w(i+1,j,k))*0.5 &
               * (axz(3,k)*u(i,j,k+2) + axz(2,k)*u(i,j,k+1) + axz(1,k)*u(i,j,k) + axz(0,k)*u(i,j,k-1)) &
               * (dsz(i,j)+dsz(i+1,j))*0.5
       END DO
       END DO

       IF (k==kstr-1) CYCLE

       DO j=1, jsize
       DO i=0, isize+1
          fx(i,j) = imask3d(i-1,j,km0)*imask3d(i,j,km0)*imask3d(i+1,j,km0) * (u(i-1,j,k)+u(i,j,k))*0.5 &
               * (axx(3,i,j)*u(i+1,j,k) + axx(2,i,j)*u(i,j,k) + axx(1,i,j)*u(i-1,j,k) + axx(0,i,j)*u(i-2,j,k)) &
               * dsx(i,j,kl0)
       END DO
       END DO

       DO j=0, jsize
       DO i=0, isize
          fy(i,j) = imask3d(i,j,km0)*imask3d(i+1,j,km0)*imask3d(i,j+1,km0)*imask3d(i+1,j+1,km0) * (v(i,j,k)+v(i+1,j,k))*0.5 &
               * (axy(3,i,j)*u(i,j+2,k) + axy(2,i,j)*u(i,j+1,k) + axy(1,i,j)*u(i,j,k) + axy(0,i,j)*u(i,j-1,k))  &
               * (dsy(i,j,kl0)+dsy(i+1,j,kl0))*0.5
       END DO
       END DO

       DO j=1, jsize
       DO i=0, isize
          gx(i,j,k) = gx(i,j,k) + ( fx(i,j)   - fx(i+1,j)   &
                                  + fy(i,j-1) - fy(i,  j)   &
                                  + fz(k0,i,j)- fz(k1,i,j)) * 2.0/(dvol(i,j,kl0)+dvol(i+1,j,kl0))
       END DO
       END DO
    END DO

    DO k=kstr-1, kend
       kt = k0
       k0 = k1
       k1 = kt

       kl0=dzindex(k)
       km0=maskindex(k)
       km1=maskindex(k+1)

       DO j=0, jsize
       DO i=1, isize
          fz(k1,i,j) = imask3d(i,j,km0)*imask3d(i,j+1,km0)*imask3d(i,j,km1)*imask3d(i,j+1,km1) * (w(i,j,k)+w(i,j+1,k))*0.5 &
               * (ayz(3,k)*v(i,j,k+2) + ayz(2,k)*v(i,j,k+1) + ayz(1,k)*v(i,j,k) + ayz(0,k)*v(i,j,k-1)) &
               * (dsz(i,j)+dsz(i,j+1))*0.5
       END DO
       END DO

       IF (k==kstr-1) CYCLE

       DO j=0, jsize
       DO i=0, isize
          fx(i,j) = imask3d(i,j,km0)*imask3d(i,j+1,km0)*imask3d(i+1,j,km0)*imask3d(i+1,j+1,km0) * (u(i,j,k)+u(i,j+1,k))*0.5 &
               * (ayx(3,i,j)*v(i+2,j,k) + ayx(2,i,j)*v(i+1,j,k) + ayx(1,i,j)*v(i,j,k) + ayx(0,i,j)*v(i-1,j,k)) &
               * (dsx(i,j,kl0)+dsx(i,j+1,kl0))*0.5
       END DO
       END DO

       DO j=0, jsize+1
       DO i=1, isize
          fy(i,j) = imask3d(i,j-1,km0)*imask3d(i,j,km0)*imask3d(i,j+1,km0) * (v(i,j-1,k)+v(i,j,k))*0.5 &
               * (ayy(3,i,j)*v(i,j+1,k) + ayy(2,i,j)*v(i,j,k) + ayy(1,i,j)*v(i,j-1,k) + ayy(0,i,j)*v(i,j-2,k)) &
               * dsy(i,j,kl0)
       END DO
       END DO

       DO j=0, jsize
       DO i=1, isize
          gy(i,j,k) = gy(i,j,k) + ( fy(i,j)    - fy(i,j+1)  &
                                  + fz(k0,i,j) - fz(k1,i,j) &
                                  + fx(i-1,j)  - fx(i,j)) * 2.0/(dvol(i,j,kl0)+dvol(i,j+1,kl0))
       END DO
       END DO
    END DO

    DO k=kstr-2, kend
       kt = k0
       k0 = k1
       k1 = kt

       kl0=dzindex(k)
       kl1=dzindex(k+1)

       km0=maskindex(k)
       km1=maskindex(k+1)
       km2=maskindex(k+2)

!$     IF (tid/= 0 .AND. k==kstr-2) CYCLE

       DO j=1, jsize
       DO i=1, isize
             fz(k1,i,j) = imask3d(i,j,km0)*imask3d(i,j,km1)*imask3d(i,j,km2) * (w(i,j,k)+w(i,j,k+1))*0.5 &
                  * (azz(3,k+1)*w(i,j,k+2) + azz(2,k+1)*w(i,j,k+1) + azz(1,k+1)*w(i,j,k) + azz(0,k+1)*w(i,j,k-1)) &
                  * dsz(i,j)
       END DO
       END DO

       IF (k==kstr-2) CYCLE

!$     IF (tid/= 0 .AND. k==kstr-1) CYCLE

       DO j=1, jsize
       DO i=0, isize
          fx(i,j) = imask3d(i,j,km0)*imask3d(i,j,km1)*imask3d(i+1,j,km0)*imask3d(i+1,j,km1) * (u(i,j,k)+u(i,j,k+1))*0.5 &
               * (azx(3,i,j)*w(i+2,j,k) + azx(2,i,j)*w(i+1,j,k) + azx(1,i,j)*w(i,j,k) + azx(0,i,j)*w(i-1,j,k)) &
               * (dsx(i,j,kl0)+dsx(i,j,kl1))*0.5
       END DO
       END DO

       DO j=0, jsize
       DO i=1, isize
          fy(i,j) = imask3d(i,j,km0)*imask3d(i,j,km1)*imask3d(i,j+1,km0)*imask3d(i,j+1,km1) * (v(i,j,k)+v(i,j,k+1))*0.5 &
               * (azy(3,i,j)*w(i,j+2,k) + azy(2,i,j)*w(i,j+1,k) + azy(1,i,j)*w(i,j,k) + azy(0,i,j)*w(i,j-1,k)) &
               * (dsy(i,j,kl0)+dsy(i,j,kl1))*0.5
       END DO
       END DO

       DO j=1, jsize
       DO i=1, isize
          gz(i,j,k) = gz(i,j,k) + ( fz(k0,i,j)  - fz(k1,i,j)  &
                                  + fx(i-1,j)   - fx(i,j)     &
                                  + fy(i,  j-1) - fy(i,j)) * 2.0/(dvol(i,j,kl0)+dvol(i,j,kl1))
       END DO
       END DO
    END DO
!$OMP END PARALLEL

    CALL metric_term(gx, gy, gz)

  CONTAINS
!---shoud be inlined---
    SUBROUTINE o4_coefficients(d0, d1, d2, d3, a0, a1, a2, a3)
      REAL(8), INTENT(IN)  :: d0, d1, d2, d3
      REAL(8), INTENT(OUT) :: a0, a1, a2, a3
      REAL(8) :: a

      a0 = -d3*d2*d1*(d3*d2*(d2-d3) + d2*d1*(d1-d2) + d1*d3*(d3-d1))
      a1 =  d3*d2*d0*(d3*d2*(d2-d3) + d2*d0*(d0-d2) + d0*d3*(d3-d0))
      a2 = -d3*d1*d0*(d3*d1*(d1-d3) + d1*d0*(d0-d1) + d0*d3*(d3-d0))
      a3 =  d2*d1*d0*(d2*d1*(d1-d2) + d1*d0*(d0-d1) + d0*d2*(d2-d0))

      a = a3+a2+a1+a0

      a0 = a0/a
      a1 = a1/a
      a2 = a2/a
      a3 = a3/a
    END SUBROUTINE o4_coefficients

  END SUBROUTINE nonlinear_o4

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE metric_term(gx, gy, gz)
    USE parameters, ONLY: earth_radius

    REAL(8), INTENT(INOUT) :: gx(0:isize, 1:jsize, 1:ksize)
    REAL(8), INTENT(INOUT) :: gy(1:isize, 0:jsize, 1:ksize)
    REAL(8), INTENT(INOUT) :: gz(1:isize, 1:jsize, 0:ksize)

    INTEGER :: i, j, k
    INTEGER :: kl0, kl1
    INTEGER :: km0, km1

!$OMP PARALLEL DO PRIVATE(i, j, k, kl0, km0)
    DO k=1, ksize
       kl0=dzindex(k)
       km0=maskindex(k)

       DO j=1, jsize
       DO i=0, isize
          gx(i,j,k) = gx(i,j,k) - imask3d(i,j,km0)*imask3d(i+1,j,km0) &
               * 0.5 * (u(i,j,k)*((v(i,  j-1,k)*metxy(i,j-1) + v(i,  j,k)*metxy(i,j))*dvol(i,  j,kl0)  &
                                   +(v(i+1,j-1,k)*metxy(i,j-1) + v(i+1,j,k)*metxy(i,j))*dvol(i+1,j,kl0)) &
                         - ((v(i,  j-1,k)**2*metyx(i,j-1) + v(i,  j,k)**2*metyx(i,j))*dvol(i,  j,kl0)    &
                           +(v(i+1,j-1,k)**2*metyx(i,j-1) + v(i+1,j,k)**2*metyx(i,j))*dvol(i+1,j,kl0)))  &
               / (dvol(i,j,kl0)+dvol(i+1,j,kl0))
       END DO
       END DO

       DO j=0, jsize
       DO i=1, isize
          gy(i,j,k) = gy(i,j,k) - imask3d(i,j,km0)*imask3d(i,j+1,km0) &
               * 0.5 * (v(i,j,k)*((u(i-1,j,  k)*metyx(i-1,j) + u(i,j,  k)*metyx(i,j))*dvol(i,j,  kl0)  &
                                   +(u(i-1,j+1,k)*metyx(i-1,j) + u(i,j+1,k)*metyx(i,j))*dvol(i,j+1,kl0)) &
                         - ((u(i-1,j,  k)**2*metxy(i-1,j) + u(i,j,  k)**2*metxy(i,j))*dvol(i,j,  kl0)    &
                           +(u(i-1,j+1,k)**2*metxy(i-1,j) + u(i,j+1,k)**2*metxy(i,j))*dvol(i,j+1,kl0)))  &
               / (dvol(i,j,kl0)+dvol(i,j+1,kl0))
       END DO
       END DO
    END DO

#ifdef METRIC_SPHERIC_EARTH
!$OMP PARALLEL PRIVATE(i, j, k, kl0, km0)
!$OMP DO
    DO k=1, ksize
       kl0=dzindex(k)
       km0=maskindex(k)

       DO j=1, jsize
       DO i=0, isize
          gx(i,j,k) = gx(i,j,k) - imask3d(i,j,km0)*imask3d(i+1,j,km0)                &
               * 0.5 * u(i,j,k)*( (w(i,  j,k-1)+w(i,  j,k))*dvol(i,  j,kl0)  &
                                   +(w(i+1,j,k-1)+w(i+1,j,k))*dvol(i+1,j,kl0)) &
               / (earth_raduis*(dvol(i,j,kl0)+dvol(i+1,j,kl0)))
       END DO
       END DO

       DO j=0, jsize
       DO i=1, isize
          gy(i,j,k) = gy(i,j,k) - imask3d(i,j,km0)*imask3d(i,j+1,km0)                &
               * 0.5 * v(i,j,k)*( (w(i,j,  k-1)+w(i,j,  k))*dvol(i,j,  kl0)  &
                                   +(w(i,j+1,k-1)+w(i,j+1,k))*dvol(i,j+1,kl0)) &
               / (earth_raduis*(dvol(i,j,kl0)+dvol(i,j+1,kl0)))
       END DO
       END DO
    END DO

!$OMP DO
    DO k=0, ksize
       kl0=dzindex(k)
       kl1=dzindex(k+1)

       km0=maskindex(k)
       km1=maskindex(k+1)

       DO j=1, jsize
       DO i=1, isize
          gz(i,j,k) = gz(i,j,k) + imask3d(i,j,km0)*imask3d(i,j,km1) &
               * 0.5 * ( (u(i-1,j,k  )**2+u(i,j,k  )**2+v(i,j-1,k  )**2+v(i,j,k  )**2)*dvol(i,j,kl0)  &
                          +(u(i-1,j,k+1)**2+u(i,j,k+1)**2+v(i,j-1,k+1)**2+v(i,j,k+1)**2)*dvol(i,j,kl1)) &
               / (earth_radius*(dvol(i,j,kl0)+dvol(i,j,kl1)))
       END DO
       END DO
    END DO
!$OMP END PARALLEL
#endif

  END SUBROUTINE metric_term

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE nonlinear_vi(gx, gy, gz)
    REAL(8), INTENT(INOUT) :: gx(0:isize, 1:jsize, 1:ksize)
    REAL(8), INTENT(INOUT) :: gy(1:isize, 0:jsize, 1:ksize)
    REAL(8), INTENT(INOUT) :: gz(1:isize, 1:jsize, 0:ksize)

    INTEGER :: i, j, k
    INTEGER :: kl0, kl1, km0, km1

!$OMP PARALLEL PRIVATE(i, j, k, kl0, kl1, km0, km1)
!$OMP DO
    DO k=0, ksize
       kl0=dzindex(k)
       kl1=dzindex(k+1)
       km0=maskindex(k)
       km1=maskindex(k+1)

       DO j=1, jsize
       DO i=1, isize
          gz(i,j,k) = gz(i,j,k) + imask3d(i,j,km0)*imask3d(i,j,km1)*(2.0*(ke(i,j,k)-ke(i,j,k+1))*dsz(i,j) &
               + 0.5*dvol(i,j,kl0) * ( u(i-1,j,  k)  *(dudz(i-1,j,k)-dwdx(i-1,j,k)) + u(i,j,k)  *(dudz(i,j,k)-dwdx(i,j,k))   &
                                     - v(i,  j-1,k)  *(dwdy(i,j-1,k)-dvdz(i,j-1,k)) - v(i,j,k)  *(dwdy(i,j,k)-dvdz(i,j,k)))  &
               + 0.5*dvol(i,j,kl1) * ( u(i-1,j,  k+1)*(dudz(i-1,j,k)-dwdx(i-1,j,k)) + u(i,j,k+1)*(dudz(i,j,k)-dwdx(i,j,k))   &
                                     - v(i,  j-1,k+1)*(dwdy(i,j-1,k)-dvdz(i,j-1,k)) - v(i,j,k+1)*(dwdy(i,j,k)-dvdz(i,j,k)))) &
               / (dvol(i,j,kl0)+dvol(i,j,kl1))
       END DO
       END DO

       IF (k==0) CYCLE

       DO j=1, jsize
       DO i=0, isize
          gx(i,j,k) = gx(i,j,k) + imask3d(i,j,km0)*imask3d(i+1,j,km0)*(2.0*(ke(i,j,k)-ke(i+1,j,k))*dsx(i,j,kl0) &
               + 0.5*dvol(i,  j,kl0) * ( v(i,  j-1,k)  *(dvdx(i,j-1,k)-dudy(i,j-1,k)) + v(i,  j,k)*(dvdx(i,j,k)-dudy(i,j,k))   &
                                       - w(i,  j,  k-1)*(dudz(i,j,k-1)-dwdx(i,j,k-1)) - w(i,  j,k)*(dudz(i,j,k)-dwdx(i,j,k)))  &
               + 0.5*dvol(i+1,j,kl0) * ( v(i+1,j-1,k)  *(dvdx(i,j-1,k)-dudy(i,j-1,k)) + v(i+1,j,k)*(dvdx(i,j,k)-dudy(i,j,k))   &
                                       - w(i+1,j,  k-1)*(dudz(i,j,k-1)-dwdx(i,j,k-1)) - w(i+1,j,k)*(dudz(i,j,k)-dwdx(i,j,k)))) &
               / (dvol(i,j,kl0)+dvol(i+1,j,kl0))
       END DO
       END DO

       DO j=0, jsize
       DO i=1, isize
          gy(i,j,k) = gy(i,j,k) + imask3d(i,j,km0)*imask3d(i,j+1,km0)*(2.0*(ke(i,j,k)-ke(i,j+1,k))*dsy(i,j,kl0) &
               + 0.5*dvol(i,j,  kl0) * ( w(i,  j,  k-1)*(dwdy(i,j,k-1)-dvdz(i,j,k-1)) + w(i,j,  k)*(dwdy(i,j,k)-dvdz(i,j,k))   &
                                       - u(i-1,j,  k)  *(dvdx(i-1,j,k)-dudy(i-1,j,k)) - u(i,j,  k)*(dvdx(i,j,k)-dudy(i,j,k)))  &
               + 0.5*dvol(i,j+1,kl0) * ( w(i,  j+1,k-1)*(dwdy(i,j,k-1)-dvdz(i,j,k-1)) + w(i,j+1,k)*(dwdy(i,j,k)-dvdz(i,j,k))   &
                                       - u(i-1,j+1,k)  *(dvdx(i-1,j,k)-dudy(i-1,j,k)) - u(i,j+1,k)*(dvdx(i,j,k)-dudy(i,j,k)))) &
               / (dvol(i,j,kl0)+dvol(i,j+1,kl0))
       END DO
       END DO
    END DO
!$OMP END PARALLEL

  END SUBROUTINE nonlinear_vi

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE coriolis(gx, gy, gz)
    REAL(8), INTENT(INOUT) :: gx(0:isize, 1:jsize, 1:ksize)
    REAL(8), INTENT(INOUT) :: gy(1:isize, 0:jsize, 1:ksize)
    REAL(8), INTENT(INOUT) :: gz(1:isize, 1:jsize, 0:ksize)

    LOGICAL :: stat
    INTEGER :: i, j, k
    INTEGER :: kl0, kl1
    INTEGER :: km0, km1

    CALL checkin('CORIOLIS_X', corx)
    CALL checkin('CORIOLIS_Y', cory)
    CALL checkin('CORIOLIS_Z', corz)

!$OMP PARALLEL PRIVATE(i, j, k, kl0, kl1, km0, km1)
!$OMP DO
    DO j=1, jsize
    DO i=1, isize
       gz(i,j,0) = gz(i,j,0) + ( cory(i,j) * (u(i-1,j,0)+u(i,j,0))*0.5 * dvol(i,j,0)  &
                               + cory(i,j) * (u(i-1,j,1)+u(i,j,1))*0.5 * dvol(i,j,1)  &
                               - corx(i,j) * (v(i,j-1,0)+v(i,j,0))*0.5 * dvol(i,j,0)  &
                               - corx(i,j) * (v(i,j-1,1)+v(i,j,1))*0.5 * dvol(i,j,1)) &
                             * imask3d(i,j,0)*imask3d(i,j,1) / (dvol(i,j,0)+dvol(i,j,1))
    END DO
    END DO

!$OMP DO
    DO k=1, ksize
       kl0=dzindex(k)
       kl1=dzindex(k+1)

       km0=maskindex(k)
       km1=maskindex(k+1)

       DO j=1, jsize
       DO i=0, isize
          gx(i,j,k) = gx(i,j,k) + ( corz(i,  j) * (v(i,  j-1,k)+v(i,  j,k))*0.5 * dvol(i,  j,kl0)  &
                                  + corz(i+1,j) * (v(i+1,j-1,k)+v(i+1,j,k))*0.5 * dvol(i+1,j,kl0)  &
                                  - cory(i,  j) * (w(i,  j,k-1)+w(i,  j,k))*0.5 * dvol(i,  j,kl0)  &
                                  - cory(i+1,j) * (w(i+1,j,k-1)+w(i+1,j,k))*0.5 * dvol(i+1,j,kl0)) &
                                * imask3d(i,j,km0)*imask3d(i+1,j,km0) / (dvol(i,j,kl0)+dvol(i+1,j,kl0))
       END DO
       END DO

       DO j=0, jsize
       DO i=1, isize
          gy(i,j,k) = gy(i,j,k) + ( corx(i,j)   * (w(i,j,  k-1)+w(i,j,  k))*0.5 * dvol(i,j,  kl0)  &
                                  + corx(i,j+1) * (w(i,j+1,k-1)+w(i,j+1,k))*0.5 * dvol(i,j+1,kl0)  &
                                  - corz(i,j)   * (u(i-1,j,  k)+u(i,j,  k))*0.5 * dvol(i,j,  kl0)  &
                                  - corz(i,j+1) * (u(i-1,j+1,k)+u(i,j+1,k))*0.5 * dvol(i,j+1,kl0)) &
                                * imask3d(i,j,km0)*imask3d(i,j+1,km0) / (dvol(i,j,kl0)+dvol(i,j+1,kl0))
       END DO
       END DO

       DO j=1, jsize
       DO i=1, isize
          gz(i,j,k) = gz(i,j,k) + ( cory(i,j) * (u(i-1,j,k)  +u(i,j,k)  )*0.5 * dvol(i,j,kl0)  &
                                  + cory(i,j) * (u(i-1,j,k+1)+u(i,j,k+1))*0.5 * dvol(i,j,kl1)  &
                                  - corx(i,j) * (v(i,j-1,k)  +v(i,j,k)  )*0.5 * dvol(i,j,kl0)  &
                                  - corx(i,j) * (v(i,j-1,k+1)+v(i,j,k+1))*0.5 * dvol(i,j,kl1)) &
                                * imask3d(i,j,km0)*imask3d(i,j,km1) / (dvol(i,j,kl0)+dvol(i,j,kl1))
       END DO
       END DO
    END DO
!$OMP END PARALLEL

  END SUBROUTINE coriolis

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE body_forcing(gx, gy, gz)
    REAL(8), INTENT(INOUT) :: gx(0:isize, 1:jsize, 1:ksize)
    REAL(8), INTENT(INOUT) :: gy(1:isize, 0:jsize, 1:ksize)
    REAL(8), INTENT(INOUT) :: gz(1:isize, 1:jsize, 0:ksize)

    INTEGER :: i, j, k
    LOGICAL :: stat

    REAL(8) :: tmp(0:isize, 0:jsize, 0:ksize)

    CALL checkin('BODYFORCE_U', tmp, stat=stat)
    IF (stat) THEN
!$OMP PARALLEL DO
       DO k=1, ksize
       DO j=1, jsize
       DO i=0, isize
          gx(i,j,k) = gx(i,j,k) + tmp(i,j,k)*imask3d(i,j,k)*imask3d(i+1,j,k)
       END DO
       END DO
       END DO
    END IF

    CALL checkin('BODYFORCE_V', tmp, stat=stat)
    IF (stat) THEN
!$OMP PARALLEL DO
       DO k=1, ksize
       DO j=0, jsize
       DO i=1, isize
          gy(i,j,k) = gy(i,j,k) + tmp(i,j,k)*imask3d(i,j,k)*imask3d(i,j+1,k)
       END DO
       END DO
       END DO
    END IF

    CALL checkin('BODYFORCE_W', tmp, stat=stat)
    IF (stat) THEN
!$OMP PARALLEL DO
       DO k=0, ksize
       DO j=1, jsize
       DO i=1, isize
          gz(i,j,k) = gz(i,j,k) + tmp(i,j,k)*imask3d(i,j,k)*imask3d(i,j,k+1)
       END DO
       END DO
       END DO
    END IF

  END SUBROUTINE body_forcing

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE wave_forcing(gx, gy, gz)
    REAL(8), INTENT(INOUT) :: gx(0:isize, 1:jsize, 1:ksize)
    REAL(8), INTENT(INOUT) :: gy(1:isize, 0:jsize, 1:ksize)
    REAL(8), INTENT(INOUT) :: gz(1:isize, 1:jsize, 0:ksize)

    INTEGER :: i, j, k
    INTEGER :: kl0, kl1, km0, km1

    REAL(8) :: u_st(0:ksize)
    REAL(8) :: v_st(0:ksize)

    LOGICAL :: stat_u, stat_v, stat_cor

    CALL checkin('U_STOKES', u_st, axis='Z', stat=stat_u)
    CALL checkin('V_STOKES', v_st, axis='Z', stat=stat_v)

    IF (.NOT. (stat_u .OR. stat_v)) RETURN

    IF (.NOT. stat_u) u_st(:) = 0.0
    IF (.NOT. stat_v) v_st(:) = 0.0


!$OMP PARALLEL PRIVATE(i, j, k, kl0, kl1, km0, km1)
!$OMP DO
    DO k=0, ksize
       kl0=dzindex(k)
       kl1=dzindex(k+1)
       km0=maskindex(k)
       km1=maskindex(k+1)

       DO j=1, jsize
       DO i=1, isize
          gz(i,j,k) = gz(i,j,k) + imask3d(i,j,km0)*imask3d(i,j,km1) &
               * ( u_st(k) *(0.5*(dudz(i-1,j,k)-dwdx(i-1,j,k)+dudz(i,j,k)-dwdx(i,j,k)) + cory(i,j))  &
                 - v_st(k) *(0.5*(dwdy(i,j-1,k)-dvdz(i,j-1,k)+dwdy(i,j,k)-dvdz(i,j,k)) + corx(i,j)))
       END DO
       END DO

       IF (k==0) CYCLE

       DO j=1, jsize
       DO i=0, isize
          gx(i,j,k) = gx(i,j,k) + imask3d(i-1,j,km0)*imask3d(i,j,km0) &
               * 0.5*(v_st(k-1)+v_st(k)) * (dvol(i,  j,kl0)*(0.5*(dvdx(i,j-1,k)-dudy(i,j-1,k)+dvdx(i,j,k)-dudy(i,j,k)) + corz(i,  j))  &
                                           +dvol(i+1,j,kl0)*(0.5*(dvdx(i,j-1,k)-dudy(i,j-1,k)+dvdx(i,j,k)-dudy(i,j,k)) + corz(i+1,j))) &
                                           / (dvol(i,j,kl0)+dvol(i+1,j,kl0))
       END DO
       END DO

       DO j=0, jsize
       DO i=1, isize
          gy(i,j,k) = gy(i,j,k) - imask3d(i,j-1,km0)*imask3d(i,j,km0) &
               * 0.5*(u_st(k-1)+u_st(k)) * (dvol(i,j,  kl0)*(0.5*(dvdx(i-1,j,k)-dudy(i-1,j,k)+dvdx(i,j,k)-dudy(i,j,k)) + corz(i,j  ))  &
                                           +dvol(i,j+1,kl0)*(0.5*(dvdx(i-1,j,k)-dudy(i-1,j,k)+dvdx(i,j,k)-dudy(i,j,k)) + corz(i,j+1))) &
                                           / (dvol(i,j,kl0)+dvol(i,j+1,kl0))

       END DO
       END DO
    END DO
!$OMP END PARALLEL

  END SUBROUTINE wave_forcing

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE grad_slp(gx, gy)
    USE parameters, ONLY: rho_0

    REAL(8), INTENT(INOUT) :: gx(0:isize, 1:jsize, 1:ksize)
    REAL(8), INTENT(INOUT) :: gy(1:isize, 0:jsize, 1:ksize)

    REAL(8) :: slp(0:isize+1, 0:jsize+1)

    LOGICAL :: stat
    INTEGER :: i, j, k, km

    CALL checkin('SLP', slp, stat)
    IF (.NOT. stat) RETURN

    slp = slp(:,:) / rho_0

!$OMP PARALLEL DO PRIVATE(km)
    DO k=1, ksize
       km = maskindex(k)
       DO j=1, jsize
       DO i=0, isize
          gx(i,j,k) = gx(i,j,k) - imask3d(i,j,km)*imask3d(i+1,j,km) * (slp(i+1,j) - slp(i,j))*idx1(i,j)
       END DO
       END DO

       DO j=0, jsize
       DO i=1, isize
          gy(i,j,k) = gy(i,j,k) - imask3d(i,j,km)*imask3d(i,j+1,km) * (slp(i,j+1) - slp(i,j))*idy1(i,j)
       END DO
       END DO
    END DO

  END SUBROUTINE grad_slp

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE grad_ssh(gx, gy)
    USE parameters, ONLY: gravity, pgradient_scheme
    REAL(8), INTENT(INOUT) :: gx(0:isize, 1:jsize, 1:ksize)
    REAL(8), INTENT(INOUT) :: gy(1:isize, 0:jsize, 1:ksize)

    REAL(8) :: tmpx(0:isize, 1:jsize, 2)
    REAL(8) :: tmpy(1:isize, 0:jsize, 2)

    INTEGER :: i, j, k
    INTEGER :: km

    SELECT CASE (pgradient_scheme)
    CASE (1)
!$OMP PARALLEL DO PRIVATE(km)
       DO k=1, ksize
          km = maskindex(k)

          DO j=1, jsize
          DO i=0, isize
             gx(i,j,k) = gx(i,j,k) - imask3d(i,j,km)*imask3d(i+1,j,km) * gravity * (ssh(i+1,j)-ssh(i,j))*idx1(i,j)
          END DO
          END DO

          DO j=0, jsize
          DO i=1, isize
             gy(i,j,k) = gy(i,j,k) - imask3d(i,j,km)*imask3d(i,j+1,km) * gravity * (ssh(i,j+1)-ssh(i,j))*idy1(i,j)
          END DO
          END DO
       END DO

    CASE (2)
!$OMP PARALLEL DO PRIVATE(km)
       DO k=1, ksize
          km = maskindex(k)

          DO j=1, jsize
          DO i=0, isize
             gx(i,j,k) = gx(i,j,k) - imask3d(i,j,km)*imask3d(i+1,j,km) * gravity * (ssh(i+1,j)*dz_star(i+1,j,k)*dy(i+1,j) - ssh(i,j)*dz_star(i,j,k)*dy(i,j)) * 2.0/(dvol(i,j,k)+dvol(i+1,j,k))
          END DO
          END DO

          DO j=0, jsize
          DO i=1, isize
             gy(i,j,k) = gy(i,j,k) - imask3d(i,j,km)*imask3d(i,j+1,km) * gravity * (ssh(i,j+1)*dz_star(i,j+1,k)*dx(i,j+1) - ssh(i,j)*dz_star(i,j,k)*dx(i,j)) * 2.0/(dvol(i,j,k)+dvol(i,j+1,k))
          END DO
          END DO
       END DO


    CASE (3)
       tmpx(:,:,:) = 0.0
       tmpy(:,:,:) = 0.0

       DO k=1, ksize
          km = maskindex(k)

          DO j=1, jsize
          DO i=0, isize
             tmpx(i,j,1) = tmpx(i,j,1) + imask3d(i,j,km)*imask3d(i+1,j,km) * (ssh(i+1,j)*dvol(i+1,j,k)*idx0(i+1,j) - ssh(i,j)*dvol(i,j,k)*idx0(i,j)) * gravity
             tmpx(i,j,2) = tmpx(i,j,2) + imask3d(i,j,km)*imask3d(i+1,j,km) * 0.5 * (dvol(i,j,k)+dvol(i+1,j,k))
          END DO
          END DO

          DO j=0, jsize
          DO i=1, isize
             tmpy(i,j,1) = tmpy(i,j,1) + imask3d(i,j,km)*imask3d(i,j+1,km) * (ssh(i,j+1)*dvol(i,j+1,k)*idy0(i,j+1) - ssh(i,j)*dvol(i,j,k)*idy0(i,j)) * gravity
             tmpy(i,j,2) = tmpy(i,j,2) + imask3d(i,j,km)*imask3d(i,j+1,km) * 0.5 * (dvol(i,j,k)+dvol(i,j+1,k))
          END DO
          END DO
       END DO

       CALL vsum(tmpx, all=.TRUE.)
       CALL vsum(tmpy, all=.TRUE.)

!$OMP PARALLEL DO PRIVATE(i, j, km)
       DO k=1, ksize
          km = maskindex(k)

          DO j=1, jsize
          DO i=0, isize
             IF (tmpx(i,j,2)/=0.0) gx(i,j,k) = gx(i,j,k) - imask3d(i,j,km)*imask3d(i+1,j,km) * tmpx(i,j,1)/tmpx(i,j,2)
          END DO
          END DO

          DO j=0, jsize
          DO i=1, isize
             IF (tmpy(i,j,2)/=0.0) gy(i,j,k) = gy(i,j,k) - imask3d(i,j,km)*imask3d(i,j+1,km) * tmpy(i,j,1)/tmpy(i,j,2)
          END DO
          END DO
       END DO

    END SELECT

  END SUBROUTINE grad_ssh

!-----------------------------------------------------------------------------------------------------------------------


  SUBROUTINE stream_uv(psi, u, v)
    REAL(8), INTENT(IN)  :: psi(0:isize, 0:jsize, 1:ksize)
    REAL(8), INTENT(OUT) :: u( -slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
    REAL(8), INTENT(OUT) :: v(1-slv:isize+slv,  -slv:jsize+slv, 1-slv:ksize+slv)

    INTEGER :: i, j, k

    DO k=1, ksize
       DO j=1, jsize
       DO i=0, isize
          u(i,j,k) = (psi(i,j,k)-psi(i,j-1,k))*0.5*(dz_ref(i,j,k)+dz_ref(i+1,j,k))/dsx_ref(i,j,k)
       END DO
       END DO

       DO j=0, jsize
       DO i=1, isize
          v(i,j,k) = (psi(i,j,k)-psi(i-1,j,k))*0.5*(dz_ref(i,j,k)+dz_ref(i,j+1,k))/dsy_ref(i,j,k)
       END DO
       END DO
    END DO

  END SUBROUTINE stream_uv

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE update_w(u, v, w, topdown)
    REAL(8), INTENT(IN)    :: u( -slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
    REAL(8), INTENT(IN)    :: v(1-slv:isize+slv,  -slv:jsize+slv, 1-slv:ksize+slv)
    REAL(8), INTENT(INOUT) :: w(1-slv:isize+slv, 1-slv:jsize+slv,  -slv:ksize+slv)
    LOGICAL, INTENT(IN), OPTIONAL :: topdown


    INTEGER :: i, j, k
    INTEGER :: kl, km0, km1
    LOGICAL :: topdown_

    topdown_ = .FALSE.
    IF (present(topdown)) topdown_ = topdown

    IF (topdown_) THEN
       w(:,:,ksize) = 0.0
       CALL urecv(w(:,:,ksize))

!$OMP PARALLEL PRIVATE(kl, km0, km1)
       DO k=ksize-1, 0, -1
          kl=dzindex(k)
          km0=maskindex(k)
          km1=maskindex(k+1)

!$OMP DO
          DO j=1, jsize
          DO i=1, isize
             IF (k==surface_k(i,j)) THEN
                w(i,j,k) = (ssh(i,j) - ssh_old(i,j))/dtime
             ELSE
                w(i,j,k) = w(i,j,k+1) - ( u(i-1,j,k)*dsx(i-1,j,kl) - u(i,j,k)*dsx(i,j,kl) &
                                        + v(i,j-1,k)*dsy(i,j-1,kl) - v(i,j,k)*dsy(i,j,kl)) / dsz(i,j)
             END IF
          END DO
          END DO
       END DO
!$OMP END PARALLEL

       CALL lsend(w(:,:,0))
    ELSE
       w(:,:,0) = 0.0
       CALL lrecv(w(:,:,0))

!$OMP PARALLEL PRIVATE(kl, km0, km1)
       DO k=1, ksize
          kl=dzindex(k)
          km0=maskindex(k)
          km1=maskindex(k+1)

!$OMP DO
          DO j=1, jsize
          DO i=1, isize
             IF (imask3d(i,j,km0)*imask3d(i,j,km1) == 0.0) CYCLE
             w(i,j,k) = w(i,j,k-1) + ( u(i-1,j,k)*dsx(i-1,j,kl) - u(i,j,k)*dsx(i,j,kl) &
                                     + v(i,j-1,k)*dsy(i,j-1,kl) - v(i,j,k)*dsy(i,j,kl)) / dsz(i,j)
          END DO
          END DO
       END DO
!$OMP END PARALLEL

       CALL usend(w(:,:,ksize))
    END IF

    CALL update_boundary(w)

  END SUBROUTINE update_w

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE update_velocity_boundary(u, v, w)
    USE parameters
    REAL(8), INTENT(INOUT) :: u( -slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
    REAL(8), INTENT(INOUT) :: v(1-slv:isize+slv,  -slv:jsize+slv, 1-slv:ksize+slv)
    REAL(8), INTENT(INOUT) :: w(1-slv:isize+slv, 1-slv:jsize+slv,  -slv:ksize+slv)
    OPTIONAL w

    CALL update_boundary(u)
    CALL update_boundary(v)

    IF (present(w))  CALL update_boundary(w)

!$OMP PARALLEL
    IF (cycle_x .AND. (u_cycleoffset_x/=0.0 .OR. v_cycleoffset_x/=0.0)) THEN
       IF (icoord==0) THEN
!$OMP WORKSHARE
          u( -slv:0,:,:) = u( -slv:0,:,:) - u_cycleoffset_x
          v(1-slv:0,:,:) = v(1-slv:0,:,:) - v_cycleoffset_x
!$OMP END WORKSHARE
       END IF
       IF (icoord==ipes-1) THEN
!$OMP WORKSHARE
          u(isize+1:isize+slv,:,:) = u(isize+1:isize+slv,:,:) + u_cycleoffset_x
          v(isize+1:isize+slv,:,:) = v(isize+1:isize+slv,:,:) + v_cycleoffset_x
!$OMP END WORKSHARE
       END IF
    END IF

    IF (cycle_y .AND. (u_cycleoffset_y/=0.0 .OR. v_cycleoffset_y/=0.0)) THEN
       IF (jcoord==0) THEN
!$OMP WORKSHARE
          u(:,1-slv:0,:) = u(:,1-slv:0,:) - u_cycleoffset_y
          v(:, -slv:0,:) = v(:, -slv:0,:) - v_cycleoffset_y
!$OMP END WORKSHARE
       END IF
       IF (jcoord==jpes-1) THEN
!$OMP WORKSHARE
          u(:,jsize+1:jsize+slv,:) = u(:,jsize+1:jsize+slv,:) + u_cycleoffset_y
          v(:,jsize+1:jsize+slv,:) = v(:,jsize+1:jsize+slv,:) + v_cycleoffset_y
!$OMP END WORKSHARE
       END IF
    END IF

    IF (cycle_z .AND. (u_cycleoffset_z/=0.0 .OR. v_cycleoffset_z/=0.0)) THEN
       IF (kcoord==0) THEN
!$OMP WORKSHARE
          u(:,:,1-slv:0) = u(:,:,1-slv:0) - u_cycleoffset_z
          v(:,:,1-slv:0) = v(:,:,1-slv:0) - v_cycleoffset_z
!$OMP END WORKSHARE
       END IF
       IF (kcoord==kpes-1) THEN
!$OMP WORKSHARE
          u(:,:,ksize+1:ksize+slv) = u(:,:,ksize+1:ksize+slv) + u_cycleoffset_z
          v(:,:,ksize+1:ksize+slv) = v(:,:,ksize+1:ksize+slv) + v_cycleoffset_z
!$OMP END WORKSHARE
       END IF
    END IF

    IF (tripolar .AND. jcoord==jpes-1) THEN
!$OMP WORKSHARE
       u(:,jsize+1:jsize+slv,:) = -u(:,jsize+1:jsize+slv,:)
       v(:,jsize+1:jsize+slv,:) = -v(:,jsize+1:jsize+slv,:)
!$OMP END WORKSHARE
    END IF
!$OMP END PARALLEL

  END SUBROUTINE update_velocity_boundary

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE restore_velocity(gx, gy, gz)
    REAL(8), INTENT(INOUT) :: gx(0:isize, 1:jsize, 1:ksize)
    REAL(8), INTENT(INOUT) :: gy(1:isize, 0:jsize, 1:ksize)
    REAL(8), INTENT(INOUT) :: gz(1:isize, 1:jsize, 0:ksize)

    REAL(8) :: resrate(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)

    REAL(8) :: restore_u( -slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
    REAL(8) :: restore_v(1-slv:isize+slv,  -slv:jsize+slv, 1-slv:ksize+slv)
    REAL(8) :: restore_w(1-slv:isize+slv, 1-slv:jsize+slv,  -slv:ksize+slv)

    LOGICAL :: stat, stat_u, stat_v, stat_w

    INTEGER :: i, j, k

    CALL checkin('RESTORE_U', restore_u, stat_u)
    CALL checkin('RESTORE_V', restore_v, stat_v)
    CALL checkin('RESTORE_W', restore_w, stat_w)

    IF (.NOT. (stat_u .OR. stat_v .OR. stat_w)) RETURN

    CALL checkin('RESRATE_VELOCITY', resrate, stat)
    IF (.NOT. stat) CALL checkin('RESRATE', resrate, stat)
    CALL assert(stat, "restoring rate (RESRATE) for velocity is not specified.")

!$OMP PARALLEL DO
    DO k=0, ksize
       IF (stat_w) THEN
          DO j=1, jsize
          DO i=1, isize
             IF (restore_w(i,j,k)==UNDEF) CYCLE
             gz(i,j,k) = gz(i,j,k) + imask3d(i,j,k)*imask3d(i,j,k+1) &
                  * (restore_w(i,j,k) - w(i,j,k)) * (resrate(i,j,k)+resrate(i,j,k+1))*0.5
          END DO
          END DO
       END IF

       IF(k==0) CYCLE

       IF (stat_u) THEN
          DO j=1, jsize
          DO i=0, isize
             IF (restore_u(i,j,k)==UNDEF) CYCLE
             gx(i,j,k) = gx(i,j,k) + imask3d(i,j,k)*imask3d(i+1,j,k) &
                  * (restore_u(i,j,k) - u(i,j,k)) * (resrate(i,j,k)+resrate(i+1,j,k))*0.5
          END DO
          END DO
       END IF

       IF (stat_v) THEN
          DO j=0, jsize
          DO i=1, isize
             IF (restore_v(i,j,k)==UNDEF) CYCLE
             gy(i,j,k) = gy(i,j,k) + imask3d(i,j,k)*imask3d(i,j+1,k) &
                  * (restore_v(i,j,k) - v(i,j,k)) * (resrate(i,j,k)+resrate(i,j+1,k))*0.5
          END DO
          END DO
       END IF
    END DO

  END SUBROUTINE restore_velocity

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE kinetic_energy(u, v, w, ke)
    REAL(8), INTENT(IN)  :: u( -slv:isize+slv,1-slv:jsize+slv,1-slv:ksize+slv)
    REAL(8), INTENT(IN)  :: v(1-slv:isize+slv, -slv:jsize+slv,1-slv:ksize+slv)
    REAL(8), INTENT(IN)  :: w(1-slv:isize+slv,1-slv:jsize+slv, -slv:ksize+slv)
    REAL(8), INTENT(OUT) :: ke(0:isize+1, 0:jsize+1, 0:ksize+1)

    INTEGER :: i, j, k

!$OMP PARALLEL DO
    DO k=0, ksize+1
    DO j=0, jsize+1
    DO i=0, isize+1
       ke(i,j,k) = 0.25 * ( u(i-1,j,k)**2 + u(i,j,k)**2 &
                          + v(i,j-1,k)**2 + v(i,j,k)**2 &
                          + w(i,j,k-1)**2 + w(i,j,k)**2)
    END DO
    END DO
    END DO

  END SUBROUTINE kinetic_energy

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE strain_rate(dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz, srate)
    REAL(8), INTENT(IN)  :: dudx( 0:isize+1, 0:jsize+1, 0:ksize+1)
    REAL(8), INTENT(IN)  :: dudy(-1:isize+1,-1:jsize+1, 0:ksize+1)
    REAL(8), INTENT(IN)  :: dudz(-1:isize+1, 0:jsize+1,-1:ksize+1)
    REAL(8), INTENT(IN)  :: dvdx(-1:isize+1,-1:jsize+1, 0:ksize+1)
    REAL(8), INTENT(IN)  :: dvdy( 0:isize+1, 0:jsize+1, 0:ksize+1)
    REAL(8), INTENT(IN)  :: dvdz( 0:isize+1,-1:jsize+1,-1:ksize+1)
    REAL(8), INTENT(IN)  :: dwdx(-1:isize+1, 0:jsize+1,-1:ksize+1)
    REAL(8), INTENT(IN)  :: dwdy( 0:isize+1,-1:jsize+1,-1:ksize+1)
    REAL(8), INTENT(IN)  :: dwdz( 0:isize+1, 0:jsize+1, 0:ksize+1)
    REAL(8), INTENT(OUT) :: srate(0:isize+1, 0:jsize+1, 0:ksize+1)

    INTEGER :: i, j, k

!$OMP PARALLEL DO
    DO k=0, ksize+1
    DO j=0, jsize+1
    DO i=0, isize+1
       srate(i,j,k) = dudx(i,j,k)**2 + dvdy(i,j,k)**2 + dwdz(i,j,k)**2 &
            + 0.125 * ( (dudy(i-1,j-1,k)+dvdx(i-1,j-1,k))**2 + (dudy(i,j-1,k)+dvdx(i,j-1,k))**2 &
                      + (dudy(i-1,j,  k)+dvdx(i-1,j,  k))**2 + (dudy(i,j,  k)+dvdx(i,j,  k))**2 &
                      + (dvdz(i,j-1,k-1)+dwdy(i,j-1,k-1))**2 + (dvdz(i,j,k-1)+dwdy(i,j,k-1))**2 &
                      + (dvdz(i,j-1,k  )+dwdy(i,j-1,k  ))**2 + (dvdz(i,j,k  )+dwdy(i,j,k  ))**2 &
                      + (dwdx(i-1,j,k-1)+dudz(i-1,j,k-1))**2 + (dwdx(i-1,j,k)+dudz(i-1,j,k))**2 &
                      + (dwdx(i,  j,k-1)+dudz(i,  j,k-1))**2 + (dwdx(i,  j,k)+dudz(i,  j,k))**2)

       srate(i,j,k) = sqrt(srate(i,j,k))
    END DO
    END DO
    END DO

  END SUBROUTINE strain_rate

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE shear_freq2(dudz, dvdz, sfreq2)
    REAL(8), INTENT(IN)  :: dudz( -1:isize+1, 0:jsize+1,-1:ksize+1)
    REAL(8), INTENT(IN)  :: dvdz(  0:isize+1,-1:jsize+1,-1:ksize+1)
    REAL(8), INTENT(OUT) :: sfreq2(0:isize+1, 0:jsize+1,-1:ksize+1)

    INTEGER :: i, j, k

!$OMP PARALLEL DO
    DO k =-1, ksize+1
    DO j = 0, jsize+1
    DO i = 0, isize+1
       sfreq2(i,j,k) = 0.25*((dudz(i-1,j,k)+dudz(i,j,k))**2 + (dvdz(i,j-1,k)+dvdz(i,j,k))**2)
    END DO
    END DO
    END DO

  END SUBROUTINE shear_freq2

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE enstrophy(result)
    REAL(8), INTENT(OUT) :: result(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)

    INTEGER :: i, j, k

!$OMP PARALLEL DO
    DO k=1, ksize
    DO j=1, jsize
    DO i=1, isize
       result(i,j,k) = 0.125 * ( (dwdy(i,j-1,k-1)-dvdz(i,j-1,k-1))**2 + (dwdy(i,j,k-1)-dvdz(i,j,k-1))**2 &
                                 + (dwdy(i,j-1,k  )-dvdz(i,j-1,k  ))**2 + (dwdy(i,j,k  )-dvdz(i,j,k  ))**2 &
                                 + (dudz(i-1,j,k-1)-dwdx(i-1,j,k-1))**2 + (dudz(i-1,j,k)-dwdx(i-1,j,k))**2 &
                                 + (dudz(i  ,j,k-1)-dwdx(i  ,j,k-1))**2 + (dudz(i,  j,k)-dwdx(i,  j,k))**2 &
                                 + (dvdx(i-1,j-1,k)-dudy(i-1,j-1,k))**2 + (dvdx(i,j-1,k)-dudy(i,j-1,k))**2 &
                                 + (dvdx(i-1,j,  k)-dudy(i-1,j,  k))**2 + (dvdx(i,j,  k)-dudy(i,j,  k))**2)
    END DO
    END DO
    END DO

  END SUBROUTINE enstrophy

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE lambda2(result)
    REAL(8), INTENT(OUT) :: result(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)

    INTEGER :: i, j, k
    REAL(8) :: a(3,3)
    REAL(8) :: l1, l2, l3

!$OMP PARALLEL DO PRIVATE(a, l1, l2, l3)
    DO k=1, ksize
    DO j=1, jsize
    DO i=1, isize
       a(1,1) = dudx(i,j,k)
       a(1,2) = 0.25*(dudy(i-1,j-1,k)+dudy(i-1,j,k)+dudy(i,j-1,k)+dudy(i,j,k))
       a(1,3) = 0.25*(dudz(i-1,j,k-1)+dudz(i-1,j,k)+dudz(i,j,k-1)+dudz(i,j,k))

       a(2,1) = 0.25*(dvdx(i-1,j-1,k)+dvdx(i,j-1,k)+dvdx(i-1,j,k)+dvdx(i,j,k))
       a(2,2) = dvdy(i,j,k)
       a(2,3) = 0.25*(dvdz(i,j-1,k-1)+dvdz(i,j-1,k)+dvdz(i,j,k-1)+dvdz(i,j,k))

       a(3,1) = 0.25*(dwdx(i-1,j,k-1)+dwdx(i,j,k-1)+dwdx(i-1,j,k)+dwdx(i,j,k))
       a(3,2) = 0.25*(dwdy(i,j-1,k-1)+dwdy(i,j,k-1)+dwdy(i,j-1,k)+dwdy(i,j,k))
       a(3,3) = dwdz(i,j,k)

       a = matmul(a, a)
       a(1,2) = 0.5 * (a(1,2)+a(2,1))
       a(2,3) = 0.5 * (a(2,3)+a(3,2))
       a(3,1) = 0.5 * (a(3,1)+a(1,3))
       a(2,1) = a(1,2)
       a(3,2) = a(2,3)
       a(1,3) = a(3,1)

       CALL eigenvalues_sym3x3(a, l1, l2, l3)

       result(i,j,k) = l2
    END DO
    END DO
    END DO

  END SUBROUTINE lambda2

!-----------------------------------------------------------------------------------------------------------------------

!calculate Okubo-Weiss parameter that detect 2D eddy
  SUBROUTINE owparam(result)
    REAL(8), INTENT(OUT) :: result(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)

    INTEGER :: i, j, k

!$OMP PARALLEL DO
    DO k=1, ksize
    DO j=1, jsize
    DO i=1, isize
       result(i,j,k) = (dudx(i,j,k)-dvdy(i,j,k))**2 + 0.25*(dvdx(i,j,k)+dvdx(i-1,j,k)+dvdx(i,j-1,k)+dvdx(i-1,j-1,k)) &
                                                          *(dudy(i,j,k)+dudy(i-1,j,k)+dudy(i,j-1,k)+dudy(i-1,j-1,k))
    END DO
    END DO
    END DO

  END SUBROUTINE owparam

  SUBROUTINE owparam_xz(result)
    REAL(8), INTENT(OUT) :: result(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)

    INTEGER :: i, j, k

!$OMP PARALLEL DO
    DO k=1, ksize
    DO j=1, jsize
    DO i=1, isize
       result(i,j,k) = (dudx(i,j,k)-dwdz(i,j,k))**2 + 0.25*(dwdx(i,j,k)+dwdx(i-1,j,k)+dwdx(i,j,k-1)+dwdx(i-1,j,k-1)) &
                                                          *(dudz(i,j,k)+dudz(i-1,j,k)+dudz(i,j,k-1)+dudz(i-1,j,k-1))
    END DO
    END DO
    END DO

  END SUBROUTINE owparam_xz

  SUBROUTINE owparam_yz(result)
    REAL(8), INTENT(OUT) :: result(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)

    INTEGER :: i, j, k

!$OMP PARALLEL DO
    DO k=1, ksize
    DO j=1, jsize
    DO i=1, isize
       result(i,j,k) = (dvdy(i,j,k)-dwdz(i,j,k))**2 + 0.25*(dwdy(i,j,k)+dwdy(i,j-1,k)+dwdy(i,j,k-1)+dwdy(i,j-1,k-1)) &
                                                          *(dvdz(i,j,k)+dvdz(i,j-1,k)+dvdz(i,j,k-1)+dvdz(i,j-1,k-1))
    END DO
    END DO
    END DO

  END SUBROUTINE owparam_yz

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE principal_strain_rate(s, ax, ay, az)
    REAL(8), INTENT(OUT) :: s(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
    REAL(8), INTENT(OUT) :: ax(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
    REAL(8), INTENT(OUT) :: ay(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)
    REAL(8), INTENT(OUT) :: az(1-slv:isize+slv, 1-slv:jsize+slv, 1-slv:ksize+slv)

    INTEGER :: i, j, k
    REAL(8) :: a(3,3), v(3), a2, v2
    REAL(8) :: l1, l2, l3

    REAL(8), PARAMETER :: eps = 1.0D-30

!$OMP PARALLEL DO PRIVATE(a, l1, l2, l3)
    DO k=1, ksize
    DO j=1, jsize
    DO i=1, isize
       a(1,1) = dudx(i,j,k)
       a(2,2) = dvdy(i,j,k)
       a(3,3) = dwdz(i,j,k)
       a(1,2) = 0.125*( dudy(i-1,j-1,k)+dudy(i-1,j,k)+dudy(i,j-1,k)+dudy(i,j,k) &
                        + dvdx(i-1,j-1,k)+dvdx(i-1,j,k)+dvdx(i,j-1,k)+dvdx(i,j,k))
       a(2,3) = 0.125*( dvdz(i,j-1,k-1)+dvdz(i,j-1,k)+dvdz(i,j,k-1)+dvdz(i,j,k) &
                        + dwdy(i,j-1,k-1)+dwdy(i,j-1,k)+dwdy(i,j,k-1)+dwdy(i,j,k))
       a(3,1) = 0.125*( dwdx(i-1,j,k-1)+dwdx(i,j,k-1)+dwdx(i-1,j,k)+dwdx(i,j,k) &
                        + dudz(i-1,j,k-1)+dudz(i,j,k-1)+dudz(i-1,j,k)+dudz(i,j,k))
       a(2,1) = a(1,2)
       a(3,2) = a(2,3)
       a(1,3) = a(3,1)

       a2 = sum(a(:,:)**2)

       IF (maxval(abs(a(1,:))) == 0.0) THEN
          CALL eigenvalues_sym2x2(a(2:3,2:3), l1, l2)
          IF (abs(l1) > abs(l2)) THEN
             s(i,j,k) = l1
          ELSE
             s(i,j,k) = l2
          END IF
          CALL eigenvector_2x2(a(2:3,2:3), s(i,j,k), v(1:2))
          ax(i,j,k) = 0.0
          ay(i,j,k) = v(1)
          az(i,j,k) = v(2)
       ELSE IF (maxval(abs(a(2,:))) == 0.0) THEN
          a(1,2) = a(1,3)
          a(2,1) = a(3,1)
          a(2,2) = a(3,3)
          CALL eigenvalues_sym2x2(a(1:2,1:2), l1, l2)
          IF (abs(l1) > abs(l2)) THEN
             s(i,j,k) = l1
          ELSE
             s(i,j,k) = l2
          END IF
          CALL eigenvector_2x2(a(1:2,1:2), s(i,j,k), v(1:2))
          ax(i,j,k) = v(1)
          ay(i,j,k) = 0.0
          az(i,j,k) = v(2)
       ELSE IF (maxval(abs(a(3,:))) == 0.0) THEN
          CALL eigenvalues_sym2x2(a(1:2,1:2), l1, l2)
          IF (abs(l1) > abs(l2)) THEN
             s(i,j,k) = l1
          ELSE
             s(i,j,k) = l2
          END IF
          CALL eigenvector_2x2(a(1:2,1:2), s(i,j,k), v(1:2))
          ax(i,j,k) = v(1)
          ay(i,j,k) = v(2)
          az(i,j,k) = 0.0
       ELSE
          CALL eigenvalues_sym3x3(a, l1, l2, l3)

          IF (abs(l1-l2) < l3*eps) THEN
             s(i,j,k) = l3
          ELSE IF (abs(l3-l2) < l1*eps) THEN
             s(i,j,k) = l1
          ELSE IF (abs(l1) > abs(l3)) THEN
             s(i,j,k) = l1
          ELSE
             s(i,j,k) = l3
          END IF

          CALL eigenvector_3x3(a, s(i,j,k), v)
          ax(i,j,k) = v(1)
          ay(i,j,k) = v(2)
          az(i,j,k) = v(3)
       END IF

       ! v2 = ax(i,j,k)**2 + ay(i,j,k)**2 + az(i,j,k)**2
       ! IF (v2 < eps**2) THEN
       !    ax(i,j,k) = 0.0
       !    ay(i,j,k) = 0.0
       !    az(i,j,k) = 0.0
       ! ELSE
       !    ax(i,j,k) = ax(i,j,k)/sqrt(v2)
       !    ay(i,j,k) = ay(i,j,k)/sqrt(v2)
       !    az(i,j,k) = az(i,j,k)/sqrt(v2)
       ! END IF
    END DO
    END DO
    END DO
  END SUBROUTINE principal_strain_rate

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE eigenvalues_sym3x3(a, l1, l2, l3)
    REAL(8), INTENT(IN)  :: a(3,3)
    REAL(8), INTENT(OUT) :: l1
    REAL(8), INTENT(OUT) :: l2
    REAL(8), INTENT(OUT) :: l3

    REAL(8) :: p, q, r
    REAL(8) :: alpha, beta
    REAL(8) :: t1, t2, t3

    COMPLEX(8) :: c, o

    o = CMPLX(-0.5D0, 0.5D0*sqrt(3.0D0), KIND=8)

    p = a(1,1) + a(2,2) + a(3,3)

    q =   a(1,1)*a(2,2) + a(2,2)*a(3,3) + a(3,3)*a(1,1) &
        - a(1,2)*a(2,1) - a(2,3)*a(3,2) - a(3,1)*a(1,3)

    r =   a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2)) &
        + a(1,2)*(a(2,3)*a(3,1)-a(2,1)*a(3,3)) &
        + a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))

    alpha = (p**2)/9.0 - q/3.0
    beta  = r/2.0 - p*q/6.0 + (p**3)/27.0

    c = CMPLX(beta, sqrt(max(0.0D0, alpha**3 - beta**2)), KIND=8)**(1/3.0)

    t1 = 2*REAL(c,          KIND=8) + p/3.0
    t2 = 2*REAL(c*o,        KIND=8) + p/3.0
    t3 = 2*REAL(c*conjg(o), KIND=8) + p/3.0

    l1 = max(t1, t2, t3)
    l3 = min(t1, t2, t3)
    l2 = t1 + t2 + t3 - l1 - l3

  END SUBROUTINE eigenvalues_sym3x3

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE eigenvalues_sym2x2(a, l1, l2)
    REAL(8), INTENT(IN)  :: a(2,2)
    REAL(8), INTENT(OUT) :: l1
    REAL(8), INTENT(OUT) :: l2

    l1 = (a(1,1)+a(2,2) + sqrt((a(1,1)-a(2,2))**2 + 4*a(1,2)**2)) / 2.0
    l2 = (a(1,1)+a(2,2) - sqrt((a(1,1)-a(2,2))**2 + 4*a(1,2)**2)) / 2.0
  END SUBROUTINE eigenvalues_sym2x2

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE eigenvector_3x3(a, lambda, v)
    REAL(8), INTENT(IN)  :: a(3,3)
    REAL(8), INTENT(IN)  :: lambda
    REAL(8), INTENT(OUT) :: v(3)

    v(1) = (a(1,2)*a(2,3) - a(1,3)*(a(2,2)-lambda)) * (a(2,3)*a(3,1) - a(2,1)*(a(3,3)-lambda))
    v(2) = (a(2,1)*a(1,3) - a(2,3)*(a(1,1)-lambda)) * (a(2,3)*a(3,1) - a(2,1)*(a(3,3)-lambda))
    v(3) = (a(2,1)*a(1,3) - a(2,3)*(a(1,1)-lambda)) * (a(3,2)*a(2,1) - a(3,1)*(a(2,2)-lambda))

  END SUBROUTINE eigenvector_3x3

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE eigenvector_2x2(a, lambda, v)
    REAL(8), INTENT(IN)  :: a(2,2)
    REAL(8), INTENT(IN)  :: lambda
    REAL(8), INTENT(OUT) :: v(2)

    v(1) = lambda - a(2,2)
    v(2) = a(2,1)
  END SUBROUTINE eigenvector_2x2

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE checkout_velocity
    REAL(8) :: tmp_xy(isize,jsize)
    REAL(8) :: tmp_xz(isize,ksize)
    REAL(8) :: tmp_yz(jsize,ksize)

    INTEGER :: i, j, k, n

    CALL checkout('U', u) ! u, x-component of velocity [m/s]
    CALL checkout('V', v) ! v, y-component of velocity [m/s]
    CALL checkout('W', w) ! w, z-component of velocity [m/s]

    IF (open_e) THEN
       IF (icoord/=ipes-1) tmp_yz(:,:) = 0.0

       IF(require_checkout('U_EAST')) THEN
          IF (icoord==ipes-1) tmp_yz(:,:) = u(isize,1:jsize,1:ksize)
          CALL checkout('U_EAST', tmp_yz, section='YZ')
       END IF

       IF(require_checkout('V_EAST')) THEN
          IF (icoord==ipes-1)  tmp_yz(:,:) = v(isize+1,1:jsize,1:ksize)
          CALL checkout('V_EAST', tmp_yz, section='YZ')
       END IF

       IF(require_checkout('W_EAST')) THEN
          IF (icoord==ipes-1) tmp_yz(:,:) = w(isize+1,1:jsize,1:ksize)
          CALL checkout('W_EAST', tmp_yz, section='YZ')
       END IF
    END IF

    IF (open_w) THEN
       IF (icoord/=0) tmp_yz(:,:) = 0.0

       IF(require_checkout('U_WEST')) THEN
          IF (icoord==0) tmp_yz(:,:) = u(0,1:jsize,1:ksize)
          CALL checkout('U_WEST', tmp_yz, section='YZ')
       END IF

       IF(require_checkout('V_WEST')) THEN
          IF (icoord==0)  tmp_yz(:,:) = v(0,1:jsize,1:ksize)
          CALL checkout('V_WEST', tmp_yz, section='YZ')
       END IF

       IF(require_checkout('W_WEST')) THEN
          IF (icoord==0) tmp_yz(:,:) = w(0,1:jsize,1:ksize)
          CALL checkout('W_WEST', tmp_yz, section='YZ')
       END IF
    END IF

    IF (open_n) THEN
       IF (jcoord/=jpes-1) tmp_xz(:,:) = 0.0

       IF(require_checkout('U_NORTH')) THEN
          IF (jcoord==jpes-1) tmp_xz(:,:) = u(1:isize,jsize+1,1:ksize)
          CALL checkout('U_NORTH', tmp_xz, section='XZ')
       END IF

       IF(require_checkout('V_NORTH')) THEN
          IF (jcoord==jpes-1)  tmp_xz(:,:) = v(1:isize,jsize,1:ksize)
          CALL checkout('V_NORTH', tmp_xz, section='XZ')
       END IF

       IF(require_checkout('W_NORTH')) THEN
          IF (jcoord==jpes-1) tmp_xz(:,:) = w(1:isize,jsize+1,1:ksize)
          CALL checkout('W_NORTH', tmp_xz, section='XZ')
       END IF
    END IF

    IF (open_s) THEN
       IF (jcoord/=0) tmp_xz(:,:) = 0.0

       IF(require_checkout('U_SOUTH')) THEN
          IF (jcoord==0) tmp_xz(:,:) = u(1:isize,0,1:ksize)
          CALL checkout('U_SOUTH', tmp_xz, section='XZ')
       END IF

       IF(require_checkout('V_SOUTH')) THEN
          IF (jcoord==0)  tmp_xz(:,:) = v(1:isize,0,1:ksize)
          CALL checkout('V_SOUTH', tmp_xz, section='XZ')
       END IF

       IF(require_checkout('W_SOUTH')) THEN
          IF (jcoord==0) tmp_xz(:,:) = w(1:isize,0,1:ksize)
          CALL checkout('W_SOUTH', tmp_xz, section='XZ')
       END IF
    END IF

    CALL checkout('DUDX', dudx) ! \partial u / \partial x [1/s]
    CALL checkout('DUDY', dudy) ! \partial u / \partial y [1/s]
    CALL checkout('DUDZ', dudz) ! \partial u / \partial z [1/s]
    CALL checkout('DVDX', dvdx) ! \partial v / \partial x [1/s]
    CALL checkout('DVDY', dvdy) ! \partial v / \partial y [1/s]
    CALL checkout('DVDZ', dvdz) ! \partial v / \partial z [1/s]
    CALL checkout('DWDX', dwdx) ! \partial w / \partial x [1/s]
    CALL checkout('DWDY', dwdy) ! \partial w / \partial y [1/s]
    CALL checkout('DWDZ', dwdz) ! \partial w / \partial z [1/s]

    IF (require_checkout('OMEGAX')) CALL checkout('OMEGAX', dwdy - dvdz) ! omega_x, x-component of vorticity [1/s]
    IF (require_checkout('OMEGAY')) CALL checkout('OMEGAY', dudz - dwdx) ! omega_y, y-component of vorticity [1/s]
    IF (require_checkout('OMEGAZ')) CALL checkout('OMEGAZ', dvdx - dudy) ! omega_z, z-component of vorticity [1/s]

    CALL checkout('SFREQ2', sfreq2) ! square of vertical shear frequency [1/s^2]

    CALL vsum(taux_sfc)
    CALL vsum(tauy_sfc)

    CALL checkout('SFC_TAUX', taux_sfc)  ! x-component of surface stress [kg / (m s^2) = Pa]
    CALL checkout('SFC_TAUY', tauy_sfc)  ! y-component of surface stress [kg / (m s^2) = Pa]

    IF (require_checkout('SFC_TAU')) THEN
       DO j=1, jsize
       DO i=1, isize
          tmp_xy(i,j) = sqrt(0.25*(taux_sfc(i-1,j)+taux_sfc(i,j))**2 + 0.25*(tauy_sfc(i,j-1)+tauy_sfc(i,j))**2)
       END DO
       END DO
       CALL checkout('SFC_TAU', tmp_xy) ! absolute value of surface stress [kg / (m s^2) = Pa]
    END IF

    CALL vsum(taux_btm)
    CALL vsum(tauy_btm)

    CALL checkout('BTM_TAUX', taux_btm) !  x-component of bottom stress [kg / (m s^2) = Pa]
    CALL checkout('BTM_TAUY', tauy_btm) !  x-component of bottom stress [kg / (m s^2) = Pa]

    IF (require_checkout('BTM_TAU')) THEN
       DO j=1, jsize
       DO i=1, isize
          tmp_xy(i,j) = sqrt(0.25*(taux_btm(i-1,j)+taux_btm(i,j))**2 + 0.25*(tauy_btm(i,j-1)+tauy_btm(i,j))**2)
       END DO
       END DO
       CALL checkout('BTM_TAU', tmp_xy) !  absolute value of bottom stress [kg / (m s^2) = Pa]
    END IF

!$OMP PARALLEL WORKSHARE
    taux_sfc(:,:) = 0.0
    tauy_sfc(:,:) = 0.0
    taux_btm(:,:) = 0.0
    tauy_btm(:,:) = 0.0
!$OMP END PARALLEL WORKSHARE

    IF (require_checkout('STREAM'))    CALL checkout_stream
    IF (require_checkout('STREAM_XZ')) CALL checkout_stream_xz
    IF (require_checkout('STREAM_YZ')) CALL checkout_stream_yz

  CONTAINS
    SUBROUTINE checkout_stream
      REAL(8) :: uvint(2, 0:isize, 0:jsize)
      REAL(8) :: psi(0:isize, 0:jsize)

      INTEGER :: i, j, k
      INTEGER :: ierr

!$OMP PARALLEL PRIVATE(i,j)
!$OMP WORKSHARE
      uvint(:,:,:) = 0.0
!$OMP END WORKSHARE

      DO k=1, ksize
!$OMP DO COLLAPSE(2)
         DO j=1, jsize
         DO i=0, isize
            uvint(1,i,j) = uvint(1,i,j) + u(i,j,k)*dsx(i,j,k)
         END DO
         END DO

!$OMP DO COLLAPSE(2)
         DO j=0, jsize
         DO i=1, isize
            uvint(2,i,j) = uvint(2,i,j) + v(i,j,k)*dsy(i,j,k)
         END DO
         END DO
      END DO
!$OMP END PARALLEL
      CALL vsum(uvint)

      IF (vrank==0) THEN
         psi(0,0) = 0.0

#ifdef PARALLEL_MPI
         IF (icoord > 0) CALL mpi_recv(psi(0,0), 1, MPI_REAL8, rank_w, 0, comm, MPI_STATUS_IGNORE, ierr)
         IF (jcoord > 0) CALL mpi_recv(psi(0,0), 1, MPI_REAL8, rank_s, 1, comm, MPI_STATUS_IGNORE, ierr)
#endif

         DO i=1, isize
            psi(i,0) = psi(i-1,0)-uvint(2,i,0)
         END DO

         DO j=1, jsize
            DO i=0, isize
               psi(i,j) = psi(i,j-1)+uvint(1,i,j)
            END DO
         END DO

#ifdef PARALLEL_MPI
         IF (icoord < ipes-1) CALL mpi_send(psi(isize,0), 1, MPI_REAL8, rank_e, 0, comm, ierr)
         IF (jcoord < jpes-1) CALL mpi_send(psi(0,jsize), 1, MPI_REAL8, rank_n, 1, comm, ierr)
#endif
      ELSE
         psi(0,0) = 0.0
      END IF

      CALL checkout('STREAM', psi)

    END SUBROUTINE checkout_stream

    SUBROUTINE checkout_stream_xz
      REAL(8) :: uint(0:isize, 1:ksize)
      REAL(8) :: psi(0:isize, 0:ksize)

      INTEGER :: i, j, k
      INTEGER :: ierr

!$OMP PARALLEL PRIVATE(i,j)
!$OMP WORKSHARE
      uint(:,:) = 0.0
!$OMP END WORKSHARE

!$OMP DO
      DO k=1, ksize
         DO j=1, jsize
         DO i=0, isize
            uint(i,k) = uint(i,k) + u(i,j,k)*dsx(i,j,k)
         END DO
         END DO
      END DO
!$OMP END PARALLEL

      psi(:,0) = 0.0

      IF (cycle_z .AND. kcoord==0) THEN
#ifdef PARALLEL_MPI
         IF (icoord > 0) CALL mpi_recv(psi(0,0), 1, MPI_REAL8, rank_w, 0, comm, MPI_STATUS_IGNORE, ierr)
#endif
         DO i=1, isize
            psi(i,0) = psi(i-1,0) + sum(w(i,1:jsize,0)*dsz(i,1:jsize))
         END DO
#ifdef PARALLEL_MPI
         IF (icoord < ipes-1) CALL mpi_send(psi(isize,0), 1, MPI_REAL8, rank_e, 0, comm, ierr)
#endif
      ELSE
         CALL lrecv(psi(:,0))
      END IF

      DO k=1, ksize
         DO i=0, isize
            psi(i,k) = psi(i,k-1)+uint(i,k)
         END DO
      END DO

      IF (kcoord < kpes-1) CALL usend(psi(:,ksize))

      CALL checkout('STREAM_XZ', psi, section='XZ')

    END SUBROUTINE checkout_stream_xz

!-----------------------------------------------------------------------------------------------------------------------

    SUBROUTINE checkout_stream_yz
      REAL(8) :: vint(0:jsize, 1:ksize)
      REAL(8) :: psi(0:jsize, 0:ksize)

      INTEGER :: i, j, k
      INTEGER :: ierr

!$OMP PARALLEL PRIVATE(i,j)
!$OMP WORKSHARE
      vint(:,:) = 0.0
!$OMP END WORKSHARE

!$OMP DO
      DO k=1, ksize
         DO j=0, jsize
         DO i=1, isize
            vint(j,k) = vint(j,k) + v(i,j,k)*dsy(i,j,k)
         END DO
         END DO
      END DO
!$OMP END PARALLEL

      psi(:,0) = 0.0

      IF (cycle_z .AND. kcoord==0) THEN
#ifdef PARALLEL_MPI
         IF (jcoord > 0) CALL mpi_recv(psi(0,0), 1, MPI_REAL8, rank_s, 0, comm, MPI_STATUS_IGNORE, ierr)
#endif
         DO j=1, jsize
            psi(j,0) = psi(j-1,0) + sum(w(1:isize,j,0)*dsz(1:isize,j))
         END DO

#ifdef PARALLEL_MPI
         IF (jcoord < jpes-1) CALL mpi_send(psi(jsize,0), 1, MPI_REAL8, rank_n, 0, comm, ierr)
#endif
      ELSE
         CALL lrecv(psi(:,0))
      END IF

      DO k=1, ksize
         DO j=0, jsize
            psi(j,k) = psi(j,k-1)+vint(j,k)
         END DO
      END DO

      IF (kcoord < kpes-1) CALL usend(psi(:,ksize))

      CALL checkout('STREAM_YZ', psi, section='YZ')

    END SUBROUTINE checkout_stream_yz

  END SUBROUTINE checkout_velocity

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE report_cfl(force)
    USE parameters
    USE calendar

    LOGICAL, INTENT(IN), OPTIONAL :: force

    REAL(8) :: cfl(1:isize,1:jsize,1:ksize)
    REAL(8) :: maxcfl(3)
    INTEGER :: loc(3)

    INTEGER :: i, j, k

    IF (.NOT. present(force)) THEN
       IF (cfl_report_interval <= 0)              RETURN
       IF (mod(n_timestep, cfl_report_interval)/=0) RETURN
    ELSE
       IF (.NOT. force) RETURN
    END IF

!$OMP PARALLEL DO
    DO k=1, ksize
    DO j=1, jsize
    DO i=1, isize
       cfl(i,j,k) = abs(u(i,j,k)*dtime*idx1(i,j))
    END DO
    END DO
    END DO

    maxcfl(1) = maxval(cfl)

!$OMP PARALLEL DO
    DO k=1, ksize
    DO j=1, jsize
    DO i=1, isize
       cfl(i,j,k) = abs(v(i,j,k)*dtime*idy1(i,j))
    END DO
    END DO
    END DO

    maxcfl(2) = maxval(cfl)

!$OMP PARALLEL DO
    DO k=1, ksize
    DO j=1, jsize
    DO i=1, isize
       cfl(i,j,k) = abs(w(i,j,k)*dtime*idz1(k))
    END DO
    END DO
    END DO

    maxcfl(3) = maxval(cfl)
    CALL gmax(maxcfl)

    IF (rank==0) WRITE(REPORT_UNIT, '(A, ": max CFL = ", ES10.3, ", ", ES10.3, ", ", ES10.3)') trim(current_datetime), maxcfl(:)

  END SUBROUTINE report_cfl

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE report_open
    REAL(8) :: balance(4)
    INTEGER :: i, j, k

    IF (open_report_interval <= 0) RETURN

    IF (mod(n_timestep, open_report_interval) == 0) THEN
       balance(:) = 0.0

!$OMP PARALLEL DO PRIVATE(i,j,k) REDUCTION(+:balance)
       DO k=1, ksize
          DO j=1, jsize
             balance(1) = balance(1) + u(    0,j,k)*dsx(    0,j,k)
             balance(2) = balance(2) - u(isize,j,k)*dsx(isize,j,k)
          END DO
          DO i=1, isize
             balance(3) = balance(3) + v(i,    0,k)*dsy(i,    0,k)
             balance(4) = balance(4) - v(i,jsize,k)*dsy(i,jsize,k)
          END DO
       END DO

       IF (.NOT. icoord==0)      balance(1) = 0.0
       IF (.NOT. icoord==ipes-1) balance(2) = 0.0
       IF (.NOT. jcoord==0)      balance(3) = 0.0
       IF (.NOT. jcoord==jpes-1) balance(4) = 0.0

       CALL gsum(balance)

       IF (rank==0) THEN
          WRITE(REPORT_UNIT, '(A, ": net volume balance at W/E/S/N = ", ES10.3, ", ", ES10.3, ", ", ES10.3, ", ", ES10.3, " [m^3/s]")') &
               trim(current_datetime), balance
       END IF

    ENDIF

  END SUBROUTINE report_open


  SUBROUTINE report_momentum(force)
    USE parameters
    USE calendar

    LOGICAL, INTENT(IN), OPTIONAL :: force

    REAL(8) :: mom(3)

    INTEGER :: i, j, k

    IF (.NOT. present(force)) THEN
       IF (momentum_report_interval <= 0)              RETURN
       IF (mod(n_timestep, momentum_report_interval)/=0) RETURN
    ELSE
       IF (.NOT. force) RETURN
    END IF

    mom(:) = 0.0

!$OMP PARALLEL DO REDUCTION(+:mom)
    DO k=1, ksize
    DO j=1, jsize
    DO i=1, isize
       mom(1) = mom(1) + u(i,j,k)*(dvol(i,j,k)+dvol(i+1,j,k))
       mom(2) = mom(2) + v(i,j,k)*(dvol(i,j,k)+dvol(i,j+1,k))
       mom(3) = mom(3) + w(i,j,k)*(dvol(i,j,k)+dvol(i,j,k+1))
    END DO
    END DO
    END DO

    CALL gsum(mom)

    mom(:) = mom(:) * rho_0 * 0.5

    IF (rank==0) THEN
       WRITE(REPORT_UNIT, '(A, ": total momentum = ", es10.3, ", ", es10.3, ", ", es10.3, " [kg m/s]")') &
            trim(current_datetime), mom(1), mom(2), mom(3)
    ENDIF

  END SUBROUTINE report_momentum

!-----------------------------------------------------------------------------------------------------------------------

  SUBROUTINE report_energy(force)
    USE parameters
    USE calendar

    LOGICAL, INTENT(IN), OPTIONAL :: force

    REAL(8) :: tke, tpe

    INTEGER :: i, j, k

    IF (.NOT. present(force)) THEN
       IF (energy_report_interval <= 0)                RETURN
       IF (mod(n_timestep, energy_report_interval)/=0) RETURN
    ELSE
       IF (.NOT. force) RETURN
    END IF

    tke = 0.0

!$OMP PARALLEL DO REDUCTION(+:tke)
    DO k=1, ksize
    DO j=1, jsize
    DO i=1, isize
       tke = tke + ke(i,j,k)*dvol(i,j,k)*imask3d(i,j,k)
    END DO
    END DO
    END DO

    CALL gsum(tke)
    tke = tke * rho_0

    tpe = 0.0
    IF (vrank==0) THEN
       DO j=1, jsize
       DO i=1, isize
          tpe = tpe + 0.5*ssh(i,j)**2 * dsz(i,j) * imask2d(i,j)
       END DO
       END DO
    END IF

    CALL gsum(tpe)
    tpe = tpe * rho_0 * gravity

    IF (rank==0) THEN
       WRITE(REPORT_UNIT, '(A, ": total kinetic energy and free-surface potential energy =", es10.3, ", ", es10.3, " [J]")') &
            trim(current_datetime), tke, tpe
    ENDIF

  END SUBROUTINE report_energy

!-----------------------------------------------------------------------------------------------------------------------

END MODULE velocity
