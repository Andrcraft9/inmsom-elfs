!=====================================================
subroutine cyclize_x(ff,nx,ny,nz,mmm,mm)
    use mpi_parallel_tools
    implicit none
!---------------------------------------------------------------------
! adds periodically left (m=mmm-1) and right (m=mm+1) for cyclic lines
    integer :: nx, ny, nz
    integer :: mmm, mm, n, k
    real(4) :: ff(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz)

    integer, dimension(2) :: p_dist
    integer :: dist_rank
    integer :: ierr
    integer stat(MPI_status_size)

    if (p_coord(1) .eq. 0) then
!-------------- proc has mmm-1 area --------------------------------------------
        p_dist(1) = p_size(1) - 1
        p_dist(2) = p_coord(2)
        call mpi_cart_rank(cart_comm, p_dist,dist_rank,ierr)

        call mpi_sendrecv(ff(mmm, bnd_y1:bnd_y2, 1:nz),         &
                          (bnd_y2 - bnd_y1 + 1)*nz,              &
                          mpi_real4, dist_rank, 1,          &
                          ff(mmm-1, bnd_y1:bnd_y2, 1:nz),       &
                          (bnd_y2 - bnd_y1 + 1)*nz,              &
                          mpi_real4, dist_rank, 1,          &
                          cart_comm, stat, ierr)
!                write(*,*)rank,p_coord,'0',ierr
    endif

    if (p_coord(1) .eq. (p_size(1) - 1)) then
!-------------- proc has mm+1 area ---------------------------------------------
        p_dist(1) = 0
        p_dist(2) = p_coord(2)
        call mpi_cart_rank(cart_comm, p_dist,dist_rank,ierr)

        call mpi_sendrecv(ff(mm, bnd_y1:bnd_y2, 1:nz),          &
                          (bnd_y2 - bnd_y1 + 1)*nz,              &
                          mpi_real4, dist_rank, 1,          &
                          ff(mm+1, bnd_y1:bnd_y2, 1:nz),        &
                          (bnd_y2 - bnd_y1 + 1)*nz,              &
                          mpi_real4, dist_rank, 1,          &
                          cart_comm, stat, ierr)
!                write(*,*)rank,p_coord,'size-1',ierr
    endif
endsubroutine cyclize_x

!======================================================
subroutine cyclize8_x(ff,nx,ny,nz,mmm,mm)
    use mpi_parallel_tools
    implicit none
!---------------------------------------------------------------------
! adds periodically left (m=mmm-1) and right (m=mm+1) for cyclic lines
    integer :: nx, ny, nz
    integer :: mmm, mm, n, k
    real(8) :: ff(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz)

    integer, dimension(2) :: p_dist
    integer :: dist_rank
    integer :: ierr
    integer stat(MPI_status_size)

    if (p_coord(1) .eq. 0) then
!-------------- proc has mmm-1 area --------------------------------------------
        p_dist(1) = p_size(1) - 1
        p_dist(2) = p_coord(2)
        call mpi_cart_rank(cart_comm, p_dist,dist_rank,ierr)

        call mpi_sendrecv(ff(mmm, bnd_y1:bnd_y2, 1:nz),         &
                          (bnd_y2 - bnd_y1 + 1)*nz,              &
                          mpi_real8, dist_rank, 1,          &
                          ff(mmm-1, bnd_y1:bnd_y2, 1:nz),       &
                          (bnd_y2 - bnd_y1 + 1)*nz,              &
                          mpi_real8, dist_rank, 1,          &
                          cart_comm, stat, ierr)
!                write(*,*)rank,p_coord,'0',ierr
    endif

    if (p_coord(1) .eq. (p_size(1) - 1)) then
!-------------- proc has mm+1 area ---------------------------------------------
        p_dist(1) = 0
        p_dist(2) = p_coord(2)
        call mpi_cart_rank(cart_comm, p_dist,dist_rank,ierr)

        call mpi_sendrecv(ff(mm, bnd_y1:bnd_y2, 1:nz),          &
                          (bnd_y2 - bnd_y1 + 1)*nz,              &
                          mpi_real8, dist_rank, 1,          &
                          ff(mm+1, bnd_y1:bnd_y2, 1:nz),        &
                          (bnd_y2 - bnd_y1 + 1)*nz,              &
                          mpi_real8, dist_rank, 1,          &
                          cart_comm, stat, ierr)
!                write(*,*)rank,p_coord,'size-1',ierr
    endif
endsubroutine cyclize8_x

!=====================================================
subroutine cyclize_y(ff,nx,ny,nz,nnn,nn)
implicit none
!---------------------------------------------------------------------
! adds periodically bottom (n=nnn-1) and top (n=nn+1) for cyclic lines
 integer nx, ny, nz
 integer nnn, nn, m, k
 real(4) ff(nx,ny,nz)

!$omp parallel do private(m,k)
  do m=1,nx
   do k=1,nz
      ff(m,nnn-1,k) = ff(m,nn ,k)
      ff(m, nn+1,k) = ff(m,nnn,k)
   enddo
  enddo
!$omp end parallel do

endsubroutine cyclize_y

!=====================================================
subroutine cyclize8_y(ff,nx,ny,nz,nnn,nn)
implicit none
!---------------------------------------------------------------------
! adds periodically bottom (n=nnn-1) and top (n=nn+1) for cyclic lines
 integer nx, ny, nz
 integer nnn, nn, m, k
 real(8) ff(nx,ny,nz)

!$omp parallel do private(m,k)
  do m=1,nx
   do k=1,nz
      ff(m,nnn-1,k) = ff(m,nn ,k)
      ff(m, nn+1,k) = ff(m,nnn,k)
   enddo
  enddo
!$omp end parallel do

endsubroutine cyclize8_y


!=============================================================================
subroutine air_sea_turbulent_fluxes(wnd,        &   ! wind modulo, m/s
                                    slp,        &   ! sea level pressure, Pa
                                    tair,       &   ! air temperature,  �C
                                    tsea,       &   ! sea surface temp, �C
                                    qair,       &   ! air specific humidity, kg/kg
                                    u_hgt,      &   ! height of wind datasets, m
                                    t_hgt,      &   ! height of tair datasets, m
                                    q_hgt,      &   ! height of qair datasets, m
                                    tx,         &   !      zonal wind stress, Pa
                                    ty          )   ! meridional wind stress, Pa
implicit none
  real(8), parameter :: grav    = 9.80d0     !m/s^2 Gravity acceleration
  real(8), parameter :: vonkarm = 0.40d0     !      Von Karman constant
  real(8), parameter :: rgas    = 287.04d0   !J/kg/K Gas constant
  real(8), parameter :: cp_air  = 1000.5d0   !J/kg/K Air heat capacity
  real(8), parameter :: lambda_v= 2.5d6      !J/kg   Heat of evaporation

  real(8) wnd, slp, tair, tsea, qair, u_hgt, t_hgt, q_hgt       !input data
  real(8) tx, ty                !output data

  integer niter, iter
  parameter(niter=2)

  real(8) q1, q2
  parameter(q1=0.98d0*640380.0d0,       &  !kg/m^3
            q2=-5107.4d0)                  !K

  real(8) rho_a, tv, qsat
  real(8) u, u10, t_zu, q_zu, tstar, qstar, z0, xx, stab
  real(8) cd, ch, ce, ustar, bstar
  real(8) cd_n10, ce_n10, ch_n10, cd_n10_rt   ! neutral 10m drag coefficients
  real(8) cd_rt
  real(8) zeta_u, zeta_t, zeta_q, x2, x,  &
              psi_m_u, psi_m_t, psi_m_q,  &
              psi_h_u, psi_h_t, psi_h_q       ! stability parameters

      t_zu = tair
      q_zu = qair

      u = max(wnd, 0.5d0)         ! 0.5 m/s floor on wind
      u10 = u

      cd_n10 = (2.7d0/u10+0.142d0+0.0764d0*u10)/1d3           ! L-Y eqn. 6a
      cd_n10_rt = dsqrt(cd_n10)
      ce_n10 =  34.6d0 *cd_n10_rt/1d3                         ! L-Y eqn. 6b
      stab = 0.5d0 + dsign(0.5d0,t_zu-tsea)
      ch_n10 = (18.0d0*stab+32.7d0*(1.0d0-stab))      &
                                        *cd_n10_rt/1d3        ! L-Y eqn. 6c

! first guess for exchange coeff's at z
      cd    = cd_n10
      ch    = ch_n10
      ce    = ce_n10

      tv= (t_zu+273.15d0)*(1.0d0+0.608d0*q_zu)
      rho_a=slp/(rgas*tv)
      qsat=q1*dexp(q2/(tsea+273.15d0))/rho_a

!     if(nint(u_hgt)/=10.or.nint(t_hgt)/=10.or.nint(q_hgt)/=10) then
      do iter=1,niter

        cd_rt = dsqrt(cd)
        ustar = cd_rt*u                                    ! L-Y eqn. 7a
        tstar = (ch/cd_rt)*(t_zu-tsea)                     ! L-Y eqn. 7b
        qstar = (ce/cd_rt)*(q_zu-qsat)                     ! L-Y eqn. 7c
        bstar    = grav*(tstar/tv+qstar/(q_zu+1.0d0/0.608d0))

! stability for U-height
        zeta_u   = vonkarm*bstar*u_hgt/(ustar*ustar)            ! L-Y eqn. 8a
        zeta_u   = dsign(min(dabs(zeta_u),10.0d0), zeta_u )   ! undocumented NCAR
        x2 = dsqrt(dabs(1.0d0-16.0d0*zeta_u))                 ! L-Y eqn. 8b
        x2 = max(x2, 1.0d0)                                   ! undocumented NCAR
        x = dsqrt(x2);

        if (zeta_u > 0.0d0) then
          psi_m_u = -5.0d0*zeta_u                             ! L-Y eqn. 8c
          psi_h_u = -5.0d0*zeta_u                             ! L-Y eqn. 8c
        else
          psi_m_u = dlog((1.0d0+2.0d0*x+x2)*(1.0d0+x2)/8.0d0)  &
                            -2.0d0*(datan(x)-datan(1.0d0))    ! L-Y eqn. 8d
          psi_h_u = 2.0d0*dlog((1.0d0+x2)/2.0d0)              ! L-Y eqn. 8e
        end if

 ! stability for T-height
        zeta_t   = vonkarm*bstar*t_hgt/(ustar*ustar)          ! L-Y eqn. 8a
        zeta_t   = dsign(min(dabs(zeta_t),10.0d0), zeta_t )   ! undocumented NCAR
        x2 = dsqrt(dabs(1.0d0-16.0d0*zeta_t))                 ! L-Y eqn. 8b
        x2 = max(x2, 1.0d0)                                   ! undocumented NCAR
        x = dsqrt(x2);

        if (zeta_t > 0.0d0) then
          psi_m_t = -5.0d0*zeta_t                             ! L-Y eqn. 8c
          psi_h_t = -5.0d0*zeta_t                             ! L-Y eqn. 8c
        else
          psi_m_t = dlog((1.0d0+2.0d0*x+x2)*(1.0d0+x2)/8.0d0)  &
                            -2.0d0*(datan(x)-datan(1.0d0))    ! L-Y eqn. 8d
          psi_h_t = 2.0d0*dlog((1.0d0+x2)/2.0d0)              ! L-Y eqn. 8e
        end if

! stability for Q-height
        zeta_q   = vonkarm*bstar*q_hgt/(ustar*ustar)            ! L-Y eqn. 8a
        zeta_q   = dsign(min(dabs(zeta_q),10.0d0), zeta_q )   ! undocumented NCAR
        x2 = dsqrt(dabs(1.0d0-16.0d0*zeta_q))                 ! L-Y eqn. 8b
        x2 = max(x2, 1.0d0)                                   ! undocumented NCAR
        x = dsqrt(x2);

        if (zeta_q > 0.0d0) then
          psi_m_q = -5.0d0*zeta_q                             ! L-Y eqn. 8c
          psi_h_q = -5.0d0*zeta_q                             ! L-Y eqn. 8c
        else
          psi_m_q = dlog((1.0d0+2.0d0*x+x2)*(1.0d0+x2)/8.0d0)  &
                            -2.0d0*(datan(x)-datan(1.0d0))    ! L-Y eqn. 8d
          psi_h_q = 2.0d0*dlog((1.0d0+x2)/2.0d0)              ! L-Y eqn. 8e
        end if

        u10 = u/(1.0d0+cd_n10_rt*(dlog(u_hgt/10.0d0)-psi_m_u)/vonkarm) ! L-Y eqn. 9
        t_zu= tair  - tstar/vonkarm*(dlog(t_hgt/u_hgt) +psi_h_u - psi_h_t )
        q_zu= qair  - qstar/vonkarm*(dlog(q_hgt/u_hgt) +psi_h_u - psi_h_q )

        cd_n10 = (2.7d0/u10+0.142d0+0.0764d0*u10)/1d3               ! L-Y eqn. 6a again
        cd_n10_rt = dsqrt(cd_n10)
        ce_n10 = 34.6d0*cd_n10_rt/1d3                               ! L-Y eqn. 6b again
        stab = 0.5d0 + dsign(0.5d0,zeta_t)
        ch_n10 = (18.0d0*stab+32.7d0*(1.0d0-stab))*cd_n10_rt/1d3    ! L-Y eqn. 6c again

        xx = (dlog(u_hgt/10.0d0)-psi_m_u)/vonkarm
        cd = cd_n10/(1.0d0+cd_n10_rt*xx)**2                        ! L-Y 10a

        xx = (dlog(u_hgt/10.0d0)-psi_h_u)/vonkarm
        ch = ch_n10/(1.0d0+ch_n10*xx/cd_n10_rt)*dsqrt(cd/cd_n10) ! 10b (corrected code aug2007)
        ce = ce_n10/(1.0d0+ce_n10*xx/cd_n10_rt)*dsqrt(cd/cd_n10) ! 10c (corrected code aug2007)

        tv= (t_zu+273.15d0)*(1.0d0+0.608d0*q_zu)
        rho_a=slp/(rgas*tv)
        qsat=q1*dexp(q2/(tsea+273.15d0))/rho_a

      end do
!    endif

      tx = rho_a*cd*u
      ty = rho_a*cd*u

endsubroutine air_sea_turbulent_fluxes
