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

subroutine cyclize82d_x(ff, mmm,mm)
    use mpi_parallel_tools
    implicit none
!---------------------------------------------------------------------
! adds periodically left (m=mmm-1) and right (m=mm+1) for cyclic lines
    integer :: mmm, mm, n, k
    real(8) :: ff(bnd_x1:bnd_x2, bnd_y1:bnd_y2)

    integer, dimension(2) :: p_dist
    integer :: dist_rank
    integer :: ierr
    integer stat(MPI_status_size)

    if (p_coord(1) .eq. 0) then
!-------------- proc has mmm-1 area --------------------------------------------
        p_dist(1) = p_size(1) - 1
        p_dist(2) = p_coord(2)
        call mpi_cart_rank(cart_comm, p_dist,dist_rank,ierr)

        call mpi_sendrecv(ff(mmm, bnd_y1:bnd_y2),         &
                          (bnd_y2 - bnd_y1 + 1),              &
                          mpi_real8, dist_rank, 1,          &
                          ff(mmm-1, bnd_y1:bnd_y2),       &
                          (bnd_y2 - bnd_y1 + 1),              &
                          mpi_real8, dist_rank, 1,          &
                          cart_comm, stat, ierr)
                write(*,*) rank, p_coord, '0', mmm, ierr
    endif

    if (p_coord(1) .eq. (p_size(1) - 1)) then
!-------------- proc has mm+1 area ---------------------------------------------
        p_dist(1) = 0
        p_dist(2) = p_coord(2)
        call mpi_cart_rank(cart_comm, p_dist,dist_rank,ierr)

        call mpi_sendrecv(ff(mm, bnd_y1:bnd_y2),          &
                          (bnd_y2 - bnd_y1 + 1),              &
                          mpi_real8, dist_rank, 1,          &
                          ff(mm+1, bnd_y1:bnd_y2),        &
                          (bnd_y2 - bnd_y1 + 1),              &
                          mpi_real8, dist_rank, 1,          &
                          cart_comm, stat, ierr)
                write(*,*) rank, p_coord, 'size-1', mm, ierr
    endif
endsubroutine cyclize82d_x

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



!==========================================================================
function density(t,s,p_pa)
implicit none
! SEA WATER DENSITY DEVIATION FROM 1020 [kg/m**3]
! AS FUNCTION OF T[�C](potential),S[PPT],P[Pa]
! By Brydon et al.:"A new approximation of the equation of state for
!                   seawater, suitable for numerical ocean models."
! In: J.Geophis.Res.,v.104,No.C1, p.1537-1540, 1999.
! Pressure profile is subtracted
       real(8)  t,s,p,p_pa,density, den_cgs
       real(8), parameter:: prestompa=1.0d-6

!-2<t<40;0<s<42;0<p<100.

      p=p_pa*prestompa
        den_cgs   =   -2.00920601d-02                              &
!      + (             5.07043d-04*p-5.43283d-07*p*p)              &
      + ( 5.10768d-05-3.69119d-06*p+6.54837d-09*p*p)*t            &
      + ( 8.05999d-04-9.34012d-07*p+1.38777d-09*p*p)*s            &
      + (-7.40849d-06+5.33243d-08*p-1.01563d-10*p*p)*t*t          &
      + (-3.01036d-06+1.75145d-08*p-2.34892d-11*p*p)*t*s          &
      + ( 3.32267d-08-3.25887d-10*p+4.98612d-13*p*p)*t*t*t        &
      + ( 3.21931d-08-1.65849d-10*p+2.17612d-13*p*p)*t*t*s
     density=den_cgs*1000.0d0
endfunction density

!============================================================
  ! calculate a distance in 3-d between points on unit sphere
  function distance(lon1, lat1, lon2, lat2)
  implicit none
    real(8) dpip180           !for degrees to radians convers
    parameter(dpip180=3.1415926535897/180.0d0)
    real(8) lon1,lat1,lon2,lat2
    real(8) distance
    real(8) x1,y1,z1, x2,y2,z2

    z1 = dsin(lat1*dpip180)
    x1 = dcos(lat1*dpip180)*dcos(lon1*dpip180)
    y1 = dcos(lat1*dpip180)*dsin(lon1*dpip180)

    z2 = dsin(lat2*dpip180)
    x2 = dcos(lat2*dpip180)*cos(lon2*dpip180)
    y2 = dcos(lat2*dpip180)*sin(lon2*dpip180)

    distance = dsqrt((z1-z2)**2+(x1-x2)**2+(y1-y2)**2)

  end function distance

!=============================================================================
subroutine air_sea_turbulent_fluxes(wnd,        &   ! wind modulo, m/s
                                    slp,        &   ! sea level pressure, Pa
                                    tair,       &   ! air temperature,  �C
                                    tsea,       &   ! sea surface temp, �C
                                    qair,       &   ! air specific humidity, kg/kg
                                    u_hgt,      &   ! height of wind datasets, m
                                    t_hgt,      &   ! height of tair datasets, m
                                    q_hgt,      &   ! height of qair datasets, m
                                    sens_heat,  &   ! sensible heat flux, W/m^2
                                    evap_rate,  &   ! evaporation rate, kg/m^2/s
                                    lat_heat,   &   ! latent heat flux, W/m^2
                                    tx,         &   !      zonal wind stress, Pa
                                    ty          )   ! meridional wind stress, Pa
implicit none
  real(8), parameter :: grav    = 9.80d0     !m/s^2 Gravity acceleration
  real(8), parameter :: vonkarm = 0.40d0     !      Von Karman constant
  real(8), parameter :: rgas    = 287.04d0   !J/kg/K Gas constant
  real(8), parameter :: cp_air  = 1000.5d0   !J/kg/K Air heat capacity
  real(8), parameter :: lambda_v= 2.5d6      !J/kg   Heat of evaporation

  real(8) wnd, slp, tair, tsea, qair, u_hgt, t_hgt, q_hgt       !input data
  real(8) sens_heat, evap_rate, lat_heat, tx, ty                !output data

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

      sens_heat = cp_air*rho_a*ch*u*(t_zu-tsea)
      evap_rate =        rho_a*ce*u*(q_zu-qsat)
      lat_heat  = evap_rate*lambda_v
      tx = rho_a*cd*u
      ty = rho_a*cd*u

endsubroutine air_sea_turbulent_fluxes

!=============================================================================
subroutine air_ice_turbulent_fluxes(wnd,        &   ! wind modulo, m/s
                                    slp,        &   ! sea level pressure, Pa
                                    tair,       &   ! air temperature,  �C
                                    tsrf,       &   ! ice-snow surface temp, �C
                                    qair,       &   ! air specific humidity, kg/kg
                                    u_hgt,      &   ! height of wind datasets, m
                                    t_hgt,      &   ! height of tair datasets, m
                                    q_hgt,      &   ! height of qair datasets, m
                                    sens_heat,  &   ! sensible heat flux, W/m^2
                                    evap_rate,  &   ! evaporation rate, kg/m^2/s
                                    lat_heat,   &   ! latent heat flux, W/m^2
                                    tx,         &   !      zonal wind stress, Pa
                                    ty          )   ! meridional wind stress, Pa
implicit none
  real(8), parameter :: grav    = 9.80d0     !m/s^2 Gravity acceleration
  real(8), parameter :: vonkarm = 0.40d0     !      Von Karman constant
  real(8), parameter :: rgas    = 287.04d0   !J/kg/K Gas constant
  real(8), parameter :: cp_air  = 1000.5d0   !J/kg/K Air heat capacity
  real(8), parameter :: lambda_s= 2.839d6    !J/kg   Heat of sublimation

  real(8) wnd, slp, tair, tsrf, qair, u_hgt, t_hgt, q_hgt    !input data
  real(8) sens_heat, evap_rate, lat_heat, tx, ty             !output data

  integer niter, iter
  parameter(niter=5)

  real(8) q1, q2
  parameter(q1=11637800.0d0,       &  !kg/m^3
            q2=-5897.8d0)             !K

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

      cd_n10 = 1.63d-3
      ce_n10 = 1.63d-3
      ch_n10 = 1.63d-3

! first guess for exchange coeff's at z
      cd    = cd_n10
      ch    = ch_n10
      ce    = ce_n10

      tv= (t_zu+273.15d0)*(1.0d0+0.608d0*q_zu)
      rho_a=slp/(rgas*tv)
      qsat=q1*dexp(q2/(tsrf+273.15d0))/rho_a

!     if(nint(u_hgt)/=10.or.nint(t_hgt)/=10.or.nint(q_hgt)/=10) then
      do iter=1,niter

        cd_rt = dsqrt(cd)
        ustar = cd_rt*u                                    ! L-Y eqn. 7a
        tstar = (ch/cd_rt)*(t_zu-tsrf)                     ! L-Y eqn. 7b
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

!        cd_n10 = (2.7d0/u10+0.142d0+0.0764d0*u10)/1d3               ! L-Y eqn. 6a again
!        cd_n10_rt = dsqrt(cd_n10)
!        ce_n10 = 34.6d0*cd_n10_rt/1d3                               ! L-Y eqn. 6b again
!        stab = 0.5d0 + dsign(0.5d0,zeta_t)
!        ch_n10 = (18.0d0*stab+32.7d0*(1.0d0-stab))*cd_n10_rt/1d3    ! L-Y eqn. 6c again

        xx = (dlog(u_hgt/10.0d0)-psi_m_u)/vonkarm
        cd = cd_n10/(1.0d0+cd_n10_rt*xx)**2                        ! L-Y 10a

        xx = (dlog(u_hgt/10.0d0)-psi_h_u)/vonkarm
        ch = ch_n10/(1.0d0+ch_n10*xx/cd_n10_rt)*dsqrt(cd/cd_n10) ! 10b (corrected code aug2007)
        ce = ce_n10/(1.0d0+ce_n10*xx/cd_n10_rt)*dsqrt(cd/cd_n10) ! 10c (corrected code aug2007)

        tv= (t_zu+273.15d0)*(1.0d0+0.608d0*q_zu)
        rho_a=slp/(rgas*tv)
        qsat=q1*dexp(q2/(tsrf+273.15d0))/rho_a

      end do
!    endif

      sens_heat = cp_air*rho_a*ch*u*(t_zu-tsrf)
      evap_rate =        rho_a*ce*u*(q_zu-qsat)
      lat_heat  = evap_rate*lambda_s
      tx = rho_a*cd*u
      ty = rho_a*cd*u

endsubroutine air_ice_turbulent_fluxes
!==============================================================================
function freeze_temp(sal,pres)
 implicit none
 real(8) sal, freeze_temp,pres, s, p
    s=max(sal,0.0d0)
    p=pres*0.0001d0
    freeze_temp= (-0.0575d0 + 1.710523d-3 * dsqrt(s) - 2.154996d-4 * s) *s - 7.53d-4 * p
endfunction freeze_temp

!======================================================================
! module contains the sweep procedure needed for different cases
subroutine factor8(nd,a,b,c,eta,rksi,ii,jj)
implicit none
!----------------------------------------------------------------------
! up-down : a(n)*f(n-1) + b(n)*f(n) + c(n)*f(n+1) = eta(n) ; n=ii,jj
! a(ii) and c(jj) are not used, so they may be set = 0.
! optimized by dianski n.a.
      integer ii, jj, nd
      real(8) a(nd),b(nd),c(nd), eta(nd), rksi(nd)
      real(8), allocatable:: x(:), y(:)
      integer j1, jj1
      real(8) scr

      allocate(x(nd), y(nd))

      x(ii) = - c(ii) / b(ii)
      y(ii) = eta(ii) / b(ii)
      jj1 = ii+1
      do  j1=jj1,jj
       scr=1.0/(a(j1) * x(j1-1) + b(j1))
       x(j1) = -c(j1) * scr
       y(j1) = (eta(j1) - a(j1)*y(j1-1)) * scr
      end do
      rksi(jj)=y(jj)
      jj1 = jj - 1
      do j1=jj1,ii,-1
       rksi(j1) = x(j1) * rksi(j1+1) + y(j1)
      end do

      deallocate(y, x)

endsubroutine factor8

!========program for interpolation from z-horizons to s-levels========
subroutine z2s(f1,      &     !input  3d field (on z-levels)
               f2,      &     !output 3d field (on s-levels)
              hhq,      &     !2d field of bottom topography (in metres)
              msk,      &     !temperature mask
           zsigma,      &     !array of s-levels
             zlvl,      &     !array of z-levels (in metres)
               nx,      &     !dimension on x
               ny,      &     !dimension on y
               nz,      &     !number of s-levels
             nlvl,      &     !number of z-levels
             lev1,      &     !parameter of task
             over,      &     !undefined value
             ierr)            !error index

! if lev1=1 then first level of f1 is definited as fist level of f1
! if lev1=0 then first level of f1 is definited correspondly to its level position

integer nx, ny, nz, nlvl

real(4) f2(nx,ny,nz),         &
        f1(nx,ny,nlvl),       &
       hhq(nx,ny),            &
       msk(nx,ny),            &
       zsigma(nz),            &
         zlvl(nlvl)
real(4) fz1,fz2,z1,z2,zs(150),fz(150),ds,deep
real(4) over, over5
integer i,j,k,kdeep,kups,nocegrid,lev1,ierr, kr

! over5 - 5% vicinity to over
  over5=abs(over)-abs(over)/20.0
  ierr=0
  write(*,*)' interpolation from z-levels to s-levels'

  if(nz>150) then
     write(*,*)' error in routine z2s:'
     write(*,*)' number of s-levels is greater than 150'
     ierr=1
     return
  endif

  if(nlvl>150) then
     write(*,*)' error in routine z2s:'
     write(*,*)' number of z-levels is greater than 150'
     ierr=1
     return
  end if

!  interpolation
  nocegrid=0

  do j=1,ny
   do i=1,nx
    if(msk(i,j)>0.5) then
      nocegrid=nocegrid+1
      deep=hhq(i,j)

      do k=1,nz
         zs(k)=zsigma(k)*deep
      end do

      fz(1)=f1(i,j,1)

      if(abs(fz(1))>over5) then
         write(*,*)' error in routine z2s:'
         write(*,'(a,2i4,a)') ' in point ',i,j,' input value of upper level is undefined!'
         ierr=1
         return
      end if

!--------- making profile without bottom-------------
      ierr=0

      do k=2,nlvl
       fz(k)=f1(i,j,k)
       if(abs(fz(k))>over5) then
         fz(k)=fz(k-1)
         ierr=ierr+1
       end if
      end do

      if(ierr==nlvl-1) then
        write(*,*)' warning in routine z2s:'
        write(*,1000) i,j,deep,zlvl(1)
      end if

      ierr=0

!  searching number of upper sigma level
      kups=1

       do while(zs(kups)<zlvl(1).and.kups<nz)
        kups=kups+1
       end do

      kups=kups-1

      if(kups>=1) then

!   for upper sigma levels of f2
         do k=1,kups
           f2(i,j,k)=f1(i,j,1)
         enddo

      else
         if (lev1/=0) then
          f2(i,j,1)=f1(i,j,1)
          kups=1
         end if
      end if


! searching the deepest z level
      kdeep=1

      do while(zlvl(kdeep)<deep.and.kdeep<nlvl)
       kdeep=kdeep+1
      end do

      kr=1
      fz1=fz(1)
      fz2=fz1
      z1=zlvl(1)-zlvl(1)/20.0
      z2=zlvl(1)

      do k=kups+1,nz

       ds=zs(k)

       do while((ds<=z1.or.ds>z2).and.kr<kdeep)
          kr=kr+1
          fz1=fz2
          z1=z2
          fz2=fz(kr)
          z2=zlvl(kr)
       end do

       if(ds>=zlvl(kdeep)) then
         f2(i,j,k)=fz(kdeep)
       else
         f2(i,j,k)=(fz1*(z2-ds)+fz2*(ds-z1))/(z2-z1)
       end if

      end do

      end if
    enddo
 enddo

 write(*,*) ' for control: number of ocean grid in interpolation is ',nocegrid
 return

1000  format(' in point i=',i5,',',' j=',i5,' deep=',f9.2,'m'/    &
            ' there is only one level for interpolation'/         &
            ' at z-level of',f9.2,'m')

endsubroutine z2s

!==============================================================
function pot_temp(tem,sal,pres)
implicit none

! potential temperature = function of T[�C], S[PSU], p[bar]
! s - salinity deviation from 35ppt

real(8) tem, sal, pres
real(8) t,s,p, pot_temp

t=tem
s=sal-35.0d0
p=pres

pot_temp = t - p*(3.6504d-04       + 8.3198d-05*t - 5.4065d-07*t*t     &
                + 4.0274d-09*t*t*t +(1.7439d-05   - 2.9778d-07*t)*s )  &
           - p*p*(8.9309d-07       - 3.1628d-08*t + 2.1987d-10*t*t     &
                - 4.1057d-09*s     -(1.6056d-10    -5.0484d-12*t)*p )

endfunction pot_temp

!==========density for potential temperature computing=========================================
function denz(tem,sal,pres)
! tem - temperature [�C]
! sal - salinity [PSU]
! pres - pressure[bars]
! denz - density [kg/m**3]
implicit none
real(8), parameter:: rh0=1.02d0, proc=1.00005d0

! polinom coefficients for the fresh water state equation [*1e-03]
real(8), parameter:: fw0= 0.999842594d0, fw1= 6.793952d-05, fw2=-9.095290d-06,     &
                     fw3= 1.001685d-07,  fw4=-1.120083d-09, fw5= 6.536332d-12

! polinom coefficients for the sea water state equation [*1e-03]
real(8), parameter:: sw10= 0.824493d-03, sw11=-4.089900d-06, sw12= 7.643800d-08,   &
                     sw13=-8.246700d-10, sw14= 5.387500d-12, sw20=-5.724660d-06,   &
                     sw21= 1.022700d-07, sw22=-1.654600d-09, sw31= 4.831400d-07

! polinom coefficients for the fresh water Young module
real(8), parameter:: pfw0= 19652.21d0,   pfw1= 148.42060d0,  pfw2=-2.327105d0,     &
                     pfw3= 1.360477d-02, pfw4=-5.155288d-05

! polinom coefficients for the sea water Young module
real(8), parameter:: psw10= 54.67460d0,   psw11=-0.603459d0,   psw12= 1.099870d-02,    &
                     psw13=-6.167000d-05, psw20= 7.944000d-02, psw21= 1.648300d-02,    &
                     psw22=-5.300900d-04

!polinom coefficients for the full module junga
real(8), parameter:: pr10= 3.239908d0,   pr11= 1.437130d-03, pr12= 1.160920d-04,      &
                     pr13=-5.779050d-07, pr20= 2.283800d-03, pr21=-1.098100d-05,      &
                     pr22=-1.607800d-06, pr30=-9.934800d-07, pr31= 2.081600d-08,      &
                     pr32= 9.169700d-10, pr40= 1.910750d-04

! sea water density(+-0.00005kg/m**3)=function of T(�C),S(PSU)& Z(m)
real(8) t,s,pressure, rho, ppym, ss
real(8) tem, sal, pres, denz, den_cgs

t=tem
s=max(sal,0.0d0)
pres=pressure

! density of sea water (p=0) [cm/s**3]
ss = dsqrt(s)
rho =     fw0+ fw1*t + fw2*t**2  + fw3*t**3  + fw4*t**4 + fw5*t**5      &
   +s*(sw10 + sw11*t + sw12*t**2 + sw13*t**3 +sw14*t**4                 &
   + ( sw20 + sw21*t + sw22*t**2 )*ss       + sw31*s )

  if(pressure>0.0d0) then
! pressure (in bars=10**5 n/m**2) per Young's module
  ppym=pressure/( pfw0 + pfw1*t  + pfw2*t**2  + pfw3*t**3 + pfw4*t**4     &      !for fresh water(p=0)
             +s*( psw10+ psw11*t + psw12*t**2 + psw13*t**3                &      !for sea   water(p=0)
             +  ( psw20+ psw21*t + psw22*t**2 )*ss       )                &      !for sea   water(p=0)
   +pressure*(    pr10 + pr11*t  + pr12*t**2  + pr13*t**3                 &      !p>0
             +s*( pr20 + pr21*t  + pr22*t**2  + pr40*ss                   &      !p>0
             +  ( pr30 + pr31*t  + pr32*t**2 ) *pressure ) ) )                   !p>0
! density
   den_cgs=rho/(1.0d0-ppym)-rh0      ![gr/cm**3]
  else
   den_cgs=rho-rh0                 ![gr/cm**3]
 endif

   denz=den_cgs*1000.0d0 !output density is in kg/m^3

 endfunction denz

!==============density at surface========================================
function dens(tem,sal)
implicit none

real(8) t, s, dens, tem, sal, den_cgs
!c unesco81 formula for surface density deviation from 1.020[g/cm**3]
!as function of T[�C], S[PSU]

t=tem
s=max(sal, 0.0d0)
den_cgs = -0.020157406d0     + 6.793952d-05*t    - 9.095290d-06*t**2       &
         + 1.001685d-07*t**3 - 1.120083d-09*t**4 + 6.536332d-12*t**5       &
     + s*( 8.414029d-04     -  4.089900d-06*t    + 7.643800d-08*t**2       &
         + 4.831400d-07*(s-35.0d0) -8.246700d-10*t**3 + 5.387500d-12*t**4  &
       + dsqrt(s)*(-5.72466d-06+1.022700d-07*t   -1.654600d-09*t**2) )

dens=den_cgs*1000.0d0   !output density is in kg/m^3

endfunction dens

!======================================================================
subroutine t2th(tt,ss,den,tth)
! making potential temperature from temperature in situ
! input:
! tt - temperature in situ
! ss - salinity
! output:
! tth - potential temperature
! den - potential density
use main_basin_pars
use basin_grid

implicit none

! hhq - depth on t-grid in m
! lu  - t-grid mask
! nx,ny,nz - dimension on longitude, latitude, depth
 real(4) tt(nx,ny,nz),  ss(nx,ny,nz),     &
        den(nx,ny,nz), tth(nx,ny,nz)
 integer m,n,k
 real(8) dens, denz, pot_temp
 real(8) p_in_bar, d, tem, sal, tpot, den_surf

! define initial density
do n=nnn,nn
  do m=mmm,mm
     if (lu(m,n)>0.5) then
      p_in_bar= FreeFallAcc*RefDen*z(1)*hhq(m,n)*1.0d-05
      do k=1,nz-1
       tem=dble(tt(m,n,k))
       sal=dble(ss(m,n,k))
       d=denz(tem,sal,p_in_bar)
       p_in_bar=p_in_bar      &
        +FreeFallAcc*(RefDen+d)*hzt(k+1)*hhq(m,n)*1.0d-05

        tpot= pot_temp (tem,sal,p_in_bar)
        den_surf=dens(tpot,sal)
        tth(m,n,k)=sngl(tpot)
        den(m,n,k )=sngl(den_surf)
       end do

       tem=dble(tt(m,n,nz))
       sal=dble(ss(m,n,nz))
       tpot= pot_temp (tem,sal,p_in_bar)
       den_surf=dens(tpot,sal)
        tth(m,n,nz)=sngl(tpot)
        den(m,n,nz)=sngl(den_surf)
     end if
   end do
 end do

  write(*,*)' Control tests:'
  write(*,*)' 1) testing of formula for potential temperature:'
  write(*,*)'    theta(T=10�C, S=25 PSU, P=1000) = ', pot_temp(10.d0,25.0d0,1000.d0)
  write(*,*)'                    control number =        8,4678516'
  write(*,*)' 2) testing of formula for density (unesco81):'
  write(*,*)'    rho(S=0 PSU, T=5�C) = ', dens(5.0d0, 0.0d0)+RefDen
  write(*,*)'    rho(S=0 PSU, T=5�C, P=0) = ', denz(5.0d0, 0.0d0,0.0d0)+RefDen
  write(*,*)'           control number =      999,96675'
  write(*,*)'    rho(S=35 PSU, T=5�C)     = ', dens(5.0d0, 35.0d0)+RefDen
  write(*,*)'    rho(S=35 PSU, T=5�C, P=0) = ',denz(5.0d0, 35.0d0, 0.0d0)+RefDen
  write(*,*)'           control number =     1027,67547'
  return
endsubroutine t2th
