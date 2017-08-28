subroutine ocean_model_parameters(tau)
use main_basin_pars
use mpi_parallel_tools
use basin_grid
use ocean_variables

implicit none
character(256)    t_mask_file,       &  !name of file with temperature point sea-land mask
       bottom_topography_file,       &  !name of file with bottom topography
       help_string
real(4) array4(bnd_x1:bnd_x2,bnd_y1:bnd_y2)
integer  m, n, k, ierr
real(8) tau
real(8) :: hx2, hy2

! define parameters of task
! description of parameters see in file with mame filepar
 if (rank .eq. 0) then
     open (90,file='phys_proc.par',status='old')
     read(90,*) ksw_ts            !Temperature and salinity equation solving (0 - no, 1 - yes)
     read(90,*) ksw_age           !Ideal age equation solving (0 - no, 1 - yes)
     read(90,*) ksw_pt            !Passive tracer equation solving (0 - no, 1 - yes)
     read(90,*) ksw_uv            !Momentum equation solving (0 - no, 1 - yes)
     read(90,*) ksw_lat           !Lateral 2nd order mix parametrization (0 - constant coeff, 1 - Smagorinski)
     read(90,*) ksw_lat4          !Lateral 4nd order momentum mix parametrization (0 - no, 1 - yes)
     read(90,*) ksw_vert          !vertical mix parametrization (0 - constant coeff, 1 - Pacanowski&Philander)
     read(90,*) ksw_dens          !pressure gradient computation (0 - no (constant density), 1 - yes)
     read(90,*) ksw_ice_th        !sea ice thermodynamics using (0 - no, 1 - yes)
     read(90,*) ksw_ice_tran      !sea ice transport using (0 - no, 1 - yes)
     read(90,*) ksw_ice_dyn       !sea ice dynamics using (0 - no, 1 - yes)
     read(90,*) ksw_ssbc          !Type of surface boundary conditions (1 - surface T&S and wind stress are set; 2 - T&S fluxes and wind stress are set; 3 - T&S fluxes and wind stress are computed
     read(90,*) ksw_wflux         !normalize global mean salt balance (0 - no, 1 - normalize water flux, 2 - normalize salinity flux)
     read(90,*) ksw_lbc_ts        !open boundary conditions for T&S using (0 - no, 1 - yes)
     read(90,*) ksw_lbc_uv        !open boundary conditions for U&V using (0 - no, 1 - yes)
     read(90,*) ksw_lbc_ssh       !open boundary conditions for SSH using (0 - no, 1 - yes)

     read(90,*) sst_relax         !Relaxation coefficient for temperature [m/s]
     read(90,*) sss_relax         !Relaxation coefficient for salinity [m/s]
     read(90,*) ldiff_ts          !lateral diffusion for temperature [m**2/s]
     read(90,*) lvisc_2           !lateral  vicosity(2nd order)[m**2/s]
     read(90,*) lvisc_4           !lateral  vicosity(4th order) [undimensional]
     read(90,*) tsfrac_lat        !fraction of salinity lateral diffusion due to one for temperature
     read(90,*) vdiff_ts_min      !vertical background diffusion coefficient for T [m**2/s]
     read(90,*) vdiff_ts_max      !vertical max(top) diffusion coefficient for T [m**2/s]
     read(90,*) vvisc_min         !vertical background viscous coefficient [m**2/s]
     read(90,*) vvisc_max         !vertical max(top) viscous coefficient [m**2/s]
     read(90,*) tsfrac_vert       !fraction of salinity vertical diffusion due to one for temperature
     read(90,*) z_frac            !weight coefficient for lateral Z-Diffusion
     read(90,*) r_frac            !weight coefficient for lateral R-Diffusion
     read(90,*) gm_ratio          !weight coefficient for Gent&McWilliams transport

     help_string =' '
     read (90,'(a)') help_string   ! file with t-mask'
     call get_first_lexeme(help_string ,t_mask_file   )

     help_string =' '
     read (90,'(a)') help_string  ! file with bottom topography'
     call get_first_lexeme(help_string ,bottom_topography_file  )

     close(90)
 endif

 call mpi_bcast(ksw_ts      , 1, mpi_integer, 0, cart_comm, ierr)      !Temperature and salinity equation solving (0 - no, 1 - yes)
 call mpi_bcast(ksw_age     , 1, mpi_integer, 0, cart_comm, ierr)      !Ideal age equation solving (0 - no, 1 - yes)
 call mpi_bcast(ksw_pt      , 1, mpi_integer, 0, cart_comm, ierr)      !Passive tracer equation solving (0 - no, 1 - yes)
 call mpi_bcast(ksw_uv      , 1, mpi_integer, 0, cart_comm, ierr)      !Momentum equation solving (0 - no, 1 - yes)
 call mpi_bcast(ksw_lat     , 1, mpi_integer, 0, cart_comm, ierr)      !Lateral 2nd order mix parametrization (0 - constant coeff, 1 - Smagorinski)
 call mpi_bcast(ksw_lat4    , 1, mpi_integer, 0, cart_comm, ierr)      !Lateral 4nd order momentum mix parametrization (0 - no, 1 - yes)
 call mpi_bcast(ksw_vert    , 1, mpi_integer, 0, cart_comm, ierr)      !vertical mix parametrization (0 - constant coeff, 1 - Pacanowski&Philander)
 call mpi_bcast(ksw_dens    , 1, mpi_integer, 0, cart_comm, ierr)      !pressure gradient computation (0 - no (constant density), 1 - yes)
 call mpi_bcast(ksw_ice_th  , 1, mpi_integer, 0, cart_comm, ierr)      !sea ice thermodynamics using (0 - no, 1 - yes)
 call mpi_bcast(ksw_ice_tran, 1, mpi_integer, 0, cart_comm, ierr)      !sea ice transport using (0 - no, 1 - yes)
 call mpi_bcast(ksw_ice_dyn , 1, mpi_integer, 0, cart_comm, ierr)      !sea ice dynamics using (0 - no, 1 - yes)
 call mpi_bcast(ksw_ssbc    , 1, mpi_integer, 0, cart_comm, ierr)      !Type of surface boundary conditions (1 - surface T&S and wind stress are set; 2 - T&S fluxes and wind stress are set; 3 - T&S fluxes and wind stress are computed
 call mpi_bcast(ksw_wflux   , 1, mpi_integer, 0, cart_comm, ierr)      !normalize global mean salt balance (0 - no, 1 - normalize water flux, 2 - normalize salinity flux)
 call mpi_bcast(ksw_lbc_ts  , 1, mpi_integer, 0, cart_comm, ierr)      !open boundary conditions for T&S using (0 - no, 1 - yes)
 call mpi_bcast(ksw_lbc_uv  , 1, mpi_integer, 0, cart_comm, ierr)      !open boundary conditions for U&V using (0 - no, 1 - yes)
 call mpi_bcast(ksw_lbc_ssh , 1, mpi_integer, 0, cart_comm, ierr)      !open boundary conditions for SSH using (0 - no, 1 - yes)

 call mpi_bcast(sst_relax   , 1, mpi_real8,   0, cart_comm, ierr)      !Relaxation coefficient for temperature [m/s]
 call mpi_bcast(sss_relax   , 1, mpi_real8,   0, cart_comm, ierr)      !Relaxation coefficient for salinity [m/s]
 call mpi_bcast(ldiff_ts    , 1, mpi_real8,   0, cart_comm, ierr)      !lateral diffusion for temperature [m**2/s]
 call mpi_bcast(lvisc_2     , 1, mpi_real8,   0, cart_comm, ierr)      !lateral  vicosity(2nd order)[m**2/s]
 call mpi_bcast(lvisc_4     , 1, mpi_real8,   0, cart_comm, ierr)      !lateral  vicosity(4th order) [undimensional]
 call mpi_bcast(tsfrac_lat  , 1, mpi_real8,   0, cart_comm, ierr)      !fraction of salinity lateral diffusion due to one for temperature
 call mpi_bcast(vdiff_ts_min, 1, mpi_real8,   0, cart_comm, ierr)      !vertical background diffusion coefficient for T [m**2/s]
 call mpi_bcast(vdiff_ts_max, 1, mpi_real8,   0, cart_comm, ierr)      !vertical max(top) diffusion coefficient for T [m**2/s]
 call mpi_bcast(vvisc_min   , 1, mpi_real8,   0, cart_comm, ierr)      !vertical background viscous coefficient [m**2/s]
 call mpi_bcast(vvisc_max   , 1, mpi_real8,   0, cart_comm, ierr)      !vertical max(top) viscous coefficient [m**2/s]
 call mpi_bcast(tsfrac_vert , 1, mpi_real8,   0, cart_comm, ierr)      !fraction of salinity vertical diffusion due to one for temperature
 call mpi_bcast(z_frac      , 1, mpi_real8,   0, cart_comm, ierr)      !weight coefficient for lateral Z-Diffusion
 call mpi_bcast(r_frac      , 1, mpi_real8,   0, cart_comm, ierr)      !weight coefficient for lateral R-Diffusion
 call mpi_bcast(gm_ratio    , 1, mpi_real8,   0, cart_comm, ierr)      !weight coefficient for Gent&McWilliams transport

 call mpi_bcast(t_mask_file, 256, mpi_character, 0, cart_comm, ierr)
 call mpi_bcast(bottom_topography_file, 256, mpi_character, 0, cart_comm, ierr)

 if (rank .eq. 0) then
     write(*,'(i7,a)') ksw_ts,   ' - Temperature and salinity equation solving'
     write(*,'(i7,a)') ksw_age,  ' - Ideal age equation solving'
     write(*,'(i7,a)') ksw_pt,   ' - Passive tracer equation solving'
     write(*,'(i7,a)') ksw_uv,   ' - Momentum equation solving'
     write(*,'(i7,a)') ksw_lat,      ' - Lateral 2nd order mix parametrization'
     write(*,'(i7,a)') ksw_lat4,     ' - Lateral 4nd order momentum mix parametrization'
     write(*,'(i7,a)') ksw_vert,     ' - Vertical mix parametrization'
     write(*,'(i7,a)') ksw_dens,     ' - Pressure gradient computation'
     write(*,'(i7,a)') ksw_ice_th,   ' - Sea ice thermodynamics using'
     write(*,'(i7,a)') ksw_ice_tran, ' - Sea ice transport using'
     write(*,'(i7,a)') ksw_ice_dyn,  ' - Sea ice dynamics using'
     write(*,'(i7,a)') ksw_ssbc,     ' - Type of surface boundary conditions'
     write(*,'(i7,a)') ksw_wflux,    ' - Normalize global mean salt balance'
     write(*,'(i7,a)') ksw_lbc_ts,   ' - Open boundary conditions for T&S using'
     write(*,'(i7,a)') ksw_lbc_uv,   ' - Open boundary conditions for U&V using'
     write(*,'(i7,a)') ksw_lbc_ssh,  ' - Open boundary conditions for SSH using'
     write(*,'(e12.4,a)') sst_relax,     ' - Relaxation coefficient for temperature [m/s]'
     write(*,'(e12.4,a)') sss_relax,     ' - Relaxation coefficient for salinity [m/s]'
     write(*,'(e12.4,a)') ldiff_ts,      ' - Lateral diffusion for temperature [m**2/s]'
     write(*,'(e12.4,a)') lvisc_2,       ' - Lateral  vicosity(2nd order)[m**2/s]'
     write(*,'(e12.4,a)') lvisc_4,       ' - Lateral  vicosity(4th order)[undim]'
     write(*,'(e12.4,a)') tsfrac_lat,    ' - fraction of salinity lateral diffusion'
     write(*,'(e12.4,a)') vdiff_ts_min,  ' - Vertical background diffusion coefficient for T [m**2/s]'
     write(*,'(e12.4,a)') vdiff_ts_max,  ' - Vertical max(top) diffusion coefficient for T [m**2/s]'
     write(*,'(e12.4,a)') vvisc_min,     ' - Vertical background viscosity coefficient [m**2/s]'
     write(*,'(e12.4,a)') vvisc_max,     ' - Vertical max(top) viscosity coefficient [m**2/s]'
     write(*,'(e12.4,a)') tsfrac_vert,   ' - fraction of salinity vertical diffusion'
     write(*,'(e12.4,a)') z_frac,        ' - Weight coefficient for lateral Z-Diffusion'
     write(*,'(e12.4,a)') r_frac,        ' - Weight coefficient for lateral R-Diffusion'
     write(*,'(e12.4,a)') gm_ratio,      ' - Weight coefficient for Gent&McWilliams transport'
     write(*,'(a,a)')  ' file with T-point sea-land mask: ', t_mask_file(1:len_trim (t_mask_file))
     write(*,'(a,a)')  '     file with bottom topography: ', bottom_topography_file(1:len_trim (bottom_topography_file))
 endif

! igrzts_surf  = min(IABS(ksw_ssbc),2) ! type of condition for T and S on sea surface
! igrzts_bot = 2                   ! type of condition for T and S on sea bottom

 ! area mask initialization
   call gridcon(t_mask_file)
   if (rank .eq. 0) print *, "--------------------END OF GRIDCON----------------------"

!  setting vertical t-,w- grid levels
   call vgrid
   if (rank .eq. 0) print *, "--------------------END OF VGRID----------------------"

! define grid geographical coordinates, steps and coriolis parameters
   call basinpar
   if (rank .eq. 0) print *, "--------------------END OF BASINPAR----------------------"

!   hhq_rest = 3000.0d0
   array4=0.0
   call prdstd(' ',bottom_topography_file,1,array4,lu,nx,ny,1, mmm,mm,nnn,nn,1,1,ierr)
   hhq_rest=dble(array4)
   call syncborder_real8(hhq_rest, 1)

   if(periodicity_x/=0) then
       call cyclize8_x(hhq_rest,nx,ny,1,mmm,mm)
   end if

   if(periodicity_y/=0) then
       call cyclize8_y(hhq_rest,nx,ny,1,nnn,nn)
   end if

!--------------Rayleigh friction initialization
!$omp parallel do private(m, n, hx2, hy2)
   do n=ny_start,ny_end
       do m=nx_start,nx_end
           if (lu(m,n)*lu(m+1,n)>0.5 .and. lu(m,n)*lu(m,n+1)>0.5) then
               hx2= ( ((hhq_rest(m+1,n)-hhq_rest(m  ,n))/dxt(m  ,n))**2 * dble(lcu(m  ,n))    &
                    +((hhq_rest(m  ,n)-hhq_rest(m-1,n))/dxt(m-1,n))**2 * dble(lcu(m-1,n)) )/dble(lcu(m,n)+lcu(m-1,n))
               hy2= ( ((hhq_rest(m,n+1)-hhq_rest(m,n  ))/dyt(m,n  ))**2 * dble(lcv(m,n  ))    &
                    +((hhq_rest(m,n  )-hhq_rest(m,n-1))/dyt(m,n-1))**2 * dble(lcv(m,n-1)) )/dble(lcv(m,n)+lcv(m,n-1))

               r_diss(m,n)=r_fric*dsqrt(hx2+hy2)
           endif
       enddo
   enddo
!$omp end parallel do

   call syncborder_real8(r_diss, 1)

   if(periodicity_x/=0) then
       call cyclize8_x( r_diss, nx,ny,1,mmm,mm)
   end if

   if(periodicity_y/=0) then
       call cyclize8_y( r_diss, nx,ny,1,nnn,nn)
   end if

endsubroutine ocean_model_parameters

subroutine test_init
    use main_basin_pars
    use mpi_parallel_tools
    use basin_grid
    use ocean_variables
    implicit none

    integer :: m, n

    ssh_i = -1
    do m = nx_start, nx_end
        do n = ny_start, ny_end
            ssh_i(m, n) = rank
        enddo
    enddo

end subroutine

subroutine sw_test2
    use main_basin_pars
    use mpi_parallel_tools
    use basin_grid
    use ocean_variables

    implicit none
    integer :: m, n, k, ierr
    real*8 :: hx2, hy2
    real*8 :: u0
    real   :: a
    real*8 :: d2r

    d2r= Pi / 180.0d0

    a = 0.0d0
    u0 = 2d0*Pi*RadEarth / (12.0d0*24*60*60)

    ssh_i = 0.0d0
    sshp_i = 0.0d0
    ubrtr_i = 0.0d0
    vbrtr_i = 0.0d0

    ssh_e = 0.0d0
    sshp_e = 0.0d0

    ubrtr_e = 0.0d0
    ubrtrp_e = 0.0d0
    vbrtr_e = 0.d0
    vbrtrp_e = 0.d0

!---------------------Test 2: --------------------------------------------------!
    do n=ny_start, ny_end
      do m=nx_start, nx_end
          if (lcu(m, n) > 0.5) then
              ubrtr_e(m,n) = u0 * (dcos(d2r*geo_lat_u(m, n))*dcos(d2r*a)            &
                  - dcos(d2r*geo_lon_u(m,n))*dsin(d2r*geo_lat_u(m,n))*dsin(d2r*a))
          endif

          if (lcv(m, n) > 0.5) then
              vbrtr_e(m,n) = u0 * dsin(d2r*geo_lon_v(m, n))*dsin(d2r*a)
          endif

          if (lu(m ,n) > 0.5) then
              ssh_e(m,n) = -(1.0d0 / FreeFallAcc)                                 &
                * (RadEarth*EarthAngVel*u0 + 0.5d0*u0*u0)                         &
                  * (( dcos(d2r*geo_lon_t(m,n))*dcos(d2r*geo_lat_t(m,n))*dsin(d2r*a) &
                      + dsin(d2r*geo_lat_t(m,n))*dcos(d2r*a) )**2)
          endif
      enddo
    enddo

    call syncborder_real8(ubrtr_e, 1)
    call syncborder_real8(vbrtr_e, 1)
    call syncborder_real8(ssh_e, 1)

    ubrtrp_e = ubrtr_e
    vbrtrp_e = vbrtr_e
    sshp_e = ssh_e

    if(periodicity_x/=0) then
        call cyclize8_x(ssh_e ,  nx,ny,1,mmm,mm)
        call cyclize8_x(sshp_e,  nx,ny,1,mmm,mm)
        call cyclize8_x( ubrtr_e,nx,ny,1,mmm,mm)
        call cyclize8_x(ubrtrp_e,nx,ny,1,mmm,mm)
        call cyclize8_x( vbrtr_e,nx,ny,1,mmm,mm)
        call cyclize8_x(vbrtrp_e,nx,ny,1,mmm,mm)
    end if

    if(periodicity_y/=0) then
        call cyclize8_y(ssh_e ,  nx,ny,1,nnn,nn)
        call cyclize8_y(sshp_e,  nx,ny,1,nnn,nn)
        call cyclize8_y( ubrtr_e,nx,ny,1,nnn,nn)
        call cyclize8_y(ubrtrp_e,nx,ny,1,nnn,nn)
        call cyclize8_y( vbrtr_e,nx,ny,1,nnn,nn)
        call cyclize8_y(vbrtrp_e,nx,ny,1,nnn,nn)
    end if

    ssh_i = ssh_e
    sshp_i = sshp_e
    ubrtr_i = ubrtrp_e
    vbrtr_i = vbrtrp_e

!initialize depth for internal mode
    call hh_init(hhq, hhqp, hhqn,    &
                 hhu, hhup, hhun,    &
                 hhv, hhvp, hhvn,    &
                 hhh, hhhp, hhhn,    &
                 ssh_i, sshp_i, hhq_rest)
!initialize depth for external mode
    call hh_init(hhq_e, hhqp_e, hhqn_e,    &
                 hhu_e, hhup_e, hhun_e,    &
                 hhv_e, hhvp_e, hhvn_e,    &
                 hhh_e, hhhp_e, hhhn_e,    &
                 ssh_e, sshp_e, hhq_rest)

endsubroutine sw_test2

subroutine sw_only_inicond(flag_init, path2ocp)
    use main_basin_pars
    use mpi_parallel_tools
    use basin_grid
    use ocean_variables

    implicit none
    integer :: flag_init
    character*(*) path2ocp

    integer :: ierr
    real(4) array4(bnd_x1:bnd_x2, bnd_y1:bnd_y2)

    ubrtr_e = 0.0d0
    vbrtr_e = 0.d0

! Read init sea level
    if (flag_init > 0) then
        if (rank.eq.0) print *, "Read init sea level"
        array4 = 0.0
        call prdstd(path2ocp, 'slf.dat', 1, array4, lu,nx,ny,1, mmm,mm,nnn,nn,1,1,ierr)
        ssh_e = dble(array4)
    else
        if (rank.eq.0) print *, "Init sea level is zero"
        ssh_e = 0.0d0
    endif

    call syncborder_real8(ubrtr_e, 1)
    call syncborder_real8(vbrtr_e, 1)
    call syncborder_real8(ssh_e, 1)

    ubrtrp_e = ubrtr_e
    vbrtrp_e = vbrtr_e
    sshp_e = ssh_e
!    ssh_err = ssh_e

    if(periodicity_x/=0) then
        call cyclize8_x(ssh_e ,  nx,ny,1,mmm,mm)
        call cyclize8_x(sshp_e,  nx,ny,1,mmm,mm)
        call cyclize8_x( ubrtr_e,nx,ny,1,mmm,mm)
        call cyclize8_x(ubrtrp_e,nx,ny,1,mmm,mm)
        call cyclize8_x( vbrtr_e,nx,ny,1,mmm,mm)
        call cyclize8_x(vbrtrp_e,nx,ny,1,mmm,mm)
    end if

    if(periodicity_y/=0) then
        call cyclize8_y(ssh_e ,  nx,ny,1,nnn,nn)
        call cyclize8_y(sshp_e,  nx,ny,1,nnn,nn)
        call cyclize8_y( ubrtr_e,nx,ny,1,nnn,nn)
        call cyclize8_y(ubrtrp_e,nx,ny,1,nnn,nn)
        call cyclize8_y( vbrtr_e,nx,ny,1,nnn,nn)
        call cyclize8_y(vbrtrp_e,nx,ny,1,nnn,nn)
    end if

    ssh_i = ssh_e
    sshp_i = sshp_e
    ubrtr_i = ubrtrp_e
    vbrtr_i = vbrtrp_e

    !initialize depth for internal mode
    call hh_init(hhq, hhqp, hhqn,    &
               hhu, hhup, hhun,    &
               hhv, hhvp, hhvn,    &
               hhh, hhhp, hhhn,    &
               ssh_i, sshp_i, hhq_rest)
    !initialize depth for external mode
    call hh_init(hhq_e, hhqp_e, hhqn_e,    &
               hhu_e, hhup_e, hhun_e,    &
               hhv_e, hhvp_e, hhvn_e,    &
               hhh_e, hhhp_e, hhhn_e,    &
               ssh_e, sshp_e, hhq_rest)


endsubroutine sw_only_inicond
