!--------------------- SUBROUTINE FOR: -----------------------------------------!
!----------------- Only shallow water solving ----------------------------------!
!-------------------------------------------------------------------------------!
subroutine shallow_water_model_step(tau,nstep)
    use main_basin_pars
    use mpi_parallel_tools
    use basin_grid
    use ocean_variables

    use time_integration

    implicit none
    integer :: m, n, k, nstep, mark, ierr
    integer*8 :: curstep
    real(8) density, tau
    real*8 :: time_count

!------------------------ Init variables ---------------------------------------
    if (atm_forcing_on == 1) then
        !Computation of sea surface boundary conditions
        if(ksw_ssbc > 0) then
            call sea_surface_fluxes
        endif

        !Computing bottom stresses
        !if(type_fric>0) then
        !    call sea_bottom_fluxes
        !endif

        !RHSx2d = ( surf_stress_x(m,n)+bot_stress_x(m,n) )*dxt(m,n)*dyh(m,n)    &
        RHSx2d = (surf_stress_x(m,n))*dxt(m,n)*dyh(m,n)    &
                 -(slpr(m+1,n)-slpr(m,n))*hhu(m,n)*dyh(m,n)/RefDen

        !RHSy2d = ( surf_stress_y(m,n)+bot_stress_y(m,n) )*dyt(m,n)*dxh(m,n)    &
        RHSy2d = (surf_stress_y(m,n))*dyt(m,n)*dxh(m,n)    &
                 -(slpr(m,n+1)-slpr(m,n))*hhv(m,n)*dxh(m,n)/RefDen
    else
        RHSx2d = 0.0d0
        RHSy2d = 0.0d0
        wf_tot = 0.0d0
    endif

    r_diss = 0.0d0
    amuv2d  = lvisc_2
    amuv42d = lvisc_4

!---------------------- Shallow water equ solver -------------------------------
    call start_timer(time_count)
    call barotropic_dynamics(tau,     &
                          nstep,     &
                          ksw_lat4,  &
                          ubrtr_e,   &
                          ubrtrp_e,  &
                          vbrtr_e,   &
                          vbrtrp_e,  &
                          ssh_e,     &
                          sshp_e,    &
                          ubrtr_i,   &
                          vbrtr_i,   &
                             ssh_i,  &
                            RHSx2d,  &
                            RHSy2d,  &
                            wf_tot,  &
                            amuv2d,  &
                           amuv42d,  &
                          r_vort2d,  &
                        stress_t2d,  &
                        stress_s2d,  &
                               xxt,  &
                               yyt,  &
                            r_diss,  &
                  RHSx2d_tran_disp,  &
                  RHSy2d_tran_disp,  &
                  RHSx2d_diff_disp,  &
                  RHSy2d_diff_disp,  &
                  RHSx2d_bfc,        &
                  RHSy2d_bfc)
    call end_timer(time_count)
    time_barotrop = time_barotrop + time_count

    do n=ny_start,ny_end
      do m=nx_start,nx_end
          if(lu(m,n)>0.5) then
              if(ssh_i(m,n)<10000.0d0.and.ssh_i(m,n)>-10000.0d0) then
                  continue
              else
                  write(*,*) rank, 'ERROR!!! In the point m=', m, 'n=', n, 'ssh_i=', ssh_i(m,n),   &
                    'step: ', num_step, 'lon: ', geo_lon_t(m, n), 'lat: ', geo_lat_t(m, n)

                  stop
                  call mpi_finalize(ierr)
              endif
          endif
      enddo
    enddo

endsubroutine shallow_water_model_step
