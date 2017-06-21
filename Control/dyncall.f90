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

!----------------------- Set constant RHS --------------------------------!
!   do n=ny_start,ny_end
!    do m=nx_start,nx_end
!        if(lcu(m,n)>0.5) then
!          RHSx2d(m, n) = -(10.0d0)*hhu(m,n)*dyh(m,n)/RefDen
!        else
!          RHSx2d(m, n) = 0.0d0
!        endif
!
!        if(lcv(m,n)>0.5) then
!          RHSy2d(m, n) = -(10.0d0)*hhv(m,n)*dxh(m,n)/RefDen
!        else
!          RHSy2d(m, n) = 0.0d0
!        endif
!    enddo
!   enddo
!   call syncborder_real8(RHSx2d, 1)
!   call syncborder_real8(RHSy2d, 1)

!------------------------ Init variables: --------------------------------------!
   wf_tot  = 0.0d0
   r_diss  = 0.0d0
!   amuv2d  = 0.0d0
!   amuv42d = 0.0d0
   amuv2d  = lvisc_2
   amuv42d = lvisc_4

   stress_t2d = 0.0d0
   stress_s2d = 0.0d0
   xxt = 0.0d0
   yyt = 0.0d0
   RHSy2d_diff_disp = 0.0d0
   RHSx2d_diff_disp = 0.0d0

   r_vort2d = 0.0d0
   RHSx2d_tran_disp = 0.0d0
   RHSy2d_tran_disp = 0.0d0

!    do n=ny_start,ny_end
!        do m=nx_start,nx_end
!            print *, "m, n," , m, n, "rls", rlh_s(m,n)
!        enddo
!    enddo
!    stop

!    if (rank.eq.0) print *, "Begin swallow water"
!---------------------- Shallow water equ solver -------------------------------!
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
                             pgrx,   &
                             pgry,   &
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
!   if (rank .eq. 0) print *, "SW step time: ", time_count
   time_barotrop = time_barotrop + time_count

   ssh_i = pgrx
!   Computing internal SSH from continuity equation using internal 2d velocity
!    call ssh_internal(tau,  &
!                    ssh_i,  &
!                   sshp_i,  &
!                  ubrtr_i,  &
!                  vbrtr_i,  &
!                   wf_tot )

   do n=ny_start,ny_end
      do m=nx_start,nx_end
          if(lu(m,n)>0.5) then
              if(ssh_i(m,n)<10000.0d0.and.ssh_i(m,n)>-10000.0d0) then
                  continue
              else
                  write(*,*) rank, 'ERROR!!! In the point m=', m, 'n=', n, 'ssh_i=', ssh_i(m,n),   &
                    'step: ', num_step, 'lon: ', geo_lon_t(m, n), 'lat: ', geo_lat_t(m, n)

                  num_step=num_step+1
                  key_time_print=0
                  call model_time_def(   num_step,           &     !step counter,            input
                                          time_step,           &    !time step in seconds,    input
                                          ndays_in_4yr,        &    !integer day distribution in 4-years (49 months)
                                          seconds_of_day,      &    !current seconds in day  ,output
                                          m_sec_of_min,        &    !second counter in minute,output
                                          m_min_of_hour,       &    !minute counter in hour  ,output
                                          m_hour_of_day,       &    !hour counter in day     ,output
                                          m_day_of_month,      &    !day counter in month    ,output
                                          m_day_of_year,       &    !day counter in year     ,output
                                          m_day_of_4yr,        &    !day counter in 4-years  ,output
                                          m_month_of_year,     &    !mon counter in year     ,output
                                          m_month_of_4yr,      &    !mon counter in 4-years  ,output
                                          m_year_of_4yr,       &    !year counter in 4yrs    ,output
                                          m_day,               &    !model elapsed day counter starting from zero
                                          m_month,             &    !model elapsed month counter starting from zero
                                          m_year,              &    !year counter            ,output
                                          m_4yr,               &    !counter of 4-yr groups  ,output
                                          m_time_changed,      &    !change indicator of time,output
                                          key_time_print,      &    !key of printing time:0-not,1-print
                                          init_year)                !initial real-time year

                  call  parallel_local_output(path2ocp,  &
                                        nrec_loc+1,  &
                                        year_loc,  &
                                         mon_loc,  &
                                         day_loc,  &
                                        hour_loc,  &
                                         min_loc,  &
                                  loc_data_tstep,  &
                                         yr_type  )
                  stop
                  call mpi_finalize(ierr)
              endif
          endif
      enddo
   enddo


endsubroutine shallow_water_model_step
