subroutine ocean_model_step(tau,nstep)
use main_basin_pars
use mpi_parallel_tools
use basin_grid
use ocean_variables
use ocean_bc

implicit none
integer m, n, k, nstep, mark
real(8) density, tau
real(8) time_count

! vertical diffusion coefficients
      if(iabs(ksw_vert)==1) then
! philander-pacanowski mixing
          call richnum(den_pot,uu,vv,rit)
          call ppmix(rit,anzu,vvisc_max,vvisc_min,anzt,vdiff_ts_max,vdiff_ts_min)
      end if

      if(iabs(ksw_vert)==2) then
! Mellor-Yamada mixing
         call turb_gen_diss(q2,     &
                           q2p,     &
                           q2l,     &
                          q2lp,     &
                            uu,     &
                            vv,     &
                       den_pot,     &
                          anzt,     &
                          anzu,     &
                     vvisc_min,     &
                  vdiff_ts_min,     &
                        rhs_q2,     &
                        rhs_q2l     )
      endif

! Computing stress components
      call stress_components(uup,vvp,stress_t,stress_s,nz, 1, 0)
      call stress_components(uup2d,vvp2d,stress_t2d,stress_s2d,1, 1, 0)

      if(ksw_lat>0) then
       call smagorinsky_coeff(lvisc_2, ldiff_ts, stress_t, stress_s, amts, amuv)
       call depth_ave(amuv ,amuv2d ,lu,0)
      end if

      if(z_frac>0.0d0.or.r_frac>0) then
!      calculating isopycnal slopes
       call diffusion_slopes(den_pot,slrx,slry,slzx,slzy)
	end if

!Computation of sea surface boundary conditions
      if(ksw_ssbc>0) then
       call sea_surface_fluxes
      endif
!Computing bottom stresses
      if(type_fric>0) then
       call sea_bottom_fluxes
      endif

!--------------- velocity module-------------------------------

if(ksw_uv>0) then

 RHSx3d=0.0d0; RHSy3d=0.0d0
 RHSx2d_diff_disp =0.0d0; RHSy2d_diff_disp =0.0d0
!computing advective terms for 3d-velocity
 call uv_trans( uu, vv, r_vort,   &
             hhq, hhu, hhv, hhh,          &
         RHSx3d_tran, RHSy3d_tran, nz,    &
         1, 0)

!computing advective terms for 2d-velocity
 call uv_trans( uu2d, vv2d, r_vort2d,   &
                 hhq, hhu, hhv, hhh,            &
    RHSx2d_tran_disp, RHSy2d_tran_disp, 1,     &
    1, 0)

!computing viscous terms for 3d-velocity
 call uv_diff2( amuv, stress_t, stress_s,  &
               hhq, hhu, hhv, hhh,        &
               RHSx3d, RHSy3d, nz,    &
               0)

!computing viscous terms for 2d-velocity
 call uv_diff2( amuv2d, stress_t2d, stress_s2d,  &
                hhq, hhu, hhv, hhh,        &
                RHSx2d_diff_disp, RHSy2d_diff_disp, 1, &
                0)

if(ksw_lat4>0) then
 call uv_diff4( amuv4, stress_t, stress_s,        &
                xxt, yyt, hhq, hhu, hhv, hhh,     &
                RHSx3d, RHSy3d, nz     )
 call uv_diff4( amuv42d, stress_t2d, stress_s2d,  &
                xxt, yyt, hhq, hhu, hhv, hhh,     &
                RHSx2d_diff_disp, RHSy2d_diff_disp, 1 )
endif

!Computing pressure gradients
  if(ksw_dens>0) then

     call pressure_gradients(den,        &
                             RHSx3d,     &
                             RHSy3d)
  endif

   call depth_ave(RHSx3d,RHSx2d,llu,0) !computing barotrop. comp from zonal RHS
   call depth_ave(RHSy3d,RHSy2d,llv,0) !computing barotrop. comp from meridional RHS

   call depth_ave(RHSx3d_tran,RHSx2d_tran,llu,0) !computing barotrop. comp from zonal RHS
   call depth_ave(RHSy3d_tran,RHSy2d_tran,llv,0) !computing barotrop. comp from meridional RHS

!add surface and bottom stresses to 2d RHS
 !$omp parallel do private(m,n,k)
      do n=ny_start,ny_end
       do m=nx_start,nx_end

        if(lcu(m,n)>0.5) then
            RHSx2d(m,n)=RHSx2d(m,n)     -RHSx2d_diff_disp(m,n)               &
                       +RHSx2d_tran(m,n)-RHSx2d_tran_disp(m,n)               &
                       +( surf_stress_x(m,n)+bot_stress_x(m,n) )*dxt(m,n)*dyh(m,n)    &
                       -(slpr(m+1,n)-slpr(m,n))*hhu(m,n)*dyh(m,n)/RefDen
        endif

        if(lcv(m,n)>0.5) then
            RHSy2d(m,n)=RHSy2d(m,n)     -RHSy2d_diff_disp(m,n)               &
                       +RHSy2d_tran(m,n)-RHSy2d_tran_disp(m,n)               &
                       +( surf_stress_y(m,n)+bot_stress_y(m,n) )*dyt(m,n)*dxh(m,n)    &
                       -(slpr(m,n+1)-slpr(m,n))*hhv(m,n)*dxh(m,n)/RefDen
        endif

       end do
	end do
!$omp end parallel do

!    call syncborder_real8(RHSx2d, 1)
!    call syncborder_real8(RHSy2d, 1)

   call start_timer(time_count)
   !computing 2d fast gravity waves in external mode and time-mean internal characteristics
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
                  RHSy2d_diff_disp  )
  call end_timer(time_count)
  if (rank .eq. 0) print *, "Barotropic time: ", time_count
  time_barotrop = time_barotrop + time_count

! removing barotropic component from 3d velocity
!$omp parallel do	private(m,n,k)
  do n=ny_start-1,ny_end+1
   do m=nx_start-1,nx_end+1

    if(lcu(m,n)>0.5) then
      uu(m,n,1:nz)=uu(m,n,1:nz)- uu2d(m,n)
    endif

    if(lcv(m,n)>0.5) then
      vv(m,n,1:nz)=vv(m,n,1:nz)- vv2d(m,n)
    endif

   enddo
  enddo
!$omp end parallel do

  call vertical_velocity(uu,vv,ww,hhu,hhv,wf_tot)

! computing corrected 2d velocity

  !$omp parallel do private(m,n)
  do n=ny_start-1,ny_end+1
   do m=nx_start-1,nx_end+1

    if(lcu(m,n)>0.5) then
      uu2d(m,n)=ubrtr_i(m,n)/hhu(m,n)
      pgrx(m,n)=pgrx(m,n)-(slpr(m+1,n)-slpr(m,n))*hhu(m,n)*dyh(m,n)/RefDen
    endif

    if(lcv(m,n)>0.5) then
      vv2d(m,n)=vbrtr_i(m,n)/hhv(m,n)
      pgry(m,n)=pgry(m,n)-(slpr(m,n+1)-slpr(m,n))*hhv(m,n)*dxh(m,n)/RefDen
    endif

   enddo
  enddo
!$omp end parallel do

  !computing internal SSH from continuity equation using internal 2d velocity
   call ssh_internal(tau,  &
                   ssh_i,  &
                  sshp_i,  &
                 ubrtr_i,  &
                 vbrtr_i,  &
                  wf_tot )

! correcting barotropic component of 3d velocity
!$omp parallel do	private(m,n,k)
  do n=ny_start-1,ny_end+1
   do m=nx_start-1,nx_end+1

    if(lcu(m,n)>0.5) then
      uu(m,n,1:nz)=uu(m,n,1:nz)+ uu2d(m,n)
    endif

    if(lcv(m,n)>0.5) then
      vv(m,n,1:nz)=vv(m,n,1:nz)+ vv2d(m,n)
    endif

   enddo
  enddo
!$omp end parallel do
endif

if(ksw_vert>1) then
!transport-diffusion of turbulent characteristics
  call turb_tran_diff(q2,     &
                     q2p,     &
                     tau,     &
                      uu,     &
                      vv,     &
                      ww,     &
                    amuv,     &
           tur_factor_nu,     &
                    anzu,     &
           tur_factor_nu,     &
                 q2_surf,     &
                  q2_bot,     &
                   rhs_q2    )

  call turb_tran_diff(q2l,    &
                     q2lp,    &
                     tau,     &
                      uu,     &
                      vv,     &
                      ww,     &
                    amuv,     &
           tur_factor_nu,     &
                    anzu,     &
           tur_factor_nu,     &
                q2l_surf,     &
                 q2l_bot,     &
                 rhs_q2l    )
endif

if(ksw_ts>0) then
      call start_timer(time_count)
!----------temperature computation------------------------------
      call tracer_tran_diff(tt,     &
                           ttp,     &
                           tau,     &
                            uu,     &
                            vv,     &
                            ww,     &
                          amts,     &
                         1.0d0,     &
                          anzt,     &
                         1.0d0,     &
                          slrx,     &
                          slry,     &
                          slzx,     &
                          slzy,     &
                        z_frac,     &
                        r_frac,     &
                      gm_ratio,     &
                    tflux_surf,     &
                     tflux_bot,     &
                   igrzts_surf,     &
                    igrzts_bot,     &
                      flux_tem_x,   &
                      flux_tem_y,   &
                          sw_bal,   &
                        ksw_lbc_ts, &
                      numb_of_lqp,  &
                           lqpx,    &
                           lqpy,    &
                       index_of_lb, &
                          tlqbw,    &
                      divswrad,     &
               1.0d0/(HeatCapWater*RefDen))
      call end_timer(time_count)
      time_tracer_tt = time_tracer_tt + time_count

      call start_timer(time_count)
!----------salinity computation------------------------------
      call tracer_tran_diff(ss,     &
                           ssp,     &
                           tau,     &
                            uu,     &
                            vv,     &
                            ww,     &
                          amts,     &
                    tsfrac_lat,     &
                          anzt,     &
                   tsfrac_vert,     &
                          slrx,     &
                          slry,     &
                          slzx,     &
                          slzy,     &
                        z_frac,     &
                        r_frac,     &
                      gm_ratio,     &
                    sflux_surf,     &
                     sflux_bot,     &
                   igrzts_surf,     &
                    igrzts_bot,     &
                      flux_sal_x,   &
                      flux_sal_y,   &
                          sw_bal,   &
                        ksw_lbc_ts, &
                      numb_of_lqp,  &
                           lqpx,    &
                           lqpy,    &
                       index_of_lb, &
                          slqbw,    &
                      divswrad,     &
                          0.0d0)
      call end_timer(time_count)
      time_tracer_ss = time_tracer_ss + time_count
endif

if(ksw_uv>0) then
!correcting advective terms for 3d-velocity
   call uv_trans( uu, vv,  r_vort, &
             hhq, hhu, hhv, hhh,         &
             RHSx3d_tran, RHSy3d_tran, nz ,  &
             1, 0 )

   call start_timer(time_count)
  !solving full equations for 3d horizontal velocity
   call baroclinic_dynamics(tau,          &
                            uu,           &
                            uup,          &
                            vv,           &
                            vvp,          &
                            ww,           &
                            anzu,         &
                            RHSx3d,       &
                            RHSy3d,       &
                            RHSx3d_tran,  &
                            RHSy3d_tran,  &
                                pgrx,     &
                                pgry,     &
                              r_diss,     &
                      surf_stress_x,      &
                      surf_stress_y,      &
                       bot_stress_x,      &
                       bot_stress_y)
   call end_timer(time_count)
   time_baroclin = time_baroclin + time_count
endif !end of velocity block


 !Updating depth functions
 if(full_free_surface>0) then
  call hh_shift(hhq, hhqp, hhqn,   &
                hhu, hhup, hhun,   &
                hhv, hhvp, hhvn,   &
                hhh, hhhp, hhhn, 1, &
                0 )
 endif

!!$omp parallel do private(m,n,k)
!  do n=ny_start-1,ny_end+1
!   do m=nx_start-1,nx_end+1

!    if(lcu(m,n)>0.5) then
!      do k=1,nz
!        uup(m,n,k) = uup(m,n,k)/hhup(m,n)
!      enddo
!    endif

!    if(lcv(m,n)>0.5) then
!      do k=1,nz
!        vvp(m,n,k) = vvp(m,n,k)/hhvp(m,n)
!      enddo
!    endif

!    if(lu(m,n)>0.5) then
!      do k=1,nz
!        ttp(m,n,k) = ttp(m,n,k)/hhqp(m,n)
!        ssp(m,n,k) = ssp(m,n,k)/hhqp(m,n)
!      enddo
!    endif

!  enddo
! enddo
!!$omp end parallel do


!if(ksw_vert>1) then
!!$omp parallel do private(m,n,k)
!  do n=ny_start-1,ny_end+1
!   do m=nx_start-1,nx_end+1

!    if(lu(m,n)>0.5) then
!      do k=2,nz
!        q2p(m,n,k) =  q2p(m,n,k)/hhqp(m,n)
!       q2lp(m,n,k) = q2lp(m,n,k)/hhqp(m,n)
!      enddo
!    endif

!   enddo
!  enddo
!!$omp end parallel do
!endif

  !compute depth mean and vertical velocity

  xxt=uu
  yyt=vv

    call depth_ave(xxt,uu2d ,llu,1)
    call depth_ave(yyt,vv2d ,llv,1)

    call depth_ave(uup,uup2d,llu,0)
    call depth_ave(vvp,vvp2d,llv,0)

!-----------------density definition-----------------------------------
      if (ksw_dens>0) then
!$omp parallel do private(m,n,k)
       do n=ny_start-1,ny_end+1
	  do m=nx_start-1,nx_end+1
            if(lu(m,n)>0.5) then
              do k=1,nz
               den(m,n,k)=density(tt(m,n,k),ss(m,n,k),FreeFallAcc*RefDen*hhq(m,n)*z(k))
           den_pot(m,n,k)=density(tt(m,n,k),ss(m,n,k),0.0d0)
              enddo
            endif
          enddo
         enddo
!$omp end parallel do
      endif

	mark=0

!$omp parallel do	private(m,n,k)
	 do n=ny_start,ny_end
        do m=nx_start,nx_end

	   if(lu(m,n)>0.5) then

           if(ssh_i(m,n)<10.0d0.and.ssh_i(m,n)>-10.0d0) then
            continue
           else
            write(*,*) 'in the point m=',m,'n=',n,'ssh_i=',ssh_i(m,n)
	      mark=1
           endif

          do k=1,nz

           if(tt(m,n,k)<100.0d0.and.tt(m,n,k)>-100.0d0) then
            continue
           else
            write(*,*) 'in the point m=',m,'n=',n,'k=',k,'tt=',tt(m,n,k)
	      mark=1
           endif

           if(ss(m,n,k)<100.0d0.and.ss(m,n,k)>-100.0d0) then
            continue
           else
            write(*,*) 'in the point m=',m,'n=',n,'k=',k,'ss=',ss(m,n,k)
	      mark=1
           endif

          enddo
         endif

	   if(lcu(m,n)>0.5) then
          do k=1,nz

           if(uu(m,n,k)<100.0d0.and.uu(m,n,k)>-100.0d0) then
            continue
           else
            write(*,*) 'in the point m=',m,'n=',n,'k=',k,'uu=',uu(m,n,k)
	      mark=1
           endif

          enddo
         endif

	   if(lcv(m,n)>0.5) then
          do k=1,nz

           if(vv(m,n,k)<100.0d0.and.vv(m,n,k)>-100.0d0) then
            continue
           else
            write(*,*) 'in the point m=',m,'n=',n,'k=',k,'vv=',vv(m,n,k)
	      mark=1
           endif

          enddo
         endif

        enddo
	 enddo
!$omp end parallel do

	if(mark>0) stop

endsubroutine ocean_model_step
