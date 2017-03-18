!========================================================================
subroutine sea_surface_fluxes
use main_basin_pars
use mpi_parallel_tools
use basin_grid
use ocean_variables
use atm_forcing
use ocean_bc
implicit none

integer m,n,k

real(8) evap_rate
real(8) wf_ave, sf_ave, m_calc, wf, wnd

if(ksw_ssbc==1) then

!$omp parallel do private(m,n)
 do n=ny_start,ny_end
   do m=nx_start,nx_end

     if(lu(m,n)>0.5) then
      tflux_surf(m,n)=dble(sst_obs(m,n))
      sflux_surf(m,n)=dble(sss_obs(m,n))
      sw_bal(m,n)=0.0d0
     endif

     if(lcu(m,n)>0.5) then
      surf_stress_x(m,n)=(taux(m  ,n)*dx(m  ,n)*dy(m  ,n)   &
                         +taux(m+1,n)*dx(m+1,n)*dy(m+1,n))/RefDen/2.0d0/dxt(m,n)/dyh(m,n)
     endif

     if(lcv(m,n)>0.5) then
      surf_stress_y(m,n)=(tauy(m,n  )*dx(m,n  )*dy(m,n  )   &
                         +tauy(m,n+1)*dx(m,n+1)*dy(m,n+1))/RefDen/2.0d0/dxh(m,n)/dyt(m,n)
     endif

   enddo
 enddo
!$omp end parallel do

endif

if(ksw_ssbc==2) then

!$omp parallel do
 do n=ny_start,ny_end
   do m=nx_start,nx_end

     if(lu(m,n)>0.5) then
      tflux_surf(m,n)=(hf_tot(m,n)-sw_bal(m,n))/(HeatCapWater*RefDen)   &
                     +dkft(m,n)*(dble(sst_obs(m,n))-ttp(m,n,1))
      sflux_surf(m,n)=-wf_tot(m,n)*ssp(m,n,1)/RefDen                    &
                     +dkfs(m,n)*(dble(sss_obs(m,n))-ssp(m,n,1))
     endif

     if(lcu(m,n)>0.5) then
      surf_stress_x(m,n)=(taux(m  ,n)*dx(m  ,n)*dy(m  ,n)   &
                         +taux(m+1,n)*dx(m+1,n)*dy(m+1,n))/RefDen/2.0d0/dxt(m,n)/dyh(m,n)
     endif

     if(lcv(m,n)>0.5) then
      surf_stress_y(m,n)=(tauy(m,n  )*dx(m,n  )*dy(m,n  )   &
                         +tauy(m,n+1)*dx(m,n+1)*dy(m,n+1))/RefDen/2.0d0/dxh(m,n)/dyt(m,n)
     endif

   enddo
 enddo
!$omp end parallel do

endif

if(ksw_ssbc==3) then
!$omp parallel do private(m, n, evap_rate)
 do n=ny_start,ny_end
   do m=nx_start,nx_end
     if(lu(m,n)>0.5) then

      call air_sea_turbulent_fluxes(wind(m,n),         &   ! wind modulo, m/s
                                    slpr(m,n),         &   ! sea level pressure, Pa
                                    tatm(m,n),         &   ! air temperature,  �C
                                   ttp(m,n,1),         &   ! sea surface temp, �C
                                    qatm(m,n),         &   ! air specific humidity, kg/kg
                                     u_height,         &   ! height of wind datasets, m
                                     t_height,         &   ! height of tair datasets, m
                                     q_height,         &   ! height of qair datasets, m
                                    sensheat(m,n),     &   ! sensible heat flux, W/m^2
                                    evap_rate,         &   ! evaporation rate, kg/m^2/s
                                     latheat(m,n),     &   ! latent heat flux, W/m^2
                                    taux(m,n),         &   !      zonal wind stress, Pa
                                    tauy(m,n)          )   ! meridional wind stress, Pa

      lw_bal(m,n)=EmissWater*(lwr(m,n)-StephBoltz*(ttp(m,n,1)+273.15d0)**4)
      sw_bal(m,n)=swr(m,n)*(1.0d0-AlbOpWater)
      hf_tot(m,n)=sensheat(m,n)+latheat(m,n)+lw_bal(m,n)+sw_bal(m,n)-snow(m,n)*lambda_f

      wf_tot(m,n)=rain(m,n)+snow(m,n)+dble(runoff(m,n))+evap_rate

     endif
   enddo
 enddo
!$omp end parallel do

  call syncborder_real8(wf_tot, 1)

  call syncborder_real8(taux, 1)
  call syncborder_real8(tauy, 1)

  if(periodicity_x/=0) then
    call cyclize8_x(taux,nx,ny,1,mmm,mm)
  endif

  if(periodicity_y/=0) then
    call cyclize8_y(tauy,nx,ny,1,nnn,nn)
  endif

  wf_ave=0.0d0
  sf_ave=0.0d0
  m_calc=0.0d0

  if(ksw_wflux>0) then

  !$omp parallel do private(m,n) reduction(+:wf_ave, m_calc)
   do n=ny_start,ny_end
     do m=nx_start,nx_end
       if(lu(m,n)>0.5) then
         wf_ave=wf_ave+wf_tot(m,n)*dx(m,n)*dy(m,n)
         m_calc=m_calc+dx(m,n)*dy(m,n)
       endif
     enddo
   enddo
  !$omp end parallel do

   wf_tot=wf_tot-wf_ave/m_calc
  endif

  !$omp parallel do private(m, n, k) reduction(+:sf_ave)
   do n=ny_start,ny_end
     do m=nx_start,nx_end
      if(lu(m,n)>0.5) then
       sflux_surf(m,n)= dkfs(m,n)*(dble(sss_obs(m,n))-ssp(m,n,1))
       tflux_surf(m,n)=(hf_tot(m,n)-sw_bal(m,n))/(HeatCapWater*RefDen)   &
                      +dkft(m,n)*(dble(sst_obs(m,n))-ttp(m,n,1))         &
                      +wf_tot(m,n)/RefDen*tatm(m,n)
       sf_ave=sf_ave+sflux_surf(m,n)*dx(m,n)*dy(m,n)
      endif
     enddo
   enddo
  !$omp end parallel do

   if(ksw_wflux>1) then
    sflux_surf=sflux_surf-sf_ave/m_calc
   endif

!$omp parallel do private(m,n,wf,wnd)
 do n=ny_start,ny_end
   do m=nx_start,nx_end
     if(lcu(m,n)>0.5) then
      wnd=(uwnd(m  ,n)*dx(m  ,n)*dy(m  ,n)         &
         + uwnd(m+1,n)*dx(m+1,n)*dy(m+1,n))/2.0d0/dxt(m,n)/dyh(m,n)
      wf=(wf_tot(m  ,n)*uwnd(m  ,n)*dx(m  ,n)*dy(m  ,n)         &
         +wf_tot(m+1,n)*uwnd(m+1,n)*dx(m+1,n)*dy(m+1,n))/2.0d0/dxt(m,n)/dyh(m,n)
      surf_stress_x(m,n)=(taux(m,n)+taux(m+1,n))/RefDen/2.0d0     &
                       *( wnd -uup(m,n,1) )  + wf/RefDen
     endif

     if(lcv(m,n)>0.5) then
      wnd=(vwnd(m,n  )*dx(m,n  )*dy(m,n  )         &
         + vwnd(m,n+1)*dx(m,n+1)*dy(m,n+1))/2.0d0/dxh(m,n)/dyt(m,n)
      wf=(wf_tot(m,n  )*vwnd(m,n  )*dx(m,n  )*dy(m,n  )         &
         +wf_tot(m,n+1)*vwnd(m,n+1)*dx(m,n+1)*dy(m,n+1))/2.0d0/dxh(m,n)/dyt(m,n)
      surf_stress_y(m,n)=(tauy(m,n)+tauy(m,n+1))/RefDen/2.0d0     &
                       *(wnd -vvp(m,n,1) )    + wf/RefDen
     endif

   enddo
 enddo
 !$omp end parallel do

endif

!Boundary condition for turbulent kinetic energy
if(iabs(ksw_vert)>1) then
  !$omp parallel do private(m,n)
   do n=ny_start,ny_end
     do m=nx_start,nx_end
      if(lu(m,n)>0.5) then
        q2_surf(m,n)= dsqrt( ( surf_stress_x(m  ,n  )**2*dxt(m  ,n  )*dyh(m  ,n  )     &
                              +surf_stress_x(m-1,n  )**2*dxt(m-1,n  )*dyh(m-1,n  )     &
                              +surf_stress_y(m  ,n  )**2*dxh(m  ,n  )*dyt(m  ,n  )     &
                              +surf_stress_y(m  ,n-1)**2*dxh(m  ,n-1)*dyt(m  ,n-1) )/2.0d0/dx(m,n)/dy(m,n) )      &
                              *B1_t**(2.0d0/3.0d0)
      endif
     enddo
   enddo
  !$omp end parallel do
endif

endsubroutine sea_surface_fluxes
!==============================================================================================
subroutine sea_bottom_fluxes
use main_basin_pars
use mpi_parallel_tools
use basin_grid
use ocean_variables

implicit none
!----------------------------------------------------------------------
integer m, n
real(8) ut, vt

   if (type_fric==1 ) then
       BottomFriction = Cb_l

! Nonlinear bottom friction
       elseif (type_fric==2) then
!$omp parallel do private(m,n,ut,vt)
        do n=ny_start,ny_end
         do m=nx_start,nx_end
          if(lu(m,n)>0.5) then
           ut = (uu(m  ,n,nz)*dxt(m  ,n)*dyh(m  ,n)*hhu(m  ,n)     &
               + uu(m-1,n,nz)*dxt(m-1,n)*dyh(m-1,n)*hhu(m-1,n)) /2.0d0/dx(m,n)/dy(m,n)/hhq(m,n)
           vt = (vv(m,n  ,nz)*dxh(m,n  )*dyt(m,n  )*hhv(m,n  )     &
               + vv(m,n-1,nz)*dxh(m,n-1)*dyt(m,n-1)*hhv(m,n-1)) /2.0d0/dx(m,n)/dy(m,n)/hhq(m,n)
           BottomFriction(m,n) = Cb_nl * dsqrt(ut*ut + vt*vt+Ebottom)
          endif
         enddo
        enddo
!$omp end parallel do
!
! 3. E r r o r
! ------------
     else
       stop 'bottomfriction'
   endif !TYPE_FRIC

   call syncborder_real8(BottomFriction, 1)

   if(periodicity_x/=0) then
    call cyclize8_x(BottomFriction,nx,ny,1,mmm,mm)
   endif

   if(periodicity_y/=0) then
    call cyclize8_y(BottomFriction,nx,ny,1,nnn,nn)
   endif

!$omp parallel do
   do n=ny_start,ny_end
    do m=nx_start,nx_end
     if(lcu(m,n)>0.5) then
      bot_stress_x(m,n)= - (BottomFriction(m,n)+BottomFriction(m+1,n))/2.0d0 * uup(m,n,nz)
     endif
     if(lcv(m,n)>0.5) then
      bot_stress_y(m,n)= - (BottomFriction(m,n)+BottomFriction(m,n+1))/2.0d0 * vvp(m,n,nz)
     endif
    enddo
   enddo
!$omp end parallel do

!Boundary condition for turbulent kinetic energy
if(iabs(ksw_vert)>1) then
  !$omp parallel do private(m,n)
   do n=ny_start,ny_end
     do m=nx_start,nx_end
      if(lu(m,n)>0.5) then
        q2_bot(m,n)= dsqrt( ( bot_stress_x(m  ,n  )**2*dxt(m  ,n  )*dyh(m  ,n  )     &
                             +bot_stress_x(m-1,n  )**2*dxt(m-1,n  )*dyh(m-1,n  )     &
                             +bot_stress_y(m  ,n  )**2*dxh(m  ,n  )*dyt(m  ,n  )     &
                             +bot_stress_y(m  ,n-1)**2*dxh(m  ,n-1)*dyt(m  ,n-1) )/2.0d0/dx(m,n)/dy(m,n) )      &
                              *B1_t**(2.0d0/3.0d0)
      endif
     enddo
   enddo
  !$omp end parallel do
endif

endsubroutine sea_bottom_fluxes
