!===========================================================================================
subroutine uv_trans( u, v, wflux, vort,    &
                   hq, hu, hv, hh,         &
                   RHSx, RHSy, nlev    )
use main_basin_pars
use mpi_parallel_tools
use basin_grid
implicit none

 integer nlev

 real(8) u(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nlev),        & !Transporting zonal velocity
         v(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nlev)           !Transporting meridional velocity

 real(8) RHSx(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nlev ),      & !Zonal source function
         RHSy(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nlev ),      & !meridional source function
         wflux(bnd_x1:bnd_x2,bnd_y1:bnd_y2)              !source term in continuity equation

 real(8) hq(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &
         hu(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &
         hv(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &
         hh(bnd_x1:bnd_x2,bnd_y1:bnd_y2)

real(8) vort(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nlev)

real(8) fx_p,fx_m,fy_p,fy_m   !fluxes through cell edges

integer m,n,k

!$omp parallel do private(m,n,k)
 do n=ny_start, ny_end
   do m=nx_start, nx_end
    if(luu(m,n)>0.5) then
     do k=1,nlev
      vort(m,n,k)= (v(m+1,n,k)*dyt(m+1,n)-v(m,n,k)*dyt(m,n))     &
                  -(u(m,n+1,k)*dxt(m,n+1)-u(m,n,k)*dxt(m,n))     &
                  -((v(m+1,n,k)-v(m,n,k))*dyb(m,n)-(u(m,n+1,k)-u(m,n,k))*dxb(m,n))
     enddo
    endif
   enddo
 enddo
!$omp end parallel do

      call syncborder_real8(vort)

      if(periodicity_x/=0) then
       call cyclize8_x(vort,nx,ny,nlev,mmm,mm)
	end if

      if(periodicity_y/=0) then
       call cyclize8_y(vort,nx,ny,nlev,nnn,nn)
	end if

!$omp parallel do private(m,n,k,fx_p,fx_m,fy_p,fy_m)
  do n=ny_start,ny_end
    do m=nx_start,nx_end

!zonal velocity
      if(lcu(m,n)>0.5) then

        do k=1,nlev

         fx_p=(u(m  ,n  ,k)*dyh(m,n)*hu(m,n) + u(m+1,n  ,k)*dyh(m+1,n)*hu(m+1,n))/2.0d0   &
             *(u(m  ,n  ,k) + u(m+1,n  ,k))/2.0d0

         fx_m=(u(m  ,n  ,k)*dyh(m,n)*hu(m,n) + u(m-1,n  ,k)*dyh(m-1,n)*hu(m-1,n))/2.0d0   &
             *(u(m  ,n  ,k) + u(m-1,n  ,k))/2.0d0

         fy_p=(v(m  ,n  ,k)*dxh(m,n  )*hv(m,n  ) + v(m+1,n  ,k)*dxh(m+1,n  )*hv(m+1,n  ))/2.0d0   &
             *(u(m  ,n+1,k) + u(m  ,n  ,k))/2.0d0*dble(luu(m,n  ))

         fy_m=(v(m  ,n-1,k)*dxh(m,n-1)*hv(m,n-1) + v(m+1,n-1,k)*dxh(m+1,n-1)*hv(m+1,n-1))/2.0d0   &
             *(u(m  ,n-1,k) + u(m  ,n  ,k))/2.0d0*dble(luu(m,n-1))

         RHSx(m,n,k)= - (fx_p - fx_m + fy_p - fy_m)           &
                      + ( wflux(m  ,n)*dx(m  ,n)*dy(m  ,n)    &
                         +wflux(m+1,n)*dx(m+1,n)*dy(m+1,n) )*u(m,n,k)      &
                            /2.0d0/RefDen*dfloat(full_free_surface)    &
            + ( vort(m,n  ,k)*hh(m,n  )*(v(m+1,n  ,k)+v(m,n  ,k))              &
            +   vort(m,n-1,k)*hh(m,n-1)*(v(m+1,n-1,k)+v(m,n-1,k))  )/4.0d0

        end do

      end if

!meridional velocity
      if(lcv(m,n)>0.5) then

        do k=1,nlev

         fy_p=(v(m  ,n  ,k)*dxh(m,n)*hv(m,n) + v(m  ,n+1,k)*dxh(m,n+1)*hv(m,n+1))/2.0d0    &
             *(v(m  ,n  ,k) + v(m  ,n+1,k))/2.0d0

         fy_m=(v(m  ,n  ,k)*dxh(m,n)*hv(m,n) + v(m  ,n-1,k)*dxh(m,n-1)*hv(m,n-1))/2.0d0    &
             *(v(m  ,n  ,k) + v(m  ,n-1,k))/2.0d0

	   fx_p=(u(m  ,n  ,k)*dyh(m  ,n)*hu(m  ,n) + u(m  ,n+1,k)*dyh(m  ,n+1)*hu(m  ,n+1))/2.0d0    &
             *(v(m+1,n  ,k) + v(m  ,n  ,k))/2.0d0

	   fx_m=(u(m-1,n  ,k)*dyh(m-1,n)*hu(m-1,n) + u(m-1,n+1,k)*dyh(m-1,n+1)*hu(m-1,n+1))/2.0d0    &
             *(v(m-1,n  ,k) + v(m  ,n  ,k))/2.0d0

         RHSy(m,n,k)= - (fx_p - fx_m + fy_p - fy_m)          &
                      + (wflux(m,n  )*dx(m,n  )*dy(m,n  )    &
                        +wflux(m,n+1)*dx(m,n+1)*dy(m,n+1))*v(m,n,k)     &
                        /2.0d0/RefDen *dfloat(full_free_surface)        &
             - ( vort(m  ,n,k)*hh(m  ,n)*(u(m  ,n+1,k)+u(m  ,n,k))               &
             +   vort(m-1,n,k)*hh(m-1,n)*(u(m-1,n+1,k)+u(m-1,n,k))  )/4.0d0
        end do

      end if

    end do
  end do
!$omp end parallel do

endsubroutine uv_trans

!===========================================================================================
subroutine uv_diff2( mu, str_t, str_s,    &
                     hq, hu, hv, hh,      &
                     RHSx, RHSy, nlev     )
use main_basin_pars
use mpi_parallel_tools
use basin_grid
implicit none

 integer nlev
 real(8) muh_p, muh_m

 real(8) mu(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nlev ),      & !lateral viscosity coefficient
       RHSx(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nlev ),      & !Zonal source function
       RHSy(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nlev ),      & !meridional source function
      str_t(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nlev ),      & !Tension stress
      str_s(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nlev )         !Shearing stress

 real(8) hq(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &
         hu(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &
         hv(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &
         hh(bnd_x1:bnd_x2,bnd_y1:bnd_y2)

integer m,n,k

!$omp parallel do private(muh_p, muh_m)
  do n=ny_start,ny_end
    do m=nx_start,nx_end

!zonal velocity
      if(lcu(m,n)>0.5) then

        do k=1,nlev

         muh_p=(mu(m,n,k)+mu(m+1,n,k)+mu(m,n+1,k)+mu(m+1,n+1,k))/4.0d0
         muh_m=(mu(m,n,k)+mu(m+1,n,k)+mu(m,n-1,k)+mu(m+1,n-1,k))/4.0d0

         RHSx(m,n,k)=( dy(m+1,n)**2*mu(m+1,n,k)*hq(m+1,n)*str_t(m+1,n,k)             &
                      -dy(m  ,n)**2*mu(m  ,n,k)*hq(m  ,n)*str_t(m  ,n,k) )/dyh(m,n)  &
                   + (dxb(m,n  )**2*muh_p*hh(m,n  )*str_s(m,n  ,k)                   &
                     -dxb(m,n-1)**2*muh_m*hh(m,n-1)*str_s(m,n-1,k) )/dxt(m,n)
        end do

      end if

!meridional velocity
      if(lcv(m,n)>0.5) then

        do k=1,nlev

         muh_p=(mu(m,n,k)+mu(m+1,n,k)+mu(m,n+1,k)+mu(m+1,n+1,k))/4.0d0
         muh_m=(mu(m,n,k)+mu(m-1,n,k)+mu(m,n+1,k)+mu(m-1,n+1,k))/4.0d0

         RHSy(m,n,k)=-( dx(m,n+1)**2*mu(m,n+1,k)*hq(m,n+1)*str_t(m,n+1,k)              &
                       -dx(m,n  )**2*mu(m,n  ,k)*hq(m,n  )*str_t(m,n  ,k) ) /dxh(m,n)  &
                    + (dyb(m  ,n)**2*muh_p*hh(m  ,n)*str_s(m  ,n,k)                    &
                      -dyb(m-1,n)**2*muh_m*hh(m-1,n)*str_s(m-1,n,k) ) /dyt(m,n)
        end do

      end if

    end do
  end do
!$omp end parallel do

endsubroutine uv_diff2

!================================================================================
subroutine uv_diff4( mu, str_t, str_s,      &
               fx, fy, hq, hu, hv, hh,      &
               RHSx, RHSy, nlev     )
use main_basin_pars
use mpi_parallel_tools
use basin_grid
implicit none
 integer nlev
 real(8) muh_p, muh_m

 real(8) mu(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nlev ),      & !lateral viscosity coefficient

       RHSx(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nlev ),      & !Zonal source function
       RHSy(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nlev ),      & !meridional source function
         fx(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nlev ),      & !Temporary array
         fy(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nlev ),      & !Temporary array
      str_t(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nlev ),      & !Tension stress
      str_s(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nlev )         !Shearing stress

 real(8) hq(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &
         hu(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &
         hv(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &
         hh(bnd_x1:bnd_x2,bnd_y1:bnd_y2)

integer m,n,k

fx=0.0d0
fy=0.0d0

!$omp parallel do private(muh_p, muh_m)
  do n=ny_start,ny_end
    do m=nx_start,nx_end

!zonal velocity
      if(lcu(m,n)>0.5) then

        do k=1,nlev

         muh_p=(mu(m,n,k)+mu(m+1,n,k)+mu(m,n+1,k)+mu(m+1,n+1,k))/4.0d0
         muh_m=(mu(m,n,k)+mu(m+1,n,k)+mu(m,n-1,k)+mu(m+1,n-1,k))/4.0d0

           fx(m,n,k)=( dy(m+1,n)**2*mu(m+1,n,k)*hq(m+1,n)*str_t(m+1,n,k)             &
                      -dy(m  ,n)**2*mu(m  ,n,k)*hq(m  ,n)*str_t(m  ,n,k) )/dyh(m,n)  &
                   + (dxb(m,n  )**2*muh_p*hh(m,n  )*str_s(m,n  ,k)                   &
                     -dxb(m,n-1)**2*muh_m*hh(m,n-1)*str_s(m,n-1,k) )/dxt(m,n)
           fx(m,n,k)=-fx(m,n,k)/hu(m,n)/dxt(m,n)/dyh(m,n)
        end do

      end if

!meridional velocity
      if(lcv(m,n)>0.5) then

        do k=1,nlev

         muh_p=(mu(m,n,k)+mu(m+1,n,k)+mu(m,n+1,k)+mu(m+1,n+1,k))/4.0d0
         muh_m=(mu(m,n,k)+mu(m-1,n,k)+mu(m,n+1,k)+mu(m-1,n+1,k))/4.0d0

           fy(m,n,k)=-( dx(m,n+1)**2*mu(m,n+1,k)*hq(m,n+1)*str_t(m,n+1,k)              &
                       -dx(m,n  )**2*mu(m,n  ,k)*hq(m,n  )*str_t(m,n  ,k) ) /dxh(m,n)  &
                    + (dyb(m  ,n)**2*muh_p*hh(m  ,n)*str_s(m  ,n,k)                    &
                      -dyb(m-1,n)**2*muh_m*hh(m-1,n)*str_s(m-1,n,k) ) /dyt(m,n)
           fy(m,n,k)=-fy(m,n,k)/hv(m,n)/dxh(m,n)/dyt(m,n)

        end do

      end if

    end do
  end do
!$omp end parallel do

      if(periodicity_x/=0) then
       call cyclize8_x(fx,nx,ny,nlev,mmm,mm)
       call cyclize8_x(fy,nx,ny,nlev,mmm,mm)
	end if

      if(periodicity_y/=0) then
       call cyclize8_y(fx,nx,ny,nlev,nnn,nn)
       call cyclize8_y(fy,nx,ny,nlev,nnn,nn)
	end if

  call stress_components(fx,fy,str_t,str_s,nlev)

!$omp parallel do private(muh_p, muh_m)
  do n=ny_start,ny_end
    do m=nx_start,nx_end

!zonal velocity
      if(lcu(m,n)>0.5) then

        do k=1,nlev

         muh_p=(mu(m,n,k)+mu(m+1,n,k)+mu(m,n+1,k)+mu(m+1,n+1,k))/4.0d0
         muh_m=(mu(m,n,k)+mu(m+1,n,k)+mu(m,n-1,k)+mu(m+1,n-1,k))/4.0d0

         RHSx(m,n,k)=RHSx(m,n,k) + ( dy(m+1,n)**2*mu(m+1,n,k)*hq(m+1,n)*str_t(m+1,n,k)             &
                                    -dy(m  ,n)**2*mu(m  ,n,k)*hq(m  ,n)*str_t(m  ,n,k) )/dyh(m,n)  &
                                 + (dxb(m,n  )**2*muh_p*hh(m,n  )*str_s(m,n  ,k)                   &
                                   -dxb(m,n-1)**2*muh_m*hh(m,n-1)*str_s(m,n-1,k) )/dxt(m,n)
        end do

      end if

!meridional velocity
      if(lcv(m,n)>0.5) then

        do k=1,nlev

         muh_p=(mu(m,n,k)+mu(m+1,n,k)+mu(m,n+1,k)+mu(m+1,n+1,k))/4.0d0
         muh_m=(mu(m,n,k)+mu(m-1,n,k)+mu(m,n+1,k)+mu(m-1,n+1,k))/4.0d0

         RHSy(m,n,k)=RHSy(m,n,k) - ( dx(m,n+1)**2*mu(m,n+1,k)*hq(m,n+1)*str_t(m,n+1,k)              &
                                    -dx(m,n  )**2*mu(m,n  ,k)*hq(m,n  )*str_t(m,n  ,k) ) /dxh(m,n)  &
                                 + (dyb(m  ,n)**2*muh_p*hh(m  ,n)*str_s(m  ,n,k)                    &
                                   -dyb(m-1,n)**2*muh_m*hh(m-1,n)*str_s(m-1,n,k) ) /dyt(m,n)
        end do

      end if

    end do
  end do
!$omp end parallel do

endsubroutine uv_diff4

!================================================================================
subroutine pressure_gradients(den,        &
                              RHSx,       &
                              RHSy)
use main_basin_pars
use mpi_parallel_tools
use basin_grid
implicit none

 real(8) den(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),     &       !water density
        RHSx(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),     &       !RHS in zonal direction
        RHSy(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz)              !RHS in meridional direction

 real(8) p_hint, p_lat

 integer m,n,k

 !$omp parallel do private(m,n,k,p_hint, p_lat)
 do n=ny_start,ny_end
   do m=nx_start,nx_end

     if(lcu(m,n)>0.5) then
          p_hint =  hhq(m+1,n)*z(1)*den(m+1,n,1) - hhq(m,n)*z(1)*den(m  ,n,1)
          p_lat  =  hhq(m+1,n)*z(1)*den(m  ,n,1) - hhq(m,n)*z(1)*den(m+1,n,1)
          RHSx(m,n,1)=RHSx(m,n,1)-(p_hint-p_lat)*hhu(m,n)*dyh(m,n)*FreeFallAcc/RefDen/2.0d0

	   do k=2,nz

          p_hint = p_hint + ( hhq(m+1,n)*z(k)*den(m+1,n,k-1) - hhq(m+1,n)*z(k-1)*den(m+1,n,k) )     &
                           -( hhq(m  ,n)*z(k)*den(m  ,n,k-1) - hhq(m  ,n)*z(k-1)*den(m  ,n,k) )
          p_lat  = hhq(m+1,n)*z(k)*den(m,n,k)-hhq(m,n)*z(k)*den(m+1,n,k)
          RHSx(m,n,k)=RHSx(m,n,k)-(p_hint-p_lat)*hhu(m,n)*dyh(m,n)*FreeFallAcc/RefDen/2.0d0
         enddo
     endif

     if(lcv(m,n)>0.5) then
          p_hint = hhq(m,n+1)*z(1)*den(m,n+1,1) - hhq(m,n)*z(1)*den(m,n  ,1)
          p_lat  = hhq(m,n+1)*z(1)*den(m,n  ,1) - hhq(m,n)*z(1)*den(m,n+1,1)
          RHSy(m,n,1)=RHSy(m,n,1)-(p_hint-p_lat)*hhv(m,n)*dxh(m,n)*FreeFallAcc/RefDen/2.0d0

	   do k=2,nz

          p_hint = p_hint + ( hhq(m,n+1)*z(k)*den(m,n+1,k-1) - hhq(m,n+1)*z(k-1)*den(m,n+1,k) )     &
                           -( hhq(m,n  )*z(k)*den(m,n  ,k-1) - hhq(m,n  )*z(k-1)*den(m,n  ,k) )
          p_lat =  hhq(m,n+1)*z(k)*den(m,n  ,k)-hhq(m,n)*z(k)*den(m,n+1,k)
          RHSy(m,n,k)=RHSy(m,n,k)-(p_hint-p_lat)*hhv(m,n)*dxh(m,n)*FreeFallAcc/RefDen/2.0d0
         enddo

     endif

   enddo
 enddo
 !$omp end parallel do

endsubroutine pressure_gradients
!==========================================================================================================
!shallow water equation sloving
subroutine barotropic_dynamics(tau,     &
                             nstep,     &
                              ksw4,     &
                             ubrtr_e,   &
                             ubrtrp_e,  &
                             vbrtr_e,   &
                             vbrtrp_e,  &
                             ssh_e,     &
                             sshp_e,    &
                             ubrtr_i,   &
                             vbrtr_i,   &
                            ssh4gradx,  &
                            ssh4grady,  &
                                 RHSx,  &
                                 RHSy,  &
                                wflux,  &
                                   mu,  &
                                  mu4,  &
                                 vort,  &
                              str_t2d,  &
                              str_s2d,  &
                                   fx,  &
                                   fy,  &
                                 rdis,  &
                             RHSx_adv,  &
                             RHSy_adv,  &
                             RHSx_dif,  &
                             RHSy_dif)
use main_basin_pars
use mpi_parallel_tools
use basin_grid
implicit none

real(8) tau, tau_inner
integer step, nstep, ksw4, m, n

real(8) ubrtr_e(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
       ubrtrp_e(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
        vbrtr_e(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
       vbrtrp_e(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
          ssh_e(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
         sshp_e(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
        ubrtr_i(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
        vbrtr_i(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
      ssh4gradx(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
      ssh4grady(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
          wflux(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
           RHSx(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
           RHSy(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
             mu(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
            mu4(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
           vort(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
        str_t2d(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
        str_s2d(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
             fx(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
             fy(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
           rdis(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
       RHSx_adv(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
       RHSy_adv(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
       RHSx_dif(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
       RHSy_dif(bnd_x1:bnd_x2,bnd_y1:bnd_y2)

real(8), allocatable::   u(:,:),   &
                         v(:,:),   &
                       ssh(:,:),   &
                        up(:,:),   &
                        vp(:,:),   &
                      sshp(:,:),   &
                        un(:,:),   &
                        vn(:,:),   &
                      sshn(:,:)


real(8) bp, bp0, grx, gry, slx, sly, slxn, slyn

allocate(  u(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &
           v(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &
         ssh(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &
          up(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &
          vp(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &
        sshp(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &
          un(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &
          vn(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &
        sshn(bnd_x1:bnd_x2,bnd_y1:bnd_y2)    )


  tau_inner=tau/dfloat(nstep)

   un=0.0d0
   vn=0.0d0
 sshn=0.0d0

   up=ubrtrp_e
   vp=vbrtrp_e
 sshp=sshp_e

    u=ubrtr_e
    v=vbrtr_e
  ssh=ssh_e

! ubrtr_i   =0.0d0
! vbrtr_i   =0.0d0
! ssh4gradx =0.0d0
! ssh4grady =0.0d0

do step=1,2*nstep

!computing ssh
!$omp parallel do
 do n=ny_start,ny_end
  do m=nx_start,nx_end

    if(lu(m,n)>0.5) then
     sshn(m,n)=sshp(m,n) + 2.0d0*tau_inner*( wflux(m,n)/RefDen*dfloat(full_free_surface)   &
                - ( u(m,n)*hhu_e(m,n)*dyh(m,n)-u(m-1,n)*hhu_e(m-1,n)*dyh(m-1,n)            &
                  + v(m,n)*hhv_e(m,n)*dxh(m,n)-v(m,n-1)*hhv_e(m,n-1)*dxh(m,n-1) )/(dx(m,n)*dy(m,n))  )
    endif

  enddo
 enddo
!$omp end parallel do

 call syncborder_real8(sshn)

 if(periodicity_x/=0) then
   call cyclize8_x(sshn,nx,ny,1,mmm,mm)
 endif

 if(periodicity_y/=0) then
   call cyclize8_y(sshn,nx,ny,1,nnn,nn)
 endif

   if(full_free_surface>0) then
    call hh_update(hhqn_e, hhun_e, hhvn_e, hhhn_e, sshn, hhq_rest)
   endif

! computing advective and lateral-viscous terms for 2d-velocity
! call stress_components(up,vp,str_t2d,str_s2d,1)
! computing advective and lateral-viscous terms for 2d-velocity
 call uv_trans( u, v, wflux, vort,            &
              hhq_e, hhu_e, hhv_e, hhh_e,     &
              RHSx_adv, RHSy_adv, 1  )
! call uv_diff2( mu, str_t2d, str_s2d,          &
!               hhq_e, hhu_e, hhv_e, hhh_e,     &
!               RHSx_dif, RHSy_dif, 1  )

! if(ksw4>0) then
!   call uv_diff4( mu4, str_t2d, str_s2d,  &
!                  fx, fy, hhq_e, hhu_e, hhv_e, hhh_e,    &
!                  RHSx_dif, RHSy_dif, 1 )
! endif

!$omp parallel do private(bp,bp0,grx,gry, slx, sly, slxn, slyn)
 do n=ny_start,ny_end
  do m=nx_start,nx_end

!zonal flux
    if(lcu(m,n)>0.5) then

      bp =hhun_e(m,n)*dxt(m,n)*dyh(m,n)/2.0d0/tau_inner
      bp0=hhup_e(m,n)*dxt(m,n)*dyh(m,n)/2.0d0/tau_inner

     slx = - FreeFallAcc * ( ssh(m+1,n) - ssh(m,n))*dyh(m,n)* hhu_e(m,n)
     slxn= - FreeFallAcc * (sshn(m+1,n) -sshn(m,n))*dyh(m,n)*hhun_e(m,n)
     grx= RHSx(m,n) +slx  +RHSx_dif(m,n)+RHSx_adv(m,n)      &
          - (rdis(m,n)+rdis(m+1,n))/2.0d0 *up(m,n)*dxt(m,n)*dyh(m,n)*hhu_e(m,n)        &
          + ( rlh_s(m,n  )*hhh_e(m,n  )*dxb(m,n  )*dyb(m,n  )*(v(m+1,n  )+v(m,n  ))             &
            + rlh_s(m,n-1)*hhh_e(m,n-1)*dxb(m,n-1)*dyb(m,n-1)*(v(m+1,n-1)+v(m,n-1))  )/4.0d0

     un(m,n)=(up(m,n)*bp0 + grx )/bp
     !ubrtr_i(m,n)  = ubrtr_i(m,n)+ (u(m,n)*hhu_e(m,n)+un(m,n)*hhun_e(m,n))/dfloat(4*nstep)
     !ssh4gradx(m,n)= ssh4gradx(m,n)+(slx+slxn)/dfloat(4*nstep)
    endif

!meridional flux
    if(lcv(m,n)>0.5) then

      bp =hhvn_e(m,n)*dyt(m,n)*dxh(m,n)/2.0d0/tau_inner
      bp0=hhvp_e(m,n)*dyt(m,n)*dxh(m,n)/2.0d0/tau_inner

     sly = - FreeFallAcc * ( ssh(m,n+1)- ssh(m,n))*dxh(m,n)* hhv_e(m,n)
     slyn= - FreeFallAcc * (sshn(m,n+1)-sshn(m,n))*dxh(m,n)*hhvn_e(m,n)
     gry= RHSy(m,n) +sly  +RHSy_dif(m,n)+RHSy_adv(m,n)      &
          - (rdis(m,n)+rdis(m,n+1))/2.0d0 *vp(m,n)*dxh(m,n)*dyt(m,n)*hhv_e(m,n)        &
          - ( rlh_s(m  ,n)*hhh_e(m  ,n)*dxb(m  ,n)*dyb(m  ,n)*(u(m  ,n+1)+u(m  ,n))             &
            + rlh_s(m-1,n)*hhh_e(m-1,n)*dxb(m-1,n)*dyb(m-1,n)*(u(m-1,n+1)+u(m-1,n))  )/4.0d0

     vn(m,n)=(vp(m,n)*bp0 + gry )/bp
     !vbrtr_i(m,n)  = vbrtr_i(m,n)+(v(m,n)*hhv_e(m,n)+vn(m,n)*hhvn_e(m,n))/dfloat(4*nstep)
     !ssh4grady(m,n)= ssh4grady(m,n)+(sly+slyn)/dfloat(4*nstep)
    endif

  enddo
 enddo
!$omp end parallel do

 call syncborder_real8(un)
 call syncborder_real8(vn)

 if(periodicity_x/=0) then
   call cyclize8_x(  un,nx,ny,1,mmm,mm)
   call cyclize8_x(  vn,nx,ny,1,mmm,mm)
 endif

 if(periodicity_y/=0) then
   call cyclize8_y(  un,nx,ny,1,nnn,nn)
   call cyclize8_y(  vn,nx,ny,1,nnn,nn)
 endif

 !shifting time indices
  !$omp parallel do private(m,n)
      do n=ny_start-1,ny_end+1
       do m=nx_start-1,nx_end+1

        if(lu(m,n)>0.5) then
          sshp(m,n) =  ssh(m,n)+time_smooth*(sshn(m,n)-2.0d0*ssh(m,n)+sshp(m,n))/2.0d0/dfloat(nstep)
           ssh(m,n) =sshn(m,n)
        endif

        if(lcu(m,n)>0.5) then
            up(m,n) =  u(m,n)+time_smooth*(un(m,n)-2.0d0*u(m,n)+up(m,n))/2.0d0/dfloat(nstep)
             u(m,n) = un(m,n)
        endif

        if(lcv(m,n)>0.5) then
            vp(m,n) =  v(m,n)+time_smooth*(vn(m,n)-2.0d0*v(m,n)+vp(m,n))/2.0d0/dfloat(nstep)
             v(m,n) = vn(m,n)
        endif

       enddo
      enddo
!$omp end parallel do

!      sshp = ssh
!      ssh = sshn
!      up = u
!      u = un
!      vp = v
!      v = vn


  if(full_free_surface>0) then
    call hh_shift(hhq_e, hhqp_e, hhqn_e,   &
                  hhu_e, hhup_e, hhun_e,   &
                  hhv_e, hhvp_e, hhvn_e,   &
                  hhh_e, hhhp_e, hhhn_e, nstep )
  endif

 if(step==nstep) then

   ubrtrp_e=up
   vbrtrp_e=vp
     sshp_e=sshp

   ubrtr_e=u
   vbrtr_e=v
     ssh_e=ssh

   ssh4gradx = ssh
   ssh4grady = ssh
   ubrtr_i = u
   vbrtr_i = v

 endif

enddo

! if(periodicity_x/=0) then
!   call cyclize8_x(ubrtr_i,nx,ny,1,mmm,mm)
!   call cyclize8_x(vbrtr_i,nx,ny,1,mmm,mm)
!!   call cyclize8_x(ssh4gradx,nx,ny,1,mmm,mm)
!!   call cyclize8_x(ssh4grady,nx,ny,1,mmm,mm)
! endif

! if(periodicity_y/=0) then
!   call cyclize8_y(ubrtr_i,nx,ny,1,nnn,nn)
!   call cyclize8_y(vbrtr_i,nx,ny,1,nnn,nn)
!!   call cyclize8_y(ssh4gradx,nx,ny,1,nnn,nn)
!!   call cyclize8_y(ssh4grady,nx,ny,1,nnn,nn)
! endif

 if(full_free_surface>0) then
!initialize depth for external mode
  call hh_init(hhq_e, hhqp_e, hhqn_e,    &
               hhu_e, hhup_e, hhun_e,    &
               hhv_e, hhvp_e, hhvn_e,    &
               hhh_e, hhhp_e, hhhn_e,    &
               ssh_e, sshp_e, hhq_rest)
  endif

deallocate(sshn,vn,un,sshp,vp,up,ssh,v,u)

endsubroutine barotropic_dynamics

!======================================================================
subroutine ssh_internal(tau,    &
                      ssh_i,    &
                     sshp_i,    &
                          u,    &
                          v,    &
                      wflux   )
 use main_basin_pars
 use mpi_parallel_tools
 use basin_grid
 implicit none

 real(8) tau
 real(8)  ssh_i(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
         sshp_i(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
              u(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
              v(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
          wflux(bnd_x1:bnd_x2,bnd_y1:bnd_y2)

 real(8), allocatable:: ssh_new(:,:)

 integer m,n

 allocate( ssh_new(bnd_x1:bnd_x2,bnd_y1:bnd_y2)   )
 ssh_new=0.0d0

 !$omp parallel do
  do n=ny_start, ny_end
    do m=nx_start, nx_end
     if(lu(m,n)>0.5) then

      ssh_new(m,n)=sshp_i(m,n)      &
           + 2.0d0*tau*( wflux(m,n)/RefDen*dfloat(full_free_surface)    &
                       - (u(m,n)*dyh(m,n)-u(m-1,n)*dyh(m-1,n)           &
                        + v(m,n)*dxh(m,n)-v(m,n-1)*dxh(m,n-1))/(dx(m,n)*dy(m,n))  )
     endif
    enddo
  enddo
 !$omp end parallel do

  if(periodicity_x/=0) then
    call cyclize8_x(ssh_new,nx,ny,1,mmm,mm)
  endif

  if(periodicity_y/=0) then
    call cyclize8_y(ssh_new,nx,ny,1,nnn,nn)
  endif


  !computing new depths
  if(full_free_surface>0) then
    call hh_update(hhqn, hhun, hhvn, hhhn, ssh_new, hhq_rest)
  endif

   !Updating ssh function
  !$omp parallel do private(m,n)
       do n=ny_start-1,ny_end+1
        do m=nx_start-1,nx_end+1

         if(lu(m,n)>0.5) then
           sshp_i(m,n)= ssh_i(m,n) + time_smooth*(ssh_new(m,n)-2.0d0*ssh_i(m,n)+sshp_i(m,n))/2.0d0
            ssh_i(m,n)= ssh_new(m,n)
         endif

        end do
 	end do
!$omp end parallel do
 deallocate(ssh_new)
endsubroutine ssh_internal

!======================================================================
subroutine ssh_gradients(ssh,gradx,grady)
use main_basin_pars
use mpi_parallel_tools
use basin_grid
implicit none
      real(8) ssh(bnd_x1:bnd_x2,bnd_y1:bnd_y2),       &
	      gradx(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &
            grady(bnd_x1:bnd_x2,bnd_y1:bnd_y2)
	integer m,n

!$omp parallel do
      do n=ny_start,ny_end
	  do m=nx_start,nx_end

           if(lcu(m,n)>0.5) then
            gradx(m,n)=-FreeFallAcc*(ssh(m+1,n)-ssh(m,n))/dxt(m,n)
           endif

           if(lcv(m,n)>0.5) then
            grady(m,n)=-FreeFallAcc*(ssh(m,n+1)-ssh(m,n))/dyt(m,n)
           endif

        end do
	end do
!$omp end parallel do

!      if(periodicity_x/=0) then
!          call cyclize8_x(gradx,nx,ny, 1,mmm,mm)
!      end if

!      if(periodicity_y/=0) then
!          call cyclize8_y(grady,nx,ny, 1,nnn,nn)
!      end if

endsubroutine ssh_gradients

!=====================================================================
subroutine baroclinic_dynamics(tau,       &
                               uu,        &
                               uup,       &
                               vv,        &
                               vvp,       &
                               ww,        &
                               nu,        &
                               RHSx,      &
                               RHSy,      &
                               RHSx_tran, &
                               RHSy_tran, &
                                 pgrx,    &
                                 pgry,    &
                                 rdis,    &
                               tx_surf,   &
                               ty_surf,   &
                               tx_bot,    &
                               ty_bot)
use main_basin_pars
use mpi_parallel_tools
use basin_grid
implicit none

real(8) tau, bp, bp0, dp, dm, fz_p, fz_m

real(8) uu(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),       &
       uup(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),       &
        vv(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),       &
       vvp(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),       &
        ww(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1),     &
        nu(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1),     &
      RHSx(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),       &
      RHSy(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),       &
      RHSx_tran(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),  &
      RHSy_tran(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),  &
      pgrx(bnd_x1:bnd_x2,bnd_y1:bnd_y2),          &
      pgry(bnd_x1:bnd_x2,bnd_y1:bnd_y2),          &
      rdis(bnd_x1:bnd_x2,bnd_y1:bnd_y2),          &
   tx_surf(bnd_x1:bnd_x2,bnd_y1:bnd_y2),          &
   ty_surf(bnd_x1:bnd_x2,bnd_y1:bnd_y2),          &
    tx_bot(bnd_x1:bnd_x2,bnd_y1:bnd_y2),          &
    ty_bot(bnd_x1:bnd_x2,bnd_y1:bnd_y2)

integer m,n,k
real(8) a(nz),b(nz),c(nz),eta(nz),rksi(nz)
real(8), allocatable::  un(:,:,:), vn(:,:,:)

allocate(un(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),    &
         vn(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz)    )

         un=0.0d0
         vn=0.0d0

!$omp parallel do private(m,n,k,bp,bp0,dp,dm,a,b,c,eta,rksi,fz_p,fz_m)
      do n=ny_start,ny_end
       do m=nx_start,nx_end

        if (lcu(m,n)>0.5) then
         bp = hhun(m,n)*dxt(m,n)*dyh(m,n) / tau/2.0d0
         bp0= hhup(m,n)*dxt(m,n)*dyh(m,n) / tau/2.0d0

! surface point
         k=1

          dm=0.0d0
          dp =(nu(m,n,k+1)+nu(m+1,n,k+1))/2.0d0/hhun(m,n)/hzt(k+1)*dxt(m,n)*dyh(m,n)

          fz_m=0.0d0
          fz_p=(ww(m  ,n  ,k+1)*dx(m,n)*dy(m,n) + ww(m+1,n  ,k+1)*dx(m+1,n)*dy(m+1,n) ) /2.0d0    &
             * (uu(m  ,n  ,k  ) + uu(m  ,n  ,k+1) ) /2.0d0

          c(k) =  - dp/dz(k)
          a(k) =  0.0d0
          b(k) =  bp + dp/dz(k)
          eta(k) = bp0*uup(m,n,k)+RHSx(m,n,k)+RHSx_tran(m,n,k) + pgrx(m,n) - (fz_p-fz_m )/ dz(k)          &
                 + ( rlh_s(m,n  )*hhh(m,n  )*dxb(m,n  )*dyb(m,n  )*(vv(m+1,n  ,k)+vv(m,n  ,k))            &
                 +   rlh_s(m,n-1)*hhh(m,n-1)*dxb(m,n-1)*dyb(m,n-1)*(vv(m+1,n-1,k)+vv(m,n-1,k))  )/4.0d0   &
                 - (rdis(m,n)+rdis(m+1,n))/2.0d0 *uup(m,n,k)*dxt(m,n)*dyh(m,n)*hhu(m,n)                   &
                 + tx_surf(m,n)/dz(k)*dxt(m,n)*dyh(m,n)

! internal points.
         do k=2,nz-1

          dm = dp
          dp =(nu(m,n,k+1)+nu(m+1,n,k+1))/2.0d0/hhun(m,n)/hzt(k+1)*dxt(m,n)*dyh(m,n)

          fz_m=fz_p
          fz_p=(ww(m  ,n  ,k+1)*dx(m,n)*dy(m,n) + ww(m+1,n  ,k+1)*dx(m+1,n)*dy(m+1,n) ) /2.0d0    &
             * (uu(m  ,n  ,k  ) + uu(m  ,n  ,k+1) ) /2.0d0
          c(k) =  - dp/dz(k)
          a(k) =  - dm/dz(k)
          b(k) =  bp + (dp + dm)/dz(k)
          eta(k) = bp0*uup(m,n,k)+RHSx(m,n,k)+RHSx_tran(m,n,k) + pgrx(m,n) - (fz_p-fz_m )/ dz(k)          &
                 + ( rlh_s(m,n  )*hhh(m,n  )*dxb(m,n  )*dyb(m,n  )*(vv(m+1,n  ,k)+vv(m,n  ,k))            &
                 +   rlh_s(m,n-1)*hhh(m,n-1)*dxb(m,n-1)*dyb(m,n-1)*(vv(m+1,n-1,k)+vv(m,n-1,k))  )/4.0d0   &
                 - (rdis(m,n)+rdis(m+1,n))/2.0d0 *uup(m,n,k)*dxt(m,n)*dyh(m,n)*hhu(m,n)
         enddo

! bottom point
         k=nz
          dm = dp
          dp = 0.0d0

          fz_m=fz_p
          fz_p=0.0d0

          c(k) =  0.0d0
          a(k) =  - dm/dz(k)
          b(k) =  bp + dm/dz(k)
          eta(k) = bp0*uup(m,n,k)+RHSx(m,n,k)+RHSx_tran(m,n,k) + pgrx(m,n) - (fz_p-fz_m )/ dz(k)          &
                 + ( rlh_s(m,n  )*hhh(m,n  )*dxb(m,n  )*dyb(m,n  )*(vv(m+1,n  ,k)+vv(m,n  ,k))            &
                 +   rlh_s(m,n-1)*hhh(m,n-1)*dxb(m,n-1)*dyb(m,n-1)*(vv(m+1,n-1,k)+vv(m,n-1,k))  )/4.0d0   &
                 - (rdis(m,n)+rdis(m+1,n))/2.0d0 *uup(m,n,k)*dxt(m,n)*dyh(m,n)*hhu(m,n)                   &
                 + tx_bot(m,n)/dz(k)*dxt(m,n)*dyh(m,n)

         call factor8(nz,a,b,c,eta,rksi,1,nz)

         do k=1,nz
          un(m,n,k)=rksi(k)
         enddo

        endif

        if (lcv(m,n)>0.5) then
         bp =  hhvn(m,n)*dxh(m,n)*dyt(m,n) / tau/2.0d0
         bp0=  hhvp(m,n)*dxh(m,n)*dyt(m,n) / tau/2.0d0

! surface point
         k=1
          dm=0.0d0
          dp =(nu(m,n,k+1)+nu(m,n+1,k+1))/2.0d0/hhvn(m,n)/hzt(k+1)*dxh(m,n)*dyt(m,n)

          fz_m=0.0d0
          fz_p=(ww(m  ,n  ,k+1)*dx(m,n)*dy(m,n) + ww(m  ,n+1,k+1)*dx(m,n+1)*dy(m,n+1) ) /2.0d0    &
              *(vv(m  ,n  ,k  ) + vv(m  ,n  ,k+1) ) /2.0d0

          c(k) =  - dp/dz(k)
          a(k) =  0.0d0
          b(k) =  bp + dp/dz(k)
          eta(k) = bp0*vvp(m,n,k)+RHSy(m,n,k)+RHSy_tran(m,n,k) + pgry(m,n) - (fz_p-fz_m )/ dz(k)            &
                - ( rlh_s(m  ,n)*hhh(m  ,n)*dxb(m  ,n)*dyb(m  ,n)*(uu(m  ,n+1,k)+uu(m  ,n,k))               &
                +   rlh_s(m-1,n)*hhh(m-1,n)*dxb(m-1,n)*dyb(m-1,n)*(uu(m-1,n+1,k)+uu(m-1,n,k))  )/4.0d0      &
                - (rdis(m,n)+rdis(m,n+1))/2.0d0 *vvp(m,n,k)*dxh(m,n)*dyt(m,n)*hhv(m,n)                      &
                +  ty_surf(m,n)/dz(k)*dxh(m,n)*dyt(m,n)

! internal points.
         do k=2,nz-1

          dm = dp
          dp =(nu(m,n,k+1)+nu(m,n+1,k+1))/2.0d0/hhvn(m,n)/hzt(k+1)*dxh(m,n)*dyt(m,n)

          fz_m=fz_p
          fz_p=(ww(m  ,n  ,k+1)*dx(m,n)*dy(m,n) + ww(m  ,n+1,k+1)*dx(m,n+1)*dy(m,n+1) ) /2.0d0    &
              *(vv(m  ,n  ,k  ) + vv(m  ,n  ,k+1) ) /2.0d0

          c(k) =  - dp/dz(k)
          a(k) =  - dm/dz(k)
          b(k) =  bp + (dp + dm)/dz(k)
          eta(k) = bp0*vvp(m,n,k)+RHSy(m,n,k)+RHSy_tran(m,n,k) + pgry(m,n) - (fz_p-fz_m )/ dz(k)            &
                - ( rlh_s(m  ,n)*hhh(m  ,n)*dxb(m  ,n)*dyb(m  ,n)*(uu(m  ,n+1,k)+uu(m  ,n,k))               &
                +   rlh_s(m-1,n)*hhh(m-1,n)*dxb(m-1,n)*dyb(m-1,n)*(uu(m-1,n+1,k)+uu(m-1,n,k))  )/4.0d0      &
                - (rdis(m,n)+rdis(m,n+1))/2.0d0 *vvp(m,n,k)*dxh(m,n)*dyt(m,n)*hhv(m,n)
         enddo

! bottom point
         k=nz

          dm = dp
          dp = 0.0d0

          fz_m=fz_p
          fz_p=0.0d0

          c(k) =  0.0d0
          a(k) =  - dm/dz(k)
          b(k) =  bp + dm/dz(k)
          eta(k) = bp0*vvp(m,n,k)+RHSy(m,n,k)+RHSy_tran(m,n,k) + pgry(m,n) - (fz_p-fz_m )/ dz(k)            &
                - ( rlh_s(m  ,n)*hhh(m  ,n)*dxb(m  ,n)*dyb(m  ,n)*(uu(m  ,n+1,k)+uu(m  ,n,k))               &
                +   rlh_s(m-1,n)*hhh(m-1,n)*dxb(m-1,n)*dyb(m-1,n)*(uu(m-1,n+1,k)+uu(m-1,n,k))  )/4.0d0      &
                - (rdis(m,n)+rdis(m,n+1))/2.0d0 *vvp(m,n,k)*dxh(m,n)*dyt(m,n)*hhv(m,n)                      &
                + ty_bot(m,n)/dz(k)*dxh(m,n)*dyt(m,n)

         call factor8(nz,a,b,c,eta,rksi,1,nz)

         do k=1,nz
          vn(m,n,k)=rksi(k)
         enddo

        endif

       enddo
      enddo
!$omp end parallel do

 if(periodicity_x/=0)  then
  call cyclize8_x(un,nx,ny,nz,mmm,mm)
  call cyclize8_x(vn,nx,ny,nz,mmm,mm)
 end if

 if(periodicity_y/=0)  then
  call cyclize8_y(un,nx,ny,nz,nnn,nn)
  call cyclize8_y(vn,nx,ny,nz,nnn,nn)
 end if


!$omp parallel do private(m,n,k)
      do n=ny_start-1,ny_end+1
       do m=nx_start-1,nx_end+1

        if(lcu(m,n)>0.5) then
         do k=1,nz
          uup(m,n,k)=uu(m,n,k)+time_smooth*(un(m,n,k)-2.0d0*uu(m,n,k)+uup(m,n,k))/2.0d0
          uu(m,n,k)=un(m,n,k)
         enddo
        endif

        if(lcv(m,n)>0.5) then
         do k=1,nz
          vvp(m,n,k)=vv(m,n,k)+time_smooth*(vn(m,n,k)-2.0d0*vv(m,n,k)+vvp(m,n,k))/2.0d0
          vv(m,n,k)=vn(m,n,k)
         enddo
        endif

       enddo
      enddo
!$omp end parallel do

deallocate(vn,un)
endsubroutine baroclinic_dynamics

!======================================================================
subroutine vertical_velocity(u,v,w,hu,hv)
use main_basin_pars
use mpi_parallel_tools
use basin_grid

implicit none
!----------------------------------------------------------------------
! computing w from continuity equation on grid "c".
! input:
! u - zonal velocity on uc-grid
! v - meridional velocity on vc-grid
! output:
! w - vertical velocity on t-grid
real(8) u(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),        &
        v(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),        &
        w(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1),      &
        hu(bnd_x1:bnd_x2,bnd_y1:bnd_y2),          &
        hv(bnd_x1:bnd_x2,bnd_y1:bnd_y2)
integer m, n, k

!$omp parallel do private(m,n,k)
      do n=ny_start,ny_end
       do m=nx_start,nx_end
        if (lu(m,n)>0.5) then
         w(m,n,nz+1)=0.0d0
         do k=nz,1,-1
          w(m,n,k) = w(m,n,k+1) + dz(k)/(dx(m,n)*dy(m,n)) *     &
            (  u(m,n,k)*hu(m,n)*dyh(m,n) - u(m-1,n  ,k)*hu(m-1,n  )*dyh(m-1,n  )    +   &
               v(m,n,k)*hv(m,n)*dxh(m,n) - v(m  ,n-1,k)*hv(m  ,n-1)*dxh(m  ,n-1)  )
         enddo
        endif
       enddo
      enddo
!$omp end parallel do

  if(periodicity_x/=0) then
   call cyclize8_x(w,nx,ny,nz+1,mmm,mm)
  endif

  if(periodicity_y/=0) then
   call cyclize8_y(w,nx,ny,nz+1,nnn,nn)
  endif
endsubroutine vertical_velocity
