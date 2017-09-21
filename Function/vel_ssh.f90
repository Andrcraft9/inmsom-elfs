!===============================================================================
! RHS for implicit bfc scheme
subroutine uv_bfc(u, v, hq, hu, hv, hh, RHSx, RHSy)
    use main_basin_pars
    use mpi_parallel_tools
    use basin_grid
    implicit none

    real(8) u(bnd_x1:bnd_x2,bnd_y1:bnd_y2),       & !Transporting zonal velocity
            v(bnd_x1:bnd_x2,bnd_y1:bnd_y2)          !Transporting meridional velocity

    real(8) RHSx(bnd_x1:bnd_x2,bnd_y1:bnd_y2),    &
            RHSy(bnd_x1:bnd_x2,bnd_y1:bnd_y2)

    real(8) hq(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &
            hu(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &
            hv(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &
            hh(bnd_x1:bnd_x2,bnd_y1:bnd_y2)

    integer :: m, n
    real*8 :: k_bfc, s
    real*8 :: k1, k2

    !$omp parallel do private(m, n, k_bfc, s, k1, k2)
    do n=ny_start, ny_end
        do m=nx_start, nx_end
            if (lcu(m,n)>0.5) then
                ! Discretization in h-points
                k_bfc = FreeFallAcc * (nbfc**2) / (hh(m, n)**(1.0/3.0))
                s = 0.5d0 * sqrt( (u(m, n) + u(m, n+1))**2 + (v(m, n) + v(m+1, n))**2 )
                k1 = -dxb(m, n) * dyb(m, n) * k_bfc * s
                !k1 = k1 * 0.5d0*(u(m, n) + u(m, n+1))

                ! Discretization in h-points
                k_bfc = FreeFallAcc * (nbfc**2) / (hh(m, n-1)**(1.0/3.0))
                s = 0.5d0 * sqrt( (u(m, n) + u(m, n-1))**2 + (v(m, n-1) + v(m+1, n-1))**2 )
                k2 = -dxb(m, n-1) * dyb(m, n-1) * k_bfc * s
                !k2 = k2 * 0.5d0*(u(m, n) + u(m, n-1))

                ! Discretization in u-points
                RHSx(m, n) = 0.5d0 * (k1 + k2)
             endif

             if (lcv(m,n)>0.5) then
                ! Discretization in h-points
                k_bfc = FreeFallAcc * (nbfc**2) / (hh(m, n)**(1.0/3.0))
                s = 0.5d0 * sqrt( (u(m, n) + u(m, n+1))**2 + (v(m, n) + v(m+1, n))**2 )
                k1 = -dxb(m, n) * dyb(m, n) * k_bfc * s
                !k1 = k1 * 0.5d0*(v(m, n) + v(m+1, n))

                ! Discretization in h-points
                k_bfc = FreeFallAcc * (nbfc**2) / (hh(m-1, n)**(1.0/3.0))
                s = 0.5d0 * sqrt( (u(m-1, n) + u(m-1, n+1))**2 + (v(m, n) + v(m-1, n))**2 )
                k2 = -dxb(m-1, n) * dyb(m-1, n) * k_bfc * s
                !k2 = k2 * 0.5d0*(v(m, n) + v(m-1, n))

                ! Discretization in v-points
                RHSy(m, n) = 0.5d0 * (k1 + k2)
             endif
        enddo
    enddo
    !$omp end parallel do

end subroutine uv_bfc

!===========================================================================================
subroutine uv_trans( u, v, vort,    &
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
         RHSy(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nlev )         !meridional source function


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

      call syncborder_real8(vort, nlev)

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
             - ( vort(m  ,n,k)*hh(m  ,n)*(u(m  ,n+1,k)+u(m  ,n,k))               &
             +   vort(m-1,n,k)*hh(m-1,n)*(u(m-1,n+1,k)+u(m-1,n,k))  )/4.0d0
        end do

      end if

    end do
  end do
!$omp end parallel do

!  call syncborder_real8(RHSx, nlev)
!  call syncborder_real8(RHSy, nlev)

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

!  call syncborder_real8(RHSx, nlev)
!  call syncborder_real8(RHSy, nlev)

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

      call syncborder_real8(fx, nlev)
      call syncborder_real8(fy, nlev)

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

!  call syncborder_real8(RHSx, nlev)
!  call syncborder_real8(RHSy, nlev)

endsubroutine uv_diff4

!===============================================================================
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
                               ssh_i,   &
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
                             RHSy_dif,  &
                             RHSx_bfc,  &
                             RHSy_bfc)
    use main_basin_pars
    use mpi_parallel_tools
    use basin_grid
    implicit none

    real(8) tau, tau_inner
    integer step, nstep, ksw4, m, n

    real(8)  ubrtr_e(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
            ubrtrp_e(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
             vbrtr_e(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
            vbrtrp_e(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
               ssh_e(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
              sshp_e(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
             ubrtr_i(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
             vbrtr_i(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
               ssh_i(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
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
            RHSy_dif(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
            RHSx_bfc(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
            RHSy_bfc(bnd_x1:bnd_x2,bnd_y1:bnd_y2)

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

    real*8 time_count
    integer ierr

    allocate(u(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &
             v(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &
           ssh(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &
            up(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &
            vp(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &
          sshp(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &
            un(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &
            vn(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &
            sshn(bnd_x1:bnd_x2,bnd_y1:bnd_y2))

    tau_inner=tau/dfloat(nstep)

    un = 0.0d0
    vn = 0.0d0
    sshn = 0.0d0

    up=ubrtrp_e
    vp=vbrtrp_e
    sshp=sshp_e

    u = ubrtr_e
    v = vbrtr_e
    ssh = ssh_e

    ubrtr_i = 0.0d0
    vbrtr_i = 0.0d0
    ssh_i = 0.0d0

    do step = 1, nstep
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

        call syncborder_real8(sshn, 1)
        if(periodicity_x/=0) then
            call cyclize8_x(sshn,nx,ny,1,mmm,mm)
        endif
        if(periodicity_y/=0) then
            call cyclize8_y(sshn,nx,ny,1,nnn,nn)
        endif

        if(full_free_surface>0) then
            call hh_update(hhqn_e, hhun_e, hhvn_e, hhhn_e, sshn, hhq_rest)
        endif

        !computing advective and lateral-viscous terms for 2d-velocity
        ! call stress_components(up,vp,str_t2d,str_s2d,1)

        !computing advective and lateral-viscous terms for 2d-velocity
        call uv_trans( u, v, vort,            &
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

        ! compute BottomFriction (bfc)
        !call uv_bfc(up, vp, hhq_e, hhu_e, hhv_e, hhh_e, RHSx_bfc, RHSy_bfc)

        !$omp parallel do private(bp,bp0,grx,gry, slx, sly, slxn, slyn)
        do n=ny_start,ny_end
            do m=nx_start,nx_end
                !zonal flux
                if(lcu(m,n)>0.5) then

                    bp =hhun_e(m,n)*dxt(m,n)*dyh(m,n)/2.0d0/tau_inner
                    bp0=hhup_e(m,n)*dxt(m,n)*dyh(m,n)/2.0d0/tau_inner

                    slx = - FreeFallAcc * ( ssh(m+1,n) - ssh(m,n))*dyh(m,n)* hhu_e(m,n)
                    !slxn= - FreeFallAcc * (sshn(m+1,n) -sshn(m,n))*dyh(m,n)*hhun_e(m,n)
                    grx= RHSx(m,n) + slx  + RHSx_dif(m,n) + RHSx_adv(m,n)      &
                        - (rdis(m,n)+rdis(m+1,n))/2.0d0 *up(m,n)*dxt(m,n)*dyh(m,n)*hhu_e(m,n)        &
                        + ( rlh_s(m,n  )*hhh_e(m,n  )*dxb(m,n  )*dyb(m,n  )*(v(m+1,n  )+v(m,n  ))             &
                        + rlh_s(m,n-1)*hhh_e(m,n-1)*dxb(m,n-1)*dyb(m,n-1)*(v(m+1,n-1)+v(m,n-1))  )/4.0d0

                    un(m,n)=(up(m,n)*bp0 + grx )/(bp - RHSx_bfc(m, n))
                endif

                !meridional flux
                if(lcv(m,n)>0.5) then

                    bp =hhvn_e(m,n)*dyt(m,n)*dxh(m,n)/2.0d0/tau_inner
                    bp0=hhvp_e(m,n)*dyt(m,n)*dxh(m,n)/2.0d0/tau_inner

                    sly = - FreeFallAcc * ( ssh(m,n+1)- ssh(m,n))*dxh(m,n)* hhv_e(m,n)
                    !slyn= - FreeFallAcc * (sshn(m,n+1)-sshn(m,n))*dxh(m,n)*hhvn_e(m,n)
                    gry= RHSy(m,n) + sly  + RHSy_dif(m,n) + RHSy_adv(m,n)      &
                        - (rdis(m,n)+rdis(m,n+1))/2.0d0 *vp(m,n)*dxh(m,n)*dyt(m,n)*hhv_e(m,n)        &
                        - ( rlh_s(m  ,n)*hhh_e(m  ,n)*dxb(m  ,n)*dyb(m  ,n)*(u(m  ,n+1)+u(m  ,n))             &
                        + rlh_s(m-1,n)*hhh_e(m-1,n)*dxb(m-1,n)*dyb(m-1,n)*(u(m-1,n+1)+u(m-1,n))  )/4.0d0

                    vn(m,n)=(vp(m,n)*bp0 + gry )/(bp - RHSy_bfc(m, n))
                endif
            enddo
        enddo
        !$omp end parallel do

        call syncborder_real8(un ,1)
        call syncborder_real8(vn, 1)
        if(periodicity_x/=0) then
            call cyclize8_x(un,nx,ny,1,mmm,mm)
            call cyclize8_x(vn,nx,ny,1,mmm,mm)
        endif
        if(periodicity_y/=0) then
            call cyclize8_y(un,nx,ny,1,nnn,nn)
            call cyclize8_y(vn,nx,ny,1,nnn,nn)
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
                    !up(m,n) =  hhu_e(m,n)*u(m,n)+time_smooth*(hhun_e(m,n)*un(m,n)-2.0d0*hhu_e(m,n)*u(m,n)+hhup_e(m,n)*up(m,n))/2.0d0/dfloat(nstep)
                    up(m,n) =  u(m,n)+time_smooth*(un(m,n)-2.0d0*u(m,n)+up(m,n))/2.0d0/dfloat(nstep)
                    u(m,n) = un(m,n)
                endif
                if(lcv(m,n)>0.5) then
                    !vp(m,n) =  hhv_e(m,n)*v(m,n)+time_smooth*(hhvn_e(m,n)*vn(m,n)-2.0d0*hhv_e(m,n)*v(m,n)+hhvp_e(m,n)*vp(m,n))/2.0d0/dfloat(nstep)
                    vp(m,n) =  v(m,n)+time_smooth*(vn(m,n)-2.0d0*v(m,n)+vp(m,n))/2.0d0/dfloat(nstep)
                    v(m,n) = vn(m,n)
                endif

            enddo
        enddo
        !$omp end parallel do

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

            ssh_i = ssh
            ubrtr_i = u
            vbrtr_i = v
        endif

    enddo

    !call syncborder_real8(ubrtr_i, 1)
    !call syncborder_real8(vbrtr_i, 1)
    !if(periodicity_x/=0) then
    !    call cyclize8_x(ubrtr_i,nx,ny,1,mmm,mm)
    !    call cyclize8_x(vbrtr_i,nx,ny,1,mmm,mm)
    !endif
    !if(periodicity_y/=0) then
    !    call cyclize8_y(ubrtr_i,nx,ny,1,nnn,nn)
    !    call cyclize8_y(vbrtr_i,nx,ny,1,nnn,nn)
    !endif

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
