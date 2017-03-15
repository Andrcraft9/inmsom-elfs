!=====================================================================
subroutine richnum(den,uu,vv,rit)
use main_basin_pars
use mpi_parallel_tools
use basin_grid
implicit none
!----------------------------------------------------------------------
! richardson' number on u-grid. (d uv/d z .ne. 0). k=2,nz.
! richnu = g/r0*(d den/d z) / ((d u/d z)^2+(d v/d z)^2), d z = hzt*hh
! hzt(k) = z(k)-z(k-1), k=2,nz
! rit - richardson number on W-grid

      real(8) den(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),       &
               uu(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),       &
               vv(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),       &
              rit(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1)
      real(8) u1,u2,v1,v2,r1,r2
      integer m, n, k

!$omp parallel do private(m,n,k,u2,v2,u1,v1,r1,r2)
      do n=ny_start,ny_end
       do m=nx_start,nx_end
        if (lu(m,n)>0.5) then
          r2=den(m,n,1)
          u2 = (uu(m  ,n,1)*dxt(m  ,n)*dyh(m  ,n)*hhu(m  ,n)    &
              + uu(m-1,n,1)*dxt(m-1,n)*dyh(m-1,n)*hhu(m-1,n))/2.0d0/hhq(m,n)/dx(m,n)/dy(m,n)       !u on t-grid
          v2 = (vv(m,n  ,1)*hhv(m,n  )*dxh(m,n  )*dyt(m,n  )    &
              + vv(m,n-1,1)*hhv(m,n-1)*dxh(m,n-1)*dyt(m,n-1))/2.0d0/hhq(m,n)/dx(m,n)/dy(m,n)       !v on t-grid

	   do k=2,nz
          r1 = r2  
          u1 = u2
          v1 = v2

          r2 = den(m,n,k)
          u2 = (uu(m  ,n,k)*dxt(m  ,n)*dyh(m  ,n)*hhu(m  ,n)    &
              + uu(m-1,n,k)*dxt(m-1,n)*dyh(m-1,n)*hhu(m-1,n))/2.0d0/hhq(m,n)/dx(m,n)/dy(m,n)       !u on t-grid
          v2 = (vv(m,n  ,k)*hhv(m,n  )*dxh(m,n  )*dyt(m,n  )    &
              + vv(m,n-1,k)*hhv(m,n-1)*dxh(m,n-1)*dyt(m,n-1))/2.0d0/hhq(m,n)/dx(m,n)/dy(m,n)       !v on t-grid

! direct discrete form of Ri
          rit(m,n,k)=min(FreeFallAcc/RefDen*(r2-r1)*hzt(k)*hhq(m,n)     &
                           / ( max((u2-u1)**2+(v2-v1)**2,0.0025d-4) ) ,1d+05)
         enddo

        endif
       enddo
      enddo
!$omp end parallel do

endsubroutine richnum
!========================================================
subroutine ppmix(rit,anzu,anumaxu,anubgru,anzt,anumaxt,anubgrt)
use main_basin_pars
use mpi_parallel_tools
use basin_grid
implicit none
   
 real(8) dimdepth1, dimdepth2, undimdepth1, undimdepth2
 parameter(dimdepth1=3.0d0,dimdepth2=20000.0d0)

integer m,n,k
real(8) rit(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1),  &
       anzu(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1),  &
       anzt(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1)
real(8) anumaxu,anubgru,anumaxt,anubgrt

!$omp parallel do private(m,n,k,undimdepth1,undimdepth2)
      do n=ny_start,ny_end
       do m=nx_start,nx_end

        if (lu(m,n)>0.5) then
          undimdepth1 = dimdepth1/hhq(m,n)          
          undimdepth2 = dimdepth2/hhq(m,n) 
          do k=2,nz
           
           if(zw(k)<undimdepth2) then
            
            if(rit(m,n,k)>0.0.and.zw(k)>undimdepth1) then
             anzt(m,n,k) = (anumaxt-anubgrt)/(1.0+5.0*rit(m,n,k))**2 + anubgrt
             anzu(m,n,k) = (anumaxu-anubgru)/(1.0+5.0*rit(m,n,k))**3 + anubgru
            else
             anzt(m,n,k) = anumaxt
             anzu(m,n,k) = anumaxu
            end if
           
           else
             anzt(m,n,k) = anubgrt
            if(rit(m,n,k)>0.0) then
             anzu(m,n,k) = (anumaxu-anubgru)/(1.0+5.0*rit(m,n,k))**3  + anubgru
            else
             anzu(m,n,k) = anumaxu
            endif
           endif
          enddo
        endif

       enddo
      enddo
!$omp end parallel do

      if(periodicity_x/=0) then 
        call cyclize8_x(anzu,nx,ny,nz+1,mmm,mm)
      endif

      if(periodicity_y/=0) then 
        call cyclize8_y(anzu,nx,ny,nz+1,nnn,nn)
      endif

endsubroutine ppmix

!====================================================================================
subroutine stress_components(u,v,str_t,str_s,nlev)
 use main_basin_pars
 use mpi_parallel_tools
 use basin_grid
 implicit none

 integer nlev
 real(8) u(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nlev),    &
         v(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nlev),    &
     str_t(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nlev),    &
     str_s(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nlev)
 
 integer m,n,k

 !$omp parallel do private(m,n,k) 
 do n=ny_start, ny_end
   do m=nx_start, nx_end
    if(lu(m,n)>0.5) then
     do k=1,nlev 
      str_t(m,n,k)=dy(m,n)/dx(m,n)*(u(m,n,k)/dyh(m,n)-u(m-1,n,k)/dyh(m-1,n))     &
                  -dx(m,n)/dy(m,n)*(v(m,n,k)/dxh(m,n)-v(m,n-1,k)/dxh(m,n-1)) 
     enddo
    endif
   enddo
 enddo
!$omp end parallel do

!$omp parallel do private(m,n,k) 
 do n=ny_start, ny_end
   do m=nx_start, nx_end
    if(luu(m,n)>0.5) then
     do k=1,nlev
      str_s(m,n,k)=dxb(m,n)/dyb(m,n)*(u(m,n+1,k)/dxt(m,n+1)-u(m,n,k)/dxt(m,n))     &
                  +dyb(m,n)/dxb(m,n)*(v(m+1,n,k)/dyt(m+1,n)-v(m,n,k)/dyt(m,n))    
     enddo
    endif
   enddo
 enddo
!$omp end parallel do

      if(periodicity_x/=0) then
       call cyclize8_x(str_t,nx,ny,nlev,mmm,mm)
       call cyclize8_x(str_s,nx,ny,nlev,mmm,mm)
	end if

      if(periodicity_y/=0) then
       call cyclize8_y(str_t,nx,ny,nlev,nnn,nn)
       call cyclize8_y(str_s,nx,ny,nlev,nnn,nn)
	end if

endsubroutine stress_components

!======================================================================
subroutine smagorinsky_coeff(visc2, dift, str_t,str_s, amts, amuv)
 use main_basin_pars
 use mpi_parallel_tools
 use basin_grid
 implicit none

 real(8), parameter:: prandtl=0.1d0
 real(8) str_s2, mu_smag
 real(8) visc2, dift
 real(8)   str_t(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),   &
           str_s(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),   &
           amts(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),    &
           amuv(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz)

 integer m,n,k
 

 !$omp parallel do private(m,n,k,str_s2,mu_smag) 
 do n=ny_start, ny_end
   do m=nx_start, nx_end

    if(lu(m,n)>0.5) then
     do k=1,nz
      str_s2= ( str_s(m  ,n  ,k)**2 +str_s(m-1,n  ,k)**2    &
               +str_s(m  ,n-1,k)**2 +str_s(m-1,n-1,k)**2 )/4.0d0     
      mu_smag=(vonkarman/Pi)**2 *dx(m,n)*dy(m,n) *dsqrt(str_t(m,n,k)**2 + str_s2)  
      amuv(m,n,k) =visc2+ mu_smag
      amts(m,n,k) =dift + mu_smag*prandtl
     enddo
    endif

   enddo
 enddo
!$omp end parallel do

      if(periodicity_x/=0) then
       call cyclize8_x(amuv ,nx,ny,nz,mmm,mm)
       call cyclize8_x(amts ,nx,ny,nz,mmm,mm)
	end if

      if(periodicity_y/=0) then
       call cyclize8_y(amuv ,nx,ny,nz,nnn,nn)
       call cyclize8_y(amts ,nx,ny,nz,nnn,nn)
	end if

endsubroutine smagorinsky_coeff

!======================================================================
subroutine diffusion_slopes(den,slrx,slry,slzx,slzy)
 use main_basin_pars
 use mpi_parallel_tools
 use basin_grid
 implicit none
! calculating slope for isopycnal diffusion
! den  -- density
! periodicity is taken to account by masks lu, lcu and lcv
! vertical density gradient is agreement with the operator of isopycnal diffusion
! epsrho -- min vertical density gradient (commonly 1.0e-12 g/cm3)
      real(8)    epsrho, slope,angle_max,slope_bot
      parameter(epsrho=1.0d-5,angle_max=0.01d0)

      real(8) den(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),      &
             slrx(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1),      &
             slry(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1),      &
             slzx(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1),      &
             slzy(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1)
      integer m,n,k
      real(8) rhoz


!$omp parallel do private(m,n,k,rhoz,slope,slope_bot)
     do n=ny_start,ny_end
      do m=nx_start,nx_end

        if(lcu(m,n)>0.5) then
! vertical density gradient on xc-w-grid as mean on t-points

          do k=2,nz

            rhoz=(den(m,n,k  )+den(m+1,n,k  )     &
                - den(m,n,k-1)-den(m+1,n,k-1)) /hhu(m,n)/hzt(k)/2.0d0
            rhoz=max(rhoz,epsrho)
            
            slope=    ( den(m+1,n,k  ) - den(m,n,k  )     &
                      + den(m+1,n,k-1) - den(m,n,k-1)) /2.0d0/rhoz
            slope_bot=( hhq(m+1,n) - hhq(m,n))*zw(k) 

            slrx(m,n,k) = max( min(slope,slope_bot+angle_max*dxt(m,n)), slope_bot-angle_max*dxt(m,n) )
            slzx(m,n,k) = slope_bot            
          end do               
         end if

        if(lcv(m,n)>0.5) then
! vertical density gradient on yc-w-grid as mean on t-points

          do k=2,nz

            rhoz=(den(m,n,k  )+den(m,n+1,k  )     &
                - den(m,n,k-1)-den(m,n+1,k-1))/hhv(m,n)/hzt(k)/2.0d0
            rhoz=max(rhoz,epsrho)
            
            slope=    ( den(m,n+1,k  ) - den(m,n,k  )     &
                      + den(m,n+1,k-1) - den(m,n,k-1)) /2.0d0/rhoz
            slope_bot=( hhq(m,n+1) - hhq(m,n))*zw(k) 

            slry(m,n,k) = max( min(slope,slope_bot+angle_max*dyt(m,n)), slope_bot-angle_max*dyt(m,n) )
            slzy(m,n,k) = slope_bot            
          end do               
         end if
       

        end do
     end do
!$omp end parallel do
      
	if(periodicity_x/=0) then
	 call cyclize8_x(slrx,nx,ny,nz,mmm,mm)
	 call cyclize8_x(slzx,nx,ny,nz,mmm,mm)
      end if

	if(periodicity_y/=0) then
	 call cyclize8_y(slry,nx,ny,nz,nnn,nn)
	 call cyclize8_y(slzy,nx,ny,nz,nnn,nn)                
      end if

endsubroutine diffusion_slopes

!========================================================================
subroutine turb_tran_diff(ff,     &
                         ffp,     &
                         tau,     &
                          uu,     &
                          vv,     &
                          ww,     &
                          mu,     &
                   factor_mu,     &
                          nu,     &
                   factor_nu,     &
                      ff_top,     &
                      ff_bot,     &
                       gendis    )
use main_basin_pars
use mpi_parallel_tools
use basin_grid
implicit none

real(8) tau, factor_mu, factor_nu

real(8) ff(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1),    &   !Transported tracer
       ffp(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1),    &   !Transported tracer (previous step)
        uu(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),      &
        vv(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),      &
        ww(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1)         !Transporting velocities

real(8) mu(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz  ),    &
        nu(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1)

real(8) gendis(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1)       !Internal source (divergence of SW-radiation)

real(8) ff_top(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
        ff_bot(bnd_x1:bnd_x2,bnd_y1:bnd_y2)

real(8) fz_p,fz_m   !fluxes through cell edges
real(8) flux_adv, flux_diff

integer m,n,k

real(8) a(nz+1),b(nz+1),c(nz+1),eta(nz+1),rksi(nz+1) 
real(8) bp, bp0, dp, dm, nu_p, nu_m,rhs 

real(8), allocatable::  flux_x(:,:,:), flux_y(:,:,:),      &
                        fn(:,:,:)

   allocate(flux_x(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1),     &
            flux_y(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1),     &
                fn(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1) ) 

   flux_x=0.0d0
   flux_y=0.0d0
       fn=0.0d0


!$omp parallel do private(m, n, k, flux_adv, flux_diff) 
   do n=ny_start,ny_end
    do m=nx_start,nx_end    
     
     if(lcu(m,n)>0.5) then
      do k=2,nz
       flux_adv = -(uu(m,n,k)*dz(k)+uu(m,n,k-1)*dz(k-1))/2.0d0/hzt(k)    &
                   *dyh(m,n)*hhu(m,n)*(ff(m,n,k)+ff(m+1,n,k))/2.0d0
       flux_diff=(mu(m,n,k)+mu(m,n,k-1)+mu(m+1,n,k)+mu(m+1,n,k-1))/4.0d0*factor_mu  &
                 *dyh(m,n)/dxt(m,n)*hhu(m,n)*(ff(m+1,n,k)-ff(m,n,k))
       flux_x(m,n,k)=flux_adv+flux_diff
      enddo
     endif

     if(lcu(m,n)>0.5) then
      do k=2,nz
       flux_adv = -(vv(m,n,k)*dz(k)+vv(m,n,k-1)*dz(k-1))/2.0d0/hzt(k)    &
                   *dxh(m,n)*hhv(m,n)*(ff(m,n,k)+ff(m,n+1,k))/2.0d0
       flux_diff=(mu(m,n,k)+mu(m,n,k-1)+mu(m,n+1,k)+mu(m,n+1,k-1))/4.0d0*factor_mu  &
                 *dxh(m,n)/dyt(m,n)*hhv(m,n)*(ff(m,n+1,k)-ff(m,n,k))
       flux_y(m,n,k)=flux_adv+flux_diff
      enddo
     endif

    enddo
   enddo
!$omp end parallel do        

	if(periodicity_x/=0) then
	 call cyclize8_x( flux_x,nx,ny,nz+1,mmm,mm)
      end if

	if(periodicity_y/=0) then
	 call cyclize8_y( flux_y,nx,ny,nz+1,nnn,nn)            
      end if

!$omp parallel do private(m,n,k,bp,bp0,dp,dm,a,b,c,eta,rksi,fz_p,fz_m, nu_p, nu_m, rhs)
      do n=ny_start,ny_end
       do m=nx_start,nx_end
        if (lu(m,n)>0.5) then
        
         bp =  hhqn(m,n)*dx(m,n)*dy(m,n) / tau/2.0d0
         bp0=  hhqp(m,n)*dx(m,n)*dy(m,n) / tau/2.0d0

!Upper layer          
          k=2
          
          nu_m=(nu(m,n,k)+nu(m,n,k-1))/2.0d0*factor_nu
          nu_p=(nu(m,n,k)+nu(m,n,k+1))/2.0d0*factor_nu
          
          fz_m= -(ww(m,n,k-1)+ww(m,n,k))/2.0d0*dx(m,n)*dy(m,n)*(ff(m,n,k)+ff_top(m,n))/2.0d0
          fz_p= -(ww(m,n,k+1)+ww(m,n,k))/2.0d0*dx(m,n)*dy(m,n)*(ff(m,n,k)+ff(m,n,k+1))/2.0d0

          dm = nu_m*dx(m,n)*dy(m,n) /hhqn(m,n)/dz(k-1)
          dp = nu_p*dx(m,n)*dy(m,n) /hhqn(m,n)/dz(k  )

          rhs =  flux_x(m,n,k) - flux_x(m-1,n,k) + flux_y(m,n,k) - flux_y(m,n-1,k) + (fz_p - fz_m)/hzt(k)    &           
                              + gendis(m,n,k)*dx(m,n)*dy(m,n) 

! "2." above is due to aproximation.
          c(k) =  - dp/hzt(k)
          a(k) =  0.0d0
          
          b(k) =  bp + (dp + dm)/hzt(k)
          eta(k) = bp0*ffp(m,n,k) + rhs + dm*ff_top(m,n)/hzt(k) 

! internal points.
         do k=3,nz-1
          
          nu_p=(nu(m,n,k)+nu(m,n,k+1))/2.0d0*factor_nu
          
          fz_m=fz_p
          fz_p= -(ww(m,n,k+1)+ww(m,n,k))/2.0d0*dx(m,n)*dy(m,n)*(ff(m,n,k)+ff(m,n,k+1))/2.0d0

          dm = dp
          dp = nu_p*dx(m,n)*dy(m,n) /hhqn(m,n)/dz(k  )       
          
          rhs =  flux_x(m,n,k) - flux_x(m-1,n,k) + flux_y(m,n,k) - flux_y(m,n-1,k) + (fz_p - fz_m)/hzt(k)    &           
                              + gendis(m,n,k)*dx(m,n)*dy(m,n)

! "2." above is due to aproximation.
          c(k) =  - dp/hzt(k)
          a(k) =  - dm/hzt(k)
          b(k) =  bp + (dp + dm)/hzt(k)

          eta(k) = bp0*ffp(m,n,k) + rhs
         enddo

!Bottom layer
          k=nz

          nu_p=nu(m,n,k+1)*factor_nu

          fz_m= fz_p
          fz_p= -(ww(m,n,k+1)+ww(m,n,k))/2.0d0*dx(m,n)*dy(m,n)*(ff(m,n,k)+ff_bot(m,n))/2.0d0

          dm = dp
          dp = nu_p*dx(m,n)*dy(m,n) /hhqn(m,n)/dz(k  )  
          
          rhs =  flux_x(m,n,k) - flux_x(m-1,n,k) + flux_y(m,n,k) - flux_y(m,n-1,k) + (fz_p - fz_m)/hzt(k)    &           
                              + gendis(m,n,k)*dx(m,n)*dy(m,n)

! "2." above is due to aproximation.
          c(k) =  0.0d0
          a(k) =  - dm/hzt(k)
          
            b(k) =  bp + (dp + dm)/hzt(k)
          eta(k) = bp0*ffp(m,n,k) + rhs + dp*ff_bot(m,n)/hzt(k) 
         
         call factor8(nz+1,a,b,c,eta,rksi,2,nz)
         do k=2,nz
          fn(m,n,k)=rksi(k)
         enddo
        endif
       enddo
      enddo
!$omp end parallel do

	if(periodicity_x/=0) then
	 call cyclize8_x(fn,nx,ny,nz,mmm,mm)
      end if

	if(periodicity_y/=0) then
	 call cyclize8_y(fn,nx,ny,nz,nnn,nn)
      end if


!$omp parallel do private(m,n,k)
      do n=ny_start-1,ny_end+1
       do m=nx_start-1,nx_end+1
        
        if(lu(m,n)>0.5) then
         do k=2,nz
!          ffp(m,n,k)=hhq(m,n)*ff(m,n,k)+time_smooth*(hhqn(m,n)*fn(m,n,k)-2.0d0*hhq(m,n)*ff(m,n,k)+hhqp(m,n)*ffp(m,n,k))/2.0d0
          ffp(m,n,k)=max(ff(m,n,k)+time_smooth*(fn(m,n,k)-2.0d0*ff(m,n,k)+ffp(m,n,k))/2.0d0,tur_var_min)
           ff(m,n,k)=max(fn(m,n,k),tur_var_min)
         enddo
        endif

       enddo
      enddo
!$omp end parallel do

   deallocate(fn,flux_y,flux_x)

endsubroutine turb_tran_diff

!========================================================================
subroutine turb_gen_diss(q2,     &
                        q2p,     &
                        q2l,     &
                       q2lp,     &
                         uu,     &
                         vv,     &
                        den,     &
                       anzt,     &
                       anzu,     &
                    anubgru,     &
                    anubgrt,     &
                     rhs_q2,     &
                     rhs_q2l     )
use main_basin_pars
use mpi_parallel_tools
use basin_grid
implicit none

real(8)    q2(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1),     &
          q2p(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1),     &
          q2l(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1),     &
         q2lp(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1),     &
           uu(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz  ),     &
           vv(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz  ),     &
          den(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz  ),     &
         anzt(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1),     &
         anzu(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1),     &
       rhs_q2(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1),     &
       rhs_q2l(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1)     
real(8) anubgrt, anubgru

real(8) shf2, vbf2, q, lscl, rich
real(8) stab_t, stab_u
real(8) wall, lm1, dpth
real(8) l_max
real(8) u1,u2,v1,v2,r1,r2
integer m, n, k


!$omp parallel do private(m, n, k, u1, u2, v1, v2, r1, r2, shf2, vbf2, q, l_max, lscl, rich, stab_t, stab_u,    &
!$omp                     wall, lm1, dpth ) 
   do n=ny_start,ny_end
    do m=nx_start,nx_end    
     
     if(lu(m,n)>0.5) then
     
        r2=den(m,n,1)
        u2 = (uu(m  ,n,1)*dxt(m  ,n)*dyh(m  ,n)*hhu(m  ,n)    &
            + uu(m-1,n,1)*dxt(m-1,n)*dyh(m-1,n)*hhu(m-1,n))/2.0d0/hhq(m,n)/dx(m,n)/dy(m,n)       !u on t-grid
        v2 = (vv(m,n  ,1)*hhv(m,n  )*dxh(m,n  )*dyt(m,n  )    &
            + vv(m,n-1,1)*hhv(m,n-1)*dxh(m,n-1)*dyt(m,n-1))/2.0d0/hhq(m,n)/dx(m,n)/dy(m,n)       !v on t-grid
     
      do k=2,nz
        r1 = r2  
        u1 = u2
        v1 = v2

        r2 = den(m,n,k)
        u2 = (uu(m  ,n,k)*dxt(m  ,n)*dyh(m  ,n)*hhu(m  ,n)    &
            + uu(m-1,n,k)*dxt(m-1,n)*dyh(m-1,n)*hhu(m-1,n))/2.0d0/hhq(m,n)/dx(m,n)/dy(m,n)       !u on t-grid
        v2 = (vv(m,n  ,k)*hhv(m,n  )*dxh(m,n  )*dyt(m,n  )    &
            + vv(m,n-1,k)*hhv(m,n-1)*dxh(m,n-1)*dyt(m,n-1))/2.0d0/hhq(m,n)/dx(m,n)/dy(m,n)       !v on t-grid
      
        shf2= ( (u2-u1)**2+(v2-v1)**2 )/(hzt(k)*hhq(m,n))**2
        
        vbf2= FreeFallAcc/RefDen *(r2-r1)/(hzt(k)*hhq(m,n))

        q=dsqrt(q2p(m,n,k))
        
        l_max= 0.53d0 * q/ dsqrt(max(vbf2,1.0d-10))

        lscl=min(q2lp(m,n,k)/q2p(m,n,k),l_max)
        
        rich=min(-vbf2*lscl**2/q2p(m,n,k),0.028d0)

        stab_t=A2_t*(1.0d0-6.0d0*A1_t/B1_t)     &
               / (1.0d0 - (3.0d0*A2_t*B2_t+18.0d0*A1_t*A2_t)*rich)
        
        stab_u=( A1_t*(1.0d0-3.0d0*C1_t-6.0d0*A1_t/B1_t)          &
                +stab_t*((18.0d0*A1_t**2+9.0d0*A1_t*A2_t)*rich) ) &
                /(1.0d0-9.0d0*A1_t*A2_t*rich)

        anzt(m,n,k)=anubgrt + q*lscl*stab_t
        anzu(m,n,k)=anubgru + q*lscl*stab_u

        dpth=zw(k)*hhq(m,n)

        lm1=1.0d0/dpth + 1.0d0/(hhq(m,n)-dpth)

        wall=1.0d0+E2_t*(lscl*lm1/vonkarman)**2

        rhs_q2(m,n,k) =    2.0d0* hhq(m,n)*(anzu(m,n,k)*shf2 - anzt(m,n,k)*vbf2      &
                                - q2p(m,n,k)*q/B1_t/lscl)

        rhs_q2l(m,n,k)= E1_t*lscl*hhq(m,n)*(anzu(m,n,k)*shf2 - anzt(m,n,k)*vbf2)     &
                                - hhq(m,n)*q2p(m,n,k)*q/B1_t*wall

      enddo
     
     endif

    enddo
   enddo
!$omp end parallel do        

      if(periodicity_x/=0) then 
        call cyclize8_x(anzu,nx,ny,nz+1,mmm,mm)
      endif

      if(periodicity_y/=0) then 
        call cyclize8_y(anzu,nx,ny,nz+1,nnn,nn)
      endif

endsubroutine turb_gen_diss
