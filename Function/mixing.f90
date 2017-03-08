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
! rit - richardson number on t-grid

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

 real(8), parameter:: vonkarman=0.4d0, prandtl=0.1d0
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
