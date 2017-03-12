!============================================================
subroutine hh_init(hq, hqp, hqn,    &
                   hu, hup, hun,    &
                   hv, hvp, hvn,    &
                   hh, hhp, hhn,    &
                   sh, shp, h_r)
 use main_basin_pars
 use mpi_parallel_tools
 use basin_grid
 implicit none

 real(8) hq(bnd_x1:bnd_x2, bnd_y1:bnd_y2), hqp(bnd_x1:bnd_x2, bnd_y1:bnd_y2), hqn(bnd_x1:bnd_x2, bnd_y1:bnd_y2),    &
         hu(bnd_x1:bnd_x2, bnd_y1:bnd_y2), hup(bnd_x1:bnd_x2, bnd_y1:bnd_y2), hun(bnd_x1:bnd_x2, bnd_y1:bnd_y2),    &
         hv(bnd_x1:bnd_x2, bnd_y1:bnd_y2), hvp(bnd_x1:bnd_x2, bnd_y1:bnd_y2), hvn(bnd_x1:bnd_x2, bnd_y1:bnd_y2),    &
         hh(bnd_x1:bnd_x2, bnd_y1:bnd_y2), hhp(bnd_x1:bnd_x2, bnd_y1:bnd_y2), hhn(bnd_x1:bnd_x2, bnd_y1:bnd_y2),    &
         sh(bnd_x1:bnd_x2, bnd_y1:bnd_y2), shp(bnd_x1:bnd_x2, bnd_y1:bnd_y2), h_r(bnd_x1:bnd_x2, bnd_y1:bnd_y2)

 real(8) slu
 integer m,n

       hq =h_r + sh *dfloat(full_free_surface)
       hqp=h_r + shp*dfloat(full_free_surface)
       hqn=h_r

!$omp parallel do private(m,n,slu)
      do n=ny_start-1,ny_end
       do m=nx_start-1,nx_end

        if(llu(m,n)>0.5) then
! interpolating hhq given on T-grid(lu) to hhu given on u-grid(lcu).
          slu=dble(lu(m,n)+lu(m+1,n))
          hu(m,n)=( hq(m  ,n)*dx(m  ,n)*dy(m  ,n)*dble(lu(m  ,n))   &
                  + hq(m+1,n)*dx(m+1,n)*dy(m+1,n)*dble(lu(m+1,n)) )/slu/dxt(m,n)/dyh(m,n)
         hup(m,n)=( hqp(m  ,n)*dx(m  ,n)*dy(m  ,n)*dble(lu(m  ,n))   &
                  + hqp(m+1,n)*dx(m+1,n)*dy(m+1,n)*dble(lu(m+1,n)) )/slu/dxt(m,n)/dyh(m,n)
         hun(m,n)=( hqn(m  ,n)*dx(m  ,n)*dy(m  ,n)*dble(lu(m  ,n))   &
                  + hqn(m+1,n)*dx(m+1,n)*dy(m+1,n)*dble(lu(m+1,n)) )/slu/dxt(m,n)/dyh(m,n)
        endif

        if(llv(m,n)>0.5) then
! interpolating hhq given on T-grid(lu) to hhv given on v-grid(lcv).
          slu=dble(lu(m,n)+lu(m,n+1))
          hv(m,n)=( hq(m,n  )*dx(m,n  )*dy(m,n  )*dble(lu(m,n  ))       &
                  + hq(m,n+1)*dx(m,n+1)*dy(m,n+1)*dble(lu(m,n+1)) )/slu/dxh(m,n)/dyt(m,n)
         hvp(m,n)=( hqp(m,n  )*dx(m,n  )*dy(m,n  )*dble(lu(m,n  ))       &
                  + hqp(m,n+1)*dx(m,n+1)*dy(m,n+1)*dble(lu(m,n+1)) )/slu/dxh(m,n)/dyt(m,n)
         hvn(m,n)=( hqn(m,n  )*dx(m,n  )*dy(m,n  )*dble(lu(m,n  ))       &
                  + hqn(m,n+1)*dx(m,n+1)*dy(m,n+1)*dble(lu(m,n+1)) )/slu/dxh(m,n)/dyt(m,n)
        endif

        if(luh(m,n)>0.5) then
! interpolating hhq given on T-grid(lu) to hhh given on h-grid(luu).
          slu=dble(lu(m,n)+lu(m+1,n)+lu(m,n+1)+lu(m+1,n+1))
          hh(m,n)=( hq(m  ,n  )*dx(m  ,n  )*dy(m  ,n  )*dble(lu(m  ,n  ))       &
                  + hq(m+1,n  )*dx(m+1,n  )*dy(m+1,n  )*dble(lu(m+1,n  ))       &
                   +hq(m  ,n+1)*dx(m  ,n+1)*dy(m  ,n+1)*dble(lu(m  ,n+1))       &
                  + hq(m+1,n+1)*dx(m+1,n+1)*dy(m+1,n+1)*dble(lu(m+1,n+1)) )/slu/dxb(m,n)/dyb(m,n)
         hhp(m,n)=( hqp(m  ,n  )*dx(m  ,n  )*dy(m  ,n  )*dble(lu(m  ,n  ))       &
                  + hqp(m+1,n  )*dx(m+1,n  )*dy(m+1,n  )*dble(lu(m+1,n  ))       &
                   +hqp(m  ,n+1)*dx(m  ,n+1)*dy(m  ,n+1)*dble(lu(m  ,n+1))       &
                  + hqp(m+1,n+1)*dx(m+1,n+1)*dy(m+1,n+1)*dble(lu(m+1,n+1)) )/slu/dxb(m,n)/dyb(m,n)
         hhn(m,n)=( hqn(m  ,n  )*dx(m  ,n  )*dy(m  ,n  )*dble(lu(m  ,n  ))       &
                  + hqn(m+1,n  )*dx(m+1,n  )*dy(m+1,n  )*dble(lu(m+1,n  ))       &
                   +hqn(m  ,n+1)*dx(m  ,n+1)*dy(m  ,n+1)*dble(lu(m  ,n+1))       &
                  + hqn(m+1,n+1)*dx(m+1,n+1)*dy(m+1,n+1)*dble(lu(m+1,n+1)) )/slu/dxb(m,n)/dyb(m,n)
        endif

       end do
	end do
!$omp end parallel do

      call syncborder_real8(hu)
      call syncborder_real8(hup)
      call syncborder_real8(hun)
      call syncborder_real8(hv)
      call syncborder_real8(hvp)
      call syncborder_real8(hvn)
      call syncborder_real8(hh)
      call syncborder_real8(hhp)
      call syncborder_real8(hhn)

      if(periodicity_x/=0) then
        call cyclize8_x(hu, nx,ny,1,mmm,mm)
        call cyclize8_x(hup,nx,ny,1,mmm,mm)
        call cyclize8_x(hun,nx,ny,1,mmm,mm)
        call cyclize8_x(hv, nx,ny,1,mmm,mm)
        call cyclize8_x(hvp,nx,ny,1,mmm,mm)
        call cyclize8_x(hvn,nx,ny,1,mmm,mm)
        call cyclize8_x(hh, nx,ny,1,mmm,mm)
        call cyclize8_x(hhp,nx,ny,1,mmm,mm)
        call cyclize8_x(hhn,nx,ny,1,mmm,mm)
      end if

      if(periodicity_y/=0) then
        call cyclize8_y(hu, nx,ny,1,nnn,nn)
        call cyclize8_y(hup,nx,ny,1,nnn,nn)
        call cyclize8_y(hun,nx,ny,1,nnn,nn)
        call cyclize8_y(hv, nx,ny,1,nnn,nn)
        call cyclize8_y(hvp,nx,ny,1,nnn,nn)
        call cyclize8_y(hvn,nx,ny,1,nnn,nn)
        call cyclize8_y(hh, nx,ny,1,nnn,nn)
        call cyclize8_y(hhp,nx,ny,1,nnn,nn)
        call cyclize8_y(hhn,nx,ny,1,nnn,nn)
      end if
endsubroutine hh_init

!============================================================
subroutine hh_update(hqn, hun, hvn, hhn, sh,h_r)
 use main_basin_pars
 use mpi_parallel_tools
 use basin_grid
 implicit none

 real(8) hqn(bnd_x1:bnd_x2, bnd_y1:bnd_y2), hun(bnd_x1:bnd_x2, bnd_y1:bnd_y2),    &
         hvn(bnd_x1:bnd_x2, bnd_y1:bnd_y2), hhn(bnd_x1:bnd_x2, bnd_y1:bnd_y2),    &
         sh(bnd_x1:bnd_x2, bnd_y1:bnd_y2),  h_r(bnd_x1:bnd_x2, bnd_y1:bnd_y2)

 integer m,n
 real(8) slu

      hqn =h_r + sh

!$omp parallel do private(m,n,slu)
      do n=ny_start-1,ny_end
       do m=nx_start-1,nx_end

        if(llu(m,n)>0.5) then
! interpolating hhq given on T-grid(lu) to hhu given on u-grid(lcu).
          slu=dble(lu(m,n)+lu(m+1,n))
          hun(m,n)=( hqn(m  ,n)*dx(m  ,n)*dy(m  ,n)*dble(lu(m  ,n))   &
                   + hqn(m+1,n)*dx(m+1,n)*dy(m+1,n)*dble(lu(m+1,n)) )/slu/dxt(m,n)/dyh(m,n)
        endif

        if(llv(m,n)>0.5) then
! interpolating hhq given on T-grid(lu) to hhv given on v-grid(lcv).
          slu=dble(lu(m,n)+lu(m,n+1))
          hvn(m,n)=( hqn(m,n  )*dx(m,n  )*dy(m,n  )*dble(lu(m,n  ))       &
                   + hqn(m,n+1)*dx(m,n+1)*dy(m,n+1)*dble(lu(m,n+1)) )/slu/dxh(m,n)/dyt(m,n)
        endif

        if(luh(m,n)>0.5) then
! interpolating hhq given on T-grid(lu) to hhh given on h-grid(luu).
          slu=dble(lu(m,n)+lu(m+1,n)+lu(m,n+1)+lu(m+1,n+1))
          hhn(m,n)=( hqn(m  ,n  )*dx(m  ,n  )*dy(m  ,n  )*dble(lu(m  ,n  ))       &
                   + hqn(m+1,n  )*dx(m+1,n  )*dy(m+1,n  )*dble(lu(m+1,n  ))       &
                    +hqn(m  ,n+1)*dx(m  ,n+1)*dy(m  ,n+1)*dble(lu(m  ,n+1))       &
                   + hqn(m+1,n+1)*dx(m+1,n+1)*dy(m+1,n+1)*dble(lu(m+1,n+1)) )/slu/dxb(m,n)/dyb(m,n)
        endif

       end do
	end do
!$omp end parallel do

      call syncborder_real8(hun)
      call syncborder_real8(hvn)
      call syncborder_real8(hhn)

      if(periodicity_x/=0) then
        call cyclize8_x(hun, nx,ny,1,mmm,mm)
        call cyclize8_x(hvn, nx,ny,1,mmm,mm)
        call cyclize8_x(hhn, nx,ny,1,mmm,mm)
      end if

      if(periodicity_y/=0) then
        call cyclize8_y(hun, nx,ny,1,nnn,nn)
        call cyclize8_y(hvn, nx,ny,1,nnn,nn)
        call cyclize8_y(hhn, nx,ny,1,nnn,nn)
      end if

endsubroutine hh_update


!============================================================
subroutine hh_shift(hq, hqp, hqn,   &
                    hu, hup, hun,   &
                    hv, hvp, hvn,   &
                    hh, hhp, hhn, nstep  )
 use main_basin_pars
 use mpi_parallel_tools
 use basin_grid
 implicit none

 real(8) hq(bnd_x1:bnd_x2, bnd_y1:bnd_y2), hqp(bnd_x1:bnd_x2, bnd_y1:bnd_y2), hqn(bnd_x1:bnd_x2, bnd_y1:bnd_y2),    &
         hu(bnd_x1:bnd_x2, bnd_y1:bnd_y2), hup(bnd_x1:bnd_x2, bnd_y1:bnd_y2), hun(bnd_x1:bnd_x2, bnd_y1:bnd_y2),    &
         hv(bnd_x1:bnd_x2, bnd_y1:bnd_y2), hvp(bnd_x1:bnd_x2, bnd_y1:bnd_y2), hvn(bnd_x1:bnd_x2, bnd_y1:bnd_y2),    &
         hh(bnd_x1:bnd_x2, bnd_y1:bnd_y2), hhp(bnd_x1:bnd_x2, bnd_y1:bnd_y2), hhn(bnd_x1:bnd_x2, bnd_y1:bnd_y2)

 integer m,n, nstep


!$omp parallel do private(m,n)
      do n=ny_start-1,ny_end+1
       do m=nx_start-1,nx_end+1

        if(llu(m,n)>0.5) then
          hup(m,n)= hu(m,n) + time_smooth*(hun(m,n)-2.0d0*hu(m,n)+hup(m,n))/2.0d0/dfloat(nstep)
           hu(m,n)= hun(m,n)
        endif

        if(llv(m,n)>0.5) then
          hvp(m,n)= hv(m,n) + time_smooth*(hvn(m,n)-2.0d0*hv(m,n)+hvp(m,n))/2.0d0/dfloat(nstep)
           hv(m,n)= hvn(m,n)
        endif

        if(lu(m,n)>0.5) then
          hqp(m,n)= hq(m,n) + time_smooth*(hqn(m,n)-2.0d0*hq(m,n)+hqp(m,n))/2.0d0/dfloat(nstep)
           hq(m,n)= hqn(m,n)
        endif

        if(luh(m,n)>0.5) then
          hhp(m,n)= hh(m,n) + time_smooth*(hhn(m,n)-2.0d0*hh(m,n)+hhp(m,n))/2.0d0/dfloat(nstep)
           hh(m,n)= hhn(m,n)
        endif

       end do
	end do
!$omp end parallel do


endsubroutine hh_shift

!======================================================================
subroutine depth_ave(u,utr,mask,kbcl)
use main_basin_pars
use mpi_parallel_tools
use basin_grid
implicit none
!-----------------------------------------------------------------------
! makes baroclinic and barotropic velocity component
! utr= <u>, u' = u - <u>.
 real(8) u(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz)
 real(4) mask(bnd_x1:bnd_x2,bnd_y1:bnd_y2)
 real(8) utr(bnd_x1:bnd_x2,bnd_y1:bnd_y2),ff0
 integer m, n, k
 integer kbcl      ! key to remove barotropic component from the field: 0 - not, 1 - yes

!$omp parallel do private(m,n,k,ff0)
      do n=bnd_y1+1,bnd_y2-1
       do m=bnd_x1+1,bnd_x2-1
        if (mask(m,n)>0.5) then

           ff0 = 0.0d0
           do k=1,nz
             ff0 = ff0 + u(m,n,k) * dz(k)
           enddo
           utr(m,n) = ff0
           u(m,n,1:nz) = u(m,n,1:nz)-dfloat(kbcl)*ff0
        endif
       enddo
      enddo
!$omp end parallel do


endsubroutine depth_ave
