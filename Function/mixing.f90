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

      call syncborder_real8(str_t, nlev)
      call syncborder_real8(str_s, nlev)

      if(periodicity_x/=0) then
       call cyclize8_x(str_t,nx,ny,nlev,mmm,mm)
       call cyclize8_x(str_s,nx,ny,nlev,mmm,mm)
	  end if

      if(periodicity_y/=0) then
       call cyclize8_y(str_t,nx,ny,nlev,nnn,nn)
       call cyclize8_y(str_s,nx,ny,nlev,nnn,nn)
	  end if

endsubroutine stress_components

