subroutine tracer_tran_diff(ff,     &
                           ffp,     &
                           tau,     &
                            uu,     &
                            vv,     &
                            ww,     &
                            mu,     &
                     factor_mu,     &
                            nu,     &
                     factor_nu,     &
                          slrx,     &
                          slry,     &
                          slzx,     &
                          slzy,     &
                            az,     &
                            ar,     &
                      alpha_gm,     &
                        ff_top,     &
                        ff_bot,     &
                        ig_top,     &
                        ig_bot,     &
                           RHS_f,   &
                           swbal,   &
                           wflux,   &
                           ksw_lbc, &
                           numlbc,  &
                           lqpx,    &
                           lqpy,    &
                           ind_lbc, &
                           flbc,    &
                           src,     &
                           src_factor)
use main_basin_pars
use mpi_parallel_tools
use basin_grid
implicit none
 
 real(8) tau,az,ar, alpha_gm, src_factor, factor_mu, factor_nu

 real(8) ff(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),      &   !Transported tracer
        ffp(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),      &   !Transported tracer (previous step)
         uu(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),      &
         vv(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),      &
         ww(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1)         !Transporting velocities

 real(8)  slrx(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1),  &
          slry(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1),  &
          slzx(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1),  &
          slzy(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1),  &
            mu(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz  ),  &
            nu(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1)
 
 real(8) ff_top(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
         ff_bot(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
          swbal(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
          wflux(bnd_x1:bnd_x2,bnd_y1:bnd_y2)  
 
 integer ig_top(bnd_x1:bnd_x2,bnd_y1:bnd_y2),    &  !type of sea surface boundary condition (1/2)
         ig_bot(bnd_x1:bnd_x2,bnd_y1:bnd_y2)        !type of sea bottom  boundary condition (1/2)

 real(8) RHS_f(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),      &   !RHS of transport and lareral diffusion
           src(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz)           !Internal source (divergence of SW-radiation)
 integer numlbc, ksw_lbc, lqp
 integer lqpx(numlbc), lqpy(numlbc)
 real(4) flbc(numlbc,nz)
 character ind_lbc(numlbc)
real(8) fz_p,fz_m   !fluxes through cell edges

integer m,n,k
real(8) a(nz),b(nz),c(nz),eta(nz),rksi(nz) 
real(8) dfdx,dfdy,dfdz, mu_u, mu_v, slope_x, slope_y,       &
        bp, bp0, rhs1, rhs2, dp, dm, nu_p, nu_m 


real(8), allocatable::  flux_x(:,:,:), flux_y(:,:,:),      &
                       flux_zx(:,:,:),flux_zy(:,:,:),      &
                       coef_zx(:,:,:),coef_zy(:,:,:),      &
                       fn(:,:,:)

real(8) flux_x1d(nz+1), flux_y1d(nz+1)

    allocate(flux_x(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),     &
             flux_y(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),     &
            flux_zx(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1),   &
            flux_zy(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1),   &
            coef_zx(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1), &
            coef_zy(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1), &
                 fn(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz) ) 

   flux_x=0.0d0
   flux_y=0.0d0
   flux_zx=0.0d0
   flux_zy=0.0d0
   coef_zx=0.0d0
   coef_zy=0.0d0
       fn=0.0d0



!$omp parallel do private(m, n, k, dfdx,dfdy,dfdz, slope_x, slope_y, mu_u, mu_v, flux_x1d, flux_y1d) 
   do n=ny_start,ny_end
    do m=nx_start,nx_end    
     
     if(lcu(m,n)>0.5) then
      
      flux_x1d=0.0d0
      
      do k=2,nz
       dfdz=( ffp(m  ,n,k  )-ffp(m  ,n,k-1)      &
            + ffp(m+1,n,k  )-ffp(m+1,n,k-1)  )/2.0d0
       dfdx=( ffp(m+1,n,k  )-ffp(m  ,n,k  )      &
            + ffp(m+1,n,k-1)-ffp(m  ,n,k-1)  )/2.0d0
       mu_u=(mu(m,n,k)+mu(m,n,k-1)+mu(m+1,n,k)+mu(m+1,n,k-1))/4.0d0*factor_mu
       
       slope_x=ar*( slrx(m,n,k) - alpha_gm*(slrx(m,n,k)-slzx(m,n,k)) ) + az*slzx(m,n,k)
       flux_x1d(k)=mu_u*dyh(m,n)/dxt(m,n)* ( hhu(m,n)*dfdx*hzt(k) - slope_x*dfdz)

       slope_x=ar*( slrx(m,n,k) + alpha_gm*(slrx(m,n,k)-slzx(m,n,k)) ) + az*slzx(m,n,k)
       flux_zx(m,n,k)=mu_u*dyh(m,n)/dxt(m,n) *slope_x*dfdx

       slope_x=ar*slrx(m,n,k) + az*slzx(m,n,k)
       coef_zx(m,n,k)=mu_u*dyh(m,n)/dxt(m,n) *slope_x**2/hhu(m,n)/hzt(k)
      enddo
       
        k=1
       flux_x(m,n,k)=flux_x1d(k+1)/hzt(k+1)        &
                    -uu(m,n,k)*hhu(m,n)*dyh(m,n)*(ff(m,n,k)+ff(m+1,n,k))/2.0d0
        
        do k=2,nz-1
       flux_x(m,n,k)=(flux_x1d(k)+flux_x1d(k+1))/dz(k)/2.0d0        &
                    -uu(m,n,k)*hhu(m,n)*dyh(m,n)*(ff(m,n,k)+ff(m+1,n,k))/2.0d0 
        enddo
        
        k=nz
       flux_x(m,n,k)=flux_x1d(k)/hzt(k)           &
                    -uu(m,n,k)*hhu(m,n)*dyh(m,n)*(ff(m,n,k)+ff(m+1,n,k))/2.0d0
     endif

     if(lcv(m,n)>0.5) then

      flux_y1d=0.0d0

      do k=2,nz
       dfdz=( ffp(m,n  ,k  )-ffp(m,n  ,k-1)      &
            + ffp(m,n+1,k  )-ffp(m,n+1,k-1)  )/2.0d0
       dfdy=( ffp(m,n+1,k  )-ffp(m,n  ,k  )      &
            + ffp(m,n+1,k-1)-ffp(m,n  ,k-1)  )/2.0d0
       mu_v=(mu(m,n,k)+mu(m,n,k-1)+mu(m,n+1,k)+mu(m,n+1,k-1))/4.0d0*factor_mu
       
       slope_y=ar*( slry(m,n,k) - alpha_gm*(slry(m,n,k)-slzy(m,n,k)) ) + az*slzy(m,n,k)
       flux_y1d(k)=mu_v*dxh(m,n)/dyt(m,n)* ( hhv(m,n)*dfdy*hzt(k) - slope_y*dfdz)

       slope_y=ar*( slry(m,n,k) + alpha_gm*(slry(m,n,k)-slzy(m,n,k)) ) + az*slzy(m,n,k)
       flux_zy(m,n,k)=mu_v*dxh(m,n)/dyt(m,n) *slope_y*dfdy

       slope_y=ar*slry(m,n,k) + az*slzy(m,n,k)
       coef_zy(m,n,k)=mu_v*dxh(m,n)/dyt(m,n) *slope_y**2 /hhv(m,n)/hzt(k)       
      enddo

        k=1
       flux_y(m,n,k)=flux_y1d(k+1)/hzt(k+1)       &
                    -vv(m,n,k)*hhv(m,n)*dxh(m,n)*(ff(m,n,k)+ff(m,n+1,k))/2.0d0
        
        do k=2,nz-1
       flux_y(m,n,k)=(flux_y1d(k)+flux_y1d(k+1))/dz(k)/2.0d0        &
                    -vv(m,n,k)*hhv(m,n)*dxh(m,n)*(ff(m,n,k)+ff(m,n+1,k))/2.0d0
        enddo
        
        k=nz
       flux_y(m,n,k)=flux_y1d(k)/hzt(k)           &
                    -vv(m,n,k)*hhv(m,n)*dxh(m,n)*(ff(m,n,k)+ff(m,n+1,k))/2.0d0

     endif  

    enddo
   enddo
!$omp end parallel do        

	if(periodicity_x/=0) then
	 call cyclize8_x( flux_x,nx,ny,nz  ,mmm,mm)
	 call cyclize8_x(flux_zx,nx,ny,nz+1,mmm,mm)
	 call cyclize8_x(coef_zx,nx,ny,nz+1,mmm,mm)
      end if

	if(periodicity_y/=0) then
	 call cyclize8_y( flux_y,nx,ny,nz  ,nnn,nn)            
	 call cyclize8_y(flux_zy,nx,ny,nz+1,nnn,nn)
	 call cyclize8_y(coef_zy,nx,ny,nz+1,nnn,nn)
      end if

!$omp parallel do private(m,n,k,bp,bp0,dp,dm,a,b,c,eta,rksi,fz_p,fz_m, nu_p, nu_m)
      do n=ny_start,ny_end
       do m=nx_start,nx_end
        if (lu(m,n)>0.5) then
        
         bp =  hhqn(m,n)*dx(m,n)*dy(m,n) / tau/2.0d0
         bp0=  hhqp(m,n)*dx(m,n)*dy(m,n) / tau/2.0d0

!Upper layer          
          k=1
          
          nu_m=nu(m,n,k  )*factor_nu
          nu_p=nu(m,n,k+1)*factor_nu

          fz_m=0.0d0
          fz_p=-(flux_zx(m,n,k+1)+flux_zx(m-1,n,k+1))/dble(max(lcu(m,n)+lcu(m-1,n),1e-10))        &
               -(flux_zy(m,n,k+1)+flux_zy(m,n-1,k+1))/dble(max(lcv(m,n)+lcv(m,n-1),1e-10))        &
               + ( coef_zx(m  ,n,k+1)*(ffp(m+1,n,k+1)-ffp(m+1,n,k))                                               &
                  +coef_zx(m-1,n,k+1)*(ffp(m-1,n,k+1)-ffp(m-1,n,k)) )/dble(max(lcu(m,n)+lcu(m-1,n),1e-10))/2.0d0  &
               + ( coef_zy(m,n  ,k+1)*(ffp(m,n+1,k+1)-ffp(m,n+1,k))                                               &
                  +coef_zy(m,n-1,k+1)*(ffp(m,n-1,k+1)-ffp(m,n-1,k)) )/dble(max(lcv(m,n)+lcv(m,n-1),1e-10))/2.0d0  &                  
              -ww(m,n,k+1)*dx(m,n)*dy(m,n)*(ff(m,n,k)+ff(m,n,k+1))/2.0d0

          RHS_f(m,n,k)=( flux_x(m,n,k) - flux_x(m-1,n,k)    &
                       + flux_y(m,n,k) - flux_y(m,n-1,k) )  + (fz_p - fz_m)/dz(k)          &
             + wflux(m,n)*ff(m,n,k)*dx(m,n)*dy(m,n)/RefDen*dfloat(full_free_surface)

          dm = nu_m*dx(m,n)*dy(m,n) /hhqn(m,n)/hzt(k  )
          dp =   (coef_zx(m,n,k+1)+coef_zx(m-1,n,k+1))/dble(max(lcu(m,n)+lcu(m-1,n),1e-10))/2.0d0         &
               + (coef_zy(m,n,k+1)+coef_zy(m,n-1,k+1))/dble(max(lcv(m,n)+lcv(m,n-1),1e-10))/2.0d0         &          
               + nu_p*dx(m,n)*dy(m,n) /hhqn(m,n)/hzt(k+1)

! "2." above is due to aproximation.
          c(k) =  - dp/dz(k)
          a(k) =  0.0d0
         if(ig_top(m,n)==1) then
          b(k) =  bp + (dp + dm)/dz(k)
          eta(k) = bp0*ffp(m,n,k)+dm*ff_top(m,n)/dz(k) +RHS_f(m,n,k)                     &
                 + src(m,n,k)*hhq_rest(m,n)*swbal(m,n)*src_factor*dx(m,n)*dy(m,n)  
         elseif(ig_top(m,n)==2) then
          b(k) =  bp + dp/dz(k) 
          eta(k) = bp0*ffp(m,n,k)+ff_top(m,n)/dz(k)*dx(m,n)*dy(m,n) +RHS_f(m,n,k)  &
                 + src(m,n,k)*hhq_rest(m,n)*swbal(m,n)*src_factor*dx(m,n)*dy(m,n)
         endif

! internal points.
         do k=2,nz-1
          
          nu_m=nu_p
          nu_p=nu(m,n,k+1)*factor_nu
          
          fz_m=fz_p
          fz_p=-(flux_zx(m,n,k+1)+flux_zx(m-1,n,k+1))/dble(max(lcu(m,n)+lcu(m-1,n),1e-10))        &
               -(flux_zy(m,n,k+1)+flux_zy(m,n-1,k+1))/dble(max(lcv(m,n)+lcv(m,n-1),1e-10))        &
               + ( coef_zx(m  ,n,k+1)*(ffp(m+1,n,k+1)-ffp(m+1,n,k))                                               &
                  +coef_zx(m-1,n,k+1)*(ffp(m-1,n,k+1)-ffp(m-1,n,k)) )/dble(max(lcu(m,n)+lcu(m-1,n),1e-10))/2.0d0  &
               + ( coef_zy(m,n  ,k+1)*(ffp(m,n+1,k+1)-ffp(m,n+1,k))                                               &
                  +coef_zy(m,n-1,k+1)*(ffp(m,n-1,k+1)-ffp(m,n-1,k)) )/dble(max(lcv(m,n)+lcv(m,n-1),1e-10))/2.0d0  &
              -ww(m,n,k+1)*dx(m,n)*dy(m,n)*(ff(m,n,k)+ff(m,n,k+1))/2.0d0

          RHS_f(m,n,k)=( flux_x(m,n,k) - flux_x(m-1,n,k)    &
                       + flux_y(m,n,k) - flux_y(m,n-1,k) )  + (fz_p - fz_m)/dz(k)          &
             + wflux(m,n)*ff(m,n,k)*dx(m,n)*dy(m,n)/RefDen*dfloat(full_free_surface)
          
          dm = dp
          dp =  (coef_zx(m,n,k+1)+coef_zx(m-1,n,k+1))/dble(max(lcu(m,n)+lcu(m-1,n),1e-10))/2.0d0     &
               +(coef_zy(m,n,k+1)+coef_zy(m,n-1,k+1))/dble(max(lcv(m,n)+lcv(m,n-1),1e-10))/2.0d0     &
               + nu_p*dx(m,n)*dy(m,n) /hhqn(m,n)/hzt(k+1)         
          
! "2." above is due to aproximation.
          c(k) =  - dp/dz(k)
          a(k) =  - dm/dz(k)
          b(k) =  bp + (dp + dm)/dz(k)
          eta(k) = bp0*ffp(m,n,k) +RHS_f(m,n,k)                     &
                 + src(m,n,k)*hhq_rest(m,n)*swbal(m,n)*src_factor*dx(m,n)*dy(m,n)  
         enddo

!Bottom layer
          k=nz

          nu_m=nu_p
          nu_p=nu(m,n,k+1)*factor_nu

          fz_m= fz_p
          fz_p= 0.0d0

          RHS_f(m,n,k)=( flux_x(m,n,k) - flux_x(m-1,n,k)    &
                       + flux_y(m,n,k) - flux_y(m,n-1,k) )  + (fz_p - fz_m)/dz(k)          &
             + wflux(m,n)*ff(m,n,k)*dx(m,n)*dy(m,n)/RefDen*dfloat(full_free_surface)

          dm = dp
          dp = nu_p*dx(m,n)*dy(m,n) /hhqn(m,n)/hzt(k+1)

! "2." above is due to aproximation.
          c(k) =  0.0d0
          a(k) =  - dm/dz(k)
         
         if(ig_bot(m,n)==1) then
          b(k) =  bp + (dp + dm)/dz(k)
          eta(k) = bp0*ffp(m,n,k) + dp*ff_bot(m,n)/dz(k)+RHS_f(m,n,k)                     &
                 + src(m,n,k)*hhq_rest(m,n)*swbal(m,n)*src_factor*dx(m,n)*dy(m,n) 
         elseif(ig_bot(m,n)==2) then
          b(k) =  bp + dm/dz(k)
          eta(k) = bp0*ffp(m,n,k) + ff_bot(m,n)/dz(k)*dx(m,n)*dy(m,n)+RHS_f(m,n,k)                     &
                 + src(m,n,k)*hhq_rest(m,n)*swbal(m,n)*src_factor*dx(m,n)*dy(m,n) 
         endif

         call factor8(nz,a,b,c,eta,rksi,1,nz)
         do k=1,nz
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
         do k=1,nz
          ffp(m,n,k)=ff(m,n,k)+time_smooth*(fn(m,n,k)-2.0d0*ff(m,n,k)+ffp(m,n,k))/2.0d0
          ff(m,n,k)=fn(m,n,k)
         enddo
        endif

       enddo
      enddo
!$omp end parallel do

   deallocate(fn,coef_zy,coef_zx,flux_zy,flux_zx,flux_y,flux_x)
endsubroutine tracer_tran_diff
