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
                        flux_f_x,   &
                        flux_f_y,   &
                           swbal,   &
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
          swbal(bnd_x1:bnd_x2,bnd_y1:bnd_y2)  
 
 integer ig_top(bnd_x1:bnd_x2,bnd_y1:bnd_y2),    &  !type of sea surface boundary condition (1/2)
         ig_bot(bnd_x1:bnd_x2,bnd_y1:bnd_y2)        !type of sea bottom  boundary condition (1/2)

 real(8) flux_f_x(bnd_x1:bnd_x2,bnd_y1:bnd_y2,4),      &   !total flux on x-direction(1-total, 2-advection, 3-diffusion, 4-GM)
         flux_f_y(bnd_x1:bnd_x2,bnd_y1:bnd_y2,4),      &   !total flux on y-direction(1-total, 2-advection, 3-diffusion, 4-GM)
             src(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz)           !Internal source (divergence of SW-radiation)

 integer numlbc, ksw_lbc, lqp
 integer lqpx(numlbc), lqpy(numlbc)
 real(4) flbc(numlbc,nz)
 character ind_lbc(numlbc)
real(8) fz_p,fz_m   !fluxes through cell edges

integer m,n,k
real(8) a(nz),b(nz),c(nz),eta(nz),rksi(nz) 
real(8) dfdx,dfdy,dfdz, slope_x, slope_y,       &
        bp, bp0, rhs1, rhs2, dp, dm, nu_p, nu_m 

real(8), allocatable::  flux_x(:,:,:), flux_y(:,:,:),      &
                       flux_zx(:,:,:),flux_zy(:,:,:),      &
                       coef_zx(:,:,:),coef_zy(:,:,:),      &
                       sfgm_zx(:,:,:),sfgm_zy(:,:,:),      &
                       fn(:,:,:)

real(8) flux_1d(nz+1), mu_1d(nz+1)

real(8) flux_adv, flux_diff, flux_gm, vel_gm, rhs

    allocate(flux_x(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),     &
             flux_y(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),     &
            flux_zx(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1),   &
            flux_zy(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1),   &
            coef_zx(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1), &
            coef_zy(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1), &
            sfgm_zx(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1), &
            sfgm_zy(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1), &
                 fn(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz) ) 

   flux_x=0.0d0
   flux_y=0.0d0
   flux_zx=0.0d0
   flux_zy=0.0d0
   coef_zx=0.0d0
   coef_zy=0.0d0
   sfgm_zx=0.0d0
   sfgm_zy=0.0d0
       fn=0.0d0

   flux_f_x=0.0d0
   flux_f_y=0.0d0

!$omp parallel do private(m, n, k, dfdx,dfdy,dfdz, slope_x, slope_y, mu_1d,    &
!$omp                     flux_1d, flux_adv, flux_diff, flux_gm,vel_gm) 
   do n=ny_start,ny_end
    do m=nx_start,nx_end    
     
     if(lcu(m,n)>0.5) then
      
      flux_1d=0.0d0
        mu_1d=0.0d0
      
      do k=2,nz
       
       dfdz=( ffp(m  ,n,k  )-ffp(m  ,n,k-1)      &
            + ffp(m+1,n,k  )-ffp(m+1,n,k-1)  )/2.0d0
       dfdx=( ffp(m+1,n,k  )-ffp(m  ,n,k  )      &
            + ffp(m+1,n,k-1)-ffp(m  ,n,k-1)  )/2.0d0
       
       mu_1d(k)=(mu(m,n,k)+mu(m,n,k-1)+mu(m+1,n,k)+mu(m+1,n,k-1))/4.0d0*factor_mu*dyh(m,n)/dxt(m,n)
       slope_x=ar*slrx(m,n,k) + az*slzx(m,n,k)

       flux_1d(k)     = mu_1d(k) * ( hhu(m,n)*dfdx*hzt(k) - slope_x*dfdz)
       flux_zx(m,n,k) = mu_1d(k) * slope_x*dfdx
       coef_zx(m,n,k) = mu_1d(k) * slope_x**2/hhu(m,n)/hzt(k)
       sfgm_zx(m,n,k) = mu_1d(k) * (slrx(m,n,k)-slzx(m,n,k)) * alpha_gm
      
      enddo
       
        k=1
      
       flux_diff =  flux_1d(k+1)/hzt(k+1)
       flux_adv  = -uu(m,n,k)*hhu(m,n)*dyh(m,n)*(ff(m,n,k)+ff(m+1,n,k))/2.0d0 
        vel_gm   = (sfgm_zx(m,n,k+1)-sfgm_zx(m,n,k))/dz(k)
       flux_gm   = -vel_gm * (ff(m,n,k)+ff(m+1,n,k))/2.0d0

       flux_x(m,n,k)=flux_adv+flux_diff+flux_gm
       
       flux_f_x(m,n,1)=flux_f_x(m,n,1)+dz(k)*(flux_adv+flux_diff+flux_gm)
       flux_f_x(m,n,2)=flux_f_x(m,n,2)+dz(k)*flux_adv
       flux_f_x(m,n,3)=flux_f_x(m,n,3)+dz(k)*flux_diff
       flux_f_x(m,n,4)=flux_f_x(m,n,4)+dz(k)*flux_gm
       
        do k=2,nz-1
       
       flux_diff = (flux_1d(k)+flux_1d(k+1))/dz(k)/2.0d0 
       flux_adv  = -uu(m,n,k)*hhu(m,n)*dyh(m,n)*(ff(m,n,k)+ff(m+1,n,k))/2.0d0 
        vel_gm   = (sfgm_zx(m,n,k+1)-sfgm_zx(m,n,k))/dz(k)
       flux_gm   = -vel_gm * (ff(m,n,k)+ff(m+1,n,k))/2.0d0

       flux_x(m,n,k)=flux_adv+flux_diff+flux_gm
       
       flux_f_x(m,n,1)=flux_f_x(m,n,1)+dz(k)*(flux_adv+flux_diff+flux_gm)
       flux_f_x(m,n,2)=flux_f_x(m,n,2)+dz(k)*flux_adv
       flux_f_x(m,n,3)=flux_f_x(m,n,3)+dz(k)*flux_diff
       flux_f_x(m,n,4)=flux_f_x(m,n,4)+dz(k)*flux_gm
        
        enddo
        
        k=nz
       
       flux_diff = flux_1d(k)/hzt(k) 
       flux_adv  = -uu(m,n,k)*hhu(m,n)*dyh(m,n)*(ff(m,n,k)+ff(m+1,n,k))/2.0d0 
        vel_gm   = (sfgm_zx(m,n,k+1)-sfgm_zx(m,n,k))/dz(k)
       flux_gm   = -vel_gm * (ff(m,n,k)+ff(m+1,n,k))/2.0d0

       flux_x(m,n,k)=flux_adv+flux_diff+flux_gm
       
       flux_f_x(m,n,1)=flux_f_x(m,n,1)+dz(k)*(flux_adv+flux_diff+flux_gm)
       flux_f_x(m,n,2)=flux_f_x(m,n,2)+dz(k)*flux_adv
       flux_f_x(m,n,3)=flux_f_x(m,n,3)+dz(k)*flux_diff
       flux_f_x(m,n,4)=flux_f_x(m,n,4)+dz(k)*flux_gm
     
     endif

     if(lcv(m,n)>0.5) then

      flux_1d=0.0d0
        mu_1d=0.0d0

      do k=2,nz
       
       dfdz=( ffp(m,n  ,k  )-ffp(m,n  ,k-1)      &
            + ffp(m,n+1,k  )-ffp(m,n+1,k-1)  )/2.0d0
       dfdy=( ffp(m,n+1,k  )-ffp(m,n  ,k  )      &
            + ffp(m,n+1,k-1)-ffp(m,n  ,k-1)  )/2.0d0
       mu_1d(k)=(mu(m,n,k)+mu(m,n,k-1)+mu(m,n+1,k)+mu(m,n+1,k-1))/4.0d0*factor_mu*dxh(m,n)/dyt(m,n)
       
       slope_y=ar*slry(m,n,k) + az*slzy(m,n,k)
              
       flux_1d(k)     = mu_1d(k) * ( hhv(m,n)*dfdy*hzt(k) - slope_y*dfdz)
       flux_zy(m,n,k) = mu_1d(k) * slope_y*dfdy
       coef_zy(m,n,k) = mu_1d(k) * slope_y**2 /hhv(m,n)/hzt(k)
       sfgm_zy(m,n,k) = mu_1d(k) * (slry(m,n,k)-slzy(m,n,k)) * alpha_gm              
      
      enddo

        k=1
       
       flux_diff=flux_1d(k+1)/hzt(k+1)
       flux_adv  = -vv(m,n,k)*hhv(m,n)*dxh(m,n)*(ff(m,n,k)+ff(m,n+1,k))/2.0d0 
        vel_gm   = (sfgm_zy(m,n,k+1)-sfgm_zy(m,n,k))/dz(k)
       flux_gm   = -vel_gm * (ff(m,n,k)+ff(m,n+1,k))/2.0d0

       flux_y(m,n,k)=flux_adv+flux_diff+flux_gm
       
       flux_f_y(m,n,1)=flux_f_y(m,n,1)+dz(k)*(flux_adv+flux_diff+flux_gm)
       flux_f_y(m,n,2)=flux_f_y(m,n,2)+dz(k)*flux_adv
       flux_f_y(m,n,3)=flux_f_y(m,n,3)+dz(k)*flux_diff
       flux_f_y(m,n,4)=flux_f_y(m,n,4)+dz(k)*flux_gm
        
        do k=2,nz-1
       
       flux_diff=(flux_1d(k)+flux_1d(k+1))/dz(k)/2.0d0  
       flux_adv  = -vv(m,n,k)*hhv(m,n)*dxh(m,n)*(ff(m,n,k)+ff(m,n+1,k))/2.0d0 
        vel_gm   = (sfgm_zy(m,n,k+1)-sfgm_zy(m,n,k))/dz(k)
       flux_gm   = -vel_gm * (ff(m,n,k)+ff(m,n+1,k))/2.0d0

       flux_y(m,n,k)=flux_adv+flux_diff+flux_gm
       
       flux_f_y(m,n,1)=flux_f_y(m,n,1)+dz(k)*(flux_adv+flux_diff+flux_gm)
       flux_f_y(m,n,2)=flux_f_y(m,n,2)+dz(k)*flux_adv
       flux_f_y(m,n,3)=flux_f_y(m,n,3)+dz(k)*flux_diff
       flux_f_y(m,n,4)=flux_f_y(m,n,4)+dz(k)*flux_gm
       
        enddo
        
        k=nz
       
       flux_diff=flux_1d(k)/hzt(k) 
       flux_adv  = -vv(m,n,k)*hhv(m,n)*dxh(m,n)*(ff(m,n,k)+ff(m,n+1,k))/2.0d0 
        vel_gm   = (sfgm_zy(m,n,k+1)-sfgm_zy(m,n,k))/dz(k)
       flux_gm   = -vel_gm * (ff(m,n,k)+ff(m,n+1,k))/2.0d0

       flux_y(m,n,k)=flux_adv+flux_diff+flux_gm
       
       flux_f_y(m,n,1)=flux_f_y(m,n,1)+dz(k)*(flux_adv+flux_diff+flux_gm)
       flux_f_y(m,n,2)=flux_f_y(m,n,2)+dz(k)*flux_adv
       flux_f_y(m,n,3)=flux_f_y(m,n,3)+dz(k)*flux_diff
       flux_f_y(m,n,4)=flux_f_y(m,n,4)+dz(k)*flux_gm

     endif  

    enddo
   enddo
!$omp end parallel do        

	if(periodicity_x/=0) then
	 call cyclize8_x( flux_x,nx,ny,nz  ,mmm,mm)
	 call cyclize8_x(flux_zx,nx,ny,nz+1,mmm,mm)
	 call cyclize8_x(coef_zx,nx,ny,nz+1,mmm,mm)
	 call cyclize8_x(sfgm_zx,nx,ny,nz+1,mmm,mm)
      end if

	if(periodicity_y/=0) then
	 call cyclize8_y( flux_y,nx,ny,nz  ,nnn,nn)            
	 call cyclize8_y(flux_zy,nx,ny,nz+1,nnn,nn)
	 call cyclize8_y(coef_zy,nx,ny,nz+1,nnn,nn)
	 call cyclize8_y(sfgm_zy,nx,ny,nz+1,nnn,nn)
      end if

!$omp parallel do private(m,n,k,bp,bp0,dp,dm,a,b,c,eta,rksi,fz_p,fz_m, nu_p, nu_m, vel_gm,rhs)
      do n=ny_start,ny_end
       do m=nx_start,nx_end
        if (lu(m,n)>0.5) then
        
         bp =  hhqn(m,n)*dx(m,n)*dy(m,n) / tau/2.0d0
         bp0=  hhqp(m,n)*dx(m,n)*dy(m,n) / tau/2.0d0

!Upper layer          
          k=1
          
          nu_m=nu(m,n,k  )*factor_nu
          nu_p=nu(m,n,k+1)*factor_nu
          
          vel_gm=-(sfgm_zx(m,n,k+1)-sfgm_zx(m-1,n,k+1)+sfgm_zy(m,n,k+1)-sfgm_zy(m,n-1,k+1))
          
          fz_m=0.0d0
          fz_p=-(flux_zx(m,n,k+1)+flux_zx(m-1,n,k+1))/dble(max(lcu(m,n)+lcu(m-1,n),1e-10))        &
               -(flux_zy(m,n,k+1)+flux_zy(m,n-1,k+1))/dble(max(lcv(m,n)+lcv(m,n-1),1e-10))        &
               + ( coef_zx(m  ,n,k+1)*(ffp(m+1,n,k+1)-ffp(m+1,n,k))                                               &
                  +coef_zx(m-1,n,k+1)*(ffp(m-1,n,k+1)-ffp(m-1,n,k)) )/dble(max(lcu(m,n)+lcu(m-1,n),1e-10))/2.0d0  &
               + ( coef_zy(m,n  ,k+1)*(ffp(m,n+1,k+1)-ffp(m,n+1,k))                                               &
                  +coef_zy(m,n-1,k+1)*(ffp(m,n-1,k+1)-ffp(m,n-1,k)) )/dble(max(lcv(m,n)+lcv(m,n-1),1e-10))/2.0d0  &                  
              -(ww(m,n,k+1)*dx(m,n)*dy(m,n)+vel_gm)*(ff(m,n,k)+ff(m,n,k+1))/2.0d0

          dm = nu_m*dx(m,n)*dy(m,n) /hhqn(m,n)/hzt(k  )
          dp =   (coef_zx(m,n,k+1)+coef_zx(m-1,n,k+1))/dble(max(lcu(m,n)+lcu(m-1,n),1e-10))/2.0d0         &
               + (coef_zy(m,n,k+1)+coef_zy(m,n-1,k+1))/dble(max(lcv(m,n)+lcv(m,n-1),1e-10))/2.0d0         &          
               + nu_p*dx(m,n)*dy(m,n) /hhqn(m,n)/hzt(k+1)

          rhs =  flux_x(m,n,k) - flux_x(m-1,n,k) + flux_y(m,n,k) - flux_y(m,n-1,k) + (fz_p - fz_m)/dz(k)    &           
                              + src(m,n,k)*hhq_rest(m,n)*swbal(m,n)*src_factor*dx(m,n)*dy(m,n) 

! "2." above is due to aproximation.
          c(k) =  - dp/dz(k)
          a(k) =  0.0d0
         if(ig_top(m,n)==1) then
          
          b(k) =  bp + (dp + dm)/dz(k)
          eta(k) = bp0*ffp(m,n,k) + rhs + dm*ff_top(m,n)/dz(k) 

         elseif(ig_top(m,n)==2) then
          
          b(k) =  bp + dp/dz(k) 
          eta(k) = bp0*ffp(m,n,k) + rhs + ff_top(m,n)/dz(k)*dx(m,n)*dy(m,n)
         endif

! internal points.
         do k=2,nz-1
          
          nu_p=nu(m,n,k+1)*factor_nu
          
          vel_gm=-(sfgm_zx(m,n,k+1)-sfgm_zx(m-1,n,k+1)+sfgm_zy(m,n,k+1)-sfgm_zy(m,n-1,k+1))

          fz_m=fz_p
          fz_p=-(flux_zx(m,n,k+1)+flux_zx(m-1,n,k+1))/dble(max(lcu(m,n)+lcu(m-1,n),1e-10))        &
               -(flux_zy(m,n,k+1)+flux_zy(m,n-1,k+1))/dble(max(lcv(m,n)+lcv(m,n-1),1e-10))        &
               + ( coef_zx(m  ,n,k+1)*(ffp(m+1,n,k+1)-ffp(m+1,n,k))                                               &
                  +coef_zx(m-1,n,k+1)*(ffp(m-1,n,k+1)-ffp(m-1,n,k)) )/dble(max(lcu(m,n)+lcu(m-1,n),1e-10))/2.0d0  &
               + ( coef_zy(m,n  ,k+1)*(ffp(m,n+1,k+1)-ffp(m,n+1,k))                                               &
                  +coef_zy(m,n-1,k+1)*(ffp(m,n-1,k+1)-ffp(m,n-1,k)) )/dble(max(lcv(m,n)+lcv(m,n-1),1e-10))/2.0d0  &
              -(ww(m,n,k+1)*dx(m,n)*dy(m,n)+vel_gm)*(ff(m,n,k)+ff(m,n,k+1))/2.0d0

          dm = dp
          dp =  (coef_zx(m,n,k+1)+coef_zx(m-1,n,k+1))/dble(max(lcu(m,n)+lcu(m-1,n),1e-10))/2.0d0     &
               +(coef_zy(m,n,k+1)+coef_zy(m,n-1,k+1))/dble(max(lcv(m,n)+lcv(m,n-1),1e-10))/2.0d0     &
               + nu_p*dx(m,n)*dy(m,n) /hhqn(m,n)/hzt(k+1)         
          
          rhs =  flux_x(m,n,k) - flux_x(m-1,n,k) + flux_y(m,n,k) - flux_y(m,n-1,k) + (fz_p - fz_m)/dz(k)    &     
                              + src(m,n,k)*hhq_rest(m,n)*swbal(m,n)*src_factor*dx(m,n)*dy(m,n) 

! "2." above is due to aproximation.
          c(k) =  - dp/dz(k)
          a(k) =  - dm/dz(k)
          b(k) =  bp + (dp + dm)/dz(k)

          eta(k) = bp0*ffp(m,n,k) + rhs
         enddo

!Bottom layer
          k=nz

          nu_p=nu(m,n,k+1)*factor_nu

          fz_m= fz_p
          fz_p= 0.0d0

          dm = dp
          dp = nu_p*dx(m,n)*dy(m,n) /hhqn(m,n)/hzt(k+1)
          
          rhs =  flux_x(m,n,k) - flux_x(m-1,n,k) + flux_y(m,n,k) - flux_y(m,n-1,k) + (fz_p - fz_m)/dz(k)    &           
                              + src(m,n,k)*hhq_rest(m,n)*swbal(m,n)*src_factor*dx(m,n)*dy(m,n) 

! "2." above is due to aproximation.
          c(k) =  0.0d0
          a(k) =  - dm/dz(k)
         
         if(ig_bot(m,n)==1) then
          
            b(k) =  bp + (dp + dm)/dz(k)
          eta(k) = bp0*ffp(m,n,k) + rhs + dp*ff_bot(m,n)/dz(k) 
         
         elseif(ig_bot(m,n)==2) then

          b(k) =  bp + dm/dz(k)
          eta(k) = bp0*ffp(m,n,k) + rhs + ff_bot(m,n)/dz(k)*dx(m,n)*dy(m,n)
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
!          ffp(m,n,k)=hhq(m,n)*ff(m,n,k)+time_smooth*(hhqn(m,n)*fn(m,n,k)-2.0d0*hhq(m,n)*ff(m,n,k)+hhqp(m,n)*ffp(m,n,k))/2.0d0
          ffp(m,n,k)=ff(m,n,k)+time_smooth*(fn(m,n,k)-2.0d0*ff(m,n,k)+ffp(m,n,k))/2.0d0
           ff(m,n,k)=fn(m,n,k)
         enddo
        endif

       enddo
      enddo
!$omp end parallel do

   deallocate(fn,sfgm_zy,sfgm_zx,coef_zy,coef_zx,flux_zy,flux_zx,flux_y,flux_x)
endsubroutine tracer_tran_diff
