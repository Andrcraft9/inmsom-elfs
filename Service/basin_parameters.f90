!======================================================================
subroutine vgrid
use main_basin_pars
use basin_grid
use mpi_parallel_tools
implicit none
  ! for setting vertical t-,w- grid levels
  ! zw(1) is surface, zw(nz+1) is bottom w-levels.
  ! version with calibration on levitus levels.
  ! hh0 is reference depth of ocean.

include 'vgrid.fi'
real(8), parameter:: hh0=3500.0d0 !mean depth of the World Ocean
integer, parameter:: nlev=33
real(8) dlev(nlev)         !levitus horizonts in meters for analitical set
data dlev/0.0d0,10.0d0,20.0d0,30.0d0,50.0d0,75.0d0,100.0d0,125.0d0,150.0d0,200.0d0,250.0d0,   &
  300.0d0,400.0d0,500.0d0,600.0d0,700.0d0,800.0d0,900.0d0,1000.0d0,1100.0d0,1200.0d0,         &
 1300.0d0,1400.0d0,1500.0d0,1750.0d0,2000.0d0,2500.0d0,3000.0d0,3500.0d0,4000.0d0,            &
 4500.0d0,5000.0d0,5500.0d0/
integer k, nrefl
real(8) bottom, devih, a, b, unidepth

  ! finding levitus number of reference depth
   nrefl=1
   devih=dlev(33)

 ! analytical set of vertical grid
    do k=2,33
       if (abs(dlev(k)-hh0) <= devih) then
            devih=abs(dlev(k)-hh0)
            nrefl=k
       endif
    enddo
    if (rank .eq. 0) then
        write(*,'(a,i3,a,f8.2)')  ' number of levitus horizonts:',nrefl,' for h=',hh0
    endif
    a=10.0d0*dfloat(nrefl-1)/dlev(nrefl)   ! gradient in upper ocean
    b=exp(1.0)
 ! may use b1 for more slightly levels in upper ocean
 !     b1=sqrt((float(nrefl)/float(nz)))

    zw(1)=0.0d0                  !sea surface
    zw(nz+1)=1.0d0               !bottom

    if (wgr_in_tgr) then
 ! w-levels are arranged in the middles of t-layers
      if(analytical_set) then
 ! analitical t-levels setting
       do k=2,nz-1
        z(k)=unidepth((dfloat(k)-0.5d0)/(dfloat(nz)-0.5d0),a,b)
       enddo
      else
 ! non-analitical t-levels setting
 ! prove the levels
      do k = 2,nz
       if(z_manual(k) <= z_manual(k-1)) then
        write(*,'(a,i4,f10.5)')  '  error in seting z-levels in 0vgrid.fi. horizont #',k,z_manual(k)
        stop 1
       end if
      enddo

 ! correct the levels for 1-depth
      bottom=z_manual(nz)+(z_manual(nz)-z_manual(nz-1))/2.0d0
      do k=1,nz
       z(k)=z_manual(k)/bottom
      enddo

    end if
 ! regulating top and bottom t-levels
    z( 1) =        z(2)   /3.0d0
    z(nz) =2.0d0/3.0d0 + z(nz-1)/3.0d0
 ! w-levels setting in the middles of t-layers
    do k =2,nz
      zw(k  )= (z (k)  + z (k-1))/2.0d0
    enddo

    else

 ! t-levels are arranged in the middles of w-layers
    if(analytical_set) then
 ! analitical w-level setting
     do k=3,nz
      zw(k)=unidepth((dfloat(k-1))/(dfloat(nz)-0.5d0),a,b)
     enddo

    else

 ! non-analytical w-levels setting
 ! prove the levels
    do k = 2,nz
     if(z_manual(k) <= z_manual(k-1)) then
     write(*,'(a,i4,f10.5)') '  error in setting z-levels in 1vgrid.fi. horizon �',k,z_manual(k)
     stop 1
     end if
    enddo

 ! correct the levels for 1-depth
          bottom=z_manual(nz)+(z_manual(nz)-z_manual(nz-1))
          do k=2,nz
           zw(k)=z_manual(k)/bottom
          enddo

         endif
 ! regulating top and bottom levels
           zw( 2) =  zw(3)/2.0d0
           zw(nz) = (zw(nz-1)+zw(nz+1))/2.0d0

 ! t-level setting in the middle of w-layer
          do k =1,nz
           z(k )  = (zw(k+1) + zw(k))/2.0d0
          enddo
       endif

 ! t and w -grid steps:
        hzt(1) = z (1)
         dz(1) = zw(2)
       do  k=2,nz
        hzt(k) = z (k)   - z (k-1)
         dz(k) = zw(k+1) - zw(k)
       enddo
        hzt(nz+1) = 1.0d0-z(nz)

 !     bottom=hh0
       bottom=1000.0d0

       if (rank .eq. 0) then
           if (wgr_in_tgr) then
               write(*,*)'  w-levels are arranged in the middles of t-layers.'
           else
               write(*,*)'  t-levels are arranged in the middles of w-layers.'
           end if

           write(*,110) bottom
   110     format('  w-levels w-steps  t-levels t-steps *',f7.2)
           do k=1,nz
               write(*,111)zw(k+1)*bottom,dz(k)*bottom,z(k)*bottom,hzt(k)*bottom
           enddo
   111     format(2(2x,2f8.2))
       endif

endsubroutine vgrid
!============================================================================================
      function unidepth(x,a,b)
        implicit none
      !  universal dimensionless function of non-uniform oceanographic horizons
      !  constructed on levitus oceanographic horizons
      !  x-dimensionless level value from [0,1]
        real(8) unidepth, x, a, b
        unidepth=(2.0d0-a)**(x**b)+a*x-1.0d0
      endfunction unidepth
!===========================================================================================
!initializing basin grid parameters
 subroutine basinpar
 use main_basin_pars
 use mpi_parallel_tools
 use basin_grid
 use rec_length
 implicit none
  integer m,n,ierr
  real(8), parameter:: dpip180 = 3.1415926535897/180.0d0   !for degrees to radians convers
  real(4) array4(bnd_x1:bnd_x2,bnd_y1:bnd_y2)

! temperature grid initialization

! x-coordinate (in degrees)
! in case of regular grid
 if(xgr_type==0) then
  do m=bnd_x1,bnd_x2
   xt(m)=rlon+dfloat(m-mmm)*dxst
  end do
 else !in case of irregular grid
  do m=bnd_x1,bnd_x2
   xt(m)=x_levels(m)
  end do
 endif

! y-coordinate (in degrees)
! in case of regular grid

 if(ygr_type==0) then
  do n=bnd_y1,bnd_y2
   yt(n)=rlat+dfloat(n-nnn)*dyst
  end do
 else !in case of irregular grid
  do n=bnd_y1,bnd_y2
   yt(n)=y_levels(n)
  end do
 endif

! parameters:
 if (rank .eq. 0) then
     write(*,'(2x,a)')' Basin parameters from 1basinpar.inc:'

     if(curve_grid==0) then        !Carthesian coordinates
         write(*,*) 'Coordinate system is carthesian'
     elseif(curve_grid==1) then
         write(*,*) 'Coordinate system is undistorted sphere'
         write(*,'(a,f10.3)') ' rotation angle on longitude is =',rotation_on_lon,    &
                           ' rotation angle on  latitude is =',rotation_on_lat
     elseif(curve_grid==2) then
         write(*,*) 'Coordinate system is distorted sphere'
         write(*,'(a,f10.3)') ' geo longitude of new north pole is =',x_pole,    &
                           ' geo  latitude of new north pole is =',y_pole,    &
                           ' geo longitude of new south pole is =',p_pole,    &
                           ' geo  latitude of new south pole is =',q_pole
     endif

     if(xgr_type==0) then
         write(*,*) 'X-grid is uniform'
         write(*,'(2(a,f10.3),a)') ' initial x-coordinate (m=mmm) =',rlon,' step on x =',dxst,'[dgr] '
     else
         write(*,*) 'X-grid is non-uniform'
         write(*,'(a,f10.3)') ' minimal x-coordinate (m=mmm) =',xt(mmm),    &
                               ' maximal x-coordinate (m=mm ) =',xt(mm)
     endif

     if(ygr_type==0) then
         write(*,*) 'Y-grid is uniform'
         write(*,'(2(a,f10.3),a)') ' initial y-coordinate (n=nnn) =',rlat,' step on y =',dyst,'[dgr] '
     else
         write(*,*) 'Y-grid is non-uniform'
         write(*,'(a,f10.3)') ' minimal y-coordinate (n=nnn) =',yt(nnn),    &
                               ' maximal y-coordinate (n=nn ) =',yt(nn)
     endif

     write(*,'(2(a,i2))') 'Periodicity on X =', periodicity_x,', Periodicity on Y =', periodicity_y
     write(*,'(4(a,i4))') '  nx=',nx, ';  ny=',ny,';  nz=',nz
     write(*,'(4(a,i4))') ' mmm=',mmm,';  mm=',mm,'; nnn=',nnn,';  nn=',nn

     write(*,'(2x,a,g14.7,a)')' Earth radius =',RadEarth,'(m)'
     write(*,'(2x,a,g14.7,a)') 'Earth angular velocity(omega) =',EarthAngVel,'[rad/sec]'
     write(*,'(2x,a,g14.7,a)') 'Heat capacity of water =',HeatCapWater,'[J/kg/�C] for 35%. sal'
     write(*,'(2x,a,g14.7,a)') 'reference density =',RefDen,'[kg/m**3]'
     write(*,'(2x,a,f10.3,a)') 'free fall acceleration(grv)=',FreeFallAcc,'[m/s**2]'
 endif
! velocity grid initialization

! x-coordinate (in degrees)
      do m=bnd_x1,bnd_x2-1
	  xu(m)=(xt(m)+xt(m+1))/2.0d0
      end do

! y-coordinate (in degrees)
      do n=bnd_y1,bnd_y2-1
	  yv(n)=(yt(n)+yt(n+1))/2.0d0
      end do

!Initialization of x-steps
if(xgr_type>0) then
!      do n=ny_start-1,ny_end+1
!       do m=nx_start-1,nx_end+1
      do n=bnd_y1+1, bnd_y2-1
       do m=bnd_x1-1, bnd_x2-1
!-----initialization of t- and v-grid x-steps in centimeters
	    dxt(m,n)=(xt(m+1)-xt(m))*dpip180*RadEarth
	    dxb(m,n)=(xt(m+1)-xt(m))*dpip180*RadEarth
!-----initialization of u- and h-grid x-steps in centimeters
	    dx(m,n)=(xu(m)-xu(m-1))*dpip180*RadEarth
	   dxh(m,n)=(xu(m)-xu(m-1))*dpip180*RadEarth
        end do
       end do
else
!      do n=ny_start-1,ny_end+1
!       do m=nx_start-1,nx_end+1
      do n=bnd_y1+1, bnd_y2-1
       do m=bnd_x1-1, bnd_x2-1
!-----initialization of t- and v-grid x-steps in centimeters
	    dxt(m,n)=dxst*dpip180*RadEarth
	    dxb(m,n)=dxst*dpip180*RadEarth
!-----initialization of u- and h-grid x-steps in centimeters
	    dx(m,n)=dxst*dpip180*RadEarth
	   dxh(m,n)=dxst*dpip180*RadEarth
        end do
       end do
endif

!Initialization of y-steps
if(ygr_type>0) then
!      do n=ny_start-1,ny_end+1
!       do m=nx_start-1,nx_end+1
      do n=bnd_y1+1, bnd_y2-1
       do m=bnd_x1-1, bnd_x2-1
!-----initialization of t- and u-grid y-steps in centimeters
	    dyt(m,n)=(yt(n+1)-yt(n))*dpip180*RadEarth
	    dyb(m,n)=(yt(n+1)-yt(n))*dpip180*RadEarth
!-----initialization of v- and h-grid y-steps in centimeters
	    dy(m,n)=(yv(n)-yv(n-1))*dpip180*RadEarth
	   dyh(m,n)=(yv(n)-yv(n-1))*dpip180*RadEarth
        end do
       end do
else
!      do n=ny_start-1,ny_end+1
!       do m=nx_start-1,nx_end+1
      do n=bnd_y1+1, bnd_y2-1
       do m=bnd_x1-1, bnd_x2-1
!-----initialization of t- and u-grid y-steps in centimeters
	    dyt(m,n)=dyst*dpip180*RadEarth
	    dyb(m,n)=dyst*dpip180*RadEarth
!-----initialization of v- and h-grid y-steps in centimeters
	    dy(m,n)=dyst*dpip180*RadEarth
	   dyh(m,n)=dyst*dpip180*RadEarth
        end do
       end do
endif

!-----initialization of Coriolis terms--------------------------
       rlh_s= 2.0d0*EarthAngVel
       rlh_c=-2.0d0*EarthAngVel

!-----metric initialization--------------------------------------------------------------
      if(curve_grid==0) then   !in case of carthesian grid

!On T-grid
       call grid_parameters_carthesian(xt,   &   !model x-coordinate in degrees
                                       yt,   &   !model y-coordinate in degrees
                                   bnd_x1,   &   !left   boundary of arrays
                                   bnd_x2,   &   !right  boundary of arrays
                                   bnd_y1,   &   !lower  boundary of arrays
                                   bnd_y2,   &   !upper  boundary of arrays
                                geo_lon_t,   &   !geographical longitude in degrees
                                geo_lat_t,   &   !geographical latitude  in degrees
                                      dx,    &   !metrical coefficient on x
                                      dy,    &   !metrical coefficient on x
                              rotvec_coeff,  &   !rotation coefficients for vector transform
                                    rlh_s,   &   !coriolis main term (sin)
                                    rlh_c,   &   !coriolis second term (cos)
                                        1,   &   !key to compute rotation coefficients (0/1)
                                        0,   &   !key to compute coriolis terms (0/1)
                                 nx_start-1, &   !first significant point in x-direction (output)
                                 nx_end+1,   &   ! last significant point in x-direction (output)
                                 ny_start-1, &   !first significant point in y-direction (output)
                                 ny_end+1    )   ! last significant point in y-direction (output)

!On U-grid
       call grid_parameters_carthesian(xu,   &   !model x-coordinate in degrees
                                       yt,   &   !model y-coordinate in degrees
                                   bnd_x1,   &   !left   boundary of arrays
                                   bnd_x2,   &   !right  boundary of arrays
                                   bnd_y1,   &   !lower  boundary of arrays
                                   bnd_y2,   &   !upper  boundary of arrays
                                geo_lon_u,   &   !geographical longitude in degrees
                                geo_lat_u,   &   !geographical latitude  in degrees
                                     dxt,    &   !metrical coefficient on x
                                     dyh,    &   !metrical coefficient on x
                              rotvec_coeff,  &   !rotation coefficients for vector transform
                                    rlh_s,   &   !coriolis main term (sin)
                                    rlh_c,   &   !coriolis second term (cos)
                                        0,   &   !key to compute rotation coefficients (0/1)
                                        0,   &   !key to compute coriolis terms (0/1)
                                 nx_start-1, &   !first significant point in x-direction (output)
                                 nx_end+1,   &   ! last significant point in x-direction (output)
                                 ny_start-1, &   !first significant point in y-direction (output)
                                 ny_end+1    )   ! last significant point in y-direction (output)

!On V-grid
       call grid_parameters_carthesian(xt,   &   !model x-coordinate in degrees
                                       yv,   &   !model y-coordinate in degrees
                                   bnd_x1,   &   !left   boundary of arrays
                                   bnd_x2,   &   !right  boundary of arrays
                                   bnd_y1,   &   !lower  boundary of arrays
                                   bnd_y2,   &   !upper  boundary of arrays
                                geo_lon_v,   &   !geographical longitude in degrees
                                geo_lat_v,   &   !geographical latitude  in degrees
                                     dxh,    &   !metrical coefficient on x
                                     dyt,    &   !metrical coefficient on x
                              rotvec_coeff,  &   !rotation coefficients for vector transform
                                    rlh_s,   &   !coriolis main term (sin)
                                    rlh_c,   &   !coriolis second term (cos)
                                        0,   &   !key to compute rotation coefficients (0/1)
                                        0,   &   !key to compute coriolis terms (0/1)
                                 nx_start-1, &   !first significant point in x-direction (output)
                                 nx_end+1,   &   ! last significant point in x-direction (output)
                                 ny_start-1, &   !first significant point in y-direction (output)
                                 ny_end+1    )   ! last significant point in y-direction (output)

!On H-grid
       call grid_parameters_carthesian(xu,   &   !model x-coordinate in degrees
                                       yv,   &   !model y-coordinate in degrees
                                   bnd_x1,   &   !left   boundary of arrays
                                   bnd_x2,   &   !right  boundary of arrays
                                   bnd_y1,   &   !lower  boundary of arrays
                                   bnd_y2,   &   !upper  boundary of arrays
                                geo_lon_h,   &   !geographical longitude in degrees
                                geo_lat_h,   &   !geographical latitude  in degrees
                                     dxb,    &   !metrical coefficient on x
                                     dyb,    &   !metrical coefficient on x
                              rotvec_coeff,  &   !rotation coefficients for vector transform
                                    rlh_s,   &   !coriolis main term (sin)
                                    rlh_c,   &   !coriolis second term (cos)
                                        0,   &   !key to compute rotation coefficients (0/1)
                                        1,   &   !key to compute coriolis terms (0/1)
                                 nx_start-1, &   !first significant point in x-direction (output)
                                 nx_end+1,   &   ! last significant point in x-direction (output)
                                 ny_start-1, &   !first significant point in y-direction (output)
                                 ny_end+1    )   ! last significant point in y-direction (output)

      elseif(curve_grid==1) then !in case of spherical grid

!On T-grid
       call grid_parameters_spherical (xt,   &   !model x-coordinate in degrees
                                       yt,   &   !model y-coordinate in degrees
                                   bnd_x1,   &   !left   boundary of arrays
                                   bnd_x2,   &   !right  boundary of arrays
                                   bnd_y1,   &   !lower  boundary of arrays
                                   bnd_y2,   &   !upper  boundary of arrays
                          rotation_on_lon,   &   !euler angle on longitude
                          rotation_on_lat,   &   !euler angle on latitude
                                geo_lon_t,   &   !geographical longitude in degrees
                                geo_lat_t,   &   !geographical latitude  in degrees
                                      dx,    &   !metrical coefficient on x
                                      dy,    &   !metrical coefficient on x
                              rotvec_coeff,  &   !rotation coefficients for vector transform
                                    rlh_s,   &   !coriolis main term (sin)
                                    rlh_c,   &   !coriolis second term (cos)
                                        1,   &   !key to compute rotation coefficients (0/1)
                                        0,   &   !key to compute coriolis terms (0/1)
                                 bnd_x1+1,   &   !first significant point in x-direction (output)
                                 bnd_x2-1,   &   ! last significant point in x-direction (output)
                                 bnd_y1+1,   &   !first significant point in y-direction (output)
                                 bnd_y2-1    )   ! last significant point in y-direction (output)

!On U-grid
       call grid_parameters_spherical (xu,   &   !model x-coordinate in degrees
                                       yt,   &   !model y-coordinate in degrees
                                   bnd_x1,   &   !left   boundary of arrays
                                   bnd_x2,   &   !right  boundary of arrays
                                   bnd_y1,   &   !lower  boundary of arrays
                                   bnd_y2,   &   !upper  boundary of arrays
                          rotation_on_lon,   &   !euler angle on longitude
                          rotation_on_lat,   &   !euler angle on latitude
                                geo_lon_u,   &   !geographical longitude in degrees
                                geo_lat_u,   &   !geographical latitude  in degrees
                                     dxt,    &   !metrical coefficient on x
                                     dyh,    &   !metrical coefficient on x
                              rotvec_coeff,  &   !rotation coefficients for vector transform
                                    rlh_s,   &   !coriolis main term (sin)
                                    rlh_c,   &   !coriolis second term (cos)
                                        0,   &   !key to compute rotation coefficients (0/1)
                                        0,   &   !key to compute coriolis terms (0/1)
                                 bnd_x1+1,   &   !first significant point in x-direction (output)
                                 bnd_x2-1,   &   ! last significant point in x-direction (output)
                                 bnd_y1+1,   &   !first significant point in y-direction (output)
                                 bnd_y2-1    )   ! last significant point in y-direction (output)

!On V-grid
       call grid_parameters_spherical (xt,   &   !model x-coordinate in degrees
                                       yv,   &   !model y-coordinate in degrees
                                   bnd_x1,   &   !left   boundary of arrays
                                   bnd_x2,   &   !right  boundary of arrays
                                   bnd_y1,   &   !lower  boundary of arrays
                                   bnd_y2,   &   !upper  boundary of arrays
                          rotation_on_lon,   &   !euler angle on longitude
                          rotation_on_lat,   &   !euler angle on latitude
                                geo_lon_v,   &   !geographical longitude in degrees
                                geo_lat_v,   &   !geographical latitude  in degrees
                                     dxh,    &   !metrical coefficient on x
                                     dyt,    &   !metrical coefficient on x
                              rotvec_coeff,  &   !rotation coefficients for vector transform
                                    rlh_s,   &   !coriolis main term (sin)
                                    rlh_c,   &   !coriolis second term (cos)
                                        0,   &   !key to compute rotation coefficients (0/1)
                                        0,   &   !key to compute coriolis terms (0/1)
                                 bnd_x1+1,   &   !first significant point in x-direction (output)
                                 bnd_x2-1,   &   ! last significant point in x-direction (output)
                                 bnd_y1+1,   &   !first significant point in y-direction (output)
                                 bnd_y2-1    )   ! last significant point in y-direction (output)

!On H-grid
       call grid_parameters_spherical (xu,   &   !model x-coordinate in degrees
                                       yv,   &   !model y-coordinate in degrees
                                   bnd_x1,   &   !left   boundary of arrays
                                   bnd_x2,   &   !right  boundary of arrays
                                   bnd_y1,   &   !lower  boundary of arrays
                                   bnd_y2,   &   !upper  boundary of arrays
                          rotation_on_lon,   &   !euler angle on longitude
                          rotation_on_lat,   &   !euler angle on latitude
                                geo_lon_h,   &   !geographical longitude in degrees
                                geo_lat_h,   &   !geographical latitude  in degrees
                                     dxb,    &   !metrical coefficient on x
                                     dyb,    &   !metrical coefficient on x
                              rotvec_coeff,  &   !rotation coefficients for vector transform
                                    rlh_s,   &   !coriolis main term (sin)
                                    rlh_c,   &   !coriolis second term (cos)
                                        0,   &   !key to compute rotation coefficients (0/1)
                                        1,   &   !key to compute coriolis terms (0/1)
                                 bnd_x1+1,   &   !first significant point in x-direction (output)
                                 bnd_x2-1,   &   ! last significant point in x-direction (output)
                                 bnd_y1+1,   &   !first significant point in y-direction (output)
                                 bnd_y2-1    )   ! last significant point in y-direction (output)

      elseif(curve_grid==2) then   !in case of curve grid

!On T-grid
       call grid_parameters_curvilinear (xt,   &   !model x-coordinate in degrees
                                         yt,   &   !model y-coordinate in degrees
                                     bnd_x1,   &   !left   boundary of arrays
                                     bnd_x2,   &   !right  boundary of arrays
                                     bnd_y1,   &   !lower  boundary of arrays
                                     bnd_y2,   &   !upper  boundary of arrays
                                     x_pole,   &   !geo longitude of new north pole
                                     y_pole,   &   !geo latitude  of new north pole
                                     p_pole,   &   !geo longitude of new south pole
                                     q_pole,   &   !geo latitude  of new south pole
                                  geo_lon_t,   &   !geographical longitude in degrees
                                  geo_lat_t,   &   !geographical latitude  in degrees
                                        dx,    &   !metrical coefficient on x
                                        dy,    &   !metrical coefficient on x
                                rotvec_coeff,  &   !rotation coefficients for vector transform
                                      rlh_s,   &   !coriolis main term (sin)
                                      rlh_c,   &   !coriolis second term (cos)
                                          1,   &   !key to compute rotation coefficients (0/1)
                                          0,   &   !key to compute coriolis terms (0/1)
                                   nx_start-1, &   !first significant point in x-direction (output)
                                   nx_end+1,   &   ! last significant point in x-direction (output)
                                   ny_start-1, &   !first significant point in y-direction (output)
                                   ny_end+1    )   ! last significant point in y-direction (output)

!On U-grid
       call grid_parameters_curvilinear(xu,   &   !model x-coordinate in degrees
                                        yt,   &   !model y-coordinate in degrees
                                    bnd_x1,   &   !left   boundary of arrays
                                    bnd_x2,   &   !right  boundary of arrays
                                    bnd_y1,   &   !lower  boundary of arrays
                                    bnd_y2,   &   !upper  boundary of arrays
                                    x_pole,   &   !geo longitude of new north pole
                                    y_pole,   &   !geo latitude  of new north pole
                                    p_pole,   &   !geo longitude of new south pole
                                    q_pole,   &   !geo latitude  of new south pole
                                 geo_lon_u,   &   !geographical longitude in degrees
                                 geo_lat_u,   &   !geographical latitude  in degrees
                                      dxt,    &   !metrical coefficient on x
                                      dyh,    &   !metrical coefficient on x
                               rotvec_coeff,  &   !rotation coefficients for vector transform
                                     rlh_s,   &   !coriolis main term (sin)
                                     rlh_c,   &   !coriolis second term (cos)
                                         0,   &   !key to compute rotation coefficients (0/1)
                                         0,   &   !key to compute coriolis terms (0/1)
                                  nx_start-1, &   !first significant point in x-direction (output)
                                  nx_end+1,   &   ! last significant point in x-direction (output)
                                  ny_start-1, &   !first significant point in y-direction (output)
                                  ny_end+1    )   ! last significant point in y-direction (output)

!On V-grid
       call grid_parameters_curvilinear(xt,   &   !model x-coordinate in degrees
                                        yv,   &   !model y-coordinate in degrees
                                    bnd_x1,   &   !left   boundary of arrays
                                    bnd_x2,   &   !right  boundary of arrays
                                    bnd_y1,   &   !lower  boundary of arrays
                                    bnd_y2,   &   !upper  boundary of arrays
                                    x_pole,   &   !geo longitude of new north pole
                                    y_pole,   &   !geo latitude  of new north pole
                                    p_pole,   &   !geo longitude of new south pole
                                    q_pole,   &   !geo latitude  of new south pole
                                 geo_lon_v,   &   !geographical longitude in degrees
                                 geo_lat_v,   &   !geographical latitude  in degrees
                                      dxh,    &   !metrical coefficient on x
                                      dyt,    &   !metrical coefficient on x
                               rotvec_coeff,  &   !rotation coefficients for vector transform
                                     rlh_s,   &   !coriolis main term (sin)
                                     rlh_c,   &   !coriolis second term (cos)
                                         0,   &   !key to compute rotation coefficients (0/1)
                                         0,   &   !key to compute coriolis terms (0/1)
                                  nx_start-1, &   !first significant point in x-direction (output)
                                  nx_end+1,   &   ! last significant point in x-direction (output)
                                  ny_start-1, &   !first significant point in y-direction (output)
                                  ny_end+1    )   ! last significant point in y-direction (output)

!On H-grid
       call grid_parameters_curvilinear(xu,   &   !model x-coordinate in degrees
                                        yv,   &   !model y-coordinate in degrees
                                    bnd_x1,   &   !left   boundary of arrays
                                    bnd_x2,   &   !right  boundary of arrays
                                    bnd_y1,   &   !lower  boundary of arrays
                                    bnd_y2,   &   !upper  boundary of arrays
                                    x_pole,   &   !geo longitude of new north pole
                                    y_pole,   &   !geo latitude  of new north pole
                                    p_pole,   &   !geo longitude of new south pole
                                    q_pole,   &   !geo latitude  of new south pole
                                 geo_lon_h,   &   !geographical longitude in degrees
                                 geo_lat_h,   &   !geographical latitude  in degrees
                                      dxb,    &   !metrical coefficient on x
                                      dyb,    &   !metrical coefficient on x
                               rotvec_coeff,  &   !rotation coefficients for vector transform
                                     rlh_s,   &   !coriolis main term (sin)
                                     rlh_c,   &   !coriolis second term (cos)
                                         0,   &   !key to compute rotation coefficients (0/1)
                                         1,   &   !key to compute coriolis terms (0/1)
                                  nx_start-1, &   !first significant point in x-direction (output)
                                  nx_end+1,   &   ! last significant point in x-direction (output)
                                  ny_start-1, &   !first significant point in y-direction (output)
                                  ny_end+1    )   ! last significant point in y-direction (output)

      end if
!-----end of metric initialization------------------------------------------

!array4=sngl(dx)
!call wdstd(' ','dx.dat',1,array4,lu1,nx,ny,1,mmm,mm,nnn,nn,1,1,ierr)
!call ctl_file_write('dx.dat',    &     !file name
!                    undef,    &     !value for undefined points
!                    mm-mmm+1,       &     !x-dimension
!                    nn-nnn+1,       &     !y-dimension
!                    1,        &     !z-dimension
!                    1,        &     !t-dimension
!                    xgr_type,       &     !x-grid type (0 - linear, 1 - levels)
!                    xt(nx_start),   &     !first x-value (if linear) or x-array (if levels)
!                    dxst,       &     !x-step (if linear)
!                    ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
!                    yt(ny_start),       &     !first y-value (if linear) or x-array (if levels)
!                    dyst,       &     !y-step (if linear)
!                    0,          &     !z-grid type (0 - linear, 1 - levels)
!                    0.0d0,        &     !first z-value (if linear) or x-array (if levels)
!                    1.0d0,        &     !z-step (if linear)
!                    0,          &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
!                   2000,       &     !year   of the first field
!                    1,          &     !month  of the first field
!                    1,          &     !day    of the first field
!                    0,          &     !hour   of the first field
!                    0,          &     !minute of the first field
!                    3600.0,     &     !time step (in seconds)
!                    'Grid step in X-direction at T-grid, m',    &     !title of dataset
!                    'dx'   )     !variable name
!
!
!array4=sngl(dy)
!call wdstd(' ','dy.dat',1,array4,lu1,nx,ny,1,mmm,mm,nnn,nn,1,1,ierr)
!call ctl_file_write('dy.dat',    &     !file name
!                    undef,    &     !value for undefined points
!                    mm-mmm+1,       &     !x-dimension
!                    nn-nnn+1,       &     !y-dimension
!                    1,        &     !z-dimension
!                    1,        &     !t-dimension
!                    xgr_type,       &     !x-grid type (0 - linear, 1 - levels)
!                    xt(nx_start),   &     !first x-value (if linear) or x-array (if levels)
!                    dxst,       &     !x-step (if linear)
!                    ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
!                    yt(ny_start),       &     !first y-value (if linear) or x-array (if levels)
!                    dyst,       &     !y-step (if linear)
!                    0,          &     !z-grid type (0 - linear, 1 - levels)
!                    0.0d0,        &     !first z-value (if linear) or x-array (if levels)
!                    1.0d0,        &     !z-step (if linear)
!                    0,          &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
!                    2000,       &     !year   of the first field
!                    1,          &     !month  of the first field
!                    1,          &     !day    of the first field
!                    0,          &     !hour   of the first field
!                    0,          &     !minute of the first field
!                    3600.0,     &     !time step (in seconds)
!                    'Grid step in Y-direction at T-grid, m',    &     !title of dataset
!                    'dy'   )     !variable name
!
!array4=sngl(geo_lon_t)
!call wdstd(' ','geolonT.dat',1,array4,lu1,nx,ny,1,mmm,mm,nnn,nn,1,1,ierr)
!call ctl_file_write('geolonT.dat',    &     !file name
!                    undef,    &     !value for undefined points
!                    mm-mmm+1,       &     !x-dimension
!                    nn-nnn+1,       &     !y-dimension
!                    1,        &     !z-dimension
!                    1,        &     !t-dimension
!                    xgr_type,       &     !x-grid type (0 - linear, 1 - levels)
!                    xt(nx_start),   &     !first x-value (if linear) or x-array (if levels)
!                    dxst,       &     !x-step (if linear)
!                    ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
!                    yt(ny_start),       &     !first y-value (if linear) or x-array (if levels)
!                    dyst,       &     !y-step (if linear)
!                    0,          &     !z-grid type (0 - linear, 1 - levels)
!                    0.0d0,        &     !first z-value (if linear) or x-array (if levels)
!                    1.0d0,        &     !z-step (if linear)
!                    0,          &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
!                    2000,       &     !year   of the first field
!                    1,          &     !month  of the first field
!                    1,          &     !day    of the first field
!                    0,          &     !hour   of the first field
!                    0,          &     !minute of the first field
!                    3600.0,     &     !time step (in seconds)
!                    'geographical longitude, degrees',    &     !title of dataset
!                    'glon'   )     !variable name
!
!array4=sngl(geo_lat_t)
!call wdstd(' ','geolatT.dat',1,array4,lu1,nx,ny,1,mmm,mm,nnn,nn,1,1,ierr)
!call ctl_file_write('geolatT.dat',    &     !file name
!                    undef,    &     !value for undefined points
!                    mm-mmm+1,       &     !x-dimension
!                    nn-nnn+1,       &     !y-dimension
!                    1,        &     !z-dimension
!                    1,        &     !t-dimension
!                    xgr_type,       &     !x-grid type (0 - linear, 1 - levels)
!                    xt(nx_start),   &     !first x-value (if linear) or x-array (if levels)
!                    dxst,       &     !x-step (if linear)
!                    ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
!                    yt(ny_start),       &     !first y-value (if linear) or x-array (if levels)
!                    dyst,       &     !y-step (if linear)
!                    0,          &     !z-grid type (0 - linear, 1 - levels)
!                    0.0d0,        &     !first z-value (if linear) or x-array (if levels)
!                    1.0d0,        &     !z-step (if linear)
!                    0,          &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
!                    2000,       &     !year   of the first field
!                    1,          &     !month  of the first field
!                    1,          &     !day    of the first field
!                    0,          &     !hour   of the first field
!                    0,          &     !minute of the first field
!                    3600.0,     &     !time step (in seconds)
!                    'geographical latitude, degrees',    &     !title of dataset
!                    'glat'   )     !variable name
!
!array4=sngl(rlh_s)
!call wdstd(' ','rlh_s.dat',1,array4,lu1,nx,ny,1,mmm-1,mm,nnn-1,nn,1,1,ierr)
!call ctl_file_write('rlh_s.dat',    &     !file name
!                    undef,    &     !value for undefined points
!                    mm-mmm+2,       &     !x-dimension
!                    nn-nnn+2,       &     !y-dimension
!                    1,        &     !z-dimension
!                    1,        &     !t-dimension
!                    xgr_type,       &     !x-grid type (0 - linear, 1 - levels)
!                    xu(nx_start),   &     !first x-value (if linear) or x-array (if levels)
!                    dxst,       &     !x-step (if linear)
!                    ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
!                    yv(ny_start),       &     !first y-value (if linear) or x-array (if levels)
!                    dyst,       &     !y-step (if linear)
!                    0,          &     !z-grid type (0 - linear, 1 - levels)
!                    0.0d0,      &     !first z-value (if linear) or x-array (if levels)
!                    1.0d0,      &     !z-step (if linear)
!                    0,          &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
!                    2000,       &     !year   of the first field
!                    1,          &     !month  of the first field
!                    1,          &     !day    of the first field
!                    0,          &     !hour   of the first field
!                    0,          &     !minute of the first field
!                    3600.0,     &     !time step (in seconds)
!                    'main Coriolis parameter, rad/s',    &     !title of dataset
!                    'rlh_s'   )     !variable name
!
!array4=sngl(rlh_c)
!call wdstd(' ','rlh_c.dat',1,array4,lu1,nx,ny,1,mmm-1,mm,nnn-1,nn,1,1,ierr)
!call ctl_file_write('rlh_c.dat',    &     !file name
!                    undef,    &     !value for undefined points
!                    mm-mmm+2,       &     !x-dimension
!                    nn-nnn+2,       &     !y-dimension
!                    1,        &     !z-dimension
!                    1,        &     !t-dimension
!                    xgr_type,       &     !x-grid type (0 - linear, 1 - levels)
!                    xu(nx_start),   &     !first x-value (if linear) or x-array (if levels)
!                    dxst,       &     !x-step (if linear)
!                    ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
!                    yv(ny_start),       &     !first y-value (if linear) or x-array (if levels)
!                    dyst,       &     !y-step (if linear)
!                   0,          &     !z-grid type (0 - linear, 1 - levels)
!                    0.0d0,        &     !first z-value (if linear) or x-array (if levels)
!                    1.0d0,        &     !z-step (if linear)
!                    0,          &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
!                    2000,       &     !year   of the first field
!                    1,          &     !month  of the first field
!                    1,          &     !day    of the first field
!                    0,          &     !hour   of the first field
!                    0,          &     !minute of the first field
!                    3600.0,     &     !time step (in seconds)
!                    'second Coriolis parameter, rad/s',    &     !title of dataset
!                    'rlh_c'   )     !variable name

!open(35,file='geo_coordinates.txt')
!write(35,*) 'Number on x ; Number on y ;   geolon  ;  geolat    ;  cos of rot;  sin of rot'
! do n=nnn,nn
!  do m=mmm,mm
!    write(35,'(i6,a,i6,a,f13.6,a,f13.6,a,f13.6,a,f13.6)') m-mmm+1,';', n-nnn+1,';', geo_lon_t(m,n),';',    &
!                                            geo_lat_t(m,n),';', rotvec_coeff(m,n,1),';', rotvec_coeff(m,n,2)
!  enddo
! enddo
!close(35)

 endsubroutine basinpar

!====================================================================
 subroutine grid_parameters_carthesian(x_mod,   &   !model x-coordinate in degrees
                                       y_mod,   &   !model y-coordinate in degrees
                                      bnd_x1,   &   !left   boundary of arrays
                                      bnd_x2,   &   !right  boundary of arrays
                                      bnd_y1,   &   !lower  boundary of arrays
                                      bnd_y2,   &   !upper  boundary of arrays
                                     geo_lon,   &   !geographical longitude in degrees
                                     geo_lat,   &   !geographical latitude  in degrees
                                     metr_x,    &   !metrical coefficient on x
                                     metr_y,    &   !metrical coefficient on x
                                     rot_coef,  &   !rotation coefficients for vector transform
                                     cor_sin,   &   !coriolis main term (sin)
                                     cor_cos,   &   !coriolis second term (cos)
                                     key_rot,   &   !key to compute rotation coefficients (0/1)
                                     key_cor,   &   !key to compute coriolis terms (0/1)
                                     mmm_out,   &   !first significant point in x-direction (output)
                                     mm_out,    &   ! last significant point in x-direction (output)
                                     nnn_out,   &   !first significant point in y-direction (output)
                                     nn_out)        !last significant point in y-direction (output)

 implicit none
 integer bnd_x1,bnd_x2,bnd_y1,bnd_y2
 integer mmm_out, mm_out, nnn_out, nn_out

 real(8), parameter :: dpip180=3.1415926535897/180.0d0 !for degrees to radians convers
 real(8)          x_mod(bnd_x1:bnd_x2),   &
                  y_mod(bnd_y1:bnd_y2),   &
  geo_lon(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &
  geo_lat(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &
   metr_x(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &
   metr_y(bnd_x1:bnd_x2,bnd_y1:bnd_y2)

 real(8) rot_coef(bnd_x1:bnd_x2,bnd_y1:bnd_y2,4),     &
          cor_sin(bnd_x1:bnd_x2,bnd_y1:bnd_y2),       &
          cor_cos(bnd_x1:bnd_x2,bnd_y1:bnd_y2)

 integer key_rot, key_cor
 integer m,n

 !$omp parallel do private (m,n)
      do n=nnn_out,nn_out
       do m=mmm_out,mm_out

!         necessary latitude
	    geo_lat(m,n) = y_mod(n)

!         necessary longitude
          geo_lon(m,n)= x_mod(m)

          metr_x(m,n)=metr_x(m,n)*1.0d0
          metr_y(m,n)=metr_y(m,n)*1.0d0

       if(key_rot==1) then
           !--------definition of angles between parallels-----------------------
	  rot_coef(m,n,1) = 1.0d0
        rot_coef(m,n,2) = 0.0d0
        rot_coef(m,n,3)=-rot_coef(m,n,2)
        rot_coef(m,n,4)= rot_coef(m,n,1)
       endif

       if(key_cor==1) then
        cor_sin(m,n)=cor_sin(m,n)/dsqrt(2.0d0)
        cor_cos(m,n)=cor_cos(m,n)/dsqrt(2.0d0)
       endif

       enddo
      enddo
 !$omp end parallel do

 endsubroutine grid_parameters_carthesian

!====================================================================
 subroutine grid_parameters_spherical (x_mod,   &   !model x-coordinate in degrees
                                       y_mod,   &   !model y-coordinate in degrees
                                      bnd_x1,   &   !left   boundary of arrays
                                      bnd_x2,   &   !right  boundary of arrays
                                      bnd_y1,   &   !lower  boundary of arrays
                                      bnd_y2,   &   !upper  boundary of arrays
                             rotation_on_lon,   &   !euler angle on longitude
                             rotation_on_lat,   &   !euler angle on latitude
                                     geo_lon,   &   !geographical longitude in degrees
                                     geo_lat,   &   !geographical latitude  in degrees
                                     metr_x,    &   !metrical coefficient on x
                                     metr_y,    &   !metrical coefficient on x
                                     rot_coef,  &   !rotation coefficients for vector transform
                                     cor_sin,   &   !coriolis main term (sin)
                                     cor_cos,   &   !coriolis second term (cos)
                                     key_rot,   &   !key to compute rotation coefficients (0/1)
                                     key_cor,   &   !key to compute coriolis terms (0/1)
                                     mmm_out,   &   !first significant point in x-direction (output)
                                     mm_out,    &   ! last significant point in x-direction (output)
                                     nnn_out,   &   !first significant point in y-direction (output)
                                     nn_out)        !last significant point in y-direction (output)

 implicit none
 integer bnd_x1,bnd_x2,bnd_y1,bnd_y2
 integer mmm_out, mm_out, nnn_out, nn_out

 real(8), parameter :: dpip180=3.1415926535897/180.0d0 !for degrees to radians convers
 real(8)          x_mod(bnd_x1:bnd_x2),   &
                  y_mod(bnd_y1:bnd_y2),   &
  geo_lon(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &
  geo_lat(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &
  metr_x(bnd_x1:bnd_x2,bnd_y1:bnd_y2),    &
  metr_y(bnd_x1:bnd_x2,bnd_y1:bnd_y2)

 real(8) rot_coef(bnd_x1:bnd_x2,bnd_y1:bnd_y2,4),     &
          cor_sin(bnd_x1:bnd_x2,bnd_y1:bnd_y2),       &
          cor_cos(bnd_x1:bnd_x2,bnd_y1:bnd_y2)

 integer key_rot, key_cor

 real(8) rotation_on_lon, rotation_on_lat

 real(8) sin_lon, sin_lat, cos_lon, cos_lat, lat_mod   	!auxilary variables
 real(8) free_term_coslon, free_term_sinlon

 real(8), parameter:: lat_extr= 89.99999d0
 real(8) sinlat_extr, coslat_extr, sum_rot_coef, sum_sincos

 integer m,n

  coslat_extr= dcos(dpip180*lat_extr)
  sinlat_extr= dsin(dpip180*lat_extr)

 !$omp parallel do private (m,n,sin_lat,cos_lat,free_term_coslon,free_term_sinlon,cos_lon,sin_lon,  &
 !$omp                      sum_rot_coef, sum_sincos, lat_mod)
      do n=nnn_out,nn_out
          lat_mod=max(min(y_mod(n),lat_extr),-lat_extr)
       do m=mmm_out,mm_out
          sin_lat = dsin(dpip180*y_mod(n)) * dcos(dpip180*rotation_on_lat)                           &
                  + dcos(dpip180*x_mod(m)) * dcos(dpip180*y_mod(n)) * dsin(dpip180*rotation_on_lat)

          sin_lat=min(max(sin_lat, -sinlat_extr),sinlat_extr)
          cos_lat = dsqrt(1d0-sin_lat**2)

!         necessary latitude
	    geo_lat(m,n) = dasin(sin_lat)/dpip180

	    free_term_coslon =(  dcos(dpip180*x_mod(m)) * dcos(dpip180*y_mod(n)) * dcos(dpip180*rotation_on_lat)   &
                            -  dsin(dpip180*y_mod(n)) * dsin(dpip180*rotation_on_lat)  )  / cos_lat

	    free_term_sinlon =(  dsin(dpip180*x_mod(m)) * dcos(dpip180*y_mod(n))  ) / cos_lat

		cos_lon=free_term_coslon*dcos(dpip180*rotation_on_lon)     &
                   -free_term_sinlon*dsin(dpip180*rotation_on_lon)

		sin_lon=free_term_sinlon*dcos(dpip180*rotation_on_lon)     &
                   +free_term_coslon*dsin(dpip180*rotation_on_lon)

          sum_sincos=max(dsqrt(cos_lon**2+sin_lon**2),1d-10)
          cos_lon=cos_lon/sum_sincos
          sin_lon=sin_lon/sum_sincos

!         necessary longitude
          geo_lon(m,n)=dsign(dacos(cos_lon)/dpip180,sin_lon)

          metr_x(m,n)=metr_x(m,n)*dcos(dpip180*lat_mod)
          metr_y(m,n)=metr_y(m,n)*1.0d0

       if(key_rot==1) then
           !--------definition of angles between parallels-----------------------
	  rot_coef(m,n,1) = (  cos_lat*dcos(dpip180*rotation_on_lat) + sin_lat*dsin(dpip180*rotation_on_lat)    &
                         * ( cos_lon*dcos(dpip180*rotation_on_lon) + sin_lon*dsin(dpip180*rotation_on_lon) )  &
                                            )  / dcos(dpip180*lat_mod)

        rot_coef(m,n,2) = (  -dsin(dpip180*rotation_on_lat)                                                     &
                         * ( sin_lon*dcos(dpip180*rotation_on_lon) - cos_lon*dsin(dpip180*rotation_on_lon) )  )      &
                                     / dcos(dpip180*lat_mod)

	  rot_coef(m,n,3)=-rot_coef(m,n,2)
        rot_coef(m,n,4)= rot_coef(m,n,1)
        sum_rot_coef=   max(dsqrt(rot_coef(m,n,1)*rot_coef(m,n,4)-rot_coef(m,n,2)*rot_coef(m,n,3)),1d-10)
        rot_coef(m,n,:)=rot_coef(m,n,:)/sum_rot_coef
       endif

       if(key_cor==1) then
        cor_sin(m,n)=cor_sin(m,n)*sin_lat
        cor_cos(m,n)=cor_cos(m,n)*cos_lat
       endif

       enddo
      enddo
 !$omp end parallel do

 endsubroutine grid_parameters_spherical

 !====================================================================
 subroutine grid_parameters_curvilinear (x_mod,   &   !model x-coordinate in degrees
                                         y_mod,   &   !model y-coordinate in degrees
                                        bnd_x1,   &   !left   boundary of arrays
                                        bnd_x2,   &   !right  boundary of arrays
                                        bnd_y1,   &   !lower  boundary of arrays
                                        bnd_y2,   &   !upper  boundary of arrays
                                        x_pole,   &   !geo longitude of new north pole
                                        y_pole,   &   !geo latitude  of new north pole
                                        p_pole,   &   !geo longitude of new south pole
                                        q_pole,   &   !geo latitude  of new south pole
                                       geo_lon,   &   !geographical longitude in degrees
                                       geo_lat,   &   !geographical latitude  in degrees
                                       metr_x,    &   !metrical coefficient on x
                                       metr_y,    &   !metrical coefficient on x
                                       rot_coef,  &   !rotation coefficients for vector transform
                                       cor_sin,   &   !coriolis main term (sin)
                                       cor_cos,   &   !coriolis second term (cos)
                                       key_rot,   &   !key to compute rotation coefficients (0/1)
                                       key_cor,   &   !key to compute coriolis terms (0/1)
                                        mmm_out,  &   !first significant point in x-direction (output)
                                        mm_out,   &   ! last significant point in x-direction (output)
                                        nnn_out,  &   !first significant point in y-direction (output)
                                        nn_out)       !last significant point in y-direction (output)
 implicit none

 integer bnd_x1,bnd_x2,bnd_y1,bnd_y2
 integer mmm_out, mm_out, nnn_out, nn_out

 real(8), parameter:: dpip180=3.1415926535897/180.0d00 !for degrees to radians convers

 real(8)          x_mod(bnd_x1:bnd_x2),   &
                  y_mod(bnd_y1:bnd_y2),   &
  geo_lon(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &
  geo_lat(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &
  metr_x(bnd_x1:bnd_x2,bnd_y1:bnd_y2),    &
  metr_y(bnd_x1:bnd_x2,bnd_y1:bnd_y2)

 real(8) rot_coef(bnd_x1:bnd_x2,bnd_y1:bnd_y2,4),   &
          cor_sin(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &
          cor_cos(bnd_x1:bnd_x2,bnd_y1:bnd_y2)

 integer key_rot, key_cor

 real(8) x_pole, y_pole, p_pole, q_pole
 real(8) y_pole1, q_pole1, r3d,r2d

 real(8) sin_lon, sin_lat, cos_lon, cos_lat, lat_mod	!auxilary variables
 real(8) a,b,s,t,a0,b0,s0,t0                    !auxilary variables

 real(8) num1,num2,numa,numb,denom1
 real(8) numd1,numd2,numd3,numd4,numas,numat,numbs,numbt
 real(8) alpha_scale
 real(8) xn,yn,zn,xs,ys,zs,xm,ym,zm,lm,phm, sinphm,coslm,sinlm, phm1

 real(8) dx_da, dx_db, dy_da, dy_db, da_ds, da_dt,   &
         db_ds, db_dt, ds_dp, ds_dq, dt_dp, dt_dq,   &
         da_dp, da_dq, db_dp, db_dq, dx_dp, dx_dq,   &
         dy_dp, dy_dq, det, hp_divide_r, hq_divide_r
 real(8) df(2,2), dfm1(2,2)

 real(8), parameter:: lat_extr= 89.99999d0
 real(8) sinlat_extr, coslat_extr, sum_rot_coef, sum_sincos

 integer m,n

  y_pole1=min(max(y_pole, -lat_extr),lat_extr)
  q_pole1=min(max(q_pole, -lat_extr),lat_extr)

  coslat_extr= dcos(dpip180*lat_extr)
  sinlat_extr= dsin(dpip180*lat_extr)

      xn=dcos(x_pole*dpip180)*dcos(y_pole*dpip180)
      yn=dsin(x_pole*dpip180)*dcos(y_pole*dpip180)
      zn=dsin(y_pole*dpip180)

      xs=dcos(p_pole*dpip180)*dcos(q_pole*dpip180)
      ys=dsin(p_pole*dpip180)*dcos(q_pole*dpip180)
      zs=dsin(q_pole*dpip180)

      xm=(xn+xs)/2.0d0
      ym=(yn+ys)/2.0d0
      zm=(zn+zs)/2.0d0

      r3d=max(dsqrt(xm**2+ym**2+zm**2),1d-10)
      r2d=max(dsqrt(xm**2+ym**2),1d-10)

      sinphm=zm/r3d
      sinlm= ym/r2d
      coslm= xm/r2d

      sinphm=min(max(sinphm, -sinlat_extr),sinlat_extr)
      phm=dasin(sinphm)/dpip180

      sum_sincos=max(dsqrt(coslm**2+sinlm**2),1d-10)
      coslm=coslm/sum_sincos
      sinlm=sinlm/sum_sincos

      lm=dsign(dacos(coslm)/dpip180,sinlm)

      s0 = 2.0d0*dtan((45.0d0 + y_pole1/2.0d0)*dpip180) *dcos(x_pole*dpip180)
      t0 = 2.0d0*dtan((45.0d0 + y_pole1/2.0d0)*dpip180) *dsin(x_pole*dpip180)
      a0 = 2.0d0*dtan((45.0d0 + q_pole1/2.0d0)*dpip180) *dcos(p_pole*dpip180)
      b0 = 2.0d0*dtan((45.0d0 + q_pole1/2.0d0)*dpip180) *dsin(p_pole*dpip180)

      alpha_scale=1.0d0

      phm1=min(max(phm, -lat_extr),lat_extr)

      S = 2.0d0*dtan((45.0d0 + phm1/2.0d0)*dpip180) *dcos(lm*dpip180)
      T = 2.0d0*dtan((45.0d0 + phm1/2.0d0)*dpip180) *dsin(lm*dpip180)

      num1=(S-alpha_scale*A0)*(S-alpha_scale*S0) + (T-alpha_scale*B0)*(T-alpha_scale*T0)
      num2=(T-alpha_scale*B0)*(S-alpha_scale*S0) - (S-alpha_scale*A0)*(T-alpha_scale*T0)

      numa=s0*num1-t0*num2
      numb=s0*num2+t0*num1

      denom1=(S-alpha_scale*S0)**2+(T-alpha_scale*T0)**2

      a=numa/denom1
      b=numb/denom1

      alpha_scale=2./dsqrt(a**2+b**2)

      write(*,*) 'alpha-scale is ', alpha_scale

!$omp parallel do private(m,n,s,t,num1,num2,numa,numb,denom1,  &
!$omp      a,b,cos_lon,sin_lon,cos_lat,sin_lat, &
!$omp      dx_da,dx_db,dy_da,dy_db,numd1,numd2, &
!$omp      numd3,numd4,numas,numat,numbs,numbt, &
!$omp      da_ds,da_dt,db_ds,db_dt,ds_dp,ds_dq,dt_dp,dt_dq,  &
!$omp      da_dp,da_dq,db_dp,db_dq,dx_dp,dx_dq,dy_dp,dy_dq,  &
!$omp      dfm1,det,df,hp_divide_r,hq_divide_r, sum_rot_coef, sum_sincos, lat_mod)
      do n=nnn_out,nn_out
          lat_mod=max(min(y_mod(n),lat_extr),-lat_extr)
       do m=mmm_out,mm_out

!trasformation from new to old grid
      s = 2.0d0*dtan((45.0d0 + lat_mod/2.0d0)*dpip180) *dcos(x_mod(m)*dpip180)
      t = 2.0d0*dtan((45.0d0 + lat_mod/2.0d0)*dpip180) *dsin(x_mod(m)*dpip180)

      num1=(S-alpha_scale*A0)*(S-alpha_scale*S0) + (T-alpha_scale*B0)*(T-alpha_scale*T0)
      num2=(T-alpha_scale*B0)*(S-alpha_scale*S0) - (S-alpha_scale*A0)*(T-alpha_scale*T0)

      numa=s0*num1-t0*num2
      numb=s0*num2+t0*num1

      denom1=(S-alpha_scale*S0)**2+(T-alpha_scale*T0)**2

      a=numa/denom1
      b=numb/denom1

!     necessary latitude

      sin_lat=(a**2+b**2-4.0d0)/(a**2+b**2+4.0d0)
      sin_lat=min(max(sin_lat, -sinlat_extr),sinlat_extr)
	cos_lat=dsqrt(1.0d0-sin_lat**2)

      geo_lat(m,n)=dasin(sin_lat)/dpip180

!     necessary longitude

      cos_lon=a/dsqrt(a**2+b**2)
	sin_lon=b/dsqrt(a**2+b**2)

      sum_sincos=max(dsqrt(cos_lon**2+sin_lon**2),1d-10)
      cos_lon=cos_lon/sum_sincos
      sin_lon=sin_lon/sum_sincos

      geo_lon(m,n)=dsign(dacos(cos_lon)/dpip180,sin_lon)


!     differential of transformation

	dx_da = -b / (a**2 + b**2)
	dx_db =  a / (a**2 + b**2)

	dy_da = a / ( dsqrt(a**2 + b**2) * (1.0d0+ (a**2 + b**2)/4.0d0))
	dy_db = b / ( dsqrt(a**2 + b**2) * (1.0d0+ (a**2 + b**2)/4.0d0))


      numd1=s-alpha_scale*s0+s-alpha_scale*a0
      numd2=t-alpha_scale*t0+t-alpha_scale*b0
      numd3=alpha_scale*(t0-b0)
      numd4=alpha_scale*(a0-s0)

      numas=s0*numd1-t0*numd3
      numat=s0*numd2-t0*numd4
      numbs=t0*numd1+s0*numd3
      numbt=t0*numd2+s0*numd4

      da_ds=numas/denom1-numa*2.0d0*(s-alpha_scale*s0)/(denom1*denom1)
      da_dt=numat/denom1-numa*2.0d0*(t-alpha_scale*t0)/(denom1*denom1)
      db_ds=numbs/denom1-numb*2.0d0*(s-alpha_scale*s0)/(denom1*denom1)
      db_dt=numbt/denom1-numb*2.0d0*(t-alpha_scale*t0)/(denom1*denom1)


      ds_dp = -2.0d0*dtan((45.0d0 + lat_mod/2.0d0)*dpip180) *dsin(x_mod(m)*dpip180)
      ds_dq = dcos(x_mod(m)*dpip180) /(dcos((45.0d0 + lat_mod/2.0d0)*dpip180))**2
      dt_dp =  2.0d0*dtan((45.0d0 + lat_mod/2.0d0)*dpip180) *dcos(x_mod(m)*dpip180)
      dt_dq = dsin(x_mod(m)*dpip180) /(dcos((45.0d0 + lat_mod/2.0d0)*dpip180))**2

	da_dp = da_ds*ds_dp + da_dt*dt_dp
	da_dq = da_ds*ds_dq + da_dt*dt_dq
	db_dp = db_ds*ds_dp + db_dt*dt_dp
	db_dq = db_ds*ds_dq + db_dt*dt_dq

	dx_dp = dx_da*da_Dp + dx_db*db_dp
	dx_dq = dx_da*da_Dq + dx_db*db_dq
	dy_dp = dy_da*da_Dp + dy_db*db_dp
	dy_dq = dy_da*da_Dq + dy_db*db_dq

	dfm1(1,1)=dx_dp   !*dcos(ret_lat*dpip180)
	dfm1(1,2)=dx_dq   !*dcos(ret_lat*dpip180)
	dfm1(2,1)=dy_dp
	dfm1(2,2)=dy_dq


	det=dfm1(2,2)*dfm1(1,1)-dfm1(1,2)*dfm1(2,1)

	df(1,1)= dfm1(2,2)/det
	df(1,2)=-dfm1(1,2)/det
	df(2,1)=-dfm1(2,1)/det
	df(2,2)= dfm1(1,1)/det

	hp_divide_r = dsqrt((dx_dp*cos_lat)**2 + (dy_dp)**2)
	hq_divide_r = dsqrt((dx_dq*cos_lat)**2 + (dy_dq)**2)

       metr_x(m,n)=metr_x(m,n)*hp_divide_r
       metr_y(m,n)=metr_y(m,n)*hq_divide_r

      if(key_rot==1) then
     !--------definition of angles between parallels-----------------------
	  rot_coef(m,n,1)=df(1,1)*hp_divide_r/cos_lat
        rot_coef(m,n,2)=df(1,2)*hp_divide_r
	  rot_coef(m,n,3)=df(2,1)*hq_divide_r/cos_lat
        rot_coef(m,n,4)=df(2,2)*hq_divide_r

        sum_rot_coef=  max(dsqrt(rot_coef(m,n,1)*rot_coef(m,n,4)-rot_coef(m,n,2)*rot_coef(m,n,3)),1d-10)
        rot_coef(m,n,:)=rot_coef(m,n,:)/sum_rot_coef
      endif

      if(key_cor==1) then
        cor_sin(m,n)=cor_sin(m,n)*sin_lat
        cor_cos(m,n)=cor_cos(m,n)*cos_lat
      endif

       enddo
      enddo
 !$omp end parallel do

 endsubroutine grid_parameters_curvilinear
