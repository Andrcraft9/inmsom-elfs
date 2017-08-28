Program make_potential_temperature
use main_basin_pars
use basin_grid
implicit none

! making potential temperature from temperature in situ
! tt - temperature in situ
! ss - salinity
! tth - potential temperature
! den - potential density
 real(4) tt(nx,ny,nz), ss(nx,ny,nz), tth(nx,ny,nz), den(nx,ny,nz)
 character pathdat*32,fileh*32,filet*32,files*32,fileth*32, ftemask*32,yno
 character*128 fname
 integer n3da, m, n, k, l, ierr
 real(4) array4(nx,ny)

 call non_mpi_array_boundary_definition
 call model_grid_allocate

 tth=0.0
 den=0.0
  
 open(90,file='mkpot.par',status='old',err=89)
 read(90,'(a32)',err=90) pathdat !path to data(input)
 read(90,'(a32)',err=90) ftemask !file with temperature mask(input)
 read(90,'(a32)',err=90) fileh   !file with bottom data(input)
 read(90,'(a32)',err=90) filet   !file with sigma level temperature(in situ)
 read(90,'(a32)',err=90) files   !file with sigma level salinity
 read(90,'(a32)',err=90) fileth  !file with sigma level potential temperature
 read(90,*,err=90)     n3da      !number of 3d arrays:
 close(90)

 write(*,*)' path to data(input):                        ',pathdat
 write(*,*)' file with temperature mask:                 ',ftemask
 write(*,*)' file with bottom data(input):               ',fileh
 write(*,*)' file with sigma level temperature (in situ):',filet
 write(*,*)' file with sigma level salinity:             ',files
 write(*,*)' file with sigma level potential temperature:',fileth
 write(*,*)' number of 3d arrays:                        ',n3da
 write(*,*)' is it ok?'
 write(*,'(1x,a)')' (input if yes - y, if no - n)'
 read(*,'(a)') yno
 if(yno/='y'.and.yno/='Y') stop

! setting vertical t-,w- grid levels
 call vgrid

! define lame coefficients and coriolis parameter
 call basinpar

 ! grid parameter setting
 call gridcon(ftemask)

 array4=0.0
 call rdstd(' ',fileh,1,array4,lu,nx,ny,1, mmm,mm,nnn,nn,1,1,ierr)
 hhq=dble(array4)
   
   if(periodicity_x/=0) then
     call cyclize8_x(hhq,nx,ny,1,mmm,mm)
   end if

   if(periodicity_y/=0) then
     call cyclize8_y(hhq,nx,ny,1,nnn,nn)
   end if

!$omp parallel do private(m,n,k) 
   do n=1,ny
    do m=1,nx
     
     if(luu(m,n)>0.5) then
! interpolating hhq given on T-grid(lu) to hhh given on h-grid(luu).
       hhh(m,n)=( hhq(m,n)+ hhq(m+1,n)+ hhq(m,n+1)+ hhq(m+1,n+1) )/4.0d0
     endif

    end do
   end do
!$omp end parallel do

   if(periodicity_x/=0) then   
     call cyclize8_x(hhh, nx,ny,1,mmm,mm)
   end if
  
   if(periodicity_y/=0) then  
     call cyclize8_y(hhh, nx,ny,1,nnn,nn)
   end if


   do l = 1,n3da
    
    write(*,'(a,i4)')' Number of 3d array:',l
!read
    call rdstd(pathdat,filet,l,tt,lu,nx,ny,nz,mmm,mm,nnn,nn,1,nz,ierr)
    call rdstd(pathdat,files,l,ss,lu,nx,ny,nz,mmm,mm,nnn,nn,1,nz,ierr)

!----making potential temperature from temperature in situ------

    call t2th(tt,ss,den,tth)

!write

    call   wdstd(pathdat,fileth,   l,tth,lu,nx,ny,nz,mmm,mm,nnn,nn, 1,nz,ierr)
    call   wdstd(pathdat,'den.dat',l,den,lu,nx,ny,nz,mmm,mm,nnn,nn, 1,nz,ierr)

   enddo

   call fulfname(fname,pathdat,fileth,ierr)
   call ctl_file_write(fname,   &     !file name
                      -1e+32,   &     !value for undefined points
                    mm-mmm+1,   &     !x-dimension
                    nn-nnn+1,   &     !y-dimension
                          nz,   &     !z-dimension
                        n3da,   &     !t-dimension
                    xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                     xt(mmm),   &     !first x-value (if linear) or x-array (if levels)
                        dxst,   &     !x-step (if linear)
                    ygr_type,   &     !y-grid type (0 - linear, 1 - levels)
                     yt(nnn),   &     !first y-value (if linear) or x-array (if levels)
                        dyst,   &     !y-step (if linear)
                          1,    &     !z-grid type (0 - linear, 1 - levels)
                 z*1000.0d0,    &     !first z-value (if linear) or x-array (if levels)
                      0.0d0,    &     !z-step (if linear)
                          1,    &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                       2000,    &     !year   of the first field
                          1,    &     !month  of the first field
                          1,    &     !day    of the first field
                          0,    &     !hour   of the first field
                          0,    &     !minute of the first field
               86400.0*30.0,    &     !time step (in seconds)
 'Potential temperature on s-levels',    &     !title of dataset
                        'tp'   )     !variable name
   
   call fulfname(fname,pathdat,'den.dat',ierr)
   call ctl_file_write(fname,   &     !file name
                      -1e+32,   &     !value for undefined points
                    mm-mmm+1,   &     !x-dimension
                    nn-nnn+1,   &     !y-dimension
                          nz,   &     !z-dimension
                        n3da,   &     !t-dimension
                    xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                     xt(mmm),   &     !first x-value (if linear) or x-array (if levels)
                        dxst,   &     !x-step (if linear)
                    ygr_type,   &     !y-grid type (0 - linear, 1 - levels)
                     yt(nnn),   &     !first y-value (if linear) or x-array (if levels)
                        dyst,   &     !y-step (if linear)
                          1,    &     !z-grid type (0 - linear, 1 - levels)
                 z*1000.0d0,    &     !first z-value (if linear) or x-array (if levels)
                      0.0d0,    &     !z-step (if linear)
                          1,    &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                       2000,    &     !year   of the first field
                          1,    &     !month  of the first field
                          1,    &     !day    of the first field
                          0,    &     !hour   of the first field
                          0,    &     !minute of the first field
               86400.0*30.0,    &     !time step (in seconds)
 'Dinsity in-situ on s-levels',    &     !title of dataset
                        'den'   )     !variable name

 write(*,*) 'Procedure is over successfully!'
 
 call model_grid_deallocate

 stop
89 write(*,*) 'error in open: mkpot.par !'
 stop
90 write(*,*) 'error in reading file: mkpot.par !'
 stop

endprogram make_potential_temperature