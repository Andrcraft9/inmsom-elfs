!programme for interpolation from z to s-levels
program z2s_interpolation
use main_basin_pars
use basin_grid
implicit none

real(4)  tt(nx,ny,nz),    &
         uu(nx,ny,nz),    &
         ww(nx,ny,nz+1),  &
       mask(nx,ny),       &
      depth(nx,ny),       &
            zs(nz),       &
            undef

real(4) array4(nx,ny)

integer nlvl, lev1, ierr, ihtask, iztask,         &
        nt, lerr, nb, mb, k, l, m, n

real(4), allocatable::  zlvl(:),ttn(:,:,:)

character*16 name(2)
integer*2 nx1,nx2,nx3,nx4
character pathdat*32,fileh*32,file1*32,file2*32,ftemask*32,yno
character*128  fname

call non_mpi_array_boundary_definition
call model_grid_allocate

open(90,file='z2s.par',status='old',err=90)
read(90,'(a32)',err=90) pathdat !path to data(input)
read(90,'(a32)',err=90) ftemask !file with temperature mask(input)
read(90,'(a32)',err=90) fileh   !file with bottom data(input)
read(90,'(a32)',err=90) file1   !file with z level data(input)
read(90,'(a32)',err=90) file2   !file with sigma level data(output)
read(90,*,err=90) ihtask        !horizontal grid: t/u/h=1/2/3
read(90,*,err=90) iztask        !vertical grid: t/w=1/2
read(90,*      ,err=90) lev1    !task for first level
read(90,*      ,err=90) nt      !number of 3d fields
read(90,*      ,err=90) nlvl    !number of normal levels

allocate ( zlvl(nlvl),ttn(nx,ny,nlvl) )

read(90,*      ,err=90)(zlvl(k),k=1,nlvl) !levels in meter
read(90,*      ,err=90) undef !undefined value on bottom in file1
close(90)

write(*,'(2x,a,a)')'path to data(input):                    ',pathdat
write(*,'(2x,a,a)')'file with temperature mask:             ',ftemask
write(*,'(2x,a,a)')'file with bottom data(input):           ',fileh
write(*,'(2x,a,a)') 'file with data on z level (input):      ',file1
write(*,'(2x,a,a)') 'file with data on sigma level (output): ',file2
write(*,'(2x,a,i3)')' horizontal grid: t/u/h=1/2/3           ',ihtask
write(*,'(2x,a,i3)')' vertical grid: t/w=1/2                 ',iztask
write(*,'(2x,a,i3)')'1-st level = 1-st z level(1/0=y/n):     ',lev1
write(*,'(2x,a,i3)') 'number of 3d fields ',nt
write(*,'(2x,a,i3,a)')'values of ',nlvl,' levels:'
write(*,'(1x,7f8.1)')(zlvl(k),k=1,nlvl)
write(*,'(2x,a,g14.3)')'undefined value on bottom in file1 :',undef
write(*,'(2x,a)')'is it ok?'
write(*,'(2x,a)')'(input if yes - y, if no - n)'
read(*,'(a)') yno
if(yno/='y'.and.yno/='Y') stop

!setting vertical t-,w- grid levels
call vgrid

!vertical grid
if(iztask==1) then
 lerr=lerr-1
 zs(1)=0.0
 do k=1,nz
  zs(k)=sngl(z(k))
 enddo
endif
 
if(iztask==2) then
 lerr=lerr-1
 zs(1)=0.0
 do k=1,nz
  zs(k)=sngl(zw(k))
 enddo
endif

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


! task definition:
! horizontal grid
  lerr=2
  if(ihtask==1) then
   lerr=lerr-1
   mb=mmm
   nb=nnn
   
   do n = 1,ny
     do m = 1,nx
        mask(m,n) =lu(m,n)
        depth(m,n)=sngl(hhq(m,n))
     enddo
   enddo

   call fulfname(fname,pathdat,file2,ierr)
   call ctl_file_write(fname,   &     !file name
                      -1e+32,   &     !value for undefined points
                    mm-mmm+1,   &     !x-dimension
                    nn-nnn+1,   &     !y-dimension
                          nz,   &     !z-dimension
                          nt,   &     !t-dimension
                    xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                     xt(mmm),   &     !first x-value (if linear) or x-array (if levels)
                        dxst,   &     !x-step (if linear)
                    ygr_type,   &     !y-grid type (0 - linear, 1 - levels)
                     yt(nnn),   &     !first y-value (if linear) or x-array (if levels)
                        dyst,   &     !y-step (if linear)
                          1,    &     !z-grid type (0 - linear, 1 - levels)
          dble(zs)*1000.0d0,    &     !first z-value (if linear) or x-array (if levels)
                      0.0d0,    &     !z-step (if linear)
                          1,    &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                       2000,    &     !year   of the first field
                          1,    &     !month  of the first field
                          1,    &     !day    of the first field
                          0,    &     !hour   of the first field
                          0,    &     !minute of the first field
               86400.0*30.0,    &     !time step (in seconds)
         'Data on s-levels',    &     !title of dataset
                       'data'   )     !variable name

  end if

  if(ihtask==2) then
    lerr=lerr-1
    mb=mmm
    nb=nnn
    
    do n = 1,ny
     do m = 1,nx
       mask(m,n) =luu(m,n)
       depth(m,n)=sngl(hhh(m,n))
     enddo
    enddo

   call fulfname(fname,pathdat,file2,ierr)
   call ctl_file_write(fname,   &     !file name
                      -1e+32,   &     !value for undefined points
                    mm-mmm+1,   &     !x-dimension
                    nn-nnn+1,   &     !y-dimension
                          nz,   &     !z-dimension
                          nt,   &     !t-dimension
                    xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                     xu(mmm),   &     !first x-value (if linear) or x-array (if levels)
                        dxst,   &     !x-step (if linear)
                    ygr_type,   &     !y-grid type (0 - linear, 1 - levels)
                     yv(nnn),   &     !first y-value (if linear) or x-array (if levels)
                        dyst,   &     !y-step (if linear)
                          1,    &     !z-grid type (0 - linear, 1 - levels)
          dble(zs)*1000.0d0,    &     !first z-value (if linear) or x-array (if levels)
                      0.0d0,    &     !z-step (if linear)
                          1,    &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                       2000,    &     !year   of the first field
                          1,    &     !month  of the first field
                          1,    &     !day    of the first field
                          0,    &     !hour   of the first field
                          0,    &     !minute of the first field
               86400.0*30.0,    &     !time step (in seconds)
         'Data on s-levels',    &     !title of dataset
                       'data'   )     !variable name

  end if

  if(ihtask==3) then
    lerr=lerr-1
    mb=mmm-1
    nb=nnn-1
    do n = 1,ny
     do m = 1,nx
      mask(m,n) =luh(m,n)
      depth(m,n)=sngl(hhh(m,n))
     enddo
    enddo
    
   call fulfname(fname,pathdat,file2,ierr)
   call ctl_file_write(fname,   &     !file name
                      -1e+32,   &     !value for undefined points
                    mm-mmm+2,   &     !x-dimension
                    nn-nnn+2,   &     !y-dimension
                          nz,   &     !z-dimension
                          nt,   &     !t-dimension
                    xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                   xu(mmm-1),   &     !first x-value (if linear) or x-array (if levels)
                        dxst,   &     !x-step (if linear)
                    ygr_type,   &     !y-grid type (0 - linear, 1 - levels)
                   yv(nnn-1),   &     !first y-value (if linear) or x-array (if levels)
                        dyst,   &     !y-step (if linear)
                          1,    &     !z-grid type (0 - linear, 1 - levels)
          dble(zs)*1000.0d0,    &     !first z-value (if linear) or x-array (if levels)
                      0.0d0,    &     !z-step (if linear)
                          1,    &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                       2000,    &     !year   of the first field
                          1,    &     !month  of the first field
                          1,    &     !day    of the first field
                          0,    &     !hour   of the first field
                          0,    &     !minute of the first field
               86400.0*30.0,    &     !time step (in seconds)
         'Data on s-levels',    &     !title of dataset
                       'data'   )     !variable name

  end if


do l = 1,nt
 
 write(*,'(a,i4)')'  number of 3d array:',l

 !read
 call rdstd(pathdat,file1,l,ttn,mask,nx,ny,nlvl,mb,mm,nb,nn,1,nlvl,ierr)
 
 !interpolation
 call z2s(ttn,tt,depth,mask,zs,zlvl,nx,ny,nz,nlvl,lev1,undef,ierr)
 
 !write
 call wdstd(pathdat,file2,l,tt,mask,nx,ny,nz,mb,mm,nb,nn,1,nz,ierr)

enddo


deallocate ( zlvl, ttn )
call model_grid_deallocate

write(*,*) ' Procedure is over successfully!'

stop

90  write(*,*) ' Error in open or read file: s2z.par!'

endprogram z2s_interpolation
