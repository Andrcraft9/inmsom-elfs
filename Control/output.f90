subroutine create_hh(path2data,     &
                        year,       &
                       month,       &
                         day,       &
                        hour,       &
                      minute,       &
                      tstep,        &
                      calendar  )
use main_basin_pars
use mpi_parallel_tools
use basin_grid
use ocean_variables
use ocean_bc
use rec_length

implicit none
include 'locout.fi'

integer year, month, day, hour, minute, calendar, ierr
character fname*256
character*(*) path2data
real(4) array4_2d(bnd_x1:bnd_x2,bnd_y1:bnd_y2)

real(4) tstep
integer m,n,k

! area mask initialization
  call gridcon('mask_planet.txt')
! setting vertical t-,w- grid levels
  call vgrid
! define grid geographical coordinates, steps and coriolis parameters
  call basinpar


! writing HH
  ierr=0
  array4_2d = sngl( 2.94d+4 / FreeFallAcc )
  call wdstd(path2data,'LOCAL/hh.dat',1,array4_2d,lu,nx,ny,1,m1loc,m2loc,n1loc,n2loc,1,1,ierr)
  call fulfname(fname,path2data,'LOCAL/hh.dat',ierr)
  call ctl_file_write(fname,    &     !file name
                     undef,    &     !value for undefined points
                     nx_loc,   &     !x-dimension
                     ny_loc,   &     !y-dimension
                          1,   &     !z-dimension
                          1,   &     !t-dimension
                   xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                  xt(m1loc),   &     !first x-value (if linear) or x-array (if levels)
                      dxst,    &     !x-step (if linear)
                  ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                  yt(n1loc),   &     !first y-value (if linear) or x-array (if levels)
                      dyst,    &     !y-step (if linear)
                      0,       &     !z-grid type (0 - linear, 1 - levels)
                      0.0d0,   &     !first z-value (if linear) or x-array (if levels)
                      1.0d0,   &     !z-step (if linear)
                   calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                       year,   &     !year   of the first field
                      month,   &     !month  of the first field
                        day,   &     !day    of the first field
                       hour,   &     !hour   of the first field
                     minute,   &     !minute of the first field
                      tstep,   &     !time step (in seconds)
                   'Topography',   &     !title of dataset
                      'hh'   )      !variable name

end subroutine



subroutine local_output(path2data,  &
                        nrec,       &
                        year,       &
                       month,       &
                         day,       &
                        hour,       &
                      minute,       &
                      tstep,        &
                      calendar  )
use main_basin_pars
use mpi_parallel_tools
use basin_grid
use ocean_variables
use ocean_bc
use rec_length

implicit none
include 'locout.fi'

integer nrec, year, month, day, hour, minute, calendar, ierr
character fname*256
character*(*) path2data
real(4) array4_2d(bnd_x1:bnd_x2,bnd_y1:bnd_y2),       &
        array4_3d(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz)
real(4) tstep
integer m,n,k

write(*,*) 'Writing local output, record number ', nrec

if(nrec==1) then
!writing HHQ
 ierr=0
 array4_2d=sngl(hhq_rest)
 call wdstd(path2data,'LOCAL/hhq.dat',nrec,array4_2d,lu,nx,ny,1,m1loc,m2loc,n1loc,n2loc,1,1,ierr)
 call fulfname(fname,path2data,'LOCAL/hhq.dat',ierr)
 call ctl_file_write(fname,    &     !file name
                     undef,    &     !value for undefined points
                     nx_loc,   &     !x-dimension
                     ny_loc,   &     !y-dimension
                          1,   &     !z-dimension
                       nrec,   &     !t-dimension
                   xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                  xt(m1loc),   &     !first x-value (if linear) or x-array (if levels)
                      dxst,    &     !x-step (if linear)
                  ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                  yt(n1loc),   &     !first y-value (if linear) or x-array (if levels)
                      dyst,    &     !y-step (if linear)
                      0,       &     !z-grid type (0 - linear, 1 - levels)
                      0.0d0,   &     !first z-value (if linear) or x-array (if levels)
                      1.0d0,   &     !z-step (if linear)
                   calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                       year,   &     !year   of the first field
                      month,   &     !month  of the first field
                        day,   &     !day    of the first field
                       hour,   &     !hour   of the first field
                     minute,   &     !minute of the first field
                      tstep,   &     !time step (in seconds)
                   'HHQ, m',   &     !title of dataset
                      'hhq'   )      !variable name
endif

if(ssh_output>0) then
!writing SSH
 ierr=0
 array4_2d=sngl(ssh_i)
 call wdstd(path2data,'LOCAL/ssh.dat',nrec,array4_2d,lu,nx,ny,1,m1loc,m2loc,n1loc,n2loc,1,1,ierr)
 call fulfname(fname,path2data,'LOCAL/ssh.dat',ierr)
 call ctl_file_write(fname,    &     !file name
                     undef,    &     !value for undefined points
                     nx_loc,   &     !x-dimension
                     ny_loc,   &     !y-dimension
                          1,   &     !z-dimension
                       nrec,   &     !t-dimension
                   xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                  xt(m1loc),   &     !first x-value (if linear) or x-array (if levels)
                      dxst,    &     !x-step (if linear)
                  ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                  yt(n1loc),   &     !first y-value (if linear) or x-array (if levels)
                      dyst,    &     !y-step (if linear)
                      0,       &     !z-grid type (0 - linear, 1 - levels)
                      0.0d0,   &     !first z-value (if linear) or x-array (if levels)
                      1.0d0,   &     !z-step (if linear)
                   calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                       year,   &     !year   of the first field
                      month,   &     !month  of the first field
                        day,   &     !day    of the first field
                       hour,   &     !hour   of the first field
                     minute,   &     !minute of the first field
                      tstep,   &     !time step (in seconds)
                   'SSH, m',   &     !title of dataset
                      'ssh'   )      !variable name

!writing SSH error
ierr=0
array4_2d=sngl(ssh_i - ssh_err)
call wdstd(path2data,'LOCAL/ssh_err.dat',nrec,array4_2d,lu,nx,ny,1,m1loc,m2loc,n1loc,n2loc,1,1,ierr)
call fulfname(fname,path2data,'LOCAL/ssh_err.dat',ierr)
call ctl_file_write(fname,    &     !file name
                   undef,    &     !value for undefined points
                   nx_loc,   &     !x-dimension
                   ny_loc,   &     !y-dimension
                        1,   &     !z-dimension
                     nrec,   &     !t-dimension
                 xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                xt(m1loc),   &     !first x-value (if linear) or x-array (if levels)
                    dxst,    &     !x-step (if linear)
                ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                yt(n1loc),   &     !first y-value (if linear) or x-array (if levels)
                    dyst,    &     !y-step (if linear)
                    0,       &     !z-grid type (0 - linear, 1 - levels)
                    0.0d0,   &     !first z-value (if linear) or x-array (if levels)
                    1.0d0,   &     !z-step (if linear)
                 calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                     year,   &     !year   of the first field
                    month,   &     !month  of the first field
                      day,   &     !day    of the first field
                     hour,   &     !hour   of the first field
                   minute,   &     !minute of the first field
                    tstep,   &     !time step (in seconds)
                 'SSH_error, m',   &     !title of dataset
                    'ssh_err'   )      !variable name
endif

if(uv_output>0) then
!-----------------------------------------------------------------------------------------------------
!writing zonal velocity
 ierr=0
 if(grid_shift==0) then !writing on the model grid
  array4_3d=sngl(uu)
  call wdstd(path2data,'LOCAL/uu.dat',nrec,array4_3d,llu,nx,ny,nz,m1loc-1,m2loc,n1loc,n2loc,1,nz,ierr)
  call fulfname(fname,path2data,'LOCAL/vv.dat',ierr)
  call ctl_file_write(fname,    &     !file name
                      undef,    &     !value for undefined points
                    nx_loc+1,   &     !x-dimension
                      ny_loc,   &     !y-dimension
                          nz,   &     !z-dimension
                        nrec,   &     !t-dimension
                    xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                 xu(m1loc-1),   &     !first x-value (if linear) or x-array (if levels)
                       dxst,    &     !x-step (if linear)
                   ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                   yt(n1loc),   &     !first y-value (if linear) or x-array (if levels)
                       dyst,    &     !y-step (if linear)
                       1,       &     !z-grid type (0 - linear, 1 - levels)
                  z*1000.0d0,   &     !first z-value (if linear) or x-array (if levels)
                       1.0d0,   &     !z-step (if linear)
                    calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                        year,   &     !year   of the first field
                       month,   &     !month  of the first field
                         day,   &     !day    of the first field
                        hour,   &     !hour   of the first field
                      minute,   &     !minute of the first field
                       tstep,   &     !time step (in seconds
        'zonal velocity, m/s',  &     !title of dataset
                         'u'   )      !variable name

 else !writing on T-grid

 !$omp parallel do
  do n=ny_start, ny_end
   do m=nx_start, nx_end
    if(lu(m,n)>0.5) then
     do k=1,nz
     array4_3d(m,n,k)=sngl( (uu(m  ,n,k)*dxt(m  ,n)*dyh(m  ,n)*hhu(m  ,n)   &
                            +uu(m-1,n,k)*dxt(m-1,n)*dyh(m-1,n)*hhu(m-1,n))/2.0d0/hhq(m,n)/dx(m,n)/dy(m,n) )
     enddo
    endif
   enddo
  enddo
 !$omp end parallel do

  call wdstd(path2data,'LOCAL/uu.dat',nrec,array4_3d,lu,nx,ny,nz,m1loc,m2loc,n1loc,n2loc,1,nz,ierr)
  call fulfname(fname,path2data,'LOCAL/uu.dat',ierr)
  call ctl_file_write(fname,    &     !file name
                      undef,    &     !value for undefined points
                      nx_loc,   &     !x-dimension
                      ny_loc,   &     !y-dimension
                          nz,   &     !z-dimension
                        nrec,   &     !t-dimension
                    xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                   xt(m1loc),   &     !first x-value (if linear) or x-array (if levels)
                       dxst,    &     !x-step (if linear)
                   ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                   yt(n1loc),   &     !first y-value (if linear) or x-array (if levels)
                       dyst,    &     !y-step (if linear)
                       1,       &     !z-grid type (0 - linear, 1 - levels)
                  z*1000.0d0,   &     !first z-value (if linear) or x-array (if levels)
                       1.0d0,   &     !z-step (if linear)
                    calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                        year,   &     !year   of the first field
                       month,   &     !month  of the first field
                         day,   &     !day    of the first field
                        hour,   &     !hour   of the first field
                      minute,   &     !minute of the first field
                       tstep,   &     !time step (in seconds
        'zonal velocity, m/s',  &     !title of dataset
                         'u'   )      !variable name
 endif

!------------------------------------------------------------------------------------------------------
!writing meridional velocity
 ierr=0
 if(grid_shift==0) then !writing on the model grid
  array4_3d=sngl(vv)
  call wdstd(path2data,'LOCAL/vv.dat',nrec,array4_3d,llv,nx,ny,nz,m1loc,m2loc,n1loc-1,n2loc,1,nz,ierr)
  call fulfname(fname,path2data,'LOCAL/vv.dat',ierr)
  call ctl_file_write(fname,    &     !file name
                      undef,    &     !value for undefined points
                      nx_loc,   &     !x-dimension
                    ny_loc+1,   &     !y-dimension
                          nz,   &     !z-dimension
                        nrec,   &     !t-dimension
                    xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                   xt(m1loc),   &     !first x-value (if linear) or x-array (if levels)
                       dxst,    &     !x-step (if linear)
                   ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                 yv(n1loc-1),   &     !first y-value (if linear) or x-array (if levels)
                       dyst,    &     !y-step (if linear)
                       1,       &     !z-grid type (0 - linear, 1 - levels)
                  z*1000.0d0,   &     !first z-value (if linear) or x-array (if levels)
                       1.0d0,   &     !z-step (if linear)
                    calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                        year,   &     !year   of the first field
                       month,   &     !month  of the first field
                         day,   &     !day    of the first field
                        hour,   &     !hour   of the first field
                      minute,   &     !minute of the first field
                       tstep,   &     !time step (in seconds
   'meridional velocity, m/s',  &     !title of dataset
                         'v'   )      !variable name

 else !writing on T-grid

 !$omp parallel do
  do n=ny_start, ny_end
   do m=nx_start, nx_end
    if(lu(m,n)>0.5) then
     do k=1,nz
     array4_3d(m,n,k)=sngl( (vv(m,n  ,k)*dxh(m,n  )*dyt(m,n  )*hhv(m,n  )    &
                            +vv(m,n-1,k)*dxh(m,n-1)*dyt(m,n-1)*hhv(m,n-1))/2.0d0/hhq(m,n)/dx(m,n)/dy(m,n) )
     enddo
    endif
   enddo
  enddo
 !$omp end parallel do

  call wdstd(path2data,'LOCAL/vv.dat',nrec,array4_3d,lu,nx,ny,nz,m1loc,m2loc,n1loc,n2loc,1,nz,ierr)
  call fulfname(fname,path2data,'LOCAL/vv.dat',ierr)
  call ctl_file_write(fname,    &     !file name
                      undef,    &     !value for undefined points
                      nx_loc,   &     !x-dimension
                      ny_loc,   &     !y-dimension
                          nz,   &     !z-dimension
                        nrec,   &     !t-dimension
                    xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                   xt(m1loc),   &     !first x-value (if linear) or x-array (if levels)
                       dxst,    &     !x-step (if linear)
                   ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                   yt(n1loc),   &     !first y-value (if linear) or x-array (if levels)
                       dyst,    &     !y-step (if linear)
                       1,       &     !z-grid type (0 - linear, 1 - levels)
                  z*1000.0d0,   &     !first z-value (if linear) or x-array (if levels)
                       1.0d0,   &     !z-step (if linear)
                    calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                        year,   &     !year   of the first field
                       month,   &     !month  of the first field
                         day,   &     !day    of the first field
                        hour,   &     !hour   of the first field
                      minute,   &     !minute of the first field
                       tstep,   &     !time step (in seconds
   'meridional velocity, m/s',  &     !title of dataset
                         'v'   )      !variable name
 endif
endif

!-----------------------------------------------------------------------------------------------------
if(wstr_output>0) then
!writing zonal wind stress
 ierr=0
 if(grid_shift==0) then !writing on the model grid
  array4_2d=sngl(surf_stress_x)
  call wdstd(path2data,'LOCAL/tx.dat',nrec,array4_2d,llu,nx,ny,1,m1loc-1,m2loc,n1loc,n2loc,1,1,ierr)
  call fulfname(fname,path2data,'LOCAL/tx.dat',ierr)
  call ctl_file_write(fname,    &     !file name
                      undef,    &     !value for undefined points
                    nx_loc+1,   &     !x-dimension
                      ny_loc,   &     !y-dimension
                           1,   &     !z-dimension
                        nrec,   &     !t-dimension
                    xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                 xu(m1loc-1),   &     !first x-value (if linear) or x-array (if levels)
                       dxst,    &     !x-step (if linear)
                   ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                   yt(n1loc),   &     !first y-value (if linear) or x-array (if levels)
                       dyst,    &     !y-step (if linear)
                       1,       &     !z-grid type (0 - linear, 1 - levels)
                       0.0d0,   &     !first z-value (if linear) or x-array (if levels)
                       1.0d0,   &     !z-step (if linear)
                    calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                        year,   &     !year   of the first field
                       month,   &     !month  of the first field
                         day,   &     !day    of the first field
                        hour,   &     !hour   of the first field
                      minute,   &     !minute of the first field
                       tstep,   &     !time step (in seconds
 'zonal wind stress, (m/s)^2',  &     !title of dataset
                         'tx'   )      !variable name

 else !writing on T-grid

 !$omp parallel do
  do n=ny_start, ny_end
   do m=nx_start, nx_end
    if(lu(m,n)>0.5) then
     array4_2d(m,n)=sngl( (surf_stress_x(m  ,n)*dxt(m  ,n)*dyh(m  ,n)      &
                          +surf_stress_x(m-1,n)*dxt(m-1,n)*dyh(m-1,n))/2.0d0/dx(m,n)/dy(m,n) )
    endif
   enddo
  enddo
 !$omp end parallel do

  call wdstd(path2data,'LOCAL/tx.dat',nrec,array4_2d,lu,nx,ny,1,m1loc,m2loc,n1loc,n2loc,1,1,ierr)
  call fulfname(fname,path2data,'LOCAL/tx.dat',ierr)
  call ctl_file_write(fname,    &     !file name
                      undef,    &     !value for undefined points
                      nx_loc,   &     !x-dimension
                      ny_loc,   &     !y-dimension
                           1,   &     !z-dimension
                        nrec,   &     !t-dimension
                    xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                   xt(m1loc),   &     !first x-value (if linear) or x-array (if levels)
                       dxst,    &     !x-step (if linear)
                   ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                   yt(n1loc),   &     !first y-value (if linear) or x-array (if levels)
                       dyst,    &     !y-step (if linear)
                       1,       &     !z-grid type (0 - linear, 1 - levels)
                       0.0d0,   &     !first z-value (if linear) or x-array (if levels)
                       1.0d0,   &     !z-step (if linear)
                    calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                        year,   &     !year   of the first field
                       month,   &     !month  of the first field
                         day,   &     !day    of the first field
                        hour,   &     !hour   of the first field
                      minute,   &     !minute of the first field
                       tstep,   &     !time step (in seconds
 'zonal wind stress, (m/s)^2',  &     !title of dataset
                        'tx'   )      !variable name
 endif

!------------------------------------------------------------------------------------------------------
!writing meridional wind stress
 ierr=0
 if(grid_shift==0) then !writing on the model grid
  array4_2d=sngl(surf_stress_y)
  call wdstd(path2data,'LOCAL/ty.dat',nrec,array4_2d,llv,nx,ny,1,m1loc,m2loc,n1loc-1,n2loc,1,1,ierr)
  call fulfname(fname,path2data,'LOCAL/ty.dat',ierr)
  call ctl_file_write(fname,    &     !file name
                      undef,    &     !value for undefined points
                      nx_loc,   &     !x-dimension
                    ny_loc+1,   &     !y-dimension
                           1,   &     !z-dimension
                        nrec,   &     !t-dimension
                    xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                   xt(m1loc),   &     !first x-value (if linear) or x-array (if levels)
                       dxst,    &     !x-step (if linear)
                   ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                 yv(n1loc-1),   &     !first y-value (if linear) or x-array (if levels)
                       dyst,    &     !y-step (if linear)
                       1,       &     !z-grid type (0 - linear, 1 - levels)
                       0.0d0,   &     !first z-value (if linear) or x-array (if levels)
                       1.0d0,   &     !z-step (if linear)
                    calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                        year,   &     !year   of the first field
                       month,   &     !month  of the first field
                         day,   &     !day    of the first field
                        hour,   &     !hour   of the first field
                      minute,   &     !minute of the first field
                       tstep,   &     !time step (in seconds
 'meridional wind stress, (m/s)^2',  &     !title of dataset
                         'ty'   )      !variable name

 else !writing on T-grid

 !$omp parallel do
  do n=ny_start, ny_end
   do m=nx_start, nx_end
    if(lu(m,n)>0.5) then
     array4_2d(m,n)=sngl( (surf_stress_y(m,n  )*dxh(m,n  )*dyt(m,n  )  &
                          +surf_stress_y(m,n-1)*dxh(m,n-1)*dyt(m,n-1))/2.0d0/dx(m,n)/dy(m,n) )
    endif
   enddo
  enddo
 !$omp end parallel do

  call wdstd(path2data,'LOCAL/ty.dat',nrec,array4_2d,lu,nx,ny,1,m1loc,m2loc,n1loc,n2loc,1,1,ierr)
  call fulfname(fname,path2data,'LOCAL/ty.dat',ierr)
  call ctl_file_write(fname,    &     !file name
                      undef,    &     !value for undefined points
                      nx_loc,   &     !x-dimension
                      ny_loc,   &     !y-dimension
                           1,   &     !z-dimension
                        nrec,   &     !t-dimension
                    xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                   xt(m1loc),   &     !first x-value (if linear) or x-array (if levels)
                       dxst,    &     !x-step (if linear)
                   ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                   yt(n1loc),   &     !first y-value (if linear) or x-array (if levels)
                       dyst,    &     !y-step (if linear)
                       1,       &     !z-grid type (0 - linear, 1 - levels)
                       0.0d0,   &     !first z-value (if linear) or x-array (if levels)
                       1.0d0,   &     !z-step (if linear)
                    calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                        year,   &     !year   of the first field
                       month,   &     !month  of the first field
                         day,   &     !day    of the first field
                        hour,   &     !hour   of the first field
                      minute,   &     !minute of the first field
                       tstep,   &     !time step (in seconds
 'meridional wind stress, (m/s)^2',  &     !title of dataset
                         'ty'   )      !variable name
 endif
endif

!--------------------------------------------------------------------------------
if(tt_output>0) then
!writing temperature
 ierr=0
 array4_3d=sngl(tt)
 call wdstd(path2data,'LOCAL/tt.dat',nrec,array4_3d,lu,nx,ny,nz,m1loc,m2loc,n1loc,n2loc,1,nz,ierr)
 call fulfname(fname,path2data,'LOCAL/tt.dat',ierr)
 call ctl_file_write(fname,    &     !file name
                     undef,    &     !value for undefined points
                     nx_loc,   &     !x-dimension
                     ny_loc,   &     !y-dimension
                         nz,   &     !z-dimension
                       nrec,   &     !t-dimension
                   xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                  xt(m1loc),   &     !first x-value (if linear) or x-array (if levels)
                      dxst,    &     !x-step (if linear)
                  ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                  yt(n1loc),   &     !first y-value (if linear) or x-array (if levels)
                      dyst,    &     !y-step (if linear)
                      1,       &     !z-grid type (0 - linear, 1 - levels)
                 z*1000.0d0,   &     !first z-value (if linear) or x-array (if levels)
                      1.0d0,   &     !z-step (if linear)
                   calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                       year,   &     !year   of the first field
                      month,   &     !month  of the first field
                        day,   &     !day    of the first field
                       hour,   &     !hour   of the first field
                     minute,   &     !minute of the first field
                      tstep,   &     !time step (in seconds)
          'Temperature, �C',   &     !title of dataset
                       'tt'   )      !variable name
endif

!--------------------------------------------------------------------------------
if(ss_output>0) then
!writing temperature
 ierr=0
 array4_3d=sngl(ss)
 call wdstd(path2data,'LOCAL/ss.dat',nrec,array4_3d,lu,nx,ny,nz,m1loc,m2loc,n1loc,n2loc,1,nz,ierr)
 call fulfname(fname,path2data,'LOCAL/ss.dat',ierr)
 call ctl_file_write(fname,    &     !file name
                     undef,    &     !value for undefined points
                     nx_loc,   &     !x-dimension
                     ny_loc,   &     !y-dimension
                         nz,   &     !z-dimension
                       nrec,   &     !t-dimension
                   xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                  xt(m1loc),   &     !first x-value (if linear) or x-array (if levels)
                      dxst,    &     !x-step (if linear)
                  ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                  yt(n1loc),   &     !first y-value (if linear) or x-array (if levels)
                      dyst,    &     !y-step (if linear)
                      1,       &     !z-grid type (0 - linear, 1 - levels)
                 z*1000.0d0,   &     !first z-value (if linear) or x-array (if levels)
                      1.0d0,   &     !z-step (if linear)
                   calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                       year,   &     !year   of the first field
                      month,   &     !month  of the first field
                        day,   &     !day    of the first field
                       hour,   &     !hour   of the first field
                     minute,   &     !minute of the first field
                      tstep,   &     !time step (in seconds)
            'Salinity, PSU',   &     !title of dataset
                       'ss'   )      !variable name
endif

!--------------------------------------------------------------------------------------------------
!writing surface fluxes
if(sfl_output>0) then
 !$omp parallel do
  do n=ny_start, ny_end
   do m=nx_start, nx_end
    if(lu(m,n)>0.5) then
     array4_3d(m,n,1)=sngl(tatm(m,n))
     array4_3d(m,n,2)=sngl(qatm(m,n))
     array4_3d(m,n,3)=sngl( lwr(m,n))
     array4_3d(m,n,4)=sngl( swr(m,n))
     array4_3d(m,n,5)=sngl(slpr(m,n))
     array4_3d(m,n,6)=sngl(sensheat(m,n))
     array4_3d(m,n,7)=sngl( latheat(m,n))
     array4_3d(m,n,8)=sngl(   lw_bal(m,n))
     array4_3d(m,n,9)=sngl(   sw_bal(m,n))
     array4_3d(m,n,10)=sngl(  hf_tot(m,n))
     array4_3d(m,n,11)=sngl(rain(m,n))
     array4_3d(m,n,12)=sngl(snow(m,n))
     array4_3d(m,n,13)=runoff(m,n)
     array4_3d(m,n,14)=sngl(wf_tot(m,n))
     array4_3d(m,n,15)=sngl(tflux_surf(m,n))
     array4_3d(m,n,16)=sngl(sflux_surf(m,n))
    endif
   enddo
  enddo
 !$omp end parallel do
ierr=0
 call wdstd(path2data,'LOCAL/sfl.dat',nrec,array4_3d,lu,nx,ny,16,m1loc,m2loc,n1loc,n2loc,1,16,ierr)
 call fulfname(fname,path2data,'LOCAL/sfl.dat',ierr)
 call ctl_file_write(fname,    &     !file name
                     undef,    &     !value for undefined points
                     nx_loc,   &     !x-dimension
                     ny_loc,   &     !y-dimension
                         16,   &     !z-dimension
                       nrec,   &     !t-dimension
                   xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                  xt(m1loc),   &     !first x-value (if linear) or x-array (if levels)
                      dxst,    &     !x-step (if linear)
                  ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                  yt(n1loc),   &     !first y-value (if linear) or x-array (if levels)
                      dyst,    &     !y-step (if linear)
                      0,       &     !z-grid type (0 - linear, 1 - levels)
                      1.0d0,   &     !first z-value (if linear) or x-array (if levels)
                      1.0d0,   &     !z-step (if linear)
                   calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                       year,   &     !year   of the first field
                      month,   &     !month  of the first field
                        day,   &     !day    of the first field
                       hour,   &     !hour   of the first field
                     minute,   &     !minute of the first field
                      tstep,   &     !time step (in seconds)
'1-TA, 2-QA, 3-LWdw, 4-SWdw, 5-SLP, 6-SensBal, 7-LatBal, 8-LWbal, 9-SWbal, 10-HFtot, 11-rain, 12-snow, 13-roff, 14-WFtot, 15-Tflux, 16-Sflux',   &     !title of dataset
                      'sfl'   )      !variable name

endif

!-----------------------------------------------------------------------------------------------
!write wind speed
if(wind_output>0) then

  array4_2d=sngl(uwnd)

  call wdstd(path2data,'LOCAL/uwnd.dat',nrec,array4_2d,lu,nx,ny,1,m1loc,m2loc,n1loc,n2loc,1,1,ierr)
  call fulfname(fname,path2data,'LOCAL/uwnd.dat',ierr)
  call ctl_file_write(fname,    &     !file name
                      undef,    &     !value for undefined points
                      nx_loc,   &     !x-dimension
                      ny_loc,   &     !y-dimension
                           1,   &     !z-dimension
                        nrec,   &     !t-dimension
                    xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                   xt(m1loc),   &     !first x-value (if linear) or x-array (if levels)
                       dxst,    &     !x-step (if linear)
                   ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                   yt(n1loc),   &     !first y-value (if linear) or x-array (if levels)
                       dyst,    &     !y-step (if linear)
                       1,       &     !z-grid type (0 - linear, 1 - levels)
                       0.0d0,   &     !first z-value (if linear) or x-array (if levels)
                       1.0d0,   &     !z-step (if linear)
                    calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                        year,   &     !year   of the first field
                       month,   &     !month  of the first field
                         day,   &     !day    of the first field
                        hour,   &     !hour   of the first field
                      minute,   &     !minute of the first field
                       tstep,   &     !time step (in seconds
      'zonal wind speed, m/s',  &     !title of dataset
                         'u'   )      !variable name

  array4_2d=sngl(vwnd)

  call wdstd(path2data,'LOCAL/vwnd.dat',nrec,array4_2d,lu,nx,ny,1,m1loc,m2loc,n1loc,n2loc,1,1,ierr)
  call fulfname(fname,path2data,'LOCAL/vwnd.dat',ierr)
  call ctl_file_write(fname,    &     !file name
                      undef,    &     !value for undefined points
                      nx_loc,   &     !x-dimension
                      ny_loc,   &     !y-dimension
                           1,   &     !z-dimension
                        nrec,   &     !t-dimension
                    xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                   xt(m1loc),   &     !first x-value (if linear) or x-array (if levels)
                       dxst,    &     !x-step (if linear)
                   ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                   yt(n1loc),   &     !first y-value (if linear) or x-array (if levels)
                       dyst,    &     !y-step (if linear)
                       1,       &     !z-grid type (0 - linear, 1 - levels)
                       0.0d0,   &     !first z-value (if linear) or x-array (if levels)
                       1.0d0,   &     !z-step (if linear)
                    calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                        year,   &     !year   of the first field
                       month,   &     !month  of the first field
                         day,   &     !day    of the first field
                        hour,   &     !hour   of the first field
                      minute,   &     !minute of the first field
                       tstep,   &     !time step (in seconds
 'meridional wind speed, m/s',  &     !title of dataset
                         'v'   )      !variable name
endif

!--------------------------------------------------------------------------------
if(amuv_output>0) then
!writing lateral viscocity
 ierr=0
 array4_3d=sngl(amuv)
 call wdstd(path2data,'LOCAL/amuv.dat',nrec,array4_3d,lu,nx,ny,nz,m1loc,m2loc,n1loc,n2loc,1,nz,ierr)
 call fulfname(fname,path2data,'LOCAL/amuv.dat',ierr)
 call ctl_file_write(fname,    &     !file name
                     undef,    &     !value for undefined points
                     nx_loc,   &     !x-dimension
                     ny_loc,   &     !y-dimension
                         nz,   &     !z-dimension
                       nrec,   &     !t-dimension
                   xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                  xt(m1loc),   &     !first x-value (if linear) or x-array (if levels)
                      dxst,    &     !x-step (if linear)
                  ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                  yt(n1loc),   &     !first y-value (if linear) or x-array (if levels)
                      dyst,    &     !y-step (if linear)
                      1,       &     !z-grid type (0 - linear, 1 - levels)
                 z*1000.0d0,   &     !first z-value (if linear) or x-array (if levels)
                      1.0d0,   &     !z-step (if linear)
                   calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                       year,   &     !year   of the first field
                      month,   &     !month  of the first field
                        day,   &     !day    of the first field
                       hour,   &     !hour   of the first field
                     minute,   &     !minute of the first field
                      tstep,   &     !time step (in seconds)
 'Lateral viscosity, m^2/s',   &     !title of dataset
                       'mu'   )      !variable name
endif

!--------------------------------------------------------------------------------
if(anzu_output>0) then
!writing vertical viscosity
 ierr=0
 array4_3d(:,:,1:nz)=sngl(anzu(:,:,1:nz))
 call wdstd(path2data,'LOCAL/anzu.dat',nrec,array4_3d,lu,nx,ny,nz,m1loc,m2loc,n1loc,n2loc,1,nz,ierr)
 call fulfname(fname,path2data,'LOCAL/anzu.dat',ierr)
 call ctl_file_write(fname,    &     !file name
                     undef,    &     !value for undefined points
                     nx_loc,   &     !x-dimension
                     ny_loc,   &     !y-dimension
                         nz,   &     !z-dimension
                       nrec,   &     !t-dimension
                   xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                  xt(m1loc),   &     !first x-value (if linear) or x-array (if levels)
                      dxst,    &     !x-step (if linear)
                  ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                  yt(n1loc),   &     !first y-value (if linear) or x-array (if levels)
                      dyst,    &     !y-step (if linear)
                      1,       &     !z-grid type (0 - linear, 1 - levels)
                zw*1000.0d0,   &     !first z-value (if linear) or x-array (if levels)
                      1.0d0,   &     !z-step (if linear)
                   calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                       year,   &     !year   of the first field
                      month,   &     !month  of the first field
                        day,   &     !day    of the first field
                       hour,   &     !hour   of the first field
                     minute,   &     !minute of the first field
                      tstep,   &     !time step (in seconds)
'vertical viscosity, m^2/s',   &     !title of dataset
                       'nu'   )      !variable name
endif

!--------------------------------------------------------------------------------
if(ww_output>0) then
!writing vertical velocity
 ierr=0
 array4_3d(:,:,1:nz)=sngl(ww(:,:,1:nz))
 call wdstd(path2data,'LOCAL/ww.dat',nrec,array4_3d,lu,nx,ny,nz,m1loc,m2loc,n1loc,n2loc,1,nz,ierr)
 call fulfname(fname,path2data,'LOCAL/ww.dat',ierr)
 call ctl_file_write(fname,    &     !file name
                     undef,    &     !value for undefined points
                     nx_loc,   &     !x-dimension
                     ny_loc,   &     !y-dimension
                         nz,   &     !z-dimension
                       nrec,   &     !t-dimension
                   xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                  xt(m1loc),   &     !first x-value (if linear) or x-array (if levels)
                      dxst,    &     !x-step (if linear)
                  ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                  yt(n1loc),   &     !first y-value (if linear) or x-array (if levels)
                      dyst,    &     !y-step (if linear)
                      1,       &     !z-grid type (0 - linear, 1 - levels)
                zw*1000.0d0,   &     !first z-value (if linear) or x-array (if levels)
                      1.0d0,   &     !z-step (if linear)
                   calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                       year,   &     !year   of the first field
                      month,   &     !month  of the first field
                        day,   &     !day    of the first field
                       hour,   &     !hour   of the first field
                     minute,   &     !minute of the first field
                      tstep,   &     !time step (in seconds)
   'vertical velocity, m/s',   &     !title of dataset
                        'w'   )      !variable name
endif

endsubroutine local_output
!===================================================================================================================
subroutine ocpwrite(path2data,  &
                        year,       &
                       month,       &
                         day,       &
                        hour,       &
                      minute,       &
                      tstep,        &
                      calendar,     &
                      nstep_btr     )
use main_basin_pars
use mpi_parallel_tools
use basin_grid
use ocean_variables
use rec_length

implicit none

integer year, month, day, hour, minute, calendar, ierr
integer nstep_btr
character fname*256
character*(*) path2data
real(4) array4_2d(bnd_x1:bnd_x2,bnd_y1:bnd_y2),       &
        array4_3d(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz)
real(4) tstep

ierr=0
call wdstd8(path2data,'cptt8.dat',1,tt ,lu,nx,ny,nz,mmm,mm,nnn,nn,1,nz,ierr)
call wdstd8(path2data,'cptt8.dat',2,ttp,lu,nx,ny,nz,mmm,mm,nnn,nn,1,nz,ierr)

ierr=0
array4_3d=sngl(tt)
call wdstd(path2data,'cptt.dat',1,array4_3d,lu,nx,ny,nz,mmm,mm,nnn,nn,1,nz,ierr)
array4_3d=sngl(ttp)
call wdstd(path2data,'cptt.dat',2,array4_3d,lu,nx,ny,nz,mmm,mm,nnn,nn,1,nz,ierr)
call fulfname(fname,path2data,'cptt.dat',ierr)
call ctl_file_write(fname,    &     !file name
                    undef,    &     !value for undefined points
                   nx_glob,   &     !x-dimension
                   ny_glob,   &     !y-dimension
                        nz,   &     !z-dimension
                         2,   &     !t-dimension
                  xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                   xt(mmm),   &     !first x-value (if linear) or x-array (if levels)
                      dxst,   &     !x-step (if linear)
                  ygr_type,   &     !y-grid type (0 - linear, 1 - levels)
                   yt(nnn),   &     !first y-value (if linear) or x-array (if levels)
                      dyst,   &     !y-step (if linear)
                         1,   &     !z-grid type (0 - linear, 1 - levels)
                z*1000.0d0,   &     !first z-value (if linear) or x-array (if levels)
                     1.0d0,   &     !z-step (if linear)
                  calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                      year,   &     !year   of the first field
                     month,   &     !month  of the first field
                       day,   &     !day    of the first field
                      hour,   &     !hour   of the first field
                    minute,   &     !minute of the first field
                     tstep,   &     !time step (in seconds
'Potential temperature, C',   &     !title of dataset
                       'tt'   )     !variable name

!-------------------------------------------------------------------------------------
ierr=0
call wdstd8(path2data,'cpss8.dat',1,ss ,lu,nx,ny,nz,mmm,mm,nnn,nn,1,nz,ierr)
call wdstd8(path2data,'cpss8.dat',2,ssp,lu,nx,ny,nz,mmm,mm,nnn,nn,1,nz,ierr)

ierr=0
array4_3d=sngl(ss)
call wdstd(path2data,'cpss.dat',1,array4_3d,lu,nx,ny,nz,mmm,mm,nnn,nn,1,nz,ierr)
array4_3d=sngl(ssp)
call wdstd(path2data,'cpss.dat',2,array4_3d,lu,nx,ny,nz,mmm,mm,nnn,nn,1,nz,ierr)

call fulfname(fname,path2data,'cpss.dat',ierr)
call ctl_file_write(fname,    &     !file name
                    undef,    &     !value for undefined points
                   nx_glob,   &     !x-dimension
                   ny_glob,   &     !y-dimension
                        nz,   &     !z-dimension
                         2,   &     !t-dimension
                  xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                   xt(mmm),   &     !first x-value (if linear) or x-array (if levels)
                      dxst,   &     !x-step (if linear)
                  ygr_type,   &     !y-grid type (0 - linear, 1 - levels)
                   yt(nnn),   &     !first y-value (if linear) or x-array (if levels)
                      dyst,   &     !y-step (if linear)
                         1,   &     !z-grid type (0 - linear, 1 - levels)
                z*1000.0d0,   &     !first z-value (if linear) or x-array (if levels)
                     1.0d0,   &     !z-step (if linear)
                  calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                      year,   &     !year   of the first field
                     month,   &     !month  of the first field
                       day,   &     !day    of the first field
                      hour,   &     !hour   of the first field
                    minute,   &     !minute of the first field
                     tstep,   &     !time step (in seconds
           'Salinity, PSU',   &     !title of dataset
                       'ss'   )     !variable name

!-----------------------------------------------------------------------------------
ierr=0
call wdstd8(path2data,'cpuu8.dat',1,uu ,llu,nx,ny,nz,mmm-1,mm,nnn,nn,1,nz,ierr)
call wdstd8(path2data,'cpuu8.dat',2,uup,llu,nx,ny,nz,mmm-1,mm,nnn,nn,1,nz,ierr)

ierr=0
array4_3d=sngl(uu)
call wdstd(path2data,'cpuu.dat',1,array4_3d,llu,nx,ny,nz,mmm-1,mm,nnn,nn,1,nz,ierr)
array4_3d=sngl(uup)
call wdstd(path2data,'cpuu.dat',2,array4_3d,llu,nx,ny,nz,mmm-1,mm,nnn,nn,1,nz,ierr)

call fulfname(fname,path2data,'cpuu.dat',ierr)
call ctl_file_write(fname,    &     !file name
                    undef,    &     !value for undefined points
                   nx_glob+1, &     !x-dimension
                   ny_glob,   &     !y-dimension
                        nz,   &     !z-dimension
                         2,   &     !t-dimension
                  xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                 xu(mmm-1),   &     !first x-value (if linear) or x-array (if levels)
                      dxst,   &     !x-step (if linear)
                  ygr_type,   &     !y-grid type (0 - linear, 1 - levels)
                   yt(nnn),   &     !first y-value (if linear) or x-array (if levels)
                      dyst,   &     !y-step (if linear)
                         1,   &     !z-grid type (0 - linear, 1 - levels)
                z*1000.0d0,   &     !first z-value (if linear) or x-array (if levels)
                     1.0d0,   &     !z-step (if linear)
                  calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                      year,   &     !year   of the first field
                     month,   &     !month  of the first field
                       day,   &     !day    of the first field
                      hour,   &     !hour   of the first field
                    minute,   &     !minute of the first field
                     tstep,   &     !time step (in seconds
     'Zonal velocity, m/s',   &     !title of dataset
                       'u'   )     !variable name

!----------------------------------------------------------------------------------
ierr=0
call wdstd8(path2data,'cpvv8.dat',1,vv ,llv,nx,ny,nz,mmm,mm,nnn-1,nn,1,nz,ierr)
call wdstd8(path2data,'cpvv8.dat',2,vvp,llv,nx,ny,nz,mmm,mm,nnn-1,nn,1,nz,ierr)

ierr=0
array4_3d=sngl(vv)
call wdstd(path2data,'cpvv.dat',1,array4_3d,llv,nx,ny,nz,mmm,mm,nnn-1,nn,1,nz,ierr)
array4_3d=sngl(vvp)
call wdstd(path2data,'cpvv.dat',2,array4_3d,llv,nx,ny,nz,mmm,mm,nnn-1,nn,1,nz,ierr)

call fulfname(fname,path2data,'cpvv.dat',ierr)
call ctl_file_write(fname,    &     !file name
                    undef,    &     !value for undefined points
                   nx_glob,   &     !x-dimension
                   ny_glob+1, &     !y-dimension
                        nz,   &     !z-dimension
                         2,   &     !t-dimension
                  xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                   xt(mmm),   &     !first x-value (if linear) or x-array (if levels)
                      dxst,   &     !x-step (if linear)
                  ygr_type,   &     !y-grid type (0 - linear, 1 - levels)
                 yv(nnn-1),   &     !first y-value (if linear) or x-array (if levels)
                      dyst,   &     !y-step (if linear)
                         1,   &     !z-grid type (0 - linear, 1 - levels)
                z*1000.0d0,   &     !first z-value (if linear) or x-array (if levels)
                     1.0d0,   &     !z-step (if linear)
                  calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                      year,   &     !year   of the first field
                     month,   &     !month  of the first field
                       day,   &     !day    of the first field
                      hour,   &     !hour   of the first field
                    minute,   &     !minute of the first field
                     tstep,   &     !time step (in seconds
'Meridional velocity, m/s',   &     !title of dataset
                       'v'   )     !variable name

!-----------------------------------------------------------------------------------------
ierr=0
call wdstd8(path2data,'cpsshi8.dat',1, ssh_i, lu,nx,ny,1,mmm,mm,nnn,nn,1,1,ierr)
call wdstd8(path2data,'cpsshi8.dat',2,sshp_i, lu,nx,ny,1,mmm,mm,nnn,nn,1,1,ierr)

ierr=0
array4_2d=sngl(ssh_i)
call wdstd(path2data,'cpsshi.dat',1,array4_2d,lu,nx,ny,1,mmm,mm,nnn,nn,1,1,ierr)
array4_2d=sngl(sshp_i)
call wdstd(path2data,'cpsshi.dat',2,array4_2d,lu,nx,ny,1,mmm,mm,nnn,nn,1,1,ierr)

call fulfname(fname,path2data,'cpsshi.dat',ierr)
call ctl_file_write(fname,    &     !file name
                    undef,    &     !value for undefined points
                   nx_glob,   &     !x-dimension
                   ny_glob,   &     !y-dimension
                         1,   &     !z-dimension
                         2,   &     !t-dimension
                  xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                   xt(mmm),   &     !first x-value (if linear) or x-array (if levels)
                      dxst,   &     !x-step (if linear)
                  ygr_type,   &     !y-grid type (0 - linear, 1 - levels)
                   yt(nnn),   &     !first y-value (if linear) or x-array (if levels)
                      dyst,   &     !y-step (if linear)
                         0,   &     !z-grid type (0 - linear, 1 - levels)
                     0.0d0,   &     !first z-value (if linear) or x-array (if levels)
                     1.0d0,   &     !z-step (if linear)
                  calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                      year,   &     !year   of the first field
                     month,   &     !month  of the first field
                       day,   &     !day    of the first field
                      hour,   &     !hour   of the first field
                    minute,   &     !minute of the first field
                     tstep,   &     !time step (in seconds
         'Internal SSH, m',   &     !title of dataset
                      'ssh'   )     !variable name

!-----------------------------------------------------------------------------------------
ierr=0
call wdstd8(path2data,'cpsshe8.dat',1, ssh_e, lu,nx,ny,1,mmm,mm,nnn,nn,1,1,ierr)
call wdstd8(path2data,'cpsshe8.dat',2,sshp_e, lu,nx,ny,1,mmm,mm,nnn,nn,1,1,ierr)

ierr=0
array4_2d=sngl(ssh_e)
call wdstd(path2data,'cpsshe.dat',1,array4_2d,lu,nx,ny,1,mmm,mm,nnn,nn,1,1,ierr)
array4_2d=sngl(sshp_e)
call wdstd(path2data,'cpsshe.dat',2,array4_2d,lu,nx,ny,1,mmm,mm,nnn,nn,1,1,ierr)

call fulfname(fname,path2data,'cpsshe.dat',ierr)
call ctl_file_write(fname,    &     !file name
                    undef,    &     !value for undefined points
                   nx_glob,   &     !x-dimension
                   ny_glob,   &     !y-dimension
                         1,   &     !z-dimension
                         2,   &     !t-dimension
                  xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                   xt(mmm),   &     !first x-value (if linear) or x-array (if levels)
                      dxst,   &     !x-step (if linear)
                  ygr_type,   &     !y-grid type (0 - linear, 1 - levels)
                   yt(nnn),   &     !first y-value (if linear) or x-array (if levels)
                      dyst,   &     !y-step (if linear)
                         0,   &     !z-grid type (0 - linear, 1 - levels)
                     0.0d0,   &     !first z-value (if linear) or x-array (if levels)
                     1.0d0,   &     !z-step (if linear)
                  calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                      year,   &     !year   of the first field
                     month,   &     !month  of the first field
                       day,   &     !day    of the first field
                      hour,   &     !hour   of the first field
                    minute,   &     !minute of the first field
                     tstep,   &     !time step (in seconds)
         'External SSH, m',   &     !title of dataset
                      'ssh'   )     !variable name

!----------------------------------------------------------------------------------------
ierr=0
call wdstd8(path2data,'cpube8.dat',1, ubrtr_e,  llu,nx,ny,1,mmm-1,mm,nnn,nn,1,1,ierr)
call wdstd8(path2data,'cpube8.dat',2, ubrtrp_e, llu,nx,ny,1,mmm-1,mm,nnn,nn,1,1,ierr)

ierr=0
array4_2d=sngl(ubrtr_e)
call wdstd(path2data,'cpubi.dat',1,array4_2d,llu,nx,ny,1,mmm-1,mm,nnn,nn,1,1,ierr)
array4_2d=sngl(ubrtrp_e)
call wdstd(path2data,'cpubi.dat',2,array4_2d,llu,nx,ny,1,mmm-1,mm,nnn,nn,1,1,ierr)

call fulfname(fname,path2data,'cpube.dat',ierr)
call ctl_file_write(fname,    &     !file name
                    undef,    &     !value for undefined points
                   nx_glob+1, &     !x-dimension
                   ny_glob,   &     !y-dimension
                         1,   &     !z-dimension
                         2,   &     !t-dimension
                  xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                   xu(mmm-1), &     !first x-value (if linear) or x-array (if levels)
                      dxst,   &     !x-step (if linear)
                  ygr_type,   &     !y-grid type (0 - linear, 1 - levels)
                   yt(nnn),   &     !first y-value (if linear) or x-array (if levels)
                      dyst,   &     !y-step (if linear)
                         0,   &     !z-grid type (0 - linear, 1 - levels)
                     0.0d0,   &     !first z-value (if linear) or x-array (if levels)
                     1.0d0,   &     !z-step (if linear)
                  calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                      year,   &     !year   of the first field
                     month,   &     !month  of the first field
                       day,   &     !day    of the first field
                      hour,   &     !hour   of the first field
                    minute,   &     !minute of the first field
                     tstep,   &     !time step (in seconds)
        'External UB, m/s',   &     !title of dataset
                      'ub'   )      !variable name

!--------------------------------------------------------------------------------------------------
ierr=0
call wdstd8(path2data,'cpvbe8.dat',1,  vbrtr_e, llv,nx,ny,1,mmm,mm,nnn-1,nn,1,1,ierr)
call wdstd8(path2data,'cpvbe8.dat',2, vbrtrp_e, llv,nx,ny,1,mmm,mm,nnn-1,nn,1,1,ierr)

ierr=0
array4_2d=sngl(vbrtr_e)
call wdstd(path2data,'cpvbe.dat',1,array4_2d,llv,nx,ny,1,mmm,mm,nnn-1,nn,1,1,ierr)
array4_2d=sngl(vbrtrp_e)
call wdstd(path2data,'cpvbe.dat',2,array4_2d,llv,nx,ny,1,mmm,mm,nnn-1,nn,1,1,ierr)

call fulfname(fname,path2data,'cpvbe.dat',ierr)
call ctl_file_write(fname,    &     !file name
                    undef,    &     !value for undefined points
                   nx_glob,   &     !x-dimension
                   ny_glob+1, &     !y-dimension
                         1,   &     !z-dimension
                         2,   &     !t-dimension
                  xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                   xt(mmm),   &     !first x-value (if linear) or x-array (if levels)
                      dxst,   &     !x-step (if linear)
                  ygr_type,   &     !y-grid type (0 - linear, 1 - levels)
                   yv(nnn-1), &     !first y-value (if linear) or x-array (if levels)
                      dyst,   &     !y-step (if linear)
                         0,   &     !z-grid type (0 - linear, 1 - levels)
                     0.0d0,   &     !first z-value (if linear) or x-array (if levels)
                     1.0d0,   &     !z-step (if linear)
                  calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                      year,   &     !year   of the first field
                     month,   &     !month  of the first field
                       day,   &     !day    of the first field
                      hour,   &     !hour   of the first field
                    minute,   &     !minute of the first field
                     tstep,   &     !time step (in seconds)
        'External VB, m/s',   &     !title of dataset
                      'vb'   )      !variable name

endsubroutine ocpwrite

!============================================================================
subroutine data_calc
use main_basin_pars
use mpi_parallel_tools
use basin_grid
use ocean_variables
use ocean_bc
implicit none
integer m,n,k

 meancalc=meancalc+1

!$omp parallel do
 do n=ny_start-1,ny_end+1
  do m=nx_start-1,nx_end+1
    do k=1,nz
      tt_calc(m,n,k)=  tt_calc(m,n,k)+tt(m,n,k)
      ss_calc(m,n,k)=  ss_calc(m,n,k)+ss(m,n,k)
      uu_calc(m,n,k)=  uu_calc(m,n,k)+uu(m,n,k)
      vv_calc(m,n,k)=  vv_calc(m,n,k)+vv(m,n,k)
    RHSt_calc(m,n)  =RHSt_calc(m,n) + RHS_tem(m,n,k)*dz(k)/dx(m,n)/dy(m,n)
    RHSs_calc(m,n)  =RHSs_calc(m,n) + RHS_sal(m,n,k)*dz(k)/dx(m,n)/dy(m,n)
    enddo
      ssh_calc(m,n)=  ssh_calc(m,n)+ ssh_i(m,n)
     uwnd_calc(m,n)= uwnd_calc(m,n)+uwnd(m,n)
     vwnd_calc(m,n)= vwnd_calc(m,n)+vwnd(m,n)
      txo_calc(m,n)=  txo_calc(m,n)+surf_stress_x(m,n)
      tyo_calc(m,n)=  tyo_calc(m,n)+surf_stress_y(m,n)
      sfl_calc(m,n,1) = sfl_calc(m,n,1) +  tatm(m,n)
      sfl_calc(m,n,2) = sfl_calc(m,n,2) +  qatm(m,n)
      sfl_calc(m,n,3) = sfl_calc(m,n,3) +  lwr(m,n)
      sfl_calc(m,n,4) = sfl_calc(m,n,4) +  swr(m,n)
      sfl_calc(m,n,5) = sfl_calc(m,n,5) + slpr(m,n)
      sfl_calc(m,n,6) = sfl_calc(m,n,6) + sensheat(m,n)
      sfl_calc(m,n,7) = sfl_calc(m,n,7) + latheat(m,n)
      sfl_calc(m,n,8) = sfl_calc(m,n,8) +  lw_bal(m,n)
      sfl_calc(m,n,9) = sfl_calc(m,n,9) +  sw_bal(m,n)
      sfl_calc(m,n,10)= sfl_calc(m,n,10)+  hf_tot(m,n)
      sfl_calc(m,n,11)= sfl_calc(m,n,11)+    rain(m,n)
      sfl_calc(m,n,12)= sfl_calc(m,n,12)+    snow(m,n)
      sfl_calc(m,n,13)= sfl_calc(m,n,13)+  runoff(m,n)
      sfl_calc(m,n,14)= sfl_calc(m,n,14)+  wf_tot(m,n)
      sfl_calc(m,n,15)= sfl_calc(m,n,15)+  tflux_surf(m,n)
      sfl_calc(m,n,16)= sfl_calc(m,n,16)+  sflux_surf(m,n)
  enddo
 enddo
!$omp end parallel do

endsubroutine data_calc

!==============================================================================================
subroutine global_output(path2data,  &
                         nrec,       &
                         year,       &
                        month,       &
                          day,       &
                         hour,       &
                       minute,       &
                       tstep,        &
                       calendar  )
use main_basin_pars
use mpi_parallel_tools
use basin_grid
use ocean_variables
use ocean_bc
use rec_length

implicit none
include 'globout.fi'

integer nrec, year, month, day, hour, minute, calendar, ierr
character fname*256
character*(*) path2data
real(4) array4_2d(bnd_x1:bnd_x2,bnd_y1:bnd_y2),       &
        array4_3d(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz)
real(4) tstep
integer m,n,k

write(*,*) 'Writing global output, record number ', nrec

if(nrec==1) then
!writing HHQ
 ierr=0
 array4_2d=sngl(hhq_rest)
 call wdstd(path2data,'GLOBAL/hhq.dat',nrec,array4_2d,lu,nx,ny,1,mmm,mm,nnn,nn,1,1,ierr)
 call fulfname(fname,path2data,'GLOBAL/hhq.dat',ierr)
 call ctl_file_write(fname,    &     !file name
                     undef,    &     !value for undefined points
                    nx_glob,   &     !x-dimension
                    ny_glob,   &     !y-dimension
                          1,   &     !z-dimension
                       nrec,   &     !t-dimension
                   xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                    xt(mmm),   &     !first x-value (if linear) or x-array (if levels)
                      dxst,    &     !x-step (if linear)
                  ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                    yt(nnn),   &     !first y-value (if linear) or x-array (if levels)
                      dyst,    &     !y-step (if linear)
                      0,       &     !z-grid type (0 - linear, 1 - levels)
                      0.0d0,   &     !first z-value (if linear) or x-array (if levels)
                      1.0d0,   &     !z-step (if linear)
                   calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                       year,   &     !year   of the first field
                      month,   &     !month  of the first field
                        day,   &     !day    of the first field
                       hour,   &     !hour   of the first field
                     minute,   &     !minute of the first field
                      tstep,   &     !time step (in seconds)
                   'HHQ, m',   &     !title of dataset
                      'hhq'   )      !variable name
endif

if(ssh_output>0) then
!writing SSH
 ierr=0
 array4_2d=sngl(ssh_calc/dfloat(meancalc))
 call wdstd(path2data,'GLOBAL/ssh.dat',nrec,array4_2d,lu,nx,ny,1,mmm,mm,nnn,nn,1,1,ierr)
 call fulfname(fname,path2data,'GLOBAL/ssh.dat',ierr)
 call ctl_file_write(fname,    &     !file name
                     undef,    &     !value for undefined points
                    nx_glob,   &     !x-dimension
                    ny_glob,   &     !y-dimension
                          1,   &     !z-dimension
                       nrec,   &     !t-dimension
                   xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                    xt(mmm),   &     !first x-value (if linear) or x-array (if levels)
                      dxst,    &     !x-step (if linear)
                  ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                    yt(nnn),   &     !first y-value (if linear) or x-array (if levels)
                      dyst,    &     !y-step (if linear)
                      0,       &     !z-grid type (0 - linear, 1 - levels)
                      0.0d0,   &     !first z-value (if linear) or x-array (if levels)
                      1.0d0,   &     !z-step (if linear)
                   calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                       year,   &     !year   of the first field
                      month,   &     !month  of the first field
                        day,   &     !day    of the first field
                       hour,   &     !hour   of the first field
                     minute,   &     !minute of the first field
                      tstep,   &     !time step (in seconds)
                   'SSH, m',   &     !title of dataset
                      'ssh'   )      !variable name
endif
ssh_calc=0.0d0

if(uv_output>0) then
!-----------------------------------------------------------------------------------------------------
!writing zonal velocity
 ierr=0
 if(grid_shift==0) then !writing on the model grid
  array4_3d=sngl(uu_calc/dfloat(meancalc))
  call wdstd(path2data,'GLOBAL/uu.dat',nrec,array4_3d,llu,nx,ny,nz,mmm-1,mm,nnn,nn,1,nz,ierr)
  call fulfname(fname,path2data,'GLOBAL/vv.dat',ierr)
  call ctl_file_write(fname,    &     !file name
                      undef,    &     !value for undefined points
                    nx_glob+1,  &     !x-dimension
                    ny_glob,    &     !y-dimension
                          nz,   &     !z-dimension
                        nrec,   &     !t-dimension
                    xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                   xu(mmm-1),   &     !first x-value (if linear) or x-array (if levels)
                       dxst,    &     !x-step (if linear)
                   ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                     yt(nnn),   &     !first y-value (if linear) or x-array (if levels)
                       dyst,    &     !y-step (if linear)
                       1,       &     !z-grid type (0 - linear, 1 - levels)
                  z*1000.0d0,   &     !first z-value (if linear) or x-array (if levels)
                       1.0d0,   &     !z-step (if linear)
                    calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                        year,   &     !year   of the first field
                       month,   &     !month  of the first field
                         day,   &     !day    of the first field
                        hour,   &     !hour   of the first field
                      minute,   &     !minute of the first field
                       tstep,   &     !time step (in seconds
        'zonal velocity, m/s',  &     !title of dataset
                         'u'   )      !variable name

 else !writing on T-grid

 !$omp parallel do
  do n=ny_start, ny_end
   do m=nx_start, nx_end
    if(lu(m,n)>0.5) then
     do k=1,nz
     array4_3d(m,n,k)=sngl( (uu_calc(m  ,n,k)*dxt(m  ,n)*dyh(m  ,n)*hhu(m  ,n)   &
                            +uu_calc(m-1,n,k)*dxt(m-1,n)*dyh(m-1,n)*hhu(m-1,n))/2.0d0/hhq(m,n)/dx(m,n)/dy(m,n)/dfloat(meancalc) )
     enddo
    endif
   enddo
  enddo
 !$omp end parallel do

  call wdstd(path2data,'GLOBAL/uu.dat',nrec,array4_3d,lu,nx,ny,nz,mmm,mm,nnn,nn,1,nz,ierr)
  call fulfname(fname,path2data,'GLOBAL/uu.dat',ierr)
  call ctl_file_write(fname,    &     !file name
                      undef,    &     !value for undefined points
                     nx_glob,   &     !x-dimension
                     ny_glob,   &     !y-dimension
                          nz,   &     !z-dimension
                        nrec,   &     !t-dimension
                    xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                     xt(mmm),   &     !first x-value (if linear) or x-array (if levels)
                       dxst,    &     !x-step (if linear)
                   ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                     yt(nnn),   &     !first y-value (if linear) or x-array (if levels)
                       dyst,    &     !y-step (if linear)
                       1,       &     !z-grid type (0 - linear, 1 - levels)
                  z*1000.0d0,   &     !first z-value (if linear) or x-array (if levels)
                       1.0d0,   &     !z-step (if linear)
                    calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                        year,   &     !year   of the first field
                       month,   &     !month  of the first field
                         day,   &     !day    of the first field
                        hour,   &     !hour   of the first field
                      minute,   &     !minute of the first field
                       tstep,   &     !time step (in seconds
        'zonal velocity, m/s',  &     !title of dataset
                         'u'   )      !variable name
 endif
 uu_calc=0.0d0

!------------------------------------------------------------------------------------------------------
!writing meridional velocity
 ierr=0
 if(grid_shift==0) then !writing on the model grid
  array4_3d=sngl(vv_calc/dfloat(meancalc))
  call wdstd(path2data,'GLOBAL/vv.dat',nrec,array4_3d,llv,nx,ny,nz,mmm,mm,nnn-1,nn,1,nz,ierr)
  call fulfname(fname,path2data,'GLOBAL/vv.dat',ierr)
  call ctl_file_write(fname,    &     !file name
                      undef,    &     !value for undefined points
                     nx_glob,   &     !x-dimension
                   ny_glob+1,   &     !y-dimension
                          nz,   &     !z-dimension
                        nrec,   &     !t-dimension
                    xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                     xt(mmm),   &     !first x-value (if linear) or x-array (if levels)
                       dxst,    &     !x-step (if linear)
                   ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                   yv(nnn-1),   &     !first y-value (if linear) or x-array (if levels)
                       dyst,    &     !y-step (if linear)
                       1,       &     !z-grid type (0 - linear, 1 - levels)
                  z*1000.0d0,   &     !first z-value (if linear) or x-array (if levels)
                       1.0d0,   &     !z-step (if linear)
                    calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                        year,   &     !year   of the first field
                       month,   &     !month  of the first field
                         day,   &     !day    of the first field
                        hour,   &     !hour   of the first field
                      minute,   &     !minute of the first field
                       tstep,   &     !time step (in seconds
   'meridional velocity, m/s',  &     !title of dataset
                         'v'   )      !variable name

 else !writing on T-grid

 !$omp parallel do
  do n=ny_start, ny_end
   do m=nx_start, nx_end
    if(lu(m,n)>0.5) then
     do k=1,nz
     array4_3d(m,n,k)=sngl( (vv_calc(m,n  ,k)*dxh(m,n  )*dyt(m,n  )*hhv(m,n  )    &
                            +vv_calc(m,n-1,k)*dxh(m,n-1)*dyt(m,n-1)*hhv(m,n-1))/2.0d0/hhq(m,n)/dx(m,n)/dy(m,n)/dfloat(meancalc) )
     enddo
    endif
   enddo
  enddo
 !$omp end parallel do

  call wdstd(path2data,'GLOBAL/vv.dat',nrec,array4_3d,lu,nx,ny,nz,mmm,mm,nnn,nn,1,nz,ierr)
  call fulfname(fname,path2data,'GLOBAL/vv.dat',ierr)
  call ctl_file_write(fname,    &     !file name
                      undef,    &     !value for undefined points
                     nx_glob,   &     !x-dimension
                     ny_glob,   &     !y-dimension
                          nz,   &     !z-dimension
                        nrec,   &     !t-dimension
                    xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                     xt(mmm),   &     !first x-value (if linear) or x-array (if levels)
                       dxst,    &     !x-step (if linear)
                   ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                     yt(nnn),   &     !first y-value (if linear) or x-array (if levels)
                       dyst,    &     !y-step (if linear)
                       1,       &     !z-grid type (0 - linear, 1 - levels)
                  z*1000.0d0,   &     !first z-value (if linear) or x-array (if levels)
                       1.0d0,   &     !z-step (if linear)
                    calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                        year,   &     !year   of the first field
                       month,   &     !month  of the first field
                         day,   &     !day    of the first field
                        hour,   &     !hour   of the first field
                      minute,   &     !minute of the first field
                       tstep,   &     !time step (in seconds
   'meridional velocity, m/s',  &     !title of dataset
                         'v'   )      !variable name
 endif
endif
vv_calc=0.0d0

!-----------------------------------------------------------------------------------------------------
if(wstr_output>0) then
!writing zonal wind stress
 ierr=0
 if(grid_shift==0) then !writing on the model grid
  array4_2d=sngl(txo_calc/dfloat(meancalc))
  call wdstd(path2data,'GLOBAL/tx.dat',nrec,array4_2d,llu,nx,ny,1,mmm-1,mm,nnn,nn,1,1,ierr)
  call fulfname(fname,path2data,'GLOBAL/tx.dat',ierr)
  call ctl_file_write(fname,    &     !file name
                      undef,    &     !value for undefined points
                   nx_glob+1,   &     !x-dimension
                     ny_glob,   &     !y-dimension
                           1,   &     !z-dimension
                        nrec,   &     !t-dimension
                    xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                   xu(mmm-1),   &     !first x-value (if linear) or x-array (if levels)
                       dxst,    &     !x-step (if linear)
                   ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                     yt(nnn),   &     !first y-value (if linear) or x-array (if levels)
                       dyst,    &     !y-step (if linear)
                       1,       &     !z-grid type (0 - linear, 1 - levels)
                       0.0d0,   &     !first z-value (if linear) or x-array (if levels)
                       1.0d0,   &     !z-step (if linear)
                    calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                        year,   &     !year   of the first field
                       month,   &     !month  of the first field
                         day,   &     !day    of the first field
                        hour,   &     !hour   of the first field
                      minute,   &     !minute of the first field
                       tstep,   &     !time step (in seconds
 'zonal wind stress, (m/s)^2',  &     !title of dataset
                         'tx'   )      !variable name

 else !writing on T-grid

 !$omp parallel do
  do n=ny_start, ny_end
   do m=nx_start, nx_end
    if(lu(m,n)>0.5) then
     array4_2d(m,n)=sngl( (txo_calc(m  ,n)*dxt(m  ,n)*dyh(m  ,n)      &
                          +txo_calc(m-1,n)*dxt(m-1,n)*dyh(m-1,n))/2.0d0/dx(m,n)/dy(m,n)/dfloat(meancalc) )
    endif
   enddo
  enddo
 !$omp end parallel do

  call wdstd(path2data,'GLOBAL/tx.dat',nrec,array4_2d,lu,nx,ny,1,mmm,mm,nnn,nn,1,1,ierr)
  call fulfname(fname,path2data,'GLOBAL/tx.dat',ierr)
  call ctl_file_write(fname,    &     !file name
                      undef,    &     !value for undefined points
                     nx_glob,   &     !x-dimension
                     ny_glob,   &     !y-dimension
                           1,   &     !z-dimension
                        nrec,   &     !t-dimension
                    xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                     xt(mmm),   &     !first x-value (if linear) or x-array (if levels)
                       dxst,    &     !x-step (if linear)
                   ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                     yt(nnn),   &     !first y-value (if linear) or x-array (if levels)
                       dyst,    &     !y-step (if linear)
                       1,       &     !z-grid type (0 - linear, 1 - levels)
                       0.0d0,   &     !first z-value (if linear) or x-array (if levels)
                       1.0d0,   &     !z-step (if linear)
                    calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                        year,   &     !year   of the first field
                       month,   &     !month  of the first field
                         day,   &     !day    of the first field
                        hour,   &     !hour   of the first field
                      minute,   &     !minute of the first field
                       tstep,   &     !time step (in seconds
 'zonal wind stress, (m/s)^2',  &     !title of dataset
                        'tx'   )      !variable name
 endif
 txo_calc=0.0d0

!------------------------------------------------------------------------------------------------------
!writing meridional wind stress
 ierr=0
 if(grid_shift==0) then !writing on the model grid
  array4_2d=sngl(tyo_calc/dfloat(meancalc))
  call wdstd(path2data,'GLOBAL/ty.dat',nrec,array4_2d,llv,nx,ny,1,mmm,mm,nnn-1,nn,1,1,ierr)
  call fulfname(fname,path2data,'GLOBAL/ty.dat',ierr)
  call ctl_file_write(fname,    &     !file name
                      undef,    &     !value for undefined points
                     nx_glob,   &     !x-dimension
                   ny_glob+1,   &     !y-dimension
                           1,   &     !z-dimension
                        nrec,   &     !t-dimension
                    xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                     xt(mmm),   &     !first x-value (if linear) or x-array (if levels)
                       dxst,    &     !x-step (if linear)
                   ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                   yv(nnn-1),   &     !first y-value (if linear) or x-array (if levels)
                       dyst,    &     !y-step (if linear)
                       1,       &     !z-grid type (0 - linear, 1 - levels)
                       0.0d0,   &     !first z-value (if linear) or x-array (if levels)
                       1.0d0,   &     !z-step (if linear)
                    calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                        year,   &     !year   of the first field
                       month,   &     !month  of the first field
                         day,   &     !day    of the first field
                        hour,   &     !hour   of the first field
                      minute,   &     !minute of the first field
                       tstep,   &     !time step (in seconds
 'meridional wind stress, (m/s)^2',  &     !title of dataset
                         'ty'   )      !variable name

 else !writing on T-grid

 !$omp parallel do
  do n=ny_start, ny_end
   do m=nx_start, nx_end
    if(lu(m,n)>0.5) then
     array4_2d(m,n)=sngl( (tyo_calc(m,n  )*dxh(m,n  )*dyt(m,n  )  &
                          +tyo_calc(m,n-1)*dxh(m,n-1)*dyt(m,n-1))/2.0d0/dx(m,n)/dy(m,n)/dfloat(meancalc) )
    endif
   enddo
  enddo
 !$omp end parallel do

  call wdstd(path2data,'GLOBAL/ty.dat',nrec,array4_2d,lu,nx,ny,1,mmm,mm,nnn,nn,1,1,ierr)
  call fulfname(fname,path2data,'GLOBAL/ty.dat',ierr)
  call ctl_file_write(fname,    &     !file name
                      undef,    &     !value for undefined points
                     nx_glob,   &     !x-dimension
                     ny_glob,   &     !y-dimension
                           1,   &     !z-dimension
                        nrec,   &     !t-dimension
                    xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                     xt(mmm),   &     !first x-value (if linear) or x-array (if levels)
                       dxst,    &     !x-step (if linear)
                   ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                     yt(nnn),   &     !first y-value (if linear) or x-array (if levels)
                       dyst,    &     !y-step (if linear)
                       1,       &     !z-grid type (0 - linear, 1 - levels)
                       0.0d0,   &     !first z-value (if linear) or x-array (if levels)
                       1.0d0,   &     !z-step (if linear)
                    calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                        year,   &     !year   of the first field
                       month,   &     !month  of the first field
                         day,   &     !day    of the first field
                        hour,   &     !hour   of the first field
                      minute,   &     !minute of the first field
                       tstep,   &     !time step (in seconds
 'meridional wind stress, (m/s)^2',  &     !title of dataset
                         'ty'   )      !variable name
 endif
endif
tyo_calc=0.0d0

!--------------------------------------------------------------------------------
if(tt_output>0) then
!writing temperature
 ierr=0
 array4_3d=sngl(tt_calc/dfloat(meancalc))
 call wdstd(path2data,'GLOBAL/tt.dat',nrec,array4_3d,lu,nx,ny,nz,mmm,mm,nnn,nn,1,nz,ierr)
 call fulfname(fname,path2data,'GLOBAL/tt.dat',ierr)
 call ctl_file_write(fname,    &     !file name
                     undef,    &     !value for undefined points
                    nx_glob,   &     !x-dimension
                    ny_glob,   &     !y-dimension
                         nz,   &     !z-dimension
                       nrec,   &     !t-dimension
                   xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                    xt(mmm),   &     !first x-value (if linear) or x-array (if levels)
                      dxst,    &     !x-step (if linear)
                  ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                    yt(nnn),   &     !first y-value (if linear) or x-array (if levels)
                      dyst,    &     !y-step (if linear)
                      1,       &     !z-grid type (0 - linear, 1 - levels)
                 z*1000.0d0,   &     !first z-value (if linear) or x-array (if levels)
                      1.0d0,   &     !z-step (if linear)
                   calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                       year,   &     !year   of the first field
                      month,   &     !month  of the first field
                        day,   &     !day    of the first field
                       hour,   &     !hour   of the first field
                     minute,   &     !minute of the first field
                      tstep,   &     !time step (in seconds)
          'Temperature, �C',   &     !title of dataset
                       'tt'   )      !variable name
endif
tt_calc=0.0d0

!--------------------------------------------------------------------------------
if(ss_output>0) then
!writing temperature
 ierr=0
 array4_3d=sngl(ss_calc/dfloat(meancalc))
 call wdstd(path2data,'GLOBAL/ss.dat',nrec,array4_3d,lu,nx,ny,nz,mmm,mm,nnn,nn,1,nz,ierr)
 call fulfname(fname,path2data,'GLOBAL/ss.dat',ierr)
 call ctl_file_write(fname,    &     !file name
                     undef,    &     !value for undefined points
                    nx_glob,   &     !x-dimension
                    ny_glob,   &     !y-dimension
                         nz,   &     !z-dimension
                       nrec,   &     !t-dimension
                   xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                    xt(mmm),   &     !first x-value (if linear) or x-array (if levels)
                      dxst,    &     !x-step (if linear)
                  ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                    yt(nnn),   &     !first y-value (if linear) or x-array (if levels)
                      dyst,    &     !y-step (if linear)
                      1,       &     !z-grid type (0 - linear, 1 - levels)
                 z*1000.0d0,   &     !first z-value (if linear) or x-array (if levels)
                      1.0d0,   &     !z-step (if linear)
                   calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                       year,   &     !year   of the first field
                      month,   &     !month  of the first field
                        day,   &     !day    of the first field
                       hour,   &     !hour   of the first field
                     minute,   &     !minute of the first field
                      tstep,   &     !time step (in seconds)
            'Salinity, PSU',   &     !title of dataset
                       'ss'   )      !variable name
endif
ss_calc=0.0d0

!--------------------------------------------------------------------------------------------------
!writing surface fluxes
if(sfl_output>0) then
 !$omp parallel do
  do n=ny_start, ny_end
   do m=nx_start, nx_end
    if(lu(m,n)>0.5) then
     do k=1,16
      array4_3d(m,n,k)=sngl(sfl_calc(m,n,k)/dfloat(meancalc))
     enddo
    endif
   enddo
  enddo
 !$omp end parallel do
ierr=0
 call wdstd(path2data,'GLOBAL/sfl.dat',nrec,array4_3d,lu,nx,ny,16,mmm,mm,nnn,nn,1,16,ierr)
 call fulfname(fname,path2data,'GLOBAL/sfl.dat',ierr)
 call ctl_file_write(fname,    &     !file name
                     undef,    &     !value for undefined points
                    nx_glob,   &     !x-dimension
                    ny_glob,   &     !y-dimension
                         16,   &     !z-dimension
                       nrec,   &     !t-dimension
                   xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                    xt(mmm),   &     !first x-value (if linear) or x-array (if levels)
                      dxst,    &     !x-step (if linear)
                  ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                    yt(nnn),   &     !first y-value (if linear) or x-array (if levels)
                      dyst,    &     !y-step (if linear)
                      0,       &     !z-grid type (0 - linear, 1 - levels)
                      1.0d0,   &     !first z-value (if linear) or x-array (if levels)
                      1.0d0,   &     !z-step (if linear)
                   calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                       year,   &     !year   of the first field
                      month,   &     !month  of the first field
                        day,   &     !day    of the first field
                       hour,   &     !hour   of the first field
                     minute,   &     !minute of the first field
                      tstep,   &     !time step (in seconds)
'1-TA, 2-QA, 3-LWdw, 4-SWdw, 5-SLP, 6-SensBal, 7-LatBal, 8-LWbal, 9-SWbal, 10-HFtot, 11-rain, 12-snow, 13-roff, 14-WFtot, 15-Tflux, 16-Sflux',   &     !title of dataset
                      'sfl'   )      !variable name

endif
sfl_calc=0.0d0

!-----------------------------------------------------------------------------------------------
!write wind speed
if(wind_output>0) then

  array4_2d=sngl(uwnd_calc/dfloat(meancalc))

  call wdstd(path2data,'GLOBAL/uwnd.dat',nrec,array4_2d,lu,nx,ny,1,mmm,mm,nnn,nn,1,1,ierr)
  call fulfname(fname,path2data,'GLOBAL/uwnd.dat',ierr)
  call ctl_file_write(fname,    &     !file name
                      undef,    &     !value for undefined points
                     nx_glob,   &     !x-dimension
                     ny_glob,   &     !y-dimension
                           1,   &     !z-dimension
                        nrec,   &     !t-dimension
                    xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                     xt(mmm),   &     !first x-value (if linear) or x-array (if levels)
                       dxst,    &     !x-step (if linear)
                   ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                     yt(nnn),   &     !first y-value (if linear) or x-array (if levels)
                       dyst,    &     !y-step (if linear)
                       1,       &     !z-grid type (0 - linear, 1 - levels)
                       0.0d0,   &     !first z-value (if linear) or x-array (if levels)
                       1.0d0,   &     !z-step (if linear)
                    calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                        year,   &     !year   of the first field
                       month,   &     !month  of the first field
                         day,   &     !day    of the first field
                        hour,   &     !hour   of the first field
                      minute,   &     !minute of the first field
                       tstep,   &     !time step (in seconds
      'zonal wind speed, m/s',  &     !title of dataset
                         'u'   )      !variable name

  array4_2d=sngl(vwnd_calc/dfloat(meancalc))

  call wdstd(path2data,'GLOBAL/vwnd.dat',nrec,array4_2d,lu,nx,ny,1,mmm,mm,nnn,nn,1,1,ierr)
  call fulfname(fname,path2data,'GLOBAL/vwnd.dat',ierr)
  call ctl_file_write(fname,    &     !file name
                      undef,    &     !value for undefined points
                     nx_glob,   &     !x-dimension
                     ny_glob,   &     !y-dimension
                           1,   &     !z-dimension
                        nrec,   &     !t-dimension
                    xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                     xt(mmm),   &     !first x-value (if linear) or x-array (if levels)
                       dxst,    &     !x-step (if linear)
                   ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                     yt(nnn),   &     !first y-value (if linear) or x-array (if levels)
                       dyst,    &     !y-step (if linear)
                       1,       &     !z-grid type (0 - linear, 1 - levels)
                       0.0d0,   &     !first z-value (if linear) or x-array (if levels)
                       1.0d0,   &     !z-step (if linear)
                    calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                        year,   &     !year   of the first field
                       month,   &     !month  of the first field
                         day,   &     !day    of the first field
                        hour,   &     !hour   of the first field
                      minute,   &     !minute of the first field
                       tstep,   &     !time step (in seconds
 'meridional wind speed, m/s',  &     !title of dataset
                         'v'   )      !variable name
endif
uwnd_calc=0.0d0
vwnd_calc=0.0d0

meancalc=0


endsubroutine global_output



subroutine parallel_local_output(path2data,  &
                                 nrec,       &
                                 year,       &
                                month,       &
                                  day,       &
                                 hour,       &
                               minute,       &
                               tstep,        &
                               calendar  )
use main_basin_pars
use mpi_parallel_tools
use basin_grid
use ocean_variables
use ocean_bc
use rec_length

implicit none
include 'locout.fi'

integer nrec, year, month, day, hour, minute, calendar, ierr
character fname*256
character*(*) path2data
real(4) array4_2d(bnd_x1:bnd_x2,bnd_y1:bnd_y2),       &
        array4_3d(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz) !,    &
!        array4_2dn(nx_start:nx_end, ny_start:ny_end)
real(4) tstep
integer m,n,k

if (rank .eq. 0) write(*,*) 'Writing local output, record number ', nrec

if(nrec==1) then
!writing HHQ
 ierr=0
 array4_2d=sngl(hhq_rest)
! array4_2dn(nx_start:nx_end, ny_start:ny_end) = sngl(hhq_rest(nx_start:nx_end, ny_start:ny_end))
! call wdstd(path2data,'LOCAL/hhq.dat',nrec,array4_2d,lu,nx,ny,1,m1loc,m2loc,n1loc,n2loc,1,1,ierr)
 call pwdstd(path2data, 'LOCAL/hhq.dat', nrec,array4_2d, lu,         &
            nx, nx_start, nx_end,                                    &
            ny, ny_start, ny_end, 1,                                 &
            m1loc, m2loc, n1loc, n2loc, 1, 1, ierr)

 call fulfname(fname,path2data,'LOCAL/hhq.dat',ierr)
 if (rank .eq. 0) then
     call ctl_file_write(fname,    &     !file name
                     undef,    &     !value for undefined points
                     nx_loc,   &     !x-dimension
                     ny_loc,   &     !y-dimension
                          1,   &     !z-dimension
                       nrec,   &     !t-dimension
                   xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                  xt(m1loc),   &     !first x-value (if linear) or x-array (if levels)
                      dxst,    &     !x-step (if linear)
                  ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                  yt(n1loc),   &     !first y-value (if linear) or x-array (if levels)
                      dyst,    &     !y-step (if linear)
                      0,       &     !z-grid type (0 - linear, 1 - levels)
                      0.0d0,   &     !first z-value (if linear) or x-array (if levels)
                      1.0d0,   &     !z-step (if linear)
                   calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                       year,   &     !year   of the first field
                      month,   &     !month  of the first field
                        day,   &     !day    of the first field
                       hour,   &     !hour   of the first field
                     minute,   &     !minute of the first field
                      tstep,   &     !time step (in seconds)
                   'HHQ, m',   &     !title of dataset
                      'hhq'   )      !variable name
 endif
endif

if(ssh_output>0) then
!if (0 .eq. 1) then
! writing SSH
 ierr=0
 array4_2d=sngl(ssh_i)
! call wdstd(path2data,'LOCAL/ssh.dat',nrec,array4_2d,lu,nx,ny,1,m1loc,m2loc,n1loc,n2loc,1,1,ierr)
 call pwdstd(path2data, 'LOCAL/ssh.dat', nrec,array4_2d, lu,         &
            nx, nx_start, nx_end,                                    &
            ny, ny_start, ny_end, 1,                                 &
            m1loc, m2loc, n1loc, n2loc, 1, 1, ierr)

 call fulfname(fname,path2data,'LOCAL/ssh.dat',ierr)
 if (rank .eq. 0) then
     call ctl_file_write(fname,    &     !file name
                     undef,    &     !value for undefined points
                     nx_loc,   &     !x-dimension
                     ny_loc,   &     !y-dimension
                          1,   &     !z-dimension
                       nrec,   &     !t-dimension
                   xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                  xt(m1loc),   &     !first x-value (if linear) or x-array (if levels)
                      dxst,    &     !x-step (if linear)
                  ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                  yt(n1loc),   &     !first y-value (if linear) or x-array (if levels)
                      dyst,    &     !y-step (if linear)
                      0,       &     !z-grid type (0 - linear, 1 - levels)
                      0.0d0,   &     !first z-value (if linear) or x-array (if levels)
                      1.0d0,   &     !z-step (if linear)
                   calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                       year,   &     !year   of the first field
                      month,   &     !month  of the first field
                        day,   &     !day    of the first field
                       hour,   &     !hour   of the first field
                     minute,   &     !minute of the first field
                      tstep,   &     !time step (in seconds)
                   'SSH, m',   &     !title of dataset
                      'ssh'   )      !variable name
 endif

! writing SSH error
 ierr=0
 array4_2d=sngl(ssh_i - ssh_err)
! array4_2d=sngl(ssh_err)
! call wdstd(path2data,'LOCAL/ssh_err.dat',nrec,array4_2d,lu,nx,ny,1,m1loc,m2loc,n1loc,n2loc,1,1,ierr)
 call pwdstd(path2data, 'LOCAL/ssh_err.dat', nrec,array4_2d, lu,         &
            nx, nx_start, nx_end,                                        &
            ny, ny_start, ny_end, 1,                                     &
            m1loc, m2loc, n1loc, n2loc, 1, 1, ierr)

 call fulfname(fname,path2data,'LOCAL/ssh_err.dat',ierr)
 if (rank .eq. 0) then
     call ctl_file_write(fname,    &     !file name
                   undef,    &     !value for undefined points
                   nx_loc,   &     !x-dimension
                   ny_loc,   &     !y-dimension
                        1,   &     !z-dimension
                     nrec,   &     !t-dimension
                 xgr_type,   &     !x-grid type (0 - linear, 1 - levels)
                xt(m1loc),   &     !first x-value (if linear) or x-array (if levels)
                    dxst,    &     !x-step (if linear)
                ygr_type,    &     !y-grid type (0 - linear, 1 - levels)
                yt(n1loc),   &     !first y-value (if linear) or x-array (if levels)
                    dyst,    &     !y-step (if linear)
                    0,       &     !z-grid type (0 - linear, 1 - levels)
                    0.0d0,   &     !first z-value (if linear) or x-array (if levels)
                    1.0d0,   &     !z-step (if linear)
                 calendar,   &     !type   of calendar (0 - without leap-year, 1 - with leap-year)
                     year,   &     !year   of the first field
                    month,   &     !month  of the first field
                      day,   &     !day    of the first field
                     hour,   &     !hour   of the first field
                   minute,   &     !minute of the first field
                    tstep,   &     !time step (in seconds)
                 'SSH_error, m',   &     !title of dataset
                    'ssh_err'   )      !variable name
 endif
endif

endsubroutine parallel_local_output
