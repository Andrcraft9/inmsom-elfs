subroutine print_basin_grid
    use main_basin_pars
    use mpi_parallel_tools
    use basin_grid
    use ocean_variables

    implicit none

    integer :: m, n, k, procn, ierr

    call mpi_comm_size(cart_comm, procn, ierr)
    do k = 0, procn-1
        if (rank .eq. k) then
            print *, "RANK IS ", rank
            ! LU mask
            print *, "|------------------------- LU MASK: --------------------------|"
            print *, "m ", "n "
            do m=bnd_x1, bnd_x2
              do n=bnd_y1, bnd_y2
                  print *, m, n, lu(m, n)
              enddo
            enddo

            ! LCU mask
            print *, "|------------------------- LCU MASK: --------------------------|"
            print *, "m ", "n "
            do m=bnd_x1, bnd_x2
              do n=bnd_y1, bnd_y2
                  print *, m, n, lcu(m, n)
              enddo
            enddo

            ! LCV mask
            print *, "|------------------------- LCV MASK: --------------------------|"
            print *, "m ", "n "
            do m=bnd_x1, bnd_x2
              do n=bnd_y1, bnd_y2
                  print *, m, n, lcv(m, n)
              enddo
            enddo
        endif
        call mpi_barrier(cart_comm, ierr)
    enddo

end subroutine print_basin_grid

subroutine parallel_point_output(path2data, nstep, lon, lat, name)
    use main_basin_pars
    use mpi_parallel_tools
    use basin_grid
    use ocean_variables

    implicit none

    character fname*256
    character(20) name
    character*(*) path2data
    integer :: nstep, m, n, r, ierr
    real*8 :: lon, lat

!------------------------------- First point -----------------------------------!
    m = floor((lon - rlon) / dxst) + mmm
    n = floor((lat - rlat) / dyst) + nnn

    r = get_rank_by_point(m, n)

    if (rank .eq. r) then
!        print *, rank, lon, lat, geo_lon_t(m, n), geo_lat_t(m, n)
        call fulfname(fname, path2data, name, ierr)
        open(40, file=fname, status='unknown', position='append')
        write(40, *) nstep, ssh_i(m, n)
        close(40)
    endif

    return
end subroutine parallel_point_output

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
 call pwdstd(path2data,'LOCAL/hhq.dat',nrec,array4_2d,lu,nx,ny,1,m1loc,m2loc,n1loc,n2loc,1,1,ierr)

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
! writing SSH
 ierr=0
 array4_2d=sngl(ssh_i)
! call wdstd(path2data,'LOCAL/ssh.dat',nrec,array4_2d,lu,nx,ny,1,m1loc,m2loc,n1loc,n2loc,1,1,ierr)
 call pwdstd(path2data,'LOCAL/ssh.dat',nrec,array4_2d,lu,nx,ny,1,m1loc,m2loc,n1loc,n2loc,1,1,ierr)
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

 ierr = 0
 array4_2d=sngl(ubrtr_i)
 call pwdstd(path2data,'LOCAL/ubrtr.dat',nrec,array4_2d,llu,nx,ny,1,m1loc-1,m2loc,n1loc,n2loc,1,1,ierr)
 call fulfname(fname,path2data,'LOCAL/ubrtr.dat',ierr)
 if (rank .eq. 0) then
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
                      0,       &     !z-grid type (0 - linear, 1 - levels)
                      0.0d0,   &     !first z-value (if linear) or x-array (if levels)
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

 ierr=0
 array4_2d=sngl(vbrtr_i)
 call pwdstd(path2data,'LOCAL/vbrtr.dat',nrec,array4_2d,llv,nx,ny,1,m1loc,m2loc,n1loc-1,n2loc,1,1,ierr)
 call fulfname(fname,path2data,'LOCAL/vbrtr.dat',ierr)
 if (rank .eq. 0) then
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
                       0,       &     !z-grid type (0 - linear, 1 - levels)
                       0.0d0,   &     !first z-value (if linear) or x-array (if levels)
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
!if (rank .eq. 0) print *, "---------------------SSH OUT-------------------------"

endsubroutine parallel_local_output
