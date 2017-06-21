!----------------------array boundary definition for non-mpi arrays
subroutine non_mpi_array_boundary_definition
 use main_basin_pars
 use mpi_parallel_tools
 implicit none

       nx_start=mmm
       nx_end  =mm
       ny_start =nnn
       ny_end  =nn

       bnd_x1=nx_start-2
       bnd_x2=nx_end  +2
       bnd_y1=ny_start-2
       bnd_y2=ny_end  +2

endsubroutine non_mpi_array_boundary_definition
!end of array boundary definition for non-mpi arrays

!----------------------array boundary definition for mpi arrays
subroutine mpi_array_boundary_definition
    use main_basin_pars
    use mpi_parallel_tools
    implicit none

    include "omp_lib.h"

    integer :: ierr, procn, i
    integer :: locn
    integer :: count_threads, num_thread

    call mpi_init(ierr)

    period = (/1,1/)
    p_size = (/0,0/)
    ierr = 0

    call mpi_comm_rank(mpi_comm_world, rank, ierr)
    call mpi_comm_size(mpi_comm_world, procs, ierr)
    call mpi_dims_create(procs, 2, p_size, ierr)
    call mpi_cart_create(mpi_comm_world, 2, p_size, period, 0, cart_comm, ierr)
    call mpi_cart_coords(cart_comm, rank, 2, p_coord, ierr)

!-----------------------------------NX------------------------------------------!
    locn = floor(real(nx - 4)/real(p_size(1)))
    nx_start = locn*p_coord(1) + 1 + 2
    if ( p_coord(1) .EQ. p_size(1) - 1 ) then
        locn = (nx - 2) - nx_start + 1
    endif
    nx_end = nx_start + locn - 1
    nx_start = nx_start
!   border area
    bnd_x1 = nx_start - 2
!    if (bnd_x1 < 1) bnd_x1 = 1
    bnd_x2 = nx_end + 2
!    if (bnd_x2 > nx) bnd_x2 = nx

!-----------------------------------NY------------------------------------------!
    locn = floor(real(ny - 4)/real(p_size(2)))
    ny_start = locn*p_coord(2) + 1 + 2
    if ( p_coord(2) .EQ. p_size(2) - 1 ) then
        locn = (ny - 2) - ny_start + 1
    endif
    ny_end = ny_start + locn - 1
    ny_start = ny_start
!   border area
    bnd_y1 = ny_start - 2
!    if (bnd_y1 < 1) bnd_y1 = 1
    bnd_y2 = ny_end + 2
!    if (bnd_y2 > ny) bnd_y2 = ny

    call mpi_comm_size(cart_comm, procn, ierr)
    if (rank .eq. 0) print *, "MPI pocs: ", procn, " Domain decomposition:"
    do i = 0, procn-1
        if (rank .eq. i) then
            print *, "nx ", rank, p_coord, nx_start, nx_end, ny_start, ny_end
            print *, "bnd", rank, p_coord, bnd_x1, bnd_x2, bnd_y1, bnd_y2
        endif
        call mpi_barrier(cart_comm, ierr)
    enddo

!$omp parallel
    count_threads = omp_get_num_threads()
    num_thread = omp_get_thread_num()
    if (num_thread .eq. 0) print *, "OMP Threads: ", count_threads
!$omp end parallel

endsubroutine mpi_array_boundary_definition


!-----------------------------------------allocation of arrays
   subroutine model_grid_allocate
   use main_basin_pars
   use mpi_parallel_tools
   use basin_grid
   implicit none

      allocate  ( lu(bnd_x1:bnd_x2,bnd_y1:bnd_y2),             &  !mask of t-grid
                 lu1(bnd_x1:bnd_x2,bnd_y1:bnd_y2),             &  !mask of t-grid (1 everywhere)
                 luu(bnd_x1:bnd_x2,bnd_y1:bnd_y2),             &  !mask of h-grid (0 on boundary)
                 luh(bnd_x1:bnd_x2,bnd_y1:bnd_y2),             &  !mask of h-grid (1 on boundary)
                 lcu(bnd_x1:bnd_x2,bnd_y1:bnd_y2),             &  !mask of u-grid (0 on boundary)
                 lcv(bnd_x1:bnd_x2,bnd_y1:bnd_y2),             &  !mask of v-grid (0 on boundary)
                 llu(bnd_x1:bnd_x2,bnd_y1:bnd_y2),             &  !mask of u-grid (0 on boundary)
                 llv(bnd_x1:bnd_x2,bnd_y1:bnd_y2) )                !mask of v-grid (0 on boundary)
      allocate   (lbasins(nx,ny))       !integer masks of regional basins
      allocate   (      hhh(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &  !ocean depth on luh (h-points)
                       hhhp(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &  !ocean depth on luh (h-points) at previous step
                       hhhn(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &  !ocean depth on luh (h-points) at pre-previous step
                   hhq_rest(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &  !ocean depth on lu  (t-points) at rest state
                        hhq(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &  !ocean depth on lu  (t-points)
                       hhqp(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &  !ocean depth on lu  (t-points) at previous step
                       hhqn(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &  !ocean depth on lu  (t-points) at pre-previous step
                        hhu(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &  !ocean depth on lcu (u-points)
                       hhup(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &  !ocean depth on lcu (u-points) at previous step
                       hhun(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &  !ocean depth on lcu (u-points) at pre-previous step
                        hhv(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &  !ocean depth on lcv (v-points)
                       hhvp(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &  !ocean depth on lcv (v-points) at previous step
                       hhvn(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &  !ocean depth on lcv (v-points) at pre-previous step
                        rlh_s(bnd_x1:bnd_x2,bnd_y1:bnd_y2),    &  !coriolis-1 parameter on edge (t-centers) points
                        rlh_c(bnd_x1:bnd_x2,bnd_y1:bnd_y2),    &  !coriolis-2 parameter on edge (t-centers) points
                     z(nz), zw(nz+1),    &  !vertical sigma-levels (t-points and w-points)
                 hzt(nz+1), dz(nz),      &  !steps between t-levels and w-levels
             dxt(bnd_x1:bnd_x2,bnd_y1:bnd_y2),dyt(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &  !horizontal grid steps between   t-points (in radians or meters)
             dx (bnd_x1:bnd_x2,bnd_y1:bnd_y2),dy (bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &  !horizontal grid steps between u,v-points (in radians or meters)
             dxh(bnd_x1:bnd_x2,bnd_y1:bnd_y2),dyh(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &  !horizontal grid steps between   h-points (in radians or meters)
             dxb(bnd_x1:bnd_x2,bnd_y1:bnd_y2),dyb(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &  !horizontal grid steps between v,u-points (in radians or meters)
                   xt(bnd_x1:bnd_x2),yt(bnd_y1:bnd_y2),        &  !horizontal t-grid            x- and y-coordinates (in degrees)
                   xu(bnd_x1:bnd_x2),yv(bnd_y1:bnd_y2)    )       !horizontal u-grid and v-grid x- and y-coordinates (in degrees)

      allocate  (  geo_lon_t(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &    !geographical longitudes of T-points
                   geo_lat_t(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &    !geographical latitudes  of T-points
                   geo_lon_u(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &    !geographical longitudes of U-points
                   geo_lat_u(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &    !geographical latitudes  of U-points
                   geo_lon_v(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &    !geographical longitudes of V-points
                   geo_lat_v(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &    !geographical latitudes  of V-points
                   geo_lon_h(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &    !geographical longitudes of H-points
                   geo_lat_h(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &    !geographical latitudes  of H-points
                rotvec_coeff(bnd_x1:bnd_x2,bnd_y1:bnd_y2,4)    )

      allocate   (      hhh_e(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &  !ocean depth on luh (h-points)
                       hhhp_e(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &  !ocean depth on luh (h-points) at previous step
                       hhhn_e(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &  !ocean depth on luh (h-points) at pre-previous step
                        hhq_e(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &  !ocean depth on lu  (t-points)
                       hhqp_e(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &  !ocean depth on lu  (t-points) at previous step
                       hhqn_e(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &  !ocean depth on lu  (t-points) at pre-previous step
                        hhu_e(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &  !ocean depth on lcu (u-points)
                       hhup_e(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &  !ocean depth on lcu (u-points) at previous step
                       hhun_e(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &  !ocean depth on lcu (u-points) at pre-previous step
                        hhv_e(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &  !ocean depth on lcv (v-points)
                       hhvp_e(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &  !ocean depth on lcv (v-points) at previous step
                       hhvn_e(bnd_x1:bnd_x2,bnd_y1:bnd_y2)    )     !ocean depth on lcv (v-points) at pre-previous step

         lu=0.0; lu1=0.0; luu=0.0; luh=0.0; lcu=0.0; lcv=0.0; llu=0.0; llv=0.0; lbasins=0
         hhh=0.0d0; hhhp=0.0d0; hhhn=0.0d0; hhq_rest=0.0d0
         hhq=0.0d0; hhqp=0.0d0; hhqn=0.0d0
         hhu=0.0d0; hhup=0.0d0; hhun=0.0d0
         hhv=0.0d0; hhvp=0.0d0; hhvn=0.0d0
         rlh_s=0.0d0; rlh_c=0.0d0
         z=0.0d0; zw=0.0d0; hzt=0.0d0; dz=0.0d0
         dxt=0.0d0; dyt=0.0d0; dx=0.0d0; dy=0.0d0; dxh=0.0d0; dyh=0.0d0; dxb=0.0d0; dyb=0.0d0
         xt=0.0d0; yt=0.0d0; xu=0.0d0; yv=0.0d0

         geo_lon_t=0.0d0; geo_lat_t=0.0d0; geo_lon_u=0.0d0; geo_lat_u=0.0d0
         geo_lon_v=0.0d0; geo_lat_v=0.0d0; geo_lon_h=0.0d0; geo_lat_h=0.0d0
         rotvec_coeff=0.0d0

         hhh_e=0.0d0; hhhp_e=0.0d0; hhhn_e=0.0d0
         hhq_e=0.0d0; hhqp_e=0.0d0; hhqn_e=0.0d0
         hhu_e=0.0d0; hhup_e=0.0d0; hhun_e=0.0d0
         hhv_e=0.0d0; hhvp_e=0.0d0; hhvn_e=0.0d0

  endsubroutine model_grid_allocate

!deallocation of arrays
  subroutine model_grid_deallocate
  use basin_grid
  implicit none
           deallocate(hhvn_e,hhvp_e,hhv_e,hhun_e,hhup_e,hhu_e,hhqn_e,hhqp_e,hhq_e,hhhn_e,hhhp_e,hhh_e)
           deallocate(rotvec_coeff)
           deallocate(geo_lat_h,geo_lon_h,geo_lat_v,geo_lon_v,geo_lat_u,geo_lon_u,geo_lat_t,geo_lon_t)
           deallocate(yv,xu,yt,xt,dyb,dxb,dyh,dxh,dy,dx,dyt,dxt,dz,hzt,zw,z,rlh_c,rlh_s)
           deallocate(hhvn,hhvp,hhv,hhun,hhup,hhu,hhqn,hhqp,hhq,hhq_rest,hhhn,hhhp,hhh)
           deallocate(lbasins)
           deallocate(llv,llu,lcv,lcu,luh,luu,lu1,lu)

  endsubroutine model_grid_deallocate
!-------------------------------------------------------------------------------------------

!allocation of arrays
!--------------------------------------------------------------------------------------------
   subroutine ocean_variables_allocate
   use main_basin_pars
   use mpi_parallel_tools
   use ocean_variables
   implicit none

      allocate (ssh_i (bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &  !sea surface height (SSH) at current  time step [m] (internal mode)
                sshp_i(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &  !sea surface height (SSH) at previous time step [m] (internal mode)
                  pgrx(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &  !pressure gradient x-component for RHS
                  pgry(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &  !pressure gradient y-component for RHS
               ubrtr_i(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &  !barotropic velocity      zonal[m/s] at current  time step [m] (internal mode)
               vbrtr_i(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &  !barotropic velocity meridional[m/s] at current  time step [m] (internal mode)
                RHSx2d(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &  !x-component of external force(barotropic)
                RHSy2d(bnd_x1:bnd_x2,bnd_y1:bnd_y2))        !y-component of external force(barotropic)

      allocate (ssh_e(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &  !sea surface height (SSH) at current  time step [m] (external mode)
               sshp_e(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &  !sea surface height (SSH) at previous time step [m] (external mode)
              ubrtr_e(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &  !barotropic velocity      zonal[m/s]                           (external mode)
             ubrtrp_e(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &  !barotropic velocity      zonal[m/s] at previous time step [m] (external mode)
              vbrtr_e(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &  !barotropic velocity meridional[m/s]                           (external mode)
             vbrtrp_e(bnd_x1:bnd_x2,bnd_y1:bnd_y2))         !barotropic velocity meridional[m/s] at previous time step [m] (external mode)

      allocate (xxt(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),    &  !auxiliary array 1
                yyt(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz))       !auxiliary array 2

      allocate (wf_tot(bnd_x1:bnd_x2,bnd_y1:bnd_y2))       !total water flux

      allocate(BottomFriction(bnd_x1:bnd_x2,bnd_y1:bnd_y2),        &
                       r_diss(bnd_x1:bnd_x2,bnd_y1:bnd_y2))

      allocate(amuv2d(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &    !depth mean lateral viscosity
               amuv42d(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &    !depth mean lateral viscosity
              r_vort2d(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &    !relative vorticity of depth mean velocity
            stress_t2d(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &    !Horizontal tension tensor component (barotropic)
            stress_s2d(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &    !Horizontal shearing tensor component(barotropic)
      RHSx2d_tran_disp(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &    !dispersion x-component of external force(barotropic)
      RHSy2d_tran_disp(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &    !dispersion x-component of external force(barotropic)
      RHSx2d_diff_disp(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &    !dispersion x-component of external force(barotropic)
      RHSy2d_diff_disp(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &    !dispersion x-component of external force(barotropic)
      RHSx2d_bfc(bnd_x1:bnd_x2,bnd_y1:bnd_y2),           &
      RHSy2d_bfc(bnd_x1:bnd_x2,bnd_y1:bnd_y2))


     ssh_i =0.0d0; sshp_i=0.0d0
     pgrx=0.0d0; pgry=0.0d0
     ubrtr_i=0.0d0;  vbrtr_i=0.0d0
     RHSx2d=0.0d0; RHSy2d=0.0d0

     ssh_e=0.0d0; sshp_e=0.0d0
     ubrtr_e=0.0d0; ubrtrp_e=0.0d0
     vbrtr_e=0.0d0; vbrtrp_e=0.0d0

     xxt=0.0d0; yyt=0.0d0

     wf_tot=0.0d0

     BottomFriction=0.0d0; r_diss=0.0d0

     amuv2d=0.0d0; amuv42d=0.0d0; r_vort2d=0.0d0
     stress_t2d=0.0d0; stress_s2d=0.0d0
     RHSx2d_tran_disp=0.0d0; RHSy2d_tran_disp=0.0d0
     RHSx2d_diff_disp=0.0d0; RHSy2d_diff_disp=0.0d0

     RHSx2d_bfc = 0.0d0; RHSy2d_bfc = 0.0d0

   endsubroutine ocean_variables_allocate

!deallocation of arrays
   subroutine ocean_variables_deallocate
   use ocean_variables
   implicit none

   deallocate(RHSx2d_bfc, RHSy2d_bfc)
   deallocate(RHSy2d_diff_disp,RHSx2d_diff_disp,RHSy2d_tran_disp,RHSx2d_tran_disp,     &
              stress_s2d,stress_t2d,r_vort2d,amuv42d,amuv2d)
   deallocate(r_diss, BottomFriction)
   deallocate(wf_tot)
   deallocate(xxt, yyt)
   deallocate(vbrtrp_e,vbrtr_e,ubrtrp_e,ubrtr_e,sshp_e,ssh_e)
   deallocate(RHSy2d,RHSx2d,vbrtr_i,ubrtr_i,pgry,pgrx,sshp_i,ssh_i)

   endsubroutine ocean_variables_deallocate
