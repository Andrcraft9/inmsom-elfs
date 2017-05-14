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

    bnd_length = 10 ! Need set bnd_length=2k, k >= 2

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
    bnd_x1 = nx_start - bnd_length - 1
    if (bnd_x1 < 1) bnd_x1 = 1
    bnd_x2 = nx_end + bnd_length + 1
    if (bnd_x2 > nx) bnd_x2 = nx

!-----------------------------------NY------------------------------------------!
    locn = floor(real(ny - 4)/real(p_size(2)))
    ny_start = locn*p_coord(2) + 1 + 2
    if ( p_coord(2) .EQ. p_size(2) - 1 ) then
        locn = (ny - 2) - ny_start + 1
    endif
    ny_end = ny_start + locn - 1
    ny_start = ny_start
!   border area
    bnd_y1 = ny_start - bnd_length - 1
    if (bnd_y1 < 1) bnd_y1 = 1
    bnd_y2 = ny_end + bnd_length + 1
    if (bnd_y2 > ny) bnd_y2 = ny

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
                RHSy2d(bnd_x1:bnd_x2,bnd_y1:bnd_y2)      )  !y-component of external force(barotropic)

      allocate (ssh_e(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &  !sea surface height (SSH) at current  time step [m] (external mode)
               sshp_e(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &  !sea surface height (SSH) at previous time step [m] (external mode)
              ubrtr_e(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &  !barotropic velocity      zonal[m/s]                           (external mode)
             ubrtrp_e(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &  !barotropic velocity      zonal[m/s] at previous time step [m] (external mode)
              vbrtr_e(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &  !barotropic velocity meridional[m/s]                           (external mode)
             vbrtrp_e(bnd_x1:bnd_x2,bnd_y1:bnd_y2)  )       !barotropic velocity meridional[m/s] at previous time step [m] (external mode)

      allocate (uu(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),    &  !     zonal velocity [m/s]
               uup(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),    &  !     zonal velocity [m/s]
                vv(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),    &  !meridional velocity [m/s]
               vvp(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),    &  !meridional velocity [m/s]
                ww(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1),  &  !  vertical velocity in sigma-coord [m/s]
                tt(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),    &  !potential temperature[�C]
               ttp(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),    &  !potential temperature[�C]
                ss(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),    &  !salinity [psu-35ppt]
               ssp(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),    &  !salinity [psu-35ppt]
               den(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),    &  !in-situ density  [kg/m^3]
           den_pot(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),    &  !potential density  [kg/m^3]
            RHSx3d(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),    &  !x-component of external force(baroclinic)
            RHSy3d(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),    &  !y-component of external force(baroclinic)
       RHSx3d_tran(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),    &  !x-component of external force(baroclinic)
       RHSy3d_tran(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),    &  !y-component of external force(baroclinic)
               age(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),    &  !water ideal age [days]
              agep(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),    &  !water ideal age [days]
               xxt(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),    &  !auxiliary array 1
               yyt(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz)     )  !auxiliary array 2

     allocate ( q2(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1),   &  !turbulent kinetic energy (m/s)^2
               q2p(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1),   &  !turbulent kinetic energy at the previous time step (m/s)^2
               q2l(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1),   &  !turbulent kinetic energy by turbulent length scale (m/s)^2*m
              q2lp(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1)    )  !turbulent kinetic energy by turbulent length scale at the previous time step (m/s)^2*m

     allocate ( RHS_q2 (bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1),    &  !RHS for turbulent kinetic energy
                RHS_q2l(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1) )      !RHS for turbulent kinetic energy by turbulent length scale

      allocate (uu2d(bnd_x1:bnd_x2,bnd_y1:bnd_y2),    &
                vv2d(bnd_x1:bnd_x2,bnd_y1:bnd_y2),    &
               uup2d(bnd_x1:bnd_x2,bnd_y1:bnd_y2),    &
               vvp2d(bnd_x1:bnd_x2,bnd_y1:bnd_y2),    &
         RHSx2d_tran(bnd_x1:bnd_x2,bnd_y1:bnd_y2),    &  !x-component of external force(baroclinic)
         RHSy2d_tran(bnd_x1:bnd_x2,bnd_y1:bnd_y2) )      !y-component of external force(baroclinic) )

      allocate (stress_t(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),   &      !Horizontal tension tensor component
                stress_s(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),   &      !Horizontal shearing tensor component
                  r_vort(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz) )

      allocate (igrzts_surf(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &     ! igrz[T,S]= 1 : f = f0
                 igrzts_bot(bnd_x1:bnd_x2,bnd_y1:bnd_y2)  )       !          = 2 : df/dz = f0

      allocate (amts(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),     &   !T lateral diffusion in T-points [m^2/s]
                amuv(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),     &   !U and V 2th order lateral diffusion in T-points[m^2/s]
               amuv4(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz)  )       !U and V 4th order lateral viscosity in T-points[m^4/s]^(1/2)

      allocate (slrx(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1),   &   !isopycnal diffusion slope in x-direction
                slry(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1),   &   !isopycnal diffusion slope in y-direction
                slzx(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1),   &   !horizontal diffusion slope in x-direction
                slzy(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1)     )  !horizontal diffusion slope in y-direction

      allocate (rit(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1),     &     !Richardson number
               anzt(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1),     &     !T vertical diffusion [m^2/s]
               anzu(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz+1)      )     !U and V vertical viscosity [m^2/s]

      allocate (tflux_surf(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &       !total surface heat flux [�C*m/s]
                tflux_bot(bnd_x1:bnd_x2,bnd_y1:bnd_y2),       &       !total bottom heat flux [�C*m/s]
                sflux_surf(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &       !total surface salt flux [   m/s]
                sflux_bot(bnd_x1:bnd_x2,bnd_y1:bnd_y2),       &       !total bottom salt flux [   m/s]
            surf_stress_x(bnd_x1:bnd_x2,bnd_y1:bnd_y2),       &       !wind      zonal stress per water density [m^2/s^2]
            surf_stress_y(bnd_x1:bnd_x2,bnd_y1:bnd_y2),       &       !wind meridional stress per water density [m^2/s^2]
             bot_stress_x(bnd_x1:bnd_x2,bnd_y1:bnd_y2),       &       !bottom      zonal stress per water density [m^2/s^2]
             bot_stress_y(bnd_x1:bnd_x2,bnd_y1:bnd_y2),       &       !bottom meridional stress per water density [m^2/s^2]
                 divswrad(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),    &       !shortwave radiation divergence coefficients
                     dkft(bnd_x1:bnd_x2,bnd_y1:bnd_y2),       &       !relaxation coefficient for SST, [m/s]
                     dkfs(bnd_x1:bnd_x2,bnd_y1:bnd_y2),       &       !relaxation coefficient for SSS, [m/s]
                 sensheat(bnd_x1:bnd_x2,bnd_y1:bnd_y2),       &       !sensible heat flux
                  latheat(bnd_x1:bnd_x2,bnd_y1:bnd_y2),       &       !latent heat flux
                   lw_bal(bnd_x1:bnd_x2,bnd_y1:bnd_y2),       &       !longwave radiation balance
                   sw_bal(bnd_x1:bnd_x2,bnd_y1:bnd_y2),       &       !shortwave radiation balance
                   hf_tot(bnd_x1:bnd_x2,bnd_y1:bnd_y2),       &       !total heat flux
                   wf_tot(bnd_x1:bnd_x2,bnd_y1:bnd_y2) )              !total water flux

      allocate (q2_surf(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &     !surface boundary condition for q2
                 q2_bot(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &     !bottom  boundary condition for q2
               q2l_surf(bnd_x1:bnd_x2,bnd_y1:bnd_y2),      &     !surface boundary condition for q2l
                q2l_bot(bnd_x1:bnd_x2,bnd_y1:bnd_y2) )           !bottom  boundary condition for q2l

      allocate (tatm(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &    !Air temperature, [�C]
                qatm(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &    !Air humidity, [kg/kg]
                rain(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &    !rain, [m/s]
                snow(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &    !snow, [m/s]
                wind(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &    !Wind speed module, [m/s]
                 lwr(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &    !Downward  longwave radiation, [W/m^2]
                 swr(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &    !Downward shortwave radiation, [W/m^2]
                slpr(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &    !Sea level pressure, [Pa]
                uwnd(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &    !Zonal      wind speed, [m/s]
                vwnd(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &    !Meridional wind speed, [m/s]
                taux(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &    !Zonal      wind speed, [m/s]
                tauy(bnd_x1:bnd_x2,bnd_y1:bnd_y2) )       !Meridional wind speed, [m/s]


      allocate( BottomFriction(bnd_x1:bnd_x2,bnd_y1:bnd_y2),        &
                       r_diss(bnd_x1:bnd_x2,bnd_y1:bnd_y2) )

      allocate ( hice(bnd_x1:bnd_x2,bnd_y1:bnd_y2,mgrad),    &
                 aice(bnd_x1:bnd_x2,bnd_y1:bnd_y2,mgrad),    &
                aice0(bnd_x1:bnd_x2,bnd_y1:bnd_y2)      ,    &
                hsnow(bnd_x1:bnd_x2,bnd_y1:bnd_y2,mgrad),    &
                 tice(bnd_x1:bnd_x2,bnd_y1:bnd_y2,mgrad),    &
                tsnow(bnd_x1:bnd_x2,bnd_y1:bnd_y2,mgrad),    &
              dhsnowt(bnd_x1:bnd_x2,bnd_y1:bnd_y2),          &
               dhicet(bnd_x1:bnd_x2,bnd_y1:bnd_y2),          &
                swice(bnd_x1:bnd_x2,bnd_y1:bnd_y2),          &
           heatice2oc(bnd_x1:bnd_x2,bnd_y1:bnd_y2),          &
           ice_stress11(bnd_x1:bnd_x2,bnd_y1:bnd_y2),        &
           ice_stress22(bnd_x1:bnd_x2,bnd_y1:bnd_y2),        &
           ice_stress12(bnd_x1:bnd_x2,bnd_y1:bnd_y2),        &
                   uice(bnd_x1:bnd_x2,bnd_y1:bnd_y2),        &
                   vice(bnd_x1:bnd_x2,bnd_y1:bnd_y2)       )

      allocate (Flux_tem_x(bnd_x1:bnd_x2,bnd_y1:bnd_y2,4), &       !Total temperature flux along x-direction
                Flux_tem_y(bnd_x1:bnd_x2,bnd_y1:bnd_y2,4), &       !Total temperature flux along y-direction
                Flux_sal_x(bnd_x1:bnd_x2,bnd_y1:bnd_y2,4), &       !Total   salinity  flux along x-direction
                Flux_sal_y(bnd_x1:bnd_x2,bnd_y1:bnd_y2,4)  )       !Total   salinity  flux along y-direction

      allocate( amuv2d(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &    !depth mean lateral viscosity
               amuv42d(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &    !depth mean lateral viscosity
              r_vort2d(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &    !relative vorticity of depth mean velocity
            stress_t2d(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &    !Horizontal tension tensor component (barotropic)
            stress_s2d(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &    !Horizontal shearing tensor component(barotropic)
           RHSx2d_tran_disp(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &    !dispersion x-component of external force(barotropic)
           RHSy2d_tran_disp(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &    !dispersion x-component of external force(barotropic)
           RHSx2d_diff_disp(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &    !dispersion x-component of external force(barotropic)
           RHSy2d_diff_disp(bnd_x1:bnd_x2,bnd_y1:bnd_y2)  )        !dispersion x-component of external force(barotropic)

      allocate ( tt_calc(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),     &
                 ss_calc(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),     &
                 uu_calc(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),     &
                 vv_calc(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),     &
                sfl_calc(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),     &
                ssh_calc(bnd_x1:bnd_x2,bnd_y1:bnd_y2),        &
                txo_calc(bnd_x1:bnd_x2,bnd_y1:bnd_y2),        &
                tyo_calc(bnd_x1:bnd_x2,bnd_y1:bnd_y2),        &
               uwnd_calc(bnd_x1:bnd_x2,bnd_y1:bnd_y2),        &
               vwnd_calc(bnd_x1:bnd_x2,bnd_y1:bnd_y2),        &
             Fltx_calc(bnd_x1:bnd_x2,bnd_y1:bnd_y2,4),        &
             Flty_calc(bnd_x1:bnd_x2,bnd_y1:bnd_y2,4),        &
             Flsx_calc(bnd_x1:bnd_x2,bnd_y1:bnd_y2,4),        &
             Flsy_calc(bnd_x1:bnd_x2,bnd_y1:bnd_y2,4) )

     ssh_i =0.0d0; sshp_i=0.0d0
     pgrx=0.0d0; pgry=0.0d0
     ubrtr_i=0.0d0;  vbrtr_i=0.0d0
     RHSx2d=0.0d0; RHSy2d=0.0d0

     ssh_e=0.0d0; sshp_e=0.0d0
     ubrtr_e=0.0d0; ubrtrp_e=0.0d0
     vbrtr_e=0.0d0; vbrtrp_e=0.0d0

     uu=0.0d0; uup=0.0d0; vv=0.0d0; vvp=0.0d0; ww=0.0d0
     tt=0.0d0; ttp=0.0d0; ss=0.0d0; ssp=0.0d0
     den=0.0d0; den_pot=0.0d0
     RHSx3d     =0.0d0; RHSy3d     =0.0d0
     RHSx3d_tran=0.0d0; RHSy3d_tran=0.0d0
     age=0.0d0; agep=0.0d0
     xxt=0.0d0; yyt=0.0d0

     q2=0.0d0; q2p=0.0d0; q2l=0.0d0; q2lp=0.0d0
     RHS_q2=0.0d0; RHS_q2l=0.0d0

     uu2d=0.0d0; vv2d=0.0d0; uup2d=0.0d0; vvp2d=0.0d0
     RHSx2d_tran =0.0d0; RHSy2d_tran=0.0d0

     stress_t=0.0d0; stress_s=0.0d0; r_vort=0.0d0

     igrzts_surf=0; igrzts_bot=0

     amts=0.0d0; amuv=0.0d0; amuv4=0.0d0
     slrx=0.0d0; slry=0.0d0; slzx=0.0d0; slzy=0.0d0

     rit=0.0d0; anzt=0.0d0; anzu=0.0d0

     tflux_surf=0.0d0; tflux_bot=0.0d0
     sflux_surf=0.0d0; sflux_bot=0.0d0
     surf_stress_x=0.0d0; surf_stress_y=0.0d0
      bot_stress_x=0.0d0;  bot_stress_y=0.0d0
     divswrad=0.0d0; dkft=0.0d0; dkfs=0.0d0
     sensheat=0.0d0; latheat=0.0d0; lw_bal=0.0d0; sw_bal=0.0d0
     hf_tot=0.0d0; wf_tot=0.0d0

     q2_surf=0.0d0; q2_bot=0.0d0; q2l_surf=0.0d0; q2l_bot=0.0d0

     tatm=0.0d0; qatm=0.0d0; rain=0.0d0; snow=0.0d0; wind=0.0d0
     lwr=0.0d0; swr=0.0d0; slpr=0.0d0; uwnd=0.0d0; vwnd=0.0d0
     taux=0.0d0; tauy=0.0d0
     BottomFriction=0.0d0; r_diss=0.0d0

     hice=0.0d0; aice=0.0d0; aice0=1.0d0; hsnow=0.0d0
     ice_stress11=0.0d0; ice_stress22=0.0d0; ice_stress12=0.0d0
     uice=0.0d0; vice=0.0d0; tice=0.0d0; tsnow=0.0d0
     dhsnowt=0.0d0; dhicet=0.0d0; swice=0.0d0; heatice2oc=0.0d0

     Flux_tem_x=0.0d0; Flux_tem_y=0.0d0; Flux_sal_x=0.0d0; Flux_sal_y=0.0d0

     amuv2d=0.0d0; amuv42d=0.0d0; r_vort2d=0.0d0
     stress_t2d=0.0d0; stress_s2d=0.0d0
     RHSx2d_tran_disp=0.0d0; RHSy2d_tran_disp=0.0d0
     RHSx2d_diff_disp=0.0d0; RHSy2d_diff_disp=0.0d0

       tt_calc=0.0d0;   ss_calc=0.0d0;   uu_calc=0.0d0;   vv_calc=0.0d0
      sfl_calc=0.0d0;  ssh_calc=0.0d0;  txo_calc=0.0d0;  tyo_calc=0.0d0
     uwnd_calc=0.0d0; vwnd_calc=0.0d0
     Fltx_calc=0.0d0; Flty_calc=0.0d0; Flsx_calc=0.0d0; Flsy_calc=0.0d0

     meancalc=0

   endsubroutine ocean_variables_allocate

!deallocation of arrays
   subroutine ocean_variables_deallocate
   use ocean_variables
   implicit none

      deallocate(Flsy_calc, Flsx_calc, Flty_calc, Fltx_calc,      &
                  vwnd_calc, uwnd_calc,                           &
                  tyo_calc,  txo_calc,  ssh_calc,  sfl_calc,      &
                   vv_calc,   uu_calc,   ss_calc,   tt_calc)
      deallocate(RHSy2d_diff_disp,RHSx2d_diff_disp,RHSy2d_tran_disp,RHSx2d_tran_disp,     &
                 stress_s2d,stress_t2d,r_vort2d,amuv42d,amuv2d)
      deallocate(Flux_sal_y, Flux_sal_x, Flux_tem_y, Flux_tem_x)
      deallocate(vice, uice, ice_stress12, ice_stress22, ice_stress11,       &
                 heatice2oc, swice, dhicet, dhsnowt,                         &
                 tsnow, tice, hsnow, aice0, aice, hice)
      deallocate(r_diss, BottomFriction)
      deallocate(tauy,taux,vwnd,uwnd,slpr,swr,lwr,wind,snow,rain,qatm,tatm)

      deallocate(q2l_bot,q2l_surf,q2_bot,q2_surf)

      deallocate(wf_tot,hf_tot,sw_bal,lw_bal,latheat,sensheat,dkfs,dkft,           &
                 divswrad,bot_stress_y,bot_stress_x,surf_stress_y,surf_stress_x,   &
                 sflux_bot,sflux_surf,tflux_bot,tflux_surf)
      deallocate(anzu,anzt,rit,slzy,slzx,slry,slrx,amuv4,amuv,amts)
      deallocate(igrzts_bot,igrzts_surf)
      deallocate(r_vort,stress_s, stress_t)
      deallocate(RHSy2d_tran,RHSx2d_tran,vvp2d,uup2d,vv2d,uu2d)

      deallocate(RHS_q2l,RHS_q2)
      deallocate(q2lp,q2l,q2p,q2)

      deallocate(yyt,xxt,agep,age,RHSy3d_tran,RHSx3d_tran,RHSy3d,RHSx3d,den_pot,den,     &
                 ssp,ss,ttp,tt,ww,vvp,vv,uup,uu)
      deallocate(vbrtrp_e,vbrtr_e,ubrtrp_e,ubrtr_e,sshp_e,ssh_e)
      deallocate(RHSy2d,RHSx2d,vbrtr_i,ubrtr_i,pgry,pgrx,sshp_i,ssh_i)

   endsubroutine ocean_variables_deallocate
   !Allocation of arrays
!---------------------------------------------------------------------------------
   subroutine pass_tracer_allocate
   use main_basin_pars
   use mpi_parallel_tools
   use ocean_variables
   implicit none

      allocate (pass_tracer(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz), &     !Passive tracer
               pass_tracerp(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz), &     !Passive tracer
                pt_forc_surf(bnd_x1:bnd_x2,bnd_y1:bnd_y2),   &     !Passive tracer surface forcing
                pt_forc_bot(bnd_x1:bnd_x2,bnd_y1:bnd_y2),    &     !Passive tracer bottom forcing
                pt_diff_x(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),   &     !Passive x-diffusion [m^2/s]
                pt_diff_y(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz),   &     !Passive y-diffusion [m^2/s]
                igrzpt_surf(bnd_x1:bnd_x2,bnd_y1:bnd_y2),    &     !Types of boundary condition
                igrzpt_bot(bnd_x1:bnd_x2,bnd_y1:bnd_y2)      )

          pass_tracer=0.0d0; pass_tracerp=0.0d0
          pt_forc_surf=0.0d0; pt_forc_bot=0.0d0
          pt_diff_x=0.0d0; pt_diff_y=0.0d0
          igrzpt_surf=0; igrzpt_bot=0

   endsubroutine pass_tracer_allocate

!deallocation of arrays
   subroutine pass_tracer_deallocate
   use ocean_variables
   implicit none

      deallocate(igrzpt_bot,igrzpt_surf,pt_diff_y,pt_diff_x,pt_forc_bot,pt_forc_surf,pass_tracer)

   endsubroutine pass_tracer_deallocate
!========================================================================
  subroutine oceanbc_arrays_allocate
  use main_basin_pars
  use mpi_parallel_tools
  use ocean_bc
  implicit none

!   allocate (sst_obs(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &       !Observed SST [�C]
!             sss_obs(bnd_x1:bnd_x2,bnd_y1:bnd_y2),     &       !Observed SSS [psu-35ppt])
!              runoff(bnd_x1:bnd_x2,bnd_y1:bnd_y2)   )          !River runoff, [m/s]

   allocate (sst_obs(nx, ny),     &       !Observed SST [�C]
             sss_obs(nx, ny),     &       !Observed SSS [psu-35ppt])
              runoff(nx, ny)   )          !River runoff, [m/s]


   allocate(lqpx(numb_of_lqp_max),        &
            lqpy(numb_of_lqp_max),        &
            tlqbw(numb_of_lqp_max,nz),    &
            slqbw(numb_of_lqp_max,nz),    &
            ulqbw(numb_of_lqp_max,nz),    &
            vlqbw(numb_of_lqp_max,nz),    &
          sshlqbw(numb_of_lqp_max))

   allocate(index_of_lb(numb_of_lqp_max))

      sst_obs=0.0; sss_obs=0.0; runoff=0.0
      lqpx=0; lqpy=0; tlqbw=0.0; slqbw=0.0
      ulqbw=0.0; vlqbw=0.0; sshlqbw=0.0
      index_of_lb=' '
  endsubroutine oceanbc_arrays_allocate

  subroutine oceanbc_arrays_deallocate
  use ocean_bc
  implicit none
   deallocate(index_of_lb)
   deallocate(sshlqbw,vlqbw,ulqbw,slqbw,tlqbw,lqpy,lqpx)
   deallocate(runoff, sss_obs, sst_obs)
  endsubroutine oceanbc_arrays_deallocate
!--------------------------------------------------------------------------------------
subroutine atm_arrays_allocate
use atm_pars
use atm_forcing
implicit none

 allocate(xa(nxa),ya(nya))

 allocate( a_hflux(nxa,nya),       &   !heat balance [w/m**2]
           a_swrad(nxa,nya),       &   !sw radiation balance[w/m**2]
           a_wflux(nxa,nya),       &   !precipitation-evaporation[m/s]
           a_stress_x(nxa,nya),    &   !zonal wind stress[pA=n/m**2]
           a_stress_y(nxa,nya),    &   !meridional wind stress[pA=n/m**2]
           a_slpr(nxa,nya),        &   !pressure at sea surface
           a_lwr(nxa,nya),         &   !dw-lw-rad[w/m**2]
           a_swr(nxa,nya),         &   !dw-sw-rad[w/m**2]
           a_rain(nxa,nya),        &   !precipit[m/s]
           a_snow(nxa,nya),        &   !precipit[m/s]
           a_tatm(nxa,nya),        &   !temp of atmosphere[�c]
           a_qatm(nxa,nya),        &   !humidity [g/kg]
           a_uwnd(nxa,nya),        &   !u-wind speed[m/s]
           a_vwnd(nxa,nya)  )          !v-wind speed[m/s]

      xa=0.0d0; ya=0.0d0

      a_hflux=0.0; a_swrad=0.0; a_wflux=0.0; a_stress_x=0.0; a_stress_y=0.0
      a_slpr=0.0;  a_lwr=0.0;   a_swr=0.0; a_rain=0.0; a_snow=0.0
      a_tatm=0.0;  a_qatm=0.0; a_uwnd=0.0; a_vwnd=0.0

        ind_change_heat =0
        ind_change_water=0
        ind_change_stress=0
        ind_change_rad  =0
        ind_change_prec =0
        ind_change_tatm =0
        ind_change_qatm =0
        ind_change_wind =0
        ind_change_slpr =0

        num_rec_heat =0
        num_rec_water=0
        num_rec_stress=0
        num_rec_rad  =0
        num_rec_prec =0
        num_rec_tatm =0
        num_rec_qatm =0
        num_rec_wind =0
        num_rec_slpr =0

endsubroutine atm_arrays_allocate
!-------------------------------------------------------------------------------
subroutine atm_arrays_deallocate
use atm_pars
use atm_forcing
implicit none

 deallocate(a_vwnd,a_uwnd,a_qatm,a_tatm,a_snow,a_rain,a_swr,a_lwr,     &
             a_slpr,a_stress_y,a_stress_x,a_wflux,a_swrad,a_hflux)
 deallocate(ya,xa)

endsubroutine atm_arrays_deallocate
!=====================================================================================
 subroutine atm2oc_allocate
 use mpi_parallel_tools
 use atm2oc_interpol
 implicit none

      allocate(   wght_mtrx_a2o(bnd_x1:bnd_x2,bnd_y1:bnd_y2,4),      &
                    i_input_a2o(bnd_x1:bnd_x2,bnd_y1:bnd_y2,4),      &
                    j_input_a2o(bnd_x1:bnd_x2,bnd_y1:bnd_y2,4)    )
      wght_mtrx_a2o=0.0d0; i_input_a2o=0; j_input_a2o=0

 endsubroutine atm2oc_allocate
!------------------------------------------------------------------------------------
 subroutine atm2oc_deallocate
 use atm2oc_interpol
 implicit none
      deallocate(j_input_a2o  , i_input_a2o  , wght_mtrx_a2o  )
 endsubroutine atm2oc_deallocate
