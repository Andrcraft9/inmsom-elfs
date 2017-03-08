!---------------------module for definition of basin grid arrays-----------------
module basin_grid
implicit none

!arrays of basin grids
real(4), allocatable:: lu(:,:),        &  !mask of t-grid
                       lu1(:,:),       &  !mask of t-grid (1 everywhere)
                       luu(:,:),       &  !mask of h-grid (0 on boundary)
                       luh(:,:),       &  !mask of h-grid (1 on boundary)
                       lcu(:,:),       &  !mask of u-grid (0 on boundary)
                       lcv(:,:),       &  !mask of v-grid (0 on boundary)
                       llu(:,:),       &  !mask of u-grid (1 on boundary)
                       llv(:,:)           !mask of v-grid (1 on boundary)

integer, allocatable:: lbasins(:,:)       !integer mask

real(8), allocatable:: hhh(:,:),      &  !ocean depth on luh (h-points)
                       hhhp(:,:),     &  !ocean depth on luh (h-points) at previous step
                       hhhn(:,:),     &  !ocean depth on luh (h-points) at pre-previous step
                   hhq_rest(:,:),     &  !ocean depth on lu  (t-points) at rest state
                       hhq(:,:),      &  !ocean depth on lu  (t-points) 
                       hhqp(:,:),     &  !ocean depth on lu  (t-points) at previous step
                       hhqn(:,:),     &  !ocean depth on lu  (t-points) at pre-previous step
                       hhu(:,:),      &  !ocean depth on lcu (u-points)
                       hhup(:,:),     &  !ocean depth on lcu (u-points) at previous step
                       hhun(:,:),     &  !ocean depth on lcu (u-points) at pre-previous step
                       hhv(:,:),      &  !ocean depth on lcv (v-points)
                       hhvp(:,:),     &  !ocean depth on lcv (v-points) at pre-previous step
                       hhvn(:,:),     &  !ocean depth on lcv (v-points) at pre-previous step
                       rlh_s(:,:),    &  !main (sin) coriolis parameter on h-points
                       rlh_c(:,:),    &  !2-nd (cos) coriolis parameter on h-points
                    z(:), zw(:),      &  !vertical sigma-levels (t-points and w-points)
                  hzt(:), dz(:),      &  !steps between t-levels and w-levels
            dxt(:,:),dyt(:,:),        &  !horizontal grid steps between   t-points (in meters)
            dx (:,:),dy (:,:),        &  !horizontal grid steps between u,v-points (in meters)
            dxh(:,:),dyh(:,:),        &  !horizontal grid steps between   h-points (in meters)
            dxb(:,:),dyb(:,:)            !horizontal grid steps between v,u-points (in meters)
 real(8), allocatable::     xt(:),yt(:),        &  !horizontal t-grid            x- and y-coordinates (in degrees)
                            xu(:),yv(:)            !horizontal u-grid and v-grid x- and y-coordinates (in degrees)


 real(8), allocatable::   geo_lon_t(:,:),   &    !geographical longitudes of T-points
                          geo_lat_t(:,:),   &    !geographical latitudes  of T-points
                          geo_lon_u(:,:),   &    !geographical longitudes of U-points
                          geo_lat_u(:,:),   &    !geographical latitudes  of U-points
                          geo_lon_v(:,:),   &    !geographical longitudes of V-points
                          geo_lat_v(:,:),   &    !geographical latitudes  of V-points
                          geo_lon_h(:,:),   &    !geographical longitudes of H-points
                          geo_lat_h(:,:),   &    !geographical latitudes  of H-points
                       rotvec_coeff(:,:,:)       !cos and sin of angles between coordinate lines
 
 real(8), allocatable:: hhh_e(:,:),      &  !ocean depth on luh (h-points)
                        hhhp_e(:,:),     &  !ocean depth on luh (h-points) at previous step
                        hhhn_e(:,:),     &  !ocean depth on luh (h-points) at pre-previous step
                        hhq_e(:,:),      &  !ocean depth on lu  (t-points) 
                        hhqp_e(:,:),     &  !ocean depth on lu  (t-points) at previous step
                        hhqn_e(:,:),     &  !ocean depth on lu  (t-points) at pre-previous step
                        hhu_e(:,:),      &  !ocean depth on lcu (u-points)
                        hhup_e(:,:),     &  !ocean depth on lcu (u-points) at previous step
                        hhun_e(:,:),     &  !ocean depth on lcu (u-points) at pre-previous step
                        hhv_e(:,:),      &  !ocean depth on lcv (v-points)
                        hhvp_e(:,:),     &  !ocean depth on lcv (v-points) at pre-previous step
                        hhvn_e(:,:)         !ocean depth on lcv (v-points) at pre-previous step

endmodule basin_grid
!---------------------end module for definition of basin grid arrays-----------------
