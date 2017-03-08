!-----------module for definition of array dimensions and boundaries------------
module mpi_parallel_tools
    implicit none

    include 'mpif.h'

    integer :: nx_start,    &   !first significant point in x-direction
               nx_end,      &   !last  significant point in x-direction
               ny_start,    &   !first significant point in y-direction
               ny_end           !last  significant point in y-direction

    integer :: len_x1,     &   !length of array boundary
               len_x2,     &
               len_y1,     &
               len_y2

    integer :: bnd_x1,      &   !left   array boundary in x-direction
               bnd_x2,      &   !right  array boundary in x-direction
               bnd_y1,      &   !bottom array boundary in y-direction
               bnd_y2           !top    array boundary in y-direction

    integer :: rank, procs
    integer :: cart_comm
    integer, dimension(2) :: p_size, period, p_coord

    contains

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
    integer function check_p_coord(coord)
        implicit none
        integer, dimension(2), intent(in) :: coord

        check_p_coord = 0
!            write(*,*) coord,all(coord.ge.0),all((p_size-coord).ge.1)
!            print *, coord, p_size - coord, all((p_size-coord).ge.1)
        if (all(coord.ge.0) .and. all((p_size-coord).ge.1)) then
            check_p_coord = 1
        endif
        return
    end function

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
    subroutine directsync_real8(field, p_dist, xd1, xd2, yd1, yd2,             &
                                       p_src,  xs1, xs2, ys1, ys2)
        implicit none
        real*8, intent(in out) :: field(bnd_x1:bnd_x2, bnd_y1:bnd_y2)
        integer, dimension(2), intent(in) :: p_dist, p_src
        integer :: xd1, xd2, yd1, yd2 ! bound of array which sending to p_dist
        integer :: xs1, xs2, ys1, ys2 ! bound of array which recieving from p_src

        integer :: dist_rank, src_rank
        integer :: flag_dist, flag_src
        integer :: ierr, debg
        integer :: stat(mpi_status_size)

        debg = 0

        if ( ((xd1-xd2+1)*(yd1-yd2+1)) .ne. (xs1-xs2+1)*(ys1-ys2+1) ) then
            print *, "Error in sync arrays size!"
        endif

        flag_dist = check_p_coord(p_dist)
        flag_src = check_p_coord(p_src)

        if ( (flag_src .eq. 1) .and. (flag_dist .eq. 1) ) then
            call mpi_cart_rank(cart_comm, p_dist,dist_rank,ierr)
            call mpi_cart_rank(cart_comm, p_src, src_rank, ierr)

            call mpi_sendrecv(field(xd1:xd2, yd1:yd2),                          &
                              (xd2 - xd1 + 1)*(yd2 - yd1 + 1),                 &
                              mpi_real8, dist_rank, 1,                         &
                              field(xs1:xs2, ys1:ys2),                          &
                              (xs2 - xs1 + 1)*(ys2 - ys1 + 1),                 &
                              mpi_real8, src_rank, 1,                          &
                              cart_comm, stat, ierr)
!            print *, rank, "rsendecv", ierr
        else
            if (flag_src .eq. 1) then
                call mpi_cart_rank(cart_comm,p_src,src_rank,ierr)

                call mpi_recv(field(xs1:xs2, ys1:ys2),                          &
                              (xs2 - xs1 + 1)*(ys2 - ys1 + 1),                 &
                              mpi_real8, src_rank, 1,                          &
                              cart_comm, stat, ierr)
!                print *, rank, src_rank, "recv", xs1, xs2, ys1, ys2, field(xs1:xs2, ys1:ys2)
            endif

            if (flag_dist .eq. 1) then
                call mpi_cart_rank(cart_comm,p_dist,dist_rank,ierr)

                call mpi_send(field(xd1:xd2, yd1:yd2),                          &
                             (xd2 - xd1 + 1)*(yd2 - yd1 + 1),                  &
                             mpi_real8, dist_rank, 1,                          &
                             cart_comm, stat, ierr)
!                print *, rank, dist_rank, "send", xd1, xd2, yd1, yd2, field(xd1:xd2, yd1:yd2)
            endif
        endif

    end subroutine directsync_real8

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
    subroutine syncborder_real8(field)
        implicit none
        real*8, intent(in out) :: field(bnd_x1:bnd_x2, bnd_y1:bnd_y2)

        integer, dimension(2) :: p_dist, p_src

!------------------ send-recv in ny+ -------------------------------------------
        p_dist(1) = p_coord(1)
        p_dist(2) = p_coord(2) + 1
        p_src(1) = p_coord(1)
        p_src(2) = p_coord(2) - 1
        call directsync_real8(field, p_dist, nx_start, nx_end, ny_end - len_y1 + 1, ny_end,       &
                                     p_src,  nx_start, nx_end, bnd_y1, bnd_y1 + len_y1 - 1)
!------------------ send-recv in nx+ -------------------------------------------
        p_dist(1) = p_coord(1) + 1
        p_dist(2) = p_coord(2)
        p_src(1) = p_coord(1) - 1
        p_src(2) = p_coord(2)
        call directsync_real8(field, p_dist, nx_end - len_x1 + 1, nx_end, ny_start, ny_end,       &
                                     p_src,  bnd_x1, bnd_x1 + len_x1 - 1, ny_start, ny_end)
!------------------ send-recv in ny- -------------------------------------------
        p_dist(1) = p_coord(1)
        p_dist(2) = p_coord(2) - 1
        p_src(1) = p_coord(1)
        p_src(2) = p_coord(2) + 1
        call directsync_real8(field, p_dist, nx_start, nx_end, ny_start, ny_start + len_y1 - 1,   &
                                     p_src,  nx_start, nx_end, bnd_y2 - len_y2 + 1, bnd_y2)
!------------------ send-recv in nx- -------------------------------------------
        p_dist(1) = p_coord(1) - 1
        p_dist(2) = p_coord(2)
        p_src(1) = p_coord(1) + 1
        p_src(2) = p_coord(2)
        call directsync_real8(field, p_dist, nx_start, nx_start + len_x2 - 1, ny_start, ny_end,   &
                                     p_src,  bnd_x2 - len_x2 + 1, bnd_x2, ny_start, ny_end)


!------------------ Sync edge points (EP) --------------------------------------
!------------------ send-recv EP in nx+,ny+ ------------------------------------
         p_dist(1) = p_coord(1) + 1
         p_dist(2) = p_coord(2) + 1
         p_src(1) = p_coord(1) - 1
         p_src(2) = p_coord(2) - 1
         call directsync_real8(field, p_dist, nx_end - len_x1 + 1, nx_end, ny_end - len_y1 + 1, ny_end,   &
                                      p_src,  bnd_x1, bnd_x1 + len_x1 - 1, bnd_y1, bnd_y1 + len_y1 - 1)
!------------------ send-recv EP in nx+,ny- ------------------------------------
         p_dist(1) = p_coord(1) + 1
         p_dist(2) = p_coord(2) - 1
         p_src(1) = p_coord(1) - 1
         p_src(2) = p_coord(2) + 1
         call directsync_real8(field, p_dist, nx_end - len_x1 + 1, nx_end, ny_start, ny_start + len_y2 - 1,   &
                                      p_src,  bnd_x1, bnd_x1 + len_x1 - 1, bnd_y2 - len_y2 + 1 , bnd_y2)
!------------------ send-recv EP in nx-,ny- ------------------------------------
         p_dist(1) = p_coord(1) - 1
         p_dist(2) = p_coord(2) - 1
         p_src(1) = p_coord(1) + 1
         p_src(2) = p_coord(2) + 1
         call directsync_real8(field, p_dist, nx_start, nx_start + len_x2 - 1, ny_start, ny_start + len_y2 - 1,   &
                                      p_src,  bnd_x2 - len_x2 + 1, bnd_x2, bnd_y2 - len_y2 + 1, bnd_y2)

!------------------ send-recv EP in nx-,ny+ ------------------------------------
         p_dist(1) = p_coord(1) - 1
         p_dist(2) = p_coord(2) + 1
         p_src(1) = p_coord(1) + 1
         p_src(2) = p_coord(2) - 1
         call directsync_real8(field, p_dist, nx_start, nx_start + len_x2 - 1, ny_end - len_y1 + 1, ny_end,  &
                                      p_src,  bnd_x2 - len_x2 + 1, bnd_x2, bnd_y1, bnd_y1 + len_y1 - 1)

        return

    end subroutine syncborder_real8


endmodule mpi_parallel_tools
!---------------------end module for definition of array dimensions and boundaries-----------------
