!=====================================================
subroutine cyclize_x(ff,nx,ny,nz,mmm,mm)
    use mpi_parallel_tools
    implicit none
!---------------------------------------------------------------------
! adds periodically left (m=mmm-1) and right (m=mm+1) for cyclic lines
    integer :: nx, ny, nz
    integer :: mmm, mm, n, k
    real(4) :: ff(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz)

    integer, dimension(2) :: p_dist
    integer :: dist_rank
    integer :: ierr
    integer stat(MPI_status_size)

    if (p_coord(1) .eq. 0) then
!-------------- proc has mmm-1 area --------------------------------------------
        p_dist(1) = p_size(1) - 1
        p_dist(2) = p_coord(2)
        call mpi_cart_rank(cart_comm, p_dist,dist_rank,ierr)

        call mpi_sendrecv(ff(mmm, bnd_y1:bnd_y2, 1:nz),         &
                          (bnd_y2 - bnd_y1 + 1)*nz,              &
                          mpi_real4, dist_rank, 1,          &
                          ff(mmm-1, bnd_y1:bnd_y2, 1:nz),       &
                          (bnd_y2 - bnd_y1 + 1)*nz,              &
                          mpi_real4, dist_rank, 1,          &
                          cart_comm, stat, ierr)
!                write(*,*)rank,p_coord,'0',ierr
    endif

    if (p_coord(1) .eq. (p_size(1) - 1)) then
!-------------- proc has mm+1 area ---------------------------------------------
        p_dist(1) = 0
        p_dist(2) = p_coord(2)
        call mpi_cart_rank(cart_comm, p_dist,dist_rank,ierr)

        call mpi_sendrecv(ff(mm, bnd_y1:bnd_y2, 1:nz),          &
                          (bnd_y2 - bnd_y1 + 1)*nz,              &
                          mpi_real4, dist_rank, 1,          &
                          ff(mm+1, bnd_y1:bnd_y2, 1:nz),        &
                          (bnd_y2 - bnd_y1 + 1)*nz,              &
                          mpi_real4, dist_rank, 1,          &
                          cart_comm, stat, ierr)
!                write(*,*)rank,p_coord,'size-1',ierr
    endif
endsubroutine cyclize_x

!======================================================
subroutine cyclize8_x(ff,nx,ny,nz,mmm,mm)
    use mpi_parallel_tools
    implicit none
!---------------------------------------------------------------------
! adds periodically left (m=mmm-1) and right (m=mm+1) for cyclic lines
    integer :: nx, ny, nz
    integer :: mmm, mm, n, k
    real(8) :: ff(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz)

    integer, dimension(2) :: p_dist
    integer :: dist_rank
    integer :: ierr
    integer stat(MPI_status_size)

    if (p_coord(1) .eq. 0) then
!-------------- proc has mmm-1 area --------------------------------------------
        p_dist(1) = p_size(1) - 1
        p_dist(2) = p_coord(2)
        call mpi_cart_rank(cart_comm, p_dist,dist_rank,ierr)

        call mpi_sendrecv(ff(mmm, bnd_y1:bnd_y2, 1:nz),         &
                          (bnd_y2 - bnd_y1 + 1)*nz,              &
                          mpi_real8, dist_rank, 1,          &
                          ff(mmm-1, bnd_y1:bnd_y2, 1:nz),       &
                          (bnd_y2 - bnd_y1 + 1)*nz,              &
                          mpi_real8, dist_rank, 1,          &
                          cart_comm, stat, ierr)
!                write(*,*)rank,p_coord,'0',ierr
    endif

    if (p_coord(1) .eq. (p_size(1) - 1)) then
!-------------- proc has mm+1 area ---------------------------------------------
        p_dist(1) = 0
        p_dist(2) = p_coord(2)
        call mpi_cart_rank(cart_comm, p_dist,dist_rank,ierr)

        call mpi_sendrecv(ff(mm, bnd_y1:bnd_y2, 1:nz),          &
                          (bnd_y2 - bnd_y1 + 1)*nz,              &
                          mpi_real8, dist_rank, 1,          &
                          ff(mm+1, bnd_y1:bnd_y2, 1:nz),        &
                          (bnd_y2 - bnd_y1 + 1)*nz,              &
                          mpi_real8, dist_rank, 1,          &
                          cart_comm, stat, ierr)
!                write(*,*)rank,p_coord,'size-1',ierr
    endif
endsubroutine cyclize8_x

!=====================================================
subroutine cyclize_y(ff,nx,ny,nz,nnn,nn)
implicit none
!---------------------------------------------------------------------
! adds periodically bottom (n=nnn-1) and top (n=nn+1) for cyclic lines
 integer nx, ny, nz
 integer nnn, nn, m, k
 real(4) ff(nx,ny,nz)

!$omp parallel do private(m,k)
  do m=1,nx
   do k=1,nz
      ff(m,nnn-1,k) = ff(m,nn ,k)
      ff(m, nn+1,k) = ff(m,nnn,k)
   enddo
  enddo
!$omp end parallel do

endsubroutine cyclize_y

!=====================================================
subroutine cyclize8_y(ff,nx,ny,nz,nnn,nn)
implicit none
!---------------------------------------------------------------------
! adds periodically bottom (n=nnn-1) and top (n=nn+1) for cyclic lines
 integer nx, ny, nz
 integer nnn, nn, m, k
 real(8) ff(nx,ny,nz)

!$omp parallel do private(m,k)
  do m=1,nx
   do k=1,nz
      ff(m,nnn-1,k) = ff(m,nn ,k)
      ff(m, nn+1,k) = ff(m,nnn,k)
   enddo
  enddo
!$omp end parallel do

endsubroutine cyclize8_y
