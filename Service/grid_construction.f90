!======================================================================
! grid consruction module by temperature mask.
subroutine gridcon(ftemask)
use main_basin_pars
use mpi_parallel_tools
use basin_grid
use rec_length
implicit none
! subroutin for construction pass boundary, velosity and bottom masks
! using temperature mask in diogin standart
!--------------------------------------------------------------------
character*(*) ftemask
character frmt*16,comment*80
! temporary integer indexes
integer m, n, ierr

write(frmt,1000) nx
1000  format('(',i9,'i1)')

! reading mask from:
open (11,file=ftemask,status='old',recl=nx*lrecl)
 read (11,  '(a)') comment(1:min(80,nx))
 if (rank .eq. 0) write(*,'(1x,a)') comment
 do n=ny,1,-1
  read(11,frmt,end=99) (lbasins(m,n),m=1,nx)
 enddo
close(11)

! conversion integer diogin mask to real model mask
do n=bnd_y1,bnd_y2
   do m=bnd_x1,bnd_x2
    if(lbasins(m,n)==0) then
     lu(m,n)=1.0
    endif
   end do
end do


     lu1=1.0

if(periodicity_x/=0) then
  call cyclize_x(lu,nx,ny,1,mmm,mm)
endif

if(periodicity_y/=0) then
  call cyclize_y(lu,nx,ny,1,nnn,nn)
endif

!  forming mask for depth grid points
!  forming luh from lu, which have land neibours in luh.
! constructing array luh for relief hh.

if (rank .eq. 0) then
   write(*,*) 'Construction of H-grid masks: '
   write(*,*) 'LUH (includes boundary) and LUU (does not include boundary)'
endif
do n=ny_start-1,ny_end
   do m=nx_start-1,nx_end
     if(lu(m,n)+lu(m+1,n)+lu(m,n+1)+lu(m+1,n+1)>0.5)  then
      luh(m,n)=1.0
     endif
   enddo
enddo

do n=ny_start-1,ny_end
   do m=nx_start-1,nx_end
     if(lu(m,n)*lu(m+1,n)*lu(m,n+1)*lu(m+1,n+1)>0.5)  then
      luu(m,n)=1.0
     endif
   enddo
enddo

call syncborder_real(luh)
call syncborder_real(luu)

if(periodicity_x/=0) then
  call cyclize_x(luh,nx,ny,1,mmm,mm)
  call cyclize_x(luu,nx,ny,1,mmm,mm)
endif

if(periodicity_y/=0) then
  call cyclize_y(luh,nx,ny,1,nnn,nn)
  call cyclize_y(luu,nx,ny,1,nnn,nn)
endif

if (rank .eq. 0) then
    write(*,*) 'Construction of U- and V-grid masks: '
    write(*,*) 'LCU and LCV (do not include boundary) and LLU and LLV (include boundary)'
endif
do n=ny_start-1,ny_end
   do m=nx_start-1,nx_end

    if(lu(m,n)+lu(m+1,n)>0.5)  then
     llu(m,n)=1.0
    endif

    if(lu(m,n)+lu(m,n+1)>0.5)  then
     llv(m,n)=1.0
    endif

    if(lu(m,n)*lu(m+1,n)>0.5)  then
     lcu(m,n)=1.0
    endif

    if(lu(m,n)*lu(m,n+1)>0.5)  then
     lcv(m,n)=1.0
    endif

   enddo
enddo

call syncborder_real(llu)
call syncborder_real(llv)
call syncborder_real(lcu)
call syncborder_real(lcv)

if (periodicity_x/=0) then
 if (rank .eq. 0) write(*,*)'  set periodicity to u-grid mask(lcu,llu).'
 call cyclize_x(lcu,nx,ny,1,mmm,mm)
 call cyclize_x(llu,nx,ny,1,mmm,mm)
 if (rank .eq. 0) write(*,*)'  set periodicity to v-grid mask(lcv,llv).'
 call cyclize_x(lcv,nx,ny,1,mmm,mm)
 call cyclize_x(llv,nx,ny,1,mmm,mm)
endif

if (periodicity_y/=0) then
 call cyclize_y(lcu,nx,ny,1,nnn,nn)
 call cyclize_y(llu,nx,ny,1,nnn,nn)
 call cyclize_y(lcv,nx,ny,1,nnn,nn)
 call cyclize_y(llv,nx,ny,1,nnn,nn)
endif

return
99    write(*,*)'  error in reading file ',ftemask(1:len_trim(ftemask))
stop 1
endsubroutine gridcon
