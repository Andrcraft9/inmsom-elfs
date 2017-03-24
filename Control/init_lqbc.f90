!======================================================================
	subroutine lqpcoordinates(path2oceandata)    !path to ocean data
! (c) diansky n.a. DINAR@INM.RAS.RU, 2008
      use mpi_parallel_tools
      use main_basin_pars
      use ocean_bc
      use rec_length

	implicit none
	character*(*) path2oceandata !path to ocean data
	character(1024) name_of_lqpmask   !name of file with lqpmask coded as
	                                  !english letters

!------working variables--------------------------------------------------------
	integer m,n,l,ierr, mark
    integer numb_of_lqp_real

	character frmt*16,comment*1024

	character*1 alphabet(4)
	data alphabet/'N','S','W','E'/

	character, allocatable:: lqbmask(:,:)

!--------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------!
    if (rank .eq. 0) print *, "Error! Lqpcoordinates is not done yet."
    call mpi_finalize(ierr)
    stop
!--------------------------------------------------------------------------------------!
!--------------------------------------------------------------------------------------!

	allocate (lqbmask(nx,ny))

      write(frmt,1000) nx
1000  format('(',i9,'a1)')

! full file name with liquid boundary temperature
      call fulfname(name_of_lqpmask,path2oceandata,'LQBMASK.TXT',ierr)
      if (rank .eq. 0) then
          write(*,*) ' file with mask with liquid boundary pointes marked'
          write(*,'(5x,a)') name_of_lqpmask(1:len_trim (name_of_lqpmask))
      endif


! reading diogin mask from:
      open (11,file=name_of_lqpmask,status='old',recl=nx*lrecl)
      read (11,  '(a)') comment(1:min(1024,nx))
      if (rank .eq. 0) write(*,'(1x,a)') comment(1:80)

      do n=ny,1,-1
         read(11,frmt,end=99) (lqbmask(m,n),m=1,nx)
      enddo

	numb_of_lqp=0
      numb_of_lqp_real=0

! Computing number of open boundares (0 - 4)
	do l=1,4

         mark=0

          do n=1,ny
	     do m=1,nx

            if(lqbmask(m,n)==alphabet(l)) then
              mark=mark+1
		end if

           end do
          end do

         if(mark>0) then
          numb_of_lb=numb_of_lb+1
         endif

	end do

! Computing number of open boundary points
      do l=1,4

         do n=1,ny
	     do m=1,nx

		if(lqbmask(m,n)==alphabet(l)) then

              numb_of_lqp_real=numb_of_lqp_real+1
              lqpx(numb_of_lqp_real)=m
              lqpy(numb_of_lqp_real)=n
              index_of_lb(numb_of_lqp_real)=alphabet(l)

            end if

           end do
         end do

	end do

	numb_of_lqp=numb_of_lqp_real

    if (rank .eq. 0) then
	       write(*,'(a,i4)') ' number of liquid boundaries =',numb_of_lb
	       write(*,'(a,i4)') ' number of liquid pointes    =',numb_of_lqp
    endif

	deallocate (lqbmask)
	return
99    write(*,*)'  error in reading file ', name_of_lqpmask(1:len_trim (name_of_lqpmask))
	stop 1
  end subroutine lqpcoordinates
