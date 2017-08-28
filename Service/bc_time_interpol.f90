!======================================================================
! atmospheric data time interpolation on atmospheric grid
subroutine atm_data_time_interpol
use atm_pars
use atm_forcing
use time_integration
use key_switches
implicit none

! Forming full file names once a year

   if(m_time_changed(6)>0) then ! if new year comes
       !  forming full file name for freshwater balance
       call bc_full_name(ss_atmfiles(5),          & ! original file name
                         ss_atmfiles_fname(5),    & ! full file name
                              path2atmssdata,     & ! path to flux data
                                 ind_clim_water,  & ! index of climaticity (0-climatic data, 1 - real year data)
                                       m_year)      ! current year number


       !  forming full file name for SLP
       call bc_full_name(ss_atmfiles(10),         & ! original file name
                         ss_atmfiles_fname(10),   & ! full file name
                              path2atmssdata,     & ! path to flux data
                                 ind_clim_slpr,   & ! index of climaticity (0-climatic data, 1 - real year data)
                                       m_year)      ! current year number

   endif


!indicating of atm array change and reading data from files

    if(iabs(ksw_ssbc)==2) then

       ! fresh water balance
       call check_bc_field_change(time_resolution_water,    &
                                            ftype_water,    &
                                        m_month_of_year,    &
                                          m_day_of_year,    &
                                          m_hour_of_day,    &
                                          num_rec_water,    &
                                       ind_change_water    )
       if(ind_change_water>0) then
           call bc_data_read(ss_atmfiles_fname(5),      &
                                          a_wflux,      &
                                              nxa,      &
                                              nya,      &
                                            1,nxa,      &
                                            1,nya,      &
                                      num_rec_water  )

          ! converting fresh water balance to model units (kg/m**2/s)
          a_wflux=a_wflux*confac_water
       endif

	end if

!-------------------------------------------------------------------------------
    if(iabs(ksw_ssbc)>=3) then
       ! SLP
       call check_bc_field_change(time_resolution_slpr,    &
                                            ftype_slpr,    &
                                       m_month_of_year,    &
                                         m_day_of_year,    &
                                         m_hour_of_day,    &
                                          num_rec_slpr,    &
                                       ind_change_slpr    )
       if(ind_change_slpr>0) then
           call bc_data_read(ss_atmfiles_fname(10),     &
                                           a_slpr,      &
                                              nxa,      &
                                              nya,      &
                                            1,nxa,      &
                                            1,nya,      &
                                     num_rec_slpr  )
          ! converting SLP to model units (Pa)
          a_slpr=a_slpr*confac_slpr
       endif

	end if

endsubroutine atm_data_time_interpol

!===============================================================================
! subroutine for forming full file name for atmospheric forcing
subroutine bc_full_name(filename,    & ! original file name
                    full_filename,    & ! full file name
                        path2flux,    & ! path to flux data
                         ind_clim,    & ! index of climaticity (0-climatic data, 1 - real year data)
                           m_year)      ! current year number

    implicit none

    character*(*) path2flux,filename,full_filename
    character fname*64,chyear*4
    integer ind_clim, m_year, m_time_changed(7), ierr

    if(ind_clim==0) then
        call fulfname(fname,'CLIM',filename,ierr)
        call fulfname(full_filename,path2flux,fname,ierr)
    else
        write(chyear(1:4),'(i4.4)') m_year
        call fulfname(fname,chyear,filename,ierr)
        call fulfname(full_filename,path2flux,fname,ierr)
    end if

endsubroutine bc_full_name
!-----------------------------------------------------------
subroutine check_bc_field_change(time_resolution,    &
                                           ftype,    &
                                 m_month_of_year,    &
                                   m_day_of_year,    &
                                   m_hour_of_day,    &
                                   num_rec,          &
                                   ind_change      )
    implicit none

    integer m_month_of_year,           &
            m_day_of_year,             &
            m_hour_of_day

    integer ftype, time_resolution, ind_change
    integer num_rec, num_rec_old

    num_rec_old = num_rec
    ind_change=0

    if(ftype.eq.0) then
        num_rec=m_month_of_year
    else
        num_rec=(m_day_of_year-1)*24/time_resolution + m_hour_of_day/time_resolution +1
    end if

    if(num_rec/=num_rec_old) then
        ind_change=1
    endif

endsubroutine check_bc_field_change

!===============================================================================
subroutine bc_data_read(filename,      &
                           field,      &
                              nx,      &
                              ny,      &
                            m1,m2,     &
                            n1,n2,     &
                             nrec  )
    use rec_length
    use mpi_parallel_tools
    implicit none

    ! external parameters
    character*(*) filename              !forcing data file

    integer nx, ny, m1, m2, n1, n2, nrec      !its dimensions
    real(4) field(nx,ny) !data array
    integer m, n, ierr

    if (rank .eq. 0) then
        open (17,file=filename,status='old',access='direct', form='unformatted',     &
                      recl=(m2-m1+1)*(n2-n1+1)*lrecl,err=117)

        if (rank .eq. 0) then
           write(*,*) 'read boundary data: ', filename(1:len_trim (filename))
           write(*,'(2(a,i4),a,i6)') 'dimension is ',m2-m1+1,' x',n2-n1+1, ', number of record is ',nrec
        endif

        read(17,rec=nrec,err=119) ((field(m,n),m=m1,m2),n=n1,n2)

        close(17)
    endif

    call mpi_bcast(field, nx*ny, mpi_real4, 0, cart_comm, ierr)

    return

117    write(*,'(1x,a)')'   error in opening file:'
       write(*,'(5x,a)') filename(1:len_trim(filename))
       stop
119    write(*,'(1x,a)')'   error in reading file:'
       write(*,'(5x,a)') filename(1:len_trim(filename))
       stop

endsubroutine bc_data_read
