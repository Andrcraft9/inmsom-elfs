!======================================================================
subroutine build_intrp_mtrx(path2flux,atmask)
use main_basin_pars
use mpi_parallel_tools
use basin_grid
use atm_pars
use atm2oc_interpol
use rec_length
implicit none

character(256) filename
character comment*80
character*(*) path2flux , atmask

 integer,allocatable:: atm_mask(:,:)  !sea-land mask for atm grid

 character frmt*16
 integer i,j,m,n,ierr

 allocate(atm_mask(nxa,nya))

 ! atmospheric grid initialize
 ! case of regular grid on longitude
 if(xagr_type==0) then
      do i=1,nxa
	  xa(i)=x0a+dfloat(i-1)*dxa
	end do
 else !  case of irregular grid on longitude
      do i=1,nxa
	  xa(i)=xa_levels(i)
	end do
 endif

 ! case of regular grid on latitude
 if(yagr_type==0) then
      do j=1,nya
	  ya(j)=y0a+dfloat(j-1)*dya
	end do
 else !  case of irregular grid on latitude
      do j=1,nya
	  ya(j)=ya_levels(j)
	end do
 endif

 write(frmt,1000) nxa
1000  format('(',i9,'i1)')

!     initialization of atmospheric sea-land mask
      if(atmask.eq.'NONE'.or.atmask.eq.'none') then
            atm_mask=0
      else
	      atm_mask=1
! full file name of atmospheric mask on atmospheric grid
        call fulfname(filename,path2flux,atmask,ierr)

        open (11,file=filename,status='old',recl=nxa*lrecl)
        if (rank .eq. 0) then
           write(*,'(a,a)') ' read file with atmospheric mask: ', filename(1:len_trim(filename))
        endif
        do n=nya,1,-1
            read(11,frmt,end=99) (atm_mask(m,n),m=1,nxa)
        enddo
        close(11)

      end if

	   if (rank .eq. 0) write(*,*) 'building matrix for interpolation from atm to ocean'
!  calculating interpolation matrix elements

 call weight_matrix_intrp_next(xa,       & !array(1-D) of input x-grid values (input)
                               ya,       & !array(1-D) of input y-grid values (input)
                              nxa,       & !number of input x-grid points(input)
                              nya,       & !number of input y-grid points(input)
                        geo_lon_t,       & !array(2D) of output x-grid values in input coordinate system (input)
                        geo_lat_t,       & !array(2D) of output y-grid values in input coordinate system (input)
                        bnd_x1,bnd_x2,   & !boundaries of arrays for output x-grid coordinates (input)
                        bnd_y1,bnd_y2,   & !boundaries of arrays for output y-grid coordinates (input)
                        atm_mask,        & !input sea-land mask (input)
                             lu1,        & !output sea-land mask (input)
                     i_input_a2o,        & !x-grid numbers of input grid for output grid (output)
                     j_input_a2o,        & !y-grid numbers of input grid for output grid (output)
                   wght_mtrx_a2o,        & !nonzero matrix elements of interpolation (output)
                           indper,       & !index of input grid periodicity:=0 -nonperiodic,=1 -periodic case
                                1,       & !filling missed values (0-no, 1-yes)
                                1,       & !first significant point in x-direction for input grid (input)
                              nxa,       & ! last significant point in x-direction for input grid  (input)
                                1,       & !first significant point in y-direction for input grid  (input)
                              nya,       & ! last significant point in y-direction for input grid  (input)
                         nx_start-1,     & !first significant point in x-direction for output grid  (input)
                         nx_end+1,       & ! last significant point in x-direction for output grid  (input)
                         ny_start-1,     & !first significant point in y-direction for output grid  (input)
                         ny_end+1    )     ! last significant point in y-direction for output grid  (input)

      deallocate(atm_mask)

	return

99    write(*,*)' error in reading file atmmask.arr ', filename(1:len_trim(filename))
      stop 1

endsubroutine build_intrp_mtrx

!======================================================================
! atmospheric data spatial interpolation fron atm to ocean grid
subroutine atm_data_spatial_interpol
use main_basin_pars
use mpi_parallel_tools
use basin_grid
use ocean_variables
use atm_pars
use atm_forcing
use atm2oc_interpol
!use seaice
implicit none

 integer m, n
!-------------------------------------------------------------------------------------------------
 if(iabs(ksw_ssbc)<=2) then

   if(ind_change_stress>0) then
        if (rank .eq. 0) then
            write(*,'(4(a,i4))') 'Spatial interpolation of wind stress:     from ',nxa,' x',nya,' to ',mm-mmm+1,' x',nn-nnn+1
        endif
       !spatial interpolation of wind stress
        call interpolrot_vec(nxa,         &  !number of x-grid points(input)
                             nya,         &  !number of y-grid points(input)
                          bnd_x1,bnd_x2,  &  !number of x-grid points(output)
                          bnd_y1,bnd_y2,  &  !number of y-grid points(output)
                       a_stress_x,        &  !input zonal data array
                       a_stress_y,        &  !input meridional data array
                             taux,        &  !output zonal data array
                             tauy,        &  !output meridional data array
                      i_input_a2o,        &  !x-grid numbers of input grid for output grid
                      j_input_a2o,        &  !y-grid numbers of input grid for output grid
                      wght_mtrx_a2o,      &  !nonzero matrix elements of interpolation
                        rotvec_coeff,     &  !angles between parallels
                                  lu,     &  !sea-land mask for output data
                               0.0d0,     &
                          nx_start-1,     &  !first significant point in x-direction (output)
                          nx_end+1,       &  ! last significant point in x-direction (output)
                          ny_start-1,     &  !first significant point in y-direction (output)
                          ny_end+1 )        ! last significant point in y-direction (output)

   endif

 endif
!-------------------------------------------------------------------------------------------
 if(iabs(ksw_ssbc)==2) then

   if(ind_change_heat>0) then
       if (rank .eq. 0) then
           write(*,'(4(a,i4))') 'Spatial interpolation of heat balance:    from ',nxa,' x',nya,' to ',mm-mmm+1,' x',nn-nnn+1
       endif
      !spatial interpolation of heat balance
       call interpolrot_scal(nxa,         &  !number of x-grid points(input)
                            nya,         &  !number of y-grid points(input)
                         bnd_x1,bnd_x2,  &  !number of x-grid points(output)
                         bnd_y1,bnd_y2,  &  !number of y-grid points(output)
                               a_hflux,  &  !input data array
                                hf_tot,  &  !output data array
                           i_input_a2o,  &  !x-grid numbers of input grid for output grid
                           j_input_a2o,  &  !y-grid numbers of input grid for output grid
                         wght_mtrx_a2o,  &  !nonzero matrix elements of interpolation
                                    lu,  &  !sea-land mask for output data
                                 0.0d0,  &
                            nx_start-1,  &  !first significant point in x-direction (output)
                              nx_end+1,  &  ! last significant point in x-direction (output)
                            ny_start-1,  &  !first significant point in y-direction (output)
                              ny_end+1 )    ! last significant point in y-direction (output)

       if (rank .eq. 0) then
           write(*,'(4(a,i4))') 'Spatial interpolation of SW-balance:      from ',nxa,' x',nya,' to ',mm-mmm+1,' x',nn-nnn+1
       endif
      !spatial interpolation of SW-rad balance
       call interpolrot_scal(nxa,         &  !number of x-grid points(input)
                            nya,         &  !number of y-grid points(input)
                         bnd_x1,bnd_x2,  &  !number of x-grid points(output)
                         bnd_y1,bnd_y2,  &  !number of y-grid points(output)
                               a_swrad,  &  !input data array
                                sw_bal,  &  !output data array
                           i_input_a2o,  &  !x-grid numbers of input grid for output grid
                           j_input_a2o,  &  !y-grid numbers of input grid for output grid
                         wght_mtrx_a2o,  &  !nonzero matrix elements of interpolation
                                    lu,  &  !sea-land mask for output data
                                 0.0d0,  &
                            nx_start-1,  &  !first significant point in x-direction (output)
                              nx_end+1,  &  ! last significant point in x-direction (output)
                            ny_start-1,  &  !first significant point in y-direction (output)
                              ny_end+1 )    ! last significant point in y-direction (output)
   endif


   if(ind_change_water>0) then
       if (rank .eq. 0) then
           write(*,'(4(a,i4))') 'Spatial interpolation of freshwater balance: from ',nxa,' x',nya,' to ',mm-mmm+1,' x',nn-nnn+1
       endif
      !spatial interpolation of fresh water balance
       call interpolrot_scal(nxa,         &  !number of x-grid points(input)
                            nya,         &  !number of y-grid points(input)
                         bnd_x1,bnd_x2,  &  !number of x-grid points(output)
                         bnd_y1,bnd_y2,  &  !number of y-grid points(output)
                               a_wflux,  &  !input data array
                                wf_tot,  &  !output data array
                           i_input_a2o,  &  !x-grid numbers of input grid for output grid
                           j_input_a2o,  &  !y-grid numbers of input grid for output grid
                         wght_mtrx_a2o,  &  !nonzero matrix elements of interpolation
                                    lu,  &  !sea-land mask for output data
                                 0.0d0,  &
                            nx_start-1,  &  !first significant point in x-direction (output)
                              nx_end+1,  &  ! last significant point in x-direction (output)
                            ny_start-1,  &  !first significant point in y-direction (output)
                              ny_end+1 )    ! last significant point in y-direction (output)
   endif

 endif

!----------------------------------------------------------------------------------------
 if(iabs(ksw_ssbc)>=3) then

  if(ind_change_slpr>0) then
     !spatial interpolation of SLP
     if (rank .eq. 0) then
         write(*,'(4(a,i4))') 'Spatial interpolation of slp:             from ',nxa,' x',nya,' to ',mm-mmm+1,' x',nn-nnn+1
     endif
     call interpolrot_scal(nxa,         &  !number of x-grid points(input)
                           nya,         &  !number of y-grid points(input)
                        bnd_x1,bnd_x2,  &  !number of x-grid points(output)
                        bnd_y1,bnd_y2,  &  !number of y-grid points(output)
                               a_slpr,  &  !input data array
                                 slpr,  &  !output data array
                          i_input_a2o,  &  !x-grid numbers of input grid for output grid
                          j_input_a2o,  &  !y-grid numbers of input grid for output grid
                        wght_mtrx_a2o,  &  !nonzero matrix elements of interpolation
                                   lu,  &  !sea-land mask for output data
                                0.0d0,  &
                           nx_start-1,  &  !first significant point in x-direction (output)
                             nx_end+1,  &  ! last significant point in x-direction (output)
                           ny_start-1,  &  !first significant point in y-direction (output)
                             ny_end+1 )    ! last significant point in y-direction (output)
  endif


  if(ind_change_rad>0) then
       if (rank .eq. 0) then
          write(*,'(4(a,i4))') 'Spatial interpolation of DW-LW-rad:       from ',nxa,' x',nya,' to ',mm-mmm+1,' x',nn-nnn+1
       endif
       !spatial interpolation of downwelling LW-radiation
       call interpolrot_scal(nxa,         &  !number of x-grid points(input)
                             nya,         &  !number of y-grid points(input)
                          bnd_x1,bnd_x2,  &  !number of x-grid points(output)
                          bnd_y1,bnd_y2,  &  !number of y-grid points(output)
                                  a_lwr,  &  !input data array
                                    lwr,  &  !output data array
                            i_input_a2o,  &  !x-grid numbers of input grid for output grid
                            j_input_a2o,  &  !y-grid numbers of input grid for output grid
                          wght_mtrx_a2o,  &  !nonzero matrix elements of interpolation
                                     lu,  &  !sea-land mask for output data
                                  0.0d0,  &
                             nx_start-1,  &  !first significant point in x-direction (output)
                               nx_end+1,  &  ! last significant point in x-direction (output)
                             ny_start-1,  &  !first significant point in y-direction (output)
                               ny_end+1 )    ! last significant point in y-direction (output)

       if (rank .eq. 0) then
           write(*,'(4(a,i4))') 'Spatial interpolation of DW-SW-rad:       from ',nxa,' x',nya,' to ',mm-mmm+1,' x',nn-nnn+1
       endif
       !spatial interpolation of downwelling SW-radiation
       call interpolrot_scal(nxa,         &  !number of x-grid points(input)
                             nya,         &  !number of y-grid points(input)
                          bnd_x1,bnd_x2,  &  !number of x-grid points(output)
                          bnd_y1,bnd_y2,  &  !number of y-grid points(output)
                                  a_swr,  &  !input data array
                                    swr,  &  !output data array
                            i_input_a2o,  &  !x-grid numbers of input grid for output grid
                            j_input_a2o,  &  !y-grid numbers of input grid for output grid
                          wght_mtrx_a2o,  &  !nonzero matrix elements of interpolation
                                     lu,  &  !sea-land mask for output data
                                  0.0d0,  &
                             nx_start-1,  &  !first significant point in x-direction (output)
                               nx_end+1,  &  ! last significant point in x-direction (output)
                             ny_start-1,  &  !first significant point in y-direction (output)
                               ny_end+1 )    ! last significant point in y-direction (output)
  endif

  if(ind_change_prec>0) then
        if (rank .eq. 0) then
          write(*,'(4(a,i4))') 'Spatial interpolation of rain:            from ',nxa,' x',nya,' to ',mm-mmm+1,' x',nn-nnn+1
        endif
       !spatial interpolation of rain precipitation
        call interpolrot_scal(nxa,         &  !number of x-grid points(input)
                              nya,         &  !number of y-grid points(input)
                           bnd_x1,bnd_x2,  &  !number of x-grid points(output)
                           bnd_y1,bnd_y2,  &  !number of y-grid points(output)
                                  a_rain,  &  !input data array
                                    rain,  &  !output data array
                             i_input_a2o,  &  !x-grid numbers of input grid for output grid
                             j_input_a2o,  &  !y-grid numbers of input grid for output grid
                           wght_mtrx_a2o,  &  !nonzero matrix elements of interpolation
                                      lu,  &  !sea-land mask for output data
                                   0.0d0,  &
                              nx_start-1,  &  !first significant point in x-direction (output)
                                nx_end+1,  &  ! last significant point in x-direction (output)
                              ny_start-1,  &  !first significant point in y-direction (output)
                                ny_end+1 )    ! last significant point in y-direction (output)

   if(prec_split>0) then
       if (rank .eq. 0) then
           write(*,'(4(a,i4))') 'Spatial interpolation of snow:            from ',nxa,' x',nya,' to ',mm-mmm+1,' x',nn-nnn+1
       endif
       !spatial interpolation of snow precipitation
       call interpolrot_scal(nxa,         &  !number of x-grid points(input)
                             nya,         &  !number of y-grid points(input)
                          bnd_x1,bnd_x2,  &  !number of x-grid points(output)
                          bnd_y1,bnd_y2,  &  !number of y-grid points(output)
                                 a_snow,  &  !input data array
                                   snow,  &  !output data array
                            i_input_a2o,  &  !x-grid numbers of input grid for output grid
                            j_input_a2o,  &  !y-grid numbers of input grid for output grid
                          wght_mtrx_a2o,  &  !nonzero matrix elements of interpolation
                                     lu,  &  !sea-land mask for output data
                                  0.0d0,  &
                             nx_start-1,  &  !first significant point in x-direction (output)
                               nx_end+1,  &  ! last significant point in x-direction (output)
                             ny_start-1,  &  !first significant point in y-direction (output)
                               ny_end+1 )    ! last significant point in y-direction (output)
    else
       snow=0.0d0
    endif

  endif

  if(ind_change_tatm>0) then
       if (rank .eq. 0) then
           write(*,'(4(a,i4))') 'Spatial interpolation of air temperature: from ',nxa,' x',nya,' to ',mm-mmm+1,' x',nn-nnn+1
       endif
       !spatial interpolation of air temperature
       call interpolrot_scal(nxa,         &  !number of x-grid points(input)
                             nya,         &  !number of y-grid points(input)
                          bnd_x1,bnd_x2,  &  !number of x-grid points(output)
                          bnd_y1,bnd_y2,  &  !number of y-grid points(output)
                                 a_tatm,  &  !input data array
                                   tatm,  &  !output data array
                            i_input_a2o,  &  !x-grid numbers of input grid for output grid
                            j_input_a2o,  &  !y-grid numbers of input grid for output grid
                          wght_mtrx_a2o,  &  !nonzero matrix elements of interpolation
                                     lu,  &  !sea-land mask for output data
                                  0.0d0,  &
                             nx_start-1,  &  !first significant point in x-direction (output)
                               nx_end+1,  &  ! last significant point in x-direction (output)
                             ny_start-1,  &  !first significant point in y-direction (output)
                               ny_end+1 )    ! last significant point in y-direction (output)
  endif

  if(ind_change_qatm>0) then
       if (rank .eq. 0) then
          write(*,'(4(a,i4))') 'Spatial interpolation of air humidity:    from ',nxa,' x',nya,' to ',mm-mmm+1,' x',nn-nnn+1
       endif
       !spatial interpolation of air humidity
       call interpolrot_scal(nxa,         &  !number of x-grid points(input)
                             nya,         &  !number of y-grid points(input)
                          bnd_x1,bnd_x2,  &  !number of x-grid points(output)
                          bnd_y1,bnd_y2,  &  !number of y-grid points(output)
                                 a_qatm,  &  !input data array
                                   qatm,  &  !output data array
                            i_input_a2o,  &  !x-grid numbers of input grid for output grid
                            j_input_a2o,  &  !y-grid numbers of input grid for output grid
                          wght_mtrx_a2o,  &  !nonzero matrix elements of interpolation
                                     lu,  &  !sea-land mask for output data
                                  0.0d0,  &
                             nx_start-1,  &  !first significant point in x-direction (output)
                               nx_end+1,  &  ! last significant point in x-direction (output)
                             ny_start-1,  &  !first significant point in y-direction (output)
                               ny_end+1 )    ! last significant point in y-direction (output)
  endif

  if(ind_change_wind>0) then
    if (rank .eq. 0) then
        write(*,'(4(a,i4))') 'Spatial interpolation of wind speed:      from ',nxa,' x',nya,' to ',mm-mmm+1,' x',nn-nnn+1
    endif
    !spatial interpolation of wind speed
    call interpolrot_vec(nxa,         &  !number of x-grid points(input)
                         nya,         &  !number of y-grid points(input)
                      bnd_x1,bnd_x2,  &  !number of x-grid points(output)
                      bnd_y1,bnd_y2,  &  !number of y-grid points(output)
                       a_uwnd,        &  !input zonal data array
                       a_vwnd,        &  !input meridional data array
                         uwnd,        &  !output zonal data array
                         vwnd,        &  !output meridional data array
                  i_input_a2o,        &  !x-grid numbers of input grid for output grid
                  j_input_a2o,        &  !y-grid numbers of input grid for output grid
                  wght_mtrx_a2o,      &  !nonzero matrix elements of interpolation
                    rotvec_coeff,     &  !angles between parallels
                              lu,     &  !sea-land mask for output data
                           0.0d0,     &
                    nx_start-1,       &  !first significant point in x-direction (output)
                    nx_end+1,         &  ! last significant point in x-direction (output)
                    ny_start-1,       &  !first significant point in y-direction (output)
                    ny_end+1 )           ! last significant point in y-direction (output)

!$omp parallel do private(m,n)
  do n=ny_start-1,ny_end+1
    do m=nx_start-1, nx_end+1
     if(lu(m,n)>0.5) then
      wind(m,n)=dsqrt(uwnd(m,n)*uwnd(m,n)+vwnd(m,n)*vwnd(m,n))
     endif
    enddo
  enddo
!$omp end parallel do

  endif

  if(prec_split==0.and.(ind_change_tatm>0.or.ind_change_prec>0)) then
   if (rank .eq. 0) write(*,*) 'Getting snow data from total precipitation'
!if rain and snow are mixed
!$omp parallel do private(m,n)
   do n=ny_start-1, ny_end+1
    do m=nx_start-1, nx_end+1
     if(lu(m,n)>0.5) then
      rain(m,n)=rain(m,n)+snow(m,n)
      snow(m,n)=0.0d0
      if(tatm(m,n)<0.0d0) then
       snow(m,n)=rain(m,n)
       rain(m,n)=0.0d0
      endif
     endif
    enddo
   enddo
!$omp end parallel do
    endif

 endif

endsubroutine atm_data_spatial_interpol
