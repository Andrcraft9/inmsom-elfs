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

!  forming full file name for zonal wind stress
       call bc_full_name(ss_atmfiles(1),          & ! original file name
                         ss_atmfiles_fname(1),    & ! full file name
                              path2atmssdata,     & ! path to flux data
                              ind_clim_stress,    & ! index of climaticity (0-climatic data, 1 - real year data)
                                       m_year)      ! current year number
!  forming full file name for meridional wind stress
       call bc_full_name(ss_atmfiles(2),          & ! original file name
                         ss_atmfiles_fname(2),    & ! full file name
                              path2atmssdata,     & ! path to flux data
                              ind_clim_stress,    & ! index of climaticity (0-climatic data, 1 - real year data)
                                       m_year)      ! current year number
!  forming full file name for heat balance
       call bc_full_name(ss_atmfiles(3),          & ! original file name
                         ss_atmfiles_fname(3),    & ! full file name
                              path2atmssdata,     & ! path to flux data
                                 ind_clim_heat,   & ! index of climaticity (0-climatic data, 1 - real year data)
                                       m_year)      ! current year number
!  forming full file name for shortwave rad balance
       call bc_full_name(ss_atmfiles(4),          & ! original file name
                         ss_atmfiles_fname(4),    & ! full file name
                              path2atmssdata,     & ! path to flux data
                                 ind_clim_heat,   & ! index of climaticity (0-climatic data, 1 - real year data)
                                       m_year)      ! current year number
!  forming full file name for freshwater balance
       call bc_full_name(ss_atmfiles(5),          & ! original file name
                         ss_atmfiles_fname(5),    & ! full file name
                              path2atmssdata,     & ! path to flux data
                                 ind_clim_water,  & ! index of climaticity (0-climatic data, 1 - real year data)
                                       m_year)      ! current year number
!  forming full file name for air temperature
       call bc_full_name(ss_atmfiles(6),          & ! original file name
                         ss_atmfiles_fname(6),    & ! full file name
                              path2atmssdata,     & ! path to flux data
                                 ind_clim_tatm,   & ! index of climaticity (0-climatic data, 1 - real year data)
                                       m_year)      ! current year number
!  forming full file name for air humidity
       call bc_full_name(ss_atmfiles(7),          & ! original file name
                         ss_atmfiles_fname(7),    & ! full file name
                              path2atmssdata,     & ! path to flux data
                                 ind_clim_qatm,   & ! index of climaticity (0-climatic data, 1 - real year data)
                                       m_year)      ! current year number
!  forming full file name for wind zonal speed
       call bc_full_name(ss_atmfiles(8),          & ! original file name
                         ss_atmfiles_fname(8),    & ! full file name
                              path2atmssdata,     & ! path to flux data
                                 ind_clim_wind,   & ! index of climaticity (0-climatic data, 1 - real year data)
                                       m_year)      ! current year number
!  forming full file name for wind meridional speed
       call bc_full_name(ss_atmfiles(9),          & ! original file name
                         ss_atmfiles_fname(9),    & ! full file name
                              path2atmssdata,     & ! path to flux data
                                 ind_clim_wind,   & ! index of climaticity (0-climatic data, 1 - real year data)
                                       m_year)      ! current year number
!  forming full file name for SLP
       call bc_full_name(ss_atmfiles(10),         & ! original file name
                         ss_atmfiles_fname(10),   & ! full file name
                              path2atmssdata,     & ! path to flux data
                                 ind_clim_slpr,   & ! index of climaticity (0-climatic data, 1 - real year data)
                                       m_year)      ! current year number
!  forming full file name for downwelling longwave radiation
       call bc_full_name(ss_atmfiles(11),         & ! original file name
                         ss_atmfiles_fname(11),   & ! full file name
                              path2atmssdata,     & ! path to flux data
                                 ind_clim_rad,    & ! index of climaticity (0-climatic data, 1 - real year data)
                                       m_year)      ! current year number
!  forming full file name for  downwelling shortwave radiation
       call bc_full_name(ss_atmfiles(12),         & ! original file name
                         ss_atmfiles_fname(12),   & ! full file name
                              path2atmssdata,     & ! path to flux data
                                 ind_clim_rad,    & ! index of climaticity (0-climatic data, 1 - real year data)
                                       m_year)      ! current year number
!  forming full file name for liquid precipitation (rain)
       call bc_full_name(ss_atmfiles(13),         & ! original file name
                         ss_atmfiles_fname(13),   & ! full file name
                              path2atmssdata,     & ! path to flux data
                                 ind_clim_prec,   & ! index of climaticity (0-climatic data, 1 - real year data)
                                       m_year)      ! current year number
!  forming full file name for solid precipitation (snow)
       call bc_full_name(ss_atmfiles(14),         & ! original file name
                         ss_atmfiles_fname(14),   & ! full file name
                              path2atmssdata,     & ! path to flux data
                                ind_clim_prec,    & ! index of climaticity (0-climatic data, 1 - real year data)
                                      m_year)      ! current year number
   endif


!indicating of atm array change and reading data from files

      if(iabs(ksw_ssbc)<=2) then

! wind stress
       call check_bc_field_change(time_resolution_stress,    &
                                             ftype_stress,    &
                                          m_month_of_year,    &
                                            m_day_of_year,    &
                                            m_hour_of_day,    &
                                           num_rec_stress,    &
                                         ind_change_stress    )
       if(ind_change_stress>0) then

           call bc_data_read(ss_atmfiles_fname(1),      &
                                       a_stress_x,      &
                                              nxa,      &
                                              nya,      &
                                            1,nxa,      &
                                            1,nya,      &
                                   num_rec_stress  )
           call bc_data_read(ss_atmfiles_fname(2),      &
                                       a_stress_y,      &
                                              nxa,      &
                                              nya,      &
                                            1,nxa,      &
                                            1,nya,      &
                                   num_rec_stress  )

         ! converting wind stress to model units (m**2/s**2)
        a_stress_x=a_stress_x*confac_stress
        a_stress_y=a_stress_y*confac_stress

       endif

	end if
!------------------------------------------------------------------------
      if(iabs(ksw_ssbc)==2) then
! heat and SW-rad balance
       call check_bc_field_change(time_resolution_heat,    &
                                            ftype_heat,    &
                                       m_month_of_year,    &
                                         m_day_of_year,    &
                                         m_hour_of_day,    &
                                          num_rec_heat,    &
                                       ind_change_heat    )
       if(ind_change_heat>0) then
           call bc_data_read(ss_atmfiles_fname(3),      &
                                          a_hflux,      &
                                              nxa,      &
                                              nya,      &
                                            1,nxa,      &
                                            1,nya,      &
                                      num_rec_heat  )
           call bc_data_read(ss_atmfiles_fname(4),      &
                                          a_swrad,      &
                                              nxa,      &
                                              nya,      &
                                            1,nxa,      &
                                            1,nya,      &
                                      num_rec_heat  )
        ! converting heat fluxes to model units (W/m^2)
          a_hflux=a_hflux*confac_heat
          a_swrad=a_swrad*confac_heat

       endif
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

!------------------------------------------------------------------------
      if(iabs(ksw_ssbc)>=3) then

! Air temperature
       call check_bc_field_change(time_resolution_tatm,    &
                                            ftype_tatm,    &
                                       m_month_of_year,    &
                                         m_day_of_year,    &
                                         m_hour_of_day,    &
                                          num_rec_tatm,    &
                                       ind_change_tatm    )
       if(ind_change_tatm>0) then
           call bc_data_read(ss_atmfiles_fname(6),      &
                                           a_tatm,      &
                                              nxa,      &
                                              nya,      &
                                            1,nxa,      &
                                            1,nya,      &
                                     num_rec_tatm  )
       ! converting air temperature to model units (�C)
        a_tatm=a_tatm*confac_tatm

       endif

! Air humidity
       call check_bc_field_change(time_resolution_qatm,    &
                                            ftype_qatm,    &
                                       m_month_of_year,    &
                                         m_day_of_year,    &
                                         m_hour_of_day,    &
                                          num_rec_qatm,    &
                                       ind_change_qatm    )
       if(ind_change_qatm>0) then
           call bc_data_read(ss_atmfiles_fname(7),      &
                                           a_qatm,      &
                                              nxa,      &
                                              nya,      &
                                            1,nxa,      &
                                            1,nya,      &
                                     num_rec_qatm  )

       ! converting air specific humidity to model units (kg/kg)
        a_qatm=a_qatm*confac_qatm
       endif

! Wind speed
       call check_bc_field_change(time_resolution_wind,    &
                                            ftype_wind,    &
                                       m_month_of_year,    &
                                         m_day_of_year,    &
                                         m_hour_of_day,    &
                                          num_rec_wind,    &
                                       ind_change_wind    )
       if(ind_change_wind>0) then
           call bc_data_read(ss_atmfiles_fname(8),      &
                                           a_uwnd,      &
                                              nxa,      &
                                              nya,      &
                                            1,nxa,      &
                                            1,nya,      &
                                     num_rec_wind     )
           call bc_data_read(ss_atmfiles_fname(9),      &
                                           a_vwnd,      &
                                              nxa,      &
                                              nya,      &
                                            1,nxa,      &
                                            1,nya,      &
                                     num_rec_wind     )
         ! converting wind speed to model units (m/s)
         a_uwnd=a_uwnd*confac_wind
         a_vwnd=a_vwnd*confac_wind
       endif

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

! Downwelling radiation
       call check_bc_field_change(time_resolution_rad,    &
                                            ftype_rad,    &
                                      m_month_of_year,    &
                                        m_day_of_year,    &
                                        m_hour_of_day,    &
                                          num_rec_rad,    &
                                       ind_change_rad    )
       if(ind_change_rad>0) then
          call bc_data_read(ss_atmfiles_fname(11),      &
                                            a_lwr,      &
                                              nxa,      &
                                              nya,      &
                                            1,nxa,      &
                                            1,nya,      &
                                       num_rec_rad  )
          call bc_data_read(ss_atmfiles_fname(12),      &
                                            a_swr,      &
                                              nxa,      &
                                              nya,      &
                                            1,nxa,      &
                                            1,nya,      &
                                       num_rec_rad  )
       ! converting downwelling radiation to model units (W/m^2)
        a_lwr=a_lwr*confac_rad
        a_swr=a_swr*confac_rad
       endif

! Precipitation
       call check_bc_field_change(time_resolution_prec,    &
                                            ftype_prec,    &
                                       m_month_of_year,    &
                                         m_day_of_year,    &
                                         m_hour_of_day,    &
                                          num_rec_prec,    &
                                       ind_change_prec    )
       if(ind_change_prec>0) then
          call bc_data_read(ss_atmfiles_fname(13),      &
                                           a_rain,      &
                                              nxa,      &
                                              nya,      &
                                            1,nxa,      &
                                            1,nya,      &
                                      num_rec_prec  )
       ! converting rain precipitation to model units (kg/m**2/s)
        a_rain=a_rain*confac_prec

         if(prec_split>0) then
          call bc_data_read(ss_atmfiles_fname(14),      &
                                           a_snow,      &
                                              nxa,      &
                                              nya,      &
                                            1,nxa,      &
                                            1,nya,      &
                                      num_rec_prec  )
         ! converting rain precipitation to model units (kg/m**2/s)
          a_snow=a_snow*confac_prec
         endif

       endif

	end if

endsubroutine atm_data_time_interpol
!======================================================================
! oceanic data time interpolation
subroutine oc_data_time_interpol
use main_basin_pars
use mpi_parallel_tools
use ocean_bc
use time_integration
use key_switches
implicit none

! Forming full file names once a year

   if(m_time_changed(6)>0) then ! if new year comes

!  forming full file name for SST
       call bc_full_name(ss_ocfiles(1),          & ! original file name
                         ss_ocfiles_fname(1),    & ! full file name
                              path2ocssdata,     & ! path to flux data
                              ind_clim_sst,      & ! index of climaticity (0-climatic data, 1 - real year data)
                                       m_year)     ! current year number
!  forming full file name for SSS
       call bc_full_name(ss_ocfiles(2),          & ! original file name
                         ss_ocfiles_fname(2),    & ! full file name
                              path2ocssdata,     & ! path to flux data
                              ind_clim_sss,      & ! index of climaticity (0-climatic data, 1 - real year data)
                                       m_year)     ! current year number
!  forming full file name for runoff
       call bc_full_name(ss_ocfiles(3),          & ! original file name
                         ss_ocfiles_fname(3),    & ! full file name
                              path2ocssdata,     & ! path to flux data
                              ind_clim_runoff,   & ! index of climaticity (0-climatic data, 1 - real year data)
                                       m_year)     ! current year number
!  forming full file name for TLBC
       call bc_full_name(ss_ocfiles(4),          & ! original file name
                         ss_ocfiles_fname(4),    & ! full file name
                              path2ocssdata,     & ! path to flux data
                              ind_clim_tlbc,     & ! index of climaticity (0-climatic data, 1 - real year data)
                                       m_year)     ! current year number
!  forming full file name for SLBC
       call bc_full_name(ss_ocfiles(5),          & ! original file name
                         ss_ocfiles_fname(5),    & ! full file name
                              path2ocssdata,     & ! path to flux data
                              ind_clim_slbc,     & ! index of climaticity (0-climatic data, 1 - real year data)
                                       m_year)     ! current year number
!  forming full file name for ULBC
       call bc_full_name(ss_ocfiles(6),          & ! original file name
                         ss_ocfiles_fname(6),    & ! full file name
                              path2ocssdata,     & ! path to flux data
                              ind_clim_ulbc,     & ! index of climaticity (0-climatic data, 1 - real year data)
                                       m_year)     ! current year number
!  forming full file name for VLBC
       call bc_full_name(ss_ocfiles(7),          & ! original file name
                         ss_ocfiles_fname(7),    & ! full file name
                              path2ocssdata,     & ! path to flux data
                              ind_clim_ulbc,     & ! index of climaticity (0-climatic data, 1 - real year data)
                                       m_year)     ! current year number
!  forming full file name for HLBC
       call bc_full_name(ss_ocfiles(8),          & ! original file name
                         ss_ocfiles_fname(8),    & ! full file name
                              path2ocssdata,     & ! path to flux data
                              ind_clim_hlbc,     & ! index of climaticity (0-climatic data, 1 - real year data)
                                       m_year)     ! current year number

   endif


!indicating of atm array change and reading data from files

! SST
       call check_bc_field_change(time_resolution_sst,    &
                                            ftype_sst,    &
                                      m_month_of_year,    &
                                        m_day_of_year,    &
                                        m_hour_of_day,    &
                                          num_rec_sst,    &
                                       ind_change_sst    )
       if(ind_change_sst>0) then
           call bc_data_read(ss_ocfiles_fname(1),      &
                                         sst_obs,      &
                                              nx,      &
                                              ny,      &
                                          mmm,mm,      &
                                          nnn,nn,      &
                                      num_rec_sst  )
           sst_obs=sst_obs*confac_sst
       endif

! SSS
       call check_bc_field_change(time_resolution_sss,    &
                                            ftype_sss,    &
                                      m_month_of_year,    &
                                        m_day_of_year,    &
                                        m_hour_of_day,    &
                                          num_rec_sss,    &
                                       ind_change_sss    )
       if(ind_change_sss>0) then
           call bc_data_read(ss_ocfiles_fname(2),      &
                                         sss_obs,      &
                                              nx,      &
                                              ny,      &
                                          mmm,mm,      &
                                          nnn,nn,      &
                                      num_rec_sss  )
           sss_obs=sss_obs*confac_sss
       endif
! RUNOFF
     if(ksw_ssbc>=3) then
       call check_bc_field_change(time_resolution_runoff,    &
                                            ftype_runoff,    &
                                         m_month_of_year,    &
                                           m_day_of_year,    &
                                           m_hour_of_day,    &
                                          num_rec_runoff,    &
                                       ind_change_runoff    )
       if(ind_change_runoff>0) then
           call bc_data_read(ss_ocfiles_fname(3),      &
                                          runoff,      &
                                              nx,      &
                                              ny,      &
                                          mmm,mm,      &
                                          nnn,nn,      &
                                      num_rec_runoff  )
           runoff=runoff*confac_runoff
       endif
      endif

! T and S at open boundary
     if(ksw_lbc_ts>0) then
! TLBC
       call check_bc_field_change(time_resolution_tlbc,    &
                                            ftype_tlbc,    &
                                       m_month_of_year,    &
                                         m_day_of_year,    &
                                         m_hour_of_day,    &
                                          num_rec_tlbc,    &
                                       ind_change_tlbc    )
       if(ind_change_tlbc>0) then
           call bc_data_read(ss_ocfiles_fname(4),      &
                                           tlqbw,      &
                                     numb_of_lqp,      &
                                              nz,      &
                                   1,numb_of_lqp,      &
                                             1,nz,     &
                                      num_rec_tlbc  )
           tlqbw=tlqbw*confac_tlbc
       endif
! SLBC
       call check_bc_field_change(time_resolution_slbc,    &
                                            ftype_slbc,    &
                                       m_month_of_year,    &
                                         m_day_of_year,    &
                                         m_hour_of_day,    &
                                          num_rec_slbc,    &
                                       ind_change_slbc    )
       if(ind_change_slbc>0) then
           call bc_data_read(ss_ocfiles_fname(5),      &
                                           slqbw,      &
                                     numb_of_lqp,      &
                                              nz,      &
                                   1,numb_of_lqp,      &
                                             1,nz,     &
                                      num_rec_slbc  )
           slqbw=slqbw*confac_slbc
       endif

     endif

! U and V at open boundary
     if(ksw_lbc_uv>0) then
! ULBC
       call check_bc_field_change(time_resolution_ulbc,    &
                                            ftype_ulbc,    &
                                       m_month_of_year,    &
                                         m_day_of_year,    &
                                         m_hour_of_day,    &
                                          num_rec_ulbc,    &
                                       ind_change_ulbc    )
       if(ind_change_ulbc>0) then
           call bc_data_read(ss_ocfiles_fname(6),      &
                                           ulqbw,      &
                                     numb_of_lqp,      &
                                              nz,      &
                                   1,numb_of_lqp,      &
                                             1,nz,     &
                                      num_rec_ulbc  )
           call bc_data_read(ss_ocfiles_fname(7),      &
                                           vlqbw,      &
                                     numb_of_lqp,      &
                                              nz,      &
                                   1,numb_of_lqp,      &
                                             1,nz,     &
                                      num_rec_ulbc  )
           ulqbw=ulqbw*confac_ulbc
           vlqbw=vlqbw*confac_ulbc
       endif


     endif

! SSH at open boundary
     if(ksw_lbc_ssh>0) then
! HLBC
       call check_bc_field_change(time_resolution_hlbc,    &
                                            ftype_hlbc,    &
                                       m_month_of_year,    &
                                         m_day_of_year,    &
                                         m_hour_of_day,    &
                                          num_rec_hlbc,    &
                                       ind_change_hlbc    )
       if(ind_change_hlbc>0) then
           call bc_data_read(ss_ocfiles_fname(8),      &
                                         sshlqbw,      &
                                     numb_of_lqp,      &
                                               1,      &
                                   1,numb_of_lqp,      &
                                              1,1,     &
                                      num_rec_hlbc  )
           sshlqbw=sshlqbw*confac_hlbc
       endif


     endif

endsubroutine oc_data_time_interpol

!============================================================================
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
!=================================================================================
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
    integer m,n

    open (17,file=filename,status='old',access='direct', form='unformatted',     &
                  recl=(m2-m1+1)*(n2-n1+1)*lrecl,err=117)

    if (rank .eq. 0) then
       write(*,*) 'read boundary data: ', filename(1:len_trim (filename))
       write(*,'(2(a,i4),a,i6)') 'dimension is ',m2-m1+1,' x',n2-n1+1, ', number of record is ',nrec
    endif

    read(17,rec=nrec,err=119) ((field(m,n),m=m1,m2),n=n1,n2)

    close(17)

	return
    
117    write(*,'(1x,a)')'   error in opening file:'
      write(*,'(5x,a)') filename(1:len_trim(filename))
      stop
119    write(*,'(1x,a)')'   error in reading file:'
      write(*,'(5x,a)') filename(1:len_trim(filename))
      stop

endsubroutine bc_data_read
