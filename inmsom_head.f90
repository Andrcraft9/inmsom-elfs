program INMSOM
use time_integration
use key_switches
use mpi_parallel_tools
implicit none

character(128) fname
integer m, ierr
real*8 t_global, t_local

!call non_mpi_array_boundary_definition
call mpi_array_boundary_definition

m_sec_of_min   = 0     !second counter in minute
m_min_of_hour  = 0     !minute counter in hour
m_hour_of_day  = 0     !hour counter in day
m_day_of_month = 0     !day counter in month
m_day_of_year  = 0     !day counter in year
m_day_of_4yr   = 0     !day counter in 4-years
m_month_of_year= 0     !mon counter in year
m_month_of_4yr = 0     !mon counter in 4-years
m_year_of_4yr  = 0     !year counter in 4yrs
m_day          = 0     !model elapsed day counter starting from zero
m_month        = 0     !model elapsed month counter starting from zero
m_year         = 0     !year counter
m_4yr          = 0     !counter of 4-yr groups
m_time_changed = 0     !change indicator of time
key_time_print = 0     !key of printing time:0-not,1-print

seconds_of_day = 0.0d0 !Current seconds of the day

if(yr_type==0) then  !day distribution without leap-year
      ndays_in_4yr=ndays_noleap
elseif(yr_type==1) then   !day distribution with leap-year
      ndays_in_4yr=ndays_leap
endif

blank=' '
filepar='ocean_run.par'

!reading parameters from file
if (rank .eq. 0) then
    call readpar(filepar,comments,nofcom)
endif
call mpi_bcast(comments, 256*256, mpi_character, 0, cart_comm, ierr)

read(comments( 1),*) start_type          !Type of starting run (0 - from TS only, 1 - frim the full checkpoint)
read(comments( 2),*) time_step           !Model time step (in seconds)
read(comments( 3),*) run_duration        !Duration of the run in days
read(comments( 4),*) num_step            !Number of time step during the run
read(comments( 5),*) init_year           !Initial year number for the run
read(comments( 6),*) loc_data_wr_period  !Period for writing local instantaneous data (minutes)
read(comments( 7),*) glob_data_wr_period !Period for writing global time average data (minutes)
read(comments( 8),*) nstep_icedyn        !Number of internal time steps for ice dynamics
read(comments( 9),*) nstep_barotrop      !Number of internal time steps for barotropic task
call get_first_lexeme(comments(10), path2ocp    )  !path to checkpoints(results)
! Files with data on oceanic grid:
call get_first_lexeme(comments(11), path2ocssdata)   !path to ocean SS data
call get_first_lexeme(comments(12), ss_ocfiles(1)  )  !file with SST
call get_first_lexeme(comments(13), ss_ocfiles(2)  )  !file with SSS
call get_first_lexeme(comments(14), ss_ocfiles(3)  )  !file with river runoff
call get_first_lexeme(comments(15), ss_ocfiles(4)  )  !file with TLBC
call get_first_lexeme(comments(16), ss_ocfiles(5)  )  !file with SLBC
call get_first_lexeme(comments(17), ss_ocfiles(6)  )  !file with ULBC
call get_first_lexeme(comments(18), ss_ocfiles(7)  )  !file with VLBC
call get_first_lexeme(comments(19), ss_ocfiles(8)  )  !file with SSHLBC
! Files with data on atmospheric grid:
call get_first_lexeme(comments(20), path2atmssdata )  !path to atmospheric data
call get_first_lexeme(comments(21), ss_atmfiles(1) )  !file with      zonal wind stress (1 and 2 condition)
call get_first_lexeme(comments(22), ss_atmfiles(2) )  !file with meridional wind stress (1 and 2 condition)
call get_first_lexeme(comments(23), ss_atmfiles(3) )  !file with          heat balance (2 condition)
call get_first_lexeme(comments(24), ss_atmfiles(4) )  !file with shortwave rad balance (2 condition)
call get_first_lexeme(comments(25), ss_atmfiles(5) )  !file with    freshwater balance (2 condition)
call get_first_lexeme(comments(26), ss_atmfiles(6) )  !file with       air temperature (3 condition)
call get_first_lexeme(comments(27), ss_atmfiles(7) )  !file with          air humidity (3 condition)
call get_first_lexeme(comments(28), ss_atmfiles(8) )  !file with wind      zonal speed (3 condition)
call get_first_lexeme(comments(29), ss_atmfiles(9) )  !file with wind meridional speed (3 condition)
call get_first_lexeme(comments(30), ss_atmfiles(10) )  !file with SLP (3 condition)
call get_first_lexeme(comments(31), ss_atmfiles(11) )  !file with downwelling longwave radiation (3 condition)
call get_first_lexeme(comments(32), ss_atmfiles(12) )  !file with downwelling shortwave radiation (3 condition)
call get_first_lexeme(comments(33), ss_atmfiles(13) )  !file with wind liquid precipitation (rain)
call get_first_lexeme(comments(34), ss_atmfiles(14) )  !file with wind solid precipitation (snow)
call get_first_lexeme(comments(35), atmask)

if(loc_data_wr_period>0.0001) then
 key_write_local=1
else
 key_write_local=0
endif

if(glob_data_wr_period>0.0001) then
 key_write_global=1
else
 key_write_global=0
endif

time_write_global=0

!computing some time parameters
time_step_m = time_step/60.00d0
time_step_h = time_step/3600.00d0
time_step_d = time_step/86400.00d0
nstep_per_day=nint(86400.0d0/time_step)       !NUMBER OF STEP PER DAY

!for local output
loc_data_wr_period = max(min(loc_data_wr_period,1440.0),sngl(time_step_m)) !The maximum local output period is 1 day, minimum is time step
loc_data_wr_period_step = nint(loc_data_wr_period/time_step_m)   !period in steps to write to write local data
loc_data_tstep = loc_data_wr_period * 60.0                      !Time step in seconds for writing local data

year_loc=init_year
 mon_loc=1
 day_loc=int(loc_data_wr_period/1440.0)+1
 hour_loc=mod(int(loc_data_wr_period/60.0),24)
  min_loc=mod(int(loc_data_wr_period),60)

!for global output
if(abs(glob_data_wr_period)>1440.5) then
 monthly_output=1
 glob_data_tstep = 86400.0*30.0

year_glob=init_year
 mon_glob=1
 day_glob=15
hour_glob=0
 min_glob=0
else
 monthly_output=0
 glob_data_wr_period = max(min(abs(glob_data_wr_period),1440.0),sngl(time_step_m))
 glob_data_wr_period_step = nint(glob_data_wr_period/time_step_m)   !period in steps to write to write global data
 glob_data_tstep = glob_data_wr_period * 60.0                        !Time step in seconds for writing global data

 year_glob=init_year
  mon_glob=1
  day_glob=1
 hour_glob=mod(int(glob_data_wr_period/2.0/60.0),24)
  min_glob=mod(int(glob_data_wr_period/2.0),60)
endif

!Allocating main arrays
call model_grid_allocate
call ocean_variables_allocate

!Initializing ocean model parameters
call ocean_model_parameters(time_step)
if (rank .eq. 0) print *, "--------------------END OF OCEAN MODEL PARAMETERS----------------------"

!Initializing SW init conditions

call sw_only_inicond(1, path2ocp)
!call sw_test2

!call  parallel_local_output(path2ocp,  &
!                          1,  &
!                   year_loc,  &
!                    mon_loc,  &
!                    day_loc,  &
!                   hour_loc,  &
!                    min_loc,  &
!             loc_data_tstep,  &
!                    yr_type  )
!call mpi_finalize(ierr)
!stop

!-----------------------------------------------------------------------
if (rank .eq. 0) then
    write(*,'(2x,5hstep:,f7.2,4hhrs;,f9.2,4hsec;,f9.5,4hday.)')     &
             time_step_h,   time_step,    time_step_d
    write(*,'(2x,15hduration of run:,f9.2,6h days.)') run_duration
endif

key_time_print=1
call model_time_def(   num_step,            &     !step counter,            input
                       time_step,           &    !time step in seconds,    input
                       ndays_in_4yr,        &    !integer day distribution in 4-years (49 months)
                       seconds_of_day,      &    !current seconds in day  ,output
                       m_sec_of_min,        &    !second counter in minute,output
                       m_min_of_hour,       &    !minute counter in hour  ,output
                       m_hour_of_day,       &    !hour counter in day     ,output
                       m_day_of_month,      &    !day counter in month    ,output
                       m_day_of_year,       &    !day counter in year     ,output
                       m_day_of_4yr,        &    !day counter in 4-years  ,output
                       m_month_of_year,     &    !mon counter in year     ,output
                       m_month_of_4yr,      &    !mon counter in 4-years  ,output
                       m_year_of_4yr,       &    !year counter in 4yrs    ,output
                       m_day,               &    !model elapsed day counter starting from zero
                       m_month,             &    !model elapsed month counter starting from zero
                       m_year,              &    !year counter            ,output
                       m_4yr,               &    !counter of 4-yr groups  ,output
                       m_time_changed,      &    !change indicator of time,output
                       key_time_print,      &    !key of printing time:0-not,1-print
                       init_year)                !initial real-time year

num_step_max=int8(run_duration*nstep_per_day)

if (rank .eq. 0) then
    write(*,*)'=================================================================='
    write(*,*)'------------------Starting model time integration-----------------'
    write(*,*)'=================================================================='
endif

call init_times
call start_timer(t_global)
do while(num_step<num_step_max)

  call start_timer(t_local)
!computing one step of ocean dynamics
  call shallow_water_model_step(time_step, nstep_barotrop)
  call end_timer(t_local)
  time_model_step = time_model_step + t_local

!-----------------------------------------------------------
!moving to the next time step
 num_step=num_step+1
 key_time_print=0
 call model_time_def(   num_step,           &     !step counter,            input
                       time_step,           &    !time step in seconds,    input
                       ndays_in_4yr,        &    !integer day distribution in 4-years (49 months)
                       seconds_of_day,      &    !current seconds in day  ,output
                       m_sec_of_min,        &    !second counter in minute,output
                       m_min_of_hour,       &    !minute counter in hour  ,output
                       m_hour_of_day,       &    !hour counter in day     ,output
                       m_day_of_month,      &    !day counter in month    ,output
                       m_day_of_year,       &    !day counter in year     ,output
                       m_day_of_4yr,        &    !day counter in 4-years  ,output
                       m_month_of_year,     &    !mon counter in year     ,output
                       m_month_of_4yr,      &    !mon counter in 4-years  ,output
                       m_year_of_4yr,       &    !year counter in 4yrs    ,output
                       m_day,               &    !model elapsed day counter starting from zero
                       m_month,             &    !model elapsed month counter starting from zero
                       m_year,              &    !year counter            ,output
                       m_4yr,               &    !counter of 4-yr groups  ,output
                       m_time_changed,      &    !change indicator of time,output
                       key_time_print,      &    !key of printing time:0-not,1-print
                       init_year)                !initial real-time year

!Write local data

if( key_write_local>0) then

 if(mod(num_step,loc_data_wr_period_step)==0) then

  nrec_loc=num_step/loc_data_wr_period_step

  call start_timer(t_local)
  call  parallel_local_output(path2ocp,  &
                     nrec_loc,  &
                     year_loc,  &
                      mon_loc,  &
                      day_loc,  &
                     hour_loc,  &
                      min_loc,  &
               loc_data_tstep,  &
                      yr_type  )
  call end_timer(t_local)
  time_output = time_output + t_local

  call model_time_print(num_step,         &
                        m_sec_of_min,     &    !second counter in minute,output
                        m_min_of_hour,    &    !minute counter in hour  ,output
                        m_hour_of_day,    &    !hour counter in day     ,output
                        m_day_of_month,   &    !day counter in month    ,output
                        m_day_of_year,    &    !day counter in year     ,output
                        m_day_of_4yr,     &    !day counter in 4-years  ,output
                        m_month_of_year,  &    !mon counter in year     ,output
                        m_month,          &    !model elapsed month counter starting from zero
                        m_year )               !year counter            ,output
 endif

endif

enddo
!End of time cycle
call end_timer(t_global)
if (rank .eq. 0) print *, "Global time: ", t_global
if (rank .eq. 0) print *, "Steps: ", num_step
call print_times

!deallocating the arrays
call ocean_variables_deallocate
call model_grid_deallocate

call mpi_finalize(ierr)

endprogram INMSOM
