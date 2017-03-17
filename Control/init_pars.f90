subroutine ocean_model_parameters(tau)
use main_basin_pars
use mpi_parallel_tools
use basin_grid
use ocean_variables

implicit none
character(256)    t_mask_file,       &  !name of file with temperature point sea-land mask
       bottom_topography_file,       &  !name of file with bottom topography
       help_string
real(4) array4(bnd_x1:bnd_x2,bnd_y1:bnd_y2)
integer  m, n, k, ierr
real(8) tau

! define parameters of task
! description of parameters see in file with mame filepar

open (90,file='phys_proc.par',status='old')
 read(90,*) ksw_ts            !Temperature and salinity equation solving (0 - no, 1 - yes)
 read(90,*) ksw_age           !Ideal age equation solving (0 - no, 1 - yes)
 read(90,*) ksw_pt            !Passive tracer equation solving (0 - no, 1 - yes)
 read(90,*) ksw_uv            !Momentum equation solving (0 - no, 1 - yes)
 read(90,*) ksw_lat           !Lateral 2nd order mix parametrization (0 - constant coeff, 1 - Smagorinski)
 read(90,*) ksw_lat4          !Lateral 4nd order momentum mix parametrization (0 - no, 1 - yes)
 read(90,*) ksw_vert          !vertical mix parametrization (0 - constant coeff, 1 - Pacanowski&Philander)
 read(90,*) ksw_dens          !pressure gradient computation (0 - no (constant density), 1 - yes)
 read(90,*) ksw_ice_th        !sea ice thermodynamics using (0 - no, 1 - yes)
 read(90,*) ksw_ice_tran      !sea ice transport using (0 - no, 1 - yes)
 read(90,*) ksw_ice_dyn       !sea ice dynamics using (0 - no, 1 - yes)
 read(90,*) ksw_ssbc          !Type of surface boundary conditions (1 - surface T&S and wind stress are set; 2 - T&S fluxes and wind stress are set; 3 - T&S fluxes and wind stress are computed
 read(90,*) ksw_wflux         !normalize global mean salt balance (0 - no, 1 - normalize water flux, 2 - normalize salinity flux)
 read(90,*) ksw_lbc_ts        !open boundary conditions for T&S using (0 - no, 1 - yes)
 read(90,*) ksw_lbc_uv        !open boundary conditions for U&V using (0 - no, 1 - yes)
 read(90,*) ksw_lbc_ssh       !open boundary conditions for SSH using (0 - no, 1 - yes)
 read(90,*) sst_relax         !Relaxation coefficient for temperature [m/s]
 read(90,*) sss_relax         !Relaxation coefficient for salinity [m/s]
 read(90,*) ldiff_ts          !lateral diffusion for temperature [m**2/s]
 read(90,*) lvisc_2           !lateral  vicosity(2nd order)[m**2/s]
 read(90,*) lvisc_4           !lateral  vicosity(4th order) [undimensional]
 read(90,*) tsfrac_lat        !fraction of salinity lateral diffusion due to one for temperature
 read(90,*) vdiff_ts_min      !vertical background diffusion coefficient for T [m**2/s]
 read(90,*) vdiff_ts_max      !vertical max(top) diffusion coefficient for T [m**2/s]
 read(90,*) vvisc_min         !vertical background viscous coefficient [m**2/s]
 read(90,*) vvisc_max         !vertical max(top) viscous coefficient [m**2/s]
 read(90,*) tsfrac_vert       !fraction of salinity vertical diffusion due to one for temperature
 read(90,*) z_frac            !weight coefficient for lateral Z-Diffusion
 read(90,*) r_frac            !weight coefficient for lateral R-Diffusion
 read(90,*) gm_ratio          !weight coefficient for Gent&McWilliams transport

 help_string =' '
 read (90,'(a)') help_string   ! file with t-mask'
 call get_first_lexeme(help_string ,t_mask_file   )

 help_string =' '
 read (90,'(a)') help_string  ! file with bottom topography'
 call get_first_lexeme(help_string ,bottom_topography_file  )

close(90)

 if (rank .eq. 0) then
     write(*,'(i7,a)') ksw_ts,   ' - Temperature and salinity equation solving'
     write(*,'(i7,a)') ksw_age,  ' - Ideal age equation solving'
     write(*,'(i7,a)') ksw_pt,   ' - Passive tracer equation solving'
     write(*,'(i7,a)') ksw_uv,   ' - Momentum equation solving'
     write(*,'(i7,a)') ksw_lat,      ' - Lateral 2nd order mix parametrization'
     write(*,'(i7,a)') ksw_lat4,     ' - Lateral 4nd order momentum mix parametrization'
     write(*,'(i7,a)') ksw_vert,     ' - Vertical mix parametrization'
     write(*,'(i7,a)') ksw_dens,     ' - Pressure gradient computation'
     write(*,'(i7,a)') ksw_ice_th,   ' - Sea ice thermodynamics using'
     write(*,'(i7,a)') ksw_ice_tran, ' - Sea ice transport using'
     write(*,'(i7,a)') ksw_ice_dyn,  ' - Sea ice dynamics using'
     write(*,'(i7,a)') ksw_ssbc,     ' - Type of surface boundary conditions'
     write(*,'(i7,a)') ksw_wflux,    ' - Normalize global mean salt balance'
     write(*,'(i7,a)') ksw_lbc_ts,   ' - Open boundary conditions for T&S using'
     write(*,'(i7,a)') ksw_lbc_uv,   ' - Open boundary conditions for U&V using'
     write(*,'(i7,a)') ksw_lbc_ssh,  ' - Open boundary conditions for SSH using'
     write(*,'(e12.4,a)') sst_relax,     ' - Relaxation coefficient for temperature [m/s]'
     write(*,'(e12.4,a)') sss_relax,     ' - Relaxation coefficient for salinity [m/s]'
     write(*,'(e12.4,a)') ldiff_ts,      ' - Lateral diffusion for temperature [m**2/s]'
     write(*,'(e12.4,a)') lvisc_2,       ' - Lateral  vicosity(2nd order)[m**2/s]'
     write(*,'(e12.4,a)') lvisc_4,       ' - Lateral  vicosity(4th order)[undim]'
     write(*,'(e12.4,a)') tsfrac_lat,    ' - fraction of salinity lateral diffusion'
     write(*,'(e12.4,a)') vdiff_ts_min,  ' - Vertical background diffusion coefficient for T [m**2/s]'
     write(*,'(e12.4,a)') vdiff_ts_max,  ' - Vertical max(top) diffusion coefficient for T [m**2/s]'
     write(*,'(e12.4,a)') vvisc_min,     ' - Vertical background viscosity coefficient [m**2/s]'
     write(*,'(e12.4,a)') vvisc_max,     ' - Vertical max(top) viscosity coefficient [m**2/s]'
     write(*,'(e12.4,a)') tsfrac_vert,   ' - fraction of salinity vertical diffusion'
     write(*,'(e12.4,a)') z_frac,        ' - Weight coefficient for lateral Z-Diffusion'
     write(*,'(e12.4,a)') r_frac,        ' - Weight coefficient for lateral R-Diffusion'
     write(*,'(e12.4,a)') gm_ratio,      ' - Weight coefficient for Gent&McWilliams transport'
     write(*,'(a,a)')  ' file with T-point sea-land mask: ', t_mask_file(1:len_trim (t_mask_file))
     write(*,'(a,a)')  '     file with bottom topography: ', bottom_topography_file(1:len_trim (bottom_topography_file))
 endif

 igrzts_surf  = min(IABS(ksw_ssbc),2) ! type of condition for T and S on sea surface
 igrzts_bot = 2                   ! type of condition for T and S on sea bottom

 dkft=sst_relax
 dkfs=sss_relax

 ! area mask initialization
   call gridcon(t_mask_file)
   if (rank .eq. 0) print *, "--------------------END OF GRIDCON----------------------"

!  setting vertical t-,w- grid levels
   call vgrid
   if (rank .eq. 0) print *, "--------------------END OF VGRID----------------------"

! define grid geographical coordinates, steps and coriolis parameters
   call basinpar
   if (rank .eq. 0) print *, "--------------------END OF BASINPAR----------------------"

   array4=0.0
   call prdstd(' ',bottom_topography_file,1,array4,lu,nx,ny,1, mmm,mm,nnn,nn,1,1,ierr)
   hhq_rest=dble(array4)
   call syncborder_real8(hhq_rest, 1)

   if(periodicity_x/=0) then
       call cyclize8_x(hhq_rest,nx,ny,1,mmm,mm)
   end if

   if(periodicity_y/=0) then
       call cyclize8_y(hhq_rest,nx,ny,1,nnn,nn)
   end if

!Initializing lateral diffusion
include 'latdiff.fi'

call depth_ave(amuv ,amuv2d, lu,0)
call depth_ave(amuv4,amuv42d,lu,0)

! solar shortwave heating penetrate divergence function.
call shortwave_divergence

!$omp parallel do private(m,n)
    do n=ny_start-1,ny_end+1
     do m=nx_start-1,nx_end+1
! set upper coefficients for temperature vertical diffusion
            anzt(m,n,1) = vdiff_ts_max*lu(m,n)
! set upper coefficients for momentum vertical viskosiity
            anzu(m,n,1) = vvisc_max*lu(m,n)
! set coefficients for temperature vertical diffusion
            anzt(m,n,2:nz) = vdiff_ts_min*lu(m,n)
!c set coefficients for velocity vertical viskosiity
            anzu(m,n,2:nz) = vvisc_min*lu(m,n)
     enddo
  enddo
!$omp end parallel do

endsubroutine ocean_model_parameters

!====================================================================
subroutine shortwave_divergence
use main_basin_pars
use mpi_parallel_tools
use basin_grid
use ocean_variables
implicit none
integer m,n,k
!-----------------------------------------------------------------------
!     Solar Shortwave energy penetrates below the ocean surface. Clear
!     water assumes energy partitions between two exponentials as
!     follows:
!
!     58% of the energy decays with a 35 cm e-folding scale
!     42% of the energy decays with a 23 m  e-folding scale
!
!     Paulson and Simpson (1977 Irradiance measurements in the upper
!                               ocean JPO 7, 952-956)
!     Also see ... Jerlov (1968 Optical oceanography. Elsevier)
!                  A General Circulation Model for Upper Ocean
!                  Simulaton (Rosati and Miyakoda JPO vol 18,Nov 1988)
!-----------------------------------------------------------------------
!
! => Shortwave penetration is a double exponential as follows:
!     parameter(rpart = 0.58,efold1 = 0.35,efold2 = 23.0)
! new approuch of using short wave radiation: upper part
! of about 60% added to heat flux and residual part
! of 40% of the energy decays with a 20 m e-folding scale
      real(8) rpart, efold1, efold2
      parameter(rpart = 0.58d0, efold1 = 0.35d0, efold2 = 23.0d0)
      real(8) swarg1,swarg2,pen1,pen2,sumk,spacesum,pointsum

      spacesum=0.0d0
      pointsum=0.0d0

      do n = bnd_y1+1, bnd_y2-1
       do m = bnd_x1+1, bnd_x2-1

        if(lu(m,n)>0.5) then

          pointsum=pointsum+1.0d0

          pen2=1.0d0
          sumk=0.0d0

          do k=1,nz
           swarg1 = -zw(k+1)*hhq_rest(m,n)/efold1
           swarg2 = -zw(k+1)*hhq_rest(m,n)/efold2
           pen1 = pen2
            if(k==nz) then
              pen2=0.0d0
            else
              pen2 = rpart*exp(swarg1) + (1.0d0-rpart)*exp(swarg2)
            endif
           divswrad(m,n,k) =(pen1 - pen2)/dz(k)/hhq_rest(m,n)
           sumk=sumk+pen1-pen2
          enddo

          do k=1,nz
!           divswrad(m,n,k) = divswrad(m,n,k)/sumk
           spacesum=spacesum+divswrad(m,n,k)*dz(k)*hhq_rest(m,n)
          enddo

        end if
       enddo
      enddo

      if (rank .eq. 0) then
          write(*,'(2x,a,f10.2)') 'sum of SW divergence coefficient =',spacesum
          write(*,'(2x,a,f10.2,a)')'for ',pointsum,' t-grid points'
      endif
endsubroutine shortwave_divergence
!=======================================================
subroutine ocinicond(start_type,path2ocp)
use main_basin_pars
use mpi_parallel_tools
use basin_grid
use ocean_variables
use ocean_bc

implicit none
character(256) fname
character*(*) path2ocp
integer m, n, k, ierr, lqp
integer start_type
real(8) density, hx2, hy2
real(4) array4(bnd_x1:bnd_x2,bnd_y1:bnd_y2,nz)

! initial conditions for temperature, salinity
      if (start_type==0) then

!read potential temperature
       array4=0.0
       call prdstd(path2ocp,'cptt.dat', 1,array4,lu,nx,ny,nz, mmm,mm,nnn,nn,1,nz,ierr)
       tt=dble(array4)

       call syncborder_real8(tt, nz)
!read salinity
       array4=0.0
       call prdstd(path2ocp,'cpss.dat', 1,array4,lu,nx,ny,nz, mmm,mm,nnn,nn,1,nz,ierr)
       ss=dble(array4)
       call syncborder_real8(ss, nz)

	 if(periodicity_x/=0) then
          call cyclize8_x(tt, nx,ny,nz,mmm,mm)
          call cyclize8_x(ss, nx,ny,nz,mmm,mm)
       end if

	 if(periodicity_y/=0) then
          call cyclize8_y(tt, nx,ny,nz,nnn,nn)
          call cyclize8_y(ss, nx,ny,nz,nnn,nn)
       end if

       ttp=tt
       ssp=ss

      else

!read potential temperature
       call rdstd8(path2ocp,'cptt8.dat', 1,tt ,lu,nx,ny,nz, mmm,mm,nnn,nn,1,nz,ierr)
       call rdstd8(path2ocp,'cptt8.dat', 2,ttp,lu,nx,ny,nz, mmm,mm,nnn,nn,1,nz,ierr)
!read salinity
       call rdstd8(path2ocp,'cpss8.dat', 1,ss ,lu,nx,ny,nz, mmm,mm,nnn,nn,1,nz,ierr)
       call rdstd8(path2ocp,'cpss8.dat', 2,ssp,lu,nx,ny,nz, mmm,mm,nnn,nn,1,nz,ierr)
!read zonal velocity
       call rdstd8(path2ocp,'cpuu8.dat', 1,uu ,llu,nx,ny,nz, mmm-1,mm,nnn  ,nn,1,nz,ierr)
       call rdstd8(path2ocp,'cpuu8.dat', 2,uup,llu,nx,ny,nz, mmm-1,mm,nnn  ,nn,1,nz,ierr)
!read meridional velocity
       call rdstd8(path2ocp,'cpvv8.dat', 1,vv ,llv,nx,ny,nz, mmm  ,mm,nnn-1,nn,1,nz,ierr)
       call rdstd8(path2ocp,'cpvv8.dat', 2,vvp,llv,nx,ny,nz, mmm  ,mm,nnn-1,nn,1,nz,ierr)

	 if(periodicity_x/=0) then
        call cyclize8_x(tt ,nx,ny,nz,mmm,mm)
        call cyclize8_x(ss ,nx,ny,nz,mmm,mm)
        call cyclize8_x(uu ,nx,ny,nz,mmm,mm)
        call cyclize8_x(vv ,nx,ny,nz,mmm,mm)
        call cyclize8_x(ttp,nx,ny,nz,mmm,mm)
        call cyclize8_x(ssp,nx,ny,nz,mmm,mm)
        call cyclize8_x(uup,nx,ny,nz,mmm,mm)
        call cyclize8_x(vvp,nx,ny,nz,mmm,mm)
	 end if

	 if(periodicity_y/=0) then
        call cyclize8_y(tt ,nx,ny,nz,nnn,nn)
        call cyclize8_y(ss ,nx,ny,nz,nnn,nn)
        call cyclize8_y(uu ,nx,ny,nz,nnn,nn)
        call cyclize8_y(vv ,nx,ny,nz,nnn,nn)
        call cyclize8_y(ttp,nx,ny,nz,nnn,nn)
        call cyclize8_y(ssp,nx,ny,nz,nnn,nn)
        call cyclize8_y(uup,nx,ny,nz,nnn,nn)
        call cyclize8_y(vvp,nx,ny,nz,nnn,nn)
	 end if

!  read turbulent kinetic energy and length scale

        call rdstd8(path2ocp,'cpq28.dat' , 1,q2  ,lu,nx,ny,nz+1, mmm,mm,nnn,nn,1,nz+1,ierr)
        call rdstd8(path2ocp,'cpq28.dat' , 2,q2p ,lu,nx,ny,nz+1, mmm,mm,nnn,nn,1,nz+1,ierr)
        call rdstd8(path2ocp,'cpq2l8.dat', 1,q2l ,lu,nx,ny,nz+1, mmm,mm,nnn,nn,1,nz+1,ierr)
        call rdstd8(path2ocp,'cpq2l8.dat', 2,q2lp,lu,nx,ny,nz+1, mmm,mm,nnn,nn,1,nz+1,ierr)

	 if(periodicity_x/=0) then
        call cyclize8_x(q2  ,nx,ny,nz+1,mmm,mm)
        call cyclize8_x(q2p ,nx,ny,nz+1,mmm,mm)
        call cyclize8_x(q2l ,nx,ny,nz+1,mmm,mm)
        call cyclize8_x(q2lp,nx,ny,nz+1,mmm,mm)
	 end if

	 if(periodicity_y/=0) then
        call cyclize8_y(q2  ,nx,ny,nz,nnn,nn)
        call cyclize8_y(q2p ,nx,ny,nz,nnn,nn)
        call cyclize8_y(q2l ,nx,ny,nz,nnn,nn)
        call cyclize8_y(q2lp,nx,ny,nz,nnn,nn)
	 end if

! read sea surface heights (internal mode) for present(1) and previous(2) time steps
       call rdstd8(path2ocp,'cpsshi8.dat', 1, ssh_i, lu,nx,ny,1, mmm,mm,nnn,nn,1,1,ierr)
       call rdstd8(path2ocp,'cpsshi8.dat', 2,sshp_i, lu,nx,ny,1, mmm,mm,nnn,nn,1,1,ierr)

       if(periodicity_x/=0) then
           call cyclize8_x( ssh_i, nx,ny,1,mmm,mm)
	       call cyclize8_x(sshp_i, nx,ny,1,mmm,mm)
       end if

       if(periodicity_y/=0) then
           call cyclize8_y( ssh_i, nx,ny,1,nnn,nn)
           call cyclize8_y(sshp_i, nx,ny,1,nnn,nn)
	   end if

! read sea surface heights (external mode) for present(1) previous(2) and pre_previous(3) time steps
       call rdstd8(path2ocp,'cpsshe8.dat', 1, ssh_e, lu,nx,ny,1, mmm,mm,nnn,nn,1,1,ierr)
       call rdstd8(path2ocp,'cpsshe8.dat', 2,sshp_e, lu,nx,ny,1, mmm,mm,nnn,nn,1,1,ierr)

! read barotropic velocity (external mode) for present(1) and previous(2) time steps
       call rdstd8(path2ocp,'cpube8.dat', 1, ubrtr_e,  llu,nx,ny,1, mmm-1,mm,nnn,nn,1,1,ierr)
       call rdstd8(path2ocp,'cpube8.dat', 2, ubrtrp_e, llu,nx,ny,1, mmm-1,mm,nnn,nn,1,1,ierr)
       call rdstd8(path2ocp,'cpvbe8.dat', 1, vbrtr_e,  llv,nx,ny,1, mmm,mm,nnn-1,nn,1,1,ierr)
       call rdstd8(path2ocp,'cpvbe8.dat', 2, vbrtrp_e, llv,nx,ny,1, mmm,mm,nnn-1,nn,1,1,ierr)

       if(periodicity_x/=0) then
           call cyclize8_x(ssh_e ,  nx,ny,1,mmm,mm)
	       call cyclize8_x(sshp_e,  nx,ny,1,mmm,mm)
           call cyclize8_x( ubrtr_e,nx,ny,1,mmm,mm)
           call cyclize8_x(ubrtrp_e,nx,ny,1,mmm,mm)
           call cyclize8_x( vbrtr_e,nx,ny,1,mmm,mm)
           call cyclize8_x(vbrtrp_e,nx,ny,1,mmm,mm)
	   end if

       if(periodicity_y/=0) then
           call cyclize8_y(ssh_e ,  nx,ny,1,nnn,nn)
           call cyclize8_y(sshp_e,  nx,ny,1,nnn,nn)
           call cyclize8_y( ubrtr_e,nx,ny,1,nnn,nn)
           call cyclize8_y(ubrtrp_e,nx,ny,1,nnn,nn)
           call cyclize8_y( vbrtr_e,nx,ny,1,nnn,nn)
           call cyclize8_y(vbrtrp_e,nx,ny,1,nnn,nn)
       end if

      endif

       q2   =max(q2   ,tur_var_min)
       q2p  =max(q2p  ,tur_var_min)
       q2l  =max(q2l  ,tur_var_min)
       q2lp =max(q2lp ,tur_var_min)

!initialize depth for internal mode
      call hh_init(hhq, hhqp, hhqn,    &
                   hhu, hhup, hhun,    &
                   hhv, hhvp, hhvn,    &
                   hhh, hhhp, hhhn,    &
                   ssh_i, sshp_i, hhq_rest)
!initialize depth for external mode
      call hh_init(hhq_e, hhqp_e, hhqn_e,    &
                   hhu_e, hhup_e, hhun_e,    &
                   hhv_e, hhvp_e, hhvn_e,    &
                   hhh_e, hhhp_e, hhhn_e,    &
                   ssh_e, sshp_e, hhq_rest)

  !compute depth mean and vertical velocity

    xxt=uu
    yyt=vv

    call depth_ave(xxt,uu2d ,llu,1)
    call depth_ave(yyt,vv2d ,llv,1)

    call depth_ave(uup,uup2d,llu,0)
    call depth_ave(vvp,vvp2d,llv,0)

!-----------------density definition-----------------------------------
      if (ksw_dens>0) then
!$omp parallel do private(m,n,k)
       do n=ny_start-1,ny_end+1
	  do m=nx_start-1,nx_end+1
            if(lu(m,n)>0.5) then
              do k=1,nz
               den(m,n,k)=density(tt(m,n,k),ss(m,n,k),FreeFallAcc*RefDen*hhq(m,n)*z(k))
           den_pot(m,n,k)=density(tt(m,n,k),ss(m,n,k),0.0d0)
              enddo
            endif
          enddo
         enddo
!$omp end parallel do
      endif

!--------------Rayleigh friction initialization
!$omp parallel do private(m,n,k, hx2, hy2)
       do n=ny_start,ny_end
	  do m=nx_start,nx_end
            if(lu(m,n)>0.5) then
              hx2= ( ((hhq_rest(m+1,n)-hhq_rest(m  ,n))/dxt(m  ,n))**2 * dble(lcu(m  ,n))    &
                    +((hhq_rest(m  ,n)-hhq_rest(m-1,n))/dxt(m-1,n))**2 * dble(lcu(m-1,n)) )/dble(lcu(m,n)+lcu(m-1,n))
              hy2= ( ((hhq_rest(m,n+1)-hhq_rest(m,n  ))/dyt(m,n  ))**2 * dble(lcv(m,n  ))    &
                    +((hhq_rest(m,n  )-hhq_rest(m,n-1))/dyt(m,n-1))**2 * dble(lcv(m,n-1)) )/dble(lcv(m,n)+lcv(m,n-1))
              r_diss(m,n)=r_fric*dsqrt(hx2+hy2)
            endif
          enddo
         enddo
!$omp end parallel do
       call syncborder_real8(r_diss, 1)

       if(periodicity_x/=0) then
           call cyclize8_x( r_diss, nx,ny,1,mmm,mm)
	   end if

       if(periodicity_y/=0) then
           call cyclize8_y( r_diss, nx,ny,1,nnn,nn)
	   end if

endsubroutine ocinicond
!=======================================================
! reading control point for ice model
 subroutine icinicond(start_type,path2ocp,nstep)
 use main_basin_pars
 use mpi_parallel_tools
 use basin_grid
 use ocean_variables

 implicit none

  integer  start_type, k, ierr, nstep
  character*(*) path2ocp

      if(start_type>0.and.nstep>0) then    !if model is running from control point
! reading hice - ice layer thicknesses=>
       call rdstd8(path2ocp,'cphice8.dat',1,hice,lu,nx,ny,mgrad, mmm,mm,nnn,nn,1,mgrad,ierr)
       call syncborder_real8(hice, mgrad)

! reading aice - ice square
       call rdstd8(path2ocp,'cpaice8.dat',1,aice,lu,nx,ny,mgrad, mmm,mm,nnn,nn,1,mgrad,ierr)
       call syncborder_real8(aice, mgrad)

! reading hsnow - snow thickness
       call rdstd8(path2ocp,'cphsnow8.dat',1,hsnow,lu,nx,ny,mgrad, mmm,mm,nnn,nn,1,mgrad,ierr)
       call syncborder_real8(hsnow, mgrad)

! read zonal ice velocity
       call rdstd8(path2ocp,'cpuice8.dat', 1,uice, llu,nx,ny,1, mmm-1,mm,nnn  ,nn,1,1,ierr)
       call syncborder_real8(uice, 1)

! read meridional ice velocity
       call rdstd8(path2ocp,'cpvice8.dat', 1,vice, llv,nx,ny,1, mmm  ,mm,nnn-1,nn,1,1,ierr)
       call syncborder_real8(vice, 1)

! read sigma1 ice stress
       call rdstd8(path2ocp,'cpsig18.dat', 1,ice_stress11,lu,nx,ny,1, mmm  ,mm,nnn,nn,1,1,ierr)
       call syncborder_real8(ice_stress11, 1)

! read sigma2 ice deformation rate
       call rdstd8(path2ocp,'cpsig28.dat', 1,ice_stress22,lu,nx,ny,1, mmm  ,mm,nnn,nn,1,1,ierr)
       call syncborder_real8(ice_stress22, 1)

! read sigma12 ice deformation rate
       call rdstd8(path2ocp,'cpsig128.dat', 1,ice_stress12,luh,nx,ny,1,  mmm-1,mm,nnn-1,nn,1,1,ierr)
       call syncborder_real8(ice_stress12, 1)

      end if

	  aice0=1.0d0

       do k=1,mgrad
        aice0(:,:)=aice0(:,:)-aice(:,:,k)
	   end do

	   aice0=min(max(aice0,0.0d0),1.0d0)
       call syncborder_real8(aice0, 1)

       if(periodicity_x/=0) then
        call cyclize8_x( aice ,nx,ny,mgrad,mmm,mm)
        call cyclize8_x( aice0,nx,ny,1,mmm,mm)
        call cyclize8_x(  hice,nx,ny,mgrad,mmm,mm)
        call cyclize8_x( hsnow,nx,ny,mgrad,mmm,mm)
        call cyclize8_x(  uice,nx,ny,1,mmm,mm)
        call cyclize8_x(  vice,nx,ny,1,mmm,mm)
        call cyclize8_x(ice_stress11,nx,ny,1,mmm,mm)
        call cyclize8_x(ice_stress22,nx,ny,1,mmm,mm)
        call cyclize8_x(ice_stress12,nx,ny,1,mmm,mm)
       end if

       if(periodicity_y/=0) then
        call cyclize8_y( aice ,nx,ny,mgrad,nnn,nn)
        call cyclize8_y( aice0,nx,ny,1,nnn,nn)
        call cyclize8_y(  hice,nx,ny,mgrad,nnn,nn)
        call cyclize8_y( hsnow,nx,ny,mgrad,nnn,nn)
        call cyclize8_y(  uice,nx,ny,1,nnn,nn)
        call cyclize8_y(  vice,nx,ny,1,nnn,nn)
        call cyclize8_y(ice_stress11,nx,ny,1,nnn,nn)
        call cyclize8_y(ice_stress22,nx,ny,1,nnn,nn)
        call cyclize8_y(ice_stress12,nx,ny,1,nnn,nn)
       end if
endsubroutine icinicond
!======================================================================
subroutine ptinicond(path2ocp)
use main_basin_pars
use mpi_parallel_tools
use basin_grid
use ocean_variables
implicit none

integer ierr
character*(*) path2ocp

! initial conditions for pass_tracer

      call rdstd8(path2ocp,'cppt8.dat', 1, pass_tracer,lu ,nx,ny,nz, mmm,mm,nnn,nn,1,nz,ierr)
      call syncborder_real8(pass_tracer, nz)

	if(periodicity_x/=0) then
       call cyclize8_x(pass_tracer,nx,ny,nz,mmm,mm)
      endif

	if(periodicity_y/=0) then
       call cyclize8_y(pass_tracer,nx,ny,nz,nnn,nn)
      endif

      igrzpt_surf=2
      igrzpt_bot =2

endsubroutine ptinicond

subroutine test_init
    use main_basin_pars
    use mpi_parallel_tools
    use basin_grid
    use ocean_variables
    implicit none

    integer :: m, n

    ssh_i = -1
    do m = nx_start, nx_end
        do n = ny_start, ny_end
            ssh_i(m, n) = rank
        enddo
    enddo


end subroutine
