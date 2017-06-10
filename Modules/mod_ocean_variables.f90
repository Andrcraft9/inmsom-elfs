!input physical parameters of the run
module key_switches
implicit none

integer ksw_ts,         &     !Temperature and salinity equation solving (0 - no, 1 - yes)
        ksw_age,        &     !Ideal age equation solving (0 - no, 1 - yes)
        ksw_pt,         &     !Passive tracer equation solving (0 - no, 1 - yes)
        ksw_uv,         &     !Momentum equation solving (0 - no, 1 - yes)
        ksw_lat,        &     !Lateral 2nd order mix parametrization (0 - constant coeff, 1 - Smagorinski)
        ksw_lat4,       &     !Lateral 4nd order momentum mix parametrization (0 - no, 1 - yes)
        ksw_vert,       &     !vertical mix parametrization (0 - constant coeff, 1 - Pacanowski&Philander)
        ksw_dens,       &     !pressure gradient computation (0 - no (constant density), 1 - yes)
        ksw_ice_th,     &     !sea ice thermodynamics using (0 - no, 1 - yes)
        ksw_ice_tran,   &     !sea ice transport using (0 - no, 1 - yes)
        ksw_ice_dyn,    &     !sea ice dynamics using (0 - no, 1 - yes)
        ksw_ssbc,       &     !Type of surface boundary conditions (1 - surface T&S and wind stress are set; 2 - T&S fluxes and wind stress are set; 3 - T&S fluxes and wind stress are computed
        ksw_wflux,      &     !normalize global mean salt balance (0 - no, 1 - normalize water flux, 2 - normalize salinity flux)
        ksw_lbc_ts,     &     !open boundary conditions for T&S using (0 - no, 1 - yes)
        ksw_lbc_uv,     &     !open boundary conditions for U&V using (0 - no, 1 - yes)
        ksw_lbc_ssh           !open boundary conditions for SSH using (0 - no, 1 - yes)

real(8) sst_relax,      &     !Relaxation coefficient for temperature [m/s]
        sss_relax,      &     !Relaxation coefficient for salinity [m/s]
        ldiff_ts,       &     !lateral diffusion for temperature [m**2/s]
        lvisc_2,        &     !lateral  vicosity(2nd order)[m**2/s]
        lvisc_4,        &     !lateral  vicosity(4th order) [undim]
     tsfrac_lat,        &     !fraction of salinity lateral diffusion due to one for temperature
        vdiff_ts_min,   &     !vertical background diffusion coefficient for T [m**2/s]
        vdiff_ts_max,   &     !vertical max(top) diffusion coefficient for T [m**2/s]
        vvisc_min,      &     !vertical background viscous coefficient [m**2/s]
        vvisc_max,      &     !vertical max(top) viscous coefficient [m**2/s]
      tsfrac_vert,      &     !fraction of salinity lateral diffusion due to one for temperature
        z_frac,         &     !weight coefficient for lateral Z-Diffusion
        r_frac,         &     !weight coefficient for lateral R-Diffusion
        gm_ratio              !weight coefficient for Gent&McWilliams transport

endmodule key_switches

!-------------module for description common ogcm variables and task control parameters---------
module ocean_variables
use key_switches
implicit none
!barotropic dynamics arrays
real(8),allocatable:: ssh_i(:,:),     &  !sea surface height (SSH) at current  time step [m] (internal mode)
                      sshp_i(:,:),    &  !sea surface height (SSH) at previous time step [m] (internal mode)
                       pgrx(:,:),     &  !pressure gradient x-component for RHS
                       pgry(:,:),     &  !pressure gradient y-component for RHS
                    ubrtr_i(:,:),     &  !barotropic velocity      zonal[m/s] at current time step (internal mode)
                    vbrtr_i(:,:),     &  !barotropic velocity meridional[m/s] at current time step (internal mode)
                     RHSx2d(:,:),     &  !x-component of external force(barotropic)
                     RHSy2d(:,:)         !y-component of external force(barotropic)

real(8), allocatable:: ssh_e(:,:),   &  !sea surface height (SSH) at current  time step [m] (external mode)
                      sshp_e(:,:),   &  !sea surface height (SSH) at previous time step [m] (external mode)
                     ubrtr_e(:,:),   &  !barotropic velocity      zonal[m/s]                       (external mode)
                    ubrtrp_e(:,:),   &  !barotropic velocity      zonal[m/s] at previous time step (external mode)
                     vbrtr_e(:,:),   &  !barotropic velocity meridional[m/s]                       (external mode)
                    vbrtrp_e(:,:)       !barotropic velocity meridional[m/s] at previous time step (external mode)

!3d dynamics arrays
real(8),allocatable:: xxt(:,:,:),   &  !auxiliary array 1
                      yyt(:,:,:)       !auxiliary array 2

! sea surface boundary condition
real(8), allocatable:: wf_tot(:,:)                !total water flux

real(8), allocatable:: BottomFriction(:,:),    &    !Bottom friction rate (m/s)
                               r_diss(:,:)          !Rayleigh friction scale (1/s)

real(8), allocatable:: amuv2d(:,:),     &    !depth mean lateral viscosity
                      amuv42d(:,:),     &    !depth mean lateral viscosity
                     r_vort2d(:,:),     &    !relative vorticity of depth mean velocity
                   stress_t2d(:,:),     &    !Horizontal tension tensor component (barotropic)
                   stress_s2d(:,:),     &    !Horizontal shearing tensor component(barotropic)
             RHSx2d_tran_disp(:,:),     &    !dispersion x-component of external force(barotropic)
             RHSy2d_tran_disp(:,:),     &    !dispersion y-component of external force(barotropic)
             RHSx2d_diff_disp(:,:),     &    !dispersion x-component of external force(barotropic)
             RHSy2d_diff_disp(:,:)           !dispersion y-component of external force(barotropic)

endmodule ocean_variables
!-------------end module for description common ogcm variables and task control parameters---------
