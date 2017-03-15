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
real(8),allocatable::  uu(:,:,:),   &  !     zonal velocity [m/s]
                       uup(:,:,:),  &  !     zonal velocity [m/s] at previous time step
                       vv(:,:,:),   &  !meridional velocity [m/s]
                       vvp(:,:,:),  &  !meridional velocity [m/s] at previous time step
                       ww(:,:,:),   &  !  vertical velocity in sigma-coord [m/s]
                       tt(:,:,:),   &  !potential temperature[°C]
                       ttp(:,:,:),  &  !potential temperature[°C] at previous time step
                       ss(:,:,:),   &  !salinity [psu]
                       ssp(:,:,:),  &  !salinity [psu] at previous time step
                      den(:,:,:),   &  !in-situ density  [kg/m^3]
                  den_pot(:,:,:),   &  !potential density  [kg/m^3]
                   RHSx3d(:,:,:),   &  !x-component of external force(baroclinic)
                   RHSy3d(:,:,:),   &  !y-component of external force(baroclinic)
              RHSx3d_tran(:,:,:),   &  !x-component of external force(baroclinic)
              RHSy3d_tran(:,:,:),   &  !y-component of external force(baroclinic)
                      age(:,:,:),   &  !water ideal age [days]
                      agep(:,:,:),  &  !water ideal age [days] at previous time step
                      xxt(:,:,:),   &  !auxiliary array 1
                      yyt(:,:,:)       !auxiliary array 2

real(8),allocatable::  q2(:,:,:),   &  !turbulent kinetic energy (m/s)^2
                      q2p(:,:,:),   &  !turbulent kinetic energy at the previous time step (m/s)^2
                      q2l(:,:,:),   &  !turbulent kinetic energy by turbulent length scale (m/s)^2*m
                     q2lp(:,:,:)       !turbulent kinetic energy by turbulent length scale at the previous time step (m/s)^2*m

real(8),allocatable:: RHS_q2 (:,:,:),   &  !turbulent kinetic energy (m/s)^2
                      RHS_q2l(:,:,:)       !turbulent kinetic energy by turbulent length scale (m/s)^2*m

real(8),allocatable:: uu2d(:,:),    &
                      vv2d(:,:),    &
                     uup2d(:,:),    &
                     vvp2d(:,:),    &
               RHSx2d_tran(:,:),    &  !x-component of external force(baroclinic)
               RHSy2d_tran(:,:)        !y-component of external force(baroclinic)

real(8),allocatable:: stress_t(:,:,:),    &     !Horizontal tension tensor component
                      stress_s(:,:,:),    &     !Horizontal shearing tensor component
                        r_vort(:,:,:)           !Part of relative vorticity

!surface and bottom boundary condition types                ! igrz[T,S]= 1 : f = f0 
integer,allocatable:: igrzts_surf(:,:), igrzts_bot(:,:)     !          = 2 : df/dz = f0

!Passive tracer arrays
real(8),allocatable:: pass_tracer(:,:,:),   &     !Passive tracer
                      pass_tracerp(:,:,:),  &     !Passive tracer at previous time step
                        pt_forc_surf(:,:),  &     !Passive tracer surface forcing
                        pt_forc_bot(:,:),   &     !Passive tracer bottom forcing
                        pt_diff_x(:,:,:),   &     !Passive x-diffusion [m^2/s]
                        pt_diff_y(:,:,:)          !Passive y-diffusion [m^2/s]

integer,allocatable:: igrzpt_surf(:,:), igrzpt_bot(:,:)  !Passive tracer boundary condition type

! coefficients of viscosity  and diffusion
! lateral viscosity and diffusion coefficients for t, s, u, v:
! horizontal:
real(8),allocatable::  amts(:,:,:),     &   !T lateral diffusion in T-points [m^2/s]
                       amuv(:,:,:),     &   !U and V lateral viscosity in T-points [m^2/s]
                       amuv4(:,:,:)         !U and V 4-th order lateral viscosity in T-points [m^4/s]^(1/2)

! lateral diffusion function for t, s:
real(8),allocatable::  slrx(:,:,:),   &     !universal isopycnal diffusion slope in x-direction
                       slry(:,:,:),   &     !universal isopycnal diffusion slope in y-direction
                       slzx(:,:,:),   &     !universal horizontal diffusion slope in x-direction
                       slzy(:,:,:)          !universal horizontal diffusion slope in y-direction 

! vertical viscous and diffusion functions
real(8),allocatable:: rit(:,:,:),     &     !Richardson number
                     anzt(:,:,:),     &     !T vertical diffusion [m^2/s]
                     anzu(:,:,:)            !U and V vertical viscosity [m^2/s]

! sea surface boundary condition
real(8), allocatable:: tflux_surf(:,:),      &       !total surface heat flux [°C*m/s]
                       tflux_bot(:,:),       &       !total bottom heat flux [°C*m/s] 
                       sflux_surf(:,:),      &       !total surface salt flux [psu*m/s]
                       sflux_bot(:,:),       &       !total bottom salt flux [psu*m/s]
                   surf_stress_x(:,:),       &       !wind      zonal stress per water density [m^2/s^2]
                   surf_stress_y(:,:),       &       !wind meridional stress per water density [m^2/s^2]
                    bot_stress_x(:,:),       &       !bottom    zonal stress per water density [m^2/s^2]
                    bot_stress_y(:,:),       &       !bottom meridional stress per water density [m^2/s^2]
                      divswrad(:,:,:),       &       !shortwave radiation divergence coefficients
                            dkft(:,:),       &       !relaxation coefficient for SST, [m/s]
                            dkfs(:,:),       &       !relaxation coefficient for SSS, [m/s]
                        sensheat(:,:),       &       !sensible heat flux
                         latheat(:,:),       &       !latent heat flux
                          lw_bal(:,:),       &       !longwave radiation balance
                          sw_bal(:,:),       &       !shortwave radiation balance
                          hf_tot(:,:),       &       !total heat flux
                          wf_tot(:,:)                !total water flux

real(8), allocatable:: q2_surf(:,:),      &     !surface boundary condition for q2
                        q2_bot(:,:),      &     !bottom  boundary condition for q2
                      q2l_surf(:,:),      &     !surface boundary condition for q2l
                       q2l_bot(:,:)             !bottom  boundary condition for q2l

!Atmospheric arrays for bulk-formulae
real(8),allocatable:: tatm(:,:),   &    !Air temperature, [°C]
                      qatm(:,:),   &    !Air humidity, [kg/kg]
                      rain(:,:),   &    !rain, [kg/m^2/s]
                      snow(:,:),   &    !snow, [kg/m^2/s]
                      wind(:,:),   &    !Wind speed module, [m/s]
                       lwr(:,:),   &    !Downward  longwave radiation, [W/m^2]
                       swr(:,:),   &    !Downward shortwave radiation, [W/m^2]
                      slpr(:,:),   &    !Sea level pressure, [Pa]
                      uwnd(:,:),   &    !Zonal      wind speed, [m/s]
                      vwnd(:,:),   &    !Meridional wind speed, [m/s]
                      taux(:,:),   &    !Zonal      wind stress, [Pa]
                      tauy(:,:)         !Meridional wind stress, [Pa]
real(8), allocatable::  BottomFriction(:,:),    &    !Bottom friction rate (m/s)
                                r_diss(:,:)          !Rayleigh friction scale (1/s)

! array description for ice.
real(8), allocatable:: hice(:,:,:),    &       !Ice mass, [m]
                       aice(:,:,:),    &       !Ice compactness, [0-1]
                      aice0(:,:),      &       !Open water compactness, [0-1]
                      hsnow(:,:,:),    &       !Snow mass, [m]
                       tice(:,:,:),    &       !Ice temperature, [C]
                      tsnow(:,:,:),    &       !Snow temperature, [C]
                    dhsnowt(:,:),      &       !Total snow mass change, [m]
                     dhicet(:,:),      &       !Total ice mass change, [m]
                      swice(:,:),      &       !under-ice shortwave radiation, [W/m^2]
                 heatice2oc(:,:),      &       !heat disbalance returning to ocean, [W/m^2]
                 ice_stress11(:,:),    &       !Stress tensor 1-1 component
                 ice_stress22(:,:),    &       !Stress tensor 2-2 component
                 ice_stress12(:,:),    &       !Stress tensor 1-2 component
                         uice(:,:),    &       !Ice zonal velocity, [m/s]
                         vice(:,:)             !Ice meridional velocity, [m/s]

real(8), allocatable:: Flux_tem_x(:,:,:), &       !Total temperature flux along x-direction
                       Flux_tem_y(:,:,:), &       !Total temperature flux along y-direction
                       Flux_sal_x(:,:,:), &       !Total   salinity  flux along x-direction
                       Flux_sal_y(:,:,:)          !Total   salinity  flux along y-direction


real(8), allocatable::   amuv2d(:,:),     &    !depth mean lateral viscosity
                        amuv42d(:,:),     &    !depth mean lateral viscosity
                       r_vort2d(:,:),     &    !relative vorticity of depth mean velocity
                     stress_t2d(:,:),     &    !Horizontal tension tensor component (barotropic)
                     stress_s2d(:,:),     &    !Horizontal shearing tensor component(barotropic)
                    RHSx2d_tran_disp(:,:),     &    !dispersion x-component of external force(barotropic)
                    RHSy2d_tran_disp(:,:),     &    !dispersion y-component of external force(barotropic)
                    RHSx2d_diff_disp(:,:),     &    !dispersion x-component of external force(barotropic)
                    RHSy2d_diff_disp(:,:)           !dispersion y-component of external force(barotropic)

real(8),allocatable:: tt_calc(:,:,:),     &
                      ss_calc(:,:,:),     &
                      uu_calc(:,:,:),     &
                      vv_calc(:,:,:),     &
                     sfl_calc(:,:,:),     &
                     ssh_calc(:,:),       &
                     txo_calc(:,:),       &
                     tyo_calc(:,:),       &
                    uwnd_calc(:,:),       &
                    vwnd_calc(:,:),       &
                  Fltx_calc(:,:,:),       &
                  Flty_calc(:,:,:),       &
                  Flsx_calc(:,:,:),       &
                  Flsy_calc(:,:,:)

integer meancalc             !calculator for time mean output

endmodule ocean_variables
!-------------end module for description common ogcm variables and task control parameters---------
