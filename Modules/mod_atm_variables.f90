module atm_forcing
    implicit none
    include 'atmforcing.fi'

    integer ind_change_heat,   &    !  index of data number change
            ind_change_water,  &
            ind_change_stress, &
            ind_change_rad,    &
            ind_change_prec,   &
            ind_change_tatm,   &
            ind_change_qatm,   &
            ind_change_wind,   &
            ind_change_slpr

    integer num_rec_heat,   &    !  number of record in file
            num_rec_water,  &
            num_rec_stress, &
            num_rec_rad,    &
            num_rec_prec,   &
            num_rec_tatm,   &
            num_rec_qatm,   &
            num_rec_wind,   &
            num_rec_slpr

    real(4),allocatable:: a_hflux(:,:),       &   !heat balance [w/m**2]
                          a_swrad(:,:),       &   !sw radiation balance[w/m**2]
                          a_wflux(:,:),       &   !precipitation-evaporation[m/s]
                          a_stress_x(:,:),    &   !zonal wind stress[pA=n/m**2]
                          a_stress_y(:,:),    &   !meridional wind stress[pA=n/m**2]
                          a_slpr(:,:),        &   !pressure at sea surface
                          a_lwr(:,:),         &   !dw-lw-rad[w/m**2]
                          a_swr(:,:),         &   !dw-sw-rad[w/m**2]
                          a_rain(:,:),        &   !precipit[m/s]
                          a_snow(:,:),        &   !precipit[m/s]
                          a_tatm(:,:),        &   !temp of atmosphere[�c]
                          a_qatm(:,:),        &   !humidity [g/kg]
                          a_uwnd(:,:),        &   !u-wind speed[m/s]
                          a_vwnd(:,:)             !v-wind speed[m/s]
endmodule atm_forcing
