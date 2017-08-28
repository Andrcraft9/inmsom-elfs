module atm_forcing
    implicit none
    include 'atmforcing.fi'

    integer ind_change_water,   &    !  index of data number change
            ind_change_slpr

    integer num_rec_water,      &    !  number of record in file
            num_rec_slpr

    real(4),allocatable:: a_wflux(:,:),       &   !precipitation-evaporation[m/s]
                          a_slpr(:,:)             !pressure at sea surface

endmodule atm_forcing
