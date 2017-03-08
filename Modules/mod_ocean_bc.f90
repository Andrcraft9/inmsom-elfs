module ocean_bc
implicit none

include 'oceanbc.fi'

integer numb_of_lb,        &   !number of liquid boundaries (calculated)   
        numb_of_lqp            !number of liquid points on t-grid(calculated)
	  
integer ind_change_sst,     &    !  index of data number change
        ind_change_sss,     &           
        ind_change_runoff,  &    
        ind_change_tlbc,    &      
        ind_change_slbc,    &   
        ind_change_ulbc,    &   
        ind_change_hlbc

integer num_rec_sst,     &    !  index of data number change
        num_rec_sss,     &           
        num_rec_runoff,  &    
        num_rec_tlbc,    &      
        num_rec_slbc,    &   
        num_rec_ulbc,    &   
        num_rec_hlbc

real(4), allocatable:: tlqbw(:,:),slqbw(:,:),      &         !t and s boundary values on liquid walls
                       ulqbw(:,:),vlqbw(:,:),      &         !u and v boundary values on liquid walls
                       sshlqbw(:)
real(4), allocatable:: sst_obs(:,:),       &       !Observed SST [°C]
                       sss_obs(:,:),       &       !Observed SSS [psu]
                       runoff(:,:)                 !River runoff [kg/m^2/s]

integer, allocatable:: lqpx(:),lqpy(:)                       !grid coordinates of liquid pointes(calculated)
character, allocatable:: index_of_lb(:)


endmodule ocean_bc
