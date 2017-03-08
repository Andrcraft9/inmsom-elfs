module atm2oc_interpol
implicit none

!parameters for matrix interpolation from atm to ocean grid
 real(8), allocatable::   wght_mtrx_a2o(:,:,:)  !nonzero matrix elements for a2o interpolation

 integer, allocatable::   i_input_a2o(:,:,:),    &    !i-numbers of matrix elements
                          j_input_a2o(:,:,:)          !j-numbers of matrix elements

endmodule atm2oc_interpol