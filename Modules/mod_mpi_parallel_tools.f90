!---------------------module for definition of array dimensions and boundaries-----------------
module mpi_parallel_tools
implicit none

integer nx_start,    &   !first significant point in x-direction
        nx_end,      &   !last  significant point in x-direction
        ny_start,    &   !first significant point in y-direction       
        ny_end           !last  significant point in y-direction

integer bnd_x1,      &   !left   array boundary in x-direction
        bnd_x2,      &   !right  array boundary in x-direction
        bnd_y1,      &   !bottom array boundary in y-direction
        bnd_y2           !top    array boundary in y-direction

endmodule mpi_parallel_tools
!---------------------end module for definition of array dimensions and boundaries-----------------
