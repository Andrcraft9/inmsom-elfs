include Makefile.inc


SRCCONTROL =	\
	Control/dyncall.f90      \
	Control/init_arrays.f90  \
	Control/init_pars.f90    \
	Control/output.f90   

SRCSERVICE =	    \
	Service/basin_parameters.f90      \
	Service/grid_construction.f90     \
	Service/input_output_data.f90     \
	Service/oc_algorithms.f90         \
	Service/read_write_parameters.f90 \
	Service/rw_ctl_file.f90           \
	Service/time_tools.f90                 
         

SRCFUNCTION =	\
	Function/depth.f90   \
	Function/mixing.f90  \
	Function/vel_ssh.f90   


SRCMODULES = 	\
	Modules/mod_basin_grid.f90           \
	Modules/mod_main_basin_pars.f90      \
	Modules/mod_mpi_parallel_tools.f90   \
	Modules/mod_ocean_variables.f90      \
	Modules/mod_rec_length.f90           \
	Modules/mod_time_integration.f90 

all: inmsom clean

inmsom:
#order is important
	$(FC) -o inmsom $(FCFLAGS) $(SRCMODULES) inmsom_head.f90  $(SRCFUNCTION) $(SRCCONTROL)  $(SRCSERVICE) 	
clean:
	$(RM) *.o *.mod 
