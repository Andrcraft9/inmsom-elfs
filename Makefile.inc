MAKE = make

#Translator for Cluster
#FC        = mpif90

FC        = ifort 

# Sequental options
#FCFLAGS    = -O3 -I./Inc

#OMP parrallel options
FCFLAGS    = -O3 -openmp -auto -stack_temps -I./Inc

RM        = rm -rf
