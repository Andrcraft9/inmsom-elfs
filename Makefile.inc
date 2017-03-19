MAKE = make

#Translator for Cluster
FC        = mpif90

#FC        = gfortran

# Sequental options
FCFLAGS    = -w -ffree-line-length-0 -O3 -I./Inc
#FCFLAGS    = -w -ffree-line-length-0 -I./Inc

#OMP parrallel options
#FCFLAGS    = -O3 -openmp -auto -stack_temps -I./Inc

RM        = rm -rf
