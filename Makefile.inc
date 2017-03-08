MAKE = make

#Translator for Cluster
FC        = mpif90

#FC        = gfortran

# Sequental options
#FCFLAGS    = -O3 -I./Inc

#OMP parrallel options
FCFLAGS    = -w -ffree-line-length-0 -I./Inc

RM        = rm -rf
