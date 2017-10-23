#
# Include file for MIPSpro Fortran compiler on IRIX
#


F90C := f90
F90FLAGS := 
LD := $(F90C)
LDFLAGS := 

ifeq ($(FORMAT),little_endian)
error:
	echo "Error: machine format little_endian not available."; exit 1
endif  

ifdef DEBUG
  F90FLAGS += -g
else
  F90FLAGS += -O3
endif

#
# Library locations
#

INCLUDES =

DINEOF_LIBRARIES =  -L$(HOME)/lib -larpack_x86_64  -llapack  -lblas  -lnetcdf

CROSSVAL_LIBRARIES = -lnetcdf
