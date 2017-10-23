#
# Include file for g95 Fortran compiler on Linux
#


F90C ?= g95
F90FLAGS ?=
LD ?= $(F90C)
LDFLAGS ?= 

# OpenMP is not availble



ifdef DEBUG
  F90FLAGS += -g -fbounds-check -ftrace=full
else
  F90FLAGS += -O3
endif

ifeq ($(PRECISION),double)
  F90FLAGS += -fdefault-real-8
endif

ifeq ($(FORMAT),big_endian)
  F90FLAGS += -fendian=big
endif  

