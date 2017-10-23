#
# Include file for Intel Fortran compiler on Linux
#


F90C ?= ifort
F90FLAGS ?= -implicitnone
LD ?= $(F90C)
LDFLAGS ?= 

# avoid stack size problem
# http://software.intel.com/en-us/articles/intel-fortran-compiler-increased-stack-usage-of-80-or-higher-compilers-causes-segmentation-fault/

F90FLAGS += -heap-arrays 
PROFILING_F90FLAGS ?= -p
PROFILING_LDFLAGS ?= -p

PIC_F90FLAGS=-fPIC
PIC_CFLAGS=-fPIC

ifdef OPENMP
  F90FLAGS += -openmp
  LDFLAGS += -openmp
endif

ifdef DEBUG
  F90FLAGS += -g -check all -traceback
else
  F90FLAGS += -vec-report0 -O3 
endif

ifeq ($(PRECISION),double)
  F90FLAGS += -r8
endif

ifeq ($(FORMAT),big_endian)
  F90FLAGS += -convert big_endian
endif  

