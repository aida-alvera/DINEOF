#==========================================================================
#
# Include file for gfortran Fortran compiler on Windows (mingw)
# For more information about mingw:
# http://www.mingw.org
# http://gcc.gnu.org/wiki/GFortranBinariesWindows
#
#==========================================================================


F90C := gfortran
F90FLAGS :=
LD := $(F90C)
LDFLAGS := 

# necessary for gfortran 4.1, but default since 4.2
F90FLAGS += -frecord-marker=4


ifdef DEBUG
  F90FLAGS += -g -fbounds-check
else
  F90FLAGS += -O3 -ffast-math
endif

ifeq ($(PRECISION),double)
  F90FLAGS += -fdefault-real-8
endif

ifeq ($(FORMAT),big_endian)
  F90FLAGS += -fconvert=big-endian -frecord-marker=4
endif  

ifdef STATIC
  LDFLAGS += -static
endif

