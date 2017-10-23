#
# Include file for PGI Fortran compiler on Linux
#


F90C ?= pgf90
F90FLAGS ?=
LD ?= $(F90C)
LDFLAGS ?= 

#PROFILING_F90FLAGS ?= -Mprof=dwarf
#PROFILING_LDFLAGS ?= -Mprof=dwarf
# -Mprof=hwcts  -Mprof=func -Mprof=lines -Mprof=time
PROFILING_F90FLAGS ?= -Minfo=ccff
PROFILING_LDFLAGS ?= $(PROFILING_F90FLAGS)

PIC_F90FLAGS=-fPIC
PIC_CFLAGS=-fPIC

ifdef OPENMP
  F90FLAGS += -mp
  LDFLAGS += -mp
endif

ifdef DEBUG
  F90FLAGS += -g -C -Mchkptr
else
#  F90FLAGS += -u -fastsse -Mipa=fast
  F90FLAGS += -O3 -Mflushz
endif

ifeq ($(PRECISION),double)
  F90FLAGS += -r8
endif

ifeq ($(FORMAT),big_endian)
  F90FLAGS += -byteswapio
endif  

ifdef STATIC
  F90FLAGS += -Bstatic
endif


