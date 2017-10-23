
ifdef PIC
  F90FLAGS += $(PIC_F90FLAGS)
  CFLAGS += $(PIC_CFLAGS)
endif


#
# Library locations
#

# If all libraries are in one folder

INCDIR ?=
LIBDIR ?=

# netCDF configuration

NETCDF_CONFIG ?= nf-config
NETCDF_INCDIR ?= $(INCDIR)
NETCDF_LIBDIR ?= $(LIBDIR)
NETCDF_LIB ?= -lnetcdf

# LAPACK configuration

LAPACK_LIBDIR ?= $(LIBDIR)
LAPACK_LIB ?= -llapack

# BLAS configuration

BLAS_LIBDIR ?= $(LIBDIR)
BLAS_LIB ?= -lblas

# ARPACK configuration

ARPACK_LIBDIR ?= $(LIBDIR)
ARPACK_LIB ?= -larpack

# Extra parameters

EXTRA_F90FLAGS ?=
EXTRA_LDFLAGS ?=



# LAPACK
# use LAPACK_LIBDIR only if it is non-empty

ifneq ($(strip $(LAPACK_LIBDIR)),)
  LIBS += -L$(LAPACK_LIBDIR)
endif
LIBS += $(LAPACK_LIB)

# BLAS
# use BLAS_LIBDIR only if it is non-empty

ifneq ($(strip $(BLAS_LIBDIR)),)
  LIBS += -L$(BLAS_LIBDIR)
endif
LIBS += $(BLAS_LIB)

# ARPACK
# use ARPACK_LIBDIR only if it is non-empty

ifneq ($(strip $(ARPACK_LIBDIR)),)
  LIBS += -L$(ARPACK_LIBDIR)
endif
LIBS += $(ARPACK_LIB)


ifeq ($(PRECISION),double)
  F90FLAGS += -DDOUBLE_PRECISION
endif


# netCDF library
# * use nc-config script if present (full path can be specified with 
# NETCDF_CONFIG environement variable)
# * if not use variables NETCDF_LIBDIR, NETCDF_INCDIR and NETCDF_LIB

NETCDF_VERSION := $(shell $(NETCDF_CONFIG) --version 2> /dev/null)

### check presense of nc-config script
ifeq ($(NETCDF_VERSION),)
  ifneq ($(strip $(NETCDF_INCDIR)),)
    F90FLAGS += -I$(NETCDF_INCDIR)
  endif

  # use NETCDF_LIBDIR only if it is non-empty
  ifneq ($(strip $(NETCDF_LIBDIR)),)
    LIBS += -L$(NETCDF_LIBDIR)
  endif
  LIBS += $(NETCDF_LIB)
else
  F90FLAGS += -I$(shell $(NETCDF_CONFIG) --includedir)
  LIBS += $(shell $(NETCDF_CONFIG) --flibs)
endif


F90FLAGS += $(EXTRA_F90FLAGS)
LDFLAGS += $(EXTRA_LDFLAGS)
