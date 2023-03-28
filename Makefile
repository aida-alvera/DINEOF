#==========================================================================
#
# Please do not modify this file
# Make your changes to config.mk instead
#
#==========================================================================

include config.mk

# default settings
# If you need to adapt the variable, then make the corresponding changes in config.mk

OS ?= Linux
FORT ?= gfortran
DEBUG ?= on
OPENMP ?= on

# DINEOF version
VERSION ?= 4.0


# Example files are provided in the format "big endian". 
# Do not change the FORMAT variable, if you want to run the test case "SmallExample".
FORMAT ?= big_endian
PRECISION ?= double

#==========================================================================
# Include compiler specific options
#==========================================================================

include Compilers/$(OS)-$(strip $(FORT)).mk
include Compilers/libs.mk

F90FLAGS += $(INCLUDES)

#==========================================================================
# All Source files
#==========================================================================

SOURCE = ReadMatrix.F90 ufileformat.F90 initfile.F90 stat.F90 norm.F90 dineof_utils.F90 \
  smeanToZero.F smeanByRow.F svariExp.F ssvd_lancz.F dineof.F90
  

#==========================================================================
# Object files for "dineof"
#==========================================================================


OBJECTS =  ReadMatrix.o ufileformat.o initfile.o stat.o norm.o dineof_utils.o \
  smeanToZero.o smeanByRow.o svariExp.o ssvd_lancz.o dineof.o 
  

#==========================================================================
# Object files for "crossval"
#==========================================================================


OBJECTSC =  ReadMatrix.o ufileformat.o initfile.o stat.o norm.o crossval.o

#==========================================================================
# Declare .F90 a valid suffix
#==========================================================================

.SUFFIXES: .F90 .F

#==========================================================================
# How to compile F77 and F90 programs?
#==========================================================================

%.o: %.F
	$(F90C) $(F90FLAGS) -c $<

%.o: %.F90
	$(F90C) $(F90FLAGS) -c $<

all:  dineof

dineof: $(OBJECTS)
	$(F90C) $(LDFLAGS) -o $@ $(OBJECTS) $(LIBS)

#crossval: $(OBJECTSC)
#	$(F90C) $(LDFLAGS) -o $@ $(OBJECTSC) $(LIBS)

clean:
	rm -f $(OBJECTS) $(OBJECTSC) *.mod 

distclean:
	rm -f  $(OBJECTS) $(OBJECTSC) *.mod dineof #crossval

tarfile:
	if [ -e dineof-$(VERSION).tar.gz ]; then  rm -i dineof-$(VERSION).tar.gz; fi
	tar -C ../ --exclude-vcs -zcvf dineof-$(VERSION).tar.gz dineof

print:
	echo $(OBJECTS) 

#==========================================================================
# Dependencies
#==========================================================================

ufileformat.o: ppdef.h
initfile.o: ppdef.h
#crossval.o: ufileformat.o initfile.o
ReadMatrix.o: ufileformat.o initfile.o
writeMatrix.o: ufileformat.o
writeSVD.o:  ufileformat.o
dineof.o: ReadMatrix.h ufileformat.o initfile.o  dineof_utils.o
ssvd_lancz.f: includes/debug.h 
stat.o: ufileformat.o
siterativeEof.o:  ufileformat.o
ssvd_lancz.o: ufileformat.o dineof_utils.o
dineof_utils.o: ufileformat.o 
