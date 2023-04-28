# Targets: all, clean, tidy and help.

################################################################################

# Override with e.g. "make F90=ifort BUILD=opt OPENMP=no"

# Fortran compiler.  Options: gfortran, ifort, nagfor.
F90=gfortran

# Optimisation / debugging.  Options: opt, debug, none.
BUILD=opt

# Use OpenMP.  Options: yes, no.
OPENMP=yes

# Any other compiler flags the user wishes to specify.
# E.g., on HEC: FLAGS="-diag-disable 5140 -diag-disable 10010 -xSSE4.2
#                      -axCORE-AVX-I,CORE-AVX2"
FLAGS=

################################################################################

# Hopefully you will not need to edit the file below this point.

SHELL=/bin/sh

# Compiler options.
include makedefs.mk

# Check the compiler and build options.
$(if $(filter-out $(F90LIST),$(F90)), $(error F90 should be one of: $(F90LIST)))
$(if $(filter-out $(BUILDLIST),$(BUILD)), \
  $(error BUILD should be one of: $(BUILDLIST)))
$(if $(filter-out $(OPENMPLIST),$(OPENMP)), \
  $(error OPENMP should be one of: $(OPENMPLIST)))

# Sort out the compiler flags, etc.
FFLAGS=$(FFLAGS_$(F90)_$(BUILD)) $(FLAGS)
FFLAGS_OPT=$(FFLAGS_$(F90)_opt) $(FLAGS)
LDFLAGS_LIB=$(LDFLAGS_LIB_$(F90))
FFLAGS_OPENMP=$(FFLAGS_OPENMP_$(OPENMP)_$(F90))

# The binary files to produce.
BINFILE=fitter_mc

# Object files.
CANNED_F77_OBJFILES=ranlux.o
CANNED_F90_OBJFILES=machine_constants.o mrgrnk.o toms573.o
F90_OBJFILES=bootstrap.o fitter_mc.o utils.o
ALL_OBJFILES=$(CANNED_F77_OBJFILES) $(CANNED_F90_OBJFILES) $(F90_OBJFILES)

.SUFFIXES:

# Main target.
.PHONY: all
all: $(BINFILE)

# Rules for making object files.

# Canned routines (always using opt flags):
$(CANNED_F77_OBJFILES): %.o: %.f
	$(F90) $(FFLAGS_OPT) $(FFLAGS_OPENMP) -c -o $@ $<

$(CANNED_F90_OBJFILES): %.o: %.f90
	$(F90) $(FFLAGS) $(FFLAGS_OPENMP) -c -o $@ $<

# Other source files.
$(F90_OBJFILES): %.o: %.f90
	$(F90) $(FFLAGS) $(FFLAGS_OPENMP) -c -o $@ $<

# Module dependencies.
bootstrap.o: machine_constants.o mrgrnk.o toms573.o utils.o
fitter_mc.o: bootstrap.o utils.o
ranlux.o: machine_constants.o
toms573.o: machine_constants.o
utils.o: machine_constants.o ranlux.o

# Rules for linking object files.
$(BINFILE): $(ALL_OBJFILES)
	$(F90) $(FFLAGS) $(FFLAGS_OPENMP) -o $@ $(ALL_OBJFILES) $(LDFLAGS_LIB)

# Get rid of all results of the make process.
.PHONY: clean
clean:
	-rm -f $(BINFILE) *.o *.mod *.gcda *.gcno *~ histogram_*.dat

# Tidy up the object files, etc.
.PHONY: tidy
tidy:
	-rm -f *.o *.mod *.gcda *.gcno *~ histogram_*.dat

# Report list of targets.
.PHONY: help
help:
	@echo
	@echo "Targets: all, clean, tidy and help."
	@echo
	@echo "Type \"make F90=xxx BUILD=yyy OPENMP=zzz\", where"
	@echo "  xxx is one of: $(F90LIST),"
	@echo "  yyy is one of: $(BUILDLIST),"
	@echo "  zzz is one of: $(OPENMPLIST)."
	@echo
