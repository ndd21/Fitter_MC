# Definitions of compiler flags, etc.

# List of allowed values of F90.
F90LIST=gfortran ifort nagfor
# List of allowed values of BUILD.
BUILDLIST=debug none opt
# List of allowed values of OPENMP.
OPENMPLIST=yes no

# Compiler flags, etc., for gfortran (GCC).
MPIF90_gfortran=mpif90
FFLAGS_gfortran_debug=-Og -std=f2003 -pedantic -Wall -Wextra\
 -Wimplicit-interface -Wimplicit-procedure -Wunderflow -fimplicit-none\
 -fbacktrace -fcheck=all -g
FFLAGS_gfortran_none=
FFLAGS_gfortran_opt=-Ofast -fpeel-loops -fprotect-parens -march=native\
 -mtune=native
FFLAGS_OPENMP_yes_gfortran=-fopenmp
FFLAGS_OPENMP_no_gfortran=
LDFLAGS_LIB_gfortran= #-lopenblas

# Compiler flags, etc., for ifort.
MPIF90_ifort=mpifort
FFLAGS_ifort_debug=-O0 -check all -check noarg_temp_created -fp-stack-check\
 -g -implicitnone -std03 -traceback -warn all,nounused -debug all -ftrapuv\
 -assume buffered_io -Vaxlib
FFLAGS_ifort_none=-Vaxlib
FFLAGS_ifort_opt=-O3 -ip -no-prec-div -no-prec-sqrt -funroll-loops\
 -no-fp-port -complex-limited-range -assume protect_parens\
 -assume buffered_io -Vaxlib
FFLAGS_OPENMP_yes_ifort=-qopenmp
FFLAGS_OPENMP_no_ifort=
LDFLAGS_LIB_ifort= #-mkl=sequential

# Compiler flags, etc., for NAG.  Only debugging flags are available.
MPIF90_nagfor=nagfor
FFLAGS_nagfor_debug=-C=all -colour -f2003 -gline -info -mtrace=off -nan -O0\
 -strict95 -u -v
FFLAGS_nagfor_none=
FFLAGS_nagfor_opt=$(FFLAGS_nagfor_debug)
FFLAGS_OPENMP_yes_nagfor=-openmp
FFLAGS_OPENMP_no_nagfor=
LDFLAGS_LIB_nagfor= #-lblas -llapack
