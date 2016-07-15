# This Makefile compiles the required auxiliary functions for
# ARPACK as well as a Runge-Kutta integrator into a static
# library lauxiliary.a.
#
# This software is part of the PEPS Library.
# This software was programmed by Michael Lubasch when he was a
# PhD student in the Theory Division of the Max Planck Institute
# of Quantum Optics where his work was supervised by
# Prof. Dr. J. Ignacio Cirac and Dr. Mari-Carmen Ba\~{n}uls.
# This software is distributed under the
# PEPS Library License v1.0 (see accompanying file LICENSE.txt).

# The Fortran compiler and optimization:
FC = gfortran
FFLAGS = -O
#FFLAGS = -O -fopenmp

# The archiver:
AR = ar
ARFLAGS = rv
RANLIB = ranlib
LIB = lauxiliary.a

# Rule to obtain .o files from .f files:
.f.o:
	$(FC) $(FFLAGS) -c $<

# The object files:
OBJ = ssyevar.o ssyevar2.o dsyevar.o dsyevar2.o zrk4.o

all: lauxiliary

lauxiliary: $(OBJ)
	$(AR) $(ARFLAGS) $(LIB) $(OBJ)
	$(RANLIB) $(LIB)

# Clean:

clean:
	rm -f *.o
