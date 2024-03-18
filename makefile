# Compiles main.f

# compiler flags
FC = gfortran
LDR = gfortran

# debugging options
FFLAGS  = -ffixed-line-length-136 -fdefault-real-8 -Wall -fno-automatic -Ofast -g -fbounds-check -fbacktrace -finit-real=nan

# MacOSX libraries
GLIB = -L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib

# source files
SRCS = element_table_mod flyv leibniz find_params read_xdat factorial polint polcof calc_forces find_forces metric gdot calc_1dforcedens calc_frdf main
OBJS = $(SRCS:=.o)

# executable
MAIN = main

# compile project
all : $(MAIN)
	@echo Model compiled

$(MAIN) : $(OBJS)
	$(FC) $(FFLAGS) $(GLIB) -o $(MAIN) $(OBJS)

.SUFFIXES : .o .f .f90

.f90.o : 
	$(FC) $(FFLAGS) $(GLIB) -c $<

clean:
	$(RM) *.o *.mod $(MAIN) 
