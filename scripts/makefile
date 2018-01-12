#
# Makefile 
#
# Usage:  make all
#  will compile code and place the object files into 
# ./obj

# modifications to have different options on linux and MAC
UNAME := $(shell uname)

ifeq ($(UNAME), Linux)
  # Set flags to point the compiler to the include files
  FF           = gfortran
  CC           = cc 
  FFLAGS       = -O2 -Wall #O2
  CCOMFLG      = -O3 #O3 
  LIB          = -L/usr/libx -lnetcdff -lnetcdf
  INCLUDE      = -I/usr/include 
  EXE          = vectri.exe

endif
ifeq ($(UNAME), Darwin)
  FF           = gfortran-mp-4.8
  CC           = cc 
  FFLAGS       = -O2 -Wall #O2
  CCOMFLG      = -O3 #O3 
  LIB          = -L/opt/local/lib -lnetcdff -lnetcdf
  INCLUDE      = -I/opt/local/include 
  EXE          = vectri.exe
endif

.SUFFIXES: .f90 .c

.f90.o :
	$(FF) $(INCLUDE) -c $(FFLAGS)  $*.f90 -o $*.o $(LIB)
.c.o:
	$(CC) -c $(CCOMFLG) $*.c   -o $*.o

FORT  = mo_climate.f90 mo_constants.f90 mo_control.f90  mo_vectri.f90 mo_namelist.f90 mo_ncdf_tools.f90 open_input.f90 open_output.f90 initialize.f90 climate_pop_trends.f90 getdata.f90 writedata.f90 setup.f90 setdown.f90 advection.f90 transmission.f90 vectri.f90

CPP =\

OBJECT=$(FORT:.f90=.o) $(CPP:.c=.o)


vectri:$(OBJECT)
	$(FF) $(INLCUDE) $(OBJECT) -o $(EXE)  $(LIB)


all:	vectri.exe	

clean:
	rm -f *.o *.mod $(EXE)

