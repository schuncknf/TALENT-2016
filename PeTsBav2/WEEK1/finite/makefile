COMP = gfortran 

 
CFLAGS = -c  -g -C -fbounds-check -fcheck=all -fbacktrace -pedantic -Wextra -Wall
#CFLAGS = -c  -O3 -funroll-loops -msse4a -mtune=native -march=native -ffast-math -floop-parallelize-all -ftree-vectorize

LD_FLAGS = -O3 


PROGNAME=main


#   Fortran objects
objs= global.o main.o	\
	Simpson.o \
#
.SUFFIXES: .f90 .f	.o
$(PROGNAME):	$(objs)
	$(COMP)	-o $(PROGNAME) $(objs)   $(LD_FLAGS)
.f90.o:
	$(COMP) $(CFLAGS)	-o $(@) $<;
.f.o:
	$(COMP) $(CFLAGS)	-o $(@) $<;
modules:	
	$(COMP) $(CFLAGS) $(mods)

clear:
	rm -f *.o *.bak *.lst *.mod $(PROGNAME);
clean:
	rm -f *.o *.bak *.lst *.mod $(PROGNAME);
