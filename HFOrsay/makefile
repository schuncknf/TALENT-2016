blddir=./bld
$(shell mkdir $(blddir) 2> /dev/null)

# For Serial Version
#FC = ifort -132 -O2 -fno-alias -vec_report=0 -fomit-frame-pointer -g -traceback
FC1 = ifort -132 -O2 -fno-alias -vec_report=0 -fomit-frame-pointer
FC = ifort -132 -O2 -debug all -traceback -check all -warn usage -warn errors -warn general
#FC = ifort -debug all -traceback -check all -warn usage -warn errors -warn general
LK = ifort -O2 
#FC = gfortran -O3
#LK = gfortran -O3
#-cpp
#  FFLAGS=-O5 -align dcommons
# FFLAGS=-O4 -align dcommons
# FFLAGS=-O5 -align dcommons -tune host

OBJ =  constants.o maths.o laguerre_polynomial.o lag_pol.o ho.o pot.o matrices.o \
	 hfsolver.o gausslag.o tbme.o hforsay.o\
	
       
BLDOBJ = $(foreach obj,$(OBJ),$(blddir)/$(obj))

$(blddir)/%.o: %.f
	$(FC) -c -o $@ $< -L Nag/ -lnag

$(blddir)/%.o:%.f90
	$(FC1) -c -o $@ $< -llapack -lblas

$(blddir)/%.o:%.c
	$(CC) -c -o $@ $< -llapack -lblas

all: run

run: $(BLDOBJ)
	$(LK) -O2 -o run $(BLDOBJ) -llapack -lblas

clean:
	rm -rf $(blddir)/* run
