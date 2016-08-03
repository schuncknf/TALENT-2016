blddir=./bld
$(shell mkdir $(blddir) 2> /dev/null)

#FC1 = ifort -132 -O2 -fno-alias -fomit-frame-pointer -g -traceback -warn errors
#LK = ifort -O2 

FC1 = gfortran -O2
LK = gfortran -O2

OBJ =  constants.o maths.o base.o laguerre_polynomial.o lag_pol.o ho.o pot.o matrices.o \
	 bcs.o solver.o printer.o gausslag.o spatial.o tbme.o hforsay.o\
	
       
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
