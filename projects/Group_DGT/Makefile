SHELL= /bin/bash
cc=ifort
OPTS = -autodouble

exec=exe
SRC= radial.f90
OBJ=$(SRC:.f90=.o)

file : $(exec)

$(exec): $(OBJ)
	$(cc) -o $(exec) $(OBJ) 

radial.o:  
	$(cc) -o radial.o variables.f90 -c radial.f90
clean : 
	rm -rf *.o

mrproper : clean
	rm -rf $(exec)
