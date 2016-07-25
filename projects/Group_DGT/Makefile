#SHELL= /bin/bash
cc=ifort
OPTS = -autodouble

exec=test
SRC= radial.f90
OBJ=$(SRC:.f90=.o)

all : $(exec)

$(exec): $(OBJ)
	$(cc) -o $(exec) $(OBJ) 

radial.o:  
	$(cc) -o radial.o  -c radial.f90
clean : 
	rm -rf *.o

mrproper : clean
	rm -rf $(exec)
