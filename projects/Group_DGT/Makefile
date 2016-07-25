SHELL= /bin/bash
cc=ifort
OPTS = -autodouble

exec=test
SRC= leftright1.f90
OBJ=$(SRC:.f90=.o)

all : $(exec)

$(exec): $(OBJ)
	$(cc) -o $(exec) $(OBJ) 

leftright1.o:  
	$(cc) -o leftright1.o  -c leftright1.f90
clean : 
	rm -rf *.o

mrproper : clean
	rm -rf $(exec)
