FF = gfortran
EXEC = main
DEPS = hello.h

all: $(EXEC) 

$(EXEC): main.o variables.o potentials.o calcwaves.o
	$(FF) -fno-whole-file main.f90 variables.f90 potentials.f90 calcwaves.f90 -o main

variables.o: variables.f90 
	$(FF) -o variables.o -c variables.f90 $(CFLAGS)

potentials.o: potentials.f90 variables.f90 
	$(FF) -o variables.o -c variables.f90 potentials.f90 $(CFLAGS)

calcwaves.o: calcwaves.f90 variables.f90 
	$(FF) -o calcwaves.o -c variables.f90 calcwaves.f90 $(CFLAGS)

main.o : main.f90
	$(FF) -o main.o -c main.f90 $(CFLAGS)

clean:
	rm -rf *.o

mrproper: clean
	rm -rf $(EXEC)
