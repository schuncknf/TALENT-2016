CC	= gfortran
CFLAGS	= -fbounds-check -Wall

LDFLAGS	=
EXEC	= main

OBJS = main.o grid.o solver.o

all: $(EXEC)

$(EXEC): $(OBJS)
	$(CC) $(LDFLAGS) -o $@ $(OBJS)

grid.o: grid.f90

solver.o: solver.f90 grid.o

main.o: main.f90 grid.o solver.o

.SUFFIXES : .f90 .o

%.o : %.mod

.f90.o:
	$(CC) $(CFLAGS) -c $<

clean:
	rm -rf *.o

mrproper: clean
	rm -rf $(EXEC) *.mod
