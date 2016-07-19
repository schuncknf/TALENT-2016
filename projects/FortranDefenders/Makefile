CC	= gfortran
CFLAGS	= -g -fbounds-check -Wall
LDFLAGS	=
EXEC	= main

OBJS = main.o grid.o solver.o

SHELL =/bin/sh

all: $(EXEC)

$(EXEC): $(OBJS)
	$(CC) $(LDFLAGS) -o $@ $(OBJS)

grid.o: grid.f95

solver.o: solver.f95 grid.o

main.o: main.f95 grid.o solver.o

%.o : %.mod

.SUFFIXES : .f95 .f .o

.f95.o:
	$(CC) $(CFLAGS) -c $<

clean:
	rm -rf *.o

mrproper: clean
	rm -rf $(EXEC) *.mod
