CC	= gfortran
CFLAGS	= 
LDFLAGS	=
EXEC	= main

OBJS = main.o grid.o

all: $(EXEC)

$(EXEC): $(OBJS)
	$(CC) $(LDFLAGS) -o $@ $(OBJS)

grid.o: grid.f95

main.o: main.f95 grid.o

%.o : %.mod

.SUFFIXES : .f95 .f .o

.f95.o:
	$(CC) $(CFLAGS) -c $<

clean: 
	rm -rf *.o

mrproper: clean
	rm -rf $(EXEC) *.mod
