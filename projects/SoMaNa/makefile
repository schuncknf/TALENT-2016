CC = gcc
CFLAGS = -W -Wall -I/usr/userfs/t/tmpu3150/include/include -L/usr/userfs/t/tmpu3150/include/lib -o
LDFLAGS = -lgsl -lgslcblas -lm
PROGRAM = main
SOURCE = gauss-laguerre.c harmon.c tbmeoutput.c main.c read_file.c
OBJS = $(SOURCE: .c=.o)

all: $(PROGRAM)

$(PROGRAM): $(OBJS)
	$(CC) $(CFLAGS) $@ $^ $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) $@ $^ $(LDFLAGS)

clean: 
	rm -r *.o

rebuild: clean
	rm project
