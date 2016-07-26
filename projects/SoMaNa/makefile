CC = gcc
CFLAGS = -W -Wall -I/usr/userfs/t/tmpu3150/include/include -L/usr/userfs/t/tmpu3150/include/lib -o
LDFLAGS = -lgsl -lgslcblas -lm
PROGRAM = gauss-laguerre harmon test tbmeoutput hermv read_file
SOURCE = gauss-laguerre.c harmon.c test.c tbmeoutput.c hermv.c read_file.c
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
