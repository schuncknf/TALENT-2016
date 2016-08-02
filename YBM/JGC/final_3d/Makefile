CC	= g++
CFLAGS	= -O3 -g -Wall
LDFLAGS	=
EXEC	= exe
SRC	= spherical.c #hello.c
OBJ	= $(SRC:.c=.o)

all: $(EXEC)

$(EXEC): $(OBJ)
	$(CC) -o $@ $^ $(LDFLAGS)

%.o: %.c
	$(CC) -o $@ -c $< $(CFLAGS)

clean: 
	rm -rf *.o

mrproper: clean
	rm -rf $(EXEC)
