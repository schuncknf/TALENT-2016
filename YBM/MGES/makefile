CC	= g++
CFLAGS	= -W -Wall -ansi -pedantic
LDFLAGS	= -lstdc++
EXEC	= squareWell
SRC	= squareWell.c integration.c
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