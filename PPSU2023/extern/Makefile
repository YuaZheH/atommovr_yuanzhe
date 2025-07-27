include ../Makefile.inc

all:  cheap.o matching.o

cheap.o: cheap.c
matching.o: matching.c

%.o: %.c
	$(CC) -c $< $(CFLAGS)

clean:
	rm -f *.o

