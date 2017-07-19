CC=gcc
#CCFLAGS=-Wall -std=gnu11    -O
CCFLAGS=-Wall -std=gnu11    -g

HEADERS = infdecomp.h infdecomp_p.h
SOURCES = infdecomp.c
MAINS   =
OBJECTS = $(SOURCES:.c=.o)

all: infdecomp.o

#blah: $(OBJECTS) blah.c
#	$(CC) $(CCFLAGS) blah.c $(OBJECTS) -lm -o blah

%.o: %.c
	$(CC) $(CCFLAGS) -c $<

clean:
	rm $(OBJECTS) # $(PSLISTINGS) $(LIBRARY)

depend:
	makedepend -o.o -f Makefile.depend -- $(CCCFLAGS) -- $(SOURCES)

include Makefile.depend
