CFLAGS += -Wall -pedantic

all: RIsearch RIsearch.dbg dsm.so

RIsearch: *.c *.h
	$(CC) $(CFLAGS) -O3 $^ -DRISVERSION=1 -lm -o $@

RIsearch.dbg: *.c *.h
	$(CC) $(CFLAGS) -g -O0 -DRISVERSION=1 -DDEBUG -DVERBOSE=2 $^ -lm -o $@

dsm.so: dsm.c
	$(CC) -fPIC -shared $(CFLAGS) dsm.c -o dsm.so

.PHONY: clean

clean:
	rm -f RIsearch  RIsearch.dbg *.o dsm.so
