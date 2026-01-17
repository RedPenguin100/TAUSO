CFLAGS += -fopenmp -Wall -I../libdivsufsort-2.0.1/build/include -DRISEARCH2  -DRISVERSION=2
CLIBS += -lz -lpcre -lm

all: ../bin/risearch2.x ../bin/risearch2.dbg.x

../bin/risearch2.x: *.c *.h
	$(CC) $(CFLAGS) -O3 $^ ../libdivsufsort-2.0.1/build/lib/libdivsufsort64.a $(CLIBS) -o $@

../bin/risearch2.dbg.x: *.c *.h
	$(CC) $(CFLAGS) -ggdb -g3 -O0 -DDEBUG -DVERBOSE=2 $^ ../libdivsufsort-2.0.1/build/lib/libdivsufsort64.a $(CLIBS) -o $@

.PHONY: clean

clean:
	rm -f ../bin/risearch2.x ../bin/risearch2.dbg.x conv.x
