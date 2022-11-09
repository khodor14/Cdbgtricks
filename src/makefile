CC=g++
CFLAGS=  -Wall  -Wextra -Ofast -std=c++11  -pthread -fmax-errors=3  -pipe
LDFLAGS=-pthread

ifeq ($(gprof),1)
CFLAGS=-std=c++11 -pg -O3  -march=native
LDFLAGS=-pg
endif

ifeq ($(valgrind),1)
CFLAGS=-std=c++0x -O3 -g
LDFLAGS=-g
endif

EXEC=main

all: $(EXEC)

main:	main.o	ParseGFA.o	CommonUtils.o	unitigmodification.o index_kmer.o
	$(CC) -o $@ $^ $(LDFLAGS)

main.o: 	main.cpp	ParseGFA.h	CommonUtils.h	unitigmodification.h
	$(CC) -o $@ -c $< $(CFLAGS)

ParseGFA.o: 	ParseGFA.cpp	ParseGFA.h
	$(CC) -o $@ -c $< $(CFLAGS)

index_kmer.o:	index_kmer.cpp	index_kmer.h	CommonUtils.h
	$(CC) -o $@ -c $< $(CFLAGS)

unitigmodification.o:	unitigmodification.cpp unitigmodification.h ParseGFA.h index_kmer.h
	$(CC) -o $@ -c $< $(CFLAGS)

CommonUtils.o:	CommonUtils.cpp CommonUtils.h
	$(CC) -o $@ -c $< $(CFLAGS)

clean:
	rm -r *.o
	rm -rf $(EXEC)


rebuild: clean $(EXEC)