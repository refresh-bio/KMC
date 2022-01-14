CC = g++
CFLAGS = -std=c++20 -O3 -Wall -pedantic-errors

all: bin/counter

main.o: main.cpp
	$(CC) $(CFLAGS) -c -o $@ main.cpp

bin/counter: main.o
	mkdir -p bin
	$(CC) $(CFLAGS) -o $@ $^ -lz 


clean:
	rm -f *.o
	rm -rf bin	