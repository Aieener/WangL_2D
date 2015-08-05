# declare the variable
# CC = g++ -std=gnu++0x
CC = g++

CFLAGS = -c -Wall -pedantic -Ofast -march=native -std=gnu++0x

PROG1 = MC

all: runit


runit: $(PROG1).o cells.o square.o hardrods.o histogram.o
		$(CC) $(PROG1).o cells.o square.o hardrods.o histogram.o -o runit

$(PROG1).o: $(PROG1).cpp
		$(CC) $(CFLAGS) $(PROG1).cpp

cells.o: cells.cpp
		$(CC) $(CFLAGS) cells.cpp

square.o: square.cpp
		$(CC) $(CFLAGS) square.cpp

hardrods.o: hardrods.cpp
		$(CC) $(CFLAGS) hardrods.cpp

histogram.o: histogram.cpp
		$(CC) $(CFLAGS) histogram.cpp

clean:
		rm -rf *o run
