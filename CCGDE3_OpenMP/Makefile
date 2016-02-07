# Makefile for compiling CCGDE3 source code
CC=mpicxx
LD=mpicxx
RM=rm -f
CFLAGS=-Wall -pedantic -O3 -std=c++0x  -ggdb -DDEBUG -g -pg
LDFLAGS=-Wall -pedantic -O3 -std=c++0x  -ggdb -DDEBUG -g -pg
#CFLAGS=-fopenmp -Wall -ansi -pedantic -ggdb -DDEBUG -pg
OBJS:=$(patsubst %.cpp,%.o,$(wildcard *.cpp))
MAIN=CCGDE3
all:$(MAIN)
$(MAIN):$(OBJS)
	$(LD) $(LDFLAGS) -o $(MAIN) $(OBJS)
%.o: %.cpp CCGDE3.h EMO_test_suite.h rand.h list.h ZDT.h
	$(CC) $(CFLAGS) -c $<
clean:
	$(RM) $(OBJS)

