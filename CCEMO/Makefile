# Makefile for compiling CCEMO_MPI source code
CC=g++
LD=g++
RM=rm -f
CFLAGS=-Wall -ansi -pedantic -g
OBJS:=$(patsubst %.cpp,%.o,$(wildcard *.cpp))
MAIN=main
all:$(MAIN)
$(MAIN):$(OBJS)
	$(LD) $(LDFLAGS) $(OBJS) -o $(MAIN) -lm
%.o: %.cpp CCEMO_MPI.h EMO_test_suite.h list.h rand.h
	$(CC) $(CFLAGS) -c $<
clean:
	$(RM) $(OBJS)

