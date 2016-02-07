# Makefile for compiling CMODE source code
CC=g++
LD=g++
RM=rm -f
CFLAGS=-Wall -ansi -pedantic -g
OBJS=EMO_test_suite.o cec09.o common.o \
	 individual.o ran.o rcmode.o main.o
MAIN=main
all:$(MAIN)
$(MAIN):$(OBJS)
	$(LD) $(LDFLAGS) $(OBJS) -o $(MAIN) -lm
	
EMO_test_suite.o: EMO_test_suite.cpp \
				  EMO_test_suite.h
	$(CC) $(CFLAGS) -c $<
cec09.o: cec09.cpp cec09.h
	$(CC) $(CFLAGS) -c $<
common.o: common.h
	$(CC) $(CFLAGS) -c $<
individual.o: individual.h main.cpp \
			  EMO_test_suite.h
	$(CC) $(CFLAGS) -c $<
ran.o: ran.h main.cpp
	$(CC) $(CFLAGS) -c $<
rcmode.o: rcmode.h main.cpp
	$(CC) $(CFLAGS) -c $<
main.o: main.cpp EMO_test_suite.h \
		ran.h common.h individual.h \
		rcmode.h
	$(CC) $(CFLAGS) -c $<
clean:
	$(RM) $(OBJS)

