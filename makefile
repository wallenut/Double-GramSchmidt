CC=gcc
CFLAGS= -std=c99 -pedantic -Wall -g3

#####
# Instructions to make hw1
#####

hw1: hw1.o
	${CC} ${CFLAGS} -o hw1 hw1.o
