#Makefile
CC=gcc
CFLAGS= -lmpfr -lgmp -lm -fopenmp -std=c11 -ggdb3 -Wall 

B_QSieve: B_QSieve.c
	$(CC) -o B_QSieve B_QSieve.c sieve.c polynomial.c factoring.c $(CFLAGS) 
