#Makefile
CC=gcc
CFLAGS= -lmpfr -lgmp -lm -std=c11 -ggdb3 -Wall 

prsa: prsa.c
	$(CC) $(CFLAGS) -o prsa prsa.c split.c
prsaBlocks: prsaBlocks.c
	$(CC) $(CFLAGS) -o prsaBlocks prsaBlocks.c split.c
