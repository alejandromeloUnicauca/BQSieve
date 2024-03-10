#Makefile
CC=gcc
CFLAGS= -lmpfr -lgmp -lm -std=c11 -ggdb3 -Wall 

# Lista de archivos fuente
SRCS=B_QSieve.c sieve.c polynomial.c factoring.c

# Lista de objetos generados a partir de los archivos fuente
OBJS=$(SRCS:.c=.o)

# Nombre del ejecutable final
TARGET=B_QSieve

# Regla para compilar el ejecutable
$(TARGET): $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS)

# Regla genérica para compilar archivos fuente a objetos
%.o: %.c
	$(CC) -c $< -o $@ $(CFLAGS)

.PHONY: clean

# Regla para limpiar archivos generados
clean:
	rm -f $(TARGET) $(OBJS)