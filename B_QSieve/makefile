#Makefile
# Variables comunes
CC=gcc
CFLAGS= -fopenmp -std=c11 -ggdb3 -Wall 
LDFLAGS = -lmpfr -lgmp -lm

# Lista de archivos fuente
SRCS = B_QSieve.c sieve.c polynomial.c factoring.c 

# Lista de objetos generados a partir de los archivos fuente
OBJS = $(SRCS:.c=.o)

# Nombre del ejecutable final
TARGET = B_QSieve
DEBUG_TARGET = B_QSieve_debug

# Regla por defecto para compilar en modo normal
all: $(TARGET)

# Regla para compilar el ejecutable en modo normal
$(TARGET): $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)

# Regla genérica para compilar archivos fuente a objetos en modo normal
%.o: %.c
	$(CC) -c $< -o $@ $(CFLAGS)

# Regla para compilar en modo depuración
debug: CFLAGS += -ggdb3
debug: $(DEBUG_TARGET)

$(DEBUG_TARGET): $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)

# Regla para limpiar archivos generados
clean:
	rm -f $(TARGET) $(DEBUG_TARGET) $(OBJS)

.PHONY: clean