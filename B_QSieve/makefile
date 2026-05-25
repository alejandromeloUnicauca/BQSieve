#Makefile
# Variables comunes
CC=gcc
# -MMD -MP: gcc genera archivos .d con las dependencias de headers de cada .o,
# de modo que editar un .h recompila todos los .c que lo incluyen.
CFLAGS= -O2 -mavx2 -fopenmp -std=c11 -Wall -MMD -MP
LDFLAGS = -lmpfr -lgmp -lm

# Lista de archivos fuente
SRCS = B_QSieve.c sieve.c polynomial.c factoring.c linalg.c

# Lista de objetos generados a partir de los archivos fuente
OBJS = $(SRCS:.c=.o)

# Nombre del ejecutable final
TARGET = B_QSieve
DEBUG_TARGET = B_QSieve_debug

# Solver GF(2) nativo (standalone, útil para depuración)
SOLVE_SRCS = solve_matrix.c linalg.c
SOLVE_OBJS = solve_matrix.o linalg.o
SOLVE_TARGET = solve_matrix

# Regla por defecto
all: $(TARGET) $(SOLVE_TARGET)

# Regla para compilar el ejecutable en modo normal
$(TARGET): $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)

$(SOLVE_TARGET): $(SOLVE_OBJS)
	$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)

solve_matrix.o: solve_matrix.c linalg.h
	$(CC) -c $< -o $@ $(CFLAGS)

linalg.o: linalg.c linalg.h
	$(CC) -c $< -o $@ $(CFLAGS)

# Regla genérica para compilar archivos fuente a objetos en modo normal
%.o: %.c
	$(CC) -c $< -o $@ $(CFLAGS)

# Regla para compilar en modo depuración (sin optimización, con símbolos)
debug: CFLAGS += -O0 -ggdb3
debug: $(DEBUG_TARGET) $(SOLVE_TARGET)

$(DEBUG_TARGET): $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)

# Incluir los archivos de dependencias generados por -MMD (si existen)
-include $(OBJS:.o=.d) $(SOLVE_OBJS:.o=.d)

# Regla para limpiar archivos generados
clean:
	rm -f $(TARGET) $(DEBUG_TARGET) $(SOLVE_TARGET) \
	      $(OBJS) $(SOLVE_OBJS) \
	      $(OBJS:.o=.d) $(SOLVE_OBJS:.o=.d)

.PHONY: clean
