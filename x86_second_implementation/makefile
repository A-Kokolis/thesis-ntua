CSHELL=sh

# --- MACROS
# define program name
MAIN= infoli_rst

# define C source files
SRC = infoli_mpi.c

# define C header files
HDRS= infoli_mpi.h

#define C object files
OBJ= infoli_mpi.o

# --- TARGETS
all:	infoli main

main:	${OBJ}
	@echo #
	@echo "-- LINKING PROGRAM --"
	mpicc ${OBJ} -o ${MAIN} -O3 -lm

infoli: ${SRC} ${HDRS}
	@echo #
	@echo "-- BUILDING MAIN INFOLI --"
	mpicc -w -c -O3 -lm ${SRC}

clean:
	@echo #
	@echo "-- CLEANING PROJECT FILES --"
	rm -f *.o ${MAIN}

