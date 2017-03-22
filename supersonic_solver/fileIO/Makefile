VPATH = ./src:

.SUFFIXES : .F90 .f .f90 .c .o

# Compiler types
FC = ifort
CC = icc
CXX = icpc

# Optimization level
OPT    = -O2

# Relevant flags
FFLAGS = ${OPT} 
CFLAGS = ${OPT} 
CXXFLAGS = ${OPT} 

# Name of the program
PROG  = fileIO.x

# List of object files
OBJ = tecio.o\
      fileIO.o

.f90.o: 
	${FC} -c -o $@ ${FFLAGS} $<

.F90.o: 
	${FC} -c -o $@ ${FFLAGS} $<

.c.o: 
	${CC} -c -o $@ ${FFLAGS} $<

all:${PROG}

$(PROG): ${OBJ}
	$(FC) ${FFLAGS} $(OBJ) -o ${PROG}
	mv ${PROG} ./bin/

clean:
	/bin/rm -f *.mod *.o *.x *~ ./src/*.o ./src/*~

uninstall: 
	rm -f bin/${PROG} 

renew:
	touch *