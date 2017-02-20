alg_gridgen: alg_gridgen.o main.o
	gfortran -o alg_gridgen.exe main.o alg_gridgen.o
main.o: main.f90 alg_gridgen.o
	gfortran -c main.f90
alg_gridgen.o: alg_gridgen.f90
	gfortran -c alg_gridgen.f90
clean:
	rm -f *.o *.mod alg_gridgen.exe
