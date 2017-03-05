cf: variables_cf.o couette_main.o
	gfortran -o couette_main.exe variables_cf.o couette_main.o

couette_main.o: couette_main.f90 variables_cf.o
	gfortran -c couette_main.f90

variables_cf.o: variables_cf.f90
	gfortran -c variables_cf.f90

clean:
	rm -f *.o *.mod