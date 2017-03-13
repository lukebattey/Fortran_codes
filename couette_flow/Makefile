cf: variables_cf.o CMA_heatEq_sub.o EXPL_heatEq_sub.o couette_main.o
	gfortran -o couette_main.exe variables_cf.o CMA_heatEq_sub.o EXPL_heatEq_sub.o couette_main.o

wi: variables_cf.o write_initial.o
	gfortran -o write_initial.exe variables_cf.o write_initial.o

couette_main.o: couette_main.f90 variables_cf.o
	gfortran -c couette_main.f90

write_initial.o: write_initial.f90 variables_cf.o
	gfortran -c write_initial.f90

CMA_heatEq_sub.o: CMA_heatEq_sub.f90 variables_cf.o
	gfortran -c CMA_heatEq_sub.f90

EXPL_heatEq_sub.o: EXPL_heatEq_sub.f90 variables_cf.o
	gfortran -c EXPL_heatEq_sub.f90

variables_cf.o: variables_cf.f90
	gfortran -c variables_cf.f90

clean:
	rm -f *.o *.mod