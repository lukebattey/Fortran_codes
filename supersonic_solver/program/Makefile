ss: variables_ss.o grid_read_metrics.o init_L_BC.o s_sonic_main.o
	gfortran -o s_sonic_main.exe variables_ss.o grid_read_metrics.o \
	init_L_BC.o s_sonic_main.o 

s_sonic_main.o: s_sonic_main.f90 variables_ss.o
	gfortran -c s_sonic_main.f90

init_L_BC.o: init_L_BC.f90 variables_ss.o
	gfortran -c init_L_BC.f90

grid_read_metrics.o: grid_read_metrics.f90 variables_ss.o
	gfortran -c grid_read_metrics.f90

variables_ss.o: variables_ss.f90
	gfortran -c variables_ss.f90

clean:
	rm -f *.o *.mod