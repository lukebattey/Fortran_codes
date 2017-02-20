FC    = ifort 	 
FLAGS = -cpp -O3
EXEC  = GT_induced_v1.0

in:
	make FLAGS="$(FLAGS) -DLOOSE -DPARALLEL" FC="/usr/lib64/openmpi/bin/mpif90" EXEC="$(EXEC)" GT_induced

test:
	make FLAGS="$(FLAGS) -DLOOSE -DPARALLEL" FC="/usr/lib64/openmpi/bin/mpif90" EXEC="$(EXEC)" GT_induced_test

GT_induced_test : clean modules.o readin.o main.o
	$(FC) -o $(EXEC) \
            	  modules.o readin.o main.o

GT_induced : modules.o readin.o main.o induced_vel_at_P.o
	$(FC) -o $(EXEC) \
            modules.o readin.o main.o induced_vel_at_P.o

modules.o: modules.f90
	$(FC) $(FLAGS)  -c modules.f90

readin.o: readin.f90
	$(FC) $(FLAGS)  -c readin.f90

main.o: main.f90
	$(FC) $(FLAGS)  -c main.f90

induced_vel_at_P.o: induced_vel_at_P.f90
	$(FC) $(FLAGS)  -c induced_vel_at_P.f90
	
clean:
	rm *.o *.mod 
