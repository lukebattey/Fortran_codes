ss: variables_ss.o get_misc.o grid_read_metrics.o init_cond.o \
	boundary_conds.o extr_Ustate_o1.o upwind_U_LR_o2.o fluxlim_U_LR_o2.o \
	ViscFlux.o AUSMPWpF.o AUSMPWpG.o update_state.o s_sonic_main.o
	gfortran -o s_sonic_main.exe variables_ss.o get_misc.o \
	grid_read_metrics.o init_cond.o boundary_conds.o extr_Ustate_o1.o \
	upwind_U_LR_o2.o fluxlim_U_LR_o2.o ViscFlux.o AUSMPWpF.o AUSMPWpG.o \
	update_state.o s_sonic_main.o 

s_sonic_main.o: s_sonic_main.f90 variables_ss.o
	gfortran -c s_sonic_main.f90

update_state.o: update_state.f90 variables_ss.o
	gfortran -c update_state.f90

AUSMPWpG.o: AUSMPWpG.f90 variables_ss.o
	gfortran -c AUSMPWpG.f90

AUSMPWpF.o: AUSMPWpF.f90 variables_ss.o
	gfortran -c AUSMPWpF.f90

ViscFlux.o: ViscFlux.f90 variables_ss.o get_misc.o
	gfortran -c ViscFlux.f90

fluxlim_U_LR_o2.o: fluxlim_U_LR_o2.f90 
	gfortran -c fluxlim_U_LR_o2.f90

upwind_U_LR_o2.o: upwind_U_LR_o2.f90 variables_ss.o
	gfortran -c upwind_U_LR_o2.f90

extr_Ustate_o1.o: extr_Ustate_o1.f90 variables_ss.o
	gfortran -c extr_Ustate_o1.f90

boundary_conds.o: boundary_conds.f90 variables_ss.o
	gfortran -c boundary_conds.f90

init_cond.o: init_cond.f90 variables_ss.o
	gfortran -c init_cond.f90

grid_read_metrics.o: grid_read_metrics.f90 variables_ss.o
	gfortran -c grid_read_metrics.f90

get_misc.o: get_misc.f90 variables_ss.o
	gfortran -c get_misc.f90

variables_ss.o: variables_ss.f90
	gfortran -c variables_ss.f90

clean:
	rm -f *.o *.mod