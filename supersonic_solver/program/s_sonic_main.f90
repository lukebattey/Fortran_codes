PROGRAM s_sonic_main

USE variables_ss
USE boundary_conds
USE get_misc

    IMPLICIT NONE

!======== READ INPUT FILE (all global variables) ============
OPEN(16,FILE = 'inputfile.dat', FORM = 'FORMATTED')
READ(16,*) infile 
READ(16,*) outfile
READ(16,*) ui
READ(16,*) rhoi
READ(16,*) gama
READ(16,*) CFL
READ(16,*) order
READ(16,*) fluxlim
CLOSE(16) 

!======== READ IN GRID, GET METRICS, AND JACOBIANS ==========
CALL grid_read_metrics

!======== SET INITIAL CONDITIONS AND BCs ====================
ALLOCATE(Util(imax,jmax),Vtil(imax,jmax), &
         del2(imax),del3(imax),deln1(imax),deln2(imax), &
         h0(imax,jmax))

ALLOCATE(URG(imax-2,jmax-1,4),ULG(imax-2,jmax-1,4), &
         URF(imax-1,jmax-2,4),ULF(imax-1,jmax-2,4), &
         Fpr(imax-1,jmax,4),Gpr(imax-1,jmax,4)) 

CALL init_cond
CALL get_dt !<---- Requires init_cond beforehand..

CALL stag_enthalpy

CALL outflow_BC
CALL lower_BC
CALL upper_BC

IF (order == 1) THEN
    CALL extr_Ustate_o1
ELSE IF (order == 2) THEN
    WRITE(6,*) 'UNDER CONSTRUCTION QUITTING'
    STOP
ELSE 
    WRITE(6,*) '1st or 2nd ORDER ONLY PLEASE!'
    STOP
END IF



CALL write_sol

END PROGRAM s_sonic_main