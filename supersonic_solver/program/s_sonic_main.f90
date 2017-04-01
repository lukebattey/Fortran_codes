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
         h0(imax,jmax))

ALLOCATE(del2(imax),del3(imax),deln1(imax),deln2(imax)) 

ALLOCATE(URG(2:imax-1,jmax-1,4),ULG(2:imax-1,jmax-1,4), &
         URF(imax-1,2:jmax-1,4),ULF(imax-1,2:jmax-1,4), &
         Fpr(imax-1,2:jmax-1,4),Gpr(2:imax-1,jmax-1,4))   

CALL init_cond


!======= Numerical Scheme ===================================

IF (order == 1) THEN
DO n = 1,2

    CALL get_primitive
    CALL get_dt
    
    CALL lower_BC
    CALL upper_BC
    CALL outflow_BC

    CALL extr_Ustate_o1

    CALL AUSMPWpF
    CALL AUSMPWpG

    CALL update_state

END DO 


ELSE IF (order == 2) THEN
    WRITE(6,*) 'UNDER CONSTRUCTION, QUITTING..'
    STOP
ELSE 
    WRITE(6,*) '1st or 2nd ORDER ONLY PLEASE!'
    STOP
END IF

CALL write_sol

END PROGRAM s_sonic_main