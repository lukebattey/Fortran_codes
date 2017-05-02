PROGRAM s_sonic_main

USE variables_ss
USE boundary_conds
USE get_misc

    IMPLICIT NONE

call cpu_time(StartTime)

IF (order == 1 .and. fluxlim .eqv. .true.) THEN
    WRITE(6,*) "Flux limiters aren't used for 1st order cases..."
    WRITE(6,*) 'QUITTING'
    STOP
END IF

wfrqRMSe = 1
LastRMSe = 12345678.9
converged = .FALSE.

!======== READ INPUT FILE (all global variables) ============
OPEN(16,FILE = 'inputfile.dat', FORM = 'FORMATTED')
READ(16,*) infile 
READ(16,*) outfile
READ(16,*) Minf
READ(16,*) Tinf
READ(16,*) rhoinf
READ(16,*) gama
READ(16,*) Rgas
READ(16,*) C1
READ(16,*) C2
READ(16,*) C3
READ(16,*) C4
READ(16,*) CFL
READ(16,*) order
READ(16,*) viscousFlux
READ(16,*) viscousWall
READ(16,*) isoTwall
READ(16,*) fluxlim
READ(16,*) nmax
READ(16,*) passconverge
READ(16,*) convCrit
READ(16,*) wfrqRMSe
CLOSE(16) 

!======== Get some freestream constants ! ====================
ui = 1.0   !  It's unsafe to delete these..
rhoi = 1.0 !

Call get_freestream


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

!======================== TESTING LOOP ===========================
         DO n = 1,nmax
            CALL get_primitive
            CALL get_dt
            CALL stag_enthalpy

            CALL lower_BC
            CALL upper_BC
            CALL outflow_BC

            CALL upwind_U_LR_o2

            CALL AUSMPWpF
            CALL AUSMPWpG

            IF (viscousFlux) THEN
                CALL ViscFlux 
            END IF

            CALL update_state
            ! CALL check_converge

            ! if (converged) THEN
            !     IF (passconverge .eqv. .FALSE.) THEN 
            !         WRITE(6,*) n,ABS(RMSe-lastRMSe)
            !         CALL write_sol
            !         stop
            !     END IF
            ! end if 
            WRITE(6,*) 'woah'

            Ust(2:imax-1,2:jmax-1,:) = UstNEW(2:imax-1,2:jmax-1,:)
            lastRMSe = RMSe
        END DO
        

        call write_sol
STOP
!======================================================================


!======= Numerical Scheme ===================================
IF (order == 1) THEN
    DO n = 1,nmax
        CALL get_primitive
        CALL get_dt
        CALL stag_enthalpy
        CALL lower_BC
        CALL upper_BC
        CALL outflow_BC
        CALL extr_Ustate_o1 !<-- The only reason it's o1..
        CALL AUSMPWpF
        CALL AUSMPWpG

        IF (viscousFlux) THEN
            CALL ViscFlux 
        END IF

        CALL update_state
        CALL check_converge

        if (converged) THEN
            IF (passconverge .eqv. .FALSE.) THEN 
                WRITE(6,*) n,ABS(RMSe-lastRMSe)
                CALL write_sol
                stop
            END IF
        end if 

        Ust(2:imax-1,2:jmax-1,:) = UstNEW(2:imax-1,2:jmax-1,:)
        lastRMSe = RMSe
    END DO 

ELSE IF (order == 2) THEN

    IF (fluxlim .eqv. .TRUE.) THEN
        DO n = 1,nmax
            CALL get_primitive
            CALL get_dt
            CALL stag_enthalpy
            CALL lower_BC
            CALL upper_BC
            CALL outflow_BC
            CALL fluxlim_U_LR_o2 !<-- extrapo state vec w/ limiter
            CALL AUSMPWpF
            CALL AUSMPWpG

            IF (viscousFlux) THEN
                CALL ViscFlux 
            END IF

            CALL update_state
            CALL check_converge

            if (converged) THEN
                IF (passconverge .eqv. .FALSE.) THEN 
                    WRITE(6,*) n,ABS(RMSe-lastRMSe)
                    CALL write_sol
                    stop
                END IF
            end if 

            Ust(2:imax-1,2:jmax-1,:) = UstNEW(2:imax-1,2:jmax-1,:)
            lastRMSe = RMSe
        END DO
    ELSE  
        DO n = 1,nmax
            CALL get_primitive
            CALL get_dt
            CALL stag_enthalpy
            CALL lower_BC
            CALL upper_BC
            CALL outflow_BC
            CALL upwind_U_LR_o2 !<-- The only reason it's o2..
            CALL AUSMPWpF
            CALL AUSMPWpG

            IF (viscousFlux) THEN
                CALL ViscFlux 
            END IF

            CALL update_state
            CALL check_converge

            if (converged) THEN
                IF (passconverge .eqv. .FALSE.) THEN 
                    WRITE(6,*) n,ABS(RMSe-lastRMSe)
                    CALL write_sol
                    stop
                END IF
            end if 

            Ust(2:imax-1,2:jmax-1,:) = UstNEW(2:imax-1,2:jmax-1,:)
            lastRMSe = RMSe
        END DO
    END IF !<---fluxlim check
ELSE
    WRITE(6,*) '1st or 2nd order only please! Quitting..'
    STOP
END IF

!============== END OF MAIN PROGRAM ====================================

IF (passconverge) THEN
    IF (converged) THEN
        call write_sol
        WRITE(6,*) "IT DID CONVERGE BUT THEN WENT TO nmax THEN WROTE"
        STOP
    ELSE
        call write_sol
        WRITE(6,*) "DIDN'T ACTUALLY CONVERGE!!! Still wrote"
        STOP
    END IF
END IF

IF (n == nmax .AND. converged .eqv. .FALSE.) THEN
    WRITE(6,*) 'NOT REALLY CONVERGED, DONT TRUST THE FOLLOWING LINE!'
    CALL write_sol
    WRITE(6,*) "SOLUTION WRITTEN BUT IT IS NOOOOOT CONVERGED"
    STOP
END IF



END PROGRAM s_sonic_main