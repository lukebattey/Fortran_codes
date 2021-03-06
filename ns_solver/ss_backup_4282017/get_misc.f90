MODULE get_misc

USE variables_ss
USE boundary_conds

IMPLICIT NONE

CONTAINS 
!=================== Get Freestream =======================================
! This subroutine gets some Freestream Parameters, should only be called 
! one time...

SUBROUTINE get_freestream

muinf = C1*((Tinf)**(1.5))/(Tinf + C2)
kinf  = C3*((Tinf)**(1.5))/(Tinf + C4)

Cp = gama*Rgas / (gama-1.0)
Cv = Rgas / (gama-1.0)

Prndlinf = Cp * muinf / kinf

Vinf = Minf*SQRT(gama*Rgas*Tinf)
ReL = rhoinf*Vinf*1.0 / muinf

END SUBROUTINE get_freestream



!=================== STAGNATION ENTHALPY =======================================
! Finds the stag enthalpy (h0) whenever called instant for the whole grid.
! Why did I even make this a subroutine????

SUBROUTINE stag_enthalpy
h0(:,:) = ((p(:,:)*gama) / (rho(:,:)*(gama-1.0))) + 0.5*(u(:,:)**2 + v(:,:)**2)
END SUBROUTINE stag_enthalpy



!=================== WRITE "SOLUTION" ==========================================
! Writes the solution for tecplot whenever called...
! this may need modifying to write it at a certain time step frequency (or something)
! FOR THE TIME BEING, THIS WRITES EACH ELEMENT OF THE STATE VECTOR AS Ust1,2,3,4...

SUBROUTINE write_sol
    ! OPEN(36,FILE = outfile, FORM = 'FORMATTED')
    ! WRITE(36,'(A)') 'VARIABLES = "X" "Y" "Ust1" "Ust2" "Ust3" "Ust4"'
    ! WRITE(36,'(A,I3,A,I3)') 'ZONE I= ',imax,'    J=  ',jmax
    ! DO j=1,jmax
    !     DO i=1,imax
    !         WRITE(36,*) X(i,j),Y(i,j),Ust(i,j,1),Ust(i,j,2),Ust(i,j,3),Ust(i,j,4)
    !     END DO
    ! END DO
    ! CLOSE(36)

    ! CALL cpu_time(EndTime)

    ! RunTime = EndTime - StartTime

    ! WRITE(6,*) ' '
    ! WRITE(6,'(A,I7,A,F9.4,A)') 'Converged with',n,' iterations in',RunTime," seconds" 
    ! WRITE(6,'(A)') 'Iteration vs. Residual RMSe history printed above'
    ! WRITE(6,*) 'Last state wrote to:  ',outfile


    OPEN(36,FILE = "primative_check.dat", FORM = 'FORMATTED')
    WRITE(36,'(A)') 'VARIABLES = "X" "Y" "p"' !"rho" "u" "v"
    WRITE(36,'(A,I3,A,I3)') 'ZONE I= ',imax,'    J=  ',jmax
    DO j=1,jmax
        DO i=1,imax
            WRITE(36,*) X(i,j),Y(i,j),p(i,j)
        END DO
    END DO
    CLOSE(36)

    WRITE(6,*) ' '
    WRITE(6,*) 'primitives wrote to:  primative_check.dat'

END SUBROUTINE write_sol



!=================== FIND TIME STEP dt ===================================
!This finds the time step to use as seen in eq 13 in the handout for Project #3

SUBROUTINE get_dt
    
    ALLOCATE(dti(imax,jmax),csqrt(imax,jmax),ci(imax,jmax)) ! Deallocates After..
    Util(:,:) = siX(:,:)*u(:,:) + siY(:,:)*v(:,:) 
    Vtil(:,:) = etaX(:,:)*u(:,:) + etaY(:,:)*v(:,:) 

    CALL stag_enthalpy 

    ci(:,:) = SQRT(2*h0(:,:)*(gama-1)/(gama+1))

    csqrt(:,:) = ci(:,:)*SQRT(siX(:,:)**2 + siY(:,:)**2 + etaX(:,:)**2 + etaY(:,:)**2 &
                              + 2.0*ABS(siX(:,:)*etaX(:,:) + siY(:,:)*etaY(:,:)))

    dti(:,:) = CFL / (ABS(Util(:,:)) + ABS(Vtil(:,:)) + csqrt(:,:))

    dt = 999999.99 !initiallize to a large number...
    DO j = 1,jmax
        DO i = 1,imax
            IF (dti(i,j) < dt) THEN  !dti becomes dt if dti is smaller
                dt = dti(i,j)        ! so it finds the minimum this way...
            END IF
        END DO
    END DO
    DEALLOCATE(dti,csqrt,ci) 

END SUBROUTINE get_dt



!=================== Pressure! ======================================
! Finds the pressure (p) whenever called for the whole grid.

SUBROUTINE get_pressure
    p(:,:) = (gama-1.0)*(Ust(:,:,4) - (Ust(:,:,2)**2 + Ust(:,:,3)**2) / &
             (2.0*Ust(:,:,1)))
END SUBROUTINE get_pressure



!=================== primitive variables! ======================================
! Finds the primitive variables (p,rho,u,v) whenever called for the whole grid.
! no need for temperature right now...

SUBROUTINE get_primitive
    
    call get_pressure

    rho(:,:) = Ust(:,:,1)
    u(:,:)   = Ust(:,:,2) / Ust(:,:,1)
    v(:,:)   = Ust(:,:,3) / Ust(:,:,1)

    call lower_BC

END SUBROUTINE get_primitive



!=================== minmod! ======================================
! ! excecutes the minmod function given as Eq. 4.220, pg 205 in Anderson (2013)

FUNCTION minmod(xm, ym)
   
   IMPLICIT NONE
   DOUBLE PRECISION, INTENT(IN) :: xm(4), ym(4)
   DOUBLE PRECISION, PARAMETER :: oe = 1.00
   DOUBLE PRECISION :: minmod(4)
   
   ! DO stind = 1,4
   !  minmod(stind) = sign(xm(stind),oe)* &
   !  max(0.,min(abs(xm(stind)),ym(stind)*sign(xm(stind),oe)))
   ! END DO
   
    DO stind = 1,4
        IF (xm(stind)*ym(stind) > 0.0) THEN
            IF (ABS(xm(stind)) < ABS(ym(stind))) THEN
                minmod(stind) = xm(stind)
            ELSE
                minmod(stind) = ym(stind)
            END IF
        ELSE IF (xm(stind)*ym(stind) <= 0.0) THEN
            minmod(stind) = 0.00
        END IF
    END DO

END FUNCTION minmod



!=================== Check convergence! ======================================
! Figure it out...

SUBROUTINE check_converge

DOUBLE PRECISION :: eSoS
eSoS = 0.000

DO j=1,jmax
    DO i=1,imax
        DO stind = 1,4
            eSoS = eSoS + (UstNew(i,j,1) - Ust(i,j,1))**2
        END DO
    END DO
END DO

    RMSe = SQRT(eSoS / (imax*jmax*4))

    IF (mod(n,wfrqRMSe) == 0) THEN
        WRITE(6,*) n,ABS(RMSe-lastRMSe)
    END IF

    converged = abs(RMSe - lastRMSe) <= convCrit


END SUBROUTINE check_converge




END MODULE get_misc