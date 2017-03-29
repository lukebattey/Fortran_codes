MODULE get_misc

USE variables_ss

IMPLICIT NONE

CONTAINS 
!=================== STAGNATION ENTHALPY ======================================
! Finds the stag enthalpy (h0) whenever called instant for the whole grid.
! Why did I even make this a subroutine????

SUBROUTINE stag_enthalpy
    h0(:,:) = p(:,:)*gama/(rho(:,:)*(gama-1.0)) + &
                                              (u(:,:)**2 + v(:,:)**2)/2.0
END SUBROUTINE stag_enthalpy



!=================== WRITE "SOLUTION" ==========================================
! Writes the solution for tecplot whenever called...
! this may need modifying to write it at a certain time step frequency (or something)
! FOR THE TIME BEING, THIS WRITES EACH ELEMENT OF THE STATE VECTOR AS Ust1,2,3,4...

SUBROUTINE write_sol
    OPEN(36,FILE = outfile, FORM = 'FORMATTED')
    WRITE(36,'(A)') 'VARIABLES = "X" "Y" "Ust1" "Ust2" "Ust3" "Ust4"'
    WRITE(36,'(A,I3,A,I3)') 'ZONE I= ',imax,'    J=  ',jmax
    DO j=1,jmax
        DO i=1,imax
            WRITE(36,*) X(i,j),Y(i,j),Ust(i,j,1),Ust(i,j,2),Ust(i,j,3),Ust(i,j,4)
        END DO
    END DO
    CLOSE(36)

    WRITE(6,'(A)') ''
    WRITE(6,*) 'Last state wrote to:  ',outfile
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
                              + 2*ABS(siX(:,:)*etaX(:,:) + siY(:,:)*etaY(:,:)))

    dti(:,:) = CFL / (Util(:,:) + Vtil(:,:) + csqrt(:,:))

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






END MODULE get_misc