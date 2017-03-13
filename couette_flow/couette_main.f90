PROGRAM couette_main

USE variables_cf

    IMPLICIT NONE

INTEGER :: n

OPEN(16,FILE = 'inputfile.dat', FORM = 'FORMATTED')
READ(16,*) infile 
READ(16,*) outfile
READ(16,*) dt0
READ(16,*) jmax   ! always gotta have an input file
READ(16,*) ccRMSeSS
READ(16,*) theta
READ(16,*) nmax

OPEN(26,FILE = infile, FORM = 'FORMATTED')
    READ(26,*) jmax
    ALLOCATE(Ynd(jmax),U(jmax),Unext(jmax),Uss(jmax),Uex(jmax))

DO j = 1,jmax
    READ(26,*) U(j),Ynd(j)
END DO

!---------------- SS values ! Stay the same... -------------------
Uss(1) = 0
Uss(jmax) = 1.00 ! set to Plate speed??
DO j = 2,jmax-1
  Uss(j) = (Uss(jmax)-Uss(1))*(j-1)/(jmax-1)
END DO 
!----------------------------------------------------------------

dy0 = ABS(Ynd(3) - Ynd(2))   !THIS ASSUMES CONSTANT SPACING.. ok for now
r = dt0/(dy0**2)  ! coef of "scheme" term and stability indicator

WRITE(6,*) '------------- INFO ABOUT THIS CASE -------------------'
WRITE(6,*) '    dy0 =          dt0 =          r =          theta ='
WRITE(6,'(f15.7,f15.7,f15.7,f15.7)') dy0,dt0,r,theta     
WRITE(6,*) '------------------------------------------------------'
! Doing this can really help with de-bugging......

!---------------------------- LOOPS! -----------------------------
IF (theta == 0) THEN 
!---------------------- EXPLICIT SOLUTION! -----------------------
!=================================================================

    DO n = 1,nmax
        CALL EXPL_heatEq_sub(n)

        ! CALL convergence_check(n)

        IF (RMSeSS <= ccRMSeSS) THEN
        WRITE(6,'(A)') ' '
        WRITE(6,'(A,I8,A,F8.6)') 'Converged in ',n,&
        ' iterations. time step = ',dt0
        WRITE(6,'(A)') 'RMS error history above: iteration,ssRMSer,ExRMSer'
        WRITE(6,*) 'Solution wrote to:  ',outfile
            EXIT
        END IF  ! 142878 steps at 0.0001 
    END DO
!--------------------- COMBINED METHOD SOLUTION -----------------
ELSE
    
    ALLOCATE(a1(jmax),a2(jmax),a3(jmax),b(jmax))

    DO n = 1,nmax
        CALL CMA_heatEq_sub(n)

        IF (RMSeSS <= ccRMSeSS) THEN
        WRITE(6,'(A)') ' '
        WRITE(6,'(A,I8,A,F8.6)') 'Converged in ',n,&
        ' iterations. time step = ',dt0
        WRITE(6,'(A)') 'RMS error history above: iteration,ssRMSer,ExRMSer'
        WRITE(6,*) 'Solution wrote to:  ',outfile
            EXIT
        END IF  ! 142878 steps at 0.0001 for theta = 0.1
    END DO
END IF
!-------------------- WRITE RESULTS -----------------------------
IF (RMSeSS > ccRMSeSS) THEN
    WRITE(6,'(A)') ''
    WRITE(6,*) "DID NOT CONVERGE to Steady state!"
    WRITE(6,*) "Solution still wrote to:",outfile
END IF

OPEN(36,FILE = outfile, FORM = 'FORMATTED')
WRITE(36,*) jmax
DO j = 1,jmax
   WRITE(36,*) U(j),Ynd(j)
END DO

END PROGRAM couette_main