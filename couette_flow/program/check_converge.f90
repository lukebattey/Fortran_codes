SUBROUTINE check_converge(n)

USE variables_cf

IMPLICIT NONE

INTEGER,INTENT(IN) :: n

seSS = 0.0    !
RMSeSS = 0.0  !<-- not making THIS mistake again...
seEX = 0.0    !
RMSeEX =0.0   !

t0 = dt0*(n)  !<-- counter for exact sol

Uex(1) = 0.00
Uex(jmax) = 1.00 

DO j = 2,jmax-1
    Uex(j) = Ynd(j) + SIN(pi*Ynd(j))*EXP(-(pi**2)*t0)
    seSS = seSS + (U(j) - Uss(j))**2
    seEX = seEX + (U(j) - Uex(j))**2
END DO

RMSeSS = SQRT(seSS/(jmax-2))
RMSeEX = SQRT(seEX/(jmax-2))

IF (RMSeEXmax < RMSeEX) THEN
RMSeEXmax = RMSeEX
END IF

IF(MOD(n,efreq) == 0 .or. n == 1) THEN
    WRITE(6,*) n,RMSeSS,RMSeEX
END IF

IF (RMSeSS <= ccRMSeSS) THEN
    WRITE(6,'(A)') ' '
    WRITE(6,'(A,I8,A,F8.6)') 'Converged in ',n,' iterations. time step = ',dt0
    WRITE(6,'(A)') 'RMS error history above: iteration,ssRMSer,ExRMSer'
    WRITE(6,*) 'Solution wrote to:  ',outfile
    WRITE(6,*) 'Maximum Exact RMS error below',RMSeEXmax
END IF  

END SUBROUTINE check_converge