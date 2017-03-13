SUBROUTINE EXPL_heatEq_sub(n)

USE variables_cf

    IMPLICIT NONE

INTEGER,INTENT(IN) :: n

seSS = 0.0    !
RMSeSS = 0.0  !<-- not making THIS mistake again...
seEX = 0.0    !
RMSeEX =0.0   !

t0 = dt0*(n)  !<-- counter for excact sol

DO j = 2,jmax-1
    Unext(j) = U(j) + r*(U(j+1)-2*U(j)+U(j-1))
    Uex(j) = Ynd(j) + SIN(pi*Ynd(j))*EXP(-(pi**2)*t0)
END DO

DO j = 2,jmax-1    !<-- needs to do this in another loop!
    U(j) = Unext(j) ! found this out the hard way.... :(
    seSS = seSS + (U(j) - Uss(j))**2
    seEX = seEX + (U(j) - Uex(j))**2
END DO 

RMSeSS = SQRT(seSS/(jmax-2))
RMSeEX = SQRT(seEX/(jmax-2))

WRITE(6,*) n,RMSeSS,RMSeEX

END SUBROUTINE EXPL_heatEq_sub