SUBROUTINE CMA_heatEq_sub(n)

USE variables_cf

    IMPLICIT NONE

INTEGER,INTENT(IN) :: n

INTEGER :: i

DOUBLE PRECISION :: RTA
!------------------- Excact solution ---------------------
t0 = dt0*(n)

DO j = 1,jmax
    Uex(j) = Ynd(j) + SIN(pi*Ynd(j))*EXP(-(pi**2)*t0)
END DO 

!------------------- Initial Values! ---------------------
seSS   = 0.0    
RMSeSS = 0.0  
seEX   = 0.0    
RMSeEX = 0.0   

a1(:) = -1.0*theta*r
a2(:) =  1.0 + 2.0*theta*r  ! cool trick I learned..
a3(:) = -1.0*theta*r

a1(1)    = 0.0    !                           ^  
a1(jmax) = 0.0    !                           |
a2(1)    = 1.0    ! THESE NEED TO COME AFTER (:) 
a2(jmax) = 1.0    ! THE FULL ARRAY ASSIGNMENT ABOVE
a3(1)    = 0.0    !
a3(jmax) = 0.0    !

b(1)    = 0.0
b(jmax) = 1.0   !<--- upper plate speed = 1.0


DO j = 2,jmax-1
b(j) = (1.0 - theta)*r*(U(j+1)-2.0*U(j)+U(j-1)) + U(j)
END DO 

!-------------- Forward Thomas Sweep --------------------
DO j = 2,jmax
    RTA   = a1(j)/a2(j-1)
    a2(j) = a2(j) - RTA*a3(j-1)
    b(j)  = b(j) - RTA*b(j-1)
END DO

b(jmax) = b(jmax)/a2(jmax)

!-------------- Backward Thomas Sweep --------------------
DO j = jmax-1,1,-1
    b(j) = (b(j) - a3(j)*b(j+1))/a2(j) ! DONT FORGET IT'S (i) NOW...
END DO 

DO j = 2,jmax-1

    U(j) = b(j) !<-- Solution "stored" in b(j)...

    seSS = seSS + (U(j) - Uss(j))**2
    seEX = seEX + (U(j) - Uex(j))**2
END DO 

RMSeSS = SQRT(seSS/(jmax-2))
RMSeEX = SQRT(seEX/(jmax-2))

WRITE(6,*) n,RMSeSS,RMSeEX

END SUBROUTINE CMA_heatEq_sub