SUBROUTINE EXPL_heatEq_sub(n)

USE variables_cf

    IMPLICIT NONE

INTEGER,INTENT(IN) :: n

DO j = 2,jmax-1
    Unext(j) = U(j) + r*(U(j+1)-2*U(j)+U(j-1))
END DO

DO j = 2,jmax-1   !<-- needs to do this in another loop!
    U(j) = Unext(j)   ! (was a bug fix)
END DO 

END SUBROUTINE EXPL_heatEq_sub