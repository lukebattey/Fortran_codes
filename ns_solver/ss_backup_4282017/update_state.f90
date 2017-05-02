SUBROUTINE update_state

USE variables_ss

IMPLICIT NONE

DO j = 2,jmax-1
    DO i = 2,imax-1

    UstNEW(i,j,:) = Ust(i,j,:) - dt*IJa(i,j)* &
                    ((Fpr(i,j,:) - Fpr(i-1,j,:)) + &
                     (Gpr(i,j,:) - Gpr(i,j-1,:)))

    END DO
END DO

END SUBROUTINE update_state