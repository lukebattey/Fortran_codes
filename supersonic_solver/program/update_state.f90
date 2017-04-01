SUBROUTINE update_state

USE variables_ss

IMPLICIT NONE

UstNEW(:,:,:) = 0.00

DO j = 2,jmax-1
    DO i = 2,imax-1

    UstNEW(i,j,:) = Ust(i,j,:) - dt*Ja(i,j)*( &
                    (Fpr(i,j,:) - Fpr(i-1,j,:)) + &
                    (Gpr(i,j,:) - Fpr(i,j-1,:)))

    END DO
END DO

Ust(:,:,:) = UstNEW(:,:,:)

END SUBROUTINE update_state