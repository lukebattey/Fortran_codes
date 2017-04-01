SUBROUTINE update_state

USE variables_ss

IMPLICIT NONE

UstNEW(:,:,:) = 0.00

DO j = 2,jmax-1
    DO i = 2,imax-1

    NaNcheck = Gpr(i,j,1)
    IF (NaNcheck /= NaNcheck) THEN
        WRITE(6,*) i,j
    END IF

    UstNEW(i,j,:) = Ust(i,j,:) - dt*Ja(i,j)*( &
                    (Fpr(i,j,:) - Fpr(i-1,j,:)) + &
                    (Gpr(i,j,:) - Gpr(i,j-1,:)))

    Ust(i,j,:) = UstNEW(i,j,:)

    END DO
END DO

END SUBROUTINE update_state