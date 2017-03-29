SUBROUTINE extr_Ustate_o1

    USE variables_ss

IMPLICIT NONE

! get URF and ULF (F is i+1/2,j)
DO j = 2,jmax-1
    DO i = 1,imax-1
        ULF(i,j,:) = Ust(i,j,:)
        URF(i,j,:) = Ust(i+1,j,:)
    END DO
END DO

! get URG and ULG (G is i,j+1/2). "LEFT" IS "ABOVE" NOW 
DO j = 1,jmax-1
    DO i = 2,imax-1
        URG(i,j,:) = Ust(i,j,:)    ! R is for "below"..
        ULG(i,j,:) = Ust(i,j+1,:)  ! L is for "above"..
    END DO
END DO
 
END SUBROUTINE