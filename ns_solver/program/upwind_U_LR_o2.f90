SUBROUTINE upwind_U_LR_o2

    USE variables_ss

IMPLICIT NONE

DOUBLE PRECISION,DIMENSION(4) :: delUp3h,delUm1h,delUp1h,UstGhost

! get URF and ULF (F is i+1/2,j).=======================================
DO j = 2,jmax-1
    DO i = 1,imax-1
        IF (i == 1 .or. i == 2) THEN 
            ULF(i,j,:) = Ust(i,j,:)
            delUp3h(:) = Ust(i+2,j,:) - Ust(i+1,j,:)
            URF(i,j,:) = Ust(i+1,j,:) - 0.5*delUp3h(:)
        ELSE IF (i == imax-1 .or. i == imax-2) THEN
            delUm1h(:) = Ust(i,j,:) - Ust(i-1,j,:)
            ULF(i,j,:) = Ust(i,j,:) + 0.5*delUm1h(:)
            URF(i,j,:) = Ust(i+1,j,:) 
        ELSE
            delUm1h(:) = Ust(i,j,:) - Ust(i-1,j,:)
            ULF(i,j,:) = Ust(i,j,:) + 0.5*delUm1h(:) 
            delUp3h(:) = Ust(i+2,j,:) - Ust(i+1,j,:)
            URF(i,j,:) = Ust(i+1,j,:) - 0.5*delUp3h(:)
        END IF
    END DO
END DO
                                   
! get URG and ULG (G is i,j+1/2).=======================================  
DO j = 1,jmax-1                    
    DO i = 2,imax-1
        IF (j == 1 .or. j == 2) THEN
            ULG(i,j,:) = Ust(i,j,:)  
            delUp3h(:) = Ust(i,j+2,:) - Ust(i,j+1,:)
            URG(i,j,:) = Ust(i,j+1,:) - 0.5*delUp3h(:)
        ELSE IF (j == jmax-1 .or. j == jmax-2) THEN
            URG(i,j,:) = Ust(i,j+1,:)
            delUm1h(:) = Ust(i,j,:) - Ust(i,j-1,:)
            ULG(i,j,:) = Ust(i,j,:) + 0.5*delUm1h(:) 
        ELSE
            delUm1h(:) = Ust(i,j,:) - Ust(i,j-1,:)
            ULG(i,j,:) = Ust(i,j,:) + 0.5*delUm1h(:) 
            delUp3h(:) = Ust(i,j+2,:) - Ust(i,j+1,:)
            URG(i,j,:) = Ust(i,j+1,:) - 0.5*delUp3h(:)
        END IF
    END DO
END DO
 
END SUBROUTINE upwind_U_LR_o2