SUBROUTINE upwind_U_LR_o2

    USE variables_ss

IMPLICIT NONE

DOUBLE PRECISION,DIMENSION(4) :: delUp3h,delUm1h,delUp1h

! get URF and ULF (F is i+1/2,j)
DO j = 2,jmax-1
    DO i = 1,imax-1

        IF (i == 1) THEN ! this is probably a slow way of doing this... 

            delUp1h(:) = Ust(i+1,j,:) - Ust(i,j,:)  
            ULF(i,j,:) = Ust(i,j,:) + 0.5*delUp1h(:) 

            delUp3h(:) = Ust(i+2,j,:) - Ust(i+1,j,:)
            URF(i,j,:) = Ust(i+1,j,:) - 0.5*delUp3h(:)

        ELSE IF (i == imax-1) THEN

            delUm1h(:) = Ust(i,j,:) - Ust(i-1,j,:)
            ULF(i,j,:) = Ust(i,j,:) + 0.5*delUm1h(:) 

            delUp1h(:) = Ust(i+1,j,:) - Ust(i,j,:)
            URF(i,j,:) = Ust(i+1,j,:) + 0.5*delUp1h(:)

        ELSE

            delUm1h(:) = Ust(i,j,:) - Ust(i-1,j,:)
            ULF(i,j,:) = Ust(i,j,:) + 0.5*delUm1h(:) 

            delUp3h(:) = Ust(i+2,j,:) - Ust(i+1,j,:)
            URF(i,j,:) = Ust(i+1,j,:) - 0.5*delUp3h(:)

        END IF

    END DO
END DO
                                    !================================
! get URG and ULG (G is i,j+1/2).   !===== UPDATE THIS TO BE O2 =====
DO j = 1,jmax-1                     !================================
    DO i = 2,imax-1

        IF (j == 1) THEN

            delUp1h(:) = Ust(i,j+1,:) - Ust(i,j,:)  
            ULG(i,j,:) = Ust(i,j,:) + 0.5*delUp1h(:) 

            delUp3h(:) = Ust(i,j+2,:) - Ust(i,j+1,:)
            URG(i,j,:) = Ust(i,j+1,:) - 0.5*delUp3h(:)

        ELSE IF (j == jmax-1) THEN

            delUm1h(:) = Ust(i,j,:) - Ust(i,j-1,:)
            ULG(i,j,:) = Ust(i,j,:) + 0.5*delUm1h(:) 

            delUp1h(:) = Ust(i,j+1,:) - Ust(i,j,:)
            URG(i,j,:) = Ust(i,j+1,:) + 0.5*delUp1h(:)

        ELSE

            delUm1h(:) = Ust(i,j,:) - Ust(i,j-1,:)
            ULG(i,j,:) = Ust(i,j,:) + 0.5*delUm1h(:) 

            delUp3h(:) = Ust(i,j+2,:) - Ust(i,j+1,:)
            URG(i,j,:) = Ust(i,j+1,:) - 0.5*delUp3h(:)

        END IF

        ! URG(i,j,:) = Ust(i,j+1,:)    ! R is for "above"...
        ! ULG(i,j,:) = Ust(i,j,:)  ! L is for "below"...
    END DO
END DO
 
END SUBROUTINE upwind_U_LR_o2