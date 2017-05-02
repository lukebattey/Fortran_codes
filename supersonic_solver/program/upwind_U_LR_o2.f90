SUBROUTINE upwind_U_LR_o2

    USE variables_ss

IMPLICIT NONE

DOUBLE PRECISION,DIMENSION(4) :: delUp3h,delUm1h,delUp1h,UstGhost, &
                                 Up2ave,Um2ave

! get URF and ULF (F is i+1/2,j).=======================================
DO j = 2,jmax-1
    DO i = 1,imax-1

        IF (i == 1) THEN ! this is probably a slow way of doing this... 

            ! UstGhost(:) = 2*Ust(i,j,:) - Ust(i+1,j,:) ! Extrapolates U(-1,j,:)

            ! delUm1h(:) = Ust(i,j,:) - UstGhost(:)
            ! ULF(i,j,:) = Ust(i,j,:) + 0.5*delUm1h(:)

            ! delUp3h(:) = Ust(i+2,j,:) - Ust(i+1,j,:)
            ! URF(i,j,:) = Ust(i+1,j,:) - 0.5*delUp3h(:)

            ULF(i,j,:) = Ust(i,j,:)
            URF(i,j,:) = Ust(i+1,j,:)

        ELSE IF (i == imax-1) THEN
           
            ! delUm1h(:) = Ust(i,j,:) - Ust(i-1,j,:)
            ! ULF(i,j,:) = Ust(i,j,:) + 0.5*delUm1h(:) 

            ! UstGhost(:) = 2*Ust(i,j,:) - Ust(i-1,j,:) ! Extrapolates U(imax+1,:)

            ! delUp3h(:) = UstGhost(:) - Ust(i+1,j,:)
            ! URF(i,j,:) = Ust(i+1,j,:) - 0.5*delUp3h(:)

            ULF(i,j,:) = Ust(i,j,:)
            URF(i,j,:) = Ust(i+1,j,:)

        ELSE

            Um2ave(:) = 0.5*(Ust(i,j,:) + Ust(i-1,j,:))
            ULF(i,j,:) = 2.0*Ust(i,j,:) - Um2ave(:) 

            Up2ave(:) = 0.5*(Ust(i+1,j,:) + Ust(i+2,j,:))
            URF(i,j,:) = 2.0*Ust(i+1,j,:) - Up2ave(:)

        END IF

    END DO
END DO
                                   
! get URG and ULG (G is i,j+1/2).=======================================  
DO j = 1,jmax-1                    
    DO i = 2,imax-1

        IF (j == 1) THEN

            ! UstGhost(:) = 2.0*Ust(i,j,:) - Ust(i,j+1,:) !Extrapolates U(i,-1,:)

            ! delUm1h(:) = Ust(i,j,:) - UstGhost(:)
            ! ULG(i,j,:) = Ust(i,j,:) + 0.5*delUm1h(:) 

            ! delUp3h(:) = Ust(i,j+2,:) - Ust(i,j+1,:)
            ! URG(i,j,:) = Ust(i,j+1,:) - 0.5*delUp3h(:)

            URG(i,j,:) = Ust(i,j+1,:)
            ULG(i,j,:) = Ust(i,j,:)  

        ELSE IF (j == jmax-1) THEN

            ! delUm1h(:) = Ust(i,j,:) - Ust(i,j-1,:)
            ! ULG(i,j,:) = Ust(i,j,:) + 0.5*delUm1h(:) 

            ! UstGhost(:) = 2.0*Ust(i,j,:) - Ust(i,j-1,:) !Extrapolates U(i,jmax+1,:)

            ! delUp3h(:) = UstGhost(:) - Ust(i,j+1,:)
            ! URG(i,j,:) = Ust(i,j+1,:) - 0.5*delUp3h(:)

            URG(i,j,:) = Ust(i,j+1,:)
            ULG(i,j,:) = Ust(i,j,:) 

        ELSE

            ! delUm1h(:) = Ust(i,j,:) - Ust(i,j-1,:)
            ! ULG(i,j,:) = Ust(i,j,:) + 0.5*delUm1h(:) 

            ! delUp3h(:) = Ust(i,j+2,:) - Ust(i,j+1,:)
            ! URG(i,j,:) = Ust(i,j+1,:) - 0.5*delUp3h(:)

            Um2ave(:) = 0.5*(Ust(i,j,:) + Ust(i,j-1,:))
            ULG(i,j,:) = 2.0*Ust(i,j,:) - Um2ave(:) 

            Up2ave(:) = 0.5*(Ust(i,j-1,:) + Ust(i,j-2,:))
            URG(i,j,:) = 2.0*Ust(i,j-1,:) - Up2ave(:)

        END IF
    END DO
END DO
 
END SUBROUTINE upwind_U_LR_o2