SUBROUTINE fluxlim_U_LR_o2

    USE variables_ss
    USE get_misc

IMPLICIT NONE

DOUBLE PRECISION,DIMENSION(4) :: delUp3h,delUm1h,delUp1h, &
                                 delBPm1h,delBMp3h, &
                                 xm,ym,UstGhost

DOUBLE PRECISION :: wm,oe,Ftest1,Gtest1,RMSt1,RMSt2
  
RMSt1 = 0.00
RMSt2 = 0.00 
  
wm = 1.3 ! hard coding this. I don't want to mess it up..
oe = 1.0000 ! this is just 1, so the compiler is happy...

! WRITE(6,*) "FLUX LIMITING EXTRAPOLATION CALLED!!!"
  
!=======================================================================
!=========== Extrapolate F fluxes! =====================================
!=======================================================================
DO j = 2,jmax-1
    DO i = 1,imax-1

        IF (i == 1) THEN ! this is probably a slow way of doing this... 

            UstGhost(:) = 2.0*Ust(i,j,:) - Ust(i+1,j,:)

            delUp1h(:) = Ust(i+1,j,:) - Ust(i,j,:)  
            delUm1h(:) = Ust(i,j,:) - UstGhost(:)
            delUp3h(:) = Ust(i+2,j,:) - Ust(i+1,j,:)  

        ELSE IF (i == imax-1) THEN

            UstGhost(:) = 2.0*Ust(i,j,:) - Ust(i-1,j,:)

            delUp1h(:) = Ust(i+1,j,:) - Ust(i,j,:)  
            delUm1h(:) = Ust(i,j,:) - Ust(i-1,j,:)
            delUp3h(:) = UstGhost(:) - Ust(i+1,j,:)

        ELSE

            delUp1h(:) = Ust(i+1,j,:) - Ust(i,j,:)  
            delUm1h(:) = Ust(i,j,:) - Ust(i-1,j,:)
            delUp3h(:) = Ust(i+2,j,:) - Ust(i+1,j,:)
            
        END IF

            xm(:) = delUm1h(:)
            ym(:) = wm*delUp1h(:)
            delBPm1h(:) = minmod(xm, ym)

            xm(:) = delUp3h(:)
            ym(:) = wm*delUp1h(:)
            delBMp3h(:) = minmod(xm, ym)

            ULF(i,j,:) = Ust(i,j,:) + 0.5*delBPm1h(:)
            URF(i,j,:) = Ust(i+1,j,:) - 0.5*delBMp3h(:)

    END DO
END DO

!=======================================================================
!=========== Extrapolate G fluxes! =====================================
!=======================================================================

DO j = 1,jmax-1                  
    DO i = 2,imax-1    ! UPDATE BELOW THIS TO HAVE FLUX LIMITER

        IF (j == 1) THEN

            UstGhost(:) = 2.0*Ust(i,j,:) - Ust(i,j+1,:)

            delUp1h(:) = Ust(i,j+1,:) - Ust(i,j,:)  
            delUm1h(:) = Ust(i,j,:) - UstGhost(:)
            delUp3h(:) = Ust(i,j+2,:) - Ust(i,j+1,:)

        ELSE IF (j == jmax-1) THEN

            UstGhost(:) = 2.0*Ust(i,j,:) - Ust(i,j-1,:)

            delUp1h(:) = Ust(i,j+1,:) - Ust(i,j,:)  
            delUm1h(:) = Ust(i,j,:) - Ust(i,j-1,:)
            delUp3h(:) = UstGhost(:) - Ust(i,j+1,:)

        ELSE

            delUp1h(:) = Ust(i,j+1,:) - Ust(i,j,:)  
            delUm1h(:) = Ust(i,j,:) - Ust(i,j-1,:)
            delUp3h(:) = Ust(i,j+2,:) - Ust(i,j+1,:)

        END IF

            xm(:) = delUm1h(:)
            ym(:) = wm*delUp1h(:)  
            delBPm1h(:) = minmod(xm, ym)

            xm(:) = delUp3h(:)
            ym(:) = wm*delUp1h(:)
            delBMp3h(:) = minmod(xm, ym)

            ULG(i,j,:) = Ust(i,j,:) + 0.5*delBPm1h(:)
            URG(i,j,:) = Ust(i,j+1,:) - 0.5*delBMp3h(:)

    END DO
END DO
 
END SUBROUTINE fluxlim_U_LR_o2