SUBROUTINE fluxlim_U_LR_o2

    USE variables_ss

IMPLICIT NONE

DOUBLE PRECISION,DIMENSION(4) :: delUp3h,delUm1h,delUp1h, &
                                 delBPm1h,delBMp3h, &
                                 xm,ym

DOUBLE PRECISION :: wm,oe,Ftest1,Gtest1,RMSt1,RMSt2

RMSt1 = 0.00
RMSt2 = 0.00

wm = 2.000 ! hard coding this. I don't want to mess it up..
oe = 1.0000 ! this is just 1, so the compiler is happy...

!=======================================================================
!=========== Extrapolate F fluxes! =====================================
!=======================================================================
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

            delUp1h(:) = Ust(i+1,j,:) - Ust(i,j,:)  
            delUm1h(:) = Ust(i,j,:) - Ust(i-1,j,:)
            delUp3h(:) = Ust(i+2,j,:) - Ust(i+1,j,:)

            xm(:) = delUm1h(:)
            ym(:) = wm*delUp1h(:)
            delBPm1h(:) = sign(xm(:),oe)*max(0.,min(abs(xm(:)),ym(:)*sign(xm(:),oe)))

            xm(:) = delUp3h(:)
            ym(:) = wm*delUp1h(:)
            delBMp3h(:) = sign(xm(:),oe)*max(0.,min(abs(xm(:)),ym(:)*sign(xm(:),oe)))

            ULF(i,j,:) = Ust(i,j,:) + 0.5*delBPm1h(:)

            URF(i,j,:) = Ust(i+1,j,:) + 0.5*delBMp3h(:)

            ! RMSt1 = RMSt1 + Ftest1**2
            ! RMSt2 = RMSt2 + ULF(i,j,1)**2

        END IF
    END DO
END DO
                                 
! RMSt1 = SQRT(RMSt1/((imax-3)*(jmax-2)))
! RMSt2 = SQRT(RMSt2/((imax-3)*(jmax-2)))

! WRITE(6,*) RMSt1,RMSt2,RMSt2-RMSt1

!=======================================================================
!=========== Extrapolate G fluxes! =====================================
!=======================================================================

DO j = 1,jmax-1                  
    DO i = 2,imax-1         ! UPDATE BELOW THIS TO HAVE FLUX LIMITER

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

            delUp1h(:) = Ust(i,j+1,:) - Ust(i,j,:)  
            delUm1h(:) = Ust(i,j,:) - Ust(i,j-1,:)
            delUp3h(:) = Ust(i,j+2,:) - Ust(i,j+1,:)

            xm(:) = delUm1h(:)
            ym(:) = wm*delUp1h(:)
            delBPm1h(1) = sign(xm(1),oe)*max(0.,min(abs(xm(1)),ym(1)*sign(xm(1),oe)))

            xm(:) = delUp3h(:)
            ym(:) = wm*delUp1h(:)
            delBMp3h(1) = sign(xm(1),oe)*max(0.,min(abs(xm(1)),ym(1)*sign(xm(1),oe)))

            ULG(i,j,:) = Ust(i,j,:) + 0.5*delBPm1h(:)
            URG(i,j,:) = Ust(i,j+1,:) + 0.5*delBMp3h(:)

        END IF
    END DO
END DO
 
END SUBROUTINE fluxlim_U_LR_o2