SUBROUTINE AUSMPWpF

USE variables_ss
USE get_misc

IMPLICIT NONE

DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: UNtilL,UNtilR,VNtilL,VNtilR, &
                                               h0L,h0R,h0norm,Cs,Cave, &
                                               MtilL,MtilR,MLp,MRm,Pp,Pm,pL,pR, &
                                               fitalL,fitalR,MbtpL,MbtmR

DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: Ft1,Ft2,Ft3,Ft4                         

DOUBLE PRECISION :: UNtilAve,ps,pmin,oneDP,wPw,Mave,T1sc,T2sc,T3sc,T4sc

oneDP = 1.00 ! For the intrinsic sign function because it's dumb...

!===================================================================================
!================= F' vectors: Fiphj = F'(i+1/2,j)) ================================
!===================================================================================

ALLOCATE(UNtilL(imax-1,2:jmax-1),UNtilR(imax-1,2:jmax-1), &     
         VNtilL(imax-1,2:jmax-1),VNtilR(imax-1,2:jmax-1), &
         h0L(imax-1,2:jmax-1),h0R(imax-1,2:jmax-1),h0norm(imax-1,2:jmax-1), &
         Cs(imax-1,2:jmax-1),Cave(imax-1,2:jmax-1), &
         MtilL(imax-1,2:jmax-1),MtilR(imax-1,2:jmax-1), &
         MLp(imax-1,2:jmax-1),MRm(imax-1,2:jmax-1), &
         Pp(imax-1,2:jmax-1),Pm(imax-1,2:jmax-1), &
         pL(imax-1,2:jmax-1),pR(imax-1,2:jmax-1), &
         fitalL(imax-1,2:jmax-1),fitalR(imax-1,2:jmax-1), &
         MbtpL(imax-1,2:jmax-1),MbtmR(imax-1,2:jmax-1))
! THE ALLOCATION IS DIFFERENT FOR THE G FLUXES!!!!! <-----^^^

ALLOCATE(Ft1(imax-1,2:jmax-1,4),Ft2(imax-1,2:jmax-1,4), &
         Ft3(imax-1,2:jmax-1,4),Ft4(imax-1,2:jmax-1,4))

DO j = 2,jmax-1    !<--- LOOPS ARE ALSO DIFFERENT FOR G FLUXES...
    DO i = 1,imax-1

        ! The following 4 variables: normalized contravarient
        ! velocities from left and right extrapolations. (Eq. 14)

        UNtilL(i,j) = ((siX(i,j)*ULF(i,j,2) / ULF(i,j,1)) + &   
                      (siY(i,j)*ULF(i,j,3) / ULF(i,j,1))) / &
                      SQRT(siX(i,j)**2 + siY(i,j)**2)

        UNtilR(i,j) = ((siX(i+1,j)*URF(i,j,2) / URF(i,j,1)) + & 
                      (siY(i+1,j)*URF(i,j,3) / URF(i,j,1))) / &
                      SQRT(siX(i+1,j)**2 + siY(i+1,j)**2)   

        VNtilL(i,j) = ((etaX(i,j)*ULF(i,j,2) / ULF(i,j,1)) + & 
                      (etaY(i,j)*ULF(i,j,3) / ULF(i,j,1))) / &
                      SQRT(etaX(i,j)**2 + etaY(i,j)**2)         

        VNtilR(i,j) = ((etaX(i+1,j)*URF(i,j,2) / URF(i,j,1)) + & 
                      (etaY(i+1,j)*URF(i,j,3) / URF(i,j,1))) / &
                      SQRT(etaX(i+1,j)**2 + etaY(i+1,j)**2)


        ! Next: stagnation enthalpy normal to the interface (Eq. 15)
        h0norm(i,j) = 0.5*(ULF(i,j,4)/ULF(i,j,1) - 0.5*VNtilL(i,j)**2 + &
            URF(i,j,4)/URF(i,j,1) - 0.5*VNtilR(i,j)**2)

        IF (h0norm(i,j) <= 0.00) THEN
            WRITE(6,*) h0norm(i,j),' <--- h0norm < 0.. Quitting'
            STOP
        END IF 

        ! Next: cell-averaged speed of sound Cave (Eq. 16)
        Cs(i,j) = SQRT(2.0*h0norm(i,j)*(gama-1) / (gama + 1))

        UNtilAve = 0.5*(UNtilL(i,j) + UNtilR(i,j)) !<-- NOT AN ARRAY...

        ! This IF statement is Eq. 16 in the handout
        IF (UNtilAve >= 0.00) THEN
            Cave(i,j) = Cs(i,j)**2 / MAX(cs(i,j),ABS(UNtilL(i,j)))
        ELSE
            Cave(i,j) = Cs(i,j)**2 / MAX(cs(i,j),ABS(UNtilR(i,j)))
        END IF


        ! Next: cell-face Mach numbers from L and R (Eq. 17) 
        MtilL(i,j) = UNtilL(i,j) / Cave(i,j)
        MtilR(i,j) = UNtilR(i,j) / Cave(i,j)


        ! Next: split Mach numbers MLp and MRm (Eqs. 18a and 19a)
        IF (ABS(MtilL(i,j)) <= 1.0) THEN
            MLp(i,j) = 0.25*(MtilL(i,j) + 1.0)**2
        ELSE                                                !(18a)
            MLp(i,j) = 0.5*(MtilL(i,j) + ABS(MtilL(i,j)))
        END IF
        !--------------------------------------------------
        IF (ABS(MtilR(i,j)) <= 1.0) THEN
            MRm(i,j) = -0.25*(MtilR(i,j) - 1.0)**2
        ELSE                                                !(19a)
            MRm(i,j) = 0.5*(MtilR(i,j) - ABS(MtilR(i,j)))
        END IF


        ! Next: split Pressures Pp and Pm (plus and minus, Eqs. 18b and 19b)
        IF (ABS(MtilL(i,j)) <= 1.0) THEN
            Pp(i,j) = 0.25*((MtilL(i,j) + 1.0)**2)*(2.0 - MtilL(i,j))
        ELSE                                                       !(18b)
            Pp(i,j) = 0.5*(1.0 + sign(oneDP,MtilL(i,j)))
        END IF
        !--------------------------------------------------
        IF (ABS(MtilR(i,j)) <= 1.0) THEN
            Pm(i,j) = 0.25*((MtilR(i,j) - 1.0)**2)*(2 + MtilR(i,j))
        ELSE                                                !(19b)
            Pm(i,j) = 0.5*(1.0 - sign(oneDP,MtilR(i,j)))
        END IF       

    END DO 
END DO      

CALL get_pressure  ! This gets pressure for entire grid, no loop needed..

! Next: The pressures from the left and right extrapolated state vectors (for Eq. 20)

pL(:,:) = (gama-1.0)*(ULF(:,:,4) - (ULF(:,:,2)**2 + ULF(:,:,3)**2) / &
    (2.0*ULF(:,:,1)))

pR(:,:) = (gama-1.0)*(URF(:,:,4) - (URF(:,:,2)**2 + URF(:,:,3)**2) / &
    (2.0*URF(:,:,1)))

! Next: pressure weighing "italic f" terms are found 
DO j = 2,jmax-1    
    DO i = 1,imax-1   

        ps = Pp(i,j)*pL(i,j) + Pm(i,j)*pR(i,j)  !<-- Eq. 20 here

        pmin = min(p(i,j-1), p(i,j+1), p(i+1,j-1), p(i+1,j+1))

        ! The following IF statement is from eq. 21
        IF (ps <= 0)  THEN
            fitalL(i,j) = 0.00
            fitalR(i,j) = 0.00
        ELSE
            fitalL(i,j) = (pL(i,j)/ps - 1)*(min(1.0,(pmin/min(pL(i,j), pR(i,j))))**2)
            fitalR(i,j) = (pR(i,j)/ps - 1)*(min(1.0,(pmin/min(pL(i,j), pR(i,j))))**2)
        END IF

        

        ! Omega parameter from the pressure weighing terms (below Eq. 21)
        wPw = 1.0 - (min((pL(i,j)/pR(i,j)), (pR(i,j)/pL(i,j)) ))**3

        Mave = MLp(i,j) + MRm(i,j)

        IF (Mave >= 0.00) THEN
            MbtpL(i,j) = MLp(i,j) + MRm(i,j)*((1.0-wPw)*(1.0+fitalR(i,j))-fitalL(i,j))
            MbtmR(i,j) = MRm(i,j)*wPw*(1+fitalR(i,j))
        ELSE
            MbtpL(i,j) = MLp(i,j)*wPw*(1+fitalL(i,j))
            MbtmR(i,j) = MRm(i,j) + MLp(i,j)*((1.0-wPw)*(1.0+fitalL(i,j))-fitalR(i,j))
        END IF

    END DO
END DO


!--------------------------------------------------------------------------------------
!=================== ACTUALLY FINDING THE F FLUX's (Eq. 23) ===========================
!--------------------------------------------------------------------------------------

DO j = 2,jmax-1    
    DO i = 1,imax-1

        T1sc = MbtpL(i,j)*Cave(i,j)*SQRT(siX(i,j)**2 + siY(i,j)**2)*Ja(i,j)

        T2sc = MbtmR(i,j)*Cave(i,j)*SQRT(siX(i+1,j)**2 + siY(i+1,j)**2)*Ja(i+1,j)

        T3sc = Pp(i,j)*Ja(i,j)

        T4sc = Pm(i,j)*Ja(i+1,j)

        Fpr(i,j,1) = T1sc*ULF(i,j,1) + T2sc*URF(i,j,1) 

        Fpr(i,j,2) = T1sc*ULF(i,j,2) + T2sc*URF(i,j,2) + &
                     T3sc*six(i,j)*pL(i,j) + T4sc*six(i+1,j)*pR(i,j) 
        
        Fpr(i,j,3) = T1sc*ULF(i,j,3) + T2sc*URF(i,j,3) + &
                     T3sc*siY(i,j)*pL(i,j) + T4sc*siY(i+1,j)*pR(i,j)  
        
        Fpr(i,j,4) = T1sc*(ULF(i,j,4)+pL(i,j)) + T2sc*(URF(i,j,4)+pR(i,j))


        DO stind = 1,4
            NaNcheck = Fpr(i,j,stind)
            IF (NaNcheck /= NaNcheck) THEN
            WRITE(6,*) "F FLUX GONE BAD AT:",n,stind,i,j
            STOP
            END IF
        END DO

    END DO
END DO


DEALLOCATE(UNtilL,UNtilR,VNtilL,VNtilR,h0L,h0R,h0norm,Cs,Cave, &
         MtilL,MtilR,MLp,MRm,Pp,Pm,pL,pR,fitalL,fitalR,MbtpL,MbtmR)
DEALLOCATE(Ft1,Ft2,Ft3,Ft4)


END SUBROUTINE AUSMPWpF