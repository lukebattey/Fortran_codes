SUBROUTINE AUSMPWpG

USE variables_ss
USE get_misc

IMPLICIT NONE

DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: UNtilL,UNtilR,VNtilL,VNtilR, &
                                               h0L,h0R,h0norm,Cs,Cave, &
                                               MtilL,MtilR,MLp,MRm,Pp,Pm,pL,pR, &
                                               fitalL,fitalR,MbtpL,MbtmR

DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: Gt1,Gt2,Gt3,Gt4                         

DOUBLE PRECISION :: VNtilAve,ps,pmin,oneDP,wPw,Mave,T1sc,T2sc,T3sc,T4sc

oneDP = 1.00 ! For the intrinsic sign function because it's dumb...

!===================================================================================
!================= G' vectors: G'(i,j+1/2)) ================================
!===================================================================================

ALLOCATE(UNtilL(2:imax-1,jmax-1),UNtilR(2:imax-1,jmax-1), &     
         VNtilL(2:imax-1,jmax-1),VNtilR(2:imax-1,jmax-1), &
         h0L(2:imax-1,jmax-1),h0R(2:imax-1,jmax-1),h0norm(2:imax-1,jmax-1), &
         Cs(2:imax-1,jmax-1),Cave(2:imax-1,jmax-1), &
         MtilL(2:imax-1,jmax-1),MtilR(2:imax-1,jmax-1), &
         MLp(2:imax-1,jmax-1),MRm(2:imax-1,jmax-1), &
         Pp(2:imax-1,jmax-1),Pm(2:imax-1,jmax-1), &
         pL(2:imax-1,jmax-1),pR(2:imax-1,jmax-1), &
         fitalL(2:imax-1,jmax-1),fitalR(2:imax-1,jmax-1), &
         MbtpL(2:imax-1,jmax-1),MbtmR(2:imax-1,jmax-1))
! THE ALLOCATION IS DIFFERENT FOR THE F FLUXES!!!!! <-----^^^

ALLOCATE(Gt1(2:imax-1,jmax-1,4),Gt2(2:imax-1,jmax-1,4), &
         Gt3(2:imax-1,jmax-1,4),Gt4(2:imax-1,jmax-1,4))

DO j = 1,jmax-1    !<--- LOOPS ARE ALSO DIFFERENT FOR F FLUXES...
    DO i = 2,imax-1

        ! The following 4 variables: normalized contravarient
        ! velocities from left and right extrapolations. (Eq. 14)

        UNtilL(i,j) = ((siX(i,j)*ULG(i,j,2) / ULG(i,j,1)) + &   
                      (siY(i,j)*ULG(i,j,3) / ULG(i,j,1))) / &
                      SQRT(siX(i,j)**2 + siY(i,j)**2)

        UNtilR(i,j) = ((siX(i,j+1)*URG(i,j,2) / URG(i,j,1)) + & 
                      (siY(i,j+1)*URG(i,j,3) / URG(i,j,1))) / &
                      SQRT(siX(i,j+1)**2 + siY(i,j+1)**2)   

        VNtilL(i,j) = ((etaX(i,j)*ULG(i,j,2) / ULG(i,j,1)) + & 
                      (etaY(i,j)*ULG(i,j,3) / ULG(i,j,1))) / &
                      SQRT(etaX(i,j)**2 + etaY(i,j)**2)         

        VNtilR(i,j) = ((etaX(i,j+1)*URG(i,j,2) / URG(i,j,1)) + & 
                      (etaY(i,j+1)*URG(i,j,3) / URG(i,j,1))) / &
                      SQRT(etaX(i,j+1)**2 + etaY(i,j+1)**2)


        ! Next: stagnation enthalpy normal to the interface (Eq. 15)
        pL(i,j) = (gama-1.0)*(ULG(i,j,4) - (ULG(i,j,2)**2 + ULG(i,j,3)**2) / &
                  (2.0*ULG(i,j,1)))

        pR(i,j) = (gama-1.0)*(URG(i,j,4) - (URG(i,j,2)**2 + URG(i,j,3)**2) / &
                  (2.0*URG(i,j,1)))

        h0L(i,j) = pL(i,j)*gama/(ULG(i,j,1)*(gama-1.0)) + &
        0.5*((ULG(i,j,2)/ULG(i,j,1))**2 + (ULG(i,j,3)/ULG(i,j,1))**2)

        h0R(i,j) = pR(i,j)*gama/(URG(i,j,1)*(gama-1.0)) + &
        0.5*((URG(i,j,2)/URG(i,j,1))**2 + (URG(i,j,3)/URG(i,j,1))**2)


        h0norm(i,j) = 0.5*(h0L(i,j) - 0.5*UNtilL(i,j)**2 + &
                           h0R(i,j) - 0.5*UNtilR(i,j)**2)
        

        ! Next: cell-averaged speed of sound, Cave (Eq. 16)
        Cs(i,j) = SQRT(2.0*h0norm(i,j)*(gama-1) / (gama + 1))

        VNtilAve = 0.5*(VNtilL(i,j) + VNtilR(i,j)) !<-- NOT AN ARRAY...

        ! This IF statement is Eq. 16 in the handout
        IF (VNtilAve >= 0.00) THEN
            Cave(i,j) = Cs(i,j)**2 / MAX(cs(i,j),ABS(VNtilL(i,j)))
        ELSE
            Cave(i,j) = Cs(i,j)**2 / MAX(cs(i,j),ABS(VNtilR(i,j)))
        END IF

        ! NaNcheck = h0norm(i,j)
        IF (ULG(i,j,4) < 0.00) THEN
          WRITE(6,*) i,j,"ULG4 is < 0.00 (AUSMPWpG)"
        END IF

        ! Next: cell-face Mach numbers from L and R (Eq. 17) 
        MtilL(i,j) = VNtilL(i,j) / Cave(i,j)
        MtilR(i,j) = VNtilR(i,j) / Cave(i,j)


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

CALL get_pressure  ! This gets pressure for entire grid, no loop please..

! Next: The pressures from the left and right extrapolated state vectors (for Eq. 20)

pL(:,:) = (gama-1.0)*(ULG(:,:,4) - (ULG(:,:,2)**2 + ULG(:,:,3)**2) / &
    (2.0*ULG(:,:,1)))

pR(:,:) = (gama-1.0)*(URG(:,:,4) - (URG(:,:,2)**2 + URG(:,:,3)**2) / &
    (2.0*URG(:,:,1)))


! Next: pressure weighing "italic f" terms are found 
DO j = 1,jmax-1     !<================================= Loops are different too...
    DO i = 2,imax-1   

        ps = Pp(i,j)*pL(i,j) + Pm(i,j)*pR(i,j)  !<-- Eq. 20 here
                                                                 ! ============================
        pmin = min(p(i-1,j), p(i+1,j), p(i-1,j+1), p(i+1,j+1))   !<=== THIS IS DIFFERENT!!!====
                                                                 ! ============================
        ! The following IF statement is from e1                  !Everything else around here 
        IF (ps <= 0)  THEN                                       ! is same (optimize?)
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
!=================== ACTUALLY FINDING THE G FLUX's (Eq. 23) ===========================
!--------------------------------------------------------------------------------------

DO j = 1,jmax-1    
    DO i = 2,imax-1

        T1sc = MbtpL(i,j)*Cave(i,j)*SQRT(etaX(i,j)**2 + etaY(i,j)**2)*Ja(i,j)

        T2sc = MbtmR(i,j)*Cave(i,j)*SQRT(etaX(i,j+1)**2 + etaY(i,j+1)**2)*Ja(i,j+1)

        T3sc = Pp(i,j)*Ja(i,j)

        T4sc = Pm(i,j)*Ja(i,j+1)

        Gpr(i,j,1) = T1sc*ULG(i,j,1) + T2sc*URG(i,j,1) 

        Gpr(i,j,2) = T1sc*ULG(i,j,2) + T2sc*URG(i,j,2) + &
                     T3sc*etax(i,j)*pL(i,j) + T4sc*etax(i,j+1)*pR(i,j) 
        
        Gpr(i,j,3) = T1sc*ULG(i,j,3) + T2sc*URG(i,j,3) + &
                     T3sc*etaY(i,j)*pL(i,j) + T4sc*etaY(i,j+1)*pR(i,j)  
        
        Gpr(i,j,4) = T1sc*(ULG(i,j,4)+pL(i,j)) + T2sc*(URG(i,j,4)+pR(i,j))

    END DO
END DO


DEALLOCATE(UNtilL,UNtilR,VNtilL,VNtilR,h0L,h0R,h0norm,Cs,Cave, &
         MtilL,MtilR,MLp,MRm,Pp,Pm,pL,pR,fitalL,fitalR,MbtpL,MbtmR)
DEALLOCATE(Gt1,Gt2,Gt3,Gt4)


END SUBROUTINE AUSMPWpG