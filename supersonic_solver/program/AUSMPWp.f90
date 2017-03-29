SUBROUTINE AUSMPWp

USE variables_ss
USE get_misc

IMPLICIT NONE

DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: UNtilL,UNtilR,VNtilL,VNtilR, &
                                               h0L,h0R,h0norm,Cs,Cave, &
                                               MtilL,MtilR,MLp,MRm,Pp,Pm,pL,pR, &
                                               fitalL,fitalR

DOUBLE PRECISION :: UNtilAve,ps,pmin,oneDP

oneDP = 1.00 ! For the intrinsic sign function because it's dumb...

!==========================================================================
!================= F' vectors: Fiphj = F'(i+1/2,j)) =======================
!==========================================================================

ALLOCATE(UNtilL(imax-1,jmax-2),UNtilR(imax-1,jmax-2), &     
         VNtilL(imax-1,jmax-2),VNtilR(imax-1,jmax-2), &
         h0L(imax-1,jmax-2),h0R(imax-1,jmax-2),h0norm(imax-1,jmax-2), &
         Cs(imax-1,jmax-2),Cave(imax-1,jmax-2), &
         MtilL(imax-1,jmax-2),MtilR(imax-1,jmax-2), &
         MLp(imax-1,jmax-2),MRm(imax-1,jmax-2), &
         Pp(imax-1,jmax-2),Pm(imax-1,jmax-2),pL(imax-1,jmax-2),pR(imax-1,jmax-2), &
         fitalL(imax-1,jmax-2),fitalR(imax-1,jmax-2))
         
! THE ALLOCATION IS DIFFERENT FOR THE G FLUXES!!!!! <-----

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

! Next: The pressures from the left and right extrapolated state vectors

pL(:,:) = (gama-1.0)*(ULF(:,:,4) - (ULF(:,:,2)**2 + ULF(:,:,3)**2) / &
          (2.0*ULF(:,:,1)))

pR(:,:) = (gama-1.0)*(ULF(:,:,4) - (ULF(:,:,2)**2 + ULF(:,:,3)**2) / &
          (2.0*ULF(:,:,1)))


DO j = 2,jmax-1    
    DO i = 1,imax-1   

      ps = Pp(i,j)*pL(i,j) + Pm(i,j)*pR(i,j)

      pmin = min(p(i,j-1), p(i,j+1), p(i+1,j-1), p(i+1,j+1))
      
    END DO
END DO


END SUBROUTINE AUSMPWp