MODULE boundary_conds

USE variables_ss
USE get_misc

IMPLICIT NONE

CONTAINS 
!================ INFLOW BOUNDARY CONDITION ===================================
! finds the inflow "wall" values (i = imax). Called before each time step!

SUBROUTINE inflow_BC
    
   WRITE(6,*) '====================================   Attempted using inflow sub:'   
   WRITE(6,*) '=== UNDER CONSTRUCTION, QUITTING ===    I might not '
   WRITE(6,*) '====================================    even need it ever..'
    
    STOP
END SUBROUTINE inflow_BC



!================ OUTFLOW BOUNDARY CONDITION ==================================
! finds the outflow "wall" values (i = imax). Called before each time step!  
SUBROUTINE outflow_BC
    i = imax
        Ust(i,:,:) =  Ust(i-1,:,:)  ! i = imax !!!!! (right edge)
        rho(i,:)   =  rho(i-1,:)
        u(i,:)     =  u(i-1,:)
        v(i,:)     =  v(i-1,:)  
        p(i,:)     =  p(i-1,:)
END SUBROUTINE outflow_BC



!================ LOWER WALL BOUNDARY CONDITION ===============================
! finds the lower wall values (j = 1). Called before each time step!
SUBROUTINE lower_BC
        
    Util(:,:) = siX(:,:)*u(:,:) + siY(:,:)*v(:,:) 
    Vtil(:,:) = etaX(:,:)*u(:,:) + etaY(:,:)*v(:,:)  

    del2(:) = SQRT(siX(:,1)**2 + siY(:,1)**2) / &
              SQRT(siX(:,2)**2 + siY(:,2)**2)
    
    del3(:) = SQRT(siX(:,1)**2 + siY(:,1)**2) / &
              SQRT(siX(:,3)**2 + siY(:,3)**2)
!------------------------------------------------------------------------------

    u(:,1) = etaY(:,1)* &
             (2.0*Util(:,2)*del2(:) - Util(:,3)*del3(:)) / Ja(:,1)

    v(:,1) = -etaX(:,1)* &
             (2.0*Util(:,2)*del2(:) - Util(:,3)*del3(:)) / Ja(:,1)

    p(:,1) = p(:,2)

    rho(:,1) = p(:,1)*gama / &
               ((gama-1.0)*(h0(:,2) - 0.5*(u(:,1)**2 + v(:,1)**2))) 

    Ust(:,1,1) = rho(:,1)
    Ust(:,1,2) = rho(:,1)*u(:,1)
    Ust(:,1,3) = rho(:,1)*v(:,1)

    CALL stag_enthalpy  ! for Etot (4th el.), uses h0 based on the new values!

    Ust(:,1,4) = h0(:,1)*rho(:,1)

END SUBROUTINE lower_BC



!================ UPPER WALL BOUNDARY CONDITION ===============================
! finds the upper wall values (j = jmax). To be called before each time step!
SUBROUTINE upper_BC
   
    Util(:,:) = siX(:,:)*u(:,:) + siY(:,:)*v(:,:) 
    Vtil(:,:) = etaX(:,:)*u(:,:) + etaY(:,:)*v(:,:)

    deln1(:) = SQRT(siX(:,jmax)**2 + siY(:,jmax)**2) / &
               SQRT(siX(:,jmax-1)**2 + siY(:,jmax-1)**2)   
    
    deln2(:) = SQRT(siX(:,jmax)**2 + siY(:,jmax)**2) / &
               SQRT(siX(:,jmax-2)**2 + siY(:,jmax-2)**2)
!---------------------------------------------------------------------
    
    u(:,jmax) = etaY(:,jmax)* &
                (2.0*Util(:,jmax-1)*deln1(:) - Util(:,jmax-2)*deln2(:)) / &
                                                                     Ja(:,jmax)

    v(:,jmax) = -etaX(:,jmax)* &
                (2.0*Util(:,jmax-1)*deln1(:) - Util(:,jmax-2)*deln2(:)) / &
                                                                     Ja(:,jmax)

    p(:,jmax) = p(:,jmax-1)

    rho(:,jmax) = p(:,jmax)*gama / &
                  ((gama-1.0)*(h0(:,jmax-1) - 0.5*(u(:,jmax)**2 + v(:,jmax)**2)))

    Ust(:,jmax,1) = rho(:,jmax)
    Ust(:,jmax,2) = rho(:,jmax)*u(:,jmax)
    Ust(:,jmax,3) = rho(:,jmax)*v(:,jmax)

    CALL stag_enthalpy 
!   ^ room for optimization here, feed in the i or j line you want or something...

    Ust(:,jmax,4) = h0(:,jmax)*rho(:,jmax)

END SUBROUTINE upper_BC



END MODULE boundary_conds