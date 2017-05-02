SUBROUTINE init_cond

USE variables_ss

IMPLICIT NONE

ALLOCATE(u(imax,jmax),v(imax,jmax),p(imax,jmax),rho(imax,jmax), &
         T(imax,jmax))

ALLOCATE(Ust(imax,jmax,4), &  
         UstNEW(imax,jmax,4)) 
!============ SET INITIAL PRIMITIVE VARIABLES =======================

        u(:,:)   = ui
        v(:,:)   = 0.00 ! always zero here, so hard-coding is okay..
        rho(:,:) = rhoi
        p(:,:)   = 1.0/(gama*(Minf**2))
        T(:,:)   = 1.0

!============ SET INITIAL STATE VECTOR ==============================
Eti = p(1,1)/(gama-1) + rhoi*(ui**2)/2

        Ust(:,:,1) = rhoi
        Ust(:,:,2) = rhoi*ui
        Ust(:,:,3) = 0.00    ! same reason as above...
        Ust(:,:,4) = Eti


END SUBROUTINE init_cond