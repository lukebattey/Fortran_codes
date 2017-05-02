SUBROUTINE init_cond

USE variables_ss

IMPLICIT NONE

ALLOCATE(u(imax,jmax),v(imax,jmax),p(imax,jmax),rho(imax,jmax), &
         T(imax,jmax),M(imax,jmax))

ALLOCATE(Ust(imax,jmax,4), &  
         UstNEW(imax,jmax,4)) 

!============ SET INITIAL PRIMITIVE VARIABLES =======================
        ui = 1.0
        rhoi = 1.0
        pinit = 1.0/(gama*(Minf**2)) 
        Ti = 1.0/(gama*(Minf**2))

        u(:,:)   = ui
        v(:,:)   = 0.00 ! always zero here, so hard-coding is okay..
        rho(:,:) = rhoi
        p(:,:)   = pinit
        T(:,:)   = Ti           

!============ SET INITIAL STATE VECTOR ==============================
Eti = p(1,1)/(gama-1) + 1.0*(1.0**2)/2

        Ust(:,:,1) = ui
        Ust(:,:,2) = ui*rhoi
        Ust(:,:,3) = 0.00    ! same reason as above...
        Ust(:,:,4) = Eti


END SUBROUTINE init_cond