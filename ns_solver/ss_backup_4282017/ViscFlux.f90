SUBROUTINE ViscFlux

USE variables_ss
USE get_misc

IMPLICIT NONE

DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: M,mu,Prndl,k, &
                                               usi,vsi,ueta,veta, &
                                               ux,vx,uy,vy, &
                                               Tsi,Teta,Tx,Ty, &
                                               qx,qy,tauxx,tauyy,tauxy

DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: Fvisc,Gvisc                        

DOUBLE PRECISION :: UNtilAve,ps,pmin,oneDP,wPw,Mave,T1sc,T2sc,T3sc,T4sc


ALLOCATE(M(imax,jmax),mu(imax,jmax),Prndl(imax,jmax), &
         k(imax,jmax), &
         usi(imax,jmax),vsi(imax,jmax), &
         ueta(imax,jmax),veta(imax,jmax), &
         ux(imax,jmax),vx(imax,jmax), &
         uy(imax,jmax),vy(imax,jmax), &
         Tsi(imax,jmax),Teta(imax,jmax), &
         Tx(imax,jmax),Ty(imax,jmax), &
         qx(imax,jmax),qy(imax,jmax))

ALLOCATE(Fvisc(1:imax-1,2:jmax-1,4), &
         Gvisc(2:imax-1,1:jmax-1,4))

!============= Get viscous stuff =======================================

Call get_primitive 

T(:,:)     = gama*(Minf**2)*p(:,:) / rho(:,:)

mu(:,:)    = (C1*((Tinf*T(:,:))**(1.5)) / (C2 + Tinf*T(:,:))) / muinf

k(:,:)     = (C3*((Tinf*T(:,:))**(1.5)) / (C4 + Tinf*T(:,:))) / kinf

Prndl(:,:) = Cp * mu(:,:)*muinf / (k(:,:)*kinf)


!=======================================================================
!============= Terms used for either Fvisc or Gvisc ====================
!=======================================================================

DO j=1,jmax  
  DO i=1,imax 
    ueta(i,j) = (u(i,j+1)-u(i,j-1))/2
    veta(i,j) = (v(i,j+1)-v(i,j-1))/2
    usi(i,j) = (u(i+1,j)-u(i-1,j))/2
    vsi(i,j) = (v(i+1,j)-v(i-1,j))/2
    Tsi(i,j) = (T(i+1,j)-T(i-1,j))/2
    Teta(i,j) = (T(i,j+1)-T(i,j-1))/2
    IF (j == jmax) THEN
      ueta(i,j) = (3*u(i,j)-4*u(i,j-1)+u(i,j-2))/2
      veta(i,j) = (3*v(i,j)-4*v(i,j-1)+v(i,j-2))/2
      Teta(i,j) = (3*T(i,j)-4*T(i,j-1)+T(i,j-2))/2
    END IF
    IF (j == 1) THEN
      ueta(i,j) = -1*(3*u(i,j)-4*u(i,j+1)+u(i,j+2))/2
      veta(i,j) = -1*(3*v(i,j)-4*v(i,j+1)+v(i,j+2))/2
      Teta(i,j) = -1*(3*T(i,j)-4*T(i,j+1)+T(i,j+2))/2
    END IF
    IF (i == imax) THEN
      usi(i,j) = (3*u(i,j)-4*u(i-1,j)+u(i-2,j))/2
      vsi(i,j) = (3*v(i,j)-4*v(i-1,j)+v(i-2,j))/2
      Tsi(i,j) = (3*T(i,j)-4*T(i-1,j)+T(i-2,j))/2
    END IF
    IF (i == 1) THEN
      usi(i,j) = -1*(3*u(i,j)-4*u(i+1,j)+u(i+2,j))/2
      vsi(i,j) = -1*(3*v(i,j)-4*v(i+1,j)+v(i+2,j))/2
      Tsi(i,j) = -1*(3*T(i,j)-4*T(i+1,j)+T(i+2,j))/2
    END IF
  END DO
END DO

ux = usi*siX + ueta*etaX 
vx = vsi*siX + veta*etaX 
uy = usi*siY + ueta*etaY 
vy = vsi*siY + veta*etaY 

tauxx = 2.0*mu*(2.0*ux - vy) / (3.0*ReL)

tauyy = 2.0*mu*(2.0*vy - ux) / (3.0*ReL)

tauxy = mu*(uy + vx) / ReL


!=======================================================================
!================= Fvisc' vectors: Fvisc'(i+1/2,j)) ====================
!=======================================================================

Tx = Tsi*siX + Teta*etaX

qx = -1.0*mu*Tx / ((gama-1.0)*(Minf**2)*ReL*Prndl)


DO j = 2,jmax-1
  DO i = 1,imax-1

    Fvisc(i,j,1) = 0.00

    Fvisc(i,j,2) =-0.5*(tauxx(i,j)+tauxx(i+1,j))

    Fvisc(i,j,3) =-0.5*(tauxy(i,j)+tauxy(i+1,j))

    Fvisc(i,j,4) =-0.5*(u(i,j)+u(i+1,j))*0.5*(tauxx(i,j)+tauxx(i+1,j)) &
                  -0.5*(v(i,j)+v(i+1,j))*0.5*(tauxy(i,j)+tauxy(i+1,j)) &
                  +0.5*(qx(i,j)+qx(i+1,j))
  END DO
END DO

!=======================================================================
!================= Gvisc' vectors: Gvisc'(i+1/2,j)) ====================
!=======================================================================

Ty = Tsi*siY + Teta*etaY

qy = -1.0*mu*Ty / ((gama-1.0)*(Minf**2)*ReL*Prndl)


DO j = 1,jmax-1
  DO i = 2,imax-1

    Gvisc(i,j,1) = 0.00

    Gvisc(i,j,2) =-0.5*(tauxy(i,j)+tauxy(i,j+1))

    Gvisc(i,j,3) =-0.5*(tauyy(i,j)+tauyy(i,j+1))

    Gvisc(i,j,4) =-0.5*(u(i,j)+u(i,j+1))*0.5*(tauxy(i,j)+tauxy(i,j+1)) &
                  -0.5*(v(i,j)+v(i,j+1))*0.5*(tauyy(i,j)+tauyy(i,j+1)) &
                  +0.5*(qy(i,j)+qy(i,j+1))
  END DO
END DO

!=======================================================================
!============== Updating Fluxes to include viscosity ===================
!=======================================================================

Fpr = Fpr + Fvisc
Gpr = Gpr + Gvisc

!================= Deallocate stuff to save some memory..===============
DEALLOCATE(M,mu,Prndl,k,usi,vsi,ueta,veta,ux,vx,uy,vy, &
           Tsi,Teta,Tx,Ty,qx,qy,Fvisc,Gvisc)

END SUBROUTINE ViscFlux