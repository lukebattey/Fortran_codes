MODULE variables_ss 

    IMPLICIT NONE

INTEGER,PARAMETER :: rDef=SELECTED_REAL_KIND(10)
REAL(KIND=rDef),PARAMETER :: pi=4.0*ATAN(1.0),root2=2.0**0.5
REAL(kind=rDef),ALLOCATABLE,DIMENSION(:,:) :: X,Y,Ja,Xeta,Xsi,Yeta,Ysi, &
                                              etaX,siX,etaY,siY,u,v,p,rho

REAL(kind=rDef),ALLOCATABLE,DIMENSION(:,:,:) :: Ust,Fpr,Gpr
                                                ! ^ State and Flux Vectors

CHARACTER(len=30) :: infile,outfile
CHARACTER(len=8) :: junk8 
                                   
REAL(kind=rDef) :: dt0,dy0,theta,L,r,t0,y0,seSS,RMSeSS,seEX, &
                   RMSeEX,ccRMSeSS,RMSeEXmax,ui,rhoi,gama

INTEGER :: i,j,jmax,imax

END MODULE variables_ss