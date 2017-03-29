MODULE variables_ss 

    IMPLICIT NONE   

INTEGER,PARAMETER :: rDef=SELECTED_REAL_KIND(10)
DOUBLE PRECISION,PARAMETER :: pi=4.0*ATAN(1.0),root2=2.0**0.5
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: X,Y,Ja,Xeta,Xsi,Yeta,Ysi, &
                                              etaX,siX,etaY,siY,u,v,p,rho, &
                                              Util,Vtil,h0,dti,csqrt,ci
                                             

DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: del2,del3,deln1,deln2

DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: Ust,Fpr,Gpr, &
                                                 URG,ULG, &
                                                 URF,ULF
                                                
CHARACTER(len=30) :: infile,outfile
CHARACTER(len=8) :: junk8 
                                   
DOUBLE PRECISION :: dt0,dy0,theta,L,r,t0,y0,seSS,RMSeSS,seEX, &
                   RMSeEX,ccRMSeSS,RMSeEXmax,ui,rhoi,gama,Eti, &
                   CFL,dt

LOGICAL :: fluxlim

INTEGER :: i,j,jmax,imax,n,order

END MODULE variables_ss