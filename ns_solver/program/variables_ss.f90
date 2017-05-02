MODULE variables_ss 

    IMPLICIT NONE   

INTEGER,PARAMETER :: rDef=SELECTED_REAL_KIND(10)
DOUBLE PRECISION,PARAMETER :: pi=4.0*ATAN(1.0),root2=2.0**0.5
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: X,Y,Ja,IJa,Xeta,Xsi,Yeta,Ysi, &
                                              etaX,siX,etaY,siY,u,v,p,rho, &
                                              Util,Vtil,h0,dti,csqrt,ci,T,e,M
                                             
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: del2,del3,deln1,deln2,pcheck

DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: Ust,Fpr,Gpr, &
                                                 URG,ULG, &
                                                 URF,ULF,UstNEW
                                                
CHARACTER(len=30) :: infile,outfile
CHARACTER(len=8) :: junk8 
                                   
DOUBLE PRECISION :: dt0,dy0,theta,L,r,t0,y0,seSS,RMSeSS,seEX, &
                    RMSeEX,ccRMSeSS,RMSeEXmax,ui,rhoi,gama,Eti, &
                    CFL,dt,NaNcheck,TempoRMS,lastRMSe,convCrit, &
                    RMSe,StartTime,EndTime,RunTime,ei

DOUBLE PRECISION :: Minf,Tinf,rhoinf,C1,C2,C3,C4,pinf,Ti,pinit, &
                    ReL,Rgas,muinf,kinf,Vinf,Cp,Cv,Prndlinf,RMSe1

LOGICAL :: fluxlim,converged,isoTwall,viscousFlux,viscousWall, &
            passconverge

INTEGER :: i,j,jmax,imax,n,order,stind,nmax,wfrqRMSe


END MODULE variables_ss