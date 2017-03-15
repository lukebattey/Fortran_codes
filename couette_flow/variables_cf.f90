MODULE variables_cf 
    IMPLICIT NONE

INTEGER,PARAMETER :: rDef=SELECTED_REAL_KIND(10)
REAL(KIND=rDef),PARAMETER :: pi=4.0*ATAN(1.0),root2=2.0**0.5
REAL(kind=rDef),ALLOCATABLE,DIMENSION(:) :: U,Unext,Ynd,Uss,Uex,a1,a2,a3,b
CHARACTER(len=30) :: infile,outfile
                                   
REAL(kind=rDef) :: dt0,dy0,theta,L,r,t0,y0,seSS,RMSeSS,seEX, &
                    RMSeEX,ccRMSeSS

INTEGER :: j,jmax,nmax,efreq

END MODULE variables_cf