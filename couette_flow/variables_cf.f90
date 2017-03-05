MODULE variables_cf 
    IMPLICIT NONE

INTEGER,PARAMETER :: rDef=SELECTED_REAL_KIND(15)
REAL(KIND=rDef),PARAMETER :: pi=4.0*ATAN(1.0),root2=2.0**0.5
REAL(kind=rDef),ALLOCATABLE,DIMENSION(:) :: U,Unext
REAL(kind=rDef),ALLOCATABLE,DIMENSION(:) :: Ynd,Uss
                                                 
REAL(kind=rDef) :: nu,dt,t,dy,theta,tau,L,c,t0,y0,RMSerr,errS,SerrS

INTEGER :: j,jmax,n,nmax

!------ imin,jmin,imaxt,jmaxt  <-- these are just for testing


END MODULE variables_cf