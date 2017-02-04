PROGRAM findThoCoefs_0psiphi
	IMPLICIT NONE
!
!      
!----------------------------Variables ---------------------------------------------------
INTEGER,PARAMETER::rDef=SELECTED_REAL_KIND(12)
INTEGER::i,j,n,k,imax,jmax,nmax
REAL(KIND=rDef),ALLOCATABLE,DIMENSION(:,:,:) :: X,Xsi,Xsisi,Xeta,Xetaeta,Xetasi, &
                                                Y,Ysi,Ysisi,Yeta,Yetaeta,Yetasi, &
                                                alpha,beta,gama,psi,phi, &
                                                bbx,ddx,aax,ccx



!----------------- Open files and set stipulations ---------------
!------ I really want from her to XXX line to come from elliptic_grid_gen.f90 -------maybe the above stuff too?

OPEN(26,FILE = 'clust_alg_grid.xy', FORM = 'FORMATTED')

READ(26,*) imax,jmax
nmax = 50 ! assumes n wont exceed nmax...

ALLOCATE(X(imax,jmax,nmax),Xsi(imax,jmax,nmax),Xsisi(imax,jmax,nmax),Xeta(imax,jmax,nmax), &
         Xetaeta(imax,jmax,nmax),Xetasi(imax,jmax,nmax), &
         Y(imax,jmax,nmax),Ysi(imax,jmax,nmax),Ysisi(imax,jmax,nmax),Yeta(imax,jmax,nmax), &
         Yetaeta(imax,jmax,nmax),Yetasi(imax,jmax,nmax), &
         alpha(imax,jmax,nmax),beta(imax,jmax,nmax), &
         gama(imax,jmax,nmax),psi(imax,jmax,nmax),phi(imax,jmax,nmax))

n = 1                                                            
DO j=1,jmax                                                      
	DO i=1,imax                                                  
        READ(26,*) X(i,j,n),Y(i,j,n)
    END DO 
END DO
CLOSE(26) 

n=1
j=2 

!--------- This would already be set -----------------------------------------------------
DO n=2,3
! set boundaries to be constant
    DO i=1,imax
    X(i,1,n) = X(i,1,1)
    X(i,jmax,n) = X(i,jmax,1)
    Y(i,1,n) = Y(i,1,1)
    Y(i,jmax,n) = Y(i,jmax,1)
    END DO 
END DO

!---------------------------------------- XXX --------------------------------------------

DO i=2,imax-1
    Xeta(i,j,n) = (X(i,j+1,n)-X(i,j-1,n))/2
    Yeta(i,j,n) = (Y(i,j+1,n)-Y(i,j-1,n))/2
    alpha(i,j,n) = (Xeta(i,j,n)**2)+(Yeta(i,j,n)**2)

    Xsi(i,j,n) = (X(i+1,j,n)-X(i-1,j,n))/2
    Ysi(i,j,n) = (Y(i-1,j,n)-Y(i-1,j,n))/2
    gama(i,j,n) = (Xsi(i,j,n)**2)+(Ysi(i,j,n)**2)

    beta(i,j,n) = Xsi(i,j,n)*Xeta(i,j,n) + Ysi(i,j,n)*Yeta(i,j,n)

    bbx(i,j,n) = alpha(i,j,n) + gama(i,j,n)
    aax(i,j,n) = bbx(i,j,k)
    ccx(i,j,n) = beta(i,j,n)*(X(i+1,j+1,n)-X(i+1,j-1,n+1)-X(i-1,j+1,n)+X(i-1,j-1,n+1))
WRITE(6,*) bbx(i,j,n),ddx(i,j,n),aax(i,j,n),aax(i,j,n) 
END DO







END PROGRAM findThoCoefs_0psiphi



