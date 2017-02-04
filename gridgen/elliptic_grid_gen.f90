PROGRAM elliptic_grid_gen
	IMPLICIT NONE
!
!      
!----------------------------Variables -----------------------------
INTEGER,PARAMETER::rDef=SELECTED_REAL_KIND(12)
INTEGER::i,j,n,k,imax,jmax,nmax
REAL(KIND=rDef),ALLOCATABLE,DIMENSION(:,:,:) :: X,Xsi,Xsisi,Xeta,Xetaeta,Xetasi, &
                                                Y,Ysi,Ysisi,Yeta,Yetaeta,Yetasi, &
                                                alpha,beta,gama,psi,phi
                                                aa,ab,ac,ca,cb,cc,b,d

!----------------- Open files and set stipulations ---------------

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

DO n=2,3
! set boundaries to be constant
    DO i=1,imax
    X(i,1,n) = X(i,1,1)
    X(i,jmax,n) = X(i,jmax,1)
    Y(i,1,n) = Y(i,1,1)
    Y(i,jmax,n) = Y(i,jmax,1)
    END DO 
! loop over j (excluding boundaries)
        DO j=2,3
        




END DO

!DO j = 2,3
  



WRITE(6,*) "1st x and y are:", X(1,1,1),Y(1,1,1)

END PROGRAM elliptic_grid_gen
                   



