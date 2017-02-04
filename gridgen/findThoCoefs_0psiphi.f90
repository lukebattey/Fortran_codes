PROGRAM findThoCoefs_0psiphi
    IMPLICIT NONE
!      
!----------------------------Variables ---------------------------------------------------
INTEGER,PARAMETER::rDef=SELECTED_REAL_KIND(12)
INTEGER::i,j,n,k,p,imax,jmax,nmax,itr
REAL(KIND=rDef),ALLOCATABLE,DIMENSION(:,:,:) :: X,Xsi,Xsisi,Xeta,Xetaeta, &
                                                Y,Ysi,Ysisi,Yeta,Yetaeta, &
                                                alpha,beta,gama,psi,phi, &
                                                bbx,ddx,aax,ccx, &
                                                bby,ddy,aay,ccy


!----------------- Open files and set stipulations ---------------
!------ I really want from her to XXX line to come from elliptic_grid_gen.f90 -------maybe the above stuff too?

OPEN(26,FILE = 'clust_alg_grid.xy', FORM = 'FORMATTED')

READ(26,*) imax,jmax
nmax = 50 ! assumes n wont exceed nmax...

ALLOCATE(X(imax,jmax,nmax),Xsi(imax,jmax,nmax),Xsisi(imax,jmax,nmax),Xeta(imax,jmax,nmax), &
         Xetaeta(imax,jmax,nmax), &
         Y(imax,jmax,nmax),Ysi(imax,jmax,nmax),Ysisi(imax,jmax,nmax),Yeta(imax,jmax,nmax), &
         Yetaeta(imax,jmax,nmax), &
         alpha(imax,jmax,nmax),beta(imax,jmax,nmax), &
         gama(imax,jmax,nmax),psi(imax,jmax,nmax),phi(imax,jmax,nmax), &
         bbx(imax,jmax,nmax),ddx(imax,jmax,nmax),aax(imax,jmax,nmax),ccx(imax,jmax,nmax), &
         bby(imax,jmax,nmax),ddy(imax,jmax,nmax),aay(imax,jmax,nmax),ccy(imax,jmax,nmax))

n = 1                                                            
DO j=1,jmax                                                      
	DO i=1,imax                                                  
        READ(26,*) X(i,j,n),Y(i,j,n)
    !WRITE(6,*) X(i,j,n),Y(i,j,n)
    END DO 
END DO
CLOSE(26) 

!----------------------- set boundaries to be constant -----------------------------------
DO n=1,15
	DO i=1,imax
    X(i,1,n) = X(i,1,1)
    X(i,jmax,n) = X(i,jmax,1)
    Y(i,1,n) = Y(i,1,1)
    Y(i,jmax,n) = Y(i,jmax,1)
    !WRITE(6,*) X(i,1,n),X(i,jmax,n),Y(i,1,n),Y(i,jmax,n)
    END DO 

    DO j=1,jmax
    X(1,j,n) = X(1,j,1)
    X(imax,j,n) = X(imax,j,1)
    Y(1,j,n) = Y(1,j,1)
    Y(imax,j,n) = Y(imax,j,1)
    !WRITE(6,*) X(i,1,n),X(i,jmax,n),Y(i,1,n),Y(i,jmax,n)
    END DO
END DO



DO n = 1,10
DO j=2,jmax-1
!----------------------- get alpha, beta, and gama ---------------------------------------
	DO i=2,imax-1
    Xeta(i,j,n) = (X(i,j+1,n)-X(i,j-1,n))/2
    Yeta(i,j,n) = (Y(i,j+1,n)-Y(i,j-1,n))/2
    alpha(i,j,n) = ((Xeta(i,j,n))**2)+((Yeta(i,j,n))**2)

    Xsi(i,j,n) = (X(i+1,j,n)-X(i-1,j,n))/2
    Ysi(i,j,n) = (Y(i+1,j,n)-Y(i-1,j,n))/2 
    gama(i,j,n) = ((Xsi(i,j,n))**2)+((Ysi(i,j,n))**2)
    
    beta(i,j,n) = Xsi(i,j,n)*Xeta(i,j,n) + Ysi(i,j,n)*Yeta(i,j,n)

!-------------------Finding Thomas array values --------------------------------------------
    bbx(i,j,n) = alpha(i,j,n)
    ddx(i,j,n) = (-2)*(alpha(i,j,n) + gama(i,j,n))
    aax(i,j,n) = bbx(i,j,n) 
    ccx(i,j,n) = beta(i,j,n)*(X(i+1,j+1,n)-X(i+1,j-1,n+1)-X(i-1,j+1,n)+X(i-1,j-1,n+1))/2 &
                 - gama(i,j,n)*(X(i,j+1,n)+X(i,j-1,n+1))

    bby(i,j,n) = alpha(i,j,n)
    ddy(i,j,n) = (-2)*(alpha(i,j,n) + gama(i,j,n))
    aay(i,j,n) = bby(i,j,n) 
    ccy(i,j,n) = beta(i,j,n)*(Y(i+1,j+1,n)-Y(i+1,j-1,n+1)-Y(i-1,j+1,n)+Y(i-1,j-1,n+1))/2 &
                 - gama(i,j,n)*(Y(i,j+1,n)+Y(i,j-1,n+1))
	END DO

!----------------- Setting Thomas array boundaries -----------------------------------------
bbx(1,j,n) = 0
bbx(imax,j,n) = 0
ddx(1,j,n) = 1
ddx(imax,j,n) = 1        ! X array boundaries
aax(1,j,n) = 0
aax(imax,j,n) = 0
ccx(1,j,n) = x(1,j,n)
ccx(imax,j,n) = x(imax,j,n)

bby(1,j,n) = 0
bby(imax,j,n) = 0
ddy(1,j,n) = 1
ddy(imax,j,n) = 1        ! Y array boundaries
aay(1,j,n) = 0
aay(imax,j,n) = 0
ccy(1,j,n) = y(1,j,n)
ccy(imax,j,n) = y(imax,j,n)

!----------------------- Thomas algorithm for x values --------------------------------------

	DO i=2,imax
        ddx(i,j,n) = ddx(i,j,n) - (bbx(i,j,n)/ddx(i-1,j,n))*aax(i-1,j,n)  ! X forward sweep
        ccx(i,j,n) = ccx(i,j,n) - (bbx(i,j,n)/ddx(i-1,j,n))*ccx(i-1,j,n)

        ddy(i,j,n) = ddy(i,j,n) - (bby(i,j,n)/ddy(i-1,j,n))*aay(i-1,j,n)  ! Y forward sweep
        ccy(i,j,n) = ccy(i,j,n) - (bby(i,j,n)/ddy(i-1,j,n))*ccy(i-1,j,n)
	END DO

    ccx(imax,j,n) = ccx(imax,j,n)/ddx(imax,j,n)  
    ccy(imax,j,n) = ccy(imax,j,n)/ddy(imax,j,n)  

	DO i=imax,1,-1
        ccx(i,j,n) = (ccx(i,j,n) - aax(i,j,n)*ccx(i+1,j,n))/ddx(i,j,n)   ! X forward sweep
        X(i,j,n+1) = ccx(i,j,n)

        ccy(i,j,n) = (ccy(i,j,n) - aay(i,j,n)*ccy(i+1,j,n))/ddy(i,j,n)   ! Y forward sweep
        Y(i,j,n+1) = ccy(i,j,n)
	END DO

END DO !----j loop
END DO !----n loop

!--------------------------------- Writing Results ----------------------------------------

DO j=1,jmax
    DO i=1,imax
    WRITE(6,*) X(i,j,10),Y(i,j,10)
!bbx(i,j,k)*X(i-1,j,n)+ddx(i,j,n)*X(i,j,n)+aax(i,j,n)*X(i+1,j,n)-ccx(i,j,n) !,MAX((x(i,j,n)**2),y(i,j,n)**2)
    END DO
END DO

END PROGRAM findThoCoefs_0psiphi



