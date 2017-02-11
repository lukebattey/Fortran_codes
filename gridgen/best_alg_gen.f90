PROGRAM best_alg_grid
    IMPLICIT NONE   

!----------------------------Variables -----------------------------
INTEGER,PARAMETER::rDef=SELECTED_REAL_KIND(10)
INTEGER,PARAMETER:: imax=41,jmax=19
INTEGER::i,j,ei,di,g,n
REAL(KIND=rDef),DIMENSION(1:imax,1:jmax):: x,y
REAL(KIND=rDef),PARAMETER:: fx = -.8, &
                            dx = 1.00, &
                            cx = 1.80, &
                            t = 0.15, &
                            xint = 1.008930411365
REAL(KIND=rDef):: cy,cxa,s

!----------------- Open files and set stipulations ---------------

OPEN(16,FILE = '.dat', FORM = 'FORMATTED')
OPEN(26,FILE = 'clust_alg_grid.dat', FORM = 'FORMATTED')

!----------------- WRITE LOWER BOUNDARY (j=1) ---------------------

ei = 10
di = 36
cxa = 3.0

x(1,1) = fx
y(1,1) = 0.00
x(ei,1) = 0.00 
y(ei,1) = 0.00
x(di,1) = 1.00 
y(di,1) = 0.00

! DO i=2,ei-1
    !WRITE(6,*) 
! END DO 
DO i=ei+1,ei+6
    n = i - ei
    X(i,1) = 0.0024*(n**2)+.0031*n-0.0006

    y(i,1) = 5*t*(0.2969*sqrt(xint*x(i,1)) &
        -0.126*xint*x(i,1) &
        -0.3516*((xint*x(i,1))**2) &
        +0.2843*((xint*x(i,1))**3) &
        -0.1015*((xint*x(i,1))**4))
END DO


! DO i=ei+1,di-1
!     x(i,1) = x(ei,1) &
!     -(((x(di,1)-x(ei,1))/cxa))*log(1+(exp(-cxa)-1)*((i*1.0-ei)/(di*1.0-ei)))

!     y(i,1) = 5*t*(0.2969*sqrt(xint*x(i,1)) &
!         -0.126*xint*x(i,1) &
!         -0.3516*((xint*x(i,1))**2) &
!         +0.2843*((xint*x(i,1))**3) &
!         -0.1015*((xint*x(i,1))**4)) 
! END DO


DO i=ei,ei+6
     WRITE(6,*) x(i,1), y(i,1)
END DO


! DO i=di+1,ci
!     x(i,1) = dx+((cx-dx)/(ci-di))*(i-di)
!     y(i,1) = 0
! END DO

! !----------------- WRITE UPPER BOUNDARY (j=jmax) -------------------

! DO i=1,ci
!     x(i,jmax) = fx+((cx-fx)/(imax-1))*(i-1)
!     y(i,jmax) = 1.0
! END DO

! !----------------- WRITE INTERIOR X POINTS  ------------------------

! DO j=2,jmax-1
! 	DO i=1,imax
!         x(i,j) = x(i,1) + ((j*1.0-1.0)/(jmax*1.0-1.0))*(x(i,jmax)-x(i,1))
!         y(i,j) = y(i,1) &
!         -(((y(i,jmax)-y(i,1))/cy))*log(1+(exp(-cy)-1)*((j*1.0-1.0)/(jmax*1.0-1.0)))
!     END DO
! END DO

! WRITE(F,'(A)') 'VARIABLES = "X" "Y" '
! WRITE(F,'(A,I3,A,I3)') 'ZONE I= ',imax,'    J=  ',jmax
! DO j=1,jmax
! 	DO i=1,imax
!         WRITE(F,*) x(i,j),y(i,j)
! 	END DO
! END DO

END PROGRAM best_alg_grid

