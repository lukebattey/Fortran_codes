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
REAL(KIND=rDef):: cy,cxa,s,cxl,cx1up,cx2up,fxup,exup,cy0,cyslope,dxu

!---------- assign point indecies & boundary clustering -----------------
ei = 15
di = 35
cxa = 2.0
cxl = 1.5

!----------------- WRITE LOWER BOUNDARY (j=1) ---------------------
x(1,1) = fx
y(1,1) = 0.00
x(ei,1) = 0.00 
y(ei,1) = 0.00
x(di,1) = 1.00 
y(di,1) = 0.00

DO i = 2,ei-1
    n = ei-i+1
    x(i,1) = x(ei,1)-((x(1,1)-x(ei,1)))* &
    log(1+(exp(-cxl)-1)*((n*1.0-1.0)/(ei-1.0)))/cxl
    y(i,1) = 0
END DO

DO i = ei+1,di-1
    x(i,1) = x(ei,1) &
    -(((x(di,1)-x(ei,1))/cxa))* &
    log(1+(exp(-cxa)-1)*((i*1.0-ei)/(di*1.0-ei)))
END DO 

DO i = ei+1,di-1
    y(i,1) = 5*t*(0.2969*sqrt(xint*x(i,1)) &
        -0.126*xint*x(i,1) &
        -0.3516*((xint*x(i,1))**2) &
        +0.2843*((xint*x(i,1))**3) &
        -0.1015*((xint*x(i,1))**4))
END DO 

dxu = 1.15
DO i=di+1,imax
    x(i,1) = dx+((cx-dx)/(imax-di))*(i-di)
    x(i,jmax) = dxu+((cx-dxu)/(imax-di))*(i-di)
    y(i,1) = 0
    y(i,jmax) = 1.00
END DO

!----------------- WRITE UPPER BOUNDARY (j=jmax) -------------------
exup = -.25
X(1,jmax) = fx
y(1,jmax) = 1.00
x(di,jmax) = dxu 
y(di,jmax) = 1.00

cx1up = 0.001
cx2up = 1.1

DO i = 2,di-1
    x(i,jmax) = x(1,jmax) &
      -(((x(di,jmax)-x(1,jmax))/cx2up))* &
    log(1+(exp(-cx2up)-1)*((i*1.0-1)/(di*1.0-1)))

    y(i,jmax) = 1.000
END DO
!----------------- WRITE INTERIOR X POINTS  ------------------------
cy0 = 2.00
cyslope = .9
DO j=2,jmax-1
  DO i=1,imax
        x(i,j) = x(i,1) + ((j*1.0-1.0)/(jmax*1.0-1.0))*(x(i,jmax)-x(i,1))
        IF (x(i,j) <= 0.00) THEN
        cy = cy0
        ELSE
        cy = cy0 - ABS(cyslope*x(i,j)) 
        END IF
        y(i,j) = y(i,1) &
        -(((y(i,jmax)-y(i,1))/cy))* &
        log(1+(exp(-cy)-1)*((j*1.0-1.0)/(jmax*1.0-1.0)))
    END DO
END DO
!--------------- WRITE FILE -----------------------------------------
OPEN(16,FILE = 'Best_alg.dat', FORM = 'FORMATTED')
WRITE(16,'(A)') 'VARIABLES = "X" "Y"'
WRITE(16,'(A,I3,A,I3)') 'ZONE I= ',imax,'    J=  ',jmax
DO j=1,jmax
  DO i=1,imax
        WRITE(16,*) x(i,j),y(i,j)
  END DO
END DO

END PROGRAM best_alg_grid

