PROGRAM alg_gridgen
    IMPLICIT NONE     
!
!----------------------------Variables -----------------------------
INTEGER,PARAMETER::rDef=SELECTED_REAL_KIND(10)
INTEGER,PARAMETER:: imax=41,jmax=19,fi=11,di=31,ci=41
INTEGER::i,j,g,F
REAL(KIND=rDef),DIMENSION(1:imax,1:jmax):: x,y,Ja
REAL(KIND=rDef),PARAMETER:: fx = -.8, &
                dx = 1.00, &
                cx = 1.80, &
                t = 0.15, &
                xint = 1.008930411365
REAL(KIND=rDef):: cy,Xeta,Xsi,Yeta,Ysi

!----------------- Open files and set stipulations ---------------

OPEN(16,FILE = 'unclust_alg_grid.dat', FORM = 'FORMATTED')
OPEN(26,FILE = 'clust_alg_grid.dat', FORM = 'FORMATTED')

DO g=1,2

IF (g == 1) THEN
    cy = 0.001    ! cy : y-direction streching factor
    F = 16       ! F : file index to write to
ELSE
    cy = 2.0
    F = 26
END IF 

!----------------- WRITE LOWER BOUNDARY (j=1) ---------------------

DO i=1,fi-1
    x(i,1) = (-fx/(fi-1))*(i-1) + fx
    y(i,1) = 0
END DO

DO i=fi,di
    x(i,1) = (dx/(di-fi))*(i-fi)
    y(i,1) = 5*t*(0.2969*sqrt(xint*x(i,1)) &
        -0.126*xint*x(i,1) &
        -0.3516*((xint*x(i,1))**2) &
        +0.2843*((xint*x(i,1))**3) &
        -0.1015*((xint*x(i,1))**4)) 
END DO

DO i=di+1,ci
    x(i,1) = dx+((cx-dx)/(ci-di))*(i-di)
    y(i,1) = 0
END DO

!----------------- WRITE UPPER BOUNDARY (j=jmax) -------------------

DO i=1,ci
    x(i,jmax) = fx+((cx-fx)/(imax-1))*(i-1)
    y(i,jmax) = 1.0
END DO

!----------------- WRITE INTERIOR X POINTS  ------------------------

DO j=2,jmax-1
    DO i=1,imax
        x(i,j) = x(i,1) + ((j*1.0-1.0)/(jmax*1.0-1.0))*(x(i,jmax)-x(i,1))
        y(i,j) = y(i,1) &
        -(((y(i,jmax)-y(i,1))/cy))*log(1+(exp(-cy)-1)*((j*1.0-1.0)/(jmax*1.0-1.0)))
    END DO
END DO

DO j=1,jmax  
    DO i=1,imax
        Xeta = (X(i,j+1)-X(i,j-1))/2
        Yeta = (Y(i,j+1)-Y(i,j-1))/2
        Xsi = (X(i+1,j)-X(i-1,j))/2
        Ysi = (Y(i+1,j)-Y(i-1,j))/2
        IF (j == jmax) THEN
            Xeta = (3*X(i,j)-4*X(i,j-1)+X(i,j-2))/2
            Yeta = (3*Y(i,j)-4*Y(i,j-1)+Y(i,j-2))/2
        END IF
        IF (j == 1) THEN
             Xeta = -1*(3*X(i,j)-4*X(i,j+1)+X(i,j+2))/2
             Yeta = -1*(3*Y(i,j)-4*Y(i,j+1)+Y(i,j+2))/2
        END IF
        IF (i == imax) THEN
             Xsi = (3*X(i,j)-4*X(i-1,j)+X(i-2,j))/2
             Ysi = (3*Y(i,j)-4*Y(i-1,j)+Y(i-2,j))/2
        END IF
        IF (i == 1) THEN
             Xsi = -1*(3*X(i,j)-4*X(i+1,j)+X(i+2,j))/2
             Ysi = -1*(3*Y(i,j)-4*Y(i+1,j)+Y(i+2,j))/2
        END IF
        Ja(i,j) = Xsi*Yeta - Xeta*Ysi
    END DO
END DO

WRITE(F,'(A)') 'VARIABLES = "X" "Y" "Ja"'
WRITE(F,'(A,I3,A,I3)') 'ZONE I= ',imax,'    J=  ',jmax
    DO j=1,jmax
        DO i=1,imax
            WRITE(F,*) x(i,j),y(i,j),Ja(i,j)
        END DO
    END DO
    
END DO


END PROGRAM alg_gridgen

