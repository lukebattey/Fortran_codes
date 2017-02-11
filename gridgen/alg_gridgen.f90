PROGRAM alg_gridgen
	IMPLICIT NONE
!
!   This program writes 3 grids, the first and second grid are written excactly as specified in
!   the "Project #1" handout. The first grid does not cluster the grid points in the y direction
!   (Cy = 0.001) and is written to the "unclust_alg_grid.xy" output file. The second grid does 
!   cluster in the only the y direction (Cy = 2.0) and is written to the      
!
!----------------------------Variables -----------------------------
INTEGER,PARAMETER::rDef=SELECTED_REAL_KIND(10)
INTEGER,PARAMETER:: imax=41,jmax=19,fi=11,di=31,ci=41
INTEGER::i,j,g,F
REAL(KIND=rDef),DIMENSION(1:imax,1:jmax):: x,y
REAL(KIND=rDef),PARAMETER:: fx = -.8, &
                dx = 1.00, &
                cx = 1.80, &
                t = 0.15, &
                xint = 1.008930411365
REAL(KIND=rDef):: cy

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

WRITE(F,'(A)') 'VARIABLES = "X" "Y" '
WRITE(F,'(A,I3,A,I3)') 'ZONE I= ',imax,'    J=  ',jmax
DO j=1,jmax
	DO i=1,imax
        WRITE(F,*) x(i,j),y(i,j)
	END DO
END DO
!----------------- WRITE ALL POINTS ------------------------
END DO


END PROGRAM alg_gridgen

