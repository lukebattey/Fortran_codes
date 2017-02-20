PROGRAM find_jacobians
    IMPLICIT NONE

INTEGER,PARAMETER:: rDef=SELECTED_REAL_KIND(12)
INTEGER::i,j,imax,jmax
REAL(KIND=rDef),ALLOCATABLE,DIMENSION(:,:) :: X,Y,Ja
REAL(KIND=rDef) :: Xsi,Xeta,Ysi,Yeta
CHARACTER(len=8) :: junk 

OPEN(16,FILE = 'elliptic_grid.dat', FORM = 'FORMATTED')

READ(16,'(A)') junk 
READ(16,'(A,I3,A,I3)') junk,imax,junk,jmax

ALLOCATE(X(imax,jmax),Y(imax,jmax),Ja(imax,jmax))

DO j=1,jmax
    DO i=1,imax 
        READ(16,*) X(i,j),Y(i,j)
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

WRITE(6,*) 'VARIABLES = "X" "Y" "Jac"'
WRITE(6,'(A,I3,A,I3)') 'ZONE I= ',imax,'    J=  ',jmax
DO j=1,jmax
    DO i=1,imax 
         IF (ja(i,j) <= 0) THEN
         WRITE(6,*) 'ALERT! BAD JACOBIAN SPOTTED! PROGRAM TERMINATED'
         STOP
         ELSE
         WRITE(6,*) X(i,j),Y(i,j),Ja(i,j)
         END IF 
    END DO
END DO 

END PROGRAM find_jacobians