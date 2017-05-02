SUBROUTINE grid_read_metrics

USE variables_ss

    IMPLICIT NONE

!================== READ STUFF =============================================
OPEN(16,FILE = 'inputfile.dat', FORM = 'FORMATTED')
READ(16,*) infile 
READ(16,*) outfile
CLOSE(16)

OPEN(26,FILE = infile, FORM = 'FORMATTED')
    READ(26,'(A)') junk8 
    READ(26,'(A,I3,A,I3)') junk8,imax,junk8,jmax

    ALLOCATE(X(imax,jmax),Y(imax,jmax),Ja(imax,jmax),IJa(imax,jmax), &
            Xeta(imax,jmax),Xsi(imax,jmax),Yeta(imax,jmax),Ysi(imax,jmax), &
            etaX(imax,jmax),siX(imax,jmax),etaY(imax,jmax),siY(imax,jmax))

DO j = 1,jmax
    DO i = 1,imax
        READ(26,*) X(i,j),Y(i,j)
    END DO
END DO
CLOSE(26)

!============== Inverse Grid Metrics and Jacobians =========================
DO j=1,jmax  
  DO i=1,imax
    Xeta(i,j) = (X(i,j+1)-X(i,j-1))/2
    Yeta(i,j) = (Y(i,j+1)-Y(i,j-1))/2
    Xsi(i,j) = (X(i+1,j)-X(i-1,j))/2
    Ysi(i,j) = (Y(i+1,j)-Y(i-1,j))/2
    IF (j == jmax) THEN
      Xeta(i,j) = (3*X(i,j)-4*X(i,j-1)+X(i,j-2))/2
      Yeta(i,j) = (3*Y(i,j)-4*Y(i,j-1)+Y(i,j-2))/2
    END IF
    IF (j == 1) THEN
      Xeta(i,j) = -1*(3*X(i,j)-4*X(i,j+1)+X(i,j+2))/2
      Yeta(i,j) = -1*(3*Y(i,j)-4*Y(i,j+1)+Y(i,j+2))/2
    END IF
    IF (i == imax) THEN
      Xsi(i,j) = (3*X(i,j)-4*X(i-1,j)+X(i-2,j))/2
      Ysi(i,j) = (3*Y(i,j)-4*Y(i-1,j)+Y(i-2,j))/2
    END IF
    IF (i == 1) THEN
      Xsi(i,j) = -1*(3*X(i,j)-4*X(i+1,j)+X(i+2,j))/2
      Ysi(i,j) = -1*(3*Y(i,j)-4*Y(i+1,j)+Y(i+2,j))/2
    END IF
    Ja(i,j) = Xsi(i,j)*Yeta(i,j) - Xeta(i,j)*Ysi(i,j)
  END DO
END DO

DO j=1,jmax
    DO i=1,imax
        siX(i,j)  = Yeta(i,j) / Ja(i,j)
        siY(i,j)  = -1.0*Xeta(i,j) / Ja(i,j)
        etaX(i,j) = -1.0*Ysi(i,j) / Ja(i,j)
        etaY(i,j) = Xsi(i,j) / Ja(i,j)
    END DO
END DO

IJa(:,:) = 1.0 / Ja(:,:)  !<- This was a bug for me...

!================ WRITING A FILE TO CHECK STUFF=============================
OPEN(36,FILE = 'Jacobi_grid_plot.dat', FORM = 'FORMATTED')
WRITE(36,'(A)') 'VARIABLES = "X" "Y" "siX" "siY" "etaX" "etaY"'
WRITE(36,'(A,I3,A,I3)') 'ZONE I= ',imax,'    J=  ',jmax
DO j=1,jmax
  DO i=1,imax
    WRITE(36,*) X(i,j),Y(i,j),siX(i,j),siY(i,j),etaX(i,j),etaY(i,j)   ! comment this section if you want..
  END DO
END DO
CLOSE(36)
!===========================================================================

END SUBROUTINE grid_read_metrics