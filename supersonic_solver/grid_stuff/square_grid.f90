PROGRAM square_grid
  IMPLICIT NONE
!----------------------------Variables ----------------------------------------
INTEGER,PARAMETER:: rDef=SELECTED_REAL_KIND(12)
INTEGER::i,j,n,k,p,imax,jmax,nmax,itr
REAL(KIND=rDef),ALLOCATABLE,DIMENSION(:,:) :: X,Y,Ja, &
alpha,beta,gama,bbx,ddx,aax,ccx,bby,ddy,aay,ccy,XNext,YNext,psi,phi
CHARACTER(len=30) :: infile,outfile
CHARACTER(len=8) :: junk
LOGICAL :: srctrms
REAL(KIND=rDef) :: cRMS,eRMS,XeSumS,YeSumS,Xsi,Xeta,Ysi,Yeta,&
XeRMS,YeRMS,Xes,Yes

!----------------------- Open files and set stipulations ----------------------
imax = 26
jmax = 26

ALLOCATE(X(imax,jmax),Y(imax,jmax),Ja(imax,jmax))
 
DO j=1,jmax                                                      
  DO i=1,imax                                                  

    X(i,j) = 1.0*i
    Y(i,j) = 1.0*j

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

  !--------------------- Writing Results to file for TecPlot360ex------------------------
OPEN(36,FILE = 'uniform_grid.dat', FORM = 'FORMATTED')
WRITE(36,'(A)') 'VARIABLES = "X" "Y" "J"'
WRITE(36,'(A,I3,A,I3)') 'ZONE I= ',imax,'    J=  ',jmax
DO j=1,jmax
  DO i=1,imax
    WRITE(36,*) X(i,j),Y(i,j),Ja(i,j)
  END DO
END DO
CLOSE(36)

END PROGRAM square_grid