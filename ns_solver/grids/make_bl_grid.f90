PROGRAM bl_grid

IMPLICIT NONE 

DOUBLE PRECISION,PARAMETER :: pi=4.0*ATAN(1.0),root2=2.0**0.5
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: X,Y,p,rho,T,u,v,M
                                                
CHARACTER(len=30) :: infile,outfile
CHARACTER(len=8) :: junk8 
                        
DOUBLE PRECISION :: gama,dx
INTEGER :: i,j,jmax,imax,imaxplot,jmaxplot


!======== READ INPUT FILE (all global variables) =======================
OPEN(16,FILE = 'inputfile.dat', FORM = 'FORMATTED')
READ(16,*) infile 
READ(16,*) outfile


!========== READ THE INPUT FILE WITH THE DATA ======================
OPEN(26,FILE = infile, FORM = 'FORMATTED')
READ(26,'(A)') junk8 
READ(26,'(A,I3,A,I3)') junk8,imax,junk8,jmax

ALLOCATE(X(imax,jmax),Y(imax,jmax))

DO j=1,jmax
    DO i=1,imax
        READ(26,*) X(i,j),Y(i,j) !,p(i,j),rho(i,j),T(i,j),u(i,j),v(i,j),M(i,j)
    END DO
END DO
CLOSE(26)


!========== WRITE bl grid! ========================
imaxplot = 10
jmaxplot = 30

dx = 0.05

DO j = 1,jmaxplot
	DO i = 1,imaxplot
		X(i,j) = (i-1.0)*dx
		Y(i,j) = Y(1,j)
	END DO
END DO


OPEN(36,FILE = outfile, FORM = 'FORMATTED')
WRITE(36,'(A)') 'VARIABLES = "X" "Y"'
WRITE(36,'(A,I3,A,I3)') 'ZONE I= ',imaxplot,'    J=  ',jmaxplot
DO j=1,jmaxplot
  DO i=1,imaxplot
    WRITE(36,*) X(i,j),Y(i,j)
  END DO
END DO
CLOSE(36)



END PROGRAM bl_grid
