PROGRAM report_plots

IMPLICIT NONE 

DOUBLE PRECISION,PARAMETER :: pi=4.0*ATAN(1.0),root2=2.0**0.5
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: X,Y,p,M,h0,c,u,v
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: Ust
                                                
CHARACTER(len=30) :: infile,outfile
CHARACTER(len=8) :: junk8 
                        
DOUBLE PRECISION :: gama        
INTEGER :: i,j,jmax,imax


!======== READ INPUT FILE (all global variables) =======================
OPEN(16,FILE = 'inputfile.dat', FORM = 'FORMATTED')
READ(16,*) infile 
READ(16,*) outfile
READ(16,*) gama
CLOSE(16)



!========== READ THE INPUT FILE WITH STATE VECTOR ======================
OPEN(26,FILE = infile, FORM = 'FORMATTED')
READ(26,'(A)') junk8 
READ(26,'(A,I3,A,I3)') junk8,imax,junk8,jmax

ALLOCATE(X(imax,jmax),Y(imax,jmax),p(imax,jmax),M(imax,jmax), &
         h0(imax,jmax),c(imax,jmax),u(imax,jmax),v(imax,jmax))
ALLOCATE(Ust(imax,jmax,4))

DO j=1,jmax
    DO i=1,imax
        READ(26,*) X(i,j),Y(i,j),Ust(i,j,1),Ust(i,j,2),Ust(i,j,3),Ust(i,j,4)
    END DO
END DO
CLOSE(26)

p(:,:) = (gama-1.0)*(Ust(:,:,4) - (Ust(:,:,2)**2 + Ust(:,:,3)**2) / &
             (2.0*Ust(:,:,1)))

u(:,:) = Ust(:,:,2) / Ust(:,:,1)
v(:,:) = Ust(:,:,3) / Ust(:,:,1)

h0(:,:) = ((p(:,:)*gama) / (Ust(:,:,1)*(gama-1.0))) + &
           0.5*(u(:,:)**2 + v(:,:)**2)

c(:,:) = SQRT(2.0*h0(:,:)*(gama-1.0)/(gama+1.0))

M(:,:) = SQRT(u(:,:)**2 + v(:,:)**2) / c(:,:)

j = 1
DO i = 1,imax
	WRITE(6,*) X(i,j),p(i,j)
END DO

!========== WRITE THE FINAL PLOT! ======================================
OPEN(36,FILE = outfile, FORM = 'FORMATTED')
WRITE(36,'(A)') 'VARIABLES = "X" "Y" "p" "M"'
WRITE(36,'(A,I3,A,I3)') 'ZONE I= ',imax,'    J=  ',jmax
DO j=1,jmax
    DO i=1,imax
        WRITE(36,*) X(i,j),Y(i,j),p(i,j),M(i,j)
    END DO
END DO
CLOSE(36)

WRITE(6,*) ' '
WRITE(6,*) 'Done for: ',infile
WRITE(6,*) 'Wrote: ',outfile



END PROGRAM report_plots
