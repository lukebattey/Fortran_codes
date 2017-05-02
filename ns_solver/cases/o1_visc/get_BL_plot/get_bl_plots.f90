PROGRAM report_plots

IMPLICIT NONE 

DOUBLE PRECISION,PARAMETER :: pi=4.0*ATAN(1.0),root2=2.0**0.5
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:) :: X,Y,p,rho,T,u,v,M
                                                
CHARACTER(len=30) :: infile,outfile
CHARACTER(len=8) :: junk8 
                        
DOUBLE PRECISION :: gama        
INTEGER :: i,j,jmax,imax,iplot,jmaxplot


!======== READ INPUT FILE (all global variables) =======================
OPEN(16,FILE = 'inputfile.dat', FORM = 'FORMATTED')
READ(16,*) infile 
READ(16,*) outfile
READ(16,*) iplot
READ(16,*) jmaxplot
READ(16,*) gama
CLOSE(16)


!========== READ THE INPUT FILE WITH THE DATA ======================
OPEN(26,FILE = infile, FORM = 'FORMATTED')
READ(26,'(A)') junk8 
READ(26,'(A,I3,A,I3)') junk8,imax,junk8,jmax

ALLOCATE(X(imax,jmax),Y(imax,jmax),p(imax,jmax),rho(imax,jmax), &
         T(imax,jmax),u(imax,jmax),v(imax,jmax),M(imax,jmax))

DO j=1,jmax
    DO i=1,imax
        READ(26,*) X(i,j),Y(i,j),p(i,j),rho(i,j),T(i,j),u(i,j),v(i,j),M(i,j)
    END DO
END DO
CLOSE(26)


!========== WRITE THE FINAL PLOT! ======================================
OPEN(36,FILE = outfile, FORM = 'FORMATTED')

i = iplot
DO j = 1,jmaxplot
    WRITE(36,*) T(i,j),u(i,j),Y(i,j) 
END DO

CLOSE(36)





END PROGRAM report_plots
