PROGRAM Quality_cluster_check
    IMPLICIT NONE

INTEGER,PARAMETER:: rDef=SELECTED_REAL_KIND(12)
INTEGER::i,j,n,k,p,imax,jmax,nmax,itr,npj,np
REAL(KIND=rDef),ALLOCATABLE,DIMENSION(:,:) :: X,Y,Ja
CHARACTER(len=30) :: infile,outfile
CHARACTER(len=8) :: junk
LOGICAL :: srctrms
REAL(KIND=rDef) :: qcR,qcoX,qcoY,cRMS,Rc,Rcxs,Rcys,perc

OPEN(16,FILE = 'inputfile.dat', FORM = 'FORMATTED')
READ(16,*) infile 
READ(16,*) outfile
READ(16,*) srctrms
READ(16,*) nmax
READ(16,*) cRMS
READ(16,*) qcR
READ(16,*) qcoX,qcoY

OPEN(26,FILE = outfile, FORM = 'FORMATTED')

READ(26,'(A)') junk 
READ(26,'(A,I3,A,I3)') junk,imax,junk,jmax

ALLOCATE(X(imax,jmax),Y(imax,jmax),Ja(imax,jmax))

DO j = 1,jmax
    DO i = 1,imax
        READ(26,*) X(i,j),Y(i,j),Ja(i,j)
    END DO 
END DO 

  np = 0

DO j = 1,jmax
    npj = 0
    DO i = 1,imax
        Rcxs = (X(i,j) - qcoX)**2 
        Rcys = (Y(i,j) - qcoY)**2
        Rc = SQRT(Rcxs + Rcys)
        IF (Rc <= qcR) THEN 
            np = np + 1 
            npj = npj + 1
        END IF
    END DO
    IF (npj <= 0) THEN
         EXIT
    END IF
END DO 

perc = np*100.0/(imax*jmax)

WRITE(6,'(I4,A)') np,' grid points in the circular vicinity specified.'
WRITE(6,'(A,F7.3,A)') 'This is ',perc,' % of the total grid points.'
WRITE(6,'(A,A)') 'Check done for: ',outfile

END PROGRAM Quality_cluster_check