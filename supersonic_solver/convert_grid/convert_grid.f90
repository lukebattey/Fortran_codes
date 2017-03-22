PROGRAM convert_grid

    IMPLICIT NONE

INTEGER,PARAMETER :: rDef=SELECTED_REAL_KIND(10)
REAL(KIND=rDef),PARAMETER :: pi=4.0*ATAN(1.0),root2=2.0**0.5
REAL(kind=rDef),ALLOCATABLE,DIMENSION(:,:) :: X,Y,Ja,Xeta,Xsi,Yeta,Ysi, &
                                              etaX,siX,etaY,siY
CHARACTER(len=30) :: infile,outfile
CHARACTER(len=8) :: junk8 
REAL(kind=rDef) :: dt0,dy0,theta,L,r,t0,y0,seSS,RMSeSS,seEX, &
                   RMSeEX,ccRMSeSS,RMSeEXmax
INTEGER :: i,j,jmax,imax

!========== INPUT FILE =================================================
OPEN(16,FILE = 'inputfile.dat', FORM = 'FORMATTED')
READ(16,*) infile 
READ(16,*) outfile

!================ READ GIVEN GRID ======================================
 OPEN(26,FILE = infile, FORM = 'FORMATTED')

    imax = 71   ! You should never hard-code numbers...
    jmax = 48   ! This is the reason for converting it,
                ! and also to view in Tecplot!

    ALLOCATE(X(imax,jmax),Y(imax,jmax),Ja(imax,jmax),Xeta(imax,jmax), &
        Xsi(imax,jmax),Yeta(imax,jmax),Ysi(imax,jmax), &
        etaX(imax,jmax),siX(imax,jmax),etaY(imax,jmax),siY(imax,jmax))

DO i = 1,imax      ! I normally would do this the other way around, 
    DO j = 1,jmax  ! But this is how the grid is given....
        READ(26,*) X(i,j),Y(i,j)
    END DO
END DO
CLOSE(26)

!============= WRITE GRID IN "BETTER" FORMAT ===========================
OPEN(36,FILE = outfile, FORM = 'FORMATTED')
WRITE(36,'(A)') 'VARIABLES = "X" "Y"'
WRITE(36,'(A,I3,A,I3)') 'ZONE I= ',imax,'    J=  ',jmax
DO j=1,jmax
  DO i=1,imax
    WRITE(36,*) X(i,j),Y(i,j)
  END DO
END DO
CLOSE(36)
!=======================================================================

END PROGRAM convert_grid