PROGRAM wT_ae6701
    IMPLICIT NONE

INTEGER,PARAMETER::rDef=SELECTED_REAL_KIND(12)
REAL(KIND=rDef):: R,rho,rpm,pitch
INTEGER :: vmax,vmin,vstep,V,B
INTEGER :: i,af,a,max,imax,numaf,nAlpha,nAlphaMax,n
CHARACTER(len=90),ALLOCATABLE,DIMENSION(:) :: airfile 
REAL(KIND=rDef),ALLOCATABLE,DIMENSION(:) :: rR,cR,tw
REAL,ALLOCATABLE,DIMENSION(:,:) :: alpha,cl,cd
INTEGER,ALLOCATABLE,DIMENSION(:) :: ai 
CHARACTER(len=1) :: junk
REAL :: numjunk,n1,n2,n3

OPEN(16,FILE = 'ae6701_WT.inp', FORM = 'FORMATTED')

READ(16,*) R
READ(16,*) B
READ(16,*) rho
READ(16,*) pitch
READ(16,*) rpm
READ(16,*) vmin,vmax,vstep

READ(16,*) numaf
ALLOCATE(airfile(numaf))

DO af=1,numaf
    READ(16,*) airfile(af)
END DO

READ(16,*) imax 
ALLOCATE(rR(imax),cR(imax),tw(imax),ai(imax))

DO i=1,imax
    READ(16,*) rR(i),cR(i),tw(i),ai(i)
END DO

! THIS FILE NEEDS TO BE IN CURRENT FOLDER, FORTRAN WONT READ "/", CAN THIS BE FIXED??? (below)
DO af=1,numaf
    OPEN(af*10,FILE = airfile(af), FORM = 'FORMATTED') !  arbitrary file handle....
!-------------------- read in junk to pass it --------------------------------------
    DO a = 1,3
        READ(af*10,*) junk ! this works 
    END DO 
    DO a = 1,9
       READ(af*10,*) numjunk ! this works too
    END DO

!------------------- read in alpha,cl, and cd (58 rows is standard..) ---------------  
    nAlphaMax = 58
    ALLOCATE(alpha(nAlphaMax,numaf),cl(nAlphaMax,numaf),cd(nAlphaMax,numaf))
    DO nAlpha = 1,nAlphaMax
        READ(af*10,*) alpha(nAlpha,af),cl(nAlpha,af),cd(nAlpha,af) 
        ! WRITE(6,*) alpha(nAlpha,af),cl(nAlpha,af),cd(nAlpha,af)
    END DO
END DO

V = 10
! DO V = vmin,vmax,vstep
    i = 1
    ! DO i = 1,imax
    
    ! END DO 
  
! END DO 

! DO i=1,imax
!     READ(16,*) r(i),c(i),tw(i)
!     rR(i) = r(i)/Rmax
!     cR(i) = c(i)/Rmax
! END DO

!--- CHECK IF INTERPOLATION IS POSSIBLE --------

! IF (rR(1)<rRmin .AND. rR(imax)>rRmax) THEN
! ELSE
!     WRITE(6,*) 'ERROR! WEE-WOO (siren sound) PROGRAM TERMINATED!'
!     STOP
! END IF

! DO i=1,numrR
!     rRi(i) = rRmin + rRinc*(i*1.0-1.0)
!     DO il=1,imax-1
!         IF (rR(il) < rRi(i) .AND. rR(il+1) > rRi THEN
!             cRi(i) = cR(il) + ((cR(il+1)-cR(il))/(rR(il+1)-rR(il)))*(rRi(i)-rR(il))
!             twi(i) = tw(il) + ((tw(il+1)-tw(il))/(rR(il+1)-rR(il)))*(rRi(i)-rR(il))
!         END IF
!     END DO
! END DO

END PROGRAM wT_ae6701