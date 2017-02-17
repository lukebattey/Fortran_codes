PROGRAM wT_ae6701
    IMPLICIT NONE

INTEGER,PARAMETER :: rDef=SELECTED_REAL_KIND(10)
INTEGER :: i,imax,numai
REAL(KIND=rDef),ALLOCATABLE,DIMENSION(:) :: rR
INTEGER,ALLOCATABLE,DIMENSION(:) :: ai  

OPEN(36,FILE = 'ae6071_WT.inp', FORM = 'FORMATTED')


DO i=1,imax
    READ(16,*) r(i),c(i),tw(i)
    rR(i) = r(i)/Rmax
    cR(i) = c(i)/Rmax
END DO

!--- Turn this into a general interpolation routine soon.... --------

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