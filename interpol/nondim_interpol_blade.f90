PROGRAM nondim_interpol_blade
	IMPLICIT NONE

!------------------------------------------------------------
!			DEFINE VARIABLES
!------------------------------------------------------------

INTEGER,PARAMETER::rDef=SELECTED_REAL_KIND(12)
INTEGER,PARAMETER:: imax=16
INTEGER::i,numrR,il
REAL(KIND=rDef),DIMENSION(1:imax) :: r, rR, rRi, c, cR, cRi, tw, twi
REAL(KIND=rDef),PARAMETER:: Rmax = 5.03, &
							rRmin = 0.275, &
							rRmax = 0.975
REAL(KIND=rDef):: rRinc
numrR = 15

rRinc = (rRmax - rRmin)/(numrR - 1)


OPEN(16,FILE = 'Blade_geo_tab.txt', FORM = 'FORMATTED')
OPEN(26,FILE = 'Nondim_blade.txt', FORM = 'FORMATTED')

DO i=1,imax
	READ(16,*) r(i),c(i),tw(i)
	rR(i) = r(i)/Rmax
	cR(i) = c(i)/Rmax
	WRITE(26,*) rR(i),cR(i)
END DO

! make sure rR(1) is less than your first radial location (0.275)
! make sure rR(iMax) is greater than your last radial location (1)
! if not, print error and quit
IF (rR(1)<rRmin .AND. rR(imax)>rRmax) THEN
ELSE
WRITE(6,*) 'ERROR ERROR ERROR'
END IF

!go to 100

DO i=1,numrR
    rRi(i) = rRmin + rRinc*(i*1.0-1.0)
!WRITE(6,*) rRi(i)
! Increment il until rR(iLower)
	    DO WHILE (rRi(i) < rR(il))
        il = il + 1
	    END DO 
	    WRITE(6,*) il
	! Interpolate between rR(il) and rR(il + 1)
	!	cRi(i) = cR(i) + ((cR(i+1)-cR(i))/(rR(il+1)-rR(il)))*(rRi(i)-rR(il))
	!	WRITE(6,*) rRi(i),cRi(i)
END DO

!100 continue

END PROGRAM nondim_interpol_blade